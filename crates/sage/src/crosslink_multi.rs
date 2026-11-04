//! Multivalent crosslinks: species held together by
//! two up to `multivalent` linkers, beyond the pairs and loop-links of the
//! main search.
//!
//! # Species
//!
//! A species is a connected set of peptides joined by `k` links (no link from
//! a peptide to itself; loop-links stay in the main search). Copies of one
//! peptide count as separate peptides (`a2a`: two copies of `a`, as at a
//! homodimer interface). For `k` linkers there are 2 to `k + 1` peptides, in
//! every arrangement: several links between two peptides (`a2b`), chains
//! (`a1b1c`) and branches (`a1b,b1c,b1d`). The precursor is
//! `M = sum of peptide masses + k * linker`. A residue carries one linker.
//!
//! # Search
//!
//! 1. The best seeds of the main search (linear peptides explaining the most
//!    peaks, [`MULTI_SEEDS`] of them) are combined, with repetition, into all
//!    sets of `p - 1` peptides; the precursor then fixes the mass of the last
//!    peptide, looked up in the mass-sorted database.
//! 2. Every connected arrangement of `k` links on the `p` peptides is scored.
//!
//! # Fragment ions
//!
//! Remove peptide `X` from the species: the rest falls apart into groups, each
//! attached to `X` through one or more of its link residues. A backbone
//! fragment of `X` carries every group whose residues it holds all of (with
//! the group's peptides and linkers); one that holds only some of a group's
//! residues stays attached to the rest through the group and does not come
//! apart (as inside a loop-link), so it is not observed. Each group's mass
//! depends only on the arrangement, not on where its links sit on the other
//! peptides, so each peptide's residues are chosen on their own.
//!
//! Every peptide needs `min_chain_peaks` matched ions, `min_linked_peaks` of
//! them carrying a linker. The score is the hyperscore over all peptides.
//!
//! # MS-cleavable linkers (DSSO, DSBU, cliXlink)
//!
//! In the collision cell every linker of the species can break, freeing each
//! peptide with one stub on each of its link residues. A peptide's ions carry
//! one stub (any of its arm's) per link residue they hold; the free peptide
//! with its stubs gives doublets (triplets with two links), and every peptide
//! needs one such peak, one peptide two. Where peptides sit in the
//! arrangement is then invisible in the fragments, so arrangements are told
//! apart by the parts freed when a single linker breaks: a part of two or more
//! peptides with a stub (`b1c` + stub from `a1b1c`), counted as
//! [`Species::topology_peaks`]. cliXlink's arms differ, so each link is tried
//! both ways round.

use super::*;

/// Seeds combined into multivalent species, per spectrum and precursor charge
const MULTI_SEEDS: usize = 10;

/// Link-end role requirement: either end of the linker
const ANY: u8 = SITE1 | SITE2;

/// Every connected arrangement of `k` links (no self-links; repeats allowed)
/// on `p` numbered peptides, as edge lists `(i, j)` with `i < j`
pub(super) fn arrangements(p: usize, k: usize) -> Vec<Vec<(u8, u8)>> {
    let pairs: Vec<(u8, u8)> = (0..p as u8).flat_map(|a| (a + 1..p as u8).map(move |b| (a, b))).collect();
    let mut out = Vec::new();
    let mut edges = Vec::with_capacity(k);
    fn walk(pairs: &[(u8, u8)], from: usize, k: usize, p: usize, edges: &mut Vec<(u8, u8)>, out: &mut Vec<Vec<(u8, u8)>>) {
        if edges.len() == k {
            if connected(p, edges) {
                out.push(edges.clone());
            }
            return;
        }
        for i in from..pairs.len() {
            edges.push(pairs[i]);
            walk(pairs, i, k, p, edges, out);
            edges.pop();
        }
    }
    walk(&pairs, 0, k, p, &mut edges, &mut out);
    out
}

fn connected(p: usize, edges: &[(u8, u8)]) -> bool {
    let mut seen = vec![false; p];
    let mut stack = vec![0u8];
    seen[0] = true;
    while let Some(v) = stack.pop() {
        for &(a, b) in edges {
            let w = if a == v { b } else if b == v { a } else { continue };
            if !seen[w as usize] {
                seen[w as usize] = true;
                stack.push(w);
            }
        }
    }
    seen.iter().all(|&s| s)
}

/// The groups hanging off peptide `x` once it is removed from the
/// arrangement: per group, its mass (peptides, the linkers inside it and the
/// linkers joining it to `x`) and the links (edge indices) that join it to `x`
fn groups_of(x: u8, edges: &[(u8, u8)], masses: &[f32], linker: f32) -> Vec<(f32, Vec<usize>)> {
    let p = masses.len();
    // components of the arrangement without `x`
    let mut comp: Vec<usize> = (0..p).collect();
    fn find(comp: &mut [usize], v: usize) -> usize {
        let mut r = v;
        while comp[r] != r {
            r = comp[r];
        }
        comp[v] = r;
        r
    }
    for &(a, b) in edges {
        if a != x && b != x {
            let (ra, rb) = (find(&mut comp, a as usize), find(&mut comp, b as usize));
            comp[ra] = rb;
        }
    }
    let mut groups: Vec<(usize, f32, Vec<usize>)> = Vec::new();
    for (e, &(a, b)) in edges.iter().enumerate() {
        let other = match (a == x, b == x) {
            (true, false) => b,
            (false, true) => a,
            _ => continue,
        };
        let root = find(&mut comp, other as usize);
        match groups.iter_mut().find(|g| g.0 == root) {
            Some(g) => g.2.push(e),
            None => groups.push((root, 0.0, vec![e])),
        }
    }
    for g in &mut groups {
        let root = g.0;
        let members: Vec<usize> = (0..p).filter(|&v| v != x as usize && find(&mut comp, v) == root).collect();
        let inside = edges
            .iter()
            .filter(|&&(a, b)| members.contains(&(a as usize)) && members.contains(&(b as usize)))
            .count();
        g.1 = members.iter().map(|&v| masses[v]).sum::<f32>() + linker * (inside + g.2.len()) as f32;
    }
    groups.into_iter().map(|(_, mass, links)| (mass, links)).collect()
}

/// Best residues for one peptide of a species: `ends` are its link ends as
/// `(group, required role)`, `shifts` each group's mass. Returns the tally,
/// the matched ions carrying a linker and the residue of each end.
#[derive(Clone)]
struct ChainFit {
    tally: Tally,
    linked: u16,
    sites: Vec<usize>,
}

/// One peptide of a species whose MS-cleavable linkers all broke: its tally,
/// matched stub-carrying ions, the residue of each end, and how many of its
/// peptide + stubs peaks (doublets, triplets) are in the spectrum
#[derive(Clone)]
struct CleavedFit {
    tally: Tally,
    linked: u16,
    sites: Vec<usize>,
    doublet: u16,
}

/// Every injective choice of residues for `ends` (`(key, required role)`),
/// roles permitting; ends with the same key and role are interchangeable, so
/// their residues are kept in order
fn choose(
    depth: usize,
    ends: &[(usize, u8)],
    sites: &[(usize, u8)],
    used: &mut [bool],
    chosen: &mut Vec<usize>,
    visit: &mut dyn FnMut(&[usize]),
) {
    if depth == ends.len() {
        visit(chosen);
        return;
    }
    for (i, &(site, roles)) in sites.iter().enumerate() {
        let same_before = (0..depth).rev().find(|&d| ends[d].0 == ends[depth].0 && ends[d].1 == ends[depth].1);
        if used[i] || roles & ends[depth].1 == 0 || same_before.is_some_and(|d| chosen[d] > site) {
            continue;
        }
        used[i] = true;
        chosen.push(site);
        choose(depth + 1, ends, sites, used, chosen, visit);
        chosen.pop();
        used[i] = false;
    }
}

/// Masses of one stub on each of `arms` (any of an arm's short, long and alt
/// stubs), summed; distinct sums only
fn stub_sums(arms: &[&ArmStubs]) -> Vec<f32> {
    let mut sums = vec![0.0f32];
    for arm in arms {
        let mut next = Vec::with_capacity(sums.len() * 3);
        for s in &sums {
            for stub in arm.all() {
                let v = s + stub;
                if next.iter().all(|x: &f32| (x - v).abs() > 1e-3) {
                    next.push(v);
                }
            }
        }
        sums = next;
    }
    sums
}

/// Is there a (neutral) peak at `mass`?
fn has_peak(peaks: &[f32], tol: Tolerance, mass: f32) -> bool {
    let (lo, hi) = tol.bounds(mass);
    let i = peaks.partition_point(|&x| x < lo);
    i < peaks.len() && peaks[i] <= hi
}

/// Link ends of peptide `x` in an arrangement, as `(arm, required role)`
/// with the edge of each: `orient` bit `e` turns link `e` round (its first
/// peptide takes the second arm / end of the linker)
fn ends_of(x: u8, edges: &[(u8, u8)], orient: u32, asymmetric: bool, hetero: bool) -> (Vec<(usize, u8)>, Vec<usize>) {
    let (mut ends, mut edge_of_end) = (Vec::new(), Vec::new());
    for (e, &(a, b)) in edges.iter().enumerate() {
        if a != x && b != x {
            continue;
        }
        let flipped = orient & (1 << e) != 0;
        // this peptide holds the linker's first end?
        let first = (a == x) != flipped;
        let arm = if asymmetric && !first { 1 } else { 0 };
        let role = match (hetero, first) {
            (false, _) => ANY,
            (true, true) => SITE1,
            (true, false) => SITE2,
        };
        ends.push((arm, role));
        edge_of_end.push(e);
    }
    (ends, edge_of_end)
}

impl<'a, 'db> CrosslinkScorer<'a, 'db> {
    fn fit_chain(
        &self,
        query: &ProcessedSpectrum,
        peptide: &Peptide,
        shifts: &[f32],
        ends: &[(usize, u8)],
        max_fragment_charge: u8,
    ) -> Option<ChainFit> {
        let db = self.scorer.db;
        let tol = self.scorer.fragment_tol;
        let sites = self.settings.link_site_roles(peptide);
        if sites.len() < ends.len() {
            return None;
        }
        let mut best: Option<(f64, ChainFit)> = None;
        let mut chosen: Vec<usize> = Vec::with_capacity(ends.len());
        let mut used = vec![false; sites.len()];
        let min_peaks = self.settings.min_chain_peaks;
        let min_linked = self.settings.min_linked_peaks;
        let mut visit = |chosen: &[usize]| {
            let mut all: Vec<(Kind, f32)> = Vec::new();
            let mut shifted: Vec<(Kind, f32)> = Vec::new();
            for &kind in &db.ion_kinds {
                let n_terminal = matches!(kind, Kind::A | Kind::B | Kind::C);
                let contains = |site: usize, idx: usize| if n_terminal { idx >= site } else { idx < site };
                'ion: for (idx, ion) in IonSeries::new(peptide, kind, &db.modifications).enumerate() {
                    let mut mass = ion.monoisotopic_mass;
                    let mut carries = false;
                    for (g, &shift) in shifts.iter().enumerate() {
                        let mine: Vec<usize> = ends.iter().zip(chosen).filter(|(e, _)| e.0 == g).map(|(_, &s)| s).collect();
                        let held = mine.iter().filter(|&&s| contains(s, idx)).count();
                        if held == mine.len() {
                            mass += shift;
                            carries = true;
                        } else if held > 0 {
                            // tied to the rest through this group: does not come apart
                            continue 'ion;
                        }
                    }
                    all.push((kind, mass));
                    if carries {
                        shifted.push((kind, mass));
                    }
                }
            }
            let tally = match_ions(query, &all, max_fragment_charge, tol);
            if tally.matched() < min_peaks {
                return;
            }
            let linked = match_ions(query, &shifted, max_fragment_charge, tol).matched();
            if linked < min_linked {
                return;
            }
            let score = tally.score(self.scorer.score_type);
            if best.as_ref().is_none_or(|b| score > b.0) {
                best = Some((score, ChainFit { tally, linked, sites: chosen.to_vec() }));
            }
        };
        choose(0, ends, &sites, &mut used, &mut chosen, &mut visit);
        best.map(|(_, fit)| fit)
    }

    /// MS-cleavable linkers: one peptide of a species whose linkers all broke.
    /// Each of its link residues keeps one stub of its arm, so an ion holding
    /// some of them carries one stub per residue held. The free peptide with
    /// its stubs (`peaks`: neutral masses at reduced charge) gives doublets,
    /// or triplets with two links; at least one such peak is required.
    fn fit_cleaved_chain(
        &self,
        query: &ProcessedSpectrum,
        peptide: &Peptide,
        ends: &[(usize, u8)],
        arms: &[ArmStubs],
        peaks: &[f32],
        max_fragment_charge: u8,
    ) -> Option<CleavedFit> {
        let db = self.scorer.db;
        let tol = self.scorer.fragment_tol;
        let sites = self.settings.link_site_roles(peptide);
        if sites.len() < ends.len() {
            return None;
        }
        let end_arms: Vec<&ArmStubs> = ends.iter().map(|&(arm, _)| &arms[arm]).collect();
        let doublet = stub_sums(&end_arms)
            .into_iter()
            .filter(|&s| has_peak(peaks, tol, peptide.monoisotopic + s))
            .count() as u16;
        if doublet == 0 {
            return None;
        }
        let min_peaks = self.settings.min_chain_peaks;
        let mut best: Option<(f64, CleavedFit)> = None;
        let mut visit = |chosen: &[usize]| {
            let mut all: Vec<(Kind, f32)> = Vec::new();
            let mut stubbed: Vec<(Kind, f32)> = Vec::new();
            for &kind in &db.ion_kinds {
                let n_terminal = matches!(kind, Kind::A | Kind::B | Kind::C);
                let contains = |site: usize, idx: usize| if n_terminal { idx >= site } else { idx < site };
                for (idx, ion) in IonSeries::new(peptide, kind, &db.modifications).enumerate() {
                    let held: Vec<&ArmStubs> = chosen
                        .iter()
                        .zip(&end_arms)
                        .filter(|(&s, _)| contains(s, idx))
                        .map(|(_, &arm)| arm)
                        .collect();
                    if held.is_empty() {
                        all.push((kind, ion.monoisotopic_mass));
                        continue;
                    }
                    for sum in stub_sums(&held) {
                        all.push((kind, ion.monoisotopic_mass + sum));
                        stubbed.push((kind, ion.monoisotopic_mass + sum));
                    }
                }
            }
            let tally = match_ions(query, &all, max_fragment_charge, tol);
            if tally.matched() < min_peaks {
                return;
            }
            let score = tally.score(self.scorer.score_type);
            if best.as_ref().is_none_or(|b| score > b.0) {
                let linked = match_ions(query, &stubbed, max_fragment_charge, tol).matched();
                best = Some((score, CleavedFit { tally, linked, sites: chosen.to_vec(), doublet }));
            }
        };
        let (mut chosen, mut used) = (Vec::with_capacity(ends.len()), vec![false; sites.len()]);
        choose(0, ends, &sites, &mut used, &mut chosen, &mut visit);
        best.map(|(_, fit)| fit)
    }

    /// MS-cleavable linkers: score one arrangement of `peptides`. Every peptide
    /// needs a peptide + stubs peak and one of them a second (a doublet). Once
    /// every linker is broken the fragments no longer tell which peptides were
    /// linked; pieces where only one linker broke do: cutting a link that holds
    /// two parts of the species together frees each part with a stub
    /// (`b1c` + stub from `a1b1c`). Those peaks are counted as
    /// [`Species::topology_peaks`] and pick between arrangements.
    #[allow(clippy::too_many_arguments)]
    fn score_species_cleaved(
        &self,
        query: &ProcessedSpectrum,
        peptides: &[PeptideIx],
        masses: &[f32],
        edges: &[(u8, u8)],
        orient: u32,
        stubs: &CleavableStubs,
        peaks: &[f32],
        max_fc: u8,
        cache: &mut HashMap<(u32, Vec<(usize, u8)>), Option<CleavedFit>>,
    ) -> Option<CrosslinkMatch> {
        let db = self.scorer.db;
        let tol = self.scorer.fragment_tol;
        let linker = self.settings.linker_mass();
        let arms = stubs.arms();
        let asymmetric = arms.len() > 1;
        let hetero = self.settings.is_heterobifunctional();
        let mut fits: Vec<CleavedFit> = Vec::with_capacity(peptides.len());
        let mut end_edges: Vec<Vec<usize>> = Vec::with_capacity(peptides.len());
        for x in 0..peptides.len() as u8 {
            let (ends, edge_of_end) = ends_of(x, edges, orient, asymmetric, hetero);
            let fit = cache
                .entry((peptides[x as usize].0, ends.clone()))
                .or_insert_with(|| self.fit_cleaved_chain(query, &db[peptides[x as usize]], &ends, &arms, peaks, max_fc))
                .clone()?;
            fits.push(fit);
            end_edges.push(edge_of_end);
        }
        if fits.iter().all(|f| f.doublet < 2) {
            return None;
        }

        // parts freed by breaking one linker: link e joins two parts only if
        // nothing else does; a part of two or more peptides keeps its inner linkers
        let mut topology_peaks = 0u16;
        for (e, &(a, b)) in edges.iter().enumerate() {
            let rest: Vec<(u8, u8)> = edges.iter().enumerate().filter(|&(i, _)| i != e).map(|(_, &l)| l).collect();
            let mut side = vec![a];
            let mut grow = true;
            while grow {
                grow = false;
                for &(u, v) in &rest {
                    for (from, to) in [(u, v), (v, u)] {
                        if side.contains(&from) && !side.contains(&to) {
                            side.push(to);
                            grow = true;
                        }
                    }
                }
            }
            if side.contains(&b) {
                continue;
            }
            let flipped = orient & (1 << e) != 0;
            for (end, first) in [(a, !flipped), (b, flipped)] {
                let part: Vec<u8> = (0..peptides.len() as u8).filter(|v| side.contains(v) == side.contains(&end)).collect();
                if part.len() < 2 {
                    // a single freed peptide: its doublet is already counted
                    continue;
                }
                let inner = rest.iter().filter(|(u, v)| part.contains(u) && part.contains(v)).count();
                let mass = part.iter().map(|&v| masses[v as usize]).sum::<f32>() + linker * inner as f32;
                let arm = &arms[if asymmetric && !first { 1 } else { 0 }];
                topology_peaks += arm.all().into_iter().filter(|&s| has_peak(peaks, tol, mass + s)).count() as u16;
            }
        }

        let mut links: Vec<(u8, u16, u8, u16)> = edges.iter().map(|&(a, b)| (a, 0, b, 0)).collect();
        for (x, (fit, ends)) in fits.iter().zip(&end_edges).enumerate() {
            for (&site, &e) in fit.sites.iter().zip(ends) {
                if links[e].0 == x as u8 {
                    links[e].1 = site as u16;
                } else {
                    links[e].3 = site as u16;
                }
            }
        }
        let total = fits.iter().fold(Tally::default(), |acc, f| acc.combine(&f.tally));
        let decoys = peptides.iter().filter(|&&ix| db[ix].decoy).count() as u8;
        let class = match decoys {
            0 => PairClass::TT,
            d if d as usize == peptides.len() => PairClass::DD,
            _ => PairClass::TD,
        };
        let intra = edges.iter().all(|&(a, b)| {
            let (pa, pb) = (&db[peptides[a as usize]], &db[peptides[b as usize]]);
            pa.proteins.iter().any(|p| pb.proteins.contains(p))
        });
        let (a, sa, b, sb) = links[0];
        let score = |f: &CleavedFit| f.tally.score(self.scorer.score_type);
        Some(CrosslinkMatch {
            psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
            spec_id: query.id.clone(),
            file_id: query.file_id,
            kind: LinkKind::Multi,
            class,
            alpha: peptides[a as usize],
            beta: Some(peptides[b as usize]),
            alpha_site: sa,
            beta_site: sb,
            intra_protein: intra,
            hyperscore: total.score(self.scorer.score_type),
            alpha_score: score(&fits[a as usize]),
            beta_score: score(&fits[b as usize]),
            alpha_matched: fits[a as usize].tally.matched(),
            beta_matched: fits[b as usize].tally.matched(),
            alpha_linked: fits[a as usize].linked,
            beta_linked: fits[b as usize].linked,
            stub_peaks: fits.iter().map(|f| f.doublet).sum::<u16>().min(u8::MAX as u16) as u8,
            species: Some(Box::new(Species {
                peptides: peptides.to_vec(),
                links,
                matched: fits.iter().map(|f| f.tally.matched()).collect(),
                linked: fits.iter().map(|f| f.linked).collect(),
                decoys,
                topology_peaks,
            })),
            ..Default::default()
        })
    }

    /// Multivalent species for one spectrum: `(matches, species scored)`
    pub(super) fn score_multivalent(&self, query: &ProcessedSpectrum) -> (Vec<CrosslinkMatch>, u32) {
        let sc = self.scorer;
        let cfg = self.settings;
        let db = sc.db;
        let Some(valency) = cfg.multivalent else { return (Vec::new(), 0) };
        let Some(precursor) = query.precursors.first() else { return (Vec::new(), 0) };
        if db.peptides.is_empty() {
            return (Vec::new(), 0);
        }
        let min_peptide_mass = db.peptides[0].monoisotopic;
        let linker = cfg.linker_mass();
        let hetero = cfg.is_heterobifunctional();
        // MS-cleavable linkers: scored as broken, with their stubs
        let stubs = cfg.stubs();
        let asymmetric = stubs.as_ref().is_some_and(|s| s.other.is_some());
        let mz = precursor.mz - PROTON;
        let charges: Vec<u8> = match precursor.charge {
            Some(z) if !sc.override_precursor_charge => vec![z],
            _ => (sc.min_precursor_charge..=sc.max_precursor_charge).collect(),
        };

        let mut found: Vec<CrosslinkMatch> = Vec::new();
        let mut scored = 0u32;
        for charge in charges {
            let mass = mz * charge as f32;
            let max_fc = frag_charge_bound(sc.max_fragment_charge, charge);
            let largest = sc.precursor_tol.bounds(mass - sc.min_isotope_err as f32 * NEUTRON).1;
            // a seed is at most the precursor less two linkers and one more peptide
            let seed_limit = largest - 2.0 * linker - min_peptide_mass;
            if seed_limit < min_peptide_mass {
                continue;
            }
            let seeds: Vec<PeptideIx> = self
                .seeds(query, seed_limit, max_fc)
                .into_iter()
                .filter(|&ix| !cfg.link_sites(&db[ix]).is_empty())
                .take(MULTI_SEEDS)
                .collect();
            if seeds.is_empty() {
                continue;
            }
            // (peptide, groups and roles) -> best residues
            let mut cache: HashMap<(u32, Vec<(u32, Vec<u8>)>), Option<ChainFit>> = HashMap::new();
            // cleavable: (peptide, arms and roles of its ends) -> best residues
            let mut cleaved: HashMap<(u32, Vec<(usize, u8)>), Option<CleavedFit>> = HashMap::new();
            // peptides with stubs are charge-reduced: neutral masses up to charge - 1
            let peaks = match stubs {
                Some(_) => neutral_candidates(query, charge.saturating_sub(1).max(1)),
                None => Vec::new(),
            };
            // each species (peptides as a sorted list, its links) once per charge
            let mut seen: HashSet<(Vec<u32>, Vec<(u8, u8)>)> = HashSet::new();

            for iso in sc.min_isotope_err..=sc.max_isotope_err {
                let m = mass - iso as f32 * NEUTRON;
                let (lo, hi) = sc.precursor_tol.bounds(m);
                for k in 2..=valency as usize {
                    for p in 2..=k + 1 {
                        let shapes = arrangements(p, k);
                        // multisets of p - 1 seeds
                        let mut pick = vec![0usize; p - 1];
                        loop {
                            let chosen: Vec<PeptideIx> = pick.iter().map(|&i| seeds[i]).collect();
                            let known: f32 = chosen.iter().map(|&ix| db[ix].monoisotopic).sum();
                            let (rest_lo, rest_hi) = (lo - known - k as f32 * linker, hi - known - k as f32 * linker);
                            if rest_hi >= min_peptide_mass {
                                let start = db.peptides.partition_point(|q| q.monoisotopic < rest_lo);
                                for (offset, last) in db.peptides[start..]
                                    .iter()
                                    .take_while(|q| q.monoisotopic <= rest_hi)
                                    .take(MAX_PARTNERS_PER_SEED)
                                    .enumerate()
                                {
                                    if cfg.link_sites(last).is_empty() {
                                        continue;
                                    }
                                    let mut peptides = chosen.clone();
                                    peptides.push(PeptideIx((start + offset) as u32));
                                    let masses: Vec<f32> = peptides.iter().map(|&ix| db[ix].monoisotopic).collect();
                                    let mut key: Vec<u32> = peptides.iter().map(|ix| ix.0).collect();
                                    key.sort_unstable();
                                    for edges in &shapes {
                                        if !seen.insert((key.clone(), edges.clone())) {
                                            continue;
                                        }
                                        scored += 1;
                                        // each link both ways round when its ends differ
                                        let orientations = if hetero || asymmetric { 1u32 << k } else { 1 };
                                        for orient in 0..orientations {
                                            let csm = match &stubs {
                                                None => self.score_species(
                                                    query, &peptides, &masses, edges, orient, hetero, max_fc, &mut cache,
                                                ),
                                                Some(stubs) => self.score_species_cleaved(
                                                    query, &peptides, &masses, edges, orient, stubs, &peaks, max_fc,
                                                    &mut cleaved,
                                                ),
                                            };
                                            if let Some(csm) = csm {
                                                let calc = masses.iter().sum::<f32>() + k as f32 * linker;
                                                found.push(CrosslinkMatch {
                                                    charge,
                                                    rt: query.scan_start_time,
                                                    expmass: mass,
                                                    calcmass: calc,
                                                    precursor_ppm: ppm_delta(mass, iso, calc),
                                                    isotope_error: iso as f32 * NEUTRON,
                                                    ..csm
                                                });
                                            }
                                        }
                                    }
                                }
                            }
                            // next multiset (non-decreasing indices)
                            let Some(i) = (0..pick.len()).rev().find(|&i| pick[i] + 1 < seeds.len()) else { break };
                            pick[i] += 1;
                            let v = pick[i];
                            for j in i + 1..pick.len() {
                                pick[j] = v;
                            }
                        }
                    }
                }
            }
        }
        // best first; among equal scores (cleavable: arrangements of the same
        // peptides) the one with more partly cleaved pieces. Later sorts are
        // stable, so this order holds for ties.
        let support = |m: &CrosslinkMatch| m.species.as_ref().map_or(0, |s| s.topology_peaks);
        found.sort_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore).then(support(b).cmp(&support(a))));
        (found, scored)
    }

    /// Score one arrangement of `peptides`: `None` if a peptide lacks
    /// residues or evidence. Precursor fields are filled in by the caller.
    #[allow(clippy::too_many_arguments)]
    fn score_species(
        &self,
        query: &ProcessedSpectrum,
        peptides: &[PeptideIx],
        masses: &[f32],
        edges: &[(u8, u8)],
        orient: u32,
        hetero: bool,
        max_fc: u8,
        cache: &mut HashMap<(u32, Vec<(u32, Vec<u8>)>), Option<ChainFit>>,
    ) -> Option<CrosslinkMatch> {
        let db = self.scorer.db;
        let linker = self.settings.linker_mass();
        // role each end of each link must take: (first end, second end)
        let role = |e: usize| -> (u8, u8) {
            if !hetero {
                (ANY, ANY)
            } else if orient & (1 << e) == 0 {
                (SITE1, SITE2)
            } else {
                (SITE2, SITE1)
            }
        };
        let mut fits: Vec<ChainFit> = Vec::with_capacity(peptides.len());
        // per peptide: the edge index of each of its ends
        let mut end_edges: Vec<Vec<usize>> = Vec::with_capacity(peptides.len());
        for x in 0..peptides.len() as u8 {
            let groups = groups_of(x, edges, masses, linker);
            let mut ends: Vec<(usize, u8)> = Vec::new();
            let mut edge_of_end: Vec<usize> = Vec::new();
            for (g, (_, links)) in groups.iter().enumerate() {
                for &e in links {
                    let (a, _) = edges[e];
                    let r = if a == x { role(e).0 } else { role(e).1 };
                    ends.push((g, r));
                    edge_of_end.push(e);
                }
            }
            let shifts: Vec<f32> = groups.iter().map(|g| g.0).collect();
            let key = (
                peptides[x as usize].0,
                groups
                    .iter()
                    .enumerate()
                    .map(|(g, (mass, _))| (mass.to_bits(), ends.iter().filter(|e| e.0 == g).map(|e| e.1).collect()))
                    .collect(),
            );
            let fit = cache
                .entry(key)
                .or_insert_with(|| self.fit_chain(query, &db[peptides[x as usize]], &shifts, &ends, max_fc))
                .clone()?;
            fits.push(fit);
            end_edges.push(edge_of_end);
        }

        // residues of each link, from both ends
        let mut links: Vec<(u8, u16, u8, u16)> = edges.iter().map(|&(a, b)| (a, 0, b, 0)).collect();
        for (x, (fit, ends)) in fits.iter().zip(&end_edges).enumerate() {
            for (&site, &e) in fit.sites.iter().zip(ends) {
                if links[e].0 == x as u8 {
                    links[e].1 = site as u16;
                } else {
                    links[e].3 = site as u16;
                }
            }
        }
        let total = fits.iter().fold(Tally::default(), |acc, f| acc.combine(&f.tally));
        let score = total.score(self.scorer.score_type);
        let decoys = peptides.iter().filter(|&&ix| db[ix].decoy).count() as u8;
        let class = match decoys {
            0 => PairClass::TT,
            d if d as usize == peptides.len() => PairClass::DD,
            _ => PairClass::TD,
        };
        // every linked pair of peptides from one protein?
        let intra = edges.iter().all(|&(a, b)| {
            let (pa, pb) = (&db[peptides[a as usize]], &db[peptides[b as usize]]);
            pa.proteins.iter().any(|p| pb.proteins.contains(p))
        });
        let (a, sa, b, sb) = links[0];
        Some(CrosslinkMatch {
            psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
            spec_id: query.id.clone(),
            file_id: query.file_id,
            kind: LinkKind::Multi,
            class,
            alpha: peptides[a as usize],
            beta: Some(peptides[b as usize]),
            alpha_site: sa,
            beta_site: sb,
            intra_protein: intra,
            hyperscore: score,
            alpha_score: fits[a as usize].tally.score(self.scorer.score_type),
            beta_score: fits[b as usize].tally.score(self.scorer.score_type),
            alpha_matched: fits[a as usize].tally.matched(),
            beta_matched: fits[b as usize].tally.matched(),
            alpha_linked: fits[a as usize].linked,
            beta_linked: fits[b as usize].linked,
            species: Some(Box::new(Species {
                peptides: peptides.to_vec(),
                links,
                matched: fits.iter().map(|f| f.tally.matched()).collect(),
                linked: fits.iter().map(|f| f.linked).collect(),
                decoys,
                topology_peaks: 0,
            })),
            ..Default::default()
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn arrangements_are_connected_multigraphs() {
        // two peptides, two links: only a2b
        assert_eq!(arrangements(2, 2), vec![vec![(0, 1), (0, 1)]]);
        // three peptides, two links: the three chains
        assert_eq!(arrangements(3, 2).len(), 3);
        // four peptides, three links: the 16 trees on 4 labelled nodes
        assert_eq!(arrangements(4, 3).len(), 16);
        // five peptides, four links: Cayley, 5^3 = 125 trees
        assert_eq!(arrangements(5, 4).len(), 125);
        assert!(arrangements(3, 2).iter().all(|e| connected(3, e)));
    }

    #[test]
    fn groups_hang_off_a_peptide() {
        let masses = [1000.0, 2000.0, 3000.0];
        // chain a-b-c: removing b leaves a and c, one link each
        let g = groups_of(1, &[(0, 1), (1, 2)], &masses, 100.0);
        assert_eq!(g, vec![(1100.0, vec![0]), (3100.0, vec![1])]);
        // removing a leaves b-c as one group, with both linkers
        let g = groups_of(0, &[(0, 1), (1, 2)], &masses, 100.0);
        assert_eq!(g, vec![(5200.0, vec![0])]);
        // a2b: b hangs off a by two links
        let g = groups_of(0, &[(0, 1), (0, 1)], &masses[..2], 100.0);
        assert_eq!(g, vec![(2200.0, vec![0, 1])]);
    }
}
