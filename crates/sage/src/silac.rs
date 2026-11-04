//! SILAC: stable isotope labelling by amino acids in cell culture.
//!
//! Each configured channel gives residues (usually `K` and `R`) a heavy
//! isotope label. Every peptide is searched light and in each channel, with
//! *all* of its labelled residues carrying that channel's label: a peptide
//! comes from one population of cells, so mixed label states are not
//! generated. Labels are registered modifications (UNIMOD accessions for the
//! usual labels), so precursor and fragment masses, every output and every
//! viewer carry them; a label on a residue that already has a modification
//! (acetyl-K, a static TMT) is one combined modification. [`Labels`] knows
//! which modifications are labels, the channel of a peptide and its light
//! form.
//!
//! [`quantify`] traces each identified species (a peptide, or a crosslinked
//! species) in MS1 in every channel at once, with the isotope grids, smoothing
//! and time warping of LFQ, and integrates all channels over one shared peak,
//! so their ratios compare like with like. Species identified in one channel
//! are quantified in the others too.

use crate::lfq::{Grid, LfqSettings, Peak, GRID_SIZE, N_ISOTOPES, RT_TOL};
use crate::mass::{Tolerance, NEUTRON};
use crate::ml::retention_alignment::Alignment;
use crate::modification::{ModId, ModName, ModificationError, ModificationSpecificity, Modifications};
use crate::peptide::Peptide;
use crate::spectrum::ProcessedSpectrum;
use dashmap::DashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// One labelled channel (light, unlabelled, is implicit)
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SilacChannel {
    /// `heavy`, `medium`
    pub name: String,
    /// Residue (one letter) → label mass: `{"K": 8.014199, "R": 10.008269}`
    pub labels: BTreeMap<String, f32>,
}

/// `database.silac`: the labelled channels
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct SilacSettings {
    pub channels: Vec<SilacChannel>,
    /// Crosslinks between peptides of different channels (populations mixed
    /// before crosslinking). Off: a crosslink's peptides share a channel.
    #[serde(default)]
    pub mixed_crosslinks: bool,
}

/// At most this many labelled channels (plus light)
pub const MAX_CHANNELS: usize = 4;

impl SilacSettings {
    pub fn validate(&self) -> Result<(), String> {
        if self.channels.is_empty() {
            return Err("silac: give at least one labelled channel".into());
        }
        if self.channels.len() > MAX_CHANNELS {
            return Err(format!("silac: at most {MAX_CHANNELS} labelled channels"));
        }
        let mut names = vec!["light".to_string()];
        let mut seen: Vec<&BTreeMap<String, f32>> = Vec::new();
        for c in &self.channels {
            let name = c.name.trim().to_ascii_lowercase();
            if name.is_empty() || names.contains(&name) {
                return Err(format!("silac: channel names must be distinct and not \"light\" ({:?})", c.name));
            }
            names.push(name);
            if c.labels.is_empty() {
                return Err(format!("silac: channel {} labels no residue", c.name));
            }
            for (residue, &mass) in &c.labels {
                let ok = residue.len() == 1 && residue.bytes().all(|b| crate::mass::VALID_AA.contains(&b));
                if !ok {
                    return Err(format!("silac: {residue:?} is not a residue (one letter, e.g. \"K\")"));
                }
                if !mass.is_finite() || mass <= 0.0 {
                    return Err(format!("silac: label mass of {residue} in {} must be positive", c.name));
                }
            }
            if seen.contains(&&c.labels) {
                return Err(format!("silac: channel {} has the same labels as another", c.name));
            }
            seen.push(&c.labels);
        }
        Ok(())
    }

    /// Channel names, light first
    pub fn names(&self) -> Vec<String> {
        std::iter::once("light".to_string()).chain(self.channels.iter().map(|c| c.name.clone())).collect()
    }
}

/// UNIMOD accessions of the usual labels, which ProForma readers (Annotator)
/// know
fn accession(mass: f32) -> Option<&'static str> {
    const KNOWN: [(f32, &str); 4] = [
        // Label:13C(6)15N(2): Lys8
        (8.014199, "UNIMOD:259"),
        // Label:13C(6)15N(4): Arg10
        (10.008269, "UNIMOD:267"),
        // Label:13C(6): Lys6, Arg6
        (6.020129, "UNIMOD:188"),
        // Label:2H(4): Lys4
        (4.025107, "UNIMOD:481"),
    ];
    KNOWN.iter().find(|(m, _)| (m - mass).abs() < 1e-3).map(|(_, n)| *n)
}

/// Can a modification of `specificity` sit on residue `residue`?
fn on_residue(specificity: ModificationSpecificity, residue: u8) -> bool {
    use ModificationSpecificity::*;
    match specificity {
        Residue(r) => r == residue,
        PeptideN(Some(r)) | PeptideC(Some(r)) | ProteinN(Some(r)) | ProteinC(Some(r)) => r == residue,
        _ => false,
    }
}

/// The labels of a database, as modification ids
#[derive(Clone, Debug, Default)]
pub struct Labels {
    /// `light`, then each channel's name
    pub names: Vec<String>,
    /// Per labelled channel (channel `c` at `c - 1`): residue → (label, mass)
    channels: Vec<HashMap<u8, (ModId, f32)>>,
    /// (channel, residue, modification already there) → both together
    combined: HashMap<(u8, u8, ModId), ModId>,
    /// (modification, residue) → (channel, the modification without the label)
    reverse: HashMap<(ModId, u8), (u8, Option<ModId>)>,
    /// Residues labelled in some channel
    residues: Vec<u8>,
    pub mixed_crosslinks: bool,
}

impl Labels {
    /// Register the labels of `settings` in `registry`, and each configured
    /// modification (`configured`) of a labelled residue combined with its label
    pub fn register(
        settings: &SilacSettings,
        registry: &mut Modifications,
        configured: &[(ModificationSpecificity, ModId)],
    ) -> Result<Labels, ModificationError> {
        let mut labels = Labels {
            names: settings.names(),
            mixed_crosslinks: settings.mixed_crosslinks,
            ..Default::default()
        };
        for (k, channel) in settings.channels.iter().enumerate() {
            let c = k as u8 + 1;
            let mut map = HashMap::new();
            for (residue, &mass) in &channel.labels {
                let r = residue.as_bytes()[0];
                let name: Option<ModName> = accession(mass).and_then(|a| a.parse().ok());
                let id = registry.insert(mass, name.as_ref())?;
                map.insert(r, (id, mass));
                labels.reverse.insert((id, r), (c, None));
                if !labels.residues.contains(&r) {
                    labels.residues.push(r);
                }
                for &(specificity, base) in configured {
                    if !on_residue(specificity, r) {
                        continue;
                    }
                    let both = registry.insert(registry.get(base).mass + mass, None)?;
                    labels.combined.insert((c, r, base), both);
                    labels.reverse.insert((both, r), (c, Some(base)));
                }
            }
            labels.channels.push(map);
        }
        Ok(labels)
    }

    /// Light channel plus labelled ones
    pub fn len(&self) -> usize {
        self.channels.len() + 1
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// The label modifications themselves (not combined with another)
    pub fn label_ids(&self) -> Vec<ModId> {
        let mut ids: Vec<ModId> = self.channels.iter().flat_map(|m| m.values().map(|v| v.0)).collect();
        ids.sort();
        ids.dedup();
        ids
    }

    /// `peptide` (light) and the peptide in each channel that labels any of
    /// its residues, every labelled residue carrying the label
    pub fn states(&self, peptide: Peptide) -> Vec<Peptide> {
        let mut out = Vec::with_capacity(self.len());
        for (k, labels) in self.channels.iter().enumerate() {
            let c = k as u8 + 1;
            let mut p = peptide.clone();
            let (mut any, mut ok) = (false, true);
            for i in 0..p.sequence.len() {
                let aa = p.sequence[i];
                let Some(&(id, mass)) = labels.get(&aa) else { continue };
                let labelled = match p.modifications[i] {
                    None => Some(id),
                    Some(m) => self.combined.get(&(c, aa, m)).copied(),
                };
                match labelled {
                    Some(x) => {
                        p.modifications[i] = Some(x);
                        p.monoisotopic += mass;
                        any = true;
                    }
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if any && ok {
                out.push(p);
            }
        }
        out.insert(0, peptide);
        out
    }

    /// The channel of a peptide: `0` light, `None` if it has no residue that
    /// any channel labels (all channels look the same)
    pub fn channel(&self, peptide: &Peptide) -> Option<u8> {
        let mut labellable = false;
        for (i, aa) in peptide.sequence.iter().enumerate() {
            if !self.residues.contains(aa) {
                continue;
            }
            labellable = true;
            if let Some(m) = peptide.modifications[i] {
                if let Some(&(c, _)) = self.reverse.get(&(m, *aa)) {
                    return Some(c);
                }
            }
        }
        labellable.then_some(0)
    }

    /// The light form of a peptide (labels removed)
    pub fn unlabel(&self, peptide: &Peptide) -> Peptide {
        let mut p = peptide.clone();
        for i in 0..p.sequence.len() {
            let aa = p.sequence[i];
            let Some(m) = p.modifications[i] else { continue };
            if let Some(&(c, base)) = self.reverse.get(&(m, aa)) {
                p.modifications[i] = base;
                p.monoisotopic -= self.channels[c as usize - 1][&aa].1;
            }
        }
        p
    }

    /// Mass added to the light form of `peptide` by channel `channel`
    pub fn shift(&self, peptide: &Peptide, channel: u8) -> f32 {
        if channel == 0 {
            return 0.0;
        }
        let labels = &self.channels[channel as usize - 1];
        peptide.sequence.iter().filter_map(|aa| labels.get(aa).map(|l| l.1)).sum()
    }
}

/// One identification to quantify
#[derive(Clone, Debug)]
pub struct Identified {
    /// The species across channels: e.g. its light form
    pub key: String,
    pub channel: u8,
    pub file_id: usize,
    pub rt: f32,
    /// Precursor inverse ion mobility, 0 without
    pub ims: f32,
    pub score: f32,
    /// Neutral mass of the species in each channel, light first
    pub masses: Vec<f32>,
    /// Theoretical isotope envelope
    pub envelope: [f32; N_ISOTOPES],
}

/// One quantified species
#[derive(Clone, Debug)]
pub struct SilacQuant {
    /// Its identifications (indices into the input), best first
    pub identifications: Vec<usize>,
    /// The shared MS1 peak: score and spectral angle
    pub peak: Peak,
    /// Integrated intensity, `areas[run][channel]`
    pub areas: Vec<Vec<f64>>,
}

/// Runs need at least this many shared species for a retention-time fit
const MIN_ALIGNMENT_POINTS: usize = 3;

fn aligned(a: &Alignment, rt: f32) -> f32 {
    (rt / a.max_rt) * a.slope + a.intercept
}

/// Linear retention-time alignment of each run onto the average run, fitted
/// on species identified in at least two runs (as crosslink LFQ)
fn align(ids: &[Identified], species: &[Vec<usize>], ms1: &[ProcessedSpectrum], n_files: usize) -> Vec<Alignment> {
    let mut max_rt = vec![0.0f32; n_files];
    for s in ms1 {
        if let Some(x) = max_rt.get_mut(s.file_id) {
            *x = x.max(s.scan_start_time);
        }
    }
    for id in ids {
        if let Some(x) = max_rt.get_mut(id.file_id) {
            *x = x.max(id.rt);
        }
    }
    for x in &mut max_rt {
        if *x <= 0.0 {
            *x = 1.0;
        }
    }
    let mut points: Vec<Vec<(f32, f32)>> = vec![Vec::new(); n_files];
    for members in species {
        let mut per_file: HashMap<usize, f32> = HashMap::new();
        for &i in members {
            per_file.entry(ids[i].file_id).or_insert(ids[i].rt / max_rt[ids[i].file_id]);
        }
        if per_file.len() < 2 {
            continue;
        }
        let mean = per_file.values().sum::<f32>() / per_file.len() as f32;
        for (&file, &rt) in &per_file {
            points[file].push((rt, mean));
        }
    }
    (0..n_files)
        .map(|file_id| {
            let identity = Alignment { file_id, max_rt: max_rt[file_id], slope: 1.0, intercept: 0.0 };
            let pts = &points[file_id];
            if pts.len() < MIN_ALIGNMENT_POINTS {
                return identity;
            }
            let n = pts.len() as f32;
            let mx = pts.iter().map(|p| p.0).sum::<f32>() / n;
            let my = pts.iter().map(|p| p.1).sum::<f32>() / n;
            let sxx = pts.iter().map(|p| (p.0 - mx).powi(2)).sum::<f32>();
            let sxy = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f32>();
            if sxx <= f32::EPSILON {
                return identity;
            }
            let slope = sxy / sxx;
            if !(0.5..=2.0).contains(&slope) {
                return identity;
            }
            Alignment { slope, intercept: my - slope * mx, ..identity }
        })
        .collect()
}

/// One isotope of one charge state of one channel of a species
#[derive(Copy, Clone)]
struct Range {
    mz_lo: f32,
    mz_hi: f32,
    rt: f32,
    mobility_lo: f32,
    mobility_hi: f32,
    species: usize,
    channel: usize,
    isotope: usize,
}

/// Quantify every identified species in all channels. Each species is traced
/// at every precursor charge and the first isotopes of each channel's mass,
/// around the (aligned) retention time of its best identification, in every
/// run; all channels and runs are integrated over one peak.
pub fn quantify(
    ids: &[Identified],
    lfq: &LfqSettings,
    ms1: &[ProcessedSpectrum],
    n_files: usize,
    precursor_charge: (u8, u8),
) -> Vec<SilacQuant> {
    if ids.is_empty() || ms1.is_empty() || n_files == 0 {
        return Vec::new();
    }
    let n_ch = ids[0].masses.len();
    let mut order: Vec<usize> = (0..ids.len()).filter(|&i| ids[i].masses.len() == n_ch).collect();
    order.sort_by(|&a, &b| ids[b].score.total_cmp(&ids[a].score));
    let mut index: HashMap<&str, usize> = HashMap::new();
    let mut species: Vec<Vec<usize>> = Vec::new();
    for i in order {
        let slot = *index.entry(ids[i].key.as_str()).or_insert_with(|| {
            species.push(Vec::new());
            species.len() - 1
        });
        species[slot].push(i);
    }
    let alignments = align(ids, &species, ms1, n_files);

    let ppm = Tolerance::Ppm(-lfq.ppm_tolerance, lfq.ppm_tolerance);
    let mut ranges = Vec::new();
    for (s, members) in species.iter().enumerate() {
        let reference = &ids[members[0]];
        let rt = aligned(&alignments[reference.file_id], reference.rt);
        let (mobility_lo, mobility_hi) = if reference.ims > 0.0 {
            Tolerance::Pct(-lfq.mobility_pct_tolerance, lfq.mobility_pct_tolerance).bounds(reference.ims)
        } else {
            (f32::MIN, f32::MAX)
        };
        for (channel, &mass) in reference.masses.iter().enumerate() {
            for charge in precursor_charge.0.max(1)..=precursor_charge.1.max(1) {
                for isotope in 0..N_ISOTOPES {
                    let mz = (mass + isotope as f32 * NEUTRON) / charge as f32;
                    let (mz_lo, mz_hi) = ppm.bounds(mz);
                    ranges.push(Range { mz_lo, mz_hi, rt, mobility_lo, mobility_hi, species: s, channel, isotope });
                }
            }
        }
    }
    ranges.sort_by(|a, b| a.mz_lo.total_cmp(&b.mz_lo));
    let max_width = ranges.iter().map(|r| r.mz_hi - r.mz_lo).fold(0.0, f32::max);

    log::info!("SILAC: tracing {} species in {} channels", species.len(), n_ch);
    // one grid per species; its "files" are (run, channel)
    let grids: DashMap<usize, Grid, fnv::FnvBuildHasher> = DashMap::default();
    ms1.par_iter().for_each(|spectrum| {
        let Some(a) = alignments.get(spectrum.file_id) else { return };
        let rt = aligned(a, spectrum.scan_start_time);
        for (k, (&mass, &intensity)) in spectrum.masses.iter().zip(spectrum.intensities.iter()).enumerate() {
            let mobility = spectrum.mobilities.get(k).copied();
            let end = ranges.partition_point(|r| r.mz_lo <= mass);
            for r in ranges[..end].iter().rev() {
                if r.mz_lo < mass - max_width {
                    break;
                }
                if mass > r.mz_hi
                    || (rt - r.rt).abs() > RT_TOL
                    || mobility.is_some_and(|x| x < r.mobility_lo || x > r.mobility_hi)
                {
                    continue;
                }
                let reference = &ids[species[r.species][0]];
                grids
                    .entry(r.species)
                    .or_insert_with(|| {
                        Grid::centered(
                            r.rt,
                            reference.file_id * n_ch + reference.channel as usize,
                            RT_TOL,
                            reference.envelope,
                            n_files * n_ch,
                            GRID_SIZE,
                        )
                    })
                    .add_entry(rt, r.isotope, spectrum.file_id * n_ch + r.channel, intensity);
            }
        }
    });

    let mut found: Vec<SilacQuant> = grids
        .into_par_iter()
        .filter_map(|(s, mut grid)| {
            let (peak, flat) = grid.summarize_traces().integrate(lfq)?;
            let areas = flat.chunks(n_ch).map(|c| c.to_vec()).collect();
            Some(SilacQuant { identifications: species[s].clone(), peak, areas })
        })
        .collect();
    found.sort_by(|a, b| b.peak.score.total_cmp(&a.peak.score));
    found
}

/// `a / b`, or `None` when either is missing
pub fn ratio(a: f64, b: f64) -> Option<f64> {
    (a > 0.0 && b > 0.0).then(|| a / b)
}

/// Median of the ratios (in log space, the same for an odd count)
pub fn median_ratio(ratios: &mut [f64]) -> Option<f64> {
    if ratios.is_empty() {
        return None;
    }
    ratios.sort_by(f64::total_cmp);
    let n = ratios.len();
    Some(if n % 2 == 1 {
        ratios[n / 2]
    } else {
        // geometric mean of the middle two (symmetric for H/L and L/H)
        (ratios[n / 2 - 1] * ratios[n / 2]).sqrt()
    })
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::Builder;

    fn duplex() -> SilacSettings {
        SilacSettings {
            channels: vec![SilacChannel {
                name: "heavy".into(),
                labels: [("K".to_string(), 8.014199), ("R".to_string(), 10.008269)].into(),
            }],
            mixed_crosslinks: false,
        }
    }

    #[test]
    fn settings_are_checked() {
        assert!(duplex().validate().is_ok());
        assert!(SilacSettings::default().validate().is_err());
        let mut bad = duplex();
        bad.channels[0].name = "Light".into();
        assert!(bad.validate().is_err());
        let mut bad = duplex();
        bad.channels[0].labels.insert("KR".into(), 1.0);
        assert!(bad.validate().is_err());
        let mut twice = duplex();
        twice.channels.push(SilacChannel { name: "again".into(), ..duplex().channels[0].clone() });
        assert!(twice.validate().is_err(), "two channels with the same labels");
    }

    /// Peptides are generated light and heavy, every K/R labelled; a label on
    /// a modified residue is one combined modification
    #[test]
    fn peptides_in_label_states() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": {"C": 57.021465},
            "variable_mods": {"K": [42.010565]},
            "silac": {"channels": [{"name": "heavy", "labels": {"K": 8.014199, "R": 10.008269}}]}
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let labels = params.labels.clone().unwrap();
        assert_eq!(labels.names, ["light", "heavy"]);
        let fasta = crate::fasta::Fasta::parse(">p\nMACDKEFGHRLMNPQRSTWYAG\n".into(), "rev_", false);
        let peptides = params.digest(&fasta);
        let show = |p: &Peptide| p.display(&params.modifications.registry).to_string();
        let names: Vec<String> = peptides.iter().map(show).collect();
        // ACDK (light, acetyl, heavy, heavy + acetyl) ...
        let light = peptides.iter().find(|p| show(p) == "EFGHR").expect("EFGHR");
        let heavy = peptides.iter().find(|p| show(p) == "EFGHR[UNIMOD:267]").expect("EFGHR heavy");
        assert!((heavy.monoisotopic - light.monoisotopic - 10.008269).abs() < 1e-3);
        assert_eq!(labels.channel(light), Some(0));
        assert_eq!(labels.channel(heavy), Some(1));
        assert_eq!(show(&labels.unlabel(heavy)), "EFGHR");
        assert!((labels.shift(light, 1) - 10.008269).abs() < 1e-4);
        // the acetylated lysine, heavy: one combined modification of both masses
        let acetyl_heavy = peptides
            .iter()
            .find(|p| p.sequence.as_ref() == b"MACDK" && labels.channel(p) == Some(1) && p.modifications[4].is_some_and(|m| {
                (params.modifications.registry.get(m).mass - 50.024764).abs() < 1e-3
            }));
        assert!(acetyl_heavy.is_some(), "{names:?}");
        let acetyl_heavy = acetyl_heavy.unwrap();
        let back = labels.unlabel(acetyl_heavy);
        assert!((params.modifications.registry.mass(back.modifications[4]) - 42.010565).abs() < 1e-4);
        // the protein C-terminal peptide has no K or R: no channel, no heavy copy
        let plain: Vec<&Peptide> = peptides.iter().filter(|p| p.sequence.as_ref() == b"STWYAG").collect();
        assert_eq!(plain.len(), 1, "{names:?}");
        assert_eq!(labels.channel(plain[0]), None);
        assert_eq!(labels.label_ids().len(), 2);
    }

    /// Light and heavy traces at 1 : 2 in MS1 give a heavy / light ratio of 2
    #[test]
    fn ratios_from_ms1() {
        let light_mass = 1500.0f32;
        let heavy_mass = light_mass + 8.014199;
        let charge = 2u8;
        let mut ms1 = Vec::new();
        for scan in 0..200 {
            let rt = scan as f32 * 0.05; // 0-10 min
            let profile = (-((rt - 5.0) / 0.15).powi(2)).exp();
            let mut masses = Vec::new();
            let mut intensities = Vec::new();
            for (mass, scale) in [(light_mass, 1000.0f32), (heavy_mass, 2000.0f32)] {
                for (iso, abundance) in [1.0f32, 0.7, 0.3].iter().enumerate() {
                    // processed MS1 peaks are m/z less a proton
                    masses.push((mass + iso as f32 * NEUTRON) / charge as f32);
                    intensities.push(scale * abundance * profile + 1.0);
                }
            }
            let mut pairs: Vec<(f32, f32)> = masses.into_iter().zip(intensities).collect();
            pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
            ms1.push(ProcessedSpectrum {
                level: 1,
                id: scan.to_string(),
                file_id: 0,
                scan_start_time: rt,
                ion_injection_time: 0.0,
                precursors: Vec::new(),
                masses: pairs.iter().map(|p| p.0).collect(),
                intensities: pairs.iter().map(|p| p.1).collect(),
                charges: vec![1; pairs.len()],
                mobilities: Vec::new(),
                total_ion_current: 0.0,
            });
        }
        let id = Identified {
            key: "PEPTIDEK".into(),
            channel: 0,
            file_id: 0,
            rt: 5.0,
            ims: 0.0,
            score: 30.0,
            masses: vec![light_mass, heavy_mass],
            envelope: [1.0, 0.7, 0.3],
        };
        let lfq = LfqSettings { spectral_angle: 0.5, ..Default::default() };
        let found = quantify(&[id], &lfq, &ms1, 1, (charge, charge));
        assert_eq!(found.len(), 1);
        let a = &found[0].areas[0];
        let r = ratio(a[1], a[0]).unwrap();
        assert!((r - 2.0).abs() < 0.1, "heavy / light {r} ({a:?})");
        assert_eq!(median_ratio(&mut [4.0, 1.0]).unwrap(), 2.0);
    }
}
