//! Label-free (MS1) quantification of crosslinked peptides.
//!
//! Sage's linear LFQ ([`crate::lfq`]) is keyed by single peptides and aligned
//! with retention times learned from linear PSMs. A crosslink search has
//! neither, so this module quantifies *crosslinked species* instead:
//!
//! * A species is a peptide pair (or a loop-linked peptide). Link-site isomers
//!   of the same pair have identical MS1 signals, so they share one species.
//! * Confidently identified species (target-target, `crosslink_q` at or below
//!   `lfq_settings.peptide_q_value`) are traced in the MS1 spectra of every run,
//!   at their charge states and first isotopes, with the same grids, smoothing,
//!   time warping and integration as linear LFQ.
//! * Runs are aligned with a linear retention-time fit on species identified in
//!   several runs (identity when there are too few).
//! * Each species gets a decoy trace at a shifted m/z and earlier retention
//!   time; the MS1 peak q-value is `decoys / targets` over the peak scores.

use crate::crosslink::{CrosslinkMatch, CrosslinkSettings, LinkKind, PairClass};
use crate::database::IndexedDatabase;
use crate::lfq::{Grid, LfqSettings, Peak, GRID_SIZE, N_ISOTOPES, RT_TOL};
use crate::mass::{composition, Composition, Tolerance, NEUTRON};
use crate::ml::retention_alignment::Alignment;
use crate::spectrum::ProcessedSpectrum;
use dashmap::DashMap;
use rayon::prelude::*;
use std::collections::HashMap;

/// m/z offset of decoy traces, as in linear LFQ
const DECOY_MZ_SHIFT: f32 = 11.06;
/// Runs need at least this many shared species for a retention-time fit
const MIN_ALIGNMENT_POINTS: usize = 3;

/// MS1 quantification of one crosslinked species
#[derive(Clone, Debug)]
pub struct CrosslinkLfq {
    /// Index (into the searched matches) of the best match of this species,
    /// which names the peptides and link sites
    pub csm: usize,
    /// `None` when charge states are combined
    pub charge: Option<u8>,
    pub decoy: bool,
    /// MS1 peak: score, spectral angle, q-value
    pub peak: Peak,
    /// Integrated intensity per run (file id order)
    pub areas: Vec<f64>,
}

/// One isotope of one charge state of a species, as an m/z window
#[derive(Copy, Clone, Debug)]
struct Range {
    mz_lo: f32,
    mz_hi: f32,
    /// Aligned retention time
    rt: f32,
    mobility_lo: f32,
    mobility_hi: f32,
    species: usize,
    charge: u8,
    isotope: usize,
    decoy: bool,
}

/// Peptide pair identity, independent of link sites and alpha/beta order
fn species_key(m: &CrosslinkMatch) -> (u32, u32) {
    match (m.kind, m.beta) {
        (LinkKind::Inter, Some(b)) => (m.alpha.0.min(b.0), m.alpha.0.max(b.0)),
        _ => (m.alpha.0, u32::MAX),
    }
}

fn aligned(a: &Alignment, rt: f32) -> f32 {
    (rt / a.max_rt) * a.slope + a.intercept
}

/// Linear retention-time alignment of each run onto the average run, fitted on
/// species identified in at least two runs. `species` holds each species' match
/// indices, best first.
pub fn align_runs(
    matches: &[CrosslinkMatch],
    species: &[Vec<usize>],
    ms1: &[ProcessedSpectrum],
    n_files: usize,
) -> Vec<Alignment> {
    let mut max_rt = vec![0.0f32; n_files];
    for s in ms1 {
        if let Some(x) = max_rt.get_mut(s.file_id) {
            *x = x.max(s.scan_start_time);
        }
    }
    for m in matches {
        if let Some(x) = max_rt.get_mut(m.file_id) {
            *x = x.max(m.rt);
        }
    }
    for x in &mut max_rt {
        if *x <= 0.0 {
            *x = 1.0;
        }
    }

    // Normalized RT of each species in each run (from its best match there)
    let mut points: Vec<Vec<(f32, f32)>> = vec![Vec::new(); n_files];
    for csms in species {
        let mut per_file: HashMap<usize, f32> = HashMap::new();
        for &i in csms {
            let m = &matches[i];
            per_file
                .entry(m.file_id)
                .or_insert(m.rt / max_rt[m.file_id]);
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
            let identity = Alignment {
                file_id,
                max_rt: max_rt[file_id],
                slope: 1.0,
                intercept: 0.0,
            };
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
            // A wild fit from a handful of points is worse than none
            if !(0.5..=2.0).contains(&slope) {
                return identity;
            }
            Alignment {
                slope,
                intercept: my - slope * mx,
                ..identity
            }
        })
        .collect()
}

/// Quantify confidently identified crosslinked species in the MS1 spectra of
/// all runs. `matches` must already carry q-values. Returns targets and decoys.
#[allow(clippy::too_many_arguments)]
pub fn quantify_lfq(
    matches: &[CrosslinkMatch],
    settings: &CrosslinkSettings,
    lfq: &LfqSettings,
    db: &IndexedDatabase,
    ms1: &[ProcessedSpectrum],
    n_files: usize,
    precursor_charge: (u8, u8),
) -> Vec<CrosslinkLfq> {
    // Group confident target matches by species, best first
    let mut order = (0..matches.len())
        .filter(|&i| {
            let m = &matches[i];
            m.class == PairClass::TT && m.crosslink_q <= lfq.peptide_q_value
        })
        .collect::<Vec<_>>();
    order.sort_by(|&a, &b| matches[b].hyperscore.total_cmp(&matches[a].hyperscore));
    let mut index: HashMap<(u32, u32), usize> = HashMap::new();
    let mut species: Vec<Vec<usize>> = Vec::new();
    for i in order {
        let slot = *index.entry(species_key(&matches[i])).or_insert_with(|| {
            species.push(Vec::new());
            species.len() - 1
        });
        species[slot].push(i);
    }
    if species.is_empty() {
        log::info!("crosslink LFQ: no crosslinks pass the q-value cutoff");
        return Vec::new();
    }
    if ms1.is_empty() {
        log::warn!("crosslink LFQ: no MS1 spectra found for quantification");
        return Vec::new();
    }

    let alignments = align_runs(matches, &species, ms1, n_files);

    // Theoretical isotope envelope of each species: both peptides plus linker
    let (linker_c, linker_s) = settings.linker_composition();
    let envelopes = species
        .iter()
        .map(|csms| {
            let m = &matches[csms[0]];
            let partner = match (m.kind, m.beta) {
                (LinkKind::Inter, Some(b)) => &db[b].sequence[..],
                _ => &[],
            };
            let c = db[m.alpha]
                .sequence
                .iter()
                .chain(partner)
                .map(|r| composition(*r))
                .sum::<Composition>();
            crate::isotopes::peptide_isotopes(c.carbon + linker_c, c.sulfur + linker_s)
        })
        .collect::<Vec<_>>();

    let ppm = Tolerance::Ppm(-lfq.ppm_tolerance, lfq.ppm_tolerance);
    let mut ranges = Vec::new();
    for (s, csms) in species.iter().enumerate() {
        let m = &matches[csms[0]];
        let rt = aligned(&alignments[m.file_id], m.rt);
        let (mobility_lo, mobility_hi) = if m.ims > 0.0 {
            Tolerance::Pct(-lfq.mobility_pct_tolerance, lfq.mobility_pct_tolerance).bounds(m.ims)
        } else {
            (f32::MIN, f32::MAX)
        };
        for charge in precursor_charge.0..=precursor_charge.1 {
            for isotope in 0..N_ISOTOPES {
                let mz = (m.calcmass + isotope as f32 * NEUTRON) / charge as f32;
                let fwd = Range {
                    mz_lo: ppm.bounds(mz).0,
                    mz_hi: ppm.bounds(mz).1,
                    rt,
                    mobility_lo,
                    mobility_hi,
                    species: s,
                    charge,
                    isotope,
                    decoy: false,
                };
                let (mz_lo, mz_hi) = ppm.bounds(mz + DECOY_MZ_SHIFT);
                let rev = Range {
                    mz_lo,
                    mz_hi,
                    rt: (rt - RT_TOL * 2.0).max(0.0),
                    decoy: true,
                    ..fwd
                };
                ranges.push(fwd);
                ranges.push(rev);
            }
        }
    }
    ranges.sort_by(|a, b| a.mz_lo.total_cmp(&b.mz_lo));
    let max_width = ranges.iter().map(|r| r.mz_hi - r.mz_lo).fold(0.0, f32::max);

    log::info!("tracing MS1 features of {} crosslinked species", species.len());
    // Grid per (species, charge or 0 when combined, decoy)
    let grids: DashMap<(usize, u8, bool), Grid, fnv::FnvBuildHasher> = DashMap::default();
    ms1.par_iter().for_each(|spectrum| {
        let Some(a) = alignments.get(spectrum.file_id) else {
            return;
        };
        let rt = aligned(a, spectrum.scan_start_time);
        for (k, (&mass, &intensity)) in spectrum
            .masses
            .iter()
            .zip(spectrum.intensities.iter())
            .enumerate()
        {
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
                let charge = if lfq.combine_charge_states { 0 } else { r.charge };
                let reference = &matches[species[r.species][0]];
                grids
                    .entry((r.species, charge, r.decoy))
                    .or_insert_with(|| {
                        Grid::centered(
                            r.rt,
                            reference.file_id,
                            RT_TOL,
                            envelopes[r.species],
                            n_files,
                            GRID_SIZE,
                        )
                    })
                    .add_entry(rt, r.isotope, spectrum.file_id, intensity);
            }
        }
    });

    log::info!("integrating crosslink MS1 features");
    let mut found = grids
        .into_par_iter()
        .filter_map(|((s, charge, decoy), mut grid)| {
            let (peak, areas) = grid.summarize_traces().integrate(lfq)?;
            Some(CrosslinkLfq {
                csm: species[s][0],
                charge: (charge > 0).then_some(charge),
                decoy,
                peak,
                areas,
            })
        })
        .collect::<Vec<_>>();

    assign_peak_q_values(&mut found);
    found.sort_by(|a, b| b.peak.score.total_cmp(&a.peak.score));
    found
}

/// Target-decoy q-values over MS1 peak scores (`decoys / targets`, as in
/// linear LFQ), made monotone
fn assign_peak_q_values(found: &mut [CrosslinkLfq]) {
    let mut order = (0..found.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| found[b].peak.score.total_cmp(&found[a].peak.score));
    let (mut decoys, mut targets) = (1.0f32, 0.0f32);
    let mut q = vec![1.0f32; found.len()];
    for (rank, &i) in order.iter().enumerate() {
        match found[i].decoy {
            true => decoys += 1.0,
            false => targets += 1.0,
        }
        q[rank] = (decoys / targets.max(1.0)).min(1.0);
    }
    let mut q_min = 1.0f32;
    for rank in (0..order.len()).rev() {
        q_min = q_min.min(q[rank]);
        found[order[rank]].peak.q_value = q_min;
    }
}
