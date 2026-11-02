//! Data-parallel chain scoring for non-cleavable crosslink searches
//! (the `reflect` and `gpu` backends).
//!
//! # The idea
//!
//! For a crosslinked pair `A-B` with precursor mass `M = m(A) + m(B) + linker`,
//! every fragment of `A` that contains the linked residue carries the whole
//! other chain plus the linker: a shift of `M - m(A)`. That shift depends only
//! on `A` and the precursor, **not on B**. So each peptide can be scored on its
//! own against the spectrum, with its linked-residue ions shifted by
//! `M - m(A)`; the pair is then recovered by a join on mass (`m(B) = M -
//! linker - m(A)`). This replaces the quadratic pair enumeration with one
//! independent job per peptide, which is exactly the shape a GPU wants.
//!
//! The per-peptide job ([`score_peptide`]) is the *specification*. The CPU
//! kernel runs it directly; a GPU kernel must reproduce it.
//!
//! # Precision note
//!
//! The shift uses the *observed* precursor mass, so shifted ions carry the
//! precursor error (e.g. 10 ppm of 3 kDa = 0.03 Da). Shifted-ion windows are
//! therefore widened by that error (`slack`). The result is a coarse,
//! surrogate score used only to choose candidates; the final score is always
//! recomputed exactly with the real partner mass.

use crate::database::IndexedDatabase;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::peptide::Peptide;
use crate::spectrum::ProcessedSpectrum;

/// Peptides with more candidate link sites than this only use the first ones
pub const MAX_KERNEL_SITES: usize = 16;

/// Per-peptide ion data, flattened for upload to a device
pub struct IonTable {
    /// Neutral monoisotopic mass of every peptide in the database
    pub peptide_mass: Vec<f32>,
    /// `ion_off[p]..ion_off[p+1]` are peptide `p`'s ions (empty if it has no link site)
    pub ion_off: Vec<u32>,
    /// Neutral fragment masses (unshifted)
    pub ion_mass: Vec<f32>,
    /// bit 0: N-terminal series (a/b/c); bits 1..: ion index within its series
    pub ion_meta: Vec<u32>,
    /// `site_off[p]..site_off[p+1]` index into `sites`
    pub site_off: Vec<u32>,
    /// 0-based link residues
    pub sites: Vec<u32>,
}

impl IonTable {
    /// `link_sites` decides which residues can carry a link
    pub fn build<F>(db: &IndexedDatabase, link_sites: F) -> Self
    where
        F: Fn(&Peptide) -> Vec<usize>,
    {
        let n = db.peptides.len();
        let mut t = IonTable {
            peptide_mass: Vec::with_capacity(n),
            ion_off: Vec::with_capacity(n + 1),
            ion_mass: Vec::new(),
            ion_meta: Vec::new(),
            site_off: Vec::with_capacity(n + 1),
            sites: Vec::new(),
        };
        t.ion_off.push(0);
        t.site_off.push(0);
        for peptide in &db.peptides {
            t.peptide_mass.push(peptide.monoisotopic);
            let sites = link_sites(peptide);
            if !sites.is_empty() {
                for &kind in &db.ion_kinds {
                    let n_term = matches!(kind, Kind::A | Kind::B | Kind::C) as u32;
                    for (idx, ion) in IonSeries::new(peptide, kind, &db.modifications).enumerate() {
                        t.ion_mass.push(ion.monoisotopic_mass);
                        t.ion_meta.push(n_term | ((idx as u32) << 1));
                    }
                }
                t.sites
                    .extend(sites.iter().take(MAX_KERNEL_SITES).map(|&s| s as u32));
            }
            t.ion_off.push(t.ion_mass.len() as u32);
            t.site_off.push(t.sites.len() as u32);
        }
        t
    }

    pub fn len(&self) -> usize {
        self.peptide_mass.len()
    }

    pub fn is_empty(&self) -> bool {
        self.peptide_mass.is_empty()
    }

    /// Approximate size on a device, in bytes
    pub fn bytes(&self) -> usize {
        4 * (self.peptide_mass.len()
            + self.ion_off.len()
            + self.ion_mass.len()
            + self.ion_meta.len()
            + self.site_off.len()
            + self.sites.len())
    }
}

/// A tolerance as `center * rel + abs`, so that a kernel needs no enum
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LinearTol {
    pub rel_lo: f32,
    pub rel_hi: f32,
    pub abs_lo: f32,
    pub abs_hi: f32,
}

impl From<Tolerance> for LinearTol {
    fn from(t: Tolerance) -> Self {
        match t {
            Tolerance::Ppm(lo, hi) => Self {
                rel_lo: lo / 1_000_000.0,
                rel_hi: hi / 1_000_000.0,
                abs_lo: 0.0,
                abs_hi: 0.0,
            },
            Tolerance::Pct(lo, hi) => Self {
                rel_lo: lo / 100.0,
                rel_hi: hi / 100.0,
                abs_lo: 0.0,
                abs_hi: 0.0,
            },
            Tolerance::Da(lo, hi) => Self {
                rel_lo: 0.0,
                rel_hi: 0.0,
                abs_lo: lo,
                abs_hi: hi,
            },
        }
    }
}

impl LinearTol {
    /// Window multipliers `1 + rel`, rounded once to f32
    #[inline]
    pub fn factors(&self) -> (f32, f32) {
        ((1.0 + self.rel_lo as f64) as f32, (1.0 + self.rel_hi as f64) as f32)
    }

    /// `[center * (1 + rel_lo) + abs_lo, center * (1 + rel_hi) + abs_hi]`.
    ///
    /// Every kernel must compute window edges exactly like this so that a peak
    /// on an edge is in or out on every device. The form is immune to fused
    /// multiply-add, which GPU drivers apply at will: a tolerance is either
    /// relative (`abs` = 0, so the product is rounded once either way) or
    /// absolute (factor = 1, so the product is exact).
    #[inline]
    pub fn bounds(&self, center: f32) -> (f32, f32) {
        let (lo, hi) = self.factors();
        (center * lo + self.abs_lo, center * hi + self.abs_hi)
    }
}

/// Everything a kernel needs to know about one (spectrum, precursor hypothesis)
#[derive(Copy, Clone, Debug)]
pub struct KernelParams {
    /// Neutral precursor mass `M` of this hypothesis
    pub precursor_mass: f32,
    /// Fragment charges tried are `1..max_fragment_charge` (exclusive, as in Sage)
    pub max_fragment_charge: u8,
    pub fragment_tol: LinearTol,
    /// Extra half-window (Da, neutral mass) for shifted ions: the precursor error
    pub slack: f32,
    /// Keep only peptides whose best site matches at least this many ions
    pub min_matched: u16,
    /// Keep at most this many survivors: the highest-scoring (ties: lowest index)
    pub max_survivors: usize,
}

impl KernelParams {
    /// `(1/z, slack/z)` for each fragment charge `z` in `1..max_fragment_charge`.
    ///
    /// Rounded once here and shared by every kernel: GPU division is not
    /// correctly rounded (WGSL allows 2.5 ulp), multiplication is.
    pub fn charge_terms(&self) -> impl Iterator<Item = (f32, f32)> + '_ {
        (1..self.max_fragment_charge).map(|z| {
            let z = z as f64;
            ((1.0 / z) as f32, (self.slack as f64 / z) as f32)
        })
    }
}

/// A peptide that passed the chain filter
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Survivor {
    pub peptide: u32,
    /// Surrogate hyperscore of the best link site
    pub score: f32,
    pub site: u32,
    pub matched: u16,
}

/// Something that can score every peptide of a database against a spectrum.
///
/// `survivors` returns the peptides whose best site matches at least
/// `params.min_matched` ions, limited to the best `params.max_survivors` by
/// score (see [`top_survivors`]), in ascending peptide order.
pub trait ChainKernel: Send + Sync {
    fn name(&self) -> &str;
    fn survivors(&self, query: &ProcessedSpectrum, params: &KernelParams) -> Vec<Survivor>;
}

/// Apply the `max_survivors` limit and put the result in peptide order
pub fn top_survivors(mut survivors: Vec<Survivor>, max: usize) -> Vec<Survivor> {
    if survivors.len() > max {
        survivors.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.peptide.cmp(&b.peptide)));
        survivors.truncate(max);
    }
    survivors.sort_by_key(|s| s.peptide);
    survivors
}

/// ln(n!) table shared by every kernel so that they agree
pub fn lnfact_table() -> Vec<f32> {
    let mut table = vec![0f32; 256];
    let mut acc = 0f64;
    for (n, slot) in table.iter_mut().enumerate().skip(1) {
        acc += (n as f64).ln();
        *slot = acc as f32;
    }
    table
}

/// Surrogate hyperscore (the Sage hyperscore formula)
#[inline]
pub fn surrogate_score(lnfact: &[f32], b: u32, y: u32, sum_b: f32, sum_y: f32) -> f32 {
    ((sum_b + 1.0) * (sum_y + 1.0)).ln() + lnfact[b.min(255) as usize] + lnfact[y.min(255) as usize]
}

/// Most intense peak in `[lo, hi]`; ties go to the later peak (as in Sage)
#[inline]
fn best_peak(masses: &[f32], intensities: &[f32], lo: f32, hi: f32) -> Option<f32> {
    let mut i = masses.partition_point(|&m| m < lo);
    let mut best: Option<f32> = None;
    while i < masses.len() && masses[i] <= hi {
        if best.map_or(true, |b| intensities[i] >= b) {
            best = Some(intensities[i]);
        }
        i += 1;
    }
    best
}

/// The per-peptide chain job: the specification every kernel implements.
pub fn score_peptide(
    table: &IonTable,
    peptide: usize,
    masses: &[f32],
    intensities: &[f32],
    params: &KernelParams,
    lnfact: &[f32],
) -> Option<Survivor> {
    let (i0, i1) = (
        table.ion_off[peptide] as usize,
        table.ion_off[peptide + 1] as usize,
    );
    let (s0, s1) = (
        table.site_off[peptide] as usize,
        table.site_off[peptide + 1] as usize,
    );
    if i0 == i1 || s0 == s1 {
        return None;
    }
    let shift = params.precursor_mass - table.peptide_mass[peptide];
    if shift <= 0.0 {
        return None;
    }

    // Per ion: (hits, intensity) when unshifted and when shifted
    let mut tallies: Vec<(u32, f32, u32, f32)> = Vec::with_capacity(i1 - i0);
    for ion in i0..i1 {
        let mass = table.ion_mass[ion];
        let (mut un, mut ui, mut sn, mut si) = (0u32, 0f32, 0u32, 0f32);
        for (inv_z, slack_z) in params.charge_terms() {
            let (lo, hi) = params.fragment_tol.bounds(mass * inv_z);
            if let Some(i) = best_peak(masses, intensities, lo, hi) {
                un += 1;
                ui += i;
            }
            let (lo, hi) = params.fragment_tol.bounds((mass + shift) * inv_z);
            let (lo, hi) = (lo - slack_z, hi + slack_z);
            if let Some(i) = best_peak(masses, intensities, lo, hi) {
                sn += 1;
                si += i;
            }
        }
        tallies.push((un, ui, sn, si));
    }

    let mut best: Option<Survivor> = None;
    for &site in &table.sites[s0..s1] {
        let (mut b, mut y, mut sum_b, mut sum_y) = (0u32, 0u32, 0f32, 0f32);
        for (k, ion) in (i0..i1).enumerate() {
            let meta = table.ion_meta[ion];
            let n_term = meta & 1 == 1;
            let idx = meta >> 1;
            let contains = if n_term { idx >= site } else { idx < site };
            let (n, i) = if contains {
                (tallies[k].2, tallies[k].3)
            } else {
                (tallies[k].0, tallies[k].1)
            };
            if n_term {
                b += n;
                sum_b += i;
            } else {
                y += n;
                sum_y += i;
            }
        }
        let matched = (b + y).min(u16::MAX as u32) as u16;
        let score = surrogate_score(lnfact, b, y, sum_b, sum_y);
        if matched >= params.min_matched && best.map_or(true, |x| score > x.score) {
            best = Some(Survivor {
                peptide: peptide as u32,
                score,
                site,
                matched,
            });
        }
    }
    best
}

/// Reference kernel: runs [`score_peptide`] over every peptide, serially
/// (callers already parallelise across spectra).
pub struct CpuKernel {
    pub table: IonTable,
    lnfact: Vec<f32>,
}

impl CpuKernel {
    pub fn new(table: IonTable) -> Self {
        Self {
            table,
            lnfact: lnfact_table(),
        }
    }
}

impl ChainKernel for CpuKernel {
    fn name(&self) -> &str {
        "reflect (cpu)"
    }

    fn survivors(&self, query: &ProcessedSpectrum, params: &KernelParams) -> Vec<Survivor> {
        let all = (0..self.table.len())
            .filter_map(|p| {
                score_peptide(
                    &self.table,
                    p,
                    &query.masses,
                    &query.intensities,
                    params,
                    &self.lnfact,
                )
            })
            .collect();
        top_survivors(all, params.max_survivors)
    }
}
