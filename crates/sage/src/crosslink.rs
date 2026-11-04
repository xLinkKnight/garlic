//! Chemical crosslinking (XL-MS) search for **non-cleavable** reagents
//! (DSS, BS3, DSG, BS2G, ...).
//!
//! # Approach
//!
//! A crosslinked spectrum-match (CSM) explains one MS2 spectrum with *two*
//! peptides joined by a linker (or one peptide joined to itself, a "loop-link").
//! Enumerating all peptide pairs is quadratic in database size, so the search
//! is split in two stages, both of which reuse Sage's existing fragment index:
//!
//! 1. **Seed search.** The spectrum is searched against the fragment index with
//!    a very wide precursor window (every peptide that could be one half of the
//!    pair). Each chain of a crosslinked peptide still produces many ordinary,
//!    unshifted b/y ions, so the linear peptides that explain the most peaks
//!    are good seeds. The top `candidates` seeds are kept.
//! 2. **Partner lookup.** For each seed `A`, the partner mass is fixed by the
//!    precursor: `m(B) = M - m(A) - m(linker)`. Peptides are stored sorted by
//!    mass, so partners are found with a binary search. Every (seed, partner)
//!    pair, and every pair of reactive residues, is then scored.
//!
//! # Fragment ion model
//!
//! For chain `A` linked at residue `i` to chain `B`, every fragment ion of `A`
//! that *contains* residue `i` carries the whole other chain plus the linker:
//! `shift = m(B) + m(linker)`. Ions that don't contain residue `i` are ordinary
//! linear ions. The same holds for `B`. (Loop-links: ions between the two
//! linked residues can't be observed with a single backbone cleavage, and ions
//! containing both carry only the linker mass.)
//!
//! # Scope / limitations
//!
//! * Non-cleavable linkers only. MS-cleavable linkers (DSSO, DSBU) need a
//!   different search strategy (diagnostic-ion doublets).
//! * Mono-links (dead-end/hydrolyzed linker) are *not* handled here; search for
//!   them as an ordinary variable modification (e.g. `K: 156.0786`).
//! * The score is the joint X!Tandem-style hyperscore over both chains.

use crate::crosslink_kernel::{ChainKernel, KernelParams};
use crate::database::PeptideIx;
use crate::enzyme::Position;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::modification::Modifications;
use crate::peptide::Peptide;
use crate::scoring::{ScoreType, Scorer};
use crate::spectrum::{select_most_intense_peak, ProcessedSpectrum};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Peptides that explain fewer than this many index fragments are never seeds
const MIN_SEED_PEAKS: u16 = 2;
/// Cleavable mode: chains kept per peptide mass, per side, before pairing
const MAX_CHAINS_PER_SIDE: usize = 10;
/// Safety valve for absurdly wide precursor tolerances
const MAX_PARTNERS_PER_SEED: usize = 5_000;
/// Multivalent search: most linkers per species
pub const MAX_VALENCY: u8 = 4;

#[path = "crosslink_multi.rs"]
mod multi;

fn default_residues() -> String {
    "K".into()
}
fn default_true() -> bool {
    true
}
fn default_candidates() -> usize {
    30
}
fn default_linked_peaks() -> u16 {
    2
}

fn default_linear_search() -> bool {
    true
}

fn default_chain_survivors() -> usize {
    20_000
}
fn default_chain_peaks() -> u16 {
    3
}

/// Mass remnants ("stubs") left on each peptide when an MS-cleavable linker
/// breaks in the collision cell. One cleavage leaves the short stub on one
/// peptide and the long stub on the other, so a short and a long stub sum to
/// the linker mass.
///
/// DSSO: short = alkene 54.01056, long = sulfenic acid 103.99320 (doublet
/// spacing 49.98264); `alt` may hold the unsaturated thiol 85.98264 (sulfenic
/// acid - H2O). DSBU: 85.05276 / 111.03203.
///
/// Asymmetric linkers leave different stubs depending on the arm: `short`,
/// `long` and `alt` then describe the first arm and `other` the second. A
/// cleavage leaves one arm's short stub and the other arm's long stub, so
/// `short + other.long` and `other.short + long` both equal the linker mass
/// (cliXlink: 54.01056 + 213.04596 = 163.06333 + 103.99320 = 267.05653).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CleavableStubs {
    pub short: f32,
    pub long: f32,
    /// Additional observed stub masses (e.g. DSSO's 85.98264)
    #[serde(default)]
    pub alt: Vec<f32>,
    /// Asymmetric linkers: the stubs of the second arm; absent when both arms
    /// leave the same stubs
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub other: Option<ArmStubs>,
}

/// The stubs one arm of an MS-cleavable linker can leave on its peptide
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ArmStubs {
    pub short: f32,
    pub long: f32,
    #[serde(default)]
    pub alt: Vec<f32>,
}

impl ArmStubs {
    /// Spacing of the doublet (short / long stub) of a peptide on this arm
    pub fn delta(&self) -> f32 {
        self.long - self.short
    }

    /// Every stub mass of this arm
    pub fn all(&self) -> Vec<f32> {
        [self.short, self.long].into_iter().chain(self.alt.iter().copied()).collect()
    }
}

impl CleavableStubs {
    /// The stubs of each arm: one entry when both arms are alike, two for an
    /// asymmetric linker (first arm, second arm)
    pub fn arms(&self) -> Vec<ArmStubs> {
        let first = ArmStubs { short: self.short, long: self.long, alt: self.alt.clone() };
        match &self.other {
            None => vec![first],
            Some(second) => vec![first, second.clone()],
        }
    }

    /// Ways to place the arms on (peptide a, peptide b), as indices into
    /// [`CleavableStubs::arms`]: an asymmetric linker may sit either way round
    pub fn arm_pairs(&self) -> Vec<(usize, usize)> {
        if self.other.is_some() {
            vec![(0, 1), (1, 0)]
        } else {
            vec![(0, 0)]
        }
    }
}

/// Common reagents. Selecting one fixes the linker mass (and, for cleavable
/// reagents, the stub masses; for heterobifunctional ones, the two sites).
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Reagent {
    Dss,
    Bs3,
    Dsg,
    Bs2g,
    Dsso,
    Dsbu,
    /// Zero-length carbodiimide: amide bond between an amine and a carboxyl
    Edc,
    /// Zero-length triazine coupling reagent, same chemistry as EDC
    Dmtmm,
    /// Asymmetric MS-cleavable sulfoxide with an alkyne handle for click
    /// enrichment (XLMOD:02253; PMID 31593346)
    Clixlink,
}

/// Zero-length amide bond formation releases one water
const ZERO_LENGTH: f32 = -18.010565;

/// Link-site role: can react with the linker's first end ([`CrosslinkSettings::site1`])
pub const SITE1: u8 = 1;
/// Link-site role: can react with the linker's second end ([`CrosslinkSettings::site2`])
pub const SITE2: u8 = 2;

/// Can a site with roles `a` be linked to a site with roles `b`? One must take
/// the first end of the linker and the other the second.
pub fn pairable(a: u8, b: u8) -> bool {
    (a & SITE1 != 0 && b & SITE2 != 0) || (a & SITE2 != 0 && b & SITE1 != 0)
}

impl Reagent {
    pub fn mass(self) -> f32 {
        match self {
            Reagent::Dss | Reagent::Bs3 => 138.06808,
            Reagent::Dsg | Reagent::Bs2g => 96.02113,
            Reagent::Dsso => 158.00376,
            Reagent::Dsbu => 196.08479,
            Reagent::Edc | Reagent::Dmtmm => ZERO_LENGTH,
            // C12H13NO4S
            Reagent::Clixlink => 267.05653,
        }
    }

    /// `(site1, site2)` for reagents whose two ends react with different groups
    pub fn sites(self) -> Option<(&'static str, &'static str)> {
        match self {
            // amine (Lys, protein N-terminus) to carboxyl (Asp, Glu, protein C-terminus)
            Reagent::Edc | Reagent::Dmtmm => Some(("K[", "DE]")),
            _ => None,
        }
    }

    pub fn stubs(self) -> Option<CleavableStubs> {
        match self {
            Reagent::Dsso => Some(CleavableStubs {
                short: 54.01056,
                long: 103.99320,
                alt: vec![85.98264],
                other: None,
            }),
            Reagent::Dsbu => Some(CleavableStubs {
                short: 85.05276,
                long: 111.03203,
                alt: vec![],
                other: None,
            }),
            // Sulfoxide cleavage on either side, as for DSSO: one arm keeps an
            // alkene, the other a sulfenic acid (or the thiol, after losing water).
            // Plain arm: C3H2O / C3H4O2S / C3H2OS. Alkyne-tagged arm:
            // C9H9NO2 / C9H11NO3S / C9H9NO2S.
            Reagent::Clixlink => Some(CleavableStubs {
                short: 54.01056,
                long: 103.99320,
                alt: vec![85.98264],
                other: Some(ArmStubs {
                    short: 163.06333,
                    long: 213.04596,
                    alt: vec![195.03540],
                }),
            }),
            _ => None,
        }
    }
}

/// Which search algorithm / hardware to use for non-cleavable reagents.
///
/// * `cpu`: seed search + partner lookup (default; the validated algorithm)
/// * `reflect`: the data-parallel chain-scoring algorithm, executed on the CPU.
///   A *reference* implementation of the GPU algorithm: slow at proteome scale,
///   but it defines the results the GPU kernel must reproduce
/// * `gpu`: the same algorithm with chain scoring on a GPU (requires building
///   with `--features gpu` and a Vulkan/DX12/Metal device; an error otherwise)
/// * `auto`: `gpu` if available, else `cpu`
///
/// MS-cleavable reagents always use the doublet algorithm on the CPU: it is
/// already cheap and has nothing worth offloading.
#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Backend {
    #[default]
    Cpu,
    Reflect,
    Gpu,
    Auto,
}

/// Crosslinker chemistry and search settings (the `crosslink` block of the
/// Sage JSON configuration)
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CrosslinkSettings {
    /// Reagent preset (`dss`, `bs3`, `dsg`, `bs2g`, `dsso`, `dsbu`). Sets the
    /// linker mass and, for cleavable reagents, the stub masses. Default DSS/BS3.
    #[serde(default)]
    pub reagent: Option<Reagent>,
    /// Mass added to the *pair* of peptides when the linker bridges them.
    /// Only for custom reagents; cannot be combined with `reagent`.
    #[serde(default)]
    pub mass: Option<f32>,
    /// Custom MS-cleavable stub masses. With a cleavable `reagent` (DSSO, DSBU)
    /// they replace its preset stubs; not allowed with a non-cleavable one.
    #[serde(default)]
    pub cleavable: Option<CleavableStubs>,
    /// Residues the linker can react with (both ends alike)
    #[serde(default = "default_residues")]
    pub residues: String,
    /// Heterobifunctional chemistry: what the *first* end reacts with, as
    /// residue letters plus `[` for the protein N-terminus and `]` for the
    /// protein C-terminus (EDC/DMTMM: `"K["`). Needs `site2`; replaces
    /// `residues` and `protein_nterm`. A link always joins a site1 group to a
    /// site2 group.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site1: Option<String>,
    /// Heterobifunctional chemistry: what the *second* end reacts with
    /// (EDC/DMTMM: `"DE]"`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site2: Option<String>,
    /// Allow the linker to react with the protein N-terminal amine
    #[serde(default = "default_true")]
    pub protein_nterm: bool,
    /// Allow linking the C-terminal residue of a peptide that is not at the
    /// protein C-terminus. Normally false: a crosslinked K/R is not cleaved
    /// by trypsin, so it can't be the C-terminal residue of an internal peptide.
    #[serde(default)]
    pub allow_cterm_residue: bool,
    /// Also search for loop-links (both ends of the linker on one peptide)
    #[serde(default = "default_true")]
    pub loop_links: bool,
    /// Multivalent crosslinks: also search species held together by 2 up to
    /// this many linkers (2–4): several links between two peptides (`a2b`),
    /// chains and branches of three or more peptides (`a1b1c`), copies of one
    /// peptide (`a2a`, `a1a1b`). Non-cleavable and MS-cleavable linkers. The search space
    /// grows steeply; meant for small protein sets. `None`: off.
    #[serde(default)]
    pub multivalent: Option<u8>,
    /// Number of top-scoring linear peptides per spectrum used as seeds
    #[serde(default = "default_candidates")]
    pub candidates: usize,
    /// Minimum matched fragment ions required on *each* chain
    #[serde(default = "default_chain_peaks")]
    pub min_chain_peaks: u16,
    /// MS-cleavable linkers only: demand both stub peaks of *both* peptides
    /// (default: both of one peptide, at least one of the partner)
    #[serde(default)]
    pub require_complete_doublets: bool,
    /// Search algorithm / hardware, see [`Backend`]
    #[serde(default)]
    pub backend: Backend,
    /// `reflect`/`gpu` backends: peptides kept after the per-chain filter, per
    /// spectrum and precursor hypothesis (best by chain score). Larger is more
    /// sensitive and slower; a true chain usually ranks far inside the top few hundred.
    #[serde(default = "default_chain_survivors")]
    pub chain_survivors: usize,
    /// Also run an ordinary linear peptide search on the same spectra
    /// (`results.sage.tsv`, `peptides.sage.tsv` and `proteins.sage.tsv` next
    /// to `crosslinks.sage.tsv`). On unless set to `false`.
    #[serde(default = "default_linear_search")]
    pub linear_search: bool,
    /// Non-cleavable linkers: matched fragment ions that carry the linker
    /// (contain the link site; loop-links: both sites) required on *each*
    /// chain. These are the only direct evidence of the link: without them a
    /// correct peptide with any modification of about the linker's mass, or a
    /// strong peptide with a barely matched partner, passes as a crosslink.
    #[serde(default = "default_linked_peaks")]
    pub min_linked_peaks: u16,
}

impl Default for CrosslinkSettings {
    fn default() -> Self {
        Self {
            reagent: None,
            mass: None,
            cleavable: None,
            residues: default_residues(),
            site1: None,
            site2: None,
            protein_nterm: true,
            allow_cterm_residue: false,
            loop_links: true,
            multivalent: None,
            candidates: default_candidates(),
            min_chain_peaks: default_chain_peaks(),
            require_complete_doublets: false,
            backend: Backend::Cpu,
            chain_survivors: default_chain_survivors(),
            linear_search: default_linear_search(),
            min_linked_peaks: default_linked_peaks(),
        }
    }
}

impl CrosslinkSettings {
    /// Mass added to the peptide pair by the linker (DSS/BS3 if nothing is set)
    pub fn linker_mass(&self) -> f32 {
        self.mass
            .or(self.reagent.map(Reagent::mass))
            .unwrap_or(Reagent::Dss.mass())
    }

    /// Stub masses if the linker is MS-cleavable
    pub fn stubs(&self) -> Option<CleavableStubs> {
        self.cleavable
            .clone()
            .or(self.reagent.and_then(Reagent::stubs))
    }

    pub fn is_cleavable(&self) -> bool {
        self.stubs().is_some()
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.reagent.is_some() && self.mass.is_some() {
            return Err(
                "crosslink: set either `reagent`, or a custom `mass` (+ `cleavable`), not both"
                    .into(),
            );
        }
        // `cleavable` may override the stub masses of an MS-cleavable preset
        if self.cleavable.is_some() && self.reagent.is_some_and(|r| r.stubs().is_none()) {
            return Err(
                "crosslink.cleavable can only be combined with an MS-cleavable reagent (dsso, dsbu)"
                    .into(),
            );
        }
        let mass = self.linker_mass();
        // Zero-length reagents (EDC, DMTMM) *remove* a water: -18.010565
        if !mass.is_finite() || mass < -20.0 {
            return Err(format!(
                "crosslink.mass must be positive, or about -18.0106 for a zero-length linker; got {}",
                mass
            ));
        }
        if let Some(stubs) = self.stubs() {
            let arms = stubs.arms();
            for arm in &arms {
                if !(arm.short > 0.0 && arm.long > arm.short) {
                    return Err("crosslink.cleavable requires 0 < short < long (on each arm)".into());
                }
                if arm.alt.iter().any(|&m| !(m > 0.0)) {
                    return Err("crosslink.cleavable.alt masses must be positive".into());
                }
            }
            // a cleavage leaves one arm's short stub and the other arm's long stub
            for (a, b) in stubs.arm_pairs() {
                let sum = arms[a].short + arms[b].long;
                if (sum - mass).abs() > 0.02 {
                    let what = if stubs.other.is_some() {
                        "one arm's short stub + the other arm's long stub"
                    } else {
                        "short + long stubs"
                    };
                    return Err(format!(
                        "crosslink.cleavable: {what} ({sum}) must equal the linker mass ({mass})"
                    ));
                }
            }
        }
        match (&self.site1, &self.site2) {
            (None, None) => {
                if self.residues.is_empty()
                    || !self.residues.bytes().all(|b| b.is_ascii_uppercase())
                {
                    return Err(format!(
                        "crosslink.residues must be uppercase amino acid letters, got `{}`",
                        self.residues
                    ));
                }
            }
            (Some(_), None) | (None, Some(_)) => {
                return Err("crosslink: set both `site1` and `site2`, or neither".into());
            }
            (Some(s1), Some(s2)) => {
                for (key, site) in [("site1", s1), ("site2", s2)] {
                    if site.is_empty()
                        || !site
                            .bytes()
                            .all(|b| b.is_ascii_uppercase() || b == b'[' || b == b']')
                    {
                        return Err(format!(
                            "crosslink.{key} must be uppercase amino acid letters, plus `[` (protein \
                             N-terminus) or `]` (protein C-terminus), got `{site}`"
                        ));
                    }
                }
            }
        }
        if self.candidates == 0 {
            return Err("crosslink.candidates must be at least 1".into());
        }
        if self.chain_survivors == 0 {
            return Err("crosslink.chain_survivors must be at least 1".into());
        }
        if self.min_chain_peaks == 0 {
            return Err("crosslink.min_chain_peaks must be at least 1".into());
        }
        if let Some(n) = self.multivalent {
            if !(2..=MAX_VALENCY).contains(&n) {
                return Err(format!("crosslink.multivalent must be 2 to {MAX_VALENCY} (linkers per species), got {n}"));
            }
        }
        Ok(())
    }

    /// Carbon and sulfur atoms the linker adds to a crosslinked species, for MS1
    /// isotope envelopes. Custom reagents are estimated from their mass.
    pub fn linker_composition(&self) -> (u16, u16) {
        if let Some(mass) = self.mass {
            return ((mass / 14.0).round().max(0.0) as u16, 0);
        }
        match self.reagent.unwrap_or(Reagent::Dss) {
            Reagent::Dss | Reagent::Bs3 => (8, 0),   // C8H10O2
            Reagent::Dsg | Reagent::Bs2g => (5, 0),  // C5H4O2
            Reagent::Dsso => (6, 1),                 // C6H6O3S
            Reagent::Dsbu => (9, 0),                 // C9H14N2O3
            Reagent::Edc | Reagent::Dmtmm => (0, 0), // -H2O
            Reagent::Clixlink => (12, 1),            // C12H13NO4S
        }
    }

    /// The crosslinker as a ProForma 2.0 modification (without its `#XL`
    /// label): its XLMOD term where ProForma readers such as Annotator know the
    /// term's mass (cleavable ones then also bring their stub forms), otherwise
    /// its formula or mass delta
    pub fn proforma_linker(&self) -> String {
        match (self.reagent, self.mass) {
            (Some(Reagent::Dss), _) | (None, None) => "XLMOD:02001".into(),
            (Some(Reagent::Bs3), _) => "XLMOD:02000".into(),
            (Some(Reagent::Dsg), _) => "XLMOD:02006".into(),
            (Some(Reagent::Bs2g), _) => "XLMOD:02005".into(),
            (Some(Reagent::Dsso), _) => "XLMOD:02126".into(),
            (Some(Reagent::Dsbu), _) => "XLMOD:02120".into(),
            (Some(Reagent::Dmtmm), _) => "XLMOD:02208".into(),
            // zero-length amide bond: one water lost
            (Some(Reagent::Edc), _) => "Formula:H-2O-1".into(),
            // XLMOD:02253 carries no mass in Annotator's copy of XLMOD; the
            // exact mass, as f32 would print +267.05652
            (Some(Reagent::Clixlink), _) => "+267.05653".into(),
            (None, Some(mass)) => format!("{mass:+}"),
        }
    }

    /// The two site specifications of a heterobifunctional reagent, if any
    pub fn hetero_sites(&self) -> Option<(String, String)> {
        match (&self.site1, &self.site2) {
            (Some(s1), Some(s2)) => Some((s1.clone(), s2.clone())),
            _ => self
                .reagent
                .and_then(Reagent::sites)
                .map(|(s1, s2)| (s1.into(), s2.into())),
        }
    }

    /// Do the two ends of the linker react with different groups?
    pub fn is_heterobifunctional(&self) -> bool {
        self.hetero_sites().is_some()
    }

    /// 0-based residue indices of `peptide` that can carry a crosslink (either
    /// end of the linker, see [`CrosslinkSettings::link_site_roles`]).
    pub fn link_sites(&self, peptide: &Peptide) -> Vec<usize> {
        self.link_site_roles(peptide)
            .into_iter()
            .map(|(site, _)| site)
            .collect()
    }

    /// Link sites of `peptide` with the linker end(s) each can take: a bitmask
    /// of [`SITE1`] and [`SITE2`]. Both bits are set for homobifunctional
    /// reagents. Use [`pairable`] to check two sites against each other.
    ///
    /// A residue is excluded if it is already modified (this also covers
    /// static modifications that block the amine, e.g. TMT), and a peptide's
    /// C-terminal residue is excluded unless it is the protein C-terminus
    /// (see [`CrosslinkSettings::allow_cterm_residue`]).
    pub fn link_site_roles(&self, peptide: &Peptide) -> Vec<(usize, u8)> {
        let Some((site1, site2)) = self.hetero_sites() else {
            let sites = self.sites_of(peptide, self.residues.as_bytes(), self.protein_nterm, false);
            return sites.into_iter().map(|s| (s, SITE1 | SITE2)).collect();
        };
        let mut roles = std::collections::BTreeMap::<usize, u8>::new();
        for (spec, role) in [(site1, SITE1), (site2, SITE2)] {
            let residues = spec.bytes().filter(u8::is_ascii_uppercase).collect::<Vec<_>>();
            let sites = self.sites_of(peptide, &residues, spec.contains('['), spec.contains(']'));
            for site in sites {
                *roles.entry(site).or_default() |= role;
            }
        }
        roles.into_iter().collect()
    }

    /// Sites reacting with one group: `residues`, optionally the protein
    /// N-terminal amine and the protein C-terminal carboxyl
    fn sites_of(&self, peptide: &Peptide, residues: &[u8], nterm: bool, cterm: bool) -> Vec<usize> {
        let n = peptide.sequence.len();
        let mut sites = Vec::new();
        for (i, aa) in peptide.sequence.iter().enumerate() {
            if !residues.contains(aa) || peptide.modifications[i].is_some() {
                continue;
            }
            if i + 1 == n
                && !self.allow_cterm_residue
                && !matches!(peptide.position, Position::Cterm | Position::Full)
            {
                continue;
            }
            sites.push(i);
        }
        // The protein N-terminal amine reacts like a residue-0 site: every b-ion
        // contains it, so it shares residue 0's fragment-ion model.
        if nterm
            && matches!(peptide.position, Position::Nterm | Position::Full)
            && peptide.nterm.is_none()
            && !sites.contains(&0)
        {
            sites.insert(0, 0);
        }
        // Likewise the protein C-terminal carboxyl and the last residue: every
        // y-ion contains it
        if cterm
            && n > 0
            && matches!(peptide.position, Position::Cterm | Position::Full)
            && peptide.cterm.is_none()
            && !sites.contains(&(n - 1))
        {
            sites.push(n - 1);
        }
        sites
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum LinkKind {
    /// Two peptides joined by the linker
    #[default]
    Inter,
    /// One peptide linked to itself
    Loop,
    /// Multivalent: two or more linkers holding peptides together (several
    /// links between two peptides, or three or more peptides); see [`Species`]
    Multi,
}

/// A multivalent crosslinked species: peptides (copies of one peptide count
/// separately) and the links between them
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Species {
    pub peptides: Vec<PeptideIx>,
    /// `(peptide, site, peptide, site)`: indices into `peptides` and 0-based
    /// link residues
    pub links: Vec<(u8, u16, u8, u16)>,
    /// Matched fragment ions per peptide, and those carrying a linker
    pub matched: Vec<u16>,
    pub linked: Vec<u16>,
    /// Decoy peptides in the species
    pub decoys: u8,
    /// MS-cleavable linkers: peaks of parts of the species freed by breaking a
    /// single linker (`b1c` + stub from `a1b1c`), the evidence for which
    /// peptides are linked to which; 0 otherwise
    pub topology_peaks: u16,
}

impl Species {
    /// Topology in letters: distinct peptides `a`, `b`, … in order, links as
    /// their count: `a2b`, `a1b1c`, `a1a1b`; a branched species lists its
    /// edges, `a1b,b1c,b1d`
    pub fn topology(&self) -> String {
        let mut letters: Vec<PeptideIx> = Vec::new();
        let letter = |ix: PeptideIx, letters: &mut Vec<PeptideIx>| {
            let i = letters.iter().position(|&l| l == ix).unwrap_or_else(|| {
                letters.push(ix);
                letters.len() - 1
            });
            (b'a' + i as u8) as char
        };
        let n = self.peptides.len();
        // edge multiplicities between nodes
        let mut count: Vec<((u8, u8), usize)> = Vec::new();
        for &(a, _, b, _) in &self.links {
            let key = (a.min(b), a.max(b));
            match count.iter_mut().find(|(k, _)| *k == key) {
                Some((_, c)) => *c += 1,
                None => count.push((key, 1)),
            }
        }
        count.sort();
        let degree = |v: u8| count.iter().filter(|((a, b), _)| *a == v || *b == v).count();
        // a chain: every peptide meets at most two others; walk it from an end
        let is_chain = (0..n as u8).all(|v| degree(v) <= 2) && count.len() + 1 == n;
        if is_chain {
            let start = (0..n as u8).find(|&v| degree(v) <= 1).unwrap_or(0);
            let mut out = String::new();
            let (mut prev, mut cur) = (u8::MAX, start);
            for _ in 0..n {
                out.push(letter(self.peptides[cur as usize], &mut letters));
                let next = count
                    .iter()
                    .find(|((a, b), _)| (*a == cur && *b != prev) || (*b == cur && *a != prev));
                let Some(&((a, b), k)) = next else { break };
                out.push_str(&k.to_string());
                prev = cur;
                cur = if a == cur { b } else { a };
            }
            return out;
        }
        let parts: Vec<String> = count
            .iter()
            .map(|&((a, b), k)| {
                let la = letter(self.peptides[a as usize], &mut letters);
                let lb = letter(self.peptides[b as usize], &mut letters);
                format!("{la}{k}{lb}")
            })
            .collect();
        parts.join(",")
    }
}

/// Target/decoy composition of a match. For loop-links (one peptide) the class
/// is `TT` for a target peptide and `DD` for a decoy.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub enum PairClass {
    #[default]
    TT,
    TD,
    DD,
}

impl PairClass {
    fn of(a_decoy: bool, b_decoy: bool) -> Self {
        match (a_decoy, b_decoy) {
            (false, false) => PairClass::TT,
            (true, true) => PairClass::DD,
            _ => PairClass::TD,
        }
    }
}

/// A crosslink-spectrum match (CSM)
#[derive(Clone, Debug, Default, Serialize)]
pub struct CrosslinkMatch {
    pub psm_id: usize,
    pub spec_id: String,
    pub file_id: usize,
    pub rank: u32,
    pub kind: LinkKind,
    pub class: PairClass,
    /// Better-supported peptide of the pair
    pub alpha: PeptideIx,
    /// `None` for loop-links
    pub beta: Option<PeptideIx>,
    /// 0-based link residue on `alpha`
    pub alpha_site: u16,
    /// 0-based link residue on `beta`; for loop-links, the second link residue on `alpha`
    pub beta_site: u16,
    /// Do the two peptides share a protein? (always true for loop-links)
    pub intra_protein: bool,
    pub charge: u8,
    pub rt: f32,
    pub expmass: f32,
    pub calcmass: f32,
    pub precursor_ppm: f32,
    pub isotope_error: f32,
    /// Joint hyperscore over both chains
    pub hyperscore: f64,
    pub alpha_score: f64,
    pub beta_score: f64,
    pub alpha_matched: u16,
    pub beta_matched: u16,
    /// Non-cleavable linkers: matched ions carrying the linker, per chain
    /// (loop-links: `alpha_linked` only). Cleavable loop-links: matched
    /// stub-carrying ions inside the opened ring.
    pub alpha_linked: u16,
    pub beta_linked: u16,
    pub delta_next: f64,
    /// Number of candidate pairs scored for this spectrum
    pub pairs_scored: u32,
    /// MS-cleavable linkers: number of stub-peptide (doublet) peaks supporting the pair (0-4)
    pub stub_peaks: u8,
    /// Precursor inverse ion mobility (0 if the instrument doesn't report one)
    pub ims: f32,
    pub spectrum_q: f32,
    pub crosslink_q: f32,
    /// Multivalent species (`kind == Multi`): every peptide and link. `alpha`
    /// / `beta` and their sites are then the species' first link.
    pub species: Option<Box<Species>>,
}

static CSM_COUNTER: AtomicUsize = AtomicUsize::new(1);

/// Matched-fragment tally, same quantities that feed Sage's hyperscore
#[derive(Copy, Clone, Default, Debug)]
struct Tally {
    b: u16,
    y: u16,
    sum_b: f32,
    sum_y: f32,
}

impl Tally {
    fn matched(&self) -> u16 {
        self.b + self.y
    }

    fn combine(&self, other: &Tally) -> Tally {
        Tally {
            b: self.b + other.b,
            y: self.y + other.y,
            sum_b: self.sum_b + other.sum_b,
            sum_y: self.sum_y + other.sum_y,
        }
    }

    fn score(&self, score_type: ScoreType) -> f64 {
        score_type.score(self.b, self.y, self.sum_b, self.sum_y)
    }
}

/// How a chain is attached to the rest of the crosslinked species
#[derive(Copy, Clone, Debug)]
pub enum Link {
    /// Linked at residue `site` to another chain; ions containing `site` carry `shift`
    /// (= mass of the other chain + linker)
    Single { site: usize, shift: f32 },
    /// Loop-link between residues `first < second` of the same chain
    Loop {
        first: usize,
        second: usize,
        shift: f32,
    },
}

/// Theoretical (neutral) fragment ions of one chain of a crosslinked species.
///
/// Ion `idx` of an N-terminal series (a/b/c) contains residues `0..=idx`; ion
/// `idx` of a C-terminal series (x/y/z) contains residues `idx+1..`.
pub fn chain_ions(
    peptide: &Peptide,
    mods: &Modifications,
    kinds: &[Kind],
    link: Link,
) -> Vec<(Kind, f32)> {
    let mut ions = Vec::with_capacity(kinds.len() * peptide.sequence.len());
    for &kind in kinds {
        let n_terminal = matches!(kind, Kind::A | Kind::B | Kind::C);
        let contains = |site: usize, idx: usize| {
            if n_terminal {
                idx >= site
            } else {
                idx < site
            }
        };
        for (idx, ion) in IonSeries::new(peptide, kind, mods).enumerate() {
            let mass = match link {
                Link::Single { site, shift } => {
                    ion.monoisotopic_mass + if contains(site, idx) { shift } else { 0.0 }
                }
                Link::Loop {
                    first,
                    second,
                    shift,
                } => match (contains(first, idx), contains(second, idx)) {
                    (false, false) => ion.monoisotopic_mass,
                    (true, true) => ion.monoisotopic_mass + shift,
                    // Cleaving inside the ring doesn't separate fragments
                    _ => continue,
                },
            };
            ions.push((kind, mass));
        }
    }
    ions
}

/// The ions of [`chain_ions`] that carry the link's `shift`: the ones that
/// contain the link site (loop-links: both sites)
pub fn linked_ions(peptide: &Peptide, mods: &Modifications, kinds: &[Kind], link: Link) -> Vec<(Kind, f32)> {
    let unshifted = match link {
        Link::Single { site, .. } => Link::Single { site, shift: 0.0 },
        Link::Loop { first, second, .. } => Link::Loop { first, second, shift: 0.0 },
    };
    chain_ions(peptide, mods, kinds, link)
        .into_iter()
        .zip(chain_ions(peptide, mods, kinds, unshifted))
        .filter(|(shifted, plain)| shifted.1 != plain.1)
        .map(|(ion, _)| ion)
        .collect()
}

/// Theoretical (neutral) fragment ions of a loop-link of an MS-cleavable
/// linker whose ring opened in the collision cell: `(all ions, ring ions)`.
///
/// Residue `first < second` carry the arms `at_first` / `at_second`. Once the
/// linker is cleaved, a backbone cleavage *inside* the ring separates
/// fragments too: those contain one linked residue and carry one of its arm's
/// stubs (any of short, long, alt). These ring ions are the evidence of the
/// loop-link and are returned separately as well. Ions containing both
/// residues carry both stubs, one arm's short stub and the other arm's long
/// (= the linker mass) or alt stub, or the intact linker; ions containing
/// neither are ordinary linear ions.
#[allow(clippy::too_many_arguments)]
pub fn cleaved_loop_ions(
    peptide: &Peptide,
    mods: &Modifications,
    kinds: &[Kind],
    first: usize,
    second: usize,
    at_first: &ArmStubs,
    at_second: &ArmStubs,
    linker: f32,
) -> (Vec<(Kind, f32)>, Vec<(Kind, f32)>) {
    // stubs summed over both sites: exactly one of them the short stub
    let mut both = vec![linker];
    for (x, y) in [(at_first, at_second), (at_second, at_first)] {
        for long in std::iter::once(y.long).chain(y.alt.iter().copied()) {
            let sum = x.short + long;
            if both.iter().all(|b| (b - sum).abs() > 1e-3) {
                both.push(sum);
            }
        }
    }
    let (mut all, mut ring) = (Vec::new(), Vec::new());
    for &kind in kinds {
        let n_terminal = matches!(kind, Kind::A | Kind::B | Kind::C);
        let contains = |site: usize, idx: usize| if n_terminal { idx >= site } else { idx < site };
        for (idx, ion) in IonSeries::new(peptide, kind, mods).enumerate() {
            let mass = ion.monoisotopic_mass;
            match (contains(first, idx), contains(second, idx)) {
                (false, false) => all.push((kind, mass)),
                (true, true) => all.extend(both.iter().map(|s| (kind, mass + s))),
                (true, false) => ring.extend(at_first.all().into_iter().map(|s| (kind, mass + s))),
                (false, true) => ring.extend(at_second.all().into_iter().map(|s| (kind, mass + s))),
            }
        }
    }
    all.extend(ring.iter().copied());
    (all, ring)
}

fn frag_charge_bound(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    // Same rule as the linear search in `scoring.rs`
    precursor_charge
        .min(max_fragment_charge.map(|c| c + 1).unwrap_or(precursor_charge))
        .max(2)
}

fn match_ions(
    query: &ProcessedSpectrum,
    ions: &[(Kind, f32)],
    max_fragment_charge: u8,
    tolerance: Tolerance,
) -> Tally {
    let mut tally = Tally::default();
    for &(kind, mass) in ions {
        for charge in 1..max_fragment_charge {
            let mz = mass / charge as f32;
            if let Some(idx) = select_most_intense_peak(
                &query.masses,
                &query.intensities,
                mz,
                tolerance,
                None,
            ) {
                let intensity = query.intensities[idx];
                match kind {
                    Kind::A | Kind::B | Kind::C => {
                        tally.b += 1;
                        tally.sum_b += intensity;
                    }
                    _ => {
                        tally.y += 1;
                        tally.sum_y += intensity;
                    }
                }
            }
        }
    }
    tally
}

struct PairEval {
    score: f64,
    site_a: usize,
    site_b: usize,
    tally_a: Tally,
    tally_b: Tally,
}

/// Searches spectra for crosslinked peptides. Wraps a linear [`Scorer`] and
/// reuses its tolerances, charge/isotope settings, ion types and score type.
pub struct CrosslinkScorer<'a, 'db> {
    pub scorer: &'a Scorer<'db>,
    pub settings: &'a CrosslinkSettings,
    /// Chain-scoring kernel for the `reflect`/`gpu` backends (non-cleavable only)
    pub kernel: Option<&'a dyn ChainKernel>,
}

impl<'a, 'db> CrosslinkScorer<'a, 'db> {
    pub fn new(scorer: &'a Scorer<'db>, settings: &'a CrosslinkSettings) -> Self {
        Self {
            scorer,
            settings,
            kernel: None,
        }
    }

    /// Use the data-parallel reflection algorithm with this chain kernel
    pub fn with_kernel(mut self, kernel: Option<&'a dyn ChainKernel>) -> Self {
        self.kernel = kernel;
        self
    }

    /// Stage 1: linear peptides that explain the most peaks, over every peptide
    /// light enough to be one half of the pair
    fn seeds(
        &self,
        query: &ProcessedSpectrum,
        max_mass: f32,
        max_fragment_charge: u8,
    ) -> Vec<PeptideIx> {
        let db = self.scorer.db;
        let candidates = db.query(
            max_mass,
            Tolerance::Da(-max_mass, 0.0),
            self.scorer.fragment_tol,
        );
        let potential = candidates.pre_idx_hi.saturating_sub(candidates.pre_idx_lo) + 1;
        let mut counts = vec![0u16; potential];

        for &peak in query.masses.iter() {
            for charge in 1..max_fragment_charge {
                for frag in candidates.page_search(peak * charge as f32) {
                    let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                    counts[idx] = counts[idx].saturating_add(1);
                }
            }
        }

        let mut hits = counts
            .iter()
            .enumerate()
            .filter(|(_, &c)| c >= MIN_SEED_PEAKS)
            .map(|(i, &c)| (c, (i + candidates.pre_idx_lo) as u32))
            .collect::<Vec<_>>();

        let k = self.settings.candidates.min(hits.len());
        if k == 0 {
            return Vec::new();
        }
        if k < hits.len() {
            hits.select_nth_unstable_by(k - 1, |a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            hits.truncate(k);
        }
        // best first (the multivalent search takes the first few)
        hits.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        hits.into_iter().map(|(_, ix)| PeptideIx(ix)).collect()
    }

    /// Best site combination for a peptide pair.
    ///
    /// The chain tally at a given site doesn't depend on where the *other*
    /// chain is linked, so tallies are computed once per site and combined.
    fn score_pair(
        &self,
        query: &ProcessedSpectrum,
        a: &Peptide,
        b: &Peptide,
        sites_a: &[(usize, u8)],
        sites_b: &[(usize, u8)],
        max_fragment_charge: u8,
    ) -> Option<PairEval> {
        let db = self.scorer.db;
        let kinds = &db.ion_kinds;
        let linker = self.settings.linker_mass();
        let tol = self.scorer.fragment_tol;
        let min_peaks = self.settings.min_chain_peaks;

        let tallies = |p: &Peptide, other: &Peptide, sites: &[(usize, u8)]| {
            sites
                .iter()
                .map(|&(site, roles)| {
                    let ions = chain_ions(
                        p,
                        &db.modifications,
                        kinds,
                        Link::Single {
                            site,
                            shift: other.monoisotopic + linker,
                        },
                    );
                    (site, roles, match_ions(query, &ions, max_fragment_charge, tol))
                })
                .filter(|(_, _, t)| t.matched() >= min_peaks)
                .collect::<Vec<_>>()
        };

        let ta = tallies(a, b, sites_a);
        let tb = tallies(b, a, sites_b);

        let mut best: Option<PairEval> = None;
        for &(site_a, roles_a, tally_a) in &ta {
            for &(site_b, roles_b, tally_b) in &tb {
                if !pairable(roles_a, roles_b) {
                    continue;
                }
                let score = tally_a.combine(&tally_b).score(self.scorer.score_type);
                if best.as_ref().map_or(true, |x| score > x.score) {
                    best = Some(PairEval {
                        score,
                        site_a,
                        site_b,
                        tally_a,
                        tally_b,
                    });
                }
            }
        }
        best
    }

    /// Best pair of link sites for a loop-link on a single peptide
    fn score_loop(
        &self,
        query: &ProcessedSpectrum,
        a: &Peptide,
        sites: &[(usize, u8)],
        max_fragment_charge: u8,
    ) -> Option<(f64, usize, usize, Tally)> {
        let db = self.scorer.db;
        let mut best: Option<(f64, usize, usize, Tally)> = None;
        for (x, &(first, roles_first)) in sites.iter().enumerate() {
            for &(second, roles_second) in &sites[x + 1..] {
                if !pairable(roles_first, roles_second) {
                    continue;
                }
                let ions = chain_ions(
                    a,
                    &db.modifications,
                    &db.ion_kinds,
                    Link::Loop {
                        first,
                        second,
                        shift: self.settings.linker_mass(),
                    },
                );
                let tally = match_ions(query, &ions, max_fragment_charge, self.scorer.fragment_tol);
                if tally.matched() < self.settings.min_chain_peaks {
                    continue;
                }
                let score = tally.score(self.scorer.score_type);
                if best.as_ref().map_or(true, |x| score > x.0) {
                    best = Some((score, first, second, tally));
                }
            }
        }
        best
    }

    /// Search one MS2 spectrum. Returns up to `scorer.report_psms` matches,
    /// best first.
    pub fn score(&self, query: &ProcessedSpectrum) -> Vec<CrosslinkMatch> {
        let mut found = self.score_spectrum(query);
        let ims = query
            .precursors
            .first()
            .and_then(|p| p.inverse_ion_mobility)
            .unwrap_or_default();
        for m in &mut found {
            m.ims = ims;
        }
        found
    }

    fn score_spectrum(&self, query: &ProcessedSpectrum) -> Vec<CrosslinkMatch> {
        let sc = self.scorer;
        let cfg = self.settings;
        let db = sc.db;

        if cfg.is_cleavable() {
            return self.score_cleavable(query);
        }
        if let Some(kernel) = self.kernel {
            return self.score_reflect(query, kernel);
        }
        let Some(precursor) = query.precursors.first() else {
            return Vec::new();
        };
        if db.peptides.is_empty() {
            return Vec::new();
        }
        let min_peptide_mass = db.peptides[0].monoisotopic;
        let linker_mass = cfg.linker_mass();
        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        let charges: Vec<u8> = match precursor.charge {
            Some(z) if !sc.override_precursor_charge => vec![z],
            _ => (sc.min_precursor_charge..=sc.max_precursor_charge).collect(),
        };

        let mut found: Vec<CrosslinkMatch> = Vec::new();
        let mut pairs_scored = 0u32;

        for charge in charges {
            let mass = mz * charge as f32;
            let max_fc = frag_charge_bound(sc.max_fragment_charge, charge);

            let largest = sc
                .precursor_tol
                .bounds(mass - sc.min_isotope_err as f32 * NEUTRON)
                .1;
            let seed_limit = if cfg.loop_links {
                largest - linker_mass
            } else {
                largest - linker_mass - min_peptide_mass
            };
            if seed_limit < min_peptide_mass {
                continue;
            }

            let mut seen: HashSet<(u32, u32)> = HashSet::new();
            for seed in self.seeds(query, seed_limit, max_fc) {
                let a = &db[seed];
                let sites_a = cfg.link_site_roles(a);
                if sites_a.is_empty() {
                    continue;
                }

                for iso in sc.min_isotope_err..=sc.max_isotope_err {
                    let m = mass - iso as f32 * NEUTRON;
                    let (lo, hi) = sc.precursor_tol.bounds(m);

                    // Loop-link: one peptide plus the linker explains the precursor
                    let loop_mass = a.monoisotopic + linker_mass;
                    if cfg.loop_links
                        && sites_a.len() >= 2
                        && loop_mass >= lo
                        && loop_mass <= hi
                        && seen.insert((seed.0, u32::MAX))
                    {
                        pairs_scored += 1;
                        if let Some((score, first, second, tally)) =
                            self.score_loop(query, a, &sites_a, max_fc)
                        {
                            found.push(CrosslinkMatch {
                                psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
                                spec_id: query.id.clone(),
                                file_id: query.file_id,
                                kind: LinkKind::Loop,
                                class: PairClass::of(a.decoy, a.decoy),
                                alpha: seed,
                                beta: None,
                                alpha_site: first as u16,
                                beta_site: second as u16,
                                intra_protein: true,
                                charge,
                                rt: query.scan_start_time,
                                expmass: mass,
                                calcmass: loop_mass,
                                precursor_ppm: ppm_delta(mass, iso, loop_mass),
                                isotope_error: iso as f32 * NEUTRON,
                                hyperscore: score,
                                alpha_score: score,
                                alpha_matched: tally.matched(),
                                ..Default::default()
                            });
                        }
                    }

                    // Inter-peptide link: the partner's mass is fixed by the precursor
                    let target_lo = lo - linker_mass - a.monoisotopic;
                    let target_hi = hi - linker_mass - a.monoisotopic;
                    if target_hi < min_peptide_mass {
                        continue;
                    }
                    let start = db.peptides.partition_point(|p| p.monoisotopic < target_lo);
                    for (offset, b) in db.peptides[start..]
                        .iter()
                        .take_while(|p| p.monoisotopic <= target_hi)
                        .take(MAX_PARTNERS_PER_SEED)
                        .enumerate()
                    {
                        let partner = PeptideIx((start + offset) as u32);
                        let key = (seed.0.min(partner.0), seed.0.max(partner.0));
                        if !seen.insert(key) {
                            continue;
                        }
                        let sites_b = cfg.link_site_roles(b);
                        if sites_b.is_empty() {
                            continue;
                        }
                        pairs_scored += 1;
                        let Some(eval) =
                            self.score_pair(query, a, b, &sites_a, &sites_b, max_fc)
                        else {
                            continue;
                        };

                        let score_a = eval.tally_a.score(sc.score_type);
                        let score_b = eval.tally_b.score(sc.score_type);
                        // alpha is the better-supported chain
                        let swap = score_b > score_a;
                        let (alpha, beta, a_site, b_site, a_score, b_score, a_n, b_n) = if swap {
                            (
                                partner,
                                seed,
                                eval.site_b,
                                eval.site_a,
                                score_b,
                                score_a,
                                eval.tally_b.matched(),
                                eval.tally_a.matched(),
                            )
                        } else {
                            (
                                seed,
                                partner,
                                eval.site_a,
                                eval.site_b,
                                score_a,
                                score_b,
                                eval.tally_a.matched(),
                                eval.tally_b.matched(),
                            )
                        };
                        let calc = a.monoisotopic + b.monoisotopic + linker_mass;
                        found.push(CrosslinkMatch {
                            psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
                            spec_id: query.id.clone(),
                            file_id: query.file_id,
                            kind: LinkKind::Inter,
                            class: PairClass::of(a.decoy, b.decoy),
                            alpha,
                            beta: Some(beta),
                            alpha_site: a_site as u16,
                            beta_site: b_site as u16,
                            intra_protein: a.proteins.iter().any(|p| b.proteins.contains(p)),
                            charge,
                            rt: query.scan_start_time,
                            expmass: mass,
                            calcmass: calc,
                            precursor_ppm: ppm_delta(mass, iso, calc),
                            isotope_error: iso as f32 * NEUTRON,
                            hyperscore: eval.score,
                            alpha_score: a_score,
                            beta_score: b_score,
                            alpha_matched: a_n,
                            beta_matched: b_n,
                            ..Default::default()
                        });
                    }
                }
            }
        }

        self.finalize(query, found, pairs_scored)
    }

    /// Assemble the match for an inter-peptide pair; `a`/`b` are database indices
    #[allow(clippy::too_many_arguments)]
    fn inter_match(
        &self,
        query: &ProcessedSpectrum,
        a: PeptideIx,
        b: PeptideIx,
        eval: &PairEval,
        charge: u8,
        mass: f32,
        iso: i8,
    ) -> CrosslinkMatch {
        let sc = self.scorer;
        let db = sc.db;
        let (pa, pb) = (&db[a], &db[b]);
        let score_a = eval.tally_a.score(sc.score_type);
        let score_b = eval.tally_b.score(sc.score_type);
        // alpha is the better-supported chain
        let swap = score_b > score_a;
        let (alpha, beta, a_site, b_site, a_score, b_score, a_n, b_n) = if swap {
            (b, a, eval.site_b, eval.site_a, score_b, score_a, eval.tally_b.matched(), eval.tally_a.matched())
        } else {
            (a, b, eval.site_a, eval.site_b, score_a, score_b, eval.tally_a.matched(), eval.tally_b.matched())
        };
        let calc = pa.monoisotopic + pb.monoisotopic + self.settings.linker_mass();
        CrosslinkMatch {
            psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
            spec_id: query.id.clone(),
            file_id: query.file_id,
            kind: LinkKind::Inter,
            class: PairClass::of(pa.decoy, pb.decoy),
            alpha,
            beta: Some(beta),
            alpha_site: a_site as u16,
            beta_site: b_site as u16,
            intra_protein: pa.proteins.iter().any(|p| pb.proteins.contains(p)),
            charge,
            rt: query.scan_start_time,
            expmass: mass,
            calcmass: calc,
            precursor_ppm: ppm_delta(mass, iso, calc),
            isotope_error: iso as f32 * NEUTRON,
            hyperscore: eval.score,
            alpha_score: a_score,
            beta_score: b_score,
            alpha_matched: a_n,
            beta_matched: b_n,
            ..Default::default()
        }
    }

    /// Non-cleavable search by chain reflection (see `crosslink_kernel`).
    ///
    /// 1. A [`ChainKernel`] scores *every* peptide on its own, with its
    ///    link-residue ions shifted by `M - m(peptide)`, and returns survivors.
    /// 2. Survivors are joined on mass: `m(B) = M - linker - m(A)`.
    /// 3. The best `candidates` pairs by surrogate score are re-scored exactly
    ///    with the real partner mass (the same scoring as the `cpu` backend).
    /// 4. Loop-links: peptides whose mass plus linker equals the precursor are
    ///    looked up directly, since no partner is involved.
    fn score_reflect(&self, query: &ProcessedSpectrum, kernel: &dyn ChainKernel) -> Vec<CrosslinkMatch> {
        let sc = self.scorer;
        let cfg = self.settings;
        let db = sc.db;
        let Some(precursor) = query.precursors.first() else {
            return Vec::new();
        };
        if db.peptides.is_empty() {
            return Vec::new();
        }
        let linker = cfg.linker_mass();
        let mz = precursor.mz - PROTON;
        let charges: Vec<u8> = match precursor.charge {
            Some(z) if !sc.override_precursor_charge => vec![z],
            _ => (sc.min_precursor_charge..=sc.max_precursor_charge).collect(),
        };

        let mut found: Vec<CrosslinkMatch> = Vec::new();
        let mut pairs_scored = 0u32;
        let mut seen: HashSet<(u32, u32)> = HashSet::new();

        for charge in charges {
            let mass = mz * charge as f32;
            let max_fc = frag_charge_bound(sc.max_fragment_charge, charge);
            for iso in sc.min_isotope_err..=sc.max_isotope_err {
                let m = mass - iso as f32 * NEUTRON;
                let (lo, hi) = sc.precursor_tol.bounds(m);

                // ---- loop-links: one peptide + linker explains the precursor
                if cfg.loop_links {
                    let start = db.peptides.partition_point(|p| p.monoisotopic < lo - linker);
                    for (offset, a) in db.peptides[start..]
                        .iter()
                        .take_while(|p| p.monoisotopic <= hi - linker)
                        .take(MAX_PARTNERS_PER_SEED)
                        .enumerate()
                    {
                        let ix = (start + offset) as u32;
                        let sites = cfg.link_site_roles(a);
                        if sites.len() < 2 || !seen.insert((ix, u32::MAX)) {
                            continue;
                        }
                        pairs_scored += 1;
                        if let Some((score, first, second, tally)) = self.score_loop(query, a, &sites, max_fc) {
                            let calc = a.monoisotopic + linker;
                            found.push(CrosslinkMatch {
                                psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
                                spec_id: query.id.clone(),
                                file_id: query.file_id,
                                kind: LinkKind::Loop,
                                class: PairClass::of(a.decoy, a.decoy),
                                alpha: PeptideIx(ix),
                                beta: None,
                                alpha_site: first as u16,
                                beta_site: second as u16,
                                intra_protein: true,
                                charge,
                                rt: query.scan_start_time,
                                expmass: mass,
                                calcmass: calc,
                                precursor_ppm: ppm_delta(mass, iso, calc),
                                isotope_error: iso as f32 * NEUTRON,
                                hyperscore: score,
                                alpha_score: score,
                                alpha_matched: tally.matched(),
                                ..Default::default()
                            });
                        }
                    }
                }

                // ---- inter-peptide links: chain filter, then join on mass
                let params = KernelParams {
                    precursor_mass: m,
                    max_fragment_charge: max_fc,
                    fragment_tol: sc.fragment_tol.into(),
                    slack: (m - lo).abs().max((hi - m).abs()),
                    min_matched: cfg.min_chain_peaks,
                    max_survivors: cfg.chain_survivors,
                };
                let survivors = kernel.survivors(query, &params);
                let masses = survivors
                    .iter()
                    .map(|s| db.peptides[s.peptide as usize].monoisotopic)
                    .collect::<Vec<_>>();

                let mut candidates: Vec<(f32, usize, usize)> = Vec::new();
                for (i, sa) in survivors.iter().enumerate() {
                    let target = masses[i];
                    let (t_lo, t_hi) = (lo - linker - target, hi - linker - target);
                    let first = masses.partition_point(|&x| x < t_lo).max(i);
                    for j in first..survivors.len() {
                        if masses[j] > t_hi {
                            break;
                        }
                        candidates.push((sa.score + survivors[j].score, i, j));
                    }
                }
                candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
                candidates.truncate(cfg.candidates);

                for (_, i, j) in candidates {
                    let (ix_a, ix_b) = (survivors[i].peptide, survivors[j].peptide);
                    if !seen.insert((ix_a, ix_b)) {
                        continue;
                    }
                    let (a, b) = (&db.peptides[ix_a as usize], &db.peptides[ix_b as usize]);
                    let (sites_a, sites_b) = (cfg.link_site_roles(a), cfg.link_site_roles(b));
                    if sites_a.is_empty() || sites_b.is_empty() {
                        continue;
                    }
                    pairs_scored += 1;
                    if let Some(eval) = self.score_pair(query, a, b, &sites_a, &sites_b, max_fc) {
                        found.push(self.inter_match(query, PeptideIx(ix_a), PeptideIx(ix_b), &eval, charge, mass, iso));
                    }
                }
            }
        }
        self.finalize(query, found, pairs_scored)
    }

    /// Rank candidate matches for one spectrum and keep the best `report_psms`
    fn finalize(&self, query: &ProcessedSpectrum, mut found: Vec<CrosslinkMatch>, mut pairs_scored: u32) -> Vec<CrosslinkMatch> {
        found.retain(|m| m.hyperscore.is_finite());
        if !self.settings.is_cleavable() {
            // cleavable links are evidenced by their stub doublets instead
            for m in &mut found {
                (m.alpha_linked, m.beta_linked) = self.linked_peaks(query, m);
            }
            let min = self.settings.min_linked_peaks;
            found.retain(|m| m.alpha_linked >= min && (m.beta.is_none() || m.beta_linked >= min));
        }
        // multivalent species compete with the pairs and loop-links for the spectrum
        // (their linked ions are checked per peptide while scoring)
        if self.settings.multivalent.is_some() {
            let (multi, scored) = self.score_multivalent(query);
            found.extend(multi.into_iter().filter(|m| m.hyperscore.is_finite()));
            pairs_scored += scored;
        }
        found.sort_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));
        found.truncate(self.scorer.report_psms.max(1));
        for idx in 0..found.len() {
            let next = found.get(idx + 1).map(|m| m.hyperscore).unwrap_or_default();
            found[idx].rank = idx as u32 + 1;
            found[idx].delta_next = found[idx].hyperscore - next;
            found[idx].pairs_scored = pairs_scored;
        }
        found
    }

    /// Matched ions carrying the linker on each chain of a non-cleavable match
    fn linked_peaks(&self, query: &ProcessedSpectrum, m: &CrosslinkMatch) -> (u16, u16) {
        let db = self.scorer.db;
        let max_fc = frag_charge_bound(self.scorer.max_fragment_charge, m.charge);
        let count = |peptide: &Peptide, link: Link| {
            let ions = linked_ions(peptide, &db.modifications, &db.ion_kinds, link);
            match_ions(query, &ions, max_fc, self.scorer.fragment_tol).matched()
        };
        let linker = self.settings.linker_mass();
        let a = &db[m.alpha];
        match m.beta {
            None => {
                let link = Link::Loop { first: m.alpha_site as usize, second: m.beta_site as usize, shift: linker };
                (count(a, link), 0)
            }
            Some(beta) => {
                let b = &db[beta];
                let on_a = Link::Single { site: m.alpha_site as usize, shift: b.monoisotopic + linker };
                let on_b = Link::Single { site: m.beta_site as usize, shift: a.monoisotopic + linker };
                (count(a, on_a), count(b, on_b))
            }
        }
    }

    /// Peptides whose mass lies within the fragment tolerance of `mass`
    fn peptides_with_mass(&self, mass: f32) -> std::ops::Range<usize> {
        let peptides = &self.scorer.db.peptides;
        let (lo, hi) = self.scorer.fragment_tol.bounds(mass);
        let start = peptides.partition_point(|p| p.monoisotopic < lo);
        let len = peptides[start..]
            .iter()
            .take_while(|p| p.monoisotopic <= hi)
            .take(MAX_PARTNERS_PER_SEED)
            .count();
        start..start + len
    }

    /// Best link site for one peptide carrying a single cleavable-linker stub,
    /// separately for each end of the linker: `[best SITE1 site, best SITE2 site]`
    /// (identical for homobifunctional reagents).
    ///
    /// Fragment ions that contain the linked residue carry *a* stub of the
    /// linker arm on this peptide, but which one (short, long, or alt) is not
    /// known, so each is a candidate; ions outside the linked residue are
    /// ordinary linear ions.
    fn score_stubbed_chain(
        &self,
        query: &ProcessedSpectrum,
        peptide: &Peptide,
        arm: &ArmStubs,
        max_fragment_charge: u8,
    ) -> [Option<(usize, Tally)>; 2] {
        let db = self.scorer.db;
        let sites = self.settings.link_site_roles(peptide);
        if sites.is_empty() {
            return [None, None];
        }
        let all_stubs = arm.all();
        let mut best: [Option<(usize, Tally, f64)>; 2] = [None, None];
        for &(site, roles) in &sites {
            let with_stub = |shift: f32| {
                chain_ions(
                    peptide,
                    &db.modifications,
                    &db.ion_kinds,
                    Link::Single { site, shift },
                )
            };
            let unshifted = with_stub(0.0);
            let mut ions: Vec<(Kind, f32)> = Vec::with_capacity(unshifted.len() * 2);
            // Ions without the linked residue are the same whatever the stub
            let first = with_stub(all_stubs[0]);
            for (plain, shifted) in unshifted.iter().zip(&first) {
                if plain.1 == shifted.1 {
                    ions.push(*plain);
                }
            }
            for &stub in &all_stubs {
                for (plain, shifted) in unshifted.iter().zip(with_stub(stub)) {
                    if plain.1 != shifted.1 {
                        ions.push(shifted);
                    }
                }
            }
            let tally = match_ions(query, &ions, max_fragment_charge, self.scorer.fragment_tol);
            if tally.matched() < self.settings.min_chain_peaks {
                continue;
            }
            let score = tally.score(self.scorer.score_type);
            for (end, role) in [SITE1, SITE2].into_iter().enumerate() {
                if roles & role != 0 && best[end].as_ref().map_or(true, |b| score > b.2) {
                    best[end] = Some((site, tally, score));
                }
            }
        }
        best.map(|b| b.map(|(site, tally, _)| (site, tally)))
    }

    /// Best pair of link sites for a loop-link of an MS-cleavable linker on one
    /// peptide (see [`cleaved_loop_ions`]): `(score, first, second, tally,
    /// matched ring ions)`. A site pair needs `min_linked_peaks` matched ring
    /// ions, the ions only an opened loop-link gives.
    fn score_cleavable_loop(
        &self,
        query: &ProcessedSpectrum,
        peptide: &Peptide,
        stubs: &CleavableStubs,
        max_fragment_charge: u8,
    ) -> Option<(f64, usize, usize, Tally, u16)> {
        let db = self.scorer.db;
        let tol = self.scorer.fragment_tol;
        let sites = self.settings.link_site_roles(peptide);
        let arms = stubs.arms();
        let mut best: Option<(f64, usize, usize, Tally, u16)> = None;
        for (x, &(first, roles_first)) in sites.iter().enumerate() {
            for &(second, roles_second) in &sites[x + 1..] {
                if !pairable(roles_first, roles_second) {
                    continue;
                }
                for (arm_first, arm_second) in stubs.arm_pairs() {
                    let (all, ring) = cleaved_loop_ions(
                        peptide,
                        &db.modifications,
                        &db.ion_kinds,
                        first,
                        second,
                        &arms[arm_first],
                        &arms[arm_second],
                        self.settings.linker_mass(),
                    );
                    let tally = match_ions(query, &all, max_fragment_charge, tol);
                    if tally.matched() < self.settings.min_chain_peaks {
                        continue;
                    }
                    let ring_matched = match_ions(query, &ring, max_fragment_charge, tol).matched();
                    if ring_matched < self.settings.min_linked_peaks {
                        continue;
                    }
                    let score = tally.score(self.scorer.score_type);
                    if best.as_ref().is_none_or(|b| score > b.0) {
                        best = Some((score, first, second, tally, ring_matched));
                    }
                }
            }
        }
        best
    }

    /// Search for an MS-cleavable crosslink (DSSO, DSBU, ...).
    ///
    /// 1. A stub-bearing peptide appears in the spectrum as a *doublet* (short
    ///    and long stub of its linker arm), spaced `long - short` apart.
    /// 2. If `a = m(A) + short` is the first peak of A's doublet, the partner
    ///    peptide B (with mass `M - m(A) - linker`) shows up at `m(B) + short`
    ///    and `m(B) + long` of the other arm. For symmetric linkers both arms
    ///    are alike; an asymmetric linker (cliXlink) is tried both ways round.
    /// 3. Both peptide masses are therefore known directly; each is looked up in
    ///    the (mass-sorted) database with a narrow window and scored on its own
    ///    as a linear peptide bearing a stub. No pairwise enumeration is needed.
    fn score_cleavable(&self, query: &ProcessedSpectrum) -> Vec<CrosslinkMatch> {
        let sc = self.scorer;
        let cfg = self.settings;
        let db = sc.db;
        let Some(stubs) = cfg.stubs() else {
            return Vec::new();
        };
        let Some(precursor) = query.precursors.first() else {
            return Vec::new();
        };
        if db.peptides.is_empty() {
            return Vec::new();
        }
        let min_peptide_mass = db.peptides[0].monoisotopic;
        let arms = stubs.arms();
        let linker = cfg.linker_mass();
        let mz = precursor.mz - PROTON;

        let charges: Vec<u8> = match precursor.charge {
            Some(z) if !sc.override_precursor_charge => vec![z],
            _ => (sc.min_precursor_charge..=sc.max_precursor_charge).collect(),
        };

        let mut found: Vec<CrosslinkMatch> = Vec::new();
        let mut pairs_scored = 0u32;
        // (peptide, arm) -> best site per linker end
        let mut chain_cache: HashMap<(u32, usize), [Option<(usize, Tally)>; 2]> = HashMap::new();
        // Homobifunctional: both ends are alike, so one entry per peptide suffices
        let ends = if cfg.is_heterobifunctional() { 2 } else { 1 };

        for charge in charges {
            let mass = mz * charge as f32;
            let max_fc = frag_charge_bound(sc.max_fragment_charge, charge);
            // Stub-peptide ions are charge-reduced precursors
            let peaks = neutral_candidates(query, charge.saturating_sub(1).max(1));
            let has_peak = |target: f32| -> bool {
                let (lo, hi) = sc.fragment_tol.bounds(target);
                let i = peaks.partition_point(|&x| x < lo);
                i < peaks.len() && peaks[i] <= hi
            };

            // each peptide pair once per precursor charge: the best-scoring way round
            let mut seen: HashMap<(u32, u32), usize> = HashMap::new();
            let mut loops_seen: HashSet<usize> = HashSet::new();
            for iso in sc.min_isotope_err..=sc.max_isotope_err {
                let m = mass - iso as f32 * NEUTRON;

                // Loop-links: one peptide plus the linker explains the precursor;
                // the opened ring shows as stub-carrying fragment ions
                if cfg.loop_links {
                    let (lo, hi) = sc.precursor_tol.bounds(m);
                    let start = db.peptides.partition_point(|p| p.monoisotopic < lo - linker);
                    let end = db.peptides[start..]
                        .iter()
                        .take_while(|p| p.monoisotopic <= hi - linker)
                        .take(MAX_PARTNERS_PER_SEED)
                        .count();
                    for ix in start..start + end {
                        let p = &db.peptides[ix];
                        if !loops_seen.insert(ix) {
                            continue;
                        }
                        pairs_scored += 1;
                        let Some((score, first, second, tally, ring)) =
                            self.score_cleavable_loop(query, p, &stubs, max_fc)
                        else {
                            continue;
                        };
                        let calc = p.monoisotopic + linker;
                        found.push(CrosslinkMatch {
                            psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
                            spec_id: query.id.clone(),
                            file_id: query.file_id,
                            kind: LinkKind::Loop,
                            class: PairClass::of(p.decoy, p.decoy),
                            alpha: PeptideIx(ix as u32),
                            beta: None,
                            alpha_site: first as u16,
                            beta_site: second as u16,
                            intra_protein: true,
                            charge,
                            rt: query.scan_start_time,
                            expmass: mass,
                            calcmass: calc,
                            precursor_ppm: ppm_delta(mass, iso, calc),
                            isotope_error: iso as f32 * NEUTRON,
                            hyperscore: score,
                            alpha_score: score,
                            alpha_matched: tally.matched(),
                            alpha_linked: ring,
                            ..Default::default()
                        });
                    }
                }

                for &(arm_a, arm_b) in &stubs.arm_pairs() {
                let (x, y) = (&arms[arm_a], &arms[arm_b]);
                for &a in &peaks {
                    // A's complete doublet on its arm
                    if !has_peak(a + x.delta()) {
                        continue;
                    }
                    let mass_a = a - x.short;
                    let mass_b = m - mass_a - linker;
                    // B's doublet on the other arm
                    let comp_long = has_peak(mass_b + y.long);
                    let comp_short = has_peak(mass_b + y.short);
                    let complement = comp_long as u8 + comp_short as u8;
                    let needed = if cfg.require_complete_doublets { 2 } else { 1 };
                    if complement < needed {
                        continue;
                    }
                    if mass_a < min_peptide_mass || mass_b < min_peptide_mass {
                        continue;
                    }
                    // (peptide, linker end, site, tally, chain score)
                    let mut side = |mass: f32, arm: usize| -> Vec<(u32, usize, usize, Tally, f64)> {
                        let mut v = Vec::new();
                        for ix in self.peptides_with_mass(mass) {
                            let entry = chain_cache.entry((ix as u32, arm)).or_insert_with(|| {
                                self.score_stubbed_chain(query, &db.peptides[ix], &arms[arm], max_fc)
                            });
                            for (end, best) in entry.iter().enumerate().take(ends) {
                                if let Some((site, tally)) = *best {
                                    v.push((ix as u32, end, site, tally, tally.score(sc.score_type)));
                                }
                            }
                        }
                        v.sort_by(|a, b| b.4.total_cmp(&a.4));
                        v.truncate(ends * MAX_CHAINS_PER_SIDE.min(cfg.candidates));
                        v
                    };
                    let side_a = side(mass_a, arm_a);
                    let side_b = side(mass_b, arm_b);

                    // Best end assignment per peptide pair: the two peptides must
                    // carry different ends of a heterobifunctional linker
                    let mut pairs: Vec<(u32, usize, Tally, f64, u32, usize, Tally, f64)> = Vec::new();
                    for &(ix_a, end_a, site_a, tally_a, score_a) in &side_a {
                        for &(ix_b, end_b, site_b, tally_b, score_b) in &side_b {
                            if ends == 2 && end_a == end_b {
                                continue;
                            }
                            let pair = (ix_a, site_a, tally_a, score_a, ix_b, site_b, tally_b, score_b);
                            let joint = |p: &(u32, usize, Tally, f64, u32, usize, Tally, f64)| {
                                p.2.combine(&p.6).score(sc.score_type)
                            };
                            match pairs.iter_mut().find(|p| p.0 == ix_a && p.4 == ix_b) {
                                Some(p) if joint(&pair) > joint(p) => *p = pair,
                                Some(_) => {}
                                None => pairs.push(pair),
                            }
                        }
                    }

                    for (ix_a, site_a, tally_a, score_a, ix_b, site_b, tally_b, score_b) in pairs {
                        let key = (ix_a.min(ix_b), ix_a.max(ix_b));
                        let joint = tally_a.combine(&tally_b).score(sc.score_type);
                        let replace = match seen.get(&key) {
                            Some(&i) if found[i].hyperscore >= joint => continue,
                            Some(&i) => Some(i),
                            None => None,
                        };
                        if replace.is_none() {
                            pairs_scored += 1;
                        }
                        let (pa, pb) = (&db.peptides[ix_a as usize], &db.peptides[ix_b as usize]);
                        let calc = pa.monoisotopic + pb.monoisotopic + linker;
                        let swap = score_b > score_a;
                        let (alpha, beta, s_a, s_b, sc_a, sc_b, n_a, n_b) = if swap {
                            (ix_b, ix_a, site_b, site_a, score_b, score_a, tally_b, tally_a)
                        } else {
                            (ix_a, ix_b, site_a, site_b, score_a, score_b, tally_a, tally_b)
                        };
                        let csm = CrosslinkMatch {
                            psm_id: CSM_COUNTER.fetch_add(1, Ordering::Relaxed),
                            spec_id: query.id.clone(),
                            file_id: query.file_id,
                            kind: LinkKind::Inter,
                            class: PairClass::of(pa.decoy, pb.decoy),
                            alpha: PeptideIx(alpha),
                            beta: Some(PeptideIx(beta)),
                            alpha_site: s_a as u16,
                            beta_site: s_b as u16,
                            intra_protein: pa.proteins.iter().any(|p| pb.proteins.contains(p)),
                            charge,
                            rt: query.scan_start_time,
                            expmass: mass,
                            calcmass: calc,
                            precursor_ppm: ppm_delta(mass, iso, calc),
                            isotope_error: iso as f32 * NEUTRON,
                            hyperscore: n_a.combine(&n_b).score(sc.score_type),
                            alpha_score: sc_a,
                            beta_score: sc_b,
                            alpha_matched: n_a.matched(),
                            beta_matched: n_b.matched(),
                            stub_peaks: 2 + complement,
                            ..Default::default()
                        };
                        match replace {
                            Some(i) => found[i] = csm,
                            None => {
                                seen.insert(key, found.len());
                                found.push(csm);
                            }
                        }
                    }
                }
                }
            }
        }
        self.finalize(query, found, pairs_scored)
    }
}

/// Candidate *neutral* masses for every peak, sorted ascending.
///
/// Deisotoped peaks with an assigned charge already hold a neutral mass. Peaks
/// without one are stored on a z=1 scale, so (for multiply charged stub-peptide
/// ions) every charge up to `max_charge` is a hypothesis.
fn neutral_candidates(query: &ProcessedSpectrum, max_charge: u8) -> Vec<f32> {
    let mut out = Vec::with_capacity(query.masses.len() * max_charge as usize);
    for (idx, &mass) in query.masses.iter().enumerate() {
        match query.charges.get(idx).copied().unwrap_or(1) {
            z if z > 1 => out.push(mass),
            _ => out.extend((1..=max_charge).map(|z| mass * z as f32)),
        }
    }
    out.sort_by(f32::total_cmp);
    out
}

/// Build the chain kernel for the configured backend.
///
/// `Ok(None)` means "use the seed/partner algorithm" (the `cpu` backend, any
/// MS-cleavable reagent, or `auto` without a usable GPU).
pub fn make_kernel(
    settings: &CrosslinkSettings,
    db: &crate::database::IndexedDatabase,
) -> Result<Option<Box<dyn ChainKernel>>, String> {
    use crate::crosslink_kernel::{CpuKernel, IonTable};
    if settings.is_cleavable() {
        if settings.backend != Backend::Cpu {
            log::info!(
                "crosslink backend `{:?}` does not apply to MS-cleavable reagents; \
                 using the CPU doublet search",
                settings.backend
            );
        }
        return Ok(None);
    }
    let table = || IonTable::build(db, |p| settings.link_sites(p));
    match settings.backend {
        Backend::Cpu => Ok(None),
        Backend::Reflect => {
            let table = table();
            log::info!("crosslink backend: reflect (cpu), ion table {} MB", table.bytes() >> 20);
            Ok(Some(Box::new(CpuKernel::new(table))))
        }
        Backend::Gpu | Backend::Auto => {
            let explicit = settings.backend == Backend::Gpu;
            #[cfg(feature = "gpu")]
            {
                match crate::crosslink_gpu::GpuKernel::new(table(), explicit) {
                    Ok(k) => return Ok(Some(Box::new(k))),
                    Err(e) if explicit => return Err(format!("crosslink backend `gpu`: {e}")),
                    Err(e) => log::info!("no usable GPU ({e}); using the cpu backend"),
                }
            }
            #[cfg(not(feature = "gpu"))]
            {
                if explicit {
                    return Err("crosslink backend `gpu` requires building Sage with \
                                `--features gpu`"
                        .into());
                }
                log::info!("built without the `gpu` feature; using the cpu backend");
            }
            Ok(None)
        }
    }
}

fn ppm_delta(observed: f32, isotope_error: i8, calculated: f32) -> f32 {
    let observed = observed - isotope_error as f32 * NEUTRON;
    (observed - calculated) * 2E6 / (observed + calculated)
}

// ---------------------------------------------------------------------------
// FDR
// ---------------------------------------------------------------------------

/// q-values for items sorted by score (best first).
///
/// Inter-peptide links use the standard target/decoy crosslink estimate
/// `(TD - DD) / TT` (Walzthoeni et al., Nat. Methods 2012): a random false pair
/// is TD or DD, and the TD count over-counts by the DD count. Loop-links are
/// single peptides, so they use the ordinary `decoys / targets`.
fn q_from_sorted(classes: &[PairClass], loops: bool) -> Vec<f32> {
    let (mut tt, mut td, mut dd) = (0f64, 0f64, 0f64);
    let mut q = Vec::with_capacity(classes.len());
    for class in classes {
        match class {
            PairClass::TT => tt += 1.0,
            PairClass::TD => td += 1.0,
            PairClass::DD => dd += 1.0,
        }
        let false_hits = if loops { td + dd } else { td - dd };
        q.push((false_hits / tt.max(1.0)).clamp(0.0, 1.0) as f32);
    }
    // q-values are monotone: the best FDR at which a hit would be accepted
    let mut running = 1.0f32;
    for v in q.iter_mut().rev() {
        running = running.min(*v);
        *v = running;
    }
    q
}

/// q-values of multivalent species sorted by score (best first), from
/// `(any decoy peptide, number of distinct peptides)` of each.
///
/// A random species of `n` distinct peptides is all-target with chance
/// `1 / 2^n`, so each species with a decoy stands for `1 / (2^n - 1)` false
/// targets. For pairs that is `(TD + DD) / 3`, which `TD - DD` also estimates;
/// unlike `N1 - N2 + …` it never turns negative, and copies of one peptide
/// (`a2a`: one distinct peptide, decoy or not) count like a single peptide.
fn q_multivalent(species: &[(bool, u8)]) -> Vec<f32> {
    let (mut targets, mut false_hits) = (0f64, 0f64);
    let mut q = Vec::with_capacity(species.len());
    for &(decoy, distinct) in species {
        if decoy {
            false_hits += 1.0 / ((1u32 << distinct.clamp(1, 16)) - 1) as f64;
        } else {
            targets += 1.0;
        }
        q.push((false_hits / targets.max(1.0)).clamp(0.0, 1.0) as f32);
    }
    let mut running = 1.0f32;
    for v in q.iter_mut().rev() {
        running = running.min(*v);
        *v = running;
    }
    q
}

/// FDR groups: loop-links, intra-protein links, inter-protein links,
/// multivalent species. Intra-protein links are far more likely to be true
/// than inter-protein links, so pooling them would hide inter-protein false
/// positives (Fischer & Rappsilber, Anal. Chem. 2017).
fn fdr_group(m: &CrosslinkMatch) -> usize {
    match (m.kind, m.intra_protein) {
        (LinkKind::Loop, _) => 0,
        (LinkKind::Inter, true) => 1,
        (LinkKind::Inter, false) => 2,
        (LinkKind::Multi, _) => 3,
    }
}

/// Canonical identity of a unique crosslink (peptide pair + link sites),
/// independent of which chain is called alpha; multivalent species add a
/// hash of all their peptides and links
fn link_key(m: &CrosslinkMatch) -> (u32, u16, u32, u16, u64) {
    if let Some(s) = &m.species {
        use std::hash::{Hash, Hasher};
        // links by (peptide, residue) at both ends, ends and links in order
        let mut links: Vec<((u32, u16), (u32, u16))> = s
            .links
            .iter()
            .map(|&(a, sa, b, sb)| {
                let x = (s.peptides[a as usize].0, sa);
                let y = (s.peptides[b as usize].0, sb);
                if x <= y { (x, y) } else { (y, x) }
            })
            .collect();
        links.sort_unstable();
        let mut h = std::collections::hash_map::DefaultHasher::new();
        links.hash(&mut h);
        return (u32::MAX, 0, u32::MAX, 0, h.finish());
    }
    match m.beta {
        None => (m.alpha.0, m.alpha_site, u32::MAX, m.beta_site, 0),
        Some(b) => {
            let x = (m.alpha.0, m.alpha_site);
            let y = (b.0, m.beta_site);
            let (lo, hi) = if x <= y { (x, y) } else { (y, x) };
            (lo.0, lo.1, hi.0, hi.1, 0)
        }
    }
}

/// Assign `spectrum_q` (per CSM) and `crosslink_q` (per unique peptide pair +
/// sites, scored by its best CSM), separately within each FDR group.
///
/// Returns `(target CSMs, unique target crosslinks)` at 1% FDR.
pub fn assign_q_values(matches: &mut [CrosslinkMatch]) -> (usize, usize) {
    let (mut csm_count, mut link_count) = (0, 0);

    for group in 0..4 {
        let loops = group == 0;
        // q-values of matches sorted best first, by the group's estimate
        let q_of = |sorted: &[usize], matches: &[CrosslinkMatch]| {
            if group == 3 {
                let species: Vec<(bool, u8)> = sorted
                    .iter()
                    .map(|&i| match &matches[i].species {
                        Some(s) => {
                            let mut distinct = s.peptides.clone();
                            distinct.sort_unstable();
                            distinct.dedup();
                            (s.decoys > 0, distinct.len() as u8)
                        }
                        None => (matches[i].class != PairClass::TT, 2),
                    })
                    .collect();
                q_multivalent(&species)
            } else {
                let classes: Vec<PairClass> = sorted.iter().map(|&i| matches[i].class).collect();
                q_from_sorted(&classes, loops)
            }
        };

        // CSM level
        let mut idx = (0..matches.len())
            .filter(|&i| fdr_group(&matches[i]) == group)
            .collect::<Vec<_>>();
        idx.sort_by(|&a, &b| matches[b].hyperscore.total_cmp(&matches[a].hyperscore));
        for (&i, q) in idx.iter().zip(q_of(&idx, matches)) {
            matches[i].spectrum_q = q;
            if q <= 0.01 && matches[i].class == PairClass::TT {
                csm_count += 1;
            }
        }

        // Unique-crosslink level: best CSM per link
        let mut best: HashMap<(u32, u16, u32, u16, u64), usize> = HashMap::new();
        for &i in &idx {
            best.entry(link_key(&matches[i]))
                .and_modify(|j| {
                    if matches[i].hyperscore > matches[*j].hyperscore {
                        *j = i
                    }
                })
                .or_insert(i);
        }
        let mut reps = best.values().copied().collect::<Vec<_>>();
        reps.sort_by(|&a, &b| {
            matches[b]
                .hyperscore
                .total_cmp(&matches[a].hyperscore)
                .then(link_key(&matches[a]).cmp(&link_key(&matches[b])))
        });
        let mut link_q = HashMap::new();
        for (&i, q) in reps.iter().zip(q_of(&reps, matches)) {
            link_q.insert(link_key(&matches[i]), q);
            if q <= 0.01 && matches[i].class == PairClass::TT {
                link_count += 1;
            }
        }
        for &i in &idx {
            matches[i].crosslink_q = link_q[&link_key(&matches[i])];
        }
    }
    (csm_count, link_count)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::{Builder, IndexedDatabase, Parameters};
    use crate::enzyme::Digest;
    use crate::mass::H2O;
    use crate::spectrum::Precursor;

    const DSS: f32 = 138.06808;

    fn peptide(seq: &str, protein: &str, position: Position) -> Peptide {
        Peptide::try_from(Digest {
            sequence: seq.into(),
            protein: protein.into(),
            position,
            ..Default::default()
        })
        .unwrap()
    }

    /// Deterministic pseudo-random tryptic-looking peptides
    fn distractors(n: usize) -> Vec<Peptide> {
        let alphabet = b"ADEFGHILMNPQSTVWYKR";
        let mut state = 0x2545F4914F6CDD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..n)
            .map(|i| {
                let len = 8 + (next() % 12) as usize;
                let mut s = (0..len)
                    .map(|_| alphabet[(next() % alphabet.len() as u64) as usize])
                    .collect::<Vec<_>>();
                s.push(b'R');
                peptide(
                    std::str::from_utf8(&s).unwrap(),
                    &format!("decoyprot{i}"),
                    Position::Internal,
                )
            })
            .collect()
    }

    fn database(targets: Vec<Peptide>) -> IndexedDatabase {
        let mut all = Vec::new();
        for t in targets {
            all.push(t.reverse());
            all.push(t);
        }
        Parameters::reorder_peptides(&mut all);
        Builder {
            fasta: Some("unused.fasta".into()),
            ..Default::default()
        }
        .make_parameters()
        .unwrap()
        .build_from_peptides(all)
    }

    fn scorer(db: &IndexedDatabase) -> Scorer<'_> {
        Scorer {
            db,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            min_matched_peaks: 4,
            min_isotope_err: 0,
            max_isotope_err: 0,
            min_precursor_charge: 2,
            max_precursor_charge: 4,
            override_precursor_charge: false,
            max_fragment_charge: None,
            chimera: false,
            report_psms: 1,
            wide_window: false,
            annotate_matches: false,
            score_type: ScoreType::SageHyperScore,
        }
    }

    /// Build an MS2 spectrum from theoretical ions (charges 1 and 2) plus noise
    fn spectrum(neutral_mass: f32, charge: u8, ions: &[(Kind, f32)]) -> ProcessedSpectrum {
        let mut peaks: Vec<(f32, f32)> = Vec::new();
        for (i, &(_, m)) in ions.iter().enumerate() {
            for z in 1..=2 {
                peaks.push((m / z as f32, 100.0 + (i % 7) as f32 * 10.0));
            }
        }
        // deterministic noise peaks
        for k in 0..60 {
            peaks.push((120.31 + k as f32 * 17.137, 15.0));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        ProcessedSpectrum {
            level: 2,
            id: "scan=1".into(),
            precursors: vec![Precursor {
                mz: neutral_mass / charge as f32 + PROTON,
                charge: Some(charge),
                ..Default::default()
            }],
            total_ion_current: peaks.iter().map(|p| p.1).sum(),
            masses: peaks.iter().map(|p| p.0).collect(),
            intensities: peaks.iter().map(|p| p.1).collect(),
            charges: vec![1; peaks.len()],
            ..Default::default()
        }
    }

    #[test]
    fn ion_masses_are_complementary() {
        // For any link site, b_k + y_(n-k) must equal M(chain) + shift: the
        // shift is carried by exactly one of the two complementary ions.
        let pep = peptide("GLEKVSTAAMQR", "P1", Position::Internal);
        let mods = Modifications::default();
        let shift = 1234.5;
        for site in [0usize, 3, 7, 10] {
            let b = chain_ions(&pep, &mods, &[Kind::B], Link::Single { site, shift });
            let y = chain_ions(&pep, &mods, &[Kind::Y], Link::Single { site, shift });
            let n = pep.sequence.len();
            assert_eq!(b.len(), n - 1);
            for k in 1..n {
                let bk = b[k - 1].1;
                let y_comp = y[k - 1].1; // the y ion holding residues k..n
                assert!(
                    (bk + y_comp - (pep.monoisotopic + shift)).abs() < 1e-3,
                    "site {site} k {k}: {} vs {}",
                    bk + y_comp,
                    pep.monoisotopic + shift
                );
            }
        }
    }

    #[test]
    fn unshifted_ions_match_linear_series() {
        let pep = peptide("GLEKVSTAAMQR", "P1", Position::Internal);
        let mods = Modifications::default();
        // link at residue 3: b ions from idx 3 on and y ions before idx 3 carry the shift
        let shifted = chain_ions(&pep, &mods, &[Kind::B, Kind::Y], Link::Single { site: 3, shift: 500.0 });
        let linear = IonSeries::new(&pep, Kind::B, &mods)
            .chain(IonSeries::new(&pep, Kind::Y, &mods))
            .collect::<Vec<_>>();
        let n = pep.sequence.len() - 1;
        // b ions before the site (idx 0..3) and y ions after it (idx 3..) are unchanged
        for idx in 0..3 {
            assert_eq!(shifted[idx].1, linear[idx].monoisotopic_mass);
        }
        for idx in 3..n {
            assert!((shifted[n + idx].1 - linear[n + idx].monoisotopic_mass).abs() < 1e-4);
        }
        for idx in 3..n {
            assert!((shifted[idx].1 - linear[idx].monoisotopic_mass - 500.0).abs() < 1e-3);
        }
    }

    #[test]
    fn loop_link_ions_skip_the_ring() {
        let pep = peptide("AEKLMNPQSTVKGR", "P1", Position::Internal);
        let mods = Modifications::default();
        let ions = chain_ions(
            &pep,
            &mods,
            &[Kind::B],
            Link::Loop {
                first: 2,
                second: 11,
                shift: DSS,
            },
        );
        // b ions with idx in [2, 11) lie inside the ring and are unobservable
        assert_eq!(ions.len(), (pep.sequence.len() - 1) - 9);
    }

    #[test]
    fn link_site_rules() {
        let cfg = CrosslinkSettings::default();
        // internal K only; C-terminal K excluded
        let p = peptide("AKCDKEFK", "P", Position::Internal);
        assert_eq!(cfg.link_sites(&p), vec![1, 4]);
        // protein C-terminus: terminal K allowed
        let p = peptide("AKCDKEFK", "P", Position::Cterm);
        assert_eq!(cfg.link_sites(&p), vec![1, 4, 7]);
        // protein N-terminus adds the amine at residue 0
        let p = peptide("ACDEKFGR", "P", Position::Nterm);
        assert_eq!(cfg.link_sites(&p), vec![0, 4]);
        // N-terminal K isn't listed twice
        let p = peptide("KCDEKFGR", "P", Position::Nterm);
        assert_eq!(cfg.link_sites(&p), vec![0, 4]);
        // a modified lysine is blocked
        let mut p = peptide("AKCDKEFR", "P", Position::Internal);
        let mut mods = Modifications::default();
        let id = mods.insert(42.0106, None).unwrap();
        p.modifications[1] = Some(id);
        assert_eq!(cfg.link_sites(&p), vec![4]);
        // N-terminal amine disabled
        let cfg = CrosslinkSettings {
            protein_nterm: false,
            ..Default::default()
        };
        let p = peptide("ACDEKFGR", "P", Position::Nterm);
        assert_eq!(cfg.link_sites(&p), vec![4]);
    }

    #[test]
    fn settings_deserialize_and_validate() {
        let cfg: CrosslinkSettings = serde_json::from_str("{}").unwrap();
        assert_eq!(cfg, CrosslinkSettings::default());
        assert!(cfg.validate().is_ok());

        let cfg: CrosslinkSettings =
            serde_json::from_str(r#"{"mass": 96.02113, "residues": "KST"}"#).unwrap();
        assert_eq!(cfg.linker_mass(), 96.02113);
        assert!(serde_json::from_str::<CrosslinkSettings>(r#"{"typo": 1}"#).is_err());

        let bad = CrosslinkSettings {
            mass: Some(-50.0),
            ..Default::default()
        };
        assert!(bad.validate().is_err());
        // zero-length linkers remove a water
        let zero_length = CrosslinkSettings {
            mass: Some(-18.010565),
            ..Default::default()
        };
        assert!(zero_length.validate().is_ok());
    }

    fn edc() -> CrosslinkSettings {
        CrosslinkSettings {
            reagent: Some(Reagent::Edc),
            ..Default::default()
        }
    }

    /// Spectrum of `a` linked at `sa` to `b` linked at `sb` by a non-cleavable linker
    fn linked_spectrum(db: &IndexedDatabase, a: &Peptide, b: &Peptide, sa: usize, sb: usize, linker: f32) -> ProcessedSpectrum {
        let kinds = [Kind::B, Kind::Y];
        let mut ions = chain_ions(a, &db.modifications, &kinds, Link::Single { site: sa, shift: b.monoisotopic + linker });
        ions.extend(chain_ions(b, &db.modifications, &kinds, Link::Single { site: sb, shift: a.monoisotopic + linker }));
        spectrum(a.monoisotopic + b.monoisotopic + linker, 3, &ions)
    }

    #[test]
    fn heterobifunctional_site_roles() {
        let cfg = edc();
        assert!(cfg.validate().is_ok());
        assert!(cfg.is_heterobifunctional());
        assert!((cfg.linker_mass() + 18.010565).abs() < 1e-5);

        // internal peptide: K is site1, D/E are site2; the C-terminal K is excluded
        let p = peptide("ADKEFK", "P", Position::Internal);
        assert_eq!(cfg.link_site_roles(&p), vec![(1, SITE2), (2, SITE1), (3, SITE2)]);
        // protein N-terminus (amine, site1) on residue 0, which is also an Asp (site2)
        let p = peptide("DAKFGR", "P", Position::Nterm);
        assert_eq!(cfg.link_site_roles(&p), vec![(0, SITE1 | SITE2), (2, SITE1)]);
        // protein C-terminal carboxyl (site2) on the last residue
        let p = peptide("AKFGR", "P", Position::Cterm);
        assert_eq!(cfg.link_site_roles(&p), vec![(1, SITE1), (4, SITE2)]);

        assert!(pairable(SITE1, SITE2) && pairable(SITE2, SITE1));
        assert!(!pairable(SITE1, SITE1) && !pairable(SITE2, SITE2));
        assert!(pairable(SITE1 | SITE2, SITE1));

        // explicit sites override the preset; both are needed
        let custom: CrosslinkSettings =
            serde_json::from_str(r#"{"mass": 50.0, "site1": "C", "site2": "K["}"#).unwrap();
        assert!(custom.validate().is_ok());
        assert_eq!(custom.hetero_sites(), Some(("C".into(), "K[".into())));
        let half: CrosslinkSettings = serde_json::from_str(r#"{"site1": "K"}"#).unwrap();
        assert!(half.validate().is_err());
        let bad: CrosslinkSettings =
            serde_json::from_str(r#"{"site1": "K", "site2": "d"}"#).unwrap();
        assert!(bad.validate().is_err());
        // homobifunctional reagents: every site takes either end
        let p = peptide("AKCDKEFK", "P", Position::Internal);
        assert_eq!(
            CrosslinkSettings::default().link_site_roles(&p),
            vec![(1, SITE1 | SITE2), (4, SITE1 | SITE2)]
        );
    }

    #[test]
    fn finds_zero_length_link_between_lys_and_glu() {
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let cfg = edc();
        // GLE[K]VSTAAMQR (K, site1) to TDYHFPN[E]WKLIR (E, site2)
        let query = linked_spectrum(&db, &a, &b, 3, 7, cfg.linker_mass());
        for kernel in [None, Some(kernel_for(&db, &cfg))] {
            let xl = CrosslinkScorer::new(&sc, &cfg).with_kernel(kernel.as_ref().map(|k| k as &dyn ChainKernel));
            let hits = xl.score(&query);
            assert_eq!(hits.len(), 1, "{hits:?}");
            assert!(hits[0].precursor_ppm.abs() < 1.0);
            assert_eq!(
                pair_set(&db, &hits[0]),
                vec![("GLEKVSTAAMQR".to_string(), 3), ("TDYHFPNEWKLIR".to_string(), 7)]
            );
        }
    }

    #[test]
    fn heterobifunctional_never_links_like_groups() {
        // A Lys-Lys spectrum searched with EDC chemistry: whatever is reported
        // must still join an amine to a carboxyl
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let cfg = edc();
        let query = linked_spectrum(&db, &a, &b, 3, 9, cfg.linker_mass());
        for hit in CrosslinkScorer::new(&sc, &cfg).score(&query) {
            let roles = |ix: PeptideIx, site: u16| {
                cfg.link_site_roles(&db[ix])
                    .into_iter()
                    .find(|&(s, _)| s == site as usize)
                    .map(|(_, r)| r)
                    .unwrap()
            };
            let ra = roles(hit.alpha, hit.alpha_site);
            let rb = roles(hit.beta.unwrap(), hit.beta_site);
            assert!(pairable(ra, rb), "{hit:?}");
        }
    }

    #[test]
    fn finds_heterobifunctional_cleavable_link() {
        // A hypothetical DSSO-like linker reacting with K at one end and D/E at the other
        let cfg: CrosslinkSettings = serde_json::from_str(
            r#"{"mass": 158.00376, "cleavable": {"short": 54.01056, "long": 103.9932},
                "site1": "K", "site2": "DE"}"#,
        )
        .unwrap();
        assert!(cfg.validate().is_ok());
        let stubs = cfg.stubs().unwrap();
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let xl = CrosslinkScorer::new(&sc, &cfg);

        let mut ions = stubbed_ions(&a, &db.modifications, 3, &stubs);
        ions.extend(stubbed_ions(&b, &db.modifications, 7, &stubs));
        for p in [&a, &b] {
            ions.push((Kind::B, p.monoisotopic + stubs.short));
            ions.push((Kind::B, p.monoisotopic + stubs.long));
        }
        let query = spectrum(a.monoisotopic + b.monoisotopic + cfg.linker_mass(), 4, &ions);
        let hits = xl.score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(
            pair_set(&db, &hits[0]),
            vec![("GLEKVSTAAMQR".to_string(), 3), ("TDYHFPNEWKLIR".to_string(), 7)]
        );
    }

    #[test]
    fn finds_inter_peptide_link() {
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let (site_a, site_b) = (3usize, 9usize);

        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let xl = CrosslinkScorer::new(&sc, &cfg);

        let kinds = [Kind::B, Kind::Y];
        let mut ions = chain_ions(
            &a,
            &db.modifications,
            &kinds,
            Link::Single {
                site: site_a,
                shift: b.monoisotopic + DSS,
            },
        );
        ions.extend(chain_ions(
            &b,
            &db.modifications,
            &kinds,
            Link::Single {
                site: site_b,
                shift: a.monoisotopic + DSS,
            },
        ));
        let precursor = a.monoisotopic + b.monoisotopic + DSS;
        let query = spectrum(precursor, 3, &ions);

        let hits = xl.score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        let hit = &hits[0];
        assert_eq!(hit.kind, LinkKind::Inter);
        assert_eq!(hit.class, PairClass::TT);
        assert!(!hit.intra_protein);
        assert!(hit.precursor_ppm.abs() < 1.0, "ppm {}", hit.precursor_ppm);
        assert!((hit.calcmass - precursor).abs() < 1e-3);

        // Alpha/beta ordering is by chain score, so check the pair as a set
        let seq = |ix: PeptideIx| String::from_utf8(db[ix].sequence.to_vec()).unwrap();
        let mut found = vec![
            (seq(hit.alpha), hit.alpha_site as usize),
            (seq(hit.beta.unwrap()), hit.beta_site as usize),
        ];
        found.sort();
        let mut expected = vec![("GLEKVSTAAMQR".to_string(), site_a), ("TDYHFPNEWKLIR".to_string(), site_b)];
        expected.sort();
        assert_eq!(found, expected);
        assert!(hit.alpha_matched >= 3 && hit.beta_matched >= 3);
    }

    #[test]
    fn finds_loop_link() {
        // Loop-links are found through the unshifted ions *outside* the ring, so
        // the linked residues need flanking sequence (here K10 .. K19, 1-based)
        let l = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_L", Position::Internal);
        let mut targets = distractors(400);
        targets.push(l.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let xl = CrosslinkScorer::new(&sc, &cfg);

        let ions = chain_ions(
            &l,
            &db.modifications,
            &[Kind::B, Kind::Y],
            Link::Loop {
                first: 9,
                second: 18,
                shift: DSS,
            },
        );
        let query = spectrum(l.monoisotopic + DSS, 3, &ions);
        let hits = xl.score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].kind, LinkKind::Loop);
        assert_eq!((hits[0].alpha_site, hits[0].beta_site), (9, 18));
        assert!(hits[0].intra_protein);
        assert!(hits[0].precursor_ppm.abs() < 1.0);
        assert!(hits[0].alpha_linked >= cfg.min_linked_peaks, "{:?}", hits[0]);
    }

    #[test]
    fn loop_link_needs_linker_carrying_ions() {
        // The right peptide with some other modification of about the linker's
        // mass: its spectrum has the ions outside the ring but none that carry
        // the linker, so only the precursor mass says "loop-link"
        let l = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_L", Position::Internal);
        let mut targets = distractors(400);
        targets.push(l.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let link = Link::Loop { first: 9, second: 18, shift: DSS };
        let all = chain_ions(&l, &db.modifications, &[Kind::B, Kind::Y], link);
        let linked = linked_ions(&l, &db.modifications, &[Kind::B, Kind::Y], link);
        // ions that span the whole ring: b19 and longer, and the y ions that
        // still contain K10
        assert_eq!(linked.len(), (l.sequence.len() - 1 - 18) + 9);
        let unlinked = all.iter().filter(|ion| !linked.contains(ion)).copied().collect::<Vec<_>>();
        let query = spectrum(l.monoisotopic + DSS, 3, &unlinked);

        let cfg = CrosslinkSettings::default();
        assert!(CrosslinkScorer::new(&sc, &cfg).score(&query).is_empty());
        // the check can be switched off
        let off = CrosslinkSettings { min_linked_peaks: 0, ..Default::default() };
        let hits = CrosslinkScorer::new(&sc, &off).score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].alpha_linked, 0);
    }

    #[test]
    fn linear_peptide_spectrum_yields_no_crosslink() {
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let mut targets = distractors(400);
        targets.push(a.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let xl = CrosslinkScorer::new(&sc, &cfg);

        // An ordinary linear spectrum: precursor is just the peptide
        let ions = IonSeries::new(&a, Kind::B, &db.modifications)
            .chain(IonSeries::new(&a, Kind::Y, &db.modifications))
            .map(|i| (i.kind, i.monoisotopic_mass))
            .collect::<Vec<_>>();
        let query = spectrum(a.monoisotopic, 2, &ions);
        assert!(xl.score(&query).is_empty());
    }

    #[test]
    fn precursor_mass_matches_formula() {
        // M = m(A) + m(B) + linker, where each m already includes one water
        let a = peptide("GLEKVSTAAMQR", "A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "B", Position::Internal);
        let residues = |p: &Peptide| p.sequence.iter().map(|&c| crate::mass::monoisotopic(c)).sum::<f32>();
        let expected = residues(&a) + residues(&b) + 2.0 * H2O + DSS;
        assert!((a.monoisotopic + b.monoisotopic + DSS - expected).abs() < 1e-3);
    }

    fn dsso() -> CrosslinkSettings {
        CrosslinkSettings {
            reagent: Some(Reagent::Dsso),
            ..Default::default()
        }
    }

    /// Ions of one DSSO-cleaved chain: linear ions outside the linked residue,
    /// and stub-bearing ions (alternating short/long stub) containing it
    fn stubbed_ions(p: &Peptide, mods: &Modifications, site: usize, stubs: &CleavableStubs) -> Vec<(Kind, f32)> {
        let kinds = [Kind::B, Kind::Y];
        let plain = chain_ions(p, mods, &kinds, Link::Single { site, shift: 0.0 });
        let short = chain_ions(p, mods, &kinds, Link::Single { site, shift: stubs.short });
        let long = chain_ions(p, mods, &kinds, Link::Single { site, shift: stubs.long });
        plain
            .iter()
            .enumerate()
            .map(|(i, &(kind, m))| {
                if short[i].1 == m {
                    (kind, m)
                } else if i % 2 == 0 {
                    short[i]
                } else {
                    long[i]
                }
            })
            .collect()
    }

    #[test]
    fn reagent_presets_and_validation() {
        let cfg: CrosslinkSettings = serde_json::from_str(r#"{"reagent": "dsso"}"#).unwrap();
        assert!(cfg.validate().is_ok());
        assert!(cfg.is_cleavable());
        assert!((cfg.linker_mass() - 158.00376).abs() < 1e-4);
        let stubs = cfg.stubs().unwrap();
        assert!((stubs.long - stubs.short - 49.98264).abs() < 1e-3);

        let cfg: CrosslinkSettings = serde_json::from_str(r#"{"reagent": "dsbu"}"#).unwrap();
        assert!(cfg.validate().is_ok());
        let stubs = cfg.stubs().unwrap();
        assert!((stubs.long - stubs.short - 25.97927).abs() < 1e-3);
        assert!((stubs.short + stubs.long - cfg.linker_mass()).abs() < 1e-3);

        // non-cleavable presets and the default
        let cfg: CrosslinkSettings = serde_json::from_str(r#"{"reagent": "bs3"}"#).unwrap();
        assert!(!cfg.is_cleavable());
        assert!(!CrosslinkSettings::default().is_cleavable());

        // a preset can't be mixed with custom chemistry
        let both: CrosslinkSettings =
            serde_json::from_str(r#"{"reagent": "dsso", "mass": 100.0}"#).unwrap();
        assert!(both.validate().is_err());

        // a cleavable preset's stubs can be overridden, a non-cleavable one has none
        let dsso: CrosslinkSettings = serde_json::from_str(
            r#"{"reagent": "dsso", "cleavable": {"short": 54.01056, "long": 103.9932, "alt": []}}"#,
        )
        .unwrap();
        assert!(dsso.validate().is_ok());
        assert!(dsso.stubs().unwrap().alt.is_empty());
        let bs3: CrosslinkSettings = serde_json::from_str(
            r#"{"reagent": "bs3", "cleavable": {"short": 54.01056, "long": 84.0575}}"#,
        )
        .unwrap();
        assert!(bs3.validate().is_err());

        // custom cleavable: stubs must add up to the linker mass
        let bad: CrosslinkSettings = serde_json::from_str(
            r#"{"mass": 158.0, "cleavable": {"short": 54.0, "long": 90.0}}"#,
        )
        .unwrap();
        assert!(bad.validate().unwrap_err().contains("must equal"));
        let good: CrosslinkSettings = serde_json::from_str(
            r#"{"mass": 158.00376, "cleavable": {"short": 54.01056, "long": 103.9932}}"#,
        )
        .unwrap();
        assert!(good.validate().is_ok());
    }

    /// A confident target crosslink match of `a` and `b`, as `assign_q_values` leaves it
    fn confident_csm(db: &IndexedDatabase, a: &str, b: &str, file_id: usize, rt: f32) -> CrosslinkMatch {
        let ix = |seq: &str| {
            PeptideIx(db.peptides.iter().position(|p| p.sequence.as_ref() == seq.as_bytes()).unwrap() as u32)
        };
        let (a, b) = (ix(a), ix(b));
        CrosslinkMatch {
            kind: LinkKind::Inter,
            class: PairClass::TT,
            alpha: a,
            beta: Some(b),
            alpha_site: 3,
            beta_site: 9,
            file_id,
            rt,
            calcmass: db[a].monoisotopic + db[b].monoisotopic + DSS,
            hyperscore: 30.0,
            rank: 1,
            ..Default::default()
        }
    }

    /// MS1 scans around `apex` with a Gaussian elution profile of the species'
    /// isotope envelope at charge `z`, scaled by `height`
    fn ms1_run(file_id: usize, mass: f32, z: u8, envelope: [f32; 3], apex: f32, height: f32) -> Vec<ProcessedSpectrum> {
        (0..=50)
            .map(|i| {
                let rt = apex - 0.5 + i as f32 * 0.02;
                let profile = height * (-0.5 * ((rt - apex) / 0.1).powi(2)).exp();
                // the species plus an unrelated peak
                let mut peaks = (0..3)
                    .map(|iso| ((mass + iso as f32 * NEUTRON) / z as f32, envelope[iso] * profile))
                    .collect::<Vec<_>>();
                peaks.push((523.27, 1000.0));
                peaks.sort_by(|x, y| x.0.total_cmp(&y.0));
                ProcessedSpectrum {
                    level: 1,
                    id: format!("ms1_{i}"),
                    file_id,
                    scan_start_time: rt,
                    masses: peaks.iter().map(|p| p.0).collect(),
                    intensities: peaks.iter().map(|p| p.1).collect(),
                    ..Default::default()
                }
            })
            .chain(std::iter::once(ProcessedSpectrum {
                level: 1,
                file_id,
                scan_start_time: 100.0,
                ..Default::default()
            }))
            .collect()
    }

    #[test]
    fn crosslink_lfq_measures_relative_abundance() {
        use crate::crosslink_quant::quantify_lfq;
        use crate::lfq::LfqSettings;

        let (db, _, _) = pair_fixture();
        let cfg = CrosslinkSettings::default();
        let csm = confident_csm(&db, "GLEKVSTAAMQR", "TDYHFPNEWKLIR", 0, 50.0);
        let (pa, pb) = (&db[csm.alpha], &db[csm.beta.unwrap()]);
        let composition = pa
            .sequence
            .iter()
            .chain(pb.sequence.iter())
            .map(|r| crate::mass::composition(*r))
            .sum::<crate::mass::Composition>();
        let (lc, ls) = cfg.linker_composition();
        let envelope = crate::isotopes::peptide_isotopes(composition.carbon + lc, composition.sulfur + ls);

        // Run 1 holds three times as much of the crosslink as run 0
        let mut ms1 = ms1_run(0, csm.calcmass, 3, envelope, 50.0, 1e6);
        ms1.extend(ms1_run(1, csm.calcmass, 3, envelope, 50.0, 3e6));

        let found = quantify_lfq(&[csm], &cfg, &LfqSettings::default(), &db, &ms1, 2, (3, 3));
        let targets = found.iter().filter(|q| !q.decoy).collect::<Vec<_>>();
        assert_eq!(targets.len(), 1, "{found:?}");
        let areas = &targets[0].areas;
        assert!(areas[0] > 0.0);
        let ratio = areas[1] / areas[0];
        assert!((ratio - 3.0).abs() < 0.3, "ratio {ratio}, areas {areas:?}");
        assert!(targets[0].peak.spectral_angle > 0.9);
        assert_eq!(targets[0].csm, 0);

        // nothing passes the q-value cutoff: nothing to quantify
        let mut doubtful = confident_csm(&db, "GLEKVSTAAMQR", "TDYHFPNEWKLIR", 0, 50.0);
        doubtful.crosslink_q = 0.5;
        assert!(quantify_lfq(&[doubtful], &cfg, &LfqSettings::default(), &db, &ms1, 2, (3, 3)).is_empty());
    }

    #[test]
    fn crosslink_runs_are_aligned_on_shared_species() {
        use crate::crosslink_quant::align_runs;
        let (db, _, _) = pair_fixture();
        // run 1 elutes everything 10% later
        let mut matches = Vec::new();
        for (k, rt) in [20.0f32, 40.0, 60.0, 80.0].into_iter().enumerate() {
            let mut a = confident_csm(&db, "GLEKVSTAAMQR", "TDYHFPNEWKLIR", 0, rt);
            let mut b = confident_csm(&db, "GLEKVSTAAMQR", "TDYHFPNEWKLIR", 1, rt * 1.1);
            // distinct species per retention time
            a.alpha = PeptideIx(k as u32);
            b.alpha = PeptideIx(k as u32);
            matches.push(a);
            matches.push(b);
        }
        let species = (0..4).map(|k| vec![2 * k, 2 * k + 1]).collect::<Vec<_>>();
        let ms1 = vec![
            ProcessedSpectrum { file_id: 0, scan_start_time: 100.0, ..Default::default() },
            ProcessedSpectrum { file_id: 1, scan_start_time: 100.0, ..Default::default() },
        ];
        let al = align_runs(&matches, &species, &ms1, 2);
        let aligned = |f: usize, rt: f32| (rt / al[f].max_rt) * al[f].slope + al[f].intercept;
        for rt in [20.0f32, 40.0, 60.0, 80.0] {
            assert!((aligned(0, rt) - aligned(1, rt * 1.1)).abs() < 1e-3);
        }
        // a single shared species is too little to fit: identity
        let al = align_runs(&matches[..2], &[vec![0, 1]], &ms1, 2);
        assert_eq!((al[1].slope, al[1].intercept), (1.0, 0.0));
    }

    #[test]
    fn neutral_candidates_handle_unknown_charge() {
        let mut q = ProcessedSpectrum {
            masses: vec![500.0, 1000.0],
            intensities: vec![1.0, 1.0],
            charges: vec![1, 3],
            ..Default::default()
        };
        // charge-1 peak: hypotheses z=1,2; deisotoped 3+ peak already neutral
        assert_eq!(neutral_candidates(&q, 2), vec![500.0, 1000.0, 1000.0]);
        q.charges = vec![1, 1];
        assert_eq!(neutral_candidates(&q, 2), vec![500.0, 1000.0, 1000.0, 2000.0]);
    }

    #[test]
    fn finds_cleavable_crosslink() {
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let (site_a, site_b) = (3usize, 9usize);
        let cfg = dsso();
        let stubs = cfg.stubs().unwrap();

        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let xl = CrosslinkScorer::new(&sc, &cfg);

        let mut ions = stubbed_ions(&a, &db.modifications, site_a, &stubs);
        ions.extend(stubbed_ions(&b, &db.modifications, site_b, &stubs));
        // the four doublet peaks: each peptide with its short and long stub
        for p in [&a, &b] {
            ions.push((Kind::B, p.monoisotopic + stubs.short));
            ions.push((Kind::B, p.monoisotopic + stubs.long));
        }
        let precursor = a.monoisotopic + b.monoisotopic + cfg.linker_mass();
        let query = spectrum(precursor, 4, &ions);

        let hits = xl.score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        let hit = &hits[0];
        assert_eq!(hit.class, PairClass::TT);
        assert_eq!(hit.stub_peaks, 4);
        assert!(hit.precursor_ppm.abs() < 1.0, "ppm {}", hit.precursor_ppm);

        let seq = |ix: PeptideIx| String::from_utf8(db[ix].sequence.to_vec()).unwrap();
        let mut found = vec![
            (seq(hit.alpha), hit.alpha_site as usize),
            (seq(hit.beta.unwrap()), hit.beta_site as usize),
        ];
        found.sort();
        let mut expected = vec![
            ("GLEKVSTAAMQR".to_string(), site_a),
            ("TDYHFPNEWKLIR".to_string(), site_b),
        ];
        expected.sort();
        assert_eq!(found, expected);
    }

    /// Stub-bearing ions of a chain on one arm of an asymmetric linker
    fn arm_ions(p: &Peptide, mods: &Modifications, site: usize, arm: &ArmStubs) -> Vec<(Kind, f32)> {
        let stubs = CleavableStubs { short: arm.short, long: arm.long, alt: arm.alt.clone(), other: None };
        stubbed_ions(p, mods, site, &stubs)
    }

    #[test]
    fn clixlink_stubs_and_validation() {
        let cfg: CrosslinkSettings = serde_json::from_str(r#"{"reagent": "clixlink"}"#).unwrap();
        assert!(cfg.validate().is_ok());
        assert!((cfg.linker_mass() - 267.05653).abs() < 1e-3);
        let stubs = cfg.stubs().unwrap();
        let arms = stubs.arms();
        assert_eq!(arms.len(), 2);
        // a cleavage leaves one arm's alkene and the other arm's sulfenic acid
        assert!((arms[0].short + arms[1].long - cfg.linker_mass()).abs() < 1e-3);
        assert!((arms[1].short + arms[0].long - cfg.linker_mass()).abs() < 1e-3);
        // both arms: alkene / sulfenic acid doublets 49.98 apart, thiol = sulfenic - H2O
        for arm in &arms {
            assert!((arm.delta() - 49.98264).abs() < 1e-3);
            assert!((arm.long - arm.alt[0] - 18.010565).abs() < 1e-3);
        }

        // custom asymmetric linker: the cross sums must match the linker mass
        let bad: CrosslinkSettings = serde_json::from_str(
            r#"{"mass": 267.05653, "cleavable": {"short": 54.01056, "long": 103.9932,
                "other": {"short": 160.0, "long": 213.04596}}}"#,
        )
        .unwrap();
        assert!(bad.validate().unwrap_err().contains("other arm's long stub"));
        let good: CrosslinkSettings = serde_json::from_str(
            r#"{"mass": 267.05653, "cleavable": {"short": 54.01056, "long": 103.9932, "alt": [85.98264],
                "other": {"short": 163.06333, "long": 213.04596, "alt": [195.0354]}}}"#,
        )
        .unwrap();
        assert!(good.validate().is_ok());
        assert_eq!(good.stubs(), cfg.stubs());
    }

    /// cliXlink leaves different stubs on its two arms; either arm may sit on
    /// either peptide, and both ways round must be found
    #[test]
    fn finds_asymmetric_cleavable_crosslink_both_ways() {
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let (site_a, site_b) = (3usize, 9usize);
        let cfg = CrosslinkSettings { reagent: Some(Reagent::Clixlink), ..Default::default() };
        let arms = cfg.stubs().unwrap().arms();

        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let xl = CrosslinkScorer::new(&sc, &cfg);
        let precursor = a.monoisotopic + b.monoisotopic + cfg.linker_mass();

        for (arm_a, arm_b) in [(0usize, 1usize), (1, 0)] {
            let (x, y) = (&arms[arm_a], &arms[arm_b]);
            let mut ions = arm_ions(&a, &db.modifications, site_a, x);
            ions.extend(arm_ions(&b, &db.modifications, site_b, y));
            // each peptide's doublet at its own arm's stubs
            ions.push((Kind::B, a.monoisotopic + x.short));
            ions.push((Kind::B, a.monoisotopic + x.long));
            ions.push((Kind::B, b.monoisotopic + y.short));
            ions.push((Kind::B, b.monoisotopic + y.long));
            let hits = xl.score(&spectrum(precursor, 4, &ions));
            assert_eq!(hits.len(), 1, "arms {arm_a}/{arm_b}: {hits:?}");
            let hit = &hits[0];
            assert_eq!(hit.stub_peaks, 4);
            assert!(hit.precursor_ppm.abs() < 1.0, "ppm {}", hit.precursor_ppm);
            assert_eq!(
                pair_set(&db, hit),
                vec![("GLEKVSTAAMQR".to_string(), site_a), ("TDYHFPNEWKLIR".to_string(), site_b)]
            );
        }

        // DSSO-like doublets on both peptides are not cliXlink: the tagged arm
        // (+163.06 / +213.05) is missing
        let mut ions = arm_ions(&a, &db.modifications, site_a, &arms[0]);
        ions.extend(arm_ions(&b, &db.modifications, site_b, &arms[0]));
        for p in [&a, &b] {
            ions.push((Kind::B, p.monoisotopic + arms[0].short));
            ions.push((Kind::B, p.monoisotopic + arms[0].long));
        }
        assert!(xl.score(&spectrum(precursor, 4, &ions)).is_empty());
    }

    /// Ions of a DSSO loop-link opened in the collision cell: the first site
    /// keeps the short stub, the second the long one
    fn opened_loop(p: &Peptide, mods: &Modifications, first: usize, second: usize, cfg: &CrosslinkSettings) -> Vec<(Kind, f32)> {
        let s = cfg.stubs().unwrap();
        let short = ArmStubs { short: s.short, long: s.short, alt: vec![] };
        let long = ArmStubs { short: s.long, long: s.long, alt: vec![] };
        let (mut all, _) = cleaved_loop_ions(p, mods, &[Kind::B, Kind::Y], first, second, &short, &long, cfg.linker_mass());
        all.dedup();
        all
    }

    #[test]
    fn cleaved_loop_ions_follow_the_ring() {
        let l = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_L", Position::Internal);
        let db = database(vec![l.clone()]);
        let cfg = dsso();
        let arm = &cfg.stubs().unwrap().arms()[0];
        let (all, ring) = cleaved_loop_ions(&l, &db.modifications, &[Kind::B, Kind::Y], 9, 18, arm, arm, cfg.linker_mass());
        // b10..b18 and y9..y17 lie inside the ring (9 + 9 positions), three stubs each
        assert_eq!(ring.len(), 18 * 3);
        // spanning ions: the linker (alkene + sulfenic acid) or alkene + thiol
        let b19 = IonSeries::new(&l, Kind::B, &db.modifications).nth(18).unwrap().monoisotopic_mass;
        // b19 with each stub combination (b20 = b19 + G is 57 Da further out)
        let spanning: Vec<f32> = all
            .iter()
            .filter(|(k, m)| *k == Kind::B && *m - b19 > 100.0 && *m - b19 < 160.0)
            .map(|(_, m)| m - b19)
            .collect();
        assert_eq!(spanning.len(), 2, "{spanning:?}");
        assert!(spanning.iter().any(|d| (d - cfg.linker_mass()).abs() < 1e-3));
        assert!(spanning.iter().any(|d| (d - (54.01056 + 85.98264)).abs() < 1e-3));
    }

    #[test]
    fn finds_cleavable_loop_link() {
        let l = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_L", Position::Internal);
        let mut targets = distractors(400);
        targets.push(l.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = dsso();
        let xl = CrosslinkScorer::new(&sc, &cfg);

        let ions = opened_loop(&l, &db.modifications, 9, 18, &cfg);
        let hits = xl.score(&spectrum(l.monoisotopic + cfg.linker_mass(), 3, &ions));
        assert_eq!(hits.len(), 1, "{hits:?}");
        let hit = &hits[0];
        assert_eq!(hit.kind, LinkKind::Loop);
        assert_eq!(hit.class, PairClass::TT);
        assert_eq!((hit.alpha_site, hit.beta_site), (9, 18));
        assert_eq!(String::from_utf8(db[hit.alpha].sequence.to_vec()).unwrap(), "AEFGHDSTWKLMNPQSTVKGYEDFAR");
        assert!(hit.alpha_linked >= 9, "ring ions {}", hit.alpha_linked);
        assert!(hit.precursor_ppm.abs() < 1.0);

        // switched off: no loop-links
        let off = CrosslinkSettings { loop_links: false, ..dsso() };
        assert!(CrosslinkScorer::new(&sc, &off)
            .score(&spectrum(l.monoisotopic + cfg.linker_mass(), 3, &ions))
            .is_empty());
    }

    #[test]
    fn cleavable_loop_link_needs_ring_ions() {
        // the ring never opened (or another +158 Da modification): only ions
        // outside the ring and spanning it, none that only an opened loop gives
        let l = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_L", Position::Internal);
        let mut targets = distractors(400);
        targets.push(l.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = dsso();
        let closed = chain_ions(
            &l,
            &db.modifications,
            &[Kind::B, Kind::Y],
            Link::Loop { first: 9, second: 18, shift: cfg.linker_mass() },
        );
        let query = spectrum(l.monoisotopic + cfg.linker_mass(), 3, &closed);
        assert!(CrosslinkScorer::new(&sc, &cfg).score(&query).is_empty());

        // nor is a linear peptide spectrum at the loop-link's precursor mass
        let linear = chain_ions(&l, &db.modifications, &[Kind::B, Kind::Y], Link::Single { site: 0, shift: 0.0 });
        let query = spectrum(l.monoisotopic + cfg.linker_mass(), 3, &linear);
        assert!(CrosslinkScorer::new(&sc, &cfg).score(&query).is_empty());
    }

    /// Ions of one peptide of a multivalent species: `groups` are the
    /// residues each group of other peptides hangs on, and the group's mass.
    /// An ion holding only some of a group's residues does not come apart.
    fn species_ions(p: &Peptide, mods: &Modifications, groups: &[(&[usize], f32)]) -> Vec<(Kind, f32)> {
        let mut out = Vec::new();
        for kind in [Kind::B, Kind::Y] {
            for (idx, ion) in IonSeries::new(p, kind, mods).enumerate() {
                let contains = |s: usize| if kind == Kind::B { idx >= s } else { idx < s };
                let mut mass = ion.monoisotopic_mass;
                let mut apart = true;
                for (sites, shift) in groups {
                    let held = sites.iter().filter(|&&s| contains(s)).count();
                    if held == sites.len() {
                        mass += shift;
                    } else if held > 0 {
                        apart = false;
                    }
                }
                if apart {
                    out.push((kind, mass));
                }
            }
        }
        out
    }

    fn multivalent(n: u8) -> CrosslinkSettings {
        CrosslinkSettings { multivalent: Some(n), ..Default::default() }
    }

    fn sites_of(db: &IndexedDatabase, s: &Species) -> Vec<(String, Vec<u16>)> {
        let mut out: Vec<(String, Vec<u16>)> = (0..s.peptides.len())
            .map(|x| {
                let mut sites: Vec<u16> = s
                    .links
                    .iter()
                    .flat_map(|&(a, sa, b, sb)| [(a, sa), (b, sb)])
                    .filter(|&(p, _)| p as usize == x)
                    .map(|(_, site)| site)
                    .collect();
                sites.sort_unstable();
                (String::from_utf8(db[s.peptides[x]].sequence.to_vec()).unwrap(), sites)
            })
            .collect();
        out.sort();
        out
    }

    #[test]
    fn finds_doubly_linked_pair() {
        // a2b: two DSS links between the same two peptides (K10, K19 on a; K9, K14 on b)
        let a = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEKWLIGKSTPR", "PROT_B", Position::Internal);
        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = multivalent(2);
        let mods = &db.modifications;
        let mut ions = species_ions(&a, mods, &[(&[9, 18], b.monoisotopic + 2.0 * DSS)]);
        ions.extend(species_ions(&b, mods, &[(&[8, 13], a.monoisotopic + 2.0 * DSS)]));
        let query = spectrum(a.monoisotopic + b.monoisotopic + 2.0 * DSS, 4, &ions);

        let hits = CrosslinkScorer::new(&sc, &cfg).score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        let hit = &hits[0];
        assert_eq!(hit.kind, LinkKind::Multi);
        assert_eq!(hit.class, PairClass::TT);
        assert!(hit.precursor_ppm.abs() < 1.0);
        let s = hit.species.as_ref().unwrap();
        assert_eq!(s.topology(), "a2b");
        assert_eq!(
            sites_of(&db, s),
            vec![("AEFGHDSTWKLMNPQSTVKGYEDFAR".into(), vec![9, 18]), ("TDYHFPNEKWLIGKSTPR".into(), vec![8, 13])]
        );
        assert!(s.linked.iter().all(|&l| l >= cfg.min_linked_peaks));

        // off by default: a two-linker species is not looked for
        assert!(CrosslinkScorer::new(&sc, &CrosslinkSettings::default()).score(&query).is_empty());
    }

    #[test]
    fn finds_three_peptide_chain() {
        // a1b1c: a (K4) – b (K10) and b (K19) – c (K10)
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_B", Position::Internal);
        let c = peptide("TDYHFPNEWKLIR", "PROT_C", Position::Internal);
        let mut targets = distractors(400);
        targets.extend([a.clone(), b.clone(), c.clone()]);
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = multivalent(2);
        let mods = &db.modifications;
        let (ma, mb, mc) = (a.monoisotopic, b.monoisotopic, c.monoisotopic);
        let mut ions = species_ions(&a, mods, &[(&[3], mb + mc + 2.0 * DSS)]);
        ions.extend(species_ions(&b, mods, &[(&[9], ma + DSS), (&[18], mc + DSS)]));
        ions.extend(species_ions(&c, mods, &[(&[9], ma + mb + 2.0 * DSS)]));
        let query = spectrum(ma + mb + mc + 2.0 * DSS, 4, &ions);

        let hits = CrosslinkScorer::new(&sc, &cfg).score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        let s = hits[0].species.as_ref().unwrap();
        assert_eq!(hits[0].kind, LinkKind::Multi);
        assert_eq!(s.peptides.len(), 3);
        assert_eq!(s.links.len(), 2);
        // b in the middle: a chain of three distinct peptides
        assert_eq!(s.topology(), "a1b1c");
        assert_eq!(
            sites_of(&db, s),
            vec![
                ("AEFGHDSTWKLMNPQSTVKGYEDFAR".into(), vec![9, 18]),
                ("GLEKVSTAAMQR".into(), vec![3]),
                ("TDYHFPNEWKLIR".into(), vec![9]),
            ]
        );
    }

    #[test]
    fn multivalent_search_keeps_plain_pairs() {
        // an ordinary crosslink stays a pair with the multivalent search on
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let cfg = multivalent(3);
        let query = pair_spectrum(&db, &a, &b, 3, 9);
        let hits = CrosslinkScorer::new(&sc, &cfg).score(&query);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].kind, LinkKind::Inter);
    }

    #[test]
    fn species_topology_names() {
        let s = |peptides: &[u32], links: &[(u8, u8)]| Species {
            peptides: peptides.iter().map(|&p| PeptideIx(p)).collect(),
            links: links.iter().map(|&(a, b)| (a, 0, b, 0)).collect(),
            ..Default::default()
        };
        assert_eq!(s(&[1, 2], &[(0, 1), (0, 1)]).topology(), "a2b");
        assert_eq!(s(&[1, 1], &[(0, 1), (0, 1)]).topology(), "a2a");
        assert_eq!(s(&[1, 1, 2], &[(0, 1), (1, 2)]).topology(), "a1a1b");
        assert_eq!(s(&[1, 2, 2], &[(0, 1), (1, 2)]).topology(), "a1b1b");
        assert_eq!(s(&[5, 6, 7], &[(0, 1), (1, 2), (1, 2)]).topology(), "a1b2c");
        // the chain is walked from an end and lettered in that order, whatever
        // the order of the peptides
        assert_eq!(s(&[5, 6, 7], &[(0, 2), (1, 2)]).topology(), "a1b1c");
        assert_eq!(s(&[1, 2, 3, 4], &[(0, 1), (1, 2), (1, 3)]).topology(), "a1b,b1c,b1d");
    }

    #[test]
    fn multivalent_fdr_and_settings() {
        // two distinct peptides: each decoy species is 1/3 false target, as
        // (TD + DD) / 3 for pairs
        let (t, d) = ((false, 2), (true, 2));
        let q = q_multivalent(&[t, t, t, d, t, d, t, d]);
        assert_eq!(q[0], 0.0);
        assert!((q[7] - 1.0 / 5.0).abs() < 1e-6, "{q:?}");
        // three distinct peptides: 1/7 each; copies of one peptide (a2a): 1 each
        assert!((q_multivalent(&[(false, 3), (true, 3)])[1] - 1.0 / 7.0).abs() < 1e-6);
        assert_eq!(q_multivalent(&[(false, 1), (true, 1)])[1], 1.0);
        // many all-decoy species never push the estimate below zero
        let q = q_multivalent(&[t, (true, 2), (true, 2), (true, 2), t]);
        assert!(q.iter().all(|&v| v > 0.0 || v == 0.0) && q[4] > 0.4, "{q:?}");
        assert!(multivalent(2).validate().is_ok());
        assert!(multivalent(4).validate().is_ok());
        assert!(multivalent(5).validate().unwrap_err().contains("2 to 4"));
        // MS-cleavable linkers too
        assert!(CrosslinkSettings { multivalent: Some(3), ..dsso() }.validate().is_ok());
        let cfg: CrosslinkSettings = serde_json::from_str(r#"{"multivalent": 3}"#).unwrap();
        assert_eq!(cfg.multivalent, Some(3));
        let off: CrosslinkSettings = serde_json::from_str("{}").unwrap();
        assert_eq!(off.multivalent, None);
    }

    /// Ions of one peptide of a DSSO species whose linkers all broke: each
    /// link residue keeps the stub given for it, and an ion carries the stubs
    /// of the residues it holds
    fn cleaved_species_ions(p: &Peptide, mods: &Modifications, stubs: &[(usize, f32)]) -> Vec<(Kind, f32)> {
        let mut out = Vec::new();
        for kind in [Kind::B, Kind::Y] {
            for (idx, ion) in IonSeries::new(p, kind, mods).enumerate() {
                let contains = |s: usize| if kind == Kind::B { idx >= s } else { idx < s };
                let shift: f32 = stubs.iter().filter(|(s, _)| contains(*s)).map(|(_, m)| m).sum();
                out.push((kind, ion.monoisotopic_mass + shift));
            }
        }
        out
    }

    #[test]
    fn finds_cleavable_doubly_linked_pair() {
        // a2b with DSSO: each link breaks into alkene (A) on one side, sulfenic
        // acid (S) on the other
        let a = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEKWLIGKSTPR", "PROT_B", Position::Internal);
        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings { multivalent: Some(2), ..dsso() };
        let s = cfg.stubs().unwrap();
        let mods = &db.modifications;
        let mut ions = cleaved_species_ions(&a, mods, &[(9, s.short), (18, s.long)]);
        ions.extend(cleaved_species_ions(&b, mods, &[(8, s.long), (13, s.short)]));
        let fragments = ions.clone();
        // the free peptides with both their stubs: a triplet member pair for a, one for b
        ions.push((Kind::B, a.monoisotopic + s.short + s.long));
        ions.push((Kind::B, a.monoisotopic + 2.0 * s.long));
        ions.push((Kind::B, b.monoisotopic + s.short + s.long));
        let precursor = a.monoisotopic + b.monoisotopic + 2.0 * cfg.linker_mass();

        let hits = CrosslinkScorer::new(&sc, &cfg).score(&spectrum(precursor, 4, &ions));
        assert_eq!(hits.len(), 1, "{hits:?}");
        let hit = &hits[0];
        assert_eq!(hit.kind, LinkKind::Multi);
        assert_eq!(hit.class, PairClass::TT);
        assert!(hit.stub_peaks >= 3, "stub peaks {}", hit.stub_peaks);
        let sp = hit.species.as_ref().unwrap();
        assert_eq!(sp.topology(), "a2b");
        assert_eq!(
            sites_of(&db, sp),
            vec![("AEFGHDSTWKLMNPQSTVKGYEDFAR".into(), vec![9, 18]), ("TDYHFPNEKWLIGKSTPR".into(), vec![8, 13])]
        );

        // the same fragments without any peptide + stubs peak: not reported
        assert!(CrosslinkScorer::new(&sc, &cfg).score(&spectrum(precursor, 4, &fragments)).is_empty());
    }

    #[test]
    fn finds_cleavable_three_peptide_chain() {
        // a1b1c with DSSO: a (K4) – b (K10), b (K19) – c (K10)
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_B", Position::Internal);
        let c = peptide("TDYHFPNEWKLIR", "PROT_C", Position::Internal);
        let mut targets = distractors(400);
        targets.extend([a.clone(), b.clone(), c.clone()]);
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings { multivalent: Some(2), ..dsso() };
        let s = cfg.stubs().unwrap();
        let mods = &db.modifications;
        let (ma, mb, mc, l) = (a.monoisotopic, b.monoisotopic, c.monoisotopic, cfg.linker_mass());
        let mut ions = cleaved_species_ions(&a, mods, &[(3, s.short)]);
        ions.extend(cleaved_species_ions(&b, mods, &[(9, s.long), (18, s.short)]));
        ions.extend(cleaved_species_ions(&c, mods, &[(9, s.long)]));
        // free peptides with their stubs (a doublet for a)
        ions.push((Kind::B, ma + s.short));
        ions.push((Kind::B, ma + s.long));
        ions.push((Kind::B, mb + s.long + s.short));
        ions.push((Kind::B, mc + s.long));
        // only one linker broken: b1c with b's stub, a1b with b's stub
        ions.push((Kind::B, mb + mc + l + s.long));
        ions.push((Kind::B, ma + mb + l + s.short));

        let hits = CrosslinkScorer::new(&sc, &cfg).score(&spectrum(ma + mb + mc + 2.0 * l, 4, &ions));
        assert_eq!(hits.len(), 1, "{hits:?}");
        let sp = hits[0].species.as_ref().unwrap();
        assert_eq!(sp.topology(), "a1b1c");
        assert!(sp.topology_peaks >= 2, "topology peaks {}", sp.topology_peaks);
        assert_eq!(
            sites_of(&db, sp),
            vec![
                ("AEFGHDSTWKLMNPQSTVKGYEDFAR".into(), vec![9, 18]),
                ("GLEKVSTAAMQR".into(), vec![3]),
                ("TDYHFPNEWKLIR".into(), vec![9]),
            ]
        );
        // b sits in the middle: it is linked to both a and c
        let middle = sp.peptides.iter().position(|&ix| db[ix].sequence.len() == 26).unwrap() as u8;
        assert!(sp.links.iter().all(|&(x, _, y, _)| x == middle || y == middle));
    }

    #[test]
    fn cleavable_search_needs_the_doublets() {
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let cfg = dsso();
        let stubs = cfg.stubs().unwrap();
        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let xl = CrosslinkScorer::new(&sc, &cfg);

        // good fragment ions, but no stub-peptide doublets
        let mut ions = stubbed_ions(&a, &db.modifications, 3, &stubs);
        ions.extend(stubbed_ions(&b, &db.modifications, 9, &stubs));
        let query = spectrum(a.monoisotopic + b.monoisotopic + cfg.linker_mass(), 4, &ions);
        assert!(xl.score(&query).is_empty());

        // one peptide's doublet only (no complement) is also rejected
        let mut partial = ions.clone();
        partial.push((Kind::B, a.monoisotopic + stubs.short));
        partial.push((Kind::B, a.monoisotopic + stubs.long));
        let query = spectrum(a.monoisotopic + b.monoisotopic + cfg.linker_mass(), 4, &partial);
        assert!(xl.score(&query).is_empty());
    }

    fn kernel_for(db: &IndexedDatabase, cfg: &CrosslinkSettings) -> crate::crosslink_kernel::CpuKernel {
        crate::crosslink_kernel::CpuKernel::new(crate::crosslink_kernel::IonTable::build(db, |p| cfg.link_sites(p)))
    }

    fn pair_set(db: &IndexedDatabase, hit: &CrosslinkMatch) -> Vec<(String, usize)> {
        let seq = |ix: PeptideIx| String::from_utf8(db[ix].sequence.to_vec()).unwrap();
        let mut v = vec![(seq(hit.alpha), hit.alpha_site as usize)];
        if let Some(b) = hit.beta {
            v.push((seq(b), hit.beta_site as usize));
        }
        v.sort();
        v
    }

    /// Fixture: a pair-link spectrum over a 400-peptide distractor database
    fn pair_fixture() -> (IndexedDatabase, Peptide, Peptide) {
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let mut targets = distractors(400);
        targets.push(a.clone());
        targets.push(b.clone());
        (database(targets), a, b)
    }

    fn pair_spectrum(db: &IndexedDatabase, a: &Peptide, b: &Peptide, sa: usize, sb: usize) -> ProcessedSpectrum {
        let kinds = [Kind::B, Kind::Y];
        let mut ions = chain_ions(a, &db.modifications, &kinds, Link::Single { site: sa, shift: b.monoisotopic + DSS });
        ions.extend(chain_ions(b, &db.modifications, &kinds, Link::Single { site: sb, shift: a.monoisotopic + DSS }));
        spectrum(a.monoisotopic + b.monoisotopic + DSS, 3, &ions)
    }

    #[test]
    fn reflect_finds_inter_peptide_link() {
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let kernel = kernel_for(&db, &cfg);
        let xl = CrosslinkScorer::new(&sc, &cfg).with_kernel(Some(&kernel));
        let hits = xl.score(&pair_spectrum(&db, &a, &b, 3, 9));
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].kind, LinkKind::Inter);
        assert_eq!(hits[0].class, PairClass::TT);
        assert!(hits[0].precursor_ppm.abs() < 1.0);
        assert_eq!(
            pair_set(&db, &hits[0]),
            vec![("GLEKVSTAAMQR".to_string(), 3), ("TDYHFPNEWKLIR".to_string(), 9)]
        );
    }

    #[test]
    fn reflect_agrees_with_cpu_backend() {
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let kernel = kernel_for(&db, &cfg);
        let cpu = CrosslinkScorer::new(&sc, &cfg);
        let reflect = CrosslinkScorer::new(&sc, &cfg).with_kernel(Some(&kernel));

        // a few different true links, including a homodimer-style self pair
        let c = peptide("AKCDKEFRQLMNR", "PROT_C", Position::Internal);
        let mut targets = distractors(400);
        for p in [&a, &b, &c] {
            targets.push(p.clone());
        }
        let db2 = database(targets);
        let sc2 = scorer(&db2);
        let kernel2 = kernel_for(&db2, &cfg);
        let cpu2 = CrosslinkScorer::new(&sc2, &cfg);
        let reflect2 = CrosslinkScorer::new(&sc2, &cfg).with_kernel(Some(&kernel2));
        let _ = (&cpu, &reflect);

        for (x, y, sx, sy) in [(&a, &b, 3, 9), (&b, &c, 9, 1), (&a, &c, 3, 4), (&c, &c, 1, 4)] {
            let q = pair_spectrum(&db2, x, y, sx, sy);
            let (h_cpu, h_ref) = (cpu2.score(&q), reflect2.score(&q));
            assert_eq!(h_cpu.len(), 1, "cpu found nothing for {sx},{sy}");
            assert_eq!(h_ref.len(), 1, "reflect found nothing for {sx},{sy}");
            assert_eq!(pair_set(&db2, &h_cpu[0]), pair_set(&db2, &h_ref[0]));
            assert!((h_cpu[0].hyperscore - h_ref[0].hyperscore).abs() < 1e-6);
            assert_eq!(h_cpu[0].class, h_ref[0].class);
        }
    }

    #[test]
    fn reflect_finds_loop_link_and_rejects_linear() {
        let l = peptide("AEFGHDSTWKLMNPQSTVKGYEDFAR", "PROT_L", Position::Internal);
        let mut targets = distractors(400);
        targets.push(l.clone());
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let kernel = kernel_for(&db, &cfg);
        let xl = CrosslinkScorer::new(&sc, &cfg).with_kernel(Some(&kernel));

        let ions = chain_ions(&l, &db.modifications, &[Kind::B, Kind::Y], Link::Loop { first: 9, second: 18, shift: DSS });
        let hits = xl.score(&spectrum(l.monoisotopic + DSS, 3, &ions));
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].kind, LinkKind::Loop);
        assert_eq!((hits[0].alpha_site, hits[0].beta_site), (9, 18));

        let linear = IonSeries::new(&l, Kind::B, &db.modifications)
            .chain(IonSeries::new(&l, Kind::Y, &db.modifications))
            .map(|i| (i.kind, i.monoisotopic_mass))
            .collect::<Vec<_>>();
        assert!(xl.score(&spectrum(l.monoisotopic, 2, &linear)).is_empty());
    }

    #[test]
    fn chain_filter_keeps_the_true_chains() {
        use crate::crosslink_kernel::{ChainKernel, KernelParams};
        let (db, a, b) = pair_fixture();
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let kernel = kernel_for(&db, &cfg);
        let q = pair_spectrum(&db, &a, &b, 3, 9);
        let m = a.monoisotopic + b.monoisotopic + DSS;
        let params = KernelParams {
            precursor_mass: m,
            max_fragment_charge: 3,
            fragment_tol: sc.fragment_tol.into(),
            slack: 0.05,
            min_matched: 3,
            max_survivors: 40,
        };
        let survivors = kernel.survivors(&q, &params);
        assert!(survivors.windows(2).all(|w| w[0].peptide < w[1].peptide));
        let find = |p: &Peptide| {
            survivors
                .iter()
                .find(|s| db[PeptideIx(s.peptide)].sequence == p.sequence && !db[PeptideIx(s.peptide)].decoy)
                .unwrap_or_else(|| panic!("{:?} missing", String::from_utf8_lossy(&p.sequence)))
        };
        assert_eq!(find(&a).site, 3);
        assert_eq!(find(&b).site, 9);
        // the top-K limit is honoured, and the true chains rank inside it
        assert!(survivors.len() <= 40, "{} survivors", survivors.len());
    }

    #[test]
    fn linear_tolerance_matches_sage_bounds() {
        use crate::crosslink_kernel::LinearTol;
        for tol in [Tolerance::Ppm(-10.0, 12.0), Tolerance::Da(-0.02, 0.03), Tolerance::Pct(-0.001, 0.002)] {
            let lin: LinearTol = tol.into();
            for c in [150.0f32, 987.654, 2500.0] {
                let (a, b) = (tol.bounds(c), lin.bounds(c));
                assert!((a.0 - b.0).abs() < 1e-3 && (a.1 - b.1).abs() < 1e-3, "{tol:?} {c}");
            }
        }
    }

    #[test]
    fn backend_selection_rules() {
        let (db, _, _) = pair_fixture();
        let cfg = |json: &str| serde_json::from_str::<CrosslinkSettings>(json).unwrap();
        assert!(make_kernel(&cfg("{}"), &db).unwrap().is_none());
        assert!(make_kernel(&cfg(r#"{"backend":"reflect"}"#), &db).unwrap().is_some());
        // cleavable reagents ignore the backend
        assert!(make_kernel(&cfg(r#"{"reagent":"dsso","backend":"reflect"}"#), &db).unwrap().is_none());
        // `auto` never fails; explicit `gpu` fails when no GPU backend is available
        assert!(make_kernel(&cfg(r#"{"backend":"auto"}"#), &db).is_ok());
        #[cfg(not(feature = "gpu"))]
        assert!(make_kernel(&cfg(r#"{"backend":"gpu"}"#), &db).err().unwrap().contains("--features gpu"));
        assert!(serde_json::from_str::<CrosslinkSettings>(r#"{"backend":"cuda"}"#).is_err());
    }

    /// The GPU kernel on every adapter of this machine (hardware and software,
    /// each driver API), or none when there is no GPU support
    #[cfg(feature = "gpu")]
    fn gpu_kernels(db: &IndexedDatabase, cfg: &CrosslinkSettings) -> Vec<crate::crosslink_gpu::GpuKernel> {
        use crate::crosslink_kernel::{ChainKernel, IonTable};
        crate::crosslink_gpu::GpuKernel::on_every_adapter(|| IonTable::build(db, |p| cfg.link_sites(p)))
            .into_iter()
            .filter_map(|k| match k {
                Ok(k) => {
                    eprintln!("GPU kernel on: {}", k.name());
                    Some(k)
                }
                Err(e) => {
                    eprintln!("adapter skipped: {e}");
                    None
                }
            })
            .collect()
    }

    /// Survivors of `gpu` must equal the CPU reference exactly (score to float accuracy)
    #[cfg(feature = "gpu")]
    fn assert_same_survivors(
        want: &[crate::crosslink_kernel::Survivor],
        got: &[crate::crosslink_kernel::Survivor],
        what: &str,
    ) {
        assert_eq!(
            want.iter().map(|s| s.peptide).collect::<Vec<_>>(),
            got.iter().map(|s| s.peptide).collect::<Vec<_>>(),
            "survivor sets differ ({what})"
        );
        for (w, g) in want.iter().zip(got) {
            assert_eq!((w.site, w.matched), (g.site, g.matched), "peptide {} ({what})", w.peptide);
            assert!((w.score - g.score).abs() < 1e-3, "peptide {}: {} vs {} ({what})", w.peptide, w.score, g.score);
        }
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn gpu_kernel_matches_cpu_reference() {
        use crate::crosslink_kernel::{ChainKernel, KernelParams};
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let c = peptide("AKCDKEFRQLMNR", "PROT_C", Position::Internal);
        let mut targets = distractors(400);
        for p in [&a, &b, &c] {
            targets.push(p.clone());
        }
        let db = database(targets);
        let sc = scorer(&db);
        let cfg = CrosslinkSettings::default();
        let cpu = kernel_for(&db, &cfg);
        let gpus = gpu_kernels(&db, &cfg);
        if gpus.is_empty() {
            eprintln!("skipping GPU parity test: no adapter");
            return;
        }

        for gpu in &gpus {
            let mut compared = 0usize;
            for (x, y, sx, sy) in [(&a, &b, 3, 9), (&b, &c, 9, 1), (&a, &c, 3, 4), (&c, &c, 1, 4)] {
                let q = pair_spectrum(&db, x, y, sx, sy);
                for (max_charge, slack, min_matched) in [(2u8, 0.0f32, 1u16), (3, 0.05, 3), (4, 0.03, 2)] {
                    let params = KernelParams {
                        precursor_mass: x.monoisotopic + y.monoisotopic + DSS,
                        max_fragment_charge: max_charge,
                        fragment_tol: sc.fragment_tol.into(),
                        slack,
                        min_matched,
                        max_survivors: usize::MAX,
                    };
                    let want = cpu.survivors(&q, &params);
                    let what = format!("{}, charge<{max_charge}, slack {slack}, min {min_matched}", gpu.name());
                    assert_same_survivors(&want, &gpu.survivors(&q, &params), &what);
                    compared += want.len();
                }
            }
            assert!(compared > 100, "parity test compared too little: {compared}");
            eprintln!("{}: GPU and CPU agree on {compared} survivors", gpu.name());
        }
    }

    /// Peaks exactly on a tolerance window's edges count on every device:
    /// rounding differences (GPU division, fused multiply-add) once dropped them
    #[cfg(feature = "gpu")]
    #[test]
    fn gpu_kernel_agrees_on_window_edges() {
        use crate::crosslink_kernel::{ChainKernel, IonTable, KernelParams};
        let a = peptide("GLEKVSTAAMQR", "PROT_A", Position::Internal);
        let b = peptide("TDYHFPNEWKLIR", "PROT_B", Position::Internal);
        let mut targets = distractors(200);
        targets.push(a.clone());
        targets.push(b.clone());
        let db = database(targets);
        let cfg = CrosslinkSettings::default();
        let cpu = kernel_for(&db, &cfg);
        let gpus = gpu_kernels(&db, &cfg);
        if gpus.is_empty() {
            eprintln!("skipping GPU window-edge test: no adapter");
            return;
        }
        let table = IonTable::build(&db, |p| cfg.link_sites(p));
        let ia = db.peptides.iter().position(|p| p.sequence == a.sequence).unwrap();

        for tol in [Tolerance::Ppm(-10.0, 10.0), Tolerance::Ppm(-7.0, 13.0), Tolerance::Da(-0.02, 0.02)] {
            for slack in [0.0f32, 0.03] {
                let params = KernelParams {
                    precursor_mass: a.monoisotopic + b.monoisotopic + DSS,
                    max_fragment_charge: 4,
                    fragment_tol: tol.into(),
                    slack,
                    min_matched: 1,
                    max_survivors: usize::MAX,
                };
                // a peak on each edge of every window of peptide `a`
                let shift = params.precursor_mass - table.peptide_mass[ia];
                let mut masses = Vec::new();
                for &mass in &table.ion_mass[table.ion_off[ia] as usize..table.ion_off[ia + 1] as usize] {
                    for (inv_z, slack_z) in params.charge_terms() {
                        let (lo, hi) = params.fragment_tol.bounds(mass * inv_z);
                        let (slo, shi) = params.fragment_tol.bounds((mass + shift) * inv_z);
                        masses.extend([lo, hi, slo - slack_z, shi + slack_z]);
                    }
                }
                masses.sort_by(f32::total_cmp);
                masses.dedup();
                let q = ProcessedSpectrum {
                    level: 2,
                    id: "edges".into(),
                    intensities: (0..masses.len()).map(|i| 50.0 + (i % 11) as f32).collect(),
                    charges: vec![1; masses.len()],
                    masses,
                    ..Default::default()
                };

                let want = cpu.survivors(&q, &params);
                let on_a = want.iter().find(|s| s.peptide as usize == ia).expect("peptide a survives");
                assert!(on_a.matched >= 20, "edge peaks should count: {on_a:?}");
                for gpu in &gpus {
                    assert_same_survivors(&want, &gpu.survivors(&q, &params), &format!("{}, {tol:?}, slack {slack}", gpu.name()));
                }
            }
        }
    }

    fn csm(score: f64, class: PairClass, kind: LinkKind, alpha: u32, beta: Option<u32>) -> CrosslinkMatch {
        CrosslinkMatch {
            hyperscore: score,
            class,
            kind,
            intra_protein: true,
            alpha: PeptideIx(alpha),
            beta: beta.map(PeptideIx),
            ..Default::default()
        }
    }

    #[test]
    fn q_values_use_td_minus_dd() {
        use PairClass::*;
        // sorted best-first: 8 TT, then TD, DD, TT ...
        let classes = [TT, TT, TT, TT, TT, TT, TT, TT, TD, DD, TT, TD, TD];
        let q = q_from_sorted(&classes, false);
        // After 8 TT: fdr 0. After TD: 1/8. After DD: (1-1)/8 = 0. After TT: 0.
        assert_eq!(q[7], 0.0);
        assert_eq!(q[9], 0.0);
        // After the next two TD: td=3, dd=1, tt=9 -> 2/9
        assert!((q[12] - 2.0 / 9.0).abs() < 1e-6);
        // monotone non-decreasing along the list
        assert!(q.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn q_values_loops_use_decoys_over_targets() {
        use PairClass::*;
        let q = q_from_sorted(&[TT, TT, TT, TT, DD, TT], true);
        // raw FDR is 1/4 after the decoy and 1/5 after the next target; a
        // q-value is the minimum FDR at or beyond a rank, so both become 0.2
        assert!((q[4] - 0.2).abs() < 1e-6);
        assert!((q[5] - 0.2).abs() < 1e-6);
        assert_eq!(q[3], 0.0);
    }

    #[test]
    fn assign_q_values_groups_and_collapses_links() {
        use PairClass::*;
        let mut m = vec![
            csm(50.0, TT, LinkKind::Inter, 1, Some(2)),
            csm(48.0, TT, LinkKind::Inter, 2, Some(1)), // same link, other chain order
            csm(40.0, TT, LinkKind::Inter, 3, Some(4)),
            csm(10.0, TD, LinkKind::Inter, 5, Some(6)),
        ];
        for x in &mut m {
            x.beta_site = 1;
            x.alpha_site = 1;
        }
        let (csms, links) = assign_q_values(&mut m);
        assert_eq!(csms, 3);
        assert_eq!(links, 2, "first two CSMs describe the same crosslink");
        assert_eq!(m[0].crosslink_q, m[1].crosslink_q);
        assert!(m[3].spectrum_q >= m[0].spectrum_q);
    }
}
