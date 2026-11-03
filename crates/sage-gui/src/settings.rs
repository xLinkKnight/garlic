//! The GUI's model of a Sage search, converted to and from the JSON
//! configuration file that `sage.exe` reads.
//!
//! Keys the GUI does not edit (e.g. `bruker_config`, LFQ peak scoring) are kept
//! when a configuration is loaded and written back unchanged.

use serde_json::{json, Map, Value};

/// Residues Sage accepts in modification targets and enzyme rules
pub const VALID_AA: &str = "ACDEFGHIKLMNPQRSTVWYUO";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TolUnit {
    Ppm,
    Da,
}

impl TolUnit {
    pub fn key(self) -> &'static str {
        match self {
            TolUnit::Ppm => "ppm",
            TolUnit::Da => "da",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Tolerance {
    pub unit: TolUnit,
    pub lo: f64,
    pub hi: f64,
}

impl Tolerance {
    fn ppm(lo: f64, hi: f64) -> Self {
        Self {
            unit: TolUnit::Ppm,
            lo,
            hi,
        }
    }

    fn to_json(&self) -> Value {
        let mut map = Map::new();
        map.insert(self.unit.key().into(), json!([self.lo, self.hi]));
        Value::Object(map)
    }

    fn from_json(v: &Value) -> Option<Self> {
        let (key, range) = v.as_object()?.iter().next()?;
        let unit = match key.as_str() {
            "ppm" => TolUnit::Ppm,
            "da" => TolUnit::Da,
            _ => return None,
        };
        let range = range.as_array()?;
        Some(Self {
            unit,
            lo: range.first()?.as_f64()?,
            hi: range.get(1)?.as_f64()?,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ModKind {
    Fixed,
    Variable,
}

#[derive(Clone, Debug)]
pub struct Modification {
    pub target: String,
    pub mass: f64,
    /// Optional ProForma name, e.g. `UNIMOD:35`
    pub name: String,
    pub kind: ModKind,
    /// Variable mods only: maximum occurrences per peptide, 0 = no limit
    pub max_count: u32,
}

/// A common modification offered in the "Add" menu
pub struct Preset {
    pub label: &'static str,
    pub target: &'static str,
    pub mass: f64,
    pub name: &'static str,
    pub kind: ModKind,
}

pub const PRESETS: &[Preset] = &[
    Preset { label: "Carbamidomethyl (C)", target: "C", mass: 57.021464, name: "UNIMOD:4", kind: ModKind::Fixed },
    Preset { label: "Propionamide (C)", target: "C", mass: 71.037114, name: "UNIMOD:24", kind: ModKind::Fixed },
    Preset { label: "Oxidation (M)", target: "M", mass: 15.994915, name: "UNIMOD:35", kind: ModKind::Variable },
    Preset { label: "Acetyl (protein N-term)", target: "[", mass: 42.010565, name: "UNIMOD:1", kind: ModKind::Variable },
    Preset { label: "Phospho (S)", target: "S", mass: 79.966331, name: "UNIMOD:21", kind: ModKind::Variable },
    Preset { label: "Phospho (T)", target: "T", mass: 79.966331, name: "UNIMOD:21", kind: ModKind::Variable },
    Preset { label: "Phospho (Y)", target: "Y", mass: 79.966331, name: "UNIMOD:21", kind: ModKind::Variable },
    Preset { label: "Deamidated (N)", target: "N", mass: 0.984016, name: "UNIMOD:7", kind: ModKind::Variable },
    Preset { label: "Deamidated (Q)", target: "Q", mass: 0.984016, name: "UNIMOD:7", kind: ModKind::Variable },
    Preset { label: "Gln->pyro-Glu (peptide N-term Q)", target: "^Q", mass: -17.026549, name: "UNIMOD:28", kind: ModKind::Variable },
    Preset { label: "Glu->pyro-Glu (peptide N-term E)", target: "^E", mass: -18.010565, name: "UNIMOD:27", kind: ModKind::Variable },
    Preset { label: "GlyGly / ubiquitin remnant (K)", target: "K", mass: 114.042927, name: "UNIMOD:121", kind: ModKind::Variable },
    Preset { label: "Methyl (K)", target: "K", mass: 14.01565, name: "UNIMOD:34", kind: ModKind::Variable },
    Preset { label: "Methyl (R)", target: "R", mass: 14.01565, name: "UNIMOD:34", kind: ModKind::Variable },
    Preset { label: "Dimethyl (K)", target: "K", mass: 28.0313, name: "UNIMOD:36", kind: ModKind::Variable },
    Preset { label: "Acetyl (K)", target: "K", mass: 42.010565, name: "UNIMOD:1", kind: ModKind::Variable },
    Preset { label: "TMT 6/10/11-plex (K)", target: "K", mass: 229.162932, name: "UNIMOD:737", kind: ModKind::Fixed },
    Preset { label: "TMT 6/10/11-plex (peptide N-term)", target: "^", mass: 229.162932, name: "UNIMOD:737", kind: ModKind::Fixed },
    Preset { label: "TMTpro 16/18-plex (K)", target: "K", mass: 304.207146, name: "UNIMOD:2016", kind: ModKind::Fixed },
    Preset { label: "TMTpro 16/18-plex (peptide N-term)", target: "^", mass: 304.207146, name: "UNIMOD:2016", kind: ModKind::Fixed },
    Preset { label: "DSS/BS3 mono-link, hydrolyzed (K)", target: "K", mass: 156.078644, name: "", kind: ModKind::Variable },
];

impl From<&Preset> for Modification {
    fn from(p: &Preset) -> Self {
        Self {
            target: p.target.into(),
            mass: p.mass,
            name: p.name.into(),
            kind: p.kind,
            max_count: 0,
        }
    }
}

/// Check a modification target such as `C`, `^`, `^Q`, `$`, `[`, `]M`
pub fn check_target(target: &str) -> Result<(), String> {
    let mut chars = target.chars();
    let first = chars.next().ok_or("target is empty")?;
    let rest: Vec<char> = chars.collect();
    let aa_ok = |c: char| VALID_AA.contains(c);
    match first {
        '^' | '$' | '[' | ']' => match rest[..] {
            [] => Ok(()),
            [c] if aa_ok(c) => Ok(()),
            _ => Err(format!("`{target}`: a terminus may be followed by at most one amino acid")),
        },
        c if rest.is_empty() && aa_ok(c) => Ok(()),
        _ => Err(format!(
            "`{target}`: use one amino acid ({VALID_AA}) or a terminus (^ $ [ ]) with an optional amino acid"
        )),
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Reagent {
    Dss,
    Bs3,
    Dsg,
    Bs2g,
    Dsso,
    Dsbu,
    Edc,
    Dmtmm,
    Clixlink,
    Custom,
}

impl Reagent {
    pub const ALL: [Reagent; 10] = [
        Reagent::Dss,
        Reagent::Bs3,
        Reagent::Dsg,
        Reagent::Bs2g,
        Reagent::Dsso,
        Reagent::Dsbu,
        Reagent::Edc,
        Reagent::Dmtmm,
        Reagent::Clixlink,
        Reagent::Custom,
    ];

    pub fn key(self) -> &'static str {
        match self {
            Reagent::Dss => "dss",
            Reagent::Bs3 => "bs3",
            Reagent::Dsg => "dsg",
            Reagent::Bs2g => "bs2g",
            Reagent::Dsso => "dsso",
            Reagent::Dsbu => "dsbu",
            Reagent::Edc => "edc",
            Reagent::Dmtmm => "dmtmm",
            Reagent::Clixlink => "clixlink",
            Reagent::Custom => "custom",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Reagent::Dss => "DSS (138.068 Da, non-cleavable)",
            Reagent::Bs3 => "BS3 (138.068 Da, non-cleavable)",
            Reagent::Dsg => "DSG (96.021 Da, non-cleavable)",
            Reagent::Bs2g => "BS2G (96.021 Da, non-cleavable)",
            Reagent::Dsso => "DSSO (158.004 Da, MS-cleavable)",
            Reagent::Dsbu => "DSBU (196.085 Da, MS-cleavable)",
            Reagent::Edc => "EDC (zero-length, -18.011 Da)",
            Reagent::Dmtmm => "DMTMM (zero-length, -18.011 Da)",
            Reagent::Clixlink => "cliXlink (267.057 Da, MS-cleavable, asymmetric)",
            Reagent::Custom => "Custom mass",
        }
    }

    pub fn is_cleavable(self) -> bool {
        self.stubs().is_some()
    }

    /// Linker mass of a preset (Da), as in `sage_core::crosslink::Reagent::mass`
    pub fn mass(self) -> Option<f64> {
        match self {
            Reagent::Dss | Reagent::Bs3 => Some(138.06808),
            Reagent::Dsg | Reagent::Bs2g => Some(96.02113),
            Reagent::Dsso => Some(158.00376),
            Reagent::Dsbu => Some(196.08479),
            Reagent::Edc | Reagent::Dmtmm => Some(-18.010565),
            Reagent::Clixlink => Some(267.05653),
            Reagent::Custom => None,
        }
    }

    /// Preset `(short, long, alt)` stub masses of an MS-cleavable reagent
    /// (for an asymmetric one, its first arm)
    pub fn stubs(self) -> Option<(f64, f64, &'static [f64])> {
        match self {
            Reagent::Dsso | Reagent::Clixlink => Some((54.01056, 103.9932, &[85.98264])),
            Reagent::Dsbu => Some((85.05276, 111.03203, &[])),
            _ => None,
        }
    }

    /// Asymmetric cleavable reagents: the stubs of the second arm
    pub fn second_arm(self) -> Option<(f64, f64, &'static [f64])> {
        match self {
            // the alkyne-tagged arm: alkene / sulfenic acid / thiol
            Reagent::Clixlink => Some((163.06333, 213.04596, &[195.0354])),
            _ => None,
        }
    }

    /// Default `(site1, site2)` for reagents with two different reactive ends
    pub fn sites(self) -> Option<(Site, Site)> {
        match self {
            Reagent::Edc | Reagent::Dmtmm => Some((
                Site { residues: "K".into(), protein_nterm: true, protein_cterm: false },
                Site { residues: "DE".into(), protein_nterm: false, protein_cterm: true },
            )),
            _ => None,
        }
    }
}

/// What one end of a heterobifunctional linker reacts with
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub residues: String,
    pub protein_nterm: bool,
    pub protein_cterm: bool,
}

impl Site {
    /// Sage's notation: residue letters, `[` protein N-terminus, `]` protein C-terminus
    pub fn to_spec(&self) -> String {
        let mut s = self.residues.clone();
        if self.protein_nterm {
            s.push('[');
        }
        if self.protein_cterm {
            s.push(']');
        }
        s
    }

    fn from_spec(spec: &str) -> Self {
        Site {
            residues: spec.chars().filter(char::is_ascii_alphabetic).collect(),
            protein_nterm: spec.contains('['),
            protein_cterm: spec.contains(']'),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Backend {
    Cpu,
    Reflect,
    Gpu,
    Auto,
}

impl Backend {
    pub const ALL: [Backend; 4] = [Backend::Cpu, Backend::Reflect, Backend::Gpu, Backend::Auto];

    pub fn key(self) -> &'static str {
        match self {
            Backend::Cpu => "cpu",
            Backend::Reflect => "reflect",
            Backend::Gpu => "gpu",
            Backend::Auto => "auto",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Crosslink {
    pub enabled: bool,
    pub reagent: Reagent,
    pub custom_mass: f64,
    pub custom_cleavable: bool,
    /// Fragment (stub) masses of an MS-cleavable linker, preset or custom
    pub stub_short: f64,
    pub stub_long: f64,
    /// Extra stub masses, comma separated
    pub stub_alt: String,
    /// Custom cleavable linkers: the two arms leave different stubs
    pub asymmetric: bool,
    /// Asymmetric linkers: stubs of the second arm
    pub stub2_short: f64,
    pub stub2_long: f64,
    pub stub2_alt: String,
    /// Homobifunctional: residues both ends react with
    pub residues: String,
    pub protein_nterm: bool,
    /// Heterobifunctional: each end reacts with a different group
    pub two_sites: bool,
    pub site1: Site,
    pub site2: Site,
    pub allow_cterm_residue: bool,
    pub loop_links: bool,
    pub candidates: u32,
    pub min_chain_peaks: u32,
    pub require_complete_doublets: bool,
    pub backend: Backend,
    pub chain_survivors: u32,
    /// Also run a linear peptide search (results.sage.tsv)
    pub linear_search: bool,
    /// Non-cleavable: matched linker-carrying ions required per chain
    pub min_linked_peaks: u32,
}

impl Default for Crosslink {
    fn default() -> Self {
        let (site1, site2) = Reagent::Edc.sites().unwrap();
        Self {
            enabled: false,
            reagent: Reagent::Dss,
            custom_mass: 138.06808,
            custom_cleavable: false,
            stub_short: 54.01056,
            stub_long: 103.9932,
            stub_alt: String::new(),
            asymmetric: false,
            stub2_short: 163.06333,
            stub2_long: 213.04596,
            stub2_alt: String::new(),
            residues: "K".into(),
            two_sites: false,
            site1,
            site2,
            protein_nterm: true,
            allow_cterm_residue: false,
            loop_links: true,
            candidates: 30,
            min_chain_peaks: 3,
            require_complete_doublets: false,
            backend: Backend::Cpu,
            chain_survivors: 20000,
            linear_search: false,
            min_linked_peaks: 2,
        }
    }
}

impl Crosslink {
    /// Is the selected linker MS-cleavable (a cleavable preset, or custom with stubs)?
    pub fn is_cleavable(&self) -> bool {
        self.reagent.is_cleavable() || (self.reagent == Reagent::Custom && self.custom_cleavable)
    }

    /// Mass the linker adds to the peptide pair (Da)
    pub fn linker_mass(&self) -> f64 {
        self.reagent.mass().unwrap_or(self.custom_mass)
    }

    /// Does the linker leave different stubs on its two arms?
    pub fn has_second_arm(&self) -> bool {
        match self.reagent {
            Reagent::Custom => self.custom_cleavable && self.asymmetric,
            r => r.second_arm().is_some(),
        }
    }

    /// Load the selected preset's stub masses (no change for non-cleavable / custom)
    pub fn reset_stubs(&mut self) {
        let list = |alt: &[f64]| alt.iter().map(f64::to_string).collect::<Vec<_>>().join(", ");
        if let Some((short, long, alt)) = self.reagent.stubs() {
            self.stub_short = short;
            self.stub_long = long;
            self.stub_alt = list(alt);
        }
        if let Some((short, long, alt)) = self.reagent.second_arm() {
            self.stub2_short = short;
            self.stub2_long = long;
            self.stub2_alt = list(alt);
        }
    }

    /// Have the stub masses of a cleavable preset been changed from its defaults?
    pub fn stubs_edited(&self) -> bool {
        let Some(first) = self.reagent.stubs() else {
            return false;
        };
        let same = |a: f64, b: f64| (a - b).abs() < 1e-6;
        let arm_same = |(short, long, alt): (f64, f64, &[f64]), (s, l, a): (f64, f64, &str)| {
            same(s, short)
                && same(l, long)
                && parse_list(a).is_ok_and(|v| v.len() == alt.len() && v.iter().zip(alt).all(|(x, y)| same(*x, *y)))
        };
        let first_same = arm_same(first, (self.stub_short, self.stub_long, &self.stub_alt));
        let second_same = self
            .reagent
            .second_arm()
            .is_none_or(|second| arm_same(second, (self.stub2_short, self.stub2_long, &self.stub2_alt)));
        !(first_same && second_same)
    }

    fn stubs_json(&self) -> Value {
        let mut stubs = json!({
            "short": self.stub_short,
            "long": self.stub_long,
            "alt": parse_list(&self.stub_alt).unwrap_or_default(),
        });
        if self.has_second_arm() {
            stubs["other"] = json!({
                "short": self.stub2_short,
                "long": self.stub2_long,
                "alt": parse_list(&self.stub2_alt).unwrap_or_default(),
            });
        }
        stubs
    }

    fn to_json(&self) -> Value {
        let mut m = Map::new();
        if self.reagent == Reagent::Custom {
            m.insert("mass".into(), json!(self.custom_mass));
            if self.custom_cleavable {
                m.insert("cleavable".into(), self.stubs_json());
            }
        } else {
            m.insert("reagent".into(), json!(self.reagent.key()));
            // Preset stubs are implied by the reagent; write them only when edited
            if self.stubs_edited() {
                m.insert("cleavable".into(), self.stubs_json());
            }
        }
        if self.two_sites {
            m.insert("site1".into(), json!(self.site1.to_spec()));
            m.insert("site2".into(), json!(self.site2.to_spec()));
        } else {
            m.insert("residues".into(), json!(self.residues));
            m.insert("protein_nterm".into(), json!(self.protein_nterm));
        }
        m.insert("allow_cterm_residue".into(), json!(self.allow_cterm_residue));
        m.insert("loop_links".into(), json!(self.loop_links));
        m.insert("candidates".into(), json!(self.candidates));
        m.insert("min_chain_peaks".into(), json!(self.min_chain_peaks));
        m.insert(
            "require_complete_doublets".into(),
            json!(self.require_complete_doublets),
        );
        m.insert("backend".into(), json!(self.backend.key()));
        m.insert("chain_survivors".into(), json!(self.chain_survivors));
        // only when on, so configs still run in a sage.exe without this option
        if self.linear_search {
            m.insert("linear_search".into(), json!(true));
        }
        // Sage's default (2) is left out, so configs still run in an older sage.exe
        if self.min_linked_peaks != 2 {
            m.insert("min_linked_peaks".into(), json!(self.min_linked_peaks));
        }
        Value::Object(m)
    }

    fn from_json(v: &Map<String, Value>) -> Self {
        let mut c = Crosslink {
            enabled: true,
            ..Default::default()
        };
        if let Some(r) = v.get("reagent").and_then(Value::as_str) {
            if let Some(r) = Reagent::ALL.iter().find(|x| x.key() == r) {
                c.reagent = *r;
            }
        } else if let Some(mass) = v.get("mass").and_then(Value::as_f64) {
            c.reagent = Reagent::Custom;
            c.custom_mass = mass;
        }
        c.reset_stubs();
        if let Some(stubs) = v.get("cleavable").and_then(Value::as_object) {
            c.custom_cleavable = c.reagent == Reagent::Custom;
            c.stub_short = f64_or(stubs.get("short"), c.stub_short);
            c.stub_long = f64_or(stubs.get("long"), c.stub_long);
            let list = |v: Option<&Value>| {
                v.and_then(Value::as_array)
                    .map(|a| a.iter().filter_map(Value::as_f64).map(|x| x.to_string()).collect::<Vec<_>>().join(", "))
                    .unwrap_or_default()
            };
            // Sage's default `alt` is empty, not the preset's
            c.stub_alt = list(stubs.get("alt"));
            if let Some(other) = stubs.get("other").and_then(Value::as_object) {
                c.asymmetric = true;
                c.stub2_short = f64_or(other.get("short"), c.stub2_short);
                c.stub2_long = f64_or(other.get("long"), c.stub2_long);
                c.stub2_alt = list(other.get("alt"));
            }
        }
        if let Some(r) = v.get("residues").and_then(Value::as_str) {
            c.residues = r.into();
        }
        c.protein_nterm = bool_or(v.get("protein_nterm"), c.protein_nterm);
        if let (Some(s1), Some(s2)) = (
            v.get("site1").and_then(Value::as_str),
            v.get("site2").and_then(Value::as_str),
        ) {
            c.two_sites = true;
            c.site1 = Site::from_spec(s1);
            c.site2 = Site::from_spec(s2);
        } else if let Some((s1, s2)) = c.reagent.sites() {
            c.two_sites = true;
            c.site1 = s1;
            c.site2 = s2;
        }
        c.allow_cterm_residue = bool_or(v.get("allow_cterm_residue"), c.allow_cterm_residue);
        c.loop_links = bool_or(v.get("loop_links"), c.loop_links);
        c.candidates = u32_or(v.get("candidates"), c.candidates);
        c.min_chain_peaks = u32_or(v.get("min_chain_peaks"), c.min_chain_peaks);
        c.require_complete_doublets =
            bool_or(v.get("require_complete_doublets"), c.require_complete_doublets);
        if let Some(b) = v.get("backend").and_then(Value::as_str) {
            if let Some(b) = Backend::ALL.iter().find(|x| x.key() == b) {
                c.backend = *b;
            }
        }
        c.chain_survivors = u32_or(v.get("chain_survivors"), c.chain_survivors);
        c.linear_search = bool_or(v.get("linear_search"), c.linear_search);
        c.min_linked_peaks = u32_or(v.get("min_linked_peaks"), c.min_linked_peaks);
        c
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tmt {
    None,
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
}

impl Tmt {
    pub const ALL: [Tmt; 6] = [Tmt::None, Tmt::Tmt6, Tmt::Tmt10, Tmt::Tmt11, Tmt::Tmt16, Tmt::Tmt18];

    /// Name used in the Sage JSON, `None` for no TMT quantification
    pub fn key(self) -> Option<&'static str> {
        match self {
            Tmt::None => None,
            Tmt::Tmt6 => Some("Tmt6"),
            Tmt::Tmt10 => Some("Tmt10"),
            Tmt::Tmt11 => Some("Tmt11"),
            Tmt::Tmt16 => Some("Tmt16"),
            Tmt::Tmt18 => Some("Tmt18"),
        }
    }

    pub fn label(self) -> &'static str {
        self.key().unwrap_or("None")
    }
}

pub const ION_KINDS: [&str; 7] = ["a", "b", "c", "x", "y", "z", "z_radical"];

/// A protease Sage knows by name (`sage_core::enzyme::protease`)
pub struct Protease {
    pub name: &'static str,
    pub label: &'static str,
    pub cleave_at: &'static str,
    pub restrict: &'static str,
    pub c_terminal: bool,
    pub help: &'static str,
}

pub const PROTEASES: &[Protease] = &[
    Protease { name: "trypsin", label: "Trypsin", cleave_at: "KR", restrict: "P", c_terminal: true,
        help: "Cleaves after K and R, not before P. The standard enzyme for bottom-up proteomics." },
    Protease { name: "trypsin/p", label: "Trypsin/P", cleave_at: "KR", restrict: "", c_terminal: true,
        help: "Cleaves after K and R, also before P." },
    Protease { name: "lys-c", label: "Lys-C", cleave_at: "K", restrict: "P", c_terminal: true,
        help: "Cleaves after K, not before P. Often combined with trypsin for more complete digestion." },
    Protease { name: "lys-n", label: "Lys-N", cleave_at: "K", restrict: "", c_terminal: false,
        help: "Cleaves before K (N-terminal side)." },
    Protease { name: "arg-c", label: "Arg-C", cleave_at: "R", restrict: "P", c_terminal: true,
        help: "Cleaves after R, not before P." },
    Protease { name: "asp-n", label: "Asp-N", cleave_at: "D", restrict: "", c_terminal: false,
        help: "Cleaves before D (N-terminal side). Combined with trypsin for sequential digests, e.g. in XL-MS." },
    Protease { name: "glu-c", label: "Glu-C", cleave_at: "E", restrict: "P", c_terminal: true,
        help: "Cleaves after E, not before P (ammonium bicarbonate buffer). In phosphate buffer it also \
               cuts after D: add a custom rule cleave at D." },
    Protease { name: "chymotrypsin", label: "Chymotrypsin", cleave_at: "FYWL", restrict: "P", c_terminal: true,
        help: "Cleaves after F, Y, W and L, not before P (low specificity)." },
    Protease { name: "chymotrypsin-high", label: "Chymotrypsin (high specificity)", cleave_at: "FYW", restrict: "P", c_terminal: true,
        help: "Cleaves after F, Y and W only, not before P." },
    Protease { name: "proalanase", label: "ProAlanase", cleave_at: "PA", restrict: "", c_terminal: true,
        help: "Cleaves after P and A. Useful for proline-rich proteins and histones." },
];

/// A hand-written cleavage rule
#[derive(Clone, Debug, PartialEq)]
pub struct EnzymeRule {
    pub cleave_at: String,
    pub restrict: String,
    pub c_terminal: bool,
}

/// Lower-case alphanumerics only, so `Asp-N` matches `aspn`
fn normalize(name: &str) -> String {
    name.chars()
        .filter(char::is_ascii_alphanumeric)
        .collect::<String>()
        .to_ascii_lowercase()
}

/// Sort a rule into a named protease when one matches it exactly
fn classify(rule: EnzymeRule, names: &mut Vec<String>, custom: &mut Vec<EnzymeRule>) {
    let preset = PROTEASES.iter().find(|p| {
        p.cleave_at == rule.cleave_at && p.restrict == rule.restrict && p.c_terminal == rule.c_terminal
    });
    match preset {
        Some(p) if !names.iter().any(|n| n == p.name) => names.push(p.name.into()),
        Some(_) => {}
        None => custom.push(rule),
    }
}

/// `(named proteases, custom rules, non-specific)` from a Sage `enzyme` block
fn parse_enzymes(enzyme: Option<&Map<String, Value>>) -> (Vec<String>, Vec<EnzymeRule>, bool) {
    let (mut names, mut custom) = (Vec::new(), Vec::new());
    let Some(e) = enzyme else {
        // Sage's default enzyme is trypsin
        return (vec!["trypsin".into()], custom, false);
    };
    if let Some(list) = e.get("enzymes").and_then(Value::as_array) {
        for item in list {
            match item {
                Value::String(name) => {
                    if let Some(p) = PROTEASES.iter().find(|p| normalize(p.name) == normalize(name)) {
                        if !names.iter().any(|n| n == p.name) {
                            names.push(p.name.into());
                        }
                    }
                }
                Value::Object(o) => classify(
                    EnzymeRule {
                        cleave_at: str_or(o.get("cleave_at"), ""),
                        restrict: str_or(o.get("restrict"), ""),
                        c_terminal: bool_or(o.get("c_terminal"), true),
                    },
                    &mut names,
                    &mut custom,
                ),
                _ => {}
            }
        }
        return (names, custom, false);
    }
    // Single enzyme, Sage's defaults for missing keys
    let rule = EnzymeRule {
        cleave_at: str_or(e.get("cleave_at"), "KR"),
        restrict: str_or(e.get("restrict"), ""),
        c_terminal: bool_or(e.get("c_terminal"), true),
    };
    if rule.cleave_at.is_empty() {
        return (names, custom, true);
    }
    classify(rule, &mut names, &mut custom);
    (names, custom, false)
}

#[derive(Clone, Debug)]
pub struct Settings {
    // Files
    pub fasta: String,
    pub spectra: Vec<String>,
    pub output_directory: String,

    // Database & enzyme
    pub missed_cleavages: u32,
    pub min_len: u32,
    pub max_len: u32,
    /// Named proteases used together (co-digestion), see [`PROTEASES`]
    pub enzymes: Vec<String>,
    /// Further hand-written cleavage rules
    pub custom_enzymes: Vec<EnzymeRule>,
    /// Non-specific search: every peptide bond
    pub nonspecific: bool,
    pub semi_enzymatic: bool,
    /// `None` lets Sage choose (500 Da, or 350 Da for crosslink searches)
    pub peptide_min_mass: Option<f64>,
    pub peptide_max_mass: f64,
    pub ion_kinds: [bool; 7],
    pub min_ion_index: u32,
    pub max_variable_mods: u32,
    pub max_combinations: Option<u32>,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    pub bucket_size: u32,
    pub prefilter: bool,
    pub prefilter_chunk_size: u32,
    pub prefilter_low_memory: bool,

    // Modifications
    pub mods: Vec<Modification>,

    // Spectra & tolerances
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: u32,
    /// `None` lets Sage choose (150, or 300 for cleavable crosslinkers)
    pub max_peaks: Option<u32>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u32,
    pub report_psms: u32,

    // Scoring, quantification & output
    pub openms_score: bool,
    pub predict_rt: bool,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f64,
    pub write_pin: bool,
    pub write_report: bool,
    /// `results.annotator.csv` / `crosslinks.annotator.csv` for Annotator
    pub write_annotator: bool,
    pub annotate_matches: bool,
    pub lfq: bool,
    /// MS1 mass tolerance (ppm) for LFQ
    pub lfq_ppm: f64,
    /// Only peptides / crosslinks at or below this q-value are quantified
    pub lfq_q: f64,
    pub lfq_combine_charges: bool,
    /// Minimum isotope-envelope spectral angle of an integrated peak
    pub lfq_spectral_angle: f64,
    pub tmt: Tmt,
    pub tmt_level: u8,
    pub tmt_sn: bool,

    pub crosslink: Crosslink,

    /// Keys from a loaded file the GUI does not edit, written back unchanged
    extra: Map<String, Value>,
    extra_database: Map<String, Value>,
    extra_enzyme: Map<String, Value>,
    extra_quant: Map<String, Value>,
    extra_lfq: Map<String, Value>,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            fasta: String::new(),
            spectra: Vec::new(),
            output_directory: String::new(),
            missed_cleavages: 2,
            min_len: 5,
            max_len: 50,
            enzymes: vec!["trypsin".into()],
            custom_enzymes: Vec::new(),
            nonspecific: false,
            semi_enzymatic: false,
            peptide_min_mass: None,
            peptide_max_mass: 5000.0,
            ion_kinds: [false, true, false, false, true, false, false],
            min_ion_index: 2,
            max_variable_mods: 2,
            max_combinations: None,
            generate_decoys: true,
            decoy_tag: "rev_".into(),
            bucket_size: 8192,
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
            mods: vec![(&PRESETS[0]).into(), (&PRESETS[2]).into()],
            precursor_tol: Tolerance::ppm(-10.0, 10.0),
            fragment_tol: Tolerance::ppm(-10.0, 10.0),
            precursor_charge: (2, 4),
            override_precursor_charge: false,
            isotope_errors: (0, 0),
            deisotope: true,
            chimera: false,
            wide_window: false,
            min_peaks: 15,
            max_peaks: None,
            max_fragment_charge: None,
            min_matched_peaks: 4,
            report_psms: 1,
            openms_score: false,
            predict_rt: true,
            protein_grouping: true,
            protein_grouping_peptide_fdr: 0.01,
            write_pin: false,
            write_report: false,
            write_annotator: false,
            annotate_matches: false,
            lfq: false,
            lfq_ppm: 5.0,
            lfq_q: 0.01,
            lfq_combine_charges: true,
            lfq_spectral_angle: 0.7,
            tmt: Tmt::None,
            tmt_level: 3,
            tmt_sn: false,
            crosslink: Crosslink::default(),
            extra: Map::new(),
            extra_database: Map::new(),
            extra_enzyme: Map::new(),
            extra_quant: Map::new(),
            extra_lfq: Map::new(),
        }
    }
}

const TOP_KEYS: &[&str] = &[
    "database", "precursor_tol", "fragment_tol", "precursor_charge", "override_precursor_charge",
    "isotope_errors", "deisotope", "chimera", "wide_window", "min_peaks", "max_peaks",
    "max_fragment_charge", "min_matched_peaks", "report_psms", "score_type", "predict_rt",
    "protein_grouping", "protein_grouping_peptide_fdr", "write_pin", "write_report",
    "write_annotator", "annotate_matches", "quant", "crosslink", "mzml_paths", "output_directory",
];
const DATABASE_KEYS: &[&str] = &[
    "fasta", "enzyme", "peptide_min_mass", "peptide_max_mass", "ion_kinds", "min_ion_index",
    "static_mods", "variable_mods", "max_variable_mods", "max_combinations", "generate_decoys",
    "decoy_tag", "bucket_size", "prefilter", "prefilter_chunk_size", "prefilter_low_memory",
];
const ENZYME_KEYS: &[&str] = &[
    "missed_cleavages", "min_len", "max_len", "cleave_at", "restrict", "c_terminal", "enzymes",
    "semi_enzymatic",
];
const QUANT_KEYS: &[&str] = &["lfq", "tmt", "tmt_settings", "lfq_settings"];
const LFQ_KEYS: &[&str] = &["ppm_tolerance", "peptide_q_value", "combine_charge_states", "spectral_angle"];

fn f64_or(v: Option<&Value>, default: f64) -> f64 {
    v.and_then(Value::as_f64).unwrap_or(default)
}

fn u32_or(v: Option<&Value>, default: u32) -> u32 {
    v.and_then(Value::as_u64).map(|x| x as u32).unwrap_or(default)
}

fn bool_or(v: Option<&Value>, default: bool) -> bool {
    v.and_then(Value::as_bool).unwrap_or(default)
}

fn str_or(v: Option<&Value>, default: &str) -> String {
    v.and_then(Value::as_str).unwrap_or(default).to_string()
}

fn pair<T: TryFrom<i64>>(v: Option<&Value>, default: (T, T)) -> (T, T) {
    let get = |i: usize| {
        v.and_then(|v| v.get(i))
            .and_then(Value::as_i64)
            .and_then(|x| T::try_from(x).ok())
    };
    match (get(0), get(1)) {
        (Some(a), Some(b)) => (a, b),
        _ => default,
    }
}

/// Keys of `obj` that are not in `known`
fn unknown_keys(obj: Option<&Map<String, Value>>, known: &[&str]) -> Map<String, Value> {
    obj.map(|o| {
        o.iter()
            .filter(|(k, _)| !known.contains(&k.as_str()))
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect()
    })
    .unwrap_or_default()
}

/// Parse a comma separated list of numbers
pub fn parse_list(s: &str) -> Result<Vec<f64>, String> {
    s.split(',')
        .map(str::trim)
        .filter(|x| !x.is_empty())
        .map(|x| x.parse::<f64>().map_err(|_| format!("`{x}` is not a number")))
        .collect()
}

fn mod_entry(m: &Modification) -> Value {
    let mut e = Map::new();
    e.insert("mass".into(), json!(m.mass));
    if m.kind == ModKind::Variable && m.max_count > 0 {
        e.insert("max_count".into(), json!(m.max_count));
    }
    if !m.name.trim().is_empty() {
        e.insert("name".into(), json!(m.name.trim()));
    }
    if e.len() == 1 {
        json!(m.mass)
    } else {
        Value::Object(e)
    }
}

fn parse_mod_entry(target: &str, v: &Value, kind: ModKind) -> Option<Modification> {
    let (mass, name, max_count) = match v {
        Value::Number(n) => (n.as_f64()?, String::new(), 0),
        Value::Object(o) => (
            o.get("mass")?.as_f64()?,
            str_or(o.get("name"), ""),
            u32_or(o.get("max_count"), 0),
        ),
        _ => return None,
    };
    Some(Modification {
        target: target.into(),
        mass,
        name,
        kind,
        max_count,
    })
}

impl Settings {
    pub fn is_thermo_raw(path: &str) -> bool {
        path.to_ascii_lowercase().ends_with(".raw") && std::path::Path::new(path).is_file()
    }

    /// Build the Sage JSON configuration. `spectra` replaces the file list
    /// (used to substitute converted mzML files for Thermo .raw files).
    pub fn to_json(&self, spectra: &[String]) -> Value {
        let mut enzyme = self.extra_enzyme.clone();
        enzyme.insert("missed_cleavages".into(), json!(self.missed_cleavages));
        enzyme.insert("min_len".into(), json!(self.min_len));
        enzyme.insert("max_len".into(), json!(self.max_len));
        match (self.nonspecific, &self.enzymes[..], &self.custom_enzymes[..]) {
            (true, _, _) => {
                enzyme.insert("cleave_at".into(), json!(""));
            }
            // A single hand-written rule keeps the classic form, which also
            // allows `$` (no cleavage)
            (false, [], [rule]) => {
                enzyme.insert("cleave_at".into(), json!(rule.cleave_at));
                enzyme.insert("restrict".into(), json!(rule.restrict));
                enzyme.insert("c_terminal".into(), json!(rule.c_terminal));
            }
            (false, names, rules) => {
                let list = names
                    .iter()
                    .map(|n| json!(n))
                    .chain(rules.iter().map(|r| {
                        json!({
                            "cleave_at": r.cleave_at,
                            "restrict": r.restrict,
                            "c_terminal": r.c_terminal,
                        })
                    }))
                    .collect::<Vec<_>>();
                enzyme.insert("enzymes".into(), Value::Array(list));
            }
        }
        enzyme.insert("semi_enzymatic".into(), json!(self.semi_enzymatic));

        let mut static_mods = Map::new();
        let mut variable_mods = Map::new();
        for m in &self.mods {
            match m.kind {
                ModKind::Fixed => {
                    static_mods.insert(m.target.clone(), mod_entry(m));
                }
                ModKind::Variable => {
                    let list = variable_mods
                        .entry(m.target.clone())
                        .or_insert_with(|| json!([]));
                    list.as_array_mut().unwrap().push(mod_entry(m));
                }
            }
        }

        let ion_kinds: Vec<&str> = ION_KINDS
            .iter()
            .zip(self.ion_kinds)
            .filter(|(_, on)| *on)
            .map(|(k, _)| *k)
            .collect();

        let mut db = self.extra_database.clone();
        db.insert("fasta".into(), json!(self.fasta));
        db.insert("enzyme".into(), Value::Object(enzyme));
        if let Some(mass) = self.peptide_min_mass {
            db.insert("peptide_min_mass".into(), json!(mass));
        }
        db.insert("peptide_max_mass".into(), json!(self.peptide_max_mass));
        db.insert("ion_kinds".into(), json!(ion_kinds));
        db.insert("min_ion_index".into(), json!(self.min_ion_index));
        db.insert("static_mods".into(), Value::Object(static_mods));
        db.insert("variable_mods".into(), Value::Object(variable_mods));
        db.insert("max_variable_mods".into(), json!(self.max_variable_mods));
        if let Some(n) = self.max_combinations {
            db.insert("max_combinations".into(), json!(n));
        }
        db.insert("generate_decoys".into(), json!(self.generate_decoys));
        db.insert("decoy_tag".into(), json!(self.decoy_tag));
        db.insert("bucket_size".into(), json!(self.bucket_size));
        db.insert("prefilter".into(), json!(self.prefilter));
        if self.prefilter {
            db.insert("prefilter_chunk_size".into(), json!(self.prefilter_chunk_size));
            db.insert("prefilter_low_memory".into(), json!(self.prefilter_low_memory));
        }

        let mut quant = self.extra_quant.clone();
        quant.insert("lfq".into(), json!(self.lfq));
        let mut lfq = self.extra_lfq.clone();
        lfq.insert("ppm_tolerance".into(), json!(self.lfq_ppm));
        lfq.insert("peptide_q_value".into(), json!(self.lfq_q));
        lfq.insert("combine_charge_states".into(), json!(self.lfq_combine_charges));
        lfq.insert("spectral_angle".into(), json!(self.lfq_spectral_angle));
        quant.insert("lfq_settings".into(), Value::Object(lfq));
        if let Some(tmt) = self.tmt.key() {
            quant.insert("tmt".into(), json!(tmt));
            quant.insert(
                "tmt_settings".into(),
                json!({ "level": self.tmt_level, "sn": self.tmt_sn }),
            );
        }

        let mut top = self.extra.clone();
        top.insert("database".into(), Value::Object(db));
        top.insert("precursor_tol".into(), self.precursor_tol.to_json());
        top.insert("fragment_tol".into(), self.fragment_tol.to_json());
        top.insert(
            "precursor_charge".into(),
            json!([self.precursor_charge.0, self.precursor_charge.1]),
        );
        top.insert(
            "override_precursor_charge".into(),
            json!(self.override_precursor_charge),
        );
        top.insert(
            "isotope_errors".into(),
            json!([self.isotope_errors.0, self.isotope_errors.1]),
        );
        top.insert("deisotope".into(), json!(self.deisotope));
        top.insert("chimera".into(), json!(self.chimera));
        top.insert("wide_window".into(), json!(self.wide_window));
        top.insert("min_peaks".into(), json!(self.min_peaks));
        if let Some(n) = self.max_peaks {
            top.insert("max_peaks".into(), json!(n));
        }
        if let Some(z) = self.max_fragment_charge {
            top.insert("max_fragment_charge".into(), json!(z));
        }
        top.insert("min_matched_peaks".into(), json!(self.min_matched_peaks));
        top.insert("report_psms".into(), json!(self.report_psms));
        top.insert(
            "score_type".into(),
            json!(if self.openms_score { "OpenMSHyperScore" } else { "SageHyperScore" }),
        );
        top.insert("predict_rt".into(), json!(self.predict_rt));
        top.insert("protein_grouping".into(), json!(self.protein_grouping));
        top.insert(
            "protein_grouping_peptide_fdr".into(),
            json!(self.protein_grouping_peptide_fdr),
        );
        top.insert("write_pin".into(), json!(self.write_pin));
        top.insert("write_report".into(), json!(self.write_report));
        top.insert("write_annotator".into(), json!(self.write_annotator));
        top.insert("annotate_matches".into(), json!(self.annotate_matches));
        top.insert("quant".into(), Value::Object(quant));
        if self.crosslink.enabled {
            top.insert("crosslink".into(), self.crosslink.to_json());
        }
        top.insert("mzml_paths".into(), json!(spectra));
        top.insert("output_directory".into(), json!(self.output_directory));
        Value::Object(top)
    }

    /// Read a Sage JSON configuration. Missing keys take Sage's defaults.
    pub fn from_json(v: &Value) -> Result<Self, String> {
        let top = v.as_object().ok_or("the file is not a JSON object")?;
        let db = top.get("database").and_then(Value::as_object);
        let enzyme = db.and_then(|d| d.get("enzyme")).and_then(Value::as_object);
        let quant = top.get("quant").and_then(Value::as_object);
        let d = Settings::default();
        let dbv = |k: &str| db.and_then(|d| d.get(k));
        let ev = |k: &str| enzyme.and_then(|e| e.get(k));
        let lfq_settings = quant.and_then(|q| q.get("lfq_settings")).and_then(Value::as_object);
        let lv = |k: &str| lfq_settings.and_then(|l| l.get(k));

        let mut mods = Vec::new();
        if let Some(fixed) = dbv("static_mods").and_then(Value::as_object) {
            for (target, entry) in fixed {
                mods.extend(parse_mod_entry(target, entry, ModKind::Fixed));
            }
        }
        if let Some(var) = dbv("variable_mods").and_then(Value::as_object) {
            for (target, entries) in var {
                match entries {
                    Value::Array(list) => mods.extend(
                        list.iter()
                            .filter_map(|e| parse_mod_entry(target, e, ModKind::Variable)),
                    ),
                    e => mods.extend(parse_mod_entry(target, e, ModKind::Variable)),
                }
            }
        }

        let mut ion_kinds = d.ion_kinds;
        if let Some(kinds) = dbv("ion_kinds").and_then(Value::as_array) {
            ion_kinds = [false; 7];
            for k in kinds.iter().filter_map(Value::as_str) {
                if let Some(i) = ION_KINDS.iter().position(|x| *x == k) {
                    ion_kinds[i] = true;
                }
            }
        }

        let (enzymes, custom_enzymes, nonspecific) = parse_enzymes(enzyme);

        let tmt = quant
            .and_then(|q| q.get("tmt"))
            .and_then(Value::as_str)
            .and_then(|t| Tmt::ALL.iter().find(|x| x.key() == Some(t)).copied())
            .unwrap_or(Tmt::None);
        let tmt_settings = quant
            .and_then(|q| q.get("tmt_settings"))
            .and_then(Value::as_object);

        Ok(Settings {
            fasta: str_or(dbv("fasta"), ""),
            spectra: top
                .get("mzml_paths")
                .and_then(Value::as_array)
                .map(|a| a.iter().filter_map(Value::as_str).map(String::from).collect())
                .unwrap_or_default(),
            output_directory: str_or(top.get("output_directory"), ""),
            missed_cleavages: u32_or(ev("missed_cleavages"), 0),
            min_len: u32_or(ev("min_len"), d.min_len),
            max_len: u32_or(ev("max_len"), d.max_len),
            enzymes,
            custom_enzymes,
            nonspecific,
            semi_enzymatic: bool_or(ev("semi_enzymatic"), d.semi_enzymatic),
            peptide_min_mass: dbv("peptide_min_mass").and_then(Value::as_f64),
            peptide_max_mass: f64_or(dbv("peptide_max_mass"), d.peptide_max_mass),
            ion_kinds,
            min_ion_index: u32_or(dbv("min_ion_index"), d.min_ion_index),
            max_variable_mods: u32_or(dbv("max_variable_mods"), d.max_variable_mods),
            max_combinations: dbv("max_combinations").and_then(Value::as_u64).map(|x| x as u32),
            generate_decoys: bool_or(dbv("generate_decoys"), d.generate_decoys),
            decoy_tag: str_or(dbv("decoy_tag"), &d.decoy_tag),
            bucket_size: u32_or(dbv("bucket_size"), d.bucket_size),
            prefilter: bool_or(dbv("prefilter"), d.prefilter),
            prefilter_chunk_size: u32_or(dbv("prefilter_chunk_size"), d.prefilter_chunk_size),
            prefilter_low_memory: bool_or(dbv("prefilter_low_memory"), d.prefilter_low_memory),
            mods,
            precursor_tol: top
                .get("precursor_tol")
                .and_then(Tolerance::from_json)
                .unwrap_or(d.precursor_tol),
            fragment_tol: top
                .get("fragment_tol")
                .and_then(Tolerance::from_json)
                .unwrap_or(d.fragment_tol),
            precursor_charge: pair(top.get("precursor_charge"), d.precursor_charge),
            override_precursor_charge: bool_or(
                top.get("override_precursor_charge"),
                d.override_precursor_charge,
            ),
            isotope_errors: pair(top.get("isotope_errors"), d.isotope_errors),
            deisotope: bool_or(top.get("deisotope"), d.deisotope),
            chimera: bool_or(top.get("chimera"), d.chimera),
            wide_window: bool_or(top.get("wide_window"), d.wide_window),
            min_peaks: u32_or(top.get("min_peaks"), d.min_peaks),
            max_peaks: top.get("max_peaks").and_then(Value::as_u64).map(|x| x as u32),
            max_fragment_charge: top
                .get("max_fragment_charge")
                .and_then(Value::as_u64)
                .map(|x| x as u8),
            min_matched_peaks: u32_or(top.get("min_matched_peaks"), d.min_matched_peaks),
            report_psms: u32_or(top.get("report_psms"), d.report_psms),
            openms_score: top.get("score_type").and_then(Value::as_str)
                == Some("OpenMSHyperScore"),
            predict_rt: bool_or(top.get("predict_rt"), d.predict_rt),
            protein_grouping: bool_or(top.get("protein_grouping"), d.protein_grouping),
            protein_grouping_peptide_fdr: f64_or(
                top.get("protein_grouping_peptide_fdr"),
                d.protein_grouping_peptide_fdr,
            ),
            write_pin: bool_or(top.get("write_pin"), d.write_pin),
            write_report: bool_or(top.get("write_report"), d.write_report),
            // absent from a file means off, as in Sage
            write_annotator: bool_or(top.get("write_annotator"), false),
            annotate_matches: bool_or(top.get("annotate_matches"), d.annotate_matches),
            lfq: bool_or(quant.and_then(|q| q.get("lfq")), false),
            lfq_ppm: f64_or(lv("ppm_tolerance"), d.lfq_ppm),
            lfq_q: f64_or(lv("peptide_q_value"), d.lfq_q),
            lfq_combine_charges: bool_or(lv("combine_charge_states"), d.lfq_combine_charges),
            lfq_spectral_angle: f64_or(lv("spectral_angle"), d.lfq_spectral_angle),
            tmt,
            tmt_level: tmt_settings
                .and_then(|t| t.get("level"))
                .and_then(Value::as_u64)
                .map(|x| x as u8)
                .unwrap_or(d.tmt_level),
            tmt_sn: bool_or(tmt_settings.and_then(|t| t.get("sn")), d.tmt_sn),
            crosslink: top
                .get("crosslink")
                .and_then(Value::as_object)
                .map(Crosslink::from_json)
                .unwrap_or_default(),
            extra: unknown_keys(Some(top), TOP_KEYS),
            extra_database: unknown_keys(db, DATABASE_KEYS),
            extra_enzyme: unknown_keys(enzyme, ENZYME_KEYS),
            extra_quant: unknown_keys(quant, QUANT_KEYS),
            extra_lfq: unknown_keys(lfq_settings, LFQ_KEYS),
        })
    }

    /// Problems that stop a search (`errors`) and ones worth a look (`warnings`)
    pub fn check(&self) -> (Vec<String>, Vec<String>) {
        let mut errors = Vec::new();
        let mut warnings = Vec::new();

        if self.fasta.trim().is_empty() {
            errors.push("Select a FASTA database (Files tab).".into());
        } else if !std::path::Path::new(&self.fasta).is_file() {
            errors.push(format!("FASTA file not found: {}", self.fasta));
        }
        if self.spectra.is_empty() {
            errors.push("Add at least one spectrum file or Bruker .d folder (Files tab).".into());
        }
        for f in &self.spectra {
            if !std::path::Path::new(f).exists() {
                errors.push(format!("Spectrum file not found: {f}"));
            }
        }
        if self.output_directory.trim().is_empty() {
            errors.push("Select an output folder for the results (Files tab).".into());
        }

        if self.min_len > self.max_len {
            errors.push("Minimum peptide length is larger than the maximum.".into());
        }
        if !self.nonspecific {
            if self.enzymes.is_empty() && self.custom_enzymes.is_empty() {
                errors.push("Select at least one enzyme, or choose a non-specific search (Database tab).".into());
            }
            let single = self.enzymes.is_empty() && self.custom_enzymes.len() == 1;
            for rule in &self.custom_enzymes {
                let no_cleavage = single && rule.cleave_at == "$";
                if !no_cleavage
                    && (rule.cleave_at.is_empty() || !rule.cleave_at.chars().all(|c| VALID_AA.contains(c)))
                {
                    errors.push(format!(
                        "Custom enzyme: `cleave at` must be amino acid letters from {VALID_AA} \
                         (or `$` for no cleavage, as the only enzyme), got `{}`.",
                        rule.cleave_at
                    ));
                }
                if !rule.restrict.chars().all(|c| VALID_AA.contains(c)) {
                    errors.push(format!("Custom enzyme: `restrict` may only contain {VALID_AA}."));
                }
            }
        }
        if !self.ion_kinds.iter().any(|x| *x) {
            errors.push("Select at least one fragment ion type.".into());
        }

        let mut fixed_targets = std::collections::HashSet::new();
        for m in &self.mods {
            if let Err(e) = check_target(&m.target) {
                errors.push(format!("Modification target {e}"));
            }
            if m.kind == ModKind::Fixed && !fixed_targets.insert(m.target.clone()) {
                errors.push(format!(
                    "Two fixed modifications on `{}`: only one fixed modification per target is allowed.",
                    m.target
                ));
            }
        }

        for (name, tol) in [("Precursor", &self.precursor_tol), ("Fragment", &self.fragment_tol)] {
            if tol.lo > tol.hi {
                errors.push(format!("{name} tolerance: the lower bound is above the upper bound."));
            }
        }
        if self.precursor_charge.0 > self.precursor_charge.1 {
            errors.push("Precursor charge: minimum is above maximum.".into());
        }
        if self.isotope_errors.0 > self.isotope_errors.1 {
            errors.push("Isotope errors: minimum is above maximum.".into());
        }

        if self.crosslink.enabled {
            let xl = &self.crosslink;
            if self.prefilter {
                errors.push("Database prefilter cannot be used with crosslink searches.".into());
            }
            if self.wide_window {
                errors.push("Wide-window (DIA) mode is not supported for crosslink searches.".into());
            }
            if self.missed_cleavages == 0 {
                warnings.push(
                    "Crosslink search with 0 missed cleavages finds no linkable peptides; use 2.".into(),
                );
            }
            // A fixed modification on a reactive residue (typically TMT on K) blocks every link there
            let reactive = if xl.two_sites {
                format!("{}{}", xl.site1.residues, xl.site2.residues)
            } else {
                xl.residues.clone()
            };
            let blocked = self
                .mods
                .iter()
                .filter(|m| m.kind == ModKind::Fixed && m.target.len() == 1 && reactive.contains(m.target.as_str()))
                .map(|m| m.target.clone())
                .collect::<Vec<_>>();
            if !blocked.is_empty() {
                warnings.push(format!(
                    "Fixed modification on crosslink-reactive residue {}: those residues can't carry a link. \
                     For TMT applied after crosslinking, make TMT on K a Variable modification \
                     (keep it Fixed on the peptide N-terminus, ^).",
                    blocked.join(", ")
                ));
            }
            if xl.two_sites {
                for (name, site) in [("Site 1", &xl.site1), ("Site 2", &xl.site2)] {
                    if !site.residues.chars().all(|c| VALID_AA.contains(c)) {
                        errors.push(format!("Crosslink {name}: residues must be amino acids from {VALID_AA}."));
                    }
                    if site.residues.is_empty() && !site.protein_nterm && !site.protein_cterm {
                        errors.push(format!("Crosslink {name}: choose residues or a protein terminus."));
                    }
                }
            } else if xl.residues.is_empty() || !xl.residues.chars().all(|c| VALID_AA.contains(c)) {
                errors.push(format!("Crosslink residues must be amino acids from {VALID_AA}."));
            }
            if xl.is_cleavable() {
                let mut arms = vec![("", xl.stub_short, xl.stub_long, &xl.stub_alt)];
                if xl.has_second_arm() {
                    arms[0].0 = "arm 1 ";
                    arms.push(("arm 2 ", xl.stub2_short, xl.stub2_long, &xl.stub2_alt));
                }
                for &(arm, short, long, alt) in &arms {
                    if !(short > 0.0 && long > short) {
                        errors.push(format!("Cleavable linker: {arm}stub masses need 0 < short < long."));
                    }
                    match parse_list(alt) {
                        Err(e) => errors.push(format!("Additional {arm}stub masses: {e}")),
                        Ok(alt) if alt.iter().any(|&m| m <= 0.0) => {
                            errors.push(format!("Additional {arm}stub masses must be positive."))
                        }
                        Ok(_) => {}
                    }
                }
                // a cleavage leaves one arm's short stub and the other arm's long stub
                let last = arms.len() - 1;
                for (i, &(arm, short, ..)) in arms.iter().enumerate() {
                    let (other, _, long, _) = arms[last - i];
                    if (short + long - xl.linker_mass()).abs() > 1e-3 {
                        errors.push(format!(
                            "Cleavable linker: {arm}short + {other}long stub ({:.5} Da) must equal the linker mass ({:.5} Da).",
                            short + long,
                            xl.linker_mass()
                        ));
                    }
                }
            }
        }

        (errors, warnings)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn mod_targets() {
        for ok in ["C", "M", "^", "^Q", "$", "$K", "[", "[M", "]"] {
            assert!(check_target(ok).is_ok(), "{ok}");
        }
        for bad in ["", "B", "CC", "^QQ", "x", "^1"] {
            assert!(check_target(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn json_round_trip_keeps_settings_and_unknown_keys() {
        let mut s = Settings::default();
        s.mods.push((&PRESETS[4]).into());
        s.mods.last_mut().unwrap().max_count = 2;
        s.max_peaks = Some(250);
        s.fragment_tol = Tolerance { unit: TolUnit::Da, lo: -0.5, hi: 0.5 };
        s.tmt = Tmt::Tmt16;
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Custom;
        s.crosslink.custom_mass = 158.00376;
        s.crosslink.custom_cleavable = true;
        s.crosslink.stub_alt = "85.98264".into();
        s.crosslink.backend = Backend::Reflect;
        let mut v = s.to_json(&["a.mzML".to_string()]);
        v["bruker_config"] = json!({"ms2": {}});
        v["database"]["enzyme"]["future_key"] = json!(1);

        let back = Settings::from_json(&v).unwrap();
        let v2 = back.to_json(&back.spectra);
        // everything, including keys the GUI does not edit, must survive
        assert_eq!(v, v2);
        assert_eq!(v["database"]["variable_mods"]["S"], json!([{"mass": 79.966331, "max_count": 2, "name": "UNIMOD:21"}]));
        assert_eq!(v["database"]["static_mods"]["C"], json!({"mass": 57.021464, "name": "UNIMOD:4"}));
        assert!(v["database"].get("peptide_min_mass").is_none());
        assert_eq!(v["crosslink"]["cleavable"]["alt"], json!([85.98264]));
    }

    #[test]
    fn enzymes_and_two_site_chemistry_round_trip() {
        let mut s = Settings::default();
        s.enzymes = vec!["trypsin".into(), "asp-n".into()];
        s.custom_enzymes.push(EnzymeRule { cleave_at: "DE".into(), restrict: "P".into(), c_terminal: true });
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Edc;
        (s.crosslink.site1, s.crosslink.site2) = Reagent::Edc.sites().unwrap();
        s.crosslink.two_sites = true;
        let v = s.to_json(&[]);
        assert_eq!(
            v["database"]["enzyme"]["enzymes"],
            json!(["trypsin", "asp-n", {"cleave_at": "DE", "restrict": "P", "c_terminal": true}])
        );
        assert!(v["database"]["enzyme"].get("cleave_at").is_none());
        assert_eq!(v["crosslink"]["reagent"], json!("edc"));
        assert_eq!(v["crosslink"]["site1"], json!("K["));
        assert_eq!(v["crosslink"]["site2"], json!("DE]"));
        assert!(v["crosslink"].get("residues").is_none());
        assert_eq!(Settings::from_json(&v).unwrap().to_json(&[]), v);

        // a preset reagent alone implies its sites
        let edc = Settings::from_json(&json!({"crosslink": {"reagent": "dmtmm"}})).unwrap();
        assert!(edc.crosslink.two_sites);
        assert_eq!(edc.crosslink.site2.residues, "DE");
    }

    #[test]
    fn classic_enzyme_blocks_map_to_presets() {
        let enzymes = |e: Value| {
            let s = Settings::from_json(&json!({"database": {"enzyme": e}})).unwrap();
            (s.enzymes, s.custom_enzymes, s.nonspecific)
        };
        assert_eq!(enzymes(json!({"cleave_at": "KR", "restrict": "P"})).0, vec!["trypsin"]);
        assert_eq!(enzymes(json!({"cleave_at": "D", "c_terminal": false})).0, vec!["asp-n"]);
        // Sage's defaults inside an enzyme block: KR without restriction
        assert_eq!(enzymes(json!({"missed_cleavages": 2})).0, vec!["trypsin/p"]);
        assert!(enzymes(json!({"cleave_at": ""})).2);
        let (names, custom, _) = enzymes(json!({"cleave_at": "$"}));
        assert!(names.is_empty() && custom[0].cleave_at == "$");
        assert_eq!(enzymes(json!({"enzymes": ["Lys-C", "GluC"]})).0, vec!["lys-c", "glu-c"]);

        // a single custom rule is written back in the classic form (so `$` works)
        let s = Settings::from_json(&json!({"database": {"enzyme": {"cleave_at": "$"}}})).unwrap();
        assert_eq!(s.to_json(&[])["database"]["enzyme"]["cleave_at"], json!("$"));
        assert!(s.check().0.iter().all(|e| !e.contains("Custom enzyme")));
    }

    #[test]
    fn preset_stub_masses_are_editable() {
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Dsso;
        s.crosslink.reset_stubs();
        // unedited preset stubs are implied by the reagent
        let v = s.to_json(&[]);
        assert!(v["crosslink"].get("cleavable").is_none());
        let back = Settings::from_json(&v).unwrap();
        assert_eq!((back.crosslink.stub_short, back.crosslink.stub_long), (54.01056, 103.9932));
        assert_eq!(back.crosslink.stub_alt, "85.98264");
        assert!(!back.crosslink.stubs_edited());

        // dropping the sulfenic acid stub writes an override next to the reagent
        s.crosslink.stub_alt.clear();
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["reagent"], json!("dsso"));
        assert_eq!(v["crosslink"]["cleavable"], json!({"short": 54.01056, "long": 103.9932, "alt": []}));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.crosslink.reagent, Reagent::Dsso);
        assert!(back.crosslink.stub_alt.is_empty());
        assert_eq!(back.to_json(&[]), v);

        // stubs must still add up to the linker mass
        s.crosslink.stub_long = 100.0;
        assert!(s.check().0.iter().any(|e| e.contains("must equal the linker mass")));
    }

    #[test]
    fn asymmetric_linkers_have_two_arms() {
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Clixlink;
        s.crosslink.reset_stubs();
        assert!(s.crosslink.has_second_arm());
        assert!(s.check().0.iter().all(|e| !e.contains("Cleavable")), "{:?}", s.check().0);
        // the preset is implied by the reagent
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["reagent"], json!("clixlink"));
        assert!(v["crosslink"].get("cleavable").is_none());

        // an edited second arm writes both arms
        s.crosslink.stub2_alt.clear();
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["cleavable"]["other"], json!({"short": 163.06333, "long": 213.04596, "alt": []}));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.crosslink.reagent, Reagent::Clixlink);
        assert_eq!(back.crosslink.stub_alt, "85.98264");
        assert!(back.crosslink.stub2_alt.is_empty());
        assert_eq!(back.to_json(&[]), v);

        // each arm's short stub + the other arm's long stub must make the linker
        s.crosslink.stub2_long = 200.0;
        assert!(s.check().0.iter().any(|e| e.contains("arm 1 short + arm 2 long")));

        // custom asymmetric linkers
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Custom;
        s.crosslink.custom_mass = 267.05653;
        s.crosslink.custom_cleavable = true;
        s.crosslink.asymmetric = true;
        s.crosslink.stub_short = 54.01056;
        s.crosslink.stub_long = 103.9932;
        let v = s.to_json(&[]);
        assert!(v["crosslink"]["cleavable"].get("other").is_some());
        let back = Settings::from_json(&v).unwrap();
        assert!(back.crosslink.custom_cleavable && back.crosslink.asymmetric);
        assert_eq!(back.to_json(&[]), v);
    }

    #[test]
    fn missing_keys_take_sage_defaults() {
        let s = Settings::from_json(&json!({"database": {"fasta": "x.fasta"}})).unwrap();
        assert_eq!(s.missed_cleavages, 0);
        assert!(s.mods.is_empty());
        assert!(!s.crosslink.enabled);
        assert_eq!(s.precursor_charge, (2, 4));
    }

    /// Run the release `sage.exe` on configurations written by the GUI, with
    /// native (e.g. `C:\...`) paths. Skipped when sage has not been built.
    #[test]
    fn generated_configs_run_in_sage() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..").canonicalize().unwrap();
        let sage = root.join("target").join("release").join(if cfg!(windows) { "sage.exe" } else { "sage" });
        if !sage.is_file() {
            eprintln!("skipped: {} not built", sage.display());
            return;
        }
        let data = root.join("tests");
        let mut cases = vec![("linear", Settings::default())];
        let mut dsso = Settings::default();
        dsso.crosslink.enabled = true;
        dsso.crosslink.reagent = Reagent::Dsso;
        cases.push(("dsso", dsso.clone()));
        // edited preset stubs, wider charge range
        dsso.crosslink.stub_alt.clear();
        dsso.precursor_charge = (3, 6);
        cases.push(("dsso_edited_stubs", dsso));
        let mut reflect = Settings::default();
        reflect.crosslink.enabled = true;
        reflect.crosslink.backend = Backend::Reflect;
        cases.push(("reflect", reflect));
        let mut edc = Settings::default();
        edc.enzymes = vec!["trypsin".into(), "asp-n".into(), "glu-c".into()];
        edc.crosslink.enabled = true;
        edc.crosslink.reagent = Reagent::Edc;
        (edc.crosslink.site1, edc.crosslink.site2) = Reagent::Edc.sites().unwrap();
        edc.crosslink.two_sites = true;
        cases.push(("edc_codigest", edc));
        let mut multi = Settings::default();
        multi.enzymes = vec!["trypsin".into(), "lys-c".into(), "chymotrypsin".into(), "proalanase".into()];
        cases.push(("multi_enzyme", multi));
        let mut quant = Settings::default();
        quant.crosslink.enabled = true;
        quant.lfq = true;
        quant.tmt = Tmt::Tmt6;
        quant.tmt_level = 2;
        cases.push(("xl_quant", quant));
        let mut clixlink = Settings::default();
        clixlink.crosslink.enabled = true;
        clixlink.crosslink.reagent = Reagent::Clixlink;
        clixlink.crosslink.reset_stubs();
        clixlink.crosslink.stub2_alt.clear();
        cases.push(("clixlink_edited_stubs", clixlink));
        let mut with_linear = Settings::default();
        with_linear.crosslink.enabled = true;
        with_linear.crosslink.linear_search = true;
        cases.push(("xl_with_linear", with_linear));

        for (name, mut s) in cases {
            let out = std::env::temp_dir().join(format!("sage-gui-test-{name}"));
            let _ = std::fs::remove_dir_all(&out);
            s.fasta = data.join("Q99536.fasta").display().to_string();
            s.spectra = vec![data.join("LQSRPAAPPAPGPGQLTLR.mzML").display().to_string()];
            s.output_directory = out.display().to_string();
            let (errors, _) = s.check();
            assert!(errors.is_empty(), "{name}: {errors:?}");

            std::fs::create_dir_all(&out).unwrap();
            let config = out.join("sage_search.json");
            std::fs::write(&config, serde_json::to_string_pretty(&s.to_json(&s.spectra)).unwrap()).unwrap();
            let run = std::process::Command::new(&sage)
                .arg(&config)
                .args(["--threads", "2"])
                .arg("--disable-telemetry-i-dont-want-to-improve-sage")
                .output()
                .unwrap();
            assert!(run.status.success(), "{name}: {}", String::from_utf8_lossy(&run.stderr));
            let results = if s.crosslink.enabled { "crosslinks.sage.tsv" } else { "results.sage.tsv" };
            let tsv = std::fs::read_to_string(out.join(results)).unwrap();
            if name == "linear" || name == "xl_with_linear" {
                let linear = std::fs::read_to_string(out.join("results.sage.tsv")).unwrap();
                assert!(linear.contains("LQSRPAAPPAPGPGQLTLR"), "{name}: linear search should identify the peptide");
            }
            if name == "xl_quant" {
                // reporter columns on every crosslink match, and an LFQ table
                let header = tsv.lines().next().unwrap();
                assert!(header.ends_with("ion_injection_time\ttmt_1\ttmt_2\ttmt_3\ttmt_4\ttmt_5\ttmt_6"), "{header}");
                let lfq = std::fs::read_to_string(out.join("crosslinks.lfq.tsv")).unwrap();
                assert!(lfq.starts_with("peptide_a\tpeptide_b"));
            }
            let _ = std::fs::remove_dir_all(&out);
        }
    }
}
