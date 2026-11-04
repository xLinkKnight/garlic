//! Data behind the Spectra tab: the matches in a results folder, the spectra
//! they were found in (read exactly as Sage read them, by Sage's spectrum ids),
//! and their annotation by Annotator's engine (github.com/snijderlab/annotator):
//! `mzcore` parses the ProForma peptidoform, crosslinks included, `mzannotate`
//! generates its fragments and matches them to the peaks.

use mzannotate::fragment::FragmentType;
use mzannotate::prelude::*;
use mzcore::chemistry::OutputMolecularFormula;
use mzcore::prelude::*;
use mzcore::system::{isize::Charge, thomson, MassOverCharge};
use mzdata::mzpeaks::{CentroidPeak, PeakSet};
use mzdata::spectrum::{CentroidSpectrum, SignalContinuity, SpectrumDescription};
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::Url;
use sage_core::crosslink::CrosslinkSettings;
use sage_core::mass::Tolerance;
use sage_core::modification::is_proforma_accession;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Linear,
    Crosslink,
    Loop,
    /// Multivalent: two or more linkers (see [`Match::species`])
    Multi,
}

/// One match from `results.sage.tsv` or `crosslinks.sage.tsv`
#[derive(Clone, Debug)]
pub struct Match {
    pub kind: Kind,
    /// Spectrum file name and Sage's spectrum id
    pub file: String,
    pub scan: String,
    pub charge: u8,
    /// Sage's notation (linear peptide, or the crosslink's first peptide)
    pub peptide: String,
    /// Crosslink partner (empty for linear and loop-links)
    pub partner: String,
    /// 1-based link residues within the peptides (loop: both on `peptide`)
    pub sites: (usize, usize),
    pub proteins: String,
    pub score: f64,
    pub spectrum_q: f32,
    /// Peptide q-value (linear) or crosslink q-value
    pub q: f32,
    /// Target (crosslinks: target-target)
    pub target: bool,
    pub rank: u32,
    pub expmass: f32,
    pub calcmass: f32,
    pub rt: f32,
    /// Precursor inverse ion mobility 1/K0 (V·s/cm²), timsTOF data
    pub ion_mobility: Option<f32>,
    /// Crosslinks: both peptides from the same protein
    pub intra_protein: bool,
    /// Multivalent species: topology (`a2b`, `a1b1c`) and the `species`
    /// column (`[1] PEPK [2] KPEP | 1:K4-2:K1; …`); empty otherwise
    pub topology: String,
    pub species: String,
}

/// The peptides of a `species` column and its links as 0-based
/// `(peptide, residue, peptide, residue)`
pub fn parse_species(text: &str) -> (Vec<String>, Vec<(usize, usize, usize, usize)>) {
    let (peptides, links) = text.split_once(" | ").unwrap_or((text, ""));
    let is_number = |t: &str| t.len() > 2 && t.starts_with('[') && t.ends_with(']') && t[1..t.len() - 1].bytes().all(|b| b.is_ascii_digit());
    let peptides = peptides.split_whitespace().filter(|t| !is_number(t)).map(String::from).collect();
    // `1:K10` → (peptide 0, residue 9)
    let end = |s: &str| -> Option<(usize, usize)> {
        let (p, residue) = s.split_once(':')?;
        let digits: String = residue.chars().skip_while(|c| !c.is_ascii_digit()).collect();
        Some((p.trim().parse::<usize>().ok()? - 1, digits.parse::<usize>().ok()? - 1))
    };
    let links = links
        .split(';')
        .filter_map(|l| {
            let (a, b) = l.trim().split_once('-')?;
            let ((pa, sa), (pb, sb)) = (end(a)?, end(b)?);
            Some((pa, sa, pb, sb))
        })
        .collect();
    (peptides, links)
}

impl Match {
    /// Short description for lists: `PEPTIDE` or `A (K5) – B (K3)`
    pub fn title(&self) -> String {
        match self.kind {
            Kind::Linear => self.peptide.clone(),
            Kind::Crosslink => format!("{} ({}) – {} ({})", self.peptide, self.sites.0, self.partner, self.sites.1),
            Kind::Loop => format!("{} (loop {}–{})", self.peptide, self.sites.0, self.sites.1),
            Kind::Multi => format!("{} {}", self.topology, parse_species(&self.species).0.join(" – ")),
        }
    }

    pub fn key(&self) -> SpectrumKey {
        (self.file.clone(), self.scan.clone())
    }
}

/// (spectrum file name, Sage spectrum id)
pub type SpectrumKey = (String, String);

/// A results folder written by Sage
pub struct Results {
    pub matches: Vec<Match>,
    /// Spectrum file name → where Sage read it
    files: HashMap<String, Url>,
    bruker: BrukerProcessingConfig,
    /// Crosslink searches: the linker chemistry
    crosslink: Option<CrosslinkSettings>,
    /// Fragment tolerance of the search
    pub fragment_tol: Tolerance,
    /// Configured modification names → mass, for names ProForma readers do not know
    mod_masses: HashMap<String, f32>,
}

/// Header → column index of a TSV
fn columns(header: &str) -> HashMap<&str, usize> {
    header.split('\t').enumerate().map(|(i, h)| (h, i)).collect()
}

fn read_tsv(path: &Path) -> Option<String> {
    std::fs::read_to_string(path).ok()
}

/// The value in column `name` of a split TSV line ("" if absent)
fn field<'a>(line: &[&'a str], columns: &HashMap<&str, usize>, name: &str) -> &'a str {
    columns.get(name).and_then(|&i| line.get(i)).copied().unwrap_or("")
}

impl Results {
    /// Read `results.json` and the match tables of a results folder
    pub fn load(folder: &Path) -> Result<Results, String> {
        let json_path = folder.join("results.json");
        let text = std::fs::read_to_string(&json_path)
            .map_err(|e| format!("{}: {e} (run a search into this folder first)", json_path.display()))?;
        let json: Value = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", json_path.display()))?;

        let mut files = HashMap::new();
        for path in json.get("mzml_paths").and_then(Value::as_array).into_iter().flatten() {
            if let Some(url) = path.as_str().and_then(|p| sage_cloudpath::to_url(p).ok()) {
                let name = sage_cloudpath::filename(&url).unwrap_or(url.as_str()).to_string();
                files.insert(name, url);
            }
        }
        let bruker = json
            .get("bruker_config")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .unwrap_or_default();
        let crosslink: Option<CrosslinkSettings> = json
            .get("crosslink")
            .filter(|v| !v.is_null())
            .and_then(|v| serde_json::from_value(v.clone()).ok());
        let fragment_tol = json
            .get("fragment_tol")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .unwrap_or(Tolerance::Ppm(-20.0, 20.0));
        let mut mod_masses = HashMap::new();
        let database = json.get("database");
        for key in ["static_mods", "variable_mods"] {
            let Some(mods) = database.and_then(|d| d.get(key)).and_then(Value::as_object) else {
                continue;
            };
            for entry in mods.values().flat_map(|v| match v {
                Value::Array(list) => list.clone(),
                other => vec![other.clone()],
            }) {
                if let (Some(name), Some(mass)) = (
                    entry.get("name").and_then(Value::as_str),
                    entry.get("mass").and_then(Value::as_f64),
                ) {
                    mod_masses.insert(name.to_string(), mass as f32);
                }
            }
        }

        let mut matches = Vec::new();
        if let Some(text) = read_tsv(&folder.join("results.sage.tsv")) {
            matches.extend(parse_linear(&text));
        }
        if let Some(text) = read_tsv(&folder.join("crosslinks.sage.tsv")) {
            matches.extend(parse_crosslinks(&text));
        }
        if matches.is_empty() {
            return Err(format!("no results.sage.tsv or crosslinks.sage.tsv in {}", folder.display()));
        }
        Ok(Results {
            matches,
            files,
            bruker,
            crosslink,
            fragment_tol,
            mod_masses,
        })
    }

    /// The link of a crosslinked match: sites, linker and, for MS-cleavable
    /// linkers, each peptide with each stub on its linked residue (ProForma)
    pub fn link(&self, m: &Match) -> Option<Link> {
        let site = |s: usize| s.saturating_sub(1);
        let settings = self.crosslink.clone().unwrap_or_default();
        let sites = match m.kind {
            // multivalent species are non-cleavable: no stubs, peptides stacked in the map
            Kind::Linear | Kind::Multi => return None,
            Kind::Crosslink => vec![(0, site(m.sites.0)), (1, site(m.sites.1))],
            Kind::Loop => vec![(0, site(m.sites.0)), (0, site(m.sites.1))],
        };
        // A cleaved loop-link stays one peptide: its ring ions carry the stub of
        // one linked residue, so each residue's stubs are looked for in turn
        let stubs = stubs_of(&settings);
        let second = if m.kind == Kind::Loop { &m.peptide } else { &m.partner };
        let stub_peptides = [(&m.peptide, sites[0].1), (second, sites[1].1)]
            .iter()
            .filter(|_| !stubs.is_empty())
            .map(|(peptide, s)| {
                stubs
                    .iter()
                    .map(|stub| rewrite(peptide, &[(*s, format!("{:+}", stub.mass))], &self.mod_masses))
                    .collect()
            })
            .collect();
        Some(Link {
            kind: m.kind,
            sites,
            linker: linker_name(&settings),
            stubs,
            stub_peptides,
        })
    }

    /// ProForma 2.0 of a match, as Annotator reads it: crosslinks as
    /// `A[linker#XL1]..//B[#XL1]..`, modifications by accession or mass
    pub fn proforma(&self, m: &Match) -> String {
        let term = self
            .crosslink
            .as_ref()
            .map_or_else(|| "XLMOD:02001".to_string(), CrosslinkSettings::proforma_linker);
        let linker = format!("{term}#XL1");
        let site = |s: usize| s.saturating_sub(1);
        match m.kind {
            // every peptide, link n as #XLn on both its residues
            Kind::Multi => {
                let (peptides, links) = parse_species(&m.species);
                peptides
                    .iter()
                    .enumerate()
                    .map(|(x, peptide)| {
                        let marks: Vec<(usize, String)> = links
                            .iter()
                            .enumerate()
                            .flat_map(|(n, &(a, sa, b, sb))| {
                                [(a, sa, format!("{term}#XL{}", n + 1)), (b, sb, format!("#XL{}", n + 1))]
                            })
                            .filter(|&(p, _, _)| p == x)
                            .map(|(_, s, mark)| (s, mark))
                            .collect();
                        rewrite(peptide, &marks, &self.mod_masses)
                    })
                    .collect::<Vec<_>>()
                    .join("//")
            }
            Kind::Linear => rewrite(&m.peptide, &[], &self.mod_masses),
            Kind::Crosslink => format!(
                "{}//{}",
                rewrite(&m.peptide, &[(site(m.sites.0), linker)], &self.mod_masses),
                rewrite(&m.partner, &[(site(m.sites.1), "#XL1".into())], &self.mod_masses),
            ),
            Kind::Loop => rewrite(
                &m.peptide,
                &[(site(m.sites.0), linker), (site(m.sites.1), "#XL1".into())],
                &self.mod_masses,
            ),
        }
    }

    /// Read the spectra of `wanted` from their files, as Sage read them.
    /// `progress` gets a message per file.
    pub fn read_spectra(
        &self,
        wanted: &HashSet<SpectrumKey>,
        progress: &dyn Fn(String),
    ) -> Result<HashMap<SpectrumKey, Peaks>, String> {
        let mut by_file: HashMap<&str, HashSet<&str>> = HashMap::new();
        for (file, scan) in wanted {
            by_file.entry(file.as_str()).or_default().insert(scan.as_str());
        }
        let mut out = HashMap::new();
        for (file, ids) in by_file {
            let url = self
                .files
                .get(file)
                .ok_or_else(|| format!("{file} is not among the searched files in results.json"))?;
            progress(format!("reading {file} ({} spectra wanted)", ids.len()));
            let spectra = sage_cloudpath::util::read_spectra(url, 0, None, self.bruker, false)
                .map_err(|e| format!("{file}: {e}"))?;
            for s in spectra.into_iter().filter(|s| ids.contains(s.id.as_str())) {
                let precursor = s.precursors.first();
                let peaks = Peaks {
                    mz: s.mz.iter().map(|&x| x as f64).collect(),
                    intensity: s.intensity,
                    precursor_mz: precursor.map(|p| p.mz as f64),
                    ion_mobility: precursor.and_then(|p| p.inverse_ion_mobility),
                };
                out.insert((file.to_string(), s.id), peaks);
            }
        }
        Ok(out)
    }
}

/// A remnant an MS-cleavable linker leaves on its peptide
#[derive(Clone, Debug, PartialEq)]
pub struct Stub {
    /// Short name for labels: `A`, `T`, `S` (DSSO), `Bu`, `BuUr` (DSBU)
    pub name: String,
    pub description: String,
    pub mass: f64,
}

/// How a crosslinked match is linked
#[derive(Clone, Debug)]
pub struct Link {
    pub kind: Kind,
    /// (peptide, 0-based residue) of both ends; a loop-link has both on peptide 0
    pub sites: Vec<(usize, usize)>,
    /// Display name of the linker (`DSSO`)
    pub linker: String,
    /// Stubs of an MS-cleavable linker (both arms), empty otherwise
    pub stubs: Vec<Stub>,
    /// Peptide → stub → ProForma of the peptide alone with that stub on its
    /// linked residue
    pub stub_peptides: Vec<Vec<String>>,
}

fn linker_name(settings: &CrosslinkSettings) -> String {
    use sage_core::crosslink::Reagent;
    match (settings.reagent, settings.mass) {
        (Some(Reagent::Dss), _) | (None, None) => "DSS".into(),
        (Some(Reagent::Bs3), _) => "BS3".into(),
        (Some(Reagent::Dsg), _) => "DSG".into(),
        (Some(Reagent::Bs2g), _) => "BS2G".into(),
        (Some(Reagent::Dsso), _) => "DSSO".into(),
        (Some(Reagent::Dsbu), _) => "DSBU".into(),
        (Some(Reagent::Edc), _) => "EDC".into(),
        (Some(Reagent::Dmtmm), _) => "DMTMM".into(),
        (Some(Reagent::Clixlink), _) => "cliXlink".into(),
        (None, Some(mass)) => format!("{mass:+.4} Da"),
    }
}

/// The stubs of an MS-cleavable linker, named for labels. Asymmetric linkers
/// (cliXlink) list both arms; the second arm's names carry a prime.
pub fn stubs_of(settings: &CrosslinkSettings) -> Vec<Stub> {
    use sage_core::crosslink::Reagent;
    let Some(stubs) = settings.stubs() else { return Vec::new() };
    let stub = |name: &str, description: &str, mass: f32| Stub {
        name: name.into(),
        description: description.into(),
        mass: mass as f64,
    };
    // (short, long, alt...) names
    let arm = |short: f32, long: f32, alt: &[f32], prime: &str, out: &mut Vec<Stub>| match settings.reagent {
        Some(Reagent::Dsso | Reagent::Clixlink) => {
            out.push(stub(&format!("A{prime}"), "alkene", short));
            for (i, &m) in alt.iter().enumerate() {
                let name = if i == 0 { format!("T{prime}") } else { format!("T{}{prime}", i + 1) };
                out.push(stub(&name, "thiol", m));
            }
            out.push(stub(&format!("S{prime}"), "sulfenic acid", long));
        }
        Some(Reagent::Dsbu) => {
            out.push(stub(&format!("Bu{prime}"), "butyl amine", short));
            for (i, &m) in alt.iter().enumerate() {
                out.push(stub(&format!("x{}{prime}", i + 1), "stub", m));
            }
            out.push(stub(&format!("BuUr{prime}"), "butyl urea", long));
        }
        _ => {
            out.push(stub(&format!("s{prime}"), "short stub", short));
            for (i, &m) in alt.iter().enumerate() {
                out.push(stub(&format!("x{}{prime}", i + 1), "stub", m));
            }
            out.push(stub(&format!("l{prime}"), "long stub", long));
        }
    };
    let mut out = Vec::new();
    arm(stubs.short, stubs.long, &stubs.alt, "", &mut out);
    if let Some(other) = &stubs.other {
        arm(other.short, other.long, &other.alt, "'", &mut out);
    }
    for s in &mut out {
        if settings.reagent == Some(Reagent::Clixlink) && s.name.ends_with('\'') {
            s.description = format!("{} (alkyne arm)", s.description);
        }
    }
    out
}

/// Sage's peptide notation → ProForma: modification names that are not
/// ontology accessions become their mass; `extra` adds bracketed text after
/// residues (0-based), after the residue's own modification
pub fn rewrite(display: &str, extra: &[(usize, String)], mod_masses: &HashMap<String, f32>) -> String {
    let mod_text = |content: &str| -> String {
        if content.parse::<f64>().is_ok() || is_proforma_accession(content) {
            content.to_string()
        } else if let Some(mass) = mod_masses.get(content) {
            format!("{mass:+}")
        } else {
            content.to_string()
        }
    };
    let mut out = String::with_capacity(display.len() + 32);
    let mut residue: Option<usize> = None;
    let mut flushed = true;
    let flush = |out: &mut String, residue: Option<usize>| {
        if let Some(r) = residue {
            for (_, text) in extra.iter().filter(|(site, _)| *site == r) {
                out.push('[');
                out.push_str(text);
                out.push(']');
            }
        }
    };
    let mut chars = display.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '[' => {
                let mut depth = 1;
                let mut content = String::new();
                for c in chars.by_ref() {
                    match c {
                        '[' => depth += 1,
                        ']' => {
                            depth -= 1;
                            if depth == 0 {
                                break;
                            }
                        }
                        _ => {}
                    }
                    content.push(c);
                }
                out.push('[');
                out.push_str(&mod_text(&content));
                out.push(']');
            }
            c if c.is_ascii_alphabetic() => {
                if !flushed {
                    flush(&mut out, residue);
                }
                residue = Some(residue.map_or(0, |r| r + 1));
                flushed = false;
                out.push(c);
            }
            other => {
                // `-` before a C-terminal modification
                if !flushed {
                    flush(&mut out, residue);
                    flushed = true;
                }
                out.push(other);
            }
        }
    }
    if !flushed {
        flush(&mut out, residue);
    }
    out
}

fn parse_linear(text: &str) -> Vec<Match> {
    let mut lines = text.lines();
    let Some(header) = lines.next() else { return Vec::new() };
    let col = columns(header);
    lines
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .map(|f| Match {
            kind: Kind::Linear,
            file: field(&f, &col, "filename").to_string(),
            scan: field(&f, &col, "scannr").to_string(),
            charge: field(&f, &col, "charge").parse().unwrap_or(2),
            peptide: field(&f, &col, "peptide").to_string(),
            partner: String::new(),
            sites: (0, 0),
            proteins: field(&f, &col, "proteins").to_string(),
            score: field(&f, &col, "hyperscore").parse().unwrap_or(0.0),
            spectrum_q: field(&f, &col, "spectrum_q").parse().unwrap_or(1.0),
            q: field(&f, &col, "peptide_q").parse().unwrap_or(1.0),
            target: field(&f, &col, "label") == "1",
            rank: field(&f, &col, "rank").parse().unwrap_or(1),
            expmass: field(&f, &col, "expmass").parse().unwrap_or(0.0),
            calcmass: field(&f, &col, "calcmass").parse().unwrap_or(0.0),
            rt: field(&f, &col, "rt").parse().unwrap_or(0.0),
            ion_mobility: field(&f, &col, "ion_mobility").parse().ok().filter(|&x: &f32| x > 0.0),
            intra_protein: false,
            topology: String::new(),
            species: String::new(),
        })
        .collect()
}

fn parse_crosslinks(text: &str) -> Vec<Match> {
    let mut lines = text.lines();
    let Some(header) = lines.next() else { return Vec::new() };
    let col = columns(header);
    lines
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .map(|f| {
            let proteins = [field(&f, &col, "proteins_a"), field(&f, &col, "proteins_b")]
                .into_iter()
                .filter(|p| !p.is_empty())
                .collect::<Vec<_>>()
                .join(" – ");
            let kind = match field(&f, &col, "link_type") {
                "loop" => Kind::Loop,
                "multivalent" => Kind::Multi,
                _ => Kind::Crosslink,
            };
            Match {
                kind,
                file: field(&f, &col, "filename").to_string(),
                scan: field(&f, &col, "scannr").to_string(),
                charge: field(&f, &col, "charge").parse().unwrap_or(3),
                peptide: field(&f, &col, "peptide_a").to_string(),
                partner: field(&f, &col, "peptide_b").to_string(),
                sites: (
                    field(&f, &col, "link_a").parse().unwrap_or(1),
                    field(&f, &col, "link_b").parse().unwrap_or(1),
                ),
                proteins,
                score: field(&f, &col, "hyperscore").parse().unwrap_or(0.0),
                spectrum_q: field(&f, &col, "spectrum_q").parse().unwrap_or(1.0),
                q: field(&f, &col, "crosslink_q").parse().unwrap_or(1.0),
                target: field(&f, &col, "class") == "TT" || field(&f, &col, "label") == "1",
                rank: field(&f, &col, "rank").parse().unwrap_or(1),
                expmass: field(&f, &col, "expmass").parse().unwrap_or(0.0),
                calcmass: field(&f, &col, "calcmass").parse().unwrap_or(0.0),
                rt: field(&f, &col, "rt").parse().unwrap_or(0.0),
                ion_mobility: field(&f, &col, "ion_mobility").parse().ok().filter(|&x: &f32| x > 0.0),
                intra_protein: field(&f, &col, "intra_protein") == "true",
                topology: if kind == Kind::Multi { field(&f, &col, "topology").to_string() } else { String::new() },
                species: if kind == Kind::Multi { field(&f, &col, "species").to_string() } else { String::new() },
            }
        })
        .collect()
}

/// The peaks of one spectrum
#[derive(Clone, Debug, Default)]
pub struct Peaks {
    pub mz: Vec<f64>,
    pub intensity: Vec<f32>,
    pub precursor_mz: Option<f64>,
    /// Precursor inverse ion mobility 1/K0 (V·s/cm²), timsTOF data
    pub ion_mobility: Option<f32>,
}

/// Annotator's built-in fragmentation models
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    CidHcd,
    EThcd,
    Etd,
    Ead,
    All,
}

impl Model {
    pub const ALL: [Model; 5] = [Model::CidHcd, Model::EThcd, Model::Etd, Model::Ead, Model::All];

    pub fn label(self) -> &'static str {
        match self {
            Model::CidHcd => "CID / HCD",
            Model::EThcd => "EThcD",
            Model::Etd => "ETD",
            Model::Ead => "EAD",
            Model::All => "All ion types",
        }
    }

    fn model(self) -> &'static FragmentationModel {
        match self {
            Model::CidHcd => FragmentationModel::cid(),
            Model::EThcd => FragmentationModel::etcid(),
            Model::Etd => FragmentationModel::etd(),
            Model::Ead => FragmentationModel::ead(),
            Model::All => FragmentationModel::all(),
        }
    }
}

/// What kind of ion a peak was annotated as, for colours
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Series {
    A,
    B,
    C,
    X,
    Y,
    Z,
    Precursor,
    Immonium,
    Internal,
    Other,
    /// Carries an MS-cleavable linker stub (index into [`Link::stubs`])
    Stub(u8),
}

impl Series {
    fn of(ion: &FragmentType) -> Series {
        match ion {
            FragmentType::a(..) => Series::A,
            FragmentType::b(..) => Series::B,
            FragmentType::c(..) => Series::C,
            FragmentType::x(..) => Series::X,
            FragmentType::y(..) => Series::Y,
            FragmentType::z(..) => Series::Z,
            FragmentType::Precursor => Series::Precursor,
            FragmentType::Immonium(..) => Series::Immonium,
            FragmentType::Internal(..) => Series::Internal,
            _ => Series::Other,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Series::A => "a",
            Series::B => "b",
            Series::C => "c",
            Series::X => "x",
            Series::Y => "y",
            Series::Z => "z",
            Series::Precursor => "precursor",
            Series::Immonium => "immonium",
            Series::Internal => "internal",
            Series::Other => "other",
            Series::Stub(_) => "stub",
        }
    }
}

/// One interpretation of a peak
#[derive(Clone, Debug)]
pub struct IonLabel {
    /// e.g. `y7`, `b5²+`, `y4-H2O`, `y6+S` (with the S stub)
    pub text: String,
    pub series: Series,
    /// Which peptide of a crosslink (0 = first / α, 1 = second / β)
    pub peptide: usize,
    pub theoretical_mz: f64,
    /// Backbone ions: (N-terminal ion?, residue it ends after / starts at)
    backbone: Option<(bool, usize)>,
}

/// Where a stub-carrying ion lies on its peptide
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StubIonKind {
    /// The whole peptide with the stub (the doublet peaks)
    Peptide,
    /// N-terminal (a/b/c) ion ending after this residue, holding the linked one
    NTerm(usize),
    /// C-terminal (x/y/z) ion starting at this residue, holding the linked one
    CTerm(usize),
}

/// An ion with a linker stub found in the spectrum
#[derive(Clone, Debug)]
pub struct StubIon {
    pub peptide: usize,
    /// Index into [`Link::stubs`]
    pub stub: usize,
    pub kind: StubIonKind,
    pub charge: u8,
    pub theoretical_mz: f64,
    /// The peak it matched (index into [`Annotation::peaks`])
    pub peak: usize,
    /// e.g. `y6+S`, `P+A²+` (peptide + stub)
    pub text: String,
}

#[derive(Clone, Debug)]
pub struct PeakView {
    pub mz: f64,
    pub intensity: f32,
    pub labels: Vec<IonLabel>,
}

/// Backbone ions seen per residue boundary of each peptide
#[derive(Clone, Debug, Default)]
pub struct Coverage {
    pub residues: Vec<String>,
    /// N-terminal ion (a/b/c) ending after residue i
    pub n_term: Vec<bool>,
    /// C-terminal ion (x/y/z) starting at residue i
    pub c_term: Vec<bool>,
}

pub struct Annotation {
    pub proforma: String,
    pub peaks: Vec<PeakView>,
    pub coverage: Vec<Coverage>,
    pub explained_peaks: usize,
    pub explained_intensity: f32,
    pub max_intensity: f32,
    /// Crosslinked matches: the link, and the stub ions found
    pub link: Option<Link>,
    pub stub_ions: Vec<StubIon>,
}

/// `y7`, `b5²+`, `y4-H2O`: the label of a fragment
fn fragment_text(f: &Fragment<OutputMolecularFormula>) -> String {
    let (prefix, main) = f.ion.label();
    let mut text = format!("{}{}{}", prefix.unwrap_or_default(), main, f.ion.position_label().unwrap_or_default());
    // charge as `y9²+`: the GUI font has superscript digits but no ⁺
    let z = f.charge.value.unsigned_abs();
    if z > 1 {
        text.push_str(&superscript(z));
        text.push('+');
    }
    for loss in &f.neutral_loss {
        text.push_str(&tidy_formula(&loss.to_string()));
    }
    text
}

/// Is `observed` within the search's fragment tolerance of `theoretical`?
fn within(tolerance: Tolerance, theoretical: f64, observed: f64) -> bool {
    match tolerance {
        Tolerance::Ppm(lo, hi) => {
            ((observed - theoretical) / theoretical * 1e6).abs() <= lo.abs().max(hi.abs()) as f64
        }
        Tolerance::Da(lo, hi) | Tolerance::Pct(lo, hi) => {
            (observed - theoretical).abs() <= lo.abs().max(hi.abs()) as f64
        }
    }
}

/// The peak (index) nearest to `mz` within tolerance; `peaks` sorted by m/z
fn find_peak(peaks: &[PeakView], mz: f64, tolerance: Tolerance) -> Option<usize> {
    let i = peaks.partition_point(|p| p.mz < mz);
    [i.checked_sub(1), Some(i)]
        .into_iter()
        .flatten()
        .filter(|&j| j < peaks.len() && within(tolerance, mz, peaks[j].mz))
        .min_by(|&a, &b| (peaks[a].mz - mz).abs().total_cmp(&(peaks[b].mz - mz).abs()))
}

/// Find the ions that carry a linker stub: each peptide alone, with each stub
/// on its linked residue, fragmented by Annotator's engine. Kept: the whole
/// peptide and the backbone ions that hold the linked residue, without
/// neutral losses. Their peaks get the stub's label (replacing the engine's
/// own interpretation of the same mass).
fn find_stub_ions(views: &mut [PeakView], link: &Link, charge: u8, model: Model, tolerance: Tolerance) -> Vec<StubIon> {
    let mut found: Vec<StubIon> = Vec::new();
    let is_loop = link.kind == Kind::Loop;
    for (side, per_stub) in link.stub_peptides.iter().enumerate() {
        let (k, site) = link.sites.get(side).copied().unwrap_or((side, 0));
        // loop-links: ring ions contain this residue but not the other one
        let other = link.sites.get(1 - side.min(1)).map(|s| s.1).filter(|_| is_loop);
        for (j, proforma) in per_stub.iter().enumerate() {
            let Ok((ions, _)) = PeptidoformIonSet::pro_forma(proforma, &mzcore::ontology::STATIC_ONTOLOGIES) else {
                continue;
            };
            let n = ions
                .peptidoform_ions()
                .first()
                .and_then(|i| i.peptidoforms().first())
                .map_or(0, |p| p.len());
            // a stub peptide is charge-reduced: at least one charge went with the other half
            let max_charge = charge.saturating_sub(1).max(1);
            let fragments = ions.generate_theoretical_fragments::<OutputMolecularFormula>(
                Charge::new::<mzcore::system::e>(max_charge as isize),
                model.model(),
            );
            // The whole peptide + stub (the doublet peaks), at every charge up to
            // `max_charge`: the model only gives the precursor at its own charge
            let peptide_mass = ions.formulas().first().map(|f| f.monoisotopic_mass().value);
            // (not for loop-links: their peptide carries both stubs)
            for z in (1..=max_charge).filter(|_| !is_loop) {
                let Some(mass) = peptide_mass else { break };
                let mz = (mass + z as f64 * 1.007_276_466_6) / z as f64;
                let Some(peak) = find_peak(views, mz, tolerance) else { continue };
                let stub = &link.stubs[j];
                let text = if z > 1 {
                    format!("P+{}{}+", stub.name, superscript(z as usize))
                } else {
                    format!("P+{}", stub.name)
                };
                add_stub_label(views, peak, mz, k, j, &text);
                found.push(StubIon { peptide: k, stub: j, kind: StubIonKind::Peptide, charge: z, theoretical_mz: mz, peak, text });
            }
            for f in fragments.iter().filter(|f| f.neutral_loss.is_empty()) {
                // N-terminal ion `bN` holds residues 0..N, C-terminal `yN` n-N..n.
                // Loop-links: the ion must hold this residue and not the other,
                // so N-terminal ring ions belong to the first residue and
                // C-terminal ones to the second
                let n_ok = |series: usize| series > site && other.is_none_or(|o| o > site && series <= o);
                let c_ok = |series: usize| {
                    series >= n.saturating_sub(site) && series <= n && other.is_none_or(|o| o < site && n - series > o)
                };
                let kind = match &f.ion {
                    FragmentType::a(p, _) | FragmentType::b(p, _) | FragmentType::c(p, _) if n_ok(p.series_number) => {
                        StubIonKind::NTerm(p.series_number - 1)
                    }
                    FragmentType::x(p, _) | FragmentType::y(p, _) | FragmentType::z(p, _) if c_ok(p.series_number) => {
                        StubIonKind::CTerm(n - p.series_number)
                    }
                    _ => continue,
                };
                let Some(mz) = f.mz(MassMode::Monoisotopic).map(|m| m.value) else { continue };
                let Some(peak) = find_peak(views, mz, tolerance) else { continue };
                let z = f.charge.value.unsigned_abs() as u8;
                if found
                    .iter()
                    .any(|s| s.peptide == k && s.stub == j && s.kind == kind && s.charge == z)
                {
                    continue;
                }
                // `y6²+` → `y6+S²+`: the stub before the charge
                let base = fragment_text(f);
                let (ion, charge_text) = match base.find(|c: char| "⁰¹²³⁴⁵⁶⁷⁸⁹".contains(c)) {
                    Some(i) => (base[..i].to_string(), base[i..].to_string()),
                    None => (base.clone(), String::new()),
                };
                let text = format!("{ion}+{}{charge_text}", link.stubs[j].name);
                add_stub_label(views, peak, mz, k, j, &text);
                found.push(StubIon { peptide: k, stub: j, kind, charge: z, theoretical_mz: mz, peak, text });
            }
        }
    }
    found
}

/// Label a peak as a stub ion, replacing the engine's interpretations of the
/// same mass (which do not say which stub)
fn add_stub_label(views: &mut [PeakView], peak: usize, mz: f64, peptide: usize, stub: usize, text: &str) {
    let labels = &mut views[peak].labels;
    labels.retain(|l| (l.theoretical_mz - mz).abs() > 2e-3 || matches!(l.series, Series::Stub(_)));
    labels.insert(
        0,
        IonLabel {
            text: text.to_string(),
            series: Series::Stub(stub as u8),
            peptide,
            theoretical_mz: mz,
            backbone: None,
        },
    );
}

/// `-H2O1` → `-H2O`: drop element counts of one
fn tidy_formula(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::with_capacity(text.len());
    for (i, &c) in chars.iter().enumerate() {
        let after_element = i > 0 && chars[i - 1].is_ascii_alphabetic();
        let ends_number = chars.get(i + 1).is_none_or(|n| !n.is_ascii_digit());
        if c == '1' && after_element && ends_number {
            continue;
        }
        out.push(c);
    }
    out
}

fn superscript(n: usize) -> String {
    const DIGITS: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
    n.to_string().chars().map(|c| DIGITS[c.to_digit(10).unwrap_or(0) as usize]).collect()
}

/// Annotate `peaks` with the peptidoform `proforma` at precursor `charge`;
/// for a crosslink, `link` adds its stub ions (MS-cleavable linkers)
pub fn annotate(
    peaks: &Peaks,
    proforma: &str,
    charge: u8,
    model: Model,
    tolerance: Tolerance,
    link: Option<&Link>,
) -> Result<Annotation, String> {
    let ions = PeptidoformIonSet::pro_forma(proforma, &mzcore::ontology::STATIC_ONTOLOGIES)
        .map_err(|errors| {
            errors
                .iter()
                .map(|e| e.to_string().lines().filter(|l| !l.trim().is_empty()).last().unwrap_or("").to_string())
                .collect::<Vec<_>>()
                .join("; ")
        })?
        .0;
    let fragments = ions.generate_theoretical_fragments::<OutputMolecularFormula>(
        Charge::new::<mzcore::system::e>(charge.max(1) as isize),
        model.model(),
    );

    let mut description = SpectrumDescription::default();
    description.ms_level = 2;
    description.signal_continuity = SignalContinuity::Centroid;
    let centroids = peaks
        .mz
        .iter()
        .zip(&peaks.intensity)
        .enumerate()
        .map(|(i, (&mz, &intensity))| CentroidPeak::new(mz, intensity, i as u32))
        .collect::<Vec<_>>();
    let spectrum = CentroidSpectrum::new(description, PeakSet::new(centroids));
    let width = match tolerance {
        Tolerance::Ppm(lo, hi) => mzcore::quantities::Tolerance::new_ppm(lo.abs().max(hi.abs()) as f64),
        Tolerance::Da(lo, hi) | Tolerance::Pct(lo, hi) => mzcore::quantities::Tolerance::new_absolute(
            MassOverCharge::new::<thomson>(lo.abs().max(hi.abs()) as f64),
        ),
    };
    let parameters = MatchingParameters::default().tolerance(width);
    let annotated = spectrum.annotate(ions.clone(), &fragments, &parameters, MassMode::Monoisotopic);

    // Residues of each peptide, for the coverage map
    let mut coverage: Vec<Coverage> = ions
        .peptidoform_ions()
        .first()
        .map(|ion| {
            ion.peptidoforms()
                .iter()
                .map(|p| {
                    let residues = p.sequence().iter().map(|s| s.aminoacid.to_string()).collect::<Vec<_>>();
                    let n = residues.len();
                    Coverage { residues, n_term: vec![false; n], c_term: vec![false; n] }
                })
                .collect()
        })
        .unwrap_or_default();

    let mut views = Vec::new();
    for peak in annotated.peaks.iter() {
        let mut labels = Vec::new();
        for f in &peak.annotations {
            let peptide = f.peptidoform_index.unwrap_or(0);
            let length = coverage.get(peptide).map_or(0, |c| c.residues.len());
            let backbone = match &f.ion {
                FragmentType::a(p, _) | FragmentType::b(p, _) | FragmentType::c(p, _) => {
                    Some((true, p.series_number.saturating_sub(1)))
                }
                FragmentType::x(p, _) | FragmentType::y(p, _) | FragmentType::z(p, _) => {
                    length.checked_sub(p.series_number).map(|i| (false, i))
                }
                _ => None,
            };
            labels.push(IonLabel {
                text: fragment_text(f),
                series: Series::of(&f.ion),
                peptide,
                theoretical_mz: f.mz(MassMode::Monoisotopic).map_or(peak.mz.value, |mz| mz.value),
                backbone,
            });
        }
        views.push(PeakView { mz: peak.mz.value, intensity: peak.intensity, labels });
    }
    views.sort_by(|a, b| a.mz.total_cmp(&b.mz));

    let stub_ions = match link {
        Some(link) if !link.stubs.is_empty() => find_stub_ions(&mut views, link, charge, model, tolerance),
        _ => Vec::new(),
    };

    // Backbone coverage from the plain ions (stub ions are drawn on their own)
    let (mut explained_peaks, mut explained_intensity, mut total, mut max_intensity) = (0usize, 0.0f32, 0.0f32, 0.0f32);
    for peak in &views {
        total += peak.intensity;
        max_intensity = max_intensity.max(peak.intensity);
        if !peak.labels.is_empty() {
            explained_peaks += 1;
            explained_intensity += peak.intensity;
        }
        for l in &peak.labels {
            let (Some((n_term, residue)), Some(c)) = (l.backbone, coverage.get_mut(l.peptide)) else { continue };
            let slots = if n_term { &mut c.n_term } else { &mut c.c_term };
            if let Some(slot) = slots.get_mut(residue) {
                *slot = true;
            }
        }
    }
    Ok(Annotation {
        proforma: proforma.to_string(),
        peaks: views,
        coverage,
        explained_peaks,
        explained_intensity: if total > 0.0 { explained_intensity / total } else { 0.0 },
        max_intensity,
        link: link.cloned(),
        stub_ions,
    })
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;

    #[test]
    fn sage_notation_to_proforma() {
        let none = HashMap::new();
        assert_eq!(rewrite("PEPTIDE", &[], &none), "PEPTIDE");
        assert_eq!(rewrite("[+42.0106]-PEPTM[+15.9949]IDE", &[], &none), "[+42.0106]-PEPTM[+15.9949]IDE");
        // a link after the residue's own modification
        assert_eq!(
            rewrite("M[+15.9949]AGLEK", &[(0, "XLMOD:02001#XL1".into())], &none),
            "M[+15.9949][XLMOD:02001#XL1]AGLEK"
        );
        assert_eq!(rewrite("STKPLEK", &[(2, "#XL1".into()), (6, "X".into())], &none), "STK[#XL1]PLEK[X]");
        assert_eq!(rewrite("PEPK-[+0.984]", &[(3, "#XL1".into())], &none), "PEPK[#XL1]-[+0.984]");
        // free text names become masses, accessions stay
        let named = HashMap::from([("my ox".to_string(), 15.9949f32)]);
        assert_eq!(rewrite("PEPM[my ox]K", &[], &named), "PEPM[+15.9949]K");
        assert_eq!(rewrite("PEPM[UNIMOD:35]K", &[], &named), "PEPM[UNIMOD:35]K");
    }

    /// Peaks at the b/y ions of PEPTIDE (singly charged) and some noise
    fn peptide_peaks() -> Peaks {
        let ions = PeptidoformIonSet::pro_forma("PEPTIDE", &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0;
        let fragments = ions.generate_theoretical_fragments::<OutputMolecularFormula>(
            Charge::new::<mzcore::system::e>(1),
            FragmentationModel::cid(),
        );
        let mut mz: Vec<f64> = fragments
            .iter()
            .filter(|f| matches!(f.ion, FragmentType::b(..) | FragmentType::y(..)) && f.neutral_loss.is_empty())
            .filter_map(|f| f.mz(MassMode::Monoisotopic).map(|m| m.value))
            .collect();
        mz.extend([150.5, 333.3, 555.5]);
        mz.sort_by(f64::total_cmp);
        mz.dedup();
        Peaks { intensity: vec![100.0; mz.len()], mz, precursor_mz: None, ion_mobility: None }
    }

    #[test]
    fn annotates_linear_peptide() {
        let peaks = peptide_peaks();
        let a = annotate(&peaks, "PEPTIDE", 1, Model::CidHcd, Tolerance::Ppm(-10.0, 10.0), None).unwrap();
        assert_eq!(a.peaks.len(), peaks.mz.len());
        assert_eq!(a.explained_peaks, peaks.mz.len() - 3, "noise stays unexplained");
        assert!(a.peaks.iter().any(|p| p.labels.iter().any(|l| l.text == "y2" && l.series == Series::Y)));
        assert!(a.peaks.iter().any(|p| p.labels.iter().any(|l| l.text == "b3" && l.series == Series::B)));
        assert_eq!(a.coverage.len(), 1);
        assert_eq!(a.coverage[0].residues.join(""), "PEPTIDE");
        assert!(a.coverage[0].c_term[5], "y2 starts at residue 6 (D)");
    }

    #[test]
    fn crosslinks_parse_and_errors_are_readable() {
        let peaks = peptide_peaks();
        let a = annotate(
            &peaks,
            "GLEK[XLMOD:02126#XL1]VSTAAMQR//TDYHFPNEWK[#XL1]LIR",
            4,
            Model::All,
            Tolerance::Ppm(-20.0, 20.0),
            None,
        )
        .unwrap();
        assert_eq!(a.coverage.len(), 2, "both peptides");
        let error = annotate(&peaks, "PEPT[nonsense]IDE", 2, Model::CidHcd, Tolerance::Ppm(-20.0, 20.0), None).err().unwrap();
        assert!(!error.is_empty());
    }

    fn dsso_link() -> Link {
        let settings = CrosslinkSettings {
            reagent: Some(sage_core::crosslink::Reagent::Dsso),
            ..Default::default()
        };
        let stubs = stubs_of(&settings);
        let none = HashMap::new();
        let sites = vec![(0, 3), (1, 9)];
        let stub_peptides = [("GLEKVSTAAMQR", 3), ("TDYHFPNEWKLIR", 9)]
            .iter()
            .map(|(p, s)| stubs.iter().map(|st| rewrite(p, &[(*s, format!("{:+}", st.mass))], &none)).collect())
            .collect();
        Link { kind: Kind::Crosslink, sites, linker: "DSSO".into(), stubs, stub_peptides }
    }

    /// A multivalent chain from `crosslinks.sage.tsv`: read, written as
    /// ProForma with one `#XLn` per link, and annotated as three peptides
    #[test]
    fn multivalent_species_annotate() {
        let species = "[1] GLEKVSTAAMQR [2] AEFGHDSTWKLMNPQSTVKGYEDFAR [3] TDYHFPNEWKLIR | 1:K4-2:K10; 2:K19-3:K10";
        let (peptides, links) = parse_species(species);
        assert_eq!(peptides, ["GLEKVSTAAMQR", "AEFGHDSTWKLMNPQSTVKGYEDFAR", "TDYHFPNEWKLIR"]);
        assert_eq!(links, [(0, 3, 1, 9), (1, 18, 2, 9)]);
        // peptides with modifications keep their brackets apart from the [n] numbers
        assert_eq!(parse_species("[1] M[+15.9949]K [2] KR | 1:K2-2:K1").0, ["M[+15.9949]K", "KR"]);

        let results = Results {
            matches: Vec::new(),
            files: HashMap::new(),
            bruker: Default::default(),
            crosslink: Some(CrosslinkSettings::default()),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            mod_masses: HashMap::new(),
        };
        let m = Match {
            kind: Kind::Multi,
            file: "a.d".into(),
            scan: "1".into(),
            charge: 4,
            peptide: peptides[0].clone(),
            partner: peptides[1].clone(),
            sites: (4, 10),
            proteins: String::new(),
            score: 50.0,
            spectrum_q: 0.0,
            q: 0.0,
            target: true,
            rank: 1,
            expmass: 0.0,
            calcmass: 0.0,
            rt: 10.0,
            ion_mobility: Some(1.0),
            intra_protein: false,
            topology: "a1b1c".into(),
            species: species.into(),
        };
        assert_eq!(m.title(), "a1b1c GLEKVSTAAMQR – AEFGHDSTWKLMNPQSTVKGYEDFAR – TDYHFPNEWKLIR");
        let proforma = results.proforma(&m);
        assert_eq!(
            proforma,
            "GLEK[XLMOD:02001#XL1]VSTAAMQR//AEFGHDSTWK[#XL1]LMNPQSTVK[XLMOD:02001#XL2]GYEDFAR//TDYHFPNEWK[#XL2]LIR"
        );
        assert!(results.link(&m).is_none());
        let peaks = Peaks { mz: vec![300.0, 500.0, 700.0], intensity: vec![10.0; 3], precursor_mz: None, ion_mobility: None };
        let a = annotate(&peaks, &proforma, 4, Model::CidHcd, Tolerance::Ppm(-10.0, 10.0), None).unwrap();
        assert_eq!(a.coverage.len(), 3);
    }

    /// An opened DSSO loop-link (NEKIQQKLQNER, K3–K7): ring ions carry the stub
    /// of one residue: b3–b6 the first, y6–y9 the second
    #[test]
    fn loop_link_ring_ions_carry_one_stub() {
        let settings = CrosslinkSettings { reagent: Some(sage_core::crosslink::Reagent::Dsso), ..Default::default() };
        let stubs = stubs_of(&settings);
        let none = HashMap::new();
        let (peptide, first, second) = ("NEKIQQKLQNER", 2usize, 6usize);
        let stub_peptides = [first, second]
            .iter()
            .map(|&s| stubs.iter().map(|st| rewrite(peptide, &[(s, format!("{:+}", st.mass))], &none)).collect())
            .collect();
        let link = Link { kind: Kind::Loop, sites: vec![(0, first), (0, second)], linker: "DSSO".into(), stubs, stub_peptides };

        // b4 with the alkene on K3 and y7 with the sulfenic acid on K7 (ring ions),
        // and y9 with the alkene on K7 (it holds K3 too: not a ring ion of K7 alone)
        let fragments = |proforma: &str| {
            PeptidoformIonSet::pro_forma(proforma, &mzcore::ontology::STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .generate_theoretical_fragments::<OutputMolecularFormula>(Charge::new::<mzcore::system::e>(1), FragmentationModel::cid())
        };
        let pick = |proforma: &str, n_term: bool, number: usize| {
            fragments(proforma)
                .into_iter()
                .find(|f| {
                    f.neutral_loss.is_empty()
                        && f.charge.value == 1
                        && match &f.ion {
                            FragmentType::b(p, _) => n_term && p.series_number == number,
                            FragmentType::y(p, _) => !n_term && p.series_number == number,
                            _ => false,
                        }
                })
                .and_then(|f| f.mz(MassMode::Monoisotopic))
                .unwrap()
                .value
        };
        let mut mz = vec![
            pick(&link.stub_peptides[0][0], true, 4),
            pick(&link.stub_peptides[1][2], false, 7),
            pick(&link.stub_peptides[1][0], false, 11),
        ];
        mz.sort_by(f64::total_cmp);
        let peaks = Peaks { intensity: vec![50.0; mz.len()], mz, precursor_mz: None, ion_mobility: None };
        let a = annotate(&peaks, "NEK[XLMOD:02126#XL1]IQQK[#XL1]LQNER", 3, Model::CidHcd, Tolerance::Ppm(-10.0, 10.0), Some(&link))
            .unwrap();
        let mut found: Vec<(StubIonKind, usize)> = a.stub_ions.iter().map(|s| (s.kind, s.stub)).collect();
        found.sort_by_key(|(k, s)| (format!("{k:?}"), *s));
        assert_eq!(found, vec![(StubIonKind::CTerm(5), 2), (StubIonKind::NTerm(3), 0)], "{:?}", a.stub_ions);
        assert!(a.stub_ions.iter().all(|s| s.peptide == 0));
    }

    #[test]
    fn dsso_stubs_are_named() {
        let names: Vec<_> = dsso_link().stubs.iter().map(|s| (s.name.clone(), (s.mass * 1000.0).round() as i64)).collect();
        assert_eq!(names, [("A".to_string(), 54011), ("T".to_string(), 85983), ("S".to_string(), 103993)]);
        let clix = stubs_of(&CrosslinkSettings {
            reagent: Some(sage_core::crosslink::Reagent::Clixlink),
            ..Default::default()
        });
        assert_eq!(clix.len(), 6);
        assert!(clix[3].name == "A'" && clix[3].description.contains("alkyne"));
        assert!(stubs_of(&CrosslinkSettings::default()).is_empty(), "DSS is not cleavable");
    }

    #[test]
    fn finds_stub_ions_of_a_cleaved_crosslink() {
        let (link, a) = cleaved_dsso();
        let has = |peptide: usize, stub: &str, kind: StubIonKind| {
            a.stub_ions
                .iter()
                .any(|s| s.peptide == peptide && link.stubs[s.stub].name == stub && s.kind == kind)
        };
        for peptide in [0, 1] {
            assert!(has(peptide, "A", StubIonKind::Peptide), "peptide {peptide} + A: {:?}", a.stub_ions);
            assert!(has(peptide, "S", StubIonKind::Peptide), "peptide {peptide} + S");
        }
        // y4..y6 of TDYHFPNEWKLIR start at residues 9, 8, 7 (0-based)
        for start in [9, 8, 7] {
            assert!(has(1, "S", StubIonKind::CTerm(start)), "y from residue {start}: {:?}", a.stub_ions);
        }
        assert!(a.stub_ions.iter().all(|s| s.charge == 1));
        // every peak is a stub ion and labelled as such
        assert_eq!(a.explained_peaks, a.peaks.len());
        assert!(a.peaks.iter().all(|p| matches!(p.labels[0].series, Series::Stub(_))));
        assert!(a.peaks.iter().any(|p| p.labels[0].text == "y4+S"));
        assert!(a.peaks.iter().any(|p| p.labels[0].text == "P+A"));
        // stub ions are not counted as plain backbone coverage
        assert!(a.coverage[1].c_term.iter().all(|c| !c));
    }

    /// A DSSO crosslink spectrum after linker cleavage: each peptide with the
    /// A and S stubs (singly charged), and y4-y6+S of the second peptide
    pub(crate) fn cleaved_dsso() -> (Link, Annotation) {
        let link = dsso_link();
        let proton = 1.007276;
        // The doublets: each peptide with the A and S stubs, singly charged
        let mass = |p: &str| {
            PeptidoformIonSet::pro_forma(p, &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0.formulas()[0]
                .monoisotopic_mass()
                .value
        };
        let mut mz = Vec::new();
        for (p, stub) in [("GLEKVSTAAMQR", 0), ("GLEKVSTAAMQR", 2), ("TDYHFPNEWKLIR", 0), ("TDYHFPNEWKLIR", 2)] {
            mz.push(mass(p) + link.stubs[stub].mass + proton);
        }
        // and y ions of the second peptide with the S stub (y4 = KLIR holds K10)
        let s_peptide = &link.stub_peptides[1][2];
        let ions = PeptidoformIonSet::pro_forma(s_peptide, &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0;
        for f in ions.generate_theoretical_fragments::<OutputMolecularFormula>(Charge::new::<mzcore::system::e>(1), FragmentationModel::cid()) {
            if let FragmentType::y(p, _) = &f.ion {
                if f.neutral_loss.is_empty() && (4..=6).contains(&p.series_number) && f.charge.value == 1 {
                    mz.push(f.mz(MassMode::Monoisotopic).unwrap().value);
                }
            }
        }
        mz.sort_by(f64::total_cmp);
        let peaks = Peaks { intensity: vec![50.0; mz.len()], mz, precursor_mz: None, ion_mobility: None };
        let a = annotate(
            &peaks,
            "GLEK[XLMOD:02126#XL1]VSTAAMQR//TDYHFPNEWK[#XL1]LIR",
            4,
            Model::CidHcd,
            Tolerance::Ppm(-10.0, 10.0),
            Some(&link),
        )
        .unwrap();
        (link, a)
    }
}

/// End to end on a real results folder: `SAGE_GUI_RESULTS=<folder> cargo test
/// -p sage-gui real_results -- --ignored --nocapture`
#[cfg(test)]
mod real_data {
    use super::*;

    #[test]
    #[ignore]
    fn real_results() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let results = Results::load(Path::new(&folder)).unwrap();
        let pick = |kind: Kind| {
            results
                .matches
                .iter()
                .filter(move |m| m.kind == kind && m.target && m.rank == 1 && m.spectrum_q <= 0.01)
                .take(40)
                .collect::<Vec<_>>()
        };
        let chosen = [pick(Kind::Linear), pick(Kind::Crosslink), pick(Kind::Loop)];
        let wanted = chosen.iter().flatten().map(|m| m.key()).collect::<HashSet<_>>();
        let start = std::time::Instant::now();
        let spectra = results.read_spectra(&wanted, &|msg| println!("{msg}")).unwrap();
        println!("read {} spectra in {:.1} s", spectra.len(), start.elapsed().as_secs_f32());
        for (name, list) in ["linear", "crosslink", "loop-link"].iter().zip(&chosen) {
            for model in [Model::CidHcd, Model::All] {
                let (mut right, mut wrong, mut n) = (0.0, 0.0, 0);
                let keys: Vec<_> = list.iter().map(|m| m.key()).collect();
                for (i, m) in list.iter().enumerate() {
                    let proforma = results.proforma(m);
                    let annotate_at = |key: &SpectrumKey| {
                        let a = annotate(&spectra[key], &proforma, m.charge, model, results.fragment_tol, results.link(m).as_ref())
                            .unwrap_or_else(|e| panic!("{proforma}: {e}"));
                        a.explained_intensity
                    };
                    right += annotate_at(&keys[i]);
                    // control: the same peptidoform on another match's spectrum
                    wrong += annotate_at(&keys[(i + 1) % keys.len()]);
                    n += 1;
                    if i < 2 {
                        println!("  {proforma}");
                    }
                }
                println!(
                    "{name:9} {:28} intensity explained: own spectrum {:5.1}%, other spectrum {:5.1}% ({n} matches)",
                    model.label(),
                    100.0 * right / n as f32,
                    100.0 * wrong / n as f32
                );
            }
        }

        // Stub ions of the crosslinks, against the same peptidoforms on other spectra
        let list = &chosen[1];
        let keys: Vec<_> = list.iter().map(|m| m.key()).collect();
        let stats = |offset: usize| {
            let (mut both_doublets, mut fragments, mut matches) = (0, 0, 0);
            for (i, m) in list.iter().enumerate() {
                let Some(link) = results.link(m) else { continue };
                let key = &keys[(i + offset) % keys.len()];
                let a = annotate(&spectra[key], &results.proforma(m), m.charge, Model::CidHcd, results.fragment_tol, Some(&link)).unwrap();
                let on = |k: usize| a.stub_ions.iter().any(|s| s.peptide == k && s.kind == StubIonKind::Peptide);
                both_doublets += (on(0) && on(1)) as usize;
                fragments += a.stub_ions.iter().filter(|s| s.kind != StubIonKind::Peptide).count();
                matches += 1;
                if offset == 0 && i < 3 {
                    println!("  {}: {:?}", results.proforma(m), a.stub_ions.iter().map(|s| (s.peptide + 1, s.text.clone())).collect::<Vec<_>>());
                }
            }
            (both_doublets, fragments, matches)
        };
        let (own, other) = (stats(0), stats(1));
        println!(
            "crosslink stub ions: own spectrum: peptide+stub on both peptides in {}/{}, {} stub fragment ions; \
             other spectrum: {}/{}, {} stub fragment ions",
            own.0, own.2, own.1, other.0, other.2, other.1
        );

        // Ring ions of the loop-links (stub on one residue), against other spectra
        let list = &chosen[2];
        if !list.is_empty() {
            let keys: Vec<_> = list.iter().map(|m| m.key()).collect();
            let ring = |offset: usize| {
                let (mut with, mut ions) = (0, 0);
                for (i, m) in list.iter().enumerate() {
                    let Some(link) = results.link(m) else { continue };
                    let key = &keys[(i + offset) % keys.len()];
                    let a = annotate(&spectra[key], &results.proforma(m), m.charge, Model::CidHcd, results.fragment_tol, Some(&link))
                        .unwrap();
                    with += (!a.stub_ions.is_empty()) as usize;
                    ions += a.stub_ions.len();
                    if offset == 0 && i < 2 {
                        println!("  {}: {:?}", results.proforma(m), a.stub_ions.iter().map(|s| s.text.clone()).collect::<Vec<_>>());
                    }
                }
                (with, ions)
            };
            let (own, other) = (ring(0), ring(1));
            println!(
                "loop-link ring ions: own spectrum: in {}/{} spectra, {} ions; other spectrum: in {}/{}, {} ions",
                own.0,
                list.len(),
                own.1,
                other.0,
                list.len(),
                other.1
            );
        }

        // Ion mobility of the timsTOF precursors, and example figures when
        // SAGE_GUI_FIGURES names a folder
        let with_mobility = spectra.values().filter(|p| p.ion_mobility.is_some()).count();
        println!("{with_mobility} of {} spectra carry a precursor ion mobility", spectra.len());
        if let Ok(dir) = std::env::var("SAGE_GUI_FIGURES") {
            use crate::viewer::figure::{export, svg, FigureOptions, Format};
            for (name, m) in [("crosslink", chosen[1][1]), ("linear", chosen[0][0])] {
                let peaks = &spectra[&m.key()];
                let a = annotate(peaks, &results.proforma(m), m.charge, Model::CidHcd, results.fragment_tol, results.link(m).as_ref()).unwrap();
                let title = format!(
                    "{} · {}+ · scan {} · RT {:.2} min · 1/K0 {:.3} V·s/cm²",
                    m.title(), m.charge, m.scan, m.rt, peaks.ion_mobility.unwrap_or(0.0)
                );
                for (background, suffix) in [(true, "white"), (false, "transparent")] {
                    let opts = FigureOptions { background, title: Some(title.clone()), label_all: false, mz_range: None, coverage: true, dpi: 300 };
                    let text = svg(&a, &opts);
                    for format in [Format::Svg, Format::Png, Format::Tiff] {
                        let path = Path::new(&dir).join(format!("{name}_{suffix}.{}", format.extension()));
                        export(&text, format, 300, background, &path).unwrap();
                    }
                }
                println!("figures for {} written", m.title());
            }
        }
    }
}

#[cfg(test)]
mod labels {
    #[test]
    fn neutral_losses_read_naturally() {
        assert_eq!(super::tidy_formula("-H2O1"), "-H2O");
        assert_eq!(super::tidy_formula("-H3N1"), "-H3N");
        assert_eq!(super::tidy_formula("-C1H4O1S1"), "-CH4OS");
        assert_eq!(super::tidy_formula("-C10H12"), "-C10H12");
    }
}
