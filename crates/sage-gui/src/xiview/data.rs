//! An XNova results folder as the PRIDE crosslinking API serves an uploaded
//! mzIdentML file to xiVIEW: the same JSON rows (`get_xiview_matches`,
//! `get_xiview_peptides`, `get_xiview_proteins` and the search metadata),
//! built from `results.json` and `crosslinks.sage.tsv` through the viewer's
//! [`Results`].
//!
//! * Crosslinks are donor / acceptor peptide pairs (`MS:1002509` /
//!   `MS:1002510` on the link residues, the linker mass on the donor), as in
//!   XNova's `crosslinks.mzid`; loop-links one peptide with both sites.
//! * A multivalent species (xiVIEW knows pairs only) is one match per link.
//! * Decoy proteins are `DECOY_<accession>`, a prefix xiVIEW maps back to
//!   the target protein.

use crate::viewer::model::{parse_species, Kind, Match, Results};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::path::Path;

const PROTON: f64 = 1.007_276_466_6;
/// XNova's crosslink FDR at which matches are served to xiVIEW: target
/// matches (rank 1) with crosslink q-value ≤ this. xiVIEW shows them as
/// passing ("Stored validation"); its own decoy-based FDR is not used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fdr(pub f32);

impl Default for Fdr {
    /// 1%, as `crosslinks_FDRfiltered.sage.tsv`
    fn default() -> Self {
        Fdr(0.01)
    }
}

impl Fdr {
    pub const LEVELS: [Fdr; 3] = [Fdr(0.01), Fdr(0.05), Fdr(0.10)];

    /// `1%`
    pub fn percent(self) -> String {
        format!("{}%", (self.0 * 100.0).round())
    }
}

/// A served match: the XNova match and, for a multivalent species, how
/// xiVIEW's spectrum viewer shows the whole species
#[derive(Clone, Debug)]
pub struct Served {
    pub index: usize,
    pub layout: Option<Layout>,
    /// Its peptides as xiSPEC sends them back (residues, each followed by
    /// its modification names), to tell an unchanged request: as xiVIEW
    /// asks (the two peptides of the link), and for a multivalent species
    /// as its spectrum viewer shows it (every peptide)
    pub signatures: Vec<Vec<String>>,
}

/// A multivalent species in xiVIEW's spectrum viewer (xiSPEC): every
/// peptide as a row, the served link between the first two
#[derive(Clone, Debug, Default)]
pub struct Layout {
    /// xiSPEC row → peptide of the species (its index in `species`)
    pub order: Vec<usize>,
    /// xiSPEC's `Peptides` and `LinkSite` (one `id` per link)
    pub peptides: Value,
    pub link_sites: Value,
    /// Its modifications as xiSPEC names them: `(oxidation)` → mass
    pub modifications: Vec<(String, f64)>,
}

/// One results folder, served as one xiVIEW upload
pub struct Upload {
    pub id: u32,
    /// The results folder and the FDR its matches are served at
    pub folder: std::path::PathBuf,
    pub fdr: Fdr,
    /// The folder name: xiVIEW's `file` URL parameter
    pub name: String,
    pub results: Results,
    /// Spectrum file names; the index is xiVIEW's spectra data id
    pub files: Vec<String>,
    /// xiVIEW match id → XNova match
    pub served: HashMap<String, Served>,
    /// Endpoint name (`get_xiview_matches`) → JSON
    endpoints: HashMap<&'static str, Vec<u8>>,
}

impl Upload {
    pub fn endpoint(&self, name: &str) -> Option<&[u8]> {
        self.endpoints.get(name).map(Vec::as_slice)
    }

    /// The XNova match behind a xiVIEW match id
    pub fn served(&self, id: &str) -> Option<(&Match, &Served)> {
        let s = self.served.get(id)?;
        Some((self.results.matches.get(s.index)?, s))
    }

    pub fn load(folder: &Path, id: u32, fdr: Fdr) -> Result<Upload, String> {
        let results = Results::load(folder)?;
        let json: Value = std::fs::read_to_string(folder.join("results.json"))
            .ok()
            .and_then(|t| serde_json::from_str(&t).ok())
            .unwrap_or(Value::Null);
        let name = folder.file_name().map_or_else(|| folder.display().to_string(), |n| n.to_string_lossy().into_owned());
        let files: Vec<(String, String)> = results.spectrum_files();
        let mut builder = Builder::new(id, &results, &json);
        let file_index: HashMap<&str, usize> = files.iter().enumerate().map(|(i, (n, _))| (n.as_str(), i)).collect();

        let mut matches = Vec::new();
        let mut served = HashMap::new();
        // multivalent links: (XNova match, link, donor and acceptor peptide rows)
        let mut multivalent: Vec<(usize, usize, usize, usize)> = Vec::new();
        for (index, m) in results.matches.iter().enumerate() {
            if m.kind == Kind::Linear || m.rank != 1 {
                continue;
            }
            // XNova's FDR: targets within it
            if !(m.target && m.q <= fdr.0) {
                continue;
            }
            let sd = file_index.get(m.file.as_str()).copied().unwrap_or(0);
            match m.kind {
                Kind::Crosslink => {
                    let a = builder.peptide(&m.peptide, &m.protein_lists[0], Role::Donor(m.sites.0), m.target);
                    let b = builder.peptide(&m.partner, m.protein_lists.get(1).map_or("", String::as_str), Role::Acceptor(m.sites.1), m.target);
                    let mid = format!("{id}_{index}");
                    matches.push(match_row(&mid, m, id, a, Some(b), sd));
                    let signature = vec![builder.signature(&m.peptide), builder.signature(&m.partner)];
                    served.insert(mid, Served { index, layout: None, signatures: vec![signature] });
                }
                Kind::Loop => {
                    let a = builder.peptide(&m.peptide, &m.protein_lists[0], Role::Loop(m.sites.0, m.sites.1), m.target);
                    let mid = format!("{id}_{index}");
                    let mut row = match_row(&mid, m, id, a, None, sd);
                    row["pi2"] = Value::Null;
                    matches.push(row);
                    let signature = vec![builder.signature(&m.peptide)];
                    served.insert(mid, Served { index, layout: None, signatures: vec![signature] });
                }
                Kind::Multi => {
                    // one match per link: xiVIEW shows pairs
                    let (peptides, links) = parse_species(&m.species);
                    for (n, &(pa, sa, pb, sb)) in links.iter().enumerate() {
                        let (Some(x), Some(y)) = (peptides.get(pa), peptides.get(pb)) else { continue };
                        let px = builder.located(x, m.target);
                        let py = builder.located(y, m.target);
                        let a = builder.peptide(x, &px, Role::Donor(sa + 1), m.target);
                        let b = builder.peptide(y, &py, Role::Acceptor(sb + 1), m.target);
                        let mid = format!("{id}_{index}_XL{}", n + 1);
                        matches.push(match_row(&mid, m, id, a, Some(b), sd));
                        let layout = builder.layout(&peptides, &links, n);
                        let species: Vec<String> = layout.order.iter().map(|&p| builder.signature(&peptides[p])).collect();
                        let signatures = vec![vec![builder.signature(x), builder.signature(y)], species];
                        served.insert(mid.clone(), Served { index, layout: Some(layout), signatures });
                        multivalent.push((index, n, a, b));
                    }
                }
                Kind::Linear => {}
            }
        }

        let mut endpoints: HashMap<&'static str, Vec<u8>> = HashMap::new();
        let mut put = |name: &'static str, v: Value| {
            endpoints.insert(name, serde_json::to_vec(&v).unwrap_or_else(|_| b"[]".to_vec()));
        };
        put("get_xiview_matches", Value::Array(matches));
        put("get_xiview_peptides", Value::Array(builder.peptides.clone()));
        put("get_xiview_proteins", Value::Array(builder.protein_rows()));
        let meta = Metadata { id, name: &name, json: &json, results: &results, files: &files, fdr };
        put("get_xiview_mzidentml_files", meta.mzidentml_files());
        put("get_xiview_analysis_collection_spectrum_identifications", meta.analysis_collection());
        put("get_xiview_spectrum_identification_protocols", meta.protocols());
        put("get_xiview_spectra_data", meta.spectra_data());
        put("get_xiview_enzymes", meta.enzymes());
        put("get_xiview_search_modifications", meta.search_modifications());
        // the links of multivalent species, as link metadata xiVIEW loads (CSV)
        let metadata = builder.link_metadata(&results, &multivalent);
        endpoints.insert("get_xnova_link_metadata", metadata.into_bytes());

        Ok(Upload {
            id,
            folder: folder.to_path_buf(),
            fdr,
            name,
            files: files.into_iter().map(|(n, _)| n).collect(),
            results,
            served,
            endpoints,
        })
    }
}

fn match_row(id: &str, m: &Match, upload: u32, a: usize, b: Option<usize>, sd: usize) -> Value {
    let z = m.charge.max(1) as f64;
    json!({
        "id": id,
        "pi1": a,
        "pi2": b,
        "sc": {"XNova:hyperscore": m.score, "XNova:q-value": m.q},
        "ui": upload,
        "c_mz": (m.calcmass as f64 + z * PROTON) / z,
        "pc_c": m.charge,
        "pc_mz": (m.expmass as f64 + z * PROTON) / z,
        "sp": m.scan,
        "sd": sd,
        // every match served passes XNova's FDR (as the PRIDE API serves
        // only passing ones)
        "p": true,
        "r": m.rank,
        "sip": 0,
        "msi_id": Value::Null,
        "msi_pc": Value::Null,
    })
}

/// How a peptide is linked (1-based residues)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Role {
    /// Carries the linker mass (`MS:1002509`)
    Donor(usize),
    /// The other end (`MS:1002510`)
    Acceptor(usize),
    /// Both ends on one peptide
    Loop(usize, usize),
}

/// A peptide in Sage's notation: its residues and modifications as
/// (mzIdentML location: 0 N-terminus, i residue i (1-based), n + 1
/// C-terminus, the bracketed text)
pub fn parse_peptide(display: &str) -> (String, Vec<(usize, String)>) {
    let mut residues = String::new();
    let mut mods = Vec::new();
    let mut chars = display.chars().peekable();
    let mut after_c_dash = false;
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
                let location = if residues.is_empty() {
                    0
                } else if after_c_dash {
                    residues.len() + 1
                } else {
                    residues.len()
                };
                mods.push((location, content));
            }
            '-' => after_c_dash = !residues.is_empty(),
            c if c.is_ascii_alphabetic() => residues.push(c),
            _ => {}
        }
    }
    (residues, mods)
}

/// A modification name as xiVIEW writes it in peptides (`(oxidation)`)
pub fn xmod(name: &str) -> String {
    format!("({})", name.to_lowercase().split_whitespace().collect::<String>())
}

/// A request peptide (xiSPEC's residues and modifications) as a signature
pub fn request_signature(peptide: &Value) -> String {
    peptide["sequence"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|r| format!("{}{}", r["aminoAcid"].as_str().unwrap_or(""), r["Modification"].as_str().unwrap_or("")))
        .collect()
}

/// A decoy peptide back to its target: Sage reverses all but the ends
fn unreverse(sequence: &str) -> String {
    let mut r: Vec<char> = sequence.chars().collect();
    let n = r.len();
    if n > 2 {
        r[1..n - 1].reverse();
    }
    r.into_iter().collect()
}

/// UniProt accession of a FASTA accession (`sp|P04050|RPB1_YEAST` → `P04050`)
fn uniprot_accession(id: &str) -> &str {
    let parts: Vec<&str> = id.split('|').collect();
    if parts.len() >= 3 && matches!(parts[0], "sp" | "tr") {
        parts[1]
    } else {
        id
    }
}

/// Name and mass of common UNIMOD modifications Sage writes by accession
fn unimod_name(accession: &str) -> Option<&'static str> {
    Some(match accession {
        "UNIMOD:1" => "Acetyl",
        "UNIMOD:4" => "Carbamidomethyl",
        "UNIMOD:5" => "Carbamyl",
        "UNIMOD:7" => "Deamidated",
        "UNIMOD:21" => "Phospho",
        "UNIMOD:27" => "Glu->pyro-Glu",
        "UNIMOD:28" => "Gln->pyro-Glu",
        "UNIMOD:34" => "Methyl",
        "UNIMOD:35" => "Oxidation",
        "UNIMOD:36" => "Dimethyl",
        "UNIMOD:121" => "GG",
        "UNIMOD:188" => "Label:13C(6)",
        "UNIMOD:259" => "Label:13C(6)15N(2)",
        "UNIMOD:267" => "Label:13C(6)15N(4)",
        "UNIMOD:481" => "Label:2H(4)",
        "UNIMOD:737" => "TMT6plex",
        "UNIMOD:2016" => "TMTpro",
        _ => return None,
    })
}

/// Mass of a modification written as `content` in Sage's notation
fn modification_mass(content: &str, named: &HashMap<String, f64>, cache: &mut HashMap<String, f64>) -> f64 {
    if let Ok(mass) = content.parse::<f64>() {
        return mass;
    }
    if let Some(&mass) = named.get(content) {
        return mass;
    }
    if let Some(&mass) = cache.get(content) {
        return mass;
    }
    // an ontology accession: what it adds to a glycine
    use mzcore::prelude::*;
    let mass_of = |p: &str| {
        PeptidoformIonSet::pro_forma(p, &mzcore::ontology::STATIC_ONTOLOGIES)
            .ok()
            .and_then(|(set, _)| set.formulas().first().map(|f| f.monoisotopic_mass().value))
    };
    let mass = match (mass_of(&format!("G[{content}]")), mass_of("G")) {
        (Some(a), Some(b)) => a - b,
        _ => 0.0,
    };
    cache.insert(content.to_string(), mass);
    mass
}

/// Builds the peptide and protein rows, each peptide once
struct Builder<'a> {
    upload: u32,
    results: &'a Results,
    decoy_tag: String,
    /// Linker: ontology term (XLMOD) and name, mass
    linker_term: Option<(String, String)>,
    linker_mass: f64,
    /// Configured modification names → mass
    named: HashMap<String, f64>,
    mass_cache: HashMap<String, f64>,
    peptides: Vec<Value>,
    peptide_ids: HashMap<String, usize>,
    /// xiVIEW protein id → (FASTA accession, decoy)
    proteins: HashMap<String, (String, bool)>,
}

impl<'a> Builder<'a> {
    fn new(upload: u32, results: &'a Results, json: &Value) -> Builder<'a> {
        let settings = results.crosslink_settings().cloned().unwrap_or_default();
        let term = settings.proforma_linker();
        let linker_name = results
            .matches
            .iter()
            .find(|m| m.kind != Kind::Linear)
            .and_then(|m| results.link(m))
            .map_or_else(|| "crosslinker".to_string(), |l| l.linker);
        let mut named = HashMap::new();
        for key in ["static_mods", "variable_mods"] {
            for entry in json["database"][key].as_object().into_iter().flat_map(|o| o.values()).flat_map(|v| match v {
                Value::Array(list) => list.clone(),
                other => vec![other.clone()],
            }) {
                if let (Some(name), Some(mass)) = (entry["name"].as_str(), entry["mass"].as_f64()) {
                    named.insert(name.to_string(), mass);
                }
            }
        }
        Builder {
            upload,
            results,
            decoy_tag: results.decoy_tag().to_string(),
            linker_term: term.starts_with("XLMOD:").then(|| (term.clone(), linker_name)),
            linker_mass: settings.linker_mass() as f64,
            named,
            mass_cache: HashMap::new(),
            peptides: Vec::new(),
            peptide_ids: HashMap::new(),
            proteins: HashMap::new(),
        }
    }

    /// The name xiVIEW shows for a modification written `content`
    fn modification_name(&mut self, content: &str) -> (String, f64) {
        let mass = modification_mass(content, &self.named, &mut self.mass_cache);
        let name = if content.parse::<f64>().is_ok() {
            format!("{mass:+.4}")
        } else {
            unimod_name(content).map_or_else(|| content.to_string(), String::from)
        };
        (name, mass)
    }

    /// A peptide as xiSPEC sends it back: each residue followed by its
    /// modifications as `(name)` (terminal ones on the end residues)
    fn signature(&mut self, display: &str) -> String {
        let (bare, mods) = parse_peptide(display);
        let n = bare.len();
        let mut out = String::new();
        for (i, aa) in bare.chars().enumerate() {
            out.push(aa);
            for (location, content) in &mods {
                let here = *location == i + 1 || (*location == 0 && i == 0) || (*location == n + 1 && i + 1 == n);
                if here {
                    let (name, _) = self.modification_name(content);
                    out.push_str(&xmod(&name));
                }
            }
        }
        out
    }

    /// A peptide as xiSPEC's `Peptides` holds it (residues and modification
    /// names) and its modifications' masses
    fn xispec_peptide(&mut self, display: &str) -> (Value, Vec<(String, f64)>) {
        let (bare, mods) = parse_peptide(display);
        let n = bare.len();
        let mut residues = Vec::new();
        let mut masses = Vec::new();
        for (i, aa) in bare.chars().enumerate() {
            let mut names = String::new();
            for (location, content) in &mods {
                if *location == i + 1 || (*location == 0 && i == 0) || (*location == n + 1 && i + 1 == n) {
                    let (name, mass) = self.modification_name(content);
                    names.push_str(&xmod(&name));
                    masses.push((xmod(&name), mass));
                }
            }
            residues.push(json!({"aminoAcid": aa.to_string(), "Modification": names}));
        }
        (json!({"sequence": residues}), masses)
    }

    /// A multivalent species in xiSPEC, its link `shown` between the first
    /// two rows and the other peptides in the order the links reach them
    fn layout(&mut self, peptides: &[String], links: &[(usize, usize, usize, usize)], shown: usize) -> Layout {
        let (a, _, b, _) = links[shown];
        let mut order = vec![a];
        if b != a {
            order.push(b);
        }
        let mut grew = true;
        while grew {
            grew = false;
            for &(x, _, y, _) in links {
                for (from, to) in [(x, y), (y, x)] {
                    if order.contains(&from) && !order.contains(&to) {
                        order.push(to);
                        grew = true;
                    }
                }
            }
        }
        for p in 0..peptides.len() {
            if !order.contains(&p) {
                order.push(p);
            }
        }
        let row = |p: usize| order.iter().position(|&q| q == p).unwrap_or(0);
        // the shown link first (id 0)
        let mut sites = Vec::new();
        let mut ids = vec![shown];
        ids.extend((0..links.len()).filter(|&n| n != shown));
        for (id, &n) in ids.iter().enumerate() {
            let (x, sx, y, sy) = links[n];
            sites.push(json!({"id": id, "peptideId": row(x), "linkSite": sx}));
            sites.push(json!({"id": id, "peptideId": row(y), "linkSite": sy}));
        }
        let mut rows = Vec::new();
        let mut modifications: Vec<(String, f64)> = Vec::new();
        for &p in &order {
            let (peptide, masses) = self.xispec_peptide(&peptides[p]);
            rows.push(peptide);
            for (name, mass) in masses {
                if !modifications.iter().any(|(n, _)| *n == name) {
                    modifications.push((name, mass));
                }
            }
        }
        Layout { order, peptides: Value::Array(rows), link_sites: Value::Array(sites), modifications }
    }

    /// Link metadata for xiVIEW (its "link metadata" CSV): each link of a
    /// multivalent species gets its topology (`Multivalent`) and the species
    /// (`Multivalent species`: its links, the same for every link of it)
    fn link_metadata(&self, results: &Results, multivalent: &[(usize, usize, usize, usize)]) -> String {
        // (protein, residue) ends of a served peptide pair
        let ends = |a: usize, b: usize| -> Vec<((String, usize), (String, usize))> {
            let (pa, pb) = (&self.peptides[a], &self.peptides[b]);
            let site = |p: &Value, i: usize| -> Option<(String, usize)> {
                Some((p["prt"][i].as_str()?.to_string(), p["pos"][i].as_u64()? as usize + p["ls1"].as_u64()? as usize - 1))
            };
            let mut out = Vec::new();
            for i in 0..pa["prt"].as_array().map_or(0, Vec::len) {
                for j in 0..pb["prt"].as_array().map_or(0, Vec::len) {
                    if let (Some(x), Some(y)) = (site(pa, i), site(pb, j)) {
                        out.push((x, y));
                    }
                }
            }
            out
        };
        let label = |(protein, residue): &(String, usize)| {
            let accession = protein.strip_prefix("DECOY_").unwrap_or(protein);
            let p = results.protein(accession);
            let entry = accession.split('|').nth(2).unwrap_or(accession);
            let name = p.map(|p| p.gene.as_str()).filter(|g| !g.is_empty()).unwrap_or(entry);
            let aa = p.and_then(|p| p.sequence.chars().nth(residue.saturating_sub(1))).unwrap_or('X');
            format!("{name} {aa}{residue}")
        };
        // each species: its links, written once
        let mut species: HashMap<usize, Vec<String>> = HashMap::new();
        for &(index, _, a, b) in multivalent {
            if let Some((x, y)) = ends(a, b).first() {
                species.entry(index).or_default().push(format!("{}–{}", label(x), label(y)));
            }
        }
        // link (both orders) → (topologies, species labels)
        let mut links: std::collections::BTreeMap<(String, usize, String, usize), (Vec<String>, Vec<String>)> = Default::default();
        for &(index, _, a, b) in multivalent {
            let topology = results.matches[index].topology.clone();
            let mut name = species.get(&index).cloned().unwrap_or_default();
            name.sort();
            let name = format!("{topology}: {}", name.join(" + "));
            for (x, y) in ends(a, b) {
                let entry = links.entry((x.0, x.1, y.0, y.1)).or_default();
                if !entry.0.contains(&topology) {
                    entry.0.push(topology.clone());
                }
                if !entry.1.contains(&name) {
                    entry.1.push(name.clone());
                }
            }
        }
        if links.is_empty() {
            return String::new();
        }
        let quote = |s: &str| format!("\"{}\"", s.replace('"', "\"\""));
        let mut csv = String::from("Protein1,SeqPos1,Protein2,SeqPos2,Multivalent,Multivalent species\n");
        for ((p1, r1, p2, r2), (topologies, names)) in links {
            csv.push_str(&format!(
                "{},{r1},{},{r2},{},{}\n",
                quote(&p1),
                quote(&p2),
                quote(&topologies.join(" | ")),
                quote(&names.join(" | "))
            ));
        }
        csv
    }

    /// The proteins of a multivalent species' peptide (not listed by Sage):
    /// found in the FASTA, as a `;`-list like Sage's
    fn located(&self, display: &str, target: bool) -> String {
        let (bare, _) = parse_peptide(display);
        let mut found: Vec<String> = self
            .results
            .proteins()
            .filter(|(_, p)| p.sequence.contains(&bare))
            .map(|(a, _)| a.clone())
            .collect();
        if found.is_empty() && !target {
            let t = unreverse(&bare);
            found = self
                .results
                .proteins()
                .filter(|(_, p)| p.sequence.contains(&t))
                .map(|(a, _)| format!("{}{a}", self.decoy_tag))
                .collect();
        }
        found.sort();
        found.join(";")
    }

    /// The row of a peptide (once per sequence, modifications, link and proteins)
    fn peptide(&mut self, display: &str, proteins: &str, role: Role, _target: bool) -> usize {
        let key = format!("{display}|{role:?}|{proteins}");
        if let Some(&id) = self.peptide_ids.get(&key) {
            return id;
        }
        let (bare, mods) = parse_peptide(display);
        let (mut m_as, mut m_ps, mut m_ms) = (Vec::new(), Vec::new(), Vec::new());
        for (location, content) in &mods {
            let (name, mass) = self.modification_name(content);
            let accession = if content.contains(':') && content.parse::<f64>().is_err() {
                content.clone()
            } else {
                "MS:1001460".to_string()
            };
            m_as.push(json!({ accession: name }));
            m_ps.push(*location);
            m_ms.push(mass);
        }
        let donor = |pair: usize| {
            let mut cv = serde_json::Map::new();
            if let Some((term, name)) = &self.linker_term {
                cv.insert(term.clone(), json!(name));
            }
            cv.insert("MS:1002509".into(), json!(pair));
            Value::Object(cv)
        };
        let (ls1, ls2, cl_m) = match role {
            Role::Donor(s) => {
                m_as.push(donor(0));
                m_ps.push(s);
                m_ms.push(self.linker_mass);
                (s, None, self.linker_mass)
            }
            Role::Acceptor(s) => {
                m_as.push(json!({"MS:1002510": 0}));
                m_ps.push(s);
                m_ms.push(0.0);
                (s, None, 0.0)
            }
            Role::Loop(s1, s2) => {
                m_as.push(donor(0));
                m_ps.push(s1);
                m_ms.push(self.linker_mass);
                m_as.push(json!({"MS:1002510": 0}));
                m_ps.push(s2);
                m_ms.push(0.0);
                (s1, Some(s2), self.linker_mass)
            }
        };
        // mzIdentML orders modifications by location
        let mut order: Vec<usize> = (0..m_ps.len()).collect();
        order.sort_by_key(|&i| m_ps[i]);
        let m_as: Vec<Value> = order.iter().map(|&i| m_as[i].clone()).collect();
        let m_ps: Vec<usize> = order.iter().map(|&i| m_ps[i]).collect();
        let m_ms: Vec<f64> = order.iter().map(|&i| m_ms[i]).collect();

        let (mut prt, mut pos, mut dec) = (Vec::new(), Vec::new(), Vec::new());
        for accession in proteins.split(';').map(str::trim).filter(|a| !a.is_empty()) {
            let (target, decoy) = match accession.strip_prefix(&self.decoy_tag) {
                Some(t) if !self.decoy_tag.is_empty() => (t, true),
                _ => (accession, false),
            };
            let sequence = if decoy { unreverse(&bare) } else { bare.clone() };
            let start = self
                .results
                .protein(target)
                .and_then(|p| p.sequence.find(&sequence).or_else(|| p.sequence.replace('I', "L").find(&sequence.replace('I', "L"))))
                .map_or(1, |i| i + 1);
            let id = if decoy { format!("DECOY_{target}") } else { target.to_string() };
            self.proteins.insert(id.clone(), (target.to_string(), decoy));
            prt.push(id);
            pos.push(start);
            dec.push(decoy);
        }
        let id = self.peptides.len();
        self.peptides.push(json!({
            "id": id,
            "u_id": self.upload,
            "seq": bare,
            "prt": prt,
            "pos": pos,
            "dec": dec,
            "ls1": ls1,
            "ls2": ls2,
            "m_as": m_as,
            "m_ps": m_ps,
            "m_ms": m_ms,
            "cl_m": cl_m,
        }));
        self.peptide_ids.insert(key, id);
        id
    }

    fn protein_rows(&self) -> Vec<Value> {
        let mut ids: Vec<&String> = self.proteins.keys().collect();
        ids.sort();
        ids.into_iter()
            .map(|id| {
                let (accession, decoy) = &self.proteins[id];
                let protein = self.results.protein(accession);
                let entry = accession.split('|').nth(2).unwrap_or(accession);
                let name = protein.map(|p| p.gene.as_str()).filter(|g| !g.is_empty()).unwrap_or(entry);
                let prefix = if *decoy { "DECOY_" } else { "" };
                json!({
                    "id": id,
                    "name": format!("{prefix}{name}"),
                    "accession": format!("{prefix}{}", uniprot_accession(accession)),
                    "sequence": protein.map(|p| p.sequence.clone()).unwrap_or_default(),
                    "search_id": self.upload.to_string(),
                    "description": protein.map(|p| p.name.clone()).unwrap_or_default(),
                })
            })
            .collect()
    }
}

/// The search metadata rows
struct Metadata<'a> {
    id: u32,
    name: &'a str,
    json: &'a Value,
    results: &'a Results,
    files: &'a [(String, String)],
    fdr: Fdr,
}

impl Metadata<'_> {
    fn software(&self) -> Value {
        json!({
            "id": "xnova",
            "name": "XNova (Sage)",
            "version": self.json["version"].as_str().unwrap_or(env!("CARGO_PKG_VERSION")),
            "uri": "https://github.com/xLinkKnight/XNova",
        })
    }

    fn file_format(name: &str) -> &'static str {
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".d") || lower.ends_with(".d/") {
            "MS:1002817"
        } else if lower.contains(".mgf") {
            "MS:1001062"
        } else {
            "MS:1000584"
        }
    }

    fn mzidentml_files(&self) -> Value {
        json!([{
            "id": self.id,
            "project_id": "xnova",
            "identification_file_name": self.name,
            "analysis_software_list": {"AnalysisSoftware": [self.software()]},
            "provider": {},
            "audit_collection": {},
            "analysis_sample_collection": {},
            "bib": [{
                "id": "10.1021/acs.jproteome.3c00486",
                "title": "Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale",
                "authors": "Michael R. Lazear",
                "publication": "Journal of Proteome Research",
                "doi": "10.1021/acs.jproteome.3c00486",
            }],
            "spectra_formats": self.files.iter().enumerate().map(|(i, (name, location))| json!({
                "id": i,
                "location": location,
                "name": name,
                "FileFormat": {Self::file_format(name): ""},
                "SpectrumIDFormat": {"MS:1000774": "multiple peak list nativeID format"},
            })).collect::<Vec<_>>(),
        }])
    }

    fn analysis_collection(&self) -> Value {
        json!([{
            "upload_id": self.id,
            "spectrum_identification_list_ref": "SIL_1",
            "spectrum_identification_protocol_ref": "SIP_1",
            "protocol_id": 0,
            "spectra_data_refs": (0..self.files.len()).collect::<Vec<_>>(),
            "search_database_refs": ["SDB_1"],
        }])
    }

    fn protocols(&self) -> Value {
        let (tol, unit) = match self.results.fragment_tol {
            sage_core::mass::Tolerance::Ppm(lo, hi) => (lo.abs().max(hi.abs()) as f64, "ppm"),
            sage_core::mass::Tolerance::Da(lo, hi) | sage_core::mass::Tolerance::Pct(lo, hi) => (lo.abs().max(hi.abs()) as f64, "Da"),
        };
        json!([{
            "id": 0,
            "sip_ref": "SIP_1",
            "upload_id": self.id,
            "frag_tol": tol,
            "frag_tol_unit": unit,
            "additional_search_params": {
                "MS:1002494": "cross-linking search",
                "MS:1001211": "parent mass type mono",
                "MS:1001256": "fragment mass type mono",
            },
            "analysis_software": self.software(),
            "threshold": {"XNova:crosslink q-value": self.fdr.0},
        }])
    }

    fn spectra_data(&self) -> Value {
        Value::Array(
            self.files
                .iter()
                .enumerate()
                .map(|(i, (name, location))| {
                    json!({
                        "id": i,
                        "upload_id": self.id,
                        "location": location,
                        "name": name,
                        "external_format_documentation": Value::Null,
                        "file_format": Self::file_format(name),
                        "spectrum_id_format": "MS:1000774",
                    })
                })
                .collect(),
        )
    }

    fn enzymes(&self) -> Value {
        let e = &self.json["database"]["enzyme"];
        let missed = e["missed_cleavages"].as_u64().unwrap_or(2);
        let semi = e["semi_enzymatic"].as_bool().unwrap_or(false);
        let mut names: Vec<String> = e["enzymes"].as_array().into_iter().flatten().filter_map(Value::as_str).map(String::from).collect();
        if names.is_empty() {
            names.push(match (e["cleave_at"].as_str(), e["restrict"].as_str(), e["c_terminal"].as_bool()) {
                (Some(at), restrict, c) => format!("custom:{at}:{}:{}", restrict.unwrap_or(""), c.unwrap_or(true)),
                _ => "trypsin".into(),
            });
        }
        let rows: Vec<Value> = names
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let (regexp, accession, display) = enzyme_regexp(name);
                json!({
                    "id": format!("ENZ_{i}"),
                    "upload_id": self.id,
                    "protocol_id": 0,
                    "c_term_gain": "OH",
                    "min_distance": Value::Null,
                    "missed_cleavages": missed,
                    "n_term_gain": "H",
                    "name": display,
                    "semi_specific": semi,
                    "site_regexp": regexp,
                    "accession": accession,
                })
            })
            .collect();
        Value::Array(rows)
    }

    fn search_modifications(&self) -> Value {
        let mut rows = Vec::new();
        let db = &self.json["database"];
        for (key, fixed) in [("static_mods", true), ("variable_mods", false)] {
            let Some(mods) = db[key].as_object() else { continue };
            for (site, entries) in mods {
                let list = match entries {
                    Value::Array(list) => list.clone(),
                    other => vec![other.clone()],
                };
                for entry in list {
                    let mass = entry["mass"].as_f64().or_else(|| entry.as_f64()).unwrap_or(0.0);
                    let name = entry["name"].as_str();
                    let residues: String = site.chars().filter(|c| c.is_ascii_uppercase()).collect();
                    let mut accessions = serde_json::Map::new();
                    match name {
                        Some(n) if n.contains(':') => {
                            accessions.insert(n.to_string(), json!(unimod_name(n).unwrap_or(n)));
                        }
                        Some(n) => {
                            accessions.insert("MS:1001460".into(), json!(n));
                        }
                        None => {
                            accessions.insert("MS:1001460".into(), json!(format!("{mass:+.4}")));
                        }
                    }
                    if site.starts_with('^') || site.starts_with('[') {
                        accessions.insert("MS:1001189".into(), json!("modification specificity peptide N-term"));
                    }
                    if site.starts_with('$') || site.starts_with(']') {
                        accessions.insert("MS:1001190".into(), json!("modification specificity peptide C-term"));
                    }
                    rows.push(json!({
                        "id": rows.len(),
                        "upload_id": self.id,
                        "protocol_id": 0,
                        "mass": mass,
                        "residues": if residues.is_empty() { ".".to_string() } else { residues },
                        "fixed_mod": fixed,
                        "accessions": accessions,
                        "crosslinker_id": Value::Null,
                    }));
                }
            }
        }
        // the linker's two ends
        if let Some(settings) = self.results.crosslink_settings() {
            let linker = self
                .results
                .matches
                .iter()
                .find(|m| m.kind != Kind::Linear)
                .and_then(|m| self.results.link(m))
                .map_or_else(|| "crosslinker".to_string(), |l| l.linker);
            let term = settings.proforma_linker();
            let (first, second) = settings
                .hetero_sites()
                .unwrap_or_else(|| (settings.residues.clone(), settings.residues.clone()));
            for (donor, sites) in [(true, first), (false, second)] {
                let residues: String = sites.chars().filter(|c| c.is_ascii_uppercase()).collect();
                let mut accessions = serde_json::Map::new();
                if donor {
                    if term.starts_with("XLMOD:") {
                        accessions.insert(term.clone(), json!(linker));
                    }
                    accessions.insert("MS:1002509".into(), json!(linker));
                } else {
                    accessions.insert("MS:1002510".into(), json!("cross-link receiver"));
                }
                rows.push(json!({
                    "id": rows.len(),
                    "upload_id": self.id,
                    "protocol_id": 0,
                    "mass": if donor { settings.linker_mass() as f64 } else { 0.0 },
                    "residues": if residues.is_empty() { ".".to_string() } else { residues },
                    "fixed_mod": false,
                    "accessions": accessions,
                    "crosslinker_id": linker,
                }));
            }
        }
        Value::Array(rows)
    }
}

/// Cleavage rule of an enzyme as xiVIEW reads it (`(?<=[KR])(?!P)`), its
/// PSI-MS term and name
fn enzyme_regexp(name: &str) -> (String, &'static str, String) {
    let (r, acc, display): (&str, &'static str, &str) = match name.to_ascii_lowercase().as_str() {
        "trypsin" => ("(?<=[KR])(?!P)", "MS:1001251", "Trypsin"),
        "trypsin/p" => ("(?<=[KR])", "MS:1001313", "Trypsin/P"),
        "lys-c" => ("(?<=K)", "MS:1001309", "Lys-C"),
        "lys-n" => ("(?=K)", "MS:1003093", "Lys-N"),
        "arg-c" => ("(?<=R)(?!P)", "MS:1001303", "Arg-C"),
        "asp-n" => ("(?=D)", "MS:1001304", "Asp-N"),
        "glu-c" => ("(?<=E)(?!P)", "MS:1001917", "Glu-C"),
        "chymotrypsin" => ("(?<=[FWYL])(?!P)", "MS:1001306", "Chymotrypsin"),
        "chymotrypsin-high" => ("(?<=[FWY])(?!P)", "MS:1001306", "Chymotrypsin (high specificity)"),
        "proalanase" => ("(?<=[PA])", "MS:1001955", "ProAlanase"),
        other => {
            // custom:<sites>:<restrict>:<c-terminal>
            let parts: Vec<&str> = other.split(':').collect();
            if let [_, at, restrict, c] = parts.as_slice() {
                let at = at.to_ascii_uppercase();
                let restrict = restrict.to_ascii_uppercase();
                let rule = if *c == "false" {
                    format!("(?=[{at}])")
                } else if restrict.is_empty() {
                    format!("(?<=[{at}])")
                } else {
                    format!("(?<=[{at}])(?![{restrict}])")
                };
                return (rule, "MS:1001957", format!("custom ({at})"));
            }
            return ("(?<=[KR])(?!P)".into(), "MS:1001251", name.to_string());
        }
    };
    (r.to_string(), acc, display.to_string())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn sage_notation_to_mzidentml_locations() {
        assert_eq!(parse_peptide("PEPTIDE"), ("PEPTIDE".into(), vec![]));
        assert_eq!(parse_peptide("PEPM[+15.9949]K"), ("PEPMK".into(), vec![(4, "+15.9949".into())]));
        assert_eq!(
            parse_peptide("[+42.0106]-MPEK[UNIMOD:259]-[+0.984]"),
            ("MPEK".into(), vec![(0, "+42.0106".into()), (4, "UNIMOD:259".into()), (5, "+0.984".into())])
        );
        assert_eq!(unreverse("TDKPEK"), "TEPKDK");
        assert_eq!(uniprot_accession("sp|P04050|RPB1_YEAST"), "P04050");
        assert_eq!(uniprot_accession("my_protein"), "my_protein");
    }

    #[test]
    fn modification_masses() {
        let named = HashMap::from([("my ox".to_string(), 15.9949)]);
        let mut cache = HashMap::new();
        assert_eq!(modification_mass("+57.021465", &named, &mut cache), 57.021465);
        assert_eq!(modification_mass("my ox", &named, &mut cache), 15.9949);
        assert!((modification_mass("UNIMOD:35", &named, &mut cache) - 15.994915).abs() < 1e-4);
        assert!((modification_mass("UNIMOD:259", &named, &mut cache) - 8.014199).abs() < 1e-4);
    }

    #[test]
    fn enzyme_rules() {
        assert_eq!(enzyme_regexp("trypsin").0, "(?<=[KR])(?!P)");
        assert_eq!(enzyme_regexp("asp-n").0, "(?=D)");
        assert_eq!(enzyme_regexp("custom:DE:P:true").0, "(?<=[DE])(?![P])");
    }

    /// A small crosslink results folder: rows for xiVIEW
    #[test]
    fn serves_a_results_folder() {
        let dir = std::env::temp_dir().join(format!("xnova-xiview-data-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("p.fasta");
        std::fs::write(&fasta, ">sp|P11111|ONE_HUMAN One protein OS=Homo sapiens GN=ONE\nMKGLEKVSTAAMQRK\n>sp|P22222|TWO_HUMAN Two\nMTDYHFPNEWKLIRGG\n").unwrap();
        let config = json!({
            "version": "0.23.0",
            "database": {
                "fasta": fasta.display().to_string(), "decoy_tag": "rev_",
                "enzyme": {"missed_cleavages": 2, "enzymes": ["trypsin"]},
                "static_mods": {"C": {"mass": 57.021465, "name": "UNIMOD:4"}},
                "variable_mods": {"M": [{"mass": 15.994915, "max_count": 2, "name": "UNIMOD:35"}]}
            },
            "crosslink": {"reagent": "dss", "residues": "K"},
            "fragment_tol": {"ppm": [-20.0, 20.0]},
            "mzml_paths": [dir.join("run.mzML").display().to_string()]
        });
        std::fs::write(dir.join("results.json"), config.to_string()).unwrap();
        // (spectrum files are listed when they exist)
        std::fs::write(dir.join("run.mzML"), "").unwrap();
        let header = "psm_id\tfilename\tscannr\trank\tlabel\tclass\tlink_type\tpeptide_a\tpeptide_b\tlink_a\tlink_b\tproteins_a\tproteins_b\tintra_protein\texpmass\tcalcmass\tcharge\thyperscore\tspectrum_q\tcrosslink_q";
        let rows = [
            "1\trun.mzML\tscan=5\t1\t1\tTT\tcrosslink\tGLEKVSTAAM[UNIMOD:35]QR\tTDYHFPNEWKLIR\t4\t10\tsp|P11111|ONE_HUMAN\tsp|P22222|TWO_HUMAN\tfalse\t3000.0\t3000.001\t3\t40.0\t0.0\t0.001",
            "2\trun.mzML\tscan=6\t1\t-1\tTD\tcrosslink\tGLEKVSTAAMQR\tRILKWENPFHYDT\t4\t4\tsp|P11111|ONE_HUMAN\trev_sp|P22222|TWO_HUMAN\tfalse\t3000.0\t3000.0\t3\t10.0\t0.5\t0.5",
            "3\trun.mzML\tscan=7\t1\t1\tTT\tloop\tGLEKVSTAAMQRK\t\t4\t13\tsp|P11111|ONE_HUMAN\t\ttrue\t1500.0\t1500.0\t2\t30.0\t0.0\t0.002",
        ];
        std::fs::write(dir.join("crosslinks.sage.tsv"), format!("{header}\n{}\n", rows.join("\n"))).unwrap();
        let get = |u: &Upload, e: &str| -> Value { serde_json::from_slice(u.endpoint(e).unwrap()).unwrap() };

        let passing = Upload::load(&dir, 1, Fdr::default()).unwrap();
        let matches = get(&passing, "get_xiview_matches");
        assert_eq!(matches.as_array().unwrap().len(), 2, "the decoy pair is left out: {matches}");
        let peptides = get(&passing, "get_xiview_peptides");
        let first = &matches[0];
        let donor = &peptides[first["pi1"].as_u64().unwrap() as usize];
        let acceptor = &peptides[first["pi2"].as_u64().unwrap() as usize];
        assert_eq!(donor["seq"], "GLEKVSTAAMQR");
        assert_eq!((donor["ls1"].as_u64(), donor["pos"][0].as_u64()), (Some(4), Some(3)));
        assert_eq!(donor["prt"][0], "sp|P11111|ONE_HUMAN");
        // the oxidation at M10 and the linker (donor) at K4, in location order
        assert_eq!(donor["m_ps"], json!([4, 10]));
        assert_eq!(donor["m_as"][0]["MS:1002509"], 0);
        assert_eq!(donor["m_as"][0]["XLMOD:02001"], "DSS");
        assert_eq!(donor["m_as"][1]["UNIMOD:35"], "Oxidation");
        assert!((donor["cl_m"].as_f64().unwrap() - 138.068).abs() < 1e-3);
        assert_eq!(acceptor["m_as"][0], json!({"MS:1002510": 0}));
        assert_eq!(acceptor["cl_m"], 0.0);
        assert_eq!((acceptor["pos"][0].as_u64(), acceptor["ls1"].as_u64()), (Some(2), Some(10)));
        // m/z from neutral masses
        assert!((first["pc_mz"].as_f64().unwrap() - (3000.0 + 3.0 * PROTON) / 3.0).abs() < 1e-6);
        assert_eq!(first["sc"]["XNova:hyperscore"], 40.0);
        assert_eq!(first["sp"], "scan=5");
        // the loop-link: one peptide, both sites, no second peptide
        let lp = &matches[1];
        assert!(lp["pi2"].is_null());
        let loop_peptide = &peptides[lp["pi1"].as_u64().unwrap() as usize];
        assert_eq!((loop_peptide["ls1"].as_u64(), loop_peptide["ls2"].as_u64()), (Some(4), Some(13)));
        // proteins: gene names, UniProt accessions for xiVIEW's UniProt lookup
        let proteins = get(&passing, "get_xiview_proteins");
        let one = proteins.as_array().unwrap().iter().find(|p| p["id"] == "sp|P11111|ONE_HUMAN").unwrap();
        assert_eq!((one["name"].as_str(), one["accession"].as_str()), (Some("ONE"), Some("P11111")));
        assert_eq!(one["description"], "One protein");

        // XNova's FDR: targets only, at any level
        let all = Upload::load(&dir, 2, Fdr(1.0)).unwrap();
        let matches = get(&all, "get_xiview_matches");
        assert_eq!(matches.as_array().unwrap().len(), 2, "never the decoy pair: {matches}");
        assert!(matches.as_array().unwrap().iter().all(|m| m["p"] == true));
        // the served match leads back to XNova's, and knows how xiSPEC writes its peptides
        let (m, s) = all.served(matches[0]["id"].as_str().unwrap()).unwrap();
        assert_eq!((m.scan.as_str(), s.layout.is_none()), ("scan=5", true));
        assert_eq!(s.signatures[0], vec!["GLEKVSTAAM(oxidation)QR".to_string(), "TDYHFPNEWKLIR".to_string()]);
        // no multivalent species: no link metadata
        assert!(all.endpoint("get_xnova_link_metadata").unwrap().is_empty());

        // metadata
        let sip = get(&all, "get_xiview_spectrum_identification_protocols");
        assert_eq!((sip[0]["frag_tol"].as_f64(), sip[0]["frag_tol_unit"].as_str()), (Some(20.0), Some("ppm")));
        let enzymes = get(&all, "get_xiview_enzymes");
        assert_eq!(enzymes[0]["site_regexp"], "(?<=[KR])(?!P)");
        let mods = get(&all, "get_xiview_search_modifications");
        assert!(mods.as_array().unwrap().iter().any(|m| m["accessions"]["UNIMOD:4"] == "Carbamidomethyl" && m["fixed_mod"] == true));
        assert!(mods.as_array().unwrap().iter().any(|m| m["accessions"].get("MS:1002509").is_some() && m["residues"] == "K"));
        let sd = get(&all, "get_xiview_spectra_data");
        assert_eq!(sd[0]["name"], "run.mzML");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// A multivalent species (a1b1c, DSS): one xiVIEW match per link, its
    /// links as link metadata, and the whole species in the spectrum viewer
    #[test]
    fn serves_a_multivalent_species() {
        use super::super::annotate::{annotate, Known};
        use crate::viewer::model::Peaks;
        let dir = std::env::temp_dir().join(format!("xnova-xiview-multi-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("p.fasta");
        std::fs::write(
            &fasta,
            ">sp|P11111|AAA_HUMAN A OS=Homo sapiens GN=AAA\nMKGLEKVSTAAMQRK\n\
             >sp|P22222|BBB_HUMAN B OS=Homo sapiens GN=BBB\nMRAEFGHDSTWKLMNPQSTVKGYEDFARGG\n\
             >sp|P33333|CCC_HUMAN C OS=Homo sapiens GN=CCC\nMTDYHFPNEWKLIRGG\n",
        )
        .unwrap();
        let config = json!({
            "version": "0.23.0",
            "database": {"fasta": fasta.display().to_string(), "decoy_tag": "rev_", "enzyme": {"enzymes": ["trypsin"]}},
            "crosslink": {"reagent": "dss", "residues": "K"},
            "fragment_tol": {"ppm": [-20.0, 20.0]},
            "mzml_paths": [dir.join("run.mzML").display().to_string()]
        });
        std::fs::write(dir.join("results.json"), config.to_string()).unwrap();
        std::fs::write(dir.join("run.mzML"), "").unwrap();
        let header = "psm_id\tfilename\tscannr\trank\tlabel\tclass\tlink_type\tpeptide_a\tpeptide_b\tlink_a\tlink_b\tproteins_a\tproteins_b\tintra_protein\texpmass\tcalcmass\tcharge\thyperscore\tspectrum_q\tcrosslink_q\ttopology\tspecies";
        let row = "1\trun.mzML\tscan=9\t1\t1\tTT\tmultivalent\tGLEKVSTAAMQR\tAEFGHDSTWKLMNPQSTVKGYEDFAR\t4\t10\tsp|P11111|AAA_HUMAN\tsp|P22222|BBB_HUMAN\tfalse\t5000.0\t5000.0\t4\t50.0\t0.0\t0.001\ta1b1c\t[1] GLEKVSTAAMQR [2] AEFGHDSTWKLMNPQSTVKGYEDFAR [3] TDYHFPNEWKLIR | 1:K4-2:K10; 2:K19-3:K10";
        std::fs::write(dir.join("crosslinks.sage.tsv"), format!("{header}\n{row}\n")).unwrap();

        let upload = Upload::load(&dir, 1, Fdr::default()).unwrap();
        let matches: Value = serde_json::from_slice(upload.endpoint("get_xiview_matches").unwrap()).unwrap();
        assert_eq!(matches.as_array().unwrap().len(), 2, "one match per link");
        // the link metadata: both links, one species label
        let csv = String::from_utf8(upload.endpoint("get_xnova_link_metadata").unwrap().to_vec()).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], "Protein1,SeqPos1,Protein2,SeqPos2,Multivalent,Multivalent species");
        assert_eq!(lines.len(), 3, "{csv}");
        assert!(csv.contains("\"sp|P11111|AAA_HUMAN\",6,\"sp|P22222|BBB_HUMAN\",12"), "A K6 – B K12: {csv}");
        assert!(csv.contains("\"sp|P22222|BBB_HUMAN\",21,\"sp|P33333|CCC_HUMAN\",11"), "B K21 – C K11: {csv}");
        let species = "\"a1b1c: AAA K6–BBB K12 + BBB K21–CCC K11\"";
        assert_eq!(csv.matches(species).count(), 2, "the same species label on both links: {csv}");

        // the second link's match: its link first, the whole species shown
        let (m, s) = upload.served("1_0_XL2").unwrap();
        let layout = s.layout.as_ref().unwrap();
        assert_eq!(layout.order, vec![1, 2, 0], "B and C (the link shown), then A");
        assert_eq!(layout.peptides.as_array().unwrap().len(), 3);
        assert_eq!(layout.link_sites[0], json!({"id": 0, "peptideId": 0, "linkSite": 18}));
        assert_eq!(layout.link_sites[1], json!({"id": 0, "peptideId": 1, "linkSite": 9}));
        assert_eq!(layout.link_sites[2], json!({"id": 1, "peptideId": 2, "linkSite": 3}));

        // annotated: peaks at b ions of A, one past its link (carries B and C)
        let proforma = upload.results.proforma(m);
        let set = mzcore::prelude::PeptidoformIonSet::pro_forma(&proforma, &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0;
        use mzannotate::prelude::*;
        use mzcore::prelude::*;
        let fragments = set.generate_theoretical_fragments::<mzcore::chemistry::OutputMolecularFormula>(
            mzcore::system::isize::Charge::new::<mzcore::system::e>(4),
            FragmentationModel::cid(),
        );
        let mut mz: Vec<f64> = fragments
            .iter()
            .filter(|f| f.peptidoform_index == Some(0) && f.neutral_loss.is_empty() && f.charge.value == 2)
            .filter(|f| matches!(&f.ion, mzannotate::fragment::FragmentType::b(p, _) if p.series_number == 3 || p.series_number == 4))
            .filter_map(|f| f.mz(MassMode::Monoisotopic).map(|m| m.value))
            .collect();
        mz.sort_by(f64::total_cmp);
        assert_eq!(mz.len(), 2);
        let peaks = Peaks { intensity: vec![100.0; mz.len()], mz, precursor_mz: None, ion_mobility: None };
        let request = json!({
            "Peptides": [], "LinkSite": [], "peaks": [],
            "annotation": {
                "requestID": "1_0_XL2", "precursorCharge": 4, "modifications": [],
                "fragmentTolerance": {"tolerance": 10.0, "unit": "ppm"},
                "ions": [{"type": "BIon"}, {"type": "YIon"}], "crosslinker": {"modMass": 138.06808}
            }
        });
        let out = annotate(&request, &peaks, Some(Known { results: &upload.results, m, layout: Some(layout) })).unwrap();
        assert_eq!(out["Peptides"].as_array().unwrap().len(), 3);
        assert_eq!(out["LinkSite"].as_array().unwrap().len(), 4);
        let b = |name: &str| out["fragments"].as_array().unwrap().iter().find(|f| f["name"] == name).cloned();
        let b3 = b("b3").expect("b3 of A");
        assert_eq!(b3["peptideId"], 2, "A is row 2");
        let b4 = b("b4+2P").expect("b4 of A holds K4: B and, through B, C");
        let ranges: Vec<(u64, u64, u64)> = b4["range"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| (r["peptideId"].as_u64().unwrap(), r["from"].as_u64().unwrap(), r["to"].as_u64().unwrap()))
            .collect();
        assert_eq!(ranges, vec![(0, 0, 25), (1, 0, 12), (2, 0, 3)]);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
