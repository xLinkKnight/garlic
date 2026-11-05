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
    /// Sage's protein lists (`;`-separated): of the peptide (linear), or of
    /// `peptide` and `partner` (crosslinks, loop-links, multivalent species)
    pub protein_lists: Vec<String>,
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
    /// The proteins of the search's FASTA, by accession (as Sage names them:
    /// the header's first word); empty when the FASTA cannot be read
    proteins: HashMap<String, Protein>,
    /// Prefix of decoy protein accessions (`rev_`)
    decoy_tag: String,
}

/// A protein of the search's FASTA
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Protein {
    /// `DNA-directed RNA polymerase II subunit RPB1` (UniProt headers: up to `OS=`)
    pub name: String,
    /// `RPO21` (UniProt `GN=`), empty if the header has none
    pub gene: String,
    pub sequence: String,
}

/// The FASTA file at `path` by accession (the header's first word)
fn read_fasta(path: &str) -> HashMap<String, Protein> {
    let Ok(text) = std::fs::read_to_string(path) else { return HashMap::new() };
    let mut proteins = HashMap::new();
    let mut current: Option<(String, Protein)> = None;
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty()) {
        if let Some(header) = line.strip_prefix('>') {
            proteins.extend(current.take());
            let (accession, rest) = header.split_once(char::is_whitespace).unwrap_or((header, ""));
            let name = rest[..rest.find(" OS=").unwrap_or(rest.len())].trim().to_string();
            let gene = rest.split_whitespace().find_map(|w| w.strip_prefix("GN=")).unwrap_or("").to_string();
            current = Some((accession.to_string(), Protein { name, gene, sequence: String::new() }));
        } else if let Some((_, protein)) = &mut current {
            protein.sequence.push_str(line);
        }
    }
    proteins.extend(current);
    proteins
}

/// A peptide in Sage's notation without its modifications: `M[+15.9949]K` → `MK`
fn bare(peptide: &str) -> String {
    let mut out = String::with_capacity(peptide.len());
    let mut depth = 0;
    for c in peptide.chars() {
        match c {
            '[' => depth += 1,
            ']' => depth -= 1,
            c if depth == 0 && c.is_ascii_alphabetic() => out.push(c),
            _ => {}
        }
    }
    out
}

/// What the peptide map says left of a peptide: its protein and where its
/// link residues are in it
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeptideLabel {
    /// Gene (`RPO21`), else the entry name (`RPB1_YEAST`) or accession;
    /// `+ 1` when the peptide is in more proteins; empty if not known
    pub protein: String,
    /// `DNA-directed RNA polymerase II subunit RPB1`
    pub name: String,
    /// The link residues as positions in the protein: `K1246`, `K18, K22`
    pub sites: String,
    /// Every protein and position, for the tooltip
    pub detail: String,
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

        let decoy_tag = database.and_then(|d| d.get("decoy_tag")).and_then(Value::as_str).unwrap_or("rev_").to_string();
        // the FASTA, or both crosslink databases
        let two = database.and_then(|d| d.get("crosslink_databases"));
        let fastas: Vec<&str> = match two {
            Some(two) => ["database_1", "database_2"].iter().filter_map(|k| two.get(*k).and_then(Value::as_str)).collect(),
            None => database.and_then(|d| d.get("fasta")).and_then(Value::as_str).into_iter().collect(),
        };
        let mut proteins = HashMap::new();
        for f in fastas.into_iter().filter(|f| !f.is_empty()) {
            for (accession, protein) in read_fasta(f) {
                proteins.entry(accession).or_insert(protein);
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
            proteins,
            decoy_tag,
        })
    }

    /// The label of each peptide of a match in the peptide map: its protein
    /// (gene and name, from the FASTA) and the positions of its link residues
    /// in the protein. The proteins are the ones Sage listed for the peptide;
    /// a multivalent species' other peptides are looked up in the FASTA.
    pub fn peptide_labels(&self, m: &Match) -> Vec<PeptideLabel> {
        let link = self.link(m);
        let peptides: Vec<String> = match m.kind {
            Kind::Linear | Kind::Loop => vec![m.peptide.clone()],
            Kind::Crosslink => vec![m.peptide.clone(), m.partner.clone()],
            Kind::Multi => parse_species(&m.species).0,
        };
        peptides
            .iter()
            .enumerate()
            .map(|(k, peptide)| {
                let sequence = bare(peptide);
                let listed = [&m.peptide, &m.partner]
                    .iter()
                    .position(|p| bare(p) == sequence)
                    .and_then(|i| m.protein_lists.get(i))
                    .filter(|list| !list.is_empty());
                let accessions: Vec<String> = match listed {
                    Some(list) => list.split(';').map(str::trim).filter(|a| !a.is_empty()).map(String::from).collect(),
                    None => self.proteins_with(&sequence),
                };
                let ends = link.as_ref().map(|l| l.ends_of(k)).unwrap_or_default();
                self.label(&sequence, &accessions, &ends)
            })
            .collect()
    }

    /// Accessions of the FASTA proteins holding `sequence` (a decoy, reversed
    /// as Sage reverses peptides, when no target holds it)
    fn proteins_with(&self, sequence: &str) -> Vec<String> {
        let mut found: Vec<String> =
            self.proteins.iter().filter(|(_, p)| p.sequence.contains(sequence)).map(|(a, _)| a.clone()).collect();
        if found.is_empty() {
            let target = unreverse(sequence);
            found = self
                .proteins
                .iter()
                .filter(|(_, p)| p.sequence.contains(&target))
                .map(|(a, _)| format!("{}{a}", self.decoy_tag))
                .collect();
        }
        found.sort();
        found
    }

    /// The label of a peptide (`sequence`, no modifications) in `accessions`,
    /// with link residues `ends` (0-based in the peptide)
    fn label(&self, sequence: &str, accessions: &[String], ends: &[usize]) -> PeptideLabel {
        let n = sequence.len();
        let residues: Vec<char> = sequence.chars().collect();
        // (accession, protein, decoy, 0-based start of the peptide in the protein)
        let mut places: Vec<(&str, Option<&Protein>, bool, Option<usize>)> = Vec::new();
        for accession in accessions {
            let (target, decoy) = match accession.strip_prefix(&self.decoy_tag) {
                Some(target) => (target, true),
                None => (accession.as_str(), false),
            };
            let protein = self.proteins.get(target);
            // a decoy peptide is its target peptide reversed but for its ends
            let query = if decoy { unreverse(sequence) } else { sequence.to_string() };
            let starts: Vec<usize> = protein.map(|p| p.sequence.match_indices(&query).map(|(i, _)| i).collect()).unwrap_or_default();
            if starts.is_empty() {
                places.push((accession, protein, decoy, None));
            }
            places.extend(starts.into_iter().map(|s| (accession.as_str(), protein, decoy, Some(s))));
        }
        // 1-based position in the protein of residue r of the peptide
        let position = |start: usize, decoy: bool, r: usize| {
            let i = if decoy && r > 0 && r + 1 < n { n - 1 - r } else { r };
            start + i + 1
        };
        let sites_at = |start: usize, decoy: bool| -> String {
            ends.iter().map(|&r| format!("{}{}", residues.get(r).copied().unwrap_or('?'), position(start, decoy, r))).collect::<Vec<_>>().join(", ")
        };
        let short = |accession: &str, protein: Option<&Protein>| -> String {
            match protein {
                Some(p) if !p.gene.is_empty() => p.gene.clone(),
                // UniProt `sp|P04050|RPB1_YEAST` → `RPB1_YEAST`
                _ => accession.rsplit('|').next().unwrap_or(accession).to_string(),
            }
        };
        let Some(&(accession, protein, decoy, start)) = places.first() else { return PeptideLabel::default() };
        let others = places.iter().map(|p| p.0).collect::<HashSet<_>>().len() - 1;
        let mut short_name = short(accession.strip_prefix(&self.decoy_tag).unwrap_or(accession), protein);
        if decoy {
            short_name.push_str(" (decoy)");
        }
        if others > 0 {
            short_name.push_str(&format!(" + {others}"));
        }
        let detail = places
            .iter()
            .map(|&(accession, protein, decoy, start)| {
                let what = match protein {
                    Some(p) if !p.gene.is_empty() => format!("{} ({})", p.name, p.gene),
                    Some(p) => p.name.clone(),
                    None => "not in the FASTA".into(),
                };
                let note = if decoy { ", decoy (positions of its target)" } else { "" };
                match start {
                    Some(s) if ends.is_empty() => format!("{accession}: {what}, residues {}–{}{note}", s + 1, s + n),
                    Some(s) => format!("{accession}: {what}, residues {}–{}, link {}{note}", s + 1, s + n, sites_at(s, decoy)),
                    None => format!("{accession}: {what}"),
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        PeptideLabel {
            protein: short_name,
            name: protein.map(|p| p.name.clone()).unwrap_or_default(),
            sites: start.map(|s| sites_at(s, decoy)).unwrap_or_default(),
            detail,
        }
    }

    /// The link of a crosslinked match: sites, linker and, for MS-cleavable
    /// linkers, each peptide with each stub on its linked residue (ProForma)
    pub fn link(&self, m: &Match) -> Option<Link> {
        let site = |s: usize| s.saturating_sub(1);
        let settings = self.crosslink.clone().unwrap_or_default();
        let sites = match m.kind {
            Kind::Linear => return None,
            Kind::Multi => {
                // every link of the species; MS-cleavable linkers: each peptide
                // alone, its ions with a stub on each link residue they hold
                let (peptides, links) = parse_species(&m.species);
                return Some(Link {
                    kind: Kind::Multi,
                    sites: links.iter().flat_map(|&(a, sa, b, sb)| [(a, sa), (b, sb)]).collect(),
                    links,
                    linker: linker_name(&settings),
                    linker_mass: settings.linker_mass() as f64,
                    stubs: stubs_of(&settings),
                    stub_peptides: Vec::new(),
                    peptides: peptides.iter().map(|p| rewrite(p, &[], &self.mod_masses)).collect(),
                });
            }
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
            links: vec![(sites[0].0, sites[0].1, sites[1].0, sites[1].1)],
            sites,
            linker: linker_name(&settings),
            linker_mass: settings.linker_mass() as f64,
            stubs,
            stub_peptides,
            peptides: Vec::new(),
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
    /// `progress` gets a message per file. (The viewer reads them in the
    /// background, see [`super::loader`].)
    #[cfg(test)]
    pub fn read_spectra(
        &self,
        wanted: &HashSet<SpectrumKey>,
        progress: &dyn Fn(String),
    ) -> Result<HashMap<SpectrumKey, Peaks>, String> {
        let mut by_file: HashMap<&str, HashSet<String>> = HashMap::new();
        for (file, scan) in wanted {
            by_file.entry(file.as_str()).or_default().insert(scan.clone());
        }
        let mut out = HashMap::new();
        for (file, ids) in by_file {
            progress(format!("reading {file} ({} spectra wanted)", ids.len()));
            let source = self.open_file(file)?;
            for (id, peaks) in source.read(&ids).map_err(|e| format!("{file}: {e}"))? {
                out.insert((file.to_string(), id), peaks);
            }
        }
        Ok(out)
    }

    /// Open a searched spectrum file to read spectra from it
    pub fn open_file(&self, file: &str) -> Result<SpectrumFile, String> {
        let url = self
            .files
            .get(file)
            .ok_or_else(|| format!("{file} is not among the searched files in results.json"))?;
        let bruker = sage_cloudpath::util::FileFormat::from(url.as_str()) == sage_cloudpath::util::FileFormat::TDF;
        if bruker && url.scheme() == "file" {
            let path = url.to_file_path().map_err(|_| format!("{file}: not a local path"))?;
            let run = sage_cloudpath::tdf::TdfSpectra::open(&path, 0, self.bruker).map_err(|e| format!("{file}: {e}"))?;
            return Ok(SpectrumFile::Bruker(run));
        }
        Ok(SpectrumFile::Whole(url.clone(), self.bruker))
    }

    /// Is `file` read whole (mzML, MGF), not spectrum by spectrum (Bruker)?
    pub fn read_whole(&self, file: &str) -> bool {
        self.files
            .get(file)
            .is_none_or(|url| sage_cloudpath::util::FileFormat::from(url.as_str()) != sage_cloudpath::util::FileFormat::TDF)
    }
}

/// A spectrum file opened for reading the spectra of matches
pub enum SpectrumFile {
    /// Bruker .d: each spectrum read on its own, by its id, exactly as the
    /// search read it, without decoding the rest of the run
    Bruker(sage_cloudpath::tdf::TdfSpectra),
    /// mzML, MGF: read whole (they have no index of spectra by id)
    Whole(Url, BrukerProcessingConfig),
}

impl std::fmt::Debug for SpectrumFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpectrumFile::Bruker(_) => f.write_str("timsTOF run"),
            SpectrumFile::Whole(url, _) => write!(f, "{url}"),
        }
    }
}

impl SpectrumFile {
    /// The spectra with these ids (Sage's spectrum ids), as `(id, peaks)`
    pub fn read(&self, ids: &HashSet<String>) -> Result<Vec<(String, Peaks)>, String> {
        let spectra = match self {
            SpectrumFile::Bruker(run) => {
                let numbers: Vec<usize> = ids.iter().filter_map(|id| id.parse().ok()).collect();
                run.read(&numbers)
            }
            SpectrumFile::Whole(url, bruker) => sage_cloudpath::util::read_spectra_with_ids(url, 0, None, *bruker, ids)
                .map_err(|e| e.to_string())?,
        };
        Ok(spectra
            .into_iter()
            .map(|s| {
                let precursor = s.precursors.first();
                let peaks = Peaks {
                    mz: s.mz.iter().map(|&x| x as f64).collect(),
                    precursor_mz: precursor.map(|p| p.mz as f64),
                    ion_mobility: precursor.and_then(|p| p.inverse_ion_mobility),
                    intensity: s.intensity,
                };
                (s.id, peaks)
            })
            .collect())
    }
}

/// A decoy peptide back to its target: Sage reverses all but the ends
fn unreverse(sequence: &str) -> String {
    let mut residues: Vec<char> = sequence.chars().collect();
    let n = residues.len();
    if n > 2 {
        residues[1..n - 1].reverse();
    }
    residues.into_iter().collect()
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
    /// (peptide, 0-based residue) of both ends; a loop-link has both on
    /// peptide 0; a multivalent species lists both ends of every link
    pub sites: Vec<(usize, usize)>,
    /// Every link as `(peptide, residue, peptide, residue)`, 0-based: one for
    /// a crosslink or loop-link, two or more for a multivalent species
    /// (link n is `#XLn` in the ProForma)
    pub links: Vec<(usize, usize, usize, usize)>,
    /// Display name of the linker (`DSSO`)
    pub linker: String,
    /// Mass of the whole linker (Da)
    pub linker_mass: f64,
    /// Stubs of an MS-cleavable linker (both arms), empty otherwise
    pub stubs: Vec<Stub>,
    /// Peptide → stub → ProForma of the peptide alone with that stub on its
    /// linked residue (crosslinks and loop-links)
    pub stub_peptides: Vec<Vec<String>>,
    /// Multivalent species: each peptide alone, as ProForma
    pub peptides: Vec<String>,
}

impl Link {
    /// The link residues of `peptide`, in sequence order
    pub fn ends_of(&self, peptide: usize) -> Vec<usize> {
        let mut ends: Vec<usize> = self.sites.iter().filter(|s| s.0 == peptide).map(|s| s.1).collect();
        ends.sort_unstable();
        ends
    }
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
            protein_lists: vec![field(&f, &col, "proteins").to_string()],
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
                protein_lists: vec![field(&f, &col, "proteins_a").to_string(), field(&f, &col, "proteins_b").to_string()],
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
    /// The whole peptide with the stub (the doublet peaks); in a multivalent
    /// species, with a stub on each of its link residues (triplets with two)
    Peptide,
    /// N-terminal (a/b/c) ion ending after this residue, holding the linked one
    NTerm(usize),
    /// C-terminal (x/y/z) ion starting at this residue, holding the linked one
    CTerm(usize),
    /// Multivalent species: the peptides (bit per peptide) still held
    /// together after only link `link` broke, with its stub on the end of
    /// [`StubIon::peptide`] (`{2,3}+A` from `a1b1c`). These show which
    /// peptides were linked to which.
    Part { link: u8, peptides: u8 },
}

/// An ion with a linker stub found in the spectrum
#[derive(Clone, Debug)]
pub struct StubIon {
    pub peptide: usize,
    /// Index into [`Link::stubs`] (the first of [`StubIon::stubs`])
    pub stub: usize,
    /// Every stub it carries, sorted: one, or in a multivalent species one
    /// per link residue it holds (`y9+A+S`)
    pub stubs: Vec<usize>,
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
                found.push(StubIon { peptide: k, stub: j, stubs: vec![j], kind: StubIonKind::Peptide, charge: z, theoretical_mz: mz, peak, text });
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
                found.push(StubIon { peptide: k, stub: j, stubs: vec![j], kind, charge: z, theoretical_mz: mz, peak, text });
            }
        }
    }
    found
}

const PROTON: f64 = 1.007_276_466_6;

/// `²+` for charge 2, nothing for 1
fn charge_text(z: u8) -> String {
    if z > 1 {
        format!("{}+", superscript(z as usize))
    } else {
        String::new()
    }
}

/// Every choice of `k` stubs out of `n` kinds, repeats allowed, sorted
/// (`A+A`, `A+S`, `S+S`...)
fn stub_choices(n: usize, k: usize) -> Vec<Vec<usize>> {
    fn walk(n: usize, k: usize, from: usize, current: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if current.len() == k {
            out.push(current.clone());
            return;
        }
        for j in from..n {
            current.push(j);
            walk(n, k, j, current, out);
            current.pop();
        }
    }
    let mut out = Vec::new();
    walk(n, k, 0, &mut Vec::new(), &mut out);
    out
}

/// The peptides joined to `start` by the links other than `skip`
fn joined(start: usize, links: &[(usize, usize, usize, usize)], skip: usize) -> Vec<usize> {
    let mut side = vec![start];
    let mut grow = true;
    while grow {
        grow = false;
        for (_, &(a, _, b, _)) in links.iter().enumerate().filter(|(e, _)| *e != skip) {
            for (from, to) in [(a, b), (b, a)] {
                if side.contains(&from) && !side.contains(&to) {
                    side.push(to);
                    grow = true;
                }
            }
        }
    }
    side
}

/// Stub ions of a multivalent species whose MS-cleavable linkers broke, as
/// the search scores them. Once every linker is broken each peptide is free
/// and keeps one stub on each of its link residues: the whole peptide with
/// its stubs (doublets, or triplets with two links), and its backbone ions
/// with one stub per link residue they hold. Parts held together by the
/// other linkers when only one broke (`b1c` + stub from `a1b1c`) show which
/// peptides were linked to which. Asymmetric linkers: either arm's stubs on
/// any end (the species does not say which way round each linker sits).
fn find_multi_stub_ions(views: &mut [PeakView], link: &Link, charge: u8, model: Model, tolerance: Tolerance) -> Vec<StubIon> {
    let mut found: Vec<StubIon> = Vec::new();
    // each piece holds at least one charge less than the species
    let max_charge = charge.saturating_sub(1).max(1);
    let names = |stubs: &[usize]| stubs.iter().map(|&j| link.stubs[j].name.as_str()).collect::<Vec<_>>().join("+");
    let added = |stubs: &[usize]| stubs.iter().map(|&j| link.stubs[j].mass).sum::<f64>();
    let mut masses: Vec<Option<f64>> = Vec::new();
    for (x, proforma) in link.peptides.iter().enumerate() {
        let Ok((ions, _)) = PeptidoformIonSet::pro_forma(proforma, &mzcore::ontology::STATIC_ONTOLOGIES) else {
            masses.push(None);
            continue;
        };
        let mass = ions.formulas().first().map(|f| f.monoisotopic_mass().value);
        masses.push(mass);
        let ends = link.ends_of(x);
        if ends.is_empty() {
            continue;
        }
        let n = ions
            .peptidoform_ions()
            .first()
            .and_then(|i| i.peptidoforms().first())
            .map_or(0, |p| p.len());

        // the whole peptide, a stub on each link residue
        if let Some(mass) = mass {
            for stubs in stub_choices(link.stubs.len(), ends.len()) {
                for z in 1..=max_charge {
                    let mz = (mass + added(&stubs) + z as f64 * PROTON) / z as f64;
                    let Some(peak) = find_peak(views, mz, tolerance) else { continue };
                    let text = format!("P+{}{}", names(&stubs), charge_text(z));
                    add_stub_label(views, peak, mz, x, stubs[0], &text);
                    found.push(StubIon {
                        peptide: x,
                        stub: stubs[0],
                        stubs: stubs.clone(),
                        kind: StubIonKind::Peptide,
                        charge: z,
                        theoretical_mz: mz,
                        peak,
                        text,
                    });
                }
            }
        }

        // backbone ions holding link residues: a stub on each one held
        let fragments = ions.generate_theoretical_fragments::<OutputMolecularFormula>(
            Charge::new::<mzcore::system::e>(max_charge as isize),
            model.model(),
        );
        for f in fragments.iter().filter(|f| f.neutral_loss.is_empty()) {
            // `bN` holds residues 0..N, `yN` n-N..n
            let (kind, held) = match &f.ion {
                FragmentType::a(p, _) | FragmentType::b(p, _) | FragmentType::c(p, _) => (
                    StubIonKind::NTerm(p.series_number.saturating_sub(1)),
                    ends.iter().filter(|&&r| r < p.series_number).count(),
                ),
                FragmentType::x(p, _) | FragmentType::y(p, _) | FragmentType::z(p, _) if p.series_number <= n => (
                    StubIonKind::CTerm(n - p.series_number),
                    ends.iter().filter(|&&r| r >= n - p.series_number).count(),
                ),
                _ => continue,
            };
            if held == 0 {
                continue;
            }
            let Some(base) = f.mz(MassMode::Monoisotopic).map(|m| m.value) else { continue };
            let z = f.charge.value.unsigned_abs().max(1) as u8;
            // `y6²+` → `y6+A+S²+`: the stubs before the charge
            let text = fragment_text(f);
            let (ion, charge_part) = match text.find(|c: char| "⁰¹²³⁴⁵⁶⁷⁸⁹".contains(c)) {
                Some(i) => (text[..i].to_string(), text[i..].to_string()),
                None => (text.clone(), String::new()),
            };
            for stubs in stub_choices(link.stubs.len(), held) {
                let mz = base + added(&stubs) / z as f64;
                let Some(peak) = find_peak(views, mz, tolerance) else { continue };
                if found.iter().any(|s| s.peptide == x && s.kind == kind && s.charge == z && s.stubs == stubs) {
                    continue;
                }
                let text = format!("{ion}+{}{charge_part}", names(&stubs));
                add_stub_label(views, peak, mz, x, stubs[0], &text);
                found.push(StubIon {
                    peptide: x,
                    stub: stubs[0],
                    stubs: stubs.clone(),
                    kind,
                    charge: z,
                    theoretical_mz: mz,
                    peak,
                    text,
                });
            }
        }
    }

    // parts freed when only one linker broke: link e holds two parts together
    // when nothing else does; a part of two or more peptides keeps the
    // linkers inside it and gets one stub where e was
    let peptides = link.peptides.len();
    if masses.iter().all(Option::is_some) && peptides <= 8 {
        for (e, &(a, _, b, _)) in link.links.iter().enumerate() {
            let side = joined(a, &link.links, e);
            if side.contains(&b) {
                continue;
            }
            for end in [a, b] {
                let part: Vec<usize> = (0..peptides).filter(|v| side.contains(v) == side.contains(&end)).collect();
                if part.len() < 2 {
                    // one peptide on its own: its peptide + stub peaks, above
                    continue;
                }
                let inner = link
                    .links
                    .iter()
                    .enumerate()
                    .filter(|&(i, l)| i != e && part.contains(&l.0) && part.contains(&l.2))
                    .count();
                let mass = part.iter().filter_map(|&v| masses[v]).sum::<f64>() + link.linker_mass * inner as f64;
                let mask = part.iter().fold(0u8, |m, &v| m | (1 << v));
                let members = part.iter().map(|v| (v + 1).to_string()).collect::<Vec<_>>().join(",");
                for (j, stub) in link.stubs.iter().enumerate() {
                    for z in 1..=max_charge {
                        let mz = (mass + stub.mass + z as f64 * PROTON) / z as f64;
                        let Some(peak) = find_peak(views, mz, tolerance) else { continue };
                        let text = format!("{{{members}}}+{}{}", stub.name, charge_text(z));
                        add_stub_label(views, peak, mz, end, j, &text);
                        found.push(StubIon {
                            peptide: end,
                            stub: j,
                            stubs: vec![j],
                            kind: StubIonKind::Part { link: e as u8, peptides: mask },
                            charge: z,
                            theoretical_mz: mz,
                            peak,
                            text,
                        });
                    }
                }
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
        Some(link) if !link.stubs.is_empty() && link.kind == Kind::Multi => {
            find_multi_stub_ions(&mut views, link, charge, model, tolerance)
        }
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

    /// SILAC peptides carry their labels as UNIMOD accessions (a label on a
    /// modified residue as a mass), which Annotator's engine resolves: the
    /// heavy peptide is 8.0142 + 10.0083 Da heavier, and its y ions carry the
    /// label of the C-terminal R
    #[test]
    fn silac_labels_annotate_at_their_masses() {
        let mass = |p: &str| {
            PeptidoformIonSet::pro_forma(&rewrite(p, &[], &HashMap::new()), &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0.formulas()[0]
                .monoisotopic_mass()
                .value
        };
        let light = mass("PEPTIDEKAR");
        assert!((mass("PEPTIDEK[UNIMOD:259]AR[UNIMOD:267]") - light - 18.022468).abs() < 1e-3);
        assert!((mass("PEPTIDEK[UNIMOD:481]AR[UNIMOD:188]") - light - 10.045236).abs() < 1e-3);
        // acetyl-K, heavy: one combined mass
        assert!((mass("PEPTIDEK[+50.024764]AR") - light - 50.024764).abs() < 1e-3);
        // annotated: the heavy y1 (R + 10.0083) is matched
        let y1 = 175.118952 + 10.008269;
        let peaks = Peaks { mz: vec![y1], intensity: vec![100.0], precursor_mz: None, ion_mobility: None };
        let a = annotate(&peaks, "PEPTIDEK[UNIMOD:259]AR[UNIMOD:267]", 2, Model::CidHcd, Tolerance::Ppm(-10.0, 10.0), None).unwrap();
        assert!(a.peaks[0].labels.iter().any(|l| l.text == "y1"), "{:?}", a.peaks[0].labels);
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
        Link {
            kind: Kind::Crosslink,
            links: vec![(0, 3, 1, 9)],
            sites,
            linker: "DSSO".into(),
            linker_mass: 158.0038,
            stubs,
            stub_peptides,
            peptides: Vec::new(),
        }
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
            proteins: HashMap::new(),
            decoy_tag: "rev_".into(),
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
            protein_lists: Vec::new(),
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
        // DSS: the links, no stubs
        let link = results.link(&m).unwrap();
        assert_eq!(link.links, [(0, 3, 1, 9), (1, 18, 2, 9)]);
        assert_eq!(link.ends_of(1), [9, 18]);
        assert!(link.stubs.is_empty());
        let peaks = Peaks { mz: vec![300.0, 500.0, 700.0], intensity: vec![10.0; 3], precursor_mz: None, ion_mobility: None };
        let a = annotate(&peaks, &proforma, 4, Model::CidHcd, Tolerance::Ppm(-10.0, 10.0), Some(&link)).unwrap();
        assert_eq!(a.coverage.len(), 3);
        assert!(a.stub_ions.is_empty());
    }

    /// A crosslink match for the label tests
    fn crosslink(kind: Kind, peptide: &str, partner: &str, sites: (usize, usize), lists: [&str; 2], species: &str) -> Match {
        Match {
            kind,
            file: "a.d".into(),
            scan: "1".into(),
            charge: 4,
            peptide: peptide.into(),
            partner: partner.into(),
            sites,
            proteins: String::new(),
            protein_lists: lists.iter().map(|l| l.to_string()).collect(),
            score: 50.0,
            spectrum_q: 0.0,
            q: 0.0,
            target: true,
            rank: 1,
            expmass: 0.0,
            calcmass: 0.0,
            rt: 10.0,
            ion_mobility: None,
            intra_protein: false,
            topology: String::new(),
            species: species.into(),
        }
    }

    /// The peptide map's labels: gene, protein name and the link residues'
    /// positions in the protein, from the search's FASTA
    #[test]
    fn peptide_labels_name_proteins_and_link_positions() {
        let dir = std::env::temp_dir().join(format!("sage-gui-labels-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("two.fasta");
        std::fs::write(
            &fasta,
            ">sp|P04050|RPB1_YEAST DNA-directed RNA polymerase II subunit RPB1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=RPO21 PE=1 SV=2\n\
             MAAAKGLEKVSTAAMQR\nCCCPPPPKR\n\n\
             >sp|Q00001|TEST_YEAST Test protein OS=Yeast\nWWTDYHFPNEWKLIRGLEKVSTAAMQR\n",
        )
        .unwrap();
        let proteins = read_fasta(&fasta.display().to_string());
        let rpb1 = &proteins["sp|P04050|RPB1_YEAST"];
        assert_eq!(rpb1.name, "DNA-directed RNA polymerase II subunit RPB1");
        assert_eq!(rpb1.gene, "RPO21");
        assert_eq!(rpb1.sequence, "MAAAKGLEKVSTAAMQRCCCPPPPKR");
        assert_eq!(proteins["sp|Q00001|TEST_YEAST"].gene, "");
        std::fs::remove_dir_all(&dir).unwrap();

        let results = Results {
            matches: Vec::new(),
            files: HashMap::new(),
            bruker: Default::default(),
            crosslink: Some(CrosslinkSettings::default()),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            mod_masses: HashMap::new(),
            proteins,
            decoy_tag: "rev_".into(),
        };
        // GLEKVSTAAMQR (K4) is residues 6-17 of RPB1 and 16-27 of the test
        // protein; TDYHFPNEWKLIR (K10) residues 3-15 of the test protein
        let m = crosslink(
            Kind::Crosslink,
            "GLEKVSTAAMQR",
            "TDYHFPNEWKLIR",
            (4, 10),
            ["sp|P04050|RPB1_YEAST;sp|Q00001|TEST_YEAST", "sp|Q00001|TEST_YEAST"],
            "",
        );
        let labels = results.peptide_labels(&m);
        assert_eq!(labels.len(), 2);
        assert_eq!(labels[0].protein, "RPO21 + 1", "gene, and one more protein");
        assert_eq!(labels[0].name, "DNA-directed RNA polymerase II subunit RPB1");
        assert_eq!(labels[0].sites, "K9");
        assert!(labels[0].detail.contains("sp|Q00001|TEST_YEAST: Test protein, residues 16–27, link K19"), "{}", labels[0].detail);
        // no gene name: the UniProt entry name
        assert_eq!((labels[1].protein.as_str(), labels[1].sites.as_str()), ("TEST_YEAST", "K12"));

        // a decoy (Sage reverses all but the ends): its target's positions
        let decoy = crosslink(Kind::Crosslink, "GQMAATSVKELR", "TDYHFPNEWKLIR", (9, 10), ["rev_sp|P04050|RPB1_YEAST", "sp|Q00001|TEST_YEAST"], "");
        let labels = results.peptide_labels(&decoy);
        assert_eq!((labels[0].protein.as_str(), labels[0].sites.as_str()), ("RPO21 (decoy)", "K9"));

        // a loop-link: both residues
        let looped = crosslink(Kind::Loop, "AAAKGLEKVSTAAMQR", "", (4, 8), ["sp|P04050|RPB1_YEAST", ""], "");
        assert_eq!(results.peptide_labels(&looped)[0].sites, "K5, K9");

        // multivalent: the third peptide is not in the protein columns, so it
        // is looked up in the FASTA (residues 18-26 of RPB1)
        let species = "[1] GLEKVSTAAMQR [2] TDYHFPNEWKLIR [3] CCCPPPPKR | 1:K4-2:K10; 2:T1-3:K8";
        let multi = crosslink(Kind::Multi, "GLEKVSTAAMQR", "TDYHFPNEWKLIR", (4, 10), ["sp|P04050|RPB1_YEAST", "sp|Q00001|TEST_YEAST"], species);
        let labels = results.peptide_labels(&multi);
        let sites: Vec<&str> = labels.iter().map(|l| l.sites.as_str()).collect();
        assert_eq!(sites, ["K9", "T3, K12", "K25"]);
        assert_eq!(labels[2].protein, "RPO21");
        // without the FASTA: no label (the map shows "peptide 3" alone)
        let bare_results = Results { proteins: HashMap::new(), ..results };
        assert_eq!(bare_results.peptide_labels(&multi)[2], PeptideLabel::default());
    }

    pub(crate) const MULTI_SPECIES: &str =
        "[1] GLEKVSTAAMQR [2] AEFGHDSTWKLMNPQSTVKGYEDFAR [3] TDYHFPNEWKLIR | 1:K4-2:K10; 2:K19-3:K10";

    /// A DSSO `a1b1c` species (peptide 2 holds both links) after its linkers
    /// broke: peptide 1 + A / S (doublet), peptide 2 + A+S and S+S (triplet,
    /// 2+), peptide 3 + A; y8 of peptide 2 with S (on K19 only) and y17 with
    /// A+S (both link residues, 2+); and the parts freed when one linker
    /// broke: {2,3} + A (XL1) and {1,2} + S (XL2), both 2+
    pub(crate) fn cleaved_multi() -> (Link, Annotation) {
        let settings = CrosslinkSettings { reagent: Some(sage_core::crosslink::Reagent::Dsso), ..Default::default() };
        let results = Results {
            matches: Vec::new(),
            files: HashMap::new(),
            bruker: Default::default(),
            crosslink: Some(settings.clone()),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            mod_masses: HashMap::new(),
            proteins: HashMap::new(),
            decoy_tag: "rev_".into(),
        };
        let (peptides, _) = parse_species(MULTI_SPECIES);
        let m = Match {
            kind: Kind::Multi,
            file: "a.d".into(),
            scan: "1".into(),
            charge: 4,
            peptide: peptides[0].clone(),
            partner: peptides[1].clone(),
            sites: (4, 10),
            proteins: String::new(),
            protein_lists: Vec::new(),
            score: 50.0,
            spectrum_q: 0.0,
            q: 0.0,
            target: true,
            rank: 1,
            expmass: 0.0,
            calcmass: 0.0,
            rt: 10.0,
            ion_mobility: None,
            intra_protein: false,
            topology: "a1b1c".into(),
            species: MULTI_SPECIES.into(),
        };
        let link = results.link(&m).unwrap();
        let stub = |name: &str| link.stubs.iter().find(|s| s.name == name).unwrap().mass;
        let mass = |p: &str| {
            PeptidoformIonSet::pro_forma(p, &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0.formulas()[0]
                .monoisotopic_mass()
                .value
        };
        let mz = |mass: f64, z: f64| (mass + z * PROTON) / z;
        let (m1, m2, m3) = (mass(&peptides[0]), mass(&peptides[1]), mass(&peptides[2]));
        let linker = settings.linker_mass() as f64;
        let y = |number: usize, z: isize| {
            PeptidoformIonSet::pro_forma(&peptides[1], &mzcore::ontology::STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .generate_theoretical_fragments::<OutputMolecularFormula>(Charge::new::<mzcore::system::e>(2), FragmentationModel::cid())
                .into_iter()
                .find(|f| {
                    f.neutral_loss.is_empty()
                        && f.charge.value == z
                        && matches!(&f.ion, FragmentType::y(p, _) if p.series_number == number)
                })
                .and_then(|f| f.mz(MassMode::Monoisotopic))
                .unwrap()
                .value
        };
        let mut peaks = vec![
            mz(m1 + stub("A"), 1.0),
            mz(m1 + stub("S"), 1.0),
            mz(m2 + stub("A") + stub("S"), 2.0),
            mz(m2 + 2.0 * stub("S"), 2.0),
            mz(m3 + stub("A"), 1.0),
            y(8, 1) + stub("S"),
            y(17, 2) + (stub("A") + stub("S")) / 2.0,
            mz(m2 + m3 + linker + stub("A"), 2.0),
            mz(m1 + m2 + linker + stub("S"), 2.0),
        ];
        peaks.sort_by(f64::total_cmp);
        let peaks = Peaks { intensity: vec![50.0; peaks.len()], mz: peaks, precursor_mz: None, ion_mobility: None };
        let a = annotate(&peaks, &results.proforma(&m), 4, Model::CidHcd, Tolerance::Ppm(-10.0, 10.0), Some(&link)).unwrap();
        (link, a)
    }

    #[test]
    fn multivalent_stub_ions_of_a_cleaved_species() {
        let (link, a) = cleaved_multi();
        assert_eq!(link.kind, Kind::Multi);
        let names = |stubs: &[usize]| stubs.iter().map(|&j| link.stubs[j].name.clone()).collect::<Vec<_>>().join("+");
        let found: HashSet<(usize, String, String)> = a
            .stub_ions
            .iter()
            .map(|s| (s.peptide, format!("{:?}", s.kind), names(&s.stubs)))
            .collect();
        let has = |peptide: usize, kind: StubIonKind, stubs: &str| {
            found.contains(&(peptide, format!("{kind:?}"), stubs.to_string()))
        };
        // doublets of the peptides with one link, triplet peaks of the one with two
        assert!(has(0, StubIonKind::Peptide, "A") && has(0, StubIonKind::Peptide, "S"), "{found:?}");
        assert!(has(1, StubIonKind::Peptide, "A+S") && has(1, StubIonKind::Peptide, "S+S"), "{found:?}");
        assert!(has(2, StubIonKind::Peptide, "A"), "{found:?}");
        // y8 holds K19 only (one stub), y17 holds K10 and K19 (two)
        assert!(has(1, StubIonKind::CTerm(18), "S"), "{found:?}");
        assert!(has(1, StubIonKind::CTerm(9), "A+S"), "{found:?}");
        // the parts freed when one linker broke
        assert!(has(1, StubIonKind::Part { link: 0, peptides: 0b110 }, "A"), "{found:?}");
        assert!(has(1, StubIonKind::Part { link: 1, peptides: 0b011 }, "S"), "{found:?}");
        // the labels say it all
        let labels: Vec<&str> = a.peaks.iter().map(|p| p.labels[0].text.as_str()).collect();
        for expected in ["P+A+S²+", "P+S+S²+", "y8+S", "y17+A+S²+", "{2,3}+A²+", "{1,2}+S²+"] {
            assert!(labels.contains(&expected), "{expected} in {labels:?}");
        }
        assert!(a.peaks.iter().all(|p| matches!(p.labels[0].series, Series::Stub(_))), "{labels:?}");
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
        let link = Link {
            kind: Kind::Loop,
            sites: vec![(0, first), (0, second)],
            links: vec![(0, first, 0, second)],
            linker: "DSSO".into(),
            linker_mass: 158.0038,
            stubs,
            stub_peptides,
            peptides: Vec::new(),
        };

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

    /// Where the time of annotating goes, for the matches whose title starts
    /// with one of `SAGE_GUI_TITLES` (`;`-separated): `SAGE_GUI_RESULTS=<folder>`
    #[test]
    #[ignore]
    fn annotation_timing() {
        use std::time::Instant;
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let titles = std::env::var("SAGE_GUI_TITLES").unwrap_or_default();
        let results = Results::load(Path::new(&folder)).unwrap();
        for wanted in titles.split(';').filter(|t| !t.is_empty()) {
            let Some(m) = results.matches.iter().find(|m| m.target && m.title().starts_with(wanted)) else {
                println!("not found: {wanted}");
                continue;
            };
            let ids: HashSet<String> = [m.scan.clone()].into();
            let t = Instant::now();
            let source = results.open_file(&m.file).unwrap();
            let opened = t.elapsed();
            let t = Instant::now();
            let peaks = source.read(&ids).unwrap().remove(0).1;
            let read = t.elapsed();
            let proforma = results.proforma(m);
            let link = results.link(m);
            let t = Instant::now();
            let ions = PeptidoformIonSet::pro_forma(&proforma, &mzcore::ontology::STATIC_ONTOLOGIES).unwrap().0;
            let parsed = t.elapsed();
            let t = Instant::now();
            let fragments = ions.generate_theoretical_fragments::<OutputMolecularFormula>(
                Charge::new::<mzcore::system::e>(m.charge as isize),
                Model::CidHcd.model(),
            );
            let generated = t.elapsed();
            let t = Instant::now();
            let a = annotate(&peaks, &proforma, m.charge, Model::CidHcd, results.fragment_tol, link.as_ref()).unwrap();
            let whole = t.elapsed();
            let t = Instant::now();
            let mut views: Vec<PeakView> =
                peaks.mz.iter().zip(&peaks.intensity).map(|(&mz, &intensity)| PeakView { mz, intensity, labels: Vec::new() }).collect();
            let stubs = match &link {
                Some(l) if l.kind == Kind::Multi && !l.stubs.is_empty() => {
                    find_multi_stub_ions(&mut views, l, m.charge, Model::CidHcd, results.fragment_tol).len()
                }
                Some(l) if !l.stubs.is_empty() => find_stub_ions(&mut views, l, m.charge, Model::CidHcd, results.fragment_tol).len(),
                _ => 0,
            };
            let stub_time = t.elapsed();
            println!(
                "{} ({}+, {} peaks)\n  open run {opened:.2?}, read {read:.2?}, parse {parsed:.2?}, fragments {} in {generated:.2?}, \
                 whole annotate {whole:.2?}, stub ions {stubs} in {stub_time:.2?}, {} of {} peaks annotated",
                m.title(),
                m.charge,
                peaks.mz.len(),
                fragments.len(),
                a.explained_peaks,
                a.peaks.len()
            );
        }
    }

    /// Memory and time of keeping every run open (`SAGE_GUI_RESULTS`)
    #[test]
    #[ignore]
    fn every_run_open() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let results = Results::load(Path::new(&folder)).unwrap();
        let memory = || {
            let out = std::process::Command::new("tasklist")
                .args(["/FI", &format!("PID eq {}", std::process::id()), "/FO", "CSV", "/NH"])
                .output()
                .unwrap();
            // the last CSV field, `"1,234,567 K"`
            let text = String::from_utf8_lossy(&out.stdout).trim().trim_end_matches('"').to_string();
            text.rsplit("\",\"").next().unwrap_or("").to_string()
        };
        let mut files: Vec<&String> = results.files.keys().collect();
        files.sort();
        println!("before: {}", memory());
        let t = std::time::Instant::now();
        let mut open = Vec::new();
        for f in &files {
            let t1 = std::time::Instant::now();
            open.push(results.open_file(f).unwrap());
            println!("{f}: {:.2?}, process {}", t1.elapsed(), memory());
        }
        println!("{} runs open in {:.2?}: {}", open.len(), t.elapsed(), memory());
        drop(open);
        // two at a time
        let t = std::time::Instant::now();
        std::thread::scope(|s| {
            for half in files.chunks(files.len().div_ceil(2)) {
                s.spawn(|| half.iter().for_each(|f| drop(results.open_file(f).unwrap())));
            }
        });
        println!("opened again, two at a time: {:.2?}", t.elapsed());
    }

    /// Bruker spectra read by id are the ones the whole-run read gives, peak
    /// for peak: `SAGE_GUI_RESULTS=<folder>` (a search of timsTOF data)
    #[test]
    #[ignore]
    fn spectra_read_by_id_are_the_whole_run_spectra() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let results = Results::load(Path::new(&folder)).unwrap();
        let file = results.matches.iter().find(|m| m.target && m.spectrum_q <= 0.01).unwrap().file.clone();
        let ids: HashSet<String> = results
            .matches
            .iter()
            .filter(|m| m.file == file && m.target && m.spectrum_q <= 0.01)
            .map(|m| m.scan.clone())
            .collect();
        let t = std::time::Instant::now();
        let source = results.open_file(&file).unwrap();
        assert!(matches!(source, SpectrumFile::Bruker(_)));
        let opened = t.elapsed();
        let by_id: HashMap<String, Peaks> = source.read(&ids).unwrap().into_iter().collect();
        let by_id_time = t.elapsed();
        let one: HashSet<String> = ids.iter().take(1).cloned().collect();
        let t1 = std::time::Instant::now();
        let single = source.read(&one).unwrap();
        println!(
            "{file}: opened in {opened:.2?}, {} spectra by id in {by_id_time:.2?}, one spectrum in {:.2?}",
            by_id.len(),
            t1.elapsed()
        );
        assert_eq!(single.len(), 1);

        let url = &results.files[&file];
        let t = std::time::Instant::now();
        let whole: HashMap<String, Peaks> = sage_cloudpath::util::read_spectra(url, 0, None, results.bruker, false)
            .unwrap()
            .into_iter()
            .filter(|s| ids.contains(&s.id))
            .map(|s| {
                let precursor = s.precursors.first();
                let peaks = Peaks {
                    mz: s.mz.iter().map(|&x| x as f64).collect(),
                    precursor_mz: precursor.map(|p| p.mz as f64),
                    ion_mobility: precursor.and_then(|p| p.inverse_ion_mobility),
                    intensity: s.intensity,
                };
                (s.id, peaks)
            })
            .collect();
        println!("whole run read in {:.2?}", t.elapsed());
        assert_eq!(by_id.len(), ids.len(), "every spectrum found by its id");
        assert_eq!(whole.len(), by_id.len());
        let mut peaks = 0;
        for (id, p) in &whole {
            let q = &by_id[id];
            assert_eq!(p.mz, q.mz, "spectrum {id}: m/z");
            assert_eq!(p.intensity, q.intensity, "spectrum {id}: intensities");
            assert_eq!(p.precursor_mz, q.precursor_mz);
            assert_eq!(p.ion_mobility, q.ion_mobility);
            peaks += p.mz.len();
        }
        println!("{} spectra identical, {:.0} peaks per spectrum", whole.len(), peaks as f64 / whole.len() as f64);
    }

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
