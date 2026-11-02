//! Crosslink results for downstream viewers, written next to `crosslinks.sage.tsv`:
//!
//! * `crosslinks_xiView.csv` for xiVIEW (CLMS-UI), in the CSV layout read by
//!   its importer (Rappsilber-Laboratory/mzidentml-reader, `csv_parser`).
//! * `crosslinks_proXL.xml` for ProXL, valid against the ProXL XML schema
//!   (yeastrc/proxl-import-api, `xsd/proxl-xml.xsd`).
//!
//! Neither Sage's database nor its matches know where a peptide starts in its
//! protein, so the protein sequences are taken from the FASTA.

use sage_core::crosslink::{CrosslinkMatch, CrosslinkSettings, LinkKind, PairClass, Reagent};
use sage_core::database::IndexedDatabase;
use sage_core::ion_series::Kind;
use sage_core::mass::{Tolerance, PROTON};
use sage_core::modification::{ModId, ModificationSpecificity};
use sage_core::peptide::Peptide;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write;

/// Matches at or below this crosslink q-value pass (xiVIEW `PassThreshold`,
/// ProXL default filter)
pub const Q_CUTOFF: f32 = 0.01;

/// Spectrum file formats, as far as the viewers care
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectrumKind {
    /// Read by both xiVIEW and ProXL
    MzML,
    /// Read by xiVIEW only, which finds spectra by position (`index=N`)
    Mgf,
    /// Bruker .d and anything else: neither viewer reads it
    Other,
}

/// How the viewers find a match's spectrum in its file
#[derive(Clone, Debug)]
pub struct SpectrumFile {
    /// Name to upload the file under: `.gz` is removed, as both viewers read
    /// the uncompressed file (xiVIEW also tries `name.gz`)
    pub name: String,
    pub kind: SpectrumKind,
    /// Spectrum id → 0-based position in the file, when it had to be read
    /// (MGF files; mzML whose ids carry no scan number)
    pub positions: Option<HashMap<String, usize>>,
}

impl SpectrumFile {
    /// Name and format from a spectrum file name (`run.mzML.gz`, `run.mgf`, `run.d`)
    pub fn from_name(file: &str) -> Self {
        let lower = file.to_ascii_lowercase();
        let name = if lower.ends_with(".gz") { &file[..file.len() - 3] } else { file };
        let lower = name.to_ascii_lowercase();
        let kind = if lower.ends_with(".mzml") {
            SpectrumKind::MzML
        } else if lower.ends_with(".mgf") {
            SpectrumKind::Mgf
        } else {
            SpectrumKind::Other
        };
        SpectrumFile { name: name.to_string(), kind, positions: None }
    }

    /// `ScanId` for xiVIEW, which reads mzML spectra by their id and MGF
    /// spectra by position (`index=N`)
    pub fn xiview_scan_id(&self, spec_id: &str) -> String {
        match (self.kind, &self.positions) {
            (SpectrumKind::Mgf, Some(pos)) => pos
                .get(spec_id)
                .map_or_else(|| spec_id.to_string(), |i| format!("index={i}")),
            _ => spec_id.to_string(),
        }
    }

    /// `(scan_file_name, scan_number)` for ProXL, which reads mzML / mzXML only
    pub fn proxl_scan(&self, spec_id: &str) -> Option<(&str, Option<u64>)> {
        (self.kind == SpectrumKind::MzML).then(|| {
            let position = self.positions.as_ref().and_then(|p| p.get(spec_id).copied());
            (self.name.as_str(), scan_number(spec_id, position))
        })
    }

    /// Does finding these spectra by position need the file read again?
    pub fn needs_positions<'m>(&self, mut ids: impl Iterator<Item = &'m str>) -> bool {
        match self.kind {
            SpectrumKind::Mgf => true,
            SpectrumKind::MzML => ids.any(|id| id_scan_number(id).is_none()),
            SpectrumKind::Other => false,
        }
    }
}

/// Everything the exporters need besides the matches
pub struct ExportContext<'a> {
    pub db: &'a IndexedDatabase,
    pub settings: &'a CrosslinkSettings,
    /// Spectrum file of each `file_id`
    pub files: &'a [SpectrumFile],
    /// Target protein accession (as in the FASTA) → sequence
    pub proteins: &'a HashMap<String, String>,
    pub fasta_name: &'a str,
    pub fragment_tol: Tolerance,
    pub ion_kinds: &'a [Kind],
    pub static_mods: &'a [(ModificationSpecificity, ModId)],
    /// `(site, modification, per-peptide limit)`
    pub variable_mods: &'a [(ModificationSpecificity, ModId, Option<usize>)],
    pub precursor_tol: Tolerance,
    pub version: &'a str,
}

/// The target sequence a peptide came from: decoys are targets with their
/// interior residues reversed (`Peptide::reverse`), which undoes itself
pub(super) fn target_sequence(p: &Peptide) -> Vec<u8> {
    let mut s = p.sequence.to_vec();
    let n = s.len().saturating_sub(1);
    if p.decoy && n > 1 {
        s[1..n].reverse();
    }
    s
}

/// `(protein name as reported, 1-based start of the peptide in it)` for each
/// protein of `p`. Decoys take their target's position; `None` if the
/// sequence is not found (e.g. the protein is missing from the FASTA).
pub(super) fn protein_starts(ctx: &ExportContext, p: &Peptide) -> Vec<(String, Option<usize>)> {
    let target = target_sequence(p);
    p.proteins
        .iter()
        .map(|acc| {
            let start = ctx.proteins.get(acc.as_ref()).and_then(|seq| {
                seq.as_bytes()
                    .windows(target.len())
                    .position(|w| w == target.as_slice())
                    .map(|i| i + 1)
            });
            let name = if p.decoy && ctx.db.generate_decoys {
                format!("{}{}", ctx.db.decoy_tag, acc)
            } else {
                acc.to_string()
            };
            (name, start)
        })
        .collect()
}

pub(super) fn mod_mass(ctx: &ExportContext, id: Option<ModId>) -> Option<f32> {
    id.map(|id| ctx.db.modifications.get(id).mass)
}

/// Mass shift as xiVIEW accepts it after a residue: `(15.9949)`, `(-17.0265)`
fn xi_mod(mass: f32) -> String {
    format!("({:.4})", mass)
}

/// Peptide sequence for xiVIEW: uppercase residues, each modification written
/// after its residue; terminal modifications go with the first / last residue
fn xi_sequence(ctx: &ExportContext, p: &Peptide) -> String {
    let mut out = String::new();
    let last = p.sequence.len().saturating_sub(1);
    for (i, (&aa, &m)) in p.sequence.iter().zip(p.modifications.iter()).enumerate() {
        out.push(aa as char);
        if i == 0 {
            if let Some(mass) = mod_mass(ctx, p.nterm) {
                out.push_str(&xi_mod(mass));
            }
        }
        if let Some(mass) = mod_mass(ctx, m) {
            out.push_str(&xi_mod(mass));
        }
        if i == last {
            if let Some(mass) = mod_mass(ctx, p.cterm) {
                out.push_str(&xi_mod(mass));
            }
        }
    }
    out
}

fn tolerance_text(t: Tolerance) -> String {
    match t {
        Tolerance::Ppm(lo, hi) => format!("{} ppm", lo.abs().max(hi.abs())),
        Tolerance::Da(lo, hi) => format!("{} Da", lo.abs().max(hi.abs())),
        Tolerance::Pct(lo, hi) => format!("{} ppm", 1e4 * lo.abs().max(hi.abs())),
    }
}

pub(super) fn mz(mass: f32, charge: u8) -> f32 {
    mass / charge.max(1) as f32 + PROTON
}

/// `crosslinks_xiView.csv`. Loop-links are left out: xiVIEW's CSV cannot
/// describe them. Returns the file and the number of loop-links skipped.
pub fn xiview_csv(ctx: &ExportContext, matches: &[CrosslinkMatch]) -> anyhow::Result<(Vec<u8>, usize)> {
    let mut wtr = csv::WriterBuilder::new().from_writer(vec![]);
    wtr.write_record([
        "PepSeq1", "PepPos1", "LinkPos1", "Protein1", "Decoy1",
        "PepSeq2", "PepPos2", "LinkPos2", "Protein2", "Decoy2",
        "Score", "Charge", "ExpMZ", "CalcMZ", "CrossLinkerModMass",
        "PassThreshold", "Rank", "PeakListFileName", "ScanId",
        "FragmentTolerance", "IonTypes", "MetaCrosslinkQ", "MetaSpectrumQ",
    ])?;
    let ion_types = std::iter::once("peptide".to_string())
        .chain(ctx.ion_kinds.iter().map(|k| match k {
            Kind::ZRadical => "z".to_string(),
            k => format!("{:?}", k).to_ascii_lowercase(),
        }))
        .collect::<Vec<_>>()
        .join(";");
    let tolerance = tolerance_text(ctx.fragment_tol);
    let linker = ctx.settings.linker_mass();
    let mut loops = 0;
    for m in matches {
        let Some(beta) = m.beta.filter(|_| m.kind == LinkKind::Inter) else {
            loops += 1;
            continue;
        };
        let side = |p: &Peptide| {
            let starts = protein_starts(ctx, p);
            let names = starts.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>().join(";");
            let positions = starts
                .iter()
                .map(|(_, s)| s.map_or("-1".to_string(), |s| s.to_string()))
                .collect::<Vec<_>>()
                .join(";");
            let decoy = if p.decoy { "True" } else { "False" };
            (xi_sequence(ctx, p), positions, names, decoy)
        };
        let (a, b) = (&ctx.db[m.alpha], &ctx.db[beta]);
        let (seq1, pos1, prot1, decoy1) = side(a);
        let (seq2, pos2, prot2, decoy2) = side(b);
        let (peak_list, scan_id) = match ctx.files.get(m.file_id) {
            Some(file) => (file.name.clone(), file.xiview_scan_id(&m.spec_id)),
            None => (String::new(), m.spec_id.clone()),
        };
        wtr.write_record([
            seq1,
            pos1,
            (m.alpha_site + 1).to_string(),
            prot1,
            decoy1.to_string(),
            seq2,
            pos2,
            (m.beta_site + 1).to_string(),
            prot2,
            decoy2.to_string(),
            m.hyperscore.to_string(),
            m.charge.to_string(),
            mz(m.expmass, m.charge).to_string(),
            mz(m.calcmass, m.charge).to_string(),
            linker.to_string(),
            (if m.crosslink_q <= Q_CUTOFF { "True" } else { "False" }).to_string(),
            m.rank.to_string(),
            peak_list,
            scan_id,
            tolerance.clone(),
            ion_types.clone(),
            m.crosslink_q.to_string(),
            m.spectrum_q.to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok((wtr.into_inner()?, loops))
}

// ---------------------------------------------------------------------------
// ProXL

pub(super) fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c if (c as u32) < 0x20 && !matches!(c, '\t' | '\n' | '\r') => {}
            c => out.push(c),
        }
    }
    out
}

/// An `xs:decimal` (no exponent, no NaN): Rust's float Display never uses an exponent
pub(super) fn decimal(x: f64) -> String {
    if x.is_finite() { format!("{}", x) } else { "0".into() }
}

/// Sage keeps masses as f32: write their shortest form (15.9949, not 15.99489974975586)
pub(super) fn decimal32(x: f32) -> String {
    if x.is_finite() { format!("{}", x) } else { "0".into() }
}

/// The scan number in an mzML spectrum id, read the way ProXL's mzML reader
/// (JRAP `MLScanAndHeaderParser`) does: the value of `scan=` or `scanId=` up to
/// the next space, or the whole id if it is a number
fn id_scan_number(spec_id: &str) -> Option<u64> {
    let mut s = spec_id;
    if let Some(i) = s.find("scan=") {
        s = &s[i + "scan=".len()..];
    }
    if let Some(i) = s.find("scanId=") {
        s = &s[i + "scanId=".len()..];
    }
    if let Some(i) = s.find(' ') {
        s = &s[..i];
    }
    s.parse::<u64>().ok()
}

/// ProXL's scan number of a spectrum: from its id, otherwise (like JRAP) its
/// position in the file + 1
fn scan_number(spec_id: &str, position: Option<usize>) -> Option<u64> {
    id_scan_number(spec_id)
        .or_else(|| position.map(|i| i as u64 + 1))
        .filter(|&n| n > 0)
}

/// Name, chemical formula and spacer arm (Å) ProXL knows reagents by
fn proxl_linker(settings: &CrosslinkSettings) -> (&'static str, Option<&'static str>, Option<f32>) {
    match settings.reagent {
        Some(Reagent::Dss) => ("dss", Some("C8H10O2"), Some(11.4)),
        Some(Reagent::Bs3) => ("bs3", Some("C8H10O2"), Some(11.4)),
        Some(Reagent::Dsg) => ("dsg", Some("C5H4O2"), Some(7.7)),
        Some(Reagent::Bs2g) => ("bs2g", Some("C5H4O2"), Some(7.7)),
        Some(Reagent::Dsso) => ("dsso", Some("C6H6O3S"), Some(10.1)),
        Some(Reagent::Dsbu) => ("dsbu", Some("C9H14N2O3"), Some(12.5)),
        Some(Reagent::Edc) => ("edc", None, Some(0.0)),
        Some(Reagent::Dmtmm) => ("dmtmm", None, Some(0.0)),
        Some(Reagent::Clixlink) => ("clixlink", Some("C12H13NO4S"), Some(9.2)),
        None if settings.mass.is_some() => ("custom", None, None),
        None => ("dss", Some("C8H10O2"), Some(11.4)),
    }
}

/// A `<linked_end>`: residues plus `[`/`]` protein termini, Sage's site notation
fn linked_end(spec: &str) -> String {
    let mut x = String::from("<linked_end>");
    let residues: Vec<char> = spec.chars().filter(|c| c.is_ascii_uppercase()).collect();
    if !residues.is_empty() {
        x.push_str("<residues>");
        for r in residues {
            let _ = write!(x, "<residue>{}</residue>", r);
        }
        x.push_str("</residues>");
    }
    if spec.contains('[') || spec.contains(']') {
        x.push_str("<protein_termini>");
        if spec.contains('[') {
            // the initiator methionine is often removed: allow positions 0 and 1
            x.push_str(r#"<protein_terminus terminus_end="n" distance_from_terminus="0"/>"#);
            x.push_str(r#"<protein_terminus terminus_end="n" distance_from_terminus="1"/>"#);
        }
        if spec.contains(']') {
            x.push_str(r#"<protein_terminus terminus_end="c" distance_from_terminus="0"/>"#);
        }
        x.push_str("</protein_termini>");
    }
    x.push_str("</linked_end>");
    x
}

/// One peptide of a reported peptide: key text, sequence and `<modifications>`
struct ProxlPeptide {
    key: String,
    sequence: String,
    mods: String,
    link: Vec<u16>,
}

fn proxl_peptide(ctx: &ExportContext, p: &Peptide, links: Vec<u16>) -> ProxlPeptide {
    let is_static = |aa: u8, id: ModId| {
        ctx.static_mods
            .iter()
            .any(|&(spec, sid)| spec == ModificationSpecificity::Residue(aa) && sid == id)
    };
    let mut mods = String::new();
    if let Some(m) = mod_mass(ctx, p.nterm) {
        let _ = write!(mods, r#"<modification mass="{}" is_n_terminal="true"/>"#, decimal32(m));
    }
    for (i, (&aa, &m)) in p.sequence.iter().zip(p.modifications.iter()).enumerate() {
        // static residue modifications are declared once in <static_modifications>
        if let Some(id) = m.filter(|&id| !is_static(aa, id)) {
            let _ = write!(
                mods,
                r#"<modification mass="{}" position="{}"/>"#,
                decimal32(ctx.db.modifications.get(id).mass),
                i + 1
            );
        }
    }
    if let Some(m) = mod_mass(ctx, p.cterm) {
        let _ = write!(mods, r#"<modification mass="{}" is_c_terminal="true"/>"#, decimal32(m));
    }
    let positions = links.iter().map(|l| (l + 1).to_string()).collect::<Vec<_>>().join(",");
    ProxlPeptide {
        key: format!("{}({})", p.display(&ctx.db.modifications), positions),
        sequence: String::from_utf8_lossy(&p.sequence).into_owned(),
        mods,
        link: links,
    }
}

const PROGRAM: &str = "Sage";

/// `crosslinks_proXL.xml`: rank-1 target-target matches, grouped by reported
/// peptide, with Sage's q-values and scores as filterable annotations.
/// `None` when there is nothing to report (ProXL needs at least one peptide).
pub fn proxl_xml(ctx: &ExportContext, matches: &[CrosslinkMatch]) -> Option<String> {
    // reported peptide string → (type, peptides, psms)
    let mut reported: BTreeMap<String, (&str, Vec<ProxlPeptide>, Vec<&CrosslinkMatch>)> = BTreeMap::new();
    let mut accessions = BTreeSet::new();
    for m in matches.iter().filter(|m| m.rank == 1 && m.class == PairClass::TT) {
        let alpha = &ctx.db[m.alpha];
        let (kind, mut peptides) = match (m.kind, m.beta) {
            (LinkKind::Inter, Some(beta)) => (
                "crosslink",
                vec![
                    proxl_peptide(ctx, alpha, vec![m.alpha_site]),
                    proxl_peptide(ctx, &ctx.db[beta], vec![m.beta_site]),
                ],
            ),
            _ => {
                let mut sites = vec![m.alpha_site, m.beta_site];
                sites.sort_unstable();
                ("looplink", vec![proxl_peptide(ctx, alpha, sites)])
            }
        };
        // the same pair in either order is one reported peptide
        peptides.sort_by(|a, b| a.key.cmp(&b.key));
        let key = peptides.iter().map(|p| p.key.as_str()).collect::<Vec<_>>().join("-");
        for ix in std::iter::once(m.alpha).chain(m.beta) {
            accessions.extend(ctx.db[ix].proteins.iter().map(|a| a.to_string()));
        }
        reported.entry(key).or_insert((kind, peptides, Vec::new())).2.push(m);
    }
    if reported.is_empty() {
        return None;
    }

    let mut x = String::new();
    x.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(x, r#"<proxl_input fasta_filename="{}">"#, xml_escape(ctx.fasta_name));

    // search program and its annotation types
    let filterable = [
        ("crosslink q-value", "q-value of the crosslinked peptide pair (Sage target-decoy FDR)", "below", Some(Q_CUTOFF)),
        ("PSM q-value", "q-value of the crosslink-spectrum match", "below", None),
        ("hyperscore", "Sage hyperscore over both peptides", "above", None),
        ("delta next", "hyperscore difference to the next best candidate", "above", None),
    ];
    let descriptive = [
        ("precursor ppm", "precursor mass error (ppm)"),
        ("matched ions A", "fragment ions matched on the first peptide"),
        ("matched ions B", "fragment ions matched on the second peptide"),
        ("stub peaks", "MS-cleavable linkers: stub doublet peaks supporting the pair"),
    ];
    x.push_str("<search_program_info><search_programs>");
    let _ = write!(
        x,
        r#"<search_program name="{p}" display_name="{p}" version="{}" description="Sage crosslink search">"#,
        xml_escape(ctx.version),
        p = PROGRAM
    );
    x.push_str("<psm_annotation_types><filterable_psm_annotation_types>");
    for (name, desc, dir, default) in filterable {
        let _ = write!(
            x,
            r#"<filterable_psm_annotation_type name="{name}" description="{desc}" filter_direction="{dir}""#
        );
        if let Some(v) = default {
            let _ = write!(x, r#" default_filter="true" default_filter_value="{}""#, decimal32(v));
        } else {
            x.push_str(r#" default_filter="false""#);
        }
        x.push_str("/>");
    }
    x.push_str("</filterable_psm_annotation_types><descriptive_psm_annotation_types>");
    for (name, desc) in descriptive {
        let _ = write!(x, r#"<descriptive_psm_annotation_type name="{name}" description="{desc}"/>"#);
    }
    x.push_str("</descriptive_psm_annotation_types></psm_annotation_types></search_program></search_programs>");
    x.push_str("<default_visible_annotations><visible_psm_annotations>");
    for name in ["crosslink q-value", "PSM q-value", "hyperscore", "precursor ppm"] {
        let _ = write!(x, r#"<search_annotation search_program="{PROGRAM}" annotation_name="{name}"/>"#);
    }
    x.push_str("</visible_psm_annotations></default_visible_annotations>");
    x.push_str("<annotation_sort_order><psm_annotation_sort_order>");
    for name in ["crosslink q-value", "hyperscore"] {
        let _ = write!(x, r#"<search_annotation search_program="{PROGRAM}" annotation_name="{name}"/>"#);
    }
    x.push_str("</psm_annotation_sort_order></annotation_sort_order></search_program_info>\n");

    // linker
    let (name, formula, spacer) = proxl_linker(ctx.settings);
    let _ = write!(x, r#"<linkers><linker name="{name}""#);
    if let Some(s) = spacer {
        let _ = write!(x, r#" spacer_arm_length="{}""#, decimal32(s));
    }
    x.push_str("><crosslink_masses>");
    let _ = write!(x, r#"<crosslink_mass mass="{}""#, decimal32(ctx.settings.linker_mass()));
    if let Some(f) = formula {
        let _ = write!(x, r#" chemical_formula="{f}""#);
    }
    x.push_str("/>");
    if let Some(stubs) = ctx.settings.stubs() {
        // every stub of every arm (asymmetric linkers such as cliXlink have two sets)
        for m in stubs.arms().iter().flat_map(|arm| arm.all()) {
            let _ = write!(x, r#"<cleaved_crosslink_mass mass="{}"/>"#, decimal32(m));
        }
    }
    x.push_str("</crosslink_masses><linked_ends>");
    let (site1, site2) = ctx.settings.hetero_sites().unwrap_or_else(|| {
        let mut s = ctx.settings.residues.clone();
        if ctx.settings.protein_nterm {
            s.push('[');
        }
        (s.clone(), s)
    });
    x.push_str(&linked_end(&site1));
    x.push_str(&linked_end(&site2));
    x.push_str("</linked_ends></linker></linkers>\n");

    // reported peptides and their PSMs
    x.push_str("<reported_peptides>\n");
    for (key, (kind, peptides, psms)) in &reported {
        let _ = write!(
            x,
            r#"<reported_peptide reported_peptide_string="{}" type="{kind}"><peptides>"#,
            xml_escape(key)
        );
        for p in peptides {
            let _ = write!(x, r#"<peptide sequence="{}">"#, p.sequence);
            if !p.mods.is_empty() {
                let _ = write!(x, "<modifications>{}</modifications>", p.mods);
            }
            x.push_str("<linked_positions>");
            for l in &p.link {
                let _ = write!(x, r#"<linked_position position="{}"/>"#, l + 1);
            }
            x.push_str("</linked_positions></peptide>");
        }
        x.push_str("</peptides><psms>");
        for m in psms {
            x.push_str("<psm");
            if let Some((name, scan)) = ctx.files.get(m.file_id).and_then(|f| f.proxl_scan(&m.spec_id)) {
                let _ = write!(x, r#" scan_file_name="{}""#, xml_escape(name));
                if let Some(n) = scan {
                    let _ = write!(x, r#" scan_number="{n}""#);
                }
            }
            let _ = write!(
                x,
                r#" precursor_charge="{}" precursor_retention_time="{}" precursor_m_z="{}" linker_mass="{}">"#,
                m.charge.max(1),
                decimal32(m.rt * 60.0),
                decimal32(mz(m.expmass, m.charge)),
                decimal32(ctx.settings.linker_mass()),
            );
            x.push_str("<filterable_psm_annotations>");
            for (name, value) in [
                ("crosslink q-value", decimal32(m.crosslink_q)),
                ("PSM q-value", decimal32(m.spectrum_q)),
                ("hyperscore", decimal(m.hyperscore)),
                ("delta next", decimal(m.delta_next)),
            ] {
                let _ = write!(
                    x,
                    r#"<filterable_psm_annotation search_program="{PROGRAM}" annotation_name="{name}" value="{value}"/>"#
                );
            }
            x.push_str("</filterable_psm_annotations><descriptive_psm_annotations>");
            for (name, value) in [
                ("precursor ppm", decimal32(m.precursor_ppm)),
                ("matched ions A", m.alpha_matched.to_string()),
                ("matched ions B", m.beta_matched.to_string()),
                ("stub peaks", m.stub_peaks.to_string()),
            ] {
                let _ = write!(
                    x,
                    r#"<descriptive_psm_annotation search_program="{PROGRAM}" annotation_name="{name}" value="{value}"/>"#
                );
            }
            x.push_str("</descriptive_psm_annotations></psm>");
        }
        x.push_str("</psms></reported_peptide>\n");
    }
    x.push_str("</reported_peptides>\n");

    // proteins the reported peptides were matched to
    let found: Vec<_> = accessions.iter().filter_map(|a| ctx.proteins.get(a).map(|s| (a, s))).collect();
    if !found.is_empty() {
        x.push_str("<matched_proteins>\n");
        for (acc, seq) in found {
            let _ = writeln!(
                x,
                r#"<protein sequence="{}"><protein_annotation name="{}"/></protein>"#,
                xml_escape(seq),
                xml_escape(acc)
            );
        }
        x.push_str("</matched_proteins>\n");
    }

    let statics: Vec<_> = ctx
        .static_mods
        .iter()
        .filter_map(|&(spec, id)| match spec {
            ModificationSpecificity::Residue(aa) => Some((aa as char, ctx.db.modifications.get(id).mass)),
            _ => None,
        })
        .collect();
    if !statics.is_empty() {
        x.push_str("<static_modifications>");
        for (aa, mass) in statics {
            let _ = write!(x, r#"<static_modification amino_acid="{aa}" mass_change="{}"/>"#, decimal32(mass));
        }
        x.push_str("</static_modifications>\n");
    }
    if ctx.db.generate_decoys {
        let _ = writeln!(
            x,
            r#"<decoy_labels><decoy_label prefix="{}"/></decoy_labels>"#,
            xml_escape(&ctx.db.decoy_tag)
        );
    }
    x.push_str("</proxl_input>\n");
    Some(x)
}

#[cfg(test)]
pub(super) mod test {
    use super::*;
    use sage_core::database::{Builder, PeptideIx};
    use sage_core::fasta::Fasta;

    const FASTA: &str = ">sp|P11111|ONE_HUMAN first protein\n\
        MAGLEKVSTAAMQRCTDYHFPNEWKLIRSSPEKAQLR\n\
        >sp|P22222|TWO_HUMAN second protein\n\
        GGAKCDKEFRQLMNRSTKPLEK\n";

    /// A real Sage database (static carbamidomethyl C, variable oxidised M)
    #[allow(clippy::type_complexity)]
    pub(in crate::runner) fn fixture() -> (
        IndexedDatabase,
        Vec<(ModificationSpecificity, ModId)>,
        Vec<(ModificationSpecificity, ModId, Option<usize>)>,
        HashMap<String, String>,
    ) {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "test.fasta",
            "enzyme": {"missed_cleavages": 2, "min_len": 4, "max_len": 30},
            "peptide_min_mass": 300.0,
            "static_mods": {"C": 57.021464},
            "variable_mods": {"M": [15.9949]},
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let statics = params.modifications.static_mods.clone();
        let variables = params.modifications.variable_mods.clone();
        let fasta = Fasta::parse(FASTA.into(), "rev_", true);
        let proteins = fasta.targets.iter().map(|(a, s)| (a.to_string(), s.clone())).collect();
        (params.build(fasta), statics, variables, proteins)
    }

    /// Index of a peptide by its modified sequence as Sage displays it
    pub(in crate::runner) fn find(db: &IndexedDatabase, shown: &str, decoy: bool) -> PeptideIx {
        let ix = db
            .peptides
            .iter()
            .position(|p| p.decoy == decoy && p.display(&db.modifications).to_string() == shown)
            .unwrap_or_else(|| {
                let have: Vec<String> = db
                    .peptides
                    .iter()
                    .map(|p| format!("{}{}", p.display(&db.modifications), if p.decoy { " (decoy)" } else { "" }))
                    .collect();
                panic!("{shown} (decoy {decoy}) not in the database: {have:?}")
            });
        PeptideIx(ix as u32)
    }

    pub(in crate::runner) fn csm(alpha: PeptideIx, beta: Option<PeptideIx>, sites: (u16, u16), class: PairClass, rank: u32, q: f32) -> CrosslinkMatch {
        CrosslinkMatch {
            psm_id: 1,
            spec_id: "controllerType=0 controllerNumber=1 scan=4321".into(),
            file_id: 0,
            rank,
            kind: if beta.is_some() { LinkKind::Inter } else { LinkKind::Loop },
            class,
            alpha,
            beta,
            alpha_site: sites.0,
            beta_site: sites.1,
            intra_protein: false,
            charge: 4,
            rt: 12.5,
            expmass: 2500.0,
            calcmass: 2500.001,
            precursor_ppm: 0.4,
            isotope_error: 0.0,
            hyperscore: 42.5,
            alpha_score: 30.0,
            beta_score: 20.0,
            alpha_matched: 9,
            beta_matched: 6,
            alpha_linked: 4,
            beta_linked: 3,
            delta_next: 3.25,
            pairs_scored: 100,
            stub_peaks: 0,
            ims: 0.0,
            spectrum_q: 0.001,
            crosslink_q: q,
        }
    }

    /// Writes both exports; with `SAGE_EXPORT_DIR` set they are also saved there
    /// for checking with the real tools (ProXL XSD, xiVIEW CSV reader)
    #[test]
    fn exports_real_matches() {
        let (db, statics, variables, proteins) = fixture();
        let settings = CrosslinkSettings::default();
        // K6 of MAGLEK... linked to K3 of CDKEFR (whose decoy is CFEKDR)
        let a = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        let b = find(&db, "C[+57.021465]DKEFR", false);
        let decoy = find(&db, "C[+57.021465]FEKDR", true);
        let loop_pep = find(&db, "STKPLEK", false);
        let matches = vec![
            csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.001),
            csm(a, Some(decoy), (5, 3), PairClass::TD, 1, 0.2),
            csm(loop_pep, None, (2, 6), PairClass::TT, 1, 0.004),
            // rank 2 matches are kept for xiVIEW, not for ProXL
            csm(b, Some(a), (2, 5), PairClass::TT, 2, 0.5),
        ];
        let files = vec![SpectrumFile::from_name("run01.mzML.gz")];
        let ctx = ExportContext {
            db: &db,
            settings: &settings,
            files: &files,
            proteins: &proteins,
            fasta_name: "test.fasta",
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            ion_kinds: &[Kind::B, Kind::Y],
            static_mods: &statics,
            variable_mods: &variables,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            version: "test",
        };

        // xiVIEW: crosslinks only, decoys flagged, positions in the protein
        let (csv, loops) = xiview_csv(&ctx, &matches).unwrap();
        let text = String::from_utf8(csv.clone()).unwrap();
        let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split(',').collect()).collect();
        assert_eq!(loops, 1);
        assert_eq!(rows.len(), 4, "{text}");
        assert_eq!(rows[0][..4], ["PepSeq1", "PepPos1", "LinkPos1", "Protein1"]);
        // MAGLEK... starts at residue 1 of P11111, CDKEFR at residue 5 of P22222
        assert_eq!(rows[1][..6], ["MAGLEKVSTAAM(15.9949)QR", "1", "6", "sp|P11111|ONE_HUMAN", "False", "C(57.0215)DKEFR"]);
        assert_eq!(rows[1][6..10], ["5", "3", "sp|P22222|TWO_HUMAN", "False"]);
        assert_eq!(rows[1][15], "True");
        // the decoy keeps its target's position under the decoy protein name
        assert_eq!(rows[2][5..10], ["C(57.0215)FEKDR", "5", "4", "rev_sp|P22222|TWO_HUMAN", "True"]);
        assert_eq!(rows[2][15], "False");
        assert_eq!(rows[3][16], "2");
        assert!(text.contains("controllerType=0 controllerNumber=1 scan=4321"));
        // the uncompressed name: what xiVIEW and ProXL read
        assert!(text.contains(",run01.mzML,") && !text.contains(".gz"));

        // ProXL: rank-1 target matches only
        let xml = proxl_xml(&ctx, &matches).unwrap();
        assert_eq!(xml.matches("<reported_peptide ").count(), 2, "{xml}");
        assert_eq!(xml.matches("<psm ").count(), 2);
        assert!(xml.contains(r#"type="crosslink""#) && xml.contains(r#"type="looplink""#));
        assert!(xml.contains(r#"<linked_position position="3"/><linked_position position="7"/>"#));
        assert!(xml.contains(r#"<modification mass="15.9949" position="12"/>"#));
        // carbamidomethyl is static: declared once, not on each peptide
        assert!(xml.contains(r#"<static_modification amino_acid="C" mass_change="57.021465"/>"#));
        assert!(!xml.contains(r#"<modification mass="57.021465""#));
        assert!(xml.contains(r#"scan_file_name="run01.mzML" scan_number="4321""#));
        assert!(xml.contains(r#"precursor_retention_time="750""#));
        assert!(xml.contains(r#"<decoy_label prefix="rev_"/>"#));
        assert!(!xml.contains("rev_sp"), "decoys are not reported to ProXL");

        if let Some(dir) = std::env::var_os("SAGE_EXPORT_DIR") {
            let dir = std::path::PathBuf::from(dir);
            std::fs::write(dir.join("crosslinks_xiView.csv"), &csv).unwrap();
            std::fs::write(dir.join("crosslinks_proXL.xml"), &xml).unwrap();
        }
    }

    /// The references written for every spectrum of real files, read through
    /// Sage's own readers, for checking with the viewers' readers (pyteomics for
    /// xiVIEW, ProXL's JRAP numbering). Runs when `SAGE_SPECTRA_DIR` is set.
    #[test]
    fn spectrum_references_of_real_files() {
        let Some(dir) = std::env::var_os("SAGE_SPECTRA_DIR") else { return };
        let dir = std::path::PathBuf::from(dir);
        let mut out = String::from("file\tsage_id\txiview_peaklist\txiview_scanid\tproxl_file\tproxl_scan\n");
        let mut names: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
            .filter(|n| n.contains(".mzML") || n.ends_with(".mgf"))
            .collect();
        names.sort();
        for name in names {
            let url = sage_cloudpath::to_url(dir.join(&name).to_str().unwrap()).unwrap();
            let spectra = if name.ends_with(".mgf") {
                sage_cloudpath::util::read_mgf(&url, 0).unwrap()
            } else {
                sage_cloudpath::util::read_mzml(&url, 0, None).unwrap()
            };
            let mut file = SpectrumFile::from_name(&name);
            if file.needs_positions(spectra.iter().map(|s| s.id.as_str())) {
                file.positions = Some(sage_cloudpath::util::read_spectrum_positions(&url).unwrap());
            }
            for s in &spectra {
                let (pfile, pscan) = file.proxl_scan(&s.id).map_or(("", None), |(f, n)| (f, n));
                out.push_str(&format!(
                    "{name}\t{}\t{}\t{}\t{pfile}\t{}\n",
                    s.id,
                    file.name,
                    file.xiview_scan_id(&s.id),
                    pscan.map_or(String::new(), |n| n.to_string())
                ));
            }
        }
        std::fs::write(dir.join("references.tsv"), out).unwrap();
    }

    #[test]
    fn scan_numbers_follow_proxl() {
        assert_eq!(scan_number("controllerType=0 controllerNumber=1 scan=1234", None), Some(1234));
        assert_eq!(scan_number("frame=12 scan=3", Some(99)), Some(3));
        assert_eq!(scan_number("scanId=17", None), Some(17));
        assert_eq!(scan_number("2442", Some(5)), Some(2442));
        // no scan number in the id: position in the file + 1, as ProXL's mzML reader does
        assert_eq!(scan_number("spectrum=2442", Some(286)), Some(287));
        assert_eq!(scan_number("index=0", Some(0)), Some(1));
        assert_eq!(scan_number("spectrum=2442", None), None);
    }

    #[test]
    fn spectrum_files_as_viewers_see_them() {
        let f = SpectrumFile::from_name("Run_01.mzML.gz");
        assert_eq!((f.name.as_str(), f.kind), ("Run_01.mzML", SpectrumKind::MzML));
        assert_eq!(SpectrumFile::from_name("a.MGF").kind, SpectrumKind::Mgf);
        assert_eq!(SpectrumFile::from_name("sample.d").kind, SpectrumKind::Other);
        // positions are read only where the ids do not give the scan number
        let mzml = SpectrumFile::from_name("a.mzML");
        assert!(!mzml.needs_positions(["controllerType=0 controllerNumber=1 scan=5"].into_iter()));
        assert!(mzml.needs_positions(["spectrum=5"].into_iter()));
        assert!(SpectrumFile::from_name("a.mgf").needs_positions(std::iter::empty()));
    }

    /// MGF spectra go to xiVIEW as `index=N`; ProXL gets no scan reference for them
    #[test]
    fn mgf_spectra_by_position() {
        let (db, statics, variables, proteins) = fixture();
        let settings = CrosslinkSettings::default();
        let a = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        let b = find(&db, "C[+57.021465]DKEFR", false);
        let mut m = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.001);
        m.spec_id = "Elution from: 12.3 to 12.3 period: 0 experiment: 1 cycles: 1".into();
        let mut file = SpectrumFile::from_name("run02.mgf");
        file.positions = Some(HashMap::from([(m.spec_id.clone(), 41)]));
        let files = vec![file];
        let ctx = ExportContext {
            db: &db,
            settings: &settings,
            files: &files,
            proteins: &proteins,
            fasta_name: "test.fasta",
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            ion_kinds: &[Kind::B, Kind::Y],
            static_mods: &statics,
            variable_mods: &variables,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            version: "test",
        };
        let (csv, _) = xiview_csv(&ctx, std::slice::from_ref(&m)).unwrap();
        assert!(String::from_utf8(csv).unwrap().contains(",run02.mgf,index=41,"));
        let xml = proxl_xml(&ctx, &[m]).unwrap();
        assert!(!xml.contains("scan_file_name") && !xml.contains("scan_number"));
    }

    #[test]
    fn xml_values_are_safe() {
        assert_eq!(xml_escape(r#"a<b & "c">'"#), "a&lt;b &amp; &quot;c&quot;&gt;&apos;");
        assert_eq!(decimal(0.00001), "0.00001");
        assert_eq!(decimal(f64::NAN), "0");
        assert_eq!(xi_mod(-17.026549), "(-17.0265)");
    }

    #[test]
    fn linked_ends_follow_sage_site_notation() {
        assert_eq!(
            linked_end("K["),
            "<linked_end><residues><residue>K</residue></residues><protein_termini>\
             <protein_terminus terminus_end=\"n\" distance_from_terminus=\"0\"/>\
             <protein_terminus terminus_end=\"n\" distance_from_terminus=\"1\"/>\
             </protein_termini></linked_end>"
        );
        assert_eq!(linked_end("DE]").matches("<residue>").count(), 2);
    }
}
