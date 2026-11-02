//! `crosslinks.mzid`: crosslink results as mzIdentML 1.3.0, following the PSI
//! "mzIdentML extension for crosslinking data" 1.0.0 (HUPO-PSI/mzIdentML,
//! `specification_document/specdoc1_3/asciidoc/crosslinking_ext.adoc`).
//!
//! This is the format the current xiVIEW (including a local installation, via
//! Rappsilber-Laboratory/mzidentml-reader) and PRIDE take for crosslinking data.
//!
//! * Each crosslinked species is a pair of `<Peptide>` elements joined by
//!   "crosslink donor" / "crosslink acceptor" (MS:1002509 / MS:1002510) with a
//!   shared value; the donor carries the linker mass and the reagent term, the
//!   acceptor a zero mass and no reagent term.
//! * A crosslink-spectrum match is two `<SpectrumIdentificationItem>`s sharing a
//!   "crosslink spectrum identification item" (MS:1002511) value; a loop-link is
//!   one item marked "looplink spectrum identification item" (MS:1003329) on a
//!   peptide carrying both ends.
//! * Decoys are included (`isDecoy`) and `passThreshold` is set at crosslink
//!   q ≤ 0.01, so viewers can show the FDR.

use super::crosslink_export::{
    decimal, decimal32, mz, protein_starts, target_sequence, xml_escape, ExportContext,
    SpectrumKind, Q_CUTOFF,
};
use sage_core::crosslink::{CrosslinkMatch, LinkKind, Reagent, SITE1, SITE2};
use sage_core::ion_series::Kind;
use sage_core::mass::Tolerance;
use sage_core::modification::{ModId, Modification, ModificationSpecificity};
use sage_core::peptide::Peptide;
use std::collections::HashMap;
use std::fmt::Write;

fn cv(cv_ref: &str, accession: &str, name: &str, value: Option<&str>) -> String {
    match value {
        Some(v) => format!(
            r#"<cvParam cvRef="{cv_ref}" accession="{accession}" name="{}" value="{}"/>"#,
            xml_escape(name),
            xml_escape(v)
        ),
        None => format!(r#"<cvParam cvRef="{cv_ref}" accession="{accession}" name="{}"/>"#, xml_escape(name)),
    }
}

fn ms(accession: &str, name: &str, value: Option<&str>) -> String {
    cv("PSI-MS", accession, name, value)
}

fn user(name: &str, value: &str) -> String {
    format!(r#"<userParam name="{}" value="{}"/>"#, xml_escape(name), xml_escape(value))
}

/// The reagent's term: XLMOD where it has one, UNIMOD for EDC (as the
/// crosslinking specification's own EDC example), "unknown modification" for
/// a custom mass
fn reagent_cv(ctx: &ExportContext) -> String {
    let (cv_ref, acc, name) = match (ctx.settings.reagent, ctx.settings.mass) {
        (Some(Reagent::Dss), _) | (None, None) => ("XLMOD", "XLMOD:02001", "DSS"),
        (Some(Reagent::Bs3), _) => ("XLMOD", "XLMOD:02000", "BS3"),
        (Some(Reagent::Dsg), _) => ("XLMOD", "XLMOD:02006", "DSG"),
        (Some(Reagent::Bs2g), _) => ("XLMOD", "XLMOD:02005", "BS2G"),
        (Some(Reagent::Dsso), _) => ("XLMOD", "XLMOD:02126", "DSSO"),
        (Some(Reagent::Dsbu), _) => ("XLMOD", "XLMOD:02120", "BuUrBu"),
        (Some(Reagent::Dmtmm), _) => ("XLMOD", "XLMOD:02208", "DMTMM"),
        (Some(Reagent::Clixlink), _) => ("XLMOD", "XLMOD:02253", "cliXLink"),
        (Some(Reagent::Edc), _) => ("UNIMOD", "UNIMOD:2018", "Xlink:EDC"),
        (None, Some(_)) => return ms("MS:1001460", "unknown modification", Some("custom crosslinker")),
    };
    cv(cv_ref, acc, name, None)
}

/// "crosslinker cleavage characteristics" (MS:1003390), `name:mass:pairs with`.
/// Short stub `A`, long stub `S`, further stubs `T`, `U`, ...: for DSSO this
/// is exactly the specification's encoding (A:54.0105647:ST, S:103.9932001:A,
/// T:85.9826354:A). A short stub pairs with the long and further stubs left on
/// the other peptide; for an asymmetric linker (cliXlink) those are the other
/// arm's, whose stubs are named `B`, `X`, `Y`, ...
fn cleavage_cvs(ctx: &ExportContext) -> String {
    let Some(stubs) = ctx.settings.stubs() else { return String::new() };
    let arms = stubs.arms();
    let letters: [&[u8]; 2] = [b"ASTUVW", b"BXYZQR"];
    let named: Vec<Vec<(char, f32)>> = arms
        .iter()
        .zip(letters)
        .map(|(arm, l)| arm.all().into_iter().zip(l.iter().map(|&c| c as char)).map(|(m, c)| (c, m)).collect())
        .collect();
    let mut out = String::new();
    for (i, own) in named.iter().enumerate() {
        // the stubs a cleavage leaves on the other peptide
        let partner = &named[if named.len() == 2 { 1 - i } else { i }];
        for (k, &(c, m)) in own.iter().enumerate() {
            let with: String = if k == 0 {
                partner[1..].iter().map(|(c, _)| *c).collect()
            } else {
                partner[0].0.to_string()
            };
            out += &ms("MS:1003390", "crosslinker cleavage characteristics", Some(&format!("{c}:{m}:{with}")));
        }
    }
    out
}

/// Common UNIMOD entries, for modifications named `UNIMOD:n` in the search settings
const UNIMOD_NAMES: &[(u32, &str)] = &[
    (1, "Acetyl"),
    (4, "Carbamidomethyl"),
    (7, "Deamidated"),
    (21, "Phospho"),
    (24, "Propionamide"),
    (27, "Glu->pyro-Glu"),
    (28, "Gln->pyro-Glu"),
    (34, "Methyl"),
    (35, "Oxidation"),
    (36, "Dimethyl"),
    (121, "GG"),
    (737, "TMT6plex"),
    (2016, "TMTpro"),
];

/// Term for an ordinary modification: its UNIMOD entry when named after one we
/// know, otherwise "unknown modification" with the name or mass as value
fn modification_cv(m: &Modification) -> String {
    let name = m.name.as_ref().map(|n| n.as_str());
    let unimod = name
        .and_then(|n| n.strip_prefix("UNIMOD:"))
        .and_then(|n| n.parse::<u32>().ok())
        .and_then(|n| UNIMOD_NAMES.iter().find(|(id, _)| *id == n));
    match unimod {
        Some((id, label)) => cv("UNIMOD", &format!("UNIMOD:{id}"), label, None),
        None => ms("MS:1001460", "unknown modification", Some(name.unwrap_or(&decimal32(m.mass)))),
    }
}

/// `residues` and specificity rule of a modification site
fn specificity(spec: ModificationSpecificity) -> (String, Option<(&'static str, &'static str)>) {
    let residue = |r: Option<u8>| r.map_or(".".to_string(), |r| (r as char).to_string());
    match spec {
        ModificationSpecificity::Residue(r) => ((r as char).to_string(), None),
        ModificationSpecificity::PeptideN(r) => (residue(r), Some(("MS:1001189", "modification specificity peptide N-term"))),
        ModificationSpecificity::PeptideC(r) => (residue(r), Some(("MS:1001190", "modification specificity peptide C-term"))),
        ModificationSpecificity::ProteinN(r) => (residue(r), Some(("MS:1002057", "modification specificity protein N-term"))),
        ModificationSpecificity::ProteinC(r) => (residue(r), Some(("MS:1002058", "modification specificity protein C-term"))),
    }
}

fn search_modification(fixed: bool, mass: f32, residues: &str, rule: Option<(&str, &str)>, params: &str) -> String {
    let mut x = format!(
        r#"<SearchModification fixedMod="{fixed}" massDelta="{}" residues="{residues}">"#,
        decimal32(mass)
    );
    if let Some((acc, name)) = rule {
        let _ = write!(x, "<SpecificityRules>{}</SpecificityRules>", ms(acc, name, None));
    }
    let _ = write!(x, "{params}</SearchModification>");
    x
}

/// `<SearchModification>`s of one linker end (Sage site notation: residue
/// letters, `[` protein N-terminus, `]` protein C-terminus)
fn linker_end(ctx: &ExportContext, donor: bool, sites: &str) -> String {
    let (role, mass, term) = if donor {
        ("donor", ctx.settings.linker_mass(), ("MS:1002509", "crosslink donor"))
    } else {
        ("acceptor", 0.0, ("MS:1002510", "crosslink acceptor"))
    };
    let params = |id: &str| {
        format!(
            "{}{}{}{}",
            ms("MS:1003392", "search modification id", Some(id)),
            reagent_cv(ctx),
            ms(term.0, term.1, Some("0")),
            cleavage_cvs(ctx)
        )
    };
    let residues: Vec<String> = sites.chars().filter(|c| c.is_ascii_uppercase()).map(String::from).collect();
    let mut x = String::new();
    if !residues.is_empty() {
        x += &search_modification(false, mass, &residues.join(" "), None, &params(&format!("xl_{role}")));
    }
    if sites.contains('[') {
        x += &search_modification(
            false,
            mass,
            ".",
            Some(("MS:1002057", "modification specificity protein N-term")),
            &params(&format!("xl_{role}_protein_n_term")),
        );
    }
    if sites.contains(']') {
        x += &search_modification(
            false,
            mass,
            ".",
            Some(("MS:1002058", "modification specificity protein C-term")),
            &params(&format!("xl_{role}_protein_c_term")),
        );
    }
    x
}

/// A `<FragmentTolerance>` / `<ParentTolerance>`. Plus and minus are written
/// equal (the wider side): xiVIEW rejects unequal values.
fn tolerance(tag: &str, t: Tolerance) -> String {
    let (value, unit_acc, unit) = match t {
        Tolerance::Ppm(lo, hi) => (lo.abs().max(hi.abs()), "UO:0000169", "parts per million"),
        Tolerance::Pct(lo, hi) => (1e4 * lo.abs().max(hi.abs()), "UO:0000169", "parts per million"),
        Tolerance::Da(lo, hi) => (lo.abs().max(hi.abs()), "UO:0000221", "dalton"),
    };
    let side = |acc: &str, name: &str| {
        format!(
            r#"<cvParam cvRef="PSI-MS" accession="{acc}" name="{name}" value="{}" unitCvRef="UO" unitAccession="{unit_acc}" unitName="{unit}"/>"#,
            decimal32(value)
        )
    };
    format!(
        "<{tag}>{}{}</{tag}>",
        side("MS:1001412", "search tolerance plus value"),
        side("MS:1001413", "search tolerance minus value")
    )
}

/// `xsd:anyURI` for a plain file name: percent-encode everything unusual
fn uri(name: &str) -> String {
    name.bytes()
        .map(|b| match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'.' | b'_' | b'~' => (b as char).to_string(),
            b => format!("%{b:02X}"),
        })
        .collect()
}

/// Current UTC time as `xsd:dateTime`
fn now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // civil-from-days (H. Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z", rem / 3600, rem / 60 % 60, rem % 60)
}

/// What the first (donor) and second (acceptor) linker ends react with, in
/// Sage site notation
fn end_sites(ctx: &ExportContext) -> (String, String) {
    ctx.settings.hetero_sites().unwrap_or_else(|| {
        let mut s = ctx.settings.residues.clone();
        if ctx.settings.protein_nterm {
            s.push('[');
        }
        (s.clone(), s)
    })
}

/// The end of the linker attached at `site` of `p`, for heterobifunctional
/// reagents: does it react as the first end (the donor)?
fn reacts_as_first_end(ctx: &ExportContext, p: &Peptide, site: u16) -> Option<bool> {
    ctx.settings.is_heterobifunctional().then(|| {
        ctx.settings
            .link_site_roles(p)
            .iter()
            .find(|(s, _)| *s == site as usize)
            .is_some_and(|(_, roles)| roles & SITE1 != 0 && roles & SITE2 == 0)
    })
}

/// Peptide elements written so far, keyed by their content
struct Peptides {
    xml: String,
    ids: HashMap<String, String>,
    /// peptide element id → its `<PeptideEvidence>` ids
    evidence: HashMap<String, Vec<String>>,
    evidence_xml: String,
    /// protein name → (DBSequence id, sequence if a target)
    proteins: Vec<(String, Option<String>)>,
    protein_ids: HashMap<String, usize>,
}

/// One end of the link on a peptide: (0-based residue, donor?, pair id)
type LinkEnd = (u16, bool, usize);

impl Peptides {
    fn new() -> Self {
        Peptides {
            xml: String::new(),
            ids: HashMap::new(),
            evidence: HashMap::new(),
            evidence_xml: String::new(),
            proteins: Vec::new(),
            protein_ids: HashMap::new(),
        }
    }

    /// The `<Peptide>` for `p` carrying the given link ends (and its evidence)
    fn get(&mut self, ctx: &ExportContext, p: &Peptide, ends: &[LinkEnd]) -> String {
        let key = format!("{}|{:?}", p.display(&ctx.db.modifications), ends);
        if let Some(id) = self.ids.get(&key) {
            return id.clone();
        }
        let id = format!("PEP_{}", self.ids.len() + 1);
        let len = p.sequence.len();
        let _ = write!(
            self.xml,
            r#"<Peptide id="{id}"><PeptideSequence>{}</PeptideSequence>"#,
            String::from_utf8_lossy(&p.sequence)
        );
        let modification = |location: usize, residue: Option<u8>, id: ModId| {
            let m = ctx.db.modifications.get(id);
            let residues = residue.map_or(String::new(), |r| format!(r#" residues="{}""#, r as char));
            format!(
                r#"<Modification location="{location}"{residues} monoisotopicMassDelta="{}">{}</Modification>"#,
                decimal32(m.mass),
                modification_cv(m)
            )
        };
        if let Some(n) = p.nterm {
            self.xml += &modification(0, None, n);
        }
        for (i, (&aa, &m)) in p.sequence.iter().zip(p.modifications.iter()).enumerate() {
            if let Some(m) = m {
                self.xml += &modification(i + 1, Some(aa), m);
            }
            for &(site, donor, pair) in ends.iter().filter(|e| e.0 as usize == i) {
                let (role, mass, term) = if donor {
                    ("donor", ctx.settings.linker_mass(), ("MS:1002509", "crosslink donor"))
                } else {
                    ("acceptor", 0.0, ("MS:1002510", "crosslink acceptor"))
                };
                let _ = write!(
                    self.xml,
                    r#"<Modification location="{}" residues="{}" monoisotopicMassDelta="{}">"#,
                    site + 1,
                    aa as char,
                    decimal32(mass)
                );
                // the search modification this end matches: the residue, or a
                // protein terminus the linker also reacts with
                let (first, second) = end_sites(ctx);
                let sites = if donor { first } else { second };
                let target = if sites.contains(aa as char) {
                    format!("xl_{role}")
                } else if i == 0 && sites.contains('[') {
                    format!("xl_{role}_protein_n_term")
                } else if i + 1 == len && sites.contains(']') {
                    format!("xl_{role}_protein_c_term")
                } else {
                    format!("xl_{role}")
                };
                self.xml += &ms("MS:1003393", "search modification id ref", Some(&target));
                // only the donor names the reagent
                if donor {
                    self.xml += &reagent_cv(ctx);
                }
                self.xml += &ms(term.0, term.1, Some(&pair.to_string()));
                self.xml += "</Modification>";
            }
        }
        if let Some(c) = p.cterm {
            self.xml += &modification(len + 1, None, c);
        }
        self.xml += "</Peptide>\n";

        // where the peptide is found
        let target = target_sequence(p);
        let mut refs = Vec::new();
        for (name, start) in protein_starts(ctx, p) {
            let dbs = match self.protein_ids.get(&name) {
                Some(&i) => i,
                None => {
                    let seq = if p.decoy { None } else { ctx.proteins.get(&name).cloned() };
                    self.proteins.push((name.clone(), seq));
                    self.protein_ids.insert(name.clone(), self.proteins.len());
                    self.proteins.len()
                }
            };
            let pe = format!("PE_{}", self.evidence.values().map(Vec::len).sum::<usize>() + refs.len() + 1);
            let _ = write!(self.evidence_xml, r#"<PeptideEvidence id="{pe}" peptide_ref="{id}" dBSequence_ref="DBS_{dbs}""#);
            if let Some(s) = start {
                let _ = write!(self.evidence_xml, r#" start="{s}" end="{}""#, s + len - 1);
                // flanking residues of the target protein
                let acc = name.strip_prefix(&ctx.db.decoy_tag).filter(|_| p.decoy).unwrap_or(&name);
                if let Some(seq) = ctx.proteins.get(acc).map(|s| s.as_bytes()) {
                    let pre = if s >= 2 { seq[s - 2] as char } else { '-' };
                    let post = seq.get(s - 1 + target.len()).map_or('-', |&c| c as char);
                    let _ = write!(self.evidence_xml, r#" pre="{pre}" post="{post}""#);
                }
            }
            let _ = writeln!(self.evidence_xml, r#" isDecoy="{}"/>"#, p.decoy);
            refs.push(pe);
        }
        self.evidence.insert(id.clone(), refs);
        self.ids.insert(key, id.clone());
        id
    }
}

/// `crosslinks.mzid`. `None` when there are no matches: mzIdentML requires at
/// least one spectrum identification.
pub fn mzid_xml(ctx: &ExportContext, matches: &[CrosslinkMatch]) -> Option<String> {
    if matches.is_empty() {
        return None;
    }
    let mut peptides = Peptides::new();
    // crosslinked species → pair id
    let mut pairs: HashMap<String, usize> = HashMap::new();
    // spectrum → identification items, in the order spectra appear
    let mut results: Vec<((usize, String), String)> = Vec::new();
    let mut result_of: HashMap<(usize, String), usize> = HashMap::new();
    let mut local: HashMap<usize, usize> = HashMap::new();

    for m in matches {
        let alpha = &ctx.db[m.alpha];
        let mut species = |a: String, b: String| {
            let n = pairs.len() + 1;
            *pairs.entry(format!("{a}|{b}")).or_insert(n)
        };
        let key = (m.file_id, m.spec_id.clone());
        let r = *result_of.entry(key.clone()).or_insert_with(|| {
            results.push((key.clone(), String::new()));
            results.len() - 1
        });
        let k = {
            let c = local.entry(r).or_insert(0);
            *c += 1;
            *c
        };
        let common = format!(
            r#"rank="{}" chargeState="{}" experimentalMassToCharge="{}" calculatedMassToCharge="{}" passThreshold="{}""#,
            m.rank,
            m.charge,
            decimal32(mz(m.expmass, m.charge)),
            decimal32(mz(m.calcmass, m.charge)),
            m.crosslink_q <= Q_CUTOFF
        );
        let scores = format!(
            "{}{}{}{}",
            ms("MS:1002354", "PSM-level q-value", Some(&decimal32(m.spectrum_q))),
            user("Sage:crosslink q-value", &decimal32(m.crosslink_q)),
            user("Sage:hyperscore", &decimal(m.hyperscore)),
            user("Sage:delta next", &decimal(m.delta_next))
        );
        let item = |id: String, pep: &str, peptides: &Peptides, marker: &str| {
            let refs: String = peptides.evidence[pep]
                .iter()
                .map(|pe| format!(r#"<PeptideEvidenceRef peptideEvidence_ref="{pe}"/>"#))
                .collect();
            format!(
                r#"<SpectrumIdentificationItem id="{id}" peptide_ref="{pep}" {common}>{refs}{marker}{scores}</SpectrumIdentificationItem>"#
            )
        };
        let sii = format!("SII_{}_{k}", r + 1);
        let xml = match (m.kind, m.beta) {
            (LinkKind::Inter, Some(beta_ix)) => {
                let beta = &ctx.db[beta_ix];
                // the donor matches the first linker end; otherwise the longer
                // peptide (then alphabetical), as the specification recommends
                let alpha_donor = match reacts_as_first_end(ctx, alpha, m.alpha_site) {
                    Some(first) => first,
                    None => (alpha.sequence.len(), beta.sequence.as_ref()) > (beta.sequence.len(), alpha.sequence.as_ref()),
                };
                let (d, ds, a, as_) = if alpha_donor {
                    (alpha, m.alpha_site, beta, m.beta_site)
                } else {
                    (beta, m.beta_site, alpha, m.alpha_site)
                };
                let pair = species(
                    format!("{}({})", d.display(&ctx.db.modifications), ds),
                    format!("{}({})", a.display(&ctx.db.modifications), as_),
                );
                let donor = peptides.get(ctx, d, &[(ds, true, pair)]);
                let acceptor = peptides.get(ctx, a, &[(as_, false, pair)]);
                let marker = ms("MS:1002511", "crosslink spectrum identification item", Some(&k.to_string()));
                item(format!("{sii}_d"), &donor, &peptides, &marker)
                    + &item(format!("{sii}_a"), &acceptor, &peptides, &marker)
            }
            _ => {
                let (mut first, mut second) = (m.alpha_site.min(m.beta_site), m.alpha_site.max(m.beta_site));
                if reacts_as_first_end(ctx, alpha, second) == Some(true) {
                    std::mem::swap(&mut first, &mut second);
                }
                let pair = species(format!("{}({first},{second})", alpha.display(&ctx.db.modifications)), String::new());
                let pep = peptides.get(ctx, alpha, &[(first, true, pair), (second, false, pair)]);
                item(sii, &pep, &peptides, &ms("MS:1003329", "looplink spectrum identification item", None))
            }
        };
        results[r].1 += &xml;
    }

    let mut x = String::with_capacity(1 << 20);
    x.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(
        x,
        r#"<MzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psidev.info/psi/pi/mzIdentML/1.3 https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.3.0.xsd" id="sage_crosslinks" version="1.3.0" creationDate="{}">"#,
        now()
    );
    x.push_str(
        r#"<cvList><cv id="PSI-MS" fullName="Proteomics Standards Initiative Mass Spectrometry Vocabularies" uri="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/><cv id="UNIMOD" fullName="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo"/><cv id="UO" fullName="UNIT-ONTOLOGY" uri="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/><cv id="XLMOD" fullName="PSI cross-link modifications" uri="https://raw.githubusercontent.com/HUPO-PSI/xlmod-CV/main/XLMOD.obo"/></cvList>"#,
    );
    x.push('\n');
    x += &ms("MS:1003385", "mzIdentML crosslinking extension document version", Some("1.0.0"));
    let _ = writeln!(
        x,
        r#"
<AnalysisSoftwareList><AnalysisSoftware id="AS_sage" name="Sage" version="{}" uri="https://github.com/lazear/sage"><SoftwareName>{}</SoftwareName></AnalysisSoftware></AnalysisSoftwareList>"#,
        xml_escape(ctx.version),
        ms("MS:1004007", "Sage", None)
    );

    // proteins, peptides and where they are found
    x.push_str("<SequenceCollection>\n");
    for (i, (name, seq)) in peptides.proteins.iter().enumerate() {
        let _ = write!(
            x,
            r#"<DBSequence id="DBS_{}" accession="{}" searchDatabase_ref="SDB_1""#,
            i + 1,
            xml_escape(name)
        );
        match seq {
            Some(s) => {
                let _ = writeln!(x, r#" length="{}"><Seq>{}</Seq></DBSequence>"#, s.len(), xml_escape(s));
            }
            None => x.push_str("/>\n"),
        }
    }
    x += &peptides.xml;
    x += &peptides.evidence_xml;
    x.push_str("</SequenceCollection>\n");

    // what was searched, and how
    x.push_str(r#"<AnalysisCollection><SpectrumIdentification id="SI_1" spectrumIdentificationProtocol_ref="SIP_1" spectrumIdentificationList_ref="SIL_1">"#);
    for i in 0..ctx.files.len().max(1) {
        let _ = write!(x, r#"<InputSpectra spectraData_ref="SD_{}"/>"#, i + 1);
    }
    x.push_str("<SearchDatabaseRef searchDatabase_ref=\"SDB_1\"/></SpectrumIdentification></AnalysisCollection>\n");

    x.push_str(r#"<AnalysisProtocolCollection><SpectrumIdentificationProtocol id="SIP_1" analysisSoftware_ref="AS_sage">"#);
    let _ = write!(x, "<SearchType>{}</SearchType><AdditionalSearchParams>", ms("MS:1001083", "ms-ms search", None));
    x += &ms("MS:1001211", "parent mass type mono", None);
    x += &ms("MS:1001256", "fragment mass type mono", None);
    x += &ms("MS:1002494", "crosslinking search", None);
    for kind in ctx.ion_kinds {
        let (acc, name) = match kind {
            Kind::A => ("MS:1001108", "param: a ion"),
            Kind::B => ("MS:1001118", "param: b ion"),
            Kind::C => ("MS:1001119", "param: c ion"),
            Kind::X => ("MS:1001261", "param: x ion"),
            Kind::Y => ("MS:1001262", "param: y ion"),
            Kind::Z | Kind::ZRadical => ("MS:1001263", "param: z ion"),
        };
        x += &ms(acc, name, None);
    }
    x.push_str("</AdditionalSearchParams><ModificationParams>");
    // the linker: first end (donor) and second end (acceptor)
    let (site1, site2) = end_sites(ctx);
    x += &linker_end(ctx, true, &site1);
    x += &linker_end(ctx, false, &site2);
    for &(spec, id) in ctx.static_mods {
        let (residues, rule) = specificity(spec);
        let m = ctx.db.modifications.get(id);
        x += &search_modification(true, m.mass, &residues, rule, &modification_cv(m));
    }
    for &(spec, id, _) in ctx.variable_mods {
        let (residues, rule) = specificity(spec);
        let m = ctx.db.modifications.get(id);
        x += &search_modification(false, m.mass, &residues, rule, &modification_cv(m));
    }
    x.push_str("</ModificationParams>");
    x += &tolerance("FragmentTolerance", ctx.fragment_tol);
    x += &tolerance("ParentTolerance", ctx.precursor_tol);
    let _ = write!(x, "<Threshold>{}</Threshold>", user("Sage:crosslink q-value", &decimal32(Q_CUTOFF)));
    x.push_str("</SpectrumIdentificationProtocol></AnalysisProtocolCollection>\n");

    // input files and the identifications
    x.push_str("<DataCollection><Inputs>");
    let _ = write!(
        x,
        r#"<SearchDatabase id="SDB_1" location="{}"><FileFormat>{}</FileFormat><DatabaseName>{}</DatabaseName>{}</SearchDatabase>"#,
        uri(ctx.fasta_name),
        ms("MS:1001348", "FASTA format", None),
        format_args!(r#"<userParam name="{}"/>"#, xml_escape(ctx.fasta_name)),
        ms("MS:1001073", "database type amino acid", None)
    );
    for (i, file) in ctx.files.iter().enumerate() {
        let (format, id_format) = match file.kind {
            SpectrumKind::MzML => (ms("MS:1000584", "mzML format", None), ms("MS:1001530", "mzML unique identifier", None)),
            SpectrumKind::Mgf => (
                ms("MS:1001062", "Mascot MGF format", None),
                ms("MS:1000774", "multiple peak list nativeID format", None),
            ),
            SpectrumKind::Other => (ms("MS:1002817", "Bruker TDF format", None), ms("MS:1000824", "no nativeID format", None)),
        };
        let _ = write!(
            x,
            r#"<SpectraData id="SD_{}" location="{}" name="{}"><FileFormat>{format}</FileFormat><SpectrumIDFormat>{id_format}</SpectrumIDFormat></SpectraData>"#,
            i + 1,
            uri(&file.name),
            xml_escape(&file.name)
        );
    }
    x.push_str("</Inputs>\n<AnalysisData><SpectrumIdentificationList id=\"SIL_1\">\n");
    for (r, ((file_id, spec_id), items)) in results.iter().enumerate() {
        let file = ctx.files.get(*file_id);
        let spectrum_id = file.map_or_else(|| spec_id.clone(), |f| f.xiview_scan_id(spec_id));
        let rt = matches
            .iter()
            .find(|m| m.file_id == *file_id && &m.spec_id == spec_id)
            .map_or(0.0, |m| m.rt);
        let _ = writeln!(
            x,
            r#"<SpectrumIdentificationResult id="SIR_{}" spectrumID="{}" spectraData_ref="SD_{}">{items}<cvParam cvRef="PSI-MS" accession="MS:1000016" name="scan start time" value="{}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/></SpectrumIdentificationResult>"#,
            r + 1,
            xml_escape(&spectrum_id),
            file_id + 1,
            decimal32(rt)
        );
    }
    x.push_str("</SpectrumIdentificationList></AnalysisData></DataCollection>\n</MzIdentML>\n");
    Some(x)
}

#[cfg(test)]
mod test {
    use super::super::crosslink_export::test::{csm, find, fixture};
    use super::super::crosslink_export::SpectrumFile;
    use super::*;
    use sage_core::crosslink::{CrosslinkSettings, PairClass};

    #[test]
    fn dates_are_xsd_datetimes() {
        let d = now();
        assert_eq!(d.len(), 20, "{d}");
        assert!(d.starts_with("20") && d.ends_with('Z') && &d[10..11] == "T");
    }

    #[test]
    fn cleavable_stubs_as_in_the_specification() {
        let settings: CrosslinkSettings = serde_json::from_str(r#"{"reagent": "dsso"}"#).unwrap();
        let (db, statics, variables, proteins) = fixture();
        let files = vec![SpectrumFile::from_name("a.mzML")];
        let ctx = context(&db, &settings, &files, &proteins, &statics, &variables);
        let x = cleavage_cvs(&ctx);
        for v in ["A:54.01056:ST", "S:103.9932:A", "T:85.98264:A"] {
            assert!(x.contains(&format!(r#"value="{v}""#)), "{x}");
        }
    }

    fn context<'a>(
        db: &'a sage_core::database::IndexedDatabase,
        settings: &'a CrosslinkSettings,
        files: &'a [SpectrumFile],
        proteins: &'a HashMap<String, String>,
        statics: &'a [(ModificationSpecificity, ModId)],
        variables: &'a [(ModificationSpecificity, ModId, Option<usize>)],
    ) -> ExportContext<'a> {
        ExportContext {
            db,
            settings,
            files,
            proteins,
            fasta_name: "test proteins.fasta",
            fragment_tol: Tolerance::Ppm(-10.0, 15.0),
            ion_kinds: &[Kind::B, Kind::Y],
            static_mods: statics,
            variable_mods: variables,
            precursor_tol: Tolerance::Ppm(-5.0, 5.0),
            version: "test",
        }
    }

    /// Cleavable (DSSO, stubs) and zero-length heterobifunctional (EDC: amine
    /// to carboxyl, protein termini) linkers; saved with `SAGE_EXPORT_DIR` too
    #[test]
    fn writes_other_reagents() {
        let (db, statics, variables, proteins) = fixture();
        let a = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        let b = find(&db, "C[+57.021465]DKEFR", false);
        let files = vec![SpectrumFile::from_name("run.mzML")];
        for (name, json) in [
            ("dsso", r#"{"reagent": "dsso"}"#),
            ("edc", r#"{"reagent": "edc"}"#),
            ("clixlink", r#"{"reagent": "clixlink"}"#),
        ] {
            let settings: CrosslinkSettings = serde_json::from_str(json).unwrap();
            let ctx = context(&db, &settings, &files, &proteins, &statics, &variables);
            // EDC: K of the first peptide (amine end) to D of CDKEFR (carboxyl end)
            let m = csm(a, Some(b), (5, 1), PairClass::TT, 1, 0.001);
            let x = mzid_xml(&ctx, &[m]).unwrap();
            if name == "dsso" {
                assert!(x.contains(r#"accession="XLMOD:02126" name="DSSO""#));
                assert!(x.contains(r#"value="A:54.01056:ST""#));
            } else if name == "clixlink" {
                assert!(x.contains(r#"accession="XLMOD:02253" name="cliXLink""#));
                // 267.05653 as a 32-bit float
                assert!(x.contains(r#"monoisotopicMassDelta="267.05652""#));
                // each arm's alkene pairs with the other arm's sulfenic acid / thiol
                for v in ["A:54.01056:XY", "S:103.9932:B", "T:85.98264:B", "B:163.06332:ST", "X:213.04596:A", "Y:195.0354:A"] {
                    assert!(x.contains(&format!(r#"value="{v}""#)), "{v}");
                }
            } else {
                assert!(x.contains(r#"accession="UNIMOD:2018" name="Xlink:EDC""#));
                assert!(x.contains(r#"monoisotopicMassDelta="-18.010565""#));
                // the amine end (K) is the donor; carboxyl end D E + protein C-terminus
                assert!(x.contains(r#"massDelta="-18.010565" residues="K">"#));
                assert!(x.contains(r#"massDelta="0" residues="D E">"#));
                assert!(x.contains("xl_acceptor_protein_c_term"));
                assert!(x.contains(r#"<Modification location="6" residues="K" monoisotopicMassDelta="-18.010565">"#));
                assert!(x.contains(r#"<Modification location="2" residues="D" monoisotopicMassDelta="0">"#));
            }
            if let Some(dir) = std::env::var_os("SAGE_EXPORT_DIR") {
                std::fs::write(std::path::PathBuf::from(dir).join(format!("crosslinks_{name}.mzid")), &x).unwrap();
            }
        }
    }

    /// Crosslinks, a decoy pair, a loop-link and a rank-2 match; with
    /// `SAGE_EXPORT_DIR` set the file is saved for the schema and reader checks
    #[test]
    fn writes_crosslinks_loops_and_decoys() {
        let (db, statics, variables, proteins) = fixture();
        let settings = CrosslinkSettings::default();
        let a = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        let b = find(&db, "C[+57.021465]DKEFR", false);
        let decoy = find(&db, "C[+57.021465]FEKDR", true);
        let loop_pep = find(&db, "STKPLEK", false);
        let mut matches = vec![
            csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.001),
            csm(a, Some(decoy), (5, 3), PairClass::TD, 1, 0.2),
            csm(loop_pep, None, (2, 6), PairClass::TT, 1, 0.004),
            csm(b, Some(a), (2, 5), PairClass::TT, 2, 0.5),
        ];
        // two spectra: the first holds rank 1 and 2 of one precursor
        matches[1].spec_id = "controllerType=0 controllerNumber=1 scan=4400".into();
        matches[2].spec_id = "controllerType=0 controllerNumber=1 scan=4500".into();
        let files = vec![SpectrumFile::from_name("run 01.mzML.gz")];
        let ctx = context(&db, &settings, &files, &proteins, &statics, &variables);
        let x = mzid_xml(&ctx, &matches).unwrap();
        if let Some(dir) = std::env::var_os("SAGE_EXPORT_DIR") {
            std::fs::write(std::path::PathBuf::from(dir).join("crosslinks.mzid"), &x).unwrap();
        }

        // the same species (rank 1 and rank 2 in either order) shares one pair of peptides
        assert_eq!(x.matches("<Peptide ").count(), 5, "{x}");
        // search modifications carry the reagent-local value 0, the peptides their pair
        assert_eq!(x.matches(r#"name="crosslink donor" value="0""#).count(), 2, "residues + protein N-term");
        assert_eq!(x.matches(r#"name="crosslink donor" value="1""#).count(), 1);
        assert_eq!(x.matches(r#"name="crosslink acceptor" value="1""#).count(), 1);
        // donor: the longer peptide, with the linker mass and the reagent
        assert!(x.contains(r#"<Peptide id="PEP_1"><PeptideSequence>MAGLEKVSTAAMQR</PeptideSequence>"#));
        assert!(x.contains(r#"<Modification location="12" residues="M" monoisotopicMassDelta="15.9949">"#));
        assert!(x.contains(r#"<Modification location="6" residues="K" monoisotopicMassDelta="138.06808">"#));
        assert!(x.contains(r#"<Modification location="3" residues="K" monoisotopicMassDelta="0">"#));
        // reagent on both ends of the search modifications (residues + protein N-term), and on donor peptides only
        assert_eq!(x.matches(r#"accession="XLMOD:02001""#).count(), 4 + 3);
        // loop-link: one item, both ends on one peptide
        assert_eq!(x.matches("looplink spectrum identification item").count(), 1);
        // a crosslink is two items with the same local id and identical values
        assert!(x.contains(r#"<SpectrumIdentificationItem id="SII_1_1_d""#) && x.contains(r#"<SpectrumIdentificationItem id="SII_1_1_a""#));
        assert_eq!(x.matches(r#"name="crosslink spectrum identification item" value="2""#).count(), 2);
        // decoys flagged, positions from the FASTA, flanking residues
        assert!(x.contains(r#"start="1" end="14" pre="-" post="C" isDecoy="false""#), "{x}");
        assert!(x.contains(r#"isDecoy="true""#) && x.contains(r#"accession="rev_sp|P22222|TWO_HUMAN""#));
        // spectra: the uncompressed name, percent-encoded as a URI
        assert!(x.contains(r#"location="run%2001.mzML" name="run 01.mzML""#));
        assert!(x.contains(r#"spectrumID="controllerType=0 controllerNumber=1 scan=4321""#));
        // plus and minus tolerance equal (xiVIEW requires it)
        assert!(x.contains(r#"name="search tolerance plus value" value="15""#) && x.contains(r#"name="search tolerance minus value" value="15""#));
        assert!(x.contains(r#"<cvParam cvRef="UNIMOD" accession="UNIMOD:4"#) || x.contains("unknown modification"));

    }
}
