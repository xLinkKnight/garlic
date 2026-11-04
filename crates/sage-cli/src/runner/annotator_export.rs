//! Identifications for [Annotator](https://github.com/snijderlab/annotator),
//! which draws annotated spectra of linear and crosslinked peptidoforms:
//! `results.annotator.csv` (linear PSMs) and `crosslinks.annotator.csv`
//! (crosslink-spectrum matches), both target matches at 1% spectrum FDR.
//!
//! The files use the "basic CSV" PSM format of Annotator's `mzident` library:
//! `scan_index,sequence,raw_file,z,mode`. `sequence` is ProForma 2.0; a
//! crosslink is `PEPK[XLMOD:02001#XL1]R//PEPK[#XL1]IDE`, a loop-link puts both
//! ends on one peptide. `scan_index` is the 0-based position of the spectrum in
//! its file. Annotator matches `raw_file` against the spectrum files it has
//! open, so open the spectrum files together with the CSV:
//! `annotator run1.mzML run2.mzML crosslinks.annotator.csv`.

use super::Runner;
use sage_cloudpath::Url;
use sage_core::crosslink::{CrosslinkMatch, CrosslinkSettings, LinkKind, PairClass};
use sage_core::database::IndexedDatabase;
use sage_core::peptide::Peptide;
use sage_core::scoring::Feature;
use std::collections::{BTreeSet, HashMap};

/// Matches at or below this spectrum q-value are exported
const Q_CUTOFF: f32 = 0.01;

/// One row of an Annotator CSV
#[derive(Clone, Debug, PartialEq)]
pub struct AnnotatorRow {
    pub scan_index: usize,
    pub sequence: String,
    pub raw_file: String,
    pub charge: u8,
}

/// ProForma of one peptide (modifications as
/// [`sage_core::modification::Modification::proforma_label`]). `extra` adds
/// bracketed text after residues (by 0-based index), e.g. a crosslinker and
/// its `#XL1` label.
pub fn peptide_proforma(db: &IndexedDatabase, p: &Peptide, extra: &[(usize, String)]) -> String {
    let mods = &db.modifications;
    let mut s = String::new();
    if let Some(id) = p.nterm {
        s.push_str(&format!("[{}]-", mods.get(id).proforma_label()));
    }
    for (i, (&aa, m)) in p.sequence.iter().zip(p.modifications.iter()).enumerate() {
        s.push(aa as char);
        if let Some(id) = m {
            s.push_str(&format!("[{}]", mods.get(*id).proforma_label()));
        }
        for (_, text) in extra.iter().filter(|(site, _)| *site == i) {
            s.push_str(&format!("[{text}]"));
        }
    }
    if let Some(id) = p.cterm {
        s.push_str(&format!("-[{}]", mods.get(id).proforma_label()));
    }
    s
}

/// ProForma of a crosslinked match. A protein N-terminal link is written on
/// the first residue (Annotator rejects crosslinks on a terminus); every b-ion
/// contains both, so the annotation is the same.
pub fn crosslink_proforma(db: &IndexedDatabase, settings: &CrosslinkSettings, m: &CrosslinkMatch) -> String {
    let linker = format!("{}#XL1", settings.proforma_linker());
    let alpha = &db[m.alpha];
    // multivalent: every peptide, link `n` as `#XLn` on both of its residues
    if let Some(s) = &m.species {
        let term = settings.proforma_linker();
        return s
            .peptides
            .iter()
            .enumerate()
            .map(|(x, &ix)| {
                let marks: Vec<(usize, String)> = s
                    .links
                    .iter()
                    .enumerate()
                    .flat_map(|(n, &(a, sa, b, sb))| {
                        [(a, sa, format!("{term}#XL{}", n + 1)), (b, sb, format!("#XL{}", n + 1))]
                    })
                    .filter(|(p, _, _)| *p as usize == x)
                    .map(|(_, site, mark)| (site as usize, mark))
                    .collect();
                peptide_proforma(db, &db[ix], &marks)
            })
            .collect::<Vec<_>>()
            .join("//");
    }
    match (m.kind, m.beta) {
        (LinkKind::Inter, Some(beta)) => format!(
            "{}//{}",
            peptide_proforma(db, alpha, &[(m.alpha_site as usize, linker)]),
            peptide_proforma(db, &db[beta], &[(m.beta_site as usize, "#XL1".into())]),
        ),
        _ => peptide_proforma(
            db,
            alpha,
            &[(m.alpha_site as usize, linker), (m.beta_site as usize, "#XL1".into())],
        ),
    }
}

/// The CSV text. Columns are those of `mzident`'s basic CSV format; `mode`
/// (fragmentation model) is left empty, so it is chosen in Annotator.
pub fn annotator_csv(rows: &[AnnotatorRow]) -> anyhow::Result<Vec<u8>> {
    let mut wtr = csv::WriterBuilder::new().from_writer(vec![]);
    wtr.write_record(["scan_index", "sequence", "raw_file", "z", "mode"])?;
    for r in rows {
        wtr.write_record([
            r.scan_index.to_string(),
            r.sequence.clone(),
            r.raw_file.clone(),
            r.charge.to_string(),
            String::new(),
        ])?;
    }
    wtr.flush()?;
    Ok(wtr.into_inner()?)
}

/// Bruker timsTOF `.d` (DDA-PASEF): Sage's spectrum id (the `Precursors.Id` of
/// `analysis.tdf`) → the 0-based position of that precursor's first PASEF
/// window, as Annotator's reader (mzdata) numbers the spectra.
///
/// mzdata walks the frames in order, skipping frames without peaks; an MS1
/// frame is one spectrum, a DDA-PASEF frame one spectrum per row of
/// `PasefFrameMsMsInfo` (in primary key order: frame, then first scan). Sage
/// merges all frames of a precursor into one spectrum, so a precursor
/// fragmented in several frames shows Annotator only its first frame's.
pub fn bruker_positions(tdf: &std::path::Path) -> anyhow::Result<HashMap<String, usize>> {
    let conn = rusqlite::Connection::open_with_flags(tdf, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let precursors: usize = conn.query_row("SELECT COUNT(*) FROM Precursors", [], |r| r.get(0))?;
    let mut windows: HashMap<i64, Vec<i64>> = HashMap::new();
    let mut q = conn.prepare("SELECT Frame, Precursor FROM PasefFrameMsMsInfo ORDER BY Frame, ScanNumBegin")?;
    for row in q.query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, Option<i64>>(1)?)))? {
        let (frame, precursor) = row?;
        windows.entry(frame).or_default().push(precursor.unwrap_or(0));
    }
    let mut positions = HashMap::new();
    let mut position = 0usize;
    let mut q = conn.prepare("SELECT Id, MsMsType, NumPeaks FROM Frames ORDER BY Id")?;
    for row in q.query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, i64>(1)?, r.get::<_, i64>(2)?)))? {
        let (frame, msms_type, peaks) = row?;
        if peaks == 0 {
            continue;
        }
        match msms_type {
            // MS1
            0 => position += 1,
            // DDA-PASEF: mzdata keeps a window if `precursors[Precursor - 1]`
            // exists (saturating at 0)
            8 => {
                for &precursor in windows.get(&frame).into_iter().flatten() {
                    if ((precursor - 1).max(0) as usize) < precursors {
                        if precursor >= 1 {
                            positions.entry(precursor.to_string()).or_insert(position);
                        }
                        position += 1;
                    }
                }
            }
            // DIA-PASEF and others: Sage's ids are not precursors there
            _ => anyhow::bail!("frame {frame} is not MS1 or DDA-PASEF (MsMsType {msms_type})"),
        }
    }
    Ok(positions)
}

/// Where Annotator finds one spectrum file's spectra
struct SpectrumFile {
    /// Local path (or URL) as passed to Annotator
    path: String,
    /// Sage spectrum id → 0-based position in the file
    positions: HashMap<String, usize>,
}

impl Runner {
    /// Path and spectrum positions of each spectrum file in `file_ids`
    /// (mzML, MGF, and local Bruker `.d` DDA-PASEF folders)
    fn annotator_files(&self, file_ids: &BTreeSet<usize>) -> HashMap<usize, SpectrumFile> {
        let mut files = HashMap::new();
        for &file_id in file_ids {
            let Some(url) = self.parameters.mzml_paths.get(file_id) else {
                continue;
            };
            let local = (url.scheme() == "file").then(|| url.to_file_path().ok()).flatten();
            let path = local
                .as_ref()
                .map_or_else(|| url.to_string(), |p| p.display().to_string());
            let lower = url.path().trim_end_matches('/').to_ascii_lowercase();
            let positions = if lower.ends_with(".d") {
                match &local {
                    Some(dir) => bruker_positions(&dir.join("analysis.tdf")),
                    None => Err(anyhow::anyhow!("Bruker .d folders must be local")),
                }
            } else if [".mzml", ".mzml.gz", ".mgf", ".mgf.gz"].iter().any(|ext| lower.ends_with(ext)) {
                sage_cloudpath::util::read_spectrum_positions(url).map_err(anyhow::Error::from)
            } else {
                Err(anyhow::anyhow!("Annotator reads mzML, MGF and Bruker .d here"))
            };
            match positions {
                Ok(positions) => {
                    files.insert(file_id, SpectrumFile { path, positions });
                }
                Err(e) => log::warn!("Annotator export: {url} left out ({e})"),
            }
        }
        files
    }

    fn write_annotator(&self, name: &str, rows: Vec<(usize, String, String, u8)>) -> anyhow::Result<Url> {
        let file_ids = rows.iter().map(|r| r.0).collect::<BTreeSet<_>>();
        let files = self.annotator_files(&file_ids);
        let mut out = Vec::with_capacity(rows.len());
        let mut missing = 0usize;
        for (file_id, spec_id, sequence, charge) in rows {
            let Some(file) = files.get(&file_id) else {
                continue;
            };
            match file.positions.get(&spec_id) {
                Some(&scan_index) => out.push(AnnotatorRow {
                    scan_index,
                    sequence,
                    raw_file: file.path.clone(),
                    charge,
                }),
                None => missing += 1,
            }
        }
        if missing > 0 {
            log::warn!("Annotator export: {missing} spectra not found in their files were left out");
        }
        let path = self.make_path(name);
        sage_cloudpath::write_bytes_sync(&path, annotator_csv(&out)?)?;
        log::info!("wrote {name} ({} matches for Annotator)", out.len());
        Ok(path)
    }

    /// `results.annotator.csv`: target PSMs (rank 1) at 1% spectrum FDR
    pub(super) fn write_annotator_linear(&self, features: &[Feature]) -> anyhow::Result<Url> {
        let rows = features
            .iter()
            .filter(|f| f.label == 1 && f.rank == 1 && f.spectrum_q <= Q_CUTOFF)
            .map(|f| {
                let sequence = peptide_proforma(&self.database, &self.database[f.peptide_idx], &[]);
                (f.file_id, f.spec_id.clone(), sequence, f.charge)
            })
            .collect();
        self.write_annotator("results.annotator.csv", rows)
    }

    /// `crosslinks.annotator.csv`: target-target crosslink matches at 1% spectrum FDR
    pub(super) fn write_annotator_crosslinks(
        &self,
        matches: &[CrosslinkMatch],
        settings: &CrosslinkSettings,
    ) -> anyhow::Result<Url> {
        let rows = matches
            .iter()
            .filter(|m| m.class == PairClass::TT && m.spectrum_q <= Q_CUTOFF)
            .map(|m| {
                let sequence = crosslink_proforma(&self.database, settings, m);
                (m.file_id, m.spec_id.clone(), sequence, m.charge)
            })
            .collect();
        self.write_annotator("crosslinks.annotator.csv", rows)
    }
}

#[cfg(test)]
mod test {
    use super::super::crosslink_export::test::{csm, find, fixture};
    use sage_core::crosslink::Reagent;
    use sage_core::modification::Modification;
    use super::*;

    #[test]
    fn modification_labels_keep_accessions_only() {
        let m = |name: Option<&str>, mass: f32| Modification {
            name: name.map(|n| n.parse().unwrap()),
            mass,
        };
        assert_eq!(m(None, 15.9949).proforma_label(), "+15.9949");
        assert_eq!(m(None, -17.0265).proforma_label(), "-17.0265");
        assert_eq!(m(Some("UNIMOD:35"), 15.9949).proforma_label(), "UNIMOD:35");
        assert_eq!(m(Some("UniMod:4"), 57.0215).proforma_label(), "UniMod:4");
        assert_eq!(m(Some("U:Oxidation"), 15.9949).proforma_label(), "U:Oxidation");
        // free text names would not parse in Annotator
        assert_eq!(m(Some("my oxidation"), 15.9949).proforma_label(), "+15.9949");
        assert_eq!(m(Some("UNIMOD:"), 1.0).proforma_label(), "+1");
    }

    #[test]
    fn linker_labels() {
        let with = |reagent: Option<Reagent>, mass: Option<f32>| CrosslinkSettings {
            reagent,
            mass,
            ..Default::default()
        };
        assert_eq!(CrosslinkSettings::default().proforma_linker(), "XLMOD:02001");
        assert_eq!(with(Some(Reagent::Dsso), None).proforma_linker(), "XLMOD:02126");
        assert_eq!(with(Some(Reagent::Edc), None).proforma_linker(), "Formula:H-2O-1");
        assert_eq!(with(Some(Reagent::Clixlink), None).proforma_linker(), "+267.05653");
        assert_eq!(with(None, Some(96.5)).proforma_linker(), "+96.5");
    }

    #[test]
    fn linear_peptides_as_proforma() {
        // without named modifications, ProForma is exactly Sage's own notation
        let (db, ..) = fixture();
        for p in &db.peptides {
            assert_eq!(peptide_proforma(&db, p, &[]), p.display(&db.modifications).to_string());
        }
        let oxidised = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        assert_eq!(peptide_proforma(&db, &db[oxidised], &[]), "MAGLEKVSTAAM[+15.9949]QR");
    }

    #[test]
    fn crosslinks_as_proforma() {
        let (db, ..) = fixture();
        let settings = CrosslinkSettings::default();
        let a = find(&db, "MAGLEKVSTAAMQR", false);
        let b = find(&db, "STKPLEK", false);
        let inter = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.0);
        assert_eq!(
            crosslink_proforma(&db, &settings, &inter),
            "MAGLEK[XLMOD:02001#XL1]VSTAAMQR//STK[#XL1]PLEK"
        );
        let loop_link = csm(b, None, (2, 6), PairClass::TT, 1, 0.0);
        assert_eq!(loop_link.kind, LinkKind::Loop);
        assert_eq!(
            crosslink_proforma(&db, &settings, &loop_link),
            "STK[XLMOD:02001#XL1]PLEK[#XL1]"
        );
        // a link on a modified residue keeps both: M[+15.9949][XLMOD:02001#XL1]
        let oxidised = find(&db, "M[+15.9949]AGLEKVSTAAMQR", false);
        let n_term = csm(oxidised, Some(b), (0, 2), PairClass::TT, 1, 0.0);
        assert!(crosslink_proforma(&db, &settings, &n_term)
            .starts_with("M[+15.9949][XLMOD:02001#XL1]AGLEK"));

        // multivalent a1b1a: two copies of a, each linked to b; link n is #XLn
        let mut chain = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.0);
        chain.kind = LinkKind::Multi;
        chain.species = Some(Box::new(sage_core::crosslink::Species {
            peptides: vec![a, b, a],
            links: vec![(0, 5, 1, 2), (1, 6, 2, 5)],
            ..Default::default()
        }));
        assert_eq!(
            crosslink_proforma(&db, &settings, &chain),
            "MAGLEK[XLMOD:02001#XL1]VSTAAMQR//STK[#XL1]PLEK[XLMOD:02001#XL2]//MAGLEK[#XL2]VSTAAMQR"
        );
    }

    #[test]
    fn csv_layout() {
        let rows = [AnnotatorRow {
            scan_index: 41,
            sequence: "PEPK[XLMOD:02001#XL1]R//AK[#XL1]R".into(),
            raw_file: r"C:\data\run, 1.mzML".into(),
            charge: 4,
        }];
        let text = String::from_utf8(annotator_csv(&rows).unwrap()).unwrap();
        assert_eq!(
            text,
            "scan_index,sequence,raw_file,z,mode\n41,PEPK[XLMOD:02001#XL1]R//AK[#XL1]R,\"C:\\data\\run, 1.mzML\",4,\n"
        );
    }
}
