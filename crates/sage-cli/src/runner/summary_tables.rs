//! Peptide- and protein-level tables of a linear search, the usual reports of a
//! bottom-up analysis next to the PSM table `results.sage.tsv`:
//!
//! - `peptides.sage.tsv`: one row per target peptide (modified sequence) at 1%
//!   peptide FDR, with its PSM count per file, best scores and q-values and, with
//!   LFQ, its MS1 intensity per file (all charge states summed).
//! - `proteins.sage.tsv`: one row per target protein group with a unique
//!   peptide in `peptides.sage.tsv`, with its q-value (Sage's picked group FDR on
//!   proteotypic peptides), unique and shared peptides, PSM counts per file and,
//!   with LFQ, its intensity per file (sum of its unique peptides).
//!
//! Protein groups are not cut at 1% group FDR: the picked FDR counts one decoy
//! from the start, so its q-values cannot fall below 1 / (number of target
//! groups). A purified complex of a dozen proteins, the usual crosslinking
//! sample, would never pass. Filter on `protein_group_q` for large databases.
//!
//! Only rank-1 PSMs count. A peptide shared between protein groups counts for
//! each of them as a shared peptide; group q-values come from unique peptides
//! alone, so a group seen only through shared peptides is not listed.

use super::Runner;
use sage_cloudpath::Url;
use sage_core::database::{IndexedDatabase, PeptideIx};
use sage_core::lfq::{Peak, PrecursorId};
use sage_core::scoring::Feature;
use std::collections::{BTreeSet, HashMap};

/// Peptides at or below this q-value are listed (and protein groups counted as passing)
const Q_CUTOFF: f32 = 0.01;

/// LFQ areas as `sage_core::lfq` reports them: per precursor (and decoy flag),
/// the peak and its intensity in each file
pub type Areas = HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>;

/// One row of `peptides.sage.tsv`
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideRow {
    /// Modified sequence, as in `results.sage.tsv`
    pub peptide: String,
    /// Plain amino-acid sequence
    pub sequence: String,
    pub proteins: String,
    pub protein_groups: String,
    pub num_proteins: usize,
    pub num_protein_groups: u32,
    pub missed_cleavages: u8,
    pub calcmass: f32,
    pub charges: BTreeSet<u8>,
    /// PSMs in each file
    pub psms: Vec<usize>,
    /// Best (lowest) spectrum q-value of its PSMs
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,
    pub protein_group_q: f32,
    /// Best discriminant score of its PSMs, and the hyperscore and retention
    /// time of that PSM
    pub score: f32,
    pub hyperscore: f64,
    pub rt: f32,
    /// Summed MS2 intensity of its PSMs
    pub ms2_intensity: f32,
    /// LFQ: MS1 intensity in each file
    pub intensity: Option<Vec<f64>>,
}

/// One row of `proteins.sage.tsv`
#[derive(Clone, Debug, PartialEq)]
pub struct ProteinRow {
    /// Protein group, as in `results.sage.tsv` (members joined by `/`)
    pub protein_group: String,
    pub num_proteins: usize,
    pub protein_group_q: f32,
    pub unique_peptides: usize,
    pub shared_peptides: usize,
    /// PSMs of all its peptides (unique and shared) in each file
    pub psms: Vec<usize>,
    /// Best discriminant score of its unique peptides
    pub score: f32,
    /// LFQ: summed intensity of its unique peptides in each file
    pub intensity: Option<Vec<f64>>,
}

/// Peptides at 1% FDR from the rank-1 target PSMs, best first
pub fn peptide_table(
    db: &IndexedDatabase,
    features: &[Feature],
    files: usize,
    areas: Option<&Areas>,
) -> Vec<PeptideRow> {
    // LFQ intensity per peptide and file, charge states summed
    let lfq = areas.map(|areas| {
        let mut sums: HashMap<PeptideIx, Vec<f64>> = HashMap::new();
        for ((id, decoy), (_, data)) in areas {
            if *decoy {
                continue;
            }
            let ix = match id {
                PrecursorId::Combined(ix) | PrecursorId::Charged((ix, _)) => *ix,
            };
            let sum = sums.entry(ix).or_insert_with(|| vec![0.0; files]);
            for (s, x) in sum.iter_mut().zip(data) {
                *s += x;
            }
        }
        sums
    });

    let mut rows: HashMap<String, PeptideRow> = HashMap::new();
    let mut ixs: HashMap<String, BTreeSet<PeptideIx>> = HashMap::new();
    for f in features
        .iter()
        .filter(|f| f.label == 1 && f.rank == 1 && f.peptide_q <= Q_CUTOFF)
    {
        let p = &db[f.peptide_idx];
        let peptide = p.display(&db.modifications).to_string();
        ixs.entry(peptide.clone()).or_default().insert(f.peptide_idx);
        let row = rows.entry(peptide.clone()).or_insert_with(|| PeptideRow {
            peptide,
            sequence: String::from_utf8_lossy(&p.sequence).into_owned(),
            proteins: p.proteins(&db.decoy_tag, db.generate_decoys),
            protein_groups: f.protein_groups.clone().unwrap_or_default(),
            num_proteins: p.proteins.len(),
            num_protein_groups: f.num_protein_groups,
            missed_cleavages: p.missed_cleavages,
            calcmass: f.calcmass,
            charges: BTreeSet::new(),
            psms: vec![0; files],
            spectrum_q: f.spectrum_q,
            peptide_q: f.peptide_q,
            protein_q: f.protein_q,
            protein_group_q: f.protein_group_q,
            score: f.discriminant_score,
            hyperscore: f.hyperscore,
            rt: f.rt,
            ms2_intensity: 0.0,
            intensity: None,
        });
        row.charges.insert(f.charge);
        if let Some(n) = row.psms.get_mut(f.file_id) {
            *n += 1;
        }
        row.spectrum_q = row.spectrum_q.min(f.spectrum_q);
        row.peptide_q = row.peptide_q.min(f.peptide_q);
        row.protein_q = row.protein_q.min(f.protein_q);
        row.protein_group_q = row.protein_group_q.min(f.protein_group_q);
        row.ms2_intensity += f.ms2_intensity;
        if f.discriminant_score > row.score {
            row.score = f.discriminant_score;
            row.hyperscore = f.hyperscore;
            row.rt = f.rt;
        }
    }

    let mut rows: Vec<PeptideRow> = rows
        .into_values()
        .map(|mut row| {
            if let Some(lfq) = &lfq {
                let mut sum = vec![0.0; files];
                for ix in &ixs[&row.peptide] {
                    for (s, x) in sum.iter_mut().zip(lfq.get(ix).into_iter().flatten()) {
                        *s += x;
                    }
                }
                row.intensity = Some(sum);
            }
            row
        })
        .collect();
    rows.sort_by(|a, b| b.score.total_cmp(&a.score).then_with(|| a.peptide.cmp(&b.peptide)));
    rows
}

/// Protein groups with a unique peptide in the peptide table, best first
pub fn protein_table(peptides: &[PeptideRow], files: usize, lfq: bool) -> Vec<ProteinRow> {
    let mut rows: HashMap<&str, ProteinRow> = HashMap::new();
    for p in peptides.iter().filter(|p| !p.protein_groups.is_empty()) {
        let unique = p.num_protein_groups == 1;
        for group in p.protein_groups.split(';') {
            let row = rows.entry(group).or_insert_with(|| ProteinRow {
                protein_group: group.to_string(),
                num_proteins: group.split('/').count(),
                protein_group_q: 1.0,
                unique_peptides: 0,
                shared_peptides: 0,
                psms: vec![0; files],
                score: f32::MIN,
                intensity: lfq.then(|| vec![0.0; files]),
            });
            for (n, x) in row.psms.iter_mut().zip(&p.psms) {
                *n += x;
            }
            if !unique {
                row.shared_peptides += 1;
                continue;
            }
            row.unique_peptides += 1;
            row.protein_group_q = row.protein_group_q.min(p.protein_group_q);
            row.score = row.score.max(p.score);
            if let (Some(sum), Some(x)) = (&mut row.intensity, &p.intensity) {
                for (s, x) in sum.iter_mut().zip(x) {
                    *s += x;
                }
            }
        }
    }
    let mut rows: Vec<ProteinRow> = rows
        .into_values()
        .filter(|r| r.unique_peptides > 0)
        .collect();
    rows.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.protein_group.cmp(&b.protein_group))
    });
    rows
}

fn float<F: ryu::Float>(x: F) -> String {
    ryu::Buffer::new().format(x).to_string()
}

/// `peptides.sage.tsv`
pub fn peptides_tsv(rows: &[PeptideRow], filenames: &[String]) -> anyhow::Result<Vec<u8>> {
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
    let mut headers: Vec<String> = [
        "peptide",
        "sequence",
        "proteins",
        "protein_groups",
        "num_proteins",
        "num_protein_groups",
        "missed_cleavages",
        "calcmass",
        "charges",
        "psms",
        "spectrum_q",
        "peptide_q",
        "protein_q",
        "protein_group_q",
        "sage_discriminant_score",
        "hyperscore",
        "rt",
        "ms2_intensity",
    ]
    .map(String::from)
    .to_vec();
    headers.extend(filenames.iter().map(|f| format!("psms:{f}")));
    let lfq = rows.first().is_some_and(|r| r.intensity.is_some());
    if lfq {
        headers.extend(filenames.iter().map(|f| format!("intensity:{f}")));
    }
    wtr.write_record(&headers)?;
    for r in rows {
        let charges = r.charges.iter().map(u8::to_string).collect::<Vec<_>>().join(";");
        let mut record = vec![
            r.peptide.clone(),
            r.sequence.clone(),
            r.proteins.clone(),
            r.protein_groups.clone(),
            r.num_proteins.to_string(),
            r.num_protein_groups.to_string(),
            r.missed_cleavages.to_string(),
            float(r.calcmass),
            charges,
            r.psms.iter().sum::<usize>().to_string(),
            float(r.spectrum_q),
            float(r.peptide_q),
            float(r.protein_q),
            float(r.protein_group_q),
            float(r.score),
            float(r.hyperscore),
            float(r.rt),
            float(r.ms2_intensity),
        ];
        record.extend(r.psms.iter().map(usize::to_string));
        if lfq {
            record.extend(r.intensity.iter().flatten().map(|&x| float(x)));
        }
        wtr.write_record(&record)?;
    }
    Ok(wtr.into_inner()?)
}

/// `proteins.sage.tsv`
pub fn proteins_tsv(rows: &[ProteinRow], filenames: &[String]) -> anyhow::Result<Vec<u8>> {
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
    let mut headers: Vec<String> = [
        "protein_group",
        "num_proteins",
        "protein_group_q",
        "peptides",
        "unique_peptides",
        "shared_peptides",
        "psms",
        "sage_discriminant_score",
    ]
    .map(String::from)
    .to_vec();
    headers.extend(filenames.iter().map(|f| format!("psms:{f}")));
    let lfq = rows.first().is_some_and(|r| r.intensity.is_some());
    if lfq {
        headers.extend(filenames.iter().map(|f| format!("intensity:{f}")));
    }
    wtr.write_record(&headers)?;
    for r in rows {
        let mut record = vec![
            r.protein_group.clone(),
            r.num_proteins.to_string(),
            float(r.protein_group_q),
            (r.unique_peptides + r.shared_peptides).to_string(),
            r.unique_peptides.to_string(),
            r.shared_peptides.to_string(),
            r.psms.iter().sum::<usize>().to_string(),
            float(r.score),
        ];
        record.extend(r.psms.iter().map(usize::to_string));
        if lfq {
            record.extend(r.intensity.iter().flatten().map(|&x| float(x)));
        }
        wtr.write_record(&record)?;
    }
    Ok(wtr.into_inner()?)
}

impl Runner {
    /// `peptides.sage.tsv` and `proteins.sage.tsv` of a linear search
    pub(super) fn write_summary_tables(
        &self,
        features: &[Feature],
        filenames: &[String],
        areas: Option<&Areas>,
    ) -> anyhow::Result<[Url; 2]> {
        let peptides = peptide_table(&self.database, features, filenames.len(), areas);
        let proteins = protein_table(&peptides, filenames.len(), areas.is_some());

        let peptides_path = self.make_path("peptides.sage.tsv");
        sage_cloudpath::write_bytes_sync(&peptides_path, peptides_tsv(&peptides, filenames)?)?;
        let proteins_path = self.make_path("proteins.sage.tsv");
        sage_cloudpath::write_bytes_sync(&proteins_path, proteins_tsv(&proteins, filenames)?)?;
        log::info!(
            "wrote peptides.sage.tsv ({} peptides at 1% FDR) and proteins.sage.tsv ({} protein groups, {} at 1% group FDR)",
            peptides.len(),
            proteins.len(),
            proteins.iter().filter(|p| p.protein_group_q <= Q_CUTOFF).count()
        );
        Ok([peptides_path, proteins_path])
    }
}

#[cfg(test)]
mod test {
    use super::super::crosslink_export::test::{find, fixture};
    use super::*;

    fn psm(db: &IndexedDatabase, shown: &str, file_id: usize, score: f32, groups: &str) -> Feature {
        let ix = find(db, shown, false);
        Feature {
            peptide_idx: ix,
            file_id,
            rank: 1,
            label: 1,
            charge: 2,
            calcmass: db[ix].monoisotopic,
            discriminant_score: score,
            hyperscore: score as f64 * 10.0,
            rt: score,
            ms2_intensity: 100.0,
            spectrum_q: 0.001,
            peptide_q: 0.001,
            protein_q: 0.001,
            protein_group_q: if groups.contains(';') { 1.0 } else { 0.002 },
            protein_groups: Some(groups.into()),
            num_protein_groups: groups.split(';').count() as u32,
            ..Default::default()
        }
    }

    #[test]
    fn peptides_collect_their_psms() {
        let (db, ..) = fixture();
        let one = "sp|P11111|ONE_HUMAN";
        let mut features = vec![
            psm(&db, "VSTAAMQR", 0, 2.0, one),
            psm(&db, "VSTAAMQR", 1, 3.0, one),
            psm(&db, "VSTAAM[+15.9949]QR", 0, 1.0, one),
        ];
        features[1].charge = 3;
        // not listed: a decoy, a rank-2 match and a peptide above 1% FDR
        let mut decoy = psm(&db, "VSTAAMQR", 0, 9.0, one);
        decoy.label = -1;
        let mut second = psm(&db, "VSTAAMQR", 0, 9.0, one);
        second.rank = 2;
        let mut weak = psm(&db, "AQLR", 0, 0.5, one);
        weak.peptide_q = 0.2;
        features.extend([decoy, second, weak]);

        let rows = peptide_table(&db, &features, 2, None);
        assert_eq!(rows.len(), 2, "{rows:?}");
        let top = &rows[0];
        assert_eq!(top.peptide, "VSTAAMQR");
        assert_eq!(top.sequence, "VSTAAMQR");
        assert_eq!(top.psms, vec![1, 1]);
        assert_eq!(top.charges, BTreeSet::from([2, 3]));
        assert_eq!((top.score, top.rt), (3.0, 3.0));
        assert_eq!(top.ms2_intensity, 200.0);
        assert_eq!(rows[1].peptide, "VSTAAM[+15.9949]QR");
        assert_eq!(rows[1].sequence, "VSTAAMQR");
    }

    #[test]
    fn protein_groups_count_unique_and_shared_peptides() {
        let (db, ..) = fixture();
        let (one, two) = ("sp|P11111|ONE_HUMAN", "sp|P22222|TWO_HUMAN");
        let shared = format!("{one};{two}");
        let features = vec![
            psm(&db, "VSTAAMQR", 0, 3.0, one),
            psm(&db, "AQLR", 0, 2.0, one),
            psm(&db, "AQLR", 0, 2.5, one),
            psm(&db, "QLMNR", 0, 1.0, &shared),
        ];
        let peptides = peptide_table(&db, &features, 1, None);
        let proteins = protein_table(&peptides, 1, false);
        // TWO_HUMAN has only a shared peptide: no group q-value, not listed
        assert_eq!(proteins.len(), 1, "{proteins:?}");
        let p = &proteins[0];
        assert_eq!(p.protein_group, one);
        assert_eq!((p.unique_peptides, p.shared_peptides), (2, 1));
        assert_eq!(p.psms, vec![4]);
        assert_eq!(p.score, 3.0);
        assert_eq!(p.protein_group_q, 0.002);

        // a group above 1% group FDR (a small database cannot get lower) is listed with its q
        let mut small = psm(&db, "STKPLEK", 0, 0.5, two);
        small.protein_group_q = 0.08;
        let peptides = peptide_table(&db, &[features, vec![small]].concat(), 1, None);
        let proteins = protein_table(&peptides, 1, false);
        assert_eq!(proteins.len(), 2, "{proteins:?}");
        assert_eq!((proteins[1].protein_group.as_str(), proteins[1].protein_group_q), (two, 0.08));
        assert_eq!((proteins[1].unique_peptides, proteins[1].shared_peptides), (1, 1));
    }

    #[test]
    fn lfq_intensities_sum_charges_and_unique_peptides() {
        let (db, ..) = fixture();
        let one = "sp|P11111|ONE_HUMAN";
        let features = vec![psm(&db, "VSTAAMQR", 0, 3.0, one), psm(&db, "AQLR", 1, 2.0, one)];
        let peak = || Peak { rt: 0, spectral_angle: 1.0, score: 1.0, q_value: 0.0 };
        let (a, b) = (find(&db, "VSTAAMQR", false), find(&db, "AQLR", false));
        let mut areas = Areas::default();
        areas.insert((PrecursorId::Charged((a, 2)), false), (peak(), vec![10.0, 20.0]));
        areas.insert((PrecursorId::Charged((a, 3)), false), (peak(), vec![1.0, 2.0]));
        areas.insert((PrecursorId::Charged((a, 2)), true), (peak(), vec![500.0, 500.0]));
        areas.insert((PrecursorId::Combined(b), false), (peak(), vec![0.0, 7.0]));

        let peptides = peptide_table(&db, &features, 2, Some(&areas));
        let vst = peptides.iter().find(|p| p.peptide == "VSTAAMQR").unwrap();
        assert_eq!(vst.intensity, Some(vec![11.0, 22.0]));
        let proteins = protein_table(&peptides, 2, true);
        assert_eq!(proteins[0].intensity, Some(vec![11.0, 29.0]));

        let names = ["a.mzML".to_string(), "b.mzML".to_string()];
        let text = String::from_utf8(proteins_tsv(&proteins, &names).unwrap()).unwrap();
        let mut lines = text.lines();
        assert_eq!(
            lines.next().unwrap(),
            "protein_group\tnum_proteins\tprotein_group_q\tpeptides\tunique_peptides\tshared_peptides\tpsms\t\
             sage_discriminant_score\tpsms:a.mzML\tpsms:b.mzML\tintensity:a.mzML\tintensity:b.mzML"
        );
        assert_eq!(lines.next().unwrap(), format!("{one}\t1\t0.002\t2\t2\t0\t2\t3.0\t1\t1\t11.0\t29.0"));
        let text = String::from_utf8(peptides_tsv(&peptides, &names).unwrap()).unwrap();
        assert!(text.lines().next().unwrap().ends_with("psms:a.mzML\tpsms:b.mzML\tintensity:a.mzML\tintensity:b.mzML"));
        assert_eq!(text.lines().count(), 3);
    }
}
