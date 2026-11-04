//! SILAC ratio tables, written when the database has SILAC channels:
//!
//! - `silac.sage.tsv`: one row per peptide (its light form: labels removed,
//!   other modifications kept) identified in any channel, with its MS1
//!   intensity in every channel and run, and each labelled channel's ratio to
//!   light per run.
//! - `silac_proteins.sage.tsv`: per protein, the median ratio of its unique
//!   peptides (peptides of one protein only) per run, and how many there are.
//! - `crosslinks.silac.tsv`: the same per crosslinked species (peptide pair,
//!   loop-linked peptide or multivalent species; link-site isomers have the
//!   same MS1 signal and share a row).
//!
//! Identifications pass at the LFQ `peptide_q_value` (spectrum and peptide
//! q-values; crosslinks: crosslink q-value), rank 1, targets only. All
//! channels of a species are integrated over one MS1 peak
//! ([`sage_core::silac::quantify`]).

use super::Runner;
use sage_cloudpath::Url;
use sage_core::crosslink::{CrosslinkMatch, CrosslinkSettings, LinkKind, PairClass};
use sage_core::database::{IndexedDatabase, PeptideIx};
use sage_core::mass::{composition, Composition};
use sage_core::peptide::Peptide;
use sage_core::scoring::Feature;
use sage_core::silac::{self, Identified, Labels, SilacQuant};
use sage_core::spectrum::ProcessedSpectrum;
use std::collections::{BTreeMap, BTreeSet};

/// What a row names: the species and its proteins
#[derive(Clone, Debug)]
struct Named {
    /// `linear`, `crosslink`, `loop`, `multivalent`
    kind: &'static str,
    /// The light form(s), as in `results.sage.tsv`
    name: String,
    proteins: String,
}

fn envelope<'a>(sequences: impl Iterator<Item = &'a [u8]>, extra: (u16, u16)) -> [f32; 3] {
    let c: Composition = sequences.flat_map(|s| s.iter()).map(|&r| composition(r)).sum();
    sage_core::isotopes::peptide_isotopes(c.carbon + extra.0, c.sulfur + extra.1)
}

/// Mass of the light form in every channel
fn masses(labels: &Labels, lights: &[Peptide], extra: f32) -> Vec<f32> {
    let light = lights.iter().map(|p| p.monoisotopic).sum::<f32>() + extra;
    (0..labels.len() as u8).map(|c| light + lights.iter().map(|p| labels.shift(p, c)).sum::<f32>()).collect()
}

/// The channel name of a peptide (`""` without labelled residues)
pub(super) fn channel_name(db: &IndexedDatabase, ix: PeptideIx) -> String {
    match (&db.silac, db.silac_channel(ix)) {
        (Some(labels), Some(c)) => labels.names[c as usize].clone(),
        _ => String::new(),
    }
}

/// The channel name of a crosslinked species: its peptides' channel, `mixed`
/// if they differ
pub(super) fn species_channel(db: &IndexedDatabase, peptides: &[PeptideIx]) -> String {
    let channels: BTreeSet<String> = peptides.iter().map(|&p| channel_name(db, p)).filter(|c| !c.is_empty()).collect();
    match channels.len() {
        0 => String::new(),
        1 => channels.into_iter().next().unwrap(),
        _ => "mixed".into(),
    }
}

/// Linear PSMs to quantify
fn linear_ids(db: &IndexedDatabase, labels: &Labels, features: &[Feature], q: f32) -> (Vec<Identified>, Vec<Named>) {
    let mut ids = Vec::new();
    let mut named = Vec::new();
    for f in features.iter().filter(|f| f.label == 1 && f.rank == 1 && f.spectrum_q <= q && f.peptide_q <= q) {
        let peptide = &db[f.peptide_idx];
        let Some(channel) = labels.channel(peptide) else { continue };
        let light = labels.unlabel(peptide);
        let name = light.display(&db.modifications).to_string();
        ids.push(Identified {
            key: name.clone(),
            channel,
            file_id: f.file_id,
            rt: f.rt,
            ims: f.ims,
            score: f.hyperscore as f32,
            masses: masses(labels, std::slice::from_ref(&light), 0.0),
            envelope: envelope(std::iter::once(&light.sequence[..]), (0, 0)),
        });
        named.push(Named { kind: "linear", name, proteins: peptide.proteins(&db.decoy_tag, db.generate_decoys) });
    }
    (ids, named)
}

/// Crosslink matches to quantify (a species whose peptides are in different
/// channels is not)
fn crosslink_ids(
    db: &IndexedDatabase,
    labels: &Labels,
    matches: &[CrosslinkMatch],
    settings: &CrosslinkSettings,
    q: f32,
) -> (Vec<Identified>, Vec<Named>) {
    let (linker_c, linker_s) = settings.linker_composition();
    let mut ids = Vec::new();
    let mut named = Vec::new();
    for m in matches.iter().filter(|m| m.class == PairClass::TT && m.rank == 1 && m.crosslink_q <= q) {
        let peptides: Vec<PeptideIx> = match &m.species {
            Some(s) => s.peptides.clone(),
            None => std::iter::once(m.alpha).chain(m.beta).collect(),
        };
        let channels: BTreeSet<u8> = peptides.iter().filter_map(|&p| db.silac_channel(p)).collect();
        let (Some(&channel), 1) = (channels.iter().next(), channels.len()) else { continue };
        let lights: Vec<Peptide> = peptides.iter().map(|&p| labels.unlabel(&db[p])).collect();
        let links = m.species.as_ref().map_or(1, |s| s.links.len()) as u16;
        let mut shown: Vec<String> = lights.iter().map(|p| p.display(&db.modifications).to_string()).collect();
        let kind = match m.kind {
            LinkKind::Inter => "crosslink",
            LinkKind::Loop => "loop",
            LinkKind::Multi => "multivalent",
        };
        let proteins = peptides
            .iter()
            .map(|&p| db[p].proteins(&db.decoy_tag, db.generate_decoys))
            .collect::<Vec<_>>()
            .join(" – ");
        let name = shown.join(" – ");
        shown.sort();
        ids.push(Identified {
            key: format!("{kind}:{}", shown.join("|")),
            channel,
            file_id: m.file_id,
            rt: m.rt,
            ims: m.ims,
            score: m.hyperscore as f32,
            masses: masses(labels, &lights, settings.linker_mass() * links as f32),
            envelope: envelope(lights.iter().map(|p| &p.sequence[..]), (linker_c * links, linker_s * links)),
        });
        named.push(Named { kind, name, proteins });
    }
    (ids, named)
}

/// A float cell, empty when there is none
fn cell(x: Option<f64>) -> String {
    x.map(|x| x.to_string()).unwrap_or_default()
}

/// The species table: who it is, where it was identified, and its
/// intensities and ratios per run and channel
fn species_tsv(
    quants: &[SilacQuant],
    ids: &[Identified],
    named: &[Named],
    names: &[String],
    filenames: &[String],
    crosslinks: bool,
) -> anyhow::Result<Vec<u8>> {
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
    let mut headers: Vec<String> = if crosslinks {
        vec!["link_type".into(), "peptides".into(), "proteins".into(), "identified_in".into(), "csms".into()]
    } else {
        vec!["peptide".into(), "proteins".into(), "identified_in".into(), "psms".into()]
    };
    headers.extend(["peak_score".into(), "spectral_angle".into()]);
    for file in filenames {
        for name in names {
            headers.push(format!("{file} {name}"));
        }
    }
    for file in filenames {
        for name in &names[1..] {
            headers.push(format!("{file} {name}/light"));
        }
    }
    wtr.write_record(&headers)?;
    for q in quants {
        let first = &named[q.identifications[0]];
        let seen: BTreeSet<u8> = q.identifications.iter().map(|&i| ids[i].channel).collect();
        let identified_in = seen.iter().map(|&c| names[c as usize].as_str()).collect::<Vec<_>>().join(";");
        let mut record: Vec<String> = Vec::new();
        if crosslinks {
            record.push(first.kind.into());
        }
        record.extend([
            first.name.clone(),
            first.proteins.clone(),
            identified_in,
            q.identifications.len().to_string(),
            q.peak.score.to_string(),
            q.peak.spectral_angle.to_string(),
        ]);
        for run in &q.areas {
            record.extend(run.iter().map(|a| a.to_string()));
        }
        for run in &q.areas {
            record.extend(run[1..].iter().map(|&a| cell(silac::ratio(a, run[0]))));
        }
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    Ok(wtr.into_inner()?)
}

/// Per protein: the median ratio of its unique peptides per run and channel
fn proteins_tsv(
    quants: &[SilacQuant],
    named: &[Named],
    names: &[String],
    filenames: &[String],
    decoy_tag: &str,
) -> anyhow::Result<Vec<u8>> {
    // protein → (peptides, ratios per (run, channel))
    type Ratios = BTreeMap<(usize, usize), Vec<f64>>;
    let mut by_protein: BTreeMap<String, (usize, Ratios)> = BTreeMap::new();
    for q in quants {
        let proteins = &named[q.identifications[0]].proteins;
        if proteins.contains(';') || proteins.starts_with(decoy_tag) || proteins.is_empty() {
            continue;
        }
        let entry = by_protein.entry(proteins.clone()).or_default();
        entry.0 += 1;
        for (run, areas) in q.areas.iter().enumerate() {
            for c in 1..areas.len() {
                if let Some(r) = silac::ratio(areas[c], areas[0]) {
                    entry.1.entry((run, c)).or_default().push(r);
                }
            }
        }
    }
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
    let mut headers = vec!["protein".to_string(), "peptides".to_string()];
    for file in filenames {
        for name in &names[1..] {
            headers.push(format!("{file} {name}/light"));
            headers.push(format!("{file} {name}/light peptides"));
        }
    }
    wtr.write_record(&headers)?;
    for (protein, (n, ratios)) in &mut by_protein {
        let mut record = vec![protein.clone(), n.to_string()];
        for run in 0..filenames.len() {
            for c in 1..names.len() {
                let list = ratios.get_mut(&(run, c));
                let count = list.as_ref().map_or(0, |l| l.len());
                record.push(cell(list.and_then(|l| silac::median_ratio(l))));
                record.push(count.to_string());
            }
        }
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    Ok(wtr.into_inner()?)
}

impl Runner {
    /// `silac.sage.tsv` and `silac_proteins.sage.tsv` of a linear search
    pub(super) fn write_silac_linear(
        &self,
        features: &[Feature],
        ms1: &[ProcessedSpectrum],
        filenames: &[String],
    ) -> anyhow::Result<Vec<Url>> {
        let Some(labels) = self.database.silac.clone() else { return Ok(Vec::new()) };
        let lfq = &self.parameters.quant.lfq_settings;
        let (ids, named) = linear_ids(&self.database, &labels, features, lfq.peptide_q_value);
        if ms1.is_empty() {
            log::warn!("SILAC: no MS1 spectra, so no ratios");
        }
        let quants = silac::quantify(&ids, lfq, ms1, filenames.len(), self.parameters.precursor_charge);
        let path = self.make_path("silac.sage.tsv");
        sage_cloudpath::write_bytes_sync(&path, species_tsv(&quants, &ids, &named, &labels.names, filenames, false)?)?;
        let proteins_path = self.make_path("silac_proteins.sage.tsv");
        let proteins = proteins_tsv(&quants, &named, &labels.names, filenames, &self.database.decoy_tag)?;
        sage_cloudpath::write_bytes_sync(&proteins_path, proteins)?;
        log::info!(
            "SILAC: {} peptides quantified in {} channels ({} identifications)",
            quants.len(),
            labels.len(),
            ids.len()
        );
        Ok(vec![path, proteins_path])
    }

    /// `crosslinks.silac.tsv`
    pub(super) fn write_silac_crosslinks(
        &self,
        matches: &[CrosslinkMatch],
        settings: &CrosslinkSettings,
        ms1: &[ProcessedSpectrum],
        filenames: &[String],
    ) -> anyhow::Result<Option<Url>> {
        let Some(labels) = self.database.silac.clone() else { return Ok(None) };
        let lfq = &self.parameters.quant.lfq_settings;
        let (ids, named) = crosslink_ids(&self.database, &labels, matches, settings, lfq.peptide_q_value);
        let quants = silac::quantify(&ids, lfq, ms1, filenames.len(), self.parameters.precursor_charge);
        let path = self.make_path("crosslinks.silac.tsv");
        sage_cloudpath::write_bytes_sync(&path, species_tsv(&quants, &ids, &named, &labels.names, filenames, true)?)?;
        log::info!("SILAC: {} crosslinked species quantified ({} identifications)", quants.len(), ids.len());
        Ok(Some(path))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn tables_name_channels_and_ratios() {
        let names = vec!["light".to_string(), "heavy".to_string()];
        let ids = vec![Identified {
            key: "PEPTIDEK".into(),
            channel: 1,
            file_id: 0,
            rt: 1.0,
            ims: 0.0,
            score: 10.0,
            masses: vec![1000.0, 1008.0],
            envelope: [1.0, 0.5, 0.2],
        }];
        let named = vec![Named { kind: "linear", name: "PEPTIDEK".into(), proteins: "sp|P1|A".into() }];
        let quants = vec![SilacQuant {
            identifications: vec![0],
            peak: Default::default(),
            areas: vec![vec![100.0, 300.0], vec![0.0, 50.0]],
        }];
        let files = vec!["a.d".to_string(), "b.d".to_string()];
        let text = String::from_utf8(species_tsv(&quants, &ids, &named, &names, &files, false).unwrap()).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "peptide\tproteins\tidentified_in\tpsms\tpeak_score\tspectral_angle\ta.d light\ta.d heavy\tb.d light\tb.d heavy\ta.d heavy/light\tb.d heavy/light"
        );
        // 300 / 100 = 3 in a.d; no light in b.d, so no ratio
        assert!(lines[1].starts_with("PEPTIDEK\tsp|P1|A\theavy\t1\t"), "{}", lines[1]);
        assert!(lines[1].ends_with("\t3\t"), "{}", lines[1]);
        let proteins = String::from_utf8(proteins_tsv(&quants, &named, &names, &files, "rev_").unwrap()).unwrap();
        assert!(proteins.lines().nth(1).unwrap().starts_with("sp|P1|A\t1\t3\t1\t\t0"), "{proteins}");
    }
}
