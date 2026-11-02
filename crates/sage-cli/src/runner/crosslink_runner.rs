//! Crosslinking (XL-MS) search workflow: read spectra, search for crosslinked
//! peptide pairs, apply target/decoy FDR, and write `crosslinks.sage.tsv`.
//!
//! This is a child module of `runner` so that it can reuse `Runner`'s private
//! spectrum-reading and path helpers.

use super::Runner;
use crate::telemetry;
use log::info;
use rayon::prelude::*;
use sage_core::crosslink::{
    assign_q_values, CrosslinkMatch, CrosslinkScorer, LinkKind, PairClass,
};
use sage_core::crosslink_quant::{quantify_lfq, CrosslinkLfq};
use sage_core::mass::Tolerance;
use sage_core::scoring::Scorer;
use std::collections::HashMap;
use std::sync::atomic::Ordering;
use std::time::Instant;

/// Isobaric reporter intensities of one MS2 spectrum: (ion injection time,
/// reporter intensities), keyed by (file id, MS2 spectrum id)
/// CSM-level q-value cutoff of the filtered report
const FDR_CUTOFF: f32 = 0.01;
/// Target-target CSMs that pass [`FDR_CUTOFF`], same columns as `crosslinks.sage.tsv`
const FDR_FILTERED: &str = "crosslinks_FDRfiltered.sage.tsv";

type Reporters = HashMap<(usize, String), (f32, Vec<f32>)>;

impl Runner {
    pub(super) fn run_crosslink(
        mut self,
        parallel: usize,
        parquet: bool,
    ) -> anyhow::Result<telemetry::Telemetry> {
        let settings = self
            .parameters
            .crosslink
            .clone()
            .expect("run_crosslink requires crosslink settings");
        if parquet {
            log::warn!("parquet output is not available for crosslink searches, writing tsv");
        }

        let scorer = Scorer {
            db: &self.database,
            precursor_tol: self.parameters.precursor_tol,
            fragment_tol: self.parameters.fragment_tol,
            min_matched_peaks: self.parameters.min_matched_peaks,
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            min_precursor_charge: self.parameters.precursor_charge.0,
            max_precursor_charge: self.parameters.precursor_charge.1,
            override_precursor_charge: self.parameters.override_precursor_charge,
            max_fragment_charge: self.parameters.max_fragment_charge,
            chimera: false,
            report_psms: self.parameters.report_psms,
            wide_window: false,
            annotate_matches: false,
            score_type: self.parameters.score_type,
        };
        let kernel = sage_core::crosslink::make_kernel(&settings, &self.database)
            .map_err(anyhow::Error::msg)?;
        let xl = CrosslinkScorer::new(&scorer, &settings).with_kernel(kernel.as_deref());
        if let Some(kernel) = &kernel {
            info!("crosslink chain scoring backend: {}", kernel.name());
        }

        let lfq = self.parameters.quant.lfq;
        let tmt = self.parameters.quant.tmt.clone();
        let tmt_level = self.parameters.quant.tmt_settings.level;

        let batch_size = parallel.max(1);
        let mut matches: Vec<CrosslinkMatch> = Vec::new();
        let mut ms1_spectra = Vec::new();
        let mut reporters = Reporters::new();
        for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
            let (ms1, msn) = self.read_processed_spectra(chunk, chunk_idx, batch_size);
            if lfq {
                ms1_spectra.extend(ms1);
            }
            if let Some(isobaric) = &tmt {
                // MS3 reporters are keyed by their MS2 precursor scan, like the CSMs
                let quant = sage_core::tmt::quantify(
                    &msn,
                    isobaric,
                    Tolerance::Ppm(-20.0, 20.0),
                    tmt_level,
                );
                reporters.extend(
                    quant
                        .into_iter()
                        .map(|q| ((q.file_id, q.spec_id), (q.ion_injection_time, q.peaks))),
                );
            }

            let start = Instant::now();
            let searched = msn
                .par_iter()
                .filter(|spec| spec.level == 2 && spec.masses.len() >= self.parameters.min_peaks)
                .count();
            let found: Vec<CrosslinkMatch> = msn
                .par_iter()
                .filter(|spec| spec.level == 2 && spec.masses.len() >= self.parameters.min_peaks)
                .flat_map(|spec| xl.score(spec))
                .collect();
            info!(
                "- crosslink search: {:8} ms ({} spectra, {} candidate matches)",
                start.elapsed().as_millis(),
                searched,
                found.len()
            );
            self.spectra_searched.fetch_add(searched, Ordering::Relaxed);
            matches.extend(found);
        }

        let (csms, links) = assign_q_values(&mut matches);
        matches.sort_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));

        info!(
            "discovered {} target crosslink-spectrum matches at 1% FDR",
            csms
        );
        info!("discovered {} unique target crosslinks at 1% FDR", links);

        let filenames = self
            .parameters
            .mzml_paths
            .iter()
            .map(|url| {
                sage_cloudpath::filename(url)
                    .unwrap_or_else(|| url.as_str())
                    .to_string()
            })
            .collect::<Vec<_>>();

        let reporters_out = tmt.as_ref().map(|_| &reporters);
        let path = self.write_crosslinks("crosslinks.sage.tsv", &matches, &filenames, reporters_out)?;
        self.parameters.output_paths.push(path);
        // Only the accepted results: target-target CSMs at 1% FDR
        let accepted = matches
            .iter()
            .filter(|m| m.class == PairClass::TT && m.spectrum_q <= FDR_CUTOFF)
            .cloned()
            .collect::<Vec<_>>();
        let path = self.write_crosslinks(FDR_FILTERED, &accepted, &filenames, reporters_out)?;
        info!("wrote {} ({} target crosslink-spectrum matches at 1% FDR)", FDR_FILTERED, accepted.len());
        self.parameters.output_paths.push(path);
        self.write_viewer_exports(&matches, &settings, &filenames)?;
        if tmt.is_some() {
            let quantified = matches
                .iter()
                .filter(|m| reporters.contains_key(&(m.file_id, m.spec_id.clone())))
                .count();
            info!("isobaric reporter ions found for {} of {} crosslink-spectrum matches", quantified, matches.len());
        }

        if lfq {
            let areas = quantify_lfq(
                &matches,
                &settings,
                &self.parameters.quant.lfq_settings,
                &self.database,
                &ms1_spectra,
                self.parameters.mzml_paths.len(),
                self.parameters.precursor_charge,
            );
            let quantified = areas.iter().filter(|a| !a.decoy && a.peak.q_value <= 0.05).count();
            info!("discovered {} target crosslink MS1 peaks at 5% FDR", quantified);
            let path = self.write_crosslink_lfq(&matches, &areas, &filenames)?;
            self.parameters.output_paths.push(path);
        }

        let path = self.make_path("results.json");
        println!("{}", serde_json::to_string_pretty(&self.parameters)?);
        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);

        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);

        Ok(telemetry::Telemetry::new(
            self.parameters,
            self.database.peptides.len(),
            self.database.fragments.len(),
            self.spectra_searched.load(Ordering::Relaxed),
            false,
            run_time,
        ))
    }

    fn write_crosslinks(
        &self,
        name: &str,
        matches: &[CrosslinkMatch],
        filenames: &[String],
        reporters: Option<&Reporters>,
    ) -> anyhow::Result<sage_cloudpath::Url> {
        let path = self.make_path(name);
        let db = &self.database;

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        // Isobaric reporter columns follow the standard ones when TMT is set
        let reporter_headers = match (reporters, &self.parameters.quant.tmt) {
            (Some(_), Some(tmt)) => tmt.headers(),
            _ => Vec::new(),
        };
        let mut headers = vec![
            "psm_id",
            "filename",
            "scannr",
            "rank",
            "label",
            "class",
            "link_type",
            "peptide_a",
            "peptide_b",
            "link_a",
            "link_b",
            "proteins_a",
            "proteins_b",
            "intra_protein",
            "expmass",
            "calcmass",
            "charge",
            "isotope_error",
            "precursor_ppm",
            "hyperscore",
            "score_a",
            "score_b",
            "matched_a",
            "matched_b",
            "linked_a",
            "linked_b",
            "delta_next",
            "pairs_scored",
            "stub_peaks",
            "rt",
            "spectrum_q",
            "crosslink_q",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
        if !reporter_headers.is_empty() {
            headers.push("ion_injection_time".into());
            headers.extend(reporter_headers.iter().cloned());
        }
        wtr.write_record(&headers)?;

        for m in matches {
            let alpha = &db[m.alpha];
            let (peptide_b, proteins_b) = match m.beta {
                Some(ix) => (
                    db[ix].display(&db.modifications).to_string(),
                    db[ix].proteins(&db.decoy_tag, db.generate_decoys),
                ),
                None => (String::new(), String::new()),
            };
            // `link_a`/`link_b` are 1-based residue positions *within each peptide*.
            // For a loop-link both are positions on peptide_a.
            let mut record = vec![
                m.psm_id.to_string(),
                filenames
                    .get(m.file_id)
                    .cloned()
                    .unwrap_or_else(|| m.file_id.to_string()),
                m.spec_id.clone(),
                m.rank.to_string(),
                (if m.class == PairClass::TT { 1 } else { -1 }).to_string(),
                format!("{:?}", m.class),
                match m.kind {
                    LinkKind::Inter => "crosslink",
                    LinkKind::Loop => "loop",
                }
                .to_string(),
                alpha.display(&db.modifications).to_string(),
                peptide_b,
                (m.alpha_site + 1).to_string(),
                (m.beta_site + 1).to_string(),
                alpha.proteins(&db.decoy_tag, db.generate_decoys),
                proteins_b,
                m.intra_protein.to_string(),
                m.expmass.to_string(),
                m.calcmass.to_string(),
                m.charge.to_string(),
                m.isotope_error.to_string(),
                m.precursor_ppm.to_string(),
                m.hyperscore.to_string(),
                m.alpha_score.to_string(),
                m.beta_score.to_string(),
                m.alpha_matched.to_string(),
                m.beta_matched.to_string(),
                m.alpha_linked.to_string(),
                m.beta_linked.to_string(),
                m.delta_next.to_string(),
                m.pairs_scored.to_string(),
                m.stub_peaks.to_string(),
                m.rt.to_string(),
                m.spectrum_q.to_string(),
                m.crosslink_q.to_string(),
            ];
            if let Some(reporters) = reporters.filter(|_| !reporter_headers.is_empty()) {
                // Spectra without reporter ions (e.g. no matching MS3) get empty cells
                match reporters.get(&(m.file_id, m.spec_id.clone())) {
                    Some((injection_time, peaks)) => {
                        record.push(injection_time.to_string());
                        record.extend(peaks.iter().map(|x| x.to_string()));
                    }
                    None => record.extend(std::iter::repeat_n(String::new(), reporter_headers.len() + 1)),
                }
            }
            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// `crosslinks_xiView.csv` and `crosslinks_proXL.xml` for the xiVIEW and
    /// ProXL viewers (see `crosslink_export`)
    fn write_viewer_exports(
        &mut self,
        matches: &[CrosslinkMatch],
        settings: &sage_core::crosslink::CrosslinkSettings,
        filenames: &[String],
    ) -> anyhow::Result<()> {
        use super::crosslink_export::{proxl_xml, xiview_csv, ExportContext, SpectrumFile};
        // Protein sequences locate each peptide in its protein
        let fasta_url = sage_cloudpath::to_url(&self.parameters.database.fasta)?;
        let fasta = sage_cloudpath::util::read_fasta(
            &fasta_url,
            &self.parameters.database.decoy_tag,
            self.parameters.database.generate_decoys,
        )?;
        let wanted: std::collections::HashSet<&str> = matches
            .iter()
            .flat_map(|m| std::iter::once(m.alpha).chain(m.beta))
            .flat_map(|ix| self.database[ix].proteins.iter().map(|p| p.as_ref()))
            .collect();
        let proteins: HashMap<String, String> = fasta
            .targets
            .iter()
            .filter(|(acc, _)| wanted.contains(acc.as_ref()))
            .map(|(acc, seq)| (acc.to_string(), seq.clone()))
            .collect();
        let fasta_name = sage_cloudpath::filename(&fasta_url)
            .unwrap_or(&self.parameters.database.fasta)
            .to_string();

        // How the viewers find each match's spectrum; spectrum positions are
        // read again only where ids alone are not enough (MGF, mzML without scan=)
        let files: Vec<SpectrumFile> = filenames
            .iter()
            .enumerate()
            .map(|(file_id, name)| {
                let mut file = SpectrumFile::from_name(name);
                let ids = matches.iter().filter(|m| m.file_id == file_id).map(|m| m.spec_id.as_str());
                if file.needs_positions(ids) {
                    match self
                        .parameters
                        .mzml_paths
                        .get(file_id)
                        .ok_or_else(|| anyhow::anyhow!("no path"))
                        .and_then(|url| Ok(sage_cloudpath::util::read_spectrum_positions(url)?))
                    {
                        Ok(positions) => file.positions = Some(positions),
                        Err(e) => log::warn!(
                            "could not read spectrum positions of {name} ({e}): xiVIEW/ProXL may not find its spectra"
                        ),
                    }
                }
                file
            })
            .collect();
        let ctx = ExportContext {
            db: &self.database,
            settings,
            files: &files,
            proteins: &proteins,
            fasta_name: &fasta_name,
            fragment_tol: self.parameters.fragment_tol,
            ion_kinds: &self.database.ion_kinds,
            static_mods: &self.parameters.database.modifications.static_mods,
            variable_mods: &self.parameters.database.modifications.variable_mods,
            precursor_tol: self.parameters.precursor_tol,
            version: &self.parameters.version,
        };

        let (csv, loops) = xiview_csv(&ctx, matches)?;
        let path = self.make_path("crosslinks_xiView.csv");
        sage_cloudpath::write_bytes_sync(&path, csv)?;
        if loops > 0 {
            info!("{} loop-link matches are not in crosslinks_xiView.csv: xiVIEW's CSV format cannot describe them", loops);
        }
        let mut written = vec![path];

        match proxl_xml(&ctx, matches) {
            Some(xml) => {
                let path = self.make_path("crosslinks_proXL.xml");
                sage_cloudpath::write_bytes_sync(&path, xml.into_bytes())?;
                written.push(path);
            }
            None => log::warn!("no target crosslink matches: crosslinks_proXL.xml was not written (ProXL needs at least one)"),
        }

        // mzIdentML 1.3.0 for xiVIEW (also a local installation) and PRIDE
        match super::crosslink_mzid::mzid_xml(&ctx, matches) {
            Some(xml) => {
                let path = self.make_path("crosslinks.mzid");
                sage_cloudpath::write_bytes_sync(&path, xml.into_bytes())?;
                written.push(path);
            }
            None => log::warn!("no crosslink matches: crosslinks.mzid was not written (mzIdentML needs at least one)"),
        }
        let names = written
            .iter()
            .filter_map(|p| sage_cloudpath::filename(p).map(String::from))
            .collect::<Vec<_>>();
        info!("wrote {} for crosslink viewers", names.join(" and "));
        self.parameters.output_paths.extend(written);
        Ok(())
    }

    /// `crosslinks.lfq.tsv`: MS1 intensities of each quantified (target)
    /// crosslinked species in every run, named by its best match
    fn write_crosslink_lfq(
        &self,
        matches: &[CrosslinkMatch],
        areas: &[CrosslinkLfq],
        filenames: &[String],
    ) -> anyhow::Result<sage_cloudpath::Url> {
        let path = self.make_path("crosslinks.lfq.tsv");
        let db = &self.database;

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = [
            "peptide_a",
            "peptide_b",
            "link_type",
            "link_a",
            "link_b",
            "proteins_a",
            "proteins_b",
            "intra_protein",
            "calcmass",
            "charge",
            "crosslink_q",
            "q_value",
            "score",
            "spectral_angle",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
        headers.extend(filenames.iter().cloned());
        wtr.write_record(&headers)?;

        for quant in areas.iter().filter(|a| !a.decoy) {
            let m = &matches[quant.csm];
            let alpha = &db[m.alpha];
            let (peptide_b, proteins_b) = match m.beta {
                Some(ix) => (
                    db[ix].display(&db.modifications).to_string(),
                    db[ix].proteins(&db.decoy_tag, db.generate_decoys),
                ),
                None => (String::new(), String::new()),
            };
            let mut record = vec![
                alpha.display(&db.modifications).to_string(),
                peptide_b,
                match m.kind {
                    LinkKind::Inter => "crosslink",
                    LinkKind::Loop => "loop",
                }
                .to_string(),
                (m.alpha_site + 1).to_string(),
                (m.beta_site + 1).to_string(),
                alpha.proteins(&db.decoy_tag, db.generate_decoys),
                proteins_b,
                m.intra_protein.to_string(),
                m.calcmass.to_string(),
                quant.charge.map_or(-1, i32::from).to_string(),
                m.crosslink_q.to_string(),
                quant.peak.q_value.to_string(),
                quant.peak.score.to_string(),
                quant.peak.spectral_angle.to_string(),
            ];
            record.extend(quant.areas.iter().map(|x| x.to_string()));
            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }
}
