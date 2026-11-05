use anyhow::{ensure, Context};
use clap::ArgMatches;
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::Url;
use sage_core::scoring::ScoreType;
use sage_core::{
    crosslink::CrosslinkSettings,
    database::{Builder, Parameters},
    lfq::LfqSettings,
    mass::Tolerance,
    tmt::Isobaric,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Clone)]
/// Actual search parameters - may include overrides or default values not set by user
pub struct Search {
    pub version: String,
    pub database: Parameters,
    pub quant: QuantSettings,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: usize,
    pub max_peaks: usize,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub report_psms: usize,
    pub predict_rt: bool,
    pub mzml_paths: Vec<Url>,
    pub output_paths: Vec<Url>,
    pub bruker_config: BrukerProcessingConfig,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f32,
    /// Present only for crosslinking (XL-MS) searches
    #[serde(skip_serializing_if = "Option::is_none")]
    pub crosslink: Option<CrosslinkSettings>,

    #[serde(skip_serializing)]
    pub output_directory: Url,

    #[serde(skip_serializing)]
    pub write_pin: bool,

    #[serde(skip_serializing)]
    pub write_report: bool,

    /// Write `results.annotator.csv` / `crosslinks.annotator.csv` for Annotator
    #[serde(skip_serializing)]
    pub write_annotator: bool,

    #[serde(skip_serializing)]
    pub annotate_matches: bool,

    pub score_type: ScoreType,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    pub database: Builder,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub report_psms: Option<usize>,
    pub chimera: Option<bool>,
    pub wide_window: Option<bool>,
    pub min_peaks: Option<usize>,
    pub max_peaks: Option<usize>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: Option<u16>,
    pub precursor_charge: Option<(u8, u8)>,
    pub override_precursor_charge: Option<bool>,
    pub isotope_errors: Option<(i8, i8)>,
    pub deisotope: Option<bool>,
    pub quant: Option<QuantOptions>,
    pub predict_rt: Option<bool>,
    pub output_directory: Option<String>,
    pub mzml_paths: Option<Vec<String>>,
    pub bruker_config: Option<BrukerProcessingConfig>,
    pub protein_grouping: Option<bool>,
    pub protein_grouping_peptide_fdr: Option<f32>,
    /// Enables a chemical crosslinking search (non-cleavable linkers)
    pub crosslink: Option<CrosslinkSettings>,

    pub annotate_matches: Option<bool>,
    pub write_pin: Option<bool>,
    pub write_report: Option<bool>,
    pub write_annotator: Option<bool>,
    pub score_type: Option<ScoreType>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LfqOptions {
    pub peak_scoring: Option<sage_core::lfq::PeakScoringStrategy>,
    pub integration: Option<sage_core::lfq::IntegrationStrategy>,
    pub spectral_angle: Option<f64>,
    pub ppm_tolerance: Option<f32>,
    pub mobility_pct_tolerance: Option<f32>,
    pub combine_charge_states: Option<bool>,
    pub peptide_q_value: Option<f32>,
}

impl From<LfqOptions> for LfqSettings {
    fn from(value: LfqOptions) -> LfqSettings {
        let default = LfqSettings::default();
        let settings = LfqSettings {
            peak_scoring: value.peak_scoring.unwrap_or(default.peak_scoring),
            integration: value.integration.unwrap_or(default.integration),
            spectral_angle: value.spectral_angle.unwrap_or(default.spectral_angle).abs(),
            ppm_tolerance: value.ppm_tolerance.unwrap_or(default.ppm_tolerance).abs(),
            peptide_q_value: value.peptide_q_value.unwrap_or(default.peptide_q_value),
            mobility_pct_tolerance: value
                .mobility_pct_tolerance
                .unwrap_or(default.mobility_pct_tolerance),
            combine_charge_states: value
                .combine_charge_states
                .unwrap_or(default.combine_charge_states),
        };
        if settings.ppm_tolerance > 20.0 {
            log::warn!("lfq_settings.ppm_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance > 4.0 {
            log::warn!("lfq_settings.mobility_pct_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance < 0.05 {
            log::warn!("lfq_settings.mobility_pct_tolerance is smaller than expected");
        }
        if settings.spectral_angle < 0.50 {
            log::warn!("lfq_settings.spectral_angle is lower than expected");
        }
        if settings.peptide_q_value > 0.01 {
            log::info!("lfq_settings.peptide_q_value is higher than expected, expect increased runtime and memory usage");
        }
        if settings.peptide_q_value < 0.01 {
            log::warn!("lfq_settings.peptide_q_value is lower than expected, not all identified peptides will have MS1 intensities extracted");
        }

        settings
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TmtOptions {
    pub level: Option<u8>,
    pub sn: Option<bool>,
}

#[derive(Copy, Clone, Serialize, Debug)]
pub struct TmtSettings {
    pub level: u8,
    pub sn: bool,
}

impl From<TmtOptions> for TmtSettings {
    fn from(value: TmtOptions) -> Self {
        let default = Self::default();
        Self {
            level: value.level.unwrap_or(default.level),
            sn: value.sn.unwrap_or(default.sn),
        }
    }
}

impl Default for TmtSettings {
    fn default() -> Self {
        Self {
            level: 3,
            sn: false,
        }
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct QuantOptions {
    pub tmt: Option<Isobaric>,
    #[serde(rename = "tmt_settings")]
    pub tmt_options: Option<TmtOptions>,

    pub lfq: Option<bool>,
    #[serde(rename = "lfq_settings")]
    pub lfq_options: Option<LfqOptions>,
}

#[derive(Serialize, Default, Clone)]
pub struct QuantSettings {
    pub tmt: Option<Isobaric>,
    pub tmt_settings: TmtSettings,
    pub lfq: bool,
    pub lfq_settings: LfqSettings,
}

impl From<QuantOptions> for QuantSettings {
    fn from(value: QuantOptions) -> Self {
        Self {
            tmt: value.tmt,
            tmt_settings: value.tmt_options.map(Into::into).unwrap_or_default(),

            lfq: value.lfq.unwrap_or(false),
            lfq_settings: value.lfq_options.map(Into::into).unwrap_or_default(),
        }
    }
}

impl Input {
    pub fn from_arguments(matches: ArgMatches) -> anyhow::Result<Self> {
        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");
        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Handle JSON configuration overrides
        if let Some(output_directory) = matches.get_one::<String>("output_directory") {
            input.output_directory = Some(output_directory.into());
        }
        if let Some(fasta) = matches.get_one::<String>("fasta") {
            input.database.fasta = Some(fasta.into());
        }
        if let Some(mzml_paths) = matches.get_many::<String>("mzml_paths") {
            input.mzml_paths = Some(mzml_paths.into_iter().map(|p| p.into()).collect());
        }

        // Output switches: a flag turns an output on; without the flag the
        // configuration file decides. (`SetTrue` flags read as `false` when
        // absent, which used to override `"write_pin": true` in the file.)
        if matches.get_flag("write-pin") {
            input.write_pin = Some(true);
        }
        if matches.get_flag("write-report") {
            input.write_report = Some(true);
        }
        if matches.get_flag("write-annotator") {
            input.write_annotator = Some(true);
        }
        if matches.get_flag("annotate-matches") {
            input.annotate_matches = Some(true);
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)

        ensure!(
            input.database.fasta.is_some() || input.database.crosslink_databases.is_some(),
            "`database.fasta` (or `database.crosslink_databases`) must be set. For more information try '--help'"
        );
        if input.database.crosslink_databases.is_some() && input.database.fasta.is_some() {
            log::warn!("`database.crosslink_databases` is set: `database.fasta` is ignored");
        }
        ensure!(
            input
                .mzml_paths
                .as_ref()
                .map(|p| p.len())
                .unwrap_or_default()
                > 0,
            "`mzml_paths` must be set. For more information try '--help'"
        );

        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        sage_cloudpath::util::read_json(path).map_err(anyhow::Error::from)
    }

    fn check_mass_tolerances(tolerance: &Tolerance) {
        let (lo, hi) = match tolerance {
            Tolerance::Ppm(lo, hi) => (*lo, *hi),
            Tolerance::Pct(lo, hi) => {
                log::warn!(
                    "Pct tolerances are very rarely used for mass tolerances, did you mean ppm?"
                );
                (*lo, *hi)
            }
            Tolerance::Da(lo, hi) => (*lo, *hi),
        };
        if hi.abs() > lo.abs() {
            log::warn!(
                "Tolerances are applied to experimental masses, not theoretical: [{}, {}]",
                lo,
                hi
            );
        }
        if lo > 0.0 {
            log::warn!(
                "The `left` tolerance should probably be negative, for example: [{}, {}]",
                -lo,
                hi.abs()
            )
        }
        if hi < 0.0 {
            log::warn!(
                "The `right` tolerance should probably be positive, for example: [{}, {}]",
                -lo.abs(),
                hi
            )
        }
    }

    /// Checks and adjusts settings that only matter for crosslink searches
    fn prepare_crosslink(&mut self) -> anyhow::Result<()> {
        let Some(crosslink) = &self.crosslink else {
            return Ok(());
        };
        crosslink.validate().map_err(anyhow::Error::msg)?;

        // The prefilter keeps only peptides that identify *linearly*. A crosslink
        // partner often can't, so prefiltering would silently drop true links.
        ensure!(
            !self.database.prefilter.unwrap_or(false),
            "`database.prefilter` cannot be combined with `crosslink`: partner peptides of a \
             crosslink usually don't identify as linear peptides and would be discarded"
        );
        ensure!(
            !self.wide_window.unwrap_or(false),
            "`wide_window` is not supported for crosslink searches"
        );
        if self.chimera.unwrap_or(false) {
            log::warn!("`chimera` is ignored for crosslink searches");
        }
        // A residue carrying a static modification can't be linked: with TMT on K
        // as a static mod, no lysine could carry a crosslink. Labelling after
        // crosslinking leaves linked lysines unlabelled, which a *variable* label
        // on K models (Sage then only links the unlabelled form).
        let reactive = match crosslink.hetero_sites() {
            Some((site1, site2)) => site1 + &site2,
            None => crosslink.residues.clone(),
        };
        let blocked = self
            .database
            .static_mods
            .iter()
            .flatten()
            .filter_map(|(target, _)| match target {
                sage_core::modification::ModificationSpecificity::Residue(aa)
                    if reactive.as_bytes().contains(aa) =>
                {
                    Some(*aa as char)
                }
                _ => None,
            })
            .collect::<String>();
        if !blocked.is_empty() {
            log::warn!(
                "static modification(s) on crosslink-reactive residue(s) `{}`: these residues can't \
                 carry a crosslink. For isobaric labels (TMT) applied after crosslinking, make the \
                 label on these residues a variable modification instead (keep it static on the \
                 peptide N-terminus, `^`)",
                blocked
            );
        }
        if self.quant.as_ref().and_then(|q| q.lfq).unwrap_or(false) {
            log::info!(
                "crosslink LFQ: quantifying crosslinked peptide pairs in MS1 (crosslinks.lfq.tsv)"
            );
        }

        if crosslink.is_cleavable() {
            if crosslink.loop_links {
                log::info!("MS-cleavable crosslinker: loop-links are not searched");
            }
            // Stub-peptide (doublet) peaks are often weak; keep more peaks than the
            // linear default of 150 unless the user chose a value
            if self.max_peaks.is_none() {
                log::info!("cleavable crosslink search: setting `max_peaks` to 300");
                self.max_peaks = Some(300);
            }
            if !self.deisotope.unwrap_or(true) {
                log::warn!(
                    "cleavable crosslink search works best with `deisotope: true`, which \
                     assigns charge states to the stub-peptide ions"
                );
            }
        }

        // Crosslinked partners can be small; the linear default (500 Da) would hide them
        if self.database.peptide_min_mass.is_none() {
            log::info!("crosslink search: setting `database.peptide_min_mass` to 350 Da");
            self.database.peptide_min_mass = Some(350.0);
        }

        // A crosslinked lysine is not cleaved by trypsin, so linked peptides contain
        // an *internal* K. With 0 missed cleavages there is no such peptide.
        let missed_cleavages = match &self.database.enzyme {
            None => 0, // Sage's default enzyme allows no missed cleavages
            Some(enzyme) => enzyme.missed_cleavages.unwrap_or(1),
        };
        if missed_cleavages == 0 {
            log::warn!(
                "crosslink search with `database.enzyme.missed_cleavages` = 0: crosslinked \
                 residues block cleavage, so no linkable peptides will be generated. \
                 Set `missed_cleavages` to at least 1 (typically 2)."
            );
        }
        Ok(())
    }

    pub fn build(mut self) -> anyhow::Result<Search> {
        self.prepare_crosslink()?;

        if let Some(enzyme) = &self.database.enzyme {
            enzyme.validate().map_err(anyhow::Error::msg)?;
            if enzyme.enzymes.is_some() && enzyme.cleave_at.is_some() {
                log::warn!("`database.enzyme.enzymes` is set: `cleave_at`, `restrict` and `c_terminal` are ignored");
            }
        }

        let database = self
            .database
            .make_parameters()
            .context("invalid modification configuration")?;

        Self::check_mass_tolerances(&self.fragment_tol);
        Self::check_mass_tolerances(&self.precursor_tol);

        if let Some(isotope_errors) = self.isotope_errors {
            if isotope_errors.0 > isotope_errors.1 {
                log::error!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
                std::process::exit(1);
            }
        }
        if let Some(charges) = self.precursor_charge {
            if charges.0 > charges.1 {
                log::error!(
                    "Precursor charges should be specified [low, high], user provided: [{}, {}]",
                    charges.0,
                    charges.1
                );
                std::process::exit(1);
            }
        }

        if !self.predict_rt.unwrap_or(true)
            && self.quant.as_ref().and_then(|q| q.lfq).unwrap_or(false)
        {
            log::warn!(
                "`predict_rt: false` and `lfq: true` are incompatible. Setting `predict_rt: true`"
            );
            self.predict_rt = Some(true);
        }

        let mzml_paths = self
            .mzml_paths
            .expect("'mzml_paths' must be provided!")
            .iter()
            .map(|s| sage_cloudpath::to_url(s))
            .collect::<Result<Vec<_>, _>>()?;

        let output_directory = match self.output_directory {
            Some(path) => {
                match sage_cloudpath::try_parse_url(&path) {
                    Some(mut url) => {
                        // Valid URL, might still be a local directory that doesn't exist
                        if url.scheme() == "file" {
                            let path = url.to_file_path().expect("url scheme is file");
                            std::fs::create_dir_all(path)?;
                        }

                        if !url.path().ends_with("/") {
                            url.set_path(&format!("{}/", url.path()));
                        }
                        url
                    }
                    None => {
                        // Treat as a local path (covers Windows `C:\...` which
                        // otherwise parses as a URL with scheme `c`).
                        let path = std::path::Path::new(&path);
                        std::fs::create_dir_all(path)?;
                        Url::from_directory_path(path.canonicalize()?).expect("valid path")
                    }
                }
            }
            None => {
                let dir = std::env::current_dir()?;
                Url::from_directory_path(dir).expect("valid path")
            }
        };

        let score_type = self.score_type.unwrap_or(ScoreType::SageHyperScore);

        Ok(Search {
            version: clap::crate_version!().into(),
            database,
            quant: self.quant.map(Into::into).unwrap_or_default(),
            mzml_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            report_psms: self.report_psms.unwrap_or(1),
            max_peaks: self.max_peaks.unwrap_or(150),
            min_peaks: self.min_peaks.unwrap_or(15),
            min_matched_peaks: self.min_matched_peaks.unwrap_or(4),
            max_fragment_charge: self.max_fragment_charge,
            annotate_matches: self.annotate_matches.unwrap_or(false),
            precursor_charge: self.precursor_charge.unwrap_or((2, 4)),
            override_precursor_charge: self.override_precursor_charge.unwrap_or(false),
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            wide_window: self.wide_window.unwrap_or(false),
            predict_rt: self.predict_rt.unwrap_or(true),
            output_paths: Vec::new(),
            write_pin: self.write_pin.unwrap_or(false),
            bruker_config: self.bruker_config.unwrap_or_default(),
            write_report: self.write_report.unwrap_or(false),
            write_annotator: self.write_annotator.unwrap_or(false),
            protein_grouping: self.protein_grouping.unwrap_or(true),
            protein_grouping_peptide_fdr: self.protein_grouping_peptide_fdr.unwrap_or(0.01),
            crosslink: self.crosslink,
            score_type,
        })
    }
}

#[cfg(test)]
mod test {

    use sage_core::{database::EnzymeBuilder, enzyme::EnzymeParameters};

    #[test]
    fn deserialize_enzyme_builder() -> Result<(), serde_json::Error> {
        let a: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
        }))?;
        let b: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "P",
        }))?;
        let c: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "",
        }))?;

        let a: EnzymeParameters = a.into();
        let b: EnzymeParameters = b.into();
        let c: EnzymeParameters = c.into();

        assert_eq!(a.enzyme.map(|e| e.skip_suffix), Some([false; 26]));
        {
            let mut expected = [false; 26];
            expected[(b'P' - b'A') as usize] = true;
            assert_eq!(b.enzyme.map(|e| e.skip_suffix), Some(expected));
        }
        assert_eq!(c.enzyme.map(|e| e.skip_suffix), Some([false; 26]));

        Ok(())
    }

    #[test]
    fn deserialize_enzyme_list() -> Result<(), serde_json::Error> {
        let builder: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "missed_cleavages": 1,
            "min_len": 2,
            "enzymes": ["Trypsin", {"cleave_at": "D", "c_terminal": false}],
        }))?;
        assert!(builder.validate().is_ok());
        let params: EnzymeParameters = builder.into();
        let peptides = params
            .digest("MADEEKLPPGWEKR", std::sync::Arc::default())
            .into_iter()
            .map(|d| d.sequence)
            .collect::<Vec<_>>();
        assert!(peptides.contains(&"DEEK".to_string()));

        let unknown: EnzymeBuilder =
            serde_json::from_value(serde_json::json!({ "enzymes": ["pepsin"] }))?;
        assert!(unknown.validate().unwrap_err().contains("unknown enzyme"));
        let empty: EnzymeBuilder = serde_json::from_value(serde_json::json!({ "enzymes": [] }))?;
        assert!(empty.validate().is_err());
        Ok(())
    }
}
