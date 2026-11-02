//! Follows a running search through its log: the stage it is in, how many
//! files have been searched, the metrics Sage reports, warnings and errors,
//! and, when it fails, what probably went wrong and how to fix it.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// Phases of a search, in the order they happen
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Stage {
    Convert,
    Start,
    Database,
    Linear,
    Search,
    Rescore,
    Proteins,
    Quant,
    Results,
}

impl Stage {
    pub fn label(self) -> &'static str {
        match self {
            Stage::Convert => "Convert Thermo .raw files to mzML",
            Stage::Start => "Start sage.exe and check settings",
            Stage::Database => "Digest proteins and build the fragment index",
            Stage::Linear => "Linear peptide search (results.sage.tsv)",
            Stage::Search => "Search spectrum files",
            Stage::Rescore => "Align retention times and rescore matches",
            Stage::Proteins => "Infer proteins",
            Stage::Quant => "Quantify",
            Stage::Results => "Estimate FDR and write results",
        }
    }

    /// Share of the progress bar (relative)
    fn weight(self) -> f32 {
        match self {
            Stage::Convert => 3.0,
            Stage::Start => 0.2,
            Stage::Database => 0.8,
            Stage::Linear => 4.0,
            Stage::Search => 6.0,
            Stage::Rescore => 0.8,
            Stage::Proteins => 0.3,
            Stage::Quant => 1.0,
            Stage::Results => 0.5,
        }
    }
}

/// What the search is set up to do, to know its stages in advance
#[derive(Clone, Debug)]
pub struct Plan {
    pub files: usize,
    pub raw_files: usize,
    pub crosslink: bool,
    /// Crosslink searches: a linear peptide search runs first
    pub linear: bool,
    pub quant: bool,
    pub gpu_backend: bool,
    pub output_dir: PathBuf,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum State {
    Running,
    Finished,
    Failed,
    Cancelled,
}

/// Explanation of a failed search
#[derive(Clone, Debug, PartialEq)]
pub struct Diagnosis {
    pub title: String,
    pub cause: String,
    pub fixes: Vec<String>,
    /// The log lines that show the error
    pub evidence: Vec<String>,
}

pub struct Tracker {
    pub plan: Plan,
    pub stages: Vec<Stage>,
    /// Index into `stages` of the stage in progress (or last reached)
    pub current: usize,
    /// Start of each stage, and how long finished ones took
    started: Vec<Option<Instant>>,
    pub took: Vec<Option<Duration>>,
    /// Log lines (by `Job::total` count) already read
    seen: usize,
    pub converted: usize,
    /// A conversion is in progress (its header was seen)
    converting: bool,
    pub files_done: usize,
    pub files_started: usize,
    /// In the linear search that precedes a crosslink search
    in_linear: bool,
    /// The linear search reported its results; the crosslink search is next
    linear_done: bool,
    /// Metrics in the order Sage reports them: (label, value)
    pub metrics: Vec<(String, String)>,
    search_ms: u64,
    io_ms: u64,
    spectra: u64,
    candidates: u64,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub start: Instant,
    pub end: Option<Duration>,
    pub state: State,
    pub diagnosis: Option<Diagnosis>,
    /// Output files listed in results.json
    pub outputs: Vec<(PathBuf, u64)>,
}

/// The message of a Sage log line `[time LEVEL module] message`, and its level
fn parse_log(line: &str) -> Option<(&str, &str)> {
    let rest = line.strip_prefix('[')?;
    let (head, msg) = rest.split_once("] ")?;
    let level = head.split_whitespace().nth(1)?;
    Some((level, msg))
}

fn number(text: &str) -> Option<u64> {
    text.trim().split_whitespace().next()?.parse().ok()
}

/// Friendly names for Sage's "discovered N <what> at X% FDR" lines
fn result_label(what: &str) -> String {
    let what = what.trim();
    let name = match what {
        "target peptide-spectrum matches" => "Peptide-spectrum matches (PSMs)",
        "target peptides" => "Peptides",
        "target proteins (supported by proteotypic peptides only)" => "Proteins (proteotypic peptides)",
        "target protein groups (supported by proteotypic peptides only)" => "Protein groups (proteotypic peptides)",
        "target MS1 peaks" => "Quantified precursors (LFQ)",
        "target crosslink-spectrum matches" => "Crosslink-spectrum matches (CSMs)",
        "unique target crosslinks" => "Unique crosslinks",
        "target crosslink MS1 peaks" => "Quantified crosslinks (LFQ)",
        other => other,
    };
    name.to_string()
}

fn format_ms(ms: u64) -> String {
    if ms >= 60_000 {
        format!("{}m {:02}s", ms / 60_000, (ms % 60_000) / 1000)
    } else if ms >= 1000 {
        format!("{:.1} s", ms as f64 / 1000.0)
    } else {
        format!("{ms} ms")
    }
}

/// Rust's `{:?}` of a Duration ("5.5147ms", "1.234567s") rounded for people
fn tidy_duration(text: &str) -> String {
    let split = text.find(|c: char| !(c.is_ascii_digit() || c == '.')).unwrap_or(text.len());
    let (value, unit) = text.split_at(split);
    let Ok(value) = value.parse::<f64>() else { return text.to_string() };
    let ms = match unit {
        "s" => value * 1000.0,
        "ms" => value,
        "µs" | "us" => value / 1000.0,
        "ns" => value / 1e6,
        _ => return text.to_string(),
    };
    if ms < 1.0 {
        "< 1 ms".into()
    } else {
        format_ms(ms.round() as u64)
    }
}

pub fn format_duration(d: Duration) -> String {
    format_ms(d.as_millis() as u64)
}

impl Tracker {
    pub fn new(plan: Plan) -> Self {
        let mut stages = Vec::new();
        if plan.raw_files > 0 {
            stages.push(Stage::Convert);
        }
        stages.extend([Stage::Start, Stage::Database]);
        if plan.crosslink && plan.linear {
            stages.push(Stage::Linear);
        }
        stages.push(Stage::Search);
        if plan.crosslink {
            // crosslinks: FDR right after the search, then quantification
            stages.push(Stage::Results);
            if plan.quant {
                stages.push(Stage::Quant);
            }
        } else {
            stages.extend([Stage::Rescore, Stage::Proteins]);
            if plan.quant {
                stages.push(Stage::Quant);
            }
            stages.push(Stage::Results);
        }
        let n = stages.len();
        let mut t = Tracker {
            plan,
            stages,
            current: 0,
            started: vec![None; n],
            took: vec![None; n],
            seen: 0,
            converted: 0,
            converting: false,
            files_done: 0,
            files_started: 0,
            in_linear: false,
            linear_done: false,
            metrics: Vec::new(),
            search_ms: 0,
            io_ms: 0,
            spectra: 0,
            candidates: 0,
            warnings: Vec::new(),
            errors: Vec::new(),
            start: Instant::now(),
            end: None,
            state: State::Running,
            diagnosis: None,
            outputs: Vec::new(),
        };
        t.started[0] = Some(t.start);
        t
    }

    pub fn stage(&self) -> Stage {
        self.stages[self.current]
    }

    pub fn elapsed(&self) -> Duration {
        self.end.unwrap_or_else(|| self.start.elapsed())
    }

    /// Move forward to `stage` (never back); stages in between count as done
    fn reach(&mut self, stage: Stage) {
        let Some(i) = self.stages.iter().position(|&s| s == stage) else { return };
        if i <= self.current {
            return;
        }
        let now = Instant::now();
        for k in self.current..i {
            if self.took[k].is_none() {
                self.took[k] = Some(self.started[k].map_or(Duration::ZERO, |s| now - s));
            }
        }
        self.current = i;
        self.started[i] = Some(now);
    }

    fn set_metric(&mut self, label: &str, value: String) {
        match self.metrics.iter_mut().find(|(l, _)| l == label) {
            Some(slot) => slot.1 = value,
            None => self.metrics.push((label.to_string(), value)),
        }
    }

    /// Read the lines added to the job's log since the last call.
    /// `total` is the number of lines ever written (`Job::total`).
    pub fn update(&mut self, log: &[String], total: usize) {
        let new = total.saturating_sub(self.seen).min(log.len());
        for line in &log[log.len() - new..] {
            self.line(line);
        }
        self.seen = total;
    }

    fn line(&mut self, line: &str) {
        // Lines from the GUI's runner
        if let Some(rest) = line.strip_prefix("--- Convert ") {
            self.reach(Stage::Convert);
            // a new header means the previous conversion finished
            if self.converting {
                self.converted += 1;
            }
            self.converting = !rest.contains("already exists, skipped");
            if !self.converting {
                self.converted += 1;
            }
            return;
        }
        if line.starts_with("--- Sage search ---") {
            if self.stage() == Stage::Convert {
                self.converted = self.plan.raw_files;
                self.converting = false;
            }
            self.reach(Stage::Start);
            return;
        }
        if line.starts_with("RAM limit: ") {
            self.set_metric("RAM limit", line["RAM limit: ".len()..].to_string());
            return;
        }
        if line.starts_with("=== ") {
            return;
        }

        let Some((level, msg)) = parse_log(line) else {
            // Not a log line: errors sage prints on exit ("Error: ...",
            // "Caused by:", panics, allocation failures)
            let t = line.trim();
            if t.starts_with("Error:") || t.contains("panicked at") || t.starts_with("memory allocation of") {
                self.errors.push(t.to_string());
            } else if !self.errors.is_empty() && !t.is_empty() && !t.starts_with(['{', '}', '[', ']', '"']) {
                // continuation of an error message ("Caused by:", its lines)
                if self.errors.len() < 30 {
                    self.errors.push(t.to_string());
                }
            }
            return;
        };
        match level {
            "WARN" => self.warnings.push(msg.to_string()),
            "ERROR" => self.errors.push(msg.to_string()),
            _ => {}
        }

        if msg.starts_with("linear peptide search (results.sage.tsv)") {
            self.reach(Stage::Linear);
            self.in_linear = true;
            return;
        }
        if self.in_linear {
            // the linear search's own stages (rescoring, proteins, quantification)
            // all belong to this one step; its first new batch of files after
            // the results starts the crosslink search
            if let Some(rest) = msg.strip_prefix("processing files ") {
                if self.linear_done {
                    self.in_linear = false;
                    self.files_done = 0;
                    self.files_started = 0;
                } else {
                    if let Some(n) = rest.split_once("..").and_then(|(_, end)| number(end)) {
                        self.files_started = n as usize;
                    }
                    return;
                }
            } else {
                if let Some(rest) = msg.strip_prefix("- file IO:") {
                    self.io_ms += number(rest).unwrap_or(0);
                } else if let Some(rest) = msg.strip_prefix("- search:") {
                    self.search_ms += number(rest).unwrap_or(0);
                    self.files_done = self.files_started;
                } else if let Some(rest) = msg.strip_prefix("discovered ") {
                    self.linear_done = true;
                    if let Some((count, what)) = rest.split_once(' ') {
                        if let Some((what, fdr)) = what.rsplit_once(" at ") {
                            self.set_metric(&format!("{} at {}", result_label(what), fdr), count.to_string());
                        }
                    }
                }
                return;
            }
        }

        if let Some(rest) = msg.strip_prefix("using ").filter(|r| r.contains("CPU thread")) {
            self.reach(Stage::Database);
            if let Some(n) = number(rest) {
                self.set_metric("CPU threads", n.to_string());
            }
        } else if let Some(rest) = msg.strip_prefix("generated ") {
            // "generated 5500 fragments, 168 peptides in 2.4083ms"
            self.reach(Stage::Database);
            let mut parts = rest.split(", ");
            if let (Some(frag), Some(pep)) = (parts.next(), parts.next()) {
                if let Some(n) = number(pep) {
                    self.set_metric("Peptides in database", n.to_string());
                }
                if let Some(n) = number(frag) {
                    self.set_metric("Fragments indexed", n.to_string());
                }
                if let Some((_, t)) = pep.split_once(" in ") {
                    self.set_metric("Database build time", tidy_duration(t.trim()));
                }
            }
        } else if let Some(rest) = msg.strip_prefix("processing files ") {
            // "processing files 0 .. 4 " (a batch searched in parallel)
            self.reach(Stage::Search);
            if let Some((_, end)) = rest.split_once("..") {
                if let Some(n) = number(end) {
                    self.files_started = n as usize;
                }
            }
        } else if let Some(rest) = msg.strip_prefix("- file IO:") {
            self.io_ms += number(rest).unwrap_or(0);
        } else if let Some(rest) = msg.strip_prefix("- search:") {
            // "- search:  1234 ms (1000 spectra/s)"
            self.search_ms += number(rest).unwrap_or(0);
            self.files_done = self.files_started;
        } else if let Some(rest) = msg.strip_prefix("- crosslink search:") {
            // "- crosslink search:  12 ms (1 spectra, 0 candidate matches)"
            self.search_ms += number(rest).unwrap_or(0);
            if let Some((_, inner)) = rest.split_once('(') {
                let mut parts = inner.split(", ");
                self.spectra += parts.next().and_then(number).unwrap_or(0);
                self.candidates += parts.next().and_then(number).unwrap_or(0);
            }
            self.files_done = self.files_started;
        } else if let Some(rest) = msg.strip_prefix("crosslink chain scoring backend: ") {
            self.set_metric("Crosslink backend", rest.to_string());
        } else if msg.starts_with("aligning file") || msg.starts_with("aligned retention times")
            || msg.contains("model") || msg.contains("discriminant")
        {
            self.reach(Stage::Rescore);
        } else if msg.starts_with("Protein grouping") || msg.starts_with("Grouped and inferred") {
            self.reach(Stage::Proteins);
        } else if msg.starts_with("tracing MS1") || msg.starts_with("integrating MS1") || msg.starts_with("isobaric reporter") {
            self.reach(Stage::Quant);
            if let Some(rest) = msg.strip_prefix("isobaric reporter ions found for ") {
                self.set_metric("TMT-quantified CSMs", rest.replace(" crosslink-spectrum matches", ""));
            }
        } else if let Some(rest) = msg.strip_prefix("discovered ") {
            // "discovered 1234 target peptides at 1% FDR"
            let quant = rest.contains("MS1 peaks");
            self.reach(if quant { Stage::Quant } else { Stage::Results });
            if let Some((count, what)) = rest.split_once(' ') {
                if let Some((what, fdr)) = what.rsplit_once(" at ") {
                    self.set_metric(&format!("{} at {}", result_label(what), fdr), count.to_string());
                }
            }
        } else if let Some(rest) = msg.strip_prefix("finished in ") {
            self.set_metric("Sage run time", rest.to_string());
        }
    }

    /// 0..1, from the stage and, while searching, the files done
    pub fn fraction(&self) -> f32 {
        if self.state == State::Finished {
            return 1.0;
        }
        let total: f32 = self.stages.iter().map(|s| s.weight()).sum();
        let done: f32 = self.stages[..self.current].iter().map(|s| s.weight()).sum();
        let within = match self.stage() {
            Stage::Search | Stage::Linear if self.plan.files > 0 => {
                // half credit for the batch in progress
                (self.files_done as f32 + 0.5 * (self.files_started - self.files_done) as f32) / self.plan.files as f32
            }
            Stage::Convert if self.plan.raw_files > 0 => self.converted as f32 / self.plan.raw_files as f32,
            _ => 0.3,
        };
        ((done + self.stage().weight() * within.clamp(0.0, 1.0)) / total).clamp(0.0, 0.99)
    }

    /// One line on what is happening now
    pub fn activity(&self) -> String {
        match self.stage() {
            Stage::Convert => format!(
                "Converting Thermo .raw file {} of {}",
                (self.converted + 1).min(self.plan.raw_files),
                self.plan.raw_files
            ),
            Stage::Search if self.files_started == self.files_done + 1 => {
                format!("Searching file {} of {}", self.files_started, self.plan.files)
            }
            Stage::Search if self.files_started > self.files_done => format!(
                "Searching files {}–{} of {}",
                self.files_done + 1,
                self.files_started,
                self.plan.files
            ),
            Stage::Search => format!("Searched {} of {} files", self.files_done, self.plan.files),
            Stage::Linear if self.linear_done => "Linear peptide search: FDR and results".into(),
            Stage::Linear => format!(
                "Linear peptide search: file {} of {}",
                self.files_started.max(1).min(self.plan.files),
                self.plan.files
            ),
            s => s.label().to_string(),
        }
    }

    /// The job ended: record how, finish the metrics, explain a failure
    pub fn finish(&mut self, state: State, failure: Option<&str>, log: &[String]) {
        if self.state != State::Running {
            return;
        }
        self.state = state;
        let end = self.start.elapsed();
        self.end = Some(end);
        let now = Instant::now();
        if state == State::Finished {
            let last = self.stages.len() - 1;
            self.reach(self.stages[last]);
            self.took[last] = Some(self.started[last].map_or(Duration::ZERO, |s| now - s));
        }
        self.set_metric("Spectrum files", self.plan.files.to_string());
        if self.spectra > 0 {
            self.set_metric("Spectra searched", self.spectra.to_string());
        }
        if self.candidates > 0 {
            self.set_metric("Candidate crosslink matches", self.candidates.to_string());
        }
        if self.search_ms > 0 {
            self.set_metric("Search time", format_ms(self.search_ms));
        }
        if self.io_ms > 0 {
            self.set_metric("File reading time", format_ms(self.io_ms));
        }
        self.set_metric("Total time", format_duration(end));
        if state == State::Finished {
            self.outputs = read_outputs(&self.plan.output_dir);
        }
        if state == State::Failed {
            self.diagnosis = Some(diagnose(failure.unwrap_or(""), &self.errors, log, &self.plan, self.stage()));
        }
    }

    /// Result files a finished crosslink search was expected to write but did
    /// not, with the reason (e.g. no ProXL file when nothing was found)
    pub fn missing_exports(&self) -> Vec<String> {
        if !self.plan.crosslink || self.state != State::Finished {
            return Vec::new();
        }
        let written = |name: &str| self.outputs.iter().any(|(p, _)| p.file_name().is_some_and(|n| n == name));
        let mut notes = Vec::new();
        if !written("crosslinks_proXL.xml") {
            let why = if self.warnings.iter().any(|w| w.contains("crosslinks_proXL.xml was not written")) {
                "no target crosslinks were found, and ProXL rejects a file without any"
            } else {
                "this sage.exe may be older than the ProXL export; rebuild it"
            };
            notes.push(format!("crosslinks_proXL.xml was not written: {why}."));
        }
        if !written("crosslinks_xiView.csv") {
            notes.push("crosslinks_xiView.csv was not written: this sage.exe may be older than the xiVIEW export; rebuild it.".into());
        }
        if !written("crosslinks.mzid") {
            let why = if self.warnings.iter().any(|w| w.contains("crosslinks.mzid was not written")) {
                "no crosslinks were found, and mzIdentML requires at least one identification"
            } else {
                "this sage.exe may be older than the mzIdentML export; rebuild it"
            };
            notes.push(format!("crosslinks.mzid was not written: {why}."));
        }
        notes
    }

    /// Main identification count, if Sage reported it
    pub fn identifications(&self) -> Option<u64> {
        let key = if self.plan.crosslink { "Crosslink-spectrum matches (CSMs)" } else { "Peptide-spectrum matches (PSMs)" };
        self.metrics.iter().find(|(l, _)| l.starts_with(key)).and_then(|(_, v)| v.parse().ok())
    }
}

/// Output files listed in `results.json`, with their sizes
fn read_outputs(dir: &Path) -> Vec<(PathBuf, u64)> {
    let Ok(text) = std::fs::read_to_string(dir.join("results.json")) else {
        return Vec::new();
    };
    let Ok(json) = serde_json::from_str::<serde_json::Value>(&text) else {
        return Vec::new();
    };
    let mut out: Vec<(PathBuf, u64)> = json
        .get("output_paths")
        .and_then(|v| v.as_array())
        .into_iter()
        .flatten()
        .filter_map(|v| v.as_str())
        .map(|p| {
            // paths may be file:// URLs
            let p = p.strip_prefix("file:///").unwrap_or(p);
            let path = PathBuf::from(p);
            let size = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
            (path, size)
        })
        .collect();
    out.dedup_by(|a, b| a.0 == b.0);
    out
}

/// Explain a failure from the runner's message, Sage's error lines and the log
fn diagnose(failure: &str, errors: &[String], log: &[String], plan: &Plan, stage: Stage) -> Diagnosis {
    let tail: Vec<&String> = log.iter().rev().take(40).collect();
    let text = std::iter::once(failure.to_string())
        .chain(errors.iter().cloned())
        .chain(tail.iter().rev().map(|s| s.to_string()))
        .collect::<Vec<_>>()
        .join("\n");
    let lower = text.to_ascii_lowercase();
    let has = |needles: &[&str]| needles.iter().any(|n| lower.contains(&n.to_ascii_lowercase()));
    let mut evidence: Vec<String> = errors.iter().take(12).cloned().collect();
    if evidence.is_empty() && !failure.is_empty() {
        evidence.push(failure.to_string());
    }
    let d = |title: &str, cause: &str, fixes: &[&str]| Diagnosis {
        title: title.into(),
        cause: cause.into(),
        fixes: fixes.iter().map(|f| f.to_string()).collect(),
        evidence: Vec::new(),
    };

    let mut found = if has(&["ran out of memory: it needed more than"]) {
        d("Out of memory (RAM limit reached)",
          "sage.exe needed more memory than the Max RAM limit set under Files > Resources, so Windows stopped it.",
          &["Raise Max RAM under Files > Resources, or set it to 0 to allow all available memory.",
            "Search fewer files at once: fewer CPU threads also means fewer files are loaded together.",
            "Reduce the search space: fewer variable modifications, fewer missed cleavages, or a smaller FASTA."])
    } else if has(&["memory allocation of", "out of memory"]) {
        d("Out of memory",
          "The computer ran out of memory during the search.",
          &["Close other programs and run again.",
            "Use fewer CPU threads under Files > Resources: fewer files are then searched at the same time.",
            "Reduce the search space: fewer variable modifications, fewer missed cleavages, a narrower precursor tolerance or a smaller FASTA.",
            "For crosslink searches, an MS-cleavable reagent (DSSO, DSBU) needs far less memory than a non-cleavable one."])
    } else if has(&["could not start"]) {
        d("A program could not be started",
          "The GUI could not run one of the programs it needs.",
          &["Check the paths under Files > Programs (sage.exe and the RAW converter).",
            "If antivirus software quarantined sage.exe, restore it or allow it."])
    } else if stage == Stage::Convert || failure.starts_with("Convert ") {
        d("Converting a Thermo .raw file failed",
          "The RAW converter (msconvert or ThermoRawFileParser) stopped with an error, so the search did not start.",
          &["Check that the .raw file opens in Thermo FreeStyle or Xcalibur, and is not still being copied or acquired.",
            "Use a current ProteoWizard (msconvert) installation; the Thermo vendor reader needs the 64-bit version.",
            "Try the other converter: ThermoRawFileParser instead of msconvert, or the reverse (Files > Programs).",
            "Delete any half-written .mzML in <output>\\converted_mzML and run again."])
    } else if has(&["unexpected argument", "unrecognized"]) {
        d("sage.exe does not understand the GUI's options",
          "This sage.exe is older than the GUI and does not know an option it passes (such as --threads).",
          &["Rebuild sage.exe from this repository (cargo build --release -p sage-cli) or select the matching sage.exe under Files > Programs."])
    } else if has(&["crosslink backend `gpu`", "requires building sage with", "no gpu adapter", "could not open gpu device", "only a software adapter"]) {
        d("The GPU crosslink backend is not available",
          "The crosslink search was set to the gpu backend, but sage.exe could not use a GPU.",
          &["Choose cpu or auto as the backend on the Crosslinking tab (auto falls back to the CPU).",
            "Check Files > Resources > GPU support: it needs a sage.exe built with --features gpu and a DirectX 12 capable graphics driver.",
            "Update the graphics driver."])
    } else if has(&["0xc0000005", "access_violation", "3221225477"]) {
        if plan.gpu_backend {
            d("sage.exe crashed (graphics driver)",
              "sage.exe stopped abruptly while the GPU crosslink backend was selected; this usually comes from the graphics driver.",
              &["Choose the cpu backend on the Crosslinking tab and run again.",
                "Update the graphics driver.",
                "If it happens with the cpu backend too, report it with the log (Copy log)."])
        } else {
            d("sage.exe crashed",
              "sage.exe stopped abruptly (an access violation).",
              &["Run again: if it happens once, it may have been a transient problem.",
                "Check that the spectrum files are complete and not damaged.",
                "If it keeps happening, report it to the Sage developers with the log (Copy log)."])
        }
    } else if has(&["failed to read parameters", "missing field", "unknown field", "invalid type", "invalid value"]) {
        d("Sage rejected the search settings",
          "sage_search.json contains a setting this sage.exe does not accept.",
          &["Read the error below: it names the setting.",
            "Use 'Show JSON' to look at the configuration, and 'Reset to defaults' if a loaded file contains old settings.",
            "Make sure sage.exe was built from the same version as the GUI (Files > Programs)."])
    } else if has(&["os error 112", "not enough space", "no space left", "disk full"]) {
        d("The disk is full",
          "There was not enough free disk space to write the results.",
          &["Free up space on the output drive, or choose an output folder on another drive (Files > Results)."])
    } else if has(&["os error 32", "being used by another process"]) {
        d("A file is in use by another program",
          "Sage could not write a result file because another program has it open.",
          &["Close result files that are open in Excel or another program, then run again.",
            "Or choose a new output folder (Files > Results)."])
    } else if has(&["os error 5", "access is denied", "permission denied"]) {
        d("Access denied",
          "Sage was not allowed to read an input file or write to the output folder.",
          &["Choose an output folder you can write to, such as one in Documents (Files > Results).",
            "Close result files that are open in Excel.",
            "Check that the FASTA and spectrum files are readable (not on a disconnected network drive)."])
    } else if has(&["os error 2", "os error 3", "no such file", "cannot find the file", "cannot find the path", "system cannot find"]) {
        d("A file or folder was not found",
          "Sage could not find a file it needs.",
          &["Check that the FASTA file and every spectrum file still exist at the paths on the Files tab.",
            "Network drives must be connected; paths with special characters can fail on some systems.",
            "Bruker data must be added as the .d folder, not a file inside it."])
    } else if has(&["fasta"]) && has(&["error", "failed", "invalid"]) {
        d("The FASTA database could not be read",
          "Sage could not read the protein database.",
          &["Check that the FASTA file is plain text with '>' header lines (not compressed or a Word document).",
            "Download the database again if it may be incomplete."])
    } else if has(&["mzml", "xml", "mgf", "tdf"]) && has(&["error", "failed", "invalid", "unexpected"]) {
        d("A spectrum file could not be read",
          "Sage could not read one of the spectrum files.",
          &["Check that the file is complete (not still being copied or converted).",
            "Convert it again with peak picking (centroiding) turned on.",
            "Remove the file from the list to see whether the others search fine."])
    } else if has(&["panicked at"]) {
        d("Sage hit an internal error",
          "sage.exe stopped because of an unexpected condition inside Sage.",
          &["Run again with the same settings to see whether it is reproducible.",
            "Try a different combination of settings near the one that failed.",
            "Report it to the Sage developers with the log (Copy log) and the settings (Show JSON)."])
    } else {
        d("The search stopped with an error",
          "sage.exe exited with an error that the GUI does not recognise.",
          &["Read the error lines below and the end of the log.",
            "Check the settings for unusual values (Show JSON), and that all input files are readable.",
            "If the problem persists, report it with the log (Copy log)."])
    };
    if evidence.is_empty() {
        evidence = tail.iter().rev().filter(|l| !l.trim().is_empty()).rev().take(8).map(|s| s.to_string()).collect();
    }
    found.evidence = evidence;
    found
}

#[cfg(test)]
mod test {
    use super::*;

    fn plan(crosslink: bool, quant: bool) -> Plan {
        Plan { files: 4, raw_files: 0, crosslink, linear: false, quant, gpu_backend: false, output_dir: PathBuf::new() }
    }

    fn feed(t: &mut Tracker, lines: &[&str]) -> Vec<String> {
        let log: Vec<String> = lines.iter().map(|s| s.to_string()).collect();
        t.update(&log, log.len());
        log
    }

    /// Lines as sage.exe prints them (from a real run)
    const LINEAR: &[&str] = &[
        "--- Sage search ---",
        "[2026-11-02T17:38:49Z INFO  sage] using 16 CPU thread(s)",
        "[2026-11-02T17:38:49Z INFO  sage_cli::runner] generated 5500 fragments, 168 peptides in 2.4083ms",
        "[2026-11-02T17:38:49Z INFO  sage_cli::runner] processing files 0 .. 2 ",
        "[2026-11-02T17:38:49Z INFO  sage_cli::runner] - file IO:        3 ms",
        "[2026-11-02T17:38:49Z INFO  sage_cli::runner] - search:        10 ms (1000 spectra/s)",
        "[2026-11-02T17:38:49Z INFO  sage_cli::runner] processing files 2 .. 4 ",
    ];

    #[test]
    fn follows_a_linear_search() {
        let mut t = Tracker::new(plan(false, true));
        feed(&mut t, LINEAR);
        assert_eq!(t.stage(), Stage::Search);
        assert_eq!((t.files_done, t.files_started), (2, 4));
        assert_eq!(t.activity(), "Searching files 3–4 of 4");
        let mid = t.fraction();
        assert!(mid > 0.1 && mid < 0.9, "{mid}");

        let mut all: Vec<&str> = LINEAR.to_vec();
        all.extend([
            "[2026-11-02T17:38:49Z INFO  sage_cli::runner] - search:        12 ms (1000 spectra/s)",
            "[2026-11-02T17:38:49Z INFO  sage_core::ml::retention_alignment] aligning file #0: y = 1.0000x + 0.0000",
            "[2026-11-02T17:38:49Z WARN  sage_cli::runner] linear model fitting failed, falling back to heuristic discriminant score",
            "[2026-11-02T17:38:49Z INFO  sage_core::protein_grouping] Protein grouping: 0 unique peptides (threshold=0.01) in 0ms",
            "[2026-11-02T17:38:49Z INFO  sage_core::lfq] tracing MS1 features",
            "[2026-11-02T17:38:49Z INFO  sage_cli::runner] discovered 7 target MS1 peaks at 5% FDR",
            "[2026-11-02T17:38:49Z INFO  sage_cli::runner] discovered 1234 target peptide-spectrum matches at 1% FDR",
            "[2026-11-02T17:38:49Z INFO  sage_cli::runner] discovered 456 target peptides at 1% FDR",
            "[2026-11-02T17:38:49Z INFO  sage_cli::runner] finished in 3s",
            "=== Search finished ===",
        ]);
        let mut t = Tracker::new(plan(false, true));
        let log = feed(&mut t, &all);
        assert_eq!(t.stage(), Stage::Results);
        t.finish(State::Finished, None, &log);
        assert_eq!(t.fraction(), 1.0);
        assert!(t.took.iter().all(Option::is_some), "{:?}", t.took);
        let m = |k: &str| t.metrics.iter().find(|(l, _)| l == k).map(|(_, v)| v.clone());
        assert_eq!(m("Peptides in database").as_deref(), Some("168"));
        assert_eq!(m("Peptide-spectrum matches (PSMs) at 1% FDR").as_deref(), Some("1234"));
        assert_eq!(m("Quantified precursors (LFQ) at 5% FDR").as_deref(), Some("7"));
        assert_eq!(m("Search time").as_deref(), Some("22 ms"));
        assert_eq!(t.identifications(), Some(1234));
        assert_eq!(t.warnings.len(), 1);
    }

    #[test]
    fn follows_a_linear_search_before_the_crosslinks() {
        let mut t = Tracker::new(Plan { linear: true, ..plan(true, false) });
        assert_eq!(t.stages, [Stage::Start, Stage::Database, Stage::Linear, Stage::Search, Stage::Results]);
        let mut log: Vec<String> = Vec::new();
        let mut step = |t: &mut Tracker, msgs: &[&str]| {
            log.extend(msgs.iter().map(|m| format!("[2026-11-02T22:00:07Z INFO  sage_cli::runner] {m}")));
            t.update(&log, log.len());
        };
        step(&mut t, &["using 16 CPU thread(s)",
            "linear peptide search (results.sage.tsv), then the crosslink search",
            "processing files 0 .. 2 "]);
        assert_eq!(t.stage(), Stage::Linear);
        assert_eq!(t.activity(), "Linear peptide search: file 2 of 4");
        step(&mut t, &["- file IO:        4 ms", "- search:         1 ms (500 spectra/s)",
            "aligning file #0: y = 1.0000x + 0.0000", "Protein grouping: 0 unique peptides (threshold=0.01) in 0ms",
            "discovered 12 target peptide-spectrum matches at 1% FDR",
            "discovered 10 target peptides at 1% FDR"]);
        // rescoring and protein inference stay inside the linear step
        assert_eq!(t.stage(), Stage::Linear);
        step(&mut t, &["processing files 0 .. 4 "]);
        assert_eq!(t.stage(), Stage::Search);
        assert_eq!((t.files_started, t.files_done), (4, 0));
        step(&mut t, &["- crosslink search:        1 ms (1 spectra, 0 candidate matches)",
            "discovered 3 target crosslink-spectrum matches at 1% FDR"]);
        assert_eq!(t.stage(), Stage::Results);
        let m = |k: &str| t.metrics.iter().find(|(l, _)| l == k).map(|(_, v)| v.clone());
        assert_eq!(m("Peptide-spectrum matches (PSMs) at 1% FDR").as_deref(), Some("12"));
        assert_eq!(t.identifications(), Some(3));
    }

    #[test]
    fn follows_a_crosslink_search_and_reads_new_lines_only() {
        let mut t = Tracker::new(plan(true, false));
        assert_eq!(t.stages, [Stage::Start, Stage::Database, Stage::Search, Stage::Results]);
        let mut log: Vec<String> = Vec::new();
        for (i, line) in [
            "--- Sage search ---",
            "[2026-11-02T17:38:50Z INFO  sage_cli::runner] processing files 0 .. 4 ",
            "[2026-11-02T17:38:50Z INFO  sage_cli::runner::crosslink_runner] - crosslink search:        5 ms (120 spectra, 9 candidate matches)",
            "[2026-11-02T17:38:50Z INFO  sage_cli::runner::crosslink_runner] discovered 3 target crosslink-spectrum matches at 1% FDR",
        ]
        .iter()
        .enumerate()
        {
            log.push(line.to_string());
            // reading after every line must not count any twice
            t.update(&log, i + 1);
            t.update(&log, i + 1);
        }
        t.finish(State::Finished, None, &log);
        let m = |k: &str| t.metrics.iter().find(|(l, _)| l == k).map(|(_, v)| v.clone());
        assert_eq!(m("Spectra searched").as_deref(), Some("120"));
        assert_eq!(m("Candidate crosslink matches").as_deref(), Some("9"));
        assert_eq!(t.identifications(), Some(3));
    }

    #[test]
    fn notes_skipped_viewer_exports() {
        let mut t = Tracker::new(plan(true, false));
        let log = feed(&mut t, &[
            "--- Sage search ---",
            "[2026-11-02T19:27:00Z WARN  sage_cli::runner::crosslink_runner] no target crosslink matches: crosslinks_proXL.xml was not written (ProXL needs at least one)",
        ]);
        t.finish(State::Finished, None, &log);
        // read_outputs found no results.json here, so both look missing
        t.outputs = vec![
            (PathBuf::from("out/crosslinks_xiView.csv"), 10),
            (PathBuf::from("out/crosslinks.mzid"), 10),
        ];
        let notes = t.missing_exports();
        assert_eq!(notes.len(), 1, "{notes:?}");
        assert!(notes[0].contains("no target crosslinks were found"));

        t.outputs.push((PathBuf::from("out/crosslinks_proXL.xml"), 10));
        assert!(t.missing_exports().is_empty());
        // linear searches have no viewer exports
        let mut linear = Tracker::new(plan(false, false));
        linear.finish(State::Finished, None, &[]);
        assert!(linear.missing_exports().is_empty());
    }

    #[test]
    fn durations_read_well() {
        assert_eq!(tidy_duration("5.5147ms"), "6 ms");
        assert_eq!(tidy_duration("2.5s"), "2.5 s");
        assert_eq!(tidy_duration("310µs"), "< 1 ms");
        assert_eq!(format_ms(125_000), "2m 05s");
    }

    #[test]
    fn diagnoses_common_failures() {
        let p = plan(false, false);
        let cases = [
            ("Sage search ran out of memory: it needed more than the 4.0 GB RAM limit (Files > Resources); raise the limit", "RAM limit"),
            ("memory allocation of 1073741824 bytes failed", "Out of memory"),
            ("Error: Failed to read parameters from `x.json`\nCaused by:\n    missing field `precursor_tol`", "rejected the search settings"),
            ("Error: The system cannot find the file specified. (os error 2)", "not found"),
            ("Error: The process cannot access the file because it is being used by another process. (os error 32)", "in use"),
            ("Error: crosslink backend `gpu` requires building Sage with `--features gpu`", "GPU crosslink backend"),
            ("error: unexpected argument '--threads' found", "does not understand"),
            ("thread 'main' panicked at src/lib.rs:1:1", "internal error"),
            ("something odd", "stopped with an error"),
        ];
        for (text, title) in cases {
            let errors: Vec<String> = text.lines().map(String::from).collect();
            let d = diagnose("Sage search exited with exit code: 1", &errors, &errors, &p, Stage::Search);
            assert!(d.title.contains(title), "{text}: {}", d.title);
            assert!(!d.fixes.is_empty() && !d.evidence.is_empty());
        }
        let d = diagnose("Convert a.raw exited with exit code: 1", &[], &[], &p, Stage::Convert);
        assert!(d.title.contains("Thermo"));
    }

    #[test]
    fn sage_errors_after_the_log_are_collected() {
        let mut t = Tracker::new(plan(false, false));
        let log = feed(&mut t, &[
            "--- Sage search ---",
            "Error: Failed to read parameters from `C:\\x\\sage_search.json`",
            "",
            "Caused by:",
            "    missing field `precursor_tol` at line 1 column 474",
            "=== Failed: Sage search exited with exit code: 1 ===",
        ]);
        t.finish(State::Failed, Some("Sage search exited with exit code: 1"), &log);
        let d = t.diagnosis.as_ref().unwrap();
        assert!(d.title.contains("settings"), "{}", d.title);
        assert!(d.evidence.iter().any(|l| l.contains("missing field")), "{:?}", d.evidence);
    }
}
