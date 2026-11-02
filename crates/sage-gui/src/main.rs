//! `sage-gui`: a desktop window for editing Sage search settings, choosing
//! input files and running `sage.exe`. Every field explains itself on hover.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod runner;
mod settings;

use eframe::egui::{self, Color32, ComboBox, DragValue, Response, RichText, TextEdit, Ui};
use runner::{Job, Status, Step};
use settings::{Backend, ModKind, Modification, Reagent, Settings, TolUnit, Tolerance, Tmt};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Sage search")
            .with_inner_size([1040.0, 780.0])
            .with_min_inner_size([720.0, 480.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Sage search",
        options,
        Box::new(|_cc| Ok(Box::new(SageGui::new()))),
    )
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tab {
    Files,
    Database,
    Modifications,
    Spectra,
    Scoring,
    Crosslink,
}

impl Tab {
    const ALL: [(Tab, &'static str); 6] = [
        (Tab::Files, "Files"),
        (Tab::Database, "Database & enzyme"),
        (Tab::Modifications, "Modifications"),
        (Tab::Spectra, "Tolerances & spectra"),
        (Tab::Scoring, "Scoring & output"),
        (Tab::Crosslink, "Crosslinking"),
    ];
}

struct SageGui {
    s: Settings,
    tab: Tab,
    sage_exe: String,
    converter_exe: String,
    parquet: bool,
    telemetry: bool,
    job: Option<Job>,
    /// Errors and warnings from the last check, shown above the log
    problems: Vec<(bool, String)>,
    show_json: bool,
    last_saved: Option<PathBuf>,
}

/// Where the GUI remembers its last settings and tool paths between runs
fn state_path() -> Option<PathBuf> {
    let base = std::env::var_os("APPDATA").or_else(|| std::env::var_os("HOME"))?;
    Some(PathBuf::from(base).join("sage-gui").join("state.json"))
}

impl SageGui {
    fn new() -> Self {
        let mut app = SageGui {
            s: Settings::default(),
            tab: Tab::Files,
            sage_exe: String::new(),
            converter_exe: String::new(),
            parquet: false,
            telemetry: true,
            job: None,
            problems: Vec::new(),
            show_json: false,
            last_saved: None,
        };
        if let Some(state) = state_path()
            .and_then(|p| std::fs::read_to_string(p).ok())
            .and_then(|s| serde_json::from_str::<serde_json::Value>(&s).ok())
        {
            let text = |k: &str| state.get(k).and_then(|v| v.as_str()).unwrap_or("").to_string();
            app.sage_exe = text("sage_exe");
            app.converter_exe = text("converter_exe");
            app.parquet = state.get("parquet").and_then(|v| v.as_bool()).unwrap_or(false);
            app.telemetry = state.get("telemetry").and_then(|v| v.as_bool()).unwrap_or(true);
            if let Some(Ok(s)) = state.get("config").map(Settings::from_json) {
                app.s = s;
            }
        }
        if !Path::new(&app.sage_exe).is_file() {
            app.sage_exe = runner::find_sage().map(|p| p.display().to_string()).unwrap_or_default();
        }
        if app.converter_exe.is_empty() {
            app.converter_exe = runner::find_msconvert()
                .map(|p| p.display().to_string())
                .unwrap_or_default();
        }
        app
    }

    fn save_state(&self) {
        let Some(path) = state_path() else { return };
        let state = serde_json::json!({
            "sage_exe": self.sage_exe,
            "converter_exe": self.converter_exe,
            "parquet": self.parquet,
            "telemetry": self.telemetry,
            "config": self.s.to_json(&self.s.spectra),
        });
        if let Some(dir) = path.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let _ = std::fs::write(path, serde_json::to_string_pretty(&state).unwrap_or_default());
    }

    fn running(&self) -> bool {
        self.job.as_ref().is_some_and(Job::is_running)
    }

    fn load_config(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .set_title("Load Sage configuration")
            .add_filter("Sage configuration (JSON)", &["json"])
            .pick_file()
        else {
            return;
        };
        let loaded = std::fs::read_to_string(&path)
            .map_err(|e| e.to_string())
            .and_then(|t| serde_json::from_str(&t).map_err(|e| e.to_string()))
            .and_then(|v| Settings::from_json(&v));
        match loaded {
            Ok(s) => {
                self.s = s;
                self.problems = vec![(false, format!("Loaded {}", path.display()))];
            }
            Err(e) => self.problems = vec![(true, format!("Could not load {}: {e}", path.display()))],
        }
    }

    fn save_config(&mut self) {
        let mut dialog = rfd::FileDialog::new()
            .set_title("Save Sage configuration")
            .add_filter("Sage configuration (JSON)", &["json"])
            .set_file_name("sage_search.json");
        if !self.s.output_directory.is_empty() {
            dialog = dialog.set_directory(&self.s.output_directory);
        }
        let Some(path) = dialog.save_file() else { return };
        let text = serde_json::to_string_pretty(&self.s.to_json(&self.s.spectra)).unwrap_or_default();
        match std::fs::write(&path, text) {
            Ok(()) => {
                self.problems = vec![(false, format!("Saved {}", path.display()))];
                self.last_saved = Some(path);
            }
            Err(e) => self.problems = vec![(true, format!("Could not save {}: {e}", path.display()))],
        }
    }

    /// Check the settings, write the configuration into the output folder
    /// and start converting / searching
    fn start_run(&mut self, ctx: &egui::Context) {
        let (mut errors, warnings) = self.s.check();
        if !Path::new(&self.sage_exe).is_file() {
            errors.push("sage.exe not found: set its location under Files > Programs.".into());
        }
        let has_raw = self.s.spectra.iter().any(|f| Settings::is_thermo_raw(f));
        if has_raw && !Path::new(&self.converter_exe).is_file() {
            errors.push(
                "Thermo .raw files need a converter: install ProteoWizard (msconvert) or \
                 ThermoRawFileParser and set it under Files > Programs."
                    .into(),
            );
        }
        self.problems = errors
            .iter()
            .map(|e| (true, e.clone()))
            .chain(warnings.iter().map(|w| (false, w.clone())))
            .collect();
        if !errors.is_empty() {
            return;
        }

        let out = PathBuf::from(&self.s.output_directory);
        if let Err(e) = std::fs::create_dir_all(&out) {
            self.problems.push((true, format!("Could not create {}: {e}", out.display())));
            return;
        }
        let converted_dir = out.join("converted_mzML");
        let mut steps = Vec::new();
        let mut spectra = Vec::new();
        for f in &self.s.spectra {
            if Settings::is_thermo_raw(f) {
                let raw = Path::new(f);
                let _ = std::fs::create_dir_all(&converted_dir);
                steps.push(runner::conversion_step(
                    Path::new(&self.converter_exe),
                    raw,
                    &converted_dir,
                ));
                spectra.push(runner::converted_path(raw, &converted_dir).display().to_string());
            } else {
                spectra.push(f.clone());
            }
        }

        let config_path = out.join("sage_search.json");
        let text = serde_json::to_string_pretty(&self.s.to_json(&spectra)).unwrap_or_default();
        if let Err(e) = std::fs::write(&config_path, text) {
            self.problems.push((true, format!("Could not write {}: {e}", config_path.display())));
            return;
        }

        let mut args = vec![config_path.into_os_string()];
        if self.parquet && !self.s.crosslink.enabled {
            args.push("--parquet".into());
        }
        if !self.telemetry {
            args.push("--disable-telemetry-i-dont-want-to-improve-sage".into());
        }
        steps.push(Step {
            title: "Sage search".into(),
            program: PathBuf::from(&self.sage_exe),
            args,
            skip_if_exists: None,
        });
        self.save_state();
        self.job = Some(Job::start(steps, ctx.clone()));
    }

    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw.dropped_files.iter().filter_map(|f| f.path.clone()).collect()
        });
        for path in dropped {
            let lower = path.to_string_lossy().to_ascii_lowercase();
            let text = path.display().to_string();
            if [".fasta", ".fa", ".faa", ".fas"].iter().any(|e| lower.ends_with(e)) {
                self.s.fasta = text;
            } else if lower.ends_with(".json") {
                // ignored: use "Load settings" so a drop never silently replaces everything
            } else if !self.s.spectra.contains(&text) {
                self.s.spectra.push(text);
            }
        }
    }
}

impl eframe::App for SageGui {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.handle_dropped_files(ctx);
        if ctx.input(|i| i.viewport().close_requested()) {
            if let Some(job) = &self.job {
                job.cancel();
            }
            self.save_state();
        }

        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui
                    .button("Load settings…")
                    .on_hover_text("Open an existing Sage JSON configuration and fill in every field from it.")
                    .clicked()
                {
                    self.load_config();
                }
                if ui
                    .button("Save settings…")
                    .on_hover_text(
                        "Save the current settings as a Sage JSON configuration. It can be run \
                         later from the command line:  sage.exe settings.json",
                    )
                    .clicked()
                {
                    self.save_config();
                }
                if ui
                    .button("Reset to defaults")
                    .on_hover_text("Restore default search settings. Your file selections are kept.")
                    .clicked()
                {
                    let files = (
                        std::mem::take(&mut self.s.fasta),
                        std::mem::take(&mut self.s.spectra),
                        std::mem::take(&mut self.s.output_directory),
                    );
                    self.s = Settings::default();
                    (self.s.fasta, self.s.spectra, self.s.output_directory) = files;
                }
                ui.toggle_value(&mut self.show_json, "Show JSON")
                    .on_hover_text("Show the configuration file that will be passed to sage.exe.");
            });
            ui.add_space(2.0);
            ui.horizontal_wrapped(|ui| {
                for (tab, name) in Tab::ALL {
                    let label = if tab == Tab::Crosslink && self.s.crosslink.enabled {
                        format!("{name} ●")
                    } else {
                        name.to_string()
                    };
                    ui.selectable_value(&mut self.tab, tab, label);
                }
            });
            ui.add_space(2.0);
        });

        egui::TopBottomPanel::bottom("actions").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                let running = self.running();
                let run = ui
                    .add_enabled(!running, egui::Button::new(RichText::new("▶ Run search").strong()))
                    .on_hover_text(
                        "Check the settings, convert any Thermo .raw files to mzML, write \
                         sage_search.json into the output folder and run sage.exe.",
                    );
                if run.clicked() {
                    self.start_run(ctx);
                }
                if ui.add_enabled(running, egui::Button::new("■ Cancel")).clicked() {
                    if let Some(job) = &self.job {
                        job.cancel();
                    }
                }
                let out = self.s.output_directory.clone();
                if ui
                    .add_enabled(Path::new(&out).is_dir(), egui::Button::new("Open results folder"))
                    .on_hover_text("Open the output folder in Explorer.")
                    .clicked()
                {
                    let _ = std::process::Command::new("explorer").arg(&out).spawn();
                }
                if let Some(job) = &self.job {
                    let (text, color) = match job.status() {
                        Status::Running(step) => (format!("Running: {step}"), Color32::from_rgb(70, 130, 220)),
                        Status::Finished => ("Finished".to_string(), Color32::from_rgb(40, 160, 70)),
                        Status::Failed(e) => (format!("Failed: {e}"), Color32::from_rgb(210, 60, 60)),
                        Status::Cancelled => ("Cancelled".to_string(), Color32::GRAY),
                    };
                    if job.is_running() {
                        ui.spinner();
                    }
                    ui.label(RichText::new(text).color(color));
                }
            });
            ui.add_space(4.0);
        });

        if self.job.is_some() || !self.problems.is_empty() {
            egui::TopBottomPanel::bottom("log")
                .resizable(true)
                .default_height(220.0)
                .min_height(80.0)
                .show(ctx, |ui| {
                    ui.add_space(4.0);
                    for (is_error, msg) in &self.problems {
                        let (icon, color) = if *is_error {
                            ("✖", Color32::from_rgb(210, 60, 60))
                        } else {
                            ("ℹ", Color32::from_rgb(200, 150, 30))
                        };
                        ui.label(RichText::new(format!("{icon} {msg}")).color(color));
                    }
                    if let Some(job) = &self.job {
                        ui.separator();
                        let lines = job.log.lock().unwrap().clone();
                        ui.horizontal(|ui| {
                            ui.strong("Log");
                            if ui.small_button("Copy").clicked() {
                                ui.ctx().copy_text(lines.join("\n"));
                            }
                        });
                        egui::ScrollArea::vertical()
                            .auto_shrink([false, false])
                            .stick_to_bottom(true)
                            .show(ui, |ui| {
                                for line in &lines {
                                    ui.label(RichText::new(line).monospace().size(12.0));
                                }
                            });
                    }
                });
        }

        if self.show_json {
            let text = serde_json::to_string_pretty(&self.s.to_json(&self.s.spectra)).unwrap_or_default();
            egui::Window::new("Configuration JSON")
                .open(&mut self.show_json)
                .default_size([520.0, 560.0])
                .show(ctx, |ui| {
                    if ui.button("Copy").clicked() {
                        ui.ctx().copy_text(text.clone());
                    }
                    egui::ScrollArea::both().show(ui, |ui| {
                        ui.label(RichText::new(&text).monospace().size(12.0));
                    });
                });
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                ui.add_space(6.0);
                let running = self.running();
                ui.add_enabled_ui(!running, |ui| match self.tab {
                    Tab::Files => self.files_tab(ui),
                    Tab::Database => database_tab(ui, &mut self.s),
                    Tab::Modifications => modifications_tab(ui, &mut self.s),
                    Tab::Spectra => spectra_tab(ui, &mut self.s),
                    Tab::Scoring => scoring_tab(ui, &mut self.s),
                    Tab::Crosslink => crosslink_tab(ui, &mut self.s),
                });
            });
        });
    }
}

// ---------------------------------------------------------------------------
// Widgets

/// One labelled row in a two-column grid; `tip` shows on hover over both
fn row(ui: &mut Ui, label: &str, tip: &str, add: impl FnOnce(&mut Ui) -> Response) {
    ui.label(label).on_hover_text(tip);
    ui.horizontal(|ui| add(ui).on_hover_text(tip));
    ui.end_row();
}

fn grid(ui: &mut Ui, id: &str, add: impl FnOnce(&mut Ui)) {
    egui::Grid::new(id)
        .num_columns(2)
        .spacing([18.0, 8.0])
        .min_col_width(200.0)
        .show(ui, add);
}

fn heading(ui: &mut Ui, text: &str) {
    ui.add_space(8.0);
    ui.heading(text);
    ui.add_space(2.0);
}

fn note(ui: &mut Ui, text: &str) {
    ui.label(RichText::new(text).weak());
    ui.add_space(4.0);
}

fn num<T: egui::emath::Numeric>(ui: &mut Ui, v: &mut T, range: RangeInclusive<T>, speed: f64) -> Response {
    ui.add(DragValue::new(v).range(range).speed(speed))
}

fn mass(ui: &mut Ui, v: &mut f64) -> Response {
    ui.add(DragValue::new(v).speed(0.01).max_decimals(6).suffix(" Da"))
}

/// A value Sage can choose itself: unchecked writes nothing to the JSON
fn optional<T: egui::emath::Numeric>(
    ui: &mut Ui,
    v: &mut Option<T>,
    fallback: T,
    auto: &str,
    range: RangeInclusive<T>,
) -> Response {
    let mut on = v.is_some();
    let mut r = ui.checkbox(&mut on, "set");
    if on != v.is_some() {
        *v = on.then_some(fallback);
    }
    match v {
        Some(x) => r |= ui.add(DragValue::new(x).range(range)),
        None => r |= ui.label(RichText::new(auto).weak()),
    }
    r
}

fn tolerance(ui: &mut Ui, id: &str, t: &mut Tolerance) -> Response {
    let speed = if t.unit == TolUnit::Ppm { 0.5 } else { 0.005 };
    let mut r = ui.add(DragValue::new(&mut t.lo).speed(speed).max_decimals(4).prefix("from "));
    r |= ui.add(DragValue::new(&mut t.hi).speed(speed).max_decimals(4).prefix("to "));
    r |= ComboBox::from_id_salt(id)
        .width(70.0)
        .selected_text(t.unit.key())
        .show_ui(ui, |ui| {
            ui.selectable_value(&mut t.unit, TolUnit::Ppm, "ppm");
            ui.selectable_value(&mut t.unit, TolUnit::Da, "da");
        })
        .response;
    r
}

fn path_field(ui: &mut Ui, value: &mut String, hint: &str, pick: impl FnOnce() -> Option<PathBuf>) -> Response {
    let mut r = ui.add(TextEdit::singleline(value).hint_text(hint).desired_width(520.0));
    let browse = ui.button("Browse…");
    if browse.clicked() {
        if let Some(p) = pick() {
            *value = p.display().to_string();
        }
    }
    r |= browse;
    r
}

// ---------------------------------------------------------------------------
// Tabs

const SPECTRA_TIP: &str = "Mass spectrometry runs to search.\n\
    • mzML / mzML.gz and MGF / MGF.gz: read directly\n\
    • Bruker timsTOF .d: add the .d FOLDER (\"Add Bruker .d folder…\")\n\
    • Thermo .raw: converted to mzML with msconvert or ThermoRawFileParser before the \
    search, into <output>\\converted_mzML (already converted files are reused)\n\
    You can also drag files onto this window.";

impl SageGui {
    fn files_tab(&mut self, ui: &mut Ui) {
        heading(ui, "Protein database");
        grid(ui, "fasta", |ui| {
            row(
                ui,
                "FASTA file",
                "Protein sequences to search against, e.g. C:\\db\\human_uniprot_2026.fasta.\n\
                 Use a database WITHOUT decoys: Sage generates reversed decoys itself \
                 (Database tab > Generate decoys).",
                |ui| {
                    path_field(ui, &mut self.s.fasta, "C:\\db\\proteome.fasta", || {
                        rfd::FileDialog::new()
                            .add_filter("FASTA", &["fasta", "fa", "faa", "fas"])
                            .add_filter("All files", &["*"])
                            .pick_file()
                    })
                },
            );
        });

        heading(ui, "Spectrum files");
        ui.horizontal(|ui| {
            if ui.button("Add files…").on_hover_text(SPECTRA_TIP).clicked() {
                if let Some(files) = rfd::FileDialog::new()
                    .add_filter("Spectra (mzML, MGF, Thermo RAW)", &["mzML", "mzml", "gz", "mgf", "raw"])
                    .add_filter("All files", &["*"])
                    .pick_files()
                {
                    for f in files {
                        let f = f.display().to_string();
                        if !self.s.spectra.contains(&f) {
                            self.s.spectra.push(f);
                        }
                    }
                }
            }
            if ui
                .button("Add Bruker .d folder…")
                .on_hover_text("Select one or more Bruker timsTOF acquisition folders ending in .d (e.g. sample_01.d).")
                .clicked()
            {
                if let Some(dirs) = rfd::FileDialog::new().pick_folders() {
                    for d in dirs {
                        let d = d.display().to_string();
                        if !self.s.spectra.contains(&d) {
                            self.s.spectra.push(d);
                        }
                    }
                }
            }
            if ui.button("Clear list").clicked() {
                self.s.spectra.clear();
            }
            ui.label(format!("{} file(s)", self.s.spectra.len()));
        });
        ui.add_space(4.0);
        if self.s.spectra.is_empty() {
            note(ui, "No files yet. Use the buttons above or drag files and .d folders onto the window.");
        }
        let mut remove = None;
        egui::Grid::new("spectra").num_columns(3).striped(true).show(ui, |ui| {
            for (i, f) in self.s.spectra.iter().enumerate() {
                if ui.small_button("✖").on_hover_text("Remove from the list").clicked() {
                    remove = Some(i);
                }
                let (kind, tip) = spectrum_kind(f);
                ui.label(RichText::new(kind).strong()).on_hover_text(tip);
                ui.label(f);
                ui.end_row();
            }
        });
        if let Some(i) = remove {
            self.s.spectra.remove(i);
        }

        heading(ui, "Results");
        grid(ui, "output", |ui| {
            row(
                ui,
                "Output folder",
                "Folder where Sage writes results.sage.tsv (PSMs), results.json (the full settings \
                 used), quantification tables and, for crosslink searches, crosslinks.sage.tsv.\n\
                 Created if it does not exist. Example: C:\\results\\2026-11-01_hela",
                |ui| {
                    path_field(ui, &mut self.s.output_directory, "C:\\results\\my_search", || {
                        rfd::FileDialog::new().pick_folder()
                    })
                },
            );
            let xl = self.s.crosslink.enabled;
            row(
                ui,
                "Parquet output",
                "Also write results as Apache Parquet (results.sage.parquet, lfq.parquet) for fast \
                 loading in Python/R. Not available for crosslink searches.",
                |ui| ui.add_enabled(!xl, egui::Checkbox::new(&mut self.parquet, "write .parquet files")),
            );
        });

        heading(ui, "Programs");
        grid(ui, "programs", |ui| {
            row(
                ui,
                "sage.exe",
                "The Sage search engine. Found automatically when sage.exe sits next to sage-gui.exe \
                 (both are built into target\\release).",
                |ui| {
                    path_field(ui, &mut self.sage_exe, "C:\\Tools\\sage\\sage.exe", || {
                        rfd::FileDialog::new().add_filter("Programs", &["exe"]).pick_file()
                    })
                },
            );
            row(
                ui,
                "RAW converter",
                "Only needed for Thermo .raw files, which Sage cannot read directly.\n\
                 • ProteoWizard msconvert.exe (found automatically if installed), run with vendor \
                 peak picking: --mzML --zlib --filter \"peakPicking vendor msLevel=1-\"\n\
                 • or ThermoRawFileParser.exe, run with -f 2 (indexed mzML)",
                |ui| {
                    let mut r = path_field(ui, &mut self.converter_exe, "C:\\Program Files\\ProteoWizard\\…\\msconvert.exe", || {
                        rfd::FileDialog::new().add_filter("Programs", &["exe"]).pick_file()
                    });
                    let find = ui.button("Find msconvert");
                    if find.clicked() {
                        if let Some(p) = runner::find_msconvert() {
                            self.converter_exe = p.display().to_string();
                        }
                    }
                    r |= find;
                    r
                },
            );
            row(
                ui,
                "Usage statistics",
                "Sage sends anonymous usage statistics (run time, number of spectra, settings) to its \
                 developers by default. Untick to run with --disable-telemetry-i-dont-want-to-improve-sage.",
                |ui| ui.checkbox(&mut self.telemetry, "send anonymous usage statistics to the Sage developers"),
            );
        });
    }
}

fn spectrum_kind(path: &str) -> (&'static str, &'static str) {
    let lower = path.to_ascii_lowercase();
    let trimmed = lower.trim_end_matches(['\\', '/']);
    if trimmed.ends_with(".d") {
        ("Bruker .d", "Bruker timsTOF folder, read directly")
    } else if Settings::is_thermo_raw(path) {
        ("Thermo RAW", "Converted to mzML before the search")
    } else if lower.ends_with(".mzml") || lower.ends_with(".mzml.gz") {
        ("mzML", "Read directly")
    } else if lower.ends_with(".mgf") || lower.ends_with(".mgf.gz") {
        ("MGF", "Read directly")
    } else if Path::new(path).is_dir() {
        ("Folder?", "Sage reads folders only if they are Bruker .d acquisitions")
    } else {
        ("Unknown", "Sage reads mzML, MGF and Bruker .d data; this file may not be readable")
    }
}

fn database_tab(ui: &mut Ui, s: &mut Settings) {
    heading(ui, "Enzyme");
    note(ui, "Tick one or more enzymes. Several enzymes are searched as a co-digestion: peptides are cut \
        wherever any of them cleaves (e.g. trypsin + Asp-N, trypsin + Lys-C).");
    ui.checkbox(&mut s.nonspecific, "Non-specific search (no enzyme)")
        .on_hover_text("Consider every peptide bond as a cleavage site, e.g. for immunopeptidomics or \
            endogenous peptides. The database becomes very large: use a narrow peptide length range \
            (e.g. 8 to 12) and expect long searches. Missed cleavages do not apply.");
    ui.add_space(4.0);
    ui.add_enabled_ui(!s.nonspecific, |ui| {
        egui::Grid::new("proteases").num_columns(3).spacing([24.0, 6.0]).show(ui, |ui| {
            for (i, p) in settings::PROTEASES.iter().enumerate() {
                let mut on = s.enzymes.iter().any(|n| n == p.name);
                let side = if p.c_terminal { "after" } else { "before" };
                let restrict = if p.restrict.is_empty() {
                    String::new()
                } else {
                    format!(", not before {}", p.restrict)
                };
                let tip = format!("{}\nRule: cut {side} {}{restrict}.", p.help, p.cleave_at);
                if ui.checkbox(&mut on, p.label).on_hover_text(tip).changed() {
                    if on {
                        s.enzymes.push(p.name.into());
                        // keep the list in the displayed order
                        s.enzymes.sort_by_key(|n| settings::PROTEASES.iter().position(|x| x.name == n));
                    } else {
                        s.enzymes.retain(|n| n != p.name);
                    }
                }
                if i % 3 == 2 {
                    ui.end_row();
                }
            }
        });

        ui.add_space(6.0);
        let mut remove = None;
        if !s.custom_enzymes.is_empty() {
            egui::Grid::new("custom_enzymes").num_columns(4).spacing([12.0, 6.0]).show(ui, |ui| {
                ui.strong("Cleave at").on_hover_text("Residues at which this enzyme cuts, e.g. DE for Glu-C in phosphate buffer.");
                ui.strong("Not before").on_hover_text("Do not cut when the next residue is one of these, e.g. P. Leave empty for none.");
                ui.strong("Side").on_hover_text("After the residue (C-terminal, like trypsin) or before it (N-terminal, like Asp-N).");
                ui.label("");
                ui.end_row();
                for (i, rule) in s.custom_enzymes.iter_mut().enumerate() {
                    ui.add(TextEdit::singleline(&mut rule.cleave_at).desired_width(90.0))
                        .on_hover_text("One-letter codes, e.g. DE");
                    ui.add(TextEdit::singleline(&mut rule.restrict).desired_width(70.0))
                        .on_hover_text("One-letter codes, e.g. P");
                    ComboBox::from_id_salt(("side", i))
                        .selected_text(if rule.c_terminal { "after (C-term)" } else { "before (N-term)" })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut rule.c_terminal, true, "after (C-term)");
                            ui.selectable_value(&mut rule.c_terminal, false, "before (N-term)");
                        });
                    if ui.small_button("✖").on_hover_text("Remove this rule").clicked() {
                        remove = Some(i);
                    }
                    ui.end_row();
                }
            });
        }
        if let Some(i) = remove {
            s.custom_enzymes.remove(i);
        }
        if ui
            .button("Add custom enzyme rule")
            .on_hover_text("For an enzyme not in the list. Example: Glu-C in phosphate buffer = cleave at DE, \
                not before P, after. Pepsin-like: cleave at FL, after.")
            .clicked()
        {
            s.custom_enzymes.push(settings::EnzymeRule {
                cleave_at: "DE".into(),
                restrict: "P".into(),
                c_terminal: true,
            });
        }
    });
    ui.add_space(6.0);

    grid(ui, "enzyme", |ui| {
        row(ui, "Missed cleavages", "Maximum number of internal cleavage sites a peptide may contain.\n\
            Typical: 2 for trypsin. Crosslink searches need at least 1 (a crosslinked lysine is not cleaved).\n\
            Sage's own default when the field is absent is 0.", |ui| num(ui, &mut s.missed_cleavages, 0..=10, 0.1));
        row(ui, "Peptide length", "Shortest and longest peptide (amino acids) to search. Typical: 7 to 50.", |ui| {
            let r = ui.add(DragValue::new(&mut s.min_len).range(1..=100).prefix("min "));
            r | ui.add(DragValue::new(&mut s.max_len).range(1..=200).prefix("max "))
        });
        row(ui, "Semi-enzymatic", "Also allow peptides with only one enzymatic end. Finds in-source fragments \
            and degradation products, but the database grows ~10x and searches get slower.",
            |ui| ui.checkbox(&mut s.semi_enzymatic, "semi-specific digestion"));
    });

    heading(ui, "Peptides & fragments");
    grid(ui, "peptides", |ui| {
        row(ui, "Minimum peptide mass", "Lightest peptide (monoisotopic, Da) to search.\n\
            Unset = Sage's choice: 500 Da, or 350 Da for crosslink searches so small partner peptides are kept.",
            |ui| optional(ui, &mut s.peptide_min_mass, 500.0, "automatic (500 Da; 350 Da for crosslinks)", 0.0..=10000.0));
        row(ui, "Maximum peptide mass", "Heaviest peptide (monoisotopic, Da) to search. Default 5000.",
            |ui| ui.add(DragValue::new(&mut s.peptide_max_mass).range(100.0..=20000.0).speed(10.0).suffix(" Da")));
        row(ui, "Fragment ion types", "Ion series that are matched in MS2 spectra.\n\
            • CID/HCD: b and y (default)\n• ETD/ECD: c and z (or z_radical)\n• EThcD: b, y, c, z",
            |ui| {
                let mut r = ui.label("");
                for (i, name) in settings::ION_KINDS.iter().enumerate() {
                    r |= ui.checkbox(&mut s.ion_kinds[i], *name);
                }
                r
            });
        row(ui, "Minimum ion index", "Skip the first N ions of each series: 2 (default) ignores b1/y1, \
            which are rarely observed and match noise easily.", |ui| num(ui, &mut s.min_ion_index, 0..=10, 0.1));
    });

    heading(ui, "Decoys & indexing");
    grid(ui, "decoys", |ui| {
        row(ui, "Generate decoys", "Add reversed decoy sequences for target-decoy FDR estimation. \
            Turn off only if your FASTA already contains decoys (then set the decoy tag to match them).",
            |ui| ui.checkbox(&mut s.generate_decoys, "generate reversed decoys"));
        row(ui, "Decoy tag", "Prefix marking decoy protein accessions. Default rev_ (e.g. rev_sp|P12345|…).",
            |ui| ui.add(TextEdit::singleline(&mut s.decoy_tag).desired_width(120.0)));
        row(ui, "Bucket size", "Internal fragment-index tuning (rounded up to a power of two). \
            Default 8192; larger values (e.g. 32768) can be faster for wide precursor windows / open searches.",
            |ui| num(ui, &mut s.bucket_size, 1024..=1_048_576, 256.0));
        row(ui, "Prefilter database", "For very large databases (e.g. metaproteomics): first find peptides \
            that match at all, then search only those. Lowers memory use. Not allowed for crosslink searches.",
            |ui| ui.checkbox(&mut s.prefilter, "prefilter in chunks"));
        if s.prefilter {
            row(ui, "Prefilter chunk size", "Number of protein sequences per prefilter chunk. 0 = automatic.",
                |ui| num(ui, &mut s.prefilter_chunk_size, 0..=10_000_000, 1000.0));
            row(ui, "Low-memory prefilter", "Use as little memory as possible during prefiltering, at some cost in speed.",
                |ui| ui.checkbox(&mut s.prefilter_low_memory, "low memory"));
        }
    });
}

const TARGET_TIP: &str = "Where the modification can occur:\n\
    • one amino acid: C, M, S, K, …\n\
    • ^ peptide N-terminus (any residue), ^Q only when the N-terminal residue is Q\n\
    • $ peptide C-terminus, $K\n\
    • [ protein N-terminus, [M\n\
    • ] protein C-terminus";

fn modifications_tab(ui: &mut Ui, s: &mut Settings) {
    heading(ui, "Modifications");
    note(ui, "Fixed modifications are always present (e.g. carbamidomethyl C after iodoacetamide). \
        Variable modifications may or may not be present; each one increases search time.");

    ui.horizontal(|ui| {
        ComboBox::from_id_salt("add_preset")
            .selected_text("Add common modification…")
            .width(280.0)
            .show_ui(ui, |ui| {
                for p in settings::PRESETS {
                    let kind = if p.kind == ModKind::Fixed { "fixed" } else { "variable" };
                    if ui
                        .selectable_label(false, p.label)
                        .on_hover_text(format!("{:+.6} Da on {}, {} by default ({})", p.mass, p.target, kind, p.name))
                        .clicked()
                    {
                        s.mods.push(p.into());
                    }
                }
            })
            .response
            .on_hover_text("Pick a frequently used modification; mass and UniMod name are filled in.");
        if ui
            .button("Add custom modification")
            .on_hover_text("Add an empty row and enter the target and mass yourself.")
            .clicked()
        {
            s.mods.push(Modification {
                target: "K".into(),
                mass: 0.0,
                name: String::new(),
                kind: ModKind::Variable,
                max_count: 0,
            });
        }
    });
    ui.add_space(6.0);

    let mut remove = None;
    egui::Grid::new("mods").num_columns(6).spacing([12.0, 6.0]).striped(true).show(ui, |ui| {
        ui.strong("Type").on_hover_text("Fixed: always applied. Variable: searched with and without.");
        ui.strong("Target").on_hover_text(TARGET_TIP);
        ui.strong("Mass shift").on_hover_text("Monoisotopic mass added to the residue, in Da. Negative values allowed (e.g. -17.026549 for pyro-Glu from Q).");
        ui.strong("Max per peptide").on_hover_text("Variable modifications only: how often this modification may occur on one peptide. 0 = only limited by 'Max variable mods'.");
        ui.strong("Name (optional)").on_hover_text("ProForma name written to the results, e.g. UNIMOD:35 for oxidation. Leave empty to report the mass.");
        ui.label("");
        ui.end_row();

        for (i, m) in s.mods.iter_mut().enumerate() {
            ComboBox::from_id_salt(("kind", i))
                .width(90.0)
                .selected_text(if m.kind == ModKind::Fixed { "Fixed" } else { "Variable" })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut m.kind, ModKind::Fixed, "Fixed");
                    ui.selectable_value(&mut m.kind, ModKind::Variable, "Variable");
                });
            let bad = settings::check_target(&m.target).is_err();
            let target = ui.add(TextEdit::singleline(&mut m.target).desired_width(50.0).text_color_opt(bad.then_some(Color32::RED)));
            target.on_hover_text(TARGET_TIP);
            mass(ui, &mut m.mass).on_hover_text("Monoisotopic mass shift in Da, e.g. 15.994915 for oxidation.");
            ui.add_enabled(m.kind == ModKind::Variable, DragValue::new(&mut m.max_count).range(0..=10))
                .on_hover_text("0 = no per-modification limit.");
            ui.add(TextEdit::singleline(&mut m.name).hint_text("UNIMOD:…").desired_width(110.0))
                .on_hover_text("e.g. UNIMOD:4 (carbamidomethyl), UNIMOD:35 (oxidation), UNIMOD:21 (phospho)");
            if ui.small_button("✖").on_hover_text("Remove this modification").clicked() {
                remove = Some(i);
            }
            ui.end_row();
        }
    });
    if let Some(i) = remove {
        s.mods.remove(i);
    }

    heading(ui, "Limits");
    grid(ui, "modlimits", |ui| {
        row(ui, "Max variable mods", "Maximum variable modifications on one peptide (all kinds together). \
            Default 2; 3 is common for phospho studies. Higher values grow the search space quickly.",
            |ui| num(ui, &mut s.max_variable_mods, 1..=10, 0.1));
        row(ui, "Max combinations", "Hard cap on modified forms generated per peptide (including the unmodified \
            one). Unset = no cap. Useful against combinatorial explosion, e.g. 64.",
            |ui| optional(ui, &mut s.max_combinations, 64, "no limit", 1..=100_000));
    });
}

fn spectra_tab(ui: &mut Ui, s: &mut Settings) {
    heading(ui, "Mass tolerances");
    note(ui, "Tolerances are applied to the experimental mass: [-10, 10] ppm accepts peptides within 10 ppm either way.");
    grid(ui, "tol", |ui| {
        row(ui, "Precursor tolerance", "Allowed difference between measured and calculated precursor mass.\n\
            • Orbitrap / timsTOF: -10 to 10 ppm (or -20 to 20)\n\
            • Open search for unknown modifications: -500 to 100 da (raise the bucket size too)",
            |ui| tolerance(ui, "ptol", &mut s.precursor_tol));
        row(ui, "Fragment tolerance", "Allowed difference for fragment ions.\n\
            • Orbitrap HCD (high resolution MS2): -10 to 10 ppm (or ±20 ppm)\n\
            • timsTOF: -20 to 20 ppm\n\
            • Ion trap (low resolution MS2): -0.5 to 0.5 da",
            |ui| tolerance(ui, "ftol", &mut s.fragment_tol));
        row(ui, "Isotope errors", "Also accept precursors picked on the wrong isotope peak (C13 errors).\n\
            Example: -1 to 3 tries M-1 … M+3. Default 0 to 0 (off). Increases search time.",
            |ui| {
                let r = ui.add(DragValue::new(&mut s.isotope_errors.0).range(-5..=5).prefix("from "));
                r | ui.add(DragValue::new(&mut s.isotope_errors.1).range(-5..=5).prefix("to "))
            });
    });

    heading(ui, "Charge states");
    grid(ui, "charges", |ui| {
        row(ui, "Precursor charge", "Charge states to try when the file does not report one. Default 2 to 4; \
            crosslinked peptides are often 3+ to 6+.", |ui| {
            let r = ui.add(DragValue::new(&mut s.precursor_charge.0).range(1..=10).prefix("from "));
            r | ui.add(DragValue::new(&mut s.precursor_charge.1).range(1..=10).prefix("to "))
        });
        row(ui, "Override charge", "Ignore the charge state in the file and always try the range above.",
            |ui| ui.checkbox(&mut s.override_precursor_charge, "always use the charge range"));
        row(ui, "Max fragment charge", "Highest fragment ion charge to consider. Unset = precursor charge - 1. \
            Example: 1 for low charge peptides, 2 for most HCD data.",
            |ui| optional(ui, &mut s.max_fragment_charge, 2, "precursor charge - 1", 1..=6));
    });

    heading(ui, "Spectrum processing");
    grid(ui, "spectrum", |ui| {
        row(ui, "Deisotope", "Collapse isotope peaks and assign charges to fragment peaks. Recommended for \
            high-resolution MS2 (Orbitrap, timsTOF) and for crosslink searches.",
            |ui| ui.checkbox(&mut s.deisotope, "deisotope MS2 spectra"));
        row(ui, "Chimeric spectra", "Search for a second co-fragmented peptide in each spectrum. \
            Ignored in crosslink searches.", |ui| ui.checkbox(&mut s.chimera, "report chimeric (co-isolated) peptides"));
        row(ui, "Wide window (DIA)", "Data-independent acquisition mode: ignore the precursor mass and use the \
            isolation window. Not supported for crosslink searches.",
            |ui| ui.checkbox(&mut s.wide_window, "DIA / wide-window search"));
        row(ui, "Min peaks", "Skip spectra with fewer peaks than this. Default 15.",
            |ui| num(ui, &mut s.min_peaks, 0..=1000, 1.0));
        row(ui, "Max peaks", "Keep only the N most intense peaks per spectrum. Unset = Sage's choice: 150, \
            or 300 for MS-cleavable crosslinkers (their doublet peaks are often weak).",
            |ui| optional(ui, &mut s.max_peaks, 150, "automatic (150; 300 for DSSO/DSBU)", 10..=5000));
        row(ui, "Min matched peaks", "A PSM needs at least this many matched b/y (etc.) ions. Default 4.",
            |ui| num(ui, &mut s.min_matched_peaks, 1..=50, 0.1));
    });
}

fn scoring_tab(ui: &mut Ui, s: &mut Settings) {
    heading(ui, "Scoring & FDR");
    grid(ui, "scoring", |ui| {
        row(ui, "Score function", "Sage hyperscore (default) or the OpenMS hyperscore variant.", |ui| {
            let r = ui.radio_value(&mut s.openms_score, false, "Sage hyperscore");
            r | ui.radio_value(&mut s.openms_score, true, "OpenMS hyperscore")
        });
        row(ui, "PSMs per spectrum", "Number of peptide matches reported per spectrum. Default 1; use \
            e.g. 5 when rescoring with Percolator or mokapot.", |ui| num(ui, &mut s.report_psms, 1..=100, 0.1));
        row(ui, "Predict retention time", "Fit a retention time model and use the RT error as a feature for \
            rescoring (recommended). Required for LFQ.", |ui| ui.checkbox(&mut s.predict_rt, "predict RT"));
        row(ui, "Protein grouping", "Group proteins that share the same peptide evidence (parsimony).",
            |ui| ui.checkbox(&mut s.protein_grouping, "group proteins"));
        row(ui, "Grouping peptide FDR", "Peptide q-value cutoff used when building protein groups. Default 0.01 (1%).",
            |ui| ui.add(DragValue::new(&mut s.protein_grouping_peptide_fdr).range(0.0..=1.0).speed(0.001).max_decimals(3)));
    });

    heading(ui, "Quantification");
    if s.crosslink.enabled {
        note(ui, "Crosslink search: label-free intensities are written per crosslinked peptide pair to \
            crosslinks.lfq.tsv; TMT reporter intensities are added as columns of crosslinks.sage.tsv.");
    }
    grid(ui, "quant", |ui| {
        row(ui, "Label-free (LFQ)", "Extract MS1 peak areas across all runs. Needs MS1 spectra in the files.\n\
            • Linear search: per peptide, in lfq.tsv\n\
            • Crosslink search: per crosslinked peptide pair (link-site isomers can't be told apart in MS1), \
            in crosslinks.lfq.tsv. Runs are aligned on crosslinks identified in several runs.",
            |ui| ui.checkbox(&mut s.lfq, "label-free quantification"));
        if s.lfq {
            row(ui, "LFQ q-value cutoff", "Only peptides (or crosslinks, by crosslink q-value) at or below this \
                q-value are quantified. Default 0.01 (1% FDR).",
                |ui| ui.add(DragValue::new(&mut s.lfq_q).range(0.0..=1.0).speed(0.001).max_decimals(3)));
            row(ui, "MS1 tolerance", "Mass tolerance for tracing precursor isotopes in MS1, in ppm. Default 5; \
                10 for older or lower-resolution instruments.",
                |ui| ui.add(DragValue::new(&mut s.lfq_ppm).range(0.5..=50.0).speed(0.1).suffix(" ppm")));
            row(ui, "Min spectral angle", "Minimum agreement (0 to 1) between the observed and theoretical isotope \
                envelope for a peak to be integrated. Default 0.7.",
                |ui| ui.add(DragValue::new(&mut s.lfq_spectral_angle).range(0.0..=1.0).speed(0.01).max_decimals(2)));
            row(ui, "Charge states", "On: sum all charge states into one intensity per peptide / crosslink \
                (default). Off: report each charge state separately.",
                |ui| ui.checkbox(&mut s.lfq_combine_charges, "combine charge states"));
        }
        row(ui, "TMT", "Isobaric labelling kit used. Reporter ions are integrated into tmt.tsv \
            (crosslink searches: columns tmt_1 … of crosslinks.sage.tsv).\n\
            Add the matching TMT modifications on the Modifications tab: Fixed on the peptide N-terminus (^). \
            On K: Fixed for linear searches, but Variable for crosslink searches, since a crosslinked lysine \
            carries no label and a fixed label would block every lysine link.",
            |ui| {
                ComboBox::from_id_salt("tmt")
                    .selected_text(s.tmt.label())
                    .show_ui(ui, |ui| {
                        for t in Tmt::ALL {
                            ui.selectable_value(&mut s.tmt, t, t.label());
                        }
                    })
                    .response
            });
        if s.tmt != Tmt::None {
            row(ui, "TMT MS level", "MS level of the reporter ions: 2 for MS2 methods, 3 for SPS-MS3 (default 3).",
                |ui| num(ui, &mut s.tmt_level, 2..=3, 0.05));
            row(ui, "Signal-to-noise", "Report reporter signal-to-noise instead of intensity (Thermo files with noise data).",
                |ui| ui.checkbox(&mut s.tmt_sn, "use S/N"));
        }
    });

    heading(ui, "Extra output");
    grid(ui, "extra", |ui| {
        row(ui, "Percolator input", "Also write results.sage.pin for rescoring with Percolator or mokapot.",
            |ui| ui.checkbox(&mut s.write_pin, "write .pin file"));
        row(ui, "HTML report", "Write an HTML summary report with plots.", |ui| ui.checkbox(&mut s.write_report, "write report"));
        row(ui, "Annotate matches", "Write matched_fragments.sage.tsv listing every matched fragment ion \
            (useful for plotting annotated spectra; large file).",
            |ui| ui.checkbox(&mut s.annotate_matches, "write matched fragment ions"));
    });
}

fn crosslink_tab(ui: &mut Ui, s: &mut Settings) {
    heading(ui, "Crosslinking (XL-MS)");
    ui.checkbox(&mut s.crosslink.enabled, RichText::new("Search for crosslinked peptides").strong())
        .on_hover_text("Turns this into a crosslink search: results go to crosslinks.sage.tsv. \
            Label-free and TMT quantification work (Scoring & output tab); chimera, DIA and parquet \
            output are not used in this mode. \
            Details: CROSSLINKING.md");
    note(ui, "Recommended database settings: 2 missed cleavages, deisotoping on, precursor charge 3 to 6.");

    let xl = &mut s.crosslink;
    ui.add_enabled_ui(xl.enabled, |ui| {
        grid(ui, "xl", |ui| {
            row(ui, "Reagent", "Crosslinker used.\n\
                • DSS / BS3: 138.06808 Da, non-cleavable\n• DSG / BS2G: 96.02113 Da, non-cleavable\n\
                • DSSO: 158.00376 Da, MS-cleavable (stubs 54.011 / 103.993 / 85.983)\n\
                • DSBU: 196.08479 Da, MS-cleavable (stubs 85.053 / 111.032)\n\
                • EDC / DMTMM: zero-length (-18.01056 Da, the water released when the amide forms); \
                link an amine (K, protein N-term) to a carboxyl (D, E, protein C-term)\n\
                MS-cleavable reagents use a fast doublet search that scales to whole proteomes.",
                |ui| {
                    let before = xl.reagent;
                    let r = ComboBox::from_id_salt("reagent")
                        .width(260.0)
                        .selected_text(xl.reagent.label())
                        .show_ui(ui, |ui| {
                            for r in Reagent::ALL {
                                ui.selectable_value(&mut xl.reagent, r, r.label());
                            }
                        })
                        .response;
                    if xl.reagent != before {
                        // Presets bring their chemistry; leaving one restores amine-only sites
                        if let Some((s1, s2)) = xl.reagent.sites() {
                            xl.two_sites = true;
                            xl.site1 = s1;
                            xl.site2 = s2;
                        } else if before.sites().is_some() && xl.reagent != Reagent::Custom {
                            xl.two_sites = false;
                        }
                    }
                    r
                });
            if xl.reagent == Reagent::Custom {
                row(ui, "Linker mass", "Mass the linker adds to the pair of peptides (Da). Example: 138.06808 for DSS; \
                    -18.010565 for a zero-length coupling (one water is lost).",
                    |ui| mass(ui, &mut xl.custom_mass));
                row(ui, "MS-cleavable", "The custom linker breaks in the collision cell and leaves stubs on each peptide.",
                    |ui| ui.checkbox(&mut xl.custom_cleavable, "cleavable linker with stub masses"));
                if xl.custom_cleavable {
                    row(ui, "Stub masses", "Short and long stub remnants. short + long must equal the linker mass \
                        (DSSO: 54.01056 + 103.99320 = 158.00376).", |ui| {
                        let r = mass(ui, &mut xl.stub_short);
                        r | mass(ui, &mut xl.stub_long)
                    });
                    row(ui, "Additional stubs", "Other observed stub masses, comma separated. Example (DSSO sulfenic acid): 85.98264",
                        |ui| ui.add(TextEdit::singleline(&mut xl.stub_alt).hint_text("85.98264").desired_width(200.0)));
                }
            }
            row(ui, "Reactive groups", "Same: both ends of the linker react with the same residues \
                (homobifunctional NHS esters such as DSS, BS3, DSSO, DSBU).\n\
                Different: each end reacts with another group (heterobifunctional or zero-length reagents \
                such as EDC and DMTMM, which join an amine to a carboxyl). Every link then joins a \
                Site 1 group to a Site 2 group; K-K or D-E links are never reported.", |ui| {
                let r = ui.radio_value(&mut xl.two_sites, false, "same at both ends");
                r | ui.radio_value(&mut xl.two_sites, true, "different at each end (Site 1 / Site 2)")
            });
            if xl.two_sites {
                site_row(ui, "Site 1", "first", &mut xl.site1);
                site_row(ui, "Site 2", "second", &mut xl.site2);
            } else {
                row(ui, "Reactive residues", "Amino acids the linker reacts with. NHS esters (DSS, BS3, DSSO, DSBU): K, \
                    or KSTY to include the weaker side reactions.",
                    |ui| ui.add(TextEdit::singleline(&mut xl.residues).desired_width(100.0)));
                row(ui, "Protein N-terminus", "Allow links to the protein N-terminal amine (NHS esters react with it).",
                    |ui| ui.checkbox(&mut xl.protein_nterm, "link protein N-termini"));
            }
            row(ui, "Peptide C-terminal residue", "Allow a link on a peptide's last residue when it is not the \
                protein C-terminus. Normally off: a crosslinked lysine blocks trypsin cleavage.",
                |ui| ui.checkbox(&mut xl.allow_cterm_residue, "allow"));
            let cleavable = cleavable_selected(xl);
            row(ui, "Loop-links", "Also search links with both ends on the same peptide. Not searched for \
                MS-cleavable reagents.", |ui| ui.add_enabled(!cleavable, egui::Checkbox::new(&mut xl.loop_links, "search loop-links")));
            if cleavable {
                row(ui, "Complete doublets", "Require both stub peaks for BOTH peptides. Stricter, fewer hits. \
                    Default: both stubs of one peptide and at least one of the other.",
                    |ui| ui.checkbox(&mut xl.require_complete_doublets, "require complete doublets"));
            }
            row(ui, "Min peaks per chain", "Each of the two peptides needs at least this many matched fragments. Default 3.",
                |ui| num(ui, &mut xl.min_chain_peaks, 1..=20, 0.1));
        });

        heading(ui, "Non-cleavable search engine");
        if cleavable_selected(xl) {
            note(ui, "MS-cleavable reagents always use the CPU doublet search; these settings are ignored.");
        }
        grid(ui, "xl_backend", |ui| {
            row(ui, "Backend", "How non-cleavable crosslinks are searched:\n\
                • cpu (default): seed search + partner lookup; validated, fine for complexes and small databases\n\
                • reflect: chain reflection on the CPU; reference for the GPU, slow\n\
                • gpu: chain reflection on the graphics card, for proteome-scale searches (needs a sage.exe built with --features gpu)\n\
                • auto: gpu when a real GPU is available, otherwise cpu",
                |ui| {
                    ComboBox::from_id_salt("backend")
                        .selected_text(xl.backend.key())
                        .show_ui(ui, |ui| {
                            for b in Backend::ALL {
                                ui.selectable_value(&mut xl.backend, b, b.key());
                            }
                        })
                        .response
                });
            row(ui, "Seed candidates", "cpu backend: linear peptides per spectrum used as seeds to find partners. \
                Default 30; also the number of pairs rescored exactly by reflect/gpu.",
                |ui| num(ui, &mut xl.candidates, 1..=1000, 1.0));
            row(ui, "Chain survivors", "reflect/gpu backends: peptides kept per spectrum after the per-chain filter. \
                Default 20000; raise it if you suspect missed links, lower it for speed.",
                |ui| num(ui, &mut xl.chain_survivors, 100..=1_000_000, 100.0));
        });
    });
}

/// Residues and protein termini one end of a heterobifunctional linker reacts with
fn site_row(ui: &mut Ui, label: &str, end: &str, site: &mut settings::Site) {
    let tip = format!(
        "What the {end} end of the linker reacts with: one-letter residue codes plus, optionally, the \
         protein termini.\nEDC / DMTMM: Site 1 = K + protein N-terminus (amines), \
         Site 2 = DE + protein C-terminus (carboxyls)."
    );
    row(ui, label, &tip, |ui| {
        let mut r = ui.add(TextEdit::singleline(&mut site.residues).hint_text("e.g. K").desired_width(80.0));
        r |= ui.checkbox(&mut site.protein_nterm, "protein N-terminus");
        r | ui.checkbox(&mut site.protein_cterm, "protein C-terminus")
    });
}

fn cleavable_selected(xl: &settings::Crosslink) -> bool {
    xl.reagent.is_cleavable() || (xl.reagent == Reagent::Custom && xl.custom_cleavable)
}
