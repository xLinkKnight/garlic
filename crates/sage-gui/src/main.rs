//! `sage-gui`: a desktop window for editing Sage search settings, choosing
//! input files and running `sage.exe`. Every field explains itself on hover.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod progress;
mod runner;
mod settings;
mod system;
mod usage;
mod theme;
mod unimod;
mod viewer;

use eframe::egui::{self, Color32, ComboBox, Response, RichText, TextEdit, Ui};
use runner::{Job, Status, Step};
use settings::{Backend, ModKind, Modification, Reagent, Settings, Silac, TolUnit, Tolerance, Tmt};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};

/// Start-up frames used to measure the tabs: two per tab
const FIT_STEPS: usize = 2 * Tab::ALL.len();

/// The project on GitHub (the link next to the light/dark switch)
const REPO_URL: &str = "https://github.com/xLinkKnight/garlic";

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Sage search")
            // the smiling garlic, drawn by examples/make_icon.rs
            .with_icon(
                eframe::icon_data::from_png_bytes(include_bytes!("../assets/garlic.png"))
                    .expect("the bundled icon is a valid PNG"),
            )
            .with_inner_size([1040.0, 780.0])
            .with_min_inner_size([720.0, 480.0]),
        ..Default::default()
    };
    // `sage-gui --view <results folder>` opens the annotated spectra of a search
    let mut args = std::env::args().skip(1);
    let view = match args.next().as_deref() {
        Some("--view") => args.next(),
        _ => None,
    };
    eframe::run_native(
        "Sage search",
        options,
        Box::new(move |cc| {
            let mut app = SageGui::new(&cc.egui_ctx);
            if let Some(folder) = view {
                app.view_spectra(Path::new(&folder));
            }
            Ok(Box::new(app))
        }),
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
    /// Shown once a search has been started
    Progress,
    /// Annotated spectra of the results
    Viewer,
}

impl Tab {
    const ALL: [(Tab, &'static str); 8] = [
        (Tab::Files, "Files"),
        (Tab::Database, "Database & enzyme"),
        (Tab::Modifications, "Modifications"),
        (Tab::Spectra, "Tolerances & spectra"),
        (Tab::Scoring, "Scoring & output"),
        (Tab::Crosslink, "Crosslinking"),
        (Tab::Progress, "Progress"),
        (Tab::Viewer, "Annotated spectra"),
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
    unimod: unimod::Lookup,
    /// Steps, metrics and errors of the current or last search
    tracker: Option<progress::Tracker>,
    /// The mobiligrams of the last search, being written in the background
    figures: Option<std::sync::mpsc::Receiver<Result<Vec<PathBuf>, String>>>,
    /// CPU / RAM / GPU / vRAM utilisation, shown on the Progress tab
    usage: usage::Monitor,
    dark: bool,
    /// Start-up frames spent measuring the tabs, see `measure_step`
    fit_step: usize,
    /// Largest height (points) by which a tab overflowed the scroll area
    overflow: f32,
    /// CPU threads sage.exe may use (`--threads`)
    threads: usize,
    /// RAM cap for sage.exe in GB, 0 = no cap
    max_ram_gb: f64,
    gpus: Vec<system::Gpu>,
    gpu_support: Option<system::GpuSupport>,
    gpu_probe: Option<std::sync::mpsc::Receiver<system::GpuSupport>>,
    /// sage.exe path the GPU support was checked for
    probed_sage: String,
    /// The "Annotated spectra" tab
    viewer: viewer::Viewer,
}

/// Where the GUI remembers its last settings and tool paths between runs
fn state_path() -> Option<PathBuf> {
    let base = std::env::var_os("APPDATA").or_else(|| std::env::var_os("HOME"))?;
    Some(PathBuf::from(base).join("sage-gui").join("state.json"))
}

impl SageGui {
    fn new(ctx: &egui::Context) -> Self {
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
            unimod: Default::default(),
            tracker: None,
            figures: None,
            usage: usage::Monitor::start(system::gpus()),
            dark: false,
            fit_step: 0,
            overflow: 0.0,
            threads: system::cpu_threads(),
            max_ram_gb: 0.0,
            gpus: system::gpus(),
            gpu_support: None,
            gpu_probe: None,
            probed_sage: String::new(),
            viewer: Default::default(),
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
            app.dark = state.get("dark").and_then(|v| v.as_bool()).unwrap_or(false);
            if let Some(n) = state.get("threads").and_then(|v| v.as_u64()) {
                app.threads = (n as usize).clamp(1, system::cpu_threads());
            }
            app.max_ram_gb = state.get("max_ram_gb").and_then(|v| v.as_f64()).unwrap_or(0.0).max(0.0);
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
        theme::install(ctx, app.dark);
        app
    }

    fn save_state(&self) {
        let Some(path) = state_path() else { return };
        let state = serde_json::json!({
            "sage_exe": self.sage_exe,
            "converter_exe": self.converter_exe,
            "parquet": self.parquet,
            "telemetry": self.telemetry,
            "dark": self.dark,
            "threads": self.threads,
            "max_ram_gb": self.max_ram_gb,
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

    /// Show a results folder on the "Annotated spectra" tab
    fn view_spectra(&mut self, out: &Path) {
        self.viewer.open(&out.display().to_string());
        self.tab = Tab::Viewer;
    }

    /// Check the settings, write the configuration into the output folder
    /// and start converting / searching
    fn start_run(&mut self, ctx: &egui::Context) {
        let (mut errors, mut warnings) = self.s.check();
        if !Path::new(&self.sage_exe).is_file() {
            errors.push("sage.exe not found: set its location under Files > Programs.".into());
        }
        // Only non-cleavable searches use the backend (see Crosslinking tab)
        let xl = &self.s.crosslink;
        if xl.enabled && !xl.is_cleavable() && xl.backend == Backend::Gpu {
            match &self.gpu_support {
                Some(system::GpuSupport::Yes(_)) => {}
                Some(system::GpuSupport::No) => errors.push(
                    "Crosslink backend gpu: this sage.exe was built without GPU support (Files > Resources). \
                     Choose cpu or auto, or rebuild sage.exe (GPU support is included by default: cargo build --release -p sage-cli)."
                        .into(),
                ),
                _ => warnings.push(
                    "Crosslink backend gpu: GPU support of sage.exe is unknown (Files > Resources); \
                     the search fails if it has none."
                        .into(),
                ),
            }
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

        let mut args = vec![config_path.into_os_string(), "--threads".into(), self.threads.to_string().into()];
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
            memory_limit: (self.max_ram_gb > 0.0).then(|| (self.max_ram_gb * system::GB) as u64),
        });
        self.save_state();
        let xl = &self.s.crosslink;
        self.tracker = Some(progress::Tracker::new(progress::Plan {
            files: spectra.len(),
            raw_files: self.s.spectra.iter().filter(|f| Settings::is_thermo_raw(f)).count(),
            crosslink: xl.enabled,
            linear: xl.enabled && xl.linear_search,
            quant: self.s.lfq || self.s.tmt != Tmt::None || self.s.silac != Silac::Off,
            gpu_backend: xl.enabled && matches!(xl.backend, Backend::Gpu | Backend::Auto),
            output_dir: out.clone(),
        }));
        self.job = Some(Job::start(steps, ctx.clone()));
        self.figures = None;
        self.tab = Tab::Progress;
    }

    /// Feed new log lines to the tracker, and close it when the job ends
    fn follow_job(&mut self, ctx: &egui::Context) {
        // the mobiligrams written in the background: listed with the search's outputs
        if let (Some(rx), Some(tracker)) = (&self.figures, &mut self.tracker) {
            match rx.try_recv() {
                Ok(Ok(written)) => {
                    for path in written {
                        let size = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
                        tracker.outputs.push((path, size));
                    }
                    self.figures = None;
                }
                Ok(Err(e)) => {
                    tracker.warnings.push(format!("mobiligram not written: {e}"));
                    self.figures = None;
                }
                Err(std::sync::mpsc::TryRecvError::Empty) => {}
                Err(std::sync::mpsc::TryRecvError::Disconnected) => self.figures = None,
            }
        }
        let (Some(job), Some(tracker)) = (&self.job, &mut self.tracker) else { return };
        let log = job.log.lock().unwrap().clone();
        tracker.update(&log, job.total.load(std::sync::atomic::Ordering::Relaxed));
        match job.status() {
            Status::Running(_) => {}
            Status::Finished if tracker.state != progress::State::Finished => {
                tracker.finish(progress::State::Finished, None, &log);
                // timsTOF data: the mobiligrams next to the results, written in
                // the background so the window stays responsive
                let out = tracker.plan.output_dir.clone();
                let (tx, rx) = std::sync::mpsc::channel();
                let ctx = ctx.clone();
                let _ = std::thread::Builder::new().name("mobiligrams".into()).spawn(move || {
                    let _ = tx.send(viewer::mobiligram::write_figures(&out));
                    ctx.request_repaint();
                });
                self.figures = Some(rx);
            }
            Status::Finished => {}
            Status::Cancelled => tracker.finish(progress::State::Cancelled, None, &log),
            Status::Failed(e) => tracker.finish(progress::State::Failed, Some(&e), &log),
        }
    }

    /// Start-up only: record how far `tab` overflows the scroll area, and
    /// once every tab has been measured, grow the window so the tallest tab
    /// needs no scrolling (limited by the screen height). Each tab is drawn
    /// invisibly for two frames, since grids size their rows from the previous frame.
    fn measure_step(&mut self, ctx: &egui::Context, overflow: f32) {
        if self.fit_step % 2 == 1 {
            self.overflow = self.overflow.max(overflow);
        }
        self.fit_step += 1;
        ctx.request_repaint();
        if self.fit_step < FIT_STEPS || self.overflow <= 0.0 {
            return;
        }
        let window = ctx.screen_rect().size();
        let mut height = window.y + self.overflow + 8.0;
        let monitor = ctx.input(|i| i.viewport().monitor_size);
        if let Some(monitor) = monitor {
            // egui-winit reports the monitor in physical pixels
            let monitor = monitor / ctx.pixels_per_point();
            // leave room for the title bar and the taskbar
            height = height.min(monitor.y - 90.0);
            ctx.send_viewport_cmd(egui::ViewportCommand::OuterPosition(egui::pos2(
                ((monitor.x - window.x) / 2.0).max(0.0),
                ((monitor.y - height) / 2.0 - 30.0).max(0.0),
            )));
        }
        ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(egui::vec2(window.x, height)));
    }

    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw.dropped_files.iter().filter_map(|f| f.path.clone()).collect()
        });
        for path in dropped {
            let lower = path.to_string_lossy().to_ascii_lowercase();
            let text = path.display().to_string();
            if [".fasta", ".fa", ".faa", ".fas"].iter().any(|e| lower.ends_with(e)) {
                // a two-database search: fill Database 1, then Database 2
                if !self.s.two_databases() {
                    self.s.fasta = text;
                } else if self.s.database_1.trim().is_empty() {
                    self.s.database_1 = text;
                } else {
                    self.s.database_2 = text;
                }
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
        self.follow_job(ctx);
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
                        std::mem::take(&mut self.s.database_1),
                        std::mem::take(&mut self.s.database_2),
                    );
                    self.s = Settings::default();
                    (self.s.fasta, self.s.spectra, self.s.output_directory, self.s.database_1, self.s.database_2) =
                        files;
                }
                ui.toggle_value(&mut self.show_json, "Show JSON")
                    .on_hover_text("Show the configuration file that will be passed to sage.exe.");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if theme::theme_switch(ui, &mut self.dark) {
                        theme::set_dark(ctx, self.dark);
                    }
                    // right to left: left of the theme switch
                    ui.add_space(12.0);
                    theme::github_link(ui, REPO_URL);
                });
            });
            ui.add_space(2.0);
            ui.horizontal_wrapped(|ui| {
                for (tab, name) in Tab::ALL {
                    if tab == Tab::Progress {
                        let Some(state) = self.tracker.as_ref().map(|t| t.state) else { continue };
                        let selected = self.tab == tab;
                        let (fill, border) = progress_tab_colors(ui, state, selected);
                        let r = ui.add(
                            egui::Button::new(theme::bold(name)).selected(selected).fill(fill).stroke(egui::Stroke::NONE),
                        );
                        if r.clicked() {
                            self.tab = tab;
                        }
                        ui.painter().rect_stroke(
                            r.rect.expand(2.0),
                            4.0,
                            egui::Stroke::new(2.0_f32, border),
                            egui::StrokeKind::Outside,
                        );
                    } else {
                        ui.selectable_value(&mut self.tab, tab, theme::bold(name));
                    }
                }
            });
            ui.add_space(2.0);
        });

        egui::TopBottomPanel::bottom("actions").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                let running = self.running();
                let run = ui
                    .add_enabled(!running, theme::tone_button(ui, "▶ Run search", theme::Tone::Blue))
                    .on_hover_text(
                        "Check the settings, convert any Thermo .raw files to mzML, write \
                         sage_search.json into the output folder and run sage.exe.",
                    );
                if run.clicked() {
                    self.start_run(ctx);
                }
                if ui.add_enabled(running, theme::tone_button(ui, "■ Cancel", theme::Tone::Tan)).clicked() {
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
                if ui
                    .add_enabled(
                        !running && Path::new(&out).join("results.json").is_file(),
                        egui::Button::new("View spectra"),
                    )
                    .on_hover_text(VIEW_TIP)
                    .clicked()
                {
                    self.view_spectra(Path::new(&out));
                }
                if let Some(job) = &self.job {
                    let (text, color) = match job.status() {
                        Status::Running(step) => (format!("Running: {step}"), Color32::from_rgb(70, 130, 220)),
                        Status::Finished => ("Finished".to_string(), Color32::from_rgb(40, 160, 70)),
                        Status::Failed(_) => (
                            "Failed: see the Progress tab for what went wrong".to_string(),
                            Color32::from_rgb(210, 60, 60),
                        ),
                        Status::Cancelled => ("Cancelled".to_string(), Color32::GRAY),
                    };
                    if job.is_running() {
                        ui.spinner();
                    }
                    ui.label(RichText::new(text).color(color));
                }
                if self.tab == Tab::Progress {
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.add_space(12.0);
                        self.usage_meters(ui);
                    });
                }
            });
            ui.add_space(4.0);
        });

        // Settings problems and load/save messages; the search log is on the Progress tab
        if !self.problems.is_empty() {
            egui::TopBottomPanel::bottom("problems").show(ctx, |ui| {
                ui.add_space(4.0);
                egui::ScrollArea::vertical().max_height(160.0).show(ui, |ui| {
                    for (is_error, msg) in &self.problems {
                        let (icon, color) = if *is_error {
                            ("✖", Color32::from_rgb(210, 60, 60))
                        } else {
                            ("ℹ", Color32::from_rgb(200, 150, 30))
                        };
                        ui.label(RichText::new(format!("{icon} {msg}")).color(color));
                    }
                });
                ui.add_space(2.0);
            });
        }

        if self.tab == Tab::Progress {
            self.progress_log(ctx);
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
            let measuring = self.fit_step < FIT_STEPS;
            let tab = if measuring { Tab::ALL[self.fit_step / 2].0 } else { self.tab };
            let out = egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                if measuring {
                    ui.set_invisible();
                }
                ui.add_space(6.0);
                let running = self.running();
                ui.add_enabled_ui(!running, |ui| match tab {
                    Tab::Files => self.files_tab(ui),
                    Tab::Database => database_tab(ui, &mut self.s),
                    Tab::Modifications => modifications_tab(ui, &mut self.s, &mut self.unimod),
                    Tab::Spectra => spectra_tab(ui, &mut self.s),
                    Tab::Scoring => scoring_tab(ui, &mut self.s),
                    Tab::Crosslink => crosslink_tab(ui, &mut self.s),
                    Tab::Progress => self.progress_tab(ui),
                    Tab::Viewer => self.viewer.ui(ui, &self.s.output_directory),
                });
            });
            if measuring {
                self.measure_step(ctx, out.content_size.y - out.inner_rect.height());
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Widgets

/// One labelled row in a two-column grid; `tip` shows on hover over both
fn row(ui: &mut Ui, label: &str, tip: &str, add: impl FnOnce(&mut Ui) -> Response) {
    let muted = theme::muted(ui);
    ui.label(RichText::new(label).color(muted)).on_hover_text(tip);
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
    theme::section(ui, text);
}

fn note(ui: &mut Ui, text: &str) {
    ui.label(RichText::new(text).weak());
    ui.add_space(4.0);
}

/// A number typed in, or stepped by `step` with the + / − buttons
fn num<T: egui::emath::Numeric>(ui: &mut Ui, v: &mut T, range: RangeInclusive<T>, step: f64) -> Response {
    theme::number(ui, v, range, step, |d| d)
}

fn mass(ui: &mut Ui, v: &mut f64) -> Response {
    theme::number(ui, v, -10_000.0..=10_000.0, 1.0, |d| d.max_decimals(6).suffix(" Da"))
}

/// A value Sage can choose itself: unchecked writes nothing to the JSON
fn optional<T: egui::emath::Numeric>(
    ui: &mut Ui,
    v: &mut Option<T>,
    fallback: T,
    auto: &str,
    range: RangeInclusive<T>,
    step: f64,
) -> Response {
    let mut on = v.is_some();
    let mut r = ui.checkbox(&mut on, "set");
    if on != v.is_some() {
        *v = on.then_some(fallback);
    }
    match v {
        Some(x) => r |= theme::number(ui, x, range, step, |d| d),
        None => r |= ui.label(RichText::new(auto).weak()),
    }
    r
}

const CHARGE_TIP: &str = "Precursor charge states to try when the file does not report one (or always, \
    with 'Override charge'). Default 2 to 4; crosslinked peptides are usually 3+ to 6+ (sometimes up to 8+).";

/// Precursor charge range; each bound is limited by the other so min <= max
fn charge_range(ui: &mut Ui, z: &mut (u8, u8)) -> Response {
    let hi = z.1;
    let r = theme::number(ui, &mut z.0, 1..=hi, 1.0, |d| d.prefix("from ").suffix("+"));
    let lo = z.0;
    r | theme::number(ui, &mut z.1, lo..=10, 1.0, |d| d.prefix("to ").suffix("+"))
}

fn tolerance(ui: &mut Ui, id: &str, t: &mut Tolerance) -> Response {
    let step = if t.unit == TolUnit::Ppm { 1.0 } else { 0.01 };
    let any = -10_000.0..=10_000.0;
    let mut r = theme::number(ui, &mut t.lo, any.clone(), step, |d| d.max_decimals(4).prefix("from "));
    r |= theme::number(ui, &mut t.hi, any, step, |d| d.max_decimals(4).prefix("to "));
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

/// Files tab: two FASTA files for a crosslink search between them, e.g. two
/// complexes. Crosslinks only between a Database 1 and a Database 2 protein;
/// linear peptides in both, as one database.
fn crosslink_databases(ui: &mut Ui, s: &mut Settings) {
    let accent = theme::accent(ui);
    egui::Frame::group(ui.style())
        .stroke(egui::Stroke::new(1.0_f32, accent.gamma_multiply(0.6)))
        .inner_margin(egui::Margin::same(10))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(theme::bold("Two-database crosslink search"));
                ui.label(
                    RichText::new(" FOR COMPLEX CROSSLINK SEARCHES ")
                        .small()
                        .strong()
                        .color(Color32::WHITE)
                        .background_color(accent),
                );
            });
            note(
                ui,
                "Crosslinks are searched only between a protein of Database 1 and a protein of Database 2 \
                 (e.g. two complexes, or a bait and its candidate partners). Crosslinks within either database, \
                 and loop-links, are not searched. Linear peptides are searched in both, as one database. Setting \
                 either file replaces the FASTA file above.",
            );
            let pick = || {
                rfd::FileDialog::new()
                    .add_filter("FASTA", &["fasta", "fa", "faa", "fas"])
                    .add_filter("All files", &["*"])
                    .pick_file()
            };
            grid(ui, "crosslink_databases", |ui| {
                row(
                    ui,
                    "Database 1",
                    "Proteins on one side of the crosslinks (e.g. complex 1). A protein in both files links to either side.",
                    |ui| path_field(ui, &mut s.database_1, "C:\\db\\complex1.fasta", pick),
                );
                row(
                    ui,
                    "Database 2",
                    "Proteins on the other side of the crosslinks (e.g. complex 2).",
                    |ui| path_field(ui, &mut s.database_2, "C:\\db\\complex2.fasta", pick),
                );
            });
            if s.two_databases() {
                ui.horizontal(|ui| {
                    if !s.crosslink.enabled {
                        ui.label(
                            RichText::new("Crosslinking is off (Crosslinking tab): the two files are searched as one FASTA.")
                                .color(Color32::from_rgb(220, 140, 30)),
                        );
                    }
                    if ui.small_button("clear both").on_hover_text("Search the FASTA file above instead.").clicked() {
                        s.database_1.clear();
                        s.database_2.clear();
                    }
                });
            }
        });
}

const SPECTRA_TIP: &str = "Mass spectrometry runs to search.\n\
    • mzML / mzML.gz and MGF / MGF.gz: read directly\n\
    • Bruker timsTOF .d: add the .d FOLDER (\"Add Bruker .d folder…\")\n\
    • Thermo .raw: converted to mzML with msconvert or ThermoRawFileParser before the \
    search, into <output>\\converted_mzML (already converted files are reused)\n\
    You can also drag files onto this window.";

impl SageGui {
    fn files_tab(&mut self, ui: &mut Ui) {
        heading(ui, "Protein database");
        // Database 1 / 2 set: the FASTA file is not searched
        let two = self.s.two_databases();
        ui.add_enabled_ui(!two, |ui| {
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
        });
        if two {
            note(ui, "Not searched: Database 1 and Database 2 below replace the FASTA file. Clear both to use it.");
        }
        ui.add_space(6.0);
        crosslink_databases(ui, &mut self.s);

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
                "Folder where Sage writes results.sage.tsv (PSMs), peptides.sage.tsv and \
                 proteins.sage.tsv (peptides at 1% FDR and their protein groups), results.json (the full settings \
                 used), quantification tables and, for crosslink searches, crosslinks.sage.tsv \
                 plus crosslinks_xiView.csv and crosslinks.mzid (for xiVIEW, PRIDE) and \
                 crosslinks_proXL.xml (for ProXL).\n\
                 Created if it does not exist. Example: C:\\results\\2026-11-01_hela",
                |ui| {
                    // dialogs open where the current output folder is, or its closest existing parent
                    let current = self.s.output_directory.clone();
                    let start = Path::new(&current)
                        .ancestors()
                        .find(|d| d.is_dir() && !d.as_os_str().is_empty())
                        .map(Path::to_path_buf);
                    let start2 = start.clone();
                    let mut r = path_field(ui, &mut self.s.output_directory, "C:\\results\\my_search", || {
                        let mut dialog = rfd::FileDialog::new().set_title("Choose the output folder");
                        if let Some(dir) = start {
                            dialog = dialog.set_directory(dir);
                        }
                        dialog.pick_folder()
                    });
                    // A Save dialog: go to the location, type the new folder's name, press Save
                    let new = ui.button("New folder…").on_hover_text(
                        "Create a new folder for the results: go to where it should be, type its name \
                         in 'File name' and press Save. The folder is created and selected.",
                    );
                    if new.clicked() {
                        let mut dialog = rfd::FileDialog::new()
                            .set_title("New output folder: go to its location, type the folder name, press Save")
                            .set_file_name("sage_results");
                        if let Some(dir) = start2 {
                            dialog = dialog.set_directory(dir);
                        }
                        if let Some(path) = dialog.save_file() {
                            match std::fs::create_dir_all(&path) {
                                Ok(()) => {
                                    self.s.output_directory = path.display().to_string();
                                    self.problems = vec![(false, format!("Created {}", path.display()))];
                                }
                                Err(e) => {
                                    self.problems = vec![(true, format!("Could not create {}: {e}", path.display()))]
                                }
                            }
                        }
                    }
                    r |= new;
                    r
                },
            );
            let out = self.s.output_directory.trim();
            if !out.is_empty() && !Path::new(out).is_dir() {
                ui.label("");
                ui.label(
                    RichText::new("This folder does not exist yet: it will be created when the search starts.")
                        .color(theme::muted(ui)),
                );
                ui.end_row();
            }
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

        self.resources(ui);
    }

    /// The Progress tab: steps, progress bar, log, and the report or error
    fn progress_tab(&mut self, ui: &mut Ui) {
        let Some(t) = &self.tracker else {
            heading(ui, "Search progress");
            note(ui, "No search has been run yet. Press ▶ Run search to start one; its progress appears here.");
            return;
        };
        let (ok, bad, warn) = (Color32::from_rgb(40, 160, 70), Color32::from_rgb(210, 60, 60), Color32::from_rgb(220, 135, 30));
        let muted = theme::muted(ui);
        // Results folder to open in Annotator, once the tracker is no longer borrowed
        let mut annotate: Option<PathBuf> = None;
        let running = t.state == progress::State::Running;
        // keep clear of the scroll bar on the right
        ui.set_max_width(ui.available_width() - 20.0);

        heading(ui, "Search progress");
        ui.horizontal(|ui| {
            let (text, color) = match t.state {
                progress::State::Running => (format!("Running: {}", t.activity()), theme::accent(ui)),
                progress::State::Finished => ("✔ Search finished successfully".to_string(), ok),
                progress::State::Failed => (
                    format!("✖ Search failed: {}", t.diagnosis.as_ref().map_or("see below", |d| d.title.as_str())),
                    bad,
                ),
                progress::State::Cancelled => ("■ Search cancelled".to_string(), muted),
            };
            if running {
                ui.spinner();
            }
            ui.label(theme::bold(text).color(color));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(RichText::new(format!("elapsed {}", progress::format_duration(t.elapsed()))).color(muted));
            });
        });
        let fill = match t.state {
            progress::State::Running => theme::accent(ui),
            progress::State::Finished => ok,
            progress::State::Failed => bad,
            progress::State::Cancelled => muted,
        };
        let bar = ui.add(
            egui::ProgressBar::new(t.fraction())
                .fill(fill)
                .animate(running),
        );
        // Percentage centred on the bar, dark blue on a pale backdrop so it reads
        // over every fill colour and over the empty track
        let navy = Color32::from_rgb(15, 40, 100);
        let galley = ui.painter().layout_no_wrap(
            format!("{:.0}%", 100.0 * t.fraction()),
            egui::FontId::proportional(14.0),
            navy,
        );
        let pos = bar.rect.center() - galley.size() / 2.0;
        let backdrop = egui::Rect::from_min_size(pos, galley.size()).expand2(egui::vec2(6.0, 0.0));
        ui.painter().rect_filled(backdrop, 4.0, Color32::from_white_alpha(190));
        ui.painter().galley(pos, galley, navy);
        if running {
            // keep the elapsed time and animation moving
            ui.ctx().request_repaint_after(std::time::Duration::from_millis(250));
        }

        heading(ui, "Steps");
        egui::Grid::new("progress_steps").num_columns(3).spacing([12.0, 6.0]).show(ui, |ui| {
            for (i, stage) in t.stages.iter().enumerate() {
                let done = t.took[i].is_some() && (i < t.current || t.state == progress::State::Finished);
                let current = i == t.current && !done;
                if done {
                    ui.label(RichText::new("✔").color(ok));
                } else if current && running {
                    ui.spinner();
                } else if current && t.state == progress::State::Failed {
                    ui.label(RichText::new("✖").color(bad));
                } else if current && t.state == progress::State::Cancelled {
                    ui.label(RichText::new("■").color(muted));
                } else {
                    ui.label(RichText::new("○").color(muted));
                }
                let label = RichText::new(stage.label());
                ui.label(if done || current { label } else { label.color(muted) });
                let detail = match (stage, done) {
                    (progress::Stage::Search, _) if t.plan.files > 0 && t.stage() != progress::Stage::Linear => {
                        format!("{} of {} files", t.files_done, t.plan.files)
                    }
                    (progress::Stage::Linear, true) => format!("{0} of {0} files", t.plan.files),
                    (progress::Stage::Linear, false) if t.plan.files > 0 => {
                        format!("{} of {} files", t.files_done, t.plan.files)
                    }
                    (progress::Stage::Convert, _) => format!("{} of {} files", t.converted, t.plan.raw_files),
                    _ => String::new(),
                };
                let time = t.took[i].filter(|_| done).map(progress::format_duration).unwrap_or_default();
                let sep = if !detail.is_empty() && !time.is_empty() { "  ·  " } else { "" };
                ui.label(RichText::new(format!("{detail}{sep}{time}")).color(muted));
                ui.end_row();
            }
        });

        if t.state == progress::State::Failed {
            if let Some(d) = &t.diagnosis {
                heading(ui, "What went wrong");
                egui::Frame::group(ui.style())
                    .stroke(egui::Stroke::new(1.5_f32, bad))
                    .inner_margin(10.0)
                    .show(ui, |ui| {
                        ui.set_width(ui.available_width());
                        ui.label(theme::bold(&d.title).color(bad).size(16.0));
                        ui.label(&d.cause);
                        ui.add_space(6.0);
                        ui.label(theme::bold("How to fix it"));
                        for fix in &d.fixes {
                            ui.label(format!("•  {fix}"));
                        }
                        if !d.evidence.is_empty() {
                            ui.add_space(6.0);
                            ui.label(theme::bold("Error message"));
                            for line in &d.evidence {
                                ui.label(RichText::new(line).monospace().size(12.0).color(muted));
                            }
                        }
                    });
            }
        }

        if t.state == progress::State::Finished {
            heading(ui, "Report");
            if t.identifications() == Some(0) {
                let hint = if t.plan.crosslink {
                    "No crosslinks were identified at 1% FDR. Check the reagent and its reactive residues, the \
                     precursor charge range (crosslinks are usually 3+ to 6+), the tolerances, and that the \
                     FASTA contains the crosslinked proteins."
                } else {
                    "No peptides were identified at 1% FDR. Check that the FASTA matches the sample's species, \
                     the enzyme, the precursor and fragment tolerances, the fixed modifications (e.g. \
                     carbamidomethyl C), and that the files contain MS2 spectra."
                };
                ui.label(RichText::new(format!("⚠ {hint}")).color(warn));
                ui.add_space(4.0);
            }
            for note in t.missing_exports() {
                ui.label(RichText::new(format!("ℹ {note}")).color(warn));
            }
            // Two columns of label / value pairs: the first half of the list on
            // the left, the second half on the right
            let half = t.metrics.len().div_ceil(2);
            egui::Grid::new("progress_metrics").num_columns(4).spacing([18.0, 6.0]).show(ui, |ui| {
                for i in 0..half {
                    for (label, value) in [t.metrics.get(i), t.metrics.get(i + half)].into_iter().flatten() {
                        ui.label(RichText::new(label).color(muted));
                        ui.label(RichText::new(value).strong());
                    }
                    ui.end_row();
                }
            });
            if !t.outputs.is_empty() {
                ui.add_space(6.0);
                ui.label(theme::bold("Output files"));
                let half = t.outputs.len().div_ceil(2);
                egui::Grid::new("progress_outputs").num_columns(4).spacing([18.0, 4.0]).show(ui, |ui| {
                    for i in 0..half {
                        for (path, size) in [t.outputs.get(i), t.outputs.get(i + half)].into_iter().flatten() {
                            ui.label(path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default())
                                .on_hover_text(path.display().to_string());
                            ui.label(RichText::new(human_size(*size)).color(muted));
                        }
                        ui.end_row();
                    }
                });
            }
            if self.figures.is_some() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(
                        RichText::new(format!(
                            "writing the mobiligrams to {}\\ (precursor m/z against ion mobility)…",
                            viewer::mobiligram::FOLDER
                        ))
                        .color(muted),
                    );
                });
            }
            let out = t.plan.output_dir.clone();
            ui.horizontal(|ui| {
                if ui.button("Open results folder").clicked() {
                    let _ = std::process::Command::new("explorer").arg(&out).spawn();
                }
                if ui.button("View spectra").on_hover_text(VIEW_TIP).clicked() {
                    annotate = Some(PathBuf::from(&out));
                }
            });
        }

        if !t.warnings.is_empty() && t.state != progress::State::Running {
            ui.add_space(6.0);
            egui::CollapsingHeader::new(RichText::new(format!("Warnings from Sage ({})", t.warnings.len())).color(warn))
                .id_salt("progress_warnings")
                .show(ui, |ui| {
                    for w in &t.warnings {
                        ui.label(format!("•  {w}"));
                    }
                });
        }

        if let Some(out) = annotate {
            self.view_spectra(&out);
        }
    }

    /// CPU, RAM, GPU and vRAM utilisation of the whole computer, right to left
    fn usage_meters(&self, ui: &mut Ui) {
        let Some(u) = self.usage.latest(ui.ctx()) else {
            ui.label(RichText::new("measuring usage…").color(theme::muted(ui)));
            return;
        };
        let gb = |b: u64| b as f64 / system::GB;
        let meter = |ui: &mut Ui, name: &str, pct: Option<f32>, tip: String| {
            let value = match pct {
                Some(p) => RichText::new(format!("{p:.0}%")).strong().color(if p >= 90.0 {
                    Color32::from_rgb(210, 60, 60)
                } else {
                    ui.visuals().text_color()
                }),
                None => RichText::new("–").color(theme::muted(ui)),
            };
            // right to left: value first, then its name
            let r = ui.label(value).on_hover_text(&tip);
            ui.label(RichText::new(name).color(theme::muted(ui))).on_hover_text(&tip);
            ui.add_space(10.0);
            r
        };
        let all_gpus = u
            .gpus
            .iter()
            .map(|g| match g.vram_pct() {
                Some(v) if g.vram_total >= 512 << 20 => format!(
                    "{}: {:.0}% busy, {:.1} of {:.1} GB vRAM ({v:.0}%)",
                    g.name, g.util, gb(g.vram_used), gb(g.vram_total)
                ),
                _ => format!("{}: {:.0}% busy (integrated, uses system RAM)", g.name, g.util),
            })
            .collect::<Vec<_>>()
            .join("
");
        let main = u.main_gpu().filter(|g| g.vram_total >= 512 << 20).or(u.main_gpu());
        let gpu_tip = format!("Whole computer, as in Task Manager (busiest engine type).
{all_gpus}");
        let vram_pct = main.filter(|g| g.vram_total >= 512 << 20).and_then(|g| g.vram_pct());
        meter(ui, "vRAM", vram_pct, match main {
            Some(g) if vram_pct.is_some() => format!(
                "{}: {:.1} of {:.1} GB dedicated video memory in use (whole computer).",
                g.name, gb(g.vram_used), gb(g.vram_total)
            ),
            _ => "No graphics card with its own video memory.".into(),
        });
        meter(ui, "GPU", main.map(|g| g.util), gpu_tip);
        meter(ui, "RAM", Some(u.ram_pct()), format!(
            "{:.1} of {:.1} GB of memory in use (whole computer).",
            gb(u.ram_used), gb(u.ram_total)
        ));
        meter(ui, "CPU", Some(u.cpu), format!(
            "Processor use of the whole computer across {} logical CPUs.",
            system::cpu_threads()
        ));
    }

    /// The search log, in a panel along the bottom of the Progress tab whose
    /// top edge can be dragged to show more or fewer lines
    fn progress_log(&self, ctx: &egui::Context) {
        let Some(job) = &self.job else { return };
        let lines = job.log.lock().unwrap().clone();
        let (bad, warn) = (Color32::from_rgb(210, 60, 60), Color32::from_rgb(220, 135, 30));
        let max = ctx.screen_rect().height() * 0.8;
        egui::TopBottomPanel::bottom("progress_log")
            .resizable(true)
            .default_height(260.0)
            .height_range(90.0..=max)
            .show(ctx, |ui| {
                // grip on the draggable top edge
                let top = ui.max_rect().top() - 3.0;
                let x = ui.max_rect().center().x;
                ui.painter().line_segment(
                    [egui::pos2(x - 22.0, top), egui::pos2(x + 22.0, top)],
                    egui::Stroke::new(3.0_f32, theme::muted(ui)),
                );
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.label(RichText::new("Log").size(17.0).color(theme::accent(ui)));
                    ui.label(RichText::new(format!("{} lines", lines.len())).color(theme::muted(ui)));
                    if ui.small_button("Copy log").on_hover_text("Copy the whole log, e.g. to report a problem.").clicked() {
                        ui.ctx().copy_text(lines.join("\n"));
                    }
                    ui.label(RichText::new("drag the top edge to resize").color(theme::muted(ui)).size(12.0));
                });
                let row_height = ui.text_style_height(&egui::TextStyle::Monospace);
                egui::Frame::group(ui.style()).show(ui, |ui| {
                    egui::ScrollArea::both()
                        .id_salt("progress_log")
                        .auto_shrink([false, false])
                        .stick_to_bottom(true)
                        .show_rows(ui, row_height, lines.len(), |ui, rows| {
                            for line in &lines[rows] {
                                let color = if line.contains(" ERROR ") || line.starts_with("Error") || line.starts_with("=== Failed") {
                                    bad
                                } else if line.contains(" WARN ") {
                                    warn
                                } else if line.starts_with("---") || line.starts_with("===") {
                                    theme::accent(ui)
                                } else {
                                    ui.visuals().text_color()
                                };
                                // one row per line (show_rows needs equal heights): scroll, don't wrap
                                ui.add(egui::Label::new(RichText::new(line).monospace().size(12.0).color(color)).extend());
                            }
                        });
                });
            });
    }

    /// CPU threads, RAM cap and GPUs available to sage.exe
    fn resources(&mut self, ui: &mut Ui) {
        // Re-check GPU support whenever sage.exe changes
        if self.probed_sage != self.sage_exe {
            self.probed_sage = self.sage_exe.clone();
            self.gpu_support = None;
            self.gpu_probe = Some(system::probe_in_background(Path::new(&self.sage_exe), ui.ctx()));
        }
        if let Some(Ok(support)) = self.gpu_probe.as_ref().map(|rx| rx.try_recv()) {
            self.gpu_support = Some(support);
            self.gpu_probe = None;
        }

        heading(ui, "Resources");
        let cpus = system::cpu_threads();
        let total_gb = system::total_ram().map(|b| b as f64 / system::GB);
        grid(ui, "resources", |ui| {
            let tip = format!(
                "Maximum CPU threads sage.exe may use (passed as --threads). This computer has {cpus}. \
                 Fewer threads leave the computer responsive for other work; searches take longer."
            );
            row(ui, "CPU threads", &tip, |ui| {
                let r = theme::number(ui, &mut self.threads, 1..=cpus, 1.0, |d| d);
                ui.label(RichText::new(format!("of {cpus}")).color(theme::muted(ui)));
                let all = ui.add_enabled(self.threads != cpus, egui::Button::new("Use all").small());
                if all.clicked() {
                    self.threads = cpus;
                }
                r | all
            });
            let installed = total_gb.map(|g| format!("{g:.1} GB")).unwrap_or("unknown".into());
            let tip = format!(
                "Most memory sage.exe may use, in GB. 0 = no limit (all available). Installed: {installed}.\n\
                 Enforced by Windows: if the search needs more, sage.exe is stopped with a \"ran out of memory\" \
                 message instead of slowing the whole computer down. Large databases, many variable \
                 modifications, crosslink searches and many files at once need the most memory."
            );
            row(ui, "Max RAM", &tip, |ui| {
                let max = total_gb.unwrap_or(1024.0).ceil();
                let r = theme::number(ui, &mut self.max_ram_gb, 0.0..=max, 1.0, |d| d.max_decimals(1).suffix(" GB"));
                let text = match (self.max_ram_gb > 0.0, total_gb) {
                    (false, Some(t)) => format!("0 = all available ({t:.1} GB installed)"),
                    (false, None) => "0 = all available".into(),
                    (true, Some(t)) => format!("of {t:.1} GB installed"),
                    (true, None) => String::new(),
                };
                ui.label(RichText::new(text).color(theme::muted(ui)));
                r
            });
            row(ui, "GPU", "Graphics adapters Windows reports, with their dedicated video memory (vRAM). \
                Integrated graphics have no vRAM of their own and borrow system RAM.", |ui| {
                ui.vertical(|ui| {
                    if self.gpus.is_empty() {
                        ui.label(RichText::new("none found").color(theme::muted(ui)));
                    }
                    for gpu in &self.gpus {
                        let memory = if gpu.vram >= 512 << 20 {
                            format!("{:.1} GB vRAM", gpu.vram as f64 / system::GB)
                        } else {
                            format!("integrated, shares up to {:.1} GB of system RAM", gpu.shared as f64 / system::GB)
                        };
                        ui.horizontal(|ui| {
                            ui.label(&gpu.name);
                            ui.label(RichText::new(memory).color(theme::muted(ui)));
                        });
                    }
                })
                .response
            });
            row(ui, "GPU support", "Whether this sage.exe can use a GPU. Only the gpu / auto backend of \
                non-cleavable crosslink searches uses one (Crosslinking tab); everything else runs on the CPU.\n\
                It needs a sage.exe with GPU support (included by default: cargo build --release -p sage-cli), and a Vulkan or \
                DirectX 12 graphics driver.", |ui| {
                let (ok, warn) = (Color32::from_rgb(40, 160, 70), Color32::from_rgb(200, 130, 30));
                let r = match &self.gpu_support {
                    None => {
                        ui.spinner();
                        ui.label(RichText::new("checking sage.exe…").color(theme::muted(ui)))
                    }
                    Some(system::GpuSupport::Yes(adapters)) => {
                        // hardware adapters only, each name once (it may appear per backend)
                        let mut names: Vec<&str> = Vec::new();
                        for name in adapters
                            .iter()
                            .filter(|a| !a.ends_with("| Cpu"))
                            .map(|a| a.split(" | ").next().unwrap_or(a))
                        {
                            if !names.contains(&name) {
                                names.push(name);
                            }
                        }
                        if names.is_empty() {
                            ui.label(
                                RichText::new("✖ sage.exe has GPU support, but no Vulkan / DirectX 12 GPU is usable \
                                    (update the graphics driver)")
                                    .color(warn),
                            )
                        } else {
                            ui.label(RichText::new(format!("✔ ready: sage.exe can use {}", names.join(", "))).color(ok))
                        }
                    }
                    Some(system::GpuSupport::No) => ui.label(
                        RichText::new("✖ this sage.exe was built without GPU support (rebuild it: cargo build --release -p sage-cli)")
                            .color(warn),
                    ),
                    Some(system::GpuSupport::Unknown(why)) => ui.label(RichText::new(format!("? {why}")).color(warn)),
                };
                let again = ui.add_enabled(self.gpu_probe.is_none(), egui::Button::new("Check again").small());
                if again.clicked() {
                    self.probed_sage.clear();
                    self.gpus = system::gpus();
                }
                r | again
            });
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
                ui.label(theme::bold("Cleave at")).on_hover_text("Residues at which this enzyme cuts, e.g. DE for Glu-C in phosphate buffer.");
                ui.label(theme::bold("Not before")).on_hover_text("Do not cut when the next residue is one of these, e.g. P. Leave empty for none.");
                ui.label(theme::bold("Side")).on_hover_text("After the residue (C-terminal, like trypsin) or before it (N-terminal, like Asp-N).");
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
            Sage's own default when the field is absent is 0.", |ui| num(ui, &mut s.missed_cleavages, 0..=10, 1.0));
        row(ui, "Peptide length", "Shortest and longest peptide (amino acids) to search. Typical: 7 to 50.", |ui| {
            let r = theme::number(ui, &mut s.min_len, 1..=100, 1.0, |d| d.prefix("min "));
            r | theme::number(ui, &mut s.max_len, 1..=200, 1.0, |d| d.prefix("max "))
        });
        row(ui, "Semi-enzymatic", "Also allow peptides with only one enzymatic end. Finds in-source fragments \
            and degradation products, but the database grows ~10x and searches get slower.",
            |ui| ui.checkbox(&mut s.semi_enzymatic, "semi-specific digestion"));
    });

    heading(ui, "Precursor charge states");
    note(ui, "Used by linear and crosslink searches alike.");
    grid(ui, "charges", |ui| {
        row(ui, "Precursor charge", CHARGE_TIP, |ui| {
            let r = charge_range(ui, &mut s.precursor_charge);
            let presets = [((2, 4), "Linear: 2+ to 4+"), ((3, 6), "Crosslinks: 3+ to 6+")];
            presets.into_iter().fold(r, |r, (z, label)| {
                let b = ui
                    .add_enabled(s.precursor_charge != z, egui::Button::new(label).small())
                    .on_hover_text("Typical range for this kind of search.");
                if b.clicked() {
                    s.precursor_charge = z;
                }
                r | b
            })
        });
        row(ui, "Override charge", "Ignore the charge state in the file and always try the range above.",
            |ui| ui.checkbox(&mut s.override_precursor_charge, "always use the charge range"));
    });

    heading(ui, "Peptides & fragments");
    grid(ui, "peptides", |ui| {
        row(ui, "Minimum peptide mass", "Lightest peptide (monoisotopic, Da) to search.\n\
            Unset = Sage's choice: 500 Da, or 350 Da for crosslink searches so small partner peptides are kept.",
            |ui| optional(ui, &mut s.peptide_min_mass, 500.0, "automatic (500 Da; 350 Da for crosslinks)", 0.0..=10000.0, 50.0));
        row(ui, "Maximum peptide mass", "Heaviest peptide (monoisotopic, Da) to search. Default 5000.",
            |ui| theme::number(ui, &mut s.peptide_max_mass, 100.0..=20000.0, 100.0, |d| d.suffix(" Da")));
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
            which are rarely observed and match noise easily.", |ui| num(ui, &mut s.min_ion_index, 0..=10, 1.0));
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
            |ui| num(ui, &mut s.bucket_size, 1024..=1_048_576, 1024.0));
        row(ui, "Prefilter database", "For very large databases (e.g. metaproteomics): first find peptides \
            that match at all, then search only those. Lowers memory use. Not allowed for crosslink searches.",
            |ui| ui.checkbox(&mut s.prefilter, "prefilter in chunks"));
        if s.prefilter {
            row(ui, "Prefilter chunk size", "Number of protein sequences per prefilter chunk. 0 = automatic.",
                |ui| num(ui, &mut s.prefilter_chunk_size, 0..=10_000_000, 10_000.0));
            row(ui, "Low-memory prefilter", "Use as little memory as possible during prefiltering, at some cost in speed.",
                |ui| ui.checkbox(&mut s.prefilter_low_memory, "low memory"));
        }
    });
    labelled_proteins(ui, s);
}

thread_local! {
    /// The proteins of the FASTA file(s) last shown: (each path with its
    /// modification time, proteins)
    static FASTA_PROTEINS: std::cell::RefCell<Option<(Vec<(String, Option<std::time::SystemTime>)>, std::rc::Rc<Vec<(String, String)>>)>> =
        const { std::cell::RefCell::new(None) };
}

/// The proteins of the FASTA file(s), read again only when a file changes;
/// a protein in several files once
fn fasta_proteins(files: &[&str]) -> Option<std::rc::Rc<Vec<(String, String)>>> {
    let key: Vec<(String, Option<std::time::SystemTime>)> = files
        .iter()
        .map(|f| (f.to_string(), std::fs::metadata(f).and_then(|m| m.modified()).ok()))
        .collect();
    FASTA_PROTEINS.with(|cache| {
        let mut cache = cache.borrow_mut();
        if let Some((k, proteins)) = cache.as_ref() {
            if *k == key {
                return Some(proteins.clone());
            }
        }
        let mut all: Vec<(String, String)> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for f in files {
            for p in settings::fasta_proteins(f)? {
                if seen.insert(p.0.clone()) {
                    all.push(p);
                }
            }
        }
        let proteins = std::rc::Rc::new(all);
        *cache = Some((key, proteins.clone()));
        Some(proteins)
    })
}

/// Proteins of the FASTA expressed in 15N and/or 13C media: pick each one's
/// label from the list
fn labelled_proteins(ui: &mut Ui, s: &mut Settings) {
    use sage_core::labelled::Isotope;
    heading(ui, "Isotope-labelled proteins");
    note(ui, "Proteins expressed in 15N and/or 13C media, every nitrogen and/or carbon atom heavy. Their peptides \
        are searched labelled (in linear and crosslink searches); the labels do not count towards the variable \
        modification limit, and a labelled residue keeps its fixed or variable modifications and still crosslinks.");
    let label_text = |l: Option<Isotope>| match l {
        None => "unlabelled",
        Some(Isotope::N15) => "15N",
        Some(Isotope::C13) => "13C",
        Some(Isotope::C13N15) => "13C + 15N",
    };
    let files = s.fasta_files();
    // one FASTA, or both crosslink databases
    let Some(proteins) = (!files.is_empty()).then(|| fasta_proteins(&files)).flatten() else {
        note(ui, "Select a FASTA file (Files tab) to pick proteins from it.");
        return;
    };
    // one protein: its label, accession and name
    let protein_row = |ui: &mut Ui, s: &mut Settings, salt: &str, accession: &str, name: &str| {
        ui.horizontal(|ui| {
            let mut label = s.labelled_proteins.get(accession).copied();
            let before = label;
            ComboBox::from_id_salt((salt, accession))
                .selected_text(label_text(label))
                .width(110.0)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut label, None, label_text(None));
                    for i in Isotope::ALL {
                        ui.selectable_value(&mut label, Some(i), label_text(Some(i)));
                    }
                });
            if label != before {
                match label {
                    Some(l) => {
                        s.labelled_proteins.insert(accession.to_string(), l);
                    }
                    None => {
                        s.labelled_proteins.remove(accession);
                    }
                }
            }
            let text = RichText::new(accession).monospace();
            ui.label(if label.is_some() { text.strong() } else { text });
            ui.label(RichText::new(name).color(theme::muted(ui)));
        });
    };
    let row_height = ui.spacing().interact_size.y + 6.0;

    // the chosen proteins, always shown (a protein may be named by a part of
    // its accession, e.g. in a configuration file)
    let chosen: Vec<(String, String)> = s
        .labelled_proteins
        .keys()
        .map(|k| {
            let name = proteins
                .iter()
                .find(|(a, _)| a == k || a.split('|').skip(1).any(|p| p == k))
                .map(|(_, n)| n.clone())
                .unwrap_or_default();
            (k.clone(), name)
        })
        .collect();
    if !chosen.is_empty() {
        egui::ScrollArea::vertical()
            .id_salt("labelled_proteins_chosen")
            .max_height(160.0)
            .auto_shrink([false, true])
            .show_rows(ui, row_height, chosen.len(), |ui, rows| {
                for (accession, name) in &chosen[rows] {
                    protein_row(ui, s, "labelled_chosen", accession, name);
                }
            });
        ui.add_space(4.0);
    }
    ui.horizontal(|ui| {
        let n = s.labelled_proteins.len();
        ui.label(theme::bold(format!("{n} labelled")));
        ui.add_enabled(n > 0, egui::Checkbox::new(&mut s.include_unlabelled, "also search them unlabelled"))
            .on_hover_text("On: the labelled proteins are searched both labelled and unlabelled (e.g. a labelled \
                protein mixed with its unlabelled form). Off: labelled only. Peptides they share with other proteins \
                are searched unlabelled for those.");
        if n > 0 && ui.small_button("clear").clicked() {
            s.labelled_proteins.clear();
        }
    });

    // every protein of the FASTA, hidden until opened: a large FASTA would
    // otherwise fill the tab
    egui::CollapsingHeader::new(format!("Choose proteins ({} in the FASTA)", proteins.len()))
        .id_salt("labelled_proteins_list")
        .default_open(false)
        .show(ui, |ui| {
            let filter_id = ui.id().with("labelled_filter");
            let mut filter: String = ui.data_mut(|d| d.get_temp::<String>(filter_id).unwrap_or_default());
            ui.horizontal(|ui| {
                ui.label("Search");
                ui.add(TextEdit::singleline(&mut filter).desired_width(240.0).hint_text("accession, gene or protein name"));
            });
            ui.data_mut(|d| d.insert_temp(filter_id, filter.clone()));
            let needle = filter.to_ascii_lowercase();
            let shown: Vec<&(String, String)> = proteins
                .iter()
                .filter(|(a, n)| {
                    needle.is_empty() || a.to_ascii_lowercase().contains(&needle) || n.to_ascii_lowercase().contains(&needle)
                })
                .collect();
            if shown.is_empty() {
                note(ui, "No protein matches the search.");
                return;
            }
            egui::ScrollArea::vertical()
                .id_salt("labelled_proteins")
                .max_height(220.0)
                .auto_shrink([false, true])
                .show_rows(ui, row_height, shown.len(), |ui, rows| {
                    for (accession, name) in &shown[rows] {
                        protein_row(ui, s, "labelled", accession, name);
                    }
                });
        });
}

const TARGET_TIP: &str = "Where the modification can occur:\n\
    • one amino acid: C, M, S, K, …\n\
    • ^ peptide N-terminus (any residue), ^Q only when the N-terminal residue is Q\n\
    • $ peptide C-terminus, $K\n\
    • [ protein N-terminus, [M\n\
    • ] protein C-terminus";

/// Enter a UniMod accession, fetch the record and add it at one of its sites
fn unimod_panel(ui: &mut Ui, s: &mut Settings, lk: &mut unimod::Lookup) {
    let tip = "Type a UniMod accession number (the Accession # column of unimod.org's modification list), \
        e.g. 35 for oxidation or UNIMOD:21 for phospho, and press Enter. The name, monoisotopic mass and \
        allowed sites are fetched from unimod.org (downloaded once, then cached).";
    let busy = lk.busy();
    ui.horizontal(|ui| {
        ui.label("UniMod accession #").on_hover_text(tip);
        let edit = ui
            .add(TextEdit::singleline(&mut lk.input).hint_text("e.g. 35").desired_width(100.0))
            .on_hover_text(tip);
        let enter = edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        let look = ui.add_enabled(!busy, egui::Button::new("Look up")).on_hover_text(tip);
        if !busy && (look.clicked() || enter) {
            lk.start(ui.ctx());
        }
        if busy {
            ui.spinner();
            ui.label(RichText::new("Fetching from unimod.org…").weak());
        }
    });

    let entry = match &lk.result {
        None => return,
        Some(Err(e)) => {
            ui.label(RichText::new(format!("✖ {e}")).color(Color32::from_rgb(210, 60, 60)));
            return;
        }
        Some(Ok(entry)) => entry.clone(),
    };
    ui.group(|ui| {
        ui.horizontal(|ui| {
            ui.label(theme::bold(format!("{}  {}", entry.proforma(), entry.title)));
            ui.label(format!("({})", entry.full_name));
            ui.hyperlink_to(
                "view on unimod.org",
                format!("https://www.unimod.org/modifications_view.php?editid1={}", entry.accession),
            );
        });
        ui.label(format!("Monoisotopic mass {:+.6} Da   ·   composition {}", entry.mono_mass, entry.composition));
        ui.add_space(4.0);
        egui::Grid::new("unimod_sites").num_columns(5).spacing([14.0, 4.0]).striped(true).show(ui, |ui| {
            ui.label(theme::bold("Site"));
            ui.label(theme::bold("Position"));
            ui.label(theme::bold("Classification"));
            ui.label(theme::bold("Sage target")).on_hover_text(TARGET_TIP);
            ui.label("");
            ui.end_row();
            for spec in &entry.specificities {
                ui.label(&spec.site);
                ui.label(&spec.position);
                ui.label(&spec.classification);
                let Some(target) = &spec.target else {
                    ui.label(RichText::new("not supported").weak())
                        .on_hover_text("Sage can't express this site (e.g. a residue that is not a standard amino acid).");
                    ui.label("");
                    ui.end_row();
                    continue;
                };
                ui.monospace(target);
                ui.horizontal(|ui| {
                    for (kind, label) in [(ModKind::Fixed, "Add fixed"), (ModKind::Variable, "Add variable")] {
                        let added = s.mods.iter().any(|m| {
                            m.kind == kind && m.target == *target && (m.mass - entry.mono_mass).abs() < 1e-6
                        });
                        let b = ui
                            .add_enabled(!added, egui::Button::new(if added { "✔ added" } else { label }).small())
                            .on_hover_text(format!(
                                "Add {} ({:+.6} Da) on {target} as a {} modification.",
                                entry.title,
                                entry.mono_mass,
                                if kind == ModKind::Fixed { "fixed" } else { "variable" }
                            ));
                        if b.clicked() {
                            s.mods.push(Modification {
                                target: target.clone(),
                                mass: entry.mono_mass,
                                name: entry.proforma(),
                                kind,
                                max_count: 0,
                            });
                        }
                    }
                });
                ui.end_row();
            }
        });
    });
}

fn modifications_tab(ui: &mut Ui, s: &mut Settings, lk: &mut unimod::Lookup) {
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
    ui.add_space(4.0);
    unimod_panel(ui, s, lk);
    ui.add_space(6.0);

    let mut remove = None;
    egui::Grid::new("mods").num_columns(6).spacing([12.0, 6.0]).striped(true).show(ui, |ui| {
        ui.label(theme::bold("Type")).on_hover_text("Fixed: always applied. Variable: searched with and without.");
        ui.label(theme::bold("Target")).on_hover_text(TARGET_TIP);
        ui.label(theme::bold("Mass shift")).on_hover_text("Monoisotopic mass added to the residue, in Da. Negative values allowed (e.g. -17.026549 for pyro-Glu from Q).");
        ui.label(theme::bold("Max per peptide")).on_hover_text("Variable modifications only: how often this modification may occur on one peptide. 0 = only limited by 'Max variable mods'.");
        ui.label(theme::bold("Name (optional)")).on_hover_text("ProForma name written to the results, e.g. UNIMOD:35 for oxidation. Leave empty to report the mass.");
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
            ui.add_enabled_ui(m.kind == ModKind::Variable, |ui| theme::number(ui, &mut m.max_count, 0..=10, 1.0, |d| d))
                .inner
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
            |ui| num(ui, &mut s.max_variable_mods, 1..=10, 1.0));
        row(ui, "Max combinations", "Hard cap on modified forms generated per peptide (including the unmodified \
            one). Unset = no cap. Useful against combinatorial explosion, e.g. 64.",
            |ui| optional(ui, &mut s.max_combinations, 64, "no limit", 1..=100_000, 16.0));
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
                let r = theme::number(ui, &mut s.isotope_errors.0, -5..=5, 1.0, |d| d.prefix("from "));
                r | theme::number(ui, &mut s.isotope_errors.1, -5..=5, 1.0, |d| d.prefix("to "))
            });
    });

    heading(ui, "Fragment charge");
    note(ui, "The precursor charge range is on the Database & enzyme tab.");
    grid(ui, "fragment_charge", |ui| {
        row(ui, "Max fragment charge", "Highest fragment ion charge to consider. Unset = precursor charge - 1. \
            Example: 1 for low charge peptides, 2 for most HCD data.",
            |ui| optional(ui, &mut s.max_fragment_charge, 2, "precursor charge - 1", 1..=6, 1.0));
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
            |ui| optional(ui, &mut s.max_peaks, 150, "automatic (150; 300 for DSSO/DSBU)", 10..=5000, 10.0));
        row(ui, "Min matched peaks", "A PSM needs at least this many matched b/y (etc.) ions. Default 4.",
            |ui| num(ui, &mut s.min_matched_peaks, 1..=50, 1.0));
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
            e.g. 5 when rescoring with Percolator or mokapot.", |ui| num(ui, &mut s.report_psms, 1..=100, 1.0));
        row(ui, "Predict retention time", "Fit a retention time model and use the RT error as a feature for \
            rescoring (recommended). Required for LFQ.", |ui| ui.checkbox(&mut s.predict_rt, "predict RT"));
        row(ui, "Protein grouping", "Group proteins that share the same peptide evidence (parsimony).",
            |ui| ui.checkbox(&mut s.protein_grouping, "group proteins"));
        row(ui, "Grouping peptide FDR", "Peptide q-value cutoff used when building protein groups. Default 0.01 (1%).",
            |ui| theme::number(ui, &mut s.protein_grouping_peptide_fdr, 0.0..=1.0, 0.01, |d| d.max_decimals(3)));
    });

    heading(ui, "Quantification");
    if s.crosslink.enabled {
        note(ui, "Crosslink search: label-free intensities are written per crosslinked peptide pair to \
            crosslinks.lfq.tsv; TMT reporter intensities are added as columns of crosslinks.sage.tsv; \
            SILAC ratios per crosslinked species go to crosslinks.silac.tsv.");
    }
    grid(ui, "quant", |ui| {
        row(ui, "Label-free (LFQ)", "Extract MS1 peak areas across all runs. Needs MS1 spectra in the files.\n\
            • Linear search: per peptide, in lfq.tsv\n\
            • Crosslink search: per crosslinked peptide pair (link-site isomers can't be told apart in MS1), \
            in crosslinks.lfq.tsv. Runs are aligned on crosslinks identified in several runs.",
            |ui| ui.checkbox(&mut s.lfq, "label-free quantification"));
        if s.lfq || s.silac != Silac::Off {
            row(ui, "LFQ q-value cutoff", "Only peptides (or crosslinks, by crosslink q-value) at or below this \
                q-value are quantified, by LFQ and SILAC. Default 0.01 (1% FDR).",
                |ui| theme::number(ui, &mut s.lfq_q, 0.0..=1.0, 0.01, |d| d.max_decimals(3)));
            row(ui, "MS1 tolerance", "Mass tolerance for tracing precursor isotopes in MS1, in ppm. Default 5; \
                10 for older or lower-resolution instruments.",
                |ui| theme::number(ui, &mut s.lfq_ppm, 0.5..=50.0, 1.0, |d| d.suffix(" ppm")));
            row(ui, "Min spectral angle", "Minimum agreement (0 to 1) between the observed and theoretical isotope \
                envelope for a peak to be integrated. Default 0.7.",
                |ui| theme::number(ui, &mut s.lfq_spectral_angle, 0.0..=1.0, 0.05, |d| d.max_decimals(2)));
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
                |ui| num(ui, &mut s.tmt_level, 2..=3, 1.0));
            row(ui, "Signal-to-noise", "Report reporter signal-to-noise instead of intensity (Thermo files with noise data).",
                |ui| ui.checkbox(&mut s.tmt_sn, "use S/N"));
        }
        row(ui, "SILAC", "Stable isotope labels in cell culture. Every peptide is searched light and in each \
            labelled channel, with all of its K and R labelled (no mixed states), so the search takes 2 or 3 times \
            as long. A labelled lysine still crosslinks.\n\
            Each identified peptide (and crosslinked species) is traced in MS1 in every channel at once, over one \
            shared peak, with the MS1 tolerance and spectral angle above: silac.sage.tsv (intensities and ratios \
            to light per run), silac_proteins.sage.tsv (median ratio of each protein's unique peptides) and \
            crosslinks.silac.tsv. results.sage.tsv and crosslinks.sage.tsv get a silac_channel column.",
            |ui| {
                ComboBox::from_id_salt("silac")
                    .selected_text(s.silac.label())
                    .show_ui(ui, |ui| {
                        for x in Silac::ALL {
                            ui.selectable_value(&mut s.silac, x, x.label());
                        }
                    })
                    .response
            });
        if s.silac == Silac::Triplex {
            row(ui, "Medium labels", "Label masses of the medium channel on K and R (Da); 0: not labelled.",
                |ui| silac_labels(ui, "silac_medium", &mut s.silac_medium));
        }
        if s.silac != Silac::Off {
            row(ui, "Heavy labels", "Label masses of the heavy channel on K and R (Da); 0: not labelled.",
                |ui| silac_labels(ui, "silac_heavy", &mut s.silac_heavy));
            if s.crosslink.enabled {
                row(ui, "Mixed crosslinks", "Off (default): a crosslink joins peptides of one channel, since they \
                    come from the same cells. On: also light-heavy crosslinks, e.g. when lysates were mixed before \
                    crosslinking.",
                    |ui| ui.checkbox(&mut s.silac_mixed, "crosslinks between channels"));
            }
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
        row(ui, "Annotator app files", "Only for the separate Annotator app (github.com/snijderlab/annotator); \
            garlic's \"Annotated spectra\" tab needs nothing extra.\n\
            Writes results.annotator.csv (linear matches) and crosslinks.annotator.csv: target \
            matches at 1% FDR as ProForma (a crosslink is PEPK[XLMOD:02126#XL1]R//AK[#XL1]R) with \
            each spectrum's position in its file. Open them in Annotator together with the \
            spectrum files.",
            |ui| ui.checkbox(&mut s.write_annotator, "write .annotator.csv files"));
    });
}

/// K and R label masses of a SILAC channel: a preset or any mass
fn silac_labels(ui: &mut Ui, id: &str, masses: &mut (f64, f64)) -> egui::Response {
    ui.horizontal(|ui| {
        for (residue, mass) in [('K', &mut masses.0), ('R', &mut masses.1)] {
            ui.label(theme::bold(residue.to_string()));
            let shown = settings::SILAC_LABELS
                .iter()
                .find(|(_, r, m)| *r == residue && (*m - *mass).abs() < 1e-4)
                .map_or_else(|| if *mass > 0.0 { "custom".to_string() } else { "none".to_string() }, |(n, _, _)| n.to_string());
            ComboBox::from_id_salt((id, residue))
                .selected_text(shown)
                .width(150.0)
                .show_ui(ui, |ui| {
                    for (name, _, m) in settings::SILAC_LABELS.iter().filter(|(_, r, _)| *r == residue) {
                        ui.selectable_value(mass, *m, *name);
                    }
                    ui.selectable_value(mass, 0.0, "none");
                });
            theme::number(ui, mass, 0.0..=50.0, 1.0, |d| d.max_decimals(6).suffix(" Da"));
            ui.add_space(10.0);
        }
    })
    .response
}

const VIEW_TIP: &str = "Show the matches of this results folder on the \"Annotated spectra\" tab: \
    each spectrum annotated with Annotator's engine (built in), linear and crosslinked.";

fn crosslink_tab(ui: &mut Ui, s: &mut Settings) {
    heading(ui, "Crosslinking (XL-MS)");
    ui.checkbox(&mut s.crosslink.enabled, theme::bold("Search for crosslinked peptides"))
        .on_hover_text("Turns this into a crosslink search: results go to crosslinks.sage.tsv \
            (and, with 'Linear search' below, on by default, the usual bottom-up results: \
            results.sage.tsv, peptides.sage.tsv and proteins.sage.tsv). \
            Label-free and TMT quantification work (Scoring & output tab); chimera, DIA and parquet \
            output are not used in this mode. \
            Details: CROSSLINKING.md");
    note(ui, "Recommended settings: 2 missed cleavages, deisotoping on, precursor charge 3 to 6 \
        (Database & enzyme tab).");

    // two crosslink databases: no loop-links (each is within one protein)
    let two_databases = s.two_databases();
    let xl = &mut s.crosslink;
    ui.add_enabled_ui(xl.enabled, |ui| {
        grid(ui, "xl", |ui| {
            row(ui, "Reagent", "Crosslinker used.\n\
                • DSS / BS3: 138.06808 Da, non-cleavable\n• DSG / BS2G: 96.02113 Da, non-cleavable\n\
                • DSSO: 158.00376 Da, MS-cleavable (stubs 54.011 / 103.993 / 85.983)\n\
                • DSBU: 196.08479 Da, MS-cleavable (stubs 85.053 / 111.032)\n\
                • cliXlink: 267.05653 Da, MS-cleavable and asymmetric: the plain arm leaves DSSO-like stubs \
                (54.011 / 103.993 / 85.983), the alkyne-tagged arm 163.063 / 213.046 / 195.035\n\
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
                        xl.reset_stubs();
                    }
                    r
                });
            if xl.reagent == Reagent::Custom {
                row(ui, "Linker mass", "Mass the linker adds to the pair of peptides (Da). Example: 138.06808 for DSS; \
                    -18.010565 for a zero-length coupling (one water is lost).",
                    |ui| mass(ui, &mut xl.custom_mass));
                row(ui, "MS-cleavable", "The custom linker breaks in the collision cell and leaves stubs on each peptide. \
                    Enter the stub masses under 'Cleavable linker fragments' below.",
                    |ui| ui.checkbox(&mut xl.custom_cleavable, "cleavable linker with stub masses"));
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
            let cleavable = xl.is_cleavable();
            row(ui, "Loop-links", "Also search links with both ends on the same peptide. Not searched for \
                MS-cleavable reagents, nor in a two-database crosslink search (Files tab).", |ui| {
                let r = ui.add_enabled(
                    !cleavable && !two_databases,
                    egui::Checkbox::new(&mut xl.loop_links, "search loop-links"),
                );
                if two_databases {
                    ui.label(RichText::new("off: two-database search").weak());
                }
                r
            });
            if cleavable {
                row(ui, "Complete doublets", "Require both stub peaks for BOTH peptides. Stricter, fewer hits. \
                    Default: both stubs of one peptide and at least one of the other.",
                    |ui| ui.checkbox(&mut xl.require_complete_doublets, "require complete doublets"));
            }
            row(ui, "Min peaks per chain", "Each of the two peptides needs at least this many matched fragments. Default 3.",
                |ui| num(ui, &mut xl.min_chain_peaks, 1..=20, 1.0));
            row(ui, "Min linked ions per chain", "Non-cleavable linkers: each peptide needs at least this many matched                 fragments that carry the linker (ions containing the link site; for loop-links, both sites). They are                 the direct evidence of the link: without them, a correct peptide carrying any modification of about                 the linker's mass, or a strong peptide paired with a barely matched partner, is reported as a                 crosslink. Default 2; 0 turns the check off. MS-cleavable linkers rely on their stub doublets instead.",
                |ui| ui.add_enabled_ui(!cleavable, |ui| num(ui, &mut xl.min_linked_peaks, 0..=20, 1.0)).inner);
            row(ui, "Linear search", "Before the crosslink search, also run an ordinary linear peptide search on \
                the same spectra and database, giving the usual bottom-up results next to crosslinks.sage.tsv: \
                results.sage.tsv (PSMs, with peptide, protein and protein-group FDR), peptides.sage.tsv and \
                proteins.sage.tsv (peptides at 1% FDR and their protein groups); LFQ / TMT tables, pin file and HTML \
                report follow the Scoring & output tab. It uses the crosslink search's database settings (e.g. \
                minimum peptide mass 350 Da). Unlinked peptides, mono-links and the seeds of crosslinks all show \
                up here. On by default; roughly doubles the file reading time.",
                |ui| ui.checkbox(&mut xl.linear_search, "also run a linear peptide search (peptide and protein reports)"));
            row(ui, "Multivalent crosslinks", "Also search species held together by two or more linkers, up to the \
                number set here (2-4 linkers per species): several links between two peptides (a2b), chains and \
                branches of three or more peptides (a1b1c, a1b,b1c,b1d), and copies of one peptide (a2a, a1a1b, as \
                at a homodimer interface). Non-cleavable linkers: each peptide's ions carry whatever hangs off its \
                link residues, and every peptide needs its own linker-carrying ions. MS-cleavable linkers (DSSO, \
                DSBU, cliXlink): every linker breaks, each peptide keeps a stub on each link residue and must show \
                its peptide + stubs peaks (doublets, triplets); pieces freed by a single broken linker (b1c + stub) \
                tell the arrangement (topology_peaks). Written to crosslinks.sage.tsv as link_type 'multivalent' \
                with their topology and species, with their own FDR. The search space grows steeply: only for \
                small protein sets.",
                |ui| {
                    ui.horizontal(|ui| {
                        let mut r = ui.checkbox(&mut xl.multivalent, "search multivalent crosslinks");
                        if xl.multivalent {
                            ui.label("max linkers per species");
                            r |= num(ui, &mut xl.valency, 2..=sage_core::crosslink::MAX_VALENCY, 1.0);
                        }
                        r
                    })
                    .inner
                });
        });
        if xl.multivalent {
            let warn = Color32::from_rgb(220, 135, 30);
            ui.label(RichText::new(format!("⚠ {}", settings::multivalent_warning(xl.valency))).color(warn));
            ui.add_space(4.0);
        }

        if xl.is_cleavable() {
            stub_masses(ui, xl);
        }

        heading(ui, "Non-cleavable search engine");
        if xl.is_cleavable() {
            note(ui, "MS-cleavable reagents always use the CPU doublet search; these settings are ignored.");
        }
        grid(ui, "xl_backend", |ui| {
            row(ui, "Backend", "How non-cleavable crosslinks are searched:\n\
                • cpu (default): seed search + partner lookup; validated, fine for complexes and small databases\n\
                • reflect: chain reflection on the CPU; reference for the GPU, slow\n\
                • gpu: chain reflection on the graphics card, for proteome-scale searches (needs a sage.exe with GPU support, the default)\n\
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
                |ui| num(ui, &mut xl.chain_survivors, 100..=1_000_000, 1000.0));
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

/// Fragment (stub) masses an MS-cleavable linker leaves on each peptide, with
/// the doublet spacing and mass balance they imply
fn stub_masses(ui: &mut Ui, xl: &mut settings::Crosslink) {
    heading(ui, "Cleavable linker fragments");
    note(ui, "When the linker breaks in the collision cell, each peptide keeps a short or a long stub, so it \
        appears as a doublet. Preset values are filled in; edit them for a modified reagent or another stub form.");
    grid(ui, "xl_stubs", |ui| {
        if xl.reagent == Reagent::Custom {
            row(ui, "Asymmetric", "The two ends of the linker differ (e.g. cliXlink, where one arm carries an \
                alkyne tag), so each arm leaves its own short / long stubs. A cleavage leaves one arm's short stub \
                on one peptide and the other arm's long stub on its partner.",
                |ui| ui.checkbox(&mut xl.asymmetric, "the two arms leave different stubs"));
        }
        let two = xl.has_second_arm();
        let arm1 = if two { "Arm 1: " } else { "" };
        row(ui, &format!("{arm1}Short stub"), "Lighter remnant left on a peptide (Da). DSSO: alkene 54.01056; \
            DSBU: butyl amine 85.05276; cliXlink arm 1: alkene 54.01056.",
            |ui| mass(ui, &mut xl.stub_short));
        row(ui, &format!("{arm1}Long stub"), "Heavier remnant left on the other side of the cleaved bond (Da). \
            DSSO: sulfenic acid 103.99320; DSBU: butyl isocyanate 111.03203; cliXlink arm 1: sulfenic acid 103.99320.",
            |ui| mass(ui, &mut xl.stub_long));
        row(ui, &format!("{arm1}Additional stubs"), "Other observed stub forms, comma separated, used when scoring \
            each peptide. DSSO / cliXlink arm 1: thiol 85.98264 (the sulfenic acid minus water). Leave empty for none.",
            |ui| ui.add(TextEdit::singleline(&mut xl.stub_alt).hint_text("e.g. 85.98264").desired_width(200.0)));
        row(ui, &format!("{arm1}Doublet spacing"), "long - short: the mass gap between the two peaks of every \
            stub-peptide doublet on this arm.",
            |ui| ui.label(format!("{:.5} Da", xl.stub_long - xl.stub_short)));
        if two {
            row(ui, "Arm 2: Short stub", "Lighter remnant of the second arm (Da). cliXlink: alkyne-tagged alkene 163.06333.",
                |ui| mass(ui, &mut xl.stub2_short));
            row(ui, "Arm 2: Long stub", "Heavier remnant of the second arm (Da). cliXlink: alkyne-tagged sulfenic acid 213.04596.",
                |ui| mass(ui, &mut xl.stub2_long));
            row(ui, "Arm 2: Additional stubs", "Other stub forms of the second arm, comma separated. \
                cliXlink: alkyne-tagged thiol 195.03540. Leave empty for none.",
                |ui| ui.add(TextEdit::singleline(&mut xl.stub2_alt).hint_text("e.g. 195.0354").desired_width(200.0)));
            row(ui, "Arm 2: Doublet spacing", "long - short of the second arm.",
                |ui| ui.label(format!("{:.5} Da", xl.stub2_long - xl.stub2_short)));
        }
        let linker = xl.linker_mass();
        let ok = Color32::from_rgb(40, 160, 70);
        let bad = Color32::from_rgb(210, 60, 60);
        let sums = if two {
            vec![
                ("arm 1 short + arm 2 long", xl.stub_short + xl.stub2_long),
                ("arm 2 short + arm 1 long", xl.stub2_short + xl.stub_long),
            ]
        } else {
            vec![("short + long", xl.stub_short + xl.stub_long)]
        };
        row(ui, "Mass balance", "A cleavage splits the linker into a short and a long stub, so short + long \
            (for an asymmetric linker: each arm's short stub + the other arm's long stub) must equal the linker \
            mass, or the search will not start.", |ui| {
            ui.vertical(|ui| {
                for (what, sum) in sums {
                    let (text, color) = if (sum - linker).abs() <= 1e-3 {
                        (format!("{what} = {sum:.5} Da = linker mass ✔"), ok)
                    } else {
                        (format!("{what} = {sum:.5} Da ≠ linker mass {linker:.5} Da"), bad)
                    };
                    ui.label(RichText::new(text).color(color));
                }
            })
            .response
        });
        if xl.reagent.is_cleavable() {
            let name = match xl.reagent {
                Reagent::Clixlink => "cliXlink".to_string(),
                r => r.key().to_ascii_uppercase(),
            };
            row(ui, "Defaults", "Restore the preset stub masses of the selected reagent.", |ui| {
                let r = ui.add_enabled(xl.stubs_edited(), egui::Button::new(format!("Reset to {name} values")));
                if r.clicked() {
                    xl.reset_stubs();
                }
                r
            });
        }
    });
}

/// `(fill, border)` of the Progress tab: pulsing orange while the search runs,
/// green when it finished, red when it failed, grey when cancelled. The fill
/// is a translucent tint of the border colour (stronger when the tab is open),
/// so the label stays readable in both themes.
fn progress_tab_colors(ui: &Ui, state: progress::State, selected: bool) -> (Color32, Color32) {
    let (color, pulse) = match state {
        progress::State::Running => {
            let t = ui.ctx().input(|i| i.time) as f32;
            ui.ctx().request_repaint_after(std::time::Duration::from_millis(33));
            (Color32::from_rgb(235, 140, 30), 0.5 + 0.5 * (t * std::f32::consts::TAU / 1.4).sin())
        }
        progress::State::Finished => (Color32::from_rgb(40, 160, 70), 1.0),
        progress::State::Failed => (Color32::from_rgb(210, 60, 60), 1.0),
        progress::State::Cancelled => (Color32::GRAY, 1.0),
    };
    let tint = if selected { 0.38 } else { 0.24 };
    (
        color.gamma_multiply(tint * (0.45 + 0.55 * pulse)),
        color.gamma_multiply(0.3 + 0.7 * pulse),
    )
}

fn human_size(bytes: u64) -> String {
    match bytes {
        b if b >= 1 << 30 => format!("{:.1} GB", b as f64 / (1u64 << 30) as f64),
        b if b >= 1 << 20 => format!("{:.1} MB", b as f64 / (1u64 << 20) as f64),
        b if b >= 1 << 10 => format!("{:.0} KB", b as f64 / 1024.0),
        b => format!("{b} B"),
    }
}

/// Screenshots of settings tabs, rendered off-screen:
/// `SAGE_GUI_FIGURES=<folder> cargo test -p sage-gui screenshots -- --ignored`
#[cfg(test)]
mod screenshots {
    use super::*;
    use egui_kittest::Harness;

    #[test]
    #[ignore]
    fn screenshots() {
        let out = std::path::PathBuf::from(std::env::var("SAGE_GUI_FIGURES").expect("set SAGE_GUI_FIGURES"));
        std::fs::create_dir_all(&out).unwrap();
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.multivalent = true;
        s.crosslink.valency = 3;
        s.silac = Silac::Triplex;
        s.fasta = "C:/Users/vahen/source/repos/garlic/tests/PolII_Complex.fasta".into();
        s.labelled_proteins.insert("sp|P08518|RPB2_YEAST".into(), sage_core::labelled::Isotope::N15);
        s.labelled_proteins.insert("sp|P04050|RPB1_YEAST".into(), sage_core::labelled::Isotope::C13N15);
        let tabs: [(&str, fn(&mut Ui, &mut Settings)); 3] =
            [("database", database_tab), ("scoring", scoring_tab), ("crosslinking", crosslink_tab)];
        for (name, tab) in tabs {
            let mut settings = s.clone();
            // fonts take effect on the frame after they are set: draw from the second
            let mut installed = false;
            let mut harness = Harness::builder()
                .with_size(egui::vec2(980.0, 1900.0))
                .build_ui(move |ui| {
                    if !installed {
                        theme::install(ui.ctx(), false);
                        installed = true;
                        return;
                    }
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.horizontal(|ui| {
                            for t in ["Files", "Database & enzyme", "Crosslinking"] {
                                let _ = ui.selectable_label(t == "Crosslinking", theme::bold(t));
                            }
                        });
                        tab(ui, &mut settings)
                    });
                });
            harness.run();
            let image = harness.render().expect("render");
            let path = out.join(format!("gui_{name}.png"));
            image.save(&path).unwrap();
            println!("{}", path.display());
        }

        // the Files tab: with and without two crosslink databases
        for (name, two) in [("files", false), ("files_two_databases", true)] {
            let mut settings = s.clone();
            if two {
                settings.database_1 = "C:\\db\\PolII_core.fasta".into();
                settings.database_2 = "C:\\db\\TFIIH.fasta".into();
            }
            let mut app: Option<SageGui> = None;
            let mut harness = Harness::builder()
                .with_size(egui::vec2(980.0, 700.0))
                .build_ui(move |ui| {
                    let Some(app) = app.as_mut() else {
                        theme::install(ui.ctx(), false);
                        let mut new = SageGui::new(ui.ctx());
                        new.s = settings.clone();
                        app = Some(new);
                        return;
                    };
                    egui::ScrollArea::vertical().show(ui, |ui| app.files_tab(ui));
                });
            harness.run();
            let path = out.join(format!("gui_{name}.png"));
            harness.render().expect("render").save(&path).unwrap();
            println!("{}", path.display());
        }
    }

    /// The list of every FASTA protein starts hidden; the chosen ones are
    /// always shown, and the list opens on a click
    #[test]
    fn labelled_protein_list_starts_hidden() {
        use egui_kittest::kittest::Queryable;
        let fasta = std::env::temp_dir().join(format!("sage_gui_labelled_{}.fasta", std::process::id()));
        std::fs::write(
            &fasta,
            ">sp|P08518|RPB2_YEAST RNA polymerase II subunit RPB2 OS=Yeast\nMSDLANSEK\n\
             >sp|P16370|RPB3_YEAST RNA polymerase II subunit RPB3 OS=Yeast\nMSEEGPQVK\n\
             >sp|P34087|RPB7_YEAST RNA polymerase II subunit RPB7 OS=Yeast\nMFFIKDLSK\n",
        )
        .unwrap();
        let mut s = Settings::default();
        s.fasta = fasta.to_string_lossy().into_owned();
        // named by a part of the accession, as a configuration file may
        s.labelled_proteins.insert("P08518".into(), sage_core::labelled::Isotope::N15);
        // fonts take effect on the frame after they are set
        let mut installed = false;
        let mut harness = Harness::builder()
            .with_size(egui::vec2(900.0, 700.0))
            .build_ui(move |ui| {
                if !installed {
                    theme::install(ui.ctx(), false);
                    installed = true;
                    return;
                }
                labelled_proteins(ui, &mut s)
            });
        harness.run();
        assert!(harness.query_by_label("P08518").is_some(), "the chosen protein is shown");
        assert!(harness.query_by_label("RNA polymerase II subunit RPB2").is_some(), "with its name");
        assert!(harness.query_by_label("sp|P16370|RPB3_YEAST").is_none(), "the full list starts hidden");
        harness.get_by_label("Choose proteins (3 in the FASTA)").click();
        harness.run();
        assert!(harness.query_by_label("sp|P16370|RPB3_YEAST").is_some(), "the list opens");
        assert!(harness.query_by_label("sp|P34087|RPB7_YEAST").is_some());
        let _ = std::fs::remove_file(&fasta);
    }
}
