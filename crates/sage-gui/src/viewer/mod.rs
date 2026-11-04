//! "Annotated spectra" tab: the matches of a results folder and, for the one
//! picked, its spectrum annotated by Annotator's engine (see [`model`]).

pub mod figure;
pub mod loader;
pub mod mobiligram;
pub mod model;
pub mod species;

use eframe::egui::{self, Color32, RichText, Ui};
use egui_plot::{Bar, BarChart, Plot, PlotPoint, Text};
use model::{annotate, Annotation, Kind, Model, Peaks, Results, Series, SpectrumKey, StubIonKind};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver, TryRecvError};
use std::sync::Arc;

pub struct Viewer {
    pub folder: String,
    results: Option<Arc<Results>>,
    message: Option<(bool, String)>,
    spectra: HashMap<SpectrumKey, Peaks>,
    /// The results folder being read (in the background)
    opening: Option<Receiver<Result<Results, String>>>,
    /// Reads spectra in the background: the match picked, then the rows in view
    loader: Option<loader::Loader>,
    /// Spectra asked for (each once), and those not in their file
    requested: HashSet<SpectrumKey>,
    not_found: HashSet<SpectrumKey>,
    /// The spectrum file being read
    reading: Option<String>,
    show_linear: bool,
    show_crosslinks: bool,
    targets_only: bool,
    max_q: f32,
    search: String,
    /// Indices into `results.matches` passing the filters, best first
    shown: Vec<usize>,
    filter_key: String,
    selected: Option<usize>,
    model: Model,
    label_all: bool,
    export: ExportSettings,
    annotation: Option<Result<Annotation, String>>,
    annotated: Option<(usize, Model)>,
    mobiligram: MobiligramView,
    /// Width of the match list; the divider next to it is dragged to change it
    list_width: f32,
}

/// Match list width: at start (and on double-clicking the divider), and its limits
const LIST_WIDTH: f32 = 380.0;
const MIN_LIST_WIDTH: f32 = 220.0;
/// Room the spectrum keeps when the list is widened
const MIN_SPECTRUM_WIDTH: f32 = 420.0;
/// Height of the match list and its divider
const LIST_HEIGHT: f32 = 560.0;

impl Default for Viewer {
    fn default() -> Self {
        Viewer {
            folder: String::new(),
            results: None,
            message: None,
            spectra: HashMap::new(),
            opening: None,
            loader: None,
            requested: HashSet::new(),
            not_found: HashSet::new(),
            reading: None,
            show_linear: true,
            show_crosslinks: true,
            targets_only: true,
            max_q: 0.01,
            search: String::new(),
            shown: Vec::new(),
            filter_key: String::new(),
            selected: None,
            model: Model::CidHcd,
            label_all: false,
            export: ExportSettings::default(),
            annotation: None,
            annotated: None,
            mobiligram: MobiligramView::default(),
            list_width: LIST_WIDTH,
        }
    }
}

/// Colour of an ion series (unannotated peaks are grey)
fn series_color(series: Series) -> Color32 {
    match series {
        Series::B => Color32::from_rgb(31, 119, 180),
        Series::Y => Color32::from_rgb(214, 39, 40),
        Series::C => Color32::from_rgb(44, 160, 44),
        Series::Z => Color32::from_rgb(255, 127, 14),
        Series::A => Color32::from_rgb(148, 103, 189),
        Series::X => Color32::from_rgb(140, 86, 75),
        Series::Precursor => Color32::from_rgb(90, 90, 90),
        Series::Immonium => Color32::from_rgb(23, 190, 207),
        Series::Internal => Color32::from_rgb(227, 119, 194),
        Series::Other => Color32::from_rgb(188, 189, 34),
        Series::Stub(j) => stub_color(j as usize),
    }
}

/// Colours of linker stubs (DSSO: A teal, T amber, S magenta), apart from the
/// ion series colours
const STUB_COLORS: [Color32; 6] = [
    Color32::from_rgb(0, 150, 136),
    Color32::from_rgb(230, 159, 0),
    Color32::from_rgb(204, 0, 153),
    Color32::from_rgb(86, 180, 233),
    Color32::from_rgb(120, 160, 0),
    Color32::from_rgb(160, 82, 45),
];

fn stub_color(stub: usize) -> Color32 {
    STUB_COLORS[stub % STUB_COLORS.len()]
}

const UNMATCHED: Color32 = Color32::from_rgb(175, 175, 175);

impl Viewer {
    /// Read the results in `folder` (in the background); spectra are read
    /// as matches are picked
    pub fn open(&mut self, folder: &str) {
        self.folder = folder.to_string();
        self.results = None;
        self.spectra.clear();
        self.loader = None;
        self.requested.clear();
        self.not_found.clear();
        self.reading = None;
        self.selected = None;
        self.annotation = None;
        self.annotated = None;
        self.filter_key.clear();
        self.mobiligram.points = None;
        self.mobiligram.file = None;
        self.message = Some((false, "reading the results…".into()));
        let (tx, rx) = channel();
        let path = PathBuf::from(folder);
        let _ = std::thread::Builder::new().name("results".into()).spawn(move || {
            let _ = tx.send(Results::load(&path));
        });
        self.opening = Some(rx);
    }

    /// The results are read: list the matches, start the spectrum reader
    /// and show the best match
    fn loaded(&mut self, results: Results, ctx: &egui::Context) {
        let (linear, crosslinks) = results.matches.iter().fold((0, 0), |(l, x), m| match m.kind {
            Kind::Linear => (l + 1, x),
            _ => (l, x + 1),
        });
        self.message = Some((false, format!("{linear} linear and {crosslinks} crosslink matches")));
        // Annotator's CID model includes MS-cleavable linker (DSSO / DSBU) stub ions
        self.model = Model::CidHcd;
        let results = Arc::new(results);
        self.loader = Some(loader::Loader::start(results.clone(), ctx.clone()));
        self.results = Some(results);
        self.refilter();
        self.selected = self.shown.first().copied();
    }

    /// Ask for a spectrum, unless it was asked for already: `now` for the
    /// match picked, otherwise after it (rows in view)
    fn request(&mut self, keys: Vec<SpectrumKey>, now: bool) {
        let Some(loader) = &self.loader else { return };
        let keys: Vec<SpectrumKey> = keys
            .into_iter()
            .filter(|k| !self.spectra.contains_key(k) && self.requested.insert(k.clone()))
            .collect();
        if keys.is_empty() {
            return;
        }
        if now {
            keys.into_iter().for_each(|k| loader.send(loader::Request::Now(k)));
        } else {
            loader.send(loader::Request::Soon(keys));
        }
    }

    fn poll(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.opening {
            match rx.try_recv() {
                Ok(result) => {
                    self.opening = None;
                    match result {
                        Ok(results) => self.loaded(results, ctx),
                        Err(e) => self.message = Some((true, e)),
                    }
                }
                Err(TryRecvError::Empty) => ctx.request_repaint_after(std::time::Duration::from_millis(50)),
                Err(TryRecvError::Disconnected) => {
                    self.opening = None;
                    self.message = Some((true, "could not read the results".into()));
                }
            }
        }
        let Some(loader) = &self.loader else { return };
        while let Ok(event) = loader.events.try_recv() {
            match event {
                loader::Event::Spectra(spectra) => self.spectra.extend(spectra),
                loader::Event::NotFound(keys) => self.not_found.extend(keys),
                loader::Event::Failed(e) => self.message = Some((true, e)),
                loader::Event::Busy(file) => self.reading = file,
            }
        }
    }

    /// Wait for [`Viewer::open`] to finish (tests)
    #[cfg(test)]
    fn wait_until_open(&mut self) {
        if let Some(rx) = self.opening.take() {
            match rx.recv().unwrap() {
                Ok(results) => self.loaded(results, &egui::Context::default()),
                Err(e) => self.message = Some((true, e)),
            }
        }
    }

    /// Recompute the shown matches when a filter changed
    fn refilter(&mut self) {
        let key = format!(
            "{}{}{}{}{}",
            self.show_linear, self.show_crosslinks, self.targets_only, self.max_q, self.search
        );
        if key == self.filter_key {
            return;
        }
        self.filter_key = key;
        let Some(results) = &self.results else { return };
        let needle = self.search.to_ascii_uppercase();
        let mut shown: Vec<usize> = results
            .matches
            .iter()
            .enumerate()
            .filter(|(_, m)| match m.kind {
                Kind::Linear => self.show_linear,
                _ => self.show_crosslinks,
            })
            .filter(|(_, m)| !self.targets_only || m.target)
            .filter(|(_, m)| m.spectrum_q <= self.max_q)
            .filter(|(_, m)| {
                needle.is_empty()
                    || m.peptide.to_ascii_uppercase().contains(&needle)
                    || m.partner.to_ascii_uppercase().contains(&needle)
                    || m.proteins.to_ascii_uppercase().contains(&needle)
                    || m.scan.contains(&self.search)
            })
            .map(|(i, _)| i)
            .collect();
        shown.sort_by(|&a, &b| results.matches[b].score.total_cmp(&results.matches[a].score));
        // files read whole (mzML, MGF) keep the spectra of every shown match
        if let Some(loader) = &self.loader {
            loader.send(loader::Request::Shown(shown.iter().map(|&i| results.matches[i].key()).collect()));
        }
        self.shown = shown;
    }

    pub fn ui(&mut self, ui: &mut Ui, output_directory: &str) {
        self.poll(ui.ctx());
        ui.heading("Annotated spectra");
        ui.label(
            RichText::new(
                "Pick a match to see its spectrum annotated by Annotator's engine (mzcore / mzannotate, \
                 github.com/snijderlab/annotator), built into garlic. Spectra are read as Sage read them; \
                 for timsTOF that is the precursor's merged spectrum.",
            )
            .weak(),
        );
        ui.horizontal_wrapped(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            ui.label(
                RichText::new(
                    "When using these spectra or figures, please cite Annotator: Schulte et al., \
                     Anal. Chem. 2025, 97, 23120–23130, ",
                )
                .weak(),
            );
            ui.hyperlink_to("doi:10.1021/acs.analchem.5c02832", "https://doi.org/10.1021/acs.analchem.5c02832");
        });
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.label("Results folder");
            ui.add(egui::TextEdit::singleline(&mut self.folder).desired_width(460.0).hint_text("folder with results.json"));
            if ui.button("Browse…").clicked() {
                if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                    self.folder = dir.display().to_string();
                }
            }
            if ui
                .add_enabled(!output_directory.is_empty(), egui::Button::new("This search's folder"))
                .on_hover_text(output_directory)
                .clicked()
            {
                self.folder = output_directory.to_string();
            }
            if ui.button(RichText::new("Load").strong()).clicked() {
                let folder = self.folder.clone();
                self.open(&folder);
            }
        });
        if let Some((error, text)) = &self.message {
            let color = if *error { Color32::from_rgb(210, 60, 60) } else { ui.visuals().weak_text_color() };
            ui.horizontal(|ui| {
                if self.opening.is_some() {
                    ui.spinner();
                }
                ui.label(RichText::new(text).color(color));
            });
        }
        let Some(results) = self.results.clone() else { return };

        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut self.show_linear, "linear");
            ui.checkbox(&mut self.show_crosslinks, "crosslinks");
            ui.checkbox(&mut self.targets_only, "targets only")
                .on_hover_text("Crosslinks: target-target matches only.");
            ui.label("spectrum q ≤");
            crate::theme::number(ui, &mut self.max_q, 0.0..=1.0, 0.01, |d| d.max_decimals(3));
            ui.label("search");
            ui.add(egui::TextEdit::singleline(&mut self.search).desired_width(160.0).hint_text("peptide, protein, scan"));
            self.refilter();
            ui.label(RichText::new(format!("{} shown", self.shown.len())).weak());
            // spectra are read as matches are picked or come into view
            if let Some(file) = &self.reading {
                ui.spinner();
                ui.label(RichText::new(format!("reading spectra from {file}")).weak());
            } else if !self.spectra.is_empty() {
                ui.label(RichText::new(format!("{} spectra read", self.spectra.len())).weak()).on_hover_text(
                    "Spectra are read as you pick matches or scroll to them: timsTOF spectra one by one \
                     (as the search read them), mzML / MGF files whole, once.",
                );
            }
            ui.separator();
            let has_mobility = results.matches.iter().any(|m| m.ion_mobility.is_some());
            if ui
                .add_enabled(has_mobility, egui::Button::new("Mobiligram…"))
                .on_hover_text(
                    "Every precursor of the run by m/z and ion mobility (1/K0): marker shape by charge, \
                     colour by identification (crosslinks, linear peptides, unidentified). Export as SVG or PDF.",
                )
                .on_disabled_hover_text(
                    "Needs ion mobility: timsTOF (Bruker .d) data, searched with garlic 0.18 or later \
                     (crosslinks.sage.tsv with an ion_mobility column).",
                )
                .clicked()
            {
                self.mobiligram.open = true;
            }
        });
        ui.add_space(4.0);

        ui.horizontal_top(|ui| {
            // the list may not squeeze the spectrum below its minimum (e.g. after the window shrank)
            let widest = (ui.available_width() - MIN_SPECTRUM_WIDTH).max(MIN_LIST_WIDTH);
            self.list_width = self.list_width.clamp(MIN_LIST_WIDTH, widest);
            self.match_list(ui, &results);
            self.divider(ui, widest);
            ui.vertical(|ui| self.spectrum_view(ui, &results));
        });

        if self.mobiligram.open {
            let view = &mut self.mobiligram;
            let points = view.points.get_or_insert_with(|| Arc::new(mobiligram::points(&results))).clone();
            let mut open = true;
            egui::Window::new("Mobiligram")
                .open(&mut open)
                .default_size([920.0, 680.0])
                .resizable(true)
                .show(ui.ctx(), |ui| view.ui(ui, &points));
            view.open = open;
        }
    }

    /// The border between the match list and the spectrum: drag it to widen
    /// either side (up to `widest` for the list), double-click for the default
    fn divider(&mut self, ui: &mut Ui, widest: f32) {
        let (rect, response) = ui.allocate_exact_size(egui::vec2(10.0, LIST_HEIGHT), egui::Sense::click_and_drag());
        response.widget_info(|| egui::WidgetInfo::labeled(egui::WidgetType::Other, true, "Resize the match list"));
        let response = response
            .on_hover_cursor(egui::CursorIcon::ResizeHorizontal)
            .on_hover_text("Drag to resize the list; double-click to reset");
        if response.dragged() {
            ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeHorizontal);
            self.list_width = (self.list_width + response.drag_delta().x).clamp(MIN_LIST_WIDTH, widest);
        }
        if response.double_clicked() {
            self.list_width = LIST_WIDTH.min(widest);
        }
        let active = response.hovered() || response.dragged();
        let stroke = if active {
            egui::Stroke::new(2.0f32, crate::theme::accent(ui))
        } else {
            ui.visuals().widgets.noninteractive.bg_stroke
        };
        let x = rect.center().x;
        ui.painter().vline(x, rect.y_range(), stroke);
        // a grip in the middle, so the border reads as something to grab
        let grip = if active { crate::theme::accent(ui) } else { ui.visuals().weak_text_color() };
        for dy in [-8.0f32, 0.0, 8.0] {
            ui.painter().circle_filled(egui::pos2(x, rect.center().y + dy), 1.8, grip);
        }
    }

    fn match_list(&mut self, ui: &mut Ui, results: &Results) {
        let row_height = ui.text_style_height(&egui::TextStyle::Body) + 6.0;
        let list = egui::Layout::top_down(egui::Align::LEFT);
        let width = self.list_width;
        let mut in_view: Vec<SpectrumKey> = Vec::new();
        ui.allocate_ui_with_layout(egui::vec2(width, LIST_HEIGHT), list, |ui| {
            ui.set_width(width);
            // one line per match (rows must be of equal height): long titles end in "…",
            // shown whole on hover or by widening the list
            ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Truncate);
            egui::ScrollArea::vertical()
                .id_salt("viewer_matches")
                .auto_shrink([false, false])
                .max_height(LIST_HEIGHT)
                .show_rows(ui, row_height, self.shown.len(), |ui, rows| {
                    for &i in &self.shown[rows] {
                        let m = &results.matches[i];
                        in_view.push(m.key());
                        let tag = match m.kind {
                            Kind::Linear => "",
                            Kind::Crosslink => "XL ",
                            Kind::Loop => "loop ",
                            Kind::Multi => "multi ",
                        };
                        let text = format!("{tag}{}  {}+  q {:.3}", m.title(), m.charge, m.spectrum_q);
                        let response = ui
                            .selectable_label(self.selected == Some(i), text)
                            .on_hover_text(format!(
                                "{}\n{}\n{} scan {}\nscore {:.2}, spectrum q {:.4}, {} q {:.4}",
                                m.title(),
                                m.proteins,
                                m.file,
                                m.scan,
                                m.score,
                                m.spectrum_q,
                                if m.kind == Kind::Linear { "peptide" } else { "crosslink" },
                                m.q
                            ));
                        if response.clicked() {
                            self.selected = Some(i);
                        }
                    }
                });
        });
        // the spectra of the rows in view, after the one picked
        if let Some(i) = self.selected {
            self.request(vec![results.matches[i].key()], true);
        }
        self.request(in_view, false);
    }

    fn spectrum_view(&mut self, ui: &mut Ui, results: &Results) {
        let Some(i) = self.selected else {
            ui.label(RichText::new("Pick a match on the left.").weak());
            return;
        };
        let m = &results.matches[i];
        ui.horizontal_wrapped(|ui| {
            ui.label(crate::theme::bold(m.title()));
            let rank = if m.rank > 1 { format!(" · rank {}", m.rank) } else { String::new() };
            ui.label(
                RichText::new(format!(
                    "{}+ · {} scan {} · score {:.2} · spectrum q {:.4}{rank} · RT {:.2} min · \
                     mass {:.4} Da, {:+.1} ppm",
                    m.charge,
                    m.file,
                    m.scan,
                    m.score,
                    m.spectrum_q,
                    m.rt,
                    m.calcmass,
                    (m.expmass - m.calcmass) as f64 / m.calcmass.max(1.0) as f64 * 1e6
                ))
                .weak(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Fragmentation");
            egui::ComboBox::from_id_salt("viewer_model")
                .selected_text(self.model.label())
                .show_ui(ui, |ui| {
                    for model in Model::ALL {
                        ui.selectable_value(&mut self.model, model, model.label());
                    }
                })
                .response
                .on_hover_text(
                    "Annotator's fragmentation models. CID / HCD: b and y ions with neutral losses, and \
                     for MS-cleavable linkers (DSSO, DSBU) the cleaved stub ions. ETD, EThcD, EAD: c and z \
                     ions as well. All: every ion type (many coincidental matches).",
                );
            ui.checkbox(&mut self.label_all, "label all matched peaks");
            ui.label(RichText::new(format!("tolerance {}", tolerance_text(results.fragment_tol))).weak());
        });

        let key = m.key();
        let Some(peaks) = self.spectra.get(&key) else {
            if self.not_found.contains(&key) {
                ui.label(
                    RichText::new(format!("Spectrum {} is not in {} (was the file changed after the search?).", m.scan, m.file))
                        .color(Color32::from_rgb(210, 60, 60)),
                );
            } else {
                ui.horizontal(|ui| {
                    ui.spinner();
                    // the first spectrum of a run waits for the run to be opened
                    let text = if !self.spectra.keys().any(|(f, _)| f == &m.file) {
                        format!("Opening {} (about a second for a timsTOF run)…", m.file)
                    } else {
                        "Reading the spectrum…".to_string()
                    };
                    ui.label(RichText::new(text).weak());
                });
            }
            return;
        };
        if self.annotated != Some((i, self.model)) {
            let proforma = results.proforma(m);
            self.annotation = Some(annotate(peaks, &proforma, m.charge, self.model, results.fragment_tol, results.link(m).as_ref()));
            self.annotated = Some((i, self.model));
        }
        let Some(annotation) = &self.annotation else { return };
        let a = match annotation {
            Ok(a) => a,
            Err(e) => {
                ui.label(
                    RichText::new(format!("Could not annotate {}: {e}", results.proforma(m)))
                        .color(Color32::from_rgb(210, 60, 60)),
                );
                return;
            }
        };
        let mut precursor = peaks
            .precursor_mz
            .map_or(String::new(), |mz| format!(" · precursor m/z {mz:.4}"));
        if let Some(im) = peaks.ion_mobility {
            precursor.push_str(&format!(" · ion mobility 1/K0 {im:.4} V·s/cm²"));
        }
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new(&a.proforma).monospace());
            ui.label(
                RichText::new(format!(
                    "· {} of {} peaks annotated · {:.1}% of the intensity{precursor}",
                    a.explained_peaks,
                    a.peaks.len(),
                    100.0 * a.explained_intensity
                ))
                .weak(),
            );
        });
        let range = spectrum_plot(ui, a, self.label_all);
        coverage_view(ui, a);
        legend(ui, a);

        // Figure export
        ui.add_space(4.0);
        let mut export = false;
        ui.horizontal_wrapped(|ui| {
            ui.label(crate::theme::bold("Export figure"));
            egui::ComboBox::from_id_salt("viewer_export_format")
                .selected_text(self.export.format.label())
                .show_ui(ui, |ui| {
                    for f in figure::Format::ALL {
                        ui.selectable_value(&mut self.export.format, f, f.label());
                    }
                })
                .response
                .on_hover_text(
                    "SVG / PDF: vector graphics, editable in Illustrator / Inkscape, any size. \
                     PNG / TIFF: rendered at the chosen resolution, 7 inches wide.",
                );
            ui.add_enabled_ui(!self.export.format.is_vector(), |ui| {
                egui::ComboBox::from_id_salt("viewer_export_dpi")
                    .selected_text(format!("{} dpi", self.export.dpi))
                    .show_ui(ui, |ui| {
                        for dpi in [300, 600, 1200] {
                            ui.selectable_value(&mut self.export.dpi, dpi, format!("{dpi} dpi"));
                        }
                    });
            });
            ui.checkbox(&mut self.export.background, "white background")
                .on_hover_text("Off: transparent background (PNG and TIFF with an alpha channel).");
            ui.checkbox(&mut self.export.title, "title")
                .on_hover_text("Peptides, charge, scan, retention time and ion mobility above the spectrum.");
            ui.checkbox(&mut self.export.map, "peptide map");
            ui.checkbox(&mut self.export.zoom, "current zoom")
                .on_hover_text("Export the m/z range shown above; off: the whole spectrum.");
            export = ui.button("Save…").clicked();
        });
        if let Some((error, text)) = &self.export.message {
            let color = if *error { Color32::from_rgb(210, 60, 60) } else { ui.visuals().weak_text_color() };
            ui.label(RichText::new(text).color(color));
        }
        if export {
            let mut title = format!("{} · {}+ · {} scan {} · RT {:.2} min", m.title(), m.charge, m.file, m.scan, m.rt);
            if let Some(im) = peaks.ion_mobility {
                title.push_str(&format!(" · 1/K0 {im:.3} V·s/cm²"));
            }
            let options = figure::FigureOptions {
                background: self.export.background,
                title: self.export.title.then_some(title),
                label_all: self.label_all,
                mz_range: if self.export.zoom { range } else { None },
                coverage: self.export.map,
                dpi: self.export.dpi,
            };
            let svg = figure::svg(a, &options);
            let format = self.export.format;
            let name = figure_name(m, format.extension());
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name(&name)
                .add_filter(format.label(), &[format.extension()])
                .save_file()
            {
                self.export.message = Some(match figure::export(&svg, format, options.dpi, options.background, &path) {
                    Ok(()) => (false, format!("Saved {}", path.display())),
                    Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                });
            }
        }
        fragment_table(ui, a);
    }
}

/// Settings of the figure export
pub struct ExportSettings {
    format: figure::Format,
    dpi: u32,
    background: bool,
    title: bool,
    map: bool,
    zoom: bool,
    message: Option<(bool, String)>,
}

impl Default for ExportSettings {
    fn default() -> Self {
        ExportSettings {
            format: figure::Format::Svg,
            dpi: 600,
            background: true,
            title: true,
            map: true,
            zoom: true,
            message: None,
        }
    }
}

/// The mobiligram window: its precursors and export settings
struct MobiligramView {
    open: bool,
    /// Computed when the window first opens
    points: Option<Arc<Vec<mobiligram::Point>>>,
    /// Only this spectrum file (`None`: all)
    file: Option<String>,
    unidentified: bool,
    format: figure::Format,
    background: bool,
    title: bool,
    message: Option<(bool, String)>,
}

impl Default for MobiligramView {
    fn default() -> Self {
        MobiligramView {
            open: false,
            points: None,
            file: None,
            unidentified: true,
            format: figure::Format::Svg,
            background: true,
            title: true,
            message: None,
        }
    }
}

/// egui_plot's marker for a charge state's shape
fn plot_marker(shape: mobiligram::Shape) -> egui_plot::MarkerShape {
    use egui_plot::MarkerShape as M;
    use mobiligram::Shape;
    match shape {
        Shape::Cross => M::Plus,
        Shape::Circle => M::Circle,
        Shape::Triangle => M::Up,
        Shape::Square => M::Square,
        Shape::Diamond => M::Diamond,
        Shape::TriangleDown => M::Down,
    }
}

/// A charge-state marker drawn with the painter (for the key under the plot)
fn paint_marker(painter: &egui::Painter, c: egui::Pos2, shape: mobiligram::Shape, r: f32, stroke: egui::Stroke) {
    use mobiligram::Shape;
    let p = |x: f32, y: f32| egui::pos2(c.x + x * r, c.y + y * r);
    match shape {
        Shape::Circle => {
            painter.circle_stroke(c, r, stroke);
        }
        Shape::Square => {
            painter.rect_stroke(egui::Rect::from_center_size(c, egui::vec2(1.7 * r, 1.7 * r)), 0.0, stroke, egui::StrokeKind::Middle);
        }
        Shape::Diamond => {
            painter.add(egui::Shape::closed_line(vec![p(0.0, -1.15), p(1.15, 0.0), p(0.0, 1.15), p(-1.15, 0.0)], stroke));
        }
        Shape::Triangle => {
            painter.add(egui::Shape::closed_line(vec![p(0.0, -1.2), p(1.04, 0.6), p(-1.04, 0.6)], stroke));
        }
        Shape::TriangleDown => {
            painter.add(egui::Shape::closed_line(vec![p(0.0, 1.2), p(1.04, -0.6), p(-1.04, -0.6)], stroke));
        }
        Shape::Cross => {
            painter.line_segment([p(-1.0, 0.0), p(1.0, 0.0)], stroke);
            painter.line_segment([p(0.0, -1.0), p(0.0, 1.0)], stroke);
        }
    }
}

impl MobiligramView {
    fn ui(&mut self, ui: &mut Ui, points: &[mobiligram::Point]) {
        use mobiligram::{Class, Shape};
        let mut files: Vec<&str> = points.iter().map(|p| p.file.as_str()).collect();
        files.sort_unstable();
        files.dedup();
        let options = mobiligram::MobiligramOptions {
            background: self.background,
            title: None,
            unidentified: self.unidentified,
            file: self.file.clone(),
        };
        let shown = mobiligram::selected(points, &options);

        ui.horizontal_wrapped(|ui| {
            if files.len() > 1 {
                ui.label("file");
                egui::ComboBox::from_id_salt("mobiligram_file")
                    .selected_text(self.file.as_deref().unwrap_or("all files"))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.file, None, "all files");
                        for f in &files {
                            ui.selectable_value(&mut self.file, Some(f.to_string()), *f);
                        }
                    });
            }
            ui.checkbox(&mut self.unidentified, "unidentified precursors")
                .on_hover_text("Also show the precursors with no match at 1% FDR (grey).");
            let crosslinks = shown.iter().filter(|p| p.class.is_crosslink()).count();
            let linear = shown.iter().filter(|p| p.class == Class::Linear).count();
            ui.label(
                RichText::new(format!("{} precursors · {crosslinks} crosslinks · {linear} linear peptides", shown.len())).weak(),
            );
        });

        let plot_height = (ui.available_height() - 70.0).max(320.0);
        Plot::new("mobiligram_plot")
            .height(plot_height)
            .legend(egui_plot::Legend::default().position(egui_plot::Corner::RightTop))
            .x_axis_label("m/z")
            .y_axis_label("1/K0 (V·s/cm²)")
            .label_formatter(|name, p| {
                let name = if name.is_empty() { String::new() } else { format!("{name}\n") };
                format!("{name}m/z {:.3}\n1/K0 {:.3}", p.x, p.y)
            })
            .show(ui, |plot_ui| {
                for class in Class::ALL {
                    let color = match class {
                        Class::Unidentified => class.color().gamma_multiply(0.6),
                        Class::Linear => class.color().gamma_multiply(0.5),
                        _ => class.color(),
                    };
                    for charge in 1..=6u8 {
                        let xy: Vec<[f64; 2]> = shown
                            .iter()
                            .filter(|p| p.class == class && p.charge.clamp(1, 6) == charge)
                            .map(|p| [p.mz, p.mobility as f64])
                            .collect();
                        if xy.is_empty() {
                            continue;
                        }
                        plot_ui.points(
                            egui_plot::Points::new(xy)
                                .name(class.label())
                                .color(color)
                                .shape(plot_marker(Shape::of(charge)))
                                .radius(if class.is_crosslink() { 4.0f32 } else { 2.5 }),
                        );
                    }
                }
            });

        // charge key: the shapes in use
        let mut charges: Vec<u8> = shown.iter().map(|p| p.charge.clamp(1, 6)).collect();
        charges.sort_unstable();
        charges.dedup();
        ui.horizontal(|ui| {
            ui.label("charge:");
            let stroke = egui::Stroke::new(1.3f32, ui.visuals().text_color());
            for charge in charges {
                let (rect, _) = ui.allocate_exact_size(egui::vec2(14.0, 14.0), egui::Sense::hover());
                paint_marker(ui.painter(), rect.center(), Shape::of(charge), 5.0, stroke);
                ui.label(if charge == 6 { "6+ and up".to_string() } else { format!("{charge}+") });
                ui.add_space(6.0);
            }
        });

        // export
        let mut save = false;
        ui.horizontal_wrapped(|ui| {
            ui.label(crate::theme::bold("Export"));
            egui::ComboBox::from_id_salt("mobiligram_format")
                .selected_text(self.format.label())
                .show_ui(ui, |ui| {
                    for f in [figure::Format::Svg, figure::Format::Pdf] {
                        ui.selectable_value(&mut self.format, f, f.label());
                    }
                });
            ui.checkbox(&mut self.background, "white background").on_hover_text("Off: transparent background.");
            ui.checkbox(&mut self.title, "title");
            save = ui.button("Save…").clicked();
        });
        if let Some((error, text)) = &self.message {
            let color = if *error { Color32::from_rgb(210, 60, 60) } else { ui.visuals().weak_text_color() };
            ui.label(RichText::new(text).color(color));
        }
        if save {
            // the names a search uses: `<MS file>_MG` (also for all files when
            // there is only one), or `Combined_Mobiligram` for several
            let file = mobiligram::figure_file(points, self.file.as_deref());
            let name = mobiligram::figure_name(file.as_deref());
            let options = mobiligram::MobiligramOptions { file, ..options };
            let options = mobiligram::MobiligramOptions {
                title: self.title.then(|| mobiligram::title(points, &options)),
                ..options
            };
            let svg = mobiligram::svg(points, &options);
            let format = self.format;
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name(format!("{name}.{}", format.extension()))
                .add_filter(format.label(), &[format.extension()])
                .save_file()
            {
                self.message = Some(match figure::export(&svg, format, 300, self.background, &path) {
                    Ok(()) => (false, format!("Saved {}", path.display())),
                    Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                });
            }
        }
    }
}

/// File name for a figure: `GLEKVSTAAMQR-TDYHFPNEWKLIR_scan4321.svg`
fn figure_name(m: &model::Match, extension: &str) -> String {
    let plain = |p: &str| {
        let mut depth = 0;
        p.chars()
            .filter(|c| {
                match c {
                    '[' => depth += 1,
                    ']' => depth -= 1,
                    _ => return depth == 0 && c.is_ascii_alphabetic(),
                }
                false
            })
            .collect::<String>()
    };
    let peptides = if m.partner.is_empty() { plain(&m.peptide) } else { format!("{}-{}", plain(&m.peptide), plain(&m.partner)) };
    let scan: String = m.scan.chars().filter(|c| c.is_ascii_alphanumeric()).collect();
    let scan = if scan.len() > 24 { scan[scan.len() - 24..].to_string() } else { scan };
    format!("{peptides}_scan{scan}.{extension}")
}

fn tolerance_text(t: sage_core::mass::Tolerance) -> String {
    match t {
        sage_core::mass::Tolerance::Ppm(lo, hi) => format!("±{} ppm", lo.abs().max(hi.abs())),
        sage_core::mass::Tolerance::Da(lo, hi) => format!("±{} Da", lo.abs().max(hi.abs())),
        sage_core::mass::Tolerance::Pct(lo, hi) => format!("±{}%", lo.abs().max(hi.abs())),
    }
}

/// Peak label as shown: `y7`, or `2:b5` for the second peptide of a crosslink
fn label_text(a: &Annotation, label: &model::IonLabel) -> String {
    // parts of a multivalent species (`{2,3}+A`) name their peptides
    if a.coverage.len() > 1 && !label.text.starts_with('{') {
        format!("{}:{}", label.peptide + 1, label.text)
    } else {
        label.text.clone()
    }
}

/// Draw the annotated spectrum; returns the m/z range shown (the zoom)
fn spectrum_plot(ui: &mut Ui, a: &Annotation, label_all: bool) -> Option<(f64, f64)> {
    let max_mz = a.peaks.iter().map(|p| p.mz).fold(0.0, f64::max).max(1.0);
    let width = max_mz / 1500.0;
    let mut by_series: HashMap<Option<Series>, Vec<Bar>> = HashMap::new();
    for p in &a.peaks {
        let series = p.labels.first().map(|l| l.series);
        by_series
            .entry(series)
            .or_default()
            .push(Bar::new(p.mz, p.intensity as f64).width(width));
    }
    // Labels: matched peaks above 2% of the base peak, or all
    let threshold = if label_all { 0.0 } else { 0.02 * a.max_intensity };
    let peaks = &a.peaks;
    let hover = |_: &str, point: &PlotPoint| -> String {
        // the peak nearest to the pointer
        let nearest = peaks
            .iter()
            .min_by(|x, y| (x.mz - point.x).abs().total_cmp(&(y.mz - point.x).abs()));
        match nearest {
            Some(p) if (p.mz - point.x).abs() < max_mz / 200.0 => {
                let mut text = format!("m/z {:.4}\nintensity {:.3e}", p.mz, p.intensity);
                for l in &p.labels {
                    let ppm = (p.mz - l.theoretical_mz) / l.theoretical_mz * 1e6;
                    text.push_str(&format!("\n{}  ({:.4}, {:+.1} ppm)", label_text(a, l), l.theoretical_mz, ppm));
                }
                text
            }
            _ => format!("m/z {:.2}", point.x),
        }
    };
    let mut series: Vec<(Option<Series>, Vec<Bar>)> = by_series.into_iter().collect();
    // unmatched first, so matched peaks are drawn on top
    series.sort_by_key(|(s, _)| s.map(|s| s.name()));
    let shown = Plot::new("viewer_spectrum")
        .height(360.0)
        .allow_scroll(false)
        .x_axis_label("m/z")
        .y_axis_label("intensity")
        .include_y(0.0)
        .include_y(a.max_intensity as f64 * 1.12)
        .label_formatter(hover)
        .show(ui, |plot| {
            for (s, bars) in series {
                let (color, name) = s.map_or((UNMATCHED, "unmatched"), |s| (series_color(s), s.name()));
                plot.bar_chart(BarChart::new(bars).color(color).name(name));
            }
            for p in a.peaks.iter().filter(|p| !p.labels.is_empty() && p.intensity >= threshold) {
                let l = &p.labels[0];
                let mut text = label_text(a, l);
                if p.labels.len() > 1 {
                    text.push('*');
                }
                plot.text(
                    Text::new(PlotPoint::new(p.mz, p.intensity as f64), RichText::new(text).size(11.0))
                        .color(series_color(l.series))
                        .anchor(egui::Align2::CENTER_BOTTOM),
                );
            }
        });
    ui.label(
        RichText::new("Drag to pan, ctrl + scroll to zoom, double-click to reset. * = more than one interpretation (hover).")
            .weak()
            .small(),
    );
    let bounds = shown.transform.bounds();
    Some((bounds.min()[0].max(0.0), bounds.max()[0]))
}

/// Width and height of one residue in the coverage map
const CELL: egui::Vec2 = egui::vec2(18.0, 30.0);
/// Room above and below a peptide row for nested stub hooks
const PAD: f32 = 12.0;
/// Width of the "peptide 1" column
const LABEL_W: f32 = 76.0;

/// A hook at a residue corner: N-terminal ions over its right (ending after
/// it), C-terminal ions under its left (starting at it). `level` 0 is a plain
/// ion; stub ions nest outward at levels 1, 2, ...
fn hook(painter: &egui::Painter, r: egui::Rect, n_term: bool, level: usize, color: Color32) {
    let out = 3.0 * level as f32;
    let stroke = egui::Stroke::new(2.0_f32, color);
    let (corner, down, back) = if n_term {
        (egui::pos2(r.right() + out * 0.6, r.top() + 3.0 - out), egui::vec2(0.0, 9.0 + out), egui::vec2(-6.0 - out * 0.6, 0.0))
    } else {
        (egui::pos2(r.left() - out * 0.6, r.bottom() - 3.0 + out), egui::vec2(0.0, -9.0 - out), egui::vec2(6.0 + out * 0.6, 0.0))
    };
    painter.line_segment([corner, corner + down], stroke);
    painter.line_segment([corner, corner + back], stroke);
}

/// One peptide row starting at `origin` (top left of its first residue):
/// residues, plain backbone hooks and the hooks of stub-carrying ions
fn draw_peptide(painter: &egui::Painter, origin: egui::Pos2, a: &Annotation, peptide: usize, linked: &[usize], text_color: Color32) {
    let c = &a.coverage[peptide];
    let font = egui::FontId::monospace(16.0);
    for (i, residue) in c.residues.iter().enumerate() {
        let r = egui::Rect::from_min_size(origin + egui::vec2(CELL.x * i as f32, 0.0), CELL);
        if linked.contains(&i) {
            painter.rect_filled(r.shrink2(egui::vec2(1.0, 6.0)), 3.0, text_color.gamma_multiply(0.18));
        }
        painter.text(r.center(), egui::Align2::CENTER_CENTER, residue, font.clone(), text_color);
        if c.n_term[i] {
            hook(painter, r, true, 0, series_color(Series::B));
        }
        if c.c_term[i] {
            hook(painter, r, false, 0, series_color(Series::Y));
        }
        // stub-carrying ions, one nested hook per stub
        for n_term in [true, false] {
            let mut stubs: Vec<usize> = a
                .stub_ions
                .iter()
                .filter(|s| s.peptide == peptide)
                .filter(|s| match s.kind {
                    StubIonKind::NTerm(j) => n_term && j == i,
                    StubIonKind::CTerm(j) => !n_term && j == i,
                    StubIonKind::Peptide | StubIonKind::Part { .. } => false,
                })
                .flat_map(|s| s.stubs.iter().copied())
                .collect();
            stubs.sort_unstable();
            stubs.dedup();
            for (level, stub) in stubs.into_iter().enumerate() {
                hook(painter, r, n_term, level + 1, stub_color(stub));
            }
        }
    }
}

/// Filled triangle with its tip toward `up` (the peptide carrying the stub)
fn triangle(painter: &egui::Painter, center: egui::Pos2, up: bool, color: Color32) {
    let (tip, base) = if up { (-5.0, 4.0) } else { (5.0, -4.0) };
    let points = vec![
        center + egui::vec2(0.0, tip),
        center + egui::vec2(-5.0, base),
        center + egui::vec2(5.0, base),
    ];
    painter.add(egui::Shape::convex_polygon(points, color, egui::Stroke::NONE));
}

fn backbone_text(a: &Annotation, peptide: usize) -> String {
    let c = &a.coverage[peptide];
    let seen = c.n_term.iter().zip(&c.c_term).filter(|(n, c)| **n || **c).count();
    format!("{seen} of {} positions with a backbone ion", c.residues.len())
}

/// The stub ions found, per peptide and stubs carried, for the tooltip of
/// the map; then the parts of a multivalent species freed by one linker
fn stub_summary(a: &Annotation, link: &model::Link) -> String {
    let ion_text = |s: &model::StubIon| {
        let observed = a.peaks[s.peak].mz;
        let ppm = (observed - s.theoretical_mz) / s.theoretical_mz * 1e6;
        format!("{} {observed:.4} ({ppm:+.1} ppm)", s.text)
    };
    let mut lines = Vec::new();
    for peptide in 0..a.coverage.len() {
        let mut groups: Vec<(&Vec<usize>, Vec<String>)> = Vec::new();
        for s in a.stub_ions.iter().filter(|s| s.peptide == peptide && !matches!(s.kind, StubIonKind::Part { .. })) {
            match groups.iter_mut().find(|(stubs, _)| *stubs == &s.stubs) {
                Some((_, ions)) => ions.push(ion_text(s)),
                None => groups.push((&s.stubs, vec![ion_text(s)])),
            }
        }
        groups.sort_by(|x, y| x.0.cmp(y.0));
        for (stubs, ions) in groups {
            let what = match stubs.as_slice() {
                [j] => {
                    let stub = &link.stubs[*j];
                    format!("{} ({}, {:+.4} Da)", stub.name, stub.description, stub.mass)
                }
                many => many.iter().map(|&j| link.stubs[j].name.as_str()).collect::<Vec<_>>().join("+"),
            };
            lines.push(format!("peptide {} + {what}: {}", peptide + 1, ions.join(", ")));
        }
    }
    for s in &a.stub_ions {
        if let StubIonKind::Part { link: e, .. } = s.kind {
            lines.push(format!("XL{} broken, the rest held together: {}", e + 1, ion_text(s)));
        }
    }
    if lines.is_empty() {
        "No stub ions found in this spectrum.".into()
    } else {
        lines.join("
")
    }
}

/// Residues of each peptide with the backbone ions seen, drawn as in
/// Annotator: a hook over the right of a residue for an N-terminal (a/b/c)
/// ion ending there, one under the left for a C-terminal (x/y/z) ion starting
/// there. Crosslinked pairs are joined by a bridge between their linked
/// residues; loop-links by a bracket.
fn coverage_view(ui: &mut Ui, a: &Annotation) {
    match &a.link {
        Some(link) if link.kind == Kind::Crosslink && a.coverage.len() == 2 => crosslink_map(ui, a, link),
        Some(link) if link.kind == Kind::Loop && a.coverage.len() == 1 => loop_map(ui, a, link),
        Some(link) if link.kind == Kind::Multi && a.coverage.len() >= 2 && !link.links.is_empty() => {
            species_map(ui, a, link)
        }
        _ => {
            for k in 0..a.coverage.len() {
                ui.horizontal(|ui| {
                    let width = CELL.x * a.coverage[k].residues.len() as f32;
                    let (rect, _) = ui.allocate_exact_size(egui::vec2(width, CELL.y + 2.0 * PAD), egui::Sense::hover());
                    let painter = ui.painter_at(rect.expand(4.0));
                    draw_peptide(&painter, rect.min + egui::vec2(0.0, PAD), a, k, &[], ui.visuals().text_color());
                    ui.label(RichText::new(format!("  {}", backbone_text(a, k))).weak());
                });
            }
        }
    }
}

/// Two peptides, their linked residues aligned and joined by a bridge. On the
/// bridge, a triangle per stub found on the whole peptide (the doublet
/// peaks), pointing to the peptide that carries it; around the residues,
/// stub-carrying fragment ions as nested hooks in the stub's colour.
fn crosslink_map(ui: &mut Ui, a: &Annotation, link: &model::Link) {
    let (s0, s1) = (link.sites[0].1, link.sites[1].1);
    let shift = s0.max(s1);
    let (o0, o1) = ((shift - s0) as f32, (shift - s1) as f32);
    let (n0, n1) = (a.coverage[0].residues.len() as f32, a.coverage[1].residues.len() as f32);
    let cells = (o0 + n0).max(o1 + n1);
    let row_h = CELL.y + 2.0 * PAD;
    let bridge_h = 54.0;
    let size = egui::vec2(LABEL_W + CELL.x * cells + 330.0, 2.0 * row_h + bridge_h);
    let (rect, response) = ui.allocate_exact_size(size, egui::Sense::hover());
    let painter = ui.painter_at(rect.expand(4.0));
    let text_color = ui.visuals().text_color();
    let weak = ui.visuals().weak_text_color();
    let small = egui::FontId::proportional(12.0);

    let x0 = rect.left() + LABEL_W;
    let top0 = rect.top() + PAD;
    let top1 = rect.top() + row_h + bridge_h + PAD;
    for (k, top, offset) in [(0, top0, o0), (1, top1, o1)] {
        painter.text(
            egui::pos2(rect.left(), top + CELL.y / 2.0),
            egui::Align2::LEFT_CENTER,
            format!("peptide {}", k + 1),
            egui::FontId::proportional(14.0),
            text_color,
        );
        let origin = egui::pos2(x0 + CELL.x * offset, top);
        draw_peptide(&painter, origin, a, k, &[link.sites[k].1], text_color);
        let end = x0 + CELL.x * cells + 14.0;
        painter.text(egui::pos2(end, top + CELL.y / 2.0), egui::Align2::LEFT_CENTER, backbone_text(a, k), small.clone(), weak);
    }

    // the bridge between the linked residues
    let bx = x0 + CELL.x * (shift as f32 + 0.5);
    let (y_top, y_bottom) = (top0 + CELL.y - 2.0, top1 + 2.0);
    painter.line_segment([egui::pos2(bx, y_top), egui::pos2(bx, y_bottom)], egui::Stroke::new(3.0_f32, weak));
    painter.text(
        egui::pos2(bx - 10.0, (y_top + y_bottom) / 2.0),
        egui::Align2::RIGHT_CENTER,
        &link.linker,
        egui::FontId::proportional(13.0),
        text_color,
    );

    // stubs seen on the whole peptides: near each end of the bridge
    for (k, y, up) in [(0, y_top + 13.0, true), (1, y_bottom - 13.0, false)] {
        let mut stubs: Vec<usize> = a
            .stub_ions
            .iter()
            .filter(|s| s.peptide == k && s.kind == StubIonKind::Peptide)
            .map(|s| s.stub)
            .collect();
        stubs.sort_unstable();
        stubs.dedup();
        let mut x = bx + 12.0;
        for stub in stubs {
            let color = stub_color(stub);
            triangle(&painter, egui::pos2(x, y), up, color);
            let name = painter.text(egui::pos2(x + 8.0, y), egui::Align2::LEFT_CENTER, &link.stubs[stub].name, small.clone(), color);
            x = name.right() + 12.0;
        }
    }
    if !link.stubs.is_empty() {
        response.on_hover_text(stub_summary(a, link));
    }
}

/// Multivalent species: the stubs found at the `peptide` end of link `e`:
/// the parts freed when only that linker broke, and a peptide with no other
/// link with its stub (its doublet)
pub(super) fn bridge_end_stubs(a: &Annotation, link: &model::Link, e: usize, peptide: usize) -> Vec<usize> {
    let single = link.ends_of(peptide).len() == 1;
    let mut stubs: Vec<usize> = a
        .stub_ions
        .iter()
        .filter(|s| s.peptide == peptide)
        .filter(|s| match s.kind {
            StubIonKind::Part { link, .. } => link as usize == e,
            StubIonKind::Peptide => single,
            _ => false,
        })
        .map(|s| s.stub)
        .collect();
    stubs.sort_unstable();
    stubs.dedup();
    stubs
}

/// Multivalent species: the stub sets found on a peptide with two or more
/// link residues, freed whole (its triplet peaks)
pub(super) fn peptide_stub_sets(a: &Annotation, link: &model::Link, peptide: usize) -> Vec<Vec<usize>> {
    if link.ends_of(peptide).len() < 2 {
        return Vec::new();
    }
    let mut sets: Vec<Vec<usize>> = a
        .stub_ions
        .iter()
        .filter(|s| s.peptide == peptide && s.kind == StubIonKind::Peptide)
        .map(|s| s.stubs.clone())
        .collect();
    sets.sort();
    sets.dedup();
    sets
}

/// Height of the gap between two peptide rows (the bridges)
const GAP: f32 = 54.0;
/// Spacing of the lanes of bridges between rows that are not neighbours
const LANE: f32 = 12.0;

/// A multivalent species: its peptides as rows, each link residue shaded
/// and joined to its partner by a bridge (`XLn`, as `#XLn` in the ProForma).
/// MS-cleavable linkers: on each bridge, a triangle per stub found where
/// that linker broke (pointing at the side carrying it), after a peptide
/// with two or more links its peptide + stubs peaks, and around the residues
/// the stub-carrying fragment ions as nested hooks.
fn species_map(ui: &mut Ui, a: &Annotation, link: &model::Link) {
    let lengths: Vec<usize> = a.coverage.iter().map(|c| c.residues.len()).collect();
    let l = species::layout(&lengths, &link.links);
    let rows = l.rows.len() as f32;
    let row_h = CELL.y + 2.0 * PAD;
    let lanes_w = if l.lanes > 0 { 14.0 + LANE * l.lanes as f32 } else { 0.0 };
    let size = egui::vec2(LABEL_W + CELL.x * l.cells as f32 + lanes_w + 380.0, rows * row_h + (rows - 1.0) * GAP);
    ui.label(
        RichText::new(format!(
            "{} {} linkers: XL1–XL{} as #XL1… in the ProForma above",
            link.links.len(),
            link.linker,
            link.links.len()
        ))
        .weak()
        .small(),
    );
    let (rect, response) = ui.allocate_exact_size(size, egui::Sense::hover());
    let painter = ui.painter_at(rect.expand(4.0));
    let text_color = ui.visuals().text_color();
    let weak = ui.visuals().weak_text_color();
    let small = egui::FontId::proportional(12.0);

    let x0 = rect.left() + LABEL_W;
    let top = |row: usize| rect.top() + PAD + row as f32 * (row_h + GAP);
    let x_of = |row: usize, residue: usize| x0 + CELL.x * (l.column(row, residue) as f32 + 0.5);
    let cells_end = x0 + CELL.x * l.cells as f32;
    let text_x = cells_end + lanes_w + 14.0;
    for (row, &k) in l.rows.iter().enumerate() {
        let y = top(row);
        painter.text(
            egui::pos2(rect.left(), y + CELL.y / 2.0),
            egui::Align2::LEFT_CENTER,
            format!("peptide {}", k + 1),
            egui::FontId::proportional(14.0),
            text_color,
        );
        let origin = egui::pos2(x0 + CELL.x * l.offsets[row] as f32, y);
        draw_peptide(&painter, origin, a, k, &link.ends_of(k), text_color);
        // a peptide with two or more links, freed whole with its stubs
        let mut x = text_x;
        for set in peptide_stub_sets(a, link, k) {
            for &stub in &set {
                triangle(&painter, egui::pos2(x + 5.0, y + CELL.y / 2.0), true, stub_color(stub));
                x += 11.0;
            }
            let names = set.iter().map(|&j| link.stubs[j].name.as_str()).collect::<Vec<_>>().join("+");
            let r = painter.text(egui::pos2(x + 2.0, y + CELL.y / 2.0), egui::Align2::LEFT_CENTER, names, small.clone(), weak);
            x = r.right() + 12.0;
        }
        painter.text(egui::pos2(x, y + CELL.y / 2.0), egui::Align2::LEFT_CENTER, backbone_text(a, k), small.clone(), weak);
    }

    let stroke = egui::Stroke::new(3.0_f32, weak);
    for b in &l.bridges {
        let (xu, yu) = (x_of(b.upper.0, b.upper.1), top(b.upper.0) + CELL.y - 2.0);
        let (xl, yl) = (x_of(b.lower.0, b.lower.1), top(b.lower.0) + 2.0);
        let name = format!("XL{}", b.link + 1);
        // where the bridge is near each end
        let (x_upper_end, x_lower_end) = match b.lane {
            None => {
                painter.line_segment([egui::pos2(xu, yu), egui::pos2(xl, yl)], stroke);
                let mid = egui::pos2((xu + xl) / 2.0, (yu + yl) / 2.0);
                painter.text(mid - egui::vec2(9.0, 0.0), egui::Align2::RIGHT_CENTER, &name, egui::FontId::proportional(13.0), text_color);
                let at = |y: f32| xu + (xl - xu) * (y - yu) / (yl - yu);
                (at(yu + 13.0), at(yl - 13.0))
            }
            Some(lane) => {
                // down into the gap, along a lane right of the rows, up into the gap above
                let x_lane = cells_end + 14.0 + LANE * lane as f32;
                let y1 = top(b.upper.0) + CELL.y + PAD + GAP / 2.0 + 4.0 * lane as f32;
                let y2 = top(b.lower.0) - PAD - GAP / 2.0 - 4.0 * lane as f32;
                let points = [(xu, yu), (xu, y1), (x_lane, y1), (x_lane, y2), (xl, y2), (xl, yl)];
                for w in points.windows(2) {
                    painter.line_segment([egui::pos2(w[0].0, w[0].1), egui::pos2(w[1].0, w[1].1)], stroke);
                }
                painter.text(
                    egui::pos2((xu + x_lane) / 2.0, y1 - 3.0),
                    egui::Align2::CENTER_BOTTOM,
                    &name,
                    egui::FontId::proportional(13.0),
                    text_color,
                );
                (xu, xl)
            }
        };
        // the stubs found where this linker broke, at each end
        for ((row, _), x, y, up) in [(b.upper, x_upper_end, yu + 13.0, true), (b.lower, x_lower_end, yl - 13.0, false)] {
            let mut x = x + 12.0;
            for stub in bridge_end_stubs(a, link, b.link, l.rows[row]) {
                let color = stub_color(stub);
                triangle(&painter, egui::pos2(x, y), up, color);
                let r = painter.text(egui::pos2(x + 8.0, y), egui::Align2::LEFT_CENTER, &link.stubs[stub].name, small.clone(), color);
                x = r.right() + 12.0;
            }
        }
    }
    if !link.stubs.is_empty() {
        response.on_hover_text(stub_summary(a, link));
    }
}

/// One peptide with a bracket under it joining its two linked residues
fn loop_map(ui: &mut Ui, a: &Annotation, link: &model::Link) {
    let (i0, i1) = (link.sites[0].1.min(link.sites[1].1), link.sites[0].1.max(link.sites[1].1));
    let n = a.coverage[0].residues.len() as f32;
    let row_h = CELL.y + 2.0 * PAD;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(LABEL_W + CELL.x * n + 330.0, row_h + 24.0), egui::Sense::hover());
    let painter = ui.painter_at(rect.expand(4.0));
    let text_color = ui.visuals().text_color();
    let weak = ui.visuals().weak_text_color();
    let x0 = rect.left() + LABEL_W;
    let top = rect.top() + PAD;
    painter.text(egui::pos2(rect.left(), top + CELL.y / 2.0), egui::Align2::LEFT_CENTER, "loop-link", egui::FontId::proportional(14.0), text_color);
    draw_peptide(&painter, egui::pos2(x0, top), a, 0, &[i0, i1], text_color);
    let (xa, xb) = (x0 + CELL.x * (i0 as f32 + 0.5), x0 + CELL.x * (i1 as f32 + 0.5));
    let (y0, y1) = (top + CELL.y + 4.0, top + CELL.y + PAD + 8.0);
    let stroke = egui::Stroke::new(3.0_f32, weak);
    painter.line_segment([egui::pos2(xa, y0), egui::pos2(xa, y1)], stroke);
    painter.line_segment([egui::pos2(xa, y1), egui::pos2(xb, y1)], stroke);
    painter.line_segment([egui::pos2(xb, y1), egui::pos2(xb, y0)], stroke);
    painter.text(egui::pos2((xa + xb) / 2.0, y1 + 3.0), egui::Align2::CENTER_TOP, &link.linker, egui::FontId::proportional(12.0), text_color);
    painter.text(
        egui::pos2(x0 + CELL.x * n + 14.0, top + CELL.y / 2.0),
        egui::Align2::LEFT_CENTER,
        backbone_text(a, 0),
        egui::FontId::proportional(12.0),
        weak,
    );
}

fn legend(ui: &mut Ui, a: &Annotation) {
    let mut present: Vec<Series> = a
        .peaks
        .iter()
        .flat_map(|p| p.labels.iter().map(|l| l.series))
        .filter(|s| !matches!(s, Series::Stub(_)))
        .collect();
    present.sort_by_key(|s| s.name());
    present.dedup();
    ui.horizontal_wrapped(|ui| {
        for s in present {
            ui.label(RichText::new(format!("■ {}", s.name())).color(series_color(s)));
        }
        ui.label(RichText::new("■ unmatched").color(UNMATCHED));
        // the linker's stubs, dimmed when none of its ions was found
        if let Some(link) = &a.link {
            if !link.stubs.is_empty() {
                ui.separator();
                ui.label("stubs:");
                for (j, stub) in link.stubs.iter().enumerate() {
                    let found = a.stub_ions.iter().filter(|s| s.stubs.contains(&j)).count();
                    let text = format!("■ {} {} {:+.3}", stub.name, stub.description, stub.mass);
                    let label = if found > 0 {
                        RichText::new(text).color(stub_color(j))
                    } else {
                        RichText::new(text).color(stub_color(j).gamma_multiply(0.35))
                    };
                    ui.label(label).on_hover_text(format!("{found} ions with this stub"));
                }
            }
        }
    });
}

fn fragment_table(ui: &mut Ui, a: &Annotation) {
    egui::CollapsingHeader::new("Matched fragments")
        .id_salt("viewer_fragments")
        .show(ui, |ui| {
            egui::ScrollArea::vertical().id_salt("viewer_fragment_rows").max_height(240.0).show(ui, |ui| {
                egui::Grid::new("viewer_fragment_grid").striped(true).num_columns(5).show(ui, |ui| {
                    for h in ["ion", "theoretical m/z", "observed m/z", "error (ppm)", "intensity"] {
                        ui.label(crate::theme::bold(h));
                    }
                    ui.end_row();
                    for p in &a.peaks {
                        for l in &p.labels {
                            ui.label(RichText::new(label_text(a, l)).color(series_color(l.series)));
                            ui.label(format!("{:.4}", l.theoretical_mz));
                            ui.label(format!("{:.4}", p.mz));
                            ui.label(format!("{:+.1}", (p.mz - l.theoretical_mz) / l.theoretical_mz * 1e6));
                            ui.label(format!("{:.3e}", p.intensity));
                            ui.end_row();
                        }
                    }
                });
            });
        });
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::Path;

    /// Lay out the coverage map of a cleaved DSSO crosslink headlessly and
    /// look at what it draws
    #[test]
    fn crosslink_map_draws_bridge_and_stub_glyphs() {
        let (_, a) = model::test::cleaved_dsso();
        let ctx = egui::Context::default();
        let output = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| coverage_view(ui, &a));
        });
        let (mut triangles, mut segments, mut texts) = (Vec::new(), Vec::new(), Vec::new());
        for clipped in &output.shapes {
            match &clipped.shape {
                egui::Shape::Path(path) if path.closed && path.points.len() == 3 => triangles.push(path.fill),
                egui::Shape::LineSegment { stroke, .. } => segments.push(stroke.color),
                egui::Shape::Text(text) => texts.push((text.galley.text().to_string(), text.fallback_color)),
                _ => {}
            }
        }
        // P+A and P+S on both peptides: four triangles in the A and S colours
        assert_eq!(triangles.len(), 4, "{triangles:?}");
        assert_eq!(triangles.iter().filter(|c| **c == stub_color(0)).count(), 2);
        assert_eq!(triangles.iter().filter(|c| **c == stub_color(2)).count(), 2);
        // the S stub hooks (y4-y6+S of peptide 2: two segments each)
        assert_eq!(segments.iter().filter(|c| **c == stub_color(2)).count(), 6);
        // the bridge, labelled with the linker; both peptides and the stub names
        assert!(texts.iter().any(|(t, _)| t == "DSSO"), "{texts:?}");
        assert!(texts.iter().any(|(t, _)| t == "peptide 1") && texts.iter().any(|(t, _)| t == "peptide 2"));
        // the stub names next to the triangles, in the stub colour
        assert_eq!(texts.iter().filter(|(t, c)| t == "A" && *c == stub_color(0)).count(), 2);
        assert_eq!(texts.iter().filter(|(t, c)| t == "S" && *c == stub_color(2)).count(), 2);
        assert!(segments.len() >= 7, "bridge and hooks");
    }

    /// The map of a cleaved DSSO `a1b1c` species: both links drawn and
    /// labelled, the stubs where each linker broke, the triplet peaks of the
    /// middle peptide and its stub-carrying y ions as hooks
    #[test]
    fn species_map_draws_every_link_and_its_stubs() {
        let (_, a) = model::test::cleaved_multi();
        let ctx = egui::Context::default();
        let output = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| coverage_view(ui, &a));
        });
        let (mut triangles, mut segments, mut texts) = (Vec::new(), Vec::new(), Vec::new());
        for clipped in &output.shapes {
            match &clipped.shape {
                egui::Shape::Path(path) if path.closed && path.points.len() == 3 => triangles.push(path.fill),
                egui::Shape::LineSegment { stroke, .. } => segments.push(stroke.color),
                egui::Shape::Text(text) => texts.push(text.galley.text().to_string()),
                _ => {}
            }
        }
        for t in ["peptide 1", "peptide 2", "peptide 3", "XL1", "XL2", "A+S", "S+S"] {
            assert!(texts.iter().any(|x| x == t), "{t} in {texts:?}");
        }
        // XL1: peptide 1 + A, + S and {2,3} + A; XL2: {1,2} + S and peptide 3 + A;
        // after peptide 2: A+S and S+S
        assert_eq!(triangles.len(), 9, "{triangles:?}");
        // y8+S (K19) and y17+A+S (K10 and K19): S hooks at two residues, an A hook at one
        assert_eq!(segments.iter().filter(|c| **c == stub_color(2)).count(), 4);
        assert_eq!(segments.iter().filter(|c| **c == stub_color(0)).count(), 2);
    }

    /// Timings on a real results folder (`SAGE_GUI_RESULTS=<folder>`): the
    /// mobiligrams written after a search (into a copy of the folder), and
    /// the viewer from "Load" to the first spectrum shown. `SAGE_GUI_FIGURES`
    /// saves a render of a multivalent species.
    #[test]
    #[ignore]
    fn real_folder_timing() {
        use egui_kittest::Harness;
        use std::time::Instant;
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let copy = std::env::temp_dir().join(format!("sage-gui-timing-{}", std::process::id()));
        std::fs::create_dir_all(&copy).unwrap();
        for name in ["results.json", "results.sage.tsv", "crosslinks.sage.tsv"] {
            let _ = std::fs::copy(Path::new(&folder).join(name), copy.join(name));
        }
        let t = Instant::now();
        let written = mobiligram::write_figures(&copy).unwrap();
        println!("mobiligrams: {} files in {:.2?}", written.len(), t.elapsed());
        std::fs::remove_dir_all(&copy).unwrap();

        let mut viewer = Viewer::default();
        let t = Instant::now();
        viewer.open(&folder);
        viewer.wait_until_open();
        println!("results read in {:.2?}: {} shown", t.elapsed(), viewer.shown.len());
        let results = viewer.results.clone().unwrap();
        // the best multivalent species, picked as a user would
        let multi = viewer.shown.iter().copied().find(|&i| results.matches[i].kind == Kind::Multi);
        viewer.selected = multi.or(viewer.selected);
        let mut installed = false;
        let mut harness = Harness::builder().with_size(egui::vec2(1500.0, 1500.0)).build_ui_state(
            move |ui, v: &mut Viewer| {
                if !installed {
                    crate::theme::install(ui.ctx(), false);
                    installed = true;
                    return;
                }
                egui::ScrollArea::vertical().show(ui, |ui| v.ui(ui, ""));
            },
            viewer,
        );
        let t = Instant::now();
        let key = results.matches[harness.state().selected.unwrap()].key();
        while !harness.state().spectra.contains_key(&key) && t.elapsed().as_secs() < 60 {
            harness.step();
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        println!("first spectrum ({} scan {}) shown after {:.2?}", key.0, key.1, t.elapsed());
        // another spectrum of the same run, now open
        let other = harness
            .state()
            .shown
            .iter()
            .map(|&i| results.matches[i].key())
            .find(|k| k.0 == key.0 && !harness.state().spectra.contains_key(k));
        if let Some(other) = other {
            let i = harness.state().shown.iter().copied().find(|&i| results.matches[i].key() == other).unwrap();
            harness.state_mut().selected = Some(i);
            let t = Instant::now();
            while !harness.state().spectra.contains_key(&other) && t.elapsed().as_secs() < 60 {
                harness.step();
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
            println!("next spectrum of the run shown after {:.2?}", t.elapsed());
            harness.state_mut().selected = multi;
        }
        harness.run_steps(4);
        let state = harness.state();
        if let Some(Ok(a)) = &state.annotation {
            println!("{}: {} stub ions: {:?}", a.proforma, a.stub_ions.len(), a.stub_ions.iter().map(|s| &s.text).collect::<Vec<_>>());
        }
        if let Ok(out) = std::env::var("SAGE_GUI_FIGURES") {
            std::fs::create_dir_all(&out).unwrap();
            harness.render().unwrap().save(Path::new(&out).join("viewer_multivalent.png")).unwrap();
        }
    }

    /// A results folder with three crosslinks, two with long titles
    fn long_titles() -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("sage-gui-divider-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("results.json"), r#"{"mzml_paths": ["a.d"]}"#).unwrap();
        let mut tsv = String::from(
            "filename\tscannr\trank\tlabel\tclass\tlink_type\tpeptide_a\tpeptide_b\tlink_a\tlink_b\t\
             proteins_a\tproteins_b\tcharge\thyperscore\tspectrum_q\tcrosslink_q\trt\n",
        );
        for (i, (a, b)) in [
            ("KVPSPTGEYVKEPEWVLETDGVNLSEVMTVPGIDPTR", "IQYAKDILQKEFLPHITQLEGFESR"),
            ("TVISGDPNLELDQVGVPKSIAKTLTYPEVVTPYNIDR", "SAGNIVQVFKKAAPSPISHVAEIR"),
            ("GLEKVSTAAMQR", "TDYHFPNEWKLIR"),
        ]
        .iter()
        .enumerate()
        {
            tsv.push_str(&format!(
                "a.d\t{i}\t1\t1\tTT\tcrosslink\t{a}\t{b}\t1\t5\tsp|P1|A\tsp|P2|B\t4\t{}\t0.001\t0.001\t30\n",
                100 - i
            ));
        }
        std::fs::write(dir.join("crosslinks.sage.tsv"), tsv).unwrap();
        dir
    }

    /// Drag the divider between the match list and the spectrum with the
    /// pointer, as a user would. `SAGE_GUI_FIGURES` saves before / after images.
    #[test]
    fn divider_drags_to_resize_the_list() {
        use egui_kittest::{kittest::Queryable, Harness};
        let dir = long_titles();
        let mut viewer = Viewer::default();
        viewer.open(&dir.display().to_string());
        viewer.wait_until_open();
        let mut installed = false;
        // 60 frames a second, so two clicks fall within egui's double-click time
        let mut harness = Harness::builder().with_size(egui::vec2(1400.0, 900.0)).with_step_dt(1.0 / 60.0).build_ui_state(
            move |ui, v: &mut Viewer| {
                // fonts take effect on the frame after they are set
                if !installed {
                    crate::theme::install(ui.ctx(), false);
                    installed = true;
                    return;
                }
                v.ui(ui, "");
            },
            viewer,
        );
        harness.run_steps(6);
        let save = |harness: &mut Harness<'_, Viewer>, name: &str| {
            if let Ok(out) = std::env::var("SAGE_GUI_FIGURES") {
                std::fs::create_dir_all(&out).unwrap();
                harness.render().unwrap().save(std::path::Path::new(&out).join(name)).unwrap();
            }
        };
        save(&mut harness, "viewer_divider_before.png");
        assert_eq!(harness.state().list_width, LIST_WIDTH);

        let drag = |harness: &mut Harness<'_, Viewer>, dx: f32| {
            let r = harness.get_by_label("Resize the match list").bounding_box().unwrap();
            let start = egui::pos2(((r.x0 + r.x1) / 2.0) as f32, ((r.y0 + r.y1) / 2.0) as f32);
            let button = |pressed| egui::Event::PointerButton {
                pos: if pressed { start } else { start + egui::vec2(dx, 0.0) },
                button: egui::PointerButton::Primary,
                pressed,
                modifiers: Default::default(),
            };
            harness.input_mut().events.push(egui::Event::PointerMoved(start));
            harness.run_steps(6);
            harness.input_mut().events.push(button(true));
            harness.run_steps(6);
            for step in 1..=5 {
                harness.input_mut().events.push(egui::Event::PointerMoved(start + egui::vec2(dx * step as f32 / 5.0, 0.0)));
                harness.run_steps(6);
            }
            harness.input_mut().events.push(button(false));
            harness.run_steps(6);
        };
        // wider by the distance dragged
        drag(&mut harness, 250.0);
        let wide = harness.state().list_width;
        assert!((wide - (LIST_WIDTH + 250.0)).abs() < 2.0, "list width {wide}");
        save(&mut harness, "viewer_divider_after.png");
        // never so wide that the spectrum loses its room, never narrower than the minimum
        drag(&mut harness, 2000.0);
        assert!(harness.state().list_width <= 1400.0 - MIN_SPECTRUM_WIDTH, "{}", harness.state().list_width);
        drag(&mut harness, -3000.0);
        assert_eq!(harness.state().list_width, MIN_LIST_WIDTH);
        // double-click: back to the default width
        let r = harness.get_by_label("Resize the match list").bounding_box().unwrap();
        let at = egui::pos2(((r.x0 + r.x1) / 2.0) as f32, ((r.y0 + r.y1) / 2.0) as f32);
        harness.input_mut().events.push(egui::Event::PointerMoved(at));
        harness.run_steps(6);
        for _ in 0..2 {
            for pressed in [true, false] {
                harness.input_mut().events.push(egui::Event::PointerButton {
                    pos: at,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                });
                harness.step();
            }
        }
        harness.run_steps(6);
        assert_eq!(harness.state().list_width, LIST_WIDTH);
        // rows stay one line: every match label is as tall as the shortest
        let heights: Vec<f64> = harness
            .get_all_by_role(egui::accesskit::Role::Button)
            .filter(|n| n.label().is_some_and(|l| l.starts_with("XL ")))
            .filter_map(|n| n.bounding_box())
            .map(|b| b.y1 - b.y0)
            .collect();
        assert_eq!(heights.len(), 3, "{heights:?}");
        assert!(heights.iter().all(|h| (h - heights[0]).abs() < 0.5), "{heights:?}");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
