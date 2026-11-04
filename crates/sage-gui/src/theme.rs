//! Look of the window: a clean light theme (white panels, thin grey borders,
//! blue accents) and a matching dark theme, plus the few widgets that carry it
//! (section headings with a rule, outlined coloured buttons, a light/dark switch,
//! number fields with + / − buttons, bold titles).

use eframe::egui::{
    self, emath::Numeric, pos2, vec2, Button, Color32, CornerRadius, DragValue, FontFamily, FontId, Rect,
    Response, RichText, Sense, Stroke, StrokeKind, TextStyle, Theme, Ui, Visuals,
};
use std::ops::RangeInclusive;

/// Font family of titles: a bold system font, or the regular one when none is found
pub const BOLD: &str = "bold";

fn bold_family() -> FontFamily {
    FontFamily::Name(BOLD.into())
}

/// Title text (tabs, sections, column headers) in the bold font
pub fn bold(text: impl Into<String>) -> RichText {
    RichText::new(text).family(bold_family())
}

/// Register the bold font family: Segoe UI Bold or Arial Bold on Windows
/// (DejaVu Sans Bold, Arial Bold elsewhere), then the regular fonts for
/// anything it lacks
fn install_fonts(ctx: &egui::Context) {
    let mut fonts = egui::FontDefinitions::default();
    let candidates = [
        r"C:\Windows\Fonts\segoeuib.ttf",
        r"C:\Windows\Fonts\arialbd.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ];
    let mut family = Vec::new();
    if let Some(bytes) = candidates.iter().find_map(|p| std::fs::read(p).ok()) {
        fonts.font_data.insert(BOLD.into(), std::sync::Arc::new(egui::FontData::from_owned(bytes)));
        family.push(BOLD.to_string());
    }
    family.extend(fonts.families.get(&FontFamily::Proportional).cloned().unwrap_or_default());
    fonts.families.insert(bold_family(), family);
    ctx.set_fonts(fonts);
}

/// A number field: the value is typed in (dragging does nothing), and two
/// stacked buttons to its right step it by `step`, + above −, within `range`.
/// `style` adds prefix, suffix or decimals to the field.
pub fn number<T: Numeric>(
    ui: &mut Ui,
    value: &mut T,
    range: RangeInclusive<T>,
    step: f64,
    style: impl FnOnce(DragValue) -> DragValue,
) -> Response {
    let (lo, hi) = (range.start().to_f64(), range.end().to_f64());
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 1.0;
        // 2 px taller than other fields, so each stacked button has room for its symbol
        ui.spacing_mut().interact_size.y += 2.0;
        // without a drag speed egui would show every decimal: up to 6 (or the
        // field's own maximum), trailing zeros trimmed
        let field = DragValue::new(value)
            .range(range)
            .speed(0.0)
            .custom_formatter(|v, decimals| plain_number(v, (*decimals.end()).min(6)));
        let mut field = ui.add(style(field)).on_hover_cursor(egui::CursorIcon::Text);
        let (area, _) = ui.allocate_exact_size(vec2(16.0, field.rect.height()), Sense::hover());
        let mut buttons = Vec::with_capacity(2);
        for (up, rect) in [
            (true, Rect::from_min_max(area.min, pos2(area.max.x, area.center().y))),
            (false, Rect::from_min_max(pos2(area.min.x, area.center().y), area.max)),
        ] {
            let rect = rect.shrink2(vec2(0.0, 0.5));
            let button = ui.interact(rect, field.id.with(up), Sense::click());
            let look = ui.style().interact(&button);
            let painter = ui.painter();
            painter.rect(rect, 2.0, look.weak_bg_fill, look.bg_stroke, StrokeKind::Inside);
            // the symbol fits its half-height box with 2 px clear of the border,
            // on whole pixels so it stays crisp
            // pixel centre at or below the middle (rounding up would push the + into the top border)
            let ppp = ui.pixels_per_point();
            let mid = rect.center();
            let c = pos2(
                ((mid.x * ppp - 0.5).round() + 0.5) / ppp,
                ((mid.y * ppp - 0.5).ceil() + 0.5) / ppp,
            );
            let arm = 2.5f32;
            let rise = (rect.height() / 2.0 - 3.0).clamp(1.5, arm);
            let stroke = Stroke::new(1.2f32, look.fg_stroke.color);
            painter.line_segment([c - vec2(arm, 0.0), c + vec2(arm, 0.0)], stroke);
            if up {
                painter.line_segment([c - vec2(0.0, rise), c + vec2(0.0, rise)], stroke);
            }
            if button.clicked() {
                *value = T::from_f64(stepped(value.to_f64(), step, up, lo, hi));
                field.mark_changed();
            }
            buttons.push(button.on_hover_text(if up { format!("+{step}") } else { format!("−{step}") }));
        }
        buttons.into_iter().fold(field, |acc, b| acc | b)
    })
    .inner
}

/// `v` with at most `decimals` decimals, trailing zeros trimmed: 5000, 0.01, 15.9949
fn plain_number(v: f64, decimals: usize) -> String {
    let text = format!("{v:.decimals$}");
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        text
    }
}

/// `value` one `step` up or down, within `lo..=hi`. Rounded to the decimals of
/// the step or of the value, whichever has more (up to 6), so 0.03 + 0.01 is
/// 0.04, not 0.04000000000000001, and 15.9949 + 1 stays 16.9949.
fn stepped(value: f64, step: f64, up: bool, lo: f64, hi: f64) -> f64 {
    let decimals_of = |x: f64| {
        (0..=6)
            .find(|&d| {
                let s = 10f64.powi(d);
                ((x * s).round() / s - x).abs() < 1e-9
            })
            .unwrap_or(6)
    };
    let scale = 10f64.powi(decimals_of(step).max(decimals_of(value)));
    let next = value + if up { step } else { -step };
    ((next * scale).round() / scale).clamp(lo, hi)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn steps_round_and_stay_in_range() {
        assert_eq!(stepped(3.0, 1.0, true, 0.0, 10.0), 4.0);
        assert_eq!(stepped(0.0, 1.0, false, 0.0, 10.0), 0.0);
        assert_eq!(stepped(10.0, 1.0, true, 0.0, 10.0), 10.0);
        assert_eq!(stepped(0.03, 0.01, true, 0.0, 1.0), 0.04);
        assert_eq!(stepped(0.7, 0.05, false, 0.0, 1.0), 0.65);
        // a typed value keeps its own decimals
        assert_eq!(stepped(20.37, 1.0, true, -1e4, 1e4), 21.37);
        assert_eq!(stepped(15.9949, 1.0, true, -1e4, 1e4), 16.9949);
        assert_eq!(stepped(-17.026549, 1.0, false, -1e4, 1e4), -18.026549);
        assert_eq!(stepped(4096.0, 1024.0, true, 1024.0, 1_048_576.0), 5120.0);
    }

    #[test]
    fn numbers_show_without_trailing_zeros() {
        assert_eq!(plain_number(5000.0, 6), "5000");
        assert_eq!(plain_number(0.01, 6), "0.01");
        assert_eq!(plain_number(15.994915, 6), "15.994915");
        assert_eq!(plain_number(0.7, 6), "0.7");
        assert_eq!(plain_number(-17.026549, 6), "-17.026549");
        assert_eq!(plain_number(0.0125, 3), "0.013");
        assert_eq!(plain_number(2.0, 0), "2");
    }
}

const fn rgb(r: u8, g: u8, b: u8) -> Color32 {
    Color32::from_rgb(r, g, b)
}

/// Colours of one theme
struct Palette {
    panel: Color32,
    field: Color32,
    faint: Color32,
    text: Color32,
    border: Color32,
    rule: Color32,
    hover: Color32,
    pressed: Color32,
    accent: Color32,
    selection: Color32,
}

const LIGHT: Palette = Palette {
    panel: Color32::WHITE,
    field: Color32::WHITE,
    faint: rgb(246, 248, 251),
    text: rgb(45, 48, 54),
    border: rgb(200, 205, 212),
    rule: rgb(225, 228, 233),
    hover: rgb(244, 247, 252),
    pressed: rgb(232, 240, 252),
    accent: rgb(47, 111, 207),
    selection: rgb(214, 228, 248),
};

const DARK: Palette = Palette {
    panel: rgb(30, 31, 34),
    field: rgb(24, 25, 27),
    faint: rgb(38, 40, 44),
    text: rgb(220, 222, 226),
    border: rgb(78, 82, 90),
    rule: rgb(56, 59, 65),
    hover: rgb(40, 44, 52),
    pressed: rgb(40, 56, 84),
    accent: rgb(110, 160, 240),
    selection: rgb(40, 62, 100),
};

/// Thin border / rule of colour `c`
fn line(c: Color32) -> Stroke {
    Stroke::new(1.0_f32, c)
}

fn palette(dark: bool) -> &'static Palette {
    if dark {
        &DARK
    } else {
        &LIGHT
    }
}

fn visuals(dark: bool) -> Visuals {
    let p = palette(dark);
    let mut v = if dark { Visuals::dark() } else { Visuals::light() };
    v.panel_fill = p.panel;
    v.window_fill = p.panel;
    v.extreme_bg_color = p.field;
    v.faint_bg_color = p.faint;
    v.window_stroke = line(p.border);
    v.window_corner_radius = CornerRadius::same(4);
    v.hyperlink_color = p.accent;
    v.selection.bg_fill = p.selection;
    v.selection.stroke = line(p.accent);

    let w = &mut v.widgets;
    w.noninteractive.bg_fill = p.panel;
    w.noninteractive.weak_bg_fill = p.panel;
    w.noninteractive.bg_stroke = line(p.rule);
    w.noninteractive.fg_stroke = line(p.text);
    for (state, fill, stroke) in [
        (&mut w.inactive, p.field, p.border),
        (&mut w.hovered, p.hover, p.accent),
        (&mut w.active, p.pressed, p.accent),
        (&mut w.open, p.pressed, p.accent),
    ] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.bg_stroke = line(stroke);
        state.fg_stroke = line(p.text);
        state.corner_radius = CornerRadius::same(2);
        state.expansion = 0.0;
    }
    w.noninteractive.corner_radius = CornerRadius::same(2);
    v
}

/// Install both themes and pick one
pub fn install(ctx: &egui::Context, dark: bool) {
    install_fonts(ctx);
    for theme in [Theme::Light, Theme::Dark] {
        ctx.style_mut_of(theme, |s| {
            s.visuals = visuals(theme == Theme::Dark);
            s.text_styles = [
                (TextStyle::Heading, FontId::new(18.0, bold_family())),
                (TextStyle::Body, FontId::proportional(14.0)),
                (TextStyle::Button, FontId::proportional(14.0)),
                (TextStyle::Monospace, FontId::monospace(13.0)),
                (TextStyle::Small, FontId::proportional(11.0)),
            ]
            .into();
            s.spacing.item_spacing = vec2(8.0, 6.0);
            s.spacing.button_padding = vec2(10.0, 3.0);
            s.spacing.interact_size.y = 22.0;
        });
    }
    set_dark(ctx, dark);
}

pub fn set_dark(ctx: &egui::Context, dark: bool) {
    ctx.set_theme(if dark { Theme::Dark } else { Theme::Light });
}

pub fn accent(ui: &Ui) -> Color32 {
    palette(ui.visuals().dark_mode).accent
}

/// Grey of field labels, a step lighter than values
pub fn muted(ui: &Ui) -> Color32 {
    if ui.visuals().dark_mode {
        rgb(160, 165, 175)
    } else {
        rgb(110, 115, 125)
    }
}

/// Blue section title followed by a thin rule to the right edge
pub fn section(ui: &mut Ui, text: &str) {
    ui.add_space(10.0);
    ui.horizontal(|ui| {
        let title = ui.label(bold(text).size(17.0).color(accent(ui)));
        let right = ui.available_rect_before_wrap().right();
        let left = title.rect.right() + 6.0;
        if right > left {
            let rule = ui.visuals().widgets.noninteractive.bg_stroke;
            ui.painter().hline(left..=right, title.rect.center().y, rule);
        }
    });
    ui.add_space(2.0);
}

#[derive(Clone, Copy)]
pub enum Tone {
    /// Primary action (Run)
    Blue,
    /// Stopping / cancelling
    Tan,
}

/// Outlined button tinted by `tone`
pub fn tone_button(ui: &Ui, text: &str, tone: Tone) -> Button<'static> {
    let dark = ui.visuals().dark_mode;
    let (stroke, fill, label) = match (tone, dark) {
        (Tone::Blue, false) => (LIGHT.accent, rgb(232, 240, 252), LIGHT.accent),
        (Tone::Blue, true) => (DARK.accent, rgb(30, 45, 70), DARK.accent),
        (Tone::Tan, false) => (rgb(190, 150, 100), rgb(250, 243, 232), rgb(160, 110, 60)),
        (Tone::Tan, true) => (rgb(200, 160, 110), rgb(60, 48, 35), rgb(220, 180, 130)),
    };
    Button::new(RichText::new(text).color(label))
        .fill(fill)
        .stroke(line(stroke))
        .min_size(vec2(84.0, 26.0))
}

/// ☀ [switch] 🌙 light/dark toggle; returns true when flipped
pub fn theme_switch(ui: &mut Ui, dark: &mut bool) -> bool {
    let muted = muted(ui);
    // right-to-left layout: the moon is placed first
    ui.label(RichText::new("🌙").color(muted));
    let size = vec2(34.0, 18.0);
    let (rect, response): (_, Response) = ui.allocate_exact_size(size, Sense::click());
    let clicked = response.clicked();
    if clicked {
        *dark = !*dark;
    }
    let t = ui.ctx().animate_bool_responsive(response.id, *dark);
    let v = ui.style().interact_selectable(&response, *dark);
    let radius = 0.5 * rect.height();
    ui.painter().rect(rect, radius, v.bg_fill, v.bg_stroke, StrokeKind::Inside);
    let x = egui::lerp((rect.left() + radius)..=(rect.right() - radius), t);
    ui.painter().circle(pos2(x, rect.center().y), 0.7 * radius, accent(ui), Stroke::NONE);
    response.on_hover_text(if *dark { "Switch to the light theme" } else { "Switch to the dark theme" });
    ui.label(RichText::new("☀").color(muted));
    clicked
}
