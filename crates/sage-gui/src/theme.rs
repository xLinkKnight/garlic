//! Look of the window: a clean light theme (white panels, thin grey borders,
//! blue accents) and a matching dark theme, plus the few widgets that carry it
//! (section headings with a rule, outlined coloured buttons, a light/dark switch).

use eframe::egui::{
    self, pos2, vec2, Button, Color32, CornerRadius, FontId, Response, RichText, Sense, Stroke,
    StrokeKind, TextStyle, Theme, Ui, Visuals,
};

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
    for theme in [Theme::Light, Theme::Dark] {
        ctx.style_mut_of(theme, |s| {
            s.visuals = visuals(theme == Theme::Dark);
            s.text_styles = [
                (TextStyle::Heading, FontId::proportional(18.0)),
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
        let title = ui.label(RichText::new(text).size(17.0).color(accent(ui)));
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
