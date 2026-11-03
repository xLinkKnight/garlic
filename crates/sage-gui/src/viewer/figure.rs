//! Publication figures of an annotated spectrum: an SVG (vector, editable in
//! Illustrator / Inkscape) drawn from the annotation, rendered to PNG or TIFF
//! at print resolution. White or transparent background.
//!
//! The figure is 7 inches wide (a two-column journal figure): the SVG's user
//! units are 150 per inch, so 1050 wide.

use super::model::{Annotation, Kind, Link, PeakView, Series, StubIonKind};
use super::{label_text, series_color, stub_color, UNMATCHED};
use eframe::egui::Color32;
use std::fmt::Write;
use std::path::Path;

/// Figure width in inches, and SVG user units per inch
const WIDTH_IN: f64 = 7.0;
const UNITS_PER_IN: f64 = 150.0;
const W: f64 = WIDTH_IN * UNITS_PER_IN;
const FONT: &str = "Arial, Helvetica, sans-serif";
const MONO: &str = "Consolas, 'Courier New', monospace";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Svg,
    Png,
    Tiff,
}

impl Format {
    pub const ALL: [Format; 3] = [Format::Svg, Format::Png, Format::Tiff];

    pub fn label(self) -> &'static str {
        match self {
            Format::Svg => "SVG (vector)",
            Format::Png => "PNG",
            Format::Tiff => "TIFF",
        }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Format::Svg => "svg",
            Format::Png => "png",
            Format::Tiff => "tif",
        }
    }
}

#[derive(Clone, Debug)]
pub struct FigureOptions {
    /// White background; otherwise transparent
    pub background: bool,
    /// Title line (peptides, charge, scan) above the spectrum
    pub title: Option<String>,
    /// Label every matched peak, not only those above 2% of the base peak
    pub label_all: bool,
    /// m/z range shown (the viewer's zoom); whole spectrum when `None`
    pub mz_range: Option<(f64, f64)>,
    /// Peptide map (with the crosslink bridge) under the spectrum
    pub coverage: bool,
    /// Raster resolution (PNG, TIFF)
    pub dpi: u32,
}

fn hex(c: Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// Label text as SVG: charges written as superscript digits become a raised
/// `<tspan>` (`y9²+` → y9 with a raised `2+`)
fn label_svg(text: &str) -> String {
    const SUP: &str = "⁰¹²³⁴⁵⁶⁷⁸⁹";
    let mut out = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if let Some(d) = SUP.chars().position(|s| s == c) {
            let mut raised = d.to_string();
            while let Some(&n) = chars.peek() {
                match SUP.chars().position(|s| s == n) {
                    Some(d) => raised.push_str(&d.to_string()),
                    None if n == '+' => raised.push('+'),
                    None => break,
                }
                chars.next();
                if raised.ends_with('+') {
                    break;
                }
            }
            let _ = write!(out, r#"<tspan baseline-shift="super" font-size="75%">{}</tspan>"#, esc(&raised));
        } else {
            out.push_str(&esc(&c.to_string()));
        }
    }
    out
}

/// A round step giving about `n` ticks over `span`
fn tick_step(span: f64, n: f64) -> f64 {
    let raw = (span / n).max(f64::MIN_POSITIVE);
    let magnitude = 10f64.powf(raw.log10().floor());
    [1.0, 2.0, 2.5, 5.0, 10.0]
        .into_iter()
        .map(|f| f * magnitude)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * magnitude)
}

/// The figure as SVG
pub fn svg(a: &Annotation, opts: &FigureOptions) -> String {
    let (lo, hi) = opts.mz_range.unwrap_or_else(|| {
        // the data, with a little room, on round hundreds
        let (min, max) = a.peaks.iter().fold((f64::MAX, 0.0f64), |(lo, hi), p| (lo.min(p.mz), hi.max(p.mz)));
        if min > max {
            return (0.0, 100.0);
        }
        (((min * 0.95) / 100.0).floor() * 100.0, (max * 1.03).max(min + 100.0))
    });
    let peaks: Vec<_> = a.peaks.iter().filter(|p| p.mz >= lo && p.mz <= hi).collect();
    let base = peaks.iter().map(|p| p.intensity).fold(0.0f32, f32::max).max(f32::MIN_POSITIVE);

    let mut body = String::new();
    let mut y = 16.0;
    if let Some(title) = &opts.title {
        let _ = write!(body, r#"<text x="20" y="{:.1}" font-family="{FONT}" font-size="15" font-weight="bold">{}</text>"#, y + 12.0, esc(title));
        y += 26.0;
    }

    // ---- the spectrum
    let (x0, x1) = (78.0, W - 24.0);
    let (top, bottom) = (y + 46.0, y + 360.0);
    let px = |mz: f64| x0 + (mz - lo) / (hi - lo) * (x1 - x0);
    let py = |rel: f64| bottom - rel * (bottom - top);
    let axis = "#333333";
    let _ = write!(body, r#"<g font-family="{FONT}" font-size="12" fill="{axis}">"#);
    // axes and ticks
    let _ = write!(body, r#"<path d="M{x0:.1},{top:.1} V{bottom:.1} H{x1:.1}" fill="none" stroke="{axis}" stroke-width="1.2"/>"#);
    let step = tick_step(hi - lo, 8.0);
    let mut t = (lo / step).ceil() * step;
    while t <= hi + 1e-9 {
        let x = px(t);
        let _ = write!(body, r#"<line x1="{x:.1}" y1="{bottom:.1}" x2="{x:.1}" y2="{:.1}" stroke="{axis}" stroke-width="1"/>"#, bottom + 5.0);
        let _ = write!(body, r#"<text x="{x:.1}" y="{:.1}" text-anchor="middle">{}</text>"#, bottom + 19.0, t.round() as i64);
        t += step;
    }
    for pct in [0, 25, 50, 75, 100] {
        let yy = py(pct as f64 / 100.0);
        let _ = write!(body, r#"<line x1="{:.1}" y1="{yy:.1}" x2="{x0:.1}" y2="{yy:.1}" stroke="{axis}" stroke-width="1"/>"#, x0 - 5.0);
        let _ = write!(body, r#"<text x="{:.1}" y="{:.1}" text-anchor="end">{pct}</text>"#, x0 - 8.0, yy + 4.0);
    }
    let _ = write!(body, r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-style="italic">m/z</text>"#, (x0 + x1) / 2.0, bottom + 38.0);
    let _ = write!(
        body,
        r#"<text transform="translate(26,{:.1}) rotate(-90)" text-anchor="middle">Relative intensity (%)</text>"#,
        (top + bottom) / 2.0
    );
    body.push_str("</g>");

    // peaks: unmatched first, matched on top
    for matched in [false, true] {
        for p in peaks.iter().filter(|p| p.labels.is_empty() != matched) {
            let (color, width) = match p.labels.first() {
                Some(l) => (hex(series_color(l.series)), 1.6),
                None => (hex(UNMATCHED), 1.0),
            };
            let x = px(p.mz);
            let _ = write!(
                body,
                r#"<line x1="{x:.2}" y1="{bottom:.1}" x2="{x:.2}" y2="{:.2}" stroke="{color}" stroke-width="{width}"/>"#,
                py((p.intensity / base) as f64)
            );
        }
    }

    // labels, strongest first, moved up (with a leader) to avoid overlaps
    let threshold = if opts.label_all { 0.0 } else { 0.02 * base };
    // internal fragments only with "label all": their labels (`my:b7:10`) crowd a figure
    let shown = |p: &&&PeakView| match p.labels.first() {
        Some(l) => opts.label_all || !matches!(l.series, Series::Internal | Series::Other),
        None => false,
    };
    let mut order: Vec<_> = peaks.iter().filter(|p| p.intensity >= threshold).filter(shown).collect();
    order.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
    let mut placed: Vec<(f64, f64, f64, f64)> = Vec::new();
    let _ = write!(body, r#"<g font-family="{FONT}" font-size="11" text-anchor="middle">"#);
    for p in order {
        let l = &p.labels[0];
        let mut text = label_text(a, l);
        if p.labels.len() > 1 {
            text.push('*');
        }
        let width = text.chars().count() as f64 * 6.2 + 4.0;
        let x = px(p.mz);
        let peak_top = py((p.intensity / base) as f64);
        let mut baseline = peak_top - 4.0;
        let overlaps = |b: f64, placed: &[(f64, f64, f64, f64)]| {
            placed
                .iter()
                .any(|&(l, r, t, btm)| x - width / 2.0 < r && x + width / 2.0 > l && b - 11.0 < btm && b > t)
        };
        let mut tries = 0;
        while overlaps(baseline, &placed) && tries < 8 {
            baseline -= 12.0;
            tries += 1;
        }
        if tries == 8 || baseline - 11.0 < top - 40.0 {
            continue;
        }
        let color = hex(series_color(l.series));
        if tries > 0 {
            let _ = write!(
                body,
                r#"<line x1="{x:.1}" y1="{:.1}" x2="{x:.1}" y2="{:.1}" stroke="{color}" stroke-width="0.6" stroke-dasharray="2,2"/>"#,
                peak_top - 1.0,
                baseline + 2.0
            );
        }
        let _ = write!(body, r#"<text x="{x:.1}" y="{baseline:.1}" fill="{color}">{}</text>"#, label_svg(&text));
        placed.push((x - width / 2.0, x + width / 2.0, baseline - 11.0, baseline + 2.0));
    }
    body.push_str("</g>");
    y = bottom + 56.0;

    // ---- peptide map
    if opts.coverage {
        y = coverage_svg(&mut body, a, y);
    }

    // ---- legend
    y = legend_svg(&mut body, a, y + 6.0);

    let height = y + 12.0;
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_IN}in" height="{:.3}in" viewBox="0 0 {W} {height:.1}">"#,
        height / UNITS_PER_IN
    );
    if opts.background {
        let _ = write!(svg, r#"<rect width="100%" height="100%" fill="white"/>"#);
    }
    svg.push_str(&body);
    svg.push_str("</svg>\n");
    svg
}

const CELL_W: f64 = 21.0;
const CELL_H: f64 = 26.0;

/// Hook at a residue corner, as in the viewer: N-terminal ions over its
/// right, C-terminal under its left; stubs nest outward
fn hook_svg(out: &mut String, left: f64, top: f64, n_term: bool, level: usize, color: Color32) {
    let o = 3.0 * level as f64;
    let (cx, cy, dy, dx) = if n_term {
        (left + CELL_W + o * 0.6, top + 3.0 - o, 9.0 + o, -6.0 - o * 0.6)
    } else {
        (left - o * 0.6, top + CELL_H - 3.0 + o, -9.0 - o, 6.0 + o * 0.6)
    };
    let _ = write!(
        out,
        r#"<path d="M{:.1},{:.1} L{cx:.1},{cy:.1} L{cx:.1},{:.1}" fill="none" stroke="{}" stroke-width="1.8"/>"#,
        cx + dx,
        cy,
        cy + dy,
        hex(color)
    );
}

fn peptide_svg(out: &mut String, a: &Annotation, peptide: usize, x: f64, top: f64, linked: &[usize]) {
    let c = &a.coverage[peptide];
    for (i, residue) in c.residues.iter().enumerate() {
        let left = x + CELL_W * i as f64;
        if linked.contains(&i) {
            let _ = write!(out, r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" rx="3" fill="#000000" fill-opacity="0.12"/>"##, left + 1.0, top + 5.0, CELL_W - 2.0, CELL_H - 10.0);
        }
        let _ = write!(out, r#"<text x="{:.1}" y="{:.1}" font-family="{MONO}" font-size="16" text-anchor="middle">{}</text>"#, left + CELL_W / 2.0, top + CELL_H / 2.0 + 6.0, esc(residue));
        if c.n_term[i] {
            hook_svg(out, left, top, true, 0, series_color(Series::B));
        }
        if c.c_term[i] {
            hook_svg(out, left, top, false, 0, series_color(Series::Y));
        }
        for n_term in [true, false] {
            let mut stubs: Vec<usize> = a
                .stub_ions
                .iter()
                .filter(|s| s.peptide == peptide)
                .filter(|s| match s.kind {
                    StubIonKind::NTerm(j) => n_term && j == i,
                    StubIonKind::CTerm(j) => !n_term && j == i,
                    StubIonKind::Peptide => false,
                })
                .map(|s| s.stub)
                .collect();
            stubs.sort_unstable();
            stubs.dedup();
            for (level, stub) in stubs.into_iter().enumerate() {
                hook_svg(out, left, top, n_term, level + 1, stub_color(stub));
            }
        }
    }
}

fn coverage_svg(out: &mut String, a: &Annotation, y: f64) -> f64 {
    let x0 = 110.0;
    let label = |out: &mut String, text: &str, top: f64| {
        let _ = write!(out, r#"<text x="24" y="{:.1}" font-family="{FONT}" font-size="13">{}</text>"#, top + CELL_H / 2.0 + 5.0, esc(text));
    };
    let pad = 12.0;
    match &a.link {
        Some(link) if link.kind == Kind::Crosslink && a.coverage.len() == 2 => {
            let (s0, s1) = (link.sites[0].1, link.sites[1].1);
            let shift = s0.max(s1);
            let top0 = y + pad;
            let bridge = 52.0;
            let top1 = top0 + CELL_H + 2.0 * pad + bridge;
            label(out, "peptide 1", top0);
            label(out, "peptide 2", top1);
            peptide_svg(out, a, 0, x0 + CELL_W * (shift - s0) as f64, top0, &[s0]);
            peptide_svg(out, a, 1, x0 + CELL_W * (shift - s1) as f64, top1, &[s1]);
            bridge_svg(out, a, link, x0 + CELL_W * (shift as f64 + 0.5), top0 + CELL_H - 2.0, top1 + 2.0);
            top1 + CELL_H + pad
        }
        Some(link) if link.kind == Kind::Loop && a.coverage.len() == 1 => {
            let (i0, i1) = (link.sites[0].1.min(link.sites[1].1), link.sites[0].1.max(link.sites[1].1));
            let top = y + pad;
            label(out, "loop-link", top);
            peptide_svg(out, a, 0, x0, top, &[i0, i1]);
            let (xa, xb) = (x0 + CELL_W * (i0 as f64 + 0.5), x0 + CELL_W * (i1 as f64 + 0.5));
            let (ya, yb) = (top + CELL_H + 4.0, top + CELL_H + pad + 8.0);
            let _ = write!(out, r##"<path d="M{xa:.1},{ya:.1} V{yb:.1} H{xb:.1} V{ya:.1}" fill="none" stroke="#888888" stroke-width="2.5"/>"##);
            let _ = write!(out, r#"<text x="{:.1}" y="{:.1}" font-family="{FONT}" font-size="12" text-anchor="middle">{}</text>"#, (xa + xb) / 2.0, yb + 14.0, esc(&link.linker));
            yb + 22.0
        }
        _ => {
            let mut top = y + pad;
            for k in 0..a.coverage.len() {
                if a.coverage.len() > 1 {
                    label(out, &format!("peptide {}", k + 1), top);
                }
                peptide_svg(out, a, k, x0, top, &[]);
                top += CELL_H + 2.0 * pad;
            }
            top
        }
    }
}

/// The bridge between the linked residues, labelled with the linker; a
/// triangle per peptide + stub ion found, pointing at its peptide
fn bridge_svg(out: &mut String, a: &Annotation, link: &Link, x: f64, y_top: f64, y_bottom: f64) {
    let _ = write!(out, r##"<line x1="{x:.1}" y1="{y_top:.1}" x2="{x:.1}" y2="{y_bottom:.1}" stroke="#888888" stroke-width="2.5"/>"##);
    let _ = write!(
        out,
        r#"<text x="{:.1}" y="{:.1}" font-family="{FONT}" font-size="12" text-anchor="end">{}</text>"#,
        x - 9.0,
        (y_top + y_bottom) / 2.0 + 4.0,
        esc(&link.linker)
    );
    for (k, cy, up) in [(0, y_top + 13.0, true), (1, y_bottom - 13.0, false)] {
        let mut stubs: Vec<usize> = a
            .stub_ions
            .iter()
            .filter(|s| s.peptide == k && s.kind == StubIonKind::Peptide)
            .map(|s| s.stub)
            .collect();
        stubs.sort_unstable();
        stubs.dedup();
        let mut cx = x + 12.0;
        for stub in stubs {
            let color = hex(stub_color(stub));
            let (tip, base) = if up { (-5.0, 4.0) } else { (5.0, -4.0) };
            let _ = write!(
                out,
                r#"<polygon points="{cx:.1},{:.1} {:.1},{:.1} {:.1},{:.1}" fill="{color}"/>"#,
                cy + tip,
                cx - 5.0,
                cy + base,
                cx + 5.0,
                cy + base
            );
            let name = &link.stubs[stub].name;
            let _ = write!(out, r#"<text x="{:.1}" y="{:.1}" font-family="{FONT}" font-size="12" fill="{color}">{}</text>"#, cx + 8.0, cy + 4.0, esc(name));
            cx += 8.0 + name.chars().count() as f64 * 7.0 + 10.0;
        }
    }
}

fn legend_svg(out: &mut String, a: &Annotation, y: f64) -> f64 {
    let mut items: Vec<(Color32, String)> = Vec::new();
    let mut present: Vec<Series> = a
        .peaks
        .iter()
        .flat_map(|p| p.labels.iter().map(|l| l.series))
        .filter(|s| !matches!(s, Series::Stub(_)))
        .collect();
    present.sort_by_key(|s| s.name());
    present.dedup();
    for s in present {
        items.push((series_color(s), s.name().to_string()));
    }
    if let Some(link) = &a.link {
        for (j, stub) in link.stubs.iter().enumerate() {
            if a.stub_ions.iter().any(|s| s.stub == j) {
                items.push((stub_color(j), format!("{} {} ({:+.3})", stub.name, stub.description, stub.mass)));
            }
        }
    }
    items.push((UNMATCHED, "unmatched".into()));
    let (mut x, mut y) = (24.0, y + 14.0);
    for (color, text) in items {
        let width = 16.0 + text.chars().count() as f64 * 6.6 + 18.0;
        if x + width > W - 20.0 {
            x = 24.0;
            y += 18.0;
        }
        let _ = write!(out, r#"<rect x="{x:.1}" y="{:.1}" width="10" height="10" fill="{}"/>"#, y - 9.0, hex(color));
        let _ = write!(out, r#"<text x="{:.1}" y="{y:.1}" font-family="{FONT}" font-size="12">{}</text>"#, x + 15.0, esc(&text));
        x += width;
    }
    y
}

/// Write the figure: the SVG itself, or rendered at `dpi` to PNG / TIFF
pub fn export(svg: &str, format: Format, dpi: u32, background: bool, path: &Path) -> Result<(), String> {
    if format == Format::Svg {
        return std::fs::write(path, svg).map_err(|e| e.to_string());
    }
    let image = rasterize(svg, dpi)?;
    match (format, background) {
        // opaque figures without an alpha channel, as journals ask for
        (_, true) => image::DynamicImage::ImageRgba8(image).to_rgb8().save(path),
        _ => image.save(path),
    }
    .map_err(|e| e.to_string())
}

/// Render the SVG at `dpi` (straight, not premultiplied, RGBA)
pub fn rasterize(svg: &str, dpi: u32) -> Result<image::RgbaImage, String> {
    use resvg::{tiny_skia, usvg};
    let mut options = usvg::Options::default();
    options.fontdb_mut().load_system_fonts();
    let tree = usvg::Tree::from_str(svg, &options).map_err(|e| e.to_string())?;
    // usvg sizes the figure (given in inches) in CSS pixels, `options.dpi` per inch
    let scale = dpi as f32 / options.dpi;
    let size = tree.size();
    let (w, h) = ((size.width() * scale).ceil() as u32, (size.height() * scale).ceil() as u32);
    let mut pixmap = tiny_skia::Pixmap::new(w, h).ok_or("figure too large")?;
    resvg::render(&tree, tiny_skia::Transform::from_scale(scale, scale), &mut pixmap.as_mut());
    let mut rgba = Vec::with_capacity((w * h * 4) as usize);
    for p in pixmap.pixels() {
        let c = p.demultiply();
        rgba.extend_from_slice(&[c.red(), c.green(), c.blue(), c.alpha()]);
    }
    image::RgbaImage::from_raw(w, h, rgba).ok_or_else(|| "could not build the image".into())
}

#[cfg(test)]
mod test {
    use super::*;

    fn options(background: bool) -> FigureOptions {
        FigureOptions {
            background,
            title: Some("GLEKVSTAAMQR (4) – TDYHFPNEWKLIR (10)".into()),
            label_all: true,
            mz_range: None,
            coverage: true,
            dpi: 300,
        }
    }

    #[test]
    fn charges_are_raised() {
        assert_eq!(label_svg("y9"), "y9");
        assert_eq!(
            label_svg("1:y9²+-H2O"),
            r#"1:y9<tspan baseline-shift="super" font-size="75%">2+</tspan>-H2O"#
        );
        assert_eq!(label_svg("P+A²+"), r#"P+A<tspan baseline-shift="super" font-size="75%">2+</tspan>"#);
        assert_eq!(label_svg("a<b"), "a&lt;b");
    }

    #[test]
    fn ticks_are_round() {
        assert_eq!(tick_step(1400.0, 8.0), 200.0);
        assert_eq!(tick_step(95.0, 8.0), 20.0);
        assert_eq!(tick_step(12.0, 8.0), 2.0);
    }

    #[test]
    fn crosslink_figure_holds_bridge_stubs_and_background_choice() {
        let (_, a) = super::super::model::test::cleaved_dsso();
        let with = svg(&a, &options(true));
        let without = svg(&a, &options(false));
        assert!(with.contains(r#"fill="white""#) && !without.contains(r#"fill="white""#));
        assert!(with.starts_with("<svg") && with.contains(r#"width="7in""#));
        // the bridge, the linker, the stub triangles and the stub-labelled peaks
        assert!(with.contains(">DSSO</text>"));
        assert_eq!(with.matches("<polygon").count(), 4);
        assert!(with.contains("y4+S") && with.contains("P+A"));
        assert!(with.contains("A alkene") && with.contains("S sulfenic acid"));
        // well-formed: rendering parses it
        let image = rasterize(&with, 150).unwrap();
        assert_eq!(image.width(), 1050);
        // transparent figure: the corner pixel is see-through, opaque otherwise
        assert_eq!(rasterize(&without, 150).unwrap().get_pixel(0, 0)[3], 0);
        assert_eq!(image.get_pixel(0, 0).0, [255, 255, 255, 255]);
    }

    #[test]
    fn exports_every_format() {
        let (_, a) = super::super::model::test::cleaved_dsso();
        let dir = std::env::temp_dir().join(format!("sage-gui-figure-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for format in Format::ALL {
            for background in [true, false] {
                let path = dir.join(format!("figure_{background}.{}", format.extension()));
                export(&svg(&a, &options(background)), format, 300, background, &path).unwrap();
                assert!(std::fs::metadata(&path).unwrap().len() > 1000, "{}", path.display());
                if format != Format::Svg {
                    let img = image::open(&path).unwrap();
                    assert_eq!(img.width(), 2100, "7 in at 300 dpi");
                    assert_eq!(img.color().has_alpha(), !background);
                }
            }
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
