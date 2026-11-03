//! Mobiligram of a timsTOF (Bruker `.d`) run: every precursor Sage searched,
//! placed by m/z and inverse ion mobility (1/K0). The marker shape gives the
//! charge state, the colour what the precursor was identified as: crosslinks
//! (between proteins, within a protein, loop-links), linear peptides, or
//! nothing at 1% FDR. Drawn in the viewer and exported as SVG or PDF.

use super::figure;
use super::model::{Kind, Match, Results};
use eframe::egui::Color32;
use std::collections::HashMap;
use std::fmt::Write;

/// Proton mass, for m/z from Sage's neutral precursor mass
const PROTON: f64 = 1.007_276_466_6;
/// Matches at or below this spectrum q-value count as identified
pub const Q_CUTOFF: f32 = 0.01;

const WIDTH_IN: f64 = 7.0;
const UNITS_PER_IN: f64 = 150.0;
const W: f64 = WIDTH_IN * UNITS_PER_IN;
const FONT: &str = "Arial, Helvetica, sans-serif";

/// What a precursor was identified as, in drawing order (later on top)
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Class {
    Unidentified,
    Linear,
    Loop,
    IntraProtein,
    InterProtein,
}

impl Class {
    pub const ALL: [Class; 5] = [Class::Unidentified, Class::Linear, Class::Loop, Class::IntraProtein, Class::InterProtein];

    pub fn label(self) -> &'static str {
        match self {
            Class::InterProtein => "crosslink, between proteins",
            Class::IntraProtein => "crosslink, within a protein",
            Class::Loop => "loop-link",
            Class::Linear => "linear peptide",
            Class::Unidentified => "not identified at 1% FDR",
        }
    }

    /// Fill colour: crosslinks between proteins red (#FF523B), within a
    /// protein orange (#FF9D35), loop-links pink; linear peptides blue,
    /// unidentified precursors grey
    pub fn color(self) -> Color32 {
        match self {
            Class::InterProtein => Color32::from_rgb(0xFF, 0x52, 0x3B),
            Class::IntraProtein => Color32::from_rgb(0xFF, 0x9D, 0x35),
            Class::Loop => Color32::from_rgb(204, 121, 167),
            Class::Linear => Color32::from_rgb(0, 114, 178),
            Class::Unidentified => Color32::from_rgb(190, 190, 190),
        }
    }

    /// Fill and outline as SVG attributes: crosslinks solid with a dark
    /// outline, linear peptides blue at 50% opacity, unidentified light grey
    fn svg_style(self) -> String {
        let color = hex(self.color());
        match self {
            Class::Unidentified => format!(r#"fill="{color}" fill-opacity="0.55""#),
            Class::Linear => format!(r#"fill="{color}" fill-opacity="0.5""#),
            _ => format!(r##"fill="{color}" stroke="#222222" stroke-width="0.6""##),
        }
    }

    pub fn is_crosslink(self) -> bool {
        matches!(self, Class::InterProtein | Class::IntraProtein | Class::Loop)
    }
}

/// Marker shape of a charge state
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    Cross,
    Circle,
    Triangle,
    Square,
    Diamond,
    TriangleDown,
}

impl Shape {
    /// 1+ cross, 2+ circle, 3+ triangle, 4+ square, 5+ diamond, 6+ and up
    /// triangle pointing down
    pub fn of(charge: u8) -> Shape {
        match charge {
            0 | 1 => Shape::Cross,
            2 => Shape::Circle,
            3 => Shape::Triangle,
            4 => Shape::Square,
            5 => Shape::Diamond,
            _ => Shape::TriangleDown,
        }
    }

    /// SVG path of the marker centred on 0,0 with radius `r`
    fn path(self, r: f64) -> String {
        match self {
            Shape::Circle => format!("M{r:.2},0 A{r:.2},{r:.2} 0 1 1 {:.2},0 A{r:.2},{r:.2} 0 1 1 {r:.2},0 Z", -r),
            Shape::Square => {
                let s = r * 0.85;
                format!("M{:.2},{:.2} H{s:.2} V{s:.2} H{:.2} Z", -s, -s, -s)
            }
            Shape::Diamond => {
                let s = r * 1.15;
                format!("M0,{:.2} L{s:.2},0 L0,{s:.2} L{:.2},0 Z", -s, -s)
            }
            Shape::Triangle => {
                let s = r * 1.2;
                format!("M0,{:.2} L{:.2},{:.2} L{:.2},{:.2} Z", -s, s * 0.87, s * 0.5, -s * 0.87, s * 0.5)
            }
            Shape::TriangleDown => {
                let s = r * 1.2;
                format!("M0,{s:.2} L{:.2},{:.2} L{:.2},{:.2} Z", s * 0.87, -s * 0.5, -s * 0.87, -s * 0.5)
            }
            Shape::Cross => {
                // a plus sign: arms `a` long, `b` half thick
                let (a, b, na, nb) = (r, r * 0.3, -r, -r * 0.3);
                format!(
                    "M{nb:.2},{na:.2} H{b:.2} V{nb:.2} H{a:.2} V{b:.2} H{b:.2} V{a:.2} H{nb:.2} V{b:.2} H{na:.2} V{nb:.2} H{nb:.2} Z"
                )
            }
        }
    }
}

/// One precursor of the mobiligram
#[derive(Clone, Debug, PartialEq)]
pub struct Point {
    pub file: String,
    pub mz: f64,
    /// 1/K0, V·s/cm²
    pub mobility: f32,
    pub charge: u8,
    pub class: Class,
    pub rt: f32,
    /// What it was identified as (empty when not identified)
    pub title: String,
}

fn class_of(m: &Match) -> Class {
    let identified = m.target && m.spectrum_q <= Q_CUTOFF;
    match m.kind {
        _ if !identified => Class::Unidentified,
        Kind::Linear => Class::Linear,
        Kind::Loop => Class::Loop,
        Kind::Crosslink if m.intra_protein => Class::IntraProtein,
        Kind::Crosslink => Class::InterProtein,
    }
}

/// Every precursor with an ion mobility, one per spectrum: the best
/// identification wins (a crosslink over a linear peptide, the rank-1 match of
/// each search)
pub fn points(results: &Results) -> Vec<Point> {
    let mut best: HashMap<(&str, &str), (Class, &Match)> = HashMap::new();
    for m in results.matches.iter().filter(|m| m.rank == 1) {
        let Some(mobility) = m.ion_mobility else { continue };
        if mobility <= 0.0 || m.expmass <= 0.0 || m.charge == 0 {
            continue;
        }
        let class = class_of(m);
        let entry = best.entry((m.file.as_str(), m.scan.as_str())).or_insert((class, m));
        if class > entry.0 {
            *entry = (class, m);
        }
    }
    let mut points: Vec<Point> = best
        .into_values()
        .map(|(class, m)| Point {
            file: m.file.clone(),
            mz: m.expmass as f64 / m.charge as f64 + PROTON,
            mobility: m.ion_mobility.unwrap_or_default(),
            charge: m.charge,
            class,
            rt: m.rt,
            title: if class == Class::Unidentified { String::new() } else { m.title() },
        })
        .collect();
    points.sort_by(|a, b| a.class.cmp(&b.class).then(a.mz.total_cmp(&b.mz)));
    points
}

#[derive(Clone, Debug)]
pub struct MobiligramOptions {
    /// White background; otherwise transparent
    pub background: bool,
    pub title: Option<String>,
    /// Also draw the precursors not identified at 1% FDR
    pub unidentified: bool,
    /// Only this spectrum file; all files when `None`
    pub file: Option<String>,
}

/// The precursors `opts` keeps
pub fn selected<'a>(points: &'a [Point], opts: &MobiligramOptions) -> Vec<&'a Point> {
    points
        .iter()
        .filter(|p| opts.unidentified || p.class != Class::Unidentified)
        .filter(|p| opts.file.as_ref().is_none_or(|f| *f == p.file))
        .collect()
}

/// A round range on `step` around `lo..hi`, with 2% room so no point sits on the frame
fn round_range(lo: f64, hi: f64, step: f64) -> (f64, f64) {
    let pad = (hi - lo) * 0.02;
    let a = ((lo - pad) / step + 1e-9).floor() * step;
    let b = ((hi + pad) / step - 1e-9).ceil() * step;
    (a, b.max(a + step))
}

/// The mobiligram as SVG, 7 inches wide
pub fn svg(points: &[Point], opts: &MobiligramOptions) -> String {
    let shown = selected(points, opts);
    let (mz_lo, mz_hi, im_lo, im_hi) = shown.iter().fold((f64::MAX, f64::MIN, f64::MAX, f64::MIN), |(a, b, c, d), p| {
        (a.min(p.mz), b.max(p.mz), c.min(p.mobility as f64), d.max(p.mobility as f64))
    });
    let (mz_lo, mz_hi) = if shown.is_empty() { (100.0, 1700.0) } else { round_range(mz_lo, mz_hi, 100.0) };
    let (im_lo, im_hi) = if shown.is_empty() { (0.6, 1.6) } else { round_range(im_lo, im_hi, 0.1) };

    let mut body = String::new();
    let mut y = 16.0;
    if let Some(title) = &opts.title {
        let _ = write!(body, r#"<text x="20" y="{:.1}" font-family="{FONT}" font-size="15" font-weight="bold">{}</text>"#, y + 12.0, esc(title));
        y += 26.0;
    }

    // plot area on the left, legend on the right
    let (x0, x1) = (84.0, 770.0);
    let (top, bottom) = (y + 14.0, y + 560.0);
    let px = |mz: f64| x0 + (mz - mz_lo) / (mz_hi - mz_lo) * (x1 - x0);
    let py = |im: f64| bottom - (im - im_lo) / (im_hi - im_lo) * (bottom - top);
    let axis = "#333333";

    // marker shapes, reused by every point
    body.push_str("<defs>");
    for charge in 1..=6u8 {
        let shape = Shape::of(charge);
        let _ = write!(body, r#"<path id="z{charge}" d="{}"/>"#, shape.path(3.4));
        let _ = write!(body, r#"<path id="zl{charge}" d="{}"/>"#, shape.path(5.0));
    }
    body.push_str("</defs>");

    // frame, grid and ticks
    let _ = write!(body, r#"<g font-family="{FONT}" font-size="12" fill="{axis}">"#);
    let step = figure::tick_step(mz_hi - mz_lo, 8.0);
    let mut t = (mz_lo / step).ceil() * step;
    while t <= mz_hi + 1e-9 {
        let x = px(t);
        let _ = write!(body, r##"<line x1="{x:.1}" y1="{top:.1}" x2="{x:.1}" y2="{bottom:.1}" stroke="#e6e6e6" stroke-width="0.8"/>"##);
        let _ = write!(body, r#"<line x1="{x:.1}" y1="{bottom:.1}" x2="{x:.1}" y2="{:.1}" stroke="{axis}" stroke-width="1"/>"#, bottom + 5.0);
        let _ = write!(body, r#"<text x="{x:.1}" y="{:.1}" text-anchor="middle">{}</text>"#, bottom + 19.0, t.round() as i64);
        t += step;
    }
    let step = figure::tick_step(im_hi - im_lo, 6.0);
    let mut t = (im_lo / step).ceil() * step;
    while t <= im_hi + 1e-9 {
        let yy = py(t);
        let _ = write!(body, r##"<line x1="{x0:.1}" y1="{yy:.1}" x2="{x1:.1}" y2="{yy:.1}" stroke="#e6e6e6" stroke-width="0.8"/>"##);
        let _ = write!(body, r#"<line x1="{:.1}" y1="{yy:.1}" x2="{x0:.1}" y2="{yy:.1}" stroke="{axis}" stroke-width="1"/>"#, x0 - 5.0);
        let _ = write!(body, r#"<text x="{:.1}" y="{:.1}" text-anchor="end">{t:.2}</text>"#, x0 - 8.0, yy + 4.0);
        t += step;
    }
    let _ = write!(body, r#"<rect x="{x0:.1}" y="{top:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="{axis}" stroke-width="1.2"/>"#, x1 - x0, bottom - top);
    let _ = write!(body, r#"<text x="{:.1}" y="{:.1}" text-anchor="middle"><tspan font-style="italic">m/z</tspan></text>"#, (x0 + x1) / 2.0, bottom + 40.0);
    let _ = write!(
        body,
        r#"<text transform="translate(28,{:.1}) rotate(-90)" text-anchor="middle">1/K<tspan baseline-shift="sub" font-size="75%">0</tspan> (V·s/cm<tspan baseline-shift="super" font-size="75%">2</tspan>)</text>"#,
        (top + bottom) / 2.0
    );
    body.push_str("</g>");

    // points, by class: unidentified underneath, crosslinks on top
    let _ = write!(body, r#"<g clip-path="url(#plot)">"#);
    for class in Class::ALL {
        let style = class.svg_style();
        let large = if class.is_crosslink() { "l" } else { "" };
        let _ = write!(body, r#"<g {style}>"#);
        for p in shown.iter().filter(|p| p.class == class) {
            let _ = write!(
                body,
                r##"<use xlink:href="#z{large}{}" x="{:.1}" y="{:.1}"/>"##,
                p.charge.clamp(1, 6),
                px(p.mz),
                py(p.mobility as f64)
            );
        }
        body.push_str("</g>");
    }
    body.push_str("</g>");

    // legend: colours (with counts), then charge shapes
    let lx = x1 + 26.0;
    let mut ly = top + 6.0;
    let _ = write!(body, r#"<g font-family="{FONT}" font-size="12" fill="{axis}">"#);
    let _ = write!(body, r#"<text x="{lx:.1}" y="{ly:.1}" font-weight="bold">Precursor</text>"#);
    ly += 20.0;
    for class in Class::ALL.into_iter().rev() {
        if class == Class::Unidentified && !opts.unidentified {
            continue;
        }
        // every class is listed, with n = 0 when none was found (e.g. no
        // loop-links), so an empty class does not look left out
        let n = shown.iter().filter(|p| p.class == class).count();
        let _ = write!(body, r##"<use xlink:href="#zl2" x="{:.1}" y="{:.1}" {}/>"##, lx + 6.0, ly - 4.0, class.svg_style());
        let (first, second) = match class.label().split_once(", ") {
            Some((a, b)) => (format!("{a},"), Some(b)),
            None => (class.label().to_string(), None),
        };
        let _ = write!(body, r#"<text x="{:.1}" y="{ly:.1}">{}</text>"#, lx + 18.0, esc(&first));
        if let Some(second) = second {
            ly += 15.0;
            let _ = write!(body, r#"<text x="{:.1}" y="{ly:.1}">{}</text>"#, lx + 18.0, esc(second));
        }
        let _ = write!(body, r##"<text x="{:.1}" y="{:.1}" fill="#777777" font-size="11">n = {n}</text>"##, lx + 18.0, ly + 14.0);
        ly += 34.0;
    }
    ly += 8.0;
    let _ = write!(body, r#"<text x="{lx:.1}" y="{ly:.1}" font-weight="bold">Charge</text>"#);
    ly += 20.0;
    let charges: Vec<u8> = {
        let mut c: Vec<u8> = shown.iter().map(|p| p.charge.clamp(1, 6)).collect();
        c.sort_unstable();
        c.dedup();
        c
    };
    for charge in charges {
        let label = if charge == 6 { "6+ and up".to_string() } else { format!("{charge}+") };
        let _ = write!(body, r##"<use xlink:href="#zl{charge}" x="{:.1}" y="{:.1}" fill="none" stroke="#333333" stroke-width="1.2"/>"##, lx + 6.0, ly - 4.0);
        let _ = write!(body, r#"<text x="{:.1}" y="{ly:.1}">{label}</text>"#, lx + 18.0);
        ly += 20.0;
    }
    body.push_str("</g>");

    let height = (bottom + 56.0).max(ly + 12.0);
    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{WIDTH_IN}in" height="{:.3}in" viewBox="0 0 {W} {height:.1}">"#,
        height / UNITS_PER_IN
    );
    let _ = write!(svg, r#"<clipPath id="plot"><rect x="{x0:.1}" y="{top:.1}" width="{:.1}" height="{:.1}"/></clipPath>"#, x1 - x0, bottom - top);
    if opts.background {
        svg.push_str(r#"<rect width="100%" height="100%" fill="white"/>"#);
    }
    svg.push_str(&body);
    svg.push_str("</svg>\n");
    svg
}

fn hex(c: Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

#[cfg(test)]
mod test {
    use super::*;

    fn point(mz: f64, mobility: f32, charge: u8, class: Class) -> Point {
        Point { file: "run.d".into(), mz, mobility, charge, class, rt: 10.0, title: String::new() }
    }

    fn options(unidentified: bool) -> MobiligramOptions {
        MobiligramOptions { background: true, title: Some("Pol II".into()), unidentified, file: None }
    }

    #[test]
    fn charges_get_their_own_shapes() {
        let shapes: Vec<Shape> = (1..=8).map(Shape::of).collect();
        assert_eq!(
            shapes,
            [
                Shape::Cross,
                Shape::Circle,
                Shape::Triangle,
                Shape::Square,
                Shape::Diamond,
                Shape::TriangleDown,
                Shape::TriangleDown,
                Shape::TriangleDown
            ]
        );
    }

    #[test]
    fn figure_has_points_legend_and_renders() {
        let points = vec![
            point(512.3, 0.81, 2, Class::Unidentified),
            point(640.8, 0.92, 2, Class::Linear),
            point(733.1, 0.98, 3, Class::Linear),
            point(855.4, 1.02, 4, Class::InterProtein),
            point(901.9, 1.05, 4, Class::IntraProtein),
            point(1002.2, 1.10, 5, Class::Loop),
        ];
        let with = svg(&points, &options(true));
        assert_eq!(with.matches("<use xlink:href=\"#z").count(), 6 + 5 + 4, "points, colour and charge keys");
        for class in Class::ALL {
            assert!(with.contains(&hex(class.color())), "{class:?}");
        }
        assert!(with.contains("crosslink,") && with.contains("between proteins") && with.contains("n = 2"));
        assert!(with.contains(">2+<") && with.contains(">5+<") && !with.contains(">1+<"));
        // linear peptides: blue at 50%, no outline
        assert!(with.contains(r##"<g fill="#0072b2" fill-opacity="0.5">"##));
        // a class with no precursors is still listed, with n = 0
        let no_loops: Vec<Point> = points.iter().filter(|p| p.class != Class::Loop).cloned().collect();
        let legend = svg(&no_loops, &options(true));
        assert!(legend.contains(">loop-link</text>") && legend.contains("n = 0"));
        let without = svg(&points, &options(false));
        assert!(!without.contains("not identified"));
        assert_eq!(without.matches("<use xlink:href=\"#z").count(), 5 + 4 + 4);

        // renders, and converts to a one-page PDF
        let image = figure::rasterize(&with, 150).unwrap();
        assert_eq!(image.width(), 1050);
        // a crosslink marker is drawn in its colour where the point is
        assert!(image.pixels().any(|p| p.0[..3] == [0xFF, 0x52, 0x3B]));
        assert!(figure::pdf(&with).unwrap().starts_with(b"%PDF"));
    }

    /// A real results folder (`SAGE_GUI_RESULTS`): the mobiligram as SVG, PDF
    /// and a PNG preview in `SAGE_GUI_FIGURES`
    #[test]
    #[ignore]
    fn real_mobiligram() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let out = std::path::PathBuf::from(std::env::var("SAGE_GUI_FIGURES").expect("set SAGE_GUI_FIGURES"));
        let results = Results::load(std::path::Path::new(&folder)).unwrap();
        let points = points(&results);
        for class in Class::ALL {
            let of: Vec<_> = points.iter().filter(|p| p.class == class).collect();
            let (lo, hi) = of.iter().fold((f32::MAX, 0.0f32), |(a, b), p| (a.min(p.mobility), b.max(p.mobility)));
            println!("{:28} {:6} precursors, 1/K0 {lo:.3}–{hi:.3}", class.label(), of.len());
        }
        std::fs::create_dir_all(&out).unwrap();
        for (name, unidentified, background) in [("all", true, true), ("identified", false, true), ("transparent", true, false)] {
            let opts = MobiligramOptions {
                background,
                title: Some(format!("Mobiligram · Pol II DSSO · {name}")),
                unidentified,
                file: None,
            };
            let svg = svg(&points, &opts);
            for format in [figure::Format::Svg, figure::Format::Pdf, figure::Format::Png] {
                let path = out.join(format!("mobiligram_{name}.{}", format.extension()));
                figure::export(&svg, format, 150, background, &path).unwrap();
                println!("{} ({} kB)", path.display(), std::fs::metadata(&path).unwrap().len() / 1024);
            }
        }
    }

    #[test]
    fn axes_round_to_the_data() {
        assert_eq!(round_range(512.3, 1002.2, 100.0), (500.0, 1100.0));
        assert_eq!(round_range(520.0, 980.0, 100.0), (500.0, 1000.0));
        let (lo, hi) = round_range(0.81, 1.10, 0.1);
        assert!((lo - 0.8).abs() < 1e-9 && (hi - 1.2).abs() < 1e-9, "{lo} {hi}");
    }
}
