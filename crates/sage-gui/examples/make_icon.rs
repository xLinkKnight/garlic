//! Draws the application icon, a smiling cartoon garlic bulb, from simple shapes
//! (signed distance fields, so every size is rendered sharp rather than scaled)
//! and writes `assets/garlic.png` (window icon) and `assets/garlic.ico` (exe icon).
//!
//!     cargo run -p sage-gui --example make_icon

use image::{ImageEncoder, Rgba, RgbaImage};
use std::path::Path;

type Color = [f32; 4];

const BACKGROUND_TOP: Color = [1.0, 0.70, 0.30, 1.0];
const BACKGROUND_BOTTOM: Color = [0.96, 0.52, 0.18, 1.0];
const OUTLINE: Color = [0.33, 0.20, 0.13, 1.0];
const BULB_LIGHT: Color = [1.0, 0.97, 0.91, 1.0];
const BULB_SHADE: Color = [0.95, 0.85, 0.72, 1.0];
const CLOVE_LINE: Color = [0.88, 0.74, 0.58, 1.0];
const LEAF: Color = [0.42, 0.74, 0.27, 1.0];
const LEAF_DARK: Color = [0.29, 0.56, 0.18, 1.0];
const CHEEK: Color = [1.0, 0.55, 0.55, 0.55];
const WHITE: Color = [1.0, 1.0, 1.0, 1.0];

fn mix(a: Color, b: Color, t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    [0, 1, 2, 3].map(|i| a[i] + (b[i] - a[i]) * t)
}

/// Paint `color` over `dst` with coverage `cov` (0..1)
fn over(dst: Color, color: Color, cov: f32) -> Color {
    let a = color[3] * cov.clamp(0.0, 1.0);
    let out_a = a + dst[3] * (1.0 - a);
    if out_a <= 0.0 {
        return [0.0; 4];
    }
    let mut c = [0.0; 4];
    for i in 0..3 {
        c[i] = (color[i] * a + dst[i] * dst[3] * (1.0 - a)) / out_a;
    }
    c[3] = out_a;
    c
}

fn circle(x: f32, y: f32, cx: f32, cy: f32, r: f32) -> f32 {
    ((x - cx).powi(2) + (y - cy).powi(2)).sqrt() - r
}

fn ellipse(x: f32, y: f32, cx: f32, cy: f32, rx: f32, ry: f32) -> f32 {
    // approximate distance, good enough for soft shapes
    let k = (((x - cx) / rx).powi(2) + ((y - cy) / ry).powi(2)).sqrt();
    (k - 1.0) * rx.min(ry)
}

fn rounded_rect(x: f32, y: f32, half: f32, r: f32) -> f32 {
    let qx = (x - 0.5).abs() - half + r;
    let qy = (y - 0.5).abs() - half + r;
    (qx.max(0.0).powi(2) + qy.max(0.0).powi(2)).sqrt() + qx.max(qy).min(0.0) - r
}

/// Smooth union: blends two shapes into one rounded silhouette
fn smooth_min(a: f32, b: f32, k: f32) -> f32 {
    let h = (0.5 + 0.5 * (b - a) / k).clamp(0.0, 1.0);
    b + (a - b) * h - k * h * (1.0 - h)
}

/// Distance to a line segment
fn segment(x: f32, y: f32, ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    let (px, py, dx, dy) = (x - ax, y - ay, bx - ax, by - ay);
    let t = ((px * dx + py * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
    ((px - dx * t).powi(2) + (py - dy * t).powi(2)).sqrt()
}

/// A leaf: the lens between two circles, rotated by `angle` around (cx, cy)
fn leaf(x: f32, y: f32, cx: f32, cy: f32, len: f32, width: f32, angle: f32) -> f32 {
    let (s, c) = angle.sin_cos();
    let (u, v) = ((x - cx) * c + (y - cy) * s, -(x - cx) * s + (y - cy) * c);
    // lens of length `len` along u, `width` wide
    let r = (len * len / 4.0 + width * width / 4.0) / width;
    let off = r - width / 2.0;
    circle(u, v, len / 2.0, off, r).max(circle(u, v, len / 2.0, -off, r))
}

/// The garlic bulb's outline: a round body that tapers into a pointed top
fn bulb(x: f32, y: f32) -> f32 {
    let body = ellipse(x, y, 0.5, 0.62, 0.285, 0.255);
    let neck = ellipse(x, y, 0.5, 0.36, 0.085, 0.13);
    smooth_min(body, neck, 0.11)
}

/// Colour of the icon at (x, y) in 0..1, for a pixel of size `px`
fn shade(x: f32, y: f32, px: f32) -> Color {
    let fill = |d: f32| (0.5 - d / px).clamp(0.0, 1.0);
    let stroke = |d: f32, w: f32| (0.5 - (d.abs() - w / 2.0) / px).clamp(0.0, 1.0);
    let line = 0.022_f32.max(1.1 * px);

    // rounded square background with a warm gradient
    let mut c: Color = [0.0; 4];
    c = over(c, mix(BACKGROUND_TOP, BACKGROUND_BOTTOM, y), fill(rounded_rect(x, y, 0.47, 0.11)));

    // sprout: a short stem and two leaves
    let stem = segment(x, y, 0.5, 0.25, 0.5, 0.17) - 0.016;
    let leaf_l = leaf(x, y, 0.5, 0.21, 0.17, 0.085, -2.35);
    let leaf_r = leaf(x, y, 0.5, 0.215, 0.15, 0.075, -0.75);
    for (d, col) in [(stem, LEAF_DARK), (leaf_l, LEAF), (leaf_r, LEAF)] {
        c = over(c, col, fill(d));
    }
    // leaf veins and outlines
    let vein_l = segment(x, y, 0.5, 0.21, 0.405, 0.13);
    let vein_r = segment(x, y, 0.5, 0.215, 0.588, 0.145);
    c = over(c, LEAF_DARK, fill(vein_l.min(vein_r) - 0.005));
    c = over(c, OUTLINE, stroke(leaf_l.min(leaf_r).min(stem), line * 0.8));

    // roots below the bulb
    for dx in [-0.05_f32, 0.0, 0.05] {
        let d = segment(x, y, 0.5 + dx * 0.6, 0.86, 0.5 + dx, 0.92) - line * 0.35;
        c = over(c, OUTLINE, fill(d));
    }

    // the bulb: cream, shaded towards the lower right
    let b = bulb(x, y);
    let light = ((x - 0.40).powi(2) + (y - 0.48).powi(2)).sqrt() / 0.42;
    c = over(c, mix(BULB_LIGHT, BULB_SHADE, light * light), fill(b));
    // clove lines: arcs from the tip to the root, bowing outwards around the
    // face (circles through both points, bulging by `sag`), inside the bulb only
    if b < 0.0 {
        let (top, bottom, sag) = (0.30_f32, 0.875_f32, 0.17_f32);
        let half = (bottom - top) / 2.0;
        let r = (half * half + sag * sag) / (2.0 * sag);
        for side in [-1.0_f32, 1.0] {
            // centre on the far side, so the arc bows towards `side`
            let arc = circle(x, y, 0.5 - side * (r - sag), top + half, r);
            let on_side = if (x - 0.5) * side > 0.0 { 1.0 } else { 0.0 };
            let fade = ((-b) / 0.03).clamp(0.0, 1.0);
            c = over(c, CLOVE_LINE, stroke(arc, line * 0.55) * fade * on_side);
        }
    }
    c = over(c, OUTLINE, stroke(b, line));

    // face: rosy cheeks, eyes with a highlight, and a smile
    for cx in [0.355_f32, 0.645] {
        c = over(c, CHEEK, fill(ellipse(x, y, cx, 0.665, 0.045, 0.027)));
    }
    for cx in [0.415_f32, 0.585] {
        c = over(c, OUTLINE, fill(ellipse(x, y, cx, 0.6, 0.026, 0.032)));
        c = over(c, WHITE, fill(circle(x, y, cx + 0.009, 0.588, 0.009)));
    }
    // smile: an arc of the lower part of a circle, with round ends
    c = over(c, OUTLINE, fill(arc(x, y, 0.5, 0.615, 0.068, 0.45, 2.69) - line * 0.45));
    c
}

/// Distance to the arc of circle (cx, cy, r) between angles `from`..`to`
/// (radians, y pointing down), with round ends
fn arc(x: f32, y: f32, cx: f32, cy: f32, r: f32, from: f32, to: f32) -> f32 {
    let angle = (y - cy).atan2(x - cx);
    if (from..=to).contains(&angle) {
        (circle(x, y, cx, cy, r)).abs()
    } else {
        let end = |a: f32| ((x - cx - r * a.cos()).powi(2) + (y - cy - r * a.sin()).powi(2)).sqrt();
        end(from).min(end(to))
    }
}

fn render(size: u32) -> RgbaImage {
    // 4x4 supersampling on top of the analytic edges
    const SS: u32 = 4;
    let px = 1.0 / size as f32;
    RgbaImage::from_fn(size, size, |ix, iy| {
        let mut acc = [0.0_f32; 4];
        for sy in 0..SS {
            for sx in 0..SS {
                let x = (ix as f32 + (sx as f32 + 0.5) / SS as f32) * px;
                let y = (iy as f32 + (sy as f32 + 0.5) / SS as f32) * px;
                let c = shade(x, y, px);
                // premultiplied accumulation
                for i in 0..3 {
                    acc[i] += c[i] * c[3];
                }
                acc[3] += c[3];
            }
        }
        let n = (SS * SS) as f32;
        let a = acc[3] / n;
        let ch = |v: f32| ((if a > 0.0 { v / acc[3] } else { 0.0 }) * 255.0).round().clamp(0.0, 255.0) as u8;
        Rgba([ch(acc[0]), ch(acc[1]), ch(acc[2]), (a * 255.0).round() as u8])
    })
}

fn png(img: &RgbaImage) -> Vec<u8> {
    let mut out = Vec::new();
    image::codecs::png::PngEncoder::new(&mut out)
        .write_image(img.as_raw(), img.width(), img.height(), image::ExtendedColorType::Rgba8)
        .unwrap();
    out
}

/// An .ico holding PNG images (supported since Windows Vista)
fn ico(sizes: &[u32]) -> Vec<u8> {
    let images: Vec<Vec<u8>> = sizes.iter().map(|&s| png(&render(s))).collect();
    let mut out = vec![0, 0, 1, 0];
    out.extend((sizes.len() as u16).to_le_bytes());
    let mut offset = 6 + 16 * sizes.len() as u32;
    for (&s, data) in sizes.iter().zip(&images) {
        let dim = if s >= 256 { 0 } else { s as u8 };
        out.extend([dim, dim, 0, 0]);
        out.extend(1u16.to_le_bytes()); // colour planes
        out.extend(32u16.to_le_bytes()); // bits per pixel
        out.extend((data.len() as u32).to_le_bytes());
        out.extend(offset.to_le_bytes());
        offset += data.len() as u32;
    }
    for data in images {
        out.extend(data);
    }
    out
}

fn main() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("assets");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("garlic.png"), png(&render(256))).unwrap();
    std::fs::write(dir.join("garlic.ico"), ico(&[16, 20, 24, 32, 40, 48, 64, 128, 256])).unwrap();
    // a large copy to look at
    std::fs::write(dir.join("garlic-512.png"), png(&render(512))).unwrap();
    println!("wrote {}", dir.display());
}
