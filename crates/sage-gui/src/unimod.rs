//! Look up modifications in UniMod (https://www.unimod.org) by accession number.
//!
//! The record pages behind unimod.org's modification list are read from
//! UniMod's own XML export, downloaded once and cached next to the GUI state.
//! A missing accession triggers one fresh download, for records added since.

use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver};

pub const XML_URL: &str = "https://www.unimod.org/xml/unimod.xml";

/// One allowed position of a modification
#[derive(Clone, Debug, PartialEq)]
pub struct Specificity {
    /// Residue letter, or `N-term` / `C-term`
    pub site: String,
    /// `Anywhere`, `Any N-term`, `Any C-term`, `Protein N-term`, `Protein C-term`
    pub position: String,
    pub classification: String,
    /// UniMod hides rare or unverified specificities from its default view
    pub hidden: bool,
    /// The same position as a Sage modification target, if Sage can express it
    pub target: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub accession: u32,
    /// PSI-MS name, e.g. `Oxidation`
    pub title: String,
    pub full_name: String,
    pub mono_mass: f64,
    pub composition: String,
    pub specificities: Vec<Specificity>,
}

impl Entry {
    /// ProForma name written into the Sage configuration
    pub fn proforma(&self) -> String {
        format!("UNIMOD:{}", self.accession)
    }
}

/// Accept `35`, `#35`, `UNIMOD:35`, `unimod 35`
pub fn parse_accession(text: &str) -> Option<u32> {
    let t = text.trim().to_ascii_lowercase();
    let t = t.strip_prefix("unimod").unwrap_or(&t);
    let t = t.trim_start_matches([':', '#', ' ']);
    t.parse().ok().filter(|&n| n > 0)
}

/// UniMod site + position as a Sage target: `K`, `^`, `^Q`, `$`, `[`, `[M`, `]`
fn sage_target(site: &str, position: &str) -> Option<String> {
    let terminus = match position {
        "Anywhere" => "",
        "Any N-term" => "^",
        "Any C-term" => "$",
        "Protein N-term" => "[",
        "Protein C-term" => "]",
        _ => return None,
    };
    match site {
        "N-term" | "C-term" if !terminus.is_empty() => Some(terminus.into()),
        aa if aa.len() == 1 && crate::settings::VALID_AA.contains(aa) => Some(format!("{terminus}{aa}")),
        _ => None,
    }
}

/// Find one record in the UniMod XML. `Ok(None)` when it is not there.
pub fn find(xml: &str, accession: u32) -> Result<Option<Entry>, String> {
    let doc = roxmltree::Document::parse(xml).map_err(|e| format!("UniMod XML is unreadable: {e}"))?;
    let id = accession.to_string();
    let Some(node) = doc
        .descendants()
        .find(|n| n.has_tag_name("mod") && n.attribute("record_id") == Some(id.as_str()))
    else {
        return Ok(None);
    };
    let attr = |n: roxmltree::Node, k: &str| n.attribute(k).unwrap_or("").to_string();
    let delta = node
        .children()
        .find(|n| n.has_tag_name("delta"))
        .ok_or(format!("UniMod record {accession} has no mass"))?;
    let mono_mass = delta
        .attribute("mono_mass")
        .and_then(|m| m.parse().ok())
        .ok_or(format!("UniMod record {accession} has no monoisotopic mass"))?;
    let mut specificities: Vec<Specificity> = node
        .children()
        .filter(|n| n.has_tag_name("specificity"))
        .map(|s| {
            let (site, position) = (attr(s, "site"), attr(s, "position"));
            Specificity {
                target: sage_target(&site, &position),
                site,
                position,
                classification: attr(s, "classification"),
                hidden: s.attribute("hidden") == Some("1"),
            }
        })
        .collect();
    // Shown specificities first, then by site, as on unimod.org
    specificities.sort_by(|a, b| (a.hidden, &a.site, &a.position).cmp(&(b.hidden, &b.site, &b.position)));
    Ok(Some(Entry {
        accession,
        title: attr(node, "title"),
        full_name: attr(node, "full_name"),
        mono_mass,
        composition: attr(delta, "composition"),
        specificities,
    }))
}

fn cache_path() -> Option<PathBuf> {
    let base = std::env::var_os("APPDATA").or_else(|| std::env::var_os("HOME"))?;
    Some(PathBuf::from(base).join("sage-gui").join("unimod.xml"))
}

fn download() -> Result<String, String> {
    let xml = ureq::get(XML_URL)
        .timeout(std::time::Duration::from_secs(60))
        .call()
        .map_err(|e| format!("Could not download UniMod ({XML_URL}): {e}"))?
        .into_string()
        .map_err(|e| format!("Could not read the UniMod download: {e}"))?;
    if !xml.contains("<umod:unimod") {
        return Err(format!("{XML_URL} did not return the UniMod database"));
    }
    if let Some(path) = cache_path() {
        if let Some(dir) = path.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let _ = std::fs::write(path, &xml);
    }
    Ok(xml)
}

/// Look up an accession: cached copy first, then a fresh download
pub fn lookup(accession: u32) -> Result<Entry, String> {
    if let Some(xml) = cache_path().and_then(|p| std::fs::read_to_string(p).ok()) {
        if let Ok(Some(entry)) = find(&xml, accession) {
            return Ok(entry);
        }
    }
    find(&download()?, accession)?.ok_or(format!("UniMod has no record with accession #{accession}"))
}

/// State of the "Add from UniMod" panel on the Modifications tab
#[derive(Default)]
pub struct Lookup {
    pub input: String,
    pub pending: Option<Receiver<Result<Entry, String>>>,
    pub result: Option<Result<Entry, String>>,
}

impl Lookup {
    /// Fetch on a background thread so the window stays responsive
    pub fn start(&mut self, ctx: &eframe::egui::Context) {
        let Some(accession) = parse_accession(&self.input) else {
            self.result = Some(Err(format!(
                "`{}` is not a UniMod accession number; enter e.g. 35 or UNIMOD:35",
                self.input.trim()
            )));
            return;
        };
        let (tx, rx) = channel();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let _ = tx.send(lookup(accession));
            ctx.request_repaint();
        });
        self.pending = Some(rx);
        self.result = None;
    }

    pub fn busy(&mut self) -> bool {
        if let Some(rx) = &self.pending {
            if let Ok(result) = rx.try_recv() {
                self.result = Some(result);
                self.pending = None;
            }
        }
        self.pending.is_some()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const XML: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">
  <umod:modifications>
    <umod:mod title="Oxidation" full_name="Oxidation or Hydroxylation" record_id="35">
      <umod:specificity hidden="1" site="C" position="Anywhere" classification="Post-translational" spec_group="9"/>
      <umod:specificity hidden="0" site="M" position="Anywhere" classification="Artefact" spec_group="8"/>
      <umod:specificity hidden="1" site="G" position="Any C-term" classification="Pre-translational" spec_group="11"/>
      <umod:delta mono_mass="15.994915" avge_mass="15.9994" composition="O">
        <umod:element symbol="O" number="1"/>
      </umod:delta>
    </umod:mod>
    <umod:mod title="Acetyl" full_name="Acetylation" record_id="1">
      <umod:specificity hidden="0" site="N-term" position="Protein N-term" classification="Post-translational" spec_group="5"/>
      <umod:specificity hidden="0" site="N-term" position="Any N-term" classification="Chemical derivative" spec_group="2"/>
      <umod:delta mono_mass="42.010565" avge_mass="42.0367" composition="H(2) C(2) O"/>
    </umod:mod>
  </umod:modifications>
</umod:unimod>"#;

    #[test]
    fn accessions() {
        for (text, n) in [("35", Some(35)), (" #35 ", Some(35)), ("UNIMOD:35", Some(35)), ("unimod 4", Some(4))] {
            assert_eq!(parse_accession(text), n, "{text}");
        }
        for bad in ["", "0", "ox", "UNIMOD:", "-3"] {
            assert_eq!(parse_accession(bad), None, "{bad}");
        }
    }

    #[test]
    fn records_map_to_sage_targets() {
        let ox = find(XML, 35).unwrap().unwrap();
        assert_eq!((ox.title.as_str(), ox.mono_mass, ox.composition.as_str()), ("Oxidation", 15.994915, "O"));
        assert_eq!(ox.proforma(), "UNIMOD:35");
        // shown specificity first
        let targets: Vec<_> = ox.specificities.iter().map(|s| (s.target.clone().unwrap(), s.hidden)).collect();
        assert_eq!(targets, [("M".into(), false), ("C".into(), true), ("$G".into(), true)]);

        let ac = find(XML, 1).unwrap().unwrap();
        let targets: Vec<_> = ac.specificities.iter().filter_map(|s| s.target.clone()).collect();
        assert_eq!(targets, ["^", "["]);

        assert_eq!(find(XML, 999).unwrap(), None);
    }

    /// Against the live unimod.org export; skipped when offline
    #[test]
    fn live_unimod() {
        let Ok(xml) = download() else {
            eprintln!("skipped: unimod.org not reachable");
            return;
        };
        let cam = find(&xml, 4).unwrap().unwrap();
        assert_eq!(cam.title, "Carbamidomethyl");
        assert!((cam.mono_mass - 57.021464).abs() < 1e-6);
        assert!(cam.specificities.iter().any(|s| s.target.as_deref() == Some("C") && !s.hidden));

        // common modifications reach the targets the GUI presets use
        let shown = |acc: u32| -> Vec<String> {
            let e = find(&xml, acc).unwrap().unwrap();
            e.specificities.into_iter().filter(|s| !s.hidden).filter_map(|s| s.target).collect()
        };
        for (acc, targets) in [(21, &["S", "T", "Y"][..]), (1, &["[", "K"]), (737, &["K", "^"]), (28, &["^Q"])] {
            let got = shown(acc);
            for t in targets {
                assert!(got.contains(&t.to_string()), "UNIMOD:{acc}: {t} not in {got:?}");
            }
        }
    }
}
