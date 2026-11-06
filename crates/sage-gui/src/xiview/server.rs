//! The local xiVIEW server: 127.0.0.1 only, on a free port, for as long as
//! sage-gui runs. It serves
//!
//! * the xiVIEW web app (the `xiview` folder next to sage-gui.exe),
//! * `/data/get_xiview_*`: the PRIDE crosslinking API's answers for each
//!   opened results folder ([`super::data`]), chosen by the `file` URL
//!   parameter,
//! * `/data/get_peaklist` and `/data/get_annotated_peaklist`: a match's
//!   spectrum, read as Sage read it, plain or annotated
//!   ([`super::annotate`]),
//! * `/xiAnnotator/annotate/FULL`: xiSPEC's re-annotation (its settings)
//!   with the peaks in the request.

use super::annotate::{self, Known};
use super::data::{request_signature, Fdr, Upload};
use crate::viewer::model::{Kind, Peaks, SpectrumFile};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, RwLock};

/// The files of the xiVIEW web app: the `xiview` folder next to sage-gui.exe
/// (or `XNOVA_XIVIEW_DIR`)
pub fn app_folder() -> Option<PathBuf> {
    let candidates = std::env::var_os("XNOVA_XIVIEW_DIR")
        .map(PathBuf::from)
        .into_iter()
        .chain(std::env::current_exe().ok().and_then(|e| e.parent().map(|p| p.join("xiview"))));
    candidates.into_iter().find(|p| p.join("network.html").is_file() && p.join("xiview.js").is_file())
}

/// (upload, spectrum file index, spectrum id) → peaks
type SpectrumKey = (u32, usize, String);

struct State {
    root: PathBuf,
    uploads: RwLock<Vec<Arc<Upload>>>,
    spectra: Mutex<HashMap<SpectrumKey, Arc<Peaks>>>,
    /// Spectrum files read whole (mzML, MGF): one read at a time, once
    reading: Mutex<HashSet<(u32, usize)>>,
    /// The match whose spectrum was opened last: xiSPEC's settings re-annotate
    /// it without saying which match it is
    last: Mutex<Option<(u32, String)>>,
    /// timsTOF runs kept open (opening one takes about a second), read one
    /// request at a time each
    opened: Mutex<HashMap<(u32, usize), Arc<Mutex<SpectrumFile>>>>,
}

pub struct Server {
    port: u16,
    state: Arc<State>,
    http: Arc<tiny_http::Server>,
    thread: Option<std::thread::JoinHandle<()>>,
}

impl Server {
    /// Start serving the xiVIEW files in `root` on 127.0.0.1
    pub fn start(root: PathBuf) -> Result<Server, String> {
        let http = Arc::new(tiny_http::Server::http("127.0.0.1:0").map_err(|e| format!("local xiVIEW server: {e}"))?);
        let port = http.server_addr().to_ip().map(|a| a.port()).ok_or("local xiVIEW server: no port")?;
        let state = Arc::new(State {
            root,
            uploads: RwLock::new(Vec::new()),
            spectra: Mutex::new(HashMap::new()),
            reading: Mutex::new(HashSet::new()),
            last: Mutex::new(None),
            opened: Mutex::new(HashMap::new()),
        });
        let (h, s) = (http.clone(), state.clone());
        let thread = std::thread::Builder::new()
            .name("xiview-server".into())
            .spawn(move || {
                for request in h.incoming_requests() {
                    let s = s.clone();
                    // xiVIEW asks for its data in parallel
                    let _ = std::thread::Builder::new().name("xiview-request".into()).spawn(move || handle(&s, request));
                }
            })
            .map_err(|e| e.to_string())?;
        Ok(Server { port, state, http, thread: Some(thread) })
    }

    /// The page of an opened results folder
    pub fn url(&self, name: &str) -> String {
        format!("http://127.0.0.1:{}/network.html?project=xnova&file={}", self.port, encode(name))
    }

    /// Is this results folder (with these matches) open?
    pub fn has(&self, folder: &Path, fdr: Fdr) -> Option<String> {
        let uploads = self.state.uploads.read().ok()?;
        uploads.iter().find(|u| u.folder == folder && u.fdr == fdr).map(|u| u.name.clone())
    }

    /// Load a results folder (slow: reads its tables and FASTA) and serve it;
    /// returns its xiVIEW name (`file` parameter)
    pub fn add(&self, folder: &Path, fdr: Fdr) -> Result<String, String> {
        let id = {
            let uploads = self.state.uploads.read().map_err(|e| e.to_string())?;
            uploads.iter().map(|u| u.id).max().unwrap_or(0) + 1
        };
        let mut upload = Upload::load(folder, id, fdr)?;
        if upload.served.is_empty() {
            return Err(format!(
                "{}: no crosslinks pass XNova's {} FDR (or there is no crosslinks.sage.tsv)",
                folder.display(),
                fdr.percent()
            ));
        }
        // one name per open folder and choice of matches
        let mut uploads = self.state.uploads.write().map_err(|e| e.to_string())?;
        uploads.retain(|u| !(u.folder == folder && u.fdr == fdr));
        let base = if fdr == Fdr::default() { upload.name.clone() } else { format!("{} ({} FDR)", upload.name, fdr.percent()) };
        let mut name = base.clone();
        let mut n = 2;
        while uploads.iter().any(|u| u.name == name) {
            name = format!("{base} ({n})");
            n += 1;
        }
        upload.name = name.clone();
        let upload = Arc::new(upload);
        uploads.push(upload.clone());
        drop(uploads);
        // in the background: read the spectra of mzML / MGF files (once), and
        // open the timsTOF runs, so the first spectrum of each comes quickly
        let state = self.state.clone();
        let _ = std::thread::Builder::new().name("xiview-spectra".into()).spawn(move || {
            for file in 0..upload.files.len() {
                if upload.results.read_whole(&upload.files[file]) {
                    let _ = read_whole_file(&state, &upload, file);
                } else {
                    let _ = opened_run(&state, &upload, file);
                }
            }
        });
        Ok(name)
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        self.http.unblock();
        if let Some(t) = self.thread.take() {
            let _ = t.join();
        }
    }
}

/// `%`-encoding for a URL parameter
pub fn encode(text: &str) -> String {
    let mut out = String::new();
    for b in text.bytes() {
        if b.is_ascii_alphanumeric() || b"-_.~".contains(&b) {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

fn decode(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        let hex = |j: usize| bytes.get(j).and_then(|&c| (c as char).to_digit(16));
        match bytes[i] {
            b'+' => out.push(b' '),
            b'%' => match (hex(i + 1), hex(i + 2)) {
                (Some(h), Some(l)) => {
                    out.push((h * 16 + l) as u8);
                    i += 2;
                }
                _ => out.push(b'%'),
            },
            b => out.push(b),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// Path and query parameters of a request URL
fn split_url(url: &str) -> (String, HashMap<String, String>) {
    let (path, query) = url.split_once('?').unwrap_or((url, ""));
    let params = query
        .split('&')
        .filter(|p| !p.is_empty())
        .map(|p| {
            let (k, v) = p.split_once('=').unwrap_or((p, ""));
            (decode(k), decode(v))
        })
        .collect();
    (decode(path), params)
}

fn content_type(path: &Path) -> &'static str {
    match path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase().as_str() {
        "html" | "htm" => "text/html; charset=utf-8",
        "js" => "text/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "json" => "application/json",
        "svg" => "image/svg+xml",
        "png" => "image/png",
        "gif" => "image/gif",
        "jpg" | "jpeg" => "image/jpeg",
        "woff" => "font/woff",
        "woff2" => "font/woff2",
        "ttf" => "font/ttf",
        "eot" => "application/vnd.ms-fontobject",
        "mp4" => "video/mp4",
        "webm" => "video/webm",
        "txt" | "obo" => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}

fn header(name: &str, value: &str) -> tiny_http::Header {
    tiny_http::Header::from_bytes(name.as_bytes(), value.as_bytes()).expect("valid header")
}

fn respond(request: tiny_http::Request, status: u16, kind: &str, body: Vec<u8>) {
    let response = tiny_http::Response::from_data(body)
        .with_status_code(status)
        .with_header(header("Content-Type", kind))
        .with_header(header("Cache-Control", "no-store"));
    let _ = request.respond(response);
}

fn respond_json(request: tiny_http::Request, value: &Value) {
    respond(request, 200, "application/json", serde_json::to_vec(value).unwrap_or_default());
}

fn not_found(request: tiny_http::Request, what: &str) {
    respond(request, 404, "application/json", serde_json::to_vec(&json!({"detail": what})).unwrap_or_default());
}

fn handle(state: &State, mut request: tiny_http::Request) {
    let (path, params) = split_url(request.url());
    let post = *request.method() == tiny_http::Method::Post;
    let mut body = String::new();
    if post {
        let _ = request.as_reader().take(64 << 20).read_to_string(&mut body);
    }
    let path = path.trim_start_matches('/').to_string();

    if let Some(endpoint) = path.strip_prefix("data/") {
        return data(state, request, endpoint, &params, &body);
    }
    if path.starts_with("xiAnnotator/") {
        // xiSPEC's settings: annotate the peaks it sends
        let Ok(req) = serde_json::from_str::<Value>(&body) else { return not_found(request, "no annotation request") };
        let peaks = annotate::request_peaks(&req);
        let uploads = state.uploads.read().map(|u| u.clone()).unwrap_or_default();
        let id = req["annotation"]["requestID"].as_str().map(String::from).unwrap_or_else(|| req["annotation"]["requestID"].to_string());
        // the match it is for: named, else the one opened last
        let last = state.last.lock().ok().and_then(|l| l.clone());
        let found = uploads.iter().find_map(|u| u.served(&id).map(|(m, s)| (u.clone(), m.clone(), s.clone()))).or_else(|| {
            let (upload, mid) = last?;
            let u = uploads.iter().find(|u| u.id == upload)?;
            u.served(&mid).map(|(m, s)| (u.clone(), m.clone(), s.clone()))
        });
        let result = match &found {
            Some((u, m, s)) if unchanged(&req, &s.signatures) => {
                annotate::annotate(&req, &peaks, Some(Known { results: &u.results, m, layout: s.layout.as_ref() }))
            }
            _ => annotate::annotate(&req, &peaks, None),
        };
        return match result {
            Ok(v) => respond_json(request, &v),
            Err(e) => respond(request, 400, "application/json", serde_json::to_vec(&json!({"error": e})).unwrap_or_default()),
        };
    }

    // the web app's files; `docs/` are xiVIEW's help pages
    let relative = match path.as_str() {
        "" => "network.html".to_string(),
        p => p.strip_prefix("docs/").map_or_else(|| p.to_string(), |d| format!("xidocs/{d}")),
    };
    if relative.split(['/', '\\']).any(|part| part == ".." || part.contains(':')) {
        return not_found(request, "not found");
    }
    let file = state.root.join(&relative);
    match std::fs::read(&file) {
        Ok(bytes) => respond(request, 200, content_type(&file), bytes),
        Err(_) => not_found(request, "not found"),
    }
}

/// Did xiSPEC send back the peptides as served (not edited in its settings)?
fn unchanged(request: &Value, signatures: &[Vec<String>]) -> bool {
    let peptides: Vec<String> = request["Peptides"].as_array().into_iter().flatten().map(request_signature).collect();
    signatures.iter().any(|s| *s == peptides)
}

/// The uploads a request is for: its `file` parameter, else all
fn chosen(state: &State, params: &HashMap<String, String>) -> Vec<Arc<Upload>> {
    let uploads = state.uploads.read().map(|u| u.clone()).unwrap_or_default();
    match params.get("file") {
        Some(f) if !f.is_empty() => uploads.into_iter().filter(|u| &u.name == f).collect(),
        _ => uploads,
    }
}

fn data(state: &State, request: tiny_http::Request, endpoint: &str, params: &HashMap<String, String>, body: &str) {
    match endpoint {
        "get_peaklist" | "get_annotated_peaklist" => {
            let upload_id: u32 = params.get("upload_id").and_then(|u| u.parse().ok()).unwrap_or(0);
            let file: usize = params.get("sd_ref").and_then(|s| s.parse().ok()).unwrap_or(0);
            let Some(id) = params.get("id") else { return not_found(request, "no spectrum id") };
            let upload = state.uploads.read().ok().and_then(|u| u.iter().find(|u| u.id == upload_id).cloned());
            let Some(upload) = upload else { return not_found(request, "no such results folder") };
            let peaks = match spectrum(state, &upload, file, id) {
                Ok(p) => p,
                Err(e) => return not_found(request, &e),
            };
            if endpoint == "get_peaklist" {
                return respond_json(request, &json!({"mz": peaks.mz, "intensity": peaks.intensity}));
            }
            let Ok(req) = serde_json::from_str::<Value>(body) else { return not_found(request, "no annotation request") };
            let match_id = req["annotation"]["requestID"].as_str().unwrap_or("").to_string();
            if let Ok(mut last) = state.last.lock() {
                *last = Some((upload.id, match_id.clone()));
            }
            // xiVIEW asks for its own match (as served): XNova's peptidoform
            let known = upload.served(&match_id).filter(|(_, s)| unchanged(&req, &s.signatures));
            let result = match known {
                Some((m, s)) => {
                    annotate::annotate(&req, &peaks, Some(Known { results: &upload.results, m, layout: s.layout.as_ref() }))
                }
                None => annotate::annotate(&req, &peaks, None),
            };
            match result {
                Ok(v) => respond_json(request, &v),
                Err(e) => respond(request, 400, "application/json", serde_json::to_vec(&json!({"error": e})).unwrap_or_default()),
            }
        }
        "get_matches_by_multiple_spectra_id" => respond_json(request, &json!([])),
        "get_datasets" => {
            let uploads = state.uploads.read().map(|u| u.clone()).unwrap_or_default();
            let rows: Vec<Value> =
                uploads.iter().map(|u| json!({"project_id": "xnova", "identification_file_name": u.name})).collect();
            respond_json(request, &Value::Array(rows))
        }
        name => {
            let uploads = chosen(state, params);
            if uploads.is_empty() {
                return not_found(request, "No upload found for the given project/file");
            }
            // one upload: its answer as is; several: the arrays joined
            if uploads.len() == 1 {
                return match uploads[0].endpoint(name) {
                    Some(bytes) => respond(request, 200, "application/json", bytes.to_vec()),
                    None => not_found(request, "not found"),
                };
            }
            let mut rows = Vec::new();
            for u in &uploads {
                if let Some(Value::Array(a)) = u.endpoint(name).and_then(|b| serde_json::from_slice(b).ok()) {
                    rows.extend(a);
                }
            }
            respond_json(request, &Value::Array(rows))
        }
    }
}

/// The peaks of a spectrum: from the cache, else read (timsTOF: this one
/// spectrum; mzML / MGF: every served spectrum of the file at once)
fn spectrum(state: &State, upload: &Upload, file: usize, id: &str) -> Result<Arc<Peaks>, String> {
    let key = (upload.id, file, id.to_string());
    if let Some(p) = state.spectra.lock().map_err(|e| e.to_string())?.get(&key) {
        return Ok(p.clone());
    }
    let name = upload.files.get(file).ok_or("no such spectrum file")?;
    if upload.results.read_whole(name) {
        read_whole_file(state, upload, file)?;
    } else {
        let run = opened_run(state, upload, file)?;
        let ids = HashSet::from([id.to_string()]);
        let spectra = run.lock().map_err(|e| e.to_string())?.read(&ids)?;
        let mut cache = state.spectra.lock().map_err(|e| e.to_string())?;
        for (sid, peaks) in spectra {
            cache.insert((upload.id, file, sid), Arc::new(peaks));
        }
    }
    state.spectra.lock().map_err(|e| e.to_string())?.get(&key).cloned().ok_or_else(|| format!("spectrum {id} not found in {name}"))
}

/// A timsTOF run, opened once and kept open
fn opened_run(state: &State, upload: &Upload, file: usize) -> Result<Arc<Mutex<SpectrumFile>>, String> {
    let key = (upload.id, file);
    if let Some(run) = state.opened.lock().map_err(|e| e.to_string())?.get(&key) {
        return Ok(run.clone());
    }
    // (opened without holding the lock: other runs and spectra are not held up)
    let name = upload.files.get(file).ok_or("no such spectrum file")?;
    let run = Arc::new(Mutex::new(upload.results.open_file(name)?));
    Ok(state.opened.lock().map_err(|e| e.to_string())?.entry(key).or_insert(run).clone())
}

/// Read the served spectra of an mzML / MGF file (once)
fn read_whole_file(state: &State, upload: &Upload, file: usize) -> Result<(), String> {
    // one reader per file; others wait for it
    loop {
        {
            let mut reading = state.reading.lock().map_err(|e| e.to_string())?;
            if !reading.contains(&(upload.id, file)) {
                reading.insert((upload.id, file));
                break;
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    let done = |state: &State| {
        if let Ok(mut r) = state.reading.lock() {
            r.remove(&(upload.id, file));
        }
    };
    let name = &upload.files[file];
    let wanted: HashSet<String> = upload
        .served
        .values()
        .filter_map(|s| upload.results.matches.get(s.index))
        .filter(|m| &m.file == name && m.kind != Kind::Linear)
        .map(|m| m.scan.clone())
        .collect();
    let missing: HashSet<String> = {
        let cache = state.spectra.lock().map_err(|e| e.to_string())?;
        wanted.into_iter().filter(|id| !cache.contains_key(&(upload.id, file, id.clone()))).collect()
    };
    if missing.is_empty() {
        done(state);
        return Ok(());
    }
    let result = upload.results.open_file(name).and_then(|source| source.read(&missing));
    match result {
        Ok(spectra) => {
            let mut cache = state.spectra.lock().map_err(|e| e.to_string())?;
            for (sid, peaks) in spectra {
                cache.insert((upload.id, file, sid), Arc::new(peaks));
            }
            drop(cache);
            done(state);
            Ok(())
        }
        Err(e) => {
            done(state);
            Err(e)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Serve a results folder to try xiVIEW in a browser:
    /// `XNOVA_XIVIEW_DIR=<xiview folder> XNOVA_XIVIEW_RESULTS=<results folder>
    /// cargo test -p sage-gui serve_a_results_folder -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn serve_a_results_folder() {
        let root = app_folder().expect("XNOVA_XIVIEW_DIR: the xiview folder");
        let folder = PathBuf::from(std::env::var("XNOVA_XIVIEW_RESULTS").expect("XNOVA_XIVIEW_RESULTS"));
        let minutes: u64 = std::env::var("XNOVA_XIVIEW_MINUTES").ok().and_then(|m| m.parse().ok()).unwrap_or(30);
        let server = Server::start(root).unwrap();
        for fdr in [Fdr::default(), Fdr(0.05)] {
            let start = std::time::Instant::now();
            let name = server.add(&folder, fdr).unwrap();
            println!("URL {} (loaded in {:.1} s)", server.url(&name), start.elapsed().as_secs_f32());
        }
        std::thread::sleep(std::time::Duration::from_secs(minutes * 60));
    }

    /// xiSPEC's annotation against XNova's viewer on real spectra:
    /// `XNOVA_XIVIEW_RESULTS=<results folder> cargo test -p sage-gui
    /// annotations_agree_with_the_viewer -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn annotations_agree_with_the_viewer() {
        use crate::viewer::model::{annotate as viewer_annotate, Model};
        let folder = PathBuf::from(std::env::var("XNOVA_XIVIEW_RESULTS").expect("XNOVA_XIVIEW_RESULTS"));
        let upload = Upload::load(&folder, 1, Fdr(1.0)).unwrap();
        let mut ids: Vec<&String> = upload.served.keys().collect();
        ids.sort();
        // a spread of matches, every kind
        let step = (ids.len() / 40).max(1);
        let mut picked: Vec<&String> = ids.iter().step_by(step).copied().collect();
        for kind in [Kind::Loop, Kind::Multi] {
            picked.extend(ids.iter().filter(|id| upload.served(id).is_some_and(|(m, _)| m.kind == kind)).take(4).copied());
        }
        let state = State {
            root: PathBuf::new(),
            uploads: RwLock::new(Vec::new()),
            spectra: Mutex::new(HashMap::new()),
            reading: Mutex::new(HashSet::new()),
            last: Mutex::new(None),
            opened: Mutex::new(HashMap::new()),
        };
        let (mut same, mut only_viewer, mut only_xi, mut seconds) = (0usize, 0usize, 0usize, 0f64);
        let mut by_kind: HashMap<String, usize> = HashMap::new();
        for id in picked {
            let (m, s) = upload.served(id).unwrap();
            let file = upload.files.iter().position(|f| f == &m.file).unwrap();
            let peaks = spectrum(&state, &upload, file, &m.scan).unwrap();
            let mut request = json!({
                "Peptides": [], "LinkSite": [], "peaks": [],
                "annotation": {
                    "requestID": id, "precursorCharge": m.charge, "modifications": [],
                    "fragmentTolerance": {"tolerance": 20.0, "unit": "ppm"},
                    "ions": [{"type": "PeptideIon"}, {"type": "AIon"}, {"type": "BIon"}, {"type": "YIon"}],
                    "crosslinker": {"modMass": 0}
                }
            });
            request["Peptides"] = json!(s.signatures[0].iter().map(|p| json!({"sequence": p.chars().filter(|c| c.is_ascii_uppercase()).map(|c| json!({"aminoAcid": c.to_string(), "Modification": ""})).collect::<Vec<_>>()})).collect::<Vec<_>>());
            let start = std::time::Instant::now();
            let xi = annotate::annotate(&request, &peaks, Some(Known { results: &upload.results, m, layout: s.layout.as_ref() })).unwrap();
            seconds += start.elapsed().as_secs_f64();
            // peaks explained in xiSPEC: those of a cluster a fragment matched
            let clusters: HashSet<u64> = xi["fragments"].as_array().unwrap().iter().flat_map(|f| f["clusterIds"].as_array().unwrap().iter().filter_map(Value::as_u64)).collect();
            let xi_mz: HashSet<u64> = xi["peaks"].as_array().unwrap().iter()
                .filter(|p| p["clusterIds"].as_array().unwrap().iter().any(|c| c.as_u64().is_some_and(|c| clusters.contains(&c))))
                .map(|p| (p["mz"].as_f64().unwrap() * 1e5) as u64).collect();
            let viewer = viewer_annotate(&peaks, &upload.results.proforma(m), m.charge, Model::CidHcd, upload.results.fragment_tol, upload.results.link(m).as_ref()).unwrap();
            // the viewer's b/y/a/precursor and stub ions (not immonium or internal fragments)
            use crate::viewer::model::Series;
            let viewer_mz: HashSet<u64> = viewer.peaks.iter()
                .filter(|p| p.labels.iter().any(|l| matches!(l.series, Series::A | Series::B | Series::Y | Series::Precursor | Series::Stub(_))))
                .map(|p| (p.mz * 1e5) as u64).collect();
            let shown = if s.layout.is_some() { "multivalent" } else { match m.kind { Kind::Loop => "loop", _ => "crosslink" } };
            *by_kind.entry(shown.to_string()).or_default() += 1;
            {
                same += xi_mz.intersection(&viewer_mz).count();
                only_viewer += viewer_mz.difference(&xi_mz).count();
                only_xi += xi_mz.difference(&viewer_mz).count();
            }
            println!("{id} {:?} z{} {}: xiSPEC {} peaks, viewer {} peaks, shared {}", m.kind, m.charge, m.title(), xi_mz.len(), viewer_mz.len(), xi_mz.intersection(&viewer_mz).count());
        }
        println!("{by_kind:?}");
        println!("crosslinks, loop-links and multivalent species: shared {same}, viewer only {only_viewer}, xiSPEC only {only_xi}; annotation {:.1} ms each", seconds * 1000.0 / by_kind.values().sum::<usize>() as f64);
        assert!(same > 0 && only_viewer * 10 <= same, "xiSPEC misses peaks the viewer explains");
    }

    #[test]
    fn urls() {
        assert_eq!(encode("PolII DSSO/run 1"), "PolII%20DSSO%2Frun%201");
        assert_eq!(decode("PolII%20DSSO%2Frun+1"), "PolII DSSO/run 1");
        let (path, params) = split_url("/data/get_xiview_matches?project=xnova&file=PolII%20DSSO");
        assert_eq!(path, "/data/get_xiview_matches");
        assert_eq!(params.get("file").map(String::as_str), Some("PolII DSSO"));
        assert_eq!(decode("100%"), "100%");
    }
}
