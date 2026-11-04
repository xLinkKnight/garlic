//! Reads the spectra of matches in the background, so the viewer never waits
//! for a whole run: the match picked first, then the rows in view of runs
//! already open. Bruker runs are read spectrum by spectrum, by id, exactly as
//! the search read them (a run is opened, in about a second, when one of its
//! matches is picked). Files without an index of spectra (mzML, MGF) are read
//! whole, once, keeping the spectra of every shown match, one after another.

use super::model::{Peaks, Results, SpectrumFile, SpectrumKey};
use eframe::egui;
use std::collections::{HashMap, HashSet};
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::Arc;

pub enum Request {
    /// The match picked: read before anything else
    Now(SpectrumKey),
    /// Rows in view: read next
    Soon(Vec<SpectrumKey>),
    /// Every shown match (replacing the last list): files read whole keep
    /// these, and are read in the background
    Shown(Vec<SpectrumKey>),
}

pub enum Event {
    Spectra(Vec<(SpectrumKey, Peaks)>),
    /// Asked for but not in their file
    NotFound(Vec<SpectrumKey>),
    /// A file could not be opened or read
    Failed(String),
    /// The file being read, or `None` once everything asked for is read
    Busy(Option<String>),
}

pub struct Loader {
    requests: Sender<Request>,
    pub events: Receiver<Event>,
}

impl Loader {
    pub fn start(results: Arc<Results>, ctx: egui::Context) -> Loader {
        let (requests, incoming) = channel();
        let (events, received) = channel();
        let _ = std::thread::Builder::new()
            .name("spectra".into())
            .spawn(move || serve(&results, &incoming, &events, &ctx));
        Loader { requests, events: received }
    }

    pub fn send(&self, request: Request) {
        let _ = self.requests.send(request);
    }
}

/// Up to this many rows in view are read in one go
const SOON_BATCH: usize = 64;

/// file → spectrum ids
type Ids = HashMap<String, HashSet<String>>;

fn is_done(done: &Ids, file: &str, id: &str) -> bool {
    done.get(file).is_some_and(|ids| ids.contains(id))
}

/// The first file read whole with a shown match not read yet
fn next_whole(shown: &Ids, done: &Ids) -> Option<String> {
    let mut files: Vec<&String> =
        shown.iter().filter(|(f, ids)| ids.iter().any(|id| !is_done(done, f, id))).map(|(f, _)| f).collect();
    files.sort();
    files.first().map(|f| f.to_string())
}

/// Run a reader; a damaged file that makes it panic gives an error instead
/// of ending the reading thread
fn guarded<T>(file: &str, read: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(read)).unwrap_or_else(|_| Err(format!("{file}: could not be read")))
}

/// The reading thread: ends when the viewer drops its [`Loader`]
fn serve(results: &Results, incoming: &Receiver<Request>, events: &Sender<Event>, ctx: &egui::Context) {
    let mut files: HashMap<String, Result<SpectrumFile, String>> = HashMap::new();
    let mut done: Ids = HashMap::new();
    let mut now: Vec<SpectrumKey> = Vec::new();
    let mut soon: Vec<SpectrumKey> = Vec::new();
    // files read whole: the ids of their shown matches
    let mut shown: Ids = HashMap::new();
    let mut busy = false;
    let send = |event: Event| {
        let sent = events.send(event).is_ok();
        ctx.request_repaint();
        sent
    };
    // Rows in view are read only from files already open: opening a timsTOF
    // run takes about a second, which a match picked meanwhile would wait
    // for. A file is opened when one of its matches is picked.
    let first_open = |soon: &[SpectrumKey], files: &HashMap<String, Result<SpectrumFile, String>>| {
        soon.iter().find(|(f, _)| matches!(files.get(f), Some(Ok(_)))).map(|(f, _)| f.clone())
    };
    loop {
        soon.retain(|(f, _)| !matches!(files.get(f), Some(Err(_))));
        let idle = now.is_empty() && first_open(&soon, &files).is_none() && next_whole(&shown, &done).is_none();
        if idle && busy {
            busy = false;
            if !send(Event::Busy(None)) {
                return;
            }
        }
        // wait for work when there is none; take every request waiting
        let first = if idle {
            match incoming.recv() {
                Ok(r) => Some(r),
                Err(_) => return,
            }
        } else {
            None
        };
        for request in first.into_iter().chain(incoming.try_iter()) {
            match request {
                Request::Now(key) => now.push(key),
                Request::Soon(keys) => soon.extend(keys),
                Request::Shown(keys) => {
                    shown.clear();
                    for (file, id) in keys {
                        if results.read_whole(&file) {
                            shown.entry(file).or_default().insert(id);
                        }
                    }
                }
            }
        }

        // One file at a time, so that a match picked meanwhile goes next: the
        // match picked, else rows in view (of open files), else a file read whole
        let take = |queue: &mut Vec<SpectrumKey>, file: &str, limit: usize| -> HashSet<String> {
            let mut ids = HashSet::new();
            queue.retain(|(f, id)| {
                if f == file && ids.len() < limit {
                    ids.insert(id.clone());
                    false
                } else {
                    true
                }
            });
            ids
        };
        let batch: Option<(String, HashSet<String>)> = if let Some((file, _)) = now.first().cloned() {
            Some((file.clone(), take(&mut now, &file, usize::MAX)))
        } else if let Some(file) = first_open(&soon, &files) {
            let ids = take(&mut soon, &file, SOON_BATCH);
            Some((file, ids))
        } else {
            next_whole(&shown, &done).map(|file| (file, HashSet::new()))
        };
        if let Some((file, mut ids)) = batch {
            ids.retain(|id| !is_done(&done, &file, id));
            // a file read whole: every shown match of it at once
            if let Some(more) = shown.get(&file) {
                ids.extend(more.iter().filter(|id| !is_done(&done, &file, id)).cloned());
            }
            if ids.is_empty() {
                continue;
            }
            busy = true;
            if !send(Event::Busy(Some(file.clone()))) {
                return;
            }
            let source = files.entry(file.clone()).or_insert_with(|| guarded(&file, || results.open_file(&file)));
            let read = match source {
                Ok(source) => guarded(&file, || source.read(&ids).map_err(|e| format!("{file}: {e}"))),
                Err(e) => Err(e.clone()),
            };
            let sent = match read {
                Ok(spectra) => {
                    let found: HashSet<&String> = spectra.iter().map(|(id, _)| id).collect();
                    let missing: Vec<SpectrumKey> =
                        ids.iter().filter(|id| !found.contains(id)).map(|id| (file.clone(), id.clone())).collect();
                    let spectra = spectra.into_iter().map(|(id, p)| ((file.clone(), id), p)).collect();
                    send(Event::Spectra(spectra)) && (missing.is_empty() || send(Event::NotFound(missing)))
                }
                Err(e) => send(Event::Failed(e)),
            };
            done.entry(file).or_default().extend(ids);
            if !sent {
                return;
            }
        }
    }
}
