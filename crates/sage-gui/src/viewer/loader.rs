//! Reads the spectra of matches in the background, so the viewer never waits
//! for a whole run: the match picked first, then the rows in view of runs
//! already open. Bruker runs are read spectrum by spectrum, by id, exactly as
//! the search read them. Opening a run takes about a second, so the runs of
//! the shown matches are opened in the background, the best match's run
//! first, by threads of their own (a run picked before its turn is opened at
//! once). Files without an index of spectra (mzML, MGF) are read whole, once,
//! keeping the spectra of every shown match, one after another.

use super::model::{Peaks, Results, SpectrumFile, SpectrumKey};
use eframe::egui;
use std::collections::{HashMap, HashSet};
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::{Arc, Mutex};

#[derive(Debug)]
pub enum Request {
    /// The match picked: read before anything else
    Now(SpectrumKey),
    /// Rows in view: read next
    Soon(Vec<SpectrumKey>),
    /// Every shown match (replacing the last list): files read whole keep
    /// these, and are read in the background; their timsTOF runs are opened
    Shown(Vec<SpectrumKey>),
    /// A run opened in the background
    Opened(String, Result<SpectrumFile, String>),
    /// The viewer is done with the reader
    Stop,
}

pub enum Event {
    Spectra(Vec<(SpectrumKey, Peaks)>),
    /// Asked for but not in their file
    NotFound(Vec<SpectrumKey>),
    /// A file could not be opened or read
    Failed(String),
    /// The file being read, or `None` once everything asked for is read
    Busy(Option<String>),
    /// timsTOF runs of the shown matches open, of all of them
    Runs { open: usize, total: usize },
}

pub struct Loader {
    requests: Sender<Request>,
    pub events: Receiver<Event>,
}

/// Threads opening runs in the background: two at a time take half as long
/// as one, and leave the rest of the computer free
const OPENERS: usize = 2;

impl Loader {
    pub fn start(results: Arc<Results>, ctx: egui::Context) -> Loader {
        let (requests, incoming) = channel();
        let (events, received) = channel();
        let (to_open, queue) = channel::<String>();
        let queue = Arc::new(Mutex::new(queue));
        // runs opened, or being opened, by either side
        let claimed: Arc<Mutex<HashSet<String>>> = Arc::default();
        for n in 0..OPENERS {
            let (results, queue, claimed, back) = (results.clone(), queue.clone(), claimed.clone(), requests.clone());
            let _ = std::thread::Builder::new()
                .name(format!("open runs {n}"))
                .spawn(move || open_runs(&results, &queue, &claimed, &back));
        }
        let _ = std::thread::Builder::new()
            .name("spectra".into())
            .spawn(move || serve(&results, &incoming, &events, &ctx, &to_open, &claimed));
        Loader { requests, events: received }
    }

    pub fn send(&self, request: Request) {
        let _ = self.requests.send(request);
    }

    /// A loader without a reading thread: the test reads its requests
    #[cfg(test)]
    pub fn fake() -> (Loader, Receiver<Request>) {
        let (requests, incoming) = channel();
        let (_, events) = channel();
        (Loader { requests, events }, incoming)
    }
}

impl Drop for Loader {
    /// The opener threads hold the request channel open too: stop the reader
    /// (which lets the openers end) when the viewer lets go
    fn drop(&mut self) {
        let _ = self.requests.send(Request::Stop);
    }
}

/// An opener thread: opens the runs queued, unless the reader has opened
/// them meanwhile, and hands them to the reader. Ends with the reader.
fn open_runs(results: &Results, queue: &Mutex<Receiver<String>>, claimed: &Mutex<HashSet<String>>, back: &Sender<Request>) {
    loop {
        let next = queue.lock().map_err(|_| ()).and_then(|q| q.recv().map_err(|_| ()));
        let Ok(file) = next else { return };
        if !claimed.lock().is_ok_and(|mut c| c.insert(file.clone())) {
            continue;
        }
        let run = guarded(&file, || results.open_file(&file));
        if back.send(Request::Opened(file, run)).is_err() {
            return;
        }
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
fn serve(
    results: &Results,
    incoming: &Receiver<Request>,
    events: &Sender<Event>,
    ctx: &egui::Context,
    to_open: &Sender<String>,
    claimed: &Mutex<HashSet<String>>,
) {
    let mut files: HashMap<String, Result<SpectrumFile, String>> = HashMap::new();
    let mut done: Ids = HashMap::new();
    let mut now: Vec<SpectrumKey> = Vec::new();
    let mut soon: Vec<SpectrumKey> = Vec::new();
    // files read whole: the ids of their shown matches
    let mut shown: Ids = HashMap::new();
    // timsTOF runs queued for opening in the background
    let mut queued: HashSet<String> = HashSet::new();
    let mut busy = false;
    let send = |event: Event| {
        let sent = events.send(event).is_ok();
        ctx.request_repaint();
        sent
    };
    let runs = |files: &HashMap<String, Result<SpectrumFile, String>>, queued: &HashSet<String>| Event::Runs {
        open: queued.iter().filter(|f| files.contains_key(*f)).count(),
        total: queued.len(),
    };
    // Rows in view are read only from files already open: opening a timsTOF
    // run takes about a second, which a match picked meanwhile would wait
    // for. Runs are opened in the background, or when one of their matches
    // is picked.
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
                    let before = queued.len();
                    for (file, id) in keys {
                        if results.read_whole(&file) {
                            shown.entry(file).or_default().insert(id);
                        } else if !files.contains_key(&file) && queued.insert(file.clone()) {
                            // the shown matches' runs, the best match's first
                            let _ = to_open.send(file);
                        }
                    }
                    if queued.len() != before && !send(runs(&files, &queued)) {
                        return;
                    }
                }
                Request::Opened(file, run) => {
                    // (unless a match of it was picked and the run opened meanwhile)
                    files.entry(file).or_insert(run);
                    if !send(runs(&files, &queued)) {
                        return;
                    }
                }
                Request::Stop => return,
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
            let opening = !files.contains_key(&file);
            let source = files.entry(file.clone()).or_insert_with(|| {
                // the openers skip a run opened here
                if let Ok(mut c) = claimed.lock() {
                    c.insert(file.clone());
                }
                guarded(&file, || results.open_file(&file))
            });
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
            if !sent || (opening && !queued.is_empty() && !send(runs(&files, &queued))) {
                return;
            }
        }
    }
}
