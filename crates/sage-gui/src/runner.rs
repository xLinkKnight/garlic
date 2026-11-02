//! Runs external programs (RAW conversion, then `sage.exe`) on a background
//! thread, collecting their console output for the log panel.

use std::ffi::OsString;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

const MAX_LOG_LINES: usize = 20_000;

/// One program invocation in a job
pub struct Step {
    pub title: String,
    pub program: PathBuf,
    pub args: Vec<OsString>,
    /// Skip the step when this file already exists (converted RAW files)
    pub skip_if_exists: Option<PathBuf>,
    /// Cap on the program's memory in bytes (see [`crate::system::MemoryCap`])
    pub memory_limit: Option<u64>,
}

#[derive(Clone, PartialEq, Eq)]
pub enum Status {
    Running(String),
    Finished,
    Failed(String),
    Cancelled,
}

pub struct Job {
    /// The latest output lines (older ones are dropped beyond `MAX_LOG_LINES`)
    pub log: Arc<Mutex<Vec<String>>>,
    /// Lines ever written to `log`, so readers can tell which ones are new
    pub total: Arc<AtomicUsize>,
    pub status: Arc<Mutex<Status>>,
    cancel: Arc<AtomicBool>,
}

impl Job {
    pub fn start(steps: Vec<Step>, ctx: eframe::egui::Context) -> Job {
        let job = Job {
            log: Arc::new(Mutex::new(Vec::new())),
            total: Arc::new(AtomicUsize::new(0)),
            status: Arc::new(Mutex::new(Status::Running("starting".into()))),
            cancel: Arc::new(AtomicBool::new(false)),
        };
        let log = job.log.clone();
        let total = job.total.clone();
        let status = job.status.clone();
        let cancel = job.cancel.clone();

        std::thread::spawn(move || {
            let push = |line: String| push_line(&log, &total, &ctx, line);
            let result = run_steps(&steps, &push, &status, &cancel);
            let end = match result {
                Ok(()) => Status::Finished,
                Err(_) if cancel.load(Ordering::Relaxed) => Status::Cancelled,
                Err(e) => Status::Failed(e),
            };
            match &end {
                Status::Finished => push("=== Search finished ===".into()),
                Status::Cancelled => push("=== Cancelled ===".into()),
                Status::Failed(e) => push(format!("=== Failed: {e} ===")),
                Status::Running(_) => {}
            }
            *status.lock().unwrap() = end;
            ctx.request_repaint();
        });
        job
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    pub fn status(&self) -> Status {
        self.status.lock().unwrap().clone()
    }

    pub fn is_running(&self) -> bool {
        matches!(self.status(), Status::Running(_))
    }
}

fn push_line(log: &Mutex<Vec<String>>, total: &AtomicUsize, ctx: &eframe::egui::Context, line: String) {
    let mut log = log.lock().unwrap();
    log.push(line);
    total.fetch_add(1, Ordering::Relaxed);
    if log.len() > MAX_LOG_LINES {
        let excess = log.len() - MAX_LOG_LINES;
        log.drain(..excess);
    }
    ctx.request_repaint();
}

fn run_steps(
    steps: &[Step],
    push: &(dyn Fn(String) + Sync),
    status: &Mutex<Status>,
    cancel: &AtomicBool,
) -> Result<(), String> {
    for step in steps {
        if cancel.load(Ordering::Relaxed) {
            return Err("cancelled".into());
        }
        if let Some(out) = &step.skip_if_exists {
            if out.exists() {
                push(format!("--- {}: {} already exists, skipped", step.title, out.display()));
                continue;
            }
        }
        *status.lock().unwrap() = Status::Running(step.title.clone());
        push(format!("--- {} ---", step.title));
        push(format!("> {}", command_line(&step.program, &step.args)));
        run_one(step, push, cancel)?;
    }
    Ok(())
}

fn command_line(program: &Path, args: &[OsString]) -> String {
    let quote = |s: String| if s.contains(' ') { format!("\"{s}\"") } else { s };
    std::iter::once(program.as_os_str())
        .chain(args.iter().map(OsString::as_os_str))
        .map(|s| quote(s.to_string_lossy().into_owned()))
        .collect::<Vec<_>>()
        .join(" ")
}

fn run_one(step: &Step, push: &(dyn Fn(String) + Sync), cancel: &AtomicBool) -> Result<(), String> {
    let mut cmd = Command::new(&step.program);
    cmd.args(&step.args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }
    let mut child: Child = cmd
        .spawn()
        .map_err(|e| format!("could not start {}: {e}", step.program.display()))?;
    let cap = match step.memory_limit {
        Some(limit) => match crate::system::MemoryCap::apply(&child, limit) {
            Ok(cap) => {
                push(format!("RAM limit: {}", crate::system::format_memory(limit)));
                Some(cap)
            }
            Err(e) => {
                let _ = child.kill();
                return Err(e);
            }
        },
        None => None,
    };
    let over_limit = |code| match &cap {
        Some(cap) if cap.reached() => format!(
            "{} ran out of memory: it needed more than the {} RAM limit (Files > Resources); \
             raise the limit, or search fewer files at once",
            step.title,
            crate::system::format_memory(cap.limit())
        ),
        _ => format!("{} exited with {code}", step.title),
    };

    let stdout = child.stdout.take().unwrap();
    let stderr = child.stderr.take().unwrap();

    std::thread::scope(|scope| {
        scope.spawn(|| read_lines(stdout, push));
        scope.spawn(|| read_lines(stderr, push));
        loop {
            if cancel.load(Ordering::Relaxed) {
                let _ = child.kill();
                let _ = child.wait();
                return Err("cancelled".to_string());
            }
            match child.try_wait() {
                Ok(Some(code)) if code.success() => return Ok(()),
                Ok(Some(code)) => return Err(over_limit(code)),
                Ok(None) => std::thread::sleep(Duration::from_millis(100)),
                Err(e) => return Err(e.to_string()),
            }
        }
    })
}

/// Forward output line by line. Progress lines that redraw themselves with
/// `\r` (msconvert) are reduced to their final state.
fn read_lines(source: impl Read, push: &(dyn Fn(String) + Sync)) {
    for line in BufReader::new(source).split(b'\n').map_while(Result::ok) {
        let text = String::from_utf8_lossy(&line);
        let last = text
            .trim_end_matches('\r')
            .rsplit('\r')
            .next()
            .unwrap_or_default()
            .to_string();
        if !last.trim().is_empty() {
            push(last);
        }
    }
}

/// `sage.exe` next to this program (both are built into `target\release`)
pub fn find_sage() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    let name = if cfg!(windows) { "sage.exe" } else { "sage" };
    let candidate = exe.parent()?.join(name);
    candidate.is_file().then_some(candidate)
}

/// The newest ProteoWizard msconvert in the usual install locations
pub fn find_msconvert() -> Option<PathBuf> {
    let mut roots: Vec<PathBuf> = Vec::new();
    for var in ["ProgramFiles", "ProgramFiles(x86)"] {
        if let Ok(dir) = std::env::var(var) {
            roots.push(PathBuf::from(dir).join("ProteoWizard"));
        }
    }
    if let Ok(dir) = std::env::var("LOCALAPPDATA") {
        roots.push(PathBuf::from(dir).join("Apps"));
    }
    let mut found: Vec<PathBuf> = Vec::new();
    for root in roots {
        let direct = root.join("msconvert.exe");
        if direct.is_file() {
            found.push(direct);
        }
        let Ok(entries) = std::fs::read_dir(&root) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path().join("msconvert.exe");
            let is_pwiz = entry.file_name().to_string_lossy().contains("ProteoWizard");
            if (is_pwiz || root.ends_with("ProteoWizard")) && path.is_file() {
                found.push(path);
            }
        }
    }
    // Install folders carry the version, so the last name is the newest
    found.sort();
    found.pop()
}

/// The command that converts one Thermo .raw file to mzML in `out_dir`.
/// Supports ProteoWizard msconvert and ThermoRawFileParser.
pub fn conversion_step(converter: &Path, raw: &Path, out_dir: &Path) -> Step {
    let output = converted_path(raw, out_dir);
    let exe = converter
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_ascii_lowercase();
    let args: Vec<OsString> = if exe.contains("thermorawfileparser") {
        vec![
            "-i".into(),
            raw.into(),
            "-o".into(),
            out_dir.into(),
            "-f".into(),
            "2".into(),
        ]
    } else {
        vec![
            raw.into(),
            "--mzML".into(),
            "--zlib".into(),
            "--filter".into(),
            "peakPicking vendor msLevel=1-".into(),
            "-o".into(),
            out_dir.into(),
        ]
    };
    Step {
        title: format!("Convert {}", raw.file_name().unwrap_or_default().to_string_lossy()),
        program: converter.to_path_buf(),
        args,
        skip_if_exists: Some(output),
        memory_limit: None,
    }
}

pub fn converted_path(raw: &Path, out_dir: &Path) -> PathBuf {
    let stem = raw.file_stem().unwrap_or_default().to_string_lossy();
    out_dir.join(format!("{stem}.mzML"))
}
