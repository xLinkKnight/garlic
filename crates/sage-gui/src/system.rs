//! What this computer offers a search (CPU threads, RAM, graphics adapters),
//! whether `sage.exe` can use a GPU, and the RAM cap placed on `sage.exe`.

use std::path::Path;
use std::process::{Child, Command};
use std::sync::mpsc::{channel, Receiver};

/// Logical CPUs (hardware threads)
pub fn cpu_threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

pub const GB: f64 = 1024.0 * 1024.0 * 1024.0;

/// A memory amount for people: "512 MB", "4.0 GB"
pub fn format_memory(bytes: u64) -> String {
    if (bytes as f64) < GB {
        format!("{:.0} MB", bytes as f64 / (1024.0 * 1024.0))
    } else {
        format!("{:.1} GB", bytes as f64 / GB)
    }
}

/// Installed RAM in bytes
#[cfg(windows)]
pub fn total_ram() -> Option<u64> {
    use windows::Win32::System::SystemInformation::{GlobalMemoryStatusEx, MEMORYSTATUSEX};
    let mut status = MEMORYSTATUSEX {
        dwLength: std::mem::size_of::<MEMORYSTATUSEX>() as u32,
        ..Default::default()
    };
    unsafe { GlobalMemoryStatusEx(&mut status) }.ok()?;
    Some(status.ullTotalPhys)
}

#[cfg(not(windows))]
pub fn total_ram() -> Option<u64> {
    None
}

/// A graphics adapter as Windows reports it
#[derive(Clone, Debug)]
pub struct Gpu {
    pub name: String,
    /// Dedicated video memory in bytes (0 for integrated GPUs that share system RAM)
    pub vram: u64,
    /// Integrated GPUs: system RAM Windows lets the GPU borrow
    pub shared: u64,
    /// Adapter id as Windows performance counters spell it: `luid_0x00000000_0x0000D1C2`
    pub luid: String,
}

/// Hardware adapters (the Microsoft Basic Render Driver is left out)
#[cfg(windows)]
pub fn gpus() -> Vec<Gpu> {
    use windows::Win32::Graphics::Dxgi::{CreateDXGIFactory1, IDXGIFactory1, DXGI_ADAPTER_FLAG_SOFTWARE};
    let Ok(factory) = (unsafe { CreateDXGIFactory1::<IDXGIFactory1>() }) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for i in 0.. {
        let Ok(adapter) = (unsafe { factory.EnumAdapters1(i) }) else { break };
        let Ok(desc) = (unsafe { adapter.GetDesc1() }) else { continue };
        let software = desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE.0 as u32 != 0;
        // 0x1414 = Microsoft: the Basic Render Driver is a CPU fallback
        if software || desc.VendorId == 0x1414 {
            continue;
        }
        let len = desc.Description.iter().position(|&c| c == 0).unwrap_or(desc.Description.len());
        out.push(Gpu {
            name: String::from_utf16_lossy(&desc.Description[..len]).trim().to_string(),
            vram: desc.DedicatedVideoMemory as u64,
            shared: desc.SharedSystemMemory as u64,
            luid: format!("luid_0x{:08X}_0x{:08X}", desc.AdapterLuid.HighPart as u32, desc.AdapterLuid.LowPart),
        });
    }
    out
}

#[cfg(not(windows))]
pub fn gpus() -> Vec<Gpu> {
    Vec::new()
}

/// What `sage.exe --gpu-info` says about GPU support
#[derive(Clone, Debug, PartialEq)]
pub enum GpuSupport {
    /// Built with `--features gpu`; the adapters it can use (`name | backend | type`)
    Yes(Vec<String>),
    /// Built without GPU support
    No,
    /// sage.exe could not be asked (missing, or too old to know `--gpu-info`)
    Unknown(String),
}

/// Ask `sage.exe --gpu-info`
pub fn probe_sage_gpu(sage: &Path) -> GpuSupport {
    if !sage.is_file() {
        return GpuSupport::Unknown("sage.exe not found".into());
    }
    let mut cmd = Command::new(sage);
    cmd.arg("--gpu-info");
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }
    let out = match cmd.output() {
        Ok(out) => out,
        Err(e) => return GpuSupport::Unknown(format!("could not run sage.exe: {e}")),
    };
    parse_gpu_info(&String::from_utf8_lossy(&out.stdout))
}

fn parse_gpu_info(text: &str) -> GpuSupport {
    let value = |key: &str| text.lines().find_map(|l| l.strip_prefix(key)).map(str::trim);
    match value("gpu-support:") {
        Some(v) if v.starts_with("yes") => GpuSupport::Yes(
            text.lines()
                .filter_map(|l| l.strip_prefix("gpu-adapter:"))
                .map(|l| l.trim().to_string())
                .collect(),
        ),
        Some(_) => GpuSupport::No,
        None => GpuSupport::Unknown("this sage.exe is too old to report GPU support; rebuild it".into()),
    }
}

/// Probe `sage.exe` on a background thread
pub fn probe_in_background(sage: &Path, ctx: &eframe::egui::Context) -> Receiver<GpuSupport> {
    let (tx, rx) = channel();
    let (sage, ctx) = (sage.to_path_buf(), ctx.clone());
    std::thread::spawn(move || {
        let _ = tx.send(probe_sage_gpu(&sage));
        ctx.request_repaint();
    });
    rx
}

/// A Windows job object that caps the memory of the process placed in it
#[cfg(windows)]
pub struct MemoryCap {
    job: windows::Win32::Foundation::HANDLE,
    /// Receives the job's notifications, e.g. that the memory limit was hit
    port: windows::Win32::Foundation::HANDLE,
    limit: u64,
}

/// Job notification: a process in the job exceeded its memory limit
#[cfg(windows)]
const JOB_OBJECT_MSG_PROCESS_MEMORY_LIMIT: u32 = 9;

#[cfg(windows)]
impl MemoryCap {
    /// Limit `child` to `limit` bytes of committed memory. An allocation
    /// beyond it fails, which ends `sage.exe` instead of exhausting the system.
    pub fn apply(child: &Child, limit: u64) -> Result<Self, String> {
        use std::os::windows::io::AsRawHandle;
        use windows::Win32::Foundation::{HANDLE, INVALID_HANDLE_VALUE};
        use windows::Win32::System::JobObjects::*;
        use windows::Win32::System::IO::CreateIoCompletionPort;
        unsafe {
            let job = CreateJobObjectW(None, None).map_err(|e| format!("could not create a job object: {e}"))?;
            let port = CreateIoCompletionPort(INVALID_HANDLE_VALUE, None, 0, 1).unwrap_or_default();
            let cap = MemoryCap { job, port, limit };
            if !port.is_invalid() {
                let assoc = JOBOBJECT_ASSOCIATE_COMPLETION_PORT {
                    CompletionKey: std::ptr::null_mut(),
                    CompletionPort: port,
                };
                let _ = SetInformationJobObject(
                    job,
                    JobObjectAssociateCompletionPortInformation,
                    &assoc as *const _ as *const _,
                    std::mem::size_of_val(&assoc) as u32,
                );
            }
            let mut info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION::default();
            info.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_PROCESS_MEMORY;
            info.ProcessMemoryLimit = limit as usize;
            SetInformationJobObject(
                job,
                JobObjectExtendedLimitInformation,
                &info as *const _ as *const _,
                std::mem::size_of_val(&info) as u32,
            )
            .map_err(|e| format!("could not set the memory limit: {e}"))?;
            AssignProcessToJobObject(job, HANDLE(child.as_raw_handle()))
                .map_err(|e| format!("could not apply the memory limit: {e}"))?;
            Ok(cap)
        }
    }

    /// Was an allocation refused because of the limit? (Reads the job's
    /// pending notifications, so ask once, after the process has exited.)
    pub fn reached(&self) -> bool {
        use windows::Win32::System::IO::{GetQueuedCompletionStatus, OVERLAPPED};
        if self.port.is_invalid() {
            return false;
        }
        let mut hit = false;
        loop {
            let (mut msg, mut key, mut overlapped) = (0u32, 0usize, std::ptr::null_mut::<OVERLAPPED>());
            let got = unsafe { GetQueuedCompletionStatus(self.port, &mut msg, &mut key, &mut overlapped, 0) };
            if got.is_err() {
                return hit;
            }
            hit |= msg == JOB_OBJECT_MSG_PROCESS_MEMORY_LIMIT;
        }
    }

    pub fn limit(&self) -> u64 {
        self.limit
    }
}

#[cfg(windows)]
impl Drop for MemoryCap {
    fn drop(&mut self) {
        unsafe {
            let _ = windows::Win32::Foundation::CloseHandle(self.job);
            if !self.port.is_invalid() {
                let _ = windows::Win32::Foundation::CloseHandle(self.port);
            }
        }
    }
}

/// No memory cap outside Windows: the limit is ignored
#[cfg(not(windows))]
pub struct MemoryCap {
    limit: u64,
}

#[cfg(not(windows))]
impl MemoryCap {
    pub fn apply(_child: &Child, limit: u64) -> Result<Self, String> {
        Ok(MemoryCap { limit })
    }

    pub fn reached(&self) -> bool {
        false
    }

    pub fn limit(&self) -> u64 {
        self.limit
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn gpu_info_lines() {
        assert_eq!(
            parse_gpu_info("gpu-support: yes\ngpu-adapter: NVIDIA RTX 4070 | Vulkan | DiscreteGpu\n"),
            GpuSupport::Yes(vec!["NVIDIA RTX 4070 | Vulkan | DiscreteGpu".into()])
        );
        assert_eq!(parse_gpu_info("gpu-support: no (rebuild ...)"), GpuSupport::No);
        assert!(matches!(parse_gpu_info("error: unexpected argument"), GpuSupport::Unknown(_)));
    }

    #[test]
    fn system_queries() {
        assert!(cpu_threads() >= 1);
        if cfg!(windows) {
            assert!(total_ram().unwrap() > 0);
        }
        // must not panic on machines without a GPU
        let _ = gpus();
    }

    /// Child half of `memory_cap_stops_a_process`: allocate 1 GB
    #[test]
    fn allocate_one_gb() {
        if std::env::var_os("SAGE_GUI_ALLOCATE").is_some() {
            std::thread::sleep(std::time::Duration::from_millis(300));
            let v = vec![1u8; 1 << 30];
            assert_eq!(std::hint::black_box(&v)[v.len() - 1], 1);
        }
    }

    /// A process that allocates past its cap is stopped, and the cap reports it
    #[cfg(windows)]
    #[test]
    fn memory_cap_stops_a_process() {
        let run = |cap_bytes: u64| {
            let mut child = Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "system::test::allocate_one_gb", "--test-threads", "1"])
                .env("SAGE_GUI_ALLOCATE", "1")
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .spawn()
                .unwrap();
            let cap = MemoryCap::apply(&child, cap_bytes).unwrap();
            (child.wait().unwrap().success(), cap.reached())
        };
        // 300 MB cap: the allocation fails
        assert_eq!(run(300 << 20), (false, true));
        // 4 GB cap: room enough
        assert_eq!(run(4 << 30), (true, false));
    }
}
