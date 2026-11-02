//! Live CPU, RAM, GPU and vRAM utilisation of the whole computer, sampled
//! once a second on a background thread while someone is looking.
//!
//! Sources match Task Manager: CPU from `GetSystemTimes`, RAM from
//! `GlobalMemoryStatusEx`, GPU from the "GPU Engine" and "GPU Adapter Memory"
//! performance counters, vRAM totals from DXGI.

use crate::system::Gpu;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GpuUsage {
    pub name: String,
    /// Busiest engine type (3D, compute, copy, video…), 0–100, as Task Manager reports
    pub util: f32,
    pub vram_used: u64,
    pub vram_total: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Usage {
    pub cpu: f32,
    pub ram_used: u64,
    pub ram_total: u64,
    pub gpus: Vec<GpuUsage>,
}

impl Usage {
    pub fn ram_pct(&self) -> f32 {
        if self.ram_total == 0 {
            0.0
        } else {
            100.0 * self.ram_used as f32 / self.ram_total as f32
        }
    }

    /// The GPU worth showing: the one with the most dedicated memory (the discrete card)
    pub fn main_gpu(&self) -> Option<&GpuUsage> {
        self.gpus.iter().max_by_key(|g| g.vram_total)
    }
}

impl GpuUsage {
    pub fn vram_pct(&self) -> Option<f32> {
        (self.vram_total > 0).then(|| 100.0 * self.vram_used as f32 / self.vram_total as f32)
    }
}

/// Samples in the background; `latest` returns the newest sample
pub struct Monitor {
    latest: Arc<Mutex<Option<Usage>>>,
    /// Set by the GUI each frame it shows the numbers; the sampler idles otherwise
    wanted: Arc<AtomicBool>,
}

impl Monitor {
    pub fn start(gpus: Vec<Gpu>) -> Self {
        let monitor = Monitor { latest: Arc::default(), wanted: Arc::default() };
        let (latest, wanted) = (monitor.latest.clone(), monitor.wanted.clone());
        std::thread::spawn(move || {
            let mut sampler = Sampler::new(gpus);
            loop {
                if wanted.swap(false, Ordering::Relaxed) {
                    *latest.lock().unwrap() = sampler.sample();
                    std::thread::sleep(Duration::from_millis(1000));
                } else {
                    std::thread::sleep(Duration::from_millis(250));
                }
            }
        });
        monitor
    }

    /// The newest sample; also keeps sampling going and the screen updating
    pub fn latest(&self, ctx: &eframe::egui::Context) -> Option<Usage> {
        self.wanted.store(true, Ordering::Relaxed);
        ctx.request_repaint_after(Duration::from_millis(1000));
        self.latest.lock().unwrap().clone()
    }
}

/// Task Manager's GPU figure from per-engine counter instances
/// (`pid_123_luid_0x…_0x…_phys_0_eng_0_engtype_3D` → value): per adapter, sum
/// each engine type over all processes, then take the busiest type.
fn gpu_utilisation(items: &[(String, f64)], luid: &str) -> f32 {
    let mut by_type: std::collections::HashMap<&str, f64> = Default::default();
    for (name, value) in items {
        if !name.to_ascii_lowercase().contains(&luid.to_ascii_lowercase()) {
            continue;
        }
        let engtype = name.rsplit_once("engtype_").map_or("", |(_, t)| t);
        *by_type.entry(engtype).or_default() += value;
    }
    by_type.values().copied().fold(0.0, f64::max).min(100.0) as f32
}

/// Dedicated memory in use on one adapter, from "GPU Adapter Memory" instances
fn vram_used(items: &[(String, f64)], luid: &str) -> u64 {
    items
        .iter()
        .filter(|(name, _)| name.to_ascii_lowercase().contains(&luid.to_ascii_lowercase()))
        .map(|(_, v)| *v)
        .sum::<f64>() as u64
}

#[cfg(windows)]
struct Sampler {
    gpus: Vec<Gpu>,
    /// Previous (idle, busy) CPU time
    cpu_prev: Option<(u64, u64)>,
    query: isize,
    engine: isize,
    memory: isize,
}

#[cfg(windows)]
impl Sampler {
    fn new(gpus: Vec<Gpu>) -> Self {
        use windows::core::w;
        use windows::Win32::System::Performance::*;
        let (mut query, mut engine, mut memory) = (0isize, 0isize, 0isize);
        unsafe {
            if PdhOpenQueryW(None, 0, &mut query) == 0 {
                PdhAddEnglishCounterW(query, w!("\\GPU Engine(*)\\Utilization Percentage"), 0, &mut engine);
                PdhAddEnglishCounterW(query, w!("\\GPU Adapter Memory(*)\\Dedicated Usage"), 0, &mut memory);
                // rates need two collections; the first one primes them
                PdhCollectQueryData(query);
            } else {
                query = 0;
            }
        }
        Sampler { gpus, cpu_prev: None, query, engine, memory }
    }

    fn cpu(&mut self) -> f32 {
        use windows::Win32::Foundation::FILETIME;
        use windows::Win32::System::Threading::GetSystemTimes;
        let (mut idle, mut kernel, mut user) = (FILETIME::default(), FILETIME::default(), FILETIME::default());
        if unsafe { GetSystemTimes(Some(&mut idle), Some(&mut kernel), Some(&mut user)) }.is_err() {
            return 0.0;
        }
        let ticks = |t: FILETIME| ((t.dwHighDateTime as u64) << 32) | t.dwLowDateTime as u64;
        // kernel time includes idle time
        let (idle, busy) = (ticks(idle), ticks(kernel) + ticks(user) - ticks(idle));
        let pct = match self.cpu_prev {
            Some((pi, pb)) if (idle - pi) + (busy - pb) > 0 => {
                100.0 * (busy - pb) as f32 / ((idle - pi) + (busy - pb)) as f32
            }
            _ => 0.0,
        };
        self.cpu_prev = Some((idle, busy));
        pct
    }

    /// Every instance of a wildcard counter as (instance name, value)
    fn counter_items(&self, counter: isize) -> Vec<(String, f64)> {
        use windows::Win32::System::Performance::*;
        if counter == 0 {
            return Vec::new();
        }
        unsafe {
            let (mut size, mut count) = (0u32, 0u32);
            let status = PdhGetFormattedCounterArrayW(counter, PDH_FMT_DOUBLE, &mut size, &mut count, None);
            if status != PDH_MORE_DATA || size == 0 {
                return Vec::new();
            }
            // the buffer holds the items followed by their names
            let mut buf = vec![0u64; (size as usize).div_ceil(8)];
            let items = buf.as_mut_ptr() as *mut PDH_FMT_COUNTERVALUE_ITEM_W;
            if PdhGetFormattedCounterArrayW(counter, PDH_FMT_DOUBLE, &mut size, &mut count, Some(items)) != 0 {
                return Vec::new();
            }
            std::slice::from_raw_parts(items, count as usize)
                .iter()
                .filter(|item| item.FmtValue.CStatus == 0)
                .map(|item| (item.szName.to_string().unwrap_or_default(), item.FmtValue.Anonymous.doubleValue))
                .collect()
        }
    }

    fn sample(&mut self) -> Option<Usage> {
        use windows::Win32::System::Performance::PdhCollectQueryData;
        use windows::Win32::System::SystemInformation::{GlobalMemoryStatusEx, MEMORYSTATUSEX};
        let cpu = self.cpu();
        let mut mem = MEMORYSTATUSEX { dwLength: std::mem::size_of::<MEMORYSTATUSEX>() as u32, ..Default::default() };
        unsafe { GlobalMemoryStatusEx(&mut mem) }.ok()?;
        let (engines, memory) = if self.query != 0 && unsafe { PdhCollectQueryData(self.query) } == 0 {
            (self.counter_items(self.engine), self.counter_items(self.memory))
        } else {
            (Vec::new(), Vec::new())
        };
        let gpus = self
            .gpus
            .iter()
            .map(|g| GpuUsage {
                name: g.name.clone(),
                util: gpu_utilisation(&engines, &g.luid),
                vram_used: vram_used(&memory, &g.luid),
                vram_total: g.vram,
            })
            .collect();
        Some(Usage {
            cpu,
            ram_used: mem.ullTotalPhys - mem.ullAvailPhys,
            ram_total: mem.ullTotalPhys,
            gpus,
        })
    }
}

#[cfg(windows)]
impl Drop for Sampler {
    fn drop(&mut self) {
        if self.query != 0 {
            unsafe {
                windows::Win32::System::Performance::PdhCloseQuery(self.query);
            }
        }
    }
}

/// No usage figures outside Windows
#[cfg(not(windows))]
struct Sampler;

#[cfg(not(windows))]
impl Sampler {
    fn new(_gpus: Vec<Gpu>) -> Self {
        Sampler
    }

    fn sample(&mut self) -> Option<Usage> {
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn gpu_figures_follow_task_manager() {
        let luid = "luid_0x00000000_0x0000D1C2";
        let items: Vec<(String, f64)> = [
            ("pid_10_luid_0x00000000_0x0000D1C2_phys_0_eng_0_engtype_3D", 30.0),
            ("pid_11_luid_0x00000000_0x0000D1C2_phys_0_eng_0_engtype_3D", 25.0),
            ("pid_11_luid_0x00000000_0x0000D1C2_phys_0_eng_5_engtype_Compute_0", 40.0),
            // another adapter
            ("pid_12_luid_0x00000000_0x0000AAAA_phys_0_eng_0_engtype_3D", 90.0),
        ]
        .iter()
        .map(|(n, v)| (n.to_string(), *v))
        .collect();
        // 3D: 30 + 25 = 55 beats Compute 40; the other adapter is ignored
        assert_eq!(gpu_utilisation(&items, luid), 55.0);
        // case of the hex digits does not matter
        assert_eq!(gpu_utilisation(&items, &luid.to_ascii_lowercase()), 55.0);

        let mem = vec![
            ("luid_0x00000000_0x0000D1C2_phys_0".to_string(), 1.5e9),
            ("luid_0x00000000_0x0000AAAA_phys_0".to_string(), 9.0e9),
        ];
        assert_eq!(vram_used(&mem, luid), 1_500_000_000);
    }

    /// A real sample on this machine: plausible numbers, matching the GPU list
    #[cfg(windows)]
    #[test]
    fn samples_this_computer() {
        let gpus = crate::system::gpus();
        let mut s = Sampler::new(gpus.clone());
        s.sample();
        std::thread::sleep(Duration::from_millis(600));
        let u = s.sample().unwrap();
        assert!((0.0..=100.0).contains(&u.cpu), "{u:?}");
        assert!(u.ram_total > 0 && u.ram_used <= u.ram_total);
        assert_eq!(u.gpus.len(), gpus.len());
        for g in &u.gpus {
            assert!((0.0..=100.0).contains(&g.util), "{g:?}");
        }
        eprintln!("{u:#?}");
    }
}
