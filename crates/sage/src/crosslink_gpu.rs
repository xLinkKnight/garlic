//! GPU implementation of the chain-scoring kernel (see `crosslink_kernel`).
//!
//! Built with `--features gpu`; uses `wgpu`, so the same code runs on DirectX 12
//! (Windows), Vulkan (Linux/Android) and Metal (macOS).
//!
//! One GPU thread scores one peptide, exactly as
//! [`crate::crosslink_kernel::score_peptide`] does on the CPU. The ion table is
//! uploaded once and stays resident; per spectrum only the peaks are uploaded
//! and two small result arrays (score, packed site/matches) are read back.
//! Tables larger than the device's buffer limits are split into chunks of
//! consecutive peptides and dispatched one after another.

use crate::crosslink_kernel::{
    lnfact_table, top_survivors, ChainKernel, IonTable, KernelParams, Survivor,
};
use crate::spectrum::ProcessedSpectrum;
use std::borrow::Cow;
use std::sync::Mutex;
use wgpu::util::DeviceExt;

const WORKGROUP: u32 = 64;
const MAX_WORKGROUPS: u32 = 65_535;

const SHADER: &str = r#"
struct Params {
    n_peaks: u32,
    n_pep: u32,
    max_charge: u32,
    min_matched: u32,
    ion_base: u32,
    site_base: u32,
    pad0: u32,
    pad1: u32,
    precursor_mass: f32,
    pad2: f32,
    // 1 + rel of the fragment tolerance (LinearTol::factors)
    lo_factor: f32,
    hi_factor: f32,
    abs_lo: f32,
    abs_hi: f32,
    pad3: f32,
    pad4: f32,
    // per fragment charge z = 1, 2, ...: x = 1/z, y = slack/z (KernelParams::charge_terms)
    charges: array<vec4<f32>, 32>,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> peaks_mass: array<f32>;
@group(0) @binding(2) var<storage, read> peaks_int: array<f32>;
@group(0) @binding(3) var<storage, read> ion_off: array<u32>;
@group(0) @binding(4) var<storage, read> ion_mass: array<f32>;
@group(0) @binding(5) var<storage, read> ion_meta: array<u32>;
@group(0) @binding(6) var<storage, read> site_off: array<u32>;
@group(0) @binding(7) var<storage, read> sites: array<u32>;
@group(0) @binding(8) var<storage, read> pep_mass: array<f32>;
@group(0) @binding(9) var<storage, read> lnfact: array<f32>;
@group(0) @binding(10) var<storage, read_write> out_score: array<f32>;
@group(0) @binding(11) var<storage, read_write> out_info: array<u32>;

// Intensity of the most intense peak in [lo, hi] (the later peak wins ties),
// or -1.0 when there is none. Peaks are sorted by mass.
fn best_peak(lo: f32, hi: f32) -> f32 {
    var l = 0u;
    var r = params.n_peaks;
    while (l < r) {
        let mid = (l + r) / 2u;
        if (peaks_mass[mid] < lo) {
            l = mid + 1u;
        } else {
            r = mid;
        }
    }
    var best = -1.0;
    var i = l;
    while (i < params.n_peaks && peaks_mass[i] <= hi) {
        if (peaks_int[i] >= best) {
            best = peaks_int[i];
        }
        i = i + 1u;
    }
    return best;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let p = gid.x;
    if (p >= params.n_pep) {
        return;
    }
    out_score[p] = -1.0;
    out_info[p] = 0u;

    let i0 = ion_off[p] - params.ion_base;
    let i1 = ion_off[p + 1u] - params.ion_base;
    let s0 = site_off[p] - params.site_base;
    let s1 = site_off[p + 1u] - params.site_base;
    if (i0 == i1 || s0 == s1) {
        return;
    }
    let shift = params.precursor_mass - pep_mass[p];
    if (shift <= 0.0) {
        return;
    }
    let n_sites = min(s1 - s0, 16u);

    var cnt_b: array<u32, 16>;
    var cnt_y: array<u32, 16>;
    var int_b: array<f32, 16>;
    var int_y: array<f32, 16>;
    for (var k = 0u; k < 16u; k = k + 1u) {
        cnt_b[k] = 0u;
        cnt_y[k] = 0u;
        int_b[k] = 0.0;
        int_y[k] = 0.0;
    }

    for (var ion = i0; ion < i1; ion = ion + 1u) {
        let mass = ion_mass[ion];
        var un = 0u;
        var ui = 0.0;
        var sn = 0u;
        var si = 0.0;
        // Window edges exactly as LinearTol::bounds computes them, so that a
        // peak on an edge is in or out on every device (see its docs)
        for (var c = 1u; c < params.max_charge; c = c + 1u) {
            let k = params.charges[c - 1u];
            let t = mass * k.x;
            let lo = t * params.lo_factor + params.abs_lo;
            let hi = t * params.hi_factor + params.abs_hi;
            let a = best_peak(lo, hi);
            if (a >= 0.0) {
                un = un + 1u;
                ui = ui + a;
            }
            let ts = (mass + shift) * k.x;
            let slo = ts * params.lo_factor + params.abs_lo - k.y;
            let shi = ts * params.hi_factor + params.abs_hi + k.y;
            let b = best_peak(slo, shi);
            if (b >= 0.0) {
                sn = sn + 1u;
                si = si + b;
            }
        }
        let ion_bits = ion_meta[ion];
        let n_term = (ion_bits & 1u) == 1u;
        let idx = ion_bits >> 1u;
        for (var k = 0u; k < n_sites; k = k + 1u) {
            let site = sites[s0 + k];
            var contains = idx < site;
            if (n_term) {
                contains = idx >= site;
            }
            var n = un;
            var inten = ui;
            if (contains) {
                n = sn;
                inten = si;
            }
            if (n_term) {
                cnt_b[k] = cnt_b[k] + n;
                int_b[k] = int_b[k] + inten;
            } else {
                cnt_y[k] = cnt_y[k] + n;
                int_y[k] = int_y[k] + inten;
            }
        }
    }

    var best_score = -1.0;
    var best_info = 0u;
    for (var k = 0u; k < n_sites; k = k + 1u) {
        let matched = cnt_b[k] + cnt_y[k];
        let score = log((int_b[k] + 1.0) * (int_y[k] + 1.0))
            + lnfact[min(cnt_b[k], 255u)] + lnfact[min(cnt_y[k], 255u)];
        if (matched >= params.min_matched && (best_score < 0.0 || score > best_score)) {
            best_score = score;
            best_info = (sites[s0 + k] << 16u) | min(matched, 65535u);
        }
    }
    if (best_score >= 0.0) {
        out_score[p] = best_score;
        out_info[p] = best_info;
    }
}
"#;

/// Fragment charges the shader's `charges` table holds
const MAX_CHARGE_TERMS: usize = 32;

/// Uniform block, must match `Params` in the shader (16 four-byte words,
/// then the per-charge table)
#[repr(C)]
#[derive(Copy, Clone)]
struct GpuParams {
    n_peaks: u32,
    n_pep: u32,
    max_charge: u32,
    min_matched: u32,
    ion_base: u32,
    site_base: u32,
    pad0: u32,
    pad1: u32,
    precursor_mass: f32,
    pad2: f32,
    lo_factor: f32,
    hi_factor: f32,
    abs_lo: f32,
    abs_hi: f32,
    pad3: f32,
    pad4: f32,
    charges: [[f32; 4]; MAX_CHARGE_TERMS],
}

fn bytes_of<T: Copy>(v: &[T]) -> &[u8] {
    // Safety: plain-old-data slices of u32/f32/#[repr(C)] structs of the same
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}

/// A run of consecutive peptides small enough for the device's buffer limits
struct Chunk {
    /// first peptide (database index) in this chunk, and how many
    first: usize,
    count: usize,
    ion_base: u32,
    site_base: u32,
    ion_off: wgpu::Buffer,
    ion_mass: wgpu::Buffer,
    ion_meta: wgpu::Buffer,
    site_off: wgpu::Buffer,
    sites: wgpu::Buffer,
    pep_mass: wgpu::Buffer,
    out_score: wgpu::Buffer,
    out_info: wgpu::Buffer,
    read_score: wgpu::Buffer,
    read_info: wgpu::Buffer,
}

/// Adapters on the backends the kernel uses, as `name | backend | device type`.
/// Software adapters are listed too; the kernel skips them unless asked.
pub fn adapters() -> Vec<String> {
    all_adapters()
        .into_iter()
        .map(|a| {
            let info = a.get_info();
            format!("{} | {:?} | {:?}", info.name, info.backend, info.device_type)
        })
        .collect()
}

/// Graphics APIs used: DirectX 12 on Windows, where it is the native API (some
/// Vulkan drivers there, e.g. Intel UHD 31.0.101.2145, crash opening a device);
/// Vulkan and Metal elsewhere. Never GL: slow and noisy.
fn backends() -> wgpu::Backends {
    if cfg!(windows) {
        wgpu::Backends::DX12
    } else {
        wgpu::Backends::PRIMARY
    }
}

/// One instance for the whole program: adapters and devices must not outlive it
fn instance() -> &'static wgpu::Instance {
    static INSTANCE: std::sync::OnceLock<wgpu::Instance> = std::sync::OnceLock::new();
    INSTANCE.get_or_init(|| {
        let mut descriptor = wgpu::InstanceDescriptor::new_without_display_handle();
        descriptor.backends = backends();
        wgpu::Instance::new(descriptor)
    })
}

fn all_adapters() -> Vec<wgpu::Adapter> {
    pollster::block_on(instance().enumerate_adapters(backends()))
}

pub struct GpuKernel {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    chunks: Vec<Chunk>,
    lnfact: wgpu::Buffer,
    adapter_name: String,
    label: String,
    /// The device is used by one spectrum at a time
    lock: Mutex<()>,
}

impl GpuKernel {
    /// Upload `table` to the best available adapter. Fails (rather than
    /// silently using a CPU adapter) when there is none.
    ///
    /// Software adapters (llvmpipe, WARP, ...) are CPUs in disguise and much
    /// slower than the `cpu` backend, so they are used only when `allow_software`
    /// is set (explicit `"backend": "gpu"`, or testing).
    pub fn new(table: IonTable, allow_software: bool) -> Result<Self, String> {
        let adapter = pollster::block_on(instance().request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .map_err(|e| format!("no GPU adapter found: {e}"))?;
        let info = adapter.get_info();
        if info.device_type == wgpu::DeviceType::Cpu && !allow_software {
            return Err(format!("only a software adapter ({}) is available", info.name));
        }
        Self::on_adapter(table, adapter)
    }

    /// One kernel per adapter on this machine, software ones included, so
    /// that tests can check every device and driver
    #[cfg(test)]
    pub(crate) fn on_every_adapter(table: impl Fn() -> IonTable) -> Vec<Result<Self, String>> {
        all_adapters().into_iter().map(|a| Self::on_adapter(table(), a)).collect()
    }

    fn on_adapter(table: IonTable, adapter: wgpu::Adapter) -> Result<Self, String> {
        let info = adapter.get_info();
        let limits = adapter.limits();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("sage-crosslink"),
            required_limits: wgpu::Limits {
                max_storage_buffer_binding_size: limits.max_storage_buffer_binding_size,
                max_buffer_size: limits.max_buffer_size,
                max_storage_buffers_per_shader_stage: 12.max(
                    wgpu::Limits::default().max_storage_buffers_per_shader_stage,
                ).min(limits.max_storage_buffers_per_shader_stage),
                ..wgpu::Limits::default()
            },
            ..Default::default()
        }))
        .map_err(|e| format!("could not open GPU device: {e}"))?;

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("crosslink-chain"),
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(SHADER)),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("crosslink-chain"),
            layout: None,
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });

        let max_binding = limits.max_storage_buffer_binding_size.min(limits.max_buffer_size) as usize;
        let chunks = build_chunks(&device, &table, max_binding)?;
        let lnfact = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("lnfact"),
            contents: bytes_of(&lnfact_table()),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let label = format!("gpu ({}, {:?}, {} chunk(s))", info.name, info.backend, chunks.len());
        Ok(Self {
            device,
            queue,
            pipeline,
            chunks,
            lnfact,
            adapter_name: info.name,
            label,
            lock: Mutex::new(()),
        })
    }

    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }
}

fn make_buffer(device: &wgpu::Device, label: &str, data: &[u8], usage: wgpu::BufferUsages) -> wgpu::Buffer {
    // Storage bindings may not be empty; pad to one word
    let padded;
    let data = if data.is_empty() {
        padded = [0u8; 4];
        &padded[..]
    } else {
        data
    };
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label),
        contents: data,
        usage,
    })
}

fn build_chunks(device: &wgpu::Device, table: &IonTable, max_binding: usize) -> Result<Vec<Chunk>, String> {
    let n = table.len();
    let max_ions = max_binding / 4;
    let max_peps = (MAX_WORKGROUPS * WORKGROUP) as usize;
    let mut chunks = Vec::new();
    let mut first = 0usize;
    while first < n {
        let mut end = first + 1;
        let ion_start = table.ion_off[first] as usize;
        if (table.ion_off[end] as usize) - ion_start > max_ions {
            return Err(format!("peptide {first} alone exceeds the GPU buffer limit"));
        }
        while end < n
            && end - first < max_peps
            && (table.ion_off[end + 1] as usize) - ion_start <= max_ions
            && (end + 1 - first + 1) * 4 <= max_binding
        {
            end += 1;
        }
        let (i0, i1) = (table.ion_off[first] as usize, table.ion_off[end] as usize);
        let (s0, s1) = (table.site_off[first] as usize, table.site_off[end] as usize);
        let st = wgpu::BufferUsages::STORAGE;
        let count = end - first;
        let out = |label| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: (count * 4) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            })
        };
        let read = |label| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: (count * 4) as u64,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            })
        };
        chunks.push(Chunk {
            first,
            count,
            ion_base: i0 as u32,
            site_base: s0 as u32,
            ion_off: make_buffer(device, "ion_off", bytes_of(&table.ion_off[first..=end]), st),
            ion_mass: make_buffer(device, "ion_mass", bytes_of(&table.ion_mass[i0..i1]), st),
            ion_meta: make_buffer(device, "ion_meta", bytes_of(&table.ion_meta[i0..i1]), st),
            site_off: make_buffer(device, "site_off", bytes_of(&table.site_off[first..=end]), st),
            sites: make_buffer(device, "sites", bytes_of(&table.sites[s0..s1]), st),
            pep_mass: make_buffer(device, "pep_mass", bytes_of(&table.peptide_mass[first..end]), st),
            out_score: out("out_score"),
            out_info: out("out_info"),
            read_score: read("read_score"),
            read_info: read("read_info"),
        });
        first = end;
    }
    Ok(chunks)
}

impl ChainKernel for GpuKernel {
    fn name(&self) -> &str {
        &self.label
    }

    fn survivors(&self, query: &ProcessedSpectrum, params: &KernelParams) -> Vec<Survivor> {
        let _guard = self.lock.lock().unwrap_or_else(|e| e.into_inner());
        let st = wgpu::BufferUsages::STORAGE;
        let peaks_mass = make_buffer(&self.device, "peaks_mass", bytes_of(&query.masses), st);
        let peaks_int = make_buffer(&self.device, "peaks_int", bytes_of(&query.intensities), st);

        // Fragment charges beyond the table (precursors above 33+) are not tried
        let mut charges = [[0f32; 4]; MAX_CHARGE_TERMS];
        let mut n_charges = 0;
        for (slot, (inv_z, slack_z)) in charges.iter_mut().zip(params.charge_terms()) {
            *slot = [inv_z, slack_z, 0.0, 0.0];
            n_charges += 1;
        }
        let (lo_factor, hi_factor) = params.fragment_tol.factors();

        let mut all = Vec::new();
        for chunk in &self.chunks {
            let gp = GpuParams {
                n_peaks: query.masses.len() as u32,
                n_pep: chunk.count as u32,
                max_charge: n_charges as u32 + 1,
                min_matched: params.min_matched as u32,
                ion_base: chunk.ion_base,
                site_base: chunk.site_base,
                pad0: 0,
                pad1: 0,
                precursor_mass: params.precursor_mass,
                pad2: 0.0,
                lo_factor,
                hi_factor,
                abs_lo: params.fragment_tol.abs_lo,
                abs_hi: params.fragment_tol.abs_hi,
                pad3: 0.0,
                pad4: 0.0,
                charges,
            };
            let uniform = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("params"),
                contents: bytes_of(std::slice::from_ref(&gp)),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            fn entry(binding: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
                wgpu::BindGroupEntry {
                    binding,
                    resource: buffer.as_entire_binding(),
                }
            }
            let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("crosslink-chain"),
                layout: &self.pipeline.get_bind_group_layout(0),
                entries: &[
                    entry(0, &uniform),
                    entry(1, &peaks_mass),
                    entry(2, &peaks_int),
                    entry(3, &chunk.ion_off),
                    entry(4, &chunk.ion_mass),
                    entry(5, &chunk.ion_meta),
                    entry(6, &chunk.site_off),
                    entry(7, &chunk.sites),
                    entry(8, &chunk.pep_mass),
                    entry(9, &self.lnfact),
                    entry(10, &chunk.out_score),
                    entry(11, &chunk.out_info),
                ],
            });

            let mut encoder = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                pass.set_pipeline(&self.pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups((chunk.count as u32).div_ceil(WORKGROUP), 1, 1);
            }
            let bytes = (chunk.count * 4) as u64;
            encoder.copy_buffer_to_buffer(&chunk.out_score, 0, &chunk.read_score, 0, bytes);
            encoder.copy_buffer_to_buffer(&chunk.out_info, 0, &chunk.read_info, 0, bytes);
            self.queue.submit(Some(encoder.finish()));

            let (tx, rx) = std::sync::mpsc::channel();
            for buffer in [&chunk.read_score, &chunk.read_info] {
                let tx = tx.clone();
                buffer.slice(..).map_async(wgpu::MapMode::Read, move |r| {
                    let _ = tx.send(r);
                });
            }
            let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
            for _ in 0..2 {
                rx.recv().expect("gpu map channel").expect("gpu buffer map failed");
            }
            {
                let score_view = chunk.read_score.slice(..).get_mapped_range().expect("mapped range");
                let info_view = chunk.read_info.slice(..).get_mapped_range().expect("mapped range");
                let scores: &[f32] = cast(&score_view);
                let infos: &[u32] = cast(&info_view);
                for p in 0..chunk.count {
                    if scores[p] >= 0.0 {
                        all.push(Survivor {
                            peptide: (chunk.first + p) as u32,
                            score: scores[p],
                            site: infos[p] >> 16,
                            matched: (infos[p] & 0xFFFF) as u16,
                        });
                    }
                }
            }
            chunk.read_score.unmap();
            chunk.read_info.unmap();
        }
        top_survivors(all, params.max_survivors)
    }
}

fn cast<T: Copy>(bytes: &[u8]) -> &[T] {
    // Safety: mapped ranges are 4-byte aligned and sized in whole words
    unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const T, bytes.len() / std::mem::size_of::<T>()) }
}

