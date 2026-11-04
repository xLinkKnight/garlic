use crate::{read_and_execute, tdf::BrukerProcessingConfig, Error};
use sage_core::spectrum::RawSpectrum;
use serde::Serialize;
use tokio::io::AsyncReadExt;
use url::Url;

#[derive(Debug, PartialEq, Eq)]
pub enum FileFormat {
    MzML,
    MGF,
    TDF,
    Unidentified,
}

impl FileFormat {
    /// Does this file format support parallel reading?
    /// By this I mean that there is 'within' file parallelism that
    /// would make it faster to read than reading mutiple files in
    /// parallel. (is giving 4 cores to read 1 file 4 times, faster
    /// than giving 1 cores to read 1 file and read 4 at the same time)
    pub fn within_file_parallel(&self) -> bool {
        match self {
            FileFormat::MzML => false,
            FileFormat::MGF => false,
            FileFormat::TDF => true,
            FileFormat::Unidentified => false,
        }
    }
}

impl From<&str> for FileFormat {
    fn from(s: &str) -> Self {
        let path_lower = s.to_lowercase();
        if path_lower.ends_with(".mgf.gz") || path_lower.ends_with(".mgf") {
            FileFormat::MGF
        } else if is_bruker(&path_lower) {
            FileFormat::TDF
        } else if path_lower.ends_with(".mzml.gz") || path_lower.ends_with(".mzml") {
            FileFormat::MzML
        } else {
            FileFormat::Unidentified
        }
    }
}

const BRUKER_EXTENSIONS: [&str; 5] = [".d", ".tdf", ".tdf_bin", "ms2", "raw"];

fn is_bruker(path: &str) -> bool {
    BRUKER_EXTENSIONS.iter().any(|ext| {
        if path.ends_with(std::path::MAIN_SEPARATOR) {
            path.strip_suffix(std::path::MAIN_SEPARATOR)
                .unwrap()
                .ends_with(ext)
        } else {
            path.ends_with(ext)
        }
    })
}

pub fn read_spectra(
    url: &Url,
    file_id: usize,
    sn: Option<u8>,
    bruker_processor: BrukerProcessingConfig,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    match FileFormat::from(url.as_ref()) {
        FileFormat::MzML => read_mzml(url, file_id, sn),
        FileFormat::MGF => read_mgf(url, file_id),
        FileFormat::TDF => read_tdf(url, file_id, bruker_processor, requires_ms1),
        FileFormat::Unidentified => panic!("Unable to get type for '{}'", url), // read_mzml(path, file_id, sn),
    }
}

/// Only the spectra of one file with these ids (Sage's spectrum ids), as
/// [`read_spectra`] reads them. Bruker data is read spectrum by spectrum,
/// without decoding the rest of the run; other formats are read whole and
/// filtered.
pub fn read_spectra_with_ids(
    url: &Url,
    file_id: usize,
    sn: Option<u8>,
    bruker_processor: BrukerProcessingConfig,
    ids: &std::collections::HashSet<String>,
) -> Result<Vec<RawSpectrum>, Error> {
    if FileFormat::from(url.as_ref()) == FileFormat::TDF {
        let numbers: Vec<usize> = ids.iter().filter_map(|id| id.parse().ok()).collect();
        if numbers.len() == ids.len() {
            if url.scheme() != "file" {
                return Err(Error::InvalidUri);
            }
            let path = url.to_file_path().map_err(|_| Error::InvalidUri)?;
            return crate::tdf::TdfReader
                .read_ids(&path, file_id, bruker_processor, &numbers)
                .map_err(Error::TDF);
        }
    }
    Ok(read_spectra(url, file_id, sn, bruker_processor, false)?
        .into_iter()
        .filter(|s| ids.contains(&s.id))
        .collect())
}

pub fn read_mzml(
    url: &Url,
    file_id: usize,
    signal_to_noise: Option<u8>,
) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(url, |bf| async move {
        Ok(crate::mzml::MzMLReader::with_file_id(file_id)
            .set_signal_to_noise(signal_to_noise)
            .parse(bf)
            .await?)
    })
}

pub fn read_tdf(
    url: &Url,
    file_id: usize,
    bruker_spectrum_processor: BrukerProcessingConfig,
    requires_ms1: bool,
) -> Result<Vec<RawSpectrum>, Error> {
    if url.scheme() != "file" {
        log::error!("Bruker files must be local: {}", url);
        return Err(Error::InvalidUri);
    }

    let path = url.to_file_path().map_err(|_| Error::InvalidUri)?;
    let res = crate::tdf::TdfReader.parse(&path, file_id, bruker_spectrum_processor, requires_ms1);
    match res {
        Ok(t) => Ok(t),
        Err(e) => Err(Error::TDF(e)),
    }
}

pub fn read_mgf(url: &Url, file_id: usize) -> Result<Vec<RawSpectrum>, Error> {
    read_and_execute(url, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        let res = crate::mgf::MgfReader::with_file_id(file_id).parse(contents);
        match res {
            Ok(m) => Ok(m),
            Err(e) => Err(Error::MGF(e)),
        }
    })
}

/// Position (0-based) of every spectrum in a spectrum file, by the id Sage gives
/// it: the mzML `spectrum` `id` → its `index` attribute, or the MGF `TITLE` →
/// its order in the file. Viewers such as xiVIEW (MGF `index=N`) and ProXL
/// (mzML scan numbers without `scan=`) locate spectra by position.
pub fn read_spectrum_positions(url: &Url) -> Result<std::collections::HashMap<String, usize>, Error> {
    use quick_xml::events::Event;
    use tokio::io::AsyncBufReadExt;
    let lower = url.path().to_ascii_lowercase();
    let mgf = lower.ends_with(".mgf") || lower.ends_with(".mgf.gz");
    read_and_execute(url, |reader| async move {
        let mut positions = std::collections::HashMap::new();
        if mgf {
            let mut lines = reader.lines();
            let mut order = 0usize;
            let mut in_spectrum = false;
            while let Some(line) = lines.next_line().await? {
                let line = line.trim();
                if line == "BEGIN IONS" {
                    in_spectrum = true;
                } else if line == "END IONS" {
                    in_spectrum = false;
                    order += 1;
                } else if let Some(title) = line.strip_prefix("TITLE=").filter(|_| in_spectrum) {
                    positions.entry(title.to_string()).or_insert(order);
                }
            }
            return Ok(positions);
        }
        let mut xml = quick_xml::Reader::from_reader(reader);
        let mut buf = Vec::new();
        loop {
            match xml.read_event_into_async(&mut buf).await {
                Ok(Event::Start(ref ev)) | Ok(Event::Empty(ref ev)) if ev.name().into_inner() == b"spectrum" => {
                    let attr = |key: &[u8]| {
                        ev.try_get_attribute(key)
                            .ok()
                            .flatten()
                            .and_then(|a| String::from_utf8(a.value.into_owned()).ok())
                    };
                    if let (Some(id), Some(index)) = (attr(b"id"), attr(b"index").and_then(|i| i.parse().ok())) {
                        positions.insert(id, index);
                    }
                }
                Ok(Event::Eof) => break,
                Err(e) => return Err(Error::MzML(crate::mzml::MzMLError::XMLError(e))),
                _ => {}
            }
            buf.clear();
        }
        Ok(positions)
    })
}

pub fn read_fasta<S>(
    url: &Url,
    decoy_tag: S,
    generate_decoys: bool,
) -> Result<sage_core::fasta::Fasta, Error>
where
    S: AsRef<str>,
{
    read_and_execute(url, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents)
            .await
            .map_err(crate::Error::IO)?;
        Ok(sage_core::fasta::Fasta::parse(
            contents,
            decoy_tag.as_ref(),
            generate_decoys,
        ))
    })
}

pub fn read_json<S, T>(path: S) -> Result<T, Error>
where
    S: AsRef<str>,
    T: for<'de> serde::Deserialize<'de>,
{
    read_and_execute(path, |mut bf| async move {
        let mut contents = String::new();
        bf.read_to_string(&mut contents).await?;
        Ok(serde_json::from_str(&contents)?)
    })
}

/// Send telemetry data
pub fn send_data<T>(url: &str, data: &T) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
where
    T: Serialize,
{
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    rt.block_on(async {
        let client = reqwest::ClientBuilder::default().https_only(true).build()?;
        let res = client.post(url).json(data).send().await?;
        res.error_for_status()?;
        Ok(())
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn spectrum_positions_by_id() {
        let mzml = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/LQSRPAAPPAPGPGQLTLR.mzML");
        let url = crate::to_url(mzml.to_str().unwrap()).unwrap();
        let positions = read_spectrum_positions(&url).unwrap();
        assert_eq!(positions.get("controllerType=0 controllerNumber=1 scan=30069"), Some(&0));

        // MGF: titles in file order, also gzipped
        let dir = std::env::temp_dir().join(format!("sage-positions-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mgf = "BEGIN IONS\nTITLE=first\nPEPMASS=500\n100 1\nEND IONS\n\
                   BEGIN IONS\nTITLE=second spectrum\nPEPMASS=600\n100 1\nEND IONS\n";
        let plain = dir.join("x.mgf");
        std::fs::write(&plain, mgf).unwrap();
        let p = read_spectrum_positions(&crate::to_url(plain.to_str().unwrap()).unwrap()).unwrap();
        assert_eq!((p["first"], p["second spectrum"]), (0, 1));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_identify_format() {
        assert_eq!(FileFormat::from("foo.mzml"), FileFormat::MzML);
        assert_eq!(FileFormat::from("foo.mzML"), FileFormat::MzML);
        assert_eq!(FileFormat::from("foo.mgf"), FileFormat::MGF);
        assert_eq!(FileFormat::from("foo.mgf.gz"), FileFormat::MGF);
        assert_eq!(FileFormat::from("foo.tdf"), FileFormat::TDF);
        assert_eq!(FileFormat::from("./tomato/foo.d"), FileFormat::TDF);
        assert_eq!(FileFormat::from("./tomato/foo.d/"), FileFormat::TDF);
    }
}
