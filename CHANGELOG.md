# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Fixed
- CI: the Rust workflow runs its test searches with `cargo run --bin sage`. A bare `cargo run` had failed on every push since `sage-gui` joined the workspace, as Cargo could not choose between the two programs.

### Removed
- The Docker image workflow (`.github/workflows/docker.yml`). It never published an image (the workflow token cannot write packages) and failed on every release; XNova releases are Windows zips.

## [v1.0.1] - 2026-11-06

### Added
- `sage.exe` carries the XNova icon too (Explorer, shortcuts), like `sage-gui.exe`.

### Changed
- XNova branding: the application icon (window, taskbar and `sage-gui.exe`) is the XNova icon, and the README shows the XNova logo, with a new title, introduction and list of major features.
- Renamed from garlic to XNova throughout: links point to https://github.com/xLinkKnight/XNova; the development variables are `XNOVA_XIVIEW_DIR`, `XNOVA_XIVIEW_RESULTS` and `XNOVA_XIVIEW_MINUTES` (were `GARLIC_…`); xiVIEW shows the scores as `XNova:hyperscore` and `XNova:q-value`. The `xiview` folder is rebuilt to match (it loads multivalent link metadata from `get_xnova_link_metadata`), so use it with this sage-gui, not an older one.
- `sage-gui`: the theme's blue accents are now the orange of the XNova logo: section titles, links, selected tabs, the Run search button, hover and selection tints, the text cursor and the running status. Orange (#e8650a) on the light theme and orange-gold (#f59e1f) on the dark one; borders and greys are warmer to match. Data colours (ion series, linker stubs, mobiligram classes) are unchanged.

## [v1.0.0] - 2026-11-05

### Added
- `sage-gui`: xiVIEW on this computer. **Open in xiVIEW** (bottom bar, Progress tab, Annotated spectra tab) and `sage-gui --xiview <results folder>` show a crosslink search's results in [xiVIEW](https://github.com/Rappsilber-Laboratory/xiview) in the web browser, with nothing uploaded and no database: sage-gui serves the xiVIEW web app (the `xiview` folder next to sage-gui.exe) on 127.0.0.1 and answers xiVIEW's data API (the PRIDE crosslinking API's `get_xiview_*` endpoints) from the results folder. Crosslinks within XNova's crosslink FDR (1%, 5% or 10%; xiVIEW's own decoy FDR is not used); proteins by gene name. Multivalent crosslinks: each link a crosslink, drawn dash-dotted in the network and circle views, with link metadata (`Multivalent`, `Multivalent species`; colour by species), and the whole species in the spectrum viewer (every peptide and link; fragments named by the peptides they carry, `b14+2P`). Annotated spectra for xiVIEW's spectrum viewer (xiSPEC) from Annotator's engine and XNova's MS-cleavable stub ions, converted to xiSPEC's format; xiSPEC's settings re-annotate locally. Spectra read as Sage read them (timsTOF runs kept open: a spectrum opens in a fraction of a second); a spectrum that cannot be fetched (sage-gui closed) or annotated says so. xiVIEW, xiNET and xiSPEC are the work of the Rappsilber Laboratory: README (Citation, Acknowledgments) and CITATION.cff cite them (Combe et al. 2024, Kolbowski et al. 2018, Combe et al. 2015), the repositories (xiview, xiview-server, crosslinking-api and mzidentml-reader with the PRIDE team at EMBL-EBI, xiAnnotator) and PRIDE (Perez-Riverol et al. 2025). `xiview/build_xiview.py` builds xiVIEW (GPL-3.0) from source at a pinned commit with the changes in `xiview/source_fixes.py` (multivalent display; the spectrum settings' request id; loop-link spectra, which did not open). See CROSSLINKING.md.

## [v0.23.0] - 2026-11-04

### Added
- Two-database crosslink searches, for complex crosslink searches (`database.crosslink_databases` with `database_1` and `database_2`, replacing `fasta`; `sage-gui`: Files tab → Two-database crosslink search, which greys out the FASTA field). Linear peptides are searched in both files as one database; crosslinks only between a Database 1 and a Database 2 protein (cpu, reflect, gpu and MS-cleavable searches), so pairs within one database are never scored or counted in the FDR. Loop-links are not searched; a multivalent species needs at least one link between the databases. A peptide of proteins in both files counts as either; decoys follow their target. `crosslinks.sage.tsv` gets `database_a` / `database_b` columns (`1`, `2`, `1+2`). The viewer, exports and labelled-protein list read both files. See CROSSLINKING.md.
- Uniformly isotope-labelled proteins (`database.labelled_proteins`; `sage-gui`: Database & enzyme tab → Isotope-labelled proteins: the chosen proteins, and a searchable list of the FASTA's proteins, closed at first): chosen proteins, expressed in ¹⁵N, ¹³C or ¹³C¹⁵N media, are searched labelled in linear and crosslink searches. Every residue carries its carbon and/or nitrogen atoms' mass shift, added after the fixed and variable modifications, so labels do not count towards `max_variable_mods` and a residue keeps its modifications (a modification's own atoms stay light; deamidation and pyro-Glu from Q lose a ¹⁵N). `include_unlabelled` (default on) also searches them unlabelled. Peptides shared with other proteins are searched both ways, each credited to its proteins. A labelled residue still crosslinks, and crosslinks between labelled and unlabelled proteins are found. `results.sage.tsv` gets an `isotope_label` column, `crosslinks.sage.tsv` `isotope_label_a` / `isotope_label_b`. Combines with SILAC.

## [v0.22.0] - 2026-11-04

### Added
- SILAC, duplex or triplex, in linear and crosslink searches (`database.silac`: labelled channels with a label mass per residue; `sage-gui`: Scoring & output tab → Quantification → SILAC, with Lys4/Lys6/Lys8/Arg6/Arg10 presets). Every peptide is searched light and in each channel with all its labelled residues labelled (no mixed states; 2–3× the search). Labels are modifications (`UNIMOD:259` Lys8, `UNIMOD:267` Arg10, `UNIMOD:188` Lys6/Arg6, `UNIMOD:481` Lys4), so fragments, exports and the viewer carry them; a label on a modified residue is one combined modification. Crosslinks: a labelled lysine still links, and a crosslink (or multivalent species) joins peptides of one channel unless `mixed_crosslinks` is set, in every backend. Quantification: each identified peptide and crosslinked species is traced in MS1 in all channels at once and integrated over one shared peak (identified in one channel, quantified in all): `silac.sage.tsv` (intensities and ratios to light per run), `silac_proteins.sage.tsv` (median ratio of unique peptides) and `crosslinks.silac.tsv`; `results.sage.tsv` and `crosslinks.sage.tsv` get a `silac_channel` column. Tested on synthetic data and on an unlabelled run as a negative control; not yet on a real SILAC sample. See SILAC.md.

## [v0.21.0] - 2026-11-04

### Added
- `sage-gui`: in the Annotated spectra peptide map, each crosslinked peptide's label names its protein and where it is linked: `peptide 1 · RPO21` (the gene, else the UniProt entry name), the protein's name, and the link residues' positions in the protein (`K1246`; `T3, K12` for a peptide with two links; both residues of a loop-link). Read from the search's FASTA (UniProt headers: name up to `OS=`, gene `GN=`). A peptide in more than one protein shows `+ n`; hover the label for every protein, the residues the peptide spans and its links there. Decoys show their target's positions. Multivalent species' peptides that are not in the protein columns are looked up in the FASTA.

### Changed
- `sage-gui`: a peptide map wider than the window (long peptides, protein labels) scrolls sideways instead of being cut off.

## [v0.20.1] - 2026-11-04

### Changed
- `sage-gui`: the Annotated spectra tab opens the timsTOF runs of the shown matches in the background when a results folder is loaded (two at a time, the best match's run first; about a second each, and 30 MB each while open), so the first match picked from a run no longer waits for its run to open. The tab shows "opening runs 5 of 21" meanwhile; a run picked before its turn is opened right away.

### Fixed
- `sage-gui`: in the Annotated spectra tab, some matches stayed on "Reading the spectrum…" (in v0.20.0, e.g. multivalent species of another timsTOF run than the first one shown). Rows in view are read ahead only from runs already open, and a match whose row had been in view was then never asked for when picked, so its run was not opened. A picked match is now always read first, opening its run (about a second the first time). Annotator's ontologies also load in the background when a results folder is opened, instead of on the first match shown.

## [v0.20.0] - 2026-11-04

### Added
- `sage-gui`: the GitHub mark with the program's version (e.g. v0.20.0) next to the light/dark switch, a link to the project on GitHub (https://github.com/xLinkKnight/XNova); it takes the theme's text colour.
- `sage-gui`: the Annotated spectra map of multivalent crosslinks shows every link: each peptide as a row with its link residues shaded, a bridge per link labelled `XL1`, `XL2`, … as in the ProForma (linked peptides stacked next to each other where possible; a link skipping a row runs along a lane). MS-cleavable linkers (DSSO, DSBU, cliXlink): the stub ions of the cleaved species as the search scores them: each peptide's peptide + stubs peaks (doublets; `P+A+S` triplets with two links), fragment ions with one stub per link residue held (`y17+A+S²+`, hooks in each stub's colour), and the parts held together when only one linker broke (`{2,3}+A`), as triangles on that linker's bridge. Also in exported figures.

### Changed
- `sage-gui`: the Annotated spectra tab no longer waits for whole runs. It used to read every spectrum of every searched file, one file after another, before showing any (about 7.5 s per timsTOF run, minutes for a 21-run search). Spectra are now read in the background as needed, the picked match first: timsTOF spectra one by one by their id (0.4 ms each, after about a second to open a run; identical peaks to the whole-run read), mzML / MGF files whole once each. Results load in the background as well.
- `sage-gui`: the mobiligrams of a finished search are written in the background (the report shows at once, the files are listed when ready), several at a time, with the system fonts loaded once instead of for every PDF.

## [v0.19.0] - 2026-11-03

### Added
- Crosslinking: multivalent crosslinks (`"multivalent": 2`–`4`, non-cleavable and MS-cleavable linkers): species held together by two up to that many linkers, with several links between two peptides (`a2b`, `a2a`), chains (`a1b1c`, `a1a1b`) and branches (`a1b,b1c,b1d`) of up to five peptides. Non-cleavable: each peptide is scored with the groups hanging off its residues (fragments tied through a group are not counted) and needs its own linker-carrying ions. MS-cleavable (DSSO, DSBU, cliXlink): every linker breaks, each peptide keeps a stub per link residue and needs its peptide + stubs peaks (doublets, triplets); the arrangement is told by parts freed when a single linker breaks (`b1c` + stub), reported as `topology_peaks`; multivalent species have their own FDR group (each species with a decoy peptide counts as 1 / (2^n − 1) false targets, n = distinct peptides; (TD + DD) / 3 for pairs). `crosslinks.sage.tsv`: `link_type` `multivalent` and new `topology`, `species` and `topology_peaks` columns. The Annotator CSV and the viewer include them (ProForma with `#XL1`, `#XL2`, …); xiVIEW, ProXL, mzIdentML and crosslink LFQ leave them out. `sage-gui`: "search multivalent crosslinks" with the number of linkers, and a warning that the search space is large and meant for small protein sets.
- `sage-gui`: multivalent crosslinks in the mobiligram, in green (#6FB120).
- `sage-gui`: a finished search on timsTOF (Bruker `.d`) data writes mobiligrams into a `Mobiligrams` folder inside the results folder, as SVG and PDF (all precursors, white background, titled): one per MS file, named after it with `_MG` (`Bruker.d` → `Bruker_MG.svg`), and, when two or more MS files were searched, `Combined_Mobiligram` with all of them together (a single MS file gets only its `_MG` figure). They are listed with the outputs in the Progress tab, and the viewer's Save… suggests the same names. MS files without ion mobility get none.

### Changed
- `sage-gui`: number fields are typed in, or stepped with stacked + / − buttons to their right; dragging no longer changes them (the pointer shows a text cursor). Values keep their decimals when stepped (15.9949 + 1 = 16.9949) and show without trailing zeros.
- `sage-gui`: titles (tabs, section headings, column headers, the viewer heading and match title) are bold (Segoe UI Bold on Windows).
- `sage-gui`: in the Annotated spectra tab the border between the match list and the spectrum can be dragged to give either side more room (double-click resets it); the list keeps at least 220 px and the spectrum 420 px. Matches stay on one line, long ones ending in "…", shown whole on hover.
- GPU support is built into `sage` by default (the `gpu` feature of `sage-cli` is a default feature): a plain `cargo build --release -p sage-cli` includes it, and test or other builds no longer replace `sage.exe` with one without it. `--no-default-features` builds without it. Searches still run on the CPU unless the crosslink backend is `gpu` or `auto`.

## [v0.18.0] - 2026-11-03

### Added
- Peptide- and protein-level reports from every linear search: `peptides.sage.tsv` (one row per target peptide at 1% peptide FDR: proteins and groups, charges, PSM count per file, best spectrum/peptide/protein/group q-values and scores, summed MS2 intensity and, with LFQ, MS1 intensity per file) and `proteins.sage.tsv` (one row per target protein group with a unique peptide in `peptides.sage.tsv`: its group q-value, unique and shared peptides, PSM count per file and, with LFQ, the summed intensity of its unique peptides per file), next to `results.sage.tsv`. Protein groups are listed with their q-value rather than cut at 1%: Sage's picked group FDR cannot go below 1 / (number of protein groups), about 0.08 for a 13-protein complex, so filter on `protein_group_q` for proteome-wide databases.
- `sage-gui`: mobiligram of timsTOF (Bruker `.d`) runs ("Mobiligram…" in the Annotated spectra tab): every searched precursor by m/z and ion mobility (1/K0), marker shape by charge state (1+ to 6+ and up), colour by identification at 1% FDR (crosslinks between proteins red #FF523B, within a protein orange #FF9D35, loop-links pink; linear peptides blue at 50% transparency; unidentified grey), with a legend of colours (with counts) and charge shapes. Interactive in the viewer; exported as SVG or PDF, white or transparent background.
- `sage-gui`: PDF (vector, fonts embedded) as an export format for annotated spectra too.
- `crosslinks.sage.tsv`: `ion_mobility` column (precursor 1/K0, timsTOF data; 0 otherwise), after `rt`.
- Crosslinking: loop-links of MS-cleavable linkers (DSSO, DSBU, cliXlink). The opened ring shows as fragment ions that carry the stub of one linked residue (ring ions); peptides whose mass plus the linker matches the precursor are scored with these ions, the ions spanning both residues (both stubs) and the plain ones outside, and need `min_linked_peaks` matched ring ions. On the DSSO Pol II run: 154 loop-link CSMs at 1% FDR (1 decoy), with 12.4 ring ions on average against 4.2 for decoys; spectra of a PRP40 loop-link that had been the best-scoring decoy crosslinks are now explained, and intra-protein crosslinks at 1% FDR rose from 74 to 148. `sage-gui` marks the ring ions as stub hooks in the peptide map and labels them (`b4+A`).
- `sage-gui`: the mobiligram legend always lists every class, with n = 0 when none was found.

### Changed
- Crosslinking: the linear peptide search (`linear_search`) is on by default, so crosslink searches also give the usual bottom-up results (`results.sage.tsv`, `peptides.sage.tsv`, `proteins.sage.tsv`). Set `"linear_search": false` (or untick it in `sage-gui`) to skip it.
- `sage-gui`: the HTML report (`results.sage.report.html`) is on by default for new settings.

## [v0.17.0] - 2026-11-03

### Added
- `sage-gui`: built-in annotated spectrum viewer, the "Annotated spectra" tab, using [Annotator](https://github.com/snijderlab/annotator)'s engine (`mzcore` / `mzannotate`) with no separate program. Lists the linear and crosslink matches of a results folder (filters for kind, targets, spectrum q-value, text search); for the one picked, draws its spectrum with ion-coloured, labelled peaks (zoom, hover for m/z, error and every interpretation), the backbone coverage of each peptide, and a table of matched fragments. Crosslinks and loop-links are ProForma 2.0 crosslinked peptidoforms; MS-cleavable linkers show their stub ions. Spectra are read with Sage's own reader by Sage's spectrum ids (mzML, MGF, Bruker `.d`: the precursor's merged spectrum, as searched). "View spectra" buttons (bottom bar, Progress tab) and `sage-gui --view <results folder>` open it. Crosslinked pairs are drawn joined by a bridge between their linked residues (loop-links by a bracket); for MS-cleavable linkers, ions carrying a stub are found per peptide and shown in the stub's colour: peptide + stub doublets as triangles on the bridge pointing at their peptide, stub-carrying fragment ions as hooks with the N-/C-terminal directionality of the plain series, and `+A` / `+S` / `+T` in the peak labels. The header shows retention time, mass error and, for timsTOF data, the precursor ion mobility (1/K0). Publication figures: "Export figure" saves the annotated spectrum (current zoom), peptide map with bridge and stubs, legend and an optional title (with RT and ion mobility) as SVG (vector) or PNG / TIFF at 300, 600 or 1200 dpi, 7 in wide, on a white or transparent background.
- Files for the separate Annotator app: `"write_annotator": true` (or `--write-annotator`) writes `results.annotator.csv` (linear PSMs) and `crosslinks.annotator.csv` (crosslink-spectrum matches) at 1% FDR in Annotator's basic CSV format, with ProForma 2.0 crosslinks (`A[XLMOD:02126#XL1]..//B[#XL1]..`) and spectrum positions for mzML, MGF and Bruker `.d` (DDA-PASEF, from `analysis.tdf`). Off by default. See `CROSSLINKING.md`.
- `sage-core`: `CrosslinkSettings::proforma_linker` and `Modification::proforma_label`, the ProForma of a crosslinker and a modification (shared by the export and the viewer).
- Citation and acknowledgments for Annotator (Schulte et al., *Anal. Chem.* 2025, doi:10.1021/acs.analchem.5c02832) and its `mzcore` / `mzannotate` libraries in `README.md` and `CITATION.cff`.

### Changed
- The Windows release `sage.exe` is built with `--features gpu`, so the crosslink `gpu` backend works out of the box (DirectX 12).

### Fixed
- `write_pin`, `write_report` and `annotate_matches` in the configuration file were ignored: the command-line flags, absent, overrode them with `false`. A flag now only turns its output on. This also made the corresponding `sage-gui` checkboxes ineffective.

## [v0.16.0] - 2026-11-02

### Added
- Crosslinking: asymmetric MS-cleavable linkers. A cleavable linker can give its second arm its own stubs with `"cleavable": {..., "other": {"short": .., "long": .., "alt": [..]}}`; the doublet search then tries each peptide pair both ways round. New `clixlink` preset (cliXlink, 267.05653 Da; arm 1 54.01056 / 103.99320 / 85.98264, alkyne-tagged arm 163.06333 / 213.04596 / 195.03540).
- Crosslinking: `crosslinks_FDRfiltered.sage.tsv`, the accepted results only (target-target matches with `spectrum_q` <= 0.01, the CSM count in the log), next to the full `crosslinks.sage.tsv`.
- Crosslinking: `"linear_search": true` also runs a linear peptide search on the same spectra first and writes `results.sage.tsv`.
- Crosslinking: `min_linked_peaks` (default 2), the matched fragment ions that carry the linker required on each chain of a non-cleavable match, reported in the new `linked_a` / `linked_b` columns.
- Crosslinking: exports for viewers: `crosslinks_xiView.csv` (xiVIEW), `crosslinks_proXL.xml` (ProXL) and `crosslinks.mzid` (mzIdentML 1.3 with the crosslinking extension, for xiVIEW and PRIDE), with references to the spectrum files.
- Crosslinking: the `gpu` backend can be built on Windows (`--features gpu`, DirectX 12).
- `sage-gui`: precursor charge range (Database & enzyme tab); editable stub masses of cleavable linkers with their doublet spacing and mass balance, including a second arm for asymmetric linkers and the cliXlink preset; "also run a linear peptide search" and "min linked ions per chain" options.
- `sage-gui`: UniMod accession lookup on the Modifications tab (fills in name, mass and sites from unimod.org).
- `sage-gui`: Resources section (CPU threads, RAM limit, GPU name, VRAM and backend support).
- `sage-gui`: Progress tab while a search runs: status colour that pulses while running, step list, log with a draggable border, two-column report and an explanation of common failures; CPU, RAM, GPU and VRAM usage next to Run search.
- `sage-gui`: New folder button when choosing the results folder; new theme, window sized to fit the tabs, and an app icon.
- Crosslinking (XL-MS) search: enable with a `crosslink` config block; writes `crosslinks.sage.tsv`. Supports non-cleavable (DSS/BS3/DSG/BS2G, incl. loop-links) and MS-cleavable (DSSO/DSBU, doublet-based) linkers. See `CROSSLINKING.md`.
- Crosslinking: heterobifunctional chemistry with `site1` / `site2` (a link always joins a site1 group to a site2 group), and the zero-length reagents `edc` and `dmtmm` (amine to carboxyl, -18.010565 Da).
- Multiple enzymes (co-digestion): `"enzyme": {"enzymes": ["trypsin", "asp-n"]}` cleaves wherever any listed enzyme does. Named proteases: trypsin, trypsin/p, lys-c, lys-n, arg-c, asp-n, glu-c, chymotrypsin, chymotrypsin-high, proalanase; custom rules `{"cleave_at": "DE", "restrict": "P", "c_terminal": true}` can be mixed in.
- Crosslinking: quantification. `quant.tmt` adds reporter-ion columns to `crosslinks.sage.tsv` (MS2 or SPS-MS3); `quant.lfq` writes `crosslinks.lfq.tsv` with MS1 intensities per crosslinked peptide pair across runs (alignment on shared crosslinks, decoy-based peak q-values). Sage warns when a static modification (e.g. TMT on K) blocks crosslink-reactive residues.
- `sage-gui`: a desktop window for editing search settings (with hover help on every field), choosing FASTA/mzML/MGF/Bruker `.d`/Thermo `.raw` files (RAW converted with msconvert or ThermoRawFileParser) and running searches.

### Changed
- `sage-gui`: modification sites are always shown (the "show hidden sites" switch is gone).

### Fixed
- Crosslinking: non-cleavable searches reported false crosslinks when the wrong linker was searched. A loop-link was accepted on the unmodified ions outside its ring alone, so a correct peptide carrying any modification of about the linker's mass passed; a strong peptide paired with a barely matched partner passed as a crosslink. Each chain now needs matched linker-carrying ions (`min_linked_peaks`). On a DSSO-crosslinked Pol II run searched as DSS, CSMs at 1% FDR fell from 55 to 13 (4 with isotope error 0, against 155 for the correct DSSO search).
- Crosslinking: DSSO stub names were swapped in the documentation and the GUI: 103.99320 Da is the sulfenic acid, 85.98264 Da the thiol.
- Crosslinking: the BS2G preset used 114.03170 Da (the hydrolysed mono-link); the crosslink adds 96.02113 Da, like DSG.

## [v0.15.0]
### Added
- Support for named modifications using `{"mass": 15.9949, "name": "UniMod:35"}`
- Radical z ions for ETD/ECD: set `database.ion_kinds` to `["c", "z_radical"]`. The existing `"z"` option retains the non-radical ion type.
- IDPicker-based protein grouping with picked group FDR control (`protein_grouping` setting, enabled by default). Proteins are grouped using a bipartite graph greedy set cover approach, and protein group-level q-values are reported via target-decoy competition. New output columns: `protein_groups`, `num_protein_groups`, `protein_group_q`.
- `protein_grouping_peptide_fdr` parameter to control the peptide FDR threshold used for confident peptides during protein grouping (default: 0.01)
- Initial support for LFQ on data with ion mobility.
- Speedup on the generation of databases when large number of peptides are redundant.
- Database splitting and prefiltering to reduce memory usage (`database.prefilter`, `database.prefilter_chunk_size`, and `database.prefilter_low_memory` settings)
- Per-modification occurrence limits using `{"mass": <mass>, "max_count": <limit>}` entries in `database.variable_mods`; existing bare-mass entries remain supported.
- `database.max_combinations` to cap the number of peptide variants (including the unmodified form) generated from variable modifications, preferring variants with fewer modifications.
- Initial support for searching diaPASEF data
- `override_precursor_charge` setting that forces multiple charge states to be searched
- Cross-cloud storage support: Replaced AWS SDK with the `object_store` crate. Sage now natively supports reading/writing from Amazon S3, Google Cloud Storage, and Azure Blob Storage using `s3://`, `gs://`, and `az://` URL schemes.
- HTML QC report generation (`--write-report`)
- `lfq_settings.peptide_q_value` parameter for controlling which peptides are quantified
- Stack size configuration for Rayon threads
- Allow zero or multiple amino-acids as cleavage restrictions
- CITATION.cff file

### Fixed
- Use neutral NH3 mass in c, x, and non-radical z fragment calculations.
- Handle negative mass errors for .pin files
- Extract precursor m/z from 'isolation window target m/z' if missing
- Selected ion m/z of 0.0 was overwriting precursor.mz
- C/N-term mixup in modification handling
- Bruker `.tdf` filename handling (use parent directory name)
- Performance optimizations on prefiltering
- Picked protein FDR now correctly uses only proteotypic peptides for competition
- Fixed bug in MS1 picked-precursor FDR where target and decoy q-values were not distinguished when assigning them to peaks
- Fixed bug in candidate score ordering that could retain lower-scoring peptide matches during prefiltering
- Track charges on deisotoped ions to correctly report fragment charges and m/z in matched fragment annotations

### Breaking Changes
- Variable modifications must be specified as lists (e.g. `[15.9949]`), rather than single values (e.g. `15.9949`)
- `precursor_ppm` field reports the non-absoluted average mass error, rather than the absoluted average mass error.
- Don't deisotope reporter ion regions if MS2-based TMT/iTRAQ is used
- Removed `fragment_min_mz` and `fragment_max_mz` parameters. These were decreasing the accuracy of preliminary scoring estimation when attempting to annotate multiply-charged, high-m/z ions.
- `sage-cloudpath` no longer exposes the `CloudPath` type. All paths are represented as URLs (`url::Url`). Local paths are converted to `file://` URLs internally.

## [v0.14.7]
### Added
- Added columns missing from parquet output: `semi_enzymatic` and `missed_cleavages`
### Changed
- Fixed ion mobility parsing from some mzMLs
- MGF paths were being lowercased prior to parsing

## [v0.14.6]
### Added
- Support for MGF files
- Support for writing ion mobility measurements to output files: `ion_mobility`, `predicted_mobility`, `delta_mobility` added to primary tsv and parquet reports. Ion mobility is predicted in a similar manner to RT, using a linear model trained on the data from the search.

## [v0.14.5]
### Added
- Support for semi-enzymatic digests (`database.enzyme.semi_enzymatic` parameter)
- Ability to directly export matched fragment ions (e.g. for spectral library or rescoring) with the `--annotate-matches` CLI option. This is compatible with the `--parquet` CLI option as well. Annotations will be written to `matched_fragments.sage.tsv` or `matched_fragments.sage.parquet`
- Sage sends basic telemetry data (version of Sage, run time, OS, # of CPU cores, # of peptides in database, whether LFQ is used) to a remote server. No information about your actual data is sent - e.g. identifications, quantities, organism, or modifications are NOT tracked or reported.  This data will be used to help focus efforts on improving Sage and figuring which features are most used. Please take a look at `crates/sage-cli/src/telemetry.rs` to see exactly what is sent! You can disable sending telemetry data  by using the `--disable-telemetry-i-dont-want-to-improve-sage` CLI flag.
### Changed
- Modified visibility on some crate internals to support the [sagepy project](https://github.com/theGreatHerrLebert/sagepy)
- Added `psm_id` field to various output files to match the new `--annotate-matches` option.
### Removed
- Removed the `ms1_intensity` field from CSV output, since it is essentially useless


## [v0.14.4]
### Added
- **Unstable feature**: Preliminary support for reading Bruker .d folders (ddaPASEF; no MS1/LFQ support yet)
### Changed
- Retention times are converted to minutes
### Fixed
- Fixed bug where charge state 1 would never be searched

## [v0.14.3]
### Fixed
- Hotfix for bug in parquet LFQ writer

## [v0.14.2]
### Added
- `quant.lfq_settings.combine_charge_state` boolean option. By default this is set to `true`, and LFQ is performed on the peptide-level, where all charge states are treated as the same precursor. Setting this to `false` performs LFQ on the peptide-charge-level, where each charge state will be treated separately.
### Changed
- Percolator output format now contains the integer-valued charge state encoded in the `z=other` column, if the charge state is outside the range 2-6 (e.g. a value of 7 will appear in the `z=other` column, rather than it being one-hot encoded)
- LFQ uses the the charge state range from the `precursor_charge` configuration option for tracing MS1 peaks

## [v0.14.1]
### Added
- Added additional output showing search progress if `SAGE_LOG=trace` environment variable is set
- Added additional warnings about precursor tolerances
- Added configuration option `precursor_charge` to make it explicit what charge states are being searched in the case where the mzML does not contain charge state information, or where `wide_window` is turned on.
### Changed
- Added a warning message if variable modifications are specified as single values (e.g. `15.9949`) instead of lists of values (e.g. `[15.9949]`). By v0.15 this will become a hard error and will not parse, to simply some of the internal logic.

## [v0.14.0]
### Added
- Support for parquet file format output. Search results and reporter ion quantification will be written to one file (`results.sage.parquet`) and label-free quant will be written to another (`lfq.parquet`). Parquet files tend to be significantly smaller than TSV files, faster to parse, and are compatible with a variety of distributed SQL engines.
### Changed
- Implement heapselect algorithm for faster sorting of candidate matches (#80). This is a backwards-incompatible change with respect to output - small changes in PSM ranks will be present between v0.13.4 and v0.14.0

## [v0.13.4]
### Fixed
- Bug in mzML parser, where some older specification-compliant mzMLs would not parse. If your mzMLs previously parsed, then there will be no change in behavior. Added a test case

## [v0.13.3]
### Fixed
- Bug in `database.enzyme.restrict` parameter, where `null` values were being overriden with "P" (causing Trypsin/P to behave like Trypsin)

## [v0.13.2]
### Changed
- Subtle change to TMT integration tolerance, and selection of which ion to quantify (most intense). As a result, TMT integration should be more in agreement (if not 100% so) with ProteomeDiscover/FragPipe/etc
- Remove `delta_mass` (precursor ppm) LDA feature - instead, build a delta mass (or ppm) profile using KDE/posterior error calculation code, and use the P(decoy) as a feature for LDA.

## [v0.13.1]
### Changed
- Internal performance and stability improvements for RT prediction & LDA

## [v0.13.0]
### Added
- Better error reporting thanks to @Elendol
- Added support for multiple variable mods for the same amino acid
- Added support for N/C-terminal modifications specific to an individual amino acid

New syntax:
```json
"variable_mods": {
    "M": [15.9949],
    "^Q": -17.026549,
    "^E": -18.010565,
    "[": 42.010565
}
```

Either a single floating point number (-18.0) or a list of floating point numbers ([-18.0, -15.2]) can be supplied as modifications. Support for single values may eventually be phased out to simplify the parser.

### Changed
- Changed "_fdr" columns to "_q" (e.g. "spectrum_q") in "results.sage.tsv" file
- Changed internal data representation of `Peptide` struct to allow for sharing of sequences (using `Arc`) among modified peptides
- Fragment index creation should now be faster

## [0.12.0]
### Added
- Add `wide_window` option to configuration file. This option turns off `precursor_tol`, instead using the isolation window written in the mzML file.
### Changed
- Changed internal calculation of precursor tolerances when searching with `isotope_errors`. The new version should be more accurate. This change also enables a significant boost to search speed for open searches.

## [0.11.2]
### Added
- Add rank & charge features to LDA
### Changed
- One-hot encode charge state information for percolator `.pin` files
- Change PSMId -> SpecId for Mokapot compatibility with `.pin` files

## [0.11.1]
### Added
- Support for additional fragment ion types, via the "database.ion_kinds" configuration option. Valid values are "a", "b", "c", "x", "y", "z"
### Changed
- Sort protein names alphanumerically for each peptide entry. This should enhance stability across runs, and fixes a bug with picked-protein group FDR
- Fix another bug where picked-FDR approaches assume internal decoy generation

### Changed
- Modify order of operations during deisotoping. Deisotoped peaks can contribute intensity to only 1 parent peak now, rather than potentially multiple parent peaks

## [0.11.0]
### Added
- Support for percolator output files (`--write-pin` CLI flag)
- Support for modifying file batch size (`--batch-size N` CLI flag)
- Add `delta_best` feature, which reports the delta hyperscore from the best match to current ranked PSM
- Add Sage version to `results.json` files

### Changed
- Breaking changes to `quant` section of the configuration file format
- Rename `delta_hyperscore` to `delta_next`
- Altered internal scoring algorithm. Rather than consider all MS2 peaks within a m/z tolerance window to be matches to a theoretical spectrum, consider only the closest peak. This should increase the accuracy of # of matched peaks, and subsequent scores
- Overhaul of chimeric scoring, `report_psms` can now be used to search for multiple chimeric spectra
- Completely overhauled the LFQ algorithm: added match-between runs, peak scoring using normalized spectral angle relative to theoretical isotopic envelope, target decoy scoring of MS1 integration
- Fixed bug in picked-peptide FDR that could lead to liberal FDR
- Fixed bug in picked-protein FDR that could lead to conservative FDR
- Fixed bug where using variable protein terminal (e.g. protein N-terminal acetylation) modifications could cause some determinism. This also improves the accuracy of peptide => protein assignment. Unfortunately this fix has performance implications, causing creation of the fragment index to take up to ~2x as long.

### Removed
- Remove `no-parallel` CLI flag, and `parallel` configuration file entry

## [0.10.0]
### Added
- Retention times are now globally aligned across files
- RT prediction is then performed on all files at once (on aligned RTs), rather than one file at a time - previously, there were many instances where some files in a search could not have RTs predicted, decreasing the effectiveness of delta_rt as a feature for LDA.

### Changed
- Peptide sequences within a protein are now deduplicated - previously, repeated peptides would be called multiple times for the same protein (e.g. num_proteins > 1 even if the peptide was unique)

## [0.9.4]
### Changed
- Fix issues with RT prediction (and occasionally LDA) that arise from 0's being present on the diagonals of the covariance matrix (small amount of regularization added)

## [0.9.3]
### Added
- Allow users to set minimum number of matched b+y ions for reporting PSMs (`min_matched_peaks`)

### Changed
- Internal code for calculating factorials

## [0.9.2]
### Added
- Added option for TMT signal/noise quantification, if noise values are present in mzML

## [0.9.1]
### Changed
- FASTA file path, JSON configuration file can now be specified as "s3://" paths, allowing Sage to run completely disk-free

## [0.9.0]
### Added
- Support for non-specific digests, N-terminal enzymatic digestion

## [0.8.1]
### Added
- `quant.tmt_level` configuration option to enable MS2 (or MSn) isobaric quantification

## [0.8.0]
### Added
- Support for protein N-terminal ('['), C-terminal (']') as well as peptide C-terminal ('$') modifications
- Support for k-combinations of variable modifications. This can be specified with the `database.max_variable_mods` parameter

## [0.7.1] - 2022-11-04
### Changed
- Fix bug with in silico digest: Logic around overwriting decoys with target sequences was incorrect peptides shared between targets/decoys were being annotated as decoy peptides but assigned to non-decoy proteins. We now make sure that they are assigned to non-decoy proteins and also annotated as target sequences.

## [0.7.0] - 2022-11-03
### Added
- Add support for user-specified enzymes to JSON file.  `database.enzyme.sites` and `database.enzyme.restrict` are limited to valid amino acids
- Sage can now search MS2 spectra without annotated precursor charge states. Default behavior is to search with z=2, z=3, z=4, and then merge the PSMs for scoring

### Changed
- Configuration file schema changed. `peptide_min_len`, `peptide_max_len`, `missed_cleavages` are now specified under `database.enzyme` in the JSON file
- Internal behavior of Sage was changed to enable deterministic searching
- Docker file changed from Alpine to Debian


## [0.6.0] - 2022-11-01
### Added
- Changelog
- `rank` column added to output file
- `database.generate_decoys` parameter, which turns off internal decoy generation. This enables the use of FASTA databases for SearchGUI/PeptideShaker

### Changed
- Base ProForma v2 notation is used for peptide modifications, i.e. "\[+304.2071\]-PEPTIDEM\[+15.9949\]AAC\[+57.0214\]H"
- `scannr` column now contains the full nativeID/spectrum title from the mzML file, i.e. "controllerType=0 controllerNumber=1 scan=30069"
- `discriminant_score` column renamed to `sage_discriminant_score` for PeptideShaker recognition
- `database.decoy_prefix` JSON option changed to `database.decoy_tag`. This allows decoy tagging to occur anywhere within the accession: "sp|P01234_REVERSED|HUMAN"
- Output file renamed:  `results.pin` to `results.sage.tsv`
- Output file renamed: `quant.csv` to `quant.tsv`
- Rename `pin_paths` to `output_paths` in results.json file


## [0.5.1] - 2022-10-31
### Added
- Support for selenocysteine and pyrrolysine amino acids

## [0.5.0] - 2022-10-28
### Added
- Ability to directly read/write files from AWS S3

### Changed
- Processing files in parallel processes them in batches of `num_cpus / 2` to avoid memory issues
- Fixed bug where `protein_fdr` was erroneously assigned to `peptide_fdr` output field
- Additional parallelization for assignment of PEP, FDR, writing output files

## [0.4.0] - 2022-10-18
### Added
- Label free quantification can be enabled by turning on `quant.lfq` JSON parameter 
- Commmand line arguments can be used to override configuration file

## [0.3.1] - 2022-10-06
### Added
- Workflow contributions from [@wfondrie](https://github.com/wfondrie).

### Changed
- Don't parse empty MS2 spectra

## [0.3.0] - 2015-09-15
### Added
- Retention time prediction
- Ability to filter low-number b/y-ions for faster preliminary scoring (`database.min_ion_index` option)
- Ability to toggle retention time prediction (`predict_rt`)
