<img src="assets/xnova-logo.png" alt="XNova" width="300">

# XNova: Ultra-Fast Crosslinking Search Engine, Powered by Sage

[![Rust](https://github.com/lazear/sage/actions/workflows/rust.yml/badge.svg)](https://github.com/lazear/sage/actions/workflows/rust.yml)

XNova is a search engine for crosslinking mass spectrometry (XL-MS), built on [Sage](https://github.com/lazear/sage) by Michael Lazear. It keeps everything Sage does for ordinary (linear) peptides and adds crosslink searches for non-cleavable, MS-cleavable and zero-length reagents, from purified complexes to whole proteomes. It runs as the `sage` command-line program and as `sage-gui`, a desktop app for setting up searches, following them and exploring the results.

Documentation: [CROSSLINKING.md](CROSSLINKING.md) (crosslink searches, backends, outputs and viewers), [SILAC.md](SILAC.md) (SILAC and isotope-labelled proteins), the [changelog](CHANGELOG.md), and [Sage's online documentation](https://sage-docs.vercel.app/docs) for the core search engine and its configuration file.

## Introduction

### Searching whole proteomes with a non-cleavable crosslinker

A non-cleavable crosslinker such as DSS, BS3, DSG or BS2G stays intact in the mass spectrometer. The MS2 spectrum of a crosslink therefore holds the fragments of two peptides at once, and nothing in it says which two. Three difficulties follow, and all of them grow with the size of the database:

* **The search space is quadratic.** A linear search compares each spectrum with the peptides of matching mass. A crosslink search must consider every *pair* of peptides whose masses, plus the linker, add up to the precursor mass, and every pair of reactive residues in them. A human proteome digested with trypsin and two missed cleavages gives millions of peptides, and so trillions of possible pairs. Doubling the database doubles a linear search but quadruples the pairs.
* **The linked fragment ions move.** Every fragment that contains the linked residue carries the whole partner peptide plus the linker, so its mass depends on a partner that is not yet known. These ions cannot be precomputed in a fragment index the way a linear peptide's ions are. Only the fragments that do not contain the link site look like those of an ordinary peptide.
* **Chance matches grow with the search space.** The more pairs are tried, the more random pairs score well. A well-matched peptide can carry a barely matched partner, and a peptide with an unrelated modification of about the linker's mass can pass as a loop-link. Inter-protein crosslinks, whose candidate space grows fastest, are affected most, so careful scoring and a separate FDR for each kind of link matter as much as speed.

MS-cleavable linkers (DSSO, DSBU) sidestep much of this: they break in the collision cell, and each peptide's mass shows up as a stub doublet. With a non-cleavable reagent, the pairing has to be solved by the search engine.

### How Sage and Rust make it fast

Sage is written in [Rust](https://www.rust-lang.org/), a compiled language that runs as native machine code, with no garbage collector or interpreter in the way. Rust's ownership rules let Sage search spectra on all CPU cores in parallel (with [rayon](https://github.com/rayon-rs/rayon)) without data races.

Sage's central data structure is a fragment-ion index: the theoretical fragments of every peptide, sorted by m/z in buckets, each bucket sorted by precursor mass. Every peak of a spectrum is looked up in the index by binary search, so a spectrum is compared with every candidate peptide at once, and the cost follows the number of peaks rather than the number of peptides. The same index makes very wide precursor windows (open searches over more than 500 Da) cheap. That is the property XNova builds on to search crosslinks: see [the Sage algorithm, adapted to crosslinks](#the-sage-algorithm-adapted-to-crosslinks).

## Major features

### The Sage algorithm, adapted to crosslinks

The non-cleavable search reuses Sage's index instead of enumerating peptide pairs:

1. **Seed search.** The fragments of each chain that do not contain the link site are ordinary, unshifted b/y ions. XNova searches the spectrum against Sage's fragment index with a very wide precursor window, so the partner plus linker behaves like a large unknown modification, and keeps the best-matching peptides as seeds (30 per spectrum by default).
2. **Partner lookup.** For a seed `A`, the precursor fixes the partner's mass: `m(B) = M − m(A) − linker`. Peptides are sorted by mass, so the partners are found by binary search.
3. **Crosslink scoring.** Each pair, and each pair of reactive residues, is scored with a crosslink fragment model (fragments that contain the link site carry the partner and the linker) as a joint hyperscore. Each chain needs its own matched ions (`min_chain_peaks`) and its own linker-carrying ions (`min_linked_peaks`), which rejects strong-peptide-plus-weak-partner matches and look-alike loop-links.
4. **FDR.** Target–decoy, `(TD − DD) / TT` for crosslinks, computed separately for loop-links, intra-protein and inter-protein links, per spectrum (`spectrum_q`) and per unique crosslink (`crosslink_q`).

A pair search thus costs one fragment-index search and a few binary searches per seed. Every crosslink search also runs Sage's linear search on the same spectra first (on by default), so one run gives the crosslinks and the usual peptide (`results.sage.tsv`, `peptides.sage.tsv`) and protein (`proteins.sage.tsv`) results. Several enzymes can be combined for co-digestion (e.g. trypsin + Asp-N).

### Thermo and Bruker data, with ion mobility

* **Thermo:** mzML (also gzipped) and MGF files are read directly. In `sage-gui`, Thermo `.raw` files can be added as they are: they are converted to mzML with ProteoWizard msconvert or ThermoRawFileParser before the search. Isobaric reporter ions are read from MS2 or SPS-MS3 scans, optionally as signal-to-noise for Thermo files with noise data.
* **Bruker timsTOF:** `.d` folders are read directly, each precursor's spectrum merged over all its frames. The precursor's ion mobility (1/K0) is reported for every crosslink match (`ion_mobility`), used to trace crosslinks and SILAC channels in MS1 (`mobility_pct_tolerance`), and shown in the spectrum viewer and its figures.
* **Mobiligrams:** for timsTOF runs, `sage-gui` draws every searched precursor by m/z and 1/K0, shaped by charge and coloured by identification (inter-protein, intra-protein, loop-link, multivalent, linear, unidentified). A search run from `sage-gui` writes them per run and combined, as SVG and PDF.

### CPU, reflect and GPU backends

Non-cleavable searches can run on three interchangeable backends, chosen with `"backend"`:

| backend | what it does | use it for |
|---|---|---|
| `cpu` (default) | seed search and partner lookup (above) | the validated default; purified complexes and small databases |
| `reflect` | chain reflection on the CPU | a reference for the GPU algorithm, to check it |
| `gpu` | chain reflection, with the per-peptide scoring on a GPU | proteome-scale non-cleavable searches |
| `auto` | `gpu` if a real GPU is present, else `cpu` | convenience |

**Chain reflection** takes the partner out of the first step. For a crosslink `A–B`, every fragment of `A` that contains the link site is shifted by `M − m(A)`, which depends only on `A` and the precursor, not on `B`. Every peptide in the database can therefore be scored on its own, with its link-site ions shifted, which suits a GPU: the work per spectrum grows with the number of peptides, not the number of pairs. The best peptides (`chain_survivors`, 20,000 by default) are paired by mass, and the best pairs are rescored exactly on the CPU, with the same scoring and FDR as the `cpu` backend.

The GPU code uses [wgpu](https://wgpu.rs), so it runs on DirectX 12 (Windows), Vulkan (Linux) and Metal (macOS) with NVIDIA, AMD, Intel or Apple GPUs, with no CUDA dependency. The peptide ion table is uploaded once (about 1 GB for a human proteome with two missed cleavages) and split into chunks when it exceeds the device's limits. The GPU kernel is tested against the CPU reference on every adapter of a machine, `sage --gpu-info` lists the adapters, and `sage-gui` shows the GPU, its memory and whether the GPU backend can use it. MS-cleavable searches always use the doublet search on the CPU, which is already cheap.

### Label-free and isobaric quantification

Sage's quantification works on crosslinks too, and the linear search quantifies peptides and proteins as usual:

* **Label-free (MS1):** `crosslinks.lfq.tsv` gives each crosslinked species (peptide pair or loop-linked peptide) its intensity in every run. Species are traced in MS1 at every charge state over their first three isotopes, using the isotope envelope of both peptides plus the linker (and ion mobility on timsTOF data). Runs are aligned on shared crosslinks, peaks are integrated as in Sage's linear LFQ, decoy traces give each MS1 peak a q-value, and species identified in one run are quantified in the others (match-between-runs).
* **Isobaric (TMT):** reporter ions are read from each crosslink match's MS2 spectrum, or from the SPS-MS3 scan that names it as precursor, and added to `crosslinks.sage.tsv`. Labelling usually follows crosslinking, so a linked lysine carries no label: set TMT on lysine as a variable modification, and XNova links only the unlabelled form (it warns when a static modification blocks the reactive residues).
* **SILAC:** see [Isotope-incorporation searches](#isotope-incorporation-searches).

### Database vs. database searches

To study the interface between two complexes, or a bait against its candidate partners, give two FASTA files (`database_1`, `database_2`) instead of one. Crosslinks are searched only between a Database 1 protein and a Database 2 protein: pairs within one database are never scored, so the search is smaller and faster, and they do not dilute the FDR. Decoys come from the same database as their targets, so the target–decoy competition is restricted the same way. Linear peptides are searched in both files as one database, with LFQ, TMT, SILAC and labelled proteins as usual. It works with every backend and with MS-cleavable linkers; a multivalent species needs at least one link between the databases, and `database_a` / `database_b` columns say which database each peptide comes from.

### Isotope-incorporation searches

* **SILAC**, duplex or triplex (Lys4, Lys6, Lys8, Arg6, Arg10 or any mass), in linear and crosslink searches. Every peptide is searched light and in each labelled channel. A labelled lysine still links, and a crosslink joins peptides of one channel unless `mixed_crosslinks` is set (for lysates mixed before crosslinking). Each identified peptide and crosslinked species is traced in MS1 in all channels at once over one shared peak, giving intensities and ratios to light per run (`silac.sage.tsv`, `silac_proteins.sage.tsv`, `crosslinks.silac.tsv`).
* **Uniformly ¹⁵N-, ¹³C- or ¹³C¹⁵N-labelled proteins**, chosen from the FASTA, are searched labelled (and optionally unlabelled) next to the rest. Each residue carries the mass shift of its own carbon and nitrogen atoms, on top of its modifications. Labelled residues still crosslink, and crosslinks between labelled and unlabelled proteins are found (`isotope_label_a` / `isotope_label_b`). It combines with SILAC.

### Non-cleavable and MS-cleavable crosslinkers, with asymmetric handling

| kind | presets | how XNova searches them |
|---|---|---|
| non-cleavable | `dss`, `bs3`, `dsg`, `bs2g`, or any mass | seed search and partner lookup, or chain reflection (CPU or GPU) |
| MS-cleavable | `dsso`, `dsbu`, `clixlink`, or custom stubs | stub-doublet search |
| zero-length / heterobifunctional | `edc`, `dmtmm`, or custom `site1` / `site2` | as non-cleavable, each link joining a `site1` group to a `site2` group |

* **MS-cleavable linkers** break in the collision cell and leave a stub on each peptide. The two stubs of a reagent differ by a fixed mass, so each peptide appears as a doublet. XNova finds peak pairs at the doublet spacing, which gives both peptide masses, and looks each one up directly in a narrow mass window: no pairwise enumeration, and a cost that grows roughly with the size of the database rather than its square, which makes proteome-scale searches practical. Loop-links of cleavable linkers are found through the ring ions of the opened linker.
* **Asymmetric cleavable linkers** have two different arms that leave different stubs, so one peptide carries one arm's short stub and its partner the other arm's long stub. cliXlink is built in, any other can be described with a second arm (`"other"`), and every pair is tried both ways round, keeping the better orientation.
* **Heterobifunctional and zero-length reagents** react with different groups at each end: EDC and DMTMM join an amine (K, protein N-terminus) to a carboxyl (D, E, protein C-terminus), at −18.0106 Da. Every reported link joins a `site1` group to a `site2` group (never K–K or D–E), in every backend and for loop-links, and the exports name the donor and acceptor ends accordingly.

### Multivalent crosslinks

With `"multivalent": 2` to `4`, XNova also searches species held together by more than one linker: several links between two peptides (`a2b`, `a2a`), chains of three or more peptides (`a1b1c`, `a1b1c1d`) and branches (`a1b,b1c,b1d`), up to five peptides. Each peptide's fragments are modelled with the groups hanging off its link residues. For MS-cleavable linkers every linker breaks, leaving one stub per link residue (doublets and triplets), and the parts freed when only one linker breaks show which peptides were linked to which (`topology_peaks`). Multivalent species have their own FDR group, with a decoy model for any number of peptides. The search space grows steeply, so this is meant for a complex of a few dozen proteins rather than a proteome.

### Annotated spectra

`sage-gui` has the engine of [Annotator](https://github.com/snijderlab/annotator) built in. For any linear or crosslink match it draws the spectrum with each peak coloured and labelled by ion series; a peptide map with the two peptides joined by a bridge at their linked residues (loop-links as a bracket, multivalent species as one row per peptide), naming each peptide's protein and link positions; and, for MS-cleavable linkers, every stub ion in its stub's colour. Spectra are read exactly as searched (mzML, MGF and timsTOF `.d`, with the precursor's 1/K0) and open in a fraction of a second. **Export figure** saves publication figures as SVG or PDF (vector), or PNG and TIFF at 300 to 1200 dpi. Matches can also be written for the separate Annotator app (`write_annotator`), with ProForma 2.0 crosslinks.

### xiVIEW on your computer, and exports for ProXL, xiVIEW and PRIDE

* **xiVIEW on this computer:** **Open in xiVIEW** in `sage-gui` (or `sage-gui --xiview <results folder>`) shows a search's crosslinks in the Rappsilber Laboratory's [xiVIEW](https://github.com/Rappsilber-Laboratory/xiview) in your web browser, with nothing uploaded and no database or Python to install. `sage-gui` serves xiVIEW locally and answers its data requests from the results folder, in the form of the PRIDE crosslinking API. All of xiVIEW works (network, matrix, circular and 3D structure views, filters, and annotated spectra in xiSPEC) at XNova's crosslink FDR of 1%, 5% or 10%. Multivalent species are drawn dash-dotted and can be coloured by species, and timsTOF spectra are shown too.
* **ProXL:** `crosslinks_proXL.xml`, valid against the ProXL XML schema, with XNova's crosslink q-value, PSM q-value, hyperscore and delta-next score as filterable annotations, the linker (mass, stubs, reactive ends), modifications, protein sequences, and scan numbers as ProXL's reader counts them. Upload it to [ProXL](https://proxl-ms.org/) with the FASTA and the mzML files.
* **xiVIEW servers and PRIDE:** `crosslinks.mzid` (mzIdentML 1.3 with the PSI crosslinking extension, as PRIDE expects for crosslinking submissions) and `crosslinks_xiView.csv`.

See [CROSSLINKING.md](CROSSLINKING.md) for every option and output, and its [Limitations](CROSSLINKING.md#limitations) section for what has been validated so far.

## Citation

If you use Sage in a scientific publication, please cite the following paper:

[Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale](https://doi.org/10.1021/acs.jproteome.3c00486)

If you use the annotated spectrum viewer or its exported figures (`sage-gui`, "Annotated spectra" tab), please also cite Annotator, whose engine draws them:

Schulte, D.; Leuvenink, R. W.; Jager, S.; Heck, A. J. R.; Snijder, J. [A Universal Spectrum Annotator for Complex Peptidoforms in Mass Spectrometry-Based Proteomics](https://doi.org/10.1021/acs.analchem.5c02832). *Anal. Chem.* **2025**, 97 (42), 23120–23130. doi:10.1021/acs.analchem.5c02832

If you use xiVIEW from `sage-gui` ("Open in xiVIEW"), please cite xiVIEW, and xiSPEC (its spectrum viewer) and xiNET (its network view) when you use them, all from the [Rappsilber Laboratory](https://www.rappsilberlab.org/):

Combe, C. W.; Graham, M.; Kolbowski, L.; Fischer, L.; Rappsilber, J. [xiVIEW: Visualisation of Crosslinking Mass Spectrometry Data](https://doi.org/10.1016/j.jmb.2024.168656). *J. Mol. Biol.* **2024**, 436 (17), 168656. doi:10.1016/j.jmb.2024.168656

Kolbowski, L.; Combe, C.; Rappsilber, J. [xiSPEC: web-based visualization, analysis and sharing of proteomics data](https://doi.org/10.1093/nar/gky353). *Nucleic Acids Res.* **2018**, 46 (W1), W473–W478. doi:10.1093/nar/gky353

Combe, C. W.; Fischer, L.; Rappsilber, J. [xiNET: Cross-link Network Maps With Residue Resolution](https://doi.org/10.1074/mcp.O114.042259). *Mol. Cell. Proteomics* **2015**, 14 (4), 1137–1147. doi:10.1074/mcp.O114.042259

XNova serves xiVIEW the data model of the PRIDE Crosslinking Resource's API; for the PRIDE database:

Perez-Riverol, Y.; Bandla, C.; Kundu, D. J.; Kamatchinathan, S.; Bai, J.; Hewapathirana, S.; John, N. S.; Prakash, A.; Walzer, M.; Wang, S.; Vizcaíno, J. A. [The PRIDE database at 20 years: 2025 update](https://doi.org/10.1093/nar/gkae1011). *Nucleic Acids Res.* **2025**, 53 (D1), D543–D553. doi:10.1093/nar/gkae1011

## Acknowledgments

The spectrum annotation in `sage-gui` (ProForma 2.0 peptidoforms including crosslinks, theoretical fragmentation and peak matching) uses the [`mzcore`](https://crates.io/crates/mzcore), [`mzannotate`](https://crates.io/crates/mzannotate) and [`mzident`](https://crates.io/crates/mzident) libraries by Douwe Schulte, Ralf Gabriels, Auke Heerdink and contributors ([rusteomics/mzcore](https://github.com/rusteomics/mzcore)), the engine of [Annotator](https://github.com/snijderlab/annotator) from the Snijder lab at Utrecht University. They are used under the MIT license (they are dual-licensed MIT OR Apache-2.0). Spectra in the viewer are read with [`mzdata`](https://github.com/mobiusklein/mzdata) by Joshua Klein, and figures are rendered with [`resvg`](https://github.com/linebender/resvg).

xiVIEW in `sage-gui` ("Open in xiVIEW") is the work of the [Rappsilber Laboratory](https://www.rappsilberlab.org/) (Technische Universität Berlin and the University of Edinburgh; Colin Combe, Martin Graham, Lars Kolbowski, Lutz Fischer, Juri Rappsilber and contributors; funded by the Wellcome Trust). XNova builds it from their repositories, with only the changes listed in [xiview/source_fixes.py](xiview/source_fixes.py):
- [Rappsilber-Laboratory/xiview](https://github.com/Rappsilber-Laboratory/xiview) (GPL-3.0-or-later): the xiVIEW web app, with xiNET, xiSPEC and its other views; it grew out of [CLMS-UI](https://github.com/colin-combe/CLMS-UI) by Colin Combe and Martin Graham. Shipped as a separate program (the `xiview` folder next to sage-gui.exe), with its license and source references.
- [Rappsilber-Laboratory/xiview-server](https://github.com/Rappsilber-Laboratory/xiview-server) (Apache-2.0): the page shells, styles, GO terms and help pages.
- [PRIDE-Archive/crosslinking-api](https://github.com/PRIDE-Archive/crosslinking-api) and its [Rappsilber-Laboratory fork](https://github.com/Rappsilber-Laboratory/crosslinking-api) (Apache-2.0), by the PRIDE team at EMBL-EBI and the Rappsilber Laboratory: the data API of the [PRIDE Crosslinking Resource](https://www.ebi.ac.uk/pride/markdownpage/crosslinking), which XNova's local server answers in the same form; its data model comes from [mzidentml-reader](https://github.com/PRIDE-Archive/mzidentml-reader) ([Rappsilber-Laboratory/mzidentml-reader](https://github.com/Rappsilber-Laboratory/mzidentml-reader)).
- [xiAnnotator](https://github.com/Rappsilber-Laboratory/xiAnnotator) and [xi2annotator](https://github.com/Rappsilber-Laboratory/xi2annotator) (Rappsilber Laboratory): the annotation format xiSPEC reads, which XNova writes from Annotator's engine.


## Features from Sage

- Incredible performance out of the box
- [Effortlessly cross-platform](https://sage-docs.vercel.app/docs/started#download-the-latest-binary-release) (Linux/MacOS/Windows), effortlessly parallel (uses all of your CPU cores)
- [Fragment indexing strategy](https://sage-docs.vercel.app/docs/how_it_works) allows for blazing fast narrow and open searches (> 500 Da precursor tolerance)
- [Isobaric quantification](https://sage-docs.vercel.app/docs/how_it_works#tmt-based) (MS2/MS3-TMT, or custom reporter ions)
- [Label-free quantification](https://sage-docs.vercel.app/docs/how_it_works#label-free): consider all charge states & isotopologues *a la* FlashLFQ
- Capable of searching for [chimeric/co-fragmenting spectra](https://sage-docs.vercel.app/docs/configuration/additional)
- Wide-window (dynamic precursor tolerance) search mode - [enables WWA/PRM/DIA searches](https://sage-docs.vercel.app/docs/configuration/tolerance#wide-window-mode)
- Retention time prediction models fit to each LC/MS run
- [PSM rescoring](https://sage-docs.vercel.app/docs/how_it_works#machine-learning-for-psm-rescoring) using built-in linear discriminant analysis (LDA)
- PEP calculation using a non-parametric model (KDE)
- FDR calculation using target-decoy competition and picked-peptide & picked-protein approaches
- Percolator/Mokapot [compatible output](https://sage-docs.vercel.app/docs/configuration#env)
- Configuration by [JSON file](https://sage-docs.vercel.app/docs/configuration#file)
- Built-in support for reading gzipped-mzML files
- Support for reading/writing directly from [AWS S3](https://sage-docs.vercel.app/docs/configuration/aws), Google Cloud, or Azure.
