<img src="assets/xnova-logo.png" alt="XNova" width="300">

# Sage: proteomics searching so fast it seems like magic

[![Rust](https://github.com/lazear/sage/actions/workflows/rust.yml/badge.svg)](https://github.com/lazear/sage/actions/workflows/rust.yml) [![Anaconda-Server Badge](https://anaconda.org/bioconda/sage-proteomics/badges/version.svg)](https://anaconda.org/bioconda/sage-proteomics)


For more information please read [the online documentation!](https://sage-docs.vercel.app/docs)


# Introduction
 
Sage is, at it's core, a proteomics database search engine - 
    a tool that transforms raw mass spectra from proteomics experiments into peptide identifications 
    via database searching & spectral matching. 

However, Sage includes a variety of advanced features that make it a one-stop shop: retention time prediction, quantification (both isobaric & LFQ), peptide-spectrum match rescoring, and FDR control. You can directly use results from Sage without needing to use other tools for these tasks.

Additionally, Sage was designed with cloud computing in mind - massively parallel processing and the ability to directly stream compressed mass spectrometry data to/from AWS S3 enables unprecedented search speeds with minimal cost. 

 Sage also runs just as well reading local files from your Mac/PC/Linux device!

## Why use Sage instead of other tools?

Sage is **simple to configure**, **powerful** and **flexible**. 
It also happens to be well-tested, **mind-boggingly fast**, open-source (MIT-licensed) and free.

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


## Features

- Incredible performance out of the box
- [Effortlessly cross-platform](https://sage-docs.vercel.app/docs/started#download-the-latest-binary-release) (Linux/MacOS/Windows), effortlessly parallel (uses all of your CPU cores)
- [Fragment indexing strategy](https://sage-docs.vercel.app/docs/how_it_works) allows for blazing fast narrow and open searches (> 500 Da precursor tolerance)
- [Isobaric quantification](https://sage-docs.vercel.app/docs/how_it_works#tmt-based) (MS2/MS3-TMT, or custom reporter ions)
- [Label-free quantification](https://sage-docs.vercel.app/docs/how_it_works#label-free): consider all charge states & isotopologues *a la* FlashLFQ
- SILAC (duplex or triplex) in linear and crosslink searches, with MS1 ratios of every channel to light: see [SILAC.md](SILAC.md)
- Uniformly ¹⁵N / ¹³C-labelled proteins searched labelled (and optionally unlabelled) next to the rest of the FASTA: see [SILAC.md](SILAC.md#isotope-labelled-proteins-15n-13c)
- xiVIEW on this computer: `sage-gui` opens the crosslinks of a search in [xiVIEW](https://github.com/Rappsilber-Laboratory/xiview) (network, matrix, circular and 3D views, FDR, annotated spectra) in your web browser, served locally: see [CROSSLINKING.md](CROSSLINKING.md#xiview-on-this-computer)
- Two-database crosslink searches (e.g. between two complexes): crosslinks only between a Database 1 and a Database 2 protein, linear peptides in both: see [CROSSLINKING.md](CROSSLINKING.md#two-databases-crosslinks-only-between-them)
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

## Interoperability

Sage is well-integrated into the open-source proteomics ecosystem. The following projects support analyzing results from Sage (typically in addition to other tools), or redistribute Sage binaries for use in their pipelines. 

- [SearchGUI](http://compomics.github.io/projects/searchgui): a graphical user interface for running searches
- [PeptideShaker](http://compomics.github.io/projects/peptide-shaker): visualize peptide-spectrum matches
- [MS2Rescore](http://compomics.github.io/projects/ms2rescore): AI-assisted rescoring of results
- [Picked group FDR](https://github.com/kusterlab/picked_group_fdr): scalable protein group FDR for large-scale experiments
- [sagepy](https://github.com/theGreatHerrLebert/sagepy): Python bindings to the sage-core library
- [quantms](https://github.com/bigbio/quantms): nextflow pipeline for running searches with Sage
- [OpenMS](https://github.com/OpenMS/OpenMS): Sage is included as a "TOPP" tool in OpenMS
- [sager](https://github.com/UCLouvain-CBIO/sager): R package for analyzing results from Sage searches
- [Sage results to mzIdentML](https://github.com/magnuspalmblad/shic/blob/main/shims/Peptide_identification_in_TSV_to_Peptide_identification_in_mzIdentML.sh): Bash script to convert `results.sage.tsv` files to mzIdentML
- [i2MassChroQ](http://pappso.inrae.fr/bioinfo/i2masschroq/): a graphical user interface for proteomics analysis
- [annotator](https://github.com/snijderlab/annotator): a graphical user interface for visualizing peptide-spectrum matches
- [rustyms](https://github.com/snijderlab/rustyms): a Rust library (with Python bindings) to handle peptides and identified peptide files
- If your project supports Sage and it's not listed, please open a pull request! If you need help integrating or interfacing with Sage in some way, please reach out.

Check out the (now outdated) [blog post introducing the first version of Sage](https://lazear.github.io/sage/) for more information and full benchmarks!
