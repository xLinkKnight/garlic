# xiVIEW for XNova

`sage-gui` shows crosslink results in [xiVIEW](https://github.com/Rappsilber-Laboratory/xiview)
on this computer ("Open in xiVIEW"): it serves the xiVIEW web app on 127.0.0.1
and answers xiVIEW's data API from the results folder (`crates/sage-gui/src/xiview`).

The web app is not in this repository. `build_xiview.py` builds it from source
at pinned commits of xiVIEW (GPL-3.0-or-later) and xiview-server (Apache-2.0),
with Node.js >= 24:

    python xiview/build_xiview.py --node <folder with node.exe> --work <work folder> --out <release folder>/xiview

and the output folder goes next to `sage-gui.exe` as `xiview`. Its `SOURCE.md`
lists the commits and the changes: the page shells point at XNova's server,
and `source_fixes.py` changes xiVIEW's source (multivalent crosslinks in the
spectrum, network and circle views; xiSPEC's settings request id; loop-link
spectra). For development, `XNOVA_XIVIEW_DIR`
points sage-gui at another copy.

See CROSSLINKING.md, "xiVIEW on this computer".

## Credit

xiVIEW, xiNET and xiSPEC are the work of the Rappsilber Laboratory
(https://www.rappsilberlab.org/; Colin Combe, Martin Graham, Lars Kolbowski,
Lutz Fischer, Juri Rappsilber and contributors). Please cite:

- Combe, C. W.; Graham, M.; Kolbowski, L.; Fischer, L.; Rappsilber, J. xiVIEW:
  Visualisation of Crosslinking Mass Spectrometry Data. *J. Mol. Biol.* 2024,
  436 (17), 168656. https://doi.org/10.1016/j.jmb.2024.168656
- Kolbowski, L.; Combe, C.; Rappsilber, J. xiSPEC: web-based visualization,
  analysis and sharing of proteomics data. *Nucleic Acids Res.* 2018, 46 (W1),
  W473–W478. https://doi.org/10.1093/nar/gky353
- Combe, C. W.; Fischer, L.; Rappsilber, J. xiNET: Cross-link Network Maps With
  Residue Resolution. *Mol. Cell. Proteomics* 2015, 14 (4), 1137–1147.
  https://doi.org/10.1074/mcp.O114.042259

Repositories: https://github.com/Rappsilber-Laboratory/xiview (GPL-3.0),
https://github.com/Rappsilber-Laboratory/xiview-server (Apache-2.0); the data
API XNova answers is the PRIDE Crosslinking Resource's,
https://github.com/PRIDE-Archive/crosslinking-api (PRIDE team, EMBL-EBI, and
the Rappsilber Laboratory; data model from
https://github.com/PRIDE-Archive/mzidentml-reader); the annotation format is
xiAnnotator's, https://github.com/Rappsilber-Laboratory/xiAnnotator.
