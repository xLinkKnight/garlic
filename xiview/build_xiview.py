"""Build the xiVIEW web app that sage-gui serves locally ("Open in xiVIEW").

    python xiview/build_xiview.py --node <folder with node.exe> --work <work folder> --out <output folder>

The output folder goes next to sage-gui.exe, named `xiview`.

* xiVIEW (https://github.com/Rappsilber-Laboratory/xiview, GPL-3.0-or-later)
  is built from source at a pinned commit with Node.js >= 24
  (`npm ci --ignore-scripts`, then webpack in production mode).
* The page shells, fonts, images, GO terms and help pages come from
  xiview-server (https://github.com/Rappsilber-Laboratory/xiview-server,
  Apache-2.0), also pinned; the help videos (37 MB) are left out.
* The two page shells (network.html, spectra.html) are changed to load their
  data from garlic's local server (`./data/`, `./xiAnnotator/`) instead of
  the PRIDE web service, and the link to the xiSEARCH server is removed.
* Changes to xiVIEW's source before the build (source_fixes.py): the
  spectrum settings' request id; multivalent crosslinks in the spectrum
  viewer (any number of peptides and links), the network and circle views
  (garlic's link metadata loaded with the data, multivalent links
  dash-dotted).

xiVIEW stays a separate program: sage-gui talks to it over HTTP on
127.0.0.1, it is not compiled into sage-gui.
"""

import argparse
import importlib.util
import os
import re
import shutil
import subprocess
import sys

XIVIEW = ('https://github.com/Rappsilber-Laboratory/xiview.git', '03a49f75d1f832b1e103dd78e318fe49d5354afb')
XIVIEW_SERVER = ('https://github.com/Rappsilber-Laboratory/xiview-server.git', '9258369aa4e75f50f7ebac687dc7a4b89bb943e5')

# garlic's changes to xiVIEW's source (next to this script; loaded by path,
# as `python -I` leaves the script's folder off the module path)
_spec = importlib.util.spec_from_file_location(
    'source_fixes', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'source_fixes.py'))
source_fixes = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(source_fixes)


def run(cmd, cwd=None, env=None):
    print('>', ' '.join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def checkout(url, commit, folder):
    if not os.path.isdir(os.path.join(folder, '.git')):
        run(['git', 'clone', '-q', url, folder])
    run(['git', '-C', folder, 'fetch', '-q', 'origin', commit])
    run(['git', '-C', folder, 'checkout', '-q', commit])
    # the pinned sources as they are (the fixes are applied to them afresh)
    run(['git', '-C', folder, 'checkout', '-q', '--', '.'])


def patch_shell(path, entry):
    """Point a page shell at garlic's local server"""
    text = open(path, encoding='utf-8').read()
    api = re.compile(r"const APIS = \{.*?xiview\.%s\([^;]*\);" % entry, re.S)
    assert api.search(text), f'{path}: API block not found'
    text = api.sub(
        "// garlic serves xiVIEW's data API and its annotator on this computer\n"
        "    xiview.%s(\"./data/\", \"./xiAnnotator/\");" % entry,
        text,
    )
    # the xiSEARCH server link (home icon)
    text = re.sub(r'\s*<i class="fa fa-home fa-xi"[^>]*>\s*</i>', '', text, flags=re.S)
    open(path, 'w', encoding='utf-8', newline='\n').write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node', required=True, help='folder holding node.exe (Node.js >= 24)')
    ap.add_argument('--work', required=True, help='folder for the sources and build')
    ap.add_argument('--out', required=True, help='output folder (the `xiview` folder next to sage-gui.exe)')
    a = ap.parse_args()

    os.makedirs(a.work, exist_ok=True)
    src = os.path.join(a.work, 'xiview')
    server = os.path.join(a.work, 'xiview-server')
    checkout(*XIVIEW, src)
    checkout(*XIVIEW_SERVER, server)
    source_fixes.apply(src)

    env = dict(os.environ)
    env['PATH'] = a.node + os.pathsep + env.get('PATH', '')
    env['PUPPETEER_SKIP_DOWNLOAD'] = '1'
    node = os.path.join(a.node, 'node.exe' if os.name == 'nt' else 'node')
    npm = os.path.join(a.node, 'npm.cmd' if os.name == 'nt' else 'npm')
    # --ignore-scripts: no install hooks (xiVIEW's `prepare` installs git hooks)
    run([npm, 'ci', '--ignore-scripts', '--no-audit', '--no-fund'], cwd=src, env=env)
    run([node, os.path.join('node_modules', 'webpack', 'bin', 'webpack.js'), '--mode', 'production',
         '--config', 'webpack.prod.js'], cwd=src, env=env)

    if os.path.isdir(a.out):
        shutil.rmtree(a.out)
    static = os.path.join(server, 'static')
    shutil.copytree(static, a.out, ignore=shutil.ignore_patterns('vid', 'xiview.js', 'xiview.js.LICENSE.txt',
                                                                  'vendors.js', 'vendors.js.LICENSE.txt'))
    dist = os.path.join(src, 'dist')
    for name in os.listdir(dist):
        if not name.startswith('tests_'):
            shutil.copy2(os.path.join(dist, name), os.path.join(a.out, name))
    patch_shell(os.path.join(a.out, 'network.html'), 'main')
    patch_shell(os.path.join(a.out, 'spectra.html'), 'validationPage')

    shutil.copy2(os.path.join(src, 'LICENSE'), os.path.join(a.out, 'LICENSE-xiview.txt'))
    if os.path.isfile(os.path.join(src, 'NOTICE')):
        shutil.copy2(os.path.join(src, 'NOTICE'), os.path.join(a.out, 'NOTICE-xiview.txt'))
    shutil.copy2(os.path.join(server, 'LICENSE.md'), os.path.join(a.out, 'LICENSE-xiview-server.md'))
    with open(os.path.join(a.out, 'SOURCE.md'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(
            '# xiVIEW in garlic\n\n'
            'This folder is the xiVIEW web app, served by sage-gui on 127.0.0.1 ("Open in xiVIEW").\n\n'
            f'* xiVIEW: {XIVIEW[0]} at commit {XIVIEW[1]} (GPL-3.0-or-later, LICENSE-xiview.txt),\n'
            '  built with `npm ci --ignore-scripts` and webpack (production): xiview.js, vendors.js and fonts.\n'
            f'* xiview-server: {XIVIEW_SERVER[0]} at commit {XIVIEW_SERVER[1]} (Apache-2.0,\n'
            '  LICENSE-xiview-server.md): network.html, spectra.html, styles, images, GO terms, help pages\n'
            '  (without the help videos).\n'
            '* Changed: network.html and spectra.html load their data from garlic\'s local server\n'
            '  (`./data/`, `./xiAnnotator/`) instead of the PRIDE web service; the xiSEARCH server link\n'
            '  is removed. xiVIEW\'s source is changed before the build as listed in garlic\'s\n'
            '  xiview/source_fixes.py: the spectrum settings\' request id (was NaN, so a re-annotated\n'
            '  spectrum was dropped); multivalent crosslinks in the spectrum viewer (any number of\n'
            '  peptides and links, one linker mass per link) and in the network and circle views\n'
            '  (garlic\'s link metadata loaded with the data; multivalent links dash-dotted).\n\n'
            '## Credit and citation\n\n'
            'xiVIEW, xiNET and xiSPEC are the work of the Rappsilber Laboratory (https://www.rappsilberlab.org/;\n'
            'Colin Combe, Martin Graham, Lars Kolbowski, Lutz Fischer, Juri Rappsilber and contributors;\n'
            'funded by the Wellcome Trust). Please cite:\n\n'
            '* Combe, C. W.; Graham, M.; Kolbowski, L.; Fischer, L.; Rappsilber, J. xiVIEW: Visualisation of\n'
            '  Crosslinking Mass Spectrometry Data. J. Mol. Biol. 2024, 436 (17), 168656.\n'
            '  https://doi.org/10.1016/j.jmb.2024.168656\n'
            '* Kolbowski, L.; Combe, C.; Rappsilber, J. xiSPEC: web-based visualization, analysis and sharing\n'
            '  of proteomics data. Nucleic Acids Res. 2018, 46 (W1), W473-W478. https://doi.org/10.1093/nar/gky353\n'
            '* Combe, C. W.; Fischer, L.; Rappsilber, J. xiNET: Cross-link Network Maps With Residue\n'
            '  Resolution. Mol. Cell. Proteomics 2015, 14 (4), 1137-1147. https://doi.org/10.1074/mcp.O114.042259\n\n'
            'The data API garlic answers is the PRIDE Crosslinking Resource\'s\n'
            '(https://github.com/PRIDE-Archive/crosslinking-api; PRIDE team, EMBL-EBI, and the Rappsilber\n'
            'Laboratory; data model from https://github.com/PRIDE-Archive/mzidentml-reader); the annotation\n'
            'format xiSPEC reads is xiAnnotator\'s (https://github.com/Rappsilber-Laboratory/xiAnnotator).\n'
            'PRIDE: Perez-Riverol, Y. et al. The PRIDE database at 20 years: 2025 update. Nucleic Acids Res.\n'
            '2025, 53 (D1), D543-D553. https://doi.org/10.1093/nar/gkae1011\n\n'
            'Rebuild with garlic\'s xiview/build_xiview.py.\n'
        )
    print('xiVIEW ready in', a.out)


if __name__ == '__main__':
    sys.exit(main())
