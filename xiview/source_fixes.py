"""Changes XNova makes to xiVIEW's source before building it.

Each is (file, text, replacement); the text must occur exactly once, so a
change in xiVIEW's source stops the build instead of being missed.

1. xiSPEC's settings: the request id was NaN when the spectrum was opened
   from the match table, so the re-annotated spectrum was dropped. A
   loop-link's spectrum did not open at all (no second peptide).
2. Multivalent crosslinks, spectrum viewer (xiSPEC): any number of peptides
   (a row each, peptides 3 and 4 in their own colours) and of links (a line
   each; a loop-link as a bracket), and the precursor mass with one linker
   per link. Two peptides with one link are drawn exactly as before.
3. Multivalent crosslinks, network (xiNET) and circle views: XNova's link
   metadata (`get_xnova_link_metadata`: `Multivalent`, `Multivalent
   species`) is loaded with the data, and links of multivalent species are
   drawn dash-dotted. "Multivalent species" is also a link colour scheme
   (one colour per species).
"""

import os

J = os.path.join
XISPEC = J('js', 'views', 'xispec')

SOURCE_FIXES = [
    # --- 1. settings request id
    (
        J(XISPEC, 'DataSettingsView.js'),
        'json_request["annotation"]["requestID"] = xiSPECUI.lastRequestedID + Date.now();',
        'json_request["annotation"]["requestID"] = String(xiSPECUI.lastRequestedID || "settings") + "_" + Date.now();',
    ),

    # --- 1b. loop-links: their spectrum did not open (one peptide: no stubs2)
    (
        J('js', 'load-spectrum.js'),
        '    formatted_data.stubs2 = match.matchedPeptides[1].stubs;',
        '    formatted_data.stubs2 = match.matchedPeptides[1] ? match.matchedPeptides[1].stubs : undefined;',
    ),

    # --- 1c. a spectrum that cannot be fetched or annotated: say so (it was
    # silent; e.g. sage-gui, which serves the data, was closed)
    (
        J('js', 'load-spectrum.js'),
        '/**\n * Load spectrum data from PRIDE for a given spectrum match',
        'import {showMessage} from "./ui-utils/message-popup";\n\n/**\n * Load spectrum data from PRIDE for a given spectrum match',
    ),
    (
        J('js', 'load-spectrum.js'),
        '''        .then(response => {
            if (!response.ok) {
                throw new Error("Network response was not ok");
            }
            return response.json();
        })''',
        '''        .then(response => {
            if (!response.ok) {
                // XNova: the server's reason
                return response.text().then(text => {
                    let detail = text;
                    try {
                        const json = JSON.parse(text);
                        detail = json.error || json.detail || text;
                    } catch (e) { /* plain text */ }
                    throw new Error("XNova could not show this spectrum: " + detail);
                });
            }
            return response.json();
        })''',
    ),
    (
        J('js', 'load-spectrum.js'),
        '''        .catch(error => {
            console.log("error getting annotated peak list", error);
        });''',
        '''        .catch(error => {
            console.log("error getting annotated peak list", error);
            showMessage(error instanceof TypeError ? XNOVA_NOT_RUNNING : String(error.message || error), "error");
        });''',
    ),
    (
        J('js', 'load-spectrum.js'),
        'export const loadSpectrum = function (match) {',
        '''// XNova: a request that reached no server
const XNOVA_NOT_RUNNING = "The spectrum could not be fetched: sage-gui serves this page's data and spectra, "
    + "so keep sage-gui open while you use xiVIEW (if it was closed, open the results again with Open in xiVIEW).";

export const loadSpectrum = function (match) {''',
    ),
    (
        J(XISPEC, 'spectrum-wrapper.js'),
        'import d3 from "d3";',
        'import d3 from "d3";\nimport {showMessage} from "../../ui-utils/message-popup";',
    ),
    (
        J(XISPEC, 'spectrum-wrapper.js'),
        '''            .catch(error => {
                console.error("Error during annotation request:", error);
            });''',
        '''            .catch(error => {
                console.error("Error during annotation request:", error);
                // XNova: say so (it was silent)
                showMessage(error instanceof TypeError
                    ? "The spectrum could not be re-annotated: keep sage-gui open while you use xiVIEW."
                    : "The spectrum could not be re-annotated: " + (error.message || error), "error");
            });''',
    ),

    # --- 2. xiSPEC: peptide colours for any number of peptides
    (
        J(XISPEC, 'annotated-spectrum-model.js'),
        '''    updateColors: function () {
        switch (this.get("visFragments")) {''',
        '''    // XNova: the colour of peptide p (multivalent crosslinks have more than
    // two); kind "", "_loss" or "_cluster"
    pepColor: function (p, kind) {
        if (p === 0) return this["p1color" + kind];
        if (p === 1) return this["p2color" + kind];
        // gold, grey, cyan: apart from every xiSPEC colour scheme's two colours
        const extra = [["#c49a00", "#e0c040", "#f0dc90"], ["#4d4d4d", "#8c8c8c", "#cccccc"], ["#17becf", "#7fdbe6", "#c4f0f5"]];
        const c = extra[(p - 2) % extra.length];
        return kind === "_loss" ? c[1] : (kind === "_cluster" ? c[2] : c[0]);
    },

    updateColors: function () {
        switch (this.get("visFragments")) {''',
    ),
    (
        J(XISPEC, 'annotated-spectrum-model.js'),
        '''        // NOT Multilink future proof
        if (JSONdata.LinkSite.length > 1) {
            if (JSONdata.LinkSite[0].linkSite !== -1 && JSONdata.LinkSite[1].linkSite !== -1)
                totalMass += clModMass;
        }''',
        '''        // XNova: one linker per link (both ends set; multivalent crosslinks have several)
        const linkEnds = {};
        JSONdata.LinkSite.forEach(function (l) {
            if (l.linkSite !== -1) linkEnds[l.id] = (linkEnds[l.id] || 0) + 1;
        });
        const linkCount = Object.keys(linkEnds).filter(function (id) { return linkEnds[id] >= 2; }).length;
        totalMass += clModMass * linkCount;''',
    ),
    (
        J(XISPEC, 'graph', 'Peak.js'),
        'this.colour = (non_lossy[0].peptideId === 0) ? model.p1color : model.p2color;',
        'this.colour = model.pepColor(non_lossy[0].peptideId, "");',
    ),
    (
        J(XISPEC, 'graph', 'Peak.js'),
        'this.colour = (this.fragments[0].peptideId === 0) ? model.p1color_loss : model.p2color_loss;',
        'this.colour = model.pepColor(this.fragments[0].peptideId, "_loss");',
    ),
    (
        J(XISPEC, 'graph', 'Peak.js'),
        '''        if (this.isotopes[0].peptideId === 0)
            this.colour = model.p1color_cluster;
        if (this.isotopes[0].peptideId === 1)
            this.colour = model.p2color_cluster;''',
        '''        this.colour = model.pepColor(this.isotopes[0].peptideId, "_cluster");''',
    ),
    (
        J(XISPEC, 'ErrorPlotView.js'),
        '''            if (d.lossy) {
                if (d["peptideId"] == 0) return self.model.p1color_loss;
                else return self.model.p2color_loss;
            } else {
                if (d["peptideId"] == 0) return self.model.p1color;
                else return self.model.p2color;
            }''',
        '''            return self.model.pepColor(+d["peptideId"], d.lossy ? "_loss" : "");''',
    ),
    (
        J(XISPEC, 'FragKey', 'KeyFragment.js'),
        '''        this.peptideId = peptideId;
        if (this.peptideId === 0) {
            this.color = this.FragKey.model.p1color;
        } else if (this.peptideId === 1) {
            this.color = this.FragKey.model.p2color;
        }''',
        '''        this.peptideId = peptideId;
        this.color = this.FragKey.model.pepColor(this.peptideId, "");''',
    ),
    (
        J(XISPEC, 'FragKey', 'KeyFragment.js'),
        '''        if (this.peptideId === 0)
            y_coord = 25;
        if (this.peptideId === 1)
            y_coord = 75;''',
        '''        // XNova: a row per peptide
        y_coord = 25 + 50 * this.peptideId;''',
    ),

    # --- 2. xiSPEC: the fragmentation key with a row per peptide and a line per link
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '        this.pepoffset = [0, 0];\n        for (let p = 0; p < pepCount; p++) {',
        '        this.pepoffset = this.model.peptides.map(function () { return 0; });\n        for (let p = 0; p < pepCount; p++) {',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''        if (this.linkPos.length > 0) {
            this.CLpos = this.linkPos[0].linkSite;
            if (this.linkPos[1].linkSite > this.linkPos[0].linkSite)
                this.CLpos = this.linkPos[1].linkSite;
        } else {
            this.CLpos = -1;
        }

        this.tooltip.style("opacity", 0);

        this.align_peptides_to_CL();''',
        '''        this.tooltip.style("opacity", 0);

        // XNova: any number of peptides and links (multivalent crosslinks)
        this.links = this.linksOf(this.linkPos);
        this.align_peptides_to_CL();
        this.CLpos = this.links.length > 0 ? this.pepoffset[this.links[0].a] + this.links[0].sa : -1;''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''            this.drawFragmentationEvents(0);
            if (this.peptides[1])
                this.drawFragmentationEvents(1);
            this.drawCleavedCrosslinkerFragmentationEvents();''',
        '''            for (let p = 0; p < this.peptides.length; p++)
                this.drawFragmentationEvents(p);
            this.drawCleavedCrosslinkerFragmentationEvents();''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''        // CL line svg elements
        if (this.linkPos.length > 0) {''',
        '''        // CL line svg elements (XNova: the first link between two peptides)
        if (this.links.length > 0 && this.links[0].a !== this.links[0].b) {''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''        // change-mod svg element
        let changeModLetterG = this.scaleSvgGroup.append("g");''',
        '''        // XNova: every link where it is (multivalent crosslinks, loop-links)
        this.drawLinks();

        // change-mod svg element
        let changeModLetterG = this.scaleSvgGroup.append("g");''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''    align_peptides_to_CL: function () {
        if (this.linkPos.length > 0) {

            // #==========================================================================
            // #	account for crosslink shift
            // #	this aligns the peptide sequences at the crosslink site
            // #==========================================================================
            let shift = this.linkPos[0].linkSite - this.linkPos[1].linkSite;
            if (shift <= 0) {
                this.peptides[0] = Array(Math.abs(shift) + 1).join("#") + this.peptideStrs[0];
                this.pepoffset[0] = Math.abs(shift);
            } else {
                this.peptides[1] = Array(shift + 1).join("#") + this.peptideStrs[1];
                this.pepoffset[1] = shift;
            }

            let diff = this.peptideStrs[0].length - this.peptideStrs[1].length;
            if (diff <= 0) {
                this.peptides[0] = this.peptides[0] + Array(Math.abs(diff) + 1).join("#");
            } else {
                this.peptides[1] = this.peptides[1] + Array(diff + 1).join("#");
            }
        }
    },''',
        '''    // XNova: the links as {a, sa, b, sb} (peptides and residues), the two ends
    // of each link sharing an id
    linksOf: function (linkPos) {
        const byId = new Map();
        (linkPos || []).forEach(function (l) {
            if (l.linkSite === undefined || l.linkSite < 0) return;
            if (!byId.has(l.id)) byId.set(l.id, []);
            byId.get(l.id).push(l);
        });
        const links = [];
        byId.forEach(function (ends) {
            if (ends.length >= 2)
                links.push({a: ends[0].peptideId, sa: ends[0].linkSite, b: ends[1].peptideId, sb: ends[1].linkSite});
        });
        return links;
    },

    // XNova: line the peptides up at their linked residues (the first link
    // first, then the peptides the links reach), any number of peptides
    align_peptides_to_CL: function () {
        const n = this.peptideStrs.length;
        const offset = this.peptideStrs.map(function () { return 0; });
        const placed = this.peptideStrs.map(function () { return false; });
        if (this.links.length > 0) {
            placed[this.links[0].a] = true;
            let grew = true;
            while (grew) {
                grew = false;
                for (const l of this.links) {
                    if (placed[l.a] && !placed[l.b]) {
                        offset[l.b] = offset[l.a] + l.sa - l.sb;
                        placed[l.b] = true;
                        grew = true;
                    } else if (placed[l.b] && !placed[l.a]) {
                        offset[l.a] = offset[l.b] + l.sb - l.sa;
                        placed[l.a] = true;
                        grew = true;
                    }
                }
            }
            const min = Math.min.apply(null, offset);
            for (let p = 0; p < n; p++) offset[p] -= min;
        }
        let width = 0;
        for (let p = 0; p < n; p++) width = Math.max(width, offset[p] + this.peptideStrs[p].length);
        for (let p = 0; p < n; p++) {
            this.pepoffset[p] = offset[p];
            if (this.links.length > 0)
                this.peptides[p] = Array(offset[p] + 1).join("#") + this.peptideStrs[p]
                    + Array(width - offset[p] - this.peptideStrs[p].length + 1).join("#");
        }
    },

    // XNova: a line per link between the rows of its peptides (the first
    // one is xiSPEC's own, clickable line); a loop-link as a bracket
    drawLinks: function () {
        const self = this;
        const x = function (p, s) { return self.xStep * (self.pepoffset[p] + s); };
        this.links.forEach(function (l, i) {
            if (l.a === l.b) {
                const y = 25 + 51 * l.a;
                self.scaleSvgGroup.append("path")
                    .attr("d", "M" + x(l.a, l.sa) + "," + y + " L" + x(l.a, l.sa) + "," + (y + 8)
                        + " L" + x(l.b, l.sb) + "," + (y + 8) + " L" + x(l.b, l.sb) + "," + y)
                    .attr("fill", "none")
                    .attr("stroke", "black")
                    .attr("stroke-width", 2.3);
                return;
            }
            const up = l.a < l.b ? l : {a: l.b, sa: l.sb, b: l.a, sb: l.sa};
            const coords = {x1: x(up.a, up.sa), y1: 25 + 51 * up.a, x2: x(up.b, up.sb), y2: 55 + 51 * (up.b - 1)};
            if (i === 0 && self.CL) {
                [self.CLlineHighlight, self.CLline, self.changeCLline].forEach(function (line) {
                    line.attr(coords);
                });
            } else {
                self.scaleSvgGroup.append("line")
                    .attr(coords)
                    .attr("stroke", "black")
                    .attr("stroke-width", 2.3);
            }
        });
    },''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''        let pep1SvgGroup = self.scaleSvgGroup
            .append("g")
            .attr("class", "peptide");
        let peptides = [{
            sequence: this.peptides[0],
            color: this.model.p1color,
            y: [20, 5],
            group: pep1SvgGroup
        }];

        if (this.peptides.length > 1) {
            let pep2SvgGroup = self.scaleSvgGroup
                .append("g")
                .attr("class", "peptide");
            peptides.push({
                sequence: this.peptides[1],
                color: this.model.p2color,
                y: [71, 83],
                group: pep2SvgGroup
            });
        }''',
        '''        // XNova: a row per peptide (multivalent crosslinks have more than two)
        let peptides = this.peptides.map(function (sequence, p) {
            return {
                sequence: sequence,
                color: self.model.pepColor(p, ""),
                y: p === 0 ? [20, 5] : [20 + 51 * p, 32 + 51 * p],
                group: self.scaleSvgGroup.append("g").attr("class", "peptide"),
            };
        });''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''                if (pepLetterData.pepIndex === 0) {
                    y = 5;
                    color = self.model.p1color;
                } else if (pepLetterData.pepIndex === 1) {
                    y = 83;
                    color = self.model.p2color;
                }''',
        '''                y = pepLetterData.pepIndex === 0 ? 5 : 32 + 51 * pepLetterData.pepIndex;
                color = self.model.pepColor(pepLetterData.pepIndex, "");''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''            const arrays = [[], []];''',
        '''            const arrays = this.peptides.map(function () { return []; });''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''            if (!duplicates(arrays[0]) && !duplicates(arrays[1]))''',
        '''            if (arrays.every(function (a) { return !duplicates(a); }))''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''        if (fragments === "all") {
            color(0, this.model.p1color, 0, this.pepLetters[0][0].length);
            if (this.peptides[1])
                color(1, this.model.p2color, 0, this.pepLetters[1][0].length);
        } else {
            for (let f = 0; f < fragments.length; f++) {
                for (let i = 0; i < fragments[f].range.length; i++) {
                    if (fragments[f].range[i].peptideId === 0)
                        color(0, this.model.p1color, fragments[f].range[i].from, fragments[f].range[i].to + 1);
                    if (fragments[f].range[i].peptideId === 1)
                        color(1, this.model.p2color, fragments[f].range[i].from, fragments[f].range[i].to + 1);
                }
            }
        }''',
        '''        if (fragments === "all") {
            for (let p = 0; p < this.pepLetters.length; p++)
                color(p, this.model.pepColor(p, ""), 0, this.pepLetters[p][0].length);
        } else {
            for (let f = 0; f < fragments.length; f++) {
                for (let i = 0; i < fragments[f].range.length; i++) {
                    const p = fragments[f].range[i].peptideId;
                    if (this.pepLetters[p])
                        color(p, this.model.pepColor(p, ""), fragments[f].range[i].from, fragments[f].range[i].to + 1);
                }
            }
        }''',
    ),
    (
        J(XISPEC, 'FragmentationKeyView.js'),
        '''        let scale = 1;
        if (parentWidth < fragKeyWidth + this.margin.left + this.margin.right) {
            scale = parentWidth / (fragKeyWidth + this.margin.left + this.margin.right);
        }''',
        '''        let scale = 1;
        if (parentWidth < fragKeyWidth + this.margin.left + this.margin.right) {
            scale = parentWidth / (fragKeyWidth + this.margin.left + this.margin.right);
        }
        // XNova: more than two peptides: about the height of two
        if (this.peptides && this.peptides.length > 2 && fragKeyHeight > 0)
            scale = Math.min(scale, 95 / fragKeyHeight);''',
    ),

    # --- 3. network and circle views: XNova's link metadata, multivalent links dash-dotted
    (
        J('js', 'main.js'),
        '''import vent from "./vent";''',
        '''import vent from "./vent";
import {updateLinkMetadata} from "./modelUtils";''',
    ),
    (
        J('js', 'main.js'),
        '''    // Add the blosum promise to the tasks array
    const allTasks = [...tasks, blosumPromise];''',
        '''    // XNova: multivalent crosslinks as link metadata (CSV; empty without them)
    let xnovaLinkMetadata = "";
    const xnovaTask = fetch(`${apiBase}get_xnova_link_metadata${window.location.search}`)
        .then(response => response.ok ? response.text() : "")
        .then(text => { xnovaLinkMetadata = text; })
        .catch(() => {});

    // Add the blosum promise to the tasks array
    const allTasks = [...tasks, blosumPromise, xnovaTask];''',
    ),
    (
        J('js', 'main.js'),
        '''            const compositeModelInst = models({}, clmsModel);
            compositeModelInst.set("blosumColl", blosumCollInst);''',
        '''            const compositeModelInst = models({}, clmsModel);
            if (xnovaLinkMetadata.trim()) {
                updateLinkMetadata(xnovaLinkMetadata, clmsModel);
            }
            compositeModelInst.set("blosumColl", blosumCollInst);''',
    ),
    (
        J('js', 'views', 'xinet', 'link', 'rendered-crosslink.js'),
        '''        this.dashedLine(this.crosslink.ambiguous && this.crosslink.isMonoLink() === false);
''',
        '''        this.dashedLine(this.crosslink.ambiguous && this.crosslink.isMonoLink() === false);
        // XNova: links of multivalent crosslinks, dash-dotted
        if (this.shown && !this.crosslink.ambiguous && this.crosslink.getMeta("Multivalent")) {
            this.line.setAttribute("stroke-dasharray", "10, 3, 2, 3");
        }
''',
    ),
    (
        J('js', 'views', 'circle', 'circularViewBB.js'),
        '''            .classed("ambiguous", function (d) {
                return crosslinks.get(d.id).ambiguous;
            });''',
        '''            .classed("ambiguous", function (d) {
                return crosslinks.get(d.id).ambiguous;
            })
            // XNova: links of multivalent crosslinks, dash-dotted
            .classed("xnovaMultivalent", function (d) {
                const crosslink = crosslinks.get(d.id);
                return !crosslink.ambiguous && !!crosslink.getMeta("Multivalent");
            });''',
    ),
    (
        J('css', 'circularViewBB.css'),
        '''.circleLink {''',
        '''/* XNova: links of multivalent crosslinks */
.circleLink.xnovaMultivalent {
    stroke-dasharray: 10, 3, 2, 3;
}

.circleLink {''',
    ),
]


def apply(src):
    """Apply the fixes to an xiVIEW checkout (already applied ones are left)"""
    for name, old, new in SOURCE_FIXES:
        path = os.path.join(src, name)
        text = open(path, encoding='utf-8').read()
        if new in text:
            continue
        assert text.count(old) == 1, f'{name}: fix does not apply ({text.count(old)} matches): {old[:60]!r}'
        open(path, 'w', encoding='utf-8', newline='\n').write(text.replace(old, new))
