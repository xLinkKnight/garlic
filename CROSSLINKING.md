# Crosslinking (XL-MS) searches

Sage can search for chemically crosslinked peptides. Add a `crosslink` block to
the usual JSON configuration to enable it; omit it and Sage behaves exactly as
before. Two strategies are used, chosen by the reagent:

* **MS-cleavable** reagents (**DSSO**, **DSBU**, **cliXlink**): `"crosslink": {"reagent": "dsso"}`.
  Recommended for proteome-scale work - see [Cleavable linkers](#ms-cleavable-linkers-dsso-dsbu-clixlink).
* **Non-cleavable** reagents (DSS, BS3, DSG, BS2G, custom): the default.
* **Zero-length / heterobifunctional** reagents (**EDC**, **DMTMM**, custom
  `site1`/`site2`): ends that react with different groups - see
  [Heterobifunctional and zero-length reagents](#heterobifunctional-and-zero-length-reagents-edc-dmtmm).

```json
{
  "database": {
    "fasta": "proteins.fasta",
    "enzyme": { "missed_cleavages": 2, "min_len": 6, "max_len": 40 },
    "generate_decoys": true
  },
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol":  { "ppm": [-10, 10] },
  "precursor_charge": [3, 6],
  "deisotope": true,
  "crosslink": {}
}
```

`"crosslink": {}` means DSS/BS3 on lysine and protein N-termini. All options:

| key | default | meaning |
|---|---|---|
| `reagent` | none | preset: `dss`, `bs3`, `dsg`, `bs2g`, `dsso`, `dsbu`, `clixlink`, `edc`, `dmtmm`. Sets the linker mass (and stub masses for cleavable reagents, sites for EDC/DMTMM) |
| `mass` | `138.06808` | custom reagents only (not with `reagent`): mass added to the *pair* of peptides; may be negative down to -20 for zero-length couplings (-18.010565) |
| `cleavable` | none | MS-cleavable stub masses: `{"short": .., "long": .., "alt": [..]}`; `short + long` must equal the linker mass. An asymmetric linker adds its second arm as `"other": {"short": .., "long": .., "alt": [..]}`; then each arm's short stub + the other arm's long stub must equal the linker mass. Required for custom cleavable reagents; with a cleavable preset it replaces the preset stubs |
| `backend` | `"cpu"` | non-cleavable search algorithm / hardware: `cpu`, `reflect`, `gpu`, `auto` - see [Choosing a backend](#choosing-a-backend) |
| `chain_survivors` | `20000` | `reflect`/`gpu` only: peptides kept per spectrum after the per-chain filter (best by chain score) |
| `min_linked_peaks` | `2` | matched fragment ions that carry the linker required on *each* chain, reported as `linked_a` / `linked_b`. Non-cleavable: ions that contain the link site (loop-links: both sites). Cleavable loop-links: stub-carrying ions inside the opened ring. They are the only direct evidence of the link; without them a correct peptide with any modification of about the linker's mass passes as a loop-link, and a strong peptide with a barely matched partner as a crosslink. `0` turns the check off |
| `linear_search` | `true` | also run an ordinary linear peptide search first, written to `results.sage.tsv`, `peptides.sage.tsv` and `proteins.sage.tsv` (plus the usual LFQ/TMT tables, `.pin` and HTML report when enabled). It uses the crosslink search's database settings, e.g. `peptide_min_mass` 350 Da. `false` skips it |
| `require_complete_doublets` | `false` | cleavable only: require both stub peaks of *both* peptides (default: both of one, at least one of the other) |
| `residues` | `"K"` | residues the linker reacts with (e.g. `"KST"`) |
| `protein_nterm` | `true` | allow linking the protein N-terminal amine |
| `site1`, `site2` | none | heterobifunctional chemistry, set both: what each end reacts with, as residue letters plus `[` (protein N-terminus) and `]` (protein C-terminus). Replaces `residues` / `protein_nterm`; overrides a reagent's sites |
| `allow_cterm_residue` | `false` | allow linking a peptide's C-terminal residue when it is not the protein C-terminus |
| `loop_links` | `true` | also search loop-links (both ends of the linker on one peptide) |
| `multivalent` | off | also search species held together by 2 up to this many linkers (`2`–`4`), see [Multivalent crosslinks](#multivalent-crosslinks). Non-cleavable and MS-cleavable linkers; the search space grows steeply, so use it with small protein sets |
| `candidates` | `30` | linear "seed" peptides per spectrum used to find partners |
| `min_chain_peaks` | `3` | minimum matched fragments required on *each* peptide |

Output: `crosslinks.sage.tsv` (+ `results.json`): every candidate match, one per
spectrum, targets and decoys, with their q-values (decoys are needed to estimate
the FDR, and for viewers that filter themselves). `crosslinks_FDRfiltered.sage.tsv`
has the same columns but only the accepted results: target-target matches with
`spectrum_q` <= 0.01 (the CSM count reported in the log). A linear peptide search of the same
spectra also runs first (`"linear_search"`, on unless set to `false`) and writes the usual bottom-up
results: `results.sage.tsv` (PSMs), `peptides.sage.tsv` (peptides at 1% FDR) and `proteins.sage.tsv`
(their protein groups, each with its `protein_group_q`), plus the LFQ/TMT tables and HTML report when
enabled. Protein groups are not cut at 1%: Sage's picked protein FDR cannot go below 1 / (number of
protein groups), about 0.08 for a 13-protein complex. Run as normal: `sage config.json`.

## Required database settings

* **`enzyme.missed_cleavages` must be >= 1 (typically 2).** A crosslinked lysine
  blocks trypsin, so linked peptides contain an *internal* K. Sage's default
  enzyme allows 0 missed cleavages; Sage warns if you leave it that way.
* `database.prefilter` is rejected (partners rarely identify as linear peptides).
* `peptide_min_mass` defaults to 350 Da in crosslink mode (500 Da otherwise) so
  small partner peptides are kept.
* Several enzymes can be combined (co-digestion, e.g. trypsin + Asp-N, common in
  XL-MS): `"enzyme": {"missed_cleavages": 2, "enzymes": ["trypsin", "asp-n"]}`.
  Peptides are cut wherever any listed enzyme cleaves. Names: `trypsin`,
  `trypsin/p`, `lys-c`, `lys-n`, `arg-c`, `asp-n`, `glu-c`, `chymotrypsin`,
  `chymotrypsin-high`, `proalanase`; a custom rule
  `{"cleave_at": "DE", "restrict": "P", "c_terminal": true}` may be listed too.
  The C-terminal-residue rule (`allow_cterm_residue`) assumes enzymes that cut
  *after* the linked residue; with only N-terminal cutters (Asp-N, Lys-N) consider
  setting `allow_cterm_residue: true`.

## Two databases: crosslinks only between them

For complex crosslink searches, e.g. the interface between two complexes or a
bait and its candidate partners, give two FASTA files instead of `fasta`:

```json
"database": {
  "crosslink_databases": {
    "database_1": "C:/db/complex1.fasta",
    "database_2": "C:/db/complex2.fasta"
  }
}
```

* **Linear peptides** are searched in both files as one database (identical
  proteins once), as with one FASTA, including LFQ, TMT, SILAC and labelled
  proteins.
* **Crosslinks** are searched only between a peptide of a Database 1 protein and
  a peptide of a Database 2 protein. Pairs within one database are never scored,
  so they take no part in the FDR, and the search is faster. A peptide of
  proteins in both files counts as either side (so it can also link to itself).
* **Loop-links** (both ends on one peptide) are not searched.
* **Multivalent crosslinks:** a species is searched if at least one of its links
  joins Database 1 to Database 2; its other links may stay within one database
  (e.g. two subunits of complex 1 linked to each other, one of them also to
  complex 2).
* Decoys come from the same database as their target, so the target-decoy
  competition is restricted in the same way.
* `crosslinks.sage.tsv` gets `database_a` / `database_b` columns: `1`, `2` or
  `1+2`.
* With `crosslink_databases` set, `fasta` is ignored (a warning is logged when
  both are given). In `sage-gui`: Files tab → Two-database crosslink search;
  setting either file greys out the FASTA field.

## Heterobifunctional and zero-length reagents (EDC, DMTMM)

```json
"crosslink": { "reagent": "edc" }
```

EDC and DMTMM couple a primary amine to a carboxyl directly, forming an amide
bond and releasing one water: the linker mass is **-18.010565 Da** and nothing
of the reagent remains. The two ends react with different groups:

| | reacts with | spec |
|---|---|---|
| `site1` | Lys, protein N-terminus (amines) | `"K["` |
| `site2` | Asp, Glu, protein C-terminus (carboxyls) | `"DE]"` |

Every reported link joins a `site1` group to a `site2` group (K-D, K-E,
N-term-E, ...); K-K or D-E pairs are never formed. The same holds for loop-links
and for every backend. The protein C-terminal carboxyl is modelled on the last
residue (every y-ion contains it), as the N-terminal amine is on the first.

Any other two-group chemistry is a custom reagent with both sites, e.g. a
hypothetical amine-to-sulfhydryl linker:

```json
"crosslink": { "mass": 96.02113, "site1": "K[", "site2": "C" }
```

Explicit `site1`/`site2` may also be combined with `"reagent": "edc"` to narrow
or extend its chemistry (e.g. `"site1": "KST["`).

## MS-cleavable linkers (DSSO, DSBU, cliXlink)

```json
"crosslink": { "reagent": "dsso" }
```

A cleavable linker breaks in the collision cell, so each peptide appears in the
MS2 spectrum on its own carrying a small remnant (a *stub*). The two stubs of one
reagent differ by a fixed mass, so every peptide shows up as a **doublet**:

| reagent | linker | stubs | doublet spacing |
|---|---|---|---|
| DSSO | 158.00376 | 54.01056 (alkene), 103.99320 (sulfenic acid); also 85.98264 (thiol) | 49.98264 |
| DSBU | 196.08479 | 85.05276, 111.03203 | 25.97927 |
| cliXlink | 267.05653 | arm 1: 54.01056, 103.99320; also 85.98264 / arm 2 (alkyne tag): 163.06333, 213.04596; also 195.03540 | 49.98264 (both arms) |

cliXlink is **asymmetric**: its two arms are different, so a cleavage leaves one
arm's short stub on one peptide and the *other* arm's long stub on its partner.
Its preset is the same as

```json
"crosslink": { "mass": 267.05653, "cleavable": {
  "short": 54.01056, "long": 103.9932, "alt": [85.98264],
  "other": { "short": 163.06333, "long": 213.04596, "alt": [195.0354] } } }
```

and any other asymmetric linker can be described the same way. Each peptide pair
is tried both ways round (A on arm 1 and B on arm 2, and vice versa) and the
better-scoring orientation is kept.

The search exploits this:

1. Find peak pairs spaced by the doublet mass. If `a = m(A) + short`, the partner
   peptide B (mass `M - m(A) - linker`, precursor mass `M`) must show its doublet at
   `m(B) + short` and `m(B) + long` (of the other arm, for an asymmetric linker).
2. Both peptide masses are now known, so each is looked up directly in the
   mass-sorted peptide list with a narrow window. **No pairwise enumeration and no
   wide seed search.**
3. Each peptide is scored independently as a linear peptide with a stub on one
   residue (any of the stubs may be on fragment ions that contain the link). The
   pair's score is the joint hyperscore.

Because candidates come from narrow mass windows, cost grows roughly with database
size rather than quadratically, which is what makes proteome-scale searching practical.

Recommended: `deisotope: true` (assigns charges to the stub-peptide ions; without
it the search tries every charge up to precursor-1 for each peak) and a
generous `max_peaks` (defaults to 300 in cleavable mode, since doublets are
often weak). Loop-links are searched too (`loop_links`), see below.

### Loop-links of MS-cleavable linkers

A loop-link joins two residues of one peptide, so cleaving the linker does not
split the peptide and gives no stub-peptide doublets. It opens the ring
instead, and that shows in the fragment ions:

* ions that contain **one** linked residue (cleaved inside the former ring) carry
  one stub of that residue's arm (DSSO: alkene +54.011, sulfenic acid +103.993 or
  thiol +85.983). A closed ring never gives these ions, which makes them the
  evidence of an opened loop-link;
* ions that contain **both** linked residues carry both stubs: one arm's alkene
  and the other's sulfenic acid (= the linker mass, as when intact) or thiol;
* ions that contain neither are ordinary linear ions.

Peptides whose mass plus the linker matches the precursor are scored with
these ions for every pair of reactive residues (both arm orientations for
cliXlink), and at least `min_linked_peaks` ring ions must match: without them a
linear peptide with any modification of about the linker's mass would pass.
Loop-links get their own FDR group (`decoys / targets`), as in non-cleavable
mode; `linked_a` reports the matched ring ions.

## Multivalent crosslinks

With `"multivalent": n` (2–4; in `sage-gui`, "search multivalent crosslinks"
on the Crosslinking tab) the search also looks for species held together by
two up to `n` linkers, next to the ordinary pairs and loop-links:

* several links between two peptides: `a2b`, and between two copies of one
  peptide, `a2a`;
* chains of three or more peptides: `a1b1c`, `a2b1c`, `a1b1c1d`, copies
  included (`a1a1b`, `a1b1b`);
* branches: one peptide linked to several others, written as its links,
  `a1b,b1c,b1d`.

Letters are distinct peptide sequences, numbers the links between neighbours.
`n` counts the linkers of a species, so a species has at most `n + 1`
peptides. A residue carries one linker; loop-links stay in the main search.

**Search.** The precursor is `M = sum of the peptides + k × linker` for `k`
linkers. The best seeds of the ordinary search (the 10 linear peptides that
explain the most peaks) are combined, with repetition, into every set of
`p - 1` peptides; the last peptide's mass is then fixed and looked up in the
mass-sorted database. Every connected arrangement of the `k` links on the `p`
peptides is scored (3 for three peptides and two links, 16 trees for four
peptides and three links, 125 for five peptides and four links), and for
heterobifunctional reagents each link both ways round.

**Fragment ions.** Take one peptide away from the species: the rest falls into
groups, each hanging on one or more of its residues. A backbone fragment
carries every group whose residues it holds all of, with that group's peptides
and linkers. A fragment that holds only some of a group's residues stays tied
to the rest through the group (as inside a loop-link) and is not observed.
Every peptide needs `min_chain_peaks` matched ions, `min_linked_peaks` of them
carrying a linker; the score is the hyperscore over all peptides.

**MS-cleavable linkers (DSSO, DSBU, cliXlink).** In the collision cell every
linker of the species can break, which frees each peptide with one stub on
each of its link residues (DSSO: alkene, sulfenic acid or thiol). Each
peptide is scored that way: an ion carries one stub (any of its arm's) for
every link residue it holds. The free peptide with its stubs gives a doublet,
or with two links a triplet (`P+A+A`, `P+A+S`, `P+S+S`, spaced 49.98 Da), as in
the pair search: every peptide needs at least one such peak, and one peptide
two. Once every linker has broken, the fragments no longer show which peptides
were linked to which, only how many links each carries. Pieces where only one
linker broke do: in `a1b1c`, breaking the a–b linker frees `b1c` with a stub,
and breaking b–c frees `a1b` with one. Those peaks are counted per arrangement
as `topology_peaks`, and among equally scoring arrangements the one with most
is reported. A `topology_peaks` of 0 means the arrangement is not shown by the
spectrum (always so for `a2b`, where no single broken linker frees anything).
cliXlink's two arms leave different stubs, so each link is tried both ways
round. On the DSSO Pol II run (2 or 3 linkers) 16 spectra were best explained
by a multivalent species and one passed 1% FDR: an RPB2 `a2b`,
SAGNIVQVFKKAAPSPISHVAEIR (K10, K11) with DGKLAKPR (K3, K6), 29 and 13 matched
ions and 6 doublet / triplet peaks, on a spectrum the pair search had left
unexplained; pairs and loop-links were unchanged.

**FDR.** Multivalent species are their own FDR group. A species with any decoy
peptide is a decoy. A random species of `n` distinct peptides is all-target
with chance `1 / 2^n`, so each decoy species stands for `1 / (2^n − 1)` false
targets: `(TD + DD) / 3` for two peptides, which the usual `TD − DD` also
estimates. Copies of one peptide count once (`a2a` has one distinct peptide, so
a decoy `a2a` counts fully, like a decoy loop-link). Unlike an
inclusion–exclusion estimate (`N1 − N2 + …`) this never turns negative when
all-decoy species are common, as they are with copies. The species compete
with pairs and loop-links for each spectrum.

**Output.** `crosslinks.sage.tsv` lists them with `link_type` `multivalent`;
`peptide_a` / `peptide_b` and `link_a` / `link_b` hold the first link, and two
columns describe the whole species:

* `topology`: `a2b`, `a1b1c`, … (`a1b` or `a1a` for ordinary crosslinks, `loop`
  for loop-links);
* `species`: the peptides, numbered, and each link by peptide and residue,
  `[1] GLEKVSTAAMQR [2] AEFGHDSTWKLMNPQSTVKGYEDFAR [3] TDYHFPNEWKLIR | 1:K4-2:K10; 2:K19-3:K10`;
* `topology_peaks`: MS-cleavable linkers, the peaks of parts freed by one broken
  linker that support the arrangement (above). `stub_peaks` holds the free
  peptides' doublet and triplet peaks.

The Annotator CSV writes them as ProForma with one label per link (`#XL1`,
`#XL2`, …), and the viewer annotates them with every peptide and link in the
map, and the stub ions of cleaved linkers (see Annotated spectra). They
are left out of the xiVIEW, ProXL and mzIdentML files, which describe pairs
only, and out of crosslink LFQ. The mobiligram shows them in green (#6FB120).

**Search space.** Each extra peptide multiplies the candidates by about the
number of seeds, and each extra linker adds arrangements, so a search with
`n = 3` or `4` takes far longer than an ordinary one. Use a FASTA of the
complex under study (a few dozen proteins), not a proteome; `sage-gui` warns
while the option is on.

## Choosing a backend

Non-cleavable reagents can be searched with three interchangeable strategies,
selected with `"backend"`. MS-cleavable reagents always use the doublet search on
the CPU (it is already cheap and has nothing worth offloading); the setting is
ignored for them, with a log message.

| backend | what it does | when to use it |
|---|---|---|
| `cpu` (default) | seed search + partner lookup (below) | the validated default; cheap enough for purified complexes and small databases |
| `reflect` | chain reflection (below) on the CPU | a slow *reference* of the GPU algorithm; use it to check the GPU, not for large databases |
| `gpu` | chain reflection with the per-peptide scoring on a GPU | proteome-scale non-cleavable searches. Errors if no GPU backend is available |
| `auto` | `gpu` if a real GPU is available, else `cpu` | convenience. Software adapters (llvmpipe, WARP) are never chosen by `auto`, since they are slower than `cpu`; backends can give slightly different results |

### Building with GPU support

GPU support is built into `sage` by default:

```
cargo build --release -p sage-cli
```

(`--no-default-features` builds without it.) This uses [wgpu](https://wgpu.rs),
so it runs on DirectX 12 (Windows), Vulkan (Linux) and Metal (macOS) with NVIDIA,
AMD, Intel or Apple GPUs; there is no CUDA dependency. Windows uses DirectX 12
only: some Windows Vulkan drivers (Intel UHD 31.0.101.2145) crash when a device
is opened. `sage --gpu-info` lists the adapters a build can use. It needs Rust
1.87 or newer. In a build without GPU support, `"backend": "gpu"` is an error and
`"auto"` falls back to `cpu`. Searches run on the CPU unless the backend is set
to `gpu` or `auto`, so the GPU support costs nothing when it is not used.

The peptide ion table is uploaded once and stays resident; for a human proteome at
2 missed cleavages expect on the order of 1 GB of device memory (Sage reports the
table size at start-up). Tables larger than the device's buffer limits are split
into chunks automatically.

### Checking the GPU on your machine

The GPU kernel is tested against the CPU reference on every adapter of the machine
(hardware and software). Run this on the machine with the GPU; it prints each
adapter and fails if any survivor differs (peptide set, link site, matched ions, or
score beyond 1e-3), including for peaks placed exactly on tolerance-window edges:

```
cargo test -p sage-core --features gpu gpu_kernel -- --nocapture
```

Window edges are computed so that every device rounds them identically: charges
enter as reciprocals rounded once on the CPU (GPU division is not correctly
rounded), and edges are `center * (1 + rel) + abs`, which a fused multiply-add
cannot change (`abs` is 0 for ppm, the factor is 1 for Da). Only the survivor
score passes through the GPU's `log`, hence the 1e-3 allowance. The final scores
are always recomputed exactly on the CPU.

### How chain reflection works

For a pair `A-B` with precursor mass `M = m(A) + m(B) + linker`, every fragment of
`A` containing the linked residue carries the other chain plus the linker, a shift
of `M - m(A)`. That depends only on `A` and the precursor, not on `B`, so:

1. **Chain filter (the GPU part).** Every peptide is scored independently, with
   its link-residue ions shifted by `M - m(A)`; the best link site is kept. The
   `chain_survivors` best peptides go on.
2. **Join on mass.** Survivors are paired where `m(B) = M - linker - m(A)`.
3. **Exact rescoring.** The best `candidates` pairs are rescored with the real
   partner mass, using the same scoring (and FDR) as the `cpu` backend.
4. **Loop-links** need no partner: peptides with `m(A) + linker = M` are looked up
   directly.

The shift uses the *observed* precursor mass, so shifted ions carry the precursor
error; their windows are widened by it. The filter score is therefore a coarse
surrogate (always the Sage hyperscore formula); only the final score is exact.

The top-`chain_survivors` cut makes the join a heuristic, not exact: a true chain
that scores below the cut is lost. True chains usually score far above random
peptides, but raise `chain_survivors` if you suspect misses.

## How the non-cleavable `cpu` search works

1. **Seed search** - the spectrum is searched against the fragment index with a
   very wide precursor window. Each chain of a crosslinked pair still yields many
   ordinary, unshifted b/y ions, so the best-matching linear peptides are seeds.
2. **Partner lookup** - for seed `A`, the partner mass is fixed by the precursor:
   `m(B) = M - m(A) - linker`. Peptides are mass-sorted, so partners are found by
   binary search; each pair, and each pair of reactive residues, is scored.
3. **Fragment model** - for chain `A` linked at residue `i`, every fragment that
   *contains* residue `i` carries the whole other chain plus the linker.
4. **Score** - joint hyperscore over both chains, requiring `min_chain_peaks`
   matched ions on each. `alpha` is the better-supported chain.
5. **FDR** - target-decoy, `(TD - DD) / TT` for crosslinks (Walzthoeni 2012);
   `decoys / targets` for loop-links. Computed **separately** for loop-links,
   intra-protein and inter-protein links. `spectrum_q` is per CSM;
   `crosslink_q` is per unique (peptide pair + sites).

## Output columns

`link_a` / `link_b` are 1-based residue positions *within* `peptide_a` /
`peptide_b`. For a loop-link, both are positions on `peptide_a`. `label` is 1 for
target-target, -1 otherwise; `class` is `TT`/`TD`/`DD` (loop-links: `TT`/`DD`).
`matched_a` / `matched_b` count all matched fragment ions of each peptide;
`linked_a` / `linked_b` count those that carry the linker (non-cleavable linkers
only, see `min_linked_peaks`). With two crosslink databases, `database_a` /
`database_b` say which one each peptide is from (`1`, `2` or `1+2`).

## Files for xiVIEW, ProXL and PRIDE

Every crosslink search also writes three files for viewers and repositories, next
to `crosslinks.sage.tsv`:

* **`crosslinks.mzid`**: the standard format, mzIdentML 1.3.0 following the PSI
  [crosslinking extension 1.0.0](https://github.com/HUPO-PSI/mzIdentML). This is
  what the current xiVIEW takes (including a local installation, whose loader
  reads mzIdentML only) and what PRIDE expects for crosslinking submissions.
  * Each crosslinked species is a donor/acceptor pair of `<Peptide>`s; the donor
    carries the linker mass and the reagent term (XLMOD; UNIMOD `Xlink:EDC` for
    EDC), the acceptor a zero mass. The donor is the end that reacts with the
    first site (EDC: the amine), otherwise the longer peptide.
  * A crosslink-spectrum match is two `<SpectrumIdentificationItem>`s sharing
    "crosslink spectrum identification item"; loop-links are included as
    "looplink spectrum identification item" on a peptide with both ends.
  * All matches and decoys are included (`isDecoy`), `passThreshold` is set at
    crosslink q ≤ 0.01. Scores: `PSM-level q-value` (MS:1002354), and
    `Sage:crosslink q-value`, `Sage:hyperscore`, `Sage:delta next` as user
    parameters (PSI-MS has no Sage score terms).
  * Cleavable linkers carry "crosslinker cleavage characteristics" for their
    stubs (DSSO: `A:54.01056:ST`, `S:103.9932:A`, `T:85.98264:A`, as in the
    specification).
  * Fragment and precursor tolerances are written symmetric (the wider side):
    xiVIEW rejects unequal plus and minus values.
  * Modifications named `UNIMOD:n` for common entries get their UNIMOD term;
    others are "unknown modification" with the name or mass.
  * Spectra are referenced as in the table below (xiVIEW column); upload the
    mzML/MGF files and the FASTA with it. It is not written when there are no
    matches (mzIdentML needs at least one identification).
* **`crosslinks_xiView.csv`** for [xiVIEW](https://github.com/colin-combe/CLMS-UI)
  (to upload to a xiVIEW server; to view results on this computer, see
  [xiVIEW on this computer](#xiview-on-this-computer)),
  in the CSV layout its importer reads (`PepSeq1`, `PepPos1`, `LinkPos1`,
  `Protein1`, `Decoy1`, ... `PeakListFileName`, `ScanId`). Upload it with the
  FASTA (and, to see annotated spectra, the mzML files: `ScanId` is the mzML
  spectrum id). All matches are included; decoys are flagged in `Decoy1`/`Decoy2`
  and `PassThreshold` is `True` at crosslink q ≤ 0.01. Modifications follow their
  residue as a mass, e.g. `M(15.9949)`. Loop-links are left out: the CSV format
  cannot describe them (the log says how many).
* **`crosslinks_proXL.xml`** for [ProXL](https://proxl-ms.org/), valid against the
  ProXL XML schema. It holds rank-1 target-target crosslinks and loop-links grouped
  into reported peptides, with Sage's `crosslink q-value`, `PSM q-value`,
  `hyperscore` and `delta next` as filterable annotations (default filter:
  crosslink q-value ≤ 0.01), the linker (mass, cleaved stub masses, reactive
  ends), static and variable modifications, and the matched protein sequences.
  Upload it with the same FASTA and mzML files. It is not written when there are
  no target matches (ProXL rejects an empty file). The linker is named after the
  reagent (`dss`, `bs3`, `dsg`, `bs2g`, `dsso`, `dsbu`, `clixlink`, `edc`, `dmtmm`; a custom
  mass is `custom`) and always carries its mass and reactive ends.

Peptide positions in proteins come from the FASTA; a decoy peptide is placed
where its target peptide is.

### Viewing spectra alongside the crosslinks

Both viewers read the spectra from the files you searched; upload them together
with the export and the FASTA, **under the names written in the export**:

| Searched file | xiVIEW (`PeakListFileName`, `ScanId`) | ProXL (`scan_file_name`, `scan_number`) |
|---|---|---|
| `run.mzML` | `run.mzML`, the mzML spectrum id | `run.mzML`, the scan number (below) |
| `run.mzML.gz` | `run.mzML` (xiVIEW also accepts the `.gz`) | `run.mzML`: upload it **uncompressed** |
| `run.mgf` | `run.mgf`, `index=N` (position in the file) | none: ProXL does not read MGF |
| Bruker `.d` | neither viewer reads `.d`; convert to mzML for spectra | none |

ProXL numbers mzML scans like its reader (JRAP): the number after `scan=` (Thermo
ids) or `scanId=`, or the id itself if it is a number, otherwise the spectrum's
position in the file + 1. Where the ids alone are not enough (MGF; mzML ids such
as `spectrum=2442`), Sage reads the spectrum positions from the file again after
the search.

These references were checked against the viewers' own readers: pyteomics
`get_by_id` / MGF index lookup as xiVIEW's importer uses them, and JRAP's scan
numbering from ProXL's importer.

## xiVIEW on this computer

`sage-gui` runs [xiVIEW](https://github.com/Rappsilber-Laboratory/xiview) (the
Rappsilber lab's crosslink viewer) on this computer with garlic's results: **Open
in xiVIEW** (bottom bar, Progress tab, Annotated spectra tab) or
`sage-gui --xiview <results folder>` opens it in your web browser. Nothing is
uploaded and no database or Python is needed.

* sage-gui serves the xiVIEW web app (the `xiview` folder next to sage-gui.exe)
  on `127.0.0.1` while it is open, and answers xiVIEW's data requests itself from
  the results folder, in the form of the PRIDE crosslinking API that xiVIEW
  normally reads (matches, peptides, proteins, search settings).
* All of xiVIEW works: network (xiNET), matrix, circular, alignment, 3D structure
  (load a PDB file or fetch one), distance histograms, filters, selection
  tables, exports and the annotated spectra (xiSPEC), including **View Spectra
  Only**. UniProt features and PDB structures are downloaded when those views
  use them; the rest works offline.
* **FDR: garlic's.** xiVIEW gets the target crosslinks within garlic's crosslink
  FDR (crosslink q-value; 1% by default, 5% or 10% chosen next to the button),
  all shown as passing ("Stored validation"). Decoys are not sent, so xiVIEW's
  own decoy-based "In-browser FDR" is not used. Scores: `garlic:hyperscore`
  (the default) and `garlic:q-value`.
* Loop-links are shown on their peptide (in the spectrum viewer, the link as a
  bracket).
* **Multivalent crosslinks:**
  * Network, circle (and matrix, 3D) views: each link of a species is a
    crosslink between its residues (a homodimer's same-residue links as
    homomultimer links), drawn **dash-dotted**. garlic loads link metadata with
    the data: `Multivalent` (the topology, `a2a`, `a1b1a`) and `Multivalent
    species` (the species by its links, e.g. `a1b1a: RPO21 K695–RPB3 K161 +
    RPO21 K695–RPB3 K165`, the same on every link of it). Both are link colour
    schemes: **colour by "Multivalent species"** to give each species its own
    colour (other links grey); they are also shown with the links and can be
    filtered on.
  * Spectrum viewer: every link of a species opens the whole species, every
    peptide a row (peptides 3 and 4 in gold and grey) lined up at its linked
    residues, every link a line (the link of the selected match first). Fragments
    are named by the peptides they carry through the links: `b4+P`, `b14+2P`
    (two more peptides), `P+P+A` (two peptides still joined after one linker
    broke, with the stub). The precursor counts one linker per link.
* Proteins are named by gene (from the FASTA) and looked up in UniProt by their
  accession; decoys are `DECOY_<accession>`.
* **Annotated spectra:** the fragments come from Annotator's engine as in
  garlic's Annotated spectra tab, and MS-cleavable linker stub ions (`P+A`,
  `y6+S`, `b26+T+S`) from garlic's own stub search; they are converted to the
  format xiSPEC reads, with the isotope peaks of each matched ion grouped. In a
  sample of 49 spectra from the Pol II DSSO runs (crosslinks, loop-links and
  multivalent species), every peak garlic's viewer explains (a/b/y, precursor
  and stub ions) was annotated in xiSPEC too (xiSPEC also marks their isotope
  peaks). Changing xiSPEC's settings
  (tolerance, ion types, modifications) re-annotates the spectrum locally.
* Spectra are read as Sage read them (timsTOF `.d` included, by spectrum id):
  timsTOF runs are opened in the background and kept open, mzML / MGF files read
  once, so a spectrum opens in a fraction of a second. **Keep sage-gui open**
  while you use xiVIEW: it serves the data and spectra; if it was closed, a
  spectrum click says so (open the results again with Open in xiVIEW).
* **Credit and citation.** xiVIEW, with xiNET and xiSPEC, is the work of the
  [Rappsilber Laboratory](https://www.rappsilberlab.org/) (Colin Combe, Martin
  Graham, Lars Kolbowski, Lutz Fischer, Juri Rappsilber): please cite xiVIEW
  (Combe et al., *J. Mol. Biol.* 2024, doi:10.1016/j.jmb.2024.168656), and
  xiSPEC (Kolbowski et al., *Nucleic Acids Res.* 2018, doi:10.1093/nar/gky353)
  and xiNET (Combe et al., *Mol. Cell. Proteomics* 2015,
  doi:10.1074/mcp.O114.042259) when you use them. Repositories:
  [xiview](https://github.com/Rappsilber-Laboratory/xiview) (GPL-3.0),
  [xiview-server](https://github.com/Rappsilber-Laboratory/xiview-server)
  (Apache-2.0); the data API garlic answers is the PRIDE Crosslinking
  Resource's, [crosslinking-api](https://github.com/PRIDE-Archive/crosslinking-api)
  (PRIDE team, EMBL-EBI, and the Rappsilber Laboratory; data model from
  [mzidentml-reader](https://github.com/PRIDE-Archive/mzidentml-reader)), and the
  annotation format xiSPEC reads is
  [xiAnnotator](https://github.com/Rappsilber-Laboratory/xiAnnotator)'s. See the
  README's Citation and Acknowledgments.
* xiVIEW (GPL-3.0) stays a separate program served to the browser; it is built
  from source at a pinned commit by `xiview/build_xiview.py`, with the changes
  listed in `xiview/source_fixes.py`: the spectrum viewer's settings (their
  re-annotation was dropped) and loop-links (their spectrum did not open), and
  the multivalent display above (any number of peptides and links in the
  spectrum viewer, garlic's link metadata, dash-dotted links). See the
  `SOURCE.md` in the `xiview` folder.

## Annotated spectra

`sage-gui` shows the annotated spectrum of any match on its **Annotated
spectra** tab, with the engine of [Annotator](https://github.com/snijderlab/annotator)
(Snijder lab) built in: `mzcore` parses the ProForma peptidoform, `mzannotate`
generates and matches its fragments. Nothing else needs to be installed.

* Open it with **View spectra** after a search (bottom bar or Progress tab), by
  loading any results folder on the tab, or with `sage-gui --view <folder>`.
* The list holds the matches of `results.sage.tsv` and `crosslinks.sage.tsv`,
  filtered by kind, targets, spectrum q-value and a text search (peptide,
  protein, scan), best score first.
* For the match picked: its spectrum with each peak coloured by ion series and
  labelled (`y7`, `b5²+`, `y4-H2O`; in a crosslink `1:` / `2:` name the
  peptide), zoom and pan, hover for m/z, mass error and every interpretation; the
  backbone ions seen on each peptide; and a table of matched fragments. The
  fragmentation model (CID / HCD by default, with DSSO / DSBU stub ions; ETD,
  EThcD, EAD, All) and the search's fragment tolerance are Annotator's.
* Under the spectrum, each peptide's residues with the backbone ions seen: a hook
  over the right of a residue for an N-terminal (a/b/c) ion ending there, one
  under its left for a C-terminal (x/y/z) ion starting there (b blue, y red).
  A crosslinked pair is drawn as two rows with the linked residues aligned and
  joined by a **bridge** labelled with the linker; a loop-link as a bracket
  under its peptide.
  Left of each crosslinked peptide: `peptide 1 · RPO21` (the gene, else the
  UniProt entry name), the protein's name, and the link residues' positions in
  the protein (`K1246`), read from the search's FASTA. A peptide in more
  proteins shows `+ n`; hover for every protein and position. Decoys show
  their target's positions.
* MS-cleavable linkers: ions carrying a **stub** are found for each peptide
  (the peptide alone with the stub on its linked residue, fragmented by
  Annotator's engine) and shown in the stub's colour (DSSO: A alkene teal, T
  thiol amber, S sulfenic acid magenta; DSBU: Bu, BuUr; cliXlink's alkyne arm
  A', T', S'):
  * the whole peptide + stub (the doublet peaks, at every charge below the
    precursor's) as a triangle on the bridge pointing at that peptide;
  * fragment ions holding the linked residue + stub as hooks with the same
    directionality as the plain series, nested outside them, one per stub;
  * on the peaks as `y6+S`, `P+A²+`; hover the map for every stub ion with its
    observed m/z and error.

  On the DSSO Pol II run, both peptides showed peptide + stub ions in 39 of 40
  crosslink spectra (9 of 40 with the same peptidoforms on other spectra).
* Multivalent species: every peptide is a row, its link residues shaded, and
  each link a bridge labelled `XL1`, `XL2`, … as in the ProForma (rows are
  ordered so that linked peptides are neighbours where possible; a link that
  skips a row runs along a lane right of the rows). MS-cleavable linkers, as
  the search scores them: once every linker broke each peptide keeps one stub
  per link residue, so
  * a peptide with one link shows its peptide + stub peaks as triangles at its
    end of that bridge, and one with two or more links its peptide + stubs
    peaks (`P+A+S`, the triplets) after its row;
  * fragment ions carry one stub per link residue they hold (`y8+S`,
    `y17+A+S²+`), drawn as hooks in each stub's colour;
  * the parts still held together when only one linker broke (`{2,3}+A` from
    `a1b1c`), which show which peptides were linked to which, are triangles at
    that end of that bridge.

  Hover the map for every stub ion. With an asymmetric linker (cliXlink)
  either arm's stubs are looked for at every end, as the species does not
  record which way round each linker sits.
* The header gives charge, scan, score, q-value, retention time, mass and mass
  error, the precursor m/z and, for timsTOF (Bruker `.d`) data, the precursor's
  ion mobility 1/K0 (V·s/cm²).
* **Export figure** saves the shown match as a publication figure: the spectrum
  (the current zoom, or the whole spectrum), labelled as in the viewer
  (internal fragments only with "label all matched peaks"), the peptide map with
  the bridge and stub glyphs, a legend, and optionally a title with the
  peptides, charge, scan, retention time and ion mobility. Formats: SVG and PDF
  (vector, editable in Illustrator / Inkscape), PNG and TIFF at 300, 600 or
  1200 dpi; 7 inches wide (two-column width); white or transparent background
  (PNG and TIFF then carry an alpha channel; opaque figures have none).
* **Mobiligram…** (timsTOF data) opens the whole run as a mobiligram: every
  precursor Sage searched, by m/z and ion mobility 1/K0. The marker shape gives
  the charge (1+ plus, 2+ circle, 3+ triangle, 4+ square, 5+ diamond, 6+ and up
  downward triangle) and the colour the identification at 1% spectrum FDR:
  crosslinks between proteins (red, #FF523B), within a protein (orange, #FF9D35) and loop-links
  (pink) drawn on top with a dark outline, linear peptides (blue, 50%
  transparent) and precursors not identified (grey, can be hidden). One spectrum is one
  point: a crosslink beats a linear match of the same spectrum. Hover for m/z
  and 1/K0; click a legend entry to hide that class; with several runs, pick
  one or show all. **Save…** writes it as SVG or PDF (7 inches wide, white or
  transparent, optional title) with a legend of the colours (with counts) and
  the charge shapes. Needs `crosslinks.sage.tsv` from garlic 0.18 or later,
  which has an `ion_mobility` column. Precursors where Sage found no candidate
  peptide at all are not in Sage's result tables and so not shown.
  A search run from `sage-gui` on timsTOF data also writes mobiligrams by
  itself, into a `Mobiligrams` folder inside the results folder, as SVG and
  PDF (all precursors, white background, titled): one per MS file, named
  after it with `_MG` (`Bruker.d` → `Bruker_MG.svg`), and, when two or more
  MS files were searched, `Combined_Mobiligram.svg` / `.pdf` with all of them
  together (a single MS file gets only its `_MG` figure, which would show the
  same). Save… in the viewer suggests the same names. MS files without ion
  mobility (anything but timsTOF) get none.
* Spectra are read with Sage's own reader by Sage's spectrum ids, exactly as
  searched: mzML, MGF and, for Bruker `.d`, the precursor's spectrum merged over
  all its frames. The matches' peptidoforms are written as for the Annotator
  files below.
* Spectra are read in the background as they are needed: the match picked
  first, then the rows in view of runs already open. Bruker `.d` spectra are
  read one by one by their id, the same way the whole-run read gets each of
  them (identical peaks), so a spectrum shows within milliseconds. Opening a
  run takes about a second, so the runs of the shown matches are opened in
  the background when a results folder is loaded, two at a time, the best
  match's run first ("opening runs 5 of 21"); a run picked before its turn is
  opened right away. Each open run holds about 30 MB. mzML and MGF files have no
  index by spectrum id: each is read whole once, keeping the spectra of every
  match shown, the picked match's file first. The results folder itself is
  read in the background too.

On the DSSO Pol II timsTOF run, the annotated ions carry 31% of the intensity of
the matched spectra of linear PSMs and 18 to 22% for crosslinks, against 6 to 9%
when the same peptidoforms are put on other spectra.

### Files for the Annotator app

For the separate Annotator program, `"write_annotator": true` (or
`--write-annotator`; "write .annotator.csv files" in `sage-gui`) writes the
matches at 1% spectrum FDR (targets, rank 1):

* `results.annotator.csv`: linear PSMs (also from `"linear_search": true`)
* `crosslinks.annotator.csv`: target-target crosslink-spectrum matches

Both are Annotator's "basic CSV" PSM format, `scan_index,sequence,raw_file,z,mode`:

```
scan_index,sequence,raw_file,z,mode
134649,FISTLQVK[XLMOD:02126#XL1]LYGR//K[#XL1]AAPSPISHVAEIR,C:\data\PolII.d,4,
```

* `sequence` is ProForma 2.0: a crosslink is `A[linker#XL1]...//B[#XL1]...`, a
  loop-link puts both ends on one peptide. The linker is its XLMOD term (DSS
  `XLMOD:02001`, BS3 `02000`, DSG `02006`, BS2G `02005`, DSSO `02126`, DSBU
  `02120`, DMTMM `02208`); EDC is `Formula:H-2O-1`, cliXlink and custom linkers
  their mass (Annotator's copy of XLMOD has no mass for cliXlink). A protein
  N-terminal link is written on the first residue: Annotator rejects crosslinks
  on a terminus, and every b ion contains both. Modifications keep their name
  when it is an ontology accession (`UNIMOD:35`), otherwise their mass.
* `scan_index` is the spectrum's 0-based position in its file as Annotator's
  reader (mzdata) numbers it: the mzML `index`, the MGF order, or for Bruker
  `.d` (DDA-PASEF) the position of the precursor's first PASEF window, computed
  from `analysis.tdf` (MS1 frames and every PASEF window count). Sage merges all
  frames of a precursor, so for a precursor fragmented in several frames
  Annotator shows the first frame's spectrum only.
* `mode` is empty: choose the fragmentation model in Annotator (its CID model
  includes DSSO / DSBU stub ions).

Open the spectrum files together with the CSV; Annotator matches each row's
`raw_file` against the open files by name:

```
annotator.exe run1.mzML run2.mzML results.annotator.csv crosslinks.annotator.csv
```

Annotator does not open Bruker `.d` folders from its command line: load the `.d`
with its raw file button, then the matches find their spectra.

Checked with Annotator's own libraries (mzident 0.2.1, mzcore, mzannotate,
mzdata 0.65): every row of a DSSO search of a timsTOF Pol II run (104 crosslinks,
9,618 PSMs) parses, its position holds the spectrum with Sage's precursor m/z,
and the expected fragments annotate (about twice the peaks of the neighbouring
spectrum); synthetic DSS, loop-link and EDC spectra from Sage's fragment model
annotate 98 to 100%.

## Quantification

The usual `quant` block works for crosslink searches:

```json
"quant": { "tmt": "Tmt16", "tmt_settings": { "level": 3 }, "lfq": true,
           "lfq_settings": { "peptide_q_value": 0.01, "ppm_tolerance": 5.0 } }
```

**Isobaric labels (TMT).** Reporter ions are read from each crosslink match's
spectrum (MS2, or with `"level": 3` the SPS-MS3 scan that names it as precursor)
and appended to `crosslinks.sage.tsv` as `ion_injection_time`, `tmt_1` ... (empty
when a spectrum has no reporter scan). `tmt_settings.sn` works as for linear
searches. Sum reporters per unique crosslink (or per protein pair) downstream.

Labelling usually follows crosslinking, so a **crosslinked lysine carries no
label**. A residue with a static modification can't be linked, so a static TMT on
K would block every lysine link (Sage warns). Instead:

```json
"static_mods":   { "^": 304.2071 },
"variable_mods": { "K": [304.2071] }
```

Sage then links only the unlabelled form of a lysine while the other lysines of
the peptide may carry the label.

**SILAC.** With SILAC channels (`database.silac`) a labelled lysine still links,
a crosslink joins peptides of one channel (unless `mixed_crosslinks`), and
`crosslinks.silac.tsv` gives each crosslinked species' intensity in every
channel and its ratios to light per run. See [SILAC.md](SILAC.md).

**Isotope-labelled proteins.** Proteins expressed in ¹⁵N and/or ¹³C media
(`database.labelled_proteins`) are searched labelled; a labelled lysine still
links, and crosslinks between labelled and unlabelled proteins are found
(`isotope_label_a` / `_b` in `crosslinks.sage.tsv`). See [SILAC.md](SILAC.md#isotope-labelled-proteins-15n-13c).

**Label-free (MS1).** With `"lfq": true` Sage writes `crosslinks.lfq.tsv`, one
row per crosslinked *species* (peptide pair or loop-linked peptide; link-site
isomers have the same MS1 signal and are not separated) with its intensity in
every run:

1. Species with a target-target match at `crosslink_q <= lfq_settings.peptide_q_value`
   are traced in MS1 at every precursor charge and the first three isotopes,
   within `ppm_tolerance`, using the isotope envelope of both peptides plus the
   linker (`mobility_pct_tolerance` on ion-mobility data).
2. Runs are aligned with a linear retention-time fit on species identified in at
   least two runs (at least three such species per run; otherwise unaligned).
   Peaks are then smoothed, time-warped and integrated exactly as in linear LFQ
   (`peak_scoring`, `integration`, `spectral_angle`, `combine_charge_states`).
3. Each species also gets a decoy trace (shifted m/z, earlier retention time);
   `q_value` is the MS1 peak q-value from decoys / targets. Rows are named by the
   species' best match (`peptide_a`, `link_a`, ... as above); `charge` is -1 when
   charge states are combined.

`predict_rt` is not used: it needs linear PSMs.

## Limitations

* `reflect`/`gpu`: results are close to, but not bit-identical with, the `cpu`
  backend (different candidate generation). Report which backend you used.
* Cleavable mode needs the stub doublets to be present in the MS2 spectrum
  (i.e. HCD/stepped-HCD data). It is not an MS2-MS3 workflow.
* Cleavable mode: peptides carrying a *mono-link* are not searched. Loop-links are
  found through the ring ions of the opened linker (see "Loop-links of
  MS-cleavable linkers"), so the two linked residues need sequence between them.
* **Mono-links** (hydrolyzed/amidated dead-ends) are not part of this search;
  search for them as an ordinary variable modification (`K: 156.0786`).
* Loop-links are found through unshifted ions *outside* the ring, so the linked
  residues need flanking sequence.
* Link positions are reported within peptides, not protein coordinates.
* No machine-learning rescoring yet (Sage's LDA is tied to linear PSMs); FDR is
  on the hyperscore. Adjacent reactive residues (e.g. `KK`) can be ambiguous.
* `wide_window`/DIA, `chimera`, and parquet output are not supported in
  crosslink mode.
* Crosslink LFQ aligns runs only on crosslinks (a few per run in small studies);
  with few shared identifications runs stay unaligned and rely on the local time
  warping within the +/-0.5% retention-time window. Match-between-runs works for
  species identified in at least one run.
* **Validated on synthetic data only** (unit tests in `crosslink.rs`, plus an
  independent generator). On synthetic DSSO/DSBU data an entrapment test (spectra
  of peptides absent from the FASTA) gave an observed FDR close to the estimate,
  but with very few false hits. Real spectra have chimeric precursors, missing and
  mis-assigned peaks, and non-uniform fragmentation that synthetic data lack.
  Benchmark on real data before trusting the FDR.
