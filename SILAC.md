# SILAC

garlic searches and quantifies SILAC samples (stable isotope labelling by amino
acids in cell culture), duplex or triplex, in linear and in crosslink searches.

## Configuration

Labelled channels go under `database` (light, unlabelled, is implicit):

```json
"database": {
  "silac": {
    "channels": [
      { "name": "medium", "labels": { "K": 4.025107, "R": 6.020129 } },
      { "name": "heavy",  "labels": { "K": 8.014199, "R": 10.008269 } }
    ],
    "mixed_crosslinks": false
  }
}
```

* `channels`: one to four labelled channels, each a name and its label mass per
  residue (one-letter codes; any residue, usually K and R). Usual labels:

  | Label | Mass (Da) | Written as |
  |---|---|---|
  | Lys4 (²H₄) | 4.025107 | `UNIMOD:481` |
  | Lys6 / Arg6 (¹³C₆) | 6.020129 | `UNIMOD:188` |
  | Lys8 (¹³C₆ ¹⁵N₂) | 8.014199 | `UNIMOD:259` |
  | Arg10 (¹³C₆ ¹⁵N₄) | 10.008269 | `UNIMOD:267` |

  Other masses are written as the mass (`K[+6.0138]`).
* `mixed_crosslinks` (crosslink searches): `false` (default) pairs only peptides
  of one channel; `true` also crosslinks between channels, e.g. when lysates
  were mixed before crosslinking.
* MS1 tracing uses the LFQ settings under `quant.lfq_settings`: `ppm_tolerance`,
  `mobility_pct_tolerance`, `spectral_angle`, and `peptide_q_value` as the
  identification cutoff. `"lfq": true` is not needed.

In `sage-gui`: Scoring & output tab → Quantification → SILAC (off, duplex,
triplex), with presets for each label or any mass.

## Search

* Every peptide is searched light and in each channel, with **all** of its
  labelled residues carrying that channel's label. A peptide comes from one
  population of cells, so mixed label states (one K heavy, another light) are
  not generated. The search space, and time, grows 2× (duplex) or 3× (triplex).
* Labels are registered modifications, so precursor and fragment masses, the
  exports (mzIdentML, xiVIEW, ProXL, Annotator CSV) and the Annotated spectra
  viewer all carry them: `PEPTIDEK[UNIMOD:259]`. A label on a residue that has
  another modification (acetyl-K, a static TMT) is one combined modification of
  both masses.
* A peptide without labelled residues (e.g. a protein C-terminal peptide with no
  K or R) looks the same in every channel; it has no channel and is not
  quantified.
* `results.sage.tsv` and `crosslinks.sage.tsv` get a `silac_channel` column
  (`light`, `heavy`, …; empty without labelled residues; `mixed` for a crosslink
  between channels).

### Crosslinks

* A SILAC label is an isotope substitution: a labelled lysine still reacts and is
  a link site (a *modified* lysine is not).
* A crosslink (and a multivalent species) joins peptides of one channel, unless
  `mixed_crosslinks` is on. All backends (cpu, reflect, gpu) and the cleavable
  and multivalent searches apply this.

## Quantification

Every identified peptide (target, rank 1, spectrum and peptide q-value at or
below `peptide_q_value`) and crosslinked species (target-target, crosslink
q-value at or below it) is traced in MS1 **in every channel at once**: its
first three isotopes at every precursor charge, around the aligned retention
time of its best identification, in every run. The traces are smoothed and
time-warped as in LFQ, and all channels and runs are integrated over **one
shared peak**, so the ratios compare like with like. A species identified in
one channel only is quantified in the others too. Runs are aligned on species
identified in several runs.

Tables (needs MS1 spectra in the files):

* `silac.sage.tsv`: one row per peptide in its light form (labels removed, other
  modifications kept): `proteins`, `identified_in` (channels it was identified
  in), `psms`, `peak_score`, `spectral_angle`, then the intensity of each channel
  in each run (`<run> light`, `<run> heavy`) and each labelled channel's ratio to
  light per run (`<run> heavy/light`; empty when either is zero).
* `silac_proteins.sage.tsv`: per protein, the median ratio of its unique
  peptides (peptides of that protein only) per run and channel, and how many
  peptides gave a ratio. Shared peptides are left out.
* `crosslinks.silac.tsv`: the same per crosslinked species (`link_type`,
  `peptides` and `proteins` of each peptide); link-site isomers have the same MS1
  signal and share a row.

Ratios are not normalised: apply your own normalisation (e.g. median-centring
on all peptides) for unequal mixing.

## Checks and limitations

* Tested on synthetic data: label states and masses, labelled lysines linking
  and same-channel pairing (cpu and reflect backends), a 1 : 2 light : heavy
  MS1 trace read as a ratio of 2.0, and Annotator resolving the labels.
* Negative control: one unlabelled Pol II DSSO run (timsTOF) searched as a
  duplex. Crosslink matches at 1% were all light; heavy / light peptide ratios
  had a median of 0.014 (protein medians 0.001 to 0.10) and crosslinked species
  about 0.001, as they should without a heavy sample. 2.9% of the linear PSMs at
  1% were called heavy; 44 of those 46 were long peptides whose light form was
  identified at the same time, typical of a precursor mass read whole isotopes
  off (their b ions are the same in both states). The MS1 ratios of those
  peptides were unaffected. Expect such heavy calls in light-only data and
  check heavy identifications against their MS1 ratio.
* **Not yet tested on a real SILAC sample.** Check the ratios of a known mixture
  before relying on them.
* Arginine-to-proline conversion and incomplete labelling are not modelled.
