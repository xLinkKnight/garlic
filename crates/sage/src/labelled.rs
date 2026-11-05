//! Uniformly isotope-labelled proteins: proteins expressed in ¹⁵N and/or ¹³C
//! media, every nitrogen and/or carbon atom of their residues heavy.
//!
//! The chosen proteins' peptides are searched labelled: each residue carries
//! its own mass shift (its carbon and nitrogen atoms times the isotope
//! difference) as a registered modification, applied after the fixed and
//! variable modifications, so labels never count towards `max_variable_mods`
//! and a residue keeps its fixed or variable modification (label and
//! modification together are one combined modification; a modification's own
//! atoms come from unlabelled reagents and stay light). Precursor and
//! fragment masses, crosslink candidates, every export and the viewer carry
//! the labels. `include_unlabelled` also searches the chosen proteins
//! unlabelled; peptides of other proteins are searched as usual. A labelled
//! residue is still free to crosslink.

use crate::modification::{ModId, ModificationError, ModificationSpecificity, Modifications};
use crate::peptide::Peptide;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

/// Mass of ¹³C less ¹²C, and of ¹⁵N less ¹⁴N
pub const C13: f32 = 1.003_354_8;
pub const N15: f32 = 0.997_034_9;

/// The atoms a protein is labelled in
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Isotope {
    #[serde(rename = "15N", alias = "N15")]
    N15,
    #[serde(rename = "13C", alias = "C13")]
    C13,
    #[serde(rename = "13C15N", alias = "15N13C", alias = "13C+15N", alias = "C13N15")]
    C13N15,
}

impl Isotope {
    pub const ALL: [Isotope; 3] = [Isotope::N15, Isotope::C13, Isotope::C13N15];

    pub fn name(self) -> &'static str {
        match self {
            Isotope::N15 => "15N",
            Isotope::C13 => "13C",
            Isotope::C13N15 => "13C15N",
        }
    }

    fn labels_nitrogen(self) -> bool {
        matches!(self, Isotope::N15 | Isotope::C13N15)
    }

    /// Mass shift of residue `residue` labelled uniformly
    pub fn shift(self, residue: u8) -> Option<f32> {
        let (c, n) = atoms(residue)?;
        Some(match self {
            Isotope::N15 => n as f32 * N15,
            Isotope::C13 => c as f32 * C13,
            Isotope::C13N15 => c as f32 * C13 + n as f32 * N15,
        })
    }
}

/// Carbon and nitrogen atoms of a residue (an amino acid less water)
pub fn atoms(residue: u8) -> Option<(u8, u8)> {
    Some(match residue {
        b'A' => (3, 1),
        b'R' => (6, 4),
        b'N' => (4, 2),
        b'D' => (4, 1),
        b'C' => (3, 1),
        b'E' => (5, 1),
        b'Q' => (5, 2),
        b'G' => (2, 1),
        b'H' => (6, 3),
        b'I' => (6, 1),
        b'L' => (6, 1),
        b'K' => (6, 2),
        b'M' => (5, 1),
        b'F' => (9, 1),
        b'P' => (5, 1),
        b'S' => (3, 1),
        b'T' => (4, 1),
        b'W' => (11, 2),
        b'Y' => (9, 1),
        b'V' => (5, 1),
        b'U' => (3, 1),
        b'O' => (12, 3),
        _ => return None,
    })
}

/// `database.labelled_proteins`
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct LabelledProteins {
    /// Protein → label. A protein is named by its accession as Sage reads it
    /// (the header's first word, `sp|P04050|RPB1_YEAST`) or by one of its
    /// `|`-separated parts (`P04050`, `RPB1_YEAST`)
    pub proteins: BTreeMap<String, Isotope>,
    /// Also search the chosen proteins unlabelled
    #[serde(default = "default_true")]
    pub include_unlabelled: bool,
}

fn default_true() -> bool {
    true
}

impl LabelledProteins {
    pub fn validate(&self) -> Result<(), String> {
        if self.proteins.is_empty() {
            return Err("labelled_proteins: name at least one protein".into());
        }
        if self.proteins.keys().any(|k| k.trim().is_empty()) {
            return Err("labelled_proteins: a protein name is empty".into());
        }
        Ok(())
    }
}

/// Mass a modification loses from a labelled residue: modifications that
/// remove one of its labelled atoms (deamidation of N/Q, pyro-Glu from Q
/// lose a nitrogen) are lighter by its label
fn removed(mass: f32, residue: u8, isotope: Isotope) -> f32 {
    let near = |x: f32| (mass - x).abs() < 0.002;
    let loses_nitrogen = (matches!(residue, b'N' | b'Q') && near(0.984_016)) || (residue == b'Q' && near(-17.026_548));
    if isotope.labels_nitrogen() && loses_nitrogen {
        N15
    } else {
        0.0
    }
}

/// Can a modification of `specificity` sit on residue `residue`?
fn on_residue(specificity: ModificationSpecificity, residue: u8) -> bool {
    use ModificationSpecificity::*;
    match specificity {
        Residue(r) => r == residue,
        PeptideN(Some(r)) | PeptideC(Some(r)) | ProteinN(Some(r)) | ProteinC(Some(r)) => r == residue,
        _ => false,
    }
}

/// The labels as modification ids
#[derive(Clone, Debug, Default)]
pub struct ProteinLabels {
    pub settings: Option<LabelledProteins>,
    /// (isotope, residue) → (label, mass)
    pure: HashMap<(Isotope, u8), (ModId, f32)>,
    /// (isotope, residue, modification already there) → (both, mass added)
    combined: HashMap<(Isotope, u8, ModId), (ModId, f32)>,
    /// (modification, residue) → isotope
    reverse: HashMap<(ModId, u8), Isotope>,
}

impl ProteinLabels {
    /// Register the labels of every isotope in `settings` on every residue,
    /// and each configured modification of a residue combined with them
    pub fn register(
        settings: &LabelledProteins,
        registry: &mut Modifications,
        configured: &[(ModificationSpecificity, ModId)],
    ) -> Result<ProteinLabels, ModificationError> {
        let mut labels = ProteinLabels { settings: Some(settings.clone()), ..Default::default() };
        let mut isotopes: Vec<Isotope> = settings.proteins.values().copied().collect();
        isotopes.sort();
        isotopes.dedup();
        for isotope in isotopes {
            for &residue in crate::mass::VALID_AA.iter() {
                let Some(shift) = isotope.shift(residue) else { continue };
                let id = registry.insert(shift, None)?;
                labels.pure.insert((isotope, residue), (id, shift));
                labels.reverse.insert((id, residue), isotope);
                for &(specificity, base) in configured {
                    if !on_residue(specificity, residue) {
                        continue;
                    }
                    let base_mass = registry.get(base).mass;
                    let added = shift - removed(base_mass, residue, isotope);
                    let both = registry.insert(base_mass + added, None)?;
                    labels.combined.insert((isotope, residue, base), (both, added));
                    labels.reverse.insert((both, residue), isotope);
                }
            }
        }
        Ok(labels)
    }

    /// The label of a protein, by its accession or a `|`-separated part of it
    pub fn label_of(&self, accession: &str) -> Option<Isotope> {
        let proteins = &self.settings.as_ref()?.proteins;
        proteins
            .get(accession)
            .or_else(|| accession.split('|').skip(1).find_map(|part| proteins.get(part)))
            .copied()
    }

    /// The label modifications themselves (a labelled residue without
    /// another modification), and the residue each sits on
    pub fn label_ids(&self) -> Vec<(u8, ModId)> {
        let mut ids: Vec<(u8, ModId)> = self.pure.iter().map(|(&(_, r), &(id, _))| (r, id)).collect();
        ids.sort();
        ids.dedup();
        ids
    }

    /// Every label modification, with the residues it sits on (for
    /// registering SILAC labels on top of them)
    pub fn all_ids(&self) -> Vec<(u8, ModId)> {
        let mut ids: Vec<(u8, ModId)> = self.reverse.keys().map(|&(id, r)| (r, id)).collect();
        ids.sort();
        ids.dedup();
        ids
    }

    /// `peptide` labelled with `isotope` on every residue; `None` if a
    /// residue carries a modification without a combined form
    fn label(&self, peptide: &Peptide, isotope: Isotope) -> Option<Peptide> {
        let mut p = peptide.clone();
        for i in 0..p.sequence.len() {
            let aa = p.sequence[i];
            let (id, shift) = match p.modifications[i] {
                None => *self.pure.get(&(isotope, aa))?,
                Some(m) => *self.combined.get(&(isotope, aa, m))?,
            };
            p.modifications[i] = Some(id);
            p.monoisotopic += shift;
        }
        Some(p)
    }

    /// The forms of `peptide` to search: labelled with each isotope of the
    /// chosen proteins it comes from (credited to those proteins), and
    /// unlabelled for its other proteins (and, with `include_unlabelled`, the
    /// chosen ones too)
    pub fn states(&self, peptide: Peptide) -> Vec<Peptide> {
        let include_unlabelled = self.settings.as_ref().is_some_and(|s| s.include_unlabelled);
        let mut plain: Vec<Arc<str>> = Vec::new();
        let mut labelled: BTreeMap<Isotope, Vec<Arc<str>>> = BTreeMap::new();
        for protein in &peptide.proteins {
            match self.label_of(protein) {
                Some(isotope) => labelled.entry(isotope).or_default().push(protein.clone()),
                None => plain.push(protein.clone()),
            }
        }
        if labelled.is_empty() {
            return vec![peptide];
        }
        let mut out = Vec::with_capacity(labelled.len() + 1);
        for (&isotope, proteins) in &labelled {
            if let Some(mut p) = self.label(&peptide, isotope) {
                p.proteins = proteins.clone();
                out.push(p);
            }
        }
        if !plain.is_empty() || include_unlabelled {
            let mut p = peptide;
            if !include_unlabelled {
                p.proteins = plain;
            }
            out.insert(0, p);
        }
        out
    }

    /// The label of a peptide: the isotope of its residues' labels
    pub fn isotope(&self, peptide: &Peptide) -> Option<Isotope> {
        peptide
            .sequence
            .iter()
            .zip(peptide.modifications.iter())
            .find_map(|(aa, m)| m.and_then(|m| self.reverse.get(&(m, *aa)).copied()))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::Builder;

    #[test]
    fn residue_shifts() {
        // a uniformly 13C15N lysine is Lys8, arginine Arg10
        assert!((Isotope::C13N15.shift(b'K').unwrap() - 8.014199).abs() < 1e-4);
        assert!((Isotope::C13N15.shift(b'R').unwrap() - 10.008269).abs() < 1e-4);
        assert!((Isotope::C13.shift(b'K').unwrap() - 6.020129).abs() < 1e-4);
        assert!((Isotope::N15.shift(b'W').unwrap() - 2.0 * N15).abs() < 1e-5);
        assert_eq!(Isotope::N15.shift(b'X'), None);
        let parsed: LabelledProteins = serde_json::from_value(serde_json::json!({"proteins": {"P1": "15N", "P2": "13C+15N"}})).unwrap();
        assert_eq!(parsed.proteins["P2"], Isotope::C13N15);
        assert!(parsed.include_unlabelled, "on unless turned off");
    }

    fn digested(include_unlabelled: bool) -> (crate::database::Parameters, Vec<Peptide>) {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": {"C": 57.021465},
            "variable_mods": {"M": [15.994915], "N": [0.984016]},
            "max_variable_mods": 1,
            "labelled_proteins": {"proteins": {"P1": "15N"}, "include_unlabelled": include_unlabelled}
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        // LGCMNPEK is in both proteins; AAGHWLK only in the labelled one, DDGHWLK only in the other
        let fasta = crate::fasta::Fasta::parse(
            ">sp|P1|ONE_HUMAN labelled\nLGCMNPEKAAGHWLK\n>sp|P2|TWO_HUMAN plain\nLGCMNPEKDDGHWLK\n".into(),
            "rev_",
            false,
        );
        let peptides = params.digest(&fasta);
        (params, peptides)
    }

    /// The chosen protein's peptides are labelled on every residue, credited
    /// to it; a shared peptide is also unlabelled for the other protein;
    /// fixed and variable modifications stay, and labels do not count
    /// towards max_variable_mods
    #[test]
    fn labelled_forms_of_peptides() {
        let (params, peptides) = digested(false);
        let labels = params.protein_labels.clone().unwrap();
        let forms = |seq: &[u8]| -> Vec<&Peptide> { peptides.iter().filter(|p| p.sequence.as_ref() == seq).collect() };
        let proteins = |p: &Peptide| p.proteins.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(";");

        // only in the labelled protein: labelled only (unlabelled form off)
        let aag = forms(b"AAGHWLK");
        assert_eq!(aag.len(), 1);
        assert_eq!(labels.isotope(aag[0]), Some(Isotope::N15));
        let n_atoms: u32 = b"AAGHWLK".iter().map(|&r| atoms(r).unwrap().1 as u32).sum();
        let plain_mass: f32 = b"AAGHWLK".iter().map(|&r| crate::mass::monoisotopic(r)).sum::<f32>() + crate::mass::H2O;
        assert!((aag[0].monoisotopic - plain_mass - n_atoms as f32 * N15).abs() < 2e-3);
        assert_eq!(proteins(aag[0]), "sp|P1|ONE_HUMAN");
        // only in the other protein: unlabelled
        assert!(forms(b"DDGHWLK").iter().all(|p| labels.isotope(p).is_none()));

        // shared: unlabelled for P2, labelled for P1, with each variable modification
        let shared = forms(b"LGCMNPEK");
        let labelled: Vec<&&Peptide> = shared.iter().filter(|p| labels.isotope(p) == Some(Isotope::N15)).collect();
        let plain: Vec<&&Peptide> = shared.iter().filter(|p| labels.isotope(p).is_none()).collect();
        assert!(labelled.iter().all(|p| proteins(p) == "sp|P1|ONE_HUMAN"));
        assert!(plain.iter().all(|p| proteins(p) == "sp|P2|TWO_HUMAN"));
        // unmodified, M oxidised, N deamidated: three of each (max one variable mod)
        assert_eq!((labelled.len(), plain.len()), (3, 3), "{shared:?}");
        let registry = &params.modifications.registry;
        let mass_at = |p: &Peptide, i: usize| registry.mass(p.modifications[i]);
        let oxidised = labelled.iter().find(|p| mass_at(p, 3) > 10.0).unwrap();
        // C: carbamidomethyl (its atoms light) + the cysteine's one 15N
        assert!((mass_at(oxidised, 2) - (57.021465 + N15)).abs() < 1e-3);
        // M: oxidation + one 15N
        assert!((mass_at(oxidised, 3) - (15.994915 + N15)).abs() < 1e-3);
        // a deamidated labelled N loses one of its two 15N with the NH: deamidation
        // + 1 x 15N (1.981), not + 2 x 15N (2.978); an unmodified labelled N is 2 x 15N
        let at_n: Vec<f32> = labelled.iter().map(|p| mass_at(p, 4)).collect();
        assert!(at_n.iter().any(|m| (m - (0.984016 + N15)).abs() < 1e-3), "{at_n:?}");
        assert!(at_n.iter().any(|m| (m - 2.0 * N15).abs() < 1e-3), "{at_n:?}");
        assert!(!at_n.iter().any(|m| (m - (0.984016 + 2.0 * N15)).abs() < 1e-3), "{at_n:?}");
        // mass: the unlabelled form plus every atom's shift
        let plain_ox = plain.iter().find(|p| mass_at(p, 3) > 10.0).unwrap();
        let n_atoms: u32 = b"LGCMNPEK".iter().map(|&r| atoms(r).unwrap().1 as u32).sum();
        assert!((oxidised.monoisotopic - plain_ox.monoisotopic - n_atoms as f32 * N15).abs() < 2e-3);

        // with the unlabelled form: P1 is also searched unlabelled
        let (params, peptides) = digested(true);
        let labels = params.protein_labels.clone().unwrap();
        let aag: Vec<&Peptide> = peptides.iter().filter(|p| p.sequence.as_ref() == b"AAGHWLK").collect();
        assert_eq!(aag.len(), 2);
        assert!(aag.iter().any(|p| labels.isotope(p).is_none() && proteins(p) == "sp|P1|ONE_HUMAN"));
        let shared_plain = peptides
            .iter()
            .find(|p| p.sequence.as_ref() == b"LGCMNPEK" && labels.isotope(p).is_none() && p.modifications[3].is_none())
            .unwrap();
        assert_eq!(proteins(shared_plain), "sp|P1|ONE_HUMAN;sp|P2|TWO_HUMAN");
    }

    /// SILAC on top of a uniformly labelled protein: its peptides are also
    /// searched in the SILAC heavy channel, the K carrying both labels; and
    /// the peptides are written with each residue's mass
    #[test]
    fn silac_on_a_labelled_protein() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "labelled_proteins": {"proteins": {"P1": "15N"}, "include_unlabelled": false},
            "silac": {"channels": [{"name": "heavy", "labels": {"K": 8.014199, "R": 10.008269}}]}
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let fasta = crate::fasta::Fasta::parse(">sp|P1|ONE_HUMAN labelled\nMAAGHWLK\n".into(), "rev_", false);
        let peptides = params.digest(&fasta);
        let silac = params.labels.clone().unwrap();
        let uniform = params.protein_labels.clone().unwrap();
        let registry = &params.modifications.registry;
        let forms: Vec<&Peptide> = peptides.iter().filter(|p| p.sequence.as_ref() == b"MAAGHWLK").collect();
        // 15N light and 15N + SILAC heavy (no unlabelled form)
        assert_eq!(forms.len(), 2, "{:?}", forms.iter().map(|p| p.display(registry).to_string()).collect::<Vec<_>>());
        assert!(forms.iter().all(|p| uniform.isotope(p) == Some(Isotope::N15)));
        let heavy = forms.iter().find(|p| silac.channel(p) == Some(1)).unwrap();
        let light = forms.iter().find(|p| silac.channel(p) == Some(0)).unwrap();
        assert!((heavy.monoisotopic - light.monoisotopic - 8.014199).abs() < 1e-3);
        // the heavy K: two 15N plus Lys8
        assert!((registry.mass(heavy.modifications[7]) - (2.0 * N15 + 8.014199)).abs() < 1e-3);
        // SILAC's light form of it is the 15N peptide
        assert_eq!(silac.unlabel(heavy).modifications, light.modifications);
        // written with each residue's mass: one N (0.997), two (K), three (H)
        let shown = light.display(registry).to_string();
        assert!(shown.starts_with("M[+0.997") && shown.contains("H[+2.991") && shown.contains("K[+1.994"), "{shown}");
        println!("{shown}");
    }

    #[test]
    fn proteins_are_named_by_accession_or_its_parts() {
        let labels = ProteinLabels {
            settings: Some(LabelledProteins {
                proteins: [("RPB1_YEAST".to_string(), Isotope::C13), ("sp|P08518|RPB2_YEAST".to_string(), Isotope::N15)].into(),
                include_unlabelled: true,
            }),
            ..Default::default()
        };
        assert_eq!(labels.label_of("sp|P04050|RPB1_YEAST"), Some(Isotope::C13));
        assert_eq!(labels.label_of("sp|P08518|RPB2_YEAST"), Some(Isotope::N15));
        assert_eq!(labels.label_of("sp|P08519|RPB3_YEAST"), None);
    }
}
