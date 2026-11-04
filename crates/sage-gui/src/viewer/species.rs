//! Layout of the peptide map of a multivalent species, shared by the viewer
//! and the figure export: the peptides as rows, stacked so that linked
//! peptides are neighbours where possible, each row shifted to line its link
//! residue up with its partner's, and one bridge (linker) per link. A bridge
//! between rows that are not neighbours runs along a lane right of the rows.

/// Where the rows and bridges go, in rows and residue cells
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesLayout {
    /// The peptide shown in each row, top to bottom
    pub rows: Vec<usize>,
    /// Cell of each row's first residue (by row)
    pub offsets: Vec<usize>,
    /// Cells of the widest row
    pub cells: usize,
    pub bridges: Vec<Bridge>,
    /// Lanes taken by bridges between rows that are not neighbours
    pub lanes: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bridge {
    /// Index into the links
    pub link: usize,
    /// `(row, residue)` of the end in the upper row and of the one in the lower row
    pub upper: (usize, usize),
    pub lower: (usize, usize),
    /// Rows apart: the lane it runs along
    pub lane: Option<usize>,
}

impl SpeciesLayout {
    /// The cell (column) of residue `residue` in row `row`
    pub fn column(&self, row: usize, residue: usize) -> usize {
        self.offsets[row] + residue
    }
}

/// Every ordering of `0..n`
fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn walk(rest: &mut Vec<usize>, current: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if rest.is_empty() {
            out.push(current.clone());
            return;
        }
        for i in 0..rest.len() {
            let x = rest.remove(i);
            current.push(x);
            walk(rest, current, out);
            current.pop();
            rest.insert(i, x);
        }
    }
    let mut out = Vec::new();
    walk(&mut (0..n).collect(), &mut Vec::new(), &mut out);
    out
}

/// Lay out peptides of `lengths` residues joined by `links`
/// (`(peptide, residue, peptide, residue)`, 0-based)
pub fn layout(lengths: &[usize], links: &[(usize, usize, usize, usize)]) -> SpeciesLayout {
    let n = lengths.len();
    // the row order with the fewest rows skipped by bridges; first such
    // order (peptide order when it is as good)
    let skipped = |order: &[usize]| -> usize {
        let mut row_of = vec![0; n];
        for (row, &p) in order.iter().enumerate() {
            row_of[p] = row;
        }
        links
            .iter()
            .filter(|l| l.0 < n && l.2 < n)
            .map(|&(a, _, b, _)| row_of[a].abs_diff(row_of[b]).saturating_sub(1))
            .sum()
    };
    let rows = if n <= 6 {
        permutations(n).into_iter().min_by_key(|o| skipped(o)).unwrap_or_default()
    } else {
        (0..n).collect()
    };
    let mut row_of = vec![0; n];
    for (row, &p) in rows.iter().enumerate() {
        row_of[p] = row;
    }

    // each row lines up a link residue with its partner in the row above
    // (else in any row above); the first row starts at 0
    let mut offset: Vec<i64> = vec![0; n];
    for row in 1..n {
        let p = rows[row];
        let ends = |l: &(usize, usize, usize, usize)| -> Option<(usize, usize, usize)> {
            // (partner row, partner residue, residue of p)
            match *l {
                (a, ra, b, rb) if a == p && b < n && row_of[b] < row => Some((row_of[b], rb, ra)),
                (a, ra, b, rb) if b == p && a < n && row_of[a] < row => Some((row_of[a], ra, rb)),
                _ => None,
            }
        };
        let partner = links
            .iter()
            .filter_map(ends)
            .min_by_key(|&(r, _, _)| row - r);
        if let Some((r, partner_residue, residue)) = partner {
            offset[row] = offset[r] + partner_residue as i64 - residue as i64;
        }
    }
    let least = offset.iter().copied().min().unwrap_or(0);
    let offsets: Vec<usize> = offset.iter().map(|o| (o - least) as usize).collect();
    let cells = (0..n).map(|row| offsets[row] + lengths[rows[row]]).max().unwrap_or(0);

    let mut lanes = 0;
    let bridges = links
        .iter()
        .enumerate()
        .filter(|(_, l)| l.0 < n && l.2 < n)
        .map(|(e, &(a, ra, b, rb))| {
            let (ends_a, ends_b) = ((row_of[a], ra), (row_of[b], rb));
            let (upper, lower) = if ends_a.0 <= ends_b.0 { (ends_a, ends_b) } else { (ends_b, ends_a) };
            let lane = (lower.0 - upper.0 > 1).then(|| {
                lanes += 1;
                lanes - 1
            });
            Bridge { link: e, upper, lower, lane }
        })
        .collect();
    SpeciesLayout { rows, offsets, cells, bridges, lanes }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn copies_linked_twice_line_up() {
        // a2a: two copies of one peptide, K18–K18 and K22–K22
        let l = layout(&[37, 37], &[(0, 17, 1, 17), (0, 21, 1, 21)]);
        assert_eq!(l.rows, [0, 1]);
        assert_eq!(l.offsets, [0, 0]);
        assert!(l.bridges.iter().all(|b| b.lane.is_none() && b.upper.1 == b.lower.1));
    }

    #[test]
    fn a_chain_puts_the_middle_peptide_between() {
        // a1b1a: peptides 1 and 2 both linked to peptide 3
        let l = layout(&[16, 16, 9], &[(0, 5, 2, 0), (1, 5, 2, 4)]);
        assert_eq!(l.rows[1], 2, "peptide 3 in the middle: {:?}", l.rows);
        assert_eq!(l.lanes, 0);
        // every bridge joins residues in the same column
        for b in &l.bridges {
            let (cu, cl) = (l.column(b.upper.0, b.upper.1), l.column(b.lower.0, b.lower.1));
            assert_eq!(cu, cl, "{b:?} in {l:?}");
        }
    }

    #[test]
    fn a_branch_routes_the_bridge_that_skips_a_row() {
        // a1b,b1c,b1d: b linked to three peptides; one link has to skip a row
        let l = layout(&[10, 12, 8, 9], &[(0, 2, 1, 3), (1, 5, 2, 1), (1, 8, 3, 2)]);
        assert_eq!(l.lanes, 1);
        assert_eq!(l.bridges.iter().filter(|b| b.lane.is_some()).count(), 1);
        assert!(l.offsets.contains(&0));
        assert!(l.cells >= 12);
    }

    #[test]
    fn two_links_between_two_peptides() {
        // a2b: only one of the two links can be lined up
        let l = layout(&[12, 14], &[(0, 3, 1, 9), (0, 7, 1, 2)]);
        assert_eq!(l.lanes, 0);
        let aligned = l
            .bridges
            .iter()
            .filter(|b| l.column(b.upper.0, b.upper.1) == l.column(b.lower.0, b.lower.1))
            .count();
        assert_eq!(aligned, 1);
    }
}
