//! Tree molecules, their canonical form, and one-atom edits.
//!
//! Trees rather than general graphs, deliberately.
//!
//! Canonical labelling of a tree is exact and cheap (AHU), so the measurement
//! is never confounded by a canonicalisation that is itself approximate. G2
//! asks whether *signature* space has locality; a heuristic canonical form
//! would leave it ambiguous whether a null result came from the signature or
//! from the labelling.
//!
//! Elements are indices into an invented radius table. There is no element
//! data here and no correspondence to anything real (§5, G1/G3) — the atoms
//! differ only in how much space they take up, which is all a shape-based
//! chemistry needs them to differ in.

#![allow(
    clippy::indexing_slicing,
    reason = "adjacency and distance indices are bounded by len(), checked at construction"
)]

use crate::rng::Stream;

/// FNV-1a's 64-bit offset basis — the starting value of the hash.
///
/// Both this and [`FNV_PRIME`] are the published 64-bit FNV parameters
/// (Fowler–Noll–Vo, 1991). The prime is chosen so that repeated multiplication
/// spreads byte values across the whole register; the basis is arbitrary but
/// fixed, and exists so that hashing the empty string is not zero.
const FNV_OFFSET_BASIS: u64 = 0xCBF2_9CE4_8422_2325;

/// FNV-1a's 64-bit prime. See [`FNV_OFFSET_BASIS`].
const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;

/// How many element types the harness draws from.
pub const N_ELEMENTS: usize = 4;

/// Invented atomic radii, in arbitrary spans.
///
/// Spread widely enough that element identity has a visible effect on shape —
/// if every atom were the same size, a retype would be invisible to the
/// signature and the G2 question would be answered by construction rather
/// than by measurement.
pub const ELEMENT_RADII: [f64; N_ELEMENTS] = [0.55, 0.80, 1.05, 1.40];

/// A molecule: typed atoms joined into a tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Molecule {
    /// Element index per atom.
    pub elements: Vec<u8>,
    edges: Vec<(usize, usize)>,
}

impl Molecule {
    /// Assemble from atom types and bonds. Bonds are stored as given; the
    /// caller is responsible for them forming a tree.
    #[must_use]
    pub const fn from_parts(elements: Vec<u8>, edges: Vec<(usize, usize)>) -> Self {
        Self { elements, edges }
    }

    /// Number of atoms.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.elements.len()
    }

    /// Whether the molecule has no atoms.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// The bonds, in construction order.
    #[must_use]
    #[allow(
        clippy::missing_const_for_fn,
        reason = "clippy suggests `const` here and the suggestion does not compile: coercing \
                  Vec to a slice needs Deref, which is not const-stable. A nursery false positive"
    )]
    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }

    /// A uniformly-grown random tree on `n` atoms.
    ///
    /// Each new atom attaches to a uniformly chosen existing one, which
    /// produces a mix of paths and bushes rather than the near-stars a
    /// preferential-attachment rule would give.
    #[must_use]
    pub fn random_tree(rng: &mut Stream, n: usize) -> Self {
        let mut elements = Vec::with_capacity(n);
        let mut edges = Vec::with_capacity(n.saturating_sub(1));
        for i in 0..n {
            elements.push(Self::draw_element(rng));
            if i > 0 {
                let parent = rng.index(i);
                edges.push((parent, i));
            }
        }
        Self { elements, edges }
    }

    fn draw_element(rng: &mut Stream) -> u8 {
        u8::try_from(rng.index(N_ELEMENTS)).unwrap_or(0)
    }

    /// The same molecule with atom `i` renamed to `perm[i]`.
    ///
    /// A test fixture in spirit, but public because relabelling invariance is
    /// the property the canonical form exists to have, and a test that cannot
    /// construct a relabelling cannot check it.
    #[must_use]
    pub fn relabelled(&self, perm: &[usize]) -> Self {
        let mut elements = vec![0u8; self.len()];
        for (old, &new) in perm.iter().enumerate() {
            elements[new] = self.elements[old];
        }
        let edges = self
            .edges
            .iter()
            .map(|&(a, b)| (perm[a], perm[b]))
            .collect();
        Self { elements, edges }
    }

    /// Adjacency lists, each sorted ascending so traversal order is fixed.
    #[must_use]
    pub fn adjacency(&self) -> Vec<Vec<usize>> {
        let mut adj = vec![Vec::new(); self.len()];
        for &(a, b) in &self.edges {
            adj[a].push(b);
            adj[b].push(a);
        }
        for row in &mut adj {
            row.sort_unstable();
        }
        adj
    }

    /// All-pairs shortest path lengths in bonds, by breadth-first search from
    /// every atom. `f64::INFINITY` for an unreachable pair, which cannot
    /// happen for a tree but is not silently papered over.
    #[must_use]
    pub fn graph_distances(&self) -> Vec<Vec<f64>> {
        let adj = self.adjacency();
        let n = self.len();
        let mut out = Vec::with_capacity(n);
        for src in 0..n {
            let mut d = vec![f64::INFINITY; n];
            d[src] = 0.0;
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(src);
            while let Some(v) = queue.pop_front() {
                for &w in &adj[v] {
                    if d[w].is_infinite() {
                        d[w] = d[v] + 1.0;
                        queue.push_back(w);
                    }
                }
            }
            out.push(d);
        }
        out
    }

    /// The tree's centroid atoms — the one or two atoms minimising the largest
    /// subtree hanging off them.
    ///
    /// Rooting at the centroid is what makes the canonical form canonical: a
    /// tree has at most two centroids and they are determined by its shape
    /// alone, so no atom numbering can influence the choice.
    fn centroids(&self) -> Vec<usize> {
        let n = self.len();
        if n <= 2 {
            return (0..n).collect();
        }
        let adj = self.adjacency();
        let mut degree: Vec<usize> = adj.iter().map(Vec::len).collect();
        let mut removed = vec![false; n];
        let mut remaining = n;
        let mut layer: Vec<usize> = (0..n).filter(|&v| degree[v] <= 1).collect();

        while remaining > 2 {
            let mut next = Vec::new();
            for &v in &layer {
                removed[v] = true;
                remaining -= 1;
                for &w in &adj[v] {
                    if !removed[w] {
                        degree[w] -= 1;
                        if degree[w] == 1 {
                            next.push(w);
                        }
                    }
                }
            }
            next.sort_unstable();
            layer = next;
        }
        (0..n).filter(|&v| !removed[v]).collect()
    }

    /// AHU canonical form: each subtree becomes `<element>(<children...>)`
    /// with the children's own forms sorted, so sibling order — the only place
    /// atom numbering could enter — is normalised away.
    ///
    /// The whole-tree form is the lexicographically smaller of the rootings at
    /// each centroid, which settles the two-centroid case without a tie-break
    /// that could depend on numbering.
    #[must_use]
    pub fn canonical_form(&self) -> String {
        let adj = self.adjacency();
        self.centroids()
            .into_iter()
            .map(|root| Self::subtree_form(&adj, &self.elements, root, usize::MAX))
            .min()
            .unwrap_or_default()
    }

    fn subtree_form(adj: &[Vec<usize>], elements: &[u8], v: usize, parent: usize) -> String {
        let mut child_forms: Vec<String> = adj[v]
            .iter()
            .filter(|&&w| w != parent)
            .map(|&w| Self::subtree_form(adj, elements, w, v))
            .collect();
        child_forms.sort_unstable();
        format!("{}({})", elements[v], child_forms.concat())
    }

    /// Atoms in canonical order: a depth-first walk from the canonical root,
    /// visiting children in order of their own canonical forms.
    ///
    /// This is the ordering a label-anchored frame would be built from, and it
    /// is the thing G2 puts on trial. It is canonical — two numberings of the
    /// same molecule give the same sequence — but it is not *continuous*:
    /// adding one atom can change which centroid wins and reorder everything
    /// below it. `D_frame` measures what that costs.
    #[must_use]
    pub fn canonical_order(&self) -> Vec<usize> {
        let adj = self.adjacency();
        let Some(root) = self
            .centroids()
            .into_iter()
            .min_by_key(|&r| Self::subtree_form(&adj, &self.elements, r, usize::MAX))
        else {
            return Vec::new();
        };
        let mut out = Vec::with_capacity(self.len());
        self.walk(&adj, root, usize::MAX, &mut out);
        out
    }

    fn walk(&self, adj: &[Vec<usize>], v: usize, parent: usize, out: &mut Vec<usize>) {
        out.push(v);
        let mut children: Vec<(String, usize)> = adj[v]
            .iter()
            .filter(|&&w| w != parent)
            .map(|&w| (Self::subtree_form(adj, &self.elements, w, v), w))
            .collect();
        // Sorted by canonical form, with the atom index as the tie-break.
        // Ties are isomorphic subtrees, so the choice between them cannot
        // change the geometry — but pinning it keeps the walk deterministic.
        children.sort();
        for (_, w) in children {
            self.walk(adj, w, v, out);
        }
    }

    /// A 64-bit digest of the canonical form.
    ///
    /// This is the **negative control**: it identifies the species exactly, so
    /// it carries all the information a descriptor could want, yet an
    /// avalanching hash destroys every trace of similarity. A descriptor built
    /// on it must show no locality whatsoever. If the harness reports locality
    /// here, the harness is broken and no other number it produced means
    /// anything.
    #[must_use]
    pub fn form_hash(&self) -> u64 {
        // FNV-1a over the canonical form's bytes, then SplitMix64's finaliser
        // to avalanche. FNV alone leaves low-order structure that would let
        // similar forms land near each other — which for a negative control
        // would be a false result rather than a mild flaw.
        let mut h: u64 = FNV_OFFSET_BASIS;
        for b in self.canonical_form().bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(FNV_PRIME);
        }
        crate::rng::mix64(h)
    }

    /// Change one atom's element to a different one.
    #[must_use]
    pub fn retype_one(&self, rng: &mut Stream) -> Self {
        let mut child = self.clone();
        if self.is_empty() {
            return child;
        }
        let atom = rng.index(self.len());
        // Draw from the other N-1 elements so the edit always bites; a retype
        // to the same element is a no-op that would dilute the mutant sample
        // with exact copies of the parent.
        let shift = rng.index(N_ELEMENTS - 1) + 1;
        let e = (usize::from(child.elements[atom]) + shift) % N_ELEMENTS;
        child.elements[atom] = u8::try_from(e).unwrap_or(0);
        child
    }

    /// Attach one new atom to a uniformly chosen existing one.
    #[must_use]
    pub fn add_leaf(&self, rng: &mut Stream) -> Self {
        let mut child = self.clone();
        if self.is_empty() {
            child.elements.push(Self::draw_element(rng));
            return child;
        }
        let anchor = rng.index(self.len());
        let new = child.len();
        child.elements.push(Self::draw_element(rng));
        child.edges.push((anchor, new));
        child
    }

    /// Remove one uniformly chosen leaf, renumbering the atoms above it.
    ///
    /// Leaves only, because removing an interior atom would disconnect the
    /// tree — and reconnecting the fragments would be a multi-atom edit
    /// wearing a one-atom label.
    #[must_use]
    pub fn remove_leaf(&self, rng: &mut Stream) -> Self {
        let adj = self.adjacency();
        let leaves: Vec<usize> = (0..self.len()).filter(|&v| adj[v].len() <= 1).collect();
        if leaves.is_empty() || self.len() <= 1 {
            return self.clone();
        }
        let victim = leaves[rng.index(leaves.len())];

        let renumber = |v: usize| if v > victim { v - 1 } else { v };
        let elements = (0..self.len())
            .filter(|&v| v != victim)
            .map(|v| self.elements[v])
            .collect();
        let edges = self
            .edges
            .iter()
            .filter(|&&(a, b)| a != victim && b != victim)
            .map(|&(a, b)| (renumber(a), renumber(b)))
            .collect();
        Self { elements, edges }
    }

    /// One uniformly chosen single-atom edit: retype, grow, or shrink.
    #[must_use]
    pub fn mutate(&self, rng: &mut Stream) -> Self {
        match rng.below(3) {
            0 => self.retype_one(rng),
            1 => self.add_leaf(rng),
            _ => self.remove_leaf(rng),
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::needless_range_loop,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds; casts are bounded by \
              N_ELEMENTS; the float comparisons are against exact 0.0 on a diagonal that is set, \
              not computed; the range loops index two parallel structures at once"
)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    /// A path `0-1-2-...-(n-1)`, every atom the same element.
    fn path(n: usize) -> Molecule {
        Molecule::from_parts(vec![0; n], (1..n).map(|i| (i - 1, i)).collect())
    }

    /// A star: atom 0 joined to every other, every atom the same element.
    fn star(n: usize) -> Molecule {
        Molecule::from_parts(vec![0; n], (1..n).map(|i| (0, i)).collect())
    }

    #[test]
    fn random_molecules_are_trees() {
        let mut rng = Stream::new(1);
        for n in 2..30 {
            let m = Molecule::random_tree(&mut rng, n);
            assert_eq!(m.len(), n);
            assert_eq!(m.edges().len(), n - 1, "n={n}: not n-1 edges");
            // Connected: BFS from atom 0 reaches everything. With n-1 edges
            // that is exactly the definition of a tree.
            let d = m.graph_distances();
            for j in 0..n {
                assert!(d[0][j].is_finite(), "n={n}: atom {j} unreachable");
            }
        }
    }

    #[test]
    fn graph_distances_are_a_metric() {
        let mut rng = Stream::new(2);
        let m = Molecule::random_tree(&mut rng, 20);
        let d = m.graph_distances();
        for i in 0..m.len() {
            assert_eq!(d[i][i], 0.0, "diagonal at {i}");
            for j in 0..m.len() {
                assert_eq!(d[i][j], d[j][i], "asymmetric at {i},{j}");
                assert!(d[i][j] >= 0.0);
                for k in 0..m.len() {
                    assert!(
                        d[i][j] <= d[i][k] + d[k][j] + 1e-12,
                        "triangle violated at {i},{j} via {k}"
                    );
                }
            }
        }
    }

    /// The property the whole canonicalisation exists for. If atom indices
    /// leak into the canonical form then every downstream identity — species
    /// interning, the negative control, the label-anchored frame — is keyed on
    /// an accident of construction order.
    #[test]
    fn the_canonical_form_is_invariant_under_relabelling() {
        let mut rng = Stream::new(3);
        for trial in 0..200 {
            let m = Molecule::random_tree(&mut rng, 4 + (trial % 16));
            let mut perm: Vec<usize> = (0..m.len()).collect();
            // Fisher-Yates, drawn from the same stream so the trial is
            // reproducible from the seed alone.
            for i in (1..perm.len()).rev() {
                let j = rng.index(i + 1);
                perm.swap(i, j);
            }
            let relabelled = m.relabelled(&perm);
            assert_eq!(
                m.canonical_form(),
                relabelled.canonical_form(),
                "trial {trial}: canonical form depends on atom numbering"
            );
        }
    }

    #[test]
    fn the_canonical_form_separates_different_shapes() {
        assert_ne!(path(4).canonical_form(), star(4).canonical_form());
        assert_ne!(path(5).canonical_form(), star(5).canonical_form());
        // Same shape, same size — must agree.
        assert_eq!(path(5).canonical_form(), path(5).canonical_form());
    }

    #[test]
    fn the_canonical_form_notices_a_retyped_atom() {
        let mut rng = Stream::new(4);
        let m = Molecule::random_tree(&mut rng, 12);
        for atom in 0..m.len() {
            for e in 0..N_ELEMENTS as u8 {
                if e == m.elements[atom] {
                    continue;
                }
                let mut other = m.clone();
                other.elements[atom] = e;
                assert_ne!(
                    m.canonical_form(),
                    other.canonical_form(),
                    "retyping atom {atom} to element {e} left the form unchanged"
                );
            }
        }
    }

    #[test]
    fn the_canonical_order_is_a_permutation() {
        let mut rng = Stream::new(5);
        for n in 2..25 {
            let m = Molecule::random_tree(&mut rng, n);
            let order = m.canonical_order();
            assert_eq!(order.len(), n);
            let mut seen = vec![false; n];
            for &a in &order {
                assert!(!seen[a], "atom {a} appears twice in the canonical order");
                seen[a] = true;
            }
        }
    }

    /// The canonical *order* must be relabelling-covariant: it has to select
    /// the same physical atoms in the same sequence, whatever they are called.
    /// Symmetric subtrees are genuinely interchangeable, so the test compares
    /// the element sequence rather than demanding identical indices.
    #[test]
    fn the_canonical_order_selects_the_same_atoms_under_relabelling() {
        let mut rng = Stream::new(6);
        for trial in 0..100 {
            let m = Molecule::random_tree(&mut rng, 5 + (trial % 14));
            let mut perm: Vec<usize> = (0..m.len()).collect();
            for i in (1..perm.len()).rev() {
                let j = rng.index(i + 1);
                perm.swap(i, j);
            }
            let relabelled = m.relabelled(&perm);

            let a: Vec<u8> = m.canonical_order().iter().map(|&i| m.elements[i]).collect();
            let b: Vec<u8> = relabelled
                .canonical_order()
                .iter()
                .map(|&i| relabelled.elements[i])
                .collect();
            assert_eq!(
                a, b,
                "trial {trial}: canonical order is numbering-dependent"
            );
        }
    }

    #[test]
    fn retyping_changes_exactly_one_element_and_no_bonds() {
        let mut rng = Stream::new(7);
        let m = Molecule::random_tree(&mut rng, 15);
        let child = m.retype_one(&mut rng);
        assert_eq!(child.len(), m.len());
        assert_eq!(child.edges(), m.edges());
        let differing = (0..m.len())
            .filter(|&i| m.elements[i] != child.elements[i])
            .count();
        assert_eq!(differing, 1, "a retype changed {differing} atoms");
    }

    #[test]
    fn adding_a_leaf_grows_the_molecule_by_exactly_one_atom() {
        let mut rng = Stream::new(8);
        let m = Molecule::random_tree(&mut rng, 15);
        let child = m.add_leaf(&mut rng);
        assert_eq!(child.len(), m.len() + 1);
        assert_eq!(child.edges().len(), m.edges().len() + 1);
        // The parent's atoms and bonds survive untouched: an edit is local.
        assert_eq!(child.elements[..m.len()], m.elements[..]);
        assert_eq!(child.edges()[..m.edges().len()], m.edges()[..]);
    }

    #[test]
    fn removing_a_leaf_shrinks_the_molecule_by_exactly_one_atom() {
        let mut rng = Stream::new(9);
        let m = Molecule::random_tree(&mut rng, 15);
        let child = m.remove_leaf(&mut rng);
        assert_eq!(child.len(), m.len() - 1);
        assert_eq!(child.edges().len(), m.edges().len() - 1);
        let d = child.graph_distances();
        for j in 0..child.len() {
            assert!(d[0][j].is_finite(), "removing a leaf disconnected the tree");
        }
    }

    #[test]
    fn every_mutation_yields_a_tree_of_nearly_the_same_size() {
        let mut rng = Stream::new(10);
        for _ in 0..300 {
            let m = Molecule::random_tree(&mut rng, 12);
            let child = m.mutate(&mut rng);
            assert_eq!(child.edges().len(), child.len() - 1, "mutant is not a tree");
            let delta = child.len().abs_diff(m.len());
            assert!(delta <= 1, "a one-atom edit changed the size by {delta}");
            let d = child.graph_distances();
            for j in 0..child.len() {
                assert!(d[0][j].is_finite(), "mutant is disconnected");
            }
        }
    }

    /// Two molecules that differ must not collide in the negative control's
    /// hash, or the control is measuring hash collisions rather than the
    /// absence of locality.
    #[test]
    fn distinct_molecules_get_distinct_form_hashes() {
        let mut rng = Stream::new(11);
        let mut forms = std::collections::BTreeMap::new();
        for _ in 0..2000 {
            let m = Molecule::random_tree(&mut rng, 10);
            let f = m.canonical_form();
            let h = m.form_hash();
            if let Some(prev) = forms.insert(h, f.clone()) {
                assert_eq!(prev, f, "hash collision between distinct canonical forms");
            }
        }
    }
}
