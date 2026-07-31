//! The shape signature (spec §8.2) — the heart of the engine.
//!
//! A molecule reduces to two numbers per sample direction: how far its surface
//! extends that way ([`Signature::extents`]), and what surface character faces
//! that way ([`Signature::characters`]). Everything downstream — binding,
//! catalysis, membranes, permeability — is a comparison between two of these.
//!
//! # Every comparison is group-minimised, and that is enforced
//!
//! [`Signature::group_distance`] is the **only** distance this module offers,
//! and `xtask`'s `check_no_ungrouped_signature_comparison` is what keeps it
//! that way. The reason is measured rather than stylistic. G2 scored five ways
//! of comparing signatures; ranked by how well a one-atom edit stays near its
//! parent, size-matched over the sizes [`Mol12`](crate::graph::Mol12) admits:
//!
//! | descriptor | size-matched |
//! |---|---|
//! | `D_group` — minimise over the 60 rotations | 0.9230–0.9415 |
//! | `D_raw` — no alignment at all | 0.9075–0.9355 |
//! | `D_sorted` — orientation discarded by sorting | 0.8840–0.9075 |
//! | `D_frame` — deleted at Task 8 for being worse than nothing | 0.8585–0.8885 |
//! | **`D_lexmin` — §8.2's canonicalisation, as written** | **0.7510–0.7745** |
//!
//! §8.2 specifies that a signature is *stored* as the lexicographically
//! smallest of its 60 rotations, and that form is the **worst of the five** —
//! below the frame `layout` deleted for scoring below no-alignment at all.
//!
//! **It is nonetheless safe, and exactly conditionally so.** Lex-min is itself
//! a group element, so `D_group(lexmin(A), lexmin(B))` is bit-identical to
//! `D_group(A, B)`: storing the canonical form costs nothing *provided every
//! consumer compares with the group-minimised distance*, and is severe for any
//! consumer that subtracts two stored signatures directly. §8.2 names three
//! consumers. Species identity and the fold cache read the stored form for
//! **equality**, which is what [`Signature::canonicalise`] is for and is fine.
//! The third is §15.2's shape-space novelty histogram, which takes
//! **distances** and is read against a neutral shadow — where a descriptor
//! with poor locality inflates apparent novelty as a pure representation
//! artefact, in the one metric whose job is to say which differences are not
//! significant.
//!
//! So §8.2 needs no amendment. The storage form is absorbed exactly by the
//! comparison rule, and the comparison rule is the thing with teeth.
//!
//! ## What the enforcement decides, and what is outside it
//!
//! A group-minimised comparison cannot be written without a [`Geodesic`] — it
//! is the only source of `rotation_perms`. A plain Euclidean one needs no such
//! argument. The check therefore reads this file's public surface and requires
//! that **every `pub fn` mentioning `Signature` in two or more of its inputs
//! also takes a `Geodesic`**. `Signature::distance(&self, &Self) -> f64` has
//! no spelling that passes.
//!
//! Stated as a bound rather than as "by construction", because three guards in
//! this project have shipped documenting a requirement without enforcing it:
//!
//! - a function that takes a `Geodesic` and then ignores it is **not** caught;
//!   `syn` reads signatures, not bodies;
//! - a consumer crate can fold over [`Signature::extents`] and subtract. Those
//!   accessors exist because §8.3's kernel needs the components, so this
//!   cannot be closed by privacy without breaking Task 10;
//! - a comparison against something that is not a `Signature` — a bare
//!   `[Span; D]` — is outside the framing entirely.
//!
//! What it does cover is the place a `Signature::distance` would actually be
//! added, which is here.
//!
//! ## Why there is no `sorted_extents`
//!
//! The plan's interface list carries one and it is deliberately not
//! implemented. Both halves of the argument are worth keeping, because the
//! first is a real point in its favour:
//!
//! **For.** As a *filter*, sorting is sound. Two signatures in the same proper
//! orbit have equal sorted extents, so it is a genuine necessary condition,
//! and being invariant under a *larger* group than the 60 only makes it more
//! permissive — which is the safe direction for a filter.
//!
//! **Against, and it decides.** §8.3 does not ask for a filter, it asks for "a
//! rotation-invariant **upper bound on the score**", and adds that "a filter
//! that cannot state the property it guarantees is not conservative, it is
//! merely untested". A sorted extent vector is not that bound; constructing it
//! is Task 10's job and it needs the binding constants to do it. Meanwhile
//! `D_sorted` is invariant under the **full** icosahedral group including the
//! improper elements — a configuration and its point inversion have identical
//! sorted extents — so any distance built on it makes enantiomers
//! indistinguishable and §22.8 impossible by construction, the same shape as
//! the dropped `ANTI`. It currently ranks third of five, so a future
//! measurement where it edges ahead is a live trap.
//!
//! Shipping it now would publish a rotation-invariant-looking handle that is
//! neither the bound §8.3 needs nor safe as a distance, guarded only by a doc
//! comment saying "do not compare with this". Routed to Task 10.

use crate::geodesic::{Geodesic, Rotation};
use crate::layout::Embedding;
use borbax_units::{Span, canonical_cmp};
use core::cmp::Ordering;

/// Width of the surface shell that contributes to surface character: atoms
/// within this of the supporting plane are "facing" the direction.
///
/// **Chosen, not derived**, and named as such for the reason `bonds.rs`
/// records about the strain exponent — a constant whose justification is "it
/// looked right" invites a later reader to treat it as load-bearing physics.
/// It sets how much of the molecule the character channel can see: at zero
/// only the single extremal atom contributes, and above the molecular diameter
/// every atom does and the channel becomes a composition average with no
/// geometry in it. Drawn element radii span roughly `[0.3, 2.04]`, so 0.75
/// admits the extremal atom's immediate neighbours and little else.
///
/// It has no measurement behind it and no sweep yet. Task 20 owns that; the
/// test `the_character_averages_only_the_atoms_facing_that_way` at least
/// asserts the shell is excluding *something*, so a value that quietly admits
/// everything fails rather than passing as a plain average.
const SHELL: Span = Span(0.75);

/// A molecule's shape, as §8.3's binding kernel sees it.
///
/// Not `Copy`. At `D = 162` this is 2 592 bytes and an implicit copy per use
/// is not something to hand out silently; `permuted` and `canonicalise` return
/// owned values because they genuinely build new ones.
#[derive(Debug, Clone, PartialEq)]
pub struct Signature<const D: usize> {
    /// Surface extent per direction — the support function of the union of
    /// atom spheres.
    r: [Span; D],
    /// Surface character per direction, in `[-1, 1]`.
    a: [f64; D],
}

impl<const D: usize> Signature<D> {
    const fn zeroed() -> Self {
        Self {
            r: [Span::ZERO; D],
            a: [0.0; D],
        }
    }

    /// Surface extent per direction.
    ///
    /// [`Span`], not `f64`: §8.3 compares these against `ideal_gap`, which is
    /// also a length, and an untyped extent is where G1's mix-up would enter.
    #[must_use]
    pub const fn extents(&self) -> &[Span; D] {
        &self.r
    }

    /// Surface character per direction, in `[-1, 1]`.
    ///
    /// Dimensionless — §7.1's `affinity` is an exposed *fraction*, so there is
    /// no unit for it to carry.
    #[must_use]
    pub const fn characters(&self) -> &[f64; D] {
        &self.a
    }

    /// This signature with rotation `r` applied, given as its permutation of
    /// directions.
    ///
    /// **The convention is `out[perm[i]] = self[i]`, and it is pinned by
    /// measurement rather than by this sentence.** Rotating an embedding
    /// bodily by the matrix for `r` and recomputing the signature from scratch
    /// gives exactly this — `rotating_the_embedding_permutes_the_signature`.
    /// The inverse convention is equally self-consistent and would leave every
    /// invariance test in this file green, because a group is closed under
    /// inverse; only a test that ties the table to real geometry tells them
    /// apart. That is the defect `contact_perms`' docstring shipped with at
    /// Task 7, where the prose named the wrong pose while the formula beside
    /// it was right.
    ///
    /// Passing [`Geodesic::anti`] here gives the **mirror image**, since
    /// `anti` is an involution: `out[anti[i]] = self[i]` is `out[j] =
    /// self[anti[j]]`.
    #[must_use]
    pub fn permuted(&self, perm: &[u8; D]) -> Self {
        let mut out = Self::zeroed();
        for ((&to, extent), character) in perm.iter().zip(self.r.iter()).zip(self.a.iter()) {
            let j = usize::from(to);
            if let Some(slot) = out.r.get_mut(j) {
                *slot = *extent;
            }
            if let Some(slot) = out.a.get_mut(j) {
                *slot = *character;
            }
        }
        out
    }

    /// The lexicographically smallest of this signature's 60 rotations (§8.2).
    ///
    /// This is the **equality** form: species identity and the fold cache need
    /// two embeddings of one molecule to intern as one species, and this is
    /// what makes that true regardless of how the embedding happened to land.
    /// It is *not* a distance and nothing may subtract two of these — see the
    /// module docs for the measurement behind that.
    ///
    /// **The search is 60 proper rotations, never 120.** The improper coset is
    /// what would identify a shape with its mirror image, and §22.8 makes
    /// handedness a result the simulation may produce rather than one excluded
    /// by construction. `a_chiral_signature_and_its_mirror_are_different_species`
    /// asserts both halves: that the 60-element search keeps them apart, and
    /// that extending it to 120 collapses them.
    ///
    /// The tie-break runs over the character channel as well as the extents.
    /// **It is necessary in principle and has never once decided anything, and
    /// both halves of that are worth stating.**
    ///
    /// In principle: with an extents-only order, "first among equals in
    /// `Rotation::all()` order" is *not* a function of the orbit. `S` and `S∘g`
    /// present the same candidate multiset indexed differently, so if the
    /// minimal set held two candidates with different characters, the two
    /// enumerations would return different answers and `canonicalise` would
    /// stop being well-defined. A total order removes that.
    ///
    /// In practice: measured over 600 molecules at n = 1..12, all 1 770 pose
    /// pairs each, extent ties are **common** — 91 500 of them, concentrated at
    /// n <= 2 where the shape has a continuous symmetry the discrete group
    /// cannot break — and in **0** of them do the characters differ. So the
    /// character comparison has never been reached as a tie-break on any input
    /// this crate can produce.
    ///
    /// An earlier version of this paragraph called it "load-bearing", which the
    /// measurement refutes. It is retained because 91 500 measured ties are
    /// evidence about a corpus and not a proof that a character-differing tie
    /// is unreachable, and because the cost of keeping it is one comparison
    /// that never runs. Do not delete it on the strength of the zero.
    #[must_use]
    pub fn canonicalise(&self, g: &Geodesic<D>) -> Self {
        // `Rotation::all()` is never empty, but the type does not say so, so
        // fold rather than seed-then-loop. No `unwrap`: the workspace denies it.
        let mut best: Option<Self> = None;
        for rot in Rotation::all() {
            let cand = self.permuted(g.rotation_perms(rot));
            if best
                .as_ref()
                .is_none_or(|held| cand.lex_cmp(held) == Ordering::Less)
            {
                best = Some(cand);
            }
        }
        best.unwrap_or_else(|| unreachable!("Rotation::all() yields 60 elements"))
    }

    /// Total order on signatures, canonicalising every float first.
    ///
    /// `Span::canonical_cmp` and `borbax_units::canonical_cmp` rather than
    /// `f64::total_cmp`, which orders on the sign bit — and a runtime NaN's
    /// sign bit is clear on aarch64 and set on x86-64, so the plain spelling
    /// reorders the whole search between the macOS leg of the matrix and the
    /// other two (§13.4).
    fn lex_cmp(&self, other: &Self) -> Ordering {
        for (mine, theirs) in self.r.iter().zip(other.r.iter()) {
            match mine.canonical_cmp(theirs) {
                Ordering::Equal => {}
                ord => return ord,
            }
        }
        for (mine, theirs) in self.a.iter().zip(other.a.iter()) {
            match canonical_cmp(*mine, *theirs) {
                Ordering::Equal => {}
                ord => return ord,
            }
        }
        Ordering::Equal
    }

    /// Distance between two shapes, minimised over the 60 rotations.
    ///
    /// **The only comparison this module offers**, for the reason in the
    /// module docs. It is the quantity §8.3's `affinity` maximises over and
    /// the only one §8.2's storage decision is safe under.
    ///
    /// Returns `f64` and **not** [`Span`], deliberately: it mixes an extent
    /// channel that is a length with a character channel that is a
    /// dimensionless ratio, so the sum is not a length and typing it as one
    /// would launder exactly the confusion [`Signature::extents`] carries
    /// `Span` to prevent. It is a descriptor distance; nothing may compare it
    /// against a physical length.
    ///
    /// Accumulated in direction-index order and never through `.sum()` over an
    /// unordered iterator, so the summation order is fixed on every platform
    /// (§13.1).
    #[must_use]
    pub fn group_distance(&self, other: &Self, g: &Geodesic<D>) -> f64 {
        let mut best = f64::INFINITY;
        for rot in Rotation::all() {
            let posed = other.permuted(g.rotation_perms(rot));
            let mut acc = 0.0f64;
            for ((mine, theirs), (my_a, their_a)) in self
                .r
                .iter()
                .zip(posed.r.iter())
                .zip(self.a.iter().zip(posed.a.iter()))
            {
                let dr = (*mine - *theirs).get();
                let da = *my_a - *their_a;
                acc += dr * dr + da * da;
            }
            // The comparison form, not `f64::min`: §13.1 bans the method
            // because it is non-deterministic for ±0.0. `acc` is a sum of
            // squares and cannot be NaN here.
            if acc < best {
                best = acc;
            }
        }
        best.sqrt()
    }
}

/// Sample an embedded molecule along the geodesic directions (§8.2).
///
/// `r` is the support function of the union of atom spheres: how far the
/// surface reaches each way. `a` is the surface character averaged over the
/// atoms that actually face that way, weighted linearly by how close they sit
/// to the supporting plane.
///
/// **Takes neither a `Mol12` nor a `Universe`, and that is the point** (Task
/// 9's routed requirement 3). The planned `signature(&Mol12, &Embedding, ..)`
/// read `m.elem[i]` at canonical position `i` while `e` came from a
/// `CanonMol`, with nothing typing the agreement — any `Mol12` compiled and
/// returned a plausible signature, which is the defect `CanonMol` closed one
/// task earlier. [`Embedding`] now carries the per-atom radius and character
/// it was already fetching, so the mismatched call has no spelling. It also
/// stops the projection being computed twice per (direction, atom), which the
/// planned two-pass form did — 3 888 duplicated dot products at `D = 162`.
///
/// The character weighting is polynomial rather than an exponential falloff on
/// purpose: `exp` is not portable between platform libm implementations
/// (§13.1), and a linear shell behaves the same qualitatively with no
/// reproducibility risk.
#[must_use]
pub fn signature<const D: usize>(e: &Embedding, g: &Geodesic<D>) -> Signature<D> {
    let mut out = Signature::<D>::zeroed();
    let n = e.len();
    if n == 0 {
        return out;
    }
    let characters = e.characters();

    for (d, dir) in g.dirs().iter().enumerate() {
        let reaches = e.reaches(*dir);
        let live = reaches.get(..n).unwrap_or(&[]);

        // Support function: the furthest any atom's surface reaches this way.
        // The comparison form rather than a `max_by`, so the traversal order is
        // fixed and no float comparison depends on iterator internals (§13.1).
        let mut extent = Span(f64::NEG_INFINITY);
        for reach in live {
            if *reach > extent {
                extent = *reach;
            }
        }

        // Surface character: atoms within SHELL of the supporting plane,
        // weighted by how close they are to it. Fixed iteration order, so the
        // summation order is fixed (§13.1).
        let floor = extent - SHELL;
        let (mut num, mut den) = (0.0f64, 0.0f64);
        for (reach, character) in live.iter().zip(characters) {
            let above = (*reach - floor).get();
            // `if a > b { a } else { b }`, not `f64::max` — §13.1 bans the
            // method as non-deterministic for ±0.0.
            let w = if above > 0.0 { above } else { 0.0 };
            num += w * character;
            den += w;
        }

        if let (Some(extent_slot), Some(character_slot)) = (out.r.get_mut(d), out.a.get_mut(d)) {
            *extent_slot = extent;
            // `den > 0.0` cannot be false for `n >= 1`: the extremal atom has
            // `reach == extent`, so its weight is SHELL. Kept as a guard rather
            // than an assertion because the alternative is a division whose
            // failure mode is a silent NaN that then wins every downstream
            // selection.
            //
            // **There is deliberately no clamp to `[-1, 1]` here**, and the
            // reason is measured rather than stylistic. `num / den` is a convex
            // combination of values §7.1 draws in `[-1, 1]`, so it cannot leave
            // the range mathematically; the residual worry is a one-ulp
            // overshoot from the division. Measured over 600 molecules at
            // n = 1..12 including 300 maximally symmetric ones — 25 200
            // (direction, molecule) samples — the worst overshoot is **exactly
            // 0.0**, so a clamp would never have fired.
            //
            // The decisive argument is not that it is dead, though. A clamp
            // *masks* a §7.1 violation: an element whose affinity left `[-1, 1]`
            // would be silently corrected into range and §8.3's kernel would
            // receive a wrong-but-plausible value. The range is a contract this
            // module *depends* on, not one it should quietly repair, so it is
            // asserted in `extents_are_positive_and_characters_are_in_range`
            // over a real corpus instead — where a violation fails loudly.
            *character_slot = if den > 0.0 { num / den } else { 0.0 };
        }
    }
    out
}

#[cfg(test)]
#[expect(
    clippy::indexing_slicing,
    clippy::unwrap_used,
    reason = "CLAUDE.md permits both inside #[cfg(test)] with a stated reason; every index \
              here is a loop bound over D or over an atom count already established"
)]
mod tests {
    use super::*;
    use crate::canonical::{CanonMol, canonicalise};
    use crate::geodesic::{Geodesic, Rotation, is_identity, rotation_matrices};
    use crate::graph::Mol12;
    use crate::layout::embed;
    use borbax_rng::{Domain, Stream};
    use borbax_units::Span;
    use borbax_universe::{
        BondOrder, ElementId, PeriodicTable, Universe, element::generate_elements,
    };

    /// §22.2's prior. The ladder is {12, 42, 162} and nothing between them keeps
    /// rotation an exact permutation.
    const D: usize = 42;

    /// Point inversion, `p -> -p`. Improper (determinant -1), so it is not one
    /// of the 60 and no `permuted` can reproduce it.
    const INVERT: [[f64; 3]; 3] = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];

    /// A periodic table and a universe **from the same seed**. `table(a)` with
    /// `Universe::generate(b)` for `a != b` is a chimera — `layout.rs`'s
    /// `fixture` records what that cost there, and the same helper exists here
    /// so the mismatch is unspellable rather than merely avoided.
    fn fixture(seed: u64) -> (PeriodicTable, Universe) {
        (generate_elements(seed), Universe::generate(seed))
    }

    fn chain_capable(tbl: &PeriodicTable) -> Vec<ElementId> {
        (0..120usize)
            .filter_map(ElementId::from_index)
            .filter(|id| tbl.get(*id).is_some_and(|el| el.valence >= 2))
            .collect()
    }

    fn canon(mol: &Mol12) -> CanonMol {
        canonicalise(mol)
            .unwrap_or_else(|err| unreachable!("fixture capped: {err}"))
            .0
    }

    fn geo() -> Geodesic<D> {
        Geodesic::<D>::build().unwrap_or_else(|_| unreachable!("D=42 builds"))
    }

    /// The whole pipeline for a raw graph: canonicalise, embed, sample.
    fn sig(mol: &Mol12, uni: &Universe, g: &Geodesic<D>) -> Signature<D> {
        signature(&embed(&canon(mol), uni), g)
    }

    fn chain(n: u8, tbl: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = Mol12::new();
        for i in 0..n {
            assert!(mol.add_atom(ids[usize::from(i) % ids.len()]).is_some());
        }
        for i in 0..n.saturating_sub(1) {
            assert!(mol.add_bond(i, i + 1, BondOrder::SINGLE, tbl).is_ok());
        }
        mol
    }

    fn random_tree(rng: &mut Stream, n: u8, tbl: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = Mol12::new();
        assert!(mol.add_atom(ids[0]).is_some());
        for i in 1..n {
            let pick = usize::try_from(rng.next_range(u64::try_from(ids.len()).unwrap())).unwrap();
            let Some(child) = mol.add_atom(ids[pick]) else {
                break;
            };
            let parent = u8::try_from(rng.next_range(u64::from(i))).unwrap();
            if mol.add_bond(parent, child, BondOrder::SINGLE, tbl).is_err() {
                for other in 0..i {
                    if mol.add_bond(other, child, BondOrder::SINGLE, tbl).is_ok() {
                        break;
                    }
                }
            }
        }
        mol
    }

    /// A ring-bearing molecule built from **one element throughout**, which
    /// maximises the automorphism group and is what reaches `embed`'s
    /// coincidence branch — automorphic same-element atoms are the ones with
    /// identical target rows. Deterministic given `n`, so it takes no `Stream`.
    fn symmetric_molecule(n: u8, tbl: &PeriodicTable, one: ElementId) -> Mol12 {
        let mut mol = Mol12::new();
        for _ in 0..n {
            if mol.add_atom(one).is_none() {
                break;
            }
        }
        let size = u8::try_from(mol.len()).unwrap_or(0);
        for i in 1..size {
            let _ = mol.add_bond(i - 1, i, BondOrder::SINGLE, tbl);
        }
        if size > 2 {
            let _ = mol.add_bond(0, size - 1, BondOrder::SINGLE, tbl);
        }
        mol
    }

    fn random_molecule(rng: &mut Stream, n: u8, tbl: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = random_tree(rng, n, tbl, ids);
        let size = u8::try_from(mol.len()).unwrap();
        for _ in 0..rng.next_range(4) {
            let a = u8::try_from(rng.next_range(u64::from(size))).unwrap();
            let b = u8::try_from(rng.next_range(u64::from(size))).unwrap();
            let _ = mol.add_bond(a, b, BondOrder::SINGLE, tbl);
        }
        mol
    }

    /// A molecule with `from`'s **exact element multiset**, rewired from
    /// scratch. Removes the composition confound from the locality test by
    /// construction rather than measuring against it — `layout.rs` records why
    /// that replacement was forced.
    fn rewired(rng: &mut Stream, from: &Mol12, tbl: &PeriodicTable) -> Option<Mol12> {
        let mut out = Mol12::new();
        for i in 0..u8::try_from(from.len()).ok()? {
            out.add_atom(from.element(i)?)?;
        }
        let size = u8::try_from(out.len()).ok()?;
        let mut order: Vec<u8> = (0..size).collect();
        for i in (1..order.len()).rev() {
            let j = usize::try_from(rng.next_range(u64::try_from(i + 1).ok()?)).ok()?;
            order.swap(i, j);
        }
        for pair in order.windows(2) {
            let _ = out.add_bond(pair[0], pair[1], BondOrder::SINGLE, tbl);
        }
        for _ in 0..rng.next_range(4) {
            let a = u8::try_from(rng.next_range(u64::from(size))).ok()?;
            let b = u8::try_from(rng.next_range(u64::from(size))).ok()?;
            let _ = out.add_bond(a, b, BondOrder::SINGLE, tbl);
        }
        Some(out)
    }

    /// The **plain** Euclidean distance, which the public API deliberately does
    /// not offer (routed requirement 1).
    ///
    /// It exists here and only here, because the property "the group-minimised
    /// distance is invariant under a group element and the plain one is not"
    /// cannot be stated without both halves. A test that measures the contrast
    /// is the reason the public surface may omit it; deleting this would make
    /// the omission unfalsifiable.
    fn plain_distance<const N: usize>(a: &Signature<N>, b: &Signature<N>) -> f64 {
        let mut acc = 0.0f64;
        for i in 0..N {
            let dr = (a.extents()[i] - b.extents()[i]).get();
            let da = a.characters()[i] - b.characters()[i];
            acc += dr * dr + da * da;
        }
        acc.sqrt()
    }

    // ---------------------------------------------------------------------
    // The support function, tied to geometry rather than to itself
    // ---------------------------------------------------------------------

    /// One atom sits at the origin, so its surface reaches exactly its radius
    /// in **every** direction. Checked against the element table, not against
    /// another call to `signature`.
    #[test]
    fn one_atom_reaches_exactly_its_radius_in_every_direction() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        for id in ids.iter().take(8) {
            let mut mol = Mol12::new();
            assert!(mol.add_atom(*id).is_some());
            let want = uni.element(*id).unwrap().radius;
            let s = sig(&mol, &uni, &g);
            for i in 0..D {
                assert!(
                    (s.extents()[i] - want).get().abs() < 1e-12,
                    "direction {i}: {:?} != {want:?}",
                    s.extents()[i]
                );
            }
        }
    }

    /// The extent is the support function of the **union of atom spheres**:
    /// the furthest any atom's surface reaches. Hand-computed from the two
    /// atoms' positions and radii rather than compared against another
    /// signature.
    #[test]
    fn the_extent_is_the_support_function_of_the_atom_spheres() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mol = chain(2, &tbl, &ids);
        let species = canon(&mol);
        let emb = embed(&species, &uni);
        let s = signature(&emb, &g);

        for i in 0..D {
            let dir = g.dirs()[i];
            let mut want = f64::NEG_INFINITY;
            for a in 0..emb.len() {
                let p = emb.coords()[a];
                let radius = uni
                    .element(species.mol().element(u8::try_from(a).unwrap()).unwrap())
                    .unwrap()
                    .radius
                    .get();
                let reach = p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2] + radius;
                if reach > want {
                    want = reach;
                }
            }
            assert!(
                (s.extents()[i].get() - want).abs() < 1e-12,
                "direction {i}: {} != {want}",
                s.extents()[i].get()
            );
        }
    }

    /// **The load-bearing link between the geometry and the lookup table, and
    /// the one test that decides `permuted`'s convention.**
    ///
    /// Rotating the embedded molecule bodily by the matrix for `r` and
    /// recomputing the signature from scratch must give exactly
    /// `original.permuted(rotation_perms(r))`. Every other rotation test in
    /// this file compares tables against tables, which is the right invariant
    /// and says nothing about what the composition *means* — the failure
    /// `contact_perms`' docstring shipped with at Task 7.
    ///
    /// If this fails, `canonicalise` is minimising over a table that does not
    /// correspond to rotating anything, and its answer is meaningless while
    /// looking entirely plausible.
    #[test]
    fn rotating_the_embedding_permutes_the_signature() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(901, Domain::Molecule, 0);
        let species = canon(&random_tree(&mut rng, 9, &tbl, &ids));
        let emb = embed(&species, &uni);
        let base = signature(&emb, &g);
        let mats = rotation_matrices().unwrap_or_else(|_| unreachable!("the matrices build"));

        for r in Rotation::all() {
            let turned = emb.rotated(&mats[r.index()]);
            let got = signature(&turned, &g);
            let want = base.permuted(g.rotation_perms(r));
            for i in 0..D {
                assert!(
                    (got.extents()[i] - want.extents()[i]).get().abs() < 1e-9,
                    "rotation {r}, direction {i}: extent {} != {}",
                    got.extents()[i].get(),
                    want.extents()[i].get()
                );
                assert!(
                    (got.characters()[i] - want.characters()[i]).abs() < 1e-9,
                    "rotation {r}, direction {i}: character {} != {}",
                    got.characters()[i],
                    want.characters()[i]
                );
            }
        }
    }

    /// **`anti` tied to the geometry, which is what makes the mirror in
    /// §22.8's test a mirror.**
    ///
    /// Point-inverting the embedding (`p -> -p`) must give exactly the original
    /// signature read through the antipodal permutation. Without this, the
    /// enantiomer test below is mirroring an index table and calling it a
    /// reflection.
    #[test]
    fn point_inverting_the_embedding_is_the_antipodal_permutation() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(902, Domain::Molecule, 0);
        let species = canon(&random_tree(&mut rng, 9, &tbl, &ids));
        let emb = embed(&species, &uni);
        let base = signature(&emb, &g);

        let inverted = signature(&emb.rotated(&INVERT), &g);
        let want = base.permuted(g.anti());
        for i in 0..D {
            assert!(
                (inverted.extents()[i] - want.extents()[i]).get().abs() < 1e-9,
                "direction {i}: extent {} != {}",
                inverted.extents()[i].get(),
                want.extents()[i].get()
            );
            assert!(
                (inverted.characters()[i] - want.characters()[i]).abs() < 1e-9,
                "direction {i}: character {} != {}",
                inverted.characters()[i],
                want.characters()[i]
            );
        }
        // The mutation must be real: a point inversion that changed nothing
        // would make this pass while proving nothing.
        assert!(
            plain_distance(&base, &inverted) > 1e-6,
            "the fixture is centrosymmetric, so this test proves nothing"
        );
    }

    /// The surface character reads **only** the atoms within [`SHELL`] of the
    /// supporting plane, weighted linearly by how close they sit to it.
    /// Hand-computed against the same weights rather than against another call.
    #[test]
    fn the_character_averages_only_the_atoms_facing_that_way() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(903, Domain::Molecule, 0);
        let species = canon(&random_tree(&mut rng, 10, &tbl, &ids));
        let emb = embed(&species, &uni);
        let s = signature(&emb, &g);

        let mut excluded_seen = 0u32;
        for i in 0..D {
            let dir = g.dirs()[i];
            let reaches = emb.reaches(dir);
            let extent = s.extents()[i];
            let (mut num, mut den) = (0.0f64, 0.0f64);
            for (reach, character) in reaches.iter().zip(emb.characters()) {
                let raw = (*reach - (extent - SHELL)).get();
                let w = if raw > 0.0 { raw } else { 0.0 };
                if w == 0.0 {
                    excluded_seen += 1;
                }
                num += w * character;
                den += w;
            }
            assert!(den > 0.0, "direction {i}: no atom faces it at all");
            assert!(
                (s.characters()[i] - num / den).abs() < 1e-12,
                "direction {i}: {} != {}",
                s.characters()[i],
                num / den
            );
        }
        // A shell that admitted every atom would make this a plain average and
        // the weighting inert. Counted over (direction, atom) pairs, which is
        // what the loop above increments.
        // **The bar is set from a sweep, not by eye.** Exclusions on this
        // fixture against the shell width: 0.75 -> 355 of 420 (shipped),
        // 1.5 -> 307, 2.0 -> 257, 3.0 -> 180, 6.0 -> 12, 1000 -> 0. So 250
        // trips at roughly SHELL >= 2.1 — the shell having grown past the
        // scale it was chosen at — while leaving 30% headroom below the
        // measured value. The corpus is seed-fixed, so this is exact and no
        // flakiness budget is being spent.
        assert!(
            excluded_seen > 250,
            "only {excluded_seen} of 420 (direction, atom) pairs fell outside the shell; \
             SHELL is admitting effectively everything, so the character channel has \
             become a composition average with no geometry in it"
        );
    }

    // ---------------------------------------------------------------------
    // Determinism and range
    // ---------------------------------------------------------------------

    #[test]
    fn is_deterministic() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mol = chain(6, &tbl, &ids);
        assert_eq!(sig(&mol, &uni, &g), sig(&mol, &uni, &g));
    }

    /// **This is the only guard on the `[-1, 1]` contract §8.3's kernel
    /// depends on, so it runs over a corpus rather than over one molecule per
    /// size.** `signature` deliberately does not clamp — a clamp would mask a
    /// §7.1 violation rather than surface it — which puts the whole weight of
    /// the contract here.
    ///
    /// The corpus alternates ordinary molecules with maximally symmetric ones
    /// (a single element throughout, ring-closed), because those are what reach
    /// `embed`'s coincidence branch and so are where a degenerate character
    /// average would appear.
    #[test]
    fn extents_are_positive_and_characters_are_in_range() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(904, Domain::Molecule, 0);
        let mut examined = 0u32;
        let mut symmetric_seen = 0u32;
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for trial in 0..300u32 {
            let atoms = 1 + u8::try_from(trial % 12).unwrap();
            let mol = if trial % 2 == 0 {
                random_molecule(&mut rng, atoms, &tbl, &ids)
            } else {
                symmetric_seen += 1;
                symmetric_molecule(atoms, &tbl, ids[0])
            };
            let shape = sig(&mol, &uni, &g);
            for i in 0..D {
                let extent = shape.extents()[i];
                assert!(
                    extent.is_finite() && extent.get() > 0.0,
                    "trial {trial} n={atoms} direction {i}: extent {extent:?}"
                );
                let character = shape.characters()[i];
                assert!(
                    character.is_finite() && (-1.0..=1.0).contains(&character),
                    "trial {trial} n={atoms} direction {i}: character {character} left \
                     [-1, 1]. `signature` does not clamp on purpose, so this is a real \
                     §7.1 range violation and not a rounding artefact — measured worst \
                     overshoot over this corpus is exactly 0.0"
                );
                if character < lo {
                    lo = character;
                }
                if character > hi {
                    hi = character;
                }
                examined += 1;
            }
        }
        // Counts (molecule, direction) samples, which is what the inner loop
        // increments: 300 x 42.
        assert_eq!(examined, 300 * 42, "the corpus did not run to completion");
        assert_eq!(
            symmetric_seen, 150,
            "the degenerate half of the corpus is not being built"
        );

        // **The attained range, pinned — because containment in `[-1, 1]`
        // cannot fire and Task 10's carried-in note says exactly why.**
        //
        // Measured on **this** corpus: `[0.477560, 0.888889]`. The broader
        // evidence is a separate sweep — 25 universes x 40 molecules x 42
        // directions — where the channel attains `[0.214442, 0.916667]`,
        // exactly the element-affinity range over those seeds, and **never
        // straddles zero**. (The two are different corpora and the first draft
        // of this block pinned the sweep's numbers against this test, which
        // failed immediately and correctly.) That is Task 4's finding arriving
        // downstream: §7.1's `affinity` is exhaustively `[0.140476, 1.0]` over
        // the `(k, units)` grid, and a weighted average of positive values is
        // positive. §8.3's charge term `-(a_A + a_B)^2` is maximised at zero
        // when the two are *opposite*, so an all-positive channel degenerates
        // it from a complementarity test into a monotone penalty on total
        // surface affinity, and §5's membrane mechanism — one flank positive,
        // the opposite negative — is unreachable.
        //
        // **This assertion pins the bounds; it deliberately does not assert
        // that they fail to straddle zero.** A test asserting the defect would
        // publish it as the contract and make the fix look like a regression.
        // Pinning the numbers instead makes the remedy *visible*: centring or
        // mean-subtracting the channel fails this line and sends the reader to
        // the note. The remedy is Task 10's to take, with its own discriminator
        // and against a binding kernel that does not exist yet; taking it here,
        // blind, would be tuning the physics to make a behaviour appear.
        assert!(
            (lo - 0.477_560).abs() < 5e-6 && (hi - 0.888_889).abs() < 5e-6,
            "the character channel attained [{lo:.6}, {hi:.6}], not the pinned \
             [0.477560, 0.888889]. If you have just centred or mean-subtracted it per \
             Task 10's carried-in note on §8.3's charge term, this is the expected \
             failure and the bounds want updating — check the new range straddles zero. \
             If you have not, something upstream in §7.1's affinity draw has moved."
        );
    }

    // ---------------------------------------------------------------------
    // Routed requirement 1 — every comparison is group-minimised
    // ---------------------------------------------------------------------

    /// The canonical form must be invariant under the rotation group. If it is
    /// not, the same molecule can intern as two species depending on how its
    /// embedding happened to land (§8.2).
    #[test]
    fn the_canonical_form_is_rotation_invariant() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(905, Domain::Molecule, 0);
        let base = sig(&random_tree(&mut rng, 9, &tbl, &ids), &uni, &g);
        let want = base.canonicalise(&g);
        for r in Rotation::all() {
            assert_eq!(
                base.permuted(g.rotation_perms(r)).canonicalise(&g),
                want,
                "rotation {r} broke invariance"
            );
        }
    }

    /// **The wiring proof, and it needs no corpus** (routed requirement 1).
    ///
    /// [`Signature::group_distance`] is invariant under applying any group
    /// element to either argument; the plain Euclidean distance is not. The
    /// withdrawn discriminator asked for a concordance comparison against a
    /// fixed rotation — measured 2-2-1 across five seeds at the sizes `Mol12`
    /// admits, so it would have reported a correct implementation as unwired.
    ///
    /// Three preconditions, and omitting any one makes it vacuous: the rotation
    /// must be **non-identity**, the fixture must be **asymmetric**, and the
    /// rotation must be asserted to have **actually changed** the signature.
    #[test]
    fn the_group_distance_is_rotation_invariant_and_a_plain_one_is_not() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(906, Domain::Molecule, 0);
        let base = sig(&random_tree(&mut rng, 9, &tbl, &ids), &uni, &g);
        let mats = rotation_matrices().unwrap_or_else(|_| unreachable!("the matrices build"));

        let (mut checked, mut symmetric) = (0u32, 0u32);
        for r in Rotation::all() {
            if mats.get(r.index()).is_some_and(is_identity) {
                continue;
            }
            let turned = base.permuted(g.rotation_perms(r));
            // The mutation must be a real change, or everything below is
            // trivially true.
            let plain = plain_distance(&base, &turned);
            if plain < 1e-9 {
                // This rotation is in the fixture's own symmetry group, so it
                // moved neither distance and proves nothing either way.
                symmetric += 1;
                continue;
            }
            checked += 1;
            assert!(
                base.group_distance(&turned, &g) < 1e-9,
                "rotation {r} moved the group distance by {}, plain {plain}",
                base.group_distance(&turned, &g)
            );
        }
        // **One counter, not two.** An earlier version incremented a second
        // `plain_moved` on the line below this one and asserted on it
        // separately, which made it identical to `checked` by construction — a
        // bar that read as an independent check and measured a quantity its own
        // message did not describe. `checked` already carries both facts: each
        // counted rotation moved the plain distance (it survived the `continue`)
        // and left the group distance at zero (it survived the assertion).
        // Measured here: 59 of 59 non-identity rotations, 0 symmetric.
        assert!(
            checked > 55,
            "only {checked} of 59 non-identity rotations genuinely moved the fixture \
             ({symmetric} were in its own symmetry group); it is too symmetric to prove \
             anything, and the contrast this test rests on is not there"
        );
    }

    // ---------------------------------------------------------------------
    // Routed requirement 2 — §22.8, and it has had no test anywhere in V0
    // ---------------------------------------------------------------------

    /// **A chiral shape and its mirror image are different species (§22.8).**
    ///
    /// This is the spec's actual handedness claim and it has had no test
    /// anywhere in V0. It cannot be made at Task 8 — a Borbax graph has exactly
    /// one embedding, so the enantiomer is not constructible from a graph — so
    /// it is a statement about *signatures*, and it lands here because
    /// [`Signature::canonicalise`] minimises over the 60 **proper** rotations,
    /// which do not contain `-I`.
    ///
    /// Three assertions, because the outcome alone would pass for the wrong
    /// reason. A `canonicalise` that returned `self` unchanged would separate
    /// *every* pair, chiral or not:
    ///
    /// 1. the fixture is chiral **by definition** — no proper rotation carries
    ///    the mirror onto the original — established without `canonicalise`;
    /// 2. the canonical forms differ, so the quotient preserves handedness;
    /// 3. extending the search to all 120 elements **collapses them**, which is
    ///    the mechanism rather than the symptom. That is the plan's probe made
    ///    permanent instead of run once.
    #[test]
    fn a_chiral_signature_and_its_mirror_are_different_species() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(907, Domain::Molecule, 0);

        let mut chiral_fixtures = 0u32;
        for trial in 0..40u32 {
            let n = 8 + u8::try_from(trial % 5).unwrap();
            let s = sig(&random_tree(&mut rng, n, &tbl, &ids), &uni, &g);
            let mirrored = s.permuted(g.anti());

            // (1) Chiral by definition: no proper rotation carries the mirror
            // onto the original. Computed from `permuted` alone.
            let achiral = Rotation::all()
                .any(|r| plain_distance(&s, &mirrored.permuted(g.rotation_perms(r))) < 1e-9);
            if achiral {
                continue;
            }
            chiral_fixtures += 1;

            // (2) The canonical form keeps them apart.
            assert_ne!(
                s.canonicalise(&g),
                mirrored.canonicalise(&g),
                "trial {trial}: a chiral signature and its mirror interned as one species"
            );

            // (3) The mechanism: the 120-element search is exactly what would
            // destroy §22.8. Minimising over the improper coset as well brings
            // the two to the same representative.
            let full = |x: &Signature<D>| {
                let direct = x.canonicalise(&g);
                let flipped = x.permuted(g.anti()).canonicalise(&g);
                if flipped.lex_cmp(&direct) == core::cmp::Ordering::Less {
                    flipped
                } else {
                    direct
                }
            };
            assert_eq!(
                full(&s),
                full(&mirrored),
                "trial {trial}: the 120-element search did not collapse the enantiomers, \
                 so assertion (2) is not evidence that the search is 60"
            );
        }
        // Measured: 40 of 40 fixtures are chiral, so this is a corpus
        // self-check rather than a threshold to tune.
        assert!(
            chiral_fixtures > 35,
            "only {chiral_fixtures} of 40 fixtures were chiral; §22.8 is being asserted \
             on a corpus that mostly cannot express it"
        );
    }

    // ---------------------------------------------------------------------
    // Routed requirement 4 — locality, now on the real signature
    // ---------------------------------------------------------------------

    /// **§8.2's founding property, measured on the signature itself.**
    ///
    /// `layout.rs` measures this on a proxy — atom *centres*, no radii and no
    /// affinity channel — which its own doc calls deliberately less informative
    /// so that it understates locality rather than flattering it. This is the
    /// real descriptor, and routed requirement 4 makes it Task 9's to own.
    ///
    /// Scored on [`Signature::group_distance`]. The stranger carries the
    /// mutant's **exact element multiset**, so a composition-only statistic is
    /// 0.5 structurally and the margin below is graph locality and nothing else.
    ///
    /// **The group minimisation is not what makes this number good, and the
    /// measurement says so plainly.** On this corpus the plain distance scores
    /// **0.9194** against `group_distance`'s **0.8896** — the minimisation
    /// *loses*. That is not a defect and it is not new: the routed requirement
    /// records `D_group` against `D_raw` at the sizes `Mol12` admits as 2 wins,
    /// 2 losses and 1 exact tie across five seeds, which is why the withdrawn
    /// discriminator for requirement 1 — "the same test with a fixed rotation
    /// must score measurably worse" — would have rejected a correct
    /// implementation. The reading is that the landmark start already lands
    /// parent and mutant in near-identical frames, so minimising over 60 poses
    /// cannot help and occasionally hands the *stranger* a flattering pose.
    ///
    /// So `group_distance` is justified by three things that are not this one:
    /// §8.2's stored lex-min form is bit-identically absorbed only under it,
    /// §8.3's `affinity` maximises over the group so the orbit *is* the
    /// consumer's equivalence class, and a gate finer than the consumer's
    /// equivalence class ranks implementations backwards — which `layout.rs`
    /// measured and recorded. Requirement 1's discriminator is therefore the
    /// invariance property in
    /// `the_group_distance_is_rotation_invariant_and_a_plain_one_is_not`, not
    /// this concordance. **Do not add a bar here asserting group beats plain.**
    #[test]
    fn a_one_atom_edit_moves_the_signature_less_than_an_unrelated_molecule_does() {
        let (tbl, uni) = fixture(17);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(908, Domain::Molecule, 0);

        let (mut nearer, mut trials, mut cyclic) = (0u32, 0u32, 0u32);
        let (mut null_hits, mut null_trials) = (0u32, 0u32);
        for _ in 0..1200 {
            let n = 6 + u8::try_from(rng.next_range(5)).unwrap();
            let parent = random_molecule(&mut rng, n, &tbl, &ids);
            let mut mutant = parent;
            let pick = usize::try_from(rng.next_range(u64::try_from(ids.len()).unwrap())).unwrap();
            let Some(fresh) = mutant.add_atom(ids[pick]) else {
                continue;
            };
            let anchor = u8::try_from(rng.next_range(u64::from(n))).unwrap();
            if mutant
                .add_bond(anchor, fresh, BondOrder::SINGLE, &tbl)
                .is_err()
            {
                continue;
            }
            let Some(stranger) = rewired(&mut rng, &mutant, &tbl) else {
                continue;
            };
            let (cp, cm, cs) = (canon(&parent), canon(&mutant), canon(&stranger));
            if cp.len() < 3 {
                continue;
            }
            let bonds: u32 = (0..u8::try_from(cp.len()).unwrap_or(0))
                .map(|a| cp.mol().degree(a))
                .sum::<u32>()
                / 2;
            if bonds > u32::try_from(cp.len()).unwrap_or(0).saturating_sub(1) {
                cyclic += 1;
            }

            let base = signature(&embed(&cp, &uni), &g);
            let dm = base.group_distance(&signature(&embed(&cm, &uni), &g), &g);
            let ds = base.group_distance(&signature(&embed(&cs, &uni), &g), &g);
            if dm < ds {
                nearer += 1;
            }
            trials += 1;

            // The composition null, read only from the element multiset. With a
            // composition-matched stranger it is 0.5 *structurally*; half credit
            // for a tie is the AUC convention and is load-bearing, because the
            // two sums are exactly equal here and scoring only strict `<` would
            // make the null 0.0 — a margin over which is no bar at all.
            let comp = |sp: &CanonMol| -> f64 {
                let mut total = 0.0f64;
                for i in 0..sp.len() {
                    total += u8::try_from(i)
                        .ok()
                        .and_then(|k| sp.mol().element(k))
                        .and_then(|id| uni.element(id))
                        .map_or(0.0, |el| el.radius.get());
                }
                total
            };
            let (bc, mc, sc) = (comp(&cp), comp(&cm), comp(&cs));
            let (nm, ns) = ((bc - mc).abs(), (bc - sc).abs());
            if nm < ns {
                null_hits += 2;
            } else if (nm - ns).abs() <= f64::EPSILON {
                null_hits += 1;
            }
            null_trials += 1;
        }

        assert!(trials > 600, "only {trials} usable trials");
        // Measured: 315 of 806.
        assert!(
            cyclic > 150,
            "only {cyclic} of {trials} fixtures had a cycle; the corpus is effectively \
             tree-only and this test is not measuring what it claims"
        );
        let concordance = f64::from(nearer) / f64::from(trials);
        let null = f64::from(null_hits) / (2.0 * f64::from(null_trials));
        assert!(
            (null - 0.5).abs() < 1e-12,
            "the composition null is {null:.4}, not 0.5 — the stranger is no longer \
             composition-matched, so the margin below is measuring the confound"
        );
        assert!(
            concordance - null > 0.20,
            "signature locality is not beating a composition-only null: concordance \
             {concordance:.4} ({nearer}/{trials}), null {null:.4}, margin {:.4}. Before \
             blaming the signature, check whether the losing pairs are the same \
             molecules `layout.rs` reports a 1.1-5.5% R_max convergence tail on at \
             n >= 11 (routed requirement 4).",
            concordance - null
        );
        // A regression guard, pinned against a measured value rather than a
        // tuned threshold: the shipped signature scores 0.8896 here and the
        // corpus is seed-fixed, so this is exact. It defends against the
        // signature or the embedding degrading; it does **not** defend the
        // group minimisation, which would score 0.9194 if removed and sail
        // past. The invariance test is what defends that.
        assert!(
            concordance > 0.85,
            "signature locality regressed: concordance {concordance:.4} \
             ({nearer}/{trials}), against 0.8896 measured for the shipped signature"
        );
    }

    // ---------------------------------------------------------------------
    // Span-typing (routed requirement 5)
    // ---------------------------------------------------------------------

    /// The extent is a length and carries [`Span`], so the G1 mix-up
    /// `coords()` leaves open cannot be spelled through the signature: §8.3
    /// compares `r` against `ideal_gap`, also a `Span`.
    #[test]
    fn the_extent_is_a_span_and_scales_with_the_molecule() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let small = sig(&chain(3, &tbl, &ids), &uni, &g);
        let large = sig(&chain(11, &tbl, &ids), &uni, &g);
        let biggest = |s: &Signature<D>| {
            let mut best = Span(f64::NEG_INFINITY);
            for i in 0..D {
                if s.extents()[i] > best {
                    best = s.extents()[i];
                }
            }
            best
        };
        assert!(
            biggest(&large) > biggest(&small),
            "an 11-atom chain does not reach further than a 3-atom one: {:?} vs {:?}",
            biggest(&large),
            biggest(&small)
        );
    }
}
