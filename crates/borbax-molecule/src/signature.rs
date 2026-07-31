//! The shape signature (spec §8.2) — the heart of the engine.
//!
//! A molecule reduces to two numbers per sample direction: how far its surface
//! extends that way ([`Signature::extents`]), and what surface character faces
//! that way ([`Signature::characters`]). Everything downstream — binding,
//! catalysis, membranes, permeability — is a comparison between two of these.
//!
//! # Every comparison is group-minimised, and nothing enforces it
//!
//! [`Signature::group_distance`] is the **only** distance this module offers.
//! That is a convention, not a guarantee, and the distinction is this section's
//! whole point.
//!
//! **Why it matters.** §8.2 stores a signature as the lexicographically
//! smallest of its 60 rotations, and that form is the **worst** of five
//! descriptors G2 compared — 0.7510–0.7745 size-matched, against `D_group`'s
//! 0.9230–0.9415. Re-measured on the shipped signature (seed 17, 806 trials,
//! composition-matched stranger): `D_lexmin` **0.7953**, `D_group` **0.8896**,
//! `D_raw` 0.9194, `D_sorted` 0.8102. The ordering that matters survives —
//! lex-min is clearly worst — but note that those two rows are **different
//! descriptors**: G2's are the atom-centre proxy with no radii and no character
//! channel, so the bands are not this module's numbers and are quoted here only
//! for provenance.
//!
//! Storing lex-min is nonetheless safe, and exactly conditionally so. Lex-min is
//! itself a group element, so `D_group(lexmin(A), lexmin(B))` equals
//! `D_group(A, B)` — measured: the worst `group_distance` between a molecule and
//! 60 real rotations of its own embedding is **1.13e-14**. It costs nothing
//! *provided every consumer compares group-minimised*, and is severe for any
//! consumer that subtracts two stored signatures directly. §8.2 names three:
//! species identity and the fold cache read the stored form for **equality**,
//! which [`Signature::canonicalise`] serves and which is fine; §15.2's novelty
//! histogram takes **distances**, and is the one that would be hurt.
//!
//! **How it would be hurt matters, because the obvious answer is wrong.** The
//! damage is *not* that a poor descriptor inflates distances — that is a bias,
//! and §15.3's neutral shadow is exactly a bias-removal device, so a reviewer
//! would rightly answer "the shadow handles that". What survives `run - shadow`
//! is the loss of **discriminability**: a representative-choosing map does not
//! merely raise distances, it maps some genuinely-near pairs to far while
//! leaving far pairs far, destroying the *ordering*. A shadow subtracts a
//! location parameter; it cannot restore ordering lost to noise. The second
//! surviving mechanism is a within-run temporal flip — a persistent species
//! whose representative changes between timesteps reports novelty for something
//! that has not changed, and the shadow does not reproduce that flip pattern so
//! the residue does not cancel. §15.3's own stated limits are about persistence
//! and RAF, not about representation, so it was never offered as a control here.
//!
//! So §8.2 needs no amendment. The storage form is absorbed by the comparison
//! rule.
//!
//! ## There is no guard, and that is deliberate rather than an oversight
//!
//! An `xtask` check shipped here briefly and was **deleted**. It read
//! `signature.rs`'s public surface and required that any `pub fn` mentioning
//! `Signature` twice also take a [`Geodesic`]. Review found **four independent
//! escapes** — an out-of-line `mod`, an `impl` in any other file of the crate,
//! an `impl` nested in a function body, and `Self` counted per-argument where a
//! named `Signature` was counted per-mention — and **two classes of false
//! positive**, one of which was a `Geodesic` *receiver*, i.e. `g.affinity(&a,
//! &b)`, which is how §8.3's kernel reads most naturally.
//!
//! The rule this project applies to a repaired guard is: state the class of
//! thing it enumerates and name what is outside it; if that sentence cannot be
//! written, the repair is another enumeration. It could not be written. The
//! guard enumerated AST shapes.
//!
//! It was also **pointed the wrong way**. Its corpus caught `impl PartialOrd for
//! Signature` — but §8.2's own consumers need an ordering, and a
//! determinism-safe intern table is a `BTreeMap`, which needs `Ord`. So the
//! guard forbade what the spec requires, and the cheapest way past it was to
//! move the `impl` to another file, which it could not see. A guard that fires
//! on correct code gets deleted, and the pressure to evade this one was
//! structural.
//!
//! **What holds the property today is this paragraph and code review.** That is
//! weaker than a check and is stated plainly rather than dressed up: a green
//! build is *not* evidence that no ungrouped comparison exists. Task 10 owns the
//! replacement, where the binding kernel makes the right shape visible — an
//! API-surface snapshot (`public-api`-style) is crate-wide by construction and
//! would close all four escapes, which an AST walker over one file cannot.
//!
//! ## Why there is no `sorted_extents`
//!
//! The plan's interface list carries one; it is deliberately not implemented.
//! The point in its favour is real: as a *filter* it is sound, since two
//! signatures in the same proper orbit have equal sorted extents, and being
//! invariant under a larger group only makes a filter more permissive — the safe
//! direction.
//!
//! What decides against it is that §8.3 does not ask for a filter. It asks for
//! "a rotation-invariant **upper bound on the score**", and adds that "a filter
//! that cannot state the property it guarantees is not conservative, it is
//! merely untested". A sorted extent vector is not that bound, and building it
//! needs Task 10's binding constants. Meanwhile a *distance* on sorted extents
//! identifies a shape with its mirror — in Kendall's terms it is the reflection
//! shape space `RΣ` (quotient by `O(m)`) where `D_group` is the shape space `Σ`
//! (quotient by `SO(m)`), so §22.8 would be impossible by construction. Routed
//! to Task 10.

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
/// It sets how much of the molecule the character channel can see. Measured
/// mean atoms admitted per direction (mean n = 8.37): `SHELL = 0` -> **0.000**,
/// `0.1` -> 1.056, `0.75` -> **1.572**, `3.0` -> 5.214.
///
/// **Three statements that stood here are false and are corrected rather than
/// deleted, because each would mislead a sweep.** (1) "At zero only the single
/// extremal atom contributes" — at zero the extremal atom's own weight is
/// `extent - extent = 0`, so `den == 0` and the channel is identically zero in
/// every direction. (2) "0.75 admits the extremal atom's immediate neighbours"
/// overstates by about two: it admits the extremal atom plus ~0.57 of one
/// neighbour, against a 10th-percentile interatomic distance of 1.656. (3)
/// "Above the molecular diameter ... the channel becomes a composition average
/// with no geometry in it" — the weight is `reach - extent + SHELL`, linear in
/// `reach` with slope 1 for *every* `SHELL`, so it never becomes uniform; at
/// `SHELL = 50` the character channel still scores 0.7343 on the locality
/// statistic, which a composition average could not.
///
/// **It is chosen, not derived**, and named as such for the reason `bonds.rs`
/// records about the strain exponent. Two things Task 20 needs before sweeping
/// it. First, a **prior**: the natural scale is angular, not a length — a
/// direction owns a patch of half-angle 15.86° at D = 42 (nearest-neighbour
/// spacing 31.717°, uniform over all 42 vertices), whose depth on a body of
/// radius `R` is `R(1 - cos 15.86°) = 0.0381 R`, giving ~0.114 for V0 extents
/// of ~3. The shipped 0.75 corresponds to `R ~ 19.7`, roughly 6x larger than
/// `Mol12` can build. That derivation is in-house, not a literature value.
///
/// Second, an **instrument**, because the obvious one is blind: whole-signature
/// concordance varies by 2 trials in 399 across `SHELL` in `[0, 50]`. The
/// character-channel-only concordance does vary (0.727-0.792, peaking at 3.0),
/// so the sweep needs a channel-isolating statistic — and note that peak sits
/// *outside* the bracket `the_character_averages_only_the_atoms_facing_that_way`
/// admits, so the two criteria disagree and Task 20 must say which wins.
///
/// A further property worth keeping: `SHELL` is an absolute length while the
/// embedding scales with the drawn radius series, so it means somewhat
/// different things in different universes (1.24 atoms admitted at seed 10
/// against 1.78 at seed 11, a 1.43x spread). Within a universe it is
/// size-*stable* — 1.64 atoms at n = 2 against 1.57 at n = 12 while max extent
/// grows 1.72 -> 7.90 — so it does not leak molecule size into the descriptor.
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
    /// The tie-break runs over the character channel as well as the extents,
    /// and **it fires** — which an earlier version of this comment denied on the
    /// strength of a corpus that could not see it.
    ///
    /// Why it is needed: with an extents-only order, "first among equals in
    /// `Rotation::all()` order" is *not* a function of the orbit. `S` and `S∘g`
    /// present the same candidate multiset indexed differently, so when the
    /// minimal set holds two candidates with different characters the two
    /// enumerations return different answers and `canonicalise` stops being
    /// well-defined.
    ///
    /// That is not hypothetical. Measured over **all 900 dimers** on the first
    /// 30 chain-capable elements of seed 6, an extents-only `lex_cmp` gives a
    /// canonical form that is **not orbit-invariant on 40 of them** (26 at seed
    /// 0, 25 at seed 17); over 1 680 random molecules across 14 universes the
    /// hits are 9, **all at n = 2**. In a pose-pair census, 360 of ~90 000
    /// extent ties differ in the character channel.
    ///
    /// **The withdrawn claim, kept because the way it failed is instructive.**
    /// This comment previously read "it has never once decided anything ... in
    /// **0** of them do the characters differ", from a 600-molecule corpus
    /// drawing `n = 1 + trial % 12`. That gives ~50 dimers, and at a 4.4% hit
    /// rate a zero is entirely likely — and ~97% of the ties counted were
    /// monomers, where all 60 poses are trivially identical. It was a sampling
    /// artefact reported as a property, and it argued for deleting code that
    /// keeps species identity well-defined on the size class a beaker is mostly
    /// made of. `the_character_tie_break_is_reached` now pins it with an
    /// exhaustive dimer enumeration rather than a random draw.
    ///
    /// The mechanism is not understood: the failing dimers do **not** have equal
    /// radii (0.485 against 1.123 in one case), so the obvious explanation is
    /// wrong and no guess is recorded here in its place.
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
/// planned two-pass form did — at `D = 162` that is 162 x 12 = 1 944 pairs, so
/// **1 944 duplicated** dot products of 3 888 total. (An earlier version of this
/// sentence, and the plan it came from, called 3 888 the duplicate count, which
/// doubles the saving.)
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
        // surface affinity. Measured: `corr(min_R C, D(m_A+m_B)^2) = 0.99991`,
        // so 99.98% of the between-pair variance in the charge term is just the
        // two molecules' mean affinities — almost *only* the monotone penalty.
        // §10.1's membrane mechanism (one flank positive, the opposite
        // negative) is therefore unreachable. Cited as §5 until now, which is
        // the fiction guarantees; §10.1 puts it more strongly than the
        // paraphrase — amphiphile-analogues "are a region of signature space,
        // and the sim discovers which of its molecules live there", and with an
        // all-positive channel that region is provably empty.
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
             [0.477560, 0.888889]. If you have just taken Task 10's carried-in remedy \
             on §8.3's charge term, this is the expected failure and the bounds want \
             updating — but check WHICH remedy: centring the channel here is measurably \
             inert (an exact rotation-invariant constant shift, 2400/2400) and the real \
             edit is the affine remap at borbax-universe's affinity draw. If you have \
             changed nothing, something upstream in §7.1 has moved."
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
        let mats = rotation_matrices().unwrap_or_else(|_| unreachable!("the matrices build"));

        // **n >= 3 only, and the exclusion of n = 2 is a measured finding, not
        // a convenience.** At n = 2 the canonical form is genuinely unstable
        // under rotating the embedding: 24 of 60 poses give a different
        // representative, worst extent movement **0.352** — not float noise, a
        // third of a Span. At n >= 3 the movement is 8.9e-16 to 1.8e-15, i.e.
        // rounding and nothing else. The dimer case is asserted separately
        // below with its consequence bounded rather than pinned green here.
        for n in 3..=10u8 {
            let species = canon(&random_tree(&mut rng, n, &tbl, &ids));
            let emb = embed(&species, &uni);
            let want = signature(&emb, &g).canonicalise(&g);
            for r in Rotation::all() {
                // **Rotate the embedding and recompute**, do not permute the
                // table. This is the whole content of the test and an earlier
                // version got it wrong: `base.permuted(rotation_perms(r))` is a
                // bit-exact relabelling, so `canonicalise` of it equals
                // `canonicalise` of the original *algebraically*, for any
                // implementation that minimises over a closed group. No
                // arithmetic on coordinates is involved, so nothing could fail
                // it — while its own doc claimed to guard against an embedding
                // "happening to land" differently, which is precisely a claim
                // about that arithmetic.
                let got = signature(&emb.rotated(&mats[r.index()]), &g).canonicalise(&g);
                for i in 0..D {
                    assert!(
                        (got.extents()[i] - want.extents()[i]).get().abs() < 1e-9
                            && (got.characters()[i] - want.characters()[i]).abs() < 1e-9,
                        "n={n} rotation {r} direction {i}: the canonical form moved when \
                         the molecule was rotated, so one species can intern as two"
                    );
                }
            }
        }
    }

    /// **A dimer's canonical form is not stable under rotation, and this pins
    /// what that does and does not cost.**
    ///
    /// Measured: rotating a dimer's embedding by the 60 matrices and
    /// recanonicalising gives a different representative on **24 of 60** poses,
    /// worst extent movement **0.352**. The cause is §8.2's own hazard —
    /// choosing a representative jumps wherever the choice is degenerate — and
    /// it is the reason Task 8 deleted `canonicalise_frame`. A dimer's sampled
    /// signature has near-ties in the lex order between poses that are genuinely
    /// different signatures, so last-bit noise selects between them.
    ///
    /// **What it costs, bounded here rather than left as a worry.**
    /// [`Signature::group_distance`] absorbs the representative choice exactly,
    /// so binding, novelty and every distance-taking consumer are unaffected —
    /// that is what this test asserts. What is *not* absorbed is **equality**,
    /// so species identity for a dimer is well-defined only because `embed` is
    /// deterministic for a given species. Two consequences follow and both are
    /// routed rather than fixed here: a cross-platform last-bit difference in
    /// `embed` could flip the representative and mint a second species from one
    /// molecule (Task 20's golden matrix would report that with no visible
    /// cause), and Task 12's folding is the first thing that holds two
    /// conformations of one chain at once.
    ///
    /// A dimer is the first product of every condensation, so this is not an
    /// exotic size class.
    #[test]
    fn a_dimer_canonical_form_is_unstable_but_its_distance_is_not() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mats = rotation_matrices().unwrap_or_else(|_| unreachable!("the matrices build"));
        let mut rng = Stream::new(905, Domain::Molecule, 0);
        let species = canon(&random_tree(&mut rng, 2, &tbl, &ids));
        let emb = embed(&species, &uni);
        let base = signature(&emb, &g);

        let (mut moved, mut worst_form, mut worst_distance) = (0u32, 0.0f64, 0.0f64);
        for r in Rotation::all() {
            let turned = signature(&emb.rotated(&mats[r.index()]), &g);
            // The representative may move...
            let (a, b) = (turned.canonicalise(&g), base.canonicalise(&g));
            let mut form = 0.0f64;
            for i in 0..D {
                let d = (a.extents()[i] - b.extents()[i]).get().abs();
                if d > form {
                    form = d;
                }
            }
            if form > 1e-9 {
                moved += 1;
            }
            if form > worst_form {
                worst_form = form;
            }
            // ...but the distance must not.
            let d = base.group_distance(&turned, &g);
            if d > worst_distance {
                worst_distance = d;
            }
        }
        assert!(
            moved > 10 && worst_form > 0.1,
            "the dimer discontinuity has gone away ({moved} of 60 poses moved, worst \
             {worst_form}); if `canonicalise` was made continuous this test should be \
             deleted rather than relaxed — measured 24 of 60, worst 0.352"
        );
        assert!(
            worst_distance < 1e-9,
            "group_distance did NOT absorb the representative choice ({worst_distance}) — \
             the discontinuity has escaped into every distance-taking consumer, which is \
             the whole reason §8.2's storage form is considered safe"
        );
    }

    /// **The character tie-break is reached, and a random corpus cannot show
    /// it.** Enumerating every dimer over the first 30 chain-capable elements,
    /// an extents-only ordering gives a canonical form that is *not* a function
    /// of the orbit on 40 of 900 — while the shipped `lex_cmp`, which continues
    /// into the character channel, is orbit-invariant on all of them.
    ///
    /// This exists because the doc on `canonicalise` previously asserted the
    /// opposite from a corpus containing ~50 dimers, and that assertion argued
    /// for deleting the comparison. The discriminator is not "the shipped form
    /// is well-defined" — it always is — but "the extents-only form is not".
    #[test]
    fn the_character_tie_break_is_reached() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let ids = chain_capable(&tbl);
        let take = 30.min(ids.len());

        // The extents-only order the shipped one extends, built here so the
        // difference between them is the thing under test.
        let extents_only = |a: &Signature<D>, b: &Signature<D>| {
            for (mine, theirs) in a.extents().iter().zip(b.extents().iter()) {
                match mine.canonical_cmp(theirs) {
                    Ordering::Equal => {}
                    ord => return ord,
                }
            }
            Ordering::Equal
        };
        let canon_by =
            |x: &Signature<D>, cmp: &dyn Fn(&Signature<D>, &Signature<D>) -> Ordering| {
                let mut best: Option<Signature<D>> = None;
                for rot in Rotation::all() {
                    let cand = x.permuted(g.rotation_perms(rot));
                    if best
                        .as_ref()
                        .is_none_or(|held| cmp(&cand, held) == Ordering::Less)
                    {
                        best = Some(cand);
                    }
                }
                best.unwrap_or_else(|| unreachable!("Rotation::all() yields 60 elements"))
            };

        let (mut split, mut examined) = (0u32, 0u32);
        for a in 0..take {
            for b in 0..take {
                let mut mol = Mol12::new();
                if mol.add_atom(ids[a]).is_none() || mol.add_atom(ids[b]).is_none() {
                    continue;
                }
                if mol.add_bond(0, 1, BondOrder::SINGLE, &tbl).is_err() {
                    continue;
                }
                let base = sig(&mol, &uni, &g);
                examined += 1;
                // **Every pose, not one.** Orbit-invariance is a statement
                // about the whole orbit, and probing a single rotation finds 2
                // of 900 where the full check finds 40 — which is how a version
                // of this test that looked at `Rotation::new(7)` alone came to
                // report a bar of 20 as unreachable.
                let from_base = canon_by(&base, &extents_only);
                let shipped = canon_by(&base, &|x, y| x.lex_cmp(y));
                let mut any = false;
                for rot in Rotation::all() {
                    let moved = base.permuted(g.rotation_perms(rot));
                    if canon_by(&moved, &extents_only) != from_base {
                        any = true;
                    }
                    // The shipped order must always be orbit-invariant.
                    assert_eq!(
                        canon_by(&moved, &|x, y| x.lex_cmp(y)),
                        shipped,
                        "the shipped canonicalise is not a function of the orbit for \
                         elements {a}/{b} under rotation {rot}"
                    );
                }
                if any {
                    split += 1;
                }
            }
        }
        assert!(
            examined > 800,
            "only {examined} dimers built, of a possible 900"
        );
        assert!(
            split > 20,
            "the extents-only ordering split the orbit on only {split} of {examined} \
             dimers; if this is 0 the corpus has stopped reaching the case and \
             `canonicalise`'s character tie-break is again undefended (measured: 40)"
        );
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

    /// **A chiral shape and its mirror image are different species — with
    /// respect to the 60-element quotient (§22.8).**
    ///
    /// This is the spec's handedness claim and it has had no test anywhere in
    /// V0. It cannot be made at Task 8 — a Borbax graph has exactly one
    /// embedding, so the enantiomer is not constructible from a graph — so it is
    /// a statement about *signatures*, and it lands here because
    /// [`Signature::canonicalise`] minimises over the 60 **proper** rotations,
    /// which do not contain `-I`.
    ///
    /// **The scope of "chiral" here is narrower than §22.8's prose and the gap
    /// is measured.** §22.8 makes a continuous claim — a left-handed helix
    /// "genuinely cannot be rotated onto a right-handed one, *no matter how it
    /// is turned*". What this test establishes is chirality with respect to the
    /// icosahedral quotient: that no member of the **60** carries the mirror
    /// onto the original. Those differ, and the discrete version
    /// **over-separates**. Measured on seed 6 over 1 500 molecules: 513
    /// embeddings are planar, and **368 of them — 71.7% — are classified
    /// chiral**, at margins of 0.36 to 0.66, eight orders above this test's
    /// tolerance. Planarity implies achirality by proof, not by measurement:
    /// reflection in the occupied plane fixes every atom centre, hence every
    /// atom sphere, hence the support function and the character channel. The
    /// 60-element search sees it only when that reflection happens to coincide
    /// with an improper icosahedral element on the sampled directions, which is
    /// generically false. The rate falls with size — 155/155 at n = 3, 1/119 at
    /// n = 10 — and this corpus is n = 8..12 trees with **0 of 40** planar, so
    /// the assertions below are sound; it is the word "chiral" that needs the
    /// qualifier.
    ///
    /// **Consequence, routed rather than fixed.** Two mirror-image planar
    /// conformations would intern as two species when they are one shape. That
    /// is unreachable in V0 — a graph has one embedding, and over 300 molecules
    /// there are **0 pairs** where one molecule's mirror signature equals
    /// another's — and becomes live at Task 12, where folding first holds two
    /// conformations of one chain at once. Making the quotient continuous is a
    /// design change §8.2 does not currently ask for.
    ///
    /// In Kendall's terms this is the shape space `Σ = S/SO(m)`; the rejected
    /// `D_sorted` is the reflection shape space `RΣ = S/O(m)`, which identifies
    /// an object with its mirror. That is the standard distinction and is why
    /// the rejection is not a matter of taste.
    ///
    /// Two assertions, not three. A third — that widening the search to 120
    /// elements collapses the pair — **was a tautology and has been deleted**:
    /// with `full(x) = min_lex(canon(x), canon(mirror(x)))` and `anti` an
    /// involution, `full(mirror(x)) = full(x)` for *any* `canon` whatsoever,
    /// including one that returns `self`. It tested only that `anti` is an
    /// involution, which `geodesic.rs` already tests, and its failure message
    /// described a state that cannot occur. Verified by mutation: widening
    /// `canonicalise` to 120 elements fails assertion (2) on 40 of 40 fixtures
    /// and passes the deleted assertion on all 40. That the search is 60 and not
    /// 120 is a theorem about the group, not something a test can discover, and
    /// `geodesic.rs::minus_identity_is_not_in_the_group` is where it belongs.
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

            // (1) Chiral with respect to the 60: no proper rotation carries the
            // mirror onto the original. Computed from `permuted` alone, with no
            // `canonicalise` in the loop, so it cannot inherit that function's
            // verdict. Conservative in the safe direction — a near-achiral
            // fixture is skipped, not asserted on.
            let achiral = Rotation::all()
                .any(|r| plain_distance(&s, &mirrored.permuted(g.rotation_perms(r))) < 1e-9);
            if achiral {
                continue;
            }
            chiral_fixtures += 1;

            // (2) The canonical form keeps them apart. This is what carries the
            // test: it fails on 40 of 40 if `canonicalise` searches all 120.
            assert_ne!(
                s.canonicalise(&g),
                mirrored.canonicalise(&g),
                "trial {trial}: a chiral signature and its mirror interned as one species"
            );
        }
        // Measured: 40 of 40, and 0 of 40 planar.
        assert!(
            chiral_fixtures > 35,
            "only {chiral_fixtures} of 40 fixtures were chiral under the 60-element \
             search; §22.8 is being asserted on a corpus that mostly cannot express it"
        );
    }

    // ---------------------------------------------------------------------
    // Routed requirement 4 — locality, now on the real signature
    // ---------------------------------------------------------------------

    /// **§8.2's founding property, measured on the signature itself.**
    ///
    /// Concordance is the fraction of trials where a one-atom mutant is nearer
    /// its parent than an unrelated molecule is: 1.0 perfect, **0.5 none**.
    /// Measured **0.8896** over 806 trials against a structural 0.5, which is
    /// `z = 22.1` against `Binomial(806, 0.5)`. §8.2's property is not in doubt.
    ///
    /// **The baseline is confounded and the margin is smaller than it looks.**
    /// `rewired` builds the stranger as a random Hamiltonian **path**, while
    /// parents come from `random_molecule`, a random **tree** plus extras. Mean
    /// graph diameter: mutant **5.34**, stranger **6.85**, with the stranger
    /// longer in 561 trials and shorter in 126. So the statistic partly rewards
    /// a descriptor for noticing elongation. Two independent controls: on the
    /// diameter-matched subset concordance is **0.7479** (n = 119), and with the
    /// stranger drawn from the parent's own generator it is **0.7747** (n = 821).
    /// The property survives both comfortably — the `> 0.20` margin holds — but
    /// roughly 0.11 to 0.14 of the headline is the stranger generator. The
    /// regression bar below therefore has far less slack than 0.8896 against
    /// 0.85 suggests. Inherited from `layout.rs`, which uses the same helper.
    ///
    /// **The composition confound is removed by construction**: the stranger
    /// carries the mutant's exact element multiset, so the composition null is
    /// 0.5 **bitwise** in 806 of 806 trials. That makes `concordance - null`
    /// algebraically `concordance > 0.70`, so the subtraction is decorative and
    /// the null's real job is the self-check on the line above it — one
    /// mis-scored trial moves it by ~6e-4 and fires the assertion.
    ///
    /// # The plain distance beats the group-minimised one, and by how much
    ///
    /// **0.9194 against 0.8896 on this corpus, and plain wins on 8 of 8 universe
    /// seeds** (sign test, two-sided `p = 0.0078`). Paired on the shipped 806
    /// trials by **`McNemar`**'s test: plain-wins-group-loses `b = 34`,
    /// group-wins-plain-loses `c = 10`, `z = +3.618`, exact two-sided
    /// `p = 0.000388`. The unpaired two-proportion `z` gives 2.034, `p = 0.042` —
    /// so the routed requirement's warning that an unpaired test understates is
    /// confirmed on this very comparison.
    ///
    /// Two corrections to how this was previously written up. The plan's
    /// "2 wins, 2 losses, 1 exact tie" is **G2's atom-centre proxy**, not this
    /// descriptor, and quoting it here made a systematic effect look like a
    /// wash. And the paired-test citation should be **`McNemar` (1947),
    /// Psychometrika 12(2):153-157** — `DeLong` et al. (1988) is a U-statistic
    /// covariance estimator for *AUCs on shared individuals*, which is the right
    /// family and the wrong tool for a per-trial binary outcome under two
    /// descriptors.
    ///
    /// **The mechanism, measured.** The minimisation removes 3.28% of the mean
    /// parent-to-mutant distance and **8.89%** of the parent-to-stranger
    /// distance; on the 34 trials where it loses, the mutant shrinks 0.71% and
    /// the stranger **27.76%** — a 39x asymmetry. A quotient metric collapses
    /// more of a large distance than a small one, so it is systematically
    /// generous to the more distant candidate. (It is a *min* over 60, not a
    /// max; an earlier version of this note said max.)
    ///
    /// **This is not an argument for switching to the plain distance.**
    /// `group_distance` is justified by three things that are not locality:
    /// §8.2's stored lex-min form is absorbed bit-identically only under it
    /// (measured 1.13e-14), §8.3's `affinity` maximises over the group so the
    /// orbit *is* the consumer's equivalence class, and a gate finer than the
    /// consumer's equivalence class ranks implementations backwards — which
    /// `layout.rs` measured and recorded. Requirement 1's discriminator is the
    /// invariance property in
    /// `the_group_distance_is_rotation_invariant_and_a_plain_one_is_not`, not
    /// this concordance. **Do not add a bar here asserting group beats plain.**
    ///
    /// # On `layout.rs`'s proxy, whose numbers may not be quoted beside these
    ///
    /// `layout.rs` scores 0.8787 with an atom-centre proxy and documents it as
    /// "deliberately *less* informative ... so it understates locality rather
    /// than flattering it". **Measured paired on identical trials, that is
    /// false**: proxy **0.8908** against this signature's **0.8896** (`McNemar`
    /// `z = -0.19`), and 0.8871 against 0.8797 on a second stream. The proxy is
    /// level or marginally ahead. The 0.0109 gap between the two shipped figures
    /// is corpus difference (different stream, 1 327 trials against 806); the
    /// unpaired comparison gives `p = 0.45`. So radii and the character channel
    /// buy nothing measurable *on this statistic* — which is consistent with the
    /// character channel being 2.08% of the distance, and is a reason to want a
    /// channel-isolating statistic rather than a reason to drop either.
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

    /// **A bit-level digest of the shipped signature, and the reason it is bits
    /// and not a tolerance.**
    ///
    /// Three arithmetic shapes in this module are pinned in prose with nothing
    /// enforcing them: `reaches`' `((p0*d0 + p1*d1) + p2*d2) + radius`,
    /// `targets`' `h * (r_i + r_j)`, and `group_distance`'s direction-order
    /// accumulation. Every one of them is *algebraically* invariant under the
    /// tidy-up a reader would reach for, so the suite stays green while every
    /// number in the system moves in its last bits.
    ///
    /// **Measured, which is why a tolerance assertion would be worthless here.**
    /// Reassociating `reaches` to `radius + p0*d0 + p1*d1 + p2*d2` moves
    /// **10 747 of 42 000** samples in `to_bits()`, affecting 915 of 1 000
    /// species — and **0 of 42 000** by more than 1e-12. Indexing
    /// `group_distance` forward instead of through the inverse permutation moves
    /// **433 of 1 000** pairs in bits and **0** by more than 1e-12. A test with
    /// any tolerance at all passes both wrong versions 100% of the time.
    ///
    /// So this exists before any optimisation is attempted, not after. Two are
    /// known available and deliberately not taken here — reindexing
    /// `group_distance` through the inverse permutation (2.3x release, 4.4x dev)
    /// and giving `reaches` a caller-owned buffer (2.15x) — both bit-identical
    /// by construction and measured so. This digest is what would make taking
    /// them safe; there is no consumer yet that makes them worth taking.
    ///
    /// It hashes through `canonical_bits`, not `to_bits`, because a runtime
    /// NaN's sign and `-0.0` are architecture-dependent and §13.6's matrix would
    /// report that as a simulation divergence (§13.4).
    ///
    /// **If this fails and you did not mean to change the physics**, the cause
    /// is almost certainly an accumulation order, not a formula.
    #[test]
    fn the_signature_digest_is_pinned() {
        let g = geo();
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let mut values = 0u64;
        for seed in 0..8u64 {
            let (tbl, uni) = fixture(seed);
            let ids = chain_capable(&tbl);
            if ids.is_empty() {
                continue;
            }
            let mut rng = Stream::new(4242 + seed, Domain::Molecule, 0);
            for trial in 0..25u32 {
                let atoms = 1 + u8::try_from(trial % 12).unwrap();
                let mol = if trial % 2 == 0 {
                    random_molecule(&mut rng, atoms, &tbl, &ids)
                } else {
                    symmetric_molecule(atoms, &tbl, ids[0])
                };
                // The canonical form, because that is what §8.2 stores and what
                // Task 20's state hash will read.
                let shape = sig(&mol, &uni, &g).canonicalise(&g);
                for i in 0..D {
                    for bits in [
                        borbax_units::canonical_bits(shape.extents()[i].get()),
                        borbax_units::canonical_bits(shape.characters()[i]),
                    ] {
                        hash ^= bits;
                        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                        values += 1;
                    }
                }
            }
        }
        // Counts the f64 values hashed, which is what the inner loop
        // increments: 8 seeds x 25 molecules x 42 directions x 2 channels.
        assert_eq!(
            values,
            8 * 25 * 42 * 2,
            "the corpus did not run to completion"
        );
        assert_eq!(
            hash, 0xb303_62ed_8102_4ae2,
            "the shape signature moved. If you meant to change the physics, \
             regenerate this constant and say so in the commit message; if you did \
             not, suspect an accumulation order — `reaches`, `targets` and \
             `group_distance` each carry a pinned association that is \
             algebraically invariant under the obvious tidy-up and moves 11-22% of \
             values in their last bits."
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
