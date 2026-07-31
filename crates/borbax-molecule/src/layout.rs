//! 3D embedding by stress majorization (spec §8.2).
//!
//! SMACOF rather than Laplacian eigenvectors: eigenvector sign is arbitrary and
//! degenerate eigenvalues are common in the symmetric graphs small molecules
//! actually are, so that arbitrariness would propagate straight into species
//! identity and the fold cache. Stress majorization from a deterministic start
//! sidesteps the problem rather than managing it.
//!
//! **No randomness at all**, and no transcendental. The Guttman transform needs
//! only `sqrt` and the four arithmetic operations, and the starting layout is a
//! fixed rational table rather than a trigonometric spiral — precisely so that
//! stays true (§13.1). The iteration count is fixed rather than
//! convergence-tested, because a tolerance makes the *number* of iterations
//! depend on float details that vary by platform, and a configuration that
//! stopped one iteration earlier is a different configuration (§13.4). That is
//! also why no solver crate fits: none offers "stop after exactly n iterations
//! regardless of convergence".
//!
//! # `embed` takes a [`CanonMol`], and that is the whole design
//!
//! Everything from here on is a pure function of the *species* (§8.6), but the
//! obvious signature — `embed(&Mol12, &Universe)` — is a function of the
//! caller's atom numbering, which `Mol12` does not constrain. It would compile,
//! run, and return a plausible shape that differs between two callers who built
//! the same molecule by different routes. Taking a [`CanonMol`] leaves that call
//! no spelling; see the type's own docs.
//!
//! # What this deliberately does *not* do: canonicalise the frame
//!
//! The V0 plan's Step 3 for this task ends by rotating the result into a frame
//! built from the lowest-indexed non-degenerate atoms, with a comment arguing
//! that canonicalising the signature over the 60 icosahedral rotations "cannot
//! absorb" frame drift, so without the frame "mutation becomes a random walk in
//! signature space \[and\] neutral networks never form".
//!
//! **That argument is falsified by measurement, and the measurement is checked
//! in.** `experiments/`'s G2 harness scores four descriptors on how well
//! signature distance tracks graph mutation — concordance 1.0 is perfect
//! locality, 0.5 is none — over 2000 trials × 5 seeds with a positive control
//! (composition) and a negative control (hash of the canonical form), both
//! passing. `experiments/src/signature.rs`'s `canonical_frame` is the identical
//! construction to the plan's, so `D_frame` measures exactly this proposal:
//!
//! ```text
//!               free size          size-matched (shape only)
//! D_group     0.9660 - 0.9730      0.9230 - 0.9415   best
//! D_raw       0.9560 - 0.9665      0.9075 - 0.9355   no alignment at all
//! D_sorted    0.9450 - 0.9535      0.8840 - 0.9075
//! D_frame     0.9195 - 0.9380      0.8585 - 0.8885   worst
//! ```
//!
//! The frame is **worse than doing nothing**, on every seed, in both regimes,
//! and it is not marginal: it inflates the median parent-to-mutant distance from
//! 1.698 to 2.494, because the frame jumps whenever the canonical order does.
//! `D_group` — minimise over the 60 rotations — is best, which is also *already
//! what the binding kernel does*, so relying on it removes a mechanism instead of
//! adding one (§8.3's "one mechanism, reused"). Re-derived rather than taken from
//! the note: `cargo run --locked -p borbax-experiments --release --bin g2`.
//!
//! So the locality property lives at Task 9, as a `D_group` requirement on
//! signatures, and is not restated here as a weaker proxy over atom coordinates
//! that no longer correspond across an edit. Do **not** reintroduce a frame here
//! without re-running G2 and beating 0.9560 / 0.9075.

use crate::canonical::CanonMol;
use crate::graph::{MAX_ATOMS, Mol12};
use borbax_units::Span;
use borbax_universe::Universe;

/// Iterations of the Guttman transform. Fixed for the reason in the module docs;
/// changing it is a physics change that moves every downstream number.
pub const ITERATIONS: usize = 240;

/// A molecule's atoms placed in space, indexed by **canonical position**.
///
/// Lengths are in [`Span`] units — the targets are built from
/// `Element::radius`, which is a `Span`. The coordinate array itself stays `f64`
/// because it is the solver's workspace and a squared length has no `Span`
/// spelling; the typed boundary is [`Self::distance`] and
/// [`Self::radius_of_gyration`], which is where a length is compared against
/// another length and therefore where the G1 class of mix-up actually happens.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Embedding {
    n: usize,
    /// Slots past `n` are zero padding and are never read — [`Self::coords`]
    /// hands out the live prefix so a caller cannot index into them by accident.
    pos: [[f64; 3]; MAX_ATOMS],
}

impl Embedding {
    /// How many atoms are placed.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.n
    }

    /// Whether nothing is placed.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The placed coordinates, without the padding.
    #[must_use]
    pub fn coords(&self) -> &[[f64; 3]] {
        self.pos.get(..self.n).unwrap_or(&[])
    }

    /// Distance between two atoms, `None` if either index is unplaced.
    #[must_use]
    pub fn distance(&self, i: usize, j: usize) -> Option<Span> {
        let (left, right) = (self.pos.get(i)?, self.pos.get(j)?);
        if i >= self.n || j >= self.n {
            return None;
        }
        let sep = [left[0] - right[0], left[1] - right[1], left[2] - right[2]];
        Some(Span(
            (sep[0] * sep[0] + sep[1] * sep[1] + sep[2] * sep[2]).sqrt(),
        ))
    }

    /// RMS distance from the centroid — the molecule's own length scale.
    ///
    /// The embedding is centred, so this is taken about the origin.
    #[must_use]
    pub fn radius_of_gyration(&self) -> Span {
        if self.n == 0 {
            return Span::ZERO;
        }
        let mut sum = 0.0f64;
        for point in self.coords() {
            sum += point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
        }
        #[expect(
            clippy::cast_precision_loss,
            clippy::as_conversions,
            reason = "n <= MAX_ATOMS = 12, so the widening is exact"
        )]
        let inv_n = 1.0 / self.n as f64;
        Span((sum * inv_n).sqrt())
    }
}

/// All-pairs shortest path by BFS, in hops. `n <= 12`, so this is trivially
/// cheap and there is no reason to reach for anything cleverer.
///
/// **Disconnected pairs come back as one hop beyond the graph's diameter, not as
/// an infinity**, which is the Task 2 review's routed hazard. `experiments`'
/// version returns `f64::INFINITY` and is safe there only because every fixture
/// is a tree; an infinity here would take stress to `inf - inf` and the whole
/// embedding to NaN, which then wins every downstream selection silently.
fn hop_distances(mol: &Mol12) -> [[u8; MAX_ATOMS]; MAX_ATOMS] {
    let n = mol.len();
    let mut hops = [[u8::MAX; MAX_ATOMS]; MAX_ATOMS];
    for src in 0..n {
        if let Some(cell) = hops.get_mut(src).and_then(|row| row.get_mut(src)) {
            *cell = 0;
        }
        let mut frontier: u16 = 1u16 << u32::try_from(src).unwrap_or(0);
        let mut seen: u16 = frontier;
        let mut depth = 0u8;
        while frontier != 0 {
            depth = depth.saturating_add(1);
            let mut next: u16 = 0;
            for node in 0..n {
                let Ok(node8) = u8::try_from(node) else {
                    continue;
                };
                if frontier & (1u16 << u32::from(node8)) != 0 {
                    next |= mol.neighbours(node8);
                }
            }
            next &= !seen;
            seen |= next;
            for node in 0..n {
                if next & (1u16 << u32::try_from(node).unwrap_or(0)) != 0
                    && let Some(cell) = hops.get_mut(src).and_then(|row| row.get_mut(node))
                {
                    *cell = depth;
                }
            }
            frontier = next;
        }
    }

    let diameter = hops
        .iter()
        .flatten()
        .filter(|x| **x != u8::MAX)
        .copied()
        .max()
        .unwrap_or(1);
    for row in hops.iter_mut().take(n) {
        for cell in row.iter_mut().take(n) {
            if *cell == u8::MAX {
                *cell = diameter.saturating_add(1);
            }
        }
    }
    hops
}

/// Target distance matrix: hop count scaled by the two atoms' radii, so bigger
/// atoms genuinely take up more room.
fn targets(species: &CanonMol, universe: &Universe) -> [[f64; MAX_ATOMS]; MAX_ATOMS] {
    let mol = species.mol();
    let n = mol.len();
    let hops = hop_distances(mol);
    let radius = |i: usize| -> f64 {
        u8::try_from(i)
            .ok()
            .and_then(|k| mol.element(k))
            .and_then(|id| universe.element(id))
            .map_or(1.0, |el| el.radius.get())
    };
    let mut target = [[0.0f64; MAX_ATOMS]; MAX_ATOMS];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let h = hops.get(i).and_then(|r| r.get(j)).copied().unwrap_or(1);
            if let Some(cell) = target.get_mut(i).and_then(|r| r.get_mut(j)) {
                *cell = f64::from(h) * (radius(i) + radius(j));
            }
        }
    }
    target
}

/// The objective [`embed`] minimises: `Σ_{i<j} w_ij (δ_ij − d_ij)²`, with the
/// same `w_ij = 1/d_ij²` weighting the solver uses.
///
/// **This exists because it is the one guarantee stress majorization makes.** A
/// misderived Guttman transform still converges to *something*, so no test of
/// the output alone would notice; only a monotonicity check distinguishes "the
/// optimiser runs" from "the optimiser works".
///
/// **It must be the *weighted* sum, and getting that wrong produces a
/// convincing false positive.** Written first as the raw sum
/// `Σ (δ_ij − d_ij)²`, it reported stress rising — trial 2, step 101,
/// `12.780925441019596 → 12.780936033504108` — which reads exactly like a
/// misderived transform and is not one. Weighted SMACOF decreases the
/// *weighted* objective; the raw sum is a different function and is under no
/// obligation to fall with it. The tell was that substituting uniform weights,
/// where the two objectives are proportional, made the "defect" vanish — which
/// is evidence about the objective, not about the solver. Raw stress is
/// available as [`raw_stress`] for reading off how far the distances actually
/// are; it is the interpretable one and the wrong one to assert monotonicity on.
///
/// Summed in `i < j` index order and never through `.sum()` over an unordered
/// iterator, so the accumulation order is fixed on every platform (§13.1).
#[must_use]
pub fn stress(species: &CanonMol, universe: &Universe, emb: &Embedding) -> f64 {
    let target = targets(species, universe);
    let n = emb.len();
    let mut total = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let Some(val) = emb.distance(i, j) else {
                continue;
            };
            let tgt = target.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
            if tgt <= 0.0 {
                continue;
            }
            let err = val.get() - tgt;
            total += (err * err) / (tgt * tgt);
        }
    }
    total
}

/// The unweighted mismatch `Σ_{i<j} (δ_ij − d_ij)²`, in squared [`Span`] units.
///
/// The interpretable measure of how far the realised distances are from the
/// targets. **Not** the quantity the solver minimises — see [`stress`], which is.
#[must_use]
pub fn raw_stress(species: &CanonMol, universe: &Universe, emb: &Embedding) -> f64 {
    let target = targets(species, universe);
    let n = emb.len();
    let mut total = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let Some(val) = emb.distance(i, j) else {
                continue;
            };
            let tgt = target.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
            let err = val.get() - tgt;
            total += err * err;
        }
    }
    total
}

/// Deterministic starting configuration.
///
/// These rational offsets are neither collinear nor coplanar, which is all the
/// solver requires — it does not need a good sphere, only a spread-out and
/// reproducible one. Spelling it as a table rather than a spiral is what keeps
/// `sin`/`cos` out of this file (§13.1), and using a table at all rather than
/// drawing from a `borbax_rng::Stream` is what makes the embedding a function
/// of the species alone.
const START: [[f64; 3]; MAX_ATOMS] = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.25, 0.5],
    [0.25, 1.0, -0.5],
    [-0.5, -1.0, 0.25],
    [0.5, 0.25, 1.0],
    [-0.25, 0.5, -1.0],
    [1.0, -0.75, -0.25],
    [-1.0, -0.5, 0.75],
    [0.75, 1.0, 0.5],
    [-0.75, 0.5, -1.0],
    [0.5, -1.0, 0.75],
    [-0.5, 0.75, 1.0],
];

/// [`START`] scaled to the molecule's own size.
///
/// Scaled by the mean over *all* pairs, not by row 0: one atom's row is not the
/// molecule's length scale, and using it would make the start depend on which
/// atom happens to be first.
fn initial_layout(target: &[[f64; MAX_ATOMS]; MAX_ATOMS], n: usize) -> [[f64; 3]; MAX_ATOMS] {
    let mut sum = 0.0f64;
    let mut cnt = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                sum += target.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                cnt += 1.0;
            }
        }
    }
    let mean_target = if cnt > 0.0 { sum / cnt } else { 1.0 };
    let mut pos = [[0.0f64; 3]; MAX_ATOMS];
    for i in 0..n {
        for k in 0..3 {
            if let (Some(dst), Some(src)) = (
                pos.get_mut(i).and_then(|p| p.get_mut(k)),
                START.get(i).and_then(|p| p.get(k)),
            ) {
                *dst = src * mean_target;
            }
        }
    }
    pos
}

/// Lay a species out in three dimensions.
///
/// Takes a [`CanonMol`] rather than a [`Mol12`] so the result cannot depend on
/// the caller's labelling — see the module docs.
#[must_use]
pub fn embed(species: &CanonMol, universe: &Universe) -> Embedding {
    embed_with_budget(species, universe, ITERATIONS)
}

/// [`embed`] with the iteration count exposed.
///
/// Only the fixed-budget [`embed`] is public: a caller who can choose the budget
/// can produce two different configurations for one species, which is the thing
/// this module exists to prevent. Tests use it to show that 240 is early
/// termination of a converging sequence rather than convergence to something
/// else — a claim that cannot be made without varying the budget.
fn embed_with_budget(species: &CanonMol, universe: &Universe, iterations: usize) -> Embedding {
    let n = species.len();
    if n <= 1 {
        let pos = [[0.0f64; 3]; MAX_ATOMS];
        // One atom sits at the origin, which is already centred; zero atoms have
        // nowhere to sit. Both skip the solver rather than being special-cased
        // inside it.
        return Embedding { n, pos };
    }
    let target = targets(species, universe);

    let mut pos = initial_layout(&target, n);

    // SMACOF. Weight `1/d²` makes short (bonded) distances dominate, which keeps
    // local structure faithful while letting the global shape relax.
    //
    // **This row-wise `Σw(...)/Σw` update is the *diagonal approximation* to the
    // weighted Guttman transform `X' = V_w⁺ B(X) X`, not the transform**, and the
    // distinction is worth stating because it is invisible: with uniform weights
    // `V_w⁺` collapses to `(1/n)(I − J)` and the two coincide, which is the case
    // `experiments/src/embed.rs` implements. With `w = 1/d²` they differ, so
    // SMACOF's monotonicity *theorem* does not cover this iteration.
    //
    // Monotonicity here is therefore **measured, not proven** —
    // `stress_never_increases` is what stands behind it, and it is the guard to
    // run first if this is ever changed. The exact form was implemented and
    // measured against this one: both pass, so the ~80 lines of Cholesky solve
    // it needs buy a guarantee the corpus does not currently distinguish. If the
    // test ever fires, `X' = (V_w + J)⁻¹ B(X) X` by Cholesky is the known repair —
    // `V_w + J` is positive definite, so it needs no pivoting and therefore no
    // float comparison (§13.1 bans `max`, which a pivot search would want).
    //
    // **Jacobi, not Gauss-Seidel**, and that part is load-bearing regardless:
    // every `next[i]` is computed from the *previous* sweep's `pos`, so the
    // result does not depend on the order atoms are visited in. Reading
    // partially-updated coordinates would make the answer a function of the loop
    // order — a §13.1 hazard, and an invisible one, because a Gauss-Seidel sweep
    // also converges.
    let mut next = pos;
    for _ in 0..iterations {
        for i in 0..n {
            let mut num = [0.0f64; 3];
            let mut den = 0.0f64;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dij = target.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                if dij <= 0.0 {
                    continue;
                }
                let (Some(pi), Some(pj)) = (pos.get(i), pos.get(j)) else {
                    continue;
                };
                let diff = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                let raw = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
                #[expect(
                    clippy::neg_cmp_op_on_partial_ord,
                    reason = "the negation is the point: `!(raw > eps)` also skips a NaN \
                              separation, where `raw <= eps` would admit it and take the \
                              whole embedding to NaN. procedural_review records this exact \
                              rewrite being what traps NaN, not the finiteness test beside it"
                )]
                if !(raw > 1e-9) {
                    continue;
                }
                let wt = 1.0 / (dij * dij);
                for k in 0..3 {
                    if let (Some(acc), Some(&pjk), Some(&dk)) =
                        (num.get_mut(k), pj.get(k), diff.get(k))
                    {
                        *acc += wt * (pjk + dij * dk / raw);
                    }
                }
                den += wt;
            }
            if den > 0.0 {
                for k in 0..3 {
                    if let (Some(dst), Some(&v)) =
                        (next.get_mut(i).and_then(|p| p.get_mut(k)), num.get(k))
                    {
                        *dst = v / den;
                    }
                }
            }
        }
        pos = next;
    }

    // Centre on the centroid. Stress is translation-invariant so the optimiser
    // will not do this itself, and signatures are support functions measured from
    // the origin — without it every extent would carry an arbitrary offset.
    // Fixed summation order (§13.1).
    let mut centroid = [0.0f64; 3];
    for p in pos.iter().take(n) {
        for k in 0..3 {
            if let (Some(acc), Some(&v)) = (centroid.get_mut(k), p.get(k)) {
                *acc += v;
            }
        }
    }
    #[expect(
        clippy::cast_precision_loss,
        clippy::as_conversions,
        reason = "n <= MAX_ATOMS = 12, so the widening is exact"
    )]
    let inv_n = 1.0 / n as f64;
    for p in pos.iter_mut().take(n) {
        for k in 0..3 {
            if let (Some(v), Some(&species)) = (p.get_mut(k), centroid.get(k)) {
                *v -= species * inv_n;
            }
        }
    }

    Embedding { n, pos }
}

#[cfg(test)]
#[expect(
    clippy::indexing_slicing,
    clippy::unwrap_used,
    reason = "CLAUDE.md permits both inside #[cfg(test)] with a stated reason; every index \
              here is a loop bound over an atom count already established"
)]
mod tests {
    use super::*;
    use crate::canonical::canonicalise;
    use borbax_rng::{Domain, Stream};
    use borbax_universe::{
        BondOrder, ElementId, PeriodicTable, Universe, element::generate_elements,
    };

    fn table(seed: u64) -> PeriodicTable {
        generate_elements(seed)
    }

    /// Elements that can carry two bonds, so a chain is buildable.
    fn chain_capable(tbl: &PeriodicTable) -> Vec<ElementId> {
        (0..120usize)
            .filter_map(ElementId::from_index)
            .filter(|id| tbl.get(*id).is_some_and(|el| el.valence >= 2))
            .collect()
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

    fn canon(mol: &Mol12) -> CanonMol {
        canonicalise(mol)
            .unwrap_or_else(|err| unreachable!("fixture capped: {err}"))
            .0
    }

    #[test]
    fn embedding_is_deterministic() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let species = canon(&chain(7, &tbl, &ids));
        assert_eq!(embed(&species, &uni), embed(&species, &uni));
    }

    #[test]
    fn positions_are_finite_and_centred() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let emb = embed(&canon(&chain(9, &tbl, &ids)), &uni);
        assert!(
            emb.coords()
                .iter()
                .all(|point| point.iter().all(|coord| coord.is_finite()))
        );
        let mut sum = [0.0f64; 3];
        for perm in emb.coords() {
            for k in 0..3 {
                sum[k] += perm[k];
            }
        }
        assert!(
            sum.iter().all(|axis| axis.abs() < 1e-9),
            "not centred: {sum:?}"
        );
    }

    #[test]
    fn bonded_atoms_end_up_closer_than_distant_ones() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        // A chain's canonical order is not the build order, so this walks the
        // canonical molecule's own adjacency to find an end rather than
        // assuming atom 0 is one.
        let species = canon(&chain(8, &tbl, &ids));
        let emb = embed(&species, &uni);
        let end = (0..8u8)
            .find(|&i| species.mol().degree(i) == 1)
            .unwrap_or_else(|| unreachable!("a chain has two ends"));
        let hops = hop_distances(species.mol());
        let mut by_hop: Vec<(u8, usize)> = (0..8).map(|j| (hops[usize::from(end)][j], j)).collect();
        by_hop.sort_unstable();
        let dist_to = |j: usize| emb.distance(usize::from(end), j).unwrap().get();
        assert!(
            dist_to(by_hop[1].1) < dist_to(by_hop[4].1),
            "graph distance not reflected in geometry"
        );
        assert!(dist_to(by_hop[4].1) < dist_to(by_hop[7].1));
    }

    #[test]
    fn handles_degenerate_inputs() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        assert_eq!(embed(&canon(&Mol12::new()), &uni).len(), 0);
        assert!(embed(&canon(&Mol12::new()), &uni).is_empty());
        assert_eq!(embed(&canon(&chain(1, &tbl, &ids)), &uni).len(), 1);
        // Two atoms and a straight path are the other degenerate SMACOF inputs:
        // both are exactly realisable, and both were named in the Task 2 review.
        assert_eq!(embed(&canon(&chain(2, &tbl, &ids)), &uni).len(), 2);
        assert!(
            embed(&canon(&chain(2, &tbl, &ids)), &uni)
                .coords()
                .iter()
                .all(|point| point.iter().all(|coord| coord.is_finite()))
        );
    }

    /// **The Task 2 review's routed hazard.** `experiments`' `graph_distances`
    /// returns `f64::INFINITY` for an unreachable atom, which is safe there only
    /// because every fixture is a tree. `Mol12` admits disconnected molecules —
    /// two fragments in one graph — and an infinity in the target matrix takes
    /// stress to `inf - inf`, the whole embedding to NaN, and NaN then wins every
    /// downstream selection: a species with a garbage signature and a garbage
    /// affinity, with nothing failing anywhere.
    #[test]
    fn a_disconnected_molecule_embeds_finitely() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut mol = Mol12::new();
        for i in 0..6u8 {
            assert!(mol.add_atom(ids[usize::from(i) % ids.len()]).is_some());
        }
        // Two disjoint triangles-of-a-path: 0-1-2 and 3-4-5. Nothing joins them.
        for (a, b) in [(0, 1), (1, 2), (3, 4), (4, 5)] {
            assert!(mol.add_bond(a, b, BondOrder::SINGLE, &tbl).is_ok());
        }
        assert!(
            !mol.is_connected(),
            "the fixture must actually be disconnected"
        );

        let emb = embed(&canon(&mol), &uni);
        assert_eq!(emb.len(), 6);
        assert!(
            emb.coords()
                .iter()
                .all(|point| point.iter().all(|coord| coord.is_finite())),
            "a disconnected molecule produced a non-finite embedding: {:?}",
            emb.coords()
        );
        assert!(
            stress(&canon(&mol), &uni, &emb).is_finite(),
            "stress is not finite on a disconnected molecule"
        );
    }

    /// The property [`CanonMol`] was introduced for, asserted at the strength the
    /// Task 8 preamble specifies for canonical-position indexing: **bit-identical**,
    /// not merely equal up to a rigid motion.
    ///
    /// With `embed` taking a `CanonMol` this holds by construction and the test is
    /// proving the construction rather than the solver. That is the point — the
    /// caller-dependent call has no spelling. Probed by giving `embed` the raw
    /// `Mol12` instead: fails on the first relabelling.
    #[test]
    fn the_embedding_is_a_function_of_the_species_not_the_labelling() {
        let (tbl, uni) = (table(17), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(81, Domain::Molecule, 0);
        let mol = chain(7, &tbl, &ids);
        let want = embed(&canon(&mol), &uni);
        for _ in 0..200 {
            let mut perm: Vec<u8> = (0..7).collect();
            for i in (1..perm.len()).rev() {
                let j = usize::try_from(rng.next_range(u64::try_from(i + 1).unwrap())).unwrap();
                perm.swap(i, j);
            }
            let relabelled = mol
                .relabelled(&perm)
                .unwrap_or_else(|| unreachable!("perm permutes 0..7"));
            assert_eq!(
                embed(&canon(&relabelled), &uni),
                want,
                "relabelling changed the embedding"
            );
        }
    }

    /// **The one guarantee stress majorization actually makes.** If stress ever
    /// rises the Guttman transform is misderived — and a misderived transform
    /// still converges to *something*, so nothing else here would notice.
    ///
    /// Lifted from `experiments/src/embed.rs`, which is the half of that file the
    /// plan's Step 3 does not have. It is also what decides the weighting
    /// question: the row-wise weighted update is the diagonal approximation to
    /// `V_w⁺ B(X) X`, and the textbook monotonicity proof is for the exact form.
    #[test]
    fn stress_never_increases() {
        let (tbl, uni) = (table(17), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(82, Domain::Molecule, 0);
        for trial in 0..60u32 {
            // Random *molecules*, not only trees: rings are where the target
            // matrix stops being a tree metric, and they are exactly the graphs
            // the plan's own preamble flags as untested by `experiments`.
            let n = 3 + u8::try_from(trial % 10).unwrap();
            let mol = random_molecule(&mut rng, n, &tbl, &ids);
            let species = canon(&mol);
            if species.len() < 2 {
                continue;
            }
            let first = stress(&species, &uni, &embed_with_budget(&species, &uni, 1));
            let mut previous = f64::INFINITY;
            for step in 1..=ITERATIONS {
                let current = stress(&species, &uni, &embed_with_budget(&species, &uni, step));
                assert!(
                    current <= previous + 1e-12,
                    "trial {trial} step {step}: stress rose from {previous} to {current}"
                );
                previous = current;
            }
            // **"Never increases" is also true of an iteration that does
            // nothing**, so the descent has to be asserted too, or a transform
            // that returned its input unchanged would pass. Probed: returning
            // `pos` untouched fails here and not above.
            assert!(
                previous < first,
                "trial {trial}: stress never fell — {first} to {previous}"
            );
        }
    }

    /// What makes 240 a *budget choice* rather than a wrong answer: the residual
    /// keeps falling as the budget grows, so the fixed stop is early termination
    /// of a converging sequence and not convergence to something else.
    #[test]
    fn stress_falls_toward_zero_as_the_budget_grows() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        // A path is the one family whose graph distances are realisable exactly,
        // so the optimum is known to be zero and the limit is not a guess.
        let species = canon(&chain(6, &tbl, &ids));
        let budgets = [60usize, 240, 960, 3840];
        let stresses: Vec<f64> = budgets
            .iter()
            .map(|&n| stress(&species, &uni, &embed_with_budget(&species, &uni, n)))
            .collect();
        for w in stresses.windows(2) {
            assert!(
                w[0] > w[1] * 4.0,
                "quadrupling the budget only moved stress {} -> {}",
                w[0],
                w[1]
            );
        }
        assert!(
            stresses.last().is_some_and(|&s| s < 1e-5),
            "stress is not heading to zero: {stresses:?}"
        );
    }

    /// Bonded atoms must not land on top of each other, or the support function
    /// collapses and every signature reports the same shape.
    #[test]
    fn bonded_atoms_stay_apart() {
        let (tbl, uni) = (table(17), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(83, Domain::Molecule, 0);
        for _ in 0..20 {
            let species = canon(&random_tree(&mut rng, 10, &tbl, &ids));
            let emb = embed(&species, &uni);
            for a in 0..species.len() {
                for b in (a + 1)..species.len() {
                    let bonded = species
                        .mol()
                        .bond_order(u8::try_from(a).unwrap(), u8::try_from(b).unwrap())
                        .is_some();
                    if bonded {
                        let dist_to = emb.distance(a, b).unwrap().get();
                        assert!(dist_to > 0.3, "bonded atoms {a},{b} only {dist_to} apart");
                    }
                }
            }
        }
    }

    /// A random molecule: a spanning tree plus a few extra edges, so rings and
    /// branches appear rather than trees only.
    fn random_molecule(rng: &mut Stream, n: u8, tbl: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = random_tree(rng, n, tbl, ids);
        let size = u8::try_from(mol.len()).unwrap();
        for _ in 0..rng.next_range(4) {
            let a = u8::try_from(rng.next_range(u64::from(size))).unwrap();
            let b = u8::try_from(rng.next_range(u64::from(size))).unwrap();
            // A refused bond is not a mutation; the fixture is whatever survives.
            let _ = mol.add_bond(a, b, BondOrder::SINGLE, tbl);
        }
        mol
    }

    /// A random tree, built by attaching each new atom to an existing one.
    fn random_tree(rng: &mut Stream, n: u8, tbl: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = Mol12::new();
        assert!(mol.add_atom(ids[0]).is_some());
        for i in 1..n {
            let pick = usize::try_from(rng.next_range(u64::try_from(ids.len()).unwrap())).unwrap();
            let Some(child) = mol.add_atom(ids[pick]) else {
                break;
            };
            let parent = u8::try_from(rng.next_range(u64::from(i))).unwrap();
            // A refused bond leaves the molecule alone, which would give a
            // disconnected fixture — fine for `embed`, but this helper promises a
            // tree, so it retries down the orders rather than shrugging.
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
}
