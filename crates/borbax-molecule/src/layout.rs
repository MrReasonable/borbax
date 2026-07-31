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
//! The V0 plan's Step 3 ends by rotating the result into a frame built from the
//! lowest-indexed non-degenerate atoms. That is removed. An atom-anchored frame
//! is itself a canonical-order discontinuity — the very thing it was meant to
//! cure — and `experiments/`'s G2 harness scores it worst of four descriptors,
//! below doing nothing at all, on every seed in both regimes.
//!
//! **Its doc comment's *mechanism* was right and only its remedy was wrong, and
//! conflating those cost a round.** It argued that canonicalising over the 60
//! icosahedral rotations "cannot absorb" frame drift because a generic SO(3)
//! rotation sits up to ~44° from the nearest group element. True. Deleting the
//! bad remedy did not address it; what addresses it is `initial_layout` being a
//! continuous function of the graph, which is where that argument now lives.
//!
//! **Quote G2's numbers with their provenance or not at all.** An earlier version
//! of this doc gave `D_frame 0.9195` against `D_raw 0.9560` as though they
//! described this code. They do not: they were measured on `experiments`' own
//! embedder, which differs in initialisation, weighting and radius source, over a
//! corpus of **8–20 atoms** when `MAX_ATOMS` is 12 — so more than half of it is
//! unbuildable here — and containing **no rings**, since its only generator is
//! `random_tree`. Re-measured at 4–12 atoms the ordering of `D_frame` against
//! `D_raw` is unchanged and unanimous across five seeds (0.7360–0.7530 against
//! 0.8147–0.8310, size-matched), so **the deletion is robust to the corpus**. The
//! *margins* are not, and they are not this crate's numbers.
//!
//! What this crate measures on its own embedder is in
//! `a_one_atom_edit_moves_the_shape_less_than_an_unrelated_molecule_does`.
//!
//! Two independent literatures reach the same conclusion, which matters because
//! it means the decision does not rest on one in-house harness: Kendall's shape
//! space (Bull. LMS 16, 1984) works in the **quotient** by the similarity group
//! precisely because any representative-choosing map is discontinuous wherever
//! the choice is degenerate; and 3D shape retrieval displaced PCA-based canonical
//! alignment with rotation-invariant descriptors for the same reason.
//!
//! **Do not reintroduce a frame, and do not replace it with any other
//! representative-choosing step**, which includes storing a lexicographically
//! smallest rotation — measured worse still than the frame. See Task 9's routed
//! requirement in the plan.

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
/// **The `diameter + 1` fill decides a geometry, and only its finiteness is
/// justified.** With two isolated atoms the only finite entries are the diagonal,
/// so `diameter = 0` and the fill is `1` — the same value a bond gets, placing two
/// unbonded atoms at exactly bonded distance. More generally the inter-fragment
/// target depends on the *fragments' internal diameter*, so the same two fragments
/// sit at different separations depending on how long the other one is. Unreachable
/// in V0 (nothing interns a disconnected species; complexes are a species pair with
/// an affinity), and routed to whoever implements Cleave at Task 15: either `embed`
/// refuses a disconnected molecule, or the separation is chosen with a physical
/// reason rather than inherited from a NaN guard.
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
///
/// **`h * (r_i + r_j)` is pinned as written.** Distributing it to
/// `h*r_i + h*r_j` is algebraically identical, reads as a tidy-up, and differs in
/// the last bits on **22.19%** of 74.2M sampled pairs over the measured operating
/// range — moving every coordinate of every embedding downstream with the whole
/// suite green. Measured by mutation: a corpus digest goes
/// `0xb893dc3a6ceb12a0` → `0x9027e1b01b559eb4`, 12 tests passing either way.
/// There is no golden yet, which is exactly why the shape is written down now.
fn targets(species: &CanonMol, universe: &Universe) -> [[f64; MAX_ATOMS]; MAX_ATOMS] {
    let mol = species.mol();
    let n = mol.len();
    let hops = hop_distances(mol);
    // **A missing element falls back to 1.0, which is inside the measured real
    // range of [0.3, 2.04] — so a mismatched `(CanonMol, Universe)` pair yields a
    // *plausible* shape rather than a failure.** `Mol12::mass` faces the identical
    // condition and returns `Option` instead, arguing that a silent zero would put
    // a wrong mass into a conservation check. The two decisions disagree, and this
    // one is deliberate: `embed` is infallible by design and the coupling is
    // unchecked crate-wide, so making this one site fallible would buy a `Result`
    // in every caller without closing the class. Recorded rather than fixed, and
    // the boundary to close is the pairing itself, not this lookup.
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
/// obligation to fall with it. The unweighted sum is not exported at all — see
/// the note below on why a same-signature twin was the wrong way to keep it.
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

// **`raw_stress` was removed rather than kept as a foil.** It returned the
// unweighted mismatch `Σ(δ−d)²` from a signature identical to `stress`'s, so the
// two were substitutable at any call site with nothing but the identifier to tell
// them apart — and confusing them is exactly the near-miss recorded above, which
// nearly bought an 80-line rewrite of a correct algorithm. It had zero callers. A
// documented footgun with no consumer is not an API; if Task 20's sweep wants the
// interpretable residual, give it a distinct return type rather than a twin.

/// The starting configuration: each atom placed at its target distances to three
/// landmark atoms.
///
/// **This must be a continuous function of the target matrix, not of the atom
/// index, and that is the single most consequential line in this file.** The
/// predecessor was a fixed table of twelve rational offsets indexed by canonical
/// position. Canonical order is *not* continuous under a graph edit — adding one
/// atom can permute every rank — so a parent and its one-atom mutant started from
/// unrelated configurations, and SMACOF converged to the same shape in an
/// unrelated **orientation**.
///
/// Measured on this embedder, size-matched, concordance of a rotation-sensitive
/// descriptor (1.0 = perfect locality, 0.5 = none):
///
/// ```text
///                    D_raw     D_group
/// index table       0.5104     0.5849      <- no locality at all
/// landmark start    0.7955     0.7595
/// ```
///
/// The weighting is irrelevant to this (0.5104 against 0.5166 with uniform
/// weights); the initialisation is the entire effect. And the damage was purely
/// orientational — a rotation-*invariant* shape metric was unchanged (1.365
/// against 1.361), and mean stress at 240 was 0.35572 against 0.35533. Both
/// starts found the same shapes, to the same quality, in different frames.
///
/// Why an arbitrary frame is fatal rather than untidy: §8.3's `affinity`
/// maximises over 60 discrete rotations, which sample SO(3) at roughly 44°
/// spacing and **cannot absorb** an arbitrary reorientation. So a one-atom
/// mutation moved the parent's best binding pose to an unrelated rotation index
/// and affinity jumped rather than drifted — neutral networks never form,
/// mutation is a random walk in binding space, and nothing accumulates, with
/// every test green.
///
/// **The deleted `canonicalise_frame` described this mechanism correctly and
/// proposed the wrong remedy.** An atom-anchored frame is itself a
/// canonical-order discontinuity, which is why G2 measures it as worse than no
/// alignment. Removing the bad remedy did not address the mechanism; this does.
///
/// **What carries the property is the *coordinates*, not the pivot choice**, and
/// probing separated them. Replacing `farthest_from` with fixed indices 1 and 2
/// leaves the locality test green: an atom's coordinates are still its
/// graph-distance vector, which is the intrinsic, continuous quantity. Restoring
/// the offset table fails it at exactly 0.6041. So the landmark *selection* is a
/// quality refinement and the graph-distance *coordinates* are the fix — do not
/// read the `farthest_from` call as the load-bearing part.
///
/// Landmarks are picked from the graph metric, so the construction stays a pure
/// function of the species and introduces no randomness and no transcendental —
/// `experiments/src/embed.rs` adds a `Stream` jitter here and this deliberately
/// does not, because that would make `borbax-rng` a runtime dependency in a
/// result-affecting path.
fn initial_layout(target: &[[f64; MAX_ATOMS]; MAX_ATOMS], n: usize) -> [[f64; 3]; MAX_ATOMS] {
    // Argmax by explicit comparison with a lowest-index tie-break, never
    // `max_by`/`total_cmp`: §13.4 disallows both, and an index tie-break is what
    // makes the pick deterministic when two atoms are equidistant — which is the
    // common case in the symmetric graphs small molecules actually are.
    // `targets_between_distinct_atoms_are_strictly_positive` establishes the
    // finiteness this comparison relies on.
    let farthest_from = |from: usize, second: Option<usize>| -> usize {
        let row = |a: usize| -> f64 {
            let direct = target
                .get(from)
                .and_then(|r| r.get(a))
                .copied()
                .unwrap_or(0.0);
            direct
                + second.map_or(0.0, |t| {
                    target.get(t).and_then(|r| r.get(a)).copied().unwrap_or(0.0)
                })
        };
        let mut best = 0usize;
        let mut best_val = row(0);
        for a in 1..n {
            let val = row(a);
            if val > best_val {
                best_val = val;
                best = a;
            }
        }
        best
    };

    // Three landmarks, each the atom furthest from those already chosen. The
    // first is atom 0 — and unlike the old table this is not load-bearing, since
    // what an atom's coordinates *are* is its graph-distance vector, an intrinsic
    // quantity, rather than an offset assigned to its rank.
    let p1 = farthest_from(0, None);
    let p2 = farthest_from(0, Some(p1));

    let mut pos = [[0.0f64; 3]; MAX_ATOMS];
    for i in 0..n {
        for (k, pivot) in [0usize, p1, p2].into_iter().enumerate() {
            if let (Some(dst), Some(&val)) = (
                pos.get_mut(i).and_then(|p| p.get_mut(k)),
                target.get(pivot).and_then(|r| r.get(i)),
            ) {
                *dst = val;
            }
        }
    }
    pos
}

/// Cholesky factor of `V_w + J`, the matrix the exact Guttman step solves against.
#[expect(
    clippy::indexing_slicing,
    reason = "every index is a loop bound over n <= MAX_ATOMS, the array's length"
)]
fn cholesky(w: &[[f64; MAX_ATOMS]; MAX_ATOMS], n: usize) -> Option<[[f64; MAX_ATOMS]; MAX_ATOMS]> {
    #[expect(
        clippy::cast_precision_loss,
        clippy::as_conversions,
        reason = "n <= MAX_ATOMS = 12, so the widening is exact"
    )]
    let inv_n = 1.0 / n as f64;

    // `mat = V_w + J`. `V_w` has zero row sums by construction, and `J = 11ᵀ/n`
    // adds `1/n` everywhere — the regularisation that makes the singular `V_w`
    // invertible while leaving `V_w⁺ = (V_w + J)⁻¹ − J`.
    let mut mat = [[0.0f64; MAX_ATOMS]; MAX_ATOMS];
    for row in 0..n {
        let mut diag = 0.0f64;
        for col in 0..n {
            if row != col {
                diag += w[row][col];
                mat[row][col] = -w[row][col] + inv_n;
            }
        }
        mat[row][row] = diag + inv_n;
    }

    let mut lower = [[0.0f64; MAX_ATOMS]; MAX_ATOMS];
    for col in 0..n {
        let mut pivot = mat[col][col];
        for &done in lower[col].iter().take(col) {
            pivot -= done * done;
        }
        // A non-positive pivot means the matrix was not positive definite, which
        // cannot happen for `V_w + J` with non-negative weights. Negated so a NaN
        // takes this branch rather than sliding into `sqrt`.
        #[expect(
            clippy::neg_cmp_op_on_partial_ord,
            reason = "a NaN pivot must fail the factorisation, not pass it"
        )]
        if !(pivot > 0.0) {
            return None;
        }
        let diag = pivot.sqrt();
        lower[col][col] = diag;
        for row in (col + 1)..n {
            let mut acc = mat[row][col];
            for (&a, &b) in lower[row].iter().zip(&lower[col]).take(col) {
                acc -= a * b;
            }
            lower[row][col] = acc / diag;
        }
    }
    Some(lower)
}

/// Solve `L Lᵀ x = b` by forward then back substitution.
#[expect(
    clippy::indexing_slicing,
    reason = "every index is a loop bound over n <= MAX_ATOMS, the array's length"
)]
fn solve(
    lower: &[[f64; MAX_ATOMS]; MAX_ATOMS],
    rhs: &[f64; MAX_ATOMS],
    n: usize,
) -> [f64; MAX_ATOMS] {
    let mut fwd = [0.0f64; MAX_ATOMS];
    for row in 0..n {
        let mut acc = rhs[row];
        for prev in 0..row {
            acc -= lower[row][prev] * fwd[prev];
        }
        fwd[row] = acc / lower[row][row];
    }
    let mut out = [0.0f64; MAX_ATOMS];
    for row in (0..n).rev() {
        let mut acc = fwd[row];
        for later in (row + 1)..n {
            acc -= lower[later][row] * out[later];
        }
        out[row] = acc / lower[row][row];
    }
    out
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

    // SMACOF, with the **exact** weighted Guttman transform `X' = V_w⁺ B(X) X`.
    //
    // Weight `1/d²` is **Kamada–Kawai's** choice (Inf. Process. Lett. 31(1), 1989,
    // spring constant `K/d²`), endorsed by GKN §2 — "we always take
    // `w_ij = d_ij^-2`". Not Sammon's, which is `1/d`; an earlier note here
    // mis-attributed it. The reason it is right for a *graph* target is that each
    // term becomes a squared **relative** error, so the objective is invariant
    // under uniform rescaling of the target matrix — and `targets` scales with the
    // drawn element-radius series, so that invariance is protecting something real.
    //
    // **What it costs is worth knowing before Task 9 reads this geometry.** `1/d²`
    // down-weights a 5-hop pair to 1/25 of a 1-hop pair, so long distances are
    // fitted loosely by construction — while the signature is a support function,
    // an outer-hull quantity dominated by the *largest* extents. The solver and its
    // consumer are weighted at opposite ends of the distance range. That is the
    // explanation for the measured R_max error of −2.03% at n = 12 growing with n
    // while the weighted stress looks converged, and it means weighted stress is
    // the wrong instrument for choosing `ITERATIONS`.
    //
    // **The V0 plan's row-wise `Σw(...)/Σw` form is not this, and two earlier
    // versions of this comment misdescribed it.** It is a published algorithm:
    // Gansner, Koren & North, *Graph Drawing by Stress Majorization*, GD 2004,
    // §2.3 eq. (12), which they call **localized optimization**. It is neither the
    // diagonal approximation `D⁻¹B(X)X` (those differ by O(1)) nor exact under
    // uniform weights (the first draft of this comment claimed both). Measured, it
    // is one **Jacobi sweep** of `V X' = B(X)X` — and equivalently the
    // **over-relaxed** Guttman transform `X + α(Γ(X) − X)` with `α = n/(n−1)`,
    // reproduced to 7.3e-15. That α is the whole story: 2.00 at n = 2, 1.50 at
    // n = 3, 1.09 at n = 12. The iteration is most over-relaxed exactly where the
    // molecules are smallest.
    //
    // GKN do prove it monotone — "guaranteed to strictly decrease the stress...
    // oscillations and non-convergence are impossible" — but **for the sequential
    // sweep**, moving one node at a time against already-updated neighbours. We run
    // Jacobi deliberately, for the §13.1 loop-order reason below, and the published
    // theorem does not cover that. (A separate argument does: with
    // `s = D⁻¹(VX − B(X)X)`, `τ(X';X) − τ(X;X) = −sᵀ(D + W)s`, and `D + W` is the
    // signless Laplacian, PSD for any non-negative weights — Householder–John
    // applied to a signless Laplacian, standard parts in a combination not found in
    // the MDS literature.)
    //
    // So it is monotone and shares this transform's stationary points — **and it
    // still fails at n = 2, because monotone does not mean convergent.** `D + W` is
    // singular exactly when the weight graph has a bipartite component; the weight
    // graph here is `K_n` (every target is strictly positive), and `K_n` is
    // bipartite iff n ≤ 2. At α = 2 the update `2Γ(X) − X` is a pure reflection
    // about the exact answer: an exact period-2 cycle at *constant* stress, which
    // `stress_never_increases` cannot see. A dimer then reports 92% too long at an
    // even budget and 92% too short at an odd one, making `ITERATIONS` being even
    // silently load-bearing. Measured, target 1.15625:
    //
    //     row-wise   budget 239 → 0.088633   budget 240 → 2.223867   (forever)
    //     this form  every budget → 1.156250, stress 0.000000000
    //
    // A dimer is the first product of every condensation, so that is not a corner
    // case. This form is also 6.7× better converged at n = 12 (1.18e-3 against
    // 7.98e-3 at budget 240). `a_dimer_lands_at_the_sum_of_its_radii` is the guard.
    //
    // **This form also closes a second defect the row-wise one carries**, and the
    // distinction is easy to lose: at a coincident pair the row-wise update drops
    // `w_ij` from its *denominator* as well as from the numerator, which replaces
    // `diag(V)` with something strictly smaller and breaks the descent argument
    // (measured: stress rose in 10 of 3000 artificially collapsed starts, worst
    // 6.60 → 27.00 in one step). Here `V` is built from the weights alone, before
    // the loop, and `cholesky` sees all of them — the `raw > 1e-9` skip below
    // touches only the `B(X)` term, which is the only part that depends on `X`.
    // That is the standard convention and it is correct by construction rather
    // than by luck.
    //
    // `B(X)X` has zero column sums, so `J·(B(X)X) = 0` and
    // `V_w⁺ = (V_w + J)⁻¹ − J` reduces to a solve against the Cholesky factor.
    // That factor depends only on the weights, hence only on the target matrix, so
    // it is built once here rather than once per iteration — which is also the
    // loop-invariant hoist the performance lane measured at ~10%.
    //
    // **Jacobi, not Gauss-Seidel**, load-bearing independently: `B(X)` is built
    // entirely from the previous sweep's `pos`, so the result does not depend on
    // the order atoms are visited in. Reading partially-updated coordinates would
    // make the answer a function of the loop order — a §13.1 hazard, and an
    // invisible one, because a Gauss-Seidel sweep also converges.
    let mut w = [[0.0f64; MAX_ATOMS]; MAX_ATOMS];
    for i in 0..n {
        for j in 0..n {
            let dij = target.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
            if i != j
                && dij > 0.0
                && let Some(cell) = w.get_mut(i).and_then(|r| r.get_mut(j))
            {
                *cell = 1.0 / (dij * dij);
            }
        }
    }
    let Some(chol) = cholesky(&w, n) else {
        return Embedding { n, pos };
    };

    for _ in 0..iterations {
        let mut rhs = [[0.0f64; MAX_ATOMS]; 3];
        for i in 0..n {
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
                // **Not a NaN defence for `embed`, and an earlier version of this
                // comment implied it was.** The negation does skip a NaN where
                // `raw <= eps` would admit it — but the route a non-finite value
                // actually takes is upstream: `initial_layout` reads target cells,
                // so one poisoned cell makes the whole start non-finite, this guard
                // then fires on every pair, `den` stays zero, and `embed` returns
                // the poisoned start verbatim. It does not trap the NaN; it turns
                // the solver into a no-op that launders it. Measured by poisoning.
                //
                // No finiteness check is added, deliberately: over 300 universe
                // seeds and 26,906 elements `Element::radius` is finite and lies in
                // [0.3, 2.04], so targets are provably finite and the condition is
                // unreachable. `targets_between_distinct_atoms_are_strictly_positive`
                // pins that, and adding a guard for an unreachable condition is the
                // wrong trade.
                #[expect(
                    clippy::neg_cmp_op_on_partial_ord,
                    reason = "the negation skips a coincident-or-NaN pair where `raw <= eps` \
                              would admit a NaN. `f64::max` is banned for §13.1 in \
                              clippy.toml's disallowed-methods block, which is the \
                              authority for this spelling"
                )]
                if !(raw > 1e-9) {
                    continue;
                }
                let wij = w.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                let scale = wij * dij / raw;
                for (k, row) in rhs.iter_mut().enumerate() {
                    if let (Some(acc), Some(&dk)) = (row.get_mut(i), diff.get(k)) {
                        *acc += scale * dk;
                    }
                }
            }
        }
        for (k, row) in rhs.iter().enumerate() {
            let solved = solve(&chol, row, n);
            for i in 0..n {
                if let (Some(dst), Some(&v)) =
                    (pos.get_mut(i).and_then(|p| p.get_mut(k)), solved.get(i))
                {
                    *dst = v;
                }
            }
        }
    }

    // Centre on the centroid. Stress is translation-invariant so the optimiser
    // will not do this itself, and signatures are support functions measured from
    // the origin — without it every extent would carry an arbitrary offset.
    //
    // **Two shapes here are pinned and both have a plausible rewrite that moves
    // every number** (§13.1): the `i`-then-`k` accumulation order, and the
    // reciprocal-multiply `c * inv_n` rather than `c / n`. The second differs in
    // the last bits on **21.54%** of 11.8M samples — digest
    // `0xb893dc3a6ceb12a0` → `0x4421fc86867441a5`, suite green either way.
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
            if let (Some(coord), Some(&centre)) = (p.get_mut(k), centroid.get(k)) {
                *coord -= centre * inv_n;
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

    /// **The property the whole design rests on (§8.2), pinned rather than
    /// asserted.** A one-atom edit must move the shape less than an unrelated
    /// molecule of the same size does, or mutation is a random walk and nothing
    /// accumulates.
    ///
    /// Concordance is the fraction of trials where the mutant is nearer the
    /// parent than the stranger is: 1.0 perfect locality, **0.5 none**. Measured
    /// on this embedder with everything but the initialisation held fixed:
    ///
    /// ```text
    /// index-keyed START table   0.6041   (801/1326)
    /// landmark start (shipped)  0.8303   (1101/1326)
    /// ```
    ///
    /// The descriptor is deliberately **rotation-sensitive** — a raw support
    /// function over 42 fixed directions. A rotation-invariant descriptor scores
    /// the two starts identically (measured elsewhere at 1.365 against 1.361),
    /// because both find the same *shapes*; what the index-keyed table destroyed
    /// was the *frame*. §8.3's `affinity` maximises over 60 discrete rotations at
    /// ~44° spacing and cannot absorb an arbitrary reorientation, so a frame that
    /// jumps under a one-atom edit is fatal rather than untidy.
    ///
    /// This is not Task 9's test. Task 9 will measure locality of the *signature*;
    /// if that fails and this passes, the signature is the culprit. Conflating
    /// them was the reason the first version of this task routed the property
    /// downstream and asserted nothing here.
    #[test]
    fn a_one_atom_edit_moves_the_shape_less_than_an_unrelated_molecule_does() {
        let (tbl, uni) = (table(17), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(777, Domain::Molecule, 0);
        let dirs = crate::geodesic::Geodesic::<42>::build()
            .unwrap_or_else(|_| unreachable!("D=42 builds"));
        let sig = |sp: &CanonMol| -> Vec<f64> {
            let emb = embed(sp, &uni);
            (0..42)
                .map(|d| {
                    let dir = dirs.dirs()[d];
                    let mut best = f64::NEG_INFINITY;
                    for p in emb.coords() {
                        let v = p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2];
                        if v > best {
                            best = v;
                        }
                    }
                    best
                })
                .collect()
        };
        let apart = |a: &[f64], b: &[f64]| -> f64 {
            a.iter()
                .zip(b)
                .map(|(x, y)| (x - y) * (x - y))
                .sum::<f64>()
                .sqrt()
        };

        let (mut nearer, mut trials) = (0u32, 0u32);
        for _ in 0..2000 {
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
            // Size-matched: the stranger has the mutant's atom count, so the
            // comparison is about shape rather than about size.
            let stranger =
                random_molecule(&mut rng, u8::try_from(mutant.len()).unwrap(), &tbl, &ids);
            if stranger.len() != mutant.len() {
                continue;
            }
            let (cp, cm, cs) = (canon(&parent), canon(&mutant), canon(&stranger));
            if cp.len() < 3 {
                continue;
            }
            let base = sig(&cp);
            if apart(&base, &sig(&cm)) < apart(&base, &sig(&cs)) {
                nearer += 1;
            }
            trials += 1;
        }

        assert!(trials > 1000, "only {trials} usable trials");
        let concordance = f64::from(nearer) / f64::from(trials);
        assert!(
            concordance > 0.75,
            "locality lost: concordance {concordance:.4} ({nearer}/{trials}); \
             0.5 is none, the index-keyed start scored 0.6041, this must clear 0.75"
        );
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
        for point in emb.coords() {
            for k in 0..3 {
                sum[k] += point[k];
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

        // **Finiteness alone does not pin the `diameter + 1` fill**, and that is
        // the whole point of the fill. Probed: replacing it with `0` leaves every
        // assertion above green, because a zero target is simply skipped by both
        // the weight loop and the Guttman step — the two fragments stop
        // constraining each other and drift to wherever the start put them.
        // Finite, deterministic, and physically wrong.
        //
        // What the fill buys is that unreachable pairs get the *largest* target
        // in the matrix, so fragments are pushed apart. That is observable: no
        // bonded pair may end up further apart than the closest cross-fragment
        // pair.
        let species = canon(&mol);
        let reach = |from: usize| -> u16 {
            let mut seen = 1u16 << u32::try_from(from).unwrap();
            let mut frontier = seen;
            while frontier != 0 {
                let mut next = 0u16;
                for node in 0..species.len() {
                    if frontier & (1u16 << u32::try_from(node).unwrap()) != 0 {
                        next |= species.mol().neighbours(u8::try_from(node).unwrap());
                    }
                }
                next &= !seen;
                seen |= next;
                frontier = next;
            }
            seen
        };
        let (mut longest_bond, mut closest_cross) = (0.0f64, f64::INFINITY);
        let mut cross_pairs = 0u32;
        for i in 0..species.len() {
            let component = reach(i);
            for j in (i + 1)..species.len() {
                let sep = emb.distance(i, j).unwrap().get();
                if component & (1u16 << u32::try_from(j).unwrap()) == 0 {
                    cross_pairs += 1;
                    if sep < closest_cross {
                        closest_cross = sep;
                    }
                } else if species
                    .mol()
                    .bond_order(u8::try_from(i).unwrap(), u8::try_from(j).unwrap())
                    .is_some()
                    && sep > longest_bond
                {
                    longest_bond = sep;
                }
            }
        }
        assert!(cross_pairs > 0, "the fixture has no cross-fragment pairs");
        assert!(
            closest_cross > longest_bond,
            "fragments are not pushed apart: closest cross-fragment pair {closest_cross} \
             is nearer than the longest bond {longest_bond}"
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

    /// **The guard for the defect the row-wise form shipped with.** A dimer's
    /// optimum is exact and known — two atoms at the sum of their radii — so this
    /// needs no epsilon fitted to whatever the solver happens to produce.
    ///
    /// The V0 plan's row-wise update reduces at n = 2 to `δ → |2d − δ|`, a
    /// reflection about the target with no contraction: an exact period-2 cycle at
    /// *constant* stress. `stress_never_increases` cannot see it (constant is
    /// non-increasing), and `handles_degenerate_inputs` built `chain(2)` — this
    /// exact fixture — while asserting only `len() == 2` and finiteness. So the
    /// whole suite was green with every dimer 92% too long.
    ///
    /// Asserted at **both parities** of the budget, because the defect's signature
    /// is that the answer depends on whether `ITERATIONS` is even.
    #[test]
    fn a_dimer_lands_at_the_sum_of_its_radii() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let species = canon(&chain(2, &tbl, &ids));
        let want = targets(&species, &uni)[0][1];
        assert!(want > 0.0, "the fixture has no target to hit");

        for budget in [1usize, 2, 3, 239, ITERATIONS, ITERATIONS + 1] {
            let emb = embed_with_budget(&species, &uni, budget);
            let got = emb.distance(0, 1).unwrap_or(Span::ZERO).get();
            assert!(
                (got - want).abs() < 1e-9,
                "budget {budget}: dimer at {got}, wants {want}"
            );
        }
        // Stationary, not merely non-increasing — a 2-cycle is non-increasing.
        let at = |b: usize| stress(&species, &uni, &embed_with_budget(&species, &uni, b));
        assert!(
            (at(ITERATIONS) - at(ITERATIONS + 1)).abs() < 1e-12,
            "the iteration is not stationary at the shipped budget"
        );
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
            // From 2, not 3: the period-2 cycle the exact transform fixes lives
            // at exactly n = 2, and the predecessor's corpus started at 3.
            let n = 2 + u8::try_from(trial % 11).unwrap();
            let mol = random_molecule(&mut rng, n, &tbl, &ids);
            let species = canon(&mol);
            if species.len() < 2 {
                continue;
            }
            // Budget 0 — the un-iterated start — not budget 1. With the exact
            // transform a dimer reaches its optimum in a single step, so
            // `stress(1) == stress(240)` and a descent assertion anchored at 1
            // fails on a converged fixture rather than on a broken solver.
            let first = stress(&species, &uni, &embed_with_budget(&species, &uni, 0));
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
    ///
    /// **Measured as excess above a long-budget reference, and the predecessor's
    /// assumption of a zero optimum was a fixture accident.** It used a path,
    /// documented as "the one family whose graph distances are realisable
    /// exactly". That is false in general: with `δ_ij = h_ij(r_i + r_j)` a path
    /// is exactly realisable on a line only if the radii are in arithmetic
    /// progression along the chain. The shipped fixture satisfied it by
    /// coincidence — `chain` cycles consecutive element ids and `radius` is
    /// linear in index within a shell, giving step 0.03125 exactly. Build the
    /// same chain from non-consecutive ids and the optimum is not zero, so the
    /// old assertion failed with a message naming the *solver* when the
    /// *fixture* had changed — `bonds.rs:1341`'s failure mode exactly.
    ///
    /// A star does not embed exactly, so its limit is nonzero and unknown, which
    /// is the point: measuring the excess needs no assumption about the optimum.
    #[test]
    fn stress_falls_toward_its_limit_as_the_budget_grows() {
        let (tbl, uni) = (table(4), Universe::generate(4));
        let ids = chain_capable(&tbl);
        // The centre must actually be able to carry the spokes, so it is the
        // highest-valence element in the table rather than whatever `ids[0]` is.
        // The predecessor of this fixture used `let _ = add_bond(..)` and
        // silently built a partial star.
        let hub = ids
            .iter()
            .copied()
            .max_by_key(|id| tbl.get(*id).map_or(0, |el| el.valence))
            .unwrap_or_else(|| unreachable!("the table has chain-capable elements"));
        let spokes = u8::try_from(
            tbl.get(hub)
                .map_or(2, |el| usize::from(el.valence))
                .min(MAX_ATOMS - 1),
        )
        .unwrap_or(2);
        assert!(
            spokes >= 4,
            "a {spokes}-spoke star is too small to be a hub"
        );

        let mut star = Mol12::new();
        assert!(star.add_atom(hub).is_some());
        for i in 0..spokes {
            assert!(star.add_atom(ids[usize::from(i) % ids.len()]).is_some());
            assert!(star.add_bond(0, i + 1, BondOrder::SINGLE, &tbl).is_ok());
        }
        let species = canon(&star);

        let limit = stress(&species, &uni, &embed_with_budget(&species, &uni, 20_000));
        let excess: Vec<f64> = [60usize, 240, 960, 3840]
            .iter()
            .map(|&b| stress(&species, &uni, &embed_with_budget(&species, &uni, b)) - limit)
            .collect();

        // Without this the test passes on a fixture that was already converged
        // at the shortest budget — which is what the landmark start does to a
        // path, and is how the predecessor became vacuous.
        assert!(
            excess[0] > 1e-4,
            "fixture is already converged at budget 60; it tests nothing: {excess:?}"
        );
        for w in excess.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-12,
                "excess above the limit rose: {excess:?}"
            );
        }
        assert!(
            excess[3] < excess[0] * 1e-6,
            "excess is not heading to zero: {excess:?}"
        );
    }

    /// Bonded atoms must not land on top of each other, or the support function
    /// collapses and every signature reports the same shape.
    #[test]
    fn bonded_atoms_stay_apart() {
        let (tbl, uni) = (table(17), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(83, Domain::Molecule, 0);
        // **This test iterates only bonded pairs, so a fixture with no bonds
        // passes it vacuously.** Counting what was examined is what stops a
        // generator regression from silently emptying it.
        let mut examined = 0u32;
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
                        examined += 1;
                        let dist_to = emb.distance(a, b).unwrap().get();
                        assert!(dist_to > 0.3, "bonded atoms {a},{b} only {dist_to} apart");
                    }
                }
            }
        }
        assert!(examined > 100, "only {examined} bonded pairs examined");
    }

    /// **The precondition the whole weighting rests on**, and nothing in the
    /// types states it: `w_ij = 1/d_ij²` is an infinity at a zero target, and
    /// `stress` divides by `d_ij²` outright.
    ///
    /// It holds because `hop_distances` returns at least 1 for `i != j` — a
    /// BFS depth, or `diameter + 1` for an unreachable pair — and radii are
    /// drawn strictly positive. Both halves are load-bearing and neither is
    /// local to this file, which is why this is asserted rather than commented.
    ///
    /// This also settles the asymmetry a reviewer flagged between `stress` and
    /// `raw_stress`: `stress`'s `tgt <= 0.0` guard exists for the division, and
    /// `raw_stress` deliberately has none, because if a target ever *were* zero
    /// the unweighted mismatch `δ − 0` is a real mismatch that should be
    /// counted, not skipped. Adding the guard there would silently discard it.
    #[test]
    fn targets_between_distinct_atoms_are_strictly_positive() {
        let (tbl, uni) = (table(17), Universe::generate(4));
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(84, Domain::Molecule, 0);
        let mut smallest = f64::INFINITY;
        let mut disconnected_seen = 0u32;
        for trial in 0..200u32 {
            let n = 2 + u8::try_from(trial % 11).unwrap();
            // Every third fixture is deliberately fragmented. Without this the
            // corpus is all-connected — `random_molecule` builds a spanning tree
            // first — so the `diameter + 1` fill is never reached and this test
            // says nothing about the case it exists for. Measured: it did not.
            let mol = if trial % 3 == 0 {
                disconnected_molecule(&mut rng, n, &tbl, &ids)
            } else {
                random_molecule(&mut rng, n, &tbl, &ids)
            };
            if !mol.is_connected() {
                disconnected_seen += 1;
            }
            let species = canon(&mol);
            let target = targets(&species, &uni);
            for (i, row) in target.iter().enumerate().take(species.len()) {
                for (j, &tgt) in row.iter().enumerate().take(species.len()) {
                    if i == j {
                        continue;
                    }
                    assert!(tgt > 0.0, "target[{i}][{j}] = {tgt} on a {n}-atom molecule");
                    assert!(tgt.is_finite(), "target[{i}][{j}] is not finite");
                    if tgt < smallest {
                        smallest = tgt;
                    }
                }
            }
        }
        assert!(smallest.is_finite(), "no pairs were examined");
        assert!(
            disconnected_seen > 0,
            "the corpus contained no disconnected molecule, so the fill was never reached"
        );
    }

    /// Two independently grown fragments in one molecule, so the unreachable
    /// branch of `hop_distances` is actually exercised.
    fn disconnected_molecule(
        rng: &mut Stream,
        n: u8,
        tbl: &PeriodicTable,
        ids: &[ElementId],
    ) -> Mol12 {
        let mut mol = Mol12::new();
        let split = 1 + n / 2;
        let mut prev: Option<u8> = None;
        for i in 0..n {
            let pick = usize::try_from(rng.next_range(u64::try_from(ids.len()).unwrap())).unwrap();
            let Some(atom) = mol.add_atom(ids[pick]) else {
                break;
            };
            // The break between the two fragments: no bond crosses it.
            if i != split
                && let Some(p) = prev
            {
                let _ = mol.add_bond(p, atom, BondOrder::SINGLE, tbl);
            }
            prev = Some(atom);
        }
        mol
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
            // tree, so it retries against *other parents* (all at SINGLE; there
            // is no lower order to retry down to, which an earlier version of
            // this comment claimed). If every parent refuses, the child is left
            // isolated and this returns a non-tree — measured unreachable, 0 of
            // 10,020 trees, and structurally so, since `chain_capable` elements
            // carry valence >= 2 so a tree always retains spare capacity.
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
