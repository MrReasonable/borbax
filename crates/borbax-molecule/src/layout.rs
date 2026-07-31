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
    ///
    /// **`sum * inv_n`, not `sum / n`, and that shape is deliberately *not*
    /// pinned with a measured number the way the centring block is** — because
    /// this has no caller yet. It is the natural scale normaliser for Task 9's
    /// signature, so it becomes result-affecting then. Pin it at the point it
    /// acquires a consumer, not before; the same applies to `distance`'s
    /// `(s₀² + s₁²) + s₂²` association, which *is* already on the `stress` path.
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
/// **And the reachability argument is a schedule, not a property.** Nothing in
/// V0 interns a disconnected species — but only because `ReactionKind::Cleave`
/// names no products in either plan file, which the chemistry plan itself flags
/// as a defect to close for §9.1 microscopic reversibility. Nothing *enforces*
/// it: `Mol12::clear_bond` is public, `canonicalise` accepts a disconnected
/// molecule, and `embed` takes a `CanonMol` with no connectivity check. The
/// guard expires the moment Cleave is finished.
///
/// One claim here was attacked and survived: that unreachable pairs get the
/// largest target, so fragments are pushed apart. Radii vary up to ~3.6x, so a
/// cross-fragment pair of small atoms could in principle sit below a bonded pair
/// of large ones — refuted over 600 disconnected fixtures across 10 universe
/// seeds, cross target largest in 600/600. It holds for realistic universes as a
/// coincidence of the drawn radius range rather than as a structural guarantee.
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
/// is evidence about the objective, not about the solver. The unweighted sum is
/// not exported at all — see the note below on why a same-signature twin was the
/// wrong way to keep it.
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
/// Measured by
/// `a_one_atom_edit_moves_the_shape_less_than_an_unrelated_molecule_does`, on a
/// matched fixture, scored on the group-minimised distance (1.0 = perfect
/// locality, 0.5 = none, and a composition-only null on the same corpus is
/// 0.5151):
///
/// ```text
/// index-keyed START   0.6878
/// landmark start      0.8092
/// ```
///
/// **Those are this crate's numbers.** An earlier version of this block gave
/// `0.5104 / 0.7955` here — figures from the reviewing lane's own harness, which
/// is a different descriptor on a different corpus — while the test 400 lines
/// below reported different values for the nominally same statistic. Two tables,
/// one file, one statistic, two answers, both labelled "measured on this
/// embedder". That is the exact failure the paragraph above this one is about.
///
/// The damage was purely orientational: a rotation-*invariant* shape metric is
/// unchanged between the two starts, and mean stress at 240 barely moves. Both
/// find the same shapes, to the same quality, in different frames.
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
/// **Both the coordinates and the landmark rule are load-bearing; the tie-break
/// direction is not.** An earlier version of this note said the pivot choice was
/// merely a quality refinement, on the evidence that fixing the pivots to indices
/// 1 and 2 left the locality test green. That evidence was worthless: the test
/// scored the *un-aligned* distance at the time, which ranks embeddings backwards
/// (see `group_distance`). Measured on the group-minimised distance instead, the
/// pivot rule spans **0.7624 to 0.8356** across eight variants — a 7.3-point range
/// on a statistic whose total headroom above chance is about 0.32.
///
/// What genuinely does not matter, measured: the tie-break direction (flipping
/// `>` to `>=` moves group locality by +0.0045), and making the *first* landmark
/// index-free — two structural alternatives both scored **below** canonical index
/// 0 (0.7888 and 0.7956 against 0.8213), so index 0 as the first pivot is
/// empirically fine and is not the weak point it looks like.
///
/// Landmarks are picked from the graph metric, so the construction stays a pure
/// function of the species and introduces no randomness and no transcendental —
/// `experiments/src/embed.rs` adds a `Stream` jitter here and this deliberately
/// does not, because that would make `borbax-rng` a runtime dependency in a
/// result-affecting path.
fn initial_layout(target: &[[f64; MAX_ATOMS]; MAX_ATOMS], n: usize) -> [[f64; 3]; MAX_ATOMS] {
    // Argmax by explicit comparison with a lowest-index tie-break, never
    // `max_by`/`total_cmp`: §13.4 disallows both.
    //
    // **Ties are genuinely reached and the tie-break is genuinely a choice.**
    // Censused over 5212 molecules: the first landmark ties at the maximum in
    // 4.5%, and in 2.9% of all molecules that tie is between atoms with different
    // sorted target rows — sufficient to show they are not in one automorphism
    // orbit, so symmetry does not absorb the decision. Separately, 68.3% of
    // parent/mutant pairs change the `(radius, degree)` tag of at least one pivot,
    // so a pivot's physical identity is not continuous under an edit at all.
    //
    // It nonetheless does not cost anything: flipping the tie-break to
    // highest-index moves group locality by +0.0045. The number is recorded rather
    // than the reassurance alone, because a future reader deserves to know the
    // choice was measured rather than assumed away.
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
///
/// **No pivoting, and the reason is a proof rather than a lint.**
/// `xᵀ(V_w + J)x = xᵀV_w x + (1ᵀx)²/n`. `V_w` is a weighted Laplacian with
/// non-negative weights, so it is PSD with null space exactly `span{1}` — the
/// weight graph is `K_n`, every target being strictly positive. If `1ᵀx ≠ 0` the
/// second term is positive; if `1ᵀx = 0` then `x ⊥ null(V_w)` and the first is.
/// Positive definite for every input, independent of the molecular graph, which
/// only sets hop counts and never the sparsity of `w`. Measured over 28,722 cases
/// across seven universe seeds — paths, rings, stars, all-isolated dust, every
/// two-fragment split, K(3,3), K(6,6), Petersen — zero failures, worst surviving
/// pivot 1.6e-2 at 0.119 of its original diagonal, worst κ₂ 88.2 against an
/// analytic bound of 1.2e4.
///
/// An earlier version of this note said instead that a pivot search would need
/// `max`, which §13.1 bans. **That is a doctrine error and is corrected here
/// rather than deleted**: §13.1 bans `f64::max` the *method*, for its ±0.0
/// non-determinism, and `clippy.toml`'s own prescribed replacement is
/// `if a > b { a } else { b }` — which `initial_layout` already uses for an
/// argmax. A §13.1-legal pivot search is perfectly writable. The reason none is
/// needed is above. A false determinism argument matters because it gets reused
/// later to reject something that is fine.
///
/// **Four accumulations below are pinned shapes** (§13.1). Reversing either
/// `cholesky` inner loop changes **24.43%** of entries; collapsing the pivot sum
/// to one `.sum::<f64>()` changes **28.47%**; reversing `solve`'s forward
/// substitution **23.57%** and its back substitution **41.05%**. A blocked or
/// tiled Cholesky is the same class. All four are invisible to every test here —
/// they move convergence by ~1e-16, and the suite measures convergence quality.
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
    // is one **Jacobi sweep** of `V X' = B(X)X`.
    //
    // **A third description was tried here and is also false, so it is recorded
    // rather than quietly dropped**: that it is equivalently the over-relaxed
    // Guttman transform `X + α(Γ(X) − X)` with `α = n/(n−1)`. That equivalence
    // needs `D⁻¹ = αV⁺` on the centred subspace, which forces every off-diagonal
    // weight to be equal — so it is exact at n = 2 (one pair, uniform trivially)
    // and on a uniform-weight clique, and **32% wrong at n = 12** on the `1/d²`
    // weighting this code actually uses. The reassuring "reproduced to 7.3e-15"
    // was measured on the family the code does not use.
    //
    // Three rewrites, three errors, the last in the sentence doing the persuading
    // — which is why the conclusion below rests on the n = 2 measurement and the
    // bipartite argument, neither of which needs α.
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
    // **Jacobi, not Gauss-Seidel.** `B(X)` is built entirely from the previous
    // sweep's `pos`, so reading partially-updated coordinates cannot change the
    // answer. That is a hazard worth naming because a Gauss-Seidel sweep also
    // converges, so the difference would be invisible.
    //
    // **That claim is about the mathematics, and an earlier version of it
    // overreached into the bits.** It said "the result does not depend on the
    // order atoms are visited in", sitting directly above the `j` loop — which it
    // appears to license and does not. Iterating `j` descending changes **32.99%**
    // of `rhs` entries. The `j` accumulation order is a pinned shape: `.rev()`,
    // distance-sorted traversal and chunked partials are all forbidden.
    //
    // One rewrite *is* cleared, measured rather than assumed: exploiting `B(X)`'s
    // antisymmetry to halve the loop over `i < j` is **bit-identical**, 0 of
    // 102,660 entries differing, because contributions still reach `rhs[i]` in
    // ascending `j` and both `-(a-b) == (b-a)` and `-(s*d) == s*(-d)` are exact.
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
    // **Falls through to the centring below rather than returning here**, so
    // every path out of this function yields a centred embedding. Returning early
    // gave an uncentred one, and since signatures are support functions measured
    // from the origin, an uncentred result is not a degraded shape — it is a
    // shape with an arbitrary offset baked into every extent.
    //
    // A breakdown means the weight matrix was not positive definite, which for a
    // molecule with at least two atoms and finite radii cannot happen. Carried as
    // a value rather than asserted: a panic in the profile CLAUDE.md designates
    // for beaker work is a dead simulation.
    let factor = cholesky(&w, n);

    for _ in 0..iterations {
        let Some(chol) = factor.as_ref() else { break };
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
                let wij = w.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);

                // The unit vector from `j` to `i`. **At a coincident pair it is
                // chosen deterministically rather than skipping the term, and that
                // is a correctness fix rather than a nicety.**
                //
                // **Skipping did not break descent, and saying it did would license
                // the wrong conclusion.** Measured on 267 molecules that start
                // coincident: the old form rises 0 times above 1e-12, worst 3.6e-15.
                // `w` is built independently of this skip, so `V_w` was already
                // weight-complete, and skipping is the textbook convention
                // corresponding to the trivial bound `δ >= 0`. The defect was purely
                // the saddle. This matters because "the old form broke descent"
                // would tell a future reader that any descent-preserving coincidence
                // handling is therefore correct — the old one preserved descent and
                // was still wrong.
                //
                // Skipping is wrong *here* because of how `initial_layout` works: two atoms that are
                // graph-automorphic and carry the same element have identical
                // target rows, so they start at exactly the same point — by
                // construction, not by accident. Dropping their mutual `B(X)` term
                // leaves the residual system exactly equivariant under transposing
                // them, so the solve returns equal rows and **they can never
                // separate**. An exact fixed point, in exact arithmetic.
                //
                // Measured before the fix, on an all-one-element corpus where the
                // automorphism groups are largest: **160 of 4000 molecules** placed
                // at least one pair at a single point, the worst violating a target
                // of **6.75 Span**. It surfaced as nothing: the iteration is still
                // monotone, it just converges to a saddle. `bonded_atoms_stay_apart`
                // was structurally blind — it walks bonded pairs only, and in a
                // tree two automorphic atoms are never adjacent.
                //
                // `initial_layout` cannot be repaired instead: no continuous
                // function of the graph can distinguish automorphic atoms, which is
                // precisely why it is the right initialiser.
                //
                // **The majorization survives**, so this does not trade correctness
                // for separation. At `δ_ij(Z) = 0` the bound `δ_ij(X) ≥ (x_i−x_j)·u`
                // holds for *any* unit `u` by Cauchy–Schwarz, so descent is
                // preserved for any choice — the choice only has to be deterministic
                // and antisymmetric (`u_ij = −u_ji`), or `B(X)` stops being
                // symmetric and the transform stops being the Guttman transform.
                //
                // **Antisymmetry carries a second load that is easier to miss.**
                // `cholesky` factors `V_w + J`, and `(V_w+J)⁻¹b = V_w⁺b + (1ᵀb/n)·1`,
                // so the `J` term vanishes only while `1ᵀb` is *exactly* zero. Targets
                // are bitwise symmetric (FP addition commutes), hence so are the
                // weights, and `scale·(−x) == −(scale·x)` is exact — so antisymmetry
                // is what keeps the sum exactly zero. Lose it and the step acquires a
                // spurious translation *before* the centring pass. Verified
                // exhaustively over all 132 ordered index pairs and all 4083 possible
                // coincident orbits: `‖sep‖ = 1`, both orderings pick one axis, and no
                // orbit of any size to 12 can be left with a zero differential push.
                //
                // **Antisymmetry is proven necessary and is NOT observable here, and
                // that gap is stated rather than papered over.** Substituting the
                // axis rule `(i + 2j) % 3` — deterministic and isomorphism-invariant
                // but not antisymmetric, so `u_ij` and `u_ji` land on different axes
                // — passes the entire suite, including
                // `descent_holds_on_the_coincidence_path`. The reason is that the
                // spurious translation it introduces is removed by the centring pass
                // at the end, so it cannot reach a distance. So this is a property no
                // test defends; it rests on the algebra above, and a future edit here
                // will get no signal from the suite.
                //
                // Keyed on the canonical indices, and **the reason that is sound is
                // not the one an earlier version of this comment gave.** It said the
                // atoms differ only by an automorphism, so the branches give
                // congruent shapes. Two problems: congruent is not enough — a
                // congruent-but-rotated shape is a *different* support function at
                // Task 9 — and the premise is measurably false. Brute-forcing the
                // full automorphism group over all n! permutations for n <= 9 found
                // **2 of 382** coincident pairs *not* in one orbit, equidistant from
                // all three landmarks by arithmetic coincidence rather than symmetry
                // (141 more at n = 10..12 were too expensive to test, so 2/382 is a
                // lower bound).
                //
                // The real reason needs no automorphisms and is exact: `CanonMol`'s
                // `mol()` is bit-identical for every relabelling of one molecule, so
                // `embed` receives a bit-identical input and there are no "two
                // branches" to reconcile — targets, landmarks, indices and this
                // direction are all one value. Measured: 0 mismatches over 480
                // relabellings of 60 species that start coincident, compared through
                // `to_bits()` rather than `==`.
                //
                // What that leaves open is a *locality* question, not a determinism
                // one: for those ~0.5% of pairs the axis is a genuine arbitrary
                // physical choice, which is the index-keying hazard `initial_layout`
                // was rewritten to remove, reappearing at small scale inside this
                // branch. Routed to Task 9 rather than priced here.
                //
                // A NaN `raw` takes the deterministic branch, since `NaN > 1e-9` is
                // false. **That is the *less detectable* failure, not the safe one**,
                // and an earlier version of this sentence had it backwards. Measured
                // from an all-NaN start: the old form gave an all-zero embedding that
                // `no_two_atoms_share_a_point` catches instantly; this one gives a
                // finite, well-converged, *wrong* shape (R_g 2.186 against a correct
                // 2.190), because iteration 1 makes `rhs` a pure function of indices
                // and targets and the remaining 239 are ordinary SMACOF from a
                // garbage-but-valid start. Deterministic either way, so not a §13.1
                // hazard — a laundering one. Unreachable today because targets are
                // provably finite; whoever breaks that unreachability (the routed
                // Cleave item above is the candidate) gets no signal.
                let unit = if raw > 1e-9 {
                    [diff[0] / raw, diff[1] / raw, diff[2] / raw]
                } else {
                    let mut sep = [0.0f64; 3];
                    if let Some(slot) = sep.get_mut(j.abs_diff(i) % 3) {
                        *slot = if i < j { 1.0 } else { -1.0 };
                    }
                    sep
                };
                // **`scale = wij*dij` with `unit = diff/raw` is a pinned shape, and
                // it is a second physics change this commit made.** The predecessor
                // computed `scale = wij*dij/raw` and multiplied by `diff` — the same
                // value mathematically, differing in the last bit on **35.06%** of
                // 400,000 sampled operand quadruples. Blast radius measured over 7000
                // species: **96.87%** change bits, against the 164 the coincidence fix
                // was for. So the intended change touches 164 species and the
                // incidental one touches nearly all of them.
                //
                // The rewrite to refuse is the one a performance pass would propose —
                // hoisting the division back out of the branch to save two divisions
                // per pair. Measured: **78 of 78 tests pass, debug and release**, and
                // the corpus digest moves. Nothing in the suite can see it, and it
                // would also structurally re-break the coincident branch.
                let scale = wij * dij;
                for (k, row) in rhs.iter_mut().enumerate() {
                    if let (Some(acc), Some(&dk)) = (row.get_mut(i), unit.get(k)) {
                        *acc += scale * dk;
                    }
                }
            }
        }
        for (k, row) in rhs.iter().enumerate() {
            let solved = solve(chol, row, n);
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
    use crate::geodesic::{Geodesic, Rotation, is_identity, rotation_matrices};
    use borbax_rng::{Domain, Stream};
    use borbax_universe::{
        BondOrder, ElementId, PeriodicTable, Universe, element::generate_elements,
    };

    /// A periodic table and a universe **from the same seed**.
    ///
    /// **`table(a)` with `Universe::generate(b)` for `a != b` is a chimera, and
    /// six tests here were one.** Measured on `table(17)` against universe 4: 24
    /// of the 78 chain-capable ids — **30.8%** — do not exist in that universe at
    /// all and silently take `targets`' `map_or(1.0, ..)` fallback, and of the 54
    /// that do exist, **0 of 54 agree on radius**.
    ///
    /// So roughly a third of the corpus's shape information was the constant 1.0.
    /// That manufactures exact ties in `farthest_from`, and worse, it inflates the
    /// composition confound: on the mismatched fixture a **composition-only null**
    /// — sum of atomic radii, no graph, no bonds, no solver — scores **0.7428**
    /// against the locality gate's 0.75 bar, so the gate sat 0.008 above a
    /// geometry-free null. Matched, the null collapses to ~0.47–0.52.
    ///
    /// Concordance itself is stable across 16 matched triples (0.7982–0.8533), so
    /// repairing the fixture does not rescue a number — it makes the existing one
    /// mean something.
    fn fixture(seed: u64) -> (PeriodicTable, Universe) {
        (generate_elements(seed), Universe::generate(seed))
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

    /// Signature distance minimised over the 60 proper rotations — the quantity
    /// §8.3's `affinity` actually maximises over, and the only comparison the
    /// routed Task 9 requirement permits.
    ///
    /// **The un-aligned distance is the wrong instrument here and using it made
    /// this test rank embeddings backwards.** An in-group reorientation is
    /// invisible to every named consumer — `affinity` maximises over the group and
    /// §8.2's canonicalisation quotients by it — but the raw distance sees it as
    /// total shape change. Measured: a variant that cyclically permutes the three
    /// landmark axes on `n mod 3` (a genuine frame discontinuity, but one that
    /// lands on rotation index 6 of the 60) scores **0.5083 raw — "no locality" —
    /// and 0.8213 grouped, byte-identical to the shipped embedder at 1089/1326.**
    /// Meanwhile fixing the pivots to indices 0,1,2 scores *better* raw (0.8650)
    /// while being 0.049 *worse* grouped. Raw would have failed a correct
    /// implementation and rewarded a worse one.
    fn group_distance<const D: usize>(geo: &Geodesic<D>, a: &[f64], b: &[f64]) -> f64 {
        let mut best = f64::INFINITY;
        for rot in Rotation::all() {
            let perm = geo.rotation_perms(rot);
            let mut acc = 0.0f64;
            for i in 0..D {
                let (Some(&lhs), Some(&turned_to)) = (a.get(i), perm.get(i)) else {
                    continue;
                };
                let rhs = b.get(usize::from(turned_to)).copied().unwrap_or(0.0);
                acc += (lhs - rhs) * (lhs - rhs);
            }
            if acc < best {
                best = acc;
            }
        }
        best.sqrt()
    }

    /// A support function over the geodesic directions. Task 9's `signature` does
    /// not exist yet; this is deliberately *less* informative than it will be
    /// (atom centres rather than surfaces, and no affinity channel), so it
    /// understates locality rather than flattering it.
    fn support<const D: usize>(species: &CanonMol, uni: &Universe, geo: &Geodesic<D>) -> Vec<f64> {
        let emb = embed(species, uni);
        (0..D)
            .map(|d| {
                let dir = geo.dirs().get(d).copied().unwrap_or([0.0, 0.0, 1.0]);
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
    }

    /// **The wiring proof, and it needs no corpus.** `group_distance` is invariant
    /// under applying any group element to either argument; the raw distance is
    /// not. This is what makes the concordance test below a statement about
    /// geometry rather than about frames.
    ///
    /// Three preconditions, and omitting any makes it vacuous: the rotation must
    /// be **non-identity**, the fixture must be **asymmetric**, and the rotation
    /// must be asserted to have **actually changed** the signature first.
    #[test]
    fn the_group_distance_is_rotation_invariant_and_the_raw_one_is_not() {
        let (tbl, uni) = fixture(17);
        let ids = chain_capable(&tbl);
        let geo = Geodesic::<42>::build().unwrap_or_else(|_| unreachable!("D=42 builds"));
        let mut rng = Stream::new(778, Domain::Molecule, 0);

        // An asymmetric fixture: a random tree, not a chain or a ring.
        let species = canon(&random_tree(&mut rng, 9, &tbl, &ids));
        let base = support(&species, &uni, &geo);

        let mats =
            rotation_matrices().unwrap_or_else(|_| unreachable!("the rotation matrices build"));
        let mut checked = 0u32;
        for rot in Rotation::all() {
            // Skip the identity explicitly: under it both distances are trivially
            // unchanged, so including it would let this test pass vacuously.
            if mats.get(rot.index()).is_some_and(is_identity) {
                continue;
            }
            let perm = geo.rotation_perms(rot);
            let turned: Vec<f64> = (0..42).map(|i| base[usize::from(perm[i])]).collect();
            // The mutation must be real, or everything below is trivially true.
            let raw = base
                .iter()
                .zip(&turned)
                .map(|(x, y)| (x - y) * (x - y))
                .sum::<f64>()
                .sqrt();
            if raw < 1e-9 {
                continue; // this rotation is in the fixture's own symmetry group
            }
            checked += 1;
            let grouped = group_distance(&geo, &base, &turned);
            assert!(
                grouped < 1e-9,
                "rotation {rot:?} moved the group distance by {grouped}, raw {raw}"
            );
        }
        assert!(
            checked > 30,
            "only {checked} rotations genuinely moved the fixture; it is too symmetric to prove anything"
        );
    }

    /// **The property the whole design rests on (§8.2), pinned rather than
    /// asserted.** A one-atom edit must move the shape less than an unrelated
    /// molecule of the same size does, or mutation is a random walk and nothing
    /// accumulates.
    ///
    /// Concordance is the fraction of trials where the mutant is nearer the parent
    /// than the stranger is: 1.0 perfect, **0.5 none**. Scored on
    /// [`group_distance`], not the raw one — see its docs for why the raw version
    /// ranked embeddings backwards.
    ///
    /// **Asserted as a margin over a null measured on the same corpus, not as an
    /// absolute.** An absolute bar cannot distinguish geometry from composition: a
    /// composition-only null — sum of atomic radii, no graph, no solver — scored
    /// 0.7428 on the *mismatched* fixture this test used to carry, against a 0.75
    /// bar. On a matched fixture the null collapses to ~0.5 and the margin is real.
    #[test]
    fn a_one_atom_edit_moves_the_shape_less_than_an_unrelated_molecule_does() {
        let (tbl, uni) = fixture(17);
        let ids = chain_capable(&tbl);
        let geo = Geodesic::<42>::build().unwrap_or_else(|_| unreachable!("D=42 builds"));
        let mut rng = Stream::new(777, Domain::Molecule, 0);

        let (mut nearer, mut trials) = (0u32, 0u32);
        let (mut null_nearer, mut null_trials) = (0u32, 0u32);
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
            let stranger =
                random_molecule(&mut rng, u8::try_from(mutant.len()).unwrap(), &tbl, &ids);
            if stranger.len() != mutant.len() {
                continue;
            }
            let (cp, cm, cs) = (canon(&parent), canon(&mutant), canon(&stranger));
            if cp.len() < 3 {
                continue;
            }
            let base = support(&cp, &uni, &geo);
            let dm = group_distance(&geo, &base, &support(&cm, &uni, &geo));
            let ds = group_distance(&geo, &base, &support(&cs, &uni, &geo));
            if dm < ds {
                nearer += 1;
            }
            trials += 1;

            // The null reads only the element multiset — no graph, no bonds, no
            // solver. If the gate cannot beat this, it is measuring composition.
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
            if (bc - mc).abs() < (bc - sc).abs() {
                null_nearer += 1;
            }
            null_trials += 1;
        }

        assert!(trials > 1000, "only {trials} usable trials");
        let concordance = f64::from(nearer) / f64::from(trials);
        let null = f64::from(null_nearer) / f64::from(null_trials);
        assert!(
            concordance - null > 0.20,
            "locality is not beating a composition-only null: concordance {concordance:.4} \
             ({nearer}/{trials}), null {null:.4}, margin {:.4}. 0.5 is no locality; the \
             index-keyed start scored 0.6041 raw.",
            concordance - null
        );
    }

    #[test]
    fn embedding_is_deterministic() {
        let (tbl, uni) = fixture(4);
        let ids = chain_capable(&tbl);
        let species = canon(&chain(7, &tbl, &ids));
        assert_eq!(embed(&species, &uni), embed(&species, &uni));
    }

    #[test]
    fn positions_are_finite_and_centred() {
        let (tbl, uni) = fixture(4);
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
        let (tbl, uni) = fixture(4);
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
        let (tbl, uni) = fixture(4);
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
        let (tbl, uni) = fixture(4);
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
    /// **"Bit-identical" is what is meant and slightly more than `assert_eq!`
    /// delivers**: `Embedding` derives `PartialEq`, which is `f64::eq`, so
    /// `-0.0 == +0.0` would pass. Measured live-ness: 0 of 334,734 coordinates
    /// across 15,993 species are negative zero, including a single-element corpus
    /// chosen to maximise symmetry. Latent, and it belongs with the already-open
    /// "NaN canonicalisation at the state hash" item for Task 20, where
    /// coordinates get hashed through `to_bits()`.
    ///
    /// With `embed` taking a `CanonMol` this holds by construction and the test is
    /// proving the construction rather than the solver. That is the point — the
    /// caller-dependent call has no spelling. Probed by giving `embed` the raw
    /// `Mol12` instead: fails on the first relabelling.
    #[test]
    fn the_embedding_is_a_function_of_the_species_not_the_labelling() {
        let (tbl, uni) = fixture(17);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(81, Domain::Molecule, 0);

        // **A one-element ring first, because a chain never reaches the path this
        // guarantee is newly at risk on.** The coincidence branch in
        // `embed_with_budget` picks its separation direction from *canonical
        // index*, so "identical, not merely congruent" rests on canonical order
        // being isomorphism-invariant. A chain has |Aut| = 2 and no coincident
        // atoms, so it exercises none of that. Measured over 400 symmetric
        // molecules x 5 relabellings: 2000 of 2000 bit-identical, 12 of them
        // reaching coincidence, and congruent-but-not-identical never occurred —
        // which is what licenses the `assert_eq!` rather than a shape comparison.
        let mut ring = Mol12::new();
        for _ in 0..8 {
            assert!(ring.add_atom(ids[0]).is_some());
        }
        for i in 0..8u8 {
            assert!(
                ring.add_bond(i, (i + 1) % 8, BondOrder::SINGLE, &tbl)
                    .is_ok(),
                "an eight-ring over a valence-2-capable element is legal"
            );
        }
        let ring_want = embed(&canon(&ring), &uni);
        for _ in 0..100 {
            let mut perm: Vec<u8> = (0..8).collect();
            for i in (1..perm.len()).rev() {
                let j = usize::try_from(rng.next_range(u64::try_from(i + 1).unwrap())).unwrap();
                perm.swap(i, j);
            }
            let relabelled = ring
                .relabelled(&perm)
                .unwrap_or_else(|| unreachable!("perm permutes 0..8"));
            assert_eq!(
                embed(&canon(&relabelled), &uni),
                ring_want,
                "relabelling a symmetric ring changed the embedding"
            );
        }

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
        let (tbl, uni) = fixture(4);
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
        let (tbl, uni) = fixture(17);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(82, Domain::Molecule, 0);
        for trial in 0..60u32 {
            // Random *molecules*, not only trees: rings are where the target
            // matrix stops being a tree metric, and they are exactly the graphs
            // the plan's own preamble flags as untested by `experiments`.
            // From 2, not 3: the period-2 cycle the exact transform fixes lives
            // at exactly n = 2, and the predecessor's corpus started at 3.
            //
            // Every third fixture is one element throughout, for the automorphism
            // groups. **It does not reach the coincidence branch** — measured 0 of
            // 60 at these sizes — so descent on that path is asserted by
            // `descent_holds_on_the_coincidence_path`, which filters for it rather
            // than hoping a corpus wanders in. An earlier version of this comment
            // claimed the corpus covered it; it does not.
            let n = 2 + u8::try_from(trial % 11).unwrap();
            let mol = if trial % 3 == 0 {
                symmetric_molecule(&mut rng, n, &tbl, ids[0])
            } else {
                random_molecule(&mut rng, n, &tbl, &ids)
            };
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

    /// **Descent on the coincidence path, which no other test reaches.**
    ///
    /// `stress_never_increases`' corpus starts 0 of 60 molecules coincident even
    /// with one-element fixtures at those sizes, so the branch that picks a
    /// deterministic separation direction was argued in a comment and asserted by
    /// nothing. This filters for the case instead of hoping a corpus wanders into
    /// it — the same failure `bonded_atoms_stay_apart` had.
    #[test]
    fn descent_holds_on_the_coincidence_path() {
        let (tbl, uni) = fixture(4);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(3131, Domain::Molecule, 0);
        let mut examined = 0u32;
        for _ in 0..1200 {
            let n = 4 + u8::try_from(rng.next_range(9)).unwrap();
            let mol = symmetric_molecule(&mut rng, n, &tbl, ids[0]);
            let species = canon(&mol);
            if species.len() < 3 {
                continue;
            }
            // Keep only molecules whose START is coincident — that is the branch.
            let start = initial_layout(&targets(&species, &uni), species.len());
            let mut coincident = false;
            for i in 0..species.len() {
                for j in (i + 1)..species.len() {
                    let same = (0..3).all(|k| (start[i][k] - start[j][k]).abs() < 1e-12);
                    if same {
                        coincident = true;
                    }
                }
            }
            if !coincident {
                continue;
            }
            examined += 1;

            let mut previous = f64::INFINITY;
            for step in 0..=60usize {
                let current = stress(&species, &uni, &embed_with_budget(&species, &uni, step));
                assert!(
                    current <= previous + 1e-12,
                    "step {step}: stress rose from {previous} to {current} on a \
                     molecule that starts coincident"
                );
                previous = current;
            }
        }
        assert!(
            examined > 20,
            "only {examined} molecules started coincident; this test reached nothing"
        );
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
        let (tbl, uni) = fixture(4);
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

    // **A gap this test does not close, stated rather than left implicit.** The
    // star is bounded at n = 6 by its hub's valence, and the truncation error that
    // would justify questioning 240 lives at n = 11-12. Measured against a 200,000
    // iteration reference on the current code: the *systematic* size-correlated
    // bias the review found before the landmark start is discharged — mean R_max
    // error at n = 12 is now -0.10% against the -2.03% measured on the old start,
    // an 11-20x reduction. What remains is a per-molecule *tail*: worst case
    // 1.1-5.5% at n >= 11, and it is budget-insensitive, 16x the budget buying
    // 2.7x on the worst case while p95 falls cleanly 0.90% -> 0.006%. That is a
    // slow mode in a small fraction of configurations, not truncation, so raising
    // `ITERATIONS` would not fix it and would cost a golden regeneration for
    // nothing. Routed to Task 9 with its discriminator: if signature locality
    // shows outliers, check whether they are the same molecules before blaming
    // the signature.

    /// **No two atoms may occupy one point.** Atoms are solid bodies; a collapsed
    /// pair makes an n-atom molecule present n−1 extremal points to the support
    /// function, so it reads as less bumpy than it is and two species differing
    /// only by a symmetric pair become near-indistinguishable to `affinity`. It
    /// also attacks §22.8 directly: a collapsed pair can give a chiral molecule an
    /// artificial mirror symmetry, losing homochirality candidates to an artefact.
    ///
    /// **This replaces `bonded_atoms_stay_apart`, which was structurally blind.**
    /// That test walked *bonded* pairs only, over `random_tree` fixtures — and in
    /// a tree two automorphic atoms are never adjacent, since their distances to
    /// any root differ by one. Measured on the defect it missed: 24 of 26 collapsed
    /// pairs were non-bonded, and 0 occurred in trees at all.
    ///
    /// The corpus is deliberately **one element throughout**, which maximises the
    /// automorphism group and is where the collapse is densest — but it is not the
    /// only place: a mixed-element corpus still starts 18 of 4000 molecules
    /// coincident (~0.45%) against this corpus's 171 of 4000, so "one element" is
    /// the sharpest fixture rather than the whole reach. Measured before the fix
    /// in `embed`: 160 of 4000 molecules collapsed a pair, worst violating a target
    /// of 6.75 Span. After: 0 of 4000.
    #[test]
    fn no_two_atoms_share_a_point() {
        let (tbl, uni) = fixture(4);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(4242, Domain::Molecule, 0);
        let (mut pairs, mut molecules) = (0u32, 0u32);
        let mut closest = f64::INFINITY;
        for _ in 0..600 {
            let n = 4 + u8::try_from(rng.next_range(9)).unwrap();
            let mut mol = Mol12::new();
            for _ in 0..n {
                if mol.add_atom(ids[0]).is_none() {
                    break;
                }
            }
            let size = u8::try_from(mol.len()).unwrap();
            for i in 1..size {
                let _ = mol.add_bond(i - 1, i, BondOrder::SINGLE, &tbl);
            }
            // Extra edges, so rings appear — automorphic atoms in a ring are what
            // a tree corpus can never produce.
            for _ in 0..rng.next_range(4) {
                let a = u8::try_from(rng.next_range(u64::from(size))).unwrap();
                let b = u8::try_from(rng.next_range(u64::from(size))).unwrap();
                let _ = mol.add_bond(a, b, BondOrder::SINGLE, &tbl);
            }
            let species = canon(&mol);
            if species.len() < 3 {
                continue;
            }
            molecules += 1;
            let emb = embed(&species, &uni);
            for i in 0..species.len() {
                for j in (i + 1)..species.len() {
                    let sep = emb.distance(i, j).unwrap_or(Span::ZERO).get();
                    assert!(
                        sep > 1e-9,
                        "atoms {i},{j} of a {}-atom molecule share a point (target {})",
                        species.len(),
                        targets(&species, &uni)[i][j]
                    );
                    if sep < closest {
                        closest = sep;
                    }
                }
                pairs += 1;
            }
        }
        assert!(molecules > 400, "only {molecules} usable molecules");
        assert!(pairs > 2000, "only {pairs} pairs examined");
        // A separation that merely clears 1e-9 would be a technical pass; the real
        // claim is that separations are on the scale of the molecule.
        assert!(
            closest > 0.05,
            "closest approach {closest} is not a real separation"
        );
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
    /// It is also what lets `stress` carry a `tgt <= 0.0` guard purely for the
    /// division rather than as a claim about reachable inputs.
    #[test]
    fn targets_between_distinct_atoms_are_strictly_positive() {
        let (tbl, uni) = fixture(17);
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

    /// A ring-bearing molecule built from **one element throughout**, which
    /// maximises the automorphism group and is what reaches `embed`'s coincidence
    /// branch — automorphic *same-element* atoms are the ones with identical
    /// target rows.
    fn symmetric_molecule(rng: &mut Stream, n: u8, tbl: &PeriodicTable, one: ElementId) -> Mol12 {
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
        for _ in 0..rng.next_range(3) {
            let a = u8::try_from(rng.next_range(u64::from(size))).unwrap_or(0);
            let b = u8::try_from(rng.next_range(u64::from(size))).unwrap_or(0);
            let _ = mol.add_bond(a, b, BondOrder::SINGLE, tbl);
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
