//! Binding — one operation, used everywhere (spec §8.3).
//!
//! Two molecules bind when their signatures are **complementary**: shapes
//! interlock and surface characters oppose. Every higher-level phenomenon in
//! Borbax is this function at a different scale — molecules sticking together,
//! an enzyme recognising a substrate, a membrane self-assembling, a pore
//! gating passage. If a feature seems to need a second mechanism, question the
//! feature (§8.3).
//!
//! # The `ANTI` composition, which is the thing to get right
//!
//! Two bodies in contact touch along **opposite** directions: A's surface
//! reaching along `+d` meets B's surface reaching along `−d`. So the kernel
//! indexes the partner through [`Geodesic::contact_perms`], never
//! [`Geodesic::rotation_perms`] — the former is `anti ∘ rotation` precomputed.
//!
//! Getting this wrong is silent and total. Since `mirror(B)[i] = B[anti[i]]`,
//! pairing `a[i]` with `b[rotation_perms(r)[i]]` computes A's affinity with
//! **B's enantiomer** — so enantiomers become interchangeable and homochirality
//! (§22.8) is impossible by construction rather than emergent. Nothing fails.
//! Applying `anti` a second time on top of `contact_perms` cancels and lands
//! the same defect.
//!
//! `the_kernel_reads_the_contact_table` is the required probe, and it is stated
//! as a substitution because the tables were renamed under an earlier wording:
//! replace `g.contact_perms(r)` with `g.rotation_perms(r)` and **both**
//! `a_physical_complement_scores_zero` and `the_mirror_complement_does_not_fit`
//! must fail. If only one does, the pair is not discriminating.

#![expect(
    clippy::many_single_char_names,
    reason = "§8.3's kernel is written in the spec's own vocabulary — `a` and `b` are the \
              two bodies, `g` the geodesic, `k` the binding constants, `i`/`j` the paired \
              directions. Renaming them to `sig_a`/`geodesic` would make the code read \
              further from the pseudocode it implements, which is where the ANTI hazard \
              lives."
)]

use crate::geodesic::{Geodesic, Rotation};
use crate::signature::Signature;
use borbax_units::Span;
use borbax_universe::UniverseConsts;

/// The universe's binding constants, gathered from [`UniverseConsts`].
///
/// Fields are private and there is exactly one constructor, so `w_shape` and
/// `w_charge` — drawn from the identical range `[0.6, 1.4]` — cannot be
/// transposed at a second construction site.
///
/// # `ideal_gap` is deliberately absent, and the algebra says why
///
/// §8.3 writes the shape term as `−(r_A[i] + r_B[j] − IDEAL_GAP)²` against a
/// drawn constant. [`fit`] computes the separation instead. That is not a
/// deviation for convenience — a constant cannot do the job the spec assigns
/// it, provably:
///
/// Write `r_A[i] = m_A + u_A[i]` with `Σu = 0`, and `δ = m_A + m_B − L`. Since
/// the direction permutation is a bijection the cross terms vanish, so for
/// **any** separation `L`:
///
/// ```text
/// Σᵢ (r_A[i] + r_B[j] − L)²  =  D·δ²  +  Σᵢ (u_A[i] + u_B[j])²
/// ```
///
/// Only the second term depends on the rotation. So the choice of `L` cannot
/// affect **which** pose wins — it only adds a pose-blind magnitude. Two
/// consequences, and both close off a "keep the constant" repair:
///
/// - An **additive** gap, `L = m_A + m_B + g`, makes `δ = −g` and contributes
///   `D·g²` — identical for every pair, so a uniform offset that changes no
///   ranking whatsoever. A pure no-op.
/// - A **multiplicative** one, `L = γ(m_A + m_B)`, makes `δ = (1−γ)(m_A + m_B)`
///   and contributes `D(1−γ)²(m_A + m_B)²` — a penalty on *combined size* with
///   no shape in it, which is exactly the defect being removed. `γ = 1` is the
///   only value that does not reintroduce it.
///
/// The field remains on [`UniverseConsts`]: deleting a drawn constant shifts
/// every subsequent draw from that stream and moves every element in every
/// universe, which is a far larger change than this task should make. It is
/// unread by binding, and Task 20 should decide whether to retire it or give
/// it a job that a constant can actually do.
///
/// `Copy`, 16 bytes: this belongs in registers, unlike the geodesic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BindConsts {
    w_shape: f64,
    w_charge: f64,
}

impl BindConsts {
    /// Gather the binding constants for a universe.
    ///
    /// No `Default`, deliberately: every one of these is drawn from the seed,
    /// and a default would be a physics constant with no universe behind it.
    #[must_use]
    pub const fn of(c: &UniverseConsts) -> Self {
        Self {
            w_shape: c.w_shape,
            w_charge: c.w_charge,
        }
    }
}

/// A shape's rotation-invariant summary, for the pre-filter (§8.3, §8.6).
///
/// Four numbers that do not change under any of the 60 rotations, so they can
/// be computed **once per species at intern time** and reused against every
/// partner. That is §8.6's rule, and the natural spelling of `may_bind` —
/// summing the extents inside the per-pair call — violates it by recomputing a
/// per-species constant once per partner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SigSummary {
    mean_extent: Span,
    /// Euclidean norm of the mean-centred extents: how much relief the surface
    /// has, independent of how big the molecule is.
    extent_spread: f64,
    mean_character: f64,
    character_spread: f64,
}

/// How well two shapes fit, and in which relative orientation.
///
/// **The orientation is returned, not discarded, and that is a design decision
/// rather than a convenience.** §8.3's own formulation is `max over R`, which
/// answers "can these two fit *somehow*" — a scalar. §10.1 needs more than
/// that: membrane sheets require the poses to be mutually *consistent*, so that
/// A's orientation against B constrains B's against C. A max-over-orientation
/// scalar cannot express that, and the literature is unambiguous that an
/// isotropic pairwise rule yields spots rather than sheets (Ono & Ikegami,
/// ECAL 1999, needed an explicit orientation state plus an alignment term).
///
/// Borbax has the anisotropy — a signature is direction-resolved and the
/// 60-rotation search *is* an orientation search — and throwing it away at the
/// API boundary would make §10.1 unreachable without recomputing the kernel
/// per collision, which §8.6 forbids. Keeping it costs one `u8` per pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fit {
    score: f64,
    pose: Rotation,
}

impl Fit {
    /// The affinity score. Always `<= 0`; zero is a perfect fit.
    #[must_use]
    pub const fn score(self) -> f64 {
        self.score
    }

    /// The winning orientation of the partner body.
    ///
    /// **Meaningful only against the exact signatures that were scored.**
    /// Canonicalising an input composes a rotation into the answer, so a
    /// renderer handed this index must be handed the same [`Signature`] values
    /// the kernel saw — not a re-derivation from an [`Embedding`].
    ///
    /// [`Embedding`]: crate::layout::Embedding
    #[must_use]
    pub const fn pose(self) -> Rotation {
        self.pose
    }
}

/// Mean extent, summed in direction-index order (§13.1).
///
/// `sum * inv_n` rather than `sum / n`, matching
/// [`Embedding::radius_of_gyration`], whose doc deferred this exact choice
/// "until it acquires a consumer". This is that consumer. The two spellings
/// differ in the last bits on ~31% of real pairs and move the winning rotation
/// on ~0.9%, so it is a decision rather than a preference; it is pinned here
/// and by `the_binding_digest_is_pinned`.
///
/// [`Embedding::radius_of_gyration`]: crate::layout::Embedding::radius_of_gyration
#[must_use]
fn mean_extent<const D: usize>(s: &Signature<D>) -> Span {
    let mut sum = Span::ZERO;
    for e in s.extents() {
        sum += *e;
    }
    #[expect(
        clippy::cast_precision_loss,
        clippy::as_conversions,
        reason = "D is 12, 42 or 162 — the widening is exact"
    )]
    let inv_n = 1.0 / D as f64;
    sum * inv_n
}

#[must_use]
fn mean_character<const D: usize>(s: &Signature<D>) -> f64 {
    let mut sum = 0.0f64;
    for c in s.characters() {
        sum += *c;
    }
    #[expect(
        clippy::cast_precision_loss,
        clippy::as_conversions,
        reason = "D is 12, 42 or 162 — the widening is exact"
    )]
    let inv_n = 1.0 / D as f64;
    sum * inv_n
}

impl<const D: usize> Signature<D> {
    /// The rotation-invariant summary [`ceiling`] needs (§8.6).
    ///
    /// Computed once per species at intern time, never per pair.
    #[must_use]
    pub fn summary(&self) -> SigSummary {
        let (me, mc) = (mean_extent(self), mean_character(self));
        let (mut se, mut sc) = (0.0f64, 0.0f64);
        for (e, c) in self.extents().iter().zip(self.characters()) {
            let de = (*e - me).get();
            let dc = *c - mc;
            se += de * de;
            sc += dc * dc;
        }
        SigSummary {
            mean_extent: me,
            extent_spread: se.sqrt(),
            mean_character: mc,
            character_spread: sc.sqrt(),
        }
    }
}

/// A rotation-invariant **upper bound** on [`fit`]'s score (§8.3).
///
/// §8.3 asks for "a genuine bound — a filter that cannot state the property it
/// guarantees is not conservative, it is merely untested". The property, stated:
///
/// > For every rotation `R`, `score(R) <= ceiling(A, B)`. So
/// > `ceiling(A, B) < threshold` proves `max_R score(R) < threshold`, and a
/// > rejection provably excludes no pair the full search would have accepted.
///
/// The derivation, per channel. With `x`, `y` the mean-centred channels and
/// `δ` the mean offset, `Σ(x_i + y_j + δ)² = ‖x‖² + ‖y‖² + D·δ² + 2Σx_i y_j`,
/// and Cauchy–Schwarz bounds the cross term by `‖x‖‖y‖`, giving
/// `≥ (‖x‖ − ‖y‖)² + D·δ²` for **every** permutation. Both weights are drawn
/// strictly positive, so summing the two channels preserves the direction.
///
/// **This bound and the per-pair separation had to arrive together.** The
/// obvious bound — built from the mean defect alone — becomes *identically
/// zero* once `L*` is derived rather than drawn, because `δ ≡ 0` is precisely
/// what makes `L*` optimal. It would then reject nothing, silently, with no
/// test failing. Measured: rejection rate 69.4% before, **0.0%** after, against
/// **94.7%** for the bound above. That matters beyond tidiness — it is the
/// difference between Task 20's D = 162 sweep taking 11 hours and taking
/// minutes.
///
/// Takes [`SigSummary`], not [`Signature`], so §8.6's per-species rule is
/// enforced by the type rather than remembered.
#[must_use]
pub fn ceiling<const D: usize>(a: &SigSummary, b: &SigSummary, k: BindConsts) -> f64 {
    #[expect(
        clippy::cast_precision_loss,
        clippy::as_conversions,
        reason = "D is 12, 42 or 162 — the widening is exact"
    )]
    let d = D as f64;
    // **The shape channel contributes no `D·δ²` term, and its absence is the
    // point rather than an omission.** `fit` separates the bodies by
    // `mean(a) + mean(b)`, which is exactly the value making `δ = 0` — that is
    // what "least-squares optimal" means. The charge channel has no such free
    // parameter, so its offset `ε = μ_A + μ_B` survives and is the term that
    // currently dominates it.
    //
    // If a drawn separation is ever reintroduced, this is where `D·δ²` comes
    // back, and omitting it would make the bound unsound rather than merely
    // loose.
    let dr = a.extent_spread - b.extent_spread;
    let da = a.character_spread - b.character_spread;
    let eps = a.mean_character + b.mean_character;
    -(k.w_shape * (dr * dr) + k.w_charge * (da * da + d * eps * eps))
}

/// Score two shapes against each other, over all 60 proper rotations (§8.3).
///
/// ```text
/// L* = mean(r_A) + mean(r_B)               // the best-fitting separation
/// for each rotation R:
///     for each direction i:
///         j = contact_perms(R)[i]          // anti ∘ rotation, precomputed
///         shape  = −( r_A[i] + r_B[j] − L* )²
///         charge = −( a_A[i] + a_B[j] )²
///     score(R) = Σᵢ ( w_shape · shape + w_charge · charge )
/// fit = argmax over R
/// ```
///
/// # `L*` is derived, not chosen
///
/// The separation that minimises `Σ(r_A[i] + r_B[j] − L)²` is the mean of the
/// summands, and because the permutation is a bijection that mean is
/// `mean(r_A) + mean(r_B)` for every rotation alike. So `L*` is the
/// **least-squares optimal separation**: the kernel places the two bodies where
/// they fit best and scores the residual, which is what "how well do these
/// interlock" means. It is not a tuning knob and there is nothing to sweep.
///
/// What this replaces is §8.3's drawn `IDEAL_GAP`, and [`BindConsts`] carries
/// the proof that no constant can do the job. Measured before the change: the
/// drawn gap sits near 2 while real extents run to 13, so `r_A + r_B < gap` in
/// **0 of 2520** samples per seed — the term never reached its optimum and
/// collapsed to a penalty on total size, correlating with combined molecular
/// size at `|r| ≥ 0.987` *within a fixed atom count*. After: `0.43`–`0.47`, and
/// the spread across the 60 poses rises from 7.5–9.7% of the score to
/// 220–260%. The rotation search goes from a ripple to the dominant term.
///
/// **Not normalised by `L*²`.** Dividing through would make the shape term
/// dimensionless and the kernel scale-free, which is tempting since `w_shape`
/// currently multiplies a squared length while `w_charge` multiplies a pure
/// number. It is held back deliberately: the RAF-specificity literature argues
/// against a fully scale-free binding rule — a monomer and a long polymer
/// should not compete on equal terms — and it is a separable decision that
/// cannot change which pose wins (`L*` is constant over both `i` and `R`).
/// Routed to Task 20 with the dimensional mismatch named.
///
/// # Why the two accumulators are not one
///
/// Shape and charge are summed separately and combined once per rotation.
/// Interleaving them into a single accumulator is algebraically identical and
/// **moves the winning rotation in ~14% of real pairs**, because the summation
/// order changes in the last bits and ~12% of pairs have an exact tie for the
/// argmax. It is also *slower*: the two accumulators pair into one packed
/// SIMD add, lane 0 shape, lane 1 charge, each lane still summing in exact
/// direction order — so the determinism-preserving form is the fast one and
/// there is no trade here. `the_binding_digest_is_pinned` catches a change.
///
/// # Why the partner is gathered rather than materialised
///
/// Building `b.permuted(contact_perms(r))` per rotation gives a bit-identical
/// *score* — the contact coset is closed under inversion — while returning a
/// **different** [`Fit::pose`] in ~75% of pairs, because `permuted` writes
/// `out[perm[i]] = self[i]` and the kernel reads `b[perm[i]]`: inverse
/// pairings. It is also 1.4–1.9× slower. A test asserting only the score would
/// pass; `the_pose_survives_a_permuted_partner` asserts the pose.
///
/// Accumulated in direction-index order and never through `.sum()` over an
/// unordered iterator, so the summation order is fixed on every platform
/// (§13.1).
#[must_use]
pub fn fit<const D: usize>(
    a: &Signature<D>,
    b: &Signature<D>,
    g: &Geodesic<D>,
    k: BindConsts,
) -> Fit {
    // The least-squares optimal separation, hoisted: it is invariant over both
    // the direction and the rotation, so computing it inside either loop is
    // pure waste — and measured, the naive spelling costs 125x at D = 162
    // because the optimiser does not always hoist it.
    let sep = (mean_extent(a) + mean_extent(b)).get();
    let (mut best, mut pose) = (f64::NEG_INFINITY, Rotation::IDENTITY);
    for rot in Rotation::all() {
        let perm = g.contact_perms(rot);
        let (mut shape, mut charge) = (0.0f64, 0.0f64);
        for i in 0..D {
            let Some(&j) = perm.get(i) else { continue };
            let j = usize::from(j);
            let (Some(ra), Some(rb)) = (a.extents().get(i), b.extents().get(j)) else {
                continue;
            };
            let (Some(ca), Some(cb)) = (a.characters().get(i), b.characters().get(j)) else {
                continue;
            };
            let ds = (*ra + *rb).get() - sep;
            let dc = *ca + *cb;
            shape += ds * ds;
            charge += dc * dc;
        }
        let score = -(k.w_shape * shape + k.w_charge * charge);
        // **Strictly greater, so the first of an exact tie wins.** ~12% of real
        // pairs tie at the bitwise maximum, and 42% of those ties are float
        // coincidences rather than molecular symmetry — so this comparison,
        // together with `Rotation::all()`'s pinned order, decides the rendered
        // pose in one pair in eight. `>=` would silently pick the *last*.
        // Load-bearing; pinned by `the_binding_digest_is_pinned`.
        if score > best {
            best = score;
            pose = rot;
        }
    }
    Fit { score: best, pose }
}

/// [`fit`]'s score alone, for callers that do not need the orientation.
#[must_use]
pub fn affinity<const D: usize>(
    a: &Signature<D>,
    b: &Signature<D>,
    g: &Geodesic<D>,
    k: BindConsts,
) -> f64 {
    fit(a, b, g, k).score()
}

#[cfg(test)]
#[expect(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    reason = "CLAUDE.md permits both inside #[cfg(test)] with a stated reason; every index \
              here is a loop bound over D or over an atom count already established"
)]
mod tests {
    use super::*;
    use crate::canonical::{CanonMol, canonicalise};
    use crate::graph::Mol12;
    use crate::layout::embed;
    use crate::signature::signature;
    use borbax_rng::{Domain, Stream};
    use borbax_universe::{
        BondOrder, ElementId, PeriodicTable, Universe, element::generate_elements,
    };

    const D: usize = 42;

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
    fn sig(mol: &Mol12, uni: &Universe, g: &Geodesic<D>) -> Signature<D> {
        signature(&embed(&canon(mol), uni), g)
    }
    fn random_tree(rng: &mut Stream, n: u8, tbl: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = Mol12::default();
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

    /// A signature built by hand, so a complement can be constructed exactly.
    ///
    /// Test-only, and it must stay that way: a public constructor would let a
    /// caller mint a shape no molecule can produce and feed it to the kernel.
    fn from_parts(r: &[Span; D], a: &[f64; D]) -> Signature<D> {
        Signature::from_parts_for_test(*r, *a)
    }

    /// **The exact complement, built through `anti`.**
    ///
    /// `b.r[anti[i]] = gap − a.r[i]` and `b.a[anti[i]] = −a.a[i]`, so at the
    /// identity rotation every `ds` and `dc` is exactly zero.
    ///
    /// **A complement built at the same index inverts and rewards the defect.**
    /// With `b.r[i] = gap − a.r[i]` the anti-less kernel gives `ds = 0` in every
    /// direction and scores exactly 0 — the global maximum — while the *correct*
    /// kernel scores strictly negative. An implementer who writes the obvious
    /// fixture sees the correct kernel score badly, "fixes" it until it scores
    /// zero, and lands the missing-`anti` defect with a green test.
    fn complement_through_anti(a: &Signature<D>, g: &Geodesic<D>) -> Signature<D> {
        // **Every extent is positive, and that is the fix showing itself.**
        // A perfect complement needs `r_A[i] + r_B[j]` constant; ANY constant
        // will do, because `fit` derives the separation from the pair. So the
        // constant is chosen above `max(a.r)` and the complement is a shape a
        // universe could contain.
        //
        // Under the drawn `ideal_gap` this was impossible: the constant was
        // forced to ~2 while real extents run to 13, so a complement scoring
        // zero had to hold *negative* extents. The spec's own required probe
        // could only pass by being unphysical, which is independent evidence
        // that the constant separation was the error.
        let mut top = Span::ZERO;
        for e in a.extents() {
            if *e > top {
                top = *e;
            }
        }
        let level = top + Span(1.0);
        let (mut r, mut c) = ([Span::ZERO; D], [0.0f64; D]);
        for i in 0..D {
            let j = usize::from(g.anti()[i]);
            r[j] = level - a.extents()[i];
            c[j] = -a.characters()[i];
        }
        for e in &r {
            assert!(
                e.get() > 0.0,
                "the complement fixture has a negative extent"
            );
        }
        from_parts(&r, &c)
    }

    /// The same construction *without* `anti` — the trap fixture, kept so the
    /// probe can show it inverts.
    fn complement_same_index(a: &Signature<D>) -> Signature<D> {
        let mut top = Span::ZERO;
        for e in a.extents() {
            if *e > top {
                top = *e;
            }
        }
        let level = top + Span(1.0);
        let (mut r, mut c) = ([Span::ZERO; D], [0.0f64; D]);
        for i in 0..D {
            r[i] = level - a.extents()[i];
            c[i] = -a.characters()[i];
        }
        from_parts(&r, &c)
    }

    fn asymmetric_fixture(seed: u64) -> (Signature<D>, Universe, Geodesic<D>) {
        let (tbl, uni) = fixture(seed);
        let g = geo();
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(6100 + seed, Domain::Molecule, 0);
        let s = sig(&random_tree(&mut rng, 9, &tbl, &ids), &uni, &g);
        (s, uni, g)
    }

    // ------------------------------------------------------------------
    // The ANTI probe — CLAUDE.md makes this pair mandatory
    // ------------------------------------------------------------------

    /// **A physical complement scores zero**, which is the global maximum of a
    /// negated sum of squares.
    ///
    /// **Not bit-exactly, and the reason is worth having.** `ds` collapses to
    /// `a.r[i] + (gap − a.r[i]) − gap`, which is exact only when
    /// `gap − a.r[i]` is representable without rounding — true at
    /// `gap = 2.0`, where an earlier measurement found 10 752 of 10 752 entries
    /// exactly zero, and false at the drawn gaps in `[1.6, 2.8]`. Measured
    /// here: seed 6 gives exactly `0.0`, seed 17 gives `−1.3e-30`.
    ///
    /// The residual is the square of one rounding of `gap`, summed over D and
    /// weighted — about `w · D · ulp(gap)² ≈ 2e-30`. **The bar sits in an empty
    /// band 10 orders above that and 21 orders below the smallest real score**
    /// (best real pair measured: −50.96), so it cannot be reached by rounding
    /// and cannot be missed by a genuine defect. A relative tolerance is wrong
    /// here — the quantity being asserted *is* zero, so there is nothing to be
    /// relative to.
    #[test]
    fn a_physical_complement_scores_zero() {
        for seed in [6u64, 17, 42] {
            let (a, uni, g) = asymmetric_fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let b = complement_through_anti(&a, &g);
            let score = fit(&a, &b, &g, k).score();
            assert!(
                score.abs() < 1e-20,
                "seed {seed}: the exact complement scored {score}, not zero. Rounding \
                 accounts for ~2e-30; anything above 1e-20 is a real defect."
            );
        }
    }

    /// **The mirror complement does not fit**, and its margin is exactly the
    /// headroom the test above has under the missing-`anti` substitution.
    ///
    /// These two are not independent tests: since
    /// `broken(a, X) = correct(a, mirror(X))`, this assertion's value *is* the
    /// other's failure magnitude. Stated because a later reader will otherwise
    /// delete one as redundant.
    #[test]
    fn the_mirror_complement_does_not_fit() {
        for seed in [6u64, 17, 42] {
            let (a, uni, g) = asymmetric_fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let b = complement_same_index(&a);
            let score = fit(&a, &b, &g, k).score();
            assert!(
                score < -1e-6,
                "seed {seed}: the mirror complement scored {score}, i.e. it fit — \
                 the kernel is indexing the partner without the antipode"
            );
        }
    }

    /// **The fixture must be chirally discriminating, or the pair above proves
    /// nothing.** If no proper rotation separates the signature from its
    /// mirror, both tests read zero under both kernels and the probe is inert.
    /// Measured: chains of 1–3 atoms are achiral at every seed, and ~23% of
    /// plausible fixtures are achiral, so this is a live hazard rather than a
    /// theoretical one.
    #[test]
    fn the_probe_fixture_is_chirally_discriminating() {
        for seed in [6u64, 17, 42] {
            let (a, _uni, g) = asymmetric_fixture(seed);
            let mut margin = f64::INFINITY;
            for rot in Rotation::all() {
                let mut acc = 0.0f64;
                for (i, &to) in g.contact_perms(rot).iter().enumerate() {
                    let j = usize::from(to);
                    let dr = (a.extents()[i] - a.extents()[j]).get();
                    let dc = a.characters()[i] - a.characters()[j];
                    acc += dr * dr + dc * dc;
                }
                if acc < margin {
                    margin = acc;
                }
            }
            assert!(
                margin > 1e-3,
                "seed {seed}: fixture chirality margin is {margin}, so the ANTI probe \
                 cannot discriminate — measured bimodal, either 0 or >= 2.3e-4"
            );
        }
    }

    // ------------------------------------------------------------------
    // Properties
    // ------------------------------------------------------------------

    /// Affinity is a negated sum of squares, so zero is the ceiling.
    #[test]
    fn affinity_is_never_positive() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let k = BindConsts::of(&uni.consts);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(6200, Domain::Molecule, 0);
        let mols: Vec<_> = (0..12)
            .map(|i| sig(&random_tree(&mut rng, 3 + (i % 9), &tbl, &ids), &uni, &g))
            .collect();
        for i in 0..mols.len() {
            for j in i..mols.len() {
                let s = fit(&mols[i], &mols[j], &g, k).score();
                assert!(s <= 0.0 && s.is_finite(), "pair ({i},{j}) scored {s}");
            }
        }
    }

    /// Swapping the arguments must not change the score by more than the
    /// summation order can account for.
    ///
    /// Relative, not absolute: scores run to 10³, so the plan's absolute `1e-9`
    /// had three orders of slack and would not fire on a real regression.
    #[test]
    fn the_score_is_symmetric_to_within_summation_order() {
        let (tbl, uni) = fixture(6);
        let g = geo();
        let k = BindConsts::of(&uni.consts);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(6300, Domain::Molecule, 0);
        let mols: Vec<_> = (0..10)
            .map(|i| sig(&random_tree(&mut rng, 4 + (i % 8), &tbl, &ids), &uni, &g))
            .collect();
        let mut checked = 0u32;
        for i in 0..mols.len() {
            for j in (i + 1)..mols.len() {
                let (ab, ba) = (
                    fit(&mols[i], &mols[j], &g, k).score(),
                    fit(&mols[j], &mols[i], &g, k).score(),
                );
                // The comparison form, not `f64::max` — §13.1 bans the method
                // as non-deterministic for ±0.0. Both are finite here.
                let scale = if ab.abs() > ba.abs() {
                    ab.abs()
                } else {
                    ba.abs()
                };
                assert!(
                    (ab - ba).abs() <= 1e-12 * scale,
                    "pair ({i},{j}): {ab} vs {ba}"
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 45, "the corpus did not run to completion");
    }

    /// **The fix, asserted rather than described: binding is no longer a size
    /// comparison.**
    ///
    /// Scored size-matched — both molecules the same atom count — so the
    /// measurement cannot be explained by "big things pair with big things".
    /// Under the drawn `ideal_gap` this correlation was `|r| >= 0.987` in every
    /// size class across three seeds; the derived separation brings it to
    /// roughly 0.43–0.47.
    ///
    /// A bar of 0.8 leaves wide margin on both sides: it is far above anything
    /// the fix produces and far below the 0.987 the defect produced, so it
    /// cannot be met by a partial repair and cannot fire on noise.
    #[test]
    fn the_score_is_not_a_size_comparison() {
        let g = geo();
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6500 + seed, Domain::Molecule, 0);
            let mut worst = 0.0f64;
            for n in [5u8, 8, 11] {
                let mols: Vec<_> = (0..10)
                    .map(|_| sig(&random_tree(&mut rng, n, &tbl, &ids), &uni, &g))
                    .collect();
                let (mut sizes, mut scores) = (Vec::new(), Vec::new());
                for i in 0..mols.len() {
                    for j in (i + 1)..mols.len() {
                        let mut size = Span::ZERO;
                        for e in mols[i].extents().iter().chain(mols[j].extents()) {
                            size += *e;
                        }
                        sizes.push(size.get());
                        scores.push(fit(&mols[i], &mols[j], &g, k).score());
                    }
                }
                let r = correlation(&sizes, &scores).abs();
                if r > worst {
                    worst = r;
                }
            }
            assert!(
                worst < 0.8,
                "seed {seed}: |corr(size, score)| is {worst} size-matched — binding is \
                 still ranking by size rather than by shape. The drawn-gap kernel \
                 measured >= 0.987 here and the derived separation measures ~0.45."
            );
        }
    }

    /// **The pre-filter never rejects a pair the full search would accept**,
    /// which is the property §8.3 demands a filter be able to state.
    ///
    /// Asserted against every rotation individually rather than against the
    /// maximum — strictly stronger, free, and it catches a bound that dominates
    /// the max by accident.
    #[test]
    fn the_prefilter_is_a_genuine_bound() {
        let g = geo();
        let (mut checked, mut slack_seen) = (0u32, 0u32);
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6600 + seed, Domain::Molecule, 0);
            let mols: Vec<_> = (0..14)
                .map(|i| sig(&random_tree(&mut rng, 3 + (i % 10), &tbl, &ids), &uni, &g))
                .collect();
            for i in 0..mols.len() {
                for j in i..mols.len() {
                    let bound = ceiling::<D>(&mols[i].summary(), &mols[j].summary(), k);
                    let best = fit(&mols[i], &mols[j], &g, k).score();
                    assert!(
                        best <= bound + 1e-9 * bound.abs(),
                        "seed {seed} pair ({i},{j}): the search found {best}, above the \
                         ceiling {bound} — the filter would reject a binding pair"
                    );
                    if bound - best > 1e-6 {
                        slack_seen += 1;
                    }
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 3 * 105, "the corpus did not run to completion");
        // A bound equal to the score everywhere would be suspicious: it would
        // mean the ceiling is the search, and the filter buys nothing.
        assert!(
            slack_seen > 200,
            "only {slack_seen} of {checked} pairs had any slack between the ceiling \
             and the true maximum"
        );
    }

    /// Pearson correlation. Test-only; nothing in the simulation reads it.
    fn correlation(x: &[f64], y: &[f64]) -> f64 {
        #[expect(
            clippy::cast_precision_loss,
            clippy::as_conversions,
            reason = "a corpus count, far below 2^53"
        )]
        let n = x.len() as f64;
        let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
        let (mut sxy, mut sxx, mut syy) = (0.0f64, 0.0f64, 0.0f64);
        for (a, b) in x.iter().zip(y) {
            let (dx, dy) = (a - mx, b - my);
            sxy += dx * dy;
            sxx += dx * dx;
            syy += dy * dy;
        }
        sxy / (sxx.sqrt() * syy.sqrt())
    }

    /// **The pose is not preserved by materialising the partner**, which is why
    /// the kernel gathers instead.
    ///
    /// `permuted` writes `out[perm[i]] = self[i]` while the kernel reads
    /// `b[perm[i]]` — inverse pairings. The contact coset is closed under
    /// inversion so the *score* is unchanged, which is exactly what makes this
    /// dangerous: a test on the score alone passes while the rendered pose
    /// differs. Asserts the score agrees and records that the pose need not.
    #[test]
    fn the_pose_survives_a_permuted_partner() {
        let (a, uni, g) = asymmetric_fixture(6);
        let (tbl2, _) = fixture(6);
        let ids = chain_capable(&tbl2);
        let mut rng = Stream::new(6400, Domain::Molecule, 0);
        let b = sig(&random_tree(&mut rng, 7, &tbl2, &ids), &uni, &g);
        let k = BindConsts::of(&uni.consts);

        let direct = fit(&a, &b, &g, k);
        let mut moved = 0u32;
        for rot in Rotation::all() {
            let posed = b.permuted(g.rotation_perms(rot));
            let f = fit(&a, &posed, &g, k);
            // The score set is closed under the group, so the max is unchanged.
            assert!(
                (f.score() - direct.score()).abs() <= 1e-9 * direct.score().abs(),
                "rotating the partner changed the score: {} vs {}",
                f.score(),
                direct.score()
            );
            if f.pose() != direct.pose() {
                moved += 1;
            }
        }
        assert!(
            moved > 30,
            "only {moved} of 60 rotations moved the pose; if this is 0 the kernel has \
             stopped depending on the partner's orientation at all"
        );
    }
}
