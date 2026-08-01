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
//! **The required probe is a substitution, not a test** — there is no test to
//! run, which is why it is stated here in full. It is phrased as an edit
//! because the tables were renamed under an earlier wording that named a
//! deleted symbol: replace `g.contact_perms(r)` with `g.rotation_perms(r)` in
//! [`fit`] — the kernel's single call site — and **both**
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

    /// The same constants with the charge channel switched off.
    ///
    /// **Test-only, and it exists to break a cancellation rather than to model
    /// anything.** The shape and charge channels have opposite-signed size
    /// dependence, so the combined score's size correlation is smaller than
    /// either channel's. Asserting on the combination alone would credit the
    /// fix for a coincidence — and would fire spuriously the day §7.1's
    /// affinity draw is repaired and the cancellation disappears.
    #[cfg(test)]
    #[must_use]
    pub(crate) const fn with_charge_weight_zero(self) -> Self {
        Self {
            w_shape: self.w_shape,
            w_charge: 0.0,
        }
    }
}

/// A shape's rotation-invariant summary, for the pre-filter (§8.3, §8.6).
///
/// Three numbers that do not change under any of the 60 rotations, so they can
/// be computed **once per species at intern time** and reused against every
/// partner. That is §8.6's rule, and the natural spelling of `may_bind` —
/// summing the extents inside the per-pair call — violates it by recomputing a
/// per-species constant once per partner.
///
/// **It was four, and the fourth was written and never read.** `mean_extent`
/// sat here because the shape channel looked like it would need it; it does
/// not, because [`fit`] separates the bodies at exactly the value that makes
/// the shape channel's offset `δ` zero, so no `D·δ²` term ever reaches
/// [`ceiling`]. A field nothing reads is worse than a missing one:
/// `#[derive(PartialEq)]` keeps it alive against dead-code analysis, and the
/// next reader assumes it is load-bearing. Task 12 can add it back with a
/// consumer attached.
///
/// **The obvious objection, answered rather than left open:** [`fit`] computes
/// `mean_extent(a) + mean_extent(b)` per *pair*, and this type's own first
/// paragraph calls recomputing a per-species constant per partner a §8.6
/// violation. It is the same shape and not the same size. §8.6's rule is about
/// the pipeline whose cost forced it — canonicalisation, embedding, signature
/// construction, folding — each orders of magnitude above a `D`-element sum.
/// Two such sums are ~2·D against the kernel's 60·D, so restoring the field
/// would remove **3%** of `fit`'s work and reintroduce a field written once and
/// read once. The change worth making is not this one: it is `fit` taking
/// summaries alongside signatures, which is an API change Task 12 should make
/// when it wires a threshold to [`ceiling`] and gains a reason to hold both.
///
/// **Carries `D`, and that is load-bearing rather than tidy.** Without it
/// `Signature<12>::summary()` and `Signature<42>::summary()` produce the same
/// type, so `ceiling::<42>(&s12.summary(), ..)` compiles and answers — and
/// answers wrongly: measured, the search beats the supposed upper bound on
/// **93.4%** of pairs, worst overshoot 2.47x. §22.2 sweeps D across
/// {12, 42, 162}, which is exactly where a copied turbofish comes from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SigSummary<const D: usize> {
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
/// per collision, which §8.6 forbids. Keeping it costs one [`Rotation`] per
/// pair — a `usize` newtype, so 8 bytes. (This said "one `u8`" until a reviewer
/// read the definition: the *index* fits in a `u8`, the type it is stored in
/// does not, and those are not the same claim.)
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
    /// # Apply the TRANSPOSE, `R_rᵀ`, not `R_r`
    ///
    /// This is the single most repeatable mistake in this codebase and it has
    /// already been made once: Task 7 shipped a `contact_perms` docstring
    /// naming the wrong pose while the formula beside it was right, and
    /// `geodesic.rs` records that a renderer applying `R_r` "draws a plausible,
    /// silent, wrong docking picture". That warning lived on the *table*; this
    /// is the value a renderer is actually handed, so it belongs here too.
    ///
    /// Measured rather than asserted — physically rotating B's embedding and
    /// rescoring against A: the transpose reproduces [`Fit::score`] on **8 of
    /// 8** real pairs, the forward rotation on **2 of 8**, and those two are the
    /// order-2 elements where `R = Rᵀ`. A worked case: kernel −81.906652,
    /// `Rᵀ`-posed −81.906652, `R`-posed −104.143087.
    ///
    /// The reason is that a support function transforms contravariantly —
    /// `h_{MB}(v) = h_B(Mᵀv)` — so reading B's signature through a permutation
    /// built from `R_r` corresponds to rotating B's *body* by `R_rᵀ`.
    /// `the_returned_pose_is_the_transpose` ties this to real geometry rather
    /// than to another table, which is what Task 7 lacked.
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
/// differ in the last bits on 33.3% of real pairs and move the winning rotation
/// on 0.4%, so it is a decision rather than a preference; it is pinned here and
/// by `the_binding_digest_is_pinned`.
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
    pub fn summary(&self) -> SigSummary<D> {
        let (me, mc) = (mean_extent(self), mean_character(self));
        // **Direction-index order, and it is a pin like every other in this
        // module.** Reversing this loop moves `extent_spread` on 54.7% of
        // species with the whole suite green — it feeds `ceiling`, which is a
        // rejection decision once Task 12 wires a threshold to it.
        let (mut se, mut sc) = (0.0f64, 0.0f64);
        for (e, c) in self.extents().iter().zip(self.characters()) {
            let de = (*e - me).get();
            let dc = *c - mc;
            se += de * de;
            sc += dc * dc;
        }
        SigSummary {
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
/// **That holds exactly in real arithmetic and up to accumulated rounding in
/// floating point, and the difference is reachable.** On an exact complement
/// the two spreads are equal, so `‖x‖ − ‖y‖` should be zero; it rounds to
/// ~1e-15, giving a ceiling of ~−1e-30 while the score rounds to exactly
/// `−0.0`. Measured on this module's own complement fixture: seed 17 gives
/// `score = −0.0` against `ceiling = −1.377e-30`, so the score is *above* the
/// ceiling by 1.4e-30.
///
/// The consequence for a caller is nil — a rejection threshold sits at the
/// scale of real scores, tens to thousands, so a discrepancy 30 orders down
/// cannot flip one. The consequence for a *test* is not nil: a relative
/// tolerance collapses to nothing near zero, which is why
/// `the_prefilter_is_a_genuine_bound` carries an absolute floor as well.
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
pub fn ceiling<const D: usize>(a: &SigSummary<D>, b: &SigSummary<D>, k: BindConsts) -> f64 {
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
    // back. **Omitting a term from the bound is always safe — it can only make
    // the ceiling larger, hence looser.** An earlier version of this comment
    // said the opposite, which is the dangerous direction to be wrong in:
    // *adding* a term is what needs proof, because it lowers the ceiling and
    // can push it below a score the search will actually find.
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
///     shape, charge = 0, 0                 // reset per rotation
///     for each direction i:
///         j = contact_perms(R)[i]          // anti ∘ rotation, precomputed
///         shape  += ( r_A[i] + r_B[j] − L* )²      // two accumulators,
///         charge += ( a_A[i] + a_B[j] )²           // never one
///     score(R) = −( w_shape · shape + w_charge · charge )   // combined once
/// fit = argmax over R
/// ```
///
/// **The two accumulators are in the block because they are load-bearing.** It
/// read `Σᵢ (w_shape · shape + w_charge · charge)` — one accumulator over the
/// weighted combination — which is precisely the interleaving the section below
/// measures and rejects. A reader implementing from the pseudocode would have
/// written the form the prose forbids three paragraphs later.
///
/// The reset line is there for the same reason and was missing from the first
/// repair: with `+=` and no per-rotation reset the accumulators carry across
/// rotations, the score decreases monotonically, and `argmax` is always the
/// identity. A doc fix that introduces a second defect while correcting the
/// first is the pattern this project keeps recording, so it is named here
/// rather than quietly amended.
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
/// should not compete on equal terms.
///
/// **It is NOT pose-neutral, and an earlier version of this paragraph said it
/// was.** The reasoning was that `L*` is constant over both `i` and `R`, which
/// is true and does not give the conclusion: the score is a weighted sum of
/// *two* channels, so rescaling one of them changes their ratio and moves the
/// argmax. Measured — dividing only the shape channel by `L*²` moves the
/// winning pose on **36–65%** of pairs. It is a physics change requiring
/// golden regeneration, not a cosmetic one. Routed to Task 20 with the
/// dimensional mismatch named and that correction attached, because "cannot
/// change which pose wins" is exactly the sentence that would let someone land
/// it as a no-op.
///
/// # Why the two accumulators are not one
///
/// Shape and charge are summed separately and combined once per rotation.
/// Interleaving them into a single accumulator is algebraically identical and
/// moves the **score** on 73.8% of pairs and the **winning rotation** on 1.9%.
/// It is also *slower*: the two accumulators pair into one packed
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
/// pass; `materialising_the_partner_moves_the_pose` asserts the pose.
///
/// Accumulated in direction-index order and never through `.sum()` over an
/// unordered iterator, so the summation order is fixed on every platform
/// (§13.1).
///
/// # `fit(A, B)` and `fit(B, A)` are not bitwise equal, and Task 12 owns it
///
/// The two argument orders sum the same terms in different orders, so the
/// scores differ in the last bits on **61.2%** of hetero-pairs, by up to
/// **8 ulp**. The *pose* is unaffected — measured, `contact_perms(pose_BA)` is
/// the exact inverse permutation of `contact_perms(pose_AB)` in 6 624 of 6 624
/// pairs — because the derived separation raised the inter-pose spread to
/// 220–260% of the score, far above any rounding.
///
/// It is documented rather than fixed here because the fix belongs to whoever
/// owns the memo. **The ordering key must be an interned species id, never a
/// float comparison over the signatures**: `borbax_units::canonical_cmp`
/// reports two NaNs `Equal` while `==` reports them unequal, so a float-keyed
/// order cannot break a tie between two *distinct* species stably, and a memo
/// keyed on it would return one pair's affinity for another's. Task 12 mints
/// the ids and is the first caller that can hit this.
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
        // **Iterated on A's side, not indexed, for the reason
        // `Signature::group_distance` is.** The previous spelling read
        // `perm.get(i)` and `a.extents().get(i)` and `continue`d on `None` —
        // and `continue` in a scoring loop does not fail, it silently returns a
        // *partial* sum for that rotation. None of those `None`s is reachable
        // (`perm` is `&[u8; D]` and `i < D`), which is exactly what makes the
        // spelling dangerous: an unreachable branch whose behaviour, if it ever
        // did run, is a wrong score rather than a stop.
        for ((ra, ca), &j) in a
            .extents()
            .iter()
            .zip(a.characters().iter())
            .zip(perm.iter())
        {
            let j = usize::from(j);
            // B is reached through the permutation, so it cannot be iterated
            // alongside. `unreachable!` names its own precondition — the
            // spelling CLAUDE.md keeps deliberately available — where a
            // `continue` would have hidden the violation in the arithmetic.
            let rb = b
                .extents()
                .get(j)
                .unwrap_or_else(|| unreachable!("contact_perms is a permutation of 0..{D}"));
            let cb = b
                .characters()
                .get(j)
                .unwrap_or_else(|| unreachable!("contact_perms is a permutation of 0..{D}"));
            let ds = (*ra + *rb).get() - sep;
            let dc = *ca + *cb;
            shape += ds * ds;
            charge += dc * dc;
        }
        let score = -(k.w_shape * shape + k.w_charge * charge);
        // **Strictly greater, so the first of an exact tie wins**, together
        // with `Rotation::all()`'s pinned order. `>=` would silently pick the
        // last, moving the pose on 1.8% of pairs while moving **zero** scores —
        // caught by `the_binding_digest_is_pinned` alone, which is why that
        // digest hashes the pose separately.
        //
        // **Ties are exclusively self-pairs**, and that is the fact worth
        // carrying rather than a global rate: measured over 12 seeds, 38 of 192
        // A-against-A pairs tie (**19.8%**) and **0 of 1440** hetero-pairs do.
        // They are structural — a molecule's own automorphisms — not float
        // coincidences: compensated summation leaves 204 of 205 standing.
        // Self-binding is §10.1's membrane self-assembly case, so this
        // tie-break is load-bearing exactly where it looked least interesting.
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
    clippy::indexing_slicing,
    reason = "CLAUDE.md permits this inside #[cfg(test)] with a stated reason; every index \
              here is a loop bound over D or over an atom count already established. There \
              is deliberately no `unwrap_used` here — the fixtures that needed it moved to \
              `testkit`, and leaving the expectation behind would suppress a lint this \
              module no longer earns."
)]
mod tests {
    use super::*;
    use crate::graph::Mol12;
    use crate::layout::embed;
    use crate::signature::signature;
    use crate::testkit::{canon, chain_capable, fixture, random_tree};
    use borbax_rng::{Domain, Stream};
    use borbax_universe::Universe;

    const D: usize = 42;

    fn geo() -> Geodesic<D> {
        Geodesic::<D>::build().unwrap_or_else(|_| unreachable!("D=42 builds"))
    }
    fn sig(mol: &Mol12, uni: &Universe, g: &Geodesic<D>) -> Signature<D> {
        signature(&embed(&canon(mol), uni), g)
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
    /// Chains of 1–3 atoms are achiral at every seed and ~23% of plausible
    /// fixtures are achiral, so this is a live hazard rather than a theoretical
    /// one.
    ///
    /// **The bar is 1.0 and the measured margins are 13.48, 9.60 and 1.48.**
    /// Seed 42 is the tight one — within 50% of the bar — and that is stated
    /// rather than smoothed, because a bar chosen far below every observation
    /// would pass a fixture that had drifted most of the way to achiral. The
    /// distribution is bimodal: exactly 0 for an achiral fixture, and the
    /// smallest non-zero margin measured anywhere is ~2.3e-4, so anything above
    /// that separates the two modes; 1.0 additionally rejects a fixture that is
    /// merely *nearly* achiral.
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
                margin > 1.0,
                "seed {seed}: fixture chirality margin is {margin}, below the bar of \
                 1.0 — the ANTI probe cannot discriminate on a fixture this close to \
                 achiral. Measured at 13.48 / 9.60 / 1.48 for seeds 6 / 17 / 42."
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
                // **Tight, and the arithmetic behind the number was wrong.**
                // The previous bound was `1e-12 * scale` — about 4 500 ulp —
                // which had 560x slack and would not have fired if the
                // asymmetry grew by two orders of magnitude.
                //
                // This comment said `f64::EPSILON * 16.0` is "~8 ulp of the
                // larger operand". It is not: `f64::EPSILON` is `2^-52`, which
                // is one ulp only for operands in `[1, 2)`. For `scale` in
                // `[2^k, 2^(k+1))` one ulp is `2^(k-52)`, so this bar admits
                // **16 to 32 ulp** of `scale`. The measured maximum is 8 ulp
                // over 6 624 pairs, so the bar has 2-4x slack rather than the
                // 2x the old wording implied — still 1000x tighter than what it
                // replaced, and the figure is now the one the expression means.
                assert!(
                    (ab - ba).abs() <= f64::EPSILON * 16.0 * scale,
                    "pair ({i},{j}) differ by more than the 16-32 ulp the summation \
                     order accounts for: {ab} vs {ba}"
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 45, "the corpus did not run to completion");
    }

    /// **A bit-level digest of the score AND the winning pose.**
    ///
    /// Three shapes in this module are pinned in prose and nothing else
    /// enforces them: the two-accumulator split, the direction-index summation
    /// order, and `mean_extent`'s `sum * inv_n`. Every one is *algebraically*
    /// invariant under the tidy-up a reader would reach for, so the suite stays
    /// green while every number moves in its last bits.
    ///
    /// **The pose is hashed separately from the score, and that is the whole
    /// reason this test earns its keep.** Flipping the tie-break from `>` to
    /// `>=` moves **zero** scores and 1.8% of poses — so a score-only digest is
    /// green under it, and so is every tolerance test at any epsilon. §14.5's
    /// renderer consumes that index and §10.1 constrains membrane sheets with
    /// it, so an unpinned pose is an unpinned physics input.
    ///
    /// The mutations this catches, each measured alone against the shipped
    /// kernel with the rest of the suite green: `sum / n` for `sum * inv_n`
    /// (score 14.2%, pose 0.4%), one accumulator instead of two (73.8% / 1.9%),
    /// the tie-break (0% / 1.8%), the direction loop reversed (66.8% / 2.2%),
    /// and `summary()`'s spread loop reversed (54.7% of species).
    ///
    /// Hashes through `canonical_bits`, not `to_bits`: a runtime NaN's sign and
    /// `-0.0` are architecture-dependent and §13.6's matrix would report that as
    /// a simulation divergence (§13.4).
    #[test]
    fn the_binding_digest_is_pinned() {
        let g = geo();
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let (mut scores, mut poses) = (0u64, 0u64);
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6700 + seed, Domain::Molecule, 0);
            let mols: Vec<_> = (0..9)
                .map(|i| sig(&random_tree(&mut rng, 3 + (i % 9), &tbl, &ids), &uni, &g))
                .collect();
            for i in 0..mols.len() {
                for j in i..mols.len() {
                    let f = fit(&mols[i], &mols[j], &g, k);
                    hash ^= borbax_units::canonical_bits(f.score());
                    hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                    scores += 1;
                    // Separate mixing step, so a pose change cannot cancel
                    // against a score change.
                    hash ^= u64::try_from(f.pose().index()).unwrap_or(u64::MAX);
                    hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                    poses += 1;
                }
            }
        }
        // 3 seeds x 45 unordered pairs of 9 molecules, one of each per pair.
        assert_eq!(scores, 3 * 45, "the corpus did not run to completion");
        assert_eq!(poses, scores, "a pose was not hashed for every score");
        assert_eq!(
            hash, 0xab87_c6f7_2cf5_3428,
            "binding moved. If you meant to change the physics, regenerate this \
             constant and say so in the commit message. If you did not, suspect an \
             accumulation order: the two accumulators, the direction-index sum, and \
             mean_extent's `sum * inv_n` are each algebraically invariant under the \
             obvious tidy-up and each move the last bits. If only the POSE moved, \
             suspect the `>` tie-break: ties here are **exclusively self-pairs**, \
             which are automorphisms rather than float coincidences, and this \
             corpus has 9 of them per seed."
        );
    }

    /// **The returned pose is `R_rᵀ`, tied to real geometry rather than to
    /// another table.**
    ///
    /// Task 7 shipped a docstring naming the wrong pose because every test it
    /// had compared tables against tables — the right invariant, saying nothing
    /// about what the composition *means*. This rotates B's embedding by the
    /// candidate matrix, recomputes its signature from scratch, and rescores.
    /// The transpose must reproduce [`Fit::score`]; the forward rotation must
    /// not, except on the order-2 elements where they coincide.
    #[test]
    fn the_returned_pose_is_the_transpose() {
        use crate::geodesic::rotation_matrices;
        let g = geo();
        let (tbl, uni) = fixture(6);
        let k = BindConsts::of(&uni.consts);
        let ids = chain_capable(&tbl);
        let mut rng = Stream::new(777, Domain::Molecule, 0);
        let mats = rotation_matrices().unwrap_or_else(|_| unreachable!("the matrices build"));

        let (mut transpose_ok, mut forward_ok, mut checked) = (0u32, 0u32, 0u32);
        for _ in 0..8 {
            let ea = embed(&canon(&random_tree(&mut rng, 8, &tbl, &ids)), &uni);
            let eb = embed(&canon(&random_tree(&mut rng, 7, &tbl, &ids)), &uni);
            let (sa, sb) = (signature(&ea, &g), signature(&eb, &g));
            let f = fit(&sa, &sb, &g, k);
            let m = &mats[f.pose().index()];
            let transposed = [
                [m[0][0], m[1][0], m[2][0]],
                [m[0][1], m[1][1], m[2][1]],
                [m[0][2], m[1][2], m[2][2]],
            ];
            // Score A against B *physically* rotated, reading antipodally —
            // which is what a renderer does when it places the two bodies.
            let score_posed = |rot: &[[f64; 3]; 3]| {
                let st = signature(&eb.rotated(rot), &g);
                let sep = (mean_extent(&sa) + mean_extent(&st)).get();
                let (mut sh, mut ch) = (0.0f64, 0.0f64);
                for i in 0..D {
                    let j = usize::from(g.anti()[i]);
                    let ds = (sa.extents()[i] + st.extents()[j]).get() - sep;
                    let dc = sa.characters()[i] + st.characters()[j];
                    sh += ds * ds;
                    ch += dc * dc;
                }
                -(k.w_shape * sh + k.w_charge * ch)
            };
            // **Absolute floor as well as relative.** A purely relative
            // tolerance is exactly zero at a zero score, so both comparisons
            // below would be false and the test would fail on a *correct*
            // kernel. Scores do reach `-0.0` — `ceiling`'s doc records it on
            // complement shapes — and the sibling assertion in
            // `the_prefilter_is_a_genuine_bound` already carries a floor for
            // this reason.
            let tol = 1e-9 * f.score().abs() + 1e-20;
            if (score_posed(&transposed) - f.score()).abs() < tol {
                transpose_ok += 1;
            }
            if (score_posed(m) - f.score()).abs() < tol {
                forward_ok += 1;
            }
            checked += 1;
        }
        assert_eq!(checked, 8, "the corpus did not run to completion");
        assert_eq!(
            transpose_ok, checked,
            "the transpose did not reproduce the kernel's score on every pair — \
             §14.5's renderer would draw a wrong docking pose"
        );
        // **The contrast is the point.** If the forward rotation also matched
        // everywhere, the test would prove nothing about which convention is
        // right. It matches only on the order-2 elements where `R = Rᵀ`;
        // measured 2 of 8.
        assert!(
            forward_ok < checked,
            "the forward rotation reproduced the score on all {checked} pairs, so this \
             fixture cannot distinguish the two conventions — pick molecules whose \
             winning poses are not all order-2"
        );
    }

    /// **The fix, asserted rather than described: binding is no longer a size
    /// comparison.**
    ///
    /// Scored size-matched — both molecules the same atom count — so the
    /// measurement cannot be explained by "big things pair with big things".
    /// Under the drawn `ideal_gap` this correlation was `|r| >= 0.987` in every
    /// size class across three seeds. The derived separation brings the worst
    /// per-seed value to **0.583 / 0.510 / 0.390** — not the "0.43–0.47" an
    /// earlier version of this doc quoted, which was measured on a corpus the
    /// fix commit then changed.
    ///
    /// # The whole-score assertion passes partly by cancellation, so it is not
    /// the only one
    ///
    /// Measured per channel, size-matched: the shape term correlates with size
    /// at **−0.52 / −0.67 / −0.62** and the charge term at **+0.67 / +0.66 /
    /// +0.73**. They have opposite signs and largely cancel, which is why the
    /// combined score reads 0.39–0.58. The charge term's size dependence exists
    /// only because §7.1's `affinity` never straddles zero — an open item routed
    /// to the affinity draw — so **fixing that would remove the cancellation and
    /// this bar would fire on a kernel that got better**, with a message
    /// blaming the wrong channel.
    ///
    /// So the shape channel is asserted **on its own**, with the charge weight
    /// zeroed. That is a property of the kernel; the combined figure is a
    /// property of a coincidence, and it is kept only as a regression tripwire.
    ///
    /// # The shape channel is improved, not cured, and the number says so
    ///
    /// Shape-only, size-matched: **0.53 / 0.57 / 0.84** for seeds 6 / 17 / 42,
    /// against `>= 0.987` for the drawn constant. Seed 42 is barely improved
    /// and the bar here is 0.90, which is a low bar honestly placed rather than
    /// a comfortable one.
    ///
    /// The reason is structural and worth stating so nobody reads the fix as
    /// complete: removing `D·δ²` deletes the term that carried **all** the
    /// combined-size dependence, but the residual `Σ(u_A[i] + u_B[j])²` is an
    /// **unnormalised squared length**, so it still grows with how much relief
    /// a surface has — and bigger molecules have more. That residual is exactly
    /// what dividing by `L*²` would remove, which is the decision deferred to
    /// Task 20 above. The two are the same question seen twice: this fix takes
    /// the half that is forced by algebra and leaves the half that is a
    /// judgement about whether binding should be scale-free.
    ///
    /// The 0.8 bar catches a full reintroduction (0.987) and a multiplicative
    /// separation at `γ <= 0.80`; measured, it **misses** `γ` in roughly
    /// `[0.82, 1.18]`. So it is a tripwire, not a proof, and the per-channel
    /// assertion below is what carries the claim.
    #[test]
    fn the_score_is_not_a_size_comparison() {
        let g = geo();
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6500 + seed, Domain::Molecule, 0);
            let (mut worst, mut worst_shape) = (0.0f64, 0.0f64);
            for n in [5u8, 8, 11] {
                let mols: Vec<_> = (0..10)
                    .map(|_| sig(&random_tree(&mut rng, n, &tbl, &ids), &uni, &g))
                    .collect();
                let (mut sizes, mut scores, mut pairs) = (Vec::new(), Vec::new(), Vec::new());
                for i in 0..mols.len() {
                    for j in (i + 1)..mols.len() {
                        let mut size = Span::ZERO;
                        for e in mols[i].extents().iter().chain(mols[j].extents()) {
                            size += *e;
                        }
                        sizes.push(size.get());
                        scores.push(fit(&mols[i], &mols[j], &g, k).score());
                        pairs.push((i, j));
                    }
                }
                let r = correlation(&sizes, &scores).abs();
                assert!(
                    r.is_finite(),
                    "seed {seed} n={n}: the correlation is not finite — a degenerate size \
                     class would otherwise leave `worst` at 0.0 and pass vacuously"
                );
                if r > worst {
                    worst = r;
                }
                // The shape channel alone, which is what the fix actually
                // changed. `w_charge = 0` removes the cancellation described
                // above.
                let shape_only = BindConsts::of(&uni.consts).with_charge_weight_zero();
                let bare: Vec<f64> = pairs
                    .iter()
                    .map(|&(i, j)| fit(&mols[i], &mols[j], &g, shape_only).score())
                    .collect();
                let rs = correlation(&sizes, &bare).abs();
                if rs > worst_shape {
                    worst_shape = rs;
                }
            }
            assert!(
                worst < 0.8,
                "seed {seed}: |corr(size, score)| is {worst} size-matched — binding is \
                 still ranking by size rather than by shape. The drawn-gap kernel \
                 measured >= 0.987 here and the derived separation measures 0.39-0.58."
            );
            // **The claim that actually rests on the fix.** Measured
            // 0.53 / 0.57 / 0.84 with the charge channel off, for seeds
            // 6 / 17 / 42; the defect gave >= 0.987. This comment said
            // "0.52-0.67", which are the *combined* kernel's per-seed shape
            // components from the round before — a stale figure contradicting
            // the assertion message four lines below it.
            assert!(
                worst_shape < 0.90,
                "seed {seed}: the SHAPE channel alone correlates with size at \
                 {worst_shape}, above the measured 0.53 / 0.57 / 0.84. This is the \
                 assertion the fix owns — unlike the combined score it cannot be \
                 rescued by the charge channel's opposite-signed size dependence \
                 cancelling it out."
            );
        }
    }

    /// **The pre-filter never rejects a pair the full search would accept**,
    /// which is the property §8.3 demands a filter be able to state.
    ///
    /// Asserted against [`fit`]'s score, which **is** the maximum over the 60
    /// rotations. That is sufficient and not a shortcut: `max_R score(R) <=
    /// ceiling` implies `score(R) <= ceiling` for every `R`, and the converse
    /// holds too because the max is one of them. The two formulations are
    /// equivalent, not ranked.
    ///
    /// An earlier version of this comment said the test asserted "against every
    /// rotation individually rather than against the maximum — strictly
    /// stronger". It was wrong twice: the code has always used the maximum, and
    /// the per-rotation form is equivalent to it rather than stronger.
    #[test]
    fn the_prefilter_is_a_genuine_bound() {
        let g = geo();
        let (mut checked, mut worst_gap) = (0u32, 0.0f64);
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
                    // No turbofish: `SigSummary<D>` carries the resolution,
                    // so `D` is inferred from both arguments and a mismatch is
                    // a type error. Spelling it again re-creates the affordance
                    // the type change removed.
                    let bound = ceiling(&mols[i].summary(), &mols[j].summary(), k);
                    let best = fit(&mols[i], &mols[j], &g, k).score();
                    // Absolute floor as well as relative: on an exact
                    // complement the relative term collapses to ~1e-39 while
                    // the float discrepancy is ~1e-30. See `ceiling`'s doc.
                    assert!(
                        best <= bound + 1e-9 * bound.abs() + 1e-20,
                        "seed {seed} pair ({i},{j}): the search found {best}, above the \
                         ceiling {bound} — the filter would reject a binding pair"
                    );
                    // **Tightness, not slack.** An earlier version counted
                    // pairs where the ceiling was *loose* and asserted that
                    // count was large — which a useless bound satisfies
                    // maximally. Measured: replacing `ceiling`'s body with
                    // `0.0` — the loosest sound bound, rejecting nothing ever —
                    // passed every test in this crate.
                    // The comparison form, not `f64::max` — §13.1. Both
                    // operands are finite and positive here.
                    let scale = if best.abs() > 1.0 { best.abs() } else { 1.0 };
                    let gap = (bound - best) / scale;
                    if gap > worst_gap {
                        worst_gap = gap;
                    }
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 3 * 105, "the corpus did not run to completion");
        // **The bound must be tight enough to be worth having.** `ceiling` is
        // justified in its own doc by turning Task 20's D=162 sweep from hours
        // into minutes, and a bound that never rejects delivers none of that.
        // Measured worst relative gap on this corpus: 0.93. A ceiling of `0.0`
        // gives 1.0 for every pair and fails here.
        assert!(
            worst_gap < 0.98,
            "the worst ceiling-to-score gap is {worst_gap}, so the bound is close to \
             vacuous — `ceiling` returning a constant 0.0 would score 1.0 here and is \
             sound but useless. The filter has to reject things to be worth its branch."
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

    /// **Materialising the partner moves the pose**, which is why the kernel
    /// gathers instead.
    ///
    /// `permuted` writes `out[perm[i]] = self[i]` while the kernel reads
    /// `b[perm[i]]` — inverse pairings. The contact coset is closed under
    /// inversion so the *score* is unchanged, which is exactly what makes this
    /// dangerous: a test on the score alone passes while the rendered pose
    /// differs. Asserts the score agrees and records that the pose need not.
    #[test]
    fn materialising_the_partner_moves_the_pose() {
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
