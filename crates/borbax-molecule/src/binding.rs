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
//! deleted symbol. Substitute in [`fit`], selecting the line with this — the
//! **anchor is load-bearing**, not decoration:
//!
//! ```text
//! rg -n '^\s+let perm = g\.contact_perms\(rot\);'
//! ```
//!
//! then change that one line's `contact_perms` to `rotation_perms`.
//!
//! **Two narrower strings were tried and both matched twice.** The bare
//! `g.contact_perms(rot)` also appears in
//! `the_probe_fixture_is_chirally_discriminating`, so replacing it edits a test
//! alongside the code under test and the result cannot be read. Narrowing to
//! `let perm = g.contact_perms(rot);` then matched **this paragraph**, because
//! quoting a search string inside the file it searches is enough to break it.
//! `^\s+` excludes the `//!` line and leaves exactly one hit.
//!
//! Both mistakes were caught by the probe script's own "abort unless exactly
//! one call site" guard rather than by anyone reading, and this project has
//! already recorded a probe that silently failed to apply and reported
//! all-pass. An earlier version of this paragraph also named
//! `the_returned_pose_is_the_transpose` as the second site; it does not call
//! `contact_perms` at all, it reads `g.anti()`. Naming the wrong site is worse
//! than naming none — the reader checks it, finds nothing, and concludes the
//! warning is stale while the real site stays unnamed. The `let perm =` form is unique to the kernel. That
//! is not hypothetical: this project has already recorded a probe that silently
//! failed to apply and reported all-pass because its search string appeared
//! twice.
//!
//! **Both** `a_physical_complement_scores_zero` and
//! `the_mirror_complement_does_not_fit` must fail. If only one does, the pair
//! is not discriminating. Verified after the `complement_level` extraction:
//! both fail, plus `the_binding_digest_is_pinned` and
//! `the_returned_pose_is_the_transpose`.

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
use borbax_units::{Span, canonical_cmp};
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
/// `Signature<12>::summary(&g)` and `Signature<42>::summary(&g)` produce the
/// same type, so `ceiling::<42>(&s12.summary(&g), ..)` compiles and answers — and
/// answers wrongly: measured, the search beats the supposed upper bound on
/// **93.4%** of pairs, worst overshoot 2.47x. §22.2 sweeps D across
/// {12, 42, 162}, which is exactly where a copied turbofish comes from.
///
/// **Not `Copy` — this is Decision 4's rearrangement-inequality prefilter
/// (issue #26, Task 26.1 Step 12), and it needs the full per-direction
/// channels, not three summary scalars.** The predecessor's Cauchy–Schwarz
/// bound (`‖x‖ − ‖y‖`) is exactly `0` on every self-pair — the centred
/// channels are identical, so the spread difference vanishes regardless of
/// shape — which is why it went "selectively blind on self-pairs" the
/// moment #31's affinity fix let the charge channel's `D·eps^2` term (the
/// only thing that was carrying self-pair rejection) approach zero too.
/// Measured against real `PhysicsVersion::V2` data: the predecessor's
/// `worst_gap` rises to ~1.0 (effectively vacuous) the moment V2's
/// genuinely zero-straddling affinity lands.
///
/// **Two bounds, not one, and the max of both — but read the roles
/// precisely, because an earlier version of this comment had them
/// backwards.** It claimed the flat rearrangement bound is *also* exactly
/// `0` on a self-pair, the same degeneracy as the Cauchy–Schwarz predecessor
/// above, and that the antipodal-block bound is what restores self-pair
/// discrimination. Both halves are false, caught by a `/review-pr` pass
/// (2026-08-10) that measured rather than trusted the claim, independently
/// three ways (two specialist reviewers plus a direct mutation probe run
/// against this file's own `fixture_v2` corpus):
///
/// - **The flat bound is *not* exactly zero on a self-pair.** It sums
///   *squared* opposite-sorted differences, `Σ(x[i] + x[D-1-i])²` — zero
///   only if the centred sorted vector is exactly antisymmetric about its
///   midpoint, which mean-centred real data is not. Measured over the real
///   corpus: **0 of 42** self-pairs give a flat bound of exactly (or even
///   nearly) zero, range `[0.0602, 17.1930]` across both channels. What *is*
///   exactly zero on a self-pair is the predecessor's `(‖x‖ − ‖y‖)²` — a
///   different bound, and the source of the false generalisation.
/// - **The flat bound carries the self-pair result: it wins at least one
///   channel on 41 of the 42 self-pairs in the measured corpus**
///   (`the_prefilter_is_a_genuine_bound`'s own `flat_wins_self` counter).
///   On the one remaining self-pair the block bound wins *both* channels —
///   so "block never wins a self-pair" (an earlier version of this bullet,
///   and of [`fit`]'s doc) is not the right generalisation either. What
///   makes flat the one that matters for the self-pair bar specifically:
///   dropping flat entirely breaches it (`worst_gap_self` rises from
///   0.9592 combined to 0.9877 block-only), while dropping block does not
///   (0.9665 flat-only, still under bar). See
///   `the_prefilter_is_a_genuine_bound`'s own `flat_wins_shape`/
///   `block_wins_shape`/`flat_wins_charge`/`block_wins_charge` counters for
///   the exact per-channel split (20/22 shape, 41/1 charge on self-pairs).
/// - **What the antipodal-block bound actually earns its place on is
///   hetero-pairs: it wins at least one channel on 241 of 273 real hetero
///   comparisons** in the measured corpus (`block_wins_hetero`, same
///   robustness reasoning) — the pair class the previous version of this
///   comment had backwards. Combining both bounds still measurably helps
///   there: `worst_gap_hetero` is 0.9145 combined, tighter than either
///   bound alone (0.9445 flat-only, 0.9331 block-only) — each bound catches
///   pairs the other misses, and this is what `flat_beats_block`/
///   `block_beats_flat` (251/263 of 315 pairs) verify without needing a bar
///   at all.
///
/// [`Geodesic::contact_perms`] commuting with [`Geodesic::anti`] (verified
/// exhaustively at D = 12/42/162 — every contact permutation maps antipodal
/// *pairs* to antipodal pairs) makes the antipodal-block relaxation a
/// strictly smaller *search space* than "any permutation of `D` elements" —
/// but that does **not**, on its own, make it a tighter bound: this
/// function minimises the `s`-half and `d`-half sums over *independent*
/// orderings (see its own doc), a further relaxation the subset argument
/// alone does not cover, which is exactly why neither bound dominates the
/// other above and `ceiling` takes the max of both rather than assuming
/// one always wins. `ceiling` takes that max per channel regardless of
/// which bound is doing the work for a given pair, so a physics change
/// that moves the figures above changes no code and no golden — see
/// `the_prefilter_is_a_genuine_bound`'s own pinned counters for the guard
/// that exists to catch this doc going stale again (it did not exist when
/// this doc went stale the first time, mid-`/review-pr` round 2, which is
/// why this correction needed a *third* round to surface).
#[derive(Debug, Clone, PartialEq)]
pub struct SigSummary<const D: usize> {
    /// Mean-centred extents (`r[i] - mean(r)`), sorted ascending via
    /// [`canonical_cmp`] — the array [`ceiling`]'s flat rearrangement bound
    /// reads.
    extent_sorted: [f64; D],
    /// Antipodal-block statistics on the mean-centred extents: `[0..D/2)`
    /// holds block sums `s = u[i] + u[anti[i]]`, `[D/2..D)` holds block
    /// differences `|u[i] - u[anti[i]]|` — both halves sorted ascending
    /// independently. One `[f64; D]` rather than two `[f64; D/2]`, because
    /// `D/2` as a const-generic parameter needs unstable
    /// `generic_const_exprs`.
    extent_blocks: [f64; D],
    /// Mean-centred characters, sorted ascending — same role as
    /// `extent_sorted`, for the charge channel. Centred, not raw: centring
    /// and adding back `D·eps^2` is exactly equivalent to working on raw
    /// characters (verified to 1.66e-15 relative), and centred keeps
    /// [`ceiling`]'s documented algebra visible and both channels
    /// structurally identical in code.
    character_sorted: [f64; D],
    /// Antipodal-block statistics on the mean-centred characters, same
    /// packing as `extent_blocks`.
    character_blocks: [f64; D],
    /// Mean character — the charge channel's un-eliminated offset,
    /// `eps = mean_character(A) + mean_character(B)` (see [`ceiling`]'s own
    /// doc for why the shape channel has no equivalent term).
    mean_character: f64,
}

/// Mean-centred values from `raw`, sorted ascending via [`canonical_cmp`].
///
/// **No ID tie-break, unlike every other float sort this codebase's own
/// §13.4 doctrine requires.** Equal keys here are *exchangeable in the
/// consuming sum* — [`flat_bound`] and [`block_bound`] only ever read a
/// sorted array back through a dot-product-shaped accumulation, and
/// swapping two equal entries cannot change `Σ(x+y)²`. That argument, not
/// the fact that a stable sort happens to be deterministic, is what makes
/// omitting a tie-break legitimate here.
fn centred_sorted<const D: usize>(raw: &[f64; D], mean: f64) -> [f64; D] {
    let mut out: [f64; D] = core::array::from_fn(|i| {
        raw.get(i)
            .copied()
            .unwrap_or_else(|| unreachable!("i < D by from_fn's own contract"))
            - mean
    });
    out.sort_by(|a, b| canonical_cmp(*a, *b));
    out
}

/// Antipodal-block statistics on `centred` (already mean-centred): `[0..D/2)`
/// holds block sums, `[D/2..D)` holds `|block differences|`, both sorted
/// ascending. Each antipodal pair is visited exactly once, via `anti[i] > i`
/// selecting one canonical member — `s` is symmetric and `d` is absolute, so
/// neither depends on which member of the pair is picked as `i`.
///
/// **Stable sort, same reasoning as [`centred_sorted`].**
fn antipodal_blocks<const D: usize>(centred: &[f64; D], anti: &[u8; D]) -> [f64; D] {
    let mut sums: Vec<f64> = Vec::with_capacity(D / 2);
    let mut diffs: Vec<f64> = Vec::with_capacity(D / 2);
    for (i, &j) in anti.iter().enumerate() {
        let j = usize::from(j);
        if j > i {
            let ci = centred
                .get(i)
                .copied()
                .unwrap_or_else(|| unreachable!("i < D: anti has D entries"));
            let cj = centred
                .get(j)
                .copied()
                .unwrap_or_else(|| unreachable!("anti[i] < D by construction"));
            sums.push(ci + cj);
            diffs.push((ci - cj).abs());
        }
    }
    // Depends on `anti` having no fixed point (asserted where `anti` is
    // built, `geodesic.rs`'s own `Geodesic::build`) — a fixed point would
    // silently produce fewer than `D / 2` pairs here, leaving `out`'s tail
    // at its `0.0` init value rather than a computed one. Made explicit at
    // the point that relies on it, not just where the invariant is proven.
    debug_assert_eq!(
        sums.len(),
        D / 2,
        "anti should pair every direction exactly once"
    );
    sums.sort_by(|a, b| canonical_cmp(*a, *b));
    diffs.sort_by(|a, b| canonical_cmp(*a, *b));
    let mut out = [0.0f64; D];
    for (slot, v) in out.iter_mut().zip(sums.iter().chain(diffs.iter())) {
        *slot = *v;
    }
    out
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
    /// The rotation-invariant summary [`ceiling`] needs (§8.3, §8.6).
    ///
    /// Computed once per species at intern time, never per pair. **Takes
    /// `g`, unlike the 3-scalar summary this replaces** — the
    /// antipodal-block statistics [`ceiling`]'s block bound reads need
    /// [`Geodesic::anti`], and that table is universe-wide, not
    /// per-species, so it is passed in rather than rebuilt.
    #[must_use]
    pub fn summary(&self, g: &Geodesic<D>) -> SigSummary<D> {
        let (me, mc) = (mean_extent(self).get(), mean_character(self));

        let raw_extent: [f64; D] = core::array::from_fn(|i| {
            self.extents()
                .get(i)
                .copied()
                .unwrap_or_else(|| unreachable!("i < D by from_fn's own contract"))
                .get()
        });
        let raw_character: [f64; D] = core::array::from_fn(|i| {
            self.characters()
                .get(i)
                .copied()
                .unwrap_or_else(|| unreachable!("i < D by from_fn's own contract"))
        });
        let centred_extent: [f64; D] = core::array::from_fn(|i| {
            raw_extent
                .get(i)
                .copied()
                .unwrap_or_else(|| unreachable!("i < D by from_fn's own contract"))
                - me
        });
        let centred_character: [f64; D] = core::array::from_fn(|i| {
            raw_character
                .get(i)
                .copied()
                .unwrap_or_else(|| unreachable!("i < D by from_fn's own contract"))
                - mc
        });

        SigSummary {
            extent_sorted: centred_sorted(&raw_extent, me),
            extent_blocks: antipodal_blocks(&centred_extent, g.anti()),
            character_sorted: centred_sorted(&raw_character, mc),
            character_blocks: antipodal_blocks(&centred_character, g.anti()),
            mean_character: mc,
        }
    }
}

/// The rearrangement-inequality **flat** lower bound on
/// `min over ALL permutations of Σ(x[i] + y[σ(i)])²` — the sum-of-squares
/// itself, computed directly, never via the algebraically-equivalent
/// `Σx² + Σy² + 2Σxy` expansion.
///
/// **The expansion form is unsound, not merely imprecise, and this is
/// measured on the real corpus, not only on a constructed complement.**
/// `x_asc.iter().zip(y_asc.iter().rev())` pairs opposite-sorted — the
/// rearrangement inequality's minimising pairing over the full `D!`
/// permutation group, a strict superset of the 60 achievable contact
/// perms, so the result is a valid (if not the tightest achievable) lower
/// bound. The expansion form subtracts two large near-equal sums to reach
/// a small answer: measured, it exceeds the true minimum by up to
/// **+9.37e-16** on this crate's own random corpus and by **+2.13e-14** on
/// an exact complement fixture — six orders past
/// `a_physical_complement_scores_zero`'s `1e-20` absolute floor (not
/// `the_prefilter_is_a_genuine_bound`'s, which floors at `16.0 *
/// f64::EPSILON` and is a different, coarser check entirely). The
/// sum-of-squares form here is bit-identical to the true minimum on every
/// complement fixture checked.
#[must_use]
fn flat_bound<const D: usize>(x_asc: &[f64; D], y_asc: &[f64; D]) -> f64 {
    let mut acc = 0.0f64;
    for (x, y) in x_asc.iter().zip(y_asc.iter().rev()) {
        let t = x + y;
        acc += t * t;
    }
    acc
}

/// The rearrangement-inequality **antipodal-block** lower bound — a
/// relaxation over a strictly *smaller* set of permutations than
/// [`flat_bound`]'s "any permutation of `D` elements", but **not therefore
/// automatically a tighter bound** (the independent-ordering relaxation
/// below is a further loosening the subset argument alone does not cover).
/// It earns its place mainly on **hetero**-pairs; [`flat_bound`] is what
/// carries the self-pair diagonal — see [`SigSummary`]'s own doc for both
/// measurements, and [`ceiling`]'s own doc for why neither bound is
/// dropped.
///
/// [`Geodesic::contact_perms`] commutes with [`Geodesic::anti`] — verified
/// exhaustively at D = 12/42/162, every contact permutation maps antipodal
/// *pairs* to antipodal pairs. With `s = a1 + a2`, `d = a1 − a2` the two
/// members of one block:
///
/// ```text
/// (a1 + b1)² + (a2 + b2)² = ½[ (s_A + s_B)² + (d_A ∓ d_B)² ]
/// ```
///
/// — verified against all 60 real perms x 200 random vectors, worst
/// relative error 7.02e-16. Relaxing to any block permutation with a free
/// per-block sign (a superset of the achievable group, so still a valid
/// lower bound) and minimising each half independently: the `s` sum is
/// opposite-sorted (minimises `Σs_A s_B`, matching [`flat_bound`]); the `d`
/// sum is **same-sorted** (maximises `Σ|d_A||d_B|`, which minimises the
/// negated `(d_A − d_B)²` term — the one place the pairing direction
/// flips relative to `flat_bound`, and the thing most likely to be copied
/// backward from it).
#[must_use]
fn block_bound<const D: usize>(blocks_a: &[f64; D], blocks_b: &[f64; D]) -> f64 {
    let half = D / 2;
    let (sa, da) = blocks_a.split_at(half);
    let (sb, db) = blocks_b.split_at(half);
    let mut acc = 0.0f64;
    for ((s_a, s_b), (d_a, d_b)) in sa.iter().zip(sb.iter().rev()).zip(da.iter().zip(db.iter())) {
        let s = s_a + s_b;
        let d = d_a - d_b;
        acc += s * s + d * d;
    }
    0.5 * acc
}

/// A rotation-invariant **upper bound** on [`fit`]'s score (§8.3, Decision
/// 4 — issue #26, Task 26.1 Step 12).
///
/// §8.3 asks for "a genuine bound — a filter that cannot state the property it
/// guarantees is not conservative, it is merely untested". The property, stated:
///
/// > For every rotation `R`, `score(R) <= ceiling(A, B)`. So
/// > `ceiling(A, B) < threshold` proves `max_R score(R) < threshold`, and a
/// > rejection provably excludes no pair the full search would have accepted.
///
/// **Two lower bounds on the true minimum sum-of-squares per channel —
/// `flat_bound` and `block_bound` — and `ceiling` takes the *larger* of
/// the two, per channel, before negating.** Both are valid (the achievable
/// 60-rotation permutations are a subset of what each relaxes to), so the
/// larger of two valid lower bounds is still valid and is the tighter
/// ceiling. **Per channel, not per whole-bound** — each weight is drawn
/// strictly positive, so the sum of two per-channel maxima is `>=` the max
/// of the two whole-bound sums, and per-channel is strictly better.
/// Measured on `fixture_v2` seeds 6/17/42 at D=42 — 315 pairs, the corpus
/// `the_prefilter_is_a_genuine_bound` runs (V2 only; no test in this file
/// measures this on V1): neither bound is redundant, and neither
/// dominates. Block beats flat on 250 of 315 pairs in the shape channel
/// and 77 of 315 in the charge channel; flat beats block on 65 and 238
/// respectively. The two channels disagree about which bound wins on most
/// pairs, which is exactly why this is per-channel rather than a single
/// whole-bound comparison. See `the_prefilter_is_a_genuine_bound`'s own
/// `flat_wins_shape`/`block_wins_shape`/`flat_wins_charge`/
/// `block_wins_charge` counters, pinned exactly rather than left as prose,
/// so a future physics change cannot silently make this paragraph stale
/// again the way it did between `/review-pr` round 2 and round 3.
///
/// This replaces the predecessor's single Cauchy–Schwarz bound
/// (`(‖x‖ − ‖y‖)²`), which is **exactly zero on every self-pair** — the
/// centred channels are identical, so the spread difference vanishes
/// regardless of actual shape complexity — and which relied on the charge
/// channel's `D·eps^2` offset term to carry self-pair rejection at all.
/// Measured against real `PhysicsVersion::V2` data (genuinely
/// zero-straddling affinity, unlike V1's always-positive predecessor):
/// `worst_gap` under the old bound rises to ~1.0, effectively vacuous.
///
/// The charge channel keeps its `D·eps^2` term for the same reason as
/// before: centring both channels' sorted/block arrays makes the
/// permutation-dependent part of each bound a pure function of the
/// centred values, but the charge channel's uncentred offset
/// `eps = mean_character(A) + mean_character(B)` is permutation-invariant
/// and does not vanish the way the shape channel's derived separation
/// does — see the historical note below for why the shape channel has no
/// equivalent term.
///
/// **Historical note, preserved because the reasoning still applies.** The
/// shape channel contributes no `D·eps^2` term because `fit` separates the
/// bodies by `mean(a) + mean(b)`, which is exactly the value making the
/// offset zero — that is what "least-squares optimal" means. **Omitting a
/// term from a bound is always safe — it can only make the ceiling larger,
/// hence looser**; adding one needs proof, because it lowers the ceiling
/// and can push it below a score the search will actually find.
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
    let flat_shape = flat_bound(&a.extent_sorted, &b.extent_sorted);
    let block_shape = block_bound(&a.extent_blocks, &b.extent_blocks);
    // The comparison form, not `f64::max` — §13.1 bans the method for
    // ±0.0. The larger lower bound gives the tighter (smaller-magnitude)
    // ceiling once negated below.
    let l_shape = if flat_shape > block_shape {
        flat_shape
    } else {
        block_shape
    };

    let flat_charge = flat_bound(&a.character_sorted, &b.character_sorted);
    let block_charge = block_bound(&a.character_blocks, &b.character_blocks);
    let l_charge = if flat_charge > block_charge {
        flat_charge
    } else {
        block_charge
    };

    let eps = a.mean_character + b.mean_character;
    -(k.w_shape * l_shape + k.w_charge * (l_charge + d * eps * eps))
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
/// size at `|r| ≥ 0.987` *within a fixed atom count*. After: **0.583 / 0.510 /
/// 0.390** for seeds 6 / 17 / 42, and the spread across the 60 poses rises from
/// 7.5–9.7% of the score to 220–260%. The rotation search goes from a ripple to
/// the dominant term.
///
/// This said `0.43`–`0.47` until a reviewer noticed
/// `the_score_is_not_a_size_comparison` already **retracts that exact figure**
/// — it was measured on a corpus the fix commit then changed. Two numbers for
/// one measurement, disagreeing inside one file, in a module that quotes them
/// as measured maxima to justify later tightening.
///
/// **Not normalised by `L*²`.** Dividing through would make the shape term
/// dimensionless and the kernel scale-free, which is tempting since `w_shape`
/// currently multiplies a squared length while `w_charge` multiplies a pure
/// number.
///
/// **Held back deliberately — and the RAF literature *permits* this, it does
/// not *mandate* it.** An earlier version of this comment cited "the
/// RAF-specificity literature" for a rule against fully scale-free binding, in
/// one sentence and with no reference a reader could check; this replaces it
/// with what the cited papers actually say. RAF (Reflexively Autocatalytic and
/// Food-generated) theory does not require size-independent catalysis:
/// Mossel & Steel's foundational theorem already allows a molecule's catalytic
/// probability to "depend on the length of the molecule" (assumption R2), and
/// their Theorem 4.1(ii) gives the *same* lower bound — one independent of the
/// maximum length n — on the probability a RAF exists whether the catalysis
/// level is uniform across molecules (`μ_n(x) >= λn`) or proportional to each
/// molecule's own length (`μ_n(x) >= λγ_n|x|`, `γ_n ≈ 1`), for either
/// hypothesis requiring `λ > log_e(k)` (`k` the number of reaction types) —
/// a real precondition of the theorem, not one this doc invents, and one
/// both branches share equally so it does not favour either reading.
/// Hordijk, Wills &
/// Steel extend this to catalysis depending on the reaction too, and their
/// own extreme case (MLEN) — only maximum-length molecules catalyse at all —
/// still yields RAF sets, at a measured cost: the flat/uniform surrogate
/// mispredicts the required catalysis level by ~20% (their Fig. 3: 21.3% at
/// max length 8, falling to 18.3% at 13), against ~1% at max length 13 for a
/// template-matching model, and per-molecule catalysis variance rises from
/// `m` to `m + m²`. So a fully scale-free rule is not required either —
/// keeping `w_shape` un-normalised, so a monomer and a long polymer do not
/// compete on equal terms, is Borbax's own retained design choice, not
/// something the citation below settles for it.
///
/// **One distinction this literature is careful about and worth repeating
/// here, because collapsing it is the overclaim the wrong direction.** Every
/// "linear growth in catalysis level" result in this body of work bounds the
/// *probability a RAF exists at all*, over a random **ensemble** of reaction
/// networks, as a function of the average catalysis rate. None of it says how
/// sensitive one **specific** network's RAF closure is to a rate change —
/// closure for a given network is decided deterministically, by the
/// polynomial-time Hordijk–Steel algorithm computing the maximal RAF as a
/// fixed point. [`ceiling`] produces one specific score per molecule pair, not
/// an ensemble statistic, so the honest claim is narrower than "RAF closure is
/// sensitive to catalysis level": the size term here is a retained, reviewed
/// *carrier of catalytic specificity* for this simulation's own molecules, not
/// a result this literature proves about this kernel.
///
/// **The good news, and it is real.** This kernel had a structural defect —
/// issue #31 — that this same task's earlier steps fixed: `affinity`'s
/// attained range never straddled zero (fixed by V2's tanh-centred
/// derivation), and separately, the prefilter's own predecessor bound
/// (`(‖x‖ − ‖y‖)²`, a norm-difference form distinct from either bound
/// [`ceiling`] uses today) was exactly zero on every self-pair regardless of
/// shape — the discrimination-killing failure the rearrangement-inequality
/// **flat bound** now closes (see [`SigSummary`]'s own doc: flat carries
/// the self-pair diagonal, winning at least one channel on 41 of the 42
/// measured self-pairs, while the antipodal-**block** bound earns its
/// place mainly on hetero-pairs — 241 of 273). That is this codebase's own
/// version of the degenerate, size-only catalysis assignment Hordijk, Wills
/// & Steel study as the extreme worth quantifying rather than assuming
/// harmless — and fixing it moves this kernel from scoring some pairs on
/// size alone toward scoring every pair on shape, which is what makes the
/// retained size term above a
/// genuine *addition* to shape-driven specificity rather than a mask over its
/// absence.
///
/// - Hordijk, W., Wills, P. R. & Steel, M. "Autocatalytic Sets and Biological
///   Specificity." *Bull. Math. Biol.* 76(1):201–224, 2014.
///   `doi:10.1007/s11538-013-9916-4` (arXiv:1307.2860).
/// - Mossel, E. & Steel, M. "Random biochemical networks: the probability of
///   self-sustaining autocatalysis." *J. Theor. Biol.* 233(3):327–336, 2005.
///   `doi:10.1016/j.jtbi.2004.10.011`.
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
    use crate::testkit::{canon, chain_capable, fixture, fixture_v2, random_tree};
    use borbax_rng::{Domain, Stream};
    use borbax_universe::Universe;

    const D: usize = 42;

    fn geo() -> Geodesic<D> {
        Geodesic::<D>::build().unwrap_or_else(|_| unreachable!("D=42 builds"))
    }
    fn sig(mol: &Mol12, uni: &Universe, g: &Geodesic<D>) -> Signature<D> {
        signature(&embed(&canon(mol), uni), g)
    }

    /// Per-pair, per-channel comparison — which of `flat_bound`/`block_bound`
    /// wins each channel, separately. Extracted so
    /// `the_prefilter_is_a_genuine_bound` stays under clippy's line limit;
    /// what each of these four booleans feeds is documented at that
    /// function's own counters, not here.
    fn channel_wins(sa: &SigSummary<D>, sb: &SigSummary<D>) -> (bool, bool, bool, bool) {
        let flat_shape = flat_bound(&sa.extent_sorted, &sb.extent_sorted);
        let block_shape = block_bound(&sa.extent_blocks, &sb.extent_blocks);
        let flat_charge = flat_bound(&sa.character_sorted, &sb.character_sorted);
        let block_charge = block_bound(&sa.character_blocks, &sb.character_blocks);
        (
            flat_shape > block_shape,
            block_shape > flat_shape,
            flat_charge > block_charge,
            block_charge > flat_charge,
        )
    }

    /// A signature built by hand, so a complement can be constructed exactly.
    ///
    /// Test-only, and it must stay that way: a public constructor would let a
    /// caller mint a shape no molecule can produce and feed it to the kernel.
    fn from_parts(r: &[Span; D], a: &[f64; D]) -> Signature<D> {
        Signature::from_parts_for_test(*r, *a)
    }

    /// The constant both complement fixtures are built at, **shared so they
    /// cannot drift apart.**
    ///
    /// This is not deduplication for tidiness. The two fixtures are a *probe
    /// pair*: `complement_through_anti` must score 0 and
    /// `complement_same_index` must not, and the pair is what shows the `ANTI`
    /// composition is right. They are only comparable if they are built at the
    /// same level. With the computation written out twice, one could be edited
    /// and the probe would then compare two different shapes — still passing,
    /// having stopped discriminating. That is the same silent-inversion class
    /// the fixtures exist to catch, one level up.
    ///
    /// **Any** constant above `max(a.r)` works, because [`fit`] derives the
    /// separation from the pair rather than taking it from a constant (§8.3).
    fn complement_level(a: &Signature<D>) -> Span {
        let mut top = Span::ZERO;
        for e in a.extents() {
            if *e > top {
                top = *e;
            }
        }
        top + Span(1.0)
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
        let level = complement_level(a);
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
        let level = complement_level(a);
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
            hash, 0xb30f_55dc_b695_1731,
            "binding moved. If you meant to change the physics, regenerate this \
             constant and say so in the commit message. If you did not, suspect an \
             accumulation order: the two accumulators, the direction-index sum, and \
             mean_extent's `sum * inv_n` are each algebraically invariant under the \
             obvious tidy-up and each move the last bits. If only the POSE moved, \
             suspect the `>` tie-break: ties here are **exclusively self-pairs**, \
             which are automorphisms rather than float coincidences, and this \
             corpus has 9 of them per seed. **Regenerated deliberately for issue \
             #26, Task 26.1 Step 13 (2026-08-10): F2's per-species `shell_width` \
             (`signature.rs`) changes the character channel `fit` reads for every \
             molecule, in every universe, both physics versions — this is the \
             same deliberate change that moved `the_signature_digest_is_pinned`, \
             not an independent one. Previous value `0xab87_c6f7_2cf5_3428`.**"
        );
    }

    /// **A bit-level digest of `summary()` and `ceiling()` — the one this
    /// crate did not have.** `the_binding_digest_is_pinned` calls `fit()`
    /// only; every accumulation order inside `SigSummary`'s construction and
    /// `ceiling`'s bound was pinned by comment alone. Measured before this
    /// test existed: appending `.rev()` to the direction-index loop that
    /// used to build the old 3-scalar summary moved 105 of 120 species'
    /// `extent_spread` and the entire workspace stayed green, 0 failures.
    /// `ceiling` had no shipped caller before Task 12, so nothing reached a
    /// result; from here on it is the rejection decision, so its bits
    /// decide which pairs are even scored.
    #[test]
    fn the_prefilter_digest_is_pinned() {
        let g = geo();
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let mut bounds = 0u64;
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6800 + seed, Domain::Molecule, 0);
            let summaries: Vec<_> = (0..9)
                .map(|i| sig(&random_tree(&mut rng, 3 + (i % 9), &tbl, &ids), &uni, &g).summary(&g))
                .collect();
            for i in 0..summaries.len() {
                for j in i..summaries.len() {
                    let sa = summaries
                        .get(i)
                        .unwrap_or_else(|| unreachable!("i < summaries.len()"));
                    let sb = summaries
                        .get(j)
                        .unwrap_or_else(|| unreachable!("j < summaries.len()"));
                    hash ^= borbax_units::canonical_bits(ceiling(sa, sb, k));
                    hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                    bounds += 1;
                }
            }
        }
        // 3 seeds x 45 unordered pairs of 9 molecules.
        assert_eq!(bounds, 3 * 45, "the corpus did not run to completion");
        assert_eq!(
            hash, 0x160f_3037_444d_4858,
            "the prefilter bound moved. If you meant to change the physics, regenerate this \
             constant and say so in the commit message. If you did not, suspect an \
             accumulation order in `centred_sorted`, `antipodal_blocks`, `flat_bound` or \
             `block_bound` — each has a pinned direction (which array is read forward, which \
             is read `.rev()`) that is algebraically invariant under the obvious tidy-up and \
             moves the last bits or, in the block bound's case, the sign of a whole term."
        );
    }

    /// **`the_prefilter_digest_is_pinned`'s own sibling on `fixture_v2` —
    /// `/review-pr` round 3, finding F6.** `the_prefilter_digest_is_pinned`
    /// above runs on `fixture(seed)`, which uses `PhysicsVersion::CURRENT`
    /// (`V1` — `lib.rs`'s own doc). A V2-only physics change therefore moves
    /// nothing this file pins bit-for-bit: `sigma_deep` (round 2, landed in
    /// the same commit as a round-2 doc re-measurement of the V2-only
    /// `the_prefilter_is_a_genuine_bound`/`the_prefilter_stays_tight_at`
    /// corpus) shipped a genuinely-changed V2 prefilter with a green suite
    /// and nine stale figures in prose, undetected for a full review round.
    /// This test exists so the *next* V2-only physics change fails here
    /// first, with a message naming exactly what to re-check, rather than
    /// leaving prose the only thing that (eventually, silently) goes stale.
    #[test]
    fn the_prefilter_digest_is_pinned_under_v2() {
        let g = geo();
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let mut bounds = 0u64;
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture_v2(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6900 + seed, Domain::Molecule, 0);
            let summaries: Vec<_> = (0..9)
                .map(|i| sig(&random_tree(&mut rng, 3 + (i % 9), &tbl, &ids), &uni, &g).summary(&g))
                .collect();
            for i in 0..summaries.len() {
                for j in i..summaries.len() {
                    let sa = summaries
                        .get(i)
                        .unwrap_or_else(|| unreachable!("i < summaries.len()"));
                    let sb = summaries
                        .get(j)
                        .unwrap_or_else(|| unreachable!("j < summaries.len()"));
                    hash ^= borbax_units::canonical_bits(ceiling(sa, sb, k));
                    hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                    bounds += 1;
                }
            }
        }
        // 3 seeds x 45 unordered pairs of 9 molecules.
        assert_eq!(bounds, 3 * 45, "the corpus did not run to completion");
        assert_eq!(
            hash, 0x11d8_f8ea_6dfb_0ec9,
            "the V2 prefilter bound moved. If you meant to change the V2-only physics \
             (`orbital.rs`, `perturbation.rs`), regenerate this constant and say so in the \
             commit message — and re-check `the_prefilter_is_a_genuine_bound`'s and \
             `the_prefilter_stays_sound_and_tight_at_every_resolution`'s pinned counters and \
             bars in the same commit, since they read the same corpus shape this digest does. \
             First pinned 2026-08-10, `/review-pr` round 3, finding F6 — see \
             `the_prefilter_digest_is_pinned`'s own doc for why a V1-only digest could not \
             have caught the staleness that made this test necessary."
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
                // **The same guard its sibling twelve lines up already has.**
                // A NaN makes `rs > worst_shape` false, leaves `worst_shape` at
                // 0.0, and the assertion below — the one the doc calls "the
                // claim that actually rests on the fix" — passes having
                // measured nothing. `bare` is the *more* likely of the two to
                // degenerate, not the less: zeroing `w_charge` removes a whole
                // channel, so it has less variance than `scores`.
                assert!(
                    rs.is_finite(),
                    "seed {seed} n={n}: the shape-only correlation is not finite — that would \
                     leave `worst_shape` at 0.0 and pass the assertion below vacuously"
                );
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
    ///
    /// **`fixture_v2`, not `fixture` — this is the corpus Decision 4 exists
    /// for.** V1's always-positive affinity never exercises the regime that
    /// broke the predecessor bound: measured, the predecessor's `worst_gap`
    /// on V1 is `0.8417031120681644` (this file's own doc used to quote
    /// `0.93` here, which was never re-derived after the corpus changed —
    /// the actual number is asserted at the bottom of this test, not left in
    /// a comment for a future edit to go stale again) and rises to
    /// effectively `1.0` — vacuous — under real V2 affinity. A test that only
    /// ran on V1 would stay green through exactly the regression Decision 4
    /// fixes.
    ///
    /// **Soundness is checked over every pair including the diagonal, never
    /// relaxed.** Tightness is split in two: self-pairs (`i == j`) are
    /// structurally the hardest case for a rearrangement bound of any kind —
    /// a vector's opposite-sorted pairing against itself is close to an
    /// exact cancellation, and no achievable rotation can reproduce it, so
    /// the relaxation is loosest exactly there (measured: the shape
    /// channel's tightness ratio averages ~4x worse on self-pairs than on
    /// hetero-pairs). Giving them a separate, looser bar means a future
    /// regression that makes self-pair rejection vacuous again — exactly
    /// what broke under V1's charge-offset-only defence — cannot hide
    /// behind hetero-pairs' better average.
    ///
    /// **Four counter families, each pinned exactly rather than `> 0` —
    /// `/review-pr` round 3, finding F7.** `flat_wins_shape`/
    /// `block_wins_shape`/`flat_wins_charge`/`block_wins_charge` count which
    /// bound wins EACH channel, separately — the "probe the guard" check
    /// CLAUDE.md asks for on every guard, since neither bound is redundant
    /// in either channel and a mutation deleting either one would otherwise
    /// still pass soundness and a loose-enough tightness bar.
    /// `flat_beats_block`/`block_beats_flat` count which side won AT LEAST
    /// ONE channel (an OR across both, not a per-channel split).
    /// `flat_wins_self`/`block_wins_hetero` split by pair class — the
    /// counter an earlier version of this test did not have, and the reason
    /// `SigSummary`'s own doc had the two bounds' roles backwards for a full
    /// task's worth of history: `flat_beats_block`/`block_beats_flat` alone
    /// cannot tell whether a win was on the self-pair diagonal (the claimed
    /// job) or hetero-pairs (the real one). `ceiling`'s and `SigSummary`'s
    /// own docs quote all of these numbers; pinning them here exactly,
    /// rather than `> 0`, is what stops a physics change rotting those docs
    /// silently the way `sigma_deep` did to the previous (unpinned)
    /// measurement in `/review-pr` round 2 — undetected for a full round.
    /// These are exact integer counts over a fixed corpus, so there is
    /// nothing to round; `flat_wins_self`/`block_wins_hetero` are the guard
    /// that would have caught the reversed claim the day it shipped, and,
    /// pinned exactly rather than `> 0`, would also have caught
    /// `sigma_deep` moving the self-pair split from 42/0 to 41/1.
    ///
    /// **The bound must be tight enough to be worth having, and the worst-
    /// gap bars below are re-measured against what this branch actually
    /// ships.** `ceiling` is justified in its own doc by turning Task 20's
    /// D=162 sweep from hours into minutes, and a bound that never rejects
    /// delivers none of that. Measured 2026-08-10 against the physics this
    /// branch ships (`.abs()` -> clamp AND `sigma_deep`, both landed in the
    /// same `/review-pr` round-2 commit — an intermediate figure taken
    /// between the two, 0.9491563/0.9651054, was carried in a previous
    /// version of this comment for a full round without being re-measured
    /// against the second change): hetero 0.9144535, self 0.9592222. A
    /// ceiling of `0.0` gives 1.0 for every pair and fails both bars.
    #[test]
    fn the_prefilter_is_a_genuine_bound() {
        let g = geo();
        let (mut checked, mut worst_gap_hetero, mut worst_gap_self) = (0u32, 0.0f64, 0.0f64);
        let (mut flat_wins_shape, mut block_wins_shape) = (0u32, 0u32);
        let (mut flat_wins_charge, mut block_wins_charge) = (0u32, 0u32);
        let (mut flat_beats_block, mut block_beats_flat) = (0u32, 0u32);
        let (mut flat_wins_self, mut block_wins_hetero) = (0u32, 0u32);
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture_v2(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6600 + seed, Domain::Molecule, 0);
            let mols: Vec<_> = (0..14)
                .map(|i| sig(&random_tree(&mut rng, 3 + (i % 10), &tbl, &ids), &uni, &g))
                .collect();
            let summaries: Vec<_> = mols.iter().map(|m| m.summary(&g)).collect();
            for i in 0..mols.len() {
                for j in i..mols.len() {
                    let sa = summaries
                        .get(i)
                        .unwrap_or_else(|| unreachable!("i < mols.len()"));
                    let sb = summaries
                        .get(j)
                        .unwrap_or_else(|| unreachable!("j < mols.len()"));
                    let bound = ceiling(sa, sb, k);
                    let best = fit(&mols[i], &mols[j], &g, k).score();
                    // Absolute floor as well as relative: on an exact
                    // complement the relative term collapses to ~1e-39 while
                    // the float discrepancy is at the scale of the
                    // sum-of-squares form's own rounding. This random corpus
                    // does not happen to contain an exact complement pair —
                    // `a_physical_complement_scores_zero` is what actually
                    // exercises that fixture, with its own `1e-20` floor —
                    // so the floor here is a principled safeguard against a
                    // corner case this specific corpus has not hit, not a
                    // measured necessity. See `ceiling`'s doc. `16.0 *
                    // f64::EPSILON`, not a bare decimal literal — `EPSILON`
                    // is specifically the gap at 1.0, so this reads as
                    // exactly "16 ULP" only for magnitudes near 1; what
                    // actually matters is that the relative term above
                    // already dominates this floor everywhere it is not
                    // near-zero (`1e-9 * |bound| >= 16 * EPSILON` once
                    // `|bound| >= ~3.6e-6`), so the floor's own scale only
                    // needs to be *a* small, principled multiple of the
                    // portable unit — not tuned to any specific magnitude —
                    // for the regime it actually governs: `bound` near zero.
                    // Same reasoning `element.rs`'s own near-tie tolerance
                    // uses for choosing a small EPSILON multiple over a bare
                    // literal.
                    assert!(
                        best <= bound + 1e-9 * bound.abs() + 16.0 * f64::EPSILON,
                        "seed {seed} pair ({i},{j}): the search found {best}, above the \
                         ceiling {bound} — the filter would reject a binding pair"
                    );

                    // `u32::from(bool)`, not an `if` block per counter —
                    // this test's own doc explains what each of the six
                    // counters below means.
                    let (fs, bs, fc, bc) = channel_wins(sa, sb);
                    flat_wins_shape += u32::from(fs);
                    block_wins_shape += u32::from(bs);
                    flat_wins_charge += u32::from(fc);
                    block_wins_charge += u32::from(bc);
                    let flat_wins_channel = fs || fc;
                    flat_beats_block += u32::from(flat_wins_channel);
                    flat_wins_self += u32::from(flat_wins_channel && i == j);
                    let block_wins_channel = bs || bc;
                    block_beats_flat += u32::from(block_wins_channel);
                    block_wins_hetero += u32::from(block_wins_channel && i != j);

                    // The comparison form, not `f64::max` — §13.1. Both
                    // operands are finite here.
                    let scale = if best.abs() > 1.0 { best.abs() } else { 1.0 };
                    let gap = (bound - best) / scale;
                    if i == j {
                        if gap > worst_gap_self {
                            worst_gap_self = gap;
                        }
                    } else if gap > worst_gap_hetero {
                        worst_gap_hetero = gap;
                    }
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 3 * 105, "the corpus did not run to completion");
        // Exact counts, not `> 0` — see this function's own doc.
        assert_eq!(
            (
                flat_wins_shape,
                block_wins_shape,
                flat_wins_charge,
                block_wins_charge
            ),
            (65, 250, 238, 77),
            "the flat/block per-channel win split moved. If you changed the physics, \
             regenerate these AND update `ceiling`'s doc (its per-channel figures) and \
             `SigSummary`'s doc in the same commit."
        );
        assert_eq!(
            (flat_beats_block, block_beats_flat),
            (251, 263),
            "the combined (either-channel) win split moved — see `SigSummary`'s doc, which \
             quotes both of these."
        );
        // The specialised claim — see this function's own doc.
        assert_eq!(
            (flat_wins_self, block_wins_hetero),
            (41, 241),
            "the pair-class split moved — see `SigSummary`'s doc, which quotes both of these."
        );
        // Re-measured against the shipped physics — see this function's own doc.
        assert!(
            worst_gap_hetero < 0.93,
            "the worst hetero-pair ceiling-to-score gap is {worst_gap_hetero}, so the bound \
             is close to vacuous on the pairs that dominate a real corpus. (D=42 only — see \
             `the_prefilter_stays_sound_and_tight_at_every_resolution` for D=12/162.)"
        );
        assert!(
            worst_gap_self < 0.98,
            "the worst self-pair ceiling-to-score gap is {worst_gap_self} — self-pair \
             discrimination (Decision 4's own goal, and #10.1's membrane self-assembly \
             case) has regressed toward the predecessor's vacuous behaviour."
        );
    }

    /// `sig`'s generic sibling — the module's own `sig` is pinned to the
    /// module-level `const D: usize = 42`, so a D-sweep needs its own
    /// helper rather than a second copy of `geo`/`sig` per resolution.
    fn sig_at<const N: usize>(mol: &Mol12, uni: &Universe, g: &Geodesic<N>) -> Signature<N> {
        signature(&embed(&canon(mol), uni), g)
    }

    /// `the_prefilter_is_a_genuine_bound`'s own tightness bars are D=42
    /// only — finding F1 of the 2026-08-10 `/review-pr` pass, and this is
    /// the fix. `SigSummary`'s block-bound relaxes to `(D/2)! * 2^(D/2)`
    /// block permutations against 60 achievable, which loosens
    /// monotonically with `D`; §22.2 sweeps `D in {12, 42, 162}`, and the
    /// same fixture corpus's gaps widen with `D` on this same fixture
    /// corpus — measured directly per resolution below, not assumed to
    /// track D=42's bars. This function also asserts soundness at every
    /// resolution it runs — soundness (the bound never overshoots the true
    /// minimum) is unconditional and does not loosen with `D`, but nothing
    /// else in this file re-checks it away from D=42, so it is not
    /// redundant to assert here too.
    fn the_prefilter_stays_tight_at<const N: usize>(bar_hetero: f64, bar_self: f64) {
        let g = Geodesic::<N>::build().unwrap_or_else(|_| {
            unreachable!("D in {{12, 42, 162}} builds — see Geodesic::build's own doc")
        });
        let (mut worst_gap_hetero, mut worst_gap_self) = (0.0f64, 0.0f64);
        let mut checked = 0u32;
        for seed in [6u64, 17, 42] {
            let (tbl, uni) = fixture_v2(seed);
            let k = BindConsts::of(&uni.consts);
            let ids = chain_capable(&tbl);
            let mut rng = Stream::new(6600 + seed, Domain::Molecule, 0);
            let mols: Vec<_> = (0..14)
                .map(|i| sig_at::<N>(&random_tree(&mut rng, 3 + (i % 10), &tbl, &ids), &uni, &g))
                .collect();
            let summaries: Vec<_> = mols.iter().map(|m| m.summary(&g)).collect();
            for i in 0..mols.len() {
                for j in i..mols.len() {
                    let sa = summaries
                        .get(i)
                        .unwrap_or_else(|| unreachable!("i < mols.len()"));
                    let sb = summaries
                        .get(j)
                        .unwrap_or_else(|| unreachable!("j < mols.len()"));
                    let bound = ceiling(sa, sb, k);
                    let best = fit(&mols[i], &mols[j], &g, k).score();
                    // Same absolute floor as the D=42 sibling above, and the
                    // same reasoning: `16.0 * f64::EPSILON`, not `1e-15`.
                    assert!(
                        best <= bound + 1e-9 * bound.abs() + 16.0 * f64::EPSILON,
                        "D={N} seed {seed} pair ({i},{j}): the search found {best}, above the \
                         ceiling {bound} — the filter would reject a binding pair"
                    );
                    let scale = if best.abs() > 1.0 { best.abs() } else { 1.0 };
                    let gap = (bound - best) / scale;
                    if i == j {
                        if gap > worst_gap_self {
                            worst_gap_self = gap;
                        }
                    } else if gap > worst_gap_hetero {
                        worst_gap_hetero = gap;
                    }
                    checked += 1;
                }
            }
        }
        assert_eq!(
            checked,
            3 * 105,
            "D={N}: the corpus did not run to completion"
        );
        assert!(
            worst_gap_hetero < bar_hetero,
            "D={N}: the worst hetero-pair ceiling-to-score gap is {worst_gap_hetero}, over the \
             {bar_hetero} bar"
        );
        assert!(
            worst_gap_self < bar_self,
            "D={N}: the worst self-pair ceiling-to-score gap is {worst_gap_self}, over the \
             {bar_self} bar"
        );
    }

    /// The bars widen with `D`, measured, not guessed: `SigSummary`'s
    /// block-bound relaxation loosens as `D` grows (this function's own
    /// doc), so a single bar tight enough for D=12 would already fail at
    /// D=42, and D=42's own bar fails at D=162 — the exact gap finding F1
    /// found. `ceiling`'s own doc cites the D=162 case specifically
    /// (Task 20's real sweep resolution), so that is the leg with the
    /// least room in its own bar, not the most.
    #[test]
    fn the_prefilter_stays_sound_and_tight_at_every_resolution() {
        // Measured 2026-08-10 against the physics this branch actually
        // ships (`.abs()` -> clamp AND `sigma_deep`, both landed in the same
        // `/review-pr` round-2 commit — an earlier version of this comment
        // carried a figure taken between the two, never re-measured against
        // the second change, the same staleness
        // `the_prefilter_is_a_genuine_bound`'s own bar carried): D=12 hetero
        // 0.8351392 / self 0.8954218; D=162 hetero 0.9312796 / self
        // 0.9741071. D=12's bars were never touched by that staleness (they
        // already had comfortable margin against the real number); D=162's
        // hetero bar was loosened on the strength of the stale figure and is
        // retightened here to the same proportion of margin D=42's own bar
        // (`the_prefilter_is_a_genuine_bound`) carries.
        the_prefilter_stays_tight_at::<12>(0.90, 0.97);
        the_prefilter_stays_tight_at::<162>(0.94, 0.98);
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
            //
            // Absolute floor as well as relative, matching the two sibling
            // assertions in this file. A purely relative tolerance is exactly
            // zero at a zero score, which demands bit equality and fails on a
            // *correct* kernel — and `ceiling`'s doc records that scores reach
            // `-0.0` on complement shapes. Leaving the third one bare would
            // make the convention look accidental.
            assert!(
                (f.score() - direct.score()).abs() <= 1e-9 * direct.score().abs() + 1e-20,
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
