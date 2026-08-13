//! Task 26.2's nuclear structure model (issue #26) — the SEMF-style
//! three-term per-unit binding energy, [`crate::nuclear::composition_argmax`]'s closed form
//! (Decision 11 of `docs/superpowers/plans/2026-08-07-borbax-real-physics.md`),
//! and Step 1's two-armed empirical gate.
//!
//! No production caller exists yet — Steps 2-7 wire this into [`crate::element`].
//! This module's own job is narrower: settle, with real numbers, whether the
//! mechanism produces a valley of stability at all before any of that wiring
//! happens. "If either arm fails, stop — do not proceed to wire this into
//! `Element`" (Step 1, routed requirement 4).
//!
//! # Three corrections to the plan's own text, made here rather than silently
//!
//! Landed after a three-amigos design session and two independent verification
//! passes (`geometry-numerics-reviewer`, a `physics-plausibility-reviewer`
//! dispatch) — both dispatched because this project has shipped confident,
//! internally-consistent, wrong derivations before (CLAUDE.md's own account of
//! Task 2, and this plan's own round-3/4/5 corrections of round-2's math).
//!
//! **1. [`crate::nuclear::composition_argmax`] uses the κ-free spelling, not the one Decision
//! 11 names "the correct closed form."** Decision 11 gives two algebraically
//! identical spellings of the same closed form — `(4κt + γt^(2/3)) / (8κ +
//! 2γt^(2/3))` (which it calls correct) and `(2t + c·t^(2/3)) / (4 +
//! 2c·t^(2/3))` (introduced only to show κ cancels, explicitly "not as a
//! second spelling to implement"). They are provably equal algebraically
//! (both numerator and denominator of the first factor by exactly `2κ`), but
//! **not** numerically equal in `f64`: measured over 6,951,474 sampled
//! `(κ, c, t)` cells, they differ by up to 4 ulp, and in 3 of those cells —
//! all at `c ≈ 0.0192, t = 125` — they round to a **different integer**,
//! i.e. a different element. Exactly one spelling can be shipped, and the
//! κ-free one is the right one: it has no `κ` in it at all, so a perturbation
//! of `κ` provably cannot move any isotope assignment — a physics property,
//! not just a numerical preference.
//!
//! **2. [`crate::nuclear::BASE_EPS`] is `2.625`, not `1.0`.** `1.0` was V1's own historical
//! `eps` midpoint, chosen when every cost coefficient in this expression was
//! O(0.01-0.1); Task 26.2's pre-registration pinned two of the three cost
//! coefficients (`kappa`, `c`) to real MeV-scale SEMF values and left `eps`
//! at its old, unrelated scale. Consequence, measured at the identity
//! coefficients and the shipped [`crate::nuclear::NUCLEAR_K`]` = 10`: 305 of 386
//! totals (79.0%) have **negative** on-valley per-unit binding energy — the
//! model is unbound over most of its domain, not merely peaked low. (An
//! earlier draft of this correction quoted 288/386, 74.6% — measured at the
//! superseded `NUCLEAR_K = 12` and never re-taken after correction 3 below;
//! the direction is unaffected, 79% is worse, but see
//! `the_identity_configuration_is_bound_over_all_but_two_totals` for the
//! figure as a pinned test rather than a paragraph.) Corrected value:
//! `eps = 2*a_V/z`, `a_V = 15.75` MeV (Rohlf's real SEMF volume coefficient —
//! the same citation `kappa` already uses), `z = 12` (this crate's
//! `packing::shell_size(k, 1)` at `NUCLEAR_K = 10`) — matching the real
//! limit `eps * contacts(t)/t -> a_V` as `t -> infinity` (a `/review-pr`
//! correction: an earlier draft named this `BE/A -> a_V`, which is false —
//! the real `BE/A` along the valley asymptotes to `a_V - a_A = -7.95`
//! MeV/nucleon, not `a_V`; `a_V` alone is the *volume term's own*
//! per-nucleon limit, which is what `eps * contacts(t)/t` computes, not
//! `BE/A` overall — see this module's own "Post-gate correction" for the
//! correct `BE/A` asymptote), a limit this model approaches
//! slowly and does not reach within `T_MAX` (`eps * contacts(T_MAX)/T_MAX`
//! is `13.19`, not `15.75`, at `T_MAX = 386`). This is a genuine
//! pre-registration correction, landing the same way Decision 11's own
//! round-5 correction to `base_c` did: found while verifying, not tuned to
//! Step 1's own output.
//!
//! **3. [`crate::nuclear::NUCLEAR_K`] is `10`, not `12`.** `contacts(total)` has no source
//! under the V2 electronic-configuration model (`generate_elements_v2` never
//! calls [`crate::packing`] at all), so Task 26.2 needs its own fixed
//! coordination parameter — fixed, not drawn, so the nuclear axis stays
//! independent of any table-level draw (Step 3's independence test). The
//! coordination number this crate's own shell law delivers is `z = k + 2`
//! (`packing::shell_size`'s own doc), so the FCC/icosahedral
//! coordination-12 argument that motivated fixing `k` at all actually
//! selects `k = 10`, not `12` — `k = 12` gives `z = 14`.
//!
//! # Post-gate correction, 2026-08-13: the model family is three coefficients, not four
//!
//! Found after Step 1's gate already returned "both arms pass" — by a
//! `/review-pr` follow-up question, not by verification before the gate ran,
//! unlike corrections 1-3 above. Recorded under its own heading rather than
//! folded into that list, so its own claim ("found while verifying, not
//! tuned to Step 1's own output") is not quietly extended to cover a change
//! that was made after reading the output.
//!
//! **`sigma * total^(5/3)`, and the constant that fed it, are deleted.** A
//! `physics-plausibility-reviewer` pass found the term structurally
//! unjustified, not merely miscalibrated: `eps * contacts(total)` already
//! reproduces the real SEMF's volume-and-surface behaviour, so the three
//! surviving coefficients (`eps`, `kappa`, `c`) already form the real
//! four-term SEMF at Rohlf's own values — volume, surface (emergent from
//! the packing series, not a separate term), asymmetry, Coulomb.
//! **Precisely, not "1-9% at every total" (a `/review-pr` correction — the
//! original figure was measured over the wrong range and conflated two
//! different quantities): over `t in [3, 386]`, the deficit
//! `a_V - eps*contacts(t)/t` fits `19.753 * t^-0.3445` against the real
//! surface term's own `17.8 * A^-1/3`** (exponent within 3.4% — reproduces
//! from the unrounded `0.3445`, not the rounded `0.345`, which gives 3.5%
//! — coefficient within 11%, confirming the packing series generates a
//! surface term with the right *scaling law*; the fit range matters —
//! `[1,386]` gives `19.35 * t^-0.3406`, "within 2.2%", and `[50,386]`
//! gives `19.80 * t^-0.3451`, "within 3.5%")
//! **— the combined volume-plus-surface reproduction itself is within
//! 2.63% for `t >= 50`, within 1% for `t >= 266`, degrading to 34.6% at
//! `t = 3`, where this model — like the real SEMF — has no pairing term.**
//! `sigma * total^(5/3)` was a fifth term duplicating a role already
//! filled — structurally a second, weaker Coulomb-shaped cost (its own
//! per-unit form, `sigma * t^(2/3)`, is the one *growing* shape shared
//! with a real SEMF term, and that term is Coulomb's, not surface's), not
//! a calibration problem the model's own structure already accommodated.
//!
//! **The "growth limit" argument that had justified keeping it does not
//! hold.** `docs/superpowers/plans/2026-08-07-borbax-real-physics.md`'s own
//! proof (~line 1785) established `sigma` as the model's only unbounded
//! term and concluded it was therefore required to supply a growth limit.
//! The premise is correct and the conclusion does not follow: a *bounded*
//! cost that exceeds a *saturating* gain limits growth without needing to
//! be unbounded itself, and `eps`/`kappa`'s own per-unit saturation already
//! does exactly that — the on-valley per-unit asymptote at `sigma = 0` is
//! `eps * z/2 - kappa`, which at the identity coefficients is `-7.95`
//! MeV/nucleon, **identical to** the real SEMF's own asymptote (`a_V -
//! a_A`) — by construction, not independent corroboration (`BASE_EPS` is
//! *defined* as `2*a_V/z` and `BASE_KAPPA` *is* `a_A`, so `eps*z/2 - kappa
//! == a_V - a_A` algebraically for any `z`). The substantive, actually
//! checkable claim is that `contacts(total)/total` really does converge to
//! `z/2` — verified, though slowly (`5.303` at `total=1000`, `5.684` at
//! `10^4`, `5.854` at `10^5`, against `z/2 = 6`). `sigma` was never
//! load-bearing for the model having *a* growth
//! limit — only for exactly where it fell, and where it fell was wrong:
//! at the shipped `sigma = 0.075`, the identity zero-crossing sat at
//! `total = 603`; without it, `total = 3055`, close to the real SEMF's own
//! `~3076` under the same `Z(Z-1)` Coulomb convention this model uses.
//!
//! **Deliberately not framed as a coverage-improvement fix.** Removing
//! `sigma` measurably improves Step 1's own coverage-while-bound statistic
//! (see [`crate::perturbation::MigratedConstant::C`]'s doc: 9.25% of
//! universes covering fewer than 120 elements while bound falls to 1.09%),
//! but that is a consequence of the correction, not its justification —
//! the justification is that the deleted term modelled a physical regime
//! this system does not have. Keeping the two separate matters: the next
//! coefficient that happens to improve a coverage statistic is not thereby
//! justified for removal on that basis alone.
//!
//! **Re-measured, not assumed unaffected.** Both arms were re-run over the
//! same 38,416-seed corpus. Arm 1 (drift band `[Some(16), Some(30)]`) is
//! unaffected, as expected — `sigma` never entered [`crate::nuclear::composition_argmax`].
//! Arm 2's closure moved from `38416/38416` to `38406/38416`
//! (`0.999740`) — `sigma` was, per the physics reviewer's own finding, a
//! meaningful contributor to positive shedding-`Q` at `T_MAX` for a small
//! number of universes, and losing it costs ten of them; the closure bar
//! (`0.894404`) still clears with wide margin. **A `/review-pr` follow-up
//! instrumented the split `bound_shape`/`shedding_q` alone don't
//! distinguish: all ten losses fail only the Q-value half, none the
//! interior-peak half — but the interior-peak half is passing on a
//! lattice artefact for a wider set than the ten that actually fail.**
//! [`crate::packing`]'s shell closures at `NUCLEAR_K = 10` are `1, 13, 55,
//! 147, 309`; without `sigma`'s growth-limiting term, a high-`eps`/low-`kappa`
//! universe's on-valley curve can still be rising at `T_MAX`, so the peak
//! walks up to the *last* closure rather than a real turnover from term
//! competition — 65/38,416 universes (0.17%) now peak exactly at `309`,
//! and only 10 of those 65 also fail positive-Q; the other 55 pass Arm 2
//! on margin in the same regime. Not a defect here — a coverage question
//! for Steps 2+, which own whether "interior peak" should mean more than
//! "not at the literal boundary." [`crate::nuclear::bound_shape`]'s own gate
//! still holds `38416/38416`, but the shape of what it is holding changed:
//! the last-bound-total range narrowed from `[62, 386]` to `[92, 386]`
//! (the worst universe now stays bound 30 totals further in), and
//! `98.6%` of universes (37,863/38,416) are now never unbound anywhere in
//! `4..=T_MAX` at all (most of what `bound_shape` used to discriminate
//! against no longer occurs) — the transition-count distribution that sets
//! `MAX_TRANSITIONS = 5` moved from `p50=0, p90=1, max=3` to `p50=0, p90=0,
//! max=3`, so the bound is unchanged but tighter relative to what is
//! actually observed. The `total = 309` shell-closure "island of
//! stability" (see [`crate::nuclear::bound_shape`]'s own doc) still occurs — two seeds
//! (3992, 4425) hit it in this corpus, down from three (14355, 20107,
//! 21217) before, since fewer universes reach an unbound region at all for
//! the blip to interrupt.
//!
//! # Where the numbers are recorded
//!
//! `docs/experiments/2026-08-12-nuclear-valley-gate.md` carries the full
//! corpus-sweep measurements and the regime they were taken under, updated
//! 2026-08-13 for the post-gate correction above. This module's own tests
//! carry the pass/fail verdicts and the permanent regression anchors.

#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "no production caller exists yet -- Steps 2-7 wire this into \
                  crate::element, and Step 1 (this module's own tests) is the first \
                  consumer of any of it"
    )
)]

use crate::packing::{self, PackingConsts};
use borbax_units::det_math;

/// Task 26.2's nuclear binding-energy scale — see this module's own doc,
/// correction 2, for why this is `2.625` rather than the pre-registration
/// commit's original `1.0`.
pub(crate) const BASE_EPS: f64 = 2.625;

/// Task 26.2's nuclear asymmetry coefficient, Rohlf's real `a_sym` (`MeV`).
pub(crate) const BASE_KAPPA: f64 = 23.7;

/// Task 26.2's dimensionless Coulomb/asymmetry ratio, `a_C/(2*a_sym)` at
/// Rohlf's real values (`0.711/(2*23.7)`). `gamma = 2*kappa*c` is always
/// derived from this, never independently drawn (Decision 11).
pub(crate) const BASE_C: f64 = 0.015;

/// The fixed design ceiling on `total = a + b`, an integer nucleon count
/// (Decision 11's round-5 pre-registration). Not itself P7-migrated.
pub(crate) const T_MAX: u32 = 386;

/// The fixed coordination parameter `contacts(total)` is built from — see
/// this module's own doc, correction 3, for why `10` and not `12`.
pub(crate) const NUCLEAR_K: usize = 10;

/// The three nuclear coefficients, always drawn and threaded together
/// (routed requirements 5/6: one shared rung, never independent draws) —
/// introduced 2026-08-13 during `/review-pr`'s round 3 on the "Post-gate
/// correction" above, once that correction's own diff made the case
/// concrete rather than speculative: of the nine function signatures it
/// touched, seven had zero logic change, purely a `- sigma: f64,`/
/// `- sigma,` edit repeated because these three values travel everywhere
/// together. Two prior review rounds recommended this same struct and
/// deferred it pending real production call sites (Steps 2-7); this
/// round's own diff is a measured demonstration that the churn is real
/// and entirely internal to this module, independent of what those call
/// sites end up wanting — so waiting no longer bought anything.
///
/// **Deliberately not adopted by [`composition_argmax`]/
/// [`composition_argmax_int`]/[`drift_onset`], which keep `c: f64`
/// alone** — that narrower signature statically encodes "`kappa` cancels
/// from the argmax, by construction" (also stated in prose by this
/// module's own doc and by `the_valley_is_a_one_parameter_family_in_c`,
/// this module's own test); folding `eps`/`kappa` back into these
/// signatures would smuggle them into a computation that structurally
/// does not need them, hiding the invariant the narrower signature
/// currently makes visible.
///
/// **No `PartialEq` derive, deliberately — a `/review-pr` finding
/// (`rust-developer-expert`, sharpened by `determinism-auditor`).** A
/// derived `PartialEq` on `f64` fields is invisible to `clippy::float_cmp`
/// (verified: the equivalent hand-written `a == b` on two bare `f64`s is
/// caught, the derive-generated one is not), so it would launder a float
/// comparison in a `pub(crate)` type past the exact lint this workspace
/// relies on to flag one. `the_three_nuclear_coefficients_share_one_rung`
/// (this module's own test) is the only place that needs equality — see
/// `NuclearCoeffs::is_identity_config` below, defined only under
/// `#[cfg(test)]` so a future production caller has no `==` to reach for
/// at all and has to write out what it actually means (exact match on
/// every field,
/// component-wise `.abs() < epsilon`, or something else).
#[derive(Debug, Clone, Copy)]
pub(crate) struct NuclearCoeffs {
    pub(crate) eps: f64,
    pub(crate) kappa: f64,
    pub(crate) c: f64,
}

#[cfg(test)]
impl NuclearCoeffs {
    /// Exact IEEE-754 equality against another `NuclearCoeffs` — not
    /// bit-exact (`-0.0 == 0.0` is `true` despite differing bit patterns)
    /// and not epsilon-tolerant (`NaN == NaN` is `false`). Correct for
    /// this module's one use, detecting rung 0 in
    /// `the_three_nuclear_coefficients_share_one_rung`:
    /// `draw_symmetric` returns `base` bit-for-bit there (`perturb`
    /// multiplies by exactly `1.0 + 0.0`), so exact equality is the right
    /// check for that specific comparison. Not a general "close enough"
    /// check — comparing two arbitrarily-drawn `NuclearCoeffs` with this
    /// would silently answer a question nobody asked.
    #[must_use]
    #[expect(
        clippy::float_cmp,
        reason = "the one deliberate exact-equality comparison this module makes -- see this \
                  method's own doc for why exact equality is correct here specifically"
    )]
    fn is_identity_config(self, base: Self) -> bool {
        self.eps == base.eps && self.kappa == base.kappa && self.c == base.c
    }
}

/// Contacts made by a cluster of `total` nucleons, via this crate's own
/// packing model at the fixed [`NUCLEAR_K`] (never a drawn `k` — see this
/// module's doc, correction 3, for why fixing it is required, not merely
/// convenient). **Extrapolates beyond `packing::FRONTIER_COEFF`'s own
/// calibration domain at `total` near `T_MAX`** — the constant is fitted
/// against exact counts at caps 12/42/162 (shell coordination 1-3), and
/// `total = T_MAX = 386` at `NUCLEAR_K = 10` reaches shell 5 (cap 252);
/// `packing.rs`'s own per-cap fit values (rising toward an asymptote as
/// cap grows) mean this is a genuine but bounded extrapolation, not an
/// error — a `/review-pr` finding (`geometry-numerics-reviewer`),
/// recorded here since it sits directly under this module's own headline
/// claim that `eps * contacts(total)` reproduces a real SEMF term.
#[must_use]
#[expect(
    clippy::as_conversions,
    reason = "total <= T_MAX = 386, far below usize's range on every supported target"
)]
pub(crate) fn contacts(total: u32) -> f64 {
    packing::contacts_upto(PackingConsts::new(NUCLEAR_K), total as usize)
}

/// The composition argmax at fixed `total`, Decision 11's closed form — the
/// κ-free spelling only, per this module's doc, correction 1.
///
/// `total^(2/3)` as `cbrt(total*total)`, matching `element.rs:627`'s house
/// convention (Task 26.2's own "every float spelling pinned" requirement).
/// **Not shared with [`total_binding_energy`] any more** — that function
/// computed the identical `t^(2/3)` value once, for `sigma`'s term; the
/// "Post-gate correction" (this module's own top-level doc) deleted it,
/// so this is now the only `t^(2/3)` computation in the module.
#[must_use]
pub(crate) fn composition_argmax(c: f64, total: u32) -> f64 {
    let t = f64::from(total);
    let t23 = det_math::cbrt(t * t);
    (2.0 * t + c * t23) / (4.0 + 2.0 * c * t23)
}

/// [`composition_argmax`], rounded to the nearest legal proton count.
///
/// **`round_ties_even`, not `f64::round`** — this crate's own house
/// convention with no exceptions (`Mass::new`, every mass-defect and
/// valence rounding in `element.rs`) — and **`a >= 1` is a domain
/// restriction, not a tie-break.** `composition_argmax(c, 1) == 0.5`
/// exactly for every legal `c` (`2*c` is an exact scaling, so `fl(4+2c) ==
/// 2*fl(2+c)` exactly and the quotient is exactly `0.5`), so
/// `round_ties_even` gives `a = 0` — a zero-proton "element" — in **every**
/// universe, at `total = 1`, unless this floor is applied. `f64::max` is
/// disallowed in this crate (`element.rs`'s own comment on
/// `find_peak_and_set_instability` records why); this is the explicit
/// conditional it uses instead.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    reason = "the floor at 1.0 makes this non-negative, and composition_argmax(c, total) <= \
              total <= T_MAX = 386, far below u32's range -- precondition c > 0, which every \
              call site satisfies (draw_symmetric scales a positive base by a multiplier in \
              [0.5, 1.5]) but which neither function takes as a runtime guard"
)]
pub(crate) fn composition_argmax_int(c: f64, total: u32) -> u32 {
    let rounded = composition_argmax(c, total).round_ties_even();
    if rounded < 1.0 { 1 } else { rounded as u32 }
}

/// The three-term **total** binding energy — the mechanism's own diagram,
/// stated as one signed expression; [`per_unit_binding_energy`] divides
/// this once by `total` to reach the per-unit form Task 26.2's own "state
/// the model in per-unit form" instruction actually asks for. Exactly
/// Rohlf's real SEMF's volume, asymmetry and Coulomb terms at
/// `eps`/`kappa`/`c`'s real citations. (This module's own doc, "Post-gate
/// correction," records why a fourth, `sigma * total^(5/3)`, size-cost
/// term was deleted rather than kept alongside these three: `eps *
/// contacts(total)` already supplies the real SEMF's surface term as an
/// emergent property of the packing series, so a separate size-cost term
/// duplicated a role already filled.)
///
/// `a` is the proton count, `b` the neutron count — not interchangeable:
/// the Coulomb/asymmetry term (`gamma * a*(a-1) / t13`) reads `a`
/// specifically, so a transposed call computes a different, still-plausible
/// number.
///
/// Every float spelling pinned per Task 26.2's own list: `total^(1/3)` as
/// `cbrt(total)`; `(a-b)^2` as `let d = a as f64 - b as f64; d*d`
/// (cast before subtraction, never on the unsigned difference); `a*(a-1)`
/// as `let af = a as f64; af*(af - 1.0)` (cast-first, never on the unsigned
/// product, which underflows at `a = 0` under `overflow-checks`);
/// left-to-right term association (`the_three_terms_are_summed_left_to_right`,
/// this module's own test, is the only test that can see this); the
/// division by `total` applied to the whole sum, not per term.
#[must_use]
pub(crate) fn total_binding_energy(a: u32, b: u32, coeffs: NuclearCoeffs) -> f64 {
    let NuclearCoeffs { eps, kappa, c } = coeffs;
    let total = a + b;
    let t = f64::from(total);
    let gamma = 2.0 * kappa * c;
    let t13 = det_math::cbrt(t);
    let d = f64::from(a) - f64::from(b);
    let af = f64::from(a);
    eps * contacts(total) - kappa * (d * d) / t - gamma * (af * (af - 1.0)) / t13
}

/// [`total_binding_energy`], divided once more by `total` — the per-unit
/// form the mechanism's own diagram is stated against.
#[must_use]
pub(crate) fn per_unit_binding_energy(a: u32, b: u32, coeffs: NuclearCoeffs) -> f64 {
    total_binding_energy(a, b, coeffs) / f64::from(a + b)
}

/// The on-valley total energy at a total — `total_binding_energy` at the
/// composition-changing channel's own argmax, [`composition_argmax_int`].
#[must_use]
pub(crate) fn on_valley_total_energy(total: u32, coeffs: NuclearCoeffs) -> f64 {
    let a = composition_argmax_int(coeffs.c, total);
    total_binding_energy(a, total - a, coeffs)
}

/// The on-valley per-unit energy at a total.
#[must_use]
pub(crate) fn on_valley_per_unit_energy(total: u32, coeffs: NuclearCoeffs) -> f64 {
    on_valley_total_energy(total, coeffs) / f64::from(total)
}

/// The smallest even `total` at which [`composition_argmax_int`] first
/// departs from the symmetric composition (`total / 2`) — Step 1 Arm 1's
/// scalar, per the plan's own round-4 correction ("gate on drift, not
/// interiority"). Only even totals are scanned, because "the symmetric
/// composition" is only a real thing to depart from where `total / 2` is
/// exact — at odd `total`, `total / 2` (integer division) is already the
/// floor of the true midpoint, a different and less meaningful comparison.
/// **Scanning every total instead of only even ones does not change the
/// answer** (verified at every tested `c`, including the reachable band's
/// endpoints and the identity value: the full scan's onset is identical to
/// the even-only scan's) — the correct reason to prefer even totals is
/// that they ask the physically meaningful question, not that odd totals
/// are somehow already trivially past onset (an earlier version of this
/// doc claimed the odd-total condition holds "from `total = 3` on",
/// unverified and wrong: the first odd total it actually holds at, for the
/// identity `c`, is `31`, well past the even-only onset of `22`).
///
/// `None` if no drift is ever found in `2..=T_MAX` — not expected to occur
/// for any legal `c` (`a* <= t/2` for every `t`, proved on paper and
/// pinned by `the_valley_never_reaches_past_half_of_total`, so drift is
/// not merely likely but structurally bounded within `T_MAX`).
/// `Option<u32>`, not a magic out-of-range sentinel: an earlier
/// version returned `T_MAX + 2`, and a caller comparing two such sentinels
/// against each other (as a range-containment check would) can pass
/// vacuously — the same shape of hazard `bound_shape`'s `Option<u32>`
/// already avoids.
#[must_use]
pub(crate) fn drift_onset(c: f64) -> Option<u32> {
    let mut total = 2_u32;
    while total <= T_MAX {
        if composition_argmax_int(c, total) < total / 2 {
            return Some(total);
        }
        total += 2;
    }
    None
}

/// Step-1-only definition of the size-shedding channel's Q-value: does
/// splitting an on-valley `total` into two on-valley fragments release
/// energy. **Scoped to this gate alone** — the plan's own text defines no
/// Q-value for this channel (it names only "V1's mechanism, kept" without
/// restating what that means under the corrected model), and Steps 5-7 own
/// the real per-isotope definition this module does not attempt.
///
/// **Does not conserve proton number** — both fragments are scored at
/// *their own* valley composition, not at a split of the parent's actual
/// `a`, so the returned Q is closer to a fission-plus-beta-chain total
/// than a single prompt-fission Q (measured at identity, `total = 238`:
/// parent `a = 93`, best-split fragments `93 + 101`, a `+8` proton
/// mismatch — magnitude-consistent with that reading, ~226 `MeV` against
/// real U-238's ~200 `MeV` prompt / ~215 `MeV` total). `Mass` conservation
/// (a hard invariant) is not violated anywhere in the crate by this — no
/// `Mass` value is constructed here — but a future caller promoting this
/// convention into a real reaction channel needs to know which quantity
/// it's actually computing.
#[must_use]
pub(crate) fn shedding_q(total: u32, coeffs: NuclearCoeffs) -> f64 {
    let whole = on_valley_total_energy(total, coeffs);
    (1..=total / 2)
        .map(|s| {
            on_valley_total_energy(total - s, coeffs) + on_valley_total_energy(s, coeffs) - whole
        })
        .fold(f64::MIN, |best, q| if q > best { q } else { best })
}

/// The total at which the on-valley per-unit binding energy peaks —
/// first-wins on ties by iteration order, the same idiom
/// `element.rs::find_peak_and_set_instability` already uses for the
/// identical shape of search (strict `>`, no `f64::max`). **The tie-break
/// matters here and is not free to change**: this fold returns the
/// *index* (`total`) the maximum occurred at, not just the maximum value
/// — a genuinely different total can tie on energy, and which one wins is
/// exactly "first, by iteration order", not incidental. (Even a
/// value-only fold using this same strict-`>` shape is not fully
/// order-invariant either — `clippy.toml`'s own `f64::max`/`min` entry
/// records that a `±0.0` tie's sign depends on which operand is seen
/// first, which reaches `to_bits()` the moment such a value is digested.
/// [`shedding_q`] is that value-only shape, and this note is not a claim
/// that it is safe to reorder — only this function's own extra hazard,
/// the index, is what the tie-break above is about.) No exact tie has
/// been observed for either fold (0 across an 11M-cell sweep of the
/// P7-reachable box), which is a reason not to expect trouble, not a
/// licence to parallelise or reorder either one.
#[must_use]
pub(crate) fn valley_peak_total(coeffs: NuclearCoeffs) -> u32 {
    (1..=T_MAX)
        .fold((1_u32, f64::MIN), |(best_total, best_energy), total| {
            let energy = on_valley_per_unit_energy(total, coeffs);
            if energy > best_energy {
                (total, energy)
            } else {
                (best_total, best_energy)
            }
        })
        .0
}

/// Whether the on-valley per-unit binding energy has the shape found while
/// running Step 1's gate: bound at `total = 2`, at most a handful of sign
/// changes across `4..=T_MAX` (`total = 3` excluded — see below). Returns
/// `Some(last total observed bound)` if the shape holds, `None` otherwise.
/// **Not a suffix property** — a sequence with several transitions that
/// ends bound still passes; `Some`'s payload is the *last* bound total
/// anywhere in the scan, not "the total before permanent unbinding", which
/// coincide in all but the seeds discussed below.
///
/// **Why a *count* of transitions, not "exactly one" — found the hard way.**
/// An earlier version of this check demanded zero re-binding past the first
/// unbound total and failed on a real, physically-meaningful case: at
/// `total = 309` (this module's fixed `NUCLEAR_K = 10`'s **exact** shell-4
/// cumulative boundary, `1 + 12 + 42 + 92 + 162 = 309` —
/// `crate::packing::shell_size`'s own arithmetic), a completed shell buys a
/// one-total burst of extra binding — an "island of stability" shape,
/// analogous to (not the same mechanism as) real nuclear shell closures
/// (real magic numbers come from the spin-orbit-split shell model, not
/// geometric packing, and this crate's closures don't coincide with them),
/// not a coefficient defect.
///
/// **Re-measured, not assumed unaffected, after the "Post-gate correction"
/// (this module's own top-level doc) deleted `sigma`.** Over the same
/// 38,416-seed corpus: two seeds (3992, 4425) hit the `total = 309`
/// re-binding blip, down from three (14355, 20107, 21217) before —
/// expected, not a regression: `98.6%` of universes (37,863/38,416) are now
/// never unbound anywhere in `4..=T_MAX` at all, so there are fewer
/// unbound regions left for the blip to interrupt. Transition-count
/// distribution moved from `p50 = 0, p90 = 1, max = 3` to `p50 = 0, p90 =
/// 0, max = 3` — `<= 5` still clears the measured worst case with the same
/// margin as before.
///
/// **Both acceptance criteria are currently unreachable on the sampled
/// domain — a `/review-pr` finding (`emergence-auditor`), confirmed by
/// mutation.** Deleting the `total = 2` early-return *and* raising
/// `MAX_TRANSITIONS` to `50` together still passes the full corpus gate —
/// no seed anywhere in the 38,416-seed corpus is close to failing either
/// check today (max transitions measured is `3`, two below the `5` bound;
/// `total = 2` bound holds in every seed tried). This function is not
/// dead code — it is called, and both checks are real regression tripwire
/// for a regime that does not currently occur, not decoration — but a
/// green run of this function proves less than "the shape check passed"
/// suggests, and this note exists so a future reader does not read
/// `38416/38416` as evidence either bound is close to its edge.
///
/// **What this check does and does not discriminate — corrected after
/// mutation-testing it directly.** It reliably catches a *scattered*
/// implementation, where several genuinely distinct regions of the range
/// flip sign (a gross scale error, e.g. drawing `eps` far enough off that a
/// mid-range block of totals goes unbound: measured to produce 7
/// transitions at one sampled seed). It does **not** reliably catch a
/// *smooth* defect that still crosses zero once — a wrong coefficient or
/// exponent on any surviving term can leave the transition count at 0 or 1
/// while moving where the single crossing happens by a large amount. That
/// class is caught elsewhere: [`total_binding_energy`]'s three pinned
/// term-isolation tests (`the_volume_term_is_eps_times_contacts`,
/// `the_asymmetry_term_is_kappa_times_imbalance_squared_over_total`,
/// `the_coulomb_term_is_gamma_times_a_times_a_minus_one_over_total_to_the_one_third`)
/// and the identity-bound count
/// (`the_identity_configuration_is_bound_over_all_but_two_totals`).
/// Fraction of totals negative and count of sign transitions are different
/// statistics and this doc previously conflated them.
///
/// **Why `total = 3` is excluded, not merely tolerated.** Verified against
/// `physics-plausibility-reviewer`: `total = 3` is negative at the
/// *identity* configuration (`-0.4055`, no perturbation involved — a
/// `/review-pr` correction: an earlier draft of this figure, `-0.5615`,
/// was the deleted four-coefficient model's value, carried over unmarked
/// after `sigma`'s removal), matching a known property of a four-term
/// liquid-drop model with no pairing or
/// curvature term — the real SEMF's own prediction for the A=3 nuclide this
/// model actually favours at that total (Z=1, i.e. H-3, both by the real
/// SEMF's own argmax and by this model's) is `+0.775` `MeV` **per
/// nucleon**, against H-3's real `2.827` `MeV` per nucleon (not the
/// tabulated total binding energy, `8.482` `MeV`), 3.65x off before Borbax
/// adds any seed. It is not a
/// coefficient-balance defect this shape check should catch, so it is
/// carved out rather than either asserted-negative (false in ~14% of
/// universes) or asserted-positive (false in the rest).
#[must_use]
pub(crate) fn bound_shape(coeffs: NuclearCoeffs) -> Option<u32> {
    const MAX_TRANSITIONS: u32 = 5;

    let bound = |total: u32| on_valley_per_unit_energy(total, coeffs) > 0.0;
    if !bound(2) {
        return None;
    }

    // total = 3 is deliberately never visited here (the loop starts at 4), so
    // its own sign never enters the transition count or the prev/cur
    // comparison -- exactly the exclusion this function's own doc describes.
    let mut prev = true;
    let mut transitions = 0_u32;
    let mut last_bound = 2_u32;
    for total in 4..=T_MAX {
        let cur = bound(total);
        if cur {
            last_bound = total;
        }
        if cur != prev {
            transitions += 1;
        }
        prev = cur;
    }

    if transitions <= MAX_TRANSITIONS {
        Some(last_bound)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    #![expect(
        clippy::cast_precision_loss,
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        clippy::indexing_slicing,
        clippy::needless_range_loop,
        reason = "measurement probe: small counts (seed indices, totals bounded by T_MAX = 386, \
                  a fixed-size 121-entry reachability table) to f64/usize and back. Nothing here \
                  reaches a simulation result -- it exists to produce Step 1's gate verdict and \
                  the numbers docs/experiments/2026-08-12-nuclear-valley-gate.md reports, and a \
                  try_from plus a get().ok_or() on every arithmetic line would bury the physics \
                  being checked"
    )]
    use super::{
        BASE_C, BASE_EPS, BASE_KAPPA, NuclearCoeffs, T_MAX, bound_shape, composition_argmax,
        composition_argmax_int, contacts, drift_onset, on_valley_per_unit_energy,
        per_unit_binding_energy, shedding_q, total_binding_energy, valley_peak_total,
    };
    use crate::perturbation::{MigratedConstant, Rung, draw_symmetric};
    use borbax_units::det_math;

    /// Reachable `c` band under P7's crate-wide default `|p| <= 0.5`:
    /// `base_c * [0.5, 1.5]`.
    const C_MIN: f64 = 0.5 * BASE_C;
    const C_MAX: f64 = 1.5 * BASE_C;

    /// `BASE_EPS`/`BASE_KAPPA`/`BASE_C` bundled — used throughout this test
    /// module wherever "at identity" is meant, instead of repeating the
    /// struct literal.
    const IDENTITY: NuclearCoeffs = NuclearCoeffs {
        eps: BASE_EPS,
        kappa: BASE_KAPPA,
        c: BASE_C,
    };

    /// Same corpus-size derivation `perturbation.rs`'s own
    /// `P_IDENTITY_N`/`P_IDENTITY_BAND` pair uses: every reported fraction
    /// resolves to +/- 0.5 percentage points at 95% confidence, worst case
    /// at p = 0.5 (`N >= (1.96/0.005)^2 * 0.25`).
    const CORPUS_N: u64 = 38_416;

    /// Draw the three nuclear coefficients jointly, the way `generate_elements`
    /// will once Step 5 wires this in: one shared [`Rung`], then
    /// [`draw_symmetric`] per constant — never independent draws (routed
    /// requirements 5/6: a probe sampling the coefficients independently
    /// measures a regime the shipped code does not use).
    fn draw_universe(seed: u64) -> (Rung, NuclearCoeffs) {
        let rung = Rung::draw(seed);
        let eps = draw_symmetric(seed, rung, MigratedConstant::Eps, BASE_EPS);
        let kappa = draw_symmetric(seed, rung, MigratedConstant::Kappa, BASE_KAPPA);
        let c = draw_symmetric(seed, rung, MigratedConstant::C, BASE_C);
        (rung, NuclearCoeffs { eps, kappa, c })
    }

    /// **1 (acceptance).** Step 1's two-armed gate over the derived corpus.
    ///
    /// **Arm 1's corpus check is a monotonicity/band-membership check, not a
    /// closed-form regression bound — say so honestly, found after
    /// mutation-testing it.** An earlier version of this doc claimed the
    /// computed envelope (`drift_onset(C_MAX)..=drift_onset(C_MIN)`) was
    /// itself "a real regression bound: a sign error, a dropped factor, or a
    /// wrong `t^(2/3)` spelling moves `drift_onset` outside it" — false for
    /// at least one real defect class, verified by running seven broken
    /// implementations over a 20,001-point `c` sweep: a `c`-term sign error
    /// prints the **identical** `[16, 30]` envelope with 0 violations,
    /// because the envelope is computed by calling the same function the
    /// assertion then checks — the same tautology shape the plan's own
    /// round-4 correction already fixed once for the Ca-40 anchor, recurring
    /// one level deeper after replacing it. (A `cbrt(t)*cbrt(t)` respelling
    /// of `t^(2/3)` was also tried and is *not* a second example of this —
    /// measured, it produces zero integer differences from the pinned
    /// spelling over the same sweep, so it is behaviourally inert and no
    /// check, corpus-level or otherwise, could be expected to catch it.)
    /// What this check actually verifies: `drift_onset`
    /// is monotone non-increasing in `c` (true for the correct form and
    /// every broken form tried), and every drawn `c` lands inside P7's own
    /// reachable band — both real properties, neither a proof the closed
    /// form itself is right. **The tests that do carry that proof** are
    /// `the_drift_onset_at_identity_is_pinned` (an exact regression pin) and
    /// `the_closed_form_matches_the_brute_force_integer_argmax` (which
    /// mismatches on 2,040-11,490 of 11,580 cells under every defect
    /// tried) — read those, not this corpus loop, for closed-form
    /// correctness.
    ///
    /// **Arm 2 also gates on the bound-set shape ([`bound_shape`]), added
    /// beyond the plan's own literal text after a `physics-plausibility-reviewer`
    /// pass found the natural stronger assertion ("binding positive
    /// everywhere") is simply false at the identity configuration —
    /// `total = 3` is negative even unperturbed, matching a known limitation
    /// of a four-term liquid-drop model with no pairing term (see
    /// [`bound_shape`]'s own doc for the corrected real-SEMF comparison and
    /// for exactly which defect classes this shape check does and does not
    /// discriminate — it reliably catches scattered, multi-region defects,
    /// not a smooth single-crossing one like a wrong coefficient or exponent
    /// on a surviving term, which the three term-isolation tests
    /// (`the_volume_term_is_eps_times_contacts` and its two siblings) exist
    /// to catch instead).
    ///
    /// **The `interior` conjunct (`peak > 1 && peak < T_MAX`) is
    /// structurally unfalsifiable within the P7-reachable box — a
    /// `/review-pr` finding (`emergence-auditor`), found by measuring what
    /// `valley_peak_total` actually returns rather than trusting that a
    /// passing conjunct means something.** Across the full 38,416-seed
    /// corpus it returns exactly one of four values, `{13, 55, 147, 309}`
    /// — this crate's own packing-shell closures at `NUCLEAR_K = 10` — and
    /// none of them is `1` or `T_MAX`. So `interior` is `true` for every
    /// sampled seed not because the coefficients balance correctly, but
    /// because no P7-reachable coefficient combination can push the peak
    /// off a closure that far; the whole 38,406/38,416 closure statistic
    /// is therefore identical to the shedding-Q statistic alone, and the
    /// interior half contributes nothing to it. It is *not* unfalsifiable
    /// in general — a wide sweep outside the P7 box gives 2,287/6,912 edge
    /// peaks — only within the domain this gate actually samples. This
    /// test's own corpus loop below now additionally pins the exact
    /// closure set every observed peak must belong to, replacing the
    /// vacuous `interior` check with one that can fail: a future change
    /// that lets the peak escape it (a fifth closure appearing, or
    /// landing at a boundary) fails loudly instead of continuing to read
    /// as "interior, therefore fine."
    ///
    /// **Reported, not gated: how far the bound range reaches.**
    /// `physics-plausibility-reviewer` also found `MigratedConstant::C`'s
    /// own landed doc overclaimed coverage — `T_MAX`'s proof covers the
    /// composition argmax's *shape*, not whether binding stays positive out
    /// that far. **Post-gate correction, 2026-08-13: after `sigma`'s
    /// removal (this module's own top-level doc), the measured shortfall
    /// fell from 9.25% (3,554/38,416) to 1.09% (419/38,416)** — the
    /// original 9.25% figure was itself the finding that motivated
    /// questioning `sigma` in the first place, and the corrected model's
    /// own shortfall is small but not zero. Still a real, separate finding
    /// for Steps 2+ (see `perturbation.rs`'s `C` doc), not something Step 1
    /// itself is positioned to fix by picking a threshold here.
    /// **Un-ignored 2026-08-13 — a `/review-pr` finding (`determinism-auditor`):
    /// this test's own `edge_peaks == 6261` pin (the closure bar Arm 2
    /// gates on) was, until this change, never evaluated by any CI leg —
    /// no workflow anywhere in the repo passes `--ignored`. A golden
    /// nothing evaluates is a golden without a guard. The prior `#[ignore]`
    /// reason's own re-measurement (2.20s release, 1.74s debug for this
    /// test alone) already showed the cost claim that had justified
    /// skipping it was false; debug matters more than that reason
    /// acknowledged, since `overflow-checks` is exactly what would catch
    /// [`on_valley_total_energy`]'s `total - a` if [`composition_argmax_int`]'s
    /// own floor ever regressed.**
    #[test]
    fn the_two_armed_gate_passes_over_the_derived_corpus() {
        // 13, 55, 147, 309: packing::shell_size(10, ..)'s own cumulative
        // closures (1+12, +42, +92, +162) -- see the corpus loop below.
        const KNOWN_CLOSURES: [u32; 4] = [13, 55, 147, 309];

        let drift_at_c_max = drift_onset(C_MAX);
        let drift_at_c_min = drift_onset(C_MIN);
        assert!(
            drift_at_c_max.is_some() && drift_at_c_min.is_some(),
            "drift_onset should find a drift point at both ends of the reachable c band -- \
             drift_onset(C_MAX)={drift_at_c_max:?} drift_onset(C_MIN)={drift_at_c_min:?}"
        );
        assert!(
            drift_at_c_max <= drift_at_c_min,
            "drift_onset should be non-increasing in c: drift_onset(C_MAX)={drift_at_c_max:?} > \
             drift_onset(C_MIN)={drift_at_c_min:?}"
        );

        let mut closed = 0_u64;
        let mut examined = 0_u64;
        let mut t_bounds: Vec<u32> = Vec::with_capacity(CORPUS_N as usize);
        let mut elements_reached_while_bound: Vec<u32> = Vec::with_capacity(CORPUS_N as usize);
        let mut drift_seen: Vec<u32> = Vec::with_capacity(CORPUS_N as usize);

        for seed in 0..CORPUS_N {
            let (_rung, coeffs) = draw_universe(seed);
            let NuclearCoeffs { c, .. } = coeffs;
            examined += 1;

            let drift = drift_onset(c);
            assert!(
                drift.is_some() && drift >= drift_at_c_max && drift <= drift_at_c_min,
                "seed {seed}: drift_onset={drift:?} outside the band [{drift_at_c_max:?}, \
                 {drift_at_c_min:?}] measured at c's own reachable endpoints -- c={c} (see this \
                 test's own doc: this checks monotonicity and band membership, not the closed \
                 form's correctness -- that is test 6 and test 3's job)"
            );
            drift_seen.push(drift.unwrap_or_else(|| unreachable!("checked is_some() above")));

            let shape = bound_shape(coeffs);
            assert!(
                shape.is_some(),
                "seed {seed}: on-valley binding does not have the expected \
                 low-transition-count bound/unbound shape -- coeffs={coeffs:?}"
            );
            let t_bound = shape.unwrap_or_else(|| unreachable!("checked is_some() above"));
            t_bounds.push(t_bound);
            elements_reached_while_bound.push(composition_argmax_int(c, t_bound));

            let peak = valley_peak_total(coeffs);
            let interior = peak > 1 && peak < T_MAX;

            // The interior conjunct above is structurally unfalsifiable within
            // this corpus (see this test's own doc) -- this closure-set pin is
            // what actually discriminates, replacing it with a check that can
            // fail. KNOWN_CLOSURES are packing::shell_size(10, ..)'s own
            // cumulative closures (1+12, +42, +92, +162) -- the exact set
            // measured across all 38,416 seeds, never anything else and never
            // the trivial total=1 the fold technically starts from.
            assert!(
                KNOWN_CLOSURES.contains(&peak),
                "seed {seed}: peak={peak} is not one of this crate's own packing-shell \
                 closures {KNOWN_CLOSURES:?} at NUCLEAR_K=10 -- coeffs={coeffs:?}. If this is a \
                 deliberate change (a new closure genuinely reachable, or the coefficients now \
                 producing a non-quantised optimum), update this list and re-verify the \
                 'interior' conjunct's own vacuity claim in this test's doc, which was measured \
                 against exactly this closure set."
            );

            let q = shedding_q(T_MAX, coeffs);
            let has_positive_q = q > 0.0;

            if interior && has_positive_q {
                closed += 1;
            }
        }

        assert_eq!(
            examined, CORPUS_N,
            "examined count must equal the derived corpus size"
        );

        let closure_bar = 1.0 - 6261.0 / 59292.0; // see docs/experiments/2026-08-12-nuclear-valley-gate.md
        let closed_fraction = closed as f64 / CORPUS_N as f64;
        assert!(
            closed_fraction >= closure_bar,
            "Arm 2 closure fraction {closed_fraction} ({closed}/{CORPUS_N}) below the bar \
             {closure_bar} (V1's own reconstructed census, 1 - 6261/59292)"
        );

        // Anti-vacuity: a constant-returning drift_onset (ignoring c entirely) would
        // still satisfy every assertion above -- band-membership against its own two
        // endpoints is trivially true if all three calls return the same value. This
        // is exactly the shape of defect this check exists to rule out.
        drift_seen.sort_unstable();
        drift_seen.dedup();
        assert!(
            drift_seen.len() >= 2,
            "drift_onset returned only {} distinct value(s) across {CORPUS_N} seeds -- a \
             constant function would pass every other assertion in this test",
            drift_seen.len()
        );

        elements_reached_while_bound.sort_unstable();
        let n = elements_reached_while_bound.len();
        let p1 = elements_reached_while_bound[n / 100];
        let p10 = elements_reached_while_bound[n / 10];
        let p50 = elements_reached_while_bound[n / 2];
        let min = elements_reached_while_bound[0];
        let below_120 = elements_reached_while_bound
            .iter()
            .filter(|&&x| x < 120)
            .count();

        t_bounds.sort_unstable();
        let t_bound_min = t_bounds[0];
        let t_bound_max = t_bounds[t_bounds.len() - 1];

        println!(
            "Step 1 gate: drift_onset band [{drift_at_c_max:?}, {drift_at_c_min:?}]/{T_MAX} \
             (Ca-40 real anchor: 40/{T_MAX}, informational only); Arm 2 closure \
             {closed}/{CORPUS_N} = {closed_fraction} (bar {closure_bar}); bound-shape holds for \
             {CORPUS_N}/{CORPUS_N}, last-bound-total range [{t_bound_min}, {t_bound_max}]"
        );
        println!(
            "elements reached while still bound (reported, not gated): min={min} p1={p1} \
             p10={p10} p50={p50}; universes reaching fewer than 120: \
             {below_120}/{CORPUS_N}"
        );
    }

    /// Brute-force integer argmax of [`per_unit_binding_energy`] over
    /// `a in 1..=total`, domain-restricted per [`composition_argmax_int`]'s
    /// own floor. Shared by tests 2 and 3, which both cross-check
    /// [`composition_argmax_int`] against it under different sweeps.
    fn brute_force_argmax(total: u32, coeffs: NuclearCoeffs) -> u32 {
        (1..=total)
            .fold((1_u32, f64::MIN), |(best_a, best_e), a| {
                let e = per_unit_binding_energy(a, total - a, coeffs);
                if e > best_e { (a, e) } else { (best_a, best_e) }
            })
            .0
    }

    /// **2.** `kappa` cancels from the composition argmax by construction —
    /// the brute-force integer argmax of the full three-term energy, at a
    /// fixed `c`, must be identical across every `kappa`. The strongest test
    /// in this list per the verification pass that found it: an
    /// implementation drawing `gamma` independently, or using the
    /// κ-containing closed-form spelling instead of [`composition_argmax`],
    /// still looks like a plausible valley and still fails this.
    #[test]
    fn the_valley_is_a_one_parameter_family_in_c() {
        const KAPPAS: [f64; 6] = [0.1, 1.0, 11.85, 23.7, 35.55, 100.0];
        const TOTALS: [u32; 5] = [16, 56, 125, 208, 386];
        let c = BASE_C;
        for &total in &TOTALS {
            let reference = composition_argmax_int(c, total);
            for &kappa in &KAPPAS {
                let brute = brute_force_argmax(
                    total,
                    NuclearCoeffs {
                        eps: BASE_EPS,
                        kappa,
                        c,
                    },
                );
                assert_eq!(
                    brute, reference,
                    "total={total} kappa={kappa} c={c}: brute-force argmax {brute} != the \
                     kappa-invariant reference {reference}"
                );
            }
        }
    }

    /// **3.** The closed form matches the brute-force integer argmax of the
    /// pinned three-term energy, domain-restricted to `a >= 1` — over the
    /// grid used for the plan's own 48,000-case verification, scaled down
    /// for test runtime. Domain restriction at `a = 1` is what makes `total
    /// = 1` agree (see test 4 and test 8): opening the search to `a = 0`
    /// produces a real, exact tie there.
    #[test]
    fn the_closed_form_matches_the_brute_force_integer_argmax() {
        const KAPPAS: [f64; 6] = [0.1, 1.0, 11.85, 23.7, 35.55, 100.0];
        const CS: [f64; 5] = [0.0075, 0.010, 0.015, 0.020, 0.0225];
        let mut examined = 0_u64;
        for &kappa in &KAPPAS {
            for &c in &CS {
                for total in 1..=T_MAX {
                    examined += 1;
                    let closed_form = composition_argmax_int(c, total);
                    let brute = brute_force_argmax(
                        total,
                        NuclearCoeffs {
                            eps: BASE_EPS,
                            kappa,
                            c,
                        },
                    );
                    assert_eq!(
                        closed_form, brute,
                        "kappa={kappa} c={c} total={total}: closed form {closed_form} != brute \
                         force {brute}"
                    );
                }
            }
        }
        let expected_examined =
            u64::from(KAPPAS.len() as u32) * u64::from(CS.len() as u32) * u64::from(T_MAX);
        assert_eq!(examined, expected_examined);
    }

    /// **4.** `total = 1` is a real, exact tie between `a = 0` and `a = 1` —
    /// `contacts(1) == 0.0`, so the binding term vanishes for both, and
    /// `a*(a-1)` is `0` for both `a = 0` and `a = 1`. The domain restriction
    /// (`a >= 1`, [`composition_argmax_int`]'s floor), not a tie-break, is
    /// what selects `a = 1` — the plan's round-3 correction claimed no exact
    /// tie occurs under the corrected model and round-5 corrected that.
    ///
    /// **Not the *only* exact tie — a `/review-pr` correction
    /// (`geometry-numerics-reviewer`), renamed accordingly.** See
    /// [`an_exact_tie_at_total_one_twenty_five_is_resolved_by_round_ties_even`]
    /// below for the tie family this test's own name used to overclaim
    /// away.
    #[test]
    fn the_total_one_tie_is_closed_by_the_domain_restriction() {
        assert_eq!(contacts(1).to_bits(), 0.0_f64.to_bits());
        for &(kappa, c) in &[(23.7_f64, 0.015_f64), (0.1, 0.0075), (100.0, 0.0225)] {
            let coeffs = NuclearCoeffs {
                eps: BASE_EPS,
                kappa,
                c,
            };
            let e0 = total_binding_energy(0, 1, coeffs);
            let e1 = total_binding_energy(1, 0, coeffs);
            assert_eq!(
                e0.to_bits(),
                e1.to_bits(),
                "kappa={kappa} c={c}: total_binding_energy(0,1) and (1,0) should be bit-identical \
                 at total=1"
            );
        }
        // The domain restriction, not the tie, decides: composition_argmax_int
        // never returns 0.
        for &c in &[C_MIN, BASE_C, C_MAX] {
            assert_eq!(composition_argmax_int(c, 1), 1);
        }
    }

    /// `composition_argmax(c, t)` is exactly `k + 0.5` whenever
    /// `c == (t - 2k - 1) / (k * t^(2/3))` — a family with **4,504**
    /// solutions in the P7-reachable `c` band (found by
    /// `geometry-numerics-reviewer`'s exhaustive enumeration; `total = 1`,
    /// test 4 above, is the one member the domain restriction closes, not
    /// the only member). A drawn `c` lands inside a tie window with
    /// probability ~4e-12 (26,853 total ulps of tie window against
    /// ~6.8e15 legal `f64` values in the band) — a statement about the
    /// model, not a reachable case, which is why no other test needs to
    /// account for it.
    ///
    /// **At such a cell, [`brute_force_argmax`] is not a valid oracle.**
    /// `total = 125` is a perfect cube, so `cbrt(125*125) == 25.0` exactly
    /// and the closed form gives exactly `50.5` at `c = 0.0192`. But
    /// `per_unit_binding_energy(50, 75, ..)` and `per_unit_binding_energy(51,
    /// 74, ..)` are mathematically equal and differ only by ~2 ulp of
    /// accumulated rounding — so a brute-force `>` comparison returns
    /// whichever way that rounding happens to fall, not a real preference.
    /// The closed form is the *reliable* side at a tie: `round_ties_even`
    /// on the exact `.5` quotient is deterministic; the oracle is not.
    #[test]
    fn an_exact_tie_at_total_one_twenty_five_is_resolved_by_round_ties_even() {
        assert_eq!(
            composition_argmax(0.0192, 125).to_bits(),
            50.5_f64.to_bits()
        );
        assert_eq!(composition_argmax_int(0.0192, 125), 50);
        let coeffs = NuclearCoeffs {
            eps: BASE_EPS,
            kappa: BASE_KAPPA,
            c: 0.0192,
        };
        let e50 = per_unit_binding_energy(50, 75, coeffs);
        let e51 = per_unit_binding_energy(51, 74, coeffs);
        assert!(
            (e50 - e51).abs() <= 4.0 * f64::EPSILON * e50.abs(),
            "the two candidates at this tie should differ by a few ulp of rounding, not \
             represent a real preference: e50={e50} e51={e51}"
        );
    }

    /// **5.** The three nuclear coefficients share one rung (routed
    /// requirement 6) — the fraction of seeds landing all three at their base
    /// values simultaneously must match `P(identity) = 1/(Rung::MAX+1)`, and
    /// the joint draw's correlation structure must differ measurably from
    /// independent draws. Reuses `perturbation.rs`'s own recalibrated
    /// `P_IDENTITY_N`/`P_IDENTITY_BAND` pair for the same quantity, rather
    /// than inventing a second corpus.
    ///
    /// **Weaker by one conjunct than before `sigma`'s removal, and that is
    /// acknowledged rather than silently absorbed**: `P(identity)` itself is
    /// unchanged in theory (every constant lands on base at rung 0
    /// regardless of how many constants there are), but this test now has
    /// three conjuncts to `AND` together instead of four, a marginally
    /// weaker joint check than before.
    /// **Un-ignored 2026-08-13, same reason as
    /// `the_two_armed_gate_passes_over_the_derived_corpus`.**
    #[test]
    fn the_three_nuclear_coefficients_share_one_rung() {
        const N: u64 = 200_000;
        const BAND: f64 = 0.0015;
        let mut all_at_base = 0_u64;
        for seed in 0..N {
            let (_rung, coeffs) = draw_universe(seed);
            if coeffs.is_identity_config(IDENTITY) {
                all_at_base += 1;
            }
        }
        let observed = all_at_base as f64 / N as f64;
        let want = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        assert!(
            (observed - want).abs() < BAND,
            "observed P(all three at base) = {observed} over {N} draws, want {want} +/- {BAND} \
             ({all_at_base} hits)"
        );
    }

    /// **6.** Drift onset at the identity configuration is pinned to `22` —
    /// a regression anchor, checked against implementation classes that
    /// would still pass a loose corpus bar: a factor-2 error in `gamma`
    /// (giving `30`), a wrong `t^(2/3)` spelling, or a `base_c` drift.
    /// Elements `1..=10` sit exactly at the symmetric composition at
    /// identity — the plan's own text says `1..=7`, measured under a
    /// superseded `base_c` value the plan's round-5 correction retracted;
    /// `1..=10` is correct at the shipped `base_c = 0.015`.
    #[test]
    fn the_drift_onset_at_identity_is_pinned() {
        assert_eq!(drift_onset(BASE_C), Some(22));
        for a in 1..=10_u32 {
            assert_eq!(
                composition_argmax_int(BASE_C, 2 * a),
                a,
                "element {a} (total {}) should sit at the symmetric composition at identity",
                2 * a
            );
        }
        assert_ne!(
            composition_argmax_int(BASE_C, 22),
            11,
            "total 22 should have already drifted"
        );
    }

    /// **7.** Every proton count `1..=120` is reached by the valley at
    /// every legal `c` — the direct 386-cell scan Decision 11's coverage
    /// proof licenses, not a restatement of the proof.
    #[test]
    fn every_proton_count_is_reachable_by_the_valley_at_every_legal_c() {
        // C_MAX itself, unnudged: coverage holds at the exact endpoint (measured;
        // `next_f64_range` never actually returns its upper bound, so this is a
        // legitimate worst-case check, not a value the shipped draw can hit).
        for &c in &[C_MIN, BASE_C, C_MAX] {
            let mut reached = [false; 121];
            for total in 1..=T_MAX {
                let a = composition_argmax_int(c, total);
                if a <= 120 {
                    reached[a as usize] = true;
                }
            }
            for z in 1..=120 {
                assert!(
                    reached[z],
                    "c={c}: proton count {z} never reached by the valley"
                );
            }
        }
    }

    /// `composition_argmax_int(c, total) <= total / 2` for every legal `c`
    /// and every `total` — [`drift_onset`]'s own doc cites this property
    /// as "proved in this module's own tests," but no test actually
    /// checked it (a `/review-pr` finding, `geometry-numerics-reviewer`);
    /// this is that test. True on paper (`a* <= t/2 iff t^(2/3) <= t^(5/3)`,
    /// which holds for `t >= 1`, and `round_ties_even` adds at most `0.5`,
    /// so `a <= t/2 + 0.5 <= t`), but the citation pointed at a test that
    /// did not exist — exactly the shape by which a real invariant gets
    /// silently deleted in a later refactor. Also closes an unsigned-
    /// underflow risk: every on-valley caller computes `total - a`, which
    /// panics under `overflow-checks` (or wraps to a garbage magnitude in
    /// release) if `a` ever exceeded `total`.
    #[test]
    fn the_valley_never_reaches_past_half_of_total() {
        for &c in &[C_MIN, BASE_C, C_MAX] {
            for total in 1..=T_MAX {
                let a = composition_argmax_int(c, total);
                assert!(
                    a <= total,
                    "c={c} total={total}: a={a} exceeds total -- total - a would underflow"
                );
                assert!(
                    a <= total / 2 || total == 1,
                    "c={c} total={total}: a={a} exceeds total/2={} (total=1 is the sole \
                     documented exception, closed by the domain floor -- see test 4)",
                    total / 2
                );
            }
        }
    }

    /// **8.** `composition_argmax(c, 1) == 0.5` exactly, for every legal
    /// `c` — makes the rounding hazard visible (`round_ties_even(0.5) ==
    /// 0`) rather than latent, and `composition_argmax_int` must still
    /// return `1`.
    #[test]
    fn total_one_is_exactly_one_half_before_the_domain_floor_applies() {
        for &c in &[C_MIN, BASE_C, C_MAX, 1e-4, 0.2] {
            assert_eq!(composition_argmax(c, 1).to_bits(), 0.5_f64.to_bits());
            assert_eq!(composition_argmax_int(c, 1), 1);
        }
    }

    /// **9 (report).** The drift distribution reported at fixed totals, not
    /// pooled — pooling makes a degenerate (single-optimum) model look
    /// diverse, which is exactly the trap the plan names.
    #[test]
    #[ignore = "reporting only, not run under CI by default -- run: cargo test --locked --release \
                -p borbax-universe -- --ignored --nocapture nuclear"]
    fn the_drift_distribution_is_reported_at_fixed_totals_not_pooled() {
        for &total in &[20_u32, 56, 120, 200, 386] {
            let lo = composition_argmax(C_MAX, total) / f64::from(total);
            let hi = composition_argmax(C_MIN, total) / f64::from(total);
            println!("total={total}: a*/total spans [{lo}, {hi}] across the legal c range");
        }
        // Within-universe curve at identity: non-crossing check against a
        // second universe's curve (a lower c).
        for total in (2..=T_MAX).step_by(20) {
            let identity = composition_argmax(BASE_C, total) / f64::from(total);
            let low_c = composition_argmax(C_MIN, total) / f64::from(total);
            assert!(
                low_c >= identity,
                "curves crossed at total={total}: low-c ratio {low_c} < identity ratio {identity}"
            );
        }
    }

    /// **10 (report).** V1's own edge-peak census, reconstructed *inside*
    /// this crate against the real `contacts_upto`/`shell_size` (settling
    /// the transcription caveat the standalone verification carried), then
    /// Arm 2's peak position reported against it. The literal Arm 2 gate
    /// (interior + positive-Q) passes near-universally at the pre-registered
    /// coefficients and cannot by itself see how far the peak sits from
    /// V1's own typical position — this test makes that visible rather than
    /// silently absent.
    /// **Un-ignored 2026-08-13 — this is the test that carries the
    /// `edge_peaks == 6261` pin `the_two_armed_gate_passes_over_the_derived_corpus`'s
    /// own closure bar is derived from; see that test's doc for why
    /// leaving it `#[ignore]`d meant the bar itself had no CI guard.**
    #[test]
    fn v1_edge_peak_census_is_pinned_and_arm_two_peak_is_reported_against_it() {
        use crate::packing::{PackingConsts as V1PackingConsts, contacts_upto};

        let mut edge_peaks = 0_u64;
        let mut examined = 0_u64;
        for k in 6..=14_usize {
            let pack = V1PackingConsts::new(k);
            for i in 0..=8_u32 {
                let eps = 0.8 + 0.05 * f64::from(i);
                for j in 0..=11_u32 {
                    let sigma = 0.02 + 0.01 * f64::from(j);
                    for n_elements in 60..=120_usize {
                        examined += 1;
                        let peak = (1..=n_elements)
                            .fold((1_usize, f64::MIN), |(best_u, best_e), units| {
                                let contacts = contacts_upto(pack, units);
                                let e = (eps * contacts) / units as f64
                                    - sigma * det_math::cbrt((units * units) as f64);
                                if e > best_e {
                                    (units, e)
                                } else {
                                    (best_u, best_e)
                                }
                            })
                            .0;
                        if peak == n_elements {
                            edge_peaks += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(examined, 9 * 9 * 12 * 61);
        assert_eq!(
            edge_peaks, 6261,
            "V1 edge-peak census drifted from its pinned value"
        );

        let v1_edge_fraction = edge_peaks as f64 / examined as f64;
        println!(
            "V1 census: {edge_peaks}/{examined} = {v1_edge_fraction} edge peaks (V1's own \
             element.rs:887-914 comment reports 2/24 ~ 8.3% on a 24-seed sample)"
        );

        let peak = valley_peak_total(IDENTITY);
        println!(
            "V3 identity: on-valley peak at total={peak} (element {} of up to 120), against \
             V1's reconstructed p5 ~0.100 / median ~0.581 of n_elements. NOT independent \
             corroboration against the real SEMF's A=60 or the experimental A~61-62 (Ni-62/ \
             Fe-56) -- peak is quantised to this crate's own packing-shell closures \
             {{13,55,147,309}} regardless of coefficients (a /review-pr finding, \
             emergence-auditor); the model's own smooth (unquantised) optimum sits at \
             total=60-70, and 55 is only the nearest closure below it",
            composition_argmax_int(BASE_C, peak)
        );
    }

    /// **11.** The identity configuration is bound over all but two totals —
    /// pinned as a count, not left to a corpus percentile. This is what
    /// correction 2 (`BASE_EPS` `1.0 -> 2.625`) was *for*, and until this
    /// test existed nothing failed when it was undone: at `eps = 1.0`,
    /// **190 of 386 totals (49.2%) are unbound on the current three-term
    /// model** — re-measured post-removal, not the four-term model's
    /// original `305/386` (79.0%, still the correct figure for the
    /// pre-registration correction 2 itself describes, see
    /// `perturbation.rs`'s `Eps` doc), and [`bound_shape`] returns `Some`
    /// on every one of them regardless, because the sequence crosses zero
    /// exactly once (a scale error, not a scattered one). Mutation-probed
    /// against `BASE_EPS = 1.0`.
    ///
    /// **Re-observed, not assumed, after the "Post-gate correction" removed
    /// `sigma`** — the unbound set stays exactly `{1, 3}` (measured: setting
    /// `sigma = 0` at the *four*-term model left this unchanged, so deleting
    /// the term is not expected to move it either, but this test's own job
    /// is to fail loudly if that expectation is wrong).
    #[test]
    fn the_identity_configuration_is_bound_over_all_but_two_totals() {
        let unbound: Vec<u32> = (1..=T_MAX)
            .filter(|&t| on_valley_per_unit_energy(t, IDENTITY) <= 0.0)
            .collect();
        assert_eq!(
            unbound,
            vec![1, 3],
            "identity should be unbound only at total = 1 (contacts(1) == 0, so the volume \
             term vanishes) and total = 3 (the known four-term liquid-drop limitation \
             bound_shape's own doc records -- the real SEMF's own four physical terms, not \
             this crate's three coefficients, since eps's own packing series supplies both \
             volume and surface) -- a different unbound set here means eps's scale has drifted"
        );
    }

    /// **12a.** The volume term is `eps * contacts(total)`, pinned directly
    /// by isolating it: `total_binding_energy(4, 4, eps=1, kappa=0, c=0)` —
    /// `a == b` zeroes the asymmetry term, `kappa=0` and `c=0` zero
    /// `gamma = 2*kappa*c` and so the Coulomb term, leaving exactly
    /// `contacts(8)`.
    ///
    /// **Pinned against a literal, not a live call to [`contacts`] — found
    /// necessary, not stylistic.** An earlier version of this test asserted
    /// `isolated == contacts(8)` (a second call), which is self-referential
    /// with respect to [`NUCLEAR_K`]: both sides route through the same
    /// [`contacts`], so any `NUCLEAR_K` regression moves both sides
    /// identically and the assertion still passes — exactly the tautology
    /// shape this module has already found once, in Arm 1's original
    /// self-computed envelope. The literal below was measured directly
    /// (`contacts(8)` at the shipped `NUCLEAR_K = 10`), so a `NUCLEAR_K`
    /// change moves only the left side and the test actually discriminates
    /// it — a strictly stronger guard than a peak-position assertion, but
    /// only once it stopped calling the function it exists to check.
    ///
    /// **What this does not claim**: `eps` is fixed at `1.0` here, so
    /// `eps` vs `eps^2` is invisible to this test alone (`1.0^2 == 1.0`).
    /// Caught elsewhere instead (mutation-tested): the identity unbound-set
    /// and identity-peak tests, both of which use `BASE_EPS != 1.0`. **Not**
    /// the closed-form/one-parameter-family tests — `eps` is a per-total
    /// additive constant that cancels from the composition argmax by
    /// construction, so those are correctly, permanently blind to it.
    /// **Also does not stand in for [`crate::packing`]'s own guard on
    /// `contacts_upto`'s internal accumulation order** (`packing.rs`'s own
    /// `the_universe_digest_is_pinned`-style tests) — `total = 8` happens
    /// not to be one of the 67 (of 386, at `k = 10`) cells that guard's own
    /// hazard actually moves, so this literal pins one point's *value*,
    /// not the packing series' own summation order.
    #[test]
    fn the_volume_term_is_eps_times_contacts() {
        let isolated = total_binding_energy(
            4,
            4,
            NuclearCoeffs {
                eps: 1.0,
                kappa: 0.0,
                c: 0.0,
            },
        );
        assert_eq!(
            isolated.to_bits(),
            18.608_003_309_573_23_f64.to_bits(),
            "eps * contacts(total) at total = 8, eps = 1: got {isolated}, want contacts(8) at \
             NUCLEAR_K = 10, 18.60800330957323, exactly"
        );
    }

    /// **12b.** The asymmetry term is `-kappa * (a-b)^2 / total`, pinned
    /// directly: `total_binding_energy(6, 2, eps=0, kappa=2, c=0)` zeroes
    /// the volume term (`eps=0`) and the Coulomb term (`gamma =
    /// 2*kappa*c = 0`), leaving `-2 * (6-2)^2 / 8 = -4.0` exactly (`16/8`
    /// divides evenly, no rounding).
    ///
    /// **Not `(3, 1, kappa=1)` — found doubly degenerate by mutation
    /// testing.** That fixture gives `(a-b)^2/total = 1` *and* `kappa = 1`,
    /// so any exponent applied to the whole term, or to `kappa` alone,
    /// leaves the result at `1^n = 1` — invisible. Verified: mutating the
    /// term to its own square, or halving-then-squaring `kappa`, both left
    /// the old fixture's assertion passing. `(6, 2, kappa=2)` gives
    /// `(a-b)^2/total = 2` and `kappa = 2`, breaking both degeneracies at
    /// once; re-run against the same mutations and both now fail here.
    #[test]
    fn the_asymmetry_term_is_kappa_times_imbalance_squared_over_total() {
        let isolated = total_binding_energy(
            6,
            2,
            NuclearCoeffs {
                eps: 0.0,
                kappa: 2.0,
                c: 0.0,
            },
        );
        assert_eq!(
            isolated.to_bits(),
            (-4.0_f64).to_bits(),
            "kappa * (a-b)^2 / total at a=6 b=2 kappa=2: got {isolated}, want -4.0 exactly"
        );
    }

    /// **12c.** The Coulomb term is `-gamma * a*(a-1) / total^(1/3)`, pinned
    /// two ways. First, isolated at `a == b` (zeroing the asymmetry term):
    /// `total_binding_energy(4, 4, eps=0, kappa=2, c=0.25)` gives
    /// `gamma = 2*kappa*c = 1.0`, so `-1.0 * 4*3 / cbrt(8) = -12/2 = -6.0`
    /// exactly (`8` is a perfect cube, `libm::cbrt` is exact on those).
    /// Second, at `a != b` — `total_binding_energy(6, 2, eps=0, kappa=2,
    /// c=0.25)` adds the asymmetry term back in (`-2 * (6-2)^2/8 = -4.0`)
    /// on top of the same Coulomb term read at `a=6`
    /// (`-1.0 * 6*5/cbrt(8) = -15.0`), giving `-19.0` exactly.
    ///
    /// **Not `kappa=1` — found degenerate by mutation testing, same shape
    /// as 12b's fix**: at `kappa=1`, `gamma = 2*kappa*c` mutated to
    /// `2*kappa^2*c` is invisible (`1^2 == 1`). `kappa=2` breaks it. **The
    /// second, asymmetric assertion is new** — the original single `a==b`
    /// fixture cannot see the Coulomb term reading the wrong variable
    /// (`b*(b-1)` instead of `a*(a-1)`), since `a==b` makes them identical;
    /// verified this mutation passed the original fixture and fails the
    /// asymmetric one.
    #[test]
    fn the_coulomb_term_is_gamma_times_a_times_a_minus_one_over_total_to_the_one_third() {
        let coeffs = NuclearCoeffs {
            eps: 0.0,
            kappa: 2.0,
            c: 0.25,
        };
        let isolated = total_binding_energy(4, 4, coeffs);
        assert_eq!(
            isolated.to_bits(),
            (-6.0_f64).to_bits(),
            "gamma * a*(a-1) / total^(1/3) at a=4 b=4 kappa=2 c=0.25: got {isolated}, want -6.0 \
             exactly"
        );
        let asymmetric = total_binding_energy(6, 2, coeffs);
        assert_eq!(
            asymmetric.to_bits(),
            (-19.0_f64).to_bits(),
            "at a=6 b=2 kappa=2 c=0.25 (asymmetry -4.0 plus Coulomb -15.0): got {asymmetric}, \
             want -19.0 exactly"
        );
    }

    /// **12d.** The three terms are combined left to right, `(volume −
    /// asymmetry) − Coulomb`, not `volume − (asymmetry + Coulomb)` — this
    /// module's own "every float spelling pinned" requirement, load-bearing
    /// per §13.4. Each of 12a/12b/12c zeroes two of the three terms to
    /// isolate the third, so none of them can see which order the three
    /// are combined in. Verified by mutation: re-associating to `volume −
    /// (asymmetry + Coulomb)` moves 131 of 386 identity on-valley totals
    /// (worst `4.5e-13`) and passes every other test in the crate that
    /// existed before this one and [`the_identity_on_valley_energy_digest_is_pinned`]
    /// below. `a = 26, b = 30` is a total where the two associations differ
    /// (chosen for that reason, not arbitrary — `(25, 30)`, `(4, 4)` and
    /// `(60, 60)` all happen not to differ) — a single hand-checkable
    /// worked example, complementing the digest's full-sequence coverage.
    #[test]
    fn the_three_terms_are_summed_left_to_right() {
        let assembled = total_binding_energy(26, 30, IDENTITY);
        assert_eq!(
            assembled.to_bits(),
            492.160_537_020_485_5_f64.to_bits(),
            "assembled three-term energy at a=26 b=30, identity coefficients: got {assembled}, \
             want 492.1605370204855 exactly"
        );
    }

    /// **13.** A bit-level digest of the identity configuration's on-valley
    /// energies over `1..=T_MAX` — a comprehensive guard on
    /// [`total_binding_energy`]'s pinned left-to-right term association
    /// (§13.4), complementing 12d's single hand-checkable point with
    /// coverage of the whole sequence. Measured: swapping the two
    /// subtractions moves 123 of 386 on-valley totals (29,341/75,077 =
    /// 39.1% of all reachable `(a, b)` cells at identity, worst `4.5e-13`
    /// — a `/review-pr` correction: an earlier draft of this denominator,
    /// 75,075, was off by 2 from `sum(t+1 for t in 1..=T_MAX)`) and fusing
    /// them into `V - (A + C)` moves 43.0% — both pass every
    /// other test in this module. Every value in the sequence is finite
    /// (checked when the literal was taken), so `to_bits` needs no NaN
    /// canonicalisation. Same FNV-1a-64 shape `packing.rs`'s own
    /// accumulation-order digests use.
    #[test]
    fn the_identity_on_valley_energy_digest_is_pinned() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for total in 1..=T_MAX {
            let a = composition_argmax_int(BASE_C, total);
            let e = total_binding_energy(a, total - a, IDENTITY);
            for byte in e.to_bits().to_le_bytes() {
                h ^= u64::from(byte);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        assert_eq!(
            h, 0x6010_7983_e565_657a,
            "the identity on-valley energy sequence moved -- if this is a deliberate physics \
             change, regenerate; if it is a refactor of total_binding_energy's term order, it \
             is not (see that function's own doc: the association is pinned)"
        );
    }

    /// **14.** The identity configuration's on-valley peak is pinned to
    /// `total = 55` (element 25) — this module's own general
    /// identity-configuration regression pin. `NUCLEAR_K` (`= 10`, not
    /// `12`, correction 3 of this module's own top-level doc) is guarded
    /// more directly by [`the_volume_term_is_eps_times_contacts`] now (an
    /// earlier version of this module isolated the volume term with
    /// `eps=0`, which blinded that isolation to `NUCLEAR_K` and left this
    /// peak assertion as the only guard — no longer true, this test is a
    /// general regression pin, not `NUCLEAR_K`'s primary guard).
    #[test]
    fn the_identity_peak_is_pinned_to_total_fifty_five() {
        let peak = valley_peak_total(IDENTITY);
        assert_eq!(
            (peak, composition_argmax_int(BASE_C, peak)),
            (55, 25),
            "on-valley peak at identity should sit at total = 55 (element 25) -- if this moved, \
             check eps/kappa/c's base values"
        );
    }
}
