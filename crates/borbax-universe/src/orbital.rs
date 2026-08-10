//! Screened-hydrogenic, `l`-dependent orbital filling (issue #26, Decisions
//! 3 and 8) — Task 26.1's replacement for `element.rs`'s V1 shell law.
//!
//! **Fill *order* and subshell *energy* are two different questions, answered
//! by two different mechanisms — this module's original version answered
//! both from one screened-hydrogenic energy formula, and a `/review-pr` pass
//! measured that this does not scale past `z ≈ 54`: a purely energy-driven
//! search puts 4f (and 5d) ahead of 6s, because a low-`n` candidate's
//! `1/n²` term wins the comparison however its screening is tuned — a
//! structural property of the formula, not a calibration miss (verified by
//! sweeping the screening coefficient from 0.65 to 0.99 with several
//! candidate corrections; none moved the crossover past `z ≈ 55`). This is
//! also the reason every real periodic-table teaching model uses an
//! explicit ordering rule rather than deriving order from screening alone.**
//!
//! # Fill order: Madelung's rule, not a search
//!
//! Subshells fill in [`madelung_order`]: ascending `(n + l)`, ties broken by
//! ascending `n` — the standard rule (Madelung 1936; independently
//! Klechkovsky) that correctly orders the *entire* periodic table's block
//! structure, with the same small set of individual-element exceptions
//! (Cr, Cu, and a handful of lanthanides/actinides) every simplified atomic
//! model shares with it. It needs no screening input at all: the order is
//! the same for every seed, drawn constants perturb the *energies* within
//! that fixed structure, not which subshell is next.
//!
//! # Subshell energy: Slater's rules, not a single screening constant
//!
//! For a subshell `(n, l)` given the electrons already placed,
//! `Z_eff(n, l) = z - S(n, l)` where `z` is the electron index about to be
//! placed and `S` is Slater's own two-tier screening rule (Slater 1930),
//! adapted to this module's continuous `t(l)` in place of Slater's discrete
//! `(ns,np)` grouping, and further split by the *screener's* own subshell
//! type via `sigma_deep` (this module's own extension — see below):
//!
//! ```text
//! near_sp(n) = electrons already placed with n' == n - 1, l' <= 1
//! near_df(n) = electrons already placed with n' == n - 1, l' >= 2
//! far_sp(n)  = electrons already placed with n' <= n - 2, l' <= 1
//! far_df(n)  = electrons already placed with n' <= n - 2, l' >= 2
//! same(n, l) = electrons already placed with n' == n, l' != l — the
//!              candidate's own already-placed peers are excluded
//! t(l) = min(l / 2.0, 1.0)             // 0 at l=0, 0.5 at l=1, saturated at l>=2
//!
//! near_rate(l) = if l <= 1 { sigma_near } else { 1.0 }   // by CANDIDATE's l
//! deep_rate    = min(1.0, sigma_deep)                     // by SCREENER's l
//!
//! S(n, l) = near_rate(l)          * near_sp(n)
//!         + min(near_rate(l), deep_rate) * near_df(n)
//!         +                          far_sp(n)
//!         + deep_rate              * far_df(n)
//!         + t(l) * same(n, l)
//!
//! energy(n, l) = -(Z_eff(n, l) / n)^2   // plain n, not a quantum defect
//! ```
//!
//! **The candidate-side split (`near_rate(l)`) — full screening from every
//! inner shell, not just the far ones, when the *candidate itself* is d/f —
//! is the piece the original single-constant model was missing, and it is
//! the reason 4f/5d no longer win against 6s.** Slater's own rule draws
//! exactly this distinction: an `ns`/`np` candidate is only partly screened
//! by its immediately-inner shell (electrons genuinely penetrate it to some
//! degree), but an `nd`/`nf` candidate — poorly penetrating by construction
//! — is screened at full strength by *every* group closer to the nucleus,
//! not only the ones two shells back. Verified against the real subshell
//! sequence: with `sigma_near` at Slater's own historical value (0.85), the
//! identity fold reproduces the real Madelung sequence exactly through
//! `z = 118` (this module's own tests pin it). **The candidate-side split
//! alone** (`element.rs`'s own radius formula, before `sigma_deep` below)
//! separately raised the model-vs-real covalent-radius correlation from
//! 0.50 to 0.87 and fixed the group-1 radius trend (Li through Fr) from
//! inverted to strictly increasing — neither of which this screening split
//! was tuned against; both were measured consequences of the same change,
//! at the time they were measured. **`sigma_deep` below moves both figures
//! again** (it changes `zeff_outer` for every `z >= 37`, which the group-1
//! trend reads directly) and neither has been re-measured since: `/review-pr`
//! round 3 (2026-08-10) found the group-1 trend is no longer strictly
//! increasing at HEAD (Cs-analogue to Fr-analogue now decreases), and found
//! no reproducible harness anywhere in this repository for the correlation
//! figure, so that one could not be re-checked at all rather than merely
//! being stale. Re-measuring both against the shipped `sigma_deep` physics,
//! with a real harness for the correlation figure, is follow-up work, not
//! restated here as fact.
//!
//! # The screener-side split: `sigma_deep`, and what it fixes
//!
//! **A `/review-pr` pass (2026-08-10, finding F3) measured that the
//! candidate-side split above leaves d/f-block `Z_eff` a pure integer
//! count, with zero dependence on any drawn constant, for every one of the
//! 68 d/f-block elements in a real generated table — and that this is also
//! the reason the model showed no lanthanide contraction at all: a flat,
//! bit-identical radius across the entire 14-element 4f-analogue block.**
//! This is not a bug in the candidate-side split — verified directly (see
//! this module's own doc's citation trail, and the two research passes
//! recorded in this task's own history): literal Slater 1930, applied
//! faithfully to a d/f *candidate*, gives exactly `near = far = 1.0`, no
//! partial tier of any kind. The model was already correct Slater; Slater
//! itself doesn't parameterise d/f screening the way it parameterises s/p
//! screening, and famously fails to predict any lanthanide contraction **for
//! the outermost-shell probe this module's radius reads** — the `zeff_outer`
//! quantity `element.rs`'s radius formula uses is exactly flat under literal
//! Slater, verified directly (a documented, textbook-known weakness of the
//! rule for that quantity, not an oversight in this implementation of it).
//! The claim is narrower than "Slater predicts no lanthanide contraction at
//! all": for the inner 5s/5p electrons that set the real Ln³⁺ *ionic*
//! radius `sigma_deep`'s own calibration target measures, 4f sits in
//! Slater's own `n-1` tier (partial screening, 0.85), so literal Slater
//! *does* contract that quantity — it is specifically the outer-probe
//! radius this module computes where the rule goes flat.
//!
//! **The fix reads screening imperfection off the *screener's* subshell,
//! not the candidate's — the standard atomic-physics teaching that
//! penetration ability runs symmetrically both ways.** A d/f electron is
//! not only *fully* screened *itself* (the candidate-side fact above — poor
//! penetration means every inner group screens it at full strength, not
//! partially); it is also, by that same poor-penetration property, a poor
//! screener *of everything else* — "the order of electron penetration from
//! greatest to least is s, p, d, f; the order of the amount of shielding
//! done is also in the order s, p, d, f" (standard teaching, e.g.
//! `LibreTexts`' "Penetration and Shielding"). `sigma_deep` is a ceiling
//! every d/f electron's own contribution to any *other* candidate's
//! `near_df`/`far_df` screening sum is capped at — the same-shell `t(l)`
//! term below is a separate mechanism `sigma_deep` does not reach, since a
//! same-shell peer is never routed through `near`/`far` at all. It is what
//! actually restores a
//! lanthanide-analogue contraction: during the 4f fill, the outer 6s probe
//! (see `zeff_outer`'s d/f-block fallback) reads `far(6)`, which counts the
//! newly-placed 4f electrons at `n' = 4 <= 6 - 2`; capping their
//! contribution below `1.0` stops each added proton from being exactly
//! cancelled by its accompanying 4f electron.
//!
//! **`sigma_deep` is a real, citable *qualitative* principle wearing an
//! uncitable *quantitative* coefficient — stated plainly, not disguised as
//! a Slater number.** No Slater-style rule set (checked: Slater's own 1960
//! revision, Clementi & Raimondi 1963/1967, Burns 1964, Bessis & Bessis
//! 1981) keys a screening coefficient on the screener's own angular
//! momentum; the one refinement that does key on angular momentum (Bessis &
//! Bessis) predicts the opposite sign for this effect. `BASE_SIGMA_DEEP`'s
//! own doc records how its value was found — numerical search against a
//! real, measured target (the Shannon-radius lanthanide contraction), the
//! same category of calibration `GapConsts::BASE_PROMOTION_BUDGET` already
//! uses, not a value copied from a source.
//!
//! **Independently perturbed from `sigma_near`, on its own stream
//! position** — the two constants model different physical uncertainties
//! (s/p near-shell partial screening; d/f screener penetration) with no
//! structural reason to move together. `min(near_rate(l), deep_rate)` in
//! `near_df`'s coefficient only ever differs from unconditionally applying
//! `deep_rate` there when the *candidate* is s/p (`l <= 1`): a d/f
//! candidate already has `near_rate(l) = 1.0`, and `deep_rate < 1.0`
//! always, so `min` collapses to `deep_rate` unconditionally there — the
//! same no-comparison shortcut `far_rate_df` uses below, just not written
//! that way for this term. The `min` genuinely compares two independent
//! constants only for an **s/p candidate with a d/f screener in its own
//! `near` bucket** — e.g. 6p priced with already-placed 5d electrons in
//! `near(6)` (5d fills before 6p under Madelung's `(n, n+l)` tie-break) —
//! where `sigma_near` and `sigma_deep` are unrelated and the `min` keeps
//! the tighter (smaller) of the two, never double-counting the deficit
//! upward past either bound alone.
//!
//! **Excluding the candidate's own already-placed same-subshell peers from
//! `same(n, l)` is a deliberate deviation from the most literal reading of
//! Decision 3's text, found necessary by direct numerical experiment, not
//! stated in the plan.** Including them makes a partially-filled subshell's
//! own energy rise as it fills — self-screening among its own occupants —
//! at a rate that, for the real 3d/4p transition, overtakes 3d's initial
//! energy advantage after roughly 4 of its 10 electrons, before 3d
//! finishes: 3d and 4p would trade off electron-by-electron instead of 3d
//! filling cleanly. Real atoms do not have this problem because of
//! exchange-correlation effects a screened-hydrogenic model does not
//! include; excluding same-subshell self-screening is the cheapest
//! structural fix that keeps the model's own account of *why* penetration
//! matters (`t(l)`) intact. This exclusion is no longer load-bearing for
//! *fill order* (Madelung's rule does not read subshell energy at all), but
//! it stays load-bearing for the *energy values* `e_homo`/`gap`/`zeff_outer`
//! read from — an unexcluded self-screening term would still make a
//! partially-filled subshell's own reported energy rise implausibly as it
//! fills.
//!
//! **`Z_eff >= 1` is provable, not merely observed, and the proof needs the
//! exclusion above — stated against `electronic_properties`'s real call
//! order, where the electron being placed joins `occ` *before* either probe
//! runs, not after.** Both `sigma_near` and `sigma_deep` are drawn through
//! [`crate::perturbation`] with a one-sided-downward direction, so each is
//! always in `(0, BASE] subset (0, 1]` for its own base; every other
//! coefficient in `S` (`near_rate_sp` when it resolves to `1.0`, `far_sp`'s
//! fixed `1.0`, and `t(l)`) is bounded in `[0, 1]` by construction.
//! `near_rate_df` is `min(near_rate_sp, sigma_deep)` and inherits the bound
//! the same way a `min` of two already-bounded operands always does —
//! `min` cannot exceed either. `far_rate_df` is unconditionally
//! `sigma_deep` with no runtime `min` at all (`subshell_energy`'s own
//! comment), which needs `sigma_deep < 1.0` to hold *unconditionally*
//! rather than by comparison — true because `BASE_SIGMA_DEEP < 1.0` and
//! perturbation only multiplies by a factor `<= 1`, pinned as a compile-time
//! assertion right after `OrbitalConsts`. So it is enough to show
//! `near_sp(n) + near_df(n) + far_sp(n) +
//! far_df(n) + same(n, l) <= z - 1` for whichever subshell `(n, l)` is
//! being probed — i.e. that the original, undivided `near(n) + far(n) +
//! same(n, l) <= z - 1` argument still holds, since `near_sp + near_df =
//! near` and `far_sp + far_df = far` are exact partitions of the same
//! counts by the screener's `l'`, not a different set of electrons —
//! because every term in `S` multiplies one of these five counts by a
//! coefficient `<= 1`.
//!
//! **Two probes, two different arguments for the same bound.**
//!
//! The first, `e_homo`'s probe, evaluates the electron's own subshell,
//! `(n, l) = (run.n, run.l)`. Both `near` and `far` require `n' != n`, and
//! `same` requires `l' != l`, so no electron in bucket `(run.n, run.l)` —
//! including the one just placed there — can appear in any of the three
//! counts, whatever order `occ` was updated in. The sum is therefore at
//! most `z` minus everything in that bucket, which is at least `1` (the
//! electron just placed), giving `S(n, l) <= z - 1`.
//!
//! The second, `zeff_outer`'s probe, can evaluate a *different* subshell,
//! `(outer_n, 0)`, in the d/f-block fallback (`run.n < outer_n`). There the
//! just-placed electron sits in a genuinely different bucket, and is *not*
//! excluded from `far(outer_n)` — legitimately counted as screening. The
//! bound still holds, for a Madelung-order reason rather than a screening
//! one: `(n, 0)` is always the first subshell to fill at a given `n`
//! (smallest `n + l` there), so whenever any subshell at `n = outer_n` is
//! occupied, `(outer_n, 0)` itself already holds at least one electron —
//! excluded from `near`, `far` and `same` for `(outer_n, 0)` by
//! construction, the same way the first case excludes the electron just
//! placed.
//!
//! Either way the probed bucket contributes at least `1` to `z` and `0` to
//! the sum, so `S(n, l) <= z - 1` and `Z_eff >= 1`, for both of these two
//! probes and both screening branches — not just the cases this module
//! happens to test.
//!
//! **A third probe is deliberately outside this bound, and `>= 1` is not
//! claimed for it.** `gap`'s LUMO probe, `subshell_energy(next.n, next.l, z,
//! &occ, ..)`, evaluates a subshell that has *no* electrons in `occ` yet —
//! neither probe's argument above applies, since nothing excludes it from
//! `near`/`far`/`same`, and the bound weakens to `S(n, l) <= z`, giving only
//! `Z_eff >= 0`. Measured, not merely possible: at z=2 the LUMO probe for
//! (2, 0) sees `near(2) = 2` (both 1s electrons) with nothing excluded,
//! `Z_eff = 2 - 2 * sigma_near = 0.3` at identity. Harmless where it is
//! read — `gap` only ever squares `e_lumo` inside a difference, never
//! divides by it — but a proof claiming `>= 1` "for every reachable z, n,
//! l" without this carve-out would be describing a stronger property than
//! the code actually has.
//!
//! Verified exhaustively for the frontier-subshell probe (`e_homo`'s shape)
//! by `z_eff_is_at_least_one_with_equality_at_z_one`, which mirrors
//! `electronic_properties`'s insert-then-call order rather than calling it
//! directly, and for `zeff_outer`'s own probe — including every d/f-block
//! fallback call — by `zeff_outer_is_at_least_one_through_the_real_call_path`,
//! which does call `electronic_properties` itself, across several drawn
//! `sigma_near` values, not only identity's.
//!
//! **What this costs: the 4d/5s interruption an earlier version of this
//! module treated as a genuinely emergent feature is gone.** Under
//! energy-driven ordering, a single 5s electron placed, then 4d filling
//! completely, then 5s's second electron, fell out of the formula near
//! `z = 37`. Madelung's rule never interrupts a subshell — every run in
//! [`fill`]'s output is a complete subshell, except possibly the very last
//! (`ELECTRON_CEILING` can and does cut the fold off mid-subshell: it lands
//! at `z = 250`, inside `(7, 4)`'s 18-slot capacity, with only 10 placed —
//! never observed in practice, since `n_elements`'s drawn range tops out at
//! 120, 130 electrons short of where that truncation first becomes
//! possible). That specific phenomenon did not correspond to any real
//! element's actual anomaly either (real Nb, Mo, Ru, Rh, Pd each break the
//! simple pattern their own, different way, for reasons neither version of
//! this model captures), so trading it for a fold that actually matches the
//! real table's structure through `z = 118` is the right side of that
//! trade, not a loss to mourn — but it means [`Occupancy`]'s own doc (a
//! subshell can appear as two separate, non-adjacent runs) no longer
//! describes anything this module produces within the reachable range;
//! kept accurate below, not deleted, since a future change could
//! reintroduce a fold shape where it matters again.
//!
//! This module has no production caller yet — Steps 6-9 give it one.
//! Steps 1-5 (this module) build and test the fold in isolation, per the
//! plan's own staging; `element.rs`'s `generate_elements(seed,
//! PhysicsVersion::V2)` calls [`fill`] when it lands. Every item here is
//! exercised by this module's own `#[cfg(test)]` tests, which are real
//! callers today — the blanket `#[expect]` below is scoped to non-test
//! builds for exactly that reason, matching `perturbation.rs`'s own
//! precedent.
//!
//! **`clippy::redundant_pub_crate` and rustc's own `unreachable_pub`
//! disagree about `pub(crate)` on this module's top-level items, and both
//! are enabled workspace-wide.** `redundant_pub_crate` wants `pub`
//! (`orbital` is itself `pub(crate)`, so a `pub` item inside it is already
//! capped at crate visibility); `unreachable_pub` correctly objects to that
//! spelling, since a reader seeing `pub fn fill` would reasonably expect it
//! reachable from outside the crate, which it is not. `pub(crate)` is the
//! honest spelling and matches `packing.rs`'s own established convention
//! for a `pub(crate) mod`'s items — the allow below is for
//! `redundant_pub_crate` specifically, not a blanket suppression.
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "fill(), Occupancy, OrbitalConsts, GapConsts, ElectronicProperties and \
                  electronic_properties() have no in-crate caller until the rest of Steps 6-9 \
                  wire generate_elements(seed, PhysicsVersion::V2) to this module — see this \
                  module's own doc"
    )
)]
#![allow(
    clippy::redundant_pub_crate,
    reason = "in tension with unreachable_pub for a pub(crate)-mod's own items with no \
              cross-module caller yet — see this module's own doc for which one is right and \
              why. Resolves itself once Steps 6-9 give fill() a caller in element.rs, matching \
              packing.rs's clean state today"
)]

use crate::perturbation::{Direction, MigratedConstant, Rung, draw_symmetric, perturb};
use borbax_rng::{Domain, Stream};
use std::collections::BTreeMap;

/// One contiguous run of electrons [`fill`] placed into a single subshell
/// before moving to a different one.
///
/// **A complete subshell today, except possibly the very last entry — and
/// that is new.** Under energy-driven ordering an earlier version of this
/// module could interrupt a subshell mid-fill (this module's own doc
/// records the 4d/5s case it found); [`madelung_order`] never does — every
/// subshell fills to capacity before the fold moves to the next one, so
/// `(n, l)` appears at most once in [`fill`]'s output, and every entry but
/// possibly the last is complete (`ELECTRON_CEILING` can cut the very last
/// one off mid-subshell — see this module's own doc for where, and why it
/// never surfaces in a real generated table). The field stays `count: u32`
/// rather than a fixed `capacity_of(l)`, and a reader summing every entry
/// for a given `(n, l)` still gets the right answer, because *how* this
/// module fills is not a contract the rest of the crate should have to
/// know — only a future
/// change reintroducing interruption would need this doc corrected again.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Occupancy {
    pub(crate) n: u8,
    pub(crate) l: u8,
    pub(crate) count: u32,
}

impl Occupancy {
    /// `2(2l + 1)` — the subshell's total capacity, independent of how many
    /// electrons this particular entry holds.
    #[must_use]
    pub(crate) const fn capacity(self) -> u32 {
        capacity_of(self.l)
    }

    /// The number of unpaired electrons in this entry, by Hund's rule:
    /// singly-occupying every orbital before any pairs, so `count` unpaired
    /// up to half capacity and `capacity - count` unpaired beyond it.
    ///
    /// **A pure function of `(l, count)`, spanning `l >= 1`.** `l == 0` (an
    /// s subshell, one orbital) has no room for Hund's rule to matter —
    /// `count` is already the unpaired count for every legal value — so the
    /// formula is exercised nontrivially only from `l == 1` up, which is
    /// where this module's own tests probe it.
    #[must_use]
    pub(crate) const fn unpaired_count(self) -> u32 {
        let half = capacity_of(self.l) / 2;
        if self.count <= half {
            self.count
        } else {
            capacity_of(self.l) - self.count
        }
    }
}

/// `2(2l + 1)`, as a free function so [`Occupancy::capacity`] and the fold
/// below share one definition.
#[expect(
    clippy::as_conversions,
    reason = "l is a quantum number, always small (this module's candidates go up to l=8), so \
              the widening u8 -> u32 cannot truncate or change sign"
)]
const fn capacity_of(l: u8) -> u32 {
    2 * (2 * l as u32 + 1)
}

/// `t(l) = min(l / 2, 1.0)` — see this module's own doc for why this
/// specific saturating shape, not Slater's step function.
///
/// **Written as `if`, not `f64::min`** — spec §13.1: `f64::min` is
/// non-deterministic for `±0.0` across platforms, and `l / 2.0` is never
/// negative, so this branch has no sign-of-zero hazard to inherit.
fn t_of_l(l: u8) -> f64 {
    let raw = f64::from(l) / 2.0;
    if raw < 1.0 { raw } else { 1.0 }
}

/// The physical constants [`electronic_properties`] needs to price a
/// subshell once [`madelung_order`] has chosen it — [`fill`] itself needs
/// none, see its own doc. The perturbed near-shell and deep-screener
/// screening coefficients — see [`MigratedConstant::ScreeningInner`] and
/// [`MigratedConstant::ScreeningDeep`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OrbitalConsts {
    sigma_near: f64,
    sigma_deep: f64,
}

impl OrbitalConsts {
    /// The unperturbed near-shell screening coefficient — Slater's own
    /// historical value (Slater 1930), not a value chosen for this model.
    /// See this module's own doc for the two-tier screening rule this
    /// feeds and what it was measured to reproduce.
    pub(crate) const BASE_SIGMA_NEAR: f64 = 0.85;

    /// The unperturbed deep-screener ceiling — **calibrated by numerical
    /// search, not a citable literature value** (see
    /// [`MigratedConstant::ScreeningDeep`]'s own doc for why no such value
    /// exists to adopt). Searched against the real lanthanide contraction's
    /// measured Shannon radii, La³⁺ 103.2 pm → Lu³⁺ 86.1 pm — a 16.57%
    /// contraction. `0.928` reproduces that ratio to within 0.0025
    /// percentage points across the identity fold's own lanthanide-analogue
    /// block (z = 57..=70), against a swept range `0.85..=1.00`; nearby
    /// values move the reproduced contraction by roughly 0.13 percentage
    /// points per `0.001` of `sigma_deep`, so this is a genuine optimum, not
    /// a coincidence of the sweep's granularity.
    ///
    /// **The comparison is between spans of different lengths, honestly,
    /// not silently.** The model's own 4f-filling block under this fold's
    /// strict Madelung order is 14 elements, `z = 57..=70` (5d starts at
    /// `z = 71`). The real La(57)→Lu(71) span quoted above is 15 *elements*
    /// — real La's own ground state has no 4f electron at all (`[Xe]5d¹6s²`),
    /// an anomaly this fold's idealised order does not and should not
    /// reproduce. So the target is the real span's *overall* contraction
    /// magnitude, applied as an analogue to this model's own 4f-filling
    /// block, not a claim that `z = 57` in this fold corresponds element-
    /// for-element to real La. `/review-pr` round 3 measured the
    /// alternative framings this ambiguity admits: fitting the model's
    /// `57..=71` span (extending one element to match the real range
    /// length) against the same 16.57% target gives `sigma_deep ≈ 0.945`;
    /// fitting `57..=70` against the real La→Yb span's own 15.89%
    /// contraction instead gives `≈ 0.933`. Both are within the same
    /// single-digit-percent neighbourhood as `0.928` and do not change
    /// which qualitative conclusion this constant supports (real,
    /// non-trivial lanthanide-analogue contraction where none existed
    /// before); the shipped value is not re-derived from either alternative
    /// mapping, since nothing here turns on the third decimal place.
    pub(crate) const BASE_SIGMA_DEEP: f64 = 0.928;

    /// Draw this universe's screening constants.
    ///
    /// **Both one-sided downward, per the plan's own text for
    /// `sigma_near` and the same reasoning extended to `sigma_deep`** — `p`
    /// sampled from `[-Direction::P_MAX, 0.0]`, never the symmetric range.
    /// This is a draw-site decision, not a [`MigratedConstant`]
    /// classification tag: see [`MigratedConstant::ScreeningInner`]'s own
    /// doc for why it is not generalised into a fourth tag. Two independent
    /// draws, two independent stream positions — `sigma_near` and
    /// `sigma_deep` model different physical uncertainties (s/p near-shell
    /// screening versus d/f screener penetration) and are not assumed to
    /// move together.
    #[must_use]
    pub(crate) fn draw(seed: u64, rung: Rung) -> Self {
        let mut near_stream = Stream::new(
            seed,
            Domain::Perturbation,
            MigratedConstant::ScreeningInner.index(),
        );
        let near_p = near_stream.next_f64_range(-Direction::P_MAX, 0.0);
        let near_direction = Direction::new(near_p, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("p is drawn within [-P_MAX, 0.0] by construction"));

        let mut deep_stream = Stream::new(
            seed,
            Domain::Perturbation,
            MigratedConstant::ScreeningDeep.index(),
        );
        let deep_p = deep_stream.next_f64_range(-Direction::P_MAX, 0.0);
        let deep_direction = Direction::new(deep_p, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("p is drawn within [-P_MAX, 0.0] by construction"));

        Self {
            sigma_near: perturb(Self::BASE_SIGMA_NEAR, rung, near_direction),
            sigma_deep: perturb(Self::BASE_SIGMA_DEEP, rung, deep_direction),
        }
    }

    #[cfg(test)]
    pub(crate) const fn at_identity() -> Self {
        Self {
            sigma_near: Self::BASE_SIGMA_NEAR,
            sigma_deep: Self::BASE_SIGMA_DEEP,
        }
    }
}

// `subshell_energy`'s `far_rate_df` is unconditionally `consts.sigma_deep`,
// with no runtime comparison against `1.0` — that shortcut, and the
// `Z_eff >= 1` proof this module's own doc makes, both depend on
// `BASE_SIGMA_DEEP < 1.0` holding for every possible perturbed value, which
// only holds because the *base* itself is `< 1.0` (perturbation multiplies
// by a factor `<= 1` for a one-sided-downward draw). Nothing else in the
// type system pins that; a future recalibration of `BASE_SIGMA_DEEP` to
// `>= 1.0` would silently break both without this failing to compile.
const _: () = assert!(
    OrbitalConsts::BASE_SIGMA_DEEP < 1.0,
    "far_rate_df (subshell_energy) and this module's Z_eff >= 1 proof both rely on \
     BASE_SIGMA_DEEP < 1.0 holding unconditionally"
);

/// The highest electron count [`fill`] computes to. Comfortably above the
/// naming grammar's 165-element symbol-space ceiling (`naming.rs`), so
/// every legal `n_elements` draw has a complete answer.
const ELECTRON_CEILING: u32 = 250;

/// The highest principal quantum number [`fill`] considers as a candidate.
/// `n = 9` reaches subshells well past [`ELECTRON_CEILING`]'s last electron,
/// so the fold never runs out of candidates before the ceiling.
const MAX_N: u8 = 9;

/// Every `(n, l)` subshell up to [`MAX_N`], in Madelung order: ascending
/// `(n + l)`, ties broken by ascending `n` — this module's own doc explains
/// why fill order is a fixed structural rule rather than a per-candidate
/// energy comparison. A pure function of nothing but [`MAX_N`]: the same
/// sequence for every seed and every universe.
fn madelung_order() -> Vec<(u8, u8)> {
    let mut out = Vec::new();
    for n in 1..=MAX_N {
        for l in 0..n {
            out.push((n, l));
        }
    }
    out.sort_by_key(|&(n, l)| (n + l, n));
    out
}

/// The `(energy, Z_eff)` pair for subshell `(n, l)` given electrons already
/// placed, at electron count `z` — the fold's own energy formula (this
/// module's doc), shared by [`electronic_properties`]'s own callers below
/// rather than duplicated.
///
/// **Deliberately takes `occ` as-is, with no capacity check.** Whether
/// `(n, l)` still has room is the caller's question — [`fill`] never asks
/// this function about a full subshell; [`electronic_properties`] calls it
/// *because* a subshell is full, to price the next candidate.
fn subshell_energy(
    n: u8,
    l: u8,
    z: u32,
    occ: &BTreeMap<(u8, u8), u32>,
    consts: OrbitalConsts,
) -> (f64, f64) {
    let near_sp: u32 = occ
        .iter()
        .filter(|((n2, l2), _)| n.checked_sub(1) == Some(*n2) && *l2 <= 1)
        .map(|(_, c)| c)
        .sum();
    let near_df: u32 = occ
        .iter()
        .filter(|((n2, l2), _)| n.checked_sub(1) == Some(*n2) && *l2 >= 2)
        .map(|(_, c)| c)
        .sum();
    let far_sp: u32 = occ
        .iter()
        .filter(|((n2, l2), _)| n.checked_sub(2).is_some_and(|floor| *n2 <= floor) && *l2 <= 1)
        .map(|(_, c)| c)
        .sum();
    let far_df: u32 = occ
        .iter()
        .filter(|((n2, l2), _)| n.checked_sub(2).is_some_and(|floor| *n2 <= floor) && *l2 >= 2)
        .map(|(_, c)| c)
        .sum();
    let same: u32 = occ
        .iter()
        .filter(|((n2, l2), _)| *n2 == n && *l2 != l)
        .map(|(_, c)| c)
        .sum();
    // l <= 1: near screens at sigma_near (partial); l >= 2: near screens at
    // full strength too — this module's own doc on why that split, not a
    // single constant, is what makes 6s win against 4f/5d.
    let near_rate_sp = if l <= 1 { consts.sigma_near } else { 1.0 };
    // **`sigma_deep` caps every coefficient a d/f SCREENER contributes,
    // regardless of the CANDIDATE's own l — this module's own doc on why
    // the imperfection belongs on the screener's side, not the candidate's.**
    // `sigma_deep`'s own one-sided-downward draw, base `BASE_SIGMA_DEEP <
    // 1.0`, guarantees `sigma_deep < 1.0` always (perturbation only ever
    // multiplies the base by a factor `<= 1`) — so `far_rate_df` is
    // unconditionally `sigma_deep`, not a real comparison against `1.0`.
    // `near_rate_df` genuinely does compare two independently-drawn
    // constants for an `l <= 1` candidate (`sigma_near` versus
    // `sigma_deep`, neither bounded relative to the other), which is why
    // it, unlike `far_rate_df`, needs the comparison written out. The
    // comparison form, not `f64::min` — §13.1 bans the method as
    // non-deterministic for ±0.0; not reachable here (both operands are
    // always finite and strictly positive), but the crate-wide rule
    // applies regardless of whether a given call site can hit the case it
    // guards against.
    let near_rate_df = if near_rate_sp < consts.sigma_deep {
        near_rate_sp
    } else {
        consts.sigma_deep
    };
    let far_rate_df = consts.sigma_deep;
    let s = near_rate_sp * f64::from(near_sp)
        + near_rate_df * f64::from(near_df)
        + f64::from(far_sp)
        + far_rate_df * f64::from(far_df)
        + t_of_l(l) * f64::from(same);
    let z_eff = f64::from(z) - s;
    let ratio = z_eff / f64::from(n);
    (-(ratio * ratio), z_eff)
}

/// Fill subshells electron by electron, in [`madelung_order`], up to
/// [`ELECTRON_CEILING`].
///
/// **Takes no [`OrbitalConsts`] — deliberately, and that is new.** Which
/// subshell fills next no longer depends on screening at all (this module's
/// own doc), so the run structure this returns is identical for every
/// seed; only the *energies* [`electronic_properties`] later reads from it
/// vary. Returns one [`Occupancy`] per subshell — see [`Occupancy`]'s own
/// doc for why that is now always a complete subshell, never a fragment.
#[must_use]
pub(crate) fn fill() -> Vec<Occupancy> {
    let mut out: Vec<Occupancy> = Vec::new();
    let mut placed: u32 = 0;
    for (n, l) in madelung_order() {
        if placed >= ELECTRON_CEILING {
            break;
        }
        let take = capacity_of(l).min(ELECTRON_CEILING - placed);
        placed += take;
        out.push(Occupancy { n, l, count: take });
    }
    out
}

/// Task 26.1 Steps 6-9's gap-derived constants (issue #26, Decision 8) —
/// currently just the valence-promotion energy budget. See
/// [`MigratedConstant::PromotionBudget`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct GapConsts {
    promotion_budget: f64,
}

impl GapConsts {
    /// Calibrated by direct numerical search — unlike
    /// `OrbitalConsts::BASE_SIGMA_NEAR`, which is Slater's own historical
    /// value rather than a value chosen for this model — to sit strictly
    /// between real alkaline-earth-like HOMO-LUMO gaps (~1.1-1.6 at the
    /// identity configuration, e.g. beryllium- and magnesium-analogues) and
    /// real noble-gas-like ones (~3.9 and up, e.g. helium- and
    /// neon-analogues) — a >25%-margin choice on both sides, not a boundary
    /// value.
    pub(crate) const BASE_PROMOTION_BUDGET: f64 = 2.0;

    /// Draw this universe's gap-derived constants. Symmetric excursion,
    /// unlike [`OrbitalConsts::draw`]'s one-sided-downward screening
    /// coefficient — see [`MigratedConstant::PromotionBudget`]'s own doc
    /// for why no structural reason favours one direction here.
    #[must_use]
    pub(crate) fn draw(seed: u64, rung: Rung) -> Self {
        Self {
            promotion_budget: draw_symmetric(
                seed,
                rung,
                MigratedConstant::PromotionBudget,
                Self::BASE_PROMOTION_BUDGET,
            ),
        }
    }

    #[cfg(test)]
    pub(crate) const fn at_identity() -> Self {
        Self {
            promotion_budget: Self::BASE_PROMOTION_BUDGET,
        }
    }
}

/// One element's derived electronic-structure properties at the z-electron
/// configuration (issue #26, Task 26.1 Steps 6-9, Decision 8) — the shared
/// "gap" computation valence/affinity/period-boundaries/bond-capacity all
/// read from, with one call site each.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ElectronicProperties {
    /// The outermost occupied shell — `max(n)` over every subshell with a
    /// nonzero count, **not** the last-placed electron's own `n`, which can
    /// regress inside the d/f-block (this module's own doc records why, and
    /// the bug it caused when radius was keyed on the wrong one). Monotone
    /// non-decreasing in z by construction — a max over an only-ever-growing
    /// occupied set — which is what makes `outer_n(z) > outer_n(z - 1)` a
    /// threshold-free period-boundary test.
    pub(crate) outer_n: u8,
    /// `Z_eff` of `(outer_n, l)` at this z, `l` the frontier's own angular
    /// momentum when the frontier sits in the outer shell and `0` only as
    /// the d/f-block fallback — see [`electronic_properties`]'s own probe
    /// site for why, and this module's own doc for what fixing the
    /// unconditional `l = 0` probe measurably changed. Feeds radius.
    pub(crate) zeff_outer: f64,
    /// The HOMO-LUMO gap. Zero whenever `fill()`'s own run structure shows
    /// this electron's subshell is not yet at capacity; nonzero only when
    /// this electron exactly completes its subshell **and** a different
    /// subshell follows.
    ///
    /// **Which case applies is read from the fold's actual output, not
    /// re-derived from `(z, n, l)` alone, even though `madelung_order()`
    /// alone now fixes *which* position is subshell-final — screening plays
    /// no part in that any more (this module's own doc: an earlier,
    /// energy-driven ordering could interrupt a partially-filled subshell,
    /// e.g. real 5s at z=37→38; Madelung's rule never does).** The value at
    /// a subshell-final position still cannot be re-derived independently:
    /// it is the *following* subshell's own `Z_eff`-dependent energy, which
    /// only [`subshell_energy`] can supply, so reading `gap` from the same
    /// pass that already computed it avoids a second, redundant call rather
    /// than avoiding any real unpredictability in the structure itself.
    pub(crate) gap: f64,
    /// The frontier subshell's own binding energy, `-(Z_eff/n)^2` — the
    /// electron most recently placed, evaluated at this z.
    pub(crate) e_homo: f64,
    /// Unpaired count of the frontier subshell, or the promotion bonus (2 —
    /// moving one electron out of a doubly-occupied orbital into an empty
    /// one always creates exactly two unpaired sites, regardless of which
    /// subshell) if the frontier just closed and `gap <= promotion_budget`.
    pub(crate) valence: u32,
    /// The frontier subshell's own `l` — the electron most recently
    /// placed's angular momentum, not `outer_n`'s. Element.rs's V2 `group`/
    /// `block()` read this directly; it is a structural readout of the
    /// fold's own state, not a fitted quantity.
    pub(crate) frontier_l: u8,
    /// Paired-electron count across **every** currently-occupied subshell,
    /// `(z - total_unpaired) / 2` — element.rs's V2 mass-defect "contacts"
    /// analogue. **Deliberately summed over all of `occ`, not just the
    /// frontier subshell**, even though [`fill`] never leaves more than one
    /// subshell incomplete at a time today — a closed subshell always
    /// contributes 0 unpaired electrons regardless, so this is the honest
    /// general form rather than one that happens to agree with a
    /// frontier-only count only because of how this module currently fills.
    pub(crate) paired_count: u32,
}

/// Compute [`ElectronicProperties`] for every z from 1 to `runs`' own total
/// electron count, by walking [`fill`]'s run output in one linear pass — no
/// re-search.
///
/// **Why this is a separate pass over `fill`'s output rather than a second
/// fold**: `fill` is already tested and shipped (Steps 1-5); this only
/// *reads* its run boundaries (which subshell a given z belongs to, and
/// whether it exactly completes that subshell) rather than re-deriving the
/// aufbau choice, so there remains exactly one fold implementation to keep
/// correct.
#[must_use]
pub(crate) fn electronic_properties(
    runs: &[Occupancy],
    orbital_consts: OrbitalConsts,
    gap_consts: GapConsts,
) -> Vec<ElectronicProperties> {
    let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
    let mut out = Vec::new();
    let mut z: u32 = 0;
    for (run_index, run) in runs.iter().enumerate() {
        for _ in 0..run.count {
            z += 1;
            *occ.entry((run.n, run.l)).or_insert(0) += 1;

            let outer_n = occ
                .iter()
                .filter(|&(_, &c)| c > 0)
                .map(|(&(n, _), _)| n)
                .max()
                .unwrap_or_else(|| unreachable!("(run.n, run.l) was just inserted with count > 0"));
            // **The frontier's own `l` when the frontier sits in the outer
            // shell, `l = 0` only as a fallback for the d/f-block region
            // (`run.n < outer_n`) where no better probe subshell exists.**
            // Probing `l = 0` unconditionally (an earlier version of this
            // line) zeroes same-shell screening for every candidate whose
            // real frontier is p/d/f (`t_of_l(0) == 0`), which is finding
            // F4 from the 2026-08-10 `/review-pr` pass: it over-contracted
            // radius within a period (Li->Ne 5.12x modelled against 2.21x
            // real) and could invert the group trend. This module's own
            // doc records the measured fix: correlation against real
            // covalent radii 0.50 -> 0.87, group 1 (Li..Fr) radius trend
            // inverted -> strictly increasing.
            let probe_l = if run.n == outer_n { run.l } else { 0 };
            let (_, zeff_outer) = subshell_energy(outer_n, probe_l, z, &occ, orbital_consts);
            let (e_homo, _) = subshell_energy(run.n, run.l, z, &occ, orbital_consts);

            let count_now = occ
                .get(&(run.n, run.l))
                .copied()
                .unwrap_or_else(|| unreachable!("just inserted above"));
            let subshell_full = count_now == capacity_of(run.l);
            // **Clamped at zero, not `.abs()` — `.abs()` was a `/review-pr`
            // finding (2026-08-10): it reads as a distance but is used as a
            // cost, and those disagree exactly where it matters.** At a
            // transition like 6s -> 4f (z=56), the *structurally* next
            // subshell has a more negative raw energy than the one that
            // just closed (a small n dominates `-(Z_eff/n)^2` regardless of
            // screening, this module's own doc) — not a bug in either
            // formula, but a consequence of `subshell_energy` no longer
            // being the thing that chose this order. `gap` is read
            // downstream as the energy *cost* of promoting an electron out
            // of a just-closed subshell (`valence = 2 if gap <=
            // promotion_budget`): a negative raw difference means the next
            // state is *lower*, so the promotion is free, and the correct
            // cost is `0.0`, not `|raw difference|`. `.abs()` read "free" as
            // "maximally expensive" instead, and did so at exactly five
            // positions in `1..=250` (z = 56, 88, 120, 170, 218 — measured):
            // every one is a closed `s^2` whose Madelung-next subshell is a
            // low-n, high-l state, and every one is where the fold's fixed
            // order and `subshell_energy`'s raw energies structurally
            // disagree. Barium-, radium- and (at `n_elements = 120`
            // specifically) the z=120 analogue are the three reachable in a
            // real generated table (`n_elements` draws `60..=120`, not
            // `1..=118` — 170 and 218 are the two genuinely unreachable
            // ones), and `.abs()` gave all three valence 0 instead of the
            // physically correct 2 — noble-gas-inert alkaline earths, in
            // every V2 universe, unconditionally.
            let gap = if subshell_full {
                runs.get(run_index + 1).map_or(0.0, |next| {
                    let (e_lumo, _) = subshell_energy(next.n, next.l, z, &occ, orbital_consts);
                    let d = e_lumo - e_homo;
                    if d > 0.0 { d } else { 0.0 }
                })
            } else {
                0.0
            };

            let unpaired = Occupancy {
                n: run.n,
                l: run.l,
                count: count_now,
            }
            .unpaired_count();
            let valence = if unpaired > 0 {
                unpaired
            } else if gap <= gap_consts.promotion_budget {
                2
            } else {
                0
            };

            let total_unpaired: u32 = occ
                .iter()
                .map(|(&(n, l), &count)| Occupancy { n, l, count }.unpaired_count())
                .sum();
            let paired_count = (z - total_unpaired) / 2;

            out.push(ElectronicProperties {
                outer_n,
                zeff_outer,
                gap,
                e_homo,
                valence,
                frontier_l: run.l,
                paired_count,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::{
        ELECTRON_CEILING, ElectronicProperties, GapConsts, MAX_N, Occupancy, OrbitalConsts,
        capacity_of, electronic_properties, fill, madelung_order, subshell_energy,
    };
    use std::collections::BTreeMap;

    /// The real subshell fill order through `z = 118` — Madelung's rule,
    /// which needs no screening input (this module's own doc), pinned as a
    /// computed constant (CLAUDE.md, §5-permitted under Decision 2), not a
    /// lookup table. Extended from the pre-fix `through 5p` (11 entries) to
    /// the full table once Madelung ordering replaced the energy-driven
    /// search that could not reach past it.
    const IDENTITY_FILL_ORDER: [(u8, u8); 19] = [
        (1, 0),
        (2, 0),
        (2, 1),
        (3, 0),
        (3, 1),
        (4, 0),
        (3, 2),
        (4, 1),
        (5, 0),
        (4, 2),
        (5, 1),
        (6, 0),
        (4, 3),
        (5, 2),
        (6, 1),
        (7, 0),
        (5, 3),
        (6, 2),
        (7, 1),
    ];

    #[test]
    fn subshell_capacity_is_2_times_2l_plus_1() {
        for l in 0..=6u8 {
            assert_eq!(capacity_of(l), 2 * (2 * u32::from(l) + 1), "l={l}");
        }
    }

    #[test]
    fn hunds_rule_maximizes_unpaired_across_l_geq_1() {
        for l in 1..=3u8 {
            let cap = capacity_of(l);
            for count in 0..=cap {
                let occ = Occupancy { n: 4, l, count };
                let want = if count <= cap / 2 { count } else { cap - count };
                assert_eq!(
                    occ.unpaired_count(),
                    want,
                    "l={l} count={count}: unpaired should follow Hund's rule"
                );
            }
        }
    }

    /// [`madelung_order`] is a pure, deterministic sort — the discriminator
    /// this replaces (`fill_order_is_invariant_to_candidate_shuffle`,
    /// through the pre-Madelung `best_unfilled`'s own energy search) no
    /// longer applies: there is no scan order left to be invariant to, since
    /// picking the next subshell is no longer a search at all. What is
    /// still worth pinning: the sequence is genuinely sorted by `(n + l,
    /// n)`, not merely equal to [`IDENTITY_FILL_ORDER`]'s own prefix (which
    /// the next test already checks) — this catches a `sort_by_key` swapped
    /// for a stable-looking but wrong comparator.
    #[test]
    fn madelung_order_is_sorted_by_n_plus_l_then_n() {
        let order = madelung_order();
        for pair in order.windows(2) {
            let [a, b] = pair else {
                unreachable!("windows(2) always yields a 2-element slice")
            };
            let (a, b) = (*a, *b);
            let key = |(n, l): (u8, u8)| (n + l, n);
            assert!(
                key(a) <= key(b),
                "{a:?} (key {:?}) should sort at or before {b:?} (key {:?})",
                key(a),
                key(b)
            );
        }
        // And every legal (n, l) with n <= MAX_N appears exactly once —
        // sorting cannot silently drop, duplicate or substitute an entry.
        // Comparing lengths before and after dedup only catches duplicates:
        // a *missing* entry shrinks both sides equally and passes silently.
        // Comparing the actual set against the pairs the definition of
        // legal (n, l) generates catches all three at once.
        let expected: std::collections::BTreeSet<(u8, u8)> = (1..=MAX_N)
            .flat_map(|n| (0..n).map(move |l| (n, l)))
            .collect();
        let seen: std::collections::BTreeSet<(u8, u8)> = order.iter().copied().collect();
        assert_eq!(
            seen, expected,
            "madelung_order() should contain exactly every (n, l) with 1 <= n <= MAX_N, 0 <= l < n"
        );
        assert_eq!(
            order.len(),
            seen.len(),
            "madelung_order() should contain no duplicates"
        );
    }

    #[test]
    fn subshell_closures_are_running_partial_sums() {
        // Derived from IDENTITY_FILL_ORDER's own (n, l) pairs via capacity_of,
        // never a hardcoded [2, 4, 10, 12, 18, 20, 30, 36, ...] literal.
        let expected_running: Vec<u32> = IDENTITY_FILL_ORDER
            .iter()
            .scan(0u32, |sum, &(_, l)| {
                *sum += capacity_of(l);
                Some(*sum)
            })
            .collect();

        let order = fill();
        // Every entry is a completed subshell today (madelung_order() never
        // interrupts one, this module's own doc), so all of
        // IDENTITY_FILL_ORDER lines up 1:1 with fill()'s own output, not
        // just a prefix.
        let mut running = 0u32;
        let mut checked = 0usize;
        for (occ, want) in order
            .iter()
            .take(IDENTITY_FILL_ORDER.len())
            .zip(expected_running.iter())
        {
            running += occ.count;
            assert_eq!(
                occ.count,
                occ.capacity(),
                "entry {checked} should be a completed subshell"
            );
            assert_eq!(
                running, *want,
                "entry {checked} ({}, {}): running total should match the derived closure",
                occ.n, occ.l
            );
            checked += 1;
        }
        assert_eq!(
            checked,
            IDENTITY_FILL_ORDER.len(),
            "should have checked every entry"
        );
    }

    /// **Every run but possibly the last is a complete subshell — the claim
    /// `Occupancy`'s own doc makes, settled directly rather than left as an
    /// assertion nobody checks.** `ELECTRON_CEILING` (250) lands inside
    /// `(7, 4)`'s 18-slot capacity with only 10 placed, so the fold's full
    /// output genuinely has one partial entry; this pins that it is exactly
    /// one, and exactly the last.
    #[test]
    fn every_run_but_possibly_the_last_is_a_complete_subshell() {
        let runs = fill();
        let (last, rest) = runs
            .split_last()
            .unwrap_or_else(|| unreachable!("fill() always places at least one electron"));
        for r in rest {
            assert_eq!(
                r.count,
                r.capacity(),
                "({}, {}) should be a completed subshell",
                r.n,
                r.l
            );
        }
        assert!(
            last.count <= last.capacity(),
            "({}, {}): count {} exceeds capacity {}",
            last.n,
            last.l,
            last.count,
            last.capacity()
        );
        assert!(
            last.count < last.capacity(),
            "the last run is exactly complete at this ELECTRON_CEILING — update this test's \
             own doc, which claims it is genuinely partial"
        );
        assert_eq!(
            (last.n, last.l, last.count),
            (7, 4, 10),
            "the truncated run moved — ELECTRON_CEILING or MAX_N changed, or the fold's \
             structure did; update this test's own doc to match"
        );
    }

    #[test]
    fn the_identity_fill_order_matches_the_real_periodic_table_through_z_118() {
        let order = fill();
        let touched: Vec<(u8, u8)> = order
            .iter()
            .take(IDENTITY_FILL_ORDER.len())
            .map(|o| (o.n, o.l))
            .collect();
        assert_eq!(
            touched, IDENTITY_FILL_ORDER,
            "the identity configuration's fill order should match the real periodic table \
             through z = 118"
        );
    }

    /// `t(l)` is no longer load-bearing for fill *order* — Madelung's rule
    /// reads no screening input at all (this module's own doc) — but it
    /// stays load-bearing for the *energy values* `subshell_energy` reports,
    /// both through `same(n, l)`'s own coefficient and through
    /// [`subshell_energy`]'s near/far split, which routes `l <= 1` and
    /// `l >= 2` candidates through different screening rates. This probes
    /// the value-level claim directly: a p subshell (l=1, t=0.5) and a d
    /// subshell (l=2, t=1.0, and the near-tier switch) at the same `(n, z,
    /// occ)` must report different `Z_eff`, or `t(l)` and the near/far split
    /// are both decorative.
    #[test]
    fn t_of_l_and_the_near_far_split_change_the_reported_z_eff() {
        let consts = OrbitalConsts::at_identity();
        let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        occ.insert((1, 0), 2);
        occ.insert((2, 0), 2);
        occ.insert((2, 1), 4);
        occ.insert((3, 0), 1);
        occ.insert((3, 1), 1);
        let z = 11;
        let (_, zeff_p) = subshell_energy(3, 1, z, &occ, consts);
        let (_, zeff_d) = subshell_energy(3, 2, z, &occ, consts);
        assert!(
            (zeff_p - zeff_d).abs() > 1e-9,
            "a p candidate (zeff={zeff_p}) and a d candidate (zeff={zeff_d}) at the same (n, z, \
             occ) reported the same Z_eff — t(l) and the near/far split are not distinguishing them"
        );
    }

    #[test]
    fn z_eff_is_at_least_one_with_equality_at_z_one() {
        let consts = OrbitalConsts::at_identity();
        let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        let mut z: u32 = 0;
        for run in fill() {
            for _ in 0..run.count {
                z += 1;
                // PRODUCTION'S REAL ORDER (`electronic_properties`): the
                // current electron joins `occ` *before* `subshell_energy` is
                // called, for both the frontier probe and the outer-shell
                // probe. Testing the call-before-insert order here would not
                // discriminate anything real.
                *occ.entry((run.n, run.l)).or_insert(0) += 1;
                let (_, z_eff) = subshell_energy(run.n, run.l, z, &occ, consts);
                assert!(
                    z_eff >= 1.0,
                    "z={z} ({},{}): Z_eff {z_eff} fell below 1.0",
                    run.n,
                    run.l
                );
                if z == 1 {
                    assert!(
                        (z_eff - 1.0).abs() < f64::EPSILON,
                        "Z_eff at Z=1 should be exactly 1.0"
                    );
                }
            }
        }
        assert_eq!(z, ELECTRON_CEILING, "should have walked every electron");
    }

    #[test]
    fn the_screening_accumulation_is_order_invariant_and_exact() {
        // Two different insertion orders reaching the same final occupancy
        // (through argon, z=18): the real fold's own run structure (order
        // A), and the same final per-subshell counts inserted in reverse
        // key order into a fresh map (order B) — BTreeMap's own iteration
        // order is sorted regardless of insertion order, so this exercises
        // subshell_energy's summation, not the map's.
        let consts = OrbitalConsts::at_identity();
        let mut occ_a: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        let mut z: u32 = 0;
        for run in fill() {
            for _ in 0..run.count {
                z += 1;
                *occ_a.entry((run.n, run.l)).or_insert(0) += 1;
                if z == 18 {
                    break;
                }
            }
            if z == 18 {
                break;
            }
        }
        assert_eq!(z, 18, "the fixture should stop exactly at argon");

        let mut occ_b: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        for (&key, &count) in occ_a.iter().rev() {
            occ_b.insert(key, count);
        }
        assert_eq!(
            occ_a, occ_b,
            "both occupancy maps should hold the same final counts"
        );

        // The next candidate in Madelung order after argon is potassium's
        // own 4s — read from IDENTITY_FILL_ORDER rather than hardcoded, so
        // a reordering upstream cannot silently desync this fixture from
        // what fill() actually produces.
        let (n, l) = IDENTITY_FILL_ORDER
            .get(5)
            .copied()
            .unwrap_or_else(|| unreachable!("IDENTITY_FILL_ORDER has 19 entries"));
        assert_eq!(
            (n, l),
            (4, 0),
            "sixth entry should be 4s, per Madelung order"
        );

        let (_, z_eff_a) = subshell_energy(n, l, 19, &occ_a, consts);
        let (_, z_eff_b) = subshell_energy(n, l, 19, &occ_b, consts);
        assert_eq!(
            z_eff_a.to_bits(),
            z_eff_b.to_bits(),
            "Z_eff should be bit-identical regardless of insertion order into the occupancy map"
        );
    }

    #[test]
    fn the_screening_coefficient_round_trips_through_perturb_at_identity() {
        // The new Domain::Perturbation consumer: at rung == 0, the drawn
        // coefficient equals the unperturbed base value bit-for-bit,
        // matching every other migrated constant's own required test.
        let consts = OrbitalConsts::draw(7, crate::perturbation::Rung::IDENTITY);
        assert_eq!(
            consts,
            OrbitalConsts::at_identity(),
            "at rung == 0, the drawn screening coefficient should equal the base value exactly"
        );
    }

    #[test]
    fn fill_reaches_the_electron_ceiling() {
        let total: u32 = fill().iter().map(|o| o.count).sum();
        assert_eq!(total, ELECTRON_CEILING);
    }

    #[test]
    fn the_promotion_budget_round_trips_through_perturb_at_identity() {
        let consts = GapConsts::draw(7, crate::perturbation::Rung::IDENTITY);
        assert_eq!(
            consts,
            GapConsts::at_identity(),
            "at rung == 0, the drawn promotion budget should equal the base value exactly"
        );
    }

    #[test]
    fn electronic_properties_covers_every_placed_electron() {
        let runs = fill();
        let total: u32 = runs.iter().map(|o| o.count).sum();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        assert_eq!(
            props.len(),
            usize::try_from(total).unwrap_or_else(|_| unreachable!("ELECTRON_CEILING fits usize")),
            "one ElectronicProperties entry per placed electron"
        );
    }

    /// Decision 8's headline claim, checked directly rather than trusted:
    /// "unpaired states reachable within an energy budget" reproduces real
    /// ground-state valences for the two periods this model was calibrated
    /// against (Steps 1-5's own fill-order pin goes through 5p), with zero
    /// per-element branching.
    /// `props`, indexed by z-1, with clippy's `indexing_slicing` (denied,
    /// reaching inside `#[cfg(test)]` per CLAUDE.md) satisfied by `.get`
    /// rather than a bare `[]` — every call site below passes a z already
    /// known to be in range, so the fallback is unreachable, not a real
    /// error path.
    fn at(props: &[ElectronicProperties], z: usize) -> ElectronicProperties {
        props
            .get(z - 1)
            .copied()
            .unwrap_or_else(|| unreachable!("z={z} should be within this test's own fixture"))
    }

    #[test]
    fn valence_matches_real_ground_state_valences_through_period_3() {
        let runs = fill();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        // H,He,Li,Be,B,C,N,O,F,Ne, Na,Mg,Al,Si,P,S,Cl,Ar, K,Ca -- z=1..=20.
        let want: [u32; 20] = [1, 0, 1, 2, 1, 2, 3, 2, 1, 0, 1, 2, 1, 2, 3, 2, 1, 0, 1, 2];
        for (i, &w) in want.iter().enumerate() {
            let z = i + 1;
            assert_eq!(
                at(&props, z).valence,
                w,
                "z={z}: valence should match the real ground-state count"
            );
        }
    }

    /// The mechanism behind the valence test above, checked at the level of
    /// the gap itself rather than only its consequence — a planted defect
    /// that priced the wrong subshell could still coincidentally reproduce
    /// the right valence numbers while getting the gap magnitude wrong.
    #[test]
    fn promotion_budget_separates_alkaline_earth_from_noble_gas_analogues() {
        let runs = fill();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        for &z in &[2usize, 10, 18] {
            let gap = at(&props, z).gap;
            assert!(
                gap > GapConsts::BASE_PROMOTION_BUDGET,
                "z={z}: noble-gas-like gap ({gap}) should exceed the promotion budget"
            );
            assert_eq!(at(&props, z).valence, 0, "z={z}: should not promote");
        }
        for &z in &[4usize, 12, 20] {
            let gap = at(&props, z).gap;
            assert!(
                gap <= GapConsts::BASE_PROMOTION_BUDGET,
                "z={z}: alkaline-earth-like gap ({gap}) should be within the promotion budget"
            );
            assert_eq!(
                at(&props, z).valence,
                2,
                "z={z}: should promote to valence 2"
            );
        }
    }

    #[test]
    fn outer_n_is_monotone_and_jumps_at_period_starts() {
        let runs = fill();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        let ceiling = usize::try_from(ELECTRON_CEILING)
            .unwrap_or_else(|_| unreachable!("ELECTRON_CEILING is 250, well within usize"));
        let mut prev = at(&props, 1).outer_n;
        for z in 2..=ceiling {
            let outer_n = at(&props, z).outer_n;
            assert!(outer_n >= prev, "z={z}: outer_n must never decrease");
            prev = outer_n;
        }
        // Real period starts at the identity configuration, through period
        // 7 (z=87, francium-analogue) — extended from period 5 (z=37) once
        // Madelung ordering (this module's own doc) closed the gap that
        // used to put period 6's own start at z=79 instead of z=55.
        for &z in &[3usize, 11, 19, 37, 55, 87] {
            assert!(
                at(&props, z).outer_n > at(&props, z - 1).outer_n,
                "z={z}: should start a new outer shell"
            );
        }
    }

    /// G3's real-chemistry correspondence, checked at the level radius is
    /// actually built from: within a period `zeff_outer` should rise (so
    /// `outer_n^2 / zeff_outer`, radius, falls), and jump down at a new
    /// period's start (so radius jumps up).
    #[test]
    fn zeff_outer_rises_within_a_period() {
        let runs = fill();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        // Period 2 spans z=3..=10: outer_n stays 2 throughout, and
        // zeff_outer should rise strictly.
        // **`>=`, not `>`, since the frontier-aware probe (this module's own
        // doc, finding F4) — the probe subshell itself changes at z=5 (2s
        // completes, 2p's own frontier starts, so `probe_l` switches from 0
        // to 1), and at the identity configuration this is an *exact* tie
        // (both 2.3, hand-verified: the +1 electron and the newly-counted
        // 2s `same`-shell screening term cancel precisely), never a
        // decrease.** A genuine decrease anywhere in this range is still a
        // real regression and still fails this assertion.
        for z in 4..=10 {
            assert_eq!(at(&props, z).outer_n, 2, "z={z}: should stay in period 2");
            assert!(
                at(&props, z).zeff_outer >= at(&props, z - 1).zeff_outer,
                "z={z}: zeff_outer should not fall within a period"
            );
        }
        // And the period genuinely contracts end to end, which `>=` alone
        // does not pin — a formula that returned a flat constant across the
        // whole period would pass the loop above and fail only here.
        assert!(
            at(&props, 10).zeff_outer > at(&props, 3).zeff_outer,
            "zeff_outer should be strictly higher at the end of period 2 than at its start"
        );
        // z=11 (Na) starts period 3: outer_n jumps, so zeff_outer -- and
        // therefore radius -- should NOT continue period 2's own trend.
        assert_eq!(at(&props, 11).outer_n, 3, "z=11 should be in a new period");
    }

    #[test]
    fn gap_and_e_homo_are_always_finite() {
        let runs = fill();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        for (i, p) in props.iter().enumerate() {
            assert!(p.gap.is_finite(), "z={}: gap should be finite", i + 1);
            assert!(p.e_homo.is_finite(), "z={}: e_homo should be finite", i + 1);
            assert!(p.gap >= 0.0, "z={}: gap should never be negative", i + 1);
        }
    }

    /// **`zeff_outer >= 1.0` through `electronic_properties`'s real call
    /// path, not just `e_homo`'s.** This module's own `Z_eff >= 1` proof
    /// (module doc) covers both probes with one argument, but they are not
    /// the same call: `zeff_outer` can probe a *different* subshell than
    /// the electron just placed, in the d/f-block fallback (`run.n <
    /// outer_n`, `probe_l = 0`) — `occ` already contains the just-placed
    /// electron in bucket `(run.n, run.l)` at that point, a different
    /// bucket than the one being probed, `(outer_n, 0)`. The proof's bound
    /// still holds there for a Madelung-order reason, not a screening one:
    /// `(n, 0)` always fills before `(n, l>0)` for a fixed `n` (smallest
    /// `n + l` at that `n`), so whenever any subshell at `n = outer_n` is
    /// occupied, `(outer_n, 0)` specifically already has at least one
    /// electron in it — which is exactly the placed electron the sum needs
    /// excluded to stay at or below `z - 1`. `> 0.0`, the bound an earlier
    /// version of `gap_and_e_homo_are_always_finite` checked, is weaker
    /// than what the module doc claims to prove; this checks the real one.
    /// **Swept across seeds, not just identity** — `sigma_near`'s value is
    /// what the proof bounds, not something identity happens to pick, so
    /// the check needs to see it vary. [`OrbitalConsts::draw`]'s one-sided
    /// direction keeps every drawn `sigma_near` in the same `(0,
    /// BASE_SIGMA_NEAR]` interval the proof assumes, at every rung.
    #[test]
    fn zeff_outer_is_at_least_one_through_the_real_call_path() {
        let runs = fill();
        for seed in [0u64, 1, 5, 21, 42] {
            let rung = crate::perturbation::Rung::draw(seed);
            let props = electronic_properties(
                &runs,
                OrbitalConsts::draw(seed, rung),
                GapConsts::draw(seed, rung),
            );
            for (i, p) in props.iter().enumerate() {
                assert!(
                    p.zeff_outer.is_finite() && p.zeff_outer >= 1.0,
                    "seed {seed} z={}: zeff_outer {} fell below 1.0",
                    i + 1,
                    p.zeff_outer
                );
            }
        }
    }

    /// **Every closed-`s^2` position with a same-shell `p` subshell to
    /// promote into promotes, and every closed-`p^6` position does not, for
    /// every z a real generated table can actually contain** — not the three
    /// hand-picked z values
    /// `promotion_budget_separates_alkaline_earth_from_noble_gas_analogues`
    /// checks. This is the guard the `.abs()` -> clamp fix (2026-08-10
    /// `/review-pr`) needed and did not have: with `.abs()`, z=56 (Ba
    /// analogue) and z=88 (Ra analogue) — both closed `s^2` — reported
    /// valence 0, and this property would have failed on both from the day
    /// `.abs()` shipped. Threshold-free within its range: reads the closures
    /// out of `fill()`'s own run structure, never a literal list of z.
    ///
    /// **Bounded at `z = 120`, not the fold's full reach to `z = 250`.**
    /// `n_elements` draws `60..=120` (Decision 9) — no `PeriodicTable` this
    /// crate ever generates contains an element past 120, so a property
    /// beyond it is a claim about electrons no universe can show. Narrowing
    /// the range keeps any calibration question about the fold's reach past
    /// `z = 120` — a real question, but a separate one — out of the
    /// property this test was written to guard.
    ///
    /// **`n = 1` is excluded from the `s^2` case, not an oversight.** `l < n`
    /// means `n = 1` has no `l = 1` at all — there is no same-shell `p`
    /// subshell for a 1s electron to promote into, structurally, not as a
    /// matter of energy. Real helium's valence is 0 for the same reason; see
    /// `valence_matches_real_ground_state_valences_through_period_3`'s own
    /// `z = 2` entry.
    ///
    #[expect(
        clippy::as_conversions,
        reason = "run.count is at most 18 (capacity_of(l)'s max), well within usize's exact \
                  integer range"
    )]
    #[test]
    fn closed_s2_positions_promote_and_closed_p6_positions_do_not() {
        let runs = fill();
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        let (mut s2_checked, mut p6_checked) = (0u32, 0u32);
        let mut z = 0usize;
        for run in &runs {
            z += run.count as usize;
            if z > 120 {
                break;
            }
            if run.count != run.capacity() {
                continue;
            }
            let want = match (run.n, run.l) {
                (1, 0) => continue,
                (_, 0) => 2u32,
                (_, 1) => 0u32,
                _ => continue,
            };
            assert_eq!(
                at(&props, z).valence,
                want,
                "z={z} ({},{}): a closed subshell of this type should have valence {want}",
                run.n,
                run.l
            );
            if run.l == 0 {
                s2_checked += 1;
            } else {
                p6_checked += 1;
            }
        }
        assert!(
            s2_checked >= 6 && p6_checked >= 6,
            "corpus filter selected too few positions to be a real property check: \
             s2={s2_checked} p6={p6_checked}"
        );
    }
}
