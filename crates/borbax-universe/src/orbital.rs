//! Screened-hydrogenic, `l`-dependent orbital filling (issue #26, Decisions
//! 3 and 8) — Task 26.1's replacement for `element.rs`'s V1 shell law.
//!
//! **The screening model, stated precisely because the plan states only its
//! principles.** For a candidate subshell `(n, l)` given the electrons
//! already placed, `Z_eff(n, l) = z - S(n, l)` where `z` is the electron
//! index about to be placed (the neutral atom's own atomic number under
//! construction) and:
//!
//! ```text
//! S(n, l) = sigma_inner * inner(n) + t(l) * same(n, l)
//! inner(n) = electrons already placed in any subshell with n' < n
//! same(n, l) = electrons already placed in a DIFFERENT subshell with n' == n
//!              (l' != l) — the candidate's own already-placed peers are
//!              excluded, see below
//! t(l) = min(l / 2.0, 1.0)   // 0 at l=0, 0.5 at l=1, saturated at l>=2
//! energy(n, l) = -(Z_eff(n, l) / n)^2   // plain n, not Slater's n*
//! ```
//!
//! **Excluding the candidate's own already-placed same-subshell peers is a
//! deliberate deviation from the most literal reading of Decision 3's text,
//! found necessary by direct numerical experiment, not stated in the
//! plan.** Including them (`same(n, l)` summing over every `n' == n`
//! regardless of `l'`) makes a partially-filled subshell's own energy rise
//! as it fills — self-screening among its own occupants — at a rate that,
//! for the real 3d/4p transition, overtakes 3d's initial energy advantage
//! after roughly 4 of its 10 electrons, before 3d finishes: 3d and 4p
//! trade off electron-by-electron instead of 3d filling cleanly.
//! Real atoms do not have this problem because of exchange-correlation
//! effects a screened-hydrogenic model does not include; excluding
//! same-subshell self-screening is the cheapest structural fix that keeps
//! the model's own account of *why* penetration matters (`t(l)`) intact.
//! Measured after the fix: the real subshell fill order is reproduced
//! exactly through 5p at the identity configuration, with no flicker in
//! that range.
//!
//! **`Z_eff >= 1` is provable, not merely observed, and the proof needs the
//! exclusion above.** `sigma_inner` is drawn through
//! [`crate::perturbation`] with base `OrbitalConsts::BASE_SIGMA_INNER` and a
//! one-sided-downward direction, so `sigma_inner` is always in
//! `(0, BASE_SIGMA_INNER] subset (0, 1]`; `t(l)` is bounded in `[0, 1]` by
//! construction. `inner(n) + same(n, l)` counts a subset of the electrons
//! already placed (every placed electron *except* the candidate's own
//! subshell's peers, which are excluded, and any electron in a subshell
//! with `n' > n`, which never contributes), so it is at most `z - 1`. Hence
//! `S(n, l) <= 1 * (z - 1) = z - 1`, and `Z_eff = z - S(n, l) >= z - (z - 1)
//! = 1`, for every reachable `z`, `n`, `l` and occupancy — not just the
//! cases this module happens to test.
//!
//! **A genuinely emergent, not hardcoded, consequence**: the model produces
//! its own transition-metal-region complexity near the 4d/5s boundary (a
//! single 5s electron placed, then 4d filling completely, then 5s's second
//! electron) — not identical to any specific real anomaly (real Pd, Nb, Mo,
//! Ru, Rh all break the simple pattern too, for reasons this model does not
//! capture), but the same qualitative phenomenon, falling out of the
//! formula rather than a per-element exception.
//!
//! **This module has no production caller yet — Steps 6-9 give it one.**
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
use borbax_units::canonical_cmp;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// One contiguous run of electrons [`fill`] placed into a single subshell
/// before moving to a different one.
///
/// **Not always a complete subshell.** Most entries are — the fold fills a
/// subshell to capacity before a lower-energy alternative exists — but the
/// 4d/5s region (see this module's own doc) genuinely interrupts a subshell
/// mid-fill, so the same `(n, l)` can appear as two separate, non-adjacent
/// entries. A reader wanting one element's total configuration at a given
/// electron count must sum every entry for a given `(n, l)` up to that
/// point, not assume the first entry for `(n, l)` is the only one.
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

/// The physical constants [`fill`] needs. Currently just the perturbed
/// inner-shell screening coefficient — see [`MigratedConstant::ScreeningInner`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OrbitalConsts {
    sigma_inner: f64,
}

impl OrbitalConsts {
    /// The unperturbed inner-shell screening coefficient. Calibrated by
    /// direct numerical search (not derived from Slater's empirical
    /// values, which this model does not use) to be the value at which the
    /// identity configuration's fold reproduces the real subshell fill
    /// order exactly through 5p. See this module's own doc for the fold
    /// this was calibrated against.
    pub(crate) const BASE_SIGMA_INNER: f64 = 0.65;

    /// Draw this universe's screening constants.
    ///
    /// **One-sided downward, per the plan's own text — `p` is sampled from
    /// `[-Direction::P_MAX, 0.0]`, never the symmetric range.** This is a
    /// draw-site decision, not a [`MigratedConstant`] classification tag:
    /// see [`MigratedConstant::ScreeningInner`]'s own doc for why it is not
    /// generalised into a fourth tag.
    #[must_use]
    pub(crate) fn draw(seed: u64, rung: Rung) -> Self {
        let mut stream = Stream::new(
            seed,
            Domain::Perturbation,
            MigratedConstant::ScreeningInner.index(),
        );
        let p = stream.next_f64_range(-Direction::P_MAX, 0.0);
        let direction = Direction::new(p, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("p is drawn within [-P_MAX, 0.0] by construction"));
        Self {
            sigma_inner: perturb(Self::BASE_SIGMA_INNER, rung, direction),
        }
    }

    #[cfg(test)]
    pub(crate) const fn at_identity() -> Self {
        Self {
            sigma_inner: Self::BASE_SIGMA_INNER,
        }
    }
}

/// The highest electron count [`fill`] computes to. Comfortably above the
/// naming grammar's 165-element symbol-space ceiling (`naming.rs`), so
/// every legal `n_elements` draw has a complete answer.
const ELECTRON_CEILING: u32 = 250;

/// The highest principal quantum number [`fill`] considers as a candidate.
/// `n = 9` reaches subshells well past [`ELECTRON_CEILING`]'s last electron,
/// so the fold never runs out of candidates before the ceiling.
const MAX_N: u8 = 9;

/// Every `(n, l)` subshell up to [`MAX_N`], in ascending `(n, l)` order —
/// the canonical order the fold iterates in, and the tie-break order two
/// candidates at bit-identical energy resolve through.
fn candidates() -> Vec<(u8, u8)> {
    let mut out = Vec::new();
    for n in 1..=MAX_N {
        for l in 0..n {
            out.push((n, l));
        }
    }
    out
}

/// The `(energy, Z_eff)` pair for subshell `(n, l)` given electrons already
/// placed, at electron count `z` — the fold's own energy formula (this
/// module's doc), extracted so [`best_unfilled`] (candidate selection) and
/// the derived electronic-structure properties below (Steps 6-9) share one
/// formula rather than two copies that could drift.
///
/// **Deliberately takes `occ` as-is, with no capacity check.** Whether
/// `(n, l)` still has room is the caller's question — [`best_unfilled`]
/// skips full subshells before calling this; the Steps 6-9 properties below
/// call it *because* a subshell is full, to price the next candidate.
fn subshell_energy(
    n: u8,
    l: u8,
    z: u32,
    occ: &BTreeMap<(u8, u8), u32>,
    consts: OrbitalConsts,
) -> (f64, f64) {
    let inner: u32 = occ
        .iter()
        .filter(|((n2, _), _)| *n2 < n)
        .map(|(_, c)| c)
        .sum();
    let same: u32 = occ
        .iter()
        .filter(|((n2, l2), _)| *n2 == n && *l2 != l)
        .map(|(_, c)| c)
        .sum();
    let s = consts.sigma_inner * f64::from(inner) + t_of_l(l) * f64::from(same);
    let z_eff = f64::from(z) - s;
    let ratio = z_eff / f64::from(n);
    (-(ratio * ratio), z_eff)
}

/// The lowest-energy unfilled candidate, given the electrons already placed.
///
/// **Iterates `order`, not a canonical iteration of `occ`, so the caller
/// controls the scan order — this is what lets a test probe order-
/// invariance directly, by calling this twice with two different orders
/// over the same `occ` and asserting the same winner.** The winner cannot
/// depend on `order` regardless: every candidate is still visited exactly
/// once, and the comparison is a strict total order
/// (`canonical_cmp` on energy, then `(n, l)` ascending as the tie-break),
/// so a linear scan finds the same minimum however it is traversed.
///
/// `occ`'s own iteration (to sum `inner`/`same`, inside [`subshell_energy`])
/// uses `BTreeMap`'s sorted order, never `HashMap` — CLAUDE.md's §13.1 —
/// though this is moot for the *sums themselves*, which are order-invariant
/// integer addition regardless of traversal order.
fn best_unfilled(
    order: &[(u8, u8)],
    occ: &BTreeMap<(u8, u8), u32>,
    z: u32,
    consts: OrbitalConsts,
) -> (u8, u8) {
    let mut best: Option<(f64, u8, u8)> = None;
    for &(n, l) in order {
        let count = occ.get(&(n, l)).copied().unwrap_or(0);
        if count >= capacity_of(l) {
            continue;
        }
        let (e, _z_eff) = subshell_energy(n, l, z, occ, consts);
        let better = match best {
            None => true,
            Some((best_e, best_n, best_l)) => match canonical_cmp(e, best_e) {
                Ordering::Less => true,
                Ordering::Equal => (n, l) < (best_n, best_l),
                Ordering::Greater => false,
            },
        };
        if better {
            best = Some((e, n, l));
        }
    }
    best.map_or_else(
        || unreachable!("candidates() covers {ELECTRON_CEILING} electrons' worth of subshells at MAX_N = {MAX_N}, so a candidate is always available"),
        |(_, n, l)| (n, l),
    )
}

/// Fill subshells electron by electron, in aufbau order, up to
/// [`ELECTRON_CEILING`].
///
/// Returns one [`Occupancy`] per contiguous run — see [`Occupancy`]'s own
/// doc for why that is not always one entry per subshell.
#[must_use]
pub(crate) fn fill(consts: OrbitalConsts) -> Vec<Occupancy> {
    let order = candidates();
    let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
    let mut out: Vec<Occupancy> = Vec::new();
    for z in 1..=ELECTRON_CEILING {
        let (n, l) = best_unfilled(&order, &occ, z, consts);
        *occ.entry((n, l)).or_insert(0) += 1;
        match out.last_mut() {
            Some(last) if last.n == n && last.l == l => last.count += 1,
            _ => out.push(Occupancy { n, l, count: 1 }),
        }
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
    /// Calibrated by direct numerical search (same method as
    /// [`OrbitalConsts::BASE_SIGMA_INNER`]) to sit strictly between real
    /// alkaline-earth-like HOMO-LUMO gaps (~1.1-1.6 at the identity
    /// configuration, e.g. beryllium- and magnesium-analogues) and real
    /// noble-gas-like ones (~3.9 and up, e.g. helium- and neon-analogues) —
    /// a >25%-margin choice on both sides, not a boundary value.
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
    /// `Z_eff` of `(outer_n, l = 0)` at this z — feeds radius.
    pub(crate) zeff_outer: f64,
    /// The HOMO-LUMO gap. Zero whenever `fill()`'s own run structure shows
    /// this electron's subshell is not yet at capacity; nonzero only when
    /// this electron exactly completes its subshell **and** a different
    /// subshell follows. Read from the fold's actual output, not
    /// re-derived — "gap = 0 mid-run" is not independently provable from
    /// screening alone (this module's own doc records the counterexample:
    /// the fold can and does interrupt a partially-filled subshell, e.g.
    /// 5s at z=37→38, purely because a smaller-n competitor's energy has
    /// higher curvature in `Z_eff` under the identical +1 shift).
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
    /// frontier subshell** — a closed subshell always contributes 0
    /// unpaired electrons, so this differs from a frontier-only count only
    /// during an interruption (this module's own doc: the 5s/4d region has
    /// *two* simultaneously-incomplete subshells, both contributing).
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
            let (_, zeff_outer) = subshell_energy(outer_n, 0, z, &occ, orbital_consts);
            let (e_homo, _) = subshell_energy(run.n, run.l, z, &occ, orbital_consts);

            let count_now = occ
                .get(&(run.n, run.l))
                .copied()
                .unwrap_or_else(|| unreachable!("just inserted above"));
            let subshell_full = count_now == capacity_of(run.l);
            let gap = if subshell_full {
                runs.get(run_index + 1).map_or(0.0, |next| {
                    let (e_lumo, _) = subshell_energy(next.n, next.l, z, &occ, orbital_consts);
                    e_lumo - e_homo
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
        ELECTRON_CEILING, ElectronicProperties, GapConsts, Occupancy, OrbitalConsts, best_unfilled,
        candidates, capacity_of, electronic_properties, fill,
    };
    use borbax_rng::{Domain, Stream};
    use std::collections::BTreeMap;

    /// The real subshell fill order through 5p — the output of running this
    /// module's formula at the identity configuration, pinned as a
    /// computed constant (CLAUDE.md, §5-permitted under Decision 2), not a
    /// lookup table.
    const IDENTITY_FILL_ORDER: [(u8, u8); 11] = [
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

    #[test]
    fn fill_order_is_invariant_to_candidate_shuffle() {
        let consts = OrbitalConsts::at_identity();
        let canonical = candidates();
        let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        // Build a representative mid-fill occupancy (through argon) so the
        // shuffle test exercises a nontrivial, multi-subshell state rather
        // than the trivial all-empty one.
        for z in 1..=18 {
            let (n, l) = best_unfilled(&canonical, &occ, z, consts);
            *occ.entry((n, l)).or_insert(0) += 1;
        }
        let canonical_pick = best_unfilled(&canonical, &occ, 19, consts);

        let mut stream = Stream::new(0x5EED, Domain::Reaction, 0);
        let mut shuffled = canonical.clone();
        for i in (1..shuffled.len()).rev() {
            let bound = u64::try_from(i + 1).unwrap_or_else(|_| {
                unreachable!("candidates().len() is a small constant, well within u64")
            });
            let j = usize::try_from(stream.next_range(bound))
                .unwrap_or_else(|_| unreachable!("next_range(bound) is < bound, which fits usize"));
            shuffled.swap(i, j);
        }
        assert_ne!(
            shuffled, canonical,
            "the shuffle should actually reorder the list"
        );
        let shuffled_pick = best_unfilled(&shuffled, &occ, 19, consts);

        assert_eq!(
            canonical_pick, shuffled_pick,
            "the winning subshell must not depend on scan order"
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

        let order = fill(OrbitalConsts::at_identity());
        // Through 4p (the first 8 entries of IDENTITY_FILL_ORDER), the fold
        // is uninterrupted (see this module's own doc on the 4d/5s region),
        // so fill()'s first 8 entries line up 1:1 with IDENTITY_FILL_ORDER.
        let mut running = 0u32;
        let mut checked = 0usize;
        for (occ, want) in order.iter().take(8).zip(expected_running.iter()) {
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
            checked, 8,
            "should have checked all 8 uninterrupted entries"
        );
    }

    #[test]
    fn the_identity_fill_order_matches_the_real_periodic_table_through_5p() {
        let order = fill(OrbitalConsts::at_identity());
        let touched: Vec<(u8, u8)> = order.iter().map(|o| (o.n, o.l)).collect();
        let mut first_seen: Vec<(u8, u8)> = Vec::new();
        for &(n, l) in &touched {
            if !first_seen.contains(&(n, l)) {
                first_seen.push((n, l));
            }
        }
        let first_seen_prefix: Vec<(u8, u8)> = first_seen
            .iter()
            .copied()
            .take(IDENTITY_FILL_ORDER.len())
            .collect();
        assert_eq!(
            first_seen_prefix, IDENTITY_FILL_ORDER,
            "the identity configuration's fill order should match the real periodic table \
             through 5p"
        );
    }

    /// Blocker 7's actual regression test: an `l`-blind screening
    /// implementation (screening that does not depend on the candidate's
    /// own `l`) fails to reproduce the real fill order — confirming `t(l)`
    /// is load-bearing, not decorative.
    #[test]
    fn an_l_blind_screening_implementation_fails_the_identity_order() {
        fn best_unfilled_l_blind(
            order: &[(u8, u8)],
            occ: &BTreeMap<(u8, u8), u32>,
            z: u32,
            sigma_inner: f64,
        ) -> (u8, u8) {
            let mut best: Option<(f64, u8, u8)> = None;
            for &(n, l) in order {
                let count = occ.get(&(n, l)).copied().unwrap_or(0);
                if count >= capacity_of(l) {
                    continue;
                }
                let inner: u32 = occ
                    .iter()
                    .filter(|((n2, _), _)| *n2 < n)
                    .map(|(_, c)| c)
                    .sum();
                let same: u32 = occ
                    .iter()
                    .filter(|((n2, l2), _)| *n2 == n && *l2 != l)
                    .map(|(_, c)| c)
                    .sum();
                // The planted defect: a constant 0.5 in place of t(l).
                let s = sigma_inner * f64::from(inner) + 0.5 * f64::from(same);
                let z_eff = f64::from(z) - s;
                let ratio = z_eff / f64::from(n);
                let e = -(ratio * ratio);
                let better = match best {
                    None => true,
                    Some((be, bn, bl)) => {
                        use borbax_units::canonical_cmp;
                        use std::cmp::Ordering;
                        match canonical_cmp(e, be) {
                            Ordering::Less => true,
                            Ordering::Equal => (n, l) < (bn, bl),
                            Ordering::Greater => false,
                        }
                    }
                };
                if better {
                    best = Some((e, n, l));
                }
            }
            best.map_or_else(|| unreachable!(), |(_, n, l)| (n, l))
        }

        let order = candidates();
        let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        let mut touched: Vec<(u8, u8)> = Vec::new();
        for z in 1..=54 {
            let (n, l) = best_unfilled_l_blind(&order, &occ, z, OrbitalConsts::BASE_SIGMA_INNER);
            *occ.entry((n, l)).or_insert(0) += 1;
            if !touched.contains(&(n, l)) {
                touched.push((n, l));
            }
        }
        let n = IDENTITY_FILL_ORDER.len().min(touched.len());
        let touched_prefix: Vec<(u8, u8)> = touched.iter().copied().take(n).collect();
        let target_prefix: Vec<(u8, u8)> = IDENTITY_FILL_ORDER.iter().copied().take(n).collect();
        assert_ne!(
            touched_prefix, target_prefix,
            "an l-blind screening function should NOT reproduce the real fill order — if it \
             does, t(l) is not load-bearing and this test is not discriminating anything"
        );
    }

    #[test]
    fn z_eff_is_at_least_one_with_equality_at_z_one() {
        let consts = OrbitalConsts::at_identity();
        let order = candidates();
        let mut occ: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        for z in 1..=ELECTRON_CEILING {
            // Recompute the winning candidate's own Z_eff the same way
            // best_unfilled does internally, to assert on it directly.
            let (n, l) = best_unfilled(&order, &occ, z, consts);
            let inner: u32 = occ
                .iter()
                .filter(|((n2, _), _)| *n2 < n)
                .map(|(_, c)| c)
                .sum();
            let same: u32 = occ
                .iter()
                .filter(|((n2, l2), _)| *n2 == n && *l2 != l)
                .map(|(_, c)| c)
                .sum();
            let s = consts.sigma_inner * f64::from(inner) + super::t_of_l(l) * f64::from(same);
            let z_eff = f64::from(z) - s;
            assert!(
                z_eff >= 1.0,
                "z={z} ({n},{l}): Z_eff {z_eff} fell below 1.0"
            );
            if z == 1 {
                assert!(
                    (z_eff - 1.0).abs() < f64::EPSILON,
                    "Z_eff at Z=1 should be exactly 1.0"
                );
            }
            *occ.entry((n, l)).or_insert(0) += 1;
        }
    }

    #[test]
    fn the_screening_accumulation_is_order_invariant_and_exact() {
        // Two different insertion orders reaching the same final occupancy:
        // fill argon (Z=1..=18) via the real fold (order A), and separately
        // insert the same final counts in reverse-(n,l) order (order B).
        let consts = OrbitalConsts::at_identity();
        let order = candidates();
        let mut occ_a: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        for z in 1..=18 {
            let (n, l) = best_unfilled(&order, &occ_a, z, consts);
            *occ_a.entry((n, l)).or_insert(0) += 1;
        }

        let mut occ_b: BTreeMap<(u8, u8), u32> = BTreeMap::new();
        for (&key, &count) in occ_a.iter().rev() {
            occ_b.insert(key, count);
        }
        assert_eq!(
            occ_a, occ_b,
            "both occupancy maps should hold the same final counts"
        );

        let pick_a = best_unfilled(&order, &occ_a, 19, consts);
        let pick_b = best_unfilled(&order, &occ_b, 19, consts);
        assert_eq!(pick_a, pick_b);

        // And the Z_eff bits for the winning candidate agree exactly.
        let z_eff_of = |occ: &BTreeMap<(u8, u8), u32>, n: u8, l: u8| -> f64 {
            let inner: u32 = occ
                .iter()
                .filter(|((n2, _), _)| *n2 < n)
                .map(|(_, c)| c)
                .sum();
            let same: u32 = occ
                .iter()
                .filter(|((n2, l2), _)| *n2 == n && *l2 != l)
                .map(|(_, c)| c)
                .sum();
            let s = consts.sigma_inner * f64::from(inner) + super::t_of_l(l) * f64::from(same);
            19.0 - s
        };
        let (n, l) = pick_a;
        assert_eq!(
            z_eff_of(&occ_a, n, l).to_bits(),
            z_eff_of(&occ_b, n, l).to_bits(),
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
        let total: u32 = fill(OrbitalConsts::at_identity())
            .iter()
            .map(|o| o.count)
            .sum();
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
        let runs = fill(OrbitalConsts::at_identity());
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
        let runs = fill(OrbitalConsts::at_identity());
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
        let runs = fill(OrbitalConsts::at_identity());
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
        let runs = fill(OrbitalConsts::at_identity());
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        let mut prev = at(&props, 1).outer_n;
        for z in 2..=60 {
            let outer_n = at(&props, z).outer_n;
            assert!(outer_n >= prev, "z={z}: outer_n must never decrease");
            prev = outer_n;
        }
        // Real period starts at the identity configuration, through period 5.
        for &z in &[3usize, 11, 19, 37] {
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
        let runs = fill(OrbitalConsts::at_identity());
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        // Period 2 spans z=3..=10: outer_n stays 2 throughout, and
        // zeff_outer should rise strictly.
        for z in 4..=10 {
            assert_eq!(at(&props, z).outer_n, 2, "z={z}: should stay in period 2");
            assert!(
                at(&props, z).zeff_outer > at(&props, z - 1).zeff_outer,
                "z={z}: zeff_outer should rise within a period"
            );
        }
        // z=11 (Na) starts period 3: outer_n jumps, so zeff_outer -- and
        // therefore radius -- should NOT continue period 2's own trend.
        assert_eq!(at(&props, 11).outer_n, 3, "z=11 should be in a new period");
    }

    #[test]
    fn gap_and_e_homo_are_always_finite() {
        let runs = fill(OrbitalConsts::at_identity());
        let props = electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        for (i, p) in props.iter().enumerate() {
            assert!(p.gap.is_finite(), "z={}: gap should be finite", i + 1);
            assert!(p.e_homo.is_finite(), "z={}: e_homo should be finite", i + 1);
            assert!(
                p.zeff_outer.is_finite() && p.zeff_outer > 0.0,
                "z={}: zeff_outer should be finite and positive",
                i + 1
            );
            assert!(p.gap >= 0.0, "z={}: gap should never be negative", i + 1);
        }
    }
}
