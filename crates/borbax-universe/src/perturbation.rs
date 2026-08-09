//! The shared perturbation mechanism (issue #26, Decision 10) — one integer
//! rung per universe, and a bounded per-constant direction, composed as
//! `base * (1.0 + strength * p)`.
//!
//! **Scope note, and why this module has no consumer yet.** `Rung` and
//! `Direction` are complete, independently tested types, and `perturb`
//! is the mechanism's whole formula. What this module deliberately does
//! *not* do is touch [`crate::Universe`] or [`crate::element::generate_elements`]:
//! every constant the plan actually migrates through this mechanism
//! (`w_shape`, `w_charge`, `contact_defect`, `base_mass`, `bonds.rs`'s rate
//! constants) is tagged "V2+ only" or "V3 only" in the plan's own migration
//! table, and `PhysicsVersion::V2` does not exist in this crate yet. Wiring
//! a `rung` field onto `Universe` today, ahead of a real consumer, was
//! measured (by the determinism-auditor, before this landed) to be safe
//! against `PhysicsVersion::V1`'s two pinned goldens *only* by discarding
//! the field — but Decision 9's `n_elements` gate, which the real-name guard
//! (P6, `naming.rs`) needs `Rung` for, has no `PhysicsVersion` qualifier in
//! the plan's own text, and wiring it into V1's actual generation path would
//! let *any* V1 universe — which has no real-physics correspondence at all,
//! being "shape and surface character" per `PhysicsVersion::V1`'s own doc —
//! claim `IdentityWitness` and dispense real element names. That is a G2
//! breach, not a golden-stability question, and it is the same shape of
//! problem P1's `CLOSURE_PREDICATE_WATCHED_FNS` redesign and P4's explicit
//! "lands after Task 26.1 Step 8" both exist to avoid: a mechanism whose
//! real subject does not exist yet. Task 26.1 is where `Universe` gains a
//! `rung` field and Decision 9's gate gains a consumer, together.
//!
//! **`Rung`'s own visibility, corrected in `/review-pr` — an earlier version
//! of this paragraph got both of its claims wrong.** It said `Rung` "has
//! to be" `pub` because `crate::naming::IdentityWitness::new` takes one as a
//! parameter — false: `new` is `pub(crate)`, so only in-crate callers ever
//! need to name the type, and this module (and the type) are `pub(crate)`
//! now. It also said "the trap is a future `Universe.rung` field, not
//! `Rung` itself" — also false, and demonstrated so in review:
//! `Rung::draw` is a pure function of the seed, so `Rung::draw(seed).is_identity()`
//! already recovered the identity predicate from `borbax-molecule` with no
//! field and no witness, while `Rung` was still `pub`. Narrowing this module
//! closes that route today, not merely at some future field.
//!
//! **This is narrow-API hygiene (CLAUDE.md's idiom tier), not a G2 fix, and
//! the distinction is worth keeping straight.** No fiction guarantee is at
//! stake either way: `Universe::seed` is already `pub`, so the identity
//! predicate is already reconstructible today via `borbax_rng::Domain::PerturbationRung`
//! regardless of this module's own visibility, and G2's 2026-08-07 revision
//! does not treat the identity configuration as a secret — narrowing here is
//! free (nothing outside this crate names `Rung`, `Direction` or `perturb`),
//! not load-bearing. Filing it as a G2 fix would imply the identity
//! configuration needs protecting, which contradicts the ruling that
//! unblocked issue #26 in the first place.
//!
//! When Task 26.1 gives `Universe` a `rung` field, that field's own
//! visibility is the decision that actually matters — it must be held no
//! more visible than the code that needs to read it requires, not given
//! `pub` by the struct's house-style convention.
//!
//! **Every item below is `dead_code` in a `not(test)` build, module-wide,
//! for the reason this whole module doc has already given at length: no
//! in-crate caller exists until Task 26.1.** One blanket `#[expect]` rather
//! than thirteen individually-reasoned ones, because the reason is one fact
//! ("this module's real consumer does not exist yet"), not thirteen —
//! `cfg_attr`-gated to `not(test)` for the same reason every other such
//! attribute in this crate is: this module's own tests are callers, and an
//! unconditional `#[expect(dead_code)]` would fire `unfulfilled_lint_expectations`
//! on the test target.
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "Rung, Direction and perturb have no in-crate caller until Task 26.1 gives \
                  Universe a rung field and wires the first migrated constant through this \
                  mechanism — see this module's own doc for why landing the mechanism ahead \
                  of its consumer is deliberate"
    )
)]

use borbax_rng::{Domain, Stream};

/// How far a universe's constants sit from the identity configuration —
/// `0..=MAX`, drawn once per universe.
///
/// **`MAX` is provisional.** Routed requirement 8 (issue #26's plan) calls
/// for deriving it from a target `P(identity) = 1/(MAX + 1)`, and that
/// target has not been chosen. `8` is a placeholder giving `P(identity) =
/// 1/9 ≈ 11.1%` — large enough to reach in a small seed sweep, small enough
/// that "identity" stays a minority of universes. Task 26.1 revisits this
/// value when the real target is picked; nothing here depends on the exact
/// number beyond what `the_binomial_band_pins_p_identity_as_a_number` (this
/// module's own `#[cfg(test)]` tests, not linkable from public docs) checks
/// against whatever `MAX` currently is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Rung(u8);

impl Rung {
    /// The top of the rung's legal range, `m` in the plan's notation.
    pub(crate) const MAX: u8 = 8;

    /// `Rung(0)` — the identity configuration.
    pub(crate) const IDENTITY: Self = Self(0);

    /// A rung, or `None` outside `0..=MAX`.
    ///
    /// **Fallible, not a bare `u8` field, for the reason `BondOrder::new` and
    /// `Direction::new` both are: an unbounded rung reopens the excursion
    /// hazard `Direction`'s own bound exists to close.** At `rung > MAX` the
    /// perturbation factor's infimum is `+0.0` exactly once `rung = 2*MAX`,
    /// and past `4*MAX` a quarter of `p`'s range gives a **negative**
    /// factor — see [`perturb`]'s doc for the proof this constructor makes
    /// unreachable.
    #[must_use]
    pub(crate) const fn new(r: u8) -> Option<Self> {
        if r <= Self::MAX { Some(Self(r)) } else { None }
    }

    /// This universe's rung, from its own stream — never a shared index.
    ///
    /// **`next_range(MAX + 1)`, not `next_range(2*MAX + 1)`.** The latter
    /// gives `2*MAX + 1` outcomes and halves `P(identity)` silently; the
    /// former is the spelling that actually delivers `P(identity) =
    /// 1/(MAX + 1)`, pinned by `the_binomial_band_pins_p_identity_as_a_number`
    /// (this module's own tests) as a number, not an order-of-magnitude
    /// check that both readings would satisfy.
    ///
    /// Draws from `Stream::new(seed, Domain::PerturbationRung, 0)` — its own
    /// `Domain`, never an index shared with anything else. See
    /// `Domain::PerturbationRung`'s own doc for why a shared index is a live
    /// hazard under a discrete `p` ladder, not merely a style preference.
    ///
    /// **Routes through [`Self::new`], found bypassing it in `/review-pr`.**
    /// An earlier version built `Self(r as u8)` directly — the only
    /// production producer of a `Rung` skipping the one place `r <= MAX` is
    /// enforced. Under the exact defect named above, that bypass let an
    /// out-of-range `Rung` reach [`perturb`]'s annihilation proof, which
    /// assumes every `Rung` was checked. `new`'s own fallible check makes
    /// this the same shape of fix as `Direction::new`'s `bound <= P_MAX`
    /// clause: the invariant enforced at the one place a caller could skip
    /// it, not merely at every site that happens to construct one today.
    #[must_use]
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        reason = "next_range(u64::from(Self::MAX) + 1) returns a value in 0..=MAX, and MAX is a \
                  u8 literal, so the result fits u8 by construction — proven by the domain of \
                  next_range's input, not merely assumed"
    )]
    pub(crate) fn draw(seed: u64) -> Self {
        let mut stream = Stream::new(seed, Domain::PerturbationRung, 0);
        let r = stream.next_range(u64::from(Self::MAX) + 1);
        Self::new(r as u8).unwrap_or_else(|| unreachable!("next_range(MAX + 1) cannot exceed MAX"))
    }

    /// Every rung except the identity, `1..=MAX` — for Decision 9's
    /// exhaustive test (`IdentityWitness::new(r).is_none()` for every `r`
    /// here), so the corpus is derived from `MAX` rather than a hardcoded
    /// `1..=RUNGS` literal that could drift from it.
    pub(crate) fn all_off_identity() -> impl Iterator<Item = Self> {
        (1..=Self::MAX).map(Self)
    }

    /// `true` for `Rung(0)` — the identity configuration.
    #[must_use]
    pub(crate) const fn is_identity(self) -> bool {
        self.0 == 0
    }

    /// This rung's magnitude, `rung / MAX`, in `[0.0, 1.0]` — derived, never
    /// stored. "Store the integer rung, not the derived `f64` strength"
    /// (Decision 10): the integer is what every test and every digest
    /// should key on, so a strength recomputed here can never disagree with
    /// one stored somewhere else.
    #[must_use]
    pub(crate) fn strength(self) -> f64 {
        f64::from(self.0) / f64::from(Self::MAX)
    }
}

/// A bounded, finite perturbation direction — `p` in the plan's notation.
///
/// **The bound is a constructor parameter, not a shared constant.** A
/// single global `P_MAX` cannot express `w_shape`/`w_charge`'s tighter
/// `1/3` bound (issue #26's plan, P7), so `new` takes `bound` explicitly
/// rather than reading [`Self::P_MAX`] internally — that constant is the
/// *default* most constants use, not the only legal value.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct Direction(f64);

impl Direction {
    /// The default excursion bound most migrated constants use. `w_shape`
    /// and `w_charge` are tighter, at `1/3` — see issue #26's plan, P7, for
    /// the numerical-soundness argument (a monotone Möbius transform of
    /// their ratio, degenerating at its extremes) that bound rests on.
    pub(crate) const P_MAX: f64 = 0.5;

    /// A direction, or `None` if `v` is non-finite, `bound` exceeds
    /// [`Self::P_MAX`], or `v` exceeds `bound` in magnitude.
    ///
    /// **Finite-and-bounded by construction, not by a `debug_assert!`** —
    /// which is absent from the `--release` profile `goldens --emit` runs
    /// in, and this workspace has already shipped that exact debug/release
    /// asymmetry twice (`borbax-units`, Task 2). `borbax_rng`'s own
    /// `next_f64_range` documents swallowing a non-finite bound (`lo =
    /// -inf` with finite `hi` returns `-inf` on every draw), so an infinite
    /// `p` is producible by a typo'd range and nothing upstream of this
    /// constructor rejects it.
    ///
    /// **`bound` is itself checked against `P_MAX`, found missing in
    /// review.** [`perturb`]'s "annihilation is unreachable" proof holds
    /// "for any `Direction` constructed with `bound <= 0.5`" — a caller
    /// obligation the earlier version of this constructor never enforced,
    /// so `Direction::new(-2.0, 5.0)` would have succeeded and, fed through
    /// `perturb` at `Rung::MAX`, produced a negative multiplier the proof
    /// claims impossible. Checking `bound` here makes it a structural
    /// guarantee instead, at the one call site rather than at every future
    /// one Task 26.1 adds.
    ///
    /// **`-0.0` is canonicalised to `+0.0`, per Decision 10's own text**
    /// ("canonicalise `-0.0` in `p` before any digest hashes the
    /// perturbation vector"). Adding `+0.0` rather than comparing against
    /// zero: IEEE-754 addition leaves every finite non-zero value unchanged
    /// and turns negative zero into positive zero under the default
    /// rounding mode, so this moves exactly one bit pattern and needs no
    /// `float_cmp` suppression.
    ///
    /// **Not observable through [`perturb`] at `Rung::IDENTITY`, which is
    /// why this constructor is tested directly rather than through it —
    /// found claimed the other way round in `/review-pr`.** At
    /// `strength() == 0.0`, `0.0 * -0.0` and `0.0 * 0.0` are both `0.0`,
    /// `1.0 + 0.0` is exactly `1.0` regardless of which zero sign fed it,
    /// and `base * 1.0` preserves `base` bit for bit — so the *multiplier*
    /// is `1.0` either way and `perturb`'s result cannot distinguish the two
    /// input signs at that rung. An earlier version of this module claimed
    /// the opposite, that routing `-0.0` through `perturb` "verified" the
    /// canonicalisation; it verified only that identity-rung perturbation
    /// is a no-op, which was never in question. See `a_negative_zero_direction_is_canonicalised`
    /// (this module's own tests) for the assertion that actually observes
    /// this constructor's output.
    #[must_use]
    pub(crate) fn new(v: f64, bound: f64) -> Option<Self> {
        if v.is_finite() && bound <= Self::P_MAX && v.abs() <= bound {
            Some(Self(v + 0.0))
        } else {
            None
        }
    }

    /// The raw value.
    #[must_use]
    pub(crate) const fn get(self) -> f64 {
        self.0
    }
}

/// Apply Decision 10's mechanism: `base * (1.0 + strength * p)`, never the
/// additive form.
///
/// **Proof that annihilation (`base * 0.0`) is unreachable, restated from
/// the plan rather than merely cited.** `1.0 + strength * p` is exactly
/// `0.0` iff `strength * p == -1.0`. `rung.strength()` is in `[0.0, 1.0]` by
/// [`Rung::new`]'s domain, and every [`Direction`] is constructed with
/// `bound <= 0.5` — enforced by [`Direction::new`] itself, not merely a
/// caller convention, since a review found the earlier version of that
/// constructor left `bound` unchecked — so `strength * p` is in `[-0.5,
/// 0.5]` and `1.0 + strength * p` is in `[0.5, 1.5]`, which never reaches
/// `0.0` — proven at the worst-case corner by
/// `the_multiplier_never_reaches_zero_at_the_worst_case_corner` (this
/// module's own tests), not merely asserted.
///
/// **The multiplier staying in `[0.5, 1.5]` does not by itself mean `base *
/// multiplier` cannot underflow to `0.0` — found stated too broadly in
/// `/review-pr`.** For any `base` above the subnormal range this is moot,
/// but at `base` in the smallest few subnormals (`f64::from_bits(1)` being
/// the extreme case) `base * 0.5` itself underflows to `0.0` — not a live
/// hazard for anything this mechanism migrates today (every constant listed
/// for Task 26.1's migration is `O(1)` or larger), but worth stating
/// precisely rather than as an unqualified "unreachable".
///
/// **The expression shape is pinned, not merely preferred.** `base * (1.0 +
/// strength * p)` and the algebraically-equal `base + base * strength * p`
/// are different `f64` values — floating-point multiplication does not
/// distribute over addition — so redistributing this after Task 26.1 wires
/// the first migrated constant moves every golden downstream of it, the
/// same class of silent physics change §13.1 exists to catch. Commuting
/// either multiply (`p * strength`, or `(1.0 + ..) * base`) is bit-identical
/// and therefore safe; distributing is not.
#[must_use]
pub(crate) fn perturb(base: f64, rung: Rung, p: Direction) -> f64 {
    base * (1.0 + rung.strength() * p.get())
}

#[cfg(test)]
mod tests {
    use super::{Direction, Rung, perturb};

    #[test]
    fn a_rung_past_max_is_rejected() {
        assert!(Rung::new(Rung::MAX).is_some());
        assert!(Rung::new(Rung::MAX + 1).is_none());
        assert!(Rung::new(u8::MAX).is_none());
    }

    #[test]
    fn all_off_identity_is_exactly_one_through_max() {
        // `r.0`: `mod tests` is a child of `perturbation`, so the private
        // field is visible here — a round-trip through `strength()` would
        // test the wrong thing (that `strength` inverts cleanly), not that
        // `all_off_identity` yields the right integers.
        let got: Vec<u8> = Rung::all_off_identity().map(|r| r.0).collect();
        let want: Vec<u8> = (1..=Rung::MAX).collect();
        assert_eq!(got, want);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "CLAUDE.md: tests may assert exactly. 0/MAX and MAX/MAX are exact under \
                  IEEE-754 division (numerator equals denominator or is zero), so an epsilon \
                  would defeat the point of asserting on it"
    )]
    fn identity_has_zero_strength_and_every_other_rung_is_positive() {
        assert_eq!(Rung::IDENTITY.strength(), 0.0);
        assert!(Rung::IDENTITY.is_identity());
        for r in Rung::all_off_identity() {
            assert!(r.strength() > 0.0, "rung {r:?} has non-positive strength");
            assert!(!r.is_identity());
        }
        let max_rung = Rung::new(Rung::MAX)
            .unwrap_or_else(|| unreachable!("Rung::MAX is in 0..=Rung::MAX by construction"));
        assert_eq!(max_rung.strength(), 1.0);
    }

    #[test]
    fn a_direction_rejects_non_finite_and_out_of_bound_values() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                Direction::new(v, Direction::P_MAX).is_none(),
                "{v} should be rejected"
            );
        }
        assert!(Direction::new(0.6, Direction::P_MAX).is_none());
        assert!(Direction::new(-0.6, Direction::P_MAX).is_none());
        // The boundary itself is legal: `<=`, not `<`.
        assert!(Direction::new(Direction::P_MAX, Direction::P_MAX).is_some());
        assert!(Direction::new(-Direction::P_MAX, Direction::P_MAX).is_some());
        assert!(Direction::new(0.0, Direction::P_MAX).is_some());
    }

    /// **The `bound <= P_MAX` clause itself, found untested in `/review-pr`.**
    /// Deleting it left all 96 `borbax-universe` tests green — the clause is
    /// what makes [`perturb`]'s "annihilation is unreachable" proof
    /// structural rather than a caller convention. Without it,
    /// `Direction::new(-2.0, 5.0)` succeeds, and at `Rung::MAX` `perturb`'s
    /// multiplier is exactly `-1.0` — the sign flip the proof claims
    /// impossible for every migrated constant, once Task 26.1 wires a real
    /// call site.
    #[test]
    fn a_bound_wider_than_p_max_is_rejected() {
        assert!(Direction::new(-2.0, 5.0).is_none());
        assert!(Direction::new(0.0, f64::INFINITY).is_none());
        assert!(Direction::new(0.0, f64::NAN).is_none());
        // The property, not just the arm: every constructible pair keeps
        // perturb's multiplier inside [0.5, 1.5].
        for r in (0..=Rung::MAX).filter_map(Rung::new) {
            for raw in [-0.5, -0.25, 0.0, 0.5, -2.0, 3.0, 5.0] {
                for bound in [0.5, 1.0 / 3.0, 0.0, 1.0, 5.0] {
                    if let Some(d) = Direction::new(raw, bound) {
                        let m = perturb(1.0, r, d);
                        assert!((0.5..=1.5).contains(&m), "multiplier {m} escaped");
                    }
                }
            }
        }
    }

    /// Decision 10's own required test: at `rung == 0`, every base value is
    /// bit-identical, over the full finite range of `p` at both signs —
    /// **compared by `to_bits()`, not by value**, because an earlier
    /// version of this test used `assert_eq!` on the `f64`s directly, and
    /// `assert_eq!(-0.0, 0.0)` passes: the one input in the `base` list
    /// chosen to exercise sign-of-zero preservation was exactly the input
    /// the comparison could not see a mismatch in.
    ///
    /// **`-0.0` is in the `p` list for parity with Decision 10's text ("at
    /// both zero signs of `p`"), not because this test can see
    /// `Direction::new`'s canonicalisation — corrected in `/review-pr`.**
    /// An earlier version of this doc claimed the `-0.0` entry "now
    /// exercises that canonicalisation"; it cannot, by construction, at
    /// `Rung::IDENTITY` — see [`Direction::new`]'s own doc for why.
    /// `a_negative_zero_direction_is_canonicalised` (below) is the test
    /// that actually observes it.
    #[test]
    fn rung_zero_leaves_every_base_value_unchanged() {
        for base in [0.0, 1.0, -1.0, 1e-300, 1e300, 42.5, -0.0] {
            for p in [
                0.0,
                -0.0,
                Direction::P_MAX,
                -Direction::P_MAX,
                Direction::P_MAX / 3.0,
                -Direction::P_MAX / 3.0,
            ] {
                let d = Direction::new(p, Direction::P_MAX)
                    .unwrap_or_else(|| unreachable!("p is within [-P_MAX, P_MAX] by construction"));
                assert_eq!(
                    perturb(base, Rung::IDENTITY, d).to_bits(),
                    base.to_bits(),
                    "base {base}, p {p}: identity rung changed the value's bits"
                );
            }
        }
    }

    /// **The assertion `rung_zero_leaves_every_base_value_unchanged` cannot
    /// make, added in `/review-pr`.** That test's `p = -0.0` entry passes
    /// whether or not `Direction::new` canonicalises the sign, because
    /// `Rung::IDENTITY`'s multiplier is `1.0` regardless — see
    /// [`Direction::new`]'s own doc. This asserts on the constructor's
    /// output directly, where the canonicalisation is actually observable.
    #[test]
    fn a_negative_zero_direction_is_canonicalised() {
        let d = Direction::new(-0.0, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("-0.0 is within [-P_MAX, P_MAX] by construction"));
        assert_eq!(
            d.get().to_bits(),
            0.0_f64.to_bits(),
            "-0.0 survived into a Direction uncanonicalised"
        );
    }

    /// The worst-case corner of the excursion-bound proof in [`perturb`]'s
    /// own doc: `rung = MAX` (`strength = 1.0`) and `p = -P_MAX` exactly —
    /// the closest this mechanism can come to annihilating a constant.
    ///
    /// **Calls [`perturb`] itself, found recomputing its formula inline in
    /// `/review-pr`.** An earlier version wrote `1.0 + max_rung.strength() *
    /// worst_p.get()` by hand — pinning a copy of the expression rather than
    /// the function, so a change to `perturb`'s own shape (caught elsewhere
    /// by the pinned grid) would not necessarily be caught here too.
    /// `perturb(1.0, ..)` is exact against the hand-written form:
    /// multiplication by `1.0` is the identity, so the substitution changes
    /// nothing about what either assertion checks.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "CLAUDE.md: tests may assert exactly. Every step is error-free — 0.5 and \
                  -0.5 are exactly representable, strength(MAX) is exactly 1.0 (MAX/MAX), \
                  1.0 * -0.5 is exact (multiplication by 1.0 is the identity), and 1.0 - 0.5 \
                  is exact by Sterbenz's lemma — so an epsilon here (as an earlier version of \
                  this test used) tolerates up to 16 ulp of drift in strength() while the \
                  message claims exactness; same argument as \
                  identity_has_zero_strength_and_every_other_rung_is_positive, above"
    )]
    fn the_multiplier_never_reaches_zero_at_the_worst_case_corner() {
        let max_rung = Rung::new(Rung::MAX)
            .unwrap_or_else(|| unreachable!("Rung::MAX is in 0..=Rung::MAX by construction"));
        let worst_p = Direction::new(-Direction::P_MAX, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("-P_MAX is within [-P_MAX, P_MAX] by construction"));
        let multiplier = perturb(1.0, max_rung, worst_p);
        assert_eq!(multiplier, 0.5, "the worst-case corner is not exactly 0.5");
        assert!(
            multiplier > 0.0,
            "the multiplier reached zero or went negative"
        );
        // And the positive corner, for completeness: 1.0 + 1.0*0.5 = 1.5.
        let best_p = Direction::new(Direction::P_MAX, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("P_MAX is within [-P_MAX, P_MAX] by construction"));
        let upper = perturb(1.0, max_rung, best_p);
        assert_eq!(upper, 1.5, "the best-case corner is not exactly 1.5");
    }

    /// The tolerance on `P(identity)` and the corpus size it is measured
    /// over — shared by [`the_binomial_band_pins_p_identity_as_a_number`]
    /// and by [`the_defect_reading_is_distinguishable_from_the_correct_one`],
    /// so widening the band cannot silently outrun its own discriminability
    /// check. Found decoupled in `/review-pr`: the discriminability test
    /// previously asserted against a hardcoded `0.04` unrelated to this
    /// value — mutation-verified to let a widened band (`0.06`) and the
    /// named `2*MAX+1` defect both pass together, undetected by either
    /// test, catchable only by the exact pinned sequence below.
    const P_IDENTITY_BAND: f64 = 0.02;

    /// See [`P_IDENTITY_BAND`].
    const P_IDENTITY_N: u32 = 20_000;

    /// Decision 10's other required test: `P(identity) = 1/(MAX + 1)`,
    /// pinned as a number over repeated draws — not an `O(m)`-order check,
    /// which both `next_range(m + 1)` and the `next_range(2*m + 1)` defect
    /// would satisfy identically.
    #[test]
    fn the_binomial_band_pins_p_identity_as_a_number() {
        let hits_count = (0..u64::from(P_IDENTITY_N))
            .filter(|&seed| Rung::draw(seed).is_identity())
            .count();
        let hits = u32::try_from(hits_count)
            .unwrap_or_else(|_| unreachable!("hits cannot exceed P_IDENTITY_N, which fits in u32"));
        let want_p = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        let observed_p = f64::from(hits) / f64::from(P_IDENTITY_N);
        // Binomial standard error at N=20,000: sqrt(p(1-p)/N) ~= 0.00222 at
        // p ~= 0.1111. P_IDENTITY_BAND = 0.02 is ~9.0 sigma on that scale —
        // corrected in `/review-pr`, where an earlier version of this
        // comment claimed "6-sigma (~0.013)", measurably wrong for the
        // value actually shipped. Generous enough to survive sampling
        // noise while tight enough that the 2*MAX+1 defect (p ~= 0.0588,
        // a ~23.5-sigma miss on this scale — also corrected from an
        // earlier "~5-sigma" claim) fails loudly rather than marginally.
        // See `the_defect_reading_is_distinguishable_from_the_correct_one`
        // for the exact, computed discriminability bound, and note that an
        // off-by-one dropping the `+ 1` (`next_range(MAX)`, p ~= 0.128)
        // passes this band — it is caught only by
        // `the_rung_sequence_over_0_to_64_is_pinned`, which is load-bearing
        // for that case, not belt-and-braces.
        assert!(
            (observed_p - want_p).abs() < P_IDENTITY_BAND,
            "observed P(identity) = {observed_p} over {P_IDENTITY_N} draws, want {want_p} +/- \
             {P_IDENTITY_BAND} ({hits} hits) — check next_range(MAX + 1) vs next_range(2*MAX + 1)"
        );
    }

    /// **A guard neither of `borbax-universe`'s two pinned goldens can
    /// provide, and the reason it exists rather than relying on them.**
    /// `every_domain_stream_is_pinned` (`borbax-rng`) pins three raw
    /// `next_u64` words per domain — it proves the *stream* moved, not
    /// that `next_range(MAX + 1)`'s rejection-sampling correctly turns
    /// those words into `0..=MAX`. The binomial test above is
    /// deliberately statistical, so it only catches a defect large enough
    /// to move `P(identity)` outside its band — a change to `m`, an
    /// off-by-one in the rejection zone, or a swapped comparison could
    /// leave `P(identity)` looking right by coincidence over 20,000 draws
    /// while every individual rung is wrong. An exact pin over a small,
    /// enumerable seed range closes that gap the way `the_universe_digest_is_pinned`
    /// closes the equivalent gap for element generation.
    #[test]
    fn the_rung_sequence_over_0_to_64_is_pinned() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |x: u64| {
            h ^= x;
            h = h.wrapping_mul(0x0100_0000_01b3);
        };
        for seed in 0..64 {
            mix(u64::from(Rung::draw(seed).0));
        }
        assert_eq!(
            h, 0xce60_cec5_2ce3_ee01,
            "the rung sequence over seeds 0..64 moved — say whether this is a deliberate \
             change to Rung::MAX, to the draw's rejection-sampling, or a bug, and regenerate \
             deliberately if the first"
        );
    }

    /// A guard against the specific defect Decision 10 names by name:
    /// `next_range(2*MAX + 1)` gives `2*MAX + 1` outcomes and roughly
    /// halves `P(identity)`. This measures that the two readings are
    /// actually distinguishable at the corpus size
    /// [`the_binomial_band_pins_p_identity_as_a_number`] uses, not merely
    /// that they are different formulas.
    ///
    /// **Derived from [`P_IDENTITY_BAND`]/[`P_IDENTITY_N`], not a
    /// standalone literal — the bug this closes in `/review-pr`.** An
    /// earlier version asserted `(correct_p - defect_p).abs() > 0.04`, a
    /// constant with no relationship to the band it was meant to certify.
    /// Mutation-verified: widening the band to `0.06` and applying the
    /// named defect together left *both* this test and the band test
    /// green — the exact "the counter doesn't count what its message
    /// says" shape CLAUDE.md records from Task 8. The threshold here is
    /// the band plus six binomial standard errors at the shared corpus
    /// size: if the defect's miss does not clear that, the band test
    /// cannot reliably tell the two apart from sampling noise alone.
    #[test]
    fn the_defect_reading_is_distinguishable_from_the_correct_one() {
        let correct_p = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        let defect_p = 1.0 / f64::from(2 * u32::from(Rung::MAX) + 1);
        let se = (correct_p * (1.0 - correct_p) / f64::from(P_IDENTITY_N)).sqrt();
        let threshold = P_IDENTITY_BAND + 6.0 * se;
        assert!(
            (correct_p - defect_p).abs() > threshold,
            "the band ({P_IDENTITY_BAND}) is too wide to reliably see the 2*MAX+1 defect at \
             N={P_IDENTITY_N}: gap {} does not clear threshold {threshold}",
            (correct_p - defect_p).abs()
        );
    }

    /// The grid the two tests below share: every legal rung, ten base
    /// magnitudes, eighteen directions spanning the legal band.
    ///
    /// **Every rung, not just `0` and `MAX`.** Before these tests existed,
    /// `perturb` was called by exactly one test and only ever at
    /// `Rung::IDENTITY`, where `strength` is `0.0` and every
    /// algebraically-equal reformulation is bit-identical — so nothing in
    /// the workspace could see a change to the expression, and nothing at
    /// all could see a change to `Rung::strength`'s ladder between its two
    /// pinned endpoints.
    ///
    /// **Short-mantissa values cannot discriminate, so the grid avoids
    /// them.** With `base == 1.0`, `base * (1.0 + s*p)` and the distributed
    /// `base + base*s*p` are bit-identical for every `s` and `p` — a grid of
    /// round numbers pins nothing about the shape. `1.0` is kept as a
    /// readable anchor; the discrimination comes from the other nine.
    /// Likewise `k/17` is inexact for every `k` but `0` and `17`, which is
    /// why the directions are derived rather than hand-picked.
    fn pin_grid() -> impl Iterator<Item = (f64, Rung, Direction)> {
        const BASES: [f64; 10] = [
            1.0,
            -1.0,
            42.5,
            -13.7,
            0.1,
            6.022_140_76e23,
            1.380_649e-23,
            2.997_924_58e8,
            1e300,
            1e-300,
        ];
        (0..=Rung::MAX).flat_map(move |r| {
            let rung =
                Rung::new(r).unwrap_or_else(|| unreachable!("r is 0..=Rung::MAX by construction"));
            BASES.into_iter().flat_map(move |base| {
                (0..=17).map(move |k| {
                    let p = -Direction::P_MAX + f64::from(k) / 17.0;
                    let d = Direction::new(p, Direction::P_MAX).unwrap_or_else(|| {
                        unreachable!("p is within [-P_MAX, P_MAX] by construction")
                    });
                    (base, rung, d)
                })
            })
        })
    }

    /// **What [`perturb`]'s doc calls "pinned, not merely preferred", pinned
    /// — and `Rung::strength`'s ladder with it, which nothing else covers.**
    ///
    /// Measured before this test was written, on this commit: replacing
    /// `strength`'s body with a geometric ladder `(r/MAX)^2` — which
    /// preserves `strength(0) == 0.0`, `strength(MAX) == 1.0` and
    /// `strength(r) > 0.0` between, everything the rest of this module
    /// checks — passed all six gate legs, changing the perturbation
    /// magnitude in seven of the nine universes. Replacing `perturb`'s body
    /// with the distributed `base + base * strength * p` was caught only by
    /// `rung_zero_leaves_every_base_value_unchanged`, and only through a
    /// single sign-of-zero cell (`base == -0.0`, `p < 0.0`) — a value-neutral
    /// accident, not a test of the shape.
    ///
    /// The fused form `base * strength.mul_add(p, 1.0)` is *not* covered
    /// here and does not need to be: `clippy.toml`'s `disallowed_methods`
    /// and `xtask`'s `BANNED_CALLS` both reject `mul_add` outright, verified
    /// by mutation on this commit (clippy: "use of a disallowed method";
    /// xtask: a `§13.1: platform transcendental` hit). A test could not
    /// name it anyway — the textual scan reads test code too, and an
    /// `#[expect]` silences clippy while leaving `xtask` firing.
    #[test]
    fn the_perturb_expression_and_strength_ladder_are_pinned_bit_for_bit() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut cells = 0u32;
        for (base, rung, p) in pin_grid() {
            h ^= perturb(base, rung, p).to_bits();
            h = h.wrapping_mul(0x0100_0000_01b3);
            cells += 1;
        }
        assert_eq!(cells, 1620, "the pinned grid changed size");
        assert_eq!(
            h, 0x0839_1cd4_8bbc_19c1,
            "perturb's output moved. This is a physics change, not a tidy-up: say whether \
             the expression shape, Rung::strength's ladder, Rung::MAX or Direction::P_MAX \
             changed deliberately, and regenerate every downstream golden with it"
        );
    }

    /// **The bookkeeping check for the golden above**, in the same spirit as
    /// `the_defect_reading_is_distinguishable_from_the_correct_one`: that the
    /// grid can actually *see* the reformulation [`perturb`]'s doc names,
    /// rather than merely covering a lot of cells. A grid trimmed back to
    /// rungs `0` and `MAX`, or to round base values, would pass the golden
    /// and discriminate nothing.
    ///
    /// It is a readable second guard too: rewrite `perturb` to the
    /// distributed shape and this drops to zero, firing with a message that
    /// names the shape instead of with a moved hash.
    #[test]
    fn the_pinned_grid_still_sees_the_distributed_reshaping() {
        let differing = pin_grid()
            .filter(|&(base, rung, p)| {
                let (s, pv) = (rung.strength(), p.get());
                perturb(base, rung, p).to_bits() != (base + base * s * pv).to_bits()
            })
            .count();
        // 294 of 1620 cells when this was written. The floor is 0, not 294:
        // a deliberate change to `strength`'s ladder moves the count without
        // weakening the grid, and should fire the golden above alone.
        assert!(
            differing > 0,
            "the grid no longer distinguishes `base + base * strength * p` from \
             `base * (1.0 + strength * p)` — either perturb IS the distributed form now, \
             or the grid lost the full-mantissa base values that make the two differ \
             (with base == 1.0 they are bit-identical at every rung and direction)"
        );
    }
}
