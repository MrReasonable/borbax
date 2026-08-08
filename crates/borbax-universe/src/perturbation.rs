//! The shared perturbation mechanism (issue #26, Decision 10) — one integer
//! rung per universe, and a bounded per-constant direction, composed as
//! `base * (1.0 + strength * p)`.
//!
//! **Scope note, and why this module has no consumer yet.** [`Rung`] and
//! [`Direction`] are complete, independently tested types, and [`perturb`]
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
//! **`Rung`'s own visibility, stated explicitly rather than left to fall
//! out of an unrelated struct's convention (issue #26's plan, P6).**
//! [`Rung`] is `pub` today, and has to be: `crate::naming::IdentityWitness::new`
//! takes one as a parameter, so any caller able to reach that constructor
//! must be able to name the type it takes. **The trap is a future
//! `Universe.rung` field, not `Rung` itself.** `crates/borbax-universe/src/lib.rs`'s
//! `pub struct Universe` gives every field `pub` by house style — a `rung`
//! field held to that same convention would let `universe.rung ==
//! Rung::IDENTITY` recover the identity predicate from `borbax-molecule`,
//! `borbax-reaction` or `borbax-ui` in one line, without ever calling
//! `IdentityWitness::new`. When Task 26.1 adds that field, it must be held
//! no more visible than the code that actually needs to read it requires —
//! break the house-style convention deliberately there, and say so, rather
//! than let a mechanical "every field `pub`" pass make the witness's
//! fallibility decorative.

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
pub struct Rung(u8);

impl Rung {
    /// The top of the rung's legal range, `m` in the plan's notation.
    pub const MAX: u8 = 8;

    /// `Rung(0)` — the identity configuration.
    pub const IDENTITY: Self = Self(0);

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
    pub const fn new(r: u8) -> Option<Self> {
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
    #[must_use]
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        reason = "next_range(u64::from(Self::MAX) + 1) returns a value in 0..=MAX, and MAX is a \
                  u8 literal, so the result fits u8 by construction — proven by the domain of \
                  next_range's input, not merely assumed"
    )]
    pub fn draw(seed: u64) -> Self {
        let mut stream = Stream::new(seed, Domain::PerturbationRung, 0);
        let r = stream.next_range(u64::from(Self::MAX) + 1);
        Self(r as u8)
    }

    /// Every rung except the identity, `1..=MAX` — for Decision 9's
    /// exhaustive test (`IdentityWitness::new(r).is_none()` for every `r`
    /// here), so the corpus is derived from `MAX` rather than a hardcoded
    /// `1..=RUNGS` literal that could drift from it.
    pub fn all_off_identity() -> impl Iterator<Item = Self> {
        (1..=Self::MAX).map(Self)
    }

    /// `true` for `Rung(0)` — the identity configuration.
    #[must_use]
    pub const fn is_identity(self) -> bool {
        self.0 == 0
    }

    /// This rung's magnitude, `rung / MAX`, in `[0.0, 1.0]` — derived, never
    /// stored. "Store the integer rung, not the derived `f64` strength"
    /// (Decision 10): the integer is what every test and every digest
    /// should key on, so a strength recomputed here can never disagree with
    /// one stored somewhere else.
    #[must_use]
    pub fn strength(self) -> f64 {
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
pub struct Direction(f64);

impl Direction {
    /// The default excursion bound most migrated constants use. `w_shape`
    /// and `w_charge` are tighter, at `1/3` — see issue #26's plan, P7, for
    /// the numerical-soundness argument (a monotone Möbius transform of
    /// their ratio, degenerating at its extremes) that bound rests on.
    pub const P_MAX: f64 = 0.5;

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
    /// `float_cmp` suppression. Verified value-neutral through [`perturb`]
    /// too — the perturbed result is exactly `1.0` for every constructible
    /// rung regardless of which sign of zero `p` started as, so
    /// canonicalising here moves no perturbed value, only a would-be
    /// digest's bit pattern.
    #[must_use]
    pub fn new(v: f64, bound: f64) -> Option<Self> {
        if v.is_finite() && bound <= Self::P_MAX && v.abs() <= bound {
            Some(Self(v + 0.0))
        } else {
            None
        }
    }

    /// The raw value.
    #[must_use]
    pub const fn get(self) -> f64 {
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
/// **The expression shape is pinned, not merely preferred.** `base * (1.0 +
/// strength * p)` and the algebraically-equal `base + base * strength * p`
/// are different `f64` values — floating-point multiplication does not
/// distribute over addition — so redistributing this after Task 26.1 wires
/// the first migrated constant moves every golden downstream of it, the
/// same class of silent physics change §13.1 exists to catch. Commuting
/// either multiply (`p * strength`, or `(1.0 + ..) * base`) is bit-identical
/// and therefore safe; distributing is not.
#[must_use]
pub fn perturb(base: f64, rung: Rung, p: Direction) -> f64 {
    base * (1.0 + rung.strength() * p.get())
}

#[cfg(test)]
#[expect(
    clippy::unwrap_used,
    reason = "CLAUDE.md: tests may unwrap freely — a fixture that does not construct is a broken \
              test, and panicking says so at the point of the mistake"
)]
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
        assert_eq!(Rung::new(Rung::MAX).unwrap().strength(), 1.0);
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

    /// Decision 10's own required test: at `rung == 0`, every base value is
    /// bit-identical, over the full finite range of `p` at both signs —
    /// **compared by `to_bits()`, not by value**, because an earlier
    /// version of this test used `assert_eq!` on the `f64`s directly, and
    /// `assert_eq!(-0.0, 0.0)` passes: the one input in the `base` list
    /// chosen to exercise sign-of-zero preservation was exactly the input
    /// the comparison could not see a mismatch in. `-0.0` is included in
    /// the `p` list too, per Decision 10's own text ("at both zero signs of
    /// `p`") — [`Direction::new`] canonicalises it to `+0.0` as of this
    /// review round, so this entry now exercises that canonicalisation
    /// rather than sign propagation through `perturb`, and is kept either
    /// way since both are worth covering.
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
                let d = Direction::new(p, Direction::P_MAX).unwrap();
                assert_eq!(
                    perturb(base, Rung::IDENTITY, d).to_bits(),
                    base.to_bits(),
                    "base {base}, p {p}: identity rung changed the value's bits"
                );
            }
        }
    }

    /// The worst-case corner of the excursion-bound proof in [`perturb`]'s
    /// own doc: `rung = MAX` (`strength = 1.0`) and `p = -P_MAX` exactly —
    /// the closest this mechanism can come to annihilating a constant.
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
        let max_rung = Rung::new(Rung::MAX).unwrap();
        let worst_p = Direction::new(-Direction::P_MAX, Direction::P_MAX).unwrap();
        let multiplier = 1.0 + max_rung.strength() * worst_p.get();
        assert_eq!(multiplier, 0.5, "the worst-case corner is not exactly 0.5");
        assert!(
            multiplier > 0.0,
            "the multiplier reached zero or went negative"
        );
        // And the positive corner, for completeness: 1.0 + 1.0*0.5 = 1.5.
        let best_p = Direction::new(Direction::P_MAX, Direction::P_MAX).unwrap();
        let upper = 1.0 + max_rung.strength() * best_p.get();
        assert_eq!(upper, 1.5, "the best-case corner is not exactly 1.5");
    }

    /// Decision 10's other required test: `P(identity) = 1/(MAX + 1)`,
    /// pinned as a number over repeated draws — not an `O(m)`-order check,
    /// which both `next_range(m + 1)` and the `next_range(2*m + 1)` defect
    /// would satisfy identically.
    #[test]
    fn the_binomial_band_pins_p_identity_as_a_number() {
        const N: u32 = 20_000;
        let hits_count = (0..u64::from(N))
            .filter(|&seed| Rung::draw(seed).is_identity())
            .count();
        let hits = u32::try_from(hits_count)
            .unwrap_or_else(|_| unreachable!("hits cannot exceed N, which fits in u32"));
        let want_p = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        let observed_p = f64::from(hits) / f64::from(N);
        // Binomial standard error at N=20,000: sqrt(p(1-p)/N) ~= 0.0022 at
        // p ~= 0.111. A 6-sigma band (~0.013) is generous enough to survive
        // sampling noise while tight enough that the 2*MAX+1 defect (which
        // would give p ~= 0.0588, a ~5-sigma miss on its own scale) fails
        // loudly rather than marginally.
        let band = 0.02;
        assert!(
            (observed_p - want_p).abs() < band,
            "observed P(identity) = {observed_p} over {N} draws, want {want_p} +/- {band} \
             ({hits} hits) — check next_range(MAX + 1) vs next_range(2*MAX + 1)"
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
    /// actually distinguishable at the corpus size the test above uses,
    /// not merely that they are different formulas.
    #[test]
    fn the_defect_reading_is_distinguishable_from_the_correct_one() {
        let correct_p = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        let defect_p = 1.0 / f64::from(2 * u32::from(Rung::MAX) + 1);
        assert!(
            (correct_p - defect_p).abs() > 0.04,
            "the two readings are too close to distinguish at this corpus size"
        );
    }
}
