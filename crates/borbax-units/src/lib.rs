//! Invented units (spec §5, G4).
//!
//! Borbax never uses real-world units. Temperature is measured in thermals,
//! energy in quanta, distance in spans, time in world-years, mass in mass
//! units defined by the generated periodic table. These are distinct types
//! and the compiler refuses to mix them — which is the point. It prevents
//! anyone, including us, from quietly reasoning about the simulation as
//! though it modelled something real.
//!
//! # The separation, proved rather than asserted
//!
//! These doctests live in the crate documentation and not in the `tests`
//! module below, because rustdoc compiles under `cfg(doctest)` and never
//! `cfg(test)`: a `compile_fail` block attached to a `#[test] fn` inside
//! `#[cfg(test)] mod tests` is not collected and never runs. `cargo test
//! --doc` reporting `0 tests` is the only sign, and it looks like success.
//!
//! A `compile_fail` block passes when the code fails to build for **any**
//! reason, so on its own it proves very little: a misspelled type name
//! satisfies all four below at once. `Mass` is a further trap of the same
//! shape — its field is private, so `Mass(1.0)` written out here fails with
//! `E0423` before it ever reaches the question being asked, which is why
//! these go through [`Mass::from_raw`].
//!
//! Two things guard against that. First, the block immediately below must
//! *compile*, and it exercises every name the failing blocks use. Second,
//! each error code was read off the compiler by inverting the block —
//! removing `compile_fail` and looking at what rustc actually said — rather
//! than predicted. That inversion is the real check, and it is a manual one:
//! **rustdoc does not enforce the error code.** Measured on the pinned
//! toolchain, `compile_fail,E0308` on a block that errors with `E0369`
//! passes silently, with no warning. The annotations are therefore accurate
//! documentation of the intended failure and nothing stronger; anyone
//! editing a block here should invert it and read the error again.
//!
//! The permitted operations — same unit added, scaled by a bare scalar,
//! divided by its own kind to give a ratio:
//!
//! ```
//! use borbax_units::{Mass, Quanta, Span, Thermal, WorldYear};
//! let _ = Quanta(1.0) + Quanta(2.0);
//! let _ = Span(3.0) * 2.0;
//! let _: f64 = Thermal(6.0) / Thermal(2.0);
//! let _ = WorldYear(1.0) - WorldYear(0.5);
//! let _ = Mass::from_raw(1024) + Mass::from_raw(2048);
//! ```
//!
//! Energy is not mass (G4):
//!
//! ```compile_fail,E0308
//! use borbax_units::{Mass, Quanta};
//! let _ = Quanta(1.0) + Mass::from_raw(1024);
//! ```
//!
//! Distance is not time:
//!
//! ```compile_fail,E0308
//! use borbax_units::{Span, WorldYear};
//! let _ = Span(1.0) + WorldYear(1.0);
//! ```
//!
//! Two quantities of the *same* unit cannot be multiplied either — the
//! result would have a unit this universe has no name for:
//!
//! ```compile_fail,E0308
//! use borbax_units::Thermal;
//! let _ = Thermal(2.0) * Thermal(3.0);
//! ```
//!
//! And [`Mass`] specifically refuses to be scaled by a real, because that
//! would leave the fixed-point grid that makes conservation exact:
//!
//! ```compile_fail,E0369
//! use borbax_units::Mass;
//! let _ = Mass::from_raw(1024) * 2.0;
//! ```
//!
//! `Thermal` also has no `Default`, because zero thermals is a physical
//! extreme rather than the value anyone means to leave unset:
//!
//! ```compile_fail,E0599
//! use borbax_units::Thermal;
//! let _ = Thermal::default();
//! ```

use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

pub mod det_math;

/// Generate a float-backed unit.
///
/// Extra derives are listed per unit rather than fixed here, because the
/// first one that had to differ differed immediately: zero quanta is the
/// additive identity and a sensible `Default`, whereas zero thermals is a
/// distinguished physical extreme that nobody means to ask for. A
/// `#[derive(Default)] struct BeakerConfig { temperature: Thermal, .. }`
/// filled in with `..Default::default()` reads as unremarkable and produces
/// a beaker in which nothing ever happens — and not loudly, because the
/// Arrhenius rate clamps the temperature away from zero before the exponent
/// rather than producing a NaN that would announce itself.
macro_rules! unit {
    ($name:ident, $doc:literal $(, $extra_derive:ident)*) => {
        #[doc = $doc]
        ///
        /// Scaling by a dimensionless scalar is allowed and like-over-like
        /// division yields a ratio, but there is deliberately no way to
        /// combine this with any other unit — see the crate docs.
        // repr(transparent) guarantees this is passed in a floating-point
        // register exactly as a bare `f64` is, which is what makes the
        // profile note about the newtypes being zero-cost a guarantee rather
        // than a hope. Without it we would be relying on unspecified
        // repr(Rust) behaviour that merely happens to work today: a plain
        // newtype has the same layout on this toolchain, measured, which is
        // precisely why the reliance would be invisible.
        //
        // It does *not* buy `&[Quanta] -> &[f64]` reinterpretation. That
        // needs `unsafe`, and the workspace forbids it — `bytemuck::Pod`
        // does not compile here. The constraint it does impose is one
        // non-ZST field, forever; a `Span { value, uncertainty }` proposal
        // has to meet that on arrival.
        #[repr(transparent)]
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd $(, $extra_derive)*)]
        pub struct $name(pub f64);

        impl $name {
            /// The additive identity.
            pub const ZERO: Self = Self(0.0);

            /// The underlying magnitude, with the unit discarded.
            ///
            /// Every call is a place where the type system stops helping, so
            /// they should be rare and local — display, serialisation, and
            /// arithmetic that genuinely leaves the unit system.
            #[must_use]
            pub const fn get(self) -> f64 {
                self.0
            }

            /// Whether the magnitude is neither infinite nor NaN.
            #[must_use]
            pub const fn is_finite(self) -> bool {
                self.0.is_finite()
            }

            /// A total order over every value including NaN, for sorting.
            ///
            /// Sorting by a float key is common in the simulation and §13.4
            /// requires every such sort to be reproducible, which
            /// [`PartialOrd`] cannot promise: `partial_cmp` returns `None` on
            /// a NaN, and the workspace denies `unwrap`, so the obvious
            /// spelling is not available. Without this the fallback is
            /// `sort_by(|a, b| a.get().total_cmp(&b.get()))` at every site —
            /// which works, but puts a `get()` into ordinary code and makes
            /// the NaN position a decision taken twelve times instead of
            /// once. It is taken here: NaN sorts last.
            ///
            /// A total order is necessary and not sufficient. §13.4 also
            /// requires a tie-break on a stable ID wherever two keys can be
            /// equal, and that belongs to the caller.
            #[must_use]
            pub fn total_cmp(&self, other: &Self) -> std::cmp::Ordering {
                self.0.total_cmp(&other.0)
            }
        }

        impl Add for $name {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                Self(self.0 + rhs.0)
            }
        }

        impl AddAssign for $name {
            fn add_assign(&mut self, rhs: Self) {
                self.0 += rhs.0;
            }
        }

        impl Sub for $name {
            type Output = Self;
            fn sub(self, rhs: Self) -> Self {
                Self(self.0 - rhs.0)
            }
        }

        impl Neg for $name {
            type Output = Self;
            fn neg(self) -> Self {
                Self(-self.0)
            }
        }

        /// Scaling by a dimensionless scalar is allowed; multiplying two
        /// quantities of the same unit is not, because the result would
        /// have a unit we have no name for.
        impl Mul<f64> for $name {
            type Output = Self;
            fn mul(self, rhs: f64) -> Self {
                Self(self.0 * rhs)
            }
        }

        impl Div<f64> for $name {
            type Output = Self;
            fn div(self, rhs: f64) -> Self {
                Self(self.0 / rhs)
            }
        }

        /// Dividing like by like yields a dimensionless ratio.
        ///
        /// Division by [`Self::ZERO`] follows IEEE-754 and produces an
        /// infinity or a NaN rather than panicking. That is the right
        /// behaviour for a simulation — but it also means a poisoned value
        /// can travel a long way before anyone notices, so callers dividing
        /// by something that could be zero should say what they expect.
        impl Div for $name {
            type Output = f64;
            fn div(self, rhs: Self) -> f64 {
                self.0 / rhs.0
            }
        }

        impl Sum for $name {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::ZERO, Add::add)
            }
        }
    };
}

// `Thermal` deliberately has no `Default` — see the macro's doc comment.
// Write `Thermal::ZERO` where zero really is meant, so it shows up in a diff.
unit!(
    Thermal,
    "Temperature, in thermals. Not Kelvin, not Celsius."
);
unit!(Quanta, "Energy, in quanta. Not joules.", Default);
unit!(Span, "Distance, in spans. Not metres.", Default);
unit!(
    WorldYear,
    "Time, in world-years — defined by the generated world's own orbit.",
    Default
);

/// Mass, in mass units defined by the generated periodic table.
///
/// **Fixed-point, not floating-point** — deliberately the odd one out.
///
/// The global constraints require that a reaction's products have *exactly*
/// the mass of its reactants. Naive `f64` summation over many terms cannot
/// deliver that: you get a tolerance argument instead of a guarantee, and
/// the tolerance has to grow as molecules get larger. Integer arithmetic
/// makes conservation exact by construction, so the test is `==` rather
/// than `(a - b).abs() < eps` — and it stays exact on every platform,
/// which matters for spec §13.4.
///
/// The representation is 1/1024 of a nominal mass unit. Generated atomic
/// masses are quantised onto this grid once, at universe-generation time,
/// so every molecular mass is an exact sum of exact atomic masses forever
/// after.
///
/// There is no `Mul<f64>` and no `Div<f64>`, which is the whole point:
/// scaling a mass by an arbitrary real would leave the grid and hand back
/// the tolerance argument the fixed point exists to avoid. Stoichiometry
/// goes through [`Mass::checked_mul_count`], and anything that genuinely
/// needs a real — a rate, a display string — goes through
/// [`Mass::to_f64`] and does not come back.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Mass(i64);

/// A real value could not be placed on the fixed-point mass grid.
///
/// Carries the offending value, because by the time this surfaces the
/// interesting question is which generated quantity went out of range — and
/// `inf` versus `NaN` versus `1e300` point at different bugs.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
#[error(
    "mass {value} is not on the fixed-point grid: must be finite and at most ±{}",
    Mass::MAX_MAGNITUDE
)]
pub struct MassRangeError {
    /// The value that could not be quantised.
    pub value: f64,
}

impl Mass {
    /// Sub-units per mass unit. A power of two so that quantisation is
    /// exact in binary and the grid has no representation error of its own.
    pub const SCALE: i64 = 1024;

    /// The additive identity.
    pub const ZERO: Self = Self(0);

    /// [`Self::SCALE`] as a float, for the two conversion functions.
    ///
    /// Written out rather than cast. 1024 is exactly representable either
    /// way, so this buys nothing numerically — it keeps `as` and a
    /// precision-loss suppression out of code that is otherwise about
    /// exactness, where a reader should not have to stop and check. The
    /// two spellings are kept in step by `quantises_on_the_declared_grid`.
    const SCALE_F64: f64 = 1024.0;

    /// Largest magnitude [`Self::from_f64_quantised`] can place on the grid:
    /// `i64::MAX / SCALE`, which for a 1024-sub-unit grid is `2^53 - 1` and so
    /// is itself exactly representable in `f64`. Spelled as a literal for the
    /// same reason as [`Self::SCALE_F64`]; `the_representable_range_matches_the_grid`
    /// pins it to `SCALE`.
    ///
    /// Public because [`Self::from_f64_quantised`] holds callers to it, and a
    /// bound a caller cannot read is not a contract.
    pub const MAX_MAGNITUDE: f64 = 9_007_199_254_740_991.0;

    /// Quantise a real value onto the fixed-point grid.
    ///
    /// Called only at universe-generation time, never in a hot path.
    /// `round_ties_even` is used rather than `round` because it is
    /// symmetric — `round` biases away from zero, which would make
    /// generated mass distributions very slightly heavier than intended.
    ///
    /// # Errors
    ///
    /// [`MassRangeError`] if `v` is NaN, infinite, or of greater magnitude
    /// than [`Self::MAX_MAGNITUDE`].
    ///
    /// **This is fallible on purpose, and it was a `debug_assert!` first.**
    /// The `as` conversion below saturates on an out-of-range value and maps
    /// NaN to zero, so an assertion that vanishes in release turns a loud
    /// failure into a quiet wrong answer — and `Mass(0)` is the worst kind,
    /// because a zero mass is conserved perfectly forever and every
    /// conservation test downstream still passes on a universe containing a
    /// bogus element. The profile that erases the assertion is also the one
    /// that mints the goldens (`cargo run -p borbax-cli --release --
    /// goldens --emit`), which is the whole argument in one line.
    ///
    /// The NaN input is not hypothetical: generated masses are drawn as
    /// `base + normal() * spread`, and `det_math::ln`'s own documentation
    /// warns that `ln(0)` is `-inf` — which is the Box–Muller path a normal
    /// draw takes.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "truncation to the grid is the operation, and the range is checked above"
    )]
    #[expect(
        clippy::as_conversions,
        reason = "no fallible f64 -> i64 conversion exists; the guard above is that conversion"
    )]
    pub fn from_f64_quantised(v: f64) -> Result<Self, MassRangeError> {
        if !v.is_finite() || v.abs() > Self::MAX_MAGNITUDE {
            return Err(MassRangeError { value: v });
        }
        Ok(Self((v * Self::SCALE_F64).round_ties_even() as i64))
    }

    /// Approximate real value, for display and for rate arithmetic.
    /// Never round-trip through this to perform a conservation check —
    /// compare `Mass` values directly.
    #[must_use]
    #[expect(
        clippy::cast_precision_loss,
        reason = "display and rate arithmetic only; conservation is checked on the i64"
    )]
    #[expect(
        clippy::as_conversions,
        reason = "no lossless i64 -> f64 conversion exists"
    )]
    pub fn to_f64(self) -> f64 {
        self.0 as f64 / Self::SCALE_F64
    }

    /// The raw count of sub-units. Exact, and the only thing worth
    /// hashing or serialising.
    #[must_use]
    pub const fn raw(self) -> i64 {
        self.0
    }

    /// Reconstruct from a raw sub-unit count, as produced by [`Self::raw`].
    ///
    /// The contract is `Mass::from_raw(m.raw()) == m`, and deserialisation is
    /// the use case. It takes any `i64` and validates nothing, deliberately:
    /// every `i64` is *on* the grid by construction, so unlike a `Mul<f64>`
    /// this cannot produce an off-grid value and does not weaken the reason
    /// that operator is absent. What it can do is exceed
    /// [`Self::MAX_MAGNITUDE`], so a value that did not come from
    /// [`Self::raw`] has not been range-checked by anybody.
    #[must_use]
    pub const fn from_raw(v: i64) -> Self {
        Self(v)
    }

    /// Scale by a whole number of copies — stoichiometry. Returns `None`
    /// on overflow rather than wrapping, because a silently wrapped mass
    /// would look like a conservation violation somewhere far away.
    #[must_use]
    pub fn checked_mul_count(self, n: u32) -> Option<Self> {
        self.0.checked_mul(i64::from(n)).map(Self)
    }
}

/// Unchecked, unlike [`Mass::checked_mul_count`], and the asymmetry is
/// deliberate. `checked_mul_count`'s multiplier is externally supplied and
/// unbounded; addition's operands are two masses that
/// [`Mass::from_f64_quantised`] has already bounded to ±2^53 sub-units, so a
/// single addition cannot leave the range — which makes that constructor's
/// range check load-bearing for this impl and not merely adjacent to it.
/// Making these checked would return `Option` into a workspace that denies
/// `unwrap`, poisoning every call site to guard an unreachable case.
impl Add for Mass {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl AddAssign for Mass {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl Sub for Mass {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Neg for Mass {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Sum for Mass {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, Add::add)
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    reason = "CLAUDE.md: tests may unwrap freely. Every unwrap here is on a \
              constant chosen to be on the grid, so a failure is a bug in the test"
)]
mod tests {
    use super::*;

    #[test]
    fn same_units_add() {
        assert_eq!(Quanta(2.0) + Quanta(3.0), Quanta(5.0));
    }

    /// The plan wrote this against `Mass`, which cannot hold a float and has
    /// no `Mul<f64>` by design. Scaling belongs to the float units.
    #[test]
    fn scales_by_scalar() {
        assert_eq!(Span(2.0) * 3.0, Span(6.0));
    }

    #[test]
    fn sums_over_iterator() {
        let total: Quanta = [Quanta(1.0), Quanta(2.0), Quanta(4.0)].into_iter().sum();
        assert_eq!(total, Quanta(7.0));
    }

    #[test]
    fn orders() {
        assert!(Thermal(300.0) > Thermal(100.0));
    }

    #[test]
    fn zero_is_additive_identity() {
        assert_eq!(Quanta::ZERO + Quanta(4.0), Quanta(4.0));
    }

    #[test]
    fn like_over_like_is_a_bare_ratio() {
        let ratio: f64 = Quanta(6.0) / Quanta(4.0);
        assert!((ratio - 1.5).abs() < f64::EPSILON);
    }

    /// Quantise a literal known to be on the grid. Every call site here is a
    /// constant chosen to be representable, so a failure is a bug in the test.
    fn grid(v: f64) -> Mass {
        Mass::from_f64_quantised(v).unwrap()
    }

    /// Pins the two spellings of the grid together without a cast: if
    /// `SCALE_F64` ever drifts from `SCALE`, one unit does not land on
    /// `SCALE` sub-units.
    #[test]
    fn quantises_on_the_declared_grid() {
        assert_eq!(grid(1.0).raw(), Mass::SCALE);
        assert_eq!(grid(0.5).raw(), Mass::SCALE / 2);
        assert_eq!(grid(-1.0).raw(), -Mass::SCALE);
    }

    /// The reason the range check is not a `debug_assert!`. Every one of
    /// these produces a plausible-looking `Mass` from a saturating `as`
    /// cast if the guard is compiled out, and `NaN` produces `Mass::ZERO`,
    /// which is conserved perfectly forever and breaks no other test here.
    ///
    /// Run this under `--release` too — that is the profile that mints the
    /// goldens, and the profile the assertion used to vanish in.
    #[test]
    fn values_off_the_grid_are_refused_in_every_profile() {
        for bad in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            1e300,
            Mass::MAX_MAGNITUDE * 2.0,
        ] {
            assert!(
                Mass::from_f64_quantised(bad).is_err(),
                "{bad} was quantised instead of refused"
            );
        }
        // The boundary itself is representable, and one sub-unit inside it.
        assert!(Mass::from_f64_quantised(Mass::MAX_MAGNITUDE).is_ok());
        assert!(Mass::from_f64_quantised(-Mass::MAX_MAGNITUDE).is_ok());
    }

    /// `round` would send both of these away from zero. Ties-to-even sends
    /// them to the nearer even sub-unit, which is unbiased — over a whole
    /// generated periodic table the difference is a systematically heavier
    /// universe.
    #[test]
    fn quantisation_breaks_ties_to_even_not_away_from_zero() {
        let half_a_sub_unit = 0.5 / Mass::SCALE_F64;
        assert_eq!(grid(half_a_sub_unit).raw(), 0);
        assert_eq!(grid(3.0 * half_a_sub_unit).raw(), 2);
    }

    /// The property the fixed point exists for. Summed in two different
    /// orders, because that is the failure `f64` would have.
    #[test]
    fn mass_conservation_is_exact_regardless_of_order() {
        let parts: Vec<Mass> = [0.1, 0.2, 0.3, 12.007, 0.000_976_562_5]
            .into_iter()
            .map(grid)
            .collect();
        let forwards: Mass = parts.iter().copied().sum();
        let backwards: Mass = parts.iter().rev().copied().sum();
        assert_eq!(forwards, backwards);
        // Subtracting every part back off returns to exactly zero — the
        // property `f64` cannot offer without a tolerance.
        assert_eq!(parts.iter().copied().fold(forwards, Sub::sub), Mass::ZERO);
    }

    /// Pins [`Mass::MAX_MAGNITUDE`] to the grid without an `as` cast in
    /// either direction: if `SCALE` moves, the first assertion fails.
    #[test]
    fn the_representable_range_matches_the_grid() {
        assert_eq!(i64::MAX / Mass::SCALE, 9_007_199_254_740_991);
        // Compared as bits: exact, and `clippy::float_cmp` is right that `==`
        // on floats needs a reason every other time.
        assert_eq!(
            Mass::MAX_MAGNITUDE.to_bits(),
            9_007_199_254_740_991.0_f64.to_bits()
        );
    }

    #[test]
    fn stoichiometry_scales_exactly_and_refuses_to_wrap() {
        assert_eq!(grid(2.5).checked_mul_count(4), Some(grid(10.0)));
        assert_eq!(Mass::from_raw(i64::MAX / 2).checked_mul_count(3), None);
    }

    #[test]
    fn mass_round_trips_through_raw() {
        let m = grid(37.25);
        assert_eq!(Mass::from_raw(m.raw()), m);
    }

    /// NaN sorts last and the order does not depend on where it started —
    /// which `partial_cmp` cannot promise, and §13.4 requires.
    #[test]
    fn total_cmp_is_a_reproducible_order_even_with_a_nan() {
        let mut a = [Span(2.0), Span(f64::NAN), Span(1.0)];
        let mut b = [Span(f64::NAN), Span(1.0), Span(2.0)];
        a.sort_by(Span::total_cmp);
        b.sort_by(Span::total_cmp);
        let bits: Vec<u64> = a.iter().map(|s| s.get().to_bits()).collect();
        let other: Vec<u64> = b.iter().map(|s| s.get().to_bits()).collect();
        assert_eq!(bits, other);
        assert!(a.last().is_some_and(|s| s.get().is_nan()));
    }

    /// The platform calls below are the one legitimate reason to call a
    /// banned method: proving `det_math` is wired up at all means comparing
    /// it against the platform. The suppression is a visible line rather
    /// than a whole exempted module, which is the point of enforcing §13.1
    /// with a resolved-path lint instead of a grep.
    #[test]
    #[expect(
        clippy::disallowed_methods,
        reason = "§13.1: comparing det_math against the platform is the test"
    )]
    fn det_math_agrees_with_the_platform_to_within_a_few_ulp() {
        // Not equality — the whole point is that they may differ. This asserts
        // we have not wired up something wrong, not that they match.
        for x in [-3.0, -0.5, 0.0, 0.5, 1.0, 7.5_f64] {
            assert!((det_math::exp(x) - x.exp()).abs() < 1e-12 * x.exp().max(1.0));
        }
        for x in [0.25, 1.0, 2.0, 1000.0_f64] {
            assert!((det_math::ln(x) - x.ln()).abs() < 1e-12);
        }
        for x in [-2.0, 0.0, 0.7, 3.5_f64] {
            assert!((det_math::sin(x) - x.sin()).abs() < 1e-12);
            assert!((det_math::cos(x) - x.cos()).abs() < 1e-12);
        }
        assert!((det_math::powf(2.0, 10.0) - 1024.0).abs() < 1e-9);
    }

    #[test]
    fn det_math_is_deterministic_within_a_process() {
        assert_eq!(det_math::exp(1.0).to_bits(), det_math::exp(1.0).to_bits());
    }
}
