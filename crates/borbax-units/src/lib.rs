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

/// The positive quiet NaN. Written as bits rather than `f64::NAN` because the
/// bit pattern is the whole point here, and `f64::NAN`'s is not guaranteed.
const CANONICAL_NAN_BITS: u64 = 0x7ff8_0000_0000_0000;

/// Negative zero, the other value whose sign bit changes an IEEE `totalOrder`
/// without changing what the number *is*.
const NEG_ZERO_BITS: u64 = 0x8000_0000_0000_0000;

/// `x` with any sign-bit ambiguity removed, so that ordering it is portable.
///
/// Two values need this and no others. A runtime-produced NaN carries an
/// architecture-dependent sign bit (clear on aarch64, set on x86-64), and
/// `-0.0` carries a sign bit while comparing equal to `0.0`. Both would move
/// under `f64::total_cmp` without meaning anything different.
///
/// **It normalises the whole NaN class, not only the sign** — payload and
/// signalling-ness go too, because it replaces the value outright rather than
/// masking a bit. That is deliberate and stronger than it looks: `det_math`
/// records that NaN *payloads* also differ between the architectures for
/// `0.0/0.0`, `inf - inf`, `inf * 0.0` and `(-1.0).sqrt()`. Anyone
/// "simplifying" this to a sign-bit mask would lose payload canonicalisation
/// and every current test would still pass.
///
/// What it deliberately leaves alone is what IEEE-754 already specifies
/// exactly: infinities and subnormals. The one hazard in this class it cannot
/// address is flush-to-zero — if `FTZ`/`DAZ` were ever set, a subnormal would
/// become `±0.0` on one platform and stay subnormal on the other. Rust sets
/// neither and this workspace is `forbid(unsafe_code)` with no C dependency,
/// so nothing here can set them; noted so nobody concludes this helper is
/// total.
#[inline]
const fn canonical_sign(x: f64) -> f64 {
    if x.is_nan() {
        f64::from_bits(CANONICAL_NAN_BITS)
    } else if x.to_bits() == NEG_ZERO_BITS {
        0.0
    } else {
        x
    }
}

/// The bit pattern to hash or serialise for a bare `f64`.
///
/// The method on each unit is the one to reach for when you have a unit. This
/// exists because the state hash will be dominated by values that are *not*
/// units: `LATTICE_SPAN` is a bare `f64` today, and the signature array is the
/// single most-hashed quantity in the system.
/// Without an exported form, a caller holding an `f64` writes `.to_bits()` and
/// the canonicalisation is silently skipped for exactly the data that matters
/// most.
///
/// Typing those quantities as [`Span`] shrinks the surface but does not remove
/// it. **Two of the named quantities have left the list.** `Signature.r` was on
/// it and no longer is — Task 9 typed it, along with the signature's own
/// `SHELL`. `BindConsts::ideal_gap` was on it and **never existed**:
/// `UniverseConsts::ideal_gap` has been a [`Span`] since Task 5, and Task 10's
/// `BindConsts` does not carry the field at all, because the separation is
/// derived per pair rather than drawn. `LATTICE_SPAN` (Task 13) is the one
/// still bare. What stays here regardless is the *character* channel,
/// which is a ratio in `[-1, 1]` with no unit to carry a method, and is why
/// [`canonical_cmp`] below exists.
#[must_use]
#[inline]
pub const fn canonical_bits(x: f64) -> u64 {
    canonical_sign(x).to_bits()
}

/// The reproducible sort order, for a bare `f64`.
///
/// The free/method pairing here is [`canonical_bits`]', for the same reason and
/// with the same split: reach for `Unit::canonical_cmp` when you have a unit,
/// and for this when the quantity genuinely has none. Every unit's method
/// delegates here, so there is one implementation of the decision rather than
/// one per unit plus a copy at each bare-`f64` call site.
///
/// **NaN sorts last and `-0.0` ties with `0.0`, on every platform** — see
/// the private `canonical_sign` for what that costs and what it deliberately
/// leaves alone. Every caveat on `Unit::canonical_cmp` applies unchanged, including
/// the important one: a reproducible order is necessary and not sufficient, so
/// §13.4's ID tie-break is still the caller's job, over a *larger* tie set than
/// `f64::total_cmp` would produce.
///
/// The first caller is `borbax-molecule`'s signature: §8.2 canonicalises a
/// signature to the lexicographically smallest of its 60 rotations, and the
/// tie-break on that ordering has to run over the surface-character channel,
/// which is a ratio in `[-1, 1]` and has no unit to carry the method. The
/// alternatives were wrapping a dimensionless ratio in [`Span`] — a G4 type
/// error to dodge a lint — or an `#[expect]` on `f64::total_cmp` resting on a
/// proof that the channel is never NaN, which a later edit invalidates
/// silently.
///
/// **Do not build `Ord` from this.** It reports two NaNs `Equal` while `==`
/// reports them unequal, so an `Ord` delegating here would violate the
/// `Ord`/`PartialEq` consistency contract.
#[must_use]
#[inline]
#[expect(
    clippy::disallowed_methods,
    reason = "§13.4: after Unit::canonical_cmp was made to delegate here, this is the \
              only caller of total_cmp in this crate — canonicalising the inputs first \
              is what makes it safe. Three sites elsewhere carry their own #[expect] \
              with their own reasons (geodesic.rs, packing.rs x2), so this is the \
              sanctioned wrapper rather than the only exemption in the workspace"
)]
pub fn canonical_cmp(a: f64, b: f64) -> core::cmp::Ordering {
    canonical_sign(a).total_cmp(&canonical_sign(b))
}

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
///
/// **Every operator below carries `#[inline]`, and it is load-bearing.**
/// These are non-generic, so without the attribute their bodies live only in
/// this crate's codegen unit and every downstream crate emits a call. Measured
/// at `opt-level = 1` — the dev profile — `a + b * s` was 12 instructions with
/// two calls and register spills, against 3 for bare `f64`; with `#[inline]`,
/// 3. Raising the opt-level does not fix it (opt-1 through opt-3 all sit at
/// ~9 000 ns on a proxy kernel against ~1 700 for the inlined form); only LTO
/// or this attribute does, and the dev profile has no LTO. The zero-cost claim
/// in the root `Cargo.toml` was false for the dev profile until these were
/// added. Release was always fine, because `lto = "thin"` closes it anyway.
///
/// It cannot move a golden: `#[inline]` grants no licence to reassociate float
/// arithmetic, and the release asm was verified byte-identical with and
/// without.
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
            #[inline]
            pub const fn get(self) -> f64 {
                self.0
            }

            /// Whether the magnitude is neither infinite nor NaN.
            #[must_use]
            #[inline]
            pub const fn is_finite(self) -> bool {
                self.0.is_finite()
            }

            /// A reproducible order over every value including NaN, for
            /// sorting. **NaN sorts last and `-0.0` ties with `0.0`, on every
            /// platform.**
            ///
            /// Sorting by a float key is common in the simulation and §13.4
            /// requires every such sort to be reproducible, which
            /// [`PartialOrd`] cannot promise: `partial_cmp` returns `None` on
            /// a NaN, and the workspace denies `unwrap`, so the obvious
            /// spelling is not available. Nor does it put a `get()` into
            /// ordinary code, or make the NaN position a decision taken
            /// twelve times instead of once.
            ///
            /// **This is not `f64::total_cmp`, and the difference is the
            /// reason the method exists.** `f64::total_cmp` is IEEE-754
            /// `totalOrder`, which sorts on the sign bit first — and a
            /// *runtime-produced* NaN has its sign bit clear on aarch64 and
            /// set on x86-64. Measured on the pinned toolchain:
            ///
            /// ```text
            /// aarch64   0.0/0.0 = 0x7ff8...   f64::total_cmp -> [1, 2, NaN]
            /// x86-64    0.0/0.0 = 0xfff8...   f64::total_cmp -> [NaN, 1, 2]
            /// ```
            ///
            /// So the plain delegation shipped in Task 2 reordered the entire
            /// slice between the macOS leg and the Linux and Windows legs.
            /// Its test did not catch that, because `f64::NAN` is
            /// const-evaluated to the positive pattern on every target — the
            /// one NaN whose sign is portable. Canonicalising both signs here
            /// makes the documented contract true everywhere, and has the
            /// second benefit of agreeing with `==` about signed zeros, which
            /// `f64::total_cmp` does not.
            ///
            /// A reproducible order is necessary and not sufficient, and
            /// this method makes §13.4's ID tie-break **more** necessary
            /// rather than less. Under `f64::total_cmp` two keys tie only if
            /// bit-identical; here `+0.0`/`-0.0`, two NaN signs, and a quiet
            /// against a signalling NaN all tie. Demonstrated: a stable sort
            /// over the same multiset in two input permutations gives two
            /// different outputs, because stability preserves *input* order
            /// and input order is only canonical if the upstream is. So the
            /// tie-break belongs to the caller, on a strictly larger set of
            /// inputs than before.
            ///
            /// **Do not build `Ord` from this.** It reports two NaNs `Equal`
            /// while the derived `PartialEq` reports them unequal, so an
            /// `Ord` impl delegating here would violate the `Ord`/`PartialEq`
            /// consistency contract and `BTreeMap` behaviour would follow it
            /// into the weeds.
            ///
            /// This is a *sort* order. `max_by`/`min_by` over it still
            /// propagate NaN, so a NaN in a selection wins; the invariant that
            /// keeps that from mattering is that signature and affinity
            /// values are finite, asserted where they are built.
            ///
            /// The body delegates to the free [`canonical_cmp`] so the two
            /// spellings cannot drift: this one is reached when the quantity
            /// has a unit, that one when it genuinely has none, and there is
            /// one implementation of the decision rather than one per unit.
            #[must_use]
            #[inline]
            pub fn canonical_cmp(&self, other: &Self) -> core::cmp::Ordering {
                crate::canonical_cmp(self.0, other.0)
            }

            /// The bit pattern to hash or serialise — **use this, not
            /// `get().to_bits()`.**
            ///
            /// Same hazard as [`Self::canonical_cmp`], one step further on.
            /// A runtime NaN's sign bit is architecture-dependent and `-0.0`
            /// carries one while equalling `0.0`, so two runs that agree on
            /// every value can still hash differently. §13.6's golden matrix
            /// would report that as a simulation divergence.
            ///
            /// It costs nothing to have this ready before the state hash
            /// exists (Task 20) and a great deal to discover the need
            /// afterwards, from a red matrix with no visible cause.
            #[must_use]
            #[inline]
            pub const fn canonical_bits(self) -> u64 {
                canonical_sign(self.0).to_bits()
            }
        }

        impl Add for $name {
            type Output = Self;
            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self(self.0 + rhs.0)
            }
        }

        impl AddAssign for $name {
            #[inline]
            fn add_assign(&mut self, rhs: Self) {
                self.0 += rhs.0;
            }
        }

        impl Sub for $name {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self(self.0 - rhs.0)
            }
        }

        impl Neg for $name {
            type Output = Self;
            #[inline]
            fn neg(self) -> Self {
                Self(-self.0)
            }
        }

        /// Scaling by a dimensionless scalar is allowed; multiplying two
        /// quantities of the same unit is not, because the result would
        /// have a unit we have no name for.
        impl Mul<f64> for $name {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: f64) -> Self {
                Self(self.0 * rhs)
            }
        }

        impl Div<f64> for $name {
            type Output = Self;
            #[inline]
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
            #[inline]
            fn div(self, rhs: Self) -> f64 {
                self.0 / rhs.0
            }
        }

        impl Sum for $name {
            #[inline]
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
/// No `Default`, for a sharper version of the reason `Thermal` has none.
/// `Mass::default()` would be `Mass::ZERO`, and this type's own documentation
/// calls that the worst possible wrong answer — a zero mass is conserved
/// perfectly forever, so every conservation test downstream still passes on a
/// universe containing a bogus element. With the derive present,
/// `Mass::from_f64_quantised(v).unwrap_or_default()` collapses the fallible
/// constructor straight back into that value, and no workspace lint objects:
/// `unwrap_or_default` is not `unwrap`.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Mass(i64);

/// The whole overflow argument, checked by the compiler rather than asserted
/// in prose — and this block previously said that while asserting only two
/// thirds of it. Three reviewers caught it independently and one mutation-tested
/// it: `MAX_MAGNITUDE` could be quadrupled, both assertions still passed, and
/// only a runtime test noticed. That is the same shape as the defect this
/// constant was rewritten to fix, so the missing assertion is now here.
///
/// The three claims, all now checked at compile time:
/// 1. `MAX_SUMMABLE_TERMS` values of `MAX_RAW` sub-units fit in an `i64`.
/// 2. `MAX_RAW` is a whole number of mass units.
/// 3. `MAX_MAGNITUDE` is *exactly* `MAX_RAW` sub-units — the one that was
///    asserted in prose while `from_f64_quantised` gated on it alone.
const _: () = assert!(
    Mass::MAX_RAW
        .checked_mul(Mass::MAX_SUMMABLE_TERMS)
        .is_some()
);
const _: () = assert!(Mass::MAX_RAW % Mass::SCALE == 0);
#[expect(
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    reason = "both operands are powers of two below 2^53, so the product is exact in f64               and the conversion back is lossless; this assertion is the proof of that"
)]
const _: () = assert!((Mass::MAX_MAGNITUDE * Mass::SCALE_F64) as i64 == Mass::MAX_RAW);

/// A real value could not be placed on the fixed-point mass grid.
///
/// Carries the offending value, because by the time this surfaces the
/// interesting question is which generated quantity went out of range — and
/// `inf` versus `NaN` versus `1e300` point at different bugs.
///
/// `PartialEq` is hand-written over canonicalised bits, and the route here is
/// worth recording. It was derived; a reviewer pointed out that NaN is the
/// headline input this type exists for and `NaN != NaN`, so the error was not
/// equal to itself. The derive was removed — and the same reviewer then
/// measured the cost of removing it and withdrew their own proposal: absence
/// also blocks `assert_eq!` on the **`Ok`** arm (`Result: PartialEq` needs
/// both bounds) and blocks `#[derive(PartialEq)]` on any error enum that later
/// wraps this one, which `borbax-universe` will plausibly want.
///
/// Comparing `canonical_sign(value).to_bits()` beats both. It is reflexive for
/// NaN, agrees about `±0.0`, restores `assert_eq!` on both arms, and keeps
/// derives available upstack. Since the relation is a genuine equivalence,
/// `Eq` follows honestly rather than by assertion.
#[derive(Debug, Clone, Copy, thiserror::Error)]
#[error(
    "mass {value} is not on the fixed-point grid: must be finite and at most ±{}",
    Mass::MAX_MAGNITUDE
)]
pub struct MassRangeError {
    /// The value that could not be quantised.
    pub value: f64,
}

impl PartialEq for MassRangeError {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        canonical_bits(self.value) == canonical_bits(other.value)
    }
}

impl Eq for MassRangeError {}

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

    /// How many `Mass` values may be summed before the running total can
    /// leave `i64`.
    ///
    /// This exists because [`Self::MAX_MAGNITUDE`] is derived from it, not the
    /// other way round. The first version of this type bounded a single value
    /// at `i64::MAX / SCALE` and asserted in a comment that "a single addition
    /// cannot leave the range" — which was false by a factor of `SCALE`, and
    /// wrong in kind as well, because [`Sum`] folds an unbounded iterator and
    /// a molecular mass is a sum over atoms while a beaker's mass is a sum
    /// over molecules. Two values the constructor accepted summed to **−2**
    /// under release wrapping. Three independent reviewers found it.
    ///
    /// The budget is denominated in **atoms**, not molecules — a molecular
    /// mass is itself a sum, so §17's 10^6 molecules at ~10 atoms each is 10^7
    /// against 4.3x10^9, roughly 400x. An earlier version of this line
    /// compared against molecules and claimed four thousand times; same
    /// conclusion, wrong unit.
    pub const MAX_SUMMABLE_TERMS: i64 = 1 << 32;

    /// The value bound in sub-units — where the arithmetic actually happens,
    /// and therefore where the overflow argument has to be made.
    const MAX_RAW: i64 = 1 << 30;

    /// Largest magnitude [`Self::from_f64_quantised`] will place on the grid,
    /// in mass units: `MAX_RAW / SCALE`, so that [`Self::MAX_SUMMABLE_TERMS`]
    /// of them can be added without leaving `i64`. Both are powers of two, so
    /// this is exact in `f64` and the const assertion below is the whole proof.
    ///
    /// **This bounds a single *element*, not a molecule.** `from_f64_quantised`
    /// is a universe-generation entry point for atomic masses; molecular and
    /// beaker masses arrive through [`Sum`] and are bounded by
    /// [`Self::MAX_SUMMABLE_TERMS`] instead, which is why a 40-monomer chain
    /// exceeding this number is fine and quantising one through this
    /// constructor would not be.
    ///
    /// Reviewers measured Task 4's generator at ~1700-2500 mass units, giving
    /// ~400x headroom — but that script was never checked in and Task 4 is
    /// unwritten, so the number is not reproducible from this repository. The
    /// previous version of this doc replaced an unmeasurable claim on exactly
    /// those grounds and then made two more; saying so here rather than
    /// repeating the trick. **Task 4 owns the other half**: it must assert
    /// that no generated element mass exceeds a fraction of this constant, or
    /// nothing connects the bound to the thing that justifies it.
    ///
    /// Public because [`Self::from_f64_quantised`] holds callers to it, and a
    /// bound a caller cannot read is not a contract.
    pub const MAX_MAGNITUDE: f64 = 1_048_576.0;

    /// Quantise a real value onto the fixed-point grid.
    ///
    /// Called only at universe-generation time, never in a hot path.
    ///
    /// `round_ties_even` rather than `round`, and the reason first written
    /// here was not the real one. It claimed `round`'s away-from-zero bias
    /// would make generated mass distributions "very slightly heavier" —
    /// measured over 180 000 draws from the Task 4 generator, the two agree
    /// on **every** input and not one lands on a tie, because a tie needs a
    /// value dyadic with denominator 2048 and these come out of `powf` plus a
    /// Gaussian. The claim was unmeasurable, which is the same shape of
    /// unbacked justification the `powi` note in `clippy.toml` records.
    ///
    /// The real reason is that ties-to-even is the mode every other float
    /// operation in the program already uses, so quantising a value and doing
    /// arithmetic on it round the same way. That bites exactly where ties are
    /// reachable — dyadic inputs, which is what a half-mass or a test fixture
    /// literal is, and what `quantisation_breaks_ties_to_even_not_away_from_zero`
    /// exercises.
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
    // Missed by the pass that added the other twenty-one, because its
    // multi-line `#[expect]` block broke the pattern being matched. Measured
    // cross-crate at opt-level 1: 1.49x in a loop with independent work, and
    // *no difference* in a reduction-bound one — which is why a spot-check in
    // the obvious shape would have shown nothing. Its own doc says "rate
    // arithmetic", and that is the beaker step.
    #[inline]
    pub fn to_f64(self) -> f64 {
        self.0 as f64 / Self::SCALE_F64
    }

    /// The bit pattern to hash or serialise, matching the float units' method
    /// of the same name so a state hash needs one spelling rather than two.
    ///
    /// No canonicalisation is possible or needed — `i64` has exactly one
    /// representation per value, which is the whole argument for the fixed
    /// point. `cast_unsigned` rather than `as`: it is the lossless
    /// bit-reinterpretation, and `raw()` fed to a hasher would otherwise need
    /// a sign-losing conversion and a suppression at every call site.
    #[must_use]
    #[inline]
    pub const fn canonical_bits(self) -> u64 {
        self.0.cast_unsigned()
    }

    /// The raw count of sub-units. Exact, and the only thing worth
    /// hashing or serialising.
    #[must_use]
    #[inline]
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
    #[inline]
    pub const fn from_raw(v: i64) -> Self {
        Self(v)
    }

    /// Scale by a whole number of copies — stoichiometry. Returns `None`
    /// on overflow rather than wrapping, because a silently wrapped mass
    /// would look like a conservation violation somewhere far away.
    #[must_use]
    #[inline]
    pub fn checked_mul_count(self, n: u32) -> Option<Self> {
        self.0.checked_mul(i64::from(n)).map(Self)
    }
}

/// Unchecked, unlike [`Mass::checked_mul_count`], and the asymmetry is
/// deliberate. `checked_mul_count`'s multiplier is externally supplied and
/// unbounded; addition's is not.
///
/// The invariant is **not** "both operands came from
/// [`Mass::from_f64_quantised`]" — most additions in the simulation have
/// operands that are themselves sums, since a molecular mass is a `Sum` over
/// atoms and a beaker mass is a `Sum` over molecules, and neither is ever
/// re-validated against [`Mass::MAX_MAGNITUDE`]. Stating it that way was a
/// category error and two reviewers caught it.
///
/// The invariant that actually holds: **every `Mass` in circulation is a sum
/// of at most [`Mass::MAX_SUMMABLE_TERMS`] constructor outputs**, and the
/// const assertion beside the type proves that many fit. The budget is
/// denominated in *atoms*, not molecules — §17's 10^6 molecules at ~10 atoms
/// each is 10^7 against a 4.3x10^9 budget.
///
/// One caveat the budget does not cover: a fold mixing [`Add`] and [`Sub`] can
/// reach twice the accumulator magnitude, and twice the budget is exactly
/// `i64::MAX + 1`. The conservation pattern that motivates `Sub` — sum, then
/// subtract the same parts — returns toward zero and cannot reach it, so this
/// is contrived rather than live; the honest budget for a mixed signed fold is
/// 2^31.
///
/// Making these checked would return `Option` into a workspace that denies
/// `unwrap`, poisoning every call site to guard a case the constructor has
/// already excluded.
///
/// The one route past that argument is [`Mass::from_raw`], which validates
/// nothing by design. A `Mass` that did not come from `raw()` has been
/// bounded by nobody, and this operator will wrap on it in release.
impl Add for Mass {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl AddAssign for Mass {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl Sub for Mass {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl Neg for Mass {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Sum for Mass {
    #[inline]
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

    /// Pins the two spellings of the bound together without an `as` cast:
    /// `MAX_MAGNITUDE` mass units must be exactly `MAX_RAW` sub-units.
    #[test]
    fn the_representable_range_matches_the_grid() {
        assert_eq!(grid(Mass::MAX_MAGNITUDE).raw(), Mass::MAX_RAW);
        assert_eq!(grid(-Mass::MAX_MAGNITUDE).raw(), -Mass::MAX_RAW);
    }

    /// The regression test for the claim `impl Add for Mass` makes.
    ///
    /// Before this, `MAX_MAGNITUDE` was `i64::MAX / SCALE`, so two values the
    /// constructor *accepted* summed to `-2` mass units under release wrapping
    /// and panicked in debug. A test over realistic masses passes either way;
    /// this one has to use the largest values the constructor admits, and has
    /// to run in the profile that mints goldens.
    #[test]
    fn the_largest_admissible_masses_sum_without_wrapping() {
        let m = grid(Mass::MAX_MAGNITUDE);
        let doubled = m + m;
        assert!(
            doubled > m,
            "addition wrapped: {} -> {}",
            m.raw(),
            doubled.raw()
        );
        assert_eq!(doubled.raw(), 2 * Mass::MAX_RAW);
        assert_eq!(doubled - m, m);
    }

    /// `Sum` folds an unbounded iterator, so the pairwise claim is not the one
    /// the type actually ships. `MAX_SUMMABLE_TERMS` is 2^32, which is too
    /// many to iterate; the const assertion beside `Mass` proves the whole
    /// budget, and this checks the fold agrees at a size a test can run.
    #[test]
    fn summing_many_maximal_masses_stays_positive_and_exact() {
        const N: i64 = 4096;
        let m = grid(Mass::MAX_MAGNITUDE);
        let total: Mass = std::iter::repeat_n(m, usize::try_from(N).unwrap_or(0)).sum();
        assert_eq!(total.raw(), N * Mass::MAX_RAW);
        assert!(total > m);
    }

    /// The budget is a real ceiling, not a slogan: `MAX_SUMMABLE_TERMS` values
    /// at the bound must fit, and one more than the budget allows must not.
    #[test]
    fn the_summation_budget_is_the_tightest_true_statement() {
        let ceiling = Mass::MAX_RAW.checked_mul(Mass::MAX_SUMMABLE_TERMS);
        assert!(ceiling.is_some(), "the declared budget already overflows");
        // `x2`, not `x4`. Headroom is exactly 2.0 and not one unit more
        // (2^62 x 2 = 2^63 = i64::MAX + 1), so asserting `x4` overflows would
        // pass unchanged at half the declared budget — the test permitted
        // precisely the slack its own failure message warns about.
        assert!(
            Mass::MAX_RAW
                .checked_mul(Mass::MAX_SUMMABLE_TERMS)
                .and_then(|c| c.checked_mul(2))
                .is_none(),
            "the bound is looser than it needs to be — say so, or tighten it"
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

    /// A NaN the compiler could not fold. `f64::NAN` is const-evaluated to the
    /// *positive* pattern on every target, which is the one NaN whose sign bit
    /// is portable — so a test built on it passes on all three CI legs while a
    /// sign-bit-ordering bug is live. That is exactly what the first version
    /// of this test did.
    fn runtime_nan() -> f64 {
        std::hint::black_box(0.0_f64) / std::hint::black_box(0.0_f64)
    }

    #[test]
    fn the_runtime_nan_is_the_one_that_varies_by_platform() {
        // Not an assertion about which pattern this platform produces — that
        // is the thing that differs. It asserts the test is using a NaN the
        // optimiser did not replace with the constant, because if it ever
        // does, `canonical_cmp_sorts_nan_last` silently stops testing anything.
        assert!(runtime_nan().is_nan());
        assert!(
            runtime_nan().to_bits() == CANONICAL_NAN_BITS
                || runtime_nan().to_bits() == CANONICAL_NAN_BITS | NEG_ZERO_BITS,
            "unexpected NaN payload: 0x{:016x}",
            runtime_nan().to_bits()
        );
    }

    /// The distinguishing test for the sign-bit bug. On x86-64 this fails
    /// against a plain `f64::total_cmp` delegation and passes with the
    /// canonicalisation; on aarch64 both pass, which is why the CI matrix and
    /// not the local run is what proves it.
    #[test]
    fn canonical_cmp_sorts_nan_last_on_every_platform() {
        let mut v = [Span(2.0), Span(runtime_nan()), Span(1.0)];
        v.sort_by(Span::canonical_cmp);
        assert!(
            v.last().is_some_and(|s| s.get().is_nan()),
            "NaN did not sort last: {v:?}"
        );
        assert_eq!(
            v.first().map(|s| s.get().to_bits()),
            Some(1.0_f64.to_bits())
        );
    }

    #[test]
    fn canonical_cmp_does_not_depend_on_input_order() {
        let mut a = [Span(2.0), Span(runtime_nan()), Span(1.0)];
        let mut b = [Span(runtime_nan()), Span(1.0), Span(2.0)];
        a.sort_by(Span::canonical_cmp);
        b.sort_by(Span::canonical_cmp);
        let key = |v: &[Span; 3]| v.iter().map(|s| s.get().is_nan()).collect::<Vec<_>>();
        assert_eq!(key(&a), key(&b));
    }

    /// `f64::total_cmp` calls these `Less`; `==` calls them equal. Two orders
    /// that disagree is a trap for the ID tie-break §13.4 requires, so
    /// `canonical_cmp` agrees with `==`.
    #[test]
    fn canonical_cmp_ties_the_two_signed_zeros() {
        let neg = -Quanta::ZERO;
        assert_eq!(neg.get().to_bits(), NEG_ZERO_BITS, "test needs a real -0.0");
        assert_eq!(neg.canonical_cmp(&Quanta::ZERO), std::cmp::Ordering::Equal);
        assert_eq!(Quanta::ZERO.canonical_cmp(&neg), std::cmp::Ordering::Equal);
        // ... and that agrees with what `==` already said.
        assert_eq!(neg, Quanta::ZERO);
    }

    /// The whole point of the wrapper: a rotation matrix built by accumulating
    /// products routinely has `tr(R) = 3.0000000000000004`, so the argument
    /// exceeds 1 and unclamped `acos` returns NaN — with an
    /// architecture-dependent sign bit. This is the case that would have hit
    /// the renderer's axis-angle extraction.
    #[test]
    fn acos_clamps_rather_than_returning_nan() {
        for x in [1.0, 1.000_000_000_000_000_2, 1.5, f64::INFINITY] {
            assert!(det_math::acos(x).is_finite(), "acos({x}) was not finite");
        }
        for x in [-1.0, -1.000_000_000_000_000_2, -1.5, f64::NEG_INFINITY] {
            assert!(det_math::acos(x).is_finite(), "acos({x}) was not finite");
        }
        // The trace of the identity rotation, the ill-conditioned end.
        let theta = det_math::acos((3.000_000_000_000_000_4 - 1.0) / 2.0);
        assert!(theta.is_finite() && theta.abs() < 1e-7, "theta = {theta}");
    }

    /// Hashing `get().to_bits()` would expose the sign-bit hazard that
    /// `canonical_cmp` closes for ordering. The helper exists so Task 20's
    /// state hash cannot get it wrong.
    #[test]
    fn canonical_bits_erases_the_sign_ambiguity() {
        assert_eq!(Span(runtime_nan()).canonical_bits(), CANONICAL_NAN_BITS);
        assert_eq!((-Quanta::ZERO).canonical_bits(), 0.0_f64.to_bits());
        // Everything else is untouched.
        assert_eq!(Span(1.5).canonical_bits(), 1.5_f64.to_bits());
        assert_eq!(Span(-1.5).canonical_bits(), (-1.5_f64).to_bits());
    }

    /// The free function is what a caller holding a bare `f64` reaches for,
    /// and the state hash will be mostly bare `f64` until Task 8 types the
    /// length quantities. It must agree with the method exactly.
    #[test]
    fn the_free_and_method_forms_of_canonical_bits_agree() {
        for x in [0.0, -0.0, 1.5, -1.5, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(canonical_bits(x), Span(x).canonical_bits(), "{x}");
        }
        assert_eq!(
            canonical_bits(runtime_nan()),
            Span(runtime_nan()).canonical_bits()
        );
    }

    /// `Mass` needs no canonicalisation — `i64` has one representation per
    /// value — but it needs the *spelling*, so a state hash has one method
    /// name rather than two.
    #[test]
    fn mass_canonical_bits_round_trips_through_raw() {
        for raw in [0_i64, 1, -1, i64::MAX, i64::MIN] {
            let m = Mass::from_raw(raw);
            assert_eq!(m.canonical_bits(), raw.cast_unsigned());
            assert_eq!(Mass::from_raw(m.canonical_bits().cast_signed()), m);
        }
    }

    /// The reason the derive came off and a hand-written impl went on: NaN is
    /// the input this error exists for, and a derived `PartialEq` made it
    /// unequal to itself.
    #[test]
    fn a_mass_range_error_is_equal_to_itself_for_every_input_it_exists_for() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e300] {
            let e = Mass::from_f64_quantised(bad);
            assert_eq!(e, Mass::from_f64_quantised(bad), "{bad}");
            assert!(e.is_err());
        }
        // And the Ok arm compares too, which dropping the derive would block.
        assert_eq!(Mass::from_f64_quantised(1.0), Ok(Mass::from_raw(1024)));
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
