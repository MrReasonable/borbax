//! Portable transcendentals (spec §13.1, §13.4).
//!
//! Backed by the `libm` crate — rust-lang's pure-Rust port of MUSL's libm.
//!
//! **The criterion is "is the result exactly specified by IEEE-754", not
//! "does it avoid arch dispatch"** — the same restatement `clippy.toml`
//! carries, and for the same reason. `exp`, `log`, `sin`, `cos` and `pow` are
//! portable because they have no dispatch **on the three architectures the
//! §13.4 matrix builds** — which is not the same as "no dispatch at all", and
//! the difference is what a future `libm` bump gets judged against. Measured on
//! 0.2.16: `exp`, `exp2` and `exp10` route to x87 inline assembly on 32-bit x86
//! without SSE2, via `use_arch_required`, which bypasses even
//! `force-soft-float` — and `libm`'s own `arch/i586.rs` documents the result as
//! possibly 1 ulp off and *hardware-dependent*. We do not build that target.
//! Adding a 32-bit leg to the matrix would make `exp` a portability hazard the
//! same day, and this paragraph is the reason someone would notice. `acos` is
//! portable for a *different* reason: it reaches `super::sqrt`, which **is**
//! arch-dispatched on both CI architectures, to the hardware instruction —
//! and IEEE-754 specifies `sqrt` to be correctly rounded, so the hardware and
//! generic paths are bit-identical.
//!
//! Applying the wrong test would reject `hypot` and `cbrt` for no reason, or
//! wave through something that reaches a genuinely unspecified helper. An
//! earlier version of this paragraph counted the wrappers and asserted none
//! reached a dispatched path; the count went stale the moment `acos` was
//! added, which is why the criterion is stated instead of a number.
//!
//! Apple's libm and glibc genuinely disagree in the last bits of `exp`, `ln`,
//! `sin` and `cos`, and glibc versions disagree with each other. Measured on
//! this project's own machine, `libm::exp` and the platform's differ on
//! **9.9%** of sampled inputs; `ln` 2.5%, `sin` 4.4%, `cos` 4.2%, `pow` 4.7%.
//! None of them is wrong: IEEE-754 does not specify these to be correctly
//! rounded.
//!
//! It does specify `+ - * / sqrt`, which is why those stay native — with one
//! caveat worth writing down rather than rediscovering. What is specified is
//! the *value*: a NaN's sign bit and payload are not, and were measured to
//! differ between aarch64 and x86-64 for `0.0/0.0`, `inf - inf`, `inf * 0.0`
//! and `(-1.0).sqrt()` alike. That is harmless until a NaN reaches a hashed
//! value, at which point the state serialiser has to canonicalise it or
//! assert it cannot happen (Tasks 20-21).
//!
//! **This lives in `borbax-units` because it must sit below every caller.**
//! `borbax-rng` and `borbax-universe` both need it and both sit below
//! `borbax-molecule`; anywhere higher and this file is unreachable from two of
//! the crates it exists to serve.
//!
//! A direct `.exp()` or `.cos()` on an `f64` anywhere in library code is a
//! determinism bug. It is enforced by `clippy::disallowed_methods`, whose path
//! list is in the workspace `clippy.toml`, with a textual second pass in
//! `cargo xtask`. This file needs no exemption from either: the wrappers call
//! `libm::exp` as a *free function*, and the banned path is `f64::exp`, the
//! inherent method. The chokepoint is the `libm` crate rather than one
//! privileged file, so there is no path an exemption could later widen.
//!
//! Adding a function here is cheap and is the right response to needing one.
//! The ban covers every non-IEEE-specified float method, including several
//! with no wrapper below (`atan2`, `tanh`, `cbrt`, `hypot`, …); when one of
//! those is genuinely wanted, it gets a wrapper here rather than an entry
//! removed there.

/// `e^x`.
#[must_use]
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// Natural logarithm. Callers must guarantee `x > 0` — `ln(0)` is `-inf` and
/// will silently poison whatever consumes it (see Task 18's clock).
#[must_use]
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// Sine, with `x` in radians.
#[must_use]
pub fn sin(x: f64) -> f64 {
    libm::sin(x)
}

/// Cosine, with `x` in radians.
#[must_use]
pub fn cos(x: f64) -> f64 {
    libm::cos(x)
}

/// Arc cosine, in radians. **Clamps its argument**, and that is the feature.
///
/// Every caller in this project will be extracting an angle from a quantity
/// that is mathematically in `[-1, 1]` and numerically is not: a dot product
/// of two unit vectors, or `(tr(R) - 1) / 2` for a rotation matrix built by
/// accumulating products, where `tr(R)` routinely comes out as
/// `3.0000000000000004`. Unclamped, that is `acos(1.0000000000000002)` — NaN,
/// and a NaN whose sign bit differs between the CI architectures.
///
/// So the clamp lives here rather than at the call site. The alternative is
/// `ln`'s contract one line down — "callers must guarantee" — which is only
/// safe while there is exactly one caller who remembered. A caller who does
/// not know the trace can exceed 3 cannot get this wrong.
/// It is **not** a NaN guard, and should not be. `clamp` propagates NaN, so a
/// NaN argument still gives a NaN result — which is right: that means the
/// caller already had a NaN, which is a different bug, and swallowing it here
/// would hide it. What this fixes is the rounding overshoot, which is not a
/// bug anywhere and cannot be avoided by the caller.
///
/// `clamp` rather than `max`/`min`, which §13.1 bans for their `±0.0`
/// behaviour; `clamp` uses ordinary comparisons and has no such ambiguity.
/// Measured: it lowers to `maxsd`/`minsd` on x86-64 and `fcmp`/`fcsel` on
/// aarch64 — branchless on both — and costs ~13% of a 3.4 ns call.
///
/// **The clamp is load-bearing, measured on this repo's own data.** Computing
/// `(tr(R) - 1) / 2` for all 3 600 compositions of `rotation_matrices()`:
/// **685 of them exceed 1** and would return NaN unclamped, worst overshoot
/// 8.9e-16. The identity element alone does it at D = 42.
///
/// **Accuracy limit, so a downstream tolerance can be checked against it.**
/// `acos` loses half its digits at both ends (`dθ/dx = -1/√(1-x²)`): angles
/// from this are trustworthy to **~1.4e-8 rad absolute**, and θ below ~1.5e-8
/// is lost entirely. No downstream tolerance on an angle derived this way may
/// be tighter than 1e-7. The icosahedral group's own angles are safe — the
/// two ill-conditioned ones, 0 and π, come out *exactly* — and the smallest
/// angle between sample directions is 0.277 rad, where the error is ~8e-16.
/// If an angle ever enters a result-affecting path near 0 or π, switch to
/// `atan2` on the axis magnitude, which is exact everywhere; that wrapper does
/// not exist yet and should not be added speculatively.
#[must_use]
pub fn acos(x: f64) -> f64 {
    libm::acos(x.clamp(-1.0, 1.0))
}

/// `x^y` for runtime `y`.
///
/// For small integer exponents write the multiplication out instead. The
/// honest reason is that `powi`'s multiply tree is lowered by LLVM and we have
/// not verified that lowering across all three targets and every toolchain
/// this project will pin — not that a divergence has been measured. It has
/// not: on aarch64, `powi` matched hand-written left-to-right binary
/// exponentiation on every one of ~20 000 sampled cases for exponents 2..=9.
/// The ban stands because `x * x` costs nothing and is clearer, which is a
/// good enough reason on its own; dressing a precaution up as a measurement is
/// how the `glam`/`nalgebra` claim in `CLAUDE.md` became a doctrine error.
#[must_use]
pub fn powf(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

/// Real cube root.
///
/// **Its portability argument is not `exp`'s, and assuming it is would be the
/// wrong reading.** The other wrappers here are safe because `libm` computes
/// them in pure Rust from operations IEEE-754 specifies exactly. `libm::cbrt`
/// is different: it is the core-math port and it calls `f64::fma`, which on
/// aarch64+neon lowers to hardware `fmadd` and on x86-64 is *runtime
/// dispatched* between `vfmadd213sd` and a software fallback. That looks like
/// exactly the hazard this module exists to prevent, so it was measured rather
/// than argued: over 200 000 samples in `[1e-6, 1000]` the digest is
/// **identical** on aarch64 and x86-64 — and identical with the dispatch taking
/// *different branches*, because under Rosetta `is_x86_feature_detected!("fma")`
/// reports false and the soft path runs. Which is what IEEE-754 predicts:
/// `fusedMultiplyAdd` is a specified, correctly-rounded operation. Same
/// argument [`acos`] makes for `sqrt`, one step further.
///
/// The platform `f64::cbrt` is a genuinely different function and is banned in
/// `clippy.toml`: over the same samples it disagrees with this on **15 223** of
/// them on aarch64 and 15 230 on x86-64, with *different digests* — a 7.6%
/// divergence rate between two macOS builds differing only in architecture.
///
/// Do not substitute `powf(x, 1.0 / 3.0)`. It is accurate to under an ulp, but
/// `1.0 / 3.0` is not one third, so it is exact on no perfect cube above 27 —
/// `powf(64.0, 1.0 / 3.0)` is `3.999_999_999_999_999_6` and `powf(343.0, …)` is
/// `6.999_999_999_999_999_1`. A permanent cube root in a result-affecting path
/// should be spelled as one.
#[must_use]
pub fn cbrt(x: f64) -> f64 {
    libm::cbrt(x)
}

#[cfg(test)]
#[expect(
    clippy::float_cmp,
    reason = "CLAUDE.md: tests may assert exactly; exactness on perfect cubes is the \
              property under test, so an epsilon would defeat it"
)]
mod tests {
    /// The reason [`cbrt`] exists rather than `powf(x, 1.0 / 3.0)`. `1.0 / 3.0`
    /// is not one third, so the substitute is exact on no perfect cube above
    /// 27 — a sparse, silent one-ulp discontinuity for no benefit.
    #[test]
    fn cbrt_is_exact_on_perfect_cubes_where_powf_is_not() {
        let mut powf_wrong = 0;
        for root in 1_u32..=40 {
            let cube = f64::from(root * root * root);
            assert_eq!(
                super::cbrt(cube),
                f64::from(root),
                "cbrt({cube}) should be exactly {root}"
            );
            if super::powf(cube, 1.0 / 3.0) != f64::from(root) {
                powf_wrong += 1;
            }
        }
        assert!(
            powf_wrong > 0,
            "powf(x, 1.0/3.0) was exact on every cube tested — if that is now \
             true, this wrapper's justification needs rewriting, not deleting"
        );
    }

    /// `libm`'s cube root and the platform's are different functions. This
    /// pins that we call the portable one; on a platform where they agree the
    /// test still passes, so it cannot fail spuriously.
    #[test]
    fn cbrt_matches_libm_not_the_platform() {
        for i in 1_u32..=2000 {
            let x = f64::from(i) * 0.37;
            assert_eq!(super::cbrt(x), libm::cbrt(x));
        }
    }
}
