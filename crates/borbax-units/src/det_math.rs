//! Portable transcendentals (spec §13.1, §13.4).
//!
//! Backed by the `libm` crate — rust-lang's pure-Rust port of MUSL's libm.
//! Its `exp`, `log`, `sin`, `cos` and `pow` carry no architecture dispatch;
//! the `use_arch` paths in that crate are `sqrt`, `fma`, the rounding
//! functions, and x87 variants gated on 32-bit x86 without SSE2, none of which
//! is a CI target or reachable from the five functions wrapped here.
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
