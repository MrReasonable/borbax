# `borbax-units` — units that cannot be mixed, and maths that cannot drift

Two jobs, both about preventing a whole class of mistake at compile time rather
than catching it later.

## Job 1: units the compiler enforces

Borbax measures things in invented units, never real ones (spec §5, G4):

| Quantity | Unit | Underlying type |
|---|---|---|
| Temperature | `Thermal` | `f64` |
| Energy | `Quanta` | `f64` |
| Distance | `Span` | `f64` |
| Time | `WorldYear` | `f64` |
| Mass | `Mass` | **`i64`**, fixed-point |

Each is a distinct type, so `Quanta + Thermal` is a **compile error**, not a
runtime surprise. There are tests in this crate that exist purely to confirm
certain expressions *fail* to compile.

**Why bother?** Partly to stop real-world reasoning leaking in — if you cannot
write "300 kelvin" you cannot quietly start treating the simulation as a model
of something real. But mostly because unit mix-ups are the classic way a physics
engine produces plausible, confident nonsense.

### Why `Mass` is a whole number

Every reaction must conserve mass **exactly**. Floating point cannot promise
that in general — the usual demonstration is that `0.1 + 0.2` gives
`0.30000000000000004` — because addition rounds and is not associative, so over
millions of reactions the crumbs accumulate and mass quietly appears or
vanishes.

Fixed-point integers have no crumbs. `Mass` stores a whole number of 1/1024ths,
so `a + b − b == a` — **for sums of up to 2³² constructor outputs**, which is the
budget the crate states and a compile-time assertion proves. Past it, addition
wraps silently in release. That is not hypothetical: this crate has already
shipped that bug once, with two in-range values summing to **−2**.

**And the deeper reason is not that 0.1 is inexact.** Borbax's masses are all
exact binary fractions by construction, so an `f64` would represent every one of
them perfectly and conserve just as well across the whole range. The real
argument for integers is that nothing in an `f64` keeps a value *on* the grid:
one multiplication by a non-dyadic factor and the mass is off-lattice forever,
with no error raised and no way to notice. An integer either lands on the
lattice or fails to compile.

## Job 2: maths that gives the same answer everywhere

Here is a fact that surprises most people: **`exp`, `ln`, `sin`, `cos` and
`powf` are not fully specified by the floating-point standard.**

IEEE-754 *requires* correct rounding — the last bit right — for `+`, `−`, `×`,
`÷`, `sqrt`, fused multiply-add, remainder, and rounding-to-integer (`ceil`,
`floor`, `rint`). Broadly: every operation whose exact answer is a rational
function of its inputs. Then it stops. The elementary functions are
*recommended*, not required, and are left to whoever writes the maths library.

(An earlier version of this README said the standard pins down "`+ − × ÷` and
`sqrt` … and then stops", which contradicted its own next section — that section
relies on `fma`, `rint`, `ceil` and `floor` being exactly specified, and they
are.)

So `(2.5_f64).exp()` can genuinely produce a *different final bit* on macOS than
on Linux. One bit, propagated through a few million simulation steps, is a
different universe.

### The fix: one door

Every transcendental in the entire workspace goes through `det_math.rs`, which
calls [`libm`](https://crates.io/crates/libm) — a pure-Rust implementation that
dispatches on processor architecture **only for operations IEEE-754 specifies
exactly**. Same answer everywhere, and see the subtlety below for why that is
the right way to state it.

(An earlier version said libm has "no platform-specific paths for these
functions". That is false for `exp`, and `CLAUDE.md` names it as *the wrong
reason to trust libm* — precisely because it would wave through a future version
that widened one. Stating the false version here, in the section a reader stops
at, defeated the protection.)

A direct `.exp()` or `.cos()` on an `f64` anywhere in library code is a
**determinism bug**, not a style preference. `clippy.toml` bans them by name so
the build fails rather than the goldens.

`sqrt` and the four arithmetic operations stay native — those *are* exactly
specified, so routing them through a library would cost speed and buy nothing.

### The subtlety worth knowing

`libm` **does** dispatch on processor architecture — it is not a myth, and
"it doesn't" is the wrong reason to trust it. What makes it safe is *which*
operations it dispatches on **for the three platforms Borbax builds**: `sqrt`
and `fma` on x86-64, plus `rint` on aarch64. All IEEE-required, so they give
identical results whichever path runs.

**That is a fact about our targets, not about `libm`.** Measured on 0.2.16, the
version `Cargo.lock` resolves: `exp`, `exp2` and `exp10` route to x87 inline
assembly on 32-bit x86 without SSE2, through a mechanism that bypasses even the
force-soft-float flag — and `libm`'s own source documents that path as possibly
1 ulp off and *hardware-dependent*. We do not build that target. Adding a 32-bit
leg to the golden matrix would make `exp` a portability hazard the same day.

The distinction matters because this paragraph is the rule a future `libm` bump
gets judged against. Stated as "libm never dispatches on `exp`" it is false, and
a version that widened that path would be waved through.

This was verified rather than assumed: `cbrt` hashed over 2.4 million inputs,
with and without the architecture-specific feature enabled, is bit-identical.

## Files

| File | What it holds |
|---|---|
| `lib.rs` | The unit newtypes and their arithmetic |
| `det_math.rs` | The one door for `exp`, `ln`, `sin`, `cos`, `powf` — and `acos` and `cbrt`, the latter being what `energy_per_unit` calls and what the `libm` note above rests on |

## Further reading

- Goldberg, *"What every computer scientist should know about floating-point
  arithmetic"*, ACM Computing Surveys 23(1), 1991 —
  [doi:10.1145/103162.103163](https://doi.org/10.1145/103162.103163). The
  standard reference for why decimals misbehave.
- IEEE 754-2019 — clause 5 for the required operations listed above, and clause
  9 for the elementary functions it merely *recommends*. The gap between those
  two clauses is the whole reason `det_math.rs` exists.
