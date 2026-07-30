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

Every reaction must conserve mass **exactly**. Decimals cannot promise that:
add 0.1 and 0.2 in floating point and you get 0.30000000000000004. Over
millions of reactions those crumbs accumulate and mass quietly appears or
vanishes.

Fixed-point integers have no crumbs. `Mass` stores a whole number of tiny
units, so `a + b - b == a` always, and conservation is exact by construction
rather than by tolerance.

## Job 2: maths that gives the same answer everywhere

Here is a fact that surprises most people: **`exp`, `ln`, `sin`, `cos` and
`powf` are not fully specified by the floating-point standard.** IEEE-754 pins
down `+`, `-`, `*`, `/` and `sqrt` to the last bit, and then stops. The
transcendental functions are left to whoever writes the maths library.

So `(2.5_f64).exp()` can genuinely produce a *different final bit* on macOS than
on Linux. One bit, propagated through a few million simulation steps, is a
different universe.

### The fix: one door

Every transcendental in the entire workspace goes through `det_math.rs`, which
calls [`libm`](https://crates.io/crates/libm) — a pure-Rust implementation with
no platform-specific paths for these functions. Same code, same answer,
everywhere.

A direct `.exp()` or `.cos()` on an `f64` anywhere in library code is a
**determinism bug**, not a style preference. `clippy.toml` bans them by name so
the build fails rather than the goldens.

`sqrt` and the four arithmetic operations stay native — those *are* exactly
specified, so routing them through a library would cost speed and buy nothing.

### The subtlety worth knowing

`libm` **does** dispatch on processor architecture — it is not a myth, and
"it doesn't" is the wrong reason to trust it. What makes it safe is *which*
operations it dispatches on: only ones IEEE-754 specifies exactly (`sqrt`,
`fma`, `rint`, `ceil`, `floor`). Those give identical results whichever path
runs.

This was verified rather than assumed: `cbrt` hashed over 2.4 million inputs,
with and without the architecture-specific feature enabled, is bit-identical.

## Files

| File | What it holds |
|---|---|
| `lib.rs` | The unit newtypes and their arithmetic |
| `det_math.rs` | The one door for `exp`, `ln`, `sin`, `cos`, `powf` |

## Further reading

- Goldberg, *"What every computer scientist should know about floating-point
  arithmetic"*, ACM Computing Surveys 23(1), 1991 —
  [doi:10.1145/103162.103163](https://doi.org/10.1145/103162.103163). The
  standard reference for why decimals misbehave.
- IEEE 754-2019, which specifies the five exact operations and, pointedly, not
  the others.
