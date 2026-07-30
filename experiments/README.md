# `experiments` — measuring things before committing to them

This crate is not part of the simulation. It exists so design questions get
answered with **numbers** instead of arguments.

The pattern throughout Borbax is: someone proposes a mechanism, someone else
says it will not work, and rather than debating it we build the smallest thing
that produces a measurement. Several of the project's load-bearing decisions
were settled here, and at least two proposals died on contact with their own
data.

## What lives here

| File | The question it answers |
|---|---|
| `rng.rs` | A small standalone generator, so experiments do not depend on the real one |
| `molecule.rs` | Toy molecules — build, mutate, canonicalise, hash |
| `signature.rs` | Turning a shape into a comparable descriptor |
| `geodesic.rs` | Points spread evenly over a sphere |
| `embed.rs` | Laying a graph out in 3D |
| `g2.rs` | **Locality**: does a small change to a molecule produce a small change in its shape? |
| `openended.rs` | Does the system keep producing genuinely new things, or plateau? |
| `bin/ptable.rs` | Periodic-table shape under different drawn constants |
| `bin/fusion.rs` | Cluster fusion and the binding peak |
| `bin/g2.rs` | Runs the locality measurement across several regimes |

## The two ideas worth understanding

### 1. Locality (`g2.rs`)

For evolution to work at all, **small changes must usually have small effects.**
If swapping one atom completely scrambles a molecule's shape, then every
mutation is fatal, nothing can be gradually improved, and there is no fitness
landscape to climb — just noise.

So `g2` measures exactly that: mutate a molecule slightly, measure how far its
shape descriptor moved, and check the distribution. This is a **guarantee**
(spec §2.5), not a nice-to-have.

### 2. The neutral shadow (`openended.rs`)

Suppose a run produces a steadily rising count of new molecule types. Is that
evolution, or is it just what *any* random process does?

The only way to know is to run a **control** — a shadow simulation with the same
machinery but no selection — and compare. A rising curve that the shadow also
produces tells you nothing.

This is why the neutral-shadow seed is a separate constant from the main one:
if the control accidentally replayed the treatment's random draws, it would stop
being a control.

Borbax's design doctrine (spec §2.7) refuses two claims on exactly these
grounds: never assume self-maintaining units compose into bigger ones, and never
declare a plateau by eyeballing a flattening curve — fit competing models and
compare.

## The random generator here

`rng.rs` implements **SplitMix64** — Steele, Lea & Flood, *"Fast splittable
pseudorandom number generators"*, OOPSLA 2014,
[doi:10.1145/2660193.2660195](https://doi.org/10.1145/2660193.2660195).

It works in two steps:

1. **Advance.** Add a fixed odd number (`GAMMA`, which is exactly
   `floor(2⁶⁴ ÷ φ)`) to the state. Because it is odd, repeatedly adding it visits
   every one of the 2⁶⁴ possible states before repeating — oddness is both
   necessary and sufficient here, since the additive order of `g` in `Z/2⁶⁴` is
   `2⁶⁴ / gcd(g, 2⁶⁴)`.

   **The golden ratio itself is doing nothing for the period.** *Any* odd
   constant gives the full 2⁶⁴. What it buys is that the raw, unmixed state
   sequence is low-discrepancy — a convention inherited from Fibonacci hashing.
   Since every value leaves through the mixer in step 2 before anyone sees it,
   that property is not what makes the output usable; step 2 is.
2. **Mix.** Scramble that state hard, using shift–XOR–multiply three times over.

Step 2 is what earns its keep. Without it, seeds 0, 1 and 2 would produce
nearly-identical first draws — and this harness derives *one stream per trial
index*, so correlated adjacent seeds would silently correlate the trials. The
test `adjacent_seeds_decorrelate_immediately` exists for precisely that, and
requires every pair of first draws to differ in at least 16 of their 64 bits.

The mixing constants are named (`MIX_MULTIPLIER_1`, `MIX_SHIFT_1`, …) rather
than inline, but naming them does not make them adjustable — they were found by
search to maximise avalanche, and round-looking replacements break the
decorrelation while leaving most tests green.

### Picking a number below `n`

The obvious `random % n` is **biased**. If `n` does not divide 2⁶⁴ exactly, the
low values come up slightly more often — and for a harness picking atoms and
edit sites, that quietly skews which molecules get generated.

`below()` uses **Lemire's multiply-shift with rejection** instead, which is
exactly uniform. The rejection branch fires with probability `(2⁶⁴ mod n)/2⁶⁴`,
which is under `n / 2⁶⁴` — and `n` here is an atom count or an edit-site count,
so in practice never. (The bound is only reassuring because `n` is small: at
`n ≈ 2⁶³` it would be one rejection in two.)

### Unit floats

`f64_unit()` takes the top 53 bits of a draw and divides by 2⁵³. 53 is exactly
the number of significant bits an `f64` holds, so the conversion is **exact** —
nothing is rounded, and every representable value in `[0, 1)` at that spacing is
reachable.

## A note on seeds

Every experiment names its seeds as constants. The values are arbitrary — the
requirements are only that they are *fixed* (so a measurement reproduces) and
*unrelated to each other* (so two regimes cannot share a stream prefix).

They are written as recognisable bit patterns rather than 1, 2, 3 so that a
value turning up somewhere unexpected is obviously a seed and not a count.

## Running them

```bash
cargo run --locked -p borbax-experiments --release --bin g2
cargo run --locked -p borbax-experiments --release --bin ptable
cargo run --locked -p borbax-experiments --release --bin fusion
```

Release matters — these do millions of trials and the debug profile is roughly
an order of magnitude slower.
