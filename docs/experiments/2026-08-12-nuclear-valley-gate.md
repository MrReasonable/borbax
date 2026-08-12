# Task 26.2 Step 1 — the two-armed empirical gate

Issue #26. Full mechanism, pre-registered coefficients and routed
requirements: `docs/superpowers/plans/2026-08-07-borbax-real-physics.md`,
Task 26.2 section and Decision 11. Implementation and the permanent
regression tests: `crates/borbax-universe/src/nuclear.rs`. Run the corpus
sweeps yourself:

```bash
cargo test --locked --release -p borbax-universe --lib nuclear:: -- --ignored --nocapture
```

**Verdict: both arms pass. Step 1 is cleared to proceed to Steps 2+.**

## Three corrections landed before the gate could run meaningfully

A three-amigos design session proposed the initial harness; two independent
verification passes (`geometry-numerics-reviewer`, a
`physics-plausibility-reviewer` dispatch) each found real errors in it
before any of it became code — this project's now-standard pattern of not
trusting a single derivation, however careful, on load-bearing physics math.

1. **`composition_argmax` ships the κ-free closed-form spelling**
   (`(2t + c·t^(2/3)) / (4 + 2c·t^(2/3))`), not the one Decision 11's text
   calls "the correct closed form." The two are algebraically identical
   (both scale by `2κ`) but not numerically identical in `f64`: over
   6,951,474 sampled `(κ, c, t)` cells they differ by up to 4 ulp, and in 3
   of them (`c ≈ 0.0192, t = 125`) round to a **different integer** — a
   different element. The κ-free spelling has no `κ` in it at all, so a
   perturbation of `κ` provably cannot move any isotope assignment.

2. **`BASE_EPS = 2.625`, not the pre-registration commit's original `1.0`.**
   `kappa`/`c` are pinned to Rohlf's real MeV-scale SEMF coefficients; `eps`
   was left at V1's old, unrelated Borbax-internal scale. Measured at the
   identity coefficients: 288 of 386 totals (74.6%) had negative on-valley
   per-unit binding energy — unbound over most of the domain. Corrected to
   `eps = 2·a_V/z`, `a_V = 15.75` MeV (Rohlf's real volume coefficient),
   `z = 12` (`crate::packing::shell_size(k, 1)` at `k = 10`) — the real
   limit `BE/A → a_V` as `A → ∞`.

3. **`NUCLEAR_K = 10`, not `12`.** `contacts(total)` has no source under the
   V2 electronic-configuration model (`generate_elements_v2` never calls
   `packing` at all), so Task 26.2 needs its own fixed coordination
   parameter. The FCC/icosahedral coordination-12 argument that motivated
   fixing it at all actually selects `k = 10`: this crate's shell law gives
   `z = k + 2`, so `k = 12` (the original proposal) gives coordination 14,
   not 12.

Both corrections 1 and 2 are recorded directly in `perturbation.rs`'s
`MigratedConstant::Eps`/`C` doc comments, not only here.

## Coefficients as actually run

| constant | base | source |
|---|---|---|
| `eps` | `2.625` | `2·a_V/z`, Rohlf `a_V = 15.75` MeV, `z = 12` |
| `sigma` | `0.075` | V1's own historical range midpoint |
| `kappa` | `23.7` MeV | Rohlf `a_sym` |
| `c` | `0.015` | Rohlf `a_C/(2·a_sym)` |
| `Rung::MAX` | `99` | `P(identity) = 1%` |

All four coefficients drawn from `p ∈ [-0.5, 0.5]` (`Direction::P_MAX`, the
crate-wide default), jointly through one shared rung — `Rung::draw(seed)`
once, then `draw_symmetric` per constant — never independently, per routed
requirements 5/6.

## Arm 1 — drift onset

**Gate: `drift_onset(c)` (smallest even `total` where the composition
argmax first departs from the symmetric split) lies in the analytic
envelope `[drift_onset(1.5·base_c), drift_onset(0.5·base_c)] = [16, 30]`,
for every sampled universe.**

The plan's own text proposed gating against a fraction of a real-world
anchor (Ca-40, the heaviest observationally-stable N=Z nuclide, `40/386`).
Measured: `kappa` provably cancels from the closed form, and `c` is bounded
to `[0.0075, 0.0225]` by P7's mechanism before any sampling happens — so
`drift_onset`'s reachable range is exactly `[16, 30]` for *every* legal
seed, deterministically. A bar of `40/386` can never fail regardless of
implementation correctness — the same tautology shape the plan's own
round-4 correction already fixed once for Arm 1's original "is the optimum
interior" framing. Replaced with the analytic envelope itself: a real
regression bound (a sign error, a dropped factor of 2, or a wrong
`t^(2/3)` spelling moves `drift_onset` outside it), derived rather than
guessed. Ca-40's `40/386` is reported alongside as a real-world comparison,
not gated on.

**Result: `[16, 30]` — passes, 38,416/38,416 seeds.**

## Arm 2 — closure, bound-set shape, and coverage

Three checks, of which two are hard gates and one is reported only.

**2a. The plan's own literal criteria — hard gate.** Does the on-valley
per-unit binding energy peak strictly inside `1..T_MAX`, and does a full
search find a positive-Q size-shedding transition for `total = T_MAX` (the
"heavy on-valley element" proxy, since `n_elements` doesn't exist without
Element wiring)? Bar: `≥ 1 − 6261/59292 = 0.894404` — V1's own edge-peak
census, reconstructed *inside* the crate against the real `contacts_upto`
(9 × 9 × 12 × 61 = 59,292 cells over `k ∈ 6..14`, `eps ∈ [0.8,1.2]`,
`sigma ∈ [0.02,0.13]`, `n_elements ∈ 60..120`; V1's own `element.rs:887-914`
comment reports a compatible 2/24 ≈ 8.3% on a 24-seed sample).

**Result: 38,416/38,416 = 100%, well above the 89.44% bar.**

**2b. Bound-set shape — hard gate, added beyond the plan's own text.** A
`physics-plausibility-reviewer` pass found the natural stronger check
("binding positive everywhere") is simply false at the *identity*
configuration: `total = 3` has negative on-valley per-unit binding
(`-0.5615`) with no perturbation involved, matching a known limitation of a
four-term liquid-drop model with no pairing/curvature term (the real SEMF's
own A=3 prediction is `+0.77` MeV against He-3's real `2.57` — already 3.3×
off before any seed is drawn). What survives contact with that fact: the
bound/unbound sequence over `4..=T_MAX` should have at most a handful of
sign transitions (a clean single crossover, once the model turns
permanently unbound past some total, since `sigma·total^(5/3)` is the
model's only unbounded term). A genuinely broken model (sign error, wrong
exponent) scatters negativity across a large fraction of the range — many
transitions, not a handful.

Measured over the corpus: **p50 = 0 transitions (never turns unbound at
all), p90 = 1 (the ordinary single crossover), max = 3** — the one outlier
(seed 14355) is a one-total re-binding blip at `total = 309`, which is
`crate::packing::shell_size(10, ·)`'s **exact** cumulative shell-4 boundary
(`1+12+42+92+162 = 309`) — a completed shell buying one extra total of
binding, the same "island of stability" shape shell closures produce in
real nuclei. Gate set at `≤ 5` transitions, clearing the measured worst
case with margin.

**Result: 38,416/38,416 pass (`total = 3` excluded from the count, as
designed).**

**2c. Elements reached while still bound — reported, not gated.** How many
proton counts does the valley reach before the bound-set's permanent
unbound tail begins? A `physics-plausibility-reviewer` pass found
`perturbation.rs`'s own `MigratedConstant::C` doc overclaimed here: its
coverage proof (that `T_MAX` is large enough for `composition_argmax` to
reach every integer up to 120) covers only the composition argmax's
*shape*, not whether binding stays positive that far out — a different
question once `eps`/`sigma` are perturbed independently of `c`.

Measured: **min = 27, p1 = 62, p10 = 122, p50 = 138; 3,554/38,416 (9.25%)
of universes cover fewer than 120 elements while still bound.** This is a
real, separate finding, corrected directly in `perturbation.rs`'s `C` doc
comment. It is *not* gated here because (a) the plan's own literal Arm 2
text asks only about interiority and a positive-Q transition, both of
which pass cleanly, and (b) fixing it would mean redesigning `eps`'s
excursion bound or coupling it to `kappa` — considered and rejected (real
SEMF volume/asymmetry coefficients are independently fitted, with no
derivable relationship: `a_A/a_V` spans 1.348–1.505 across five published
fits against `a_V`'s own ±5.7% spread) — a decision for Steps 2+, which own
how this bears on `n_elements`'s own draw, not something Step 1 is
positioned to resolve by picking a threshold.

## V1 census cross-check (in-crate, settles a transcription caveat)

The 59,292-cell V1 edge-peak census above was independently computed twice
— once by `geometry-numerics-reviewer` in a standalone transcription of
`packing.rs`, once here as a permanent `#[cfg(test)]` inside
`borbax-universe` against the real `contacts_upto`/`shell_size`. Both give
**6,261/59,292 = 10.559603%** exactly. At the identity coefficients, the
on-valley peak lands at `total = 55, element 25` — a strong independent
cross-check, since this exactly matches the physics-plausibility pass's own
scratch-script prediction for `k=10, eps=2.625` before any of this was
production code.

## What Step 1 does not settle

- **Whether the 9.25%-of-universes coverage shortfall (2c, above) needs a
  design response** — routed to Steps 2+, which decide how `n_elements`'s
  own draw interacts with it.
- **The size-shedding channel's Q-value** used by 2a is a Step-1-only
  definition (`crate::nuclear::shedding_q` — does splitting a heavy total
  into two on-valley fragments release energy), not the plan's own
  per-isotope definition, which Steps 5-7 still own.
- **Whether `sigma`/`gamma` are collinear on the valley** was pre-registered
  as a caveat before this gate ran (`gamma/4 = 0.17775` exceeds V1's own
  `sigma` ceiling of `0.13`) and stands unchanged by this measurement.
