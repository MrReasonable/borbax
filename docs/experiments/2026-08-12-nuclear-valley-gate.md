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

**Updated 2026-08-13: `sigma` was removed from the model after this gate
first passed — see "Post-gate correction, 2026-08-13" near the end of this
document for the re-measured numbers before trusting any `sigma`-referencing
figure elsewhere in this document as current.** Everything between here and
that section is left as the historical record of how the gate was
originally cleared (a four-coefficient model), not the model this crate
now ships.

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
   identity coefficients and the shipped `NUCLEAR_K = 10`: 305 of 386
   totals (79.0%) had negative on-valley per-unit binding energy — unbound
   over most of the domain. (An earlier draft of this correction quoted
   288/386, 74.6% — measured at the since-superseded `NUCLEAR_K = 12` and
   never re-taken; the direction is unaffected, 79% is worse — see
   `nuclear.rs`'s `the_identity_configuration_is_bound_over_all_but_two_totals`
   for the figure as a pinned test.) Corrected to `eps = 2·a_V/z`,
   `a_V = 15.75` MeV (Rohlf's real volume coefficient), `z = 12`
   (`crate::packing::shell_size(k, 1)` at `k = 10`) — the real limit
   `BE/A → a_V` as `A → ∞`, approached slowly and not attained within
   `T_MAX` (`eps · contacts(T_MAX)/T_MAX` is `13.19`, not `15.75`).

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

**Historical — `sigma` row removed 2026-08-13, see "Post-gate correction"
below.** This table reflects what the gate first ran against, not what
`crates/borbax-universe/src/nuclear.rs` ships today.

| constant | base | source |
|---|---|---|
| `eps` | `2.625` | `2·a_V/z`, Rohlf `a_V = 15.75` MeV, `z = 12` |
| `sigma` | `0.075` | V1's own historical range midpoint — **deleted 2026-08-13** |
| `kappa` | `23.7` MeV | Rohlf `a_sym` |
| `c` | `0.015` | Rohlf `a_C/(2·a_sym)` |
| `Rung::MAX` | `99` | `P(identity) = 1%` |

All coefficients drawn from `p ∈ [-0.5, 0.5]` (`Direction::P_MAX`, the
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
interior" framing. Ca-40's `40/386` is reported alongside as a real-world
comparison, not gated on.

**The replacement (the computed envelope itself) is a monotonicity/band
check, not a closed-form regression bound — corrected after
mutation-testing it.** An earlier draft of this section claimed the
envelope was itself "a real regression bound: a sign error, a dropped
factor of 2, or a wrong `t^(2/3)` spelling moves `drift_onset` outside
it" — false for at least one real defect class: a `c`-term sign error
prints the *identical* `[16, 30]` envelope with zero violations over a
20,001-point sweep, because the envelope is computed by calling the same
function the per-seed check then tests against — the same tautology shape
recurring one level deeper after replacing the Ca-40 anchor. (A
`cbrt(t)*cbrt(t)` respelling of `t^(2/3)` was also tried and is *not* a
second example — it produces zero integer differences from the pinned
spelling over the same sweep, so it is behaviourally inert and proves
nothing about the check's power either way.) What the corpus check
genuinely verifies: `drift_onset` is monotone non-increasing in `c`, and
every drawn `c` stays inside P7's own reachable band. **The tests that
actually carry closed-form correctness** are `the_drift_onset_at_identity_is_pinned`
(an exact regression pin, `== 22`) and
`the_closed_form_matches_the_brute_force_integer_argmax` (mismatches on
2,040-11,490 of 11,580 cells under every defect tried, including the two
above via the `gamma` factor-2 mutation).

**Result: `[16, 30]` — passes, 38,416/38,416 seeds.**

## Arm 2 — closure, bound-set shape, and coverage

Three checks, of which two are hard gates and one is reported only.

**2a. The plan's own literal criteria — hard gate.** Does the on-valley
per-unit binding energy peak strictly inside `1..T_MAX`, and does a full
search find a positive-Q size-shedding transition for `total = T_MAX` (the
"heavy on-valley element" proxy, since `n_elements` doesn't exist without
Element wiring)? Bar: `≥ 1 − 6261/59292 = 0.894404` — V1's own edge-peak
census, reconstructed *inside* the crate against the real `contacts_upto`
(9 × 9 × 12 × 61 = 59,292 cells over `k ∈ 6..=14`, `eps ∈ [0.8,1.2]`,
`sigma ∈ [0.02,0.13]`, `n_elements ∈ 60..=120`; V1's own `element.rs:887-914`
comment reports a compatible 2/24 ≈ 8.3% on a 24-seed sample).

**Result: 38,416/38,416 = 100%, well above the 89.44% bar.**

**2b. Bound-set shape — hard gate, added beyond the plan's own text.** A
`physics-plausibility-reviewer` pass found the natural stronger check
("binding positive everywhere") is simply false at the *identity*
configuration: `total = 3` has negative on-valley per-unit binding
(`-0.5615`) with no perturbation involved, matching a known limitation of a
four-term liquid-drop model with no pairing/curvature term. The real SEMF's
own prediction for the A=3 nuclide this model actually favours at that
total (Z=1, i.e. H-3 — both by the real SEMF's own argmax and by this
model's) is `+0.775` MeV, against H-3's real `2.827`, 3.65× off before any
seed is drawn. (An earlier draft of this comparison paired the model's Z=1
prediction against He-3's, Z=2, real value — the direction of the finding
is unaffected, but the two nuclides shouldn't be crossed.) What survives
contact with that fact: the bound/unbound sequence over `4..=T_MAX` should
have at most a handful of sign transitions (a clean single crossover, once
the model turns permanently unbound past some total, since
`sigma·total^(5/3)` is the model's only unbounded term).

**What this check does and does not discriminate — corrected after
mutation-testing it directly.** It reliably catches a *scattered* defect,
where several distinct regions of the range flip sign (a gross scale error
in `eps`, e.g., measured to produce 7 transitions at one sampled seed). It
does **not** reliably catch a *smooth* defect that still crosses zero once
— a wrong exponent on `sigma`'s term, or halving `sigma`'s coefficient,
both leave the transition count at 0 or 1 while moving where the single
crossing happens by a large amount, and one of the two (the wrong
exponent) makes the gate's own reported numbers look *better*, not worse.
That class is caught by two added pinned tests instead:
`the_identity_configuration_is_bound_over_all_but_two_totals` (the
identity's unbound set is exactly `{1, 3}`) and
~~`the_size_term_is_total_to_the_five_thirds`~~ (the term itself, plus the
identity peak position) — **deleted 2026-08-13 along with `sigma` itself
(see "Post-gate correction" below); its job is now split across the three
term-isolation tests in `nuclear.rs`** (`the_volume_term_is_eps_times_contacts`,
`the_asymmetry_term_is_kappa_times_imbalance_squared_over_total`,
`the_coulomb_term_is_gamma_times_a_times_a_minus_one_over_total_to_the_one_third`).

Measured over the corpus: **p50 = 0 transitions (never turns unbound at
all), p90 = 1 (the ordinary single crossover), max = 3** — three seeds
(14355, 20107, 21217), each a one-total re-binding blip at `total = 309`,
which is `crate::packing::shell_size(10, ·)`'s **exact** cumulative
shell-4 boundary (`1+12+42+92+162 = 309`) — a completed shell buying one
extra total of binding, an "island of stability" shape analogous to (not
the same mechanism as) real nuclear shell closures — real magic numbers
come from the spin-orbit-split shell model, not geometric packing, and
this crate's own closures don't coincide with them. Three independent
seeds landing on the same boundary is what makes it the mechanism rather
than coincidence. **Updated 2026-08-13**: the gate was originally set at
`≤ 5` transitions, "clearing the measured worst case with margin"; a
follow-up specialist cross-check (`geometry-numerics-reviewer`, after the
fourth `/review-pr` round on PR #38) found `5` was not a margin but an
off-by-one line drawn inside the nearest
out-of-box defect regime, and proved `3` is this mechanism's own
structural ceiling (one permanent crossing plus at most one closure-blip,
which costs exactly 2) — see `bound_shape`'s own doc in
`crates/borbax-universe/src/nuclear.rs` for the parity proof that `4`
would have been an empty margin. Gate now set at `≤ 3`, with zero verdicts
changed on the corpus below.

**Result: 38,416/38,416 pass (`total = 3` excluded from the count, as
designed).**

**2c. Elements reached while still bound — reported, not gated.** How many
proton counts does the valley reach before the bound-set's permanent
unbound tail begins? A `physics-plausibility-reviewer` pass found
`perturbation.rs`'s own `MigratedConstant::C` doc overclaimed here: its
coverage proof (that `T_MAX` is large enough for `composition_argmax` to
reach every integer up to 120) covers only the composition argmax's
*shape*, not whether binding stays positive that far out — a different
question once `eps`/`sigma` move too. (All four nuclear coefficients — now
three, `sigma` deleted 2026-08-13, see "Post-gate correction" below —
shared one rung, so their *magnitudes* moved together, not independently
of `c` — only each constant's *direction* is drawn independently, from its
own `Stream` index.)

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

**`sigma`, not `eps`, is the larger driver of this shortfall — a
`/review-pr` finding, also routed to Steps 2+, not resolved here.**
`eps · contacts(total)/total` reproduces the real SEMF's volume-plus-surface
term with worst-case error 2.63% (at `total = 78`), mean absolute error
1.20%, holding within 1% only from `total ≥ 266`. This model's
asymmetry/Coulomb terms are already exactly the real SEMF's own
(`gamma = 0.711 = a_C` at identity), so at identity this model is close to
the real four-term SEMF plus one extra term with no real counterpart and
the wrong exponent sign relative to the real surface term. Sweeping
`sigma` alone at the identity `kappa`/`c`: zero-crossing moves from
`total = 3055` (`sigma = 0`, close to the real `3076` under this model's
own `Z(Z-1)` Coulomb convention) to `603` at the shipped `BASE_SIGMA =
0.075`; the peak stays at `total = 55` regardless. Sharper still: **at
`BASE_SIGMA = 0.0`, every test in this module still passes and the
coverage shortfall improves from 9.25% to 1.09%** — `sigma` is not merely
unneeded for the interior peak at the pre-registered value, it actively
costs coverage for a term with no physical counterpart, within the range
this model actually samples. Whether `BASE_SIGMA` should be reduced is a
physics-plausibility question for Steps 2+.

## V1 census cross-check (in-crate, settles a transcription caveat)

The 59,292-cell V1 edge-peak census above was independently computed twice
— once by `geometry-numerics-reviewer` in a standalone transcription of
`packing.rs`, once here as a permanent `#[cfg(test)]` inside
`borbax-universe` against the real `contacts_upto`/`shell_size`. Both give
**6,261/59,292 = 10.559603%** exactly. At the identity coefficients, the
on-valley peak lands at `total = 55, element 25` — this exactly matches
the physics-plausibility pass's own scratch-script prediction for `k=10,
eps=2.625` before any of this was production code, a genuine cross-check
of the script against the implementation.

**Post-gate correction, 2026-08-13: `total=55`'s closeness to a real SEMF
peak (A=60, or the experimental A~61-62) is not itself independent
corroboration of the physics — a `/review-pr` finding (`emergence-auditor`).**
`valley_peak_total` returns exactly one of four values across the entire
38,416-seed corpus: `{13, 55, 147, 309}` — this crate's own packing-shell
closures at `NUCLEAR_K=10` (`1+12=13`, `+42=55`, `+92=147`, `+162=309`),
never anything else. The model's own *smooth* counterpart (no shell
quantisation, same asymmetry/Coulomb terms, the volume-plus-surface piece
alone) peaks at `total=70` (using this doc's own fitted `a_S=19.75`) or
`total=60` (Rohlf's real `a_S=17.8`) — the quantisation pulls the reported
peak down from a real term-competition optimum near 60-70 to the nearest
closure below it, 55. `55`'s proximity to `60` is a product of the lattice
happening to place a closure nearby, not of the coefficients balancing to
match real physics — a coefficient error large enough to move the smooth
optimum from 60 to, say, 50 would very likely still round to the same
closure, 55, and leave this "cross-check" reporting success. See Arm 2's
own doc, below, for the direct consequence.

## Post-gate correction, 2026-08-13: `sigma` removed

Found after the gate above already returned "both arms pass" — a
`/review-pr` follow-up question, not a pre-verification finding like
corrections 1-3. Full physics account:
`crates/borbax-universe/src/nuclear.rs`'s own module doc, "Post-gate
correction." Summary: a `physics-plausibility-reviewer` pass found
`sigma * total^(5/3)` structurally unjustified, not merely miscalibrated —
`eps * contacts(total)` already reproduces the real SEMF's volume-and-surface
behaviour, so the three surviving coefficients (`eps`, `kappa`, `c`) already
form the real four-term SEMF at Rohlf's own values. The "growth limit"
argument that had justified keeping `sigma` (§ Arm 2, 2b above, "since
`sigma·total^(5/3)` is the model's only unbounded term") does not hold: a
bounded cost exceeding a saturating gain limits growth without needing to
be unbounded itself, and `eps`/`kappa`'s own saturation already does that.
`BASE_SIGMA` and `MigratedConstant::Sigma` are deleted; `MigratedConstant`'s
index 14 is retired, not reused, so `kappa`/`c`'s draws are bit-identical
to before.

**Both arms re-run over the same 38,416-seed corpus, three-coefficient
model:**

- **Arm 1**: `[16, 30]` — unchanged, as expected (`sigma` never entered
  `composition_argmax`). Still 38,416/38,416.
- **Arm 2a** (interior peak + positive-Q): **38,406/38,416 = 99.974%**,
  down from 100% — `sigma` was, per the physics reviewer's own finding, a
  meaningful contributor to positive shedding-Q at `T_MAX` for a small
  number of universes, and losing it costs ten of them. Still comfortably
  above the 89.44% bar.
- **Arm 2b** (bound-set shape): still 38,416/38,416, but the shape
  underneath changed. Last-bound-total range narrowed from `[62, 386]` to
  `[92, 386]` (the worst universe now stays bound 30 totals further in);
  **98.6% of universes (37,863/38,416) are now never unbound
  anywhere in `4..=T_MAX` at all** (most of what this check used to
  discriminate against no longer occurs at all); the transition-count
  distribution moved from **p50 = 0, p90 = 1, max = 3** to **p50 = 0,
  p90 = 0, max = 3** — max unchanged at the top. **Updated 2026-08-13**:
  the bound is `≤ 3`, not `≤ 5` (see the note above this section) — `3`
  clears this measured max with zero margin remaining above it, which is
  correct: `3` is the mechanism's structural ceiling, not a margin below
  one. The `total = 309` shell-closure blip still occurs —
  two seeds (3992, 4425) hit it in this corpus, down from three (14355,
  20107, 21217), for the same reason: fewer universes reach an unbound
  region at all for the blip to interrupt.
- **Arm 2c** (elements reached while bound, reported not gated): **min = 38,
  p1 = 116, p10 = 130, p50 = 138; 419/38,416 (1.09%) of universes cover
  fewer than 120 elements while still bound**, down from 9.25%
  (3,554/38,416) — matching exactly the figure the `emergence-auditor`
  pass measured at `BASE_SIGMA = 0.0` on the four-coefficient model before
  the term was actually deleted, as it should.

**Not the justification for the removal, only its consequence** — see
`nuclear.rs`'s own module doc for why the coverage improvement is
deliberately not the stated reason.

## What Step 1 does not settle

- **Whether the residual 1.09%-of-universes coverage shortfall (2c, above,
  post-removal) needs a design response** — routed to Steps 2+, which
  decide how `n_elements`'s own draw interacts with it. Smaller than the
  pre-removal 9.25%, not zero.
- **The size-shedding channel's Q-value** used by 2a is a Step-1-only
  definition (`crate::nuclear::shedding_q` — does splitting a heavy total
  into two on-valley fragments release energy), not the plan's own
  per-isotope definition, which Steps 5-7 still own.
- **Whether `sigma`/`gamma` were collinear on the valley** was a
  pre-registered caveat before this gate first ran, and is now moot —
  `sigma` no longer exists.
- ~~The `NuclearCoeffs` struct... remains deferred to Steps 2+~~ — **landed
  2026-08-13**, once the `sigma`-removal diff itself made the churn
  concrete rather than speculative (7 of 9 changed function signatures had
  zero logic change). See `crates/borbax-universe/src/nuclear.rs`'s own
  doc on the struct for why `composition_argmax`/`composition_argmax_int`/
  `drift_onset` deliberately still take `c: f64` alone. Whether Steps 2+
  want a different shape at their own real call sites is still open.
