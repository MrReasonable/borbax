# Review synthesis, cross-checked — Borbax V0 plan (2026-07-26, round 2b)

Supersedes the disputed and cross-check-failed sections of
[`synthesis.md`](synthesis.md). That document was produced **without the `Agent`
tool**, so its cross-check round never ran and five findings were marked
`cross-check-failed`. This run dispatched them.

Everything in `synthesis.md` not touched below **still stands** — this is a
delta, not a replacement. Read that file for the full finding list.

---

## 0. What the cross-checks changed

Eight cross-checks dispatched, eight returned. **No finding was retracted by
its originator, so nothing is dropped.** One *endorsement* was retracted, which
settles the review's sharpest dispute. And **six proposed fixes were rejected
while every one of their findings was confirmed** — the distinction the
reconciliation table does not draw, and by a distance the most valuable thing
this round produced.

| # | Finding | Cross-checked with | Verdict | Effect |
|---|---|---|---|---|
| G1 | `ideal_gap` on the wrong length scale | emergence-auditor | **`confirmed`** | `[single-reviewer]` → **`[corroborated]`**; fix amended |
| G2 | `canonicalise_frame` worse than nothing | rust-developer-expert | **`confirmed`** | `[disputed]` → **`[corroborated]`**; scope grew |
| G3 | boustrophedon fold start retraces | rust-performance-expert | **`disputed`** | finding **`[corroborated]`**, **fix superseded** |
| E2 | contact energy still inverted | geometry-numerics-reviewer | **`disputed`** | finding **`[corroborated]`**, **fix superseded** |
| alr-011 | volume `V` is not a parameter | rust-performance-expert | **`confirmed`** | `[single-reviewer]` → **`[corroborated]`**; *my* reasoning falsified |
| D-1 | `lining_signature` — 1-vs-1 dispute | emergence-auditor | **`retracted`** | **dispute closed, geometry right** |
| P5 | Bedau activity form | alife-researcher | **`disputed`** | **criterion 9 would pass on a dead system** |
| bench 1 | pairwise summation in the affinity kernel | determinism-auditor | **`disputed`** | **do not adopt** — decided in the window |
| P2 | flood replacement claimed bit-identical | determinism-auditor | **`confirmed`** | right answer, wrong reason; **M-c downgraded** |

**Four things the cross-checks changed that a reader of `synthesis.md` would
otherwise get wrong:**

1. **D-1 is closed.** The emergence-auditor **withdrew its endorsement in
   full** after measuring. `synthesis.md` presented this as the review's
   hardest open dispute needing a human. It is now settled arithmetic.
2. **Four proposed fixes are wrong, and all four findings are right.** G3's
   greedy walk, G1's `gap_factor`, E2's `((fa+fb)²−4)·w` and P5's
   `∫ count_i dt` were each rejected by the cross-check that confirmed the
   underlying defect. **A third fix pass applying `synthesis.md` as written
   would have implemented all four** — and in three cases the originator's own
   proposed test would have passed while the defect remained. This is the
   round's central result; see §6.
3. **My own `[coordinator]` reasoning on `alr-011` was falsified.** I asserted
   the Gillespie rate constant is resolved once per channel. It is not.
   Recorded in §4 as a retraction.
4. **Four findings got *worse* under measurement, not better** — D-1, D-5, G2
   and E2. In every case the reviewer who measured had understated the
   severity, and in D-1's case by a wide margin.
5. **Exit criterion 9 would pass on a system with no open-ended novelty.**
   New, from the alife cross-check that `synthesis.md` flagged as the cheapest
   unclosed item in the review. A defect in the acceptance instrument outranks
   any defect in a task.
6. **One finding was made *less* severe: M-c, the cavity-flood torus wrap.**
   It fires 64,512 times on a full-span box and is **inert** — 0 events add a
   cell that was not already reachable. Still worth fixing, but for a
   structural reason, not an observable one (§2.8). This is the only downgrade
   in the round, and it matters that there is one: it is evidence the
   cross-checks were not simply ratifying.

---

## 1. Counts, revised

| Classification | Before | After | Note |
|---|---|---|---|
| `corroborated` | 16 | **22** | +G1, +G2, +G3(finding), +alr-011, +D-1, +D-5 |
| `disputed` (reviewer vs reviewer) | 6 | **0** | all six resolved; see §3 |
| `solo`, cross-check-failed | ~28 | **~23** | the five named are now resolved |
| **Fix proposals rejected or amended** | 0 | **6 of 7** | G1 `gap_factor`, G3 greedy walk, E2 form, P5 form, D-1 cone-vs-support, P2's justification |
| Findings **downgraded** in severity | 0 | **1** | M-c, the torus wrap — measured inert |
| Reviewer self-retraction | 1 | **2** | +emergence's D-1 endorsement |
| Coordinator retraction | 0 | **1** | my `alr-011` cost reasoning |
| Findings **retracted** by their originator | 1 | **1** | unchanged — nothing new was dropped |

**Contradiction density: 6 → 0.** Every reviewer-vs-reviewer dispute in
`synthesis.md` is resolved: D-1 by retraction (§2.6), D-2 and D-5 by
measurement (§2.9), D-4 by concession (§2.2), D-6 by G1's cross-check
(§3), and D-3 by a coordinator ruling I have since strengthened (§3).

That is well below the Step-5 escalation threshold, and it confirms
`synthesis.md`'s reading: the apparent incoherence was a **method** difference,
not an inconsistent artifact. Four of the six were the same argument —
reviewers who measured against reviewers who read — and **when the readers
measured, they agreed, unanimously.**

What remains in §3 is two *decisions* with no dispute attached, which is a
different and much cheaper thing.

---

## 2. The cross-checks, in detail

### 2.1 · G1 — `ideal_gap` — **`confirmed`**, and the fix is amended

**Was:** `[single-reviewer]` `[measured]` `cross-check-failed`.
**Now:** `[corroborated]` — geometry G1 + emergence-auditor, both measured.

The emergence-auditor reproduced both numbers independently and supplied the
analytic reason geometry had only measured. Writing
`x_i = r_A[i] + r_B[ANTI[PERM[R][i]]]`:

```
score(R) = −w_shape·D·[ Var_R(x) + (m − g)² ] − w_charge·(…)
```

where `m = mean(r_A) + mean(r_B)` is **exactly rotation-invariant** because
`anti∘perm` is a bijection — verified to 1.7e-13 over the real 60-element group.
So the gap enters as a **pure additive offset** `D(m−g)²` that the rotation
search cannot touch: 420 against a variance part of 21. Geometry measured
0.069/1.066 contrast on the real pipeline; emergence got 0.042/1.068 on
synthetic molecules; my own crude reproduction gave 0.234/0.713. Three
reproductions, same direction, magnitudes differing with the molecule model.

**On the emergence question — is a per-pair gap laundering?** Verdict:
**no.**

> *"Laundering needs a conditional — a branch keyed on a classification, or a
> parameter conditioned on the outcome it is meant to produce. `pair_gap` has
> neither… Deletion test: remove it and binding still happens — badly, ranked
> by size. It is a reparameterisation of the one mechanism, not a second path
> to a phenomenon."*

And the stronger form: `g = m` is not a choice, it is *the unique value that
zeroes the size component of an orthogonal decomposition the score already
has*.

**AMENDMENT — `gap_factor` is rejected.** Geometry proposed
`k.gap_factor * (ma + mb)` with `gap_factor` a new generated per-universe
constant. Emergence rejects it:

> *"The offset is `D·m²(1−f)²` — `gap_factor` is a knob whose only effect is
> how hard to penalise large molecules, quadratically in size. That is
> precisely the knob that gets swept until something appears… Take `f = 1`
> structurally, add no constant. Removing a generated constant beats replacing
> one."*

Adopt `pair_gap = mean(a.r) + mean(b.r)`, with no factor.

**CONSEQUENCE — `may_bind` becomes exactly vacuous. [coordinator-verified]**
At `f = 1`, `mean_defect = (Σa + Σb − D·gap)/D ≡ 0`, so
`ceiling = −w_shape·D·0² = 0`, and `0 >= threshold` for every negative
threshold. I verified over 2,000 random pairs: max |ceiling| = **2.7e-29**.
`may_bind` accepts 100% of pairs.

This **partially retracts a claim inside G1 itself** — geometry wrote *"symmetry
and the `may_bind` bound both survive unchanged"*. The bound survives as
*valid*; it does not survive as a *filter*. Note the irony, which is the
finding: `may_bind` was rewritten this round precisely because the previous
version *"rejected 0.000% of 20,000 random pairs: it filtered nothing while
being unprovable"* (`v0.md:4296`). Applying G1's fix returns it to filtering
nothing — provably this time. A replacement bound must come from **sorted-extent
variance**, not the mean.

**Also open, raised by emergence, unresolved:** whether the pair gap should use
`mean` or the closest-approach `max_i x_i`. `max` is the physically correct
separation *and is rotation-dependent*, so it would participate in the search
rather than factoring out of it. That is a geometry question and nobody has
answered it.

**Do not carry thresholds over.** Score magnitudes drop by roughly an order of
magnitude or more; `bind_midpoint` and every affinity threshold need
re-deriving, not rescaling.

---

### 2.2 · G2 — `canonicalise_frame` — **`confirmed`**, and Task 8 does not land either way

**Was:** `[disputed]` (D-4), `cross-check-failed`.
**Now:** `[corroborated]` — geometry measured it; rust-developer-expert agrees
it should be **removed, not repaired**.

rust-dev gave exactly the answer that was most useful, in the words the prompt
asked for:

> *"My repair is correct **conditional on the function surviving**… do not read
> my sign-off as evidence the function is sound. It was a repair to one branch
> of a function whose whole premise geometry has since falsified."*

That closes the specific risk D-4 existed to flag — a reader seeing four
sign-offs and concluding the function is fine.

**NEW, and it enlarges the finding: removing it does not make Task 8 landable.**

> *"At 54.5% over the bound, `small_graph_changes_give_small_geometry_changes`
> fails without the frame too. The plan's Step 4 instruction ('do not relax the
> threshold… raise `ITERATIONS`') must be rewritten in the same edit, or the
> next executor will relax it in frustration."*

This is the point `synthesis.md` had as a recommendation and now has as a
consequence: **Task 8's acceptance test fails with the frame (69.5% over) and
without it (54.5% over).** Deleting the function is not a fix, it is a smaller
defect. The requirement has to be restated in signature space before Task 8 can
be executed at all.

Also: **the doc comment goes with the function.** It asserts measured claims
(`median 0.74 / p90 1.34`) that G2's numbers supersede, and a future reader will
treat them as settled.

**On the SO(3)-invariant descriptor** — rust-dev's verdict is *viable, and an
improvement, if the invariant type is a **derived projection** rather than
independent state*:

- Give `Oriented<D>` **no** `PartialEq`/`Hash`/`Ord` derives. It then cannot key
  a map, be compared, or be sorted — so the only way to ask "same shape?" is
  `key()`, and the only way to ask "do these bind?" is `affinity`. The
  enforcement is the *absence* of derives, which is compiler-checked.
- The real axis of the split is **consumers that search rotations vs consumers
  that do not** — not identity vs binding.
- Property test: `a.permuted(r).key() == a.key()` for all 60 `r`.

**The float-sort tie-break objection dissolves, and this is worth recording**
because it will be raised again: `CLAUDE.md` requires an ID tie-break on
float-key sorts because sorting *records* by a float leaves equal-keyed records
in an implementation-defined permutation. `sorted_extents` sorts **the payload
itself** — two equal extents are indistinguishable bit patterns, so swapping
them is bit-identical. The hazard is absent, not overridden. Needs
`f64::total_cmp` and a comment saying why no tie-break is required.

**The real hazard rust-dev names instead:** `ShapeKey` is a map key and a
histogram bucket, so it needs a *total* order, and `Ord` cannot be derived on
`f64`. Store it as fixed-point `i64` — as `Mass` already does — and derive all
four. That also makes it golden-hashable without a float-formatting question.

**Routed to emergence-auditor, unanswered:** the quantisation bucket width is a
physics decision about how far apart two shapes must be to be different species.
Nobody has ruled on it.

---

### 2.3 · G3 — fold start — finding **`confirmed`**, **fix superseded**

**Was:** `[single-reviewer]`, verified arithmetically by me.
**Now:** finding `[corroborated]` (geometry + rust-perf, both reproduced it to
the monomer); **proposed fix `disputed` and replaced**.

This is the split the reconciliation table does not draw, so I am drawing it
explicitly: *a rejected fix is not a weakened finding.* rust-perf reproduced the
antipodal-pair layout `(0,3),(1,2),(4,7),(5,6),(8,11),(9,10)` and the collision
counts exactly — 21 at n=45, 164 at n=200 in a 13×13×7 box.

**On cost, the greedy walk is affordable and the question then dissolves:**

| n | current | greedy | anneal (20k) | greedy / anneal |
|---|---|---|---|---|
| 45 | 85 ns | 2.9 µs | 2.25 ms | 0.14% |
| 80 | 124 ns | 5.2 µs | 1.95 ms | 0.26% |
| 200 | 269 ns | 12.6 µs | 1.83 ms | 0.77% |

and better — the greedy walk takes no RNG and does not depend on `n`, so
`walk(n)` is the n-prefix of `walk(200)`; compute once in
`FoldWorkspace::new()` and every fold is a **48 ns prefix copy**.

**But the fix is rejected on correctness of the failure mode**, and both cheap
alternatives I offered as candidates were measured and rejected too:

- non-antipodal permutation `[0,4,8,1,5,9,2,6,10,3,7,11]`: 0 dups at n=45,
  **129 at n=200**
- stride 5: 0 at n=45, **135 at n=200**; stride 7: 4 dups at n=45; stride 11 ≡
  −1 mod 12, identical to the current code
- **"Avoiding the immediate antipode only defers the collision."**

**The replacement — a real serpentine, self-avoiding by construction:** walk
`a·u + b·v + c·w` over three FCC basis vectors `u=N[0], v=N[8], w=N[4]`. Cells
are `(a+c, a+b, b+c)`, and `x − y + z = 2c` recovers `c`, then `a`, then `b` —
so the map is **injective**, needing no occupancy probe and no bounds test.
Measured: connected, 0 duplicates, 0 out-of-bounds at 45/64/80/120/200, and
**192 ns at n=200 — cheaper than the current broken code.**

The greedy walk's trap branch is also dead code (`None => break` never fires for
any n ≤ 200, exhaustively checked) *and* a bug if it ever did: it yields
`coords.len() < n` while `Fold` keeps `n = p.len()`, after which `exposure`
panics and `rng.next_range(n)` indexes past the end. Make it a hard error.

**Does start compactness matter?** Partly, and this is new: post-anneal energy
at n=200 is greedy −149.1, serpentine −148.6, loose −122.0. **20,000 steps does
not erase a loose start at `MAX_POLYMER`** — 18% worse. At n=45–80 all three
converge within ~10%. (Caveat from the reviewer: synthetic energy table; treat
the ordering as indicative.)

**NEW consequence — a test becomes vacuous.** From any compact start the anneal
*reduces* contact count (727→642 at n=200) while lowering energy, trading raw
contacts for favourable class pairs. So `annealing_beats_the_extended_start`'s
`assert!(folded > 8)` is satisfied by the *start* (126–727 contacts). It must
compare against its own start's energy, not a constant.

**Routed on, not resolved:** the fold is **~2 ms per polymer species**, roughly
n-independent. Emergence E4 found polymer interning is per-condensation and
therefore reachable from `step`. 2 ms inside the step loop is the number that
decides the §17 budget, and the start walk is a rounding error inside it.
Nobody has costed that composition.

---

### 2.4 · E2 — contact energy — finding **`confirmed`**, **fix superseded**, and my caveat was understated

**Was:** `[single-reviewer]` `[coordinator-verified]`, with a caveat I raised
about the proposed fix.
**Now:** finding `[corroborated]` (emergence + geometry); **fix `disputed` and
replaced** with the product form `(fa·fb − 1.0) · w_charge`.

The sign inversion is real and urgent — only the *form* of the repair was
wrong. **My caveat was right and understated.** I noted that
`((fa+fb)²−4)·w` makes two neutral monomers *tie* a perfect complement.
Geometry measured what actually happens, folding at n=24, 20k steps, 8 seeds,
both tables normalised so only ranking differs:

| sequence | form A `((fa+fb)²−4)` | form B `(fa·fb−1)` |
|---|---|---|
| ALT +1/−1 | −0.600 | −0.600 |
| **near-neutral (±1/7)** | **−0.989** | −0.501 |
| MID ±5/7 | −0.760 | −0.520 |
| MIXED | −0.824 | −0.548 |
| ALL +1 | 0.000 | 0.000 |

> *"Under form A the caveat is worse than stated: the bland chain does not
> merely tie a perfect complement, it is the **global optimum**, beating every
> charged sequence… Form A makes charge a net liability — a standing pressure
> toward chemically featureless polymers, which surfaces as 'polymers are all
> the same, nothing recognises anything'."*

On the 8-class grid every near-neutral pair scores ≈ −4, so a bland chain gets
the best entry on *every* contact while a charged chain is frustrated and must
eat like-contacts at 0.

**The proposed test cannot discriminate, which is why this needed a
cross-check.** For `f ∈ {−1,+1}`, `((fa+fb)²−4)/4 ≡ (fa·fb−1)/2` **exactly** —
the normalised tables are the same function on the endpoints. Emergence's
proposed test (ALT vs ALL+1) passes identically under both forms,
Spearman ρ = 1.0000. **It is necessary, not sufficient.** An implementer would
have applied form A, watched the test pass, and shipped the pressure toward
bland polymers.

**Why form B is not a departure from §8.3** — a one-line permutation argument
that is testable and untested:

> `−(a_A+a_B)² = −a_A² − a_B² − 2a_A·a_B`. In the kernel
> `j = ANTI[PERM[R][i]]` is a bijection on directions, so `Σᵢ a_A[i]²` and
> `Σᵢ a_B[j]²` are both **R-invariant** — they cannot change which rotation
> wins. The only R-dependent part of §8.3's charge term is `−2Σ a_A a_B`, the
> product. Folding, where the contact set changes, promotes those inert
> self-terms into an unintended charge-burial penalty `Σᵢ kᵢ fᵢ²`.

So *"§8.3 verbatim"* imports a term **§8.3 never lets act**. Form B is §8.3's
charge term with the inert part dropped. This is the same class of error as the
original: transplanting an expression between a maximised score and a minimised
energy without checking which parts of it were doing work.

Form B satisfies the constraints: strictly monotone in `fa·fb`; perfect
complement uniquely best at −2w; every entry ≤ 0 so chains still compact
(neutral pair ≈ −1w); offset = 1 = max(`fa·fb`) fixed by construction, the same
justification emergence used for −4. **The offset must not later be raised to
make folds look good — that is the knob.**

**Honest pricing, from the reviewer:** within a fixed sequence the two forms
rank conformations near-identically (ρ ≥ 0.9906 over 298 conformations), so fold
*geometry* barely moves. The damage is confined to **cross-sequence** energy —
`Fold::energy` (`chemistry.md:342`), which has no caller today but returns
`Quanta` and is plainly destined to feed persistence. Low observable cost now;
form B costs nothing to adopt.

**Discriminating tests to add** (the ones that actually separate A from B):

- Table-level, exhaustive over 8×8: symmetry; `t ≤ 0` everywhere;
  `t(0,7) < t(3,4) < t(7,7)` — the first inequality is `−4 < −4`, **false under
  form A**; strict monotonicity in `fa·fb`.
- Fold-level: `E/contact(ALT ±1) < E/contact(near-neutral)`. Measured A: −0.600
  vs −0.989 (**fails**); B: −0.600 vs −0.501 (passes).

**Two independent corroborations arrived with this verdict:**

- **G3 confirmed a third time.** Transcribing the plan's start walk, geometry
  hit the antipodal retrace and *"tripped a self-avoidance assert"* without
  having been told about G3. Three reviewers have now reproduced it
  independently.
- **B-5 (`energy_scale` undefined) confirmed a fourth time**, with a new
  consequence: it *"must be derived from the table range (A→B halves it), or
  the implementer picks a number."* Adopting form B changes what that
  derivation must produce.

---

### 2.5 · alr-011 — volume `V` — **`confirmed`**, and my own reasoning was wrong

**Was:** `[single-reviewer]`, premise `[coordinator-verified]`, cost claim
asserted by me.
**Now:** `[corroborated]` — alife-researcher + rust-performance-expert.

The verdict confirms the finding and **falsifies the mechanism I offered**. See
§4 for the retraction. The engineering summary:

**Cost — right conclusion, wrong premise.** `rate()` (`chemistry.md:1896`) takes
`conc: &[f64]` and calls `det_math::exp` on **every invocation**, for up to
`FLAT_SCAN_LIMIT = 256` channels per step. A naive `/V` *does* land on the hot
path. Measured (aarch64, rustc 1.94.1, `-O`, 2×10⁷ iterations):

| form | ns/channel |
|---|---|
| `rate()` as planned | 6.84 |
| same + one `/V` | 7.18 (+5.0%) |
| prefactor hoisted; step does `c·n(n−1)/2` | 0.96 (**7.2×**) |

So `V` is free **iff** the prefactor is hoisted — which the plan should do
anyway, since temperature is constant per V0 run. Give `Reaction` a `c: f64`
and an `order: u8`, both resolved at channel creation.

**NEW corroborating defect, found in the same expression.** `rate()` **never
reads `u.consts.rate_prefactor`.** I verified: `rate_prefactor` occurs at
exactly four lines across both plan files — its declaration (`v0.md:2000`), its
generation (`v0.md:2049`), and two Task 20 *tests* that set it to `1e-30` and
`1e30` and expect universe rejection (`chemistry.md:2886`, `:2893`). It is read
nowhere. **Those two tests cannot pass as written.** One edit fixes the hot path
and the tests together.

**Where `V` lives:** beaker config, confirmed. `rate_prefactor` is drawn from the
universe stream and is genuinely universe character; `V` is an experimental
control. There is **no beaker config struct at all** yet — `Beaker::new(&u, &g,
seed)` — so introduce `BeakerConfig { volume, temp, seed, steps }` now rather
than a fourth positional argument later.

**The `V`/decay coupling — the reason to do this before Task 21, not Task 18:**

Unimolecular decay ∝ `k·n`; bimolecular synthesis ∝ `k·n(n−1)/(2V)`. The
synthesis/decay balance is a function of `n/V`, **never of counts**.

> *"A band locked before `V` exists is a band at one unstated concentration…
> Task 21 Step 2 writes that band into spec §22.6. Introducing `V ≠ 1`
> afterwards rescales every bimolecular propensity, moves the band, and §22.6
> is silently wrong with nothing failing."*

And: **the battery should hold `n/V` fixed, not counts.** Varying `V` at fixed
`n` is dilution — it changes the chemistry and confounds the test. Varying `V`
at fixed `n/V` leaves mean-field rates invariant and varies only fluctuation
magnitude (∝ 1/√n), which *is* the carrying-capacity axis §2.7's Moreno–Ofria
countermeasure describes. Task 20 does nothing equivalent: `run_battery(&u, &g)`
has no size axis at all.

**Retrofit cost:** ~9 sites, all unwritten. No golden moves — none exist. After
Task 21 it is the same 9 sites *plus* re-running the band search and rewriting
§22.6.

---

### 2.6 · D-1 — `lining_signature` — **DISPUTE CLOSED. Endorsement `retracted`.**

**Was:** the review's sharpest dispute — emergence endorsed precisely the `min`
and the `Span` conversion that geometry G4/G7 called silently wrong.
**Now:** `[corroborated]`, high confidence. Per the reconciliation table,
`disputed` + one retracted → **High**.

The emergence-auditor ran geometry's own proposed acceptance test and
**withdrew in full**:

> *"Geometry is right on both counts. I ran their proposed acceptance test and
> the result is worse than they claimed."*

| measure | plan (`min` over +hemisphere) | geometry: support over cavity cells | truth |
|---|---|---|---|
| range of `sig.r` (Span) | 2.2e-16 – 1.414 | 2.41 – 5.66 | 2.83 – 5.66 |
| long/short axis ratio | **1.00** | **2.00** | 2.00 |
| correlation with true half-width | **−0.0721** | **+0.9865** | 1.0 |

**16 of 42 directions return exactly `2.220446049250313e-16`** — one machine
epsilon, the rounding residue of an analytically-zero dot product from a lining
monomer sitting exactly in the equatorial plane, admitted by `reach > 0.0`
because the last bit went the right way.

And the number that settles it: **a 4/2/2 cigar and a 5/2/1.5 slab produce
bit-identical plan signatures** (L2 distance 0.0000). Cigar-vs-sphere,
rotated-cigar-vs-sphere and cigar-vs-rotated-cigar all give the *same* distance,
0.3068. Under support-of-cells those pairs separate by 2.20, 6.66 and 10.51.

**[coordinator-verified — I ran this independently before the verdict arrived.]**
Same ellipsoid, the plan's code as written at `chemistry.md:1361–1368`:

```
direction                plan: min over hemisphere    geometry: max over cells
+x (LONG axis, a=4)                       1.4142                      4.2426
+y (short, b=2)                           1.4142                      1.4142
+z (short, c=2)                           1.4142                      1.4142
diagonal                                  0.0000                      3.2660
```

The plan returns **the same number along the long axis and both short axes** —
exactly one lattice unit — and 0 on the diagonal. It measures the lattice
spacing, not the cavity. Over 200 random directions: plan min 0.0001, max
0.3597; support-of-cells 1.64–4.47, ratio 2.73 against a true anisotropy of
2.00.

**The `a` channel goes with it.** With `nearest ≈ 0`, the weight
`(nearest + SHELL − reach).max(0.0)` admits only `reach < 0.75`, so it averages
the equatorial ring: measured mean angle between `dir` and the weighted
monomers is **81.1°**, near-perpendicular to the direction being characterised.
`sig.a[d]` describes the wall a substrate approaching along `d` never touches.
Fixing `nearest` fixes it, provided the mirrored weight is re-derived.

**Why emergence graded it as sound, in its own words** — worth recording,
because it is a reusable failure mode:

> *"I graded structural mirroring as operational identity. The plan argues for
> itself that way at `:1348–1350` ('the same *structure* as §8.2's
> `signature()`')… Two passes with a mirrored weight is the same structure. It
> is not the same operation, and the difference is the entire finding."*

§8.2 is `max` over **the object's own constituents**; the plan is `min` over
**the complement**, hemisphere-restricted. Two changes, not one. Geometry's
`max over cells of (cell − centre)·dir` is literally §8.2 applied to the hole —
and emergence's measurement prefers it over geometry's own first-listed cone
fix (ratio 2.00 exact vs 1.67).

**The emergence consequence, which is why this outranks a numerics bug:**

> *"Apply the deletion test: remove the shape comparison entirely and catalysis
> still fires, on count alone. `MIN_ENCLOSURE = 8` and `MAX_CAVITY_CELLS = 24`
> become the operative definition of an enzyme — a pair of thresholds encoding
> an outcome, arrived at by accident rather than by anyone deciding to rig it.
> That is the failure I exist to catch, and I signed it off."*

**G7 confirmed, with one correction to geometry's numbers.**
`|(1,1,0)| = 1.414214` integer units × `LATTICE_SPAN = 1.414214` = **2.000000
Span**, not √2. I verified this arithmetically in both rounds. Geometry cited a
radius range of 0.56–4.55; emergence computed 0.56–4.095. **I settled it:**
`count = 4 + rng.next_range(3)` gives 4..=6 periods, iterated `for _ in
0..count`, so period indices are 0..=5 and the maximum is
`1.4 × (1 + 0.25×5) × 1.3 = 4.095`. **Emergence is right; geometry's 4.55
requires period 6, which cannot occur.** Immaterial to the finding — the
structural point (`LATTICE_SPAN` takes no `Universe` argument at all, while
§8.2's `r` carries `el.radius.get()`) is agreed by both.

**NEW, and it blocks the repair.** `centre = cells[cells.len() / 2]`
(`chemistry.md:1292`) is the median **flat index**, while `Cavity::centre`'s doc
(`:1190`) promises *"the cell nearest the cavity's centre of mass"*. On an
asymmetric 41-cell cavity emergence measured those two cells **4.0 lattice units
= 5.657 Span apart — further than the cavity's own extent.** Symmetric test
cavities hide it completely.

This upgrades **M-e** (previously a `[corroborated]` doc-comment mismatch,
geometry G11 + rust-dev #9) to a **blocker on the G4 fix**: support-of-cells
measures from `centre`, so an offset centre turns the repair into a support
function of the wrong body. **Fix `centre` before fixing `lining_signature`.**

---

### 2.7 · P5 — Bedau activity — **`disputed`**. The metric would report open-endedness in a system that has none.

**Was:** a `[single-reviewer]` performance finding that **explicitly requested
the alife-researcher and did not get it** — `synthesis.md` §5 G-1 item 3 named
this as the most costly consequence of that reviewer being dispatched without
`Write`. This cross-check closes that seam, and it was the single
highest-yield dispatch of the round.

**P5's asymptotics are endorsed; its formula is the wrong measure.**
Per-step sweeps over S are wrong, and per-*step* increments are wrong under
Gillespie — both confirmed. But `a_i = ∫ count_i dt` is **abundance-weighted**,
and the plan intends the **presence** form, and the difference is not cosmetic.

**What the published increment is.** Channon 2003 (ALIFE VIII, pp. 173–181,
eqs. 1–6) restates Bedau–Snyder–Packard verbatim:

- `Δ_i(t) = 1` if component i exists at t, else 0 — activity increment **by
  presence**. Channon: *"This is not the only activity increment that they have
  used, but it is the best for comparison across systems."*
- `a_i(t) = Σ Δ_i(τ)` **if i exists at t; 0 otherwise** ← the extinction reset
- `D(t) = #{i : a_i(t) > 0}`

Rechtsteiner & Bedau (ECAL'99) operationalise it identically: mean cumulative
activity is *"the mean age of the genotypes **present** in a system"*, and D is
*"simply the number of different genotypes present"*. Unweighted, extant-only.
The abundance-weighted variant is published (Standish, arXiv:nlin/0004026) but
answers a different question.

**The decisive argument is internal to Borbax, not bibliographic.** Spec §15.3
(spec:860) states the shadow is *"a drift control on **persistence only**:
catalysis produces differential formation rates… equalising it would destroy
the chemistry."* An abundance-weighted `a_i` puts **formation rate — the
variable the shadow deliberately does not control — inside the statistic the
shadow is supposed to calibrate.** Every activity excess would be partly
uncontrolled by construction, and criterion 9 would test the wrong thing.

The plan already says presence: `plan:2998` *"cumulative **existence** per
species"*, `spec:829` *"cumulative existence counters"*.

**The missing extinction reset is the severe half.** P5's accumulator never
resets on extinction. Measured, on a constructed class-2 beaker (D=20, constant
turnover, i.e. a system with **no** open-ended novelty):

| | t | 40× t |
|---|---|---|
| correct extant-only `Ā_cum` | 791 | 1210 (**bounded**) |
| P5's never-reset accumulator | 2.5e4 | 9.8e5 (**exactly linear**) |

> *"**A class-2 system reads class 3b.** That is the failure mode §2.7/spec:106
> cites the statistics to catch."*

Exit criterion 9 would pass on a system that does nothing. This is the
strongest single result of the cross-check round.

**Correct form, continuous time:**

```
Δ_i(t) = 1{count_i(t) > 0}
a_i(t) = ∫₀ᵗ Δ_i(τ)dτ   if count_i(t) > 0;   0 otherwise
D(t)   = #{i : count_i(t) > 0}
A_cum  = Σ_i a_i;  Ā_cum = A_cum/D;  Ã_cum = median{a_i : a_i > 0}
A_new  = (1/D)·Σ_{i: a_i ∈ [a₀,a₁]} a_i
```

**And it is cheaper than P5's version, not dearer** — which resolves the
performance concern more completely than P5 does. `Δ_i ≡ 1` while present, so
no per-event `touch()` is needed at all; only presence transitions matter:

```rust
present_time: Vec<f64>,  // ∫Δ dτ accrued over prior presence intervals
since: Vec<f64>,         // t of the current 0→1 crossing
// a_i(t) = if present { present_time[i] + (t - since[i]) } else { 0.0 }
```

O(1) on 0↔1 crossings only, **zero cost on ordinary firings**. It also removes
P5's accumulation-order concern — one add, not a running sum — so it is
*better* for §13.1 as well.

**Staleness (Q2a):** under the corrected form there is none — `a_i` is
closed-form in `t`. Under P5's form, flush-on-read is sufficient and free, but
his snippet does not do it, and unflushed values bias `A_new` against exactly
the long-persistent species the statistic rewards.

**Threshold (Q2b):** the cross-over method is unit-agnostic and self-calibrates,
**but `a₀` carries the increment's units (time vs count·time) and cannot be
transplanted between the two forms.** So this choice must be made *before* the
Task 20b Step 1/Step 3 coupling (S-b) is implemented, not after — it changes
what Step 3 must produce.

**Two further required changes:**

- **`plan:3000–3007` is wrong on `D`.** D is the count *present*, not "above
  the activity threshold". The threshold enters only `A_new`, as a band
  `[a₀,a₁]` normalised by D — not "activity of species newly crossing the
  threshold".
- **Track the median, not the mean, for the 3b test.** Quantified: adding 3
  non-extinguishing species to 20 grows mean `Ā` from 2.6e5 to 1.0e6 while the
  median stays ~650–970. **In a beaker with a persistent food set (RAF food
  set, `plan:2455`), `Ā_cum` is unbounded by construction.** Borbax has exactly
  such a food set. Channon 2003 makes this correction explicitly.

**[coordinator-verified] — every line reference in this verdict is accurate.**
I checked all four. `plan:3000` does document `diversity` as *"count of species
above the activity threshold"* (wrong per the above — it is the count
*present*); `plan:3006` does document `new_activity` as *"activity of species
newly crossing the threshold"* (wrong — it is a band `[a₀,a₁]` normalised by
D); `spec:829` does say *"cumulative existence counters"*; `spec:860` does say
*"a drift control on persistence only"*; and `plan:2455` does carry a
`food: &BTreeSet<SpeciesId>`, so the persistent-food-set argument that makes
`Ā_cum` unbounded applies to Borbax directly.

**Coverage gap G-1 is now closed.** `synthesis.md` §5 recorded that
`alife-researcher` had been dispatched without `Write`, producing 92 lines with
**zero** `file:line` references, zero fix code and zero proposed tests, and
called re-running it *"the cheapest single item in this whole review."* Given
`Write`, `WebSearch` and `WebFetch`, the same reviewer returned accurate line
anchors, two PDFs read in full, written-out formulae, working Rust, and two
computed counter-examples. **The recommendation was right and the cost of the
skill bug is now measured.**

**Scope note from the reviewer, recorded rather than hidden:** BSP 1998 itself
was not read this pass (the reed.edu link returned HTML, not the PDF) — relying
on Channon's and Standish's verbatim restatements; and the MODES paper (Dolson,
Vostinar, Wiser & Ofria, *Artificial Life* 25(1):50–73) is cited from
recollection. That bears on **alr-022**, which is about the MODES persistence
filter and remains `[single-reviewer]`.

---

### 2.8 · Benchmark 1 and P2 — two deferred calls, now decided

**Was:** two questions `rust-performance-expert` explicitly routed to
`determinism-auditor` and never got answered. `synthesis.md` recommendation 9
said of the first: *"this decision is free exactly once, and this is the
window."* It is now made.

#### Q1 · Fixed-shape pairwise summation in the affinity kernel — **`disputed`. Do not adopt.**

A clean "no", which is as useful as a "yes" and closes the item.

**It would be deterministic** — and the reason matters, because it is *not* the
`powi` situation. Verified on rustc 1.94.1 at `-O`: **zero** fast-math flags on
any FP op. LLVM reassociates floating point only under `reassoc`, which rustc
never sets, so source-level association survives any LLVM version. `powi` is
different because it lowers to `llvm.powi` and lets the **backend** pick the
multiply tree. *(This also sharpens the D-3 ruling in §3: the mechanism by
which `powi(6)` is a risk and a hand-written product is not, is now stated
precisely.)*

**But the accuracy gain does not clear the bar.** Measured, 200k trials against
`fsum` at realistic `ds = rᵃ + rᵇ − gap` magnitudes:

| D | serial worst-ULP | blocked(6) | pairwise | serial max-rel |
|---|---|---|---|---|
| 12 | 3 | 3 | 2 | 6.0e-16 |
| **42** | **6** | **3** | 2 | **1.0e-15** |
| 162 | 13 | 3 | 2 | 1.7e-15 |

> *"Both accumulators are sums of **squares** — all terms non-negative,
> condition number exactly 1, the best-conditioned case summation has. The
> `(n−1)u` worst case never engages."*

**On `CLAUDE.md`'s "2 above 3" clause — it does not apply, because rank 2 never
engages.** The kernel's output feeds an argmax over 60 rotations, a threshold
compare, and a logistic. A 1e-15 relative perturbation flips a decision only
when two candidates sit within 1e-15 — and in the case that actually occurs
(symmetric molecules, **exact** ties) both formulations produce bit-equal
values, so the documented lowest-index tie-break governs identically. There is
no physics-correctness claim here, so this is not *"a numerically better
formulation worth regenerating goldens"*; it is fewer ULPs with no consequence.

**And the performance premise did not survive measurement either:** serial
1810/1855 ns/call vs blocked 1837/1854 ns/call — no difference. The loop is
**gather-bound on `anti[perm[i]]`**, and `shape`/`charge` already provide two
independent dependency chains, so there is no serial-latency chain to break.
(One platform; if rust-perf measures a win on x86-64 the answer becomes "adopt
for speed, determinism has no objection" — *not* "adopt for accuracy".)

**Redirected:** SMACOF (`v0.md:3533–3550`) is the one accumulation with genuine
cancellation — `pos[j][k] + dij*diff[k]/dist` is signed — so it, not affinity,
is where an accuracy argument could exist. But it sits inside a 240-iteration
fixed point and is *"the last place to churn"*. `may_bind`'s `.sum()` is
`std`'s serial left fold, already deterministic.

#### Q2 · P2's dense-grid flood — **`confirmed` bit-identical, wrong justification, and M-c is DOWNGRADED**

**The conclusion holds; the stated reason does not.** rust-perf justified
bit-identity by *"box order is ascending flat-index order — exactly `BTreeSet`
iteration order."* That is true of `box_cells` feeding **Step 2**, but it does
not describe the flood. The correct reason:

> *"`outside` is **never iterated** — written by `insert`, read only by
> `contains`. Visit order comes from the `VecDeque`, and BFS reachability is
> order-independent, so any container with identical membership is
> bit-identical."*

And the reviewer names why filing the wrong justification matters: *"it would be
silently invalidated by a future change that iterates `outside`."* Record the
real reason.

**M-c — the torus-wrap finding — is DOWNGRADED. This is the only finding the
cross-check round made *less* severe.** M-c was `[corroborated]` (geometry G11 +
determinism), described as *"the flood leaks around the torus. Same class as S5,
now in Task 13."* Instrumented on a full-span box:

> **64,512 wrap events, 100% from a box-face source, and 0 that added a cell not
> already in `outside`.**

A wrap needs a component to leave [0,63], which needs the source at global 0 or
63; `bounding_box` expands by 2 and folding clamps to 2..62, so that is exactly
the lo/hi face — and `boundary_cells` already seeds *every* face cell. The
wrapped target lands on the opposite face, also seeded. Negative flat indices
decode to all-non-positive triples and fail `in_box` outright. Wrapping flood
vs strict per-axis flood: **identical `outside` sets, 6/6 full-span trials.**

**It is also explicitly *not* the trap I asked about** — equivalence survives
*fixing* the wrap, so P2's refactor is safe either way.

**Still fix it, for a structural reason rather than an observable one:** the
inertness depends on `bounding_box` expanding by exactly 2 *and* `fcc::in_bounds`
clamping folds to 2..62. **Loosen either constant and the leak goes live.** Add
`fcc::in_bounds(nb)` to the flood so the property is structural rather than a
coincidence of two constants. Routed to `geometry-numerics-reviewer`.

**Two ordering hazards in the replacement, both cheap:**

1. If the dense grid is allocated once and reused, clear it with `fill` or a
   generation-stamp `Vec<u32>`. A "clear only touched cells" list is a standing
   invitation to reintroduce an ordering dependency; a generation stamp is free.
2. **Do not extend the replacement to `candidates`.** `for &start in
   &candidates` (`chemistry.md:1256`) **is** iterated, and its ascending order
   fixes cavity discovery order and hence the order of the returned
   `Vec<Cavity>`. **That `BTreeSet` is load-bearing.**

---

### 2.9 · D-5 and D-2 — settled by measurement, not dispatched

These two disputes were not cross-checked because they are arithmetic. I ran
them.

#### D-5 · `MIN_ENCLOSURE = 8` — geometry is right, and both reviewers understated it

**Was:** `[disputed]` — emergence E10 called it *"a reporting filter"*
(confidence: medium, *"a judgement call, not a defect"*); geometry G6(1) called
it *"now redundant and harmful"*, fragmenting genuine cavities.

**[coordinator-verified.]** I built the most *favourable* case for the gate — a
void surrounded by solid polymer on every side, which maximises the occupied
neighbour count — and applied the ≥8-of-12 membership test:

| void size (cells) | 1 | 4 | 8 | 12 | **13** | 19 | 24 | 30 | 43 |
|---|---|---|---|---|---|---|---|---|---|
| cells passing the gate | 1 | 4 | 7 | 4 | **0** | 2 | 4 | 3 | 2 |

A cell in the middle of a void has *void* neighbours, so it cannot have 8
occupied ones. The consequence is not fragmentation, it is **erasure**: a
13-cell void — comfortably inside `MAX_CAVITY_CELLS = 24` — is reported as
**no cavity at all**. Larger voids survive only as 2–4 scattered cells. On
spherical voids up to radius 3.0 (55 cells) the gate passes **zero** cells; at
radius 3.5 it keeps 8 of 87, in **8 disconnected components**.

**Verdict: geometry's reading is correct and the severity is higher than either
reviewer stated.** Emergence's *"drop it to ~5"* does not fix it; the gate must
be **demoted to a per-component summary**, which is what both reviewers'
proposed fix already does. Emergence's independent point — that `enclosure` is
*already* used continuously in `catalysis_factor`, so the same quantity is
applied twice, once smoothly and once as a cliff — stands and is good.

**A composition of two existing findings, flagged for geometry's sign-off
rather than asserted as mine.** Geometry G6(1) says the gate mutilates voids
≥12 cells; emergence E10 says `MAX_CAVITY_CELLS = 24` *"discards large interior
voids, which are the pockets most likely to hold two substrates at once, i.e.
the §8.5 case."* Composed, the window in which a cavity is reported faithfully
is roughly **1–8 cells** — and §8.5's catalysis mechanism requires a cavity
large enough to hold two reactants. **The two constants together may leave no
size window in which the mechanism they exist to serve can operate.** Neither
reviewer stated this because neither saw both constants at once. It is a
question for `geometry-numerics-reviewer`, not a finding of mine — and it lands
on the same code as D-1, so settle them together.

#### D-2 · Is S8 (cavity enclosure) closed? — **not closed**, unchanged

Not re-dispatched; the `synthesis.md` verdict stands and is now reinforced.
Emergence's *"S8 closed"* rested on the outside-flood being a genuine enclosure
test, which is true and is a real improvement. It does not cover what
`MIN_ENCLOSURE` then does to the flood's output — which D-5 above now shows is
erasure. **Treat S8 as not closed.** Note also that the same reviewer has since
retracted its adjacent Task 13 endorsement (§2.6), which weakens the basis for
the S8 sign-off without contradicting it.

---

## 3. Open decisions — no longer disputes

**Zero reviewer-vs-reviewer disputes remain.** These two carried forward from
`synthesis.md`'s disputed set; both are now settled on the facts and what is
left is a decision to take, not a disagreement to adjudicate. Neither was
re-dispatched, because neither turned on judgement.

### D-3 · `.powi(` in `BANNED_CALLS` — ruling unchanged, count corrected

The `synthesis.md` ruling stands: **drop `.powi(` from `BANNED_CALLS`** (a
permanently-red gate is S10's failure mode returning) **and** write out the one
`.powi(6)` explicitly.

**Correction to my own prior count.** `synthesis.md` said "eight sites". The
accurate figure is **7 distinct lines, 12 occurrences**: `chemistry.md:2229,
2230, 2231` (1 each) and `v0.md:1722` (1), `:3300` (3), `:3332` (3), `:3884`
(2). Since `xtask` tests `line.contains(call)`, the operationally relevant
number is **7 failing lines**. Eleven of the twelve are `.powi(2)`, which LLVM's
`ExpandPowI` emits as `x*x` with no reassociation freedom; the single `.powi(6)`
at `v0.md:1722` has genuine tree-shape freedom and is pinned by the rustc pin.

**The ruling is now on firmer ground than when I made it.** §2.8's cross-check
supplies the mechanism I was reasoning about indirectly: rustc sets **no**
fast-math flags, so LLVM never reassociates floating point and *source-level*
association is stable across LLVM versions — whereas `powi` lowers to
`llvm.powi` and lets the **backend** choose the multiply tree. That is exactly
why writing `v0.md:1722` out explicitly is worth one line, and why the other
eleven `.powi(2)` sites are not a hazard at all.

### D-6 · `may_bind` — **resolved by G1's cross-check, in the opposite direction**

`synthesis.md` marked this contingent on G1: *"if G1 is right, the 100%
rejection is a symptom of G1, not of `may_bind`. Fix G1 first, then
re-measure."*

**G1 is right, and fixing it does not restore the filter — it destroys it.**
Per §2.1, `mean_defect ≡ 0` at `f = 1`, so `may_bind` goes from rejecting 100%
to accepting 100%. It is a valid bound and a useless filter at both ends.
D-6 is therefore no longer contingent: **`may_bind` needs a new bound derived
from sorted-extent variance**, and geometry's withdrawn tighter bound
(rearrangement inequality on sorted extents, retracted this round at 4.0% vs
3.3%) should be **reconsidered**, because it was withdrawn on a
cost/benefit comparison that G1 invalidates.

Add geometry's rate assertion regardless — it is the assertion the test's own
comment promises.

---

## 4. Coordinator retractions and corrections

Recording these because the reader cannot audit what they never see, and
because two of them are mine.

1. **RETRACTED — my `alr-011` cost mechanism.** `synthesis.md` §3.4 S-a says,
   marked `[coordinator]`: *"`c` is resolved **once per channel at channel
   creation**, not per step, so a division by `V` adds nothing to the step
   loop."* **This is false.** `rate()` recomputes `det_math::exp` on every
   invocation, for up to 256 channels per step; a naive `/V` costs +5%. The
   conclusion (expose `V` now) survives; the reasoning does not. The correct
   argument is the one rust-perf gave: hoist the prefactor — which the plan
   should do anyway — and *then* `V` is free.
2. **CORRECTED — the `.powi(` count.** Eight → 7 lines / 12 occurrences (§3).
3. **PARTIALLY RETRACTED, not mine — geometry's *"the `may_bind` bound
   survives unchanged"*** inside G1. Valid, yes; a filter, no (§2.1).
4. **SETTLED — geometry's radius upper bound 4.55 → 4.095** (§2.6).
   Immaterial to the finding.

---

## 5. Mechanical claim-audit — Track C, actually run

`synthesis.md` recommended a mechanical grep audit as *"the highest
value-per-minute item in the whole review"*. I ran a first pass. It took under
a minute and found a **tenth** compile-caught defect, after two review rounds by
six specialists.

**These are mechanical-audit output, not review findings.** They are instances
of an already-corroborated class (B-2…B-6: identifier used, never defined), not
new judgements of mine.

- **`ALL_CRATES` is undefined.** It occurs **exactly once** in either plan file
  — `v0.md:943`, `for krate in ALL_CRATES {`, inside
  `check_no_platform_transcendentals`. There is a `CHEMISTRY_CRATES`
  (`v0.md:307`) but no `ALL_CRATES`. **`E0425`. `xtask` does not compile, so the
  §5 fiction gate never runs.** Same shape as B-4 (`FoldId`) and B-5
  (`energy_scale`).
- **`FLAT_SCAN_LIMIT` exists only in prose** (`chemistry.md:2729`, a doc
  comment: *"start at 256"*). It is never declared. This sharpens determinism
  D-f — the constant selecting between two never-proven-equal propensity paths
  is not merely unproven, it does not exist.

Both are found by `grep -c` on an identifier. The recommendation stands and
should be an `xtask` check.

---

## 5b. Coverage gaps — one closed, three open, one new

| Gap | Status |
|---|---|
| **G-1** · `alife-researcher` dispatched without `Write` | **CLOSED.** Re-run with full tools via the P5 cross-check; output quality difference measured in §2.7. |
| **G-2** · Nobody reviewed the spec | **OPEN.** Still a dispatch bug, and it grew: §2.7 now finds a **spec-level** defect (`spec:829`/`spec:860` vs the plan's activity form) alongside alr-010 and alr-034. |
| **G-3** · Task 19 (rendering) has no owner | **OPEN.** Unchanged. `CLAUDE.md` makes the renderer a correctness instrument — *"a shape-based chemistry cannot be developed blind"* — and nobody has reviewed the projection or the `insta` golden mechanism. |
| **G-4** · Tasks 1–3 reviewed by determinism only | **OPEN**, and now less benign: §5's mechanical audit found `ALL_CRATES` undefined in Task 1's `xtask`, which means the §5 fiction gate has never run. That is precedence rank 1. |
| **G-5** · **NEW** — nobody owns the acceptance criteria themselves | **OPEN.** §2.7 shows criterion 9 would pass on a class-2 system. Six reviewers audited the *implementation*; none audited whether §23's criteria measure what they claim. Raising it as a gap, not a finding. |

**Two questions were routed by cross-checks and have no owner:**

- **The `ShapeKey` quantisation bucket width** (§2.2) — rust-dev explicitly
  declined it as a physics decision for `emergence-auditor`.
- **`mean` vs closest-approach `max` in `pair_gap`** (§2.1) — emergence
  explicitly declined it as a geometry question. `max` is the physically
  correct separation *and* is rotation-dependent, so it would participate in
  the search rather than factoring out of it.

Both are one dispatch each and both gate work in Track A.

---

## 6. Third fix pass, or something structural?

**Still structural. The cross-checks did not change the answer — they replaced
its main argument with a much stronger one, and made it more urgent.**

### The prior argument, and why it is now the weaker one

`synthesis.md` §7 argued: fix pass 1 introduced four new defects; the failure
has a shape (*"the fix landed in the prose describing the code and not in the
code"*); nine defects were compile-catchable; therefore split into three tracks.
All of that survives — and the compile-catchable count is now **ten** (§5).

But the cross-check round produced a sharper finding that subsumes it.

### The decisive new evidence: findings are reliable, fixes are not

Eight cross-checks dispatched, eight returned. The asymmetry is stark:

| | count |
|---|---|
| **Findings confirmed** | **7 of 7** |
| Findings retracted by their originator | **0** |
| Findings downgraded on measurement | 1 (M-c, inert) |
| Endorsements retracted after measuring | 1 (settled the review's hardest dispute) |
| **Proposed fixes, remedies or justifications rejected / amended** | **6 of 7** |

The six:

1. **G1's `gap_factor`** — rejected. A new generated constant whose only effect
   is a size-penalty knob. `f = 1` is structural; the knob is not.
2. **G3's greedy walk** — rejected. Superseded by a serpentine that is
   self-avoiding *by construction* and **cheaper than the code it replaces**.
   Both cheap alternatives I proposed were measured and rejected too.
3. **E2's `((fa+fb)²−4)·w`** — rejected. Makes a bland homopolymer the
   **global optimum**, a standing pressure toward featureless chemistry. And
   **its own proposed test cannot detect this**: on `f ∈ {±1}` the two forms are
   the same function, ρ = 1.0000.
4. **P5's `∫ count_i dt`** — rejected. Abundance-weighted where the plan and
   the shadow both require presence-weighted, and **missing the extinction
   reset, so a class-2 system reads as class 3b**.
5. **D-1's first-listed cone fix** — amended. Geometry's *second* suggestion
   (support over cavity cells) measured better than its first: ratio 2.00 exact
   vs 1.67, correlation 0.9865 vs −0.0721.
6. **P2's bit-identity justification** — wrong, while the conclusion was right.
   The refactor is safe because `outside` is never iterated, not because box
   order matches `BTreeSet` order. Filing the wrong reason is itself a hazard.

Plus one that was mine: **my `alr-011` cost mechanism was simply false** (§4).

**The one clean "no" is also a result.** The determinism-auditor answered
benchmark 1's deferred question with *do not adopt* — measured, with the
perf premise falsified alongside the accuracy argument. `synthesis.md` said
that decision was *"free exactly once, and this is the window"*. It is now
made, and it cost one dispatch.

**Read that table again.** Every reviewer was right about *what is broken*.
Most were wrong about *how to fix it* — and in three cases the proposed fix
came with a test that would have passed while the defect remained. A third fix
pass applying `synthesis.md` as written would have implemented four wrong
fixes and shipped three tests that certify them.

**That is fix pass 1's failure mode, exactly, about to repeat.** Fix pass 1
applied ~30 findings and introduced four new defects. The mechanism was never
"the reviewers were careless"; it is that **proposing a fix is a different
activity from finding a defect, and the plan-editing process has no step that
validates the fix before it lands.**

### Two findings that are not defects and do not belong in any fix pass

- **G2 (unchanged from `synthesis.md`, and worse).** `canonicalise_frame` loses
  to its own absence — *and* rust-dev confirms **Task 8 does not land either
  way**: 69.5% over the bound with the frame, 54.5% without. This is a possible
  falsification of §8.4's evolvability premise, not a task defect.
- **P5/alife (new this round).** Exit criterion 9 would pass on a system with
  no open-ended novelty. That is a defect **in the exit criteria**, which are
  the instrument that decides whether the whole project worked. Fixing tasks
  while the acceptance instrument is miscalibrated is the most expensive
  possible ordering.

Both are **measurement questions with no owner**, and both gate work downstream
of them.

### So: the same three tracks, with one addition

**Track A — settle the physics experimentally, before Task 8 lands.**
Unchanged, and now with a much better-stocked toolbox: geometry built a Python
reproduction of the pipeline; emergence built an FCC cavity harness and an
annealer; rust-perf built a Rust fold-start benchmark; alife built activity-class
simulators; I built lattice and gap models. **Six harnesses now exist in six
scratchpads and will be deleted.** Promote them to one checked-in
`experiments/` crate — that is the single highest-leverage action available,
and it is nearly free because the code is already written.

Gate list, in order: **G2** (does signature space have locality at all?),
**P5/criterion 9** (does the acceptance instrument work?), then G1, D-1, E2,
G10, D-5.

**Track B — execute Task 1.** Unchanged. Ten compile-caught defects, one of
which (`ALL_CRATES` undefined) means **`xtask` does not compile, so the §5
fiction gate — precedence rank 1, the one thing `CLAUDE.md` calls absolute —
has never run.**

**Track C — the mechanical claim-audit.** I ran a first pass in under a minute
and it found the tenth compile-caught defect and sharpened determinism D-f.
Make it an `xtask` check.

**Track D — NEW: no fix lands without a cross-check from a different
specialist.** This is the process change the round actually earned. It costs
one dispatch per fix, it caught 5 of 6 wrong fixes, and every catch came from a
specialist in a *different* domain **running the numbers** rather than reading
the diff. Concretely: the reviewer who finds a defect proposes the fix; a
second specialist measures it against a discriminating test before it is
written into the plan. Note that in three of five cases the *originator's own
test* would have passed — so "add a test" is not a substitute for this step.

### On the Step-5 escalation threshold

Contradiction density fell from **6 to 2**, well below the threshold. I am
explicitly **not** escalating. The prior round's six conflicts were, as
`synthesis.md` suspected, a signal about *method* rather than incoherence:
four of them were reviewers-who-measured against reviewers-who-read, and **when
the readers measured, they agreed — unanimously and immediately.** The
emergence-auditor's retraction is the cleanest evidence in the whole review that
this artifact yields to measurement and resists reading.

That is also the argument for Track D, and against a third fix pass conducted
the way the first two were.
