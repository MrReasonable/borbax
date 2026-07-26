# Emergence audit — V0 plan, post-fix re-review (2026-07-26)

Scope: both plan files. Verifying the fixes recorded in
`docs/review/2026-07-26-plan-audit.md` and judging the new code they introduced.

Headline: **two of the four physics fixes did not land in the code they claim to
fix.** Task 16's fix is prose only; Task 11's fix is sign-inverted. Both are
load-bearing for exit criterion 6.

---

## E1 · Task 16 Step 3 — the catalysis fix is a comment, not a change. SEVERE

`plans/...-v0-chemistry.md:2169-2254`.

The module header (2169–2191) describes at length how `IDEAL_SEPARATION`,
`BIND_THRESHOLD`, `MAX_ENHANCEMENT` and the `* 1.0e3` multiplier were removed.
The function body below it is **byte-for-byte the pre-fix version**:

- `2215`, `2222`: `if bind_a < BIND_THRESHOLD` — constant no longer defined.
- `2237`: `(sep - IDEAL_SEPARATION).abs() / IDEAL_SEPARATION` — still there.
- `2249`: `let factor = 1.0 + proximity * enclosure * fit * 1.0e3;` — the
  multiplier the header says moved to `UniverseConsts.catalytic_prefactor`.
- `2250`: `factor.min(MAX_ENHANCEMENT)`.
- `BIND_FRACTION` (2195) is declared and **never read**.
- `catalytic_prefactor` appears nowhere in either plan file except in that
  comment. `UniverseConsts` (v0.md:1992–2051) does not have the field.
- The Interfaces line (2096) still exports `IDEAL_SEPARATION`,
  `MAX_ENHANCEMENT`; the tests (2123, 2130, 2137, 2138) still use them.

Which principle: §3.1 "if we find ourselves writing a special case for life,
the physics is wrong"; §22.2 (the resolution sweep must measure resolution).

Detector or mechanism: **mechanism.** Delete the `1.0e3` and catalysis
disappears — `proximity·enclosure·fit` is ≲1, so the factor collapses to ~2.
That literal *is* the answer to exit criterion 6, and the test at 2125 asserts
`f > 2.0`, which is calibrated to it. `BIND_THRESHOLD = -12.0` scales with D
and with generated `w_shape`/`w_charge`/`ideal_gap`, so the §22.2 sweep over
D ∈ {12,42,162} measures the threshold moving, not the resolution.

Why it matters: criterion 6 is "the one that matters" (§23). As written it
reports a fact about two literals.

Fix: apply the header's own prescription to the body.

```rust
// Reference score: a perfect complement scores 0, so scale the cutoff by the
// worst plausible score for this pair at this D. Dimensionless => resolution
// independent (§22.2).
let ref_score = reference_affinity(&ca.signature, sub_a, g, k); // = -w_shape*D*ideal_gap^2, say
if bind_a < ref_score * (1.0 - BIND_FRACTION) { continue; }
...
let sep_ideal = k.ideal_gap + sub_a.mean_extent() + sub_b.mean_extent();
let factor = 1.0 + proximity * enclosure * fit * k.catalytic_prefactor;
best = best.max(factor);   // no clamp: the three factors already bound it
```
and add `catalytic_prefactor: rng.next_f64_range(..)` to `UniverseConsts`,
plus the field to `BindConsts`.

Is `BIND_FRACTION = 0.25` resolution-independent as claimed? Only if it is
actually used as a *fraction of a score that scales the same way*. As a bare
constant it is neither used nor defined against a reference, so the claim is
currently unverifiable. Written as above (fraction of a D-scaled reference) it
is genuinely dimensionless. Confidence: high on the finding, medium on the
exact reference-score form.

---

## E2 · Task 11 Step 3 — the contact-energy fix is inverted. SEVERE

`plans/...-v0-chemistry.md:608`, comment 589–607.

```rust
t[pair_index(a, b)] = -(fa + fb) * (fa + fb) * u.consts.w_charge;
```

`w_charge ∈ [0.6, 1.4]`, always positive (v0.md:2048). The annealer
**minimises**: `delta = trial_energy - current; accept = delta <= 0.0 || ...`
(chemistry.md:517–519).

Work it through, per contact:

| fa | fb | table entry | minimiser's preference |
|---|---|---|---|
| +1 | −1 (complementary) | `−0² · w = 0` | worst |
| +1 | +1 (like) | `−4w` | **best** |
| −1 | −1 (like) | `−4w` | **best** |

So like-attracts-like — the exact defect M1 raised, reintroduced with the
opposite sign error. §8.3 uses `−(a_A + a_B)²` inside a quantity it
**maximises** (`affinity(A,B) = max over R of score`), where it correctly
rewards `a_A ≈ −a_B`. Transplanting the expression verbatim into a minimised
energy reverses its meaning. "§8.3's charge term, verbatim" is true of the
characters and false of the physics.

Detector or mechanism: **mechanism.** Every fold in the engine is produced by
this ranking, and every cavity, every solvent-exposure figure and every
catalysis result is downstream of it.

`annealing_beats_the_extended_start` (109–116) cannot catch it: with the
current sign every table entry is ≤ 0, so *all* contacts are favourable and
the chain compacts regardless of complementarity. Same structural blind spot
as the original.

Fix — negate and add a contact-vs-no-contact reference offset, so
complementary contacts are the most favourable *and* contacts remain
favourable at all (a plain `+(fa+fb)²w` would be ≥ 0 everywhere and the chain
would refuse to compact):

```rust
// §8.3's charge term expressed as an *energy*. §8.3 maximises -(a+b)²; an
// annealer minimises, so the sign flips. The -4 is the reference state (a
// like-like contact is neutral, not attractive), which keeps compaction
// favourable while ranking complementarity best.
t[pair_index(a, b)] = ((fa + fb) * (fa + fb) - 4.0) * u.consts.w_charge;
```
gives −4w for a perfect complement, 0 for a like pair.

Add a test the current one cannot substitute for: fold an alternating
`+1,−1,+1,−1` sequence and an all-`+1` sequence of the same length, and assert
the alternating one reaches lower energy *per contact*. Confidence: high.

---

## E3 · Task 18 Step 5 — catalysis recompute is on the hot path and unmemoised. SEVERE (§8.6)

`plans/...-v0-chemistry.md:2725-2731`.

> "Recompute `catalysis` only when a *cavity-bearing* species' count crosses zero"

Three problems.

1. **Not well-defined.** "Crosses zero" does not say 0→1, 1→0 or both, and does
   not say which channels are recomputed or how a catalyst is joined to the
   channels it acts on. It also contradicts Task 14's `Reaction::catalysis` doc
   (1875–1882), which says the value is "resolved once when the channel is
   created". Two statements, two different lifetimes.
2. **It is most frequent exactly where it matters least affordably.** At small
   counts — the regime in which a RAF nucleates — a catalyst present in one or
   two copies crosses zero on most firings that touch it. So the trigger fires
   near-continuously precisely in the interesting regime.
3. **It calls `affinity_ordered` and nothing memoises it.** `catalysis_factor`
   runs `affinity_ordered(&cavity.signature, substrate, …)` — 2,520 iterations
   at D = 42 — per ordered cavity pair per channel. `AffinityMemo`
   (1799–1839) is keyed on `(SpeciesId, SpeciesId)`; a cavity signature has no
   `SpeciesId`, so the one affinity call now reachable from `step` is the one
   affinity call the memo does not cover.

Which principle: §8.6 "no code reachable from a simulation step may call any of
them"; §17 budget.

Detector or mechanism: neither — it is a scheduling defect, but it is the §8.6
regression "arriving disguised as a small feature" that CLAUDE.md warns about.

Fix: make `catalysis` genuinely per-(channel, catalyst-species) and computed
once, at the point the *catalyst species* is interned (it is a pure function of
the catalyst's cavities and the two substrate signatures — nothing per-molecule
about it). Store a `Vec<(SpeciesId, f64)>` of candidate catalysts on the
channel, sorted by id. Then a count change does no geometry at all: it selects
the best factor among catalysts with `count > 0`, which is a max over a short
sorted list. Add a `CavityAffinityMemo` keyed on
`(catalyst SpeciesId, cavity index, substrate SpeciesId)` for the intern-time
computation. Confidence: high.

---

## E4 · Task 18 Step 6 — interning a *polymer* is not a rare event. MEDIUM-SEVERE

`plans/...-v0-chemistry.md:2732-2733`: "When a genuinely novel species appears,
intern it and create its channels. **This is the only path on which
canonicalisation runs.**"

§8.6 sanctions this for canonicalisation because a novel species is rare. It is
not rare for polymers: every condensation that extends a chain produces a
sequence never seen before, by construction, so the Task 12 fold cache misses
*by construction* on exactly this path. Each such intern is `ANNEAL_STEPS`
(~20k) Metropolis steps plus the cavity flood fill plus a D-direction lining
signature — inside `step`.

The S2 fix (fold and cavities in `SpeciesRecord`) removed the *repeat* cost.
It did not bound the *first-sight* cost, and the plan says nothing about it.

Detector or mechanism: mechanism-adjacent. It does not rig an outcome, but it
is the §8.6 hazard the fix was supposed to close, still open on the polymer
path.

Fix: state the amortisation explicitly, and make `fold`/`cavities` **lazy per
species** — computed on first query rather than eagerly at intern — so a
polymer that appears once and decays never pays for a fold it is never asked
about. Lazy-per-species is still per-species and does not breach §8.6. Record
the measured novel-polymer rate per 10⁶ steps in Task 21. Confidence: medium
(the right bound depends on the measured rate).

---

## E5 · Task 14 Step 3 — solvent attack is size-blind, and branches on representation. MEDIUM

`plans/...-v0-chemistry.md:1752-1770`.

**Sign: now correct.** `fit = affinity_ordered(...)` is ≤ 0 with 0 the perfect
fit; `bind_probability` is monotone increasing in its first argument
(v0.md:4312–4320), so the rate rises with complementarity. `exposed_fraction`
∈ [0,1] falls with burial. Both directions are right. S4 is resolved.

Two residual problems, and they point the same way.

**(a) The rate is intensive.** `solvent_rate = base_solvent * exposed_fraction`
is never multiplied by `bond_count`. A 200-mer with mean exposure 0.5 is
attacked at *half* the per-molecule rate of a fully exposed dimer, despite
having ~200× the attackable bonds. §9.5 says the solvent "binds exposed *bonds*
and cleaves them" — the propensity must scale with the number of exposed bonds.
Compounding it, `fit` is the *whole-molecule* signature against a lone solvent
atom: a big folded globule has huge `r` against a tiny one, so the shape defect
is large, `fit` is very negative, and `base_solvent` is small. Size therefore
protects twice, and neither halving is derived from anything.

**(b) `exposed_fraction = 1.0` for small molecules** (1768) is a physics branch
on `SpeciesKind`. Defensible in isolation — a small molecule has no interior —
but combined with (a) it means the same physical object gets a different decay
law depending on which constructor made it, and the discontinuity lands exactly
at the small-molecule/polymer boundary where polymerisation has to compete.

Which principle: §3.6 "load the dice, don't rig the game". Neither branch is
conditioned on *being a replicator* — this is loaded dice, not a rigged game —
but the loading points at the outcome the run is trying to produce, and nothing
in the spec asks for it.

Detector or mechanism: **mechanism.** Delete the size-blindness and long
polymers become materially harder to keep; that is a different chemistry.

Fix: derive it from the one mechanism, per bond:

```rust
// One propensity per exposed bond, not one per molecule (§9.5).
let solvent_rate: f64 = (0..bond_count).map(|b| {
    let local = bond_local_signature(b, ...);      // the near-wall construction
    let fit = affinity_ordered(&local, &solvent_sig, g, k);   // already exists
    bind_probability(fit, k.reference_temp, k.bind_midpoint) * bond_exposure(b)
}).sum();
```
where `bond_exposure` is 1.0 for a small molecule and the §9.5 empty-neighbour
fraction for a folded monomer pair. This removes the `SpeciesKind` branch
(exposure is 1.0 because nothing is buried, not because of a type tag), makes
the rate extensive in bonds, and reuses the local-signature construction Task
13 already builds for cavities rather than comparing a globule to an atom.
Confidence: high on (a), medium on (b).

**Also broken mechanically here:** `intern` pushes
`SpeciesRecord { canon, sig, mass, bond_count, cleave_propensity, solvent_rate }`
(1772–1779) against a struct declaring `kind, sig, mass, bond_count,
radiogenic_rate, cleave_propensity, solvent_rate, fold, cavities`. `fold` is
read at 1762 and never bound. `intern_polymer` is used at 1551–1552 and never
defined. `u.consts.reference_temp` and `u.consts.bind_midpoint` (1754) are not
fields of `UniverseConsts`. The polymer intern path — the whole point of the S2
fix — does not exist in the plan.

---

## E6 · Task 15 — the §8.6 fix was not propagated. MEDIUM

`plans/...-v0-chemistry.md:2074-2081`. `SpeciesRecord::radiogenic_rate` was
added (1646–1651) with a doc comment explaining why per-call computation was
wrong. Task 15 still computes it per call, still carries the old
"hoisted if it shows up in a profile" comment, and still reads
`it.record(id).canon` — a field that no longer exists. `decay_channels`
(2054–2072) still allocates a `Vec` per call, on the path taken every time a
count changes.

Fix: `decay_channels` writes into a caller-supplied `&mut [(DecayKind,
SpeciesId, f64); 3]`; `radiogenic_rate` reads `record.radiogenic_rate`.
Confidence: high, low severity — it is a missed edit, not a design error.

---

## E7 · Task 4 — G3: the valence rule and the stability curve are still real. MEDIUM

`plans/...-v0.md:1706` and `1722`.

The fix landed for period lengths (drawn, 1632–1643) and for trend *directions*
(`affinity_sign`, `radius_period_sign`, `radius_group_sign`, 1669–1671). Those
are genuinely fixed and are the two that mattered most.

What remains real:

**Valence** (1706): `max_valence * (1 - |2·across - 1|)` — a triangle peaking
mid-period, falling to zero at both ends. That is real valence behaviour, and
the comment names it as such ("falls to zero at the noble end" — also a G6
framing hit). Only the *amplitude* is drawn; the *rule* is fixed across every
universe. G3 says "Valence rules are generated", and §7.1's fix note is
explicit that generated magnitude is not enough. This is the same defect as the
trend signs, one field over.

Fix: draw the profile, not just its height —
```rust
let lobes = 1 + rng.next_range(3) as f64;          // 1..=3 valence peaks
let phase = rng.next_f64();                        // peak position, not always mid-period
let v = (f64::from(max_valence)
    * (0.5 - 0.5 * det_math::cos(TAU * (lobes * across + phase)))).round() as u8;
```
A two-lobed period is a table nothing in real chemistry looks like, which is
the point.

**Stability** (1722): `(index/90.0).powi(6)` — heavier elements less stable, no
drawn sign, and `90` sits at the real onset of radioactivity. With the `>= 120`
element cap (1680) this is the most identifiable real-chemistry echo left. It
also feeds a live physics channel (radiogenic decay and mutation, §9.5).

Fix: draw the instability profile — `unstable_onset ∈ [0.3, 0.95]` as a
fraction of table size, `sharpness ∈ [2, 10]`, and optionally a drawn choice
between "rises with index" and "peaks in a mid-table band". Which elements are
your energy-and-mutation source then becomes a per-universe fact.

Also `.powi(6)` at 1722 is a native call on a path that sets `stability`; `powi`
is repeated multiplication so it is exactly specified, but it sits beside a
`powf` that was deliberately routed through `det_math` — worth a comment saying
why one is routed and the other is not.

**Checked and clean here:** `mass_base * (index+1)^mass_growth` reproduces
"mass increases down the table", but index *is* the table ordering, so this is a
convention rather than a trend claim, and §7.1 sanctions it explicitly. Not a
finding. Confidence on valence: high. On stability: medium-high.

---

## E8 · `catalytic_class` is a second catalysis mechanism waiting to be used. LOW (V0), watch for V1

`plans/...-v0.md:1735`, spec §7.1. `catalytic_class: group % n_catalytic` —
"which reaction families this element promotes when exposed on a mineral
surface". It is written and never read anywhere in V0, so no violation today.

But it is a catalysis path keyed on an *element label*, parallel to §8.5's
cavity geometry. Principle 2: "if a feature needs a second mechanism, question
the feature". When V1 adds mineral surfaces, this must be derived from
`affinity(mineral_surface_signature, molecule)` — the same operation at another
scale — not from a class integer. Flagging now because the field's existence is
what will make the shortcut look natural later. Confidence: high that it is a
smell, low urgency.

---

## E9 · Task 20b — clean instrumentation, two under-specified knobs. LOW

`plans/...-v0-chemistry.md:2968-3082`. This is measurement throughout. Nothing
in it can cause a phenomenon; deleting the whole task leaves the chemistry
identical. No special case for life. The `min`-not-mean novelty choice, the
constant-N histogram subsampling, the increments-not-cumulative fit, and the
held-out tail are all guards against measuring an answer into existence, and
the null-trajectory validation (3034, 3080) is the right acceptance test.

Two places an implementer can still build a fitness function by accident:

**(a) `ActivityStats::diversity` = "count of species above the activity
threshold" (3001).** The threshold is undefined. In Bedau–Packard it is
calibrated *against the neutral shadow* — that is the whole reason the shadow
exists. If someone picks a literal, it becomes a knob tunable until the graph
looks good. Specify: the threshold is the shadow's per-species activity
distribution at a stated quantile, computed from the shadow run, and the plan
must forbid a compile-time default.

**(b) "The common rate is set so total removal flux matches the focal run at
the fork point" (3038–3041)** states the constraint but not the formula. Write
it: `k_common = (Σ_i n_i · k_i) / (Σ_i n_i)`, count-weighted at the fork
keyframe, and **held fixed thereafter** — recomputing it as the shadow drifts
would feed the focal run's composition back into the control, which is no
longer a control. An unweighted mean over species does not match flux and is
the natural wrong guess. Confidence: high; both are one-line spec additions.

---

## E10 · Task 13 — `MAX_CAVITY_CELLS` outlived its justification. LOW

`plans/...-v0-chemistry.md:1173`: `const MAX_CAVITY_CELLS: usize = 24;` —
"Anything bigger is the outside world."

That justification was true before the fix and is not true now: the
outside-flood (1224–1239) already establishes what is outside, so anything
surviving it is interior by construction. What the constant now does is discard
large interior voids — which are the pockets most likely to hold two substrates
at once, i.e. the §8.5 case. Deletion test: delete it and *more* cavities are
reported, so it is a filter, not a mechanism. But it caps the substrate size a
cavity can serve, by fiat, with a stale reason attached.

Fix: keep a bound (an unbounded flood is a real hazard) but re-justify it —
either as a fraction of the chain's own volume (`n/4` cells, so it scales with
the molecule rather than being absolute), or state plainly that it exists to
bound the flood and is not a claim about what a cavity is.

**`MIN_ENCLOSURE = 8` survives, and is now a detector.** With a real enclosure
test alongside it, deleting it does not remove cavities — it reports more,
shallower ones. It is a reporting filter. One caveat worth a line: `enclosure`
is *already* used continuously in `catalysis_factor` (`ca.enclosure.min(...)
/ 12.0`), so the hard cutoff at 8 applies the same quantity twice, once
smoothly and once as a cliff. Preferring the smooth use and dropping the cutoff
to a low sanity bound (say 5) would let the chemistry decide where the useful
enclosure lies rather than the constant deciding it. Confidence: medium — this
is a judgement call, not a defect.

---

## Fiction guarantees, G1–G6

- **G1** clean. No data files, no tables of real values anywhere in the
  chemistry crates.
- **G2** clean. `REAL_ELEMENT_SYMBOLS` / `REAL_WORDS` (v0.md:1394–1425) are
  exclusion-only, never read as data, and `blocklist_actually_rejects`
  (1366–1372) exercises both directions. The `xtask` name mismatch from the
  previous audit is fixed — the const is `REAL_ELEMENT_SYMBOLS`.
- **G3** — see E7. Period lengths and trend signs fixed; valence rule and
  stability curve still real. `affinity` carries an explicit "not
  electronegativity" disclaimer (v0.md:1596–1598) and its derivation bears that
  out.
- **G4** — `Thermal`, `Quanta`, `Span`, `WorldYear`, `Mass` exist as newtypes
  with compile-fail tests. Leaks remain and are unchanged from the previous
  audit: `Signature.r[]`, `ideal_gap`, `Element.stability`, `cleave_propensity`,
  `solvent_rate` and `catalysis` are all bare `f64`. `bind_probability(affinity:
  f64, t: Thermal, midpoint: f64)` divides a bare `f64` by a `Thermal` — the
  one place where the type system would have caught a unit error and does not.
  Low individually; collectively this is how G4 is defeated quietly.
- **G5** clean. No format import/export anywhere.
- **G6** — small. "falls to zero at the noble end" (v0.md:1704) imports real
  vocabulary into the generator's own justification. The `enzyme` / `membrane`
  mentions in the binding module (v0.md:4159, 4353) and Task 16's header
  (2148–2161) are §8.3's own table of "what this function is at another scale",
  and Task 16's is explicitly *anti*-framing ("nothing anywhere that knows what
  an enzyme is"). Those are fine. Fix the "noble" one.

## Damage

Still implicit, everywhere. No damage field, no age, no hitpoints, no
per-molecule struct in either plan file. Task 15's module header (2013–2017)
states the coupling correctly — decay is O(1) per species *because* damage is
implicit — and Task 21 Step 3 (3123–3129) now instructs an implementer to
rework the measurement before reaching for the provenance fallback, with the
asymptotic-class reason attached. This remains the best-defended invariant in
the plan.

## Call-graph trace from `Beaker::step` (Task 18 Step 4)

| Step | Reaches | Verdict |
|---|---|---|
| 1 `tree.total()` | array read | clean |
| 2 `det_math::ln` | one transcendental | clean |
| 3 `tree.sample(u2)` | tree descent | clean |
| 4 apply channel | count arithmetic; products pre-resolved | clean — canonicalisation genuinely absent |
| 5 recompute catalysis | `catalysis_factor` → `affinity_ordered`, unmemoised | **E3 — §8.6 breach** |
| 5 propensity update | `decay_channels` → `Vec` alloc + `radiogenic_rate` recompute | **E6 — not a breach, a missed edit** |
| 6 novel species | `canonicalise`, `embed`, `signature`; for polymers also `fold` + `cavities` | sanctioned for small molecules; **E4** for polymers |

`fold`, `cavities`, `canonicalise`, `embed` and `signature` are otherwise
unreachable from `step`, which is the fix working. The `FoldCache` (Task 12)
memoises a pure function with eviction that cannot change results, and
`AffinityMemo` is proven equal to direct computation by test (1563–1577) rather
than assumed — both are §8.6-required memoisation, not per-molecule state.

## Checked and found sound

- **Catalysis really does fall out of two cavities and nothing else.** No
  `is_enzyme`, no enzyme type, no species-name branch, no rate bonus keyed on a
  classification. The only inputs are cavity separation, cavity enclosure, and
  two `affinity` scores. Subject to E1, the *shape* of the mechanism is right.
- **RAF detection is detection.** Task 17 sees a reaction list and a food set,
  returns a subset, and nothing anywhere boosts a reaction because it is in a
  RAF. `find_raf` returning `None` means "there is no RAF", which the doc
  proves via closure monotonicity rather than asserting.
- **"Chains of length ≥ 20 must form and persist"** (Task 20 battery, §7.2) is
  a universe-acceptance filter, not a stability rule. Delete it and polymer
  physics is unchanged; only the set of simulated universes changes. This is
  the constant that most *looks* like the bad example and is not one.
- **`decay_scale` is a band parameter**, tuned empirically per §22.6, not an
  outcome. Legitimate.
- **The battery rejects `decay_scale = 0.0`** with §9.4 cited — decay treated as
  a precondition rather than a balancing mechanic.
- **Exit criteria 7 and 9 are now separable**: decay-*off* for criterion 7,
  decay-*equalised* for §15.3, stated explicitly in Task 21 Step 5 and Task 20b
  Step 3. Conflating them was the S12 defect and it is properly closed.
- **The randomised-catalysis control** (Task 20b Step 4) is the right second
  control and asks the question the persistence shadow cannot.
- **Cavity detection's outside-flood** is a genuine enclosure test, and
  `an_open_groove_is_not_reported_as_a_cavity` can actually fail. S8 closed.
- **Task 13's near-wall (min, not max) lining signature in `Span`** is correct
  and is the same operation as §8.2, not a parallel one.
