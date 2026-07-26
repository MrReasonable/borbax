# V0 Plan Audit — 2026-07-26

Four specialist agents reviewed the spec and both plan files **before any code
was written**, and a fifth reviewed dependencies. This is the consolidated
finding list.

Severity is "how expensive is this to discover late", not "how wrong is it".

## Status: all findings applied

Fixed across commits `62afa68`..`HEAD`. Every severe and medium finding below
has been addressed in the spec, both plan files, or `CLAUDE.md`. The two that
changed scope rather than fixing a defect:

- **S12** added exit criterion 9 and Task 20b, because the plan could not
  otherwise reach criterion 7.
- The dependency review found the stated reason for hand-rolling vector maths
  was **factually wrong** — `glam` and `nalgebra` would be bit-identical here,
  verified in their source. The hand-roll stays, but on a size argument. A
  false determinism doctrine would have been used later to reject a dependency
  that is genuinely correct.

Kept as the record of *why* each change was made; the reasoning is not
reconstructable from the diffs alone.

---

## S1 · The binding kernel searches reflections, not rotations

**Found by:** geometry-numerics-reviewer · **Verified** by exhaustive group
computation at D = 12/42/162 · **Affects: spec §8.3, plan Task 10**

Two bodies in contact touch along *opposite* directions: A's extent along `+d`
must meet B's extent along `−d`. The pseudocode in §8.3 and the implementation
in Task 10 both pair `a.r[i]` with `b.r[perm[i]]` — **the antipodal map is
missing.**

Because `−I` is not in the 60-element rotation group (det = −1), the function
as written maximises over the 60 **improper** elements of I_h. That is the
exact inverse of §22.8 and of CLAUDE.md's "searches rotations but not
reflections". Every molecule is scored against its partner's point inversion.

Nothing fails. Binding still happens, scores still vary, chemistry still runs.
What is lost is that *bumps meet hollows* — complementarity is computed between
mirror partners. Any homochirality result reports the wrong handedness, and
cavity–substrate recognition selects the enantiomer of the molecule that fits.

**The existing test cannot catch it.** `complementary_pair` builds its
complement at the *same* index, so it was written to match the implementation.

Fix: build an antipode table in Task 7 (the vertex set is antipodally closed to
1e-31, verified), then index `b` through `anti[perm[i]]`. Symmetry survives.
Replace the test with one that builds the *physical* complement from a real
embedded molecule.

## S2 · Polymers have no species identity, so folding lands in the hot loop

**Found by:** emergence-auditor · **Affects: plan Tasks 14, 18**

`rate()` takes a `catalysis: f64` that nothing constructs. `Reaction` has no
catalyst field. `find_raf` demands a catalyst list nothing builds. The join is
missing, and every natural way to supply it violates §8.6.

Root cause: `Interner` interns `Mol12` only. There is nowhere to *put* a
polymer's fold or cavities, so an implementer computes them at the point of
need — inside `step`. On a cache miss that is 20,000 anneal steps plus a
2,520-iteration `affinity` per cavity×substrate pair. **~4 orders of
magnitude**, not 2.

Fix before Task 14 is written: give `Polymer` a `SpeciesId` from the same
interner; put `fold` and `cavities` in `SpeciesRecord`; add a `catalysis` field
to `Reaction` resolved once at channel creation.

## S3 · `det_math` is in a crate two of its callers cannot depend on

**Found by:** determinism-auditor · **Affects: File Structure, Task 10**

It lives in `borbax-molecule`; `borbax-rng` and `borbax-universe` are *below*
that in the dependency order. The one file meant to isolate every transcendental
is unreachable from two of the three crates that need it. Move it to
`borbax-units`, which has no dependencies. **This is a File Structure change, so
it is cheapest before Task 1.**

## S4 · Solvent attack is sign-inverted and `exposure` is never read

**Found by:** emergence-auditor · **Affects: plan Task 14, spec §9.5**

`affinity` is a negated sum of squares, always ≤ 0. So `-affinity ≥ 0` and is
**largest when complementarity is worst** — molecules the solvent cannot touch
decay fastest. The burial mechanism runs backwards.

Compounding it: `Fold::exposure` — the quantity §9.5 actually specifies — is
computed, tested, cached, rendered, and **read by nothing**.

Exit criterion 7 would still pass. It only checks that decay bites, not which
direction.

## S5 · The folding chain wraps around the grid for n ≥ 60

**Found by:** geometry-numerics-reviewer · **Verified arithmetically** ·
**Affects: plan Task 11**

The extended start runs `centre ± n/2` along one diagonal. At `GRID = 64` a
64-mer already wraps. The flat index stays in range, so there is no panic and
`in_bounds` cannot detect it — it decodes the *wrapped* coordinates, which pass.
The comment "cannot escape the grid even fully extended" is wrong by ~3×.

Task 13's own tests fold 70- and 80-mers. Those chains are geometrically
meaningless — self-avoidance and connectivity evaluated on a sheared torus.

Fix: compact (boustrophedon) start, plus a `debug_assert!(in_bounds)` per
monomer and a property test at `MAX_POLYMER`.

## S6 · Canonicalisation indexes 12-element arrays by raw element id

**Found by:** geometry-numerics-reviewer · **Affects: plan Task 6**

`colour[v] = m.elem[v]`, and element ids reach 119. `masks[colour[v]]` is
`[u16; 12]`; `1u16 << colour[v]` panics in debug and **silently masks the shift
in release**, so refinement terminates before it is stable — two relabellings of
the same molecule intern as different species.

Every test uses `ElementId(0..5)`, so all of Task 6 passes.

## S7 · The embedding's orientation is not stable

**Found by:** geometry-numerics-reviewer · **Measured** over 200 random trees ·
**Affects: plan Task 8**

Shape is stable (aligned RMSD/R_g median 0.13). The *frame* is not: unaligned
displacement median **0.74**, p90 1.34, max 2.59. The Fibonacci init depends on
`n`, so adding one atom reorients the whole configuration by a generic rotation
— and canonicalising over 60 rotations cannot remove a generic SO(3) rotation
(up to ~44° from the nearest icosahedral element).

**The plan's own test passes** on chain(8)→chain(9) at 0.45. Over random trees,
**135 of 200 exceed the threshold.** A textbook case of an example-based test
passing on everything except the graphs that matter.

Consequence: mutation is a random walk in signature space, neutral networks
never appear, and the battery rejects universes for reasons unrelated to the
universe.

Fix: canonicalise the embedding frame from the canonical atom order. Replace the
example with a property test asserting *both* aligned and unaligned bounds.

## S8 · Cavity detection has no enclosure test

**Found by:** geometry-numerics-reviewer · **Affects: plan Task 13**

Candidates are "empty cell with ≥ 8 of 12 neighbours occupied". Nothing tests
whether the region connects to the outside — the flood fill floods the
*candidate* set, not the empty space. A deep surface groove, freely open to
solvent, is indistinguishable from a sealed void.

`cavities_are_enclosed_by_construction` cannot fail: candidates already passed
the threshold, so the mean over them exceeds it by arithmetic.

Fix: flood empty space inward from the bounding box; interior = empty ∧
unreached. Test a hollow shell against the same shell with one monomer removed.

## S9 · The RAF closure produces non-RAFs

**Found by:** alife-researcher · **Affects: plan Task 17**

`missing[ri]` counts only *non-food* reactants, but `Incidence::build` indexes
**every** reactant including food ones. For `A + B → C` with `A ∈ F`: popping
`A` decrements to zero and produces `C` **without `B` ever existing**. The
F-generated condition is violated and the returned set need not be a RAF.

All four tests use single-reactant reactions, so none catches it.

Also: RAF is a *structural* property. `find_raf` sees no counts, so §9.6's death
criterion fires only when a channel disappears, never when a catalyst's count
hits zero.

## S10 · CI is red for twenty tasks, disabling the determinism matrix

**Found by:** determinism-auditor · **Affects: plan Task 1**

The `test` job unconditionally runs `borbax-cli`, created in Task 21. CI fails
on all three platforms from Task 1, so nobody reads it — and `determinism-matrix`,
gated on `needs: test`, never runs at all. Its inertness is accidental.

No task implements the `goldens` subcommand or defines the state-hash format,
which is the mechanism the whole §13.6 contract rests on.

## S11 · The fold cache key omits the universe

**Found by:** determinism-auditor · **Affects: plan Task 12**

The fold depends on the universe through monomer classes and `w_charge`, but
`key_of` hashes the sequence only. A hit from universe A answers for universe B.
Task 20 iterates 40 universes.

## S12 · The plan cannot reach its own exit criteria

**Found by:** alife-researcher · **Affects: spec §15, §23; plan Task 20**

§15's novelty metrics, neutral shadow, and plateau fitting are specified but
**no task builds them**. Exit criterion 7 requires a shadow comparison. Task
20's only novelty test asserts that two enum variants are unequal.

Also: "selection switched off" must mean *equalised* decay rates, not disabled
decay. §23 criterion 7 and Task 21 Step 5 say "decay disabled", which is a
different experiment.

Either add a Task 20b, or move criterion 7 to V1 explicitly.

---

## Medium

- **M1 · Fold contact energy is a second mechanism with the wrong sign.**
  `-fa*fb` favours like-attracts-like, contradicting its own comment, §8.3, and
  Task 5's bond matrix. Three sites, two conventions. Use `-(fa+fb)²`.
  *(emergence-auditor)*
- **M2 · The anneal is not Metropolis.** Uphill acceptance is independent of ΔE,
  so a catastrophic move is as likely as a marginal one. The fold becomes an RNG
  artefact rather than a function of the sequence — which destroys the
  many-to-one folding map the evolvability argument rests on.
  *(geometry-numerics-reviewer)*
- **M3 · `PERIOD_LENGTHS = [2,6,8,10,14,18]` is real shell structure** — real
  period lengths ∪ real subshell capacities — plus four real periodic trends in
  the real direction. G3 says period lengths are generated "rather than
  following any real shell structure". `xtask` cannot see it: it scans for data
  *files* and this is a `const`. *(emergence-auditor)*
- **M4 · The `1.0e3` in `catalysis_factor` decides exit criterion 6.** Every
  other magnitude is generated per universe; this one is a bare literal, and the
  test asserting `f > 2.0` is calibrated to it. `BIND_THRESHOLD = -12.0` scales
  with D — so the §22.2 resolution sweep would measure the threshold moving,
  not the resolution. *(emergence-auditor)*
- **M5 · Cavity lining signature measures the far wall, in grid units** while
  molecule signatures are in `Span`. `affinity(cavity, molecule)` adds
  incommensurable numbers. *(geometry-numerics-reviewer)*
- **M6 · Three unrouted transcendentals** — `next_normal`'s `ln`/`cos`,
  `powf` in element generation, `cos`/`sin` in the embedding init. Each feeds a
  *discrete* decision downstream, so 1 ULP flips a rotation argmax or a
  Metropolis accept. Not "last digit differs" — different species, different
  folds. *(determinism-auditor, corroborated by geometry)*
- **M7 · Bedau class mislabelled.** The failure mode described is **class 3b**
  (bounded diversity, unbounded activity per component), not class 2 (bounded
  activity). Appears in spec §2.7, CLAUDE.md, and a Task 20 test comment. The
  countermeasure is right; the label is wrong. *(alife-researcher)*
- **M8 · Gillespie: `u1` can be exactly 0.0** → `ln(0) = -∞` → the clock jumps to
  infinity. Also, same-species bimolecular propensity is `c·n(n-1)/2`, not
  `c·n²/2` — invisible at large n, wrong exactly where a RAF nucleates.
  *(alife-researcher)*
- **M9 · `Stream` derives `Copy`**, contradicting §13.1 and CLAUDE.md. The plan
  already contains a duplicate-stream instance: the shell pattern is generated
  twice from a re-derived stream. *(determinism-auditor)*
- **M10 · `may_bind` filters nothing** — rejects 0.000% of 20,000 random pairs,
  and its bound is not a bound on anything the full search computes. Replace
  with a genuine Cauchy–Schwarz bound. *(geometry-numerics-reviewer)*

## Low / mechanical

- `xtask` greps for `REAL_ELEMENT_SYMBOLS`; Task 4 names it `REAL_SYMBOLS`. CI
  goes red at Task 4 and the tempting fix is to weaken the check.
- Task 7 Step 3 does not compile — E0502, verified with rustc. Hoist the
  midpoint into a `let`.
- `bind_probability` caps at 0.5, since `affinity ≤ 0` always.
- `ideal_gap`, `Signature.r[]` and `Element.stability` are bare `f64` while
  being typed quantities — a quiet G4 defeat.
- `radiogenic_rate` and `decay_channels` allocate per call on a §8.6 path.
- Four accumulation sites deterministic but not commented as load-bearing;
  `Fold::energy` and `contact_energy` are two copies of the same loop.
- `Task 19`'s `depth_sort` will be the project's first float-keyed sort *with a
  payload* — needs an id tie-break before it exists.
- `SEARCH_LEAF_CAP` hits return a partial minimum silently.

---

## What the reviewers checked and found sound

- **The rotation table is correct**, verified exhaustively at all three
  resolutions: 60 distinct permutations, closed under composition (3600
  compositions checked per resolution), inverse-closed, identity present once,
  every element det = +1. The 1e-12 tolerance has 19 orders of margin.
- **Colour refinement's rank computation is correct** and isomorphism-invariant;
  refinement is monotone so the stopping test is valid.
- **The segment tree** is the best-defended determinism decision in the plan.
- **Damage stays implicit** — no per-molecule struct, no damage field, no age.
  The plan defends this invariant better than any other.
- **No hardcoded biology anywhere.** No code branches on a biological name.
- **The blocklist is exclusion-only.** G2 clean.
- **Exactly one `HashMap`**, and it is verified never iterated.
- **Fixed-point `Mass` quantisation is sound** — error enters once, at
  generation; molecular mass is an exact integer sum.
- **FCC neighbours are exactly the twelve `(±1,±1,0)` permutations**, all
  parity-preserving.
- **SMACOF's update is the standard Guttman transform**, and 240 iterations is
  adequate (stress changes 0.06% over the last 10).
