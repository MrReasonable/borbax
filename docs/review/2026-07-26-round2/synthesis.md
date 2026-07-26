# Review synthesis — Borbax V0 plan, post-fix round (2026-07-26)

Six specialist reviews reconciled. This is a transformation of their output plus
gap detection; **no finding here is mine** except where explicitly marked
*[coordinator]*, and those are questions, not findings.

---

## 0. Process failure to declare up front

**The cross-check round did not happen.** My tool set for this task was
`Read`, `Write`, `Bash` only — no `Agent`/`Task`, no `WebSearch`, no `WebFetch`.
The five cross-checks you asked for could not be dispatched. Per protocol I did
not block; I marked them and moved on.

I compensated where the question was factual rather than a judgement, using
`Bash` against the plan files and a Python reproduction. That turned out to be
worth more than expected: **five findings that were slated for cross-check are
now verified arithmetically**, which is a stronger result than a specialist
verdict. What remains open is the genuinely judgement-shaped half, and I say so
per finding.

`cross-check-failed` below means "no independent verdict obtained", not
"doubted".

**This is a bug in the calling skill, twice over:**
1. The coordinator was dispatched without `Agent` while being told it had one.
2. `alife-researcher` was dispatched without `Write` and its findings were
   transcribed by the parent (its own file says so, line 2:
   *"Fix its `tools:` line"*). See §5 for what that cost.

---

## 1. Counts

At the grain of *distinct issue after dedup* (raw finding count across the six
files is ~100; many are the same defect seen from different angles).

| Classification | Count | Note |
|---|---|---|
| `corroborated` (2+ reviewers) | **16** | C1–C16 |
| `disputed` | **6** | D-1 … D-6 — read these first |
| `solo`, verified by coordinator | **5** | G3, E2, G10, D2, alr-011 premise |
| `solo`, cross-check-failed | **~28** | tagged `[single-reviewer]` |
| `gap` (domain in scope, nothing produced) | **4** | §5 |
| Reviewer self-retraction | **1** | geometry withdrew the tighter `may_bind` bound (G8) — dropped, noted here as required |

**Contradiction density: 6 genuine conflicts. This is at the Step-5 escalation
threshold.** See §7 — but the meta-signal here is not "the artifact is
internally inconsistent"; it is more specific and worse than that, and I set it
out in §7 rather than §1.

---

## 2. DISPUTED — read this section first

These are where a human is genuinely needed. Everything else is a work item.

### D-1 · Task 13 `lining_signature` — endorsed by one reviewer, called silently wrong by another

`chemistry.md:1316–1368`

- **geometry G4**: *"`nearest` is a minimum over a hemisphere, so it is ≈0
  everywhere… the previous review's `max → min` change swapped one wrong
  quantity (the far wall) for another (the equator)."*
- **geometry G7**: `LATTICE_SPAN = SQRT_2` *"converts in the wrong direction and
  ignores the universe."*
- **emergence, "Checked and found sound"**: *"Task 13's near-wall (min, not max)
  lining signature in `Span` is correct and is the same operation as §8.2, not a
  parallel one."*

Emergence endorsed **precisely the two things geometry says are wrong** — the
`min`, and the `Span` conversion. This is a direct contradiction of premise, not
two reviewers looking at different properties.

**Two other reviewers touched the same code and are NOT in conflict, which the
reader must not misread as a third and fourth endorsement:**

- **rust-developer-expert #9**: *"`reach_of` reads well… The two-pass structure
  is right."* This is a **style** verdict. The same section then flags an
  inversion in the same expression — *"`0.0` means 'wall touching this
  direction', the tightest possible fit, but the case being encoded is 'no
  lining monomer that way', i.e. open space"* — and **explicitly defers**: *"the
  right physical value is geometry-numerics' call."*
- **determinism-auditor**: *"summation order is pinned and commented… This one
  is right."* A **determinism** verdict, scoped to accumulation order.

So: the two-pass *structure* is endorsed by three reviewers and disputed by
none. The *quantity it computes* is disputed 1-vs-1, with rust-dev's independent
observation of an inversion in the same expression landing on geometry's side of
the ledger without claiming to.

**Coordinator note:** I verified G7 arithmetically. `LATTICE_SPAN` is applied as
a **multiplier on integer coordinate differences** (`chemistry.md:1354–1356`).
FCC nearest neighbours are √2 integer units apart, so multiplying by √2 sets the
nearest-neighbour distance to **2.0 Span**, not √2. It is also a compile-time
constant while molecule extents scale with generated radii (`radius_base ∈
[0.8,1.4]`, `v0.md:1660`, reaching ~5.0 after the period/group factors). G7 is
correct as arithmetic. That does not settle G4, which is the harder claim.

**What settles it:** geometry's own proposed test — a synthetic ellipsoidal
cavity, semi-axes 4/2/2, assert `sig.r` is ~2× larger along the long axis. One
run. Do that before anyone edits this function again.

---

### D-2 · Is S8 (cavity enclosure) closed?

- **emergence, "Checked and found sound"**: *"Cavity detection's outside-flood
  is a genuine enclosure test, and `an_open_groove_is_not_reported_as_a_cavity`
  can actually fail. **S8 closed.**"*
- **geometry G6**: *"the enclosure gate now fragments the cavities the flood
  found, and the flood only finds sealed voids"* — `MIN_ENCLOSURE` applied
  after a real flood *"drops the interior cells of any void bigger than ~6
  cells… which splits one genuine cavity into several disconnected shells."*
  Directly threatens criterion 6: *"one void split in three reads as an enzyme."*
- **rust-developer-expert #9**: the new enclosure test passes
  `remove_one_monomer(&sealed)` together with the **original** `p`. *"If that
  helper drops a coordinate, the test is measuring a fold/polymer mismatch, not
  an open groove."*

**2 against 1**, and the two dissenters found *different* defects — one in the
mechanism, one in the test that is supposed to protect it. Emergence's "S8
closed" rests on the test being able to fail; rust-dev's point is that it can
fail *for the wrong reason*.

**Verdict:** treat S8 as **not closed**. Emergence's claim is about the flood
being a genuine enclosure test, which is true and is a real improvement; it does
not cover what the gate downstream of the flood then does to the result.

---

### D-3 · `.powi(` — same symptom, opposite remedies. **[coordinator ruling]**

- **determinism D5**: `.powi(` is a **false positive** in `BANNED_CALLS`; keep
  `powi`, drop it from the list. *"CI goes red at Task 4 and stays red. The
  tempting fix is to weaken the check — the same failure mode as S10."*
- **rust-developer-expert #6**: the three sites *"will fail the plan's own CI
  gate"*; rewrite them as explicit multiplication, citing `det_math`'s own doc
  (`v0.md:892–894`) that *"`powi`'s multiply tree is chosen by LLVM, so it can
  be re-associated by a compiler upgrade without our source changing."*
- **emergence E7** leans determinism's way: *"`powi` is repeated multiplication
  so it is exactly specified."*

This is factual, not a priority conflict, so precedence does not apply. I settled
it. There are **eight** sites, not three (`v0.md:1722, 3300×3, 3332×3, 3884×2`;
`chemistry.md:2229–2231×3`):

- **Seven of eight are `.powi(2)`.** LLVM's `ExpandPowI` emits `x*x` for
  exponent 2 with *no* reassociation freedom. There is nothing for a compiler
  upgrade to change.
- **One is `.powi(6)`** (`v0.md:1722`, `stability`) and genuinely has tree-shape
  freedom.
- The tree shape is fixed by the **rustc pin**, which `CLAUDE.md` already names
  as part of the determinism contract — *"Bumping it invalidates golden hashes
  and requires deliberate regeneration."* So even the `powi(6)` risk is
  already covered by the mechanism that exists for exactly this.

**Ruling: do both, they are not alternatives.** Drop `.powi(` from
`BANNED_CALLS` (determinism D5 is right that a permanently-red gate is the worst
available outcome — it is S10's failure mode returning), **and** write out
`v0.md:1722` explicitly as rust-dev proposes, because it feeds a discrete
propensity channel and costs one line. Then `det_math`'s doc comment must be
corrected — as written it argues for a rule the plan does not follow, which is
what generated the disagreement.

---

### D-4 · `canonicalise_frame` — repair it, or put it on trial? **[precedence ruling]**

- **geometry G2**: measured **worse than not having it**. Median unaligned
  displacement 0.87 with the function, **0.64 with no frame canonicalisation at
  all**. *"I do not have a validated [fix], and I would rather say that than
  invent one."*
- **rust-developer-expert #4**: the collinear early return *"is dishonest"* —
  proposes completing the frame from the standard basis.
- **determinism**: *"Determinism clean — its stability under near-degeneracy is
  a question for the geometry reviewer, not me."* (correctly scoped, no conflict)
- **rust-performance-expert**: *"Fine, and by a wide margin… 0.03% of `embed`.
  No action, and no need to revisit."* (cost only, no conflict)

Nobody contradicts geometry's measurement. But rust-dev proposes *repairing* a
function geometry measured as net-negative, and two other reviewers signed off
on it for properties that are not the one in question. A reader skimming four
sign-offs will conclude the function is fine.

**Applying `CLAUDE.md` Review precedence: rank 2 (correctness of the physics)
over rank 6 (idiom and readability).** rust-dev #4 is *correct* — the early
return is dishonest and the `(array, n)` signature is unenforceable — and should
be applied **if the function survives**. It should not be applied first, because
polishing the collinear branch of a function that loses to its own absence is
work aimed at the wrong question.

**This is the single most important open item in the review** and I have set out
why in §7.

---

### D-5 · `MIN_ENCLOSURE = 8` — filter or defect?

- **emergence E10**: *"`MIN_ENCLOSURE = 8` survives, and is now a detector…
  It is a reporting filter."* Suggests dropping it to ~5. *"Confidence: medium —
  this is a judgement call, not a defect."*
- **geometry G6(1)**: *"now redundant and harmful"* — as a **membership** filter
  after a real flood it fragments genuine cavities. Silently wrong.

The **fixes converge** (both want it demoted or loosened); the **severity does
not**. Emergence read it as a threshold that reports fewer cavities; geometry
read it as one that splits single cavities into several. Those are different
failure modes and only geometry's threatens criterion 6.

**Verdict:** geometry's is the stronger reading and is arithmetic — a cell in
the middle of a void has void neighbours, so it cannot pass a ≥8-of-12 occupied
test. Both reviewers' fix is the same: demote it to a per-component summary, not
a membership gate. Emergence's additional point — that `enclosure` is *already*
used continuously in `catalysis_factor`, so the hard cutoff applies the same
quantity twice, once smoothly and once as a cliff — stands independently and is
good.

---

### D-6 · `may_bind` — two sign-offs and one measurement saying it rejects everything

- **rust-developer-expert #14**: *"a genuine bound with a stated property… a
  real improvement on M10."*
- **rust-performance-expert**: *"its correctness is settled"* — dropped the
  prefilter benchmark from its list.
- **geometry G8**: *"The bound is valid, verified… 12,000 random pairs across
  three thresholds: 0 violations."* **And then:** *"Given G1, it rejects
  everything. On 1,770 pairs of real embedded 8-atom signatures it rejected
  **100%** at thresholds −40, −100, −300 and −1000."* Plus: the test's own
  comment says it *"must actually reject something"* and the assertion below it
  checks that `may_bind` **accepts** a perfect complement — M10's actual
  complaint is not retested.

Not strictly a contradiction — geometry agrees the bound is sound and is
strictly better informed about the rate. But **two reviewers closed M10 on a
filter a third measured as a "nothing binds" switch**, and rust-perf dropped the
benchmark that would have caught it. Flagging because the reader will otherwise
see M10 as closed.

**Contingent on G1.** If G1 is right, the 100% rejection is a symptom of G1, not
of `may_bind`. Fix G1 first, then re-measure. Add geometry's rate assertion
regardless — it is the assertion the test's own comment promises.

---

## 3. Findings, severity-grouped

Tags: `[corroborated]` · `[disputed]` · `[single-reviewer]` · `[measured]`
(reproduced by the reviewer or by me) · `[coordinator-verified]`

### 3.1 BLOCKING — physics silently wrong (precedence rank 2)

These compile, pass their tests, and are wrong forever. They change *what the
code is*, so they must be settled before the code that depends on them is
written.

**P-1 · `ideal_gap` is on the wrong length scale — binding is 93% size, 7% shape**
`[single-reviewer]` `[measured]` `cross-check-failed`
geometry G1 · `v0.md:4234`, `v0.md:2046`

> *"The 2,520-iteration rotation search — the mechanism the entire design rests
> on — moves the score by 7% of its magnitude. 93% of `affinity(A,B)` is 'how
> big are these two molecules'."*

Rotation contrast `(max−min)/|mean|`: **0.069** with fixed `ideal_gap = 2.0`,
**1.066** pair-scaled. *"Homochirality, cavity selectivity and neutral networks
all degrade to noise on top of a size ranking. §23 criterion 6 reads as
'catalysis did not emerge' and the decay band gets blamed."*

**[coordinator-verified — the structural half.]** `ideal_gap` is drawn from
`[1.6, 2.8]` and is a fixed per-universe scalar (`v0.md:2046`). `Signature.r` is
documented as *"the support function of the atom spheres"* (`v0.md:3812`), so it
is `pos·d + radius` and grows with molecular extent; generated radii reach ~5.0.
Comparing a size-independent constant against a size-dependent support function
is a genuine structural mismatch, independent of the exact contrast numbers,
which I could not reproduce without the full embedding.

**Open, and it is the question the cross-check was for:** is `pair_gap`
legitimate physics, or does it launder the answer in? *[coordinator]* — my read,
offered as a prior for the specialist and not as a verdict: `pair_gap` is built
from `mean(a.r)` and `mean(b.r)`, both **rotation-invariant**, so it cannot
encode a preferred orientation, which is the laundering risk. It changes binding
from "absolute size match" to "relative shape match", which is what §8.3 says
binding *is*. I lean legitimate. **This needs emergence-auditor's verdict, not
mine.**

---

**P-2 · `canonicalise_frame` measures worse than doing nothing; signature space
has no locality** `[disputed]` (D-4) `[measured]` `cross-check-failed`
geometry G2 · `v0.md:3604`

| frame | unaligned median | fraction over the plan's 0.6 bound |
|---|---|---|
| the plan's `canonicalise_frame` | **0.87** | 69.5% |
| **no frame canonicalisation at all** | **0.64** | 54.5% |

End to end, over 120 trials: *"in 32% of cases the mutant is FURTHER from its
parent than an unrelated molecule is… **signature space has no locality.**
Mutation is a random walk, neutral networks cannot form, and §8.4's evolvability
argument has nothing under it."*

Root cause 3 is the one that matters and is why none of the reviewer's five
alternative frames helped: *"Canonical atom order is not continuous under graph
edits. Adding one atom can permute the whole canonical labelling, so **any**
frame anchored to 'the first atom in canonical order' inherits a discontinuity
that no better selection rule removes."*

**No validated fix exists and the reviewer says so.** The plan's own remediation
note (*"Raise `ITERATIONS` or revisit the weighting first"*) is called
misleading — aligned RMSD is already 0.12, so iterations are not the constraint.

**This is not a defect in a task. It is a possible falsification of a design
premise.** See §7.

---

**P-3 · Contact energy is sign-inverted under a minimising annealer**
`[single-reviewer]` `[coordinator-verified]` `cross-check-failed`
emergence E2 · `chemistry.md:608`

**Verified arithmetically by me.** `fa`, `fb` are class representatives in
[−1,1]; `w_charge ∈ [0.6,1.4]`, always positive; the annealer accepts on
`delta <= 0.0`, i.e. it **minimises**:

| fa | fb | `-(fa+fb)²·w` | minimiser's preference |
|---|---|---|---|
| +1 | −1 (complementary) | **0.00** | worst |
| +1 | +1 (like) | **−4.00** | **best** |
| −1 | −1 (like) | **−4.00** | **best** |

Like-attracts-like — *"the exact defect M1 raised, reintroduced with the
opposite sign error."* The comment above it correctly diagnoses the old
`-fa*fb` bug and then reproduces it. **Introduced by fix pass 1.**

**[coordinator] — a caveat on the proposed fix that a cross-check would likely
have raised.** Emergence's `((fa+fb)²−4)·w` correctly inverts the ranking, but
it ranks on `|fa+fb|` alone. Tabulating it:

| fa | fb | `((fa+fb)²−4)·w` |
|---|---|---|
| +1 | −1 (perfect complement) | −4.00 |
| −0.5 | +0.5 (complement) | −4.00 |
| **0.0** | **0.0 (both neutral)** | **−4.00** |
| +1 | +1 (like) | 0.00 |

**A contact between two neutral monomers scores identically to a perfect
complement.** So a bland homopolymer folds as tightly as an alternating charged
one, which weakens the sequence→fold map that §8.4's evolvability argument
needs. This is inherited from §8.3's own `−(a_A+a_B)²`, so it is arguably
consistent by construction — but at monomer scale there is no shape term to
break the degeneracy, which is exactly the difference the comment at
`chemistry.md:605–607` relies on. **Worth one question to geometry-numerics
before applying.** The fix is still clearly better than what is there.

---

**P-4 · Cavity–substrate fit is inverted: a bigger cavity prefers a smaller substrate**
`[corroborated]` (geometry G5 + rust-dev #9, who found the same inversion and
deferred on the physics) · `chemistry.md:1361–1368` → §8.3

> *"`d(preferred r_sub)/dh = −1`: a roomier cavity selects a smaller substrate.
> Containment is a difference (`h − r_sub = clearance`); §8.3's form is a sum,
> because there the two extents are measured from two **different** centres."*

Also: *"the kernel pairs slot `i` with the substrate's antipodal slot, but for a
body sitting **inside** a cavity both are measured from the same centre and want
the same direction."*

rust-dev #9 independently: *"`0.0` means 'wall touching this direction', the
tightest possible fit, but the case being encoded is 'no lining monomer that
way', i.e. open space. That inversion feeds straight into
`a.r[i] + b.r[j] − ideal_gap`."*

Geometry's own confidence: *"high on the monotonicity argument (it is one line
of algebra); medium on my proposed encoding being the best repair."* Take the
finding, treat the fix as a proposal, and let the test decide.

---

**P-5 · The SMACOF target matrix is not a metric, and bond order is chemically inert**
`[single-reviewer]` `[coordinator-verified]` — **not covered by the previous
review round**
geometry G10 · `v0.md:3492`

Two defects in one expression.

*(a)* `target[i][j] = hops[i][j] * (radius(i) + radius(j))` uses the
**endpoints'** radii for a multi-hop distance. **Verified by me:** for a chain
i–j–k with `r_i = r_k = R` and a small `r_j = ρ`, the triangle inequality is
violated by `2R+2R` vs `2(R+ρ)` — a factor approaching **exactly 2× as ρ→0**,
and 1.92× at the plan's actual radius extremes (R=5.0, ρ=0.2). SMACOF converges
against whatever you give it, so this is silent. *"Shape is driven by which two
atoms are at the ends rather than by what the path is made of."*

*(b)* **Verified by me:** `hop_distances` calls `m.neighbours()`, which is
`self.adj.iter().fold(0, |acc, plane| acc | plane[a])` (`v0.md:2288`) — it ORs
all three bond-order planes together. **A single bond and a double bond between
the same atoms produce an identical embedding, hence identical signature, hence
identical binding.** Bond order changes `CanonForm.planes` (species identity)
and `valence_used`, and nothing else. A generated degree of freedom that the
chemistry cannot see.

Geometry's Floyd-Warshall replacement fixes both at once and is 8 lines for
n ≤ 12. *(b)* is arguably an **emergence** finding as much as a geometry one and
was not seen by that reviewer.

---

### 3.2 BLOCKING — the plan does not compile or CI is red

Nine defects here. **Every one is caught by `cargo check` or `cargo clippy` in
seconds.** That fact is the whole of §7.

| # | Defect | Reporters | Compiler says |
|---|---|---|---|
| B-1 | Task 16 body uses `BIND_THRESHOLD`, `IDEAL_SEPARATION`, `MAX_ENHANCEMENT`, `* 1.0e3` — none defined; the header says all four were replaced | emergence E1, rust-dev #2, rust-perf P4, geometry G11 — **4 reporters** `[corroborated]` | `E0425` |
| B-2 | Task 14 `intern` pushes `SpeciesRecord { canon, .. }`; `canon` is not a field and `kind`/`radiogenic_rate`/`fold`/`cavities` are missing | rust-dev #1, emergence E5, rust-perf P1 — **3** `[corroborated]` | `E0560`, `E0063` |
| B-3 | `intern_polymer` called (`chemistry.md:1551–1552`), defined nowhere | rust-dev #1, emergence E5 `[corroborated]` | `E0599` |
| B-4 | `FoldId` occurs exactly once — as a field type at `chemistry.md:1663`. Never defined | rust-dev #1 `[single-reviewer]` | `E0412` |
| B-5 | `energy_scale` (`chemistry.md:519`) undefined; sole occurrence in either file | determinism D7, geometry G11, rust-dev #5 — **3** `[corroborated]` | `E0425` |
| B-6 | `fold_seed` defined as `FoldWorkspace::fold_seed`, called as a free fn (`chemistry.md:418`) | determinism, rust-dev #5 `[corroborated]` | `E0425` |
| B-7 | `local_contacts`' `coords` param never read → `unused_variables` → `-D warnings` | rust-dev #5 `[single-reviewer]` | CI red |
| B-8 | `recount_contacts`: `(k * 0) + class_of(p,k)` trips `clippy::erasing_op` (**deny by default**) | rust-dev #5 `[single-reviewer]` | CI red |
| B-9 | `.powi(` in `BANNED_CALLS` against 8 library sites | determinism D5, rust-dev #6, emergence E7 — **3** `[disputed]` (D-3) | xtask red from Task 4 |

Plus, from Task 8: *"Undefined test helpers — `random_tree`, `permute`, `dist`
and `kabsch_rmsd` are all used and none is provided"* (geometry G11), and
`//!` inner doc after `use` at `chemistry.md:2169` will not compile
(rust-perf P4).

**B-1 verified by me in one grep.** `catalytic_prefactor` appears at exactly one
line in either plan file — `chemistry.md:2190`, inside the doc comment claiming
it exists. `UniverseConsts` (`v0.md:1992–2010`) has nine fields and is not one
of them. `BIND_FRACTION` is declared at `:2195` and never read.

> rust-dev #2: *"A false 'already fixed' comment is worse than no comment."*

**B-2 verified by me.** `SpeciesKind` (`chemistry.md:1634`) is declared and
never constructed; `Interner::ids` is `BTreeMap<CanonForm, SpeciesId>` and a
polymer has no `CanonForm` — *"there is no identity path for
`SpeciesKind::Polymer` at all, which is precisely the hole S2 was raised to
close"* (rust-dev #1).

---

### 3.3 SEVERE — determinism (precedence rank 3)

**D-a · A second `det_math` calling the platform libm, exempted from the one
gate that would catch it** `[corroborated]` `[coordinator-verified]`
determinism D1 + rust-dev #7 · `v0.md:4323–4338`, `4367–4370`;
`chemistry.md:1898, 2046, 2067, 3095`

Verified by me, all six sites. Task 10 Step 3 still instructs *"create
`crates/borbax-molecule/src/det_math.rs`"* with a platform-delegating body and a
comment referencing **Task 25, which does not exist** in a 21-task plan — while
the code block 6 lines above already calls `borbax_units::det_math::exp`
(`v0.md:4319`). Tasks 14 and 15 then consume the shadow module at three sites:
the Arrhenius factor and **both** thermal-decay propensities.

*"Every reaction rate and every decay propensity… One ULP flips a channel
selection and the trajectory diverges completely."*

The compounding fact, found by both reviewers: `xtask`'s exemption is
**by filename anywhere in the tree** (`v0.md:949`), so the shadow module is
exempt from the check designed to catch it. Both propose the identical fix —
exempt by path.

---

**D-b · `libm` is a caret range and `Cargo.lock` is never committed**
`[single-reviewer]` `[coordinator-verified]`
determinism D2 · `v0.md:181–182`, `v0.md:511`

Verified by me. `[workspace.dependencies.libm] version = "0.2"` is `^0.2`. Task
1 Step 7's `git add` is an explicit file list (`git add .prototools Cargo.toml
rustfmt.toml .github xtask`) and `Cargo.lock` is not in it. **`.gitignore` does
not exclude it** — so adding it is genuinely a one-word change, as claimed.

*"Two clones of the same commit at different times resolve different `libm`
patch versions and produce different goldens."* This crate **is** the
transcendental implementation. `chemistry.md:3189` says to pin it exactly; the
manifest does not.

---

**D-c · CI will report a false determinism failure on Windows, on day one**
`[single-reviewer]` · determinism D3 · `v0.md:452–454`

PowerShell's `>` is `Out-File`: it rewrites the stream with CRLF. *"The first
thing the matrix ever does is fire falsely. The predictable response is
`diff --strip-trailing-cr` or `|| true`, and the gate is dead from that day."*
Fix: `shell: bash`, or better `--out <path>` so no shell touches the stream.

---

**D-d · Task 21 Step 6 is stale and contradicts Task 2** `[corroborated]`
determinism D4 + rust-dev #7 · `chemistry.md:3162–3198`

*"It tells the implementer that the platform libm is acceptable for Tasks 2–20
and that goldens are knowingly non-portable until the last task. Combined with
D1 that is a coherent (and wrong) instruction set which survived the fix."*

---

**D-e · `BANNED_CALLS` is wrong in both directions** `[single-reviewer]`
determinism D5 · `v0.md:939–940`

False negatives: `.log(` misses `.log10(`/`.log2(`; `.exp()` misses `.exp2()`;
entirely absent — `.atan2(`, `.cbrt(`, `.hypot(`, `.asin(`, `.acos(`, `.atan(`,
`.sinh(`, `.cosh(`, `.tanh(`, `.exp_m1(`, `.ln_1p(`, `.to_degrees(`,
`.to_radians(`, `.mul_add(`, and every free-function form (`f64::exp(x)`,
`use libm::*`). `.cbrt(` matters specifically — determinism's own `libm` source
audit found `cbrt.rs` is one of only two files reaching the arch-dispatch path.

---

**D-f · Two propensity-update paths, selected at runtime, never proven equal**
`[single-reviewer]` · determinism D8 · `chemistry.md:2727–2731`

*"A run crosses 256 channels mid-trajectory and silently switches physics. Two
runs of the same seed at different `FLAT_SCAN_LIMIT` values disagree; the CI
matrix cannot see it because all three platforms take the same path."* Fix is a
`debug_assertions` bit-comparison plus one release test at `FLAT_SCAN_LIMIT`
forced to 0 and `usize::MAX`.

---

**D-g · The golden state-hash digest is still unnamed** `[single-reviewer]`
determinism D9 · `v0.md:494–506`

*"The only hashing crate `CLAUDE.md` names is `rustc-hash`, whose algorithm
**changed between 1.x and 2.0** — an implementer following the doctrine picks
the one crate that reintroduces the problem."* Also unspecified: newline/encoding
(which is what makes D-c bite), and what happens to a NaN count.

---

**D-h · Fold cache keyed on `u.seed` while the battery mutates `u.consts` in
place; §13.1's `config_hash` appears nowhere** `[single-reviewer]`
determinism D6 · `chemistry.md:922–933` vs `2886, 2893, 2900, 3120`

Latent today (nothing sweeps a fold-relevant constant), *"one line of a future
sweep away from re-becoming S11."* Separately and definitely: §13.1's
reproducibility tuple is `(universe_seed, world_seed, config_hash)` and
**no `config_hash` exists in either plan file.**

---

### 3.4 STRUCTURAL / SCOPE — interfaces that need to exist before the tasks that use them

**S-a · Volume `V` is not a parameter anywhere** `[single-reviewer]`
`[coordinator-verified]` `cross-check-failed`
alife alr-011 · Task 18

Verified by me: `grep -in 'volume'` across both plan files returns three prose
mentions of *"one well-mixed volume"* and no parameter. `rate_prefactor`
(`v0.md:2000, 2049`, drawn from `[1e3, 1e5]`) is the only prefactor.

> *"Gillespie's `c = k/V` (or `2k/V` same-species), so beaker volume is silently
> baked into `rate_prefactor`. Consequence: §2.7's Moreno-Ofria 'region size is
> a complexity budget' countermeasure is NOT implementable in V0 — a rejected
> universe cannot be distinguished from a beaker too small to nucleate."*

*[coordinator]* — the cheapness question the cross-check was for, answered as
best I can without a specialist: `c` is resolved **once per channel at channel
creation**, not per step, so a division by `V` adds nothing to the step loop. It
belongs on the beaker config rather than `UniverseConsts` (it is a beaker
property, not a universe property), which also keeps it out of the fold cache
key. No goldens exist yet. **My assessment: genuinely cheap now, expensive after
Task 18 lands, and the diagnostic value is high** — without it, "this universe
is barren" and "this beaker is too small" are the same observation. Mark this
as my reasoning, not rust-performance-expert's verdict.

---

**S-b · Task 20b Steps 1 and 3 are presented as independent and are not**
`[corroborated]` · emergence E9(a) + alife alr-012/031/032

The activity threshold *"must come from the shadow cross-over (Rechtsteiner &
Bedau 1999 — gave 50 vs the arbitrary 10 used previously). Step 1 cannot compute
`diversity` or `new_activity` without Step 3's output. As written an implementer
hard-codes a constant, which is the practice the paper was written to replace."*

Emergence puts it in emergence terms: *"If someone picks a literal, it becomes a
knob tunable until the graph looks good."* Two reviewers, two framings, same
defect. alr-032: *"Nothing computes the shadow's own activity distribution,
which is the object Rechtsteiner & Bedau's method consumes. Missing ~5 lines."*

---

**S-c · Catalysis recompute is on the hot path and unmemoised** `[single-reviewer]`
emergence E3 · `chemistry.md:2725–2731` — **§8.6 breach**

Three problems, of which the third is decisive: *"`catalysis_factor` runs
`affinity_ordered(&cavity.signature, substrate, …)` — 2,520 iterations at D=42 —
per ordered cavity pair per channel. `AffinityMemo` is keyed on
`(SpeciesId, SpeciesId)`; a cavity signature has no `SpeciesId`, so the one
affinity call now reachable from `step` is the one affinity call the memo does
not cover."*

Also: *"It is most frequent exactly where it matters least affordably. At small
counts — the regime in which a RAF nucleates — a catalyst present in one or two
copies crosses zero on most firings that touch it."*

And it contradicts Task 14's own `Reaction::catalysis` doc, which says the value
is *"resolved once when the channel is created"* — two statements, two lifetimes.

**Precedence note *[coordinator]*:** rust-performance-expert **P4** touches the
same function and proposes hoisting the loop-invariant `affinity_ordered` out of
the inner loop (O(C²)→O(C)). That is correct and is rank 5. E3 is rank 4 and its
fix — compute at **catalyst-intern time** — subsumes P4's entirely. **Apply E3's
fix; P4's optimisation then comes for free and P4 should not be applied
separately**, or someone will believe E3 is handled.

---

**S-d · Interning a polymer is not a rare event** `[single-reviewer]`
emergence E4 · `chemistry.md:2732–2733`

*"§8.6 sanctions this for canonicalisation because a novel species is rare. It
is not rare for polymers: every condensation that extends a chain produces a
sequence never seen before, by construction, so the Task 12 fold cache misses
**by construction** on exactly this path."* Each intern is ~20k Metropolis steps
plus a cavity flood plus a D-direction lining signature — inside `step`.
*"The S2 fix removed the repeat cost. It did not bound the first-sight cost."*

---

**S-e · `find_raf` returns `Option<RafSet>`, so §15.2's "RAF set count" is 0 or 1**
`[single-reviewer]` · alife alr-026

*"The literature's answer is the count of **irreducible** RAFs. One irrRAF is
polynomial (~15 lines on `find_raf`); smallest is NP-hard. Also bears on §9.6:
maxRAF can shrink drastically while a core irrRAF still closes, so 'the RAF
stopped closing' is a blunt death criterion."*

---

### 3.5 MEDIUM — corroborated

- **M-a `[corroborated]`** (emergence E6, rust-perf P1, determinism) ·
  `radiogenic_rate` is both a `SpeciesRecord` field and a free function
  recomputing it; `decay_channels` allocates a `Vec` per call on a path taken
  every count change; both still read `record.canon`, a field that no longer
  exists. *"Two sources of the same number is how they drift."* rust-perf prices
  it: *"~4×10⁵ allocations/sec/core — 2–4% of the entire budget burned on an
  allocation whose result is three tuples."*
- **M-b `[corroborated]`** (determinism D11(d), rust-dev #12) ·
  `levenberg-marquardt` pulls **`nalgebra ^0.34`** into `borbax-beaker`'s
  runtime dependency graph. Both propose the same fix — a separate
  `borbax-analysis` / `borbax-metrics` crate — and rust-dev states the test:
  *"`cargo tree -p borbax-beaker | grep nalgebra`."* rust-dev adds the doctrine
  point: `chemistry.md:3075`'s claim that routing the models' `exp`/`powf`
  through `det_math` makes the fit cross-platform is **false** (the solver's own
  QR runs through nalgebra), *"and the sentence as written will be cited later
  as evidence that a nalgebra-backed path is cross-platform-safe, which is the
  same class of doctrine error as the vector-maths determinism claim §18.2 just
  corrected."*
- **M-c `[corroborated]`** (geometry G11, determinism) · The cavity outside-flood
  does not check `fcc::in_bounds`; `in_box` decodes via `from_flat`, which is
  `% GRID`, so a cell at `x=0` stepping to `x=−1` decodes as `x=63` and the
  flood leaks around the torus. **Same class as S5, now in Task 13.**
- **M-d `[corroborated]`** (determinism, rust-dev #5) · `fold_seed` claims to
  share construction with `key_of` *"so the two cannot drift apart"*. It does
  not — different seed constants (`0xF0_1D_5EED` vs `0xB0_1BAA_5EED`), different
  stream indices (3 vs 1), 64-bit vs 128-bit. Both are pure functions of
  `(sequence, universe)` so **the invariant holds and the stated reason does
  not.**
- **M-e `[corroborated]`** (geometry G11, rust-dev #9) · `Cavity::centre`
  documented as *"nearest the cavity's centre of mass"*, implemented as
  `cells[cells.len()/2]` — the median of sorted **flat indices**, i.e. median in
  z-major order. It sets the origin for every `r[d]`.
- **M-f `[corroborated]`** (geometry G8, rust-dev #14) · `may_bind`'s
  Cauchy–Schwarz soundness rests on `w_shape ≥ 0 && w_charge ≥ 0`, which holds
  by construction and is asserted nowhere. rust-dev goes further and would make
  it a property of the type: `BindConsts::new(..) -> Result<Self, ConstsError>`.
- **M-g `[corroborated]`** (emergence E9(b), alife alr-024) · The shadow's
  equalised rate. Both give `k* = (Σ n_i k_i)/(Σ n_i)`, count-weighted at the
  fork. Both insist the held-fixed-vs-re-equalised choice *"must be a decision,
  not an accident"* — emergence: *"recomputing it as the shadow drifts would
  feed the focal run's composition back into the control, which is no longer a
  control."*
- **M-h `[corroborated]`** (determinism, rust-dev #5) · `fold_with_seed` is
  `pub`, so *"diagnostics only"* is a comment, not a constraint. rust-dev's
  answer to the question that was asked: no marker, no sealed trait, no
  `#[cfg(test)]` — *"its only caller is a test in the same file… Make it
  private and the doc comment stops being a request and becomes a fact."*
- **M-i `[corroborated]`** (determinism, rust-dev #14 deferring) · Result-
  affecting accumulations not commented as load-bearing: the SMACOF update
  (`v0.md:3533–3550`, *"the more important of the two"* — the centroid three
  lines below it **is** commented), `may_bind`'s `.sum()` (`v0.md:4299`, gates a
  discrete accept/reject), and `contact_energy`.

### 3.6 MEDIUM — solo

`[single-reviewer]` throughout; all `cross-check-failed`.

- **emergence E5** · Solvent attack is **intensive**, never multiplied by
  `bond_count`. *"A 200-mer with mean exposure 0.5 is attacked at half the
  per-molecule rate of a fully exposed dimer, despite having ~200× the
  attackable bonds… Size therefore protects twice, and neither halving is
  derived from anything."* Sign is now correct — S4 is genuinely resolved.
- **emergence E7** · G3 (fiction) partially open. The **valence rule**
  (`v0.md:1706`) is a fixed triangle peaking mid-period with only its amplitude
  drawn — *"That is real valence behaviour, and the comment names it as such."*
  The **stability curve** `(index/90.0).powi(6)` (`v0.md:1722`) — *"`90` sits at
  the real onset of radioactivity… the most identifiable real-chemistry echo
  left. It also feeds a live physics channel."* Precedence rank **1**; the
  period-length and trend-sign fixes did land.
- **emergence E10** · `MAX_CAVITY_CELLS = 24` outlived its justification —
  *"it discards large interior voids, which are the pockets most likely to hold
  two substrates at once, i.e. the §8.5 case."* See also D-5.
- **emergence E8** · `catalytic_class` is written and never read: *"a catalysis
  path keyed on an element label, parallel to §8.5's cavity geometry."* No V0
  violation; flagged because *"the field's existence is what will make the
  shortcut look natural later."*
- **rust-dev #3** · `affinity_ordered` may swap the pair and **discards that it
  swapped**, so physics and renderer can select different rotations on ties —
  *"Same defect [as the one `affinity_with_rotation` exists to prevent], arrived
  at by a different route."* Both `affinity` and `affinity_with_rotation` are
  `pub`, so the discipline is a doc comment. Proposes one entry point returning
  `Bind { score, rotation, swapped }`.
- **rust-dev #8** · `rate()` returns 0.0 for an out-of-range species id —
  *"nothing distinguishes 'no reactant present' from 'the vectors are out of
  sync'."* And `Reaction` now carries `catalysis: f64` while `rate` still takes
  it as a parameter and ignores the field.
- **rust-dev #9** · `cavities(&Fold, &Polymer)` accepts a fold and a polymer that
  never met; `lining_signature` guards the write to `pos_of` and reads it
  unguarded. Proposes a `Folded` newtype. Also `bounding_box`'s bare tuples —
  *"`in_box(cell, hi, lo)` compiles and yields an empty box, so cavity detection
  returns nothing and no test fails."*
- **rust-dev #10** · Four `unreachable!()` in library code with a note saying to
  fix them later. *"`unreachable!()` is worse than the `unwrap()` the Global
  Constraints ban — same panic, and it reads as an assertion of impossibility
  rather than an admission."* Supplies a version needing neither. Also
  `evict_one` spins forever if every slot is `None`, and `key_of` shadows its
  `u: &Universe` parameter with a `u8` loop binding.
- **rust-perf P2** · The cavity outside-flood is **B-tree-bound, not cell-bound**:
  ~8×10⁴ B-tree descents per call, *"3–8 ms per call — likely several times the
  cost of the 20,000-step anneal that produced the fold, which inverts the
  expected balance between Tasks 11 and 13."* Replacement is **bit-identical**
  and the reviewer explains why: `box_cells` iterates z-outer/x-inner, so box
  order **is** ascending flat-index order — exactly `BTreeSet` iteration order.
- **rust-perf P3** · `SpeciesRecord` is **832 bytes** at D=42 (2,752 at D=162);
  the rescan reads 32 of them from up to 256 records — *"208 KB, over M1's 128 KB
  L1D, so every rescan is an L2-latency walk."* Proposes hot/warm/cold split.
  *"Do it now rather than later: the retrofit touches every call site, and there
  are currently none."*
- **rust-perf P5** · `ActivityStats::per_species` reads as a per-step sweep:
  *"at S=10⁴ and 10⁵ steps/sec that is 10⁹ updates/sec — three orders past the
  machine, and it would make Task 20b, a metrics task, the dominant cost."* Also
  wrong on its own terms — *"Gillespie steps have variable Δt, so a per-step
  increment of 1.0 measures event count, not persistence."* **This one has a
  natural corroborator it never reached: rust-perf explicitly says *"I would
  want the alife-researcher to confirm it matches the Bedau definition"*, and
  alife's transcribed file does not address it.** See §5.
- **rust-perf P6** · `search()` allocates a `Vec` inside a loop over 256 colour
  values, in a search bounded at 50,000 leaves — *"purely to learn whether a
  colour class has more than one member. The answer is a popcount."* Loop bound
  also drops from 256 to `n` after the densification fix.
- **rust-perf P7** · `AffinityMemo` uses `Vec<Option<f64>>`; `f64` has no niche,
  so the table is **2× the needed size** — 134 MB vs 67 MB at k=4096. Affinity
  is a negated sum of squares, so `NAN` is a sound sentinel.
- **geometry G9** · `the_mirror_complement_does_not_fit` is *"the right pair of
  tests"* and each fails if `anti` is dropped — **verified numerically, roles
  swap exactly**. But *"the seed is weak and the threshold is far too loose"*:
  discrimination is 0.06% of the natural scale against a `−1e-6` assertion,
  three more orders below that. A branched seed gives ~50× the margin.
- **geometry G11** · `SHELL = 0.75` is absolute while extents scale with
  generated radii (0.56–4.55) — *"in a small-radius universe the shell spans the
  whole molecule and `a[d]` collapses to a global average."* Same class as G7.
- **determinism D10** · `depth_sort` declared, never specified. Correctly
  **downgraded** from the previous round: *"a float-keyed `sort_unstable_by` is
  still reproducible here… this is **not** a cross-platform hazard. It is a
  fragility hazard"* that churns `insta` goldens and trains the team to accept
  snapshot diffs unread.
- **determinism D12** · `Domain::Shadow` declared, never used; Task 20b Step 3
  states no RNG discipline, *"and the obvious implementation clones the focal
  run's stream, which is precisely the duplicate-replay that dropping `Copy`
  exists to prevent — the shadow would then be a correlated control."*
- **determinism D11(b)(c)** · The `ln` in AIC/BIC is unrouted (*"AICc and BIC
  routinely differ by <1 unit, which is exactly the margin where a last-bit
  difference flips the verdict"*), and Levenberg–Marquardt terminates on a
  tolerance — *"the hazard SMACOF's fixed 240 iterations exists to avoid."*
  D11(a) explicitly **defends** the headline reproducibility claim.
- **alife alr-013** · *"The null-trajectory test will fail as written."*
  Plug-in histogram-distance estimators are positively biased; subsampling to
  constant N removes the N-dependence but not the K-dependence. Change the
  criterion from *"reads ~0"* to *"constant in time, matching the analytic
  null"*.
- **alife alr-014** · Novelty increments are **counts**, so variance scales with
  mean and homoscedastic-Gaussian RSS is violated *"exactly where models differ
  most."* Use Poisson/quasi-Poisson likelihood or variance-stabilise first.
  *"Otherwise 'fit increments' swaps one wrong error assumption for another."*
- **alife alr-021** · The bounded model in the source (Wiser, Ribeck & Lenski
  2013) is **hyperbolic**, not exponential-saturating. *"Hyperbolic decays as
  t^-2 and is much harder to separate from a power law — it is the
  discrimination that actually tests the boundedness illusion. Costs one more
  model."*
- **alife alr-025** · The randomised-catalysis control's *"holding density
  fixed"* is under-specified: uniform reassignment also destroys the degree
  distribution, *"so it measures 'shape chemistry vs uniform' not 'vs
  same-density random'."* Needs a degree-preserving edge shuffle.
- **alife alr-027** · Neutral-network percolation threshold is
  `λ* = 1 − κ^(−1/(κ−1))`, depending on **alphabet size κ**, which Borbax
  generates per universe. *"A fixed LCC-fraction bar would reject universes for
  having a small alphabet."*
- **alife alr-020** · AICc's `K` must include σ², so K=3 for 2-param models.
  *"Cancels in AIC/BIC; does NOT cancel in AICc. At n=10 the ΔAICc error is
  ~1.07 and biases toward the 2-param models."*

### 3.7 LOW / mechanical

alife alr-030 (`find_raf_present` referenced, defined nowhere), alr-033
(plastogenetic congruence is partly an annealer artefact; a fixed-temperature
chain after burn-in is a faithful variant at no extra cost), alr-035 (power-law
increment singular at t=0 for b<1 — index from t=1), alr-036 (AICc/BIC on the
training half only); determinism's `CacheStats.bytes` using
`size_of::<Fold>()` (fine, but must never reach a golden), `xtask`'s
filesystem-order `read_dir`; rust-dev #11 (`GeoError` hand-rolled with no
`Display`, and `UnsupportedResolution` is a compile-time fact reported at
runtime — a sealed `Resolution` trait removes the error path entirely), #13
(`(ShellPattern, Vec<Element>)` wants a name), #14's cast/panic sweep
(`self.occ[cell as usize]` unbounded, `exposure` calling `is_free` with no
`in_bounds` guard, `f64 as usize` saturating silently, `SpeciesId(pub u32)`
letting anyone panic the interner); rust-perf P8; emergence's G4/G6 fiction
notes (*"falls to zero at the noble end"* imports real vocabulary; the bare-`f64`
leaks are *"how G4 is defeated quietly"*).

### 3.8 Retracted — dropped, recorded as required

- **geometry G8, tighter `may_bind` bound.** *"I also tested a strictly tighter
  bound (rearrangement inequality on the sorted extents, which is exact over all
  bijections rather than just Cauchy–Schwarz): 4.0% vs 3.3% rejection on random
  pairs. Not worth the complexity — **withdrawing that suggestion.**"*
  A reviewer retracting its own suggestion on measurement is the process working.

---

## 4. Two misattributions, one introduced by the fix pass

- **alr-010 `[single-reviewer]`** · §2.7's 3b/3c subdivision is **Channon 2003
  (ALIFE VIII)**, not Bedau/Snyder/Packard 1998. *"First three clauses verified
  correct; the fourth is Channon's extension."* **INTRODUCED BY THE FIX PASS** —
  fix pass 1 corrected M7's class label and imported a new citation error doing
  it.
- **alr-022 `[single-reviewer]`** · The MODES persistence filter is **lineage**
  persistence on a **phylogeny**. *"Borbax has no phylogeny (species arise by
  reaction, not descent) so it is not implementable as published. Also the plan
  calls it 'the substitute for a full shadow' while building the shadow —
  inverts the relationship."*
- **alr-023 `[single-reviewer]`** · *"'Lehman-Stanley archive measure' is wrong:
  theirs is k-NN mean with k~15, chosen because a pure minimum is too noisy.
  Min-distance is ASAL's choice."* **The plan's `min`-not-mean choice is
  correct** — it is faithful to ASAL, verified against `SakanaAI/asal`. Only the
  attribution is wrong. Note that **emergence independently praised this same
  choice** (E9: *"The `min`-not-mean novelty choice… [is a] guard against
  measuring an answer into existence"*) — so the design is corroborated and only
  the citation needs changing.
- **alr-034 `[could-not-verify]`** · The AlChemy 16%/60% figures could not be
  confirmed; one render gave ~13%/17%/70% and the paper emphasises *dominance*.
  These are in the **spec**, not the plan. *"Qualitative conclusion (never assume
  composition) is unaffected."*

---

## 5. Coverage gaps

### G-1 · `alife-researcher` produced no file — findings were transcribed by the parent. **[skill bug]**

Its own first line: *"(Transcribed by the parent: this agent has no write tools.
Fix its `tools:` line.)"* This is a bug in the calling skill, not a reviewer
failure. **Verdict: `missed`, and the loss is measurable.**

| | alife | median of the other five |
|---|---|---|
| File length | 92 lines | ~500 |
| `file:line` references | **0** | dozens each |
| Proposed fix code | **0 blocks** | 5–15 blocks each |
| Proposed tests | **0** | several each |

Specifically lossy, in descending order of what it will cost to recover:

1. **No line references anywhere.** Every other reviewer anchors findings to
   `file:line`. Acting on alr-011/012/013/014 means re-locating the code first.
2. **alr-013 and alr-014 are numerical-analysis findings with no corrected
   formula.** alr-014 says *"use Poisson/quasi-Poisson likelihood, or
   variance-stabilise (Anscombe/sqrt) first"* — that is a direction, not a fix,
   and it is the kind of thing this reviewer would normally write out.
3. **rust-performance-expert P5 explicitly asked for this reviewer and did not
   get it**: *"the exact-integration form is standard and I would want the
   alife-researcher to confirm it matches the Bedau definition they intend."*
   That confirmation is missing, and P5 is a rewrite of the activity metric.
4. **alr-034 is a spec finding with no pointer into the spec**, and it asks for
   a decision (fix the figures against the published table, or drop them).
5. The four `VERIFIED CORRECT` entries (alr-000..003) are the round's cleanest
   result — RAF closure verified against Hordijk & Steel 2004 by walking the
   algorithm — and they are compressed to four lines. As the regression-
   protection record for §6, that is thin.

**Recommendation: re-run this reviewer with a `Write` tool before acting on
Task 20b or Task 17.** It is the cheapest single item in this whole review.

### G-2 · Nobody reviewed the spec. **[dispatch gap, not a reviewer gap]**

All six reviewed the two plan files. The PRD
(`docs/superpowers/specs/2026-07-26-borbax-prd.md`) was not in scope for anyone.
alr-010 (Channon) and alr-034 (AlChemy) are **spec** defects surfaced
incidentally, and §2.7 is where fix pass 1 introduced a new citation error. Per
my Step-4 rule: **a domain in scope with no reviewer dispatched is a bug in the
calling skill.** Saying so explicitly.

### G-3 · Task 19 (rendering) has no owner. **[dispatch gap]**

Coverage is one determinism finding (D10, `depth_sort`) and one type finding
(rust-dev #3, renderer/physics rotation divergence). Nobody reviewed the
projection, the `insta` golden mechanism, or whether the diagram agrees with the
physics. `CLAUDE.md` states the reason SVG exists from the start — *"a
shape-based chemistry cannot be developed blind"* — which makes the renderer a
correctness instrument, not decoration. If it is wrong, every visual judgement
made during chemistry tuning is wrong.

### G-4 · Tasks 1–3 (workspace, CI, RNG) reviewed by determinism only.

D2/D3/D5/D9 are thorough on the determinism axis. Nobody else looked. Low risk —
the tasks are mechanical — but it is silence, not a clean bill.

### Intentional clears — genuine results, recorded as such

Not gaps. Each reviewer produced a substantial "checked and found sound"
section; see §6.

---

## 6. What the plan does well — do not let a third pass regress this

Assembled from all six "checked and found sound" sections. **Anything on this
list that changes should change deliberately.**

**Verified exhaustively, by computation:**

- **The rotation table, at all three resolutions** (geometry). 60 distinct
  permutations, closed under composition (3,600 products each), inverse-closed,
  identity present exactly once, **every element det = +1.000000000000** with
  `RᵀR − I` within 1.6e−15. *"Zero improper elements. 60, not 120 — §22.8
  holds."* `anti` is an involution, and `anti[perm[i]] == perm[anti[i]]` **as
  integers** at all three resolutions.
- **The `anti` fix genuinely searches proper rotations**, and
  `affinity_ordered`'s symmetry argument survives it — *"Verified numerically:
  the difference is exactly 0.0."* S1 is properly closed.
- **`may_bind`'s bound is real** — 12,000 random pairs, 0 violations, holds for
  every rotation rather than on average (see D-6 for the caveat).
- **The RAF closure fix is correct** (alife alr-000) — walked by hand,
  *"Matches Hordijk & Steel 2004"* — and all four doc-comment claims check out.
  *"RAF detection is detection"* (emergence): nothing anywhere boosts a reaction
  because it is in a RAF.
- **`libm`'s `arch` feature is safe for the five functions used** (determinism,
  by downloading and reading 0.2.16): arch overrides only `sqrt`, `fma`, `rint`,
  `ceil`, `floor` — all IEEE-754-exact. `exp/log/pow/sin/cos` contain no `fma`
  and no arch dispatch. **Only `cbrt` reaches the CPU-feature detection**, which
  is why D5 wants it in the banned list.
- **Task 6's densified colouring genuinely fixes S6** (geometry) — the rank is
  dense in `0..n_distinct`, isomorphism-invariant, and array bounds hold at
  every recursion depth.
- **FCC parity does not break the outside flood** (geometry) — the flood
  propagates only through `NEIGHBOURS`, which preserves parity, so the two
  sublattices flood independently and odd-parity cells can never produce
  spurious cavities.

**Structures and disciplines all reviewers want kept:**

- **The segment tree** — *"the best-defended determinism decision in the plan"*
  (determinism), *"Best-defended structure in the plan"* (rust-perf). Parents
  recomputed as `left + right`, bit-identical to a rebuild; freed slots set to
  0.0 rather than compacted, with a test.
- **Damage stays implicit** — *"the best-defended invariant in the plan"*
  (emergence). No damage field, no age, no per-molecule struct anywhere.
  Task 21 Step 3 now tells an implementer to rework the measurement before
  reaching for the provenance fallback, with the asymptotic reason attached.
- **Exactly one `HashMap`, verified never iterated** (determinism). `BTreeMap`/
  `BTreeSet`/`VecDeque` throughout; `PassThrough` replaces `RandomState`. **No
  rayon, no threads, no `Instant`, no `env::var`, no `rand::`, no `f32`** —
  grepped exhaustively.
- **`Stream` is not `Copy`**, and the M9 duplicate-stream instance is gone;
  `Universe::generate` now carries the reason inline.
- **Fixed 240 SMACOF iterations and 20,000 anneal steps** rather than
  convergence tests. Aligned RMSD/R_g is 0.12 — *"SMACOF is doing its job"*
  (geometry).
- **The fold workspace trail reset is sound** (geometry, traced) and
  **incremental contacts equal a full recount exactly** — integer histogram, so
  the assertion is exact rather than tolerant. *"Correct call."*
- **`propose` cannot break self-avoidance or connectivity** (geometry) —
  *"Sound; the problem is the start (G3), not the move set."*
- **The Metropolis fix is right** (rust-dev, emergence, determinism), and
  drawing `roll` **unconditionally** so stream position never depends on a float
  comparison is *"exactly the discipline §13.1 needs."*
- **Catalysis really does fall out of two cavities and nothing else**
  (emergence) — no `is_enzyme`, no enzyme type, no species-name branch, no rate
  bonus keyed on a classification. *"Subject to E1, the shape of the mechanism
  is right."*
- **Exit criteria 7 and 9 are now separable** — decay-*off* for criterion 7,
  decay-*equalised* for §15.3. *"Conflating them was the S12 defect and it is
  properly closed."*
- **Task 20b is measurement throughout** (emergence) — *"deleting the whole task
  leaves the chemistry identical."* The `min`-not-mean novelty choice, the
  constant-N subsampling, the increments-not-cumulative fit and the held-out
  tail are *"all guards against measuring an answer into existence."*
- **The randomised-catalysis control is the right second control** and asks the
  question the persistence shadow cannot.
- **`AffinityMemo` is proven equal to direct computation by test** rather than
  assumed; the fold cache memoises a pure function with result-neutral eviction.
- **The dev/release opt-level split is not a determinism concern** (determinism):
  LLVM does not reassociate f64 or contract to FMA without fast-math.
- **Test-pair design.** `a_physical_complement_scores_zero` +
  `the_mirror_complement_does_not_fit` together pin handedness and *"either
  alone would not"* — and geometry verified they are mutually self-checking: a
  degenerate seed makes them contradict each other, so it fails loudly.
  rust-dev: *"The comment explaining that omitting `anti` converts the search
  into the improper elements is the single most valuable comment in the plan."*
- **The dependency doctrine correction from round 1 held.** §18.2 now records
  the vector-maths hand-roll as a size judgement, not a determinism one. rust-dev
  flags one new instance of the same error class (M-b) — evidence the lesson is
  being applied.

---

## 7. The question you actually asked: third fix pass, or something structural?

**Something structural. And the reason is sharper than "fix passes introduce
defects".**

### What the evidence says

Fix pass 1 applied ~30 findings and introduced at least four new defects:
alr-010 (citation error, explicitly *"INTRODUCED BY THE FIX PASS"*), E2 (contact
energy sign-inverted — verified by me), G3 (boustrophedon retrace — verified by
me), and G2 (`canonicalise_frame` made the metric measurably **worse**: 0.87
median against 0.64 for no frame at all).

But the failure has a **shape**, and four reviewers named it independently
without seeing each other:

> determinism D1: *"the code was edited and the prose below it was not."*
> determinism, smaller items: *"the fix pass clearly edited the comment and not
> the code — **the same pattern as D1**."*
> emergence E1: *"the catalysis fix is **a comment, not a change**."*
> rust-dev #2: *"A false 'already fixed' comment is worse than no comment."*

Four instances of one mode: **the fix landed in the prose describing the code
and not in the code.** Task 16, Task 10's `det_math`, Task 21 Step 6, Task 15's
`radiogenic_rate`. In three of those the header now actively asserts the fix is
done.

That is worse than the defects themselves, because it **destroys the review
artifact's own reliability**. A reader cannot tell a fixed section from an
unfixed one by reading it.

### The decisive observation

I counted **nine defects in this round that `cargo check` or `cargo clippy`
would find in seconds** (§3.2): `E0560`, `E0063`, `E0599`, `E0412`, `E0425` ×3,
`unused_variables`, `erasing_op`. They survived **two review rounds by six
specialists** — because six specialists reading prose is a strictly worse
compiler than a compiler, at perhaps a thousand times the cost.

Meanwhile the things review is uniquely good at — G1's measured 7% rotation
contrast, G2's 32% locality failure, E2's sign table, G10's triangle-inequality
violation — are *exactly* the things a compiler cannot see, and they were found
by two reviewers who **stopped reading and started computing**.

The process is currently spending its most expensive resource on the cheapest
class of defect.

### Second decisive observation: G2 is not a fix-pass item at all

Geometry measured that `canonicalise_frame` is worse than its own absence, that
**32% of one-atom mutants land further from their parent than an unrelated
molecule does**, and — critically — that the root cause is that *canonical atom
order is not continuous under graph edits*, so **no frame anchored to canonical
order can fix it.** Five alternative frames were tried and all measured the same
or worse.

If that holds, §8.4's evolvability argument has nothing under it, neutral
networks cannot form, and no amount of fixing Tasks 10–21 reaches the exit
criteria. **That is a possible falsification of a design premise, and it does
not belong in a defect list.** It belongs in an experiment, and the reviewer was
right to say *"I do not have a validated [fix], and I would rather say that than
invent one."* A third fix pass would either skip it or invent one.

### So: split the work by what can verify it

**Track A — settle the physics experimentally, before Task 8 lands.**
G1, G2, G4, G5, G7, G10 and E2 all compile fine and are all *measurable*.
Geometry already built the harness (a Python reproduction of `smacof`,
`canonicalise_frame`, `signature`, `canonicalise`, `affinity`, `may_bind`, the
geodesic construction and the FCC walk) and it produced every number in that
review. **Promote that harness from a scratchpad artifact to a checked-in
experiment.** Each of these findings comes with the reviewer's own acceptance
test; run them. G2 is the gate — if signature space has no locality, the answer
changes the plan, not a task.

**Track B — stop editing prose about the compile-caught defects. Execute
Task 1.** Every item in §3.2 is discovered and fixed in minutes once there is a
`cargo check`. Editing them by hand in a markdown file is how B-1 through B-9
got here.

**Track C — one mechanical audit, not a fix pass.** Every claim of the form
"constant X was replaced by Y" gets a grep asserting X is absent. I ran four
such greps and they found B-1, D-a, D-d and M-a. That audit is minutes, is
repeatable, and would have caught the entire "fix landed in the comment" class
before this review round was commissioned.

### On the Step-5 escalation threshold

Six genuine conflicts is at the threshold, and I am **not** invoking the
standard escalation ("the artifact is internally inconsistent; resolve before
acting"). The disputes are not scattered — **four of the six (D-1, D-2, D-4,
D-5) are the same argument**, between the reviewer who measured the geometry and
the reviewers who read it. That is not incoherence; it is a signal about
*method*, and it points the same way as everything above: **on this artifact,
measurement beats reading, and the review process should be reorganised around
that.**

---

## 8. Recommended next action

1. **Re-run `alife-researcher` with a `Write` tool.** Cheapest item here. Its
   findings are currently a lossy transcription with no line references, and
   `rust-performance-expert` P5 is explicitly blocked on a confirmation it never
   gave. (§5, G-1)
2. **Run the mechanical claim-audit (Track C) — today, ~15 minutes.** For every
   "X was replaced by Y" in either plan file, grep that X is absent. Fix the
   prose/code divergence class in one pass, and **add this grep to `xtask`** so
   it cannot recur. This is the process fix, and it is the highest
   value-per-minute item in the whole review.
3. **Decide G2 before anything else touches Tasks 8–10.** Promote geometry's
   reproduction harness to a checked-in experiment and run its proposed test
   (`a_one_atom_mutation_stays_near_its_parent_in_signature_space`, over random
   trees with `canonicalise()` called on the *graph*). Then evaluate the
   reviewer's own proposal: an **SO(3)-invariant descriptor** (`sorted_extents`,
   or a low-order spherical-harmonic power spectrum) for species identity and
   novelty, leaving `affinity`'s rotation search to handle orientation —
   *"it already does"*. **If signature space has no locality, this is a design
   decision, not a fix.** Do not apply rust-dev #4's collinear repair until the
   function's survival is decided (§2, D-4 — precedence rank 2 over rank 6).
4. **Settle the measured physics set in the same harness, in this order:** G1
   (`pair_gap` — and get emergence-auditor's verdict on whether a per-pair gap
   launders the answer in), E2 (apply the inversion fix; ask geometry about the
   neutral-neutral degeneracy I found), G10 (Floyd-Warshall target + bond-order
   shortening — the triangle-inequality half is settled, the shape-effect size
   is not), then G4/G5/G7 together, which is dispute D-1 and needs one
   synthetic-ellipsoid test to resolve.
5. **Execute Task 1 (Track B)** with three changes folded in before the first
   commit: `libm = "=0.2.16"`, `Cargo.lock` added to the `git add` list, and
   `shell: bash` on the goldens emit step. These are D-b and D-c, both verified,
   both one line, both far more expensive after goldens exist.
6. **Delete the shadow `det_math` instruction and make the xtask exemption
   path-based** (`v0.md:4323–4338`, `4367–4370`; repoint `chemistry.md:1898,
   2046, 2067, 3095`; rewrite Task 21 Step 6). Then apply the D-3 ruling: drop
   `.powi(` from `BANNED_CALLS`, write out `v0.md:1722`'s `.powi(6)`, and
   correct `det_math`'s doc comment so it stops arguing for a rule the plan does
   not follow.
7. **Expose `V` and derive the activity threshold from the shadow** before
   Tasks 18 and 20b are written (S-a, S-b). Both are interface changes that are
   cheap now and expensive later, and both are the difference between a
   diagnostic that works and one that cannot distinguish two hypotheses.
8. **Apply E3's fix, not P4's** (S-c). Computing catalysis at catalyst-intern
   time closes the §8.6 breach *and* subsumes P4's O(C²)→O(C) hoist. Applying
   P4 alone leaves the breach open while looking like it was addressed.
9. **Before Task 10 lands, decide rust-perf benchmark 1's second variant.**
   Fixed-shape pairwise summation in the affinity kernel is deterministic (the
   tree is fixed at compile time) and *more accurate* than the serial sum, but
   it moves every golden. rust-performance-expert correctly flagged it as
   determinism-auditor's call and cited **`CLAUDE.md`'s precedence clause "2
   above 3"** — adopt a numerically better formulation and regenerate goldens
   rather than preserving a worse one to protect a hash. **There are no goldens
   yet, so this decision is free exactly once, and this is the window.**
10. **When the fix pass is done, re-run only `geometry-numerics-reviewer` and
    `emergence-auditor`** — the two whose findings this round were mostly
    *measured* rather than read, and the two `CLAUDE.md` names as the ones
    *"no test suite fully protects"*. Give them the harness from step 3. Do not
    re-run all six on prose.
