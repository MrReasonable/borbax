# Borbax Real Physics — plan (issue #26)

**Status:** specced 2026-08-07, revised twice same day after two `/review-plan`
rounds (12 reviewer-passes total). Not started. An epic under issue #26 —
prerequisite guard/seam fixes, then three sequenced tasks sharing one
`Element` struct and coordinated golden regeneration.

## Why this exists

Issue #26 (Ian, 2026-08-04): two defects that look separate and are not.
**Slot 0 has valence 0 in every universe measured** — a one-unit element can
never bond. **Every monomer has negative binding energy and instability
pegged at 1.0** (#25). Both come from applying formulas derived for large
clusters at N=1.

**What changed since the issue was filed: the project's own purpose changed,
twice, and the plan changed with it each time.** First, reviewing this
epic's draft surfaced a tension in §5's G1 (no real chemistry data) — a
calibration step that had to hit a real target was real chemistry data
however it was phrased. Ian's ruling went further than the plan asked:
*"G1 doesn't need to protect. We're not 'inventing chemistry' anymore. We're
demonstrating how it really works in an interactive tool,"* and, tying it to
precedent already in §5, *"That is why we removed G3 through G6."* G1 and G2
were revised the same day — see spec §5 for the exact wording.

Second, a full review round found that "the identity configuration reproduces
the real periodic table" was underspecified in two ways: whether it includes
real *names* (not derivable from physics, need their own guard), and how much
of the physics it must reach (only the element table, or every perturbable
constant in the simulation). Ian's rulings, both more ambitious than the
plan's draft assumed: **real names for the identity configuration**, and
**everything physics-related, table included** — every perturbable constant
in the simulation, not just `Element`'s own properties, should bottom out at
its real value when the seed is unperturbed. Both are now Decisions 9 and 10,
below, and both carry real architectural consequences worked out in this
revision.

## Scope of this epic

**In scope:**
1. **Prerequisite guard/seam fixes** (below) — land first, own commits. Now
   includes a cross-cutting constant-drawing convention (Decision 10) that
   every perturbable constant reachable from the identity configuration must
   use, and a real-name guard (Decision 9).
2. **Task 26.1 — electron configuration.** A perturbable, self-consistent,
   `l`-dependent orbital-filling generator (screened-hydrogenic), landing as
   `PhysicsVersion::V2`, additive. Produces real periods, groups, valence,
   affinity, radius, and (new, from the second review round) period
   boundaries and V2 bond capacity — all derived from one quantity the fill
   already computes (Decision 8), rather than four separately-specified
   outcomes. Folds in validated fixes to §8.3's binding kernel, the
   shape-signature character channel, and a rank-2 bond-energy formula
   (closes #31) — plus the prefilter-tightness fix #31's own closure turns out
   to require.
3. **Task 26.2 — nuclear structure and isotopes.** Lands as
   `PhysicsVersion::V3` (not V2 — see Decision 6). A second, independent axis
   giving each element a composition, with a *local* stability definition
   (Decision 7) rather than V1's global one — real nuclear physics is local
   (hydrogen-1 is the least-bound nuclide and is perfectly stable), and this
   is what actually closes #25, not the composition axis by itself. Isotopes
   are an abundance-weighted **distribution** on one `Element`, never a
   separate species identity (Decision 5) — the naive "isotope = new
   `ElementId`" design was found to panic under realistic table sizes and to
   silently corrupt canonicalisation.

**Explicitly deferred, with reasons:**
- **Flipping `PhysicsVersion::CURRENT`.** Its own commit, after §7.2's battery
  exists (it doesn't yet). Needs a review checkpoint attached when it happens
  — the binding kernel's shape/charge balance and the prefilter's tightness
  bar have never been exercised against real V2/V3 data (see Task 26.1's
  verification steps), only against a one-off cross-check measurement.
- **The deeper tension in `signature.rs` between two conflicting stated
  intents.** F2 (below) has a working, now fully-specified remedy; the file's
  own older self-contradiction is a separate `borbax-molecule` design
  question.
- **`embed`'s approximate (not exact) scale-equivariance on symmetric rings**
  at n=9,12 — found incidentally during this review, deterministic, not a
  §13.4 hazard, unrelated to this epic. Raise against `layout.rs` separately.

## Decisions settled during design (2026-08-07, two review rounds)

1. **Shape-complementarity binding (§8.3) needs no code changes** to its
   core kernel. Traced every production read of an `Element` field from
   `borbax-molecule`: `radius`/`affinity` feed the shape signature, `valence`
   gates bond legality, `mass` sums molecular mass. (The prefilter *does*
   need a change — see Decision 4.)

2. **G1 and G2 are revised, seed-scoped.** `check_no_data_files` stays
   absolute — no literal data file, no lookup table, under any seed. What's
   permitted: a generative formula's parameters may equal, or be chosen to
   approximate, their real physical counterparts, because there is one
   calculation for every seed and a seed landing on real-reproducing values is
   a point in that calculation's parameter space, not a special case. Full
   wording: spec §5.

3. **Energy form: screened-hydrogenic (`Z_eff`-based), and it needs
   `l`-dependent screening specifically, not merely self-consistency.**
   Self-consistency alone (screening that doesn't depend on the candidate
   subshell's own `l`) was measured to produce **exactly one** fill order over
   the entire coefficient space — no d-block, no periodic structure at all.
   The mechanism that actually produces ordering anomalies is **orbital
   penetration**: the candidate's own `l` must enter the screening through a
   function that is ≈0 at `l=0` and essentially saturated by `l=2` — a step
   is the simplest such function, not the only one. Measured against
   Slater's own step function: `t(l) = (0, 0.5, 1, 1)` retains a real
   d-block at 98–100% of draws out to ±50% perturbation, against the Slater
   step's 8% — a smoother, more robust, *more* physically faithful choice.
   Pin the energy denominator as `−(Z_eff/n)²` (plain `n`, not Slater's `n*`)
   — measured, plain `n` reproduces the real fill order to 5p at the identity
   point and `n*` does not, despite `n*` being individually more standard.

4. **Closing #31 (affinity's degenerate range) breaks a different proven
   bound, and both need fixing together.** The binding kernel's prefilter
   (`ceiling`) is currently tight partly *because* affinity never straddles
   zero — its dominant term is a size proxy the prefilter can reject
   confidently. Once affinity is fixed (Task 26.1), the existing prefilter's
   tightness bar (`worst_gap < 0.98`) fails outright (measured: rises to
   ≈1.0), and it goes selectively blind on self-pairs — exactly the class
   §10.1's membrane self-assembly needs to discriminate. **Fix: a
   rearrangement-inequality prefilter** (sort each molecule's shape/charge
   channels, dot the sorted sequences) — strictly tighter than the current
   bound everywhere, brings `worst_gap` back under the bar, restores
   self-pair discrimination, and costs one `O(D log D)` sort per species at
   intern time (§8.6-compliant, ~3% of the kernel's own cost). This was
   already routed to a prior task and never executed; it executes now,
   bundled with Task 26.1 since #31's fix is what activates it.

5. **Isotopes must not get their own `ElementId`.** Four independent review
   angles (type-system ceiling, canonicalisation-memo cost, intern-space
   growth, graph-automorphism breakage) converged on the same conclusion, and
   two reviewers' numbers close it decisively: the flat-`ElementId` scheme
   permits at most **2.13 isotopes/element** at a typical table size. Against
   that ceiling: because the composition cost is strictly convex, exactly
   one `a` (proton count) is composition-stable at each total under Decision
   11's corrected model, so `mean stable isotopes per element == t_max /
   n_elements` **exactly** (verified at n = 30, 60, 120) — not an
   independent "a working stability model demands" figure, a direct
   consequence of the model's own algebra once `t_max`/`n_elements` is
   chosen (Decision 11 pins that ratio to a 2.52–2.95 band). 2.13 sits below
   that floor for every legal choice, so the conclusion is unchanged; only
   its derivation is restated here, because the figure's earlier phrasing
   made it sound like an independent empirical target rather than what it
   actually is — a function of `c = gamma/(2*kappa)` and `n_elements` that
   Decision 11 now derives explicitly. The scheme would panic in library
   code on realistic seeds regardless. The subtler reason
   it's mandatory rather than merely expensive: isotope substitution with a
   separate `ElementId` silently changes canonical atom ordering, which
   changes SMACOF's embedding pivots, which reorients (not reshapes) the
   molecule in a way the 60-rotation binding search cannot absorb — isotopic
   substitution would become an undetected *chemical* change. **Fix**:
   `ElementId` stays one-per-chemical-element; isotopic composition is an
   abundance-weighted distribution attached to `Element` (mass and
   `instability` as share-weighted averages over a fixed ascending order — no
   sort, no accumulation-order question). Required test: two molecules
   differing only in isotopic composition give `group_distance == 0` and
   bit-identical `affinity` against any third species. This decision is made
   *before* Task 26.2 Step 1, not informed by it — the numbers that would
   "inform" it were actually already conclusive.

6. **Task 26.2 lands as `PhysicsVersion::V3`, not `V2`.** Landing both new
   tasks under `V2` would mean a shared universe address (`@2`, made visible
   on screen by Task 26.1's own wiring) names two different physics
   depending on which task has landed — violating `PhysicsVersion`'s own
   documented contract that a universe keeps the physics it was born with.

7. **`instability` must be redefined as a local (Q-value) quantity, and this
   — not the composition axis by itself — is what closes #25.** V1's
   `instability` is "distance below the binding-energy peak," a *global*
   comparison; under that definition exactly one element is ever "stable," by
   construction of the argmax, in every universe — before or after Task
   26.2's composition axis lands. Real nuclear stability is *local*: a
   nuclide is stable if no accessible transition has positive energy release,
   which is why hydrogen-1 (the least-bound nuclide by binding-energy
   standards) is nonetheless perfectly stable — nothing lower is reachable
   from it. **Fix**: `instability`/decay derives from a Q-value — the best
   available transition's energy release, clamped at zero (arithmetic, not a
   declaration: zero is where a transition stops being exoergic, fixed by the
   energy function's own structure) — evaluated over both a
   composition-changing channel and the existing size-shedding channel.
   Measured under this definition: 33 elements have ≥1 stable isotope (not
   1), stable isotopes/element averages 2.5–3.0 (which is also Decision 5's
   demand figure) — **round 4 correction: "with no tuning" no longer holds
   as stated, now that Decision 11 pins this same ratio deliberately via a
   chosen `t_max`/`n_elements` target rather than reading it off an
   untuned draw; see Decision 11's round-4 correction for the derivation.
   The absence of *ad hoc* tuning (no branch naming an outcome, no
   per-element special case) still holds and is the property that
   matters** — and the heaviest stable composition is a
   genuine emergent property varying across universes — not a constant. This
   also removes `element.rs`'s existing `if deficit > 1.0 { 1.0 }` cap, whose
   own comment concedes it's a category error, and its scale constant is the
   **same missing constant** Task 15's radiogenic decay channel already
   needs — one golden regeneration, shared with that known-open item, not
   two.

   **Dependency note, added at pre-registration (2026-08-12) rather than
   left to be rediscovered at Step 5: this constant must land wired, not
   stubbed, in the same commit that removes the cap.** `element.rs:641-642`
   (`find_peak_and_set_instability`) is where the cap currently lives and
   where its replacement will read the scale constant; its own comment
   already points at Task 15 rather than naming a value. If the cap is
   simply deleted (or replaced with an unbounded pass-through) without the
   scale constant actually being wired to something, `instability` becomes
   a well-defined Q-value that nothing downstream converts into a decay
   propensity — the composition-changing channel this task exists to add
   would be computed and never read, silently degrading the RESOLVED
   decay-only coupling above from "decay-only" to "decay-none." Track this
   explicitly through Step 5/6 rather than assuming Task 15 will supply it
   later; Task 15 itself is not in this task's scope, but the constant is.

8. **Four previously-unspecified V2 outcomes — valence promotion, `affinity`,
   period boundaries, V2's bond capacity — are all answered by one quantity
   the fill already produces**, found independently by three reviewers and
   assembled by the coordinator: **the energy gap between the highest
   occupied and lowest unoccupied subshell.** Period boundaries sit where
   that gap is large relative to within-period gaps (energetic, not a
   labelling convention). Valence promotion is "unpaired states reachable
   within an energy budget, where the budget is the gap and one drawn
   constant" — this gets beryllium (divalent) *and* helium (inert) right,
   where either "closed subshell ⇒ valence 0" or "closed subshell ⇒ promote"
   alone gets one of them wrong, and it survives the emergence-auditor's
   discriminator because it names no outcome, only an energy comparison.
   `affinity` is a centred function of the same gap, straddling zero by
   construction — closing #31 *structurally* rather than by an offset. V2's
   bond capacity (the thing Task 26.1's rank-2 bond formula needs and V1's
   `ShellPattern`/`contact_density` can't supply under an orbital model) is
   the same gap again. This is §8.3's one-mechanism discipline applied one
   layer down, and it's what makes Task 26.1's wiring steps startable at all.

9. **Real names for the identity configuration — with an enforceable guard,
   not a workaround.** Ian's ruling: use real names, not just real structure.
   Names aren't derivable from physics (there's no calculation whose
   parameter space contains "Carbon"), so this needs an actual string table —
   and the review found `naming.rs` is the one file the literal-scanning
   guard (`check_no_real_chemistry_in_literals`) exempts, verified by
   planting a real-name table there and confirming `cargo xtask` passes
   clean.

   **"The existing element blocklist already pairs 1:1 with what the
   real-name table needs to contain" was false, and it removed the argument
   that building on it is cheap.** Counted from source: `REAL_ELEMENT_SYMBOLS`
   has 118 entries; `REAL_WORDS`, before the separate fix below, covered only
   75 of the 118 real element names — 43 absent, the first gap at
   praseodymium (Z=59) — and neither array *pairs* a symbol with its name in
   the first place, being flat exclusion lists rather than a lookup. "Build
   the real-name table as its companion" assumed a pairing that does not
   exist, and the 43 missing names had to be written by hand regardless of
   this decision (they were, on `main`, independently — see the standing
   G2 blocklist fix below). Two further consequences the false framing hid:
   if the two arrays were merged so they genuinely paired, the merged table
   would acquire a **second** reader (`naming::is_real`, `pub`, called from
   `borbax-ui`'s acceptance suite) — directly at war with this decision's own
   "exactly one reader" requirement below. And `REAL_ELEMENT_SYMBOLS` has
   exactly 118 entries while `element.rs` draws `n_elements` up to 120
   (`60 + next_range(61)` = `60..=120`) — a table reaching its 119th or
   120th element has no real symbol to give it.

   **Fix, structure**: keep `REAL_ELEMENT_SYMBOLS` and `REAL_WORDS` exactly
   as they are — flat `&[&str]`, so both existing `xtask` guards (which
   extract by exact source text and assert a floor) keep working unchanged.
   Add a private `REAL_TABLE: &[(&str, &str)]` beside them, with exactly one
   reader (the accessor below), and a **test** enforcing its correspondence
   with `REAL_ELEMENT_SYMBOLS`/`REAL_WORDS` rather than shared storage — the
   duplication is already deliberate, documented policy (`xtask`'s
   `REAL_WORDS_COPY`: *"duplicated, never imported, and the duplication is
   the point"*), and a test-checked third copy is consistent with that
   design, not a new sin. **The overflow case is resolved by construction,
   not by handling — and not through Decision 10's multiplicative mechanism,
   which `n_elements` does not go through** (it is a structural integer with
   its own cliffs; see P7's per-constant migration list). Instead, gate
   `n_elements`'s own existing draw directly on the shared rung: **draw and
   discard.** At `rung == 0`, still call `next_range(61)` — consuming
   exactly the same word off the stream every other rung consumes — but
   ignore its result and return **118**, the real, current element count;
   at any other rung, use the drawn value exactly as today. **Round 4
   correction: an earlier revision of this fix skipped the draw entirely at
   `rung == 0` rather than discarding its result.** `element.rs` draws
   `n_bands` immediately after `n_elements` off the same sequential stream,
   so skipping the draw would offset every subsequent draw by one word at
   the identity rung relative to every perturbed universe — exactly the
   "stream position becomes a function of caller data" hazard
   `borbax_rng`'s own documentation forbids at the method level, and
   exactly the "insert a draw, shift everything after it" hazard P7's
   stream-layout fix (below) exists to prevent, arriving here from the
   opposite direction: a removed draw rather than an inserted one. It also
   confounds Task 26.2 Step 3's independence test, whose whole premise is
   that perturbing one axis leaves every other axis's draws untouched —
   under the skip, *every* draw after `n_elements`'s slot differs between
   identity and perturbed universes, for a reason that has nothing to do
   with the axis under test. That makes the identity configuration's
   table size exactly `REAL_TABLE`'s size by construction; it never needs a
   119th or 120th entry. Perturbed universes
   keep the existing `60..=120` legal range untouched — 119- and
   120-element tables stay legal, they are simply never identity
   configurations, so the accessor below is never called for them.

   **Fix, guard mechanics — a four-part stack, not one AST check.** An
   earlier version of this decision said the real-name table's one call site
   would be "provably dominated by the identity check, verified by AST."
   Measured against 12 realistic call-site spellings, a straightforward
   `syn`-based dominance check **accepts only 2 of 7 legitimate spellings**
   (rejecting early return, `let`-else, `match`, `?`, and a compound `&&`
   condition as insufficiently "dominated") **and misses 1 of 5 leaks** (a
   shadowing case: `let strength = 0.0; if strength == 0.0 { real_name(z) }`
   passes as dominated while reading the wrong variable). Worse: dominance
   at the call site does not close the hole even where it works, because the
   identity classification is reconstructible **without ever calling the
   guarded accessor at all** — under `value = base * (1.0 + strength * p)`,
   every constant is bit-identical to its base exactly when `strength == 0`,
   so `if consts.some_constant == BASE_SOME_CONSTANT` recovers the identity
   predicate from outside the guard entirely, with no call to gate.
   "Dominated by the identity check" was never the right property to prove.

   The four parts, each closing what the others don't. **Round 4 correction:
   an earlier revision of this decision gave `real_name` two different
   signatures across parts 1 and 2 — a plain `Rung` argument in part 1, a
   required `&IdentityWitness` in part 2 — which are mutually exclusive, not
   two views of the same design. Only one can be real; the witness-taking
   form is the one that actually closes the hole, so it is the only one kept
   below.**
   1. **Move the identity test inside the accessor itself — by requiring the
      witness, not by re-checking the rung.** `fn real_name(_:
      &IdentityWitness, z: usize) -> Option<(&'static str, &'static str)>` —
      no `Rung` parameter. No caller can reach the table by any spelling,
      because there is no spelling that skips the check: the table is
      unreachable without a witness, and part 2 is what makes a witness
      unforgeable off-identity.
   2. **Make the predicate unforgeable by the type system, not just by the
      accessor's own logic.** `pub struct IdentityWitness(())` with a
      private field and one **fallible** constructor, `IdentityWitness::new(
      rung: Rung) -> Option<Self>`, returning `None` unless `rung == 0`. No
      witness exists off-identity, so `real_name` is dominated by
      construction — no witness, no call.
   3. **Replace the AST dominance check with two checks that are decidable
      and complete**: `REAL_TABLE` is named in exactly one file and read by
      exactly one function (inverting the default — an unrecognised shape
      fails closed, not open); and the number of `IdentityWitness(`
      **construction expressions** in the crate is exactly one (its own
      constructor) — a textual count of the substring `IdentityWitness(`
      cannot be this check, because the struct's own definition,
      `pub struct IdentityWitness(())`, textually contains it, so a textual
      count is at least two by construction and the check as stated could
      never pass. The count must be an expression-level `syn` count of
      `IdentityWitness(...)` call expressions, not struct/tuple-field
      syntax. Add a third: a `syn` check that `real_name`'s parameter list
      contains `IdentityWitness` and does not contain `Rung` — the check
      that would have caught this decision's own earlier contradiction.
   4. **Confine `Rung`'s visibility, so the identity predicate can't be
      recovered by reading the rung directly.** This is a G1-uniformity
      guard, not a second real-name check — see P6's restatement of this
      part, below, for why and for the `BASE_*` half of the mechanism.

   Store the **integer** rung as the identity coordinate (`rung == 0` is
   exact, no `-0.0` question, no `clippy::float_cmp` argument), not the
   derived `f64` strength. Required test, exhaustive over the whole rung
   space, on the **constructor** — `real_name` no longer takes a `Rung` at
   all, so there is nothing rung-shaped left on it to probe exhaustively:
   `for r in 1..=RUNGS { assert!(IdentityWitness::new(Rung(r)).is_none()) }`
   for every off-identity `r`. The compile-fail half of the guard is a
   separate probe, stated in P6 below, because a witness constructed off a
   non-zero rung is a `None` *value* at runtime, not a compiler error — the
   only compiler error available is constructing `IdentityWitness` from
   outside its own module. `xtask`'s `syn` dependency has no `visit` feature
   enabled — say in the guard's own doc comment which hand-rolled recursion
   it uses, the same discipline every existing AST guard in this file
   already follows.

10. **The identity configuration must be reachable across every
    physics-relevant constant, not just the element table — and the current
    per-constant drawing mechanisms can't do that.** `element.rs`'s own
    constants already use discrete ladders and are seed-reachable (confirmed
    by search: three identity seeds found in 21 seconds). `UniverseConsts`
    (§8.3's `ideal_gap`, `w_shape`, `w_charge`) and `bonds.rs`'s rate
    constants draw *continuously* — and for several realistic real-valued
    targets, the identity value has **zero** exact floating-point preimages
    in the drawn range, meaning no seed, however large the search, can ever
    reach it. **Fix, adopted rather than patched per-constant**: replace
    every physics-relevant constant's independent draw with one shared
    mechanism — draw a single `perturbation_strength` rung, `Rung(0..=m)`
    (`m` a tuned resolution), and apply it multiplicatively to every such
    constant: `value = base * (1.0 + strength * p)`, where `strength =
    rung as f64 / m as f64` is a non-negative magnitude and `p` is a
    per-constant, per-seed *signed* direction/scale drawn as before —
    direction lives entirely in `p`; the rung carries only "how far". At
    `rung = 0`, **every** perturbable constant lands on `base`
    simultaneously — the identity configuration — with probability
    `1/(m+1)` regardless of how many constants exist, which removes the
    per-constant-ladder design's combinatorial reachability problem
    entirely. This is new, cross-cutting infrastructure (Prerequisite P7,
    below) that Task 26.1's and Task 26.2's own constants must be drawn
    through, and that `UniverseConsts`/`bonds.rs` need a version-gated path
    to as well — V1 keeps its existing continuous draws untouched.

    **Rung arithmetic, stated once and referenced everywhere else, because an
    earlier revision of this decision contradicted itself here and no test
    could see it.** `Rung(0..=m)` has exactly `m + 1` members. Draw it via
    `next_range(m + 1)` — **not** `next_range(2*m + 1)`, which produces
    `2m + 1` outcomes and silently makes `P(identity) = 1/(2m+1)`, half the
    stated figure, while the plan's own text kept claiming `1/(m+1)` in the
    same paragraph. `next_range(m + 1)` is what actually delivers
    `P(identity) = 1/(m+1)`. Store the integer rung, not a derived float —
    `rung == 0` is exact by construction, with no `-0.0` question and no
    `clippy::float_cmp` argument, and it gives every downstream identity
    check (including P6's, below) one named, unambiguous coordinate to key
    on. The required test is a binomial-band assertion over repeated draws
    that pins `1/(m+1)` as a *number*, not an `O(m)`-order check — the `O(m)`
    form is satisfied by both the correct reading and the contradicted one
    and therefore cannot discriminate them.

    Determinism notes, otherwise unchanged: use `base * (1.0 + strength *
    p)`, never `base + base * strength * p` (the additive form produces NaN
    at `base = ±inf` and mishandles `-0.0`); canonicalise `-0.0` in `p`
    before any digest that hashes the perturbation vector.

11. **κ and γ are not drawn as two independently perturbed coefficients —
    their dimensionless ratio is, and γ is derived.** Task 26.2's original
    coefficient list drew `kappa` and `gamma` independently through Decision
    10's mechanism and then imposed `kappa >= gamma * t_max^(2/3) /
    (4*(t_max - 1))` as a post-hoc constraint. **Round 4 correction: the
    closed form this section derives its conclusions from was the argmin of
    the superseded `a^2` spelling, not the `a(a-1)` mechanism this document
    actually ships (above).** The correct closed form, argmin of the shipped
    energy over integer-relaxed `a` at fixed `total = t`, is
    `a*(t) = (4*kappa*t + gamma*t^(2/3)) / (8*kappa + 2*gamma*t^(2/3))` —
    computed as a ratio of sums of positive terms, never as a difference, so
    no cancellation risk; the equivalent form `a*(t) = t/2 − (t−1)/2·w` with
    `w = c·t^(2/3) / (2 + c·t^(2/3))` is algebraically identical and is
    stated here only because it makes the two facts below immediate, not as
    a second spelling to implement. At `t = 1` this form gives **exactly**
    `1/2` for every `(kappa, gamma)` — matching, and now actually
    supporting, the "provably, not by tuning" claim above. The old
    `a^2`-derived form instead gives `0.496...` at realistic coefficients,
    disagreeing with the true integer argmin on order-10% of totals across
    several independently sampled boxes — not a rounding-scale error, since
    the true argmin is exact under a strictly convex integer objective and
    the wrong closed form is systematically off, not merely imprecise.

    **The `a = 0` wall is closed by the `a(a-1)` spelling itself, not by any
    bound on κ and γ, and the old post-hoc bound's "278×–1071× too
    permissive" diagnosis belongs entirely to the superseded `a^2` model.**
    Under `a(a-1)`, `f(1) − f(0) = −4·kappa·(t−1)/t` exactly — negative for
    every `t > 1` and every `kappa > 0`, with `gamma` cancelling out of the
    comparison entirely, so the integer argmin never lands at the excluded
    wall regardless of κ or γ. (The old bound, `kappa >= gamma *
    t_max^(2/3) / (4*(t_max - 1))`, is exactly the `a^2` model's own wall
    condition — the two are algebraically identical — which is why it was
    evaluated at the wrong end of a strictly monotone family and measured
    278×–1071× too permissive *for that model*; a corrected version of the
    same bound has no job under `a(a-1)` and is not needed here.) The model
    is still one-sided where it's degenerate at *both* ends (as κ/γ → ∞ it
    returns to the rejected three-term form), and no clamp, rejection
    sample, or γ-dependent draw range for κ can be expressed through
    Decision 10's mechanism without breaking either its bit-identical-at-
    `rung=0` guarantee or routed requirement 5's pre-registration — that
    residual concern is what the reparameterization below actually answers.

    **Fix, converged on independently by every lane that looked at it, and
    justified by a stronger property than wall-avoidance**: don't constrain
    κ and γ separately. Draw the dimensionless ratio `c = gamma / (2 *
    kappa)` directly through Decision 10's mechanism, with its own
    `base_c`/`p_c` pair, and derive `gamma = 2 * kappa * c` after `kappa` is
    drawn. The corrected closed form above depends on `(kappa, gamma)`
    **only through `c`** — substituting `gamma = 2*kappa*c` and cancelling
    `kappa` gives `a*(t) = (2t + c·t^(2/3)) / (4 + 2c·t^(2/3))`, with no
    `kappa` left in it at all, confirming `kappa` is pure scale for this
    shape and `c` alone determines the valley's drift. That is the reason
    to draw `c` directly, independent of any wall — the two-sided legal
    range for `c` then holds **by construction** (no clamp, no rejection, no
    conditional draw order), and the identity configuration is preserved
    exactly, because `c` is bit-identical to `base_c` at `rung = 0` on
    exactly the same terms as every other migrated constant. The four Task
    26.2 coefficients become `(eps, sigma, kappa, c)`; `gamma` is computed,
    never independently drawn. **Digest `c`, not `gamma`, in P7's
    perturbation-vector digest** — `c` is the quantity actually drawn;
    `gamma` is derived from it, and hashing a derived quantity alongside its
    own inputs is redundant. State explicitly whether this digest is
    diagnostic-only or feeds a production golden hash before implementation
    starts, since that decides how much churn a later change to the
    derivation costs.

    **`t_max` in `c`'s legal range is a fixed design constant, not the
    per-universe table's actual largest total.** Using the drawn table's own
    size made the bound circular (Task 26.2 needs `c`'s range before the
    table that would supply `t_max` exists) and coupled the nuclear
    coefficients to Task 26.1's electronic axis through `n_elements` — so
    perturbing an electronic-tagged constant would move `n_elements` → the
    table's largest total → κ's legal range → every Q-value, which directly
    contradicts Step 3's own required independence test. A fixed constant,
    chosen once as a property of the model family rather than of any one
    table, removes the coupling entirely.

    **One more ambiguity resolved, not glossed over**: "the composition
    argmax" names two different questions depending on what is held fixed,
    and the plan's own text answers each one differently in different
    places. Step 1 arm 1 asks for "the distribution of optimal `a`/total …
    at a fixed total" — that is the argmax **over `a` at fixed `total =
    a + b`**, which is what the closed form `a*(t)` above computes, and is a
    global validation of the whole `(a, b)` energy surface's shape. Decision
    5's per-`Element` isotopic distribution is a different slice of the
    *same* energy function: fixed `a` (the element's own atomic number),
    varying `b`. Both are legitimate and both are needed — Step 1 validates
    the model before Step 5 uses one axis-aligned cut of it to build the
    per-`Element` distribution — but they answer different questions, and a
    future reader conflating them will expect Step 1's closed form to
    produce Decision 5's per-element table directly, which it does not.
    `t_max` in this decision's `c`-range derivation means the fixed design
    bound on `a + b`, never the largest atomic number (`units`) reachable in
    a table — Task 26.1 Step 6 is explicit that `units` is never a
    composition sum, and this decision does not change that.

    **Round 4 correction — the `t_max`/`n_elements` identity this decision
    used to pin `c`'s band is real but answers the wrong question, and the
    band it was supposed to produce cannot be derived from it at all.**
    Because the composition cost is strictly convex, exactly one `a` is
    **composition-stable** at each total under the corrected model, so
    `mean composition-stable isotopes per element == t_max / n_elements`
    exactly (verified at n = 30, 60, 120) — but "composition-stable" (Q ≤ 0
    on the composition-changing channel alone) is **not** the same quantity
    as Decision 7's "stable" (Q ≤ 0 across **both** the composition channel
    and the existing size-shedding channel). Arm 2, below, requires a
    positive-Q size-shedding transition to exist for heavy on-valley
    elements in most universes — so wherever arm 2 passes, the full-Q
    stable count is strictly **less** than the composition-stable count,
    and the identity above holds exactly at the same time as the true
    stable-count average falls short of it. The equality cannot be what
    pins `c`'s band, because it holds identically **regardless of `c`** —
    verified across a `c` range spanning three orders of magnitude,
    including `c = 0` itself, the degenerate point this whole reparameterization
    exists to leave. `kappa` cancels out of the composition argmax
    entirely (shown above), so this identity is uninformative about both
    `kappa` and `c`, and routed requirement 5's instruction to "derive `c`'s
    band from a chosen `t_max`/`n_elements` ratio" using this identity is
    **not executable** as written.

    **Fix: pre-register `c` from coverage, not from the mean.** The real
    constraint is that no element may have an empty isotopic distribution —
    Decision 5's abundance-weighted `mass`/`instability` averages are
    undefined (0/0, in library code) for an element the composition argmax
    never reaches. Since `a*(t) ≤ t/2` for every legal `c` (proved above:
    `a* = t/2 − (t−1)/2·w` with `w ∈ (0,1)`), coverage up to `n_elements`
    requires `t_max ≥ 2 * n_elements` as a hard floor before `c` is even
    considered — state this floor against the **legal** `n_elements`
    domain (P7's per-constant migration list, above), not the drawn value,
    so it stays non-circular. Given `t_max` fixed at or above that floor,
    pre-register `base_c` from the coverage boundary itself: choose it so
    the continuous curve's endpoint at `t_max` reaches the identity
    configuration's `n_elements = 118` — `base_c = (t_max / 118 − 2) /
    t_max^(2/3)` — and derive `c`'s legal band from the `p_c` range that
    keeps coverage across the whole **legal** `n_elements` domain, not just
    the identity value.

    **Two sub-decisions this fix surfaces and does not resolve on its own —
    flagging both explicitly rather than picking silently, because each
    changes the chemistry, not just the derivation:**
    - When the composition argmax overflows `n_elements` for a legal
      `(c, t_max)` pair, does the model **clamp** `a` to `n_elements`
      (piling every overflowing total onto the last element) or **drop**
      those totals (leaving them with no stable composition)? Nothing in
      this document says which, and the two produce different chemistry —
      not just different numbers.
    - `t_max ≥ 2 * n_elements_max` (against the *legal*, not drawn, ceiling)
      may conflict with keeping the mean `composition-stable` count inside
      any particular target range at the *identity* configuration — the two
      constraints are not guaranteed compatible for every legal
      `n_elements` ceiling. **This needs a decision before Task 26.2 Step 1
      is implemented**, not an assumption baked in here: narrow the domain
      one of the constraints binds over, or accept a defined fallback for
      the empty-distribution case above, or raise `t_max` and accept
      whatever mean results at the identity configuration.

    Decision 5's `ElementId` conclusion does not change (it is
    over-determined by four independent angles that don't depend on either
    figure at all) — only the derivation of its supporting number does.
    **Discriminating test**: sweep `c` over a decade at fixed `(t_max,
    n_elements)` and assert the composition-stable mean is **invariant**
    (it is `t_max/n_elements` regardless of `c`, by the identity above) —
    then assert that **coverage** (whether every element's distribution is
    non-empty) is the quantity that actually moves with `c`. A test that
    asserts the mean *selects* `c` is measuring the boundary rule from the
    first bullet above, not the physics.

    **Round 5, 2026-08-12 — pre-registration (routed requirement 5), fixing
    two more errors in the fix immediately above and settling both
    sub-decisions it left open.** Independently verified by
    `geometry-numerics-reviewer` (0/48,000 mismatches against brute-force
    integer argmax over `κ∈[0.1,100], c∈[1e-3,0.1], t∈1..400`).

    **The `base_c` formula stated two paragraphs up is itself the
    superseded `a²`-model inversion — the exact error round 4 announced it
    had fixed, one level further in.** `base_c = (t_max/118 − 2) /
    t_max^(2/3)` is algebraically the inversion of `a* = t/(2 + c·t^(2/3))`
    (the `a²` closed form), not of this document's own shipped `a(a−1)`
    form. It is **0.4237% low at every `t_max`** — small, but the same
    class of error this document has already had to retract once. The
    correct inversion of `a*(t) = (2t + c·t^(2/3)) / (4 + 2c·t^(2/3))` for
    the `c` at which `a*(t) = N`:

    ```
    c(t, N) = 2·(t − 2N) / ((2N − 1)·t^(2/3))
    ```

    **The recipe itself (pin `base_c` to the identity value exactly, at
    `n_elements = 118`) is unsound independent of that arithmetic — it
    leaves zero perturbation headroom.** At `t_max = 386` under the stated
    formula, `base_c = 0.024038`, giving `a*(t_max, 1.5·base_c) = 98.87` —
    every perturbed universe with `n_elements ≥ 100` would have empty
    isotope distributions for its heaviest elements. **Fix, replacing the
    recipe two paragraphs up**: pin `base_c` to the real physical value
    first (Decision 9: `base_c = a_C/(2·a_sym) = 0.711/(2·23.7) = 0.015`
    exactly, Rohlf's self-consistent pair — see the citation correction
    above), then solve `c(t_max, N) = base_c` for `t_max`, using the
    corrected inversion above and `N` set to the coverage target *with
    margin* (not 118 — see below), not the other way round.

    **`T_MAX = 386` (an integer — `t = a + b` is a nucleon count, never a
    continuous quantity, which the float `t_max ≈ 385.65` an earlier draft
    of this pre-registration used got backwards).** Chosen so that the
    worst-case reachable `c` (`base_c · 1.5`, at the crate-wide default
    `|p| ≤ Direction::P_MAX = 0.5`) gives `a*(T_MAX) = 121.08` — a 1.5-unit
    buffer above the bare `119.5` threshold `round(a*) ≥ 120` needs. The
    margin's job is landing on an integer `T_MAX` and headroom for a future
    added term (Decision 9's own text already anticipates one, e.g. a
    pairing term) — **not** floating-point precision: the spread across
    `t^(2/3)`'s three candidate spellings (`cbrt(t*t)`, `cbrt(t)²`,
    `powf(t, 2/3)`) is `2.84e−14` over `t = 1..400`, thirteen orders of
    magnitude below the margin, so citing precision as the reason would be
    a false rationale for a real number. Pin `t^(2/3)` as
    `det_math::cbrt(t*t)` regardless, per this document's own house
    convention (`element.rs:627` and the "every float spelling pinned"
    section above) — consistency with the runtime energy function, not
    defensiveness against this specific 386-cell scan.

    | at `T_MAX = 386` | `0.5·base_c` | `base_c` (identity) | `1.5·base_c` (worst-case) |
    |---|---|---|---|
    | `a*(T_MAX)` | 161.077 | 138.235 | **121.083** |

    Clears the hard floor (`T_MAX ≥ 2·n_elements_max = 240`) by 146.

    **Independent cross-check on the physics, not just the arithmetic**: at
    `c = 0.015`, this model's valley line against the textbook β-stability
    line `Z = A/(1.98 + 0.0155·A^(2/3))` — `A=16`: 7.659 vs 7.698; `A=56`:
    25.280 vs 25.375; `A=208`: 82.429 vs 82.404; `A=238`: 92.501 vs 92.417.
    Maximum deviation 0.10 protons across the whole range — independent
    support for `c_real = 0.015` beyond the citation fix alone.

    **Why coverage has no gaps — a proof, not the continuity hand-wave an
    earlier draft of this pre-registration used.** Two facts, both
    verified: (1) the per-unit energy is **exactly quadratic** in `a` at
    fixed `t` (leading coefficient `4κ/t + γ/t^(1/3) > 0`), so the integer
    argmax equals `round(a*(t))` exactly, never merely `{floor, ceil}` of
    it — confirmed 0 mismatches over the same 48,000-case sweep above.
    (2) `da*/dt` is strictly positive and **strictly less than 1/2**
    everywhere (measured supremum `0.499975` over `c ∈ [1e-4, 1], t ∈
    1..400`; the analytic supremum is exactly `1/2` as `c, t → 0`). A
    slope under `1/2` means `round(a*(t))` can step by at most 1 as `t`
    increases by 1 (with margin against adversarial tie-breaking), so no
    integer in `[round(a*(1)), round(a*(T_MAX))]` is ever skipped. Step 1's
    own coverage assertion should be the direct 386-cell scan this proof
    licenses (`for t in 1..=T_MAX`, assert every element `1..=n_elements`
    is hit at least once), not a re-statement of the argument.

    **The clamp-vs-drop sub-decision (the first bullet, two paragraphs
    up) is closed: drop — and the finding that closes it changes what
    "rare" meant.** At `T_MAX = 386`, `base_c` (identity), **67 of 386
    totals (17%) give `a*(t) > 118.5`** — no representable element in the
    118-element identity table — a routine fraction, not an edge case.
    (The earlier three-amigos pass that first proposed "drop" called it "a
    should-never-happen guard," which undersold this — at this `T_MAX` it
    is neither rare nor a guard, it is the expected outcome for roughly
    one total in six, at identity, before any perturbation is even
    applied.) **Drop is still the right choice**, on a real physical
    correspondence: no real stable nuclide exists above `A ≈ 209`, so a
    model that simply has no stable composition past a certain total
    matches the real periodic table's own boundary rather than fighting
    it. **Clamp would be visibly wrong at this `T_MAX`, not merely
    inelegant**: under clamp, the identity table's heaviest element would
    absorb 67+ isotopes against 2–4 for every other element — a gross,
    physically unmotivated distortion of the abundance-weighted
    `mass`/`instability` averages Decision 5 builds. Implement as a
    well-formedness filter on the candidate daughter/isotope set (per this
    document's own "domain-restrict... never a special case in the energy
    function" rule above), not a branch naming the outcome.

    **The second sub-decision (the `t_max ≥ 2·n_elements_max` floor
    possibly conflicting with a target mean ratio at identity) is now
    moot, not merely deferred**: nothing above targets a mean ratio.
    `T_MAX/n_elements = 386/118 ≈ 3.27` sits outside the `2.52–2.95` band
    cited earlier in this document (Decision 5, Decision 7, and routed
    requirement 5's own text) — that citation predates this section's own
    finding, a few paragraphs up, that the mean-composition-stable
    identity is **uninformative about `c`** (it holds regardless of `c`,
    including at `c = 0`) and therefore was never a real constraint to
    satisfy. `T_MAX/118 ≈ 3.27` is the *composition-stable* mean, a strict
    upper bound on the true post-both-channels stable-isotope count Step 1
    Arm 2 and Decision 7 actually care about (2.5–3.0) — the two are not
    in tension, since filtering through the size-shedding channel can only
    reduce the count, and by how much is exactly what Arm 2 measures
    empirically rather than what this derivation is entitled to assume.

    **The other three coefficients, pre-registered alongside `c` per
    routed requirement 5's "three things, not one interval" instruction —
    `base_i`, each coefficient's `p` range, and `m` together:**

    | constant | `base_i` | `p` range | `Counterpart` | source |
    |---|---|---|---|---|
    | `kappa` | `23.7` (MeV) | `±0.5` (`Direction::P_MAX`) | `HasReal` | Rohlf's `a_sym`, self-consistent with `base_c` above |
    | `c` | `0.015` | `±0.5` | `HasReal` | `a_C/(2·a_sym) = 0.711/(2·23.7)`, Rohlf |
    | `sigma` | `0.075` | `±0.5` | `DefaultOnly` | V1's own historical `sigma` range midpoint (`element.rs`'s `0.02 + 0.01·next_range(12)`, i.e. `[0.02, 0.13]`) — the scale Step 1 Arm 2 already needs for a like-for-like comparison against the V1 census |
    | `eps` | `1.0` | `±0.5` | `DefaultOnly` (composite — folds the real formula's volume *and* surface terms into one packing-derived contact count, no single clean real counterpart) | V1's own historical `eps` range midpoint (`[0.80, 1.20]`) |
    | `Rung::MAX` (`m`) | — | — | — | `99` — `P(identity) = 1/100 = 1%` exactly, a product decision (Ian, 2026-08-12): the search-cost case for a common identity seed is weak (found once, shipped as a named universe, Step 1b's precedent), so real chemistry stays a genuine rarity; `m=99` also gives every other migrated constant (Task 26.1's electronic axis included) a finer 100-rung perturbation ladder than the placeholder `m=8` did |

    No constant needed a custom-tightened `p` bound (unlike `w_shape`/
    `w_charge`'s `1/3`) — the crate-wide default closes every case checked
    above with margin.

    **One consequence worth stating plainly, found while pre-registering
    `kappa`: `gamma = 2·kappa·c` is derived, and at these real values
    `gamma = 2 · 23.7 · 0.015 = 0.711` exactly — `a_C` itself, which is the
    self-consistency check working as intended, not a coincidence.**
    `gamma/4 = 0.17775` **exceeds V1's own `sigma` range ceiling (`0.13`)
    outright**, not merely its midpoint as an earlier draft of this
    section anticipated qualitatively. Since the size-cost term's actual
    governing quantity on the valley is `sigma + gamma/4` (shown above,
    Step 1), the effective value at these real coefficients is
    `0.075 + 0.17775 ≈ 0.253` — roughly double `sigma`'s own historical
    ceiling. State this explicitly per the plan's own instruction: **Arm 2
    does not measure `sigma`'s growth-limiting behaviour in isolation
    against the V1 census; it measures `sigma + gamma/4`'s**, which at the
    real-physics coefficients sits well outside V1's own drawn range. This
    is not a defect — Decision 11's per-unit bound proof (`kappa +
    kappa/4` bounding all non-`sigma` cost, `sigma·t^(5/3)` supplying the
    only unbounded growth limit) holds at any positive `sigma`, including
    this one — but a reader expecting Arm 2's comparison to land inside
    V1's historical range at the real-physics point should not be
    surprised when it does not.

12. **The plan (this document, in its first revision) contradicted itself in
    eight places** — found by two review rounds' worth of specialists reading
    it against itself and against source. Each such spot in the sections
    below has been rewritten for consistency; the two structural sources were
    (a) claiming "zero `borbax-molecule` change" for a fix that provably
    touches `borbax-molecule` (F2), and (b) using `units` to mean two
    different things across the two tasks. Both are fixed by the specific
    remedies in Decisions 4, 5, and Task 26.1's wiring section.

## Prerequisites — land first, each its own commit, in the order given

**Files:** `xtask/src/main.rs`, `crates/borbax-universe/src/element.rs`,
`crates/borbax-universe/src/lib.rs`, `crates/borbax-universe/src/naming.rs`,
`crates/borbax-universe/src/bonds.rs`, and the actual (not assumed) location
of the shared test-fixture helper used across `borbax-molecule`'s test
modules — verify the real path before writing P3; an earlier draft of this
plan named a file (`testkit.rs`) that does not exist.

- [ ] **P1 — Re-key the closure-predicate guard correctly, not just wider.**
  The straightforward fix (add `unpaired`/`occupancy` to the existing global
  vocabulary lists) was tested against the actual guard code and caught 2 of
  13 realistic valence-promotion spellings — the vocabulary is *global*
  across the whole file, so `l == 0` (needed legitimately by `2(2l+1)` and by
  `block`'s derivation) becomes an accidental exemption on the exact function
  the guard exists to watch, and the likeliest real spelling
  (`unpaired.max(1)`) has no comparison operator for any vocabulary list to
  catch. **Fix**: don't rely on re-keying alone — per Decision 8, valence is
  now derived from a *named, specified* energy-gap computation rather than an
  ad-hoc branch, so P1's real job is verifying that specified derivation has
  no unaccounted branch, via a dedicated check scoped to that one function
  rather than a file-wide vocabulary scan. Probe with the actual likeliest
  spelling (`unpaired.max(1)`, not just an equality comparison) and confirm
  it's caught before this lands.
- [ ] **P2 — Pin V1's digest by name, not by `CURRENT`; derive the
  per-version digest list from an exhaustive `match`, not a hand-written
  array.** (Unchanged from the first revision — reviewed twice, found sound
  both times.) Probe: add a version arm byte-identical to V1's derivation,
  move `CURRENT`, confirm the assembled-universe digest test fails in a way
  that discards V1's pin before implementing the fix.
- [ ] **P3 — Give `generate_elements` a `PhysicsVersion` parameter, and give
  the returned `PeriodicTable` a `physics` field.** A parameter alone isn't
  enough: `BondEnergyMatrix::generate(&table, rng)` — the exact consumer that
  needs a version-aware path — takes the table, not the generation call, and
  has nothing to dispatch on without it. `PeriodicTable::new` is already
  `pub(crate)`. Probe: construct a table under one version and a bond matrix
  under another; confirm it does **not** compile once the field exists and a
  consuming match is exhaustive over it. Fix the real shared test-fixture
  helper in the same commit (find its actual path first) — it currently
  builds a version-blind table alongside a `CURRENT`-following universe,
  which would silently pair mismatched versions the day `CURRENT` moves.
- [ ] **P4 — Land immediately after Task 26.1 Step 8, not before it,
  because its subject doesn't exist yet.** Rescope the `block`-field read
  guard to allowlist by *call site*, not by crate, and make `block` a method
  (`Element::block()`, a pure function of `l`) rather than a stored field —
  this converts a field-read scan (which can't tell `if block == Block::D`
  from the semantically identical `if l == 2`) into a call-site guard, which
  is what the property actually requires. Allowlist exactly two legitimate
  call sites: construction (deriving `block` from `l`) and the viewer's
  display read; explicitly include reading `l` directly to locate period
  boundaries, since that's the first legitimate near-miss someone will write.
  Probe **inside** the function that computes valence (not in an unrelated
  placeholder) with `if block() == Block::D { valence += 2; }` and confirm
  the guard fires — a probe placed anywhere else doesn't discriminate the
  actual risk.
- [ ] **P5 — Run the full six-leg gate and commit** the prerequisites so
  far as their own unit.
- [ ] **P6 — Build the real-name guard (Decision 9), as the four-part
  stack, not a single AST dominance check.** Private `REAL_TABLE: &[(&str,
  &str)]` beside the existing `REAL_ELEMENT_SYMBOLS`/`REAL_WORDS` (which stay
  untouched flat lists), with a test enforcing correspondence rather than
  shared storage. One accessor, `real_name(_: &IdentityWitness, z: usize) ->
  Option<(&'static str, &'static str)>` — no `Rung` parameter; see Decision
  9's round-4 correction, which removes the contradictory `Rung`-taking form
  an earlier revision also gave it. The identity test lives **inside the
  witness's own construction**, not at this accessor — no call-site
  dominance to prove. `pub struct IdentityWitness(())`, private field, one
  **fallible** constructor, `IdentityWitness::new(rung: Rung) ->
  Option<Self>`, returning `None` unless `rung == 0` — unforgeable by the
  type system, not just by the accessor's own logic. `n_elements`'s base
  value (Decision 9) is drawn and discarded, then pinned to 118, so the
  identity table's size never exceeds `REAL_TABLE`'s. Three `xtask` checks
  on the real-name accessor, all decidable and complete rather than
  dominance-based: `REAL_TABLE` named in one file, read by one function
  (unrecognised shape fails closed); exactly one `IdentityWitness(`
  **construction expression** in the crate, counted at the `syn`
  expression level rather than textually — a textual count of the
  substring `IdentityWitness(` is at least two by construction, since the
  struct's own definition (`pub struct IdentityWitness(())`) contains it,
  and could never pass; and `real_name`'s parameter list contains
  `IdentityWitness` and does not contain `Rung`.

  **Round 4 correction — the fourth part is not a real-name check, and its
  key was porous on both readings.** An earlier revision filed "no `pub fn`
  outside the constant-drawing module reads a `BASE_*` constant" as the
  real-name stack's fourth part, reasoning that it closes the
  reconstruction path (`consts.x == BASE_X` recovers `rung == 0` without
  ever calling `real_name`). It does close that path — but the predicate it
  recovers is `rung == 0`, not a real name, so what it actually protects is
  G1's revised "not a special case in the code" clause and Decision 10's
  uniformity property, not G2. Sized against the real-name hazard it looked
  adequate; sized against the uniformity hazard it is the wrong key — and
  as stated it was porous on both counts regardless: a **private** helper
  reading `BASE_X`, called from a `pub fn`, passes a `pub fn`-keyed scan
  cleanly, and a **child** module reading its parent's private items via
  `super::` passes a file-keyed scan just as cleanly. **Fix, three parts,
  replacing the single scan:** (a) make every `BASE_*` constant
  `module`-private, not merely unread by a `pub fn` — the compiler then
  closes both the `pub use` path and the cross-module read path that no
  source scan can fully cover; (b) keep a file-keyed reader scan as a
  second, independent check (naming-convention drift is still worth
  catching even once visibility is the load-bearing guard); (c) record
  explicitly, in the guard's own doc comment, the residual **no** source
  analysis can close: a branch keyed on a **locally recomputed** base value
  (`let base = draw_constant(Rung(0), p); if consts.x == base { .. }`)
  needs no `BASE_*` identifier and no literal at all, and no scan can see
  it — the differential test below is what stands in its place. **State
  `Rung`'s own visibility explicitly here, rather than letting it fall out
  of an unrelated struct's convention** — Decision 9's `n_elements` gate,
  Decision 10's `strength` application in `element.rs`, `bonds.rs` and
  `UniverseConsts`, and P7's draw all need `Rung` visible across at least
  three modules, and `crates/borbax-universe/src/lib.rs`'s
  `pub struct Universe` gives every field `pub` by house style — a `rung`
  field held to that same convention would let `universe.rung ==
  Rung(0)` recover the identity predicate from `borbax-molecule`,
  `borbax-reaction` or `borbax-ui` in one line, defeating (a) and (b) from
  outside `borbax-universe` entirely. Keep `Rung`, and any field holding
  one, no more visible than the constant-drawing module actually requires.
  **Add one behavioural test alongside the structural ones**: generate the
  identity universe and a non-zero-rung universe whose `p` values are
  chosen so every migrated constant lands bit-identical to the identity
  universe's; digest both. Everything except the two universes' names must
  match exactly — a branch keyed on the rung itself rather than on the
  constants' values diverges here even when no source scan catches it,
  because this test does not depend on how the branch was spelled.

  Store the integer rung, not the derived `f64` strength. Probe: construct
  `IdentityWitness` from outside its own module; confirm it fails to
  compile — this is the guard's actual compile-time half, since a witness
  built off a non-zero rung is a `None` **value**, produced by a fallible
  constructor evaluated at runtime, not a type error the compiler can
  reject. Confirm separately, as a *runtime* test, that
  `IdentityWitness::new(Rung(r)).is_none()` for every `r != 0` (Decision 9's
  exhaustive test, above, restated here on the constructor rather than on
  `real_name`, which no longer takes a rung to probe). Probe: read a
  `BASE_*` constant from outside its module; confirm it fails to compile.
  Probe: read a `BASE_*` constant from a **child** module of the
  constant-drawing module via `super::`; confirm it still fails to compile
  — this is the case a file-keyed scan alone would have missed.
- [ ] **P7 — Build the shared `perturbation_strength` mechanism (Decision
  10).** One integer rung, `Rung(0..=m)`, drawn once per universe via
  `next_range(m + 1)` — **not** `next_range(2*m + 1)`, which gives `2m + 1`
  outcomes and halves `P(identity)` silently; `next_range(m + 1)` is the
  spelling that actually delivers `P(identity) = 1/(m+1)`. Store the integer
  rung; derive `strength = rung as f64 / m as f64` only where a magnitude is
  needed. Applied multiplicatively (`base * (1.0 + strength * p)`, never the
  additive form) to every migrated constant. Test: at `rung == 0`, every
  migrated constant is bit-identical to its `base` value, over the full
  range of finite `p`, at both zero signs of `p`. Test: a binomial-band
  assertion over repeated draws that pins `P(identity) = 1/(m+1)` as a
  number — not an `O(m)`-order check, which is satisfied by both the correct
  reading and the `2*m+1` defect and therefore proves nothing. Canonicalise
  `-0.0` in `p` before any digest hashes the perturbation vector.

  **Stream layout, stated as a decision, not left to the implementer.**
  An earlier version of this prerequisite said only "drawn once per
  universe," naming no `(Domain, index)` — every natural place to put it
  silently moves an existing universe (`Domain::Universe, 0` shifts every
  V1 element constant; index 1 shifts every `UniverseConsts` real; index 2
  shifts every bond-matrix cell), and `borbax-rng`'s own documented rule
  (*"new physics draws from a new `Domain`"*) exists specifically to
  prevent this. **Fix**: append `PerturbationRung = 10` to `Domain`; draw
  the rung from `Stream::new(seed, Domain::PerturbationRung, 0)` — its own
  `Domain`, not index `0` of a `Domain::Perturbation` shared with the
  per-constant `p` draws below. **Round 4 correction: an earlier revision
  reserved index `0` of one shared `Domain::Perturbation` for the rung and
  started per-constant indices from `1`.** Under today's continuous `p`
  spelling that collision is a doctrine violation with no measurable
  consequence (`next_f64` keeps the top 53 bits of the word; measured
  correlation between a rung sharing an index with a constant's direction
  and that direction's value is at the noise floor). But line 654's ladder
  spelling for `p` — "drawn on a discrete ladder at realistic resolution"
  — makes the shared index load-bearing: under it, `rung ≡ p_index (mod
  gcd(m+1, L))` with probability 1, and at the single most natural choice
  (`L = m+1`) the first constant's direction becomes **literally equal to
  the rung**, every seed, forever. A separate `Domain` removes the
  obligation to reason about this at all, rather than discharging it by
  choosing indices carefully — it inherits `every_domain_stream_is_pinned`
  free (the exhaustive `match` over `Domain::ALL` forces the new variant
  in, with its own pinned golden stream), and matches `borbax-rng`'s own
  precedent of a `Domain` naming a mechanism (`Domain::Hash`) rather than a
  location. Give each perturbable constant's `p` a **permanently assigned
  per-constant index** on its own `Domain::Perturbation` — an enum with
  explicit discriminants, appended and never renumbered — never a
  sequential draw off one shared stream, which would reintroduce the exact
  "insert a draw, shift everything after it" hazard
  `the_universe_digest_is_pinned`'s own doc comment exists to describe.
  Required probe: point the rung at `Stream::new(seed, Domain::Perturbation,
  0)` instead of its own `Domain`; confirm some existing test fails — as
  stated today, neither of P7's two required tests would notice.

  **Excursion bound and finiteness, stated as a mechanism, not implied by
  a measurement.** `1.0 + strength * p` is exactly `0.0` iff
  `strength * p == -1.0`, which annihilates the constant (`base * 0.0 =
  ±0.0`). **Round 4 correction: the mechanism that makes this unreachable
  is `strength ≤ 1` (proved below), not merely `|p| ≤ 0.5`** — an earlier
  revision stated the bound three times with three different premises
  (`strength ≥ 0` only; `|p| ≤ 0.5`; routed requirement 8's contradictory
  `p ∈ [-1, 1]`, which reopens exactly the annihilation case this bound
  exists to close, at `p = -1, strength = 1`), and the one that actually
  carries the "unreachable" conclusion appeared in none of the three.
  **Fix**: parse both `p` and the rung's own domain as finite-and-bounded
  by construction, not by a `debug_assert!` — which is absent from the
  `--release` profile `goldens --emit` runs in, and this repository has
  already shipped that exact debug/release asymmetry twice. `struct
  Direction(f64)` with `fn new(v: f64, bound: f64) -> Option<Self>` gated
  on `v.is_finite() && v.abs() <= bound`, with `bound` supplied
  **per-constant** — a single global `P_MAX` cannot express
  `w_shape`/`w_charge`'s tighter `1/3` bound below, so the constructor
  takes the bound as a parameter rather than reading one shared constant.
  `borbax_rng`'s own `next_f64_range` documents swallowing a non-finite
  bound (`lo = -inf` with finite `hi` returns `-inf` on every draw), so an
  infinite `p` is producible by a typo'd range and nothing today rejects
  it; a freshly generated NaN's sign differs by architecture, which would
  be a genuine §13.4 divergence, not a plausibility concern. Give `Rung`
  the same treatment: a private field and a fallible `Rung::new(r) ->
  Option<Self>` gated on `r <= m`, with `Rung::all_off_identity() -> impl
  Iterator<Item = Rung>` for Decision 9's exhaustive test to iterate rather
  than hardcoding `1..=RUNGS` (the plan states the rung's bound as both
  `RUNGS` and `m` in different places, and a literal range is satisfied by
  any `RUNGS`, including one that has drifted from the draw's actual `m`).
  This closes a real gap in the "unreachable" proof: at `rung > m` the
  factor's infimum is `+0.0` exactly (not merely small) once `rung = 2m`,
  and past `4m` a quarter of `p`'s range gives a **negative** factor —
  unbounded `Rung` construction, not just an unbounded `p`, is part of
  what the excursion bound has to close. **Proof that `strength ≤ 1`
  together with `|p| ≤ P_MAX = 0.5` makes annihilation unreachable**:
  `strength ∈ [0, 1]` by the domain just fixed, so `strength * p ∈
  [-0.5, 0.5]` and `1.0 + strength * p ∈ [0.5, 1.5]`, which never reaches
  `0.0`. **Drop the "reachable in roughly half of top-rung universes"
  estimate for annihilation under a discrete `p` ladder — it was never
  measured against a real ladder, and attaching a specific fraction to an
  unmeasured claim manufactures provenance routed requirement 6 exists to
  prevent.** State the actual, provable mechanism instead: any `p` ladder
  written in this codebase's own idiom (`lo + step * next_range(N)`, as
  `element.rs` already uses) contains its low endpoint by construction, so
  a ladder spanning `±100%` written that way contains `p = -1.0` exactly
  with probability `1/N` per constant — not by chance, by the idiom. At
  `strength = 1`, `1.0 + p` is then computed at `p = -1.0` exactly, which
  is error-free by Sterbenz's lemma (no nearby value rounds to zero), so
  the annihilation this bound closes is a real, exactly-reachable state
  under any such ladder — which is exactly why `strength` must be bounded
  to `1` and `p` to `P_MAX`, not a reason to estimate how often it would
  otherwise occur.

  **The per-constant migration list, with a reason for every inclusion and
  every exclusion — an earlier version of this prerequisite promised
  `element.rs`'s constants would migrate *and* promised Task 26.1 Steps
  14/17 that V1's digests stay unchanged, and both cannot hold at once for
  a migration with no version gate:**

  | constant | migrated through P7? | reason |
  |---|---|---|
  | `k` | **No.** | Structural integer — the shell law forces `z = k + 2`; `PackingConsts::new(k)` derives `z = shell_size(k, 1)`, and below `k = 6` a shell cannot triangulate a sphere. A continuous perturbation rounded back to an integer reintroduces the discrete-ladder mechanism P7 exists to remove. |
  | `n_bands` | **No.** | Structural `u8` count feeding `outer_fill_band`, which Task 26.1 Step 8 deletes — migrating a field that no longer exists is dead work. |
  | `n_elements` | **No.** | Structural integer with two undocumented cliffs under continuous perturbation: past 165 the naming grammar's own collision space is exhausted (measured 6.07 ms/mint, a 1100× slowdown, falling back to `Q{n}` names that violate `Element::symbol`'s documented 1–2 character contract, and no test today would fail on it); at 256, `PeriodicTable::new`'s `assert!` panics in library code. **Round 4 correction: two different intervals were both called this constant's "legal" range, self-contradictorily.** Its **drawn range stays `60..=120`, unchanged** (`60 + next_range(61)`, exactly today's mechanism); a *separate*, wider number — its **symbol-space ceiling, asserted at 165** — is the point past which the naming grammar's own collision space is exhausted, and exists only so a future widening of the *drawn* range is caught by an assertion naming the real cliff rather than silently exceeding it. The two numbers are not interchangeable: Task 26.2's isotope machinery (Decision 5) builds an abundance-weighted average over each element's isotopic distribution, and an `n_elements` draw whose composition ceiling (`a*(t_max)`, Decision 11) falls short of the drawn value leaves that element's distribution empty — a 0/0 average in library code. At `n_elements` confined to its drawn `60..=120` range this affects a measurable but small fraction of draws; were the drawn range itself ever widened toward the 165 ceiling, the empty-distribution fraction rises sharply, so Task 26.2's feasibility check (Step 1, below) must be written against the **drawn** range, never against the wider symbol-space ceiling. Its identity-configuration value (exactly 118) is handled separately, directly on the rung — see Decision 9's fix, above — precisely because this constant does *not* go through P7's multiplicative mechanism. |
  | `contact_defect` | **V2+ only, drawn through P7.** | Continuous on a `1/512` grid already; no structural constraint blocks multiplicative perturbation. |
  | `base_mass` | **V2+ only, drawn through P7, with its exactness constraint carried forward explicitly.** | `element.rs:486-494` documents, in bold and probed both ways, that `base_mass` must be **dyadic** — `base_mass * (1.0 + strength * p)` with continuous `p` is not. Mass conservation still holds (quantisation to `base_sub` precedes the integer multiply that actually matters for `Mass`), but the file's own recorded invariant does not survive unless the perturbed value is re-quantised to the dyadic grid immediately after the multiply, before anything else reads it. `MAX_BASE_SUB` (`element.rs:1083`, pinned tight against the current ladder's top, `1792 = 1.75 × 1024`) must be re-derived in the same commit against the new legal range, or `every_mass_is_inside_the_range_from_raw_does_not_check` silently stops checking anything. |
  | `eps`, `sigma`, `base_radius` | **V1 unchanged; V2+ replaces them.** | Fusion/packing-model constants that Task 26.1's screened-hydrogenic orbital fill does not use — migrating them is either a V1-affecting change (contradicting Steps 14/17) or work on constants V2 has already discarded. Task 26.1 draws its own screening coefficients, unrelated to these. |
  | `UniverseConsts`'s `w_shape`, `w_charge` | **V2+ only, drawn through P7, ratio-bounded — see below.** | Version-gated path added alongside the existing V1 continuous draw, which stays untouched. |
  | `ideal_gap` | **V2+ only, migrated for stream-layout uniformity — not because it is used.** | Zero production consumers anywhere in the workspace, verified by search. Migrating it costs nothing it wasn't already going to cost (one more permanently-assigned index) and keeps every `UniverseConsts` field on one mechanism rather than one migrated and one not; it is not migrated because anything downstream needs it to be identity-reachable. |
  | `bonds.rs`'s rate constants | **V2+ only, drawn through P7.** | Version-gated path added alongside the existing V1 continuous draws, which stay untouched. |
  | Task 26.2's `eps`, `sigma`, `kappa`, Decision 11's `c` | **V3 only, drawn through P7.** | New constants; no V1/V2 digest to protect. `gamma` is derived from `kappa` and `c` (Decision 11) and is not itself a P7 migration target. |

  This resolves Steps 14/17's contradiction by construction: nothing that
  feeds a V1 digest migrates. State in Task 26.1's own commit message
  which constants stayed and which moved, so a reviewer checking "V1
  digest unchanged" has the list in front of them rather than having to
  re-derive it from a diff.

  **Round 4 correction: don't leave a migrated constant's classifications as
  three separate prose lists — hang all three off the per-constant index
  enum as total, wildcard-free functions.** Each migrated constant needs
  three independent tags: which axis it belongs to (electronic, drawn on by
  Task 26.1; nuclear, drawn on by Task 26.2 — the independence Step 3
  requires depends on this being right), whether it `has a real
  counterpart` or is `default only` (below), and its legal `p` excursion
  bound (`P_MAX` by default, `1/3` for `w_shape`/`w_charge`). Stated as
  prose — as an earlier revision left them — a constant added later can be
  tagged correctly in two of the three lists and omitted from the third,
  and nothing catches it: the electronic-tagged set is a live instance of
  exactly this today, since Task 26.1's own screened-hydrogenic screening
  coefficients are required to migrate through P7 but appear in no P7 list
  by that name, and two of P7's own migrated names (`eps`, `sigma`) are
  shared with existing, unrelated V1 constants the migration list itself
  says stay untouched — a name alone cannot disambiguate them; the enum
  variant can. **Fix**: generate the per-constant enum, its `ALL` array,
  and the three total functions (`axis`, `p_bound`, `counterpart`) from one
  macro-driven table, with a `const _` compile-time assertion that every
  enum variant appears in `ALL` — a hand-written `ALL` array checked
  separately from the enum can omit a variant and still compile clean, so
  the assertion has to be structural, not a second list to keep in sync.
  Step 3's independence test derives its electronic/nuclear lists from
  `axis`, not from a hand-maintained name list.

  **`w_shape`/`w_charge` bound `p` as a ratio, not per-constant — the
  kernel only sees the ratio, and today's `[0.6, 1.4]` independent draws
  hard-bound `w_shape/w_charge` to `[3/7, 7/3] = [0.4286, 2.3333]` in a way
  nobody wrote down as a requirement, but that §8.3's kernel depends on.**
  Measured against the shipped prefilter tightness bar
  (`worst_gap < 0.98`, Decision 4): the safe zone is the *middle* of the
  ratio axis, and both ends fail — at ratio 0.0625 the size correlation
  reaches 0.9610, at 16.0 it reaches 0.8260, and today's independent
  `[0.6, 1.4]` draws keep every universe in the middle by construction,
  never by design. Migrating both constants through P7 with independent,
  unstated `p` ranges reopens exactly #31's defect and the size-comparison
  pathology Task 10's thirteen commits closed — measured at `|p| ≤ 1.0`
  (the natural "±100%" reading, and Decision 3 already works at ±50%):
  13.6% of universes leave today's support, 5% land where the measured
  size-correlation is 0.79–0.96. Nothing fails when this happens — every
  six-leg gate leg stays green — the affected universes are simply the
  boring, size-ranked ones. **Fix**: bound `p` for `w_shape` and
  `w_charge` at `|p| ≤ 1/3` each. **Round 4 correction: this does not
  reproduce today's worst-case ratio "exactly."** Today's independent
  `[0.6, 1.4]` draws hard-bound the ratio to `[3/7, 7/3] = [0.4286,
  2.3333]`; `|p| ≤ 1/3` gives `[0.5, 2.0]` — **tighter than today, not a
  reproduction of it**. Say so plainly, because "exactly" is what would
  later justify quietly widening the bound back to `2.3333` on the
  strength of a false equivalence.

  **Round 4 correction: the bound's justification should rest on the
  prefilter's own soundness, not on preserving a capability — the capability
  framing this bound inherits by citing "the shipped prefilter tightness
  bar" is a Principle 6 risk (tying a physics parameter's legal range to
  preserving a specific emergent behaviour) that the numerical argument
  below makes unnecessary.** The relative gap the prefilter tests,
  `gap(ρ) = (ρ(S−A) + (C−B)) / (ρS + C)` as a function of the shape/charge
  ratio `ρ`, is a Möbius transform of `ρ` and therefore **monotone** —
  verified over 200,000 random draws of its coefficients, 0 monotonicity
  violations, 0 violations of its own bracketing limits. As `ρ → 0` it
  tends to `(C−B)/C`; as `ρ → ∞`, to `(S−A)/S` — at each end of the ratio
  axis the two-channel bound degenerates to whichever single channel
  dominates, which is a description of the prefilter losing information,
  not of any chemistry becoming boring. That is sufficient on its own to
  justify bounding the reachable ratio range: a monotone bound that
  degenerates at its extremes needs its extremes kept away from, on pure
  numerical-soundness grounds, independent of what kind of universe results.
  Add a property test checking the *far end* of the ratio's reachable range
  against the prefilter bar — P7's two stated tests (bit-identical at
  `rung == 0`; `O(m)` search cost) are both about reachability, and neither
  looks at what the far end of a legal draw actually does to the kernel.
  Extend the existing `the_score_is_not_a_size_comparison`
  (`binding.rs:1076`) test's corpus over the reachable perturbed ratio
  range in the same commit — its current corpus is fixed seeds at V1's
  `[0.6, 1.4]` draws, so no perturbed universe is in it today, which is
  *why* every six-leg gate leg stays green over the degradation this bound
  exists to prevent. Re-measure the numeric bound on a V2 table after
  Decision 4's prefilter lands — the sweep above is on V1, where `affinity`
  never straddles zero and the charge channel is what carries the size
  proxy; #31 changes the curve's shape, not the conclusion that the ratio
  is load-bearing.

  **Round 4 correction: tag every constant whose value differs between the
  identity configuration and a perturbed one — not just every
  P7-migrated one — with `has a real counterpart` or `default only`, and
  correct what the identity name gate actually keys on.** The gate
  (Decision 9, P6 above) keys mechanically on `rung == 0` via
  `IdentityWitness`, not on this tag set — this section had said the gate
  "keys on" the `has a real counterpart` set reaching their real values,
  which overstated the mechanism's dependency on documentation; the
  accurate claim is that the gate's *honesty* — whether `rung == 0`
  actually means "the real periodic table" and not merely "every constant
  at some chosen default" — is only as good as this tagging being complete
  and correct. `n_elements` belongs in this set even though it is
  explicitly **not** P7-migrated (see the migration list, above): it is
  gated on the same rung, returns exactly 118 at identity, and the size of
  the identity table depends on it — leaving it out of the tag paragraph
  because it technically isn't "P7-migrated" is exactly the seam an
  earlier revision's narrower scoping fell into.

  As stated, this tag set is **empty** — no constant here has actually been
  tagged `has a real counterpart`. Complete it: `w_shape`/`w_charge` are
  weights in §8.3's own shape-complementarity kernel with no real physical
  analogue to measure against, so their `base` is necessarily a *chosen*
  constant in the sense `signature.rs:141-144` already warns about — tag
  both `default only`, chosen against the prefilter-soundness criterion
  above rather than by eye, and record the criterion in the constant's own
  doc comment. `ideal_gap` is `default only` by construction (zero
  production consumers). `contact_defect`, `base_mass`, and `bonds.rs`'s
  rate constants are `has a real counterpart` in principle but this
  document does not pin their real values — do so, or tag them `default
  only` explicitly rather than leaving them silently untagged. `n_elements`
  is `has a real counterpart`, pinned to 118 already. Task 26.2's `sigma`
  has no real-mass-formula counterpart (`sigma * t^(5/3)` is this model's
  own packing-series size term; no term in the real formula carries that
  exponent) — tag it `default only`, base value drawn from V1's own
  historical `sigma` range (`element.rs`'s `0.02 + 0.01 * next_range(12)`,
  i.e. `[0.02, 0.13]`), which Step 1 Arm 2 already needs as a like-for-like
  comparison scale. `kappa` and `c` **do** have real counterparts and must
  be pinned to them, not left as defensible-but-arbitrary points in
  Decision 11's legal band.

  **Citation correction, 2026-08-12 — the previous version of this
  paragraph cited "Rohlf" for a pair Rohlf's own text does not contain.**
  Verified against Rohlf's actual published coefficients (J. W. Rohlf,
  *Modern Physics from α to Z⁰*, Wiley, 1994, ISBN 0-471-57270-5:
  `a_V=15.75, a_S=17.8, a_C=0.711, a_sym=23.7, a_P=11.18` MeV) and
  cross-checked against Krane's *Introductory Nuclear Physics* and Imperial
  College's Dauncey lecture notes (both agree on `a_C≈0.71–0.72,
  a_sym≈23.0–23.2` for the classical fits): `base_kappa = 23.2` and
  `base_c ≈ 0.015323` were never a real published pair. `0.711/(2×23.2) =
  0.0153233` reproduces the old number to 7 digits, but it silently mixed
  Rohlf's own `a_C = 0.711` with `a_sym = 23.2` from a *different* fit (the
  "least-squares fit (1)" set, `15.8, 18.3, 0.714, 23.2, 12`) — two
  incompatible columns, not one source. No 2019 refit was found supporting
  `0.015323` either (checked Benzaid, Bentridi, Kerraci & Amrani, *Nucl.
  Sci. Tech.* 31:9, 2020, an actual AME2016 refit: `a_c=0.64, a_a=21.07` →
  `c=0.015187`, 0.88% off, not the "0.1%" this document previously claimed
  with no source named). **Fix: `base_kappa = 23.7` (MeV, Rohlf's own
  `a_sym`) and `base_c = a_C / (2 * a_sym) = 0.711 / (2 * 23.7) = 0.015`
  exactly** — Rohlf's own self-consistent pair, so both constants come from
  one fitted set rather than two. Published `a_sym` spans 21.1–23.7 MeV
  depending on the fit and mass evaluation used; "the real `a_sym`" should
  read "Rohlf's `a_sym`" throughout this document. The physics the
  coefficient expresses is unaffected: `c = a_C/(2*a_sym)` is exactly the
  coefficient in the β-stability line `Z = A/(2 + c·A^(2/3))`, and at
  `c=0.015` it predicts `Z=92.4` for `A=238` (real 92) and `Z=82.3` for
  `A=208` (real 82) — this was a wrong coefficient, not a wrong quantity.

  **The "band" claim above it, and the reasoning behind it, are both
  superseded — not by this citation fix, but by Decision 11's own round-4
  correction (below), found independently.** `t_max/n_elements ∈
  2.52–2.95` "pinning" `c`'s band via the mean-stable-isotopes identity was
  later shown to be **uninformative about `c`** (the identity holds
  identically regardless of `c`, including at `c=0`) — see Decision 11's
  round-4 text for the proof. `base_c = 0.015` still survives on its own
  terms: solving the *coverage* requirement (not the mean-identity one)
  with `base_c` pinned to `0.015` exactly gives `t_max ≈ 385.65`, at which
  the worst-case perturbed `c` (the crate-wide default `|p| ≤ 0.5`) still
  covers every legal `n_elements` up to 120 with a safety margin — see
  Decision 11's own text for the full derivation. That the real physical
  ratio and the coverage requirement are jointly satisfiable at all, with
  no custom-tightened excursion bound needed, is worth stating plainly:
  it did not have to work out that way. `eps` is the genuinely hard one to
  tag: `eps * contacts(total)` folds the real formula's volume *and*
  surface terms into one packing-derived contact count, so it has no
  single clean real counterpart — tag it explicitly as a composite and say
  so in its doc comment, rather than forcing it into either category; base
  value drawn from V1's own historical `eps` range (`element.rs`'s
  `0.8 + 0.05 * next_range(9)`, i.e. `[0.80, 1.20]`), same reasoning as
  `sigma` above. A `default only` constant's value at `rung == 0` is *a*
  chosen default, not *the* real value, and the plan should not imply
  otherwise.
- [ ] **P8 — Run the full six-leg gate and commit P6/P7 as their own
  unit.**

## Task 26.1: electron configuration — self-consistent, `l`-dependent orbital filling

**Files:**
- Create: `crates/borbax-universe/src/orbital.rs`
- Modify: `crates/borbax-universe/src/element.rs`, `lib.rs`, `bonds.rs`
- Modify: `crates/borbax-molecule/src/binding.rs` (the rearrangement-inequality
  prefilter, Decision 4 — this is the one real kernel-adjacent change; the
  scoring kernel itself is untouched), `signature.rs` (F2, per-species `SHELL`)
- Test: all modified/created files; `crates/borbax-ui` acceptance tests

**Interfaces:**
- Consumes: `borbax_units::{Mass, Span, Quanta, canonical_cmp}`,
  `borbax_rng::{Stream, Domain}`, P7's `perturbation_strength` mechanism
- Produces: `orbital::{Occupancy, OrbitalConsts, fill(consts: OrbitalConsts) -> Vec<Occupancy>}`,
  `PhysicsVersion::V2`, `Element::block() -> Block` (method, not field,
  replacing `outer_fill_band: u8`), `PeriodicTable::physics: PhysicsVersion`

### Orbital fill, screened-hydrogenic with `l`-dependent screening

- [x] **Step 1: Write the failing tests.** Subshell capacity `2(2l+1)`,
  answering correctly for `l` up to 6 so it can't be faked as a lookup table.
  Hund's-rule unpaired count spanning `l ≥ 1`. Order-invariance under a
  shuffled candidate list, sourced from a `borbax_rng::Stream`. The full
  subshell-closure set as running partial sums, never a hardcoded noble-gas
  subset. **The screening function's `l`-dependence, tested directly**: at
  the identity configuration, assert the real fill order through at least
  5p (a single deterministic array, built as the output of running the
  formula, discussed with Ian as a §5-permitted assertion under Decision 2 —
  it is not a lookup table, it is the *expected output of a computed
  formula*, same status `packing.rs`'s Euler-forced constant already has).
  Assert this fails under an `l`-blind screening implementation — plant one,
  confirm the test catches it (this is Blocker 7's actual regression test;
  the diversity floor alone does not catch it, see below). `Z_eff >= 1`, with
  equality permitted at Z=1 (**not** `Z_eff > 1` — that's false at the
  element this epic exists to fix), provable from a one-sided-downward
  screening-coefficient draw (`σ ∈ (0, 1]`, closed upper bound) rather than
  merely observed. The screening sum's accumulation: build the same final
  occupancy via two different insertion orders, assert `Z_eff.to_bits()`
  equality (the shuffled-candidate test above does *not* exercise this,
  since the argmin is unique under the tie-break and the fill order doesn't
  change under a shuffle).
- [x] **Step 2: Verify each test fails for the stated reason.**
- [x] **Step 3: Implement the fold — explicit, not a sort — with an exact,
  order-free screening accumulation.** `Z_eff` depends on current occupancy,
  built via a repeated argmin over unfilled subshells given occupancy so far
  (`canonical_cmp(e, best_e).then_with(|| (n,l).cmp(&best))`, never
  `min_by(|a,b| a.partial_cmp(b).unwrap())`). The screening-`l` function
  enters as `t(l) = (0, 0.5, 1, 1)` (or an equivalent smooth function that is
  ≈0 at l=0 and saturated by l=2) — **not** the Slater step, which is less
  robust to perturbation (measured: 8% d-block retention at ±50% vs this
  form's 98–100%). Keep occupancy as `u32` counts (exact, order-free) per
  screening class; reduce the screening sum to a **fixed-arity float
  expression applied once**, association written out in a comment — this
  makes "accumulate from the fill's own emitted order" and "a fixed number of
  pinned float terms" the same requirement (they were stated as two
  incompatible ones in the first revision) by ensuring what's accumulated is
  an exact integer count, with the float expression applied only at the end.
  **Ban incremental `z_eff -= sigma` update.** No `HashMap` anywhere in this
  path. Draw the screening coefficients through P7's `perturbation_strength`
  mechanism, one-sided downward on the inner-shell coefficient.
- [x] **Step 4: Digest the derived quantities with the raw bit function, not
  the canonicalised one.** `Z_eff` and `radius` (both float) must be hashed
  via raw `to_bits()`, matching this crate's own existing documented
  rationale for its detector digests — canonicalisation is lossy by design
  and would collapse a genuine `-0.0`-sign regression onto `+0.0`, reporting
  "unchanged" over a real defect. `Mass`'s `canonical_bits()` is fine to use
  as-is (it's lossless by construction — `i64` has one representation per
  value). Also hash the emitted fill order itself, since `block`/`period`/
  `group` are functions of it.
- [x] **Step 5: Run all Step 1 tests.**

### Wiring: `Element`, `PeriodicTable`, `PhysicsVersion::V2`, and the shared "gap" derivation

- [x] **Step 6: Write the failing tests.** All four of `valence`,
  `affinity`, period boundaries, and V2's bond capacity derive from the
  HOMO-LUMO-analogue gap (Decision 8) — assert this as a shared computation
  with one call site each, not four independently-tuned formulas. `radius`
  falls across a period and jumps at each boundary, with boundaries read
  from the gap-derived computation, never a literal list. `affinity`
  periodic in Z with at least one within-period sign change and a pre-fixed
  correlation bound against Z. Valence promotion: beryllium's analogue
  (closed subshell, but promotable within the energy budget) gets nonzero
  valence; helium's analogue (closed subshell, no accessible promotion within
  the budget) gets zero — both from the same rule, no per-element exception.
  Bond-matrix rank-2 residual correlating with the affinity *difference*
  specifically. `n_elements`'s last element sitting at a subshell boundary.
  `block()` unreadable outside its two allowlisted call sites (P4).
  `Z_eff`/`radius`/`mass` finite and `>= 1`/positive respectively, asserted
  at creation. **`units` keeps one meaning across this task and Task 26.2 —
  atomic number, never redefined as a composition sum** (the plan's first
  revision used it both ways; Task 26.2 uses a distinct field for total
  composition size, never overloading `units`).
- [x] **Step 7: Verify failure for the stated reasons.**
- [x] **Step 8: Implement the wiring.** `Element`'s stored fields:
  `symbol, name, units, period, group, mass, valence, affinity, radius,
  energy_per_unit, instability, abundance` (dropping `outer_fill_band`
  entirely — no replacement stored field; `block()` is a method on `l`).
  `generate_elements(seed, PhysicsVersion::V2)` returns a `PeriodicTable`
  carrying `physics: PhysicsVersion` (P3). `ShellPattern`'s V1-only member
  (`k`) is split out into a version-specific type (`#[non_exhaustive] enum
  ShellLaw`) so the compiler forces every consumer to handle both, rather
  than a V2 table carrying a meaningless `k`. `BondEnergyMatrix::generate`
  gains a V2 path keyed on the gap-derived bond capacity from Decision 8 —
  this is what makes bond strength correlate with the same structure #31's
  fix uses, rather than an independently invented capacity model.
- [x] **Step 9: Run all tests from Steps 1 and 6. Add one binding test that
  actually runs on a V2 table** (not just V1, which every existing fixture
  currently pins to via P3) — `worst_shape < 0.90`, `worst_gap < 0.98`
  (post-prefilter-fix), `score <= 0 && is_finite`, so the shape/charge
  balance and the prefilter's tightness are exercised before `CURRENT` ever
  flips, not discovered afterward.

### The prefilter fix (Decision 4) and F2

- [x] **Step 10: Write the failing tests.** The rearrangement-inequality
  prefilter: sound (bound ≥ actual score) on the full existing corpus;
  strictly tighter than the current bound; `worst_gap < 0.98` restored on a
  table where affinity has been fixed to straddle zero; self-pair
  discrimination restored (a self-pair's bound must not degenerate to
  unrejectable). F2: per-species `SHELL = 0.9 * this molecule's own mean
  atom radius` (not `≈0.90` — pin the exact value, and document that `f < 2`
  is the geometric requirement "the channel never sees through the atom
  facing it"), and the property it must satisfy is *composition-within-one-
  table* invariance — not "flat across molecule size" (already true of the
  shipped constant, unrelated to F2) or "flat across seeds" (delivered by
  radius normalisation, unrelated to F2). Test: rescaling every radius in a
  universe by a constant leaves the character channel's per-direction
  weights unchanged (scale-equivariance) — this is the property F2 actually
  owns. **Before writing this test, measure the V2 table's actual radius
  max/min** — "broken under V2's wider range" was an assumption in the first
  revision, not yet a measurement.
- [x] **Step 11: Verify failure.**
- [x] **Step 12: Implement.** Rearrangement-inequality prefilter in
  `binding.rs` (sort each side's shape/charge channels via `canonical_cmp`,
  dot the sorted sequences — `SigSummary` stops being trivially `Copy` at
  this size). Per-species `SHELL` in `signature.rs`, computed once at intern
  time from the molecule's own atoms (§8.6-compliant — no per-pair or
  per-step recomputation).
- [x] **Step 13: Run all `borbax-molecule` tests. This step DOES regenerate
  `borbax-molecule` goldens — say so, and say why.** The prefilter change and
  F2 both touch `borbax-molecule` directly (F2 specifically **cannot** be
  version-gated: `signature()` takes no `Universe`/`PhysicsVersion` by a
  deliberate prior design decision, so per-species `SHELL` changes scoring
  for V1 chemistry too). This is a real, deliberate physics-adjacent change
  affecting geometry/shape scoring universally, not an artifact of V2 landing
  — measured to move the winning pose on 2.3–5.5% of pairs (against 0.7%
  for the radius-normalisation fix that *is* gatable, and the 54.9% that got
  an earlier, rejected kernel-internal approach thrown out). Name this
  explicitly as the deliberate regeneration event, distinct from Task 26.2's
  separate one.

### Verification and commit

- [x] **Step 14: Confirm `borbax-universe`'s V1 digests are unchanged**,
  under `generate_under(seed, PhysicsVersion::V1)`, using P2's name-pinned
  digest.
- [x] **Step 15: Point the viewer at a V2 seed** via a visible on-screen
  selector (the physics version already travels in the universe address as
  `@2`) — never an `env::var`. Confirm the row count matches the gap-derived
  period boundaries as a *consistency* check, not "does it look like 7 rows."
  Confirm real names appear only at the exact identity configuration
  (`rung == 0` — the integer coordinate, per Decision 9/10's fix; not a
  float `perturbation_strength` comparison), generated names everywhere
  else, and that
  `borbax-ui`'s existing "never atomic number, never Z, never protons" string
  is updated to reflect that under V2, `units` genuinely is the atomic
  number — that shipped test needs a deliberate edit, not silent survival.
- [x] **Step 16: Rewrite `binding.rs`'s RAF-specificity comment properly,
  with an actual citation, not one sentence.** **Correction to this bullet's
  own first draft, made during implementation and verified against both
  papers directly:** the "same bound whether catalysis is flat or scales
  with length" result is Mossel & Steel (2005)'s Theorem 4.1(ii), not
  Hordijk, Wills & Steel (2014)'s — HWS extends *that* result to catalysis
  depending on the reaction too, and separately studies their own extreme
  case (MLEN, only maximum-length molecules catalyse) as a still-RAF-viable
  but measurably costly degenerate case. Neither paper requires size-free
  binding; neither paper's core theorems are confined to the fully
  scale-free case — Mossel & Steel's own theorem covers both the flat and
  the length-proportional case with one bound. State what RAF closure
  actually is sensitive to (the mean catalysis level crossing a required
  threshold) and that the size signal here is a *retained, reviewed carrier
  of catalytic specificity*, not a literature mandate. State the good news
  too: closing #31 moves this codebase's binding kernel into the
  independence regime those theorems assume, out of the size-degenerate
  regime it was in before.
- [x] **Step 17: Run the full six-leg gate and commit.** State explicitly:
  `borbax-universe` V1 goldens unchanged; `borbax-molecule` goldens
  regenerated deliberately (Step 13 — affects V1 chemistry too, this is
  intentional); #31 closed structurally via the gap-derived `affinity`; the
  prefilter's tightness bound restored via the rearrangement inequality;
  #25 is **not** closed by this task — see Task 26.2.

## Task 26.2: nuclear structure and isotopes — a local valley of stability

Built from `docs/superpowers/plans/2026-07-30-borbax-isotopes.md`. Lands as
`PhysicsVersion::V3`. Its original three-term mechanism was found to be
mathematically degenerate (the imbalance term's optimum sits at exactly the
symmetric point for every drawn coefficient, every size, every universe) —
the corrected mechanism below (Variant B) was arrived at after two reviewers
proposed different fixes and a cross-check between them reversed the
initially-recommended one.

**The mechanism, corrected:**

```
binding      eps   * contacts(total)
size cost    sigma * total^(5/3)
imbalance    kappa * (a - b)^2 / total
asymmetry    -gamma * a*(a - 1) / total^(1/3)      — NEW: the term that lets the valley drift
```

**`a*(a - 1)`, not `a^2` — the pair-counting SEMF spelling, and the difference
is load-bearing, not stylistic.** An earlier revision of this mechanism used
`a^2`, under which composition `(a=1, b=0)` — the lightest element, at exactly
the size issue #25 is about — is exoergically unstable (`Q = gamma > 0`) **in
every universe, for every legal `(kappa, gamma)`**, directly contradicting
Decision 7's own motivating example (*"hydrogen-1 … is nonetheless perfectly
stable"*). `a^2` counts each unit's self-interaction with itself; `a(a-1)` is
the standard pair-counting form that removes it, and it makes the `t=1`
optimum exactly `1/2`, so `Q ≡ 0` identically there, **provably, not by
tuning** — verified over 20,000 random `(eps, sigma, kappa, gamma)` draws: 0/
20000 unstable under `a(a-1)`, against 20000/20000 under `a^2`. Valley drift
(`a*/total`) is unchanged by the substitution.

Three drawn coefficients (`eps`, `sigma`, `kappa`) plus the dimensionless
ratio `c = gamma / (2 * kappa)`, each through P7's `perturbation_strength`
mechanism (Decision 10) — `gamma` itself is **derived** as `2 * kappa * c`
after `kappa` is drawn, never independently drawn or independently
constrained (Decision 11). This replaces an earlier, defective form of this
paragraph that drew `kappa` and `gamma` separately and imposed `kappa >=
gamma * t_max^(2/3) / (4*(t_max - 1))` as a post-hoc bound — that bound was
evaluated at the wrong end of a strictly monotone family and was measured
278×–1071× too permissive at the plan's own stated value; see Decision 11
for the full derivation and the reparameterisation that replaces it. `c`'s
two-sided legal range holds the model's integer-rounding degenerate case
closed **by construction**, not by luck of the draw.

**Term names, corrected to match the physics they approximate — this costs
nothing and prevents a future "fix."** The plan's `asymmetry` term
(`gamma * a*(a-1) / total^(1/3)`) is the real Bethe–Weizsäcker formula's
**Coulomb** term; the plan's `imbalance` term (`kappa * (a-b)^2 / total`)
is the real formula's **asymmetry** term. The names in the diagram above
and everywhere this document refers to them should say so, or be renamed to
match the literature (`coulomb`/`asymmetry`) — pick one and apply it
everywhere, so a future reader citing a real SEMF coefficient maps it to
the right term on the first try.

**Every float spelling pinned here, in the plan text, not left to
implementation time — an earlier version of this section pinned the
trivial one (`(a-b)^2` as `d * d`, correct, required since `f64::powi` is
denied) and said nothing about where the bits actually move:**

- `total^(5/3)` as `total * det_math::cbrt(total * total)` — cube-root-of-
  the-square, matching `element.rs:627`'s existing house convention, not
  `det_math::powf(total, 5.0/3.0)` (differs in 247/250 sampled totals,
  worst 4 ulp, and carries 2.3× the ulp error against a reference value —
  `5.0/3.0` as `f64` is not exactly `5/3`, so `powf` evaluates a slightly
  different exponent than the one written) and not `cbrt(total).powi(5)`
  or `cbrt(total.powi(5))` (each differs from the pinned spelling on more
  than half of sampled totals). Tabulate `cbrt(total * total)` once per
  integer total and reuse it — the closed-form validation formula below
  (`a*(t) = (4*kappa*t + gamma*t^(2/3)) / (8*kappa + 2*gamma*t^(2/3))`,
  Decision 11, used by Step 1) needs the same `t^(2/3)` value at the same
  `t`, and computing it twice by different spellings would silently
  reintroduce the divergence this pin exists to close.
- `total^(1/3)` (the Coulomb/`asymmetry` term's denominator) as
  `det_math::cbrt(total)` directly — a different computation from
  `cbrt(total * total)` above (`cbrt(t*t)` and `cbrt(t)` squared differ on
  115/250 sampled totals), not interchangeable with it.
- `let d = a as f64 - b as f64;` before squaring, never
  `(a - b) as f64` — `a` and `b` are unsigned compositions, and the
  subtraction-before-cast form panics under `overflow-checks` in
  `cargo test` and silently wraps in `cargo test --release`, the exact
  debug/release asymmetry CLAUDE.md already records shipping twice in
  Task 2.
- **`a * (a - 1)` itself, added in round 4 — the section's own headline
  term was missing from a list introduced as "every float spelling
  pinned."** `let af = a as f64; af * (af - 1.0)`, never `(a * (a - 1))
  as f64` on the unsigned product — `a` is an unsigned composition and
  `a - 1` underflows at `a = 0` under the same debug-panics/release-wraps
  asymmetry as the `(a - b)` bullet above. This term is exact under every
  spelling for the domain that matters (`a * (a-1) ≤ 14,280` for `a ≤ 120`,
  far below `f64`'s 53-bit mantissa — 0 bit mismatches across the
  integer-product, cast-first, and `af*af − af` forms, checked over
  `a ∈ 1..=200`), so this bullet is about **totality and sign**, not
  accumulation order like the others above: pin the evaluation order so
  `a = 0` is filtered by `reachable_compositions` (below) **before** the
  energy function is ever called on it, not after — the shipped channel
  never evaluates this term at `a = 0`, but Decision 11's own closed-form
  validation and Step 1's own probe scans do reach `a = 0` deliberately,
  so both must use the exact, cast-first spelling above rather than the
  unsigned one. Note for the one-signed-expression form above: `af * (af -
  1.0)` at `af = 0.0` gives `-0.0`, which the surrounding negation and sum
  erase before any digest reads it — provided the bare per-term value is
  never digested on its own.
- Left-to-right term association for the four-term sum, explicitly, per
  `element.rs:624-625`'s existing precedent and vocabulary
  ("parenthesised to preserve the left-association exactly … (§13.4)") —
  grouping the two negative terms first differs from left-to-right
  association on 34.6% of sampled compositions. The **division by `total`**
  is already correctly pinned as applying to the whole sum, not per-term
  (dividing per-term differs on 65.6% of compositions) — this bullet only
  pins the sum's own internal order, which was previously unstated.
- Write the four-term binding-energy expression as **one signed
  expression**, not four rows with role labels and an implicit sign on one
  of them. An earlier version of this section's diagram carried an
  unlabelled `-gamma` on only the fourth row; read as *additional* to an
  implicit positive convention on the first three, the second derivative
  of the per-unit energy can go negative whenever `gamma * total^(2/3) >
  4 * kappa`, which is a materially **weaker** convexity condition than
  the one Decision 11 derives (119× weaker at `total = 120`) — Step 1's
  own gate would have caught the resulting degenerate model, so this is a
  wasted cycle rather than a shipped defect, but writing one signed
  expression removes the ambiguity before it costs either.

**The tie-break stays, and its justification is corrected a second time —
the round-3 correction was itself wrong, in the direction that matters
most: an exact tie does still occur, at the one composition #25 is about.**
"Near-symmetric compositions still tie exactly" was true of the rejected
three-term model, where `kappa*(a-b)^2/total` alone is symmetric under
`a <-> b`. Round 3 corrected this to "an exact tie no longer occurs" under
the asymmetric four-term model — **false, at `total = 1`.** Under
`a(a-1)`, both legal compositions at `total = 1` — `(a=1, b=0)` and
`(a=0, b=1)` — give `a*(a-1) = 0`, so the Coulomb-analogue term vanishes
for both, and `(a-b)^2 = 1` for both — the two costs are **identical**,
verified bit-identical (`0x4037333333333333` both sides) across 50,000
draws spanning six decades of `kappa` and `c`, with no other exact argmin
tie found in 116,000 further `(coefficients, total ≥ 2)` cases. **This is
also why the "endoergic" wording below (Step 7's commit-message paragraph)
is wrong**: the real transition at `total = 1` is **isoergic** — `Q = 0`
exactly, sitting on the clamp rather than under it, not "endoergic"
(`Q < 0`) as stated there. **Keep the tie-break on §13.1's float-key-argmax
discipline alone** — it costs nothing regardless of whether a tie is
reachable. **Round 4 correction: `gamma = 0` is not independently
reachable, and the citation for that claim ("H8") does not exist anywhere
in this document** — delete it. Under P7's excursion bound (`strength ≤ 1`,
`|p| ≤ P_MAX = 0.5`), both `kappa` and `c` stay within `[0.5, 1.5] ×` their
base values, so `gamma = 2 * kappa * c ≥ base_kappa * base_c / 4 > 0` in
every reachable universe — `gamma = 0` is not a live case this tie-break
needs to guard against. **What actually re-arms the exact-tie case is the
`total = 1` structural tie above, present in every universe regardless of
coefficients — and it is closed by the domain restriction (`a ≥ 1`,
below), not by the tie-break.** Without the domain restriction, the
comparator (`(a, b).cmp(&best)`, below) would select `(0, 1)` — a
zero-proton element — over the physically sensible `(1, 0)`; with it, the
tie is moot because `(0, 1)` is never a candidate. The tie-break remains
correct defensive practice for any future coefficient or term this section
does not yet have, which is reason enough to keep it even though it is not
what closes the `total = 1` case. Break ties on the composition argmax
with the same comparator rigour Task 26.1 Step 3 already specifies —
`canonical_cmp(e, best_e).then_with(|| (a, b).cmp(&best))`, never
`min_by(|x, y| x.partial_cmp(y).unwrap())`, which would be the first
`.unwrap()` this crate ships in library code. Probe: force `gamma = 0.0`
and confirm the tie-break fires — this remains a valid unit-test probe
even though real universes never reach `gamma = 0`. State the sampling
regime for the "no other exact tie" measurement above, per routed
requirement 6. State the model in **per-unit** form (the total-energy form
above is off by a factor of the total size) — divide the whole expression
by `total` before it's assigned.

**Round 5 correction, 2026-08-12 — the line above used to end "…it feeds
`energy_per_unit`", and that is wrong: fixed here rather than left to
contradict the resolution below at Step 5 implementation time.**
`energy_per_unit`/`ShellPattern::eps` is what both `bonds.rs`'s V1
(`contact_density`, "the positive half of `energy_per_unit`, divided by
`eps`") and V2 (`base = pattern.eps * scale`) paths price every bond from —
reusing it for this task's nuclear binding model would make bond strength a
function of nuclear structure, and real nuclear binding energies are
~2.25×10⁶× a real chemical bond energy (Ni-62's BE/A peak, 8.794555
MeV/nucleon, against a C–C bond at 376.98 kJ/mol) — no real channel above
~10⁻⁵ of a bond energy connects the two. Caught 2026-08-11/12 by an ad-hoc
three-amigos session after surviving four `/review-plan` rounds, then
independently confirmed by two review-agent lenses in parallel (physics
magnitude; alife emergence — a correctly-rescaled coupling would still be
wrong, because real per-nucleon binding is nearly flat above A≈12 while
real bond energies span 7–11×, so feeding one into the other would flatten
the whole bond-energy landscape and kill differential persistence, §9.4).
**This task's four-term per-unit binding model feeds a quantity Step 5
names when it lands — not `energy_per_unit`, not `ShellPattern::eps` —
which in turn feeds `instability`'s Q-value redefinition (Decision 7),
`mass`'s defect term, and eventually `abundance` (the Fe-peak in cosmic
abundance **is** the nuclear BE/A peak, NSE, Ni-56→Fe-56 — legitimate
physical backing for that route specifically).** One coupling stays legal
and is worth keeping deliberately, not rejecting wholesale: composition →
`mass` → the Arrhenius prefactor (`A ∝ 1/√μ`, the kinetic isotope effect) —
real, small (1–2%, from the H₂/D₂ dissociation-energy comparison), and
already permitted by Decision 5's independence test, which names
valence/affinity/radius/block as the electronic-side observables, not
mass. **Step 3's independence-test observable set must add
`BondEnergyMatrix::energy`** (`bonds.rs`) to what it asserts bit-identical
under a nuclear-tagged perturbation — as originally planned, the test
suite has no observable that could see this exact coupling leak back in
by accident, which is precisely the shape of miss that let the plan's own
wording say "feeds `energy_per_unit`" for two review rounds running.

**The rejected alternative, recorded so it isn't tried again**: replacing
the size-cost term with an asymmetric one instead of adding a fourth term
(fewer coefficients) was measured to delete the growth limit entirely — at
the optimum, its size cost falls to 0.3–21% of the original, and a heavy
on-valley element has no interior binding peak at all in 20–78% of draws.
It also turned out to be structurally incompatible with the Q-value
redefinition below (that redefinition's whole argument rests on having two
distinct decay channels; the rejected alternative deletes one of them).

**`instability` is redefined as a Q-value (Decision 7), not a global
distance.** Best-available-transition energy release, clamped at zero,
across both a composition-changing channel (this task) and the existing
size-shedding channel (V1's mechanism, kept). The clamp is arithmetic —
zero is where a transition stops releasing energy — not a tuned threshold.
**Round 4 addition: spell the clamp as `if q > 0.0 { q } else { 0.0 }`,
never `q.clamp(0.0, f64::MAX)`.** `f64::clamp` is not on `clippy.toml`'s
denied-methods list (unlike `f64::max`/`f64::min`, which are, precisely
for this reason) and passes a negative zero through unchanged
(`(-0.0f64).clamp(0.0, f64::MAX)` returns `-0.0`), while the `if`-form
returns `+0.0` at `q = -0.0`. Whether a `-0.0` Q is actually reachable
depends on the sign convention chosen for `Q`, which this section leaves
open below (does the channel conserve `total`, or shed an `a`-unit) — spell
the clamp safely regardless of which convention is chosen, rather than
deferring the question until after the convention is picked.
`max(Q, 0)` over the candidate transitions needs **no tie-break of its
own** — a fold-max over any finite set of `f64` values (including
duplicates, `±0.0`, and `NaN`) is order-invariant by IEEE `maxNum`
semantics, and this crate already forbids the one spelling
(`max_by(partial_cmp().unwrap())`) that would turn it into one. **Round 4
addition — say why explicitly, because the reason is stronger than "no
hazard found": at this granularity, the winning *daughter* is not even a
well-defined question.** `Element.instability` is one abundance-weighted
average over the *whole* isotopic distribution (Decision 5); the
`(daughter_composition, Q)` pair belongs to one isotope's transition, not
to the element, so there is no single "winning daughter" for an `Element`
to expose even in principle — only the scalar survives averaging. **This
constrains Task 15, which is not itself in scope here but inherits a
requirement from this decision**: Task 15's drafted `decay_channels`
mechanism needs a real state-change vector per channel (a Gillespie
consumer needs `(a_μ, ν_μ)` — a propensity *and* a product state, not a
propensity alone), and the per-isotope `(daughter_composition, Q)` set
this task computes internally cannot be reconstructed from the scalar
`instability` this task exposes. Task 15 will need the per-isotope set
itself, not a derived average — record that as an open requirement routed
forward, not solved here. The `Q → propensity` conversion needs one scale
constant, shared with Task 15's known-open "radiogenic channel has no
scale constant" item — one golden regeneration serves both. Decide
deliberately whether `instability` becomes `Quanta`-typed or stays a
Q-derived `f64` with the energy scale kept internal to the nuclear module
(the latter avoids a `Quanta` leaking into a Gillespie-propensity consumer
that doesn't want one). **State explicitly, for whichever task first wires
this into a stepping beaker**: the composition-changing channel enters
simulation only as a per-species propensity computed once at intern time —
no per-molecule isotope state, no per-molecule damage field, per §8.6.
Decay is this project's own recorded most-common disguise for exactly that
regression, and nothing in this task's scope currently states the
boundary explicitly.

**The composition-changing channel is bidirectional, and the plan's earlier
"a composition-changing channel" (singular, undirected) understated what
that requires.** Measured: element 1's stable-isotope count is 3
bidirectionally against ≥200 unidirectionally, under every variant including
the `a(a-1)` fix above — direction is not a property the energy function
has, and a unidirectional channel undercounts stable isotopes by orders of
magnitude by silently imposing one. Evaluate `Q` both ways (`(a,b) →
(a+1,b-1)` and `(a,b) → (a-1,b+1)`) and take the best-available transition
over the union. State explicitly whether the channel conserves `total`
(β-like — the size terms cancel identically) or sheds an `a`-unit; the plan
does not currently say, and the two give different `Q` formulae.

**Domain-restrict the channel to compositions the table can represent, and
spell the restriction as a well-formedness condition on the daughter set —
never as a special case in the energy function.** `for d in
reachable_compositions(parent)` derived from the table's own bounds, with
`Q = max` over that set, is a domain fact: delete it and the channel emits a
daughter the element table cannot represent. `if a == 1 { return 0 }`
computes the same number at the one composition that matters and **is
rigging** — it names an outcome (element 1 must be stable) rather than a
domain boundary (no element has `a = 0`). This distinction is why Step 7's
commit message can honestly say #25 closes because stability became local:
under `a(a-1)` the domain restriction is no longer load-bearing for the
headline result (real H-1 is stable because the relevant transition is
**isoergic** — `Q = 0` exactly, sitting *on* the clamp rather than under
it; see the tie-break paragraph above for the exact-tie this produces at
`total = 1` — not because a daughter is unlisted). This is a knife-edge
zero, not a strict negative, and has no margin against a future added term
or a changed float spelling — worth stating plainly rather than leaving
the reader to infer it from "stable." It only prevents the channel from
proposing a daughter outside the table, which it must do
regardless of what makes element 1 specifically stable.

**Round 4 addition: state, as a bound on what may be claimed rather than as
a to-do list, what the identity configuration can and cannot reproduce.**
G1/G2's 2026-08-07 revision permits real-world correspondence and
calibration at the identity configuration; the model's fidelity there is
genuinely strong — at the real semi-empirical mass formula's own
coefficients, the valley line this model produces matches the textbook
β-stability line `Z = A/(1.98 + 0.0155·A^(2/3))` to well under half a
proton across the whole table range. It does **not** reproduce every real
structural fact: the model omits the neutron–proton mass-difference and
pairing terms the real formula carries, so it does not reproduce even-`A`
double stability or the real distribution of isotope multiplicities per
element. Adding either term later is not rigging under the revised G1 — it
is real physics, exactly the distinction the domain-restriction fix above
already draws (`if a == 1 { return 0 }` names an outcome; a well-formedness
condition names a boundary) — but it is **not free**: any added term
changes the energy surface, which requires re-deriving `c`'s coverage band
(Decision 11, above) and re-running Step 1's gate, not simply appending a
correction. **One condition on this list, stated so a future change trips
it rather than passes silently**: this bound belongs in the docs, never as
a literal array in the crate. If any of these correspondence facts is ever
encoded as, say, `const REAL_STABLE_ISOTOPE_COUNTS: [u8; 118]` inside
`borbax-universe`, that is a G1 breach in kind — `check_no_data_files`
cannot see an inline source-level array any more than it can see a
checked-in file, and no existing guard would catch it.

**Files:**
- Modify: `crates/borbax-universe/src/element.rs` (new composition fields on
  `Element`, per Decision 5 — a distribution, not a new `ElementId` axis),
  the nuclear energy model (new module or extend `packing.rs`'s V3
  counterpart), `bonds.rs` if `instability`'s consumers need updating
- Test: same files, plus any `borbax-molecule` fixture that reads `mass`
  (the one field composition affects there) to confirm it's still a simple
  scalar under the distribution model, not a new species dimension

**Interfaces:** `ElementId` is **unchanged** — this was closed as a decision
before any implementation, not left for Step 1 to inform (see Decision 5).
`Element` gains composition-distribution fields; `mass`/`instability` are
computed as abundance-weighted averages over that distribution.

### Does the mechanism actually produce a valley — settle this before building anything else

- [x] **Step 1: Run the corrected model's own gating probe, with a second
  arm, sampled the way the shipped code will draw.** **Round 4 correction:
  the first arm's original framing — "is the most-bound composition
  interior, and specifically not always at the symmetric point" — is a
  tautology of the corrected model's own algebra, guaranteed for every
  positive coefficient, and Step 3 below already bans this exact assertion
  shape ("composition is interior… must not be treated as [a]
  discriminating test") while this step still asks it as its gate.**
  `a*(t) = t/2 − (t−1)/2·w` with `w = c·t^(2/3)/(2 + c·t^(2/3)) ∈ (0, 1)`
  for every `c > 0`, so `a* ∈ (1/2, t/2]` always — verified over 200,000
  sampled `(c, t)` pairs spanning nine orders of magnitude in `c`: the
  minimum observed `a*` is exactly `0.5`, the maximum observed
  `a* − t/2` is exactly `0.0`. No coefficient draw can put the optimum at
  the boundary or outside the symmetric point; the question the first arm
  literally asks cannot fail. **What the first arm should actually gate on
  is drift, not interiority**: at small `t` the symmetric composition is
  *physically correct* (light nuclei genuinely cluster near N≈Z), so a
  flat "minimum drift in integer units" criterion would fail a correct
  model — report the **drift-onset total**, the smallest `total` at which
  `round(a*(t))` first departs from the symmetric composition, as the
  scalar this arm actually gates on, and require it below a pre-registered
  fraction of `t_max`. That number moves several-fold across a plausible
  `c` range even though the continuous curve barely does — at the identity
  configuration, elements `a = 1..7` sit exactly at the symmetric
  composition, and halving `c` roughly doubles the drift-onset total; a
  gate on the continuous curve's shape alone cannot see this at all. Report
  the *distribution* of optimal `a/total` across universes at a **fixed**
  total, not pooled across totals (pooling across totals makes a degenerate
  model look diverse — the original, broken model produces 60 distinct
  ratios that way while having exactly one optimum at every individual
  size). **Report it as `a*/total` — a curve in `total`
  within each universe, comparing curve shapes across universes — not as
  the spread of a single scalar at one `total`.** One drawn ratio genuinely
  produces one curve (the real SEMF has exactly one Coulomb/asymmetry ratio
  too, so this is physically correct, not a shortfall of the model) — a
  gate that reports the spread of `a*/total` at one fixed `total` across
  universes cannot distinguish genuine valley diversity between universes
  from the single knob each universe already has; "the valley drifts" is a
  claim about a derivative with respect to `total`, and only the curve
  form can support it. **Second arm, required because the first alone
  doesn't discriminate a model that's interior but has lost its growth
  limit**: does per-unit energy along the valley still peak *inside* the
  table (compare against the existing V1 census — on a closure in most
  universes), and does a full search still find a positive-Q transition
  for a heavy on-valley element, in what fraction of universes. **`sigma`
  and `gamma` are collinear on the valley — pre-register this before
  reading arm 2's result, don't discover it after.** On the valley, `gamma`
  contributes an addition to `sigma` of `gamma/(4*(1+u)^2)` (`u` the
  imbalance ratio), which in the curvature-dominant regime — the only
  regime with realistic drift — is flat at `gamma/4`; the two coefficients
  are not separately identifiable from size-cost behaviour alone, the
  `total`-only effective coefficient is always `sigma + gamma/4`. Arm 2
  compares the binding peak against the existing V1 census, drawn under
  `sigma` alone from `0.02..0.13`; under this corrected model the *actual*
  governing quantity is `sigma + gamma/4`, and at realistic `gamma` values
  that addition can exceed `sigma`'s own range midpoint. The comparison is
  only like-for-like if `gamma/4` stays small against `sigma`'s range —
  state that explicitly as a pre-registration constraint, and do not read
  arm 2 as measuring `sigma`'s growth-limiting behaviour in isolation; it
  measures `sigma + gamma/4`'s. **The growth limit itself rests entirely on
  `sigma`, provably — stated per-unit, and completed to cover the
  Coulomb-analogue term the original version of this proof omitted.**
  **Round 4 correction: stated in *total* form, as an earlier revision did,
  the assertion is false for every sufficiently large `total`** (on-valley
  imbalance-term cost grows `∝ kappa * total`, unbounded) — the proof only
  holds **per-unit**, which this paragraph must say explicitly since it is
  the form that actually feeds `energy_per_unit`. Per-unit, the on-valley
  imbalance-term cost (`kappa * (a-b)^2 / total^2`, i.e. `kappa * u^2` in
  the imbalance ratio `u`) is bounded above by `kappa` — it approaches
  `kappa` as `total` grows, never exceeds it. **The Coulomb-analogue term
  is not automatically bounded and must be shown, not skipped**: per-unit
  on-valley Coulomb cost has the closed form `(gamma/c) * x / (2+x)^2`
  with `x = c * total^(2/3)`, and `x/(2+x)^2` is maximised at `x = 2`,
  where it equals `1/8` exactly — giving a per-unit Coulomb bound of
  `gamma / (8*c) = kappa / 4` (using `gamma = 2*kappa*c`). Verified against
  the maximum observed value at several `(kappa, c)` draws: `5.7923` vs the
  bound `5.8000`, `3.7363` vs `3.7500`, `8.7456` vs `8.7500` — the bound
  holds and is close to tight. So total non-sigma on-valley per-unit cost
  is bounded by `kappa + kappa/4`, and only `sigma * total^(5/3)` (still
  stated in total form here, since it is what feeds the size-cost term
  directly) is unbounded and can supply a growth limit against a saturating
  binding term — adopt this as the stated reason the corrected model
  cannot lose its growth limit the way the rejected alternative did,
  replacing the earlier regime-dependent "0.3–21% / 20–78%" figures, which
  needed a stated sampling regime this proof does not. **Pre-register
  `base_i` (each coefficient's base value),
  each coefficient's `p` range, **and** `m` together, as three things, not
  one interval** — P7's shared rung correlates all four Task 26.2
  coefficients through one `strength`, so `kappa`'s reachable range given
  the others is **narrower** than under independent draws and collapses
  toward `base_kappa`/`base_c` at low rungs; any degeneracy figure measured
  under independently-sampled coefficients (including this document's own
  earlier "98–100% of draws out to ±50%" and "8% d-block retention"
  figures, both measured independently of P7) does not transfer to the
  regime the shipped code will actually draw from. Sample `(rung, p_eps,
  p_sigma, p_kappa, p_c)` jointly, the way `generate_elements` will. **If
  either arm fails, stop — do not proceed to wire this into `Element`.**

  **Round 6, 2026-08-12 — both arms pass; cleared to proceed to Step 2.**
  Full account: `docs/experiments/2026-08-12-nuclear-valley-gate.md`,
  implementation: `crates/borbax-universe/src/nuclear.rs`. Three corrections
  landed before the gate could run meaningfully, each found by an
  independent verification pass rather than trusted on the first
  derivation: (1) `composition_argmax` ships the κ-free closed-form
  spelling, not the one this section's own text calls "the correct closed
  form" two paragraphs up — the two differ by up to 4 ulp and round to a
  different *element* in 3 of 6,951,474 sampled cells; (2) `BASE_EPS =
  2.625`, not `1.0` — the original value sat at an unrelated scale next to
  `kappa`/`c`'s real Rohlf MeV values, leaving 74.6% of totals unbound at
  identity; (3) the nuclear packing model's fixed coordination parameter is
  `k = 10`, not `12` — this crate's shell law gives `z = k + 2`, so `12`
  actually selects coordination 14. Both (1) and (2) are also recorded in
  `perturbation.rs`'s own `MigratedConstant` doc comments.

  Arm 1's gate is the analytic envelope `drift_onset(c) ∈ [16, 30]`,
  derived from `c`'s own P7-reachable range rather than this section's
  Ca-40 anchor (`40/386`) — measured to be structurally unfalsifiable,
  since `kappa` provably cancels and `c` is bounded before any sampling
  happens, so no implementation could ever fail it. Arm 2 passes the
  literal criteria above (100% closure against the 89.44% bar) plus an
  added bound-set-shape check (at most a handful of sign transitions in the
  on-valley binding curve, tolerating the one measured shell-closure
  "island of stability" at `total = 309`). One further finding, reported
  rather than gated: 9.25% of universes cover fewer than 120 elements while
  still bound — real, and now corrected into `MigratedConstant::C`'s own
  doc (its coverage proof covered the composition argmax's shape, not
  whether binding stays positive that far out) — routed to Steps 2+, which
  own how it bears on `n_elements`'s draw.
- [ ] **Step 2: Confirm `ElementId`'s shape is unchanged** (it is, per
  Decision 5). **Round 4 correction: "two molecules differing only in
  isotopic composition" is not constructible under Decision 5's own design,
  which makes the required cross-species test as originally stated
  tautological — it cannot fail under any implementation, correct or not.**
  Under Decision 5, isotopic composition is a property of the *universe*
  (the element's abundance distribution), never of a molecule — `Mol12`
  has one `ElementId` per atom and no isotope-shaped field for two
  molecules to differ in, so there is no constructor by which the required
  pair could exist; the comparison degenerates to a species compared with
  itself, and `group_distance(s, s) == 0` is a theorem of the signature
  construction (identity is always in the rotation group; the accumulator
  is a sum of squares whose floor is zero), not a property this task's
  design contributes. **Fix, two parts:**
  1. **Make per-atom isotope state unrepresentable**, with a `const _`
     compile-time size assertion on `Mol12` (the pattern already used
     elsewhere in this crate for a fixed-layout invariant) pinned against
     today's measured struct size — the assertion's own comment must name
     what it actually guards ("no additional per-atom field of adjacency
     size was added"), not "no isotope field," since the assertion trips on
     *any* qualifying addition and a comment naming only isotopes would
     mis-describe what future failures mean.
  2. **Rescope the behavioural test from two molecules to two universes**:
     build the *same* `Mol12` in the identity universe and in a
     non-zero-rung universe whose nuclear-tagged `p` values are chosen so
     every migrated constant still lands bit-identical to the identity
     universe's (mirroring Step 3's own vacuity fix, below — assert the
     perturbation actually landed on some other observable *before*
     asserting `group_distance == 0` and bit-identical `affinity` here,
     or this test inherits the same vacuity-at-the-wrong-rung failure mode
     for the same reason). This test is falsifiable today: it would fail
     the moment `embed` or `signature` starts reading `Element::mass` or
     `Element::instability`, which is exactly the §8.6 boundary Decision 5
     asserts and which currently has no test guarding it at all.

### Wiring

- [ ] **Step 3: Write the failing tests.** #25 actually closes: the smallest
  element in every drawn universe has `instability == 0` under the Q-value
  definition (not "a valley exists," which is true by construction and
  proves nothing — see below).

  **The composition-distribution independence test is vacuous at
  `rung == 0`, and must not run only there.** At `rung == 0`,
  `base * (1.0 + 0 * p) == base` for every `p` on every constant — P7's own
  headline test asserts exactly this — so "perturbing any nuclear-tagged
  constant leaves `(valence, affinity, radius, block())` bit-identical" is
  trivially true at the identity rung **regardless of whether the two axes
  are actually independent**, and so is its converse. A fully coupled
  implementation passes this test unchanged if it is only run at the
  identity. **Fix**: state the rung the independence test runs at — it
  must be `> 0`, and several distinct rungs, not the identity — and
  **assert the perturbation actually landed** (the perturbed constant
  differs from its base value) before asserting the other axis did not
  move; a test that only checks the second half would pass just as
  vacuously if the perturbation silently failed to apply. Perturbing any
  nuclear-tagged constant (`eps`, `sigma`, `kappa`, `c` — not `gamma`,
  which is derived from `kappa` and `c` and is not an independent
  perturbation target, per Decision 11) leaves `(valence, affinity,
  radius, block(), BondEnergyMatrix::energy)` bit-identical; perturbing
  any electronic-tagged constant leaves the composition distribution and
  Q-values bit-identical. **`BondEnergyMatrix::energy` added to this
  observable set 2026-08-12, per the RESOLVED bond-pricing decision
  above** — without it, the test suite has no way to see whether a
  nuclear-tagged perturbation leaks into bond pricing (via
  `energy_per_unit`/`ShellPattern::eps`) by accident; every other
  observable in this set is silent to exactly that failure mode.

  Mass conservation: composition changes under a reaction still sum `Mass`
  exactly. **Both of the original note's headline assertions ("composition
  is interior," "instability at the valley is less than at the extremes")
  are tautologies of the corrected model's own definition and must not be
  treated as discriminating tests** — replace them with Step 1's actual
  gating measurements (pre-registered ratio spread, pre-registered
  peak-location and Q-positivity fractions) reported as instrumentation,
  not asserted as pass/fail on values the model's algebra already
  guarantees.
- [ ] **Step 4: Verify failure for the stated reasons.**
- [ ] **Step 5: Implement**, per Decision 5's `ElementId`-unchanged
  structure and Decision 7's Q-value.
- [ ] **Step 6: Regenerate goldens deliberately** — the largest physics
  change in the epic. State in the commit message which V1/V2 behaviour is
  V3-superseded (the size-shedding decay channel gains a composition-changing
  partner) and which stays available by physics version. This is `V3`'s
  digest, separate from Task 26.1's `V2` digest (Decision 6) — confirm both
  are independently pinned.
- [ ] **Step 7: Run the full six-leg gate and commit.** State explicitly:
  #25 is closed because stability became local, not because a valley
  appeared — the composition axis is what makes the local definition
  *expressible*, not what makes it true.

## Routed requirements

1. **Every corpus-sweeping test's examined-count must be asserted and
   derived** (sum of table lengths), not a round number chosen before
   measuring.
2. **Plant the literal fixed version back into wherever the viewer reads
   the physics version, after V2/V3 land, and confirm the vacuous-today
   acceptance test now fails.**
3. **Amend §7.1 in the same PR, covering all four passages this epic
   touches**: the shell-pattern-lengths line (energy-based filling's
   reachable period lengths are a small structured set, which that line
   currently forbids), the one-axis-limitation paragraph (which explicitly
   forecloses "adding a neutron count" — Task 26.2 is exactly that, and its
   authority, G3-absolute, is revised), the `instability` row (needs to
   describe the Q-value shape), and the `outer_fill_band` row (delete —
   `block()` is a method, not a stored field with its own spec row).
4. **Task 26.2's Step 1 is a hard, two-armed gate on the rest of Task
   26.2.** Do not let "the corrected model's algebra already guarantees
   this" substitute for actually running both arms against the real V2
   electronic-side numbers.
5. **Pre-register `kappa`'s drawn range and `c = gamma/(2*kappa)`'s two-sided
   legal band, in their own commit, before Task 26.2 Step 1 runs** — `gamma`
   is derived (Decision 11) and is not independently registered. Derive
   `c`'s band from a chosen `t_max`/`n_elements` ratio (Decision 11's H7
   fix: mean stable isotopes/element equals that ratio exactly under the
   corrected model, and Decision 5 already commits to a 2.52–2.95 band), not
   from `kappa` and `gamma` separately. Every quantitative claim about the
   corrected model, in this document and in the review that produced it, is
   a function of these ranges, and none of them are specified yet. **P7
   correlates all four of Task 26.2's coefficients through one shared rung,
   which narrows `kappa`'s reachable range given the others below what an
   independent-draw pre-registration would show — so pre-register `base_i`
   (each coefficient's base value), each coefficient's `p` range, **and**
   `m` together, as three things, not one interval, and re-run Step 1's
   gating probe sampling `(rung, p_eps, p_sigma, p_kappa, p_c)` jointly, the
   way `generate_elements` will draw them** — a probe sampling the four
   coefficients independently measures a regime the shipped code does not
   use.
6. **State the measurement regime alongside any figure requoted from this
   plan's review history** (e.g. the 96.3%/31.7%/68.3% figures, any
   κ/γ-dependent number) — several are correct only under a specific,
   previously-unstated perturbation regime, and restating them without that
   context risks a future reviewer "correcting" a true claim into a false
   one. This now explicitly includes this document's own "98–100% of draws
   out to ±50%" (Decision 3) and "8% d-block retention" (Decision 3) figures
   — both measured under independent coefficient draws, before P7 existed,
   and do not transfer to P7's shared-rung regime without re-measurement.
7. **`crates/borbax-ui/tests/acceptance.rs`'s
   `no_element_reaches_the_screen_under_a_real_name_or_symbol` breaks the
   day P6/P7 land, and its replacement must not ask production what it
   thinks the identity configuration is.** It currently sweeps seeds
   `1..120` asserting *no* seed in that range shows a real name — true
   only because no seed anywhere can produce one yet. Once the identity
   configuration exists, expected identity seeds in `1..120` are roughly
   ten at `m = 10` under the corrected `P(identity) = 1/(m+1)` (item 8,
   below) — the test goes red the day P7 lands, correctly, because its
   current universal negative becomes false. **The tempting replacement —
   branch the assertion on `state.is_identity()` or an equivalent
   production predicate — is the vacuous shape this project has already
   shipped once** (`naming.rs`'s recorded 96.7%-vacuous G2 test): if a bug
   widens the identity predicate, the test agrees with the bug and stays
   green, because it is asking the thing under test what the thing under
   test believes. **Fix, three assertions, none of which asks production
   for its own opinion**: a **pinned positive** — one specific seed in
   range, found by search exactly as Decision 10 already describes
   ("three identity seeds found in 21 seconds"), hardcoded, asserted to
   show real names; a **negative** over every other seed in the same fixed
   range, asserted to show none; and an **anti-vacuity assertion on the
   *count* of identity seeds actually observed** in the fixed range — the
   one that fires if the identity predicate widens and starts matching
   seeds it should not, which the pinned-positive/negative pair alone
   would not catch (a widened predicate could still pass both if the
   pinned seed and the excluded seeds happen not to cross the new,
   incorrect boundary).
8. **`P(identity)` is a product decision, not a tuning parameter, and the
   plan must say what it decided rather than leave the number implicit in
   `m`.** State the *target* `P(identity)` explicitly and derive `m` from
   it, rather than choosing `m` and reporting whatever probability falls
   out. The search-cost argument for making identity common is weaker than
   it looks: the identity seed only has to be found **once**, and can then
   be shipped as a named universe (exactly as Step 1b already does for
   other seeds) — which is a strictly better product answer than raising
   `P(identity)` for every future search. State plainly, in words a reader
   cannot miss, that under this design **most seeds do not produce
   genuinely different chemistry from the real one — they produce a
   bounded perturbation of it**: `base * (1.0 + strength * p)` with
   `strength ∈ [0, 1]` and `|p| ≤ P_MAX` (P7, above — **round 4
   correction: not `p ∈ [-1, 1]`, which contradicts P7's own excursion
   bound and, at `p = -1, strength = 1`, reopens exactly the annihilation
   case that bound exists to close**) means no migrated constant can ever
   change sign, and the reachable parameter space is a box centred on and
   shaped by the real
   universe's own values. Under §5's revision that is the intended
   design — but a reader who inherits the pre-revision sense that
   "different seeds produce genuinely different, unrelated chemistry"
   will be wrong about this epic's output unless the document says so in
   those words. (Do not requote the specific "1 in 11" figure from this
   plan's review history without checking which rung-arithmetic reading
   it assumes — see item 6, and Decision 10's rung-arithmetic fix, above:
   the same reasoning gives materially different numbers under the two
   readings that decision closes.)

**Round 4 status, superseding the "next" paragraph below (kept for its
still-open items, which round 4 did not touch): a fourth, narrower
`/review-plan` pass targeted at exactly the scope the paragraph below
called for — Decisions 9–12, P6, P7, and Task 26.2 — has run and its
findings are applied throughout this document, each marked "round 4
correction."** Per that pass's coordinator, P1–P5 are unblocked and ready
to implement now; P6 and P7 were blocked pending fixes that are now
applied above; Task 26.2 has one design decision this pass deliberately
left open rather than resolving unilaterally — the composition-argmax
overflow boundary rule (clamp vs. drop, Decision 11's coverage-band fix,
above) and the `n_elements`-domain/`t_max`/ratio-band interaction it sits
inside — which needs a decision before Task 26.2 Step 1 is implemented.
The `the_score_is_not_a_size_comparison` test corpus (`binding.rs:1076`)
still needs extending over the reachable perturbed ratio range once P7
lands, per the w_shape/w_charge fix above — a test-suite change, not a
plan edit, and easy to lose track of since nothing in this document's own
gate gates it.

Original "next" paragraph, its own H-numbered findings from round 3's
process (distinct from round 4's H-numbered coordinator findings elsewhere
in this document — the two schemes collide by coincidence, not by
reference): **three smaller items to carry into Task 26.1's own review
rather than resolve here**, because they depend on artefacts Task 26.1
produces and this document does not yet have: the prefilter's data layout
(`SigSummary` at 672/2592 bytes, unreviewed by anyone so far);
re-measuring round 3's `w_shape`/`w_charge` ratio-to-correlation curve on
an actual V2 table once Decision 4's prefilter lands, since the measured
curve above is on V1 where `affinity` never straddles zero; and a single
stated pattern for `ShellLaw`/`Nucleus`-shaped version-tagged data,
so V2 and V3 do not each invent their own.
