# Determinism audit — V0 plan, post-fix re-review (2026-07-26)

Repo contains no code, so neither cheap reproduction was runnable
(`RAYON_NUM_THREADS` sweep, double `goldens --emit`). One claim was verified by
downloading and reading crate source; everything else is reasoned from the plan
text and is marked as such.

---

## D1 · Task 10 Step 3 still creates a second `det_math` that calls the platform libm — SEVERE

**Where:** `plans/2026-07-26-borbax-v0.md:4323-4338` and its commit message
`:4367-4370`; call sites `plans/2026-07-26-borbax-v0-chemistry.md:1898, 2046,
2067`; file list `chemistry.md:3095`.

**What:** The fix pass moved `det_math` to `borbax-units` (v0.md:847-899) but did
not delete the *old* instruction. Task 10 Step 3 still says, verbatim:

> create `crates/borbax-molecule/src/det_math.rs` with a single
> `pub fn exp(x: f64) -> f64 { x.exp() }`

…with a header comment saying "**Currently delegating to the platform**… Task 25
replaces these bodies". Task 25 does not exist; V0 has 21 tasks. The code block
40 lines above it already calls `borbax_units::det_math::exp` (`v0.md:4319`),
which is decisive evidence the code was edited and the prose below it was not.

Tasks 14 and 15 then consume that shadow module: the Arrhenius factor
(`chemistry.md:1898`) and **both** thermal-decay propensities
(`chemistry.md:2046, 2067`) call `borbax_molecule::det_math::exp`.

**Does it reach a result?** Yes — every reaction rate and every decay propensity.
Those feed the Gillespie propensity tree, which selects which reaction fires.
One ULP flips a channel selection and the trajectory diverges completely.

**How it manifests:** macOS/aarch64 vs Linux/x86-64 vs Windows. Silently — and
`xtask`'s own check cannot see it, because `check_no_platform_transcendentals`
skips **any file named `det_math.rs`** (`v0.md:948-951`), so the shadow module is
exempt from the one gate designed to catch this.

**Fix:**
- Delete `v0.md:4323-4338` down to the closing fence; keep only the
  `Signature::lex_cmp_pub` sentence.
- Delete the `det_math.rs` paragraph from the Task 10 commit message
  (`v0.md:4367-4370`).
- Replace `borbax_molecule::det_math::` with `borbax_units::det_math::` at
  `chemistry.md:1898, 2046, 2067`.
- `chemistry.md:3095` — change `crates/borbax-molecule/src/det_math.rs` to
  `crates/borbax-units/src/det_math.rs`.
- Tighten the xtask exemption so it cannot re-shield a future copy:

```rust
// Exempt by *path*, not by filename. A second file called det_math.rs
// anywhere else is exactly the bug this check exists to catch.
if entry.ends_with("borbax-units/src/det_math.rs") { continue; }
```

**Confidence:** High — read directly, contradicted by adjacent code in the same
file.

---

## D2 · `libm` is a caret range and `Cargo.lock` is never committed — SEVERE

**Where:** `v0.md:181-182` (`[workspace.dependencies.libm] version = "0.2"`),
`v0.md:511` (`git add .prototools Cargo.toml rustfmt.toml .github xtask`).

**What:** `"0.2"` is `^0.2`, so any `0.2.x` satisfies it. `chemistry.md:3189`
says "pin it exactly and treat any upgrade as a physics change" — the manifest
does not do that. Task 1's `git add` is an explicit file list and `Cargo.lock` is
not in it, so nothing pins the resolved version either. Every transcendental in
the project comes from this crate.

**Does it reach a result?** Yes — this *is* the transcendental implementation.

**How it manifests:** Two clones of the same commit at different times resolve
different `libm` patch versions and produce different goldens. Also
machine-to-machine, if one dev ran `cargo update`. Diagnosing it means noticing
that a file nobody edited changed.

**Fix:**

```toml
[workspace.dependencies.libm]
# Exact, not caret. This crate *is* the transcendental implementation
# (§13.1); a patch bump is a physics change requiring golden regeneration,
# exactly like a rustc bump (§18.1).
version = "=0.2.16"
```

…and add `Cargo.lock` to Task 1 Step 7's `git add`, with the same reason in the
commit message. `.gitignore` does not exclude it, so this is a one-word change.

**Confidence:** High — read the manifest and the `git add` line; `cargo search`
confirms 0.2.16 is current, so the range is live, not vestigial.

---

## D3 · The CI emit step will report a false determinism failure on Windows — HIGH

**Where:** `v0.md:452-454`.

```yaml
- name: Emit golden hashes
  run: cargo run -p borbax-cli --release -- goldens --emit > goldens-${{ matrix.os }}.txt
```

**What:** `run:` on `windows-latest` defaults to `pwsh`. PowerShell's `>` is
`Out-File` semantics: it splits the native command's stdout into strings and
rewrites them with `[Environment]::NewLine`, i.e. **CRLF**. The macOS and Linux
files are LF. `diff -u` in `determinism-matrix` then fails on every line of an
otherwise identical file.

**Does it reach a result?** It reaches the *verdict on* results, which is worse:
the first thing the matrix ever does is fire falsely. The predictable response is
`diff --strip-trailing-cr` or `|| true`, and the gate is dead from that day.

**How it manifests:** Windows vs the other two, on the first commit after Task 21.

**Fix — force the shell:**

```yaml
- name: Emit golden hashes
  if: hashFiles('crates/borbax-cli/Cargo.toml') != ''
  shell: bash
  run: cargo run -p borbax-cli --release -- goldens --emit > goldens-${{ matrix.os }}.txt
```

Better: give the subcommand `--out <path>` so it writes the bytes itself and no
shell touches the stream. Add `newline_style = "Unix"` equivalent to the emitter
— the CLI must write `\n`, never `writeln!` through a platform newline.

**Confidence:** High on the mechanism (documented PowerShell redirection
behaviour), not run.

---

## D4 · Task 21 Step 6 is stale and contradicts Task 2 — HIGH

**Where:** `chemistry.md:3162-3198`.

**What:** "Replace the `det_math` bodies with the `libm` crate… Goldens generated
on the M1 will not reproduce on CI runners until this lands." Task 2 Step 6
already creates `det_math` backed by `libm` (`v0.md:851-899`). The stale block
also omits `powf`, which Task 5 needs (`v0.md:1700`).

**Does it reach a result?** Indirectly and severely — it tells the implementer
that the platform libm is acceptable for Tasks 2-20 and that goldens are
knowingly non-portable until the last task. Combined with D1 that is a coherent
(and wrong) instruction set which survived the fix.

**Fix:** Replace the whole of Step 6 with:

> **Step 6: Enable the cross-platform matrix.** `det_math` has been `libm`-backed
> since Task 2 Step 6; nothing to replace. Confirm `cargo run -p xtask` reports
> zero direct transcendental calls, then enable `determinism-matrix` and confirm
> goldens are byte-identical across macOS/aarch64, Windows/x86-64,
> Linux/x86-64. If they are not, this is D1's shadow `det_math` or an unrouted
> call — check both before touching a hash.

**Confidence:** High.

---

## D5 · `xtask`'s banned-call list is wrong in both directions — HIGH

**Where:** `v0.md:939-940`.

```rust
const BANNED_CALLS: &[&str] =
    &[".exp()", ".ln()", ".log(", ".sin()", ".cos()", ".tan()", ".powf(", ".powi("];
```

**False positives.** `.powi(` is banned but the plan uses it in eight
result-affecting places: `v0.md:1722, 3300, 3332, 3884` and
`chemistry.md:2229-2231`. CI goes red at Task 4 and stays red. The tempting fix
is to weaken the check — the same failure mode as S10.

`powi` *is* fine: it lowers to LLVM's target-independent `ExpandPowI` multiply
tree, which is exact repeated multiplication under a pinned rustc. The det_math
doc comment (`v0.md:892-894`) argues otherwise; if that argument is accepted the
plan must stop using `powi`, not ban it while using it. I'd keep `powi` and drop
it from the list.

**False negatives.** `.log(` does not match `.log10(` or `.log2(`. `.exp()` does
not match `.exp2()`. Entirely absent: `.atan2(`, `.cbrt(`, `.hypot(`,
`.exp_m1(`, `.ln_1p(`, `.asin(`, `.acos(`, `.atan(`, `.sinh(`, `.cosh(`,
`.tanh(`, `.to_degrees(`, `.to_radians(`, `.mul_add(`, and every free-function
form (`f64::exp(x)`, a bare `use libm::*`).

**Fix:**

```rust
/// §13.1. `powi` is absent deliberately: it is repeated multiplication,
/// exact and target-independent. `sqrt` likewise — IEEE-754 specifies it.
const BANNED_CALLS: &[&str] = &[
    ".exp()", ".exp2()", ".exp_m1()", ".ln()", ".ln_1p()",
    ".log(", ".log2()", ".log10()", ".powf(",
    ".sin()", ".cos()", ".tan()", ".asin()", ".acos()", ".atan()", ".atan2(",
    ".sinh()", ".cosh()", ".tanh()", ".cbrt()", ".hypot(",
    ".to_degrees()", ".to_radians()", ".mul_add(",
    "f64::exp", "f64::ln", "f64::powf", "f64::sin", "f64::cos",
    "libm::",  // only det_math.rs may name it, and that file is exempt
];
```

Also add: `.mul_add` is not a determinism bug on its own (FMA is exactly
specified) but it *is* a silent accumulation-shape change, so banning it outside
a reviewed site is right.

**Confidence:** High — grepped both plan files for every listed pattern.

---

## D6 · The fold cache is keyed on `u.seed`, but the plan mutates `u.consts` in place — MEDIUM-HIGH

**Where:** `chemistry.md:922-933` (`key_of`), and the mutation sites
`chemistry.md:2886, 2893, 2900` (`u.consts.rate_prefactor = 1e-30`,
`u.consts.decay_scale = 0.0`) plus `chemistry.md:3120` (`borbax band --sweep
decay_scale` bisects a constant on a fixed universe).

**What:** S11 was fixed by folding `u.seed` into the key. But `Universe` is
demonstrably mutated *after* generation, so `u.seed` no longer identifies the
universe's constants. The fold currently depends on `u.consts.w_charge` and on
element affinities (`chemistry.md:439, 608`); nothing today sweeps `w_charge`,
so this is latent rather than live — but it is one line of a future sweep away
from re-becoming S11, and the failure mode is identical: silent, order-dependent,
survives eviction.

Separately, §13.1's reproducibility tuple is `(universe_seed, world_seed,
config_hash)` and **no `config_hash` appears anywhere in either plan file**.

**Does it reach a result?** Today: no (nothing mutates a fold-relevant constant).
If `band --sweep` is ever pointed at `w_charge`: yes, and the wrong conformation
becomes the wrong species becomes the wrong chemistry.

**What would settle it:** whether `UniverseConsts` gains a sweep over any
constant read by `fold`. I could not determine that from the plan; `band --sweep
decay_scale` is the only sweep named.

**Fix — make it structural rather than vigilant:**

```rust
/// Hash of every constant the simulation reads, recomputed whenever
/// `consts` is touched. §13.1's reproducibility tuple is
/// (universe_seed, world_seed, config_hash) — the seed alone stopped
/// identifying the universe the moment the battery started mutating
/// `consts` in place (Task 20).
pub fn config_hash(&self) -> u64 { /* mix over consts, in field order */ }
```

…then key on `(config_hash, sequence)` in both `key_of` and `fold_seed`, and emit
`config_hash` in the golden header so a mismatched config fails loudly instead of
silently.

**Confidence:** Medium-high on the latency; high that `config_hash` is missing.

---

## D7 · `energy_scale` in the Metropolis criterion is undefined — MEDIUM

**Where:** `chemistry.md:519`.

```rust
|| roll < borbax_units::det_math::exp(-delta / (temp * energy_scale).max(1e-9));
```

`energy_scale` is declared nowhere in either plan file (grepped; only occurrence
is this line). It is the sole knob setting the Metropolis accept probability, so
the implementer's invented value determines every fold, hence every polymer
species identity, hence every golden in Tasks 12-21.

**Does it reach a result?** Everything downstream of folding.

**How it manifests:** Not a platform divergence — a value pinned into a golden by
accident, unreviewed, and unrecoverable later without regenerating everything.

**Fix:** derive it from the universe so it is not a bare literal (M4's lesson):

```rust
// Metropolis energy scale, from the universe rather than a literal:
// contact energies are O(w_charge), so this sets ΔE ≈ 1 at temp = 1.
let energy_scale = u.consts.w_charge.abs().max(1e-9);
```

Whatever value is chosen, name it and comment it as pinned.

**Confidence:** High that it is undefined; the suggested derivation is a
proposal, not a finding.

---

## D8 · The beaker's two propensity-update paths must be proven bit-identical — MEDIUM

**Where:** `chemistry.md:2727-2731`.

> Below `FLAT_SCAN_LIMIT` (start at 256) rescan all channels… Above it, use the
> incidence structure from Task 17 to touch only affected channels.

**What:** Two code paths compute the same tree state, selected by a runtime
channel count. The segment tree makes the *combination* order-safe (every `set`
recomputes ancestors as `left + right`, so the tree is bit-identical to a
rebuild — this part is genuinely sound and is the plan's best determinism
decision). What is *not* guaranteed is that the incidence set is complete: any
channel whose propensity depends on something the incidence graph does not model
gets a stale leaf on the fast path and a fresh one on the slow path.

**Does it reach a result?** Yes. And it is nastier than a platform divergence: a
run crosses 256 channels mid-trajectory and silently switches physics. Two runs
of the same seed at different `FLAT_SCAN_LIMIT` values disagree; the CI matrix
cannot see it because all three platforms take the same path.

**Fix — make the equivalence a test, not a hope:**

```rust
// The fast path is only legitimate if it lands on the identical tree.
// Compared by bits: "close enough" is how a selection boundary moves.
#[cfg(debug_assertions)]
fn assert_paths_agree(&self) {
    let mut rebuilt = SegmentTree::with_capacity(self.tree.cap);
    for (slot, ch) in self.channels.iter().enumerate() {
        rebuilt.set(slot, self.propensity_of(ch));
    }
    for i in 0..self.tree.cap {
        debug_assert_eq!(
            self.tree.leaf(i).to_bits(), rebuilt.leaf(i).to_bits(),
            "incidence-updated leaf {i} diverged from a full rescan"
        );
    }
}
```

…called every step under `debug_assertions`, plus one release test that runs a
seed to 100k events with `FLAT_SCAN_LIMIT` forced to 0 and to `usize::MAX` and
asserts the two final state hashes are equal.

**Confidence:** High on the hazard, reasoned from the source text.

---

## D9 · The golden state-hash digest is still unnamed — MEDIUM

**Where:** `v0.md:494-506`.

The specification improved and is now most of the way there — ordering by
`SpeciesId`, `f64::to_bits`, `Mass::raw()`, "never `DefaultHasher`". Three gaps
remain, and the third is the one that bites:

1. **No digest is named.** "a fixed digest" is still hand-waving. The only
   hashing crate CLAUDE.md names is `rustc-hash`, whose algorithm **changed
   between 1.x and 2.0** — an implementer following the doctrine picks the one
   crate that reintroduces the problem.
2. **No newline/encoding rule for the emitted text**, which is what makes D3
   bite.
3. **No statement of what happens to absent species** (skipped, or emitted with
   count 0) or to a NaN count (`to_bits` of NaN is payload-dependent, and a NaN
   count should abort, not hash).

**Fix — use what already exists, no new dependency:**

```
The digest is `borbax_rng::mix` in Merkle-Damgård form, seeded from
`Domain::Hash`. It is already in the workspace, already pure, already
version-stable by virtue of being our own source. Do not reach for
rustc-hash: FxHash's algorithm changed between 1.x and 2.0, so it is
not stable across a dependency bump.

Emit exactly:
  <config_hash:016x> <universe_seed:016x> <step> <digest:016x>\n
one line per checkpoint, LF only, ASCII only, written as bytes.
Species with count 0 are omitted. A non-finite count aborts the run
with an error rather than being hashed.
```

**Confidence:** High.

---

## D10 · `depth_sort` is still declared and never specified — LOW-MEDIUM

**Where:** `chemistry.md:2749` (interface list). No implementation, no
requirement text, and Task 19 Step 3's bullet list does not mention it.

The previous audit's Low finding ("needs an id tie-break before it exists") was
not applied. To be precise about severity, which the earlier note was not: a
float-keyed `sort_unstable_by` is still *reproducible* here — identical input
plus identical rustc gives identical pdqsort output, and the depths come from the
bit-identical embedding. So this is **not** a cross-platform hazard. It is a
fragility hazard: equal depths reorder under any change to input permutation,
which silently churns an `insta` SVG golden and trains the team to accept
snapshot diffs unread.

**Fix — one line in Step 3's bullet list:**

```rust
// Depth ties broken by atom index. Equal depths are common (symmetric
// molecules, planar rings), and without the tie-break they reorder
// under any change to input order, churning the SVG golden for no
// visible reason.
order.sort_by(|&a, &b| depth[a].total_cmp(&depth[b]).then(a.cmp(&b)));
```

**Confidence:** High that it is unspecified.

---

## D11 · Task 20b: the plateau fitter's reproducibility claim is half-true — LOW-MEDIUM

**Where:** `chemistry.md:3056-3078`.

Three separate points, and the headline claim survives:

**(a) "does not enter golden hashes" — defensible, and I'd keep it.** The fitter
runs over already-emitted trajectories and nothing feeds back. That is a correct
reading and I am not blocking the dependency.

**(b) The `ln` in AIC/BIC is unrouted.** Line 3075 says route the models'
`exp`/`powf` through `det_math`; line 3069's `aic = n·ln(rss/n) + 2k` and
`bic = n·ln(rss/n) + k·ln(n)` are not covered by that sentence. AICc and BIC
routinely differ by <1 unit, which is exactly the margin where a last-bit
difference flips the verdict. Add `ln` to the sentence.

**(c) LM terminates on a tolerance, which is the hazard SMACOF's fixed 240
iterations exists to avoid.** `levenberg-marquardt`'s `LevenbergMarquardt`
defaults stop on gradient/step tolerances, so the iterate count — and therefore
the fitted parameters — depends on the last bits of the residual. Same input,
different platform, different iteration count, potentially different model
selected.

**Fix:** set `.with_patience(n)` / an explicit max-iteration cap and treat the
result as "the parameters after exactly N iterations", the same contract as
SMACOF. Document it in the module header alongside the `det_math` routing.

**(d) Crate placement.** `levenberg-marquardt` depends on `nalgebra`. Putting it
in `borbax-beaker`'s `[dependencies]` makes a large numerics tree reachable from
simulation code, where the next person to need a matrix will reasonably use it.
Spec §18.3 already has `borbax-metrics` for exactly this. Move
`{activity,shadow,plateau}.rs` there, or make it CLI-only.

**Confidence:** (b) and (d) high; (c) high on the mechanism, from the crate's
documented default behaviour rather than from reading its source this session.

---

## D12 · `Domain::Shadow` is declared and never used — LOW

`Domain::Shadow = 8` exists at `v0.md:1138`. Task 20b Step 3 ("fork from a
keyframe with decay rates equalised") states no RNG discipline at all — the
obvious implementation clones the focal run's stream, which is precisely the
duplicate-replay that dropping `Copy` exists to prevent, and the shadow would
then be a *correlated* control rather than an independent one.

**Fix — state it in Step 3:**

```rust
// The shadow draws from Domain::Shadow, never a clone of the focal
// run's stream. A cloned stream replays the identical sequence, which
// makes the "control" perfectly correlated with the thing it controls —
// and the comparison then measures nothing. Stream is not Copy so this
// has to be written explicitly either way (§13.1).
let shadow_rng = Stream::new(seed, Domain::Shadow, fork_step);
```

---

## Smaller items, grouped

- **`fold_seed` does not share construction with `key_of`, contrary to its own
  comment.** `chemistry.md:377-391` vs `922-933`: different seed constants
  (`0xF0_1D_5EED` vs `0xB0_1BAA_5EED`), different indices (3 vs 1), and
  `fold_seed` is 64-bit while `key_of` is 128-bit. Both are pure functions of
  `(sequence, universe)` so cached and uncached folds *do* agree — the invariant
  holds, the stated reason for it does not. Either factor out a shared
  `mix_sequence(p, u) -> u128` and take `fold_seed` as its low 64 bits, or fix
  the comment. Right now a future edit to one will not obviously require the
  other.
- **`fold_seed` is called unqualified from an inherent method** (`chemistry.md:418`,
  `self.fold_with_seed(p, u, fold_seed(p, u))`) while being defined inside
  `impl FoldWorkspace` at :383. Needs `Self::fold_seed`. Compile error, not a
  determinism bug.
- **`fold_with_seed` is reachable from production.** It is `pub`, so
  "diagnostics only" is a comment, not a constraint. Make it
  `#[cfg(any(test, feature = "diagnostics"))]`, or `pub(crate)` plus a test-only
  re-export. The comment is doing work the type system should do.
- **The SMACOF accumulation is not commented as load-bearing** (`v0.md:3533-3550`).
  `num[k] += w * (pos[j][k] + dij * diff[k] / dist)` and `den += w` accumulate in
  `j` order and feed species identity. The centroid three lines later *is*
  commented ("Fixed summation order (spec §13.1)"). This is the more important of
  the two. Same for `may_bind`'s `a.r.iter().sum()` (`v0.md:4299-4300`), which
  gates a discrete accept/reject, and `contact_energy` (`chemistry.md:614-620`).
- **`Fold::energy` (`chemistry.md:342`) and `contact_energy` (`:614`) are still two
  copies of the same loop.** Deterministic today; the hazard is that reordering
  one leaves the annealer's `current` disagreeing with the reported energy.
- **`CacheStats.bytes` uses `size_of::<Fold>()`** (`chemistry.md:969`). `repr(Rust)`
  layout is not guaranteed stable across compilations. Eviction is result-neutral
  so this is fine — but `CacheStats` must never be emitted into a golden state
  hash. Worth one line in Task 12.
- **`radiogenic_rate` is both a `SpeciesRecord` field (`chemistry.md:1651`) and a
  free function recomputing it (`:2074-2081`)**, and the `intern` constructor at
  `:1772-1779` populates neither it, `kind`, `fold`, nor `cavities`. Two sources
  of the same number is how they drift.
- **`xtask`'s `read_dir`** (`v0.md:378`) iterates in filesystem order. The set of
  failures is order-independent, so this is diagnostic-only — but sorting makes
  CI log diffs readable and costs one line.
- **The cavity outside-flood does not check `fcc::in_bounds`** (`chemistry.md:1235`)
  where step 2 does (`:1244`). `in_box` decodes via `from_flat`, which wraps. Same
  class as S5, deterministic either way, so this is a note for the geometry
  reviewer rather than a determinism finding.
- **`catalysis_factor` still contains the bare `1.0e3`** (`chemistry.md:2249`)
  despite the header at `:2188-2190` saying it moved to
  `u.consts.catalytic_prefactor`. Emergence's finding, not mine, but the fix pass
  clearly edited the comment and not the code — the same pattern as D1.

---

## Checked and found sound

Verified by reading, not assumed:

- **`libm`'s `arch` default feature is safe for the five functions used.**
  Downloaded `libm 0.2.16` and read it. `arch` overrides only `sqrt`, `fma`,
  `rint`, `ceil`, `floor` (`src/math/arch/{aarch64,x86,i586}.rs`) — all exactly
  specified by IEEE-754. `exp.rs`, `log.rs`, `pow.rs`, `sin.rs`, `cos.rs` contain
  no `fma` call and no arch dispatch; only `cbrt.rs` uses `fma`. The runtime
  CPU-feature detection in `arch/x86/detect.rs` is reachable only through `fma`
  and `cbrt`, neither of which the project calls. `default-features = false` is
  therefore not required — but add `cbrt` to the banned list (D5) so it stays
  true.
- **`Stream` is `#[derive(Debug, Clone)]`, not `Copy`** (`v0.md:1163`), and no
  by-value duplication survives in either plan file. The M9 shell-pattern
  duplicate is gone; the only `fork` uses are in Task 3's own test at
  `v0.md:1076-1079`, which forks a `base` by reference and is correct.
- **`fork`'s new odd constant does not move any pinned golden.** The only pinned
  RNG values are `GOLDEN_UNIVERSE_0` (`v0.md:1268-1273`), which Task 3 Step 5
  generates *after* implementation from the failure message. Nothing else in
  either file states a numeric RNG expectation.
- **`Domain::Hash = 9` is used consistently** and is genuinely decoupled: both
  `key_of` (`chemistry.md:923-924`) and `fold_seed` (`:384`) use it, no live
  domain is used as a hash function anywhere, and adding a draw in folding cannot
  now shift a cache key.
- **Exactly one `HashMap`, still never iterated.** `FoldCache.map`
  (`chemistry.md:899`) sees only `get`/`insert`/`remove`/`is_empty`; eviction
  walks `slots` by index. `PassThrough` replaces `RandomState`. `Interner` uses
  `BTreeMap`; Task 17 uses `BTreeSet`/`VecDeque`; Task 13 uses
  `BTreeMap`/`BTreeSet`/`VecDeque`. No `HashSet`, no `DefaultHasher`, no
  `RandomState` outside a ban comment.
- **No rayon, no `par_iter`, no threads, no `Instant`/`SystemTime`, no
  `env::var`, no `rand::`, no `f32`** anywhere in either file. Grepped
  exhaustively. Parallel-reduction hazard is currently zero because there is no
  parallelism — worth restating when Task 20's battery gets parallelised over 40
  universes, which is the obvious next step and the obvious next hazard.
- **Task 13's outside-flood is deterministic.** `boundary_cells` derives from
  `box_cells`, which iterates `z`, then `y`, then `x` in fixed ranges;
  `BTreeSet`/`VecDeque` preserve that; `fcc::NEIGHBOURS` is a fixed array.
  Component discovery iterates `candidates` (a `BTreeSet`) in ascending order.
  `cells.sort_unstable()` on distinct `i32` has no ties. `enclosure` is integer
  division. Clean.
- **Task 13's `lining_signature` summation order is pinned and commented.**
  `lining` is a `BTreeSet<usize>`, both passes iterate it in the same ascending
  order, and `:1316-1317` says so explicitly. `reach_of` is pure `+`/`*`, so the
  two calls are bit-identical. The `nearest` selection uses strict `<`, keeping
  the lowest monomer index on a tie. This one is right.
- **Task 8's `canonicalise_frame` uses only `+ - * / sqrt`.** Every branch
  (`norm(pos[i]) > 1e-9 * scale_by`, `m > 1e-9 * scale_by`) tests a bit-identical
  value, so the same atoms are selected on every platform. `e2 = e0 × e1` fixes a
  right-handed frame, preserving the handedness §22.8 depends on. No accumulation
  needs pinning beyond what SMACOF already produces. Determinism clean — its
  *stability* under near-degeneracy is a question for the geometry reviewer, not
  me.
- **Task 11's unconditional `roll` claim is correct.** `let roll = rng.next_f64()`
  precedes both float comparisons, so stream position never depends on `delta` or
  on the `exp` result. Worth narrowing the comment: position *does* depend on the
  earlier `propose`/`in_bounds`/`is_free` short-circuits, but every one of those
  is an integer test, so it stays platform-independent. The claim as stated
  ("depends only on the step count") is slightly stronger than the truth and
  should say "never on a float comparison" — which is what actually matters.
- **The anneal cannot drift.** `contacts` are `u32` and `contact_energy` is
  recomputed from them each step in fixed index order, so the incremental path is
  bit-identical to a full recount rather than merely close.
- **Task 14's `solvent_rate` is clean.** `affinity_ordered` is `+`/`*`;
  `bind_probability` routes `exp` through `borbax_units::det_math`;
  `exposed_fraction` is an integer sum divided by an exact constant. `AffinityMemo`
  memoises a pure function, which cannot change a result.
- **`Signature::canonicalise`** uses strict `<` over `lex_cmp`, so ties keep the
  lowest rotation index; `total_cmp` throughout, no `partial_cmp().unwrap()`
  anywhere in either file.
- **`SegmentTree` is correct and remains the best-defended decision in the plan.**
  `set` recomputes ancestors as `left + right`, so stored totals are bit-identical
  to a rebuild regardless of update count; freeing sets 0.0 rather than compacting;
  the zero-propensity fallback scans forward deterministically.
- **`.prototools` pins `rust = "1.97.1"` with the §18.1 reasoning inline**, and
  `[profile.release]` sets `codegen-units = 1`. The dev/release opt-level split is
  not a determinism concern: LLVM does not reassociate f64 or contract to FMA
  without fast-math, so `cargo test` (opt-level 1) and `cargo run --release`
  produce identical values.
- **`hashFiles('crates/borbax-cli/Cargo.toml') != ''`** is the correct guard idiom
  and will start emitting the moment Task 21 lands.
  `actions/download-artifact@v4` with no `name` downloads every artifact into
  `<artifact-name>/`, so the `goldens-*/goldens-*.txt` glob resolves. The job
  works — modulo D3.

## Not verifiable this session

Neither cheap reproduction could be run: the repository is plan-only. Both should
be run at Task 3 (first RNG) and again at Task 18 (first trajectory), not
deferred to Task 21:

```bash
RAYON_NUM_THREADS=1 cargo test --workspace
RAYON_NUM_THREADS=8 cargo test --workspace
cargo run -p borbax-cli --release -- goldens --emit > /tmp/a.txt
cargo run -p borbax-cli --release -- goldens --emit > /tmp/b.txt
diff /tmp/a.txt /tmp/b.txt
```
