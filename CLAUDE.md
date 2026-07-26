# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Borbax simulates chemical evolution — abiogenesis up to protocells — in a universe
whose chemistry is **invented from a 64-bit seed**, not modelled on the real world.
Molecules are shapes; binding is a geometric complementarity test; catalysis,
membranes and heredity are all downstream consequences of that single mechanism.

**The repository currently contains no code.** It holds a spec and an
implementation plan. Read them before proposing anything:

| Document | What it is |
|---|---|
| `docs/superpowers/specs/2026-07-26-borbax-prd.md` | The PRD. Numbered sections (§n) are referenced everywhere — cite them in code comments and commit messages. V0 exit criteria are §23. |
| `docs/superpowers/plans/2026-07-26-borbax-v0.md` | V0 plan, Tasks 1–10. Holds the Global Constraints and File Structure that both plan files share. |
| `docs/superpowers/plans/2026-07-26-borbax-v0-chemistry.md` | V0 plan, Tasks 11–21 (folding onward). Complete. |

Work is plan-driven: execute tasks in order with `superpowers:subagent-driven-development`
or `superpowers:executing-plans`, ticking the `- [ ]` checkboxes. Every task is TDD —
write the failing test, verify it fails for the stated reason, then implement.
Commit after every task, conventional-commit style, with the trailer
`MrReasonable <4990954+MrReasonable@users.noreply.github.com>`.

Next action if starting fresh: **Task 1** — workspace, `.prototools`, `xtask`, CI.

## Commands

These exist once Task 1 lands; none of them work before that.

```bash
proto install                                          # toolchain from .prototools (Rust only)
cargo run -p xtask -- setup                            # prek + git hooks (once per clone)
cargo test --workspace                                 # full suite
cargo test -p borbax-molecule signature                # one crate, one filter
cargo test -p borbax-units --doc                       # doc tests (compile_fail unit-mixing tests live here)
cargo run -p xtask                                     # fiction-guarantee checks (spec §5) — CI gate
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo run -p borbax-cli --release -- goldens --emit    # golden state hashes (cross-platform matrix)
UPDATE_GOLDENS=1 cargo test -p borbax-render           # regenerate SVG goldens after a deliberate render change
```

`borbax-cli` subcommands: `beaker`, `sweep`, `band`, `battery`, `render`, `goldens`.

Transcendentals route through `crates/borbax-units/src/det_math.rs`, backed by
the `libm` crate — rust-lang's pure-Rust MUSL port, no platform dispatch. It
lives in `borbax-units` and nowhere else because `borbax-rng` and
`borbax-universe` both call it and both sit below `borbax-molecule` in the crate
order; anywhere higher and the one file guaranteeing portability is unreachable
from two crates that need it.

Every `exp`, `ln`, `sin`, `cos` and `powf` in the workspace goes through it.
`sqrt` and the four arithmetic operations stay native — IEEE-754 specifies those
exactly. A direct `.exp()` or `.cos()` on an `f64` in library code is a
determinism bug, not a style choice.

The dev profile is deliberately `opt-level = 1` (deps at 2): the unit newtypes are
zero-cost only after inlining, and the beaker harness's value is that a
chemistry-tuning iteration takes seconds, not minutes.

## Architecture

Two seeds, and a strictly one-directional pipeline:

```
universe_seed → universe_gen  (elements, valences, bond energies, folding rules)
                    ↓          generated once, immutable, cached
world_seed  →  world_gen      (V1: terrain, vents, climate, abundances)
                    ↓
              simulation core  bulk 10⁶ patches ⇄ stochastic 10³ ⇄ molecular ~10  (V1)
                    ↓          V0 has none of this — one well-mixed beaker
              borbax-geometry → SVG (headless, golden-tested) | egui+wgpu (V1)
```

Crate dependency order — later depends on earlier, never the reverse:

`borbax-units` → `borbax-rng` → `borbax-universe` → `borbax-molecule` →
`borbax-reaction` → `borbax-beaker`, with `borbax-geometry` → `borbax-render`
and `borbax-cli` on top.

`borbax-molecule` carries the conceptual weight (graph → canonical form → 3D
embedding → signature → binding → fold → cavity) and is split into small files
because that is where nearly all iteration happens.

### The three ideas that everything else is a consequence of

**1. One mechanism, reused at every scale.** `affinity(A, B)` scores two shape
signatures over the 60-element icosahedral rotation group: bumps must meet
hollows, positive faces must meet negative ones. Molecules sticking together,
an enzyme recognising a substrate, a membrane self-assembling, and a pore gating
passage are all *that same function* at a different scale. If a feature seems to
need a second mechanism, question the feature (spec §8.3).

**2. Nothing about life is hardcoded.** No `struct Cell`, no replication rule, no
"if conditions are right, spawn a protocell". Catalysis is not implemented — a
folded chain that happens to have two well-placed cavities holds two reactants
adjacent, and the rate enhancement falls out of the geometry (§8.5). Protocells
and death are *detected* (via RAF closure), never declared. A special case for
life means the physics is wrong; fix the physics.

**3. Everything expensive is per-species, not per-molecule** (§8.6). Canonicalisation,
3D embedding, signature construction, folding, cavity extraction and per-bond decay
rates are pure functions of the species, computed once at intern time.
**No code reachable from a simulation step may call any of them.** Violating this
makes the performance budget unreachable by ~2 orders of magnitude, and the
regression usually arrives disguised as a small feature — most often a
per-molecule damage field, which would turn decay from one propensity channel
back into a sweep over every molecule.

### Design choices that are load-bearing, not incidental

- **Icosahedral sampling** (D ∈ {12, 42, 162} only) makes rotation an exact
  *permutation* of sample directions — precomputed, no interpolation, no
  resampling error. Intermediate resolutions destroy that and are not options.
  42 is the prior; V0 sweeps and **locks** the value (§22.2).
- **Stress majorization, not Laplacian eigenvectors**, for the 3D embedding.
  Eigenvector sign is arbitrary and degenerate eigenvalues are common in the
  symmetric graphs small molecules actually are — that arbitrariness would
  propagate straight into species identity (§8.2).
- **FCC lattice for folding** (12 neighbours, integer coords with even parity)
  and **3D throughout**: in 2D a chain has ~√n interior against √n boundary, so
  burial protects nothing and a cavity is an open notch. The catalysis and
  solvent-attack mechanisms both depend on real burial (§8.4).
- **Binding pairs each direction with its partner's *antipode*.** Two bodies in
  contact touch along opposite directions, so the kernel indexes `b` through
  `ANTI[PERM[R][i]]`, never `PERM[R][i]` directly. Dropping `ANTI` is not a
  small error: because `−I` is not in the rotation group, it silently converts
  the search into the 60 *improper* elements — reflections only — and nothing
  fails. A complement built at the same index cannot test this.
- **Binding searches rotations but not reflections** — 60 elements, not 120. This
  makes Borbax chemistry *handed*, so homochirality is an emergent result the
  simulation could produce rather than something excluded by construction (§22.8).
- **Decay is foundational, not a balancing mechanic.** Selection *is* differential
  persistence; without decay, fitness differences can never express themselves
  (§9.4). Every decay source is a physical process already being run for other
  reasons. The decay band is the most sensitive parameter in the system — suspect
  it first whenever a run produces nothing.

## Hard invariants

Violations of these are correctness bugs, not style preferences.

**Determinism** — same seeds → bit-identical output on macOS/aarch64,
Windows/x86-64 and Linux/x86-64, enforced by a CI golden matrix (§13.4, §13.6):

- All randomness comes from `borbax_rng::Stream` (counter-based, one stream per
  `Domain` × index). Never `rand::random`, never a thread-local RNG, never
  `RandomState`. `Stream` is deliberately not `Copy`.
- No `HashMap`/`HashSet` iteration in any result-affecting path. Use `BTreeMap`,
  or sort into a `Vec` first. Every float-key sort carries an ID tie-break.
- No wall-clock. No scheduler-ordered `.sum()`/`.reduce()` — fold into per-chunk
  partials indexed by chunk ID, combine in index order.
- No platform transcendentals. `exp`, `ln`, `sin`, `cos`, `powf` are not
  correctly-rounded by IEEE-754 and genuinely differ between libms; route them
  through a vendored pure-Rust implementation or a table. `+ - * /` and `sqrt`
  are exactly specified and may stay native.
- Float accumulation order is pinned and commented as load-bearing. Nobody
  "tidies" it later.
- `Mass` is fixed-point `i64`, so conservation is exact by construction. Every
  reaction conserves mass exactly, checked in tests.
- `f64` everywhere in simulation arithmetic; no `f32`.
- The `.prototools` Rust pin is part of the determinism contract. Bumping it
  invalidates golden hashes and requires deliberate regeneration, reviewed like
  a change to physics.

**Fiction guarantees** (§5, enforced by `cargo run -p xtask`):

- No real chemistry data anywhere — no element tables, reaction databases,
  molecular structures or sequence data. No data files in `borbax-universe`,
  `borbax-molecule`, `borbax-reaction`.
- Generated element names are blocklist-checked against real symbols and names
  (`REAL_ELEMENT_SYMBOLS` in `borbax-universe/src/naming.rs`).
- No real-world units. Temperature is `Thermal`, energy `Quanta`, distance
  `Span`, time `WorldYear`, mass `Mass` — distinct newtypes; cross-unit
  arithmetic must not compile.
- Never SMILES, InChI, MOL, PDB or FASTA import/export, no Borbax↔real mapping
  table, no real-world calibration targets. Feature requests implying otherwise
  are rejected on sight.

**Other:** no `unwrap()`/`expect()` in library code (`clippy::unwrap_used` and
`expect_used` are `deny` at workspace level); tests may unwrap freely.

## Scope discipline

V0 is chemistry in a test tube: no world, no spatial simulation, no live viewer —
but SVG rendering from the start, because a shape-based chemistry cannot be
developed blind. V1 adds the region, multi-scale engine, keyframes, neutral
shadow and the Chronicle, targeting protocells. Genomes are **V2** and are
explicitly deferred however tempting they become; V1 succeeds or fails on §15.1
alone.

Two claims the design deliberately refuses to make, both from the prior-art
survey (§2.7): never assume self-maintaining units compose into higher-order
ones (the field's weakest result), and never declare a plateau by eyeballing a
flattening curve — fit competing models and compare, and read every novelty
metric against a neutral shadow run.

## Agents

Six project agents in `.claude/agents/`, most mapped to a named risk in §20
rather than to a generic role. They review; none of them writes code.

| Agent | Invoke when |
|---|---|
| `rust-developer-expert` | On any new module before it settles; when a type or trait is introduced; when code works but reads badly; as the last pass before a task is committed. Also owns **dependency review** — whether a well-established crate should be doing this instead of us. |
| `rust-performance-expert` | Before committing to a data structure in a hot path; when a task specifies performance targets; to audit that an optimisation preserved determinism. |
| `determinism-auditor` | Before adding parallelism or a dependency that touches results; when a golden hash moves; on any diff touching float arithmetic, collection iteration, sorting or RNG. Every hazard here is invisible on one machine. |
| `emergence-auditor` | Before adding a threshold, special case or per-molecule field; when a behaviour is not emerging and there is a temptation to help it along; on any chemistry or simulation diff. |
| `geometry-numerics-reviewer` | On canonicalisation, geodesic construction, the rotation table, stress majorization, signatures, binding, FCC folding, cavity detection — code that is either correct or silently wrong forever. |
| `alife-researcher` | When implementing an algorithm or metric from the literature (RAF, Gillespie, Bedau activity statistics, neutral networks, plateau fitting), or when a design decision rests on a cited result. Constrained never to return real chemistry (§5). |

The auditors are independent lenses on the same diff and can run in parallel.
Run the emergence and determinism auditors before any commit that touches
physics — those are the two invariants no test suite fully protects.

### Review precedence

The agents are told to argue with each other rather than defer, so genuine
conflicts will surface. This order breaks ties. **It ranks how expensive a
mistake is to discover late, not how important each concern is** — every one
of them matters.

1. **Fiction guarantees (§5, G1–G6).** Absolute. No result is worth breaching
   them, and there is nothing to trade against.
2. **Correctness of the physics.** Geometry and numerics that are silently
   wrong poison everything downstream and are the hardest thing here to
   detect.
3. **Determinism (§13.1, §13.4).** A result nobody can reproduce cannot be
   verified, shared, or debugged.
4. **Emergence invariants (§3, §8.6).** A special case for life means the
   physics is wrong.
5. **Performance (§17).** The budget is real, but a fast wrong answer is
   worthless.
6. **Idiom and readability.** Genuinely matters, and yields to the above.

Two notes on applying it. **2 above 3 is deliberate**: when a numerically
better formulation changes bit patterns, adopt it and regenerate goldens as a
deliberate physics change — do not preserve a worse formulation to protect a
hash. And **most conflicts are not about this order at all** — they are two
correct observations pulling opposite ways, where the resolution is a
measurement or a line of the spec, not a ruling.

## Dependencies

The best code is code someone else maintains. Reach for the settled crate:
`thiserror` for library errors, `anyhow` for binaries only, `clap` for CLI,
`proptest` for property tests, `insta` for snapshot and golden tests,
`criterion` for benchmarks, `libm` for portable transcendentals, `serde` for
serialisation, `rustc-hash` for deterministic hashing.

The asymmetry that governs the decision:

- **Dev-dependencies are cheap.** `proptest`, `insta`, `criterion`, `rstest`
  cannot affect simulation output. Low bar.
- **A runtime dependency in a result-affecting path is expensive.** It must be
  version-pinned, and upgrading it is a physics change requiring deliberate
  golden regeneration — exactly like a compiler bump (§18.1). Run the
  determinism auditor on any such addition.

Hand-rolled, and the reason matters: three-component vector maths is ~20 lines
using none of what `glam` or `nalgebra` sell. **Not** for determinism — both
would be bit-identical (`glam`'s f64 module has no arch dispatch and documents
cross-platform bit-identity; `nalgebra` special-cases dim 3 to `a + b + c`).
Claiming otherwise was a doctrine error worth naming, because a false
determinism argument gets used later to reject a dependency that is correct.

SMACOF stays hand-rolled on firmer ground: its fixed 240-iteration count is a
determinism requirement, and no solver crate offers "stop after exactly n
iterations regardless of convergence".

*Not* hand-rolled: transcendentals (`libm`), plateau model fitting
(`levenberg-marquardt`), snapshot tests (`insta` — do not hand-roll
`UPDATE_GOLDENS`), benchmarks (`criterion` — a timing assertion inside
`#[test]` is a flaky test that will get deleted, taking the guarantee with it).

Held in reserve: `nauty-Traces-sys` for canonicalisation if the hand-rolled
search ever becomes a bottleneck. It would be a result-affecting runtime
dependency *and* C FFI against a `forbid(unsafe_code)` workspace, so it stays
an escape hatch rather than a plan.
