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
cargo test --workspace                                 # full suite
cargo test -p borbax-molecule signature                # one crate, one filter
cargo test -p borbax-units --doc                       # doc tests (compile_fail unit-mixing tests live here)
cargo run -p xtask                                     # fiction-guarantee checks (spec §5) — CI gate
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo run -p borbax-cli --release -- goldens --emit    # golden state hashes (cross-platform matrix)
```

`borbax-cli` subcommands: `beaker`, `sweep`, `battery`, `render`, `goldens`.

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

`.claude/agents/rust-performance-expert.md` — invoke before committing to a data
structure in a hot path, when a task specifies performance targets, and to audit
that determinism constraints survived an optimisation.
