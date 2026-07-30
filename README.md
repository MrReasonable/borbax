# Borbax

Chemical evolution — abiogenesis up to protocells — in a universe whose
chemistry is **invented from a 64-bit seed**, not modelled on the real world.

Molecules are shapes. Two of them react when their shapes fit together and
their surface characters oppose. Catalysis, membranes and heredity are not
implemented; they are downstream consequences of that one test.

## The premise

Real chemistry is, at the level that matters for life, about shape. Borbax
takes that literally and refuses to model anything else:

```
For each of the 60 rotations R in the icosahedral rotation group:
    for each direction i:
        j = ANTI[ PERM[R][i] ]
        shape_term  = −( r_A[i] + r_B[j] − IDEAL_GAP )²      # bumps meet hollows
        charge_term = −( a_A[i] + a_B[j] )²                  # + faces meet −
    score(R) = Σᵢ ( w_shape · shape_term + w_charge · charge_term )

affinity(A,B) = max over R of score(R)
```

That function is the entire mechanism. Two molecules sticking together, an
enzyme recognising a substrate, a membrane self-assembling, and a pore gating
what crosses it are all *the same call* at a different scale. If a feature
seems to need a second mechanism, the feature is wrong.

## Three things that make it unusual

**The chemistry is invented, and that is enforced.** There is no periodic
table in this repository. Elements, valences, bond energies and folding rules
are all generated from `universe_seed`, so a Borbax universe has its own
chemistry that nobody — including us — knew in advance. A CI gate rejects real
element names, real chemical file formats, and any data file in a chemistry
crate. Temperature is measured in *thermals*, energy in *quanta*, distance in
*spans*: distinct types, and mixing them does not compile.

**Nothing about life is hardcoded.** No `struct Cell`, no replication rule, no
"if conditions are right, spawn a protocell". Catalysis is not implemented — a
folded chain that happens to have two well-placed cavities holds two reactants
adjacent, and the rate enhancement falls out of the geometry. Protocells and
death are *detected*, via autocatalytic-set closure, never declared. A special
case for life would mean the physics is wrong.

**The same seed gives the same universe, everywhere.** Bit-identical on
macOS/aarch64, Linux/x86-64 and Windows/x86-64, enforced by a CI matrix that
diffs golden state hashes across all three. That means counter-based RNG, no
`HashMap` iteration in any result path, pinned float accumulation order,
fixed-point mass so conservation is exact by construction, and every
transcendental routed through one portable implementation — because `exp` and
`cos` genuinely differ between platform libms, and a world you cannot
reproduce is a world you cannot share.

## Status

**V0, in progress — Tasks 1–2 of 21 complete.** The workspace, toolchain and
fiction-guarantee gate are in; so are the invented-unit types and the portable
transcendental chokepoint. Universe generation is next.

| | | |
|---|---|---|
| **V0** | Chemistry in a test tube: universe generation, molecules, shape and binding, folding, decay, a beaker harness, SVG rendering. No world, no live viewer. | ← *here* |
| **V1** | One focused region, multi-scale engine, keyframes, neutral shadow, the Chronicle. **Target: protocells.** | The hard part |
| **V2** | Template-copying polymers, mutation, true heredity. | Life gets a genome |
| **V3** | Cells, ecology, a planetary field. | The world gets big |

Genomes are deliberately deferred: V1 succeeds or fails on its own criteria,
however tempting V2 becomes.

## Getting started

```bash
proto install          # toolchain, from .prototools
cargo xtask setup      # prek + git hooks, once per clone
cargo test --workspace
cargo xtask            # fiction-guarantee checks
```

`CLAUDE.md` has the full quality gate and explains why each command is on it.

## Layout

```
crates/            the simulation, in dependency order — ✅ exists, ○ planned
  ✅ borbax-units      invented units + the portable-transcendental chokepoint
  ✅ borbax-rng        counter-based deterministic streams
  ✅ borbax-universe   generated elements, bond energies, constants
  ○  borbax-molecule   graph → canonical form → 3D embedding → signature → binding → fold → cavity
  ○  borbax-reaction   reaction classes and rates
  ○  borbax-beaker     the well-mixed reactor
  ○  borbax-geometry   → borbax-render (headless SVG, golden-tested)
  ○  borbax-cli        beaker, sweep, band, battery, render, goldens
experiments/       measurement harnesses that answer a question and then stay as evidence
xtask/             the fiction-guarantee gate
docs/
  superpowers/specs/    the PRD — numbered sections (§n) are cited everywhere
  superpowers/plans/    the V0 implementation plan, Tasks 1–21
  experiments/          measured results, with their methods and near-misses
```

Crates later in that list depend on earlier ones and never the reverse. The
planned ones arrive in plan order; nothing is stubbed ahead of time.

## How the maths works, per crate

Each crate has a README written for a reader who wants to understand the
*ideas* rather than the API — what problem the maths solves, why that approach
and not the obvious one, and where the constants came from. They assume no
background and cite sources.

| README | What it explains |
|---|---|
| [`crates/borbax-rng`](crates/borbax-rng/README.md) | Why random numbers are *computed from a counter* rather than shuffled, and how Philox scrambles them |
| [`crates/borbax-units`](crates/borbax-units/README.md) | Units the compiler enforces, why mass is an integer, and why `exp` is banned |
| [`crates/borbax-universe`](crates/borbax-universe/README.md) | Elements as packed spheres — the shell law, the binding peak, and where valence and bond energies come from |
| [`experiments`](experiments/README.md) | The measurement harnesses, locality, and the neutral shadow |
| [`xtask`](xtask/README.md) | The fiction-guarantee gate and why it is a program rather than a checklist |

## Reading order

1. `docs/superpowers/specs/2026-07-26-borbax-prd.md` — the design, and why each
   decision is the way it is. §3 is the principles; §8 is the chemistry; §13 is
   the determinism contract.
2. `CLAUDE.md` — the working agreement: hard invariants, review precedence, and
   how work reaches `main`.
3. `docs/superpowers/plans/` — the task-by-task plan, with the "as built" notes
   recording where execution disagreed with the plan and who won.

## Two claims this design refuses to make

Both come from the prior-art survey, and both are places where the field's
consensus is weaker than it looks:

- **Never assume self-maintaining units compose into higher-order ones.** It is
  the most-cited claim in artificial chemistry and its weakest result. V2 onward
  treats composition as an open question rather than something that will fall
  out.
- **Never declare a plateau by eyeballing a flattening curve.** Fit competing
  models and compare, and read every novelty metric against a neutral shadow
  run — otherwise "look what evolved" is a person looking at a graph and feeling
  optimistic.

## Licence

MIT — see [LICENSE](LICENSE). `[workspace.package] license` already said so;
this is the file that makes it true.
