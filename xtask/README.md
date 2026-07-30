# `xtask` — the checks that protect the fiction

Borbax's central promise is that **no real chemistry is in this repository**
(spec §5). That promise is easy to state and easy to erode: one helpful pull
request adding "just a small reference table" and it is gone.

So the promise is enforced by a program, and that program is `xtask`. It runs in
CI and in the pre-commit hook.

## What it checks

| Guarantee | The check |
|---|---|
| **G1** — no real chemistry data | No data files under the chemistry crates. No element tables, reaction databases, molecular structures or sequence data anywhere. |
| **G2** — generated names cannot collide with real ones | The name generator's blocklist is present and consulted. |
| **G4** — no real-world units | Temperature is `Thermal`, energy `Quanta`, distance `Span`, time `WorldYear`, mass `Mass`. |
| **G5** — permanent non-goals | No SMILES, InChI, MOL, PDB or FASTA anywhere; no Borbax↔real mapping table. |
| §13.1 | Textual checks on determinism hazards — though `clippy::disallowed_methods` is the real authority here. |

## Why a program and not a review checklist

Because reviews get tired and rules that live only in documentation are
*advisory*. This project has already found several guards that were present,
believed, and completely inert — a lint config that silently applied to nothing,
a doc-link check nobody was running, a test asserting a tautology.

The pattern is consistent enough to be a rule: **if it is not executed, it is
not a guarantee.**

## Running it

```bash
cargo xtask          # an alias for `cargo run -p xtask --`
```

CI and the git hooks use the long form deliberately: flags belonging to
`cargo run` have to sit before `run`, and an alias cannot splice them in.

## Adding a check

The bar is whether the thing being checked would otherwise be enforced by
nothing but good intentions. If a human has to remember it, it belongs here.
