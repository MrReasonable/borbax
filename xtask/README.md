# `xtask` — the checks that protect the fiction

Borbax's central promise is that **no real chemistry is in this repository**
(spec §5). That promise is easy to state and easy to erode: one helpful pull
request adding "just a small reference table" and it is gone.

So the promise is enforced by a program, and that program is `xtask`. It runs in
CI and in the pre-commit hook.

## What it checks

**This table lists the checks `check_guarantees` actually calls, and nothing
else.** An earlier version described a gate about four times wider than the code:
it claimed G5 was enforced "anywhere" when the scan read one file, listed G4 as a
check when nothing tests it, and omitted four checks that do run. A review
planted `SMILES`, `InChI`, `FASTA` and `MOL` in `bonds.rs` and the gate reported
all checks passed. (The scan is now genuinely repository-wide — that was fixed
rather than documented around.)

| Function | Guarantee | What it does |
|---|---|---|
| `check_no_data_files` | G1 | No data files under the chemistry crates — no element tables, reaction databases, structures or sequence data. |
| `check_blocklist_present` | G2 | `naming.rs` exists and contains `REAL_ELEMENT_SYMBOLS`. Fails loudly if the file is missing, rather than reporting green for a path that moved. |
| `check_no_real_chemical_formats` | G5 | No `SMILES`, `InChI`, `FASTA`, `PDB format`, `MOL format` or `SDF format` in the **code** of any `.rs` file under the scanned roots. Comment-only lines are skipped — see below. |
| `check_toolchain_pins_agree` | §18.1 | The Rust pin is the same in every place that states it. |
| `check_no_platform_transcendentals` | §13.1 | No direct `exp`/`ln`/`sin`/`cos`/`powf` on `f64` — everything routes through `det_math`. |
| `check_no_stream_deriving_method` | §13.1 | `Stream` does not gain a method that would silently fork a sequence. |
| `check_no_closure_predicate_branch` | §13.1 | The shell-closure predicate has not sprouted a second spelling. |
| `check_packing_matches_probe` | §7.1 | The shipped packing model still agrees with the experiment that established it. |

**G4 is not in this list, deliberately.** No-real-world-units is enforced by the
newtypes in `borbax-units` — `Quanta + Thermal` is a compile error — which is
stronger than any grep. But it is not `xtask`, and a table headed "what it
checks" that lists it is claiming a guard that is not here.

**Why the format scan skips comments.** G5 forbids *importing or exporting*
these formats. A parser lives in an identifier, a match arm or a file extension —
never in a comment. Scanning prose as well made it impossible to write the
prohibition down: a doc comment reading "Borbax will never import or export
SMILES, InChI, MOL format or FASTA" failed the gate. A guard that fires on its
own doctrine being documented is one somebody eventually disables, so it now
scans code and lets the doctrine be written.

Things this gate does **not** do, stated because their absence is easy to assume
away:

- **No mapping-table check.** G5 also forbids a Borbax↔real correspondence table;
  nothing here looks for one.
- **Substring matching, case-sensitive.** `MOL` and `PDB` are spelled `MOL
  format` and `PDB format` so they do not fire on `MOLECULE` or on ordinary
  prose — and `parse_smiles` slips through where `parse_SMILES` would not.
- **It cannot scan itself.** `xtask` is outside the scan roots, necessarily: the
  token list would match.

None of that is an argument for widening it until it cries wolf. The value is
that adding a real parser becomes awkward and visible, and §5 sits at the top of
the review precedence with a human on it regardless.

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
