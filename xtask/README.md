# `xtask` — the checks that protect the fiction

Borbax's central promise is that **no real chemistry data is imported into
this repository** (spec §5, revised 2026-08-07 — the promise used to be
broader; see §5 for what changed and why). That promise is easy to state and
easy to erode: one helpful pull request adding "just a small reference table"
and it is gone.

So the promise is enforced by a program, and that program is `xtask`. It runs in
CI and in the pre-commit hook.

## What it checks

**This table covers §5 and its closely adjacent determinism/portability
checks — it is a subset of what `check_guarantees` calls, not the whole
list.** `check_guarantees` runs 17 functions; this table names 9. The other
eight (the viewer-seam checks, lint inheritance, the real-chemistry literal
scan, the signature pin, the wall-clock and `libm` single-home checks) are real
and enforced, and are
documented where they're more legible: CLAUDE.md's seam section, and the doc
comment on each function in `xtask/src/main.rs`, which remains the only
complete and current list.

An earlier version of this table claimed completeness it did not have — listed
G4 as a check when nothing tests it, and separately omitted checks that do
run. Keep this one honest about being partial rather than repeat that: a
README claiming a wider or more complete gate than the code is how a guarantee
comes to be believed rather than enforced.

**§5's G5 was withdrawn on 2026-08-06 and its two checks were deleted with it**
— the chemical-format scan and the palette signature/lookup guard. See the spec
§5 for why. Nothing here enforces G5 any more, and nothing should be added that
does without reopening that decision.

| Function | Guarantee | What it does |
|---|---|---|
| `check_no_data_files` | G1 | No data files under the chemistry crates — no element tables, reaction databases, structures or sequence data. |
| `check_blocklist_present` | G2 | `naming.rs` exists and contains `REAL_ELEMENT_SYMBOLS`. Fails loudly if the file is missing, rather than reporting green for a path that moved. |
| `check_no_crate_escapes_the_scan` | §5 + §13.1 | Every crate directory in the workspace is inside a scanned root, so no crate can be added where the textual checks do not look. |
| `check_no_unscanned_includes` | §5 | No `include!`, `include_str!`, `include_bytes!` or `#[path]` in scanned code — each reaches source the scanners never lex. |
| `check_toolchain_pins_agree` | §18.1 | The Rust pin is the same in every place that states it. |
| `check_no_platform_transcendentals` | §13.1 | No direct `exp`/`ln`/`sin`/`cos`/`powf` on `f64` — everything routes through `det_math`. |
| `check_no_stream_deriving_method` | §13.1 | `Stream` does not gain a method that would silently fork a sequence. |
| `check_no_closure_predicate_branch` | §13.1 | The shell-closure predicate has not sprouted a second spelling. |
| `check_packing_matches_probe` | §7.1 | The shipped packing model still agrees with the experiment that established it. |

**G4 is not in this list, deliberately.** No-real-world-units is enforced by the
newtypes in `borbax-units` — `Quanta + Thermal` is a compile error — which is
stronger than any grep. But it is not `xtask`, and a table headed "what it
checks" that lists it is claiming a guard that is not here.

**`check_no_unscanned_includes` (§5) bans `include!`, `include_str!`,
`include_bytes!` and `#[path]` outright**, rather than trying to scan whatever
they pull in. `include!` of a non-`.rs` file was verified to put a full importer
in the public API with no name of any kind in the `.rs` source, and `#[path]`
into a directory with no manifest walks past the crate-escape check too. Both
were zero occurrences in the tree when banned, so the ban forbids nothing that
exists — it closes a route future code could take, not one it has taken.

**The deleted format check's story is kept in git history, not here**, per the
commit that deleted it — a review once wrote `SmilesParser`, `parse_smiles`,
`read_molfile`, `read_pdb`, `write_fasta` and a real↔Borbax mapping table, and
the gate reported all checks passed with `clippy -D warnings` clean. Worth
knowing only for the transferable half, which the surviving scanners still rely
on: **lex, do not grep.** An identifier is an identifier, a comment is not a
token, and matching on whole lowercase segments is what made `SmilesParser`
findable when a case-sensitive substring scan reported zero hits — ordinary
camelCase, not obfuscation.

Two limits worth stating plainly rather than losing when the story above moved
to git history:

- **A determined obfuscator wins.** `concat!("Car", "bon")` and a string split
  across lines both pass. The job these scanners do is catching accident and
  casual addition, not sabotage.
- **`xtask` cannot scan itself.** It is necessarily outside the roots it
  walks — the vocabulary would match its own definition.

§5 sits at the top of the review precedence with a human on it regardless.

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
