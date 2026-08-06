# `xtask` — the checks that protect the fiction

Borbax's central promise is that **no real chemistry is in this repository**
(spec §5). That promise is easy to state and easy to erode: one helpful pull
request adding "just a small reference table" and it is gone.

So the promise is enforced by a program, and that program is `xtask`. It runs in
CI and in the pre-commit hook.

## What it checks

**This table lists the checks `check_guarantees` actually calls, and nothing
else.** An earlier version described a gate about four times wider than the
code: it listed G4 as a check when nothing tests it, and omitted four checks
that do run. Keep it that way — a README describing a wider gate than the code
is how a guarantee comes to be believed rather than enforced.

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

**Why the format check lexes instead of grepping, and what that cost to learn.**
Two versions of this check failed before the current one:

1. A case-sensitive substring scan for `"SMILES"`, `"InChI"`, `"PDB format"` …
   which read **one file**, so planting `SMILES` in `bonds.rs` passed.
2. The same scan widened to every file — which then fired on *doc comments
   stating the prohibition*, including a sentence lifted from `CLAUDE.md`. Fixed
   by skipping comment-only lines, which quietly **deleted half the list**:
   `"MOL format"`, `"PDB format"` and `"SDF format"` are English phrases that
   occur only in prose.

A **fourth** version followed, because the third was defeated four more times —
each defeat built, compiled and run green. The two that mattered most:

- **Acronym-cased identifiers scored zero.** `SMILESParser`, `PDBReader`,
  `FASTAWriter`, `MMCIFParser` and `InChI` itself all missed, because the
  camelCase split only fired after a lowercase letter, so an acronym absorbed the
  word following it. That made version 3 **strictly weaker than version 1** for
  all-caps spellings — the substring scan had caught them. Fixed with an
  `ACRONYMWord` boundary, plus a second vocabulary tier matched as a substring,
  because `InChI`'s own capitalisation splits to `in`-`ch`-`i` and no boundary
  rule recovers it.
- **A committed `experiments/src/target/` escaped every scan** — `.gitignore`
  anchors `/target/` to the root while the walker pruned the name at any depth.
  That is a **G1** hole, not G5: a real element table planted there passed the
  whole gate, while the identical file one directory up failed as it should.

Plus `include!` of a non-`.rs` file (verified to put a full importer in the
public API with no format name anywhere in the `.rs` source) and `#[path]` into a
directory with no manifest, which walks past the crate-escape check too. Both now
banned outright — there were zero occurrences in the tree, so the ban forbids
nothing that exists.

**The six-format importer story belonged to the deleted check** and is kept in
git history rather than here: a review wrote `SmilesParser`, `parse_smiles`,
`read_molfile`, `read_pdb`, `write_fasta` and a real↔Borbax mapping table, and
the gate reported all checks passed with `clippy -D warnings` clean. It is worth
knowing only for its transferable half, which the surviving scanners still rely
on: **lex, do not grep.** An identifier is an identifier, a comment is not a
token, and matching on whole lowercase segments is what made `SmilesParser`
findable when a case-sensitive substring scan reported zero hits.

Things this gate still does **not** do:

- **It cannot scan itself.** `xtask` is necessarily outside the roots — the
  vocabulary would match its own definition.
- **A determined obfuscator wins.** `concat!("SMI", "LES")` and a string split
  across lines both pass. The job is accident and casual addition, not sabotage.

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
