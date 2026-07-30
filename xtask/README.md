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
| `check_no_real_chemical_formats` | G5 | Lexes every `.rs` file under the scanned roots and rejects identifiers or string literals containing a real chemical-format name as a whole segment — `SmilesParser`, `parse_smiles`, `"pdb"`, `read_molfile`. Comments are not tokens, so doctrine may name them freely. |
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

Then a review wrote a working six-format importer/exporter — `SmilesParser`,
`parse_smiles`, `read_molfile`, `read_pdb`, `write_fasta`, an extension table
`["smi", "smiles", "inchi", "mol", "sdf", "pdb", "fasta"]`, and the real↔Borbax
mapping table G5 also forbids — and the gate reported **all checks passed**,
with `clippy -D warnings` clean. Nobody writing Rust types `SMILESParser`.

`xtask/Cargo.toml` had already written the verdict, about a different check:
*"Hand-rolling was tried twice and failed twice … a textual matcher covers the
shapes someone thought to probe."* `syn` and `proc-macro2` were already
dependencies for that reason.

So it lexes. An identifier is an identifier, a comment is not a token, and
matching is on whole lowercase **segments** — `SmilesParser` splits to
`["smiles", "parser"]`. That parser now produces 11 failures; the doctrine
sentence produces none. Four adversarial tests pin both directions.

Things this gate still does **not** do:

- **No mapping-table check.** G5 forbids a Borbax↔real correspondence table;
  nothing looks for one. The review's `REAL_MAPPING` const passed on its name.
- **`mol` is not in the vocabulary.** It collides with `fn embed(mol: &Molecule)`
  six times in `experiments/`. `sdf` and `molfile` cover the same import path.
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
