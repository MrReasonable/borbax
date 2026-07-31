//! Fiction-guarantee enforcement (spec §5). Run in CI and from the
//! pre-commit hook.
//!
//! These checks make "no real chemistry" an architectural property rather
//! than a promise. They are cheap and they run on every commit.
//!
//! Git hooks are not this binary's job — `prek` installs and runs them from
//! `.pre-commit-config.yaml`, which calls back into `cargo run -p xtask`.
//! `setup` only *invokes* prek; it does not reimplement it.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Extensions that would indicate bulk data smuggled into a chemistry crate.
/// The entire chemistry is generated from a seed; there is nothing legitimate
/// for a data file to be doing in these crates.
const DATA_EXTENSIONS: &[&str] = &["csv", "tsv", "json", "yaml", "yml", "parquet", "bin", "dat"];

/// Crates that must contain no data files at all (G1).
///
/// `experiments` is on the list even though it is not part of the simulation:
/// it implements a working version of Tasks 8-9, so a molecule set or a
/// signature table smuggled in there is the same breach by the same route.
/// Measurement *output* is not exempt either — it belongs in `docs/` or gets
/// regenerated, exactly as it would for the crates below.
const CHEMISTRY_CRATES: &[&str] = &[
    "crates/borbax-universe",
    "crates/borbax-molecule",
    "crates/borbax-reaction",
    "experiments",
];

/// Real chemical-format tokens that must never appear anywhere (G5).
/// §5, G5 — real chemical interchange formats, as **identifier segments**.
///
/// **Matched against lexed identifiers and string literals, case-insensitively
/// and segment by segment — not as substrings of raw source.** The predecessor
/// scanned for `"SMILES"`, `"InChI"`, `"PDB format"` and friends, and a review
/// defeated it completely: a working six-format importer/exporter —
/// `SmilesParser`, `parse_smiles`, `read_molfile`, `read_pdb`, `write_fasta`, an
/// extension table `["smi", "smiles", "inchi", "mol", "sdf", "pdb", "fasta"]`
/// and the real-to-Borbax mapping table G5 also forbids — passed with **zero**
/// hits while `clippy -D warnings` stayed clean. Nobody writing Rust types
/// `SMILESParser`.
///
/// Worse, three of the six were already *dead*: `"MOL format"`, `"PDB format"`
/// and `"SDF format"` are English phrases occurring only in prose, and the
/// commit before this one taught the scan to skip prose. A fix meant to narrow
/// the check had deleted half of it.
///
/// `xtask/Cargo.toml` already carried the verdict, about a different check:
/// *"Hand-rolling was tried twice and failed twice … a textual matcher covers
/// the shapes someone thought to probe."* `syn` and `proc-macro2` are
/// dependencies for that reason, and this check now uses them.
///
/// **`mol` is deliberately absent**: it appears six times in
/// `experiments/src/embed.rs` as `fn embed(mol: &Molecule)`, and `sdf` /
/// `molfile` cover the same import path. `xyz` and `cif` are absent for now but
/// are the obvious next entries once 3D geometry lands.
const FORMAT_SEGMENTS: &[&str] = &["smi", "sdf", "pdb"];

/// Format names long enough to match as a **substring** of a squashed
/// identifier, rather than as a whole segment.
///
/// **Segmentation alone cannot reach `InChI`**: its own capitalisation splits to
/// `in`-`ch`-`i`, so no boundary rule recovers it. Matching these against the
/// identifier with non-alphanumerics removed catches `InChIString`,
/// `MOLFile`, `SDFile` and `PDBx` — measured, 46 of 47 probe spellings, the miss
/// being `mol`, which is a deliberate exclusion.
///
/// The cost, stated because it is real: substring matching means the
/// `inching`/`pinching`/`flinching` family would hit on `inchi`. A
/// segment-prefix variant removes them and loses `InChI`, `MOLFile`, `SDFile`
/// and `PDBx` — measured — so substring is the right trade. Zero hits across
/// every `.rs` file in `crates/` and `experiments/` today.
const FORMAT_SUBSTRINGS: &[&str] = &[
    "smiles", "smarts", "inchi", "molfile", "sdfile", "mmcif", "fasta", "fastq", "pdbx",
];

/// Split an identifier into lowercase segments on `_`, `-`, `.` and camelCase
/// boundaries.
///
/// `SmilesParser` -> `["smiles", "parser"]`; `read_pdb_file` -> `["read", "pdb",
/// "file"]`; `SMILES` -> `["smiles"]`; `b"pdb"` -> `["b", "pdb"]`. Whole-segment
/// matching is what lets `smi` be in the vocabulary without firing on `smith`.
fn identifier_segments(name: &str) -> Vec<String> {
    let chars: Vec<char> = name.chars().collect();
    let mut out = Vec::new();
    let mut cur = String::new();
    for (i, ch) in chars.iter().enumerate() {
        // **Any non-alphanumeric is a boundary.** A review probe found `b"pdb"`
        // scoring zero because the byte-string prefix stayed glued on. ASCII
        // rather than Unicode, to match an ASCII-only vocabulary and to stay
        // consistent with `to_ascii_lowercase` below.
        if !ch.is_ascii_alphanumeric() {
            if !cur.is_empty() {
                out.push(std::mem::take(&mut cur));
            }
            continue;
        }
        let prev = i
            .checked_sub(1)
            .and_then(|j| chars.get(j))
            .copied()
            .unwrap_or('\0');
        let next = chars.get(i + 1).copied().unwrap_or('\0');
        // Three boundaries, and the middle one is the fix for a real defect:
        //   `fooBar`   -> foo|Bar   (lower then upper)
        //   `PDBFile`  -> PDB|File  (upper run, then an upper followed by lower)
        //   `pdb2`     -> pdb|2     (digit run starts or ends)
        // Without the ACRONYMWord case an acronym absorbs the word after it, so
        // `SMILESParser` becomes one segment `smilesparser` and matches nothing
        // — measured, and it made this check *weaker* than the substring scan it
        // replaced, which caught `SMILES` inside it.
        let boundary = !cur.is_empty()
            && ((ch.is_ascii_uppercase() && !prev.is_ascii_uppercase())
                || (ch.is_ascii_uppercase()
                    && prev.is_ascii_uppercase()
                    && next.is_ascii_lowercase())
                || (ch.is_ascii_digit() != prev.is_ascii_digit()));
        if boundary {
            out.push(std::mem::take(&mut cur));
        }
        cur.push(ch.to_ascii_lowercase());
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}
/// §13.1 — float operations that are not specified exactly by IEEE-754, so two
/// correct libm implementations may return different bits for the same input.
/// Every one of these must route through `borbax_units::det_math`.
///
/// **`clippy.toml`'s `disallowed-methods` list is the authority; this is a
/// second pass.** Clippy matches a resolved path, so `x.exp()`, `f64::exp(x)`,
/// `<f64>::sin(x)` and the point-free `.map(f64::exp)` are one rule to it. A
/// text scanner cannot manage that — the first version of this list held only
/// the `.exp()` spelling, and both Task 2 reviewers independently walked a
/// `f64::exp(x)` straight past it while `cargo xtask` printed "all checks
/// passed". The bare-path forms below close the common cases; `<f64>::sin(`
/// and `std::primitive::f64::cos(` are still invisible here and are caught
/// only by clippy.
///
/// Keeping this at all is worth the duplication for two reasons: it runs in
/// the sub-second pre-commit hook where clippy does not, and it reads files
/// as text, so code behind an inactive `cfg` — which clippy never builds and
/// therefore never lints — is still seen.
///
/// `+ - * /` and `sqrt` are absent because IEEE-754 *does* specify them to be
/// correctly rounded; they stay native. `mul_add` is present for a subtler
/// reason: it is one instruction on aarch64 and on x86-64 with FMA enabled,
/// but a dynamic call into the platform's libm on the baseline x86-64 the
/// Linux and Windows CI legs compile for. That is the same hazard as `exp`,
/// and it is why `clippy::suboptimal_flops` — which rewrites `a * b + c` into
/// exactly this — is allowed at workspace level rather than taken.
const BANNED_CALLS: &[&str] = &[
    ".exp()",
    ".exp2()",
    ".exp_m1()",
    ".ln()",
    ".ln_1p()",
    ".log(",
    ".log2()",
    ".log10()",
    ".sin()",
    ".cos()",
    ".tan()",
    ".sin_cos()",
    ".asin()",
    ".acos()",
    ".atan()",
    ".atan2(",
    ".sinh()",
    ".cosh()",
    ".tanh()",
    ".asinh()",
    ".acosh()",
    ".atanh()",
    ".powf(",
    ".powi(",
    ".cbrt()",
    ".hypot(",
    ".mul_add(",
    // Fully-qualified and point-free spellings. No trailing paren: `.map(
    // f64::exp)` has none, and that is the form that reads most naturally
    // over a slice of directions or signature bins.
    "f64::exp",
    "f64::exp2",
    "f64::exp_m1",
    "f64::ln",
    "f64::ln_1p",
    "f64::log",
    "f64::sin",
    "f64::cos",
    "f64::tan",
    "f64::asin",
    "f64::acos",
    "f64::atan",
    "f64::sinh",
    "f64::cosh",
    "f64::tanh",
    "f64::asinh",
    "f64::acosh",
    "f64::atanh",
    "f64::powf",
    "f64::powi",
    "f64::cbrt",
    "f64::hypot",
    "f64::mul_add",
    // Both float widths, matching `clippy.toml`. §13.4 says f64 everywhere in
    // simulation arithmetic, but `borbax-render` may legitimately use f32 for
    // output and its SVG goldens are compared byte-for-byte across the matrix,
    // so an `f32::sin` there diverges just as an f64 one would. The method
    // spellings above (`.exp()` and friends) are type-blind and already catch
    // f32; only the path forms need mirroring.
    "f32::exp",
    "f32::ln",
    "f32::log",
    "f32::sin",
    "f32::cos",
    "f32::tan",
    "f32::asin",
    "f32::acos",
    "f32::atan",
    "f32::sinh",
    "f32::cosh",
    "f32::tanh",
    "f32::asinh",
    "f32::acosh",
    "f32::atanh",
    "f32::powf",
    "f32::powi",
    "f32::cbrt",
    "f32::hypot",
    "f32::mul_add",
    // `max`/`min` in path form only. The method spelling `x.max(y)` is
    // deliberately absent: it is type-blind here, and `a.len().min(b.len())`
    // is `Ord::min` on a `usize`, which is perfectly deterministic. Banning
    // `.min(` as text would fire on every integer clamp in the workspace, and
    // a check that cries wolf gets relaxed. Clippy resolves the type and has
    // no such problem, which is the clearest single case for why it is the
    // authority and this is the second pass.
    "f64::max",
    "f64::min",
    "f32::max",
    "f32::min",
    // Path forms previously "covered" only by prefix accident — `f64::atan`
    // happening to be a substring of `f64::atan2`. Requiring an exact match in
    // the agreement test is what surfaced them.
    "f64::log2",
    "f64::log10",
    "f64::sin_cos",
    "f64::atan2",
    "f32::log2",
    "f32::log10",
    "f32::sin_cos",
    "f32::atan2",
    "f32::exp2",
    "f32::exp_m1",
    "f32::ln_1p",
];

/// §13.4 — parallel iteration, which reorders float reductions.
///
/// `impl Sum for Quanta` is *exactly* the bound `rayon`'s
/// `ParallelIterator::sum` requires, so `iter()` → `par_iter()` is a one-word
/// change that compiles, typechecks, passes clippy and is invisible in review.
/// Measured on the shipped `Sum` shape: five thread counts gave five different
/// answers, while the serial fold was stable at all of them.
///
/// Guarded here rather than in `clippy.toml` because rayon is not a dependency,
/// so there is no path for `disallowed-methods` to resolve — a text scan is the
/// only check available until the day it would be too late to add one.
///
/// Integer reductions are exempt in principle (`Sum for Mass` is associative
/// and exact), but not in this check: the spelling is identical, and the
/// day someone genuinely wants a parallel `Mass` sum is a day for a per-site
/// `#[allow]` and a determinism review, not a hole in a grep.
const BANNED_PARALLEL_CALLS: &[&str] = &[
    ".par_iter()",
    ".par_iter_mut()",
    ".into_par_iter()",
    ".par_bridge()",
    ".par_sort",
    ".par_chunks",
    ".par_extend(",
];

/// Directories scanned for §13.1 violations, relative to the workspace root.
///
/// `experiments` is here for the same reason it is in `CHEMISTRY_CRATES`: it
/// implements a working version of Tasks 8-9 and is meant to be lifted into
/// `borbax-molecule` more or less unchanged, so a direct `.cos()` written
/// there arrives in the simulation later having never been checked.
///
/// `xtask` is *not* scanned. It is a build tool whose output cannot reach a
/// simulation result — and it necessarily contains every banned string above
/// as a literal, so scanning it would mean exempting it, which is worse.
const TRANSCENDENTAL_SCAN_ROOTS: &[&str] = &["crates", "experiments"];

/// Workspace directories that hold crates and are deliberately **not** scanned.
///
/// Only `xtask` — it cannot scan itself, because [`FORMAT_SEGMENTS`] and
/// the banned-call lists would match their own definitions.
const UNSCANNED_CRATE_DIRS: &[&str] = &["xtask"];

fn main() -> Result<(), String> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .ok_or("no workspace root")?
        .to_path_buf();

    // `check-guarantees` is accepted as well as the bare form because CI and
    // the plan both spell it out. An unrecognised argument is an error rather
    // than being ignored: a typo silently running the checks anyway is how a
    // subcommand that never ran gets believed.
    match std::env::args().nth(1).as_deref() {
        None | Some("check-guarantees") => check_guarantees(&root),
        Some("setup") => setup(&root),
        Some(other) => Err(format!(
            "unknown subcommand {other:?} — expected `check-guarantees` or `setup`"
        )),
    }
}

/// Bring a fresh clone or worktree up to working order.
///
/// This exists because of the one real cost of choosing `prek` over a
/// build-script installer like `husky-rs`: hooks are installed by an explicit
/// command, so a fresh clone has none and nothing says so. Rather than accept
/// "remember to run prek install", this is the one command to run after
/// cloning, and it is idempotent — running it again costs nothing.
///
/// It orchestrates rather than duplicates. Hook installation is still prek's
/// job; this only makes sure prek exists and then asks it.
fn setup(root: &Path) -> Result<(), String> {
    if tool_present("prek") {
        println!("prek: already installed");
    } else {
        println!("prek: not found — installing with `cargo install prek --locked`");
        run(root, "cargo", &["install", "prek", "--locked"])?;
    }

    // Bare `install`: the stages come from `default_install_hook_types` in
    // .pre-commit-config.yaml, so the config stays the single source of truth
    // for which hooks exist rather than being restated in an argument here.
    run(root, "prek", &["install"])?;

    println!("setup complete — pre-commit and pre-push hooks active");
    Ok(())
}

/// Whether `bin` can be executed at all.
///
/// Runs `--version` rather than searching PATH by hand, so it works the same
/// on Windows and does not have to know about extensions or shims.
fn tool_present(bin: &str) -> bool {
    Command::new(bin)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn run(root: &Path, bin: &str, args: &[&str]) -> Result<(), String> {
    let status = Command::new(bin)
        .args(args)
        .current_dir(root)
        .status()
        .map_err(|e| format!("could not run `{bin}`: {e}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("`{bin} {}` failed", args.join(" ")))
    }
}

fn check_guarantees(root: &Path) -> Result<(), String> {
    let mut failures = Vec::new();
    check_no_data_files(root, &mut failures)?;
    check_blocklist_present(root, &mut failures)?;
    check_no_real_chemical_formats(root, &mut failures)?;
    check_no_crate_escapes_the_scan(root, &mut failures)?;
    check_no_unscanned_includes(root, &mut failures)?;
    check_toolchain_pins_agree(root, &mut failures)?;
    check_no_platform_transcendentals(root, &mut failures)?;
    check_no_stream_deriving_method(root, &mut failures)?;
    check_no_closure_predicate_branch(root, &mut failures)?;
    check_packing_matches_probe(root, &mut failures)?;

    if failures.is_empty() {
        // Deliberately not an unqualified "all checks passed". The §13.1 scan
        // here is textual and cannot see `<f64>::sin(x)`; its doc says an
        // absence means nothing on its own, and a doc comment on a private
        // function does not reach the person reading hook output. Naming the
        // authority in the success line is the only place that lands.
        println!(
            "repository invariants: all checks passed \
             (§13.1 textually — clippy::disallowed_methods is the authority)"
        );
        Ok(())
    } else {
        for f in &failures {
            eprintln!("FAIL: {f}");
        }
        Err(format!("{} invariant check(s) failed", failures.len()))
    }
}

/// §18.1 — the toolchain is pinned in *three* files, which must all agree.
///
/// `.prototools` installs the toolchain, `rust-toolchain.toml` selects it, and
/// `[workspace.package] rust-version` is the MSRV that makes cargo refuse to
/// build with anything older. All three are needed, for the reasons documented
/// in `rust-toolchain.toml` and `xtask/Cargo.toml`.
///
/// The MSRV was previously left out of this check, and a gap between it and
/// the pin is a hole in the guard exactly the size of that gap. The MSRV
/// exists because proto selects the toolchain with an environment variable, so
/// any context that misses the selection — a GUI git client, an agent shell, a
/// stale `RUSTUP_TOOLCHAIN` — falls back to the machine default. With the MSRV
/// equal to the pin that fails loudly. With the MSRV even one release behind,
/// that context compiles happily on the older rustc and produces different
/// floats, which is the exact failure the MSRV line was added to prevent.
///
/// Keeping the pin *equal* to current stable is the goal; tracking `stable`
/// itself is not. A floating channel would change the compiler under the
/// project without a commit, and every golden hash with it. Renovate proposes
/// bumps as labelled PRs that never automerge, so moving to a new release
/// stays a deliberate act reviewed like a change to physics.
fn check_toolchain_pins_agree(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    let sources = [
        (".prototools `rust`", root.join(".prototools"), "rust"),
        (
            "rust-toolchain.toml `channel`",
            root.join("rust-toolchain.toml"),
            "channel",
        ),
        (
            "Cargo.toml `rust-version`",
            root.join("Cargo.toml"),
            "rust-version",
        ),
    ];

    let mut found: Vec<(&str, String)> = Vec::new();
    for (label, path, key) in &sources {
        match read_pin(path, key)? {
            Some(value) => found.push((label, value)),
            None => failures.push(format!("§18.1: no pin found for {label}")),
        }
    }

    let Some((first_label, first_value)) = found.first() else {
        return Ok(());
    };
    let disagreeing: Vec<String> = found
        .iter()
        .filter(|(_, value)| value != first_value)
        .map(|(label, value)| format!("{label} has {value:?}"))
        .collect();
    if !disagreeing.is_empty() {
        failures.push(format!(
            "§18.1: toolchain pins disagree — {first_label} has {first_value:?}, but {}",
            disagreeing.join(", ")
        ));
    }
    Ok(())
}

/// The value of the first `<key> = "<value>"` line in a simple TOML file.
///
/// Deliberately not a TOML parser: xtask has no dependencies, and every pin is
/// a bare key in a file this repository owns.
///
/// `Cargo.toml` does have tables, which stretches that further than
/// `.prototools` and `rust-toolchain.toml` do. It is still safe today because
/// `rust-version` appears exactly once as a bare key, and the member crates
/// spell theirs `rust-version.workspace`, which does not match. A real parser
/// becomes the right call the moment a second `rust-version` can appear —
/// a `[package.metadata]` table would do it.
fn read_pin(path: &Path, key: &str) -> Result<Option<String>, String> {
    let src = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    for line in src.lines() {
        let line = line.trim();
        if line.starts_with('#') {
            continue;
        }
        let Some((lhs, rhs)) = line.split_once('=') else {
            continue;
        };
        if lhs.trim() != key {
            continue;
        }
        return Ok(rhs.trim().trim_matches('"').to_owned().into());
    }
    Ok(None)
}

/// §13.1 — nothing in `borbax-rng` may derive one `Stream` from another.
///
/// **This is the only enforcement of the crate's central duplicate-stream
/// guarantee.** `Stream::sub` is a constructor precisely so that
/// `base.sub(a).sub(b)` — which under the old `fork` silently returned a
/// duplicate of an unrelated sibling — has no spelling. That held only because
/// nobody had written the method: reinstating it verbatim under any other name
/// passes the whole suite and `clippy -D warnings`, measured.
///
/// A `compile_fail` doctest provably cannot do this job. `sub` cannot be both
/// an associated function and a method (`E0592`), so a regression *must* arrive
/// under a different name than any doctest could spell — and
/// `compile_fail,E0599` does not enforce the error code on the pinned
/// toolchain, so such a block passes when it fails for an unrelated reason.
/// `borbax-units` records that second measurement already.
///
/// **Why `syn` and not a matcher.** Hand-rolled textual versions were written
/// twice. The first caught one signature shape of six; the second, built from
/// the enumerated list of the first's misses, caught six of fifteen — it was
/// blind to `-> (Self, Self)` (which `clippy::use_self` actively pushes authors
/// toward), `-> Option<Self>`, `-> [Stream; 2]`, `-> impl Iterator<Item =
/// Stream>`, free functions, and trait methods (which are spelled `fn`, never
/// `pub fn`, so the scanner never even started). It also *false*-positived on
/// `self_seed` as a parameter name and on a block comment describing the
/// forbidden shape — and a guard that fires on correct code gets deleted.
///
/// The first AST version then failed a third time, and the diagnosis is worth
/// keeping: `syn` performs **no name resolution**, so matching the bare
/// identifier `Stream` was still enumeration — of *type* spellings rather than
/// of *signature* spellings. One line, `type Child = Stream;`, defeated it and
/// compiled clean under `-D warnings`. It also gated the whole check on the
/// enclosing `impl` being `Stream`, so a method on any other type taking a
/// `&Stream` and returning one was invisible; and it read `Self::Item` as a
/// stream, so `impl Iterator for Stream` — an ordinary thing to want from an
/// RNG — was a false positive.
///
/// **What is enforced now, stated as a bound rather than as "by
/// construction".** Local aliases (`type X = Stream`, `use .. as X`) are
/// resolved one level; `Self` is resolved three ways (on `Stream`, on another
/// type, or unknowable in a trait declaration); the return type counts wherever
/// a stream is mentioned, including inside `Option`, tuples, arrays and
/// `impl Trait`; free functions and trait methods are walked; items declared
/// inside function bodies are walked; and anything the AST cannot read — a
/// macro invocation, a verbatim item — is **reported**, because silence must
/// mean "looked and it was clean" rather than "did not look".
///
/// **Known open**, and named rather than implied: a function that mutates a
/// stream in place instead of returning one (`fn rekey(&mut self, sub: u64)`)
/// is outside this framing entirely, and is arguably worse than `fork` — it
/// converts a stream mid-use, so a caller that drew three words and then
/// rekeyed silently aliases a fresh sub-stream at offset 3. An alias of an
/// alias is also unresolved.
///
/// `Clone` is exempt and must be: it duplicates a stream rather than deriving a
/// different one, and `!Copy` exists to make that duplication visible at the
/// call site. A *hand-written* `impl Clone` is not exempt, because it could
/// derive rather than duplicate.
fn check_no_stream_deriving_method(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    let dir = root.join("crates/borbax-rng/src");
    if !dir.is_dir() {
        // A silent `Ok` here is a disabled guarantee that survives a rename or
        // a crate move. If the check cannot run, that is itself the failure.
        failures.push(format!(
            "§13.1: {} is missing — the no-derived-stream check cannot run",
            dir.display()
        ));
        return Ok(());
    }
    for path in walk(&dir)? {
        if path.extension().and_then(|e| e.to_str()) != Some("rs") {
            continue;
        }
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let rel = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .display()
            .to_string();
        for f in scan_for_derived_streams(&rel, &text)? {
            failures.push(f);
        }
    }
    Ok(())
}

/// The predicate behind [`check_no_stream_deriving_method`], split out so it
/// takes `(&str, &str)` and can be unit-tested against a string.
///
/// That shape is the point. `scan_rust_source` has this signature and carries
/// twenty-odd tests; the first two versions of this check took `(root,
/// failures)`, could only be exercised against the real crate — which by
/// construction produces no failures and therefore tests nothing — and each
/// shipped enforcing far less than its own doc claimed.
fn scan_for_derived_streams(rel: &str, text: &str) -> Result<Vec<String>, String> {
    let file = syn::parse_file(text).map_err(|e| format!("{rel}: {e}"))?;
    // `Stream` can be renamed. Collect every local name for it first, because
    // the check is otherwise defeated by one line — `type Child = Stream;` —
    // which compiles clean under `-D warnings`. `syn` does no name resolution,
    // so matching the bare identifier was still enumeration: of *type*
    // spellings rather than of *signature* spellings, which is the same class
    // one level up.
    let mut names: Vec<String> = alloc_stream_names(&file.items);
    names.push("Stream".to_owned());
    let mut out = Vec::new();
    collect_derived_streams(rel, &file.items, &names, SelfTy::Unknown, &mut out);
    out.sort();
    out.dedup();
    Ok(out)
}

/// Local aliases for `Stream`: `type X = Stream;` and `use ... Stream as X;`.
///
/// One pass, not transitive — an alias of an alias is rare enough that the
/// residual is worth naming rather than solving. It is named in
/// [`check_no_stream_deriving_method`]'s doc.
fn alloc_stream_names(items: &[syn::Item]) -> Vec<String> {
    alias_names(items, "Stream")
}

/// Local aliases for `target`: `type X = Target;` and `use ... Target as X;`.
///
/// Generalised from the `Stream`-only version so the signature guard shares one
/// implementation rather than growing a second copy with its own bugs — the
/// enumeration failure this file records three times, one level up.
fn alias_names(items: &[syn::Item], target: &str) -> Vec<String> {
    fn walk_use(tree: &syn::UseTree, target: &str, out: &mut Vec<String>) {
        match tree {
            syn::UseTree::Rename(r) if r.ident == target => out.push(r.rename.to_string()),
            syn::UseTree::Path(p) => walk_use(&p.tree, target, out),
            syn::UseTree::Group(g) => {
                for t in &g.items {
                    walk_use(t, target, out);
                }
            }
            _ => {}
        }
    }
    let mut out = Vec::new();
    for item in items {
        match item {
            syn::Item::Type(t) if idents_of(&t.ty).iter().any(|i| i == target) => {
                out.push(t.ident.to_string());
            }
            syn::Item::Use(u) => walk_use(&u.tree, target, &mut out),
            syn::Item::Mod(m) => {
                if let Some((_, inner)) = &m.content {
                    out.extend(alias_names(inner, target));
                }
            }
            _ => {}
        }
    }
    out
}

/// What `Self` means in the item being inspected.
///
/// Three states, not two. The previous version had `Option<&[String]>` and
/// treated "not an impl on `Stream`" as "stop looking", which missed every
/// method on *another* type that takes a `&Stream` and returns one —
/// `impl Nursery { pub fn child(&self, p: &Stream) -> Stream }` — and
/// false-positived on `impl Iterator for Stream`, where `Self::Item` is a
/// `u64`, and on any trait declaration with `-> Self`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SelfTy {
    /// Inside `impl .. for Stream`: a receiver is a stream and `Self` is one.
    Stream,
    /// Inside `impl .. for` anything else: neither is.
    Other,
    /// A trait declaration or a free function: `Self` is unknowable, so only a
    /// literal `Stream` (or alias) counts. A receiver still counts, because a
    /// trait can be implemented *for* `Stream`.
    Unknown,
}

/// Every identifier appearing in a type, via its token stream.
fn idents_of(ty: &syn::Type) -> Vec<String> {
    fn walk_tt(tt: proc_macro2::TokenTree, out: &mut Vec<String>) {
        match tt {
            proc_macro2::TokenTree::Ident(i) => out.push(i.to_string()),
            proc_macro2::TokenTree::Group(g) => {
                for inner in g.stream() {
                    walk_tt(inner, out);
                }
            }
            _ => {}
        }
    }
    let mut out = Vec::new();
    for tt in quote::ToTokens::to_token_stream(ty) {
        walk_tt(tt, &mut out);
    }
    out
}

/// Does this type hand back a stream?
///
/// `Self` is counted only where it means one. A bare `Self` inside a
/// *projection* — `Self::Item` — does not, which is what lets `Stream`
/// implement `Iterator<Item = u64>` without tripping the guard.
fn yields_stream(ty: &syn::Type, names: &[String], self_ty: SelfTy) -> bool {
    let toks = idents_of(ty);
    if toks.iter().any(|t| names.iter().any(|n| n == t)) {
        return true;
    }
    // `Self` means a stream only inside `impl .. for Stream`. In a trait
    // *declaration* it is the implementor and is unconstrained, so
    // `pub trait Reset { fn reset(&self) -> Self; }` must not fire — a literal
    // `Stream` in that position still does.
    if self_ty != SelfTy::Stream {
        return false;
    }
    // `Self` counts only where it is the WHOLE type. `Self::Item` and
    // `<Self as Keyed>::Key` are projections to something that is not a
    // stream, which is what lets `Stream` implement `Iterator<Item = u64>`.
    // Checked structurally: a string test on the rendered form got
    // `<Self as Keyed>::Key` wrong, because the `::` is not adjacent to `Self`.
    bare_self(ty)
}

/// Is `Self` used as a complete type anywhere in `ty`, rather than as the
/// qualifier of a projection?
fn bare_self(ty: &syn::Type) -> bool {
    match ty {
        syn::Type::Path(p) => {
            // A qualified path `<Self as Trait>::Assoc` is a projection.
            if p.qself.is_some() {
                return false;
            }
            if p.path.segments.len() == 1 && p.path.is_ident("Self") {
                return true;
            }
            // `Self::Item` — more than one segment — is also a projection.
            if p.path.segments.first().is_some_and(|s| s.ident == "Self") {
                return false;
            }
            p.path.segments.iter().any(|seg| match &seg.arguments {
                syn::PathArguments::AngleBracketed(a) => a.args.iter().any(|arg| match arg {
                    syn::GenericArgument::Type(t) => bare_self(t),
                    syn::GenericArgument::AssocType(t) => bare_self(&t.ty),
                    _ => false,
                }),
                syn::PathArguments::Parenthesized(a) => {
                    a.inputs.iter().any(|t| bare_self(&t.ty))
                        || match &a.output {
                            syn::ReturnType::Type(_, t) => bare_self(t),
                            syn::ReturnType::Default => false,
                        }
                }
                syn::PathArguments::None => false,
            })
        }
        syn::Type::Tuple(t) => t.elems.iter().any(bare_self),
        syn::Type::Array(a) => bare_self(&a.elem),
        syn::Type::Slice(s) => bare_self(&s.elem),
        syn::Type::Reference(r) => bare_self(&r.elem),
        syn::Type::Ptr(p) => bare_self(&p.elem),
        syn::Type::Paren(p) => bare_self(&p.elem),
        syn::Type::Group(g) => bare_self(&g.elem),
        syn::Type::ImplTrait(i) => i.bounds.iter().any(|b| match b {
            syn::TypeParamBound::Trait(t) => t.path.segments.iter().any(|seg| {
                matches!(&seg.arguments, syn::PathArguments::AngleBracketed(a)
                if a.args.iter().any(|arg| match arg {
                    syn::GenericArgument::Type(t) => bare_self(t),
                    syn::GenericArgument::AssocType(t) => bare_self(&t.ty),
                    _ => false,
                }))
            }),
            _ => false,
        }),
        _ => false,
    }
}

/// Walk every item, recursing through modules and function bodies, and flag any
/// function that both receives a `Stream` and hands one back.
fn collect_derived_streams(
    rel: &str,
    items: &[syn::Item],
    names: &[String],
    outer: SelfTy,
    out: &mut Vec<String>,
) {
    for item in items {
        match item {
            syn::Item::Mod(m) => {
                if let Some((_, inner)) = &m.content {
                    collect_derived_streams(rel, inner, names, SelfTy::Unknown, out);
                }
            }
            syn::Item::Fn(f) => {
                check_sig(rel, &f.sig, names, outer, out);
                collect_derived_streams(rel, &block_items(&f.block), names, outer, out);
            }
            syn::Item::Impl(i) => {
                let on = if idents_of(&i.self_ty)
                    .iter()
                    .any(|t| names.iter().any(|n| n == t))
                {
                    SelfTy::Stream
                } else {
                    SelfTy::Other
                };
                for it in &i.items {
                    match it {
                        syn::ImplItem::Fn(f) => {
                            check_sig(rel, &f.sig, names, on, out);
                            collect_derived_streams(rel, &block_items(&f.block), names, on, out);
                        }
                        syn::ImplItem::Macro(m) => out.push(unanalysable(rel, &m.mac.path)),
                        syn::ImplItem::Verbatim(_) => {
                            out.push(unanalysable_at(rel, "verbatim impl item"));
                        }
                        _ => {}
                    }
                }
            }
            // A trait method is spelled `fn`, never `pub fn`, so a matcher
            // keyed on visibility never saw an extension trait — the idiomatic
            // way to add a method in Rust, and so the most likely
            // reinstatement rather than the least.
            syn::Item::Trait(t) => {
                for it in &t.items {
                    match it {
                        syn::TraitItem::Fn(f) => {
                            check_sig(rel, &f.sig, names, SelfTy::Unknown, out);
                            // **A default body is a function body.** It can
                            // declare items, and those items can derive
                            // streams — measured, a `fn derive(&Stream) ->
                            // Stream` nested in one produced zero findings.
                            // `Item::Fn` and `ImplItem::Fn` were both already
                            // recursed into; this arm checked the signature
                            // and stopped, which is the same enumeration
                            // mistake one variant further along.
                            if let Some(block) = &f.default {
                                collect_derived_streams(
                                    rel,
                                    &block_items(block),
                                    names,
                                    SelfTy::Unknown,
                                    out,
                                );
                            }
                        }
                        syn::TraitItem::Macro(m) => out.push(unanalysable(rel, &m.mac.path)),
                        syn::TraitItem::Verbatim(_) => {
                            out.push(unanalysable_at(rel, "verbatim trait item"));
                        }
                        _ => {}
                    }
                }
            }
            // "There is code here I cannot read" must not be silence. The
            // missing-directory branch above already resolves that same
            // situation this way.
            syn::Item::Macro(m) => out.push(unanalysable(rel, &m.mac.path)),
            syn::Item::Verbatim(_) => out.push(unanalysable_at(rel, "verbatim item")),
            _ => {}
        }
    }
}

/// The `syn::Item`s declared directly inside a function body.
fn block_items(block: &syn::Block) -> Vec<syn::Item> {
    block
        .stmts
        .iter()
        .filter_map(|s| match s {
            syn::Stmt::Item(i) => Some(i.clone()),
            _ => None,
        })
        .collect()
}

fn unanalysable(rel: &str, path: &syn::Path) -> String {
    let name = path
        .segments
        .last()
        .map_or_else(|| "?".to_owned(), |s| s.ident.to_string());
    unanalysable_at(rel, &format!("`{name}!` expansion"))
}

fn unanalysable_at(rel: &str, what: &str) -> String {
    format!(
        "§13.1: {what} in {rel} cannot be checked for stream derivation — this guard reads the \
         AST and does not expand macros, so silence here would mean \"did not look\" rather than \
         \"looked and it was clean\". Write the impl out, or move it outside borbax-rng."
    )
}

/// Flag `sig` if a `Stream` goes in and a `Stream` comes out.
fn check_sig(
    rel: &str,
    sig: &syn::Signature,
    names: &[String],
    self_ty: SelfTy,
    out: &mut Vec<String>,
) {
    let syn::ReturnType::Type(_, ret) = &sig.output else {
        return;
    };
    if !yields_stream(ret, names, self_ty) {
        return;
    }
    let takes_stream = sig.inputs.iter().any(|arg| match arg {
        // Any receiver spelling — `self`, `&self`, `&mut self`, `self: &Self`,
        // `&'a self` — but only where the receiver *is* a stream.
        syn::FnArg::Receiver(_) => self_ty != SelfTy::Other,
        syn::FnArg::Typed(t) => idents_of(&t.ty)
            .iter()
            .any(|i| names.iter().any(|n| n == i)),
    });
    if !takes_stream {
        return;
    }
    let line = syn::spanned::Spanned::span(sig).start().line;
    let name = &sig.ident;
    out.push(format!(
        "§13.1: `fn {name}` at {rel}:{line} takes a Stream and returns one — that is how a \
         second sub-level silently returns a duplicate of an unrelated sibling. Sub-streams \
         are constructed (`Stream::sub`), not derived."
    ));
}

/// G1 — no real chemistry data enters the repository.
fn check_no_data_files(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    for krate in CHEMISTRY_CRATES {
        let dir = root.join(krate);
        if !dir.exists() {
            continue;
        }
        for entry in walk(&dir)? {
            let Some(ext) = entry.extension().and_then(|e| e.to_str()) else {
                continue;
            };
            if DATA_EXTENSIONS.contains(&ext) {
                failures.push(format!(
                    "G1: data file in chemistry crate: {}",
                    entry.display()
                ));
            }
        }
    }
    Ok(())
}

/// G2 — the blocklist must exist and be non-trivial, or name generation
/// is unguarded and could mint a real element symbol.
fn check_blocklist_present(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    let path = root.join("crates/borbax-universe/src/naming.rs");
    if !path.exists() {
        // **Loud, not `Ok(())`.** This returned success for a missing target
        // until Task 4 created the file, on the reasoning that the file did not
        // exist yet. That reasoning expired the moment it did, and what it
        // leaves behind is a check reporting green for a path that is wrong or
        // a file that has been split — the state the G2 blocklist is least able
        // to afford.
        failures.push(
            "G2: crates/borbax-universe/src/naming.rs is missing, so the blocklist \
             check did not look"
                .into(),
        );
        return Ok(());
    }
    let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
    if !src.contains("REAL_ELEMENT_SYMBOLS") {
        failures.push("G2: naming.rs has no REAL_ELEMENT_SYMBOLS blocklist".into());
    }
    Ok(())
}

/// Every crate in the workspace is inside a scanned root.
///
/// **The scan roots are hand-kept, and a crate added outside them escapes every
/// textual check at once** — G5's format scan, §13.1's transcendental scan and
/// §13.4's parallel-call scan all iterate `TRANSCENDENTAL_SCAN_ROOTS`. A review
/// verified it: a workspace member at `tools/borbax-import/` containing both
/// `pub const SMILES: &str = "SMILES";` and `(-e / t).exp()` produced
/// "all checks passed".
///
/// §13.1 keeps a second guard there — `clippy::disallowed_methods` is
/// warn-by-default and `clippy.toml` is read from the workspace root, so the
/// `exp` would still be caught. **§5 has no second guard at all**, which is what
/// makes this precedence 1 rather than housekeeping.
///
/// So rather than trusting the list, fail when a `Cargo.toml` turns up outside
/// it. That is the same self-maintaining argument
/// `check_no_platform_transcendentals` already makes against a hand-kept crate
/// list, applied one level up.
fn check_no_crate_escapes_the_scan(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    let entries = std::fs::read_dir(root).map_err(|e| e.to_string())?;
    let mut dirs: Vec<std::path::PathBuf> = Vec::new();
    for entry in entries {
        let path = entry.map_err(|e| e.to_string())?.path();
        if path.is_dir() {
            dirs.push(path);
        }
    }
    dirs.sort();
    for dir in dirs {
        let Some(name) = dir.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.starts_with('.') || name == "target" || name == "docs" {
            continue;
        }
        if TRANSCENDENTAL_SCAN_ROOTS.contains(&name) || UNSCANNED_CRATE_DIRS.contains(&name) {
            continue;
        }
        // A directory is a crate root if it, or any child, declares a manifest.
        let mut manifests = vec![dir.join("Cargo.toml")];
        if let Ok(children) = std::fs::read_dir(&dir) {
            for child in children.flatten() {
                manifests.push(child.path().join("Cargo.toml"));
            }
        }
        if manifests.iter().any(|m| m.exists()) {
            failures.push(format!(
                "§5/§13.1: crate directory {name:?} is outside every scan root, so the \
                 G5 format, transcendental and parallel-call checks do not see it. Add it \
                 to TRANSCENDENTAL_SCAN_ROOTS, or to UNSCANNED_CRATE_DIRS with a reason"
            ));
        }
    }
    Ok(())
}

/// §5 — no scanned file may pull in source the scanner cannot see.
///
/// **Two defeats of the format check went through here, both built and run
/// green.** `tables.rs` containing `include!("tables.in")` puts a full
/// six-format importer in the public API while the only literal is `"tables.in"`,
/// which names nothing — confirmed present in `cargo doc` output. And
/// `#[path = "../../tools/importer.rs"]` reaches a file in a directory with no
/// manifest, so `check_no_crate_escapes_the_scan` never fires either.
///
/// A blanket ban is exact rather than blunt: there are **zero** occurrences of
/// any of these in `crates/` or `experiments/` today, so this forbids nothing
/// that exists. Following the literal instead — resolving each path and lexing
/// it — is the better fix if one is ever wanted; it is more code for a case
/// nobody needs yet.
fn check_no_unscanned_includes(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    const FORMS: &[&str] = &["include!", "include_str!", "include_bytes!", "#[path"];
    for scan_root in TRANSCENDENTAL_SCAN_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            continue;
        }
        for entry in walk(&dir)? {
            if entry.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let src = std::fs::read_to_string(&entry).map_err(|e| e.to_string())?;
            let rel = entry.strip_prefix(root).unwrap_or(&entry);
            for (n, line) in src.lines().enumerate() {
                if line.trim_start().starts_with("//") {
                    continue;
                }
                for form in FORMS {
                    if line.contains(form) {
                        failures.push(format!(
                            "§5: {} at {}:{} pulls in source the §5 scanners do not lex. \
                             Move the content into a scanned `.rs` file, or add the target \
                             to the scan explicitly",
                            form,
                            rel.display(),
                            n + 1
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

/// §5, G5 — no real chemical interchange format is read or written anywhere.
///
/// **Lexes rather than greps, for the reason `xtask/Cargo.toml` already
/// records.** Tokenising makes an identifier an identifier and a comment not a
/// token at all, so this can be strict about code and silent about prose
/// without a hand-rolled rule for either — the previous version hand-rolled
/// `starts_with("//")` and got block comments and trailing comments wrong in
/// opposite directions.
///
/// Doc comments arrive as `#[doc = "..."]` and are skipped explicitly; `//` and
/// `/* */` never reach the token stream at all.
fn check_no_real_chemical_formats(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    for scan_root in TRANSCENDENTAL_SCAN_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            continue;
        }
        for entry in walk(&dir)? {
            if entry.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let src = std::fs::read_to_string(&entry).map_err(|e| e.to_string())?;
            let rel = entry.strip_prefix(root).unwrap_or(&entry);
            let Ok(stream) = src.parse::<proc_macro2::TokenStream>() else {
                // **Loud, not skipped.** A file this cannot lex is one it cannot
                // vouch for, and passing it silently is how a scanner reports
                // green over the only file that needed it.
                failures.push(format!("G5: {} could not be tokenised", rel.display()));
                continue;
            };
            let mut hits: Vec<(String, String)> = Vec::new();
            scan_format_tokens(stream, &mut hits);
            for (seg, ctx) in hits {
                failures.push(format!(
                    "G5: real chemical format {seg:?} in {} (as {ctx})",
                    rel.display()
                ));
            }
        }
    }
    Ok(())
}

/// Match one identifier or literal against both tiers of the vocabulary.
///
/// Short names are whole **segments** — `smi` must not fire on `smith`. Long
/// names are substrings of the squashed text, because segmentation cannot reach
/// a name whose own capitalisation splits it (`InChI` -> `in`|`ch`|`i`).
fn check_text(text: &str, ctx: &str, hits: &mut Vec<(String, String)>) {
    for seg in identifier_segments(text) {
        if FORMAT_SEGMENTS.contains(&seg.as_str()) {
            hits.push((seg, ctx.to_owned()));
        }
    }
    let squashed: String = text
        .chars()
        .filter(char::is_ascii_alphanumeric)
        .map(|c| c.to_ascii_lowercase())
        .collect();
    for name in FORMAT_SUBSTRINGS {
        if squashed.contains(name) {
            hits.push(((*name).to_owned(), ctx.to_owned()));
        }
    }
}

/// Walk a token stream, reporting identifiers and string literals that contain a
/// [`FORMAT_SEGMENTS`] entry as a whole segment.
fn scan_format_tokens(stream: proc_macro2::TokenStream, hits: &mut Vec<(String, String)>) {
    for tt in stream {
        match tt {
            proc_macro2::TokenTree::Group(g) => {
                let is_doc = matches!(g.delimiter(), proc_macro2::Delimiter::Bracket)
                    && g.stream().into_iter().next().is_some_and(
                        |t| matches!(&t, proc_macro2::TokenTree::Ident(i) if i == "doc"),
                    );
                if !is_doc {
                    scan_format_tokens(g.stream(), hits);
                }
            }
            proc_macro2::TokenTree::Ident(id) => {
                let name = id.to_string();
                check_text(&name, &format!("identifier `{name}`"), hits);
            }
            proc_macro2::TokenTree::Literal(lit) => {
                // **No `trim_matches` on the quotes.** It stripped `r` from the
                // *content* as well as the prefix, so `"pdbr"` and `"rpdbr"`
                // both scored a false hit on `pdb`. `identifier_segments`
                // already treats `"`, `#` and the `b`/`r`/`c` prefixes as
                // boundaries, so the raw token text is what to pass.
                let text = lit.to_string();
                check_text(&text, &format!("string literal {text}"), hits);
            }
            proc_macro2::TokenTree::Punct(_) => {}
        }
    }
}

/// §13.1 — every transcendental must route through `borbax_units::det_math`.
///
/// A direct call is invisible locally and only diverges on another platform,
/// so a grep at commit time is worth more than a code review. The whole check
/// is textual and knows nothing about types; that is fine, because the names
/// in `BANNED_CALLS` do not appear on anything else in this workspace.
///
/// **This is the second pass, not the authority.** `clippy::disallowed_methods`
/// is, via the path list in the workspace `clippy.toml`; it matches resolved
/// paths and so cannot be fooled by a spelling. What this adds is that it runs
/// in the sub-second pre-commit hook where clippy does not, and that it reads
/// files as text — so code behind an inactive `cfg`, which clippy never builds
/// and therefore never lints, is still seen. Treat a finding here as real and
/// an absence here as meaning nothing on its own.
///
/// **There is no path exemption**, which is a deliberate departure from the
/// obvious design of letting `det_math.rs` call the platform. It does not need
/// to: the wrappers call `libm::exp` and friends as free functions, so the
/// chokepoint is the `libm` crate rather than one privileged file. An
/// exemption that currently exempts nothing is a hole with a sign on it — and
/// the thing the exemption was there to prevent, a second `det_math.rs`
/// appearing in whichever crate happened to want one, is prevented more simply
/// by there being nothing to gain from writing one.
///
/// One known false positive, kept rather than fixed: an integration test under
/// `crates/*/tests/` has no `#[cfg(test)]` attribute, so the region skip never
/// applies to it and a legitimate platform comparison there is reported. The
/// clippy route handles that correctly with a one-line `#[expect]` at the call
/// site, which is where the exemption belongs anyway.
fn check_no_platform_transcendentals(
    root: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    for scan_root in TRANSCENDENTAL_SCAN_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            continue;
        }
        // Scanned wholesale rather than from a list of crate names. A
        // hand-kept list silently exempts every crate added after it was
        // written, which is the same failure mode as a path exemption.
        for entry in walk(&dir)? {
            if entry.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let src = std::fs::read_to_string(&entry).map_err(|e| e.to_string())?;
            let rel = entry.strip_prefix(root).unwrap_or(&entry);
            scan_rust_source(&rel.display().to_string(), &src, failures);
        }
    }
    Ok(())
}

/// Report every banned call in `src` that is not inside a `#[cfg(test)]` item.
///
/// Test code is exempt because proving `det_math` is wired up at all means
/// comparing it against the platform, and that comparison has to call the
/// platform. Tests cannot reach a simulation result, so the exemption costs
/// nothing.
///
/// Finding the end of a `#[cfg(test)]` item means counting braces, and counting
/// braces means knowing which ones are inside comments, strings and character
/// literals. Miscounting in the direction that *over*-counts opening braces
/// would leave the scanner believing it is still inside a test module for the
/// rest of the file — a silent false negative, in the one check whose whole job
/// is to not be silent. So the three pieces of state are asserted at end of
/// file: braces balanced, no string left open, no block comment left open. Any
/// of those failing is reported rather than shrugged off, because it means the
/// file was not really scanned.
fn scan_rust_source(rel: &str, src: &str, failures: &mut Vec<String>) {
    let mut lex = LexState::default();
    let mut depth: i64 = 0;
    // Depth at which the innermost enclosing `#[cfg(test)]` item opened.
    let mut test_region: Option<i64> = None;
    let mut pending_cfg_test = false;

    for (i, raw) in src.lines().enumerate() {
        let code = strip_comments_and_literals(raw, &mut lex);
        let trimmed = code.trim();

        if test_region.is_none() && trimmed.starts_with("#[cfg(test)]") {
            pending_cfg_test = true;
        }

        if test_region.is_none() {
            for call in BANNED_CALLS {
                if code.contains(call) {
                    failures.push(format!(
                        "§13.1: platform transcendental {call} at {rel}:{} — route through \
                         borbax_units::det_math",
                        i + 1
                    ));
                }
            }
            for call in BANNED_PARALLEL_CALLS {
                if code.contains(call) {
                    failures.push(format!(
                        "§13.4: parallel iteration {call} at {rel}:{} — float reductions \
                         must fold in index order",
                        i + 1
                    ));
                }
            }
        }

        let mut opens: i64 = 0;
        let mut closes: i64 = 0;
        for c in code.chars() {
            match c {
                '{' => opens += 1,
                '}' => closes += 1,
                _ => {}
            }
        }

        if pending_cfg_test {
            if opens > 0 {
                test_region = Some(depth);
                pending_cfg_test = false;
            } else if trimmed.ends_with(';') {
                // `#[cfg(test)] use ...;` or `#[cfg(test)] mod tests;` — the
                // attribute is spent without opening a block. Without this the
                // next unrelated `{` in the file would be taken for the start
                // of a test region and everything after it silently skipped.
                pending_cfg_test = false;
            }
        }

        depth += opens - closes;
        if test_region.is_some_and(|start| depth <= start) {
            test_region = None;
        }
    }

    if depth != 0 || !lex.is_clean() {
        failures.push(format!(
            "§13.1: could not reliably scan {rel} (brace depth {depth}, string open {}, block \
             comment open {}, raw string open {}) — the file was not checked",
            lex.in_string,
            lex.in_block_comment,
            lex.raw_hashes.is_some()
        ));
    }
}

/// Whether the scanner is currently inside a construct that spans lines.
#[derive(Default)]
struct LexState {
    in_string: bool,
    in_block_comment: bool,
    /// Hash count of the raw string being consumed, if any. `Some(0)` is `r"`.
    raw_hashes: Option<usize>,
}

impl LexState {
    /// Whether anything is still open. At end of file this must be false, or
    /// the scan reached the end in a state that means it was not really
    /// reading code.
    const fn is_clean(&self) -> bool {
        !self.in_string && !self.in_block_comment && self.raw_hashes.is_none()
    }
}

/// `line` with comments, string literals and character literals removed, so
/// that what remains is code whose braces can be counted and whose method
/// calls are real.
///
/// Character literals are handled only to keep `'{'` from unbalancing the
/// count, and the lookahead is written to leave lifetimes alone — `&'a str`
/// has an apostrophe that never closes, and treating it as a literal would
/// swallow the rest of the line.
///
/// Raw strings are handled because `borbax-render` will be emitting SVG from
/// templates full of `"`, and a raw string with an odd number of inner quotes
/// flips the ordinary string state on. One of those trips the end-of-file
/// assertion, which is loud and survivable; *two* flip the parity back and
/// re-balance the braces, and everything between them is skipped in silence.
fn strip_comments_and_literals(line: &str, lex: &mut LexState) -> String {
    let chars: Vec<char> = line.chars().collect();
    let mut out = String::with_capacity(line.len());
    let mut i = 0;

    while let Some(&c) = chars.get(i) {
        let next = chars.get(i + 1).copied();

        if let Some(hashes) = lex.raw_hashes {
            if c == '"' && (1..=hashes).all(|k| chars.get(i + k) == Some(&'#')) {
                lex.raw_hashes = None;
                i += hashes;
            }
            i += 1;
            continue;
        }

        if lex.in_block_comment {
            if c == '*' && next == Some('/') {
                lex.in_block_comment = false;
                i += 1;
            }
            i += 1;
            continue;
        }

        if lex.in_string {
            if c == '\\' {
                i += 1;
            } else if c == '"' {
                lex.in_string = false;
            }
            i += 1;
            continue;
        }

        match c {
            // Rest of the line is a comment — including `///` and `//!`, which
            // is why a doc example may mention `.exp()` freely.
            '/' if next == Some('/') => break,
            '/' if next == Some('*') => {
                lex.in_block_comment = true;
                i += 2;
            }
            'r' if raw_string_prefix(&chars, i).is_some() => {
                // Safe to re-evaluate: the guard just proved it is `Some`.
                let (hashes, len) = raw_string_prefix(&chars, i).unwrap_or((0, 1));
                lex.raw_hashes = Some(hashes);
                i += len;
            }
            '"' => {
                lex.in_string = true;
                i += 1;
            }
            '\'' => {
                i += char_literal_len(&chars, i).unwrap_or(1);
            }
            _ => {
                out.push(c);
                i += 1;
            }
        }
    }
    out
}

/// `(hash count, prefix length in chars)` if a raw string literal opens at
/// `start`, which must be its `r`.
///
/// The `r` has to start a token, or the `r` in `for` would open one. A `b` or
/// `c` immediately before it is allowed, for `br"..."` and `cr"..."` — Rust
/// has three raw-string prefixes and the first version of this function knew
/// about two, which left `cr#"` reproducing the exact silent-blinding bug the
/// function was written to prevent.
fn raw_string_prefix(chars: &[char], start: usize) -> Option<(usize, usize)> {
    let preceded_by = |k: usize| chars.get(start.checked_sub(k)?).copied();
    if start > 0 {
        let prev = preceded_by(1)?;
        let literal_prefix = (prev == 'b' || prev == 'c')
            && (start < 2 || preceded_by(2).is_none_or(|p| !p.is_alphanumeric() && p != '_'));
        if !literal_prefix && (prev.is_alphanumeric() || prev == '_') {
            return None;
        }
    }
    let mut hashes = 0;
    while chars.get(start + 1 + hashes) == Some(&'#') {
        hashes += 1;
    }
    if chars.get(start + 1 + hashes) == Some(&'"') {
        Some((hashes, hashes + 2))
    } else {
        None
    }
}

/// Length in `char`s of the character literal starting at `start`, or `None`
/// if what starts there is a lifetime rather than a literal.
fn char_literal_len(chars: &[char], start: usize) -> Option<usize> {
    // `'\n'`, `'\u{1F600}'` — escaped, so scan for the closing quote. The bound
    // is generous enough for the longest unicode escape and short enough that a
    // lifetime followed by an unrelated quote later on the line cannot match.
    let escaped = chars.get(start + 1) == Some(&'\\');
    let limit = if escaped { 12 } else { 3 };
    for offset in 2..=limit {
        match chars.get(start + offset) {
            Some('\'') => return Some(offset + 1),
            Some(_) if escaped || offset < 3 => {}
            _ => return None,
        }
    }
    None
}

/// Every file under `dir`, without following symlinks.
///
/// The distinction matters in both directions, which is why this uses
/// `symlink_metadata` rather than the `Path::is_dir` that follows links:
///
/// - A symlink pointing at an ancestor is an infinite descent. `is_dir()`
///   resolves the target, so the walk would recurse until it ran out of memory.
/// - A symlink is still *reported*, because the G1 check keys off the path's
///   own extension. A `crates/borbax-universe/elements.json` symlinked to a
///   real table outside the crate is exactly the smuggling route this check
///   exists to close, and skipping links entirely would have opened it.
fn walk(dir: &Path) -> Result<Vec<PathBuf>, String> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        for entry in std::fs::read_dir(&d).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            let meta = std::fs::symlink_metadata(&path).map_err(|e| e.to_string())?;
            if meta.is_symlink() {
                out.push(path);
            } else if meta.is_dir() {
                // **Prune cargo's build directory, not every directory called
                // `target`.** `.gitignore` anchors `/target/`, so a nested
                // `experiments/src/target/` is *committed* — and pruning by name
                // at any depth meant it was never scanned. Measured: a real
                // element table at `experiments/src/target/elements.json` passed
                // the whole gate, while the identical file one directory up
                // failed `check_no_data_files` as it should. That is a **G1**
                // hole, not merely G5, because `walk` is shared with the
                // data-file check.
                let is_build_dir = path.file_name().and_then(|n| n.to_str()) == Some("target")
                    && path.parent().and_then(|p| p.file_name()).is_none();
                if !is_build_dir {
                    stack.push(path);
                }
            } else {
                out.push(path);
            }
        }
    }
    // Sort so failures are reported in a stable order.
    out.sort();
    Ok(out)
}

/// Files whose frontier path must not branch on the closure predicate.
const CLOSURE_PREDICATE_FILES: &[&str] = &[
    "crates/borbax-universe/src/packing.rs",
    "crates/borbax-universe/src/element.rs",
];

/// Functions exempted from [`scan_closure_predicates`], **by name and with a
/// reason**.
///
/// `compact_bound` contains `if a == 0`, and `a = min(outer, cap - outer)`, so
/// that *is* a branch on the closure predicate. It is a legitimate guard
/// against `sqrt(12*0 - 3)`, and measured, it alone holds zero-at-closure when
/// the other factor is broken. A whole-file grep either flags it forever or is
/// switched off and stops seeing the real thing.
///
/// The exemption is bounded by *output* as well as by name:
/// `unmade_lateral_is_finite_and_non_negative_everywhere` asserts the frontier
/// count is finite and non-negative for every `(cap, outer)` the table reaches.
/// Without that, `NaN` and negatives escape the `min` in `unmade_lateral` and
/// both saturate to 0 through `round_ties_even() as u8` — so
/// `closed_shells_have_valence_zero` would read a correct zero while
/// `contacts_upto` was poisoned.
const CLOSURE_PREDICATE_EXEMPT_FNS: &[&str] = &["compact_bound"];

/// Variables that name the fill count on the frontier path.
const FRONTIER_VARS: &[&str] = &["outer", "a", "take", "fill", "filled", "f", "smaller"];

/// Values that mean "this shell is closed".
const CLOSURE_SENTINELS: &[&str] = &[
    "0", "0.0", "0_usize", "0.0_f64", "1", "cap", "capacity", "units",
];

/// Comparison operators, longest first so `<=` is never read as `<`.
const COMPARISONS: &[&str] = &["==", "!=", "<=", ">=", "<", ">"];

/// True if a comparison ending at `after` is the *whole* predicate rather than a
/// subexpression.
///
/// Without this, `if outer < cap - outer` — the legitimate `smaller` computation
/// in `unmade_lateral` — matches `outer < cap`. It was the only false positive in
/// the prototype, so it is the thing to keep a test on.
fn terminates(code: &str, after: usize) -> bool {
    code.get(after..)
        .and_then(|rest| rest.trim_start().chars().next())
        .is_none_or(|c| matches!(c, '{' | '}' | ';' | ')' | ',' | '&' | '|'))
}

/// The trailing identifier on `left`, if any.
fn last_token(left: &str) -> &str {
    let t = left.trim_end();
    let start = t
        .rfind(|c: char| !(c.is_alphanumeric() || c == '_' || c == '.'))
        .map_or(0, |i| i + 1);
    t.get(start..).unwrap_or("")
}

/// The leading identifier on `right`, if any.
fn first_token(right: &str) -> &str {
    let t = right.trim_start();
    let end = t
        .find(|c: char| !(c.is_alphanumeric() || c == '_' || c == '.'))
        .unwrap_or(t.len());
    t.get(..end).unwrap_or("")
}

/// A closure-predicate comparison on `code`, in **either operand order**.
///
/// **Operand order, not another pattern.** The fixed-string list this replaces
/// matched `outer == 0` and missed `0 == outer` — the same predicate with its
/// operands swapped — and the emergence auditor demonstrated end to end that the
/// Yoda spelling walks past every gate while a load-bearing declaration sits
/// behind it. Measured: `rustfmt` normalises the *spacing* (`outer==0` becomes
/// `outer == 0`) but never reorders operands, so that bypass was permanently
/// fmt-clean. Matching a comparison between one [`FRONTIER_VARS`] entry and one
/// [`CLOSURE_SENTINELS`] entry, in either position, closes it once rather than
/// one spelling at a time.
fn closure_predicate_hit(code: &str) -> Option<String> {
    for op in COMPARISONS {
        let mut from = 0;
        while let Some(rel) = code.get(from..).and_then(|c| c.find(op)) {
            let at = from + rel;
            let after = at + op.len();
            // `<` and `>` must not fire inside `<=`, `>=`, `==`, `!=`.
            let bytes = code.as_bytes();
            let glued = (op.len() == 1 && bytes.get(after) == Some(&b'='))
                || (at > 0 && matches!(bytes.get(at - 1), Some(b'=' | b'!' | b'<' | b'>')));
            if glued {
                from = after;
                continue;
            }
            let (lhs, rhs) = (
                last_token(code.get(..at).unwrap_or("")),
                first_token(code.get(after..).unwrap_or("")),
            );
            let pair_ok = (FRONTIER_VARS.contains(&lhs) && CLOSURE_SENTINELS.contains(&rhs))
                || (CLOSURE_SENTINELS.contains(&lhs) && FRONTIER_VARS.contains(&rhs));
            if pair_ok && terminates(code, after + rhs.len() + 1) {
                return Some(format!("{lhs} {op} {rhs}"));
            }
            from = after;
        }
    }
    None
}

/// §7.1 — no branch in the frontier path may read the closure predicate.
///
/// The claim this protects is that zero-valence-at-closure is *arithmetic*,
/// not a declaration. Neither behavioural test can establish it: deletion is
/// uninformative, because two independent factors hold the zero and mutating
/// either alone leaves a zero-at-closure assertion passing; and continuity was
/// shown defeatable by a load-bearing `if outer == 0` that passed every test in
/// both files while producing a bit-identical table. So the property is checked
/// where it lives — in the source.
///
/// **This is a tripwire, not a proof, and the class it enumerates is textual
/// shapes that name the frontier variables.** Outside that class, and therefore
/// invisible here: a predicate computed into a differently-named binding
/// (`let done = outer == 0` is caught, `let done = k_outer(); if done` is not),
/// a branch reached through a helper, a `const` comparison written as
/// `outer.eq(&0)`, and — most importantly — a change to the *formula* rather
/// than an added branch. A mutation that breaks `f(1-f)` produces no new
/// branch and this check stays green.
///
/// That paragraph is the point of the step. The previous three repairs to the
/// §13.1 guard each moved to a more principled-*looking* mechanism without
/// naming what the new mechanism decides, and each shipped enforcing less than
/// its doc claimed. Pair this with `closed_shells_have_valence_zero`, which
/// guards the output; neither substitutes for the other.
fn check_no_closure_predicate_branch(
    root: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    for rel in CLOSURE_PREDICATE_FILES {
        let path = root.join(rel);
        // **Loudly, not `Ok(())`.** `check_blocklist_present` returned success
        // for a missing target, which reports green for a path that is wrong or
        // a file that has been split — worse than no check at all.
        if !path.exists() {
            failures.push(format!(
                "§7.1: {rel} is missing, so the closure-predicate check did not look"
            ));
            continue;
        }
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        scan_closure_predicates(rel, &src, failures);
    }
    Ok(())
}

/// Report every closure-predicate branch in `src` outside an exempt function
/// and outside `#[cfg(test)]`.
///
/// Function scoping is by the innermost `fn` name at the current brace depth,
/// which is what lets `compact_bound` be exempted without switching the check
/// off for the file it lives in.
fn scan_closure_predicates(rel: &str, src: &str, failures: &mut Vec<String>) {
    let mut lex = LexState::default();
    let mut depth: i64 = 0;
    let mut test_region: Option<i64> = None;
    let mut pending_cfg_test = false;
    // (depth at which the fn body opened, name).
    let mut fn_stack: Vec<(i64, String)> = Vec::new();
    let mut pending_fn: Option<String> = None;

    for (i, raw) in src.lines().enumerate() {
        let code = strip_comments_and_literals(raw, &mut lex);
        let trimmed = code.trim();

        if test_region.is_none() && trimmed.starts_with("#[cfg(test)]") {
            pending_cfg_test = true;
        }
        if let Some(name) = fn_name_of(trimmed) {
            pending_fn = Some(name);
        }

        if test_region.is_none() {
            // A line that *declares* a function is attributed to that function,
            // not to its parent. Otherwise a single-line nested item —
            // `#[inline] fn inner(outer) { if outer == 0 { .. } }` inside
            // `compact_bound` — is scanned under the parent's frame and
            // inherits its exemption, because the frame is only pushed when the
            // brace is processed at the end of the line.
            let current = pending_fn
                .as_deref()
                .unwrap_or_else(|| fn_stack.last().map_or("", |(_, n)| n.as_str()));
            if !CLOSURE_PREDICATE_EXEMPT_FNS.contains(&current)
                && let Some(pat) = closure_predicate_hit(&code)
            {
                failures.push(format!(
                    "§7.1: {rel}:{} branches on the closure predicate ({pat:?}) \
                         inside `{}` — zero at a closure must be arithmetic, not a \
                         declaration",
                    i + 1,
                    if current.is_empty() {
                        "<file scope>"
                    } else {
                        current
                    },
                ));
            }
        }

        for c in code.chars() {
            match c {
                '{' => {
                    depth += 1;
                    if pending_cfg_test && test_region.is_none() {
                        test_region = Some(depth);
                        pending_cfg_test = false;
                    }
                    if let Some(name) = pending_fn.take() {
                        fn_stack.push((depth, name));
                    }
                }
                '}' => {
                    if test_region == Some(depth) {
                        test_region = None;
                    }
                    if fn_stack.last().is_some_and(|(d, _)| *d == depth) {
                        fn_stack.pop();
                    }
                    depth -= 1;
                }
                ';' => {
                    // `#[cfg(test)] use ..;` and a bodyless `fn f();` never
                    // open a block, so a pending marker must not survive.
                    pending_cfg_test = false;
                    pending_fn = None;
                }
                _ => {}
            }
        }
    }

    // **The end-of-scan assertion, inherited from `scan_for_derived_streams`
    // and missing here for one review round.** `strip_comments_and_literals`
    // does not handle *nested* block comments, which Rust accepts. One placed
    // inside the `#[cfg(test)]` module inflates `depth`, so the module's
    // closing brace never matches `test_region`, the region stays latched, and
    // every subsequent line goes unscanned — reporting zero findings and zero
    // failures. Same outcome for an unterminated string.
    //
    // That function's own doc claims such a construct "trips the end-of-file
    // assertion, which is loud and survivable". Without this block that
    // sentence was false for this caller.
    if depth != 0 || !lex.is_clean() {
        failures.push(format!(
            "§7.1: {rel} could not be reliably scanned (end depth {depth}, \
             lexer {}) — the closure-predicate check did not look at the whole file",
            if lex.is_clean() {
                "clean"
            } else {
                "left a string or block comment open"
            }
        ));
    }
}

/// The name in a `fn NAME(` declaration on `line`, if there is one.
///
/// **A rejection rule, not an accept-list of modifiers.** The first version
/// accepted only `pub`, `)`, `const`, `async`, `unsafe` and `extern` before the
/// `fn`, which silently declined to register anything it had not heard of —
/// `#[inline] fn` among them. Failing to register is normally safe, because an
/// unregistered function leaves the scan at file scope, which is not exempt.
/// It is unsafe in exactly one direction: a nested item inside an *exempt*
/// function inherits the exemption. Measured, `#[inline] fn inner(outer) { if
/// outer == 0 ... }` nested inside `compact_bound` produced zero findings.
///
/// So anything at a word boundary counts as a declaration unless it is plainly
/// a type (`-> fn(..)`, `: fn(..)`), and a declaration whose name cannot be
/// parsed still opens a frame under [`UNNAMED_FN`] — which is not exempt.
fn fn_name_of(line: &str) -> Option<String> {
    let mut search = 0;
    loop {
        let idx = line.get(search..)?.find("fn ")? + search;
        let prev = line.get(..idx).and_then(|s| s.chars().next_back());
        if prev.is_some_and(|c| c.is_alphanumeric() || c == '_') {
            // `sfn `, `my_fn ` — an identifier that happens to end in `fn`.
            search = idx + 3;
            continue;
        }
        let before = line.get(..idx)?.trim_end();
        if before.ends_with("->") || before.ends_with(':') {
            // A function *type*, not a declaration.
            return None;
        }
        let rest = line.get(idx + 3..)?;
        let name = rest
            .find(['(', '<', ' ', '\r'])
            .and_then(|end| rest.get(..end))
            .map(str::trim)
            .filter(|n| !n.is_empty() && n.chars().all(|c| c.is_alphanumeric() || c == '_'));
        return Some(name.map_or_else(|| UNNAMED_FN.to_owned(), str::to_owned));
    }
}

/// Frame name for a `fn` whose declaration parsed but whose name did not.
///
/// Deliberately not a valid Rust identifier, so it can never collide with an
/// entry in [`CLOSURE_PREDICATE_EXEMPT_FNS`] and therefore never inherits an
/// exemption.
const UNNAMED_FN: &str = "<unnamed fn>";

/// The probe and `packing.rs` share these function bodies verbatim.
const SHARED_PACKING_FNS: &[&str] = &[
    "shell_size",
    "lateral_coordination",
    "unmade_lateral",
    "continuum_at",
    "frontier_notches",
    "compact_bound",
    "lateral_made",
    "lateral_made_raw",
];

/// Constants the probe and `packing.rs` must hold in common.
///
/// **Bodies are not enough, and this is the gap two lanes found independently.**
/// `unmade_lateral`'s body reads `FRONTIER_COEFF` *by name*, so the two copies
/// can hold different values while their texts stay identical. Measured: with
/// `packing.rs` at 6.50 and the probe at 6.90 — a 6% divergence in the constant
/// that sets every element's valence — `cargo xtask` reported all checks passed
/// and all 28 crate tests passed. 7.20 is caught, but only incidentally, by
/// `valence_ceiling_is_four_to_six_over_the_drawn_range`; roughly [6.4, 7.1] was
/// caught by nothing at all.
const SHARED_PACKING_CONSTS: &[&str] = &["FRONTIER_COEFF"];

/// Step 8d — pin `packing.rs` to the probe it was measured from.
///
/// **Every number Task 4 asserts was measured in `experiments/src/bin/fusion.rs`,
/// and the two copies have already diverged four times.** `compact_bound`
/// carried a clamp in one and not the other under a doc saying no clamp is
/// written; `unmade_lateral` bound `compact_bound(smaller)` twice where the
/// probe binds it once; the strain comment kept a defence the probe had
/// deleted; and the boundary-contact clause said "three of six" after the probe
/// said two. None of these is visible to a compiler, a test, or a reviewer
/// reading one file.
///
/// **Text, not numbers.** Comparing outputs would pass for two implementations
/// that agree on the sampled domain and differ off it — the defect shape this
/// project keeps hitting, most recently a proposed fix whose replacement and
/// original were the same function over the tested domain (rho = 1.0000). The
/// text is the artefact that drifts, so the text is what is checked.
///
/// **What this does not see, stated rather than implied:** comments are
/// stripped before comparison, so the *third* historical divergence — a comment
/// keeping a defence the probe deleted — would not be caught. Bodies only, and
/// only for the functions named in [`SHARED_PACKING_FNS`]; `contacts_upto` is
/// excluded because the two take different constant types.
///
/// When the probe is retired, delete this check in the same commit and say so —
/// a check silently passing because its input vanished is the missing-file
/// failure in a different hat.
fn check_packing_matches_probe(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    let ours = root.join("crates/borbax-universe/src/packing.rs");
    let probe = root.join("experiments/src/bin/fusion.rs");
    for (label, path) in [("packing.rs", &ours), ("fusion.rs", &probe)] {
        if !path.exists() {
            failures.push(format!(
                "Step 8d: {label} is missing, so the probe-divergence check did not look"
            ));
            return Ok(());
        }
    }
    let a = std::fs::read_to_string(&ours).map_err(|e| e.to_string())?;
    let b = std::fs::read_to_string(&probe).map_err(|e| e.to_string())?;

    for name in SHARED_PACKING_CONSTS {
        match (extract_const_value(&a, name), extract_const_value(&b, name)) {
            (Some(x), Some(y)) if x == y => {}
            (Some(x), Some(y)) => failures.push(format!(
                "Step 8d: `{name}` has diverged from the probe.\n         \
                 packing.rs: {x}\n         fusion.rs : {y}"
            )),
            (None, _) => failures.push(format!(
                "Step 8d: const `{name}` not found in packing.rs — renamed, or inlined"
            )),
            (_, None) => failures.push(format!(
                "Step 8d: const `{name}` not found in fusion.rs — the probe has dropped it"
            )),
        }
    }

    for name in SHARED_PACKING_FNS {
        match (extract_fn_body(&a, name), extract_fn_body(&b, name)) {
            (Some(x), Some(y)) if x == y => {}
            (Some(x), Some(y)) => failures.push(format!(
                "Step 8d: `{name}` has diverged from the probe.\n         \
                 packing.rs: {x}\n         fusion.rs : {y}"
            )),
            (None, _) => failures.push(format!(
                "Step 8d: `{name}` not found in packing.rs — renamed, or the extractor is wrong"
            )),
            (_, None) => failures.push(format!(
                "Step 8d: `{name}` not found in fusion.rs — the probe has dropped it"
            )),
        }
    }
    Ok(())
}

/// The right-hand side of `const NAME: ... = ...;` in `src`.
///
/// Comment- and literal-stripped before searching, for the same reason
/// [`extract_fn_body`] is: a `const NAME` mentioned in prose must not be
/// mistaken for the declaration.
///
/// Named constants only. A value inlined as a literal into a body is compared
/// by [`extract_fn_body`] instead; a value that is neither is invisible to this
/// check, and that is the honest limit of it.
fn extract_const_value(src: &str, name: &str) -> Option<String> {
    let mut lex = LexState::default();
    for raw in src.lines() {
        let code = strip_comments_and_literals(raw, &mut lex);
        let t = code.trim();
        // Visibility-agnostic. An earlier version enumerated `const ` and
        // `pub const `, the two spellings that happened to be in front of me,
        // and broke the moment `packing` was narrowed to `pub(crate)`. It failed
        // *loudly* — "not found in packing.rs — renamed, or inlined" — which is
        // the only reason it cost two minutes instead of shipping as a silent
        // pass, and is the whole argument for the missing-target failures above.
        let Some(rest) = t
            .split_once("const ")
            .filter(|(before, _)| before.is_empty() || before.trim_end().starts_with("pub"))
            .map(|(_, rest)| rest)
        else {
            continue;
        };
        let Some(after) = rest.strip_prefix(name) else {
            continue;
        };
        if !after.trim_start().starts_with(':') {
            continue;
        }
        let value = after.split_once('=')?.1.trim().trim_end_matches(';').trim();
        return Some(value.to_owned());
    }
    None
}

/// The body of `fn name` in `src`, normalised for comparison.
///
/// Normalisation removes exactly the differences that are legitimate between a
/// library and a single-file probe — `pub`/`pub(super)`, `#[must_use]`, comments
/// and whitespace — and nothing else.
///
/// **Strip first, then search. The reverse order made this check pass
/// vacuously.** The first version located both `fn NAME(` and its opening brace
/// on *raw* source and lexed only what followed, so a `{` appearing in a
/// comment between the two silently redefined "the body". Because the two files
/// are deliberate near-copies *including their prose*, it redefined it
/// identically in both: two functions differing only in their bodies extracted
/// to the same string and compared equal. That is the fifth divergence, missed
/// by the check written to stop the first four.
///
/// **`None` on imbalance, never a truncated body.** Unbalanced braces mean the
/// extraction failed and the caller has a branch that says so; returning `Some`
/// of whatever was read makes that branch unreachable and silently compares two
/// wrong strings.
fn extract_fn_body(src: &str, name: &str) -> Option<String> {
    let mut lex = LexState::default();
    let flat: String = src
        .lines()
        .map(|raw| strip_comments_and_literals(raw, &mut lex))
        .collect::<Vec<_>>()
        .join("\n");

    let at = flat.find(&format!("fn {name}("))?;
    let open = flat.get(at..)?.find('{')? + at;

    let mut depth = 0_i64;
    let mut end = None;
    for (i, c) in flat.get(open..)?.char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    end = Some(open + i + 1);
                    break;
                }
            }
            _ => {}
        }
    }
    let body = flat.get(open..end?)?;
    Some(body.split_whitespace().collect::<Vec<_>>().join(" "))
}

/// The §13.1 scanner is the only piece of logic here that can fail *quietly* —
/// every other check either finds its target or does not exist. A miscounted
/// brace makes it skip the rest of a file while still reporting success, so
/// each way it could miscount gets a test.
#[cfg(test)]
mod tests {
    use super::scan_for_derived_streams;
    use super::{extract_const_value, extract_fn_body, scan_closure_predicates};
    use super::{identifier_segments, scan_format_tokens};

    /// A `{` inside a comment between `fn NAME(` and the real body must not be
    /// mistaken for the body brace. The two files are near-copies including
    /// their prose, so the mistake lands identically in both and two differing
    /// bodies compare equal — the check passing vacuously.
    #[test]
    fn a_brace_in_a_comment_does_not_become_the_body() {
        let lib = "/// Mirrors `fn f(a)`; see the {7.2} note.\n\
                   pub fn f(a: usize) -> f64 { 6.0 * (a as f64) }\n";
        let probe = "/// Mirrors `fn f(a)`; see the {7.2} note.\n\
                     fn f(a: usize) -> f64 { 1.0 * (a as f64) }\n";
        assert_ne!(
            extract_fn_body(lib, "f"),
            extract_fn_body(probe, "f"),
            "a comment brace was taken for the body, so two different bodies compared equal"
        );
    }

    /// Unbalanced braces mean extraction failed. Returning `Some` of whatever
    /// was read makes the caller's "the extractor is wrong" branch unreachable.
    #[test]
    fn an_unterminated_body_extracts_to_none() {
        assert_eq!(
            extract_fn_body("fn f(a: usize) -> f64 {\n    6.0\n", "f"),
            None
        );
    }

    /// ...and extraction must stop at its own closing brace, not run on.
    #[test]
    fn extraction_stops_at_the_end_of_the_function() {
        let body = extract_fn_body("fn f() { 6.0 } fn g() { 9.9 }\n", "f");
        assert_eq!(body.as_deref(), Some("{ 6.0 }"));
    }

    /// **The Yoda spelling was a permanently fmt-clean bypass**, demonstrated
    /// end to end by the emergence auditor: `if 0 == outer` walked past the
    /// fixed-string list while a load-bearing declaration sat behind it and all
    /// six gate legs stayed green. It is not another pattern, it is the same
    /// predicate with its operands swapped.
    #[test]
    fn the_predicate_is_matched_in_either_operand_order() {
        for src in [
            "fn f(cap: usize, outer: usize) -> f64 { if outer == 0 { return 0.0; } 1.0 }",
            "fn f(cap: usize, outer: usize) -> f64 { if 0 == outer { return 0.0; } 1.0 }",
            "fn f(cap: usize, outer: usize) -> f64 { if 0 >= outer { return 0.0; } 1.0 }",
            "fn f(cap: usize, outer: usize) -> f64 { if outer <= 0 { return 0.0; } 1.0 }",
            "fn f(cap: usize, outer: usize) -> f64 { if cap == outer { return 0.0; } 1.0 }",
            "fn f(cap: usize, outer: usize) -> f64 { let done = 0 == outer; 1.0 }",
        ] {
            assert_eq!(closure_hits(src).len(), 1, "missed a spelling: {src}");
        }
    }

    /// ...and must not fire on the legitimate lines of the real file. The
    /// `smaller` computation is the one that matters: without the
    /// whole-predicate check, `outer < cap - outer` matches `outer < cap`.
    #[test]
    fn the_matcher_does_not_fire_on_legitimate_frontier_code() {
        for src in [
            "fn f(cap: usize, outer: usize) -> usize { let smaller = if outer < cap - outer { outer } else { cap - outer }; smaller }",
            "fn f(remaining: usize, cap: usize) -> usize { let take = if remaining < cap { remaining } else { cap }; take }",
            "fn f(continuum: f64, discrete: f64) -> f64 { if continuum < discrete { continuum } else { discrete } }",
        ] {
            assert!(closure_hits(src).is_empty(), "false positive on: {src}");
        }
    }

    /// A nested item inside an exempt function must not inherit the exemption.
    #[test]
    fn a_nested_item_does_not_inherit_the_exemption() {
        let hits = closure_hits(
            "fn compact_bound(a: usize) -> f64 {\n\
             \x20   #[inline] fn inner(outer: usize) -> f64 { if outer == 0 { 0.0 } else { 1.0 } }\n\
             \x20   inner(a)\n}\n",
        );
        assert_eq!(
            hits.len(),
            1,
            "the exemption leaked into a nested item: {hits:?}"
        );
    }

    /// An identifier ending in `fn` is not a declaration.
    #[test]
    fn an_identifier_ending_in_fn_is_not_a_declaration() {
        assert_eq!(super::fn_name_of("let my_fn = 1;"), None);
        assert_eq!(
            super::fn_name_of("fn shell_size(k: usize)"),
            Some("shell_size".to_owned())
        );
        assert_eq!(
            super::fn_name_of("#[inline] fn inner(a: usize)"),
            Some("inner".to_owned())
        );
        assert_eq!(
            super::fn_name_of("pub(crate) fn f(a: usize)"),
            Some("f".to_owned())
        );
    }

    /// A construct the lexer cannot follow must be reported, not silently
    /// swallow the rest of the file. Nested block comments are valid Rust and
    /// `strip_comments_and_literals` does not handle them.
    #[test]
    fn an_unscannable_file_is_reported_rather_than_skipped() {
        let hits = closure_hits(
            "#[cfg(test)]\nmod tests {\n\
             \x20   /* a /* nested */ {{{ */\n\
             \x20   fn t() {}\n}\n\
             fn frontier_notches(cap: usize, outer: usize) -> f64 {\n\
             \x20   if outer == 0 { return 0.0; }\n\
             \x20   1.0\n}\n",
        );
        assert!(
            !hits.is_empty(),
            "a nested block comment latched the test region and the rest of the file went unscanned"
        );
    }

    /// A constant read by name inside a pinned body is invisible to a body
    /// comparison, so it needs its own check.
    #[test]
    fn a_divergent_constant_is_extracted_and_compared() {
        let lib = "pub const FRONTIER_COEFF: f64 = 6.90;\n";
        let probe = "const FRONTIER_COEFF: f64 = 6.50;\n";
        assert_eq!(
            extract_const_value(lib, "FRONTIER_COEFF").as_deref(),
            Some("6.90")
        );
        assert_ne!(
            extract_const_value(lib, "FRONTIER_COEFF"),
            extract_const_value(probe, "FRONTIER_COEFF")
        );
    }

    /// Every visibility spelling is a declaration. Enumerating the two I had in
    /// front of me broke the moment `packing` was narrowed to `pub(crate)`.
    #[test]
    fn a_constant_is_found_under_any_visibility() {
        for src in [
            "const FRONTIER_COEFF: f64 = 6.90;\n",
            "pub const FRONTIER_COEFF: f64 = 6.90;\n",
            "pub(crate) const FRONTIER_COEFF: f64 = 6.90;\n",
            "pub(super) const FRONTIER_COEFF: f64 = 6.90;\n",
        ] {
            assert_eq!(
                extract_const_value(src, "FRONTIER_COEFF").as_deref(),
                Some("6.90"),
                "missed a visibility spelling: {src}"
            );
        }
    }

    /// A constant named only in prose must not be taken for the declaration.
    #[test]
    fn a_constant_named_in_a_comment_is_not_the_declaration() {
        let src = "/// See FRONTIER_COEFF for why.\nconst FRONTIER_COEFF: f64 = 6.90;\n";
        assert_eq!(
            extract_const_value(src, "FRONTIER_COEFF").as_deref(),
            Some("6.90")
        );
    }

    fn closure_hits(src: &str) -> Vec<String> {
        let mut f = Vec::new();
        scan_closure_predicates("probe.rs", src, &mut f);
        f
    }

    /// The tripwire must fire on the shape it exists to stop — the exact
    /// declaration a review built, which passed every test in both files while
    /// producing a bit-identical table.
    #[test]
    fn a_planted_closure_branch_is_reported() {
        let hits = closure_hits(
            "fn unmade_lateral(cap: usize, outer: usize) -> f64 {\n\
             \x20   if outer == 0 { return 0.0; }\n\
             \x20   1.0\n}\n",
        );
        assert_eq!(hits.len(), 1, "expected one finding, got {hits:?}");
        assert!(
            hits.first().is_some_and(|h| h.contains("unmade_lateral")),
            "{hits:?}"
        );
    }

    /// ...and must NOT fire on the identical line inside the function exempted
    /// by name. An exemption that is never exercised is indistinguishable from
    /// a check that does not scope by function at all.
    #[test]
    fn the_same_branch_inside_compact_bound_is_exempt() {
        let hits = closure_hits(
            "fn compact_bound(a: usize) -> f64 {\n\
             \x20   if outer == 0 { return 0.0; }\n\
             \x20   6.0\n}\n",
        );
        assert!(hits.is_empty(), "the exemption did not apply: {hits:?}");
    }

    /// The exemption must be scoped to `compact_bound` and end with it, or it
    /// silently covers whatever function follows.
    #[test]
    fn the_exemption_ends_with_the_exempt_function() {
        let hits = closure_hits(
            "fn compact_bound(a: usize) -> f64 {\n\
             \x20   if a == 0 { return 0.0; }\n\
             \x20   6.0\n}\n\
             fn frontier_notches(cap: usize, outer: usize) -> f64 {\n\
             \x20   if outer == 0 { return 0.0; }\n\
             \x20   1.0\n}\n",
        );
        assert_eq!(hits.len(), 1, "expected one finding, got {hits:?}");
        assert!(
            hits.first().is_some_and(|h| h.contains("frontier_notches")),
            "{hits:?}"
        );
    }

    /// Test code legitimately enumerates `outer` and must not trip the check —
    /// `the_clamp_fires_only_at_a_lone_outer_site` asserts on `outer == 1`.
    #[test]
    fn cfg_test_regions_are_skipped() {
        let hits = closure_hits(
            "#[cfg(test)]\nmod tests {\n\
             \x20   fn t() { assert_eq!(raw(cap, outer) < 0.0, outer == 0); }\n}\n",
        );
        assert!(hits.is_empty(), "test region was scanned: {hits:?}");
    }

    /// Step 8d's extractor has to normalise away exactly the differences that
    /// are legitimate between a library and a single-file probe, and nothing
    /// else.
    #[test]
    fn extraction_normalises_visibility_attributes_and_comments() {
        let lib =
            "#[must_use]\npub fn f(a: usize) -> f64 {\n    // a comment\n    6.0 * a as f64\n}\n";
        let probe = "fn f(a: usize) -> f64 {\n    6.0 * a as f64\n}\n";
        assert_eq!(extract_fn_body(lib, "f"), extract_fn_body(probe, "f"));
    }

    /// ...and must still see a real divergence through that normalisation.
    #[test]
    fn extraction_still_sees_a_real_divergence() {
        let lib = "pub fn f(a: usize) -> f64 { let d = g(a); if d < 1.0 { d } else { 1.0 } }\n";
        let probe = "fn f(a: usize) -> f64 { if g(a) < 1.0 { g(a) } else { 1.0 } }\n";
        assert_ne!(extract_fn_body(lib, "f"), extract_fn_body(probe, "f"));
    }

    /// A renamed function must be reported, not silently skipped.
    #[test]
    fn extraction_reports_a_missing_function() {
        assert_eq!(extract_fn_body("fn other() {}\n", "shell_size"), None);
    }

    #[expect(
        clippy::expect_used,
        reason = "CLAUDE.md: tests may unwrap freely — a probe that does not parse is a broken \
                  test, and panicking says so at the point of the mistake"
    )]
    fn hits(src: &str) -> Vec<String> {
        scan_for_derived_streams("probe.rs", src).expect("probe must parse")
    }

    /// Every shape that reintroduces a derived stream must be caught.
    ///
    /// **This corpus is the whole point of the check having a test.** Two
    /// hand-rolled versions shipped before it existed: the first caught one of
    /// these, the second six. Each was built from the list of misses the
    /// previous round happened to enumerate, so each was blind to whatever
    /// nobody probed. Every row below was found by a reviewer, not by the
    /// author — including `-> (Self, Self)`, which `clippy::use_self` actively
    /// pushes an author toward, and the trait method, which is spelled `fn`
    /// and so was invisible to a `pub fn` matcher.
    #[test]
    fn every_derived_stream_shape_is_caught() {
        const FORBIDDEN: &[(&str, &str)] = &[
            (
                "&self -> Self",
                "impl Stream { pub fn c(&self, i: u64) -> Self { todo!() } }",
            ),
            (
                "&mut self",
                "impl Stream { pub fn c(&mut self, i: u64) -> Self { todo!() } }",
            ),
            (
                "self by value",
                "impl Stream { pub fn c(self, i: u64) -> Self { todo!() } }",
            ),
            (
                "mut self",
                "impl Stream { pub fn c(mut self, i: u64) -> Self { todo!() } }",
            ),
            (
                "self: &Self",
                "impl Stream { pub fn c(self: &Self, i: u64) -> Self { todo!() } }",
            ),
            (
                "-> Stream",
                "impl Stream { pub fn c(&self, i: u64) -> Stream { todo!() } }",
            ),
            (
                "-> crate::Stream",
                "impl Stream { pub fn c(&self) -> crate::Stream { todo!() } }",
            ),
            (
                "-> (Self, Self)",
                "impl Stream { pub fn split(&self) -> (Self, Self) { todo!() } }",
            ),
            (
                "-> Option<Self>",
                "impl Stream { pub fn c(&self) -> Option<Self> { todo!() } }",
            ),
            (
                "-> Result<Self, ()>",
                "impl Stream { pub fn c(&self) -> Result<Self, ()> { todo!() } }",
            ),
            (
                "-> [Self; 4]",
                "impl Stream { pub fn q(&self) -> [Self; 4] { todo!() } }",
            ),
            (
                "-> Vec<Self>",
                "impl Stream { pub fn c(&self) -> Vec<Self> { todo!() } }",
            ),
            (
                "-> Box<Self>",
                "impl Stream { pub fn c(&self) -> Box<Self> { todo!() } }",
            ),
            (
                "-> impl Iterator<Item = Stream>",
                "impl Stream { pub fn c(&self) -> impl Iterator<Item = Stream> { todo!() } }",
            ),
            (
                "pub(crate)",
                "impl Stream { pub(crate) fn c(&self) -> Self { todo!() } }",
            ),
            (
                "private fn",
                "impl Stream { fn c(&self) -> Self { todo!() } }",
            ),
            (
                "free fn",
                "pub fn derive_child(p: &Stream, i: u64) -> Stream { todo!() }",
            ),
            (
                "inside a module",
                "mod inner { impl Stream { pub fn c(&self) -> Self { todo!() } } }",
            ),
            (
                "wrapped signature",
                "impl Stream {\n pub fn c(\n &self,\n i: u64,\n ) -> Self { todo!() }\n}",
            ),
            (
                "generic with parens",
                "impl Stream { pub fn c<F: Fn(u64) -> u64>(&self, f: F) -> Self { todo!() } }",
            ),
            (
                "async",
                "impl Stream { pub async fn c(&self) -> Self { todo!() } }",
            ),
        ];
        for (label, src) in FORBIDDEN {
            assert_eq!(hits(src).len(), 1, "{label}: not caught\n{src}");
        }
    }

    /// The item shapes that carry a derivation, split from the signature
    /// shapes above only to stay under `clippy::too_many_lines`.
    #[test]
    fn every_derived_stream_item_shape_is_caught() {
        const FORBIDDEN: &[(&str, &str)] = &[
            (
                "free fn",
                "pub fn derive_child(p: &Stream, i: u64) -> Stream { todo!() }",
            ),
            (
                "extension trait decl",
                "pub trait Ext { fn child(&self, i: u64) -> Stream; }",
            ),
            (
                "trait impl for Stream",
                "impl Ext for Stream { fn child(&self, i: u64) -> Self { todo!() } }",
            ),
            (
                "hand-written Clone",
                "impl Clone for Stream { fn clone(&self) -> Self { todo!() } }",
            ),
            (
                "inside a module",
                "mod inner { impl Stream { pub fn c(&self) -> Self { todo!() } } }",
            ),
            (
                "pub(crate)",
                "impl Stream { pub(crate) fn c(&self) -> Self { todo!() } }",
            ),
            (
                "private fn",
                "impl Stream { fn c(&self) -> Self { todo!() } }",
            ),
            (
                "async",
                "impl Stream { pub async fn c(&self) -> Self { todo!() } }",
            ),
        ];
        for (label, src) in FORBIDDEN {
            assert_eq!(hits(src).len(), 1, "{label}: not caught\n{src}");
        }
    }

    /// Shapes found by four independent sources *after* the AST rewrite
    /// shipped: `CodeRabbit` (aliases), and three review lanes (methods on other
    /// types, macro items, items in function bodies). Every row here passed
    /// the guard when it was written.
    #[test]
    fn shapes_found_after_the_ast_rewrite_are_caught() {
        const FORBIDDEN: &[(&str, &str)] = &[
            (
                "type alias hides the return",
                "type Child = Stream;\nimpl Stream { pub fn c(&self) -> Child { todo!() } }",
            ),
            (
                "renamed import hides the return",
                "use crate::Stream as Rng;\npub fn child(p: &Rng) -> Rng { todo!() }",
            ),
            (
                "alias on both sides",
                "type Rng = Stream;\nimpl Rng { pub fn c(&self) -> Rng { todo!() } }",
            ),
            (
                "method on another type",
                "impl Nursery { pub fn child(&self, p: &Stream, i: u64) -> Stream { todo!() } }",
            ),
            (
                "method on another type, wrapped return",
                "impl Nursery { pub fn brood(&self, p: &Stream) -> Vec<Stream> { todo!() } }",
            ),
            (
                "trait impl on another type",
                "impl Spawn for Nursery { fn child(&self, p: &Stream) -> Stream { todo!() } }",
            ),
            (
                "generic type",
                "impl<T> Factory<T> { pub fn child(&self, p: &Stream) -> Stream { todo!() } }",
            ),
            (
                "item inside a function body",
                "pub fn outer() { impl Stream { pub fn c(&self) -> Self { todo!() } } }",
            ),
            (
                "free fn inside a function body",
                "pub fn outer() { pub fn c(p: &Stream) -> Stream { todo!() } }",
            ),
            // A trait *default body* is a function body, and this arm used to
            // check the signature and stop — so a derivation nested one level
            // inside produced zero findings against the real crate. `Item::Fn`
            // and `ImplItem::Fn` were both already recursed into; this is the
            // same enumeration mistake one variant further along, found by
            // CodeRabbit after three earlier rounds of repairing this guard.
            (
                "free fn inside a trait default body",
                "pub trait Ext { fn spawn(&self) { pub fn c(p: &Stream) -> Stream { todo!() } } }",
            ),
            (
                "impl block inside a trait default body",
                "pub trait Ext { fn spawn(&self) { impl Stream { pub fn c(&self) -> Self { \
                 todo!() } } } }",
            ),
        ];
        for (label, src) in FORBIDDEN {
            assert_eq!(hits(src).len(), 1, "{label}: not caught\n{src}");
        }
    }

    /// Code the AST cannot read must be reported, not skipped.
    ///
    /// The missing-directory branch already resolves the same situation the
    /// same way: silence has to mean "looked and it was clean", never "did not
    /// look".
    #[test]
    fn unreadable_code_is_reported_rather_than_skipped() {
        for (label, src) in [
            ("macro item", "define_stream_ext!();"),
            (
                "macro inside an impl",
                "impl Stream { derive_children!(); }",
            ),
            // `ImplItem::Macro` was covered and `TraitItem::Macro` was not —
            // the two sit in sibling enums and were fixed one round apart.
            (
                "macro inside a trait",
                "pub trait Ext { derive_children!(); }",
            ),
        ] {
            let h = hits(src);
            assert_eq!(h.len(), 1, "{label}: not reported\n{src}");
            assert!(
                h.first().is_some_and(|m| m.contains("cannot be checked")),
                "{label}: wrong message"
            );
        }
    }

    /// The control set, and it matters as much as the corpus above.
    ///
    /// A guard that fires on correct code gets an `#[allow]` and then a
    /// deletion, taking the guarantee with it. Two of these — `self_seed` as a
    /// parameter name, and a block comment *describing* the forbidden shape —
    /// were live false positives in the textual version, in a crate that
    /// documents its invariants at length.
    #[test]
    fn legitimate_code_is_not_flagged() {
        const ALLOWED: &[(&str, &str)] = &[
            (
                "the real constructors",
                "impl Stream {\n pub const fn new(seed: u64, domain: Domain, index: u64) -> Self { todo!() }\n pub const fn sub(seed: u64, domain: Domain, index: u64, sub: u64) -> Self { todo!() }\n}",
            ),
            (
                "a draw",
                "impl Stream { pub const fn next_u64(&mut self) -> u64 { todo!() } }",
            ),
            (
                "self_seed param",
                "impl Stream { pub fn from_parts(self_seed: u64) -> Self { todo!() } }",
            ),
            (
                "selfish param",
                "impl Stream { pub fn of(selfish: bool) -> Self { todo!() } }",
            ),
            (
                "Domain builder",
                "impl Domain { pub fn next(&self) -> Self { todo!() } }",
            ),
            (
                "other type's builder",
                "impl Key { pub fn with_round(&self, r: u64) -> Self { todo!() } }",
            ),
            (
                "fn-typed argument",
                "impl Stream { pub fn f(&mut self, make: fn(u64) -> Self) -> u64 { todo!() } }",
            ),
            (
                "block comment describing the ban",
                "/* pub fn child(&self) -> Self {} */\npub fn ok() -> u64 { 0 }",
            ),
            (
                "doc comment describing the ban",
                "/// `pub fn child(&self) -> Self` is banned.\npub fn ok() -> u64 { 0 }",
            ),
            (
                "string literal",
                "pub fn ok() -> &'static str { \"pub fn child(&self) -> Self\" }",
            ),
            (
                "derived Clone",
                "#[derive(Clone)]\npub struct Stream { k: u64 }",
            ),
            (
                "Domain::ALL",
                "impl Domain { pub const ALL: [Self; 9] = [todo!(); 9]; }",
            ),
            // Found firing wrongly by two lanes after the AST rewrite. An RNG
            // implementing Iterator<Item = u64> is entirely ordinary, and
            // `Self::Item` is a projection to something that is not a stream.
            (
                "impl Iterator for Stream",
                "impl Iterator for Stream { type Item = u64; fn next(&mut self) -> Option<Self::Item> { todo!() } }",
            ),
            (
                "associated-type projection",
                "impl Stream { pub fn k(&self) -> <Self as Keyed>::Key { todo!() } }",
            ),
            // In a trait *declaration* `Self` is the implementor, unconstrained.
            (
                "trait declaration returning Self",
                "pub trait Reset { fn reset(&self) -> Self; }",
            ),
            (
                "builder trait declaration",
                "pub trait Builder { fn with(&self, x: u64) -> Self; }",
            ),
            (
                "other type's builder returning Self",
                "impl Nursery { pub fn tuned(&self, x: u64) -> Self { todo!() } }",
            ),
            (
                "Debug impl on Stream",
                "impl Debug for Stream { fn fmt(&self, f: &mut Formatter) -> Result { todo!() } }",
            ),
        ];
        for (label, src) in ALLOWED {
            assert!(
                hits(src).is_empty(),
                "{label}: false positive\n{src}\n{:?}",
                hits(src)
            );
        }
    }

    /// Unparseable source must fail loudly rather than scan to nothing.
    #[test]
    fn a_file_that_does_not_parse_is_an_error() {
        assert!(scan_for_derived_streams("probe.rs", "fn broken( {").is_err());
    }

    use super::{BANNED_CALLS, scan_rust_source};

    fn scan(src: &str) -> Vec<String> {
        let mut failures = Vec::new();
        scan_rust_source("t.rs", src, &mut failures);
        failures
    }

    #[test]
    fn reports_a_call_in_plain_code() {
        let found = scan("fn f(x: f64) -> f64 { x.exp() }\n");
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(
            found.first().is_some_and(|f| f.contains(".exp()")),
            "{found:?}"
        );
    }

    #[test]
    fn allows_the_operations_ieee_754_specifies_exactly() {
        assert!(scan("fn f(x: f64) -> f64 { (x + 1.0).sqrt() / 2.0 }\n").is_empty());
    }

    #[test]
    fn catches_mul_add_as_well_as_the_obvious_ones() {
        assert_eq!(
            scan("fn f(a: f64) -> f64 { a.mul_add(2.0, 1.0) }\n").len(),
            1
        );
    }

    /// Both Task 2 reviewers found this hole independently: the first version
    /// of `BANNED_CALLS` held only the method spelling, so `f64::exp(x)` went
    /// past while `cargo xtask` printed "all checks passed". The point-free
    /// case is the one a naive fix misses — `.map(f64::exp)` has no trailing
    /// paren for a `"f64::exp("` pattern to match.
    #[test]
    fn catches_fully_qualified_and_point_free_spellings() {
        assert_eq!(scan("fn f(x: f64) -> f64 { f64::exp(x) }\n").len(), 1);
        let point_free =
            "fn f(v: &[f64]) -> Vec<f64> { v.iter().copied().map(f64::exp).collect() }\n";
        assert_eq!(scan(point_free).len(), 1, "{:?}", scan(point_free));
    }

    /// `libm::exp` is what every wrapper in `det_math` calls, so a pattern
    /// that caught it would make that file unscannable — which is how path
    /// exemptions get introduced.
    #[test]
    fn leaves_the_libm_free_functions_alone() {
        assert!(scan("pub fn exp(x: f64) -> f64 {\n    libm::exp(x)\n}\n").is_empty());
    }

    /// A raw string with an odd number of inner quotes flips ordinary string
    /// state on; a second one flips it back and re-balances the braces, so
    /// everything between them is skipped in silence. `borbax-render` will be
    /// full of exactly this shape.
    #[test]
    fn a_raw_string_containing_quotes_does_not_swallow_the_next_line() {
        let src = "fn open() -> &'static str { r#\"<path d=\"M0 0\"# }\n\
                   fn hidden(x: f64) -> f64 { x.exp() }\n\
                   fn close() -> &'static str { r#\"<path d=\"M1 1\"# }\n";
        let found = scan(src);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(
            found.first().is_some_and(|f| f.contains(".exp()")),
            "{found:?}"
        );
    }

    #[test]
    fn raw_strings_of_every_shape_are_consumed() {
        // Braces and banned calls inside each; none may escape, and the file
        // must still balance.
        let src = "fn f() {\n    let _a = r\"{ .exp() \";\n    let _b = r#\"} .cos()\"#;\n\
                   let _c = br##\"{{ .powf( \"##;\n}\n";
        assert!(scan(src).is_empty(), "{:?}", scan(src));
    }

    /// Rust has three raw-string prefixes and the first version of
    /// `raw_string_prefix` knew about two, so `cr#"` reproduced the exact
    /// silent-blinding bug the function exists to prevent. The original test
    /// used `r#` only, which is why it passed.
    #[test]
    fn every_raw_string_prefix_is_recognised() {
        for prefix in ["r", "br", "cr"] {
            let src = format!(
                "fn open() -> T {{ {prefix}#\"<path d=\"M0 0\"# }}\n\
                 fn hidden(x: f64) -> f64 {{ x.exp() }}\n\
                 fn close() -> T {{ {prefix}#\"<path d=\"M1 1\"# }}\n"
            );
            let found = scan(&src);
            assert_eq!(found.len(), 1, "prefix {prefix:?} -> {found:?}");
            assert!(
                found.first().is_some_and(|f| f.contains(".exp()")),
                "prefix {prefix:?} -> {found:?}"
            );
        }
    }

    /// The two ban lists are the redundancy that makes a deleted `clippy.toml`
    /// entry survivable, and redundancy only works while it agrees. Nothing
    /// else keeps them in step, and they were already one short on day one —
    /// four `clippy.toml` names were covered in `BANNED_CALLS` only by prefix
    /// containment (`f64::atan` matching `f64::atan2`), which is coverage by
    /// accident rather than by intent.
    /// `.min(` in a text scan fires on `a.len().min(b.len())`, which is
    /// `Ord::min` on a `usize` and perfectly deterministic. A check that cries
    /// wolf gets relaxed, so these keep the path form only.
    const METHOD_FORM_IS_TYPE_BLIND: &[&str] = &["max", "min"];

    /// Enforced by clippy alone, because they have legitimate library-code
    /// uses that need a per-site `#[expect]` — which the text scan cannot
    /// express, and deliberately so.
    const CLIPPY_ONLY: &[&str] = &["total_cmp"];

    #[test]
    fn the_two_ban_lists_cover_the_same_functions() {
        let clippy_toml = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .map(|p| p.join("clippy.toml"))
            .and_then(|p| std::fs::read_to_string(p).ok())
            .unwrap_or_default();
        assert!(!clippy_toml.is_empty(), "clippy.toml not found or empty");

        let mut missing = Vec::new();
        for line in clippy_toml.lines() {
            let Some(rest) = line.split_once("path = \"") else {
                continue;
            };
            let Some((path, _)) = rest.1.split_once('"') else {
                continue;
            };
            let Some((_, func)) = path.split_once("::") else {
                continue;
            };
            // Both spellings required, and the method match exact.
            //
            // This was a prefix match, so `.sinh()` "covered" `sin`, `.exp2()`
            // covered `exp` and `.ln_1p()` covered `ln`. Measured: deleting
            // *both* spellings of `exp`, `ln`, `sin` and `cos` left this test
            // green — the four functions with `det_math` wrappers and the four
            // §13.1 names. Coverage by accident, which is the thing this test
            // exists to prevent, reintroduced in the other direction.
            // The deliberate divergences, named rather
            // than hidden by a loose predicate. The method spelling `.min(` is
            // type-blind in a text scan, and `a.len().min(b.len())` is
            // `Ord::min` on a `usize` — perfectly deterministic. Banning it as
            // text would fire on every integer clamp in the workspace, and a
            // check that cries wolf gets relaxed. Clippy resolves the type and
            // has no such problem, which is the clearest single case for why it
            // is the authority and this is the second pass.
            //
            // `total_cmp` diverges for a different reason: it has legitimate
            // uses in library code, each needing a *per-site* exemption, and
            // the text scan deliberately cannot express one — an `#[expect]`
            // smuggling a call past clippy is still reported here, which is a
            // feature. So it lives in `clippy.toml` alone, where `#[expect]`
            // works and the reason travels with the call site.
            let path_exempt = CLIPPY_ONLY.contains(&func);
            let method_exempt = path_exempt || METHOD_FORM_IS_TYPE_BLIND.contains(&func);

            let has_path = BANNED_CALLS.contains(&path);
            let has_method = BANNED_CALLS
                .iter()
                .any(|b| *b == format!(".{func}()") || *b == format!(".{func}("));
            if !(has_path || path_exempt) || !(has_method || method_exempt) {
                missing.push(format!(
                    "{path} (path: {}, method: {})",
                    if has_path { "ok" } else { "MISSING" },
                    if has_method { "ok" } else { "MISSING" }
                ));
            }
        }
        assert!(
            missing.is_empty(),
            "in clippy.toml but not BANNED_CALLS: {missing:?} — the text scan is \
             the only thing that would catch one of these being deleted from \
             clippy.toml, so it has to know about them"
        );

        // The other direction, which nothing checked: `clippy.toml` is the
        // authority, so a path silently disappearing from *it* is the more
        // dangerous of the two edits.
        let mut orphaned = Vec::new();
        for banned in BANNED_CALLS {
            let Some((width, func)) = banned.split_once("::") else {
                continue;
            };
            if width != "f64" && width != "f32" {
                continue;
            }
            if !CLIPPY_ONLY.contains(&func)
                && !clippy_toml.contains(&format!("path = \"{width}::{func}\""))
            {
                orphaned.push((*banned).to_owned());
            }
        }
        assert!(
            orphaned.is_empty(),
            "in BANNED_CALLS but not clippy.toml: {orphaned:?} — clippy is the \
             authority, so a path missing there is unenforced everywhere clippy \
             can see, which is almost everywhere"
        );
    }

    /// `impl Sum for Quanta` is exactly `rayon::ParallelIterator::sum`'s
    /// bound, so this is a one-word change that compiles and passes clippy.
    /// rayon is not a dependency, so `disallowed-methods` has no path to
    /// resolve and this text scan is the only check available.
    #[test]
    fn parallel_reductions_are_reported() {
        for call in [
            "v.par_iter().sum::<Quanta>()",
            "v.into_par_iter().map(f).sum::<Quanta>()",
            "v.par_sort_by(cmp)",
            "v.par_bridge().count()",
        ] {
            let src = format!("fn f(v: Vec<Quanta>) {{ let _ = {call}; }}\n");
            let found = scan(&src);
            assert_eq!(found.len(), 1, "{call:?} -> {found:?}");
            assert!(
                found.first().is_some_and(|f| f.contains("§13.4")),
                "{found:?}"
            );
        }
    }

    /// A serial fold must not be mistaken for a parallel one.
    #[test]
    fn serial_iteration_is_left_alone() {
        assert!(scan("fn f(v: Vec<Quanta>) { let _ = v.iter().sum::<Quanta>(); }\n").is_empty());
    }

    /// The `r` in `for` must not open a raw string.
    #[test]
    fn an_r_inside_an_identifier_is_not_a_raw_string_prefix() {
        let src = "fn f(v: &[f64]) {\n    for _x in v {\n        let _ = 1.0;\n    }\n}\n\
                   fn g(y: f64) -> f64 { y.cos() }\n";
        assert_eq!(scan(src).len(), 1, "{:?}", scan(src));
    }

    #[test]
    fn skips_calls_inside_a_cfg_test_module() {
        assert!(
            scan("#[cfg(test)]\nmod tests {\n    fn t(x: f64) { let _ = x.exp(); }\n}\n")
                .is_empty()
        );
    }

    /// The false-negative direction, and the reason the brace counting has to
    /// be right: a test module that never appears to close silently exempts
    /// every line after it.
    #[test]
    fn resumes_checking_after_a_cfg_test_module_closes() {
        let src = "#[cfg(test)]\nmod tests {\n    fn t(x: f64) { let _ = x.exp(); }\n}\n\
                   fn real(y: f64) -> f64 { y.cos() }\n";
        let found = scan(src);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(
            found.first().is_some_and(|f| f.contains(".cos()")),
            "{found:?}"
        );
    }

    #[test]
    fn a_cfg_test_item_that_opens_no_block_does_not_start_a_region() {
        let src = "#[cfg(test)]\nuse std::f64;\nfn real(y: f64) -> f64 { y.cos() }\n";
        assert_eq!(scan(src).len(), 1, "{:?}", scan(src));
    }

    /// Lex a source string for G5 hits the way the real check does.
    fn g5(src: &str) -> Vec<(String, String)> {
        let mut hits = Vec::new();
        // Neither `unwrap`, `expect` nor `panic!` is available — all three are
        // denied workspace-wide and the deny reaches inside `#[cfg(test)]`. The
        // emptiness assert is not decoration: a source that failed to lex would
        // make every assertion below pass vacuously.
        let stream: proc_macro2::TokenStream = src.parse().unwrap_or_default();
        assert!(!stream.is_empty(), "test source did not lex: {src:?}");
        scan_format_tokens(stream, &mut hits);
        hits
    }

    /// **The parser that defeated the predecessor.** A case-sensitive substring
    /// scan reported zero hits on this while `clippy -D warnings` stayed clean.
    #[test]
    fn a_real_format_parser_is_caught_however_it_is_spelled() {
        let src = r#"
            pub struct SmilesParser { depth: usize }
            impl SmilesParser {
                pub fn parse_smiles(&mut self, s: &str) -> usize { s.len() }
                pub fn to_smiles(&self) -> String { String::new() }
            }
            pub const IMPORT_EXTENSIONS: &[&str] =
                &["smi", "smiles", "inchi", "sdf", "pdb", "fasta"];
            pub fn read_molfile(t: &str) -> usize { t.len() }
            pub fn write_fasta() -> String { String::new() }
        "#;
        let found = g5(src);
        assert!(found.len() >= 10, "only {} hits: {found:?}", found.len());
        for want in ["smiles", "smi", "inchi", "sdf", "pdb", "fasta", "molfile"] {
            assert!(
                found.iter().any(|(seg, _)| seg == want),
                "missed {want}: {found:?}"
            );
        }
    }

    /// Doctrine must be writable where it matters. The predecessor failed on a
    /// sentence lifted from CLAUDE.md, and a gate that fires on its own
    /// rationale is one somebody disables.
    #[test]
    fn every_comment_form_may_state_the_prohibition() {
        for src in [
            "/// Never import SMILES, InChI, MOL format or FASTA.\nfn f() {}\n",
            "//! No FASTA, no InChI, no SMILES anywhere.\nfn f() {}\n",
            "/* Borbax will never read SMILES or PDB format files. */\nfn f() {}\n",
            "fn f() {}\n// a SMILES parser does not live here\n",
            "pub const Q: u8 = 1; // no SMILES parser here\n",
        ] {
            assert!(g5(src).is_empty(), "fired on prose: {src:?}");
        }
    }

    /// `mol` is deliberately out of the vocabulary: `fn embed(mol: &Molecule)`
    /// appears six times in `experiments/`, and `sdf`/`molfile` cover the same
    /// import path without the collision.
    #[test]
    fn ordinary_molecule_vocabulary_does_not_fire() {
        for src in [
            "fn embed(mol: &Molecule) -> usize { 0 }\n",
            "fn f() { let smith = 1; let summary = 2; }\n",
            "struct Molecule { atoms: usize }\n",
        ] {
            assert!(g5(src).is_empty(), "false positive: {src:?}");
        }
    }

    /// Segment splitting is what makes `smi` safe to carry: it matches
    /// `parse_smi` and not `smith`.
    /// **The acronym defeat.** A review built a working six-format importer
    /// using only all-caps type names and scored **zero** — which made this
    /// check *weaker* than the substring scan it replaced, since that one
    /// caught `SMILES` inside `SMILESParser` and matched `InChI` exactly.
    #[test]
    fn acronym_cased_identifiers_are_caught() {
        for src in [
            "pub struct SMILESParser { d: usize }",
            "pub struct PDBReader { d: usize }",
            "pub struct InChIString(String);",
            "pub struct FASTAWriter;",
            "pub struct MMCIFParser;",
            "pub fn read_pdb2_file(t: &str) -> usize { t.len() }",
            r#"pub const K: &str = "InChI";"#,
        ] {
            assert!(!g5(src).is_empty(), "escaped: {src}");
        }
    }

    /// `trim_matches` stripped `r` from the literal's *content*, not just its
    /// prefix, so these scored a false hit on `pdb`.
    #[test]
    fn a_literal_ending_in_r_is_not_a_format() {
        for src in [
            r#"pub const A: &str = "pdbr";"#,
            r#"pub const B: &str = "rpdbr";"#,
        ] {
            assert!(g5(src).is_empty(), "false positive: {src}");
        }
    }

    /// Literal prefixes and quotes must not glue to the segment. `b"pdb"` once
    /// scored zero because the `b` stayed attached.
    #[test]
    fn literal_prefixes_and_quotes_are_segment_boundaries() {
        for src in [
            r#"pub const X: &[u8] = b"pdb";"#,
            r##"pub const X: &str = r#"smiles"#;"##,
            r#"pub const D: &str = include_str!("data.inchi");"#,
        ] {
            assert!(!g5(src).is_empty(), "escaped the scan: {src}");
        }
    }

    #[test]
    fn identifiers_split_on_underscores_and_camel_case() {
        assert_eq!(identifier_segments("SmilesParser"), ["smiles", "parser"]);
        assert_eq!(
            identifier_segments("read_pdb_file"),
            ["read", "pdb", "file"]
        );
        assert_eq!(identifier_segments("SMILES"), ["smiles"]);
        assert_eq!(identifier_segments("smith"), ["smith"]);
    }

    #[test]
    fn skips_calls_in_line_and_doc_comments() {
        assert!(
            scan("/// Routes `.exp()` through det_math.\n// and .powf(2.0) here\nfn f() {}\n")
                .is_empty()
        );
    }

    #[test]
    fn skips_calls_and_braces_inside_string_literals() {
        // The unmatched `{` in the literal would leave the file looking
        // unbalanced if strings were not stripped.
        let src = "fn f(y: f64) -> f64 {\n    let _s = \"x.exp() {\";\n    y.cos()\n}\n";
        let found = scan(src);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(
            found.first().is_some_and(|f| f.contains(".cos()")),
            "{found:?}"
        );
    }

    #[test]
    fn skips_a_block_comment_spanning_lines() {
        let src =
            "/* opens a brace {\n   and mentions .exp()\n*/\nfn f(y: f64) -> f64 { y.cos() }\n";
        let found = scan(src);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(
            found.first().is_some_and(|f| f.contains(".cos()")),
            "{found:?}"
        );
    }

    #[test]
    fn a_brace_in_a_character_literal_does_not_unbalance_the_file() {
        assert!(scan("fn f() {\n    let _open = '{';\n    let _close = '}';\n}\n").is_empty());
    }

    /// `&'a str` has an apostrophe that never closes. Treating it as a
    /// character literal would swallow the rest of the line — and with it any
    /// banned call sitting on that line.
    #[test]
    fn a_lifetime_is_not_a_character_literal() {
        let src = "fn f<'a>(s: &'a str, y: f64) -> (&'a str, f64) { (s, y.cos()) }\n";
        let found = scan(src);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(
            found.first().is_some_and(|f| f.contains(".cos()")),
            "{found:?}"
        );
    }

    #[test]
    fn an_escaped_character_literal_is_consumed_whole() {
        assert!(scan("fn f() {\n    let _nl = '\\n';\n    let _q = '\\'';\n}\n").is_empty());
    }

    /// If the scanner loses track, it says so. Reporting "all checks passed"
    /// for a file it did not actually read is the one outcome that must not
    /// be possible.
    #[test]
    fn a_file_it_cannot_parse_is_reported_rather_than_skipped() {
        for broken in [
            "fn f() {\n",
            "fn f() { }\n}\n",
            "let _s = \"never closed;\n",
        ] {
            let found = scan(broken);
            assert!(
                found.iter().any(|f| f.contains("could not reliably scan")),
                "{broken:?} -> {found:?}"
            );
        }
    }
}
