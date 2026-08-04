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

/// Directories under which **every** crate must contain no data files (G1).
///
/// **This was a hand-kept list of four crate paths and it had already gone
/// stale.** `crates/borbax-ui` was created without being added to it, so
/// `crates/borbax-ui/assets/elements.json` passed `cargo xtask` cleanly — in the
/// crate with the largest G2/G4 surface in the project, whose entire job is
/// putting generated names on a screen and which therefore has a standing motive
/// to acquire a lookup table. `borbax-units` and `borbax-rng` were never on it
/// either.
///
/// So the list now names **roots, not crates**, and applies Task 10's rule for
/// exactly this failure: **invert the default — unknown ⇒ covered.** A crate
/// added under `crates/` is guarded on the day it is created, by nobody
/// remembering anything. The predecessor's verdict was computed from the paths
/// it knew, so a path it did not know was a silent pass; this one's verdict is
/// "is this file under a scanned root", which is total.
///
/// `experiments` is a root of its own even though it is not part of the
/// simulation: it implements a working version of Tasks 8-9, so a molecule set
/// or a signature table smuggled in there is the same breach by the same route.
/// Measurement *output* is not exempt either — it belongs in `docs/` or gets
/// regenerated.
///
/// `xtask` is deliberately absent, for the same reason it is absent from
/// [`TRANSCENDENTAL_SCAN_ROOTS`]: it is a build tool that cannot reach a
/// simulation result, and `.pre-commit-config.yaml` and friends would have to be
/// exempted one by one.
///
/// **Checked against the tree before this widened**: no `.csv`, `.tsv`, `.json`,
/// `.yaml`, `.yml`, `.parquet`, `.bin` or `.dat` exists anywhere under `crates/`
/// or `experiments/` today, so this costs nothing now. It is also safe for Task
/// 19: `borbax-render`'s SVG goldens are `.svg` and `insta`'s snapshots are
/// `.snap`, and neither extension is in [`DATA_EXTENSIONS`] — a golden is
/// generated output that a reader can verify by eye, not smuggled reference
/// data.
const DATA_FREE_ROOTS: &[&str] = &["crates", "experiments"];

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
    segments_preserving_case(name)
        .iter()
        .map(|seg| seg.to_ascii_lowercase())
        .collect()
}

/// [`identifier_segments`] without the case fold.
///
/// **The G2 symbol tier needs the original case and the name tier must not
/// have it**, so the boundary rules live here once and the fold is applied on
/// top. Writing a second splitter instead would let the two drift, and the
/// boundary logic is the part that has already been fixed twice — once for
/// `b"pdb"` keeping its prefix glued on, once for `SMILESParser` absorbing the
/// word after the acronym.
///
/// Case matters for symbols and not for names because `Fe` is iron while `fe`,
/// `nO`, `aT` and `iN` are ordinary fragments of source code — measured, a
/// case-insensitive symbol rule gives **442** false positives on this tree
/// against **0** for the case-sensitive one. Case carries no meaning in a name:
/// `Carbon`, `carbon` and `CARBON` all read as carbon.
fn segments_preserving_case(name: &str) -> Vec<String> {
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
        cur.push(*ch);
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

/// §13.1 — hashers whose output is not pinned across releases or across runs.
///
/// **`clippy.toml`'s `disallowed-types` list is the authority; this is a second
/// pass**, on exactly the terms [`BANNED_CALLS`] states: this runs in the
/// sub-second pre-commit hook where clippy does not, and it reads files as
/// text, so code behind an inactive `cfg` — which clippy never builds and
/// therefore never lints — is still seen. `clippy.toml` carries the criterion
/// and the argument; this carries the spellings.
///
/// **Bare identifiers, not paths, and that is the opposite choice from
/// [`BANNED_CALLS`].** There the method form (`.exp()`) and the path form
/// (`f64::exp`) are both needed because `exp` is a name that belongs to a
/// *type*, so the bare word is meaningless on its own. A hasher is imported and
/// then written unqualified — `use std::hash::DefaultHasher;` followed by
/// `DefaultHasher::new()` — so the bare identifier is the spelling that
/// actually appears, and it is distinctive enough to match on: none of these
/// three names belongs to anything else in this workspace, and a path form
/// would miss the imported case entirely.
///
/// The class this decides, stated so the next person can see what is outside
/// it: **these three identifiers as whole tokens, in comment-stripped and
/// literal-stripped code, outside `#[cfg(test)]`**.
///
/// Whole tokens rather than substrings, which was the machine reviewer's finding and is
/// half right in a way worth recording. `DefaultHasherMetrics` did fire —
/// measured — and that is a real false positive, and a check that cries wolf
/// gets relaxed. Its second example, `make_random_state`, does **not** fire:
/// the case differs.
///
/// **"Costs nothing real" was wrong, and `FxRandomState` above is why.** Under
/// substring matching, `RandomState` caught `rustc_hash::FxRandomState`;
/// tokenising lost it. That is not a curiosity — `rustc-hash` is the crate
/// CLAUDE.md names as *the* approved deterministic hasher, and its own source
/// says of `FxRandomState`: "This mirrors what
/// `std::collections::hash_map::RandomState` does". Randomly seeded per
/// process, in the crate someone reaching for determinism would add. It is
/// behind that crate's non-default `rand` feature today, so it is a forward
/// hazard rather than a live one — which is exactly when a guard is free.
/// Verified in the vendored source; two review lanes disagreed about whether
/// the type existed and one was wrong.
///
/// `std::hash::SipHasher13` really is `#[unstable]` and unnameable on stable —
/// but `siphasher::sip::SipHasher13` and `SipHasher24` are exported on stable
/// by a mainstream crate. They are outside this class deliberately: keyed
/// explicitly they are deterministic, so banning the identifier would cry
/// wolf. A third-party hasher arriving in a result-affecting path is a
/// dependency review, which is where it belongs.
///
/// Outside the class: a renaming import's *use site* (`H::new()` after
/// `use ... as H;` — though the `use` line itself is caught, because it spells
/// the identifier), a macro expansion, and a hasher from a future third-party
/// crate. An earlier version of this paragraph drew a distinction between
/// `type H = ..` and `use .. as H` that does not exist — both are caught where
/// they are written and both escape at the use site. A determinism reviewer
/// measured it. Clippy's path resolution covers the first two; the third is a
/// dependency review, which is where it belongs.
///
/// Test code is exempt for the same reason it is exempt from the transcendental
/// scan — no more than that. An earlier version claimed the exemption was
/// "load-bearing rather than incidental" because proving `seed_from_phrase`
/// disagrees with a `DefaultHasher` would mean writing one down. **No such test
/// exists**, and the sentence was a justification invented for an exemption
/// that was simply inherited. What the exemption *is* load-bearing for is the
/// `#[expect(clippy::disallowed_types)]` liveness anchor in
/// `borbax-universe`'s phrase tests, which has to name a banned type in order
/// to prove clippy still resolves it.
const BANNED_TYPES: &[&str] = &["DefaultHasher", "RandomState", "SipHasher", "FxRandomState"];

/// Directories scanned for §13.1 violations, relative to the workspace root.
///
/// `experiments` is here for the same reason it is in [`DATA_FREE_ROOTS`]: it
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
    check_signature_surface_is_pinned(root, &mut failures)?;
    check_no_closure_predicate_branch(root, &mut failures)?;
    check_packing_matches_probe(root, &mut failures)?;
    check_every_member_inherits_the_lints(root, &mut failures)?;
    check_no_real_chemistry_in_literals(root, &mut failures)?;
    check_the_viewer_stays_a_leaf(root, &mut failures)?;
    check_the_engine_stays_in_the_viewer(root, &mut failures)?;
    check_the_viewer_seam_holds(root, &mut failures)?;
    check_wall_clock_has_one_home(root, &mut failures)?;
    check_libm_has_one_home(root, &mut failures)?;

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
/// **Generalised for a second consumer that no longer exists.** The signature
/// guard it was split out for was deleted in the same task for enumerating AST
/// shapes, so there is one caller again. Kept parameterised rather than folded
/// back: the next textual guard will want it, and the fold would be a diff with
/// no behaviour change.
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

/// What a signature-facing item can do, stated by a human and checked where a
/// structural check is possible.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum SigFacing {
    /// Can see at most one signature, so no comparison is expressible in it.
    /// Refutable: if the scan counts two, the claim fails.
    Single,
    /// Can see two or more **and** takes a `Geodesic`. This is the only
    /// category in which a comparison may be written. Checked directly.
    Minimised,
    /// Hands out the stored arrays. Sees one signature, so it cannot compare —
    /// but its *return value* can be subtracted by anyone. This is the residual
    /// the guard cannot close, and it is a category rather than a `Single` so
    /// that the residual is enumerated rather than merely admitted in prose.
    Accessor,
    /// Sees two with no `Geodesic` because it **is a step of** the minimisation
    /// — `lex_cmp` is what `canonicalise` uses to pick the representative over
    /// the 60 rotations the caller supplies. Bounded: must be private, so it
    /// cannot become a consumer-facing distance without changing category.
    Primitive,
    /// A comparison on the **stored** form for identity or ordering, not for
    /// shape distance — species interning, a `BTreeMap` key, `PartialEq`.
    ///
    /// **This category is why the guard is not the one that got deleted.**
    /// §8.2's own module doc argues `impl PartialOrd for Signature` must be
    /// *listed rather than forbidden*, because a determinism-safe intern table
    /// is a `BTreeMap` and that needs `Ord`. Without `Storage` every other
    /// category refutes it — measured, all five — so the guard forbade what the
    /// spec requires, which is precisely the structural pressure that made
    /// evading the predecessor worthwhile.
    ///
    /// Bounded: the return type must be an ordering or a `bool`. A shape
    /// *distance* returns a float, so it cannot be filed here.
    Storage,
    /// Sees two with no `Geodesic` because it exists to be **measured against**
    /// the minimised form. `plain_distance` is the only one: the property
    /// "group-minimised distance is invariant under a group element and the
    /// plain one is not" cannot be stated without both halves. Bounded: must be
    /// inside `#[cfg(test)]`, so library code can never be waved through here.
    Contrast,
}

/// Every item in the workspace that can see a `borbax_molecule::Signature`,
/// pinned.
///
/// Adding a signature-facing item anywhere fails the build until it is added
/// here with a category. That friction is the mechanism, not a side effect —
/// the plan's acceptance criterion names "adding a new signature-comparison
/// consumer" as the thing to stop.
const SIGNATURE_SURFACE: &[(&str, &str, SigFacing)] = &[
    // --- crates/borbax-molecule/src/signature.rs — §8.2 itself ---------------
    //
    // The derive list is *in the entry*, so adding `PartialOrd` — which would
    // compare the stored lex-min arrays element by element, written by the
    // compiler with no body to review — reads as an unlisted item rather than
    // a silent change to a listed one.
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature: derive(Clone, Debug, PartialEq)",
        SigFacing::Storage,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::zeroed",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::extents",
        SigFacing::Accessor,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::characters",
        SigFacing::Accessor,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::permuted",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::canonicalise",
        SigFacing::Single,
    ),
    // The lex-min ordering `canonicalise` minimises *with*. Private, which is
    // what stops it becoming a consumer-facing distance under another name.
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::lex_cmp",
        SigFacing::Primitive,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::group_distance",
        SigFacing::Minimised,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "signature",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "Signature::from_parts_for_test",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/signature.rs",
        "sig",
        SigFacing::Single,
    ),
    // The ungrouped half of the invariance contrast. `#[cfg(test)]`-bounded:
    // the property "the group-minimised distance is invariant under a group
    // element and the plain one is not" cannot be stated without both halves,
    // and deleting this would make the public API's omission unfalsifiable.
    (
        "crates/borbax-molecule/src/signature.rs",
        "plain_distance",
        SigFacing::Contrast,
    ),
    // --- crates/borbax-molecule/src/binding.rs — §8.3's kernel ---------------
    (
        "crates/borbax-molecule/src/binding.rs",
        "mean_extent",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "mean_character",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "Signature::summary",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "fit",
        SigFacing::Minimised,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "affinity",
        SigFacing::Minimised,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "sig",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "from_parts",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "complement_level",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "complement_through_anti",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "complement_same_index",
        SigFacing::Single,
    ),
    (
        "crates/borbax-molecule/src/binding.rs",
        "asymmetric_fixture",
        SigFacing::Single,
    ),
];

/// The methods that hand out a signature's stored arrays.
///
/// **This is the second dimension, and it is what makes the residual small
/// enough to state.** Outside `signature.rs` the fields are private, so *every*
/// ungrouped comparison must read one of these — which means an item that calls
/// one is a comparison site regardless of what its own declaration names. That
/// closes the shapes a declaration-only scan cannot see: a method on a type
/// that merely *holds* a signature (`impl Ord for SpeciesKey`), a closure
/// inside an unlisted function, and a function that constructs both operands
/// internally. All three were measured passing green before this existed.
///
/// Matched by method name, so a different type's `extents()` also lists. That
/// is a false positive, which is the safe direction: it forces a listing, it
/// never grants a pass.
const RAW_ACCESSORS: &[&str] = &["extents", "characters"];

/// Type constructors that can hold arbitrarily many values, so one mention of
/// `Signature` inside them still means "can see two".
const PLURAL_TYPES: &[&str] = &[
    "Vec",
    "VecDeque",
    "Iterator",
    "IntoIterator",
    "BTreeMap",
    "BTreeSet",
    "HashMap",
    "HashSet",
];

/// §8.2 — no signature comparison may skip the group minimisation.
///
/// **This replaces a guard that was deleted, and the difference is the whole
/// point.** Task 9 shipped a check that read `signature.rs`'s public surface and
/// required any `pub fn` mentioning `Signature` twice to also take a
/// `Geodesic`. Review found four independent escapes — an out-of-line `mod`,
/// an `impl` in any other file of the crate, an `impl` nested in a function
/// body, and `Self` counted per-argument where a named `Signature` was counted
/// per-mention — and two classes of false positive, one of which was a
/// `Geodesic` *receiver*, i.e. `g.affinity(&a, &b)`.
///
/// It was deleted because it **enumerated AST shapes**: its verdict was
/// pass/fail computed from the shapes it knew, so a shape it did not know was a
/// silent pass. The repair rule this project applies is that you must be able to
/// name the class a guard decides and what sits outside it.
///
/// **The class this decides**, in two dimensions, because one is not enough:
///
/// 1. **Declaration.** Any item whose own declaration *names* a signature type
///    — or a workspace-wide alias of one — in a parameter, a return type, a
///    generic bound, a struct or enum field, a `const`/`static` type, or a
///    `derive`/`cfg_attr` list.
/// 2. **Body.** Any *shipped* item whose body reads a signature's stored arrays
///    through [`RAW_ACCESSORS`], whatever its declaration says.
///
/// Every such item must appear in [`SIGNATURE_SURFACE`] with a category, and
/// every entry there must still exist.
///
/// **The second dimension is what makes the residual small enough to state.**
/// Outside `signature.rs` the fields are private, so every ungrouped comparison
/// must read an accessor — which catches the shapes a declaration-only scan
/// cannot see, and all three were measured passing green with only dimension 1:
/// `impl Ord for SpeciesKey` where `SpeciesKey` merely *holds* a signature (the
/// intern-table shape §8.2's own doc invites), a closure inside an unlisted
/// function, and an item declared inside a block *expression*.
///
/// The pinned-list decision is **total**: an item either appears or it does not,
/// and no shape analysis decides that. Shape analysis decides only *which
/// category* an item may claim, and there it can only **refute** — every
/// category check is a rejection. An unreadable construct is reported, exactly
/// as [`check_no_stream_deriving_method`] reports one, because silence must mean
/// "looked and it was clean".
///
/// **`Minimised` is the one place shape analysis grants rather than refutes, and
/// it is a smoke alarm rather than a proof.** It passes when some argument type
/// or the receiver names a `Geodesic` — which `PhantomData<Geodesic<D>>`,
/// `Option<&Geodesic<D>>` and `[Geodesic<D>; 0]` all satisfy while carrying no
/// rotation group. Tightening that would be shape enumeration again; naming it
/// is the honest option.
///
/// Each of Task 9's four escapes is closed by construction rather than by a
/// rule: the scan walks **every** `.rs` file under `crates/` and `experiments/`
/// (out-of-line `mod`, `impl` in another file), recurses through
/// [`block_items`] (`impl` in a function body), and counts *values a function
/// can see* rather than mentions of a name, with `Self` resolved from the
/// enclosing `impl` (the counting bug). Both false positives are gone: a
/// `Geodesic` receiver counts, and `impl PartialOrd for Signature` is *listed*
/// — as [`SigFacing::Storage`] — rather than forbidden, because §8.2's own
/// consumers need an ordering for a `BTreeMap` intern table and a guard that
/// fires on correct code gets deleted.
///
/// **What sits outside it, named rather than implied:**
///
/// - **The bodies of the listed items.** A `Minimised` entry that takes a
///   `Geodesic` and then forgets to use it is not caught, and an `Accessor`'s
///   returned array can be subtracted by any arithmetic anywhere.
/// - **Test bodies.** Dimension 2 is shipped-code-only: every test writing
///   `a.extents()[i]` in an assertion reads an accessor, and listing all 15 of
///   them is the noise that gets a guard switched off. A test whose
///   *declaration* names a signature is still counted — that is how
///   `plain_distance` is listed.
/// - **Cross-file reachability through a type this guard cannot resolve.** A
///   holder declared in one file and compared in another is caught by dimension
///   2 only because the comparison must read an accessor. If `Signature` ever
///   grows a third way out, it belongs in [`RAW_ACCESSORS`] the same day.
/// - **Anything a macro expands to.** Reported, never passed.
///
/// A green build is evidence that no *new* comparison site appeared unreviewed.
/// It is **not** evidence that every listed body minimises.
fn check_signature_surface_is_pinned(
    root: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    let dir = root.join("crates");
    if !dir.is_dir() {
        // Same reasoning as `check_no_stream_deriving_method`: a check that
        // cannot run is itself the failure, never a silent `Ok`.
        failures.push(format!(
            "§8.2: {} is missing — the signature-surface check cannot run",
            dir.display()
        ));
        return Ok(());
    }
    let mut sources: Vec<(String, String)> = Vec::new();
    for root_dir in [&dir, &root.join("experiments")] {
        if !root_dir.is_dir() {
            continue;
        }
        for path in walk(root_dir)? {
            if path.extension().and_then(|e| e.to_str()) != Some("rs") {
                continue;
            }
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let rel = path
                .strip_prefix(root)
                .unwrap_or(&path)
                .display()
                .to_string()
                .replace('\\', "/");
            sources.push((rel, text));
        }
    }
    sources.sort();

    // **Aliases are collected workspace-wide before anything is scanned, and
    // that is a two-pass design for a reason a one-pass version cannot fix.**
    // `alias_names` is per-file, so a rename living in a *different* file left
    // every consuming file dark. Measured against the real tree: adding
    // `pub use signature::Signature as Shape;` to `lib.rs` — which already
    // re-exports `Signature` two lines up — and a new `probe.rs` containing a
    // `pub fn probe_alias(a: &Shape<12>, b: &Shape<12>) -> f64` that never
    // writes the identifier `Signature` produced **zero failures**. That is a
    // public, ungrouped, element-wise comparison in a brand-new file, and it is
    // the same one-line alias move `scan_signature_surface`'s own comment says
    // this guard defends against.
    //
    // Parsed **once** each: the alias pass and the scan pass share one AST per
    // file, so the two passes cannot disagree about what the file contains and
    // the second is not paying to re-lex it.
    let mut parsed: Vec<(&str, syn::File)> = Vec::new();
    for (rel, text) in &sources {
        let file = syn::parse_file(text).map_err(|e| format!("{rel}: {e}"))?;
        parsed.push((rel.as_str(), file));
    }
    let mut sig_names = vec!["Signature".to_owned()];
    let mut geo_names = vec!["Geodesic".to_owned()];
    for (_, file) in &parsed {
        sig_names.extend(alias_names(&file.items, "Signature"));
        geo_names.extend(alias_names(&file.items, "Geodesic"));
    }
    sig_names.sort();
    sig_names.dedup();
    geo_names.sort();
    geo_names.dedup();

    let mut found: Vec<SigItem> = Vec::new();
    for (rel, file) in &parsed {
        let (items, unreadable) = scan_parsed_surface(rel, file, &sig_names, &geo_names);
        found.extend(items);
        failures.extend(unreadable);
    }
    found.sort_by(|a, b| (&a.file, &a.name).cmp(&(&b.file, &b.name)));
    compare_signature_surface(&found, failures);
    Ok(())
}

/// One item that can see a signature.
#[derive(Debug)]
#[expect(
    clippy::struct_excessive_bools,
    reason = "the same argument `Where` makes, for the same reason: these are independent facts \
              about one item, each bounding a different category, and the lint's remedy — a \
              state enum — would claim they are mutually exclusive when `derive` + \
              `returns_ordering` and `is_pub` + `in_test` both co-occur in the pinned list \
              today. Naming them is what stops a category's bound being checked against the \
              wrong flag."
)]
struct SigItem {
    file: String,
    name: String,
    /// How many signature *values* it can see. Two is the interesting bar; a
    /// plural container counts as two because it holds arbitrarily many.
    sees: usize,
    /// Does a `Geodesic` reach it, as a parameter or as the receiver?
    geodesic: bool,
    /// Bounds `SigFacing::Primitive` — a step of the minimisation may not be
    /// reachable from outside the module that owns the minimisation.
    is_pub: bool,
    /// Bounds `SigFacing::Contrast` — an ungrouped comparison kept for
    /// measurement may not exist in shipped code.
    in_test: bool,
    /// Does its body read a signature's stored arrays? Recorded so the
    /// failure message can say *why* an item is surface — a reviewer meeting
    /// `Species::cmp2` in the list needs to know it was the body, not the
    /// declaration.
    reads_raw: bool,
    /// Does it return an ordering or a `bool` rather than a number? Bounds
    /// `SigFacing::Storage`: an identity comparison answers "same?" or "which
    /// first?", a shape distance answers "how far?".
    returns_ordering: bool,
    /// A `#[derive(..)]` list rather than an item with a body. It has no
    /// arguments to count and no visibility of its own, so the `Primitive` and
    /// `Contrast` bounds are *unmeasurable* on it and must be refused rather
    /// than passed.
    derive: bool,
    line: usize,
}

/// Can this item be named from outside the module that declares it?
///
/// **`pub(crate)` counts, and the distinction is load-bearing.** `SigFacing::
/// Primitive` exists for a step *inside* the minimisation, and its bound is that
/// no consumer can reach it. `binding.rs` is in the same crate as
/// `signature.rs`, so `pub(crate) fn lex_cmp` is reachable by every consumer
/// that matters while `matches!(vis, Visibility::Public(_))` reads it as
/// private. Only `Inherited` — a bare `fn` — holds the bound.
const fn reachable(vis: &syn::Visibility) -> bool {
    !matches!(vis, syn::Visibility::Inherited)
}

/// Match what was found against [`SIGNATURE_SURFACE`], both ways.
///
/// Both directions matter. An unlisted item is a new comparison site nobody
/// reviewed; a listed item that no longer exists is a guard pinning a deleted
/// thing, which is how a guard silently becomes vacuous.
fn compare_signature_surface(found: &[SigItem], failures: &mut Vec<String>) {
    for item in found {
        let Some((_, _, facing)) = SIGNATURE_SURFACE
            .iter()
            .find(|(f, n, _)| *f == item.file && *n == item.name)
        else {
            // Which dimension made it surface. A reviewer meeting a new entry
            // needs to know whether the declaration or the body put it there —
            // they call for different scrutiny.
            // Both dimensions can be true at once, and saying only the second
            // sends a reviewer to the wrong place. Reported precisely.
            let why = match (item.sees > 0 || item.derive, item.reads_raw) {
                (true, true) => {
                    "its declaration names a Signature and its body reads the stored arrays"
                }
                (false, true) => "its body reads a signature's stored arrays",
                _ => "its declaration names a Signature",
            };
            failures.push(format!(
                "§8.2: `{}` at {}:{} is surface — {why} — and is not in SIGNATURE_SURFACE. Add it \
                 with a category — Single (sees at most one), Minimised (sees two and takes a \
                 Geodesic), Accessor (hands out the stored arrays), Primitive (a private step \
                 *inside* the minimisation), Storage (identity or ordering on the stored form, \
                 for interning — must return an ordering or a bool), or Contrast \
                 (`#[cfg(test)]` only, the ungrouped half of a comparison the minimised one is \
                 measured against). This list is where a reviewer sees a new comparison site \
                 appear.",
                item.name, item.file, item.line
            ));
            continue;
        };
        let at = format!("`{}` at {}:{}", item.name, item.file, item.line);
        // **Every category asserts something, and a derive entry can satisfy
        // none of them.** A derive list is not a function: it has no arguments
        // to count and no visibility of its own, so `sees` is 0 and `is_pub` is
        // `false` unconditionally. Pinning one as `Primitive` or `Contrast`
        // would satisfy those bounds *vacuously*, which is the shape this whole
        // guard exists to refuse.
        if item.derive && matches!(facing, SigFacing::Primitive | SigFacing::Contrast) {
            failures.push(format!(
                "§8.2: {at} is a derive list pinned as {facing:?}, and that category's bound is \
                 unmeasurable on a derive — it would pass vacuously. Derive lists are Single, \
                 Minimised or Accessor."
            ));
            continue;
        }
        match facing {
            // The property itself, checked directly.
            SigFacing::Minimised if !item.geodesic => failures.push(format!(
                "§8.2: {at} is pinned as Minimised but takes no Geodesic, so any comparison in \
                 it is on the stored lex-min form. §8.2's storage form is safe only because \
                 every consumer compares group-minimised."
            )),
            // A category refuted by structure. Not "this is wrong" — "this
            // claim is not the one the code supports".
            SigFacing::Single | SigFacing::Accessor if item.sees >= 2 => failures.push(format!(
                "§8.2: {at} is pinned as {facing:?} but can see {} signatures, so a comparison \
                 is expressible in it. Either it takes a Geodesic and is Minimised, or it \
                 should not see two.",
                item.sees
            )),
            // The two escape-hatch categories, each held shut by the one fact
            // that makes it not an escape hatch.
            //
            // **`pub(crate)` counts as public here, and must.** `Primitive`'s
            // bound is that a step *inside* the minimisation cannot be reached
            // by a consumer — and `binding.rs` is in the same crate as
            // `signature.rs`, so `pub(crate) fn lex_cmp` is reachable by every
            // consumer that matters. Only `Visibility::Inherited` — a bare
            // private `fn` — actually holds the bound.
            SigFacing::Primitive if item.is_pub => failures.push(format!(
                "§8.2: {at} is pinned as Primitive — a step *inside* the minimisation — but it \
                 is reachable outside its module. A reachable ungrouped comparison is the thing \
                 §8.2 forbids, whatever it is called, and `pub(crate)` reaches every consumer \
                 in this crate."
            )),
            // A distance returns a number. If this returns one, it is not an
            // identity comparison whatever it is called.
            SigFacing::Storage if !item.returns_ordering && !item.derive => {
                failures.push(format!(
                    "§8.2: {at} is pinned as Storage — an identity or ordering comparison on the \
                     stored form — but it does not return an ordering or a bool. A shape \
                     distance is not storage, whatever it is named."
                ));
            }
            SigFacing::Contrast if !item.in_test => failures.push(format!(
                "§8.2: {at} is pinned as Contrast — an ungrouped comparison kept so the \
                 minimised one can be measured against it — but it is not inside `#[cfg(test)]`. \
                 That category exists for measurement, not for shipping."
            )),
            // **The three "sees two" categories must actually see two.** Without
            // this, `Minimised`/`Primitive`/`Contrast` drift into labels for
            // items that cannot compare at all, and the list stops describing
            // where comparisons live — which is the only thing it is for.
            SigFacing::Minimised | SigFacing::Primitive | SigFacing::Contrast if item.sees < 2 => {
                failures.push(format!(
                    "§8.2: {at} is pinned as {facing:?}, which claims it can see two signatures, \
                     but it sees {}. It is Single or Accessor.",
                    item.sees
                ));
            }
            _ => {}
        }
    }
    check_surface_bookkeeping(found, failures);
}

/// Both directions of the pinned list, **counted rather than searched**.
///
/// Split out of [`compare_signature_surface`] for length only. `find`/`any`
/// answer "at least one", so two items sharing a `(file, name)` key would take
/// the first entry's category for both, and a duplicated list line would never
/// be noticed.
fn check_surface_bookkeeping(found: &[SigItem], failures: &mut Vec<String>) {
    bookkeeping_against(found, SIGNATURE_SURFACE, failures);
}

/// [`check_surface_bookkeeping`] against an arbitrary pinned list.
///
/// **Both directions are counted, not merely searched.** `find`/`any` answer
/// "at least one", so two items sharing a `(file, name)` key — two `impl
/// Signature` blocks in one file each declaring `fn cmp`, say — would take the
/// first entry's category for both, and a duplicated list line would never be
/// noticed.
///
/// Parameterised only so a test can drive the **duplicate-pin** branch, whose
/// `take(i).all(..)` guard exists to report a repeated key exactly once. No
/// real-tree plant can reach it: you cannot plant a duplicate list entry
/// without editing the shipped list, so without this split an inverted
/// condition there would stay green.
fn bookkeeping_against(
    found: &[SigItem],
    pinned: &[(&str, &str, SigFacing)],
    failures: &mut Vec<String>,
) {
    for (i, (file, name, _)) in pinned.iter().enumerate() {
        let copies = pinned
            .iter()
            .filter(|(f, n, _)| f == file && n == name)
            .count();
        if copies > 1
            && pinned
                .iter()
                .take(i)
                .all(|(f, n, _)| f != file || n != name)
        {
            failures.push(format!(
                "§8.2: SIGNATURE_SURFACE pins `{name}` in {file} {copies} times. Duplicate keys \
                 mean one category silently governs items the other entry was written for."
            ));
        }
        let hits = found
            .iter()
            .filter(|i| i.file == *file && i.name == *name)
            .count();
        if hits > 1 {
            failures.push(format!(
                "§8.2: {hits} items share the key `{name}` in {file}, so one pinned category \
                 governs all of them. Give them distinguishable names, or the list is claiming \
                 something it has not checked."
            ));
        }
        if hits == 0 {
            failures.push(format!(
                "§8.2: SIGNATURE_SURFACE pins `{name}` in {file}, which no longer exists. A \
                 guard that pins a deleted item is one entry closer to vacuous — delete the \
                 line, or fix the rename."
            ));
        }
    }
}

/// The predicate behind [`check_signature_surface_is_pinned`], split out so it
/// takes `(&str, &str)` and can be unit-tested against a string.
///
/// Same reasoning as [`scan_for_derived_streams`]: a check that can only be
/// exercised against the real tree — which by construction produces no failures
/// — tests nothing.
#[cfg(test)]
fn scan_signature_surface(
    rel: &str,
    text: &str,
    sig_names: &[String],
    geo_names: &[String],
) -> Result<(Vec<SigItem>, Vec<String>), String> {
    let file = syn::parse_file(text).map_err(|e| format!("{rel}: {e}"))?;
    Ok(scan_parsed_surface(rel, &file, sig_names, geo_names))
}

/// The signature-surface scan, on an already-parsed file.
///
/// The real check parses each file once and hands the AST to both the alias
/// pass and this. The `#[cfg(test)]` `scan_signature_surface` wrapper above
/// takes source text instead, so the corpus tests can be written as strings.
fn scan_parsed_surface(
    rel: &str,
    file: &syn::File,
    sig_names: &[String],
    geo_names: &[String],
) -> (Vec<SigItem>, Vec<String>) {
    // **A file that never names `Signature` is not scanned for unreadable
    // constructs**, and that narrowing is what keeps this guard usable. Without
    // it the check reports every `macro_rules!` in the scanned tree — five in
    // `borbax-units/src/lib.rs` alone: four `unit!` invocations plus the
    // `macro_rules! unit` definition itself — and a guard whose output is mostly
    // noise gets switched off, which is the recorded way this project loses a
    // guarantee. ("Five expansions" was the first count and it was wrong by a
    // category, not by a number.)
    //
    // The residual, named: a macro *defined* in a `Signature`-free file that
    // expands to a comparison. It would have to name the type through a path,
    // and its invocation site is then in a file that does name `Signature` and
    // is scanned. A macro reaching `Signature` while naming neither is outside
    // this guard.
    //
    // **As an identifier, not as text.** `borbax-units/src/lib.rs` mentions
    // `Signature.r` in a doc comment, and a raw `text.contains` reads that as
    // "this file could compare signatures" — which is how all five of its
    // `unit!` expansions got reported. Doc comments survive parsing as
    // `#[doc = ".."]` string literals, so an ident walk excludes them by
    // construction rather than by stripping.
    // One walk, not one per name. `mentions_ident` rendered the whole file to a
    // `TokenStream` on every call, so each file was re-tokenised once per
    // signature alias plus once per raw accessor — and `sig_names` grows with
    // every alias in the workspace.
    let wanted: Vec<&str> = sig_names
        .iter()
        .map(String::as_str)
        .chain(RAW_ACCESSORS.iter().copied())
        .collect();
    if !mentions_any_ident(file, &wanted) {
        return (Vec::new(), Vec::new());
    }
    // Aliases first, for the reason `scan_for_derived_streams` gives: matching
    // the bare identifier is enumeration of *type spellings*, defeated by one
    // line — `type Shape = Signature;` — that compiles clean under `-D warnings`.
    let scan = SigScan {
        rel,
        sig: sig_names,
        geo: geo_names,
    };
    let mut out = Vec::new();
    let mut unreadable = Vec::new();
    scan.walk(&file.items, None, false, &mut out, &mut unreadable);
    (out, unreadable)
}

/// The name lists for one file, so the walker's recursion stays readable.
struct SigScan<'a> {
    rel: &'a str,
    sig: &'a [String],
    geo: &'a [String],
}

/// Everything the enclosing item makes true of a function signature.
///
/// Bundled rather than passed as five positional arguments because four of the
/// five are `bool` — the shape in which an argument-order slip compiles clean
/// and silently inverts a category.
#[derive(Clone, Copy)]
#[expect(
    clippy::struct_excessive_bools,
    reason = "the lint's remedy — a state enum, or separate arguments — is the disease here. \
              These four travel together through a recursive walk, and the alternative spelling \
              is four positional `bool`s at each call site, which is exactly the shape in which \
              an argument-order slip compiles clean and silently inverts a category. They are \
              also genuinely independent flags, not a state machine: `pub` + `cfg(test)` and \
              `impl Signature` + `impl Geodesic` are both reachable combinations."
)]
struct Where<'a> {
    scope: Option<&'a str>,
    /// Inside `impl .. Signature`: the receiver is a signature and so is `Self`.
    self_sig: bool,
    /// Inside `impl .. Geodesic`: `g.affinity(&a, &b)` is a *minimised*
    /// comparison, and reading it as ungrouped was one of the predecessor's two
    /// false positives.
    self_geo: bool,
    is_pub: bool,
    in_test: bool,
}

/// Is this item behind `#[cfg(test)]`?
///
/// Matched on the attribute's tokens rather than parsed, because the shapes that
/// matter — `#[cfg(test)]`, `#[cfg(all(test, ..))]`, `#[cfg(any(test, ..))]` —
/// all contain `test` as a bare ident inside `cfg`.
///
/// **An earlier version of this comment argued the loose reading was safe in
/// every direction, and that was false in exactly one place.** It said
/// over-reading `cfg(test)` "can only widen a category the pinned list already
/// forces a human to justify by name". `#[cfg(not(test))]` is the counter-case:
/// it means *shipped-only*, the loose reading called it test-only, and
/// `Contrast` is the one category bounded away from shipped code. So
/// `#[cfg(not(test))] fn plain_distance(a: &Signature<D>, b: &Signature<D>)`
/// pinned as `Contrast` would ship in every release build with the guard green.
/// `has_test` therefore refuses `not` rather than descending into it.
///
/// This is CLAUDE.md's "probe the guard's bookkeeping" one level down: the
/// assertion had a mutation test, its stated justification did not.
fn cfg_test(attrs: &[syn::Attribute]) -> bool {
    fn has_test(stream: proc_macro2::TokenStream) -> bool {
        let mut it = stream.into_iter();
        while let Some(tt) = it.next() {
            match tt {
                proc_macro2::TokenTree::Ident(i) if i == "test" => return true,
                // **`not(..)` inverts the predicate, so its contents are the
                // opposite of what this function is being asked.** Descending
                // into it read `#[cfg(not(test))]` — shipped-only code — as
                // test-only, which is the one direction `Contrast`'s bound
                // cannot survive: it exists precisely to keep an ungrouped
                // comparison out of shipped code.
                proc_macro2::TokenTree::Ident(i) if i == "not" => {
                    it.next();
                }
                // `all(test, ..)` and `any(test, ..)` nest one level deeper, and
                // a `TokenTree::Group`'s `to_string` is `(test, ..)`, not `test`.
                proc_macro2::TokenTree::Group(g) if has_test(g.stream()) => return true,
                _ => {}
            }
        }
        false
    }
    attrs.iter().any(|a| {
        a.path().is_ident("cfg")
            && a.meta
                .require_list()
                .is_ok_and(|l| has_test(l.tokens.clone()))
    })
}

/// How many values of a type named in `names` a parameter of this type can see.
///
/// A slice, array or collection counts as two, because holding arbitrarily many
/// is exactly the shape that lets a comparison be written without either operand
/// appearing twice in the signature.
fn sees_count(ty: &syn::Type, names: &[String]) -> usize {
    let hits = idents_of(ty)
        .iter()
        .filter(|i| names.iter().any(|n| n == *i))
        .count();
    if hits == 0 {
        0
    } else if hits >= 2 || plural(ty) {
        2
    } else {
        1
    }
}

impl SigScan<'_> {
    /// Read one function signature into a [`SigItem`], if it sees a signature.
    ///
    /// `self_sig` and `self_geo` say what the enclosing `impl` makes the
    /// receiver. This is where Task 9's fourth escape lived: it counted mentions
    /// of a named `Signature` but resolved `Self` per-argument, so
    /// `fn d(&self, o: &Self)` inside `impl Signature` read as one.
    fn item(&self, sig: &syn::Signature, at: Where<'_>, reads_raw: bool) -> Option<SigItem> {
        let Where {
            scope,
            self_sig,
            self_geo,
            is_pub,
            in_test,
        } = at;
        // Inside `impl Signature`, `Self` *is* a signature. Adding it to the
        // name list is what makes `other: &Self` count, and it is deliberately
        // not added outside such an impl, where `Self` is something else.
        let mut names: Vec<String> = self.sig.to_vec();
        if self_sig {
            names.push("Self".to_owned());
        }
        let mut geo_names: Vec<String> = self.geo.to_vec();
        if self_geo {
            geo_names.push("Self".to_owned());
        }
        // **A generic bound is a type name one indirection out.**
        // `fn d<T: AsRef<Signature<D>>>(a: T, b: T)` puts two signatures in
        // front of the body with no `Signature` in any argument type — measured
        // passing green as a `pub fn` in `binding.rs`. Any parameter bounded by
        // a signature name *is* a signature for counting purposes.
        //
        // **Both spellings.** `fn d<T: AsRef<Signature<D>>>(..)` and
        // `fn d<T>(..) where T: AsRef<Signature<D>>` are the same bound, and the
        // first version of this loop read only `type_params()` — the corpus
        // test caught the `where` half immediately, which is what a corpus is
        // for.
        for p in sig.generics.type_params() {
            if bound_by(&p.bounds, &names) {
                names.push(p.ident.to_string());
            }
        }
        if let Some(w) = &sig.generics.where_clause {
            for pred in &w.predicates {
                if let syn::WherePredicate::Type(t) = pred
                    && bound_by(&t.bounds, &names)
                {
                    names.extend(idents_of(&t.bounded_ty));
                }
            }
        }
        let mut sees = 0;
        let mut geodesic = false;
        for arg in &sig.inputs {
            match arg {
                syn::FnArg::Receiver(_) => {
                    sees += usize::from(self_sig);
                    geodesic |= self_geo;
                }
                syn::FnArg::Typed(t) => {
                    sees += sees_count(&t.ty, &names);
                    geodesic |= sees_count(&t.ty, &geo_names) > 0;
                }
            }
        }
        // A function that only *returns* a signature can see none, so it cannot
        // compare — but it is still surface, and listing it is what makes a
        // later change to it visible.
        let returns = match &sig.output {
            syn::ReturnType::Type(_, ret) => sees_count(ret, &names) > 0,
            syn::ReturnType::Default => false,
        };
        let returns_ordering = match &sig.output {
            syn::ReturnType::Type(_, ret) => {
                let idents = idents_of(ret);
                idents.iter().any(|i| i == "Ordering" || i == "bool")
                    && !idents.iter().any(|i| i == "f64" || i == "f32")
            }
            syn::ReturnType::Default => false,
        };
        // **A body that reads the stored arrays is surface even when nothing in
        // the declaration names a signature.** That is the whole point of the
        // second dimension: `impl Ord for SpeciesKey { fn cmp(&self, o: &Self) }`
        // where `SpeciesKey` holds a `Signature` names no signature anywhere in
        // its own declaration, and measured green before this clause existed.
        // **The body dimension applies to shipped code only, and the reason is
        // the recorded one.** Every test that writes `a.extents()[i]` in an
        // assertion reads a raw accessor — 15 of them — and a guard whose output
        // is mostly noise gets switched off, which is how this project loses a
        // guarantee. The hazard `RAW_ACCESSORS` exists for is a *shipped*
        // ungrouped comparison.
        //
        // The residual this leaves, named: a test body that compares two
        // signatures it obtained without naming the type in its own
        // declaration. A test whose declaration *does* name it is still
        // counted, which is how `plain_distance` is listed — and `Contrast` is
        // the category for the one test that legitimately compares ungrouped.
        let body_surface = reads_raw && !in_test;
        if sees == 0 && !returns && !body_surface {
            return None;
        }
        Some(SigItem {
            file: self.rel.to_owned(),
            name: scope.map_or_else(|| sig.ident.to_string(), |s| format!("{s}::{}", sig.ident)),
            sees,
            geodesic,
            is_pub,
            in_test,
            derive: false,
            returns_ordering,
            reads_raw,
            line: syn::spanned::Spanned::span(sig).start().line,
        })
    }

    /// Walk every item, at every nesting depth.
    fn walk(
        &self,
        items: &[syn::Item],
        scope: Option<&str>,
        in_test: bool,
        out: &mut Vec<SigItem>,
        unreadable: &mut Vec<String>,
    ) {
        for item in items {
            match item {
                syn::Item::Mod(m) => {
                    if let Some((_, inner)) = &m.content {
                        self.walk(inner, scope, in_test | cfg_test(&m.attrs), out, unreadable);
                    }
                }
                syn::Item::Fn(f) => {
                    let test = in_test | cfg_test(&f.attrs);
                    out.extend(self.item(
                        &f.sig,
                        Where {
                            scope,
                            self_sig: false,
                            self_geo: false,
                            is_pub: reachable(&f.vis),
                            in_test: test,
                        },
                        reads_raw(&f.block),
                    ));
                    self.walk(&block_items(&f.block), scope, test, out, unreadable);
                }
                syn::Item::Impl(i) => self.walk_impl(i, in_test, out, unreadable),
                syn::Item::Trait(t) => self.walk_trait(t, in_test, out, unreadable),
                // **A field is a comparison site one step removed.** A struct
                // holding two signatures gives its methods both operands with
                // no `Signature` anywhere in their own signatures, which is the
                // one hole a function-only scan cannot see.
                syn::Item::Struct(s) => {
                    self.fields(&s.ident.to_string(), &s.fields, &s.attrs, in_test, out);
                }
                // **`fold`, not `any`.** `any` short-circuits, so every variant
                // after the first signature-bearing one went unvisited — measured
                // with `enum Holder { First(Signature<D>), Second(Signature<D>) }`,
                // where only `Holder::First.0` was reported. A walker that stops
                // walking is the exact failure mode this guard replaced.
                syn::Item::Enum(e) => {
                    // Seeded with the type-level check: `enum Signature { .. }`
                    // holds signatures by being one, and without this its
                    // derive list never reached the pinned name.
                    let named = self.sig.iter().any(|n| e.ident == n.as_str());
                    let holds = e.variants.iter().fold(named, |acc, v| {
                        let held = self.fields(
                            &format!("{}::{}", e.ident, v.ident),
                            &v.fields,
                            &v.attrs,
                            in_test,
                            out,
                        );
                        acc | held
                    });
                    self.derives(&e.ident.to_string(), &e.attrs, holds, in_test, out);
                }
                // A union's fields are reached exactly as a struct's are, and
                // omitting the arm would have been a silent pass.
                syn::Item::Union(u) => {
                    self.fields(
                        &u.ident.to_string(),
                        &syn::Fields::Named(u.fields.clone()),
                        &u.attrs,
                        in_test,
                        out,
                    );
                }
                // A `const` or `static` holding signatures is a comparison site
                // with no function and no field: `static PAIR: [Signature<12>; 2]`
                // puts both operands in scope of every body in the module.
                syn::Item::Const(c) => {
                    self.value(&c.ident.to_string(), &c.ty, &c.vis, in_test, out);
                }
                syn::Item::Static(c) => {
                    self.value(&c.ident.to_string(), &c.ty, &c.vis, in_test, out);
                }
                syn::Item::Macro(m) => unreadable.push(sig_unanalysable(self.rel, &m.mac.path)),
                syn::Item::Verbatim(ts) => unreadable.push(sig_unanalysable_at(
                    self.rel,
                    syn::spanned::Spanned::span(ts).start().line,
                    "verbatim item",
                )),
                // **Provably cannot hold a comparison**, so ignored by name
                // rather than by falling through. `Type` is an alias, which
                // `alias_names` resolves; `Use` and `ExternCrate` bring names
                // into scope without declaring a body; `TraitAlias` declares no
                // items.
                syn::Item::Type(_)
                | syn::Item::Use(_)
                | syn::Item::ExternCrate(_)
                | syn::Item::TraitAlias(_) => {}
                // **Everything else is reported.** `syn::Item` is
                // `#[non_exhaustive]`, so a bare `_ => {}` is a silent pass for
                // every variant nobody thought of — `Item::ForeignMod` today,
                // and whatever a `syn` upgrade adds tomorrow. That is precisely
                // the "a shape it did not know was a silent pass" failure this
                // guard replaced.
                other => unreadable.push(sig_unanalysable_at(
                    self.rel,
                    syn::spanned::Spanned::span(other).start().line,
                    "an item variant this walker does not name",
                )),
            }
        }
    }

    /// An `impl` block. The receiver's meaning comes from the self type, which
    /// is where Task 9's fourth escape lived: it counted mentions of a named
    /// `Signature` but resolved `Self` per-argument, so `fn d(&self, o: &Self)`
    /// inside `impl Signature` read as one signature rather than two.
    fn walk_impl(
        &self,
        i: &syn::ItemImpl,
        in_test: bool,
        out: &mut Vec<SigItem>,
        unreadable: &mut Vec<String>,
    ) {
        let idents = idents_of(&i.self_ty);
        let self_sig = idents.iter().any(|t| self.sig.iter().any(|n| n == t));
        let self_geo = idents.iter().any(|t| self.geo.iter().any(|n| n == t));
        let name = idents.first().cloned().unwrap_or_else(|| "?".to_owned());
        let test = in_test | cfg_test(&i.attrs);
        for it in &i.items {
            match it {
                syn::ImplItem::Fn(f) => {
                    let ftest = test | cfg_test(&f.attrs);
                    out.extend(self.item(
                        &f.sig,
                        Where {
                            scope: Some(&name),
                            self_sig,
                            self_geo,
                            // A method of a trait impl has no visibility of its
                            // own and is as public as the trait — never
                            // `Primitive`.
                            is_pub: i.trait_.is_some() || reachable(&f.vis),
                            in_test: ftest,
                        },
                        reads_raw(&f.block),
                    ));
                    self.walk(&block_items(&f.block), None, ftest, out, unreadable);
                }
                // The same hole `Item::Const` had, one nesting level in:
                // `impl Foo { const PAIR: [Signature<12>; 2] = ..; }` puts both
                // operands in scope of every method on `Foo`.
                syn::ImplItem::Const(c) => self.value(
                    &format!("{name}::{}", c.ident),
                    &c.ty,
                    &c.vis,
                    test | cfg_test(&c.attrs),
                    out,
                ),
                syn::ImplItem::Macro(m) => {
                    unreadable.push(sig_unanalysable(self.rel, &m.mac.path));
                }
                syn::ImplItem::Verbatim(ts) => unreadable.push(sig_unanalysable_at(
                    self.rel,
                    syn::spanned::Spanned::span(ts).start().line,
                    "verbatim impl item",
                )),
                // An associated type is an alias; nothing else in an `impl` is
                // ignorable by construction, so the rest is reported.
                syn::ImplItem::Type(_) => {}
                other => unreadable.push(sig_unanalysable_at(
                    self.rel,
                    syn::spanned::Spanned::span(other).start().line,
                    "an impl-item variant this walker does not name",
                )),
            }
        }
    }

    /// A trait declaration. `Self` is unknowable here, so only a literal
    /// `Signature` counts — a trait implemented *for* `Signature` is caught by
    /// [`SigScan::walk_impl`] instead.
    fn walk_trait(
        &self,
        t: &syn::ItemTrait,
        in_test: bool,
        out: &mut Vec<SigItem>,
        unreadable: &mut Vec<String>,
    ) {
        let name = t.ident.to_string();
        let test = in_test | cfg_test(&t.attrs);
        for it in &t.items {
            match it {
                syn::TraitItem::Fn(f) => {
                    out.extend(self.item(
                        &f.sig,
                        Where {
                            scope: Some(&name),
                            self_sig: false,
                            self_geo: false,
                            // A trait method is as reachable as the trait, and
                            // it is spelled `fn`, never `pub fn` — keying on the
                            // keyword is how the predecessor missed extension
                            // traits, the idiomatic way to add a method in Rust
                            // and so the most likely reinstatement rather than
                            // the least.
                            is_pub: true,
                            in_test: test,
                        },
                        f.default.as_ref().is_some_and(reads_raw),
                    ));
                    // **A default body is a function body.** It can declare
                    // items, and those items can compare signatures.
                    if let Some(block) = &f.default {
                        self.walk(&block_items(block), None, test, out, unreadable);
                    }
                }
                // A trait constant is as reachable as the trait, hence the
                // `Public` visibility — the same reasoning as its methods.
                syn::TraitItem::Const(c) => self.value(
                    &format!("{name}::{}", c.ident),
                    &c.ty,
                    &syn::Visibility::Public(syn::token::Pub::default()),
                    test | cfg_test(&c.attrs),
                    out,
                ),
                syn::TraitItem::Macro(m) => {
                    unreadable.push(sig_unanalysable(self.rel, &m.mac.path));
                }
                syn::TraitItem::Verbatim(ts) => unreadable.push(sig_unanalysable_at(
                    self.rel,
                    syn::spanned::Spanned::span(ts).start().line,
                    "verbatim trait item",
                )),
                syn::TraitItem::Type(_) => {}
                other => unreadable.push(sig_unanalysable_at(
                    self.rel,
                    syn::spanned::Spanned::span(other).start().line,
                    "a trait-item variant this walker does not name",
                )),
            }
        }
    }

    /// A `const` or `static` whose type holds signatures.
    ///
    /// It sees whatever its type holds and takes no `Geodesic`, so a
    /// `static PAIR: [Signature<12>; 2]` is a two-signature site available to
    /// every body in the module without appearing in any function signature.
    fn value(
        &self,
        name: &str,
        ty: &syn::Type,
        vis: &syn::Visibility,
        in_test: bool,
        out: &mut Vec<SigItem>,
    ) {
        let sees = sees_count(ty, self.sig);
        if sees == 0 {
            return;
        }
        out.push(SigItem {
            file: self.rel.to_owned(),
            name: name.to_owned(),
            sees,
            geodesic: false,
            is_pub: reachable(vis),
            in_test,
            derive: false,
            returns_ordering: false,
            reads_raw: false,
            line: syn::spanned::Spanned::span(ty).start().line,
        });
    }

    /// Record signature-typed fields, and the type's derive list if it is or
    /// holds a signature.
    fn fields(
        &self,
        owner: &str,
        fields: &syn::Fields,
        attrs: &[syn::Attribute],
        in_test: bool,
        out: &mut Vec<SigItem>,
    ) -> bool {
        let mut holds = self.sig.iter().any(|n| n == owner);
        for (i, f) in fields.iter().enumerate() {
            let sees = sees_count(&f.ty, self.sig);
            if sees == 0 {
                continue;
            }
            holds = true;
            let field = f
                .ident
                .as_ref()
                .map_or_else(|| i.to_string(), std::string::ToString::to_string);
            out.push(SigItem {
                file: self.rel.to_owned(),
                name: format!("{owner}.{field}"),
                sees,
                geodesic: false,
                is_pub: reachable(&f.vis),
                in_test,
                derive: false,
                returns_ordering: false,
                reads_raw: false,
                line: syn::spanned::Spanned::span(f).start().line,
            });
        }
        self.derives(owner, attrs, holds, in_test, out);
        holds
    }

    /// A derive list is a comparison site: `#[derive(PartialOrd)]` on a type
    /// holding signatures compares the stored lex-min arrays element by element,
    /// which is precisely the ungrouped comparison §8.2 forbids — written by the
    /// compiler, so no function body exists to review.
    fn derives(
        &self,
        owner: &str,
        attrs: &[syn::Attribute],
        holds: bool,
        in_test: bool,
        out: &mut Vec<SigItem>,
    ) {
        if !holds {
            return;
        }
        for attr in attrs {
            // **`#[cfg_attr(pred, derive(..))]` is a derive.** Keying on
            // `is_ident("derive")` alone made it invisible: measured, adding
            // `#[cfg_attr(not(test), derive(PartialOrd))]` to `Signature` left
            // the pinned name at `derive(Clone, Debug, PartialEq)` and the guard
            // green — a lexicographic ordering over the stored lex-min arrays
            // present in every *shipped* build, added with no diff to the list.
            //
            // The predicate goes into the name too, so a conditional ordering
            // reads as an unlisted item rather than a listed one.
            let cfg_attr = attr.path().is_ident("cfg_attr");
            if !attr.path().is_ident("derive") && !cfg_attr {
                continue;
            }
            if cfg_attr {
                let tokens = attr
                    .meta
                    .require_list()
                    .map_or_else(|_| "<unparsed>".to_owned(), |l| l.tokens.to_string());
                if tokens.contains("derive") {
                    out.push(SigItem {
                        file: self.rel.to_owned(),
                        name: format!("{owner}: cfg_attr({tokens})"),
                        sees: 0,
                        geodesic: false,
                        is_pub: false,
                        in_test,
                        derive: true,
                        returns_ordering: true,
                        reads_raw: false,
                        line: syn::spanned::Spanned::span(attr).start().line,
                    });
                }
                continue;
            }
            let mut traits: Vec<String> = Vec::new();
            // **Every segment, joined — not `get_ident()`.** `get_ident` returns
            // `None` for any multi-segment path, so `#[derive(Clone, Debug,
            // PartialEq, ::core::cmp::PartialOrd)]` produced the *identical*
            // pinned name as the plain three-trait list: a silent pass giving
            // `Signature` a lexicographic ordering over the stored lex-min
            // arrays, written by the compiler, with the guard green. Measured.
            // `rustc` accepts the qualified spelling without a warning, and
            // `unused_qualifications` is not in `[workspace.lints.rust]`.
            //
            // Single-segment paths join to themselves, so the pinned entries
            // are unaffected.
            let parsed = attr.parse_nested_meta(|meta| {
                traits.push(
                    meta.path
                        .segments
                        .iter()
                        .map(|seg| seg.ident.to_string())
                        .collect::<Vec<_>>()
                        .join("::"),
                );
                Ok(())
            });
            if parsed.is_err() {
                // Swallowing this truncates `traits` and silently *shrinks* the
                // pinned name — the same shape as the defect above, arrived at
                // from the other side.
                traits.push("<unparsed>".to_owned());
            }
            traits.sort();
            out.push(SigItem {
                file: self.rel.to_owned(),
                // The derived traits are *in the name*, so adding `PartialOrd`
                // is an unlisted item rather than a silent change to a listed
                // one.
                name: format!("{owner}: derive({})", traits.join(", ")),
                sees: 0,
                geodesic: false,
                is_pub: false,
                in_test,
                derive: true,
                returns_ordering: true,
                reads_raw: false,
                line: syn::spanned::Spanned::span(attr).start().line,
            });
        }
    }
}

/// Does any of `names` appear anywhere in the file as an identifier?
///
/// Deliberately not `text.contains(name)`: a doc comment parses to a `#[doc =
/// ".."]` string literal, which carries no `Ident`, so prose about a type never
/// makes its file look like a consumer of one.
///
/// Takes a slice rather than one name because the caller has several and the
/// file was otherwise re-tokenised once per name.
fn mentions_any_ident(file: &syn::File, names: &[&str]) -> bool {
    fn walk(stream: proc_macro2::TokenStream, names: &[&str]) -> bool {
        stream.into_iter().any(|tt| match tt {
            proc_macro2::TokenTree::Ident(i) => {
                let s = i.to_string();
                names.contains(&s.as_str())
            }
            proc_macro2::TokenTree::Group(g) => walk(g.stream(), names),
            _ => false,
        })
    }
    walk(quote::ToTokens::to_token_stream(file), names)
}

/// Is any of these bounds parameterised by one of `names`?
///
/// A bound is a type name one indirection out: `T: AsRef<Signature<D>>` puts a
/// signature in front of the body with no `Signature` in any argument type.
fn bound_by(
    bounds: &syn::punctuated::Punctuated<syn::TypeParamBound, syn::Token![+]>,
    names: &[String],
) -> bool {
    bounds.iter().any(|b| match b {
        syn::TypeParamBound::Trait(t) => t
            .path
            .segments
            .iter()
            .flat_map(|seg| match &seg.arguments {
                syn::PathArguments::AngleBracketed(a) => a
                    .args
                    .iter()
                    .filter_map(|arg| match arg {
                        syn::GenericArgument::Type(ty) => Some(idents_of(ty)),
                        _ => None,
                    })
                    .flatten()
                    .collect::<Vec<_>>(),
                _ => Vec::new(),
            })
            .any(|i| names.contains(&i)),
        _ => false,
    })
}

/// Does this body read a signature's stored arrays?
///
/// Matched on the method *name* anywhere in the block's tokens, including inside
/// closures and nested items — deliberately loose, because every miss is a
/// silent pass and every false positive is only a listing a human must make.
fn reads_raw(block: &syn::Block) -> bool {
    fn walk(stream: proc_macro2::TokenStream) -> bool {
        stream.into_iter().any(|tt| match tt {
            proc_macro2::TokenTree::Ident(i) => RAW_ACCESSORS.contains(&i.to_string().as_str()),
            proc_macro2::TokenTree::Group(g) => walk(g.stream()),
            _ => false,
        })
    }
    walk(quote::ToTokens::to_token_stream(block))
}

/// Can this type hold arbitrarily many values?
fn plural(ty: &syn::Type) -> bool {
    match ty {
        syn::Type::Slice(_) | syn::Type::Array(_) => true,
        syn::Type::Reference(r) => plural(&r.elem),
        syn::Type::Paren(p) => plural(&p.elem),
        syn::Type::Group(g) => plural(&g.elem),
        syn::Type::Ptr(p) => plural(&p.elem),
        // **Recurse into generic arguments.** `Box<[Signature<D>]>` and
        // `Rc<[Signature<D>]>` are neither `Type::Slice` at the top level nor
        // in `PLURAL_TYPES`, so both read as singular — measured, `fn
        // f(s: Box<[Signature<D>]>)` counted one. The name list alone is
        // enumeration; descending is structural.
        syn::Type::Path(p) => {
            idents_of(ty)
                .iter()
                .any(|i| PLURAL_TYPES.contains(&i.as_str()))
                || p.path.segments.iter().any(|seg| match &seg.arguments {
                    syn::PathArguments::AngleBracketed(a) => a.args.iter().any(|arg| match arg {
                        syn::GenericArgument::Type(t) => plural(t),
                        _ => false,
                    }),
                    _ => false,
                })
        }
        _ => idents_of(ty)
            .iter()
            .any(|i| PLURAL_TYPES.contains(&i.as_str())),
    }
}

fn sig_unanalysable(rel: &str, path: &syn::Path) -> String {
    let name = path
        .segments
        .last()
        .map_or_else(|| "?".to_owned(), |s| s.ident.to_string());
    let line = syn::spanned::Spanned::span(path).start().line;
    sig_unanalysable_at(rel, line, &format!("`{name}!` expansion"))
}

/// The line number is not decoration: without it two macros in one file produce
/// **byte-identical** failure strings, and on a three-OS matrix you cannot tell
/// "two sites" from "one message printed twice".
fn sig_unanalysable_at(rel: &str, line: usize, what: &str) -> String {
    format!(
        "§8.2: {what} at {rel}:{line} cannot be checked for signature comparisons — this guard \
         reads the AST and does not expand macros, so silence here would mean \"did not look\" \
         rather than \"looked and it was clean\". Write the impl out."
    )
}

/// G1 — no real chemistry data enters the repository.
///
/// Walks [`DATA_FREE_ROOTS`] whole rather than a list of crate paths, so a crate
/// created tomorrow is covered today. See that constant for why the hand-kept
/// version was replaced and what it had already missed.
fn check_no_data_files(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    for scan_root in DATA_FREE_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            // **Loud, not silent.** `crates/` not existing means the layout has
            // moved and this check looked at nothing — the state a §5 guard is
            // least able to afford, and the shape `check_blocklist_present`
            // already records being caught by.
            failures.push(format!(
                "G1: scan root {scan_root:?} does not exist, so the data-file check did not look"
            ));
            continue;
        }
        for entry in walk(&dir)? {
            // **Lowercased, because the extension list is exact and the
            // filesystem is not.** `elements.json` fires; `elements.JSON` did
            // not, and on macOS and Windows those are the *same file* to the OS
            // — so the rename that defeats a G1 guard is invisible to the person
            // making it, and `git mv` is not even required. The widened scan's
            // verdict is total over paths; this is what makes it total over
            // spellings too.
            let Some(ext) = entry
                .extension()
                .and_then(|e| e.to_str())
                .map(str::to_ascii_lowercase)
            else {
                continue;
            };
            if DATA_EXTENSIONS.contains(&ext.as_str()) {
                failures.push(format!(
                    "G1: data file under {scan_root}: {}",
                    entry.display()
                ));
            }
        }
    }
    Ok(())
}

/// Every workspace member opts in to `[workspace.lints]`.
///
/// **Nothing enforced this, and the trap it guards has already cost this
/// repository twice.** `[workspace.lints]` is *inert* for a member that omits
/// `[lints] workspace = true` — the crate does not inherit the workspace value,
/// it has no lints at all. `xtask` shipped that way, so `unsafe_code = "forbid"`
/// forbade nothing and `unwrap_used = "deny"` denied nothing, verified at the
/// time by planting an `.unwrap()` that drew no lint. The identical shape hit
/// `rust-version` the day before. Both were found by accident.
///
/// Four manifests carry a *comment* warning about this. A comment is not a
/// guard: `fmt`, `clippy`, `test`, `test --release` and `doc` all stay green
/// when the two lines are deleted, which is precisely what makes the omission
/// attractive. `crates/borbax-ui` is the crate with a standing motive to try it
/// — §13.1's transcendental ban reaches every member, and an orbit camera at
/// Step 3 of the viewer plan needs `sin` and `cos`. The honest answer there is a
/// per-site `#[expect(..., reason = "...")]`; the quiet one is dropping the
/// lints table, and this is what makes the quiet one fail.
///
/// Fail-closed on an unreadable manifest, for the reason
/// `check_no_real_chemical_formats` gives about untokenisable files: "could not
/// read" and "read it and it was fine" must not produce the same verdict.
fn check_every_member_inherits_the_lints(
    root: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    /// The floor for the bookkeeping check at the end. Seven members exist —
    /// five under `crates/`, plus `experiments` and `xtask` — and the bar sits
    /// one below so that adding a crate does not move it and removing one does
    /// not fire it spuriously.
    const MINIMUM_MANIFESTS: usize = 6;

    let mut examined = 0_usize;

    for scan_root in DATA_FREE_ROOTS.iter().chain(std::iter::once(&"xtask")) {
        let dir = root.join(scan_root);
        if !dir.exists() {
            // **Loud, not `continue`.** The first version of this silently
            // skipped a missing root, and its own `examined` bar could not
            // cover the gap: with `crates/` renamed, `experiments` and `xtask`
            // still supply two manifests, so `examined == 2`, `2 < 2` is false,
            // and five members went unchecked while this reported clean. Eight
            // other guards fired on that mutation; this one said nothing.
            //
            // Its sibling `check_no_data_files` was changed in this same commit
            // to fail loud here, with a comment explaining why silence is
            // unaffordable for a §5 guard — and then this was written twenty
            // lines below it with the opposite behaviour.
            failures.push(format!(
                "lints: scan root {scan_root:?} does not exist, so no manifest under it \
                 was checked for `[lints] workspace = true`"
            ));
            continue;
        }
        for manifest in walk(&dir)?
            .into_iter()
            .filter(|p| p.file_name().is_some_and(|n| n == "Cargo.toml"))
        {
            let rel = manifest
                .strip_prefix(root)
                .unwrap_or(&manifest)
                .display()
                .to_string();
            let Ok(src) = std::fs::read_to_string(&manifest) else {
                failures.push(format!(
                    "lints: {rel} could not be read, so its `[lints] workspace = true` \
                     was not checked"
                ));
                continue;
            };
            examined += 1;

            // **Section-scoped, because the substring version was defeated by
            // the manifests it was written to check.** The first attempt was
            //
            //     src.contains("[lints]") && src.contains("workspace = true")
            //
            // and the second conjunct is satisfied unconditionally by every
            // manifest in this workspace — `edition.workspace = true` is on
            // line 3 of all seven. So the predicate degenerated to
            // `contains("[lints]")`, which a *comment* satisfies, and which an
            // empty `[lints]` table satisfies. Cargo accepts an empty `[lints]`
            // table and the crate then inherits nothing: three reviewers
            // independently planted `[lints]` with the opt-in line deleted plus
            // an `.unwrap()` in shipped library code, and got a clean gate and a
            // clean `clippy -D warnings`.
            //
            // That is the precise defect this function exists to prevent,
            // reachable with this function green — and the doc comment above
            // said "deliberately strict about the spelling" while the code was
            // not strict at all. It was probed against deleting the whole table,
            // which it does catch; it was not probed against deleting one line
            // of it.
            //
            // `[lints] workspace = true` and `lints.workspace = true` are the
            // two forms cargo accepts; a `[lints.clippy]` override table is
            // *rejected by cargo itself* alongside `workspace = true`, so there
            // is no third shape to admit. Trimmed whole-line equality also
            // rejects both commented forms, since a comment line starts `#`.
            let has_table = src
                .lines()
                .map(str::trim)
                // `split('#')` so `[lints]  # inherit workspace lints` is
                // recognised. Fail-closed either way, but a valid manifest
                // reported as non-compliant is a confusing failure, and a guard
                // that cries wolf is a guard that gets deleted.
                .skip_while(|line| line.split('#').next().map(str::trim) != Some("[lints]"))
                .skip(1)
                .take_while(|line| !line.starts_with('['))
                .any(|line| line == "workspace = true");
            let has_dotted = src
                .lines()
                .map(str::trim)
                .any(|line| line == "lints.workspace = true");
            if !has_table && !has_dotted {
                failures.push(format!(
                    "lints: {rel} does not opt in to `[workspace.lints]` — without \
                     `[lints] workspace = true` the crate has NO lints rather than the \
                     workspace's, so `unsafe_code = \"forbid\"` and `unwrap_used = \"deny\"` \
                     are silently inert in it"
                ));
            }
        }
    }

    // The bookkeeping check, which is the half that does not test itself: this
    // counts **manifests read**, not crates or directories.
    //
    // **The first version of this comment was wrong in both halves and the bar
    // was decorative.** It said "six members exist"; `cargo metadata --no-deps`
    // reports **seven** — five under `crates/`, plus `experiments` and `xtask`.
    // And it justified the bar as "deleting the filter does [move it]", which
    // runs backwards: deleting the `Cargo.toml` filter makes `examined` go *up*,
    // reading every file in three trees. The bar could not fire on the mutation
    // its own comment credited it with catching.
    //
    // The missing-root case that *would* have exercised it is now caught above,
    // loudly and by name, which is the better place for it. What is left here is
    // a floor against a walk that matches nothing at all, set one below the
    // seven that exist so that adding a crate does not move it and removing one
    // does not fire it spuriously. The constant is declared at the top of this
    // function because `clippy::items_after_statements` denies it here.
    if examined < MINIMUM_MANIFESTS {
        failures.push(format!(
            "lints: only {examined} manifest(s) were read and at least \
             {MINIMUM_MANIFESTS} were expected, so this check did not look at the \
             whole workspace"
        ));
    }
    Ok(())
}

/// Source with every whole-line comment removed.
///
/// The seam checks below are about **imports and calls**, not about prose, and
/// the two are easy to confuse because the prose necessarily quotes the thing it
/// forbids: `panel.rs`'s module doc says "there is no `format!` in this file",
/// and `lib.rs`'s doc names the engine while explaining which files may.
/// A raw `grep` therefore reports more files naming it, and one `format!`,
/// in `panel.rs` — which is why an earlier version of CLAUDE.md's "check it by
/// grep" instruction was false as literally written.
fn code_only(src: &str) -> String {
    src.lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// `borbax-ui` must stay a leaf, and its per-file seam must hold.
///
/// **Both are load-bearing claims that were stated in four documents and
/// enforced in none**, which is the argument
/// [`check_every_member_inherits_the_lints`] makes about its own subject one
/// screen above: a comment is not a guard.
///
/// **Leaf-ness licenses everything else the crate is allowed to do.** `f32`,
/// `HashMap`, a wall-clock read and unordered iteration are all permitted in
/// there *because* it is a pure consumer that no result can flow out of. The day
/// something depends on it, the engine's ~500-package tree becomes reachable from
/// a result path and every one of those relaxations turns into a §13.1 breach —
/// with no diff line saying so.
///
/// **The seam is what makes the crate testable at all.** `egui` lays out on the
/// CPU and cannot open a window, so `egui_kittest` can drive `panel::draw`
/// headlessly; the engine is the shell that can, so it is confined to the files
/// no test reaches. And `status_line` lives in `state.rs` so that a test
/// asserting on a string is asserting on what the window paints — a `format!`
/// migrating into `panel.rs` leaves every such test green over a window they
/// have stopped describing, which is silent by construction.
/// Crate names that mean "an engine, a window, or a GPU".
///
/// **The list is of things that must never be reachable from a result**, not of
/// things that are large or unwelcome. Each entry either owns a frame loop,
/// opens a window, or talks to a GPU — and every one of them brings scheduling,
/// hashing or iteration order that §13.4 cannot survive.
///
/// `bevy` is the entry this list was written for. Its ECS runs systems in
/// parallel by default and its query iteration order is **not** insertion order:
/// measured on this workspace, spawning 64 entities, despawning every third and
/// spawning replacements gives back `[49, 1, 2, 62, 4, 5, 61, ...]`, because
/// entity slots are reused. That is exactly the churn a viewer rebuilding its
/// scene on a seed change produces. Anything result-affecting that read that
/// order would diverge across platforms with nothing failing.
const ENGINE_CRATES: &[&str] = &[
    "bevy", "eframe", "egui", "epaint", "wgpu", "winit", "fyrox", "three-d",
];

/// No crate but the viewer may depend on an engine, a window or a GPU.
///
/// **The physics, the chemistry and the generation are ours, and this is the
/// line that says so in a form the build can check.** The engine displays and
/// takes input; it computes nothing. Stated as a rule about *manifests* because
/// that is the only place the boundary can be crossed without writing a line of
/// code — one dependency entry is all it takes for `bevy`'s scheduler to become
/// reachable from `borbax-molecule`.
///
/// Matched as a prefix on the dependency name so `bevy_egui`, `bevy_render`,
/// `egui_kittest` and `wgpu-hal` are all caught by their base entry, and as a
/// whole token at the start so `bevyish` is not a false positive.
///
/// **Both halves are here.** The negative one refuses an engine anywhere but
/// the viewer; the positive one refuses a tree with *no* engine at all, so the
/// guard cannot pass by the dependency having been deleted. An earlier version
/// of this paragraph said the positive half was still to come — it landed with
/// the dependency, and this doc was not updated alongside it.
fn check_the_engine_stays_in_the_viewer(
    root: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    const VIEWER: &str = "crates/borbax-ui";
    let viewer_manifest = root.join(VIEWER).join("Cargo.toml");
    let mut manifests_examined = 0_usize;

    for scan_root in DATA_FREE_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            failures.push(format!(
                "§13.4: scan root {scan_root:?} does not exist, so no manifest under it \
                 was checked for an engine dependency"
            ));
            continue;
        }
        for manifest in walk(&dir)?
            .into_iter()
            .filter(|p| p.file_name().is_some_and(|n| n == "Cargo.toml"))
        {
            // Counted where a manifest is *read*, not where a directory is
            // entered — the Task 8 `pairs` defect, which counted the outer loop
            // and reported a bar four times weaker than its message claimed.
            manifests_examined += 1;
            if manifest == viewer_manifest {
                continue;
            }
            let src = std::fs::read_to_string(&manifest).map_err(|e| e.to_string())?;
            for dep in manifest_dependency_names(&code_only(&src)) {
                let Some(engine) = ENGINE_CRATES.iter().find(|e| {
                    dep == **e
                        || dep.starts_with(&format!("{e}_"))
                        || dep.starts_with(&format!("{e}-"))
                }) else {
                    continue;
                };
                failures.push(format!(
                    "§13.4: {} depends on `{dep}` ({engine}). The engine displays and \
                     takes input; it computes nothing. Bevy's ECS schedules systems in \
                     parallel and its query iteration order is not insertion order, so \
                     a result path that reaches it diverges across platforms with every \
                     test green. If this crate genuinely needs to draw, it is the \
                     viewer's job to call it — never the reverse",
                    manifest.strip_prefix(root).unwrap_or(&manifest).display()
                ));
            }
        }
    }

    // **The positive half, and it is not symmetry for its own sake.** Every
    // check above is a *negative*: it fires when a crate that should not name an
    // engine does. Negatives alone cannot tell "the engine is confined to the
    // viewer" from "there is no engine anywhere" — delete Bevy from the viewer's
    // manifest and every assertion above passes, loudly reporting success over a
    // program that no longer draws. This is `check_blocklist_present`'s "Loud,
    // not `Ok(())`" lesson, and the degeneration `check_wall_clock_has_one_home`
    // is documented as unable to see.
    let viewer_src = std::fs::read_to_string(&viewer_manifest).map_err(|e| e.to_string())?;
    let viewer_deps = manifest_dependency_names(&code_only(&viewer_src));
    if !viewer_deps.iter().any(|d| {
        ENGINE_CRATES
            .iter()
            .any(|e| d == *e || d.starts_with(&format!("{e}_")))
    }) {
        failures.push(format!(
            "§13.4: {} names no engine, so this check passed by there being nothing to \
             confine rather than by the confinement holding. If the viewer has stopped \
             drawing, say so here; if it has moved, move this check with it",
            viewer_manifest
                .strip_prefix(root)
                .unwrap_or(&viewer_manifest)
                .display()
        ));
    }

    // **A scan that read no manifests is a failure, not a pass.** The counter
    // counts manifests, which is what the message says; the workspace has one
    // per member plus the roots, so anything below 5 means the layout moved and
    // this guard checked nothing.
    if manifests_examined < 5 {
        failures.push(format!(
            "§13.4: the engine-containment scan examined only {manifests_examined} \
             manifest(s), which is fewer than this workspace has members. The layout \
             has moved and the check is blind"
        ));
    }
    Ok(())
}

/// Every dependency name declared by a manifest, including renamed entries.
///
/// A rename hides the crate behind a key of the author's choosing —
/// `renderer = { package = "bevy" }` declares `bevy` under a name no scan of
/// the keys would see — so the `package = "..."` value is read in preference to
/// the key whenever one is present.
fn manifest_dependency_names(manifest: &str) -> Vec<String> {
    /// What the current `[header]` puts us inside.
    enum Section {
        /// `[dependencies]` and friends: every `key = ..` line is a crate.
        Table,
        /// `[dependencies.foo]`: the header names one crate, and the lines
        /// inside are *its fields* — only `package` renames it.
        Entry(usize),
        Other,
    }

    let mut out: Vec<String> = Vec::new();
    let mut section = Section::Other;

    for line in manifest.lines() {
        let trimmed = line.trim();

        if let Some(header) = trimmed.strip_prefix('[').and_then(|h| h.strip_suffix(']')) {
            // `dev-dependencies` and `build-dependencies` count too: an engine
            // in a chemistry crate's test binaries is an engine a result could
            // be minted next to.
            section = match header.trim().rsplit_once('.') {
                // `[target.'cfg(..)'.dependencies]` -> a table.
                Some((_, last)) if last.ends_with("dependencies") => Section::Table,
                // `[dependencies.foo]`, `[target...dev-dependencies.foo]`.
                Some((prefix, name)) if prefix.ends_with("dependencies") => {
                    out.push(name.trim_matches('"').to_owned());
                    Section::Entry(out.len() - 1)
                }
                _ if header.trim().ends_with("dependencies") => Section::Table,
                _ => Section::Other,
            };
            continue;
        }

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let Some((key, value)) = trimmed.split_once('=') else {
            continue;
        };
        let key = key.trim();

        match section {
            Section::Table => {
                // A `package = "..."` inline renames the crate, so the key is
                // not the crate. Prefer the rename.
                out.push(inline_package(value).unwrap_or_else(|| key.trim_matches('"').to_owned()));
            }
            Section::Entry(at) => {
                if key == "package"
                    && let Some(real) = quoted(value)
                    && let Some(slot) = out.get_mut(at)
                {
                    *slot = real;
                }
            }
            Section::Other => {}
        }
    }
    out
}

/// The `package = "x"` inside an inline table such as
/// `renderer = { package = "bevy", version = "0.19" }`.
fn inline_package(value: &str) -> Option<String> {
    let at = value.find("package")?;
    quoted(value.get(at..)?)
}

/// The first double-quoted run in `s`.
fn quoted(s: &str) -> Option<String> {
    let open = s.find('"')?;
    let rest = s.get(open + 1..)?;
    let close = rest.find('"')?;
    Some(rest.get(..close)?.to_owned())
}

fn check_the_viewer_stays_a_leaf(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    const VIEWER: &str = "crates/borbax-ui";

    let viewer_src = root.join(VIEWER).join("src");
    if !viewer_src.exists() {
        // Not a failure: the viewer is one crate among several and a tree
        // without it is a legitimate state (it did not exist before Step 1).
        // Unlike a *scan root*, its absence does not mean another check was
        // blinded.
        return Ok(());
    }

    // 1. Nothing depends on the viewer.
    for scan_root in DATA_FREE_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            // Loud, for the reason its two siblings above are loud — and this
            // one is the second-order lesson: those two were *changed to be*
            // loud in this same commit, and then this function was split out
            // twenty lines below them carrying the silent version. Fixing a
            // pattern in the places you can see is not fixing the pattern.
            failures.push(format!(
                "§13.1: scan root {scan_root:?} does not exist, so no manifest under it \
                 was checked for a dependency on `borbax-ui`"
            ));
            continue;
        }
        for manifest in walk(&dir)?
            .into_iter()
            .filter(|p| p.file_name().is_some_and(|n| n == "Cargo.toml"))
            .filter(|p| !p.starts_with(root.join(VIEWER)))
        {
            let src = std::fs::read_to_string(&manifest).map_err(|e| e.to_string())?;
            if code_only(&src).contains("borbax-ui") {
                failures.push(format!(
                    "§13.1: {} depends on `borbax-ui`. Every relaxation that crate \
                     documents for itself — f32, HashMap, a wall-clock read — is \
                     licensed by nothing depending on it",
                    manifest.strip_prefix(root).unwrap_or(&manifest).display()
                ));
            }
        }
    }

    Ok(())
}

/// A whole file with comments, string literals and character literals removed.
///
/// **Stronger than [`code_only`], and the seam check needs the difference.**
/// `code_only` drops whole-line `//` comments and nothing else, so a type named
/// inside a string literal reads as code. That is not hypothetical: `lib.rs`
/// carries an `#![expect]` attribute whose `reason` string names `eframe` while
/// explaining why a lint fires, which is documentation and is
/// exactly as much prose as the `//!` block above it. Scanning it as code failed
/// the seam on correct code — and a guard that fires on correct code gets
/// deleted, which is the structural pressure that makes this worth getting
/// right rather than exempting.
///
/// Naming a type inside a string literal is definitionally not using it: there
/// is no import, no call and no path resolution through a `&str`. Whereas
/// `format!` survives, because the macro name sits *outside* the literal it
/// builds.
///
/// Returns `None` when the file ends mid-string or mid-comment, so a file that
/// cannot be lexed is reported rather than silently treated as empty.
fn code_without_prose(src: &str) -> Option<String> {
    let mut lex = LexState::default();
    let stripped = src
        .lines()
        .map(|line| strip_comments_and_literals(line, &mut lex))
        .collect::<Vec<_>>()
        .join("\n");
    lex.is_clean().then_some(stripped)
}

/// The one file under `crates/borbax-ui/src` that may open a window.
const VIEWER_SHELL_FILE: &str = "main.rs";

/// The files that may name the engine.
///
/// **`main.rs` alone was the rule and it was the wrong rule**, for a reason
/// that cost an empty window. A binary crate root cannot be imported, so wiring
/// kept in `main.rs` is wiring no test can reach — and on 2026-08-03 the app
/// shipped without the camera `bevy_egui` renders through, with 509 tests green
/// and nothing on screen. Moving the construction into `app.rs` is what lets
/// `tests/app.rs` assert about the *same* `App` the binary runs.
///
/// So the licence is two files, and `main.rs` is not one of them — **measured,
/// not assumed**. It held the licence until a review removed it and found the
/// whole gate still green: the file is one statement and names no engine, and
/// the count twelve lines below already removed it before comparing. A licence
/// nobody uses is permission for the next person to put a decision back into
/// the one file no test can import, which is the defect this branch exists
/// downstream of.
/// Everything else in the crate still names no engine at all.
///
/// **`scene.rs` was added deliberately at Step 3, which is the only way a file
/// gets here.** The default tier is the strict one, so the new file failed three
/// checks the moment it existed and had to be argued onto this list rather than
/// arriving on it. What earns it the licence: the 3D scene is entities, meshes,
/// cameras and input — all of it engine vocabulary with no honest spelling that
/// avoids the engine. What is *not* on this list, and stays off it deliberately,
/// is the arithmetic behind it: the molecule is built in `molecule.rs` and the
/// camera's angles live in `orbit.rs`, both in the strictest tier and both
/// testable with no app at all.
const VIEWER_ENGINE_FILES: &[&str] = &["app.rs", "scene.rs"];

/// The files that draw, and may therefore name `egui`.
///
/// **Membership is a licence, not a description.** Being on this list is what
/// grants a file `egui`; it costs the file `format!` in exchange, because a
/// string built where it is painted is a string no test can assert on without
/// describing a window it has stopped describing. Adding a name here is the
/// deliberate act; forgetting to is caught by [`viewer_banned_imports`]'s
/// default rather than by nobody.
const VIEWER_DRAWING_FILES: &[&str] = &["panel.rs"];

/// The files that may name `wgpu`, relative to `crates/borbax-ui/src`.
///
/// **Empty is the correct state, not an unfinished one.** No file names `wgpu`
/// today, so the licence list is empty and every file is refused. Step 3's
/// offscreen renderer adds one entry — `scene/draw.rs` — in the commit that
/// creates the file, which is what makes granting the licence a visible act
/// rather than a side effect of writing an import.
///
/// It governs [`check_the_viewer_seam_holds`]'s `wgpu` count as well as its ban,
/// so it binds [`VIEWER_SHELL_FILE`] too even though the shell's ban list is
/// empty — exactly as the engine is pinned to the files whose own ban list is
/// empty.
const VIEWER_GPU_FILES: &[&str] = &[];

/// What a file under `crates/borbax-ui/src` may not name as a whole identifier.
///
/// **Separate from [`viewer_banned_imports`] because the match rule differs, and
/// the difference is load-bearing in both directions.**
///
/// `egui` is matched as a *substring* on purpose: that is what catches
/// `egui_wgpu` and `egui_kittest` in a file with no `egui` licence, and
/// narrowing it to an identifier would open both.
///
/// These two cannot take that rule. A substring `wgpu` counts `egui_wgpu`,
/// which makes the `egui` rule and the `wgpu` rule indistinguishable and reports
/// a number one higher than the thing it names for every drawing file. And a
/// substring `epaint` fires on **`request_repaint`** — measured: `panel.rs`
/// already contains `repaints` in prose today, and `ctx.request_repaint()` is
/// ordinary correct code a continuously-animating scene will write. A guard that
/// fires on correct code gets deleted, which is the structural pressure that
/// makes the boundary worth getting right rather than exempting.
///
/// **Why `epaint` at all**, when nothing names it: it is `egui`'s painting
/// layer, re-exported through `egui` itself, so `use epaint::PaintCallbackInfo;`
/// satisfies the strictest tier's letter — it names no `egui` — while defeating
/// its purpose. It carries no licence list because no file has a reason to want
/// it; the day one does, it needs one.
///
/// **These bind [`VIEWER_SHELL_FILE`] as well**, which is the one place this
/// function deliberately disagrees with [`viewer_banned_imports`]. The shell's
/// empty ban list is justified by "no headless test reaches it, so there is
/// nothing to protect" — true of the engine, and the licence for `wgpu` is a
/// different question: it is which files were *chosen* to hold GPU code. The
/// shell may well want `wgpu` for the render state; when it does it joins
/// [`VIEWER_GPU_FILES`] like any other file.
fn viewer_banned_idents(file: &str) -> Vec<&'static str> {
    banned_idents_given(file, VIEWER_GPU_FILES)
}

/// [`viewer_banned_idents`] with the licence list passed in.
///
/// **Split out so the licensed arm is reachable from a test.**
/// [`VIEWER_GPU_FILES`] is empty today, which makes `contains` always false and
/// the exemption branch dead — so the obvious inversion of the condition
/// (`if contains` rather than `if !contains`, banning `wgpu` in exactly the
/// files allowed to have it and nowhere else) is invisible to any test written
/// against the real constant. Handing the list in is what lets both arms be
/// exercised before a file needs the licence.
fn banned_idents_given(file: &str, gpu_files: &[&str]) -> Vec<&'static str> {
    let mut banned = vec!["epaint"];
    if !gpu_files.contains(&file) {
        banned.push("wgpu");
    }
    // **`bevy` must be an identifier needle and cannot be a substring one**,
    // which is why it lives here rather than in `viewer_banned_imports`. The
    // drawing tier is *allowed* `bevy_egui` — that is where `egui` now comes
    // from — and `"bevy_egui".contains("bevy")` is true, so a substring ban
    // would fire on the one import the tier exists to permit.
    if !VIEWER_ENGINE_FILES.contains(&file) {
        banned.push("bevy");
    }
    banned
}

/// What a file under `crates/borbax-ui/src` may not name in code.
///
/// **The default is the strictest tier, and that inversion is the whole
/// design.** The shipped version of this rule was a two-entry table, so a new
/// `.rs` file — which is exactly what Step 2 wanted to add — was checked for
/// *nothing*: it could name the engine, or build strings while calling itself a
/// drawing file. Enumerating the files instead moves the hand-kept-list failure
/// mode down into the entry, where `("periodic.rs", &[])` is one word wide and
/// builds green. Defaulting to strict means the way to get a licence is to ask
/// for one.
///
/// Three tiers:
///
/// - unnamed (**default**) — names no UI type at all. This is `state.rs`'s rule
///   and it is an *identity*, not a consequence of what the file happens to
///   import, so it cannot be derived from the file's own contents.
/// - [`VIEWER_DRAWING_FILES`] — may name `egui` (which lays out on the CPU and
///   cannot open a window, so `egui_kittest` drives it headlessly), may not
///   build strings.
/// - [`VIEWER_SHELL_FILE`] — the shell that can open a window, so no test
///   reaches it and there is nothing to protect.
fn viewer_banned_imports(file: &str) -> &'static [&'static str] {
    if VIEWER_ENGINE_FILES.contains(&file) {
        // **`format!` is still forbidden, and that is not symmetry.** The
        // engine files may name `bevy` and `egui` because wiring the bridge
        // needs both; they may not build strings, for the same reason
        // `panel.rs` may not. A string built where it is painted is a string no
        // test can assert on without describing a window it has stopped
        // describing, and that rule is about *where the string is made*, not
        // about which crate the file is allowed to import.
        &["format!"]
    } else if VIEWER_DRAWING_FILES.contains(&file) {
        &["format!"]
    } else {
        &["egui"]
    }
}

/// The viewer's per-file import seam, and its single `Universe::generate` site.
///
/// Split out of [`check_the_viewer_stays_a_leaf`] because the combined function
/// crossed `clippy::too_many_lines`, which is the right instinct here — these
/// are three unrelated invariants that happen to share a subject.
fn check_the_viewer_seam_holds(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    const VIEWER: &str = "crates/borbax-ui";
    let viewer_src = root.join(VIEWER).join("src");
    // **A missing directory is a failure, not a skip**, and the first version of
    // this function returned `Ok(())` here. That is the same fail-open shape as
    // the `crates/` rename that blinded `check_every_member_inherits_the_lints`
    // while eight other guards fired: rename or move `crates/borbax-ui` and the
    // *entire* seam check — every tier, the engine count and the
    // `Universe::generate` count — passes silently, with nothing in the output
    // saying it checked nothing.
    if !viewer_src.exists() {
        failures.push(format!(
            "viewer seam: {VIEWER}/src does not exist, so the seam was not checked. \
             If the viewer has moved, move this check with it; a silent pass here \
             disables the per-file seam, the engine count and the \
             `Universe::generate` count at once"
        ));
        return Ok(());
    }

    // The per-file import seam, checked over code rather than prose, over
    // **every** `.rs` file rather than a list of two.
    let mut bevy_homes = Vec::new();
    let mut wgpu_homes = Vec::new();
    for path in walk(&viewer_src)?
        .into_iter()
        .filter(|p| p.extension().is_some_and(|e| e == "rs"))
    {
        // **Relative to `src`, not the bare file name**, and the difference is a
        // licence. `walk` recurses, so a nested `src/widgets/main.rs` has
        // `file_name() == "main.rs"` and would be handed the *unrestricted*
        // tier — the laxest rule applied to a file nobody named, which inverts
        // the unknown-is-strict rule this whole check is built on. A nested
        // `panel.rs` claims the drawing licence the same way, and the engine
        // count is fooled twice over: two files called `main.rs` count as
        // `["main.rs", "main.rs"]`, and one hidden in a subdirectory passes as a
        // legitimate `["main.rs"]`.
        let name = path
            .strip_prefix(&viewer_src)
            .unwrap_or(&path)
            .to_string_lossy()
            .replace('\\', "/");
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        let Some(code) = code_without_prose(&src) else {
            failures.push(format!(
                "viewer seam: {VIEWER}/src/{name} ended inside a string or comment, so \
                 it could not be scanned. Reported rather than skipped: a file the \
                 lexer gives up on is the one file most likely to be hiding something"
            ));
            continue;
        };

        for needle in viewer_banned_imports(&name) {
            if code.contains(needle) {
                failures.push(format!(
                    "viewer seam: {VIEWER}/src/{name} contains `{needle}` in code. \
                     Files name no UI type by default, named drawing files may name \
                     `egui` but build no strings, and the engine files may name \
                     both and build none — that split is what lets `egui_kittest` drive the \
                     real drawing code with no window, and what keeps every string \
                     assertion pointed at what the window actually paints. If this is \
                     a new drawing file, add it to `VIEWER_DRAWING_FILES` and say why \
                     it draws"
                ));
            }
        }

        for needle in viewer_banned_idents(&name) {
            if names_type(&code, needle) {
                failures.push(format!(
                    "viewer seam: {VIEWER}/src/{name} names `{needle}` in code. `wgpu` \
                     belongs only to the files listed in `VIEWER_GPU_FILES`, and \
                     `epaint` — `egui`'s painting layer, which the strictest tier's \
                     `egui` needle does not see — belongs to none. Matched at \
                     identifier boundaries, so `egui_wgpu` and `request_repaint` are \
                     not this. If this file is meant to hold GPU code, add it to \
                     `VIEWER_GPU_FILES` and say what it renders"
                ));
            }
        }

        if names_type(&code, "bevy") {
            bevy_homes.push(name.clone());
        }
        if names_type(&code, "wgpu") {
            wgpu_homes.push(name);
        }
    }

    // **The engine is named in code by exactly the licensed files.**
    //
    // This is the half a per-file ban list cannot enforce, and before this it
    // was enforced by *nobody*: "only one file names the engine" held only
    // because there happened to be three files and two of them banned the
    // string. A fourth file added with a lax entry — `("periodic.rs", &[])` —
    // reopens it with a green build, which is the hand-kept-list failure mode
    // moved from the file level down into the entry, where it is *less* visible.
    // Counting the homes is what makes the invariant independent of the list.
    bevy_homes.sort();
    let mut engine_files: Vec<&str> = VIEWER_ENGINE_FILES.to_vec();
    engine_files.sort_unstable();
    // `main.rs` is one statement long and names the engine only through
    // `borbax_ui::app`, so it does not appear here — the count is of files that
    // name `bevy` *in code*, which is `app.rs` alone.
    engine_files.retain(|f| *f != VIEWER_SHELL_FILE);
    if bevy_homes != engine_files {
        failures.push(format!(
            "viewer seam: `bevy` is named in code by {bevy_homes:?} under \
             {VIEWER}/src, expected exactly {engine_files:?}. The engine opens \
             a window and owns the frame loop, so every file that names it is a file \
             no headless test can reach; confining it to one is what keeps the drawing \
             code testable and the decisions in `state.rs` reachable with no app at \
             all. Matched as a whole identifier, so `bevy_egui` — which the drawing \
             tier is deliberately allowed — is not this"
        ));
    }

    // **`wgpu` is named by exactly the files licensed to name it**, and the
    // count is the half the per-file ban list cannot do.
    //
    // The ban list is keyed on a tier, and being in a tier is a licence for
    // *that tier's* needles: a file added to `VIEWER_DRAWING_FILES` gets
    // `["egui", "format!"]` and — before this — nothing about `wgpu` at all.
    // So the ban alone lets a lax entry reopen the rule with a green build,
    // which is precisely why `wgpu_homes` exists two lines up. Counting makes
    // the invariant independent of the tier table.
    //
    // **The equality is against the licence list, so it fails in both
    // directions.** An unlicensed file naming `wgpu` fails; a licensed file that
    // has stopped naming it also fails, because a stale entry is a licence
    // nobody is using and the next file added under that name inherits it.
    wgpu_homes.sort();
    let mut licensed: Vec<&str> = VIEWER_GPU_FILES.to_vec();
    licensed.sort_unstable();
    if wgpu_homes != licensed {
        failures.push(format!(
            "viewer seam: `wgpu` is named in code by {wgpu_homes:?} under {VIEWER}/src, \
             expected exactly {licensed:?}. Rendering code is confined to named files so \
             the rest of the crate stays drivable with no GPU — `state.rs` in particular \
             is the file whose entire identity is that a test can reach it without a \
             window or an adapter. Note this counts whole identifiers: `egui_wgpu` is \
             caught by the `egui` needle instead, in every tier that bans it"
        ));
    }

    check_the_viewer_calls_the_chemistry_once(root, &viewer_src, failures)?;
    check_the_chemistry_has_one_door(root, &viewer_src, failures)?;
    check_the_viewer_narrows_in_one_place(&viewer_src, failures)?;
    Ok(())
}

/// The viewer's call-site counts: one `Universe::generate`, no `.pattern()`.
///
/// Split from [`check_the_viewer_seam_holds`] when the per-file tiers pushed the
/// combined function past `clippy::too_many_lines` — the same split, for the
/// same reason, that separated the seam from
/// [`check_the_viewer_stays_a_leaf`].
fn check_the_viewer_calls_the_chemistry_once(
    root: &Path,
    viewer_src: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    const VIEWER: &str = "crates/borbax-ui";

    // **`PeriodicTable::pattern()` has no consumer in the viewer at Step 2.**
    //
    // The attractive wrong implementation is building the grid's rows from
    // `pattern().closures`: the closure list gives exactly the row boundaries,
    // it *works*, and it makes the viewer derive its layout from the shell law
    // instead of reading `period`/`group` off each element. That is a second
    // encoding of the same fact in the crate least entitled to hold one, and it
    // opens the door to `peak` and `k` reaching the screen as undisclosed
    // physics.
    //
    // **This check's only correct response is deletion**, and it says so rather
    // than being worked around: `peak` and `closures` are legitimate things a
    // later step will want to show.
    // **`(` without the closing paren.** Requiring `)` matched only an *empty*
    // argument list, so `x.pattern(&table)` scored nothing and the guard would
    // go silent the day `pattern` took an argument — the fail-open direction
    // this check exists to close. The `(` stays, because `.pattern` as a *field*
    // is not a call.
    let (pattern_total, pattern_sites) =
        call_sites_under(root, viewer_src, PATTERN_CALL, ".pattern()", failures)?;
    if pattern_total > 0 {
        failures.push(format!(
            "viewer seam: `.pattern()` is called in {pattern_sites:?}. Step 2 has no \
             consumer for `PeriodicTable::pattern()` — the grid reads `period` and \
             `group` off each element, and deriving the layout from the shell law \
             instead would be the viewer holding a second encoding of the physics. If \
             you are adding a legitimate consumer, delete this check and say what it \
             shows"
        ));
    }

    // Exactly one `Universe::generate` call site.
    //
    // **This is the half the `regenerations` counter cannot cover, and the
    // counter's own doc says so.** A `Universe::generate` written directly into
    // `panel::draw` never touches `reload`, so it never increments the counter:
    // `an_unchanged_seed_box_does_not_rebuild_the_universe` stays green over a
    // 141 us call running every frame, which is 0.85% of a frame budget — too
    // small to see and invisible to every test. The plan deferred this grep to
    // Step 7; the surface it guards exists now.
    let (total, where_) = call_sites_under(
        root,
        viewer_src,
        GENERATE_CALL,
        "Universe::generate",
        failures,
    )?;
    if total != 1 {
        failures.push(format!(
            "viewer seam: `Universe::generate` has {total} call site(s) in {VIEWER}/src \
             ({where_}), expected exactly 1 (in state.rs, inside `reload`). A second one \
             in the paint body would run ~141 us every frame and increment no counter, \
             so no test would see it. Note this matches the literal spelling: a \
             `use borbax_universe::Universe as U; U::generate(..)` alias passes, which is \
             a hole a symbol-aware check would close"
        ));
    }

    // Exactly one `embed` and one `canonicalise` call site.
    //
    // **The same argument as `Universe::generate`, at the price that made §8.6
    // a rule.** Both are per-species work: `canonicalise` searches for a
    // canonical labelling and `embed` runs 240 fixed stress-majorization
    // iterations. Written into a paint body they would run every frame, and
    // `re_embeds` would never see it — a counter on a state method is untouched
    // by a call that does not go through the state method, which is exactly what
    // the `regenerations` doc says about its own blind spot one function above.
    //
    // **Two counters, never one over the pair**, because they are not the same
    // quantity and will stop being called together: `canonicalise` returns a
    // `Result` and can fail, so a combined count would report "embeddings" while
    // measuring "attempts" the first time one of them is retried.
    for (needle, what, site) in [
        (EMBED_CALL, "embed", "molecule.rs, inside `build`"),
        (
            CANONICALISE_CALL,
            "canonicalise",
            "molecule.rs, inside `build`",
        ),
    ] {
        let (total, where_) = call_sites_under(root, viewer_src, needle, what, failures)?;
        if total != 1 {
            failures.push(format!(
                "viewer seam: `{what}` has {total} call site(s) in {VIEWER}/src ({where_}), \
                 expected exactly 1 (in {site}). This counts **call sites in shipped \
                 code**, which is not the same as work done — `if x {{ {what}(a) }} else \
                 {{ {what}(b) }}` is two sites and one call — so it bounds where the \
                 chemistry is reached from and says nothing about how often. That bound \
                 is `ViewerState::re_embeds`'s job, and it cannot see a call written \
                 straight into a paint body or a frame system, which is why this exists. \
                 The hole, stated so the right probe gets run: an alias hiding the *only* \
                 site gives 0 and fails, so it is not that — it is an alias hiding a \
                 *second* site (`use borbax_molecule::{what} as lay_out;`) while a \
                 literal one keeps the count at 1"
            ));
        }
    }
    Ok(())
}

/// Count `needle` across every shipped `.rs` file under `viewer_src`.
///
/// Returns the total and a human-readable site list for the failure message.
///
/// **A file that cannot be lexed is a failure, never a zero.** The `.pattern()`
/// count used to `continue` silently where the `Universe::generate` half failed
/// loudly; unifying them makes all three loud.
///
/// **That is defence in depth, not a closed live hole**, and the difference was
/// measured rather than reasoned: an unlexable file under `src` aborts the whole
/// `xtask` run at an earlier `syn::parse_file`, so the seam counting is never
/// reached and the silent `continue` was unreachable. It is worth keeping
/// because that ordering is not a property anyone declared — it changes the day
/// a check moves or narrows.
///
/// **Shipped code only.** A `#[cfg(test)]` module that legitimately calls one of
/// these would otherwise fire the guard on correct code, and a guard that fires
/// on correct code gets deleted — which is the structural pressure that made the
/// last one worth evading.
fn call_sites_under(
    root: &Path,
    viewer_src: &Path,
    needle: &[&str],
    what: &str,
    failures: &mut Vec<String>,
) -> Result<(usize, String), String> {
    let mut sites = Vec::new();
    for path in walk(viewer_src)?
        .into_iter()
        .filter(|p| p.extension().is_some_and(|e| e == "rs"))
    {
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        let Some(hits) = count_outside_tests(&src, needle) else {
            failures.push(format!(
                "viewer seam: crates/borbax-ui/src/{} could not be lexed, so its \
                 `{what}` calls were not counted",
                path.strip_prefix(viewer_src).unwrap_or(&path).display()
            ));
            continue;
        };
        if hits > 0 {
            sites.push((path, hits));
        }
    }
    let total = sites.iter().map(|(_, n)| n).sum();
    let where_ = sites
        .iter()
        .map(|(p, n)| format!("{}x{n}", p.strip_prefix(root).unwrap_or(p).display()))
        .collect::<Vec<_>>()
        .join(", ");
    Ok((total, where_))
}

/// The chemistry is reached from exactly one file in the viewer.
///
/// **This is the §8.6 guard, and it replaces an enumeration that covered two of
/// six operations.** The previous form counted call sites of `embed` and
/// `canonicalise` by name. §8.6 names six per-species operations —
/// canonicalisation, embedding, signature construction, folding, cavity
/// extraction and per-bond decay — and two of the four uncounted ones,
/// `signature` and `affinity`, are reachable *today*, because `borbax-molecule`
/// is a dependency of the viewer. Folding and cavity extraction arrive later
/// into the same crate and would have arrived unguarded.
///
/// Worse, the enumeration had a hole its own message disclosed and a review
/// then demonstrated: an aliased `embed`, called from a system that runs **every
/// frame**, passed the whole gate with `re_embeds` flat. All three §8.6 guards
/// were blind at once.
///
/// **So the default is inverted, the way the seam tiers already do it**: instead
/// of naming the functions that are forbidden, name the files that are allowed
/// to reach the chemistry crate at all. An alias still has to name the crate to
/// import from it, so the hole closes; operations that do not exist yet are
/// covered the day they are written; and the rule is one line rather than a list
/// that has to be kept in step with a spec section.
///
/// Measured when written: `borbax_molecule` appears in exactly one file under
/// `crates/borbax-ui/src`, so this needs no exemption list.
fn check_the_chemistry_has_one_door(
    root: &Path,
    viewer_src: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    const CHEMISTRY_DOOR: &str = "molecule.rs";
    /// The crate every per-species operation lives in, as one token.
    const CHEMISTRY_CRATE: &[&str] = &["borbax_molecule"];
    let mut doors = Vec::new();
    for path in walk(viewer_src)?
        .into_iter()
        .filter(|p| p.extension().is_some_and(|e| e == "rs"))
    {
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        // **Counted as a token, not as text.** `identifier_segments` splits on
        // underscores, so it never sees the crate's name whole — the first
        // version of this check found zero doors in a tree that has one.
        let Some(hits) = count_outside_tests(&src, CHEMISTRY_CRATE) else {
            failures.push(format!(
                "viewer seam: crates/borbax-ui/src/{} could not be lexed, so its \
                 chemistry imports were not counted",
                path.strip_prefix(viewer_src).unwrap_or(&path).display()
            ));
            continue;
        };
        if hits > 0 {
            doors.push(
                path.strip_prefix(viewer_src)
                    .unwrap_or(&path)
                    .display()
                    .to_string(),
            );
        }
    }
    doors.sort();
    if doors != vec![CHEMISTRY_DOOR.to_owned()] {
        failures.push(format!(
            "viewer seam: `borbax_molecule` is named in code by {doors:?} under \
             crates/borbax-ui/src, expected exactly [\"{CHEMISTRY_DOOR}\"]. Every \
             per-species operation the chemistry owns — canonicalisation, \
             embedding, signatures, folding, cavities, decay rates — is work that \
             must happen once per species and never inside a frame, and confining \
             the *import* is what covers the ones that do not exist yet. The \
             counter on `ViewerState` cannot help: a call written straight into a \
             frame system never goes through it. If a second file genuinely needs \
             the chemistry, add it here and say what it computes and when"
        ));
    }
    let _ = root;
    Ok(())
}

/// The viewer narrows `f64` to `f32` in exactly one place.
///
/// **`scene.rs`'s own doc claimed this check existed before it did.** A second
/// `const fn shrink2` beside the first passed the entire gate; what was actually
/// enforcing anything was `clippy::as_conversions`, which demands an `#[expect]`
/// per site and says nothing about how many sites there are.
///
/// The invariant is worth holding: the narrowing is a display requirement, and a
/// display requirement spreading through code that handles chemistry values is
/// how a rounding choice reaches a number that came out of the solver.
///
/// The whole-pixel narrowing for the viewport is a *different* cast (`f32` to
/// `u32`, forced by the engine's viewport API) and is deliberately not counted
/// here.
fn check_the_viewer_narrows_in_one_place(
    viewer_src: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    const NARROW_CAST: &[&str] = &["x", "as", "f32"];
    let mut sites = Vec::new();
    for path in walk(viewer_src)?
        .into_iter()
        .filter(|p| p.extension().is_some_and(|e| e == "rs"))
    {
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        let Some(hits) = count_outside_tests(&src, NARROW_CAST) else {
            failures.push(format!(
                "viewer seam: crates/borbax-ui/src/{} could not be lexed, so its \
                 `as f32` casts were not counted",
                path.strip_prefix(viewer_src).unwrap_or(&path).display()
            ));
            continue;
        };
        if hits > 0 {
            sites.push((
                path.strip_prefix(viewer_src)
                    .unwrap_or(&path)
                    .display()
                    .to_string(),
                hits,
            ));
        }
    }
    let total: usize = sites.iter().map(|(_, n)| n).sum();
    if total != 1 {
        failures.push(format!(
            "viewer seam: the viewer narrows `f64` to `f32` at {total} site(s) \
             ({sites:?}), expected exactly 1 (in scene.rs, `narrow`). Every \
             coordinate, radius and camera position crosses that one function, \
             which is what keeps a display-driven rounding from spreading into \
             code that handles a value the solver produced. Note the needle \
             matches the binding `x as f32` specifically, so it is the declared \
             boundary being counted and not every cast in the crate"
        ));
    }
    Ok(())
}

/// Wall-clock is readable in exactly one file under [`DATA_FREE_ROOTS`].
///
/// **Two limits, stated because the first version of this line said "the whole
/// workspace" and meant neither.** `xtask` is not scanned — it is a build tool
/// whose output cannot reach a simulation result, the same exclusion
/// [`TRANSCENDENTAL_SCAN_ROOTS`] makes and for the same reason. And the match is
/// textual over two spellings, so `use std::time::Instant as Clock; Clock::now()`,
/// `SystemTime::UNIX_EPOCH.elapsed()` and a bare `.elapsed()` on a held `Instant`
/// all pass. That is a real hole and the guard is still worth having: it catches
/// the spelling anyone actually writes, and it converts a second wall-clock read
/// from an invisible act into one that has to be disguised.
///
/// **§13.1 bans wall-clock and, until this, nothing enforced it.** The ban
/// sits in CLAUDE.md's hard invariants beside things that all have guards —
/// `clippy.toml`'s `disallowed-methods` list is 129 lines of `f64::`/`f32::`
/// entries and mentions no time type at all, and `xtask`'s `BANNED_CALLS`
/// mirrors it. So the invariant was held by nobody having written one, and
/// viewer Step 1 is the first read in the workspace's history.
///
/// **Why not add `SystemTime::now` to `clippy.toml` instead**, which is the
/// obvious move: a cross-check test requires everything in `clippy.toml` to
/// appear in `BANNED_CALLS`, and the transcendental scan has **no path
/// exemption by design** — its own comment calls an exemption that currently
/// exempts nothing "a hole with a sign on it". So a `#[expect]` on `moment()`
/// would silence clippy and not `xtask`, and there would be no honest way to
/// let the one sanctioned call through. A single-site rule is the shape that
/// fits: the read is permitted, in one named place, and moving or copying it
/// fails the gate.
///
/// The viewer may read the clock because it produces no results — §13.1
/// protects results and a seed the user is shown and can retype is an input.
/// What must never happen is a *second* read appearing somewhere that does
/// produce results, on the precedent of this one.
fn check_wall_clock_has_one_home(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    const VIEWER: &str = "crates/borbax-ui";
    const CLOCK_READS: &[&str] = &["SystemTime::now", "Instant::now"];
    let clock_home = root.join(VIEWER).join("src").join("state.rs");

    for scan_root in DATA_FREE_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            failures.push(format!(
                "§13.1: scan root {scan_root:?} does not exist, so no source under it was \
                 checked for a wall-clock read"
            ));
            continue;
        }
        for path in walk(&dir)?
            .into_iter()
            .filter(|p| p.extension().is_some_and(|e| e == "rs"))
            .filter(|p| *p != clock_home)
        {
            let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
            let code = code_only(&src);
            for needle in CLOCK_READS {
                if code.contains(needle) {
                    failures.push(format!(
                        "§13.1: `{needle}` in {} — wall-clock is readable in exactly one \
                         file, `{VIEWER}/src/state.rs`, where it seeds a button and \
                         reaches no result. A second read is how it gets into one",
                        path.strip_prefix(root).unwrap_or(&path).display()
                    ));
                }
            }
        }
    }

    Ok(())
}

/// The token run that spells a call into the `libm` crate.
///
/// Tokens rather than the text `"libm::"` so that `#[cfg(test)]` bodies are
/// skipped by the same machinery the viewer's call-site counts already use —
/// see [`check_libm_has_one_home`]'s source half for the probe that made that
/// necessary.
const LIBM_CALL: &[&str] = &["libm", ":", ":"];

/// The one manifest that may declare `libm`.
const LIBM_HOME_MANIFEST: &str = "crates/borbax-units/Cargo.toml";

/// The one source file that may call it.
const LIBM_HOME_SOURCE: &str = "crates/borbax-units/src/det_math.rs";

/// Strip TOML comments so prose naming `libm` is not mistaken for a dependency.
///
/// TOML has no block comments, so this is `#` to end of line — but only when
/// the `#` is outside a string, because a `reason = "..."` may contain one and
/// truncating there would hide whatever followed on that line. Both quote
/// styles and TOML's basic-string escapes are handled; multi-line strings are
/// not, and cannot matter here because a `libm` dependency cannot be written
/// inside one.
fn toml_code_only(src: &str) -> String {
    let mut out = String::with_capacity(src.len());
    for line in src.lines() {
        let mut quote: Option<char> = None;
        let mut escaped = false;
        for c in line.chars() {
            match quote {
                Some(q) => {
                    out.push(c);
                    if escaped {
                        escaped = false;
                    } else if c == '\\' && q == '"' {
                        escaped = true;
                    } else if c == q {
                        quote = None;
                    }
                }
                None if c == '#' => break,
                None => {
                    out.push(c);
                    if c == '"' || c == '\'' {
                        quote = Some(c);
                    }
                }
            }
        }
        out.push('\n');
    }
    out
}

/// §13.1 — `libm` is declared by one crate and called from one file.
///
/// **This closes a hole that passed all six gate legs, measured rather than
/// supposed.** `libm::cos(t)` written in `crates/borbax-ui/src/state.rs`, with
/// one `libm = { workspace = true }` line added to that crate's manifest, gives
/// "all checks passed" from `cargo xtask` and a clean `cargo clippy`.
/// [`BANNED_CALLS`] has no `libm::` needle and `clippy.toml` bans the *inherent
/// method* `f64::sin`, not every function named `sin` — its own header says so,
/// and `leaves_the_libm_free_functions_alone` pins that pass as deliberate.
/// So the gap is a designed property of the transcendental scan, not an
/// oversight in it, and closing it needs a separate guard: a `"libm::"` entry
/// in [`BANNED_CALLS`] would fire on `det_math.rs` itself and force exactly the
/// path exemption that scan exists without.
///
/// **This is a *scope* guard, not a divergence guard, and calling it the latter
/// would be overclaiming.** `libm::sin` written in `borbax-molecule` is
/// bit-identical to `det_math::sin`. Nothing diverges. Three other things break:
///
/// - **The libm-bump review loses its scope.** CLAUDE.md commits to Renovate
///   TIER B delivering libm bumps as `determinism-review` PRs, and that review
///   reads `det_math.rs` and checks each function's arch-dispatch status
///   against the new version. That procedure is correct *only* while
///   `det_math.rs` is the complete list of libm functions this workspace calls.
///   Today it is — by accident, one manifest line from ending, and nothing in
///   the diff that ends it would say the review's scope had just become wrong.
/// - **`libm::acos` is unclamped.** `det_math::acos` clamps, and its doc
///   records that as load-bearing: 685 of 3600 compositions of
///   `rotation_matrices()` give `(tr(R) - 1) / 2 > 1` and would return NaN
///   unclamped. Measured again for the shape a camera actually produces — 500
///   000 near-parallel unit-vector pairs, each independently normalised and
///   perturbed by ~1e-9 — **89 403 (17.88%) have `|dot| > 1`**, worst overshoot
///   4.44e-16. Every one is NaN from `libm::acos` and correct from the wrapper.
///   That is a correctness hazard rather than a determinism one: a NaN camera
///   goes black, deterministically. It is named here because an orbit camera
///   extracting a pitch angle from a dot product is precisely where a direct
///   `libm::acos` gets written.
/// - **The next case is result-affecting.** Folding and cavity extraction will
///   want trig in `borbax-molecule`, which has zero `det_math` call sites
///   today. One manifest line plus `libm::atan2(..)` there passes clippy,
///   passes this scan, is genuinely portable — and is a transcendental in a
///   result path outside the audited surface with no wrapper contract.
///
/// **Two halves, because neither covers the other.** The manifest half is
/// airtight *between* crates and blind *inside* `borbax-units`; the source half
/// is the reverse. Cross-crate airtightness is not an argument, it is a
/// compiler property verified by the negative: with the manifest line removed,
/// the probe gives `error[E0433]: cannot find module or crate 'libm' in this
/// scope`. A transitive dependency is not nameable, so no crate can call
/// `libm::` without declaring it.
///
/// **`--locked` is a first line of defence and not a substitute**, which the
/// probes had to be rewritten to see: planting the manifest line alone makes
/// every `cargo` invocation refuse to build, so the guard never speaks. That is
/// the *build* failing, not the invariant holding — anyone adding the
/// dependency runs cargo once without `--locked`, commits the lockfile, and
/// from then on the tree is quiet. The probes below therefore update the
/// lockfile first, which is what a real change looks like.
///
/// **The obvious template carries a flaw that would be fatal here, and the
/// polarity is what makes it fatal.** [`check_wall_clock_has_one_home`] only
/// ever inspects files that are *not* the sanctioned home, so it cannot
/// distinguish one home from zero. For wall-clock that degenerates harmlessly —
/// zero clock reads is a fine state. For `libm` the meaning inverts: zero
/// `libm::` inside `borbax-units` means `det_math` has stopped calling libm,
/// which is exactly the state this guard exists to make impossible, and a
/// copied template would be green for it. Hence the positive assertion below,
/// which is `check_blocklist_present`'s "loud, not `Ok(())`" lesson applied to
/// a file rather than a constant.
///
/// **Why `> 0` and not a tighter floor.** Seven wrapper call sites were
/// measured on 2026-08-03, eight once `tan` landed — but the state being
/// defended is "the chokepoint still calls libm at all". A bar of eight would
/// fire on the deliberate removal of any wrapper, and a check that cries wolf
/// gets relaxed. The count is reported in the message so a collapse from eight
/// to one is visible to a reader even though it does not fail.
///
/// **Known residual holes, written down rather than fixed**, both dependency
/// reviews rather than gate business — the same disposition [`BANNED_TYPES`]
/// gives third-party hashers: a member vendoring libm's source under another
/// crate name, and `[patch]`/`[replace]` entries redirecting `libm` itself.
fn check_libm_has_one_home(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    // **Declarations, not uses.** This reads dependency tables; a crate that
    // reaches libm transitively is invisible to it, which is correct, because
    // a transitive dependency cannot be named in code.
    let mut declaring = Vec::new();
    let mut manifests_examined = 0usize;

    for scan_root in TRANSCENDENTAL_SCAN_ROOTS {
        let dir = root.join(scan_root);
        if !dir.exists() {
            failures.push(format!(
                "§13.1: scan root {scan_root:?} does not exist, so no manifest under it \
                 was checked for a `libm` dependency. Reported rather than skipped: a \
                 fail-open here disables the whole chokepoint check"
            ));
            continue;
        }
        for path in walk(&dir)?
            .into_iter()
            .filter(|p| p.file_name().is_some_and(|n| n == "Cargo.toml"))
        {
            manifests_examined += 1;
            let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
            // Matched as a whole token anywhere in the manifest rather than as a
            // key, so `mathlib = { package = "libm", .. }` is caught — a rename
            // still has to spell the crate name somewhere.
            if names_type(&toml_code_only(&src), "libm") {
                declaring.push(
                    path.strip_prefix(root)
                        .unwrap_or(&path)
                        .to_string_lossy()
                        .replace('\\', "/"),
                );
            }
        }
    }

    // **The corpus filter, asserted rather than assumed.** If the layout moves
    // and `walk` returns no manifests, every loop body above is skipped and the
    // check is silently green. Six members were measured on 2026-08-03; the
    // floor is stated as a floor because adding a crate must not fail this.
    if manifests_examined < 6 {
        failures.push(format!(
            "§13.1: only {manifests_examined} manifest(s) were examined under \
             {TRANSCENDENTAL_SCAN_ROOTS:?}, and 6 were measured on 2026-08-03. The \
             `libm` chokepoint check found almost nothing to check, which means the \
             workspace layout moved rather than that the invariant holds"
        ));
    }

    declaring.sort();
    if declaring != [LIBM_HOME_MANIFEST] {
        failures.push(format!(
            "§13.1: `libm` is declared by {declaring:?}, expected exactly \
             [\"{LIBM_HOME_MANIFEST}\"]. Transcendentals route through \
             `borbax_units::det_math` and nowhere else — a second declaration lets a \
             crate call `libm::acos` directly, which skips the clamp that wrapper \
             exists for, and silently takes the libm-bump determinism review out of \
             scope. If `borbax-units` no longer declares it, the chokepoint itself is \
             gone. Add `borbax-units` as a dependency and call `det_math` instead; if \
             the function you want has no wrapper, add one there"
        ));
    }

    // The source half: inside `borbax-units`, only `det_math.rs` may call libm.
    let units_src = root.join("crates/borbax-units/src");
    if !units_src.exists() {
        failures.push(
            "§13.1: crates/borbax-units/src does not exist, so the `libm` call-site \
             half was not checked. If the units crate has moved, move this check with \
             it"
            .to_owned(),
        );
        return Ok(());
    }

    let mut home_calls = 0usize;
    for path in walk(&units_src)?
        .into_iter()
        .filter(|p| p.extension().is_some_and(|e| e == "rs"))
    {
        let rel = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .to_string_lossy()
            .replace('\\', "/");
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        // **Shipped code only, and the probe is why.** The first version counted
        // `libm::` over comment-stripped text, which includes `#[cfg(test)]`
        // bodies — and `det_math.rs`'s own tests call `libm::cbrt` and
        // `libm::tan` to prove the wrappers are the portable ones. So rewriting
        // every wrapper to `x.cos()` left the count at 2 and the positive
        // assertion below stayed silent, in the one state it exists to forbid.
        // Its message said "no `libm::` call at all", which was true of what it
        // counted and false of what it claimed to defend. Found by running the
        // mutation rather than by reading the code.
        let Some(calls) = count_outside_tests(&src, LIBM_CALL) else {
            failures.push(format!(
                "§13.1: {rel} could not be lexed, so its `libm::` call sites were not \
                 counted. Reported rather than skipped: a file the lexer gives up on \
                 is the one most likely to be hiding something"
            ));
            continue;
        };
        if rel == LIBM_HOME_SOURCE {
            home_calls = calls;
        } else if calls > 0 {
            failures.push(format!(
                "§13.1: {calls} `libm::` call site(s) in {rel} — inside `borbax-units`, \
                 libm is called from `{LIBM_HOME_SOURCE}` alone. That file is what the \
                 libm-bump determinism review reads, and a call outside it is a \
                 transcendental nobody audits. It is also where the wrapper contracts \
                 live: `acos` clamps, and the raw one returns NaN on 17.88% of \
                 near-parallel unit-vector dot products"
            ));
        }
    }

    // **The positive assertion the wall-clock template lacks.** Without it, the
    // whole check passes in the one state it exists to forbid: `det_math`
    // rewritten to call `x.cos()` directly. The transcendental scan would catch
    // that too, and both firing is correct — this one names the cause.
    if home_calls == 0 {
        failures.push(format!(
            "§13.1: `{LIBM_HOME_SOURCE}` contains no `libm::` call at all. It is the \
             workspace's portable-transcendental chokepoint, so either it has been \
             rewritten to call the platform (which `check_no_platform_transcendentals` \
             should also be reporting) or it has moved and this check has not"
        ));
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

/// The root the fiction scans cover. **Unknown ⇒ covered.**
///
/// `crates/` and not the workspace: `experiments/` is deliberately out of
/// scope, declared here rather than skipped silently. It holds measurement
/// probes that ship to nobody and print to no child, and it carries the only
/// two hits the one-letter symbol tier finds on this tree (`"N"` in
/// `src/bin/fusion.rs` and `src/bin/ptable.rs`, both axis labels). `xtask`
/// itself is out for the reason [`UNSCANNED_CRATE_DIRS`] gives — it holds these
/// vocabularies as data, so scanning it would fire on every entry.
const FICTION_SCAN_ROOT: &str = "crates";

/// The one file allowed to contain the blocklist, by exact path.
///
/// It **must exist**: a missing or moved `naming.rs` is reported as a failure,
/// not treated as "nothing to exempt". That is the `check_blocklist_present`
/// lesson — a guard that reports success when its subject has vanished is worse
/// than no guard, because the green tick is read as evidence.
const BLOCKLIST_FILE: &str = "borbax-universe/src/naming.rs";

/// `xtask`'s own copy of the real element symbols (G2).
///
/// **Duplicated, never imported**, and the duplication is the point. `xtask`
/// depending on `borbax-universe` to enforce §5 inverts the gate: deleting the
/// blocklist would make the gate *fail to compile* rather than fail, and a gate
/// that goes down with the thing it guards is worse than one that reports red.
/// It would also make [`check_blocklist_present`] vacuous, since its whole job
/// is asserting the identifier still exists.
///
/// Held honest by `the_blocklist_the_guard_uses_is_the_blocklist_the_generator_uses`,
/// which extracts both arrays from `naming.rs` as text and compares them in both
/// directions — the same shape as the `clippy.toml` ↔ [`BANNED_CALLS`] pair.
const REAL_SYMBOLS: &[&str] = &[
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// `xtask`'s own copy of the real element names and chemical terms (G2).
///
/// Duplicated rather than imported, for the reason [`REAL_SYMBOLS`] gives.
const REAL_WORDS_COPY: &[&str] = &[
    "hydrogen",
    "helium",
    "lithium",
    "beryllium",
    "boron",
    "carbon",
    "nitrogen",
    "oxygen",
    "fluorine",
    "neon",
    "sodium",
    "magnesium",
    "aluminium",
    "aluminum",
    "silicon",
    "phosphorus",
    "sulfur",
    "sulphur",
    "chlorine",
    "argon",
    "potassium",
    "calcium",
    "scandium",
    "titanium",
    "vanadium",
    "chromium",
    "manganese",
    "iron",
    "cobalt",
    "nickel",
    "copper",
    "zinc",
    "gallium",
    "germanium",
    "arsenic",
    "selenium",
    "bromine",
    "krypton",
    "rubidium",
    "strontium",
    "yttrium",
    "zirconium",
    "niobium",
    "molybdenum",
    "technetium",
    "ruthenium",
    "rhodium",
    "palladium",
    "silver",
    "cadmium",
    "indium",
    "tin",
    "antimony",
    "tellurium",
    "iodine",
    "xenon",
    "caesium",
    "cesium",
    "barium",
    "lanthanum",
    "cerium",
    "tungsten",
    "platinum",
    "gold",
    "mercury",
    "thallium",
    "lead",
    "bismuth",
    "polonium",
    "astatine",
    "radon",
    "francium",
    "radium",
    "actinium",
    "thorium",
    "uranium",
    "neptunium",
    "plutonium",
    "water",
    "protein",
    "enzyme",
    "peptide",
    "lipid",
    "sugar",
    "glucose",
    "amine",
    "acid",
    "alcohol",
    "ester",
    "ketone",
    "methane",
    "ammonia",
    "benzene",
    "cellulose",
    "chitin",
];

/// Real-world unit spellings that may not appear in a string literal (G4).
///
/// **A deny-list over an open vocabulary cannot invert to fail-on-unknown**, and
/// saying so is part of shipping it: there is no finite set of real-world units,
/// so this catches the spellings anyone actually writes and nothing else. The
/// fail-on-unknown half lives in `borbax-ui`'s
/// `every_unit_word_on_screen_is_one_this_universe_invented`, an allow-list over
/// the *attained* set of unit words, which also catches one assembled at runtime
/// that no source scan can see. Neither half is sufficient; anyone proposing to
/// drop the acceptance test because "xtask covers it" is proposing to make the
/// guarantee pass-on-unknown.
///
/// Matched as case-insensitive substrings, so every entry must be distinctive
/// enough to survive that. The ones that are **not** on this list are the
/// interesting part, each measured against this tree:
///
/// - `second` — 3 hits (`"{first} {second}"`, `"a second layout pass"`).
/// - `pm`, `fm` — 2 and 1, and they are also *element symbols*, so they would
///   collide with the G2 fixtures.
/// - `year` — `world-years` is a **sanctioned** unit and contains it.
/// - `mol` alone — `Molecule` and `Mol12` are everywhere; it appears here only
///   inside the compound `kj/mol`, which is a literal.
/// - `metre`, `meter`, `gram`, `mole`, `litre`, `liter` — all fire as substrings
///   of ordinary words (`parameter`, `program`/`diagram`/`histogram`,
///   `molecule`, `obliterate`). They are matched as whole segments instead.
const REAL_UNIT_SUBSTRINGS: &[&str] = &[
    "angstrom",
    "kelvin",
    "celsius",
    "fahrenheit",
    "joule",
    "newton",
    "pascal",
    "dalton",
    "kilogram",
    "kj/mol",
    "kcal",
    "g/mol",
    "j/mol",
    "hertz",
    "\u{c5}",
    "\u{212b}",
    "\u{b0}c",
    "\u{b0}f",
    "\u{b5}m",
    // **The escaped spellings, because the scan reads source text and not
    // decoded strings.** A literal written `"\u{c5}"` renders as an angstrom
    // sign on screen and contains no `\u{c5}` character for a substring match to
    // find — measured, it walked past the first version of this list. Both
    // codepoints for the sign are listed: U+00C5 (Latin capital A with ring) and
    // U+212B (the angstrom sign proper), which look identical and are
    // interchangeable in practice.
    "\\u{c5}",
    "\\u{212b}",
    "\\u{b0}c",
    "\\u{b0}f",
    "\\u{b5}m",
];

/// Real-world units that are ordinary words inside other words, so they are
/// matched as whole segments rather than substrings. See [`REAL_UNIT_SUBSTRINGS`].
const REAL_UNIT_SEGMENTS: &[&str] = &[
    "metre", "metres", "meter", "meters", "gram", "grams", "mole", "moles", "litre", "litres",
    "liter", "liters", "amu", "ev", "nm",
];

/// G2 and G4 — no real element name, symbol or real-world unit in a literal.
///
/// **Neither check existed for any crate before Step 2**, and an earlier plan
/// said otherwise. [`check_blocklist_present`] asserts only that the
/// *identifier* `REAL_ELEMENT_SYMBOLS` still appears in one file; the blocklist
/// itself is a generation-time filter. No file was scanned for real element
/// names, and nothing at all was scanned for real-world units.
///
/// **The generator cannot mint a breach** — `naming::is_real` rejects and
/// redraws — so the subject here is the *display layer inventing a name the
/// generator did not produce*: a hard-coded `"Carbon"`, a `match units { 1 => "H" }`
/// table, a `"—"`-style placeholder that reaches for a real symbol.
///
/// **The rule is a split and the split is measured, not chosen.** Over this
/// tree, string literals only, comments stripped, `naming.rs` excluded:
///
/// | tier | match | false positives |
/// |---|---|---|
/// | real names | whole segment, case-insensitive | 0 |
/// | two-letter symbols | whole segment, case-**sensitive** | 0 |
/// | one-letter symbols | whole **literal**, case-sensitive | 0 |
/// | *(rejected)* any symbol, case-insensitive | | **442** |
///
/// The case asymmetry looks like the `elements.JSON` mistake and is not, so the
/// numbers stay here: `elements.JSON` is a *filesystem* lesson — on macOS and
/// Windows that is the same file, so the rename is invisible to whoever makes
/// it. Two string literals differing in case are two different tokens. And case
/// is semantic in the symbol vocabulary while meaningless in the name
/// vocabulary.
///
/// **Residual, stated so nobody reads this as total:** a hard-coded `"FE"` or
/// `"fe"` escapes, and so does a name assembled at runtime — `&name[..2]`
/// turning "helion" into `"He"` has no literal to scan. The second is covered by
/// `no_element_reaches_the_screen_under_a_real_name_or_symbol` in `borbax-ui`,
/// which asserts against `is_real` over what is actually displayed.
fn check_no_real_chemistry_in_literals(
    root: &Path,
    failures: &mut Vec<String>,
) -> Result<(), String> {
    let scan_root = root.join(FICTION_SCAN_ROOT);
    if !scan_root.exists() {
        failures.push(format!(
            "§5: {FICTION_SCAN_ROOT}/ does not exist, so no file was scanned for real \
             element names or real-world units. Reported rather than skipped — a \
             renamed root is how a guard goes quiet while every other check stays green"
        ));
        return Ok(());
    }
    let blocklist = scan_root.join(BLOCKLIST_FILE);
    if !blocklist.exists() {
        failures.push(format!(
            "§5: the blocklist at {FICTION_SCAN_ROOT}/{BLOCKLIST_FILE} is missing, so the \
             one file exempt from the name scan cannot be identified. If it moved, move \
             this check with it"
        ));
    }

    for path in walk(&scan_root)?
        .into_iter()
        .filter(|p| p.extension().is_some_and(|e| e == "rs"))
    {
        if path == blocklist || is_test_only_tree(&path) {
            continue;
        }
        let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        let rel = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .display()
            .to_string();
        let Ok(stream) = src.parse::<proc_macro2::TokenStream>() else {
            failures.push(format!(
                "§5: {rel} could not be lexed, so it was not scanned for real element \
                 names or units. Reported rather than skipped"
            ));
            continue;
        };
        let mut literals = Vec::new();
        collect_shipped_literals(stream, &mut literals);
        for lit in literals {
            for (word, why) in fiction_breaches(&lit) {
                failures.push(format!(
                    "§5: {rel} has `{word}` in the string literal {lit} — {why}. The viewer \
                     shows *generated* names and this universe's own units; a real one in a \
                     label is exactly as much a breach as one in a data file"
                ));
            }
        }
    }
    Ok(())
}

/// Whether a bracket group is a `#[cfg(..)]` that includes `test`.
///
/// **Not an equality test against `"cfg(test)"`, which is what this was.**
/// `#[cfg(all(test, feature = "x"))]` and `#[cfg(any(test, doc))]` are ordinary
/// spellings, and against a string comparison none of them matched — so the test
/// module underneath was scanned as shipped code and its fixtures, which
/// legitimately spell out real element names, fired the guard. A check that
/// fires on correct code gets deleted.
///
/// Deliberately narrow in the other direction: `test` must appear as a bare
/// **identifier**, so `#[cfg(feature = "test-utils")]` does not qualify — there
/// the word is inside a string literal and the item is shipped code.
fn names_cfg_test(stream: &proc_macro2::TokenStream) -> bool {
    let mut tokens = stream.clone().into_iter();
    let is_cfg = matches!(tokens.next(), Some(proc_macro2::TokenTree::Ident(i)) if i == "cfg");
    is_cfg && mentions_test_ident(stream)
}

/// Whether `test` appears as a bare identifier anywhere in `stream`, **outside
/// any `not(..)`**.
///
/// **`not(..)` inverts the predicate, so descending into it answers the opposite
/// question.** `#[cfg(not(test))]` marks *shipped-only* code; reading it as
/// test-only means the block beneath it is skipped, and shipped code going
/// unscanned is the one direction this guard must never fail in.
/// `#[cfg(all(not(test), unix))]` is the same shape.
///
/// **`cfg_test` in this same file already refuses `not` for exactly this
/// reason**, with a doc comment recording that its stated justification had
/// been wrong while its mutation test passed — and this function reintroduced
/// the defect 1900 lines away. Pinned here by
/// `cfg_not_test_is_not_read_as_a_test_module`.
fn mentions_test_ident(stream: &proc_macro2::TokenStream) -> bool {
    let mut it = stream.clone().into_iter();
    while let Some(tt) = it.next() {
        match tt {
            proc_macro2::TokenTree::Ident(i) if i == "test" => return true,
            // Skip the `not` and the parenthesised group that follows it.
            proc_macro2::TokenTree::Ident(i) if i == "not" => {
                it.next();
            }
            proc_macro2::TokenTree::Group(g) => {
                if mentions_test_ident(&g.stream()) {
                    return true;
                }
            }
            proc_macro2::TokenTree::Ident(_)
            | proc_macro2::TokenTree::Literal(_)
            | proc_macro2::TokenTree::Punct(_) => {}
        }
    }
    false
}

/// The token run that spells a `.pattern(..)` call.
///
/// **Named, so the test and the check cannot pass different needles.** The first
/// version of this was `[".", "pattern", "(", ")"]` — matching only an *empty*
/// argument list, so `pattern(&table)` scored nothing. The unit test caught
/// nothing when that was reintroduced, because it supplied its own needle
/// literal rather than this constant: the test pinned `count_outside_tests`, and
/// the defect was in what the caller handed it.
const PATTERN_CALL: &[&str] = &[".", "pattern", "("];

/// The token run that spells a `Universe::generate` call.
const GENERATE_CALL: &[&str] = &["Universe", ":", ":", "generate"];

/// The token run that spells an `embed(..)` call.
///
/// **The `(` is required and it is doing two jobs.** It keeps `use
/// borbax_molecule::embed;` — which is an import, not a call — from scoring, and
/// it keeps the needle from matching a bare mention. Identifiers lex whole, so
/// `re_embeds(` is one token `re_embeds` and does not match; that matters,
/// because the counter this guard's sibling reads is spelled exactly that.
const EMBED_CALL: &[&str] = &["embed", "("];

/// The token run that spells a `canonicalise(..)` call.
///
/// Matches the free function and a `.canonicalise(..)` method equally, which is
/// deliberate: `Signature::canonicalise` is also chemistry, and the viewer has
/// no more business calling one per frame than the other.
const CANONICALISE_CALL: &[&str] = &["canonicalise", "("];

/// Count how many times `needle` — a `::`-joined path or a `.method()` — appears
/// in `src` **outside** any `#[cfg(test)]` item.
///
/// **Attribute-aware, replacing a positional cut that was fail-open.** Both call
/// counts used to truncate the file at the first `#[cfg(test)]` and scan only
/// what came before, so anything written *after* a test module was invisible.
/// That was documented as a residual on the grounds that below a test module is
/// a strange place to hide a call — true, and still a hole a rename or a
/// re-ordering could open without anyone choosing to.
///
/// It is cheap now only because this PR already built the attribute pairing for
/// the fiction scan; before that it would have been new machinery. It closes the
/// *positional* half only. The **alias** half — `use borbax_universe::Universe
/// as U; U::generate(..)` — needs real symbol resolution and is deliberately
/// still open, named in the failure message and assigned to Step 7.
fn count_outside_tests(src: &str, needle: &[&str]) -> Option<usize> {
    let stream = src.parse::<proc_macro2::TokenStream>().ok()?;
    let mut flat = Vec::new();
    collect_shipped_tokens(stream, &mut flat);
    Some(
        flat.windows(needle.len())
            .filter(|w| w.iter().zip(needle).all(|(got, want)| got == want))
            .count(),
    )
}

/// The token *text* of `stream`, skipping `#[cfg(test)]` items and `#[doc]`
/// attributes, flattened so a path or a method call is a contiguous run.
fn collect_shipped_tokens(stream: proc_macro2::TokenStream, out: &mut Vec<String>) {
    let mut skip_next_body = false;
    for tt in stream {
        match tt {
            proc_macro2::TokenTree::Group(g) => {
                let is_bracket = matches!(g.delimiter(), proc_macro2::Delimiter::Bracket);
                let is_doc = is_bracket
                    && g.stream().into_iter().next().is_some_and(
                        |t| matches!(&t, proc_macro2::TokenTree::Ident(i) if i == "doc"),
                    );
                let is_cfg_test = is_bracket && names_cfg_test(&g.stream());
                let is_body = matches!(g.delimiter(), proc_macro2::Delimiter::Brace);
                if is_cfg_test {
                    skip_next_body = true;
                } else if is_body && skip_next_body {
                    skip_next_body = false;
                } else if !is_doc {
                    // A `()` after an ident is what makes `.pattern()` a call
                    // rather than a field, so the delimiter is recorded.
                    if matches!(g.delimiter(), proc_macro2::Delimiter::Parenthesis) {
                        out.push("(".to_owned());
                    }
                    collect_shipped_tokens(g.stream(), out);
                    if matches!(g.delimiter(), proc_macro2::Delimiter::Parenthesis) {
                        out.push(")".to_owned());
                    }
                }
            }
            proc_macro2::TokenTree::Ident(i) => out.push(i.to_string()),
            proc_macro2::TokenTree::Punct(p) => {
                if p.as_char() == ';' {
                    skip_next_body = false;
                }
                out.push(p.as_char().to_string());
            }
            proc_macro2::TokenTree::Literal(_) => out.push(String::new()),
        }
    }
}

/// Whether a path lives in a crate's integration-test, benchmark or example
/// tree.
///
/// **Skipped, because nothing in them reaches a screen and the alternative
/// fires on correct code.** A file under `crates/*/tests/` carries no
/// `#[cfg(test)]` attribute — the whole file *is* the test — so the attribute
/// pairing that exempts a unit-test module cannot see it, and a fixture
/// legitimately naming a real element would fail the gate. A check that fires
/// on correct code gets deleted, and it would take the name and unit scans with
/// it.
///
/// What covers the gap: the display path is asserted directly by
/// `no_element_reaches_the_screen_under_a_real_name_or_symbol` in `borbax-ui`,
/// which runs `naming::is_real` over what the panel actually shows — and that
/// catches a name *derived* at runtime, which no source scan can see either way.
fn is_test_only_tree(path: &Path) -> bool {
    path.components()
        .filter_map(|c| c.as_os_str().to_str())
        .any(|c| matches!(c, "tests" | "benches" | "examples"))
}

/// Every string literal in `stream`, skipping `#[doc]` groups and `#[cfg(test)]`
/// items.
///
/// The `#[doc]` skip is what makes the prose in this repository legal: these
/// files necessarily quote what they forbid. It is a **Bracket** test, which is
/// why `borbax-units`' unit descriptions had to become doc comments rather than
/// macro arguments — a bare literal inside `unit!(..)` sits in a *Parenthesis*
/// group and would not be skipped.
fn collect_shipped_literals(stream: proc_macro2::TokenStream, out: &mut Vec<String>) {
    // **Consulted for the next *brace* group, which is the bug the first draft
    // shipped.** The flag was set on seeing `#[cfg(test)]` and then only ever
    // cleared in the literal arm, so the `mod tests { .. }` body it was supposed
    // to skip was recursed into anyway. Nothing on the tree fired, because no
    // test currently spells a real name outside `naming.rs` — which is exactly
    // the shape of a guard that is wrong and looks right.
    let mut skip_next_body = false;
    for tt in stream {
        match tt {
            proc_macro2::TokenTree::Group(g) => {
                let mut inner = g.stream().into_iter().peekable();
                let is_bracket = matches!(g.delimiter(), proc_macro2::Delimiter::Bracket);
                let is_doc = is_bracket
                    && inner.peek().is_some_and(
                        |t| matches!(t, proc_macro2::TokenTree::Ident(i) if i == "doc"),
                    );
                let is_cfg_test = is_bracket && names_cfg_test(&g.stream());
                let is_body = matches!(g.delimiter(), proc_macro2::Delimiter::Brace);
                if is_cfg_test {
                    // The item this attribute decorates is test code, and a test
                    // fixture legitimately spells out what it forbids.
                    skip_next_body = true;
                } else if is_body && skip_next_body {
                    skip_next_body = false;
                } else if !is_doc {
                    // **The flag deliberately survives an attribute group.**
                    // `#[cfg(test)] #[allow(..)] mod tests { .. }` is ordinary,
                    // and clearing here made that module scan as shipped code —
                    // measured, its fixtures fired the guard. The terminator is
                    // what ends an item, and that is handled on `;` below.
                    collect_shipped_literals(g.stream(), out);
                }
            }
            proc_macro2::TokenTree::Literal(lit) => out.push(lit.to_string()),
            // **A `;` ends the item, and this is the whole fix.**
            // `#[cfg(test)]` is legal on an item that opens no brace —
            // `#[cfg(test)] use super::*;` or `#[cfg(test)] mod tests;` — and
            // with the flag cleared only by the next brace it survived to
            // swallow the next brace in the file, which can be *shipped* code:
            // the scan then reported nothing about it, silently.
            // `check_no_platform_transcendentals` carries the same rule for the
            // same reason (`a_cfg_test_item_that_opens_no_block_does_not_start_a_region`).
            //
            // Idents are the item keyword between the attribute and its body
            // (`mod`, `fn`) and attribute groups may sit there too, so neither
            // clears the flag.
            proc_macro2::TokenTree::Punct(p) if p.as_char() == ';' => skip_next_body = false,
            proc_macro2::TokenTree::Ident(_) | proc_macro2::TokenTree::Punct(_) => {}
        }
    }
}

/// A literal's whole words, split on non-alphanumerics only, case preserved.
///
/// **Deliberately *not* [`segments_preserving_case`], and the difference is a
/// measured false positive.** That splitter also breaks on camelCase, so `NaN`
/// becomes `["Na", "N"]` and sodium matches — six times on this tree, in
/// `packing.rs`, `layout.rs`, `geodesic.rs` and `borbax-units` itself, every one
/// of them a §13.4 comment about float handling. `Na` inside `NaN` is not
/// sodium, and a guard that fires on correct code gets deleted.
///
/// **The residual this buys, stated rather than hidden:** a two-letter symbol
/// glued to another word inside a single token is missed — `"NaCl"` scores
/// nothing. That is a real hole and it is the smaller one. The spelling a label
/// actually uses is a word: `"Sodium"` (caught by the name tier), `"Na"`
/// (caught here), `"carbon (C)"` (caught by the name tier). A hard-coded
/// `"NaCl"` in a UI label is caught by nothing and is worth a check of its own
/// the day formulas reach the screen — which is Step 5, where `BondError` text
/// first arrives.
fn whole_words_preserving_case(name: &str) -> Vec<String> {
    name.split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|w| !w.is_empty())
        .map(str::to_owned)
        .collect()
}

/// A string literal's content, with its prefix and delimiters removed.
///
/// **Positional, because `trim_matches` eats the content.** The first version of
/// this stripped `"`, `#`, `r` and `b` from both ends with `trim_matches`, which
/// removes *every* leading and trailing character satisfying the predicate — so
/// `"Nr"` became `N`, `"Ib"` became `I` and `"bH"` became `H`, and the
/// one-letter symbol tier reported nitrogen, iodine and hydrogen on literals
/// containing no symbol at all. That is the same defect
/// `a_literal_ending_in_r_is_not_a_format` already pins one scan over, and it
/// was reintroduced here by a reviewer's own recorded failure mode: reaching for
/// the convenient trim.
///
/// Handles `"…"`, `r"…"`, `r#"…"#`, `b"…"`, `br#"…"#` and `c"…"`. The trailing
/// delimiter is as long as the leading one, so the hash count is measured once
/// and reused.
fn literal_content(lit: &str) -> &str {
    let rest = lit.trim_start_matches(['r', 'b', 'c']);
    let hashes = rest.len() - rest.trim_start_matches('#').len();
    let rest = rest.get(hashes..).unwrap_or(rest);
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    let end = rest.len().saturating_sub(hashes + 1);
    rest.get(..end).unwrap_or(rest)
}

/// The G2/G4 breaches in one string literal's raw token text, with the reason.
fn fiction_breaches(lit: &str) -> Vec<(String, &'static str)> {
    let mut hits = Vec::new();
    if !lit.starts_with(['"', 'r', 'b', 'c']) {
        // Not a string literal — a number or a char. Nothing to say about it.
        return hits;
    }
    let lower = lit.to_ascii_lowercase();
    let content = literal_content(lit);

    // **The union of both splittings, and `cArBoN` is why.** `identifier_segments`
    // breaks on camelCase, so a mixed-case spelling shatters into `c|Ar|Bo|N` and
    // matches nothing — measured, it walked straight past this guard. Whole words
    // catch that; camelCase segments catch `HydrogenPeroxide`, which whole words
    // would miss. Neither alone is enough and the union costs nothing, because
    // names are long enough that both splittings measure clean on this tree.
    let name_forms: Vec<String> = identifier_segments(lit)
        .into_iter()
        .chain(
            whole_words_preserving_case(lit)
                .iter()
                .map(|w| w.to_ascii_lowercase()),
        )
        .collect();
    for word in REAL_WORDS_COPY {
        if name_forms.iter().any(|seg| seg == word) {
            hits.push((
                (*word).to_owned(),
                "a real element name or chemical term (G2)",
            ));
        }
    }
    // **Symbols match the whole literal, at either length.** The two-letter tier
    // used to match whole *words* inside a literal, which is unshippable:
    // measured, seven ordinary sentence-initial English words are element
    // symbols, so `"No universe loaded"` fired as nobelium, `"In this period"`
    // as indium, and `"At the top"`, `"Be careful"`, `"As shown above"`,
    // `"He typed a name"` and `"Am I right"` fired too. Those are exactly the
    // strings a UI writes, so this would have fired on correct code the first
    // time anyone wrote a sentence — and a check that fires on correct code gets
    // deleted, taking the name and unit scans with it.
    //
    // **What it costs, stated rather than hidden:** a symbol embedded in a
    // sentence — `"element Fe"` — escapes. The breach that actually happens is a
    // hard-coded table (`&["H", "He", "Li"]`, each a whole literal) or a real
    // *name*; the name tier still matches by word, because names do not collide
    // with English this way.
    for sym in REAL_SYMBOLS {
        if content == *sym {
            hits.push(((*sym).to_owned(), "a real element symbol (G2)"));
        }
    }
    for unit in REAL_UNIT_SUBSTRINGS {
        if lower.contains(unit) {
            hits.push(((*unit).to_owned(), "a real-world unit (G4)"));
        }
    }
    // **The union of both splittings, for the same reason the name tier needs
    // it.** `identifier_segments` breaks on camelCase, so `eV` becomes
    // `["e", "v"]` and `nM` becomes `["n", "m"]` — and those are the spellings
    // people actually write for electronvolts and nanometres. Whole words catch
    // them; camelCase segments still catch a unit glued into an identifier.
    let unit_forms: Vec<String> = identifier_segments(lit)
        .into_iter()
        .chain(
            whole_words_preserving_case(lit)
                .iter()
                .map(|w| w.to_ascii_lowercase()),
        )
        .collect();
    for unit in REAL_UNIT_SEGMENTS {
        if unit_forms.iter().any(|seg| seg == unit) {
            hits.push(((*unit).to_owned(), "a real-world unit (G4)"));
        }
    }
    hits
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

/// Does `code` name `ty` as a whole identifier token?
///
/// A bare `contains` fires on `DefaultHasherMetrics` and on any identifier with
/// a banned name as a prefix or suffix — measured, and a check that cries wolf
/// gets relaxed.
///
/// **Boundaries are tested as `char`s, not bytes**, which is the second
/// version. Bytes looked safe because every banned name is ASCII, but a Rust
/// identifier is not: `λDefaultHasher` is a legal identifier whose preceding
/// byte is the tail of a multi-byte sequence and therefore not
/// `is_ascii_alphanumeric`, so the byte version fired on it. A false positive
/// rather than a miss, and marginal — but the fix costs four characters and it
/// lets the doc above say "identifier boundary" instead of "ASCII identifier
/// boundary", which is the difference between a class and a habit.
fn names_type(code: &str, ty: &str) -> bool {
    let is_ident = |c: char| c.is_alphanumeric() || c == '_';
    code.match_indices(ty).any(|(at, _)| {
        let before_ok = !code
            .get(..at)
            .and_then(|head| head.chars().next_back())
            .is_some_and(is_ident);
        let after_ok = !code
            .get(at + ty.len()..)
            .and_then(|tail| tail.chars().next())
            .is_some_and(is_ident);
        before_ok && after_ok
    })
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
            for ty in BANNED_TYPES {
                if names_type(&code, ty) {
                    failures.push(format!(
                        "§13.1: unpinned hasher {ty} at {rel}:{} — its output is not stable \
                         across releases or across runs; derive seeds with \
                         borbax_rng::Stream (Domain::Hash)",
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
    use super::{GENERATE_CALL, PATTERN_CALL};
    use super::{REAL_SYMBOLS, REAL_WORDS_COPY, collect_shipped_literals, fiction_breaches};
    use super::{
        SIGNATURE_SURFACE, SigFacing, SigItem, alias_names, bookkeeping_against,
        check_surface_bookkeeping, compare_signature_surface, scan_signature_surface,
    };
    use super::{
        banned_idents_given, code_without_prose, count_outside_tests, viewer_banned_idents,
        viewer_banned_imports,
    };
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

    /// What the tests read off a scanned item.
    ///
    /// **Named fields, not a 5-tuple with three `bool`s.** `Where`'s own
    /// `#[expect(clippy::struct_excessive_bools)]` in this file argues that four
    /// positional `bool`s is the shape in which a slip compiles clean and
    /// silently inverts a category — and the first version of this helper
    /// returned exactly that, with the tests telling `is_pub` from `in_test` by
    /// counting underscores in a destructuring pattern.
    #[derive(Debug, PartialEq)]
    struct Seen {
        name: String,
        sees: usize,
        geodesic: bool,
        is_pub: bool,
        in_test: bool,
    }

    /// Mirrors the real pipeline: aliases are collected first, then the scan
    /// runs with the resulting name lists. Passing a bare `["Signature"]` here
    /// would make the alias rows of the corpus test something production does
    /// not do.
    fn surface(src: &str) -> Vec<Seen> {
        let file = syn::parse_file(src).unwrap_or_else(|e| unreachable!("probe must parse: {e}"));
        let mut sig = vec!["Signature".to_owned()];
        sig.extend(alias_names(&file.items, "Signature"));
        let mut geo = vec!["Geodesic".to_owned()];
        geo.extend(alias_names(&file.items, "Geodesic"));
        let (items, _) = scan_signature_surface("probe.rs", src, &sig, &geo)
            .unwrap_or_else(|e| unreachable!("probe must parse: {e}"));
        items
            .into_iter()
            .map(|i| Seen {
                name: i.name,
                sees: i.sees,
                geodesic: i.geodesic,
                is_pub: i.is_pub,
                in_test: i.in_test,
            })
            .collect()
    }

    /// Every way of putting two signatures in front of one function.
    ///
    /// **The corpus is the point, and its rows are not hypothetical.** The
    /// predecessor guard was deleted for missing four of them, and it had
    /// string-level tests that all passed — so a row is only worth having if it
    /// is a shape someone actually reached for. The first four are the four
    /// escapes review found; the rest are the shapes that closing them opened.
    #[test]
    fn every_way_of_seeing_two_signatures_is_counted_as_two() {
        for (label, src) in [
            (
                "named type, twice",
                "fn d(a: &Signature<12>, b: &Signature<12>) -> f64 { 0.0 }",
            ),
            (
                "receiver plus `&Self` — the per-argument counting bug",
                "impl Signature<12> { fn d(&self, o: &Self) -> f64 { 0.0 } }",
            ),
            (
                "receiver plus a named type",
                "impl Signature<12> { fn d(&self, o: &Signature<12>) -> f64 { 0.0 } }",
            ),
            (
                "a tuple holding both",
                "fn d(p: (Signature<12>, Signature<12>)) -> f64 { 0.0 }",
            ),
            (
                "a slice, which holds arbitrarily many",
                "fn d(s: &[Signature<12>]) -> f64 { 0.0 }",
            ),
            (
                "an array, likewise",
                "fn d(s: [Signature<12>; 4]) -> f64 { 0.0 }",
            ),
            (
                "a Vec, likewise",
                "fn d(s: Vec<Signature<12>>) -> f64 { 0.0 }",
            ),
            (
                "an iterator, likewise",
                "fn d(s: impl Iterator<Item = Signature<12>>) -> f64 { 0.0 }",
            ),
            (
                "a generic bound, which names no signature in any argument type",
                "fn d<T: AsRef<Signature<12>>>(a: T, b: T) -> f64 { 0.0 }",
            ),
            (
                "a generic bound through a where-clause spelling",
                "fn d<T>(a: T, b: T) -> f64 where T: AsRef<Signature<12>> { 0.0 }",
            ),
            (
                "a boxed slice — plural one indirection down",
                "fn d(s: Box<[Signature<12>]>) -> f64 { 0.0 }",
            ),
            (
                "an Rc'd slice, likewise",
                "fn d(s: Rc<[Signature<12>]>) -> f64 { 0.0 }",
            ),
            (
                "an alias, which defeats matching the bare identifier",
                "type Shape = Signature; fn d(a: &Shape, b: &Shape) -> f64 { 0.0 }",
            ),
            (
                "a renaming import, likewise",
                "use crate::sig::Signature as Shape; fn d(a: &Shape, b: &Shape) -> f64 { 0.0 }",
            ),
        ] {
            let found = surface(src);
            assert!(
                found.iter().any(|f| f.sees >= 2),
                "{label}: read as fewer than two signatures, so a comparison in it would pass \
                 as Single"
            );
        }
    }

    /// A `Geodesic` *receiver* is the spelling §8.3's kernel reads most
    /// naturally — `g.affinity(&a, &b)` — and treating it as ungrouped was one
    /// of the two false positives that got the predecessor deleted. A guard
    /// that fires on correct code gets deleted, so this is a correctness test
    /// for the guard, not a nicety.
    #[test]
    fn a_geodesic_receiver_counts_as_taking_one() {
        let found = surface(
            "impl Geodesic<12> { fn affinity(&self, a: &Signature<12>, b: &Signature<12>) -> f64 { 0.0 } }",
        );
        assert!(
            found
                .iter()
                .any(|f| f.name == "Geodesic::affinity" && f.sees >= 2 && f.geodesic),
            "a Geodesic receiver was not read as supplying the rotation group: {found:?}"
        );
    }

    /// Operands that reach a body through fields appear in no function
    /// signature at all — the one hole a function-only scan cannot see.
    #[test]
    fn signature_typed_fields_are_surface() {
        let found = surface("struct Pair { a: Signature<12>, b: Signature<12> }");
        assert!(
            found.iter().any(|f| f.name == "Pair.a") && found.iter().any(|f| f.name == "Pair.b"),
            "fields holding signatures were not recorded: {found:?}"
        );
    }

    /// `#[derive(PartialOrd)]` compares the stored lex-min arrays element by
    /// element. The compiler writes it, so there is no body for a reviewer to
    /// read — which is exactly why the derive list has to be *in* the pinned
    /// name rather than beside it.
    #[test]
    fn the_derive_list_is_part_of_the_pinned_name() {
        let before = surface("#[derive(Clone, PartialEq)] struct Signature { r: u8 }");
        let after = surface("#[derive(Clone, PartialEq, PartialOrd)] struct Signature { r: u8 }");
        assert_ne!(
            before, after,
            "adding PartialOrd left the pinned name unchanged, so it would land as an edit to a \
             listed entry rather than as an unlisted item"
        );
        // `assert_ne!` alone passes for any implementation whose output merely
        // *differs* — one that put a derive *count* in the name would satisfy
        // it. The name must contain the trait.
        assert!(
            after.iter().any(|f| f.name.contains("PartialOrd"))
                && !before.iter().any(|f| f.name.contains("PartialOrd")),
            "the pinned name changed without naming the trait that changed it: {after:?}"
        );
    }

    /// A derive can be written with a qualified path, and `Path::get_ident`
    /// returns `None` for every multi-segment one — so the trait vanished from
    /// the pinned name and the entry compared **equal** to the plain list.
    /// Measured against the real tree: `#[derive(Debug, Clone, PartialEq,
    /// ::core::cmp::PartialOrd)]` on `Signature` produced zero failures.
    #[test]
    fn a_qualified_derive_path_is_not_dropped() {
        let plain = surface("#[derive(Clone, PartialEq)] struct Signature { r: u8 }");
        for spelling in [
            "core::cmp::PartialOrd",
            "::core::cmp::PartialOrd",
            "std::cmp::PartialOrd",
        ] {
            let qualified = surface(&format!(
                "#[derive(Clone, PartialEq, {spelling})] struct Signature {{ r: u8 }}"
            ));
            assert_ne!(
                plain, qualified,
                "`{spelling}` left the pinned name unchanged, so an ordering over the stored \
                 lex-min arrays lands with the guard green"
            );
        }
    }

    /// `not(..)` inverts the predicate, so `#[cfg(not(test))]` means
    /// **shipped-only** — the one direction `Contrast`'s bound cannot survive.
    #[test]
    fn cfg_not_test_is_not_read_as_test_only() {
        for src in [
            "#[cfg(not(test))] fn d(a: &Signature<12>) {}",
            "#[cfg(all(not(test), unix))] fn d(a: &Signature<12>) {}",
            "#[cfg(not(test))] mod m { fn d(a: &Signature<12>) {} }",
        ] {
            let found = surface(src);
            assert!(
                !found.is_empty() && found.iter().all(|f| !f.in_test),
                "`cfg(not(test))` read as test-only, so a Contrast entry would ship in every \
                 release build: {src}"
            );
        }
    }

    /// Build a `SigItem` by named field, so a category test says which fact it
    /// is exercising.
    fn item_with(file: &str, name: &str, f: impl FnOnce(&mut SigItem)) -> SigItem {
        let mut i = SigItem {
            file: file.to_owned(),
            name: name.to_owned(),
            sees: 0,
            geodesic: false,
            is_pub: false,
            in_test: false,
            returns_ordering: false,
            derive: false,
            reads_raw: false,
            line: 1,
        };
        f(&mut i);
        i
    }

    /// **Every category's rejection, exercised directly.**
    ///
    /// The real-tree probes drive these through `SIGNATURE_SURFACE`, which
    /// means they can only test the categories the shipped list happens to use.
    /// A category whose check is never reached by the pinned list would be dead
    /// and look fine — this is the direct test of each arm.
    #[test]
    fn every_category_rejection_fires() {
        // Only the pinned-list lookup is exercised here, so the names are the
        // real ones; the *facts* are what vary.
        let cases: Vec<(&str, SigItem, &str)> = vec![
            (
                "Minimised with no Geodesic",
                item_with("crates/borbax-molecule/src/binding.rs", "fit", |i| {
                    i.sees = 2;
                    i.geodesic = false;
                }),
                "takes no Geodesic",
            ),
            (
                "Single that can see two",
                item_with(
                    "crates/borbax-molecule/src/binding.rs",
                    "mean_extent",
                    |i| {
                        i.sees = 2;
                    },
                ),
                "can see 2 signatures",
            ),
            (
                "Primitive that is reachable",
                item_with(
                    "crates/borbax-molecule/src/signature.rs",
                    "Signature::lex_cmp",
                    |i| {
                        i.sees = 2;
                        i.is_pub = true;
                    },
                ),
                "reachable outside its module",
            ),
            (
                "Contrast outside cfg(test)",
                item_with(
                    "crates/borbax-molecule/src/signature.rs",
                    "plain_distance",
                    |i| {
                        i.sees = 2;
                        i.in_test = false;
                    },
                ),
                "not inside `#[cfg(test)]`",
            ),
            (
                "Storage returning a number",
                item_with(
                    "crates/borbax-molecule/src/signature.rs",
                    "Signature: derive(Clone, Debug, PartialEq)",
                    |i| {
                        i.sees = 2;
                        i.returns_ordering = false;
                        i.derive = false;
                    },
                ),
                "does not return an ordering",
            ),
            (
                "Minimised that cannot see two",
                item_with("crates/borbax-molecule/src/binding.rs", "affinity", |i| {
                    i.sees = 1;
                    i.geodesic = true;
                }),
                "claims it can see two signatures",
            ),
        ];
        for (label, item, wanted) in cases {
            let mut failures = Vec::new();
            compare_signature_surface(&[item], &mut failures);
            assert!(
                failures.iter().any(|f| f.contains(wanted)),
                "{label}: no failure mentioned {wanted:?} — that category's check is dead. \
                 Got: {failures:?}"
            );
        }
    }

    /// A derive list has no arguments to count and no visibility of its own, so
    /// the `Primitive` and `Contrast` bounds are *unmeasurable* on it and would
    /// pass vacuously.
    #[test]
    fn a_derive_cannot_claim_a_bounded_category() {
        let mut failures = Vec::new();
        compare_signature_surface(
            &[item_with(
                "crates/borbax-molecule/src/signature.rs",
                "Signature::lex_cmp",
                |i| {
                    i.derive = true;
                },
            )],
            &mut failures,
        );
        assert!(
            failures.iter().any(|f| f.contains("pass vacuously")),
            "a derive pinned as a bounded category was not refused: {failures:?}"
        );
    }

    /// Both bookkeeping directions, driven directly rather than through the
    /// real tree — where a duplicate key cannot be planted without also
    /// planting the duplicate item.
    #[test]
    fn the_bookkeeping_counts_rather_than_searches() {
        let mut failures = Vec::new();
        check_surface_bookkeeping(
            &[
                item_with("crates/borbax-molecule/src/binding.rs", "fit", |_| {}),
                item_with("crates/borbax-molecule/src/binding.rs", "fit", |_| {}),
            ],
            &mut failures,
        );
        assert!(
            failures.iter().any(|f| f.contains("share the key")),
            "two items sharing a pinned key went unreported — one category would silently \
             govern both: {failures:?}"
        );
        // The duplicate-PIN direction, which no real-tree plant can reach: you
        // cannot plant a repeated list entry without editing the shipped list.
        // The branch carries a `take(i).all(..)` guard so a repeated key is
        // reported exactly once, and nothing exercised it.
        let mut failures = Vec::new();
        let doubled = &[
            ("a.rs", "d", SigFacing::Single),
            ("a.rs", "d", SigFacing::Minimised),
        ];
        bookkeeping_against(&[item_with("a.rs", "d", |_| {})], doubled, &mut failures);
        let dupes = failures.iter().filter(|f| f.contains("2 times")).count();
        assert_eq!(
            dupes, 1,
            "a repeated pinned key was reported {dupes} times, not once — the `take(i)` guard \
             is the thing that makes it exactly one: {failures:?}"
        );
        // And the other direction: everything pinned but absent.
        let mut failures = Vec::new();
        check_surface_bookkeeping(&[], &mut failures);
        assert!(
            failures.len() >= SIGNATURE_SURFACE.len(),
            "an empty tree did not report every pinned entry as stale: {} of {}",
            failures.len(),
            SIGNATURE_SURFACE.len()
        );
    }

    /// The predecessor was deleted partly for firing on correct code, and the
    /// canonical example is the one §8.2's own module doc asks for: an ordering
    /// on `Signature` so a `BTreeMap` can intern species. Measured before
    /// `Storage` existed, **all five** other categories refuted it.
    #[test]
    fn an_ordering_on_the_stored_form_is_classifiable() {
        let found = surface(
            "impl Signature<12> { fn partial_cmp(&self, o: &Self) -> Option<Ordering> { None } }",
        );
        assert!(
            found.iter().any(|f| f.sees >= 2 && !f.geodesic),
            "the ordering was not even seen: {found:?}"
        );
        // Storage's bound is the return type, so a distance cannot be filed
        // under it however it is named.
        let (items, _) = scan_signature_surface(
            "probe.rs",
            "impl Signature<12> { fn partial_cmp(&self, o: &Self) -> f64 { 0.0 } }",
            &["Signature".to_owned()],
            &["Geodesic".to_owned()],
        )
        .unwrap_or_else(|e| unreachable!("probe must parse: {e}"));
        assert!(
            items.iter().all(|i| !i.returns_ordering),
            "a function returning f64 read as an ordering, so Storage's bound is inert"
        );
        let (items, _) = scan_signature_surface(
            "probe.rs",
            "impl Signature<12> { fn eq(&self, o: &Self) -> bool { true } }",
            &["Signature".to_owned()],
            &["Geodesic".to_owned()],
        )
        .unwrap_or_else(|e| unreachable!("probe must parse: {e}"));
        assert!(
            items.iter().all(|i| i.returns_ordering),
            "a `-> bool` equality did not read as an ordering, so Storage would reject the \
             `PartialEq` §8.2 says is fine"
        );
    }

    /// A body that reads the stored arrays is surface even when nothing in its
    /// declaration names a signature — the intern-key shape, measured green
    /// before this dimension existed. Shipped code only: tests reading
    /// `a.extents()[i]` in assertions are deliberately not listed.
    #[test]
    fn a_body_reading_the_stored_arrays_is_surface() {
        let holder = "struct Key { sig: Signature<12> }
             impl Key { pub fn cmp2(&self, o: &Self) -> f64 { self.sig.extents()[0].get() - o.sig.extents()[0].get() } }";
        let found = surface(holder);
        assert!(
            found.iter().any(|f| f.name == "Key::cmp2"),
            "a holder method comparing through fields went unrecorded: {found:?}"
        );
        let in_test = format!("#[cfg(test)] mod t {{ {holder} }}");
        let found = surface(&in_test);
        assert!(
            !found.iter().any(|f| f.name == "Key::cmp2"),
            "a test body was listed on the accessor dimension — that is the noise that gets a \
             guard switched off: {found:?}"
        );
    }

    /// The two bounded categories are bounded by facts the scan reads, not by
    /// the author's word: `Primitive` must be private and `Contrast` must be
    /// under `#[cfg(test)]`. Without these, either name is a way to label an
    /// ungrouped comparison and move on.
    #[test]
    fn the_escape_hatch_categories_carry_their_bounds() {
        // **Both directions for both flags.** A positive-only corpus is
        // satisfied maximally by the mutation that matters: `is_pub = true`
        // everywhere makes `Primitive` unusable but never *wrong*, and
        // `in_test = true` everywhere makes `Contrast`'s bound pass for shipped
        // code — the exact hole this category was given a bound to close.
        for (label, src, want_pub) in [
            ("`pub fn`", "mod m { pub fn d(a: &Signature<12>) {} }", true),
            (
                "`pub(crate) fn` — reaches every consumer in the crate",
                "mod m { pub(crate) fn d(a: &Signature<12>) {} }",
                true,
            ),
            (
                "`pub(super) fn`",
                "mod m { pub(super) fn d(a: &Signature<12>) {} }",
                true,
            ),
            (
                "a bare private `fn`",
                "mod m { fn d(a: &Signature<12>) {} }",
                false,
            ),
        ] {
            let found = surface(src);
            assert!(
                found.iter().all(|f| f.is_pub == want_pub),
                "{label} read as is_pub={}, wanted {want_pub} — Primitive's bound is only as \
                 good as this flag",
                !want_pub
            );
        }
        for (label, src, want_test) in [
            (
                "`#[cfg(test)] mod`",
                "#[cfg(test)] mod t { fn d(a: &Signature<12>) {} }",
                true,
            ),
            (
                "`cfg(all(test, ..))`, which nests one level deeper",
                "#[cfg(all(test, feature = \"x\"))] mod t { fn d(a: &Signature<12>) {} }",
                true,
            ),
            (
                "`#[cfg(test)]` on the fn itself",
                "#[cfg(test)] fn d(a: &Signature<12>) {}",
                true,
            ),
            (
                "a plain module",
                "mod m { fn d(a: &Signature<12>) {} }",
                false,
            ),
            (
                "`cfg(feature = ..)`, which is not `cfg(test)`",
                "#[cfg(feature = \"x\")] fn d(a: &Signature<12>) {}",
                false,
            ),
        ] {
            let found = surface(src);
            assert!(
                found.iter().all(|f| f.in_test == want_test),
                "{label} read as in_test={}, wanted {want_test} — Contrast's bound is only as \
                 good as this flag",
                !want_test
            );
        }
    }

    /// The counting must be exact in both directions. A corpus that only
    /// asserts `>= 2` is satisfied maximally by `sees_count` returning 2 for
    /// everything, which would fail every `Single` entry in the pinned list and
    /// look like the guard working.
    #[test]
    fn one_signature_is_counted_as_one() {
        for (label, src) in [
            (
                "one named parameter",
                "fn d(a: &Signature<12>) -> f64 { 0.0 }",
            ),
            (
                "a receiver alone",
                "impl Signature<12> { fn m(&self) -> f64 { 0.0 } }",
            ),
            (
                "one parameter plus a Geodesic",
                "fn d(a: &Signature<12>, g: &Geodesic<12>) -> f64 { 0.0 }",
            ),
            (
                "an array of something else",
                "fn d(a: &Signature<12>, r: &[Span; 12]) -> f64 { 0.0 }",
            ),
        ] {
            let found = surface(src);
            assert!(
                !found.is_empty() && found.iter().all(|f| f.sees <= 1),
                "{label}: read as seeing two, so a Single entry would fail and the pinned list \
                 would have to lie: {found:?}"
            );
        }
    }

    /// `any` short-circuits. The first version of the enum arm used it, so
    /// every variant after the first signature-bearing one went unvisited —
    /// measured on `enum Holder { First(Signature<D>), Second(Signature<D>) }`,
    /// where only `Holder::First.0` was reported.
    #[test]
    fn every_enum_variant_is_visited_not_just_the_first() {
        let found = surface("enum Holder { First(Signature<12>), Second(Signature<12>) }");
        assert!(
            found.iter().any(|f| f.name == "Holder::First.0")
                && found.iter().any(|f| f.name == "Holder::Second.0"),
            "a later variant's signature field went unrecorded — the walker stopped walking: \
             {found:?}"
        );
    }

    /// An *associated* const is the same hole one nesting level in, and both
    /// `impl` and `trait` had it falling to `_ => {}`.
    #[test]
    fn an_associated_const_is_surface() {
        let found = surface("impl Holder { const PAIR: [Signature<12>; 2] = []; }");
        assert!(
            found
                .iter()
                .any(|f| f.name == "Holder::PAIR" && f.sees >= 2),
            "an associated const holding signatures went unrecorded: {found:?}"
        );
        let found = surface("trait Keyed { const PAIR: [Signature<12>; 2]; }");
        assert!(
            found.iter().any(|f| f.name == "Keyed::PAIR" && f.sees >= 2),
            "a trait const holding signatures went unrecorded: {found:?}"
        );
    }

    /// A type that *is* a signature holds one by being one. Seeding `holds`
    /// from variant fields alone meant `enum Signature { .. }`'s derive list
    /// never reached the pinned name.
    #[test]
    fn a_signature_enums_own_derive_list_is_recorded() {
        let found = surface("#[derive(Clone, PartialOrd)] enum Signature { A(u8) }");
        assert!(
            found
                .iter()
                .any(|f| f.name.contains("derive") && f.name.contains("PartialOrd")),
            "an enum named Signature did not record its own derive list: {found:?}"
        );
    }

    /// A `const` or `static` puts both operands in scope of every body in the
    /// module while appearing in no function signature and no struct field.
    #[test]
    fn a_static_holding_signatures_is_surface() {
        let found = surface("static PAIR: [Signature<12>; 2] = [];");
        assert!(
            found.iter().any(|f| f.name == "PAIR" && f.sees >= 2),
            "a static array of signatures was not recorded as a two-signature site: {found:?}"
        );
    }

    /// Prose about a type must not make its file look like a consumer of one.
    /// `borbax-units/src/lib.rs` mentions `Signature.r` in a doc comment, and
    /// reading that as a mention reported five constructs there as unreadable —
    /// four `unit!` invocations and the `macro_rules!` definition — noise that
    /// gets a guard switched off.
    #[test]
    fn a_type_named_only_in_prose_does_not_make_a_file_a_consumer() {
        let (items, unreadable) = scan_signature_surface(
            "probe.rs",
            "//! `Signature.r` was on that list.\nmacro_rules! m { () => {} }\n",
            &["Signature".to_owned()],
            &["Geodesic".to_owned()],
        )
        .unwrap_or_else(|e| unreachable!("probe must parse: {e}"));
        assert!(
            items.is_empty() && unreadable.is_empty(),
            "a doc comment made the file look like a signature consumer: {items:?} {unreadable:?}"
        );
    }

    /// Silence must mean "looked and it was clean". A macro in a file that does
    /// name `Signature` is code this guard cannot read, and it says so.
    #[test]
    fn an_unreadable_construct_in_a_consumer_file_is_reported() {
        let (_, unreadable) = scan_signature_surface(
            "probe.rs",
            "fn f(a: &Signature<12>) {}\nmacro_rules! m { () => {} }\n",
            &["Signature".to_owned()],
            &["Geodesic".to_owned()],
        )
        .unwrap_or_else(|e| unreachable!("probe must parse: {e}"));
        assert!(
            !unreadable.is_empty(),
            "a macro went unreported in a file that names Signature"
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

    use super::{BANNED_CALLS, BANNED_TYPES, names_type, scan_rust_source};

    fn scan(src: &str) -> Vec<String> {
        let mut failures = Vec::new();
        scan_rust_source("t.rs", src, &mut failures);
        failures
    }

    /// The `path = ".."` values inside one named array of `clippy.toml`.
    ///
    /// **Parsed, not scanned, and that is a repair with a measurement behind
    /// it.** The first version walked lines looking for `path = "`, tracking
    /// which array it was inside by `starts_with`. A review lane fed it seven
    /// TOML shapes and it mis-read three:
    ///
    /// - a **commented-out** entry was returned as if it were active — the
    ///   silent direction, on a §13.1 guard. Someone disabling a hasher ban to
    ///   unblock a build would leave this cross-check green while clippy, the
    ///   half this test calls "the authority", enforced nothing;
    /// - a single-line `key = [ .. ]` never reset the section flag, so the
    ///   *types* array's paths were returned as the *methods* list — exactly
    ///   the failure the section-awareness was added to fix, reachable again by
    ///   a different trigger;
    /// - `disallowed-types-extra` was swallowed into `disallowed-types` by
    ///   prefix match.
    ///
    /// Each has a patch and the class does not, which is the argument
    /// `xtask/Cargo.toml` already makes for parsing Rust with `syn` rather than
    /// grepping it. The class this now decides is "whatever the TOML spec says
    /// an array of tables with a `path` key is", and nothing is outside it.
    fn clippy_paths_in(src: &str, key: &str) -> Vec<String> {
        let doc: toml::Table = src
            .parse()
            .unwrap_or_else(|e| unreachable!("clippy.toml does not parse: {e}"));
        let Some(entries) = doc.get(key).and_then(toml::Value::as_array) else {
            return Vec::new();
        };
        entries
            .iter()
            .map(|entry| {
                // **Both forms, because clippy accepts both.** A reviewer
                // measured it rather than trusting the lint-configuration page,
                // which says bare strings are invalid: clippy 1.97.1 fires on
                // `disallowed-types = [ "std::hash::RandomState" ]` identically
                // to the table form. `filter_map(|e| e.get("path"))` dropped
                // those silently, which is the direction that matters — a type
                // clippy bans would be invisible to the textual scan forever,
                // and the cross-check that exists to notice would say nothing.
                entry
                    .as_str()
                    .or_else(|| entry.get("path").and_then(toml::Value::as_str))
                    .unwrap_or_else(|| {
                        unreachable!(
                            "a `{key}` entry is neither a path string nor a table with a \
                             `path` key: {entry:?}"
                        )
                    })
                    .to_owned()
            })
            .collect()
    }

    fn workspace_cargo_toml() -> String {
        let toml = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .map(|p| p.join("Cargo.toml"))
            .and_then(|p| std::fs::read_to_string(p).ok())
            .unwrap_or_default();
        assert!(!toml.is_empty(), "workspace Cargo.toml not found or empty");
        toml
    }

    /// §13.1's two ban lists are still **denied**, not merely configured.
    ///
    /// **The third state, and neither existing guard can see it.** The
    /// `#[expect(clippy::disallowed_types)]` anchors in `borbax-universe` cannot:
    /// `#[expect]` is itself a lint-level attribute and overrides the outer
    /// level inside its own scope, so the lint still fires there and the
    /// expectation stays fulfilled. The cross-check below cannot: it compares
    /// path lists and knows nothing about levels. Three review lanes measured
    /// the same escape independently — flipping one word in the workspace
    /// manifest to `"allow"` and planting a real `DefaultHasher` gives clippy
    /// exit 0 with every anchor green.
    ///
    /// `clippy.toml` holds the *list*; the manifest holds the *level*, and a
    /// list at `warn` is advisory. `disallowed_types` in particular is
    /// warn-by-default, so before it was denied the ban rested entirely on CI
    /// passing `-D warnings`.
    #[test]
    fn the_ban_lists_are_denied_and_not_merely_configured() {
        let manifest: toml::Table = workspace_cargo_toml()
            .parse()
            .unwrap_or_else(|e| unreachable!("workspace Cargo.toml does not parse: {e}"));
        let clippy_lints = manifest
            .get("workspace")
            .and_then(|w| w.get("lints"))
            .and_then(|l| l.get("clippy"));
        for lint in ["disallowed_types", "disallowed_methods"] {
            // **Both spellings Cargo accepts.** `lint = "deny"` and
            // `lint = { level = "deny", priority = 1 }` are the same thing, and
            // the second is the form Cargo *asks* for when group priorities
            // conflict — so a legal edit that keeps the deny would have failed
            // this test. It fails closed, but it is the same one-form blindness
            // `clippy_paths_in` was repaired for, in the check written to
            // repair it.
            //
            // This landed a commit later than its own commit message said it
            // did: the scripted edit that was meant to apply it reported a miss
            // and the miss was not chased. A reviewer read the diff against the
            // claim and caught it, which is the recorded remedy working.
            let configured = clippy_lints.and_then(|c| c.get(lint)).and_then(|entry| {
                entry
                    .as_str()
                    .or_else(|| entry.get("level").and_then(toml::Value::as_str))
            });
            assert_eq!(
                configured,
                Some("deny"),
                "[workspace.lints.clippy] {lint} is not `deny` — §13.1's ban is \
                 advisory, and nothing else can see this: the #[expect] anchors \
                 set the level for their own items, and the path cross-check \
                 knows nothing about levels"
            );
        }
    }

    /// Clippy lint **group** names, as of 1.97.1.
    ///
    /// A hand-kept closed list, which this codebase distrusts elsewhere and for
    /// good reason — it covers every group clippy has today and would miss one
    /// a future release invents. The enumeration-free alternative is to stop
    /// reading the manifest and *ask clippy*, with a `trybuild`-style fixture
    /// that must be rejected. That is stronger and it is real work; this is the
    /// cheap version, and the gap is written down rather than left implied.
    const CLIPPY_GROUPS: &[&str] = &[
        "all",
        "correctness",
        "suspicious",
        "style",
        "complexity",
        "perf",
        "pedantic",
        "nursery",
        "restriction",
        "cargo",
        "deprecated",
    ];

    /// No lint **group** outranks the specific denies beneath it.
    ///
    /// **The fourth state, and the third guard in a row defeated by its
    /// neighbour rather than by anything wrong with itself.**
    /// `the_ban_lists_are_denied_and_not_merely_configured` asserts the string
    /// `"deny"` is written down. Cargo applies lower priority first, so a
    /// sibling group at priority >= 0 lands *after* every bare
    /// `lint = "deny"` (which is priority 0) and silences it — with the deny
    /// still written down verbatim and that test still green.
    ///
    /// Measured on the real tree at 1.97.1: adding
    /// `all = { level = "allow", priority = 1 }` and planting an
    /// `f64::total_cmp` and a `DefaultHasher` in shipped library code took
    /// clippy's disallowed-diagnostic count from **3 to 0**, left all 78 xtask
    /// tests green, and left all three `#[expect]` anchors fulfilled.
    ///
    /// **This one reaches a result, unlike the hasher half.** `CLIPPY_ONLY`
    /// names `total_cmp` as enforced by clippy *alone*, deliberately — the text
    /// scan cannot express a per-site `#[expect]`. Its own reason string says
    /// it "orders on the sign bit, and a runtime NaN's sign differs by
    /// architecture", which is precisely a §13.4 cross-platform divergence with
    /// no second guard behind it.
    ///
    /// `warn` counts as outranking: these bans are `deny` because a warning is
    /// advisory, and a group demoting them to `warn` at priority >= 0 is the
    /// same defect one step smaller.
    ///
    /// Measured, so the list below is a fact rather than a guess:
    /// `disallowed_types` and `disallowed_methods` live in the `style` group,
    /// so `all` and `style` reach them and `pedantic`, `nursery`,
    /// `restriction` and `cargo` do not. The manifest's own three groups sit at
    /// priority -1, which is why it is safe today.
    #[test]
    fn no_lint_group_outranks_the_specific_denies() {
        let manifest: toml::Table = workspace_cargo_toml()
            .parse()
            .unwrap_or_else(|e| unreachable!("workspace Cargo.toml does not parse: {e}"));
        let clippy = manifest
            .get("workspace")
            .and_then(|w| w.get("lints"))
            .and_then(|l| l.get("clippy"))
            .and_then(toml::Value::as_table)
            .unwrap_or_else(|| unreachable!("[workspace.lints.clippy] is missing"));

        // **The floor the specific denies actually sit at, not the literal 0.**
        // The first version compared every group against zero, which is only
        // right while the denies carry cargo's default priority. Give them a
        // negative priority and a group at a *higher* negative priority still
        // lands after them: `disallowed_types = { level = "deny", priority = -5 }`
        // beside `all = { level = "allow", priority = -1 }` silences both bans
        // while this test and its sibling stay green. Measured — clippy's
        // disallowed-diagnostic count went to 0 with an `f64::total_cmp` and a
        // `DefaultHasher` planted in shipped library code. The machine reviewer
        // found it, one step down the number line from the escape this test was
        // written for.
        let deny_priority = ["disallowed_types", "disallowed_methods"]
            .iter()
            .filter_map(|lint| clippy.get(*lint))
            .map(|value| {
                value
                    .get("priority")
                    .and_then(toml::Value::as_integer)
                    .unwrap_or(0)
            })
            .min()
            .unwrap_or_else(|| unreachable!("[workspace.lints.clippy] has neither specific deny"));

        for (name, value) in clippy {
            if !CLIPPY_GROUPS.contains(&name.as_str()) {
                continue;
            }
            // A bare string level carries cargo's default priority, which is 0.
            let priority = value
                .get("priority")
                .and_then(toml::Value::as_integer)
                .unwrap_or(0);
            let level = value
                .as_str()
                .or_else(|| value.get("level").and_then(toml::Value::as_str))
                .unwrap_or_else(|| unreachable!("[workspace.lints.clippy] {name} has no level"));
            assert!(
                level == "deny" || level == "forbid" || priority < deny_priority,
                "[workspace.lints.clippy] group `{name}` is `{level}` at priority \
                 {priority}, and the specific denies sit at {deny_priority}. Cargo \
                 applies lower priority first, so this group lands *after* them and \
                 silences them — including §13.4's `f64::total_cmp` ban, which the \
                 text scan deliberately cannot enforce. A group is only safe at a \
                 strictly lower priority than the denies it can reach"
            );
        }
    }

    /// `.cargo/config.toml` does not cap the lints the gate depends on.
    ///
    /// `[build] rustflags = ["--cap-lints=allow"]` silences every §13.1 and
    /// §13.4 clippy ban — measured. (The narrower
    /// `rustflags = ["-Aclippy::disallowed_types"]` does *not*, also measured,
    /// so this is specifically the `--cap-lints` form.)
    ///
    /// Lower severity than the group escape above because `--cap-lints=allow`
    /// also kills `unsafe_code = "forbid"` and every other workspace lint,
    /// which makes it conspicuous in a diff — where a lint-group priority reads
    /// like ordinary config tidying. Checked anyway because this file is
    /// committed, and the check costs a line.
    #[test]
    fn the_cargo_config_does_not_cap_lints() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .map(|p| p.join(".cargo/config.toml"));
        let Some(src) = path.and_then(|p| std::fs::read_to_string(p).ok()) else {
            return;
        };
        let config: toml::Table = src
            .parse()
            .unwrap_or_else(|e| unreachable!(".cargo/config.toml does not parse: {e}"));
        assert!(
            config
                .get("build")
                .and_then(|b| b.get("rustflags"))
                .is_none(),
            ".cargo/config.toml sets `build.rustflags`. `--cap-lints=allow` there \
             silences every §13.1 and §13.4 clippy ban while the whole gate stays \
             green; if a rustflag is genuinely wanted, this check needs to learn \
             which ones are safe rather than being deleted"
        );

        // **`target.<triple>.rustflags` and `target.<cfg>.rustflags` too**, and
        // this is the half worth having: the §13.4 matrix pins three targets, so
        // a target-scoped cap silences the bans on *one leg* and leaves the other
        // two green. A divergence that appears on one platform and is
        // unenforced on that platform is precisely the shape the golden matrix
        // exists to catch and would then misreport.
        if let Some(targets) = config.get("target").and_then(toml::Value::as_table) {
            for (selector, spec) in targets {
                assert!(
                    spec.get("rustflags").is_none(),
                    ".cargo/config.toml sets `target.{selector}.rustflags`, which caps \
                     lints for that target alone — leaving the other legs of the §13.4 \
                     matrix enforcing a ban this one does not"
                );
            }
        }
    }

    /// No ban entry is silenced with `allow-invalid`.
    ///
    /// Clippy warns when a configured path does not resolve, and offers
    /// `allow-invalid = true` to suppress it — measured, that is clippy's own
    /// help text. The warning is a config diagnostic with no lint level, so
    /// `-D warnings` never promoted it anyway; what `allow-invalid` removes is
    /// the last visible sign that an entry is enforcing nothing. Two reviewers
    /// reached for it independently while probing, which is how likely it is to
    /// be pasted in for real.
    #[test]
    fn no_ban_entry_is_silenced_with_allow_invalid() {
        let doc: toml::Table = clippy_toml()
            .parse()
            .unwrap_or_else(|e| unreachable!("clippy.toml does not parse: {e}"));
        for key in ["disallowed-types", "disallowed-methods"] {
            let Some(entries) = doc.get(key).and_then(toml::Value::as_array) else {
                continue;
            };
            for entry in entries {
                assert!(
                    entry.get("allow-invalid").is_none(),
                    "a `{key}` entry carries `allow-invalid`, which silences the one \
                     diagnostic that says the path stopped resolving: {entry:?}"
                );
            }
        }
    }

    fn clippy_toml() -> String {
        let toml = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .map(|p| p.join("clippy.toml"))
            .and_then(|p| std::fs::read_to_string(p).ok())
            .unwrap_or_default();
        assert!(!toml.is_empty(), "clippy.toml not found or empty");
        toml
    }

    /// §13.1's hasher ban, held in step across its two enforcement points.
    ///
    /// Same argument as `the_two_ban_lists_cover_the_same_functions`: two
    /// lists are only redundant while they agree, and nothing else keeps them
    /// in step. The match is on the **last path segment**, because
    /// `clippy.toml` names a resolvable path and [`BANNED_TYPES`] names the
    /// identifier as written — see that constant for why those differ.
    ///
    /// **The empty-list check is a better diagnostic, not a load-bearing
    /// assertion, and an earlier version of this comment claimed otherwise.**
    /// It said both directions below "are vacuously satisfied by two empty
    /// lists", so deleting `disallowed-types` wholesale would pass without it.
    /// Measured false by a review lane: the `orphaned` direction iterates
    /// `BANNED_TYPES`, which is never empty, so it fires either way. The assert
    /// earns its place by naming the actual cause instead of listing three
    /// orphans, and the commit message that shipped the false version is
    /// corrected here rather than left standing.
    #[test]
    fn the_two_type_ban_lists_cover_the_same_types() {
        let toml = clippy_toml();
        let paths = clippy_paths_in(&toml, "disallowed-types");
        assert!(
            !paths.is_empty(),
            "clippy.toml has no `disallowed-types` array — §13.1's hasher ban \
             has no authoritative half, and the textual list below cannot see \
             a `use ... as` rename or a macro expansion on its own"
        );

        let mut missing = Vec::new();
        for path in &paths {
            let leaf = path.rsplit("::").next().unwrap_or(path);
            if !BANNED_TYPES.contains(&leaf) {
                missing.push(path.clone());
            }
        }
        assert!(
            missing.is_empty(),
            "in clippy.toml but not BANNED_TYPES: {missing:?} — the text scan is \
             the only thing that would catch one of these being deleted from \
             clippy.toml, so it has to know about them"
        );

        let mut orphaned = Vec::new();
        for ty in BANNED_TYPES {
            if TEXT_SCAN_ONLY.contains(ty) {
                continue;
            }
            if !paths
                .iter()
                .any(|p| p.rsplit("::").next().unwrap_or(p) == *ty)
            {
                orphaned.push((*ty).to_owned());
            }
        }
        assert!(
            orphaned.is_empty(),
            "in BANNED_TYPES but not clippy.toml: {orphaned:?} — clippy is the \
             authority, so a type missing there is unenforced against every \
             spelling a text scan cannot see"
        );

        // An exemption that exempts something nobody bans enforces nothing, and
        // would say so to no one. Today `a_banned_hasher_is_matched_as_a_whole_token`
        // happens to cover it through a fixture line; this makes the coverage
        // intent rather than coincidence.
        for exempt in TEXT_SCAN_ONLY {
            assert!(
                BANNED_TYPES.contains(exempt),
                "{exempt:?} is exempted from the clippy half but is not banned at \
                 all — TEXT_SCAN_ONLY names types the text scan enforces alone, \
                 not types nothing enforces"
            );
        }
    }

    /// A banned name is matched as a token, not as a substring.
    ///
    /// **Both directions, or the "fix" is indistinguishable from deleting the
    /// check.** The machine reviewer found the substring behaviour; the negative cases
    /// are what stop the repair from over-correcting into silence.
    #[test]
    fn a_banned_hasher_is_matched_as_a_whole_token() {
        // Fires: the bare name, a qualified path (`::` is not an identifier
        // character), and a generic argument.
        for names in [
            "DefaultHasher::new()",
            "ahash::RandomState::default()",
            "BuildHasherDefault<DefaultHasher>",
            "let h: SipHasher;",
            "rustc_hash::FxRandomState::default()",
        ] {
            assert!(
                BANNED_TYPES.iter().any(|ty| names_type(names, ty)),
                "{names:?} should be reported"
            );
        }
        // Silent: a banned name as a prefix or a suffix of a longer identifier.
        for quiet in [
            "DefaultHasherMetrics",
            "MyDefaultHasher",
            "RandomStateBuilder",
            "a_RandomState_thing",
            "make_random_state",
            // Non-ASCII identifier characters are identifier characters. The
            // byte-boundary version fired on both of these.
            "\u{3bb}DefaultHasher",
            "DefaultHasher\u{3bb}",
        ] {
            for ty in BANNED_TYPES {
                assert!(!names_type(quiet, ty), "{quiet:?} should NOT be reported");
            }
        }
    }

    /// `clippy.toml` is parsed, so a disabled entry reads as disabled.
    ///
    /// The line-based predecessor returned a **commented-out** entry as active
    /// — the silent direction on a §13.1 guard — and lost its section on a
    /// single-line array. Both shapes are pinned here.
    #[test]
    fn a_commented_out_ban_is_not_read_as_an_active_one() {
        let src = "\
disallowed-types = [
  { path = \"a::Alpha\", reason = \"live\" },
  # { path = \"b::Beta\", reason = \"temporarily disabled\" },
]
";
        assert_eq!(clippy_paths_in(src, "disallowed-types"), vec!["a::Alpha"]);

        // A one-line array must not leak into the next key's list.
        let one_line = "\
disallowed-methods = [ { path = \"f64::exp\" } ]
disallowed-types = [ { path = \"x::Y\" } ]
";
        assert_eq!(
            clippy_paths_in(one_line, "disallowed-methods"),
            vec!["f64::exp"]
        );
        assert_eq!(clippy_paths_in(one_line, "disallowed-types"), vec!["x::Y"]);

        // A key that merely starts with the one asked for is a different key.
        let lookalike = "\
disallowed-types-extra = [ { path = \"wrong::One\" } ]
disallowed-types = [ { path = \"right::One\" } ]
";
        assert_eq!(
            clippy_paths_in(lookalike, "disallowed-types"),
            vec!["right::One"]
        );
    }

    /// The hasher ban catches the spelling anyone would actually write.
    ///
    /// Probed rather than asserted: the mutation this exists to stop is
    /// `seed_from_phrase` drafted with a `DefaultHasher`, which passed all six
    /// gate legs before this guard landed.
    #[test]
    fn reports_an_unpinned_hasher() {
        let found = scan("fn f() -> DefaultHasher { DefaultHasher::new() }\n");
        assert!(
            found.iter().any(|f| f.contains("DefaultHasher")),
            "{found:?}"
        );
        assert!(
            scan("fn f() -> RandomState { RandomState::new() }\n")
                .iter()
                .any(|f| f.contains("RandomState"))
        );
        // Imported and written unqualified — the spelling a path form misses,
        // and the reason `BANNED_TYPES` holds bare identifiers.
        assert!(
            scan("use std::hash::DefaultHasher;\n")
                .iter()
                .any(|f| f.contains("DefaultHasher")),
        );
        // The exemption that makes the golden's own negative control writable:
        // a test may name a hasher in order to prove we do not agree with it.
        assert!(scan("#[cfg(test)]\nmod t {\n    use std::hash::DefaultHasher;\n}\n").is_empty());
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

    /// Enforced by the **text scan alone**, because clippy has no path to
    /// resolve.
    ///
    /// `rustc-hash` is not a workspace dependency, so a `disallowed-types`
    /// entry for `FxRandomState` would be inert *and* would make clippy warn on
    /// every run — silenceable only with `allow-invalid`, which is the very
    /// switch `no_ban_entry_is_silenced_with_allow_invalid` forbids. Same shape
    /// as [`BANNED_PARALLEL_CALLS`] and for the reason that constant already
    /// gives: a text scan is the only check available until the day it would be
    /// too late to add one. Move it the day `rustc-hash` becomes a direct
    /// dependency.
    const TEXT_SCAN_ONLY: &[&str] = &["FxRandomState"];

    /// Enforced by clippy alone, because they have legitimate library-code
    /// uses that need a per-site `#[expect]` — which the text scan cannot
    /// express, and deliberately so.
    const CLIPPY_ONLY: &[&str] = &["total_cmp"];

    #[test]
    fn the_two_ban_lists_cover_the_same_functions() {
        let clippy_toml = clippy_toml();
        let method_paths = clippy_paths_in(&clippy_toml, "disallowed-methods");
        // As with the type cross-check: a better diagnostic, not a load-bearing
        // assertion. The `orphaned` direction iterates `BANNED_CALLS`, which is
        // never empty, so it fires either way. An earlier version of this
        // message claimed both directions were vacuous over an empty list —
        // the same false claim this commit corrected one screen above, written
        // fresh into its sibling. A reviewer caught it.
        assert!(
            !method_paths.is_empty(),
            "clippy.toml has no `disallowed-methods` array"
        );

        let mut missing = Vec::new();
        for path in &method_paths {
            let Some((_, func)) = path.split_once("::") else {
                continue;
            };
            let path = path.as_str();
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
                && !method_paths
                    .iter()
                    .any(|p| p == &format!("{width}::{func}"))
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

    /// An unnamed viewer file gets the **strictest** tier, not the laxest.
    ///
    /// This is the whole inversion. The shipped check was a two-entry table, so
    /// a new `.rs` file under `crates/borbax-ui/src` was checked for *nothing* —
    /// it could name the engine, or build strings while calling itself a drawing
    /// file. Step 2 is the first step that wanted to add a file, which is why
    /// this is the moment it was worth fixing.
    #[test]
    fn a_viewer_file_nobody_named_is_covered_rather_than_exempt() {
        for unnamed in ["periodic.rs", "elements.rs", "scene/mesh.rs", "lib.rs"] {
            assert_eq!(
                viewer_banned_imports(unnamed),
                ["egui"],
                "{unnamed} was not given the default substring tier, so a file added \
                 tomorrow is unguarded"
            );
            // **The other half of the tier, and it has to be asserted
            // separately.** `bevy` moved to the identifier-matched list when the
            // engine changed, because `bevy_egui` — which the drawing tier is
            // allowed — contains `bevy` as a substring. A test that only checked
            // `viewer_banned_imports` would report the default tier intact while
            // the engine ban had quietly left it.
            assert!(
                viewer_banned_idents(unnamed).contains(&"bevy"),
                "{unnamed} may name the engine without being the shell"
            );
        }
    }

    /// Being a drawing file is a **licence with a price**: `egui` in exchange
    /// for `format!`.
    ///
    /// The pair matters. A drawing file that could also build strings would put
    /// a string where it is painted, and every test asserting on that string
    /// would be describing a window it had stopped describing.
    #[test]
    fn a_named_drawing_file_trades_format_for_egui() {
        let banned = viewer_banned_imports("panel.rs");
        assert!(
            !banned.contains(&"egui"),
            "a drawing file must be allowed `egui` — it is what `egui_kittest` drives, \
             and since the move to Bevy it arrives as `bevy_egui::egui`"
        );
        assert!(
            banned.contains(&"format!"),
            "a drawing file must not build strings; that is the price of the licence"
        );
        assert!(
            viewer_banned_idents("panel.rs").contains(&"bevy"),
            "only the shell may name the engine — a drawing file gets `bevy_egui` and \
             not `bevy`, which is exactly why this needle is matched at identifier \
             boundaries rather than as a substring"
        );
    }

    /// The engine files trade `bevy` and `egui` for `format!`, and nothing is
    /// unrestricted.
    ///
    /// **The shell used to be exempt from everything and no longer is**, which
    /// is a tightening rather than a rename. Its old licence was justified by
    /// "no headless test reaches it, so there is nothing to protect" — and that
    /// stopped being true the moment the app construction moved into `app.rs`
    /// so `tests/app.rs` could reach it. A file a test asserts about is a file
    /// where a `format!` can hide a string the test cannot see.
    #[test]
    fn the_engine_files_are_licensed_and_not_exempt() {
        for engine in super::VIEWER_ENGINE_FILES {
            let banned = viewer_banned_imports(engine);
            assert!(
                !banned.contains(&"egui"),
                "{engine} wires the bridge and needs `egui` to build the root `Ui`"
            );
            assert!(
                !viewer_banned_idents(engine).contains(&"bevy"),
                "{engine} is where the engine lives"
            );
            assert!(
                banned.contains(&"format!"),
                "{engine} must not build strings — that rule is about where a \
                 string is made, not about which crates a file may import"
            );
        }
        // The negative arm. Without it this passes over a
        // `viewer_banned_imports` that returns `&[]` for everything.
        assert!(
            viewer_banned_idents("state.rs").contains(&"bevy"),
            "every other file is restricted; a licence everywhere is the vacuous \
             version of this guard"
        );
    }

    /// `wgpu` and `epaint` are refused by default, exactly as `egui` is.
    ///
    /// **Both were holes on `main`, found by reading the tier table rather than
    /// by a test.** `viewer_banned_imports` returns `["egui"]` for the
    /// strictest tier, and raw `wgpu` matches neither — so `state.rs`, the file
    /// whose whole identity is "a test can reach it without a window", could
    /// write `use wgpu::Device;` and pass. `egui_wgpu` *was* caught, by the
    /// substring `egui`, and that asymmetry is what made the hole look closed.
    #[test]
    fn an_unnamed_viewer_file_may_name_neither_wgpu_nor_epaint() {
        for unnamed in ["state.rs", "panel.rs", "scene/pipeline.rs", "lib.rs"] {
            let banned = viewer_banned_idents(unnamed);
            assert!(
                banned.contains(&"wgpu"),
                "{unnamed} may name `wgpu` without holding the licence"
            );
            assert!(
                banned.contains(&"epaint"),
                "{unnamed} may name `epaint`, which is `egui`'s painting layer and \
                 invisible to the `egui` needle"
            );
        }
    }

    /// The licence list is the *only* thing that lifts the `wgpu` ban.
    ///
    /// Written against a synthetic list because the real one is empty, which
    /// leaves the exemption arm unreachable — see [`super::banned_idents_given`].
    /// The negative arm is the one that carries the information: without it this
    /// passes on a function that lifts the ban for everybody.
    #[test]
    fn only_a_licensed_file_may_name_wgpu() {
        let gpu = ["scene/draw.rs"];
        assert!(
            !banned_idents_given("scene/draw.rs", &gpu).contains(&"wgpu"),
            "a listed file must be allowed `wgpu`; that is what the list is for"
        );
        assert!(
            banned_idents_given("scene/camera.rs", &gpu).contains(&"wgpu"),
            "an unlisted file next door must still be refused — unknown means strict, \
             and this arm is what catches the inverted condition"
        );
        assert!(
            banned_idents_given("scene/draw.rs", &gpu).contains(&"epaint"),
            "the GPU licence covers `wgpu` and nothing else"
        );
    }

    /// No file holds the GPU licence yet, which is why the test above is
    /// written against a synthetic list.
    ///
    /// It is not a placeholder: it fails the day a file is added to
    /// [`super::VIEWER_GPU_FILES`] without this test being read, which is the
    /// moment someone should be checking that the entry was deliberate and that
    /// `only_a_licensed_file_may_name_wgpu` still exercises both arms.
    #[test]
    fn the_gpu_licence_list_is_empty_until_a_file_earns_it() {
        assert!(
            super::VIEWER_GPU_FILES.is_empty(),
            "{:?} now holds the GPU licence — check the entry was deliberate, and that \
             the seam tests still cover both arms of the exemption",
            super::VIEWER_GPU_FILES
        );
    }

    /// `request_repaint` is not `epaint`, and `egui_wgpu` is not `wgpu`.
    ///
    /// **The near-miss is live, not hypothetical.** `panel.rs` contains the word
    /// `repaints` in prose today, and `ctx.request_repaint()` is ordinary code a
    /// continuously-animating scene writes — a substring `epaint` needle fires on
    /// it, on correct code, in the tier that file is in. Guards that cry wolf get
    /// deleted; that is why these two are matched at identifier boundaries while
    /// `egui` deliberately is not.
    #[test]
    fn the_new_needles_do_not_fire_on_the_words_that_contain_them() {
        assert!(
            !names_type("ctx.request_repaint();", "epaint"),
            "`request_repaint` is correct code and must not read as `epaint`"
        );
        assert!(
            !names_type("use egui_wgpu::Renderer;", "wgpu"),
            "`egui_wgpu` is caught by the `egui` needle, not by this one — counting it \
             twice makes the two rules indistinguishable"
        );
        // The positive arms, without which the two above pass on a matcher that
        // never fires at all.
        assert!(names_type("use epaint::Color32;", "epaint"));
        assert!(names_type("use wgpu::Device;", "wgpu"));
    }

    /// Every way a manifest can declare a dependency is seen, including the
    /// two that hide the crate's real name.
    ///
    /// **The long form and the rename are the whole reason this is parsed
    /// rather than grepped.** `borbax-ui` already uses `[dependencies.eframe]`,
    /// so a scan that only read `key = value` lines under `[dependencies]` would
    /// miss the largest dependency in the workspace. And
    /// `renderer = { package = "bevy" }` declares `bevy` under a key no amount
    /// of key-matching would catch.
    #[test]
    fn a_dependency_is_found_however_it_is_spelled() {
        let manifest = "\
[package]
name = \"borbax-nope\"
version = \"0.1.0\"

[dependencies]
thiserror = { workspace = true }
renderer = { package = \"bevy\", version = \"0.19\" }
plain = \"1.0\"

[dependencies.eframe]
version = \"0.35\"
default-features = false

[dev-dependencies.masked]
package = \"bevy_egui\"
version = \"0.41\"

[target.'cfg(unix)'.dependencies]
winit = \"0.30\"

[profile.dev]
opt-level = 1
";
        let found = super::manifest_dependency_names(manifest);
        for want in ["thiserror", "bevy", "plain", "eframe", "bevy_egui", "winit"] {
            assert!(
                found.iter().any(|d| d == want),
                "{want} was not found in {found:?}"
            );
        }
        // The negative arms. Without these the test passes on a parser that
        // returns every token in the file.
        assert!(
            !found.iter().any(|d| d == "renderer"),
            "the rename key was returned instead of the crate it renames: {found:?}"
        );
        assert!(
            !found.iter().any(|d| d == "masked"),
            "the long-form rename key survived: {found:?}"
        );
        assert!(
            !found.iter().any(|d| d == "opt-level" || d == "name"),
            "a non-dependency table was read as dependencies: {found:?}"
        );
    }

    /// A type named inside a string literal is prose, not an import.
    ///
    /// **Measured on the real tree, and it failed the seam on correct code.**
    /// `crates/borbax-ui/src/lib.rs` carries an `#![expect]` attribute whose
    /// `reason` string names `eframe` while explaining why the lint fires —
    /// documentation, and exactly
    /// as much prose as the `//!` block above it. The weaker `code_only` strip
    /// reads it as code. A guard that fires on correct code gets deleted, so
    /// this is not a nicety.
    #[test]
    fn a_type_named_in_a_string_literal_is_not_an_import() {
        let src = "#![expect(clippy::x, reason = \"`eframe` emits 31 of these\")]\nfn f() {}\n";
        let stripped = code_without_prose(src);
        // Asserted before defaulting, deliberately. `unwrap_or_default()` on a
        // `None` gives `""`, which contains no needle and would make every
        // assertion below pass over a fixture that had stopped lexing.
        assert!(stripped.is_some(), "the fixture must lex");
        let code = stripped.unwrap_or_default();
        assert!(
            !code.contains("eframe"),
            "prose inside an attribute string was read as an import: {code:?}"
        );
    }

    /// …and the macro name survives the literal it builds.
    ///
    /// The reverse arm, and it is what keeps the strip from disarming the
    /// `format!` half of the seam. `format!` sits *outside* the string, so
    /// removing string contents must not remove it — otherwise `panel.rs` could
    /// build every string it liked.
    #[test]
    fn stripping_a_literal_does_not_strip_the_macro_that_builds_it() {
        let stripped = code_without_prose("fn f() -> String { format!(\"eframe {}\", 1) }\n");
        assert!(stripped.is_some(), "the fixture must lex");
        let code = stripped.unwrap_or_default();
        assert!(
            code.contains("format!"),
            "the `format!` half of the seam was disarmed by the literal strip: {code:?}"
        );
        assert!(
            !code.contains("eframe"),
            "the literal's *contents* should still be gone: {code:?}"
        );
    }

    /// A file that ends mid-string is reported, never silently read as empty.
    ///
    /// The failure mode this forbids is the quiet one: an unlexable file
    /// stripped to nothing contains no banned needle, so it *passes* every tier
    /// while being the file most likely to be hiding something.
    #[test]
    fn a_file_that_ends_inside_a_string_is_reported_rather_than_scanned() {
        assert!(
            code_without_prose("fn f() { let _ = \"unterminated;\n").is_none(),
            "an unterminated string must not strip to a clean empty scan"
        );
        assert!(
            code_without_prose("fn f() { /* unterminated\n").is_none(),
            "an unterminated block comment must not strip to a clean empty scan"
        );
        // The positive arm, so this cannot pass by `code_without_prose` always
        // returning `None`.
        assert!(
            code_without_prose("fn f() { let _ = \"closed\"; }\n").is_some(),
            "ordinary code must lex"
        );
    }

    /// Every G2/G4 breach found in one source fixture, by the real scan path.
    fn fiction_hits(src: &str) -> Vec<String> {
        let stream = src
            .parse::<proc_macro2::TokenStream>()
            .unwrap_or_else(|_| proc_macro2::TokenStream::new());
        // **Asserted, because a fixture that does not lex strips to nothing.**
        // An empty stream contains no vocabulary entry, so every `is_empty()`
        // assertion below would pass without scanning anything — the vacuous
        // shape this repository keeps re-finding. `g5` asserts the same thing
        // for the same reason.
        assert!(
            !src.trim().is_empty() && !stream.is_empty(),
            "the fixture did not lex, so this assertion would pass over nothing: {src}"
        );
        let mut lits = Vec::new();
        collect_shipped_literals(stream, &mut lits);
        lits.iter()
            .flat_map(|l| fiction_breaches(l))
            .map(|(w, _)| w)
            .collect()
    }

    /// A real element name hard-coded into a display string is caught.
    ///
    /// The generator cannot mint one — `naming::is_real` rejects and redraws —
    /// so the subject is the display layer inventing a name the generator did
    /// not produce.
    #[test]
    fn a_hardcoded_real_element_name_in_a_display_string_is_caught() {
        for src in [
            r#"fn f() -> &'static str { "Carbon" }"#,
            r#"fn f() -> &'static str { "carbon (C)" }"#,
            r#"const T: &[&str] = &["Hydrogen", "Helium"];"#,
            r#"fn f() -> &'static str { "made of water" }"#,
        ] {
            assert!(
                !fiction_hits(src).is_empty(),
                "a real name reached the screen unnoticed: {src}"
            );
        }
    }

    /// Real symbols are caught in the spellings a person actually writes.
    ///
    /// Two-letter symbols as whole words; one-letter symbols only as the entire
    /// literal, because `"N"` is a maths label far more often than nitrogen.
    #[test]
    fn a_hardcoded_real_symbol_is_caught_in_the_spellings_a_person_writes() {
        for src in [
            r#"fn f() -> &'static str { "He" }"#,
            r#"fn f() -> &'static str { "Fe" }"#,
            r#"fn f() -> &'static str { "C" }"#,
            r#"const T: &[&str] = &["H", "He", "Li"];"#,
        ] {
            assert!(
                !fiction_hits(src).is_empty(),
                "a real symbol reached the screen unnoticed: {src}"
            );
        }
        // The negative arm, and it is what keeps the one-letter tier from being
        // the 442-false-positive version. A single letter *inside* a sentence is
        // not an element.
        assert!(
            fiction_hits(r#"fn f() -> &'static str { "N of M" }"#).is_empty(),
            "a one-letter symbol must only match as the WHOLE literal"
        );
    }

    /// Names fold case; symbols do not. Both halves asserted.
    ///
    /// **The numbers are the reason and they belong here**, because this
    /// asymmetry looks like the `elements.JSON` mistake to anyone applying that
    /// lesson uniformly: a case-insensitive symbol rule gives **442** false
    /// positives on this tree (`n`×93, `k`×50, `at`×44, `i`×40, `u`×32) against
    /// **0** for the case-sensitive one. `elements.JSON` is a *filesystem*
    /// lesson — those are the same file on macOS, so the rename is invisible.
    /// Two string literals differing in case are two different tokens.
    #[test]
    fn a_case_flipped_real_name_is_still_caught_and_a_case_flipped_symbol_is_not() {
        for spelling in ["Carbon", "carbon", "CARBON", "cArBoN"] {
            assert!(
                !fiction_hits(&format!("fn f() -> &'static str {{ \"{spelling}\" }}")).is_empty(),
                "case must not matter for a name: {spelling}"
            );
        }
        for spelling in ["fe", "FE", "nO", "aT", "iN", "no", "at", "in"] {
            assert!(
                fiction_hits(&format!("fn f() -> &'static str {{ \"{spelling} x\" }}")).is_empty(),
                "case MUST matter for a symbol, or this guard fires 442 times: {spelling}"
            );
        }
    }

    /// `Na` inside `NaN` is not sodium.
    ///
    /// **Measured, and it is why the symbol tier splits on whole words rather
    /// than reusing [`super::identifier_segments`]**, which also breaks on
    /// camelCase and turns `NaN` into `["Na", "N"]`. Six real files fired —
    /// `packing.rs`, `layout.rs`, `geodesic.rs` and `borbax-units` itself — every
    /// one of them a §13.4 comment about float handling. A guard that fires on
    /// correct code gets deleted.
    // The fixtures below quote real assertion messages verbatim, braces and all,
    // because the point is that this guard sees what is actually written in the
    // tree. They are source text handed to a lexer, never a format string.
    #[expect(
        clippy::literal_string_with_formatting_args,
        reason = "fixture text copied from the real tree, lexed rather than formatted"
    )]
    #[test]
    fn nan_is_not_sodium() {
        for src in [
            r#"fn f() { assert!(x, "NaN did not sort last: {v:?}"); }"#,
            r#"fn f() { assert!(x, "unexpected NaN payload: 0x{:016x}"); }"#,
            r#"fn f() { assert!(x, "a NaN pivot must fail the factorisation"); }"#,
        ] {
            assert!(
                fiction_hits(src).is_empty(),
                "the float spelling NaN was read as sodium: {src}"
            );
        }
    }

    /// A format placeholder and an ordinary English word are not units.
    ///
    /// Each of these was measured on the tree and each is why the token beside
    /// it is absent from the vocabulary: `second` (3 hits), `pm`/`fm` (2 and 1,
    /// and both are also element symbols), `mol` alone (`Molecule`, `Mol12`),
    /// `metre`/`gram`/`mole`/`liter` inside `parameter`, `histogram`,
    /// `molecule`, `obliterate`. And `world-years` contains `year`, so `year`
    /// can never be a bare entry.
    #[test]
    fn a_format_placeholder_is_not_a_unit() {
        for src in [
            r#"fn f() -> &'static str { "{first} {second}" }"#,
            r#"fn f() -> &'static str { "a second layout pass" }"#,
            r#"fn f() -> &'static str { "parameter histogram molecule" }"#,
            r#"fn f() -> &'static str { "obliterate the diagram" }"#,
            r#"fn f() -> &'static str { "12 world-years" }"#,
            r#"fn f() -> &'static str { "0.92 spans" }"#,
            r#"fn f() -> &'static str { "2.60 quanta per unit" }"#,
        ] {
            assert!(
                fiction_hits(src).is_empty(),
                "an invented unit or ordinary word was read as a real-world unit: {src}"
            );
        }
    }

    /// A real unit is caught in the spellings that matter.
    #[test]
    fn a_real_world_unit_is_caught_in_every_spelling_that_matters() {
        for src in [
            r#"fn f() -> &'static str { "1.15 \u{c5}" }"#,
            r#"fn f() -> &'static str { "2.4 kJ/mol" }"#,
            r#"fn f() -> &'static str { "204 Kelvin" }"#,
            r#"fn f() -> &'static str { "12 nm" }"#,
            r#"fn f() -> &'static str { "70.7 amu" }"#,
            r#"fn f() -> &'static str { "3 grams" }"#,
            r#"fn f() -> &'static str { "5 metres" }"#,
        ] {
            assert!(
                !fiction_hits(src).is_empty(),
                "a real-world unit reached a label unnoticed: {src}"
            );
        }
    }

    /// The disclaimers are doc comments, and a bare literal is **not** exempt.
    ///
    /// **This is the test that decides how G4 is enforced**, and it is silent
    /// under either fix that was rejected. `borbax-units` necessarily names the
    /// units it disclaims ("Not Kelvin, not Celsius"); as a bare `$doc:literal`
    /// those words sit in a *Parenthesis* group, which the `#[doc]` skip — a
    /// **Bracket** test — cannot reach, so the crate that exists to disclaim
    /// real units would fail the guard that enforces them.
    ///
    /// A path exemption for that file is a hand-kept list of one, and it would
    /// silence the first fixture below. A `Not `-prefix rule would silence the
    /// first fixture *and* `"Not Å"` in any UI label — a guard satisfiable by
    /// editing the string it fired on trains people to edit strings. Taking
    /// `$(#[$doc:meta])*` instead means the words arrive as `#[doc = "…"]` and
    /// the existing skip covers them, with zero new exemption surface.
    #[test]
    fn the_disclaimers_are_doc_comments_and_a_bare_literal_is_not_exempt() {
        let bare = r#"macro_rules! u { ($n:ident, $d:literal) => { pub struct $n; }; }
u!(Thermal, "Temperature, in thermals. Not Kelvin, not Celsius.");"#;
        assert!(
            !fiction_hits(bare).is_empty(),
            "a bare literal in a macro's parenthesis group must NOT be exempt — this is \
             the spelling a path exemption would have silenced"
        );

        let doc = r"macro_rules! u { ($(#[$d:meta])* $n:ident) => { $(#[$d])* pub struct $n; }; }
u!(
    /// Temperature, in thermals. Not Kelvin, not Celsius.
    Thermal
);";
        assert!(
            fiction_hits(doc).is_empty(),
            "the shipped spelling must be silent, or borbax-units cannot describe itself"
        );

        // And the rule that was rejected, shown failing: a disclaimer prefix
        // must not launder a real unit in an ordinary label.
        assert!(
            !fiction_hits(r#"fn f() -> &'static str { "Not \u{c5}" }"#).is_empty(),
            "`Not ` must not be a way to write a real unit into a label"
        );
    }

    /// `xtask`'s copy of the blocklist is the generator's blocklist.
    ///
    /// **The highest-value test of the pair.** The copy is what makes the gate
    /// independent of the crate it guards — importing would mean deleting the
    /// blocklist makes the gate *fail to compile* rather than fail. The cost of
    /// a copy is that it can silently narrow while looking maintained, and this
    /// is what stops that: both directions, with a floor on each count so the
    /// test cannot pass over two empty lists.
    #[test]
    fn the_blocklist_the_guard_uses_is_the_blocklist_the_generator_uses() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .map(|p| p.join("crates/borbax-universe/src/naming.rs"));
        let src = path.and_then(|p| std::fs::read_to_string(p).ok());
        assert!(
            src.is_some(),
            "naming.rs could not be read — reported rather than skipped, because a \
             guard that goes quiet when its subject moves is worse than no guard"
        );
        let src = src.unwrap_or_default();

        let extract = |key: &str| -> Vec<String> {
            let Some(from) = src.find(key) else {
                return Vec::new();
            };
            let tail = &src[from..];
            let Some(end) = tail.find("\n];") else {
                return Vec::new();
            };
            tail[..end]
                .split('"')
                .skip(1)
                .step_by(2)
                .map(str::to_owned)
                .collect()
        };

        let symbols = extract("const REAL_ELEMENT_SYMBOLS: &[&str] = &[");
        let words = extract("const REAL_WORDS: &[&str] = &[");

        // Floors first. Without them the whole test is satisfied by two empty
        // lists compared against two empty lists — the vacuous shape this
        // repository has shipped before.
        assert!(
            symbols.len() >= 100,
            "only {} symbols extracted from naming.rs; the extraction broke, not the list",
            symbols.len()
        );
        assert!(
            words.len() >= 80,
            "only {} words extracted from naming.rs; the extraction broke, not the list",
            words.len()
        );

        assert_eq!(
            symbols, REAL_SYMBOLS,
            "xtask's symbol copy has drifted from naming.rs"
        );
        assert_eq!(
            words, REAL_WORDS_COPY,
            "xtask's word copy has drifted from naming.rs"
        );
    }

    /// A literal ending in `r` or `b` is not an element symbol.
    ///
    /// **The regression test for a real false positive.** `trim_matches` removes
    /// *every* leading and trailing character matching its predicate, so the
    /// first version of the content extraction turned `"Nr"` into `N`, `"Ib"`
    /// into `I` and `"bH"` into `H` — reporting nitrogen, iodine and hydrogen on
    /// literals containing no symbol. `a_literal_ending_in_r_is_not_a_format`
    /// pins the identical defect one scan over; this is it arriving again in new
    /// code, which is why the fixture list here is the same shape.
    #[test]
    fn a_literal_ending_in_r_or_b_is_not_an_element_symbol() {
        for src in [
            r#"fn f() -> &'static str { "Nr" }"#,
            r#"fn f() -> &'static str { "Ib" }"#,
            r#"fn f() -> &'static str { "bH" }"#,
            r#"fn f() -> &'static str { "rNr" }"#,
        ] {
            assert!(
                fiction_hits(src).is_empty(),
                "the prefix strip ate a character out of the content: {src}"
            );
        }
        // The positive arm, so this cannot pass by the extraction returning "".
        assert!(
            !fiction_hits(r#"fn f() -> &'static str { "N" }"#).is_empty(),
            "a bare one-letter symbol must still be caught"
        );
        // And a raw string's real content is still reached.
        assert!(
            !fiction_hits("fn f() -> &'static str { r#\"Carbon\"# }").is_empty(),
            "a raw string's content must still be scanned"
        );
    }

    /// A braceless `#[cfg(test)]` does not blind the next block of shipped code.
    ///
    /// `#[cfg(test)]` is legal on an item that opens no brace — `use super::*;`,
    /// `mod tests;`. With the skip flag cleared only by the next brace group, it
    /// survived that item and swallowed the *following* block, which can be
    /// shipped code: the scan then reported nothing about it, silently.
    #[test]
    fn a_cfg_test_item_that_opens_no_block_does_not_skip_the_next_one() {
        let src = r#"
#[cfg(test)]
use super::*;

fn shipped() -> &'static str { "Carbon" }
"#;
        assert!(
            !fiction_hits(src).is_empty(),
            "a braceless #[cfg(test)] swallowed the next block, so shipped code went \
             unscanned"
        );
        // The same shape with `mod tests;` rather than `use`.
        let declared = r#"
#[cfg(test)]
mod tests;

fn shipped() -> &'static str { "Carbon" }
"#;
        assert!(
            !fiction_hits(declared).is_empty(),
            "a braceless `mod tests;` swallowed the next block"
        );

        // The pairing it is supposed to do still works — **including across an
        // intervening attribute**, which the first version of this fix broke:
        // clearing the flag on any non-brace group meant `#[allow(..)]` between
        // the `cfg` and the module made that module scan as shipped code, and
        // its fixtures fired the guard.
        for paired in [
            "#[cfg(test)]\nmod tests { const F: &str = \"Carbon\"; }",
            "#[cfg(test)]\n#[allow(clippy::all)]\nmod tests { const F: &str = \"Carbon\"; }",
        ] {
            assert!(
                fiction_hits(paired).is_empty(),
                "a real test module must still be skipped, or every fixture fires: {paired}"
            );
        }
    }

    /// `eV` and `nM` are caught, not split into single letters.
    ///
    /// `identifier_segments` breaks on camelCase, so `eV` becomes `["e", "v"]`
    /// and `nM` becomes `["n", "m"]` — and those are the spellings anyone
    /// actually writes for electronvolts and nanometres, so the lowercase-only
    /// entries in the vocabulary never matched them.
    #[test]
    fn a_camel_cased_unit_is_still_a_unit() {
        for src in [
            r#"fn f() -> &'static str { "1 eV" }"#,
            r#"fn f() -> &'static str { "12 nM" }"#,
            r#"fn f() -> &'static str { "12 NM" }"#,
        ] {
            assert!(
                !fiction_hits(src).is_empty(),
                "a real unit escaped by its capitalisation: {src}"
            );
        }
    }

    /// A `cfg` that includes `test` in any shape pairs with its module.
    ///
    /// **The string comparison this replaced matched only the literal
    /// `cfg(test)`.** `#[cfg(all(test, feature = "x"))]` is an ordinary
    /// spelling, and against that comparison the module underneath was scanned
    /// as shipped code — so its fixtures, which legitimately name real
    /// elements, fired the guard. A check that fires on correct code gets
    /// deleted.
    #[test]
    fn a_cfg_that_includes_test_in_any_shape_skips_its_module() {
        for attr in [
            "#[cfg(test)]",
            "#[cfg(all(test, feature = \"x\"))]",
            "#[cfg(any(test, doc))]",
        ] {
            let src = format!("{attr}\nmod tests {{ const F: &str = \"Carbon\"; }}");
            assert!(
                fiction_hits(&src).is_empty(),
                "{attr} did not pair with its module, so a test fixture fires the guard"
            );
        }
        // The other direction, and it is what keeps this from becoming a way to
        // hide shipped code: `test` inside a *string* is not a `cfg(test)`.
        let shipped = "#[cfg(feature = \"test-utils\")]\nmod m { const F: &str = \"Carbon\"; }";
        assert!(
            !fiction_hits(shipped).is_empty(),
            "a feature named `test-utils` must not exempt shipped code"
        );
    }

    /// `#[cfg(not(test))]` is shipped-only and must never be skipped.
    ///
    /// **The one direction this guard cannot fail in**, and it was failing:
    /// `mentions_test_ident` recursed into every group, so `not(test)` answered
    /// the opposite of the question asked and the block beneath was treated as
    /// test code. `cfg_test` in this same file already refuses `not` and its
    /// doc records the identical lesson — this function reintroduced it 1900
    /// lines away.
    #[test]
    fn cfg_not_test_is_not_read_as_a_test_module() {
        for attr in [
            "#[cfg(not(test))]",
            "#[cfg(all(not(test), unix))]",
            "#[cfg(any(not(test), doc))]",
        ] {
            let src = format!("{attr}\nmod shipped {{ const F: &str = \"Carbon\"; }}");
            assert!(
                !fiction_hits(&src).is_empty(),
                "{attr} marks SHIPPED code and it was skipped, so a real element name \
                 ships unscanned"
            );
        }
        // The other direction still holds, or the fix has simply disabled the
        // pairing.
        assert!(
            fiction_hits("#[cfg(test)]\nmod tests { const F: &str = \"Carbon\"; }").is_empty(),
            "a real test module must still be skipped"
        );
    }

    /// A call written *below* the test module is still counted.
    ///
    /// **The positional cut this replaced could not see it.** Both call counts
    /// truncated the file at the first `#[cfg(test)]` and scanned only what came
    /// before, so a `Universe::generate` or `.pattern()` after the test module
    /// was invisible — documented as a residual on the grounds that below a test
    /// module is a strange place to hide a call. True, and still a hole that a
    /// re-ordering could open without anyone choosing to.
    ///
    /// It closes the **positional** half only. The alias half —
    /// `use borbax_universe::Universe as U; U::generate(..)` — needs symbol
    /// resolution, is deliberately still open, and says so in the failure
    /// message.
    #[test]
    fn a_call_below_the_test_module_is_still_counted() {
        let below = r"
fn shipped() {}

#[cfg(test)]
mod tests {
    fn t() { let _ = Universe::generate(1); }
}

fn sneaked(u: &U) { let _ = Universe::generate(2); let _ = u.table.pattern(); }
";
        assert_eq!(
            count_outside_tests(below, GENERATE_CALL),
            Some(1),
            "a `Universe::generate` below the test module was not counted, and the one \
             inside it should not be"
        );
        assert_eq!(
            count_outside_tests(below, PATTERN_CALL),
            Some(1),
            "a `.pattern()` below the test module was not counted"
        );

        // The legitimate direction: calls inside the test module are exempt, or
        // every crate's own tests fail the gate.
        let inside = r"
#[cfg(test)]
mod tests {
    fn t() { let _ = Universe::generate(1); let _ = Universe::generate(2); }
}
";
        assert_eq!(
            count_outside_tests(inside, GENERATE_CALL),
            Some(0),
            "a test module's own calls must stay exempt"
        );

        // **A call WITH arguments still counts**, and this is the fixture the
        // needle's first form failed: requiring the closing paren matched only
        // an empty argument list, so `pattern(&table)` scored nothing and the
        // guard would go silent the day the method took an argument. The probe
        // that should have caught it did not — the fixture above calls
        // `pattern()` with no arguments, so both spellings matched it.
        assert_eq!(
            count_outside_tests(
                "fn f(u: &U) { let _ = u.table.pattern(&thing, 2); }",
                PATTERN_CALL
            ),
            Some(1),
            "a `.pattern(..)` call with arguments was not counted"
        );

        // `.pattern` as a *field* is not a call — the opening paren is what
        // makes it one, which is why the delimiter is part of the needle.
        assert_eq!(
            count_outside_tests("fn f(x: T) { let _ = x.pattern; }", &[".", "pattern", "("]),
            Some(0),
            "a field read was counted as a call"
        );

        // Unlexable input is reported, never silently counted as zero.
        assert_eq!(count_outside_tests("fn f( {", GENERATE_CALL), None);
    }

    /// An ordinary English sentence is not an element symbol.
    ///
    /// **Seven sentence-initial English words are element symbols**, and the
    /// two-letter tier matched whole *words* inside a literal, so every one of
    /// these fired on correct code: `No` (nobelium), `In` (indium), `At`
    /// (astatine), `Be` (beryllium), `As` (arsenic), `He` (helium), `Am`
    /// (americium). They are exactly what a UI writes.
    ///
    /// The earlier measurement of "0 false positives on this tree" was true and
    /// useless: none of these strings existed yet. A corpus that happens to be
    /// clean says nothing about the rule.
    #[test]
    fn an_ordinary_english_sentence_is_not_an_element_symbol() {
        for line in [
            "No universe loaded",
            "In this period",
            "At the top of the table",
            "Be careful with that seed",
            "As shown above",
            "He typed a name",
            "Am I looking at the right element",
        ] {
            let src = format!("fn f() -> &'static str {{ \"{line}\" }}");
            assert!(
                fiction_hits(&src).is_empty(),
                "an ordinary sentence fired the symbol tier: {line}"
            );
        }
        // The positive arm, so this cannot pass by the tier being switched off.
        for sym in ["He", "Fe", "No", "In", "C", "N"] {
            let src = format!("fn f() -> &'static str {{ \"{sym}\" }}");
            assert!(
                !fiction_hits(&src).is_empty(),
                "a bare symbol literal must still be caught: {sym}"
            );
        }
        // And a real *name* inside a sentence is still caught, because names do
        // not collide with English the way two-letter symbols do.
        assert!(
            !fiction_hits(r#"fn f() -> &'static str { "made of Carbon" }"#).is_empty(),
            "the name tier must still match by word"
        );
    }
}
