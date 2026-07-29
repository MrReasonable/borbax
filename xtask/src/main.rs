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
const FORBIDDEN_FORMAT_TOKENS: &[&str] = &["SMILES", "InChI", "FASTA", "PDB format"];

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
    check_toolchain_pins_agree(root, &mut failures)?;
    check_no_platform_transcendentals(root, &mut failures)?;
    check_no_stream_deriving_method(root, &mut failures)?;

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

/// §13.1 — no method in `borbax-rng` may derive one `Stream` from another.
///
/// **This is the only enforcement of the crate's central duplicate-stream
/// guarantee, and a doctest provably cannot do the job.** `Stream::sub` is a
/// constructor precisely so that `base.sub(a).sub(b)` — which silently
/// returned a duplicate of an unrelated sibling under the old `fork` — has no
/// spelling. But that held only because nobody had written the method:
/// reinstating it verbatim under any other name passes the whole suite and
/// `clippy -D warnings`, measured.
///
/// Two things rule out the obvious alternative. A `compile_fail` doctest naming
/// `sub` cannot catch a regression that arrives as `child` or `derive` — and
/// `sub` cannot be both an associated function and a method (`E0592`), so it
/// *must* arrive under a different name. And `compile_fail,E0599` does not
/// enforce the error code on the pinned toolchain, so such a block passes when
/// it fails for an unrelated reason. `borbax-units` records that measurement
/// already.
///
/// So the check is textual and structural: no `pub fn` in that crate takes
/// `&self` and returns `Self`. Same shape and same reason as the rayon guard —
/// a one-word change that compiles, typechecks, passes clippy, and is
/// invisible in review.
///
/// `Clone` is exempt and must be: it duplicates a stream rather than deriving
/// a different one, `!Copy` exists to make that duplication visible at the
/// call site, and it is a derive rather than a `pub fn`.
fn check_no_stream_deriving_method(root: &Path, failures: &mut Vec<String>) -> Result<(), String> {
    let path = root.join("crates/borbax-rng/src/lib.rs");
    if !path.exists() {
        return Ok(());
    }
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    for (i, line) in text.lines().enumerate() {
        let code = line.split("//").next().unwrap_or("").trim();
        if !code.starts_with("pub fn") && !code.starts_with("pub const fn") {
            continue;
        }
        if code.contains("&self") && (code.contains("-> Self") || code.contains("-> Stream")) {
            failures.push(format!(
                "§13.1: `{code}` at crates/borbax-rng/src/lib.rs:{} derives a Stream from a \
                 Stream — that is how a second sub-level silently returns a duplicate of an \
                 unrelated sibling. Sub-streams are constructed, not derived.",
                i + 1
            ));
        }
    }
    Ok(())
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
        // Not yet written — Task 4 creates it. Not a failure before then.
        return Ok(());
    }
    let src = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
    if !src.contains("REAL_ELEMENT_SYMBOLS") {
        failures.push("G2: naming.rs has no REAL_ELEMENT_SYMBOLS blocklist".into());
    }
    for token in FORBIDDEN_FORMAT_TOKENS {
        if src.contains(token) {
            failures.push(format!("G5: forbidden format token {token:?} in naming.rs"));
        }
    }
    Ok(())
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
                if path.file_name().and_then(|n| n.to_str()) != Some("target") {
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

/// The §13.1 scanner is the only piece of logic here that can fail *quietly* —
/// every other check either finds its target or does not exist. A miscounted
/// brace makes it skip the rest of a file while still reporting success, so
/// each way it could miscount gets a test.
#[cfg(test)]
mod tests {
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
