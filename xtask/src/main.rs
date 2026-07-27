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

    if failures.is_empty() {
        println!("repository invariants: all checks passed");
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
