//! Fiction-guarantee enforcement (spec §5). Run in CI and from the
//! pre-commit hook.
//!
//! These checks make "no real chemistry" an architectural property rather
//! than a promise. They are cheap and they run on every commit.
//!
//! Git hooks are not this binary's job — `prek` installs and runs them from
//! `.pre-commit-config.yaml`, which calls back into `cargo run -p xtask`.

use std::path::{Path, PathBuf};

/// Extensions that would indicate bulk data smuggled into a chemistry crate.
/// The entire chemistry is generated from a seed; there is nothing legitimate
/// for a data file to be doing in these crates.
const DATA_EXTENSIONS: &[&str] = &["csv", "tsv", "json", "yaml", "yml", "parquet", "bin", "dat"];

/// Crates that must contain no data files at all (G1).
const CHEMISTRY_CRATES: &[&str] = &[
    "crates/borbax-universe",
    "crates/borbax-molecule",
    "crates/borbax-reaction",
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
        Some(other) => Err(format!(
            "unknown subcommand {other:?} — expected `check-guarantees` or no argument"
        )),
    }
}

fn check_guarantees(root: &Path) -> Result<(), String> {
    let mut failures = Vec::new();
    check_no_data_files(root, &mut failures)?;
    check_blocklist_present(root, &mut failures)?;

    if failures.is_empty() {
        println!("fiction guarantees: all checks passed");
        Ok(())
    } else {
        for f in &failures {
            eprintln!("FAIL: {f}");
        }
        Err(format!(
            "{} fiction-guarantee check(s) failed",
            failures.len()
        ))
    }
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
