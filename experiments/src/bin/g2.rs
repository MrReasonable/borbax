//! Print the G2 locality measurement in both regimes.
//!
//! ```text
//! cargo run -p borbax-experiments --release --bin g2
//! ```
//!
//! Runs several independent seeds and reports the spread as well as the
//! numbers. A single-seed figure is exactly the kind of result that cannot be
//! reproduced or argued with later — the spread is what says whether a gap
//! between two descriptors is real or is one draw's noise.
//!
//! Exits non-zero if any regime's controls fail, so this cannot run in CI and
//! quietly report numbers that mean nothing.

use borbax_experiments::g2::{Config, Report, run};
use borbax_experiments::signature::Descriptor;

/// Independent seeds, one per regime.
///
/// **Arbitrary values. The only requirement is that they are fixed and
/// distinct.**
///
/// An earlier version of this doc said they must also be "not adjacent, so two
/// regimes cannot accidentally share a stream prefix". That is false, and this
/// crate proves it: `adjacent_seeds_decorrelate_immediately` shows seeds 0..63
/// give first draws differing in at least 18 of 64 bits across all 2016 pairs.
/// The finaliser is what makes adjacency harmless — asserting a requirement the
/// mixer already discharges invites someone to "fix" a perfectly good seed set.
///
/// Fixed, because a measurement reproduces only if its stream does. Written as
/// recognisable bit patterns rather than 1, 2, 3 so a value appearing somewhere
/// unexpected is obviously *this array* and not a count.
const SEEDS: [u64; 5] = [
    0x0B0B_BA05,
    0x1234_5678,
    0xFEED_FACE,
    0x0000_0001,
    0xA5A5_A5A5,
];

fn main() -> std::process::ExitCode {
    let mut healthy = true;

    for size_matched in [false, true] {
        let reports: Vec<Report> = SEEDS
            .iter()
            .map(|&seed| {
                run(&Config {
                    seed,
                    size_matched,
                    ..Config::default()
                })
            })
            .collect();

        let Some(first) = reports.first() else {
            continue;
        };
        println!("{}", first.render());

        println!("Concordance across {} seeds:", SEEDS.len());
        println!(
            "{:<24} {:>10} {:>10} {:>10}",
            "descriptor", "min", "max", "spread"
        );
        println!("{}", "-".repeat(56));
        for d in Descriptor::ALL {
            let values: Vec<f64> = reports
                .iter()
                .filter_map(|r| r.get(d))
                .map(|s| s.concordance)
                .collect();
            // Not `fold(_, f64::min)`: §13.1 bans it because std documents that
            // when the inputs compare equal — which `+0.0` and `-0.0` do —
            // either may be returned non-deterministically, and it was measured
            // returning different signs on aarch64 and x86-64. This is report
            // formatting rather than a result, but the ban is worth more than
            // the exemption would be.
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for v in values.iter().copied() {
                if v < lo {
                    lo = v;
                }
                if v > hi {
                    hi = v;
                }
            }
            println!("{:<24} {lo:>10.4} {hi:>10.4} {:>10.4}", d.name(), hi - lo);
        }
        println!();

        if reports.iter().any(|r| r.check_controls().is_err()) {
            healthy = false;
        }
    }

    if healthy {
        std::process::ExitCode::SUCCESS
    } else {
        eprintln!("controls failed on at least one seed — the numbers above are not a G2 result");
        std::process::ExitCode::FAILURE
    }
}
