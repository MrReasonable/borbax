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

/// Independent seeds per regime.
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
            let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
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
