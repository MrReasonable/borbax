//! Does a fusion-built periodic table produce a usable chemistry?
//!
//! The design question is whether elements should be *derived* from a packing
//! rule — element N is N copies of one base atom — instead of drawn from
//! fitted curves as Task 4 currently plans. The argument for it is that every
//! property acquires a cause; the argument against is that a beautiful table
//! producing a dead chemistry is still dead.
//!
//! This measures the second half, using the shape machinery that already
//! exists. It answers three questions and nothing else:
//!
//! 1. Do the derived properties span usable ranges, or does packing collapse
//!    them into a narrow band?
//! 2. Does the *radius series* — the one element property the signature
//!    actually consumes — give molecular shapes that spread as well as the
//!    harness's hand-picked ladder?
//! 3. Does the descriptor image keep its effective dimensionality? `alife`
//!    measured 2.90 of 42 under the current radii, with 57% of the variance in
//!    one mode, so this is the number a principled radius series most plausibly
//!    moves.
//!
//! Run: `cargo run -p borbax-experiments --release --bin ptable`

// Every cast below is a small count — element index, shell size, sample
// count — to `f64` for arithmetic. All far below 2^53, none reaching a
// simulation result. This is a measurement probe, not library code; the
// alternative is a `try_from` on every line of what is a spreadsheet.
#![expect(
    clippy::cast_precision_loss,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    clippy::indexing_slicing,
    clippy::needless_range_loop,
    clippy::too_many_lines,
    reason = "measurement probe: small counts to f64 and fixed-size D-by-D matrix \
              arithmetic indexed by loop bounds. Nothing here reaches a simulation \
              result — it exists to produce the numbers in Task 4's design note, and \
              a try_from plus a get().ok_or() on every line of a spreadsheet would \
              bury the arithmetic being checked"
)]

use borbax_experiments::embed::embed;
use borbax_experiments::geodesic::Geodesic;
use borbax_experiments::molecule::Molecule;
use borbax_experiments::rng::Stream;
use borbax_units::det_math;

/// Sample directions. Matches the harness so the numbers are comparable.
const D: usize = 42;

/// Units in the `n`th packing shell of an icosahedral cluster: `10n² + 2`.
///
/// Cumulative totals are the icosahedral magic numbers 13, 55, 147, 309. The
/// first two shell sizes, 12 and 42, are also the first two rungs of the
/// geodesic ladder `geodesic.rs` already builds — the same geometry that makes
/// a rotation an exact permutation would make a shell close.
const fn shell_size(n: usize) -> usize {
    10 * n * n + 2
}

/// A derived element. Nothing here is drawn; every field is a function of
/// `units` and the universe's generated constants.
#[derive(Debug, Clone, Copy)]
struct Derived {
    units: usize,
    shell: usize,
    /// Units in the incomplete outer shell.
    outer: usize,
    mass: f64,
    radius: f64,
    valence: u8,
    /// Surface-to-bulk ratio of the packed cluster.
    affinity: f64,
    /// Distance from the binding-curve peak, normalised.
    instability: f64,
}

/// Build a table by fusion: element `N` is `N` base units, packed.
fn derive_table(
    n_elements: usize,
    base_mass: f64,
    base_radius: f64,
    peak: f64,
    bulge: f64,
) -> Vec<Derived> {
    let mut out = Vec::with_capacity(n_elements);
    for units in 1..=n_elements {
        // Which packing shell is being filled, and how far into it.
        let (mut shell, mut filled) = (0_usize, 1_usize); // the core unit
        while filled + shell_size(shell + 1) <= units {
            shell += 1;
            filled += shell_size(shell);
        }
        let outer = units - filled;
        let cap = shell_size(shell + 1);

        // Mass: integer multiples of the base unit, less a binding defect that
        // grows with how well-packed the cluster is. Exact on the fixed-point
        // grid by construction, which is what `Mass` was built for.
        let closed = outer == 0;
        let defect = if closed { 0.012 } else { 0.004 };
        let mass = units as f64 * base_mass * (1.0 - defect);

        // Radius: packing gives the cube root — but that alone compresses the
        // series badly (measured: 1.59x range against the hand-picked ladder's
        // 2.55x, costing effective dimensionality). A partially-filled outer
        // shell is not smooth: it bulges. Adding an occupancy term gives a
        // sawtooth that grows through a shell and drops at closure.
        //
        // Worth noting for G3: real atomic radius *decreases* across a period
        // and jumps up at a new one. Packing gives the opposite sign, so this
        // is structurally anti-isomorphic rather than a renamed trend.
        let fill = if cap == 0 {
            0.0
        } else {
            outer as f64 / cap as f64
        };
        let radius = base_radius * det_math::powf(units as f64, 1.0 / 3.0) * (1.0 + bulge * fill);

        // Valence: bonding slots are unpaired sites in the incomplete outer
        // shell. A closed shell has none — the noble family falls out rather
        // than being placed at the end of a period.
        let valence = if closed {
            0
        } else {
            let toward = (outer.min(cap - outer)) as f64;
            (toward.sqrt().round() as u8).min(6)
        };

        // Affinity: the fraction of units exposed on the cluster surface. A
        // shape property of the nucleus, which is the kind of thing this
        // project is made of.
        let surface = if closed { cap } else { outer.max(1) };
        let ratio = surface as f64 / units as f64;
        let affinity = (2.0 * ratio - 1.0).clamp(-1.0, 1.0);

        // Instability: distance from the binding peak, with closed shells
        // stabilised. An iron-peak analogue emerges instead of being drawn.
        let d = (units as f64 - peak) / peak;
        let raw_inst = d * d + if closed { 0.0 } else { 0.15 };
        let instability = if raw_inst > 1.0 { 1.0 } else { raw_inst };

        out.push(Derived {
            units,
            shell,
            outer,
            mass,
            radius,
            valence,
            affinity,
            instability,
        });
    }
    out
}

/// Support-function signature with pluggable radii, so two radius series can
/// be compared on identical molecules and identical embeddings.
fn signature_with(
    g: &Geodesic<D>,
    coords: &[[f64; 3]],
    elements: &[u8],
    radii: &[f64],
) -> Vec<f64> {
    let mut out = vec![f64::NEG_INFINITY; D];
    for (i, u) in g.dirs.iter().enumerate() {
        let mut best = f64::NEG_INFINITY;
        for (a, p) in coords.iter().enumerate() {
            let r = radii.get(usize::from(elements[a])).copied().unwrap_or(0.0);
            let reach = p[0] * u[0] + p[1] * u[1] + p[2] * u[2] + r;
            if reach > best {
                best = reach;
            }
        }
        out[i] = best;
    }
    out
}

fn d_raw(a: &[f64], b: &[f64]) -> f64 {
    let mut acc = 0.0;
    for i in 0..a.len().min(b.len()) {
        let d = a[i] - b[i];
        acc += d * d;
    }
    acc.sqrt()
}

/// Participation ratio of the covariance spectrum — the number of directions
/// the variance is actually spread over. `alife`'s metric.
fn effective_dim(sigs: &[Vec<f64>]) -> (f64, f64) {
    let n = sigs.len() as f64;
    let mut mean = vec![0.0; D];
    for s in sigs {
        for i in 0..D {
            mean[i] += s[i] / n;
        }
    }
    // Covariance eigenvalues via the Gram matrix's trace powers is overkill;
    // the participation ratio only needs tr(C) and tr(C²), both computable
    // from the centred data without an eigendecomposition.
    let mut c = vec![vec![0.0; D]; D];
    for s in sigs {
        for i in 0..D {
            for j in 0..D {
                c[i][j] += (s[i] - mean[i]) * (s[j] - mean[j]) / n;
            }
        }
    }
    let mut tr = 0.0;
    let mut tr2 = 0.0;
    for i in 0..D {
        tr += c[i][i];
        for j in 0..D {
            tr2 += c[i][j] * c[j][i];
        }
    }
    // Largest eigenvalue by power iteration, for the "top PC" share.
    let mut v = vec![1.0 / (D as f64).sqrt(); D];
    for _ in 0..200 {
        let mut w = vec![0.0; D];
        for i in 0..D {
            for j in 0..D {
                w[i] += c[i][j] * v[j];
            }
        }
        let norm = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm <= 0.0 {
            break;
        }
        for i in 0..D {
            v[i] = w[i] / norm;
        }
    }
    let mut lam = 0.0;
    for i in 0..D {
        for j in 0..D {
            lam += v[i] * c[i][j] * v[j];
        }
    }
    (tr * tr / tr2, lam / tr)
}

fn main() {
    println!("== a fusion-built table: element N is N base units, packed ==\n");
    println!(
        "shell sizes 10n^2+2: {:?}",
        (1..=4).map(shell_size).collect::<Vec<_>>()
    );
    let cum: Vec<usize> = (1..=4)
        .scan(1, |a, n| {
            *a += shell_size(n);
            Some(*a)
        })
        .collect();
    println!("cumulative (magic numbers): {cum:?}\n");

    let table = derive_table(120, 1.7, 0.42, 62.0, 0.55);

    println!(
        "{:>4} {:>5} {:>5} {:>8} {:>7} {:>4} {:>8} {:>8}",
        "N", "shell", "outer", "mass", "radius", "val", "affin", "instab"
    );
    for e in table.iter().filter(|e| {
        // The interesting rows: shell boundaries and a few in between.
        e.outer == 0 || e.outer == 1 || e.units % 17 == 0
    }) {
        println!(
            "{:>4} {:>5} {:>5} {:>8.2} {:>7.3} {:>4} {:>8.3} {:>8.3}",
            e.units, e.shell, e.outer, e.mass, e.radius, e.valence, e.affinity, e.instability
        );
    }

    println!("\n-- property spans over 120 elements --");
    let span = |f: fn(&Derived) -> f64| {
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for e in &table {
            let v = f(e);
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        (lo, hi)
    };
    let (ml, mh) = span(|e| e.mass);
    let (rl, rh) = span(|e| e.radius);
    let (al, ah) = span(|e| e.affinity);
    let (il, ih) = span(|e| e.instability);
    println!("mass        {ml:8.2} .. {mh:8.2}");
    println!("radius      {rl:8.3} .. {rh:8.3}");
    println!("affinity    {al:8.3} .. {ah:8.3}");
    println!("instability {il:8.3} .. {ih:8.3}");
    let closed: Vec<usize> = table
        .iter()
        .filter(|e| e.outer == 0)
        .map(|e| e.units)
        .collect();
    println!("closed shells (the emergent noble family): {closed:?}");
    let vals: Vec<u8> = (0..=6)
        .map(|v| table.iter().filter(|e| e.valence == v).count() as u8)
        .collect();
    println!("valence histogram 0..=6: {vals:?}");

    // -- The question that decides it: does shape still spread? ------------
    println!("\n== shape diversity: harness ladder vs packing-derived radii ==\n");
    let Ok(g) = Geodesic::<D>::build() else {
        println!("geodesic failed to build");
        return;
    };
    let mut rng = Stream::new(0x0B0B_BA05);
    let mols: Vec<Molecule> = (0..400)
        .map(|_| Molecule::random_tree(&mut rng, 14))
        .collect();
    let coords: Vec<Vec<[f64; 3]>> = mols.iter().map(embed).collect();

    // The harness's hand-picked ladder, and the first four derived radii
    // rescaled to the same mean so only the *shape* of the series differs.
    let ladder = [0.55, 0.80, 1.05, 1.40];
    let mean_l = ladder.iter().sum::<f64>() / 4.0;
    let rescale = |v: Vec<f64>| {
        let m = v.iter().sum::<f64>() / v.len() as f64;
        v.iter().map(|r| r * mean_l / m).collect::<Vec<_>>()
    };
    // Pure cube root, and cube root with the shell-occupancy bulge. The first
    // four elements sit in shells 0 and 1, so the bulge is visible there.
    let cbrt_only = rescale(
        (1..=4)
            .map(|n| 0.42 * det_math::powf(n as f64, 1.0 / 3.0))
            .collect(),
    );
    let with_bulge = rescale(
        derive_table(4, 1.7, 0.42, 62.0, 0.55)
            .iter()
            .map(|e| e.radius)
            .collect(),
    );
    // The fair comparison. A universe draws its elements from across the whole
    // table, not from the first four. Taking N = 3, 20, 60, 110 spans three
    // shells and two closures.
    let picks = [3_usize, 20, 60, 110];
    let spread = rescale(
        picks
            .iter()
            .map(|&n| table.get(n - 1).map_or(0.0, |e| e.radius))
            .collect(),
    );
    println!("ladder    {ladder:.3?}");
    println!("cbrt      {cbrt_only:.3?}  (rescaled to the same mean)");
    println!("+bulge    {with_bulge:.3?}  (rescaled to the same mean)");
    println!("spread    {spread:.3?}  (N = 3, 20, 60, 110 — across the table)");

    for (label, radii) in [
        ("ladder", &ladder[..]),
        ("cbrt", &cbrt_only[..]),
        ("+bulge", &with_bulge[..]),
        ("spread", &spread[..]),
    ] {
        let sigs: Vec<Vec<f64>> = mols
            .iter()
            .zip(&coords)
            .map(|(m, c)| signature_with(&g, c, &m.elements, radii))
            .collect();
        let (eff, top) = effective_dim(&sigs);
        let mut ds = Vec::new();
        for i in 0..sigs.len() {
            for j in (i + 1)..sigs.len() {
                ds.push(d_raw(&sigs[i], &sigs[j]));
            }
        }
        let n = ds.len() as f64;
        let mean = ds.iter().sum::<f64>() / n;
        let sd = (ds.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / n).sqrt();
        let collapsed = ds.iter().filter(|&&d| d < 0.10 * mean).count();
        println!(
            "{label:7} eff.dim {eff:5.2}/42   top PC {:5.1}%   CV {:.4}   collapsed pairs {collapsed}",
            100.0 * top,
            sd / mean
        );
    }
}
