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
//! 3. Does the descriptor image keep its effective dimensionality?
//!
//! **Read the answer to 2 and 3 as "costs nothing", not "improves".** An
//! earlier version of this comment cited `alife`'s 2.90-of-42 next to this
//! probe's 7.59, which reads as a large improvement and is not one: 2.90 was
//! measured on molecules of *free* size (8-20 atoms), where gross size is one
//! dominant mode taking 57% of the variance, and this probe fixes every
//! molecule at 14 atoms, which removes that mode. Same metric, different
//! populations. State the population beside any figure from here.
//!
//! Two further limits on question 3, both measured after the fact:
//!
//! - **The participation ratio on raw signatures is gameable.** Spinning each
//!   molecule by a random one of the 60 rotations — a transformation `D_group`
//!   is exactly invariant to, changing no shape at all — moves the score from
//!   7.59 to 13.24. A diversity statistic for this project must be invariant
//!   under the group, because the adopted descriptor is.
//! - **It is not what §7.2 asks for.** §7.2 wants signatures spread rather than
//!   collapsed into a few clusters — a *cluster count*. The two disagree in the
//!   case that matters: uniform radii collapse 574 distinct molecules to 12
//!   distinguishable shapes while the participation ratio falls only 38%, to a
//!   value no plausible floor would reject.
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
use borbax_experiments::geodesic::{Geodesic, apply_mat, rotation_matrices};
use borbax_experiments::molecule::Molecule;
use borbax_experiments::rng::Stream;
use borbax_experiments::signature::d_group;
use borbax_units::det_math;

/// Sample directions. Matches the harness so the numbers are comparable.
const D: usize = 42;

/// Units in the `n`th packing shell of an icosahedral cluster: `10n² + 2`.
///
/// **Defined for `n >= 1`.** `shell_size(0)` returns 2, which is meaningless —
/// shell 0 is the single central unit. Nothing calls it at 0: the fill loop
/// always evaluates `shell_size(shell + 1)`, and the burial calculation
/// special-cases `shell == 0` to 1.
///
/// Cumulative totals are 1, 13, 55, 147, 309.
///
/// **On the relationship to `geodesic.rs`, which an earlier version of this
/// comment overstated.** It claimed "the same geometry that makes a rotation an
/// exact permutation would make a shell close". What the two objects genuinely
/// share is the icosahedron's face tiling and its rotation group. What they do
/// **not** share is the point set, and the ladders are different sequences:
/// geodesic is `10·4^L + 2` = 12, 42, 162, 642, packing is `10n² + 2` = 12, 42,
/// 92, 162, 252. They coincide only where `n = 2^L`, so **92 has no geodesic
/// counterpart at all**, and 162 is packing shell 4 against geodesic level 2 —
/// measured as *different* point sets, mismatch 2.29e-2, which is 8% of the
/// nearest-neighbour spacing and nowhere near `exact_index`'s 1e-12.
/// Recursive midpoint-then-normalise is not radial projection of a planar
/// triangular grid. A shell's neighbour spacing is uniform; a geodesic's is not.
///
/// So do not reuse `Geodesic::dirs` as cluster site positions. It is right at
/// 12 and 42 and silently wrong at 162.
const fn shell_size(n: usize) -> usize {
    10 * n * n + 2
}

/// Binding defect for a closed cluster, as a fraction of its unpacked mass.
///
/// Dyadic on purpose — see the derivation at the `mass` computation. `1/128`.
const DEFECT_CLOSED: f64 = 0.007_812_5;

/// Binding defect for a cluster with a partly-filled outer shell. `1/256`.
const DEFECT_OPEN: f64 = 0.003_906_25;

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
        // grows with how well-packed the cluster is.
        //
        // **The defect factors are dyadic, and that is load-bearing rather than
        // tidy.** `Mass` is fixed-point at SCALE = 1024, so a mass series is
        // exact on its grid only if every factor has a denominator dividing
        // 1024. The values here stood at 0.012 and 0.004 under a comment
        // claiming exactness — measured, **0 of 120** elements landed on the
        // grid, because 0.988 is 247/250 and 250 does not divide 1024. With
        // 1/128 and 1/256 against a 7/4 base the products are 889/512 and
        // 1785/1024, so every element is an exact multiple of 1/1024 and the
        // claim is true rather than believed.
        //
        // The magnitudes also had to move. At 0.012/0.004 the series stops
        // increasing at N = 124 — `mass(146) > mass(147)` — so the third
        // closed shell weighs less than the element below it, against §7.1's
        // "increases down the table". Monotonicity holds while
        // `N < (1 - d_closed) / (d_closed - d_open)`; these values carry it to
        // N = 307, breaking only at the fifth closure (309), far past any table
        // V0 will build. `mass_increases_with_unit_count` asserts it.
        let closed = outer == 0;
        let defect = if closed { DEFECT_CLOSED } else { DEFECT_OPEN };
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
        // `cap` is `shell_size(shell + 1)` = 10(shell+1)^2 + 2, so it is never
        // below 12. A `cap == 0` guard stood here and could not fire.
        let fill = outer as f64 / cap as f64;
        let radius = base_radius * det_math::cbrt(units as f64) * (1.0 + bulge * fill);

        // Valence: bonding slots are unpaired sites in the incomplete outer
        // shell. A closed shell has none — the closed-shell family falls out rather
        // than being placed at the end of a period.
        //
        // **There is deliberately no `if closed { 0 }` branch, and that is the
        // strongest emergence result in this file.** One stood here. Deleting
        // it changes the valence of *zero* of 120 elements, because `outer` is
        // 0 at a closure and `min(0, cap)` is 0 already. By the deletion test
        // the closed-shell family is a thing the arithmetic *detects*, not a
        // case the code *makes* — §3's distinction exactly. A redundant special
        // case on a chemically-loaded predicate is how a detector quietly
        // becomes a declaration two refactors later.
        //
        // The `sqrt` and the `min(6)` are both known-bad and are replaced
        // wholesale in Task 4: they are a fitted curve, the clamp caps this at
        // 6 where 8+ is wanted, and geometry measured that a per-site count
        // *cannot* exceed 6 (coordination 12 minus a vertex site's 6). Reaching
        // 8 needs a count over sites plus a canonical filling order.
        let toward = (outer.min(cap - outer)) as f64;
        let valence = (toward.sqrt().round() as u8).min(6);

        // Affinity: the fraction of units exposed on the cluster surface. A
        // shape property of the nucleus, which is the kind of thing this
        // project is made of.
        //
        // A unit is buried when its shell is covered by the one outside it.
        // Shells strictly inside the outermost occupied one are covered
        // completely; the outermost *complete* shell is covered in proportion
        // to how far the next shell has filled. That proportion is the whole
        // model — there is no coefficient to tune.
        //
        // It is continuous at every closure by construction: at `outer == 0`
        // it reduces to `surface == shell_units(shell)`, and at `outer == cap`
        // it reduces to the same expression for the next shell. **An earlier
        // version was not.** It read `if closed { cap }` — the capacity of the
        // shell that is still *empty* — giving N=13 a surface of 42 on a
        // 13-unit cluster, a ratio of 3.23 clamped to exactly +1.000. Every
        // closed shell scored maximum affinity and its N+1 neighbour scored
        // near minimum, a swing of 1.86 where the real exposed fraction moves
        // by less than 0.01. That made the closed-shell elements simultaneously
        // the least reactive by valence and the *stickiest* by affinity — and
        // §8.3 runs the one mechanism on affinity, so nothing would have
        // accumulated.
        let shell_units = if shell == 0 { 1 } else { shell_size(shell) };
        let inner = filled - shell_units;
        let covered = shell_units as f64 * outer as f64 / cap as f64;
        let surface = units as f64 - (inner as f64 + covered);

        // In (0, 1] by construction — `inner + covered` is a subset of `units`
        // — so `2r - 1` lands in (-1, 1] with no clamp. A clamp here would be
        // a guard that can never fire, which is the shape of defect this file
        // is partly about. `affinity_stays_in_range` asserts it instead.
        let ratio = surface / units as f64;
        let affinity = 2.0 * ratio - 1.0;

        // Instability: distance from the binding peak, with closed shells
        // stabilised.
        //
        // **`peak` is a fudge factor and this is not an emergent minimum.** With
        // the open-shell bonus at 0 the minimum sits at N = 62 — the value of
        // `peak` itself — and 0.013 is just `((55-62)/62)^2`. Sweeping `peak`
        // over 5..120 the minimum lands on a closed shell in 13 of 24 draws. It
        // is also the second encoding of "closed shells bind better", the first
        // being the mass defect above, in incommensurable units so neither can
        // check the other. Task 4 derives one binding energy per unit from the
        // packing counts and makes both a function of it, so `peak` becomes an
        // output rather than an argument.
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

/// Molecules drawn per seed in the diversity measurement.
const MOLS_PER_SEED: usize = 200;

/// Random 4-element draws used to size the *sampling* variation, which is the
/// term an earlier version of this probe left out and then reported inside.
const RANDOM_DRAWS: usize = 200;

/// Effective dimensionality that is invariant under the 60 rotations.
///
/// **Why this exists alongside [`effective_dim`].** That one takes a
/// participation ratio of the raw signature covariance, and a raw signature is
/// not rotation-invariant even though `D_group` is. Measured: spinning every
/// molecule by a group element — changing no shape at all — moves the raw score
/// from 7.59 to 13.24. A diversity statistic for this project has to be blind
/// to exactly what the adopted descriptor is blind to, or a universe can score
/// well on orientation noise.
///
/// Classical MDS on the `D_group` distance matrix, read through the same
/// participation ratio. `tr(B)^2 / tr(B^2)` needs no eigendecomposition, since
/// `tr(B) = sum(lambda)` and `tr(B^2) = sum(lambda^2) = sum of squared entries`
/// for symmetric `B`. `D_group` is not guaranteed Euclidean, so `B` may carry
/// small negative eigenvalues; they enter `tr(B^2)` positively and depress the
/// ratio slightly, which is the conservative direction.
///
/// Returns `(effective_dim, collapsed_pairs)`. The second is a control that
/// genuinely fires — uniform radii give 23 collapsed pairs on this population
/// and 3-atom trees give 4145 — so it is reported rather than dropped when the
/// raw-distance path went away.
///
/// Accumulation is in index order throughout and must stay that way (§13.4).
fn group_invariant_dim(g: &Geodesic<D>, sigs: &[Vec<f64>]) -> (f64, usize) {
    let n = sigs.len();
    if n < 2 {
        return (0.0, 0);
    }
    let mut d2 = vec![vec![0.0_f64; n]; n];
    let mut flat = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let d = d_group(g, &sigs[i], &sigs[j]);
            d2[i][j] = d * d;
            d2[j][i] = d2[i][j];
            flat.push(d);
        }
    }
    let mean_d = flat.iter().sum::<f64>() / flat.len() as f64;
    let collapsed = flat.iter().filter(|&&d| d < 0.10 * mean_d).count();
    let nf = n as f64;
    let row: Vec<f64> = (0..n).map(|i| d2[i].iter().sum::<f64>() / nf).collect();
    let grand = row.iter().sum::<f64>() / nf;
    let (mut trace, mut trace_sq) = (0.0_f64, 0.0_f64);
    for i in 0..n {
        for j in 0..n {
            let b = -0.5 * (d2[i][j] - row[i] - row[j] + grand);
            if i == j {
                trace += b;
            }
            trace_sq += b * b;
        }
    }
    if trace_sq <= 0.0 {
        return (0.0, collapsed);
    }
    (trace * trace / trace_sq, collapsed)
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
    println!("cumulative shell totals: {cum:?}\n");

    let table = derive_table(120, 1.75, 0.42, 62.0, 0.55);

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
    println!("closed shells: {closed:?}");
    let vals: Vec<u8> = (0..=6)
        .map(|v| table.iter().filter(|e| e.valence == v).count() as u8)
        .collect();
    println!("valence histogram 0..=6: {vals:?}");

    // -- Does shape still spread? -----------------------------------------
    //
    // **This section reports distributions, not single numbers, and that is a
    // repair.** An earlier version drew one molecule seed, hand-picked four
    // elements, and concluded the derived table was "slightly better on every
    // measure". Over 8 seeds the difference is smaller than the between-seed
    // spread and points the other way; the hand-picked draw sits in the top
    // few percent of what that choice could have produced. A statistic quoted
    // without its variation is an opinion with a decimal point.
    println!("\n== shape diversity: harness ladder vs packing-derived radii ==\n");
    let Ok(g) = Geodesic::<D>::build() else {
        println!("geodesic failed to build");
        return;
    };
    let Ok(mats) = rotation_matrices() else {
        println!("rotation table failed to build");
        return;
    };

    let ladder = [0.55, 0.80, 1.05, 1.40];
    let mean_l = ladder.iter().sum::<f64>() / 4.0;
    let rescale = |v: Vec<f64>| {
        let m = v.iter().sum::<f64>() / v.len() as f64;
        v.iter().map(|r| r * mean_l / m).collect::<Vec<_>>()
    };
    let radii_for = |picks: &[usize]| {
        rescale(
            picks
                .iter()
                .map(|&n| table.get(n - 1).map_or(0.0, |e| e.radius))
                .collect(),
        )
    };
    // The draw an earlier version reported. Kept so its percentile can be
    // stated rather than its value re-quoted.
    let picks = [3_usize, 20, 60, 110];
    let spread = radii_for(&picks);
    println!("ladder    {ladder:.3?}");
    println!("spread    {spread:.3?}  (N = 3, 20, 60, 110 — the earlier draft's pick)\n");

    // --- Over 8 molecule seeds, paired: same molecules feed both arms.
    let mut l_eff = Vec::new();
    let mut s_eff = Vec::new();
    let mut l_inv = Vec::new();
    let mut s_inv = Vec::new();
    let mut paired = Vec::new();
    let (mut collapsed_l, mut collapsed_s) = (0_usize, 0_usize);
    for seed_ix in 0..8_u64 {
        let mut rng = Stream::new(0x0B0B_BA05 ^ seed_ix);
        let mols: Vec<Molecule> = (0..MOLS_PER_SEED)
            .map(|_| Molecule::random_tree(&mut rng, 14))
            .collect();
        let coords: Vec<Vec<[f64; 3]>> = mols.iter().map(embed).collect();
        let sig = |radii: &[f64]| -> Vec<Vec<f64>> {
            mols.iter()
                .zip(&coords)
                .map(|(m, c)| signature_with(&g, c, &m.elements, radii))
                .collect()
        };
        let (sl, ss) = (sig(&ladder), sig(&spread));
        let (el, _) = effective_dim(&sl);
        let (es, _) = effective_dim(&ss);
        l_eff.push(el);
        s_eff.push(es);
        paired.push(es - el);
        let (li, lc) = group_invariant_dim(&g, &sl);
        let (si, sc) = group_invariant_dim(&g, &ss);
        l_inv.push(li);
        s_inv.push(si);
        collapsed_l += lc;
        collapsed_s += sc;
    }
    let stat = |v: &[f64]| -> (f64, f64) {
        let n = v.len() as f64;
        let m = v.iter().sum::<f64>() / n;
        (
            m,
            (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n).sqrt(),
        )
    };
    let (lm, ls) = stat(&l_eff);
    let (sm, ss_) = stat(&s_eff);
    let (lim, lis) = stat(&l_inv);
    let (sim, sis) = stat(&s_inv);
    let (pm, ps) = stat(&paired);
    println!("over 8 molecule seeds, {MOLS_PER_SEED} molecules of 14 atoms each:\n");
    println!("series   eff.dim (raw covariance)   eff.dim (D_group, invariant)   collapsed");
    println!(
        "ladder   {lm:6.3} +/- {ls:.3}          {lim:6.3} +/- {lis:.3}          {collapsed_l}"
    );
    println!(
        "spread   {sm:6.3} +/- {ss_:.3}          {sim:6.3} +/- {sis:.3}          {collapsed_s}"
    );
    println!(
        "\npaired difference (spread - ladder): {pm:+.3} +/- {ps:.3}  -> {}",
        if pm.abs() < ps {
            "inside its own noise"
        } else {
            "outside the noise"
        }
    );

    // --- How much of that is the *element draw* rather than the scheme?
    let mut rng = Stream::new(0x5EED_D2A4);
    let mut mols_rng = Stream::new(0x0B0B_BA05);
    let mols: Vec<Molecule> = (0..MOLS_PER_SEED)
        .map(|_| Molecule::random_tree(&mut mols_rng, 14))
        .collect();
    let coords: Vec<Vec<[f64; 3]>> = mols.iter().map(embed).collect();
    let mut draws = Vec::new();
    for _ in 0..RANDOM_DRAWS {
        let mut pick = [0_usize; 4];
        for slot in &mut pick {
            *slot = 1 + (rng.next_u64() % table.len() as u64) as usize;
        }
        let radii = radii_for(&pick);
        let sigs: Vec<Vec<f64>> = mols
            .iter()
            .zip(&coords)
            .map(|(m, c)| signature_with(&g, c, &m.elements, &radii))
            .collect();
        draws.push(effective_dim(&sigs).0);
    }
    draws.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    let pct = |q: f64| draws[((draws.len() as f64 - 1.0) * q) as usize];
    let hand = {
        let sigs: Vec<Vec<f64>> = mols
            .iter()
            .zip(&coords)
            .map(|(m, c)| signature_with(&g, c, &m.elements, &spread))
            .collect();
        effective_dim(&sigs).0
    };
    let below = draws.iter().filter(|&&d| d < hand).count();
    println!("\n{RANDOM_DRAWS} random 4-element draws from the same table:");
    println!(
        "  eff.dim  p10 {:.2}   median {:.2}   p90 {:.2}",
        pct(0.10),
        pct(0.50),
        pct(0.90)
    );
    println!(
        "  the earlier draft's hand-picked draw scores {hand:.2} — {}th percentile",
        100 * below / draws.len()
    );

    // --- Is the statistic itself measuring shape, or orientation?
    //
    // Spin every molecule by a group element. `D_group` is exactly invariant to
    // that, so no shape changes. Anything the statistic reports as a change is
    // the statistic's own artefact.
    let mut spin = Stream::new(0x5717_0001);
    let spun: Vec<Vec<[f64; 3]>> = coords
        .iter()
        .map(|c| {
            let r = (spin.next_u64() % mats.len() as u64) as usize;
            mats.get(r).map_or_else(
                || c.clone(),
                |mat| c.iter().map(|&p| apply_mat(mat, p)).collect(),
            )
        })
        .collect();
    let sig_of = |cs: &[Vec<[f64; 3]>]| -> Vec<Vec<f64>> {
        mols.iter()
            .zip(cs)
            .map(|(m, c)| signature_with(&g, c, &m.elements, &ladder))
            .collect()
    };
    let (still, spun_sigs) = (sig_of(&coords), sig_of(&spun));
    println!("\ngameability: every molecule spun by a random group element (no shape changes)");
    println!(
        "  raw covariance  {:6.2} -> {:6.2}   <- moves, so it is not measuring shape alone",
        effective_dim(&still).0,
        effective_dim(&spun_sigs).0
    );
    println!(
        "  D_group MDS     {:6.2} -> {:6.2}   <- invariant, as the adopted descriptor is",
        group_invariant_dim(&g, &still).0,
        group_invariant_dim(&g, &spun_sigs).0
    );
}

#[cfg(test)]
#[expect(
    clippy::float_cmp,
    reason = "the float comparisons are against exact multiples of 1/1024, which is \
              the property being asserted — an epsilon would defeat the test"
)]
mod tests {
    use super::*;

    /// The table Task 4's note reports on.
    fn table() -> Vec<Derived> {
        derive_table(120, 1.75, 0.42, 62.0, 0.55)
    }

    /// `affinity` is `2 * surface_fraction - 1`, and a surface fraction is a
    /// fraction — so this cannot exceed `(-1, 1]` unless burial is computed
    /// wrongly.
    ///
    /// It was. The previous formulation took a *closed* cluster's surface to be
    /// the capacity of the shell outside it, so N=13 scored 42 exposed units on
    /// 13 units total: a ratio of 3.23, saved from being visibly absurd only by
    /// a clamp.
    #[test]
    fn affinity_stays_in_range() {
        for e in table() {
            assert!(
                e.affinity > -1.0 && e.affinity <= 1.0,
                "N={}: affinity {} outside (-1, 1] — burial is miscounted",
                e.units,
                e.affinity
            );
        }
    }

    /// Adding one unit to a cluster changes its exposed fraction by about one
    /// unit's worth. The bound `2/N` is loose on purpose: the point is to catch
    /// a *discontinuity*, and the old code jumped 1.86 at every closure while
    /// the real exposed fraction moved by under 0.01.
    #[test]
    fn affinity_is_continuous_across_shell_closures() {
        let t = table();
        for pair in t.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let bound = 2.0 / a.units as f64;
            assert!(
                (b.affinity - a.affinity).abs() <= bound,
                "N={} -> {}: affinity jumped {:.4} (bound {:.4}) — closed={} -> {}",
                a.units,
                b.units,
                (b.affinity - a.affinity).abs(),
                bound,
                a.outer == 0,
                b.outer == 0
            );
        }
    }

    /// §7.1: mass increases down the table. With the original 0.012/0.004
    /// defects this failed at N=124, so the third closed shell — the element a
    /// designer would most want to include — weighed less than its predecessor.
    ///
    /// **Deliberately past 120.** Checked against the 120-element table this
    /// test passes with the broken constants, because the first inversion sits
    /// four elements beyond the end. It was written that way first; a mutation
    /// run caught it. 200 covers the fourth closure at 147 with margin.
    #[test]
    fn mass_increases_with_unit_count() {
        for pair in derive_table(200, 1.75, 0.42, 62.0, 0.55).windows(2) {
            assert!(
                pair[1].mass > pair[0].mass,
                "N={} mass {} >= N={} mass {}",
                pair[0].units,
                pair[0].mass,
                pair[1].units,
                pair[1].mass
            );
        }
    }

    /// The claim that made `Mass`'s fixed-point representation the right home
    /// for this series. It was false for **every** element before the defect
    /// factors were made dyadic, under a comment asserting it.
    #[test]
    fn every_mass_is_exact_on_the_fixed_point_grid() {
        for e in table() {
            let scaled = e.mass * 1024.0;
            assert_eq!(
                scaled,
                scaled.round(),
                "N={}: mass {} is not a multiple of 1/1024 (x1024 = {})",
                e.units,
                e.mass,
                scaled
            );
        }
    }

    /// Cumulative shell totals are 1, 13, 55, 147 — and the fill loop must
    /// agree with the closed form, since burial and valence both key on it.
    #[test]
    fn closures_land_where_the_shell_arithmetic_says() {
        let closed: Vec<usize> = table()
            .iter()
            .filter(|e| e.outer == 0)
            .map(|e| e.units)
            .collect();
        assert_eq!(closed, vec![1, 13, 55]);
    }

    /// Emergence's deletion test, pinned. The closed-shell family is produced
    /// by the general formula, not by a branch that names it — so there is
    /// nothing here to accidentally promote into a declaration.
    #[test]
    fn closed_shells_have_valence_zero_without_being_special_cased() {
        for e in table() {
            if e.outer == 0 {
                assert_eq!(e.valence, 0, "N={}", e.units);
            }
        }
    }
}
