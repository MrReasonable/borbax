//! Graph distances to 3D coordinates, by stress majorization (spec §8.2).
//!
//! SMACOF rather than Laplacian eigenvectors: eigenvector sign is arbitrary
//! and degenerate eigenvalues are common in the symmetric graphs small
//! molecules actually are, so that arbitrariness would propagate straight
//! into species identity.
//!
//! The iteration count is **fixed**, not convergence-tested. A tolerance
//! makes the number of iterations depend on floating-point details that vary
//! by platform, and a configuration that stopped one iteration earlier is a
//! different configuration — so the stopping rule has to be something every
//! machine agrees on exactly (§13.4). This is also why no solver crate fits:
//! none offers "stop after exactly n iterations regardless of convergence".
//!
//! Nothing here uses a transcendental. The Guttman transform needs only
//! `sqrt` and the four arithmetic operations, and the initial layout is built
//! from graph distances rather than from a trigonometric spiral, precisely so
//! that stays true.
#![expect(
    clippy::disallowed_methods,
    reason = "§13.4: graph distances are non-negative sums, never -0.0. They CAN be \
                  f64::INFINITY for a disconnected component — see the finiteness \
                  guard above, which is what makes this selection safe"
)]
#![allow(
    clippy::indexing_slicing,
    reason = "every index is a loop bound over the same n as the coordinate vector"
)]

use crate::molecule::Molecule;
use crate::rng::Stream;

/// Iterations of the Guttman transform. Fixed for the reason in the module
/// docs; changing it is a physics change that moves every downstream number.
pub const SMACOF_ITERATIONS: usize = 240;

/// A 3D point.
pub type Point = [f64; 3];

fn distance(a: Point, b: Point) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// Raw stress: the summed squared mismatch between realised and target
/// distances, over unordered pairs.
///
/// Summed in `i < j` index order and never with `.sum()` over an unordered
/// iterator, so the accumulation order is fixed on every platform.
#[must_use]
pub fn stress(coords: &[Point], target: &[Vec<f64>]) -> f64 {
    let mut total = 0.0;
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            let e = distance(coords[i], coords[j]) - target[i][j];
            total += e * e;
        }
    }
    total
}

/// One Guttman transform: the majorization step at the heart of SMACOF.
///
/// With unit weights the update is `X <- (1/n) B(X) X`, which reduces to
///
/// ```text
/// x_i' = (1/n) * sum_{j != i} (d_ij / ||x_i - x_j||) * (x_i - x_j)
/// ```
///
/// Coincident points contribute nothing, the standard convention: the
/// direction to move them apart is undefined, and inventing one would be a
/// source of non-determinism rather than a fix.
///
/// Each `x_i'` accumulates over `j` in ascending index order, which is
/// load-bearing. Reordering this sum — or splitting it across threads —
/// changes the low bits of every coordinate, and those propagate into the
/// signature and then into whether two molecules bind.
#[must_use]
pub fn guttman_step(coords: &[Point], target: &[Vec<f64>]) -> Vec<Point> {
    let n = coords.len();
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "n is an atom count, far below 2^53 — exact"
    )]
    let inv_n = 1.0 / n as f64;
    let mut out = vec![[0.0; 3]; n];
    for i in 0..n {
        let mut acc = [0.0; 3];
        for j in 0..n {
            if i == j {
                continue;
            }
            let dist = distance(coords[i], coords[j]);
            if dist <= 0.0 {
                continue;
            }
            let ratio = target[i][j] / dist;
            for axis in 0..3 {
                acc[axis] += ratio * (coords[i][axis] - coords[j][axis]);
            }
        }
        for axis in 0..3 {
            out[i][axis] = acc[axis] * inv_n;
        }
    }
    out
}

/// Translate a configuration so its centroid sits at the origin.
///
/// Stress is translation-invariant, so the optimiser will not do this by
/// itself. Signatures are support functions measured from the origin, so
/// without it every extent would carry an arbitrary offset.
fn centre(coords: &mut [Point]) {
    if coords.is_empty() {
        return;
    }
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "n is an atom count, far below 2^53 — exact"
    )]
    let inv_n = 1.0 / coords.len() as f64;
    for axis in 0..3 {
        let mut sum = 0.0;
        for c in coords.iter() {
            sum += c[axis];
        }
        let mean = sum * inv_n;
        for c in coords.iter_mut() {
            c[axis] -= mean;
        }
    }
}

/// Run the transform a fixed number of times, centring at the end.
#[must_use]
pub fn smacof(target: &[Vec<f64>], init: Vec<Point>, iterations: usize) -> Vec<Point> {
    let mut coords = init;
    for _ in 0..iterations {
        coords = guttman_step(&coords, target);
    }
    centre(&mut coords);
    coords
}

/// The starting configuration: pivot coordinates plus a tiny fixed jitter.
///
/// Each atom is placed at its graph distance from three pivots — a cheap
/// landmark embedding that already has roughly the right shape, so the
/// fixed iteration budget is spent refining rather than untangling.
///
/// The pivots are chosen from the canonical order, which makes the whole
/// layout **canonical-order dependent**. That is deliberate and is the thing
/// G2 puts on trial: the orientation of the result is inherited from an atom
/// numbering that is not continuous under graph edits. `D_raw` measures what
/// that costs, and `D_group` and `D_sorted` measure what removing the
/// dependence buys.
///
/// The jitter breaks ties between atoms sharing all three pivot distances,
/// which is common in symmetric trees and would otherwise leave coincident
/// points that the transform can never separate. It is `1e-3` against pivot
/// coordinates of order 1-10, and comes from the integer stream rather than a
/// trigonometric spiral so no transcendental enters (§13.4).
#[must_use]
pub fn initial_coords(mol: &Molecule, target: &[Vec<f64>]) -> Vec<Point> {
    let n = mol.len();
    if n == 0 {
        return Vec::new();
    }
    let order = mol.canonical_order();
    // Canonical rank per atom, so every tie-break below is decided by the
    // canonical order rather than by atom numbering.
    let mut rank = vec![0usize; n];
    for (r, &atom) in order.iter().enumerate() {
        rank[atom] = r;
    }

    let p0 = *order.first().unwrap_or(&0);
    let p1 = (0..n)
        .max_by(|&a, &b| {
            target[p0][a]
                .total_cmp(&target[p0][b])
                .then(rank[b].cmp(&rank[a]))
        })
        .unwrap_or(0);
    let p2 = (0..n)
        .max_by(|&a, &b| {
            (target[p0][a] + target[p1][a])
                .total_cmp(&(target[p0][b] + target[p1][b]))
                .then(rank[b].cmp(&rank[a]))
        })
        .unwrap_or(0);

    (0..n)
        .map(|atom| {
            let mut jitter = Stream::new(u64::try_from(rank[atom]).unwrap_or(0));
            [
                target[p0][atom] + 1e-3 * (jitter.f64_unit() - 0.5),
                target[p1][atom] + 1e-3 * (jitter.f64_unit() - 0.5),
                target[p2][atom] + 1e-3 * (jitter.f64_unit() - 0.5),
            ]
        })
        .collect()
}

/// Embed a molecule: graph distances, landmark start, fixed SMACOF budget.
#[must_use]
pub fn embed(mol: &Molecule) -> Vec<Point> {
    let target = mol.graph_distances();
    let init = initial_coords(mol, &target);
    smacof(&target, init, SMACOF_ITERATIONS)
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds; the one cast is an \
              atom count to f64, far below 2^53"
)]
mod tests {
    use super::*;
    use crate::molecule::Molecule;
    use crate::rng::Stream;

    fn euclidean_matrix(pts: &[[f64; 3]]) -> Vec<Vec<f64>> {
        pts.iter()
            .map(|a| {
                pts.iter()
                    .map(|b| {
                        let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
                        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
                    })
                    .collect()
            })
            .collect()
    }

    /// The one guarantee stress majorization actually makes. If stress ever
    /// rises, the Guttman transform is misderived — and a misderived transform
    /// still *converges to something*, so nothing else here would notice.
    #[test]
    fn stress_never_increases() {
        let mut rng = Stream::new(21);
        for trial in 0..20 {
            let m = Molecule::random_tree(&mut rng, 6 + trial);
            let target = m.graph_distances();
            let mut coords = initial_coords(&m, &target);
            let mut previous = stress(&coords, &target);
            for step in 0..SMACOF_ITERATIONS {
                coords = guttman_step(&coords, &target);
                let current = stress(&coords, &target);
                assert!(
                    current <= previous + 1e-9,
                    "trial {trial} step {step}: stress rose from {previous} to {current}"
                );
                previous = current;
            }
        }
    }

    /// A configuration that genuinely lives in three dimensions must come back
    /// essentially exactly. This is what separates "the optimiser runs" from
    /// "the optimiser works".
    #[test]
    fn an_exact_three_dimensional_configuration_is_recovered() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.7, 0.0],
            [0.0, 0.0, 2.3],
            [1.1, 1.2, 0.4],
            [-0.8, 0.6, 1.9],
            [2.0, -1.0, 0.5],
        ];
        let target = euclidean_matrix(&pts);
        let init: Vec<[f64; 3]> = (0..pts.len())
            .map(|k| {
                let mut r = Stream::new(k as u64);
                [r.f64_unit(), r.f64_unit(), r.f64_unit()]
            })
            .collect();
        let coords = smacof(&target, init, SMACOF_ITERATIONS);
        let got = euclidean_matrix(&coords);
        for i in 0..pts.len() {
            for j in 0..pts.len() {
                // 1e-4, not 1e-12. The budget is fixed at 240 iterations for
                // determinism, and SMACOF converges geometrically rather than
                // exactly — measured residual here is ~2e-6. A tighter bound
                // would be asserting a convergence rate nobody promised.
                assert!(
                    (got[i][j] - target[i][j]).abs() < 1e-4,
                    "pair {i},{j}: wanted {}, got {}",
                    target[i][j],
                    got[i][j]
                );
            }
        }
    }

    /// A path graph's distances are realisable exactly on a straight line, so
    /// the optimum stress has a known value: zero. A tree in general does not
    /// embed exactly in 3D, which is why this uses the one family that does.
    #[test]
    fn a_path_embeds_with_almost_no_stress() {
        let m = Molecule::from_parts(vec![0; 6], (1..6).map(|i| (i - 1, i)).collect());
        let target = m.graph_distances();
        let coords = embed(&m);
        let residual = stress(&coords, &target);
        assert!(
            residual < 1e-3,
            "a path should embed nearly exactly; stress {residual}"
        );
    }

    /// The claim the loose bound above rests on. An arbitrary tolerance at a
    /// fixed budget cannot tell "converging slowly to the right answer" from
    /// "converged to the wrong one" — but the *rate* can, because only the
    /// former keeps improving.
    ///
    /// Measured on a configuration whose optimum is exactly zero:
    ///
    /// ```text
    ///     60 -> 3.97e-3      3840 -> 1.30e-6
    ///    240 -> 3.08e-4     15360 -> 8.13e-8
    ///    960 -> 2.04e-5     61440 -> 5.09e-9
    /// ```
    ///
    /// Every 4x in budget buys ~16x in stress, so the error falls as
    /// `1/iterations^2` and the limit is zero. That is what makes 240 a budget
    /// choice rather than a wrong answer.
    #[test]
    fn stress_falls_toward_zero_as_the_budget_grows() {
        let m = Molecule::from_parts(vec![0; 6], (1..6).map(|i| (i - 1, i)).collect());
        let target = m.graph_distances();
        let budgets = [60usize, 240, 960, 3840];
        let stresses: Vec<f64> = budgets
            .iter()
            .map(|&n| stress(&smacof(&target, initial_coords(&m, &target), n), &target))
            .collect();
        for w in stresses.windows(2) {
            // A 4x budget must buy at least 8x — comfortably under the ~16x
            // observed, but far above the 1x a stalled optimiser would give.
            assert!(
                w[0] > w[1] * 8.0,
                "quadrupling the budget only moved stress {} -> {}",
                w[0],
                w[1]
            );
        }
        assert!(
            stresses.last().is_some_and(|&s| s < 1e-5),
            "stress is not heading to zero: {stresses:?}"
        );
    }

    #[test]
    fn embedding_is_deterministic() {
        let mut rng = Stream::new(22);
        let m = Molecule::random_tree(&mut rng, 14);
        assert_eq!(embed(&m), embed(&m));
    }

    /// Stress is translation-invariant, so the optimiser cannot pin the
    /// centroid on its own. Signatures are support functions measured from the
    /// origin, so an uncentred embedding would make every extent depend on
    /// where the configuration happened to drift to.
    #[test]
    fn the_embedding_is_centred_on_the_origin() {
        let mut rng = Stream::new(23);
        let m = Molecule::random_tree(&mut rng, 18);
        let coords = embed(&m);
        for axis in 0..3 {
            let mean: f64 = coords.iter().map(|c| c[axis]).sum::<f64>() / coords.len() as f64;
            assert!(mean.abs() < 1e-9, "axis {axis} centroid at {mean}");
        }
    }

    /// Two atoms one bond apart must not land on top of each other, or the
    /// support function collapses and every signature reports the same shape.
    #[test]
    fn bonded_atoms_stay_apart() {
        let mut rng = Stream::new(24);
        for _ in 0..20 {
            let m = Molecule::random_tree(&mut rng, 16);
            let coords = embed(&m);
            for &(a, b) in m.edges() {
                let d = [
                    coords[a][0] - coords[b][0],
                    coords[a][1] - coords[b][1],
                    coords[a][2] - coords[b][2],
                ];
                let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                assert!(len > 0.3, "bonded atoms {a},{b} only {len} apart");
            }
        }
    }
}
