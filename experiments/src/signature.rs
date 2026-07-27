//! Shape signatures, and the four ways of comparing two of them that G2 puts
//! side by side.
//!
//! A signature is the **support function** of the union of atom spheres,
//! sampled along the geodesic directions: along each direction it records how
//! far the molecule reaches. That is the quantity §8.2 calls radial extent.
//!
//! The four distances differ only in what they do about orientation, which is
//! the whole question:
//!
//! | Distance    | Orientation handled by                            |
//! |-------------|---------------------------------------------------|
//! | `D_frame`   | rotating into a frame anchored to canonical atoms |
//! | `D_raw`     | nothing — compared as embedded                    |
//! | `D_group`   | minimising over all 60 rotations                  |
//! | `D_sorted`  | discarding it, by sorting the extents             |
//!
//! `D_frame` is the current design. `D_raw` is the baseline it must beat.
//! `D_group` is what the binding kernel actually does, so it is the honest
//! upper bound on what a frame could ever buy. `D_sorted` is invariant under
//! *all* of SO(3) rather than just the 60, which costs information — two
//! genuinely different shapes can share a sorted extent profile — and the
//! measurement is whether that cost is worth paying.
//!
//! The two controls are not descriptors of shape at all, and that is the
//! point. They calibrate the instrument: see [`Descriptor::PositiveControl`]
//! and [`Descriptor::NegativeControl`].
#![expect(
    clippy::disallowed_methods,
    reason = "§13.4: support-function extents are sums of squares and radii — never NaN, \
                  never -0.0"
)]
#![allow(
    clippy::indexing_slicing,
    reason = "every index is a loop bound over D or over the molecule's own atom count"
)]

use crate::embed::{self, Point};
use crate::geodesic::{Geodesic, Mat3, N_ROTATIONS, apply_mat};
use crate::molecule::{ELEMENT_RADII, Molecule, N_ELEMENTS};

/// Sampling resolution. 42 is the spec's prior (§22.2); the ladder is
/// {12, 42, 162} and nothing between them preserves the exact permutation.
pub const D: usize = 42;

/// The support function of the union of atom spheres, sampled along `g.dirs`.
///
/// Written as an explicit max over atoms in index order rather than a
/// `fold`/`max_by`, so the traversal order is fixed and no float comparison
/// depends on iterator internals.
#[must_use]
pub fn signature<const N: usize>(g: &Geodesic<N>, coords: &[Point], elements: &[u8]) -> Vec<f64> {
    let mut out = vec![f64::NEG_INFINITY; N];
    for (i, u) in g.dirs.iter().enumerate() {
        let mut best = f64::NEG_INFINITY;
        for (a, p) in coords.iter().enumerate() {
            let radius = ELEMENT_RADII
                .get(usize::from(elements[a]))
                .copied()
                .unwrap_or(0.0);
            let reach = p[0] * u[0] + p[1] * u[1] + p[2] * u[2] + radius;
            if reach > best {
                best = reach;
            }
        }
        out[i] = best;
    }
    out
}

/// An orthonormal frame anchored to the canonical atom order.
///
/// The first canonical atom that is far enough from the centroid fixes the
/// first axis; the next one with a usable perpendicular component fixes the
/// second. This is the label-anchored construction G2 exists to test.
///
/// It is genuinely canonical — two numberings of the same molecule give the
/// same frame, which `the_canonical_frame_survives_relabelling` checks. What
/// it is not is *continuous*: adding one atom can change which atom comes
/// first in canonical order, and the frame then jumps. The identity matrix is
/// returned when no usable axis exists, which happens only for molecules too
/// small or too symmetric to define one.
#[must_use]
pub fn canonical_frame(coords: &[Point], order: &[usize]) -> Mat3 {
    const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let norm = |v: Point| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let unit = |v: Point| {
        let n = norm(v);
        [v[0] / n, v[1] / n, v[2] / n]
    };

    let Some(&first) = order.iter().find(|&&a| norm(coords[a]) > 1e-6) else {
        return IDENTITY;
    };
    let e0 = unit(coords[first]);

    // The second axis comes from the next canonical atom with a usable
    // component perpendicular to e0. "Usable" is a threshold, and a low one
    // would make the frame flip on numerical noise.
    let mut e1 = None;
    for &a in order {
        let p = coords[a];
        let d = p[0] * e0[0] + p[1] * e0[1] + p[2] * e0[2];
        let perp = [p[0] - d * e0[0], p[1] - d * e0[1], p[2] - d * e0[2]];
        if norm(perp) > 1e-3 {
            e1 = Some(unit(perp));
            break;
        }
    }
    let Some(e1) = e1 else { return IDENTITY };
    let e2 = [
        e0[1] * e1[2] - e0[2] * e1[1],
        e0[2] * e1[0] - e0[0] * e1[2],
        e0[0] * e1[1] - e0[1] * e1[0],
    ];
    [e0, e1, e2]
}

/// Everything a descriptor needs from one molecule, computed once.
///
/// Per-species, never per-molecule (§8.6) — the same discipline the real
/// simulation is held to, and the reason the harness can afford thousands of
/// trials.
#[derive(Debug, Clone)]
pub struct Profile {
    /// Signature from the embedding as produced, orientation and all.
    pub raw: Vec<f64>,
    /// Signature after rotating into the label-anchored canonical frame.
    pub framed: Vec<f64>,
    /// `raw`, sorted ascending — orientation discarded entirely.
    pub sorted: Vec<f64>,
    /// How many atoms of each element, for the positive control.
    pub composition: [f64; N_ELEMENTS],
    /// Atom count.
    pub atoms: f64,
    /// Digest of the canonical form, for the negative control.
    pub form_hash: u64,
}

/// Embed a molecule and derive every descriptor's view of it.
#[must_use]
pub fn profile<const N: usize>(g: &Geodesic<N>, m: &Molecule) -> Profile {
    let coords = embed::embed(m);
    let raw = signature(g, &coords, &m.elements);

    let frame = canonical_frame(&coords, &m.canonical_order());
    let framed_coords: Vec<Point> = coords.iter().map(|&p| apply_mat(&frame, p)).collect();
    let framed = signature(g, &framed_coords, &m.elements);

    let mut sorted = raw.clone();
    sorted.sort_by(f64::total_cmp);

    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "an atom count, far below 2^53 — exact"
    )]
    let atoms = m.len() as f64;
    let mut composition = [0.0; N_ELEMENTS];
    for &e in &m.elements {
        if let Some(slot) = composition.get_mut(usize::from(e)) {
            *slot += 1.0;
        }
    }
    Profile {
        raw,
        framed,
        sorted,
        composition,
        atoms,
        form_hash: m.form_hash(),
    }
}

/// Euclidean distance between two signatures, summed in index order.
#[must_use]
pub fn d_raw(a: &[f64], b: &[f64]) -> f64 {
    let mut total = 0.0;
    for i in 0..a.len().min(b.len()) {
        let e = a[i] - b[i];
        total += e * e;
    }
    total.sqrt()
}

/// Smallest Euclidean distance over all 60 rotations of `b`.
///
/// Rotation `r` carries direction `i` to `perms[r][i]`, so reading `b` through
/// the table is exactly comparing against a rotated copy. No antipode here:
/// `anti` is for binding, where two bodies touch along *opposite* directions.
/// This compares one shape against another in a different pose.
#[must_use]
pub fn d_group<const N: usize>(g: &Geodesic<N>, a: &[f64], b: &[f64]) -> f64 {
    let mut best = f64::INFINITY;
    for r in 0..N_ROTATIONS {
        let mut total = 0.0;
        for i in 0..N {
            let e = a[i] - b[usize::from(g.perms[r][i])];
            total += e * e;
        }
        if total < best {
            best = total;
        }
    }
    best.sqrt()
}

/// Euclidean distance between sorted extent profiles.
#[must_use]
pub fn d_sorted(a: &[f64], b: &[f64]) -> f64 {
    let mut x = a.to_vec();
    let mut y = b.to_vec();
    x.sort_by(f64::total_cmp);
    y.sort_by(f64::total_cmp);
    d_raw(&x, &y)
}

/// The four candidate descriptors, plus the two controls that decide whether
/// their numbers may be reported at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Descriptor {
    /// Label-anchored canonical frame — the current design.
    Frame,
    /// No alignment at all — the baseline the frame must beat.
    Raw,
    /// Minimum over the 60 rotations — what the binding kernel does.
    Group,
    /// Sorted extents: invariant under all of SO(3), not just the 60.
    Sorted,
    /// **Positive control.** L1 distance between per-element atom counts.
    ///
    /// Not a shape descriptor and not proposed as one — it uses no geometry
    /// whatsoever. It is a quantity that *must* show locality: a one-atom edit
    /// changes it by at most two, while two independently drawn molecules
    /// differ by far more. A harness that cannot see locality here cannot see
    /// locality at all, and its other numbers mean nothing.
    ///
    /// Composition rather than bare atom count, which is the obvious choice
    /// and is wrong: under the size-matched regime every molecule has the same
    /// atom count, so that control would read a flat zero and certify a broken
    /// instrument as working. Composition still discriminates when sizes are
    /// held equal.
    PositiveControl,
    /// **Negative control.** Hamming distance between digests of the canonical
    /// form. It identifies the species exactly, so it holds strictly more
    /// information than any signature — yet the hash destroys all similarity
    /// structure, so it *must* show no locality. A harness that reports
    /// locality here is reporting an artefact of its own statistics, and every
    /// other number it produced has to be thrown away.
    NegativeControl,
}

impl Descriptor {
    /// All six, in report order.
    pub const ALL: [Self; 6] = [
        Self::Frame,
        Self::Raw,
        Self::Group,
        Self::Sorted,
        Self::PositiveControl,
        Self::NegativeControl,
    ];

    /// Short label for the report.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Frame => "D_frame",
            Self::Raw => "D_raw",
            Self::Group => "D_group",
            Self::Sorted => "D_sorted",
            Self::PositiveControl => "control+ (composition)",
            Self::NegativeControl => "control- (form hash)",
        }
    }

    /// Whether this is a control rather than a candidate descriptor.
    #[must_use]
    pub const fn is_control(self) -> bool {
        matches!(self, Self::PositiveControl | Self::NegativeControl)
    }

    /// Distance between two profiles under this descriptor.
    #[must_use]
    pub fn distance<const N: usize>(self, g: &Geodesic<N>, a: &Profile, b: &Profile) -> f64 {
        match self {
            Self::Frame => d_raw(&a.framed, &b.framed),
            Self::Raw => d_raw(&a.raw, &b.raw),
            Self::Group => d_group(g, &a.raw, &b.raw),
            Self::Sorted => d_sorted(&a.sorted, &b.sorted),
            Self::PositiveControl => {
                let mut total = 0.0;
                for i in 0..N_ELEMENTS {
                    total += (a.composition[i] - b.composition[i]).abs();
                }
                total
            }
            #[allow(
                clippy::cast_precision_loss,
                reason = "a popcount of a u64: at most 64"
            )]
            Self::NegativeControl => f64::from((a.form_hash ^ b.form_hash).count_ones()),
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    clippy::many_single_char_names,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds; casts are bounded by \
              ELEMENT_RADII's length; g/m/a/b/d are this domain's vocabulary, not lazy naming"
)]
mod tests {
    use super::*;
    use crate::geodesic::{apply_mat, rotation_matrices};
    use crate::molecule::{ELEMENT_RADII, Molecule};
    use crate::rng::Stream;

    fn geo() -> Geodesic<D> {
        Geodesic::<D>::build().unwrap()
    }

    #[test]
    fn one_atom_at_the_origin_has_a_uniform_signature() {
        let g = geo();
        for (e, &radius) in ELEMENT_RADII.iter().enumerate() {
            let sig = signature(&g, &[[0.0, 0.0, 0.0]], &[e as u8]);
            for (i, &s) in sig.iter().enumerate() {
                assert!((s - radius).abs() < 1e-12, "direction {i}: {s} != {radius}");
            }
        }
    }

    /// The signature is the support function of the union of atom spheres:
    /// along each direction it reports how far the molecule reaches. Checked
    /// against a hand-computed value rather than against itself.
    #[test]
    fn the_signature_is_the_support_function() {
        let g = geo();
        let coords = [[3.0, 0.0, 0.0], [-3.0, 0.0, 0.0]];
        let sig = signature(&g, &coords, &[0, 3]);
        for (i, &s) in sig.iter().enumerate() {
            let u = g.dirs[i];
            // Spelled as a comparison, not `f64::max`, for the reason in
            // `min_pair_dist2` — and because a test that reaches for the
            // banned spelling is how the ban gets relaxed later.
            let (lo, hi) = (
                3.0 * u[0] + ELEMENT_RADII[0],
                -3.0 * u[0] + ELEMENT_RADII[3],
            );
            let want = if lo > hi { lo } else { hi };
            assert!((s - want).abs() < 1e-12, "direction {i}: {s} != {want}");
        }
    }

    /// The load-bearing link between geometry and the lookup table. Rotating
    /// the molecule by rotation `r` and recomputing the signature from scratch
    /// must give exactly the original signature read through `perms[r]`.
    ///
    /// If this fails, `D_group` is minimising over a table that does not
    /// correspond to rotating anything, and its number is meaningless while
    /// looking entirely plausible.
    #[test]
    fn rotating_the_molecule_permutes_the_signature() {
        let g = geo();
        let mut rng = Stream::new(31);
        let m = Molecule::random_tree(&mut rng, 12);
        let coords = crate::embed::embed(&m);
        let base = signature(&g, &coords, &m.elements);

        for (r, mat) in rotation_matrices().unwrap().iter().enumerate() {
            let turned: Vec<Point> = coords.iter().map(|&p| apply_mat(mat, p)).collect();
            let sig = signature(&g, &turned, &m.elements);
            for i in 0..D {
                let want = base[i];
                let got = sig[g.perms[r][i] as usize];
                assert!(
                    (got - want).abs() < 1e-9,
                    "rotation {r}, direction {i}: {got} != {want}"
                );
            }
        }
    }

    #[test]
    fn identical_molecules_are_at_distance_zero_under_every_descriptor() {
        let g = geo();
        let mut rng = Stream::new(32);
        let m = Molecule::random_tree(&mut rng, 14);
        let p = profile(&g, &m);
        for d in Descriptor::ALL {
            assert!(
                d.distance(&g, &p, &p).abs() < 1e-9,
                "{} put a molecule at nonzero distance from itself",
                d.name()
            );
        }
    }

    #[test]
    fn every_descriptor_is_symmetric() {
        let g = geo();
        let mut rng = Stream::new(33);
        let a = profile(&g, &Molecule::random_tree(&mut rng, 11));
        let b = profile(&g, &Molecule::random_tree(&mut rng, 13));
        for d in Descriptor::ALL {
            let ab = d.distance(&g, &a, &b);
            let ba = d.distance(&g, &b, &a);
            assert!(
                (ab - ba).abs() < 1e-9,
                "{} is asymmetric: {ab} vs {ba}",
                d.name()
            );
        }
    }

    /// `D_group` and `D_sorted` claim to be blind to orientation. That claim is
    /// what makes them candidates, so it is measured directly: rotate a
    /// molecule bodily and the distance to the original must stay zero.
    #[test]
    fn the_frame_free_descriptors_ignore_a_bodily_rotation() {
        let g = geo();
        let mut rng = Stream::new(34);
        let m = Molecule::random_tree(&mut rng, 15);
        let coords = crate::embed::embed(&m);
        let base = signature(&g, &coords, &m.elements);

        for (r, mat) in rotation_matrices().unwrap().iter().enumerate() {
            let turned: Vec<Point> = coords.iter().map(|&p| apply_mat(mat, p)).collect();
            let sig = signature(&g, &turned, &m.elements);
            assert!(
                d_group(&g, &base, &sig) < 1e-9,
                "D_group saw rotation {r} as a difference"
            );
            assert!(
                d_sorted(&base, &sig) < 1e-9,
                "D_sorted saw rotation {r} as a difference"
            );
        }
    }

    /// The contrast that makes the test above mean something. If `D_raw` were
    /// *also* rotation-blind, the comparison between them would be measuring
    /// nothing — so the harness needs it demonstrated that the unaligned
    /// descriptor really does see orientation.
    #[test]
    fn the_unaligned_descriptor_does_see_a_bodily_rotation() {
        let g = geo();
        let mut rng = Stream::new(35);
        let m = Molecule::random_tree(&mut rng, 15);
        let coords = crate::embed::embed(&m);
        let base = signature(&g, &coords, &m.elements);

        let moved = rotation_matrices()
            .unwrap()
            .iter()
            .map(|mat| {
                let turned: Vec<Point> = coords.iter().map(|&p| apply_mat(mat, p)).collect();
                d_raw(&base, &signature(&g, &turned, &m.elements))
            })
            .filter(|&d| d > 1e-6)
            .count();
        assert!(
            moved >= 50,
            "D_raw was unmoved by {} of 60 rotations — it is behaving as if aligned",
            60 - moved
        );
    }

    /// The canonical frame must not depend on atom numbering, or `D_frame` is
    /// measuring construction order rather than shape. This is the property
    /// the label-anchored design is *supposed* to have; G2 asks the separate
    /// question of whether having it is enough.
    #[test]
    fn the_canonical_frame_survives_relabelling() {
        let g = geo();
        let mut rng = Stream::new(36);
        for trial in 0..40 {
            let m = Molecule::random_tree(&mut rng, 6 + (trial % 12));
            let mut perm: Vec<usize> = (0..m.len()).collect();
            for i in (1..perm.len()).rev() {
                let j = rng.index(i + 1);
                perm.swap(i, j);
            }
            let relabelled = m.relabelled(&perm);
            let a = profile(&g, &m);
            let b = profile(&g, &relabelled);
            let d = d_raw(&a.framed, &b.framed);
            assert!(
                d < 1e-6,
                "trial {trial}: renumbering moved the framed signature by {d}"
            );
        }
    }

    #[test]
    fn the_positive_control_tracks_composition_and_the_negative_control_does_not() {
        let g = geo();
        let mut rng = Stream::new(37);
        let m = Molecule::random_tree(&mut rng, 10);
        let bigger = m.add_leaf(&mut rng);
        let a = profile(&g, &m);
        let b = profile(&g, &bigger);
        // One added atom moves exactly one element's count by one.
        assert!((Descriptor::PositiveControl.distance(&g, &a, &b) - 1.0).abs() < 1e-12);
        // A retype moves two counts by one each.
        let retyped = profile(&g, &m.retype_one(&mut rng));
        assert!((Descriptor::PositiveControl.distance(&g, &a, &retyped) - 2.0).abs() < 1e-12);
        // The negative control must react to *any* change, without its
        // magnitude meaning anything — that is what makes it a null.
        assert!(Descriptor::NegativeControl.distance(&g, &a, &b) > 0.0);
    }
}
