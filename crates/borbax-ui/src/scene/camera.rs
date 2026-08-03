//! The orbit camera: two angles and a distance about a fixed target.
//!
//! **No UI type, no GPU type, no physics.** This file is in the seam's strictest
//! tier — it names neither `egui`, `eframe`, `wgpu` nor `epaint` — which is what
//! lets every camera test run as plain arithmetic with no window and no adapter.
//!
//! # Why angles rather than a direction
//!
//! The camera **accumulates** yaw and pitch from drag deltas and never
//! reconstructs an orientation from a direction vector. Two reasons, and the
//! second is the one that would bite:
//!
//! - `det_math` has no `atan2`, and a camera needing one would be a design that
//!   reached for a wrapper rather than the other way round.
//! - Accumulating into the *matrix* instead of the parameter drifts, and makes
//!   every other camera test unreproducible — the camera stops being a pure
//!   function of anything a test can state.
//!
//! # §13.1
//!
//! Every transcendental here goes through [`borbax_units::det_math`]. The
//! viewer produces no simulation results, and that is deliberately **not** the
//! reason it is allowed to skip the chokepoint — it is not allowed to. The
//! camera maths is scheduled to move into `borbax-geometry`, whose bytes are
//! byte-compared across the §13.4 matrix, and trig's single-site rule already
//! exists. See the Step 3 design memo, R15.

use borbax_units::det_math;

/// A point or a direction in scene space.
///
/// The workspace spells three-component vectors this way — `borbax-molecule`'s
/// `Vec3` is the same alias, and `Embedding::coords()` hands back
/// `&[[f64; 3]]`. Repeating the shape rather than importing the alias keeps this
/// file free of a dependency it does not otherwise need.
pub type Vec3 = [f64; 3];

/// Which way is up before the camera has an opinion.
///
/// **A world constant, not a camera one.** The camera's own up vector is
/// derived; this is the reference that makes "up" mean anything at all, and the
/// pitch clamp exists to keep the camera's view direction from ever becoming
/// parallel to it.
const WORLD_UP: Vec3 = [0.0, 1.0, 0.0];

/// The camera's view-space axes, as unit vectors in scene space.
///
/// **`back`, not `forward`, and the name is doing real work.** These are the
/// rows of the view rotation in the standard right-handed convention: the
/// camera looks along **`-back`**. Spelling the field for the direction the
/// camera does *not* look reads oddly once and then stops a whole class of sign
/// error, because `[right, up, back]` is the triple whose determinant is `+1`
/// when the basis is right-handed. `[right, up, forward]` is left-handed for a
/// correct camera, so a reader checking the obvious thing about the obvious
/// name would conclude the scene was mirrored when it was not.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Basis {
    /// Screen-right, in scene space.
    pub right: Vec3,
    /// Screen-up, in scene space.
    pub up: Vec3,
    /// From the target toward the eye — the opposite of the view direction.
    pub back: Vec3,
}

/// An orbit camera: a target, a distance, and two accumulated angles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Rotation about [`WORLD_UP`], radians. Unbounded and deliberately not
    /// wrapped — see [`Camera::orbit`].
    yaw: f64,
    /// Elevation above the target's horizontal plane, radians. Clamped to
    /// [`PITCH_LIMIT`].
    pitch: f64,
    /// Distance from [`Camera::target`] to the eye. Always positive.
    distance: f64,
    /// The point the camera orbits and looks at.
    target: Vec3,
}

/// How near the pole the pitch may get, in radians.
///
/// **Strictly inside a quarter turn, and the margin is not cosmetic.** At
/// exactly `±π/2` the view direction is parallel to [`WORLD_UP`], their cross
/// product is the zero vector, and normalising it is `0/0` — a NaN camera, which
/// renders black. Clamping *inside* the pole also refuses the other failure the
/// obvious fix leaves open: a camera that guards the NaN but lets pitch pass the
/// pole flips `right`, and the whole image inverts while every vector stays
/// finite and unit-length.
///
/// The margin is `1e-3` rad — about 0.057°, far below what a drag can resolve on
/// screen and far above the `1e-8` where `cos(pitch)` starts losing the
/// precision `right` is built from.
const PITCH_LIMIT: f64 = std::f64::consts::FRAC_PI_2 - 1e-3;

/// Radians of rotation per logical point of drag.
///
/// A whole-screen drag across [`crate::WINDOW_SIZE`]'s 1200 points is about
/// twelve radians of yaw — a little under two full turns, which is enough to
/// look all the way round something without the mouse leaving the window and
/// slow enough to stop on a face.
const RADIANS_PER_POINT: f64 = 0.01;

impl Camera {
    /// A camera at `distance` from `target`, level with it and looking along
    /// `-Z`.
    #[must_use]
    pub const fn new(target: Vec3, distance: f64) -> Self {
        Self {
            yaw: 0.0,
            pitch: 0.0,
            distance,
            target,
        }
    }

    /// Turn by `dyaw` and `dpitch` **radians**, clamping pitch at the poles.
    ///
    /// The primitive [`Camera::orbit`] is written in terms of. Both exist
    /// because the pixels-to-radians conversion is a presentation choice and the
    /// rotation is not — a test that wants to state "a quarter turn" should not
    /// have to know the drag sensitivity to say it.
    ///
    /// **Yaw accumulates unwrapped.** Reducing it modulo `2π` would be free to
    /// write and would make the camera stop being a pure function of the drags
    /// applied to it: the wrap point is a discontinuity a test cannot state, and
    /// `sin`/`cos` are periodic anyway, so there is nothing to gain.
    pub fn turn(&mut self, dyaw: f64, dpitch: f64) {
        self.yaw += dyaw;
        // **`clamp` is permitted here and `min`/`max` would not be.** §13.1 bans
        // the latter two for their `±0.0` behaviour — they are free to return
        // either zero on a tie. `f64::clamp` is not on that list and is not the
        // same shape: it is specified as `if self < min { min } else if
        // self > max { max } else { self }`, which is a total order on the
        // comparisons rather than a tie-break, and it propagates NaN. It can
        // panic on `min > max` or a NaN bound; both bounds here are the same
        // positive compile-time constant negated, so neither is reachable.
        self.pitch = (self.pitch + dpitch).clamp(-PITCH_LIMIT, PITCH_LIMIT);
    }

    /// Turn by a drag of `dx`, `dy` **logical points**.
    ///
    /// Dragging right turns the molecule so its left side comes toward you, and
    /// dragging up tips its top toward you — the camera moves opposite the
    /// hand, which is what makes the object feel grabbed rather than the
    /// viewpoint pushed.
    pub fn orbit(&mut self, dx: f64, dy: f64) {
        self.turn(-dx * RADIANS_PER_POINT, dy * RADIANS_PER_POINT);
    }

    /// Where the eye is, in scene space.
    ///
    /// At `yaw = 0, pitch = 0` this is `target + [0, 0, distance]` — the camera
    /// on `+Z` looking toward `-Z`, which is the convention the projection is
    /// written against.
    ///
    /// **The accumulation order is pinned** (§13.1): each component is
    /// `target + distance * (trig * trig)`, with the two transcendentals
    /// multiplied together before the distance scales them. `mul_add` is not
    /// available here — it is a fused multiply-add, and whether it fuses depends
    /// on target features, so it is on `clippy.toml`'s ban list.
    #[must_use]
    pub fn eye(&self) -> Vec3 {
        let cos_pitch = det_math::cos(self.pitch);
        let sin_pitch = det_math::sin(self.pitch);
        let sin_yaw = det_math::sin(self.yaw);
        let cos_yaw = det_math::cos(self.yaw);
        [
            self.target[0] + self.distance * (cos_pitch * sin_yaw),
            self.target[1] + self.distance * sin_pitch,
            self.target[2] + self.distance * (cos_pitch * cos_yaw),
        ]
    }

    /// The camera's view-space axes.
    ///
    /// Orthonormal and right-handed for every reachable `(yaw, pitch)`, which is
    /// what the pitch clamp buys: `back` is never parallel to world-up, so the
    /// cross product that builds `right` is never degenerate. (Both are private
    /// constants in this module, so they are named rather than linked — a
    /// public doc may not link a private item, and `cargo doc` is a gate leg.)
    #[must_use]
    pub fn basis(&self) -> Basis {
        let eye = self.eye();
        let back = normalise(sub(eye, self.target));
        let right = normalise(cross(WORLD_UP, back));
        // `cross(back, right)`, in that order. `cross(right, back)` is the
        // negation, which leaves every vector unit-length and every pair
        // orthogonal while making the basis left-handed — the whole scene
        // mirrored, with orthonormality alone unable to see it.
        let up = cross(back, right);
        Basis { right, up, back }
    }
}

/// `a - b`, componentwise.
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The right-handed cross product `a x b`.
///
/// Written as `product - product` rather than `mul_add`: an FMA rounds once
/// where two operations round twice, so the two spellings give different bits,
/// and whether `mul_add` compiles to a fused instruction depends on target
/// features. `clippy.toml` bans it for exactly that (§13.1).
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `a . b`, summed in index order.
///
/// The order is load-bearing and pinned (§13.1). Nobody tidies it later.
fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `v` scaled to unit length.
///
/// **No zero guard, deliberately.** Every call site here is protected by
/// [`PITCH_LIMIT`], and a guard returning some arbitrary vector for the zero
/// case would convert a loud NaN into a silently wrong picture. If this ever
/// divides by zero the camera has reached a state the clamp was supposed to
/// exclude, and a black window is the correct report.
fn normalise(v: Vec3) -> Vec3 {
    let len = dot(v, v).sqrt();
    [v[0] / len, v[1] / len, v[2] / len]
}

#[cfg(test)]
mod tests {
    use super::{Camera, PITCH_LIMIT, RADIANS_PER_POINT, Vec3, cross, dot, sub};

    /// Every `(yaw, pitch)` a drag can reach, plus both clamp boundaries.
    fn poses() -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        for y in -8..=8 {
            for p in -8..=8 {
                out.push((f64::from(y) * 0.45, f64::from(p) * (PITCH_LIMIT / 8.0)));
            }
        }
        out
    }

    /// `worst = max(worst, x)`, written out because §13.1 bans `f64::max` for
    /// its `±0.0` behaviour. Every value here is a magnitude, so the tie case
    /// cannot arise — the explicit form is used anyway rather than arguing the
    /// exemption.
    fn keep_worst(worst: &mut f64, x: f64) {
        if x > *worst {
            *worst = x;
        }
    }

    /// Distance between two points.
    fn gap(a: Vec3, b: Vec3) -> f64 {
        let d = sub(a, b);
        dot(d, d).sqrt()
    }

    /// T1. The basis is orthonormal **and** right-handed.
    ///
    /// **The determinant arm is the one that carries the information.** A
    /// left-handed basis — `cross(right, back)` where the implementation writes
    /// `cross(back, right)` — is still three unit vectors, still mutually
    /// orthogonal, and mirrors the entire scene. Orthonormality cannot see it;
    /// `det == +1` is the only arm that can, and a mirrored render of a
    /// symmetric-looking molecule is exactly the defect nobody spots by looking.
    ///
    /// Both arms are kept because neither implies the other: the determinant of
    /// a *non*-orthonormal basis can still be positive.
    #[test]
    fn the_camera_basis_is_orthonormal_and_right_handed() {
        // Tolerance derived rather than chosen: each component is a handful of
        // f64 multiplies and adds off a `det_math` sine, so the error is a few
        // ulp of 1.0 — order 1e-16. 1e-14 is two orders of headroom and still
        // four orders tighter than any sign or ordering error, which lands at
        // 2.0 for the determinant and 1.0 for a dot product.
        const TOL: f64 = 1e-14;

        let mut checked = 0_u32;
        let mut worst_len = 0.0_f64;
        let mut worst_ortho = 0.0_f64;
        let mut worst_det = 0.0_f64;

        for (yaw, pitch) in poses() {
            for distance in [0.5_f64, 3.0, 40.0] {
                let mut cam = Camera::new([1.0, -2.0, 0.5], distance);
                cam.yaw = yaw;
                cam.pitch = pitch;
                let b = cam.basis();

                for v in [b.right, b.up, b.back] {
                    keep_worst(&mut worst_len, (dot(v, v).sqrt() - 1.0).abs());
                }
                for (u, v) in [(b.right, b.up), (b.up, b.back), (b.back, b.right)] {
                    keep_worst(&mut worst_ortho, dot(u, v).abs());
                }
                // `det [right, up, back]`, as the scalar triple product.
                keep_worst(
                    &mut worst_det,
                    (dot(cross(b.right, b.up), b.back) - 1.0).abs(),
                );
                checked += 1;
            }
        }

        // Bookkeeping: `checked` counts *poses*, which is what the loops
        // produce — 17 yaws x 17 pitches x 3 distances. Asserted as an equality
        // rather than a floor so a corpus that silently stops generating cases
        // fails here rather than passing three assertions over nothing.
        assert_eq!(checked, 17 * 17 * 3, "the corpus is not the one measured");

        assert!(
            worst_len <= TOL,
            "a basis vector was {worst_len} off unit length"
        );
        assert!(
            worst_ortho <= TOL,
            "two basis vectors were {worst_ortho} from orthogonal"
        );
        assert!(
            worst_det <= TOL,
            "det[right, up, back] was {worst_det} from +1. A determinant near -1 is a \
             left-handed basis: every vector is still unit-length and still \
             orthogonal, and the entire scene is mirrored"
        );
    }

    /// T2. A quarter turn moves the camera; a full turn brings it back.
    ///
    /// **The negative arm is the one that carries the information**, and the
    /// order below reflects that. "2π returns to the start" alone passes on a
    /// camera that ignores its argument entirely — the most complete failure
    /// available — so it is not a test on its own. Pairing it with "π/2 does
    /// *not* return to the start" is what pins that the angle is used, and used
    /// in radians: a camera reading its input as degrees turns 2π ≈ 6.3° (fails
    /// the first arm) and π/2 ≈ 1.6° (fails the second, which is the arm that
    /// says how far it is wrong rather than merely that it is).
    #[test]
    fn a_quarter_turn_is_not_a_full_turn() {
        const DISTANCE: f64 = 4.0;
        let target = [0.5, -1.0, 2.0];

        let start = Camera::new(target, DISTANCE).eye();

        let mut full = Camera::new(target, DISTANCE);
        full.turn(std::f64::consts::TAU, 0.0);
        assert!(
            gap(full.eye(), start) <= 1e-14,
            "a full turn left the eye {} from where it started",
            gap(full.eye(), start)
        );

        let mut quarter = Camera::new(target, DISTANCE);
        quarter.turn(std::f64::consts::FRAC_PI_2, 0.0);
        // A quarter turn about a circle of radius `DISTANCE` moves the eye by
        // the chord `r * sqrt(2)`. Asserting the *value* rather than merely
        // "it moved" is what makes this arm fail on a degrees reading, which
        // moves by 0.11 rather than 5.66.
        let expected_chord = DISTANCE * std::f64::consts::SQRT_2;
        let moved = gap(quarter.eye(), start);
        assert!(
            (moved - expected_chord).abs() <= 1e-14,
            "a quarter turn moved the eye {moved}, expected the chord \
             {expected_chord}. A camera reading radians as degrees moves ~0.11 here, \
             and passes any test that only asks whether the eye moved at all"
        );

        // And the sensitivity constant is wired the way `orbit`'s doc says: a
        // drag right turns the camera the other way, so the object follows the
        // hand. Without this, `turn` could be perfect and `orbit` inverted.
        let mut dragged = Camera::new(target, DISTANCE);
        dragged.orbit(10.0, 0.0);
        let mut turned = Camera::new(target, DISTANCE);
        turned.turn(-10.0 * RADIANS_PER_POINT, 0.0);
        assert!(
            gap(dragged.eye(), turned.eye()) <= 1e-15,
            "a rightward drag did not turn the camera the way `orbit` documents"
        );
    }

    /// T8. Pitch cannot reach the pole, and cannot pass it either.
    ///
    /// Two failures, and the second is the one the obvious fix leaves open.
    ///
    /// At exactly `±π/2` the view direction is parallel to `WORLD_UP`, the cross
    /// product that builds `right` is the zero vector, and normalising it gives
    /// NaN — a black window. That much a NaN guard catches.
    ///
    /// What it does not catch: letting pitch pass the pole. Every vector stays
    /// finite, unit-length and orthogonal, and the determinant stays `+1` —
    /// what flips is **`right`**, and with it the whole image, left for right,
    /// discontinuously and in one frame.
    ///
    /// **The obvious arm for that is vacuous, and it took a mutation to find
    /// out.** The first version of this test asserted `dot(up, WORLD_UP) > 0`,
    /// reasoning that passing the pole tips the view upside down. It cannot
    /// fire: working the construction through, `right` normalises to
    /// `sign(cos pitch) * [cos yaw, 0, -sin yaw]` and `up`'s vertical component
    /// is therefore `|cos pitch|`, which is never negative. The assertion read
    /// as coverage of the exact failure the clamp exists to prevent while being
    /// unable to fail — and it passed a `PITCH_LIMIT` widened to `π/2 + 0.5`,
    /// which is precisely the "guards the NaN, permits the flip" fix.
    ///
    /// The arm that works falls out of the same algebra: within the clamp
    /// `cos pitch > 0`, so **`right` does not depend on pitch at all**. Compare
    /// it against the same yaw at pitch zero and a pole crossing negates it.
    #[test]
    fn the_camera_does_not_degenerate_at_the_poles() {
        const TARGET: Vec3 = [0.0, 0.0, 0.0];
        let mut checked = 0_u32;

        // Far past both poles, in steps a drag could actually deliver, and from
        // both directions so a one-sided clamp fails.
        for direction in [1.0_f64, -1.0] {
            let mut cam = Camera::new(TARGET, 3.0);
            for _ in 0..40 {
                cam.turn(0.3, direction * 0.2);
                let b = cam.basis();

                assert!(
                    b.right
                        .iter()
                        .chain(&b.up)
                        .chain(&b.back)
                        .all(|c| c.is_finite()),
                    "the basis went non-finite at pitch {}: the cross product that \
                     builds `right` degenerated, which is a black window",
                    cam.pitch
                );
                assert!(
                    cam.pitch.abs() <= PITCH_LIMIT,
                    "pitch reached {}, outside the clamp at {PITCH_LIMIT}",
                    cam.pitch
                );

                // The same yaw, level with the target — `right` must be
                // identical, because within the clamp it is a function of yaw
                // alone.
                let mut level = Camera::new(TARGET, 3.0);
                level.yaw = cam.yaw;
                assert!(
                    dot(b.right, level.basis().right) > 0.0,
                    "`right` reversed at pitch {} — the view has passed the pole and \
                     the image is mirrored left for right. Every vector here is \
                     finite, unit-length and orthogonal and the determinant is still \
                     +1, so this is the only arm that can see it",
                    cam.pitch
                );
                checked += 1;
            }
        }

        // Bookkeeping: 2 directions x 40 turns, counted where a basis is
        // examined rather than at loop entry.
        assert_eq!(checked, 80, "the corpus is not the one measured");
    }
}
