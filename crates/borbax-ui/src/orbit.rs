//! The orbit camera: two angles and a distance, about the origin.
//!
//! **No engine type and no UI type**, which is what lets every test here run as
//! plain arithmetic with no window and no device.
//!
//! # What this owns, and what it deliberately does not
//!
//! It owns `(yaw, pitch, distance) -> eye position`. Everything after that —
//! the view basis, the projection matrix, the perspective divide, the depth
//! range — belongs to the engine, and re-deriving any of it here would be
//! duplicating a dependency in `f32`. `Transform::look_to` builds the same basis
//! the predecessor of this file built by hand, in the same order and with the
//! same handedness.
//!
//! # Why angles rather than a direction
//!
//! Yaw and pitch **accumulate**; an orientation is never read back out of a
//! transform to be modified. Accumulating into the matrix drifts, and worse, it
//! stops the camera being a pure function of anything a test can state.
//!
//! # The target is the origin, and that is a fact rather than a choice
//!
//! `Embedding` is centred — `layout.rs` says so where it takes the radius of
//! gyration about the origin — so the molecule's centroid *is* `[0, 0, 0]`. A
//! `target` field would be a second encoding of that, and the viewer computing a
//! centroid would be the viewer deriving a value the chemistry already gives.
//!
//! # §13.1
//!
//! Every transcendental goes through [`borbax_units::det_math`]. The viewer
//! produces no simulation results, and that is deliberately **not** why this is
//! allowed to be sloppy — it is not allowed to be. The rule is about which
//! function is called, and the workspace has exactly one place to call it from.

use borbax_units::det_math;

/// A point or a direction, in the same space `Embedding::coords` uses.
pub type Vec3 = [f64; 3];

/// Vertical field of view, radians.
///
/// **Ours rather than the engine's default, so a test can name it** without the
/// strictest seam tier having to name the engine. It is set on the projection at
/// spawn time; leaving it implicit would make every framing test a test of
/// whatever the engine's default happened to be this release.
pub const FOV_Y: f64 = std::f64::consts::FRAC_PI_4;

/// Near clipping distance, in the same units as the embedding.
///
/// Same argument as [`FOV_Y`]: pinned here so the no-clipping test is about a
/// number this crate chose.
pub const NEAR: f64 = 0.1;

/// How near the pole the pitch may get, radians.
///
/// **The reason recorded for this constant's predecessor was wrong, and the
/// correction matters because it changes what the test has to look at.** That
/// version said the pole makes the up-vector cross product `0/0` — "a NaN
/// camera, which renders black". Measured against the engine actually in use:
/// at exactly `±pi/2` the basis comes back finite and orthonormal, because the
/// normalisation falls back to an arbitrary orthogonal vector rather than
/// dividing by zero.
///
/// What actually happens past the pole is that **screen-right flips through
/// 180 degrees** and the whole image turns over between one frame and the next.
/// Measured per `5e-4` rad step of pitch: `6.7e-8` inside the clamp against
/// `1.999999762` across the pole — seven orders of magnitude, which is what the
/// continuity test discriminates on.
///
/// The margin is `1e-3` rad, about 0.057 degrees: far below what a drag can
/// resolve on screen, and far above where `cos(pitch)` starts losing precision.
const PITCH_LIMIT: f64 = std::f64::consts::FRAC_PI_2 - 1e-3;

/// Radians of rotation per logical point of drag.
const RADIANS_PER_POINT: f64 = 0.01;

/// Multiples of the radius of gyration the camera stands back by.
///
/// **Scaled from a shipped measurement rather than from a bounding radius the
/// viewer would have to compute**, and then *chosen by measurement* rather than
/// by the derivation looking right — which matters, because the derivation on
/// its own would have been wrong by a factor of two.
///
/// The molecule exactly fills the vertical field at
/// `bounding_radius / sin(FOV_Y / 2)`, which is `2.613 * bounding_radius`.
/// Measured across 500 universes, `bounding / gyration` runs **1.6491 ..
/// 1.7697**, so "just touching the edges" is **4.3092 .. 4.6243** gyration
/// radii. This constant is 6.5, which leaves the molecule spanning about
/// **66-71%** of the view's height in every universe measured: filling it edge
/// to edge leaves no room to orbit into, and the first value written here — 9.0,
/// a guess — put it twice as far away as it needed to be.
///
/// Vertical is the binding direction because the scene region is wider than it
/// is tall. If that ever stops being true, this constant is measuring the wrong
/// axis.
const FRAMING: f64 = 6.5;

/// How far in and out of the framing distance the wheel may go.
const ZOOM_IN_LIMIT: f64 = 0.35;
/// See [`ZOOM_IN_LIMIT`].
const ZOOM_OUT_LIMIT: f64 = 4.0;
/// Distance multiplier per notch of wheel.
const ZOOM_PER_NOTCH: f64 = 0.88;

/// A gyration radius this small means the molecule has no size to frame.
///
/// Reached only by a degenerate embedding — every atom on one point — which the
/// demo molecule cannot produce. It exists so the framing distance cannot come
/// out as zero and put the camera inside the scene.
const MIN_GYRATION: f64 = 1e-6;

/// Two angles and a distance, plus whether a drag currently owns the pointer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orbit {
    /// Rotation about world-up, radians. Unbounded and deliberately unwrapped.
    yaw: f64,
    /// Elevation, radians, clamped to [`PITCH_LIMIT`].
    pitch: f64,
    /// Distance from the origin to the eye. Always positive.
    distance: f64,
    /// The distance this camera was framed at, which the zoom limits scale.
    framed_at: f64,
    /// Whether the press that is currently down began over the scene.
    dragging: bool,
    /// Last frame's button state, for edge detection.
    was_pressed: bool,
}

impl Default for Orbit {
    /// A camera framing a molecule of unit gyration radius.
    ///
    /// Exists so the engine can hold this in a resource before any universe has
    /// been laid out. The first molecule replaces it.
    fn default() -> Self {
        Self::framing(1.0)
    }
}

impl Orbit {
    /// A camera standing back far enough to see a molecule of this gyration
    /// radius, level with it.
    #[must_use]
    pub fn framing(gyration: f64) -> Self {
        let size = if gyration > MIN_GYRATION {
            gyration
        } else {
            MIN_GYRATION
        };
        let distance = FRAMING * size;
        Self {
            yaw: 0.0,
            pitch: 0.0,
            distance,
            framed_at: distance,
            dragging: false,
            was_pressed: false,
        }
    }

    /// Turn by `dyaw` and `dpitch` **radians**, clamping pitch short of the
    /// poles.
    ///
    /// **Yaw accumulates unwrapped.** Reducing it modulo `2*pi` is free to write
    /// and would put a discontinuity in the middle of the parameter a test
    /// states its case in; `sin` and `cos` are periodic, so there is nothing to
    /// gain.
    pub fn turn(&mut self, dyaw: f64, dpitch: f64) {
        self.yaw += dyaw;
        // `clamp`, not `min`/`max`. §13.1 bans the latter two for being free to
        // return either zero on a `+0.0`/`-0.0` tie; `clamp` is specified as a
        // chain of comparisons and propagates NaN. It can panic on `min > max`
        // or a NaN bound, and both bounds here are one positive compile-time
        // constant and its negation.
        self.pitch = (self.pitch + dpitch).clamp(-PITCH_LIMIT, PITCH_LIMIT);
    }

    /// Turn by a drag of `dx`, `dy` **logical points**.
    ///
    /// Dragging right brings the molecule's left side toward you and dragging up
    /// tips its top toward you — the camera moves against the hand, which is
    /// what makes the object feel grabbed rather than the viewpoint shoved.
    pub fn drag(&mut self, dx: f64, dy: f64) {
        self.turn(-dx * RADIANS_PER_POINT, dy * RADIANS_PER_POINT);
    }

    /// Move `notches` of wheel closer, clamped either side of the framing
    /// distance.
    ///
    /// **Multiplicative, so a notch means the same thing at every distance.**
    /// Additive zoom crawls when far away and jumps through the molecule when
    /// near.
    pub fn zoom(&mut self, notches: f64) {
        let scaled = self.distance * det_math::powf(ZOOM_PER_NOTCH, notches);
        self.distance = scaled.clamp(
            self.framed_at * ZOOM_IN_LIMIT,
            self.framed_at * ZOOM_OUT_LIMIT,
        );
    }

    /// Feed one frame of pointer state in, and orbit if this drag is ours.
    ///
    /// **The latch is the whole point of this function.** A drag belongs to the
    /// scene if the *press* landed on the scene; once it has, it keeps the
    /// pointer until release, however far outside the cursor wanders. Deciding
    /// per frame instead — orbit whenever the cursor is not over the panel —
    /// makes the molecule freeze the instant a drag crosses onto the periodic
    /// table, which is the failure a design written from memory ships.
    ///
    /// `ui_wants_pointer` is asked only on the press edge for that reason.
    pub fn point(&mut self, pressed: bool, ui_wants_pointer: bool, dx: f64, dy: f64) {
        if pressed && !self.was_pressed {
            self.dragging = !ui_wants_pointer;
        }
        if !pressed {
            self.dragging = false;
        }
        self.was_pressed = pressed;
        if self.dragging {
            self.drag(dx, dy);
        }
    }

    /// Where the eye is.
    ///
    /// At `yaw = 0, pitch = 0` this is `[0, 0, distance]` — on `+Z`, looking
    /// toward `-Z`.
    ///
    /// **The accumulation order is pinned** (§13.1): each component is
    /// `distance * (trig * trig)`, the two transcendentals multiplied before the
    /// distance scales them. `mul_add` is not an option — whether it fuses
    /// depends on target features, which is why it is on the ban list.
    #[must_use]
    pub fn eye(&self) -> Vec3 {
        let cos_pitch = det_math::cos(self.pitch);
        let sin_pitch = det_math::sin(self.pitch);
        let sin_yaw = det_math::sin(self.yaw);
        let cos_yaw = det_math::cos(self.yaw);
        [
            self.distance * (cos_pitch * sin_yaw),
            self.distance * sin_pitch,
            self.distance * (cos_pitch * cos_yaw),
        ]
    }

    /// How far the eye is from the origin.
    #[must_use]
    pub const fn distance(&self) -> f64 {
        self.distance
    }

    /// Whether a drag currently belongs to the scene.
    #[must_use]
    pub const fn dragging(&self) -> bool {
        self.dragging
    }
}

#[cfg(test)]
mod measure {
    use crate::molecule::build;
    use borbax_universe::Universe;

    /// How far the molecule actually reaches, against the quantity the framing
    /// distance is scaled from. Printed, not asserted.
    #[test]
    #[ignore = "measurement, not a check"]
    fn bounding_against_gyration() {
        let mut ratio_lo = f64::INFINITY;
        let mut ratio_hi = 0.0f64;
        let mut bound_hi = 0.0f64;
        let seeds = 500u64;
        for seed in 0..seeds {
            let u = Universe::generate(seed);
            let Some(demo) = build(&u) else { continue };
            let mut bounding = 0.0f64;
            for (p, r) in demo.embedding.coords().iter().zip(demo.embedding.radii()) {
                let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() + r.get();
                if d > bounding {
                    bounding = d;
                }
            }
            let g = demo.embedding.radius_of_gyration().get();
            let ratio = bounding / g;
            if ratio < ratio_lo {
                ratio_lo = ratio;
            }
            if ratio > ratio_hi {
                ratio_hi = ratio;
            }
            if bounding > bound_hi {
                bound_hi = bounding;
            }
        }
        println!("bounding/gyration: {ratio_lo:.4} .. {ratio_hi:.4}");
        println!("max bounding radius: {bound_hi:.4}");
        println!("needed framing (gyration units) = 2.613 * bounding/gyration");
        println!("  => {:.4} .. {:.4}", 2.613 * ratio_lo, 2.613 * ratio_hi);
    }
}

#[cfg(test)]
mod tests {
    use super::{FOV_Y, NEAR, Orbit};
    use crate::molecule::build;
    use borbax_universe::Universe;

    /// Every claim below that ranges over universes uses this corpus.
    const SEEDS: u64 = 500;

    /// Distance from the origin to `p`.
    fn norm(p: [f64; 3]) -> f64 {
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    }

    /// A position as raw bits, for the places where the claim really is
    /// *identical* rather than *close*.
    ///
    /// Two cameras built and driven separately are two different computations,
    /// so this is not the vacuous shape where both sides of an equality read the
    /// same expression. Every caller also asserts finiteness first, because
    /// `NaN` bits compare equal to `NaN` bits and would make the comparison pass
    /// over a broken camera.
    fn bits(p: [f64; 3]) -> [u64; 3] {
        [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
    }

    /// A quarter turn moves the eye; a full turn brings it back.
    ///
    /// **Both arms, because either alone is vacuous.** "A full turn returns to
    /// the start" passes on a camera that ignores the angle entirely. "A quarter
    /// turn moves it" passes on degrees-vs-radians, which moves it by the wrong
    /// amount. Together they pin the unit: a quarter turn of a camera at
    /// distance `d` moves the eye by the chord `d * sqrt(2)`, and 0.25 *degrees*
    /// would move it by 0.0044 `d`.
    #[test]
    fn a_quarter_turn_is_not_a_full_turn() {
        let mut cam = Orbit::framing(1.0);
        let start = cam.eye();
        let d = cam.distance();

        cam.turn(std::f64::consts::FRAC_PI_2, 0.0);
        let quarter = cam.eye();
        let chord = norm([
            quarter[0] - start[0],
            quarter[1] - start[1],
            quarter[2] - start[2],
        ]);
        let expected = d * std::f64::consts::SQRT_2;
        assert!(
            (chord - expected).abs() < 1e-9,
            "a quarter turn moved the eye {chord:.6}, not the chord {expected:.6} — \
             the angle is being read in the wrong unit"
        );

        cam.turn(3.0 * std::f64::consts::FRAC_PI_2, 0.0);
        let full = cam.eye();
        assert!(
            norm([full[0] - start[0], full[1] - start[1], full[2] - start[2]]) < 1e-9,
            "a full turn did not come back to where it started"
        );
    }

    /// Pitch stops short of the pole, however far the drag goes.
    ///
    /// **Two earlier versions of this test could not fail, in different ways,
    /// and the second is the instructive one.**
    ///
    /// The first asserted on the *up* vector. That cannot fail: the engine
    /// re-derives up from the view direction, so it comes back positive on both
    /// sides of the pole. This project shipped that version once.
    ///
    /// The second asserted that the eye's *horizontal reach* never gets small.
    /// That reads as obviously right — at the pole the eye stands directly over
    /// the molecule and the reach really is zero — and it still cannot catch the
    /// defect, because **past** the pole the reach grows again. A limit widened
    /// to `pi/2 + 0.5` leaves the eye 0.479 of a radius out from the axis, which
    /// clears any floor comfortably while the image is upside down. It also
    /// stepped pitch by 0.5 rad, so the sweep jumped over the pole without ever
    /// sampling near it. Both faults found by running the mutation, not by
    /// reading it.
    ///
    /// What is claimed here is the **sign**. Yaw stays at zero, so the eye sits
    /// on `+Z` and `cos(pitch)` is what holds it there; crossing the pole makes
    /// that negative, which is precisely when screen-right flips through 180
    /// degrees and the image turns over between frames. The floor beneath it is
    /// a separate claim — that the basis never becomes degenerate — and steps
    /// are 0.05 rad so the sweep converges onto the clamp instead of stepping
    /// across it.
    #[test]
    fn the_camera_never_reaches_the_pole() {
        for direction in [1.0, -1.0] {
            let mut cam = Orbit::framing(1.0);
            let floor = cam.distance() * 5e-4;
            let mut samples = 0u32;
            for _ in 0..4000 {
                cam.drag(0.0, direction * 5.0);
                let eye = cam.eye();
                // Yaw stays at zero through this sweep, so the eye's `z` is
                // `distance * cos(pitch)` and its **sign** is the whole claim.
                assert!(
                    eye[2] > 0.0,
                    "the eye crossed the pole — at yaw 0 it is now at z = {:.6}, so \
                     screen-right has flipped and the image is upside down",
                    eye[2]
                );
                assert!(
                    eye[2] > floor,
                    "the eye came within {:.3e} of standing directly over the \
                     molecule, so the view basis is degenerate",
                    eye[2]
                );
                samples += 1;
            }
            assert_eq!(samples, 4000, "the sweep did not run");
        }
    }

    /// The camera is a pure function of the drags applied to it.
    ///
    /// **The obvious spelling of this test was wrong and passed**, which is the
    /// worse of the two ways to be wrong. It compared one whole drag against two
    /// half drags — that is an assertion that floating-point addition is
    /// *associative*, which it is not. It happened to hold for the constants in
    /// front of me, and would have broken later for a reason nobody could read
    /// off the test.
    ///
    /// What is claimed instead is two things that are actually true of a camera
    /// that accumulates into its angles:
    ///
    /// - the same sequence of drags gives the same camera, **bit for bit**;
    /// - turning by an exactly-representable angle and back returns **exactly**
    ///   to the start. This is the arm that fails on accumulating into a
    ///   rotation matrix instead — `R * R⁻¹` is not the identity in floating
    ///   point, and the drift is invisible in any single frame.
    #[test]
    fn the_camera_is_a_pure_function_of_its_parameters() {
        let steps = [(37.0, 11.0), (-91.0, 4.0), (13.0, -55.0), (220.0, 3.0)];
        let mut left = Orbit::framing(1.2);
        let mut right = Orbit::framing(1.2);
        for _ in 0..25 {
            for (dx, dy) in steps {
                left.drag(dx, dy);
                right.drag(dx, dy);
            }
        }
        assert!(
            left.eye().iter().all(|c| c.is_finite()),
            "the camera went non-finite, which a bitwise comparison would happily \
             call equal to another non-finite camera"
        );
        assert_eq!(
            bits(left.eye()),
            bits(right.eye()),
            "two cameras given the same drags disagreed"
        );

        // 0.5 and 0.25 are exact in binary, so `+x` then `-x` returns the angle
        // to exactly its starting value — for a camera that stores the angle.
        let mut there_and_back = Orbit::framing(1.2);
        let start = bits(there_and_back.eye());
        for _ in 0..64 {
            there_and_back.turn(0.5, 0.25);
            there_and_back.turn(-0.5, -0.25);
        }
        assert_eq!(
            bits(there_and_back.eye()),
            start,
            "turning back and forth drifted, so the orientation is being \
             accumulated somewhere other than the angle"
        );
    }

    /// The camera moves against the hand.
    ///
    /// Fails on an inverted sign or a swapped axis, which is the difference
    /// between the object feeling grabbed and the viewpoint feeling shoved.
    #[test]
    fn dragging_moves_the_camera_opposite_the_hand() {
        let mut right = Orbit::framing(1.0);
        right.drag(50.0, 0.0);
        assert!(
            right.eye()[0] < 0.0,
            "dragging right did not swing the eye left"
        );

        let mut up = Orbit::framing(1.0);
        up.drag(0.0, 50.0);
        assert!(up.eye()[1] > 0.0, "dragging up did not raise the eye");
    }

    /// The wheel stops, and it stops at the limit rather than wherever it
    /// happens to be.
    #[test]
    fn zooming_stops_at_both_ends() {
        let mut cam = Orbit::framing(1.0);
        let framed = cam.distance();
        for _ in 0..500 {
            cam.zoom(1.0);
        }
        let near = cam.distance();
        for _ in 0..1000 {
            cam.zoom(-1.0);
        }
        let far = cam.distance();
        assert!(near < framed, "zooming in did not come closer");
        assert!(far > framed, "zooming out did not go further");
        assert!(
            near > 0.0 && far.is_finite(),
            "the zoom ran away: {near} .. {far}"
        );
    }

    /// Zoomed all the way in, the camera is still outside the molecule and the
    /// molecule is still in front of the near plane.
    ///
    /// **This is the claim the framing constant is not allowed to make by
    /// argument.** Fails on an unclamped zoom, on a framing distance tuned to
    /// one seed, and on a near plane raised without re-checking the zoom limit —
    /// each of which shows up as the molecule vanishing, turning inside out, or
    /// swallowing the camera.
    #[test]
    fn zooming_all_the_way_in_never_crosses_the_near_plane() {
        let mut checked = 0u64;
        let mut worst = f64::INFINITY;
        for seed in 0..SEEDS {
            let universe = Universe::generate(seed);
            let demo =
                build(&universe).unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            let mut bounding = 0.0f64;
            for (p, r) in demo.embedding.coords().iter().zip(demo.embedding.radii()) {
                let reach = norm(*p) + r.get();
                if reach > bounding {
                    bounding = reach;
                }
            }
            let mut cam = Orbit::framing(demo.embedding.radius_of_gyration().get());
            for _ in 0..500 {
                cam.zoom(1.0);
            }
            let clearance = cam.distance() - bounding;
            if clearance < worst {
                worst = clearance;
            }
            assert!(
                clearance > NEAR,
                "seed {seed}: fully zoomed in, the nearest surface is {clearance:.4} \
                 from the eye against a near plane at {NEAR}, so the molecule clips"
            );
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was measured over"
        );
        assert!(
            worst.is_finite(),
            "no clearance was measured, so this test asserted nothing"
        );
    }

    /// The whole molecule is inside the field of view at the opening camera.
    ///
    /// **Fails on a framing distance tuned to one seed**, which is the shape the
    /// first value of `FRAMING` had. Checks the half-angle directly rather than
    /// through a projection matrix, so it is about the constant this crate owns
    /// and not about the engine's arithmetic.
    #[test]
    fn every_atom_is_within_the_field_of_view_at_the_opening_camera() {
        let half = FOV_Y / 2.0;
        let sin_half = borbax_units::det_math::sin(half);
        for seed in 0..SEEDS {
            let universe = Universe::generate(seed);
            let demo =
                build(&universe).unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            let mut bounding = 0.0f64;
            for (p, r) in demo.embedding.coords().iter().zip(demo.embedding.radii()) {
                let reach = norm(*p) + r.get();
                if reach > bounding {
                    bounding = reach;
                }
            }
            let cam = Orbit::framing(demo.embedding.radius_of_gyration().get());
            // The molecule fits when the bounding sphere subtends less than the
            // half-angle: bounding <= distance * sin(fov/2).
            let fits = cam.distance() * sin_half;
            assert!(
                bounding < fits,
                "seed {seed}: the molecule reaches {bounding:.4} against a view that \
                 holds {fits:.4} at the opening distance, so part of it is off screen"
            );
        }
    }

    /// A drag that began on the panel never orbits, and one that began on the
    /// scene keeps the pointer wherever it wanders.
    ///
    /// **The latch is the discriminator.** Deciding per frame — orbit whenever
    /// the pointer is not over the panel — freezes the molecule the instant a
    /// drag crosses onto the periodic table, which is a bug that only shows up
    /// with a mouse in your hand.
    #[test]
    fn a_drag_belongs_to_wherever_it_started() {
        let mut began_on_panel = Orbit::framing(1.0);
        let before = began_on_panel.eye();
        began_on_panel.point(true, true, 40.0, 0.0);
        began_on_panel.point(true, false, 40.0, 0.0);
        assert_eq!(
            bits(began_on_panel.eye()),
            bits(before),
            "a drag that started on the panel turned the molecule once it wandered off"
        );

        let mut began_on_scene = Orbit::framing(1.0);
        began_on_scene.point(true, false, 40.0, 0.0);
        let after_first = began_on_scene.eye();
        began_on_scene.point(true, true, 40.0, 0.0);
        assert_ne!(
            bits(began_on_scene.eye()),
            bits(after_first),
            "a drag that started on the scene stopped turning when it crossed the panel"
        );

        began_on_scene.point(false, false, 0.0, 0.0);
        let released = began_on_scene.eye();
        began_on_scene.point(true, true, 40.0, 0.0);
        assert_eq!(
            bits(began_on_scene.eye()),
            bits(released),
            "the latch survived the button being released"
        );
    }
}
