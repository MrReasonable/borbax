//! The camera: an orientation and a distance, about the origin.
//!
//! **No engine type and no UI type**, which is what lets every test here run as
//! plain arithmetic with no window and no device.
//!
//! # It tumbles freely, and that replaced a camera that did not
//!
//! The first version of this file was a *turntable*: a yaw about a fixed world
//! up, a pitch above the horizon, and a clamp stopping the pitch just short of
//! the poles. The clamp was there for a real reason — a turntable that passes
//! over the top turns the whole image upside down between one frame and the next
//! — but the behaviour it bought was a camera that **stops** when you drag
//! upward, which is not what something you are turning over in your hands does.
//! Found by using it, within a minute, having survived every test.
//!
//! So the fixed up is gone. The camera stores an orientation and turns it about
//! its **own** axes: dragging moves the molecule the way your hand went whichever
//! way up it currently is, rotation never stops in any direction, and the pole is
//! not a special case that has been guarded — it does not exist. The `up` the
//! engine needs comes out of the orientation, so it is perpendicular to the view
//! direction by construction rather than by a clamp.
//!
//! What that costs, stated because it is a real loss: turning by an angle and
//! back is no longer *exactly* the identity, because a rotation composed and
//! renormalised drifts in the last bits where adding and subtracting a stored
//! angle did not. Identical input still gives bit-identical output, which is the
//! claim the purity test actually needs; the round trip is now a tolerance.
//!
//! # What this owns, and what it deliberately does not
//!
//! It owns `orientation -> (eye, up)`. The view matrix, the projection, the
//! perspective divide and the depth range are the engine's, and re-deriving any
//! of them here would be duplicating a dependency in `f32`.
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
//! Every transcendental goes through [`borbax_units::det_math`].

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
pub const NEAR: f64 = 0.1;

/// Radians of rotation per logical point of drag.
const RADIANS_PER_POINT: f64 = 0.01;

/// Multiples of the radius of gyration the camera stands back by.
///
/// **Scaled from a shipped measurement rather than from a bounding radius the
/// viewer would have to compute**, and chosen by measurement rather than by the
/// derivation looking right — the derivation alone was wrong by a factor of two
/// on the first attempt, and the second attempt was still too far away when it
/// was looked at.
///
/// The molecule exactly fills the vertical field at
/// `bounding_radius / sin(FOV_Y / 2)` = `2.613 * bounding_radius`. Measured
/// across 500 universes, `bounding / gyration` runs **1.6491 .. 1.7697**, so
/// touching the edges is **4.3092 .. 4.6243** gyration radii. At 5.4 the
/// molecule spans about 80-86% of the view's height, leaving a margin to turn
/// into without any orientation pushing an atom off screen.
///
/// Vertical is the binding direction because the scene region is wider than it
/// is tall. If that stops being true, this constant is measuring the wrong axis.
const FRAMING: f64 = 5.4;

/// How far in and out of the framing distance the wheel may go.
///
/// **The inner limit is a measured floor, not a preference.** It has to keep the
/// nearest surface in front of the near plane in every universe, which
/// `zooming_all_the_way_in_never_crosses_the_near_plane` checks across 500 of
/// them.
const ZOOM_IN_LIMIT: f64 = 0.42;
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

/// A rotation, as a unit quaternion.
///
/// **Hand-rolled, and the reason is the seam rather than novelty.** This file is
/// in the strictest tier — it may name no engine type — and the engine's own
/// quaternion is `f32` anyway. Sixty lines of arithmetic keeps the camera
/// testable with no app at all and keeps the transcendentals routed through
/// `det_math`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Rotation {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl Rotation {
    /// The rotation that does nothing.
    const IDENTITY: Self = Self {
        w: 1.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// A turn of `angle` radians about a **unit** axis.
    ///
    /// Only ever called with a coordinate axis, so there is no normalisation
    /// here to forget: both callers pass a literal.
    fn about(axis: Vec3, angle: f64) -> Self {
        let half = angle * 0.5;
        let s = det_math::sin(half);
        Self {
            w: det_math::cos(half),
            x: axis[0] * s,
            y: axis[1] * s,
            z: axis[2] * s,
        }
    }

    /// `self`, then `next` **in `self`'s own frame**.
    ///
    /// Right multiplication, which is what makes a drag turn the molecule about
    /// the axes currently facing the screen rather than about the world's. Left
    /// multiplication is a turntable, and a turntable is what this replaced.
    fn then_locally(self, next: Self) -> Self {
        Self {
            w: self.w * next.w - self.x * next.x - self.y * next.y - self.z * next.z,
            x: self.w * next.x + self.x * next.w + self.y * next.z - self.z * next.y,
            y: self.w * next.y - self.x * next.z + self.y * next.w + self.z * next.x,
            z: self.w * next.z + self.x * next.y - self.y * next.x + self.z * next.w,
        }
    }

    /// Scale back to unit length.
    ///
    /// **Called after every turn, and it is not optional.** Composing rotations
    /// accumulates error in the norm, and a quaternion that has drifted off the
    /// unit sphere scales the scene as well as turning it — the molecule would
    /// slowly grow or shrink as it was played with, which is the kind of defect
    /// that takes a long session to notice and no test to hide.
    fn unit(self) -> Self {
        let len = (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if len <= 0.0 {
            return Self::IDENTITY;
        }
        Self {
            w: self.w / len,
            x: self.x / len,
            y: self.y / len,
            z: self.z / len,
        }
    }

    /// Turn `v` by this rotation.
    ///
    /// The `v + 2w(q x v) + 2(q x (q x v))` form, so there is no inverse and no
    /// division. Correct only for a unit quaternion, which [`Self::unit`] keeps
    /// true.
    fn apply(self, v: Vec3) -> Vec3 {
        let q = [self.x, self.y, self.z];
        let t = scale(cross(q, v), 2.0);
        add(add(v, scale(t, self.w)), cross(q, t))
    }
}

/// `a x b`.
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `a + b`, componentwise.
fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// `a * k`, componentwise.
fn scale(a: Vec3, k: f64) -> Vec3 {
    [a[0] * k, a[1] * k, a[2] * k]
}

/// An orientation and a distance, plus whether a drag currently owns the
/// pointer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orbit {
    /// Camera-space to world-space. Kept unit by every mutator.
    orientation: Rotation,
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
    /// radius.
    #[must_use]
    pub fn framing(gyration: f64) -> Self {
        let size = if gyration > MIN_GYRATION {
            gyration
        } else {
            MIN_GYRATION
        };
        let distance = FRAMING * size;
        Self {
            orientation: Rotation::IDENTITY,
            distance,
            framed_at: distance,
            dragging: false,
            was_pressed: false,
        }
    }

    /// Turn by `across` and `up` **radians**, about the camera's own axes.
    ///
    /// **There is no clamp and no special case at any orientation.** Turning
    /// about the axes currently facing the screen has no pole: the up vector is
    /// carried along by the rotation rather than re-derived from a fixed world
    /// axis, so there is no orientation at which it becomes ambiguous.
    pub fn turn(&mut self, across: f64, up: f64) {
        self.orientation = self
            .orientation
            .then_locally(Rotation::about([0.0, 1.0, 0.0], across))
            .then_locally(Rotation::about([1.0, 0.0, 0.0], up))
            .unit();
    }

    /// Turn by a drag of `dx`, `dy` **logical points**.
    ///
    /// Dragging right brings the molecule's left side toward you and dragging up
    /// tips its top toward you — the camera moves against the hand, which is what
    /// makes the object feel grabbed rather than the viewpoint shoved.
    pub fn drag(&mut self, dx: f64, dy: f64) {
        self.turn(-dx * RADIANS_PER_POINT, -dy * RADIANS_PER_POINT);
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

    /// Feed one frame of pointer state in, and turn if this drag is ours.
    ///
    /// **The latch is the whole point of this function**, and it is what makes
    /// the *wide* question below safe to ask. A drag belongs to the scene if the
    /// **press** landed outside the panel; once it has, it keeps the pointer
    /// until release, however far outside the cursor wanders.
    ///
    /// `over_ui` must therefore be the wide question — "is the pointer over the
    /// panel at all". The first shipped version asked the narrow one, "is the
    /// panel *using* the pointer", which is false over a panel's own background,
    /// so **the molecule turned when dragged anywhere on the interface**. Found
    /// by using it. The narrow question is only correct for a per-frame gate,
    /// which this deliberately is not.
    pub fn point(&mut self, pressed: bool, over_ui: bool, dx: f64, dy: f64) {
        if pressed && !self.was_pressed {
            self.dragging = !over_ui;
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
    /// At the identity orientation this is `[0, 0, distance]` — on `+Z`, looking
    /// toward `-Z`.
    #[must_use]
    pub fn eye(&self) -> Vec3 {
        self.orientation.apply([0.0, 0.0, self.distance])
    }

    /// Which way is up for the camera, in scene space.
    ///
    /// **Carried by the orientation rather than fixed at world `+Y`**, which is
    /// what removes the pole: it is perpendicular to the view direction at every
    /// orientation by construction, so the engine never falls back to an
    /// arbitrary axis and the image never flips.
    #[must_use]
    pub fn up(&self) -> Vec3 {
        self.orientation.apply([0.0, 1.0, 0.0])
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
        }
        println!("bounding/gyration: {ratio_lo:.4} .. {ratio_hi:.4}");
        println!(
            "just-fits framing = {:.4} .. {:.4}",
            2.613 * ratio_lo,
            2.613 * ratio_hi
        );
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

    /// `a . b`.
    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    /// A position as raw bits, for the places where the claim really is
    /// *identical* rather than *close*.
    ///
    /// Two cameras built and driven separately are two different computations,
    /// so this is not the vacuous shape where both sides read the same
    /// expression. Every caller asserts finiteness first, because `NaN` bits
    /// compare equal to `NaN` bits.
    fn bits(p: [f64; 3]) -> [u64; 3] {
        [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
    }

    /// A quarter turn moves the eye; a full turn brings it back.
    ///
    /// **Both arms, because either alone is vacuous.** "A full turn returns to
    /// the start" passes on a camera that ignores the angle entirely. "A quarter
    /// turn moves it" passes on degrees-vs-radians, which moves it by the wrong
    /// amount. Together they pin the unit: a quarter turn of a camera at
    /// distance `d` moves the eye by the chord `d * sqrt(2)`, where a quarter
    /// *degree* would move it by 0.0044 `d`.
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

    /// Turning upward never stops, and never degenerates.
    ///
    /// **This is the test for the defect that replaced the previous camera.**
    /// That one clamped its pitch just short of the poles, so dragging upward
    /// stopped dead — correct against the failure it was written for, and wrong
    /// as a thing to use. Ian found it in under a minute.
    ///
    /// Three claims, and the first is the one that fails on a clamp coming back:
    /// a sustained upward drag takes the eye all the way round, so it must spend
    /// time *behind* the molecule (`z < 0`). The others are what the old clamp
    /// existed to protect and which the free rotation gives for nothing: the up
    /// vector stays unit and stays perpendicular to the view, at every
    /// orientation, so there is no pose where the engine has to invent an axis.
    #[test]
    fn turning_never_stops_and_never_degenerates() {
        for (across, up) in [(0.0, 5.0), (0.0, -5.0), (5.0, 0.0), (3.0, 4.0)] {
            let mut cam = Orbit::framing(1.0);
            let mut went_behind = false;
            let mut steps = 0u32;
            for _ in 0..1000 {
                cam.drag(across, up);
                let eye = cam.eye();
                let up_vector = cam.up();
                assert!(
                    (norm(up_vector) - 1.0).abs() < 1e-9,
                    "the up vector stopped being a unit vector: {up_vector:?}"
                );
                assert!(
                    dot(up_vector, eye).abs() < 1e-9 * cam.distance(),
                    "up stopped being perpendicular to the view direction, so the \
                     basis is degenerate at this orientation"
                );
                assert!(
                    (norm(eye) - cam.distance()).abs() < 1e-9,
                    "the eye drifted off the sphere it orbits on, so the rotation \
                     is scaling as well as turning"
                );
                if eye[2] < 0.0 {
                    went_behind = true;
                }
                steps += 1;
            }
            assert_eq!(steps, 1000, "the sweep did not run");
            assert!(
                went_behind,
                "dragging ({across}, {up}) for 1000 steps never took the camera \
                 round the far side, so something is stopping the rotation"
            );
        }
    }

    /// The camera is a pure function of the drags applied to it.
    ///
    /// **The obvious spelling of this was wrong and passed**, which is worse than
    /// failing: it compared one whole drag against two half drags, an assertion
    /// that floating-point addition is *associative*. It held for the constants
    /// in front of me and would have broken later for a reason unreadable from
    /// the test.
    ///
    /// What is claimed instead: the same sequence of drags gives the same camera
    /// **bit for bit**, and turning by an angle and back returns to the start
    /// within a tolerance. The round trip is no longer exact — composing and
    /// renormalising a rotation drifts where adding and subtracting a stored
    /// angle did not — and that is a real cost of the free rotation, stated
    /// rather than hidden.
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

        // **One axis at a time, in matched pairs.** The first draft of this
        // arranged three mixed turns whose angles *summed* to zero and expected
        // the start back. That is only a round trip if rotations commute, and
        // they do not — which is precisely the property that lets this camera
        // tumble freely, so the test was contradicting the feature it sits
        // beside. Turns about a single axis do commute with each other, so
        // these pairs really are inverses.
        let mut there_and_back = Orbit::framing(1.2);
        let start = there_and_back.eye();
        for (across, up) in [(0.5, 0.0), (0.0, 0.25)] {
            for _ in 0..64 {
                there_and_back.turn(across, up);
            }
            for _ in 0..64 {
                there_and_back.turn(-across, -up);
            }
        }
        let end = there_and_back.eye();
        assert!(
            norm([end[0] - start[0], end[1] - start[1], end[2] - start[2]]) < 1e-9,
            "turning back and forth drifted by more than rounding explains"
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

    /// The wheel stops at both ends.
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
    /// one seed, and on a near plane raised without re-checking the zoom limit.
    #[test]
    fn zooming_all_the_way_in_never_crosses_the_near_plane() {
        let mut checked = 0u64;
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
    }

    /// The whole molecule is inside the field of view at the opening camera, in
    /// every universe and at every orientation.
    ///
    /// **Fails on a framing distance tuned to one seed**, which is the shape the
    /// first value of `FRAMING` had. Because the camera orbits at a fixed
    /// distance about a centred embedding, checking the bounding sphere covers
    /// every orientation at once.
    #[test]
    fn every_atom_is_within_the_field_of_view_at_the_opening_camera() {
        let sin_half = borbax_units::det_math::sin(FOV_Y / 2.0);
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
            let fits = cam.distance() * sin_half;
            assert!(
                bounding < fits,
                "seed {seed}: the molecule reaches {bounding:.4} against a view that \
                 holds {fits:.4} at the opening distance, so part of it is off screen"
            );
        }
    }

    /// A drag that began on the panel never turns the molecule, and one that
    /// began on the scene keeps the pointer wherever it wanders.
    ///
    /// **The latch is the discriminator, and it is what lets the caller ask the
    /// *wide* question about egui.** Deciding per frame — turn whenever the
    /// pointer is not over the panel — freezes the molecule the instant a drag
    /// crosses onto the periodic table.
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
