//! The 3D scene: one sphere per atom, a light, and a camera that orbits.
//!
//! **The engine draws; it computes nothing.** Every position and every radius
//! here comes out of `Embedding` unchanged — the only arithmetic applied to a
//! chemistry value is the narrowing to `f32` the GPU requires, and it happens in
//! exactly one function so a reviewer can find it.
//!
//! # Two cameras, and the reason is a measurement rather than tidiness
//!
//! The scene camera has a `viewport` covering only the region below the panel.
//! The egui context is on a *different*, full-window camera, because `bevy_egui`
//! takes egui's screen rectangle from the camera it is attached to: put the
//! context on the viewport-confined camera and the entire periodic table slides
//! down into the scene region under a blank bar. Nothing crashes and no test
//! notices — it was found by rendering it and looking at the pixels.
//!
//! Which camera `bevy_egui` picks is otherwise decided by spawn order, so the
//! automatic choice is switched off and the context placed deliberately.
//!
//! # Query order is never read for anything that matters
//!
//! Bevy reuses entity slots, so after the despawn/respawn a seed change causes,
//! query iteration comes back shuffled. Nothing here depends on it: atoms are
//! spawned with their transform already correct and are never zipped against
//! the coordinate list again. The hazard is designed out rather than guarded.

use bevy::camera::Viewport;
use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll, MouseScrollUnit};
use bevy::prelude::*;
use bevy::window::PrimaryWindow;
use bevy_egui::PrimaryEguiContext;

use crate::orbit::{FOV_Y, NEAR, Orbit};

/// Bevy's `SystemParam` contract hands read-only parameters by value.
///
/// `#[expect]` rather than `#[allow]`, so it fails the build if the systems here
/// ever stop taking one.
#[expect(
    clippy::needless_pass_by_value,
    reason = "a Bevy system takes `Res<T>` and `Single<T>` by value; that is the \
              framework's calling convention, not a missed borrow"
)]
mod systems {
    use super::{
        AccumulatedMouseMotion, AccumulatedMouseScroll, Atom, AtomLook, DrawnAt, FOV_Y,
        LARGEST_VIEWPORT, MouseScrollUnit, NEAR, Orbit, OrbitState, PrimaryWindow, SceneCamera,
        SceneRegion, Viewer, Viewport, WORLD_UP, narrow, pixels, place,
    };
    use bevy::prelude::*;

    /// Build the camera the egui panel is drawn through, the camera the
    /// molecule is drawn through, and the light that makes a sphere look round.
    pub(super) fn spawn_scene(
        mut commands: Commands<'_, '_>,
        mut meshes: ResMut<'_, Assets<Mesh>>,
        mut materials: ResMut<'_, Assets<StandardMaterial>>,
    ) {
        // The 3D camera. `order: 0` so it draws first and the panel's camera
        // composites over it.
        commands.spawn((
            Camera3d::default(),
            Camera {
                order: 0,
                ..default()
            },
            Projection::Perspective(PerspectiveProjection {
                fov: narrow(FOV_Y),
                near: narrow(NEAR),
                ..default()
            }),
            Transform::from_xyz(0.0, 0.0, 1.0).looking_at(Vec3::ZERO, WORLD_UP),
            // **Ambient light as well as the directional one**, so the unlit
            // side of each sphere is dark rather than absent — with a single
            // light and no ambient every atom reads as a crescent, not a ball.
            // It is a component on the camera in this engine version, not a
            // global resource.
            AmbientLight {
                brightness: 220.0,
                ..default()
            },
            SceneCamera,
        ));

        // The panel's camera: full window, no viewport, and it carries the egui
        // context. See this module's header for why it cannot be the one above.
        commands.spawn((
            Camera2d,
            Camera {
                order: 1,
                clear_color: ClearColorConfig::None,
                ..default()
            },
            super::PrimaryEguiContext,
        ));

        // **A light, because a physically-based material with nothing shining on
        // it renders black** — and black on a dark background is a scene that
        // passes every structural test and shows nothing.
        commands.spawn((
            DirectionalLight {
                illuminance: 6_000.0,
                ..default()
            },
            Transform::from_xyz(4.0, 8.0, 6.0).looking_at(Vec3::ZERO, WORLD_UP),
        ));

        // **One mesh and one material, shared by every atom**, and both halves
        // are guarded. A single material handle is what makes "no colour encodes
        // any property" checkable in one assertion rather than by reading — and
        // a single *mesh* handle closes the same door one channel along, because
        // a per-atom mesh keyed on valence (a subdivision count, a different
        // primitive) is a property map reaching the screen through geometry.
        // Measured during review: giving every atom its own mesh handle failed
        // only the material assertion.
        //
        // Atoms are told apart by radius, which is a generated quantity drawn
        // faithfully.
        commands.insert_resource(AtomLook {
            mesh: meshes.add(Mesh::from(Sphere::new(1.0))),
            material: materials.add(StandardMaterial {
                base_color: Color::srgb(0.62, 0.68, 0.78),
                perceptual_roughness: 0.45,
                ..default()
            }),
        });
    }

    /// Replace the atoms when a different molecule has been laid out.
    ///
    /// **Gated on the reload token, not on the resource being "changed".** Two
    /// separate reasons, and the second is the better one:
    ///
    /// - `ResMut<Viewer>` marks the resource changed on every frame the panel
    ///   draws — measured, 40 frames out of 40 — so `resource_changed` would
    ///   tear down and rebuild the whole scene at the refresh rate, silently,
    ///   with exactly the entity churn that shuffles query order.
    /// - Even with that fixed, `resource_changed` would still be wrong: the
    ///   viewer's state also changes on a keystroke that does not parse, on a
    ///   name-box edit and on a cell selection, none of which change the
    ///   molecule. `reloads` is the right *quantity*, not a workaround for a
    ///   framework wart — which matters, because someone who "fixes" the
    ///   `ResMut` would otherwise believe they had removed a hack.
    ///
    /// **`reloads`, not `re_embeds`.** A refused seed lays nothing out, so
    /// `re_embeds` does not move — and gating on it left the previous
    /// universe's molecule on screen underneath the refusal message, with the
    /// whole suite green.
    pub(super) fn respawn_atoms(
        mut commands: Commands<'_, '_>,
        viewer: Res<'_, Viewer>,
        look: Res<'_, AtomLook>,
        mut drawn: ResMut<'_, DrawnAt>,
        mut orbit: ResMut<'_, OrbitState>,
        existing: Query<'_, '_, Entity, With<Atom>>,
    ) {
        let generation = viewer.0.reloads();
        if drawn.0 == generation {
            return;
        }
        drawn.0 = generation;

        for entity in existing.iter() {
            commands.entity(entity).despawn();
        }

        let Some(demo) = viewer.0.demo() else {
            return;
        };
        for (position, radius) in demo.embedding.coords().iter().zip(demo.embedding.radii()) {
            commands.spawn((
                Mesh3d(look.mesh.clone()),
                MeshMaterial3d(look.material.clone()),
                // The sphere mesh has unit radius, so scaling by the atom's own
                // radius is the whole transform. No fudge, no clamp, no
                // per-atom adjustment: the picture is the chemistry.
                Transform::from_translation(place(*position))
                    .with_scale(Vec3::splat(narrow(radius.get()))),
                Atom,
            ));
        }
        orbit.0 = Orbit::framing(demo.embedding.radius_of_gyration().get());
    }

    /// Is the pointer in the region the panel left for the scene?
    ///
    /// **Refuses when it does not know**, in both arms: `SceneRegion` is `None`
    /// until the panel has drawn once, and `cursor_position` is `None` when the
    /// pointer is outside the window. A molecule that does not turn is a smaller
    /// defect than one that turns while a child is reaching for a control, so
    /// the unknown case declines rather than allows.
    fn pointer_is_in_the_scene(region: Option<Rect>, window: Option<&Window>) -> bool {
        let (Some(rect), Some(window)) = (region, window) else {
            return false;
        };
        // Logical points on both sides: `SceneRegion` is documented as logical,
        // and `cursor_position` divides the physical position by the scale
        // factor. Comparing a physical cursor against a logical rectangle would
        // put the boundary in the wrong place on every non-unity display, which
        // is most of them.
        window.cursor_position().is_some_and(|at| rect.contains(at))
    }

    /// Turn the camera with the mouse, and zoom it with the wheel.
    ///
    /// Runs after egui has decided whether it wants the pointer this frame; a
    /// system in `Update` would be reading last frame's answer.
    pub(super) fn drive_camera(
        buttons: Res<'_, ButtonInput<MouseButton>>,
        motion: Res<'_, AccumulatedMouseMotion>,
        wheel: Res<'_, AccumulatedMouseScroll>,
        region: Res<'_, SceneRegion>,
        // **`With<PrimaryWindow>`, not the first row of a bare `&Window`
        // query.** There is one window today, so both pick the same one — but
        // "whichever the query yields first" is precisely the thing this module
        // header says is never read for anything that matters, and a gate that
        // silently followed spawn order would be the same class of defect as
        // the one it replaces.
        window: Query<'_, '_, &Window, With<PrimaryWindow>>,
        mut orbit: ResMut<'_, OrbitState>,
    ) {
        // **Ask the geometry, because `egui` has no answer to give here.** Two
        // versions of this gate asked `egui` whether it wanted the pointer —
        // first the narrow question (`wants_pointer_input`), then the wide one
        // (`wants_any_pointer_input`) — and Ian reported the molecule turning
        // when dragged anywhere on the interface after *both*.
        //
        // The reason is in `app.rs`: the root `Ui` is built on
        // `LayerId::background()` and the panel is shown inside it, so the panel
        // never gets a layer of its own and `egui`'s "is the pointer over me"
        // machinery has nothing to report. Measured against that exact
        // construction with the pointer at `(20, 20)`, squarely on the panel:
        // `is_pointer_over_egui()` false, `egui_wants_pointer_input()` false,
        // `layer_id_at` the background layer. Both answers were false
        // everywhere, which is why swapping one for the other changed nothing.
        //
        // `SceneRegion` is the rectangle the panel actually left, and
        // `aim_camera` already sets the scene camera's viewport from it. Gating
        // on the same value is what stops the gate and the picture disagreeing
        // about where the scene is — the same reason the region is returned by
        // `draw_panel` rather than agreed by convention.
        //
        // The latch still does the rest: `Orbit::point` consults this only on
        // the press edge, so a drag that starts in the scene keeps turning when
        // it crosses onto the panel, which is what a person expects.
        let over_ui = !pointer_is_in_the_scene(region.0, window.iter().next());
        orbit.0.point(
            buttons.pressed(MouseButton::Left),
            over_ui,
            f64::from(motion.delta.x),
            f64::from(motion.delta.y),
        );

        if wheel.delta.y != 0.0 && !over_ui {
            // **A line of wheel and a pixel of trackpad are not the same
            // quantity.** Left unnormalised, a trackpad zooms about a hundred
            // times per notch of a mouse.
            let notches = match wheel.unit {
                MouseScrollUnit::Line => f64::from(wheel.delta.y),
                MouseScrollUnit::Pixel => f64::from(wheel.delta.y) / 50.0,
            };
            orbit.0.zoom(notches);
        }
    }

    /// Point the camera where the orbit says, and confine it to the panel's
    /// leftover space.
    pub(super) fn aim_camera(
        orbit: Res<'_, OrbitState>,
        region: Res<'_, SceneRegion>,
        window: Single<'_, '_, &Window>,
        camera: Single<'_, '_, (&mut Transform, &mut Camera), With<SceneCamera>>,
    ) {
        let (mut transform, mut camera) = camera.into_inner();
        // **The camera's own up, not the world's.** The orbit carries its up
        // vector through the rotation, so it is perpendicular to the view at
        // every orientation — which is what lets the molecule tumble freely
        // without the image ever flipping, and means `looking_at` never has to
        // fall back to an arbitrary axis.
        *transform = Transform::from_translation(place(orbit.0.eye()))
            .looking_at(Vec3::ZERO, place(orbit.0.up()));

        let Some(rect) = region.0 else { return };
        // egui works in logical points and a viewport is in physical pixels.
        // Missing this puts the scene in the wrong quarter of the window on any
        // display with a scale factor, and never on a headless test, where the
        // factor is 1.
        let scale = window.scale_factor();
        let min = rect.min * scale;
        let size = rect.size() * scale;
        // A zero-sized viewport is invalid, so a window dragged too small drops
        // the viewport rather than clamping it to a degenerate one.
        //
        // **Written as `!(a && b)` rather than `a || b`, and that is the fix
        // rather than a style choice.** The `<` form is *false* for NaN, so a
        // non-finite rect skipped the very branch that exists to reject a
        // degenerate viewport and went on to build a zero-sized one — the guard
        // shaped to catch the bad case was the one comparison the bad case
        // passed through. The upper bound is here for the other end the
        // saturating cast leaves open: a very large extent would otherwise
        // become `u32::MAX`.
        let sane = size.x >= 1.0
            && size.y >= 1.0
            && size.x <= LARGEST_VIEWPORT
            && size.y <= LARGEST_VIEWPORT;
        if !sane {
            camera.viewport = None;
            return;
        }
        camera.viewport = Some(Viewport {
            physical_position: UVec2::new(pixels(min.x), pixels(min.y)),
            physical_size: UVec2::new(pixels(size.x), pixels(size.y)),
            ..default()
        });
    }
}

/// Which way is up, before the camera has an opinion.
const WORLD_UP: Vec3 = Vec3::Y;

/// The largest viewport edge, in physical pixels, that is worth believing.
///
/// Nothing plausible reaches it; it exists so a non-finite or absurd layout rect
/// is refused before the saturating cast turns it into `u32::MAX`.
const LARGEST_VIEWPORT: f32 = 65_536.0;

/// Marks the camera the molecule is drawn through.
#[derive(Component, Debug)]
pub struct SceneCamera;

/// Marks one atom of the molecule on screen.
#[derive(Component, Debug)]
pub struct Atom;

/// The orbit camera's parameters, held where a system can reach them.
#[derive(Resource, Default, Debug)]
pub struct OrbitState(pub Orbit);

/// The region the panel left for the scene, in **logical points**.
#[derive(Resource, Default, Debug)]
pub struct SceneRegion(pub Option<Rect>);

/// Which reload the atoms currently on screen were built from.
#[derive(Resource, Default, Debug)]
pub struct DrawnAt(pub u64);

/// The one mesh and the one material every atom shares.
#[derive(Resource, Debug)]
pub struct AtomLook {
    mesh: Handle<Mesh>,
    material: Handle<StandardMaterial>,
}

/// The viewer's state, as the engine holds it.
///
/// Re-exported from `app` so the systems above can name it without `app`
/// needing to know what a scene is.
pub use crate::app::Viewer;

/// The **one** place an `f64` becomes an `f32` in this crate.
///
/// Every coordinate, radius and camera position crosses here and nowhere else,
/// which is what keeps a display-driven narrowing from spreading into a
/// chemistry value.
///
/// **`xtask` counts the spelling — and did not when this sentence was first
/// written.** The claim shipped as prose describing a check nobody had written;
/// a second `const fn shrink2` beside this one passed the whole gate. What was
/// actually enforcing anything was `clippy::as_conversions`, which requires an
/// `#[expect]` per site but does not care how many sites there are.
#[expect(
    clippy::cast_possible_truncation,
    clippy::as_conversions,
    reason = "the GPU takes f32 and the chemistry is f64; this is the single \
              declared boundary between them, and the values crossing it are \
              positions and radii of order 1"
)]
pub const fn narrow(x: f64) -> f32 {
    x as f32
}

/// A chemistry-space point, as the engine wants it.
pub const fn place(v: [f64; 3]) -> Vec3 {
    Vec3::new(narrow(v[0]), narrow(v[1]), narrow(v[2]))
}

/// A logical-point coordinate as a whole physical pixel.
///
/// **A second, different narrowing from [`narrow`], and it is named separately
/// rather than folded in.** A viewport is expressed in whole pixels and there is
/// no floating-point spelling of one, so this cast is forced by the engine's API
/// rather than by the GPU's number format.
///
/// **The guard is for legibility, not for safety, and the first version of this
/// doc claimed otherwise.** It said negatives are floored "instead of wrapping to
/// four billion" — Rust's float-to-integer casts have saturated since 1.45, so a
/// negative already gives 0 and NaN already gives 0 without any help. What the
/// explicit branch buys is a reader who can see the intent and an `#[expect]`
/// that reads as deliberate. The end that is *not* handled here is the top: a
/// very large extent saturates to `u32::MAX`, which the caller rejects before
/// building a viewport out of it.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::as_conversions,
    reason = "a viewport is whole pixels and the engine offers no float form, \
              and the explicit branch states an intent the saturating cast \
              would leave implicit"
)]
fn pixels(x: f32) -> u32 {
    if x > 0.0 { x as u32 } else { 0 }
}

/// Register everything the scene needs.
pub fn plugin(app: &mut App) {
    app.init_resource::<OrbitState>()
        .init_resource::<SceneRegion>()
        .init_resource::<DrawnAt>()
        .add_systems(Startup, systems::spawn_scene)
        // **All three after the panel has drawn, and `respawn_atoms` was in
        // `Update` until a review measured what that costs.** The panel runs in
        // `EguiPrimaryContextPass`, which the bridge schedules *inside*
        // `PostUpdate` — so a seed typed in one frame was not seen by a system
        // in `Update` until the next one, and the caption changed a frame before
        // the molecule did. Sixteen milliseconds, invisible, and exactly the
        // frame-ordering class the panel spent five review rounds on.
        //
        // Chained, because the order among the three is load-bearing:
        // `respawn_atoms` resets the orbit to frame a new molecule, and
        // `aim_camera` must see that reset in the same frame it happens.
        .add_systems(
            PostUpdate,
            (
                systems::respawn_atoms,
                systems::drive_camera,
                systems::aim_camera,
            )
                .chain()
                .after(bevy_egui::EguiPostUpdateSet::ProcessOutput),
        );
}
