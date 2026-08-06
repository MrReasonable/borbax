//! Is the app wired up such that anything can reach the screen?
//!
//! **These exist because 509 tests were green over an empty window.** On
//! 2026-08-03 the Bevy port shipped without spawning a camera. `bevy_egui`
//! paints into a camera's render graph rather than to the window, so with no
//! camera the panel was laid out perfectly and reached no pixel. Every frame
//! test in `tests/panel.rs` stayed green, because they drive `panel::draw`
//! through an `egui` context of their own and never involve the engine at all.
//! It was caught by Ian looking at the window — the fourth defect in this
//! project found that way.
//!
//! **What these tests can and cannot do, stated rather than implied.** They
//! assert on the *structure* of the app: that the resources, systems and
//! entities something has to have in order to draw are present. They do **not**
//! assert that a pixel changed colour — that needs a GPU adapter, and CI's Linux
//! runner has none. So this is the cheap guard for the specific defect that
//! shipped and for its near neighbours, not a claim that the window looks right.
//! A person still has to look; the manual checklist in the Step 3 design memo is
//! the record of that.
//!
//! They drive [`borbax_ui::app::viewer_app`], which is the value `main.rs` runs.
//! That identity is the point: the construction moved out of the binary
//! precisely because a test cannot import one, so any test would otherwise have
//! been asserting about its own copy of the wiring.

use bevy::input::mouse::{MouseButtonInput, MouseMotion, MouseWheel};
use bevy::prelude::*;
use bevy::render::renderer::RenderDevice;
use borbax_ui::app::headless_app;

/// The tests can run where there is no graphics driver.
///
/// **Not a preference — the alternative is that none of them run at all.** Bevy
/// initialises a real GPU adapter as part of building the app, and where it
/// cannot find one it calls `.expect("Unable to find a GPU! …")` inside
/// `bevy_render`. So on a runner with no driver every test in this file aborts
/// in a third-party crate, with a message that reads like a viewer bug.
///
/// **This file now runs on CI, and the correction matters more than the fact.**
/// An earlier version of this paragraph said the branch had never been through
/// CI and that the file had only ever run on a Mac with a Metal device. Both
/// stopped being true on 2026-08-04: the run went green on macOS, Windows and
/// Linux, and `the_viewer_asks_the_platform_for_no_audio_device` below was
/// written *from* a `test (windows-latest)` failure. In a file whose stated
/// purpose is that a claim must be measured rather than believed, a stale
/// provenance note is exactly the sentence the next reader acts on.
///
/// **It widens the gap between this app and the shipped one, and that is worth
/// saying plainly**, because the narrowness of that gap is what makes these
/// tests worth anything. It was exactly one plugin (`winit`); it is now one
/// plugin and the graphics backend. Nothing below asserts about pixels, so
/// nothing below is weakened by there being no device — what is given up is the
/// ability to ever add such an assertion here, which was never available on CI
/// anyway.
#[test]
fn the_headless_app_needs_no_gpu_adapter() {
    let mut app = headless_app();
    app.finish();
    app.cleanup();
    app.update();

    assert!(
        !app.world().contains_resource::<RenderDevice>(),
        "the headless app initialised a real graphics device, so on a runner \
         with no driver `bevy_render` panics with `Unable to find a GPU!` and \
         every test in this file dies inside a third-party crate"
    );
}

/// The app asks the platform for no audio device.
///
/// **This is the guard for a crash no machine here can reproduce.** On
/// 2026-08-04 `test (windows-latest)` died with `STATUS_ACCESS_VIOLATION`
/// (`0xc0000005`) after two of the fourteen tests in this file had passed — a
/// segfault in the test process, with no assertion failing and no backtrace,
/// which reads like a defect in the engine rather than in this crate.
///
/// The one thing that separated the failing run from the passing one was
/// `bevy_audio`. Counted across both logs of the same commit: `No audio device
/// found` appears **twice on Windows and not once on macOS**, while every other
/// marker — the duplicate global-logger error, the `RenderApp` warnings, the
/// missing `winit` window — appears on both, and macOS built **28** apps to
/// Windows' 4 without complaint. Windows Server runs no audio service, so
/// WASAPI enumerates nothing; each app then tears that failure down while the
/// next is building it, on a platform where that initialisation is
/// thread-affine.
///
/// **The viewer plays no sound**, so nothing is given up. The plugin is dropped
/// for the window shell too, not merely this one: if probing a device-less
/// machine can fault, it can fault for a person whose PC has no sound card, and
/// a test-only fix would leave that standing while hiding it from CI. It also
/// *narrows* the gap between the shipped app and this one rather than widening
/// it, which the note on [`borbax_ui::app::headless_app`] cares about.
///
/// This asserts on the plugin rather than on the warning text, because the
/// string belongs to a third-party crate and a quieter release of it would
/// leave this green over the crash it exists to stop. Dropping the dependency
/// outright is issue #28 and would make this test redundant — until then the
/// crate is still compiled and linked, and only this keeps it from running.
#[test]
fn the_viewer_asks_the_platform_for_no_audio_device() {
    let app = headless_app();

    assert!(
        !app.is_plugin_added::<bevy::audio::AudioPlugin>(),
        "the app builds `bevy_audio`, which probes for an output device while \
         the app is being constructed. On a machine with no audio service — \
         every Windows CI runner — that faulted the whole test process with \
         STATUS_ACCESS_VIOLATION and no failing assertion to point at it"
    );
}

/// The app has a camera after startup.
///
/// **The exact defect that shipped.** `bevy_egui` renders through a camera; no
/// camera means no panel, silently.
///
/// Built by running the app's real `Startup` schedule rather than by reading the
/// system list, because "a system is registered" and "the entity exists" are
/// different claims and only the second one is what draws.
#[test]
fn the_app_has_a_camera_to_render_through() {
    let mut app = headless_app();
    app.finish();
    app.cleanup();
    app.update();

    let world = app.world_mut();
    let mut cameras = world.query::<&Camera>();
    let found = cameras.iter(world).count();
    assert!(
        found > 0,
        "the app has no camera after startup, so `bevy_egui` has no render graph \
         to paint into and the window opens empty — with every frame test still \
         green, because they drive `panel::draw` through their own `egui` context \
         and never reach the engine"
    );
}

/// The window is the size and name the constants say.
///
/// The negative arm is the one that carries information: asserting the title is
/// non-empty passes on any string at all, and `WINDOW_TITLE` is pinned elsewhere
/// for what it may say. What this adds is that the *app* uses it — a window
/// built from a literal typed here would pass every existing test.
#[test]
fn the_window_is_built_from_the_shipped_constants() {
    let mut app = headless_app();
    app.finish();
    app.cleanup();
    app.update();

    let world = app.world_mut();
    let mut windows = world.query::<&Window>();
    let window = windows
        .iter(world)
        .next()
        .unwrap_or_else(|| unreachable!("the app always configures a primary window"));

    assert_eq!(
        window.title,
        borbax_ui::WINDOW_TITLE,
        "the window title is not the shipped constant"
    );
    assert_eq!(
        window.resolution.physical_width(),
        u32::from(borbax_ui::WINDOW_SIZE[0]),
        "the window width is not the shipped constant, so the frame test that \
         asserts every cell is reachable is describing a window that does not open"
    );
    assert_eq!(
        window.resolution.physical_height(),
        u32::from(borbax_ui::WINDOW_SIZE[1]),
        "the window height is not the shipped constant"
    );
}

/// The viewer opens on a universe rather than on an empty box.
///
/// Step 1's recorded defect, re-guarded at the app level: `ViewerState::new`
/// rather than `opening` gives a window whose first frame says "type a seed".
#[test]
fn the_app_opens_on_a_universe() {
    let mut app = headless_app();
    app.finish();
    app.cleanup();
    app.update();

    // **Read out of the app, not constructed beside it.** This assertion used
    // to build its own `ViewerState::opening()` and check that — which is a
    // test of `opening()`, not of the app, and would have passed unchanged if
    // the app had inserted `ViewerState::new()`. That is precisely the defect
    // the test is named for, and precisely the shape this file exists to stop.
    let state = &app.world().resource::<borbax_ui::app::Viewer>().0;
    assert!(
        state.status_line().contains("elements"),
        "the app's opening state does not describe a universe, so the window's \
         first frame shows no chemistry: {:?}",
        state.status_line()
    );
}

/// Build the app and run it far enough that startup and the first frames have
/// happened.
fn started() -> App {
    let mut app = headless_app();
    app.finish();
    app.cleanup();
    app.update();
    app.update();
    app
}

/// Put the pointer somewhere and drag it, driving the app's real input path.
///
/// The press and the motion arrive as they do from a device: `ButtonInput` is
/// pressed on this frame's edge, and `AccumulatedMouseMotion` carries the
/// delta. Both are written directly because there is no `winit` here to
/// synthesise them, which is the same reason the cursor is placed by hand.
fn drag_from(app: &mut App, cursor: Vec2, delta: Vec2) -> ([f64; 3], [f64; 3]) {
    let window = place_cursor(app, cursor);

    // Release first, so the next frame carries a genuine press *edge* — the
    // latch only consults the gate there, and a test arriving mid-drag would be
    // testing nothing.
    app.world_mut().write_message(MouseButtonInput {
        button: MouseButton::Left,
        state: bevy::input::ButtonState::Released,
        window,
    });
    app.update();

    let before = eye(app);
    // **Written as messages, not as the resources themselves.** Bevy rebuilds
    // `AccumulatedMouseMotion` from `MouseMotion` in `PreUpdate` every frame, so
    // a test that sets the resource directly has it zeroed before
    // `drive_camera` runs — and then *both* arms of this pair pass, the panel
    // one vacuously. Measured: that is exactly what the first version did.
    app.world_mut().write_message(MouseButtonInput {
        button: MouseButton::Left,
        state: bevy::input::ButtonState::Pressed,
        window,
    });
    app.world_mut().write_message(MouseMotion { delta });
    app.update();
    (before, eye(app))
}

/// Put the pointer at a logical position and return the window it is in.
fn place_cursor(app: &mut App, cursor: Vec2) -> Entity {
    let world = app.world_mut();
    let mut windows = world.query::<(Entity, &mut Window)>();
    let (entity, mut window) = windows
        .iter_mut(world)
        .next()
        .unwrap_or_else(|| unreachable!("the app always configures a primary window"));
    let scale = f64::from(window.scale_factor());
    window.set_physical_cursor_position(Some(bevy::math::DVec2::new(
        f64::from(cursor.x) * scale,
        f64::from(cursor.y) * scale,
    )));
    entity
}

/// Where the eye is, in the orbit's own `f64` — **never narrowed to the
/// engine's `f32`**. The scene narrows once, deliberately, in `place`; a test
/// that narrowed as well would be comparing a quantity the camera never used.
fn eye(app: &App) -> [f64; 3] {
    app.world()
        .resource::<borbax_ui::scene::OrbitState>()
        .0
        .eye()
}

/// How far the eye moved between two frames.
///
/// A magnitude rather than `assert_eq!` on the components: exact float
/// comparison is denied workspace-wide, and `1e-9` is the bar the orbit's own
/// tests already use.
fn moved(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    // Only the five operations §13.1 permits natively: `+ - * /` and `sqrt`,
    // all exactly specified by IEEE-754. Nothing here needs `det_math`.
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The region the panel leaves, as the app itself computed it this frame.
fn scene_region(app: &App) -> Rect {
    app.world()
        .resource::<borbax_ui::scene::SceneRegion>()
        .0
        .unwrap_or_else(|| unreachable!("the panel has run, so it left a region"))
}

/// A drag that starts on the panel does not turn the molecule.
///
/// **Ian found this by using the viewer, twice.** The first report was answered
/// by swapping which question the gate asks `egui` — `wants_pointer_input` for
/// `wants_any_pointer_input` — and the molecule still turned when dragged
/// anywhere on the interface, because *both* answers are `false` here and no
/// test could see it.
///
/// The reason is one line in `app.rs`: the root `Ui` is built on
/// `LayerId::background()` and the panel is shown inside it, so the panel never
/// gets a layer of its own. Measured directly against that construction, with
/// the pointer at `(20, 20)` — squarely on the panel — `is_pointer_over_egui()`
/// is `false`, `egui_wants_pointer_input()` is `false`, and `layer_id_at` gives
/// the background layer. Asking `egui` a different question could not have
/// worked; it had no answer to give.
///
/// So the gate asks the geometry instead, and takes it from
/// [`borbax_ui::scene::SceneRegion`] — the rectangle the panel actually left,
/// which the scene camera's viewport is *already* aimed at. The gate and the
/// picture therefore cannot disagree about where the scene is, which is the
/// same reason that region is returned rather than agreed by convention.
///
/// Asserted on the eye position rather than on a `dragging` flag: the flag is
/// the mechanism, and a gate that set it correctly while still turning the
/// camera would pass.
#[test]
fn a_drag_beginning_on_the_panel_does_not_turn_the_molecule() {
    let mut app = started();
    let region = scene_region(&app);
    // Above the region's top edge is the panel, by construction.
    let on_the_panel = Vec2::new(region.min.x + 20.0, region.min.y / 2.0);

    let (before, after) = drag_from(&mut app, on_the_panel, Vec2::new(40.0, 25.0));

    assert!(
        moved(before, after) < 1e-9,
        "a drag beginning at {on_the_panel:?} — on the panel, above the scene \
         region {region:?} — turned the camera by {}. Every control in the \
         viewer is up there, so the molecule spins whenever a child reaches \
         for one",
        moved(before, after)
    );
}

/// A drag that starts in the scene still turns the molecule.
///
/// **The other half, and it is what stops the fix being "never turn".** A gate
/// that refused every drag would pass the test above and take the feature with
/// it.
#[test]
fn a_drag_beginning_in_the_scene_still_turns_the_molecule() {
    let mut app = started();
    let region = scene_region(&app);
    let in_the_scene = region.center();

    let (before, after) = drag_from(&mut app, in_the_scene, Vec2::new(40.0, 25.0));

    assert!(
        moved(before, after) > 1e-9,
        "a drag beginning at {in_the_scene:?} — the middle of the scene region \
         {region:?} — did not turn the camera, so the viewer cannot be turned \
         at all"
    );
}

/// The wheel over the panel does not zoom the molecule.
///
/// **A separate route through the same gate, and the one without a latch.** The
/// drag consults `over_ui` only on the press edge; the wheel consults it every
/// frame, on the line below. So a change that fixed the drag and left the wheel
/// alone would pass both tests above while the molecule still zoomed as a child
/// scrolled the periodic table.
#[test]
fn the_wheel_over_the_panel_does_not_zoom_the_molecule() {
    let mut app = started();
    let region = scene_region(&app);
    let window = place_cursor(&mut app, Vec2::new(region.min.x + 20.0, region.min.y / 2.0));
    app.update();

    let before = distance(&app);
    app.world_mut().write_message(MouseWheel {
        unit: bevy::input::mouse::MouseScrollUnit::Line,
        x: 0.0,
        y: 3.0,
        window,
        phase: bevy::input::touch::TouchPhase::Moved,
    });
    app.update();

    assert!(
        (before - distance(&app)).abs() < 1e-9,
        "the wheel zoomed the molecule while the pointer was on the panel"
    );
}

/// The wheel over the scene still zooms.
///
/// The arm that stops the fix above being "never zoom".
#[test]
fn the_wheel_over_the_scene_still_zooms_the_molecule() {
    let mut app = started();
    let region = scene_region(&app);
    let window = place_cursor(&mut app, region.center());
    app.update();

    let before = distance(&app);
    app.world_mut().write_message(MouseWheel {
        unit: bevy::input::mouse::MouseScrollUnit::Line,
        x: 0.0,
        y: 3.0,
        window,
        phase: bevy::input::touch::TouchPhase::Moved,
    });
    app.update();

    assert!(
        (before - distance(&app)).abs() > 1e-9,
        "the wheel did not zoom with the pointer in the middle of the scene, so \
         the viewer cannot be zoomed at all"
    );
}

fn distance(app: &App) -> f64 {
    app.world()
        .resource::<borbax_ui::scene::OrbitState>()
        .0
        .distance()
}

/// There are exactly two cameras, and the panel is drawn through the one that
/// is *not* confined to the scene's region.
///
/// **This is the Step 3 form of the empty window, and it is not hypothetical.**
/// `bevy_egui` takes egui's screen rectangle from the camera its context is
/// attached to. Put the context on the 3D camera — whose viewport covers only
/// the area below the panel — and the entire periodic table slides down into
/// that area under a blank bar. Nothing crashes. Rendered and looked at, which
/// is the only way it was found.
///
/// Which camera the engine picks is otherwise spawn-order dependent, so the
/// automatic choice is switched off and this asserts the deliberate one.
#[test]
fn the_panel_is_drawn_through_the_full_window_camera() {
    let mut app = started();
    let world = app.world_mut();

    let mut all = world.query::<(Entity, &Camera)>();
    let cameras: Vec<_> = all.iter(world).map(|(e, c)| (e, c.order)).collect();
    assert_eq!(
        cameras.len(),
        2,
        "expected one camera for the scene and one for the panel, found {cameras:?}"
    );

    let mut scene = world.query_filtered::<Entity, With<borbax_ui::scene::SceneCamera>>();
    let scene_camera = scene
        .iter(world)
        .next()
        .unwrap_or_else(|| unreachable!("the scene camera is spawned at startup"));

    let mut egui_cameras = world.query_filtered::<Entity, With<bevy_egui::PrimaryEguiContext>>();
    let egui_camera = egui_cameras
        .iter(world)
        .next()
        .unwrap_or_else(|| unreachable!("some camera must carry the egui context"));

    assert_ne!(
        scene_camera, egui_camera,
        "the egui context is on the scene camera, so egui's screen rectangle is \
         the scene's viewport — the periodic table will render inside the 3D \
         region with a blank bar above it, and no test but this one can see it"
    );
}

/// The scene has a light.
///
/// **A physically-based material with nothing shining on it renders black**, and
/// black against a dark background is a scene that satisfies every other
/// assertion here and shows nothing at all.
#[test]
fn the_scene_has_a_light() {
    let mut app = started();
    let world = app.world_mut();
    let mut lights = world.query::<&DirectionalLight>();
    assert!(
        lights.iter(world).count() > 0,
        "no light, so every atom renders black"
    );
}

/// One entity per atom, each with a mesh and a material.
#[test]
fn there_is_one_atom_entity_per_atom_of_the_molecule() {
    let expected = borbax_ui::state::ViewerState::opening().demo().map_or_else(
        || unreachable!("the opening universe builds a molecule"),
        |d| d.embedding.len(),
    );

    let mut app = started();
    let world = app.world_mut();
    let mut atoms =
        world.query_filtered::<(&Mesh3d, &MeshMaterial3d<StandardMaterial>), With<borbax_ui::scene::Atom>>();
    assert_eq!(
        atoms.iter(world).count(),
        expected,
        "the scene does not hold one entity per atom"
    );
}

// **`every_atom_is_painted_with_the_same_material` was deleted here, and the
// deletion is part of Step 4 rather than a casualty of it.** It asserted that
// every atom shared one material handle, and its own doc said it held the line
// *until* Step 4 — the guarantee it published, "no colour encodes any
// property", is precisely what this step exists to retire. Amending it to
// admit colour would have left a test whose name no longer described anything.
//
// A deletion with nothing replacing it is a finding, so here is what replaces
// it, and the replacement is strictly stronger. The old guard said *no*
// information reached the colour channel. These say *this* information and no
// other:
//
// - `every_atom_is_painted_the_colour_its_own_element_earns` — the colour is
//   this element's, checked through an oracle that does not use any index;
// - `the_five_atoms_are_five_different_colours` — the palette does not collapse;
// - `every_atom_differs_from_every_other_only_in_its_base_colour` — the rest of
//   the material is still uniform, so roughness and emissive cannot become a
//   second, unguarded property map;
// - `molecule.rs`'s palette tests, which pin what the colour is a function of.
//
// `every_atom_is_drawn_with_the_same_mesh` is **kept, unchanged**. Geometry is
// now the only closed channel, which makes it more load-bearing rather than
// less, and it already queries `With<Atom>` so it excludes sticks — provided
// they carry a different marker, which `no_entity_is_both_an_atom_and_a_bond`
// is what pins.

/// The mesh handle resolves to a mesh that has vertices.
///
/// Fails on a handle to an asset that was never inserted, and on a mesh built
/// with no geometry — neither of which is an error, and both of which draw
/// nothing.
#[test]
fn the_atom_mesh_has_geometry() {
    let mut app = started();
    let world = app.world_mut();
    let mut atoms = world.query_filtered::<&Mesh3d, With<borbax_ui::scene::Atom>>();
    let handle = atoms
        .iter(world)
        .next()
        .map_or_else(|| unreachable!("the scene holds atoms"), |m| m.0.clone());
    let meshes = world.resource::<Assets<Mesh>>();
    let mesh = meshes
        .get(&handle)
        .unwrap_or_else(|| unreachable!("the atom mesh handle resolves"));
    assert!(
        mesh.count_vertices() > 0,
        "the atom mesh has no vertices, so nothing is drawn and nothing errors"
    );
}

/// An idle scene does not rebuild itself.
///
/// **The counter in `ViewerState` cannot see this**, and that is why it is a
/// separate test: rebuilding the atom entities every frame calls no chemistry at
/// all, so `re_embeds` stays flat while the scene is torn down and respawned
/// sixty times a second. It is also exactly the entity churn that shuffles Bevy's
/// query order.
///
/// The natural wrong implementation is `run_if(resource_changed::<Viewer>)` —
/// which looks right and is not, because the panel takes the viewer mutably
/// every frame and that marks it changed every frame.
#[test]
fn an_idle_scene_keeps_the_atoms_it_already_has() {
    let mut app = started();

    let before: Vec<Entity> = {
        let world = app.world_mut();
        let mut atoms = world.query_filtered::<Entity, With<borbax_ui::scene::Atom>>();
        atoms.iter(world).collect()
    };
    assert!(!before.is_empty(), "no atoms, so this asserted nothing");

    for _ in 0..30 {
        app.update();
    }

    let after: Vec<Entity> = {
        let world = app.world_mut();
        let mut atoms = world.query_filtered::<Entity, With<borbax_ui::scene::Atom>>();
        atoms.iter(world).collect()
    };
    // **Compared as sets.** `assert_eq!` on two `Vec`s compares order, and the
    // scene's own header records that query iteration is not spawn order after
    // entity churn — so the sequence form would be a test that fails
    // intermittently with a message naming the wrong cause.
    let before: std::collections::BTreeSet<Entity> = before.into_iter().collect();
    let after: std::collections::BTreeSet<Entity> = after.into_iter().collect();
    assert_eq!(
        before, after,
        "the atom entities were replaced while nothing changed, so the scene is \
         being rebuilt every frame"
    );
}

/// A refused seed takes the molecule off the screen.
///
/// **A real defect, found by review and green over 532 tests.** The scene used
/// to be gated on the layout counter, which does not move when a seed is
/// refused — so typing nonsense left the *previous* universe's molecule sitting
/// under the refusal message. The Step 3 design memo's manual checklist says in
/// as many words that this must not happen: it is `Outcome`'s whole argument —
/// a shape attributed to a universe that did not produce it — carried one layer
/// down into the entity cache, where the enum could not reach it.
#[test]
fn a_refused_seed_takes_the_molecule_off_the_screen() {
    let mut app = started();

    let before = {
        let world = app.world_mut();
        let mut atoms = world.query_filtered::<Entity, With<borbax_ui::scene::Atom>>();
        atoms.iter(world).count()
    };
    assert!(
        before > 0,
        "no molecule to begin with, so this asserts nothing"
    );

    {
        let mut viewer = app.world_mut().resource_mut::<borbax_ui::app::Viewer>();
        viewer.0.seed_text_mut().clear();
        viewer.0.seed_text_mut().push_str("not a seed");
        viewer.0.commit_typed_seed();
    }
    app.update();

    let world = app.world_mut();
    let mut atoms = world.query_filtered::<Entity, With<borbax_ui::scene::Atom>>();
    let after = atoms.iter(world).count();
    assert_eq!(
        after, 0,
        "the refusal left {after} atoms on screen — a molecule from a universe \
         that is no longer loaded, beside a message saying the seed was refused"
    );
}

/// Every atom is drawn exactly where the solver put it, at the size it said.
///
/// **This is the test that says the picture is the chemistry, and its absence
/// was the review's most damning measurement.** With no such test, atoms were
/// moved to `[x*3, y, z+9]` and scaled into ellipsoids `(1.7r, 0.4r, 1.0r)` —
/// fabricated positions, squashed spheres — and all 81 tests passed, against a
/// comment in the scene reading "No fudge, no clamp, no per-atom adjustment".
/// The guard set was aimed at the *colour* channel, which carries nothing;
/// position and radius, which carry everything, were unverified.
///
/// **Compared as a set, not a sequence.** Bevy reuses entity slots, so query
/// order after a despawn/respawn cycle is not spawn order — zipping against the
/// coordinate list would be a test that passes or fails depending on churn.
///
/// **Bit-exact, deliberately.** Both sides run the same `f64 -> f32` narrowing
/// on the same values in the same order, so any difference is a transform this
/// crate applied rather than a rounding question.
#[test]
fn every_atom_is_drawn_where_the_embedding_put_it() {
    use std::collections::BTreeSet;

    let want: BTreeSet<[u32; 6]> = {
        let state = borbax_ui::state::ViewerState::opening();
        let demo = state
            .demo()
            .unwrap_or_else(|| unreachable!("the opening universe builds a molecule"));
        // **The window opens on the sticks view, so the expected radius is the
        // drawn one — and naming it here is the point.** The obvious repair
        // when this test started failing was to delete it, which QA flagged as
        // the most expensive mistake available in this step: it is the guard
        // that exists because atoms were once moved to `[x*3, y, z+9]` and
        // squashed into ellipsoids with all 81 tests green. A display scale is
        // exactly the licence under which a per-atom fudge would re-enter, so
        // the scale is written into the expectation rather than the expectation
        // being dropped. `every_atom_is_shrunk_by_the_same_amount` in
        // `molecule.rs` is what stops the scale becoming per-atom;
        // `the_solid_view_draws_every_atom_at_its_true_size_and_no_sticks` is
        // unscaled picture honest.
        assert!(
            state.view().draws_sticks(),
            "the window no longer opens on the sticks view, so this expectation \
             names the wrong radius"
        );
        demo.embedding
            .coords()
            .iter()
            .zip(&demo.stick_radii)
            .map(|(p, r)| {
                let scale = borbax_ui::scene::narrow(*r);
                let at = borbax_ui::scene::place(*p);
                [
                    at.x.to_bits(),
                    at.y.to_bits(),
                    at.z.to_bits(),
                    scale.to_bits(),
                    scale.to_bits(),
                    scale.to_bits(),
                ]
            })
            .collect()
    };

    let mut app = started();
    let world = app.world_mut();
    let mut atoms = world.query_filtered::<&Transform, With<borbax_ui::scene::Atom>>();
    let got: BTreeSet<[u32; 6]> = atoms
        .iter(world)
        .map(|t| {
            [
                t.translation.x.to_bits(),
                t.translation.y.to_bits(),
                t.translation.z.to_bits(),
                t.scale.x.to_bits(),
                t.scale.y.to_bits(),
                t.scale.z.to_bits(),
            ]
        })
        .collect();

    assert_eq!(
        got.len(),
        want.len(),
        "the scene holds {} distinct atom transforms against {} atoms in the \
         embedding",
        got.len(),
        want.len()
    );
    assert_eq!(
        got, want,
        "an atom is not where the solver put it, or is not the size it said — \
         the picture has stopped being the chemistry"
    );
}

/// The sphere the atoms are drawn with has unit radius.
///
/// **The transform test above cannot see this**, which is why it is separate: a
/// mesh of radius 0.5 with every scale correct draws every atom at half size and
/// leaves both sets identical. Scale is only a radius if the thing being scaled
/// is a unit sphere.
#[test]
fn the_atom_mesh_is_a_unit_sphere() {
    let mut app = started();
    let world = app.world_mut();
    let mut atoms = world.query_filtered::<&Mesh3d, With<borbax_ui::scene::Atom>>();
    let handle = atoms
        .iter(world)
        .next()
        .map_or_else(|| unreachable!("the scene holds atoms"), |m| m.0.clone());
    let meshes = world.resource::<Assets<Mesh>>();
    let mesh = meshes
        .get(&handle)
        .unwrap_or_else(|| unreachable!("the atom mesh handle resolves"));
    let positions = mesh
        .attribute(Mesh::ATTRIBUTE_POSITION)
        .and_then(|values| values.as_float3())
        .unwrap_or_else(|| unreachable!("a sphere mesh has float3 positions"));
    assert!(!positions.is_empty(), "the atom mesh has no vertices");
    let mut reach = 0.0f32;
    for p in positions {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        if r > reach {
            reach = r;
        }
    }
    assert!(
        (reach - 1.0).abs() < 1e-3,
        "the atom mesh reaches {reach} from its centre, not 1.0 — so scaling it \
         by an atom's radius does not draw that radius, and the transform test \
         cannot see it because every scale would still match"
    );
}

/// Every atom shares one mesh handle, as well as one material.
///
/// The colour channel is guarded; this is the same door one channel along. A
/// per-atom mesh keyed on a generated property — subdivision count, a different
/// primitive — is a property map arriving through geometry, and the material
/// assertion cannot see it. Measured during review: giving every atom its own
/// mesh handle failed only the material test.
#[test]
fn every_atom_is_drawn_with_the_same_mesh() {
    let mut app = started();
    let world = app.world_mut();
    let mut atoms = world.query_filtered::<&Mesh3d, With<borbax_ui::scene::Atom>>();
    let handles: Vec<_> = atoms.iter(world).map(|m| m.0.id()).collect();
    assert!(!handles.is_empty(), "no atoms, so this asserted nothing");
    let first = handles
        .first()
        .copied()
        .unwrap_or_else(|| unreachable!("checked non-empty above"));
    assert!(
        handles.iter().all(|h| *h == first),
        "the atoms hold {} distinct meshes — geometry that varies per atom is a \
         property map, and Step 4 owns those with their own guard",
        handles
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    );
}

/// One stick entity per bond of the molecule.
///
/// **Weak on its own and labelled so.** It cannot separate a correct bond
/// lookup from "every pair containing the highest index" — on the demo
/// molecule those agree in every universe, which is why the discriminating
/// bond tests are the fixtures in `molecule.rs` rather than anything here. It
/// stays because a count failure names the cause where a set failure names a
/// symptom.
#[test]
fn there_is_one_stick_entity_per_bond() {
    let expected = borbax_ui::state::ViewerState::opening().demo().map_or_else(
        || unreachable!("the opening universe builds a molecule"),
        |d| d.bonds.len(),
    );
    assert!(
        expected > 0,
        "the molecule has no bonds, so this asserts nothing"
    );

    let mut app = started();
    let world = app.world_mut();
    let mut sticks = world.query_filtered::<Entity, With<borbax_ui::scene::Bond>>();
    assert_eq!(
        sticks.iter(world).count(),
        expected,
        "the scene does not hold one stick per bond"
    );
}

/// Every stick runs between the two atoms it joins.
///
/// **This is the bond channel's answer to "the picture is the chemistry".** The
/// equivalent fabrication for a stick is the right count in the right places
/// pointing the wrong way, or stretched — none of which any count test can see,
/// and all of which are exactly what the position guard was added for on the
/// atom channel.
///
/// **Compared as a set, with each pair internally ordered**, for two separate
/// reasons: Bevy query order is not spawn order after entity churn, and a bond
/// has two endpoints, so `(a, b)` and `(b, a)` are the same stick.
///
/// **Not bit-exact, unlike the atom transform test, and the difference is
/// real.** A cylinder's placement involves a rotation, so the endpoints come
/// back through `f32` trigonometry rather than through the same narrowing
/// applied to the same values. The tolerance below is the measured
/// reconstruction residue, not a number chosen to make this pass.
#[test]
fn every_stick_spans_the_two_atoms_it_bonds() {
    /// A bound on the endpoint error, comfortably above the measured residue.
    ///
    /// **It is a bound, not the measured residue, and an earlier version of
    /// this doc claimed the second.** Measured: **1.23e-7** on the opening
    /// universe and **2.98e-7** worst over 500 — so this sits ~340x above what
    /// the reconstruction actually produces. The defect it must catch moves an
    /// endpoint by order 0.05 span, which is 500x above this, so the bar sits
    /// between two populations three orders apart in each direction.
    ///
    /// **It is not safe in the other direction, and that is worth knowing
    /// before it fires.** `glam` documents `from_rotation_arc` as accurate only
    /// to about 1e-3 for near-singular inputs, and a review measured an
    /// endpoint residue of **9.3e-4** for a bond within 6e-4 radians of the Y
    /// axis — 9x this tolerance, on a completely correct implementation. That
    /// is unreachable today: the largest `|dot(bond direction, Y)|` over 500
    /// universes is **0.9543**. If a later molecule reaches the axis, this
    /// fires with a message blaming a stick that is none of the things it
    /// names, and the fix then is to compare length and midpoint tightly (they
    /// involve no rotation) and the direction against glam's documented 1e-3.
    const TOLERANCE: f32 = 1e-4;

    let want: Vec<(Vec3, Vec3)> = {
        let state = borbax_ui::state::ViewerState::opening();
        let demo = state
            .demo()
            .unwrap_or_else(|| unreachable!("the opening universe builds a molecule"));
        demo.bonds
            .iter()
            .map(|[a, b]| {
                let coords = demo.embedding.coords();
                let from = borbax_ui::scene::place(
                    *coords
                        .get(usize::from(*a))
                        .unwrap_or_else(|| unreachable!("a bond names an atom")),
                );
                let to = borbax_ui::scene::place(
                    *coords
                        .get(usize::from(*b))
                        .unwrap_or_else(|| unreachable!("a bond names an atom")),
                );
                (from, to)
            })
            .collect()
    };

    let mut app = started();
    let world = app.world_mut();
    let mut sticks = world.query_filtered::<&Transform, With<borbax_ui::scene::Bond>>();
    let got: Vec<(Vec3, Vec3)> = sticks
        .iter(world)
        .map(|t| {
            // The mesh is a unit cylinder along Y, so its ends sit at
            // ±0.5 in mesh space and the scale carries the length.
            let half = t.rotation * (Vec3::Y * (t.scale.y / 2.0));
            (t.translation - half, t.translation + half)
        })
        .collect();

    assert_eq!(
        got.len(),
        want.len(),
        "the scene holds {} sticks against {} bonds",
        got.len(),
        want.len()
    );

    let mut worst = 0.0_f32;
    let mut matched = 0_usize;
    for (from, to) in &want {
        // Each wanted bond must be matched by some drawn stick, in either
        // endpoint order.
        let mut best = f32::INFINITY;
        for (a, b) in &got {
            // Explicit comparisons: `f32::max` is disallowed under §13.1 for
            // disagreeing about NaN between platforms, and the deny reaches
            // inside integration tests too.
            let (fa, fb) = (from.distance(*a), to.distance(*b));
            let forward = if fa > fb { fa } else { fb };
            let (ra, rb) = (from.distance(*b), to.distance(*a));
            let reversed = if ra > rb { ra } else { rb };
            let error = if forward < reversed {
                forward
            } else {
                reversed
            };
            if error < best {
                best = error;
            }
        }
        assert!(
            best < TOLERANCE,
            "a bond from {from:?} to {to:?} is drawn by no stick — the nearest \
             is out by {best}, so a stick is mispointed, stretched or moved"
        );
        if best > worst {
            worst = best;
        }
        matched += 1;
    }
    // **A count, not a finiteness check, and the difference is the whole
    // point.** This assertion read `worst.is_finite()` — with `worst`
    // initialised to `0.0`, which is finite, so it passed on precisely the
    // empty corpus it exists to catch. Four independent reviewers found it. A
    // self-check that cannot fire is the shape CLAUDE.md's "probe the guard's
    // bookkeeping" section names, and it was written in the same commit as
    // three counters that get this right.
    assert_eq!(
        matched,
        want.len(),
        "{matched} of {} bonds were matched, so this test did not describe the \
         whole molecule",
        want.len()
    );
    assert!(
        matched > 0,
        "the molecule has no bonds, so this test asserted nothing at all"
    );
}

/// No entity is both an atom and a bond.
///
/// **Not independently discriminating, and saying so is the point.** Spawning
/// sticks with the `Atom` marker was measured to fail **eight** other tests. Its value is
/// that those three would report the symptom in the wrong vocabulary — a wrong
/// atom count, a wrong transform set — while this names the cause. It is also
/// what keeps `every_atom_is_drawn_with_the_same_mesh` describing atoms now
/// that a second entity kind exists.
#[test]
fn no_entity_is_both_an_atom_and_a_bond() {
    let mut app = started();
    let world = app.world_mut();

    let atoms = world
        .query_filtered::<Entity, With<borbax_ui::scene::Atom>>()
        .iter(world)
        .count();
    let bonds = world
        .query_filtered::<Entity, With<borbax_ui::scene::Bond>>()
        .iter(world)
        .count();
    // The positive arm: without it this passes on a scene holding nothing.
    assert!(
        atoms > 0 && bonds > 0,
        "the scene holds {atoms} atoms and {bonds} bonds, so the check below asserts nothing"
    );

    let both = world
        .query_filtered::<Entity, (With<borbax_ui::scene::Atom>, With<borbax_ui::scene::Bond>)>()
        .iter(world)
        .count();
    assert_eq!(
        both, 0,
        "{both} entities carry both markers, so every test querying \
         `With<Atom>` has silently started counting sticks too"
    );
}

/// Every stick shares one material.
///
/// A bond has no generated property of its own to encode — every bond this
/// crate builds is `SINGLE` — so a per-bond colour would be a property map with
/// no property, which is the shape the atom guard used to forbid and which is
/// still forbidden here.
#[test]
fn every_stick_is_painted_with_the_same_material() {
    let mut app = started();
    let world = app.world_mut();
    let mut sticks =
        world.query_filtered::<&MeshMaterial3d<StandardMaterial>, With<borbax_ui::scene::Bond>>();
    let handles: Vec<_> = sticks.iter(world).map(|m| m.0.id()).collect();
    assert!(!handles.is_empty(), "no sticks, so this asserted nothing");
    let first = handles
        .first()
        .copied()
        .unwrap_or_else(|| unreachable!("checked non-empty above"));
    assert!(
        handles.iter().all(|h| *h == first),
        "the sticks hold {} distinct materials — a per-bond colour is a property \
         map, and every bond here is the same order",
        handles
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    );
}

/// The solid view draws no sticks and shrinks nothing.
///
/// **Both halves in one test, because they are one claim**: the solid view is
/// the picture Step 3 shipped, unchanged. The radius comparison is bit-exact
/// for that reason — a display scale that was "1.0 but not quite" would be a
/// silent change to the view this exists to preserve.
#[test]
fn the_solid_view_draws_every_atom_at_its_true_size_and_no_sticks() {
    use std::collections::BTreeSet;

    let want: BTreeSet<u32> = {
        let state = borbax_ui::state::ViewerState::opening();
        let demo = state
            .demo()
            .unwrap_or_else(|| unreachable!("the opening universe builds a molecule"));
        demo.embedding
            .radii()
            .iter()
            .map(|r| borbax_ui::scene::narrow(r.get()).to_bits())
            .collect()
    };

    let mut app = started();
    {
        let mut viewer = app.world_mut().resource_mut::<borbax_ui::app::Viewer>();
        assert!(
            viewer.0.view().draws_sticks(),
            "the window no longer opens on the sticks view, so this flips the wrong way"
        );
        viewer.0.flip_view();
    }
    app.update();

    let world = app.world_mut();
    let sticks = world
        .query_filtered::<Entity, With<borbax_ui::scene::Bond>>()
        .iter(world)
        .count();
    assert_eq!(
        sticks, 0,
        "the solid view left {sticks} sticks on screen — they are buried inside \
         the overlapping atoms, so nobody would see them and nothing else would fail"
    );

    let mut atoms = world.query_filtered::<&Transform, With<borbax_ui::scene::Atom>>();
    let got: BTreeSet<u32> = atoms.iter(world).map(|t| t.scale.x.to_bits()).collect();
    assert_eq!(
        got, want,
        "the solid view is not drawing atoms at their true size, so it is no \
         longer the picture Step 3 shipped"
    );
}

/// Flipping the view redraws the molecule.
///
/// **The defect this catches shipped once already, one token along.** The scene
/// rebuilds on a change token; flipping the view generates no universe, so
/// `reloads` does not move. A token of only the reload count leaves the
/// previous picture on screen after the toggle — atoms at the wrong size and
/// sticks that should not be there — with every other test in this file green,
/// because they all inspect the opening state.
#[test]
fn flipping_the_view_redraws_the_molecule() {
    let mut app = started();

    let before: Vec<u32> = {
        let world = app.world_mut();
        let mut atoms = world.query_filtered::<&Transform, With<borbax_ui::scene::Atom>>();
        atoms.iter(world).map(|t| t.scale.x.to_bits()).collect()
    };
    assert!(!before.is_empty(), "no atoms, so this asserts nothing");

    {
        let mut viewer = app.world_mut().resource_mut::<borbax_ui::app::Viewer>();
        viewer.0.flip_view();
    }
    app.update();

    let after: Vec<u32> = {
        let world = app.world_mut();
        let mut atoms = world.query_filtered::<&Transform, With<borbax_ui::scene::Atom>>();
        atoms.iter(world).map(|t| t.scale.x.to_bits()).collect()
    };
    assert_eq!(
        after.len(),
        before.len(),
        "the flip changed how many atoms are on screen"
    );
    assert_ne!(
        before, after,
        "the atoms are the same size after flipping the view, so the scene did \
         not notice the flip and is still drawing the previous picture"
    );
}

/// A refused seed takes the bonds off too.
///
/// **The exact Step 3 defect, one entity kind along.** Despawning only
/// `With<Atom>` leaves four sticks hanging in space under the refusal message;
/// `a_refused_seed_takes_the_molecule_off_the_screen` counts atoms and stays
/// green over it.
#[test]
fn a_refused_seed_takes_the_bonds_off_the_screen_too() {
    let mut app = started();

    let before = {
        let world = app.world_mut();
        let mut sticks = world.query_filtered::<Entity, With<borbax_ui::scene::Bond>>();
        sticks.iter(world).count()
    };
    assert!(
        before > 0,
        "no sticks to begin with, so this asserts nothing"
    );

    {
        let mut viewer = app.world_mut().resource_mut::<borbax_ui::app::Viewer>();
        viewer.0.seed_text_mut().clear();
        viewer.0.seed_text_mut().push_str("not a seed");
        viewer.0.commit_typed_seed();
    }
    app.update();

    let world = app.world_mut();
    let mut sticks = world.query_filtered::<Entity, With<borbax_ui::scene::Bond>>();
    let after = sticks.iter(world).count();
    assert_eq!(
        after, 0,
        "the refusal left {after} sticks on screen — bonds of a molecule from a \
         universe that is no longer loaded"
    );
}

/// Every atom is painted the colour its own element earns.
///
/// **The oracle is the atom's drawn size, not its position in a query.** Bevy
/// query order is not spawn order, and `Demo::elements` is not canonical order,
/// so the only thing on screen that identifies an atom is how big it is —
/// which is a generated quantity, distinct across the five in every universe
/// measured. The uniqueness is asserted rather than trusted, so a corpus where
/// the oracle stopped working reports that instead of matching the wrong atom.
///
/// This is the *scene's* half of the claim: that it attaches colour `i` to the
/// atom drawn at radius `i`. That colour `i` is the right colour for the
/// element at canonical index `i` is `molecule.rs`'s
/// `every_atom_is_coloured_from_its_own_element`, which uses a wholly
/// independent oracle.
#[test]
fn every_atom_is_painted_the_colour_its_own_element_earns() {
    use std::collections::BTreeMap;

    let want: BTreeMap<u32, Color> = {
        let state = borbax_ui::state::ViewerState::opening();
        let demo = state
            .demo()
            .unwrap_or_else(|| unreachable!("the opening universe builds a molecule"));
        let map: BTreeMap<u32, Color> = demo
            .stick_radii
            .iter()
            .zip(&demo.colours)
            .map(|(r, c)| {
                (
                    borbax_ui::scene::narrow(*r).to_bits(),
                    borbax_ui::scene::srgb(*c),
                )
            })
            .collect();
        assert_eq!(
            map.len(),
            demo.stick_radii.len(),
            "two atoms of this molecule are drawn the same size, so the size \
             cannot identify which atom is which and this test's oracle has \
             stopped working"
        );
        map
    };

    let mut app = started();
    let world = app.world_mut();
    let mut atoms = world
        .query_filtered::<(&Transform, &MeshMaterial3d<StandardMaterial>), With<borbax_ui::scene::Atom>>();
    let painted: Vec<(u32, AssetId<StandardMaterial>)> = atoms
        .iter(world)
        .map(|(t, m)| (t.scale.x.to_bits(), m.0.id()))
        .collect();
    assert_eq!(
        painted.len(),
        want.len(),
        "the scene holds {} atoms against {} in the molecule",
        painted.len(),
        want.len()
    );

    let materials = world.resource::<Assets<StandardMaterial>>();
    let mut checked = 0_u32;
    for (size, handle) in &painted {
        let expected = want
            .get(size)
            .unwrap_or_else(|| unreachable!("every drawn size belongs to an atom"));
        let material = materials
            .get(*handle)
            .unwrap_or_else(|| unreachable!("an atom's material handle resolves"));
        assert_eq!(
            material.base_color, *expected,
            "the atom drawn at size {size} is painted a colour a different atom \
             earns — the scene has permuted the colours against the coordinates"
        );
        checked += 1;
    }
    assert_eq!(
        usize::try_from(checked).unwrap_or(usize::MAX),
        want.len(),
        "not every atom was checked"
    );
}

/// The five atoms are five different colours.
///
/// **The measurement behind this is why the palette is keyed the way it is.** A
/// palette on valence alone gives three of the four leaves one colour in every
/// universe — the valence multiset is `[1, 2, 2, 2, 4]` in 500 of 500 — and a
/// smooth ramp over any single channel collapses the leaves, which sit at
/// consecutive groups. Both pass every other test in this file.
#[test]
fn the_five_atoms_are_five_different_colours() {
    use std::collections::BTreeSet;

    let mut app = started();
    let world = app.world_mut();
    let mut atoms =
        world.query_filtered::<&MeshMaterial3d<StandardMaterial>, With<borbax_ui::scene::Atom>>();
    let handles: Vec<_> = atoms.iter(world).map(|m| m.0.id()).collect();
    assert!(!handles.is_empty(), "no atoms, so this asserts nothing");

    let materials = world.resource::<Assets<StandardMaterial>>();
    let colours: BTreeSet<[u32; 4]> = handles
        .iter()
        .map(|h| {
            let c = materials
                .get(*h)
                .unwrap_or_else(|| unreachable!("an atom's material handle resolves"))
                .base_color
                .to_srgba();
            [
                c.red.to_bits(),
                c.green.to_bits(),
                c.blue.to_bits(),
                c.alpha.to_bits(),
            ]
        })
        .collect();

    assert_eq!(
        colours.len(),
        handles.len(),
        "the {} atoms carry only {} distinct colours, so two of them are \
         indistinguishable on screen",
        handles.len(),
        colours.len()
    );
}

/// Atoms differ in their colour and in nothing else about their material.
///
/// **This is the part of the retired single-material guard the colour tests do
/// not recover.** That guard covered the whole material; Step 4 legitimises
/// exactly one field of it. Without this, "the palette is allowed to vary"
/// becomes the sentence under which roughness or emissive quietly becomes a
/// second, unguarded property map.
#[test]
fn every_atom_differs_from_every_other_only_in_its_base_colour() {
    use std::collections::BTreeSet;

    let mut app = started();
    let world = app.world_mut();
    let mut atoms =
        world.query_filtered::<&MeshMaterial3d<StandardMaterial>, With<borbax_ui::scene::Atom>>();
    let handles: Vec<_> = atoms.iter(world).map(|m| m.0.id()).collect();
    assert!(
        handles.len() > 1,
        "fewer than two atoms, so this asserts nothing"
    );

    let materials = world.resource::<Assets<StandardMaterial>>();
    // Clone each material and overwrite the one field that is allowed to vary.
    // **A clone-and-compare rather than an enumerated field list**, so a field
    // nobody named is covered by construction — an enumeration would go silent
    // at the next engine bump, which is the `[lints] workspace = true` shape.
    let flattened: BTreeSet<String> = handles
        .iter()
        .map(|h| {
            let mut m = materials
                .get(*h)
                .unwrap_or_else(|| unreachable!("an atom's material handle resolves"))
                .clone();
            m.base_color = Color::WHITE;
            format!("{m:?}")
        })
        .collect();

    assert_eq!(
        flattened.len(),
        1,
        "the atoms' materials differ in {} ways once colour is set aside, so \
         something other than `base_color` is carrying a per-atom property",
        flattened.len()
    );
}

/// The solid view paints the same colours the sticks view does.
///
/// **Written because Ian reported the solid view as having no colours**, and
/// because no test could have told him otherwise: every colour test in this
/// file runs on the opening state, which is the *sticks* view. The colour path
/// is shared between the two, so this should pass — and "should" is exactly the
/// word that means a test is missing.
#[test]
fn the_solid_view_paints_the_same_colours_as_the_sticks_view() {
    use std::collections::BTreeSet;

    fn colours_now(app: &mut App) -> BTreeSet<[u32; 4]> {
        let world = app.world_mut();
        let mut atoms = world
            .query_filtered::<&MeshMaterial3d<StandardMaterial>, With<borbax_ui::scene::Atom>>();
        let handles: Vec<_> = atoms.iter(world).map(|m| m.0.id()).collect();
        let materials = world.resource::<Assets<StandardMaterial>>();
        handles
            .iter()
            .map(|h| {
                let c = materials
                    .get(*h)
                    .unwrap_or_else(|| unreachable!("an atom's material handle resolves"))
                    .base_color
                    .to_srgba();
                [
                    c.red.to_bits(),
                    c.green.to_bits(),
                    c.blue.to_bits(),
                    c.alpha.to_bits(),
                ]
            })
            .collect()
    }

    let mut app = started();
    let sticks = colours_now(&mut app);
    assert_eq!(
        sticks.len(),
        5,
        "the sticks view does not show five colours"
    );

    {
        let mut viewer = app.world_mut().resource_mut::<borbax_ui::app::Viewer>();
        viewer.0.flip_view();
        assert!(
            !viewer.0.view().draws_sticks(),
            "the flip did not reach the solid view"
        );
    }
    app.update();

    let solid = colours_now(&mut app);
    assert_eq!(
        solid.len(),
        5,
        "the solid view shows {} distinct atom colours against five in the \
         sticks view — the colours are not reaching this view",
        solid.len()
    );
    assert_eq!(
        solid, sticks,
        "the solid view paints different colours from the sticks view, though \
         both read the same `Demo::colours`"
    );
}
