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
/// This branch has never been through CI — the workflow triggers on `main`
/// pushes and pull requests, and nothing here has been either — so the whole
/// file has only ever run on a Mac with a working Metal device.
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

/// Every atom shares one material handle.
///
/// **This is what makes "no colour encodes any property" checkable rather than a
/// matter of reading the code.** Step 4 owns palettes, and it lands with the
/// guard that says a palette may not be built out of an element's mass, valence
/// or affinity. Until then a per-atom colour is not merely out of scope — one
/// keyed on an index looks exactly like a property map and encodes nothing.
#[test]
fn every_atom_is_painted_with_the_same_material() {
    let mut app = started();
    let world = app.world_mut();
    let mut atoms =
        world.query_filtered::<&MeshMaterial3d<StandardMaterial>, With<borbax_ui::scene::Atom>>();
    let handles: Vec<_> = atoms.iter(world).map(|m| m.0.id()).collect();
    assert!(!handles.is_empty(), "no atoms, so this asserted nothing");
    let first = handles
        .first()
        .copied()
        .unwrap_or_else(|| unreachable!("checked non-empty above"));
    assert!(
        handles.iter().all(|h| *h == first),
        "the atoms hold {} distinct materials — a per-atom colour is a palette, \
         and a palette lands with its own G6 guard at Step 4",
        handles
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    );
}

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
        demo.embedding
            .coords()
            .iter()
            .zip(demo.embedding.radii())
            .map(|(p, r)| {
                let scale = borbax_ui::scene::narrow(r.get());
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
