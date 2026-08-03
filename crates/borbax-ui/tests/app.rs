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

    // The status line is built in `state.rs`, which is what makes asserting on
    // it an assertion about what the window paints.
    let opening = borbax_ui::state::ViewerState::opening();
    assert!(
        opening.status_line().contains("elements"),
        "the opening state does not describe a universe, so the window's first \
         frame shows no chemistry: {:?}",
        opening.status_line()
    );
}
