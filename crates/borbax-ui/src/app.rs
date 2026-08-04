//! The Bevy app: what the window is made of.
//!
//! **Split out of `main.rs` because a test cannot import a binary.** The empty
//! window that shipped on 2026-08-03 was an app missing a camera, and no test
//! could see it for a structural reason rather than an oversight: the wiring
//! lived in `main.rs`, which is a binary crate root. Any test would have had to
//! rebuild the app itself and would then have been asserting about its own copy.
//!
//! With the construction here, `viewer_app` is the *same* value the binary runs,
//! and `tests/app.rs` drives it.
//!
//! # The engine displays; it does not compute
//!
//! Ian's ruling, 2026-08-03. Bevy owns the window, the camera, the meshes and
//! the input. The physics, the chemistry and the generation are ours and stay in
//! the `borbax-*` crates. `check_the_engine_stays_in_the_viewer` fails the build
//! if any other crate so much as names an engine in its manifest.
//!
//! That is not tidiness. Measured on this workspace: Bevy reuses entity slots,
//! so after a despawn/respawn cycle — exactly what rebuilding the scene on a
//! seed change does — query iteration comes back as `[49, 1, 2, 62, 4, 5, 61,
//! ...]` from entities spawned in order. Anything result-affecting that read
//! that order would diverge across the §13.4 matrix with every test green.

use crate::panel::draw;
use crate::state::ViewerState;
use crate::{WINDOW_SIZE, WINDOW_TITLE};
use bevy::prelude::*;
use bevy::render::RenderPlugin;
use bevy::render::settings::{RenderCreation, WgpuSettings};
use bevy::window::WindowResolution;
use bevy_egui::{EguiContexts, EguiGlobalSettings, EguiPlugin, egui};

/// The viewer's whole decision-making surface, held as a Bevy resource.
///
/// **A newtype rather than `impl Resource for ViewerState`**, because the trait
/// would have to be implemented in `state.rs` — the one file whose identity is
/// that it names no engine type at all. Keeping the impl here is what lets every
/// state test run with no Bevy app at all.
#[derive(Resource, Debug)]
pub struct Viewer(pub ViewerState);

/// Hand the `egui` context to the panel.
///
/// The entire bridge between the engine and the drawing, and deliberately one
/// statement: `panel::draw` is unchanged from the `eframe` build, so the frame
/// tests that describe it still describe what ships.
fn draw_panel(
    mut contexts: EguiContexts<'_, '_>,
    mut viewer: ResMut<'_, Viewer>,
    mut region: ResMut<'_, crate::scene::SceneRegion>,
) -> Result {
    let ctx = contexts.ctx_mut()?;
    // **The root `Ui` is built here rather than handed to us**, and that is a
    // real difference from the `eframe` build rather than boilerplate. egui 0.35
    // moved panels from taking a `&Context` to taking a `&mut Ui`, so `eframe`'s
    // `fn ui` was giving us a root `Ui` it had made. Bevy's bridge gives a
    // `Context`, so the root is ours to make — this is the pattern `bevy_egui`'s
    // own examples use.
    let mut root = egui::Ui::new(
        ctx.clone(),
        "viewport".into(),
        egui::UiBuilder::new()
            .layer_id(egui::LayerId::background())
            .max_rect(ctx.viewport_rect()),
    );
    // **A top panel rather than a central one, and the whole scene depends on
    // it.** `CentralPanel` fills the viewport with an opaque background, so a 3D
    // scene rendered underneath is painted over completely — with every test
    // green, because nothing in the suite looks at a pixel. Reserving the panel
    // at the top leaves a hole, and the hole is what the scene camera is aimed
    // at.
    egui::Panel::top("controls").show(&mut root, |ui| draw(&mut viewer.0, ui));

    // **Returned rather than agreed by convention.** The camera is aimed at
    // whatever the panel actually left, so the two cannot drift apart — and a
    // test can assert on the region without hunting through the widget tree.
    let left = root.available_rect_before_wrap();
    region.0 = Some(Rect::new(left.min.x, left.min.y, left.max.x, left.max.y));
    Ok(())
}

/// The application, assembled but not running.
///
/// Returned rather than run so a test can drive it. `run()` is the binary's job.
pub fn viewer_app() -> App {
    build(Shell::Window)
}

/// The same application with the platform backends left out.
///
/// **For tests, and the difference is two named things rather than one.** An
/// earlier version of this function differed from [`viewer_app`] by exactly one
/// plugin and said so, because that narrowness is what stops a test asserting
/// about wiring the binary does not have. It now differs by two, and both are
/// forced:
///
/// - **`winit`**, because it permits one event loop per process and creates it
///   eagerly, so the second test in a file would fail with `RecreationAttempt`.
/// - **the graphics backend**, because Bevy asks the platform for a real device
///   while the app is being built and calls `.expect("Unable to find a GPU!")`
///   when there is none. A runner without a driver therefore aborts every test
///   in `tests/app.rs` inside `bevy_render`, with a message that reads like a
///   defect in this crate. Measured on this machine before the change: the
///   headless app held a live `RenderDevice`.
///
/// What that costs is stated rather than discovered later: **no assertion here
/// can ever be about a pixel.** That was never available on a driverless runner
/// anyway, so nothing is given up that CI could have had. Everything these tests
/// do assert on — the cameras, the window's configuration, the resources, the
/// schedules, entity transforms, visibility and projection — is built by the
/// same code path the binary uses and is unaffected by there being no device.
#[doc(hidden)]
pub fn headless_app() -> App {
    build(Shell::Headless)
}

/// Which platform backends the app attaches.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Shell {
    /// The shipped configuration: a real window and a real graphics device.
    Window,
    /// Neither, so many apps can exist in one process and none needs a driver.
    Headless,
}

fn build(shell: Shell) -> App {
    let mut app = App::new();
    let plugins = DefaultPlugins.set(WindowPlugin {
        primary_window: Some(Window {
            title: WINDOW_TITLE.to_owned(),
            // Read from the constant rather than typed here, so the frame test
            // asserting every cell is reachable is asserting about *this*
            // window and not one it guessed.
            resolution: WindowResolution::new(u32::from(WINDOW_SIZE[0]), u32::from(WINDOW_SIZE[1])),
            ..default()
        }),
        ..default()
    });
    match shell {
        Shell::Window => app.add_plugins(plugins),
        // `backends: None` rather than dropping `RenderPlugin`: the render
        // plugin is what registers `Mesh3d`, `MeshMaterial3d`, `Camera3d`,
        // visibility and frustum culling, so removing it would take the whole
        // scene out of reach of every test. Asking for no backend keeps all of
        // that and skips only the device.
        Shell::Headless => app.add_plugins(
            plugins
                .build()
                .disable::<bevy::winit::WinitPlugin>()
                .set(RenderPlugin {
                    render_creation: RenderCreation::Automatic(Box::new(WgpuSettings {
                        backends: None,
                        ..default()
                    })),
                    ..default()
                }),
        ),
    }
    .add_plugins(EguiPlugin::default())
    // `opening`, not `new`: the window shows a universe from the first frame
    // rather than the words "type a seed".
    .insert_resource(Viewer(ViewerState::opening()))
    // **The engine is not allowed to pick which camera egui draws through.**
    // Left automatic it attaches to whichever camera is created first, so the
    // render topology would be decided by spawn order inside a startup system.
    // `scene::plugin` places it deliberately, on the full-window camera.
    .insert_resource(EguiGlobalSettings {
        auto_create_primary_context: false,
        ..default()
    })
    .add_systems(bevy_egui::EguiPrimaryContextPass, draw_panel);
    crate::scene::plugin(&mut app);
    app
}
