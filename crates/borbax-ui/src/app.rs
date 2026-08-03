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
use bevy::window::WindowResolution;
use bevy_egui::{EguiContexts, EguiPlugin, egui};

/// The viewer's whole decision-making surface, held as a Bevy resource.
///
/// **A newtype rather than `impl Resource for ViewerState`**, because the trait
/// would have to be implemented in `state.rs` — the one file whose identity is
/// that it names no engine type at all. Keeping the impl here is what lets every
/// state test run with no Bevy app at all.
#[derive(Resource)]
struct Viewer(ViewerState);

/// Spawn the camera `egui` is drawn through.
///
/// **Without this the window opens empty and every test stays green**, which is
/// this repository's recorded failure mode for the third time — Step 1's opening
/// state, Step 2's invisible cells, and now this. `bevy_egui` paints into a
/// camera's render graph rather than to the window directly, so with no camera
/// there is nowhere for the panel to go. Nothing detects it: the frame tests
/// drive `panel::draw` through an `egui` context of their own and never involve
/// Bevy's renderer at all, so they describe a panel that is laid out correctly
/// and never reaches a pixel.
///
/// Caught by Ian looking at the window. It is the fourth defect in this project
/// found that way and the fourth that no test could have found.
fn spawn_camera(mut commands: Commands<'_, '_>) {
    commands.spawn(Camera2d);
}

/// Hand the `egui` context to the panel.
///
/// The entire bridge between the engine and the drawing, and deliberately one
/// statement: `panel::draw` is unchanged from the `eframe` build, so the frame
/// tests that describe it still describe what ships.
fn draw_panel(mut contexts: EguiContexts<'_, '_>, mut viewer: ResMut<'_, Viewer>) -> Result {
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
    egui::CentralPanel::default().show(&mut root, |ui| draw(&mut viewer.0, ui));
    Ok(())
}

/// The application, assembled but not running.
///
/// Returned rather than run so a test can drive it. `run()` is the binary's job.
pub fn viewer_app() -> App {
    build(Windowing::Real)
}

/// The same application with the windowing backend left out.
///
/// **For tests, and the difference is exactly one plugin.** `winit` permits one
/// event loop per process and creates it eagerly, so a test that built the real
/// app would fail on the second call with `RecreationAttempt` — and `cargo test`
/// runs a whole file's tests in one process.
///
/// The narrowness is the point. Everything a test asserts on — the camera, the
/// window's configuration, the resources, the schedules — is built by the same
/// code path the binary uses, so a test cannot pass over wiring the binary does
/// not have. What is *not* covered is `winit` itself: whether a real OS window
/// appears. Nothing here claims otherwise.
#[doc(hidden)]
pub fn headless_app() -> App {
    build(Windowing::None)
}

/// Whether to attach the platform's windowing backend.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Windowing {
    /// The shipped configuration: `winit` opens a real window.
    Real,
    /// `winit` omitted so more than one app can exist in a process.
    None,
}

fn build(windowing: Windowing) -> App {
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
    match windowing {
        Windowing::Real => app.add_plugins(plugins),
        Windowing::None => app.add_plugins(plugins.build().disable::<bevy::winit::WinitPlugin>()),
    }
    .add_plugins(EguiPlugin::default())
    // `opening`, not `new`: the window shows a universe from the first frame
    // rather than the words "type a seed".
    .insert_resource(Viewer(ViewerState::opening()))
    .add_systems(Startup, spawn_camera)
    .add_systems(bevy_egui::EguiPrimaryContextPass, draw_panel);
    app
}
