//! The window, which is now a Bevy app.
//!
//! **The only file in this crate that names `bevy`**, and the only one no test
//! reaches. Everything it does is wire together: it owns the [`ViewerState`] as
//! a resource and forwards an `egui` context to [`borbax_ui::panel::draw`],
//! which is the function `egui_kittest` drives headlessly. The residual untested
//! surface is that forwarding — it takes no decision, and it is visually
//! checkable.
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

use bevy::prelude::*;
use bevy::window::WindowResolution;
use bevy_egui::{EguiContexts, EguiPlugin, egui};
use borbax_ui::{WINDOW_SIZE, WINDOW_TITLE, panel::draw, state::ViewerState};

/// The viewer's whole decision-making surface, held as a Bevy resource.
///
/// **A newtype rather than `impl Resource for ViewerState`**, because the trait
/// would have to be implemented in `state.rs` — the one file whose identity is
/// that it names no engine type at all. Keeping the impl here is what lets every
/// state test run with no Bevy app at all.
#[derive(Resource)]
struct Viewer(ViewerState);

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

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: WINDOW_TITLE.to_owned(),
                // Read from the constant rather than typed here, so the frame
                // test asserting every cell is reachable is asserting about
                // *this* window and not one it guessed.
                resolution: WindowResolution::new(
                    u32::from(WINDOW_SIZE[0]),
                    u32::from(WINDOW_SIZE[1]),
                ),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin::default())
        // `opening`, not `new`: the window shows a universe from the first
        // frame rather than the words "type a seed".
        .insert_resource(Viewer(ViewerState::opening()))
        .add_systems(bevy_egui::EguiPrimaryContextPass, draw_panel)
        .run();
}
