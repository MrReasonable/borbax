//! The window.
//!
//! **The only file in this crate that names `eframe`**, and the only one no
//! test reaches. Everything it does is forward: `fn ui` hands `egui::Ui` to
//! [`borbax_ui::panel::draw`], which is the function `egui_kittest` drives
//! headlessly. The residual untested surface is that one forwarding line —
//! it takes no decision, and it is visually checkable.

use borbax_ui::{WINDOW_TITLE, panel::draw, state::ViewerState};

/// The application, which is a [`ViewerState`] and a way to draw it.
struct Viewer {
    state: ViewerState,
}

impl eframe::App for Viewer {
    // `fn ui`, not `fn update`. eframe 0.35 splits `App` into a provided
    // `fn logic(&mut self, ctx, frame)` and a required
    // `fn ui(&mut self, ui, frame)` — every tutorial, and every reflex, writes
    // `update`, which gives E0407 (`method not a member of trait`) plus E0046
    // (`missing fn ui`) and reads like a version mismatch rather than a rename.
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        draw(&mut self.state, ui);
    }
}

fn main() -> eframe::Result {
    eframe::run_native(
        WINDOW_TITLE,
        eframe::NativeOptions::default(),
        Box::new(|_cc| {
            Ok(Box::new(Viewer {
                // `opening`, not `new`: the window shows a universe from the
                // first frame rather than the words "type a seed".
                state: ViewerState::opening(),
            }))
        }),
    )
}
