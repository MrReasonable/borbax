//! The drawing, and nothing that decides anything.
//!
//! **This file imports `egui` and never `eframe`.** `egui` lays widgets out on
//! the CPU and cannot open a window, which is what lets `egui_kittest` drive
//! [`draw`] through the accessibility tree in CI with no window and no GPU.
//!
//! **There is no `format!` in this file**, deliberately, and it is the second
//! half of the seam rule: every string the window paints comes from
//! [`crate::state::ViewerState::status_line`], so a test asserting on that
//! string is asserting on what is actually shown.

use crate::state::{ViewerState, moment};

/// Draw the whole viewer into `ui`.
///
/// **The order of these four statements is load-bearing and is pinned by a
/// test.** `egui` is immediate-mode: the seed box's response for *this* frame
/// only exists after the box has been emitted, so the status line has to come
/// last. Emit it first and it paints the previous frame's state — the window
/// shows the old seed's element count for sixteen milliseconds. Nothing in
/// `tests/acceptance.rs` can see that, because a state-level test commits and
/// reads in one call; `the_element_count_appears_in_the_frame_the_seed_lands`
/// in `tests/panel.rs` is the guard, and moving `ui.label` above
/// `ui.text_edit_singleline` is the mutation it exists to catch.
pub fn draw(state: &mut ViewerState, ui: &mut egui::Ui) {
    ui.label("seed");

    let response = ui.text_edit_singleline(state.seed_text_mut());
    // Gated on `changed()`, not called unconditionally. `Universe::generate` is
    // ~146 us and this body runs at the display's refresh rate.
    if response.changed() {
        state.reload();
    }

    // Above the status line for the same reason the seed box is: whatever this
    // changes has to be visible in *this* frame, not the next one.
    //
    // `moment()` is the crate's only wall-clock read, and it is called here
    // rather than inside `randomise` so that the decision half stays a pure
    // function of the number handed to it — every test of this button passes a
    // fixed moment instead of racing a clock.
    if ui.button("surprise me").clicked() {
        state.randomise(moment());
    }

    // The only string this file paints, and it is not built here — see the
    // module note. `status_line` is in `state.rs` so that a test asserting on
    // it is asserting on what the window shows.
    ui.label(state.status_line());
}
