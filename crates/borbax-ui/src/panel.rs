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
/// **The statement order here is load-bearing and is pinned by tests.** `egui`
/// is immediate-mode: the seed box's response for *this* frame only exists after
/// the box has been emitted, so the status line has to come last. Emit it first
/// and it paints the previous frame's state — the window shows the old seed's
/// element count for sixteen milliseconds. Nothing in `tests/acceptance.rs` can
/// see that, because a state-level test commits and reads in one call; the
/// frame tests in `tests/panel.rs` are the guard, and moving `ui.label` above
/// `ui.text_edit_singleline` is the mutation they exist to catch. Measured with
/// `--no-fail-fast`: that mutation fails **three** of the four frame tests and
/// leaves all twenty state-level tests green.
///
/// **The button comes before the box, and that is the same rule again rather
/// than a layout preference.** "surprise me" writes into `seed_text`, so it has
/// to run before the widget that renders `seed_text` is emitted. With it after,
/// the frame showed the *previous* seed in the box beside the *new* universe's
/// count — the mirror image of the defect above, in the other widget. An earlier
/// version of this comment called that harmless because it self-corrects on the
/// next frame; a reviewer pointed out that the test claiming to cover the button
/// asserted on `state.seed_text()` and never on what was painted, so nothing
/// would have caught it getting worse. Both widgets are now asserted through the
/// accessibility tree.
pub fn draw(state: &mut ViewerState, ui: &mut egui::Ui) {
    // **The name box comes before the seed box, and that is the frame-order
    // rule again rather than a layout preference.** Typing a name writes into
    // `seed_text`, so it has to run before the widget that renders `seed_text`
    // is emitted — with it after, the frame shows the *previous* seed beside
    // the new name's universe. Same defect as the "surprise me" button's, in a
    // third widget, and the frame tests are what catch it.
    let name_label = ui.label("name");
    // Hinted rather than pre-filled. A pre-filled name chooses the universe
    // that name produces; an empty box with no hint teaches nothing about what
    // the box is for. A fixed literal like the labels around it — the "no
    // `format!`" rule is about strings built from state, which is what a test
    // can be wrong about.
    //
    // **`labelled_by` is load-bearing rather than politeness.** There are two
    // text boxes now, so `get_by_role(Role::TextInput)` is ambiguous and the
    // obvious repair is to index the tree — which pins the tests to the *order*
    // widgets are emitted in, in a file whose entire subject is that the order
    // is a decision that changes. Naming each box instead means the frame tests
    // address the one they mean, and a reordering breaks nothing that is not
    // actually broken. It is also the AccessKit relation a screen reader needs.
    let named = ui
        .add(egui::TextEdit::singleline(state.phrase_text_mut()).hint_text("type a name"))
        .labelled_by(name_label.id);
    if named.changed() {
        state.reload_from_phrase();
    }

    let seed_label = ui.label("seed");

    ui.horizontal(|ui| {
        // `moment()` is the crate's only wall-clock read, and it is called here
        // rather than inside `randomise` so that the decision half stays a pure
        // function of the number handed to it — every test of this button passes
        // a fixed moment instead of racing a clock.
        if ui.button("surprise me").clicked() {
            state.randomise(moment());
        }

        let response = ui
            .text_edit_singleline(state.seed_text_mut())
            .labelled_by(seed_label.id);
        // Gated on `changed()`, not called unconditionally. `Universe::generate`
        // is ~146 us and this body runs at the display's refresh rate.
        //
        // `commit_typed_seed` rather than `reload`, which is the Step 1
        // spelling: a seed typed by hand was not produced by whatever is in the
        // name box, so the name box is cleared. Calling `reload` here leaves a
        // universe on screen labelled with a name that did not make it.
        if response.changed() {
            state.commit_typed_seed();
        }
    });

    // The only string this file paints, and it is not built here — see the
    // module note. `status_line` is in `state.rs` so that a test asserting on
    // it is asserting on what the window shows.
    ui.label(state.status_line());
}
