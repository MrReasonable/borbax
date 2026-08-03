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
/// `ui.text_edit_singleline` is the mutation they exist to catch: it fails
/// several of them and leaves **every** state-level test green, which is the
/// asymmetry that justifies the file existing. (`typing_a_name_fills_the_seed_box_in_the_frame_the_name_lands`
/// is one of them — spelled on one line so it can be grepped.)
///
/// **No fraction, deliberately.** It has been written twice and gone stale
/// twice — "three of the four", then "four of the five" in the very commit that
/// added three more frame tests, while the denominator was already eight. Both
/// times the number went stale in the same edit that widened it. The named
/// mutation and the asymmetry are the durable statement; the ratio is a fact
/// about how many tests happen to exist this week.
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
    // the new name's universe, which **§6** makes the expensive kind of wrong:
    // the seed is the shareable artefact, so a stale one is a number a child
    // writes down that opens a different universe from the one on screen. `typing_a_name_fills_the_seed_box_in_the_frame_
    // the_name_lands` is the guard, and moving this block below the seed row
    // fails it and nothing else.
    //
    // **The other direction cannot be fixed by ordering, and is fixed by
    // `request_discard` instead.** The seed box and the button both write into
    // `phrase_text` and run *after* this box has been painted, so on their own
    // they would leave the old name on screen for one frame — a universe
    // labelled with a name that did not produce it, which is exactly what
    // `Outcome` was made an enum to prevent one field along.
    //
    // Reordering does not help: it swaps which direction is broken. What does
    // help is telling `egui` that the pass it just laid out is wrong, which is
    // what `Context::request_discard` is documented for — the frame is
    // discarded and re-laid before anything is presented, so the user never
    // sees the stale name at all.
    //
    // **Two earlier versions of this comment were wrong about this and the
    // second was worse.** The first implied ordering resolved all three
    // widgets. The second said the two directions "genuinely conflict" and the
    // lag was irreducible — and shipped a test asserting the stale frame, which
    // would have made removing the defect look like a regression. A reviewer
    // measured `request_discard` closing it in a single pass.
    //
    // **The discard is gated on the name box actually having something in it,
    // and that gate is load-bearing rather than an optimisation.** Ungated, it
    // fired on every seed keystroke — and because a discarded pass re-paints
    // from updated state, it *self-repaired* the stale-status defect these
    // frame tests exist to catch. Measured: with the discard ungated, moving
    // `ui.label(state.status_line())` to the top of this function failed
    // **one** frame test; before `request_discard` existed it failed four, and
    // with the gate it fails four again. A fix that quietly disarms the guard
    // around it is worse than the defect it repaired, and this one did.
    //
    // The gate also keeps the discard to the "rare occasion" the API's own doc
    // asks for. `egui` paints a red `PERF WARNING: request_discard has been
    // called N frames in a row` at three consecutive multipass frames in debug
    // builds — which is what holding a key down in the seed box would have
    // produced for anyone running `cargo run -p borbax-ui` without
    // `--release`.
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
            // Read *before* the mutation: the discard is needed only if there
            // was a name on screen for this pass to have painted stale.
            let name_was_showing = !state.phrase_text().is_empty();
            state.randomise(moment());
            if name_was_showing {
                ui.ctx().request_discard(
                    "the surprise button cleared the name box after it was painted",
                );
            }
        }

        let response = ui
            .text_edit_singleline(state.seed_text_mut())
            .labelled_by(seed_label.id);
        // Gated on `changed()`, not called unconditionally. `Universe::generate`
        // is ~146 us and this body runs at the display's refresh rate — the
        // per-frame shape **§8.6** forbids, one layer above the simulation.
        //
        // `commit_typed_seed` rather than `reload`, which is the Step 1
        // spelling: a seed typed by hand was not produced by whatever is in the
        // name box, so the name box is cleared. **§6** makes the seed a
        // universe's address, and an address labelled with a name that did not
        // produce it is the misattribution `Outcome` exists to prevent one
        // field along. Calling `reload` here reintroduces it — and did, until
        // `typing_a_seed_by_hand_empties_the_name_box_in_the_same_frame` was written:
        // the substitution failed zero tests.
        if response.changed() {
            // See the surprise button above, and the note at the top of `draw`.
            let name_was_showing = !state.phrase_text().is_empty();
            state.commit_typed_seed();
            if name_was_showing {
                ui.ctx()
                    .request_discard("a typed seed cleared the name box after it was painted");
            }
        }
    });

    // The only string this file paints, and it is not built here — see the
    // module note. `status_line` is in `state.rs` so that a test asserting on
    // it is asserting on what the window shows.
    ui.label(state.status_line());
}
