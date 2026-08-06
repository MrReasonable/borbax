//! The drawing, and nothing that decides anything.
//!
//! **This file imports `egui` and names no engine.** `egui` lays widgets out on
//! the CPU and cannot open a window, which is what lets `egui_kittest` drive
//! [`draw`] through the accessibility tree in CI with no window and no GPU.
//!
//! **There is no `format!` in this file**, deliberately, and it is the second
//! half of the seam rule: every string the window paints comes from
//! [`crate::state::ViewerState::status_line`], so a test asserting on that
//! string is asserting on what is actually shown.

use bevy_egui::egui;

use crate::state::{ViewerState, moment};

/// Make the panel readable by a child sitting next to you.
///
/// **Applied here rather than in `main.rs`, so the headless tests see the same
/// thing.** Styling in the shell would leave `egui_kittest` driving a differently
/// sized panel from the one that ships — and geometry is exactly what the
/// reachability test asserts on.
///
/// Two defaults are wrong for this audience and both were caught by *looking at
/// the window*, which is the same way Step 1's opening state and Step 2's cell
/// spacing were caught:
///
/// - `egui`'s dark theme paints normal text at `Color32::from_gray(140)`
///   (`style.rs:1679`), a grey on near-black that is genuinely hard to read.
/// - The default text size is tuned for a developer at a laptop, not a child
///   reading a table of invented elements across a desk.
///
/// **This is not a palette and does not touch Step 4's G6 guard.** A palette
/// encodes an element *property* — affinity, mass, valence — and that is the
/// thing that could smuggle real-world chemistry in through CPK colouring.
/// Contrast and text size encode nothing about any element; they are the same
/// category as the window size. The selection highlight was ruled in on exactly
/// this argument.
fn legible(ui: &mut egui::Ui) {
    let style = ui.style_mut();

    // **Pinned to the dark theme rather than inherited.** Under the previous
    // shell this followed the system appearance, so the same code produced light
    // buttons on a dark panel — a machine-dependent look no headless test can
    // see, in a program whose whole job is being looked at. Whether the current
    // bridge would inherit it is unverified, which is exactly why the pin stays
    // and why the manual checklist still has an item for running it in a light
    // theme.
    style.visuals = egui::Visuals::dark();

    // **Labels only, and NOT `override_text_color`.** That field is a blunt
    // instrument: it repaints *every* string including button text, which put
    // near-white symbols on near-white cells and made the table invisible while
    // all 58 tests stayed green. Caught by looking at the window — twice in one
    // sitting, since the first attempt at legibility is what caused it.
    style.visuals.widgets.noninteractive.fg_stroke.color = egui::Color32::from_gray(235);

    // **A box you type into has to look like one.** The dark theme fills text
    // fields with `from_gray(10)` against a near-black panel, so the name and
    // seed boxes read as floating text with no edge — and those two boxes are
    // the only things on the screen a child is asked to *do* something with.
    style.visuals.extreme_bg_color = egui::Color32::from_gray(32);

    for font in style.text_styles.values_mut() {
        font.size *= 1.3;
    }
}

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
    legible(ui);

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

    // Not built here — see the module note. `status_line` is in `state.rs` so
    // that a test asserting on it is asserting on what the window shows.
    ui.label(state.status_line());

    ui.separator();
    periodic_table(state, ui);

    // **The scene's caption, and it is emitted here rather than beside the
    // scene for the reason this whole file exists.** It was painted from the
    // engine file at first, which is outside every frame test — so the label
    // naming the entire new feature could have returned an empty string, or the
    // wrong elements, with all 532 tests green. That is this file's own seam
    // rule defeated by a different door: the string moved to `state.rs` as the
    // rule asks, and then the widget went somewhere no test looks.
    //
    // It sits outside the 3D region, because there is no text inside it.
    ui.label(state.scene_caption());

    view_switch(state, ui);
}

/// The two ways of looking at the molecule.
///
/// **Two fixed labels rather than one button whose text changes.** A button
/// reading "show sticks" has to say what you will get, while a highlighted
/// choice says what you have — and the second is the one a child can read
/// without pressing it to find out. It also keeps this file free of `format!`,
/// which the seam requires: a string built where it is painted is a string no
/// test can assert on without describing a window it has stopped describing.
///
/// `selectable_label` rather than a radio, because the choice is between two
/// pictures rather than between two settings, and the highlight is the whole
/// message.
fn view_switch(state: &mut ViewerState, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("show");
        let view = state.view();
        // Clicking the one already showing does nothing — `flip_view` takes no
        // argument precisely so the state cannot be asked to become what it
        // already is, which would move the scene's change token for no reason.
        if ui.selectable_label(!view.draws_sticks(), "solid").clicked() && view.draws_sticks() {
            state.flip_view();
        }
        if ui.selectable_label(view.draws_sticks(), "sticks").clicked() && !view.draws_sticks() {
            state.flip_view();
        }
    });
}

/// The grid of elements, and the properties of whichever one is selected.
///
/// **The grid is emitted before the properties block, and that is the frame-order
/// rule again rather than a layout preference.** Clicking a cell changes which
/// element the block describes, so the block has to be laid out *after* the
/// click has been seen. With the block first, the frame in which a cell is
/// clicked paints the newly-highlighted cell beside the **previous** element's
/// numbers — a plausible set of properties attributed to the wrong element,
/// which is the class `Outcome` was made an enum to prevent one field along.
/// `the_selected_elements_properties_appear_in_the_frame_the_cell_is_clicked`
/// is the guard, and moving `properties` above `grid` is the mutation it exists
/// to catch.
fn periodic_table(state: &mut ViewerState, ui: &mut egui::Ui) {
    // **Deferred apply, and it is required twice over.** `rows()` borrows
    // `state` immutably for the whole grid loop, so `select` cannot be called
    // inside it; and knowing whether the selection actually *changed* — which
    // the discard below is gated on — needs the comparison made after the loop
    // rather than during it. The two requirements agree, which is the sign the
    // shape is right rather than a workaround.
    let mut clicked = None;

    // **Horizontal scrolling, because the table does not fit and must not be
    // made to.** The widest measured universe is 74 columns with a content
    // extent of ~3531 points against the window's 1200. Reflowing long rows would
    // destroy the column alignment that makes group 0, group 1 and group 2 read
    // as families down the table — measured, and it is the educational payload.
    // Shrinking to fit puts 74 cells in the window's 1200 points, i.e. 16 points
    // each, which cannot
    // hold a two-character symbol, and would hide the one thing this step exists
    // to show: that the table's *shape* differs per universe.
    egui::ScrollArea::horizontal().show(ui, |ui| {
        // **Tight spacing, and this was caught by looking at the window rather
        // than by any test.** `Grid`'s default column spacing is sized for
        // prose, which gave each two-character symbol a ~55-point column: at the
        // 900-point window of the time, only **eight** of a 58-cell row were on
        // screen (`WINDOW_SIZE` is 1200 now, and widened for the same reason), so
        // the one thing this
        // step exists to show — that the table's *shape* changes with the seed —
        // was invisible behind a scroll nobody would think to drag. Every test
        // was green over it, because the accessibility tree carries all 80 cells
        // whether or not a human can see them.
        //
        // This is a layout constant and not a physical quantity, so it is
        // exempt from the no-arithmetic rule for the same reason the window size
        // is.
        egui::Grid::new("periodic")
            .spacing(egui::vec2(3.0, 2.0))
            .show(ui, |ui| {
                for row in state.rows() {
                    ui.label(&row.shell_label);
                    for cell in &row.cells {
                        // **`selectable_label`, not a `label` inside a clickable
                        // frame.** The second works when a human clicks it and
                        // carries `Role::Label` with no `Action::Click` and no
                        // `Toggled` — invisible to a screen reader and to every
                        // test, which is why
                        // `clicking_a_cell_marks_that_cell_and_only_that_cell_as_chosen`
                        // reads the toggle state out of the accessibility tree
                        // rather than asserting on `ViewerState`.
                        //
                        // Labelled by the **symbol**, which is unique within a
                        // universe by construction. Names are not — measured, 16% of
                        // universes contain a duplicate — and a duplicate label makes
                        // `egui_kittest`'s query API panic outright.
                        if ui
                            .add(egui::Button::new(cell.symbol).selected(cell.selected))
                            .clicked()
                        {
                            clicked = Some(cell.id);
                        }
                    }
                    ui.end_row();
                }
            });
    });

    if let Some(id) = clicked {
        let changed = state.selected() != Some(id);
        state.select(id);
        // **Gated on the selection actually changing, and the gate is
        // load-bearing rather than an optimisation.** This is Step 1b's recorded
        // incident pre-empted: an ungated `request_discard` fires every frame
        // something is clicked, and because a discarded pass re-paints from
        // updated state it *self-repairs* the very frame-order defects the
        // guards above exist to catch. Measured last time: a mutation went from
        // four failures to one. `selecting_an_element_does_not_cost_a_second_layout_pass`
        // asserts the pass count so that an ungated discard is a test failure
        // rather than an invisible weakening of three other tests.
        if changed {
            ui.ctx()
                .request_discard("a cell was selected after the properties block was laid out");
        }
    }

    ui.separator();
    ui.label(state.selection_heading());
    egui::Grid::new("properties").show(ui, |ui| {
        for row in state.selection_properties() {
            ui.label(row.label);
            // Value and unit as two labels rather than one string, so the unit
            // word stays a `&'static str` the G4 allow-list can enumerate. On
            // screen they read as one phrase — `0.920` `spans`.
            ui.label(row.value);
            ui.label(row.unit);
            ui.end_row();
            if let Some(gloss) = row.gloss {
                ui.label("");
                ui.label(gloss);
                ui.end_row();
            }
        }
    });
}
