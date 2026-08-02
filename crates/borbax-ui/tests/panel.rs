//! Frame-level tests, driven through the accessibility tree with no window and
//! no GPU.
//!
//! **These exist for one defect the tests in `acceptance.rs` structurally
//! cannot see.** A state-level test commits a seed and reads the line back in
//! the same call, so it is green whatever order the panel emits its widgets in.
//! An immediate-mode panel that paints the label *before* the seed box reads
//! last frame's state — the window shows the previous seed's count for one
//! frame, sixteen milliseconds, which nobody demonstrating the app will catch
//! and no state-level assertion can reach.

use borbax_ui::{panel::draw, state::ViewerState};
use borbax_universe::Universe;
use egui::accesskit::Role;
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;

/// A harness driving the real [`draw`] over a real [`ViewerState`].
fn harness<'a>() -> Harness<'a, ViewerState> {
    Harness::new_ui_state(|ui, state| draw(state, ui), ViewerState::new())
}

/// Focus a box by the label the panel gives it, and type into it.
///
/// **Addressed by label rather than by position in the tree.** Step 1b added a
/// second `TextInput`, which made the bare `get_by_role(Role::TextInput)` these
/// tests used ambiguous. The cheap repair is `.nth(1)`, and it would tie this
/// file — whose entire subject is that widget *order* is a decision that
/// changes — to the very ordering it exists to let people change. The panel
/// wires `Response::labelled_by` so each box can be named instead.
///
/// `focus()` before `type_text` is not optional: `type_text` queues an
/// `egui::Event::Text` addressed to whatever has focus, so without it the text
/// never lands and the test fails looking *exactly* like the defect it hunts.
fn type_into(harness: &Harness<'_, ViewerState>, label: &str, text: &str) {
    let field = harness.get_by_role_and_label(Role::TextInput, label);
    field.focus();
    field.type_text(text);
}

#[test]
fn the_element_count_appears_in_the_frame_the_seed_lands() {
    let mut harness = harness();

    type_into(&harness, "seed", "7");

    // **Exactly one frame.** A `run()`-shaped call steps until the UI is
    // quiescent, which hides this defect by definition — the second frame
    // repairs the lag. A harness helper that loops is the vacuity trap in this
    // test, and it is the reason `step()` is spelled out here rather than
    // wrapped in something convenient.
    harness.step();

    // Generated once and read twice, matching `expected_for` in
    // `tests/acceptance.rs`. Two `generate(7)` calls would agree — it is a pure
    // function — but they would agree by luck of purity rather than by
    // construction, and the shape is the one that stops agreeing the day
    // anything about generation becomes context-dependent.
    let universe = Universe::generate(7);
    let expected = format!(
        "{} elements · physics v{}",
        universe.table.len(),
        u8::from(universe.physics)
    );
    assert!(
        harness.query_by_label(&expected).is_some(),
        "one frame after typing `7` the panel does not show {expected:?} — \
         the label is being emitted before the seed box's response is read"
    );
}

#[test]
fn an_unchanged_seed_box_does_not_rebuild_the_universe() {
    // `Universe::generate` is ~146 us. At 60 Hz an unconditional call in the
    // paint body spends most of the frame budget regenerating a universe
    // nobody asked for again, and by Step 2 there are 120 element cells behind
    // it. This is §8.6's per-species rule one level up.
    let mut harness = harness();

    type_into(&harness, "seed", "7");
    harness.step();
    assert_eq!(
        harness.state().regenerations(),
        1,
        "committing one seed should generate exactly one universe"
    );

    for _ in 0..30 {
        harness.step();
    }
    assert_eq!(
        harness.state().regenerations(),
        1,
        "30 idle frames rebuilt the universe — `reload` is being called from \
         the paint body rather than from the seed box's `changed()` response"
    );
}

#[test]
fn the_surprise_button_fills_the_seed_box_and_loads_that_universe() {
    // The button's contract is visible-and-retypable (see
    // `the_random_button_puts_a_seed_in_the_box_that_can_be_typed_back`). This
    // is the half that can only be checked against the real widget tree: that
    // the button exists, that clicking it is what triggers this, and that the
    // count lands in the same frame rather than one behind.
    let mut harness = harness();
    assert!(
        harness.state().seed_text().is_empty(),
        "this test starts from an empty box so the button is what fills it"
    );

    harness.get_by_label("surprise me").click();
    harness.step();

    let text = harness.state().seed_text().to_owned();
    assert!(
        !text.is_empty(),
        "clicking `surprise me` left the box empty"
    );

    let expected = format!("{} elements", {
        let seed: u64 = text.parse().unwrap_or_else(|_| {
            unreachable!("`randomise` writes a plain u64, and it wrote {text}")
        });
        Universe::generate(seed).table.len()
    });
    assert!(
        harness.query_by_label_contains(&expected).is_some(),
        "the frame the button was clicked in does not show {expected:?} for seed {text}"
    );

    // **The rendered box, not the state behind it.** The first version of this
    // test asserted only on `harness.state().seed_text()`, which is true one
    // instruction after `randomise` runs and says nothing about what the window
    // painted. It missed a real inconsistency: with the button emitted *after*
    // the text box, the box had already rendered the old text for that frame,
    // so the screen showed the previous seed next to the new universe's count —
    // the mirror of the defect `the_element_count_appears_in_the_frame_the_seed_lands`
    // exists for, in the other widget. Asserting the accessibility value is what
    // makes the two agree in the frame the user is looking at.
    let rendered = harness
        .get_by_role_and_label(Role::TextInput, "seed")
        .value()
        .unwrap_or_else(|| unreachable!("a TextInput node always carries a value"));
    assert_eq!(
        rendered, text,
        "the seed box renders {rendered:?} while the state and the count say {text:?} — \
         the box was emitted before the button that changed it"
    );
}

#[test]
fn a_refused_seed_shows_the_refusal_in_the_panel_itself() {
    // The `acceptance.rs` version of this asserts on `status_line`. This one
    // asserts the panel actually paints it, which is the half that would rot
    // silently if the panel ever stopped calling `status_line`.
    let mut harness = harness();

    // Still the *seed* box, and still valid under two boxes: the seed box keeps
    // its four refusals precisely because the name box exists to absorb
    // anything that is not a number. With one control this test would have had
    // to go, and `letters_in_the_seed_box_are_refused_not_hashed` with it.
    type_into(&harness, "seed", "abc");
    harness.step();

    assert!(
        harness
            .query_by_label("seeds are whole numbers, 0 to 18446744073709551615")
            .is_some(),
        "the panel does not show the refusal for a non-numeric seed"
    );
}

/// Typing a name fills the seed box **in the same frame**.
///
/// **The third instance of this file's one defect, in a third widget.** The
/// name box writes into `seed_text`, so it has to be emitted before the widget
/// that renders `seed_text`. Emitted after, the frame shows the *previous*
/// seed beside the new name's element count — a number the user can write down
/// that opens a different universe from the one on screen, which is the exact
/// thing §6 makes the seed for.
///
/// A state-level test cannot see this: `set_phrase` commits and reads in one
/// call, and every one of those is green whatever order the panel paints in.
#[test]
fn typing_a_name_fills_the_seed_box_in_the_frame_the_name_lands() {
    let mut harness = harness();
    type_into(&harness, "name", "Emily");
    harness.step();

    let seed = borbax_universe::seed_from_phrase("Emily")
        .unwrap_or_else(|| unreachable!("`Emily` is not blank, so it names a universe"));
    let rendered = harness
        .get_by_role_and_label(Role::TextInput, "seed")
        .value()
        .unwrap_or_else(|| unreachable!("a TextInput node always carries a value"));
    assert_eq!(
        rendered,
        seed.to_string(),
        "the seed box renders {rendered:?} one frame after `Emily` was typed — \
         it was emitted before the name box that changed it"
    );

    let expected = format!("{} elements", Universe::generate(seed).table.len());
    assert!(
        harness.query_by_label_contains(&expected).is_some(),
        "the frame `Emily` landed in does not show {expected:?}"
    );
}
