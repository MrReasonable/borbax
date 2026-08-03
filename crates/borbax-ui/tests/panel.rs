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
use egui_kittest::kittest::{NodeT as _, Queryable};

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
    // The other half of the discard gate: with the name box empty there is
    // nothing stale to repaint, so no discard is wanted and one pass is right.
    // Ungated, this reads 2 — which is the assertion that fails when the gate
    // goes.
    assert_eq!(
        harness.output().platform_output.num_completed_passes,
        1,
        "a seed typed with no name showing should not cost a second layout pass"
    );

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

/// An idle name box does not rebuild the universe every frame.
///
/// **The twin of `an_unchanged_seed_box_does_not_rebuild_the_universe`, and it
/// was missing.** Measured by a review lane: replacing `if named.changed()`
/// with an unconditional `state.reload_from_phrase();` failed **zero** tests
/// across the whole crate. `naming_a_universe_generates_exactly_one` cannot see
/// it, because it calls `set_phrase` directly and never goes through `draw`.
///
/// `Universe::generate` is ~146 µs and this body runs at the display's refresh
/// rate, which is the same argument the seed box's version makes.
#[test]
fn an_unchanged_name_box_does_not_rebuild_the_universe() {
    let mut harness = harness();
    type_into(&harness, "name", "Emily");
    harness.step();
    assert_eq!(
        harness.state().regenerations(),
        1,
        "typing one name should generate exactly one universe"
    );

    for _ in 0..30 {
        harness.step();
    }
    assert_eq!(
        harness.state().regenerations(),
        1,
        "30 idle frames rebuilt the universe — `reload_from_phrase` is being \
         called from the paint body rather than from the name box's `changed()` \
         response"
    );
}

/// What the name box currently renders, read from the accessibility tree.
///
/// Reads what was *painted*, not `state.phrase_text()` — a state-level read is
/// true one instruction after the method runs and says nothing about the frame
/// the user is looking at, which is the whole distinction this file exists for.
fn name_box(harness: &Harness<'_, ViewerState>) -> String {
    harness
        .get_by_role_and_label(Role::TextInput, "name")
        .value()
        .unwrap_or_else(|| unreachable!("a TextInput node always carries a value"))
}

/// Typing a seed by hand empties the name box **in the same presented frame**.
///
/// **The wiring, which nothing tested.** Measured by a review lane: reverting
/// `panel.rs` to Step 1's `state.reload()` — one identifier's difference from
/// `state.commit_typed_seed()` — left every test that then existed green, which
/// is why this test was written. (A "refresh" of that historical figure to the
/// then-current count was wrong three ways at once, and a reviewer measured all
/// three: the measurement was against 24 tests, the crate held a different
/// number by then, and the sentence claimed this test does *not* catch the
/// mutation it demonstrably does. Numbers of this shape are dropped on sight
/// now — see `draw`'s doc for the same rule.) The
/// acceptance test that looks like it covers this calls `commit_typed_seed`
/// directly, so it tests the method and never the call.
///
/// That mutation is now a **compile error** rather than a test failure, because
/// `reload` is private: `error[E0624]: method `reload` is private`. This test
/// still earns its place — it covers the frame, which visibility cannot.
///
/// **One `step()`, and an earlier version of this test asserted two.** The name
/// box is emitted before the seed box, so clearing it mid-pass would leave the
/// stale name on screen for a frame — and the first version of this test
/// asserted that stale frame, pinning the defect as though it were a law. A
/// reviewer showed `Context::request_discard` removes it outright: the pass is
/// thrown away and re-laid before anything is presented. Asserting the lag
/// would have made fixing it look like a regression.
///
/// **`type_text` appends, so the seed typed here is deliberately a refusal.**
/// The box already holds the twenty digits `Emily` produced, so `"4"` makes a
/// twenty-one digit number — `state.rs` records this exact trap as the reason
/// `ViewerState::opening()` exists separately from `new()`. That is worth
/// asserting rather than working around: the name must stop claiming the seed
/// **even when the seed is refused**, because a refusal is precisely when a
/// stale name would be left sitting over an error message. The valid-seed frame
/// path is covered by `the_element_count_appears_in_the_frame_the_seed_lands`.
#[test]
fn typing_a_seed_by_hand_empties_the_name_box_in_the_same_frame() {
    let mut harness = harness();
    type_into(&harness, "name", "Emily");
    harness.step();
    assert_eq!(name_box(&harness), "Emily");

    type_into(&harness, "seed", "4");
    harness.step();
    assert_eq!(
        name_box(&harness),
        "",
        "the frame shows a name that did not produce the seed beside it"
    );
    // **The gate on `request_discard`, which nothing else covers.** Deleting
    // `if name_was_showing` failed 0 of the 457 tests the workspace then held —
    // and an ungated discard is
    // what round 2 shipped, so it is a demonstrated regression rather than a
    // hypothetical: it silently re-repairs the stale-status defect that three
    // of the four frame-order guards exist to catch.
    //
    // A discard costs a second layout pass, so the pass count is the one
    // observable that tells a gated discard from an ungated one. Here the name
    // box *was* showing, so the discard is correct and there are two.
    assert_eq!(
        harness.output().platform_output.num_completed_passes,
        2,
        "clearing a showing name box should discard the pass that painted it"
    );

    // The frame-level analogue of `the_boxes_agree`: whatever the seed box
    // holds, the line underneath describes *that*, with nothing left over from
    // the name.
    let seed_shown = harness
        .get_by_role_and_label(Role::TextInput, "seed")
        .value()
        .unwrap_or_else(|| unreachable!("a TextInput node always carries a value"));
    // **Derived, not written down.** The first version hard-coded the
    // twenty-one digits, which smuggles a copy of `phrase.rs`'s golden into a
    // UI frame test — a deliberate remapping would then fail here with a
    // message about frame timing. The sibling test above already shows the
    // right shape.
    let emily = borbax_universe::seed_from_phrase("Emily")
        .unwrap_or_else(|| unreachable!("`Emily` is not blank, so it names a universe"));
    assert_eq!(
        seed_shown,
        format!("{emily}4"),
        "type_text appends, so this should be Emily's seed with a `4` on the end"
    );
    assert!(
        harness
            .query_by_label("that is larger than the largest seed, 18446744073709551615")
            .is_some(),
        "the refused seed's message is not on screen"
    );
}

/// "surprise me" empties the name box in the same frame too, for the same
/// reason — one `step()`, because `request_discard` re-lays the pass.
#[test]
fn a_surprise_seed_empties_the_name_box_in_the_same_frame() {
    let mut harness = harness();
    type_into(&harness, "name", "Emily");
    harness.step();

    harness.get_by_label("surprise me").click();
    harness.step();
    assert_eq!(
        name_box(&harness),
        "",
        "the name box still claims to have produced a seed the button drew"
    );
}

/// A name is not truncated by the box it is typed into, **at any length**.
///
/// **The other half of the length-cap hazard, and it lives here rather than in
/// `borbax-universe`.** That crate's
/// `no_prefix_of_a_phrase_is_enough_to_decide_the_universe` guards the byte
/// loop; nothing guarded the widget. A reviewer measured `.char_limit(8)` on
/// the name box passing the **entire** workspace suite while collapsing two
/// different full names onto one universe, because no test above four
/// characters ever went through `draw`.
///
/// **The first version of this test compared two names sharing a 20-character
/// stem, and a reviewer showed that discriminates caps of 20 or less and
/// nothing above — `.char_limit(512)` passed the whole suite.** That is the
/// *same* defect, in the sibling file, committed in the same pass that repaired
/// it: a guard bounded by a length someone guessed.
///
/// So this asserts the property directly — whatever was typed is what the box
/// holds — over a name longer than any `char_limit` anyone would write. **The
/// class it decides, stated so the next person can see its edge:** every cap at
/// or below **2110**, one below the length typed — a cap at exactly the typed
/// length truncates nothing. A reviewer measured that edge rather than taking
/// "at or below the length typed here", which was off by one. The length is
/// chosen to sit above the caps a person reaches for (8, 32, 64, 128, 256, 512,
/// 1024) rather than above a name anyone would type. Measured at 8, 20, 512 and
/// 1024: all four fail.
#[test]
fn a_name_is_not_truncated_by_the_box_it_is_typed_into() {
    let mut harness = harness();
    let typed = "elizabeth alexandra mary windsor ".repeat(64);
    let typed = typed.trim_end();
    type_into(&harness, "name", typed);
    harness.step();

    assert_eq!(
        harness.state().phrase_text(),
        typed,
        "the name box dropped characters — it has acquired a char_limit, and \
         every name longer than it now shares a universe with its truncation"
    );
}

// ---------------------------------------------------------------------------
// Step 2 — the periodic table panel.
// ---------------------------------------------------------------------------

/// A harness wide enough that every cell of the widest table is on screen.
///
/// **3600 px, and the first figure written here was 3000, which was wrong.**
/// The worst-case content extent is ~3531 px and it does **not** depend on the
/// harness size: `egui`'s `Grid` sizes columns to their content, so there is no
/// reflow to rescue a harness that is too narrow. At 3000 the far-right cell is
/// off screen and `Node::click()` — which synthesises a pointer event at the
/// node's rect centre — is a **silent no-op**, leaving the selection unchanged
/// and the test green for the wrong reason.
///
/// This is deliberately *not* the size the window opens at. Clipping is noise
/// when the subject is placement and count; the test that cares about the real
/// window reads [`borbax_ui::WINDOW_SIZE`] instead.
fn wide_harness<'a>() -> Harness<'a, ViewerState> {
    Harness::builder()
        .with_size(egui::vec2(3600.0, 900.0))
        .build_ui_state(|ui, state| draw(state, ui), ViewerState::new())
}

/// Three distinct symbols from a seed's widest row: one to click, two to check
/// stay unselected.
///
/// Returned as a tuple rather than indexed out of the vector at the call site,
/// because `clippy::indexing_slicing` is denied workspace-wide and a panicking
/// index in a test is still a panicking index — it turns a clear assertion
/// failure into a slice-out-of-bounds with no message.
fn three_symbols(seed: u64) -> (String, String, String) {
    let symbols = widest_row_symbols(seed);
    let mut it = symbols.into_iter();
    let (a, b, c) = (it.next(), it.next(), it.next());
    assert!(
        c.is_some(),
        "seed {seed}'s widest row is too short to distinguish a selected cell from \
         its neighbours"
    );
    (
        a.unwrap_or_default(),
        b.unwrap_or_default(),
        c.unwrap_or_default(),
    )
}

/// The symbols of a seed's widest row, longest row first.
fn widest_row_symbols(seed: u64) -> Vec<String> {
    let mut state = ViewerState::new();
    seed.to_string().clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    let rows = state.rows();
    let widest = rows
        .iter()
        .max_by_key(|r| r.cells.len())
        .map(|r| r.cells.iter().map(|c| c.symbol.to_owned()).collect());
    widest.unwrap_or_default()
}

/// Clicking a cell marks that cell, and only that cell, as chosen.
///
/// **Read out of the accessibility tree, not out of `ViewerState`.** A cell
/// drawn as `ui.label` inside a clickable region works perfectly when a human
/// clicks it, and carries `Role::Label` with no `Action::Click` and no toggle
/// state — invisible to a screen reader and to every test that asks the state
/// instead of the tree. There is no state-level version of this test to write.
#[test]
fn clicking_a_cell_marks_that_cell_and_only_that_cell_as_chosen() {
    let mut harness = wide_harness();
    type_into(&harness, "seed", "7");
    harness.step();

    let (first, chosen_symbol, third) = three_symbols(7);

    harness
        .get_by_role_and_label(Role::Button, &chosen_symbol)
        .click();
    harness.step();

    // Read through `accesskit_node()`, which is the tree a screen reader sees.
    // `kittest::Node` exposes the AccessKit node rather than re-exporting every
    // accessor, so this is the route to the toggle state.
    let chosen = harness.get_by_role_and_label(Role::Button, &chosen_symbol);
    assert_eq!(
        chosen.accesskit_node().toggled(),
        Some(egui::accesskit::Toggled::True),
        "the clicked cell does not report itself as chosen — a screen reader, and \
         every test, would see nothing"
    );
    for other in [&first, &third] {
        assert_eq!(
            harness
                .get_by_role_and_label(Role::Button, other)
                .accesskit_node()
                .toggled(),
            Some(egui::accesskit::Toggled::False),
            "{other} is drawn chosen as well, so either every cell is highlighted or \
             the selection is not exclusive"
        );
    }
}

/// The selected element's properties appear in the frame its cell is clicked.
///
/// **This file's one defect, in Step 2's widget.** With the properties block
/// emitted above the grid, the frame in which a cell is clicked paints the
/// newly-highlighted cell beside the *previous* element's numbers. Mutation:
/// move the `properties` block above the `grid` block in `panel::periodic_table`.
#[test]
fn the_selected_elements_properties_appear_in_the_frame_the_cell_is_clicked() {
    let mut harness = wide_harness();
    type_into(&harness, "seed", "7");
    harness.step();

    let (_, chosen_symbol, _) = three_symbols(7);

    // Exactly one frame after the click. A second would repair the lag and hide
    // the defect, which is the whole reason `step()` is spelled out.
    harness
        .get_by_role_and_label(Role::Button, &chosen_symbol)
        .click();
    harness.step();

    let universe = Universe::generate(7);
    let element = universe
        .table
        .iter()
        .find(|(_, e)| e.symbol == chosen_symbol)
        .map(|(_, e)| e);
    assert!(element.is_some(), "the fixture symbol left the table");
    let Some(element) = element else { return };

    // The heading is the cheapest thing to read and it moves with the
    // selection, so a stale block shows the *previous* element's name here.
    harness.get_by_label(&format!("{} · {}", element.symbol, element.name));
}

/// The grid paints as many rows and columns as this universe has.
///
/// Kills `for p in 0..7 { for g in 0..18 }` in the paint body while the state's
/// own `rows()` is correct — the real periodic table's shape, which is wrong in
/// every Borbax universe: measured, periods run 2–4 and groups run up to 74.
#[test]
fn the_grid_paints_as_many_rows_and_columns_as_this_universe_has() {
    for seed in [1_u64, 7, 42] {
        let mut harness = wide_harness();
        type_into(&harness, "seed", &seed.to_string());
        harness.step();

        let universe = Universe::generate(seed);
        for (_, element) in universe.table.iter() {
            // Every element of the table has a cell that can be found by its
            // symbol. `get_by_role_and_label` panics if it finds none — or two.
            harness.get_by_role_and_label(Role::Button, &element.symbol);
        }
    }
}

/// An idle table does not rebuild the universe.
///
/// The Step 2 twin of `an_unchanged_name_box_does_not_rebuild_the_universe`,
/// whose own note records that the state-level version failed **zero** tests
/// because it never went through `draw`. The identical hole exists for the grid:
/// a `Universe::generate` in the paint body increments no counter and is
/// invisible to every state-level assertion.
#[test]
fn an_idle_table_does_not_rebuild_the_universe() {
    let mut harness = wide_harness();
    type_into(&harness, "seed", "7");
    harness.step();
    let after_commit = harness.state().regenerations();

    for _ in 0..30 {
        harness.step();
    }
    assert_eq!(
        harness.state().regenerations(),
        after_commit,
        "the table rebuilt the universe while nothing was happening"
    );
}

/// Selecting an element costs no second layout pass.
///
/// **Kills an ungated `request_discard` anywhere in the selection path**, and
/// this is Step 1b's recorded incident pre-empted rather than repeated. An
/// ungated discard re-paints from updated state, which *self-repairs* the
/// frame-order defects three other tests exist to catch — measured then, a
/// mutation went from four failures to one. Nothing but the pass count can see
/// it.
#[test]
fn selecting_an_element_does_not_cost_a_second_layout_pass() {
    let mut harness = wide_harness();
    type_into(&harness, "seed", "7");
    harness.step();

    let (only, _, _) = three_symbols(7);

    // Clicking the *already selected* cell must not discard: nothing changed,
    // so there is nothing stale to repaint.
    harness.get_by_role_and_label(Role::Button, &only).click();
    harness.step();
    harness.get_by_role_and_label(Role::Button, &only).click();
    harness.step();
    assert_eq!(
        harness.output().platform_output.num_completed_passes,
        1,
        "re-selecting the same cell asked for a second layout pass, so the discard \
         is not gated on the selection actually changing"
    );
}

/// Two cells never answer to the same label.
///
/// Kills labelling cells by group index or period, which makes every
/// `get_by_label` in this file ambiguous — and `kittest`'s query API panics on
/// an ambiguous match, so the failure would land in whichever test ran first
/// rather than where the mistake is. Symbols are unique within a universe by
/// construction; this proves the panel used them.
#[test]
fn two_cells_never_answer_to_the_same_label() {
    for seed in 1..40_u64 {
        let mut state = ViewerState::new();
        seed.to_string().clone_into(state.seed_text_mut());
        state.commit_typed_seed();
        let labels: Vec<&str> = state
            .rows()
            .iter()
            .flat_map(|r| r.cells.iter().map(|c| c.symbol))
            .collect();
        let distinct: std::collections::BTreeSet<&&str> = labels.iter().collect();
        assert_eq!(
            distinct.len(),
            labels.len(),
            "seed {seed}: two cells share a label, so the accessibility tree is ambiguous"
        );
    }
}

/// Every cell is on screen, after scrolling, at the size the window opens at.
///
/// **The only test here that says anything about the shipped program.** The
/// others use a harness sized so the whole table fits, which makes clipping a
/// non-subject; this one reads [`borbax_ui::WINDOW_SIZE`] — never a retyped
/// literal — because the defect it hunts is a function of the real window. The
/// widest measured table is ~3531 px of content against a 900 px window, so
/// most universes do not fit and a grid that does not scroll leaves the tail
/// permanently unreachable.
///
/// **This asserts geometry rather than a click, and the reason is a measured
/// harness limitation rather than a preference.** A pointer `click()` does not
/// reach widgets inside a *scrollable* `ScrollArea` under `egui_kittest`:
/// measured at 3600x900 (content fits) a cell selects, and at 900x600, 900x900
/// and 1400x600 — identical widget rect of `[200, 157]` in every case — it does
/// not. It is not input timing (a hover frame first changes nothing, and one to
/// four warm-up frames change nothing) and not drag-to-scroll (excluding
/// `ScrollSource::DRAG` changes nothing). The widget itself is correctly wired:
/// `click_accesskit()` selects it, and a pointer click on "surprise me" —
/// outside the scroll area, same harness size — works. So the app is fine and
/// the harness cannot express this one interaction.
///
/// Rect containment is not a weaker substitute for a click here; it is the
/// actual claim. The defect is "the cell can never be brought on screen", and
/// that is geometry. What is genuinely **not** covered is that a click at those
/// coordinates activates the cell — stated rather than implied, and covered at
/// a size where the harness can express it by
/// `clicking_a_cell_marks_that_cell_and_only_that_cell_as_chosen`.
#[test]
fn every_cell_is_on_screen_after_scrolling_at_the_size_the_window_opens_at() {
    let size = borbax_ui::WINDOW_SIZE;
    let mut harness = Harness::builder()
        .with_size(egui::vec2(size[0], size[1]))
        .build_ui_state(|ui, state| draw(state, ui), ViewerState::new());
    type_into(&harness, "seed", "7");
    harness.step();

    let symbols = widest_row_symbols(7);
    assert!(symbols.len() > 20, "seed 7's widest row should be long");

    // Picked from `rows()` rather than hard-coded: a literal symbol pins one
    // seed's table and stops tracking it.
    let last = symbols.last().cloned().unwrap_or_default();

    let before = harness.get_by_role_and_label(Role::Button, &last).rect();
    assert!(
        before.min.x > size[0],
        "the fixture assumes the last cell starts off screen; it is at {before:?} in a          {}-point window, so this test would pass without any scrolling at all",
        size[0]
    );

    harness
        .get_by_role_and_label(Role::Button, &last)
        .scroll_to_me();
    // **Three steps, measured.** A scroll takes effect over frames: the cell
    // sits at x = 2456 after one and two steps and reaches x = 852 on the third.
    // A fixed, explained count rather than a `run()`-shaped call, which would
    // step until quiescent and reintroduce the vacuity the rest of this file
    // forbids.
    for _ in 0..3 {
        harness.step();
    }

    let after = harness.get_by_role_and_label(Role::Button, &last).rect();
    assert!(
        after.min.x >= 0.0 && after.max.x <= size[0],
        "after scrolling, the last cell of the widest row is at {after:?}, which is          outside a {}-point window — the grid does not scroll, so the tail of every          wide table is unreachable",
        size[0]
    );
}

/// The grid is operable from the keyboard.
///
/// **A separate concern from reachability, and the one test here that
/// legitimately uses `focus()`.** Reachability asks whether a cell can be
/// brought under the pointer; this asks whether it can be reached at all without
/// one. `focus()` bypasses geometry, which is exactly why it must not be used
/// for the other question — and exactly why it is right for this one.
///
/// It matters more than politeness at 74 columns: a row that wide is unpleasant
/// to navigate by mouse, and the accessibility relations this needs are the same
/// ones every other test in this file addresses widgets through.
#[test]
fn the_grid_is_operable_from_the_keyboard() {
    let mut harness = wide_harness();
    type_into(&harness, "seed", "7");
    harness.step();

    let (only, _, _) = three_symbols(7);
    harness.get_by_role_and_label(Role::Button, &only).focus();
    harness.step();
    assert!(
        harness
            .get_by_role_and_label(Role::Button, &only)
            .is_focused(),
        "a cell cannot take keyboard focus, so the grid is mouse-only"
    );

    harness.key_combination(&[egui::Key::Enter]);
    harness.step();
    assert_eq!(
        harness
            .get_by_role_and_label(Role::Button, &only)
            .accesskit_node()
            .toggled(),
        Some(egui::accesskit::Toggled::True),
        "Enter on a focused cell does not select it"
    );
}
