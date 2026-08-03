//! The outside-in tests for the viewer: a seed or a name goes in, that
//! universe's element count comes out.
//!
//! Step 1 is the seed half; Step 1b added the name half below it.
//!
//! These are integration tests, so they compile the library **without**
//! `cfg(test)` — which is why `ViewerState::regenerations` is a shipped field
//! rather than a test-gated one. A `#[cfg(test)]` field would be invisible here.

use borbax_ui::state::{Outcome, ViewerState, parse_seed};
use borbax_universe::{PhysicsVersion, Universe};

/// Commit `seed` to a fresh state and hand back the line the window paints.
fn status_after(seed: &str) -> String {
    let mut state = ViewerState::new();
    seed.clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    state.status_line()
}

/// The expected line for `seed`, derived rather than pinned.
///
/// **Computed, never a literal.** A pinned `"87 elements · physics v1"` here
/// would be a second copy of `the_assembled_universe_digest_is_pinned`, and it
/// would make a legitimate change to the element-drawing physics fail in a
/// crate that has no opinion about physics and no business having one.
fn expected_for(seed: u64) -> String {
    let universe = Universe::generate(seed);
    format!(
        "{} elements · physics v{}",
        universe.table.len(),
        u8::from(universe.physics)
    )
}

#[test]
fn a_seed_typed_into_the_viewer_reports_that_universes_element_count() {
    assert_eq!(status_after("7"), expected_for(7));
}

/// The window opens on a universe, not on an instruction.
///
/// **Caught by looking at it rather than by a test**, which is the plan's own
/// acceptance bar working: the first build was correct on all 18 tests and
/// opened on the words "type a seed". A viewer whose first screen asks you to
/// do something has spent its opening move on a chore, and the whole argument
/// for building this before Task 11 was that there should be something to look
/// at.
///
/// `ViewerState::new` deliberately keeps the empty box — the frame tests type
/// into it, and a pre-filled box would make `type_text` append rather than
/// enter.
#[test]
fn the_viewer_opens_on_a_universe_rather_than_an_instruction() {
    let opened = ViewerState::opening();
    assert_eq!(opened.status_line(), expected_for(borbax_ui::OPENING_SEED));
    assert_eq!(
        opened.regenerations(),
        1,
        "opening the window should generate exactly one universe"
    );
}

/// **This test is vacuous today and is kept anyway. Read the reason before
/// trusting it.**
///
/// `status_line` is supposed to read the physics version out of the universe
/// (`u8::from(universe.physics)`) rather than type `v1` into the format string,
/// because §6 says a universe is `(seed, physics)` and a viewer that hardcodes
/// the version becomes a liar on the day `V2` ships — silently, on a screen
/// somebody has been reading for a year.
///
/// **Measured: substituting the literal `1` for `u8::from(universe.physics)`
/// passes every test in this crate, this one included.** The discriminator
/// asked for — "change `PhysicsVersion::CURRENT` and the line must move" — is
/// not available, because `PhysicsVersion` has exactly one variant. Every test
/// that could be written today compares 1 against 1.
///
/// **The count is deliberately gone.** It was re-measured three times as the
/// suite grew (16, 18, 24) under a sentence promising it always would be — and
/// then Step 1b took this crate to 36 and did not, in the one place that had
/// made the promise. A reviewer caught it. The mutation is the durable claim;
/// the denominator is a fact about how many tests happen to exist this week.
///
/// So what this is: a **latent** guard. It derives its expectation from
/// `PhysicsVersion::CURRENT` rather than from a literal, so it stays vacuous
/// until a second variant exists and fires the moment one does. Writing it as
/// `assert!(line.contains("v1"))` would have been the same number of characters
/// and would never fire at all.
///
/// The mechanism has a second, stronger protection that does not depend on this
/// test: `impl From<PhysicsVersion> for u8` in `borbax-universe` is an
/// exhaustive match over a `#[non_exhaustive]` enum, so adding `V2` fails to
/// compile there and forces the conversation. That protection only covers this
/// crate because the call is *made*; a literal would sail past it.
#[test]
fn the_physics_version_on_screen_is_read_from_the_universe() {
    let line = status_after("7");
    let current = u8::from(PhysicsVersion::CURRENT);
    assert!(
        line.ends_with(&format!("physics v{current}")),
        "expected the line to end with the current physics version {current}, got {line:?}"
    );
}

#[test]
fn the_random_button_puts_a_seed_in_the_box_that_can_be_typed_back() {
    // **The round-trip is the whole requirement, not a nicety.** §6 makes the
    // seed the shareable thing: a child who presses "surprise me", gets a
    // universe she likes and cannot write down what it was has lost exactly
    // what the seed is for. So the button does not secretly hold a number — it
    // *types into the box*, and what it typed must produce the universe on
    // screen when typed again by hand.
    let mut state = ViewerState::new();
    state.randomise(12_345);
    let shown = state.status_line();
    let text = state.seed_text().to_owned();

    assert!(
        !text.is_empty(),
        "the random button left the seed box empty"
    );
    assert_eq!(
        status_after(&text),
        shown,
        "retyping the seed the button produced ({text}) gives a different universe"
    );
}

#[test]
fn a_random_seed_is_short_enough_for_a_child_to_write_down() {
    // A full u64 is twenty digits. Nobody copies that off a screen, and a seed
    // nobody can copy is a seed §6's sharing promise does not reach. Typing a
    // long one by hand still works — this bounds only what the button offers.
    for moment in 0..2_000_u64 {
        let mut state = ViewerState::new();
        state.randomise(moment);
        let text = state.seed_text().to_owned();
        assert!(
            text.len() <= 6 && text.chars().all(|c| c.is_ascii_digit()),
            "the button produced {text:?}, which is not a short whole number"
        );
    }
}

#[test]
fn the_random_button_does_not_keep_giving_the_same_universe() {
    // Catches the stub (`randomise` ignoring its argument), a constant, and a
    // seed derived from something that does not move.
    let lines: std::collections::BTreeSet<String> = (0..500_u64)
        .map(|moment| {
            let mut state = ViewerState::new();
            state.randomise(moment);
            state.seed_text().to_owned()
        })
        .collect();
    assert!(
        lines.len() > 400,
        "500 presses produced only {} distinct seeds",
        lines.len()
    );
}

#[test]
fn the_same_moment_always_gives_the_same_seed() {
    // `randomise` is a pure function of the moment handed to it — the clock is
    // read at the edge, in the panel, and never inside this. That is what makes
    // every test above possible without a wall-clock in the assertions.
    let mut first = ViewerState::new();
    let mut second = ViewerState::new();
    first.randomise(99);
    second.randomise(99);
    assert_eq!(first.seed_text(), second.seed_text());
}

/// **Pins a §5 decision; has no power over the value it is pinning.**
///
/// An equality assertion cannot judge a title — it can only notice the second
/// one. So the rule the title was chosen under is written down where the next
/// person changing it will read it: **the title may name the program; it may
/// not name the program's subject.** "Borbax" is an invented proper noun and
/// asserts nothing. "Chemistry Simulator", "Molecular Explorer" and "Atomic
/// Sandbox" each name a real-world subject matter, which is the claim G6 (§5)
/// forbids — and each is the sort of thing that gets typed in without a thought
/// while making the window look finished.
///
/// A failure here is not a bug. It is a §5 conversation this test exists to
/// force.
#[test]
fn the_window_title_names_the_program_and_not_its_subject() {
    assert_eq!(borbax_ui::WINDOW_TITLE, "Borbax");
}

#[test]
fn a_refused_seed_leaves_no_number_claiming_to_be_its_universe() {
    // The defect this exists for is `parse(text).unwrap_or(0)` and its cousins
    // `unwrap_or(last_good_seed)` and "leave the old label painted" — each of
    // which puts a real, plausible element count on screen attributed to a seed
    // that did not produce it. That is the failure a human demonstrating the
    // app cannot see, because the number looks entirely fine.
    let good = expected_for(7);

    let mut state = ViewerState::new();
    "7".clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    assert_eq!(state.status_line(), good);

    "nonsense".clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    let refused = state.status_line();
    assert!(
        matches!(state.outcome(), Outcome::Rejected(_)),
        "a seed of `nonsense` left the viewer holding a universe: {refused}"
    );
    assert!(
        !refused.contains("elements"),
        "a refused seed still shows an element count: {refused}"
    );

    // And the good seed comes back exactly — no accumulated state.
    "7".clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    assert_eq!(state.status_line(), good);
}

#[test]
fn the_same_seed_gives_the_same_line_however_the_viewer_got_there() {
    // Path-independence rather than "`generate(7)` twice is equal", which is
    // vacuous — `Universe::generate` is a pure function and
    // `universes_are_deterministic` in `borbax-universe` already owns that
    // assertion. What this catches is wall-clock or `RandomState`
    // contamination in the viewer, and any state that accumulates across
    // commits.
    let fresh = status_after("7");

    let mut travelled = ViewerState::new();
    for seed in [
        "11",
        "",
        "-3",
        "not a seed",
        "99999999999999999999999",
        "11",
    ] {
        seed.clone_into(travelled.seed_text_mut());
        travelled.commit_typed_seed();
    }
    "7".clone_into(travelled.seed_text_mut());
    travelled.commit_typed_seed();

    assert_eq!(travelled.status_line(), fresh);
}

#[test]
fn two_universes_do_not_have_to_agree_on_how_many_elements_they_have() {
    // The bar is `> 1` **because its job is to catch a constant** — a stub
    // returning 100, a universe built once in `new()` and never rebuilt, a
    // `OnceLock` cache. It is deliberately NOT raised to the number actually
    // observed: pinning the distribution would make a legitimate physics change
    // fail here, and this crate must never acquire an opinion about physics.
    let counts: std::collections::BTreeSet<String> = (0..64)
        .map(|seed| status_after(&seed.to_string()))
        .collect();

    assert!(
        counts.len() > 1,
        "64 seeds produced one line, so the count is not coming from the seed"
    );
}

#[test]
fn the_reported_count_is_this_universes_table_length() {
    // Tautological against the *first* implementation — it re-derives the
    // expected value with the same call the implementation makes, so it cannot
    // fail for any honest first attempt. It is kept for what it catches later:
    // the label being re-sourced from `pattern()`, `MAX_ELEMENTS`, a filtered
    // `iter().count()`, or an off-by-one.
    //
    // It earns its place over a `60..=120` range check, which was considered
    // and dropped: a range test catches `len() + 1` only when some seed draws
    // exactly 120, which over 64 seeds is 1 - (60/61)^64 ~= 65% — a test that
    // finds a real defect two times in three is a flaky test, and flaky tests
    // get deleted taking the guarantee with them.
    for seed in 0..64_u64 {
        assert_eq!(status_after(&seed.to_string()), expected_for(seed));
    }
}

// ---------------------------------------------------------------------------
// Step 1b — a name instead of a number.
// ---------------------------------------------------------------------------

/// The whole point of the step: a child types her name and gets a universe.
///
/// **The name is a key, not a label** — Ian's decision, stated in his framing as
/// "Emily" giving Emily's universe on any machine, forever. The alternative put
/// to him was a random seed with a text box in front of it, which makes the name
/// decoration: that is the "surprise me" button, and it already exists.
#[test]
fn a_name_typed_into_the_viewer_opens_that_names_universe() {
    let mut state = ViewerState::new();
    state.set_phrase("Emily");

    let seed = borbax_universe::seed_from_phrase("Emily")
        .unwrap_or_else(|| unreachable!("`Emily` is not blank, so it names a universe"));
    assert_eq!(state.status_line(), expected_for(seed));
    assert_eq!(
        state.seed_text(),
        seed.to_string(),
        "the seed box shows the number the name produced, so it can be written down"
    );
}

/// A name and a seed cannot disagree about which universe is on screen.
///
/// **The invariant, and it is what the two-box design buys.** Whenever the name
/// box holds something that names a universe, the seed box holds that
/// universe's number — so the number can be written down, and the name never
/// labels a universe it did not produce.
/// **Both arms assert, and the `None` arm is a repair.** The helper used to be
/// a bare `if let Some(..)`, which meant it checked nothing at two of its three
/// call sites — both of them clear the name box on the line above, so
/// `seed_from_phrase("")` is `None` and the whole body was skipped, by
/// construction, every run, forever. A review lane measured that with an
/// `eprintln!`. It is the "corpus filter that silently selects nothing" class
/// CLAUDE.md names, and the fix is for the helper to have nothing to skip.
fn the_boxes_agree(state: &ViewerState) {
    match borbax_universe::seed_from_phrase(state.phrase_text()) {
        Some(seed) => {
            assert_eq!(
                state.seed_text(),
                seed.to_string(),
                "the name box says {:?} but the seed box says {:?}",
                state.phrase_text(),
                state.seed_text()
            );
            assert_eq!(state.status_line(), expected_for(seed));
        }
        // Nothing is named, so the seed box stands on its own and the line on
        // screen must describe whatever it holds — including a refusal.
        None => match parse_seed(state.seed_text()) {
            Ok(seed) => assert_eq!(state.status_line(), expected_for(seed)),
            Err(refusal) => assert_eq!(state.status_line(), refusal.to_string()),
        },
    }
}

/// A digit typed into the *name* box is a name, not a seed.
///
/// **This is the two-box decision made checkable, and it is the test that would
/// be impossible with one box.** A single control has to guess whether `"7"`
/// means seed 7 or a child called 7, and guessing means deleting
/// `letters_in_the_seed_box_are_refused_not_hashed` — the test that forbids "be
/// friendly, hash whatever was typed". With two boxes each control keeps its
/// own unambiguous job and `parse_seed` keeps all four of its refusals.
#[test]
fn a_digit_in_the_name_box_names_a_universe_rather_than_selecting_one() {
    let mut state = ViewerState::new();
    state.set_phrase("7");

    let named = borbax_universe::seed_from_phrase("7")
        .unwrap_or_else(|| unreachable!("`7` is not blank, so it names a universe"));
    assert_ne!(named, 7, "the name box fell through to seed parsing");
    assert_eq!(state.status_line(), expected_for(named));
    the_boxes_agree(&state);
}

/// Typing a seed by hand clears the name, because that name did not produce it.
///
/// Without this the window shows `Emily` beside a universe Emily did not name —
/// a plausible number attributed to the wrong thing, which is the class
/// `Outcome` was made an enum to prevent one field along.
#[test]
fn typing_a_seed_by_hand_stops_claiming_the_old_name_produced_it() {
    let mut state = ViewerState::new();
    state.set_phrase("Emily");
    assert!(!state.phrase_text().is_empty());

    "42".clone_into(state.seed_text_mut());
    state.commit_typed_seed();

    assert_eq!(
        state.phrase_text(),
        "",
        "the name outlived the seed it named"
    );
    assert_eq!(state.status_line(), expected_for(42));
    the_boxes_agree(&state);
}

/// "surprise me" clears the name too, and for the same reason.
#[test]
fn a_surprise_seed_stops_claiming_the_old_name_produced_it() {
    let mut state = ViewerState::new();
    state.set_phrase("Emily");
    state.randomise(1_234_567_890);

    assert_eq!(state.phrase_text(), "");
    the_boxes_agree(&state);
}

/// Emptying the name box leaves the universe alone rather than blanking it.
///
/// The alternative — clearing the seed box — was rejected because it
/// contradicts the state the window opens in: name box empty, a universe on
/// screen. Two rules for the same visible input is the wart; this is the
/// version where an empty name box means "not currently driving", at startup
/// and after a deletion alike.
#[test]
fn emptying_the_name_box_does_not_take_the_universe_with_it() {
    let mut state = ViewerState::new();
    state.set_phrase("Emily");
    let named = state.status_line();

    state.set_phrase("");
    assert_eq!(state.status_line(), named);
    // And blank-but-not-empty behaves the same, since it names nothing either.
    state.set_phrase("   ");
    assert_eq!(state.status_line(), named);
}

/// The same name gives the same universe however the viewer got there.
#[test]
fn a_name_is_a_key_and_not_a_label() {
    let mut typed = ViewerState::new();
    typed.set_phrase("Emily");

    let mut differently = ViewerState::new();
    differently.set_phrase("  eMiLy ");

    assert_eq!(typed.status_line(), differently.status_line());
    assert_eq!(typed.seed_text(), differently.seed_text());
}

/// Naming a universe generates exactly one, not one per keystroke path.
#[test]
fn naming_a_universe_generates_exactly_one() {
    let mut state = ViewerState::new();
    state.set_phrase("Emily");
    assert_eq!(state.regenerations(), 1);
    // A name that names nothing must not generate at all.
    state.set_phrase("  ");
    assert_eq!(state.regenerations(), 1);
}

// ---------------------------------------------------------------------------
// Step 2 — the periodic table panel.
// ---------------------------------------------------------------------------

/// A state showing `seed`'s universe.
fn loaded(seed: &str) -> ViewerState {
    let mut state = ViewerState::new();
    seed.clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    state
}

/// Every element id of a loaded state, asserted non-empty.
///
/// **The floor is the point.** If `loaded(seed)` ever returned a `Rejected`
/// outcome — a refused seed, a changed parse — `rows()` gives an empty vector,
/// the `for id in ids` loop below never runs, and the test passes having
/// asserted nothing at all. Four tests had that shape. It is the same vacuity
/// the shape-count floor further up exists for, and the same class this
/// repository has now found in a guard six times.
fn all_ids(state: &ViewerState) -> Vec<borbax_universe::ElementId> {
    let ids: Vec<_> = state
        .rows()
        .iter()
        .flat_map(|r| r.cells.iter().map(|c| c.id))
        .collect();
    assert!(
        ids.len() >= 60,
        "a loaded universe has 60..=120 elements; {} means the state was not loaded \
         and every assertion over these ids would pass vacuously",
        ids.len()
    );
    ids
}

/// The id of a loaded state's first cell.
///
/// Indexing is denied workspace-wide, and a panicking index in a test is still a
/// panicking index — it replaces a clear assertion failure with a slice
/// out-of-bounds carrying no message.
fn first_cell(state: &ViewerState) -> borbax_universe::ElementId {
    let id = state
        .rows()
        .first()
        .and_then(|row| row.cells.first())
        .map(|cell| cell.id);
    assert!(
        id.is_some(),
        "a loaded universe always has at least one cell"
    );
    id.unwrap_or(borbax_universe::ElementId::ZERO)
}

/// Every element of `seed`'s universe reaches exactly one cell, and no cell
/// holds anything the table does not.
///
/// **Both directions, and the reverse arm is the one that earns its place.**
/// Forward alone is satisfied by a grid that also carries padding cells — the
/// `ElementId::ZERO` filler somebody writes to make rows rectangular, which puts
/// a real element's symbol in a slot no element occupies. It also kills
/// `row = i / 18, col = i % 18` and every other wrap-based layout that looks
/// like a table from a distance.
#[test]
fn every_element_sits_at_its_own_period_and_group_and_no_cell_holds_two() {
    for seed in ["1", "7", "42", "15709401653729972761"] {
        let state = loaded(seed);
        let universe = Universe::generate(parse_seed(seed).unwrap_or_default());

        let placed: Vec<_> = state
            .rows()
            .iter()
            .flat_map(|row| row.cells.iter().map(|c| c.id))
            .collect();
        let mut sorted = placed.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            placed.len(),
            "seed {seed}: a cell was painted twice"
        );
        assert_eq!(
            placed.len(),
            universe.table.len(),
            "seed {seed}: the grid holds a different number of cells than the table has \
             elements — either an element is missing or a padding cell was invented"
        );
        for id in &placed {
            assert!(
                universe.table.get(*id).is_some(),
                "seed {seed}: the grid holds {id:?}, which this table has no slot for"
            );
        }
    }
}

/// The column a cell sits in is its group, and the row is its period.
///
/// **This pins a property of `generate_elements`, not of `PeriodicTable`.**
/// `rows()` walks the table in one forward pass and starts a new row when
/// `period` changes — which is only a correct layout because `group` runs
/// 0,1,2,… contiguously inside each period. Verified over 5000 universes with
/// zero exceptions, but it is nowhere documented as a contract, so if it ever
/// stops holding the grid silently transposes elements rather than failing.
#[test]
fn the_column_a_cell_sits_in_is_its_group() {
    for seed in 1..200_u64 {
        let state = loaded(&seed.to_string());
        let universe = Universe::generate(seed);
        for (r, row) in state.rows().iter().enumerate() {
            for (c, cell) in row.cells.iter().enumerate() {
                let element = universe.table.get(cell.id);
                assert!(element.is_some(), "seed {seed}: cell {c} has no element");
                let Some(element) = element else { continue };
                assert_eq!(
                    usize::from(element.group),
                    c,
                    "seed {seed}: {} sits in column {c} but its group is {}",
                    element.symbol,
                    element.group
                );
                assert_eq!(
                    usize::from(element.period),
                    r,
                    "seed {seed}: {} sits in row {r} but its period is {}",
                    element.symbol,
                    element.period
                );
            }
        }
    }
}

/// Two universes do not have to agree on the shape of the table.
///
/// Kills a shape computed once — a `OnceLock`, a `static`, or a grid built in
/// `new()` and never rebuilt by `reload`.
///
/// **Asserted on the `(rows, cols)` pair, never on rows alone.** Row counts
/// collapse to {2, 3, 4} across every universe, so "several seeds all gave the
/// same number of rows" is a 1-in-18 000 coincidence rather than a defect — a
/// test asserting on rows alone would be a flake with a plausible story.
#[test]
fn two_universes_do_not_have_to_agree_on_the_shape_of_the_table() {
    let shapes: std::collections::BTreeSet<(usize, usize)> = (1..64_u64)
        .map(|seed| {
            let state = loaded(&seed.to_string());
            let rows = state.rows();
            let widest = rows.iter().map(|r| r.cells.len()).max().unwrap_or_default();
            (rows.len(), widest)
        })
        .collect();
    assert!(
        shapes.len() > 1,
        "every seed produced the same table shape {shapes:?} — the grid is not being \
         rebuilt from the universe"
    );
}

/// Nothing is selected until something is selected.
///
/// Kills `selected: Some(ElementId::ZERO)` as a default, which opens the window
/// painting element 0's properties under a heading nobody clicked — the same
/// class as `an_empty_seed_box_does_not_quietly_mean_universe_zero` one field
/// along. Every other test here selects first, so only this one sees the
/// opening state.
#[test]
fn nothing_is_selected_until_something_is_selected() {
    let state = loaded("7");
    assert_eq!(state.selected(), None);
    assert_eq!(state.selection_heading(), "pick an element");
    assert!(
        state
            .rows()
            .iter()
            .all(|r| r.cells.iter().all(|c| !c.selected)),
        "a cell was drawn selected before anything was clicked"
    );
    // The empty block still paints all nine labels, so the grid above it does
    // not reflow on the first click.
    let rows = state.selection_properties();
    assert_eq!(rows.len(), 9);
    assert!(
        rows.iter().all(|r| r.unit.is_empty() && r.gloss.is_none()),
        "the empty state must carry no unit and no gloss: {rows:?}"
    );
}

/// A selection cannot outlive the universe it was made in.
///
/// **Two mutations, and the first is what someone writes the first hour they
/// want "keep the same element when I change the seed".** Lifting `selected` to
/// a `ViewerState` field compiles, reads as a tidy-up, and silently paints slot
/// 40 of the *new* universe under the old element's heading — both are valid
/// ids, so nothing downstream can detect the substitution. The second is
/// `Rejected { last_selected }`, the same instinct one variant over.
///
/// The property is currently held by the type rather than by this test, which
/// is the point: the test exists so that lifting the field out is a failure
/// rather than a refactor.
#[test]
fn a_selection_cannot_outlive_the_universe_it_was_made_in() {
    let mut state = loaded("7");
    let first = first_cell(&state);
    state.select(first);
    assert_eq!(state.selected(), Some(first));

    // Reload into a different, valid universe.
    "11".clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    assert_eq!(
        state.selected(),
        None,
        "the selection survived a change of universe"
    );

    // And reload into a refused seed.
    let next = first_cell(&state);
    state.select(next);
    assert!(state.selected().is_some());
    "not a seed".clone_into(state.seed_text_mut());
    state.commit_typed_seed();
    assert_eq!(
        state.selected(),
        None,
        "a selection survived into a state with no universe at all"
    );
}

/// An id this universe has no slot for selects nothing.
///
/// Kills `table.get(id).unwrap()` and its quieter cousin, clamping to the last
/// slot — which is in range, plausible, and undetectable downstream, which is
/// exactly what `PeriodicTable::get` is fallible to prevent.
#[test]
fn an_element_id_this_universe_has_no_slot_for_selects_nothing() {
    let mut state = loaded("7");
    let count = state.rows().iter().map(|r| r.cells.len()).sum::<usize>();
    let past_the_end = borbax_universe::ElementId::from_index(count + 1);
    assert!(
        past_the_end.is_some(),
        "the fixture needs a constructible id"
    );
    if let Some(id) = past_the_end {
        state.select(id);
    }
    assert_eq!(
        state.selected(),
        None,
        "an out-of-range id was stored, so the properties block will read some other \
         element's numbers"
    );
}

/// The properties shown are the selected element's and not a neighbour's.
///
/// Kills an off-by-one at the read. Derived over **every** element of several
/// universes rather than one, because an off-by-one at index 0 is invisible.
#[test]
fn the_properties_shown_are_the_selected_elements_and_not_a_neighbours() {
    for seed in ["1", "7", "42"] {
        let mut state = loaded(seed);
        let universe = Universe::generate(parse_seed(seed).unwrap_or_default());
        for id in all_ids(&state) {
            state.select(id);
            let element = universe.table.get(id);
            assert!(element.is_some());
            let Some(element) = element else { continue };
            assert_eq!(
                state.selection_heading(),
                format!("{} · {}", element.symbol, element.name),
                "seed {seed}: the heading names the wrong element"
            );
            // Looked up by label rather than by position: an index here would
            // pin the row order, which is a layout decision this test has no
            // opinion about, and it would panic rather than assert if the order
            // changed.
            let rows = state.selection_properties();
            let value_of = |label: &str| {
                rows.iter()
                    .find(|r| r.label == label)
                    .map(|r| r.value.clone())
            };
            assert_eq!(
                value_of("made of"),
                Some(element.units.to_string()),
                "seed {seed}: `made of` came from the wrong element"
            );
            assert_eq!(
                value_of("outer shell"),
                Some(element.group.to_string()),
                "seed {seed}: `outer shell` came from the wrong element"
            );
            // **All nine rows, not the two that were easy.** `selection_properties`
            // pairs `property(i, ..)` with `PROPERTY_LABELS[i]` *by position*, so
            // a reorder of either list silently re-labels every row — and a
            // mutation sweep found exactly that class already: deleting the
            // shell label's `+ 1`, and formatting mass through `Debug`, each
            // failed zero tests. Asserting two of nine left seven undefended.
            let shells = state.rows().len();
            assert_eq!(
                value_of("shell"),
                Some(format!("{} of {}", element.period + 1, shells)),
                "seed {seed}: `shell` came from the wrong element"
            );
            assert_eq!(
                value_of("bonding slots"),
                Some(if element.valence == 0 {
                    "none \u{2014} a closed shell".to_owned()
                } else {
                    element.valence.to_string()
                }),
                "seed {seed}: `bonding slots` came from the wrong element"
            );
            assert_eq!(
                value_of("size"),
                Some(format!("{:.3}", element.radius.get())),
                "seed {seed}: `size` came from the wrong element"
            );
            assert_eq!(
                value_of("surface"),
                Some(format!("{:.3}", element.affinity)),
                "seed {seed}: `surface` came from the wrong element"
            );
            assert_eq!(
                value_of("binding energy"),
                Some(format!("{:.3}", element.energy_per_unit.get())),
                "seed {seed}: `binding energy` came from the wrong element"
            );
            assert_eq!(
                value_of("instability"),
                Some(format!("{:.3}", element.instability)),
                "seed {seed}: `instability` came from the wrong element"
            );
        }
    }
}

/// Every unit word on screen is one this universe invented.
///
/// **This is the fail-on-unknown half of G4 and the only one there can be.** A
/// source deny-list over real units cannot invert — there is no finite set of
/// them — so `xtask`'s scan catches the spellings anyone writes and nothing
/// else. This is an allow-list over the *attained* set, so it also catches a
/// real unit assembled at runtime, which no source scan can see, and a
/// *missing* unit word, which is a bare number attributed to nothing.
///
/// **Asserted on the attained set, never against the declared constant** — that
/// would be a guard whose two sides are the same expression.
#[test]
fn every_unit_word_on_screen_is_one_this_universe_invented() {
    const INVENTED: [&str; 7] = [
        "",
        "base units",
        "units",
        "mass units",
        "spans",
        "quanta per unit",
        "bonding slots",
    ];
    let mut attained = std::collections::BTreeSet::new();
    for seed in 1..40_u64 {
        let mut state = loaded(&seed.to_string());
        for id in all_ids(&state) {
            state.select(id);
            for row in state.selection_properties() {
                attained.insert(row.unit);
            }
        }
    }
    for unit in &attained {
        assert!(
            INVENTED.contains(unit),
            "the word {unit:?} reached the screen and this universe did not invent it"
        );
    }
    // Negative arms: no Step 2 property carries a temperature or a time, so
    // either arriving means a property acquired a dimension it does not have —
    // and `world-years` is also the forbidden phrasing for `instability`.
    assert!(!attained.contains("thermals"));
    assert!(!attained.contains("world-years"));
    // And the block must not be all-dimensionless, which would pass the loop
    // above while showing bare numbers attributed to nothing.
    assert!(
        attained.len() > 1,
        "no property carried a unit word at all: {attained:?}"
    );
}

/// No element reaches the screen under a real name or symbol.
///
/// **Neither `xtask` guard can reach this**, and that is why it is here: a
/// display-side *derivation* has no literal to scan. `&name[..2]` for a narrow
/// cell turns a generated `helion` into `He`; `name.to_uppercase()` and
/// `symbol.chars().next()` are the same move. The mint checks `is_real` once,
/// at generation, and never again on anything the viewer computes from it.
#[test]
fn no_element_reaches_the_screen_under_a_real_name_or_symbol() {
    for seed in 1..120_u64 {
        let mut state = loaded(&seed.to_string());
        let cells: Vec<(borbax_universe::ElementId, String)> = state
            .rows()
            .iter()
            .flat_map(|r| r.cells.iter().map(|c| (c.id, c.symbol.to_owned())))
            .collect();
        assert!(
            cells.len() >= 60,
            "seed {seed}: {} cells means the state was not loaded and every assertion \
             below would pass vacuously",
            cells.len()
        );
        for (id, symbol) in cells {
            state.select(id);
            let heading = state.selection_heading();
            let name = heading.split(" · ").nth(1).unwrap_or_default().to_owned();
            assert!(
                !borbax_universe::naming::is_real(&symbol, &name),
                "seed {seed}: the screen shows {symbol}/{name}, which collides with \
                 something real"
            );
        }
    }
}

/// Selecting an element does not regenerate the universe.
///
/// Kills `select()` calling `reload()` — the plausible "refresh the view"
/// implementation, which costs a generation per click and would also reset the
/// seed box.
#[test]
fn selecting_an_element_does_not_regenerate_the_universe() {
    let mut state = loaded("7");
    let before = state.regenerations();
    for id in all_ids(&state) {
        state.select(id);
    }
    assert_eq!(
        state.regenerations(),
        before,
        "selecting rebuilt the universe"
    );
}

/// The same universe gives the same table however the viewer got there.
///
/// Kills a grid pushed to rather than rebuilt — cells accumulating across seed
/// changes.
#[test]
fn the_same_universe_gives_the_same_table_however_the_viewer_got_there() {
    let direct = loaded("7");
    let mut wandered = loaded("11");
    "42".clone_into(wandered.seed_text_mut());
    wandered.commit_typed_seed();
    "7".clone_into(wandered.seed_text_mut());
    wandered.commit_typed_seed();

    let shape = |s: &ViewerState| -> Vec<Vec<String>> {
        s.rows()
            .iter()
            .map(|r| r.cells.iter().map(|c| c.symbol.to_owned()).collect())
            .collect()
    };
    assert_eq!(shape(&direct), shape(&wandered));
}

/// No explanatory line borrows a word from real chemistry.
///
/// **The only thing standing between "surface, unlike electronegativity…" and a
/// child's screen.** These words are not element names and not units, so both
/// `xtask` guards are blind to them; and they cannot be added to the generator's
/// blocklist, which would change *generation*. A closed set of `&'static str`
/// constants is enumerable at compile time, so this needs no scan.
#[test]
fn no_explanatory_line_borrows_a_word_from_real_chemistry() {
    const FORBIDDEN: [&str; 8] = [
        "electronegativity",
        "electron",
        "half-life",
        "halflife",
        "radioactive",
        "decays at",
        "per world-year",
        "atomic number",
    ];
    let mut state = loaded("7");
    let pick = first_cell(&state);
    state.select(pick);

    let mut prose: Vec<String> = state
        .selection_properties()
        .iter()
        .flat_map(|r| {
            [
                r.label.to_owned(),
                r.unit.to_owned(),
                r.gloss.unwrap_or_default().to_owned(),
            ]
        })
        .collect();
    prose.push(state.selection_heading());

    for line in &prose {
        let lower = line.to_ascii_lowercase();
        for word in FORBIDDEN {
            assert!(
                !lower.contains(word),
                "the panel says {line:?}, which borrows {word:?} from real chemistry"
            );
        }
    }
    // The positive arm: the two glosses must actually be there, or this test
    // passes over a panel that explains nothing.
    let glossed = state
        .selection_properties()
        .iter()
        .filter(|r| r.gloss.is_some())
        .count();
    assert_eq!(glossed, 2, "the two explanatory lines are missing");
}

/// Rows are numbered from one, not from zero.
///
/// **Found by mutation probe, not by design: deleting the `+ 1` failed zero
/// tests.** `Element::period` is 0-based and the label is 1-based, so the
/// conversion is a display decision with nothing else defending it — and "shell
/// 0" is the single number a child reads wrong. The doc comment claiming this
/// was already there; the test was not, which is the difference between a
/// decision and a decoration.
#[test]
fn the_shell_label_counts_from_one_not_zero() {
    for seed in ["1", "7", "42"] {
        let state = loaded(seed);
        let universe = Universe::generate(parse_seed(seed).unwrap_or_default());
        for row in state.rows() {
            let period = row
                .cells
                .first()
                .and_then(|c| universe.table.get(c.id))
                .map(|e| e.period);
            assert!(period.is_some(), "seed {seed}: an empty row was built");
            let Some(period) = period else { continue };
            assert_eq!(
                row.shell_label,
                format!("shell {}", u16::from(period) + 1),
                "seed {seed}: the row label does not count from one"
            );
        }
        // The strongest arm, and the one the `+ 1` mutation fails: the first row
        // is always period 0 and must never read "shell 0".
        let first = state.rows().first().map(|r| r.shell_label.clone());
        assert_eq!(first, Some("shell 1".to_owned()), "seed {seed}");
    }
}

/// Mass is shown in mass units, not in raw fixed-point sub-units.
///
/// **Found by mutation probe: `format!("{:?}", element.mass)` failed zero
/// tests.** `Mass` is fixed-point over `i64`, so its `Debug` prints
/// `Mass(1792)` where the value is `1.750` — wrong by a factor of 1024, and
/// plausible enough on screen to survive a glance. Three of the other four unit
/// types `Debug`-print something close enough to the right answer to pass
/// review, which is exactly why this one needs a test rather than a comment.
#[test]
fn mass_is_shown_in_mass_units_not_raw_sub_units() {
    for seed in ["1", "7", "42"] {
        let mut state = loaded(seed);
        let universe = Universe::generate(parse_seed(seed).unwrap_or_default());
        for id in all_ids(&state) {
            state.select(id);
            let element = universe.table.get(id);
            assert!(element.is_some());
            let Some(element) = element else { continue };
            let shown = state
                .selection_properties()
                .into_iter()
                .find(|r| r.label == "mass")
                .map(|r| r.value);
            assert_eq!(
                shown,
                Some(format!("{:.3}", element.mass.to_f64())),
                "seed {seed}: mass is not the value `to_f64` gives — a `Debug` \
                 spelling prints sub-units and is wrong by 1024x"
            );
        }
    }
}
