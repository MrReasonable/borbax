//! The outside-in test for Step 1 of the viewer plan: a seed goes in, that
//! universe's element count comes out.
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
    state.reload();
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
/// passes this entire test suite — all 24 tests, this one included.** The
/// discriminator asked for — "change `PhysicsVersion::CURRENT` and the line must
/// move" — is not available, because `PhysicsVersion` has exactly one variant.
/// Every test that could be written today compares 1 against 1. The figure has
/// been re-measured each time the suite grew (16, then 18, then 24) rather than
/// carried forward, because a number quoted from an earlier run is the one kind
/// of evidence that rots without anyone noticing.
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
    state.reload();
    assert_eq!(state.status_line(), good);

    "nonsense".clone_into(state.seed_text_mut());
    state.reload();
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
    state.reload();
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
        travelled.reload();
    }
    "7".clone_into(travelled.seed_text_mut());
    travelled.reload();

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
