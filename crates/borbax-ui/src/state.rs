//! What is selected, and every decision the viewer makes.
//!
//! **This file imports no UI type**, which is the whole reason it exists
//! separately: it is the half of the crate a test can reach without a window.

use core::num::IntErrorKind;
use std::time::{SystemTime, UNIX_EPOCH};

use borbax_rng::{Domain, Stream};
use borbax_universe::{Element, ElementId, Universe, seed_from_phrase};

use crate::molecule::Demo;
use thiserror::Error;

/// What a property row shows when nothing is selected.
const EMPTY_VALUE: &str = "\u{2014}";

/// One property of `element`, as `(value, unit)`, in [`PROPERTY_LABELS`] order.
///
/// **`{:.3}` fixed, never bare `{}`.** `Display` on an `f64` prints the shortest
/// round-tripping decimal, so widths jump between `0.3` and
/// `0.30000000000000004` from row to row. Measured ranges make three decimals
/// right for all of them: mass 1…208.78, radius 0.300…2.087, affinity
/// 0.140…1.000, energy/unit −0.278…5.710, instability 0…1.
///
/// **No `is_finite` guard**, deliberately. Measured: zero non-finite values
/// across 180 864 elements × 6 fields. A branch that cannot fire cannot be
/// tested, and an untestable branch is the vacuous-guard shape this repository
/// keeps re-discovering.
///
/// **`Mass::to_f64()`, never `{:?}`.** `Mass` is fixed-point over `i64`, so its
/// `Debug` prints raw 1/1024 sub-units — `Mass(1792)` for a mass of 1.75, wrong
/// by a factor of 1024. Three of the other four units `Debug`-print something
/// close enough to survive review, which is what makes this worth stating.
fn property(i: usize, element: &Element, shells: usize) -> (String, &'static str) {
    match i {
        // The element *is* this number — element N is N copies of one base
        // unit. Never "atomic number", never "Z", never "protons": there are no
        // protons in Borbax.
        0 => (element.units.to_string(), "base units"),
        1 => (format!("{} of {}", element.period + 1, shells), ""),
        // `group` is the count of units in the incomplete outer shell, not a
        // family index. The column position already shows the group; the number
        // worth painting is the count, with its meaning attached.
        2 => (element.group.to_string(), "units"),
        // Zero is a sentence and non-zero is a bare number, deliberately: the
        // closed shell is the interesting value and the one a child can be told
        // something about.
        3 if element.valence == 0 => ("none \u{2014} a closed shell".to_owned(), ""),
        3 => (element.valence.to_string(), "bonding slots"),
        4 => (format!("{:.3}", element.mass.to_f64()), "mass units"),
        5 => (format!("{:.3}", element.radius.get()), "spans"),
        // Dimensionless, and no bar is drawn. The documented range is [-1, +1]
        // but the attained range is [0.140476, 1.0], so a bar on the documented
        // scale would put every element in the top half and teach that this
        // universe has no "negative" elements when the scale is simply not
        // reached.
        6 => (format!("{:.3}", element.affinity), ""),
        7 => (
            format!("{:.3}", element.energy_per_unit.get()),
            "quanta per unit",
        ),
        // Dimensionless. NOT a rate, and the gloss says so in words that avoid
        // "half-life", "radioactive" and "per world-year".
        _ => (format!("{:.3}", element.instability), ""),
    }
}

/// One period of the table, as a row of cells in group order.
///
/// Borrows its strings from the universe rather than owning them: the whole
/// table's cells are rebuilt on every paint, and ~120 `String` allocations a
/// frame buys nothing when the borrow is free.
#[derive(Debug)]
pub struct PeriodRow<'a> {
    /// The row's label, already built: `shell 2`.
    ///
    /// **A `String` rather than the number, because `panel.rs` has no
    /// `format!`.** That is the seam doing its job rather than an inconvenience:
    /// a row label built where it is painted is a label no test can assert on
    /// without describing a window it has stopped describing.
    ///
    /// **1-based, and the field it comes from is 0-based.** A row labelled
    /// `shell 0` is the one number a child reads wrong, and the `+ 1` is a
    /// display decision, so it belongs on this side of the seam.
    pub shell_label: String,
    /// Which period this is, 0-based as the element carries it.
    ///
    /// Kept beside the label so the row break compares numbers rather than
    /// re-parsing the string it just built.
    pub period: u8,
    /// The elements of this period, in group order.
    pub cells: Vec<Cell<'a>>,
}

/// One element's cell in the grid.
#[derive(Debug)]
pub struct Cell<'a> {
    /// Which element this cell is.
    pub id: ElementId,
    /// The generated symbol — **not** the name.
    ///
    /// **Symbols are unique within a universe and names are not.** Measured
    /// over 2000 universes: zero symbol collisions (`mint` regenerates against
    /// a `taken` set that holds symbols only) against 357 duplicate names
    /// across 321 universes — 16.05%, the first at seed 16, where `meax` names
    /// both `Mx` and `Me`. A cell labelled by name would therefore be ambiguous
    /// in one universe in six, which makes the accessibility tree ambiguous and
    /// every frame test silently address the wrong cell.
    ///
    /// It is also the only thing that fits: names run 5–9 characters and the
    /// widest measured row holds 74 of them.
    pub symbol: &'a str,
    /// Whether this is the selected cell.
    pub selected: bool,
}

/// One line of the properties block.
///
/// **Four fields rather than a `(label, value)` pair, and the split of `value`
/// from `unit` is what makes G4 checkable.** With the unit word inside the
/// formatted value, an allow-list over what reaches the screen would have to
/// parse sentences, and a unit assembled at runtime would be indistinguishable
/// from a declared one. As its own `&'static str` the attained set is
/// enumerable, which is exactly what
/// `every_unit_word_on_screen_is_one_this_universe_invented` asserts over.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PropertyRow {
    /// What the row is called.
    pub label: &'static str,
    /// The value, formatted but carrying no unit.
    pub value: String,
    /// This universe's own unit word, or `""` when the value is dimensionless.
    pub unit: &'static str,
    /// A plain-English line under the row, for the two that need one.
    ///
    /// **Attached to its row rather than collected at the bottom.** A footnote
    /// adjacent to nothing is a footnote nobody reads, and a fixed footnote
    /// position would make the panel index into the row order — the coupling
    /// `labelled_by` removed from the frame tests at Step 1b.
    pub gloss: Option<&'static str>,
}

/// The labels of the properties block, in the order they are painted.
///
/// Named separately from the values so the empty state can paint the same rows
/// with nothing in them: a block that appears on first click reflows the grid,
/// and a grid that jumps when you click it is the worst interaction on this
/// screen.
const PROPERTY_LABELS: [&str; 9] = [
    "made of",
    "shell",
    "outer shell",
    "bonding slots",
    "mass",
    "size",
    "surface",
    "binding energy",
    "instability",
];

/// The two rows that need a sentence, and the sentences.
///
/// **Both avoid a real-chemistry word, deliberately and by name.** `surface`
/// must never be glossed as electronegativity — G3 says `affinity` has a
/// different range, a different meaning and different behaviour, and the
/// measured range is [0.140476, 1.0], so it does not even reach the negative
/// half the documented scale describes. `instability` must never be glossed as
/// a half-life, a decay rate, or anything "per world-year": §7.1 calls it a
/// probability per world-year while Task 15 consumes it as an unbounded
/// Gillespie propensity, and a viewer that writes a rate settles an open
/// physics question in the crate least entitled to settle one.
///
/// Held by `no_explanatory_line_borrows_a_word_from_real_chemistry`, because
/// neither `xtask` guard can see these: they are not element names and not
/// units.
const SURFACE_GLOSS: &str = "how much of this element sits on the outside";
const INSTABILITY_GLOSS: &str = "how far this sits from the most stable element in this universe";

/// One more than the largest seed the "surprise me" button will offer.
///
/// **Six digits, because the seed has to be copyable off a screen by hand.** §6
/// makes the seed the shareable thing, and a full `u64` is twenty digits — a
/// number nobody writes down is a universe nobody can return to or pass on.
/// Typing a long seed manually still works; this bounds only what the button
/// suggests, and a million universes is not a limit anyone will reach.
const RANDOM_SEED_CEILING: u64 = 1_000_000;

/// A number from the wall clock, for seeding the "surprise me" button.
///
/// **This is the only wall-clock read in the crate, and it is deliberately not
/// in [`ViewerState::randomise`].** §13.1's ban on wall-clock exists to protect
/// simulation *results*; this produces none — it produces an *input*, which is
/// the same category as a keystroke, and the number it picks is immediately
/// written into the seed box where the user can see it, keep it, and type it
/// back. Keeping the read here rather than inside `randomise` is what leaves
/// that function a pure function of the moment, so every test of the button can
/// hand it a fixed number instead of racing a clock.
///
/// Wrapping arithmetic rather than `as_nanos()`: that returns `u128`, and the
/// conversion would need either a cast this workspace warns on or a fallible
/// path with nothing sensible to do on failure. Wrapping is exactly right for
/// an entropy source — no value of it is wrong.
#[must_use]
pub fn moment() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |since| {
            since
                .as_secs()
                .wrapping_mul(1_000_000_000)
                .wrapping_add(u64::from(since.subsec_nanos()))
        })
}

/// Why a typed seed could not be used.
///
/// **Four variants rather than one, because the messages differ and the
/// difference is the point.** This is a program a child reads, and "that is not
/// a number" in reply to `-1` teaches something false — `-1` *is* a number, it
/// is merely not a seed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum SeedError {
    /// Nothing was typed.
    #[error("type a seed")]
    Empty,
    /// Something was typed, and it was not a whole number.
    ///
    /// The bound is interpolated from [`u64::MAX`] rather than transcribed. It
    /// is the domain of `Universe::generate`, and a typed copy would stay on
    /// screen saying something false if that domain ever narrows — with the test
    /// that pins this string passing throughout, since it pins the literal to
    /// itself.
    #[error("seeds are whole numbers, 0 to {}", u64::MAX)]
    NotANumber,
    /// A negative whole number.
    #[error("seeds are not negative")]
    Negative,
    /// A whole number too large to be a seed.
    #[error("that is larger than the largest seed, {}", u64::MAX)]
    TooLarge,
}

/// Read a seed from what the user typed.
///
/// Decimal only. **Hexadecimal is deliberately refused**, tempting though it is
/// — §6 writes a universe's address as `U-7F3A21C9@1`, so `7f3a` looks like a
/// seed. That form displays eight hex digits of a **64-bit** seed, and
/// `borbax_universe`'s own note is explicit that a truncated display address is
/// never the thing to compare. Accepting hex here would silently answer an open
/// §6 question in the top crate of the workspace.
///
/// # Errors
///
/// Returns the [`SeedError`] describing what is wrong with `input`, so the
/// panel can say something useful rather than merely refuse.
pub fn parse_seed(input: &str) -> Result<u64, SeedError> {
    let trimmed = input.trim();
    trimmed.parse::<u64>().map_err(|err| match err.kind() {
        IntErrorKind::Empty => SeedError::Empty,
        IntErrorKind::PosOverflow => SeedError::TooLarge,
        // Everything else arrives as `InvalidDigit`, which is where `-1` and
        // `hello` become indistinguishable to the standard library and have to
        // be separated by hand.
        _ => sign_or_digits(trimmed),
    })
}

/// Tell a negative number from something that merely starts with a minus sign.
///
/// **The one-line version of this is wrong and it is the version anyone would
/// write.** `if s.starts_with('-') { Negative }` says "seeds are not negative"
/// about `-abc`, which is not negative and not a number — it replaces one
/// misleading message with another. The sign only means *negative* when what
/// follows it is itself a number.
///
/// Sign beats magnitude: `-99999999999999999999999` is [`SeedError::Negative`]
/// rather than [`SeedError::TooLarge`], because the sign is the first thing
/// wrong with it and "seeds are not negative" is the more useful sentence.
///
/// `-0` is [`SeedError::Negative`] too. It is arguably zero, but the sign is
/// still the thing to remove, and a separate message for it would be a rule
/// nobody could guess.
fn sign_or_digits(trimmed: &str) -> SeedError {
    match trimmed.strip_prefix('-') {
        Some(rest) if !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()) => {
            SeedError::Negative
        }
        _ => SeedError::NotANumber,
    }
}

/// The result of the last committed seed.
///
/// **An enum rather than an `Option<Universe>` beside an `Option<SeedError>`.**
/// Those two fields have four combinations and only three are meaningful; the
/// fourth — a universe *and* an error — is the one that leaves a stale element
/// count on screen underneath a red message, attributed to a seed that did not
/// produce it.
#[derive(Debug)]
pub enum Outcome {
    /// The typed seed was refused, and why.
    Rejected(SeedError),
    /// The universe the typed seed generated, and which of its elements is
    /// selected.
    ///
    /// **Boxed**, and not for taste: `size_of::<Universe>()` is 232 bytes, so
    /// without the box `clippy::large_enum_variant` fires — every `Outcome`
    /// anywhere would carry the universe's footprint to hold a four-byte error.
    Loaded {
        /// The universe on screen.
        universe: Box<Universe>,
        /// The selected element, or `None` until a cell is clicked.
        ///
        /// **Inside the variant rather than beside it, and that is the same
        /// argument this enum was made for one field along.** A
        /// `selected: Option<ElementId>` field on [`ViewerState`] would have
        /// four combinations with `outcome` and two of them are meaningless: a
        /// selection in a universe that does not exist, and — the expensive one
        /// — a selection made in a *previous* universe, which paints one
        /// element's properties under another element's heading. Ids do not
        /// carry their table with them, so slot 40 of seed 7 and slot 40 of
        /// seed 11 are both valid and entirely different elements; nothing
        /// downstream could detect the substitution.
        ///
        /// In here it is unrepresentable. `ViewerState::reload` — private,
        /// which is why this is not a link — assigns
        /// `self.outcome` wholesale, so the old selection is dropped with the
        /// old universe — there is no `self.selected = None` line to forget at
        /// the next call site, and no rule for a future reader to remember.
        /// The same reasoning that made `reload` private: make
        /// it impossible rather than merely tested.
        selected: Option<ElementId>,
        /// The molecule on screen, laid out once when this universe loaded.
        ///
        /// **Inside the variant for exactly the reason `selected` is**, and the
        /// failure it prevents is worse. A `scene: Option<Demo>` field beside
        /// `outcome` would let the previous universe's molecule stay on screen
        /// next to the new universe's element count — a shape drawn from
        /// elements that are not the ones named beside it, which nothing
        /// downstream could detect. In here there is no line to forget: `reload`
        /// replaces the whole `Outcome`.
        ///
        /// `None` means the universe has no molecule to show. Measured over 500
        /// universes it never happened; it is data rather than a panic because
        /// "measured never" is not "cannot".
        scene: Option<Box<Demo>>,
    },
}

/// Everything the window shows, and nothing about how it is drawn.
#[derive(Debug)]
pub struct ViewerState {
    phrase_text: String,
    seed_text: String,
    outcome: Outcome,
    regenerations: u64,
    re_embeds: u64,
    reloads: u64,
}

impl Default for ViewerState {
    fn default() -> Self {
        Self::new()
    }
}

impl ViewerState {
    /// A viewer with an empty seed box.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            phrase_text: String::new(),
            seed_text: String::new(),
            outcome: Outcome::Rejected(SeedError::Empty),
            regenerations: 0,
            re_embeds: 0,
            reloads: 0,
        }
    }

    /// The state the window opens in — [`crate::OPENING_SEED`], already
    /// committed.
    ///
    /// **Separate from [`Self::new`] rather than replacing it**, and the reason
    /// is a test that would otherwise be quietly wrong: the frame tests drive a
    /// real `TextEdit` with `type_text`, which *appends*. Against a pre-filled
    /// box, typing `7` produces `17` — a different universe, a green-looking
    /// harness, and a test measuring nothing it claims to.
    #[must_use]
    pub fn opening() -> Self {
        let mut state = Self::new();
        crate::OPENING_SEED
            .to_string()
            .clone_into(&mut state.seed_text);
        state.reload();
        state
    }

    /// Pick a seed from `moment`, put it in the box, and load that universe.
    ///
    /// **The button types into the box rather than holding a hidden number.**
    /// That is the whole design: whatever it picks is visible, writeable-down
    /// and retypable, so a universe someone likes is one they can come back to
    /// and pass on. §6 makes the seed the shareable thing, and a "surprise me"
    /// that surprises you with something you cannot record has taken that away.
    ///
    /// **Randomness comes from [`borbax_rng::Stream`], never `rand::random` and
    /// never `RandomState`** — a Global Constraint, and honoured here even
    /// though this value is an input rather than a result. Two things come free
    /// from obeying it: `next_range` is modulo-bias-free, which a hand-rolled
    /// `% RANDOM_SEED_CEILING` on a clock reading is not, and there is no route
    /// by which a future edit could reach for `rand` "because the viewer
    /// already does".
    ///
    /// **[`Domain::Hash`], and the first version of this said
    /// [`Domain::Universe`], which was a real collision rather than a stylistic
    /// one.** `Stream::new(m, Domain::Universe, 0)` is bit-for-bit the stream
    /// `borbax_universe`'s own `stream(seed)` helper returns — the one element
    /// generation draws its shell pattern from. Two reviewers measured the
    /// aliasing independently and it was exact on every moment tried. Nothing
    /// was wrong on screen, because the universe is keyed on the *drawn* seed
    /// and not on `moment`; what was wrong is that a future change to element
    /// generation's first draw would silently change which universes this button
    /// offers, and `Domain`'s own doc says the enum exists to stop exactly that.
    /// `Domain::Hash` is documented for this case — "using a `Stream` as a hash
    /// function" — and turning a clock reading into a seed is that.
    ///
    /// Pure in `moment`. The clock lives in [`moment()`], at the edge.
    pub fn randomise(&mut self, moment: u64) {
        let seed = Stream::new(moment, Domain::Hash, 0).next_range(RANDOM_SEED_CEILING);
        seed.to_string().clone_into(&mut self.seed_text);
        // The button produced this seed, so no name did. Leaving a name beside
        // it would attribute a universe to something that did not make it —
        // the same class of defect `Outcome` exists to prevent one field along.
        self.phrase_text.clear();
        self.reload();
    }

    /// The text in the name box, for the panel to edit in place.
    pub const fn phrase_text_mut(&mut self) -> &mut String {
        &mut self.phrase_text
    }

    /// The text in the name box.
    #[must_use]
    pub fn phrase_text(&self) -> &str {
        &self.phrase_text
    }

    /// Put `phrase` in the name box and open the universe it names.
    ///
    /// The convenience spelling for tests and for anything holding a phrase
    /// already; the panel edits [`Self::phrase_text_mut`] in place and calls
    /// [`Self::reload_from_phrase`], which is the same body.
    pub fn set_phrase(&mut self, phrase: &str) {
        phrase.clone_into(&mut self.phrase_text);
        self.reload_from_phrase();
    }

    /// Derive the seed the name box currently names, and load that universe.
    ///
    /// **An empty or blank name box deliberately leaves the seed box alone**,
    /// rather than clearing it. The alternative was considered and rejected on
    /// the state the window opens in: it opens with the name box empty and a
    /// universe on screen, because a viewer whose first screen asks you to do
    /// something has spent its opening move on a chore. Any rule under which an
    /// empty name box clears the seed contradicts that opening state, so the
    /// two would have to disagree about the same visible input. Leaving it
    /// alone means an empty name box reads as "the name box is not currently
    /// driving", which is exactly what it means at startup.
    ///
    /// **The name is a key, not a label** — see
    /// [`borbax_universe::seed_from_phrase`] for why that was the choice, and
    /// note the direction: the phrase goes *in*, a number comes back, and no
    /// `&str` crosses into a chemistry type.
    pub fn reload_from_phrase(&mut self) {
        if let Some(seed) = seed_from_phrase(&self.phrase_text) {
            seed.to_string().clone_into(&mut self.seed_text);
            self.reload();
        }
    }

    /// The seed box was edited by hand, so nothing is named any more.
    ///
    /// **Clearing the name box is the point of this method existing**, and it
    /// is why the panel does not simply call `reload` on the seed box's
    /// `changed()` as it did in Step 1. With the name left standing, typing
    /// `42` into the seed box beside a name box reading `Emily` puts a universe
    /// on screen labelled with a name that did not produce it — a plausible
    /// attribution to the wrong thing, which is precisely the state [`Outcome`]
    /// was made an enum to prevent one field along.
    ///
    /// It cannot live inside `reload`: that is also the path
    /// [`Self::reload_from_phrase`] takes, which would wipe the name a
    /// keystroke after setting it.
    pub fn commit_typed_seed(&mut self) {
        self.phrase_text.clear();
        self.reload();
    }

    /// The text in the seed box, for the panel to edit in place.
    pub const fn seed_text_mut(&mut self) -> &mut String {
        &mut self.seed_text
    }

    /// The text in the seed box.
    #[must_use]
    pub fn seed_text(&self) -> &str {
        &self.seed_text
    }

    /// Generate the universe the seed box currently names.
    ///
    /// **Private, and that is a guard rather than tidiness.** It differs from
    /// [`Self::commit_typed_seed`] by one `phrase_text.clear()`, and calling
    /// the wrong one leaves a universe on screen labelled with a name that did
    /// not produce it. While it was public, substituting it in the panel was a
    /// *test failure*; now it is a compile error — `error[E0624]: method
    /// `reload` is private` — which catches it at every future call site rather
    /// than at the three that exist today. A reviewer measured the cost:
    /// `reload` had no caller outside this file except six lines in
    /// `tests/acceptance.rs`, all of them `seed_text_mut()` + `reload()` over
    /// an already-empty phrase box, which is exactly what `commit_typed_seed`
    /// does.
    ///
    /// **Called only when the seed text actually changes**, never once per
    /// frame. `Universe::generate` is ~146 µs; at 60 Hz an unconditional call
    /// would spend most of the frame budget regenerating a universe nobody
    /// asked for again, and by Step 2 there are 120 element cells behind it.
    /// [`Self::regenerations`] is what makes that claim testable.
    /// Select the element with this id, if this universe has such a slot.
    ///
    /// **Refused rather than clamped**, and [`borbax_universe::PeriodicTable::get`]'s
    /// own doc gives the reason: clamping "silently answers with a *different
    /// element's* properties, which nothing downstream can detect". There is
    /// nothing sensible to do on failure and nothing to report — every id the
    /// panel can produce came out of [`Self::rows`], which minted it from the
    /// table itself.
    ///
    /// Selecting while the seed is refused does nothing: there is no universe
    /// to select in, which is the state the enum makes unrepresentable.
    pub fn select(&mut self, id: ElementId) {
        if let Outcome::Loaded {
            universe, selected, ..
        } = &mut self.outcome
            && universe.table.get(id).is_some()
        {
            *selected = Some(id);
        }
    }

    /// Which element is selected, if any.
    #[must_use]
    pub const fn selected(&self) -> Option<ElementId> {
        match &self.outcome {
            Outcome::Rejected(_) => None,
            Outcome::Loaded { selected, .. } => *selected,
        }
    }

    /// The table laid out by period and group.
    ///
    /// Empty when the seed was refused, so the panel paints no grid rather than
    /// an empty one.
    ///
    /// **A single forward pass with no arithmetic**, which is what keeps the
    /// "viewer computes no physics" rule from even being a question here. The
    /// grid needs no extents: `group` *is* the column index and a change of
    /// `period` *is* the row break. Verified over 5000 universes — `period`
    /// rises by exactly one at each break, `group` runs 0,1,2,… contiguously
    /// within a period, and the first element is always (0, 0), with zero
    /// exceptions. That is forced by the generator rather than lucky: `period`
    /// and `group` are written from a shell counter and a within-shell counter
    /// that increment together.
    ///
    /// It is a property of `generate_elements` and not a documented contract of
    /// `PeriodicTable`, so `the_column_a_cell_sits_in_is_its_group` pins it.
    #[must_use]
    pub fn rows(&self) -> Vec<PeriodRow<'_>> {
        let Outcome::Loaded {
            universe, selected, ..
        } = &self.outcome
        else {
            return Vec::new();
        };
        let mut rows: Vec<PeriodRow<'_>> = Vec::new();
        for (id, element) in universe.table.iter() {
            let cell = Cell {
                id,
                symbol: &element.symbol,
                selected: *selected == Some(id),
            };
            match rows.last_mut() {
                Some(row) if row.period == element.period => row.cells.push(cell),
                _ => rows.push(PeriodRow {
                    period: element.period,
                    shell_label: format!("shell {}", u16::from(element.period) + 1),
                    cells: vec![cell],
                }),
            }
        }
        rows
    }

    /// The heading over the properties block: `Mx · meax`, or the empty prompt.
    ///
    /// **Symbol first.** It is what was just clicked, so the heading confirms
    /// rather than asking the reader to look something up — and it is the half
    /// that is unique by construction. `·` because [`Self::status_line`]
    /// already separates with it.
    #[must_use]
    pub fn selection_heading(&self) -> String {
        self.selected_element().map_or_else(
            || "pick an element".to_owned(),
            // Two elements of one universe may share a name (measured: 16.05%
            // of universes), and that is deliberately not surfaced. Only one
            // heading is on screen at a time, so the collision is invisible;
            // and a badge announcing it would teach that Borbax is broken, when
            // in fact G2 never required names to be unique.
            |element| format!("{} · {}", element.symbol, element.name),
        )
    }

    /// The nine property rows, with values when something is selected.
    ///
    /// **The same nine labels either way**, so the block has a fixed footprint
    /// and the grid above it never reflows on the first click.
    #[must_use]
    pub fn selection_properties(&self) -> Vec<PropertyRow> {
        let element = self.selected_element();
        let shells = self.rows().len();
        PROPERTY_LABELS
            .iter()
            .enumerate()
            .map(|(i, label)| {
                let (value, unit) = element
                    .map_or_else(|| (EMPTY_VALUE.to_owned(), ""), |e| property(i, e, shells));
                PropertyRow {
                    label,
                    value,
                    unit: if element.is_some() { unit } else { "" },
                    gloss: match (i, element.is_some()) {
                        (6, true) => Some(SURFACE_GLOSS),
                        (8, true) => Some(INSTABILITY_GLOSS),
                        _ => None,
                    },
                }
            })
            .collect()
    }

    /// The selected element, if there is one.
    fn selected_element(&self) -> Option<&Element> {
        match &self.outcome {
            Outcome::Rejected(_) => None,
            Outcome::Loaded {
                universe, selected, ..
            } => selected.and_then(|id| universe.table.get(id)),
        }
    }

    fn reload(&mut self) {
        // **At the entry, deliberately, and this is the one counter here that
        // should be.** Its job is to tell the scene that what it is drawing is
        // out of date, and a *refused* seed makes it out of date just as surely
        // as a loaded one: the previous universe's molecule must come off the
        // screen beside the refusal, which is `Outcome`'s whole argument carried
        // one layer down into the entity cache.
        //
        // Gating the scene on `re_embeds` instead left the old atoms on screen
        // under the refusal message, with 532 tests green — exactly what the
        // Step 3 memo's manual checklist item 8 says must not happen.
        self.reloads = self.reloads.saturating_add(1);
        self.outcome = match parse_seed(&self.seed_text) {
            Ok(seed) => {
                let universe = Universe::generate(seed);
                // **Boxed**, for the reason `universe` beside it is: an
                // `Embedding` carries fixed-size arrays for the largest
                // molecule the graph type allows, so an unboxed one would make
                // every `Outcome` anywhere carry that footprint to hold a
                // four-byte refusal.
                let scene = crate::molecule::build(&universe).map(Box::new);
                // **Inside the `is_some`, and the previous version was not.**
                // `build` returns `None` without ever calling `embed`, so an
                // unconditional increment counted *attempts* while three
                // comments — including the one that used to sit here, warning
                // against this exact conflation — claimed it counted layouts.
                // Five review lanes found it independently.
                if scene.is_some() {
                    self.re_embeds = self.re_embeds.saturating_add(1);
                }
                // Incremented **here**, adjacent to the call it counts, and not
                // at the top of this function. An increment at the entry counts
                // *commits*, which is a different quantity that happens to
                // agree today — and would keep agreeing right up until parsing
                // and generation are split, at which point the counter would
                // report the wrong thing while its test stayed green.
                self.regenerations = self.regenerations.saturating_add(1);
                // **Nothing is selected in a universe nobody has looked at
                // yet**, and this is the line that makes a stale selection
                // unrepresentable rather than merely cleared: the whole
                // `Outcome` is replaced, so the previous universe's selection
                // is dropped with the universe it belonged to.
                Outcome::Loaded {
                    universe: Box::new(universe),
                    selected: None,
                    scene,
                }
            }
            Err(refusal) => Outcome::Rejected(refusal),
        };
    }

    /// The last committed seed's result.
    #[must_use]
    pub const fn outcome(&self) -> &Outcome {
        &self.outcome
    }

    /// How many universes have been generated since this state was created.
    ///
    /// **Shipped rather than `#[cfg(test)]`-gated**, which is deliberate:
    /// integration tests compile the library *without* `cfg(test)`, so a
    /// test-gated field is invisible to exactly the tests that need it. It is
    /// one `u64` in a crate that cannot affect a simulation result.
    ///
    /// **What it does and does not catch.** It catches `reload`
    /// migrating into the per-frame body. It does **not** catch a direct
    /// `Universe::generate` call added to [`crate::panel::draw`], which never
    /// touches `reload` and so never touches this counter — that variant is
    /// held by the seam rule instead: `Universe::generate` has exactly one call
    /// site in this crate, and it is in this file.
    #[must_use]
    pub const fn regenerations(&self) -> u64 {
        self.regenerations
    }

    /// How many molecules have been laid out since this state was created.
    ///
    /// **The §8.6 guard, and it now counts layouts rather than attempts.** It
    /// rises only when `build` returned a molecule, so it is `embed`'s own
    /// count: `embed` runs 240 fixed solver iterations and `canonicalise`
    /// searches for a labelling, and both are per-species work that must never
    /// reach a frame.
    ///
    /// **This is not the scene's change token** — [`Self::reloads`] is. The two
    /// were the same field until a review found that gating the scene on this
    /// one left the previous universe's molecule on screen under a refusal.
    ///
    /// **What it cannot see, stated so nobody reads it as covering more than it
    /// does.** It does not see an `embed` written straight into a paint body or
    /// a frame system — such a call never touches `reload`, so it never touches
    /// this counter. That variant is held by `cargo xtask` instead, which
    /// requires exactly one call site in this crate. It also says nothing about
    /// per-frame *scene* cost: rebuilding the atom entities every frame calls no
    /// chemistry at all and would leave this number flat.
    #[must_use]
    pub const fn re_embeds(&self) -> u64 {
        self.re_embeds
    }

    /// How many times a seed has been committed, loaded or refused.
    ///
    /// **The scene's change token, and it counts commits on purpose.** Every
    /// other counter here is adjacent to the call it counts, because counting
    /// commits would measure the wrong quantity. This one wants commits: the
    /// question it answers is "is what the scene is drawing still the right
    /// thing", and a refused seed changes that answer without laying anything
    /// out.
    #[must_use]
    pub const fn reloads(&self) -> u64 {
        self.reloads
    }

    /// The molecule to draw, if this universe has one.
    #[must_use]
    pub const fn demo(&self) -> Option<&Demo> {
        match &self.outcome {
            Outcome::Rejected(_) => None,
            Outcome::Loaded { scene, .. } => match scene {
                Some(demo) => Some(demo),
                None => None,
            },
        }
    }

    /// The line printed under the scene, saying what is being looked at.
    ///
    /// **Without it the scene reads as "the element you clicked, in 3D"**, which
    /// is false — the molecule is built from this universe's lowest elements
    /// that fit the shape, and has nothing to do with the selection. It sits
    /// outside the viewport, because there is no text inside the viewport.
    ///
    /// Built here rather than where it is painted, so a test asserting on this
    /// string is asserting on what the window shows.
    #[must_use]
    pub fn scene_caption(&self) -> String {
        let Outcome::Loaded {
            universe, scene, ..
        } = &self.outcome
        else {
            return String::new();
        };
        let Some(demo) = scene else {
            return "this universe builds no molecule".to_owned();
        };
        let mut symbols = demo
            .elements
            .iter()
            .filter_map(|id| universe.table.get(*id))
            .map(|e| e.symbol.as_str());
        let Some(centre) = symbols.next() else {
            return "this universe builds no molecule".to_owned();
        };
        let leaves = symbols.collect::<Vec<_>>().join(", ");
        format!("{centre} holding {leaves}")
    }

    /// The one line the panel paints under the seed box.
    ///
    /// **Formatting lives here rather than in [`crate::panel`]** so that a test
    /// asserting on this string is asserting on the thing the window shows. A
    /// `format!` inlined into the panel would leave every test here green over a
    /// window they had stopped describing.
    #[must_use]
    pub fn status_line(&self) -> String {
        match &self.outcome {
            Outcome::Rejected(refusal) => refusal.to_string(),
            // **`u8::from(physics)`, never the literal `"v1"`.** §6 is explicit
            // that a universe is `(seed, physics)` and not `seed` alone, so a
            // viewer showing only the seed teaches the exact misidentification
            // that field exists to prevent. Reading the discriminant is what
            // makes the label move on the day `PhysicsVersion::CURRENT` does;
            // a typed `"v1"` would satisfy every test here while quietly
            // becoming a lie, which is the whole difference between a guard and
            // a decoration.
            //
            // The full `U-7F3A21C9@N` address form is deliberately not used:
            // it shows 32 bits of a 64-bit seed, and no encoder for it exists
            // anywhere in the workspace. Inventing one here would settle an
            // open §6 question in the crate least entitled to settle it.
            Outcome::Loaded { universe, .. } => format!(
                "{} elements · physics v{}",
                universe.table.len(),
                u8::from(universe.physics)
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{SeedError, parse_seed};

    #[test]
    fn a_bare_decimal_seed_is_accepted_as_written() {
        assert_eq!(parse_seed("7"), Ok(7));
        assert_eq!(parse_seed("0"), Ok(0));
        assert_eq!(parse_seed("18446744073709551615"), Ok(u64::MAX));
        // `u64::MAX` alone does **not** catch a round-trip through `f64`:
        // Rust's saturating float-to-int cast returns `u64::MAX` for it anyway,
        // so the test would pass over the defect. These two are the inputs that
        // kill it — `2^53 + 1` is the first integer `f64` cannot represent, and
        // `u64::MAX - 1` comes back as `u64::MAX`.
        assert_eq!(parse_seed("18446744073709551614"), Ok(u64::MAX - 1));
        assert_eq!(parse_seed("9007199254740993"), Ok(9_007_199_254_740_993));
    }

    #[test]
    fn an_empty_seed_box_says_it_is_empty_rather_than_not_a_number() {
        assert_eq!(parse_seed(""), Err(SeedError::Empty));
        assert_eq!(parse_seed("   "), Err(SeedError::Empty));
    }

    #[test]
    fn an_empty_seed_box_does_not_quietly_mean_universe_zero() {
        // The defect: `if s.is_empty() { return Ok(0) }`. It puts a real
        // element count on screen for a seed nobody typed, which is the same
        // class as a refused seed leaving a stale number — a plausible number
        // attributed to nothing.
        assert!(parse_seed("").is_err());
    }

    #[test]
    fn letters_in_the_seed_box_are_refused_not_hashed() {
        assert_eq!(parse_seed("hello"), Err(SeedError::NotANumber));
        // The sharp one. A decimal parse rejects `7f3a`; `from_str_radix(s, 16)`
        // accepts it. It kills two tempting implementations at once — hex input
        // (tempting because §6 spells addresses `U-7F3A21C9@1`) and "be
        // friendly, hash whatever was typed into a seed", which would be a
        // display choice reaching a `borbax-*` call.
        assert_eq!(parse_seed("7f3a"), Err(SeedError::NotANumber));
        assert_eq!(parse_seed("0x1f"), Err(SeedError::NotANumber));
        assert_eq!(parse_seed("U-7F3A21C9"), Err(SeedError::NotANumber));
    }

    #[test]
    fn a_seed_larger_than_the_universe_can_hold_says_so() {
        assert_eq!(parse_seed("18446744073709551616"), Err(SeedError::TooLarge));
        // Catches `.parse().unwrap_or(u64::MAX)` and `.unwrap_or(0)`, both of
        // which generate a real universe for a seed the user did not type.
        assert!(parse_seed("18446744073709551616").is_err());
    }

    #[test]
    fn a_minus_sign_is_only_negative_when_what_follows_is_a_number() {
        assert_eq!(parse_seed("-1"), Err(SeedError::Negative));
        assert_eq!(parse_seed("-99"), Err(SeedError::Negative));
        // The discriminator for `sign_or_digits`: the one-line
        // `starts_with('-')` shortcut calls these negative, and "seeds are not
        // negative" is exactly as untrue of `-abc` as "that is not a number"
        // was of `-1`. Nothing else on this list catches it.
        assert_eq!(parse_seed("-abc"), Err(SeedError::NotANumber));
        assert_eq!(parse_seed("-"), Err(SeedError::NotANumber));
        assert_eq!(parse_seed("- 1"), Err(SeedError::NotANumber));
        // Sign beats magnitude.
        assert_eq!(
            parse_seed("-99999999999999999999999"),
            Err(SeedError::Negative)
        );
        // Measured, not assumed: `u64::from_str` reports `-0` as `InvalidDigit`
        // like every other signed spelling, so it lands on the sign rule.
        assert_eq!(parse_seed("-0"), Err(SeedError::Negative));
    }

    #[test]
    fn whitespace_around_a_pasted_seed_is_trimmed_and_nothing_else_is() {
        // This test pins a decision rather than discovering one, and would pass
        // for either decision once written to match it. Its value is that a
        // later `retain(|c| !c.is_whitespace())` "cleanup" turns `7 7` into
        // `77` — a different universe, silently.
        assert_eq!(parse_seed(" 7 "), Ok(7));
        assert_eq!(parse_seed("7\n"), Ok(7));
        assert_eq!(parse_seed("\t7"), Ok(7));
        assert_eq!(parse_seed("7 7"), Err(SeedError::NotANumber));
        assert_eq!(parse_seed("1_000"), Err(SeedError::NotANumber));
        assert_eq!(parse_seed("7."), Err(SeedError::NotANumber));
        // Measured: `u64::from_str` accepts a leading `+`. Left accepted — it
        // is unambiguous and refusing it would be a rule with no reason.
        assert_eq!(parse_seed("+7"), Ok(7));
    }

    #[test]
    fn every_refusal_says_something_different() {
        // A guard on the guard: four variants exist *because* the messages
        // differ, so four identical messages would make the whole enum
        // decoration. Catches a copy-paste in the `#[error(..)]` attributes.
        let messages = [
            SeedError::Empty.to_string(),
            SeedError::NotANumber.to_string(),
            SeedError::Negative.to_string(),
            SeedError::TooLarge.to_string(),
        ];
        let distinct: std::collections::BTreeSet<&String> = messages.iter().collect();
        assert_eq!(
            distinct.len(),
            messages.len(),
            "two refusals share a message: {messages:?}"
        );
    }
}
