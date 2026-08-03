//! What is selected, and every decision the viewer makes.
//!
//! **This file imports no UI type**, which is the whole reason it exists
//! separately: it is the half of the crate a test can reach without a window.

use core::num::IntErrorKind;
use std::time::{SystemTime, UNIX_EPOCH};

use borbax_rng::{Domain, Stream};
use borbax_universe::{Universe, seed_from_phrase};
use thiserror::Error;

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
    /// The universe the typed seed generated.
    ///
    /// **Boxed**, and not for taste: `size_of::<Universe>()` is 232 bytes, so
    /// without the box `clippy::large_enum_variant` fires — every `Outcome`
    /// anywhere would carry the universe's footprint to hold a four-byte error.
    Loaded(Box<Universe>),
}

/// Everything the window shows, and nothing about how it is drawn.
#[derive(Debug)]
pub struct ViewerState {
    phrase_text: String,
    seed_text: String,
    outcome: Outcome,
    regenerations: u64,
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
    fn reload(&mut self) {
        self.outcome = match parse_seed(&self.seed_text) {
            Ok(seed) => {
                let universe = Universe::generate(seed);
                // Incremented **here**, adjacent to the call it counts, and not
                // at the top of this function. An increment at the entry counts
                // *commits*, which is a different quantity that happens to
                // agree today — and would keep agreeing right up until parsing
                // and generation are split, at which point the counter would
                // report the wrong thing while its test stayed green.
                self.regenerations = self.regenerations.saturating_add(1);
                Outcome::Loaded(Box::new(universe))
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
            Outcome::Loaded(universe) => format!(
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
