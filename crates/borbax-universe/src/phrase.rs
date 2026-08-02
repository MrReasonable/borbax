//! Turning a typed phrase into a seed, so a child can type her name and get a
//! universe.
//!
//! # Why a derived seed rather than a random one with a label attached
//!
//! Ian's decision, and it was put to him as a choice: *"Emily" always gives
//! Emily's universe, on any machine, forever* — the name is a **key**, not a
//! label. The alternative was to draw a random seed and let the phrase name it
//! afterwards, which is a "surprise me" button with a text box in front of it,
//! and that button already exists. Deriving costs nothing in variety (`"Emily"`
//! and `"Emilz"` come out unrelated) and buys reproducibility and sharing,
//! which is §6's whole point.
//!
//! # Why this lives in a chemistry crate and returns a number
//!
//! Hashing a phrase inside the viewer would be a display choice reaching a
//! `borbax-*` call. It returns the **seed**, never a [`crate::Universe`], for
//! two reasons that both matter: the viewer has to show the number so it can be
//! written down (§6 makes the seed the shareable thing), and a
//! `Universe::from_phrase` would be a second `Universe::generate` call site in
//! a crate whose `xtask` guard asserts there is exactly one.
//!
//! **The phrase is destroyed here, at the crate boundary.** No `&str` and no
//! `String` reaches any chemistry type, which makes G1, G3, G5 and G6 (§5)
//! structurally unreachable from this feature rather than merely unviolated by
//! it. The check is a grep, not a reading: a `&str` field on a type in a
//! `borbax-*` crate is the defect, and there is none.
//!
//! The phrase is deliberately **not** filtered through `naming`'s G2 blocklist.
//! G2 constrains *generated* element names; a user typing "Carbon" into a text
//! box breaches nothing, and refusing it would be a program telling a child her
//! name is not allowed.
//!
//! # The normalisation rule, in one sentence
//!
//! *Capitals do not matter, and spaces at the beginning and end do not matter;
//! everything else does.*
//!
//! Spelling: trim, collapse internal whitespace runs to a single space, fold
//! ASCII case, and leave every other byte alone.
//!
//! **`to_ascii_lowercase`, never [`str::to_lowercase`], and this is a
//! determinism decision rather than a simplification.** Rust's Unicode tables
//! move between compiler releases — `UNICODE_VERSION` was 16.0.0 and is 17.0.0
//! on the pinned toolchain — so a Unicode-aware fold would let an innocent
//! toolchain bump silently reassign every universe anyone had already named.
//! ASCII folding cannot move. The accepted cost is that `É` and `é` are
//! different universes, which is consistent with the no-normalisation decision
//! below anyway.
//!
//! **No Unicode normalisation.** The NFC and NFD spellings of `café` are
//! different universes, and that is pinned in `GOLDEN` as a deliberate
//! recorded fact rather than left to be discovered. Adding
//! `unicode-normalization` would introduce a *second* Unicode clock, movable by
//! a lockfile-only dependency bump, into a result-affecting path.
//!
//! **One residual Unicode dependency is admitted rather than hidden.**
//! [`str::split_whitespace`] uses `char::is_whitespace`, which is the Unicode
//! `White_Space` property and therefore a table, not a constant. That is the
//! ruling and it is a far more stable table than the case mappings — but
//! "far more stable" is not "pinned", so the hazard is converted from silent to
//! loud instead of being argued away: `GOLDEN` carries a `U+00A0` row whose
//! value is shared with its plain-space twin, so a release that changed
//! `White_Space` for that character fails the golden rather than quietly
//! reassigning universes.
//!
//! # Why the mixing is a `Stream` and not a hasher
//!
//! §13.1 bans `DefaultHasher` (unspecified across releases) and `RandomState`
//! (randomly seeded per process); as of the commit before this one, that ban is
//! enforced by `clippy.toml` and `cargo xtask` rather than by a comment.
//! [`Domain::Hash`] exists for precisely this — its own documentation sanctions
//! "using a `Stream` as a hash function" — and it is separated from the live
//! domains so that adding a draw somewhere in the physics cannot move the
//! mapping.

use borbax_rng::{Domain, Stream};

/// The universe a phrase names, or `None` if the phrase names nothing.
///
/// Returns `None` **only** for a phrase that is empty once normalised — empty,
/// or nothing but whitespace.
///
/// # Why fallible
///
/// An infallible version would have to invent a universe for the empty string,
/// which reopens `an_empty_seed_box_does_not_quietly_mean_universe_zero` one
/// crate down: a real element count on screen, attributed to a name nobody
/// typed. The plan's sketch said `-> u64`; it was written before that test
/// existed.
///
/// # Why byte by byte, and never packed into words
///
/// Each byte is absorbed on its own, so **there is no byte order to get
/// wrong** — UTF-8 is a byte sequence and `u64::from(u8)` is a widening
/// conversion. This matters because it cannot be tested here: §13.4's matrix is
/// three platforms that are all little-endian and all 64-bit, so a
/// `to_ne_bytes` or a `usize` accumulator would be invisible to every leg of
/// it. The mitigation is structural rather than testable, which is also why the
/// accumulator's type is written out rather than inferred.
///
/// # Why `Stream::sub` rather than `Stream::new`
///
/// `Stream::new(acc, Domain::Hash, u64::from(byte))` with `byte == 0` — and a
/// `&str` may contain a NUL — is **bit-for-bit** the stream the viewer's
/// "surprise me" button draws from. `sub` offsets its coordinate by one, so
/// `counter[2]` is 1 here and 0 there, and the two constructions cannot alias
/// for *any* byte. That is a property of the construction rather than of the
/// values, which is the only kind worth relying on.
///
/// The accumulator enters as the **key**, so every byte re-keys the stream;
/// this is a chaining construction, not a sum, and a prefix therefore cannot
/// share a seed with what follows it.
///
/// The initial accumulator is `0` and is deliberately not a magic constant:
/// an unexplained number in a result-affecting path is a thing someone must
/// later justify, and every value here is equally arbitrary.
///
/// # Examples
///
/// ```
/// use borbax_universe::seed_from_phrase;
///
/// // Capitals and surrounding spaces do not matter.
/// assert_eq!(seed_from_phrase("Emily"), seed_from_phrase(" emily "));
/// // Anything else does.
/// assert_ne!(seed_from_phrase("Emily"), seed_from_phrase("Emilz"));
/// // A phrase that names nothing.
/// assert_eq!(seed_from_phrase("   "), None);
/// ```
#[must_use]
pub fn seed_from_phrase(phrase: &str) -> Option<u64> {
    let normalised = normalise(phrase);
    if normalised.is_empty() {
        return None;
    }

    let mut acc: u64 = 0;
    for &byte in normalised.as_bytes() {
        acc = absorb(acc, byte);
    }
    Some(acc)
}

/// Fold one byte into the accumulator.
///
/// **Extracted rather than written inline, and the reason is that the property
/// below is about this function and nothing downstream of it can see it.** The
/// choice of `Stream::sub` over `Stream::new` is invisible in
/// [`seed_from_phrase`]'s output: `GOLDEN` pins whatever was minted, so a
/// version built on `Stream::new` would have its own self-consistent table and
/// every test in this file would be green. Inline, the decision would be
/// defended by a comment. Here it is defended by
/// `absorbing_a_nul_does_not_draw_the_surprise_buttons_stream`.
///
/// The hazard, concretely: `Stream::new(acc, Domain::Hash, u64::from(byte))`
/// with `byte == 0` — and a `&str` may contain a NUL — is bit-for-bit the
/// stream the viewer's "surprise me" button draws from. `sub` offsets its
/// coordinate by one, so `counter[2]` is 1 here and 0 there, and the two
/// constructions cannot alias for *any* byte or *any* accumulator. That is a
/// property of the construction rather than of the values.
fn absorb(acc: u64, byte: u8) -> u64 {
    Stream::sub(acc, Domain::Hash, u64::from(byte), 0).next_u64()
}

/// Trim, collapse internal whitespace to single spaces, and fold ASCII case.
///
/// Written as an explicit fold rather than `split_whitespace().collect::<Vec<_>>()
/// .join(" ")` so there is no intermediate allocation per word, and
/// `make_ascii_lowercase` is applied to the finished string rather than per
/// word because the two orders are identical — ASCII case folding neither
/// creates nor destroys whitespace — and doing it once is the version with no
/// claim to check.
fn normalise(phrase: &str) -> String {
    let mut out = String::with_capacity(phrase.len());
    for word in phrase.split_whitespace() {
        if !out.is_empty() {
            out.push(' ');
        }
        out.push_str(word);
    }
    out.make_ascii_lowercase();
    out
}

#[cfg(test)]
mod tests {
    use super::{absorb, seed_from_phrase};
    use borbax_rng::{Domain, Stream};
    use std::collections::BTreeSet;

    /// Shorthand for a phrase that is expected to name something.
    fn seed(phrase: &str) -> u64 {
        seed_from_phrase(phrase)
            .unwrap_or_else(|| unreachable!("`{phrase:?}` should name a universe"))
    }

    #[test]
    fn a_phrase_that_is_nothing_but_space_names_no_universe() {
        assert_eq!(seed_from_phrase(""), None);
        assert_eq!(seed_from_phrase(" "), None);
        assert_eq!(seed_from_phrase("   \t\n  "), None);
        // The discriminator: an infallible `-> u64` returns *some* universe for
        // the empty string, which is a real element count on screen attributed
        // to a name nobody typed.
        assert!(seed_from_phrase("").is_none());
    }

    /// The one decision in this file that its own output cannot defend.
    ///
    /// `Stream::new(m, Domain::Hash, 0)` is exactly what the viewer's "surprise
    /// me" button draws from. Had [`absorb`] been written with `Stream::new`,
    /// a phrase containing a NUL would draw that identical stream whenever the
    /// accumulator happened to equal the clock reading the button used —
    /// nothing visibly wrong, and a live coupling between a name and a button.
    ///
    /// Checked at several accumulators rather than one, because the aliasing is
    /// a property of the *pair* `(acc, byte)` and a single value would leave a
    /// reader unable to tell a structural guarantee from a coincidence.
    #[test]
    fn absorbing_a_nul_does_not_draw_the_surprise_buttons_stream() {
        for acc in [0, 1, 7, 123_456, u64::MAX] {
            assert_ne!(
                absorb(acc, 0),
                Stream::new(acc, Domain::Hash, 0).next_u64(),
                "absorbing a NUL at acc={acc} draws the random button's stream"
            );
        }
    }

    #[test]
    fn capitals_do_not_matter() {
        assert_eq!(seed("Emily"), seed("emily"));
        assert_eq!(seed("EMILY"), seed("emily"));
        assert_eq!(seed("eMiLy"), seed("emily"));
    }

    #[test]
    fn spaces_at_the_ends_do_not_matter() {
        assert_eq!(seed(" Emily"), seed("Emily"));
        assert_eq!(seed("Emily "), seed("Emily"));
        assert_eq!(seed("\t Emily \n"), seed("Emily"));
    }

    #[test]
    fn a_run_of_spaces_inside_a_name_counts_as_one() {
        assert_eq!(seed("emily  rose"), seed("emily rose"));
        assert_eq!(seed("emily\trose"), seed("emily rose"));
        assert_eq!(seed("emily \n rose"), seed("emily rose"));
    }

    /// The other half of the rule, and the half a "be helpful" edit erodes.
    ///
    /// Collapsing a run of spaces is one step from deleting them, and
    /// `retain(|c| !c.is_whitespace())` would make `"emily rose"` and
    /// `"emilyrose"` the same universe — two different names silently sharing
    /// one chemistry.
    #[test]
    fn everything_that_is_not_case_or_an_end_space_matters() {
        assert_ne!(seed("emily rose"), seed("emilyrose"));
        assert_ne!(seed("emily"), seed("emilyy"));
        assert_ne!(seed("emily"), seed("emliy"));
        assert_ne!(seed("emily"), seed("emily!"));
        // Byte-faithful: a NUL is part of the name, not a terminator. This is
        // what a C-style `strlen` view of the string would lose, and it is why
        // `Stream::sub` rather than `Stream::new` is the constructor.
        assert_ne!(seed("a\0b"), seed("ab"));
    }

    /// Accents are not folded, and this pins the accepted cost of ASCII case.
    ///
    /// `str::to_lowercase` would make the first pair equal. It is banned here
    /// because Rust's Unicode tables move between releases, which would
    /// reassign named universes on a toolchain bump.
    #[test]
    fn accents_are_part_of_the_name() {
        assert_ne!(seed("CAFÉ"), seed("café"));
        assert_ne!(seed("café"), seed("cafe"));
    }

    /// Two Unicode spellings of one visible name are two universes.
    ///
    /// Deliberate, recorded rather than discovered: normalising would add a
    /// second Unicode clock movable by a lockfile-only dependency bump.
    #[test]
    fn the_composed_and_decomposed_spellings_differ() {
        assert_ne!(seed("caf\u{e9}"), seed("cafe\u{301}"));
    }

    /// A chained fold, not a sum: no prefix shares a seed with what follows it.
    ///
    /// The mutations this catches are an accumulator that keeps only the last
    /// byte, and any commutative combine — under which `"emil"` and `"lime"`
    /// would collide.
    #[test]
    fn no_prefix_and_no_anagram_shares_a_seed() {
        let prefixes = ["e", "em", "emi", "emil", "emily"];
        let distinct: BTreeSet<u64> = prefixes.iter().map(|p| seed(p)).collect();
        assert_eq!(distinct.len(), prefixes.len(), "a prefix shares a seed");

        assert_ne!(seed("emil"), seed("lime"));
        assert_ne!(seed("ab"), seed("ba"));
    }

    /// Names a child would plausibly type land on distinct universes.
    ///
    /// A weak-mixing implementation — folding bytes with `+` or `^` — collides
    /// short similar names constantly. Measured on this corpus: zero
    /// collisions, and the bar is exact rather than a fraction because any
    /// collision at this size is a defect, not noise.
    #[test]
    fn a_corpus_of_ordinary_names_does_not_collide() {
        let mut names = Vec::new();
        for first in [
            "emily", "olivia", "amelia", "isla", "ava", "ivy", "freya", "lily", "mia", "willow",
        ] {
            names.push(first.to_owned());
            for second in ["rose", "grace", "may", "jane"] {
                names.push(format!("{first} {second}"));
            }
        }
        let distinct: BTreeSet<u64> = names.iter().map(|n| seed(n)).collect();
        assert_eq!(
            distinct.len(),
            names.len(),
            "{} names produced {} seeds",
            names.len(),
            distinct.len()
        );
    }

    /// A one-bit change to the input changes about half the output bits.
    ///
    /// **Both bounds were set after measuring, and the measurement is recorded
    /// here so the next reader can see which is which.** My first attempt
    /// asserted that every trial landed in `24..=40` — a bar chosen before any
    /// number existed, which failed immediately on a perfectly good mixer at 23.
    /// It deserved to: the per-trial minimum over 40 trials is an *order
    /// statistic*, not an estimate of the mean, and for a good mixer
    /// `P(X <= 23)` is about 1.4% per trial, so one trial at or below 23 turns
    /// up 42% of the time.
    ///
    /// The two bars therefore measure different things:
    ///
    /// - the **mean** is the estimate, with standard error `4/sqrt(40)` = 0.63,
    ///   so `30..=34` is roughly three standard errors either side. Measured:
    ///   **32.00**, the ideal value.
    /// - each **trial** gets a deliberately wide `12..=52`, about five standard
    ///   deviations, which no good mixer will ever trip and which an identity
    ///   (0 bits moved) or an inverter (64) fails on the first sample.
    ///   Measured: min 23, max 39.
    ///
    /// **Two-sided on purpose.** A floor alone passes for the mixer's opposite:
    /// an implementation flipping *every* bit is as broken as one flipping none,
    /// and only an upper bound sees it.
    #[test]
    fn one_changed_byte_moves_about_half_the_output_bits() {
        let base = "emily";
        let mut moved_bits = Vec::new();
        for position in 0..base.len() {
            for bit in 0..8_u8 {
                let mut bytes = base.as_bytes().to_vec();
                let Some(target) = bytes.get_mut(position) else {
                    continue;
                };
                // Flip one bit within the low five, which keeps the byte inside
                // ASCII — so the result is still valid UTF-8, still not
                // whitespace, and still survives normalisation unchanged.
                *target ^= 1 << (bit % 5);
                let Ok(flipped) = std::str::from_utf8(&bytes) else {
                    continue;
                };
                if flipped == base {
                    continue;
                }
                moved_bits.push((seed(base) ^ seed(flipped)).count_ones());
            }
        }

        assert!(
            moved_bits.len() >= 30,
            "only {} usable flips — the corpus filter selected almost nothing, \
             which would make every bound below vacuous",
            moved_bits.len()
        );
        for &moved in &moved_bits {
            assert!(
                (12..=52).contains(&moved),
                "a single flip moved {moved} of 64 bits: {moved_bits:?}"
            );
        }

        let total: u32 = moved_bits.iter().sum();
        let count = u32::try_from(moved_bits.len()).unwrap_or(u32::MAX);
        let mean = f64::from(total) / f64::from(count);
        assert!(
            (30.0..=34.0).contains(&mean),
            "mean bits moved {mean:.2} over {count} flips: {moved_bits:?}"
        );
    }

    /// The seed is a full 64-bit value and is deliberately not capped.
    ///
    /// `RANDOM_SEED_CEILING` in the viewer exists because a *suggested* seed has
    /// to be copyable off a screen by hand. With a phrase, the phrase is the
    /// memorable artefact and the number is a receipt, so capping buys nothing
    /// and costs collisions — at 10^6 a thousand names collide about 40% of the
    /// time.
    #[test]
    fn the_seed_is_not_squeezed_into_six_digits() {
        assert!(
            super::GOLDEN.iter().any(|&(_, value)| value > 1_000_000),
            "no golden row exceeds the viewer's suggestion ceiling"
        );
        assert!(
            super::GOLDEN
                .iter()
                .any(|&(_, value)| value > u64::from(u32::MAX)),
            "no golden row exceeds 32 bits"
        );
    }
}

/// The mapping, pinned to literal values.
///
/// **This is the point of the whole file.** The mapping is result-affecting:
/// changing it silently reassigns every universe anyone has named, and the
/// symptom — "the universe I wrote down is not the one I get" — has no other
/// guard. A round-trip test cannot see it, because a round-trip is satisfied by
/// any self-consistent mapping.
///
/// **Its authority is the mint, not the test.** These numbers were produced by
/// running the implementation, so if the first implementation was wrong but
/// stable, this table records the wrong mapping and stays green forever. What
/// it protects against is *drift*, which is the failure that actually happens.
/// The properties above are what argue the mapping is right; this argues only
/// that it has not moved.
///
/// It is an ordinary workspace `#[test]`, deliberately, so all three CI legs
/// run it — a `goldens --emit` artefact is compared on one.
///
/// **The non-ASCII rows are load-bearing and are not decoration.** Two
/// plausible wrong implementations fail *only* on them, and would otherwise
/// ship green: absorbing `char`s instead of bytes agrees with this table on
/// every ASCII row, and swapping `to_ascii_lowercase` for `str::to_lowercase`
/// changes nothing outside ASCII either. Both are caught below.
#[cfg(test)]
const GOLDEN: &[(&str, u64)] = &[
    // Three spellings of one name, sharing one value: capitals and end spaces
    // do not matter. A shared literal is the assertion — writing three
    // different numbers here would be the defect.
    ("emily", 15_709_401_653_729_972_761),
    ("Emily", 15_709_401_653_729_972_761),
    (" Emily ", 15_709_401_653_729_972_761),
    // A run of internal whitespace collapses to one space...
    ("emily rose", 5_423_473_925_176_484_554),
    ("emily  rose", 5_423_473_925_176_484_554),
    // ...but is not deleted. `retain(|c| !c.is_whitespace())` makes this row
    // equal to the two above, which is two names sharing one chemistry.
    ("emilyrose", 8_835_804_968_681_729_229),
    // Non-ASCII, and these three rows are why the table has any. Absorbing
    // `char`s instead of bytes agrees with every ASCII row above and disagrees
    // here — `é` is one `char` and two bytes.
    ("caf\u{e9}", 13_581_006_170_393_591_422),
    // The same visible word, decomposed. Different universe, deliberately: no
    // Unicode normalisation, because it would add a second Unicode clock
    // movable by a lockfile-only dependency bump.
    ("cafe\u{301}", 1_970_033_413_177_415_331),
    ("cafe", 2_972_208_318_086_301_641),
    // `ΟΔΥΣΣΕΥΣ` and `οδυσσευς`, the second with a word-final sigma. These two
    // rows are the discriminator for `to_ascii_lowercase`: swap it for
    // `str::to_lowercase` and the first row folds onto the second — including
    // the final-sigma special case, which Rust's `to_lowercase` implements —
    // so the first row's literal stops matching. Nothing else in this file
    // catches that swap, because it is a no-op on every ASCII row.
    (
        "\u{39f}\u{394}\u{3a5}\u{3a3}\u{3a3}\u{395}\u{3a5}\u{3a3}",
        4_886_566_843_106_266_812,
    ),
    (
        "\u{3bf}\u{3b4}\u{3c5}\u{3c3}\u{3c3}\u{3b5}\u{3c5}\u{3c2}",
        11_140_088_625_042_210_881,
    ),
    // A no-break space and an ordinary one, sharing a value. This pins the one
    // residual Unicode dependency the module doc admits: `split_whitespace`
    // reads the `White_Space` property, so a release that changed it for
    // `U+00A0` fails *here*, loudly, instead of silently reassigning every
    // universe named with one.
    ("a\u{a0}b", 16_609_814_704_669_528_739),
    ("a b", 16_609_814_704_669_528_739),
    // A digit typed into the *name* box is a name. This is the two-box
    // decision made checkable: with one box, `"7"` would have to be guessed at,
    // and guessing means deleting `letters_in_the_seed_box_are_refused_not_hashed`
    // one crate down — the test that forbids "be friendly, hash whatever was
    // typed".
    ("7", 7_059_366_201_351_256_473),
    // A NUL is part of the name rather than the end of it.
    ("a\u{0}b", 15_243_831_238_367_000_186),
];

#[cfg(test)]
mod golden {
    use super::{GOLDEN, seed_from_phrase};
    use std::collections::BTreeSet;

    #[test]
    fn the_mapping_from_a_phrase_to_a_seed_is_pinned() {
        for &(phrase, expected) in GOLDEN {
            assert_eq!(
                seed_from_phrase(phrase),
                Some(expected),
                "the universe named {phrase:?} moved — every universe anyone \
                 has written down under that name is now a different one. This \
                 is a physics change, not a stale expectation: regenerate only \
                 if the mapping was *meant* to change"
            );
        }
    }

    /// Bookkeeping on the table itself.
    ///
    /// **Both counts are literals, and neither is decoration.** The table's
    /// whole method is that some rows deliberately share a value and others
    /// deliberately do not, so a copy-pasted row is invisible to the test above
    /// — it would simply assert a true thing twice. The row count catches a row
    /// being dropped; the distinct count catches one being duplicated, and
    /// catches a normalisation change that collapses two rows onto each other
    /// even when both literals still appear somewhere in the table.
    #[test]
    fn the_golden_still_has_the_shape_it_was_written_with() {
        assert_eq!(GOLDEN.len(), 15, "a golden row was added or removed");

        let distinct: BTreeSet<u64> = GOLDEN.iter().map(|&(_, value)| value).collect();
        assert_eq!(
            distinct.len(),
            11,
            "the golden has {} distinct values across {} rows — four of the \
             rows share a value on purpose (three spellings of `emily`, two of \
             `emily rose`, and the two spellings of a space); any other \
             coincidence is a normalisation change",
            distinct.len(),
            GOLDEN.len()
        );

        let phrases: BTreeSet<&str> = GOLDEN.iter().map(|&(phrase, _)| phrase).collect();
        assert_eq!(
            phrases.len(),
            GOLDEN.len(),
            "a phrase appears twice in the golden"
        );
    }
}
