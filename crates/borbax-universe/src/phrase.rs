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
//! two reasons, of which **one** is load-bearing: the viewer has to show the
//! number so it can be written down, and §6 makes the seed the shareable thing.
//!
//! The other reason usually offered — that `Universe::from_phrase` would be a
//! second `Universe::generate` call site under an `xtask` guard — is **false,
//! and is recorded here rather than deleted**, because it reads plausible.
//! `check_the_viewer_seam_holds` walks `crates/borbax-ui/src` and nothing else;
//! this crate's own `generate` call sites are counted by nothing. Measured by a
//! reviewer.
//!
//! **The mapping is frozen.** A phrase is now a human-facing address for a
//! universe and, unlike §6's `U-…@N`, it carries no version — so there is
//! nothing to bump on the day it changes, and every "Emily" ever written down
//! would silently resolve elsewhere. If the mapping must ever change, add a
//! function beside this one rather than editing this one, and let the old
//! phrases keep their universes. The mitigation that already exists is that the
//! viewer shows the *number*: the phrase is the mnemonic, the seed is the
//! receipt.
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
//! move between compiler releases — Rust reached `UNICODE_VERSION` 17.0.0 from
//! 16.0.0 during the 1.9x series — so a Unicode-aware fold would let an
//! innocent toolchain bump silently reassign every universe anyone had already
//! named. ASCII folding cannot move.
//!
//! **The hazard is prospective, and the first wording overstated it.** It said
//! the tables moved "on the pinned toolchain", which a reviewer measured false:
//! `.prototools` has a single commit and has always read `1.97.1`, whose
//! `UNICODE_VERSION` is 17.0.0, so this project's pin has never crossed that
//! boundary. That is a reason to close the hazard now rather than a reason to
//! think it already fired. The accepted cost is that `É` and `é` are
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
//! ruling, and the reason first given for it was too comfortable: "a far more
//! stable table than the case mappings". Measured against unicode.org by a
//! reviewer, Unicode's stability policy covers **`Pattern_White_Space`** and
//! grants `White_Space` no guarantee at all — and `White_Space` has moved in
//! practice, `U+180E` losing it in Unicode 6.3.
//!
//! So the whole property is pinned rather than argued about.
//! `the_set_of_characters_treated_as_whitespace_is_pinned` enumerates every
//! scalar value and asserts the exact 25-codepoint set, for about 10 ms. The
//! `U+00A0` golden row remains as a second, narrower check — it covered 2 of
//! the 25, which is why it was not enough on its own.
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
/// this is a chaining construction rather than a sum, so a prefix does not
/// share a seed with what follows it by *construction* — as it would under any
/// commutative combine.
///
/// **"Cannot" would be the wrong word and an earlier draft used it.** `absorb`
/// truncates a four-word Philox block to one word, so for a fixed byte it is a
/// pseudorandom *function*, not a bijection: collisions are improbable, not
/// impossible. `Stream::sub`'s own doc insists those are different things.
/// The honest number is P ≈ N²/2⁶⁵ — about **1.4e-8** across a million named
/// universes.
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
/// The hazard, stated **structurally**, because the tempting version of it is
/// weak enough to lose the argument. It is *not* "the accumulator might happen
/// to equal the clock reading" — that is 2⁻⁶⁴ per byte per click, and anyone
/// who wants the guard gone can say so. It is that
/// `Stream::new(acc, Domain::Hash, u64::from(byte))` at `byte == 0` claims the
/// coordinate `(Hash, index 0, root)`, which the viewer's "surprise me" button
/// **already owns** — an ownership overlap by construction, independent of any
/// value.
///
/// `sub` offsets its coordinate by one, so `counter[2]` is 1 here and 0 there,
/// and the two constructions cannot alias for *any* byte or *any* accumulator.
///
/// **This call site owns `(Domain::Hash, index 0..=255, sub 0)`**, recorded
/// because `Stream::sub`'s own doc requires it: "a sub-index range is owned by
/// exactly one call site". A future `Domain::Hash` consumer — a fold or species
/// cache key is the obvious one, and `Domain::Hash`'s doc names exactly that
/// use — must take a different `sub`, or its own `Domain`.
fn absorb(acc: u64, byte: u8) -> u64 {
    Stream::sub(acc, Domain::Hash, u64::from(byte), 0).next_u64()
}

/// Trim, collapse internal whitespace to single spaces, and fold ASCII case.
///
/// Written as an explicit fold rather than `split_whitespace().collect::<Vec<_>>()
/// .join(" ")` because that allocates a `Vec` of slices as well as the string.
/// Measured with a counting allocator: **1 allocation against 2–4**, growing
/// with word count only through `Vec` growth. (An earlier version of this
/// sentence said "no intermediate allocation per word". There is no allocation
/// per word in either version — `collect::<Vec<&str>>` copies no string data —
/// so the choice was right and the reason given for it was not.) The
/// `with_capacity` is never exceeded, since normalising only removes bytes. And
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
        // The discriminator is the assertion above rather than an extra one: an
        // infallible `-> u64` returns *some* universe for the empty string,
        // which is a real element count on screen attributed to a name nobody
        // typed. An earlier version repeated `is_none()` here under a comment
        // calling it "the discriminator", which was the same assertion twice.
        //
        // What the repeat did not cover, and this does: a phrase made only of
        // *non-ASCII* whitespace still names nothing.
        assert_eq!(seed_from_phrase("\u{2003}\u{3000}"), None);
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

    /// §13.1's hasher ban still resolves — a liveness anchor, not a test of us.
    ///
    /// **The cross-check in `xtask` compares only the last path segment**, so a
    /// wrong *module* prefix in `clippy.toml` leaves clippy resolving nothing
    /// while `the_two_type_ban_lists_cover_the_same_types` stays green. A
    /// determinism reviewer measured exactly that: changing one entry to
    /// `std::collections::hash_map::NotReallyThere::DefaultHasher` and planting
    /// a real `DefaultHasher` gave **zero** clippy output and exit 0, with the
    /// cross-check green throughout. Clippy does not report an unresolvable
    /// `disallowed-types` path, and `-D warnings` cannot promote a diagnostic
    /// that is never emitted.
    ///
    /// The older method list does not have this hole, because its cross-check
    /// compares the *function* name and so breaks on a typo anywhere after the
    /// first `::`. The type list duplicates only the leaf — the one part a
    /// prefix typo leaves intact.
    ///
    /// `#[expect]` closes it, and closes more than it: an expectation is
    /// unfulfilled if the lint fails to fire for **any** reason — an
    /// unresolvable path, a level dropped from `deny` to `allow`, or the lint
    /// leaving clippy's default group. It lives inside `#[cfg(test)]` so
    /// `xtask`'s textual scan skips it, while clippy still lints it under
    /// `--all-targets`.
    #[test]
    #[expect(
        clippy::disallowed_types,
        reason = "this test exists to make the lint fire; if it stops firing, the \
                  expectation is unfulfilled and the build fails, which is the point"
    )]
    fn the_hasher_ban_still_resolves() {
        // Named in a return type rather than a `let`: the type has to appear in
        // a resolved path for the lint to fire, and a binding to a `_`-prefixed
        // variable with no side effect is itself denied here.
        fn anchor() -> Option<std::collections::hash_map::DefaultHasher> {
            None
        }
        assert!(anchor().is_none());
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

    /// Fifty names a child would plausibly type.
    ///
    /// Shared by the collision test and the avalanche test so the two measure
    /// the same population.
    fn corpus() -> Vec<String> {
        let mut names = Vec::new();
        for first in [
            "emily", "olivia", "amelia", "isla", "ava", "ivy", "freya", "lily", "mia", "willow",
        ] {
            names.push(first.to_owned());
            for second in ["rose", "grace", "may", "jane"] {
                names.push(format!("{first} {second}"));
            }
        }
        names
    }

    /// Names a child would plausibly type land on distinct universes.
    ///
    /// A weak-mixing implementation — folding bytes with `+` or `^` — collides
    /// short similar names constantly. Measured on this corpus: zero
    /// collisions, and the bar is exact rather than a fraction because any
    /// collision at this size is a defect, not noise.
    #[test]
    fn a_corpus_of_ordinary_names_does_not_collide() {
        let names = corpus();
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
    /// **Every bar here was set after measuring, and two earlier versions of
    /// this test were wrong in ways worth recording, because both looked
    /// rigorous.**
    ///
    /// The first asserted that every trial landed in `24..=40` — a bar chosen
    /// before any number existed, which failed immediately on a perfectly good
    /// mixer at 23. The per-trial minimum is an *order statistic*, not an
    /// estimate of the mean.
    ///
    /// The second fixed that and then mis-stated its own sample. It flipped
    /// `1 << (bit % 5)` over `0..8`, which repeats shifts 0, 1 and 2 — so it
    /// pushed 40 values drawn from only **25 distinct** flips, and every
    /// statistic derived from `n = 40` was wrong: the standard error was 0.84
    /// rather than the claimed `4/sqrt(40)` = 0.63, the band was ±2.4 standard
    /// errors rather than three, and the quoted `P(X <= 23)` was 1.4% against a
    /// true 1.638%. Its self-check `moved_bits.len() >= 30` reported "the corpus
    /// filter selected almost nothing" while counting *pushes*, not flips —
    /// and all three of its `continue` arms were unreachable, so the filter it
    /// guarded could not filter. Four review lanes found this independently.
    ///
    /// This version flips **every** bit of **every** byte across the whole
    /// [`corpus`], keeps only flips that survive normalisation unchanged (a
    /// flip that merely makes a capital is not a one-bit change to the *input*
    /// of `absorb`), and deduplicates. Measured: **2590 distinct flips**, mean
    /// **31.9656**, min 17, max 44.
    ///
    /// The three bars, and what each is for:
    ///
    /// - `moved.len() >= 2000` now counts distinct flips, so it is a real
    ///   self-check rather than a restatement of the loop bounds.
    /// - the **mean** gets `31.5..=32.5`, which at `SE = 4/sqrt(2590)` = 0.0786
    ///   is ±6.4 standard errors. Measured z = −0.44.
    /// - each **trial** gets `7..=57`, about 6.25 standard deviations: exact
    ///   `P(outside)` is 9.03e-12 per trial, so 2.3e-8 across all 2590 — a
    ///   hundredfold safer than the previous `12..=52` at this sample size,
    ///   while still failing on the first sample for an identity (0 bits moved)
    ///   or an inverter (64).
    ///
    /// **Two-sided on purpose.** A floor alone passes for the mixer's opposite:
    /// flipping *every* bit is as broken as flipping none.
    #[test]
    fn one_changed_byte_moves_about_half_the_output_bits() {
        let mut seen = BTreeSet::new();
        let mut moved = Vec::new();
        for base in corpus() {
            let unflipped = seed(&base);
            for position in 0..base.len() {
                for shift in 0..8_u8 {
                    let mut bytes = base.as_bytes().to_vec();
                    let Some(target) = bytes.get_mut(position) else {
                        continue;
                    };
                    *target ^= 1 << shift;
                    let Ok(flipped) = std::str::from_utf8(&bytes) else {
                        continue;
                    };
                    // Keep only flips `normalise` leaves alone: a flip that
                    // merely makes a capital folds straight back, so counting
                    // it would measure a change that never reached `absorb`.
                    // The high bits are not uniformly rejected — measured, they
                    // contribute 465 of the 2590 survivors, because e.g. bit 5
                    // of a space gives NUL, which is not whitespace.
                    if super::normalise(flipped).as_bytes() != bytes.as_slice() {
                        continue;
                    }
                    if seen.insert(flipped.to_owned()) {
                        moved.push((unflipped ^ seed(flipped)).count_ones());
                    }
                }
            }
        }

        assert!(
            moved.len() >= 2000,
            "only {} DISTINCT flips — measured 2590, so the filter has started \
             rejecting almost everything and every bound below is vacuous",
            moved.len()
        );
        for &bits in &moved {
            assert!(
                (7..=57).contains(&bits),
                "a single flip moved {bits} of 64 bits"
            );
        }

        // `u32` throughout, so both conversions below are `f64::from` — exact,
        // and needing no `as`. The sum is bounded by 64 * 2590 = 165_760.
        let total: u32 = moved.iter().sum();
        let count = u32::try_from(moved.len()).unwrap_or(u32::MAX);
        let mean = f64::from(total) / f64::from(count);
        assert!(
            (31.5..=32.5).contains(&mean),
            "mean bits moved {mean:.4} over {count} distinct flips"
        );
    }

    /// The seed is a full 64-bit value and is deliberately not capped.
    ///
    /// `RANDOM_SEED_CEILING` in the viewer exists because a *suggested* seed has
    /// to be copyable off a screen by hand. With a phrase, the phrase is the
    /// memorable artefact and the number is a receipt, so capping buys nothing
    /// and costs collisions — at 10^6 a thousand names collide about 40% of the
    /// time.
    ///
    /// **This test asserted on `GOLDEN`'s literals and never called the
    /// function, and two review lanes measured it green over a planted
    /// `Some(acc % 1_000_000)`.** A test named for a property of the seed that
    /// can only fail when someone edits a table is not a guard on the seed. It
    /// now reads the live mapping over the whole corpus.
    #[test]
    fn the_seed_is_not_squeezed_into_six_digits() {
        let mut over_six_digits = 0_usize;
        let mut over_32_bits = 0_usize;
        for name in corpus() {
            let value = seed(&name);
            if value > 1_000_000 {
                over_six_digits += 1;
            }
            if value > u64::from(u32::MAX) {
                over_32_bits += 1;
            }
        }
        // Both bars are far below the ~50 and ~49.99 a uniform u64 gives on 50
        // names, and far above what any cap permits: a modulo into six digits
        // makes both counts exactly 0.
        assert!(
            over_six_digits >= 45,
            "only {over_six_digits} of 50 names exceed the viewer's suggestion ceiling"
        );
        assert!(
            over_32_bits >= 45,
            "only {over_32_bits} of 50 names exceed 32 bits"
        );
    }

    /// No prefix of a phrase is ever enough to decide the universe.
    ///
    /// **The whole suite was blind to a length cap, measured by a review lane:
    /// planting `.take(16)` on the byte loop passed every one of the
    /// 444 tests the workspace then held** while silently collapsing `"elizabeth alexandra mary"` and
    /// `"elizabeth alexandra rose"` onto one universe. Caps at 24 and 32 passed
    /// the golden too; the longest normalised phrase anywhere in the suite was
    /// 16 bytes.
    ///
    /// It is not a contrived mutation. It is the shape of a future "do not hash
    /// a ten-megabyte paste" guard, and `egui::TextEdit` has a `char_limit`
    /// builder someone will reach for on the panel side, where no test in this
    /// crate can see it.
    ///
    /// The discriminator is a shared *prefix* longer than any cap anyone would
    /// write: under `take(n)` these agree for every `n` up to the shared
    /// prefix's length.
    #[test]
    fn no_prefix_of_a_phrase_is_enough_to_decide_the_universe() {
        let stem = "elizabeth alexandra ";
        let distinct: BTreeSet<u64> = ["mary", "rose", "jane", "may"]
            .iter()
            .map(|tail| seed(&format!("{stem}{tail}")))
            .collect();
        assert_eq!(
            distinct.len(),
            4,
            "the tail of a long phrase was discarded — a cap at or below {} bytes",
            stem.len() + 4
        );

        // And at a length nothing plausible would cap: 512 bytes differing only
        // in the last one.
        let long_a = "a".repeat(511) + "x";
        let long_b = "a".repeat(511) + "y";
        assert_ne!(seed(&long_a), seed(&long_b));
    }

    /// The set of characters `split_whitespace` folds is pinned outright.
    ///
    /// **The `U+00A0` golden row covers 2 of the 25 codepoints, not the
    /// hazard.** `split_whitespace` reads Unicode's `White_Space` property, and
    /// — measured by a review lane against unicode.org — `White_Space` is *not*
    /// covered by Unicode's stability policy. `Pattern_White_Space` is; this is
    /// not it. `White_Space` has moved in practice: `U+180E` lost it in Unicode
    /// 6.3.
    ///
    /// So the whole table is pinned rather than one member of it. A release
    /// that adds or removes any codepoint fails *here*, with a message saying
    /// what moved, instead of silently reassigning every universe named with a
    /// phrase containing it.
    ///
    /// Enumerating the full scalar range costs about 10 ms.
    #[test]
    fn the_set_of_characters_treated_as_whitespace_is_pinned() {
        let found: Vec<u32> = (0..=0x0010_FFFF_u32)
            .filter_map(char::from_u32)
            .filter(|c| c.is_whitespace())
            .map(u32::from)
            .collect();
        let expected: &[u32] = &[
            0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x20, 0x85, 0xA0, 0x1680, 0x2000, 0x2001, 0x2002, 0x2003,
            0x2004, 0x2005, 0x2006, 0x2007, 0x2008, 0x2009, 0x200A, 0x2028, 0x2029, 0x202F, 0x205F,
            0x3000,
        ];
        assert_eq!(
            found, expected,
            "Unicode's White_Space property moved on this toolchain. Every phrase \
             containing an affected character now names a different universe. This \
             is a toolchain change, not a code change: do NOT regenerate the golden \
             without deciding that the reassignment is acceptable"
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
/// changes nothing *within* ASCII. Both are therefore invisible to an
/// ASCII-only table, and both are caught below.
///
/// (An earlier version of that sentence said `to_lowercase` "changes nothing
/// outside ASCII either", which is the opposite of the property the Greek rows
/// exist for. Three review lanes caught it.)
///
/// **What may go in this table.** Phrases only — a name, a word, a spelling
/// being pinned. Never a real element or compound name, however tempting a
/// themed row looks: this is a literal string table checked into a chemistry
/// crate, and a `("Carbon", ..)` row would put a real-chemistry name inside
/// `borbax-universe` where G1's data-file scan cannot see it.
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
    // so the first row's literal stops matching.
    //
    // They are the only *golden row* that moves under that swap, which is what
    // makes them worth their place. They are NOT the only test that catches it:
    // `accents_are_part_of_the_name` catches it independently and more cheaply,
    // via `CAFÉ` vs `café`. An earlier version of this comment claimed nothing
    // else caught it — false, measured by three lanes, and the kind of claim
    // that gets the cheaper guard deleted as redundant.
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
    /// catches a normalisation change that collapses two rows onto each other.
    ///
    /// **The distinct set is built from `seed_from_phrase`, not from the
    /// literals, and that is a repair rather than a preference.** Built from
    /// the literals it could not see a normalisation change at all — two review
    /// lanes measured it green under a `str::to_lowercase` swap, which is
    /// precisely the collapse its own comment claimed to catch. Reading the
    /// live mapping makes the claim true.
    #[test]
    fn the_golden_still_has_the_shape_it_was_written_with() {
        assert_eq!(GOLDEN.len(), 15, "a golden row was added or removed");

        let distinct: BTreeSet<u64> = GOLDEN
            .iter()
            .filter_map(|&(phrase, _)| seed_from_phrase(phrase))
            .collect();
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
