//! Procedural element naming, with a real-name blocklist (spec §5, G2).
//!
//! The blocklist below exists **solely to exclude**. It is never read as
//! data, never used to derive a property, and carries no information beyond
//! "this string is taken". Without it a generated universe could mint an
//! element called `He` or `carbon`, which would falsely imply a mapping onto
//! real chemistry — exactly the impression the whole design exists to avoid.

use borbax_rng::Stream;

/// Real element symbols. Exclusion list only (G2).
/// Name must match the string `xtask` greps for — see Task 1 Step 3.
const REAL_ELEMENT_SYMBOLS: &[&str] = &[
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// Real element names plus common real chemical and biological terms, all
/// lowercase. Exclusion list only (G2).
const REAL_WORDS: &[&str] = &[
    "hydrogen",
    "helium",
    "lithium",
    "beryllium",
    "boron",
    "carbon",
    "nitrogen",
    "oxygen",
    "fluorine",
    "neon",
    "sodium",
    "magnesium",
    "aluminium",
    "aluminum",
    "silicon",
    "phosphorus",
    "sulfur",
    "sulphur",
    "chlorine",
    "argon",
    "potassium",
    "calcium",
    "scandium",
    "titanium",
    "vanadium",
    "chromium",
    "manganese",
    "iron",
    "cobalt",
    "nickel",
    "copper",
    "zinc",
    "gallium",
    "germanium",
    "arsenic",
    "selenium",
    "bromine",
    "krypton",
    "rubidium",
    "strontium",
    "yttrium",
    "zirconium",
    "niobium",
    "molybdenum",
    "technetium",
    "ruthenium",
    "rhodium",
    "palladium",
    "silver",
    "cadmium",
    "indium",
    "tin",
    "antimony",
    "tellurium",
    "iodine",
    "xenon",
    "caesium",
    "cesium",
    "barium",
    "lanthanum",
    "cerium",
    "praseodymium",
    "neodymium",
    "promethium",
    "samarium",
    "europium",
    "gadolinium",
    "terbium",
    "dysprosium",
    "holmium",
    "erbium",
    "thulium",
    "ytterbium",
    "lutetium",
    "hafnium",
    "tantalum",
    "tungsten",
    "rhenium",
    "osmium",
    "iridium",
    "platinum",
    "gold",
    "mercury",
    "thallium",
    "lead",
    "bismuth",
    "polonium",
    "astatine",
    "radon",
    "francium",
    "radium",
    "actinium",
    "thorium",
    "protactinium",
    "uranium",
    "neptunium",
    "plutonium",
    "americium",
    "curium",
    "berkelium",
    "californium",
    "einsteinium",
    "fermium",
    "mendelevium",
    "nobelium",
    "lawrencium",
    "rutherfordium",
    "dubnium",
    "seaborgium",
    "bohrium",
    "hassium",
    "meitnerium",
    "darmstadtium",
    "roentgenium",
    "copernicium",
    "nihonium",
    "flerovium",
    "moscovium",
    "livermorium",
    "tennessine",
    "oganesson",
    // Common real terms that would imply a mapping even without a symbol clash.
    "water",
    "protein",
    "enzyme",
    "peptide",
    "lipid",
    "sugar",
    "glucose",
    "amine",
    "acid",
    "alcohol",
    "ester",
    "ketone",
    "methane",
    "ammonia",
    "benzene",
    "cellulose",
    "chitin",
];

/// True if either the symbol or the name collides with something real.
#[must_use]
pub fn is_real(symbol: &str, name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    REAL_ELEMENT_SYMBOLS
        .iter()
        .any(|s| s.eq_ignore_ascii_case(symbol))
        || REAL_WORDS.contains(&n.as_str())
}

// Phonotactics chosen to be pronounceable and to sound like they belong to
// one language. Deliberately unlike real element phonology.
const ONSETS: &[&str] = &[
    "b", "d", "f", "g", "k", "l", "m", "n", "p", "r", "s", "t", "v", "z", "br", "dr", "fl", "gl",
    "kr", "pl", "st", "tr", "vr", "zh", "th", "sk",
];
const NUCLEI: &[&str] = &["a", "e", "i", "o", "u", "ae", "ia", "au", "oi", "ei"];
const CODAS: &[&str] = &["", "", "l", "n", "r", "s", "x", "th", "rn", "st", "sk"];
const SUFFIXES: &[&str] = &["ium", "ex", "on", "ar", "yl", "is", "or", "ax"];

/// Draw one element of `xs` uniformly.
#[expect(
    clippy::indexing_slicing,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    reason = "`next_range(n)` returns `v % n`, so the index is strictly < xs.len() \
              for every non-empty slice, and every slice passed here is a non-empty \
              `const`. That bound is what makes the u64 -> usize narrowing exact on \
              a 32-bit target too — the value never exceeds a slice length"
)]
fn pick<'a>(rng: &mut Stream, xs: &[&'a str]) -> &'a str {
    xs[rng.next_range(xs.len() as u64) as usize]
}

/// Which path produced a symbol.
///
/// The `Fallback` arm is a last resort and is **not well-formed past `Q9`** —
/// three characters, and a name carrying a digit. It exists so a test meaning
/// to exercise the generator can *say so*, rather than sniffing for a `Q`
/// prefix that only works because `ONSETS` happens to contain no `q`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Provenance {
    /// From the phonotactic grammar.
    Generated,
    /// From the exhaustion fallback.
    Fallback,
}

/// Mint one element name, unique against `taken` and clear of the blocklist.
///
/// `taken` accumulates symbols already used in this universe and is mutated.
///
/// **The symbol space is 165, not the 702 an earlier version of this comment
/// implied by "26x27".** Enumerated exhaustively over the grammar: 7 one-letter
/// and 158 two-letter symbols survive the blocklist, because only the 14
/// distinct initials of `ONSETS` can start one and the second letter must
/// appear in the same generated name. Against a table that reaches 120 that is
/// 1.4x headroom. Filling 120 elements takes **480 attempts at the median**
/// (min 333, max 704, zero fallbacks over 2000 universes); the 120th element
/// alone needs a **mean of 11.5** attempts, p90 26, p95 33, p99 52, max 93.
///
/// An earlier version said "~38 retries" without saying which statistic that
/// was: it sits between p95 and p99, reached in 3.5% of universes, and reading
/// it as typical understates the headroom this paragraph is arguing for.
///
/// **The budget is a cliff, not a gradient, and it is worth naming.** Below the
/// symbol space a mint costs ~5.4 us at position 119; at 165 every mint burns
/// all 10 000 attempts at **6.07 ms**, a 1100x step. `n_elements` is drawn
/// 60..=120 so this is unreachable, and
/// `the_generator_covers_the_largest_table_without_falling_back` pins the
/// margin — but raising the draw to 200 without widening `ONSETS` turns a
/// 177 us universe into ~0.5 s, and no test would fail.
#[expect(
    clippy::indexing_slicing,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    reason = "`letters` is non-empty because every `ONSETS` entry is, so `letters[0]` \
              is in bounds; `idx` is `1 + next_range(len - 1)`, hence in `1..len`. \
              That bound is what makes the u64 -> usize narrowing exact on a 32-bit \
              target too — the value never exceeds a name length"
)]
pub fn mint(rng: &mut Stream, taken: &mut Vec<String>) -> (String, String, Provenance) {
    for _ in 0..10_000 {
        let name = format!(
            "{}{}{}{}",
            pick(rng, ONSETS),
            pick(rng, NUCLEI),
            pick(rng, CODAS),
            pick(rng, SUFFIXES)
        );

        // Symbol: first letter, plus one later letter from the name itself so
        // the symbol visibly belongs to the word.
        let letters: Vec<char> = name.chars().collect();
        let mut symbol = String::new();
        symbol.push(letters[0].to_ascii_uppercase());
        if letters.len() > 2 && rng.next_f64() < 0.85 {
            let idx = 1 + rng.next_range((letters.len() - 1) as u64) as usize;
            symbol.push(letters[idx]);
        }

        if !is_real(&symbol, &name) && !taken.contains(&symbol) {
            taken.push(symbol.clone());
            return (symbol, name, Provenance::Generated);
        }
    }
    // Reachable at 166 mints, demonstrated — not "unreachable in practice" as
    // an earlier comment said. The fallback must therefore preserve the
    // uniqueness postcondition this function's own doc promises: a draft used
    // `taken.len() % 10`, giving ten symbols and an unconditional push, which
    // made the fallback the *only* path that could emit a duplicate. It did.
    //
    // Dropping the `% 10` makes it unique for as long as `taken` grows, which
    // is the postcondition. It does **not** make it well-formed: past `Q9` the
    // symbol is three characters and the name carries a digit, so both
    // `symbols_are_well_formed` and the 1..=2 length contract fail there.
    //
    // That is stated rather than fixed, deliberately. The real guarantee is the
    // margin — `the_generator_covers_the_largest_table_without_falling_back` —
    // and a fallback elaborate enough to stay well-formed would look like a
    // guarantee while still being a last resort. If a universe ever needs more
    // than 165 elements, widen `ONSETS` and re-measure the space; do not extend
    // this.
    // Bind the index once: a draft read `taken.len()` before the push for the
    // symbol and after it for the name, so they disagreed by one in the path
    // whose whole justification is "unique and traceable".
    let idx = taken.len();
    let symbol = format!("Q{idx}");
    taken.push(symbol.clone());
    (symbol, format!("quorium{idx}"), Provenance::Fallback)
}

#[cfg(test)]
mod tests {
    use super::*;
    use borbax_rng::{Domain, Stream};

    #[test]
    fn names_are_deterministic() {
        let mut a = Stream::new(1, Domain::Naming, 0);
        let mut b = Stream::new(1, Domain::Naming, 0);
        let (mut ta, mut tb) = (Vec::new(), Vec::new());
        assert_eq!(mint(&mut a, &mut ta), mint(&mut b, &mut tb));
    }

    /// G2 (spec §5). The single most important test in this file — and a draft
    /// of it was **96.7% vacuous**.
    ///
    /// It drew 5 000 names against one shared `taken`. The space is 165, so
    /// 4 835 of those iterations exercised the `Q{n}` fallback rather than the
    /// phonotactic generator the test exists to check. A fresh `taken` per draw
    /// keeps every iteration on the generator.
    #[test]
    fn never_mints_a_real_element() {
        let mut s = Stream::new(99, Domain::Naming, 0);
        for _ in 0..5_000 {
            let mut taken = Vec::new();
            let (sym, name, _) = mint(&mut s, &mut taken);
            assert!(
                !is_real(&sym, &name),
                "minted a real element: {sym} / {name}"
            );
        }
    }

    /// The independent oracle both tests below use. **Never read `REAL_WORDS`
    /// to build this** — it is the blocklist under test, and its whole defect
    /// was a silent gap in exactly this list. Built from atomic number, zipped
    /// against `REAL_ELEMENT_SYMBOLS` (already complete and in atomic-number
    /// order), so this list's own completeness is checked structurally rather
    /// than assumed.
    ///
    /// Verified live on `main` at `cb3a28c`, before this fix: `REAL_WORDS`
    /// covered only 78 of these 118 — it jumped from `cerium` (Z=58) straight
    /// to `tungsten` (Z=74), dropping the whole Z=59..=73 lanthanide run, then
    /// dropped Z=75..=77, Z=91 and everything from Z=95 onward. 43 names
    /// missing in total, one of them `thulium`.
    const REAL_ELEMENT_NAMES: [&str; 118] = [
        "hydrogen",
        "helium",
        "lithium",
        "beryllium",
        "boron",
        "carbon",
        "nitrogen",
        "oxygen",
        "fluorine",
        "neon",
        "sodium",
        "magnesium",
        "aluminium",
        "silicon",
        "phosphorus",
        "sulfur",
        "chlorine",
        "argon",
        "potassium",
        "calcium",
        "scandium",
        "titanium",
        "vanadium",
        "chromium",
        "manganese",
        "iron",
        "cobalt",
        "nickel",
        "copper",
        "zinc",
        "gallium",
        "germanium",
        "arsenic",
        "selenium",
        "bromine",
        "krypton",
        "rubidium",
        "strontium",
        "yttrium",
        "zirconium",
        "niobium",
        "molybdenum",
        "technetium",
        "ruthenium",
        "rhodium",
        "palladium",
        "silver",
        "cadmium",
        "indium",
        "tin",
        "antimony",
        "tellurium",
        "iodine",
        "xenon",
        "caesium",
        "barium",
        "lanthanum",
        "cerium",
        "praseodymium",
        "neodymium",
        "promethium",
        "samarium",
        "europium",
        "gadolinium",
        "terbium",
        "dysprosium",
        "holmium",
        "erbium",
        "thulium",
        "ytterbium",
        "lutetium",
        "hafnium",
        "tantalum",
        "tungsten",
        "rhenium",
        "osmium",
        "iridium",
        "platinum",
        "gold",
        "mercury",
        "thallium",
        "lead",
        "bismuth",
        "polonium",
        "astatine",
        "radon",
        "francium",
        "radium",
        "actinium",
        "thorium",
        "protactinium",
        "uranium",
        "neptunium",
        "plutonium",
        "americium",
        "curium",
        "berkelium",
        "californium",
        "einsteinium",
        "fermium",
        "mendelevium",
        "nobelium",
        "lawrencium",
        "rutherfordium",
        "dubnium",
        "seaborgium",
        "bohrium",
        "hassium",
        "meitnerium",
        "darmstadtium",
        "roentgenium",
        "copernicium",
        "nihonium",
        "flerovium",
        "moscovium",
        "livermorium",
        "tennessine",
        "oganesson",
    ];

    /// **`REAL_WORDS` itself can be incomplete in a way `never_mints_a_real_element`
    /// cannot see, because that test's oracle *is* `REAL_WORDS`.** A name missing
    /// from the list is never rejected by `is_real`, however many draws are taken.
    /// This test checks the blocklist's own coverage against `REAL_ELEMENT_NAMES`,
    /// never against itself.
    #[test]
    fn the_blocklist_names_every_real_element() {
        assert_eq!(
            REAL_ELEMENT_NAMES.len(),
            REAL_ELEMENT_SYMBOLS.len(),
            "the reference list above must cover the same 118 elements \
             REAL_ELEMENT_SYMBOLS does, or the zip below is comparing the wrong pairs"
        );
        for (symbol, name) in REAL_ELEMENT_SYMBOLS.iter().zip(REAL_ELEMENT_NAMES) {
            assert!(
                REAL_WORDS.contains(&name),
                "REAL_WORDS is missing '{name}' (element {symbol}) — a generated \
                 element could be minted with that name and pass is_real() falsely"
            );
        }
    }

    /// The empirical companion to `the_blocklist_names_every_real_element`:
    /// sweep universe seeds, not just draws within one, and confirm no minted
    /// name collides with a real one under the actual generator.
    ///
    /// **This must check against `REAL_ELEMENT_NAMES`, never against
    /// `is_real`.** A first draft of this test asserted `!is_real(&sym, &name)`
    /// and passed cleanly even with the bug this test exists to catch, still
    /// present. The reason: `mint` already filters every candidate through
    /// `is_real` before returning it, so *nothing `mint` ever returns can fail
    /// an `is_real`-based check* — the assertion was checking the generator
    /// agreed with itself, not with reality. Caught by hand-tracing seed 1631
    /// before this landed, not by the test — which is exactly why the
    /// independent-oracle version below is the one that ships.
    ///
    /// **Neither pinned universe digest can see this class of regression
    /// either.** `symbol`/`name` are mixed into `the_universe_digest_is_pinned`
    /// (`element.rs`), but that digest sweeps seeds `0..64` only. The specific
    /// defect this test was added for first reproduces at seed 1631, at 43 of
    /// the first 20 000 seeds (0.215%). A digest test would have reported
    /// "unchanged", correctly, while this many universes shipped a
    /// falsely-real element name.
    #[test]
    fn no_seed_in_a_large_sweep_mints_a_real_element_name() {
        let mut examined: u64 = 0;
        for seed in 0..20_000u64 {
            let mut s = Stream::new(seed, Domain::Naming, 0);
            let mut taken = Vec::new();
            for _ in 0..120 {
                let (_, name, _) = mint(&mut s, &mut taken);
                examined += 1;
                let lower = name.to_ascii_lowercase();
                assert!(
                    !REAL_ELEMENT_NAMES.contains(&lower.as_str()),
                    "seed {seed} minted a real element name: {name}"
                );
            }
        }
        assert_eq!(
            examined,
            20_000 * 120,
            "examined count must equal every mint actually attempted"
        );
    }

    /// The margin between the symbol space and the largest table, pinned.
    ///
    /// 165 against 120 is 1.4x, not the "far larger" `mint`'s doc claimed.
    ///
    /// **Assert provenance, not uniqueness.** A draft asserted
    /// `unique.len() == taken.len()`, which **both** `mint` paths preserve —
    /// the generator only pushes unseen symbols and the fallback pushes
    /// `Q{taken.len()}`, unique because the length strictly increases. So the
    /// assertion was invariant: measured, cutting the grammar to a space of 7
    /// leaves it passing while 158 of 165 symbols come from the last resort,
    /// under a failure message reading "the generator ran out before 165".
    ///
    /// That is the same vacuity the G2 fix removed — a test counting fallback
    /// draws as generator draws — reintroduced by the fix for it. `mint`
    /// therefore returns which path produced the symbol.
    #[test]
    fn the_generator_covers_the_largest_table_without_falling_back() {
        let mut s = Stream::new(5, Domain::Naming, 0);
        let mut taken = Vec::new();
        for i in 0..165 {
            let (sym, _, provenance) = mint(&mut s, &mut taken);
            assert_eq!(
                provenance,
                Provenance::Generated,
                "the generator ran out at {i}, falling back to {sym}"
            );
        }
    }

    /// **120, not 200 — the symbol space is 165 and 200 cannot pass.**
    /// Exhaustive enumeration of the grammar (26 onsets x 10 nuclei x 11 codas
    /// x 8 suffixes = 20 800 names) yields 7 one-character and 158
    /// two-character symbols after the blocklist. `mint`'s own doc claimed
    /// "26x27", i.e. 702, a 4.3x overstatement. 120 is the largest table this
    /// task builds; `the_generator_covers_the_largest_table_without_falling_back`
    /// pins the margin separately so it cannot erode silently.
    #[test]
    fn symbols_are_unique_within_a_universe() {
        let mut s = Stream::new(7, Domain::Naming, 0);
        let mut taken = Vec::new();
        for _ in 0..120 {
            let (sym, _, _) = mint(&mut s, &mut taken);
            assert_eq!(
                taken.iter().filter(|t| **t == sym).count(),
                1,
                "duplicate {sym}"
            );
        }
    }

    /// 120 rather than 500, for the reason above: past 165 every draw is the
    /// `Q{n}` fallback, whose name carries a digit and fails the lowercase
    /// assertion.
    #[test]
    fn symbols_are_well_formed() {
        let mut s = Stream::new(3, Domain::Naming, 0);
        let mut taken = Vec::new();
        for _ in 0..120 {
            let (sym, name, _) = mint(&mut s, &mut taken);
            assert!((1..=2).contains(&sym.len()), "bad symbol length: {sym}");
            assert!(sym.chars().next().is_some_and(|c| c.is_ascii_uppercase()));
            assert!(name.len() >= 4, "name too short: {name}");
            assert!(name.chars().all(|c| c.is_ascii_lowercase()));
        }
    }

    /// The blocklist is worthless if it is never exercised. Prove it rejects.
    #[test]
    fn blocklist_actually_rejects() {
        assert!(is_real("He", "helium"));
        assert!(is_real("Xx", "carbon"));
        assert!(is_real("Fe", "vorium"));
        assert!(!is_real("Vo", "vorium"));
    }
}
