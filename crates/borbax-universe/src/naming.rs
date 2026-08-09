//! Procedural element naming, with a real-name blocklist (spec §5, G2).
//!
//! The blocklist below exists to **exclude**, for every *perturbed*
//! universe: it is never read as data, never used to derive a property, and
//! carries no information beyond "this string is taken". Without it a
//! generated universe could mint an element called `He` or `carbon`, which
//! would falsely imply a mapping onto real chemistry for a universe that was
//! never meant to be that mapping.
//!
//! **G2 was revised 2026-08-07** (spec §5): the one seed-equivalent
//! configuration whose unperturbed constants reproduce the real periodic
//! table — generated too, just unperturbed, not exempt for any other reason
//! — is now expected to use the real names rather than be blocked from
//! them. That exception is implemented in this file (`REAL_TABLE`,
//! `IdentityWitness`, `real_name`) but has no consumer yet — wiring it to a
//! live `Universe` is issue #26's Task 26.1, and changes nothing about what
//! the blocklist below does today.

use crate::perturbation::Rung;
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

/// Symbol paired with name, for the identity configuration only (Decision 9,
/// issue #26). **A third copy, deliberately** — `REAL_ELEMENT_SYMBOLS` and
/// `REAL_WORDS` stay exactly as they are (flat exclusion lists, read by two
/// `xtask` guards keyed on their exact source text), and this table exists
/// beside them rather than being built from them, matching `xtask`'s own
/// `REAL_WORDS_COPY` precedent: "duplicated, never imported, and the
/// duplication is the point." `the_real_table_corresponds_to_the_blocklists`
/// is what actually enforces the correspondence a shared-storage refactor
/// would otherwise buy for free — and buy at the cost of a second reader on
/// `REAL_ELEMENT_SYMBOLS`/`REAL_WORDS`, which the identity guard's own
/// "exactly one reader" requirement (see [`real_name`]) forbids.
///
/// Ordered by atomic number, `Z = 1..=118` — index `Z - 1`. [`real_name`] is
/// this table's only reader.
const REAL_TABLE: &[(&str, &str)] = &[
    ("H", "hydrogen"),
    ("He", "helium"),
    ("Li", "lithium"),
    ("Be", "beryllium"),
    ("B", "boron"),
    ("C", "carbon"),
    ("N", "nitrogen"),
    ("O", "oxygen"),
    ("F", "fluorine"),
    ("Ne", "neon"),
    ("Na", "sodium"),
    ("Mg", "magnesium"),
    ("Al", "aluminium"),
    ("Si", "silicon"),
    ("P", "phosphorus"),
    ("S", "sulfur"),
    ("Cl", "chlorine"),
    ("Ar", "argon"),
    ("K", "potassium"),
    ("Ca", "calcium"),
    ("Sc", "scandium"),
    ("Ti", "titanium"),
    ("V", "vanadium"),
    ("Cr", "chromium"),
    ("Mn", "manganese"),
    ("Fe", "iron"),
    ("Co", "cobalt"),
    ("Ni", "nickel"),
    ("Cu", "copper"),
    ("Zn", "zinc"),
    ("Ga", "gallium"),
    ("Ge", "germanium"),
    ("As", "arsenic"),
    ("Se", "selenium"),
    ("Br", "bromine"),
    ("Kr", "krypton"),
    ("Rb", "rubidium"),
    ("Sr", "strontium"),
    ("Y", "yttrium"),
    ("Zr", "zirconium"),
    ("Nb", "niobium"),
    ("Mo", "molybdenum"),
    ("Tc", "technetium"),
    ("Ru", "ruthenium"),
    ("Rh", "rhodium"),
    ("Pd", "palladium"),
    ("Ag", "silver"),
    ("Cd", "cadmium"),
    ("In", "indium"),
    ("Sn", "tin"),
    ("Sb", "antimony"),
    ("Te", "tellurium"),
    ("I", "iodine"),
    ("Xe", "xenon"),
    ("Cs", "caesium"),
    ("Ba", "barium"),
    ("La", "lanthanum"),
    ("Ce", "cerium"),
    ("Pr", "praseodymium"),
    ("Nd", "neodymium"),
    ("Pm", "promethium"),
    ("Sm", "samarium"),
    ("Eu", "europium"),
    ("Gd", "gadolinium"),
    ("Tb", "terbium"),
    ("Dy", "dysprosium"),
    ("Ho", "holmium"),
    ("Er", "erbium"),
    ("Tm", "thulium"),
    ("Yb", "ytterbium"),
    ("Lu", "lutetium"),
    ("Hf", "hafnium"),
    ("Ta", "tantalum"),
    ("W", "tungsten"),
    ("Re", "rhenium"),
    ("Os", "osmium"),
    ("Ir", "iridium"),
    ("Pt", "platinum"),
    ("Au", "gold"),
    ("Hg", "mercury"),
    ("Tl", "thallium"),
    ("Pb", "lead"),
    ("Bi", "bismuth"),
    ("Po", "polonium"),
    ("At", "astatine"),
    ("Rn", "radon"),
    ("Fr", "francium"),
    ("Ra", "radium"),
    ("Ac", "actinium"),
    ("Th", "thorium"),
    ("Pa", "protactinium"),
    ("U", "uranium"),
    ("Np", "neptunium"),
    ("Pu", "plutonium"),
    ("Am", "americium"),
    ("Cm", "curium"),
    ("Bk", "berkelium"),
    ("Cf", "californium"),
    ("Es", "einsteinium"),
    ("Fm", "fermium"),
    ("Md", "mendelevium"),
    ("No", "nobelium"),
    ("Lr", "lawrencium"),
    ("Rf", "rutherfordium"),
    ("Db", "dubnium"),
    ("Sg", "seaborgium"),
    ("Bh", "bohrium"),
    ("Hs", "hassium"),
    ("Mt", "meitnerium"),
    ("Ds", "darmstadtium"),
    ("Rg", "roentgenium"),
    ("Cn", "copernicium"),
    ("Nh", "nihonium"),
    ("Fl", "flerovium"),
    ("Mc", "moscovium"),
    ("Lv", "livermorium"),
    ("Ts", "tennessine"),
    ("Og", "oganesson"),
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

/// Proof that a universe is the one seed-equivalent configuration whose
/// unperturbed constants reproduce the real periodic table (spec §5, G2's
/// 2026-08-07 revision).
///
/// **Confined to this crate by the type system, not merely by convention —
/// and an earlier version of this doc comment overstated that as
/// "unforgeable by the type system" full stop, which a review falsified in
/// under a minute from a different crate.** `Rung::IDENTITY` and
/// `Rung::new(0)` are both `pub`, always succeed, and carry no information
/// about any actual universe — so if `Self::new` were `pub`, any crate in
/// the workspace could write `IdentityWitness::new(Rung::IDENTITY)` and
/// mint a witness with zero evidence that any universe exists, let alone
/// one at the identity rung. Measured: from `borbax-molecule`, with no
/// `unsafe`, three lines dumped all 118 real element names, and a
/// `Rung::draw(seed)` off a *non-identity* rung changed nothing about the
/// answer, because the witness never consulted it. `new` is `pub(crate)`
/// for exactly this reason — the field being private stops nothing on its
/// own, since a public fallible constructor is a public capability
/// regardless of what the type's internals look like.
///
/// **The type itself is `pub(crate)` too, added in `/review-pr` after that
/// same fix round.** Sealing only `new` left three routes open, all
/// reproduced end to end from `borbax-molecule` with `cargo xtask` and
/// `clippy -D warnings` both green: a `macro_rules!`-forged `Self(())`
/// (macro token streams are invisible to the `syn::visit`-based confinement
/// checks below), a `format!("{REAL_TABLE:?}")` table dump, and a plain
/// `pub fn` inside this file returning either the witness or the names
/// `real_name` reads. Sealing `IdentityWitness` and [`real_name`] to
/// `pub(crate)` closes all three the same way `new` closed the first one —
/// the compiler refuses to compile a crate-external-`pub` item whose
/// signature names a `pub(crate)` type (`E0446`), and no external crate can
/// call a `pub(crate)` function regardless of what wraps the call. The
/// three `xtask` structural checks below stay as defence-in-depth against a
/// mistake made *inside* this crate (a second reader, a second
/// construction site); they are no longer the sole enforcement of
/// confinement from outside it.
///
/// **No `PhysicsVersion` field, and that is a known, deliberate gap in this
/// prerequisite, not an oversight.** A `Rung` alone cannot express "this
/// version has real-physics constants to reproduce" — under
/// `PhysicsVersion::V1`, which has none, `rung == 0` is not meaningfully
/// "the real periodic table" at all. Wiring a version-aware witness (and the
/// `Rung` field on `Universe` this needs to be checked against a live
/// universe) is Task 26.1's job, deliberately deferred — see
/// `crate::perturbation`'s own module doc for the full reasoning and the
/// determinism-auditor finding that produced it. Task 26.1 is also where
/// this crate gains its first real caller of `new` — until then the
/// `pub(crate)` capability is unused by design, not an oversight either.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct IdentityWitness(());

impl IdentityWitness {
    /// A witness, or `None` unless `rung` is the identity rung.
    ///
    /// **`pub(crate)`, not `pub` — see [`Self`]'s own doc for why a public
    /// fallible constructor over a public, always-succeeding `Rung::IDENTITY`
    /// is not actually a capability restriction at all.** `Universe::identity_witness`
    /// (Task 26.1, not yet built) is meant to be the single public door,
    /// checking a *live* universe's own rung rather than a caller-supplied
    /// one.
    ///
    /// **No in-crate caller exists yet — deliberately, until Task 26.1 —
    /// which `dead_code` would otherwise report.** The `#[expect]` is
    /// `cfg_attr`-gated to `not(test)` because this crate's own tests *are*
    /// callers, and an unconditional `#[expect(dead_code)]` fires
    /// `unfulfilled_lint_expectations` on the test target for exactly that
    /// reason — verified.
    ///
    /// **The one place this type may be constructed in the crate — `xtask`
    /// checks it, at the `syn` level.** Written as `Self(())`, not
    /// `IdentityWitness(())`: `clippy::use_self` (workspace `-D warnings`)
    /// forbids the latter, so the plan's own framing — "a textual count of
    /// the substring `IdentityWitness(` is at least two by construction" —
    /// does not hold for this codebase's actual constructor shape; a naive
    /// textual scan would find exactly one occurrence (this struct's
    /// declaration) and could pass while counting the wrong thing entirely.
    /// The `syn`-level guard resolves this correctly: it walks
    /// `impl IdentityWitness` specifically and treats a `Self(..)` call
    /// inside it as a construction of this type, which is the one place
    /// AST-level structure can do what a textual scan cannot.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "no in-crate caller until Task 26.1 gives Universe a rung and a public \
                      identity_witness() door. Held pub(crate) rather than pub so a \
                      downstream crate cannot mint the capability from Rung::IDENTITY; this \
                      attribute pays for the deferral visibly instead of paying for it with \
                      public API surface"
        )
    )]
    pub(crate) const fn new(rung: Rung) -> Option<Self> {
        if rung.is_identity() {
            Some(Self(()))
        } else {
            None
        }
    }
}

/// The real symbol and name for atomic number `z` (`1..=118`), or `None`
/// outside that range — reachable only with proof of the identity
/// configuration in hand.
///
/// **No `Rung` parameter.** An earlier revision of this prerequisite gave
/// this accessor a `rung: Rung` parameter *and* an `IdentityWitness`
/// parameter, which reopened exactly the hazard the witness exists to
/// close: a caller could construct any `Rung` value directly and pass it
/// here without ever going through `IdentityWitness::new`'s fallible check.
/// The witness alone is the capability; `_` is deliberate — this accessor
/// reads nothing from it but the proof that constructing it succeeded.
///
/// **By value, not by reference** — `IdentityWitness` is a zero-sized
/// `Copy` type, so a `&IdentityWitness` costs a pointer to point at nothing,
/// strictly worse than moving the value itself. `clippy::trivially_copy_pass_by_ref`
/// only checks crate-internal signatures, which is why this did not fire
/// before this accessor became `pub(crate)`.
///
/// **`pub(crate)`, not `pub` — added in `/review-pr`, same round and same
/// reason as [`IdentityWitness`] itself.** See that type's doc for the
/// three bypass routes this closes. No in-crate caller exists yet, same as
/// `IdentityWitness::new` — deliberately, until Task 26.1 — so the
/// `#[expect]` is `cfg_attr`-gated to `not(test)` for the same reason: this
/// crate's own tests are callers, and an unconditional `#[expect(dead_code)]`
/// would fire `unfulfilled_lint_expectations` on the test target.
#[must_use]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "no in-crate caller until Task 26.1 gives Universe a rung and a public \
                  identity_witness() door, matching IdentityWitness::new's own reason"
    )
)]
pub(crate) fn real_name(_: IdentityWitness, z: usize) -> Option<(&'static str, &'static str)> {
    z.checked_sub(1).and_then(|i| REAL_TABLE.get(i)).copied()
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

    /// The independent oracle every test below uses. **Never read `REAL_WORDS`
    /// to build this** — it is the blocklist under test, and its whole defect
    /// was a silent gap in exactly this list.
    ///
    /// Verified live on `main` at `cb3a28c`, before this fix: `REAL_WORDS`
    /// covered 75 of these 118 (not 78 — that figure counted the three
    /// alternate spellings `aluminum`/`sulphur`/`cesium` as additional
    /// coverage, which they are not) — it jumped from `cerium` (Z=58) straight
    /// to `tungsten` (Z=74), dropping the whole Z=59..=73 lanthanide run, then
    /// dropped Z=75..=77, Z=91 and everything from Z=95 onward. 43 names
    /// missing in total.
    ///
    /// **Only one of the 43, `thulium`, was ever mint-time reachable.**
    /// `mint`'s grammar (`ONSETS x NUCLEI x CODAS x SUFFIXES`) can produce
    /// exactly six of these 118 names — `barium`, `boron`, `lithium`, `neon`,
    /// `thorium`, `thulium` — and five were already blocked. The other 42
    /// cannot be generated by the current grammar and could not have caused
    /// this specific defect; they stay in the blocklist because `REAL_WORDS`
    /// has a second consumer with no grammar constraint at all — `xtask`'s
    /// literal-content scan, which would flag any of the 118 appearing
    /// hardcoded in a label or comment regardless of whether the generator
    /// could ever mint it.
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
    ///
    /// **Length alone does not prove completeness, and a first draft of this
    /// test claimed it did.** A hand-typed 118-entry array's realistic failure
    /// is a duplicated line paired with a dropped one — the length stays 118,
    /// `REAL_WORDS.contains` still passes for the duplicate, and the oracle
    /// acquires exactly the defect it exists to catch, silently. The
    /// distinctness check below is what actually makes "checked structurally"
    /// true rather than assumed.
    #[test]
    fn the_blocklist_names_every_real_element() {
        assert_eq!(
            REAL_ELEMENT_NAMES.len(),
            REAL_ELEMENT_SYMBOLS.len(),
            "the reference list above must cover the same 118 elements \
             REAL_ELEMENT_SYMBOLS does, or the zip below is comparing the wrong pairs"
        );
        let distinct: std::collections::BTreeSet<&str> = REAL_ELEMENT_NAMES.into_iter().collect();
        assert_eq!(
            distinct.len(),
            REAL_ELEMENT_NAMES.len(),
            "REAL_ELEMENT_NAMES contains a duplicate, which means it is silently \
             missing a real element name — the length check alone cannot see this"
        );
        for (symbol, name) in REAL_ELEMENT_SYMBOLS.iter().zip(REAL_ELEMENT_NAMES) {
            assert!(
                REAL_WORDS.contains(&name),
                "REAL_WORDS is missing '{name}' (element {symbol}) — a generated \
                 element could be minted with that name and pass is_real() falsely"
            );
        }
    }

    /// **`REAL_TABLE` is a deliberate third copy (see its own doc comment) —
    /// this test is what keeps that duplication honest rather than merely
    /// declared.** Checks correspondence against both blocklists, never
    /// shared storage with either: same length, same symbols in the same
    /// order as `REAL_ELEMENT_SYMBOLS`, same names in the same order as
    /// `REAL_ELEMENT_NAMES`, and every name present in `REAL_WORDS`.
    ///
    /// **The `REAL_ELEMENT_NAMES` zip is not redundant with the
    /// `REAL_WORDS.contains` check below it, found in review.** `REAL_WORDS`
    /// is a flat, unordered 176-entry set (118 element names plus 58
    /// unrelated chemistry/biology terms) — membership in it cannot catch
    /// two adjacent names being *transposed* (`B` paired with "carbon",
    /// `C` paired with "boron"): both names are still somewhere in
    /// `REAL_WORDS`, so `contains` passes for both while the pairing itself
    /// is wrong. `REAL_ELEMENT_NAMES` is ordered by atomic number, the same
    /// axis `REAL_TABLE` is, so zipping against it is the one check that
    /// actually verifies the *pairing*, not just the *wordlist*.
    #[test]
    fn the_real_table_corresponds_to_the_blocklists() {
        assert_eq!(
            REAL_TABLE.len(),
            REAL_ELEMENT_SYMBOLS.len(),
            "REAL_TABLE must cover the same 118 elements REAL_ELEMENT_SYMBOLS does"
        );
        for ((&(table_symbol, table_name), &blocklist_symbol), reference_name) in REAL_TABLE
            .iter()
            .zip(REAL_ELEMENT_SYMBOLS)
            .zip(REAL_ELEMENT_NAMES)
        {
            assert_eq!(
                table_symbol, blocklist_symbol,
                "REAL_TABLE's symbol order has drifted from REAL_ELEMENT_SYMBOLS's"
            );
            assert_eq!(
                table_name, reference_name,
                "REAL_TABLE pairs {table_symbol} with '{table_name}', but the \
                 independently-ordered REAL_ELEMENT_NAMES says this position is \
                 '{reference_name}' — the pairing has drifted, not just the wordlist"
            );
            assert!(
                REAL_WORDS.contains(&table_name),
                "REAL_TABLE names '{table_name}' ({table_symbol}), which REAL_WORDS \
                 does not contain — the two blocklists and this table have drifted apart"
            );
        }
    }

    /// Decision 9's exhaustive test, restated on the constructor: every
    /// non-identity rung must fail to produce a witness. `Rung::all_off_identity()`
    /// is the corpus, derived from `Rung::MAX` rather than a hardcoded
    /// `1..=RUNGS` literal that could drift from it.
    #[test]
    fn identity_witness_is_none_for_every_off_identity_rung() {
        for rung in Rung::all_off_identity() {
            assert!(
                IdentityWitness::new(rung).is_none(),
                "rung {rung:?} produced a witness — only Rung::IDENTITY should"
            );
        }
    }

    #[test]
    fn identity_witness_exists_only_for_the_identity_rung() {
        assert!(IdentityWitness::new(Rung::IDENTITY).is_some());
    }

    #[test]
    fn real_name_answers_every_atomic_number_and_nothing_outside_it() {
        let witness = IdentityWitness::new(Rung::IDENTITY)
            .unwrap_or_else(|| unreachable!("Rung::IDENTITY always constructs a witness"));
        assert_eq!(
            real_name(witness, 0),
            None,
            "atomic number 0 does not exist"
        );
        assert_eq!(
            real_name(witness, 1),
            Some(("H", "hydrogen")),
            "atomic number 1 is hydrogen"
        );
        assert_eq!(
            real_name(witness, 118),
            Some(("Og", "oganesson")),
            "atomic number 118 is oganesson, the last entry in REAL_TABLE"
        );
        assert_eq!(
            real_name(witness, 119),
            None,
            "atomic number 119 is past REAL_TABLE's end"
        );
        assert_eq!(
            real_name(witness, usize::MAX),
            None,
            "an absurd atomic number must not panic or wrap into range"
        );
    }

    /// **Exhaustive, not probabilistic, and this is what makes it strictly
    /// stronger than the seed sweep below at a fraction of the cost.** The
    /// generator's name space is finite and known in advance — `ONSETS x
    /// NUCLEI x CODAS x SUFFIXES` — so every name `mint` could ever produce
    /// can be checked directly against `is_real`, with no seed, no RNG, and
    /// no probability of missing a case. See the note on `REAL_ELEMENT_NAMES`
    /// above for the six grammar-producible real names this finds.
    ///
    /// The dummy symbol is deliberately never real — `is_real` checks symbol
    /// and name independently, and the symbol side is already exhaustively
    /// safe (`REAL_ELEMENT_SYMBOLS` was complete before this fix and is
    /// unchanged by it), so this isolates the name check the fix changed.
    #[test]
    fn every_name_the_grammar_can_produce_is_checked_against_the_blocklist() {
        let mut examined: usize = 0;
        for onset in ONSETS {
            for nucleus in NUCLEI {
                for coda in CODAS {
                    for suffix in SUFFIXES {
                        let name = format!("{onset}{nucleus}{coda}{suffix}");
                        examined += 1;
                        if REAL_ELEMENT_NAMES.contains(&name.as_str()) {
                            assert!(
                                is_real("Zz", &name),
                                "the grammar can produce the real element name \
                                 '{name}', and is_real() does not catch it"
                            );
                        }
                    }
                }
            }
        }
        let expected = ONSETS.len() * NUCLEI.len() * CODAS.len() * SUFFIXES.len();
        assert_eq!(
            examined, expected,
            "examined count must equal the grammar's full production count"
        );
    }

    /// The integration companion to the exhaustive test above: exercise the
    /// real `Stream` / `taken` / retry pipeline end-to-end, not just the
    /// static grammar. 2 000 seeds, not 20 000 — measured to carry identical
    /// discriminating power to the defect this test was added for at roughly
    /// a tenth of the cost, since the grammar (checked exhaustively above) is
    /// what actually bounds which collisions exist, not the seed count.
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
    /// `name` needs no `.to_ascii_lowercase()` here: every `NUCLEI`/`ONSETS`/
    /// `CODAS`/`SUFFIXES` entry is already lowercase, so `mint`'s generated
    /// path can never produce anything else.
    ///
    /// **Neither pinned universe digest would have caught the defect this
    /// closes, and not for the same reason.** `the_universe_digest_is_pinned`
    /// (`element.rs`) does mix `symbol`/`name` into its hash, but sweeps seeds
    /// `0..64` only — the collision this fix is for first reproduces at seed
    /// 1631. `the_assembled_universe_digest_is_pinned` never mixes `symbol`/
    /// `name` at all, at any seed, so it could not have seen this regardless
    /// of range.
    #[test]
    fn no_seed_in_a_sweep_mints_a_real_element_name() {
        let mut generated: u64 = 0;
        for seed in 0..2_000u64 {
            let mut s = Stream::new(seed, Domain::Naming, 0);
            let mut taken = Vec::new();
            for _ in 0..120 {
                let (_, name, provenance) = mint(&mut s, &mut taken);
                if provenance == Provenance::Generated {
                    generated += 1;
                }
                assert!(
                    !REAL_ELEMENT_NAMES.contains(&name.as_str()),
                    "seed {seed} minted a real element name: {name}"
                );
            }
        }
        // Counting only Provenance::Generated, not every mint, is what makes
        // this sensitive to the fallback path going quiet: `Q{n}`/`quorium{n}`
        // never reaches `is_real` at all (see `mint`'s fallback arm), so a
        // sweep that counted every attempt regardless of provenance would
        // report full coverage while silently exercising a path this test
        // cannot check. Not live today — 0 fallbacks measured over 2 000
        // seeds x 120 mints — but the counter should be able to say so.
        assert_eq!(
            generated,
            2_000 * 120,
            "expected every draw in this sweep to reach the generator, not the \
             fallback — a nonzero gap means Provenance::Fallback fired and this \
             test's coverage claim no longer holds"
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
