//! The generated periodic table (spec §7.1).
//!
//! **Element `N` is `N` copies of one base unit, packed.** `mass`, `period`,
//! `group`, `valence`, `radius`, `affinity` and `energy_per_unit` are functions
//! of that packing rather than curves chosen because they looked plausible —
//! principle 2, applied to the file that defines what an atom is.
//!
//! One property is **not** a function of the packing, and an earlier version
//! of this header claimed "every property below" was, which the file's own
//! later text contradicted: `abundance` is a function of the fusion *process*
//! (one drawn constant and the unit count, with no packing quantity in it).
//! A second exception, `outer_fill_band`, existed for the same reason and was
//! removed entirely in Task 26.1 Step 8 — a rebinning of `(period, group)`
//! that never carried chemistry (see the git history for its own header,
//! which documented three successive incorrect claims about what it meant).
//!
//! See Task 4's preamble in the plan for what is claimed and, more importantly,
//! what is not: there is no shape-diversity advantage, only the measured
//! finding that a fusion-derived radius series costs nothing.
//!
//! Periodicity is a *consequence* here. A shell closes when it fills, and the
//! elements at and around a closure form families with a cause: a closed shell
//! has no frontier and so no bonding slots, and the element one unit past any
//! closure has exactly one.

use crate::PhysicsVersion;
use crate::naming;
use crate::orbital::{self, ElectronicProperties, GapConsts, OrbitalConsts};
use crate::packing::{self, PackingConsts};
use crate::perturbation::{MigratedConstant, Rung, draw_symmetric};
use borbax_rng::{Domain, Stream};
use borbax_units::{Mass, Quanta, Span, canonical_cmp, det_math};

/// Index into a universe's element table.
///
/// `u8` because tables are capped well below 256, which keeps `Mol12`'s atom
/// array one byte per slot — **by storing `[ElementId; MAX_ATOMS]`, not
/// `[u8; MAX_ATOMS]`.**
///
/// That distinction is load-bearing and a review found Task 6's draft on the
/// wrong side of it. `Self::index` is `pub(crate)`, so a downstream crate
/// holding raw `u8`s has a route *in* ([`Self::from_index`]) and none back out;
/// `self.elem[idx] = e.0` and `ElementId(self.elem[i])` are both hard compile
/// errors outside this crate. Storing the id itself costs the same byte
/// (asserted below), keeps `Copy`, and means the array cannot hold a value that
/// was never an id.
///
/// **The field is `pub(crate)`, and that is the whole point of the type.** Until
/// Task 5 an `ElementId` was write-only — minted by [`generate_elements`] and
/// read by nobody — and the shape it was heading for was `elements[id.0 as
/// usize]` at every call site: a panic site, in a workspace where
/// `clippy::indexing_slicing` is `warn` under `-D warnings`, repeated once per
/// consumer. Closing the field means the only route from an id to an element is
/// [`PeriodicTable::get`], which is fallible, so the panic has one home and that
/// home returns `None`.
///
/// [`Self::from_index`] is deliberately *not* the inverse of that rule. It goes
/// `usize -> ElementId`, which cannot index anything on its own; validity
/// against a particular table stays the table's question, because an id is a
/// bounded integer and a table is what makes one meaningful.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ElementId(pub(crate) u8);

// The doc above promises `Mol12` one byte per slot. Enforce it rather than
// asserting it: a newtype over `u8` is 1 byte with no `repr` attribute today,
// but adding a field or a niche would silently double `Mol12`'s atom array and
// the promise would rot the way a comment does.
const _: () = assert!(
    size_of::<ElementId>() == size_of::<u8>(),
    "ElementId must stay one byte — Mol12's atom array is sized on it"
);

impl ElementId {
    /// The first slot of any non-empty table.
    pub const ZERO: Self = Self(0);

    /// An id for the `i`th slot, or `None` if `i` cannot be an element index.
    ///
    /// Fails only on the `u8` bound — it says nothing about whether any
    /// particular table has such a slot. Ask [`PeriodicTable::get`] for that.
    #[must_use]
    pub fn from_index(i: usize) -> Option<Self> {
        u8::try_from(i).ok().map(Self)
    }

    /// The slot number, for indexing inside this crate only.
    pub(crate) fn index(self) -> usize {
        usize::from(self.0)
    }
}

/// One generated element: a cluster of [`Element::units`] base units, and the
/// properties that follow from how they pack.
///
/// **`#[non_exhaustive]`, so every field stays readable and none of it is
/// constructible outside this crate.** The coupling here is tighter than
/// `PackingConsts`'s — given `k`, all of `period`, `group`, `valence`, `radius`,
/// `affinity` and `mass` are functions of `units` alone, so
/// `Element { units: 3, valence: 4, .. }` is geometrically impossible and would
/// compile. Thirteen accessors on a record read field-by-field is boilerplate
/// that loses an argument later; this is the stdlib answer to "read freely,
/// construct never".
///
/// Timing is the whole argument: there is no downstream consumer today, so it
/// costs nothing. After Task 5 writes fixtures against it, it is a breaking
/// change — and a fixture the generator can never emit lets Task 5's tests pass
/// on states that do not exist.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Element {
    // **There is deliberately no `id` field.** It was here, and it encoded the
    // same fact as the element's position in [`PeriodicTable`] — which is the
    // duplicate-encoding defect `bonds.rs`'s own header rejects one level down,
    // reproduced inside the table. `PeriodicTable::get` indexes by position,
    // `BondEnergyMatrix` builds its cells by position and queries them by id,
    // and `iter` used to read the field: three consumers depending on "slot `i`
    // holds id `i`", with nothing establishing it. A table built by hand with
    // two elements transposed made `iter` and `get` disagree, and priced the
    // wrong pair — in range, plausible and silent, which `get`'s own doc calls
    // the worse of the two failures this type exists to prevent. The id now has
    // no *duplicate* source: position is the only thing that determines it, so
    // nothing can disagree with anything because there is nothing left to
    // disagree with. (Not "exactly one source" — [`ElementId::from_index`] is
    // public and `Universe::generate` mints one from a draw. Neither can
    // disagree with position, since both take a position as input, but the
    // stronger sentence stood here for two rounds after being retracted 190
    // lines below. When a claim is retracted, grep the file for the claim.)
    /// Generated symbol, 1–2 characters (§5, G2).
    pub symbol: String,
    /// Generated name (§5, G2).
    pub name: String,
    /// Base units in this element's cluster. The element *is* this number.
    pub units: usize,
    /// Which shell is filling. The period index, derived rather than drawn.
    ///
    /// **V1**: the shell index (0-based). **V2**: the outermost occupied
    /// principal quantum number `outer_n` (1-based) — the two versions do
    /// not share an indexing convention, since each is self-consistent
    /// within its own table and nothing compares a period number across
    /// versions.
    pub period: u8,
    /// Units in the incomplete outer shell. The group index.
    ///
    /// **V1**: fill-fraction position, `outer` (a function of `(period,
    /// group)`, not of `group` alone — see the historical note this field's
    /// predecessor `outer_fill_band` carried, now removed). **V2**: the
    /// frontier subshell's own angular momentum `l` directly (`0` = s, `1` =
    /// p, `2` = d, `3` = f, ...) — a structural readout, not a fitted
    /// quantity. See [`Block`].
    pub group: u8,
    /// Cluster mass: units carried, less the mass defect of made contacts.
    pub mass: Mass,
    /// Bonding slots: unmade lateral contacts on the frontier, over the
    /// lateral coordination. Zero at a closure because a closed shell has no
    /// frontier — not because a branch says so.
    pub valence: u8,
    /// Exposed fraction of the cluster, mapped to [-1, +1]. Drives
    /// complementarity in §8.3. **Not** electronegativity: different range,
    /// different meaning, different behaviour (G3).
    pub affinity: f64,
    /// Enclosing radius of the packed cluster. Rises monotonically with the
    /// occupied shell and never jumps — the opposite of a real atomic radius,
    /// which falls across a period (G3).
    pub radius: Span,
    /// Binding energy per unit: made contacts, less radial strain.
    ///
    /// [`crate::bonds::BondEnergyMatrix`] is derived from this rather than drawn
    /// independently. Two encodings of "these elements bind well" in
    /// incommensurable units is the defect that made the predecessor's `peak`
    /// unfalsifiable, and an independently drawn bond matrix reproduces it one
    /// level up.
    ///
    /// **G4, decided in Task 5 rather than inherited: this is `Quanta`.** Task 4
    /// left it a bare `f64` and flagged it, on the argument that the value is
    /// compared, summed and scaled inside this file where the newtype buys
    /// nothing. That argument is true and it stops at the crate boundary — the
    /// field is `pub` on a `pub` struct, so a bare `f64` here *is* an energy
    /// crossing the public API undimensioned, which is the "carried through
    /// three functions and quietly defeated" shape §5's G4 exists to stop.
    /// Precedence 1 has nothing to trade against, so the tie went to the type.
    ///
    /// The measured cost was lower than the flag implied, which is worth
    /// recording because the flag is what a future reader will find: the
    /// arithmetic below needed no `.get()` calls at all. `deficit` now reads
    /// `(Quanta - Quanta) / Quanta`, and `borbax-units` types that quotient as
    /// `f64` — so the dimensionless-ness of a decay rate became a fact the
    /// compiler checks instead of a convention. The universe digest did not
    /// move, which is the evidence that this was a typing change and not a
    /// physics one.
    pub energy_per_unit: Quanta,
    /// How far this element sits below the binding-energy peak, on `0..=1`.
    /// 0 at the peak, rising toward both extremes.
    ///
    /// **Renamed from `decay_rate`, because "decay" already meant something
    /// else in this project and the collision was actively misleading.** Borbax
    /// has two unrelated processes that the word covers:
    ///
    /// 1. **Bond cleavage** — heat and solvent attack breaking *chemical* bonds,
    ///    `k = A·exp(−E_bond/T)` (§9.4, §9.5). Per-bond, driven by
    ///    [`crate::bonds::BondEnergyMatrix`]. This is "the decay band", the
    ///    quantity CLAUDE.md calls the most sensitive parameter in the system.
    /// 2. **This field** — a per-*element* instability derived from
    ///    [`Self::energy_per_unit`], which is the packing/cluster model. Nothing
    ///    to do with bonds; it is the nuclear-side analogue.
    ///
    /// They meet in exactly one place: an unstable element sitting inside a
    /// molecule damages it, which is §9.5's **radiogenic** cleave channel. So
    /// this is a *source term* for one of the channels in (1), not an instance
    /// of it. Sharing a name invited exactly the conflation it produced.
    ///
    /// **And the new name deliberately does not say "rate".** §7.1 calls the old
    /// one "probability per world-year" while Task 15 consumes it as a Gillespie
    /// propensity, which is unbounded above — a recorded, unresolved
    /// contradiction, and the reason the `> 1.0` clamp below is a category error
    /// rather than a saturation. A dimensionless distance is what the value
    /// actually *is*; Task 15 decides how it becomes a rate, and now has to say
    /// so explicitly rather than inheriting an answer from a field name.
    ///
    /// (An earlier draft called this `stability`, which inverted the sense
    /// against its own description. `instability` keeps the sense — larger means
    /// less stable — without claiming a unit.)
    ///
    /// **One-axis limitation, stated because it bounds what this can model.**
    /// Real half-lives are driven by position in a *two*-dimensional (protons,
    /// neutrons) space — beta decay corrects a wrong ratio, alpha decay and
    /// spontaneous fission shed mass from nuclei too large for the strong force
    /// to hold against electrostatic repulsion. Borbax has one axis (`units`),
    /// so this can express "far from the peak is unstable" but has no second
    /// axis for a *ratio* to be wrong along. The shape it gives is
    /// alpha/fission-like, not beta-like. That is a property of invented
    /// physics, not a defect — but do not reach for a neutron count to fix it.
    pub instability: f64,
    /// Relative abundance. **The gel lever** (§7.2) — see [`generate_elements`].
    pub abundance: f64,
}

/// Which subshell block an element's frontier sits in (Task 26.1, issue #26's
/// Step 8).
///
/// **`V2` only**; V1 has no orbital picture for this to describe.
/// Derived from [`Element::group`], which under V2 stores the frontier
/// subshell's own `l` directly (a structural readout of `orbital.rs`'s own
/// fold state, not a fitted or calibrated quantity — see the free function
/// `block` in this module).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Block {
    /// `l == 0`.
    S,
    /// `l == 1`.
    P,
    /// `l == 2`.
    D,
    /// `l == 3`.
    F,
    /// `l >= 4` — reachable in principle (`orbital.rs`'s `MAX_N = 9` allows
    /// up to `l = 8`) but not observed within any table's drawn `n_elements`
    /// range (`60..=120`) at the identity configuration; kept rather than
    /// treated as unreachable because a perturbed universe's screening
    /// coefficients are not proven to keep it that way.
    Other(u8),
}

/// Which block an angular momentum `l` sits in — **a function on `l`
/// directly, not on `Element`** (Step 8's own text: "`block()` is a method
/// on `l`"), so a V1 caller cannot silently misread `Element::group` — which
/// means something else entirely for V1 — as though it were a V2 `l` value.
/// A caller must already know it holds a V2 element's `group`.
///
/// **`pub(crate)`, narrower than P4's eventual call-site guard needs to
/// restrict — safe to widen once P4 lands, never to narrow.** P4 itself
/// (the block-field read guard) lands right after this exists, not before,
/// per the plan's own sequencing; this function existing is what P4 has
/// something to guard.
///
/// **No in-crate caller yet — deliberately, until P4 gives it one.** The
/// `#[expect]` is `cfg_attr`-gated to `not(test)` for the same reason every
/// other such attribute in this crate is: this module's own tests are
/// callers, and an unconditional `#[expect(dead_code)]` fires
/// `unfulfilled_lint_expectations` on the test target.
#[must_use]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "no in-crate caller until P4's block-field read guard lands with its own \
                  allowlisted call sites — see this function's own doc"
    )
)]
pub(crate) const fn block(l: u8) -> Block {
    match l {
        0 => Block::S,
        1 => Block::P,
        2 => Block::D,
        3 => Block::F,
        other => Block::Other(other),
    }
}

/// Which shell/orbital law a table follows (Task 26.1, issue #26's Step 8).
///
/// The one member of [`ShellPattern`] that is genuinely version-specific.
/// `closures`, `fallback_symbols`, `peak` and `eps` describe both laws
/// equally, since each version populates its own values for those from its
/// own mechanism; only `k` has no V2 equivalent at all — V2's law is the
/// whole orbital fold (`orbital.rs`), not a single constant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ShellLaw {
    /// V1's packing shell law: shell `n` holds `k*n^2 + 2` units.
    V1 {
        /// The drawn shell-law constant.
        k: usize,
    },
    /// V2's screened-hydrogenic orbital fill — see `orbital.rs` (`pub(crate)`,
    /// so not linkable from this public doc comment).
    V2,
}

/// The derived shape of one universe's table.
///
/// **No `Eq`**, since the in-crate `eps` field wraps a float. That is the honest consequence
/// of carrying the energy scale here, and `PartialEq` is what the equality tests
/// actually use.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ShellPattern {
    /// Which shell/orbital law this table follows, and that law's own
    /// version-specific constant (`k` for V1; none for V2). See
    /// [`ShellLaw`]'s own doc for why only this one thing splits.
    pub law: ShellLaw,
    /// Cumulative unit counts at which a shell closes. Derived from `k` (V1)
    /// or from period boundaries (V2), never drawn directly.
    pub closures: Vec<usize>,
    /// How many symbols came from `naming::mint`'s exhaustion fallback rather
    /// than its grammar.
    ///
    /// **Zero for every table this scheme draws, and the point is that it is
    /// checkable at runtime.** `mint` computes `Provenance` and the only
    /// production call site used to drop it into `_`. A `Fallback` symbol
    /// violates [`Element::symbol`]'s own 1–2 character contract, so the day
    /// someone raises the element cap past the 165-symbol space the contract
    /// breaks *silently* — the signal existed and was discarded at the one place
    /// it could be seen. Keeping the fallback crude is right; discarding the
    /// fact that it fired is not.
    pub fallback_symbols: usize,
    /// Where the per-unit binding energy peaks. **Read off the finished
    /// series, never passed in** — that is the property distinguishing this
    /// from the predecessor, in which `peak` was an argument whose value was
    /// then reported back as an emergent minimum.
    pub peak: usize,
    /// Per-contact binding energy — **the table's own energy scale**.
    ///
    /// Carried here so [`crate::bonds`] can build the bond scale *from* it
    /// rather than drawing an independent one. A review measured that the bond
    /// matrix was exactly invariant to `energy_per_unit`'s magnitude (scale
    /// every energy by 1000, worst bond change 3.7e-16) because the capacity
    /// map is affine-invariant — so "how strongly things bind" was encoded
    /// twice, here in Quanta and there in a chosen interval, with neither able
    /// to check the other. That is the defect `bonds.rs` exists to reject,
    /// reproduced one axis over.
    ///
    /// **`pub(crate)`, not `pub`.** `#[non_exhaustive]` blocks external
    /// construction, not field *reads*, and `PeriodicTable::pattern` is public —
    /// so a `pub` field here would let `borbax-molecule` and everything above it
    /// read the raw drawn per-contact energy and build a second energy formula
    /// from it. That is the double-encoding this commit removed, made reachable
    /// one crate up. A companion `sigma` was added and then removed the same
    /// round: it had exactly one reader in the workspace, an inert line in a
    /// test's setup, which is not a reason to widen a public type.
    pub(crate) eps: Quanta,
}

/// One universe's periodic table: the elements, and the shell law they follow.
///
/// **This exists to make [`ElementId`] usable without an indexing panic**, and
/// it names the pair [`generate_elements`] used to return as a bare tuple. The
/// fields are private so that [`Self::get`] is the only route from an id to an
/// element — a `pub elements: Vec<Element>` would hand every consumer back the
/// `elements[id.0 as usize]` this type exists to prevent.
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodicTable {
    pattern: ShellPattern,
    elements: Vec<Element>,
    /// The laws this table was generated under. Private, like `pattern` and
    /// `elements` — [`Self::physics`] is the only route out, for the same
    /// reason `elements` isn't `pub`: a bare field lets a caller build a
    /// `BondEnergyMatrix` whose dispatch disagrees with the table it names,
    /// which is exactly the "in range, plausible and silent" failure this
    /// type already exists to close for element ids.
    physics: PhysicsVersion,
}

impl PeriodicTable {
    /// Assemble a table. In-crate only: outside, a table comes from
    /// [`generate_elements`] and nowhere else, so no caller can build one whose
    /// elements disagree with its shell law.
    pub(crate) fn new(
        pattern: ShellPattern,
        elements: Vec<Element>,
        physics: PhysicsVersion,
    ) -> Self {
        // **The bound belongs here, not at the read site.** An `ElementId` is a
        // `u8`, so a table longer than 256 cannot round-trip position through an
        // id: measured, a 300-slot table hands out 256 distinct ids, 44 slots
        // disagree between `iter` and `get`, and `BondEnergyMatrix::energy`
        // prices slot 299 as slot 255 — in range, plausible and silent, which is
        // the exact failure removing `Element::id` was meant to close.
        //
        // `assert!`, not `debug_assert!`: the workspace has already shipped a
        // guard present in `cargo test` and absent from the profile that mints
        // goldens — twice, in `borbax-units` (Task 2), not in this crate, which
        // an earlier version of this line misattributed. It costs one compare
        // against a 146 us generation.
        //
        // What it guards is *not* the production path: `table_size_is_workable`
        // already pins 60..=120 over 30 seeds. It guards `new`'s other callers —
        // two tests today, and any future production one — so it is not
        // redundant with that test and should not be deleted as though it were.
        //
        // Note `clippy::panic` does not see `assert!` (measured), so the green
        // gate is not clearance for this; it is a deliberate ruling.
        assert!(
            elements.len() <= Self::MAX_ELEMENTS,
            "a PeriodicTable of {} elements cannot address its own slots: an \
             ElementId is a u8, so position stops round-tripping past {}",
            elements.len(),
            Self::MAX_ELEMENTS
        );
        Self {
            pattern,
            elements,
            physics,
        }
    }

    /// The largest table an [`ElementId`] can address.
    /// The number of distinct `u8` values — which is *why* the cap is what it
    /// is, stated better than `1 + u8::MAX as usize` stated it. That form needed
    /// an `#[expect(clippy::as_conversions)]` whose reason claimed the cast
    /// "has to" stay because `usize::from` is not const (true, E0658) — but the
    /// conclusion was false, and a suppression whose reason begins "has to"
    /// should be tested against `rustc` before it is written.
    pub(crate) const MAX_ELEMENTS: usize = 1 << u8::BITS;

    /// The shell law this table was built under.
    #[must_use]
    pub const fn pattern(&self) -> &ShellPattern {
        &self.pattern
    }

    /// The laws this table was generated under.
    ///
    /// Read by `BondEnergyMatrix::generate` (`pub(crate)`, so not linkable
    /// from this public doc comment), the consumer P3 (issue #26's
    /// prerequisites) added this field for: it takes a table, not a
    /// generation call, so without this accessor it has nothing to dispatch
    /// a version-aware bond-generation path on.
    #[must_use]
    pub const fn physics(&self) -> PhysicsVersion {
        self.physics
    }

    /// Take the table apart. In-crate only, and used only by tests.
    ///
    /// An earlier version said "tests that predate this type" — round 3 then
    /// added `a_table_too_long_to_address_is_refused_at_construction`, which
    /// calls this and postdates it by 200 lines, falsifying a doc it never
    /// touched.
    ///
    /// Keeping them destructuring the same triple they always did (a pair,
    /// before P3 added `physics`) is what lets `the_universe_digest_is_pinned`
    /// keep its **destructuring** unchanged across this refactor, so a moved
    /// digest means moved physics-the-generated-quantity rather than a
    /// rewritten test. An earlier version of this sentence claimed the test
    /// stayed *textually* unchanged, which the same commit falsified: one line
    /// moved from `energy_per_unit.to_bits()` to `energy_per_unit.0.to_bits()`
    /// when the field became `Quanta`.
    #[cfg(test)]
    pub(crate) fn into_parts(self) -> (ShellPattern, Vec<Element>, PhysicsVersion) {
        (self.pattern, self.elements, self.physics)
    }

    /// The element with this id, or `None` if this table has no such slot.
    ///
    /// Fallible on purpose. Ids outlive the tables that mint them — one
    /// travelling in from another universe, or from a persisted molecule, is a
    /// number with no guarantee attached, and the two alternatives are both
    /// worse: indexing panics, and clamping silently answers with a *different
    /// element's* properties, which nothing downstream can detect.
    #[must_use]
    pub fn get(&self, id: ElementId) -> Option<&Element> {
        self.elements.get(id.index())
    }

    /// Every element with its id, in table order.
    ///
    /// **Ids come from position here**, which is what makes them agree with
    /// [`Self::get`] by construction rather than by convention.
    ///
    /// Not "the single source of an `ElementId`" — an earlier version of this
    /// line said that and it was false twice over: [`ElementId::from_index`] is
    /// public, and `Universe::generate` mints one from an RNG draw. Neither can
    /// *disagree* with position, since both take a position as input, so the
    /// guarantee holds — but the sentence a later reader would have relied on
    /// did not.
    ///
    /// The `unwrap_or` is not `saturating` and does not "keep it from aliasing
    /// slot 0"; it would alias slot **255**, which is a real slot in any table
    /// where it could fire. `Self::new` now caps the length so it cannot.
    pub fn iter(&self) -> impl Iterator<Item = (ElementId, &Element)> {
        // `map_while`, not `unwrap_or(u8::MAX)`. The old form aliased every
        // slot past 255 onto id 255 — the exact silent-wrong-answer this type
        // exists to prevent, kept alive as a second line of defence that fails
        // the same way. Stopping instead makes `len()` and `iter().count()`
        // disagree, which is visible. `Self::new`'s assert means neither fires.
        self.elements
            .iter()
            .enumerate()
            .map_while(|(i, e)| u8::try_from(i).ok().map(|b| (ElementId(b), e)))
    }

    /// How many elements this universe drew.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.elements.len()
    }

    /// Whether this table has no elements.
    ///
    /// **Not "always `false`"**, which an earlier version of this line claimed:
    /// [`generate_elements`] draws 60..=120, but the in-crate constructor takes
    /// any `Vec`, so an empty table is constructible here today. The
    /// guarantee belongs to the generator, not to the type.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

// **No `IntoIterator for &PeriodicTable`**, deliberately. It was a second copy
// of `iter`'s body, and a public `type IntoIter` has to name a concrete type —
// which pinned the `map`-over-`slice::Iter` shape that reading the id off the
// element required. Two bodies deciding what an iteration item is, where the
// public one pinned the implementation that was wrong. `for (id, e) in
// table.iter()` is one character longer than `for (id, e) in &table` and leaves
// one definition.

/// Convenience for tests and callers that need the universe-domain stream.
#[must_use]
pub const fn stream(seed: u64) -> Stream {
    Stream::new(seed, Domain::Universe, 0)
}

/// Build a universe's periodic table under stated laws.
///
/// **Dispatches on `physics`, exhaustively, the same shape
/// `BondEnergyMatrix::generate` (`pub(crate)`, so not linkable from this
/// public doc comment) dispatches on `table.physics()`.** An earlier
/// version of this function took `physics`
/// and only *stamped* it onto the returned table without branching on it —
/// which was textually compliant with P3 (issue #26's prerequisites, which
/// only named `BondEnergyMatrix::generate` as needing a version-aware path)
/// but left this `pub fn` without the one property that actually matters:
/// today, with one variant, a stamp-only parameter and an exhaustive match
/// are behaviourally identical, so this looked like unrequested structure.
/// The day Task 26.1 adds `PhysicsVersion::V2` with its own
/// electron-configuration derivation, the difference stops being cosmetic —
/// a stamp-only `generate_elements` would keep compiling and keep silently
/// returning V1's derivation for a caller that asked for V2, "in range,
/// plausible and silent," while every other version-aware site in this
/// crate (`From<PhysicsVersion> for u8`, `Universe::generate_under`,
/// `BondEnergyMatrix::generate`, `the_assembled_universe_digest_is_pinned`)
/// would already refuse to compile until V2 is handled. `generate_elements`
/// is `pub`, reachable from any downstream crate, which is exactly where
/// that gap would bite hardest.
///
/// **Draws its own `rung` for the V2 branch — this is the *external*-caller
/// door, which never already has one.** `Universe::generate_v2` (`lib.rs`)
/// draws `rung` before it has a table to build and needs that exact draw to
/// end up on `Universe.rung` and inside the table alike, so it calls
/// `generate_elements_v2` directly with its own draw rather than through
/// here — a `CodeRabbit` finding (2026-08-10) caught the earlier version of
/// this split drawing twice, harmless only because `Rung::draw` is a pure
/// function of `seed` alone; the moment it needed anything else, `Universe.
/// rung` and the table's own configuration could have silently disagreed.
#[must_use]
pub fn generate_elements(seed: u64, physics: PhysicsVersion) -> PeriodicTable {
    match physics {
        PhysicsVersion::V1 => generate_elements_v1(seed, physics),
        PhysicsVersion::V2 => generate_elements_v2(seed, physics, Rung::draw(seed)),
    }
}

/// Find the peak-energy element and set every element's `instability` as its
/// distance below it, in units of the peak's own magnitude. Returns the
/// peak's own `units`.
///
/// **Shared by both physics versions, extracted once the mechanism was
/// copied verbatim into `generate_elements_v2`.** Neither `energy_per_unit`'s
/// formula (version-specific) nor this consequence of it (version-
/// independent: "read the table's own maximum, measure distance from it")
/// belongs to one version — see each generator's own call site for what,
/// if anything, has been separately measured about *its* `energy_per_unit`
/// distribution; this function makes no claim about either.
fn find_peak_and_set_instability(out: &mut [Element]) -> usize {
    let peak = out
        .iter()
        .fold((1_usize, Quanta(f64::MIN)), |(bn, be), e| {
            if e.energy_per_unit > be {
                (e.units, e.energy_per_unit)
            } else {
                (bn, be)
            }
        })
        .0;
    let peak_energy = out
        .get(peak - 1)
        .map_or(Quanta::ZERO, |e| e.energy_per_unit);
    for e in out.iter_mut() {
        // `f64::max` is disallowed — it returns either input on a tie and
        // measured `(+0.0).max(-0.0)` differs between aarch64 and x86-64.
        // `f64::EPSILON` is the *relative* spacing at 1.0, not an absolute
        // magnitude, so using it bare as a floor on an energy is a category
        // error. Harmless here only because `peak_energy` is the series
        // maximum and is O(1) over the drawn grid — named so the guard says
        // what scale it is relative to rather than leaving the next reader
        // to assume.
        const MIN_ENERGY_SCALE: Quanta = Quanta(f64::EPSILON);
        let scale = Quanta(peak_energy.get().abs());
        let scale = if scale > MIN_ENERGY_SCALE {
            scale
        } else {
            MIN_ENERGY_SCALE
        };
        // `Quanta / Quanta` is `f64` by construction in `borbax-units`, so
        // the dimensionlessness of a decay rate is checked rather than
        // assumed. Capped at 1.0 — Task 15 owns the fix to that cap; see
        // `Element::instability`'s own doc for why it is a category error
        // against a Gillespie propensity, deferred rather than overlooked.
        let deficit = (peak_energy - e.energy_per_unit) / scale;
        e.instability = if deficit > 1.0 { 1.0 } else { deficit };
    }
    peak
}

/// **There is no `peak` parameter and there must never be one.** The
/// predecessor took one, computed instability as `((N - peak)/peak)^2`, and
/// reported the resulting minimum as emergent — it was the distance from a
/// typed-in constant, and sweeping `peak` landed the minimum on a closure in
/// 13 of 24 draws. Here one binding energy per unit is derived from the
/// packing counts and both the mass defect and the decay rate are functions
/// of it.
///
/// Infallible: the mass is built in sub-units and handed to `Mass::from_raw`,
/// so there is no range check to fail and no error to propagate. See the task
/// preamble for why a draft made it fallible and why that was wrong.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    reason = "every narrowing here is bounded by construction, and every bound is \
              *computed* by `the_bounds_named_in_the_expect_reason_are_measured` rather \
              than asserted here — three clauses of an earlier version of this string \
              were wrong, which is what a bound stated only in prose is worth. \
              Measured over k in 6..=14, units 1..=120: `period` <= 3, `group` <= 73, \
              `cap` <= 130, `valence` in 0..=6, mass sub-units \
              in 1024..=215040. `out.len()` and `units` are <= `n_elements` <= 120"
)]
fn generate_elements_v1(seed: u64, physics: PhysicsVersion) -> PeriodicTable {
    let mut rng = stream(seed);

    // `k` spans 6..=14. Below 6 a shell cannot triangulate a sphere; above 14
    // the second shell exceeds 200 units, so a 120-element table never leaves
    // it and there are no families at all.
    let k = 6 + rng.next_range(9) as usize;
    // `z = k + 2` is forced by the shell law, not chosen: `shell(1) = k + 2` and
    // shell 1 *is* the units touching the core. An earlier `6 + k/2` agreed only
    // at k = 8 and put 16 sites in shell 1 around a core with 13 neighbours.
    // Construct through `PackingConsts::new` so the relation has one home.
    let pack = PackingConsts::new(k);

    // **`base_mass` must be dyadic and `contact_defect` need not be.** `Mass`
    // is fixed-point at SCALE = 1024; `units * base_mass` is never rounded, so
    // that factor has to land on the grid by itself, while the defect is
    // quantised explicitly below. Probed both ways: drawing `base_mass` in
    // steps of 0.3 fails exactness at N = 1, changing `contact_defect`'s
    // denominator from 512 to 500 changes nothing, and deleting the
    // quantisation fails at N = 2. An earlier comment called both
    // load-bearing, which would have sent a future reader to guard the wrong
    // half.
    let base_mass = 1.0 + 0.25 * rng.next_range(4) as f64;
    let contact_defect = (1 + rng.next_range(4)) as f64 / 512.0;
    let base_radius = 0.30 + 0.02 * rng.next_range(11) as f64;
    let eps = Quanta(0.8 + 0.05 * rng.next_range(9) as f64);
    let sigma = 0.02 + 0.01 * rng.next_range(12) as f64;
    let decay = 0.04 + 0.01 * rng.next_range(13) as f64;
    let n_elements = 60 + rng.next_range(61) as usize; // 60..=120
    // `n_bands` (fed `outer_fill_band`, removed per Task 26.1 Step 8 -- see
    // this file's own header) drew here and nothing else reads further from
    // this stream in this function, so removing the draw outright (not
    // draw-and-discard) shifts nothing downstream: it was the *last* use of
    // `stream(seed)`, and `UniverseConsts`/bonds draw from their own,
    // separate `Domain::Universe` indices (1 and 2). This permanently
    // changes V1's stream shape and its pinned digest -- a deliberate,
    // one-time removal, not the conditional per-seed skip Decision 9's own
    // draw-and-discard fix exists to avoid.

    let mut naming_rng = Stream::new(seed, Domain::Naming, 0);
    let mut taken = Vec::new();
    let mut out: Vec<Element> = Vec::with_capacity(n_elements);
    let mut fallback_symbols = 0_usize;

    for units in 1..=n_elements {
        // Which shell is filling, and how far into it.
        let (mut shell, mut filled) = (0_usize, 1_usize);
        while filled + packing::shell_size(k, shell + 1) <= units {
            shell += 1;
            filled += packing::shell_size(k, shell);
        }
        let outer = units - filled;
        let cap = packing::shell_size(k, shell + 1);
        let contacts = packing::contacts_upto(pack, units);

        let (symbol, name, provenance) = naming::mint(&mut naming_rng, &mut taken);
        if provenance == naming::Provenance::Fallback {
            fallback_symbols += 1;
        }

        // Mass, in sub-units of 1/1024, so the grid is a property of the
        // construction rather than a claim about it. `base_mass` is drawn on a
        // 0.25 step so `base_mass * 1024.0` is an exact integer; the defect is
        // quantised with `round_ties_even` to match `Mass`'s own convention —
        // ties are reachable (13 in this series, at N = 2, 3, 4, 7 among
        // others), and two conventions for one grid means a later "harmonise
        // the rounding" cleanup silently moves a golden.
        let base_sub = (base_mass * 1024.0).round_ties_even() as i64;
        let defect_sub = (contacts * contact_defect * 1024.0).round_ties_even() as i64;
        let mass = Mass::from_raw(units as i64 * base_sub - defect_sub);

        // Radius: set by which shell is occupied, and it cannot shrink.
        //
        // The enclosing radius of a superset of units is never smaller, and an
        // earlier `cbrt(units) * (1 + 0.55*fill)` collapsed it 35% at every
        // closure (N=74 -> 75: 2.587 -> 1.687 on one unit added). Since G1
        // measured binding at ~93% size and ~7% shape, that cliff would have
        // read as closed-shell elements binding a different partner set — an
        // arithmetic accident wearing the look of interesting periodicity.
        //
        // For a shelled cluster the radius *is* the shell index: each shell adds
        // one lattice spacing. So `shell + fill` is both monotone and the better
        // model, and it recovers `N^(1/3)` asymptotically since N ~ k*shell^3/3.
        // Over the range a table actually spans the fitted exponent is
        // 0.29..0.32 against a true cube root's 0.333 — the `1.0 +` offset drags
        // it down. Not a defect; an unqualified range claim would be.
        //
        // G3: real atomic radius falls across a period and jumps at a new one.
        // This rises monotonically and never jumps — structurally unlike,
        // without being geometrically wrong.
        let fill = outer as f64 / cap as f64;
        let radius = base_radius * (1.0 + shell as f64 + fill);

        // Valence: docking notches on the frontier. Zero at a closure because
        // there is no frontier there.
        //
        // `round_ties_even`, for the same reason as the mass defect: `notches`
        // lands on an exact 2.5 at k=8/outer=2 and k=8/outer=8, where ties-away
        // and ties-even differ by a whole bonding slot.
        let valence = packing::frontier_notches(cap, outer).round_ties_even() as u8;

        // Affinity: exposed fraction of the cluster. Continuous across every
        // closure by construction — at `outer == 0` it reduces to the
        // outermost complete shell, and at `outer == cap` to the same
        // expression for the next one.
        let shell_units = if shell == 0 {
            1
        } else {
            packing::shell_size(k, shell)
        };
        let inner = filled - shell_units;
        let covered = shell_units as f64 * outer as f64 / cap as f64;
        let surface = units as f64 - (inner as f64 + covered);
        let affinity = 2.0 * (surface / units as f64) - 1.0;

        // Per-unit binding energy: made contacts, less radial strain.
        //
        // **The strain term is the invented destabiliser, and it was chosen to
        // avoid the semi-empirical mass formula's Coulomb term** — the familiar
        // way to bend a binding curve back down, which imports real nuclear
        // physics with no data file in sight. That reasoning was written under
        // §5's G3, which forbade resemblance to real chemistry absolutely; G3
        // was revised 2026-08-06 to "not necessarily real", so the reasoning is
        // no longer load-bearing on its own. It is kept here because nobody has
        // evaluated a Coulomb-like alternative against this formula's own
        // measured behaviour (below) — that evaluation is issue #26's, not a
        // quick swap, since #26 is the deliberate redesign toward approximated
        // real physics and a term changed in isolation would be tuned blind.
        // **The physical motivation, not the literal computation.** A lattice
        // cannot tile a sphere, so each shell is stretched over a larger
        // radius than the one below and carries a strain growing as `n^2` —
        // that picture is what motivates a strain term at all, and it is not
        // what is computed below. Read on: a per-shell sum of that quantity
        // gives a constant, not `N^(2/3)`; what is actually charged is a
        // single direct per-site penalty on the *total* cluster size, chosen
        // to reach the `2/3` exponent this scheme wants without claiming to
        // derive it from the per-shell picture above.
        //
        // (The defence "no charge, no isospin, no pairing and no asymmetry
        // term" stood here and is deleted: it enumerates the terms that are
        // *absent* and says nothing about the one that is present, which is the
        // shape of the `libm` doctrine error. §7.2's prose forbids it.)
        //
        // **The exponent 2/3 is CHOSEN, not derived, and it is load-bearing.**
        // Per-*shell* strain proportional to `n^2`, summed and divided by `N`,
        // gives `1/k` — a constant, not `N^(2/3)`. Reaching 2/3 needs
        // per-*site* strain proportional to `n^2`, which is a different claim.
        // Accumulating shell-by-shell instead — as `contacts_upto` does, for
        // the identical reason — moves the peak in **48.1% of drawn
        // universes**, so the `19 of 24` closure census is a property of this
        // choice and not of the packing. Deferred to Task 20 deliberately:
        // precedence 2, not a G3 breach, remedy known and it requotes every
        // number. Do not describe the exponent as derived.
        // Wrapped once, at the point the quantity acquires its meaning. `eps`
        // and `sigma` stay bare `f64`: they are the drawn *scale* of a contact
        // and of the strain penalty, and typing them would mint two more energy
        // units with no operations between them. The expression inside the
        // constructor is byte-for-byte the one Task 4 shipped — this is a
        // typing change, and `the_universe_digest_is_pinned` is what says so.
        // Parenthesised to preserve the left-association exactly, so typing
        // `eps` as `Quanta` moves no bits (§13.4).
        let energy_per_unit = (eps * contacts) / units as f64
            - Quanta(sigma) * det_math::cbrt((units * units) as f64);

        // Abundance: fusion builds heavy clusters from light ones, so each
        // extra unit costs a step and abundance falls geometrically. Sequential
        // addition with a constant per-step survival gives `r^N = exp(N ln r)`.
        //
        // **There is deliberately no `exp(beta * energy_per_unit)` tilt.** A
        // draft had one at beta = 0.6. Its stated mechanism compounds to
        // `E_total`, not `E_per_unit`; implementing the stated story makes the
        // heaviest element the most abundant and collapses p_c to 0.112. And
        // 0.6 was typed in where every neighbouring constant is drawn, which is
        // the same shape of defect as `peak`-as-an-argument.
        //
        // **This is the gel lever (§7.2) and Task 4 gates on nothing.** Measured
        // nominal p_c is 0.47..0.91 across 24 universes (median 0.74) — read it
        // from a probe run rather than from this comment. Do NOT justify that
        // with "burial makes realised functionality lower, so realised p_c is
        // higher" — §7.2 carried that argument and it is wrong: under
        // degree-independent thinning the beaker gels at the same bonds per atom
        // regardless, and only the coordinate was rescaled. Verified exactly by
        // binomial convolution: `f_w' = 1 + q(f_w - 1)`, `p_c' = p_c/q`, and
        // bonds-per-node at gel invariant to every printed digit.
        //
        // **Only the *mechanism* is retracted, not the conclusion** — and §7.2
        // restores the conclusion by a different route, so a reader of this file
        // alone would otherwise conclude the opposite of a reader of §7.2 alone.
        // The route that survives is the operating point: §9.1's bimolecular
        // Condense against §9.4's unimolecular Cleave, where the margin between
        // `q = 1.0` and `q = 0.5` is `F(R) = 2*p_ss(R)/p_ss(2R)`, in which `p_c`
        // and `f_w` cancel, leaving it bounded in (2, 4) and monotone. Burial
        // does move the beaker away from gel; the old argument simply was not
        // why. Report nominal f_w
        // and bonds-per-node at gel as diagnostics; the gate is Task 20's
        // battery, on `p_ss / p_c`.
        let abundance = det_math::exp(-(units as f64) * decay);

        out.push(Element {
            symbol,
            name,
            units,
            period: shell as u8,
            group: outer as u8, // `outer < cap` and the table caps at 120
            mass,
            valence,
            affinity,
            radius: Span(radius),
            energy_per_unit,
            instability: 0.0, // filled below, once the peak is known
            abundance,
        });
    }

    // **The peak is read off here, after the series exists.** This ordering is
    // the mechanism: nothing above could have consulted it.
    //
    // Every closure is a *local* maximum by mechanism — a closed shell has no
    // frontier, so it makes every lateral contact available to it and maximises
    // contacts per unit, giving a sawtooth. That is structural: zero failures
    // across the full (k, eps, sigma) product.
    //
    // The **global** max is not. A series truncated mid-shell can peak at its
    // edge, and over the drawn table size it does.
    //
    // **Measured here, over `generate_elements` seeds 0..24:** on a closure in
    // 19 of 24, at the table edge in 2, mid-shell in 3; 13 distinct peak
    // positions; table sizes 61..120.
    //
    // The plan reports the same 19 of 24 and 13 distinct, but splits the
    // remaining five as 4 at the edge and 1 mid-shell. That is not a
    // contradiction and it is worth knowing why: those figures come from
    // `experiments/src/bin/fusion.rs`, which draws its constants from a
    // *different* stream, so its "24 seeds" are 24 different universes. Any
    // per-universe census quoted in this file has to be measured from this
    // function. The aggregate shape agrees; the individual tail does not, and
    // only the aggregate was ever a property of the packing.
    //
    // Both errors have been made here — an earlier version called this a
    // coincidence rate, and its correction called it structural without
    // qualification. Neither is right, which is why the census is reported
    // rather than asserted.
    // Decay rate: distance below the peak, so the most tightly bound elements
    // persist and the extremes decay. One quantity, two consequences — which
    // is what the predecessor's two unlinked encodings could not give.
    //
    // **Cap at 1.0 — and the reason this comment used to give was wrong.
    // Left standing deliberately; Task 15 owns the fix.** It said
    // `instability` is a probability per world-year (§7.1), so 1.0 "is the
    // value the type implies and needs no defence". §7.1 does say that —
    // and Task 15 hands this quantity to the scheduler as a **Gillespie
    // propensity**, which is a rate constant: inverse time, defined on an
    // infinitesimal interval, and *unbounded above* (Gillespie 2007, Eq. 2,
    // `a_j(x) = c_j x_1` for the unimolecular case). Against that consumer a
    // ceiling is not a saturation but a category error, and 200 summed
    // per-interval probabilities for a `MAX_POLYMER` chain is not a
    // probability of anything.
    //
    // **Deferred, not overlooked.** Deleting the cap costs 0.7% of the
    // abundance-weighted mean decay and *gains* 2.1% of the table's
    // max/median ratio spread — which is the quantity §9.4's differential
    // persistence is actually about. It is not applied here because it
    // shares a fix with two larger items now written into Task 15 as
    // requirements 4 and 5: the radiogenic channel carries no scale constant
    // and so dominates cleave by 10^3..10^7, and §7.1's documented *type* is
    // what has to change. One golden regeneration, not three.
    //
    // **The cause is upstream of the cap.** `contacts_upto` returns 0 at
    // `units <= 1`, so `energy_per_unit(1) = -sigma` exactly and
    // `deficit(1) = 1 + sigma/peak_energy > 1` in **every** universe —
    // measured, the monomer's `instability` takes exactly one value across
    // 20 000 universes. The strain term charges a lone unit the full
    // per-site penalty of a lattice whose own comment above says "each shell
    // is stretched over a larger radius than the one below", and at N = 1
    // there is no shell below. Same root cause as the open `N^(2/3)` item.
    //
    // **The 0.9 this replaced was still worse, and that measurement
    // stands** —
    // measured over the 2916-universe grid it binds on 1.19% of elements
    // but on **at least one element in every universe**, and those elements
    // sit at median position N/n_elements = 0.011: it is the monomer, the
    // feedstock, the most abundant species. So a bare 0.9 systematically
    // slowed the decay of the one species whose persistence matters most,
    // by up to 26%, as a side effect of a guard nobody described. The raw
    // deficit does exceed 1 (max measured 1.218), so a cap is genuinely
    // needed — it is the value that was wrong, not the clamp.
    //
    // If a future measurement wants the monomer to persist longer, that is
    // a physics claim about feedstock and belongs in §9.4 with the
    // measurement attached, not in a magic number here. And the feedstock
    // defence does not hold in V0 regardless: there is no inflow — the
    // `Beaker` API is `{new, step, time, counts, species_count}` and §11's
    // open boundaries are V1 — so a closed beaker losing its most abundant
    // element does not turn over, it empties.
    let peak = find_peak_and_set_instability(&mut out);

    let shell = ShellPattern {
        law: ShellLaw::V1 { k },
        eps,
        closures: packing::closures(k, n_elements),
        fallback_symbols,
        peak,
    };
    PeriodicTable::new(shell, out, physics)
}

/// `tanh`, built from [`det_math::exp`] rather than a native `f64::tanh` —
/// §13.1 routes every transcendental through `det_math`, and `tanh` is not
/// one of the primitives it exposes yet, so this is the one caller building
/// it from a primitive that *is* covered rather than reaching for a
/// platform tanh. **Branches on the sign of `x` for numerical safety, not
/// style**: the naive `(exp(2x)-1)/(exp(2x)+1)` overflows to `inf/inf` =
/// `NaN` for large positive `x` (`exp` saturates to `f64::INFINITY`); using
/// `exp(-2*|x|)` on whichever side keeps the exponent non-positive, so the
/// exponential only ever underflows harmlessly toward `0.0`, never
/// overflows.
fn tanh(x: f64) -> f64 {
    if x >= 0.0 {
        let e = det_math::exp(-2.0 * x);
        (1.0 - e) / (1.0 + e)
    } else {
        let e = det_math::exp(2.0 * x);
        (e - 1.0) / (e + 1.0)
    }
}

/// Decision 8's affinity centring: the table's own median of `-e_homo` (an
/// ionisation-energy-like, always-positive quantity), read off the finished
/// series — the same "computed constant" principle `peak` already follows,
/// never passed in — plus a self-scaling denominator (the median itself,
/// floored near zero the same way `find_peak_and_set_instability` guards
/// `peak_energy`). Deterministic median: sort via `canonical_cmp`, never
/// `partial_cmp().unwrap()`, then index `n/2`.
fn median_neg_e_homo_and_scale(props: &[ElectronicProperties], n_elements: usize) -> (f64, f64) {
    const MIN_AFFINITY_SCALE: f64 = f64::EPSILON;
    let mut neg_e_homo: Vec<f64> = (0..n_elements)
        .map(|i| {
            -props
                .get(i)
                .unwrap_or_else(|| unreachable!("n_elements <= 120, well within ELECTRON_CEILING"))
                .e_homo
        })
        .collect();
    neg_e_homo.sort_by(|a, b| canonical_cmp(*a, *b));
    let median = neg_e_homo
        .get(n_elements / 2)
        .copied()
        .unwrap_or_else(|| unreachable!("n_elements >= 60, so n_elements / 2 is in bounds"));
    let scale = if median.abs() > MIN_AFFINITY_SCALE {
        median.abs()
    } else {
        MIN_AFFINITY_SCALE
    };
    (median, scale)
}

/// V2's element generator (issue #26, Task 26.1 Steps 6-9, Decision 8) —
/// screened-hydrogenic orbital fill (`orbital.rs`) replaces V1's
/// shell-packing model entirely. `symbol`/`name`/`units`/`mass`/`valence`/
/// `affinity`/`radius`/`energy_per_unit`/`instability`/`abundance` name the
/// same *kind* of quantity as V1's; `period`/`group` do not share an
/// indexing or value convention across versions — see their own field docs
/// on [`Element`].
///
/// **`n_elements` is Decision 9's draw-and-discard, not a new mechanism.**
/// The word is drawn off the stream identically at every rung; only at
/// `rung == 0` is its *result* discarded in favour of 118 (the real
/// periodic table's size), so every later draw's stream position is
/// unaffected by which rung this universe landed on.
///
/// **Not a thirteenth `MigratedConstant`, and that is a considered
/// rejection, not an oversight.** `draw_symmetric` perturbs continuously
/// around a base and lands on that base *only* at `rung == 0` — but the
/// plan's own Decision 9 fixes this constant's non-identity range at a
/// uniform `60..=120` draw, unrelated to how far the rung sits from
/// identity, and continuous perturbation around 118 cannot produce a
/// uniform distribution on that interval at every other rung. Migrating
/// this constant would either narrow the drawn range or make it
/// rung-dependent, both of which the plan explicitly forecloses: "Its
/// drawn range stays `60..=120`, unchanged... this constant does *not* go
/// through P7's multiplicative mechanism."
///
/// **`decay` (feeds `abundance`) and `radius_scale`/`base_mass`/
/// `contact_defect`'s base values are chosen, not derived — `Counterpart::DefaultOnly`
/// for all of them.** `decay` is not named in the plan's own P7 migration
/// table at all, so it stays a plain per-universe draw, same mechanism and
/// shape V1 already uses, just V2's own value on V2's own stream position.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "every narrowing here is bounded by construction: `n_elements` is 60..=120 \
              (Decision 9), `outer_n`/`frontier_l` are small quantum numbers bounded by \
              orbital.rs's own MAX_N = 9, `valence` is a Hund's-rule unpaired count or the \
              promotion bonus (2), never large, and mass sub-units follow the same \
              round_ties_even -> i64 pattern generate_elements_v1 already uses. No \
              cast_sign_loss here (unlike generate_elements_v1): every narrowing in this \
              function goes to a signed destination wide enough for its bounded input, never \
              float/signed -> unsigned"
)]
pub(crate) fn generate_elements_v2(
    seed: u64,
    physics: PhysicsVersion,
    rung: Rung,
) -> PeriodicTable {
    let mut rng = stream(seed);
    let raw_n_elements = 60 + rng.next_range(61) as usize;
    let n_elements = if rung.is_identity() {
        118
    } else {
        raw_n_elements
    };
    let decay = 0.04 + 0.01 * rng.next_range(13) as f64;

    let orbital_consts = OrbitalConsts::draw(seed, rung);
    let gap_consts = GapConsts::draw(seed, rung);
    let runs = orbital::fill();
    let props = orbital::electronic_properties(&runs, orbital_consts, gap_consts);

    let radius_scale = draw_symmetric(seed, rung, MigratedConstant::RadiusScale, 0.4);
    let base_mass = draw_symmetric(seed, rung, MigratedConstant::BaseMass, 1.375);
    // Re-quantised to the 1/1024 grid immediately after perturbation -- the
    // continuous P7 multiply does not itself land on the grid, and only the
    // quantised value may reach `Mass::from_raw` (see element.rs:486-494's
    // own note on why `base_mass` must stay dyadic, restated for V2 here).
    let base_mass_sub = (base_mass * 1024.0).round_ties_even() as i64;
    let contact_defect = draw_symmetric(seed, rung, MigratedConstant::ContactDefect, 0.005);

    let (median_neg_e_homo, affinity_scale) = median_neg_e_homo_and_scale(&props, n_elements);

    let identity_witness = naming::IdentityWitness::new(rung);
    let mut naming_rng = Stream::new(seed, Domain::Naming, 0);
    let mut taken = Vec::new();
    let mut fallback_symbols = 0_usize;

    let mut out: Vec<Element> = Vec::with_capacity(n_elements);
    let mut cumulative_e_homo = 0.0_f64;
    for units in 1..=n_elements {
        let p: ElectronicProperties = *props
            .get(units - 1)
            .unwrap_or_else(|| unreachable!("n_elements <= 120, well within ELECTRON_CEILING"));
        cumulative_e_homo += p.e_homo;

        let (symbol, name) = identity_witness.map_or_else(
            || {
                let (symbol, name, provenance) = naming::mint(&mut naming_rng, &mut taken);
                if provenance == naming::Provenance::Fallback {
                    fallback_symbols += 1;
                }
                (symbol, name)
            },
            |witness| {
                let (s, nm) = naming::real_name(witness, units).unwrap_or_else(|| {
                    unreachable!(
                        "units <= n_elements == 118 at identity, within REAL_TABLE's 1..=118 range"
                    )
                });
                (s.to_string(), nm.to_string())
            },
        );

        let period = p.outer_n;
        let group = p.frontier_l;

        let radius = radius_scale * f64::from(p.outer_n) * f64::from(p.outer_n) / p.zeff_outer;

        let defect_sub =
            (f64::from(p.paired_count) * contact_defect * 1024.0).round_ties_even() as i64;
        let mass = Mass::from_raw(units as i64 * base_mass_sub - defect_sub);

        // `cumulative_e_homo` is a running left-fold in ascending-index order --
        // bit-identical to re-summing `props[0..units]` fresh every iteration,
        // since both add the same terms in the same order, just without the
        // O(n) rewalk this loop would otherwise repeat n_elements times.
        let mean_e = cumulative_e_homo / units as f64;
        let energy_per_unit = Quanta(-mean_e);

        let affinity = tanh((-p.e_homo - median_neg_e_homo) / affinity_scale);

        let abundance = det_math::exp(-(units as f64) * decay);

        out.push(Element {
            symbol,
            name,
            units,
            period,
            group,
            mass,
            valence: p.valence as u8,
            affinity,
            radius: Span(radius),
            energy_per_unit,
            instability: 0.0, // filled below, once the peak is known
            abundance,
        });
    }

    // Peak-finding and deficit: shared with V1 via `find_peak_and_set_instability`
    // -- the "table's own energy scale, read off the finished series"
    // principle is version-independent, only what feeds `energy_per_unit`
    // differs. Neither the 19-of-24-closures census nor the monomer/0.9-cap
    // measurements V1's own call site documents have been separately
    // re-measured for V2's own `energy_per_unit` distribution.
    let peak = find_peak_and_set_instability(&mut out);
    let peak_energy = out
        .get(peak - 1)
        .map_or(Quanta::ZERO, |e| e.energy_per_unit);

    // Period boundaries: every z whose outer_n differs from z+1's -- z is
    // the LAST element of the closing period, matching V1's own "closures
    // are the last element of a completed shell" convention. No test reads
    // V2's `closures` yet (this crate's own tests only read it through
    // `table()`, which fixes `PhysicsVersion::CURRENT` == V1), so this is a
    // straightforward, not yet independently scrutinised, derivation.
    let closures: Vec<usize> = (1..n_elements)
        .filter(|&z| props.get(z - 1).map(|q| q.outer_n) != props.get(z).map(|q| q.outer_n))
        .collect();

    let shell = ShellPattern {
        law: ShellLaw::V2,
        eps: peak_energy,
        closures,
        fallback_symbols,
        peak,
    };
    PeriodicTable::new(shell, out, physics)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Named for brevity below; `generate_elements` is infallible because the
    /// mass is built in sub-units (see its doc), so there is nothing to unwrap.
    ///
    /// `PhysicsVersion::CURRENT`, not a hardcoded `V1`: this helper predates
    /// versioning and every caller below wants "the shipped table", the same
    /// intent [`crate::Universe::generate`] follows.
    fn table(seed: u64) -> (ShellPattern, Vec<Element>, PhysicsVersion) {
        generate_elements(seed, PhysicsVersion::CURRENT).into_parts()
    }

    /// V2's own table, for Task 26.1 Step 6's acceptance tests below --
    /// `table(seed)` deliberately fixes `PhysicsVersion::CURRENT`, which is
    /// V1, so every V2-specific assertion needs this instead.
    fn table_v2(seed: u64) -> (ShellPattern, Vec<Element>, PhysicsVersion) {
        generate_elements(seed, PhysicsVersion::V2).into_parts()
    }

    /// Decision 8's "one shared computation" claim, checked directly rather
    /// than trusted: `period`, `group` and `valence` must come from the
    /// *same* [`ElectronicProperties`] entry per `z`, not four
    /// independently-tuned formulas that happen to agree on this input. A
    /// version wiring `period` from a stale or re-derived quantity would
    /// pass every range/shape test above while still failing this one,
    /// since it compares against the source computation directly.
    #[test]
    fn period_group_and_valence_come_from_the_same_electronic_properties_entry() {
        for seed in [0u64, 1, 5, 21, 42] {
            let rung = Rung::draw(seed);
            let runs = orbital::fill();
            let props = orbital::electronic_properties(
                &runs,
                OrbitalConsts::draw(seed, rung),
                GapConsts::draw(seed, rung),
            );
            let (_, els, _) = table_v2(seed);
            for (i, e) in els.iter().enumerate() {
                let p = props
                    .get(i)
                    .unwrap_or_else(|| unreachable!("els.len() <= n_elements <= props.len()"));
                assert_eq!(e.period, p.outer_n, "seed {seed} z={}: period", i + 1);
                assert_eq!(e.group, p.frontier_l, "seed {seed} z={}: group", i + 1);
                assert_eq!(
                    e.valence,
                    u8::try_from(p.valence)
                        .unwrap_or_else(|_| unreachable!("valence is a small unpaired count")),
                    "seed {seed} z={}: valence",
                    i + 1
                );
            }
        }
    }

    /// `radius(z) = radius_scale * outer_n^2 / zeff_outer` falls within a
    /// period (since `zeff_outer` rises there — `orbital.rs`'s own
    /// `zeff_outer_rises_within_a_period` pins the mechanism) and jumps up
    /// at every period start, with the boundary read from `period` itself
    /// (`outer_n`), never a literal list of z values.
    #[expect(
        clippy::indexing_slicing,
        reason = "`windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn v2_radius_falls_within_a_period_and_jumps_at_period_starts() {
        for seed in [0u64, 1, 5, 21, 42] {
            let (_, els, _) = table_v2(seed);
            for w in els.windows(2) {
                let (prev, next) = (&w[0], &w[1]);
                if next.period == prev.period {
                    // **A small relative tolerance, not exact `<=`, at
                    // frontier subshell transitions.** `zeff_outer`'s probe
                    // changes subshell at an s->p (or d->f, ...) boundary
                    // within a period (this module's `orbital.rs` doc,
                    // finding F4's fix), and the +1 electron there can
                    // near-exactly cancel against the newly-counted
                    // same-shell screening term — measured, at the identity
                    // configuration the cancellation is bit-exact, but a
                    // perturbed universe's own arithmetic can land a few ULP
                    // on either side of the tie (seed 1, units 67->68:
                    // 5.0921313321890045 -> 5.092131332189016, a 1.1e-14
                    // difference). `64 * f64::EPSILON`, relative to the
                    // magnitude, is generous against that noise floor and
                    // would still fail on a genuine rise, which is orders of
                    // magnitude larger everywhere this test has measured one.
                    let tolerance = 64.0 * f64::EPSILON * prev.radius.0.abs();
                    assert!(
                        next.radius.0 <= prev.radius.0 + tolerance,
                        "seed {seed} units {}->{}: radius should not rise within a period ({} -> {}, \
                         tolerance {tolerance})",
                        prev.units,
                        next.units,
                        prev.radius.0,
                        next.radius.0
                    );
                } else {
                    assert!(
                        next.radius.0 > prev.radius.0,
                        "seed {seed} units {}->{}: radius should jump up at a new period ({} -> {})",
                        prev.units,
                        next.units,
                        prev.radius.0,
                        next.radius.0
                    );
                }
            }
        }
    }

    /// #31's closure: `affinity` must be periodic in Z, not a monotone
    /// trend. **The sign-change count is the real discriminator, and it
    /// must be a count, not a boolean.** A monotone function can cross zero
    /// at most once — so "at least one sign change" is satisfied by a
    /// monotone function that happens to cross zero, and proves nothing.
    /// "At least two" is the actual threshold no monotone function can
    /// reach, since a second crossing requires a local extremum in
    /// between. `MIN_SIGN_CHANGES = 3` sits one above that mathematical
    /// floor and a long way below **measured** counts of 10-17 within a
    /// period across `[0, 1, 5, 21, 42]` — real margin in both directions,
    /// not the bare minimum that would merely fail to be wrong. The
    /// correlation bound is supplementary and **measured, not guessed**:
    /// `-e_homo` genuinely trends upward with Z (deeper, more-screened
    /// shells still bind more tightly on average), so affinity correlates
    /// 0.86-0.91 with Z across the same seeds — real periodic structure
    /// riding on top of a real trend, not a flat oscillation. An earlier
    /// draft of this bound guessed 0.6 before measuring and failed on
    /// every seed. `0.95` is pre-fixed with real margin above the measured
    /// range and stays a long way short of the ~1.0 a literal Z-surrogate
    /// (what #31's old formula effectively was) would hit.
    #[expect(
        clippy::indexing_slicing,
        reason = "`windows(2)` yields slices of exactly length 2"
    )]
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "`els.len()` is at most 120, well within f64's exact integer range"
    )]
    #[test]
    fn v2_affinity_is_periodic_not_a_monotone_trend_in_z() {
        const CORRELATION_BOUND: f64 = 0.95;
        const MIN_SIGN_CHANGES: usize = 3;
        for seed in [0u64, 1, 5, 21, 42] {
            let (_, els, _) = table_v2(seed);
            let mut sign_changes_within_a_period = 0_usize;
            for w in els.windows(2) {
                if w[1].period == w[0].period && w[0].affinity.signum() != w[1].affinity.signum() {
                    sign_changes_within_a_period += 1;
                }
            }
            assert!(
                sign_changes_within_a_period >= MIN_SIGN_CHANGES,
                "seed {seed}: affinity changes sign only {sign_changes_within_a_period} time(s) \
                 within a period -- too few to rule out a monotone trend that happens to cross zero"
            );

            let z: Vec<f64> = (1..=els.len()).map(|u| u as f64).collect();
            let a: Vec<f64> = els.iter().map(|e| e.affinity).collect();
            let corr = pearson_correlation(&z, &a);
            assert!(
                corr.abs() < CORRELATION_BOUND,
                "seed {seed}: affinity correlates {corr:.4} with Z -- too close to a monotone \
                 trend for a periodic quantity"
            );
        }
    }

    /// Pearson correlation, hand-rolled: this file's own established
    /// precedent for small numeric helpers (`tanh`, above) is to write the
    /// ~10 lines rather than reach for a dependency, and a test-only
    /// helper is an even lower bar than that (dev-dependencies are cheap,
    /// but this needs nothing a crate would give beyond what fits here).
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "x.len() is a small test-fixture size, well within f64's exact integer range"
    )]
    fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
        let n = x.len() as f64;
        let mean_x = x.iter().sum::<f64>() / n;
        let mean_y = y.iter().sum::<f64>() / n;
        let mut cov = 0.0;
        let mut var_x = 0.0;
        let mut var_y = 0.0;
        for (&xi, &yi) in x.iter().zip(y) {
            let dx = xi - mean_x;
            let dy = yi - mean_y;
            cov += dx * dy;
            var_x += dx * dx;
            var_y += dy * dy;
        }
        cov / (var_x.sqrt() * var_y.sqrt())
    }

    /// Decision 8's headline valence-promotion claim, checked again at the
    /// level the wiring actually ships -- `orbital.rs`'s own
    /// `valence_matches_real_ground_state_valences_through_period_3` checks
    /// the *unwired* `ElectronicProperties`; this checks `Element::valence`
    /// as `generate_elements_v2` actually stores it, at the identity
    /// configuration where the real analogues are nameable.
    #[test]
    fn v2_valence_promotion_matches_beryllium_and_helium_analogues_through_the_wired_element() {
        let (_, els, _) = table_v2(21); // identity configuration (see lib.rs's own pinned seed)
        let helium = els
            .get(1)
            .unwrap_or_else(|| unreachable!("z=2 exists in any 60..=120 table"));
        assert_eq!(
            helium.valence, 0,
            "z=2 (helium-analogue) should not promote"
        );
        let beryllium = els
            .get(3)
            .unwrap_or_else(|| unreachable!("z=4 exists in any 60..=120 table"));
        assert_eq!(
            beryllium.valence, 2,
            "z=4 (beryllium-analogue) should promote to valence 2"
        );
    }

    /// The `n_elements`-at-a-boundary property this file's own design memo
    /// resolved as identity-only, not a general per-seed guarantee (every
    /// `MigratedConstant` lands on its base value at `rung == 0`, so the
    /// fold itself is seed-independent at identity): real element 118
    /// (Oganesson-analogue) completes its own subshell -- `gap > 0` at
    /// z=118 -- rather than sitting mid-fill, which is what makes 118 a
    /// structurally sound table size to have pinned rather than an
    /// arbitrary cut.
    #[test]
    fn element_118_completes_its_own_subshell_at_the_identity_configuration() {
        let runs = orbital::fill();
        let props = orbital::electronic_properties(
            &runs,
            OrbitalConsts::at_identity(),
            GapConsts::at_identity(),
        );
        let last = props
            .get(117) // z=118, 0-indexed
            .unwrap_or_else(|| unreachable!("the identity fold reaches well past z=118"));
        assert!(
            last.gap > 0.0,
            "z=118 should complete its own subshell (gap > 0), not sit mid-fill"
        );
    }

    #[test]
    fn block_maps_l_to_the_named_blocks_and_falls_back_past_f() {
        assert_eq!(block(0), Block::S);
        assert_eq!(block(1), Block::P);
        assert_eq!(block(2), Block::D);
        assert_eq!(block(3), Block::F);
        assert_eq!(block(4), Block::Other(4));
        assert_eq!(block(8), Block::Other(8));
    }

    /// The maximum valence of a table at an explicit `(k, n_elements)`.
    ///
    /// **This is the whole input space for that quantity.** `valence` is
    /// `frontier_notches(cap, outer).round_ties_even()`, and `cap`/`outer` are
    /// functions of `k` and `units` alone — none of `base_mass`,
    /// `contact_defect`, `decay`, `eps` or `sigma` enters. So enumerating
    /// `9 x 61` cells is exhaustive where sampling seeds is not.
    ///
    /// The shell walk below mirrors `generate_elements` exactly; keep them in
    /// step, or the enumeration stops measuring the thing it names. No
    /// `PackingConsts` is built because none is needed — `z` never enters a
    /// valence, and constructing one here would imply otherwise.
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "`frontier_notches` is non-negative and bounded well under 255 — \
                  `unmade_lateral_is_finite_and_non_negative_everywhere` pins the \
                  first and this test's own assertion pins the second"
    )]
    fn valence_ceiling_at(k: usize, n_elements: usize) -> u8 {
        (1..=n_elements)
            .map(|units| {
                let (mut shell, mut filled) = (0_usize, 1_usize);
                while filled + packing::shell_size(k, shell + 1) <= units {
                    shell += 1;
                    filled += packing::shell_size(k, shell);
                }
                let cap = packing::shell_size(k, shell + 1);
                packing::frontier_notches(cap, units - filled).round_ties_even() as u8
            })
            .max()
            .unwrap_or(0)
    }

    #[test]
    fn generation_is_deterministic() {
        // Purity, not determinism: `RandomState` is seeded once per process,
        // so two identically-built `HashMap`s iterate identically here and
        // would pass. Cross-platform reproducibility is the §13.6 golden
        // matrix's job and nothing in this file can stand in for it.
        assert_eq!(table(42), table(42));
    }

    /// Symbol and name must not be swapped on the way into `Element`.
    ///
    /// `let (symbol, name, _) = naming::mint(..)` feeds `Element { symbol,
    /// name, .. }` by field-init shorthand, so `let (name, symbol, _)` compiles
    /// and silently swaps them. `symbols_are_well_formed` tests `mint`
    /// directly and never sees it; nothing else checks either field at element
    /// level. Two `String`s in a tuple carry no type-level distinction between
    /// a 1–2 character symbol and a 4+ character name, so the assertion has to
    /// supply it.
    #[test]
    fn symbol_and_name_are_not_transposed() {
        for seed in 0..8 {
            for e in table(seed).1 {
                assert!(
                    (1..=2).contains(&e.symbol.len()),
                    "seed {seed}: symbol {:?} is not 1-2 chars",
                    e.symbol
                );
                assert!(
                    e.name.len() >= 4,
                    "seed {seed}: name {:?} is shorter than any minted name",
                    e.name
                );
            }
        }
    }

    #[test]
    fn different_seeds_give_different_tables() {
        assert_ne!(table(1), table(2));
    }

    /// **A cheap behavioural backstop, not the guarantee.**
    ///
    /// Zero-valence-at-closure is held by two independent factors, so the
    /// deletion test cannot tell a formula from an `if closed { 0 }`.
    /// Continuity was the proposed replacement and is **also insufficient**: it
    /// catches a declaration laid over a formula returning something *large*
    /// (`notches + 3` fails it) but not one writing 0 over a formula returning
    /// ~1, because both neighbours of a closure are pinned at valence 1 and the
    /// declaration sits inside the bound. Measured: a mutation putting every
    /// closure at valence 1 passes the entire suite.
    ///
    /// The guarantee is `xtask`'s source check; the output is guarded by
    /// `closed_shells_have_valence_zero`. This is the third line of defence,
    /// and the cheapest.
    #[expect(
        clippy::indexing_slicing,
        reason = "`windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn valence_is_continuous_across_closures() {
        for seed in 0..12 {
            let (_, els, _physics) = table(seed);
            for pair in els.windows(2) {
                let jump = i16::from(pair[1].valence) - i16::from(pair[0].valence);
                assert!(
                    jump.abs() <= 2,
                    "seed {seed} N={} -> {}: valence jumped {jump}",
                    pair[0].units,
                    pair[1].units
                );
            }
        }
    }

    /// The family §7.1 asks for, with a cause behind it. The element one unit
    /// past any closure is a closed core plus a single adhering unit, so it has
    /// exactly one docking notch, in every shell of every universe.
    ///
    /// **Its stated cause is not its cause**, and the correction is worth
    /// keeping: the valence is 1 by rounding `6/(6 - 12/cap)` in [1.14, 1.33],
    /// not because a lone site has "all six lateral contacts unmade" — the
    /// shell offers 4.5 to 5.9. Over the caps actually reached the ratio spans
    /// **[1.0156, 1.3333]**; an earlier "[1.14, 1.33]" was the range over
    /// shell-1 caps only, while this test iterates `closures`, which reaches
    /// shells 2 and 3 and caps to 130. Margin to 1.5 is comfortable either way,
    /// so the family holds.
    #[test]
    fn elements_one_past_a_closure_all_have_valence_one() {
        for seed in 0..12 {
            let (shell, els, _physics) = table(seed);
            for &c in &shell.closures {
                if let Some(e) = els.get(c) {
                    assert_eq!(
                        e.valence, 1,
                        "seed {seed}: N={} follows closure {c}",
                        e.units
                    );
                }
            }
        }
    }

    /// **Valence 8 is not reachable in a table this scheme generates, and
    /// saying so is the finding.**
    ///
    /// A per-*site* count caps at `z - 6`; summing over the frontier removes
    /// that ceiling in principle, and valence 9 does appear once a fourth shell
    /// opens past N = 200. But `n_elements` is drawn on 60..=120, and inside
    /// that range the ceiling is **4 to 6**.
    ///
    /// **Assert `4..=6`, not `4..=7`.** This assertion has been wrong three
    /// times: `5..=7` from a guess, which failed with a 4 in it; widened to
    /// `4..=7`, which was then loose enough to survive `FRONTIER_COEFF` being
    /// refitted from 7.30 to 6.90 — the refit moved the measured ceiling and
    /// the test did not notice. A range assertion that survives a change to the
    /// quantity it measures is not measuring it. Do not widen it back.
    ///
    /// A draft asserted `reached >= 12` of 24. The true rate is 273/549 = 49.7%,
    /// so the threshold sat on the population mean and the test would have
    /// failed on ~43% of seed sets — and its guidance sent the implementer to
    /// widen `k` when the controlling term is `n_elements`.
    ///
    /// **Enumerate the grid; do not sample seeds.** `valence` reads only `cap`
    /// and `outer`, both functions of `k` and `units` alone, so the ceiling is
    /// a pure function of `(k, n_elements)` — 9 x 61 = 549 cells, about a
    /// millisecond. The *third* draft of this assertion still sampled 48 seeds,
    /// and was still blind on the high side: at `FRONTIER_COEFF = 7.20`, inside
    /// the band `frontier_matches_exact_counts_on_a_real_sphere` admits, the
    /// ceiling reaches 7 in 7 of 549 cells and 48 draws miss all seven **54% of
    /// the time**. §7.1 puts valence at 0..6, so a 7 is a spec breach and this
    /// is the only test positioned to catch it. Raising 48 to 480 shrinks the
    /// probability instead of removing it and is not a fix.
    ///
    /// Assert the attained **set**, not containment — a refit collapsing every
    /// universe to 5 satisfies `.all(|v| (4..=6).contains(v))`.
    #[test]
    fn valence_ceiling_is_four_to_six_over_the_drawn_range() {
        let mut seen = std::collections::BTreeSet::new();
        for k in 6..=14 {
            for n in 60..=120 {
                seen.insert(valence_ceiling_at(k, n));
            }
        }
        let expected: std::collections::BTreeSet<u8> = [4, 5, 6].into_iter().collect();
        assert_eq!(
            seen, expected,
            "the attained valence ceiling set moved over the full drawn grid"
        );

        // **Ties `BondOrder::MAX` to the measurement that defines it.** The cap
        // was an undrawn 3 while this set was `{4,5,6}` — the emergence lane's
        // finding — and the repair raised it to 6. Nothing pointed the two at
        // each other, so a refit that moved this set would have left the cap
        // silently below the generated ceiling again, in the other file. This is
        // exhaustive over the drawn grid and costs a millisecond.
        assert!(
            seen.iter().all(|&v| v <= crate::bonds::BondOrder::MAX),
            "the valence ceiling {seen:?} now exceeds BondOrder::MAX = {} — raise \
             the cap in bonds.rs, or the one bonding rule in this crate is again \
             below the ceiling the table generates",
            crate::bonds::BondOrder::MAX
        );
    }

    /// §7.1: mass increases down the table. The fitted series this replaced
    /// stopped increasing at N = 124 — so the third closed shell weighed less
    /// than the element below it — and a 120-element table did not show it.
    #[expect(
        clippy::indexing_slicing,
        reason = "`windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn mass_increases_with_index() {
        for seed in 0..12 {
            let (_, els, _physics) = table(seed);
            for w in els.windows(2) {
                assert!(
                    w[1].mass > w[0].mass,
                    "seed {seed}: {} !> {}",
                    w[1].name,
                    w[0].name
                );
            }
        }
    }

    /// **`Mass::from_raw` validates nothing** — its own doc says so, and `impl
    /// Add` wraps in release on an out-of-range value. Since the construction is
    /// now integer, grid-exactness is a property of the type rather than a claim
    /// about the series; what still needs checking is that the drawn range
    /// cannot leave `MAX_MAGNITUDE`.
    ///
    /// Three drafts of this test were wrong in three ways. The first compared
    /// `Mass` to `Result<Mass, _>` and did not compile. The second round-tripped
    /// an already-quantised value — `to_f64` is an exact division by a power of
    /// two and `from_f64_quantised` an exact multiplication back, so it could
    /// not fail for anything the construction can produce. The third *described*
    /// the repair in this comment and left the round-trip in the body, which is
    /// how it reached review.
    ///
    /// **Assert the raw sub-unit count against the tight bound**, and the
    /// reason is *tightness* — nothing to do with overflow. A fourth draft
    /// justified it by overflow; the overflow needs `units >= 2^63/1792 ~ 5.2e15`
    /// against a loop bound of 120, unreachable by thirteen orders of magnitude.
    ///
    /// What the tight bound buys: a draw range that moves fails *here*, instead
    /// of silently spending 4993x of headroom. `MAX_MAGNITUDE` would not notice.
    #[test]
    fn every_mass_is_inside_the_range_from_raw_does_not_check() {
        // `const`, not `checked_mul(..).expect(..)`: const evaluation makes an
        // overflowing product a *compile* error (E0080), which is strictly
        // stronger than a runtime check, and avoids an `expect` that
        // `clippy::expect_used` denies — `clippy.toml` sets neither
        // `allow-expect-in-tests` nor `allow-unwrap-in-tests`, so the deny
        // reaches inside `#[cfg(test)]`.
        const MAX_UNITS: i64 = 120;
        const MAX_BASE_SUB: i64 = 1792;
        const CEILING: i64 = MAX_UNITS * MAX_BASE_SUB;
        for seed in 0..24 {
            let (_, els, _physics) = table(seed);
            for e in els {
                let raw = e.mass.raw();
                assert!(
                    (1..=CEILING).contains(&raw),
                    "seed {seed} N={}: raw {raw} outside 1..={CEILING} — a draw range moved",
                    e.units
                );
            }
        }
    }

    /// A closed shell has no frontier, so it has no unmade lateral contacts.
    /// This is the mechanism's headline output and what §7.2's "the most
    /// tightly bound element cannot bond" reading rests on.
    ///
    /// **It was in the probe and not carried across, and nothing else covers
    /// it.** Measured: a mutation breaking both zeroing factors puts every
    /// closure at valence 1 and the whole suite still passes — continuity
    /// trivially (1→1), the ceiling test, and `elements_one_past_a_closure`.
    ///
    /// It is *not* the test that distinguishes mechanism from declaration; two
    /// independent factors hold the zero, so a declaration can be laid over a
    /// broken formula and this still passes. That claim is a property of the
    /// source — see the `xtask` check.
    #[test]
    fn closed_shells_have_valence_zero() {
        for seed in 0..24 {
            let (sp, els, _physics) = table(seed);
            let ShellLaw::V1 { k } = sp.law else {
                unreachable!("table() fixes PhysicsVersion::CURRENT, which is V1")
            };
            for n in packing::closures(k, els.len()) {
                if let Some(e) = els.get(n - 1) {
                    assert_eq!(e.valence, 0, "seed {seed}: N={n} is a closure");
                }
            }
        }
    }

    /// Affinity is `2 * exposed_fraction - 1`, so it cannot leave `(-1, 1]`
    /// unless burial is miscounted. The predecessor reached exactly +1.000 at
    /// every closure by reading the capacity of the *empty* shell, making the
    /// stable elements the stickiest and stopping anything accumulating.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        clippy::indexing_slicing,
        reason = "`units` is at most 120; `windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn affinity_is_in_range_and_continuous() {
        // **Assert the attained range, both sides.** The predecessor asserted
        // `affinity > -1.0`, which has 1.14 of slack and cannot fire — the
        // defect class this file names 65 lines further down. `affinity` is a
        // pure function of `(k, units)` (no drawn constant enters it), so the
        // space is enumerable and the attained bounds are exact.
        //
        // **The lower bound being positive is a finding, not a target.** §8.3's
        // charge term `-(a_A + a_B)^2` is maximised at zero when the two
        // affinities are *opposite*, so a strictly positive range degenerates it
        // from a complementarity test into a monotone penalty on total surface
        // affinity, and §10.1's membrane mechanism — one flank positive, the
        // opposite negative — is unreachable. How `a` is summed from element
        // affinities is Task 7's code and is not written yet, so the fix may
        // belong there; what must not happen is the range being inherited by
        // default. This assertion makes it fail loudly the moment it moves.
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for k in 6..=14_usize {
            for units in 1..=120_usize {
                let (mut shell, mut filled) = (0_usize, 1_usize);
                while filled + packing::shell_size(k, shell + 1) <= units {
                    shell += 1;
                    filled += packing::shell_size(k, shell);
                }
                let cap = packing::shell_size(k, shell + 1);
                let shell_units = if shell == 0 {
                    1
                } else {
                    packing::shell_size(k, shell)
                };
                let inner = filled - shell_units;
                let outer = units - filled;
                let covered = shell_units as f64 * outer as f64 / cap as f64;
                let surface = units as f64 - (inner as f64 + covered);
                let a = 2.0 * (surface / units as f64) - 1.0;
                if a < lo {
                    lo = a;
                }
                if a > hi {
                    hi = a;
                }
            }
        }
        assert!(
            (lo - 0.140_476).abs() < 1e-6 && (hi - 1.0).abs() < 1e-12,
            "the attained affinity range moved: {lo:.6}..={hi:.6}. It has never \
             reached the negative half, which §8.3's charge term needs — see the \
             comment above before widening this."
        );

        for seed in 0..12 {
            let (_, els, _physics) = table(seed);
            assert!(els.iter().all(|e| e.affinity > -1.0 && e.affinity <= 1.0));
            for w in els.windows(2) {
                let bound = 2.0 / w[0].units as f64;
                assert!(
                    (w[1].affinity - w[0].affinity).abs() <= bound,
                    "seed {seed} N={}: affinity jumped {}",
                    w[0].units,
                    (w[1].affinity - w[0].affinity).abs()
                );
            }
        }
    }

    /// The binding peak sits at a closure by *mechanism* — a closed shell has
    /// no frontier, so it maximises contacts per unit and the energy series is
    /// a sawtooth. If a change smooths the series across closures this fails,
    /// which is the point: the predecessor's minimum was the distance from a
    /// typed-in constant and looked identical from outside.
    ///
    /// **Enumerated, not sampled.** The doc above this test used to claim "zero
    /// failures across the full `(k, eps, sigma)` product" while the body ran
    /// `table(3)` — one seed. The claim is true (measured: 2268 closure cells
    /// over 972 universes, 0 non-maximal), and the product is 9 x 9 x 12 = 972
    /// cells costing milliseconds, so there is no reason to sample it. This file
    /// argues exactly that 140 lines up, about the valence ceiling: raising a
    /// seed count shrinks a probability where enumeration removes it.
    ///
    /// `eps` and `sigma` are reconstructed from their draw expressions in
    /// `generate_elements`; keep them in step.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "loop indices reconstructing the drawn constant grid"
    )]
    #[test]
    fn energy_per_unit_is_locally_maximal_at_every_closure() {
        let mut cells = 0_usize;
        for k in 6..=14_usize {
            let pack = packing::PackingConsts::new(k);
            for ei in 0..9_usize {
                for si in 0..12_usize {
                    let eps = 0.8 + 0.05 * ei as f64;
                    let sigma = 0.02 + 0.01 * si as f64;
                    let energy = |units: usize| {
                        eps * packing::contacts_upto(pack, units) / units as f64
                            - sigma * det_math::cbrt((units * units) as f64)
                    };
                    for c in packing::closures(k, 120).into_iter().filter(|&c| c > 1) {
                        if c + 1 > 120 {
                            continue;
                        }
                        cells += 1;
                        assert!(
                            energy(c) > energy(c - 1) && energy(c) > energy(c + 1),
                            "k={k} eps={eps} sigma={sigma}: closure {c} is not a local maximum"
                        );
                    }
                }
            }
        }
        assert!(cells > 2000, "the enumeration shrank to {cells} cells");
    }

    /// `peak` is an output. If it took only one value across universes it
    /// would merely have moved from an argument to a hard-coded location.
    #[test]
    fn the_binding_peak_moves_with_the_seed() {
        let peaks: std::collections::BTreeSet<usize> = (0..24).map(|s| table(s).0.peak).collect();
        assert!(
            peaks.len() >= 6,
            "peak took only {} values: {peaks:?}",
            peaks.len()
        );
    }

    /// Different universes must get different table *shapes*, not the same
    /// shape with different numbers in it.
    #[test]
    fn shell_patterns_vary_between_universes() {
        let patterns: std::collections::BTreeSet<Vec<usize>> =
            (0..30).map(|s| table(s).0.closures).collect();
        assert!(
            patterns.len() > 5,
            "shell patterns barely vary: {}",
            patterns.len()
        );
    }

    /// **Every numeric bound named in `generate_elements`'s `#[expect]` reason,
    /// computed rather than asserted in prose.**
    ///
    /// Three of those clauses were wrong, found by three separate lanes. A
    /// bound stated only in a reason string is the stale-figure class with a
    /// compiler-shaped alibi: it is what the next reader trusts *instead of*
    /// re-deriving, and nothing fails when a draw range moves. So the reason
    /// now cites this test and this test computes the numbers.
    ///
    /// The shell walk mirrors `generate_elements` exactly; keep them in step.
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "`frontier_notches` is non-negative and this test's own assertion bounds it"
    )]
    #[test]
    fn the_bounds_named_in_the_expect_reason_are_measured() {
        let (mut max_shell, mut max_cap, mut max_outer) = (0_usize, 0_usize, 0_usize);
        let mut valences = std::collections::BTreeSet::new();
        for k in 6..=14_usize {
            for units in 1..=120_usize {
                let (mut shell, mut filled) = (0_usize, 1_usize);
                while filled + packing::shell_size(k, shell + 1) <= units {
                    shell += 1;
                    filled += packing::shell_size(k, shell);
                }
                let cap = packing::shell_size(k, shell + 1);
                let outer = units - filled;
                max_shell = max_shell.max(shell);
                max_cap = max_cap.max(cap);
                max_outer = max_outer.max(outer);
                valences.insert(packing::frontier_notches(cap, outer).round_ties_even() as u8);
            }
        }
        assert_eq!(max_shell, 3, "`period` no longer fits the stated bound");
        assert_eq!(max_cap, 130, "`cap` moved — the reason string says 130");
        assert_eq!(max_outer, 73, "`group` moved — the reason string says 73");
        assert_eq!(
            valences,
            [0, 1, 2, 3, 4, 5, 6]
                .into_iter()
                .collect::<std::collections::BTreeSet<u8>>(),
            "the attained valence set moved; note this is the set of *values*, where \
             `valence_ceiling_is_four_to_six_over_the_drawn_range` asserts the set of \
             per-table *maxima* — conflating the two is what made the reason string wrong"
        );

        // The fourth shell *completes* at 189 (k = 6) and later for larger `k`,
        // which is why `period` stops at 3 — but its *capacity* is in use from
        // N = 92, which is why `cap` reaches 130 inside the drawn range. The
        // reason string said "the fourth shell opens past N = 189" and was
        // ambiguous between the two; a reader raising the element cap would
        // have taken 189 as the safe ceiling.
        let mut fourth_completes = usize::MAX;
        let mut fourth_starts_filling = usize::MAX;
        for k in 6..=14_usize {
            let (mut shell, mut filled) = (0_usize, 1_usize);
            for units in 1..=400_usize {
                while filled + packing::shell_size(k, shell + 1) <= units {
                    shell += 1;
                    filled += packing::shell_size(k, shell);
                }
                if shell == 3 && units > filled && fourth_starts_filling == usize::MAX {
                    // `units > filled` is `outer > 0`: at N = 91 shell 3 *closes*
                    // and the fourth shell holds nothing yet.
                    fourth_starts_filling = units;
                }
                if shell == 4 {
                    fourth_completes = fourth_completes.min(units);
                    break;
                }
            }
        }
        assert_eq!(
            fourth_completes, 189,
            "the fourth shell's completion point moved"
        );
        assert_eq!(
            fourth_starts_filling, 92,
            "the fourth shell's first fill moved"
        );

        let (mut min_raw, mut max_raw) = (i64::MAX, i64::MIN);
        for seed in 0..400 {
            for e in &table(seed).1 {
                min_raw = min_raw.min(e.mass.raw());
                max_raw = max_raw.max(e.mass.raw());
            }
        }
        assert!(
            (1024..=215_040).contains(&min_raw) && (1024..=215_040).contains(&max_raw),
            "mass sub-units {min_raw}..={max_raw} left the stated range"
        );
    }

    /// The naming margin, observable from a caller rather than only from a test
    /// that calls `mint` directly.
    ///
    /// `Provenance` was computed and dropped into `_` at the only production
    /// call site, so a `Fallback` symbol — which violates `Element::symbol`'s
    /// 1–2 character contract — would have been silent. The margin is 165
    /// symbols against a cap of 120; this is what fires if someone raises the
    /// cap without widening `ONSETS`.
    #[test]
    fn no_drawn_table_falls_back_to_a_minted_symbol() {
        for seed in 0..64 {
            let (shell, els, _physics) = table(seed);
            assert_eq!(
                shell.fallback_symbols,
                0,
                "seed {seed}: {} of {} symbols came from the exhaustion fallback",
                shell.fallback_symbols,
                els.len()
            );
        }
    }

    /// **A permutation guard, and the only test in this crate that can see a
    /// bit-level change to any output.**
    ///
    /// Every constant in `generate_elements` comes from one sequential stream,
    /// so transposing two `next_range` calls — or inserting one anywhere but the
    /// end — rewrites every universe while leaving each constant inside its own
    /// drawn range. Nothing else here could see that: `generation_is_deterministic`
    /// compares a table to itself and its own comment says so, and every other
    /// assertion is a range, a shape or a set. Measured before this existed:
    /// swapping the `eps` and `sigma` draws moved seed 0's last
    /// `energy_per_unit` from 4.295755047544661 to 4.481116661947341 with 254
    /// tests and `cargo xtask` green.
    ///
    /// It matters because the failure is between *versions of this repository*,
    /// not between machines: a seed someone shared before a refactor names a
    /// different periodic table afterwards, with nothing announcing it. That is
    /// the hazard `borbax_rng`'s `Domain` doc calls grim, and the doctrine was
    /// applied to the `Domain` enum and not to the draw sequence inside one
    /// stream. `borbax-cli` does not exist yet, so §13.6's matrix has nothing to
    /// compare and the state hash is Tasks 20–21 — roughly sixteen tasks during
    /// which any silent bit-mover would be green, after which the golden is
    /// minted from whatever the code says *then* and the drift is baked in.
    ///
    /// **`to_bits()`, not a formatted float**, so it also pins the accumulation
    /// order in `contacts_upto`: fusing its two `total +=` statements into one
    /// sum moves 167 of 1080 `(k, N)` cells by up to 1.14e-13, which no printed
    /// figure would show.
    ///
    /// When this moves, say which of §18.1's three it is — an intended physics
    /// change, a `.prototools` bump, or a bug. "Regenerated it" is the failure
    /// mode.
    #[expect(
        clippy::as_conversions,
        clippy::cast_sign_loss,
        reason = "hashing raw bit patterns and small counts; every value is pinned by \
                  the assertion itself"
    )]
    #[test]
    fn the_universe_digest_is_pinned() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |x: u64| {
            h ^= x;
            h = h.wrapping_mul(0x0100_0000_01b3);
        };
        for seed in 0..64 {
            // `_physics`, deliberately: this test fixes `PhysicsVersion::
            // CURRENT` (V1) via `table()` and covers V1 only — see
            // `the_v2_universe_digest_is_pinned` for V2's own coverage.
            //
            // **An earlier version of this comment claimed mixing element
            // content here would duplicate
            // `the_assembled_universe_digest_is_pinned`'s job — false, and a
            // `/review-pr` pass (2026-08-10, determinism-auditor) measured
            // it directly: reversing this test's own `for e in &els`
            // accumulation order moved 105 of 120 species' hashes with
            // **zero** test failures anywhere in the workspace, because
            // `the_assembled_universe_digest_is_pinned` mixes `seed`, the
            // version discriminant, `UniverseConsts`'s eight fields and bond
            // cells — never an `Element` field.** This test is the *only*
            // thing pinning element content, for whichever version it
            // covers, which is exactly why a V2-covering sibling is
            // necessary rather than redundant.
            let (sp, els, _physics) = table(seed);
            // **Destructured, for the reason the assembled-universe digest is.**
            // A hand-written field list cannot see a field that is not in it: a
            // review added `pub mutant_species_field: f64` to `Element`, filled
            // it at the single construction site, and all 55 tests plus
            // `clippy -D warnings` stayed green. `Element` is *the* type §8.6
            // designates for per-species precomputed quantities, so Tasks 9 and
            // 11 will add fields to it by design — this is the digest most
            // likely to be silently outrun.
            //
            // **E0027 forces a field to be *mentioned*, not *mixed*.** Binding
            // one and never hashing it is caught only by `unused_variables` in
            // the clippy leg, and `field: _` opts out silently through all six.
            // A `_` here is therefore a deliberate act and a review finding —
            // there is exactly one below, and it carries its reason.
            let ShellPattern {
                law,
                closures,
                fallback_symbols,
                peak,
                eps,
            } = &sp;
            // `table()` fixes `PhysicsVersion::CURRENT`, which is V1 — this
            // digest predates `PhysicsVersion` and has no version axis to mix
            // `law` itself into (see this test's own doc, above), so the
            // V1-only extraction below reproduces exactly what `mix(*k as
            // u64)` mixed before `k` moved into `ShellLaw`.
            let ShellLaw::V1 { k } = law else {
                unreachable!("table() fixes PhysicsVersion::CURRENT, which is V1")
            };
            mix(*k as u64);
            mix(*peak as u64);
            mix(*fallback_symbols as u64);
            // Added when `bonds.rs` began deriving its energy scale from `eps`,
            // which gave it a second path to a result. `sigma` reaches results
            // only through `energy_per_unit`, which is already mixed below, so
            // it is not carried here.
            mix(eps.get().to_bits());
            for c in closures {
                mix(*c as u64);
            }
            for e in &els {
                let Element {
                    symbol,
                    name,
                    // The one deliberate omission: `units` is the generator's
                    // loop index (`for units in 1..=n_elements`, pushed in
                    // order), so it is exactly the element's position in this
                    // very iteration and adds no information the digest does not
                    // already carry. Mixing it would move the constant to record
                    // a value that cannot independently change.
                    units: _,
                    period,
                    group,
                    mass,
                    valence,
                    affinity,
                    radius,
                    energy_per_unit,
                    instability,
                    abundance,
                } = e;
                mix(mass.raw() as u64);
                mix(u64::from(*valence));
                mix(u64::from(*period));
                mix(u64::from(*group));
                mix(affinity.to_bits());
                mix(radius.0.to_bits());
                // `.0.to_bits()`, matching the `radius` line above rather than
                // switching to `canonical_bits()`. The point of this line
                // during Task 5 is to prove that typing the field as `Quanta`
                // moved no value, and it can only prove that if it hashes the
                // same bits it hashed before.
                mix(energy_per_unit.0.to_bits());
                mix(instability.to_bits());
                mix(abundance.to_bits());
                for b in symbol.bytes().chain(name.bytes()) {
                    mix(u64::from(b));
                }
            }
        }
        assert_eq!(
            h, 0x4e79_4d66_eee9_ca47,
            "the universe digest moved — say which of §18.1's three this is, or \
             the fourth: the digest's own input set changed shape. Recomputing \
             with the changed field distinguishes them. This constant has moved \
             three times deliberately and no *value* this digest still covers has \
             ever changed: once when `eps` joined `ShellPattern` so `bonds.rs` \
             could derive its energy scale from the table's own; once when a \
             companion `sigma` was removed again for having no reader outside a \
             test; and once when Task 26.1 Step 8 removed `outer_fill_band` \
             entirely (a rebinning of `(period, group)` that never carried \
             chemistry — see this file's own header) — the input set narrowed by \
             one field, which moves the hash on its own even though every \
             remaining mixed value is unchanged"
        );
    }

    /// V2's own element-content digest — the coverage gap a `/review-pr`
    /// determinism pass found (2026-08-10): `PhysicsVersion::ALL` grew to
    /// `[V1, V2]` and pinned an assembled-universe digest for V2, but that
    /// digest mixes no `Element` field at all (see
    /// `the_universe_digest_is_pinned`'s own updated doc), so V2's `radius`,
    /// `instability`, `period`, `group`, `abundance`, `mass` and
    /// `energy_per_unit` had no pinned value anywhere — `radius` alone
    /// reaches `embed -> signature -> binding`, a result-affecting path with
    /// zero determinism coverage. Deliberately a separate test rather than
    /// a version-parameterised `the_universe_digest_is_pinned`: `ShellLaw`'s
    /// two variants carry different fields (`V1 { k }` against `V2`'s unit
    /// variant), so there is no single destructuring that covers both
    /// without reintroducing the exact silent-narrowing risk this file's own
    /// header warns about.
    #[expect(
        clippy::as_conversions,
        clippy::cast_sign_loss,
        reason = "hashing raw bit patterns and small counts; every value is pinned by \
                  the assertion itself — same reason as the_universe_digest_is_pinned"
    )]
    #[test]
    fn the_v2_universe_digest_is_pinned() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |x: u64| {
            h ^= x;
            h = h.wrapping_mul(0x0100_0000_01b3);
        };
        for seed in 0..64 {
            let (sp, els, _physics) = table_v2(seed);
            let ShellPattern {
                law,
                closures,
                fallback_symbols,
                peak,
                eps,
            } = &sp;
            let ShellLaw::V2 = law else {
                unreachable!("table_v2() fixes PhysicsVersion::V2")
            };
            mix(*peak as u64);
            mix(*fallback_symbols as u64);
            mix(eps.get().to_bits());
            for c in closures {
                mix(*c as u64);
            }
            for e in &els {
                let Element {
                    symbol,
                    name,
                    units: _,
                    period,
                    group,
                    mass,
                    valence,
                    affinity,
                    radius,
                    energy_per_unit,
                    instability,
                    abundance,
                } = e;
                mix(mass.raw() as u64);
                mix(u64::from(*valence));
                mix(u64::from(*period));
                mix(u64::from(*group));
                mix(affinity.to_bits());
                mix(radius.0.to_bits());
                mix(energy_per_unit.0.to_bits());
                mix(instability.to_bits());
                mix(abundance.to_bits());
                for b in symbol.bytes().chain(name.bytes()) {
                    mix(u64::from(b));
                }
            }
        }
        assert_eq!(
            h, 0x74f0_d6b0_999c_25ab,
            "the V2 universe digest moved — say which of §18.1's three this is, or the \
             fourth: the digest's own input set changed shape. This golden was first \
             measured 2026-08-10, alongside the assembled-universe V2 golden — see this \
             test's own doc for what it closes"
        );
    }

    #[test]
    #[should_panic(expected = "cannot address its own slots")]
    fn a_table_too_long_to_address_is_refused_at_construction() {
        // 300 slots: measured before the guard, `iter` handed out 256 distinct
        // ids, 44 slots disagreed with `get`, and the bond matrix priced slot
        // 299 as slot 255 — in range, plausible and silent.
        let (pattern, els, physics) = generate_elements(0, PhysicsVersion::CURRENT).into_parts();
        let mut long = els.clone();
        while long.len() <= PeriodicTable::MAX_ELEMENTS {
            long.extend(els.iter().cloned());
        }
        let _ = PeriodicTable::new(pattern, long, physics);
    }

    #[test]
    fn table_size_is_workable() {
        for seed in 0..30 {
            let n = table(seed).1.len();
            // The draw is 60..=120, so the lower half of a `40..=120`
            // assertion was a guard that could not fire — the defect class this
            // task is partly about, asserted 200 lines from a comment rejecting
            // it.
            assert!((60..=120).contains(&n), "seed {seed} produced {n} elements");
        }
    }
}
