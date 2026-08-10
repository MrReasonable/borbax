//! Generated bond energies (spec §7.1).
//!
//! **Bond energy is derived from the packing — both its ordering and its
//! scale.** §7.1's `bond_energies` was the last substantive property of a
//! universe still without a cause: Task 4 derived `period`, `group`, `valence`,
//! `radius`, `affinity` and `mass` from how a cluster packs, and left this one
//! to be drawn. Drawing it reproduces exactly the defect that made the
//! predecessor's `peak` unfalsifiable — "these elements bind well" encoded
//! twice, in incommensurable units, so neither encoding can check the other.
//!
//! A bond is a *contact*, of the same kind as the ones holding a cluster
//! together. Two things follow, and the file now implements both rather than
//! only the first:
//!
//! - **What a bond is made of** is contact density, `contacts_upto(N)/N` — not
//!   `energy_per_unit`, which subtracts the strain of packing a cluster into a
//!   sphere. That strain is about a cluster's own geometry and has nothing to
//!   do with bonding to another one.
//! - **The energy scale of a bond** is `eps`, the drawn per-contact energy, so
//!   `base = eps * scale` with `scale` dimensionless.
//!
//! Only two constants are drawn here and both are dimensionless: that multiple,
//! and the order exponent. That is principle 1 — one mechanism reused — rather
//! than a second axis invented for bonding.
//!
//! # What this cost, and what a review had to find first
//!
//! An earlier version normalised `energy_per_unit` onto a floored `[w, 1]`
//! scale, with `base` drawn independently in Quanta. That map is
//! **affine-invariant**: multiply every `energy_per_unit` by 1000 and the worst
//! relative change in any bond energy was **3.7e-16**, one ulp. So the matrix
//! carried the table's *ordering* and discarded its *magnitude* — and magnitude
//! is what enters `exp(-E/T)`, hence what §9.4's differential persistence is
//! made of. The file's own defect, one axis over, with `base` as the second
//! incommensurable encoding.
//!
//! `scaling_the_tables_energy_scales_every_bond_energy` is the discriminator,
//! and it is the one a plausible fix fails: widening `base`'s interval leaves
//! the ratio at 1.0.
//!
//! Retiring the normalisation also retired `weakest_share`, whose floor existed
//! only because a min-max map puts a zero at the table minimum and zero is
//! absorbing under a geometric mean. Contact density is positive for every
//! element that can bond, so positivity is structural. Measured over 500
//! universes, the bondable single-bond ratio is now **[6.599, 10.826]** against
//! the drawn scheme's [2.301, 5.279] — wider, the direction §22.6 wants — and
//! it varies with `k`, which `1/weakest_share` could not.
//!
//! The monomer scores 0 and prices 0. That is correct by meaning: no frontier,
//! `valence == 0`, no bonds. It collides with the unknown-id sentinel, which is
//! one more reason Task 14 must retire that.
//!
//! # Two things the first draft of this file got wrong, both measured
//!
//! The version in the plan keys the matrix on `group` and derives it from
//! *affinity opposition*, `-(affinity_a * affinity_b)`, under the comment
//! "complementary affinities bond strongly; like affinities bond weakly".
//! Neither half survives contact with the table Task 4 actually generates.
//!
//! **Affinity has no negative half, so there is no opposition to read.**
//! Measured over 36 013 elements across 400 universes, `affinity` spans
//! `[0.14047619047619042, 1.0]` with **zero** negative values — it is
//! `2·(surface/units) − 1` and a cluster's exposed fraction does not fall below
//! a half. So `-(a·b)` is negative for *every* pair, the expression is monotone
//! decreasing in both arguments, and what it actually says is "high-affinity
//! elements bond weakly to everything" — not complementarity. This is the same
//! open finding Task 10 carries for §8.3's charge term, which is where the fix
//! belongs, and it is why nothing here reads `affinity` at all.
//!
//! **`group` is not a family here, and the measurement is stronger than the
//! aliasing argument.** `group` is `outer`, the occupancy of the incomplete
//! shell, so a value recurs only once per period — and `period <= 3`. At most 4
//! elements share a group. But the decisive figure is not the aliasing: against
//! a permutation null with identical class sizes, `group` explains bond
//! character with an R² *ratio of 0.296–0.893* — **below chance**, so
//! partitioning by group is worse than an arbitrary partition of the same
//! shape. `period` scores 16.3–53.4× above it. Both are sample ranges over 40
//! seeds and are quoted as such; the qualitative claim (period far above
//! chance, group below it) is what the decision rests on, and nothing in this
//! crate guards the R² itself — `the_table_still_justifies_ignoring_affinity_and_group`
//! guards the affinity floor and the multiplicity bound, which this paragraph
//! demotes.
//!
//! So families do share bond character here, and the axis is the **period, not
//! the group** — an earlier version of this header implied the column. §7.1's
//! "everything in this column makes rings" payoff is still available, through
//! `valence`, which is a function of `(cap, outer)` by construction.
//!
//! # The matrix is exactly rank 1, and that is a design position
//!
//! `bond = base·√(c_a·c_b)` factorises as `x_a·x_b`, so every 2×2 minor
//! vanishes (measured worst relative minor 6.2e-16) and `E(a,c)/E(b,c)` is
//! independent of `c`. **Every element ranks its partners in the same order**;
//! an `n`-element table stores `n²` cells holding `n+1` independent numbers.
//! There is no selective covalent bond in this chemistry and no exchange
//! reaction whose sign depends on the spectator.
//!
//! **That is a deliberate departure from §7.1, which specifies `bond_energies`
//! as a "per-partner-group bond strength matrix" — i.e. the spec asks for pair
//! structure and this does not supply it.** Recording it as a departure rather
//! than as principle 2 enforcing itself, because the stronger claim does not
//! hold: §3.2 enumerates *recognition* phenomena, and covalent dissociation
//! energy is §9.2's axis. The narrower argument does hold, and is the reason:
//! a pair-structured covalent energy would give partner preference a **second
//! source** alongside §8.3's complementarity — two independent axes scoring
//! "does A prefer B" — and §3.2's "if a feature needs a second mechanism,
//! question the feature" applies to that duplication. It is stated here because the type is
//! called a matrix, and a reader building §9.1's energetics on top will
//! otherwise assume pair structure exists. If selectivity is ever wanted here,
//! the question to answer first is what the shape model failed to do.
//!
//! # The `E/T` regime, restated after the derivation changed it
//!
//! Four review lanes converged on a finding that `E/T` sat in `exp`'s linear
//! regime, so the *rate* spread was near 1 while the energy spread looked fine.
//! **Deriving the scale substantially changed that, and the old figures are kept
//! nowhere — they were measured under a map this file no longer uses.**
//!
//! **The figures are asserted, not quoted.** Three consecutive rounds wrote a
//! measured table into this header and the next commit invalidated it — the
//! stream-index move alone shifted seed 0's `temp_min` through 206.09, 222.68
//! and 197.24, and every row here divides by `T`. Prose cannot survive that, so
//! `the_header_figures_are_current` measures each row and fails when one drifts.
//! Read the numbers there; this paragraph states only the shape.
//!
//! `E/T` rose, **and the claim needs its qualifier or it is misleading.** Across
//! the orders a pair can actually *form* — capped at `min(valence_a, valence_b)`
//! — the cleave-rate ratio spans several orders of magnitude and exceeds 2.0 in
//! every universe drawn, where before the derivation change it was near 1 in all
//! of them.
//!
//! **But that is carried by bond order, not by the energy scale.** For *single*
//! bonds alone, max `E/T` is below 1 in **218 of 500** universes (41.5% at 200
//! seeds, 43.6% at 500, 45.0% at 2000 — it is a sample statistic, not a
//! constant). Order 1 is not "the only order many pairs can form" — only 11.0%
//! of bondable pairs are capped there. The true and stronger statement is that
//! order 1 is the only order **every** pair can form, and the floor for all of
//! them, so this is a claim about the whole chemistry's baseline.
//!
//! **Do not read that as "the old diagnosis still holds in 44%".** It does not,
//! and the difference is a Task 14 sizing input. Conditioning on exactly the
//! flagged universes, at `T_min` over 2000 seeds, the single-bond rate ratio
//! there is **[1.3616, 2.4430], mean 1.9657** — against the **1.04–1.24x** that
//! prompted the original finding. The worst flagged universe now beats the best
//! pre-change one. What is true is the weaker claim: those universes sit where
//! `exp` is nearly linear, so composition alone cannot lift the ratio past `e`
//! there. (Elsewhere: [2.3625, 12.7144], mean 3.6620.)
//!
//! An earlier version of this paragraph said "no longer near zero" without the
//! order split at all, and priced every pair at order 6 to say it: measured,
//! **78.1% of universes cannot form an order-6 bond**, so the figure came from
//! an energy no bond in four universes out of five can carry.
//!
//! Task 14 must size against the formable numbers. `the_header_figures_are_current`
//! now caps at the per-pair valence bound and measures `[11.40, 2370.55]` energy,
//! `[3.58, 3.56e4]` rate ratio, max `E/T` 10.68 over 500 seeds.
//!
//! **The absolute energy scale rose about fivefold, and that was not stated
//! when it happened.** Capacity was normalised onto `[w, 1]` and is now contact
//! density, which runs toward `z/2`. On seed 0, over **formable** orders across
//! bondable pairs, energies now span `38.98..1729.49`. The new single-bond rate
//! *minimum* is roughly the old *maximum*, and orders reaching 6 multiply the
//! exponent again.
//!
//! **A third consequence rode along and is worth naming separately**: the draw
//! count on `Domain::Universe` index 1 went from four to two, so every constant
//! after the bond draws shifted — and then the bond generator moved to its own
//! index 2, shifting them again. Seed 0's `temp_min` is now **197.2437**; this
//! paragraph has quoted 206.09 and then 222.68, each true of a tree two commits
//! back. That is not "the bond derivation changed"; it is every shared seed
//! naming a different universe, which is exactly what a digest exists to make
//! visible rather than silent.
//!
//! **What that does and does not settle.** §22.6 wants "some linkages durable
//! and others fragile". From composition alone — order held at `SINGLE`, so the
//! only variation is which two elements are joined — the cleave-rate ratio over
//! 500 seeds is **[1.1886, 4.2678] at `T_mid`** and **[1.4182, 12.7144] at
//! `T_min`**. Both are a different situation from the 1.04–1.24x that prompted
//! the finding. Whether either *suffices* is a chemistry judgement for Task 20's
//! battery, not a claim this file should make.
//!
//! **Name the temperature or the number means nothing.** An earlier version of
//! this paragraph said "1.19–3.73x", which took its low end from `T_mid` and its
//! high end from no sample that reproduces, then argued in the next clause from
//! a `T_min` figure — one sentence, two temperatures, neither labelled. Every
//! rate figure in this header now carries `(T, sample)` or is not here.
//!
//! What it does settle is that any sizing done against the old numbers is wrong:
//! an implementer applying Task 14's coupling on top of these while reading
//! `1.043–1.243` would oversize it badly.
//!
//! **A STOP stood here and it was a false alarm — recorded because a Task 14
//! implementer opens this file first.** It said drawing the energy in units of
//! the universe's own temperature "would collapse a cleave-rate spread of 3 to
//! 217 758 down to 1.09x, deleting §9.4's differential persistence". Both of
//! those endpoints are stale too — measured at 2000 seeds the spread is
//! **2.4018 to 1.6094e6** — but the figures are the smaller problem. That does
//! not follow. `E = base·sqrt(c_a c_b)·n^gamma` is separable, so pinning
//! `base/T` across universes does not collapse the *within*-universe rate ratio
//! — it makes every universe have the *same* ratio, at whatever level the
//! nominal sets. The remedy was being blamed for a harm it does not do.
//!
//! The real caution is the opposite one and it is Task 20's: `log R` is exactly
//! linear in `base/T`, so a tolerance stated on the dimensionless group is a
//! tolerance on the *logarithm* of the rate ratio — "±6% on `Ea/T`" is a 3.5x
//! cross-universe spread in `R` at the measured high end, not a pin. Task 14
//! should fix the group; Task 20 chooses the value, against an observable it can
//! print rather than against the group.
//!
//! The remaining question is a *band*, not a direction: it belongs with the rate
//! law, is written into Task 14 with its discriminator, and into Task 20 for the
//! band. Picking a number here would be the typed-in constant this
//! project keeps finding reported back as an emergent result — and note that
//! deriving the *energy* scale, as this file now does, does not fix the
//! *ratio*: `eps` is still drawn independently of `temp_min`. The remedy Task 14
//! carries is to draw the energy in units of the universe's own temperature,
//! which a review measured collapses the `Ea/T` spread from 4.39x to 1.09x.

use crate::PhysicsVersion;
use crate::element::{ElementId, PeriodicTable};
use crate::perturbation::{MigratedConstant, Rung, draw_symmetric};
use borbax_rng::Stream;
use borbax_units::{Quanta, det_math};

/// How many contacts join one pair of atoms.
///
/// **A validated newtype over `1..=MAX`, not a three-variant enum, and both
/// halves of that changed for measured reasons.**
///
/// It was `Single`/`Double`/`Triple` — the real-chemistry taxonomy verbatim,
/// which is a G3/G6 concern on its own — capped at an undrawn `3`. Measured
/// over 200 universes, the per-universe ceiling of `min(valence_a, valence_b)`
/// across bondable pairs takes the values `{4, 5, 6}` and is **never 3**:
/// 25.6% of bondable pairs have both partners able to support order 4 or more.
/// So the one bonding rule in this crate that was not generated sat below the
/// generated table's own ceiling in 200 of 200 universes, while landing exactly
/// on real chemistry's main-group maximum.
///
/// The cap is now [`Self::MAX`] = 6, the valence ceiling `element.rs` already
/// asserts as the attained set. A newtype keeps illegal states unrepresentable
/// — the reason the enum was right — without inventing six names for what is a
/// count.
///
/// **Enforcing a valence budget is not this type's job**, and the routing has
/// been wrong twice. It is not Task 11's (`Polymer::push`) and it is not Task
/// 6's `add_bond`, which has no `&Universe` and stores element *indices*, not
/// valences. The per-pair bound `min(valence_a, valence_b)` — a bond of order
/// *n* consumes *n* slots at each end — is also only **necessary, not
/// sufficient**: at 12 atoms a vertex has up to 11 neighbours, so a per-pair
/// check alone admits 11x the budget at every valence. The rule that holds is
/// per *atom*: `valence_used(a) <= valence(elem[a])`. Its home is
/// `Interner::intern`, which already takes `&Universe`, is already per-species
/// (§8.6), and can return `Result<SpeciesId, OverValence>` so that no
/// `SpeciesId` reachable from a simulation step can name an over-valent
/// molecule. Nothing here forms a bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondOrder(u8);

impl BondOrder {
    /// A single contact.
    pub const SINGLE: Self = Self(1);
    /// Every representable order, ascending.
    ///
    /// **One array is the single source of the bound, its `u8` form, and the
    /// multiplier table's length.** It replaces a `MAX: u8` / `MAX_LEN: usize`
    /// pair that were two independent literals: two reviews measured that
    /// raising one without the other made order 7 price at `1.0` — *exactly*
    /// order 1's multiplier, a 4.8x drop — with the only failure being the
    /// digest, whose message then blames the order range for widening. And
    /// `element.rs` instructs a future reader to raise the cap.
    ///
    /// It also lets a consumer iterate or clamp without a fallible conversion.
    /// `MAX` as a bare `u8` beside `SINGLE: Self` forced `BondOrder::new(MAX)`
    /// plus a `None` arm nobody can `unwrap` in a workspace that denies it —
    /// visible inside this crate as an unreachable `if let` in the digest and a
    /// divisibility guard built to compensate for it.
    pub const ALL: [Self; 6] = [Self(1), Self(2), Self(3), Self(4), Self(5), Self(6)];
    /// The largest order any universe's valence ceiling admits.
    /// Derived from `ALL`'s **last element**, not its length — no `as`, no
    /// suppression, and it reads the array's content rather than its size.
    pub const MAX: u8 = match Self::ALL.last() {
        Some(o) => o.0,
        // **Not `0`.** `MAX = 0` makes `new`'s range `1..=0` empty, so every
        // `BondOrder::new` in the workspace returns `None` — silently, at every
        // call site — which is precisely the plausible-value fallback this file
        // argues against three times elsewhere. A panic in a const initializer
        // is `E0080` at compile time, not a runtime abort.
        None => panic!("BondOrder::ALL is empty"),
    };

    /// How many per-order planes a consumer must allocate to hold every order.
    ///
    /// **Naming, not a new guarantee — the doc here previously overstated it.**
    /// A consumer storing one plane per order indexes it as `order - 1`, so
    /// `[T; PLANES]` says what it means, where `[T; ALL.len()]` says an adjacent
    /// thing that coincides. It was claimed that a sparse `ALL` (`[1, 2, 3, 6]`)
    /// would make `ALL.len()` 4 against a `MAX` of 6 and put order 6 out of
    /// bounds. A review checked: **that does not compile.** The `const _` below
    /// rejects it with `E0080`, in this file, so introducing a sparse `ALL` is a
    /// deliberate act of deleting an assertion rather than an accident waiting
    /// to happen.
    ///
    /// What remains true and is the reason to keep this: `borbax-molecule` can
    /// neither see nor cite that `const _`, so `PLANES` is the name it can
    /// depend on without depending on an invariant it cannot read. Task 6's
    /// draft wrote `ALL.len()`, and so — until a review caught it — did
    /// `OrderScale::mult`, the one consumer of this pattern already in the
    /// crate.
    #[expect(
        clippy::as_conversions,
        reason = "`usize::from` is not stable as a const fn on 1.97.1 (`From` is \
                  not yet a const trait), and this must be a `const` because \
                  consumers use it as an array length. Widening u8 -> usize is \
                  lossless on every target this workspace supports."
    )]
    pub const PLANES: usize = Self::MAX as usize;

    /// This order's plane, for a consumer holding one array per order.
    ///
    /// **The arithmetic lives here so no call site does it.** Task 6's draft
    /// wrote `let plane = (order - 1) as usize` against an unvalidated `u8`,
    /// which panics at `order == 0` (subtraction overflow in debug; wraps to 255
    /// and indexes out of bounds in release, where `overflow-checks` is off) and
    /// at `order >= 7`. Both profiles verified. Taking `self` means the value
    /// was already validated by [`Self::new`] or came from [`Self::ALL`], so the
    /// result is in `0..PLANES` by construction and the call site has no `- 1`,
    /// no cast, and no `#[expect]`.
    #[must_use]
    pub fn plane_index(self) -> usize {
        // Not `const`, so this can use `usize::from` rather than `as`. `- 1`
        // cannot underflow: `self.0 >= 1` for every constructible value.
        usize::from(self.0) - 1
    }

    /// An order, or `None` outside `1..=MAX`.
    #[must_use]
    pub const fn new(n: u8) -> Option<Self> {
        if n >= 1 && n <= Self::MAX {
            Some(Self(n))
        } else {
            None
        }
    }
}

/// **`ALL` must be `1..=MAX` ascending, checked by the compiler.**
///
/// Nothing else checks its *contents*. `OrderScale` fills and reads `mult` by
/// index, so a permuted or duplicated `ALL` never reaches the multiplier table:
/// a review measured `[1, 2, 3, 4, 6, 5]` failing **only** the pinned digest,
/// 57 of 58 green, with the array contradicting its own doc line. And a digest
/// is re-pinned by hand on every deliberate change, so it is not a guard against
/// this — the argument this file makes about itself twice already.
#[expect(
    clippy::as_conversions,
    clippy::indexing_slicing,
    reason = "const evaluation: an out-of-range index is E0080 at compile time, \
              not a runtime panic, so `indexing_slicing`'s hazard cannot arise; \
              and `u8 -> usize` is widening. `.get()` and `usize::from` are both \
              unavailable in a const `while` on the pinned toolchain"
)]
const _: () = {
    let mut i = 0;
    while i < BondOrder::ALL.len() {
        assert!(
            BondOrder::ALL[i].0 as usize == i + 1,
            "BondOrder::ALL must be 1..=MAX ascending"
        );
        i += 1;
    }
};

/// The error [`BondOrder::try_from`] returns for a count it cannot represent.
///
/// `thiserror` rather than a hand-rolled `Display`/`Error` pair: it is already
/// a workspace dependency and already compiled into this crate's tree via
/// `borbax-units`, so it costs zero extra crates and zero build time.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("bond order {0} is not in 1..={max}", max = BondOrder::MAX)]
pub struct NotABondOrder(pub u8);

impl TryFrom<u8> for BondOrder {
    type Error = NotABondOrder;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        Self::new(n).ok_or(NotABondOrder(n))
    }
}

/// Closes the round-trip, so the count is never recovered by an `as`.
impl From<BondOrder> for u8 {
    fn from(order: BondOrder) -> Self {
        order.0
    }
}

/// Contacts per unit — **what a bond is made of**.
///
/// The positive half of `energy_per_unit`, divided by `eps`. The strain term
/// that `energy_per_unit` subtracts is the cost of packing a cluster into a
/// sphere, which is about that cluster's own geometry rather than about bonding
/// to another one, so it has no place in a bond.
///
/// Zero exactly at `units == 1`, where `contacts_upto` is 0 — the monomer, which
/// has no frontier and cannot bond. `u16::try_from` rather than `as`:
/// `PeriodicTable::new` bounds tables at 256 elements, so the conversion is
/// lossless with eight bits to spare, but `as_conversions` is denied
/// workspace-wide because a silent narrowing is the gotcha most likely to
/// produce a plausible wrong number.
///
/// **The unreachable fallback is `NAN`, not `INFINITY`, and the difference is a
/// defect a review caught.** `x / INFINITY` is `0.0`, and under §9.4's
/// `k = A·exp(-E_bond/T)` a zero bond energy is the **fastest-cleaving** value —
/// the exact error [`BondEnergyMatrix::energy`]'s own doc spends thirty lines
/// recording, reintroduced here by a division rather than a literal, and hidden
/// one step further for it. There is no conservative choice available: a
/// cluster above 65 535 units has a true density tending to `z/2`, the *other* end of the
/// range. So the branch holds poison instead of a plausible number.
fn contact_density(pack: crate::packing::PackingConsts, units: usize) -> f64 {
    let n = u16::try_from(units).map_or(f64::NAN, f64::from);
    crate::packing::contacts_upto(pack, units) / n
}

/// Every pairwise bond energy in one universe, computed once at generation.
///
/// Sized to the table it was generated from and indexed by [`ElementId`], so a
/// rate path can price a bond holding two ids and nothing else — no `&Element`,
/// no element table in scope, and no `String` field kept alive by a borrow.
#[derive(Debug, Clone, PartialEq)]
pub struct BondEnergyMatrix {
    /// Elements in the table this was built for.
    n: usize,
    /// Row-major single-bond energies, `n * n`, symmetric by construction.
    e: Vec<Quanta>,
    /// The order exponent. See [`OrderScale`].
    order_scale: OrderScale,
}

/// How much stronger an order-*n* contact set is than a single one.
///
/// **One drawn exponent, `n^gamma`, replacing two independently drawn
/// multipliers.** The pair had three problems the exponent removes at once.
///
/// Monotonicity was enforced by making the draw ranges disjoint (`[1.6, 2.1]`
/// and `[2.2, 3.0]`), which is a chosen constraint wearing the look of a
/// structural guarantee. Under `n^gamma` it is a theorem for every `gamma > 0`.
///
/// The *increments* were undecided and nobody had noticed: measured over
/// 2 000 000 draws, **40.04% of universes got an accelerating series** — the
/// third contact adding more energy than the second — because the two
/// multipliers were drawn independently. `gamma <= 1` makes diminishing returns
/// a single visible choice instead of a coin flip, and it has a mechanism
/// already in this codebase: `element.rs` charges strain for stacking contacts
/// into the same region.
///
/// And orders 4..=6 come free, which [`BondOrder::MAX`] needs. `gamma` is drawn
/// in `[0.68, 1.0]`, which reproduces the old ranges almost exactly (double
/// `[1.60, 2.00]` against `[1.6, 2.1]`; triple `[2.11, 3.00]` against
/// `[2.2, 3.0]`) with one constant instead of two.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderScale {
    // **`PLANES`, not `ALL.len()`.** This array is indexed at `order - 1` by
    // `of()`, so what it needs is `MAX`. They coincide only while `ALL` is
    // dense. A review pointed out that this — the one existing consumer of the
    // one-slot-per-order pattern — still said `ALL.len()` in the commit that
    // introduced `PLANES` to say otherwise; measured under a sparse `ALL` it
    // prices the top order as `NaN` and fails five tests including the digest.
    mult: [f64; BondOrder::PLANES],
}

impl OrderScale {
    /// Precomputed at generation, so `powf` never runs on a lookup.
    ///
    /// The old two-multiplier form was a field read; moving to an exponent put a
    /// transcendental in `energy()`, which is the one function a rate path will
    /// call per species. `gamma` is fixed per universe, so the table costs
    /// `ALL.len()` `powf` calls once.
    ///
    /// **Measured 32.3 ns/lookup inline against 1.0-1.3 precomputed — a 25-32x
    /// win, not the 4.6x an earlier version of this comment claimed.** That
    /// figure quoted 7.3 ns for the precomputed form, which was per-call
    /// harness overhead rather than the lookup: this repo's recorded
    /// microbenchmark incident, one level up. The decision was right and the
    /// baseline was wrong, and a per-lookup cost is exactly the number someone
    /// sizes a budget against later.
    fn new(gamma: f64) -> Self {
        let mut mult = [1.0; BondOrder::PLANES];
        for (i, m) in mult.iter_mut().enumerate() {
            // `i` runs 0..ALL.len(), so the narrowing is lossless; `as` is
            // denied. **`NAN` on the impossible branch, not `0`.** `0` there
            // saturates to order 1 and prices at `powf(1, gamma) == 1.0` —
            // *exactly* order 1's multiplier, which is the precise
            // silent-wrong-number `BondOrder::ALL`'s own doc says this type
            // exists to prevent. This file states the rule twice already, at
            // `contact_density` and at `of`: an unreachable fallback is `NAN`,
            // so it cannot be mistaken for a real answer.
            let order = u8::try_from(i).map_or(f64::NAN, |n| f64::from(n.saturating_add(1)));
            *m = det_math::powf(order, gamma);
        }
        Self { mult }
    }

    /// Total by construction, for every representable order.
    ///
    /// **`&self`, not `self`.** `OrderScale` is 48 bytes and `Copy`, and taking
    /// it by value with a dynamic index defeats SROA — LLVM materialises the
    /// copy in a stack slot and indexes *that*, surviving thin LTO at opt-level
    /// 3. Measured: the inner loop goes 7 instructions to 3, six memory ops to
    /// two, and the dev-profile frame drops from 64 bytes to 16. **-22.3%**,
    /// bit-identical over 3 087 174 cells, against a zero noise floor.
    ///
    /// The path is per-species at intern time so this buys nothing at the real
    /// call rate. It is here because the fix is one character with no downside,
    /// and because the shape — a `Copy` struct with an array field, taken by
    /// value, indexed dynamically — is what gets copy-pasted into the signature
    /// code where the rate is 10^6/s.
    ///
    /// `wrapping_sub` removes an unreachable underflow panic branch: the field
    /// is private and `new` bounds it at >= 1, but LLVM cannot see that, so a
    /// bare `- 1` compiles to a test-and-panic in dev.
    fn of(&self, order: BondOrder) -> f64 {
        let i = usize::from(u8::from(order)).wrapping_sub(1);
        // `NAN`, not `1.0`: `1.0` is exactly order 1's multiplier, so an
        // out-of-range order would price as a single bond and look plausible.
        // Same ruling `contact_density` makes about `INFINITY`; the file now
        // has one rule for unreachable fallbacks instead of three.
        self.mult.get(i).copied().unwrap_or(f64::NAN)
    }
}

impl BondEnergyMatrix {
    /// The multiple of [`crate::element::ShellPattern::eps`] that sets this bond
    /// energy scale, as a `(lo, hi)` draw range.
    ///
    /// **One definition, deliberately, because two identical ones do not pin a
    /// retune.** `both_bond_constants_are_this_universes_own_draws` replays the
    /// production stream to prove the constants are the seed's own. When that
    /// replay re-typed the literals, a single find-and-replace across the file —
    /// or a multi-cursor edit, which is how a range actually gets retuned —
    /// updated both copies and the test re-derived the *new* constant: measured,
    /// `s/84.0/84.5/` left 58 of 59 green, the digest alone firing. Sharing the
    /// range leaves the replay pinning what it can genuinely pin, the **order
    /// and count** of the draws, and leaves "the range itself moved" to the
    /// digest, which is the guard that is actually re-pinned by hand.
    pub(crate) const SCALE_RANGE: (f64, f64) = (28.0, 84.0);
    /// The bond-order exponent's draw range. See [`Self::SCALE_RANGE`] for why
    /// this is one definition rather than two.
    pub(crate) const GAMMA_RANGE: (f64, f64) = (0.68, 1.0);

    /// V2's bond-capacity floor (issue #26, Task 26.1 Steps 6-9, Decision 8)
    /// — `capacity(el) = BASE_CAPACITY_SCALE + valence(el)`, drawn through
    /// P7 as [`MigratedConstant::CapacityScale`]. Chosen, not derived, on
    /// the same order of magnitude as the valence values it sits alongside
    /// (typically 0-4 for a main-group element) — large enough that a
    /// zero-valence (noble-gas-analogue) element still prices a small
    /// positive capacity rather than exactly zero, small enough that
    /// valence differences still dominate the signal.
    pub(crate) const BASE_CAPACITY_SCALE: f64 = 0.5;

    /// V2's rank-2 ionic-excess scale (issue #26, Task 26.1 Step 6/8) — see
    /// [`MigratedConstant::IonicScale`]'s own doc for the full derivation,
    /// the worst-case positivity proof, and why this is `2.0` rather than
    /// the naive `4.0`.
    pub(crate) const BASE_IONIC_SCALE: f64 = 2.0;

    /// Derive the matrix from a table's binding energies.
    ///
    /// **Two constants are drawn and both are dimensionless** — the multiple of
    /// `eps` that sets the energy scale, and the order exponent. Everything
    /// dimensionful comes from the table. (An earlier version of this sentence
    /// said "three", naming a `weakest_share` this commit deleted and two order
    /// multipliers that became one exponent — three false clauses on the
    /// constructor's own doc, in a file whose subject is which constants get
    /// drawn.)
    ///
    /// What is deliberately **not** drawn is anything per-pair: the plan's draft
    /// added `next_normal() * 0.08` to every cell, which is an independently
    /// drawn matrix wearing a derivation, in the exact proportion of the noise.
    /// Every difference between two cells here traces to the two elements'
    /// **contact densities** — not their binding energies, which is a different
    /// quantity and the whole point of the change.
    ///
    /// **`pub(crate)` for the same reason [`PeriodicTable::new`] is.** A public
    /// constructor lets a caller build a matrix from table A and query it with
    /// ids from table B — in range, plausible and silent, which is the failure
    /// this crate keeps closing. Universes come from [`crate::Universe`], and
    /// the matrix comes with one.
    ///
    /// **Dispatches on `table.physics()`, not a second parameter.** P3
    /// (issue #26's prerequisites) added this: before it, `generate` took a
    /// table and nothing else, so there was no version to read and no way to
    /// stop a future variant's derivation from being bolted onto V1's body
    /// unnoticed. The match below is exhaustive over [`PhysicsVersion`] —
    /// despite `#[non_exhaustive]` not binding in-crate — so a second variant
    /// with no arm here is `E0004` at compile time, the same forcing function
    /// [`crate::Universe::generate_under`] already uses.
    ///
    /// **`seed` and `rung` are only consumed by the V2 arm** — V1's own
    /// energy-scale and order-exponent draws stay on `rng` unchanged,
    /// exactly as before this pair was added. Passed unconditionally rather
    /// than threaded in only for V2 so this stays one signature, not two.
    pub(crate) fn generate(table: &PeriodicTable, seed: u64, rung: Rung, rng: &mut Stream) -> Self {
        match table.physics() {
            PhysicsVersion::V1 => Self::generate_v1(table, rng),
            PhysicsVersion::V2 => Self::generate_v2(table, seed, rung, rng),
        }
    }

    fn generate_v1(table: &PeriodicTable, rng: &mut Stream) -> Self {
        let n = table.len();
        let pattern = table.pattern();

        // **Two draws, both dimensionless.** Every dimensionful constant here is
        // now the table's own, which is what makes the matrix carry the table's
        // energy *scale* rather than only its ordering.
        //
        // `scale` multiplies `eps`, the drawn per-contact energy. The header's
        // mechanism said this already: if a bond is the same kind of contact as
        // the ones holding a cluster together, its energy scale *is* the contact
        // energy. Drawing an independent `base` in Quanta encoded that quantity
        // twice, and a review measured the consequence — the matrix was exactly
        // invariant to `energy_per_unit`'s magnitude, 3.7e-16 under a x1000.
        let scale = rng.next_f64_range(Self::SCALE_RANGE.0, Self::SCALE_RANGE.1);
        let base = pattern.eps * scale;
        // The order exponent. See `OrderScale`.
        // **`gamma <= 1` is CHOSEN, not derived, and it is load-bearing.** It
        // makes the series diminishing — the third contact worth less than the
        // second. An earlier comment credited that to `element.rs` charging
        // strain for stacking contacts, which does not hold: that strain is
        // `sigma * cbrt(units^2)`, a *radial* cost of packing N units into a
        // sphere, and this same function argues 30 lines below that it "has
        // nothing to do with bonding to another one" — which is why it was
        // removed from capacity. A term cannot be irrelevant to a bond there and
        // be the mechanism for bond order here. Stated as chosen, on the
        // precedent of `element.rs`'s own `N^(2/3)` note, until something
        // derives it. `the_order_series_diminishes` asserts it.
        let order_scale =
            OrderScale::new(rng.next_f64_range(Self::GAMMA_RANGE.0, Self::GAMMA_RANGE.1));

        // **Capacity is contact density, not net binding energy**, and the
        // difference is the mechanism rather than a detail. `energy_per_unit` is
        // `eps * contacts/units - sigma * cbrt(units^2)`; the second term is the
        // strain of packing a cluster into a sphere, which is about that
        // cluster's own geometry and has nothing to do with bonding to another
        // one. A bond is a *contact*, so what a bond is made of is the contact
        // term alone.
        //
        // This retires `weakest_share`. That constant existed because a min-max
        // normalisation puts a zero at the table minimum and zero is absorbing
        // under a geometric mean; contact density is positive for every element
        // that can bond (`units >= 2` gives `contacts >= 1`), so positivity is
        // structural rather than floored. Its stated justification had also come
        // apart — "below 1/3" bought `ratio > 3`, and the only thing that ever
        // wanted `> 3` was a test deleted two rounds earlier. Measured, the
        // derived scheme's bondable ratio is 6.6-10.8 against the drawn scheme's
        // 2.3-5.3, which is wider — the direction §22.6 wants — and it varies
        // with `k`, which `1/weakest_share` cannot.
        //
        // The monomer has `contacts_upto(1) = 0` and so prices 0. That is
        // correct by meaning: no frontier, `valence == 0`, forms no bonds. It
        // does collide with the unknown-id sentinel, which is one more reason
        // Task 14 must retire that.
        let crate::element::ShellLaw::V1 { k } = pattern.law else {
            unreachable!(
                "generate_v1 is only called for a V1 table — table.physics() and \
                          the universe's own physics agree by construction, per generate_v1's \
                          own header comment in lib.rs"
            )
        };
        let pack = crate::packing::PackingConsts::new(k);
        let capacity: Vec<f64> = table
            .iter()
            .map(|(_, el)| contact_density(pack, el.units))
            .collect();

        let mut e = vec![Quanta::ZERO; n * n];
        for a in 0..n {
            // Hoisted: `a` is loop-invariant. Leaving the lookup inside made
            // LLVM version-clone the loop on an unreachable fallback — but only
            // *pre-LTO*; under `lto = "thin"`, the profile that ships, the clone
            // had already collapsed. What survives LTO is the reason to keep the
            // hoist: LLVM cannot prove `capacity: Vec<f64>` does not alias the
            // `e: Vec<Quanta>` being stored into, so it could not lift the
            // `capacity[a]` load out of the inner loop — one redundant load per
            // cell across ~4560 pairs. Legibility and one load, not speed; the
            // measured effect is nil at 4.4 us.
            let ca = capacity.get(a).copied().unwrap_or(0.0);
            for b in a..n {
                // Geometric mean, so the pair is symmetric bit-for-bit
                // (multiplication commutes exactly) and a weak partner drags the
                // bond down rather than being rescued by a strong one the way a
                // sum would allow — a bond needs both ends to hold. `sqrt` is
                // exactly specified by IEEE-754 and stays native (§13.1).
                // Both fallbacks are unreachable (`capacity.len() == n`, and
                // `a, b < n`) and both are `0.0` — which under contact density
                // is not an arbitrary sentinel but the value a non-bonding
                // element genuinely scores. An early version used `1.0`, the
                // *maximum*, thirty lines from a comment arguing the
                // conservative default at the other unreachable site.
                let cb = capacity.get(b).copied().unwrap_or(0.0);
                // **The association here is load-bearing and pinned (§13.4).**
                // It reads `base * (ca * cb).sqrt()`. Distributing it to
                // `base * ca.sqrt() * cb.sqrt()` is mathematically identical and
                // moves **47.7% of cells** by up to 5.2e-16 relative — measured
                // over 2e6 draws across the ranges these variables actually
                // occupy (`ca, cb` in [0.5, 8], `base` in [22.4, 100.8]).
                //
                // This replaces the pin the `weakest_share` normalisation
                // carried before it was retired. Deleting that comment with its
                // expression was right; leaving the file with **no** §13.4 pin
                // was not, and a review caught the gap.
                // `the_assembled_universe_digest_is_pinned` catches it.
                // Nobody tidies this later.
                let v = base * (ca * cb).sqrt();
                if let Some(cell) = e.get_mut(a * n + b) {
                    *cell = v;
                }
                if let Some(cell) = e.get_mut(b * n + a) {
                    *cell = v;
                }
            }
        }

        Self { n, e, order_scale }
    }

    /// V2's path (issue #26, Task 26.1 Steps 6-9, Decision 8): a rank-2
    /// construction, not V1's rank-1 `base * sqrt(ca * cb)`. The rank-1
    /// term is unchanged — `ca`/`cb` come from `capacity_scale +
    /// valence(el)`, the gap-derived quantity `orbital.rs` already
    /// computes and `element.rs`'s V2 path stores on
    /// [`crate::element::Element::valence`] — but a second, additive
    /// **ionic-excess** term now rides alongside it:
    /// `-ionic_base * affinity[a] * affinity[b]`, an outer product of the
    /// affinity vector with itself. See [`MigratedConstant::IonicScale`]
    /// for the full derivation (Pauling-motivated, numerically verified
    /// rank-2 and residual-correlation properties, and the worst-case
    /// positivity proof `bond_energy_stays_non_negative_at_every_reachable_rung`,
    /// below, checks directly).
    ///
    /// **`scale`/`order_scale`/`ionic_scale` are now P7-migrated too**
    /// (`MigratedConstant::Scale`/`Gamma`/`IonicScale`), so `rng` (the
    /// shared per-universe stream V1's own draws still use) has no reader
    /// left in this function — kept in the signature only because
    /// [`Self::generate`] passes it unconditionally to both version arms,
    /// per that function's own doc.
    fn generate_v2(table: &PeriodicTable, seed: u64, rung: Rung, _rng: &mut Stream) -> Self {
        let n = table.len();
        let pattern = table.pattern();

        let scale = draw_symmetric(seed, rung, MigratedConstant::Scale, 56.0);
        let base = pattern.eps * scale;
        let order_scale = OrderScale::new(draw_symmetric(seed, rung, MigratedConstant::Gamma, 0.6));

        let capacity_scale = draw_symmetric(
            seed,
            rung,
            MigratedConstant::CapacityScale,
            Self::BASE_CAPACITY_SCALE,
        );
        let ionic_scale = draw_symmetric(
            seed,
            rung,
            MigratedConstant::IonicScale,
            Self::BASE_IONIC_SCALE,
        );
        let ionic_base = pattern.eps * ionic_scale;

        let capacity: Vec<f64> = table
            .iter()
            .map(|(_, el)| capacity_scale + f64::from(el.valence))
            .collect();
        let affinity: Vec<f64> = table.iter().map(|(_, el)| el.affinity).collect();

        let mut e = vec![Quanta::ZERO; n * n];
        for a in 0..n {
            let ca = capacity.get(a).copied().unwrap_or(0.0);
            let aff_a = affinity.get(a).copied().unwrap_or(0.0);
            for b in a..n {
                let cb = capacity.get(b).copied().unwrap_or(0.0);
                let aff_b = affinity.get(b).copied().unwrap_or(0.0);
                // Rank-1 covalent term, plus the rank-2-completing ionic
                // term: opposite-signed (complementary) affinities ADD
                // energy, same-signed affinities SUBTRACT it, matching
                // Principle 1 ("bumps must meet hollows"). See
                // `MigratedConstant::IonicScale` for why this specific
                // additive form (not `(aff_a - aff_b)^2`) is what keeps
                // the whole matrix exactly rank 2.
                let v = base * (ca * cb).sqrt() - ionic_base * (aff_a * aff_b);
                if let Some(cell) = e.get_mut(a * n + b) {
                    *cell = v;
                }
                if let Some(cell) = e.get_mut(b * n + a) {
                    *cell = v;
                }
            }
        }

        Self { n, e, order_scale }
    }

    /// Bond dissociation energy for a bond of the given order.
    ///
    /// Takes ids rather than `&Element` deliberately: requiring `&Element`
    /// would couple the caller to the element table's lifetime and to a struct
    /// 48 bytes of which are `String`.
    ///
    /// **The caller is `Interner::intern` — per *species*, at intern time — and
    /// not a rate path**, which an earlier version of this comment claimed. The
    /// Gillespie loop reads a resolved `cleave_propensity`/`activation` and
    /// never reaches this method. That matters because the claim was the stated
    /// premise for the `n²` layout, and the premise was false — but so was the
    /// "0%" that replaced it, which held for one loop shape only. Measured
    /// against recomputing `base * sqrt(c_a * c_b)`, which is bit-identical over
    /// 3 087 174 cells: **+0.3%** in a single-accumulator loop (fadd-latency
    /// bound, so the `sqrt` hides entirely), **+18.6%** in the intern shape, and
    /// **+32.9%** with four independent accumulators over the whole table. The
    /// dense form is genuinely faster wherever the accumulator chain is broken,
    /// which is the shape a real consumer has.
    ///
    /// So the table stays on the measurement rather than on the extensibility
    /// story an earlier version told — that was a YAGNI argument for a design
    /// this same file rejects on §3.2 grounds. Caveat worth keeping: at a mean
    /// 67 KB the matrix fits this machine's 128 KiB L1D and does *not* fit the
    /// 32-48 KiB L1D of the x86-64 CI legs, where scattered access could shrink
    /// or invert that. Unmeasured; it changes nothing at intern-time rates.
    ///
    /// # An id this universe never minted scores [`Quanta::ZERO`], and that is
    /// the **worst** available answer, not a conservative one
    ///
    /// **Recorded as a known defect rather than defended.** An earlier version
    /// of this comment called zero "no bond, hence inert". That is backwards,
    /// and the spec says so directly: §9.4's spontaneous cleavage is
    /// `k = A·exp(−E_bond/T)` (PRD §9.1, §9.4), so `E_bond = 0` gives
    /// `exp(0) = 1` — the **maximum** of the range. A stray id does not produce
    /// an inert bond; it produces the most labile bond in the universe and the
    /// largest cleave propensity available.
    ///
    /// **This is present, not prospective, and by four orders of magnitude.**
    /// Measured over 500 universes against the slowest *formable* bond, the
    /// sentinel cleaves **3.8x to 4.36e4x** faster — so a stray `ElementId`
    /// dissolves a lineage in one step and looks like ordinary decay.
    ///
    /// Two earlier versions of this paragraph understated it. The first attached
    /// `1.23x` to the fastest bond when it was the slowest — real number, wrong
    /// referent, caught by two lanes. The second kept the seed-0 figures from
    /// before the derivation changed (`T = 206.09`, energies `6.60..42.71`) and
    /// called the hazard "entirely prospective... it becomes 10^2-10^5 once
    /// Task 14 raises `E/T`". `E/T` already rose, in the commit that derived the
    /// bond scale, and the figure is understated by roughly 5000x. A Task 14
    /// implementer reading "1.23x, prospective" would not prioritise closing
    /// this; at 4.4e4 it is not deferrable.
    ///
    /// The tell that a sentinel is the wrong shape here is that the *inert*
    /// value for `exp(−E/T)` would be `+∞`, which is a poison that propagates
    /// through any energy sum. `Option<Quanta>` is no better: the ids reaching
    /// here come from [`PeriodicTable::iter`] so `None` is unreachable in
    /// correct code, and a rate loop in a workspace that denies `unwrap` would
    /// write `.unwrap_or(Quanta::ZERO)` — reintroducing this exact value one
    /// level up, where it is harder to see.
    ///
    /// **The fix is to make the invalid case unrepresentable, and it belongs
    /// with the consumer.** An id validated against a table once, at the
    /// boundary, makes this lookup total and deletes the question.
    ///
    /// **Split the two halves of "prospective", because conflating them is what
    /// this whole comment has twice got wrong.** The *magnitude* is present
    /// today — 3.8x to 4.36e4x, measured above. What is prospective is only the
    /// *exposure*: nothing outside this crate's tests calls `energy` yet. An
    /// earlier version of this paragraph said "the hazard is entirely
    /// prospective" **four lines below** the paragraph that exists to retract
    /// exactly that phrase, so a Task 14 implementer reading to the end — which
    /// is where the what-do-I-do-about-it sentence lives — took away
    /// "deferrable" from the comment written to prevent it. It is not: the first
    /// consumer is Task 14's rate path and it inherits the full 4.4e4x on day
    /// one. That is why this is written into Task 14 (reaction rates) rather
    /// than left here. `the_zero_sentinel_is_the_fastest_cleaving_value`
    /// pins the current value *and* names it as wrong, so a future reader
    /// cannot mistake a green test for a blessing.
    #[must_use]
    pub fn energy(&self, a: ElementId, b: ElementId, order: BondOrder) -> Quanta {
        let (ia, ib) = (a.index(), b.index());
        if ia >= self.n || ib >= self.n {
            return Quanta::ZERO;
        }
        self.e
            .get(ia * self.n + ib)
            .copied()
            .unwrap_or(Quanta::ZERO)
            * self.order_scale.of(order)
    }

    /// Elements this matrix was built for.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.n
    }

    /// Whether the table this was built from was empty.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.n == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PhysicsVersion;
    use crate::element::{ElementId, PeriodicTable, generate_elements};
    use crate::perturbation::Rung;
    use borbax_rng::{Domain, Stream};

    /// V2's own table and bond matrix, for Task 26.1 Step 6's acceptance
    /// tests below -- `universe(seed)` deliberately fixes
    /// `PhysicsVersion::CURRENT`, which is V1.
    fn universe_v2(seed: u64) -> (PeriodicTable, BondEnergyMatrix) {
        let table = generate_elements(seed, PhysicsVersion::V2);
        let mut rng = Stream::new(seed, Domain::Universe, 2);
        let bonds = BondEnergyMatrix::generate(&table, seed, Rung::draw(seed), &mut rng);
        (table, bonds)
    }

    /// One universe's table and its bond matrix, drawn the way
    /// [`crate::Universe::generate`] draws them.
    fn universe(seed: u64) -> (PeriodicTable, BondEnergyMatrix) {
        let table = generate_elements(seed, PhysicsVersion::CURRENT);
        let mut rng = Stream::new(seed, Domain::Universe, 2);
        let bonds = BondEnergyMatrix::generate(&table, seed, Rung::draw(seed), &mut rng);
        (table, bonds)
    }

    fn ids(table: &PeriodicTable) -> Vec<ElementId> {
        table.iter().map(|(id, _)| id).collect()
    }

    /// Elements with at least one contact — everything but the monomer, which
    /// is the only element whose capacity is 0. This is the predicate for the
    /// **zero-cell** hazard; `bondable_ids` is the predicate for "a molecule can
    /// contain it". Conflating them excluded 404 well-behaved elements per 200
    /// seeds on a reason that was false, and a review found the comment saying
    /// so sitting directly above code that still used the wrong one.
    fn contacting_ids(table: &PeriodicTable) -> Vec<ElementId> {
        table
            .iter()
            .filter(|(_, e)| e.units > 1)
            .map(|(id, _)| id)
            .collect()
    }

    /// Elements a molecule can actually contain. `valence == 0` means a closed
    /// shell — no frontier, so no bonding slots (§7.1).
    fn bondable_ids(table: &PeriodicTable) -> Vec<ElementId> {
        table
            .iter()
            .filter(|(_, e)| e.valence >= 1)
            .map(|(id, _)| id)
            .collect()
    }

    /// Bit-exact, not approximate. `energy` reads one symmetric table, so a
    /// tolerance here would pass for an implementation that computed the two
    /// directions separately and rounded them differently.
    ///
    /// **Raw `to_bits()`, not `canonical_bits()`** — the latter is lossy by
    /// design (it ties `+0.0` with `-0.0` and every NaN with every other NaN),
    /// which is right for a product hash and wrong for a detector. Under
    /// `canonical_bits` a matrix whose two directions were `+0.0` and `-0.0`
    /// would pass a test whose own doc claims bit-exactness.
    #[test]
    fn energies_are_symmetric() {
        let (table, m) = universe(4);
        for &a in &ids(&table) {
            for &b in &ids(&table) {
                assert_eq!(
                    m.energy(a, b, BondOrder::SINGLE).get().to_bits(),
                    m.energy(b, a, BondOrder::SINGLE).get().to_bits(),
                    "{a:?} vs {b:?}"
                );
            }
        }
    }

    /// **Guards this module's reason for ignoring `affinity` and `group`, not
    /// its behaviour.** Both figures in the header are claims about the table
    /// Task 4 generates, and a claim of that shape is exactly what goes stale:
    /// the `bond_energies` derivation would still compile, still pass every
    /// other test here, and its stated justification would silently be false.
    ///
    /// If this fails, the header is wrong and the design question is genuinely
    /// reopened — an `affinity` with a negative half *would* support §8.3's
    /// opposition term, and a `group` with real multiplicity *would* be a
    /// family worth keying on.
    #[test]
    fn the_table_still_justifies_ignoring_affinity_and_group() {
        let mut aff_lo = f64::MAX;
        let mut negatives = 0_usize;
        let mut worst_multiplicity = 0_usize;
        let mut counted = 0_usize;

        for seed in 0..64 {
            let table = generate_elements(seed, PhysicsVersion::CURRENT);
            let mut per_group = [0_usize; 256];
            for (_, e) in table.iter() {
                counted += 1;
                if e.affinity < aff_lo {
                    aff_lo = e.affinity;
                }
                if e.affinity < 0.0 {
                    negatives += 1;
                }
                if let Some(slot) = per_group.get_mut(usize::from(e.group)) {
                    *slot += 1;
                }
            }
            for c in per_group {
                if c > worst_multiplicity {
                    worst_multiplicity = c;
                }
            }
        }

        assert_eq!(
            negatives, 0,
            "affinity reached the negative half over {counted} elements — §8.3's \
             opposition term is now expressible and this module's header is stale"
        );
        assert!(
            aff_lo > 0.14 && aff_lo < 0.141,
            "affinity's floor moved to {aff_lo}; the header quotes 0.14047619047619042"
        );
        assert!(
            worst_multiplicity <= 4,
            "a group now holds {worst_multiplicity} elements; the header claims at most 4, \
             and above that `group` starts to be a family worth keying on"
        );
    }

    /// **The other half of the derivation claim, and the half six mutation
    /// probes missed.**
    ///
    /// **Complementary to `perturbing_...`, not a superset of it — do not delete
    /// that one as redundant.** Gap-adaptive dilution (each element free
    /// anywhere strictly between its neighbours' midpoints, extremes pinned)
    /// *passes* this test, because ordering is preserved by construction, and is
    /// caught only by the perturbation test's untouched half: the jitter bounds
    /// are relational, so bumping the target shifts its neighbours. Measured
    /// dilution threshold with both in place: 2e-6 < eps <= 3e-6, against the
    /// 0.99 that defeated the suite before.
    ///
    /// `perturbing_one_elements_binding_energy_moves_every_bond_it_makes` says
    /// `energy_per_unit` has *non-zero* influence. It does not say it has *all*
    /// the influence, and the difference is not academic: a review measured that
    /// `capacity = 0.01·derived + 0.99·drawn`, with the two extreme elements
    /// pinned so the spread identity survives, **passes the entire suite**.
    /// Strict monotonicity in the bumped element holds at any positive weight,
    /// so the perturbation test cannot see dilution — only replacement. Every
    /// probe run before this commit tested replacement.
    ///
    /// This asserts exactness instead: self-bond energy is `base · capacity`, so
    /// ordering self-bonds must reproduce the ordering of `energy_per_unit`
    /// **exactly**, over every ordered pair. Any drawn admixture inverts roughly
    /// half of the pairs it touches.
    #[test]
    fn capacity_is_a_function_of_contact_density_and_of_nothing_else() {
        for seed in 0..24 {
            let (table, m) = universe(seed);
            // **Contact density, not `energy_per_unit`.** The derivation input
            // changed when the bond scale moved onto `eps`: a bond is a contact,
            // so what a bond is made of is the contact term, and the strain term
            // `energy_per_unit` subtracts is about a cluster's own packing. The
            // anti-dilution property this test exists for is unchanged — it
            // still fails the moment anything else reaches the capacity scale.
            let crate::element::ShellLaw::V1 { k } = table.pattern().law else {
                unreachable!("this test targets V1")
            };
            let pack = crate::packing::PackingConsts::new(k);
            let all: Vec<(ElementId, f64)> = table
                .iter()
                .map(|(id, e)| (id, contact_density(pack, e.units)))
                .collect();
            let mut inversions = 0_usize;
            for &(ia, ea) in &all {
                for &(ib, eb) in &all {
                    let self_a = m.energy(ia, ia, BondOrder::SINGLE).get();
                    let self_b = m.energy(ib, ib, BondOrder::SINGLE).get();
                    // Strict on both sides, so ties must correspond to ties.
                    if (ea < eb) != (self_a < self_b) {
                        inversions += 1;
                    }
                }
            }
            assert_eq!(
                inversions, 0,
                "seed {seed}: {inversions} pairs order differently by binding \
                 energy and by self-bond energy — something other than contact \
                 density is reaching the capacity scale"
            );
        }
    }

    /// **Restores what the deleted spread identity pinned: the combiner.**
    ///
    /// `the_spread_is_exactly_the_reciprocal_of_the_weakest_share` died with
    /// `weakest_share`, and a review measured what went with it — swapping the
    /// geometric mean for `(c_a c_b)^0.7` passes **60 of 61** tests, failing
    /// only the digest, which is re-pinned by hand on every deliberate change
    /// and so is not a guard against this.
    ///
    /// The identity survives `weakest_share`'s retirement in a base-free form.
    /// The extreme cells of the bondable submatrix are the extreme self-bonds
    /// (the test above), and a self-bond is `base * c`, so the whole-submatrix
    /// ratio is exactly `c_max / c_min` — with `base` cancelling. Under
    /// `(c_a c_b)^0.7` it is that ratio to the 1.4; under `base * c_a * c_b`, to
    /// the 2.
    #[test]
    fn the_bondable_spread_is_exactly_the_capacity_ratio() {
        for seed in 0..16 {
            let (table, m) = universe(seed);
            let crate::element::ShellLaw::V1 { k } = table.pattern().law else {
                unreachable!("this test targets V1")
            };
            let pack = crate::packing::PackingConsts::new(k);
            let bondable: Vec<(ElementId, f64)> = table
                .iter()
                .filter(|(_, e)| e.valence >= 1)
                .map(|(id, e)| (id, contact_density(pack, e.units)))
                .collect();
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            let (mut clo, mut chi) = (f64::MAX, f64::MIN);
            for &(a, ca) in &bondable {
                if ca < clo {
                    clo = ca;
                }
                if ca > chi {
                    chi = ca;
                }
                for &(b, _) in &bondable {
                    let v = m.energy(a, b, BondOrder::SINGLE).get();
                    if v < lo {
                        lo = v;
                    }
                    if v > hi {
                        hi = v;
                    }
                }
            }
            let (got, want) = (hi / lo, chi / clo);
            assert!(
                (got - want).abs() < 1e-12 * want,
                "seed {seed}: bondable spread {got} is not the capacity ratio {want} \
                 — the combiner is no longer the geometric mean"
            );
        }
    }

    /// **The order series must diminish, and `gamma <= 1` is what makes it.**
    ///
    /// A review measured that widening the draw to `[0.68, 1.6]` — admitting the
    /// accelerating series this change exists to remove — passes 60 of 61 tests.
    /// The choice had no tripwire, and its cited mechanism was one this file
    /// removed from capacity in the same commit.
    #[test]
    fn the_order_series_diminishes() {
        for seed in 0..24 {
            let (_, m) = universe(seed);
            let step = |n: u8| {
                let o = BondOrder::new(n).unwrap_or(BondOrder::SINGLE);
                m.order_scale.of(o)
            };
            for n in 2..BondOrder::MAX {
                let (prev, this) = (step(n) - step(n - 1), step(n + 1) - step(n));
                assert!(
                    this < prev,
                    "seed {seed}: order {n}->{} adds {this}, more than {}->{n}'s {prev} \
                     — the series is accelerating, which `gamma <= 1` exists to prevent",
                    n + 1,
                    n - 1
                );
            }
        }
    }

    /// **Pins the combiner, which nothing else does.**
    ///
    /// `bond = base · sqrt(c_a · c_b)` factorises as `x_a · x_b`, so the matrix
    /// is exactly rank 1 and every 2x2 minor vanishes. That identity is specific
    /// to the geometric mean: a review measured the worst relative minor at
    /// 6.21e-16 here against 6.06e-1 for the arithmetic mean, and confirmed that
    /// swapping in the arithmetic mean passes all nine of the other tests.
    /// Harmonic and `min` are likewise invisible to everything except strict
    /// monotonicity.
    ///
    /// The identity survives a change to the capacity *map* — that only moves
    /// `c` — so this pins **separability** without pinning the normalisation.
    /// The geometric mean is pinned by this together with the spread identity,
    /// and by neither alone: `base * c_a * c_b` and `(c_a c_b)^0.7` both keep
    /// rank 1 and are caught by the other test.
    #[test]
    fn the_matrix_is_exactly_rank_one() {
        for seed in 0..8 {
            let (table, m) = universe(seed);
            // **`units > 1`, not `valence >= 1`** — the zero cell belongs to the
            // *monomer*, not to every closed shell. Measured over 200 seeds: of
            // 604 `valence == 0` elements only 200 have zero capacity (one per
            // seed, `units == 1`); the other 404 have positive capacity up to
            // 5.09, and a `valence == 0` element holds the whole-table maximum
            // in 63 of 200 seeds. Filtering on valence excluded them from
            // coverage on a reason that was false.
            let ids: Vec<ElementId> = contacting_ids(&table).into_iter().take(16).collect();
            let mut worst = 0.0_f64;
            for &a in &ids {
                for &b in &ids {
                    for &c in &ids {
                        for &d in &ids {
                            let (ab, cd) = (
                                m.energy(a, b, BondOrder::SINGLE).get(),
                                m.energy(c, d, BondOrder::SINGLE).get(),
                            );
                            let (ad, cb) = (
                                m.energy(a, d, BondOrder::SINGLE).get(),
                                m.energy(c, b, BondOrder::SINGLE).get(),
                            );
                            // Guarded: `0.0 / 0.0` is NaN and `NaN > worst`
                            // is false, so an unguarded ratio would *silently
                            // drop* zero-valued cells rather than catch them —
                            // and `energies_are_positive_and_finite`, the
                            // precondition, runs on seed 9 while this runs 0..8.
                            assert!(ab * cd > 0.0, "seed {seed}: zero-valued cell");
                            let minor = (ab * cd - ad * cb).abs();
                            let rel = minor / (ab * cd);
                            if rel > worst {
                                worst = rel;
                            }
                        }
                    }
                }
            }
            assert!(
                worst < 1e-12,
                "seed {seed}: worst relative 2x2 minor {worst} — the matrix is no \
                 longer separable. This does NOT mean the combiner stopped being a \
                 geometric mean: a non-geometric separable combiner keeps rank 1, and \
                 a geometric mean plus a derived pair factor breaks it. If a pair \
                 factor was added deliberately, what rank 1 buys is that the sign of \
                 a substitution does not depend on the spectator. It is NOT what keeps \
                 a cycle from being net-favourable — bond energy is a state function of \
                 the graph for any symmetric matrix, so no cycle is net-favourable \
                 regardless, and treating this as a free-energy-pump guard would \
                 wrongly veto a legitimate change"
            );
        }
    }

    /// Over every pair, not one hand-picked pair: the order multipliers are
    /// universe-wide, so a single pair cannot distinguish "orders increase"
    /// from "these two elements happen to increase".
    #[test]
    fn higher_orders_are_stronger() {
        let (table, m) = universe(4);
        let bondable = bondable_ids(&table);
        for &a in &bondable {
            for &b in &bondable {
                let (s, d, t) = (
                    m.energy(a, b, BondOrder::SINGLE),
                    m.energy(a, b, BondOrder::new(2).unwrap_or(BondOrder::SINGLE)),
                    m.energy(a, b, BondOrder::new(3).unwrap_or(BondOrder::SINGLE)),
                );
                assert!(d > s && t > d, "{a:?}-{b:?}: {s:?} {d:?} {t:?}");
            }
        }
    }

    /// Over pairs with `units > 1`.
    ///
    /// **The monomer, and only the monomer, prices 0** — `contacts_upto(1)` is
    /// 0, so its capacity is 0 and every cell in its row and column is
    /// `Quanta::ZERO`. That is correct by meaning: no contacts, no bond.
    ///
    /// An earlier version filtered on `valence >= 1` and said a closed shell
    /// prices 0. Measured, that is false for 404 of 604 such elements, some of
    /// which carry the table's largest capacity — so the filter was excluding
    /// well-behaved cells on a wrong reason.
    #[test]
    fn energies_are_positive_and_finite() {
        let (table, m) = universe(9);
        let bondable = bondable_ids(&table);
        for &a in &bondable {
            for &b in &bondable {
                let e = m.energy(a, b, BondOrder::SINGLE);
                assert!(e.get() > 0.0 && e.is_finite(), "{a:?}-{b:?} = {e:?}");
            }
        }
    }

    /// The bondable submatrix's extremes are the self-bonds of the extreme
    /// bondable elements.
    ///
    /// **Fourth attempt, and the first that reads the bond matrix.** The three
    /// before it all tried to assert a *magnitude* — `hi/lo > 3.0` over the
    /// whole matrix (endpoints are `valence == 0` elements no molecule can
    /// contain), then the bondable ratio in `[2.4, 4.6]` (fitted to 200 seeds,
    /// false from seed 213), then a coverage fraction `> 0.6` (fitted to 2000).
    /// Round 3 replaced the third with "bondable elements straddle the median",
    /// which carries no constant and is nearly vacuous: 95–97.5% of elements are
    /// bondable with ~45% each side, so it fires only once the non-bondable
    /// fraction reaches 45% on one side. Seed 1209's entire above-median
    /// bondable set sits within 1.04% of the range and it passes.
    ///
    /// **The lesson, which is the useful part: a magnitude claim cannot be made
    /// constant-free, because a magnitude needs a scale.** Three thresholds and
    /// one vacuity is the evidence. So the magnitude is *recorded* below and
    /// asserted nowhere — measured coverage `[0.7157, 0.9018]` over 2000
    /// universes — and what is asserted instead is a *derivation* claim, which
    /// can be exact: the extreme bondable bond energies are the self-bonds of
    /// the bondable elements with extreme `energy_per_unit`. Bit-exact, 0
    /// violations over 2000 seeds.
    ///
    /// **What it does not catch: a flat matrix.** With every bond equal the
    /// extremes tie and the assertion holds trivially. That is caught twice
    /// elsewhere — by `the_single_bond_rate_spread...`'s `hi > lo` guard and by
    /// `capacity_is_a_function_of_contact_density_and_of_nothing_else` — so the
    /// suite covers it; this test does not, and saying so is the point.
    ///
    /// That also repairs a defect round 3 introduced: the previous body computed
    /// `lo`/`hi` over the submatrix and then wrote `let _ = (lo, hi);`, so the
    /// test named for bond coverage, living in `bonds.rs`, **did not read the
    /// bond matrix at all** — verified, it passed with every bond energy
    /// replaced by a constant.
    #[test]
    fn the_bondable_extremes_are_the_extreme_elements_self_bonds() {
        for seed in 0..20 {
            let (table, m) = universe(seed);
            let bondable: Vec<ElementId> = table
                .iter()
                .filter(|(_, e)| e.valence >= 1)
                .map(|(id, _)| id)
                .collect();
            assert!(
                bondable.len() >= 2,
                "seed {seed}: fewer than two bondable elements"
            );
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for &a in &bondable {
                for &b in &bondable {
                    let v = m.energy(a, b, BondOrder::SINGLE).get();
                    if v < lo {
                        lo = v;
                    }
                    if v > hi {
                        hi = v;
                    }
                }
            }
            // The extreme bondable elements by binding energy. `canonical_cmp`
            // rather than `total_cmp` under an `#[expect]`: it is the route
            // `clippy.toml`'s own reason string names, and the exemption's
            // stated justification was false anyway — no test in the workspace
            // asserts `energy_per_unit` is finite, and the nearest one asserts a
            // different quantity on seed 9 while this ran seeds 0..20.
            let crate::element::ShellLaw::V1 { k } = table.pattern().law else {
                unreachable!("this test targets V1")
            };
            let pack = crate::packing::PackingConsts::new(k);
            let density = |id: ElementId| {
                table
                    .get(id)
                    .map_or(0.0, |x| contact_density(pack, x.units))
            };
            let extreme = |pick_max: bool| {
                bondable
                    .iter()
                    .copied()
                    .fold(None::<(ElementId, f64)>, |best, id| {
                        let e = density(id);
                        match best {
                            Some((_, b)) if (e > b) != pick_max => best,
                            _ => Some((id, e)),
                        }
                    })
            };
            let (weak, strong) = (extreme(false), extreme(true));
            // `assert!` then `if let`, not `let ... else { panic! }`:
            // `clippy::panic` is denied workspace-wide and reaches inside tests.
            assert!(
                weak.is_some() && strong.is_some(),
                "seed {seed}: bondable set has no extremes"
            );
            if let (Some((weak, _)), Some((strong, _))) = (weak, strong) {
                // The derivation claim, exact: `bond = base * sqrt(c_a * c_b)` is
                // monotone in both arguments, so the extreme cells of the bondable
                // submatrix must be the self-bonds of its extreme elements. This
                // reads the matrix, which the previous version had stopped doing.
                assert_eq!(
                    m.energy(weak, weak, BondOrder::SINGLE).get().to_bits(),
                    lo.to_bits(),
                    "seed {seed}: the weakest bondable bond is not the weakest \
                     bondable element's self-bond"
                );
                assert_eq!(
                    m.energy(strong, strong, BondOrder::SINGLE).get().to_bits(),
                    hi.to_bits(),
                    "seed {seed}: the strongest bondable bond is not the strongest \
                     bondable element's self-bond"
                );
            }
        }
    }

    /// **Both bond constants are this universe's own draws, bit-exactly.**
    ///
    /// Measured: replacing `scale` with a literal passes **57 of 58** tests, and
    /// pinning `gamma` at `1 - 1e-9` passes 57 of 58 — in both cases the only
    /// failure is the digest, which is re-pinned by hand on every deliberate
    /// change and so is not a guard against this. The energy bands cannot see
    /// either: under a constant scale they read `e_lo = 22.4` against a band top
    /// of 25, and `e_hi = 1861` against a floor of 1200.
    ///
    /// `the_order_series_diminishes` cannot see the second, because it asserts
    /// only the *sign* of the second difference, which survives to
    /// `gamma = 1 - 2e-16`. The header's whole argument for `gamma <= 1` is dead
    /// at nine decimal places of linear with nothing failing.
    ///
    /// Forward direction, not a recovery, so there is no tolerance — which also
    /// makes this the §13.4 association pin for `base * (ca * cb).sqrt()`, today
    /// caught only by the digest.
    ///
    /// **It reads the real universe, not a re-typed copy of it.** An earlier
    /// version opened its own `Stream::new(seed, Domain::Universe, 2)` and so
    /// pinned "this test's index equals this test's index" — both copies inside
    /// `#[cfg(test)]`. Measured: changing the *production* index at
    /// `lib.rs`'s `generate_under` from 2 to 3 left this green and fired only
    /// the digest, while `lib.rs`'s own header explains that this index exists
    /// so an added bond constant cannot shift `temp_min` and with it the decay
    /// band — the most sensitive parameter in the system.
    #[test]
    fn both_bond_constants_are_this_universes_own_draws() {
        for seed in 0..24 {
            let u = crate::Universe::generate(seed);
            let (table, m) = (&u.table, &u.bonds);
            let crate::element::ShellLaw::V1 { k } = table.pattern().law else {
                unreachable!("this test targets V1")
            };
            let pack = crate::packing::PackingConsts::new(k);
            // Same stream, same index, same order — and the index is the
            // production one, because `Universe::generate` above is what chose
            // it. A reordered draw or an inserted one fails here; a *retuned
            // range* does not, by construction, and is the digest's job (see
            // `SCALE_RANGE`).
            let mut rng = Stream::new(seed, Domain::Universe, 2);
            let scale = rng.next_f64_range(
                BondEnergyMatrix::SCALE_RANGE.0,
                BondEnergyMatrix::SCALE_RANGE.1,
            );
            let gamma = rng.next_f64_range(
                BondEnergyMatrix::GAMMA_RANGE.0,
                BondEnergyMatrix::GAMMA_RANGE.1,
            );
            let base = table.pattern().eps * scale;
            // **Every ordered pair, not the diagonal.** On the diagonal
            // `(c * c).sqrt()` is exactly `c` — 2e6 of 2e6 samples over the
            // range contact density occupies, the standard binary-FP result
            // that `sqrt(fl(x^2)) == |x|` absent over/underflow — so a
            // diagonal-only assertion never exercises the sqrt-of-product at
            // all, which is the association this test claims to pin. Measured:
            // reassociating **only** off-diagonal cells to
            // `base * ca.sqrt() * cb.sqrt()` left this green and fired the
            // digest alone.
            let cs: Vec<(ElementId, f64)> = table
                .iter()
                .map(|(id, el)| (id, contact_density(pack, el.units)))
                .collect();
            for &(ia, ca) in &cs {
                for &(ib, cb) in &cs {
                    assert_eq!(
                        m.energy(ia, ib, BondOrder::SINGLE).get().to_bits(),
                        (base * (ca * cb).sqrt()).get().to_bits(),
                        "seed {seed}: the bond scale is no longer `eps * <the draw>`"
                    );
                }
            }
            // **Every order, not order 2.** The predecessor read
            // `BondOrder::new(2).unwrap_or(BondOrder::SINGLE)`, which would turn
            // "2 stopped being a valid order" into a *different* assertion
            // firing with a message naming the wrong cause. `unwrap` is not the
            // fix — `clippy::unwrap_used` is denied workspace-wide and reaches
            // inside `#[cfg(test)]`. Iterating `ALL` needs no fallible
            // construction at all and pins the whole series rather than one term.
            for o in BondOrder::ALL {
                assert_eq!(
                    m.order_scale.of(o).to_bits(),
                    det_math::powf(f64::from(u8::from(o)), gamma).to_bits(),
                    "seed {seed}: order {o:?}'s multiplier is no longer the drawn gamma"
                );
            }
            assert_eq!(
                m.order_scale.of(BondOrder::SINGLE).to_bits(),
                1.0_f64.to_bits(),
                "seed {seed}: the order exponent is no longer the drawn gamma"
            );
        }
    }

    /// **The header's figures, asserted rather than quoted.**
    ///
    /// Three rounds running, a measured table went into the module header and
    /// the next commit falsified it — most recently because moving the bond
    /// generator to its own stream index shifted every temperature. Each time a
    /// review had to find it, and each time the stale numbers were the sizing
    /// input for Task 14's coupling.
    ///
    /// The bands are wide on purpose: this is a drift detector for figures the
    /// header describes qualitatively, not a physics gate. A change that moves
    /// any of them by more than a factor of two is a change the header must be
    /// rewritten for, and this says so by name.
    #[test]
    fn the_header_figures_are_current() {
        let (mut e_lo, mut e_hi) = (f64::MAX, f64::MIN);
        let (mut ln_rate_hi, mut et_hi) = (f64::MIN, f64::MIN);
        let mut r_min = f64::MAX;
        // **500, because that is the sample the header quotes.** It ran 200
        // while the header attributed 500-seed figures to it, and three of the
        // five disagreed — `e_lo` 11.81 vs 11.40, `e_hi` 2361.69 vs 2370.55,
        // and the per-seed rate floor 4.15 vs 3.58, the last off by 16%. That
        // is the same stale-figure class the paragraph immediately above this
        // test declares closed, committed in the act of declaring it.
        for seed in 0..500 {
            let u = crate::Universe::generate(seed);
            let b: Vec<(ElementId, u8)> = u
                .table
                .iter()
                .filter(|(_, e)| e.valence >= 1)
                .map(|(id, e)| (id, e.valence))
                .collect();
            let t_min = u.consts.temp_min.get();
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for &(x, vx) in &b {
                for &(y, vy) in &b {
                    // **Only orders this pair can actually form.** A bond of
                    // order `n` consumes `n` slots at each end, so the cap is
                    // `min(valence_a, valence_b)` — `BondOrder`'s own doc says
                    // so, and an earlier version of this test ignored it.
                    // Measured over 20 000 seeds the per-universe ceiling is 4
                    // in 27.2%, 5 in 50.9%, 6 in 21.9%: pricing every pair at
                    // order 6 reports an energy **no bond in 78.1% of universes
                    // can carry**, overstating `e_hi` by up to **1.6310x**
                    // per-seed worst case (and `E/T` by the same factor, since
                    // `T` is fixed within a seed). An earlier version added "and
                    // the rate-ratio margin by 1.35x", which matches none of the
                    // three readings — ratio of rate ratios is 119.83x, of `ln`
                    // rate ratios 1.6533x, of `E/T` 1.6310x. An unnamed "margin"
                    // in a comment Task 14 sizes against is exactly the class of
                    // figure this file has lost four rounds to.
                    // `u8::min`, not the comparison form. The `if a < b` idiom is house style
                    // for `f64` because `f64::min` is `disallowed_methods` for a stated
                    // IEEE-754 reason; writing it on integers erodes exactly that
                    // distinction, which `clippy.toml`'s header insists on.
                    let cap = vx.min(vy);
                    for o in BondOrder::ALL.into_iter().filter(|o| u8::from(*o) <= cap) {
                        let v = u.bonds.energy(x, y, o).get();
                        // F5: `<` and `>` are both false for NaN, so a bare
                        // min/max scan *drops* a poisoned cell rather than
                        // reporting it — the same hazard
                        // `the_matrix_is_exactly_rank_one` records for zero
                        // cells. Measured: a NaN at order 5 left this test
                        // green; at order 6 it failed by a 13% margin, on luck.
                        //
                        // Verified by injection at every order — a NaN in
                        // `OrderScale::mult[k]` fails this assert for all
                        // k in 1..=6, order 6 failing on the assert itself.
                        // **Coverage note:** the formable cap means order 6 is
                        // now read in only 21.9% of universes (0.59% of pairs)
                        // and order 5 in 72.8%, where before it was every order
                        // for every pair. Sound today because every NaN source
                        // here is universe-wide — `contact_density`'s branch and
                        // a NaN `gamma` both poison all orders. If `mult` ever
                        // becomes per-universe conditional, scan all of `ALL`
                        // for finiteness and apply the cap only to the min/max.
                        assert!(v.is_finite(), "seed {seed}: {x:?}-{y:?} {o:?} is {v}");
                        if v < lo {
                            lo = v;
                        }
                        if v > hi {
                            hi = v;
                        }
                    }
                }
            }
            if lo < e_lo {
                e_lo = lo;
            }
            if hi > e_hi {
                e_hi = hi;
            }
            // Kept unexponentiated. `exp` of it is what the header quotes, but
            // the *assertion* belongs on the log — see the band block below.
            let ln_r = (hi - lo) / t_min;
            if ln_r > ln_rate_hi {
                ln_rate_hi = ln_r;
            }
            if hi / t_min > et_hi {
                et_hi = hi / t_min;
            }
            let r = borbax_units::det_math::exp(ln_r);
            if r < r_min {
                r_min = r;
            }
            // **A tail claim, so it degrades with sample size by construction**
            // — it asserts a property of *every* universe from a prefix of the
            // seed space. Measured minima: 4.15 over 200 seeds, 3.58 over 500,
            // 2.40 over 1000, 2.35 over 10 000, 2.16 over 200 000. It does not
            // cross 2.0 within 200 000, but the margin is 42% smaller at 1000
            // than at 200, so raising this loop is not a free edit and the
            // header must not state this as a property of the physics.
            assert!(
                r > 2.0,
                "seed {seed}: cleave-rate ratio {r} across formable orders is below 2.0 — \
                 the header says every drawn universe exceeds it. This is a tail \
                 statistic: check the sample size before concluding the physics moved"
            );
        }
        assert!(
            (5.0..=25.0).contains(&e_lo) && (1200.0..=5200.0).contains(&e_hi),
            "bondable energies [{e_lo}, {e_hi}] left the header's range by more than \
             a factor of two — rewrite the header, do not widen this"
        );
        assert!(
            (7.0..=30.0).contains(&et_hi),
            "max E/T is {et_hi}; the header describes a regime, and this moved out of it"
        );
        // **On the log, because the header's factor-of-two contract lives on the
        // exponent.** `rate_hi > 1e4` is `exp(ln_rate_hi)`, so it fires at a
        // 1.138x change in the energy scale, not a 2x one — bisected under a
        // uniform scaling α of every bond energy, which is exactly what
        // retuning the `scale` draw does: the composite window across all four
        // bands was α ∈ (0.8789, 2.1176), and `rate_hi` was the binding edge on
        // the low side. Stating a 2x tolerance while enforcing 1.14x is how a
        // pure-retune commit comes back as a red suite with no physics in it.
        assert!(
            (5.0..=21.0).contains(&ln_rate_hi),
            "max ln(cleave-rate ratio) {ln_rate_hi} (ratio {}) — a factor-of-two band on \
             the exponent, which is the scale the header's tolerance is stated on",
            borbax_units::det_math::exp(ln_rate_hi)
        );
        // Emitted, not asserted: the header quotes these and has now been
        // falsified four rounds running by quoting a different sample than the
        // test ran. Regenerate the header block from this line's output.
        println!(
            "header figures @500 seeds: e_lo={e_lo:.2} e_hi={e_hi:.2} \
             rate_hi={:.4e} et_hi={et_hi:.2} per_seed_r_min={r_min:.4}",
            borbax_units::det_math::exp(ln_rate_hi)
        );
    }

    /// **Records the defect four review lanes converged on — and is explicitly
    /// NOT the tripwire for its repair.**
    ///
    /// §22.6 wants "some linkages durable and others fragile", and §9.4 calls
    /// differential persistence the thing selection acts on. Both are claims
    /// about *rates*, and §9.1/§9.5 make the rate `k = A·exp(−E_bond/T)`. This
    /// crate mints both scales — `base = eps * scale`, measured `[22.5, 100.5]`
    /// Quanta, and `T ∈ [180, 760]`
    /// Thermal — so `E/T` stays in `exp`'s linear regime. The energy spread is
    /// fine; the *rate* spread is not.
    ///
    /// **A round-2 review killed the claim this doc used to make.** It said the
    /// assertion "fails when Task 14 fixes the scale", which cannot be true: the
    /// coupling belongs in `borbax-reaction`'s `rate()`, and applying it there
    /// leaves this crate's energies and temperatures untouched. The test would
    /// have stayed green with the defect repaired and a stale defect-record
    /// passing beside it — a green test whose doc describes a problem that no
    /// longer exists, which is worse than an honest one because nothing tells
    /// the next reader.
    ///
    /// So this asserts only what `borbax-universe` can see: the two scales *as
    /// minted here* are mismatched. It guards against those drifting, nothing
    /// more. Deleting it is a numbered step in Task 14, because no assertion in
    /// this crate can observe a constant in another one.
    ///
    /// The slice is stated rather than implied, which is the other round-2
    /// correction: **weakest to strongest bondable SINGLE bond at
    /// mid-temperature**. The module header used to quote this figure under the
    /// words "most durable bondable bond", which is the highest **formable**
    /// order — a different and much wider number: all formable orders at
    /// `temp_min` reaches **1.609e6** and exceeds 2.0 in **2000 of 2000**
    /// universes. Sizing Task 14's coupling from the narrow slice while reading
    /// the wide sentence would oversize it by orders of magnitude.
    ///
    /// (That parenthetical previously read "a *triple* bond … reaches 3.50 …
    /// 716 of 2000". "Triple" names the `Single/Double/Triple` enum this branch
    /// retired — the widest order is now 6 — and both figures predate the
    /// derivation change. A stale number inside the sentence warning against
    /// stale numbers.)
    #[test]
    fn the_single_bond_rate_spread_at_mid_temperature_is_narrow() {
        for seed in 0..20 {
            let u = crate::Universe::generate(seed);
            let bondable: Vec<ElementId> = u
                .table
                .iter()
                .filter(|(_, e)| e.valence >= 1)
                .map(|(id, _)| id)
                .collect();
            let t = f64::midpoint(u.consts.temp_min.get(), u.consts.temp_max.get());
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for &a in &bondable {
                for &b in &bondable {
                    let e = u.bonds.energy(a, b, BondOrder::SINGLE).get();
                    if e < lo {
                        lo = e;
                    }
                    if e > hi {
                        hi = e;
                    }
                }
            }
            // **`exp((hi-lo)/T) < hi/lo` reduces EXACTLY to `logmean(lo,hi) < T`,
            // and round 3 shipped it without noticing.** Take logs:
            // `(hi-lo)/T < ln hi - ln lo`, and `(hi-lo)/(ln hi - ln lo)` is the
            // logarithmic mean. Verified: 0 disagreements over 5000 seeds.
            //
            // So the "constant-free" form was not constant-free — it carries the
            // constant **1**, on the dimensionless group `E_typical/T`, reached
            // by algebraic coincidence rather than by choice. Worse, it is
            // *scale*-sensitive rather than *spread*-sensitive, so it merges two
            // different claims. Measured: a completely flat matrix — every bond
            // identical, rate spread exactly 1.0, the maximally dead universe —
            // **fails** it, reporting the defect as repaired and instructing the
            // reader to delete the tripwire.
            //
            // Written in the reduced form because the `hi/lo` spelling hides its
            // own threshold and reads as though it were about spread. The
            // `hi > lo` guard splits the two failure modes the old form merged.
            //
            // Honest about sensitivity, which round 3 also got wrong: this fires
            // once `E/T` rises ~5.8x, against ~2.6x for the fitted `< 1.5` it
            // replaced — constant-free bought a **2.3x less sensitive** test. It
            // is kept because it says something true about the *regime* rather
            // than about a sample, not because it is stronger. The rate spread
            // itself (0.242-0.466 of the energy spread over 500 universes) is
            // recorded here and asserted nowhere: it is a magnitude, and a
            // magnitude needs a scale.
            assert!(
                hi > lo,
                "seed {seed}: every bondable single bond has the same energy — the \
                 capacity map has collapsed. That is a bond-matrix defect, NOT an \
                 E/T one, and the `< hi/lo` form used to report it as the repair"
            );
            let logmean =
                (hi - lo) / (borbax_units::det_math::ln(hi) - borbax_units::det_math::ln(lo));
            assert!(
                logmean < t,
                "seed {seed}: typical bondable single-bond energy {logmean} has risen \
                 above mid-temperature {t} — E/T has left the linear regime. Delete \
                 this test and assert the target band in Task 14"
            );
        }
    }

    /// **The guard round 4 deleted, restored as an exact identity.**
    ///
    /// `perturbing_one_elements_binding_energy_moves_every_bond_it_makes` was
    /// removed when capacity stopped reading `energy_per_unit`, on the reasoning
    /// that its premise had gone. The premise had; **the property had not**, and
    /// a review measured the cost: gap-adaptive dilution of the capacity map —
    /// each interior element moved to a drawn point strictly between its
    /// neighbours' midpoints, ordering preserved exactly — passes **53 of 54**
    /// tests here. The one failure is the golden digest, which fires for any
    /// change and so cannot distinguish an injected drawn component from a
    /// deliberate physics edit. That distinction is what this file exists to
    /// provide.
    ///
    /// The identity needs no knowledge of `base`: since `E(a,a) = base * c_a`,
    /// any two elements satisfy `E(a,a) * c_b == E(b,b) * c_a`. Ordering is
    /// preserved by *any* monotone dilution; only a ratio pins the map.
    /// Measured worst relative violation: **3.85e-16** here against **4.91e-1**
    /// under dilution — twelve orders of separation, so the gate carries no
    /// fitted constant.
    #[test]
    fn capacity_enters_bond_energy_as_an_exact_ratio() {
        for seed in 0..16 {
            let (table, m) = universe(seed);
            let crate::element::ShellLaw::V1 { k } = table.pattern().law else {
                unreachable!("this test targets V1")
            };
            let pack = crate::packing::PackingConsts::new(k);
            let cells: Vec<(ElementId, f64, f64)> = table
                .iter()
                .map(|(id, e)| {
                    (
                        id,
                        contact_density(pack, e.units),
                        m.energy(id, id, BondOrder::SINGLE).get(),
                    )
                })
                .collect();
            let mut worst = 0.0_f64;
            for &(_, ca, ea) in &cells {
                for &(_, cb, eb) in &cells {
                    let (l, r) = (ea * cb, eb * ca);
                    let scale = if l.abs() > r.abs() { l.abs() } else { r.abs() };
                    if scale > 0.0 {
                        let rel = (l - r).abs() / scale;
                        if rel > worst {
                            worst = rel;
                        }
                    }
                }
            }
            assert!(
                worst < 1e-12,
                "seed {seed}: worst |E(a,a)*c_b - E(b,b)*c_a| relative residual is \
                 {worst} — something other than contact density is reaching the \
                 capacity scale, which ordering alone cannot see"
            );
        }
    }

    /// **The discriminator for the scale half of the derivation.**
    ///
    /// This replaces the *scale* half of
    /// `perturbing_one_elements_binding_energy_moves_every_bond_it_makes`, whose
    /// premise the round-4 change removed: capacity no longer reads
    /// `energy_per_unit`. The *dilution* half it also carried is restored above,
    /// in `capacity_enters_bond_energy_as_an_exact_ratio` — deleting it outright
    /// was wrong and a review measured the cost.
    ///
    /// What matters now is the property a review measured absent: **scale the
    /// table's energy by alpha and every bond energy must scale by alpha.**
    /// Before `base` was derived from `eps`, that measured 1.0 for any alpha —
    /// the matrix carried the table's ordering and discarded its magnitude,
    /// while magnitude is what enters `exp(-E/T)`. Widening the old `base`
    /// interval would have left it at 1.0, which is what separates this fix
    /// from the plausible one.
    #[test]
    fn scaling_the_tables_energy_scales_every_bond_energy() {
        for seed in 0..12 {
            let table = generate_elements(seed, PhysicsVersion::CURRENT);
            let alpha = 3.0;
            let mut scaled = table.pattern().clone();
            scaled.eps = scaled.eps * alpha;
            let scaled_table = PeriodicTable::new(
                scaled,
                table.iter().map(|(_, e)| e.clone()).collect(),
                table.physics(),
            );

            let mut r0 = Stream::new(seed, Domain::Universe, 2);
            let mut r1 = Stream::new(seed, Domain::Universe, 2);
            let rung = Rung::draw(seed);
            let (m0, m1) = (
                BondEnergyMatrix::generate(&table, seed, rung, &mut r0),
                BondEnergyMatrix::generate(&scaled_table, seed, rung, &mut r1),
            );
            for &a in ids(&table).iter().take(24) {
                for &b in ids(&table).iter().take(24) {
                    let (e0, e1) = (
                        m0.energy(a, b, BondOrder::SINGLE).get(),
                        m1.energy(a, b, BondOrder::SINGLE).get(),
                    );
                    if e0 == 0.0 {
                        continue;
                    }
                    assert!(
                        (e1 / e0 - alpha).abs() < 1e-12,
                        "seed {seed}: {a:?}-{b:?} scaled by {}, not {alpha} — the bond \
                         scale has come unmoored from the table's own energy",
                        e1 / e0
                    );
                }
            }
        }
    }

    /// Pins two things that must not be confused: an out-of-table id must not
    /// read a *neighbouring element's* energy (which unguarded indexing would
    /// do, since `ia * n + ib` stays in range for `ib == n`), and the value it
    /// returns instead is [`Quanta::ZERO`], which under §9.4's
    /// `k = A·exp(−E_bond/T)` is the **fastest-cleaving** bond rather than an
    /// inert one.
    ///
    /// **The second half is a defect this test records, not a guarantee it
    /// blesses.** It is asserted so the value cannot drift unnoticed before
    /// Task 14 removes the sentinel; see [`BondEnergyMatrix::energy`]. A passing
    /// test is a published guarantee, so this one says in its own name which
    /// half is which.
    #[test]
    fn the_zero_sentinel_is_the_fastest_cleaving_value() {
        let (table, m) = universe(9);
        // The first slot *past* the table, not an arbitrary large id: `n` is
        // exactly the value that makes `a * n + b` alias into a different
        // element's cell, so this is the id an unguarded lookup answers wrongly
        // rather than out of range. The fallback cannot fire for a 60..=120
        // table and is still out of range if it ever does.
        // **`len() + 1`, not `len()`.** At exactly `n` the alias `a*n + b` wraps to
        // cell `(a+1, 0)` — column 0 is the *monomer's* column, which is all
        // zeros, so the assertion passed no matter which probe element was
        // chosen. Measured: with `n`, deleting half the bounds check left all 57
        // tests green. `n + 1` aliases to `(a+1, 1)`, a real nonzero cell.
        let past_the_end = ElementId::from_index(table.len() + 1).unwrap_or(ElementId(u8::MAX));
        // **Not `ids(..).first()`** — that is element 0, the monomer, whose
        // whole row is `Quanta::ZERO` by construction. The cell an unguarded
        // read aliases into is *also* zero, so the assertion could not tell "the
        // bound check fired" from "we read a neighbouring element's energy".
        // Measured: deleting the `ib >= self.n` half of the guard left all 57
        // tests green. A probe with nonzero capacity makes it discriminate.
        let first = contacting_ids(&table)
            .first()
            .copied()
            .unwrap_or(ElementId::ZERO);
        assert_eq!(
            m.energy(past_the_end, first, BondOrder::SINGLE),
            Quanta::ZERO
        );
        assert_eq!(
            m.energy(first, past_the_end, BondOrder::SINGLE),
            Quanta::ZERO
        );
    }

    #[test]
    fn bond_order_rejects_what_it_cannot_represent() {
        assert_eq!(BondOrder::try_from(1), Ok(BondOrder::SINGLE));
        assert_eq!(u8::from(BondOrder::SINGLE), 1);
        // The ceiling is the valence ceiling — measured `{4, 5, 6}` across
        // universes and never 3, which is what the old undrawn cap asserted.
        for n in 1..=BondOrder::MAX {
            assert_eq!(BondOrder::try_from(n).map(u8::from), Ok(n));
        }
        assert!(BondOrder::try_from(0).is_err());
        assert!(BondOrder::try_from(BondOrder::MAX + 1).is_err());
    }

    /// **`PLANES` is an array length for a downstream crate, so an off-by-one
    /// here is an out-of-bounds there.**
    ///
    /// Task 6's draft sized its per-order adjacency at `ALL.len()` and indexed
    /// it at `order - 1`. Those agree only while `ALL` is dense — which is
    /// asserted by a `const _` in *this* crate that `borbax-molecule` cannot see
    /// or cite. This test states the property the consumer actually needs:
    /// every representable order has a distinct plane, and every plane in
    /// `0..PLANES` belongs to one. A sparse `ALL` fails it here rather than
    /// panicking there.
    #[test]
    fn every_order_has_its_own_plane_and_they_fill_the_array() {
        let planes: std::collections::BTreeSet<usize> =
            BondOrder::ALL.iter().map(|o| o.plane_index()).collect();
        assert_eq!(
            planes.len(),
            BondOrder::ALL.len(),
            "two orders share a plane: {planes:?}"
        );
        assert!(
            planes.iter().all(|p| *p < BondOrder::PLANES),
            "an order indexes past PLANES ({}): {planes:?}",
            BondOrder::PLANES
        );
        // **The half `ALL` cannot state.** `new` bounds on `1..=MAX` and never
        // consults `ALL`'s membership, so under a sparse `ALL` a `BondOrder`
        // can exist that `ALL` does not list — `new(4)` would still return
        // `Some`. That is the value an unvalidated `u8` actually reaches, and
        // it is a wider gap than the plane arithmetic. Assert over everything
        // `new` admits, not everything `ALL` lists.
        for n in 0..=u8::MAX {
            if let Some(o) = BondOrder::new(n) {
                assert!(
                    o.plane_index() < BondOrder::PLANES,
                    "order {n} indexes plane {} past PLANES {}",
                    o.plane_index(),
                    BondOrder::PLANES
                );
            }
        }
    }

    /// Decision 8's "one shared computation" claim, for bond capacity
    /// specifically: `capacity(el) = capacity_scale + valence(el)` must be
    /// an exact function of `Element::valence` (itself gap-derived) and
    /// nothing else -- checked by replaying the production draw and
    /// recomputing the whole V2 formula bit-for-bit, the same discipline
    /// `both_bond_constants_are_this_universes_own_draws` uses for V1.
    #[test]
    fn v2_bond_energy_is_exactly_derived_from_valence_and_affinity() {
        for seed in [0u64, 1, 5, 21, 42] {
            let (table, bonds) = universe_v2(seed);
            let rung = Rung::draw(seed);
            let pattern = table.pattern();

            let scale = draw_symmetric(seed, rung, MigratedConstant::Scale, 56.0);
            let base = pattern.eps * scale;
            let capacity_scale = draw_symmetric(
                seed,
                rung,
                MigratedConstant::CapacityScale,
                BondEnergyMatrix::BASE_CAPACITY_SCALE,
            );
            let ionic_scale = draw_symmetric(
                seed,
                rung,
                MigratedConstant::IonicScale,
                BondEnergyMatrix::BASE_IONIC_SCALE,
            );
            let ionic_base = pattern.eps * ionic_scale;

            let ids = ids(&table);
            for &a in ids.iter().step_by(7) {
                for &b in ids.iter().step_by(11) {
                    let ea = table
                        .get(a)
                        .unwrap_or_else(|| unreachable!("id from this table"));
                    let eb = table
                        .get(b)
                        .unwrap_or_else(|| unreachable!("id from this table"));
                    let ca = capacity_scale + f64::from(ea.valence);
                    let cb = capacity_scale + f64::from(eb.valence);
                    let want = base * (ca * cb).sqrt() - ionic_base * (ea.affinity * eb.affinity);
                    let got = bonds.energy(a, b, BondOrder::SINGLE);
                    assert!(
                        (got.get() - want.get()).abs() < 1e-9,
                        "seed {seed} {a:?}-{b:?}: bond energy {} != replayed formula {} -- \
                         capacity/ionic terms are not exactly derived from valence/affinity",
                        got.get(),
                        want.get()
                    );
                }
            }
        }
    }

    /// The worst-case positivity proof `MigratedConstant::IonicScale`'s own
    /// doc states, checked directly rather than trusted: the ionic term
    /// can only drive a bond negative when both elements' affinities are
    /// near-saturated same-signed *and* `Scale`/`CapacityScale` are both at
    /// their perturbed floor simultaneously. Sweep every reachable rung
    /// (not just the seed's own drawn rung) with `Universe::generate_under`
    /// unavailable at arbitrary rungs, so this drives `BondEnergyMatrix::generate`
    /// directly the way `gamma_stays_a_diminishing_series_at_every_reachable_rung`
    /// drives `draw_symmetric` directly.
    #[test]
    fn bond_energy_stays_finite_and_non_negative_at_every_reachable_rung() {
        for seed in [0u64, 1, 7, 42, 1000] {
            let table = generate_elements(seed, PhysicsVersion::V2);
            for r in 0..=Rung::MAX {
                let rung = Rung::new(r).unwrap_or_else(|| unreachable!("r <= Rung::MAX"));
                let mut rng = Stream::new(seed, Domain::Universe, 2);
                let bonds = BondEnergyMatrix::generate(&table, seed, rung, &mut rng);
                for (a, _) in table.iter() {
                    for (b, _) in table.iter() {
                        let v = bonds.energy(a, b, BondOrder::SINGLE);
                        assert!(
                            v.get().is_finite(),
                            "seed {seed} rung {r} {a:?}-{b:?}: bond energy should be finite"
                        );
                        assert!(
                            v.get() >= 0.0,
                            "seed {seed} rung {r} {a:?}-{b:?}: bond energy {} went negative -- \
                             the ionic term's worst-case margin was not what \
                             MigratedConstant::IonicScale's doc claims",
                            v.get()
                        );
                    }
                }
            }
        }
    }

    /// Power iteration for a symmetric matrix's dominant eigenpair --
    /// test-only, no external linear-algebra dependency: `M v` repeated
    /// and renormalised converges to the eigenvector of largest-magnitude
    /// eigenvalue for a symmetric matrix, which is exactly the rank-1
    /// approximation `bond_matrix_rank2_residual_correlates_with_affinity_difference_squared`
    /// needs to subtract off.
    #[expect(
        clippy::indexing_slicing,
        reason = "every index is < n by loop construction: m is n*n, v/next/mv are all built at \
                  length n"
    )]
    fn dominant_eigenpair(m: &[f64], n: usize) -> (f64, Vec<f64>) {
        let mut v = vec![1.0; n];
        for _ in 0..300 {
            let mut next = vec![0.0; n];
            for (i, slot) in next.iter_mut().enumerate() {
                *slot = (0..n).map(|j| m[i * n + j] * v[j]).sum();
            }
            let norm = next.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in &mut next {
                *x /= norm;
            }
            v = next;
        }
        let mv: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| m[i * n + j] * v[j]).sum())
            .collect();
        let lambda: f64 = (0..n).map(|i| v[i] * mv[i]).sum();
        (lambda, v)
    }

    /// Task 26.1 Step 6's headline discriminator: the bond matrix's
    /// rank-1-removed residual must correlate with the affinity
    /// *difference*, checked directly against the constructed
    /// `BondEnergyMatrix` output rather than by recomputing the source
    /// formula (which `v2_bond_energy_is_exactly_derived_from_valence_and_affinity`
    /// already does) -- this is what actually discriminates "the ionic
    /// term is present and doing its job" from "the matrix is still
    /// rank-1", the way the plan's own acceptance criterion demands.
    ///
    /// **`(affinity[a]-affinity[b])^2`, not the raw signed difference.**
    /// `BondEnergyMatrix` is symmetric by construction, so a residual
    /// correlated against an antisymmetric quantity averaged over the
    /// off-diagonal is not the meaningful comparison; the squared
    /// difference is the symmetric analogue and is what a Pauling-style
    /// ionic-excess term is actually shaped like.
    #[expect(
        clippy::indexing_slicing,
        reason = "a, b < n by loop construction; v is built at length n and m at length n*n"
    )]
    #[test]
    fn bond_matrix_rank2_residual_correlates_with_affinity_difference_squared() {
        const CORRELATION_FLOOR: f64 = 0.5;
        for seed in [0u64, 1, 5, 21, 42] {
            let (table, bonds) = universe_v2(seed);
            let n = table.len();
            let ids = ids(&table);

            let mut m: Vec<f64> = Vec::with_capacity(n * n);
            for &a in &ids {
                for &b in &ids {
                    m.push(bonds.energy(a, b, BondOrder::SINGLE).get());
                }
            }
            let (lambda, v) = dominant_eigenpair(&m, n);

            let affinity: Vec<f64> = table.iter().map(|(_, e)| e.affinity).collect();
            let mut residual = Vec::with_capacity(n * (n - 1) / 2);
            let mut diffsq = Vec::with_capacity(n * (n - 1) / 2);
            for a in 0..n {
                for b in (a + 1)..n {
                    let fit = lambda * v[a] * v[b];
                    residual.push(m[a * n + b] - fit);
                    let d = affinity.get(a).unwrap_or_else(|| unreachable!("a < n"))
                        - affinity.get(b).unwrap_or_else(|| unreachable!("b < n"));
                    diffsq.push(d * d);
                }
            }
            let corr = pearson_correlation(&residual, &diffsq);
            assert!(
                corr > CORRELATION_FLOOR,
                "seed {seed}: rank-1-removed residual correlates {corr:.4} with \
                 (affinity_diff)^2 -- below the floor, so the matrix is not \
                 discriminably rank-2"
            );
        }
    }

    /// Pearson correlation, hand-rolled -- see `element.rs`'s own copy of
    /// this helper (same reasoning: a test-only ~10-line helper is a lower
    /// bar than even a dev-dependency).
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
}
