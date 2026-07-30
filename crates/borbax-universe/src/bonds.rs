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
//! Measured on the current code over 500 universes, all bondable pairs:
//!
//! | quantity | was | is |
//! |---|---|---|
//! | single-bond energy | — | `[11.57, 501.76]` Quanta |
//! | cleave-rate ratio, single bond, mid-`T` | 1.043–1.243 | **1.185–3.728** |
//! | same, orders 1..=3, at `temp_min` | max 3.50, >2.0 in 716/2000 | **max 614**, >2.0 in **500/500** |
//! | same, orders 1..=6, at `temp_min` | not representable | **4.0 – 217 758** |
//! | `logmean(E)/T` (the assertion below) | — | `[0.083, 0.571]`, **1.75x** of headroom |
//!
//! **The absolute energy scale rose about fivefold, and that was not stated
//! when it happened.** Capacity was normalised onto `[w, 1]` and is now contact
//! density, which runs toward `z/2` — so bond energies moved from roughly
//! `6.6..42.7` to `23.9..221.0` on seed 0. The new single-bond rate *minimum* is
//! roughly the old *maximum*, and orders reaching 6 multiply the exponent again.
//!
//! **A third consequence rode along and is worth naming separately**: the draw
//! count on `Domain::Universe` index 1 went from four to two, so every constant
//! after the bond draws shifted. Seed 0's `temp_min` moved from 206.09 to
//! 222.68, and its solvent with it. That is not "the bond derivation changed";
//! it is every shared seed naming a different universe, which is exactly what a
//! digest exists to make visible rather than silent.
//!
//! **What that does and does not settle.** §22.6 wants "some linkages durable
//! and others fragile", and a 1.19–3.73x lifetime ratio from composition alone —
//! with bond order supplying five further orders of magnitude — is a different
//! situation from the 1.04–1.24x that prompted the finding. Whether it *suffices*
//! is a chemistry judgement for Task 20's battery, not a claim this file should
//! make. What it does settle is that any sizing done against the old numbers is
//! wrong: an implementer applying Task 14's coupling on top of these while
//! reading `1.043–1.243` would oversize it badly.
//!
//! `the_single_bond_rate_spread_at_mid_temperature_is_narrow` now sits **1.75x**
//! from firing, not the ~5.8x its own comment claimed under the old map. The
//! test that exists to say when the regime changed is close to saying so.
//!
//! **STOP before applying Task 14's planned remedy.** It is written as "draw the
//! energy in units of the universe's own temperature, which collapses the `Ea/T`
//! spread from 4.39x to 1.09x" — sized against the *old* map. Applied to this
//! one it would collapse a cleave-rate spread of **3 to 217 758** down to 1.09x,
//! deleting §9.4's differential persistence in order to fix a problem this
//! commit already fixed. Measured at HEAD: `E/T` reaches **12.52**, and the rate
//! ratio at `temp_min` exceeds 2.0 in **2000 of 2000** universes. Task 14 must
//! re-derive against these numbers before touching anything.
//!
//! The remaining question is a *band*, not a direction: it belongs with the rate
//! law, is written into Task 14 with its discriminator, and into Task 20 for the
//! band. Picking a number here would be the typed-in constant this
//! project keeps finding reported back as an emergent result — and note that
//! deriving the *energy* scale, as this file now does, does not fix the
//! *ratio*: `eps` is still drawn independently of `temp_min`. The remedy Task 14
//! carries is to draw the energy in units of the universe's own temperature,
//! which a review measured collapses the `Ea/T` spread from 4.39x to 1.09x.

use crate::element::{ElementId, PeriodicTable};
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
/// **The per-pair cap is Task 11's**, not this type's: a bond of order *n*
/// consumes *n* slots at each end, so the real bound is
/// `min(valence_a, valence_b)`, and nothing here forms a bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondOrder(u8);

impl BondOrder {
    /// A single contact.
    pub const SINGLE: Self = Self(1);
    /// The largest order any universe's valence ceiling admits.
    pub const MAX: u8 = 6;
    /// [`Self::MAX`] as an array length, so a multiplier table cannot be sized
    /// independently of the bound it must cover.
    const MAX_LEN: usize = 6;

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
    mult: [f64; BondOrder::MAX_LEN],
}

impl OrderScale {
    /// Precomputed at generation, so `powf` never runs on a lookup.
    ///
    /// The old two-multiplier form was a field read; moving to an exponent put a
    /// transcendental in `energy()`, which is the one function a rate path will
    /// call per species. `gamma` is fixed per universe and `MAX` is small, so
    /// the table costs six `powf` calls once.
    fn new(gamma: f64) -> Self {
        let mut mult = [1.0; BondOrder::MAX_LEN];
        for (i, m) in mult.iter_mut().enumerate() {
            // `i` runs 0..MAX_LEN, so the narrowing is lossless; `as` is denied.
            let order = u8::try_from(i).unwrap_or(0).saturating_add(1);
            *m = det_math::powf(f64::from(order), gamma);
        }
        Self { mult }
    }

    /// Total by construction, for every representable order.
    fn of(self, order: BondOrder) -> f64 {
        let i = usize::from(u8::from(order)) - 1;
        self.mult.get(i).copied().unwrap_or(1.0)
    }
}

impl BondEnergyMatrix {
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
    pub(crate) fn generate(table: &PeriodicTable, rng: &mut Stream) -> Self {
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
        let scale = rng.next_f64_range(28.0, 84.0);
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
        let order_scale = OrderScale::new(rng.next_f64_range(0.68, 1.0));

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
        let pack = crate::packing::PackingConsts::new(pattern.k);
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
    /// premise for the `n²` layout, and the premise was false: measured, the
    /// table buys **0%** over recomputing `base * sqrt(c_a * c_b)` (+0.9% in the
    /// intern loop shape, +0.3% in a hot loop, both inside a 0.02–0.08% noise
    /// floor), and both forms are bit-identical over 129 086 cells. The table is
    /// kept anyway — it is the only shape that could later hold a pair rule that
    /// does not factorise — but nobody should re-derive that from a cost story.
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
    /// largest cleave propensity available. Measured on seed 0
    /// (`T = temp_min = 206.09`, single-bond energies `6.60..42.71`): the
    /// sentinel rates **1.03x the fastest** genuine bond and 1.23x the slowest.
    /// An earlier version of this sentence attached `1.23x` to the fastest,
    /// which two review lanes caught independently — real number, wrong
    /// referent, in the one block whose stated job is recording a defect
    /// honestly. The corrected figure sharpens the point rather than weakening
    /// it: a sentinel sitting 3% outside the genuine range *is* the §22.6
    /// finding, and it becomes ~10^2-10^5 once Task 14 raises `E/T`.
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
    /// boundary, makes this lookup total and deletes the question. There is no
    /// consumer today, so nothing is presently wrong — the hazard is entirely
    /// prospective, which is exactly why it is written into Task 14 (reaction rates)
    /// rather than left in this comment. `the_zero_sentinel_is_the_fastest_cleaving_value`
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

    /// The multiplier this universe applies to an order. Test-visible so the
    /// diminishing-returns choice can be asserted rather than only documented.
    #[cfg(test)]
    fn order_multiplier(&self, order: BondOrder) -> f64 {
        self.order_scale.of(order)
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
    use crate::element::{ElementId, PeriodicTable, generate_elements};
    use borbax_rng::{Domain, Stream};

    /// One universe's table and its bond matrix, drawn the way
    /// [`crate::Universe::generate`] draws them.
    fn universe(seed: u64) -> (PeriodicTable, BondEnergyMatrix) {
        let table = generate_elements(seed);
        let mut rng = Stream::new(seed, Domain::Universe, 2);
        let bonds = BondEnergyMatrix::generate(&table, &mut rng);
        (table, bonds)
    }

    fn ids(table: &PeriodicTable) -> Vec<ElementId> {
        table.iter().map(|(id, _)| id).collect()
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
            let table = generate_elements(seed);
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
            let pack = crate::packing::PackingConsts::new(table.pattern().k);
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
            let pack = crate::packing::PackingConsts::new(table.pattern().k);
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
                m.order_multiplier(o)
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
            let ids: Vec<ElementId> = bondable_ids(&table).into_iter().take(16).collect();
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
        for &a in &bondable_ids(&table) {
            for &b in &bondable_ids(&table) {
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
        for &a in &bondable_ids(&table) {
            for &b in &bondable_ids(&table) {
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
            let pack = crate::packing::PackingConsts::new(table.pattern().k);
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
    /// words "most durable bondable bond", which is a *triple* bond — a
    /// different and wider number (all orders at `temp_min` reaches 3.50, and
    /// exceeds 2.0 in 716 of 2000 universes). Sizing Task 14's coupling from the
    /// narrow slice while reading the wide sentence would oversize it ~2.4x.
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
            let pack = crate::packing::PackingConsts::new(table.pattern().k);
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
            let table = generate_elements(seed);
            let alpha = 3.0;
            let mut scaled = table.pattern().clone();
            scaled.eps = scaled.eps * alpha;
            let scaled_table =
                PeriodicTable::new(scaled, table.iter().map(|(_, e)| e.clone()).collect());

            let mut r0 = Stream::new(seed, Domain::Universe, 2);
            let mut r1 = Stream::new(seed, Domain::Universe, 2);
            let (m0, m1) = (
                BondEnergyMatrix::generate(&table, &mut r0),
                BondEnergyMatrix::generate(&scaled_table, &mut r1),
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
        let past_the_end = ElementId::from_index(table.len()).unwrap_or(ElementId(u8::MAX));
        let first = ids(&table).first().copied().unwrap_or(ElementId::ZERO);
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
}
