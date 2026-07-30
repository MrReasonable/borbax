//! Generated bond energies (spec §7.1).
//!
//! **Bond energy is a function of the two elements' per-unit binding energies,
//! and that is this file's entire content.** §7.1's `bond_energies` was the last
//! substantive property of a universe still without a cause: Task 4 derived
//! `period`, `group`, `valence`, `radius`, `affinity` and `mass` from how a
//! cluster packs, and left this one to be drawn. Drawing it reproduces exactly
//! the defect that made the predecessor's `peak` unfalsifiable — "these
//! elements bind well" encoded twice, in incommensurable units, so neither
//! encoding can check the other.
//!
//! So: an element whose units bind strongly to each other binds strongly to
//! other elements, because a bond is the same kind of contact as the ones
//! holding the cluster together. That is principle 1 — one mechanism reused —
//! rather than a second axis invented for bonding.
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
//! character with an R² *ratio of 0.30–0.90* — **below chance**, so partitioning
//! by group is worse than an arbitrary partition of the same shape. `period`
//! scores 19.8–41.6× above the null.
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
//! That is deliberate rather than a shortfall. Principle 2 says shape is the
//! only mechanism, so partner selectivity belongs to §8.3's signature
//! complementarity; a per-pair covalent energy axis would be the second
//! mechanism principle 2 exists to refuse. It is stated here because the type is
//! called a matrix, and a reader building §9.1's energetics on top will
//! otherwise assume pair structure exists. If selectivity is ever wanted here,
//! the question to answer first is what the shape model failed to do.
//!
//! # What this file gets wrong, recorded rather than hidden
//!
//! The energy scale drawn here and the temperature scale drawn in
//! [`crate::UniverseConsts`] are independent, and §9.1/§9.5 make cleavage
//! `k = A·exp(−E_bond/T)`. With `base ∈ [30, 90]` against `T ∈ [180, 760]`,
//! `E/T` never leaves `exp`'s linear regime: the cleave-rate ratio between the
//! most fragile and most durable bondable bond is **1.04×–1.48×**, and about
//! **1.05×–1.20×** over the abundance-weighted population a beaker actually
//! contains. A molecule's spontaneous lifespan is therefore set by how many
//! bonds it has, with composition worth the equivalent of one to five extra
//! bonds — the same shape as the recorded radiogenic finding, where a
//! composition channel is swamped by an atom-count channel.
//!
//! §22.6's plan is to locate the decay band globally and then introduce
//! per-bond rates "as a redistribution rather than a fresh search". At this
//! scale that redistribution is ~1% of the range the global constant already
//! searches. **The energy spread is adequate and the rate spread is not**, and
//! no test in this crate could previously tell those apart — a review showed
//! that multiplying `base` by 0.001 leaves every spread assertion passing in a
//! universe where differential persistence is identically zero.
//!
//! The fix sets the dimensionless group `E/T`, not `E`, and it belongs with the
//! rate law: it is written into Task 14 with its discriminator, and into Task 20
//! for the band. Picking a number here would be the typed-in constant this
//! project keeps finding reported back as an emergent result.

use crate::element::{ElementId, PeriodicTable};
use borbax_rng::Stream;
use borbax_units::Quanta;

/// How many bonds join one pair of atoms.
///
/// **An enum rather than a `u8`, because the `u8` version has no honest
/// answer for zero.** The plan's signature took `order: u8` and clamped to
/// `1..=3` under the comment "a malformed molecule should be inert, not a
/// crash" — but clamping `0` up to `1` is the opposite of inert: it hands a
/// non-bond the energy of a single bond. Making the state unrepresentable
/// removes the choice, and pushes the one real question — what a caller does
/// with an out-of-range number — to [`TryFrom`], where it is visible.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondOrder {
    /// One bond.
    Single = 1,
    /// Two bonds.
    Double = 2,
    /// Three bonds.
    Triple = 3,
}

/// The error [`BondOrder::try_from`] returns for a count it cannot represent.
///
/// `thiserror` rather than a hand-rolled `Display`/`Error` pair: it is already
/// a workspace dependency and already compiled into this crate's tree via
/// `borbax-units`, so it costs zero extra crates and zero build time. It is a
/// proc macro emitting a `Display` impl and cannot reach a result-affecting
/// path, so the "runtime dependency in a result-affecting path is expensive"
/// rule does not bite.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("bond order {0} is not in 1..=3")]
pub struct NotABondOrder(pub u8);

impl TryFrom<u8> for BondOrder {
    type Error = NotABondOrder;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            1 => Ok(Self::Single),
            2 => Ok(Self::Double),
            3 => Ok(Self::Triple),
            other => Err(NotABondOrder(other)),
        }
    }
}

/// Closes the round-trip the variants' `= 1, 2, 3` discriminants promise.
///
/// Without this the only route back is `order as u8`, which `as_conversions`
/// warns on — so the discriminants would be decoration inviting exactly the
/// conversion the lint forbids, at the point `Mol12` needs to serialise a bond
/// order.
impl From<BondOrder> for u8 {
    fn from(order: BondOrder) -> Self {
        match order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        }
    }
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
    /// Multipliers for the two orders above single.
    order_scale: OrderScale,
}

/// How much stronger a double and a triple bond are than a single one.
///
/// **Named fields rather than `[f64; 3]`, because the array undid [`BondOrder`]
/// one line after it was established.** Indexing an array by a converted enum
/// reintroduces the integer the enum exists to remove, and the miss needs an
/// `unwrap_or(1.0)` fallback that would silently return the *single*-bond
/// multiplier — an unreachable branch whose value is wrong if it is ever
/// reached. Matching on the variant makes `Single => 1.0` a fact rather than an
/// array slot that happens to hold `1.0`, and makes a fourth order fail to
/// compile in one place instead of two.
#[derive(Debug, Clone, Copy, PartialEq)]
struct OrderScale {
    double: f64,
    triple: f64,
}

impl OrderScale {
    /// Total by construction — there is no order without a multiplier.
    const fn of(self, order: BondOrder) -> f64 {
        match order {
            BondOrder::Single => 1.0,
            BondOrder::Double => self.double,
            BondOrder::Triple => self.triple,
        }
    }
}

impl BondEnergyMatrix {
    /// Derive the matrix from a table's binding energies.
    ///
    /// Three constants are drawn — the energy scale, the weakest-bond fraction,
    /// and the two higher-order multipliers — because they are this universe's
    /// *character*, in the same sense as `eps` and `sigma` in
    /// [`crate::element::generate_elements`]. What is deliberately **not** drawn
    /// is anything per-pair: the plan's draft added `next_normal() * 0.08` to
    /// every cell, which is an independently drawn matrix wearing a derivation,
    /// in the exact proportion of the noise. Every difference between two cells
    /// here traces to the two elements' binding energies.
    ///
    /// **`pub(crate)` for the same reason [`PeriodicTable::new`] is.** A public
    /// constructor lets a caller build a matrix from table A and query it with
    /// ids from table B — in range, plausible and silent, which is the failure
    /// this crate keeps closing. Universes come from [`crate::Universe`], and
    /// the matrix comes with one.
    pub(crate) fn generate(table: &PeriodicTable, rng: &mut Stream) -> Self {
        let n = table.len();

        // Drawn before anything reads the table, so the draw order is a
        // property of this function rather than of the table handed to it.
        let base = rng.next_f64_range(30.0, 90.0);
        // What the *weakest* element in the table brings to a bond, as a
        // fraction of what the strongest brings.
        //
        // **This is a floor on the capacity scale, not on the finished energy.**
        // The first draft normalised binding energy onto `[0, 1]` and put the
        // floor on the bond instead. That gives the weakest element capacity
        // exactly zero — and zero is absorbing under a geometric mean, so every
        // bond to it scored identically no matter the partner.
        //
        // **The mechanism is real; the target named in the first version of this
        // comment was wrong, and the correction is worth keeping.** It said the
        // weakest element "is reliably the monomer, the feedstock, the most
        // abundant species", so the commonest bond in the chemistry would have
        // had no structure. The first clause is right — `contacts_upto(1) = 0`,
        // so the monomer is the minimum in 200 of 200 drawn universes — and the
        // conclusion does not follow, because the monomer has `outer == 0`,
        // hence no frontier, hence `valence == 0`. **It forms no bonds at all**,
        // so the row this floor was described as rescuing is never priced.
        //
        // What the floor still buys is the rare case where the minimum is *not*
        // the monomer: at large `sigma` and small `eps` the strain term
        // overtakes the contact term and the heaviest element becomes the
        // minimum, and that element can have a frontier. Keeping the floor is
        // right; believing it protects the feedstock is not.
        //
        // Below 1/3 makes the nominal whole-matrix ratio exactly
        // `1 / weakest_share` — see `bondable_energies_span_a_recorded_band`
        // for why that nominal figure is not the quantity §22.6 is about.
        let weakest_share = rng.next_f64_range(0.12, 0.28);
        // 2.2 > 2.1, so triple strictly exceeds double in every universe. The
        // ranges are adjacent on purpose; overlapping them would make
        // `higher_orders_are_stronger` a property of the draw.
        let order_scale = OrderScale {
            double: rng.next_f64_range(1.6, 2.1),
            triple: rng.next_f64_range(2.2, 3.0),
        };

        // The binding-energy range of this table. Folded in table order, and
        // `f64::min`/`max` are disallowed (§13.1 — they are non-deterministic
        // for ±0.0), so the comparisons are written out.
        let (lo, hi) = table.iter().fold((f64::MAX, f64::MIN), |(l, h), (_, el)| {
            let v = el.energy_per_unit.get();
            (if v < l { v } else { l }, if v > h { v } else { h })
        });
        let span = hi - lo;

        // What each element brings to a bond, on a `[weakest_share, 1]` scale.
        // Strictly positive at both ends, which is what keeps the weakest
        // element from annihilating its partners.
        //
        // A degenerate table — every element identically bound — would divide
        // by zero, so it maps to the top of the scale: every bond then costs
        // `base`, which is the honest answer for a universe with no binding
        // structure, and `energies_have_useful_spread` fails on it rather than
        // a fudge hiding it. Unreachable while the table draws 60..=120
        // elements, since `energy_per_unit` is strictly monotone in `units`
        // over the first shell.
        // **The association below is load-bearing and pinned (§13.4).** It reads
        // `(1-w) * (x-lo) / span`, which associates as `((1-w) * (x-lo)) / span`.
        // Rewriting it as `(1-w) * ((x-lo) / span)` is a reassociation a reviewer
        // could ask for as clarification, it is mathematically identical, and it
        // moves bond energies by an ulp: measured, the digest over 129 086 cells
        // goes 9ef10ca10dd88b5c -> 34f1ac1ed3a2e377 while all 50 tests pass.
        // `the_assembled_universe_digest_is_pinned` is what now catches it.
        // Nobody tidies this later.
        let capacity: Vec<f64> = table
            .iter()
            .map(|(_, el)| {
                if span > 0.0 {
                    weakest_share + (1.0 - weakest_share) * (el.energy_per_unit.get() - lo) / span
                } else {
                    1.0
                }
            })
            .collect();

        let mut e = vec![Quanta::ZERO; n * n];
        for a in 0..n {
            // Hoisted out of the inner loop: `a` is loop-invariant, and leaving
            // the lookup inside made LLVM version-clone the whole loop on an
            // unreachable fallback rather than lift the load. Legibility, not
            // speed — the measured effect is nil at 4.4 us.
            let ca = capacity.get(a).copied().unwrap_or(weakest_share);
            for b in a..n {
                // Geometric mean, so the pair is symmetric bit-for-bit
                // (multiplication commutes exactly) and a weak partner drags the
                // bond down rather than being rescued by a strong one the way a
                // sum would allow — a bond needs both ends to hold. `sqrt` is
                // exactly specified by IEEE-754 and stays native (§13.1).
                // Both fallbacks are unreachable (`capacity.len() == n` and
                // `a, b < n`) — but they default to `weakest_share`, not to
                // `1.0`. An earlier version defaulted to the maximum, so an
                // element that fell out of the table would have been given the
                // *strongest* bonds in the universe, thirty lines from a
                // comment arguing the conservative default at the other
                // unreachable site. Two unreachable branches in one file
                // disagreeing about which direction is safe is the finding,
                // regardless of whether either can fire.
                let cb = capacity.get(b).copied().unwrap_or(weakest_share);
                let v = Quanta(base * (ca * cb).sqrt());
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
    /// largest cleave propensity available. Measured on seed 0 (`T = 206.09`,
    /// real energies `6.60..42.71`): the zero sentinel rates `1.23x` the
    /// fastest genuine bond.
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
    /// rather than left in this comment. `the_zero_sentinel_is_known_wrong_signed`
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
    use crate::element::{ElementId, PeriodicTable, generate_elements};
    use borbax_rng::{Domain, Stream};

    /// One universe's table and its bond matrix, drawn the way
    /// [`crate::Universe::generate`] draws them.
    fn universe(seed: u64) -> (PeriodicTable, BondEnergyMatrix) {
        let table = generate_elements(seed);
        let mut rng = Stream::new(seed, Domain::Universe, 1);
        let bonds = BondEnergyMatrix::generate(&table, &mut rng);
        (table, bonds)
    }

    fn ids(table: &PeriodicTable) -> Vec<ElementId> {
        table.iter().map(|(id, _)| id).collect()
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
                    m.energy(a, b, BondOrder::Single).get().to_bits(),
                    m.energy(b, a, BondOrder::Single).get().to_bits(),
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

    /// The spread is structural, not lucky: the extremes are the weakest and
    /// strongest elements bonded to themselves, so the ratio over the whole
    /// matrix is exactly `1 / weakest_share` — which is why drawing that
    /// constant below 1/3 is what `energies_have_useful_spread` really asserts.
    #[test]
    fn the_spread_is_exactly_the_reciprocal_of_the_weakest_share() {
        for seed in 0..16 {
            let (table, m) = universe(seed);
            let all = ids(&table);
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for &a in &all {
                for &b in &all {
                    let v = m.energy(a, b, BondOrder::Single).get();
                    if v < lo {
                        lo = v;
                    }
                    if v > hi {
                        hi = v;
                    }
                }
            }
            // `weakest_share` is the second draw of the bond stream.
            let mut rng = Stream::new(seed, Domain::Universe, 1);
            let _base = rng.next_f64_range(30.0, 90.0);
            let weakest_share = rng.next_f64_range(0.12, 0.28);
            let ratio = hi / lo;
            // 1e-12, not 1e-9. Measured worst deviation over 400 universes is
            // 1.776e-15, so the old gate sat six orders above the signal — wide
            // enough to admit a one-ulp reassociation of the capacity
            // expression, which is exactly the drift a golden must catch.
            assert!(
                (ratio - 1.0 / weakest_share).abs() < 1e-12,
                "seed {seed}: spread {ratio} is not 1/{weakest_share}"
            );
            assert!(
                ratio > 3.0,
                "seed {seed}: spread {ratio} fell to the spread floor"
            );
        }
    }

    /// **The other half of the derivation claim, and the half six mutation
    /// probes missed.**
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
    fn capacity_is_a_function_of_binding_energy_and_of_nothing_else() {
        for seed in 0..24 {
            let (table, m) = universe(seed);
            let all: Vec<(ElementId, f64)> = table
                .iter()
                .map(|(id, e)| (id, e.energy_per_unit.get()))
                .collect();
            let mut inversions = 0_usize;
            for &(ia, ea) in &all {
                for &(ib, eb) in &all {
                    let self_a = m.energy(ia, ia, BondOrder::Single).get();
                    let self_b = m.energy(ib, ib, BondOrder::Single).get();
                    // Strict on both sides, so ties must correspond to ties.
                    if (ea < eb) != (self_a < self_b) {
                        inversions += 1;
                    }
                }
            }
            assert_eq!(
                inversions, 0,
                "seed {seed}: {inversions} pairs order differently by binding \
                 energy and by self-bond energy — something other than \
                 `energy_per_unit` is reaching the capacity scale"
            );
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
    /// `c` — so this pins the combiner without pinning the normalisation.
    #[test]
    fn the_matrix_is_exactly_rank_one() {
        for seed in 0..8 {
            let (table, m) = universe(seed);
            let ids: Vec<ElementId> = ids(&table).into_iter().take(16).collect();
            let mut worst = 0.0_f64;
            for &a in &ids {
                for &b in &ids {
                    for &c in &ids {
                        for &d in &ids {
                            let (ab, cd) = (
                                m.energy(a, b, BondOrder::Single).get(),
                                m.energy(c, d, BondOrder::Single).get(),
                            );
                            let (ad, cb) = (
                                m.energy(a, d, BondOrder::Single).get(),
                                m.energy(c, b, BondOrder::Single).get(),
                            );
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
                "seed {seed}: worst relative 2x2 minor {worst} — the matrix is \
                 not rank 1, so the combiner is no longer a geometric mean"
            );
        }
    }

    /// Over every pair, not one hand-picked pair: the order multipliers are
    /// universe-wide, so a single pair cannot distinguish "orders increase"
    /// from "these two elements happen to increase".
    #[test]
    fn higher_orders_are_stronger() {
        let (table, m) = universe(4);
        for &a in &ids(&table) {
            for &b in &ids(&table) {
                let (s, d, t) = (
                    m.energy(a, b, BondOrder::Single),
                    m.energy(a, b, BondOrder::Double),
                    m.energy(a, b, BondOrder::Triple),
                );
                assert!(d > s && t > d, "{a:?}-{b:?}: {s:?} {d:?} {t:?}");
            }
        }
    }

    #[test]
    fn energies_are_positive_and_finite() {
        let (table, m) = universe(9);
        for &a in &ids(&table) {
            for &b in &ids(&table) {
                let e = m.energy(a, b, BondOrder::Single);
                assert!(e.get() > 0.0 && e.is_finite(), "{a:?}-{b:?} = {e:?}");
            }
        }
    }

    /// Bond-energy spread over the elements a molecule can actually contain.
    ///
    /// **This replaces `energies_have_useful_spread`, which asserted
    /// `hi/lo > 3.0` over the whole matrix and cited §22.6 — and was measuring
    /// cells no molecule can occupy.** `valence == 0` means no frontier and so
    /// no bonding slots. Measured over 200 universes: the *weakest* element is
    /// the monomer and has valence 0 in **200 of 200**; the *strongest* sits on
    /// a binding-peak closure and has valence 0 in **170 of 200**. So both
    /// endpoints of the asserted range are unreachable, and the old threshold
    /// failed on the bondable submatrix in **86 of 200** universes while passing
    /// on all 20 seeds it ran.
    ///
    /// The band below is the measured bondable range, recorded so drift is
    /// visible. It is deliberately *not* dressed up as satisfying §22.6 — see
    /// the module header for why the quantity §22.6 is about is a **rate**
    /// ratio, which this commit's scales leave near 1.
    #[test]
    fn bondable_energies_span_a_recorded_band() {
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
                    let v = m.energy(a, b, BondOrder::Single).get();
                    if v < lo {
                        lo = v;
                    }
                    if v > hi {
                        hi = v;
                    }
                }
            }
            let ratio = hi / lo;
            assert!(
                (2.4..=4.6).contains(&ratio),
                "seed {seed}: bondable spread {ratio} left the recorded band \
                 [2.4, 4.6] — requote the module header and Task 14's note"
            );
        }
    }

    /// **Records the defect four review lanes converged on, so it cannot be
    /// mistaken for a satisfied requirement.**
    ///
    /// §22.6 wants "some linkages durable and others fragile", and §9.4 calls
    /// differential persistence the thing selection acts on. Both are claims
    /// about *rates*, and §9.1/§9.5 make the rate `k = A·exp(−E_bond/T)`. This
    /// commit mints both scales — `base ∈ [30, 90]` Quanta and `T ∈ [180, 760]`
    /// Thermal — so `E/T` lands near zero and `exp` never leaves its linear
    /// regime. The energy spread is fine; the *rate* spread is not.
    ///
    /// The assertion is deliberately an upper bound on a defect rather than a
    /// lower bound on a virtue: it fails if the ratio grows, which is what
    /// Task 14 fixing the scale will do. That is the intended way for this test
    /// to die — deliberately, in the commit that repairs the physics, rather
    /// than silently.
    #[test]
    fn the_cleave_rate_spread_is_currently_far_too_narrow() {
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
                    let e = u.bonds.energy(a, b, BondOrder::Single).get();
                    if e < lo {
                        lo = e;
                    }
                    if e > hi {
                        hi = e;
                    }
                }
            }
            // Fastest cleave is the weakest bond. Ratio of rates, not energies.
            let rate_ratio =
                borbax_units::det_math::exp(-lo / t) / borbax_units::det_math::exp(-hi / t);
            assert!(
                rate_ratio < 2.0,
                "seed {seed}: cleave-rate spread is {rate_ratio}, which is larger \
                 than when this defect was recorded. If Task 14 introduced the E/T \
                 coupling, delete this test and assert the target band instead"
            );
        }
    }

    /// **The distinguishing test for §7.1's one uncaused property.**
    ///
    /// Task 4 derives `energy_per_unit` from the packing counts and makes both
    /// the mass defect and the decay rate functions of it, because its
    /// predecessor encoded "these bind well" twice in incommensurable units and
    /// neither encoding could check the other. A `bond_energies` matrix drawn
    /// from its own stream reproduces that defect one level up, and every other
    /// test in this file passes for such a matrix.
    ///
    /// So the assertion is that the matrix *moves*: raise one element's binding
    /// energy and every bond it makes must strengthen. A mid-table element is
    /// perturbed deliberately — at either extreme the element defines the
    /// normalising range, so its own normalised value is pinned at 0 or 1 and
    /// the test would measure the denominator instead.
    #[test]
    fn perturbing_one_elements_binding_energy_moves_every_bond_it_makes() {
        let (table, before) = universe(9);
        let all = ids(&table);

        // A mid-table element by binding energy, found without sorting floats:
        // the element whose energy is nearest the midpoint of the range.
        let lo = table.iter().fold(f64::MAX, |m, (_, e)| {
            let v = e.energy_per_unit.get();
            if v < m { v } else { m }
        });
        let hi = table.iter().fold(f64::MIN, |m, (_, e)| {
            let v = e.energy_per_unit.get();
            if v > m { v } else { m }
        });
        let mid = f64::midpoint(lo, hi);
        let (target, _) = table
            .iter()
            .fold((ElementId::ZERO, f64::MAX), |(bid, bd), (id, e)| {
                let d = (e.energy_per_unit.get() - mid).abs();
                if d < bd { (id, d) } else { (bid, bd) }
            });

        // Indexed by `target.index()`, which is the invariant `PeriodicTable`
        // now establishes: ids come from position in `iter`, so position is
        // where the element lives. `get_mut`/`if let` rather than `expect` —
        // `clippy.toml` sets neither `allow-expect-in-tests` nor
        // `allow-unwrap-in-tests`, so the workspace deny reaches inside
        // `#[cfg(test)]`.
        let mut perturbed: Vec<_> = table.iter().map(|(_, e)| e.clone()).collect();
        let bump = (hi - lo) * 0.05;
        let bumped = if let Some(e) = perturbed.get_mut(target.index()) {
            e.energy_per_unit += Quanta(bump);
            true
        } else {
            false
        };
        assert!(
            bumped,
            "the target id came from this table and must be in it"
        );
        let after_table = PeriodicTable::new(table.pattern().clone(), perturbed);
        let mut rng = Stream::new(9, Domain::Universe, 1);
        let after = BondEnergyMatrix::generate(&after_table, &mut rng);

        // Every bond the perturbed element makes strengthens...
        for &other in &all {
            let b = before.energy(target, other, BondOrder::Single);
            let a = after.energy(target, other, BondOrder::Single);
            assert!(
                a > b,
                "bond {target:?}-{other:?} did not move: {b:?} -> {a:?}"
            );
        }
        // ...and bonds between two untouched elements do not, since the
        // normalising range is unchanged by a strictly interior bump.
        let untouched: Vec<ElementId> = all.iter().copied().filter(|&i| i != target).collect();
        for &a in untouched.iter().take(8) {
            for &b in untouched.iter().take(8) {
                assert_eq!(
                    before.energy(a, b, BondOrder::Single).canonical_bits(),
                    after.energy(a, b, BondOrder::Single).canonical_bits(),
                    "bond {a:?}-{b:?} moved without either element changing"
                );
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
    fn the_zero_sentinel_is_known_wrong_signed() {
        let (table, m) = universe(9);
        // The first slot *past* the table, not an arbitrary large id: `n` is
        // exactly the value that makes `a * n + b` alias into a different
        // element's cell, so this is the id an unguarded lookup answers wrongly
        // rather than out of range. The fallback cannot fire for a 60..=120
        // table and is still out of range if it ever does.
        let past_the_end = ElementId::from_index(table.len()).unwrap_or(ElementId(u8::MAX));
        let first = ids(&table).first().copied().unwrap_or(ElementId::ZERO);
        assert_eq!(
            m.energy(past_the_end, first, BondOrder::Single),
            Quanta::ZERO
        );
        assert_eq!(
            m.energy(first, past_the_end, BondOrder::Single),
            Quanta::ZERO
        );
    }

    #[test]
    fn bond_order_rejects_what_it_cannot_represent() {
        assert_eq!(BondOrder::try_from(1), Ok(BondOrder::Single));
        assert_eq!(BondOrder::try_from(3), Ok(BondOrder::Triple));
        assert!(BondOrder::try_from(0).is_err());
        assert!(BondOrder::try_from(4).is_err());
    }
}
