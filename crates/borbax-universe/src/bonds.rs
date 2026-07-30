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
//! **`group` is not a family here.** `group` is `outer`, the occupancy of the
//! incomplete shell, so a value recurs only once per period — and `period <= 3`.
//! Measured over the same 400 universes, **at most 4 elements share a group**.
//! Keying on it would not spread bond character across a family; it would alias
//! up to four unrelated elements onto one value, so element `X`'s bonds would be
//! partly set by an element it has nothing to do with. Families do share bond
//! character in what follows, but as a *consequence* — `energy_per_unit` peaks
//! at every closure by mechanism, so closure-adjacent elements land together on
//! the scale below without anything being keyed on their group.

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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NotABondOrder(pub u8);

impl std::fmt::Display for NotABondOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "bond order {} is not in 1..=3", self.0)
    }
}

impl std::error::Error for NotABondOrder {}

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

impl BondOrder {
    /// Index into [`BondEnergyMatrix`]'s multiplier array.
    ///
    /// Written out rather than `self as usize - 1`: the `as` would be a silent
    /// conversion clippy denies workspace-wide, and spelling the three cases
    /// means adding a variant fails to compile here instead of computing an
    /// index one past the array.
    const fn scale_index(self) -> usize {
        match self {
            Self::Single => 0,
            Self::Double => 1,
            Self::Triple => 2,
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
    /// Multipliers for [`BondOrder`], indexed by `scale_index`.
    order_scale: [f64; 3],
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
    #[must_use]
    pub fn generate(table: &PeriodicTable, rng: &mut Stream) -> Self {
        let n = table.len();

        // Drawn before anything reads the table, so the draw order is a
        // property of this function rather than of the table handed to it.
        let base = rng.next_f64_range(30.0, 90.0);
        // What the *weakest* element in the table brings to a bond, as a
        // fraction of what the strongest brings.
        //
        // **This is a floor on the capacity scale, not on the finished energy,
        // and the difference is a defect the discriminator test caught.** The
        // first draft normalised binding energy onto `[0, 1]` and put the floor
        // on the bond instead. That gives the weakest element capacity exactly
        // zero — and zero is absorbing under a geometric mean, so *every* bond
        // to it scored identically no matter the partner. The weakest element
        // is not an edge case: `contacts_upto(1) = 0`, so it is reliably the
        // monomer, the feedstock, the most abundant species in the beaker. The
        // most common bond in the chemistry would have had no structure at all,
        // which is the §22.6 failure this file's spread test exists to catch and
        // could not see, because the *range* over the whole matrix stayed wide.
        //
        // Below 1/3 makes `energies_have_useful_spread` structural rather than
        // lucky: the extremes are the weakest and strongest elements bonded to
        // themselves, so the ratio is exactly `1 / weakest_share`.
        let weakest_share = rng.next_f64_range(0.12, 0.28);
        // 2.2 > 2.1, so triple strictly exceeds double in every universe. The
        // ranges are adjacent on purpose; overlapping them would make
        // `higher_orders_are_stronger` a property of the draw.
        let order_scale = [
            1.0,
            rng.next_f64_range(1.6, 2.1),
            rng.next_f64_range(2.2, 3.0),
        ];

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
            for b in a..n {
                // Geometric mean, so the pair is symmetric bit-for-bit
                // (multiplication commutes exactly) and a weak partner drags the
                // bond down rather than being rescued by a strong one the way a
                // sum would allow — a bond needs both ends to hold. `sqrt` is
                // exactly specified by IEEE-754 and stays native (§13.1).
                let (ca, cb) = (
                    capacity.get(a).copied().unwrap_or(1.0),
                    capacity.get(b).copied().unwrap_or(1.0),
                );
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
    /// Takes ids rather than `&Element` deliberately: a caller in a
    /// reaction-rate path holds ids, and requiring `&Element` would couple that
    /// loop to the element table's lifetime and to a struct 48 bytes of which
    /// are `String`.
    ///
    /// **An id this universe never minted scores [`Quanta::ZERO`]** — no bond,
    /// hence inert. An `Option` was the alternative and lost on what it would
    /// cost every caller: the ids reaching here come from
    /// [`PeriodicTable::iter`], so `None` is unreachable in correct code and
    /// would be unwrapped at every site in a workspace that denies `unwrap`.
    /// The case that does reach it is an id from a *different* universe, and
    /// zero is the conservative answer there — pinned by
    /// `an_id_outside_the_table_yields_no_energy` rather than left to the
    /// indexing to decide.
    #[must_use]
    pub fn energy(&self, a: ElementId, b: ElementId, order: BondOrder) -> Quanta {
        let (ia, ib) = (a.index(), b.index());
        if ia >= self.n || ib >= self.n {
            return Quanta::ZERO;
        }
        let scale = self
            .order_scale
            .get(order.scale_index())
            .copied()
            .unwrap_or(1.0);
        self.e
            .get(ia * self.n + ib)
            .copied()
            .unwrap_or(Quanta::ZERO)
            * scale
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
    #[test]
    fn energies_are_symmetric() {
        let (table, m) = universe(4);
        for &a in &ids(&table) {
            for &b in &ids(&table) {
                assert_eq!(
                    m.energy(a, b, BondOrder::Single).canonical_bits(),
                    m.energy(b, a, BondOrder::Single).canonical_bits(),
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
            assert!(
                (ratio - 1.0 / weakest_share).abs() < 1e-9,
                "seed {seed}: spread {ratio} is not 1/{weakest_share}"
            );
            assert!(
                ratio > 3.0,
                "seed {seed}: spread {ratio} fell to the spread floor"
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

    /// Spread matters: a matrix where every bond costs the same gives a
    /// chemistry with no structure to exploit. Some linkages must be durable
    /// and others fragile, or decay has nothing to select on (spec §22.6).
    #[test]
    fn energies_have_useful_spread() {
        for seed in 0..20 {
            let (table, m) = universe(seed);
            let all: Vec<f64> = ids(&table)
                .iter()
                .flat_map(|&a| {
                    ids(&table)
                        .into_iter()
                        .map(move |b| (a, b))
                        .collect::<Vec<_>>()
                })
                .map(|(a, b)| m.energy(a, b, BondOrder::Single).get())
                .collect();
            let lo = all
                .iter()
                .copied()
                .fold(f64::MAX, |x, y| if x < y { x } else { y });
            let hi = all
                .iter()
                .copied()
                .fold(f64::MIN, |x, y| if x > y { x } else { y });
            assert!(
                hi / lo > 3.0,
                "seed {seed}: bond energies too uniform: {lo}..{hi}"
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

        // `find`/`if let` rather than `expect`: `clippy.toml` sets neither
        // `allow-expect-in-tests` nor `allow-unwrap-in-tests`, so the workspace
        // deny reaches inside `#[cfg(test)]`.
        let mut perturbed: Vec<_> = table.iter().map(|(_, e)| e.clone()).collect();
        let bump = (hi - lo) * 0.05;
        let mut bumped = false;
        for e in &mut perturbed {
            if e.id == target {
                e.energy_per_unit += Quanta(bump);
                bumped = true;
            }
        }
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

    /// An id from another universe's larger table must not silently read a
    /// neighbouring element's energy.
    #[test]
    fn an_id_outside_the_table_yields_no_energy() {
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
