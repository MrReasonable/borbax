//! Small-molecule atom graphs (spec §8.1).
//!
//! Capped at twelve atoms. Adjacency is stored as bitset rows, one plane per
//! representable bond order. That matters more than it might appear: colour
//! refinement in `canonical.rs` becomes a sequence of popcounts over masks,
//! with no allocation and no hashing anywhere in it.

use borbax_units::Mass;
use borbax_universe::{BondOrder, ElementId, PeriodicTable};

/// The atom cap. Twelve, so an adjacency row is a `u16` and every reachability
/// sweep is a handful of popcounts.
pub const MAX_ATOMS: usize = 12;

/// **The row width bound, asserted where the row lives.**
///
/// An adjacency row is a `u16`, so an atom index is shifted into it and the cap
/// must not exceed its width. `canonical.rs` asserts the *corresponding* `u128`
/// bound for its pair bitsets, and `N_PAIRS <= 128` happens to coincide with
/// `MAX_ATOMS <= 16` exactly — by accident, in another file, about another
/// array. A bound that holds by coincidence somewhere else is not enforced here.
const _: () = assert!(MAX_ATOMS <= 16, "an adjacency row is a u16");

/// One bitset plane per representable bond order.
///
/// **`BondOrder::PLANES`, not `BondOrder::ALL.len()`.** This array is indexed at
/// `order - 1`, so what it needs is `MAX` — the *count* of representable orders
/// is an adjacent quantity that currently coincides. They agree only while `ALL`
/// is dense, which is guaranteed by a `const _` inside `borbax-universe` that
/// this crate can neither see nor cite. A sparse `ALL` — `[1, 2, 3, 6]`, not
/// absurd while formability is being decided — would make `ALL.len()` 4 while
/// `MAX` is 6, and order 6 would index plane 5 of a 4-plane array.
/// `every_representable_order_has_a_plane` asserts the relation here rather than
/// trusting the invisible one.
pub const N_ORDERS: usize = BondOrder::PLANES;

/// Why a bond could not be formed.
///
/// **Every variant is a caller bug, and none of them is a silent early return.**
/// [`Mol12::add_atom`] argues the opposite for its own out-of-range case and the
/// two are not inconsistent: a molecule at capacity is chemically meaningful — it
/// simply does not react that way — whereas an out-of-range atom index or an
/// unpayable bond order has no chemical reading at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum BondError {
    /// The two indices were equal, or one of them named no atom.
    #[error("no bond between atoms {a} and {b}: a molecule of {n} cannot have both, distinctly")]
    NoSuchAtom {
        /// The first index as given.
        a: u8,
        /// The second index as given.
        b: u8,
        /// How many atoms the molecule actually has.
        n: u8,
    },
    /// The atom's element is not in the table this bond was checked against.
    #[error("atom {atom} holds an element id the table does not have")]
    NoSuchElement {
        /// The atom whose element could not be resolved.
        atom: u8,
    },
    /// The bond would spend slots the atom does not have.
    ///
    /// Carries the atom rather than the pair, because the per-atom budget is
    /// what was exceeded — `min(valence_a, valence_b)` cannot say which end ran
    /// out, and at twelve atoms it admits eleven times the budget.
    #[error("atom {atom} has valence {valence} but the bond would take it to {would_use}")]
    ValenceExceeded {
        /// The end that ran out of slots.
        atom: u8,
        /// That atom's element's valence — its whole budget.
        valence: u8,
        /// What the budget would have had to be for the bond to form.
        would_use: u32,
    },
}

/// **Fields are private, and `elem` holds [`ElementId`] rather than `u8`.**
///
/// Two findings, one shape. A public `adj` lets any write set two planes for one
/// pair, after which [`Mol12::bond_order`] returns the *lowest* order and
/// [`Mol12::valence_used`] counts both — an illegal state maintained only by
/// `add_bond` remembering to call `clear_bond` first. `pub(crate)` would not
/// close it either, because Task 11's `Polymer` will be the first thing to build
/// a `Mol12` by another route and it lives in this same crate; private fields
/// plus read-only accessors is what actually holds.
///
/// And `[u8; _]` cannot be written or read at all from here: `ElementId`'s field
/// and its `index()` are both `pub(crate)` in `borbax-universe`, deliberately, so
/// `self.elem[i] = e.0` is a hard compile error outside that crate.
/// `[ElementId; MAX_ATOMS]` is the same byte — asserted by a `const _` there —
/// stays `Copy`, and cannot hold a value that was never an id.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mol12 {
    n: u8,
    elem: [ElementId; MAX_ATOMS],
    /// `adj[order.plane_index()][i]` is a bitset of atoms bonded to `i` at that
    /// order.
    adj: [[u16; MAX_ATOMS]; N_ORDERS],
}

impl Default for Mol12 {
    fn default() -> Self {
        Self::new()
    }
}

impl Mol12 {
    /// The empty molecule: no atoms, no bonds.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            n: 0,
            elem: [ElementId::ZERO; MAX_ATOMS],
            adj: [[0; MAX_ATOMS]; N_ORDERS],
        }
    }

    /// How many atoms this molecule has.
    ///
    /// Not `const`, so this can use `usize::from` rather than `as` — the same
    /// trade `BondOrder::plane_index` makes, and for the same reason: a silent
    /// widening cast is the spelling `clippy::as_conversions` exists to refuse.
    #[must_use]
    pub fn len(&self) -> usize {
        usize::from(self.n)
    }

    /// Whether this molecule has no atoms at all.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The element at atom `i`, or `None` if there is no such atom.
    #[must_use]
    pub fn element(&self, i: u8) -> Option<ElementId> {
        if i < self.n {
            self.elem.get(usize::from(i)).copied()
        } else {
            None
        }
    }

    /// Append an atom, returning its index.
    ///
    /// `None` at capacity — a molecule that cannot grow is chemically meaningful
    /// (it simply does not react that way), whereas a panic here would take down
    /// a simulation step. The return is `#[must_use]` so "silently ignored"
    /// cannot be what a caller accidentally chooses.
    #[must_use = "at capacity this adds nothing, and the index is how you know"]
    pub fn add_atom(&mut self, e: ElementId) -> Option<u8> {
        let idx = self.n;
        let slot = self.elem.get_mut(usize::from(idx))?;
        *slot = e;
        self.n += 1;
        Some(idx)
    }

    /// The adjacency bit for atom `b`, or `None` if `b` names no atom.
    ///
    /// **Every shift on a caller-supplied index goes through here, and that is
    /// the point.** (`is_connected` shifts a loop counter it just bounded by
    /// `self.n`; the qualifier is there because an unqualified "every shift in
    /// this file" was false, and CLAUDE.md records two shipped defects that were
    /// exactly a doc asserting a bound it did not have.)
    /// An earlier version spelled `1u16 << u32::from(b)` inline in
    /// [`Self::bond_order`] and twice more in `clear_bond`, guarded in neither —
    /// so the bound was enforced by `row`/`set_row`'s `.get()` calls, which
    /// happen *after* the shift has already been evaluated. Measured on the
    /// committed code: `b = 17` panicked with "attempt to shift left with
    /// overflow" in debug, and in release the shift **masked** to `b % 16`, so
    /// `bond_order(0, 17)` answered `Some(SINGLE)` for a bond that did not exist
    /// and `clear_bond(0, 17)` cleared one row of a pair, leaving adjacency
    /// asymmetric — `degree(0) = 0` against `degree(1) = 1`, `is_connected`
    /// flipped, and a different `CanonForm`. Green in `cargo test`, wrong in the
    /// profile that mints goldens, which is the split CLAUDE.md names twice.
    ///
    /// Returning the *bit* rather than validating an index is what makes the
    /// hazard unspellable: there is no way to reach the shift without having
    /// asked whether the atom exists.
    fn bit(&self, b: u8) -> Option<u16> {
        // `b < self.n <= MAX_ATOMS <= 16` (asserted above), so the shift is in
        // range by construction rather than by inspection.
        (b < self.n).then(|| 1u16 << u32::from(b))
    }

    /// One adjacency row, or an empty one for an index that names no atom.
    ///
    /// Out-of-range reads as an empty row rather than panicking. Note this is
    /// **not** on its own what keeps the public readers total — the shift in
    /// [`Self::bit`] is the part that had to be guarded, and this `.get()` is
    /// what makes the row index safe once it has been.
    fn row(&self, plane: usize, a: u8) -> u16 {
        self.adj
            .get(plane)
            .and_then(|p| p.get(usize::from(a)))
            .copied()
            .unwrap_or(0)
    }

    fn set_row(&mut self, plane: usize, a: u8, bits: u16) {
        if let Some(slot) = self
            .adj
            .get_mut(plane)
            .and_then(|p| p.get_mut(usize::from(a)))
        {
            *slot = bits;
        }
    }

    /// Add or replace a bond, refusing any order the two atoms cannot pay for.
    ///
    /// **Takes a [`PeriodicTable`], and the reason it is not a `&Universe` is
    /// that valence is all it needs.** The bond-energy matrix and the universe
    /// constants have nothing to say about whether a bond can form. It is not
    /// two `&Element`s either: the molecule already knows which ids sit at `a`
    /// and `b`, and a caller passing the elements could pass ones that do not
    /// match those atoms — a fresh illegal state, in the function whose whole
    /// job is refusing them.
    ///
    /// **The budget is per-atom, not per-bond.** Task 5 routed
    /// `min(valence_a, valence_b)` here, which is necessary and *not sufficient*:
    /// a hub with eleven single bonds to eleven valence-1 partners satisfies it
    /// at every bond and overspends its own budget elevenfold. Maintaining
    /// `valence_used(x) + order <= valence(x)` at each addition gives the
    /// whole-molecule property by induction, and implies the pairwise check at
    /// the first bond — see `the_per_atom_budget_outlives_the_pairwise_check`.
    ///
    /// The order is a [`BondOrder`], so `0` and `7` are unrepresentable rather
    /// than guarded. A predecessor took `u8` and computed `(order - 1) as usize`,
    /// where `0` underflows in debug and wraps to 255 in release.
    ///
    /// # Errors
    ///
    /// [`BondError::NoSuchAtom`] if the indices are equal or either names no
    /// atom; [`BondError::NoSuchElement`] if an atom's element is absent from
    /// `table`; [`BondError::ValenceExceeded`] if either end lacks the slots.
    ///
    /// **A refused bond leaves the molecule exactly as it was found, and that is
    /// now structural rather than repaired.** The predecessor cleared the pair
    /// first — so a *replacement* got its old slots back — and then restored the
    /// old bond on the `ValenceExceeded` path. But the two element lookups were
    /// bare `?` sitting after that clear and restored nothing, so a molecule
    /// checked against a table that does not hold its ids returned `Err` **with
    /// the bond already deleted**. Reproduced: two atoms of one id bonded under
    /// an 89-element table, then re-stated against an 80-element one.
    /// `a_refused_bond_leaves_the_molecule_untouched` covered only the valence
    /// branch, so two of the three documented paths were unenforced.
    ///
    /// The repair is not a third restore. Refunding a replaced bond's order
    /// *arithmetically* means every fallible check happens before any mutation,
    /// so there is no error path after a write for a fourth one to be forgotten
    /// on.
    pub fn add_bond(
        &mut self,
        a: u8,
        b: u8,
        order: BondOrder,
        table: &PeriodicTable,
    ) -> Result<(), BondError> {
        let (Some(bit_a), Some(bit_b)) = (self.bit(a), self.bit(b)) else {
            return Err(BondError::NoSuchAtom { a, b, n: self.n });
        };
        if a == b {
            return Err(BondError::NoSuchAtom { a, b, n: self.n });
        }

        // A replacement must not be charged twice, so the bond being replaced is
        // refunded from the budget rather than cleared out of the molecule.
        // `saturating_sub` cannot actually saturate: a bond of order `o` between
        // `a` and `b` contributes `o` to `valence_used` at *both* ends, so the
        // refund never exceeds either end's usage. It is spelled saturating so a
        // future change cannot turn that reasoning into an underflow.
        let refund = self.bond_order(a, b).map_or(0, |o| u32::from(u8::from(o)));
        let cost = u32::from(u8::from(order));

        for atom in [a, b] {
            let id = self
                .element(atom)
                .ok_or(BondError::NoSuchElement { atom })?;
            let valence = table
                .get(id)
                .ok_or(BondError::NoSuchElement { atom })?
                .valence;
            let would_use = self.valence_used(atom).saturating_sub(refund) + cost;
            if would_use > u32::from(valence) {
                return Err(BondError::ValenceExceeded {
                    atom,
                    valence,
                    would_use,
                });
            }
        }

        // Past every fallible step, so this is a commit rather than a mutation
        // that might need undoing.
        self.write_pair(a, b, bit_a, bit_b, order);
        Ok(())
    }

    /// Set the plane bits for a pair, **clearing every other plane first**.
    ///
    /// The predecessor was documented "assumes the pair is already clear" and
    /// left the clearing to its caller. A restructure then reworded the doc to
    /// "whose bits are already known valid" — which is about the bit *masks* —
    /// and the precondition vanished with nothing in its place.
    ///
    /// Two planes set for one pair is the illegal state [`Mol12`]'s own doc
    /// names: [`Mol12::bond_order`] returns the *lowest* order while
    /// [`Mol12::valence_used`] counts them all, so `encode` and `refine`
    /// disagree about the graph. Measured: a pair written SINGLE then TRIPLE
    /// interns as the **single-bonded species** while carrying four spent
    /// valence slots — a species collision, in which a whole class of
    /// higher-order bonds silently collapses onto the single-bonded form.
    ///
    /// So the clear happens here rather than in a comment asking callers to do
    /// it. The struct doc names Task 11's `Polymer` as the next in-crate builder,
    /// which is exactly who the old precondition was waiting for. Same six
    /// `set_row`s either way, one caller obligation fewer — and a clear cannot
    /// be absent from the profile that mints goldens the way a `debug_assert`
    /// can, which is why it is not one.
    fn write_pair(&mut self, a: u8, b: u8, bit_a: u16, bit_b: u16, order: BondOrder) {
        self.clear_pair(a, b, bit_a, bit_b);
        let plane = order.plane_index();
        self.set_row(plane, a, self.row(plane, a) | bit_b);
        self.set_row(plane, b, self.row(plane, b) | bit_a);
    }

    /// Clear a pair whose bits are already known valid.
    fn clear_pair(&mut self, a: u8, b: u8, bit_a: u16, bit_b: u16) {
        for plane in 0..N_ORDERS {
            self.set_row(plane, a, self.row(plane, a) & !bit_b);
            self.set_row(plane, b, self.row(plane, b) & !bit_a);
        }
    }

    /// Remove any bond between `a` and `b`, at whatever order it was held.
    ///
    /// **Fallible, and symmetric with [`Self::add_bond`] on purpose.** It has no
    /// caller yet; the first will be Task 12/13's bond breaking, taking indices
    /// from a reaction site. A silent no-op on a bad index there would be a
    /// species split nothing reports, so the refusal is a value the caller has
    /// to handle rather than a comment promising it cannot happen.
    ///
    /// # Errors
    ///
    /// [`BondError::NoSuchAtom`] if either index names no atom. The molecule is
    /// unchanged in that case — the check precedes every write.
    pub fn clear_bond(&mut self, a: u8, b: u8) -> Result<(), BondError> {
        let (Some(bit_a), Some(bit_b)) = (self.bit(a), self.bit(b)) else {
            return Err(BondError::NoSuchAtom { a, b, n: self.n });
        };
        self.clear_pair(a, b, bit_a, bit_b);
        Ok(())
    }

    /// The order of the bond between `a` and `b`, if there is one.
    ///
    /// **`Option<BondOrder>`, not a `0`-means-absent `u8`.** The sentinel form is
    /// what let a predecessor's `permute` feed a `0` back into `add_bond`, and it
    /// is the same shape `bonds.rs` spends forty lines condemning in
    /// `Quanta::ZERO` — in a type where `0` is unrepresentable and
    /// `Option<BondOrder>` is two bytes.
    /// `None` also when either index names no atom — asking about a bond
    /// between atoms that do not exist has one honest answer, and the private
    /// `bit` helper is what makes the shift unreachable until that has been
    /// established.
    #[must_use]
    pub fn bond_order(&self, a: u8, b: u8) -> Option<BondOrder> {
        let bit_b = self.bit(b)?;
        self.bit(a)?;
        BondOrder::ALL
            .into_iter()
            .find(|o| self.row(o.plane_index(), a) & bit_b != 0)
    }

    /// Bitset of all neighbours of `a`, at any order.
    #[must_use]
    pub fn neighbours(&self, a: u8) -> u16 {
        (0..N_ORDERS).fold(0, |acc, plane| acc | self.row(plane, a))
    }

    /// Bitset of the neighbours of `a` held at exactly this order.
    ///
    /// `canonical.rs` needs the planes separately — a vertex signature counts
    /// neighbours per (order, colour) pair — and this is how it reads them
    /// without the fields being visible.
    #[must_use]
    pub fn neighbours_at(&self, order: BondOrder, a: u8) -> u16 {
        self.row(order.plane_index(), a)
    }

    /// Number of bonded neighbours, regardless of order.
    #[must_use]
    pub fn degree(&self, a: u8) -> u32 {
        self.neighbours(a).count_ones()
    }

    /// Total valence consumed by bonds at `a`, counting order.
    ///
    /// Integer addition over a fixed-length array in plane order, so there is no
    /// accumulation-order question here (§13.1 is about float sums and
    /// scheduler-ordered reductions).
    #[must_use]
    pub fn valence_used(&self, a: u8) -> u32 {
        BondOrder::ALL
            .into_iter()
            .map(|o| u32::from(u8::from(o)) * self.row(o.plane_index(), a).count_ones())
            .sum()
    }

    /// Exact total mass, or `None` if any atom's element is absent from `table`.
    ///
    /// Integer addition, so this is a guarantee rather than a tolerance (§13.1).
    /// Fallible because `add_atom` accepts any [`ElementId`] and an id is only
    /// meaningful against a particular table — answering with a silent zero for a
    /// missing element would put a wrong mass into a conservation check whose
    /// whole purpose is being exact.
    #[must_use]
    pub fn mass(&self, table: &PeriodicTable) -> Option<Mass> {
        (0..self.n)
            .map(|i| self.element(i).and_then(|id| table.get(id)).map(|e| e.mass))
            .sum()
    }

    /// This molecule with its atoms renumbered: `order[k]` is the atom that ends
    /// up at position `k`.
    ///
    /// That is the same convention `canonical.rs`'s `encode` uses, deliberately —
    /// [`canonicalise`](crate::canonicalise) applies the winning search order
    /// through here to build the canonical molecule, so the two must agree on
    /// which direction the permutation runs. Reading it backwards produces the
    /// inverse relabelling, which is *also* a legal molecule and also isomorphic,
    /// so nothing downstream fails — it silently interns a different labelling.
    ///
    /// **No [`PeriodicTable`] and no valence check, and that is sound rather than
    /// a shortcut.** Valence is a per-atom sum of bond orders; a relabelling
    /// carries every bond to a renumbered pair with its order intact, so the sum
    /// at atom `order[k]` becomes the sum at `k` unchanged. A relabelling of a
    /// legal molecule is legal, which is why this can write the rows directly.
    ///
    /// `None` if `order` is not a permutation of `0..n` — a repeated entry would
    /// duplicate one atom and drop another, giving a molecule with the right atom
    /// count and the wrong contents, which is exactly the failure that would
    /// survive every test asserting only a count.
    #[must_use]
    pub fn relabelled(&self, order: &[u8]) -> Option<Self> {
        let n = usize::from(self.n);
        if order.len() != n {
            return None;
        }

        // Permutation check and the inverse in one pass. `inv[v]` is the new
        // position of original atom `v`; `bit` rejects any `v >= n` before it is
        // used as an index, which is the same guard every public reader goes
        // through.
        let mut seen: u16 = 0;
        let mut inv = [0u8; MAX_ATOMS];
        for (k, &v) in order.iter().enumerate() {
            let bit = self.bit(v)?;
            if seen & bit != 0 {
                return None;
            }
            seen |= bit;
            *inv.get_mut(usize::from(v))? = u8::try_from(k).ok()?;
        }

        let mut out = Self::new();
        out.n = self.n;
        for (k, &v) in order.iter().enumerate() {
            *out.elem.get_mut(k)? = *self.elem.get(usize::from(v))?;
        }

        // Rows are rebuilt bit by bit rather than permuted wholesale: a row is a
        // bitset over *atom indices*, so moving the row to its new owner is only
        // half the job — every bit inside it names an atom that has also moved.
        // Permuting only the outer index leaves adjacency asymmetric, which
        // `bond_order` would then answer inconsistently depending on which end
        // it is asked from.
        //
        // The set bits are walked with `trailing_zeros` rather than by scanning
        // every column. That is not a micro-optimisation of a hot loop — this
        // runs once per `canonicalise` — but the scan visited `N_ORDERS * n²`
        // = 864 positions at n = 12 to find the `deg(v) <= valence` that are
        // actually set, and the walk is a measured **5.0x**. It is taken because
        // `relabelled_agrees_with_replaying_the_molecule_through_add_bond` checks
        // it against an independent implementation over 2000 (molecule,
        // permutation) pairs — and exhaustively by review over all 50,360
        // (row, permutation) pairs for n = 2..6 — so the equivalence is
        // established rather than argued.
        //
        // **Scoped to rows with no bits at or above `n`**, which is where the two
        // forms agree. Out of that range they diverge in 13 of 13 constructed
        // cases: the walk reads `inv[w]` for a slot `order` never wrote and sets
        // bit 0 spuriously, where the scan drops it. That range is unreachable —
        // `self.n` is only ever incremented, and every `set_row` call site derives
        // its bits from `bit()`, which refuses any index >= n — but the precondition
        // is shared rather than local, so it is stated here rather than assumed.
        for plane in 0..N_ORDERS {
            for (k, &v) in order.iter().enumerate() {
                let mut src = self.row(plane, v);
                let mut dst: u16 = 0;
                while src != 0 {
                    let w = src.trailing_zeros();
                    src &= src - 1;
                    let moved = *inv.get(usize::try_from(w).ok()?)?;
                    dst |= self.bit(moved)?;
                }
                out.set_row(plane, u8::try_from(k).ok()?, dst);
            }
        }
        Some(out)
    }

    /// Whether every atom is reachable from atom zero.
    #[must_use]
    pub fn is_connected(&self) -> bool {
        if self.n == 0 {
            return true;
        }
        let mut seen: u16 = 1;
        let mut frontier: u16 = 1;
        while frontier != 0 {
            let mut next: u16 = 0;
            for i in 0..self.n {
                if frontier & (1u16 << u32::from(i)) != 0 {
                    next |= self.neighbours(i);
                }
            }
            next &= !seen;
            seen |= next;
            frontier = next;
        }
        seen.count_ones() == u32::from(self.n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use borbax_universe::element::generate_elements;
    use borbax_universe::{BondOrder, ElementId, PeriodicTable};

    /// A real drawn element table.
    ///
    /// **`generate_elements`, not `Universe::generate`.** A graph needs element
    /// *properties* — valence to bound bond orders, mass to sum — and nothing
    /// from the bond-energy matrix or the universe constants. Building a whole
    /// `Universe` here would borrow two subsystems this file never consults and
    /// would put `Universe` in `graph.rs`'s dependency surface, where it does
    /// not belong until Task 8 embeds something.
    fn table() -> PeriodicTable {
        generate_elements(5)
    }

    /// The lowest-id element with exactly this valence.
    ///
    /// Exact rather than "at least", because the valence tests need to know the
    /// budget precisely — "at least 1" would let a valence-6 element satisfy a
    /// test meaning to exhaust a valence-1 one, and it would pass for the wrong
    /// reason.
    ///
    /// `unreachable!` rather than a silent skip: seed 5's table is fixed, so
    /// these classes provably exist, and if element generation ever moves these
    /// fixtures must fail loudly rather than quietly stop testing anything.
    /// `clippy::panic` is denied workspace-wide and reaches inside
    /// `#[cfg(test)]`, so this is the sanctioned spelling.
    fn valence_exactly(t: &PeriodicTable, v: u8) -> ElementId {
        t.iter().find(|(_, e)| e.valence == v).map_or_else(
            || unreachable!("seed 5's table has no valence-{v} element; the fixture is stale"),
            |(id, _)| id,
        )
    }

    /// The first `count` distinct elements with at least `v` valence.
    fn bondable(t: &PeriodicTable, v: u8, count: usize) -> Vec<ElementId> {
        let ids: Vec<ElementId> = t
            .iter()
            .filter(|(_, e)| e.valence >= v)
            .map(|(id, _)| id)
            .take(count)
            .collect();
        assert_eq!(
            ids.len(),
            count,
            "seed 5's table has too few valence-{v}+ elements"
        );
        ids
    }

    /// The `i`th element of a fixture list, cycling.
    fn nth(ids: &[ElementId], i: usize) -> ElementId {
        ids.get(i % ids.len()).map_or_else(
            || unreachable!("bondable() guarantees a non-empty list"),
            |id| *id,
        )
    }

    fn order(n: u8) -> BondOrder {
        BondOrder::new(n).unwrap_or_else(|| unreachable!("{n} is outside 1..=BondOrder::MAX"))
    }

    /// Append an atom, asserting the molecule had room.
    fn atom(mol: &mut Mol12, elem: ElementId) -> u8 {
        mol.add_atom(elem)
            .unwrap_or_else(|| unreachable!("the molecule was unexpectedly at capacity"))
    }

    /// Add a bond that the fixture expects to be legal.
    fn bond(mol: &mut Mol12, a: u8, b: u8, ord: BondOrder, table: &PeriodicTable) {
        assert!(
            mol.add_bond(a, b, ord, table).is_ok(),
            "fixture bond {a}-{b} at order {ord:?} was refused"
        );
    }

    /// An unbranched chain of `n` atoms, cycling through four distinct elements
    /// so the mass sum is a sum of *different* numbers rather than a multiple.
    fn chain(t: &PeriodicTable, n: u8) -> Mol12 {
        // Interior atoms carry two single bonds, so two slots is the floor.
        let ids = bondable(t, 2, 4);
        let mut m = Mol12::new();
        for i in 0..n {
            atom(&mut m, nth(&ids, usize::from(i)));
        }
        for i in 0..n.saturating_sub(1) {
            bond(&mut m, i, i + 1, BondOrder::SINGLE, t);
        }
        m
    }

    #[test]
    fn bonds_are_symmetric() {
        let t = table();
        let m = chain(&t, 4);
        for i in 0..4u8 {
            for j in 0..4u8 {
                assert_eq!(m.bond_order(i, j), m.bond_order(j, i));
            }
        }
    }

    #[test]
    fn degree_counts_all_orders() {
        let t = table();
        // Atom 0 carries a single *and* a triple, so it needs four slots.
        let hub = valence_exactly(&t, 4);
        let mut m = Mol12::new();
        for _ in 0..3 {
            atom(&mut m, hub);
        }
        bond(&mut m, 0, 1, BondOrder::SINGLE, &t);
        bond(&mut m, 0, 2, order(3), &t);
        assert_eq!(m.degree(0), 2, "two neighbours, regardless of order");
        assert_eq!(m.degree(1), 1);
        assert_eq!(m.valence_used(0), 4, "one slot plus three");
    }

    #[test]
    fn mass_is_the_exact_sum_of_atoms() {
        let t = table();
        let ids = bondable(&t, 2, 4);
        let m = chain(&t, 6);
        let expected: Option<borbax_units::Mass> =
            (0..6).map(|i| t.get(nth(&ids, i)).map(|e| e.mass)).sum();
        // Stated rather than assumed: a `None == None` comparison below would
        // pass while measuring nothing.
        assert!(
            expected.is_some(),
            "every fixture element must be in the table"
        );
        // Exact equality, not a tolerance — this is why Mass is fixed-point.
        assert_eq!(m.mass(&t), expected);
    }

    #[test]
    fn connectivity_is_detected() {
        let t = table();
        assert!(chain(&t, 5).is_connected());
        let mut m = chain(&t, 4);
        atom(&mut m, valence_exactly(&t, 2)); // isolated fifth atom
        assert!(!m.is_connected());
    }

    #[test]
    fn capacity_is_respected() {
        let t = table();
        let e = valence_exactly(&t, 2);
        let mut m = Mol12::new();
        let mut accepted = 0;
        for _ in 0..MAX_ATOMS + 3 {
            if m.add_atom(e).is_some() {
                accepted += 1;
            }
        }
        assert_eq!(
            accepted, MAX_ATOMS,
            "add_atom reports refusal past capacity"
        );
        assert_eq!(m.len(), MAX_ATOMS);
    }

    /// **A refused bond must not be a mutation.**
    ///
    /// `add_bond` clears the pair *before* checking the budget, so that replacing
    /// a bond returns its slots. That ordering means a refusal happens with the
    /// old bond already gone — so without an explicit restore, a caller that
    /// handles the error correctly still silently loses a bond it had. The
    /// `# Errors` section promises this; a promise with no test is the shape this
    /// project keeps finding.
    #[test]
    fn a_refused_bond_leaves_the_molecule_untouched() {
        let t = table();
        let two = valence_exactly(&t, 2);
        let mut m = Mol12::new();
        for _ in 0..2 {
            atom(&mut m, two);
        }
        bond(&mut m, 0, 1, order(2), &t);
        let before = m;

        // Both ends are now full, so re-stating the pair at a *higher* order is
        // refused — after the existing double has already been cleared.
        assert!(matches!(
            m.add_bond(0, 1, order(3), &t),
            Err(BondError::ValenceExceeded { .. })
        ));
        assert_eq!(
            m, before,
            "the refused bond destroyed the one that was there"
        );
        assert_eq!(m.bond_order(0, 1), Some(order(2)));
    }

    /// **All three error variants, because two of them were unenforced.**
    ///
    /// The test above covered `ValenceExceeded` only. `NoSuchElement` is
    /// reachable only with two tables of different lengths — `add_atom` accepts
    /// any `ElementId`, and an id is only meaningful against a particular table —
    /// and on that path the predecessor had already cleared the bond before the
    /// lookup failed. It returned `Err` as documented, with the molecule
    /// silently changed, which changes its `CanonForm` and interns it as a
    /// different species.
    #[test]
    fn every_refusal_leaves_the_molecule_untouched() {
        let big = table(); // seed 5
        let small = generate_elements(1);
        assert!(
            small.len() < big.len(),
            "this test needs an id valid in one table and absent from the other"
        );

        // An id the big table holds and the small one does not.
        let stranger = ElementId::from_index(small.len())
            .unwrap_or_else(|| unreachable!("a table length fits u8"));
        assert!(big.get(stranger).is_some() && small.get(stranger).is_none());

        let mut m = Mol12::new();
        atom(&mut m, stranger);
        atom(&mut m, stranger);
        bond(&mut m, 0, 1, BondOrder::SINGLE, &big);
        let before = m;

        // NoSuchElement: the ids are fine here, absent over there.
        assert!(matches!(
            m.add_bond(0, 1, order(2), &small),
            Err(BondError::NoSuchElement { .. })
        ));
        assert_eq!(m, before, "a table mismatch deleted the bond it refused");
        assert_eq!(m.bond_order(0, 1), Some(BondOrder::SINGLE));

        // NoSuchAtom, both spellings.
        assert!(matches!(
            m.add_bond(0, 0, BondOrder::SINGLE, &big),
            Err(BondError::NoSuchAtom { .. })
        ));
        assert_eq!(m, before);
        assert!(matches!(
            m.add_bond(0, 9, BondOrder::SINGLE, &big),
            Err(BondError::NoSuchAtom { .. })
        ));
        assert_eq!(m, before);
        // The comment used to say "both spellings" while testing only an
        // out-of-range `b`. `a` goes through the same let-else, so this is
        // coverage rather than a second code path — but the claim was untrue.
        assert!(matches!(
            m.add_bond(9, 0, BondOrder::SINGLE, &big),
            Err(BondError::NoSuchAtom { .. })
        ));
        assert_eq!(m, before);
    }

    /// **At most one plane may hold any pair — the invariant `Mol12`'s struct
    /// doc names and nothing asserted.**
    ///
    /// Two planes for one pair is not a slow answer, it is a species collision:
    /// `bond_order` returns the *lowest* order while `valence_used` counts all
    /// of them, so `encode` and `refine` disagree about what graph they are
    /// looking at. Measured on a build where `write_pair` did not clear first, a
    /// pair written SINGLE then TRIPLE interned as the single-bonded species
    /// while carrying four spent valence slots.
    ///
    /// Asserted over the whole random corpus rather than one fixture, because
    /// the states that reach it are the ones a *future* builder constructs.
    #[test]
    fn a_pair_never_occupies_two_planes() {
        let t = table();
        let two = valence_exactly(&t, 2);
        let four = valence_exactly(&t, 4);

        let mut m = Mol12::new();
        atom(&mut m, four);
        atom(&mut m, four);
        // Every replacement path, in both directions.
        for (first, second) in [(1u8, 3u8), (3, 1), (2, 2), (1, 4), (4, 1)] {
            bond(&mut m, 0, 1, order(first), &t);
            bond(&mut m, 0, 1, order(second), &t);
            let planes = BondOrder::ALL
                .into_iter()
                .filter(|o| m.neighbours_at(*o, 0) & (1u16 << 1) != 0)
                .count();
            assert_eq!(
                planes, 1,
                "order {first} then {second} left {planes} planes set"
            );
            assert_eq!(m.bond_order(0, 1), Some(order(second)));
            assert_eq!(
                m.valence_used(0),
                u32::from(second),
                "valence_used counts every plane, so it exposes a stale one"
            );
        }

        // And a chain, where clear_bond is the other route in.
        let mut c = Mol12::new();
        for _ in 0..3 {
            atom(&mut c, two);
        }
        bond(&mut c, 0, 1, BondOrder::SINGLE, &t);
        bond(&mut c, 1, 2, BondOrder::SINGLE, &t);
        assert!(c.clear_bond(0, 1).is_ok());
        for a in 0..3u8 {
            for b in (a + 1)..3u8 {
                let planes = BondOrder::ALL
                    .into_iter()
                    .filter(|o| c.neighbours_at(*o, a) & (1u16 << b) != 0)
                    .count();
                assert!(planes <= 1, "pair {a}-{b} is in {planes} planes");
            }
        }
    }

    /// **An index past the atom count must be refused before anything shifts.**
    ///
    /// This is a regression test for a defect that behaved differently in each
    /// profile, which is why it asserts state rather than expecting a panic.
    /// Measured on the code before the fix: `b = 17` panicked with "attempt to
    /// shift left with overflow" under `overflow-checks`, and in release the
    /// shift masked to `b % 16`, so `bond_order(0, 17)` answered `Some(SINGLE)`
    /// for a bond that did not exist and `clear_bond(0, 17)` cleared one row of
    /// a bonded pair — leaving `degree(0) = 0` against `degree(1) = 1`,
    /// `is_connected()` false for a connected dimer, and a different
    /// `CanonForm`. A test that only expected a panic would pass in debug and
    /// assert nothing in the profile that mints goldens.
    #[test]
    fn an_index_past_the_atom_count_shifts_nothing() {
        let t = table();
        let two = valence_exactly(&t, 2);
        let mut m = Mol12::new();
        atom(&mut m, two);
        atom(&mut m, two);
        bond(&mut m, 0, 1, BondOrder::SINGLE, &t);
        let before = m;

        // Past `n`, past `MAX_ATOMS`, and past the row width — the last is the
        // one that used to wrap rather than fail.
        for bad in [2u8, 12, 16, 17, 200, u8::MAX] {
            assert_eq!(
                m.bond_order(0, bad),
                None,
                "bond_order invented a bond at {bad}"
            );
            assert_eq!(
                m.bond_order(bad, 0),
                None,
                "bond_order invented a bond at {bad}"
            );
            assert!(
                matches!(m.clear_bond(0, bad), Err(BondError::NoSuchAtom { .. })),
                "clear_bond accepted {bad}"
            );
            assert_eq!(m, before, "clear_bond({bad}) mutated the molecule");
        }

        // The real bond is still there and still symmetric.
        assert_eq!(m.bond_order(0, 1), Some(BondOrder::SINGLE));
        assert_eq!(m.degree(0), m.degree(1));
        assert!(m.is_connected());
    }

    /// The requirement Task 5 routed here: `BondEnergyMatrix::energy` will price
    /// any order in `1..=BondOrder::MAX` for any pair, including orders the pair
    /// cannot physically support. Nothing in `borbax-universe` can refuse them.
    #[test]
    fn an_order_above_the_shared_valence_is_refused() {
        let t = table();
        let one = valence_exactly(&t, 1);
        let four = valence_exactly(&t, 4);
        let mut m = Mol12::new();
        atom(&mut m, one);
        atom(&mut m, four);
        // min(1, 4) == 1, so a double bond is unpayable at the valence-1 end
        // even though the other end has slots to spare.
        assert!(matches!(
            m.add_bond(0, 1, order(2), &t),
            Err(BondError::ValenceExceeded { atom: 0, .. })
        ));
        // ...and the legal order at the same pair is accepted, so the test is
        // not passing because everything is refused.
        assert!(m.add_bond(0, 1, BondOrder::SINGLE, &t).is_ok());
    }

    /// **The case `min(valence_a, valence_b)` cannot see, and the reason the
    /// check is per-atom rather than per-bond.**
    ///
    /// Every bond here is order 1 between two valence-1 elements, so every one
    /// of them satisfies `order <= min(v_a, v_b)`. The hub still runs out of
    /// slots after the first. At twelve atoms a pairwise check admits eleven
    /// times the hub's budget; this asserts the second bond is already refused.
    #[test]
    fn the_per_atom_budget_outlives_the_pairwise_check() {
        let t = table();
        let one = valence_exactly(&t, 1);
        let mut m = Mol12::new();
        for _ in 0..3 {
            atom(&mut m, one);
        }
        assert!(m.add_bond(0, 1, BondOrder::SINGLE, &t).is_ok());
        assert!(
            matches!(
                m.add_bond(0, 2, BondOrder::SINGLE, &t),
                Err(BondError::ValenceExceeded { .. })
            ),
            "the hub's one slot was spent by the first bond"
        );
    }

    /// Replacing a bond must free the slots the old one held, or a molecule
    /// becomes progressively unbondable by re-stating bonds it already has.
    #[test]
    fn replacing_a_bond_returns_its_slots_first() {
        let t = table();
        let two = valence_exactly(&t, 2);
        let mut m = Mol12::new();
        atom(&mut m, two);
        atom(&mut m, two);
        bond(&mut m, 0, 1, BondOrder::SINGLE, &t);
        // Both ends have two slots and one is spent. Upgrading the *same* pair
        // to a double is legal only if the single's slot is returned first.
        assert!(m.add_bond(0, 1, order(2), &t).is_ok());
        assert_eq!(
            m.bond_order(0, 1),
            Some(order(2)),
            "the orders did not accumulate"
        );
        assert_eq!(m.valence_used(0), 2);
    }

    #[test]
    fn a_bond_needs_two_distinct_atoms_that_exist() {
        let t = table();
        let mut m = chain(&t, 3);
        assert!(matches!(
            m.add_bond(0, 0, BondOrder::SINGLE, &t),
            Err(BondError::NoSuchAtom { .. })
        ));
        assert!(matches!(
            m.add_bond(0, 7, BondOrder::SINGLE, &t),
            Err(BondError::NoSuchAtom { .. })
        ));
    }

    /// `N_ORDERS` is the plane count and is indexed at `order - 1`, so it must
    /// track `BondOrder::MAX` rather than `BondOrder::ALL.len()`. Those coincide
    /// only while `ALL` is dense, which a `const _` inside `borbax-universe`
    /// guarantees and this crate can neither see nor cite.
    #[test]
    fn every_representable_order_has_a_plane() {
        assert_eq!(N_ORDERS, usize::from(BondOrder::MAX));
        for o in BondOrder::ALL {
            assert!(
                o.plane_index() < N_ORDERS,
                "order {o:?} indexes past the planes"
            );
        }
    }
}
