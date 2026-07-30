//! Canonical labelling, so isomorphic molecules intern to the same species.
//!
//! Colour refinement (1-WL) followed by individualization-refinement, taking the
//! lexicographically smallest encoding.
//!
//! **Completeness comes from `encode` being faithful, not from the search
//! repairing 1-WL.** An earlier version of this paragraph implied the latter.
//! `CanonForm` records the atom count, the elements in canonical order, and
//! every bonded pair at every order — it *is* the molecule, relabelled — so
//! `canon(G) == canon(H)` means G under one permutation equals H under another,
//! hence G is isomorphic to H, however weak refinement happens to be. 1-WL's
//! incompleteness costs *leaves*, not correctness. What the search is for is
//! making the leaf *set* isomorphism-invariant, so that the minimum over it is
//! a canonical form rather than merely a deterministic one.
//!
//! Verified exhaustively rather than argued: across six enumerations totalling
//! 122 251 valence-legal labelled molecules, zero collisions and zero invariance
//! failures — with the checking machinery independently validated by
//! reproducing 156 unlabelled graphs on six vertices (OEIS A000088).
//!
//! 1-WL weakness is nonetheless real inside the atom cap: `K(3,3)` and the
//! triangular prism are both 3-regular on six vertices, both collapse to one
//! colour class, and only the search separates them. Drop it and species
//! identity breaks at six atoms.
//!
//! The governing cost relation, derived and measured by review:
//! **`#leaves = |Aut(G)| x (number of distinct leaf encodings)`**. Valence
//! bounds degree; degree does **not** bound `|Aut|` — see [`SEARCH_LEAF_CAP`],
//! whose justification used to claim otherwise.

use crate::graph::{MAX_ATOMS, Mol12, N_ORDERS};
use borbax_universe::{BondOrder, ElementId};

/// Upper-triangle pair count for twelve vertices: `12 * 11 / 2 = 66`, which fits
/// one `u128` per bond order with room to spare.
const N_PAIRS: usize = MAX_ATOMS * (MAX_ATOMS - 1) / 2;
const _: () = assert!(N_PAIRS <= 128, "a plane must fit one u128");

/// How many distinct colours can exist at once. **The `+ 1` is headroom, not a
/// fix — this is written down because I first claimed the opposite.**
///
/// The worry was that [`search`] individualizes by *incrementing* every colour
/// above the target, so a colour could reach `n` and overrun an array sized on
/// `MAX_ATOMS`. Probed by setting this to `MAX_ATOMS`: the whole suite passes,
/// including the twelve-cycle written specifically to branch deeply.
///
/// The tight bound holds, and it is worth stating because it is not obvious.
/// Individualization runs **only when `classes < n`** — the discrete case
/// returns before reaching it. Colours are dense ranks, so the maximum before
/// incrementing is `classes - 1 <= n - 2`, and after incrementing it is at most
/// `n - 1 <= 11`. That indexes `[_; MAX_ATOMS]` legally at every level.
///
/// The extra slot stays because the argument above depends on the discrete-exit
/// condition in `search`, which is three lines away and not obviously coupled to
/// an array size over here. It costs six bytes per signature row.
const MAX_COLOURS: usize = MAX_ATOMS + 1;

/// A vertex signature: its own colour, then a neighbour count per
/// (bond order, colour) pair.
const SIG_LEN: usize = 1 + N_ORDERS * MAX_COLOURS;

/// Ceiling on leaves explored by the individualization search.
///
/// **This bounds search cost. It does not encode a claim about chemistry, and
/// an earlier version of this comment did.** It said valence caps degree in any
/// chemically reachable molecule, so hitting the cap "means a molecule was
/// constructed that chemistry cannot produce". That is false, and the
/// counterexample is small: K(6,6) over a valence-6 element — twelve atoms, six
/// bonds each — is **connected, valence-legal, and accepted by `add_bond`
/// without a single refusal**, and it explores 50 000 leaves without finishing.
/// Measured on seed 0, which has a valence-6 element; 22.0% of the first 500
/// seeds do. `|Aut(K₆,₆)| = 2·(6!)² = 1 036 800`, twenty times this cap.
///
/// So the cap is reachable by legal input, and [`canonicalise`] reports it as an
/// error rather than asserting it away. If it is ever hit by a molecule the
/// simulation actually produces, bind to nauty (spec §18.2) rather than growing
/// this search — that route is still right, only its trigger was misdescribed.
pub const SEARCH_LEAF_CAP: u32 = 50_000;

/// The search ran out of budget before it could finish.
///
/// **An unfinished search yields a minimum over a prefix of the leaf set, and
/// nothing establishes that it is the minimum over all of it** — so this is an
/// error rather than a flag on the result. It was a `cap_hit: bool` on
/// [`SearchStats`], documented as "sticky, so callers cannot miss it". Missing
/// it cost one character: `canonicalise(&m).0` compiled clean at `-D warnings`,
/// because `#[must_use]` on the function is discharged by touching the tuple.
/// And the `debug_assert!` meant to catch it fired on **legal** input in the dev
/// profile — the one CLAUDE.md designates for beaker work — blaming a cause that
/// had been falsified.
///
/// **What is NOT established, stated because the previous version asserted it.**
/// That doc said a partial form "interns as a second species". A review went
/// looking for that and could not produce it: 200 relabellings of K(6,6), 200 of
/// 3×K₄ and 12 of twelve isolated atoms each yielded exactly **one** canonical
/// form. The structural reason offered — that exceeding this many leaves at
/// twelve atoms forces a graph homogeneous enough for the first root branch's
/// subtree minimum to already be the global one — is on paper and untested, and
/// explicitly not expected to survive a rise in [`MAX_ATOMS`].
///
/// So the honest claim is the weaker one: a capped search has not *proved* it
/// found the minimum. That is reason enough to refuse to return it, and it is
/// less than the previous doc claimed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error(
    "canonicalisation explored {leaves} leaves without finishing, so the minimum it found is over \
     a prefix of the leaf set and is not established to be the minimum over all of it"
)]
pub struct Capped {
    /// Leaves explored before the search stopped: [`SEARCH_LEAF_CAP`] on the
    /// budget path.
    ///
    /// **Not "always the cap" — an earlier version of this line said so and the
    /// type does not agree.** [`canonicalise`]'s `ok_or` fallback constructs one
    /// with whatever count it has, which on that path is `0`, and the `Display`
    /// message then reads "explored 0 leaves without finishing". The path is
    /// unreachable — `classes <= n` always, and `classes < n` forces a
    /// non-singleton cell by pigeonhole, so a target always exists — but a doc
    /// that denies what the type admits is the shape this file has already
    /// corrected twice.
    ///
    /// Do **not** close it by seeding `Search::best` with the identity ordering:
    /// the identity is not isomorphism-invariant, so a minimum over
    /// `{identity} ∪ leaves` is not a canonical form.
    pub leaves: u32,
}

/// What the search did, so the branch rate can be measured rather than assumed.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SearchStats {
    /// Discrete colourings reached, i.e. candidate labellings encoded.
    ///
    /// For a molecule whose refinement does not discretise, this converges on
    /// the size of the automorphism group — the twelve-cycle reaches exactly 24,
    /// which is `|Aut(C₁₂)|`, the dihedral group of order 24.
    pub leaves: u32,
    /// Total refinement rounds across every node of the search tree.
    pub refinement_rounds: u32,
}

/// A molecule's canonical encoding. Comparison is a few integer compares, which
/// is what makes this usable as an intern key.
///
/// Field order defines the derived [`Ord`], and therefore defines which
/// labelling counts as "smallest". **It is part of the physics**: changing it
/// changes every species id and every golden hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonForm {
    /// Atom count. Slots past this in `elem` are [`ElementId::ZERO`] padding.
    pub n: u8,
    /// Elements in canonical position order.
    pub elem: [ElementId; MAX_ATOMS],
    /// One bitset of canonical pairs per bond order, indexed by
    /// `BondOrder::plane_index`.
    pub planes: [u128; N_ORDERS],
}

/// Canonical index for the unordered pair `{i, j}`.
#[inline]
fn pair_bit(i: usize, j: usize) -> u32 {
    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
    u32::try_from(hi * (hi - 1) / 2 + lo).unwrap_or(0)
}

/// Encode a molecule under a given vertex ordering.
/// `order[k]` is the original vertex placed at canonical position `k`.
fn encode(m: &Mol12, order: &[u8]) -> CanonForm {
    let mut elem = [ElementId::ZERO; MAX_ATOMS];
    let mut planes = [0u128; N_ORDERS];
    for (k, &v) in order.iter().enumerate() {
        if let (Some(slot), Some(id)) = (elem.get_mut(k), m.element(v)) {
            *slot = id;
        }
    }
    for a in 0..order.len() {
        for b in (a + 1)..order.len() {
            let (Some(&va), Some(&vb)) = (order.get(a), order.get(b)) else {
                continue;
            };
            // `if let Some(o)`, not `if o > 0` — there is no no-bond sentinel.
            if let Some(o) = m.bond_order(va, vb)
                && let Some(plane) = planes.get_mut(o.plane_index())
            {
                *plane |= 1u128 << pair_bit(a, b);
            }
        }
    }
    CanonForm {
        n: u8::try_from(m.len()).unwrap_or(0),
        elem,
        planes,
    }
}

/// One round-loop of colour refinement. Returns the number of colour classes.
///
/// A vertex's signature under the current colouring is its own colour plus, for
/// every (bond order, colour) pair, how many neighbours it has of that colour —
/// all computed as popcounts over colour masks. No allocation, no hashing, and
/// no map iteration whose order could vary (§13.1).
#[expect(
    clippy::indexing_slicing,
    reason = "every index here is proven in range: `v` and `w` run over `0..n` with \
              `n <= MAX_ATOMS`, and every colour is < MAX_COLOURS because `refine` \
              emits ranks in `0..n` and `search` raises them by at most one. The \
              `k` cursor is bounded by `N_ORDERS * MAX_COLOURS`, which is SIG_LEN - 1. \
              debug_assert below pins the colour bound rather than leaving it prose"
)]
fn refine(m: &Mol12, colour: &mut [u8; MAX_ATOMS], stats: &mut SearchStats) -> usize {
    let n = m.len();
    if n == 0 {
        return 0;
    }
    loop {
        stats.refinement_rounds += 1;
        debug_assert!(
            colour[..n].iter().all(|&c| usize::from(c) < MAX_COLOURS),
            "a colour escaped the array these masks are sized for"
        );
        let n_colours = usize::from(colour[..n].iter().max().copied().unwrap_or(0)) + 1;

        let mut masks = [0u16; MAX_COLOURS];
        for v in 0..n {
            masks[usize::from(colour[v])] |= 1u16 << v;
        }

        let mut sig = [[0u8; SIG_LEN]; MAX_ATOMS];
        for v in 0..n {
            // **This line is what makes the termination test sound**, and it
            // reads like an ordinary part of the signature. Because the vertex's
            // own colour is the leading component of a lexicographic array
            // compare, the new partition always *refines* the old, so an equal
            // class count forces an unchanged partition — which is what lets the
            // loop below stop on the count rather than on the partition.
            // Delete it and refinement stops refining: the class count never
            // stabilises, the recursion never bottoms out, and the observed
            // failure is a stack overflow wearing a costume. Confirmed by review
            // asserting monotone refinement over the whole suite plus 120k
            // random molecules: never fired.
            sig[v][0] = colour[v];
            let mut k = 1;
            for o in BondOrder::ALL {
                let row = m.neighbours_at(o, u8::try_from(v).unwrap_or(0));
                for mask in masks.iter().take(n_colours) {
                    sig[v][k] = u8::try_from((row & mask).count_ones()).unwrap_or(0);
                    k += 1;
                }
            }
        }

        // Rank signatures to produce new colours. Insertion-order-free: the rank
        // of a signature is how many *distinct* signatures are strictly smaller,
        // which depends only on the multiset of signatures.
        //
        // **"Is this the first occurrence of `sig[w]`" is a function of `w`
        // alone**, so it is computed once per `w` rather than once per `(v, w)`.
        // The predecessor evaluated it inside the double loop, making the round
        // O(n³·L) in 79-byte array compares where O(n²·L) does the same work.
        // Measured, bit-identical: **−36% on a realistic corpus and −47% on the
        // twelve-cycle**, against a p5 noise floor of 1.7% and 0.7%.
        //
        // The booleans are the same because `&&` is commutative over pure
        // operands — but "the suite passes" would not establish that, since the
        // relabelling properties only check that two labellings agree with *each
        // other* and would hold under a genuinely different refinement. Verified
        // instead by a digest over canonical forms across four alphabet sizes
        // (100, 3, 2 and 1 elements — the small ones branch 28% to 100% of the
        // time) plus the twelve-cycle, byte-identical either way.
        let mut first_occurrence = [false; MAX_ATOMS];
        for w in 0..n {
            first_occurrence[w] = !(0..w).any(|x| sig[x] == sig[w]);
        }
        let mut new_colour = [0u8; MAX_ATOMS];
        for v in 0..n {
            let mut rank = 0u8;
            for w in 0..n {
                if first_occurrence[w] && sig[w] < sig[v] {
                    rank += 1;
                }
            }
            new_colour[v] = rank;
        }

        let before = count_classes(colour, n);
        let after = count_classes(&new_colour, n);
        *colour = new_colour;
        if after == before {
            return after;
        }
    }
}

fn count_classes(colour: &[u8; MAX_ATOMS], n: usize) -> usize {
    let mut seen = 0u32;
    for &c in colour.iter().take(n) {
        seen |= 1u32 << c;
    }
    usize::try_from(seen.count_ones()).unwrap_or(0)
}

/// Reduce a molecule to its canonical form.
///
/// Returns the form plus search statistics, which exist so the branch rate can
/// be measured rather than assumed.
///
/// # Errors
///
/// [`Capped`] if the search exhausted [`SEARCH_LEAF_CAP`] leaves. That is
/// reachable from legal input — see the constant — so it is a value the caller
/// must handle rather than a condition asserted away.
pub fn canonicalise(m: &Mol12) -> Result<(CanonForm, SearchStats), Capped> {
    let n = m.len();
    if n == 0 {
        return Ok((
            CanonForm {
                n: 0,
                elem: [ElementId::ZERO; MAX_ATOMS],
                planes: [0; N_ORDERS],
            },
            SearchStats::default(),
        ));
    }

    // Initial colouring by element, **densified**.
    //
    // A colour indexes `masks` and shifts into `count_classes`'s bitset, both of
    // which are sized for `MAX_COLOURS` — while a table holds 60 to 120 elements.
    // So the colour cannot be the element id; it has to be the id's *rank* among
    // the distinct ids present, which is dense in `0..n`.
    //
    // The plan warned that assigning the raw id "fails silently". Worth recording
    // that it cannot be assigned at all from here: `ElementId`'s field and its
    // `index()` are both `pub(crate)` in `borbax-universe`, so there is no way to
    // get an integer out of one in this crate. The type closed that hazard; what
    // is left for this code to get right is that the rank depends only on the
    // *multiset* of ids present, which is what keeps it isomorphism-invariant —
    // and `relabelling_is_stable_across_the_full_element_range` is what checks it.
    let elems: Vec<ElementId> = (0..u8::try_from(n).unwrap_or(0))
        .filter_map(|i| m.element(i))
        .collect();
    let mut colour = [0u8; MAX_ATOMS];
    for (slot, ev) in colour.iter_mut().zip(elems.iter()) {
        let mut rank = 0u8;
        for (w, ew) in elems.iter().enumerate() {
            if ew < ev && !elems.iter().take(w).any(|seen| seen == ew) {
                rank += 1;
            }
        }
        *slot = rank;
    }

    let mut s = Search::default();
    search(m, colour, &mut s);
    if s.capped {
        return Err(Capped {
            leaves: s.stats.leaves,
        });
    }
    // For n > 0 an uncapped search always reaches a leaf: refinement either
    // discretises, or the target cell is non-empty by pigeonhole and every
    // child is explored, so depth is bounded by n.
    s.best.map(|form| (form, s.stats)).ok_or(Capped {
        leaves: s.stats.leaves,
    })
}

/// The search's working state. Bundled so `capped` cannot be dropped on the way
/// out the way `SearchStats::cap_hit` could be dropped by the caller.
#[derive(Default)]
struct Search {
    best: Option<CanonForm>,
    stats: SearchStats,
    capped: bool,
}

#[expect(
    clippy::indexing_slicing,
    reason = "three sites, needing two arguments. `colour[v]` twice: `v` runs over \
              `0..n <= MAX_ATOMS`, the array's length. `order[colour[v]]` is the one \
              that is not obvious and an earlier reason omitted it while calling \
              `order` a colour array — it is a `Vec<u8>` of length `n`, the inverse \
              permutation, indexed by a colour VALUE. That is in range only because \
              this branch runs when `classes == n` and `refine` emits dense ranks, so \
              `colour[..n]` is a permutation of `0..n`. There is no `w` in this \
              function; that was `refine`'s justification pasted onto a different \
              index shape"
)]
fn search(m: &Mol12, mut colour: [u8; MAX_ATOMS], s: &mut Search) {
    let n = m.len();
    if s.capped || s.stats.leaves >= SEARCH_LEAF_CAP {
        // **No `debug_assert!(false)` here, and its removal is the fix.** It
        // fired on K(6,6) over a valence-6 element — connected, valence-legal,
        // every one of its 36 bonds accepted — with a message asserting the
        // molecule was not valence-legal. That is a panic on legal input, in the
        // profile CLAUDE.md designates for beaker work, in a workspace where
        // `clippy::panic` is denied and `debug_assert!` routes around the deny.
        // The condition is now carried out as a value instead.
        s.capped = true;
        return;
    }
    let classes = refine(m, &mut colour, &mut s.stats);

    if classes == n {
        // Discrete: the colouring *is* an ordering.
        s.stats.leaves += 1;
        let mut order = vec![0u8; n];
        for v in 0..n {
            order[usize::from(colour[v])] = u8::try_from(v).unwrap_or(0);
        }
        let form = encode(m, &order);
        if s.best.as_ref().is_none_or(|b| form < *b) {
            s.best = Some(form);
        }
        return;
    }

    // Target cell: **the lowest colour VALUE with more than one member** — not
    // the smallest cell by cardinality.
    //
    // The distinction is a trap, and the previous wording walked into it. To a
    // graph-theory reader "smallest cell" means smallest *cardinality*, which is
    // an equally standard and equally isomorphism-invariant selector. This code
    // implements nauty's default, the first non-singleton cell by colour value.
    //
    // A review replicated the algorithm and swapped only the selector: **the two
    // give different canonical forms** — measured on the Petersen graph and on
    // two 5-rings sharing a bridging atom, an ordinary valence-legal 11-atom
    // molecule. Leaf counts are identical, so nothing looks wrong, and *both*
    // selectors pass every relabelling property here because both are invariant.
    // So somebody "correcting" the code to match a comment that said cardinality
    // would silently change every species id and every golden hash with the
    // whole suite green.
    //
    // The selector is physics in exactly the sense `CanonForm`'s field order is,
    // and unlike that one it was unlabelled. `the_target_cell_rule_is_pinned`
    // holds it to a value, because no invariance test can.
    let target = (0..MAX_COLOURS)
        .find(|&c| {
            colour
                .iter()
                .take(n)
                .filter(|&&v| usize::from(v) == c)
                .count()
                > 1
        })
        .and_then(|c| u8::try_from(c).ok());
    let Some(target) = target else {
        return;
    };

    for v in 0..n {
        if colour[v] != target {
            continue;
        }
        // Individualize v: it keeps the target colour, everyone else in the cell
        // (and every higher colour) shifts up by one.
        let mut c2 = colour;
        for (w, slot) in c2.iter_mut().enumerate().take(n) {
            if *slot > target || (*slot == target && w != v) {
                *slot += 1;
            }
        }
        search(m, c2, s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{MAX_ATOMS, Mol12};
    use borbax_rng::{Domain, Stream};
    use borbax_universe::element::generate_elements;
    use borbax_universe::{BondOrder, ElementId, PeriodicTable};

    fn table(seed: u64) -> PeriodicTable {
        generate_elements(seed)
    }

    /// Every element that can carry a chain interior.
    ///
    /// **Valence ≥ 2 is what makes the spanning tree below always completable**,
    /// and the argument is a slot count rather than a hope: a tree over `n` nodes
    /// spends `2(n - 1)` slots against a budget of at least `2n`, so at least two
    /// slots are always free and at least one earlier node always has one. Without
    /// that, `random_molecule` could emit isolated atoms — and twelve identical
    /// isolated atoms have a 12! automorphism group, which would blow
    /// `SEARCH_LEAF_CAP` for reasons that say nothing about the code under test.
    fn tree_capable(t: &PeriodicTable) -> Vec<ElementId> {
        t.iter()
            .filter(|(_, e)| e.valence >= 2)
            .map(|(id, _)| id)
            .collect()
    }

    fn draw(rng: &mut Stream, n: usize) -> usize {
        let bound = u64::try_from(n).unwrap_or_else(|_| unreachable!("a length fits u64"));
        usize::try_from(rng.next_range(bound))
            .unwrap_or_else(|_| unreachable!("next_range returns < n"))
    }

    fn draw_u8(rng: &mut Stream, n: u8) -> u8 {
        u8::try_from(rng.next_range(u64::from(n)))
            .unwrap_or_else(|_| unreachable!("next_range returns < n <= u8::MAX"))
    }

    fn nth<T: Copy>(xs: &[T], i: usize) -> T {
        xs.get(i)
            .copied()
            .unwrap_or_else(|| unreachable!("index was drawn in range"))
    }

    /// Canonicalise, asserting the search finished.
    ///
    /// Every fixture here is a molecule whose search completes; a `Capped` would
    /// mean the fixture changed into something else, so it fails loudly rather
    /// than being absorbed.
    fn canon(m: &Mol12) -> (CanonForm, SearchStats) {
        canonicalise(m).unwrap_or_else(|e| unreachable!("fixture capped: {e}"))
    }

    fn as_u8(i: usize) -> u8 {
        u8::try_from(i).unwrap_or_else(|_| unreachable!("bounded by MAX_ATOMS"))
    }

    /// Slots still unspent at atom `p`.
    fn free_slots(m: &Mol12, t: &PeriodicTable, p: u8) -> u32 {
        m.element(p).and_then(|id| t.get(id)).map_or(0, |e| {
            u32::from(e.valence).saturating_sub(m.valence_used(p))
        })
    }

    /// A random connected, valence-legal molecule.
    ///
    /// **The parent is chosen before the order, and the order is capped by the
    /// budget that actually remains.** An earlier version drew the order first
    /// and searched for a parent that could pay it, on the argument that a tree
    /// over `n` nodes spends `2(n-1)` slots against a budget of at least `2n`.
    /// That argument holds only for single bonds and the generator was drawing
    /// from all six: two valence-2 atoms joined by a double bond are both full,
    /// so the third atom has no payable parent at any order. It fired on the
    /// first run, which is the only reason the fixture is not quietly emitting
    /// disconnected graphs.
    ///
    /// Choosing the parent from those with a free slot, then taking
    /// `min(drawn, free(parent), valence(new))`, is payable by construction —
    /// every factor is at least one — and still reaches the high orders, because
    /// an early atom with a large valence and no bonds yet can pay for them.
    fn random_molecule(rng: &mut Stream, table: &PeriodicTable, ids: &[ElementId]) -> Mol12 {
        let mut mol = Mol12::new();
        let target = 2 + draw_u8(rng, as_u8(MAX_ATOMS) - 1);
        assert!(mol.add_atom(nth(ids, draw(rng, ids.len()))).is_some());

        for i in 1..target {
            let candidates: Vec<u8> = (0..i)
                .filter(|&p| free_slots(&mol, table, p) >= 1)
                .collect();
            if candidates.is_empty() {
                break; // Fully saturated: stop growing rather than disconnect.
            }
            let parent = nth(&candidates, draw(rng, candidates.len()));
            let elem = nth(ids, draw(rng, ids.len()));
            let Some(child) = mol.add_atom(elem) else {
                break;
            };

            // **Draw from `BondOrder::ALL`, not `1 + next_range(3)`.** The literal
            // 3 is the stale constant `N_ORDERS` just stopped being: with it,
            // planes 3, 4 and 5 are never set anywhere in the suite, so `encode`'s
            // writes to them and `refine`'s popcounts over them go untested while
            // both property tests report green.
            let drawn = nth(&BondOrder::ALL, draw(rng, BondOrder::ALL.len()));
            let head = table.get(elem).map_or(1, |el| u32::from(el.valence));
            let payable = u32::from(u8::from(drawn))
                .min(free_slots(&mol, table, parent))
                .min(head);
            let ord = u8::try_from(payable)
                .ok()
                .and_then(BondOrder::new)
                .unwrap_or(BondOrder::SINGLE);
            assert!(
                mol.add_bond(parent, child, ord, table).is_ok(),
                "the order was capped by both ends' remaining budget"
            );
        }

        // A few extra bonds, which is where rings and branches come from. These
        // may legitimately be refused, and a refused bond is not a mutation —
        // `a_refused_bond_leaves_the_molecule_untouched` is what licenses this.
        let size = as_u8(mol.len());
        for _ in 0..rng.next_range(4) {
            let from = draw_u8(rng, size);
            let to = draw_u8(rng, size);
            let ord = nth(&BondOrder::ALL, draw(rng, BondOrder::ALL.len()));
            let _ = mol.add_bond(from, to, ord, table);
        }
        mol
    }

    /// Relabel a molecule by a permutation. An isomorphic graph must reduce to
    /// the same canonical form — that is the whole contract.
    fn permute(m: &Mol12, p: &[u8], t: &PeriodicTable) -> Mol12 {
        let mut out = Mol12::new();
        for i in 0..m.len() {
            let e = m
                .element(nth(p, i))
                .unwrap_or_else(|| unreachable!("p permutes 0..n"));
            assert!(out.add_atom(e).is_some());
        }
        for i in 0..m.len() {
            for j in (i + 1)..m.len() {
                // `if let Some(o)`, not `if o > 0`. `bond_order` has no
                // no-bond sentinel to feed back into `add_bond`.
                if let Some(o) = m.bond_order(nth(p, i), nth(p, j)) {
                    assert!(
                        out.add_bond(as_u8(i), as_u8(j), o, t).is_ok(),
                        "a relabelling of a legal molecule is legal"
                    );
                }
            }
        }
        out
    }

    fn shuffle(rng: &mut Stream, n: u8) -> Vec<u8> {
        let mut p: Vec<u8> = (0..n).collect();
        for i in (1..p.len()).rev() {
            let j = draw(rng, i + 1);
            p.swap(i, j);
        }
        p
    }

    /// Element ids must span the real range, not `0..6`.
    ///
    /// The colouring bug this catches is invisible with a handful of low ids:
    /// colouring by raw element id indexes `refine`'s `[u16; MAX_ATOMS]` mask
    /// array at up to 119, and `count_classes` computes `1u16 << colour[v]` —
    /// which panics in debug and **masks the shift** in release, so refinement
    /// stops before it is stable and two relabellings of one molecule intern as
    /// different species.
    #[test]
    fn relabelling_is_stable_across_the_full_element_range() {
        let t = table(17);
        let ids = tree_capable(&t);
        // Stated as an assertion, because the test is worthless if the ids it
        // draws never leave the range the colour arrays are sized for.
        let past_the_colour_array =
            ElementId::from_index(MAX_ATOMS).unwrap_or_else(|| unreachable!("12 fits u8"));
        assert!(
            ids.iter().any(|id| *id >= past_the_colour_array),
            "the densification bug is only reachable with ids past MAX_ATOMS"
        );

        let mut rng = Stream::new(17, Domain::Molecule, 0);
        for _ in 0..2_000 {
            let m = random_molecule(&mut rng, &t, &ids);
            let (want, _) = canon(&m);
            let p = shuffle(&mut rng, as_u8(m.len()));
            assert_eq!(canon(&permute(&m, &p, &t)).0, want);
        }
    }

    #[test]
    fn relabelling_does_not_change_the_canonical_form() {
        let t = table(1);
        let ids = tree_capable(&t);
        let mut rng = Stream::new(1, Domain::Molecule, 0);
        for _ in 0..2_000 {
            let m = random_molecule(&mut rng, &t, &ids);
            let want = canon(&m).0;
            for _ in 0..4 {
                let p = shuffle(&mut rng, as_u8(m.len()));
                assert_eq!(
                    canon(&permute(&m, &p, &t)).0,
                    want,
                    "relabelling changed the form"
                );
            }
        }
    }

    /// **The measurement that replaces a `proptest` strategy.**
    ///
    /// The plan proposed `prop::sample::select(BondOrder::ALL)` so the generator
    /// tracks the array rather than a literal. Drawing from `ALL` does that
    /// already; what neither approach gives on its own is evidence that the high
    /// planes were ever *written*. This asserts it directly, so "planes 3, 4 and
    /// 5 are never set anywhere in the suite" is a failing test rather than a
    /// review finding.
    ///
    /// **The bar is the universe's own ceiling, not `BondOrder::MAX`.** Task 5
    /// measured order 6 to be reachable at all in **23.2%** of universes, so
    /// asserting all six planes would demand chemistry three quarters of seeds
    /// forbid — a test that fails for being right about the physics. The ceiling
    /// here is the largest valence in the drawn table, which is attainable
    /// because the generator may draw the same element at both ends.
    #[test]
    fn every_plane_the_universe_admits_is_exercised() {
        let t = table(17);
        let ids = tree_capable(&t);
        let ceiling = ids
            .iter()
            .filter_map(|id| t.get(*id))
            .map(|e| e.valence)
            .max()
            .unwrap_or(1)
            .min(as_u8(N_ORDERS));
        assert!(ceiling >= 2, "seed 17 admits no bond order above a single");

        let mut rng = Stream::new(17, Domain::Molecule, 0);
        let mut seen = [false; N_ORDERS];
        for _ in 0..2_000 {
            let m = random_molecule(&mut rng, &t, &ids);
            for i in 0..m.len() {
                for j in (i + 1)..m.len() {
                    if let Some(o) = m.bond_order(as_u8(i), as_u8(j))
                        && let Some(slot) = seen.get_mut(o.plane_index())
                    {
                        *slot = true;
                    }
                }
            }
        }
        let reached = seen.iter().take(usize::from(ceiling)).all(|s| *s);
        assert!(
            reached,
            "seed 17 admits orders 1..={ceiling} but the suite wrote {seen:?}"
        );
    }

    #[test]
    fn different_molecules_get_different_forms() {
        let t = table(5);
        let e = tree_capable(&t);
        let two = nth(&e, 0);
        let mut a = Mol12::new();
        assert!(a.add_atom(two).is_some());
        assert!(a.add_atom(two).is_some());
        assert!(a.add_bond(0, 1, BondOrder::SINGLE, &t).is_ok());

        let mut b = Mol12::new();
        assert!(b.add_atom(two).is_some());
        assert!(b.add_atom(two).is_some());
        let double = BondOrder::new(2).unwrap_or_else(|| unreachable!("2 is an order"));
        assert!(
            b.add_bond(0, 1, double, &t).is_ok(),
            "same atoms, different order"
        );

        assert_ne!(canon(&a).0, canon(&b).0);
    }

    #[test]
    fn element_labels_distinguish_otherwise_identical_graphs() {
        let t = table(5);
        let e = tree_capable(&t);
        let (first, second) = (nth(&e, 0), nth(&e, 1));

        let mut a = Mol12::new();
        assert!(a.add_atom(first).is_some());
        assert!(a.add_atom(first).is_some());
        assert!(a.add_bond(0, 1, BondOrder::SINGLE, &t).is_ok());

        let mut b = Mol12::new();
        assert!(b.add_atom(first).is_some());
        assert!(b.add_atom(second).is_some());
        assert!(b.add_bond(0, 1, BondOrder::SINGLE, &t).is_ok());

        assert_ne!(canon(&a).0, canon(&b).0);
    }

    /// A ring of identical atoms, which refinement **cannot** discretise.
    ///
    /// **This test exists because the measured branch rate is 0.8% with a worst
    /// case of 2 leaves.** That is the good news the perf review predicted, and
    /// it means the entire individualization search — target-cell selection, the
    /// colour increment, the recursion — is exercised by almost nothing in the
    /// suite. A twelve-cycle of one element is the opposite extreme: every vertex
    /// has identical colour and identical neighbourhood at every round, so 1-WL
    /// returns one class forever and the answer can only come from branching.
    ///
    /// It is also the only test that reaches high colour values, which is what
    /// `MAX_COLOURS` exists for: individualization raises colours by one per
    /// level, and a twelve-cycle nests deeply enough to matter.
    #[test]
    fn a_symmetric_ring_forces_the_search_to_branch() {
        let t = table(5);
        let ids = tree_capable(&t);
        // Valence >= 2 at every atom, since a ring gives each two bonds.
        let e = nth(&ids, 0);
        let mut m = Mol12::new();
        for _ in 0..MAX_ATOMS {
            assert!(m.add_atom(e).is_some());
        }
        for i in 0..as_u8(MAX_ATOMS) {
            let next = (i + 1) % as_u8(MAX_ATOMS);
            assert!(m.add_bond(i, next, BondOrder::SINGLE, &t).is_ok());
        }

        let (want, stats) = canon(&m);
        assert!(
            stats.leaves > 1,
            "refinement cannot discretise a vertex-transitive graph; the search must have run"
        );
        // `#leaves = |Aut(G)| x (distinct leaf encodings)`, derived and measured
        // by review. A twelve-cycle has one encoding and the dihedral group of
        // order 24, so this is an exact expected value rather than a bound — and
        // it is the cheapest available check that the search explores the orbit
        // it should, no more and no less.
        assert_eq!(stats.leaves, 24, "|Aut(C12)| is 24: the dihedral group");

        // And the whole point still holds: every relabelling agrees.
        let mut rng = Stream::new(99, Domain::Molecule, 0);
        for _ in 0..64 {
            let p = shuffle(&mut rng, as_u8(m.len()));
            assert_eq!(canon(&permute(&m, &p, &t)).0, want);
        }
    }

    /// **The target-cell selector is physics, and only a pinned value holds it.**
    ///
    /// `canonical.rs` picks the lowest colour *value* with more than one member.
    /// Choosing the smallest cell by *cardinality* is equally standard and
    /// equally isomorphism-invariant — and gives **different canonical forms**.
    /// Both selectors satisfy every relabelling property in this file, because
    /// both are invariant; invariance is what those tests check and it cannot
    /// distinguish them. So a well-meaning change to match a differently-worded
    /// comment would move every species id with the suite green.
    ///
    /// The Petersen graph is the fixture because it is where a review measured
    /// the two selectors diverging: 10 atoms, 3-regular, girth 5, one element,
    /// all single bonds. Refinement gives every vertex the same colour, so the
    /// answer comes entirely from the search.
    ///
    /// A moved value here is a deliberate physics change requiring golden
    /// regeneration, exactly like `CanonForm`'s field order. It is not a
    /// tolerance to widen.
    #[test]
    fn the_target_cell_rule_is_pinned() {
        let t = table(5);
        let three = t.iter().find(|(_, e)| e.valence >= 3).map_or_else(
            || unreachable!("seed 5 has a valence-3 element"),
            |(id, _)| id,
        );

        let mut m = Mol12::new();
        for _ in 0..10 {
            assert!(m.add_atom(three).is_some());
        }
        // Outer pentagon, spokes, inner pentagram.
        for i in 0..5u8 {
            assert!(m.add_bond(i, (i + 1) % 5, BondOrder::SINGLE, &t).is_ok());
            assert!(m.add_bond(i, i + 5, BondOrder::SINGLE, &t).is_ok());
            assert!(
                m.add_bond(i + 5, ((i + 2) % 5) + 5, BondOrder::SINGLE, &t)
                    .is_ok()
            );
        }
        assert!(m.is_connected());
        for v in 0..10u8 {
            assert_eq!(m.degree(v), 3, "Petersen is 3-regular");
        }

        let (form, stats) = canon(&m);
        // |Aut(Petersen)| = 120, and it has one leaf orbit, so
        // `#leaves = |Aut| x (distinct encodings)` predicts exactly 120.
        assert_eq!(stats.leaves, 120, "|Aut(Petersen)| is 120");
        assert_eq!(form.n, 10);
        assert_eq!(
            form.planes[0], PETERSEN_SINGLE_PLANE,
            "the canonical form moved: either the target-cell rule or the \
             encoding changed, and both are physics"
        );
    }

    /// The single-bond plane of the Petersen graph's canonical form under the
    /// lowest-colour-value target-cell rule. Regenerate deliberately, never to
    /// make a test pass.
    const PETERSEN_SINGLE_PLANE: u128 = 1_315_755_469_568;

    /// The perf review predicted refinement alone discretises for the large
    /// majority of small labelled molecules, leaving the expensive search a
    /// rarely-taken path. Measure it rather than believe it — the number decides
    /// how much this matters later.
    #[test]
    fn the_branch_rate_is_reported_against_the_palette_it_was_measured_on() {
        let t = table(2);
        let all = tree_capable(&t);
        // The palette is the dominant variable, so it is swept rather than
        // fixed. Uniform draws from 100 elements make almost every atom unique,
        // which is why the shipped figure was so low.
        let mut worst_overall = 0u32;
        for palette in [1usize, 2, 3, 6, 12, all.len()] {
            let ids: Vec<ElementId> = all.iter().copied().take(palette).collect();
            let mut rng = Stream::new(2, Domain::Molecule, 0);
            let (mut total, mut branched, mut worst) = (0u32, 0u32, 0u32);
            for _ in 0..2_000 {
                let m = random_molecule(&mut rng, &t, &ids);
                let (_, stats) = canon(&m);
                total += 1;
                if stats.leaves > 1 {
                    branched += 1;
                }
                worst = worst.max(stats.leaves);
            }
            let pct = 100.0 * f64::from(branched) / f64::from(total);
            println!("palette {palette:3}: branched {pct:5.1}%, worst {worst} leaves");
            worst_overall = worst_overall.max(worst);
        }
        // The conclusion the module doc rests on, asserted over every palette
        // rather than over the friendliest one.
        assert!(
            worst_overall < SEARCH_LEAF_CAP,
            "hit the search cap — investigate, do not raise it"
        );
        assert!(
            worst_overall <= 64,
            "worst case grew past anything previously measured: {worst_overall}"
        );
    }
}
