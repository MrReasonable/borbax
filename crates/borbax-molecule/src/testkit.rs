//! Shared fixtures for this crate's tests.
//!
//! **Why this module exists.** `random_tree`, `fixture`, `canon` and friends
//! were defined independently in `layout.rs`, `signature.rs` and `binding.rs`,
//! and by the third copy they had already drifted — `Mol12::new()` against
//! `Mol12::default()`, and one copy carrying eight lines of reasoning the
//! others had lost. Cosmetic that time. The hazard is not:
//!
//! - A measured figure is quoted in a doc comment against "the corpus". Three
//!   corpora answering to one name means the number and the reader disagree,
//!   silently, which is the stale-figure class this project has recorded
//!   repeatedly.
//! - The `random_tree` in `layout.rs` carries a *measured* justification for
//!   its retry-against-other-parents branch — "0 of 10,020 trees" left a child
//!   isolated. A copy without that comment is edited blind.
//!
//! **`fixture(seed)` returning both the table and the universe from one seed is
//! load-bearing**, not a convenience. `layout.rs` records what the mismatch
//! cost: pairing `table(17)` with `Universe::generate(4)` left 30.8% of ids
//! absent from the universe, silently taking a fallback radius, and inflated a
//! composition-only null to within 0.008 of the gate it was supposed to clear.
//! Taking one seed makes the chimera unspellable.

#![expect(
    clippy::redundant_pub_crate,
    reason = "`pub(crate)` is what these are: this module is private and `#[cfg(test)]`, so \
              the visibility that reaches the sibling test modules and no further is exactly \
              crate. `redundant_pub_crate` sees the enclosing module is private and calls it \
              redundant; `unreachable_pub` rejects the bare `pub` it suggests. Only one of \
              the two can be satisfied, and this spelling is the one that states the intent."
)]
#![expect(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    reason = "CLAUDE.md permits both inside test code with a stated reason. This module is \
              `#[cfg(test)]` in `lib.rs` but is a module rather than a `mod tests` block, so \
              it needs its own opt-in — the deny reaches here exactly as it would inside one. \
              Every index is a loop bound over a count already established, and every unwrap \
              is on a fixture the caller controls: a broken fixture should panic at the point \
              of the mistake rather than silently produce a different corpus."
)]

use crate::canonical::{CanonMol, canonicalise};
use crate::graph::Mol12;
use borbax_rng::Stream;
use borbax_universe::{BondOrder, ElementId, PeriodicTable, Universe, element::generate_elements};

/// A periodic table and a universe **from the same seed**.
#[must_use]
pub(crate) fn fixture(seed: u64) -> (PeriodicTable, Universe) {
    (generate_elements(seed), Universe::generate(seed))
}

/// Elements that can carry two bonds, so a chain is buildable.
#[must_use]
pub(crate) fn chain_capable(tbl: &PeriodicTable) -> Vec<ElementId> {
    (0..120usize)
        .filter_map(ElementId::from_index)
        .filter(|id| tbl.get(*id).is_some_and(|el| el.valence >= 2))
        .collect()
}

/// Canonicalise, or fail loudly naming the cause.
#[must_use]
pub(crate) fn canon(mol: &Mol12) -> CanonMol {
    canonicalise(mol)
        .unwrap_or_else(|err| unreachable!("fixture capped: {err}"))
        .0
}

/// A random tree, built by attaching each new atom to an existing one.
///
/// **The retry branch is measured, not defensive.** A refused bond would leave
/// the child isolated and the result would not be a tree, so it retries against
/// *other parents* — all at `SINGLE`, since there is no lower order to retry
/// down to. If every parent refuses, the child is left isolated and this
/// returns a non-tree: measured unreachable, 0 of 10,020 trees, and
/// structurally so, since `chain_capable` elements carry valence >= 2 and a
/// tree always retains spare capacity.
#[must_use]
pub(crate) fn random_tree(
    rng: &mut Stream,
    n: u8,
    tbl: &PeriodicTable,
    ids: &[ElementId],
) -> Mol12 {
    let mut mol = Mol12::new();
    assert!(mol.add_atom(ids[0]).is_some());
    for i in 1..n {
        let pick = usize::try_from(rng.next_range(u64::try_from(ids.len()).unwrap())).unwrap();
        let Some(child) = mol.add_atom(ids[pick]) else {
            break;
        };
        let parent = u8::try_from(rng.next_range(u64::from(i))).unwrap();
        if mol.add_bond(parent, child, BondOrder::SINGLE, tbl).is_err() {
            for other in 0..i {
                if mol.add_bond(other, child, BondOrder::SINGLE, tbl).is_ok() {
                    break;
                }
            }
        }
    }
    mol
}

/// A random molecule: a spanning tree plus a few extra edges, so rings and
/// branches appear rather than trees only.
#[must_use]
pub(crate) fn random_molecule(
    rng: &mut Stream,
    n: u8,
    tbl: &PeriodicTable,
    ids: &[ElementId],
) -> Mol12 {
    let mut mol = random_tree(rng, n, tbl, ids);
    let size = u8::try_from(mol.len()).unwrap();
    for _ in 0..rng.next_range(4) {
        let a = u8::try_from(rng.next_range(u64::from(size))).unwrap();
        let b = u8::try_from(rng.next_range(u64::from(size))).unwrap();
        // A refused bond is not a mutation; the fixture is whatever survives.
        let _ = mol.add_bond(a, b, BondOrder::SINGLE, tbl);
    }
    mol
}

/// A ring-bearing molecule built from **one element throughout**, which
/// maximises the automorphism group and is what reaches `embed`'s coincidence
/// branch — automorphic same-element atoms are the ones with identical target
/// rows. Deterministic given `n`, so it takes no `Stream`.
#[must_use]
pub(crate) fn symmetric_molecule(n: u8, tbl: &PeriodicTable, one: ElementId) -> Mol12 {
    let mut mol = Mol12::new();
    for _ in 0..n {
        if mol.add_atom(one).is_none() {
            break;
        }
    }
    let size = u8::try_from(mol.len()).unwrap_or(0);
    for i in 1..size {
        let _ = mol.add_bond(i - 1, i, BondOrder::SINGLE, tbl);
    }
    if size > 2 {
        let _ = mol.add_bond(0, size - 1, BondOrder::SINGLE, tbl);
    }
    mol
}

/// A one-element molecule with **randomly** placed extra bonds, on top of the
/// chain.
///
/// **Deliberately a different name from [`symmetric_molecule`], because it is a
/// different function.** Both were called `symmetric_molecule` in different
/// test modules and both documented themselves as "maximally symmetric" — one
/// deterministic with a single ring closure, one drawing 0..3 extra bonds from
/// a stream. Two corpora answering to one identifier in one crate is how a
/// measured figure quoted against "the symmetric corpus" comes to mean two
/// things, which is the drift this module exists to stop.
#[must_use]
pub(crate) fn symmetric_tangle(
    rng: &mut Stream,
    n: u8,
    tbl: &PeriodicTable,
    one: ElementId,
) -> Mol12 {
    let mut mol = Mol12::new();
    for _ in 0..n {
        if mol.add_atom(one).is_none() {
            break;
        }
    }
    let size = u8::try_from(mol.len()).unwrap_or(0);
    for i in 1..size {
        let _ = mol.add_bond(i - 1, i, BondOrder::SINGLE, tbl);
    }
    for _ in 0..rng.next_range(3) {
        let a = u8::try_from(rng.next_range(u64::from(size))).unwrap_or(0);
        let b = u8::try_from(rng.next_range(u64::from(size))).unwrap_or(0);
        let _ = mol.add_bond(a, b, BondOrder::SINGLE, tbl);
    }
    mol
}
