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
//!
//! **Spec sections, since a review asked which one governs.** No single §
//! does — this is test scaffolding, not modelled behaviour, and citing one to
//! satisfy a convention would be worse than saying so. What the fixtures *are*
//! bound by:
//!
//! - `fixture` and `chain_capable` produce elements and a table, so §7.1's
//!   element properties (valence in particular — `chain_capable` filters to
//!   `valence >= 2`) and §7.2's universe generation define what a valid draw
//!   is. The one-seed rule above is a §13.4 determinism consequence.
//! - The molecule builders serve §8.2's signature and §8.3's binding tests, and
//!   `MAX_ATOMS` is §8.1's graph bound.
//! - The rule this module actually enforces is CLAUDE.md's, not the PRD's: a
//!   measured figure quoted against "the corpus" must name one corpus.

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
              of the mistake rather than silently produce a different corpus.\
              \
              `expect` rather than `allow`, deliberately and against one review's advice. The \
              objection is that removing the last `.unwrap()` here breaks the build with \
              `unfulfilled_lint_expectation` — which is true, verified, and the point: it is a \
              LOUD failure telling you a blanket opt-in has outlived its need, and this \
              project's doctrine is that loud beats silent everywhere else too. It is also the \
              crate's convention, 17 `expect` sites to 4 `allow`, and the four are file-scope \
              opt-ins over shipped code rather than test modules."
)]

use crate::canonical::{CanonMol, canonicalise};
use crate::graph::{MAX_ATOMS, Mol12};
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
    // `ids` is `chain_capable`'s output, which already filters to elements
    // present in the table with `valence >= 2` — re-asserting that here would
    // be checking the constructor's postcondition once per fixture. What is
    // *not* guaranteed by construction is that it found any, and an empty list
    // silently yields an empty corpus.
    assert!(
        !ids.is_empty(),
        "random_tree was handed an empty element list, so the corpus it builds is empty too"
    );
    assert!(
        n >= 1 && usize::from(n) <= MAX_ATOMS,
        "random_tree({n}) is outside Mol12's capacity of {MAX_ATOMS} — it would silently \
         return a smaller molecule and every size-dependent assertion downstream would be \
         weaker than it reads"
    );
    let mut mol = Mol12::new();
    assert!(mol.add_atom(ids[0]).is_some());
    for i in 1..n {
        let pick = usize::try_from(rng.next_range(u64::try_from(ids.len()).unwrap())).unwrap();
        let Some(child) = mol.add_atom(ids[pick]) else {
            unreachable!("capacity was checked above, so add_atom cannot fail at {i} of {n}")
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
///
/// **Ring-bearing only for `n > 2`, and callers legitimately pass 1 and 2.**
/// `signature.rs` sweeps `atoms = 1..12` through here. A review proposed
/// asserting `n >= 3`; that would break those callers, and the honest fix is
/// this sentence instead. Do not write a ring-dependent assertion over a corpus
/// that includes the small sizes without filtering for them — which is the same
/// mistake as a corpus filter that silently selects nothing, arrived at from
/// the other side.
#[must_use]
pub(crate) fn symmetric_molecule(n: u8, tbl: &PeriodicTable, one: ElementId) -> Mol12 {
    assert!(
        usize::from(n) <= MAX_ATOMS,
        "symmetric_molecule({n}) is outside Mol12's capacity of {MAX_ATOMS} — it would silently build a \
         smaller fixture than asked"
    );
    let mut mol = Mol12::new();
    for _ in 0..n {
        assert!(mol.add_atom(one).is_some(), "capacity was checked above");
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
/// **The extra bonds can number zero** — `rng.next_range(3)` includes it — so a
/// given draw may be a bare chain with no ring at all. Nothing currently claims
/// otherwise; `layout.rs`'s corpus uses this for automorphism groups and states
/// separately that it does not reach the coincidence branch.
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
    assert!(
        usize::from(n) <= MAX_ATOMS,
        "symmetric_tangle({n}) is outside Mol12's capacity of {MAX_ATOMS} — it would silently build a \
         smaller fixture than asked"
    );
    let mut mol = Mol12::new();
    for _ in 0..n {
        assert!(mol.add_atom(one).is_some(), "capacity was checked above");
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
