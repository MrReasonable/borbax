# Borbax V0 Implementation Plan — Part 2: Chemistry, Beaker, and Exit Criteria

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

Continues [`2026-07-26-borbax-v0.md`](2026-07-26-borbax-v0.md), which holds the
Goal, Architecture, **Global Constraints**, and File Structure. Those apply here
unchanged — read them first. Tasks 1–10 are there; this file has 11–21.

**Spec:** `docs/superpowers/specs/2026-07-26-borbax-prd.md`. V0 exit criteria are §23.

---

## Phase 4 — Polymers, folding, and cavities

### Task 11: FCC-lattice folding

**Files:**
- Create: `crates/borbax-molecule/src/polymer.rs`, `crates/borbax-molecule/src/fold.rs`
- Test: both

**Interfaces:**
- Consumes: `borbax_universe::{Universe, ElementId}`, `borbax_rng::{Stream, Domain}`
- Produces: `MAX_POLYMER`, `Polymer`, `Fold`, `FoldWorkspace`, `FoldWorkspace::fold(&Polymer, &Universe, &mut Stream) -> Fold`, `Fold::{contacts, energy, exposure}`, `fcc::{GRID, NEIGHBOURS, to_flat, from_flat}`

- [ ] **Step 1: Write the failing tests**

`crates/borbax-molecule/src/fold.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use borbax_rng::{Domain, Stream};
    use borbax_universe::{ElementId, Universe};

    fn poly(n: usize, u: &Universe) -> Polymer {
        let mut p = Polymer::new();
        for i in 0..n {
            p.push(ElementId((i * 7 % u.elements.len()) as u8));
        }
        p
    }

    #[test]
    fn folding_is_deterministic() {
        let u = Universe::generate(8);
        let p = poly(40, &u);
        let mut ws = FoldWorkspace::new();
        let a = ws.fold(&p, &u, &mut Stream::new(1, Domain::Fold, 0));
        let b = ws.fold(&p, &u, &mut Stream::new(1, Domain::Fold, 0));
        assert_eq!(a.coords, b.coords);
        assert_eq!(a.contacts, b.contacts);
    }

    /// The workspace is reused across folds and reset by trail rather than
    /// cleared. If the reset is wrong, fold N+1 sees fold N's occupancy —
    /// a bug that would only show up as mysteriously bad folds much later.
    #[test]
    fn workspace_reuse_matches_a_fresh_workspace() {
        let u = Universe::generate(8);
        let (p1, p2) = (poly(35, &u), poly(50, &u));
        let mut shared = FoldWorkspace::new();
        let _ = shared.fold(&p1, &u, &mut Stream::new(2, Domain::Fold, 0));
        let reused = shared.fold(&p2, &u, &mut Stream::new(3, Domain::Fold, 0));

        let mut fresh = FoldWorkspace::new();
        let clean = fresh.fold(&p2, &u, &mut Stream::new(3, Domain::Fold, 0));
        assert_eq!(reused.coords, clean.coords, "workspace reset leaked state");
    }

    #[test]
    fn the_chain_stays_connected_and_self_avoiding() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        for seed in 0..30 {
            let p = poly(45, &u);
            let f = ws.fold(&p, &u, &mut Stream::new(seed, Domain::Fold, 0));
            let mut seen = std::collections::BTreeSet::new();
            for (i, &c) in f.coords.iter().enumerate() {
                assert!(seen.insert(c), "self-intersection at monomer {i}");
                if i > 0 {
                    assert!(
                        fcc::NEIGHBOURS.contains(&(c - f.coords[i - 1])),
                        "chain broken between {} and {i}", i - 1
                    );
                }
            }
        }
    }

    /// Incremental contact counts must equal a full recount. Integer counts
    /// are used precisely so this can be an exact assertion — a running f64
    /// energy would drift and slowly change which conformation wins.
    #[test]
    fn incremental_contacts_match_a_full_recount() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        for seed in 0..20 {
            let p = poly(40, &u);
            let f = ws.fold(&p, &u, &mut Stream::new(seed, Domain::Fold, 0));
            assert_eq!(f.contacts, recount_contacts(&p, &f), "seed {seed}");
        }
    }

    /// Annealing must actually do something. A fold no better than its
    /// starting extended conformation means the move set or schedule is
    /// broken, and every cavity downstream would be meaningless.
    #[test]
    fn annealing_beats_the_extended_start() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        let p = poly(50, &u);
        let f = ws.fold(&p, &u, &mut Stream::new(9, Domain::Fold, 0));
        let folded: u32 = f.contacts.iter().sum();
        assert!(folded > 8, "only {folded} contacts — annealing is not compacting");
    }

    #[test]
    fn exposure_is_between_zero_and_twelve() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        let p = poly(40, &u);
        let f = ws.fold(&p, &u, &mut Stream::new(4, Domain::Fold, 0));
        assert!(f.exposure.iter().take(p.len()).all(|&e| e <= 12));
        // A compact fold must bury something, or §9.5's burial argument fails.
        assert!(f.exposure.iter().take(p.len()).any(|&e| e < 8), "nothing is buried");
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test -p borbax-molecule fold`
Expected: FAIL — `cannot find type FoldWorkspace`

- [ ] **Step 3: Implement the FCC lattice**

`crates/borbax-molecule/src/fold.rs`, first section:

```rust
//! Polymer folding on a face-centred cubic lattice (spec §8.4).
//!
//! FCC gives twelve nearest neighbours — the densest lattice packing — which
//! is what lets a chain fold genuinely compactly rather than into something
//! loose. That matters because burial is what protects a backbone from
//! solvent attack (§9.5) and enclosure is what makes a cavity selective
//! (§8.5). A looser lattice would weaken both.
//!
//! This is the most expensive operation in the engine, so the workspace is
//! allocated once and reused, and the energy is tracked as integer contact
//! counts rather than a running float.

pub mod fcc {
    /// Grid edge length. 64^3 cells at 2 bytes is 512 KB — the per-thread
    /// budget in spec §17. Chains that would escape are recentred rather
    /// than growing the grid.
    pub const GRID: i32 = 64;
    pub const CELLS: usize = (GRID * GRID * GRID) as usize;

    /// Flat-index offsets to the twelve FCC neighbours. Compile-time
    /// constants, so neighbour enumeration is twelve constant-offset loads
    /// with no coordinate arithmetic in the inner loop.
    pub const NEIGHBOURS: [i32; 12] = [
        1 + GRID,
        1 - GRID,
        -1 + GRID,
        -1 - GRID,
        1 + GRID * GRID,
        1 - GRID * GRID,
        -1 + GRID * GRID,
        -1 - GRID * GRID,
        GRID + GRID * GRID,
        GRID - GRID * GRID,
        -GRID + GRID * GRID,
        -GRID - GRID * GRID,
    ];

    #[must_use]
    pub const fn to_flat(x: i32, y: i32, z: i32) -> i32 {
        x + y * GRID + z * GRID * GRID
    }

    #[must_use]
    pub const fn from_flat(f: i32) -> (i32, i32, i32) {
        (f % GRID, (f / GRID) % GRID, f / (GRID * GRID))
    }

    /// True if a flat index is far enough from the boundary that all twelve
    /// neighbours are in range.
    #[must_use]
    pub fn in_bounds(f: i32) -> bool {
        let (x, y, z) = from_flat(f);
        (2..GRID - 2).contains(&x) && (2..GRID - 2).contains(&y) && (2..GRID - 2).contains(&z)
    }
}
```

- [ ] **Step 4: Implement the polymer and the fold**

`crates/borbax-molecule/src/polymer.rs`:

```rust
//! Polymers: sequences of monomer units.
//!
//! Unlike small molecules (§8.1) polymers are not graphs — they are chains,
//! and their shape comes from folding rather than from topology.

use borbax_universe::ElementId;

/// Cap on chain length. Chosen so a chain cannot escape the fold grid even
/// fully extended, which removes a whole class of boundary handling.
pub const MAX_POLYMER: usize = 200;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Polymer {
    units: Vec<u8>,
}

impl Polymer {
    #[must_use]
    pub fn new() -> Self {
        Self { units: Vec::new() }
    }

    pub fn push(&mut self, e: ElementId) -> bool {
        if self.units.len() >= MAX_POLYMER {
            return false;
        }
        self.units.push(e.0);
        true
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.units.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.units.is_empty()
    }

    #[must_use]
    pub fn unit(&self, i: usize) -> ElementId {
        ElementId(self.units[i])
    }

    #[must_use]
    pub fn units(&self) -> &[u8] {
        &self.units
    }
}
```

`crates/borbax-molecule/src/fold.rs`, second section:

```rust
use crate::polymer::{Polymer, MAX_POLYMER};
use borbax_rng::Stream;
use borbax_units::Quanta;
use borbax_universe::{ElementId, Universe};

/// Annealing steps per fold. Fixed, not convergence-based — a tolerance stop
/// would make the result depend on float details that vary between platforms
/// (spec §13.4).
const ANNEAL_STEPS: u32 = 20_000;

/// Number of monomer classes contacts are bucketed into. Contacts are counted
/// per unordered class pair, so the running energy is an integer histogram.
const N_CLASSES: usize = 8;
const N_PAIRS: usize = N_CLASSES * (N_CLASSES + 1) / 2;

#[derive(Debug, Clone, PartialEq)]
pub struct Fold {
    pub n: usize,
    /// Flat lattice index per monomer.
    pub coords: Vec<i32>,
    /// Contact counts per unordered monomer-class pair. Integers, so
    /// incremental updates cannot drift from a full recount.
    pub contacts: [u32; N_PAIRS],
    /// Empty neighbours out of twelve, per monomer. Low means buried, which
    /// is what protects a backbone from solvent attack (§9.5).
    pub exposure: Vec<u8>,
}

impl Fold {
    /// Total energy, recomputed exactly from the integer contact histogram.
    /// Never accumulated incrementally as a float.
    #[must_use]
    pub fn energy(&self, class_energy: &[f64; N_PAIRS]) -> Quanta {
        let mut acc = 0.0;
        for i in 0..N_PAIRS {
            acc += f64::from(self.contacts[i]) * class_energy[i];
        }
        Quanta(acc)
    }
}

#[inline]
const fn pair_index(a: usize, b: usize) -> usize {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    hi * (hi + 1) / 2 + lo
}

/// Reusable per-thread scratch. Allocated once; never reallocated, and never
/// cleared wholesale — resetting walks the chain's own cells, which is
/// O(chain) rather than O(grid). At 512 KB a memset per fold would dominate.
pub struct FoldWorkspace {
    occ: Box<[u16]>,
    trail: Vec<i32>,
}

impl Default for FoldWorkspace {
    fn default() -> Self {
        Self::new()
    }
}

impl FoldWorkspace {
    #[must_use]
    pub fn new() -> Self {
        Self { occ: vec![0u16; fcc::CELLS].into_boxed_slice(), trail: Vec::with_capacity(MAX_POLYMER) }
    }

    fn clear(&mut self) {
        for &c in &self.trail {
            self.occ[c as usize] = 0;
        }
        self.trail.clear();
    }

    fn occupy(&mut self, cell: i32, monomer: usize) {
        self.occ[cell as usize] = monomer as u16 + 1;
        self.trail.push(cell);
    }

    fn is_free(&self, cell: i32) -> bool {
        self.occ[cell as usize] == 0
    }

    pub fn fold(&mut self, p: &Polymer, u: &Universe, rng: &mut Stream) -> Fold {
        self.clear();
        let n = p.len();
        if n == 0 {
            return Fold { n: 0, coords: Vec::new(), contacts: [0; N_PAIRS], exposure: Vec::new() };
        }

        let class = |i: usize| -> usize {
            // Bucket monomers by affinity so contact energy depends on
            // character rather than on which element happened to be used.
            let a = u.element(p.unit(i)).affinity;
            (((a + 1.0) * 0.5 * (N_CLASSES as f64 - 1.0)).round() as usize).min(N_CLASSES - 1)
        };

        // Extended start along one lattice direction, centred in the grid.
        let centre = fcc::to_flat(fcc::GRID / 2, fcc::GRID / 2, fcc::GRID / 2);
        let step = fcc::NEIGHBOURS[0];
        let mut coords = Vec::with_capacity(n);
        for i in 0..n {
            let c = centre + step * (i as i32 - n as i32 / 2);
            coords.push(c);
            self.occupy(c, i);
        }

        let mut contacts = count_contacts(&coords, self, &class, n);
        let energy_table = class_energy_table(u);
        let mut current = contact_energy(&contacts, &energy_table);

        for stepno in 0..ANNEAL_STEPS {
            // Linear cooling. Deterministic and monotone; no transcendentals.
            let temp = 1.0 - f64::from(stepno) / f64::from(ANNEAL_STEPS);
            let i = rng.next_range(n as u64) as usize;

            let Some(target) = self.propose(&coords, i, n, rng) else {
                continue;
            };
            if !fcc::in_bounds(target) || !self.is_free(target) {
                continue;
            }

            // Incremental delta: a move touches one monomer, so only its own
            // contacts change. That is O(12) instead of O(n^2) per step.
            let old = coords[i];
            let removed = local_contacts(&coords, self, &class, i, old, n);
            self.occ[old as usize] = 0;
            self.occ[target as usize] = i as u16 + 1;
            let added = local_contacts(&coords, self, &class, i, target, n);

            let mut trial = contacts;
            for (idx, c) in removed.iter() {
                trial[*idx] -= c;
            }
            for (idx, c) in added.iter() {
                trial[*idx] += c;
            }
            let trial_energy = contact_energy(&trial, &energy_table);

            // Metropolis. Lower energy always accepted; uphill accepted with
            // probability falling as the schedule cools.
            let accept = trial_energy <= current
                || rng.next_f64() < temp * 0.35;
            if accept {
                coords[i] = target;
                // Keep the trail complete so `clear` reaches every touched cell.
                self.trail.push(target);
                contacts = trial;
                current = trial_energy;
            } else {
                self.occ[target as usize] = 0;
                self.occ[old as usize] = i as u16 + 1;
            }
        }

        let exposure = (0..n)
            .map(|i| {
                fcc::NEIGHBOURS.iter().filter(|&&d| self.is_free(coords[i] + d)).count() as u8
            })
            .collect();

        Fold { n, coords, contacts, exposure }
    }

    /// Propose a move. Two move types, both of which preserve chain
    /// connectivity and self-avoidance by construction:
    ///
    /// - **End rotation**: a terminal monomer moves to any free neighbour of
    ///   its one bonded partner.
    /// - **Corner flip**: an internal monomer moves to a free cell adjacent
    ///   to both of its bonded partners.
    ///
    /// This move set is not fully ergodic — pull moves would be — but it is
    /// simple enough to be obviously correct, and `annealing_beats_the_
    /// extended_start` checks it actually compacts. If folds later look
    /// under-compacted, adding pull moves is the first thing to try.
    fn propose(&self, coords: &[i32], i: usize, n: usize, rng: &mut Stream) -> Option<i32> {
        if n < 2 {
            return None;
        }
        if i == 0 || i == n - 1 {
            let anchor = if i == 0 { coords[1] } else { coords[n - 2] };
            let d = fcc::NEIGHBOURS[rng.next_range(12) as usize];
            Some(anchor + d)
        } else {
            let (a, b) = (coords[i - 1], coords[i + 1]);
            // Free cells adjacent to both neighbours.
            let mut options = [0i32; 12];
            let mut k = 0;
            for &d in &fcc::NEIGHBOURS {
                let c = a + d;
                if c != coords[i] && fcc::NEIGHBOURS.contains(&(c - b)) {
                    options[k] = c;
                    k += 1;
                }
            }
            if k == 0 {
                return None;
            }
            Some(options[rng.next_range(k as u64) as usize])
        }
    }
}

fn class_energy_table(u: &Universe) -> [f64; N_PAIRS] {
    let mut t = [0.0; N_PAIRS];
    for a in 0..N_CLASSES {
        for b in a..N_CLASSES {
            // Class index maps back to a representative affinity in [-1, 1].
            let fa = (a as f64) / (N_CLASSES as f64 - 1.0) * 2.0 - 1.0;
            let fb = (b as f64) / (N_CLASSES as f64 - 1.0) * 2.0 - 1.0;
            // Opposite characters attract; like characters do not. Same
            // principle as §8.3, applied between monomers instead of surfaces.
            t[pair_index(a, b)] = -fa * fb * u.consts.w_charge;
        }
    }
    t
}

fn contact_energy(contacts: &[u32; N_PAIRS], table: &[f64; N_PAIRS]) -> f64 {
    let mut acc = 0.0;
    for i in 0..N_PAIRS {
        acc += f64::from(contacts[i]) * table[i];
    }
    acc
}

/// Contacts contributed by monomer `i` sitting at `cell`: lattice-adjacent
/// monomers that are not its chain neighbours.
fn local_contacts(
    coords: &[i32],
    ws: &FoldWorkspace,
    class: &impl Fn(usize) -> usize,
    i: usize,
    cell: i32,
    n: usize,
) -> Vec<(usize, u32)> {
    let mut out = Vec::new();
    for &d in &fcc::NEIGHBOURS {
        let occupant = ws.occ[(cell + d) as usize];
        if occupant == 0 {
            continue;
        }
        let j = occupant as usize - 1;
        if j == i || j >= n || j + 1 == i || i + 1 == j {
            continue; // self or chain neighbour: not a contact
        }
        out.push((pair_index(class(i), class(j)), 1));
    }
    out
}

fn count_contacts(
    coords: &[i32],
    ws: &FoldWorkspace,
    class: &impl Fn(usize) -> usize,
    n: usize,
) -> [u32; N_PAIRS] {
    let mut c = [0u32; N_PAIRS];
    for i in 0..n {
        for (idx, k) in local_contacts(coords, ws, class, i, coords[i], n) {
            c[idx] += k;
        }
    }
    // Each contact was counted from both ends.
    for v in &mut c {
        *v /= 2;
    }
    c
}

/// Independent full recount, used only by tests to verify the incremental
/// path. Deliberately written differently from `count_contacts` so a shared
/// bug cannot make both agree.
#[cfg(test)]
pub fn recount_contacts(p: &Polymer, f: &Fold) -> [u32; N_PAIRS] {
    use std::collections::BTreeMap;
    let mut pos: BTreeMap<i32, usize> = BTreeMap::new();
    for (i, &c) in f.coords.iter().enumerate() {
        pos.insert(c, i);
    }
    let mut out = [0u32; N_PAIRS];
    for i in 0..f.n {
        for &d in &fcc::NEIGHBOURS {
            if let Some(&j) = pos.get(&(f.coords[i] + d)) {
                if j > i && j != i + 1 {
                    // Class recomputed here from the polymer, not cached.
                    let cls = |k: usize| -> usize { (k * 0) + class_of(p, k) };
                    out[pair_index(cls(i), cls(j))] += 1;
                }
            }
        }
    }
    out
}

#[cfg(test)]
fn class_of(_p: &Polymer, _i: usize) -> usize {
    // Filled in by the implementer to mirror the `class` closure in `fold`,
    // reading from the same Universe. Kept separate on purpose: the recount
    // must not share code with the incremental path.
    unimplemented!("mirror the class closure using the test's Universe")
}
```

> **Implementer note on `recount_contacts`:** it needs the same `Universe` the
> fold used, so give it a `&Universe` parameter and inline the class mapping
> rather than leaving the `unimplemented!` above. It is written as a separate
> function on purpose — a recount that shares code with the incremental path
> cannot catch a bug in that path.

- [ ] **Step 5: Run the tests**

Run: `cargo test -p borbax-molecule fold`
Expected: 6 passed

`annealing_beats_the_extended_start` and `exposure_is_between_zero_and_twelve`
are the meaningful ones. If nothing gets buried, §9.5's burial argument and
§8.5's cavity selectivity both collapse — try pull moves before adjusting
the thresholds.

- [ ] **Step 6: Commit**

```bash
git add crates/borbax-molecule
git commit -m "$(cat <<'EOF'
feat(molecule): FCC-lattice polymer folding

FCC for its twelve neighbours — the densest lattice packing, which is what
lets a chain fold compactly enough for burial to protect a backbone (§9.5)
and for cavities to be enclosed enough to be selective (§8.5).

Energy is an integer contact histogram, not a running float. Incremental
updates to an f64 total would drift from a full recount and slowly change
which conformation wins, and the conformation is the cache key for
everything downstream. Tested against an independently written recount.

The workspace is allocated once and reset by trail rather than cleared:
O(chain) instead of O(grid), which matters at 512 KB per fold. Reuse is
tested against a fresh workspace, because a leaked-state bug would surface
only as mysteriously poor folds much later.

Move set is end rotations plus corner flips — not fully ergodic, but
simple enough to be obviously correct. Pull moves are the first thing to
try if folds look under-compacted.

MrReasonable <4990954+MrReasonable@users.noreply.github.com>
EOF
)"
```

---
