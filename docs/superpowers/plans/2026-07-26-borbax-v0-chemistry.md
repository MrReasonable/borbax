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

### Task 12: The fold cache

**Files:**
- Create: `crates/borbax-molecule/src/foldcache.rs`
- Test: same file

**Interfaces:**
- Consumes: `Polymer`, `Fold`, `FoldWorkspace`, `Universe`, `Stream`
- Produces: `FoldCache::{with_budget, get_or_fold, stats}`, `CacheStats { hits, misses, evictions, bytes }`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use borbax_rng::{Domain, Stream};
    use borbax_universe::{ElementId, Universe};

    fn poly(seed: u64, n: usize) -> Polymer {
        let mut p = Polymer::new();
        let mut r = Stream::new(seed, Domain::Fold, 99);
        for _ in 0..n {
            p.push(ElementId(r.next_range(20) as u8));
        }
        p
    }

    #[test]
    fn a_hit_returns_exactly_what_a_miss_computed() {
        let u = Universe::generate(2);
        let mut c = FoldCache::with_budget(8 << 20);
        let p = poly(1, 30);
        let first = c.get_or_fold(&p, &u).clone();
        let second = c.get_or_fold(&p, &u).clone();
        assert_eq!(first, second);
        assert_eq!(c.stats().hits, 1);
        assert_eq!(c.stats().misses, 1);
    }

    #[test]
    fn distinct_sequences_do_not_collide() {
        let u = Universe::generate(2);
        let mut c = FoldCache::with_budget(64 << 20);
        for i in 0..500u64 {
            let _ = c.get_or_fold(&poly(i, 24), &u);
        }
        assert_eq!(c.stats().misses, 500, "two distinct sequences shared an entry");
    }

    /// Eviction must be driven by the entries' own measured sizes at fixed
    /// step boundaries — never by available system memory, which would make
    /// results depend on what else the machine was doing (spec §13.1).
    #[test]
    fn eviction_respects_the_declared_budget() {
        let u = Universe::generate(2);
        let mut c = FoldCache::with_budget(64 << 10); // deliberately tiny
        for i in 0..400u64 {
            let _ = c.get_or_fold(&poly(i, 40), &u);
        }
        assert!(c.stats().bytes <= 64 << 10, "budget exceeded: {}", c.stats().bytes);
        assert!(c.stats().evictions > 0, "budget never bit");
    }

    #[test]
    fn eviction_is_deterministic() {
        let u = Universe::generate(2);
        let run = || {
            let mut c = FoldCache::with_budget(64 << 10);
            for i in 0..300u64 {
                let _ = c.get_or_fold(&poly(i, 40), &u);
            }
            (c.stats().hits, c.stats().misses, c.stats().evictions)
        };
        assert_eq!(run(), run());
    }
}
```

- [ ] **Step 2: Run to verify failure** — `cargo test -p borbax-molecule foldcache`, expect `cannot find type FoldCache`

- [ ] **Step 3: Implement**

```rust
//! Cache of folded conformations (spec §17).
//!
//! A fold is a pure function of the sequence, so caching cannot change any
//! result — unlike a lossy cache, this is safe under §13.1 no matter what it
//! evicts. What *would* be unsafe is letting eviction depend on anything
//! outside the simulation, so the budget is evaluated from the entries' own
//! measured sizes and never from available system memory.

use crate::fold::{Fold, FoldWorkspace};
use crate::polymer::Polymer;
use borbax_rng::{Domain, Stream};
use borbax_universe::Universe;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// Pass-through hasher: keys are already a strong 128-bit hash, so hashing
/// them again would be pure cost. `RandomState` is banned outright — it is
/// seeded nondeterministically.
#[derive(Default)]
pub struct PassThrough(u64);

impl Hasher for PassThrough {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut b = [0u8; 8];
            b[..chunk.len()].copy_from_slice(chunk);
            self.0 ^= u64::from_le_bytes(b);
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CacheStats {
    pub hits: u64,
    pub misses: u64,
    pub evictions: u64,
    pub bytes: usize,
}

struct Entry {
    fold: Fold,
    bytes: usize,
    /// CLOCK reference bit. A single hand index serialises trivially into a
    /// keyframe, which an LRU linked list would not — and CLOCK avoids an
    /// allocation and a pointer chase per access.
    referenced: bool,
}

pub struct FoldCache {
    map: HashMap<u128, usize, BuildHasherDefault<PassThrough>>,
    slots: Vec<Option<(u128, Entry)>>,
    free: Vec<usize>,
    hand: usize,
    budget: usize,
    stats: CacheStats,
    ws: FoldWorkspace,
}

/// 128-bit key from two independent passes of the RNG mixer over the
/// sequence. At a few million entries the birthday probability is ~1e-27, so
/// storing the sequence itself would be paying for a guarantee we already have.
fn key_of(p: &Polymer) -> u128 {
    let mut a = Stream::new(0xB0_1BAA_5EED, Domain::Fold, 1);
    let mut b = Stream::new(0x5EED_B0_1BAA, Domain::Fold, 2);
    let (mut ha, mut hb) = (a.next_u64(), b.next_u64());
    for (i, &u) in p.units().iter().enumerate() {
        ha ^= u64::from(u).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left((i % 61) as u32);
        ha = ha.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        hb ^= u64::from(u).wrapping_add(i as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93);
        hb = hb.rotate_left(29).wrapping_mul(0x94D0_49BB_1331_11EB);
    }
    (u128::from(ha) << 64) | u128::from(hb)
}

impl FoldCache {
    #[must_use]
    pub fn with_budget(bytes: usize) -> Self {
        Self {
            map: HashMap::default(),
            slots: Vec::new(),
            free: Vec::new(),
            hand: 0,
            budget: bytes,
            stats: CacheStats::default(),
            ws: FoldWorkspace::new(),
        }
    }

    #[must_use]
    pub fn stats(&self) -> CacheStats {
        self.stats
    }

    pub fn get_or_fold(&mut self, p: &Polymer, u: &Universe) -> &Fold {
        let k = key_of(p);
        if let Some(&slot) = self.map.get(&k) {
            self.stats.hits += 1;
            if let Some((_, e)) = &mut self.slots[slot] {
                e.referenced = true;
            }
            // Re-borrow immutably for the return.
            return &self.slots[slot].as_ref().map(|(_, e)| &e.fold).unwrap_or_else(|| unreachable!());
        }

        self.stats.misses += 1;
        // Fold seed derives from the key, so the same sequence anneals
        // identically regardless of when it is first encountered.
        let mut rng = Stream::new(k as u64, Domain::Fold, 0);
        let fold = self.ws.fold(p, u, &mut rng);
        let bytes = std::mem::size_of::<Fold>()
            + fold.coords.len() * std::mem::size_of::<i32>()
            + fold.exposure.len();

        while self.stats.bytes + bytes > self.budget && !self.map.is_empty() {
            self.evict_one();
        }

        let slot = self.free.pop().unwrap_or_else(|| {
            self.slots.push(None);
            self.slots.len() - 1
        });
        self.stats.bytes += bytes;
        self.slots[slot] = Some((k, Entry { fold, bytes, referenced: false }));
        self.map.insert(k, slot);
        &self.slots[slot].as_ref().map(|(_, e)| &e.fold).unwrap_or_else(|| unreachable!())
    }

    /// CLOCK: sweep, clearing reference bits, evict the first unreferenced
    /// entry. Deterministic given a deterministic access sequence, and the
    /// entire eviction state is one integer.
    fn evict_one(&mut self) {
        loop {
            if self.slots.is_empty() {
                return;
            }
            self.hand = (self.hand + 1) % self.slots.len();
            let idx = self.hand;
            let Some((k, e)) = &mut self.slots[idx] else { continue };
            if e.referenced {
                e.referenced = false;
                continue;
            }
            let (k, e) = (*k, self.slots[idx].take().map(|(_, e)| e).unwrap_or_else(|| unreachable!()));
            self.map.remove(&k);
            self.stats.bytes -= e.bytes;
            self.stats.evictions += 1;
            self.free.push(idx);
            return;
        }
    }
}
```

> **Implementer note:** the `unwrap_or_else(|| unreachable!())` forms above are
> placeholders for the borrow shape only — restructure with an index-then-index
> pattern or `split_at_mut` so no `unreachable!` survives into the committed
> code. The Global Constraints ban `unwrap`/`expect` in library code and this
> is the same class of thing.

- [ ] **Step 4: Run** — `cargo test -p borbax-molecule foldcache`, expect 4 passed

- [ ] **Step 5: Commit**

```bash
git add crates/borbax-molecule/src/foldcache.rs
git commit -m "$(cat <<'EOF'
feat(molecule): fold cache with CLOCK eviction

A fold is a pure function of its sequence, so caching cannot change a
result no matter what it evicts — unlike a lossy cache, this is safe under
§13.1 unconditionally. What would be unsafe is eviction depending on
anything outside the simulation, so the budget is measured from the
entries' own sizes and never from available system memory.

Keys are 128-bit hashes of the sequence rather than the sequence itself:
at a few million entries the birthday probability is around 1e-27, so
storing the sequence buys a guarantee we already have. RandomState is
banned outright since it is seeded nondeterministically.

CLOCK rather than LRU: the whole eviction state is one hand index, which
serialises into a keyframe trivially where a linked list would not, and it
avoids an allocation and a pointer chase per access.

MrReasonable <4990954+MrReasonable@users.noreply.github.com>
EOF
)"
```

---

### Task 13: Cavity detection

**Files:**
- Create: `crates/borbax-molecule/src/cavity.rs`
- Test: same file

**Interfaces:**
- Consumes: `Fold`, `Polymer`, `Geodesic<D>`, `Universe`, `fcc`
- Produces: `Cavity { cells, centre, enclosure, signature }`, `cavities<D>(&Fold, &Polymer, &Geodesic<D>, &Universe) -> Vec<Cavity<D>>`, `MIN_ENCLOSURE`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{fold::FoldWorkspace, geodesic::Geodesic, polymer::Polymer};
    use borbax_rng::{Domain, Stream};
    use borbax_universe::{ElementId, Universe};

    fn folded(seed: u64, n: usize, u: &Universe) -> (Polymer, crate::fold::Fold) {
        let mut p = Polymer::new();
        let mut r = Stream::new(seed, Domain::Fold, 7);
        for _ in 0..n {
            p.push(ElementId(r.next_range(20) as u8));
        }
        let f = FoldWorkspace::new().fold(&p, u, &mut Stream::new(seed, Domain::Fold, 0));
        (p, f)
    }

    #[test]
    fn detection_is_deterministic() {
        let (u, g) = (Universe::generate(12), Geodesic::<42>::build().unwrap());
        let (p, f) = folded(3, 60, &u);
        assert_eq!(cavities(&f, &p, &g, &u), cavities(&f, &p, &g, &u));
    }

    #[test]
    fn cavities_are_enclosed_by_construction() {
        let (u, g) = (Universe::generate(12), Geodesic::<42>::build().unwrap());
        let (p, f) = folded(4, 70, &u);
        for c in cavities(&f, &p, &g, &u) {
            assert!(c.enclosure >= MIN_ENCLOSURE, "reported an unenclosed cavity");
        }
    }

    #[test]
    fn cavity_cells_are_empty_and_connected() {
        let (u, g) = (Universe::generate(12), Geodesic::<42>::build().unwrap());
        let (p, f) = folded(5, 70, &u);
        let occupied: std::collections::BTreeSet<i32> = f.coords.iter().copied().collect();
        for c in cavities(&f, &p, &g, &u) {
            assert!(c.cells.iter().all(|x| !occupied.contains(x)), "cavity overlaps the chain");
            // Connectivity: every cell reachable from the first.
            let set: std::collections::BTreeSet<i32> = c.cells.iter().copied().collect();
            let mut seen = std::collections::BTreeSet::new();
            let mut stack = vec![c.cells[0]];
            while let Some(x) = stack.pop() {
                if !seen.insert(x) {
                    continue;
                }
                for &d in &crate::fold::fcc::NEIGHBOURS {
                    if set.contains(&(x + d)) {
                        stack.push(x + d);
                    }
                }
            }
            assert_eq!(seen.len(), c.cells.len(), "cavity is not connected");
        }
    }

    /// Longer chains fold more compactly and should produce more cavities.
    /// If they do not, either folding is not compacting or the enclosure
    /// threshold is wrong — and §8.5's catalysis has nothing to work with.
    #[test]
    fn longer_chains_produce_more_cavities() {
        let (u, g) = (Universe::generate(12), Geodesic::<42>::build().unwrap());
        let short: usize = (0..12).map(|s| { let (p, f) = folded(s, 25, &u); cavities(&f, &p, &g, &u).len() }).sum();
        let long: usize = (0..12).map(|s| { let (p, f) = folded(s, 80, &u); cavities(&f, &p, &g, &u).len() }).sum();
        assert!(long > short, "80-mers ({long}) gave no more cavities than 25-mers ({short})");
    }
}
```

- [ ] **Step 2: Run to verify failure** — expect `cannot find function cavities`

- [ ] **Step 3: Implement**

```rust
//! Cavity detection on a folded chain (spec §8.4, §8.5).
//!
//! A cavity is a connected pocket of empty lattice cells substantially
//! enclosed by the chain. This is the three-dimensional payoff: in 2D the
//! equivalent is a notch on a perimeter, open on two sides and unselective.
//! Here a cavity can be surrounded on many sides at once, which is what makes
//! a binding site *specific* rather than merely sticky — and specificity is
//! the difference between a catalyst and a patch of glue.

use crate::fold::{fcc, Fold};
use crate::geodesic::Geodesic;
use crate::polymer::Polymer;
use crate::signature::Signature;
use borbax_universe::Universe;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Minimum occupied-neighbour count for an empty cell to count as enclosed.
/// Out of twelve FCC neighbours, eight means surrounded on two-thirds of its
/// faces. Below this a "cavity" is a surface dimple, and §8.5's selectivity
/// argument does not hold for it.
pub const MIN_ENCLOSURE: u8 = 8;

/// Largest cavity worth reporting. Anything bigger is the outside world.
const MAX_CAVITY_CELLS: usize = 24;

#[derive(Debug, Clone, PartialEq)]
pub struct Cavity<const D: usize> {
    /// Empty lattice cells forming the cavity, sorted — determinism (§13.1).
    pub cells: Vec<i32>,
    /// Flat index of the cell nearest the cavity's centre of mass.
    pub centre: i32,
    /// Mean occupied-neighbour count across the cavity's cells.
    pub enclosure: u8,
    /// Local signature, built from the monomers lining the cavity. This is
    /// what a substrate's signature is compared against (§8.3).
    pub signature: Signature<D>,
}

#[must_use]
pub fn cavities<const D: usize>(
    f: &Fold,
    p: &Polymer,
    g: &Geodesic<D>,
    u: &Universe,
) -> Vec<Cavity<D>> {
    if f.n == 0 {
        return Vec::new();
    }
    let occupied: BTreeMap<i32, usize> =
        f.coords.iter().enumerate().map(|(i, &c)| (c, i)).collect();

    // Candidate cells: empty, in bounds, and enclosed enough to matter.
    let mut candidates: BTreeSet<i32> = BTreeSet::new();
    for &c in &f.coords {
        for &d in &fcc::NEIGHBOURS {
            let cell = c + d;
            if occupied.contains_key(&cell) || !fcc::in_bounds(cell) {
                continue;
            }
            let filled =
                fcc::NEIGHBOURS.iter().filter(|&&e| occupied.contains_key(&(cell + e))).count() as u8;
            if filled >= MIN_ENCLOSURE {
                candidates.insert(cell);
            }
        }
    }

    // Flood-fill candidates into connected groups. BTreeSet iteration is
    // ordered, so group discovery order is fixed.
    let mut seen: BTreeSet<i32> = BTreeSet::new();
    let mut out = Vec::new();
    for &start in &candidates {
        if seen.contains(&start) {
            continue;
        }
        let mut cells = Vec::new();
        let mut q = VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            if !seen.insert(x) {
                continue;
            }
            cells.push(x);
            if cells.len() > MAX_CAVITY_CELLS {
                break;
            }
            for &d in &fcc::NEIGHBOURS {
                if candidates.contains(&(x + d)) && !seen.contains(&(x + d)) {
                    q.push_back(x + d);
                }
            }
        }
        if cells.is_empty() || cells.len() > MAX_CAVITY_CELLS {
            continue;
        }
        cells.sort_unstable();

        let enclosure = (cells
            .iter()
            .map(|&c| fcc::NEIGHBOURS.iter().filter(|&&e| occupied.contains_key(&(c + e))).count())
            .sum::<usize>()
            / cells.len()) as u8;

        let centre = cells[cells.len() / 2];
        let signature = lining_signature(&cells, &occupied, p, g, u);
        out.push(Cavity { cells, centre, enclosure, signature });
    }
    out
}

/// Signature of a cavity, built from the monomers lining it.
///
/// Directions point *inward* from the lining monomers toward the cavity
/// centre, so the resulting signature is what a substrate sitting in the
/// cavity would have to complement — the same comparison as §8.3, applied to
/// a hole instead of a molecule.
fn lining_signature<const D: usize>(
    cells: &[i32],
    occupied: &BTreeMap<i32, usize>,
    p: &Polymer,
    g: &Geodesic<D>,
    u: &Universe,
) -> Signature<D> {
    let mut sig = Signature::<D>::zeroed();
    let (cx, cy, cz) = fcc::from_flat(cells[cells.len() / 2]);
    let centre = [f64::from(cx), f64::from(cy), f64::from(cz)];

    // Lining monomers: those adjacent to any cavity cell. BTreeSet keeps the
    // set ordered so the summation order below is fixed.
    let mut lining: BTreeSet<usize> = BTreeSet::new();
    for &c in cells {
        for &d in &fcc::NEIGHBOURS {
            if let Some(&m) = occupied.get(&(c + d)) {
                lining.insert(m);
            }
        }
    }

    for dir_idx in 0..D {
        let dir = g.dirs[dir_idx];
        let mut extent = f64::NEG_INFINITY;
        let mut num = 0.0;
        let mut den = 0.0;
        for &m in &lining {
            let (mx, my, mz) = fcc::from_flat_of(m, occupied);
            let rel = [f64::from(mx) - centre[0], f64::from(my) - centre[1], f64::from(mz) - centre[2]];
            let reach = rel[0] * dir[0] + rel[1] * dir[1] + rel[2] * dir[2];
            if reach > extent {
                extent = reach;
            }
            let w = reach.max(0.0);
            num += w * u.element(p.unit(m)).affinity;
            den += w;
        }
        sig.r[dir_idx] = if extent.is_finite() { extent } else { 0.0 };
        sig.a[dir_idx] = if den > 0.0 { (num / den).clamp(-1.0, 1.0) } else { 0.0 };
    }
    sig
}

/// Look up a monomer's lattice coordinates via the occupancy map.
fn from_flat_of(monomer: usize, occupied: &BTreeMap<i32, usize>) -> (i32, i32, i32) {
    for (&cell, &m) in occupied {
        if m == monomer {
            return fcc::from_flat(cell);
        }
    }
    (0, 0, 0)
}
```

> **Implementer note:** `from_flat_of` as written is a linear scan per
> direction per monomer, which is fine for correctness but wasteful. Build a
> `Vec<(i32,i32,i32)>` indexed by monomer once at the top of
> `lining_signature` and index into it. Kept explicit here so the intent is
> visible rather than hidden in a setup line.

- [ ] **Step 4: Run** — `cargo test -p borbax-molecule cavity`, expect 4 passed

`longer_chains_produce_more_cavities` is the one to watch. If it fails,
folding is not compacting enough and Task 11's move set needs pull moves —
adjusting `MIN_ENCLOSURE` downward would "fix" the test while destroying the
selectivity that makes catalysis possible.

- [ ] **Step 5: Commit**

```bash
git add crates/borbax-molecule/src/cavity.rs
git commit -m "$(cat <<'EOF'
feat(molecule): cavity detection on folded chains

A cavity is a connected pocket of empty lattice cells enclosed on at least
eight of twelve faces. That threshold is the 3D payoff: in 2D the
equivalent is a notch open on two sides, and §8.5's selectivity argument
does not hold for it. Specificity is the difference between a catalyst and
a patch of glue.

Cavity signatures are built from the lining monomers with directions
pointing inward, so the result is what a substrate would have to
complement — the same comparison as §8.3, applied to a hole.

BTreeSet and BTreeMap throughout rather than hash containers: group
discovery order and summation order both affect results, so they are fixed
(§13.1).

MIN_ENCLOSURE is not a tuning knob. If longer chains stop producing more
cavities, the fix is pull moves in folding, not a lower threshold.

MrReasonable <4990954+MrReasonable@users.noreply.github.com>
EOF
)"
```

---
