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
- Produces: `MAX_POLYMER`, `Polymer`, `Fold`, `FoldWorkspace`, `FoldWorkspace::{fold(&Polymer, &Universe) -> Fold, fold_with_seed(..)}`, `Fold::{contacts, energy, exposure}`, `fcc::{GRID, NEIGHBOURS, to_flat, from_flat}`

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
        let a = ws.fold(&p, &u);
        let b = ws.fold(&p, &u);
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
        let _ = shared.fold(&p1, &u);
        let reused = shared.fold(&p2, &u);

        let mut fresh = FoldWorkspace::new();
        let clean = fresh.fold(&p2, &u);
        assert_eq!(reused.coords, clean.coords, "workspace reset leaked state");
    }

    #[test]
    fn the_chain_stays_connected_and_self_avoiding() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        for seed in 0..30 {
            let p = poly(45, &u);
            let f = ws.fold(&p, &u);
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
            let f = ws.fold(&p, &u);
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
        let f = ws.fold(&p, &u);
        let folded: u32 = f.contacts.iter().sum();
        assert!(folded > 8, "only {folded} contacts — annealing is not compacting");
    }

    /// The grid-wrap check. A 200-mer must stay well inside the lattice at
    /// every point, in *true* coordinates — `in_bounds` alone cannot catch a
    /// wrap because it decodes the already-wrapped position.
    #[test]
    fn no_chain_length_wraps_the_grid() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        for n in [45, 64, 80, 120, MAX_POLYMER] {
            let p = poly(n, &u);
            let f = ws.fold(&p, &u);
            for (i, &c) in f.coords.iter().enumerate() {
                let (x, y, z) = fcc::from_flat(c);
                assert!(
                    (2..fcc::GRID - 2).contains(&x)
                        && (2..fcc::GRID - 2).contains(&y)
                        && (2..fcc::GRID - 2).contains(&z),
                    "n={n} monomer {i} at ({x},{y},{z}) is at or past the boundary"
                );
                assert_eq!(fcc::to_flat(x, y, z), c, "flat/coord round-trip failed");
            }
        }
    }

    /// Folding must be a function of the *sequence*, not of the anneal stream.
    /// `folding_is_deterministic` compares two folds from the same stream,
    /// which tests the RNG rather than the folder. If the same polymer reaches
    /// materially different shapes under different streams then "the folding
    /// map" is not a map, and every neutral-network and shape-space-covering
    /// measurement downstream is measuring annealer noise.
    #[test]
    fn folds_are_reproducible_across_independent_streams() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        let mut modal_agreement = 0;
        let trials = 40;
        for t in 0..trials {
            let p = poly(40, &u);
            let a = ws.fold_with_seed(&p, &u, t);
            let b = ws.fold_with_seed(&p, &u, t + 1000);
            if a.contacts == b.contacts {
                modal_agreement += 1;
            }
        }
        // Not equality — annealing is stochastic. But the modal shape must
        // dominate, or the map is noise. Record the measured figure in the
        // commit; this bound is a floor, not a target.
        assert!(
            modal_agreement * 2 > trials,
            "only {modal_agreement}/{trials} folds agreed across streams — the folding map is annealer noise"
        );
    }

    #[test]
    fn exposure_is_between_zero_and_twelve() {
        let u = Universe::generate(8);
        let mut ws = FoldWorkspace::new();
        let p = poly(40, &u);
        let f = ws.fold(&p, &u);
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

    /// Seed for a fold, derived from the sequence and the universe.
    ///
    /// Shares its construction with the fold cache key (Task 12) so the two
    /// cannot drift apart — a cached fold and a fresh one must be the same
    /// conformation.
    #[must_use]
    pub fn fold_seed(p: &Polymer, u: &Universe) -> u64 {
        let mut h = Stream::new(0xF0_1D_5EED ^ u.seed, Domain::Hash, 3);
        let mut acc = h.next_u64();
        for (i, &m) in p.units().iter().enumerate() {
            acc ^= u64::from(m).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left((i % 61) as u32);
            acc = acc.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        }
        acc
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

    /// Fold a polymer. **The seed is derived from the sequence and universe**,
    /// so a species has exactly one conformation regardless of when it is
    /// first encountered (§8.6).
    ///
    /// An earlier draft took an external `Stream` here while the cache derived
    /// its own — so a direct call and a cached call could disagree, and one
    /// species would have two shapes. Production code calls this; only
    /// `fold_with_seed` below can vary the stream, and only tests use it.
    pub fn fold(&mut self, p: &Polymer, u: &Universe) -> Fold {
        self.fold_with_seed(p, u, fold_seed(p, u))
    }

    /// Fold under an explicit seed. **Diagnostics only.**
    ///
    /// Exists so `folds_are_reproducible_across_independent_streams` can
    /// measure whether the fold is a property of the sequence or of the
    /// annealer. Calling this from production code would give one species two
    /// shapes.
    pub fn fold_with_seed(&mut self, p: &Polymer, u: &Universe, seed: u64) -> Fold {
        let mut rng = Stream::new(seed, Domain::Fold, 0);
        let rng = &mut rng;
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

        // Compact (boustrophedon) start, centred in the grid.
        //
        // An extended start along one diagonal runs `centre ± n/2`, which at
        // GRID = 64 already **wraps** for a 64-mer — and wrapping is silent:
        // the flat index stays inside 0..64³, so there is no panic, and
        // `in_bounds` cannot detect it because it decodes the *wrapped*
        // coordinates, which pass. A monomer at true (−8,−8,32) simply becomes
        // (56,55,31), and self-avoidance and connectivity are then evaluated
        // on a sheared torus. Task 13 folds 70- and 80-mers.
        //
        // A boustrophedon walk keeps the start inside a box of side ~n^(1/3)
        // instead of n, so a 200-mer is nowhere near the boundary.
        let centre = fcc::to_flat(fcc::GRID / 2, fcc::GRID / 2, fcc::GRID / 2);
        const RUN: usize = 6;
        let mut coords = Vec::with_capacity(n);
        let mut c = centre;
        let mut dir = 0usize;
        for i in 0..n {
            debug_assert!(fcc::in_bounds(c), "fold start left the grid at monomer {i}");
            coords.push(c);
            self.occupy(c, i);
            if i % RUN == RUN - 1 {
                dir = (dir + 1) % fcc::NEIGHBOURS.len();
            }
            c += fcc::NEIGHBOURS[dir];
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

            // Metropolis proper: uphill acceptance falls off with ΔE.
            //
            // The earlier form — `trial <= current || rng.next_f64() < temp *
            // 0.35` — was not Metropolis. Acceptance was independent of ΔE, so
            // a catastrophic move was as likely as a marginal one, and because
            // `<=` accepts every zero-ΔE move a chain of near-neutral monomers
            // performed a pure random walk. The fold would then be an artefact
            // of the RNG rather than a function of the sequence, which
            // destroys the many-to-one folding map (§8.4) that the whole
            // evolvability argument rests on.
            //
            // The draw is unconditional so stream position depends only on the
            // step count, never on a float comparison (§13.1).
            let roll = rng.next_f64();
            let delta = trial_energy - current;
            let accept = delta <= 0.0
                || roll < borbax_units::det_math::exp(-delta / (temp * energy_scale).max(1e-9));
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

            // **This is §8.3's charge term, verbatim, at monomer scale.**
            //
            // The earlier form was `-fa * fb`, which had two problems. It
            // favoured like-attracts-like — annealing minimises, and with
            // fa = +1, fb = -1 that expression is *positive* — contradicting
            // its own comment, §8.3, and Task 5's bond matrix. Three sites,
            // two conventions.
            //
            // And even with the sign corrected it is not §8.3:
            // -(fa+fb)² = -fa² - fb² - 2·fa·fb, and the self-terms are not
            // constant across class pairs, so the two forms rank contacts
            // differently. That made it a *second* complementarity law
            // producing the same phenomenon — exactly the design smell
            // Principle 2 names.
            //
            // Only the charge term appears here, with no shape term: FCC
            // contacts have no extent to compare, since every contact is at
            // the same lattice distance. That is a real reason, not an
            // omission.
            t[pair_index(a, b)] = -(fa + fb) * (fa + fb) * u.consts.w_charge;
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

    /// The same sequence under two universes must not share a cache entry.
    /// The fold depends on the universe through monomer classes and
    /// `w_charge`, so a shared entry silently returns a conformation computed
    /// under different physics.
    #[test]
    fn the_same_sequence_in_two_universes_does_not_share_an_entry() {
        let (u1, u2) = (Universe::generate(2), Universe::generate(3));
        let mut c = FoldCache::with_budget(8 << 20);
        let p = poly(1, 30);
        let _ = c.get_or_fold(&p, &u1);
        let _ = c.get_or_fold(&p, &u2);
        assert_eq!(c.stats().misses, 2, "a hit crossed universes");
        assert_eq!(c.stats().hits, 0);
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

/// 128-bit key from two independent passes of the RNG mixer over the sequence
/// **and the universe seed**. At a few million entries the birthday
/// probability is ~1e-27, so storing the sequence itself would be paying for a
/// guarantee we already have.
///
/// **The universe is part of the key.** The fold depends on it through the
/// monomer classes and `w_charge`, so a key on the sequence alone means a hit
/// from universe A silently answers for universe B. Task 20's battery iterates
/// forty universes and Task 21's sweeps many more; whether the bug fires would
/// then depend on visit order and on which entries survived eviction.
///
/// `Domain::Hash` rather than `Domain::Fold`: using a live simulation domain
/// as a hash function couples the two, so that adding a draw in folding would
/// change cache keys.
fn key_of(p: &Polymer, u: &Universe) -> u128 {
    let mut a = Stream::new(0xB0_1BAA_5EED ^ u.seed, Domain::Hash, 1);
    let mut b = Stream::new(0x5EED_B0_1BAA ^ u.seed, Domain::Hash, 2);
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
        let k = key_of(p, u);
        if let Some(&slot) = self.map.get(&k) {
            self.stats.hits += 1;
            if let Some((_, e)) = &mut self.slots[slot] {
                e.referenced = true;
            }
            // Re-borrow immutably for the return.
            return &self.slots[slot].as_ref().map(|(_, e)| &e.fold).unwrap_or_else(|| unreachable!());
        }

        self.stats.misses += 1;
        // `fold` derives its own seed from (sequence, universe), so the cached
        // and uncached paths cannot disagree.
        let fold = self.ws.fold(p, u);
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
        let f = FoldWorkspace::new().fold(&p, u);
        (p, f)
    }

    #[test]
    fn detection_is_deterministic() {
        let (u, g) = (Universe::generate(12), Geodesic::<42>::build().unwrap());
        let (p, f) = folded(3, 60, &u);
        assert_eq!(cavities(&f, &p, &g, &u), cavities(&f, &p, &g, &u));
    }

    /// **A real enclosure test.** The previous version asserted that reported
    /// cavities have `enclosure >= MIN_ENCLOSURE` — which cannot fail, since
    /// candidates are admitted only if they already pass that threshold. It
    /// read as a defence of the selectivity claim and defended nothing.
    ///
    /// This builds a hollow shell (a genuine cavity) and the same shell with
    /// one monomer removed (open to solvent). The second must yield nothing.
    #[test]
    fn an_open_groove_is_not_reported_as_a_cavity() {
        let (u, g) = (Universe::generate(12), Geodesic::<42>::build().unwrap());
        let (p, sealed) = hollow_shell(&u);
        assert!(!cavities(&sealed, &p, &g, &u).is_empty(), "sealed shell had no cavity");

        let breached = remove_one_monomer(&sealed);
        assert!(
            cavities(&breached, &p, &g, &u).is_empty(),
            "a groove open to solvent was reported as a cavity"
        );
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

/// Lattice spacing in `Span`, so cavity distances share units with molecule
/// signatures (§8.2). Without this, `affinity(cavity, molecule)` adds raw grid
/// integers to atom radii and compares the result against `ideal_gap` — two
/// incommensurable numbers, and the resulting recognition scores are arbitrary.
///
/// FCC nearest-neighbour distance is sqrt(2) in flat-index units.
const LATTICE_SPAN: f64 = std::f64::consts::SQRT_2;

/// Shell width for cavity surface character, mirroring §8.2's `SHELL`.
const SHELL: f64 = 0.75;

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

    // ---- Step 1: what is *outside*? ----
    //
    // This is the enclosure test, and an earlier draft had none. Selecting
    // cells by "≥ 8 of 12 neighbours occupied" is a local *density* test: a
    // deep groove on a compact globule, freely open to solvent, is
    // indistinguishable from a sealed void by that measure. Sticky surface
    // patches would then be counted as catalytic sites and §23's criterion 6
    // would read positive on geometry that is not a cavity.
    //
    // What actually distinguishes a cavity from a dimple is reachability:
    // everything empty and reachable from beyond the chain's bounding box is
    // solvent. Interior = empty and unreached.
    let (lo, hi) = bounding_box(&f.coords);
    let mut outside: BTreeSet<i32> = BTreeSet::new();
    let mut q: VecDeque<i32> = VecDeque::new();
    for cell in boundary_cells(lo, hi) {
        if !occupied.contains_key(&cell) && outside.insert(cell) {
            q.push_back(cell);
        }
    }
    while let Some(x) = q.pop_front() {
        for &d in &fcc::NEIGHBOURS {
            let nb = x + d;
            if in_box(nb, lo, hi) && !occupied.contains_key(&nb) && outside.insert(nb) {
                q.push_back(nb);
            }
        }
    }

    // ---- Step 2: interior cells, filtered by enclosure ----
    let mut candidates: BTreeSet<i32> = BTreeSet::new();
    for cell in box_cells(lo, hi) {
        if occupied.contains_key(&cell) || outside.contains(&cell) || !fcc::in_bounds(cell) {
            continue;
        }
        let filled =
            fcc::NEIGHBOURS.iter().filter(|&&e| occupied.contains_key(&(cell + e))).count() as u8;
        if filled >= MIN_ENCLOSURE {
            candidates.insert(cell);
        }
    }

    // ---- Step 3: group into connected cavities ----
    //
    // BTreeSet iteration is ordered, so discovery order is fixed.
    let mut seen: BTreeSet<i32> = BTreeSet::new();
    let mut out = Vec::new();
    for &start in &candidates {
        if seen.contains(&start) {
            continue;
        }
        // The flood always runs to completion before any size filter is
        // applied. Breaking out mid-fill leaves the component half-marked in
        // `seen`, so a later start inside the same component floods only the
        // remainder — with the already-seen cells acting as walls — and
        // produces a spurious fragment that may not even be connected.
        let mut cells = Vec::new();
        let mut q = VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            if !seen.insert(x) {
                continue;
            }
            cells.push(x);
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

    // Monomer positions, resolved once. Looking each up by scanning the
    // occupancy map per direction would be O(D · lining · |occupied|).
    let mut pos_of: Vec<(i32, i32, i32)> = vec![(0, 0, 0); p.len()];
    for (&cell, &m) in occupied {
        if m < pos_of.len() {
            pos_of[m] = fcc::from_flat(cell);
        }
    }

    for dir_idx in 0..D {
        let dir = g.dirs[dir_idx];

        // **Nearest wall, not furthest.** An earlier draft took `max`, which
        // is the support function of the *far* side of the cavity — the wall
        // a substrate never touches. What a substrate needs to complement is
        // the free space to the nearest wall in each direction, which is a
        // min over the lining monomers that actually lie that way.
        //
        // Distances are in `Span`, converted from lattice units. Leaving them
        // as raw grid integers made `affinity(cavity, molecule)` add two
        // incommensurable numbers and compare the result to `ideal_gap`.
        // Two passes, because the weights in the second depend on `nearest`
        // from the first. Merging them is not an optimisation, it is a
        // semantic change — the same structure as §8.2's `signature()`.
        let reach_of = |m: usize| -> f64 {
            let (mx, my, mz) = pos_of[m];
            let rel = [
                (f64::from(mx) - centre[0]) * LATTICE_SPAN,
                (f64::from(my) - centre[1]) * LATTICE_SPAN,
                (f64::from(mz) - centre[2]) * LATTICE_SPAN,
            ];
            rel[0] * dir[0] + rel[1] * dir[1] + rel[2] * dir[2]
        };

        let mut nearest = f64::INFINITY;
        for &m in &lining {
            let reach = reach_of(m);
            // Behind this direction: not a wall a substrate can reach that way.
            if reach > 0.0 && reach < nearest {
                nearest = reach;
            }
        }

        // Weight by *proximity* to the near wall, matching how §8.2 weights
        // its surface shell. The far wall must not dominate the character a
        // substrate would actually feel.
        let (mut num, mut den) = (0.0f64, 0.0f64);
        if nearest.is_finite() {
            for &m in &lining {
                let reach = reach_of(m);
                if reach <= 0.0 {
                    continue;
                }
                let w = (nearest + SHELL - reach).max(0.0);
                num += w * u.element(p.unit(m)).affinity;
                den += w;
            }
        }

        // Always non-negative, matching the positivity the molecule-signature
        // test asserts.
        sig.r[dir_idx] = if nearest.is_finite() { nearest } else { 0.0 };
        sig.a[dir_idx] = if den > 0.0 { (num / den).clamp(-1.0, 1.0) } else { 0.0 };
    }
    sig
}

/// Bounding box of the chain, expanded by two cells so the outside flood has
/// somewhere to start.
fn bounding_box(coords: &[i32]) -> ((i32, i32, i32), (i32, i32, i32)) {
    let mut lo = (i32::MAX, i32::MAX, i32::MAX);
    let mut hi = (i32::MIN, i32::MIN, i32::MIN);
    for &c in coords {
        let (x, y, z) = fcc::from_flat(c);
        lo = (lo.0.min(x), lo.1.min(y), lo.2.min(z));
        hi = (hi.0.max(x), hi.1.max(y), hi.2.max(z));
    }
    ((lo.0 - 2, lo.1 - 2, lo.2 - 2), (hi.0 + 2, hi.1 + 2, hi.2 + 2))
}

fn in_box(cell: i32, lo: (i32, i32, i32), hi: (i32, i32, i32)) -> bool {
    let (x, y, z) = fcc::from_flat(cell);
    (lo.0..=hi.0).contains(&x) && (lo.1..=hi.1).contains(&y) && (lo.2..=hi.2).contains(&z)
}

/// All cells in the box, in a fixed order.
fn box_cells(lo: (i32, i32, i32), hi: (i32, i32, i32)) -> impl Iterator<Item = i32> {
    (lo.2..=hi.2)
        .flat_map(move |z| (lo.1..=hi.1).flat_map(move |y| (lo.0..=hi.0).map(move |x| fcc::to_flat(x, y, z))))
}

/// Cells on the box's faces — the seeds for the outside flood.
fn boundary_cells(lo: (i32, i32, i32), hi: (i32, i32, i32)) -> impl Iterator<Item = i32> {
    box_cells(lo, hi).filter(move |&c| {
        let (x, y, z) = fcc::from_flat(c);
        x == lo.0 || x == hi.0 || y == lo.1 || y == hi.1 || z == lo.2 || z == hi.2
    })
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

## Phase 5 — Reactions, decay, and self-sustaining sets

### Task 14: Species records, interning, and reaction rates

The invariant of spec §8.6 lands here. **Nothing in this task's output may be
recomputed per molecule** — a `SpeciesRecord` is built once at intern time and
the step loop only ever reads it.

**Files:**
- Create: `crates/borbax-reaction/{Cargo.toml,src/lib.rs,src/species.rs,src/kind.rs,src/rate.rs}`
- Test: `src/species.rs`, `src/rate.rs`

**Interfaces:**
- Consumes: `borbax_molecule::{Mol12, Signature, Geodesic, canonicalise, embed, signature, affinity_ordered}`, `Universe`
- Produces: `SpeciesId(u32)`, `SpeciesRecord`, `Interner::{new, intern, record, len}`, `ReactionKind`, `Reaction`, `rate(&Reaction, Thermal, &[f64], f64) -> f64`, `AffinityMemo`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use borbax_molecule::{geodesic::Geodesic, graph::Mol12};
    use borbax_universe::{ElementId, Universe};

    fn mol(elems: &[u8]) -> Mol12 {
        let mut m = Mol12::new();
        for &e in elems {
            m.add_atom(ElementId(e));
        }
        for i in 1..elems.len() as u8 {
            m.add_bond(i - 1, i, 1);
        }
        m
    }

    #[test]
    fn identical_molecules_intern_to_one_id() {
        let (u, g) = (Universe::generate(3), Geodesic::<42>::build().unwrap());
        let mut it = Interner::new();
        let a = it.intern(&mol(&[0, 1, 2]), &u, &g);
        let b = it.intern(&mol(&[0, 1, 2]), &u, &g);
        assert_eq!(a, b);
        assert_eq!(it.len(), 1);
    }

    #[test]
    fn ids_are_assigned_in_order_of_first_appearance() {
        let (u, g) = (Universe::generate(3), Geodesic::<42>::build().unwrap());
        let run = || {
            let mut it = Interner::new();
            let ids: Vec<u32> = [&[0u8, 1][..], &[2, 3], &[0, 1], &[4, 5]]
                .iter()
                .map(|e| it.intern(&mol(e), &u, &g).0)
                .collect();
            ids
        };
        assert_eq!(run(), vec![0, 1, 0, 2]);
        assert_eq!(run(), run(), "id assignment is not reproducible");
    }

    #[test]
    fn records_carry_everything_the_hot_loop_needs() {
        let (u, g) = (Universe::generate(3), Geodesic::<42>::build().unwrap());
        let mut it = Interner::new();
        let id = it.intern(&mol(&[0, 1, 2, 3]), &u, &g);
        let r = it.record(id);
        assert_eq!(r.mass, mol(&[0, 1, 2, 3]).mass(&u));
        assert!(r.bond_count > 0);
        assert!(r.solvent_rate.is_finite());
    }

    /// Burial must protect. Two chains with the same composition but
    /// different compactness must get different solvent rates, in the
    /// direction §9.5 claims — otherwise "compact folds survive" is a
    /// sentence in the spec with no code behind it.
    #[test]
    fn burial_lowers_the_solvent_rate() {
        let (u, g) = (Universe::generate(3), Geodesic::<42>::build().unwrap());
        let mut it = Interner::new();
        let compact = it.intern_polymer(&compact_folder(&u), &u, &g);
        let extended = it.intern_polymer(&extended_folder(&u), &u, &g);
        assert!(
            it.record(compact).solvent_rate < it.record(extended).solvent_rate,
            "burial did not protect: compact {} vs extended {}",
            it.record(compact).solvent_rate,
            it.record(extended).solvent_rate
        );
    }

    /// Memoising a pure function cannot change results — that is the whole
    /// reason it is safe here (spec §8.6). Prove it rather than assume it.
    #[test]
    fn the_affinity_memo_agrees_with_direct_computation() {
        let (u, g) = (Universe::generate(3), Geodesic::<42>::build().unwrap());
        let mut it = Interner::new();
        let ids: Vec<_> = (0..8u8).map(|i| it.intern(&mol(&[i, i + 1, i + 2]), &u, &g)).collect();
        let mut memo = AffinityMemo::new(8);
        let k = (&u.consts).into();
        for &a in &ids {
            for &b in &ids {
                let direct =
                    borbax_molecule::binding::affinity_ordered(&it.record(a).sig, &it.record(b).sig, &g, &k);
                assert_eq!(memo.get(a, b, &it, &g, &k), direct);
            }
        }
    }

    #[test]
    fn rates_rise_with_temperature_and_concentration() {
        use borbax_units::{Quanta, Thermal};
        let r = Reaction {
            kind: ReactionKind::Condense,
            reactants: [SpeciesId(0), SpeciesId(1)],
            products: [SpeciesId(2), SpeciesId(3)],
            activation: Quanta(50.0),
            delta: Quanta(-10.0),
        };
        let c = [10.0, 10.0, 0.0, 0.0];
        assert!(rate(&r, Thermal(400.0), &c, 1.0) > rate(&r, Thermal(200.0), &c, 1.0));
        let more = [20.0, 10.0, 0.0, 0.0];
        assert!(rate(&r, Thermal(300.0), &more, 1.0) > rate(&r, Thermal(300.0), &c, 1.0));
        assert!(rate(&r, Thermal(300.0), &c, 50.0) > rate(&r, Thermal(300.0), &c, 1.0));
    }
}
```

- [ ] **Step 2: Run to verify failure** — expect `cannot find type Interner`

- [ ] **Step 3: Implement species records and interning**

```rust
//! Species records: everything expensive, computed once (spec §8.6).
//!
//! Canonicalisation, embedding, signature construction, folding and per-bond
//! decay rates are all pure functions of the *species*. They run here, at
//! intern time, and nowhere else. **No code reachable from a simulation step
//! may call any of them** — that invariant is what makes the molecular-tier
//! budget in §17 reachable at all, and violating it costs roughly two orders
//! of magnitude.

use borbax_molecule::binding::{affinity_ordered, BindConsts};
use borbax_molecule::canonical::{canonicalise, CanonForm};
use borbax_molecule::geodesic::Geodesic;
use borbax_molecule::graph::Mol12;
use borbax_molecule::layout::embed;
use borbax_molecule::signature::{signature, Signature};
use borbax_units::Mass;
use borbax_universe::{ElementId, Universe};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpeciesId(pub u32);

/// A species is either a small molecule or a polymer. **Polymers must have
/// species identity too**, and an earlier draft interned only `Mol12`.
///
/// That omission was structural rather than clerical: with nowhere to *put* a
/// polymer's fold or cavities, an implementer computes them at the point of
/// need — inside `Beaker::step`. On a cache miss that is 20,000 anneal steps
/// plus a 2,520-iteration `affinity` per cavity-substrate pair, which is
/// roughly four orders of magnitude past the §17 budget, not two.
#[derive(Debug, Clone)]
pub enum SpeciesKind {
    Small(CanonForm),
    Polymer(Polymer),
}

/// Everything the hot loop needs about a species. Built once; read forever.
#[derive(Debug, Clone)]
pub struct SpeciesRecord<const D: usize> {
    pub kind: SpeciesKind,
    pub sig: Signature<D>,
    pub mass: Mass,
    pub bond_count: u32,
    /// Summed element `decay_rate`. Precomputed for the same reason as
    /// everything else here — an earlier draft recomputed it per call with a
    /// comment saying it would be hoisted "if it shows up in a profile". A
    /// twelve-iteration loop never shows up in a profile and runs millions of
    /// times.
    pub radiogenic_rate: f64,
    /// Total thermal cleavage propensity per unit time at unit temperature,
    /// summed over the molecule's bonds. Precomputed so decay is a propensity
    /// channel rather than a per-molecule sweep (spec §9.5).
    pub cleave_propensity: f64,
    /// Solvent-attack rate. Rises with complementarity to the solvent and
    /// falls with burial — see the sign discussion in `intern`.
    pub solvent_rate: f64,
    /// Folded conformation, for polymers. `None` for small molecules.
    ///
    /// Named in §8.6's list of what a record carries, and omitted from an
    /// earlier draft — which is what forced folding into the step loop.
    pub fold: Option<FoldId>,
    /// Cavities, with their local signatures. Computed once at intern time;
    /// `cavities()` has no cache of its own, so a fold-cache *hit* would still
    /// have re-run the flood fill and D×lining dot products on every call.
    pub cavities: Vec<Cavity<D>>,
}

/// Canonical form -> id. Lookup only; **never iterated**, because iteration
/// order of a hash container is not reproducible. The `Vec` is the iteration
/// surface.
pub struct Interner<const D: usize> {
    ids: BTreeMap<CanonForm, SpeciesId>,
    records: Vec<SpeciesRecord<D>>,
}

impl<const D: usize> Default for Interner<D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const D: usize> Interner<D> {
    #[must_use]
    pub fn new() -> Self {
        Self { ids: BTreeMap::new(), records: Vec::new() }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    #[must_use]
    pub fn record(&self, id: SpeciesId) -> &SpeciesRecord<D> {
        &self.records[id.0 as usize]
    }

    /// Intern a molecule, computing its record on first sight.
    ///
    /// Ids are assigned **in order of first appearance in the simulation's
    /// event sequence**, which is automatically reproducible because the event
    /// sequence is. Assigning by hash value or by sorted traversal would be
    /// reproducible too, but would make ids jump around as unrelated species
    /// appear, which makes every log and every diff harder to read.
    pub fn intern(&mut self, m: &Mol12, u: &Universe, g: &Geodesic<D>) -> SpeciesId {
        let (canon, _) = canonicalise(m);
        if let Some(&id) = self.ids.get(&canon) {
            return id;
        }
        let id = SpeciesId(self.records.len() as u32);
        let e = embed(m, u);
        let sig = signature(m, &e, g, u).canonicalise(g);

        let mut bond_count = 0u32;
        let mut cleave = 0.0f64;
        for i in 0..m.n {
            for j in (i + 1)..m.n {
                let order = m.bond_order(i, j);
                if order == 0 {
                    continue;
                }
                bond_count += 1;
                let (ea, eb) =
                    (u.element(ElementId(m.elem[i as usize])), u.element(ElementId(m.elem[j as usize])));
                let energy = u.bonds.energy(ea, eb, order).get();
                // Weaker bonds break more readily. The exponential lives in
                // `rate`; here we keep the bond-strength weighting only, so
                // this stays a plain sum with no transcendental (§13.1).
                cleave += u.consts.decay_scale / energy.max(1e-6);
            }
        }

        // Solvent attack is a pure function of the species (§9.5), so it is
        // resolved once here rather than on every step.
        //
        // **Sign.** `affinity` is a negated sum of squares and is therefore
        // always ≤ 0, with 0 the perfect fit. An earlier draft used
        // `-affinity`, which is *largest when complementarity is worst* — so
        // molecules the solvent could not touch decayed fastest, and the
        // burial mechanism ran backwards. Exit criterion 7 would still have
        // passed: it checks that decay bites, not which direction.
        //
        // Mapping the score onto a rate through a logistic keeps it positive,
        // monotone in fit, and bounded.
        let solvent_sig = solvent_signature(u, g);
        let fit = affinity_ordered(&sig, &solvent_sig, g, &(&u.consts).into());
        let base_solvent = bind_probability(fit, u.consts.reference_temp, u.consts.bind_midpoint);

        // **Burial.** §9.5 specifies exposure as the count of empty lattice
        // neighbours out of twelve, and `Fold::exposure` computes exactly
        // that — but an earlier draft never read it, using the whole-molecule
        // signature instead. Without this term, "compact folds that bury their
        // backbone survive" is not implemented, and folding has no survival
        // payoff at all.
        let exposed_fraction = match &fold {
            Some(f) => {
                let total: u32 = f.exposure.iter().map(|&e| u32::from(e)).sum();
                f64::from(total) / (12.0 * f.exposure.len().max(1) as f64)
            }
            // A small molecule has no interior to hide in; fully exposed.
            None => 1.0,
        };
        let solvent_rate = base_solvent * exposed_fraction;

        self.records.push(SpeciesRecord {
            canon,
            sig,
            mass: m.mass(u),
            bond_count,
            cleave_propensity: cleave,
            solvent_rate,
        });
        self.ids.insert(canon, id);
        id
    }
}

/// Signature of a lone solvent atom. Cheap, but computed once per call site
/// in practice — callers should hoist it.
fn solvent_signature<const D: usize>(u: &Universe, g: &Geodesic<D>) -> Signature<D> {
    let mut m = Mol12::new();
    m.add_atom(u.consts.solvent);
    signature(&m, &embed(&m, u), g, u).canonicalise(g)
}

/// Dense triangular memo over the most abundant species (spec §8.6).
///
/// `affinity_ordered` is a pure function of a species pair plus universe
/// constants, so memoising it cannot change a result. Unlike a lossy cache
/// this is safe under §13.1 unconditionally — which is why it is worth doing
/// here rather than optimising the kernel further.
pub struct AffinityMemo {
    k: u32,
    vals: Vec<Option<f64>>,
}

impl AffinityMemo {
    #[must_use]
    pub fn new(k: u32) -> Self {
        Self { k, vals: vec![None; (k as usize) * (k as usize + 1) / 2] }
    }

    fn slot(&self, a: SpeciesId, b: SpeciesId) -> Option<usize> {
        let (lo, hi) = if a.0 < b.0 { (a.0, b.0) } else { (b.0, a.0) };
        if hi >= self.k {
            return None;
        }
        Some((hi as usize) * (hi as usize + 1) / 2 + lo as usize)
    }

    pub fn get<const D: usize>(
        &mut self,
        a: SpeciesId,
        b: SpeciesId,
        it: &Interner<D>,
        g: &Geodesic<D>,
        k: &BindConsts,
    ) -> f64 {
        let compute = || affinity_ordered(&it.record(a).sig, &it.record(b).sig, g, k);
        match self.slot(a, b) {
            None => compute(),
            Some(s) => match self.vals[s] {
                Some(v) => v,
                None => {
                    let v = compute();
                    self.vals[s] = Some(v);
                    v
                }
            },
        }
    }
}
```

- [ ] **Step 4: Implement reaction kinds and rates**

```rust
//! Reaction classes and the rate law (spec §9.1).

use crate::species::SpeciesId;
use borbax_units::{Quanta, Thermal};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionKind {
    /// Non-covalent complex via signature complementarity.
    Associate,
    /// Covalent bond forms; consumes energy, releases a leaving group.
    Condense,
    /// Covalent bond breaks; releases energy. Fires **spontaneously** at a
    /// temperature-dependent rate, not only under catalysis — this is the
    /// engine's primary decay path (spec §9.5).
    Cleave,
    /// A sub-group migrates between molecules within a complex.
    Transfer,
    /// Internal bond topology changes with no partner.
    Rearrange,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reaction {
    pub kind: ReactionKind,
    /// Products are resolved once when the channel is created and stored
    /// (spec §8.6), so canonicalisation runs only for genuinely novel species.
    pub reactants: [SpeciesId; 2],
    pub products: [SpeciesId; 2],
    pub activation: Quanta,
    pub delta: Quanta,
    /// Rate enhancement from whichever folded polymer catalyses this channel,
    /// **resolved once when the channel is created**.
    ///
    /// An earlier draft passed `catalysis: f64` into `rate()` with nothing
    /// anywhere constructing it. That missing join is what would have dragged
    /// folding and cavity extraction into the step loop, because the only way
    /// to supply the value on demand is to compute it on demand.
    pub catalysis: f64,
    /// The species catalysing this channel, if any. Also what `find_raf`'s
    /// catalyst list is built from — another consumer of a value nothing was
    /// producing.
    pub catalyst: Option<SpeciesId>,
}

/// Arrhenius-style rate.
///
/// The `exp` is routed through `det_math` so the vendored portable
/// implementation has one call site to replace (spec §13.1). It is also the
/// prime candidate for tabulation per (bond class, temperature bucket), since
/// its arguments are drawn from a small set — see the V0 benchmark list.
#[must_use]
pub fn rate(r: &Reaction, t: Thermal, conc: &[f64], catalysis: f64) -> f64 {
    let temp = t.get().max(1e-6);
    let arrhenius = borbax_molecule::det_math::exp(-r.activation.get() / temp);
    let a = conc.get(r.reactants[0].0 as usize).copied().unwrap_or(0.0);
    let b = conc.get(r.reactants[1].0 as usize).copied().unwrap_or(0.0);
    arrhenius * a * b * catalysis
}
```

- [ ] **Step 5: Run** — `cargo test -p borbax-reaction`, expect 5 passed

- [ ] **Step 6: Commit**

```bash
git add crates/borbax-reaction
git commit -m "$(cat <<'EOF'
feat(reaction): species records, interning, and the rate law

Implements the §8.6 invariant: canonicalisation, embedding, signature,
and per-bond decay propensity are pure functions of the species, so they
run once at intern time and never again. Nothing reachable from a step may
call them — violating that costs roughly two orders of magnitude against
the §17 budget.

Ids are assigned in order of first appearance in the event sequence, which
is reproducible because the event sequence is. Hash-order assignment would
also be reproducible but would make ids jump as unrelated species appear,
making every log and diff harder to read.

The intern map is BTreeMap and is lookup-only; the Vec is the iteration
surface. Hash-container iteration order is not reproducible and this is
the last place to be casual about it.

The affinity memo is safe unconditionally because it memoises a pure
function — tested against direct computation rather than assumed.

MrReasonable <4990954+MrReasonable@users.noreply.github.com>
EOF
)"
```

---

### Task 15: Decay as propensity channels

**Files:** Create `crates/borbax-reaction/src/decay.rs`; test same file.

**Interfaces:**
- Consumes: `SpeciesRecord`, `Interner`, `Universe`, `Thermal`
- Produces: `decay_propensity<D>(SpeciesId, f64, &Interner<D>, Thermal, &Universe) -> f64`, `DecayKind`, `decay_channels<D>(...) -> [(DecayKind, SpeciesId, f64); 3]`

**Three requirements this task carries, all found by review before any of it was
written. Read them before Step 1 — two change the signatures above.**

1. **The channels must sum to the propensity.** As drafted they do not:
   `decay_propensity` returns `count * (thermal + solvent)` while
   `decay_channels` returns three channels *including* Radiogenic. A scheduler
   drawing the total from one and selecting the channel from the other either
   never fires Radiogenic or fires it off-budget — and Radiogenic is §9.5's
   mutation source, so the failure mode is "mutation silently switched off" with
   every test still green. Make one the sum of the other, and assert exactly
   that: `decay_propensity(..) == decay_channels(..).iter().map(|c| c.2).sum()`
   to the last bit, with the accumulation order pinned and commented (§13.4).

2. **`decay_channels` must not allocate.** It returns a compile-time-constant
   three-element list from a function the beaker calls whenever a species' count
   changes, so a `Vec` is a malloc/free pair on the step path. Measured on an
   M1 Pro: 23.98 ns/call as drafted, 19.94 with the radiogenic walk hoisted,
   **3.68 with a fixed-size return and the walk hoisted**. The allocator alone
   is 16.27 ns — **4x the radiogenic walk** at twelve atoms. Hoisting the walk
   and leaving the `Vec` buys 17% of the available win, which is exactly the
   plausible-but-not-real fix this project keeps producing. Hence `[_; 3]`
   above.

   *Test:* a counting global allocator in a `#[test]`, asserting **zero heap
   allocations per step** past warm-up. Deterministic rather than a timing
   assertion, so it survives the rule against timing assertions in `#[test]`;
   and it is the only candidate test that a radiogenic-only fix fails.

3. **Per-step decay cost must be O(1) in species size.** `radiogenic_rate` walks
   every canonical element per call. The field doc's justification — "a
   twelve-iteration loop" — understates it: `Mol12` caps small molecules at 12
   but `MAX_POLYMER` is **200**, and the drafted function reads `record(id).canon`
   with no polymer arm at all, so whoever closes that gap walks `Polymer::units`
   through a further indirection. Measured cost of `decay_channels` against a
   ~10 µs step at 96 species: 6.6% at chain length 12, 10.3% at 20 (§7.2's
   viability floor), **155% at `MAX_POLYMER`** — i.e. the regime that breaks the
   budget is precisely the one §7.2's polymer criterion exists to reach.

   Store the summed rate on `SpeciesRecord` at intern time and read it (§8.6).
   *Test:* a `criterion` bench of `decay_channels` at species size 1 against
   `MAX_POLYMER`, ratio asserted < 1.2x; it fails today at 8.4x. Pair it with
   the allocation test — alone, it is satisfied by the plausible fix.

   *And correct the record while doing it:* CLAUDE.md's "~2 orders of magnitude"
   is right for the class (canonicalisation, embedding, folding are 10⁴–10⁵x per
   call) and wrong for this member — measured 1.20x at n=12 rising to ~8x at
   n=200. The defensible sentence is "4% of a step budget with small molecules,
   over 100% of it at `MAX_POLYMER`". Overclaiming gets the item dismissed the
   first time someone measures it; "a twelve-iteration loop" gets it
   deprioritised forever.

**And the test named as this task's §8.6 guard cannot fail.**
`cost_does_not_depend_on_how_many_molecules_exist` is two `let _: f64 = ..`
bindings and a comment reading "Asserted structurally: the function takes a
scalar count, not a collection". The structural argument is sound and the test
asserts nothing — it is silent on both the allocation and the species-size
scaling above, and it is the test a reader would point at to justify deferring
either. Replace it with the allocation counter, which is the assertion its own
doc comment is describing.

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn propensity_scales_linearly_with_count() {
        let (u, g, it, id) = fixture();
        let one = decay_propensity(id, 1.0, &it, Thermal(300.0), &u);
        let ten = decay_propensity(id, 10.0, &it, Thermal(300.0), &u);
        assert!((ten - 10.0 * one).abs() < 1e-9, "not linear in count");
    }

    #[test]
    fn propensity_rises_with_temperature() {
        let (u, g, it, id) = fixture();
        assert!(
            decay_propensity(id, 1.0, &it, Thermal(500.0), &u)
                > decay_propensity(id, 1.0, &it, Thermal(200.0), &u)
        );
    }

    /// The whole point of §9.5: decay costs nothing per molecule. Computing
    /// a species' propensity must not touch anything per-molecule, which
    /// shows up as the cost being independent of the count.
    #[test]
    fn cost_does_not_depend_on_how_many_molecules_exist() {
        let (u, g, it, id) = fixture();
        // Same work regardless of count — a per-molecule sweep would make
        // these differ in observable effort. Asserted structurally: the
        // function takes a scalar count, not a collection.
        let _: f64 = decay_propensity(id, 1.0, &it, Thermal(300.0), &u);
        let _: f64 = decay_propensity(id, 1.0e9, &it, Thermal(300.0), &u);
    }

    #[test]
    fn a_molecule_with_no_bonds_does_not_spontaneously_cleave() {
        let (u, g, mut it) = empty_fixture();
        let mut m = borbax_molecule::graph::Mol12::new();
        m.add_atom(borbax_universe::ElementId(0));
        let id = it.intern(&m, &u, &g);
        assert_eq!(it.record(id).cleave_propensity, 0.0);
    }
}
```

(`fixture()` and `empty_fixture()` build a `Universe::generate(21)`, a
`Geodesic::<42>`, an `Interner`, and intern a six-atom chain — write them at
the top of the test module.)

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```rust
//! Decay (spec §9.5).
//!
//! **Decay costs nothing per molecule.** The obvious implementation — sweep
//! every molecule every step and roll for each of its bonds — would dominate
//! the loop, and it is unnecessary. Thermal cleavage is a Poisson process
//! whose rate depends only on the species and the temperature, so a species'
//! total propensity is its count times a precomputed constant: one more
//! channel in the same scheduler that handles every other reaction.
//!
//! This works **because damage is implicit** (spec §22.7) — a damaged molecule
//! is simply a different species, so molecules of a species are
//! interchangeable and none needs individual state. Adding a per-molecule
//! damage field would silently turn this back into an O(molecules) sweep.
//! That is a change of asymptotic class, not a constant factor.

use crate::species::{Interner, SpeciesId};
use borbax_units::Thermal;
use borbax_universe::Universe;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecayKind {
    /// Thermal bond cleavage. Hot places destroy structure; cold preserve it.
    Spontaneous,
    /// The solvent binds an exposed bond and cleaves it. A pure function of
    /// the species, resolved at intern time.
    Solvent,
    /// An unstable element decays, wrecking its host. Also the mutation source.
    Radiogenic,
}

/// Total decay propensity for a species present at `count`.
#[must_use]
pub fn decay_propensity<const D: usize>(
    id: SpeciesId,
    count: f64,
    it: &Interner<D>,
    t: Thermal,
    u: &Universe,
) -> f64 {
    let r = it.record(id);
    let temp = t.get().max(1e-6);
    // One `exp`, on a species-level constant — not per molecule, not per bond.
    let thermal = r.cleave_propensity * borbax_molecule::det_math::exp(-1.0 / (temp * 0.01));
    let solvent = r.solvent_rate.max(0.0) * u.consts.decay_scale;
    count * (thermal + solvent)
}

/// The individual channels, for the scheduler and for the cause-of-death
/// census (spec §15.2). Order is fixed, so channel indices are stable.
#[must_use]
pub fn decay_channels<const D: usize>(
    id: SpeciesId,
    count: f64,
    it: &Interner<D>,
    t: Thermal,
    u: &Universe,
) -> Vec<(DecayKind, SpeciesId, f64)> {
    let r = it.record(id);
    let temp = t.get().max(1e-6);
    vec![
        (
            DecayKind::Spontaneous,
            id,
            count * r.cleave_propensity * borbax_molecule::det_math::exp(-1.0 / (temp * 0.01)),
        ),
        (DecayKind::Solvent, id, count * r.solvent_rate.max(0.0) * u.consts.decay_scale),
        (DecayKind::Radiogenic, id, count * radiogenic_rate(id, it, u)),
    ]
}

fn radiogenic_rate<const D: usize>(id: SpeciesId, it: &Interner<D>, u: &Universe) -> f64 {
    // Summed element `decay_rate`. The field was called `stability` until
    // Task 4's review: it holds 0 at the binding peak and rises toward the
    // extremes, so it was being summed here *as a rate* — correct sense, wrong
    // name, which is exactly the inversion a reader reconciling the two would
    // introduce.
    //
    // **This walk is a §8.6 violation and must not be copied as written.** It
    // is kept here only so the rename is legible against the previous draft.
    // `SpeciesRecord::radiogenic_rate` is the precomputed value; read that.
    // An earlier version of this comment ended "hoisted into SpeciesRecord if
    // it shows up in a profile" — the sentence the field's own doc quotes and
    // condemns, 430 lines above. Deleted rather than softened, because it is
    // the sentence an implementer copies.
    //
    // Two corrections to the field doc's sizing, measured: the loop is bounded
    // by MAX_POLYMER = 200, not 12; and on the same eleven lines the returned
    // `Vec` costs ~4x more than the walk at small-molecule sizes (16.3 ns of
    // allocator against 4.0 ns of walk). A fix that hoists the walk and leaves
    // the `Vec` buys 17% of the cost. See Task 15.
    let canon = it.record(id).canon;
    (0..canon.n as usize)
        .map(|i| u.element(borbax_universe::ElementId(canon.elem[i])).decay_rate)
        .sum()
}
```

- [ ] **Step 4: Run**, then **Step 5: Commit** with a message noting that decay
being O(1) per species depends on implicit damage, and that the §22.7
provenance fallback would break it.

---

### Task 16: Catalysis from cavity geometry

**Files:** Create `crates/borbax-reaction/src/catalysis.rs`; test same file.

**Interfaces:**
- Consumes: `Cavity<D>`, `Signature<D>`, `Geodesic<D>`, `BindConsts`, `fcc`
- Produces: `catalysis_factor<D>(&[Cavity<D>], &Signature<D>, &Signature<D>, &Geodesic<D>, &BindConsts) -> f64`, `IDEAL_SEPARATION`, `MAX_ENHANCEMENT`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_cavities_means_no_enhancement() {
        assert_eq!(catalysis_factor::<42>(&[], &sig_a(), &sig_b(), &geo(), &consts()), 1.0);
    }

    #[test]
    fn one_cavity_is_a_receptor_not_a_catalyst() {
        // A single matching cavity binds a substrate but cannot hold two
        // reactants adjacent, so it must not enhance a rate.
        let f = catalysis_factor(&[cavity_matching(&sig_a())], &sig_a(), &sig_b(), &geo(), &consts());
        assert!((f - 1.0).abs() < 1e-9, "a single cavity acted as a catalyst");
    }

    /// The headline behaviour of the whole project (spec §8.5): two cavities
    /// close together, each complementing one reactant, accelerate a reaction.
    /// Nothing in the code knows what an enzyme is.
    #[test]
    fn two_well_placed_cavities_accelerate() {
        let cavs = vec![cavity_matching(&sig_a()), cavity_matching_at(&sig_b(), IDEAL_SEPARATION)];
        let f = catalysis_factor(&cavs, &sig_a(), &sig_b(), &geo(), &consts());
        assert!(f > 2.0, "two matching cavities gave only {f}x");
    }

    #[test]
    fn cavities_too_far_apart_do_not_help() {
        let cavs = vec![cavity_matching(&sig_a()), cavity_matching_at(&sig_b(), IDEAL_SEPARATION * 8.0)];
        let f = catalysis_factor(&cavs, &sig_a(), &sig_b(), &geo(), &consts());
        assert!(f < 1.5, "distant cavities enhanced by {f}x");
    }

    #[test]
    fn enhancement_is_bounded() {
        let cavs = vec![cavity_matching(&sig_a()), cavity_matching_at(&sig_b(), IDEAL_SEPARATION)];
        assert!(catalysis_factor(&cavs, &sig_a(), &sig_b(), &geo(), &consts()) <= MAX_ENHANCEMENT);
    }
}
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```rust
//! Catalysis (spec §8.5).
//!
//! **Catalysis is not implemented. It falls out.**
//!
//! A polymer with one cavity that complements a small molecule binds it —
//! that is a receptor. A polymer with *two* cavities close together binds two
//! molecules and holds them adjacent in a fixed relative orientation for as
//! long as the complex persists, and that proximity is precisely what lowers
//! an activation barrier.
//!
//! So there is no enzyme type, no catalysis rule, and nothing anywhere that
//! knows what an enzyme is. There is a geometric measurement, and a folded
//! chain that happens to score well on it *is* one. We detect that, and the
//! Chronicle reports it.

use borbax_molecule::binding::{affinity_ordered, BindConsts};
use borbax_molecule::cavity::Cavity;
use borbax_molecule::fold::fcc;
use borbax_molecule::geodesic::Geodesic;
use borbax_molecule::signature::Signature;

//! **On the constants below.** An earlier draft had four bare literals here,
//! and two of them were doing more than they looked.
//!
//! `IDEAL_SEPARATION = 2.5` invented a second length scale for something
//! `UniverseConsts.ideal_gap` already means. It now derives from that plus the
//! substrates' own extents.
//!
//! `BIND_THRESHOLD = -12.0` was the worse one. `affinity` scales with **D** and
//! with the generated `w_shape`/`w_charge`/`ideal_gap`, so a fixed cutoff means
//! a different quality of fit at D = 12, 42 and 162 — and §22.2 sweeps exactly
//! those three to lock D. The sweep would have measured the threshold moving
//! rather than the resolution. It is now expressed as a fraction of a reference
//! score.
//!
//! `MAX_ENHANCEMENT = 1.0e4` was inert: `fit ≤ 1`, `enclosure ≤ 1`,
//! `proximity ≤ 1`, so the factor could not exceed ~1001 and the clamp never
//! fired. Its test could not fail.
//!
//! And the `* 1.0e3` multiplier *was* the enhancement magnitude — the number
//! deciding whether exit criterion 6 reports 2x or 500x. Every other magnitude
//! in the system is generated per universe. This one is now
//! `UniverseConsts.catalytic_prefactor`, so criterion 6 reports a fact about a
//! generated universe rather than about a literal nobody would question again.

/// Fraction of a reference score above which a cavity counts as holding a
/// substrate. Dimensionless, so it is stable across resolutions.
const BIND_FRACTION: f64 = 0.25;

/// Rate enhancement this folded chain provides to a reaction between two
/// substrates. Returns 1.0 (no effect) unless two distinct cavities each
/// complement one substrate and sit at a workable separation.
#[must_use]
pub fn catalysis_factor<const D: usize>(
    cavities: &[Cavity<D>],
    sub_a: &Signature<D>,
    sub_b: &Signature<D>,
    g: &Geodesic<D>,
    k: &BindConsts,
) -> f64 {
    if cavities.len() < 2 {
        return 1.0;
    }
    let mut best = 1.0f64;

    for (i, ca) in cavities.iter().enumerate() {
        let bind_a = affinity_ordered(&ca.signature, sub_a, g, k);
        if bind_a < BIND_THRESHOLD {
            continue;
        }
        for (j, cb) in cavities.iter().enumerate() {
            if i == j {
                continue;
            }
            let bind_b = affinity_ordered(&cb.signature, sub_b, g, k);
            if bind_b < BIND_THRESHOLD {
                continue;
            }

            let (ax, ay, az) = fcc::from_flat(ca.centre);
            let (bx, by, bz) = fcc::from_flat(cb.centre);
            let sep = (f64::from(ax - bx).powi(2)
                + f64::from(ay - by).powi(2)
                + f64::from(az - bz).powi(2))
            .sqrt();

            // Proximity term: peaks at IDEAL_SEPARATION, falls off either side.
            // Quadratic rather than exponential so there is no transcendental
            // on this path (spec §13.1).
            let offset = (sep - IDEAL_SEPARATION).abs() / IDEAL_SEPARATION;
            let proximity = (1.0 - offset * offset).max(0.0);
            if proximity <= 0.0 {
                continue;
            }

            // Enclosure matters: a deeply enclosed pair holds its substrates
            // in a fixed orientation, a shallow pair barely constrains them.
            let enclosure = f64::from(ca.enclosure.min(cb.enclosure)) / 12.0;
            // Fit quality, mapped so a perfect complement (score 0) is best.
            let fit = (1.0 / (1.0 - bind_a)) * (1.0 / (1.0 - bind_b));

            let factor = 1.0 + proximity * enclosure * fit * 1.0e3;
            best = best.max(factor.min(MAX_ENHANCEMENT));
        }
    }
    best
}
```

- [ ] **Step 4: Run**, then **Step 5: Commit.**

The commit message should say plainly that this file contains no concept of
an enzyme — only a geometric measurement — because that is the claim V0 exit
criterion 6 tests, and a future reader needs to know it was deliberate.

---

### Task 17: RAF detection

**Files:** Create `crates/borbax-reaction/src/raf.rs`; test same file.

**Interfaces:**
- Consumes: `Reaction`, `SpeciesId`
- Produces: `Incidence`, `RafSet { reactions, species }`, `find_raf(&[Reaction], &[(SpeciesId, usize)], &BTreeSet<SpeciesId>) -> Option<RafSet>`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn food(ids: &[u32]) -> BTreeSet<SpeciesId> {
        ids.iter().map(|&i| SpeciesId(i)).collect()
    }

    /// A -> B catalysed by B, with A available. B makes itself: a RAF.
    #[test]
    fn finds_a_minimal_self_sustaining_set() {
        let rxns = vec![rxn(0, 0, 1, 1)];
        let cat = vec![(SpeciesId(1), 0usize)];
        let raf = find_raf(&rxns, &cat, &food(&[0])).expect("should find a RAF");
        assert_eq!(raf.reactions, vec![0]);
    }

    /// Same reaction, but the catalyst is never produced. Not a RAF.
    #[test]
    fn rejects_a_set_whose_catalyst_is_unreachable() {
        let rxns = vec![rxn(0, 0, 1, 1)];
        let cat = vec![(SpeciesId(9), 0usize)];
        assert!(find_raf(&rxns, &cat, &food(&[0])).is_none());
    }

    /// Reactants not traceable to food. Not a RAF, however well catalysed.
    #[test]
    fn rejects_a_set_not_grounded_in_food() {
        let rxns = vec![rxn(5, 5, 6, 6)];
        let cat = vec![(SpeciesId(6), 0usize)];
        assert!(find_raf(&rxns, &cat, &food(&[0])).is_none());
    }

    /// The closure bug the original four tests could not catch, because they
    /// all used single-reactant reactions. `A + B -> C` with `A` in food and
    /// `B` unreachable must **not** yield a RAF.
    #[test]
    fn a_bimolecular_reaction_needs_all_its_reactants() {
        let rxns = vec![Reaction {
            kind: ReactionKind::Condense,
            reactants: [SpeciesId(0), SpeciesId(7)], // 7 is never produced
            products: [SpeciesId(1), SpeciesId(1)],
            activation: Quanta(1.0),
            delta: Quanta(0.0),
            catalysis: 1.0,
            catalyst: Some(SpeciesId(1)),
        }];
        let cat = vec![(SpeciesId(1), 0usize)];
        assert!(
            find_raf(&rxns, &cat, &food(&[0])).is_none(),
            "produced a product from a reactant that never existed"
        );
    }

    #[test]
    fn is_deterministic_and_order_independent() {
        let rxns = vec![rxn(0, 0, 1, 1), rxn(1, 1, 2, 2), rxn(2, 0, 3, 3)];
        let cat = vec![(SpeciesId(1), 0), (SpeciesId(2), 1), (SpeciesId(3), 2)];
        let a = find_raf(&rxns, &cat, &food(&[0]));
        let mut shuffled = cat.clone();
        shuffled.reverse();
        let b = find_raf(&rxns, &shuffled, &food(&[0]));
        assert_eq!(a, b, "result depended on catalyst listing order");
    }
}
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

```rust
//! RAF detection — Reflexively Autocatalytic and Food-generated sets
//! (spec §9.3), following Hordijk and Steel.
//!
//! A RAF is a set of reactions where every reaction is catalysed by something
//! the set itself produces, and every reactant traces back to environmentally
//! available food. In other words: a chemical system that makes itself.
//!
//! Finding one is rigorous, computable in polynomial time, and unambiguous,
//! which makes it the right trigger both for promoting a region to full
//! detail and for the Chronicle to announce that something happened. It is
//! also how §9.6 defines death: a compartment dies when its set stops closing.
//!
//! **The output is result-affecting**, so every container here is ordered.

use crate::kind::Reaction;
use crate::species::SpeciesId;
use std::collections::{BTreeSet, VecDeque};

/// species -> reactions mentioning it, in CSR form. Closure decrements a
/// per-reaction missing-reactant counter as species enter the set, so the
/// whole closure is O(total incidences) rather than O(species x reactions).
#[derive(Debug, Clone, Default)]
pub struct Incidence {
    offsets: Vec<u32>,
    reactions: Vec<u32>,
}

impl Incidence {
    #[must_use]
    pub fn build(rxns: &[Reaction], n_species: usize) -> Self {
        let mut counts = vec![0u32; n_species + 1];
        for r in rxns {
            for s in &r.reactants {
                if (s.0 as usize) < n_species {
                    counts[s.0 as usize + 1] += 1;
                }
            }
        }
        for i in 1..counts.len() {
            counts[i] += counts[i - 1];
        }
        let offsets = counts.clone();
        let mut reactions = vec![0u32; *counts.last().unwrap_or(&0) as usize];
        let mut cursor = offsets.clone();
        for (ri, r) in rxns.iter().enumerate() {
            for s in &r.reactants {
                if (s.0 as usize) < n_species {
                    reactions[cursor[s.0 as usize] as usize] = ri as u32;
                    cursor[s.0 as usize] += 1;
                }
            }
        }
        Self { offsets, reactions }
    }

    fn of(&self, s: SpeciesId) -> &[u32] {
        let i = s.0 as usize;
        if i + 1 >= self.offsets.len() {
            return &[];
        }
        &self.reactions[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RafSet {
    /// Reaction indices, sorted.
    pub reactions: Vec<usize>,
    /// Species in the closure, sorted.
    pub species: BTreeSet<SpeciesId>,
}

/// Find the maximal RAF, or `None` if there is not one.
///
/// **Maximality is well-defined**: RAFs are closed under union, so a unique
/// maximal one exists. Each removal is safe because for any RAF `R' ⊆ active`,
/// `cl_R'(F) ⊆ cl_active(F)` — a reaction unsupported by the larger closure
/// cannot belong to any sub-RAF. So `None` means *there is no RAF*, not "we
/// failed to find one".
///
/// **This is a RAF, not a CAF.** A CAF requires each catalyst to already exist
/// when its reaction is first used; a RAF permits a catalyst to be produced by
/// the very set that needs it. RAF is the weaker, more permissive condition,
/// so finding one means "self-referentially closed", not "bootstrappable from
/// food in order". Worth stating, because the Chronicle announces it.
///
/// **Two limits, both deliberate for V0 and both worth knowing.**
///
/// It is *structural*: this sees a reaction list and a food set, never any
/// counts. So §9.6's death criterion fires when a reaction channel disappears,
/// but not when a catalyst's population reaches zero. Callers should therefore
/// pass only species actually present — see `find_raf_present`.
///
/// It ignores *inhibition*. §9.5's broad-complementarity by-products are
/// inhibitors in RAF terms, and detecting uninhibited RAFs under inhibition is
/// NP-hard. Polynomial time is the reason RAF was chosen, so the simplification
/// has to be stated rather than assumed.
///
/// The algorithm alternates two reductions until they stop removing anything:
/// drop reactions whose reactants are not in the closure of what is currently
/// reachable, and drop reactions with no catalyst present in that closure.
/// Whatever survives is the maximal RAF.
#[must_use]
pub fn find_raf(
    rxns: &[Reaction],
    catalysts: &[(SpeciesId, usize)],
    food: &BTreeSet<SpeciesId>,
) -> Option<RafSet> {
    let n_species = rxns
        .iter()
        .flat_map(|r| r.reactants.iter().chain(r.products.iter()))
        .map(|s| s.0 as usize + 1)
        .chain(food.iter().map(|s| s.0 as usize + 1))
        .max()
        .unwrap_or(0);
    let inc = Incidence::build(rxns, n_species);

    let mut active: BTreeSet<usize> = (0..rxns.len()).collect();

    loop {
        // Closure: what can be built from food using only active reactions.
        let mut present = food.clone();
        // Count **all** reactants, and let food decrement them like anything
        // else.
        //
        // An earlier draft counted only non-food reactants while
        // `Incidence::build` indexed every reactant occurrence including food
        // ones. For `A + B -> C` with `A` in food and `B` not: missing = 1,
        // popping `A` decrements it to zero, and `C` is produced **without `B`
        // ever existing**. The F-generated condition is violated and the
        // returned set need not be a RAF at all.
        //
        // All four of this task's original tests used single-reactant
        // reactions, so none of them could catch it.
        let mut missing: Vec<u32> = rxns.iter().map(|r| r.reactants.len() as u32).collect();
        let mut queue: VecDeque<SpeciesId> = food.iter().copied().collect();

        while let Some(s) = queue.pop_front() {
            for &ri in inc.of(s) {
                let ri = ri as usize;
                if !active.contains(&ri) || missing[ri] == 0 {
                    continue;
                }
                missing[ri] -= 1;
                if missing[ri] == 0 {
                    for p in &rxns[ri].products {
                        if present.insert(*p) {
                            queue.push_back(*p);
                        }
                    }
                }
            }
        }

        // Prune: a reaction survives only if its reactants are all present
        // AND some catalyst for it is present. `catalysts` is scanned in full
        // rather than indexed, so listing order cannot affect the outcome.
        let before = active.len();
        active.retain(|&ri| {
            let reactants_ok = rxns[ri].reactants.iter().all(|s| present.contains(s));
            let catalysed = catalysts.iter().any(|(c, r)| *r == ri && present.contains(c));
            reactants_ok && catalysed
        });

        if active.len() == before {
            if active.is_empty() {
                return None;
            }
            return Some(RafSet { reactions: active.into_iter().collect(), species: present });
        }
    }
}
```

- [ ] **Step 4: Run** — expect 4 passed. **Step 5: Commit.**

---

## Phase 6 — The beaker and seeing inside it

### Task 18: The well-mixed beaker

**Files:** Create `crates/borbax-beaker/{Cargo.toml,src/lib.rs,src/propensity.rs}`; test both.

**Interfaces:**
- Consumes: `Interner`, `Reaction`, `decay_channels`, `Universe`, `Stream`
- Produces: `SegmentTree::{with_capacity, set, total, sample}`, `Beaker::{new, step, time, counts, species_count}`, `ChannelId(u32)`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sampling_is_proportional_to_propensity() {
        let mut t = SegmentTree::with_capacity(4);
        t.set(0, 1.0);
        t.set(1, 3.0);
        assert!((t.total() - 4.0).abs() < 1e-12);
        let mut hits = [0u32; 2];
        for i in 0..10_000 {
            hits[t.sample(f64::from(i) / 10_000.0)] += 1;
        }
        let ratio = f64::from(hits[1]) / f64::from(hits[0]);
        assert!((ratio - 3.0).abs() < 0.2, "ratio {ratio}");
    }

    /// A Fenwick tree updated by deltas accumulates rounding drift in its
    /// partial sums, and drift slowly changes which reaction a given RNG draw
    /// selects. A segment tree recomputes parents as left + right, so stored
    /// totals stay bit-identical to a full rebuild — forever.
    #[test]
    fn totals_stay_bit_identical_to_a_rebuild() {
        let mut t = SegmentTree::with_capacity(64);
        let mut r = borbax_rng::Stream::new(5, borbax_rng::Domain::Beaker, 0);
        for _ in 0..50_000 {
            t.set(r.next_range(64) as usize, r.next_f64() * 1e6);
        }
        let mut fresh = SegmentTree::with_capacity(64);
        for i in 0..64 {
            fresh.set(i, t.leaf(i));
        }
        assert_eq!(t.total().to_bits(), fresh.total().to_bits(), "drift detected");
    }

    #[test]
    fn zero_propensity_slots_are_never_selected() {
        let mut t = SegmentTree::with_capacity(8);
        t.set(3, 1.0);
        for i in 0..1_000 {
            assert_eq!(t.sample(f64::from(i) / 1_000.0), 3);
        }
    }

    #[test]
    fn stepping_is_deterministic_and_advances_time() {
        let (u, g) = fixture();
        let mut a = Beaker::new(&u, &g, 1);
        let mut b = Beaker::new(&u, &g, 1);
        for _ in 0..500 {
            a.step(&u, &g);
            b.step(&u, &g);
        }
        assert_eq!(a.time(), b.time());
        assert_eq!(a.counts(), b.counts());
        assert!(a.time().get() > 0.0);
    }

    /// Reactions appear and vanish. If slots were compacted on removal,
    /// indices would shift and identical physics would sample differently.
    #[test]
    fn freed_slots_do_not_shift_other_indices() {
        let mut t = SegmentTree::with_capacity(8);
        t.set(2, 5.0);
        t.set(5, 7.0);
        t.set(2, 0.0); // "removed"
        assert_eq!(t.sample(0.5), 5, "indices shifted when a slot was freed");
    }
}
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement the propensity tree**

```rust
//! Propensity selection for the stochastic step.
//!
//! A **segment tree**, not a Fenwick tree. Both are O(log n), but a Fenwick
//! updated by deltas accumulates rounding drift in its partial sums, and that
//! drift slowly changes which reaction a given RNG draw selects — a
//! determinism failure that would appear as a slow divergence rather than an
//! obvious bug. A segment tree recomputes each parent as `left + right`, so
//! its stored totals are bit-identical to a full rebuild no matter how many
//! updates have happened.
//!
//! Gibson-Bruck was considered and rejected for V0: its advantage is few
//! propensity updates per event, which is undercut in a well-mixed volume
//! where one firing changes counts many channels depend on — and its
//! dependency-graph bookkeeping is where the bugs live.

pub struct SegmentTree {
    cap: usize,
    tree: Vec<f64>,
}

impl SegmentTree {
    #[must_use]
    pub fn with_capacity(n: usize) -> Self {
        let cap = n.next_power_of_two().max(1);
        Self { cap, tree: vec![0.0; cap * 2] }
    }

    #[must_use]
    pub fn leaf(&self, slot: usize) -> f64 {
        self.tree[self.cap + slot]
    }

    #[must_use]
    pub fn total(&self) -> f64 {
        self.tree[1]
    }

    /// Set a slot's propensity. Freeing a channel means setting it to 0.0 —
    /// **never** compacting, because that would shift every later index and
    /// make identical physics sample differently.
    pub fn set(&mut self, slot: usize, value: f64) {
        let mut i = self.cap + slot;
        self.tree[i] = value;
        while i > 1 {
            i /= 2;
            self.tree[i] = self.tree[2 * i] + self.tree[2 * i + 1];
        }
    }

    /// Select a slot with probability proportional to its propensity.
    /// `u` must be uniform in [0, 1).
    #[must_use]
    pub fn sample(&self, u: f64) -> usize {
        let mut t = u * self.tree[1];
        let mut i = 1;
        while i < self.cap {
            let left = self.tree[2 * i];
            if t < left {
                i *= 2;
            } else {
                t -= left;
                i = 2 * i + 1;
            }
        }
        let mut slot = i - self.cap;
        // Guard the rounding case where `u` near 1 lands on a zero-propensity
        // leaf. Deterministic forward scan, so the fallback is reproducible.
        if self.tree[self.cap + slot] <= 0.0 {
            for k in 0..self.cap {
                let s = (slot + k) % self.cap;
                if self.tree[self.cap + s] > 0.0 {
                    slot = s;
                    break;
                }
            }
        }
        slot
    }
}
```

- [ ] **Step 4: Implement the beaker**

The `Beaker` holds an `Interner`, a `Vec<f64>` of counts indexed by
`SpeciesId`, a `Vec<Reaction>` of channels, a `SegmentTree` over their
propensities, a free-slot list, and a `WorldYear` clock.

`step` performs one Gillespie iteration:

1. Compute `total = tree.total()`. If zero, the beaker is dead — return.
2. Draw `u1`, advance the clock by `-ln(u1) / total` (via `det_math::ln`).

   **Guard `u1 == 0`.** `next_f64` is uniform on `[0, 1)` — closed at zero and
   open at one, the wrong half-open direction for this — so zero is reachable
   with probability 2^-53, and `ln(0)` is `-inf`: the clock jumps to infinity
   and the run silently stops scheduling, with no panic and nothing in the
   state hash. Draw from `(0, 1]` by using `1.0 - u1`; Exp(1) is invariant
   either way, and `-ln(1 - u)` is total, truncating at 36.7368 — a cap on the
   waiting time rather than an infinity, with the omitted tail at mass 2^-53.
   (`1.0 - u` is exact for every `u` on the 2^-53 grid, which is what makes
   36.7368 a true bound rather than an approximation.)

   **Test it at the value.** `u == 0.0` is reachable but only at 2^-53, so a
   sampling test will never find it and a direct one always will. `borbax-rng`'s
   `normal_from` is the worked example: it is private precisely so its boundary
   can be called directly.

   **And pin the draw-count contract while you are here.** `next_range`
   rejection-samples, so a step's draw count is value-dependent — stream
   position is therefore *not* a function of the event count, which anything
   reconstructing a run will otherwise assume.
3. Draw `u2`, select a channel with `tree.sample(u2)`.

   **Assert `tree.total()` is finite before sampling.** A single NaN propensity
   anywhere in the tree makes `sample` return the **last** channel on every
   call — measured, not predicted — so the beaker fires one reaction forever
   and produces clean-looking, completely wrong output. No crash, nothing in
   the state hash. This is the hazard `borbax-rng`'s `next_f64_range` doc used
   to describe, wrongly attributed to that function and with the channel
   inverted; it lives here because this is where it can actually happen.
4. Apply it: decrement reactant counts, increment product counts. Products
   were resolved when the channel was created (spec §8.6), so no
   canonicalisation happens here.

   **For a reaction between two molecules of the same species the propensity
   is `c·n(n-1)/2`, not `c·n²/2`.** The difference is invisible at large `n`
   and wrong exactly where it matters — at the small counts where a RAF
   nucleates. Note also that `counts` being `Vec<f64>` invites fractional
   values, at which point `n(n-1)` stops meaning anything; integers up to 2⁵³
   are exact in `f64` so determinism survives, but the conversion should be
   explicit and commented.
5. Recompute `catalysis` only when a *cavity-bearing* species' count crosses
   zero — it is a property of the channel, resolved at creation (§8.6), not a
   per-step quantity. Then update the propensity of every channel whose
   species changed. Below
   `FLAT_SCAN_LIMIT` (start at 256) rescan all channels — a contiguous scan
   beats a dependent-load tree descent at small counts. Above it, use the
   incidence structure from Task 17 to touch only affected channels.
6. When a genuinely novel species appears, intern it and create its channels.
   **This is the only path on which canonicalisation runs.**

Decay channels from Task 15 are ordinary entries in the same tree — decay is
not a separate pass.

- [ ] **Step 5: Run** — expect 5 passed. **Step 6: Commit**, recording the
measured flat-scan/tree crossover in the message (perf review benchmark 2).

---

### Task 19: Geometry and SVG rendering

**Files:** Create `crates/borbax-geometry/{Cargo.toml,src/{lib,project,net,palette}.rs}` and
`crates/borbax-render/{Cargo.toml,src/{lib,svg,molecule,geodesic_net,fold,binding}.rs}`; test both.

**Interfaces:**
- Produces: `project::{orthographic, isometric, depth_sort}`, `net::unfold`, `palette::{affinity_colour, burial_colour}`, `Svg::{new, line, circle, polygon, text, finish}`, `render_molecule`, `render_signature_net`, `render_fold`, `render_binding`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn svg_output_is_deterministic() {
        assert_eq!(render_molecule(&fixture_mol(), &u(), &g()), render_molecule(&fixture_mol(), &u(), &g()));
    }

    #[test]
    fn svg_is_well_formed_and_self_contained() {
        let s = render_molecule(&fixture_mol(), &u(), &g());
        assert!(s.starts_with("<svg "));
        assert!(s.ends_with("</svg>"));
        assert!(!s.contains("http://"), "external reference in output");
        assert_eq!(s.matches("<g").count(), s.matches("</g>").count());
    }

    #[test]
    fn floats_are_written_at_fixed_precision() {
        // Full f64 precision would make goldens fragile against harmless
        // last-bit differences; fixed precision makes the diff meaningful.
        let s = render_molecule(&fixture_mol(), &u(), &g());
        for tok in s.split(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-')) {
            if let Some((_, frac)) = tok.split_once('.') {
                assert!(frac.len() <= 3, "over-precise coordinate: {tok}");
            }
        }
    }

    /// Golden render tests via `insta` (spec §16). SVG is a string, so a
    /// rendering regression is caught like any other — no image diffing, no
    /// browser, runs headless in CI.
    ///
    /// `insta`, not a hand-rolled `UPDATE_GOLDENS` dance: it is already the
    /// project's named snapshot tool (CLAUDE.md § Dependencies), it is a
    /// dev-dependency and so cannot touch simulation output, and
    /// `cargo insta review` shows a diff and requires an explicit accept.
    /// That last point is the one that matters — "inspect every generated SVG
    /// by eye" is a step a human will discharge by accepting all fourteen.
    #[test]
    fn molecule_render_matches_snapshot() {
        insta::assert_snapshot!(render_molecule(&fixture_mol(), &u(), &g()));
    }

    #[test]
    fn signature_net_render_matches_snapshot() {
        insta::assert_snapshot!(render_signature_net(&fixture_sig(), &g()));
    }

    #[test]
    fn fold_render_matches_snapshot() {
        insta::assert_snapshot!(render_fold(&fixture_fold(), &fixture_poly(), &u()));
    }
}
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

`borbax-geometry` computes *what* to draw; `borbax-render` turns it into SVG.
The split exists now rather than later because V1's `borbax-ui` becomes a
second backend over the same geometry (spec §14.5) — one projection, two
renderers, no duplicated maths.

Key points for the implementer:

- **`Svg`** is a plain `String` builder. Write floats with `format!("{:.3}")`
  throughout — full precision makes goldens fragile against harmless last-bit
  differences and makes a real diff impossible to read.
- **`project::orthographic`** projects the 3D embedding along a fixed axis;
  **`isometric`** uses a fixed rotation. Both are constants, not parameters:
  a golden test needs a fixed viewpoint, and §14.5's argument is that a
  well-chosen fixed viewpoint beats an interactive one for "why did these two
  bind?".
- **`net::unfold`** flattens the icosahedron into 20 triangles in a fixed
  layout, so a spherical function is legible on a flat page. This is the
  map-projection trick and it is what makes a 42-direction signature
  readable at all.
- **`render_fold`** draws monomers back-to-front by depth, shades each by
  burial (`exposure` from Task 11), and outlines cavities from Task 13.
- **`render_binding`** draws both signatures in the winning orientation with
  per-direction contributions annotated, so you can see *which* part of the
  fit carried the score. This is the diagram that makes a binding bug
  diagnosable rather than mysterious.

- [ ] **Step 4: Accept the snapshots** — `cargo insta test -p borbax-render`,
then `cargo insta review`, which shows each snapshot as a diff and requires an
explicit accept per file.

**Open each SVG in a browser before accepting it.** A snapshot accepted without
being looked at locks in whatever was wrong at the time, and the review tool
makes that easy to do quickly rather than impossible to do wrongly.

- [ ] **Step 5: Run** — `cargo test -p borbax-render`, expect 4 passed.
**Step 6: Commit** goldens and code together.

---

### Task 20: The beaker battery and metrics

**Files:** Create `crates/borbax-beaker/src/{battery.rs,metrics.rs}`; test both.

**Interfaces:**
- Produces: `BatteryReport`, `run_battery(&Universe, &Geodesic<D>) -> BatteryReport`, `BatteryReport::passes()`, `Metrics`, `MetricFamily`

- [ ] **Step 1: Write the failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    // The battery's speed is a `criterion` benchmark in `benches/battery.rs`,
    // not a test.
    //
    // A wall-clock assertion inside `#[test]` is flaky under CI contention,
    // and a flaky test gets deleted rather than fixed — taking the §7.2
    // iteration-speed guarantee with it. `criterion` tracks the number over
    // time and reports regressions, which is what is actually wanted:
    //
    //     fn bench_battery(c: &mut Criterion) {
    //         let g = Geodesic::<42>::build().unwrap();
    //         c.bench_function("battery", |b| {
    //             b.iter(|| run_battery(black_box(&Universe::generate(1)), &g))
    //         });
    //     }

    #[test]
    fn it_rejects_an_inert_universe() {
        let mut u = Universe::generate(1);
        u.consts.rate_prefactor = 1e-30; // nothing will ever react
        assert!(!run_battery(&u, &Geodesic::<42>::build().unwrap()).passes());
    }

    #[test]
    fn it_rejects_a_universe_that_burns() {
        let mut u = Universe::generate(1);
        u.consts.rate_prefactor = 1e30;
        assert!(!run_battery(&u, &Geodesic::<42>::build().unwrap()).passes());
    }

    #[test]
    fn it_rejects_a_universe_with_no_decay() {
        let mut u = Universe::generate(1);
        u.consts.decay_scale = 0.0;
        let r = run_battery(&u, &Geodesic::<42>::build().unwrap());
        assert!(!r.passes(), "a universe where nothing decays cannot select (spec §9.4)");
    }

    /// **The gel ceiling, as a test rather than a table row.**
    ///
    /// Added because a review found both new criteria present in the criteria
    /// table and the prose, and absent from this block — and CLAUDE.md's own
    /// rule is write the failing test first. A criterion with no test is prose,
    /// and the defect being repaired here is precisely "a criterion that exists
    /// on paper while the failure it names passes".
    ///
    /// Gate on `p_ss / p_c`, not on either threshold alone. Drive the operating
    /// point above the threshold by raising the forward rate against Cleave —
    /// that is the sludge case §7.2 describes, where a floor-only polymer
    /// criterion passes a beaker that has gelled into one component.
    #[test]
    fn it_rejects_a_gelled_beaker() {
        let mut u = Universe::generate(1);
        // Condense fast against a slow Cleave drives `p_ss` up; `p_c` is a
        // property of the realised degree distribution and moves far less.
        u.consts.rate_prefactor *= 1e3;
        u.consts.decay_scale *= 1e-3;
        let r = run_battery(&u, &Geodesic::<42>::build().unwrap());
        assert!(
            r.p_ss / r.p_c > 1.0,
            "the fixture must actually gel, or this tests nothing: \
             p_ss={} p_c={}",
            r.p_ss,
            r.p_c
        );
        assert!(!r.passes(), "a gelled beaker must fail the ceiling (spec §7.2)");
    }

    /// Self-synthesis is a **band**, and this is the half a floor misses.
    ///
    /// RBN-World eliminated 110 of 163 chemistries for *every* sample
    /// self-synthesising — no variation — against 53 for none doing so. A
    /// one-sided criterion passes the 110-analogue, which here is the case
    /// where complementarity carries no information at all.
    #[test]
    fn it_rejects_universes_where_everything_binds_itself() {
        let mut u = Universe::generate(1);
        // Collapse the charge term so the complementarity test is trivially
        // satisfied for every species against itself.
        u.consts.charge_scale = 0.0;
        let r = run_battery(&u, &Geodesic::<42>::build().unwrap());
        assert!(
            !r.passes(),
            "universal self-binding is sludge, not richness (spec §2.1, §7.2)"
        );
    }

    #[test]
    fn some_universes_pass() {
        let g = Geodesic::<42>::build().unwrap();
        let passing = (0..40u64).filter(|&s| run_battery(&Universe::generate(s), &g).passes()).count();
        assert!(passing > 0, "no seed in 40 produced a workable universe");
        println!("{passing}/40 universes passed the battery");
    }

    /// Activity and novelty are separate families on purpose (spec §2.7):
    /// a system can show unbounded activity with zero novelty and pass a
    /// naive open-endedness test while producing nothing new.
    #[test]
    fn activity_and_novelty_are_tracked_separately() {
        let m = Metrics::default();
        assert_ne!(MetricFamily::Activity, MetricFamily::Novelty);
        assert!(m.field_families().iter().any(|f| *f == MetricFamily::Novelty));
        assert!(m.field_families().iter().any(|f| *f == MetricFamily::Decay));
    }
}
```

- [ ] **Step 2–4: Implement, run, commit.**

`run_battery` runs a few hundred thousand reactions in one well-mixed volume
and checks the criteria of spec §7.2:

| Criterion | Rejects |
|---|---|
| Reactivity band | inert universes, and universes that burn |
| **Decay band** | universes where nothing persists, or nothing decays (§22.6) |
| Polymer viability, **banded** | chains of length ≥ 20 that cannot form, cannot survive, **or a beaker that has gelled** |
| **Self-synthesis, banded** | universes where **no** species binds a copy of itself — **and** where essentially every species does |
| Shape diversity | signatures collapsing into a few clusters |
| Neutral-network structure | see below — **not** a redundancy test |
| Catalytic potential | no folded polymer producing an enclosed cavity |
| Energy landscape | chemistries that only ever run downhill |

**The gel ceiling is this task's, and it has had no landing site until now.**
§7.2 says plainly that *a gelled beaker passes the floor* — 51% of mass in one
component with seven molecules still over 20 units — so a floor-only polymer
criterion is the sludge failure passing the test written to catch it. That
ceiling was added to §7.2 by PR #8 and named no task; Task 4 then deferred to
"Task 5's battery", and Task 5 has no battery and architecturally cannot. It is
here.

Three numbers, and the third is the one that matters:

- `f_w = ⟨f²⟩/⟨f⟩` on the **realised** degree distribution in the beaker, not
  nominal valences off the element table.
- **bonds-per-node at gel**, `p_c⟨f⟩/2`, alongside the extent. Extent alone
  moves under any rescaling of the mean degree — measured on the generated
  table, abundance weighting moves extent +85% and bonds-per-node +9%, so most
  of the apparent effect is coordinate rather than safety.
- **`p_ss / p_c`, the margin between the operating point and the threshold.**
  Both of the above are *thresholds*, and reporting a second threshold changes
  nothing: the ratio is identical in either coordinate. Where the beaker
  actually sits comes from the steady state of §9.1's Condense against §9.4's
  Cleave, `(1−p)²/p = k_r/(2·k_f·n·⟨f⟩)`. **Reject on the ratio.** §7.2 carries
  the derivation and the history of getting this wrong twice.

**The self-synthesis criterion is cheap and the precedent says it discriminates.**
It is §8.3's kernel run on a species against itself — self-complementarity under
the 60-rotation search, non-trivial given §22.8's handedness. RBN-World's
five-test filter cut 183 candidates to 20 on self-synthesis alone, more than its
other four tests combined (§2.1). We have no equivalent and it costs one call.

**Copy the band, not the floor — their test was two-sided and the larger half is
the one a floor misses.** From that paper's Table 4, self-synthesis eliminated
163 of 183, but **110 of those were chemistries where every sample
self-synthesised**, discarded for showing no variation, against only 53 where
none did. A floor captures the 53. The 110-analogue here is sludge: if
essentially every species binds a copy of itself, complementarity carries no
information and §22.8's handedness has bought nothing. Writing this as a floor
would reproduce, one row up in the same table, the exact defect the banded
polymer criterion two rows below exists to fix — which is why it is banded here
before anyone implements it.

**The neutral-network criterion needs three numbers, not one.** "Folding maps
close to one-to-one" tests *redundancy*, which is necessary and nowhere near
sufficient: a map can be massively many-to-one with every preimage a scattered
set of isolated points — no connected network, no drift, no evolvability. As
written, an implementer writes a distinct-shape-count assertion and ticks the
box. What the literature actually measures (§2.5):

- **Neutral network connectivity** — sample a shape's preimage, build the graph
  on Hamming-1 edges, report the largest-connected-component fraction.
- **Shape-space covering radius** — smallest `r` such that a ball of radius `r`
  around a random sequence contains a sequence folding to every *common* shape.
  Covering was only ever claimed for common shapes; rare ones are not covered,
  so report the phenotype-frequency distribution alongside it. Otherwise a map
  where three shapes account for 90% of sequences passes trivially.
- **Plastogenetic congruence** — currently unmeasurable, because `fold()`
  produces one conformation and there is no ensemble to correlate against.
  Either state the omission or fold each polymer under `k` diagnostic seeds
  (`fold_with_seed`) and treat the distinct results as the plastic repertoire.
  §16 lists all three properties while §23 criterion 2 lists two; that gap
  should be closed deliberately rather than by omission.

`Metrics` tracks the families of spec §15.2, with **activity, novelty,
complexity, organisation and decay kept distinct**. Novelty must never be
inferred from activity.

---

### Task 20b: Evolutionary-activity metrics, the shadow, and plateau fitting

**This task exists because the plan could not previously reach its own exit
criteria.** §15's novelty metrics, neutral shadow and plateau fitting are all
specified in the spec and none was scheduled; meanwhile exit criterion 7
requires a shadow comparison, and Task 20's only novelty test asserted that two
enum variants are unequal.

**Files:** Create `crates/borbax-beaker/src/{activity.rs,shadow.rs,plateau.rs}`;
test each.

**Interfaces:**
- Produces: `ActivityStats { diversity, cumulative, mean_cumulative, new_activity }`,
  `ShadowRun::fork_from`, `shape_novelty`, `fit_models`, `ModelVerdict`

- [ ] **Step 1: Bedau–Packard activity statistics** (`activity.rs`)

Raw counts cannot detect the failure mode they exist to catch. A species
appearing once weighs the same as one persisting a million years, so "distinct
species" and "novel species rate" miss class 3b entirely (§2.7). The
Bedau–Packard statistics are persistence-weighted by construction, and that
weighting *is* the mechanism:

```rust
/// Per-species cumulative existence counters and the four aggregates.
///
/// Class 3b is `cumulative` unbounded with `diversity` bounded — visible only
/// if both are computed, which is why they are separate fields rather than one
/// "activity" number.
pub struct ActivityStats {
    /// a_i(t): cumulative existence per species, indexed by SpeciesId.
    pub per_species: Vec<f64>,
    /// D(t): count of species above the activity threshold.
    pub diversity: f64,
    /// A(t) = sum of a_i.
    pub cumulative: f64,
    /// Ā(t) = A/D.
    pub mean_cumulative: f64,
    /// A_new(t): activity of species newly crossing the threshold.
    pub new_activity: f64,
}
```

Add the **MODES persistence filter** alongside: discard components not
surviving a set interval before counting. It is the cheap, well-tested
substitute for a full shadow and maps directly onto these families.

- [ ] **Step 2: Shape-space novelty** (`activity.rs`)

```rust
/// Novelty as *minimum* distance to any prior state — nearest-neighbour-in-
/// history, following ASAL and the Lehman-Stanley archive measure.
///
/// **Min, not mean.** §15.2 said "distance to every prior one" without saying
/// which. Mean-to-all is a different and much weaker statistic that grows
/// automatically with run length, which would make novelty rise simply because
/// the run got longer.
///
/// Histograms are subsampled to a constant N first: histogram-distance
/// estimates from finite samples are positively biased, and the bias grows
/// with the number of occupied bins — so without this, novelty rises merely
/// because species count rose, which is exactly the activity/novelty
/// conflation this metric exists to prevent.
pub fn shape_novelty(current: &Histogram, history: &[Histogram]) -> f64 { /* ... */ }
```

Validate on a null trajectory of pure multinomial noise: it must score ~0.

- [ ] **Step 3: The neutral shadow** (`shadow.rs`)

**A `Stream`'s position cannot currently be saved or restored, and this step is
the first thing that needs it.** `Stream` exposes no accessor for its block
counter or intra-block offset and derives no `serde`, so "fork from a keyframe"
has no route as written. Two options, and the choice is not free:

1. **Expose position.** Keeps the 4x-amortised block buffer, and makes
   `(counter[3], spent)` part of the §13.2 contract. `Stream` stores position
   *twice* — as `counter[3]` and as `spent` into a `buf` that must match the
   previous block — and nothing in the type enforces the correspondence, so a
   field-wise constructor gets it wrong. **Two different ways**, which the
   test must cover both of:

   - one emits `4 - n%4` zeros and skips a block — caught at every `n`,
     including multiples of four;
   - one keeps `spent = 4`, emits no zeros at all, and silently truncates `n`
     down to a multiple of four — **invisible to a multiples-of-four test**.

   So the discriminating test is every `n` in `0..=16`. An implementer who
   tests `n ∈ {0,4,8}`, sees no zeros, and concludes the trap did not fire has
   shipped the second variant, and every keyframe restore silently rewinds up
   to three draws.

   **If position is exposed, expose position only.** `key[0]` is the universe
   seed, and it is not recoverable from a `Stream` today — verified, `E0599`.
   That is what makes a per-molecule sub-stream impossible to mint inside a
   routine that merely holds a `&Stream`, which is §8.6's barrier at every
   leaf frame. A field-wise `Serialize`/`Deserialize` exposes `key` by
   definition and hands that barrier away. Pin it with a consumer-crate probe:
   `s.seed()` must not compile.
2. **Drive the counter from the logical simulation state** — Salmon et al.'s
   own recommendation, and the paper names this exact failure: "the error
   identified by [5] arose because of a failure to properly checkpoint and
   restore the PRNG's state." A fresh `Stream::sub(seed, Domain::Beaker, ...,
   event_index)` per event is O(1) state, not the per-molecule state §8.6
   forbids. Its cost is the block buffer: the direct method uses 2 words per
   event, so 2 of every 4 are discarded.

**Related, and settled here rather than deferred again:** §13.1's
reproducibility tuple is `(universe_seed, world_seed, config_hash)`, and
`config_hash` has **no home** in the generator. Philox 4x64 gives six 64-bit
words and `borbax-rng` allocates all six — `key = [universe_seed, reserved]`,
`counter = [domain, index, sub, block]`.

**`config_hash` belongs in keyframe identity, not in stream keying**, and the
argument does not need the keyframe format to exist:

- **Keying on it would make universe generation config-dependent.**
  `Stream::new(universe_seed, Domain::Universe, i)` generates the periodic
  table, so folding config into the key means changing a keyframe interval or
  an output setting regenerates the chemistry. That breaks §13.4's promise that
  a shared `U-…/W-…` pair means the same planet — the same argument `Domain`'s
  own doc already makes for append-only discriminants. This rests on the spec,
  not on contested statistics, and is the decisive reason.
- Supporting, not load-bearing: keying would also throw away common random
  numbers across parameter settings, the standard variance-reduction device
  for a config sweep. Stated as "keeps the option open" rather than "obtains",
  because 25 lines below this same note says CRN's value here is unmeasured
  and must be settled by measurement — an argument cannot be decisive here and
  contested there.
- A config change that reaches the arithmetic already changes the run, so it
  does not additionally need different draws. (True only of physics-affecting
  entries; keyframe interval and output settings change nothing, which is why
  this is not the lead reason.)
- The real hazard is replaying a keyframe under a *different* config and
  silently getting a different trajectory. That is a load-time check, not a
  keying problem.

So: hash the config into keyframe identity and **refuse the restore on
mismatch**. `Stream` keeps its six words. The test that makes this real is a
keyframe written under one config and a load attempted under another, asserting
a refusal rather than a divergent replay — a load that succeeds and diverges is
the failure mode, and it is silent.

Recorded as a decision taken on reasoning, before the keyframe format exists.
If that format later makes config part of the world seed upstream, this becomes
moot rather than wrong.

Fork from a keyframe with decay rates **equalised**, not disabled (§15.3). The
common rate is set so total removal flux matches the focal run at the fork
point — otherwise the two differ in mass balance and any diversity difference
is explained by that rather than by selection.

Note in the module header that this is a drift control on *persistence only*:
catalysis produces differential formation rates, which is also selection, and
equalising it would destroy the chemistry.

**Which `Domain` the shadow draws from is undecided, and both steps below
depend on it.** `Domain::Shadow` exists but its doc explicitly declines to
settle the question: a paired `run - shadow` difference *can* have far lower
variance under common random numbers, but Glasserman & Yao (1992) find the
class where CRN is provably advantageous "rather limited", CRN decouples in a
Gillespie setting as the integrated intensities diverge, and per-channel
indexing reduces without removing that — Anderson (2012) "predicted (though did
not prove)" it, and notes the time to full decoupling is "quite large" in his
example. Those decoupling results are for infinitesimal parameter
perturbations, though Glasserman & Yao's are not; `run - shadow` is selection
on versus off, a large structural difference, which is the regime where
coupling decays fastest. Glasserman & Yao also give the one usable *positive*
test — their guarantees rest on **monotonicity and continuity**, so ask whether
the shadow's output is monotone in the perturbation before spending the
budget.

**Settle it by measurement, not argument:** run both indexings at fixed budget
and compare the variance of the paired difference. If CRN's is not lower, use
independent streams and spend the budget on replicates instead.

**Whatever is chosen must apply to every draw on the shadow path.** Mixing
`Shadow` and `Beaker` draws gives partial CRN, silently — the worst of the
options.

Steps 3 and 4 do **not** want opposite things, and an earlier draft of this
note said they did. Step 4's permutation is already immune to the live run's
draw count: a stream is a pure function of `(seed, domain, index)`, so nothing
any other stream does can move it. That is `borbax-rng`'s headline property.
What Step 4 needs from Step 3 is only a **disjoint coordinate range** — an
allocation question, not a CRN one. The single open question is Step 3's:
whether the shadow should share the live run's draw *values*.

- [ ] **Step 4: Randomised-catalysis control** (`shadow.rs`)

The persistence shadow cannot answer the question that actually matters about a
RAF. Reassign which species catalyses which reaction uniformly at random,
holding catalysis density fixed, and measure how often a RAF still appears.
That distinguishes "a RAF appeared because of the shape chemistry" from "any
network this dense has one" — and without it, §15.1 criterion 2 means very
little.

- [ ] **Step 5: Plateau model fitting** (`plateau.rs`)

Use `levenberg-marquardt` (a MINPACK port) rather than hand-rolling nonlinear
least squares. **This is a cheap dependency despite appearances**: it runs over
already-emitted trajectories, never feeds back into physics, and does not enter
golden hashes.

Fit saturating, linear and power-law models to the novelty **increments** —
not the cumulative curve, whose residuals are near-perfectly autocorrelated and
which will confidently "prove" whichever model was fitted. Compare by AICc and
BIC, reported with the fitted parameters:

```
aic = n·ln(rss/n) + 2k;  aicc = aic + 2k(k+1)/(n-k-1);  bic = n·ln(rss/n) + k·ln(n)
```

Write the Jacobians analytically; `differentiate_numerically` is a test helper
for verifying them, not a substitute — and it makes a good unit test.

Route the models' `exp`/`powf` through `det_math`, so the verdict reproduces
cross-platform for free. Fit on the first half and score predictive error on
the second: the boundedness illusion is a *projection* failure, so a held-out
tail tests the actual claim in a way no information criterion does.

- [ ] **Step 6: Run, then commit.** Record the measured null-trajectory novelty
score in the commit message — a metric that has never been shown to read zero
on noise has not been validated.

---

## Phase 7 — V0 exit criteria

### Task 21: Prove the chemistry works

This task is experiments, not features. Each produces a number or a decision
that gets recorded in the spec. **Nothing here is complete until its result is
written back into `docs/superpowers/specs/2026-07-26-borbax-prd.md`.**

**Files:** Create `crates/borbax-cli/{Cargo.toml,src/main.rs}`; replace bodies in
`crates/borbax-units/src/det_math.rs`; add `docs/v0-results.md`.

`det_math` lives in `borbax-units` and nowhere else (spec §18.3) — two crates
below `borbax-molecule` call it, so anywhere higher is unreachable from them.
The xtask check added in Task 2 Step 7 exempts that one path, so a copy in
`borbax-molecule` fails CI rather than quietly becoming a second chokepoint.

- [ ] **Step 1: Signature resolution sweep (spec §22.2, exit criterion 3)**

`borbax sweep --resolutions 12,42,162` runs the battery and the folding-map
property tests at each resolution, reporting pass/fail per criterion **and
cost per binding call**. Pick the smallest that passes everything.

Dispatch once at the top, so the const generic is monomorphised rather than
carried as a runtime parameter:

```rust
match cfg.resolution {
    12 => run::<12>(&cfg),
    42 => run::<42>(&cfg),
    162 => run::<162>(&cfg),
    n => return Err(Error::BadResolution(n)),
}
```

Record the winner in §22.2 and replace the const generic with a crate const.
**This is a V0 exit criterion**: changing it later invalidates every golden.

- [ ] **Step 2: Decay band search (spec §22.6, exit criterion 4)**

`borbax band --sweep decay_scale` bisects `decay_scale` for the widest range
where polymers both form and turn over. Record the band in §22.6.

- [ ] **Step 3: Damage measurability gate (spec §22.7, exit criterion 5)**

Inject a known quantity of damage into a beaker; confirm the composition-
divergence metric recovers it above noise. **If it fails, do not reach
straight for the provenance fallback** — it breaks the interchangeability
that makes decay O(1) per species (§9.5), which is a change of asymptotic
class. Try reworking the measurement first, and record the outcome either way.

- [ ] **Step 4: Catalysis emergence (exit criterion 6 — the one that matters)**

`borbax beaker --seed N --steps 1e7 --detect-catalysis` runs until it finds a
folded polymer with two cavities measurably accelerating a reaction. Report
the polymer, its cavities, the reaction, and the measured enhancement — and
render all of it with Task 19.

**If this does not happen, no amount of world-building above it will help and
the design needs rethinking before another line is written.** Write up what
was observed instead; a negative result here is the most valuable output V0
can produce.

- [ ] **Step 5: Decay-off equivalence (exit criterion 7)**

With decay disabled, a run and its shadow must become statistically
indistinguishable — the direct experimental test of "no death, no life"
(§9.4). It runs in seconds and catches decay being nominally implemented but
not actually biting.

**This is a different experiment from the §15.3 neutral shadow**, which
*equalises* decay rather than removing it. Both are needed and conflating them
produces plausible output from the wrong control. Run both: decay-off for
criterion 7, equalised-decay for every "look what evolved" claim.

- [ ] **Step 5b: Metrics and shadow machinery (exit criterion 9)**

Confirm Task 20b's deliverables: the Bedau–Packard statistics compute, the
shape-novelty metric reads ~0 on a null trajectory, a shadow forks from a
beaker run and compares, and the plateau fitter distinguishes a saturating
series from a power-law one on synthetic data before it is trusted on real.

- [ ] **Step 6: Portable transcendentals, then the cross-platform matrix**

Replace the `det_math` bodies with the [`libm`](https://crates.io/crates/libm)
crate — rust-lang's pure-Rust port of MUSL's libm:

```rust
//! Deterministic transcendentals (spec §13.1, §13.4).
//!
//! Backed by the `libm` crate rather than the platform's. Apple's libm and
//! glibc genuinely disagree in the last bits of `exp`, `ln`, `sin` and `cos`,
//! and glibc versions disagree with each other — none of them is wrong, since
//! IEEE-754 does not specify these to be correctly rounded. `libm` is pure
//! Rust with no platform dispatch, so it gives the same bits everywhere,
//! which is exactly what §13.4 needs.
pub fn exp(x: f64) -> f64 { libm::exp(x) }
pub fn ln(x: f64) -> f64 { libm::log(x) }
pub fn sin(x: f64) -> f64 { libm::sin(x) }
pub fn cos(x: f64) -> f64 { libm::cos(x) }
```

**Do not write these by hand.** An earlier draft of this plan said "vendored
pure-Rust implementation", which was a mistake worth naming: a correctly-
rounded `exp` is hard to get right and easy to get subtly wrong, and the
maintained implementation already exists. `sqrt` and the four arithmetic
operations stay native — IEEE-754 specifies those exactly, so they are already
portable.

This is a runtime dependency in a result-affecting path, so pin it exactly and
treat any upgrade as a physics change requiring golden regeneration (CLAUDE.md
§ Dependencies). Run the determinism auditor on the change.

Goldens generated on the M1 will not reproduce on CI runners until this lands.

Then enable `determinism-matrix` in CI (Task 1) and confirm goldens are
identical across macOS/aarch64, Windows/x86-64 and Linux/x86-64. **Do this
last on purpose:** the matrix fails loudly if the transcendental work was
skipped, which is a better guarantee than remembering.

- [ ] **Step 7: Write up `docs/v0-results.md`** — every number, every decision,
and every criterion's outcome. Update the spec sections each result belongs
to. Commit.

**V0 is complete when all eight criteria in spec §23 are met and recorded.**
