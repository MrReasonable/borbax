# rust-performance-expert — V0 plan, post-fix re-review (2026-07-26)

Verdict on the fix pass: **the eight items I was asked about are, with two
exceptions, cheap or free.** The expensive new cost is not in any of them — it
is in Task 13's cavity flood (container choice, not algorithm) and in a decay
path the fix pass half-updated.

---

## P1 · Task 15 Step 3 — `decay_channels` allocates a `Vec` per event, and
## still recomputes a value the fix pass precomputed

**Where:** chemistry plan, Task 15 Step 3, `decay_channels` / `radiogenic_rate`.

**What:** `decay_channels` returns `Vec<(DecayKind, SpeciesId, f64)>` — a heap
allocation of 72 bytes on every call. It is called whenever a species' count
changes, i.e. up to 4× per Gillespie event plus once per channel in the
`FLAT_SCAN_LIMIT` rescan. Worse, it calls `radiogenic_rate()`, which walks
`canon.elem` and does an element-table lookup per atom — while Task 14's fix
pass added `SpeciesRecord::radiogenic_rate` precisely to kill that loop. The
call site was never updated. `radiogenic_rate` also reads `it.record(id).canon`,
a field the fix pass replaced with `kind: SpeciesKind`, so it will not compile
as written.

**Why it costs:** at the §17 target of 10⁵ steps/sec/core, a 4×-per-step malloc
is ~4×10⁵ allocations/sec/core — 50-100 ns each, so 2-4% of the entire budget
burned on an allocation whose result is three tuples. The `radiogenic_rate`
loop is the exact pattern the record's own doc comment warns about: "a
twelve-iteration loop never shows up in a profile and runs millions of times."

**Fix:**

```rust
/// Channels written into a caller-owned buffer. No allocation on the step path.
pub fn decay_channels_into<const D: usize>(
    id: SpeciesId,
    count: f64,
    it: &Interner<D>,
    t: Thermal,
    u: &Universe,
    out: &mut [(DecayKind, SpeciesId, f64); 3],
) {
    let r = it.record(id);
    let thermal_k = borbax_molecule::det_math::exp(-1.0 / (t.get().max(1e-6) * 0.01));
    // Accumulation order pinned: Spontaneous, Solvent, Radiogenic (§13.1).
    out[0] = (DecayKind::Spontaneous, id, count * r.cleave_propensity * thermal_k);
    out[1] = (DecayKind::Solvent,     id, count * r.solvent_rate.max(0.0) * u.consts.decay_scale);
    // Precomputed at intern time (§8.6) — do not re-derive from `kind`.
    out[2] = (DecayKind::Radiogenic,  id, count * r.radiogenic_rate);
}
```

Delete `fn radiogenic_rate` entirely; it has a field now. Hoist `thermal_k` to
the caller if temperature is constant across a rescan (it is, in V0).

**Confidence:** high that the allocation and the duplicate loop are real; high
that the code as written does not compile against Task 14's new record.

---

## P2 · Task 13 Step 3 — the outside flood is B-tree-bound, not cell-bound

**Where:** chemistry plan, Task 13 Step 3, `cavities()`.

**What:** `occupied` is a `BTreeMap<i32, usize>`, `outside`/`candidates`/`seen`
are `BTreeSet<i32>`. Every neighbour test is a B-tree descent with pointer
chasing. The cell count is *not* the problem; the container is.

**Why it costs (magnitude):** a compact 200-mer on the FCC lattice occupies a
region roughly 8 cells on a side; `bounding_box` pads by 2, giving ~12-13 per
side, so `box_cells` enumerates ~2,200 cells (`box_cells` walks *all* integer
cells, including the odd-parity ones no monomer can occupy — half the
enumeration is dead). Per cell the candidate loop does 1 `occupied` lookup + 1
`outside` lookup + 12 neighbour `occupied` lookups = 14 descents. The flood
does ~12 lookups + 1 insert per reached empty cell. Total ≈ 8×10⁴ B-tree
operations, each a handful of cache misses. Estimate 3-8 ms per call — likely
**several times the cost of the 20,000-step anneal that produced the fold**,
which inverts the expected balance between Tasks 11 and 13.

It runs once per species at intern time, so it is not a §8.6 violation. It is
still the largest single item in the intern cost, and interning novel species
is the thing a productive run does constantly.

**Fix:** the data is already dense and the workspace already exists.

```rust
// `occupied` as an O(1) dense grid: FoldWorkspace::occ is already 512 KB of
// u16 keyed by flat cell, already reset by trail (§17). Store monomer+1.
// For the box-local flags, one byte per cell of the padded bounding box —
// ~2.2 KB, L1-resident for the whole pass.
const OUTSIDE: u8 = 1;
const CANDIDATE: u8 = 2;
const SEEN: u8 = 4;
let (nx, ny, nz) = ((hi.0 - lo.0 + 1), (hi.1 - lo.1 + 1), (hi.2 - lo.2 + 1));
let mut flag = vec![0u8; (nx * ny * nz) as usize];   // hoist to a workspace field
let bx = |c: i32| -> usize { let (x,y,z) = fcc::from_flat(c);
    ((z - lo.2) * ny * nx + (y - lo.1) * nx + (x - lo.0)) as usize };
```

`candidates` becomes a `Vec<i32>` pushed in `box_cells` order. **This is
bit-identical**, and the reason is worth writing in the comment: `to_flat` is
`x + y*GRID + z*GRID²` and `box_cells` iterates z-outer/x-inner, so box order
*is* ascending flat-index order — exactly what `BTreeSet` iteration gives. The
step-3 discovery order, which does affect results, is unchanged. `outside` and
`seen` are membership-only and never iterated.

Skip odd-parity cells in `box_cells` (no FCC site can be there) to halve the
enumeration.

**Confidence:** high that B-tree lookups dominate; high that the replacement is
bit-identical; the 3-8 ms figure is reasoned from first principles and wants
benchmark 3 below.

---

## P3 · Task 14 — `SpeciesRecord` is 832 bytes and the step loop reads four of them

**Where:** chemistry plan, Task 14 Step 3, `SpeciesRecord<D>` / `Interner::records`.

**What (the size you asked for), at D = 42:**

| Field | Bytes |
|---|---|
| `kind: SpeciesKind` (`CanonForm` is 64 B at align 16; `Polymer` is a 24-B `Vec`) | 80 |
| `sig: Signature<42>` (2 × 42 × f64) | 672 |
| `mass`, `bond_count`(+pad), `radiogenic_rate`, `cleave_propensity`, `solvent_rate` | 40 |
| `fold: Option<FoldId>` | 8 |
| `cavities: Vec<Cavity<42>>` | 24 |
| **total** | **832** |

At D = 162 it is 2,752 bytes. `Cavity<42>` is 704 bytes each and they hang off
a `Vec`, so a 6-cavity polymer is a 4.2 KB out-of-line block plus a `Vec<i32>`
per cavity — 7 allocations per polymer species.

**Why it costs:** the propensity rescan (Task 18 Step 5, `FLAT_SCAN_LIMIT` =
256) reads `cleave_propensity`, `solvent_rate`, `radiogenic_rate` and `mass` —
32 bytes — from up to 256 records per Gillespie step. Those 32 bytes sit 13
cache lines apart. Touched footprint: 256 × 832 B = 208 KB, over M1's 128 KB
L1D, so every rescan is an L2-latency walk. With the four hot fields in a
parallel array it is 256 × 32 B = 8 KB, comfortably L1-resident, and the
prefetcher gets a unit stride. At 10⁵ steps/sec that is 2.6×10⁷ record touches
per second per core, so the difference is not noise.

`sig` is cold for the step loop (only the affinity-memo *miss* path reads it)
but hot for channel creation, which does an O(species²) affinity pass — and
that pass currently strides 832 bytes to read 672 contiguous ones.

**Fix:** three parallel `Vec`s behind the existing `record()` API shape.

```rust
/// Read on every propensity update. One cache line per species (§17).
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct SpeciesHot {
    pub cleave_propensity: f64,
    pub solvent_rate: f64,
    pub radiogenic_rate: f64,
    pub mass: Mass,          // i64
}   // 32 bytes, two species per cache line

pub struct Interner<const D: usize> {
    ids: BTreeMap<CanonForm, SpeciesId>,
    hot:  Vec<SpeciesHot>,        // step loop
    sigs: Vec<Signature<D>>,      // channel creation / affinity passes
    cold: Vec<SpeciesCold<D>>,    // kind, bond_count, fold, cavities
}

impl<const D: usize> Interner<D> {
    #[inline] pub fn hot(&self, id: SpeciesId) -> &SpeciesHot { &self.hot[id.0 as usize] }
    #[inline] pub fn sig(&self, id: SpeciesId) -> &Signature<D> { &self.sigs[id.0 as usize] }
    #[inline] pub fn cold(&self, id: SpeciesId) -> &SpeciesCold<D> { &self.cold[id.0 as usize] }
}
```

Do it now rather than later: the retrofit touches every call site, and there
are currently none. It does not make the record a pointer-chasing structure —
the record was never chased; the cost is stride, not indirection.

**Confidence:** medium-high. The mechanism is certain; whether it shows up
depends on how many channels the rescan touches in a real run. This is settled
by benchmark 2, which the plan already schedules for a different reason.

---

## P4 · Task 16 — `catalysis_factor` computes O(C²) affinities where O(C) suffice

**Where:** chemistry plan, Task 16 Step 3.

**What:** `affinity_ordered(&cb.signature, sub_b, ...)` sits in the inner loop
but depends only on `j`. It is recomputed C times for each `j`.

**Why it costs:** `affinity_ordered` at D = 42 is ~2,520 inner iterations,
serially dependent on two accumulators — call it 2-3 µs (benchmark 1). With
C = 6 cavities that is 36 calls where 12 would do: ~60 µs wasted per
(polymer, substrate-pair) channel creation, and channel creation is O(species²)
in the worst case.

**Fix:**

```rust
// Per-cavity fits, computed once. Order is cavity order, which is fixed by
// Task 13's discovery order (§13.1).
let mut fits: Vec<(f64, f64)> = Vec::with_capacity(cavities.len());
for c in cavities {
    fits.push((
        affinity_ordered(&c.signature, sub_a, g, k),
        affinity_ordered(&c.signature, sub_b, g, k),
    ));
}
for (i, ca) in cavities.iter().enumerate() {
    if fits[i].0 < bind_threshold { continue; }
    for (j, cb) in cavities.iter().enumerate() {
        if i == j || fits[j].1 < bind_threshold { continue; }
        /* ... unchanged arithmetic ... */
    }
}
```

Bit-identical: same function, same arguments, evaluated once instead of C times.

**Also, not my lane but adjacent:** the module doc says `BIND_THRESHOLD`,
`IDEAL_SEPARATION`, `MAX_ENHANCEMENT` and the `* 1.0e3` were all replaced by
generated constants; the code below it still uses all four literals. The prose
was fixed and the code was not. `* 1.0e3` decides exit criterion 6. Also the
`//!` block at line 2169 sits after `use` statements and will not compile.

**Confidence:** high.

---

## P5 · Task 20b Step 1 — activity counters must be event-incremental, not a per-step sweep

**Where:** chemistry plan, Task 20b Step 1, `ActivityStats::per_species`.

**What:** `per_species: Vec<f64>` "updated over time" reads as a sweep over all
S species per step. It must not be.

**Why it costs:** at S = 10⁴ species and 10⁵ steps/sec/core that is 10⁹
updates/sec — three orders past the machine, and it would make Task 20b, a
*metrics* task, the dominant cost in the engine. It is also wrong on its own
terms: Bedau activity is the time-integral of existence, and Gillespie steps
have variable Δt, so a per-step increment of 1.0 measures event count, not
persistence.

**Fix:** integrate exactly, O(1) per event. A Gillespie firing changes at most
four species' counts.

```rust
/// a_i(t) = ∫ count_i dt, maintained exactly at O(1) per event.
/// Each species accumulates in its own event order, so the sum is
/// reproducible without pinning a global order (§13.1).
pub struct Activity {
    a: Vec<f64>,          // cumulative, indexed by SpeciesId
    last_t: Vec<f64>,     // WorldYear of that species' last count change
}

impl Activity {
    /// Call immediately *before* changing `count_i`.
    #[inline]
    pub fn touch(&mut self, i: SpeciesId, count_before: f64, now: f64) {
        let k = i.0 as usize;
        self.a[k] += count_before * (now - self.last_t[k]);
        self.last_t[k] = now;
    }
}
```

`diversity`, `cumulative`, `mean_cumulative` and `new_activity` are O(S)
reductions over `a` in `SpeciesId` order — run them on the reporting interval,
not the step loop. The shadow fork doubles step cost; §22.3 already says it is
spawned on demand, so keep it that way. The `levenberg-marquardt` fit is
offline over emitted trajectories and costs nothing in the loop.

**Confidence:** high on the asymptotics; the exact-integration form is standard
and I would want the alife-researcher to confirm it matches the Bedau
definition they intend.

---

## P6 · Task 6 Step 3 — `search()` allocates inside a loop over colour values

**Where:** main plan, Task 6 Step 3, `fn search`.

**What:** `let members: Vec<usize> = (0..n).filter(...).collect();` is inside
`for c in 0..=u8::MAX`, and `let mut order = vec![0u8; n];` runs per leaf.

**Why it costs:** up to 256 heap allocations per search node, purely to learn
whether a colour class has more than one member, in a search bounded at 50,000
leaves. Nothing about it is needed — the answer is a popcount.

**Fix:**

```rust
// Colours are dense in 0..n after the densification fix, so a u16 mask per
// colour is exact and the whole scan is register-resident.
let mut counts = [0u8; MAX_ATOMS];
for v in 0..n { counts[colour[v] as usize] += 1; }
let Some(target) = (0..n).find(|&c| counts[c] > 1) else { return };
// ... and at the leaf:
let mut order = [0u8; MAX_ATOMS];
for v in 0..n { order[colour[v] as usize] = v as u8; }
let form = encode(m, &order[..n]);
```

Note the loop bound also drops from 256 to `n`: after the densification fix
colours cannot exceed `n`, so `0..=u8::MAX` is scanning 244 impossible values.

**Confidence:** high.

---

## P7 · Task 14 — `AffinityMemo` uses `Vec<Option<f64>>`, doubling the one table §17 names

**Where:** chemistry plan, Task 14 Step 3, `AffinityMemo`.

**What:** `f64` has no niche, so `Option<f64>` is 16 bytes. The table is twice
the size it needs to be and holds half as many entries per cache line.

**Why it costs:** §17 calls the binding memo "what makes 3D binding
affordable". At k = 4,096 the triangular table is 4096·4097/2 = 8.4M slots:
134 MB with `Option`, 67 MB with a sentinel. It is also allocated and zeroed
eagerly in `new()`. Affinity is a negated sum of squares — always finite and
≤ 0 — so `NAN` is an unreachable value and a sound sentinel.

**Fix:**

```rust
pub struct AffinityMemo { k: u32, vals: Vec<f64> }   // NAN = not yet computed

impl AffinityMemo {
    #[must_use]
    pub fn new(k: u32) -> Self {
        Self { k, vals: vec![f64::NAN; (k as usize) * (k as usize + 1) / 2] }
    }
    // in get(): `let v = self.vals[s]; if !v.is_nan() { return v; }`
}
```

Determinism-neutral: memoising a pure function, exactly as the existing
`the_affinity_memo_agrees_with_direct_computation` test asserts. Keep that test.

**Confidence:** high on the size; whether 67 MB vs 134 MB matters depends on
the chosen `k`, which V0 has not fixed yet.

---

## P8 · Task 13 — `lining_signature`: hoist the invariant, pass `pos_of` in

**Where:** chemistry plan, Task 13 Step 3, `lining_signature`.

**What:** `reach_of` recomputes `rel` — three subtractions and three
multiplications by `LATTICE_SPAN` — inside both passes of the D loop, so 2 × 42
= 84 times per lining monomer. And `pos_of` is a `vec![(0,0,0); p.len()]`
(2.4 KB at `MAX_POLYMER`) allocated per *cavity*, with a full walk of the
200-entry `occupied` map behind it.

**Fix:** build `pos_of` once in `cavities()` and pass `&[(i32,i32,i32)]` in;
build the scaled relative vectors once, outside the direction loop:

```rust
// Loop-invariant: same operands, same operations, so bit-identical.
let rel: Vec<[f64; 3]> = lining.iter().map(|&m| {
    let (mx, my, mz) = pos_of[m];
    [(f64::from(mx) - centre[0]) * LATTICE_SPAN,
     (f64::from(my) - centre[1]) * LATTICE_SPAN,
     (f64::from(mz) - centre[2]) * LATTICE_SPAN]
}).collect();
for dir_idx in 0..D { /* pass 1 and 2 index `rel`, dot product only */ }
```

**Net answer to the question asked:** two passes over ~40 lining monomers × 42
directions is ~3,400 iterations — microseconds, and *better* than the previous
linear scan per monomer, which was O(D · lining · |occupied|). The two-pass
structure is semantically required (pass 2's weights depend on pass 1's
`nearest`) and must stay. Keep it; just hoist the invariant. This is a ~2× on a
function that is not the bottleneck, so it is the lowest-priority item here.

**Confidence:** high, low importance.

---

# The things you asked about that are fine

**1 · The `anti[perm[i]]` double gather (Task 10).** Free, or close enough that
I would not change it for speed. The loop is bound by its two serial FP-add
chains — `shape` and `charge` each take one dependent add per iteration, ~3
cycles' latency on M1 — so ~3-4 cycles/iteration, ~7,500 cycles per
`affinity` call at D = 42. The added gather is an L1 hit into a 42-byte table
that will never leave L1, it is independent across `i`, and M1 issues three
loads per cycle against ~1.3 loads/cycle of demand. The out-of-order window
covers the perm→anti→b.r chain entirely.

A composed `anti_perm[r][i] = anti[perms[r][i]]` is 2,520 bytes at D = 42 and
9,720 at D = 162 — no locality problem either way, since M1's L1D is 128 KB
and `perms` is already the same size. So it does not *hurt*; I just do not
expect it to help. If you build it, build it for the correctness reason rather
than the speed one: `anti_perm[r][i]` cannot be "simplified" back into the
reflection bug by deleting one index, whereas `g.anti[perm[i] as usize]` is one
deletion away from S1 all over again. That is the geometry reviewer's call, not
mine. Either way it is bit-identical — pure index composition, and
`antipode_commutes_with_every_rotation` already pins the invariant.

The rejected alternative, for the record: folding `anti` into the signature
(`r_anti[i] = r[anti[i]]`, per species at canonicalisation time) also removes
the gather, but doubles `Signature<D>` to 1,344 bytes and so makes P3 worse.
Do not do that one.

**2 · `affinity` delegating to `affinity_with_rotation`.** Free. Both are
generic over `D`, so they monomorphise into the *calling* crate and are inline
candidates cross-crate with no `#[inline]` needed; a body that is one call plus
`.0` is inlined at any optimisation level that runs the inliner, including the
`opt-level = 1` dev profile. `(f64, usize)` returns in registers under AAPCS64,
so even un-inlined the discarded half is one dead register. Inside the loop,
`best = (score, r)` adds one conditional register move per rotation — 60 out of
~7,500 cycles. No action.

**5 · Task 6's densified colouring.** Fine. n ≤ 12 so the rank loop is ≤ 1,728
`u8` comparisons with an early-exit inner scan — well under a microsecond, and
it runs **once per `canonicalise`, outside the search**. The sort-based
alternative would allocate. Leave it.

Worth knowing where the real cost in that file is: `refine()` uses the *same*
O(n³) rank pattern but over `[u8; 37]` signature rows, and it runs once per
*search node*. That is ~1,728 array comparisons of 37 bytes per refinement
round, per node, in a search capped at 50,000 leaves. Whether that matters
depends entirely on the leaf distribution over real molecules, which
`SearchStats.leaves` already exists to measure — see benchmark 4. Do not
pre-optimise it; measure it.

**6 · Task 8's `canonicalise_frame`.** Fine, and by a wide margin. ~20 flops to
find the frame plus 9 multiply-adds per atom — under 300 flops for n = 12,
against SMACOF's 240 iterations × n² × ~30 flops ≈ 10⁶. It is 0.03% of `embed`,
which itself runs once per species. No action, and no need to revisit.

**7 · Has the record become pointer-chasing?** No — see P3. `records` is a
`Vec<SpeciesRecord>` and lookups are `id.0 as usize`, one indexed load. The
`Vec<u8>` inside `Polymer` and the `Vec<Cavity>` are genuine out-of-line
allocations, but nothing on the step path dereferences them (with one
exception, which is P1's `radiogenic_rate` reading `canon.elem` — delete that
and the step path touches no pointer in the record). The problem is stride, not
indirection.

---

# Previous recommendations — survival check

| Recommendation | Status |
|---|---|
| **Affinity memo** | **Survived.** Task 14 has `AffinityMemo`, dense triangular, keyed on `SpeciesId`, with a test proving it agrees with direct computation. `affinity_ordered` gives it a canonical key. Only defect is the `Option<f64>` slot width (P7). |
| **Segment tree over Fenwick** | **Survived, and strengthened.** Parents recomputed as `left + right`; `totals_stay_bit_identical_to_a_rebuild` tests it over 50,000 updates; freed slots set to 0.0 rather than compacted, with a test. The `FLAT_SCAN_LIMIT` crossover survived too, and Task 18 Step 6 requires the measured value in the commit message. Best-defended structure in the plan. |
| **Fold workspace reset by trail** | **Survived intact.** `FoldWorkspace { occ: Box<[u16]>, trail: Vec<i32> }`, 64³ × 2 B = 512 KB, matching §17 exactly; `clear()` walks the trail, O(chain) not O(grid), with the reason in the comment. |
| **Flat arena for conformations** | **Not adopted.** `FoldCache` stores `Entry { fold: Fold, .. }` where `Fold` owns `coords: Vec<i32>` and `exposure: Vec<u8>` — two allocations per cached fold. At the §17 budget of 2 GB and ~1 KB per entry that is ~2M entries and ~4M live small allocations. Two consequences: allocator pressure and fragmentation on the miss path, and a **budget accounting gap** — `bytes` counts `size_of::<Fold>()` plus the two slices, but not the `HashMap` entry (~32 B), the `slots` tuple (~48 B), or per-allocation allocator overhead. At 2M entries the unaccounted overhead is on the order of 200 MB, so the cache will overshoot its stated budget by ~10%. Not urgent — but if the fold cache ever misses its >95% hit-rate target and gets resized, fix the accounting first, and consider a single `Vec<i32>` arena with `(offset, len)` per entry, which also makes the whole cache one allocation to checkpoint into a keyframe. |

---

# What to benchmark rather than reason about (updated)

Criterion, per `CLAUDE.md`'s dependency doctrine — a timing assertion in a
`#[test]` is a flake that gets deleted.

1. **`affinity` at D = 12 / 42 / 162**, with `black_box` on both signatures,
   real embedded molecules rather than synthetic patterns. Report cycles/call
   and variance. Settles the composed-table question, prices P4, and gives
   §22.2's resolution sweep its cost curve. *Second variant:* the same kernel
   with a fixed-shape pairwise summation (e.g. 42 = 6 blocks of 7, blocks
   combined in index order) — deterministic because the tree is fixed at
   compile time, and *more* accurate than the serial sum, but it moves every
   golden. **Decide this before Task 10 lands, while there are no goldens to
   move.** I am flagging the trade rather than recommending it; it is
   determinism-auditor's call whether the accuracy gain justifies the change,
   and the precedence order in `CLAUDE.md` (2 above 3) is the relevant clause.

2. **Flat-scan vs tree-descent crossover** for propensity updates, swept over
   channel count. Already scheduled at Task 18 Step 6. Run it twice — once with
   the monolithic `SpeciesRecord`, once with P3's split — and it settles both
   questions for the price of one harness.

3. **`cavities()` vs `fold()` on the same chain**, at n = 50 / 100 / 200,
   B-tree version against P2's dense-grid version. This is the one I most want
   measured, because my estimate says cavity detection costs more than the
   anneal that feeds it, and that would be surprising enough to be worth
   knowing before Task 16 depends on it.

4. **`canonicalise` leaf distribution** over the Task 6 random-molecule
   generator: report the `SearchStats.leaves` histogram and p99, not just the
   mean. `refine`'s O(n³) rank loop is only worth touching if the tail is fat.
   This also tells you how close the `SEARCH_LEAF_CAP` guard is to firing.

5. **Steady-state `Beaker::step` throughput** against §17's 10⁵ steps/sec/core,
   swept over species count. The single number that says whether the budget is
   met. Measure steady state, not the first thousand steps — early steps have a
   tiny channel set and an empty affinity memo, and will flatter it by an order
   of magnitude.

6. **Fold-cache hit rate at the 2 GB budget** over a Task 20 forty-universe
   battery, against §17's >95% target — plus measured RSS against the cache's
   own `stats.bytes`, which is the accounting gap above.

7. **`det_math::exp` against a tabulated Arrhenius lookup.** Task 14's `rate()`
   already names itself a tabulation candidate. `libm`'s `exp` is ~20-40 ns;
   the table is a hash or an index plus an interpolation, and if the argument
   set is not actually small the table is a pessimisation with extra
   determinism surface. Measure the *distinct argument count* in a real run
   first — that, not the timing, decides it.

Dropped from my previous list: the `may_bind` prefilter benchmark. It now has a
provable Cauchy-Schwarz bound and the plan tests the bound as a property over
100,000 pairs, so its *correctness* is settled; its rejection *rate* comes out
of benchmark 1 for free.
