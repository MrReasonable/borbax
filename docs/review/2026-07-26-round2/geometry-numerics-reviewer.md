# Geometry & numerics re-review — Tasks 6–13 (post-fix)

All figures below come from a Python reproduction of the plan's own code
(`smacof`, `canonicalise_frame`, `signature`, `canonicalise`, `affinity`,
`may_bind`, the geodesic/rotation/anti construction, the FCC start walk).
Scripts in the session scratchpad. Where I did not measure, I say so.

---

## G1 · `ideal_gap` is on the wrong length scale — binding is 93% size, 7% shape

**Where:** Task 10 Step 3 `affinity_with_rotation` (`v0.md:4234`), fed by
`UniverseConsts.ideal_gap = rng.next_f64_range(1.6, 2.8)` (`v0.md:2046`).
**Silently wrong.**

`shape_term = -(r_A[i] + r_B[j] - ideal_gap)²`. Per §8.3, `ideal_gap` is the
*sum of two extents at contact* — a centre-to-centre separation. But `r` is a
support function: `pos·d + radius`, so it grows with molecule size. Generated
radii reach ~4.5 (`v0.md:1715`) and an 8-atom embedding has R_g of several
Span, so real extents are 1.5–6 and `r_A + r_B` is 3–12 against a gap of ~2.
The squared term is then dominated by a constant size offset and the
complementarity signal is a small perturbation on it.

Measured over 30 embedded 8-atom molecules, 465 pairs, D=42:

| gap | rotation contrast `(max−min)/|mean|` over the 60 rotations |
|---|---|
| fixed `ideal_gap = 2.0` | **0.069** |
| pair-scaled `mean(r_A)+mean(r_B)` | **1.066** |

The 2,520-iteration rotation search — the mechanism the entire design rests
on — moves the score by 7% of its magnitude. 93% of `affinity(A,B)` is
"how big are these two molecules".

**How it surfaces:** not as a failure. Binding happens, scores vary, and the
ordering is monotone in *size*: small molecules bind everything, large ones
bind nothing. Homochirality, cavity selectivity and neutral networks all
degrade to noise on top of a size ranking. §23 criterion 6 reads as "catalysis
did not emerge" and the decay band gets blamed.

**Fix** (a deliberate physics change — regenerate goldens):

```rust
/// Ideal centre-to-centre separation for *this pair*.
///
/// §8.3's IDEAL_GAP is a separation, so it has to scale with the two bodies.
/// A per-universe constant on the order of one atomic radius makes the shape
/// term a size-difference term: measured, it leaves the 60-rotation search
/// contributing 6.9% of the score instead of 107%.
#[inline]
fn pair_gap<const D: usize>(a: &Signature<D>, b: &Signature<D>, k: &BindConsts) -> f64 {
    let ma: f64 = a.r.iter().sum::<f64>() / D as f64;
    let mb: f64 = b.r.iter().sum::<f64>() / D as f64;
    k.gap_factor * (ma + mb)   // gap_factor ~ 1.0, generated per universe
}
```
`ma`/`mb` are rotation-invariant, so symmetry and the `may_bind` bound both
survive unchanged.

**Test that would have caught it** (and is the right acceptance test):
```rust
#[test]
fn the_rotation_search_decides_the_score() {
    // If the best and worst orientations of the same pair score nearly the
    // same, `affinity` is not measuring complementarity (§8.3).
    for (a, b) in random_species_pairs(200) {
        let s: Vec<f64> = (0..N_ROTATIONS).map(|r| score_at(&a, &b, r, &g, &k)).collect();
        let (hi, lo) = (s.iter().copied().fold(f64::MIN, f64::max),
                        s.iter().copied().fold(f64::MAX, f64::min));
        let mean = s.iter().sum::<f64>() / N_ROTATIONS as f64;
        assert!((hi - lo) / mean.abs() > 0.25, "orientation barely matters: {}", (hi-lo)/mean.abs());
    }
}
```

**Confidence:** measured, high. The exact numbers depend on the generated
radius distribution; the *structural* mismatch (a size-independent constant
compared against a size-dependent support function) does not.

---

## G2 · `canonicalise_frame` does not fix S7, and is worse than doing nothing

**Where:** Task 8 Step 3 `canonicalise_frame` (`v0.md:3604`). **Silently wrong
in effect; the plan's own test will fail, but the remediation note points at
the wrong knob.**

Reproducing the plan's `embed` exactly (240 SMACOF iterations, `START`,
centring, then `canonicalise_frame`) over 200 random 8-atom trees with a
one-atom addition — the plan's own experiment:

| frame | unaligned median | p90 | max | fraction over the plan's 0.6 bound |
|---|---|---|---|---|
| audit's pre-fix figure | 0.74 | 1.34 | 2.59 | 67.5% |
| **plan's `canonicalise_frame`** | **0.87** | **1.79** | **3.11** | **69.5%** |
| **no frame canonicalisation at all** | **0.64** | **1.27** | 2.48 | **54.5%** |
| gyration-tensor axes, 3rd-moment signs | 1.11 | 2.96 | 3.54 | 68.5% |
| max-magnitude atom pair | 1.43 | 2.90 | 3.53 | 71.0% |
| atom 0 + max-perpendicular atom | 0.93 | 2.48 | 3.39 | 69.0% |
| weighted moments (Lipschitz in positions) | 1.06 | 2.18 | 2.94 | 83.0% |

Aligned RMSD/R_g is 0.12 median in every row — SMACOF is doing its job. The
Kabsch rotation between the two embeddings is median 21.5°, p90 53°, max 176°.

Three separate defects, and the third is why none of my alternatives helped:

1. **Conditioning.** `e0` is `pos[a0]` normalised, `a0` = *first* atom past
   `1e-9 * scale_by`. The prompt's worry (two atoms nearly equidistant flipping
   `a0`) does **not** apply — selection is by index, so `a0` is almost always 0.
   The real hazard is the opposite: `e1` is the first atom whose perpendicular
   component clears the same threshold, and that component is measured at
   **0.0013·R_g** in the worst case. The roll about `e0` is then set by 0.13% of
   the molecule and swings freely under perturbation.
2. **`scale_by = 1.0f64.max(norm(pos[0]))` is not a scale.** It is
   `max(1, |pos[0]|)`, so for a molecule smaller than 1 Span the guard is a bare
   absolute `1e-9`, and `pos[0]` is exactly the atom that may sit at the
   centroid. Every magic epsilon should name its scale; this one names the wrong
   quantity. It should be `radius_of_gyration()`.
3. **Canonical atom order is not continuous under graph edits.** Adding one atom
   can permute the whole canonical labelling, so *any* frame anchored to "the
   first atom in canonical order" inherits a discontinuity that no better
   selection rule removes. That is why max-magnitude, max-perpendicular,
   weighted-moment and principal-axis variants all measure the same or worse.
   The plan's own test does not even exercise this — it grows `m` without
   re-canonicalising.

**Handedness is fine:** `e2 = e0 × e1` and the map `v ↦ (v·e0, v·e1, v·e2)` has
det +1 (`e1 × e2 = e0`), so it is a proper rotation and §22.8's chirality
survives. Verified on paper and numerically.

**How it surfaces — measured end to end.** Building the full pipeline
(embed → frame → signature → `canonicalise` over 60 rotations) and comparing a
one-atom mutant against an unrelated random tree of the same size, 120 trials:

```
near/far signature distance ratio: median 0.90, p90 1.31
in 32% of cases the mutant is FURTHER from its parent than an unrelated molecule is
```

That is the concrete damage: **signature space has no locality.** Mutation is a
random walk, neutral networks cannot form, and §8.4's evolvability argument has
nothing under it. Task 9's `small_changes_move_the_signature_less_than_unrelated_molecules`
uses one example (chain(8)→chain(9) vs a star), which passes; the property holds
68% of the time over random trees.

**Fix:** I do not have a validated one, and I would rather say that than invent
one. What I can state:
- Do **not** relax the 0.6 threshold, and do **not** raise `ITERATIONS` — aligned
  RMSD is already 0.12, so iterations are not the constraint. The plan's Step 4
  note ("Raise `ITERATIONS` or revisit the weighting first") is misleading.
- Removing `canonicalise_frame` improves the metric. If it stays, it needs to
  earn its place against that baseline.
- The direction worth measuring next, in order of cheapness: (a) restate the
  requirement in terms of the *canonicalised signature distance* rather than
  atom displacement, and see whether an SO(3)-invariant descriptor (the
  `sorted_extents` multiset, or a low-order spherical-harmonic power spectrum)
  gives locality without any frame at all; (b) if it does, use the invariant
  descriptor for species identity and novelty, and leave `affinity`'s own
  rotation search to handle orientation — it already does.

**Test to add** (replaces the example-based one, and is the specification):
```rust
#[test]
fn a_one_atom_mutation_stays_near_its_parent_in_signature_space() {
    // Over random trees, not one chain. The parent/mutant distance must beat
    // the parent/unrelated distance far more often than chance, or mutation is
    // a random walk and §8.4's neutral networks cannot exist.
    let mut wins = 0;
    for _ in 0..200 { /* build m, m+atom, unrelated; canonicalise both graphs */
        if d(sig(&m), sig(&grown)) < d(sig(&m), sig(&unrelated)) { wins += 1; }
    }
    assert!(wins > 180, "only {wins}/200 mutants were closer than an unrelated molecule");
}
```
Note it must call `canonicalise()` on the *graph* before embedding, which the
current tests do not.

**Confidence:** measured in a faithful reproduction, high for the direction and
rough magnitude. My per-atom `affinity` values are a stand-in, but `r` dominates
the signature distance and `r` is modelled exactly.

---

## G3 · The boustrophedon fold start retraces itself from monomer 12

**Where:** Task 11 Step 4, the start walk (`chemistry.md:455–468`).
**Obviously broken, but it will not fail where the plan looks.**

`fcc::NEIGHBOURS` is laid out in **antipodal pairs**: `NEIGHBOURS[2] ==
-NEIGHBOURS[1]` and `NEIGHBOURS[3] == -NEIGHBOURS[0]`, and so on for
(4,7),(5,6),(8,11),(9,10). Cycling `dir` in index order therefore walks a
direction and then, one run later, its exact reverse. Traced:

```
i=10 (42,32,32)  dir 1
i=11 (43,31,32)  dir 1
i=12 (42,32,32)  dir 2   <-- revisits monomer 10
i=13 (41,33,32)  dir 2   <-- revisits monomer 9
...
i=22 (32,32,32)  dir 3   <-- revisits monomer 0
```

Duplicate-cell counts: n=45 → 21 collisions; n=80 → 44; n=200 → **164
collisions, 36 distinct cells**. The whole chain occupies a 13×13×7 box at every
length.

`occupy` overwrites `occ[cell]` without checking, so the collision is silent at
construction. It then corrupts the anneal: when one of two co-located monomers
moves, `self.occ[old] = 0` wipes the *other* monomer's occupancy, so
self-avoidance stops holding globally.

`the_chain_stays_connected_and_self_avoiding` (n=45) should fail — good — but
`no_chain_length_wraps_the_grid`, which is the test that goes to
`MAX_POLYMER`, only checks coordinate bounds, and Task 13 folds 60/70/80-mers
through a path with no self-avoidance assertion at all.

**Fix:**
```rust
// Greedy compact start. Cycling `dir` through NEIGHBOURS in index order
// retraces the walk, because NEIGHBOURS is stored in antipodal pairs
// ((0,3), (1,2), (4,7), (5,6), (8,11), (9,10)) — a run in direction k is
// followed one run later by a run in -k. Pick, deterministically, the free
// in-bounds neighbour that stays closest to the centre.
let mut c = centre;
for i in 0..n {
    debug_assert!(self.is_free(c), "fold start collided at monomer {i}");
    coords.push(c);
    self.occupy(c, i);
    if i + 1 == n { break; }
    let next = fcc::NEIGHBOURS.iter()
        .map(|&d| c + d)
        .filter(|&t| fcc::in_bounds(t) && self.is_free(t))
        .min_by_key(|&t| {                       // integer key: no float tie-break
            let (x, y, z) = fcc::from_flat(t);
            let (cx, cy, cz) = (fcc::GRID / 2, fcc::GRID / 2, fcc::GRID / 2);
            ((x-cx).pow(2) + (y-cy).pow(2) + (z-cz).pow(2), t)  // t breaks ties
        });
    match next { Some(t) => c = t, None => break }   // trapped: chain ends here
}
```
and `fn occupy(&mut self, cell: i32, m: usize) { debug_assert!(self.occ[cell as usize] == 0); .. }`.

**Test:** move the self-avoidance assertion into `no_chain_length_wraps_the_grid`
so it runs at 45/64/80/120/`MAX_POLYMER`, and assert it on the *initial*
conformation as well as the annealed one.

**Confidence:** verified by direct simulation of the plan's arithmetic.

---

## G4 · Cavity `lining_signature`: `nearest` is a minimum over a hemisphere, so it is ≈0 everywhere

**Where:** Task 13 Step 3, `lining_signature` (`chemistry.md:1361–1368`).
**Silently wrong.**

`nearest` = min over lining monomers of `(p_m − centre)·dir`, admitted whenever
`reach > 0.0`. "Reach > 0" is the whole `+dir` **hemisphere**, not a cone around
`dir`. A cavity is lined by a shell of monomers; their projections onto any
direction spread continuously from −R to +R, so the smallest *positive*
projection is set by whichever lining monomer happens to sit nearest the
equatorial plane — near zero, for every direction.

So `sig.r[d] ≈ 0` and nearly direction-independent: the cavity signature carries
almost no shape. The previous review's `max → min` change swapped one wrong
quantity (the far wall) for another (the equator).

**How it surfaces:** `affinity(cavity, molecule)` becomes indiscriminate — every
cavity looks the same, so every cavity "recognises" the same substrates.
Criterion 6 reports catalysis for chains with any two cavities.

**Fix** — the wall distance along the ray, not the hemisphere:
```rust
// Only monomers within one lattice unit of the ray are "that way". Without
// this the min is over the whole +dir hemisphere and is set by whichever
// lining monomer sits nearest the equatorial plane — i.e. ~0 for every
// direction, for any roughly convex cavity.
const RAY_RADIUS: f64 = 1.0 * LATTICE_UNIT;   // one lattice site, the probe width
let mut nearest = f64::INFINITY;
for &m in &lining {
    let rel = rel_of(m);
    let along = rel[0]*dir[0] + rel[1]*dir[1] + rel[2]*dir[2];
    if along <= 0.0 { continue; }
    let perp2 = rel[0]*rel[0] + rel[1]*rel[1] + rel[2]*rel[2] - along*along;
    if perp2 <= RAY_RADIUS * RAY_RADIUS && along < nearest { nearest = along; }
}
```
Cheaper and arguably better: take the support function of the **cavity cells**
themselves, `max over cells of (cell − centre)·dir`, which is the exact analogue
of §8.2's molecule support function and needs no cone test.

**Test:** build a synthetic ellipsoidal cavity (occupied shell, empty interior)
with semi-axes 4/2/2 and assert `sig.r` is ~2× larger along the long axis than
the short ones. The current code gives ~0 for both.

**Confidence:** derived, high. The test above settles it in one run.

---

## G5 · Cavity–substrate fit is inverted: a bigger cavity prefers a smaller substrate

**Where:** Task 13 `lining_signature` output fed to `affinity` (§8.3, Task 16).
**Silently wrong.**

The kernel scores `(r_cav[i] + r_sub[anti[perm[i]]] − ideal_gap)²`. With
`r_cav[i] = h(d_i)` (distance from centre to wall), the squared term is minimised
when `r_sub = ideal_gap − h`. **`d(preferred r_sub)/dh = −1`**: a roomier cavity
selects a *smaller* substrate. Containment is a difference (`h − r_sub =
clearance`); §8.3's form is a sum, because there the two extents are measured
from two *different* centres. The index is also wrong: the kernel pairs slot `i`
with the substrate's antipodal slot, but for a body sitting *inside* a cavity
both are measured from the same centre and want the same direction.

**Fix** — store the cavity so that §8.3's kernel, unchanged, computes the
clearance:
```rust
// The kernel pairs our slot i with the substrate's slot anti[perm[i]], and
// sums the two extents. For containment both quantities share a centre and
// the fit condition is a *difference*, so we store gap − h at the antipodal
// slot; the kernel then evaluates exactly (r_sub(-d) − h(-d)), the clearance.
// Storing h directly makes a roomier cavity prefer a *smaller* substrate.
let j = g.anti[dir_idx] as usize;
sig.r[j] = k.ideal_gap - nearest;   // may be negative; that is correct here
sig.a[j] = lining_character;
```

**Test:**
```rust
#[test]
fn a_bigger_cavity_prefers_a_bigger_substrate() {
    let small = sphere_signature(1.0);
    let large = sphere_signature(2.0);
    let c1 = spherical_cavity_signature(1.5);
    let c2 = spherical_cavity_signature(2.5);
    assert!(affinity(&c1, &small, &g, &k) > affinity(&c1, &large, &g, &k));
    assert!(affinity(&c2, &large, &g, &k) > affinity(&c2, &small, &g, &k));
}
```

**Confidence:** derived on paper, high on the monotonicity argument (it is one
line of algebra); medium on my proposed encoding being the best repair — the
test above is what settles it.

---

## G6 · The enclosure gate now fragments the cavities the flood found, and the flood only finds sealed voids

**Where:** Task 13 Step 3, Steps 1–3 (`chemistry.md:1224–1295`).
**Silently wrong, in the "there is no chemistry" direction.**

Two problems that compound:

1. **`MIN_ENCLOSURE` is now redundant and harmful.** With a real outside flood,
   an empty cell that is unreached is enclosed *by definition*. Applying
   "≥ 8 of 12 neighbours occupied" as a *membership* filter afterwards drops the
   interior cells of any void bigger than ~6 cells (a cell in the middle of a
   void has void neighbours, not occupied ones), which splits one genuine cavity
   into several disconnected shells of outer cells. `cavity_cells_are_empty_and_connected`
   still passes — each fragment is connected. Directly threatens criterion 6:
   §8.5 needs *two* cavities, and one void split in three reads as an enzyme.
2. **The flood is too strict in the other direction.** §8.4's own figure says a
   cavity is "reachable through a mouth — a real binding site". A hermetically
   sealed void is not reachable by any substrate. Requiring zero lattice
   connectivity to the outside means the only cavities reported are ones nothing
   can enter.

**Fix** — one change addresses both: flood with a *probe*, not with a point.
```rust
/// Width of the solvent probe, in occupied-neighbour terms. The flood may
/// pass through a cell only if the cell has room around it; a constriction
/// (many occupied neighbours) is a mouth, and the pocket behind it is
/// enclosed even though a zero-radius lattice path exists.
///
/// This is a probe size — a physical quantity — not a tuning knob: it is the
/// size of the smallest thing that counts as solvent. Deriving it from the
/// solvent species' own extent is the right long-term answer.
const PROBE_MAX_OCCUPIED: usize = 5;

// ...in the outside flood:
if in_box(nb, lo, hi) && !occupied.contains_key(&nb)
    && occupied_neighbours(&occupied, nb) <= PROBE_MAX_OCCUPIED
    && outside.insert(nb) { q.push_back(nb); }
```
then take interior components **directly** (`empty ∧ !outside`), with
`MIN_ENCLOSURE` demoted to a *per-component* summary rather than a membership
gate.

**Tests:** the plan's `an_open_groove_is_not_reported_as_a_cavity` is good and
should stay; add (a) a hollow shell with a one-cell mouth, which must still be
reported, and (b) an assertion that a 12-cell spherical void is reported as
**one** cavity, not several.

**Confidence:** derived. The fragmentation argument is arithmetic and certain;
the "sealed voids are rare in real folds" claim is a prediction — measure the
cavity-count distribution over `folded(seed, 80)` for 50 seeds before and after.

---

## G7 · `LATTICE_SPAN = SQRT_2` converts in the wrong direction and ignores the universe

**Where:** Task 13, `LATTICE_SPAN` (`chemistry.md:1181`) and its use at
`:1354–1356`; same unit reappears in Task 16's `sep` (`chemistry.md:2231`).
**Silently wrong.**

The comment is correct — the FCC nearest-neighbour separation is √2 *integer
units*. It is then used as a **multiplier on the coordinates**, which sets one
integer unit to √2 Span and the nearest-neighbour distance to **2.0** Span. To
make the NN distance equal to √2 you would multiply by 1; to make it a monomer
contact distance you divide by √2.

Worse, it is universe-independent. Molecule signatures scale with generated
radii (0.56–4.55); cavity signatures are pinned at a fixed ~2.8. `affinity(cavity,
molecule)` is therefore commensurable in *one* universe out of the 40 Task 20
iterates. M5 is only half fixed: the units now have a name, but not the right
magnitude.

```rust
// One integer lattice unit in Span. Two monomers in contact are `ideal_gap`
// apart (§8.3), and FCC nearest neighbours are sqrt(2) integer units apart,
// so one unit is ideal_gap / sqrt(2). A fixed constant here would put cavity
// distances on a different scale from molecule extents in every universe but
// one.
let lattice_unit = u.consts.ideal_gap / std::f64::consts::SQRT_2;
```

**Test:** assert `mean(cavity.signature.r)` and `mean(molecule.signature.r)` are
within a factor of ~3 of each other across 10 generated universes.
**Confidence:** derived, high.

---

## G8 · `may_bind` is a genuine bound — and rejects 100% of real pairs

**Where:** Task 10 Step 3 `may_bind` (`v0.md:4280`) and its test (`v0.md:4113`).

**The bound is valid, verified.** For any rotation, `Σ ds_i = ΣA.r + ΣB.r − D·gap`
because `anti∘perm` is a bijection, so it is rotation-invariant; Cauchy–Schwarz
gives `Σ ds² ≥ (Σ ds)²/D = D·mean_defect²`; dropping the charge term is valid
because it is a sum of squares. It therefore holds **for every rotation, not on
average**. 12,000 random pairs across three thresholds: **0 violations.**

Two caveats:

1. **The sign argument depends on `w_shape ≥ 0` and `w_charge ≥ 0`**, which
   happens to hold (`next_f64_range(0.6, 1.4)`, `v0.md:2047–2048`) but is
   asserted nowhere. Add `debug_assert!(k.w_shape >= 0.0 && k.w_charge >= 0.0,
   "the Cauchy–Schwarz ceiling in may_bind assumes both weights are non-negative")`.
2. **Given G1, it rejects everything.** On 1,770 pairs of real embedded 8-atom
   signatures it rejected **100%** at thresholds −40, −100, −300 and −1000 — with
   zero violations, because the true affinities are ~−3000. Once G1 is fixed the
   picture changes; until then `may_bind` is a "nothing binds" switch.

**The test does not test what its comment says.** `v0.md:4130-4132`:
```rust
// And it must actually reject something, or it is not a filter.
assert!(may_bind(&a, &complement_of(&a, &g, &k), &k, threshold));
```
That asserts `may_bind` **accepts** a perfect complement. M10's actual complaint
— "rejects 0.000% of pairs" — is not retested. Add a rate assertion:
```rust
let rate = pairs.iter().filter(|(x, y)| !may_bind(x, y, &k, threshold)).count() as f64
         / pairs.len() as f64;
assert!((0.02..0.98).contains(&rate), "pre-filter rejects {rate:.1%} — it is not filtering");
```

I also tested a strictly tighter bound (rearrangement inequality on the sorted
extents, which is *exact* over all bijections rather than just Cauchy–Schwarz):
4.0% vs 3.3% rejection on random pairs. Not worth the complexity — **withdrawing
that suggestion.**

**Confidence:** measured.

---

## G9 · `the_mirror_complement_does_not_fit` passes with a 4-orders-of-magnitude margin deficit

**Where:** Task 10 Step 1 (`v0.md:4064-4077`).

Answering the question directly: **yes, this is the right pair of tests, and no,
they cannot both pass with a wrong kernel.** Verified numerically on the plan's
own seed (6-atom chain, elements i%4, D=42, gap 2.0):

| | physical complement | mirror complement |
|---|---|---|
| kernel with `anti` | **−0.000e0** | −7.755e−1 |
| kernel without `anti` | −7.755e−1 | **−0.000e0** |

The roles swap exactly, so each test fails if `anti` is dropped. They are also
mutually self-checking: if the seed signature were antipodally symmetric the two
complements coincide and the two assertions contradict each other, so a
degenerate seed fails loudly rather than silently.

**But the seed is weak and the threshold is far too loose.** The mirror
complement scores −0.776 against a `self` affinity of −1395 — a discrimination
of **0.06% of the natural scale**, and the assertion threshold is `−1e-6`, three
more orders below that. A near-linear chain is nearly antipodally symmetric
(`max|r[i] − r[anti[i]]| = 0.57` against extents 1.46–5.83). Over 60 random
8-atom trees the mirror penalty ranges 3.4 to 101, median 40 — i.e. a branched
seed gives 50× the margin.

```rust
// Margin relative to the signature's own scale, not an absolute 1e-6. A
// near-linear chain is nearly antipodally symmetric: on chain(6) the mirror
// penalty is 0.06% of the self-affinity, so a 1e-6 threshold would keep
// passing long after the test stopped discriminating anything.
let self_scale = affinity(&a, &a, &g, &k).abs();
assert!(affinity(&a, &mirror, &g, &k) < -1e-3 * self_scale);
```
and seed from a branched molecule, plus a property test over 200 random graphs.

**Confidence:** measured.

---

## G10 · The SMACOF target matrix is not a metric, and ignores bond order

**Where:** Task 8 Step 3 (`v0.md:3492`). Not covered by the previous review.

`target[i][j] = hops[i][j] * (radius(i) + radius(j))` uses the **endpoints'**
radii for a multi-hop distance. The distance between atoms three bonds apart
should be the sum of the intervening bond lengths, not 3× the endpoint sum. With
`r_i = r_k = R` and a small `r_j` on the shortest path, `t(i,k) = 2R·h_ik` while
`t(i,j) + t(j,k) ≈ R·h_ik` — the triangle inequality is violated by **2×**.

SMACOF converges anyway (it minimises stress against whatever you give it), so
this is silent. The effect is that shape is driven by which two atoms are at the
ends rather than by what the path is made of, which flattens the shape diversity
the whole chemistry samples.

Also: `hop_distances` uses `m.neighbours()`, which merges all three bond-order
planes. **A double bond and a single bond produce identical geometry.** Bond
order changes species identity (`CanonForm.planes`) but not shape, signature, or
binding — so the three bond orders are chemically inert as far as §8.3 is
concerned.

```rust
// Weighted all-pairs shortest path: each bond (a,b) has length r_a + r_b,
// scaled down for higher orders. n <= 12, so Floyd-Warshall is free.
// hops * (r_i + r_j) uses the *endpoints'* radii for a multi-hop distance,
// which violates the triangle inequality by up to 2x when the intervening
// atoms are small — the target then is not a metric and SMACOF is fitting an
// inconsistent set of distances.
const ORDER_SHORTENING: [f64; 3] = [1.0, 0.88, 0.80];   // single, double, triple
let mut d = [[f64::INFINITY; MAX_ATOMS]; MAX_ATOMS];
for i in 0..n { d[i][i] = 0.0; }
for i in 0..n { for j in 0..n {
    let o = m.bond_order(i as u8, j as u8);
    if o > 0 { d[i][j] = (radius(i) + radius(j)) * ORDER_SHORTENING[(o-1) as usize]; }
}}
for k in 0..n { for i in 0..n { for j in 0..n {
    if d[i][k] + d[k][j] < d[i][j] { d[i][j] = d[i][k] + d[k][j]; }
}}}
```
A shortest-path metric satisfies the triangle inequality by construction, so the
target becomes well-posed.

**Test:** `assert!(target[i][k] <= target[i][j] + target[j][k] + 1e-12)` for all
triples over 500 random molecules; and `assert_ne!(embed(single_bonded), embed(double_bonded))`.

**Confidence:** the triangle-inequality violation is arithmetic and certain. The
*size* of the shape effect is not measured — worth measuring before treating the
`ORDER_SHORTENING` values as physics.

---

## G11 · Smaller items

- **`SHELL = 0.75` is absolute** (`v0.md:3807`) while extents scale with generated
  radii (0.56–4.55). In a small-radius universe the shell spans the whole
  molecule and `a[d]` collapses to a global average; in a large-radius one only
  the single frontmost atom contributes. Scale it: `SHELL_FRACTION * mean_radius`.
- **The cavity flood can wrap the grid.** `in_box` decodes with
  `fcc::from_flat`, which is `f % GRID` — for a cell at `x = 0` stepping to
  `x = −1`, the decode returns `x = 63` and, if the bounding box reaches 63, the
  flood leaks around the torus. Same class as S5, now in Task 13. Add
  `debug_assert!(lo.0 >= 1 && hi.0 <= fcc::GRID - 2, ..)` in `bounding_box`, or
  carry `(x,y,z)` through the flood instead of flat indices.
- **`Cavity::centre` is documented as "nearest the cavity's centre of mass" but
  is `cells[cells.len()/2]`** — the median of the *sorted flat indices*, i.e. the
  median in z-major order, which is not geometrically central. It sets the origin
  for every `r[d]`. Use the centroid and pick the cell minimising distance to it.
- **`energy_scale` in Task 11's Metropolis (`chemistry.md:519`) is never
  defined.** It sets the acceptance temperature relative to the contact energies;
  wrong by an order of magnitude in either direction and the anneal either
  freezes at step 1 or never converges. It needs a stated derivation (e.g. the
  RMS of `class_energy_table`), not a literal.
- **Task 16's prose says four constants were replaced; the code below it still
  uses `BIND_THRESHOLD`, `IDEAL_SEPARATION`, `MAX_ENHANCEMENT` and `* 1.0e3`**,
  none of which are defined in the shown code (only `BIND_FRACTION` is). The fix
  landed in the comment and not the body.
- **Undefined test helpers in Task 8**: `random_tree`, `permute`, `dist` and
  `kabsch_rmsd` are all used and none is provided. `FoldWorkspace::fold_seed` is
  an associated function but called as a free `fold_seed(p, u)`.

---

## What I checked and found sound

**The rotation table — exhaustively, at all three resolutions.** Rebuilt the
plan's `icosahedron`/`faces_of`/`intern`/`build_perms`/`build_anti` and checked:

- vertex counts exactly 12 / 42 / 162, all unit, midpoints deduplicated;
- **60 distinct permutations** at each resolution;
- **closed under composition** — all 3,600 products at each resolution land in
  the set;
- **inverse-closed**; **identity present exactly once**;
- every permutation is a genuine **proper rotation**: reconstructing the 3×3
  matrix by least squares over all D directions gives `det = +1.000000000000`
  and `RᵀR − I` within 1.6e−15, with reconstruction residual ≤ 1e−15. **Zero
  improper elements.** 60, not 120 — §22.8 holds.
- `anti` is an involution with `dirs[anti[i]] == −dirs[i]` to 1e−12 (in fact
  bit-exact, since midpoint negation and normalisation are exact under IEEE-754);
- **`anti[perm[i]] == perm[anti[i]]` for all 60×D pairs at all three
  resolutions** — exactly equal as integers, so the two forms are
  interchangeable and the choice does not matter. The plan's
  `antipode_commutes_with_every_rotation` asserts precisely this.

**The `anti` fix genuinely searches proper rotations.** `dirs[anti[perm_R[i]]] =
−R·dirs[i]`, so the term pairs A's extent along `d` with B's extent along
`−R d`; that is B rotated by `R⁻¹`, and the group is inverse-closed, so
maximising over `R` maximises over all 60 proper rotations of B. Correct.

**`affinity_ordered`'s symmetry argument survives the change.** With
`σ = anti ∘ perm_R`, the per-direction terms are symmetric, so
`score_R(a,b) = score_{σ⁻¹}(b,a)`; and `σ⁻¹ = perm_{R⁻¹} ∘ anti = anti ∘
perm_{R⁻¹}`, which is in the same set because `anti` is an involution that
commutes with every rotation. So the maxima agree. Verified numerically: the
difference is exactly 0.0 on the plan's seed.

**Task 6's densified colouring is correct.** The O(n³) loop counts *distinct*
element values strictly below `elem[v]` — the `!(0..w).any(|x| elem[x] ==
elem[w])` guard restricts the count to first occurrences. That is a dense rank in
`0..n_distinct ≤ n < MAX_ATOMS`, and it depends only on the multiset of element
ids, so it is isomorphism-invariant. It genuinely fixes S6. `refine`'s rank loop
has the identical structure and is likewise correct. Individualization preserves
density (`k` classes → `k+1`, max colour ≤ n−1), so `masks`, `sig` and
`count_classes` are all within their array bounds at every recursion depth.

**FCC neighbours and parity.** The twelve offsets are exactly the `(±1,±1,0)`
permutations; every one changes `x+y+z` by 0 or ±2, so even parity is preserved.

**The FCC parity constraint does not break Task 13's outside flood.** `box_cells`
does enumerate odd-parity triples that are not lattice sites, but the flood
propagates only through `fcc::NEIGHBOURS`, which preserves parity — so the odd
and even sublattices flood independently and neither leaks into the other. All
occupied cells are even-parity (the start is at (32,32,32), sum 96), so an
odd-parity empty cell always has **zero** occupied neighbours and can never pass
`filled >= MIN_ENCLOSURE`. Odd-parity cells cannot produce spurious cavities and
cannot leak the flood into a genuine one. The enclosure test is consistent
because occupancy, flood and candidates all live on the same sublattice.

**Task 13's two `lining_signature` passes are consistent.** `reach_of` is a pure
closure of `m` and the loop's `dir`, evaluated identically in both passes, and
`nearest = INFINITY` is handled (both `r` and `a` fall to 0.0). The defect is
*which* quantity it computes (G4), not the two-pass structure.

**The fold workspace trail reset is sound.** `clear` zeroes every cell in
`trail`; `occupy` pushes on every write; accepted moves push the new cell; and
rejected moves restore `occ[target] = 0` without ever having pushed it. No cell
can stay set across folds. `workspace_reuse_matches_a_fresh_workspace` is the
right test for it.

**Incremental contacts equal a full recount, exactly.** A move touches only
monomer `i`, `removed` is computed before `occ` is mutated and `added` after,
both count each contact once (matching `count_contacts`' divide-by-two), and
`recount_contacts` counts each pair once via `j > i`. Integer histogram, so the
assertion is exact rather than tolerant — correct call.

**`propose` cannot break self-avoidance or connectivity.** End moves land on a
neighbour of the anchor; corner flips are filtered to cells adjacent to both
chain neighbours and are checked against `coords[i]`; both go through
`is_free(target)` and `in_bounds(target)`. `target == old` is rejected by
`is_free`. Sound — the problem is the *start* (G3), not the move set.

**`Signature::permuted` is correctly equivariant**: `out.r[perm[i]] = self.r[i]`
matches `dirs[perm[i]] = R·dirs[i]`, so it implements rotating the molecule by R.
`sorted_extents` is genuinely rotation-invariant, and `total_cmp` gives a total
order with no tolerance.

**`intern`'s 1e−18 squared tolerance has ~8 orders of margin** — the closest
distinct vertices at level 2 are ~0.27 apart, and shared midpoints are in fact
bit-identical because IEEE addition is commutative. `build_anti`'s 1e−24 and
`build_perms`' 1e−12 likewise have enormous margin (the observed reconstruction
residual is ~1e−15 squared, i.e. ~1e−30).
