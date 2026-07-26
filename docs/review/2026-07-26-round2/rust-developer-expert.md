# rust-developer-expert — V0 plan re-review (post-fix pass), 2026-07-26

Scope: both plan files, concentrating on code changed by the fix pass.

---

## 1. Task 14 — `SpeciesKind` exists as a declaration and nothing else. The S2 fix is half-applied.

**Where:** chemistry.md:1634–1782 (`kind.rs` / `species.rs`)

`SpeciesKind` is declared, and then:

- `SpeciesKind` is never constructed anywhere in the workspace.
- `intern` pushes `SpeciesRecord { canon, sig, mass, bond_count, cleave_propensity, solvent_rate }` — `canon` is not a field, and `kind`, `radiogenic_rate`, `fold`, `cavities` are all missing. It does not compile.
- `let exposed_fraction = match &fold` reads a binding `fold` that is never created.
- `FoldId` occurs exactly once in the plan: as the type of the field. It is never defined.
- `intern_polymer` is called by `burial_lowers_the_solvent_rate` and never defined.
- `Interner::ids` is `BTreeMap<CanonForm, SpeciesId>`. A polymer has no `CanonForm`, so **there is no identity path for `SpeciesKind::Polymer` at all** — which is precisely the hole S2 was raised to close.

Even completed, the shape moves the problem rather than removing it: `kind: Small(_)` with `fold: Some(..)` and a non-empty `cavities` is fully representable, and `Polymer` with `fold: None` is too.

**Cost:** the §8.6 regression S2 identified is still open. An implementer reaching Task 16 needs `record.cavities` for a polymer, finds no way to get one, and computes it at the point of need — inside the step loop. That is the four-orders-of-magnitude path S2 described.

**Fix:**

```rust
pub enum SpeciesKind<const D: usize> {
    Small { canon: CanonForm },
    Polymer { units: Vec<u8>, fold: FoldId, cavities: Vec<Cavity<D>> },
}

/// Interner key. Polymers key on their sequence, small molecules on their
/// canonical form; the two spaces are disjoint by construction.
#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum SpeciesKey { Small(CanonForm), Polymer(Vec<u8>) }

pub struct SpeciesRecord<const D: usize> {
    pub kind: SpeciesKind<D>,
    pub sig: Signature<D>,
    pub mass: Mass,
    pub bond_count: u32,
    pub radiogenic_rate: f64,
    pub cleave_propensity: f64,
    pub solvent_rate: f64,
    // no `fold`, no `cavities` — they live in the variant that can have them
}

impl<const D: usize> SpeciesRecord<D> {
    #[must_use]
    pub fn cavities(&self) -> &[Cavity<D>] {
        match &self.kind {
            SpeciesKind::Polymer { cavities, .. } => cavities,
            SpeciesKind::Small { .. } => &[],
        }
    }
}
```

Folding is total over polymers, so `fold: FoldId` not `Option<FoldId>`; "small molecule has no fold" becomes unrepresentable rather than a convention. `intern` and `intern_polymer` both need writing against `SpeciesKey`.

**Severity: bug risk.** Task 14 cannot be implemented from what is written.

---

## 2. Task 16 — the module header describes a fix the code does not contain.

**Where:** chemistry.md:2180–2254

The header says `1.0e3` "is now `UniverseConsts.catalytic_prefactor`" and that `BIND_THRESHOLD` has been replaced by a dimensionless `BIND_FRACTION`. The code below it:

- still has the literal `* 1.0e3` (line 2249);
- never reads `catalytic_prefactor`;
- defines `BIND_FRACTION` and never uses it;
- uses `BIND_THRESHOLD`, which is not defined in the snippet;
- still clamps to `MAX_ENHANCEMENT`, which the header says was inert.

**Cost:** M4 is not fixed. The number deciding whether exit criterion 6 reports 2× or 500× is still a literal, and the header now says otherwise — so the next reader will not check. A false "already fixed" comment is worse than no comment.

**Fix:** `let factor = 1.0 + proximity * enclosure * fit * u.consts.catalytic_prefactor;`, thread `&UniverseConsts` in, replace `bind_a < BIND_THRESHOLD` with `bind_a < BIND_FRACTION * reference_score`, delete `MAX_ENHANCEMENT`.

**Severity: bug risk.**

---

## 3. Task 10 — the ordered/unordered split reintroduces the disagreement `affinity_with_rotation` exists to prevent.

**Where:** v0.md:4210–4270

`affinity_with_rotation`'s doc: the renderer must not re-derive the winning rotation with its own loop, because a `>=` instead of `>` picks a different winner on ties and "the diagram would then disagree with the physics".

But `affinity_ordered` may swap the pair, and **discards the fact that it swapped**. The physics goes through `affinity_ordered(a, b)`; the renderer must call `affinity_with_rotation(a, b)` to get an index. When `a.lex_cmp(b) == Greater` those two evaluate different loops and can select different rotation indices on ties. Same defect, arrived at by a different route.

Both `affinity` and `affinity_with_rotation` are `pub`, so "every caller must use `affinity_ordered`" is enforced by a doc comment, which is not an enforcement mechanism.

**Fix:** one public entry point that carries everything:

```rust
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bind {
    pub score: f64,
    /// Index into `Geodesic::perms`, valid for the pair *as evaluated*.
    pub rotation: u8,
    /// True if the pair was swapped before evaluation; the renderer must
    /// apply `rotation` to the other member when this is set.
    pub swapped: bool,
}

#[must_use]
pub fn bind<const D: usize>(a: &Signature<D>, b: &Signature<D>, g: &Geodesic<D>, k: &BindConsts) -> Bind {
    match a.lex_cmp(b) {
        std::cmp::Ordering::Greater => { let (score, r) = kernel(b, a, g, k); Bind { score, rotation: r, swapped: true } }
        _ => { let (score, r) = kernel(a, b, g, k); Bind { score, rotation: r, swapped: false } }
    }
}

/// Unordered. Private: the two argument orders sum in different sequences.
fn kernel<const D: usize>(...) -> (f64, u8) { ... }
```

Keep `affinity_ordered` as `bind(..).score` if the call sites read better. Delete `pub fn affinity` and `pub fn affinity_with_rotation`. Also rename `lex_cmp_pub` — `pub` in an identifier is a smell; make the private one private and call the public one `lex_cmp`.

**Severity: bug risk** (renderer/physics divergence), plus maintenance.

---

## 4. Task 8 — `canonicalise_frame`'s collinear early return is dishonest.

**Where:** v0.md:3604–3649

```rust
if !found {
    return; // collinear molecule; the axis alone is already canonical
}
```

It is not. No rotation has been applied at that point, so a collinear molecule keeps whatever orientation SMACOF happened to produce — exactly the instability the function exists to remove, silently, for the case where it is easiest to fix. `embed` returns the result as canonical either way. Diatomics and linear triatomics are not edge cases in this chemistry.

The `everything coincident` return is honest (there is genuinely nothing to orient). The `n < 2` return is honest.

**Fix:** complete the frame instead of bailing. The residual freedom for a collinear molecule is a rotation about `e0`, so pick a deterministic completion:

```rust
    // Collinear: the frame is determined only up to a rotation about e0, so
    // complete it deterministically from the standard basis. Returning here
    // instead would leave the molecule in SMACOF's arbitrary frame, which is
    // the exact instability this function exists to remove.
    let e1 = if !found {
        let axis = if e0[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let dot = axis[0]*e0[0] + axis[1]*e0[1] + axis[2]*e0[2];
        let perp = [axis[0]-dot*e0[0], axis[1]-dot*e0[1], axis[2]-dot*e0[2]];
        let m = norm(perp);
        [perp[0]/m, perp[1]/m, perp[2]/m]
    } else { e1 };
```

Separately: the signature should be `fn canonicalise_frame(pos: &mut [[f64; 3]])`, called as `canonicalise_frame(&mut pos[..n])`. The `(array, n)` pair is an invariant nothing enforces, and `pos[a0]`/`pos[i]` panic if `n > MAX_ATOMS`. No `Result` needed once the collinear case is handled — there is no failure left to report.

**Severity: bug risk.**

---

## 5. Task 11 — `fold` does not compile, and `fold_with_seed` wants to be private, not documented.

**Where:** chemistry.md:383–697

Compile blockers introduced by the fix pass:

- `energy_scale` (line 519, the new Metropolis criterion) is never defined.
- `fold_seed` is defined as `FoldWorkspace::fold_seed` but called as a free function (line 418).
- `Domain` is used (line 384) but only `Stream` is imported (line 311).
- `local_contacts`'s `coords` parameter is never read → `unused_variables` → `-D warnings` → CI red.
- `recount_contacts`: `(k * 0) + class_of(p, k)` trips `clippy::erasing_op`, which is **deny by default**, not pedantic.

**On the question asked:** `fold_with_seed` needs no marker, no sealed trait, no `#[cfg(test)]`. Its only caller is a test in the same file, and tests in a child module can call a private method. Make it `fn fold_with_seed(...)` — private — and the doc comment stops being a request and becomes a fact. Delete the "Diagnostics only" paragraph; the visibility says it.

The doc on `fold_seed` claims it "shares its construction with the fold cache key (Task 12) so the two cannot drift apart". It does not — `key_of` (chemistry.md:922) is an independently written copy of the same mixing loop with different constants and stream indices. Either extract `fn mix_sequence(seed: u64, units: &[u8]) -> u64` and use it in both, or delete the claim.

**Severity: bug risk** (blocks the task), maintenance (the false sharing claim).

---

## 6. Three `.powi(` calls in library code will fail the plan's own CI gate.

**Where:** v0.md:1722 (`stability`), v0.md:3884 (`Signature::distance`), chemistry.md:2229–2231 (`catalysis_factor`)

`xtask`'s `BANNED_CALLS` includes `".powi("` (v0.md:940), and `det_math.rs`'s own doc says why: "`powi`'s multiply tree is chosen by LLVM, so it can be re-associated by a compiler upgrade without our source changing". All three sites are library code outside `det_math.rs`. CI goes red at Task 4, again at Task 9, again at Task 16 — and the tempting fix each time is to weaken `BANNED_CALLS`.

`stability` is the one that matters: it feeds radiogenic decay, which is a discrete propensity channel.

**Fix:**

```rust
let s = f64::from(index) / 90.0;
let s3 = s * s * s;
let stability = (s3 * s3).min(0.9);
```

and `let (dr, da) = (self.r[i] - other.r[i], self.a[i] - other.a[i]); acc += dr * dr + da * da;`, and `let (dx, dy, dz) = (f64::from(ax - bx), ...); let sep = (dx*dx + dy*dy + dz*dz).sqrt();`.

**Severity: bug risk** (determinism) + CI.

---

## 7. `det_math` is half-moved. Tasks 10, 14, 15, 21 still point at `borbax-molecule`.

The File Structure (v0.md:51, 105) and Task 2 Step 6 are correct: `borbax-units/src/det_math.rs`, backed by `libm`. But:

- v0.md:4324 — Task 10's prose still instructs "create `crates/borbax-molecule/src/det_math.rs`", with a stub delegating to the platform and a comment referencing a **Task 25 that does not exist** in a Tasks 1–21 plan. The code six lines above already calls `borbax_units::det_math::exp`.
- chemistry.md:1898, 2046, 2067 — `borbax_molecule::det_math::exp` in `borbax-reaction`.
- chemistry.md:3095 — Task 21's file list names `crates/borbax-molecule/src/det_math.rs`.
- chemistry.md:3164 — "replace the det_math bodies with libm" at Task 21, already done at Task 2.

Compounding: the xtask exemption matches **by file name anywhere in the tree** (v0.md:949), so following Task 10's stale instruction creates a second `det_math.rs` that is silently permitted to call the platform. Make the exemption path-based:

```rust
if entry.ends_with("borbax-units/src/det_math.rs") { continue; }
```

**Severity: bug risk** (determinism, and the guard against it is bypassable).

---

## 8. Task 14 — `rate()` returns a plausible default, and `catalysis` now has two sources.

**Where:** chemistry.md:1866–1902

```rust
let a = conc.get(r.reactants[0].0 as usize).copied().unwrap_or(0.0);
```

A species id outside the concentration vector yields rate 0 — the reaction quietly never fires, and nothing distinguishes "no reactant present" from "the vectors are out of sync". Loud failure turned into a silent wrong answer.

And `Reaction` now carries `catalysis: f64` (the S2 fix) while `rate` still takes `catalysis: f64` as a fourth parameter and ignores the field. Two sources, one of which is authoritative and neither of which says which.

**Fix:** drop the parameter, read `r.catalysis`. For the lookup:

```rust
#[must_use]
pub fn rate(r: &Reaction, t: Thermal, conc: &Concentrations) -> f64
```

where `Concentrations` is a newtype over `Vec<f64>` indexed by `SpeciesId` and sized from the interner, so the out-of-range case cannot arise. Failing that, `debug_assert!(r.reactants[0].0 as usize) < conc.len())` at minimum.

Also: `Reaction` derives `PartialEq` and now contains an `f64`. Ask what equality of two channels is supposed to mean; if it is identity, compare `(kind, reactants, products, catalyst)`.

**Severity: bug risk.**

---

## 9. Task 13 — `cavities(&Fold, &Polymer)` accepts a fold and a polymer that never met, and the new test does exactly that.

**Where:** chemistry.md:1200–1205, 1329–1334; test at 1100–1104

Nothing links a `Fold` to the `Polymer` it was folded from. `lining_signature` builds `pos_of` sized `p.len()` from `occupied`, guarding the write (`if m < pos_of.len()`) and then reading unguarded (`pos_of[m]` inside `reach_of`) — so a mismatch either panics or silently reads `(0,0,0)` for the missing monomers and produces a plausible, wrong signature.

The new enclosure test passes `remove_one_monomer(&sealed)` together with the *original* `p`. If that helper drops a coordinate, the test is measuring a fold/polymer mismatch, not an open groove.

**Fix:** make the mismatch unrepresentable.

```rust
/// A fold and the polymer it was computed from. `FoldWorkspace::fold` is the
/// only constructor, so a fold can never be paired with a different chain.
pub struct Folded {
    polymer: Polymer,
    fold: Fold,
}
impl Folded {
    #[must_use] pub fn fold(&self) -> &Fold { &self.fold }
    #[must_use] pub fn polymer(&self) -> &Polymer { &self.polymer }
}

pub fn cavities<const D: usize>(f: &Folded, g: &Geodesic<D>, u: &Universe) -> Vec<Cavity<D>>
```

Then the enclosure test has to construct a breached *chain* and fold it, which is the experiment it was trying to run.

**On the rewrite:** `reach_of` reads well and captures exactly `pos_of` (by ref), `centre` and `dir` (both `Copy`) — nothing more. The two-pass structure is right and the comment explaining why merging them is a semantic change is the kind of comment worth keeping. Two smaller things:

- `Cavity::centre` is documented as "the cell nearest the cavity's centre of mass"; the code takes `cells[cells.len() / 2]`, the median in sort order. Different thing. And `lining_signature` recomputes the same expression instead of taking the field. Compute once, pass it, fix the doc.
- `sig.r[dir_idx] = if nearest.is_finite() { nearest } else { 0.0 }` — 0.0 means "wall touching this direction", the *tightest* possible fit, but the case being encoded is "no lining monomer that way", i.e. open space. That inversion feeds straight into `a.r[i] + b.r[j] - ideal_gap`. Either use the cavity's own max half-extent or make it explicit. (Flagging the type/default issue; the right physical value is geometry-numerics' call.)

**On `bounding_box` returning `((i32,i32,i32),(i32,i32,i32))`:** yes, that wants a type — mostly because of `in_box(cell, lo, hi)`, `box_cells(lo, hi)`, `boundary_cells(lo, hi)`, where `(hi, lo)` compiles and yields an empty box, so cavity detection returns nothing and no test fails.

```rust
#[derive(Clone, Copy)] pub struct Cell3 { pub x: i32, pub y: i32, pub z: i32 }
#[derive(Clone, Copy)] pub struct Bounds { lo: Cell3, hi: Cell3 }
impl Bounds {
    fn of(coords: &[i32]) -> Option<Self>   // None for an empty chain
    fn contains(self, c: Cell3) -> bool
    fn cells(self) -> impl Iterator<Item = i32>
    fn faces(self) -> impl Iterator<Item = i32>
}
```

`Option` on the constructor also removes the empty-`coords` case, which currently produces `lo = (i32::MAX-2, ..)` and relies on a guard on `f.n` in a different function.

**Severity: bug risk** (the mismatch), maintenance (the tuples).

---

## 10. Task 12 — four `unreachable!()` in library code, shipped with a note saying to fix them later.

**Where:** chemistry.md:962, 984, 1002, plus the note at 1013

The note is honest, and it will be skipped. `unreachable!()` is worse than the `unwrap()` the Global Constraints ban — same panic, and it reads as an assertion of impossibility rather than an admission. There is a version that compiles with neither:

```rust
pub fn get_or_fold(&mut self, p: &Polymer, u: &Universe) -> &Fold {
    let k = key_of(p, u);
    if let Some(&slot) = self.map.get(&k) {
        if let Some((_, e)) = self.slots.get_mut(slot).and_then(Option::as_mut) {
            self.stats.hits += 1;
            e.referenced = true;
            return &e.fold;              // NLL reborrows the &mut as &
        }
    }
    self.stats.misses += 1;
    let fold = self.ws.fold(p, u);
    let bytes = /* as written */;
    while self.stats.bytes + bytes > self.budget && !self.map.is_empty() {
        self.evict_one();
    }
    let slot = self.free.pop().unwrap_or_else(|| { self.slots.push(None); self.slots.len() - 1 });
    self.stats.bytes += bytes;
    self.map.insert(k, slot);
    // `Option::insert` hands back a &mut to what it just stored, so there is
    // no re-lookup and nothing to unwrap.
    let (_, entry) = self.slots[slot].insert((k, Entry { fold, bytes, referenced: false }));
    &entry.fold
}
```

`evict_one` similarly: `let Some((k, e)) = self.slots[idx].take() else { continue };` then use `k` and `e` directly.

Two more in the same file:

- `key_of` shadows its `u: &Universe` parameter with `for (i, &u) in p.units()` — the same identifier is a `&Universe` above the loop and a `u8` inside it. It happens to work because `u.seed` is read first. Rename the loop binding.
- `evict_one` spins forever if every slot is `None` (`else { continue }` with no progress). It is only reachable behind `!self.map.is_empty()`, so it is safe today and hangs the first time someone calls it from anywhere else. Bound the sweep to `2 * self.slots.len()` and return.

**Severity: bug risk** (the `unreachable!`s land in committed code; the hang), maintenance (the shadow).

---

## 11. Task 7 — `GeoError` is hand-rolled, and its main variant is a compile-time fact reported at runtime.

**Where:** v0.md:2936–2943, 3057–3198

CLAUDE.md's Dependencies section names `thiserror` for library errors and Task 1 already declares it at v0.md:179. `GeoError` has no `Display`, no `std::error::Error`, so it cannot be `?`-ed into `anyhow` in `borbax-cli` and cannot be printed usefully. Add:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum GeoError {
    #[error("unsupported signature resolution {0}; must be 12, 42 or 162 (spec §22.2)")]
    UnsupportedResolution(usize),
    #[error("a rotation failed to map the vertex set onto itself — this is a bug in the geodesic construction, not bad input")]
    RotationNotClosed,
}
```

Better: `UnsupportedResolution` should not exist. `D` is a const generic, so `Geodesic::<43>` is knowable at compile time, and the project already uses `compile_fail` doctests for unit mixing.

```rust
mod sealed {
    pub trait Resolution { const LEVEL: u32; }
}
pub struct Res<const D: usize>;
impl sealed::Resolution for Res<12>  { const LEVEL: u32 = 0; }
impl sealed::Resolution for Res<42>  { const LEVEL: u32 = 1; }
impl sealed::Resolution for Res<162> { const LEVEL: u32 = 2; }

impl<const D: usize> Geodesic<D> where Res<D>: sealed::Resolution {
    pub fn build() -> Result<Self, GeoError> { /* only RotationNotClosed remains */ }
}
```

`Geodesic::<43>::build()` then fails to compile, `level_for` and its `Option` go away, and the surviving `Result` means one thing: the group construction is broken. The existing runtime test becomes a `compile_fail` doctest.

**Severity: maintenance** (the error type), **polish→maintenance** (the sealed bound — it removes a whole error path).

---

## 12. Dependencies

**`insta`** (Task 19, chemistry.md:2784–2842) — correct, dev-dependency, `cargo insta review` workflow described properly. One thing to finish: `CLAUDE.md`'s Commands block still documents `UPDATE_GOLDENS=1 cargo test -p borbax-render`. Replace with `cargo insta test -p borbax-render && cargo insta review`, or the hand-rolled mechanism comes back through the door the docs left open.

**`criterion`** (Task 20, chemistry.md:2868–2873) — correct, dev-dependency, and the reasoning about a timing assertion inside `#[test]` being a flaky test that gets deleted is right.

**`levenberg-marquardt`** (Task 20b, chemistry.md:3058–3061) — right crate, **wrongly priced**. It is described as "a cheap dependency despite appearances", but it is a *runtime* dependency of `borbax-beaker`, and it pulls in **`nalgebra ^0.34`** plus `num-traits` and `cfg-if`. `borbax-beaker` is the simulation crate. Adding nalgebra to its dependency graph so that a post-hoc curve fit can run means every simulation build links it, and the "does not enter golden hashes" property is asserted rather than enforced.

Put it where it cannot reach the physics:

```
crates/borbax-analysis/   # activity.rs, shadow.rs, plateau.rs
    depends on borbax-beaker (types only), levenberg-marquardt
crates/borbax-beaker/     # unchanged dependency graph
```

`borbax-cli` depends on both. Then "cannot affect simulation output" is a fact about the dependency graph rather than a claim in a comment, and `cargo tree -p borbax-beaker | grep nalgebra` is the test.

Also correct the claim at chemistry.md:3075: routing the *models'* `exp`/`powf` through `det_math` does not make the fit reproduce cross-platform. The solver's own QR decompositions run through nalgebra and are not routed. That is fine — the verdict is not golden — but the sentence as written will be cited later as evidence that a nalgebra-backed path is cross-platform-safe, which is the same class of doctrine error as the vector-maths determinism claim §18.2 just corrected.

**Still hand-rolled and defensible:** SMACOF (fixed 240 iterations is the requirement no solver crate offers), three-component vector maths (v0.md:2981–2994 — twenty lines, and §18.2 now records this as a size judgement, correctly), colour refinement, the FCC lattice, the CLOCK cache.

---

## 13. Task 4 — `(ShellPattern, Vec<Element>)`

A two-field tuple destructured at its single production call site (`let (shell, elements) = ...`, v0.md:2030) is defensible. What makes it not quite right is that the tests reach in positionally — `generate_elements(seed).0.periods` (v0.md:1551), `generate_elements(seed).1.len()` (v0.md:1543) — and that `Universe` already has fields named `shell` and `elements`, so the tuple is a nameless copy of a struct that exists. A third element is also already implied (the per-universe constants drawn at 1657–1671 are discarded).

```rust
/// The generated table and the shell pattern that produced it. Returned
/// together because the shell must come from the same stream draw — see
/// Task 5's note on re-deriving it.
pub struct GeneratedTable {
    pub shell: ShellPattern,
    pub elements: Vec<Element>,
}
```

**Severity: polish.**

---

## 14. Sweep — casts, panics, `#[must_use]`, transposition

- **chemistry.md:406, 634** — `self.occ[cell as usize]` and `ws.occ[(cell + d) as usize]` with no bounds check. `exposure` (line 534) calls `is_free(coords[i] + d)` with **no `in_bounds` guard**, unlike the anneal loop which has one. At the grid edge that is an index panic in library code. Fix: `fn occupant(&self, cell: i32) -> Option<u16> { usize::try_from(cell).ok().and_then(|c| self.occ.get(c)).copied() }`, with out-of-grid reading as free (solvent) for exposure and as unusable for moves.
- **chemistry.md:440** — `(((a + 1.0) * 0.5 * (N_CLASSES as f64 - 1.0)).round() as usize).min(N_CLASSES - 1)`. `f64 as usize` saturates at 0 for negatives and maps NaN to 0, both silently. Clamp in float space first: `.round().clamp(0.0, (N_CLASSES - 1) as f64) as usize`.
- **chemistry.md:401, 491, 528** — `monomer as u16 + 1` truncates above 65534 and can overflow. `MAX_POLYMER` makes it safe today; `u16::try_from(monomer).ok()?` or a `debug_assert!(MAX_POLYMER < u16::MAX as usize)` next to the declaration makes it stay safe.
- **chemistry.md:1717** — `SpeciesId(self.records.len() as u32)`. Silent wrap at 4 G species is not reachable, but the cast is free to make honest: `u32::try_from(self.records.len())`.
- **chemistry.md:1623, 1702** — `pub struct SpeciesId(pub u32)` with `record()` doing `&self.records[id.0 as usize]`. A public field means anyone can build `SpeciesId(9999)` and panic the interner, and ids from two interners are silently interchangeable. Make the field private with `#[must_use] pub fn index(self) -> usize`.
- **v0.md:4312** — `bind_probability(affinity: f64, t: Thermal, midpoint: f64)`. Two bare `f64` in positions 1 and 3, both affinity-valued; `bind_probability(midpoint, t, affinity)` compiles and returns nonsense. Combined with `may_bind(.., threshold: f64)` and the audit's own note that `ideal_gap` and `Signature.r[]` are bare `f64` (a G4 defeat), the answer is one newtype: `#[derive(Clone, Copy, PartialEq, PartialOrd)] pub struct Affinity(pub f64);` used for the score, the threshold and the midpoint. `may_bind(a, b, k, threshold)` itself is symmetric in `a`/`b` and separated from `threshold` by `&BindConsts`, so its transposition risk is low — the risk is in `bind_probability`.
- **v0.md:4280–4304** — `may_bind`'s bound is genuinely a bound (QM–AM on the shape term, non-negativity of the charge term), and it is a real improvement on M10. It relies on `w_shape >= 0.0 && w_charge >= 0.0`, which nothing checks. `BindConsts` should be constructed through `BindConsts::new(..) -> Result<Self, ConstsError>` rejecting negative weights, so the property the filter's soundness rests on is a property of the type.
- **chemistry.md:1806, v0.md:4194 etc.** — `#[must_use]` coverage is good: present on `affinity*`, `may_bind`, `bind_probability`, `embed`, `signature`, `cavities`, `catalysis_factor`, `Interner::{new, len, record}`, `AffinityMemo::new`. Missing on `FoldWorkspace::fold` and `fold_with_seed` (both return a `Fold` nobody would drop by accident, but they are pure transformations), and on `Fold::energy`'s callers.
- **chemistry.md:624–645** — `local_contacts` returns `Vec<(usize, u32)>`, allocated twice per anneal step, 40,000 allocations per fold. `[(usize, u32); 12]` plus a length, or `smallvec`. Per-species so not a §8.6 violation; handing the sizing question to `rust-performance-expert`.
- **v0.md:4299** — `let sum_a: f64 = a.r.iter().sum();` sits four lines below a comment declaring the accumulation shape in `affinity` to be load-bearing physics. `Iterator::sum` on `f64` is a sequential left fold and is deterministic, so this is not a determinism bug — but it *is* result-affecting (it gates a discrete accept/reject), and the inconsistency invites someone to "tidy" the pinned loop above it to match. Either write it as a pinned loop or add one line saying why it need not be. Deferring the call to `determinism-auditor`.

---

## What I read and found good

- **The `anti` fix in Task 10 is right, and the test pair is the right test pair.** `a_physical_complement_scores_zero` plus `the_mirror_complement_does_not_fit` together pin the handedness; either alone would not. The comment explaining that omitting `anti` converts the search into the improper elements is the single most valuable comment in the plan.
- **`may_bind` is now a genuine bound with a stated property**, and the property test asserts the property rather than an example. Big improvement on M10.
- **Task 8's frame-stability test asserts both aligned and unaligned bounds over random trees**, with the note that a chain does not reveal the bug. That is the right shape of test and the right note.
- **Task 13's outside-flood enclosure test** is a real test now — sealed shell versus breached shell — and the comment explaining why the previous one could not fail is worth keeping.
- **`Fold::energy` recomputed from an integer contact histogram** rather than accumulated as a running float is the right call, and the reasoning is stated.
- **The `AffinityMemo` correctness test** (memo versus direct computation) proves the thing everyone assumes.
- **The Metropolis fix in Task 11** is right, and drawing `roll` unconditionally so stream position never depends on a float comparison is exactly the discipline §13.1 needs. (Modulo the undefined `energy_scale`.)
- **Task 2's `det_math` on `libm`**, with the note that `sqrt` and the four operations stay native because IEEE-754 specifies those exactly, and the warning about `powi`'s multiply tree — correct and well argued.

## Rules I left alone, and why

- **Pinned accumulators** in `affinity` (two named accumulators summed in direction order, combined once), `Fold::energy`, `contact_energy`, `hop_distances`, the centroid loop in `embed`, and `lining_signature`'s two passes. Each is idiomatically a `.sum()` or a `fold`. Each is result-affecting and commented as load-bearing. Left as written; the comments are doing their job.
- **`BTreeMap`/`BTreeSet`** throughout `cavity.rs`, `Interner::ids`, `recount_contacts`. `HashMap` would be the reflex and would be wrong. The one `HashMap` (`FoldCache::map`) is never iterated and is documented as such.
- **Const-generic `D`** rather than a runtime dimension, and `[[u8; D]; 60]` tables. Monomorphising three times is the point.
- **Fixed 240 SMACOF iterations and 20,000 anneal steps** rather than convergence tests.
- **`Stream` not being `Copy`.**

## One place a rule may be applied where it does not belong

`lining_signature`'s two-pass structure carries the comment "Merging them is not an optimisation, it is a semantic change". That is true — but it is a statement about *semantics* (the weights depend on `nearest`), not about accumulation order, and it reads as though it were a §13.1 pin. Someone will later cite it as a determinism constraint on a loop that is really just a data dependency. One word — "semantic change, not an accumulation-order pin" — keeps the two categories apart.
