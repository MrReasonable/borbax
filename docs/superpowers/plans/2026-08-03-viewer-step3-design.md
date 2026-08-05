# Viewer Step 3 — the 3D scene: agreed design

Three-amigos session, 2026-08-03, on branch `viewer-3-scene` (from
`viewer-2-tick`). Supersedes the viewer plan's Step 3 line and its `Structure`
block where they disagree — and they disagree in four places, listed under
*Where the plan is wrong* below.

Everything labelled *measured* was measured on this worktree, today, at the
commit named. Numbers without a measurement behind them are not in this
document.

---

## The measurements everything rests on

Harness `scratchpad/shapeprobe`, path deps on this worktree, release profile,
**seeds 1..=2000**, every row built 2000/2000.

`axis_k/axis_1` is the square root of the k-th covariance eigenvalue of
`Embedding::coords()` over the first. A curve gives `(0, 0)`; a plane gives
`(>0, 0)`; a solid gives both clear of 0.

`d/(r+r)` is `Embedding::distance(i,j)` over `radii[i] + radii[j]`. Below 1 the
two spheres interpenetrate.

### Uniform valence floor — every atom drawn from the same band

| topology | axis2/1 | axis3/1 med (p5) | bonded d/(r+r) med | bonded overlap | rmax/rmin |
|---|---|---|---|---|---|
| chain-4 (val≥2) | **0.000** | **0.000** (0.000) | 1.000 | 22.5% | 1.21 |
| ring-6 (val≥2) | 0.967 | **0.000** (0.000) | 1.117 | 0.0% | 1.36 |
| ring-6 + 2 subs (val≥3) | 0.481 | **0.000** (0.000) | 1.143 | 0.0% | 1.11 |
| star-5 (val≥4) | 0.992 | 0.984 (0.978) | 1.118 | 0.0% | 1.04 |
| branched-8 (val≥4) | 0.574 | 0.570 (0.557) | 1.161 | 0.0% | 1.06 |
| bipyramid-5 (val≥4) | 0.566 | 0.560 (0.556) | 1.029 | 33.3% | 1.04 |

### Mixed valence floor — high-valence centre(s), low-valence leaves

| topology | axis2/1 | axis3/1 med (p5) | bonded d/(r+r) med | bonded overlap | rmax/rmin med (p5..max) |
|---|---|---|---|---|---|
| **star-5 `[4,1,1,1,1]`** | 0.897 | **0.793 (0.734)** | **0.902** | **94.5%** | **2.11** (2.06..2.85) |
| branched-8 `[4,4,1,1,1,1,1,1]` | 0.595 | 0.550 (0.516) | 0.926 | 79.4% | 2.13 (2.08..2.87) |

**Non-bonded pairs never interpenetrate**, in any topology, in any of the 2000
universes: the minimum `d/(r+r)` over 12 000 non-bonded star-5 pairs is
**1.789**, and over 42 000 branched-8 pairs **1.641**. So sphere overlap
carries information — it happens across a bond and nowhere else.

### What those numbers decide

1. **A chain embeds to an exactly straight line.** `axis2/axis1 = 0.000` median
   over 2000 universes, and this is *correct* rather than a defect: a path
   graph's hop metric is exactly 1D-realisable, so SMACOF is right to flatten
   it. It is also useless as the first thing a 3D view ever shows.

2. **Every ring is exactly planar** — `axis3/axis1 = 0.000` in 2000 of 2000,
   for the bare 6-ring *and* the substituted one. Substituents do not lift it
   out of the plane; they only make the plane less round (0.967 → 0.481).

3. **`ElementId::from_index(0)` has valence 0 in every universe**, so the
   literal reading of the plan's "one hard-coded molecule" — chain the first N
   elements — is refused by `add_bond` on its first bond in 100% of universes.
   Selection must be by property. `BondError` stays the authority on legality
   and the viewer never sums a valence budget itself; that would be a second
   encoding of `add_bond`'s, which is the `PeriodicTable::pattern()` failure
   mode the seam check already bans by name.

4. **A uniform valence floor is the wrong selection rule, and this is the
   result no plan text anticipated.** Requiring `valence >= 4` of *every* atom
   draws five elements from one narrow band of the table: they come out the
   same size (rmax/rmin 1.04) and their spheres do not touch (0% overlap). The
   render is five identical balls floating apart — technically rank-3 and
   visually a screensaver.

5. **The mixed floor fixes three problems at once**, which is why it is the
   ruling rather than a preference. One `valence >= 4` centre with `valence >= 1`
   leaves gives solidly rank-3 geometry (`axis3/axis1` p5 = 0.734), a centre
   visibly ~2.1× the radius of its leaves, and 94.5% of bonded pairs
   interpenetrating — so the molecule reads as **one connected body without
   Step 4's cylinders**. That is what settles whether Step 3 is self-sufficient:
   it is, but only under this rule.

### Cost measurements

| | |
|---|---|
| `det_math::{sin,cos}` (libm) | 5.77 ns/call |
| `f64::{sin,cos}` (hardware) | 3.19 ns/call — ratio 1.81× |
| orbit camera trig | 4–6 calls/frame; at 60 fps, **2.076 µs/s**, or 0.0002% of one core |
| `embed` on 6 atoms | 79.07 µs — 0.47% of a 16.67 ms frame |
| `canonicalise` on 6 atoms | 6.42 µs |
| `Universe::generate` | 0.25 ms |
| `det_math::sin` vs platform | differs on **37 313 of 1 000 000** inputs on this machine |
| headless `RenderState` creation | 268 ms cold, 52 ms warm, no window |

---

## Scope

**In.** A 3D scene panel below the periodic table showing one molecule built
from the loaded universe's own elements, drawn as shaded spheres at
`Embedding::radii()`, with a mouse-drag orbit camera and wheel zoom. One flat
material colour. No text inside the viewport.

**Out**, each with the step that owns it: bonds and cylinders (4); colour by
any generated property, and any palette at all (4, with its G6 guard); the
molecule builder, add/remove atoms, live re-embed, `BondError` surfaced with
its reason (5); the signature overlay (6); selection-driven molecules — the
scene does not react to the clicked element; `borbax-geometry`, `depth_sort`,
fixed orthographic/isometric projections (Task 19); binding, folding,
cavities; persisting camera state across a seed change; screenshot export;
any window-title change; **issue #25**, which is a physics change that moves
every golden hash and must not be touched from the viewer.

---

## Rulings

### The molecule

**R1. The demo molecule is a `valence >= 4` centre bonded to four
`valence >= 1` leaves, resolved against the loaded table at reload time.**
Measured rank-3 in 2000/2000 universes with `axis3/axis1` p5 = 0.734. Element
*ids* are universe-dependent, so the molecule is specified by valence and
topology and never by a literal `ElementId`.

**R2. Selection is "lowest-index element meeting this slot's valence floor,
never reusing one".** Lowest-index is a rule, not a search; it takes no
opinion about which element is interesting, which would be the viewer choosing
chemistry.

**R3. Hard-coded *topology* is in scope; hard-coded *chemistry* is not.** The
line: a topology is admissible when it is justified by what it shows about
`embed`'s output — R1 is justified by a rank measurement over 2000 universes —
and is a §5 G3/G6 breach the moment it is justified by resemblance to a real
molecule, **including in a comment**. The shipped G1 guard only sees data
*files*, so a `Mol12` built in code is invisible to every guard and this is a
review obligation rather than a checkable one. Written down because that is
the only enforcement available.

**R4. The scene carries a one-line caption, built in `state.rs`.** Without
one, a scene beside a clicked element reads as "the element you clicked, in
3D", which is false. `panel.rs` is forbidden `format!` by the shipped seam
check, which is exactly what makes a string assertion describe the window.

### Architecture

> **HISTORY — this section predates the Bevy ruling of 2026-08-03 and is not
> what shipped.** `docs/superpowers/plans/2026-08-01-borbax-viewer.md` directs
> readers here as the authority on the panel, so the supersession has to be
> stated *here* rather than only there. What changed:
>
> - **R5** specifies offscreen rendering into our own texture with our own depth
>   attachment, refuses `egui_wgpu::CallbackTrait`, and defines `scene/draw.rs`
>   as a `wgpu` function. Bevy owns the render graph now; none of that exists.
>   The measurement R5 rests on is still true of `egui_wgpu` and is why the
>   `eframe` route was hard, which is part of why the engine changed.
> - **R8** locates the single `f64 -> f32` cast in the `wgpu` layer. It shipped
>   in `scene.rs` as `narrow`, and `xtask`'s
>   `check_the_viewer_narrows_in_one_place` is what holds it there.
>
> The *requirements* in this section — one narrowing site, no engine maths in
> our code, the depth-buffer hazard — all survived the change. Only the
> mechanism named for each is history.

**R5. Render offscreen into our own texture with our own depth attachment, and
paint it with `ui.image`. Do NOT use `egui_wgpu::CallbackTrait`.** This is the
session's largest deviation from the plan and it is forced by a measurement:
`WgpuTestRenderer::render` hardcodes `depth_stencil_attachment: None`, while
the real window gets a depth buffer only from `NativeOptions.depth_buffer` — a
`main.rs`-only setting no harness can replicate. A depth-buffered
`CallbackTrait` is therefore correct in the window and **untestable**, proven
by the validation error it raises under `egui_kittest`:

```text
Incompatible depth-stencil attachment format: the RenderPass uses a texture
with format None but the RenderPipeline uses an attachment with format
Some(Depth32Float)
```

Offscreen rendering makes depth and MSAA ours, makes egui's pass configuration
irrelevant, and turns `scene/draw.rs` into
`(&wgpu::Device, &wgpu::Queue, camera, instances) -> wgpu::TextureView` — a
function testable with a bare `wgpu` device and no `egui` at all. It is
smaller *and* more testable than the plan's shape.

**R6. Depth testing is required and is not an optimisation.** 94.5% of bonded
sphere pairs interpenetrate under R1, so back-to-front sorting is wrong rather
than approximate — and screen-space occlusion happens on every orbit
regardless of whether the spheres intersect in 3D.

**R7. `embed` and `canonicalise` run once, in `state.rs`'s reload path, and
the `Embedding` is cached on `ViewerState`.** `embed` is 240 fixed SMACOF
iterations at 79 µs; per-frame would *work*, which is exactly why it needs
saying. §8.6 keeps `embed` off the per-step path and the viewer's own module
doc is "the viewer computes no physics".

**R8. Exactly one `f64 -> f32` cast site**, in the function that fills the
instance and uniform buffers, carrying one `#[expect(clippy::cast_possible_truncation)]`
with a reason. Camera and embedding arithmetic stay `f64`. `crates/borbax-ui/src/lib.rs`
already licenses `f32` inside the viewer's own rendering, so this needs no new
policy — only a single boundary so it stays reviewable.

**R9. The instance buffer is built by a pure function returning
`Vec<Instance>`.** This is a testability constraint QA imposed and it is cheap:
without it, the affine cross-check (T11) cannot exist.

**R10. The scene calls `widget_info`.** Measured: a paint callback or a bare
`allocate_response` contributes **no node** to the accessibility tree, so the
scene is invisible to a screen reader *and* to every frame test unless the
implementation deliberately says what it is. Without this, test group B does
not exist.

### The rule that matters most, restated for Step 3

**R11. The viewer may *transform* a `borbax-*` value; it may not *derive*
one.** The plan's rule — "no arithmetic on a value that came out of a
`borbax-*` call" — is falsified by Step 3's first line, because every
coordinate from `embed` gets multiplied by a matrix. The surviving version
distinguishes an invertible, information-preserving map into screen space from
a new quantity that could have come from the chemistry and didn't.

*Allowed*: `model · [x,y,z,1]`, `f64 as f32`, recentring on the centroid,
scaling by a camera-chosen scalar.
*Forbidden*: computing a bond length from two coordinates when
`Embedding::distance` exists; computing a sphere radius from mass when
`Embedding::radii` exists; inferring a bond from coordinate proximity;
computing an extent that is really a signature reach when `Embedding::reaches`
exists.

### Presentation

**R12. One flat material colour plus shading that depends only on surface
normal and view direction.** No property→colour map of any kind. The precedent
is exact and already shipped: `legible()` is licensed on "a palette encodes an
element *property* — affinity, mass, valence — … Contrast and text size encode
nothing about any element". An arbitrary per-atom colour keyed on index is
**worse** than monochrome: it looks like a property map, encodes nothing, and
lands a palette ahead of Step 4's G6 guard, which is this repository's named
failure mode. Atoms are distinguished by radius, which is a real generated
property drawn faithfully.

**R13. Fixed-height scene panel below the table; both visible at once.** The
table's scarce resource is width (3531 px worst case) and the scene's is area,
so stacking costs the table nothing while a side-by-side split halves it. A
scene *inside* the table's horizontal `ScrollArea` would also put orbit-drag
and scroll-drag on one gesture, and Step 2 measured that pointer input does not
reach widgets inside a scrollable `ScrollArea` under `egui_kittest`.

**R14. No text inside the 3D viewport at all** — no labels, no per-atom
symbols, no unit words, no axis gizmo text. This keeps Step 2's
`every_unit_word_on_screen_is_one_this_universe_invented` allow-list valid
unchanged, which is the cheapest correct answer available. Stripping
`Span -> f32` at the vertex buffer is not a G4 breach: G4 governs unit *names*
reaching a screen and cross-unit arithmetic not compiling, and a vertex buffer
is neither.

### §13.1

**R15. No exception is granted. Camera trig routes through
`borbax_units::det_math`.** The exception's premise is that compliance has a
cost worth buying an exception for, and it is false: `det_math` already exports
`sin` and `cos`, `borbax-units` is already in the tree via `borbax-universe`,
and compliance costs one manifest line and the spelling `det_math::cos(t)`.
The plan's proposed resolution — "a per-site `#[expect]` plus an `xtask`
ruling" — buys **927 ns per second** in exchange for a policy change to a hard
invariant, and it does not even work: an `#[expect]` silences clippy and not
`xtask`, as the plan's own next clause says. Granting an exception would mean
*building path-exemption machinery into the transcendental scan*, whose own
comment calls an exemption that currently exempts nothing "a hole with a sign
on it". Every guard change this project has made has removed exemption
surface; this would be the first to add some, for convenience.

**Two arguments stronger than cost, and the cost argument is the weakest of
the three.** Recorded in this order because a reader who only remembers "it was
cheap to comply" will grant the exception the day it stops being cheap.

*(a) The camera maths is scheduled to move into a golden-tested crate, and this
repository already treats that as decisive.* Task 19 Step 3 and the viewer plan
both say `borbax-ui` becomes a second backend over `borbax-geometry` — "one
projection, two renderers, no duplicated maths" — and `borbax-render`'s output
is byte-compared across the §13.4 matrix. So an `#[expect]`-ed `.cos()` in
`scene/camera.rs` is code the plan says will be lifted into a crate whose bytes
are compared on three platforms; when it moves, the `#[expect]` either travels
with it into a golden path or is dropped. That is *precisely* the reasoning
already written into `xtask` for why `experiments` is scanned: "it is meant to
be lifted into `borbax-molecule` more or less unchanged, so a direct `.cos()`
written there arrives in the simulation later having never been checked." The
exception asks for the opposite ruling on the same fact pattern.

*(b) §13.1's transcendental rule is already a single-site rule, and the site
already exists.* Wall-clock got `check_wall_clock_has_one_home` because it had
no home and one had to be invented. Trig's home is `det_math`. "Grant an
exception" here means **demoting an existing single-site rule to a two-site
rule with a path exemption**.

**The best argument for the exception, stated fairly.**
`check_the_viewer_stays_a_leaf`'s doc says outright that leaf-ness "licenses
everything else the crate is allowed to do — `f32`, `HashMap`, a wall-clock
read and unordered iteration are all permitted in there because it is a pure
consumer." That is a documented, enforced licensing argument, and it is what
makes R8's `f32` vertex buffer fine. It does not carry for trig for three
reasons: everything it licenses has **no guard to tear a hole in** — and the
one that did, wall-clock, got a single-site rule rather than an exemption; the
camera maths has a documented migration into a golden path and `HashMap` does
not; and compliance costs one manifest line.

**The coupling risk this creates, named so it is refused later.** Routing the
viewer through `det_math` invites pressure to *widen* `det_math` for display
convenience — an `acos_unchecked` that skips the clamp, or f32 wrappers.
Additions to `det_math` stay `f64` and keep their contracts. An f32
transcendental is a different function, not a cheaper one: measured, `f64::sin`
rounded to f32 differs from `f32::sin` on **648 594 of 1 000 000** samples
(64.86%). There is no f32 wrapper in `det_math` and one must not be added.

**R16. Add `det_math::tan`. Do not write `cos(h) / sin(h)`.** This ruling is
the reverse of the one this session first agreed, and the reversal is the
determinism auditor's — it is recorded rather than quietly corrected because
the losing argument was mine and is the one a reader would reach for again.

*The argument that lost:* the projection's only need is `f = cot(fovy/2)`,
`sin h ∈ [0.38, 0.71]` over any plausible fov, so `cos/sin` is well
conditioned and adding a wrapper puts an entry in the determinism chokepoint
that no result-affecting code calls.

*Why it loses.* **The two spellings are different functions.** Measured over
890 000 values of `h` from 0.0001° to 89°: `1/tan(h)` and `cos(h)/sin(h)`
differ in **354 665 of 890 000 f64 results — 39.85%**, worst 3 ulp. After the
`as f32` the GPU receives, **0 of 890 000** differ, so neither choice can
produce a visible pixel. That is exactly what makes it dangerous: it is a pair
a later tidy-up can swap in either direction with no test noticing, in a
codebase whose §13.1 discipline is *"float accumulation order is pinned and
commented as load-bearing; nobody tidies it later"*.

And **both documents that own this decision already answer it.** `det_math`'s
module doc: *"Adding a function here is cheap and is the right response to
needing one… when one of those is genuinely wanted, it gets a wrapper here
rather than an entry removed there."* `clippy.toml`: *"wanting one is the
signal to add it there, not to delete a line here."* Writing `cos/sin`
deliberately to route around an absent wrapper sets the precedent that a
missing wrapper is worked around, and the next one will be `atan2`, spelled
worse.

*Portability, read from the vendored `libm` 0.2.16 source rather than assumed.*
`src/math/arch/mod.rs` exports, in full: wasm32 `{ceil, ceilf, fabs, fabsf,
floor, floorf, rint, rintf, sqrt, sqrtf, trunc, truncf}`; x86+sse2 `{sqrt,
sqrtf, fma, fmaf}`; aarch64+neon `{fma, fmaf, rint, rintf, sqrt, sqrtf}`; plus
i586 `{ceil, floor}` and the `x87_exp*` family under `x86_no_sse`. **`tan` is
in none of them.** `tan.rs` reaches `rem_pio2`, whose only arch-conditional
line is guarded on `x86` without SSE2 — **and `sin.rs` and `cos.rs` share that
same `rem_pio2`.** So `libm::tan` sits in the *identical* portability class to
the already-wrapped `libm::sin` and `libm::cos`: safe on all three §13.4
targets, and a hazard the same day a 32-bit leg is added. The only new audited
surface is `k_tan.rs`, which was read and is clean.

The wrapper's doc must say that `cos/sin` is **not** a behaviour-preserving
substitution, with the 39.85% and the 3 ulp attached, or someone will
"simplify" it back.

**R17. The camera accumulates yaw and pitch from drag deltas** and never
reconstructs an orientation from a direction vector — `det_math` has no
`atan2`, and a camera needing one would be a design that reached for a wrapper
rather than the other way round.

---

## Where the plan is wrong

1. **"Render one hard-coded molecule"** read with "the first N elements" is
   refused in 100% of universes, and the fixed version is a straight line.
   Measured; see R1 and measurement 3.
2. **The trig resolution** (plan lines 209-215) proposes an `#[expect]` that
   silences half the guard, to buy 927 ns/s. Withdrawn — R15.
3. **The plan is silent on depth**, which is the single decision that
   determines whether Step 3 is testable at all. R5/R6.
4. **`mesh.rs` "sphere and cylinder instancing"** — cylinders are bonds, which
   the plan's own step list puts in Step 4. Spheres only.
5. **`palette.rs` is listed in Step 3's structure** while the step text puts
   colour in Step 4. Not created now: an empty file in the strictest tier is a
   file the guard has nothing to say about.
6. **"`eframe` bundles `egui` and `wgpu`"** is true and misleading — it does
   not follow that we may *use* them without direct dependencies, because
   naming `eframe::wgpu` puts the file in `main.rs`'s tier and out of reach of
   every test. `egui-wgpu` and `wgpu` become direct dependencies for that
   reason alone; measured to add **zero** new packages, since both already
   resolve at those versions.
7. **`scene/camera.rs` as "the only place with its own maths"** is false for a
   UV sphere, which needs `sin`/`cos` per vertex ring. It becomes true under
   R18 below.

**R18. The sphere mesh is an icosphere by midpoint-normalise subdivision**,
which needs `sqrt` and no transcendentals at all. `sqrt` is IEEE-754-exact and
explicitly permitted. This collapses Step 3's entire workspace trig surface to
the camera basis. It also rules out using `Geodesic` as the mesh: its faces are
private and dropped after `build()`, and D is physics under sweep (§22.2), so a
mesh keyed on it would make `borbax sweep --resolutions` change how a sphere
looks.

---

## The guards

Four, and three of them close holes that exist on `main` today.

### Guard A — `libm` has one home *(lands with the trig ruling, before any camera code)*

**Measured on this worktree:** `libm::cos(t) + libm::sin(t)` written in
`crates/borbax-ui/src/state.rs`, with `libm = { workspace = true }` added to
that crate's manifest, passes **all six gate legs** — `cargo xtask` prints
"all checks passed" and `cargo clippy -p borbax-ui` is clean. `BANNED_CALLS`
has no `libm::` entry and `clippy.toml` bans the inherent method `f64::sin`,
not every function named `sin`. `clippy.toml`'s own header says so outright:
"`libm::exp` and friends are untouched by these entries".

This is not merely a doctrine gap. **`det_math::acos` clamps its argument**,
and its doc records that as load-bearing with a measurement — 685 of 3600
compositions of `rotation_matrices()` produce `(tr(R)-1)/2 > 1` and would
return NaN unclamped, and a runtime NaN's sign bit differs between aarch64 and
x86-64. An orbit camera extracting a pitch angle from a dot product is exactly
where someone writes `acos`. So the likely wrong first draft of Step 3
bypasses the one wrapper whose clamp was built for this, silently.

**Classify it correctly, because overclaiming would get it deleted.** This is a
**scope** guard, not a divergence guard. `libm::sin` written in
`borbax-molecule` is bit-identical to `det_math::sin`; nothing diverges. Three
other things break:

*(a) The libm-bump review loses its scope.* CLAUDE.md commits to Renovate
TIER B delivering libm bumps as `determinism-review`-labelled PRs, and that
review's procedure is to read `det_math.rs` and check each function's
arch-dispatch status against the new version. **That procedure is correct only
if `det_math.rs` is the complete list of libm functions this workspace calls.**
Today it is — by accident, one manifest line from ending, and nothing in the
diff that ends it would say the review's scope had just become wrong.

*(b) `libm::acos` is unclamped and a camera is exactly where that bites.*
Measured over 500 000 near-parallel unit-vector pairs — each independently
normalised, perturbed by ~1e-9, renormalised, which is the shape of an
orbit-camera angle extraction — **89 403 (17.88%) produce `|dot| > 1`**, worst
overshoot 4.44e-16. Every one is `NaN` from `libm::acos` and correct from
`det_math::acos`. Stated honestly: that is a *correctness* hazard, not a
determinism one — a NaN camera goes black, deterministically.

*(c) The future case is squarely result-affecting.* Task 11+ (folding, cavity)
will want trig in `borbax-molecule`, which has **zero** `det_math` call sites
today and no `libm`. One manifest line plus `libm::atan2(..)` there passes
clippy, passes xtask, is genuinely portable — and is a transcendental in a
result path outside the audited surface with no wrapper contract.

**The template carries a flaw that would be fatal if copied, and the polarity
is what makes it fatal here.** `check_wall_clock_has_one_home` only ever
inspects files that are *not* the sanctioned home, so it **cannot distinguish
one home from zero**. For wall-clock that degenerates harmlessly — zero clock
reads is fine. **For `libm` the meaning inverts:** zero `libm::` inside
`borbax-units` means `det_math` has stopped calling libm, which is exactly the
state this guard exists to make impossible, and the copied template would be
green for it. So the source half needs an explicit **positive** assertion that
its sanctioned site is still populated — `check_blocklist_present`'s "Loud, not
`Ok(())`" lesson, which `check_wall_clock_has_one_home` does not implement.

| Half | Mutation | Must |
|---|---|---|
| manifest | `libm = { workspace = true }` in `crates/borbax-ui/Cargo.toml` | FAIL |
| manifest | `mathlib = { package = "libm", workspace = true }` in `crates/borbax-molecule/Cargo.toml` | FAIL — the rename hole |
| manifest | `libm` **removed** from `crates/borbax-units/Cargo.toml` | FAIL — the "at most one" degeneration |
| manifest | `libm` as a `[dev-dependencies]` entry elsewhere | FAIL, or the rule is "runtime only" and the message says so |
| source | `libm::acos(x)` in `crates/borbax-units/src/lib.rs` | FAIL |
| source | `det_math.rs` renamed or deleted | FAIL — file-missing is loud |
| source | every `libm::` removed from `det_math.rs` | FAIL, from the positive assertion |

- **Bookkeeping:** it counts *manifests whose dependency tables name `libm`*,
  over comment-stripped text so the several prose mentions do not count. That
  is **declarations, not uses**; a transitive re-export is invisible to it and
  the message must say so. The counter must count **manifests, not
  directories**, and increment where a manifest is read — the Task 8 `pairs`
  defect exactly. Assert `manifests_examined >= <n members>`.
- **How it goes silent:** if `walk()` returns no `Cargo.toml` — layout moved,
  root wrong — the loop body never runs and the guard is green. A missing scan
  root must be a **failure**, not `return Ok(())`. It iterates
  `TRANSCENDENTAL_SCAN_ROOTS` so it inherits `check_no_crate_escapes_the_scan`'s
  coverage rather than a parallel constant that can drift.
- **Match the token anywhere in the manifest, not a key**, so
  `mathlib = { package = "libm" }` is caught.
- Airtight cross-crate, verified by the negative: removing the manifest line
  gives `error[E0433]: cannot find module or crate 'libm' in this scope`. A
  transitive dependency is not nameable.
- **Known residual hole, written down rather than fixed:** a member vendoring
  libm source under another crate name, and `[patch]`/`[replace]` entries. Both
  are dependency reviews, the same disposition `BANNED_TYPES` gives third-party
  hashers.
- **Why not a `BANNED_CALLS` entry instead:** a `"libm::"` needle there would
  fire on `det_math.rs` itself and force the path exemption the transcendental
  scan exists without. `leaves_the_libm_free_functions_alone` pins the current
  pass as *deliberate*, so closing the gap needs a separate guard, not a line
  in that list.

### Guard B — `wgpu` is named in code only by the scene files

**Measured hole:** `viewer_banned_imports` returns `["egui", "eframe"]` for the
default tier and `["eframe", "format!"]` for drawing files. **Raw `wgpu`
matches neither.** So `state.rs` — the file whose whole identity is "names no
UI type, so a test can reach it without a window" — may write
`use wgpu::Buffer;` and pass the seam today. `egui_wgpu` *is* caught for
default-tier files because it contains the substring `egui`; that asymmetry is
the trap.

**A second hole found in the same read and closed with it: `epaint` also passes
the strictest tier**, and it is egui's painting layer. `use epaint::PaintCallbackInfo`
would satisfy the guard's letter and defeat its purpose.

- **Mutations:** `use wgpu::Device;` in `state.rs` must fail. A new
  `scene/pipeline.rs` naming `wgpu` without being listed must fail —
  unknown ⇒ strict. `use epaint::…` in any non-drawing file must fail.
- **Bookkeeping:** count `wgpu` as an **identifier token** via `names_type`,
  not `contains`. With `contains` it counts `egui_wgpu` too, makes the two
  rules indistinguishable, and reports a number one higher than the thing it
  names for every drawing file.
- **How it goes silent:** if spelled as a per-file ban list rather than a
  **count** pinned to a named list, a new file with a lax entry reopens it —
  which is what the `eframe_homes` count already exists to prevent.

### Guard C — `embed` has exactly one call site in `borbax-ui/src`

`check_the_viewer_calls_the_chemistry_once` counts `Universe::generate` and
`.pattern(` and **nothing else**, so Step 3's new calls arrive with no guard
having an opinion about them. Extended by the same argument that put
`Universe::generate` in Step 1 rather than Step 7.

- **Mutations:** a second `embed(` in `panel.rs` must fail; **deleting the only
  call site must also fail** — `total != 1`, not `total > 1`. The existing
  `Universe::generate` guard already does this correctly; copy it rather than
  re-deriving it.
- **Bookkeeping:** it counts textual `embed(` occurrences outside
  `#[cfg(test)]`. That is **not** "embeddings computed" —
  `if x { embed(a) } else { embed(b) }` is two sites and one embedding. The
  message must say "call sites" and must not claim to bound work.
- **How it goes silent:** `use borbax_molecule::embed as lay_out;`. The
  `Universe::generate` guard already discloses this exact hole in its own
  failure text; copy the disclosure verbatim.

### Guard D — the `re_embeds` counter

- **Mutation:** move `embed(..)` into `panel::draw`'s body — T20 and T22 must
  fail, T21 must not.
- **Bookkeeping:** increment **adjacent to the `embed` call**, never at
  function entry, and count `embed` and **not** `canonicalise`. Those are
  called together today and will not always be, since `canonicalise` returns
  `Result<_, Capped>` and can fail; a counter over the pair reports
  "embeddings" while measuring "attempts".
- **How it goes silent:** exactly as `regenerations` does — an `embed` written
  straight into `panel.rs` never touches the state method. That is why Guard C
  exists, and Guard C's message must say so.

---

## The test list

Discriminator is the point. Anything without a nameable one is not here.

### Group A — pure logic, no window, no GPU

| # | Name | Fails on |
|---|---|---|
| T1 | `the_camera_basis_is_orthonormal_and_right_handed` | `cross(up, forward)` the other way — a left-handed basis mirrors the whole scene, and **orthonormality alone does not catch it**; the `det == +1` arm does |
| T2 | `a_quarter_turn_is_not_a_full_turn` | degrees-vs-radians. The **negative** arm carries the information: "2π returns to start" alone passes on a camera that ignores the angle entirely |
| T3 | `the_target_projects_to_the_centre_of_the_screen` | transposed matrix; row/column-major; eye and target swapped (`w <= 0`); translate-before-rotate |
| T4 | `right_is_right_and_up_is_up_on_screen` | a y-flip. T3 **structurally cannot see this** — the centre is the fixed point of a flip |
| T5 | `the_depth_range_is_wgpus_and_not_opengls` | any `-1..1` depth matrix. Plausible, renders, silently clips the front half |
| T6 | `the_projection_uses_the_aspect_ratio` | aspect ignored or applied to the wrong axis — atoms are ovoid and nothing else notices |
| T7 | `zooming_in_never_crosses_the_near_plane` | unclamped zoom: the molecule vanishes when you scroll |
| T8 | `the_camera_does_not_degenerate_at_the_poles` | `look_at` with fixed `up = +Y`; **and** the "fix" that guards NaN while letting the view flip |
| T9 | `the_camera_is_a_pure_function_of_its_parameters` | accumulating into the *matrix* rather than the parameter — drifts, and makes every other camera test unreproducible |
| T10 | `every_atom_is_inside_the_viewport_at_the_default_camera` | a default distance tuned to one seed; camera inside the molecule |
| T11 | `the_instances_are_embeds_coordinates_under_one_affine_map` | any per-atom fudge — an atom clamped into view, scaled by its own mass, dropped, or a 2D projection baked in. **This is the test that says the picture *is* the chemistry** |
| T12 | `one_instance_per_atom_and_no_more` | an off-by-one loop — indistinguishable on screen for a symmetric molecule |
| T13 | `sphere_radii_come_from_the_embedding_and_are_all_positive` | `radius = 0.0` — an invisible scene with every test green; radius from `mass` "because it looks better" |
| T14 | `the_demo_molecule_is_buildable_in_every_universe_it_is_shown_in` | hard-coded `ElementId`s whose valence budget refuses the bond in some universes |
| T15 | `there_is_exactly_one_f64_to_f32_cast_site` | the cast spreading, which is how a display choice reaches a chemistry value |

**T11's bookkeeping, called out because it is the Task 8 defect exactly.**
`residuals_checked` must count atoms **beyond** the four used to fit the affine
map — the fit points have zero residual by construction, so a counter over all
atoms reports a bar four higher than the evidence. Assert
`residuals_checked >= 4`, and assert the four fixture points are
**non-coplanar** before fitting, or the map is underdetermined and the test is
vacuous. The tolerance is derived, not chosen after seeing the residuals: an
f64→f32 round trip is ~6e-8 relative.

### Group B — frame tests through `egui_kittest`, CPU only

| # | Name | Fails on |
|---|---|---|
| T16 | `the_scene_has_an_accessible_identity` | the natural first implementation — measured, a bare `allocate_response` plus callback yields **no node at all** |
| T17 | `the_scene_occupies_the_space_it_claims_at_the_size_the_window_opens_at` | a scene allocated `Vec2::ZERO` or 1×1 — green everywhere else. Direct heir of Step 2's `every_cell_is_on_screen…` |
| T18 | `dragging_the_scene_orbits_the_camera_and_a_click_does_not` | drag wired to the wrong axis; sign inverted; drag read from `clicked()` so it jumps |
| T19 | `the_wheel_zooms_and_only_over_the_scene` | zoom read from global scroll state — scrolling the periodic table zooms the molecule |
| T20 | `an_idle_scene_does_not_re_embed` | `embed` in the paint body |
| T21 | `changing_the_seed_re_embeds_exactly_once` | a cached embedding that never invalidates — the previous universe's molecule beside the new universe's count |
| T22 | `orbiting_does_not_re_embed` | re-embedding "because the camera moved" — the §8.6 regression in the shape it actually arrives in |
| T23 | `the_scene_costs_no_second_layout_pass_when_nothing_changed` | an ungated `request_discard`, which Step 1b shipped and Step 2 pre-empted — it self-repairs the frame-order guards |

Input reaches a bare allocated rect, measured and **unlike** Step 2's
`ScrollArea` situation: pointer down at (100,100) through two moves accumulates
a `drag_delta` of exactly `[60.0, 50.0]` over two frames, and a line-wheel
event delivers `smooth_scroll_delta.y = 108` while hovered.

### Group C — GPU tests behind a feature, adopted deliberately

| # | Name | Fails on |
|---|---|---|
| T24 | `the_scene_paints_something_other_than_the_background` | radius 0; zero instances; black-on-black; camera inside the molecule; everything behind near; a pipeline that never binds |
| T25 | `the_molecule_is_where_the_camera_says_it_is` | a pass that paints outside its allotted rect, over the periodic table |
| T26 | `orbiting_changes_the_picture` | a camera uniform that is never uploaded — matrix correct, buffer stale, **every Group A test passes** |

**T24's bookkeeping:** it must count pixels **inside the scene rect**, not the
whole image. Counting the whole image lets the periodic table's text satisfy
the bar, and the guard then reports "the scene drew something" while measuring
"the panel drew something". Its threshold is measured across ≥20 seeds and set
below the minimum, with both numbers written down.

**Camera claims evaluated and rejected**, recorded so they are not proposed
again: *"the camera position lies on the sphere of radius r"* passes on any
swap of the two angles, either sign flip, and degrees-vs-radians — it
discriminates only "the radius is used at all". *"A full 2π orbit returns to
the start"* passes on a camera that does nothing with the angle. *"Screen size
is monotonic in zoom"* catches a sign error and nothing else; its real value is
the inverse arm, which is T7.

---

## The GPU test leg — a decision with a CI cost, taken here

**Adopted**, as a `wgpu` dev-feature on `egui_kittest`, in its own commit after
the CPU-side scene works.

It buys exactly one thing nothing else can buy: T24. "The scene is invisible
and every test is green" is this project's recorded failure mode **twice over**
— Step 1's empty opening state and Step 2's unreadable grid, then Step 2's own
fix making the cells invisible. Every one of those was caught by a person
looking at the window, and T24 is the first test that could have caught any of
them.

**What it costs, stated rather than discovered at review.** Measured working
headless on macOS/Metal: adapter `Apple M1 Pro / Metal / IntegratedGpu`,
render state in 268 ms, a real custom pipeline reaching the framebuffer
(`red_pixels=5766`). `ubuntu-latest` installs no Mesa Vulkan driver and
`create_render_state` **panics** rather than erroring, so a missing adapter
arrives as a panic inside a third-party crate and reads like a viewer bug.
Windows/WARP is untested.

**So the leg runs where an adapter exists, and says so loudly where one does
not.** A test that silently passes when it could not run is the guard-goes-
silent failure mode; a skip must be loud or it must not exist.

**A GPU snapshot is never a §13.4 golden and must never enter the determinism
matrix.** Two renders on one machine are byte-identical; Metal, lavapipe and
WARP are not, and asserting otherwise would put a false portability failure in
front of the one gate that must never cry wolf.

---

## Shipped untested, stated rather than hidden

1. **That the window shows a molecule.** Every CI leg proves only that
   `eframe` compiles on its host. `main.rs` is unreached.
2. **That the picture is legible to a person.** Both Step 2 defects lived in
   this gap, both with every test green.
3. **That depth is correct.** Nothing in the list catches a missing depth
   buffer — back atoms painting over front ones passes T24 and T25 both.
4. **That the lighting produces spheres rather than flat discs.**
5. **That the image is identical across platforms.** Deliberately not claimed;
   see the GPU leg above.
6. **That a number means what its label says** — "radius" that is actually
   mass. A reading rule, not a testable one.
7. **That the GPU leg runs on Linux CI.** Measured on macOS/Metal only. Prove
   it on a real CI run before any test depends on it.
8. **That 60 Hz holds** with N atoms. `criterion` if it ever matters; CLAUDE.md
   forbids timing assertions inside `#[test]`.

---

## The manual checklist

Step 3 is a *rendering* step, and this project has twice shipped a panel that
every test liked and no person could read. Run this in front of the window,
not instead of the suite.

1. Open on the default seed — is a molecule visible **without touching
   anything**?
2. Are the atoms round and shaded, or flat discs?
3. Drag left — does it turn the way your hand went? Drag up — same question.
4. Drag past the top — does it flip, jitter, or stop cleanly?
5. Scroll in until it stops: does the molecule pass through the camera, invert
   or vanish? Scroll out until it stops — still visible?
6. Resize narrow, then wide. Still round?
7. Type a different seed. Does the molecule change **in the same frame** as the
   element count?
8. Type a refused seed. Does the previous molecule stay on screen beside the
   refusal? It should not — that is `Outcome`'s argument, in 3D.
9. Turn it slowly and watch for atoms popping through ones that should occlude
   them.
10. Is any part of the scene painted over the periodic table, or the table over
    the scene?
11. Run it with the system theme set to light — Step 2 pinned `Visuals::dark()`
    for exactly this; check the scene did not reintroduce the dependency.
12. Show it to a child and say nothing. If she does not ask what it is, it is
    not showing anything.

---

## Carried out of the session

- **A disconnected 2-atom molecule embeds bit-identically to a bonded one.**
  `hop_distances` fills unreachable pairs with `diameter + 1`, which for two
  lone atoms is 1 — the same as a bond. Invisible at Step 3, which draws one
  connected molecule; **at Step 5 the builder will draw two unconnected atoms
  touching**, and that is where it must be handled. Routed to Step 5.
- **`n <= 1` short-circuits to the origin** with zero extent, so the camera must
  never divide by `radius_of_gyration`.
- **`palette.rs` in `borbax-ui` conflicts with §14.5/§18.3**, which put
  `palette` in `borbax-geometry` under "one projection, two renderers, no
  duplicated maths". A Step 4 conflict, flagged now because Step 3 sets the
  precedent.
- **The substring seam matcher cannot grant `egui_wgpu` without granting full
  `egui`.** Not fixed: the guard's two load-bearing halves are untouched and a
  word-boundary matcher would be new machinery guarding a distinction nobody
  has abused. Recorded so the next person knows it was seen.
- **Do not mutate the worktree while review agents are reading it.** One
  reviewer reported a probe as shipped code because it read the tree during the
  window when the probe was planted. The finding was independently real; the
  location was not.
- **`xtask`'s root is `env!("CARGO_MANIFEST_DIR")`, baked at compile time.** A
  cached `xtask` binary run inside a *copied* tree scans the **original** path
  and reports "all checks passed" for a mutation planted in the copy. One
  reviewer hit this and nearly recorded a false negative. `touch
  xtask/src/main.rs` before believing any xtask result taken in a copy.
- **`DATA_FREE_ROOTS` and `TRANSCENDENTAL_SCAN_ROOTS` are equal today**
  (`["crates", "experiments"]`) and nothing keeps them equal. Not fixed here;
  noted because two guards now depend on the coincidence.
- **The `1/tan` vs `cos/sin` conditioning near the ends of the fov range** was
  bounded at 3 ulp but not analysed. Routed to `geometry-numerics-reviewer` at
  the Step 3 review.
