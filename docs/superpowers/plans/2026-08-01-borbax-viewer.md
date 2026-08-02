# Borbax Viewer — plan

**Status:** specced 2026-08-01, not started. Begins after Task 10 (binding kernel).

**Why this exists and why it is not Task 19.** The V0 plan renders to SVG at
Task 19 — static pictures, headless, golden-tested, near the end. That is the
right tool for proving the geometry is correct and the wrong tool for the
project's actual purpose. Ian, 2026-08-01, on what Borbax is for:

> a realistic simulation of an imaginary universe that I can show my child in 3D
> — how atoms fuse to make elements, how elements bond to form molecules, how
> molecules grow and fold and create self-replicating molecules, how abiogenesis
> *could* have happened and, later, how entire cells work, replicate and apply
> natural selection + random mutation.

and, on where it started:

> I want to show how organic chemistry works, how molecular machinery operates
> to transport and transform things.

Deferring anything visible until Task 19 means ten more tasks with nothing to
look at, in a project whose point is looking at it. This brings the viewer
forward and lets the chemistry land *into* it.

**It does not replace Task 19.** SVG goldens stay: they are how §13.6 proves the
geometry is identical across platforms, and a screenshot cannot do that. Task 19
shrinks to the golden-testable subset once the viewer owns presentation.

---

## What it is

A native desktop application. `eframe`/`egui` for the interface, `wgpu` for the
3D scene, Rust throughout.

**Native, per spec §14.0**, whose reasons still hold — the geometry is 3D so a
browser needs WebGL regardless and gains nothing, and a web build adds a second
language, a second toolchain and a wire protocol for a program on the same
machine. The reason that decides it here is the fourth:

> a double-clickable application is meaningfully easier for a child than "start a
> server, then open localhost"

**The browser door stays open.** `egui` compiles to WebAssembly from the same
source. The simulation is the part that would suffer there — WASM compute is
slower than native and awkward to thread — and the *rendering* would be fine,
since WebGL and WebGPU are hardware-accelerated. If a shareable link is ever
wanted, the viewer builds for it without a rewrite.

**It grows rather than being replaced.** This is the same `egui` + `wgpu`
foundation §14 specifies for V1, so folding, docking, reactions and the
Chronicle land in this window as those tasks complete.

---

## The rule that matters most

**The viewer computes no physics. Ever.**

It calls `Universe::generate`, `canonicalise`, `embed`, `signature` and — once
they exist — `fold`, `affinity`, `cavities`, and displays what comes back. It
must never reimplement a calculation "just for display", because the moment it
does, the picture and the simulation can disagree and the picture is the thing
being trusted.

The discriminator: **every simulation-derived number on screen must be traceable
to a call into a `borbax-*` crate.** If such a value cannot be, it does not
belong on screen.

"Simulation-derived" is the load-bearing qualifier and an earlier version omitted
it, which made the rule false the moment Step 1 shipped: the seed in the box is a
number on screen and it is **user input**. Three categories are exempt and they
are exempt because they are not results — user input, the seed's own identity,
and layout values.

An earlier version added "a reviewer should be able to grep the viewer for
arithmetic and find only projection, layout and colour", and Step 1 falsified it
immediately — the crate has a `wrapping_mul`/`wrapping_add` pair (a clock
reading) and a `saturating_add` (a counter), none of which reaches the screen.
The rule that survives contact is narrower and is the one that matters: **no
arithmetic on a value that came out of a `borbax-*` call.** That is the viewer
recomputing physics; the rest is bookkeeping.

Corollary worth stating because it reads like a violation and is not: the viewer
*may* use `f32`, `HashMap`, wall-clock time and unordered iteration **inside its
own rendering**. §13.1's bans exist to protect simulation results, and the
viewer produces none — it is a pure consumer. Camera angles and frame times are
not results. What it may **not** do is feed anything back into simulation state,
or let a display choice reach a `borbax-*` call.

---

## Fiction guarantees

Everything in §5 applies unchanged, and the viewer is a new surface for
breaching them, because a UI wants labels.

- **No real element names, symbols, or data.** The viewer displays *generated*
  names. A hard-coded "Carbon" in a tooltip is exactly as much a breach as one
  in a data file.

  **There is no `xtask` G2 blocklist scan to extend, for the viewer or for any
  other crate** — an earlier version of this line said there was. Read from the
  source at Step 1: `check_blocklist_present` asserts only that the *identifier*
  `REAL_ELEMENT_SYMBOLS` still appears in `borbax-universe/src/naming.rs`, and
  the blocklist itself is a **generation-time filter** consumed by
  `is_real(symbol, name)`. No file in this repository is scanned for real
  element names, and no check exists for real-world unit strings either. Step 7
  therefore **writes three checks**; it does not add a path to an existing one.
  It also carries a structural question nobody has costed: the vocabulary lives
  in `borbax-universe`, and `xtask` depending on a chemistry crate to enforce
  §5 inverts the gate — so the blocklist is either extracted or duplicated in
  `xtask`, which is what `FORMAT_SEGMENTS` already does for G5.
- **No real-world units.** Labels read "spans", "quanta", "thermals",
  "world-years". A tooltip saying "Å" or "kJ/mol" is a G4 breach.
- **No real-world visual conventions.** CPK colouring (carbon black, oxygen
  red, nitrogen blue) is real chemistry smuggled in through a palette. Colour
  by *generated* properties — affinity, mass, valence, period — so the palette
  says something true about this universe rather than importing another one.
- The window title, the about box and any exported image must not imply the
  simulation models reality.

---

## What it shows, in the order the chemistry allows

### Now — buildable today, no new physics

**1. The periodic table of a generated universe.** Enter a seed, get that
universe's elements laid out by period and group. Each cell shows the generated
symbol and name. Select one for its properties: mass, valence, radius, affinity,
the shell that is filling. Change the seed and the whole chemistry changes.

*This is the first act — how atoms fuse to make elements — and it is the single
most striking thing the project can show with what already exists.*

**2. Molecules in 3D.** Build one by picking elements and bonding them, or
browse generated ones. Rendered from `embed`'s real coordinates: atoms as
spheres at their true radii, coloured by surface character, bonds as sticks.
Rotate, zoom, pan.

*Not a diagram — the actual output of the stress-majorization solver that the
chemistry runs on.*

**3. The molecule's surface, as binding sees it.** The signature drawn over the
molecule: how far the surface reaches in each of the 42 sampled directions and
what character faces that way. This is the part that makes "shape complementarity"
visible rather than abstract.

### As tasks land

| Lands with | Shows |
|---|---|
| Task 10 — binding | Two molecules docked at their best-scoring orientation, with the score |
| Task 11/12 — folding | A chain collapsing onto the FCC lattice; buried versus exposed atoms |
| Task 13 — cavities | Pockets highlighted on a folded chain |
| Task 14–16 — reactions, catalysis | A reaction stepping through; a cavity holding two reactants adjacent |
| Task 18 — the beaker | Species populations over time; what is actually in the pot |
| Task 17/20b — RAF, activity | Autocatalytic sets; novelty against the neutral shadow |

**The end state Ian named** — molecular machinery transporting and transforming
things — is the far end of this table. It needs folding, cavities, catalysis and
the beaker. The viewer is what makes each step visible on arrival instead of
ten tasks later.

---

## Structure

New crate `crates/borbax-ui` (library plus a `borbax` binary). Depends on
`borbax-universe`, `borbax-molecule`, and later crates as they land. **Nothing
depends on it** — it sits at the top of the crate order so it can never affect a
result.

**`borbax-ui`, not `borbax-viewer`** — corrected at Step 1, which was the only
moment it was free. Three reviewed documents already name this crate
`borbax-ui`: the PRD's §18.3 repository layout, the V0 plan's File Structure,
and Task 19 Step 3, which is written against `borbax-render` and `borbax-ui` as
**a pair of backends over one geometry** — "one projection, two renderers, no
duplicated maths". `borbax-viewer` appeared only here, in a document written in
one sitting that has never been through `/review-plan`, and it breaks that
pairing. The **binary** is `borbax`, because the double-clickable artefact is
what a child finds in a folder and `-ui` names a fact about the repository.

```text
src/
  main.rs        eframe entry and the window — the ONLY file that names eframe
  state.rs       every decision, and no UI import at all
  panel.rs       the drawing — imports egui, contains no `format!`
  periodic.rs    the element table panel
  builder.rs     assembling a molecule from elements and bonds
  scene/
    camera.rs    orbit camera — the only place with its own maths
    mesh.rs      sphere and cylinder instancing
    draw.rs      wgpu pipeline
  palette.rs     generated-property colour maps (never real-chemistry ones)
```

**The seam is per-file and greppable, which is the point.** The original note
against `state.rs` read "what is selected; nothing computed here", and Step 1
found that exactly backwards: `state.rs` is where *everything* is computed,
because it is the half a test can reach without a window. The rule that actually
holds is about imports — `state.rs` names no UI type, `panel.rs` names `egui`
(a CPU layout library that cannot open a window, so `egui_kittest` can drive it
headlessly in CI), and only `main.rs` names `eframe`. Formatting lives in
`state.rs` and `panel.rs` has no `format!`, so a test asserting on a string is
asserting on what the window paints.

**A `camera.rs` writing `.cos()` will fail the gate, and Step 3 should not be
where that is discovered.** `clippy.toml`'s §13.1 ban resolves for every
workspace member, and `xtask`'s textual scan has **no path exemption by
design** — so a clippy `#[expect]` silences only half of it. The viewer produces
no results, so routing camera trig through `det_math` is cost with no benefit;
the likely answer is a per-site `#[expect]` plus an `xtask` ruling. Settle it
before Step 3 rather than during it.

**Dependencies.** `eframe` (bundles `egui` and `wgpu`) is the only substantial
addition. It is a **binary-only, non-result-affecting** dependency, so §18's bar
is the low one — it cannot touch simulation output.

**Its cost is new and wholly the viewer's, and an earlier version of this plan
said otherwise twice.** It claimed `wgpu` "is already in the dependency tree for
§13.5's GPU screening" and "is already required by §13.5". Measured against the
tree at Step 1 rather than assumed: `wgpu`, `eframe` and `egui` appeared in
**zero** manifests and **zero** `Cargo.lock` entries, the whole workspace had
**13** lockfile packages, and §21 lists GPU screening as an explicit V1
non-goal. The PRD's §14.0 is careful — it says the tree needs `wgpu` "either
way", which is an argument about the *eventual* tree; this plan flattened
"eventually" into "already". The native decision stands on its other reasons.

The correction matters beyond tidiness, for the reason CLAUDE.md already records
about the `glam`/`nalgebra` episode: a false justification gets reused later to
wave through a dependency that has not earned it.

What it actually costs, measured on the tree as shipped: **366 packages in
`Cargo.lock`**, against 13 before. `default-features = false` is what makes that
number not 428 — it drops `accesskit → atspi → zbus` and `sctk-adwaita →
tiny-skia`, taking the resolved Linux graph from 261 to 189 and macOS from 167
to 164. ~888 s of new cold-build CPU, and a build-cache archive going from 30 MB
to several hundred per platform. Plus **31** `clippy::multiple_crate_versions`
errors on day one, held by **one** crate-level `#![expect]`, in `lib.rs`.

**Three of those five figures were wrong in the first version of this
paragraph** — 418 packages (the pre-`default-features` resolve, quoted as the
shipped cost in the same sentence that credited the trim), 34 errors, and
"a `#![expect]` in each of the two targets". The last is the worst of the three:
`lib.rs` records that exact wording as *refuted by measurement*, and it was left
standing here. Re-measure rather than carry forward; that is the whole reason
this paragraph exists.

---

## Steps

- [x] **1. The crate, the window, and the seed.** `borbax-ui` opens a window,
      takes a seed, calls `Universe::generate`, and shows the element count. Proves
      the whole chain works before any rendering exists.

      *Shipped.* The window opens on seed 1 and the line reads
      `80 elements · physics v1` — the version is read
      from the universe rather than typed, because §6 says a universe is
      `(seed, physics)` and a viewer that hardcodes the version becomes a liar
      the day `V2` ships. Seeds are decimal `u64`; hex is refused deliberately
      (`7f3a` looks like a seed because §6 spells addresses `U-7F3A21C9@1`, and
      accepting it would settle an open §6 question in the top crate of the
      workspace). A "surprise me" button draws a six-digit seed through
      `borbax_rng` and **types it into the box**, so a universe a child likes is
      one she can write down and come back to — §6 makes the seed the shareable
      thing, and a surprise you cannot record takes that away.

      24 tests, of which four drive the real panel headlessly through the
      accessibility tree — the only way to catch an immediate-mode panel
      painting the count one frame behind the seed that produced it.

      Four `xtask` guards landed with it, two closing holes that predate the
      viewer and two closing claims this crate makes about itself: G1's
      data-file check was a hand-kept list of four crate paths and is now every
      crate under `crates/`, case-insensitively (unknown ⇒ covered);
      `[lints] workspace = true` — the trap that has silently disabled
      `unsafe_code = "forbid"` in this repository once already — was enforced by
      **nothing**; and `borbax-ui`'s leaf-ness and its per-file seam, both of
      which had been asserted in four documents and checked in none.

      **The review found two of those guards did not guard.** The lints check
      read `contains("[lints]") && contains("workspace = true")` — and every
      manifest satisfies the second clause via `edition.workspace = true`, so it
      collapsed to "is the header present", which a comment satisfies and which
      an empty `[lints]` table satisfies while the crate inherits nothing. Three
      reviewers found it independently, each by planting the mutation the
      original probe had not tried: the probe deleted the whole table, which it
      caught; nobody had deleted one line of it. The G1 widening was
      case-sensitive, so `elements.JSON` passed — and on macOS and Windows that
      rename is invisible to the person making it.

- [ ] **1b. A name instead of a number.** A child types "Emily" and gets a
      universe. `seed_from_phrase(&str) -> u64` in `borbax-universe`, **not** in
      the viewer: hashing a phrase in the display layer is a display choice
      reaching a `borbax-*` call, and the obvious spelling (`DefaultHasher`,
      `RandomState`) is a flat §13.1 violation. Returns the *seed*, not a
      `Universe`, so the viewer can show the phrase **and** the number it
      produced — §6 requires the seed be the shareable thing, and a child who
      cannot write down what she got has lost the point of it.

      Needs a golden with literal expected values, not a round-trip: the mapping
      is result-affecting and changing it silently reassigns every universe
      anyone has named. Normalisation (`"Emily"` / `"emily"` / `" Emily "`) is
      pinned with it. The phrase is **not** filtered through the G2 blocklist —
      G2 constrains *generated* names, and a user typing "Carbon" breaches
      nothing.

      *Sequenced here rather than inside Step 1 because it changes a chemistry
      crate and needs a chemistry-crate review, and rather than after Step 2
      because the periodic table should arrive with a seed control a child can
      use.*

- [ ] **2. The periodic table panel.** Elements by period and group, selectable,
      with generated names and properties. First thing worth showing anyone.
      **The first step with real G2 and G4 exposure** — generated names and
      unit-bearing properties both reach the screen here, so the name and unit
      checks from Step 7 land with this step, not five steps later.
- [ ] **3. The 3D scene.** `wgpu` pipeline, orbit camera, instanced spheres. Render
      one hard-coded molecule from `embed`'s coordinates.
- [ ] **4. Bonds, and colour by generated property.** Cylinders between bonded
      atoms; palettes driven by affinity, mass and valence.
- [ ] **5. The molecule builder.** Add atoms, add bonds, re-embed live. Invalid
      bonds refused with the reason from `BondError`.
- [ ] **6. The signature overlay.** The 42 sampled directions drawn as surface
      extents, toggleable.
- [ ] **7. `xtask` gains the three §5 checks that do not exist yet** — real
      element names, real-world unit strings, and real-chemistry palettes.
      **Not "covers the viewer": these checks exist for no crate.** See the
      Fiction-guarantees section. The name and unit checks are pulled forward to
      Step 2, which is where the surface they guard first appears; a guard
      arriving after its surface has shipped is this repository's recorded
      failure mode. The palette check belongs with Step 4.

      **`Universe::generate`'s single call site landed in Step 1, not here.**
      The `regenerations` counter catches `reload` migrating into the paint
      body; it cannot see a `Universe::generate` written directly into
      `panel::draw`, which never touches `reload`. Both were mutation-probed
      before commit — a second call site, a direct call from the paint body, and
      an unconditional `reload` each fail exactly the guard named for them.

      **What the shipped check does NOT cover, so Step 7 can decide whether to
      close it.** It matches the literal spelling over comment-stripped,
      pre-`#[cfg(test)]` source. So `use borbax_universe::Universe as U;
      U::generate(..)`, a helper function wrapping the call, and a macro
      expansion all pass. Closing that means symbol-aware parsing — `xtask`
      already carries `syn` and a crate-wide walker for §13.1, so the cost is
      real but not new machinery. The judgement to make is whether a textual
      check that catches the spelling anyone actually writes is worth more than
      an AST check nobody maintains; the same question was asked and answered
      the other way for G5's format scan, which a review defeated completely at
      the string level.

Each step ends with something on screen that is better than the step before. If a
step cannot be demonstrated by looking at it, it is the wrong step.

---

## What "done" means for the first version

Ian can double-click it, type a seed, see a periodic table nobody has ever seen,
click through to a molecule, and spin it around in 3D — with his daughter, on a
laptop, without a terminal.

Everything after that is the chemistry arriving in a window that already works.
