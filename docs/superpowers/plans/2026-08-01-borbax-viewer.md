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

**Native, per spec §14.0**, whose three reasons all still hold — the geometry is
3D so a browser needs WebGL regardless and gains nothing; `wgpu` is already in
the dependency tree for §13.5's GPU screening; and a web build adds a second
language, a second toolchain and a wire protocol for a program on the same
machine. The fourth reason is the one that decides it here:

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

The discriminator: **every number on screen must be traceable to a call into a
`borbax-*` crate.** If a value cannot be, it does not belong on screen. A
reviewer should be able to grep the viewer for arithmetic and find only
projection, layout and colour.

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
  names. `xtask`'s G2 blocklist scan must cover the viewer crate — a hard-coded
  "Carbon" in a tooltip is exactly as much a breach as one in a data file.
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

New crate `crates/borbax-viewer` (binary). Depends on `borbax-universe`,
`borbax-molecule`, and later crates as they land. **Nothing depends on it** — it
sits at the top of the crate order so it can never affect a result.

```
src/
  main.rs        eframe entry, window, the seed control
  state.rs       what is selected; nothing computed here
  periodic.rs    the element table panel
  builder.rs     assembling a molecule from elements and bonds
  scene/
    camera.rs    orbit camera — the only place with its own maths
    mesh.rs      sphere and cylinder instancing
    draw.rs      wgpu pipeline
  palette.rs     generated-property colour maps (never real-chemistry ones)
```

**Dependencies.** `eframe` (bundles `egui` and `wgpu`) is the only substantial
addition. It is a **binary-only, non-result-affecting** dependency, so §18's bar
is the low one — it cannot touch simulation output. `wgpu` is already required
by §13.5.

---

## Steps

- [ ] **1. The crate, the window, and the seed.** `borbax-viewer` opens a window,
      takes a seed, calls `Universe::generate`, and shows the element count. Proves
      the whole chain works before any rendering exists.
- [ ] **2. The periodic table panel.** Elements by period and group, selectable,
      with generated names and properties. First thing worth showing anyone.
- [ ] **3. The 3D scene.** `wgpu` pipeline, orbit camera, instanced spheres. Render
      one hard-coded molecule from `embed`'s coordinates.
- [ ] **4. Bonds, and colour by generated property.** Cylinders between bonded
      atoms; palettes driven by affinity, mass and valence.
- [ ] **5. The molecule builder.** Add atoms, add bonds, re-embed live. Invalid
      bonds refused with the reason from `BondError`.
- [ ] **6. The signature overlay.** The 42 sampled directions drawn as surface
      extents, toggleable.
- [ ] **7. `xtask` covers the viewer** for G2 names, G4 units and real-chemistry
      palettes.

Each step ends with something on screen that is better than the step before. If a
step cannot be demonstrated by looking at it, it is the wrong step.

---

## What "done" means for the first version

Ian can double-click it, type a seed, see a periodic table nobody has ever seen,
click through to a molecule, and spin it around in 3D — with his daughter, on a
laptop, without a terminal.

Everything after that is the chemistry arriving in a window that already works.
