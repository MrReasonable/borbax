# Borbax Viewer — plan

**Status:** specced 2026-08-01, not started. Begins after Task 10 (binding kernel).

**Why this exists and why it is not Task 19.** The V0 plan renders to SVG at
Task 19 — static pictures, headless, golden-tested, near the end. That is the
right tool for proving the geometry is correct and the wrong tool for the
project's actual purpose. Ian, 2026-08-01, on what Borbax is for:

> a realistic simulation of an imaginary universe that I can show a child in 3D
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

A native desktop application. **Bevy** for the window, the 3D scene and input,
with `egui` (through `bevy_egui`) for the interface, Rust throughout.

**This paragraph said `eframe`/`egui` + hand-rolled `wgpu` until 2026-08-04.**
Ian's ruling on 2026-08-03 moved the viewer to a game engine — "the last thing I
want to have to focus on is the rendering, collisions etc" — under the standing
constraint that Bevy is for **displaying and interacting only**: the physics, the
chemistry and the generation stay ours. Everything below that still says
`eframe`, `wgpu` or `camera.rs` is history and is marked as such where it
matters.

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
  main.rs        one statement, and no licence to name anything
  app.rs         the engine, the plugins, the cameras — names bevy and egui
  scene.rs       the 3D entities and the input — names bevy and egui
  state.rs       every decision, and no UI import at all
  molecule.rs    the demo molecule, and the crate's only door to the chemistry
  orbit.rs       the camera's arithmetic — no engine type, testable with no app
  panel.rs       the drawing — imports egui, contains no `format!`

  -- HISTORY, superseded by the engine change; see the note below --
  builder.rs     assembling a molecule from elements and bonds
  scene/
    camera.rs    orbit camera — the only place with its own maths
    mesh.rs      sphere and cylinder instancing
    draw.rs      wgpu pipeline
  palette.rs     generated-property colour maps (never real-chemistry ones)
```

**The five entries below the marker are intent, not layout**, and line 41's
promise that anything still naming `wgpu` or `camera.rs` "is history and is
marked as such where it matters" is what this marker discharges. Specifically:

- `scene/camera.rs` and `orbit.rs` both claimed the camera. `orbit.rs` is what
  shipped; the 2026-08-03 design memo records that `scene/camera.rs`'s claim to
  be "the only place with its own maths" was false, and R18 collapses the trig
  surface to the camera basis.
- `scene/draw.rs` is the hand-rolled `wgpu` pipeline the Bevy ruling deleted.
  `xtask`'s `VIEWER_GPU_FILES` still names it as an entry a future commit adds,
  so guard and plan agree it is intent rather than layout.
- `palette.rs` is Step 4 work and deliberately not created yet: an empty file in
  the strictest seam tier is a file the guard has nothing to say about.

**The seam is per-file and greppable, which is the point.** The original note
against `state.rs` read "what is selected; nothing computed here", and Step 1
found that exactly backwards: `state.rs` is where *everything* is computed,
because it is the half a test can reach without a window. The rule that actually
holds is about imports — `state.rs`, `molecule.rs` and `orbit.rs` name no UI or
engine type; `panel.rs` names `egui` (a CPU layout library that cannot open a
window, so `egui_kittest` can drive it headlessly in CI); `app.rs` and `scene.rs`
may name the engine and may build no strings; `main.rs` names nothing at all.
Formatting lives in `state.rs`, so a test asserting on a string is asserting on
what the window paints.

**Settled at Step 3, and it went the other way.** This paragraph predicted that
routing the camera's trig through `det_math` would be "cost with no benefit" and
that the answer was "a per-site `#[expect]` plus an `xtask` ruling". It is not:
an `#[expect]` silences clippy and leaves `xtask`'s textual scan firing, so the
proposed exception would not have worked — and compliance turned out to cost one
import. The camera's trig goes through `det_math`, no exception was granted, and
`check_libm_has_one_home` landed with it.

**Dependencies.** Bevy plus `bevy_egui` are the substantial addition — a
*projected* ~530 packages against `eframe`'s 164, which is a real cost in CI time on three
platforms. They are **binary-only, non-result-affecting**, so §18's bar is the
low one: verified rather than asserted at Step 3's review, by checking that the
feature union on `libm` is unchanged, which is the way a leaf consumer *could*
have reached a lower crate.

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

What it actually costs. **Every figure here names the tree it was measured on,
because two of them did not and a reader could not tell which was current.**

- **13** — before the viewer existed.
- **366** — the `eframe` tree, measured when Step 1 shipped. This is the number
  the `default-features = false` note below is about.
- **~530** — the projection for Bevy, above.
- **590** — the Bevy tree as shipped, counted on 2026-08-05 from `Cargo.lock`.
  This is the current number; issue #28 (trimming Bevy's default features) is
  the lever on it and has not been pulled.

On the `eframe` tree: **366 packages in `Cargo.lock`**, against 13 before. `default-features = false` is what makes that
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

- [x] **1b. A name instead of a number.** A child types "Emily" and gets a
      universe. `seed_from_phrase` in `borbax-universe`, **not** in the viewer:
      hashing a phrase in the display layer is a display choice reaching a
      `borbax-*` call, and the obvious spelling (`DefaultHasher`, `RandomState`)
      is a flat §13.1 violation. Returns the *seed*, not a `Universe`, so the
      viewer can show the phrase **and** the number it produced — §6 requires the
      seed be the shareable thing, and a child who cannot write down what she got
      has lost the point of it.

      *Shipped.* Typing `emily` writes `15709401653729972761` into the seed box
      and paints `117 elements · physics v1`.

      **No test count here, deliberately.** This block argued that a stale
      denominator is the recorded failure mode and then carried one anyway,
      which a machine reviewer caught by noticing three different totals across
      the branch. The counts that remain in the source are dated measurements of
      specific probes, and each says so.

      **The signature is `-> Option<u64>`, not the `-> u64` above.** An
      infallible version has to invent a universe for the empty string, which
      reopens `an_empty_seed_box_does_not_quietly_mean_universe_zero` one crate
      down — a real element count on screen attributed to a name nobody typed.
      That test post-dates this plan step.

      **Two boxes, not one, and that is what kept every shipped test.** A single
      control has to guess whether `"7"` is seed 7 or a child called 7, and
      guessing means deleting `letters_in_the_seed_box_are_refused_not_hashed`,
      which exists specifically to forbid "be friendly, hash whatever was typed".
      **Zero test functions were removed and no assertion was weakened** —
      `parse_seed` keeps all four refusals, and
      `a_digit_in_the_name_box_names_a_universe_rather_than_selecting_one` is the
      test that one box makes unwritable. (The first version of this sentence
      said "all 12 previously shipped `borbax-ui` tests pass unchanged". A
      reviewer measured both halves wrong: `main` has **24** `borbax-ui` tests,
      not 12 — that was the acceptance file alone — and the 4 frame tests *were*
      edited, to address their text box by label rather than by tree position.
      Their assertions are untouched, which is the claim that was meant.)

      **The name box drives the seed box, and two rules follow.** Typing a seed
      by hand *clears* the name, because that name did not produce it — otherwise
      the window shows `Emily` beside a universe Emily did not name. An empty
      name box deliberately does **not** clear the seed, because the window opens
      with the name box empty and a universe on screen, so the opposite rule
      would make the opening state contradict itself.

      The golden holds literal values with its authority stated as **the mint,
      not the test**. Its non-ASCII rows are load-bearing: absorbing `char`s
      instead of bytes fails *only* the first non-ASCII **row**, and swapping
      `to_ascii_lowercase` for `str::to_lowercase` fails *only* the Greek
      **row** — both green on every ASCII row. (Row is the right unit and the
      first draft of the source comment got it wrong:
      `accents_are_part_of_the_name` catches the `to_lowercase` swap
      independently and more cheaply. The Greek rows earn their place as the
      only *golden* evidence, not as the only evidence.)

      `to_ascii_lowercase` is a determinism decision: Rust's Unicode tables move
      between compiler releases, so a Unicode-aware fold would let a bump
      silently reassign every named universe. (An earlier version said the
      tables "moved 16.0.0 → 17.0.0 inside this project's pinned toolchain". A
      reviewer measured `.prototools`: it has one commit and has always read
      `1.97.1`, so the pin has never crossed that boundary. The hazard is
      prospective, which is the whole reason to close it now.)

      The `split_whitespace` residue is closed the same way and it took a
      reviewer to show the first attempt was not enough: a single `U+00A0`
      golden row covered **2 of the 25** codepoints `char::is_whitespace`
      accepts, and Unicode's stability policy does not cover `White_Space` at
      all — it covers `Pattern_White_Space`, and `White_Space` has moved in
      practice (`U+180E`, Unicode 6.3). The whole 25-codepoint set is now
      enumerated and pinned.

      NFC and NFD are different universes, pinned deliberately. The phrase is **not** filtered through the
      G2 blocklist — G2 constrains *generated* names, and a user typing "Carbon"
      breaches nothing.

      **§13.1's hasher ban landed first, in its own commit, before any golden
      literal was minted.** `DefaultHasher`, `RandomState` and `SipHasher`
      passed all six gate legs beforehand: `clippy.toml` had no
      `disallowed-types` key at all and `xtask` had no hasher entry, so the ban
      §13.1 states was held by prose. A golden minted from a `DefaultHasher`
      draft is green forever and moves on the next rustc bump, arriving as a
      golden failure whose obvious diagnosis ("regenerate") is the wrong one.

      **One decision the golden structurally cannot defend**, which is why
      `absorb` is a named function with its own test: `Stream::sub` over
      `Stream::new` is invisible in the output, since a `Stream::new` version
      would have its own self-consistent table with every test green. The hazard
      is that `Stream::new(acc, Hash, 0)` — reachable, because a `&str` may
      contain NUL — is bit-for-bit the "surprise me" button's stream.

- [x] **2. The periodic table panel.** Elements by period and group, selectable,
      with generated names and properties. First thing worth showing anyone.
      **The first step with real G2 and G4 exposure** — generated names and
      unit-bearing properties both reach the screen here, so the name and unit
      checks from Step 7 land with this step, not five steps later.

      *Shipped as PR #24, thirteen commits, 502 tests.* The design and every
      measurement behind it are in
      `docs/superpowers/plans/2026-08-03-viewer-step2-design.md`; read that
      rather than this paragraph before touching the panel.

      **The table is 2–4 rows by up to 74 columns where a real one is 7×18**, so
      it scrolls horizontally and does not reflow — column alignment is where the
      families are. Cells show the **symbol only**, against this plan's "symbol
      and name": element names collide in 16.05% of universes and would make
      `egui_kittest`'s query API panic.

      **Two defects were caught only by looking at the window, and every test
      was green over both.** The grid's default spacing showed eight cells of a
      58-cell row, hiding the one thing the step exists to show; and the panel
      shipped grey-on-black at `Color32::from_gray(140)` until Ian said it was
      hard to read — whereupon the *first* fix made the cells invisible
      (white-on-white) with all 58 tests still passing.

      **Five review rounds, and every finding in all five was in the guards
      rather than the feature.** Four were defects introduced while fixing
      earlier ones. Worth reading before Step 7 widens any of those guards.

      **It surfaced a real physics defect** — every monomer has negative binding
      energy and instability pegged at 1.0, because the strain term charges a
      lone unit for a shell that is not there. **Issue #25**, sequenced after
      Step 3, to land with the `N^(2/3)` and Task 15 items it shares a golden
      regeneration with.
- [x] **3. The 3D scene.** Orbit camera, one sphere per atom at `embed`'s own
      coordinates, sized by `Embedding::radii()`.

      *The `wgpu` pipeline in this line is history* — Bevy owns the renderer, so
      the step is entities, a light and two cameras. The design record is
      `docs/superpowers/plans/2026-08-03-viewer-step3-design.md`, whose
      architecture sections predate the engine change; the molecule, the §8.6
      rules and the manual checklist all still bind.

      **The demo molecule is chosen by valence, never by `ElementId`** — a
      `valence >= 4` centre with four `valence >= 1` leaves, lowest-indexed
      first, no element reused. Slot 0 has valence 0 in every universe measured,
      so the plan's original "chain the first N elements" is refused on its first
      bond 100% of the time. The *mixed* floor is what makes it worth looking at:
      a uniform one gives five same-sized balls that do not touch.

      **Two defects were found by Ian using the window, both with 532 tests
      green.** The drag gate asked egui the narrow question, so the molecule
      turned when dragged anywhere on the interface; and the camera was a
      turntable with a pitch clamp, so it stopped dead when dragged upward. The
      clamp is gone — the camera stores an orientation and turns about its own
      axes, which deletes the pole rather than guarding it.

      **Nothing catches "the molecule is too far away".** Measured and stated
      rather than papered over: the field-of-view test catches too *close* only.
      That gap belongs to the manual checklist.
- [x] **4. Bonds, and colour by generated property.** Cylinders between bonded
      atoms; palettes driven by affinity, mass and valence.

      *Ian asked for both directly, 2026-08-04: "I like the sticks idea, it makes
      it easy to visualise" and "it would be nice if atoms had unique colours".*
      The G6 guard lands with the palette — a colour may encode a **generated**
      property and must never reproduce a real-world element colour scheme.

      *Shipped, and three things in the paragraph above were amended by
      execution rather than inherited.*

      **The guard is G5, not G6.** §5's G6 is about results not transferring and
      names the documentation as its enforcement, so it forbids no artefact and
      gives a reviewer nothing to check a diff against. The artefact at issue —
      a per-element colour table using conventional colours — is *literally*
      G5's "mapping table between Borbax entities and real-world entities", and
      a ramp anchored to make the picture read right against real chemistry is
      G3's forbidden transfer function. `xtask` cites G1, G2, G4 and G5 today
      and has no G6 check; this did not mint one.

      **There are two views, not one, and it is a toggle.** Ian's instruction,
      2026-08-06. Step 3 made interpenetration the way the molecule reads as a
      single body — **1893 of 2000 bonded pairs overlap** (94.65%), over 500
      universes — so a stick between two true-size atoms is buried inside them
      and shows on about **9%** of bonds. Sticks and overlap are close to
      mutually exclusive and each carries what the other loses, so `Solid` keeps
      Step 3's picture exactly and `Sticks` shrinks the atoms. The view joins the
      scene's change token rather than getting a second mechanism.

      **The drawn size is derived from the molecule, not hard-coded** — Ian's
      correction, and it was right. A first version used a constant 0.5,
      justified against a survey minimum of 0.872035 across 500 universes: that
      is correct-on-average and unjustified for any *particular* molecule, and
      wrong for anything Step 5 builds. Two drawn spheres separate exactly when
      `s < d / (rₐ + r_b)`, so `stick_geometry` takes the smallest such ratio
      over **this** molecule's bonds. Baked in `build` under §8.6, never on the
      way to a frame.

      **The palette is keyed on `group`, stepped by the golden angle, and the
      reason is a measurement that killed the obvious design.** An affinity ramp
      gives the five demo atoms a minimum pairwise separation of **0.0239** —
      four leaves collapsed into one colour — because the leaves land at
      *consecutive* groups `(0,1), (0,2), (0,3), (0,4)` in every universe
      measured. Golden-angle stepping on `group` gives **0.1932** across the
      five and **0.6748** among the leaves, and it is the coordinate that makes
      families visible, which §7.1 calls the payload. Saturation carries
      `affinity` (radius order equals mass order in 500/500, so a mass hue would
      restate the sphere sizes; it equals affinity order in **0** of 500);
      lightness carries `mass`; `valence == 0` is grey.

      **Three things nobody should re-derive.** `Demo::elements` is build order
      and disagrees with canonical order in **500 of 500** universes, in the
      permutation that paints the centre's colour onto a leaf. The demo molecule
      **cannot discriminate a bond lookup at all** — its bond set is one
      constant across the corpus, and six wrong implementations agree with the
      right one — so `mod shapes` builds a path, a ring, a tree and a two-hub
      graph inside the test. And `every_atom_is_drawn_where_the_embedding_put_it`
      was **rewritten, not deleted**, when the display scale made it fail.

      583 tests. `every_atom_is_painted_with_the_same_material` was retired
      deliberately, with its five replacements recorded at the deletion site;
      `every_atom_is_drawn_with_the_same_mesh` is kept, and is now the only
      closed channel.

      **Still manual, and stated rather than hidden.** Nothing catches "the
      molecule is too far away" (unchanged from Step 3), and nothing catches
      *colour legibility* — five distinct bit patterns is not five
      distinguishable colours on a screen, or to a colour-blind child. A
      distance bar would have to be measured against the shipped palette before
      it meant anything, which is a Step 4b job.
- [ ] **4b. Open the atoms up, and name what the generation already separates.**
      *New, from Ian on 2026-08-04, and it is mostly making visible what already
      exists rather than adding physics.*

      He asked for sub-atomic particles, for elements with distinct characters
      the way real ones have, and for something answering to the inert elements
      of the real table. Three of those are already computed and simply never
      drawn.

      **A warning about vocabulary, because the first draft of this step failed
      it.** Explaining the analogy to Ian in conversation is fine and necessary.
      Writing a real-world family name into the *plan* is how it becomes a label
      in the UI, and a label mapping an invented element onto a real chemical
      family is exactly the Borbax-to-real mapping §5 forbids. The step must name
      what the generator computes — a closed outer shell with no bonding contact
      left — and let the child discover the resemblance rather than assert it.

      - **The sub-atomic layer is the base unit.** An element *is* a cluster of
        `units` base units, and `period`/`group` are how they stack. Drawing the
        cluster — the units inside an atom, layer by layer, with the outer layer
        incomplete — shows the actual thing the generator computes. It is
        structurally the nucleus-and-shells picture without importing a single
        real-world particle.
      - **The inert elements are already there**: `valence == 0` is a closed
        outer shell with no contact left to make, so those elements bond with
        nothing. The properties block already prints "none — a closed shell".
        They need *marking in the table*, not inventing — and the mark must be
        drawn from that fact, never from a real family's name.
      - **Families should FALL OUT of the generation, never be declared.** Real
        families exist because of how the outer shell fills, and Borbax generates
        that same quantity — so the groupings are to be *found* and named, not
        imposed. A hard-coded family list would be the viewer drawing the answer
        on, which is the failure this project's whole design is against. Measure
        first; if nothing separates, say so rather than inventing a split.

      **Context worth keeping, because it changes how to read the fiction
      guarantee.** Ian's goal is showing a child how the world works. He
      chose an invented universe because another model refused to help him build
      real chemistry — not because fiction was the aim. The guarantee still holds
      for *this* project and is still enforced; but a separate, real-chemistry
      viewer is explicitly on the table as its own thing, and "we cannot build
      real chemistry" is **not** a true statement of the constraint.
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
click through to a molecule, and spin it around in 3D — with a child beside him,
on a laptop, without a terminal.

Everything after that is the chemistry arriving in a window that already works.
