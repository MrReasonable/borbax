# Borbax — Product Requirements Document

**A simulation of chemical evolution in an invented universe.**

| | |
|---|---|
| **Status** | Draft for review |
| **Date** | 2026-07-26 |
| **Scope of this document** | V1 — abiogenesis up to protocells |
| **Author** | Ian Dominey, with Claude |

---

## 1. What this is

Borbax simulates the emergence of self-replicating life from non-living chemistry, in a universe whose physics and chemistry are **invented from scratch by a generative algorithm**.

The central premise is that chemistry is about shape. Two things react when their shapes fit together and their surface properties complement each other. Everything else — catalysis, metabolism, membranes, heredity — is downstream of that one idea. Borbax takes that premise literally: molecules are shapes, binding is a geometric fit test, and nothing about replication is written into the rules. If a self-copying system appears, it appears because the shapes allowed it.

The primary audience is one child, and the primary goal is wonder: to sit with someone and watch a planet's chemistry churn for a hundred million years until something in it starts making copies of itself, and to be able to point at the exact moment it happened.

### 1.1 The invented-chemistry premise

Borbax contains no real chemistry. Not as a disclaimer — as an architectural property.

There is no periodic table in this repository. There is a *periodic table generator*: given a 64-bit seed, it invents a set of elements with their own masses, valences, bonding rules, and reactive tendencies. Different seeds produce genuinely different chemistries, with different rules about what can bond to what.

This is not a compromise made for safety reasons that costs the project something. It is a better design for what we actually want:

- **You cannot look up the answer.** Real-chemistry simulations are haunted by the fact that we already know how the story ends. In Borbax nobody knows, including the person who built it. The discovery is real.
- **The chemistry itself becomes a thing to explore.** Learning that Vorium is wildly reactive and that anything in Group 4 makes stable rings is *part of the game*, not a prerequisite you must already have.
- **It generalises.** We can ask what happens under chemistries that Earth never had, which is a far more interesting question than replaying a simulation of Earth.

Section 5 specifies the guarantees that keep this property true as the codebase grows.

---

## 2. Prior art, and what is genuinely new

Artificial life is a fifty-year-old field with a great deal of excellent work in it. This section establishes what already exists, what Borbax adds, and — most valuably — the specific ways projects in this space are known to fail.

### 2.1 Artificial chemistries

| System | Approach |
|---|---|
| **AlChemy** (Fontana & Buss) | Lambda-calculus expressions collide and apply to one another; self-maintaining "organisations" emerge |
| **Typogenetics** (Hofstadter) | Strings translate into enzymes that then edit strings |
| **Squirm3 / OrganicBuilder** (Hutton) | 2D particles with a type, a state, and typed bonds; reaction rules are explicit rewrite rules |
| **Chemlambda** | Graph rewriting on a directed graph of typed nodes |
| **Combinatory Chemistry** (Kruszewski & Mikolov) | Combinatory logic with conservation laws over primitive counts |
| **SCL** (Varela, Maturana & Uribe) | Three particle species in 2D; the original demonstration of computational autopoiesis |
| **JohnnyVon** (Smith, Turney & Ewaschuk) | Discrete finite-state machines embedded in continuous 2D physics |
| **P systems** (Păun) | Multisets rewritten inside nested membranes, with transport across boundaries |

The pattern is clear: nearly all of these are **symbolic rewriting systems** — strings, graphs, or lambda terms — or **particles with typed bonds and explicit rules**. In both cases, *which things react is written down somewhere*, either as rewrite rules or as a type-compatibility table.

None of them derives reactivity from **generated geometry**. That is the gap Borbax occupies.

> **A warning from this literature.** "Self-organization in computation and chemistry: Return to AlChemy" (Mathis, Patel, Weimer & Forrest, *Chaos* 34(9):093142, 2024) re-ran the founding experiment of the field and partially falsified it. Level-**0** organisations turn out to be more robust than reported — they resist collapse to the trivial fixed point. Level 1 is weaker than that summary suggests: the paper says stability "varies widely", and that one organisation "was not truly stable".
>
> Level 2 does not replicate. Combining organisations drawn from *different* runs across 455 trials: **mutual destruction ~68%, one dominating the other ~27%, coexistence ~5%** (Figure 4B, read from the figure itself — the table is a raster that text extraction drops).
>
> Three conditions travel with that result and are easy to lose. It combined organisations from *different* runs; the other route Fontana & Buss used — two separable organisations within a single run — was **not attempted**, so this is a non-replication of one route out of two. And the categories are a Jaccard-similarity cut at 0.1.
>
> The implication for Borbax is direct and important. **We must not assume protocells will spontaneously compose into higher-order structures.** The single most-cited claim in artificial chemistry — that self-maintaining systems stack into hierarchies — is the field's weakest result. V2 onward must treat composition as an open research question, not as something that will fall out.

### 2.2 Digital evolution

Tierra (Ray), Avida (Ofria, Adami & Lenski), Polyworld (Yaeger), Geb (Channon), Chromaria, and DISHTINY (Moreno & Ofria) all evolve **programs**, with no physics beneath the genome. Borbax inverts this entirely: in V1 there is no genome at all, only chemistry, and any heredity must arise from it.

The most relevant recent work is **BFF / "Computational Life"** (Agüera y Arcas et al., [arXiv:2406.19108](https://arxiv.org/abs/2406.19108)), in which random byte tapes are concatenated and executed as self-modifying programs with **no fitness function whatsoever** — and replicators emerge spontaneously. This is the strongest available evidence that replication can arise from a substrate rather than being selected for, which is exactly Borbax's bet.

It also carries a pointed lesson: **SUBLEQ fails at this despite being Turing-complete.** Universality of the substrate is not sufficient; its *structure* determines whether replicators are reachable. That is a strong argument for caring enormously about the details of the shape model, which is where Borbax concentrates its design effort.

### 2.3 Emergent-complexity sandboxes

Lenia and Flow-Lenia (Chan), Particle Life, Neural Cellular Automata (Mordvintsev), and SmoothLife are the visual bar to clear — they are genuinely beautiful and they demonstrate how much richness simple local rules can produce. But they have no heredity, no discrete species, and no open-ended novelty. They are pattern formation, not evolution. Borbax should look at least this good and do considerably more.

### 2.4 Games in this space

Thrive, Species: Artificial Life Real Evolution, Cell Lab, The Bibites, Evolution (Keiwan), and Spore all operate at the *organism* layer with hand-authored parts. Chemistry is decorative or absent. The consistent criticism — Spore most notoriously — is that the evolution turns out to be scripted or shallow once players look closely. Borbax's answer is that there is nothing to look behind: the mechanism *is* the chemistry, all the way down.

### 2.5 The one real precedent for the folding design

The genuine precedent for Borbax's shape map is the study of **sequence-to-secondary-structure maps** by Fontana, Schuster and colleagues, which established three properties as consequences of a many-to-one folding map:

- **Neutral networks** — many sequences fold to the same shape, so a lineage can drift extensively without losing function and then innovate from a new position
- **Shape-space covering** — every common shape is reachable within a small mutational radius of essentially any sequence
- **Plastogenetic congruence** — a sequence's alternative folds predict what its mutants will fold to

These three properties are precisely what make a shape map *evolvable* rather than merely complicated, and they are adopted directly as design targets for §8.4. Critically, they are **measurable in the beaker harness**, so we can verify the folding model has them before building anything on top of it.

But that work studied a folding map in isolation. It was never embedded in a spatial world with energy gradients, compartments, and geology. That embedding is the gap.

### 2.6 What is genuinely new here

Each ingredient exists somewhere. The combination does not:

1. **Procedurally generated chemistry** — an invented periodic table per seed, rather than a fixed rule set
2. **Geometric complementarity as the sole reaction mechanism** — no rewrite rules, no compatibility tables
3. **A multi-scale engine that moves compute toward novelty** — deep time made affordable by spending it selectively
4. **An unbroken causal chain** from element generation to protocell, with no hand-authored layer anywhere in between

**Independent corroboration of §8.3, from a substrate sharing no machinery with ours.** RBN-World (Faulconbridge, Stepney, Miller & Caves, "RBN-World: The Hunt for a Rich AChem", *Artificial Life XII*, MIT Press 2010, pp. 261–268) builds an artificial chemistry whose atoms *are* random Boolean networks. They enumerated six bonding *properties* against five *criteria* — 20 valid combinations — and put each through five tests: synthesis, self-synthesis, decomposition, substitution, catalysis. **Only the complementarity criteria survived**: `Total`/`Sum Zero`, which is `p(N_i) + p(N_j) = 0`, and `Proportion`/`Sum One`, its cycle-normalised sibling (identical only where the two networks share a cycle length, so they are two distinct rules of one family rather than a rescaling). Every similarity-based criterion — 'equal', 'similar', 'different' — died.

The often-quoted "200 chemistries" is 20 rules × 10 substrate settings, and the paper's own headline is that the substrate settings do not matter — so the honest statement of the filter's power is **2 of 2 complementarity rules survived and 0 of 18 similarity rules did**, which is still striking and is a smaller claim than "200".

`Total`/`Sum Zero` **is** our §8.3 charge term. That kernel scores `−(a_A[i] + a_B[j])²`, maximised exactly when `a_A + a_B = 0`. A five-test filter over 200 chemistries on a completely unrelated substrate selected the rule this project's central mechanism already uses, which is the strongest external support §8.3 has.

Two further transfers. Their other headline — that the substrate's own parameters "have little or no influence on the low-level properties of the chemistry", while the *bonding rule* decided richness — argues that Borbax's richness will be decided by §8.3 and not by the element table, which is worth remembering whenever the table looks like the problem. And **self-synthesis was by far their most discriminating test** (183 candidates down to 20, more than any other single test): a molecule binding a copy of itself. In Borbax that is self-complementarity under the 60-rotation search, non-trivial given §22.8's handedness — and §7.2's battery currently has no analogue of it. Adding one is cheap and the precedent says it discriminates.

**Their test was two-sided, and copying only the floor loses most of its power.** From the paper's Table 4, self-synthesis eliminated 163 of 183 — but **110 of those were chemistries where *every* sample self-synthesised**, discarded for showing no variation, against only 53 where none did. A floor ("reject if no species binds a copy of itself") captures the 53. The 110-analogue here is the sludge case: if essentially every species binds a copy of itself, complementarity carries no information and §22.8's handedness has bought nothing. Any criterion derived from this must be a **band**, not a floor — the same defect §7.2's polymer criterion had, one row above it in the same table.

**The caveat that limits it, and it lands on §8.3's boldest claim.** The paper says in its own Results: *"this work has only sampled from atomic constituents; it is not guaranteed that molecular structures will also exhibit these behaviours… molecular bRBNs of many atoms may not behave as an equivalent large bRBN atom **due to the constrictions from reciprocal bonding sites between atoms**."* §8.3's central idea is that **one affinity function works at every scale** — molecule–molecule, enzyme–substrate, membrane, pore. RBN-World tested atoms only and declines to support the lift to composites. That is the same shape as this section's governing lesson: AlChemy's level 0 and 1 replicated and the composition step did not.

**How squarely it lands, stated honestly.** An earlier version of this paragraph truncated the quote immediately before the authors' *reason*, and the reason is the load-bearing half: their obstacle is reciprocal bonding sites constraining how RBNs compose by rewiring, which has no analogue in an affinity kernel over shape signatures. So the caveat is weaker against §8.3 than the truncation made it look. It is retained rather than dropped because §2.7 asks us to err toward caution on exactly this claim — but a reader must be able to weigh it, and a quotation cut before its "because" cannot be weighed.

The authors frame their tests as *"a necessary, but not sufficient, basis for higher-order emergent properties"* — surviving the filter means "more likely to be capable of rich emergent properties", not rich. (An earlier version of this line quoted them as saying "necessary, but possibly insufficient". The hedge was ours, inside their quotation marks; they assert insufficiency flatly.) And their own previously published parameterisation failed their own battery, which is the best free argument for §7.2 existing at all.

*Two hygiene notes.* Do not repeat that paper's "less than 5%" survival figure: its Table 5 lists 19 survivors of 200, which is 9.5%, and the text contradicts the table. And this is **a single unreplicated 8-page conference paper** — in a section whose governing lesson is a widely-cited early result that failed to replicate, that has to be said out loud rather than left as an implication.

RBN-World is **not** an analogue of building atoms by packing a base unit — its atoms are Boolean networks. The resemblance is methodological, and it belongs against §7.2: it is the closest published precedent for a cheap low-level battery used to reject dead chemistries before running them.

### 2.7 Failure modes to design against

The most valuable output of the survey. These are the documented reasons projects in this space stall, and each gets a specific countermeasure.

| Failure mode | Countermeasure |
|---|---|
| **Uncreative unbounded (Bedau class 3b)** — unbounded activity concentrated in a *bounded* diversity of components, so cumulative activity climbs while nothing new ever appears. It passes a naive open-endedness test while doing nothing. (Bedau, Snyder & Packard 1998: class 1 is no adaptive activity, class 2 is *bounded* activity, class 3 is unbounded — and 3b is the subcase where diversity is bounded.) | Track novelty as a metric entirely separate from activity; never infer one from the other. Detecting 3b specifically needs the persistence-weighted Bedau–Packard statistics, not a raw species count (§15.2) |
| **The boundedness illusion** (Wiser, Dolson, Vostinar, Lenski & Ofria) — plateaus are usually declared by eyeballing a flattening curve. Judged by fitting and comparing competing models instead, the LTEE fitness trajectory is best fit by an *unbounded power law* | Never declare a plateau from a graph. Fit competing models and compare them (§15.3) |
| **No neutral baseline** — activity statistics are uninterpretable without a control | **Neutral shadow runs** (Rechtsteiner & Bedau): an identical run with selection switched off, used to calibrate thresholds and subtract drift. Adopted as a V1 deliverable (§15.3) |
| **Gameable metrics** (Hintze, "Open-Endedness for the Sake of Open-Endedness") — a trivial system can satisfy every published open-endedness criterion | No single metric is trusted. Every novelty measure is read against the neutral shadow |
| **Complexity carrying capacity is physical** (Moreno & Ofria) — it scales with space-time volume, so a too-small world caps complexity for thermodynamic reasons rather than tuning reasons | Treat region size as a complexity *budget*. If a run stalls, suspect this before touching parameters (§20) |
| **Composition is not guaranteed** — per the AlChemy correction in §2.1 | Never design a version on the assumption that lower-level units will stack |

---

## 3. Design principles

1. **Nothing about life is hardcoded.** No `struct Cell`. No replication rule. No "if conditions are right, spawn a protocell". Every biological structure must be a name we give to a pattern the chemistry produced on its own. If we find ourselves writing a special case for life, the physics is wrong and we fix the physics.

2. **Shape is the only mechanism.** Binding, catalysis, membrane formation, and selective permeability are all the same operation — a complementarity test between two surface descriptors — applied at different scales. One mechanism, reused. If a feature needs a second mechanism, question the feature.

3. **Spend compute where the surprise is.** The overwhelming majority of a planet is boring for the overwhelming majority of its history. Simulating it uniformly is the mistake every ambitious simulation makes. Borbax runs a cheap approximation almost everywhere and promotes a region to full molecular detail only when something novel starts happening there.

4. **The clock is event-driven, not tick-driven.** We never step through time waiting for something to happen. We ask "what is the next event anywhere, and when?" and jump straight to it. Quiet epochs cost almost nothing, which is what makes deep time affordable.

5. **Deterministic to the bit.** Same seeds, same configuration, same result — always. This is not a nice-to-have; it is what makes rewinding, branching, sharing, and debugging possible at all.

6. **Load the dice, don't rig the game.** The environment provides the ratchets that real prebiotic chemistry is thought to have had — energy gradients, wet/dry cycles, mineral surfaces, confined pores. These make self-organisation *probable* instead of astronomically unlikely. But they are environmental conditions, not shortcuts. The chemistry still has to find its own way through.

7. **Everything decays.** Nothing is permanent and no structure is exempt. Decay is not a balancing mechanic added afterwards to stop numbers running away — it is the precondition for selection meaning anything at all (§9.4). The region is an open, driven system with matter and energy flowing through it, because a closed system at equilibrium is by definition dead.

8. **Measure against a null model.** No claim about emergence is trusted without a neutral shadow run to compare it to (§2.7). A number that looks impressive on its own means nothing.

9. **Build the debugger before the telescope.** The schematic, clinical view of what the engine thinks is happening gets built first, because without it nothing else can be trusted. Beauty is layered on top of a simulation we already believe.

---

## 4. Locked decisions

These were settled before drafting and the rest of the document assumes them.

| Decision | Choice |
|---|---|
| **Emergence model** | Open rules, loaded dice. Nothing about replication is hardcoded; the environment supplies realistic ratchets that make it probable. |
| **V1 scope ceiling** | Stop at protocells. Genomes, cells, multicellularity and organs are V2+. |
| **Platform** | Rust throughout — simulation core plus a native `egui` + `wgpu` viewer (§14.0). No second language. |
| **Interaction model** | Observatory and time machine. Watch, zoom, scrub through deep time. Intervention tools come later. |
| **Shape model** | Hybrid — small molecules are atom graphs reduced to a directional shape signature; polymers physically fold into cavities. |
| **Geometry** | Three-dimensional. Signatures sample a subdivided icosahedron; folding runs on an FCC lattice (§8.2, §8.4). |
| **Seeding** | Two seeds. `universe_seed` generates the physics and chemistry; `world_seed` generates a planet under those laws. |
| **Run length** | Overnight. First self-replicator within roughly 8 hours on the target machine. |
| **V1 playfield** | One focused region (~256 km²) at high fidelity. The 1,000,000 km² planetary field arrives in V3 when complex life needs the room. |
| **Interface depth** | Scales with the reader — one interface, switchable between storybook, explorer, and scientist registers. |
| **Visual treatment** | Schematic scientific renderer first; naturalistic presentation layered over the same data later. |
| **Decay model** | Global cleavage rate to locate the band, then per-bond rates derived from the generated energy matrix (§22.6). |
| **Damage representation** | Implicit — a damaged molecule is simply a different molecule. Damage load measured as divergence from RAF-implied composition (§22.7). |
| **Neutral shadow** | Spawned on demand from a keyframe, not run continuously (§22.3). |
| **Platforms** | macOS, Windows and Linux, on aarch64 and x86-64, with bit-identical results everywhere (§13.4). Enforced by a CI golden matrix (§13.6). |
| **GPU** | wgpu when it lands, as a *screening* accelerator only — never a canonical result path (§13.5). |
| **Chirality** | The binding search covers rotations but **not** reflections, so Borbax chemistry is handed from the start (§22.8). |

---

## 5. Fiction guarantees

Borbax must be structurally incapable of describing real-world chemistry or biology. These are binding requirements on the implementation, testable in CI, not statements of intent.

**G1 — No real chemistry data enters the repository.** No element property tables, no reaction databases, no molecular structures, no sequence data of any kind. The entire chemistry of every universe is produced by `universe_gen` from a seed. A CI check fails the build if any data file appears under the chemistry crates.

**G2 — Generated names cannot collide with real ones.** The element name generator checks every candidate symbol and name against a blocklist of all real element symbols and names, and against a list of common real chemical and biological terms. Collisions are rejected and regenerated. The blocklist exists solely to *exclude*; it is never read as data.

**G3 — The physics is deliberately non-isomorphic to real chemistry.** Period lengths are generated per universe rather than following any real shell structure. Valence rules are generated. The `affinity` scalar governs surface complementarity and is not electronegativity — it has different units, a different range, and different behaviour. The folding model is a lattice heuristic tuned for interesting shapes, not a physical force field. There is no transfer function from a Borbax molecule to any real molecule, and constructing one would require inventing the mapping from nothing.

**G4 — No real-world units anywhere.** Temperature is measured in *thermals*, energy in *quanta*, distance in *spans*, time in *world-years* defined by the generated world's own orbit. Unit names are enforced by the type system. This is not cosmetic: it prevents anyone, including us, from quietly reasoning about the simulation as though it were a model of anything real.

**G5 — Permanent non-goals.** Borbax will never import or export real chemical formats (SMILES, InChI, MOL, PDB, FASTA or similar), will never contain a mapping table between Borbax entities and real-world entities, and will never accept real-world measurements as calibration targets. Any future feature request implying otherwise is rejected on sight.

**G6 — Outputs are not predictive of anything real, by construction.** A Borbax result is a fact about Borbax. Because the elements, bonding rules, energetics, and folding model are all invented and mutually entangled, no result transfers. The documentation states this plainly and the README leads with it.

Note that §2.5 cites prior work on folding maps for a *mathematical property* — that many-to-one maps produce neutral networks — and nothing else. No real structural data is used, referenced, or required.

Together these mean the interesting question — *what does it take for chemistry to come alive?* — can be explored freely, because the chemistry in question is one we made up.

---

## 6. Architecture

```
                        universe_seed (u64)
                               │
                    ┌──────────▼──────────┐
                    │   universe_gen      │   "the laws of physics"
                    │  elements, bonding  │   generated once, then immutable
                    │  energetics, folding│
                    └──────────┬──────────┘
                               │
        world_seed (u64) ──────┤
                               │
                    ┌──────────▼──────────┐
                    │    world_gen        │   "one planet under those laws"
                    │ terrain, vents,     │
                    │ climate, abundances │
                    └──────────┬──────────┘
                               │
   ┌───────────────────────────▼───────────────────────────┐
   │                    SIMULATION CORE                     │
   │                                                        │
   │   ┌─────────┐   promote   ┌────────────┐   promote     │
   │   │  BULK   │ ──────────► │ STOCHASTIC │ ──────────►   │
   │   │ ~10⁶    │ ◄────────── │  ~10³      │ ◄──────────   │
   │   │ patches │   demote    │  patches   │   demote      │
   │   └─────────┘             └────────────┘               │
   │        reaction-diffusion    tau-leaping                │
   │        years/step            seconds/step               │
   │                                    ┌──────────────┐    │
   │                                    │  MOLECULAR   │    │
   │                                    │  ~10 patches │    │
   │                                    └──────────────┘    │
   │                                    every molecule      │
   │                                    microseconds/step   │
   │                                                        │
   │   ┌────────────────────────────────────────────────┐   │
   │   │ Novelty detector · RAF detector · Governor     │   │
   │   └────────────────────────────────────────────────┘   │
   └───────────────────────┬────────────────────────────────┘
                           │ keyframes + event journal
                           │ in-process — no serialisation boundary
              ┌────────────▼────────────┐
              │     borbax-geometry     │  projection, layout, colour
              └──────┬───────────┬──────┘
                     │           │
        ┌────────────▼──┐   ┌────▼─────────────────┐
        │ SVG           │   │  egui + wgpu         │  Lab → Wonder skin
        │ headless,     │   │  native application  │  depth dial
        │ golden-tested │   │  The Chronicle       │
        └───────────────┘   └──────────────────────┘
```

A world is addressed by both seeds together **and the physics version that generated it**, written `U-7F3A21C9@1 / W-0004`. Share that and someone else gets a bit-identical planet.

The version is not decoration. Today the true identity of a world is `(universe_seed, world_seed, physics, compiler)` and only the seeds are written down — which is fine while there is exactly one set of laws, and stops being fine the first time a law is added. A universe carries the physics it was generated under and is replayed under those laws forever; new laws mint a new version rather than reinterpreting old worlds. The mechanism and the rules that make an added law additive rather than destructive live in the V0 plan, at Tasks 3 and 5 — an earlier version of this sentence cited §7.1, which is the generated periodic table and contains none of it.

Two things this address does **not** promise, stated here because the overclaim is tempting. It records the *physics*, not the *toolchain*: the true identity is `(universe_seed, world_seed, physics, compiler)`, and a rustc or `libm` bump still moves results under an unchanged `@1`. What versioning buys is **deliberateness** — when a replay diverges you can tell whether the laws changed — not bit-identity across time, which additionally needs the pins. And it versions the *universe*, not the world: changing world generation changes what `W-0004` produces while leaving the chemistry alone, which is a separate and currently unversioned kind of breakage. Keep the universe seed and change the world seed to get a new planet under chemistry you already understand. Change the universe seed to start over in a universe where nothing you learned applies.

---

## 7. Layer 0 — Generated physics

`universe_gen` turns a 64-bit seed into a complete, self-consistent chemistry. It runs once per universe; the output is immutable and cached.

### 7.1 The generated periodic table

The critical design choice is that the table has **genuine periodicity**. Purely random elements would give random chemistry, which is unlearnable and therefore worthless — there would be no patterns to notice, and noticing patterns is the entire educational payload.

So the generator invents a shell structure and derives properties from position:

1. Pick a **shell pattern** — a sequence of period lengths drawn from the seed, e.g. `[3, 7, 7, 15, 22]`. Different universes get genuinely different table shapes.

   **Draw the lengths, do not pick them from a menu.** G3 requires that period lengths do not follow any real shell structure, and a fixed candidate list of `{2, 8, 18}` ∪ `{6, 10, 14}` *is* real shell structure — real period lengths unioned with real subshell capacities — however randomly it is sampled. The generated-ness has to be in the numbers themselves. Note the `xtask` guarantee check cannot catch this: it scans for data *files*, and a violation of this kind lives in a `const`.
2. Lay elements out in periods and groups according to that pattern.
3. Derive each property as a smooth function of (period, group) plus a small amount of seeded noise.

The result is a table with **families**: elements sharing a group behave alike. This produces the moment we actually want — *"wait, everything in this column makes rings"* — which is the moment a person starts to understand a chemistry rather than just watch it.

**The trends must have generated *direction*, not just generated magnitude.** It is not enough that `affinity` rises across a period and `radius` shrinks across one while growing down a group — those are the real trends, in the real direction, and reproducing them makes G3 an assertion rather than a fact. The sign of each trend is drawn from the seed, so one universe has affinity rising across a period and another has it falling. This costs nothing and is the difference between "invented" and "renamed".

Each element carries:

| Property | Meaning |
|---|---|
| `symbol`, `name` | Procedurally generated, pronounceable, blocklist-checked (G2) |
| `mass` | Generative atomic weight; increases down the table with noise |
| `valence` | Bonding slots available, 0–6; a periodic function of group |
| `affinity` | Surface-interaction scalar in [−1, +1]. Drives complementarity. Trends across a period. |
| `radius` | Steric size; contributes to shape extent |
| `bond_energies` | Per-partner-group bond strength matrix |
| `decay_rate` | Decay probability per world-year. Unstable elements are both an energy source and a mutation source. **Named for what it holds**: it is 0 at the binding-energy peak and rises toward the extremes, so calling it `stability` — as an earlier version did — inverted the sense of the field against its own description, in the quantity §9.4 calls the most sensitive parameter in the system. |
| `catalytic_class` | **Reserved, and carries no chemistry.** A rebinning of `group` — `floor(fill × n_catalytic)` — with zero information beyond it. It is *not* derived, it is the one §7.1 property that is not, and **nothing downstream may branch on it as though it decided which reactions an element promotes.** An earlier version of this row said exactly that, which was a standing licence for a catalysis path keyed on an element label, parallel to §8.5's cavity geometry and forbidden by §8.3's one-mechanism rule. Catalysis is geometry (§8.5) or it is not catalysis. |

### 7.2 Validation — is this universe worth simulating?

Most random chemistries are dead: everything is inert, or everything reacts instantly into sludge. Generating those wastes an entire overnight run, so `universe_gen` self-tests before declaring a universe valid.

It runs a **beaker battery** — a few hundred thousand simulated reactions in a single well-mixed volume, taking a second or two — and measures:

- **Reactivity band.** The fraction of collisions that react must fall in a workable range. Too low and nothing happens; too high and everything burns.
- **Decay band.** Spontaneous and solvent-driven cleavage must land inside a narrow window (§9.4). Too fast and nothing survives long enough to fold; too slow and there is no selection, only accumulation. This is the most sensitive parameter in the entire system, and finding it is much of the reason the battery exists.
- **Polymer viability, banded.** Chains of length ≥ 20 must form and persist *against the prevailing decay rate* — and the beaker must **not** have gelled. Persistence is a balance between formation and destruction, never a property of the molecule alone.

  The ceiling is not symmetry for its own sake: **a gelled beaker passes the floor.** Simulated on a configuration model of 20 000 units at bond-conversion extent 0.3, with a valence distribution a packing-derived element table plausibly produces, **51% of all mass sits in a single connected component** — and seven molecules still exceed 20 units, so the floor passes. That is the sludge failure this section names in its own opening sentence, passing the test meant to catch it. The measurable form is a ceiling on the fraction of mass in the largest connected component.

  The controlling quantity is the *realised* degree distribution's weighted functionality, `f_w = ⟨f²⟩/⟨f⟩`: a giant component appears at extent `p_c = 1/(f_w − 1)`. **That threshold is Cohen, Erez, ben-Avraham & Havlin (*Phys. Rev. Lett.* 85:4626, 2000), independently Callaway et al. the same year, and in the polymer literature it is the Flory–Stockmayer gel point (Flory 1941; Stockmayer, *J. Chem. Phys.* 11:45, 1943 for the multifunctional-mixture form — the 1944 paper treats cross-linking of pre-formed chains)** — Molloy–Reed (1995) gives only the *existence* criterion `Σ i(i−2)λ_i > 0`, equivalently `f_w > 2`, which is a yes/no test rather than a threshold in an extent coordinate. An earlier version of this paragraph attributed the threshold to Molloy–Reed. The conditions travelling with it are tree-like connectivity, equal reactivity, no degree correlation, and N → ∞; the first of those cuts in our favour, since bonds spent closing a cycle are wasted for connectivity, so a finite cyclic system gels at or later than the mean-field prediction.

  So the battery must measure `f_w` in the beaker rather than reading nominal valences off the element table.

  **But it must not report `p_c` alone, and the reason overturns an earlier claim here.** This paragraph used to argue that burial and shape complementarity make realised functionality lower than nominal, so realised `p_c` is higher and any nominal figure is conservative. The arithmetic is right and the physics is not. Under degree-independent thinning at rate `q`, `f_w' = 1 + q(f_w − 1)` exactly, so `p_c' = p_c/q` — and the invariance is exact, not statistical: bonds-per-node at gel is `p_c⟨f⟩/2 = ⟨f⟩/(2(f_w−1))`, in which `q` cancels identically. (An earlier version cited a simulation giving 0.4952 / 0.4957 / 0.4928. Those sit at the Erdős–Rényi value 0.5, i.e. the run used a Poisson-like distribution — the one family closed under thinning, where the invariance is trivial. On this section's own reference distribution the invariant is 0.668, and the `q = 0.5` case cannot be run at all, because thinning takes `f_w` to 1.68, below Molloy–Reed's `f_w > 2`, so no giant component exists at any extent. The algebra holds generally; the simulation did not demonstrate it.) `p_c` rises only because the extent coordinate was rescaled. Burial acts as a hard cap exactly when `f_w' ≤ 2`, i.e. when `q ≤ p_c` — so the removal required is `1 − p_c`, which on the measured range is **8.6% to 53%, median 26%**. (This sentence has now been wrong twice in the same way. It first said 60–80%, which was `1 − p_c` computed from a `p_c` range a frontier-geometry error had depressed and not recomputed when the error was fixed. Its replacement said 11–58%, median 35% — computed from that same stale range rather than from a re-run, *inside the parenthetical explaining the defect*. The figures above are `1 − {0.914, 0.741, 0.468}`, from the probe's `gel extents` line at HEAD. Requote them from a run, never from this paragraph.) Shape complementarity moves `p_c` by ±5% with the sign depending on the split.

  **The battery reports bonds-per-node at gel alongside extent**, because that measure is invariant under exactly the rescaling — measured on the generated table, abundance weighting moves extent by +85% and bonds-per-node by +9%, so roughly 89% of the apparent effect is coordinate.

  **But neither coordinate is sufficient on its own, and this is the correction that matters.** Both are *thresholds*. Reporting a second threshold changes nothing, because the ratio of where-the-beaker-sits to where-it-gels is the same number in either coordinate. What both versions of this paragraph lacked is an **operating point**. With §9.1's bimolecular Condense against §9.4's unimolecular Cleave, steady state is `(1−p)²/p = k_r / (2·k_f·n·⟨f⟩)`, and putting that in reverses the conclusion again: burial moves the beaker *away* from gel, monotonically. The margin between `q` = 1.0 and `q` = 0.5 is `F(R) = 2·p_ss(R)/p_ss(2R)`, in which `p_c` and therefore `f_w` cancel entirely — so it depends only on the rate ratio `R = k_r/(2·k_f·n·⟨f⟩)`, and is **bounded in (2, 4)**: `F → 2` as `R → 0` (near-full conversion), `F → 4` as `R → ∞` (dilute), monotone throughout. An earlier version of this sentence said "~2.3×, robust across three orders of magnitude"; it is neither. `F = 2.3` pins `R ≈ 0.119`, and three decades centred there span 2.05 → 3.37. Worse, that `R` puts `p_ss = 0.710` against this section's own reference `p_c = 0.734` — a beaker at 96.7% of its gel threshold. The paragraph whose thesis is *"what both earlier versions lacked is an operating point"* quoted a figure that silently depended on one, and `R` is a §9.1/§9.4 quantity Task 4 cannot measure. **Quote the bound and its two limits; a single number here is always a hidden operating point.** So the original burial intuition was defensible after all — for a reason neither the original paragraph nor its first correction gave.

  **The battery must therefore report `p_ss / p_c`**, the margin between the operating point and the threshold. A universe is judged on that ratio, not on either coordinate alone.

  **`p_ss` must be *counted*, not evaluated from that formula.** The closed form assumes one forward constant for all functional groups and one reverse constant for all bonds — i.e. equal reactivity, which this design deliberately breaks: §9.1's Condense is affinity-dependent per species pair and §9.4's Cleave is `A·exp(−E_bond/T)` per bond, and broad distributions in both are the *point* (§8.3 makes binding shape-dependent; §9.4 calls the decay band the most sensitive parameter in the system). This paragraph is careful to list the four conditions travelling with `p_c` — tree-like, equal reactivity, no degree correlation, N → ∞ — and then smuggles equal reactivity back in one paragraph later, attached to the quantity it calls "the correction that matters". The bias is not obviously conservative either: weak bonds cleave preferentially, so if bond energy correlates with valence the collapse is degree-correlated as well. The remedy is not a better formula — the battery already runs the beaker, so `p_ss` is directly countable as reacted groups over total groups at steady state, unconditionally and for free. Keep the closed form as a cross-check; gate on the count. **A distribution dominated by f = 2 with f = 1 terminators and a thin high-valence tail is what produces chains; a broad or high-centred one produces a single cross-linked blob**, and the two are not distinguishable by any floor. Reference points, stated here rather than left in a merged PR body:

  | distribution | `f_w` | `p_c` |
  |---|---|---|
  | 30% f=1, 68% f=2, 2% f=8 | 2.36 | 0.734 |
  | all f = 8 | 8.00 | 0.143 |
  | uniform 6..10 | 8.25 | 0.138 |

  **PR #8 quoted "centred near 8, `f_w` = 6.01, gels at 0.200" and no such distribution exists.** Cauchy–Schwarz gives `f_w = ⟨f²⟩/⟨f⟩ ≥ ⟨f⟩`, so anything with a mean near 8 has `f_w ≥ 8` and `p_c ≤ 1/7 = 0.143`. A `f_w` of 6.01 implies a mean of at most 6.01 and is therefore not centred near 8.

  That refutes the *mean* reading of "centred"; the mode reading falls too. Maximising mass at `f = 8` subject to `f_w = 6.01` puts at most **36.2%** there (the optimum is 36.19% at f=8, 63.81% at f=3.00). No distribution with that `f_w` can put even 40% of its mass at 8 under any reading. The 0.734 row is sound.
- **Shape diversity.** The realised shape signatures must spread across signature space rather than collapsing into a few clusters.
- **Neutral-network structure.** Per §2.5 — the folding map must exhibit neutral networks and approximate shape-space covering. A universe whose folding map is close to one-to-one is rejected: it cannot support drift, and lineages in it will be brittle.
- **Catalytic potential.** At least some folded polymers must produce cavities that complement common small molecules — and those cavities must be *enclosed* enough to be selective, not merely present.
- **Energy landscape.** There must exist reactions in both directions across a usable energy range — a chemistry that only ever runs downhill cannot build anything.

Universes failing the battery are rejected and the seed is advanced. The battery is also the tuning harness: it makes the chemistry-design loop seconds long instead of hours long, which is the difference between a tractable project and an intractable one.

---

## 8. Layer 1 — Molecules and shape

### 8.1 Small molecules

An atom graph, up to about 12 atoms. Stored as an arena of atom references plus a bond adjacency list with bond orders.

Identical molecules must be recognised as identical, so every molecule is reduced to a **canonical form** by a canonical graph-labelling pass, hashed, and interned. Species identity is then a `u32` and comparison is free.

### 8.2 The shape signature

This is the heart of the engine and everything else depends on it.

A molecule's graph is embedded in **three dimensions** by stress majorization over graph-theoretic distances, initialised deterministically from the canonical atom order. This is stable under small graph changes — a small mutation produces a small shape change — and that property is what makes the chemistry *evolvable* rather than merely complicated.

Laplacian eigenvectors were the obvious first choice and are the wrong one. Eigenvectors have arbitrary sign, and degenerate eigenvalues — extremely common in the symmetric graphs that small molecules actually are — leave the eigenspace basis arbitrary. Since the signature keys species identity, the fold cache, and the shape-space novelty histogram, that arbitrariness would propagate straight into the parts of the system that most need to be well-defined. Stress majorization from a canonical starting order sidesteps the problem rather than managing it.

The embedded molecule is then sampled along **D directions distributed over a sphere** — the vertices of a subdivided icosahedron, so D ∈ {12, 42, 162}. Each direction records two numbers:

- `r` — how far the molecular surface extends that way (shape)
- `a` — the summed affinity of atoms facing that way (surface character)

```
   Sample directions = vertices of a geodesic sphere
   (icosahedron, subdivided)

            ·  ·  ·                signature =
         ·           ·               [(r₀,a₀), (r₁,a₁), … (r_D₋₁,a_D₋₁)]
       ·    ╭─────╮    ·
      ·    │  mol  │    ·          At D = 42 that is 84 floats — the whole
       ·    ╰─────╯    ·           molecule, as far as binding is concerned.
         ·           ·
            ·  ·  ·

   12 verts (icosahedron) · 42 (1 subdivision) · 162 (2 subdivisions)
```

**Why an icosahedral sampling specifically.** The rotation group of the icosahedron has 60 elements, and every one of them maps the vertex set onto itself. A rotation is therefore an exact **permutation of the sample directions** — precomputable as a lookup table, with no interpolation and no resampling error. With an arbitrary set of directions, rotating a signature would require interpolating between samples, which is both expensive and lossy, and the loss would land squarely on the comparison that decides whether two molecules bind. The geometry and the algebra line up, and that is not a coincidence worth passing up.

**The signature is then rotation-canonicalised**: of its 60 rotations, the lexicographically smallest is stored. Without this, the same molecule could produce different signatures depending on how it happened to be embedded, and species identity would stop being well-defined. Binding is unaffected, because only *relative* orientation matters there — it searches rotations regardless.

### 8.3 Binding — one operation, used everywhere

Two molecules bind when their signatures are **complementary**: their shapes interlock and their surface characters oppose.

```
For each of the 60 rotations R in the icosahedral rotation group:
    for each direction i:
        j = ANTI[ PERM[R][i] ]                     // exact, precomputed
        shape_term  = −( r_A[i] + r_B[j] − IDEAL_GAP )²
        charge_term = −( a_A[i] + a_B[j] )²
    score(R) = Σᵢ ( w_shape · shape_term + w_charge · charge_term )

affinity(A,B) = max over R of score(R)
P(bind)       = sigmoid( affinity / temperature )
```

**The `ANTI` lookup is load-bearing and easy to omit.** Two bodies in contact touch along *opposite* directions: A's surface reaching along `+d` meets B's surface reaching along `−d`. Pairing `r_A[i]` with `r_B[PERM[R][i]]` directly — the obvious form — compares each molecule against the *point inversion* of its partner.

The error does not announce itself, and it is not absorbed by maximising over rotations. The sample directions are antipodally closed, but `−I` has determinant −1 and so is **not** a member of the 60-element rotation group. Omitting `ANTI` therefore searches the 60 *improper* elements of the full icosahedral group — the exact opposite of §22.8's deliberate choice, and with nothing failing to signal it. Bumps would meet bumps, cavity recognition would select the enantiomer of the molecule that actually fits, and any homochirality result would report the wrong handedness.

A test built by constructing a complement at the *same* index cannot catch this, because it encodes the same mistake. The complement must be built through `ANTI`, from a real embedded molecule.

At D = 42 the search is 2,520 iterations per pair — around forty times the two-dimensional version, and affordable only because of §8.6: `affinity` is a pure function of a species *pair*, so it is memoised and runs once per novel pair rather than once per collision.

Where the memo misses, a rotation-invariant upper bound on the score rejects hopeless pairs before the full search runs. It must be a genuine bound — a filter that cannot state the property it guarantees is not conservative, it is merely untested.

Every higher-level phenomenon in Borbax is this function called at a different scale:

| Phenomenon | What it actually is |
|---|---|
| Two molecules sticking together | `affinity(A, B)` above threshold |
| An enzyme recognising a substrate | `affinity(cavity, molecule)` |
| A membrane self-assembling | Many molecules with mutually complementary flanks |
| A pore letting something through | `affinity(pore, molecule)` gating passage |

One mechanism. This is the design working.

### 8.4 Polymers and folding

Polymers are sequences of monomer units, and unlike small molecules they **fold**.

Folding runs on a **face-centred cubic lattice** by simulated annealing against a contact-energy model derived from monomer affinities. FCC gives each site twelve nearest neighbours — the densest packing available on a lattice — which matters because it is what lets a chain fold genuinely *compactly* rather than into something loose and stringy. Coordinates are integer triples constrained to even parity, so the twelve neighbour offsets are compile-time constants and self-avoidance is a single array lookup.

Folding is the expensive operation in the engine, so folds are aggressively cached. A fold is a pure function of the sequence, so the cache is a sequence-hash → conformation map with high hit rates once a population stabilises.

The folded conformation is then scanned for **cavities** — connected empty regions substantially enclosed by the chain — and each cavity gets its own directional signature computed from the monomers surrounding it.

```
  V K Z V V K R Z T V Q V K
              │  anneal on FCC lattice (12 neighbours per site)
              ▼
         ▓▓▓▓▓▓▓            buried core: monomers with few
        ▓▓░░░▓▓▓▓           empty neighbours, protected from
        ▓▓░ A ░▓▓▓▓         solvent attack (§9.5)
        ▓▓░░░▓▓▓▓
         ▓▓▓▓░░▓▓           A, B = cavities. Each is enclosed on
          ▓▓░ B ░▓          most sides and reachable through a
           ▓▓░░▓▓           mouth — a real binding site, not a
            ▓▓▓▓            notch on a perimeter.
```

**This is the reason for three dimensions.** In two, a chain of *n* monomers has roughly √n interior against √n boundary, so there is almost no inside to be inside of: burial barely protects anything and a "pocket" is a notch on an edge, open on two sides and correspondingly unselective. In three, the same chain has *n* interior against n^⅔ surface. Burial becomes genuinely protective, and a cavity can be enclosed on many sides at once — which is what makes a binding site *specific* rather than merely sticky. The catalysis mechanism in §8.5 depends on that specificity, and so does V0's exit criterion 6.

**The folding map is many-to-one by design.** Per §2.5, this is what produces neutral networks — many sequences folding to the same shape — and it is the single most important property for evolvability, because it lets a lineage accumulate variation without losing function and then innovate from somewhere new. The beaker battery measures it explicitly (§7.2), and a universe whose folding map lacks it is rejected before a run ever starts.

### 8.5 How catalysis emerges

Catalysis is not implemented. It falls out.

A polymer with a single cavity that complements some small molecule will bind it — that is a receptor. A polymer with **two cavities close together** will bind two molecules and hold them adjacent, in a fixed relative orientation, for as long as the complex persists. That proximity and orientation is precisely what lowers a reaction's activation barrier.

So the rate enhancement is computed geometrically from cavity separation and alignment, and an enzyme is simply a folded chain that happened to end up with two well-placed cavities. Nothing declares it an enzyme. We detect that it is one, and the Chronicle reports it.

A cavity enclosed on many sides is selective in a way a surface dimple is not — it admits molecules of roughly the right shape and excludes everything else. That selectivity is the whole difference between a catalyst and a sticky patch, and it is the concrete payoff of §8.4's three dimensions.

This is the moment the whole design pays off, and it is the first headline event of a run.

### 8.6 Everything expensive is per-species, not per-molecule

This invariant is what makes the performance budget in §17 reachable at all, and it is easy to violate by accident, so it is stated here rather than left implicit.

Canonicalisation, the 3D embedding, signature construction, folding, cavity extraction, and per-bond decay rates are **all pure functions of the species**. None depends on which particular copy of a molecule is under consideration, or on when. So each runs exactly once, when a species is first interned, and never again.

Concretely: **no code reachable from a simulation step may call any of them.** A species record computed at intern time carries everything the hot loop needs — signature, pre-rotated binding view, mass, bond inventory, precomputed solvent-attack rate, fold reference. The step loop reads that record and does arithmetic on it.

Two consequences worth stating because they are not obvious:

**Reaction products are memoisable.** The products of a given (species, species, reaction class) triple are deterministic, so they are resolved once when the reaction channel is created. Canonicalisation therefore runs only when a genuinely *novel* species appears — a rare event, and one we already want to count as the novelty metric of §15.2.

**Binding affinity is memoisable, and safely so.** `affinity(A, B)` is a pure function of two species plus universe constants; temperature enters only afterwards, at the sigmoid. Memoising a pure function cannot change results, unlike a lossy cache, so this is safe under §13.1. Because concentrations are heavily skewed, a dense table over the most abundant species captures the overwhelming majority of queries.

If a later change makes any of this per-molecule, the molecular tier's budget becomes unreachable by roughly two orders of magnitude — and that kind of regression usually arrives disguised as a small feature.

---

## 9. Layer 2 — Reactions, energy, and decay

### 9.1 Reaction classes

Five classes, all derived from shape and energetics rather than enumerated:

| Class | Description |
|---|---|
| **Associate** | Non-covalent complex forms via signature complementarity |
| **Condense** | Covalent bond forms between two molecules; consumes energy, releases a small leaving group |
| **Cleave** | Covalent bond breaks; releases energy. Fires **spontaneously** at a temperature-dependent rate, not only under catalysis — this is the engine's primary decay path (§9.5) |
| **Transfer** | A sub-group migrates between molecules within a complex |
| **Rearrange** | Internal bond topology changes without a partner |

Rates follow an Arrhenius-style form:

```
rate = A · exp( −Ea / T ) · [X] · [Y] · catalysis_factor
```

where `Ea` comes from the generated bond-energy matrix, `T` is patch temperature in thermals, and `catalysis_factor` is the geometric enhancement from §8.5.

### 9.2 Energy

Energy is conserved and tracked per patch, denominated in quanta. Sources are all environmental:

- **Thermal gradients** — hydrothermal vents; the steepest and most reliable gradient available
- **Photic** — surface illumination, cycling with day and season
- **Redox-analogue gradients** — mineral surfaces holding elements at different affinity states
- **Radiogenic** — decay of unstable elements, which also seeds structural mutation

A reaction requiring energy competes for the patch budget. When the budget runs dry, energy-hungry chemistry stalls until conditions change. This is what makes location matter, and why the vents get interesting first.

Energy also needs somewhere to **go**. Every patch is coupled to a heat sink — the deep ocean, the atmosphere — into which degraded energy drains at a rate set by local conditions. Without a sink the region would simply heat until nothing survived, and more fundamentally it would stop being an open system, which per §9.4 is the same as being dead. Sources and sinks together are what make the region a place energy flows *through* rather than a box energy accumulates in.

### 9.3 Autocatalytic set detection

The engine needs to *know* when it has found something, rather than relying on a human noticing.

Each patch maintains a reaction hypergraph. Periodically the engine runs **RAF detection** (Reflexively Autocatalytic and Food-generated set detection, following Hordijk and Steel) — a polynomial-time algorithm that finds subsets of reactions where every reaction is catalysed by a molecule the set itself produces, and every reactant traces back to environmentally available feedstock.

A RAF set is a chemical system that makes itself. Finding one is rigorous, computable, and unambiguous — which makes it the ideal trigger both for promoting a patch to full molecular detail and for the Chronicle to announce that something important just happened.

### 9.4 Why decay is mandatory

Without degradation there is no life. This is not a balancing concern, it is a precondition, and it is worth being exact about why.

**Selection *is* differential persistence.** If nothing degrades, every molecule ever formed persists forever. A superb self-replicator and an inert lump both accumulate without limit, and their relative abundance simply converges on the ratio of their formation rates. Fitness differences never get to express themselves. Selection requires things to be going away at *different rates* — otherwise the word means nothing at all.

**Replication is only worth anything in the presence of destruction.** In a world without decay, making one copy of yourself is as good as making a million: you and your copy both last forever either way. The entire advantage of self-replication is that it outruns loss. Remove loss and you have removed the reason replication matters.

**Finite matter has to be recycled.** The region has a fixed element budget. Without decay every atom eventually ends up locked inside some stable dead-end molecule, and the chemistry grinds to a halt with all its raw material sequestered in junk. A world without rot fills up with garbage and stops.

**Life is a dissipative structure.** It persists *because* matter and energy flow through it. A closed system at equilibrium is by definition dead — nothing changes, by construction. So the region must be genuinely open: energy in from vents and light and out to a heat sink, matter in from geology and out via burial and dilution. That through-flow is what makes persistent structure possible in the first place, and it is why the boundary conditions in §11 matter as much as the chemistry does.

### 9.5 Where decay comes from

Consistent with principle 1, there is no decay timer and there are no hitpoints. Every source of degradation is a physical process the engine is already running for other reasons.

| Source | Mechanism | What it means |
|---|---|---|
| **Spontaneous cleavage** | Every bond carries a thermal breaking probability `k = A·exp(−E_bond/T)`. This is the Cleave class of §9.1 firing with no catalyst present. | Hot places destroy structure quickly, cold places preserve it. Vents are energy-rich *and* corrosive — the first real environmental trade-off. |
| **Solvent attack** | The solvent binds exposed bonds by ordinary complementarity (§8.3) and cleaves them. A monomer's exposure is simply its count of empty lattice neighbours out of twelve. Same mechanism as everything else; nothing new is introduced. | **A molecule's shape determines its lifespan.** Compact folds that bury their backbone survive; sprawling ones are eaten. Folding therefore acquires a survival payoff for free, without us ever rewarding it — and in three dimensions a folded chain has a real interior for that burial to happen in (§8.4). |
| **Radiogenic damage** | Unstable elements decay according to their `decay_rate` (§7.1), destroying their host molecule and damaging neighbours. | One process serves as both the mutation source and a decay source. |
| **Reactive by-products** | Some reactions yield small molecules with unusually *broad* complementarity — they bind almost anything and cleave it. | **Metabolism produces its own poison.** The harder a system runs, the faster it damages itself. |
| **Photic damage** | High-energy input in the surface layer breaks bonds outright. | The photic zone is energy-rich and dangerous at once, which makes depth a genuine niche axis rather than just a coordinate. |
| **Burial and dilution** | Sedimentation and diffusion physically remove material from the active volume. | Not destruction, but functionally identical — it leaves the system. |

Note how little of this is new machinery. Solvent attack is the binding function from §8.3. Radiogenic damage is the `decay_rate` property from §7.1. Reactive by-products are simply molecules whose signatures happen to be broad. Only burial is genuinely new, and it belongs to the world layer regardless.

**Decay also costs nothing per molecule**, which is worth spelling out because the obvious implementation — sweep every molecule every step and roll for each of its bonds — would dominate the loop entirely. It is unnecessary. Thermal cleavage is a Poisson process whose rate depends only on the species and the temperature, so the total propensity for a species is just its count multiplied by a precomputed per-species constant: one more channel in the same scheduler that handles every other reaction. Solvent attack is likewise a pure function of the species, resolved once at intern time (§8.6).

This works *because* damage is implicit (§22.7) — a damaged molecule is simply a different species, so all molecules of a species are interchangeable and none needs individual state. The design is already consistent; the point of saying so is that a well-meant later addition of a per-molecule damage field would quietly turn decay back into a sweep over every molecule in the simulation.

### 9.6 Aging and death, without an aging mechanism

The honest answer to "what causes aging" — and the one this design commits to — is that **aging is not a mechanism at all. It is what damage accumulation looks like when repair is imperfect.**

A protocell has no lifespan counter and no age field. It has three competing rates: damage accruing continuously from §9.5, repair happening only where its own chemistry happens to rebuild what was damaged, and dilution of damage through growth and division.

The crucial asymmetry is that **a RAF set covers some of a compartment's molecules and not others.** Damage to a covered molecule is repaired, because the set remakes it by definition. Damage to an uncovered molecule accumulates with nothing to reverse it. Over time that uncovered damage load rises until it starts interfering with the covered chemistry, the autocatalytic set stops closing, and the compartment dies.

So **death is emergent and precisely defined: a compartment dies when its autocatalytic set no longer closes.** Not a threshold anybody picked — a property we compute, using the same RAF detector from §9.3 that told us it was alive in the first place. Which means the Chronicle can report a death with an actual cause: *which reaction stopped being catalysed, and what stopped catalysing it.*

Aging, in this framing, is just the observation that the damage a system cannot repair outpaces the damage it can. Which is, as far as anyone can currently tell, roughly what aging actually is.

### 9.7 Damage segregation — a prediction, not a feature

Worth flagging, because it is the class of result that would justify the whole project.

When a compartment divides, damaged material partitions stochastically. Damaged molecules tend to be larger, stickier, and slower to diffuse, so there is no particular reason the split should be even. If it turns out *uneven*, one daughter inherits most of the junk while the other gets something close to a fresh start — and the lineage as a whole becomes effectively immortal by concentrating its damage into a sacrificial line.

That is asymmetric damage segregation, and it is how single-celled organisms are understood to handle aging. **Borbax does not implement it.** If it appears, it appears because sticky things do not diffuse evenly. Detecting it is straightforward — compare damage load between sibling lineages — and it is exactly the kind of thing the Chronicle should shout about.

---

## 10. Layer 3 — Compartments and protocells

### 10.1 Membranes

Molecules with strongly **asymmetric** signatures — one flank strongly positive in affinity, the opposite flank strongly negative — will preferentially align side-by-side with their like flanks together. At sufficient concentration this produces sheets, and sheets close into vesicles because an open edge is energetically expensive.

None of this is special-cased. It is §8.3 applied to many molecules at once. Amphiphile-analogues are not a molecule type we define; they are a region of signature space, and the sim discovers which of its molecules live there.

### 10.2 Compartments

A closed vesicle is promoted to a **compartment**: a distinct reaction volume with its own concentration state, its own energy budget, and selective permeability. What crosses the membrane is decided by running the binding test between the molecule and the membrane's pore signatures. A compartment that develops a pore complementary to a useful feedstock has, in effect, evolved a transporter.

The nested-compartment topology here is close to what P systems (§2.1) formalise, and that formalism is worth borrowing from when implementing transport and containment.

### 10.3 What counts as a protocell

A protocell is a name we give to a compartment exhibiting all four of:

1. **Boundary** — a closed membrane persisting beyond a threshold duration
2. **Metabolism** — a RAF set operating inside it
3. **Growth** — net accumulation of membrane material from internal chemistry
4. **Division** — mechanical fission driven by surface-area-to-volume stress

Division is physics. A growing vesicle becomes unstable past a critical ratio and splits, partitioning its contents stochastically. Daughters inherit an imperfect sample of the parent's chemistry — and that imperfect sampling is the origin of heritable variation, which is the door into V2.

And a protocell stops being one the moment any of the four fails. In practice that is almost always metabolism: per §9.6, accumulated damage the RAF set cannot repair eventually stops the set closing, and the compartment dies with a computable cause. Sometimes it is the boundary instead — a membrane degrading faster than its chemistry replaces it, and the compartment simply comes apart. Both are ordinary consequences of §9.5 rather than special cases, and both are worth reporting: **the census of what died and why is at least as informative as the census of what lived**, because it is the clearest signal of what the environment is actually selecting for.

---

## 11. Layer 4 — World generation

`world_gen` takes a validated universe plus a `world_seed` and produces a specific place.

**V1 region: 16 km × 16 km (256 km²)**, as a 512 × 512 patch grid at ~31 m per patch, with 4 vertical layers (surface, photic, deep, sediment) — roughly 1M cells.

The region is chosen to be maximally interesting for prebiotic chemistry: a **volcanic coastal shelf** combining a hydrothermal vent field, a tidal flat that wets and dries, mineral-rich porous rock, and a shallow photic zone. That gives every energy source and every ratchet in one place.

Generated per world:

- **Terrain** — bathymetry, vent placement, mineral outcrops, porosity
- **Element abundances** — which of the universe's elements are common here, and where
- **Climate** — temperature fields, currents, seasonal and tidal cycles
- **Initial inventory** — the small-molecule feedstock produced by geology alone, computed by running the beaker battery under local conditions

**The ratchets** are environmental features, not shortcuts:

| Ratchet | Effect |
|---|---|
| Wet/dry cycling | Tidal flats concentrate solutes to the point where condensation becomes favourable |
| Mineral surfaces | Catalytic classes lower activation barriers and hold molecules in orientation |
| Thermal gradients | Vent walls provide steep, sustained energy differentials |
| Confined pores | Rock pores prevent dispersal, keeping reaction products together long enough to matter |
| Convection | Vent circulation cycles material through temperature extremes repeatedly |

Every one of these is a real feature of the environment. None of them knows what a replicator is.

---

## 12. Layer 5 — The multi-scale engine

This is the hardest engineering in the project and the reason an overnight run can cover deep time.

### 12.1 Three tiers

| Tier | Patches | Representation | Timestep | Cost |
|---|---|---|---|---|
| **Bulk** | ~10⁶ | Concentration fields, reaction-diffusion | years | very low |
| **Stochastic** | ~10³ | Discrete molecule counts, tau-leaping | seconds–minutes | moderate |
| **Molecular** | ~10 | Every molecule an object with position and shape | microseconds | high |

Bulk patches store the top 32 species by abundance plus an aggregate "other" pool — bounded memory (~256 bytes/cell, ~256 MB total) regardless of how baroque the chemistry becomes.

### 12.2 Promotion and demotion

A patch is promoted when the **novelty detector** fires on any of:

- A molecular species appearing that has never been seen in this run
- Reaction-network complexity rising sharply
- A RAF candidate partially matching
- Concentration of any species crossing a threshold
- A membrane-forming signature region becoming populated

Demotion happens when a patch goes quiet for a sustained interval. Promotion carries state up by sampling the concentration field into discrete molecules; demotion carries it down by binning molecules back into the field. Both directions are lossy, which is acceptable and must be *tested* — see §16.

This is principle 3 made concrete: the engine continuously moves its expensive attention to wherever something unprecedented is happening.

### 12.3 Event-driven time

Borbax never steps through uneventful time.

Each patch reports the time of its next scheduled event — a reaction firing, a tidal transition, a climate change, a decay. The global scheduler holds these in a priority queue and jumps directly to the earliest. Quiet regions cost nothing but a queue entry.

Two further accelerations:

- **Quasi-equilibrium jumping.** When a bulk patch's chemistry reaches a fixed point, stop integrating it. Solve for the steady state analytically and park it until an external condition changes.
- **Rare-event sampling.** The interesting events are rare by definition. Rather than waiting through 10⁶ years for an improbable condensation, sample the waiting time directly from the rate distribution and jump to it. This is Gillespie's next-reaction method, applied across the whole hierarchy.

**Consequence, stated honestly:** time compression is *variable and emergent*, not a fixed guarantee. A quiet epoch may pass at 10⁶ simulated years per wall-second; an active protocell may run at near-realtime. This is a feature — the pacing of the visualisation automatically matches where the interest is — but it means the 8-hour target cannot be guaranteed analytically. It is achieved by a governor.

### 12.4 The compute governor

The governor is a control loop targeting **10⁷–10⁸ simulated years within an 8-hour run**. It monitors the achieved rate and adjusts promotion thresholds, the number of concurrently molecular patches, and folding cache aggressiveness to steer toward the target.

It is deliberately conservative about one thing: it will never demote a patch containing an active RAF set to hit a time target. Reaching the target is worthless if we bulldoze the result on the way.

---

## 13. Layer 6 — Determinism and the time machine

### 13.1 Determinism

Bit-for-bit reproducibility given `(universe_seed, world_seed, config_hash)`.

**`config_hash` is the physics version §13.3 writes as `@n`, under a second
name, and nothing else in configuration may reach the arithmetic.** That
resolution is forced rather than chosen: §13.3 promises that sharing
`U-7F3A21C9@1 / W-0004` gets the recipient a bit-identical planet, and that
promise fails the moment a *fourth* value can change results while sitting
outside the address. Either a setting changes the physics — in which case it
mints a new version, is carried in the address, and a recipient with a
different one can see that they differ — or it cannot change results, in which
case it does not belong in a reproducibility tuple at all.

So configuration splits in two, and the split is a rule rather than a
convention:

- **Runtime knobs** — keyframe interval, output paths, tier-promotion
  thresholds, thread count. These must not reach any arithmetic that produces a
  simulation value. A knob that would is not a knob.
- **Physics** — anything that changes what the laws *are*. Sweeping such a
  value is legitimate (§22.2 sweeps and then **locks** the geodesic resolution),
  but the locked value becomes a constant under a version, not a field that
  travels with a run.

The failure this closes: a keyframe-restore check on `config_hash` only fires
on restore. Two people sharing a seed pair and starting fresh never restore
anything, so nothing would have caught a physics-affecting knob that differed
between them. Putting physics in the address rather than in a side-channel is
what makes the §13.3 promise true instead of merely stated.

This requires discipline the engine must enforce structurally:

- **Counter-based RNG**, with an independent stream derived per subsystem, per patch, per simulated time index. No shared mutable RNG state, so results never depend on thread scheduling. Streams are deliberately not `Copy`, so an accidental duplication — which would silently replay the same sequence twice — has to be written explicitly.
- **No wall-clock, no address-dependent iteration.** Hash maps are never iterated in a result-affecting path; sort into a `Vec` first. Every sort on a float key carries an ID tie-break, or equal keys reorder under a different input permutation.
- **Deterministic parallel reduction.** Fold into per-chunk partials indexed by chunk ID, then combine in index order. Never a scheduler-ordered `.sum()` or `.reduce()`.
- **No platform transcendentals.** `exp`, `ln`, `sin`, `cos` and `powf` are *not* specified by IEEE-754 to be correctly rounded, and implementations genuinely differ between Apple's libm, glibc, and glibc versions. Every transcendental routes through a single vendored pure-Rust implementation, or through a precomputed table. The four arithmetic operations and `sqrt` *are* exactly specified by IEEE-754 and can stay native.
- **Fixed summation order is physics.** `-(w₁·Σa + w₂·Σb)` and `Σ(-w₁·a − w₂·b)` are different values in floating point. Whichever shape is chosen is pinned before any golden exists, and commented as load-bearing so nobody later "tidies" it.
- **Conserved quantities are integers.** Mass is fixed-point, not `f64`, so conservation is exact by construction rather than true within a tolerance. §9's requirement that products' mass equals reactants' mass *exactly* is otherwise not satisfiable by naive floating-point summation.
- **Single canonical CPU code path.** GPU acceleration, when it lands, is an optional non-canonical path — see §13.5.

### 13.2 Keyframes and scrubbing

Full state snapshots every K simulated years, compressed, plus a complete event journal.

Scrubbing to time T loads the nearest keyframe at or before T and deterministically replays forward. Storage stays modest, seeks stay fast, and — because replay is exact — the rewound state is genuinely the state that occurred, not an approximation of it.

### 13.3 Branching

Any keyframe can be forked with a new sub-seed, producing a divergent timeline sharing history up to the fork. This is the substrate that later makes intervention possible: "what if the vent had never cooled?" is a branch, and the two timelines can be compared side by side. V1 builds the mechanism; V5 builds the controls.

It is also the mechanism behind neutral shadow runs (§15.3), which is a second reason to get it right early.

### 13.4 What "deterministic" means across machines

Borbax commits to the strong version: **the same seed pair produces a bit-identical world on macOS, Windows and Linux, on both aarch64 and x86-64.** Sharing `U-7F3A21C9 / W-0004` has to mean something, and it only does if the recipient gets the same planet.

This is achievable, but only because of a distinction that is easy to miss: **cross-platform APIs make code *run* everywhere; they do not make it produce *identical results* everywhere.** Those are different problems with different solutions, and conflating them is the standard way projects discover the issue too late.

The good news is that most of the arithmetic is already portable. IEEE-754 specifies addition, subtraction, multiplication, division and square root *exactly* — a conforming implementation must return the correctly rounded result, so those operations are bit-identical on every platform we care about. Rust also does not perform floating-point contraction by default, so `a*b + c` will not silently become a fused multiply-add on one architecture and not another. That is a significant hazard we simply do not have.

What is left is a short, closed list:

| Hazard | Resolution |
|---|---|
| Transcendentals differ between libm implementations | Vendored pure-Rust implementations, or precomputed tables (§13.1) |
| Parallel reduction order varies with scheduling | Indexed partials, combined in fixed order |
| Hash map iteration order | Never iterated in a result-affecting path |
| Float accumulation order | Pinned and commented before goldens exist |
| Exact mass conservation | Fixed-point, not floating-point |

The games industry solved this problem decades ago and its answer is worth borrowing wholesale. Lockstep real-time strategy games, rollback-netcode fighting games, and replay-driven simulations like Factorio all require every machine to reach an identical state from identical input — and they achieve it by keeping simulation state in **fixed-point integer arithmetic on the CPU**, with the GPU confined to rendering. Integer arithmetic is exact on every processor ever made, which makes the problem disappear rather than requiring it to be managed.

Borbax takes the hybrid: **fixed-point for conserved and counted quantities** where exactness is a stated requirement, and **strict IEEE-754 `f64` for the geometry**, where the operations involved are already exactly specified. That keeps the shape mathematics readable — and readability matters a great deal in the one part of the codebase that carries the project's central idea — without giving up exactness where exactness was promised.

### 13.5 GPU: what it can and cannot do

**The right tool is [wgpu](https://wgpu.rs).** It implements the WebGPU standard and compiles one WGSL compute-shader source to Metal on macOS, and DirectX 12 or Vulkan on Windows — plus a browser path we may want later for the viewer. It is mature, it is the established choice in the Rust ecosystem, and it means one shader codebase rather than a Metal backend and a DirectX backend maintained in parallel.

What it explicitly does *not* provide is bit-identical results. WebGPU permits implementations latitude in the precision of many operations, drivers differ between vendors and between driver versions, and reduction order inside a GPU dispatch is not something the API lets us pin. **A GPU result is therefore never canonical**, and no amount of care changes that.

Rather than fight this, the architecture makes it irrelevant by giving the GPU a job where it does not matter:

> **The GPU is a search accelerator, not a result generator.**
>
> It screens — running many universes and worlds at low fidelity to find the ones worth looking at. Screening does not need reproducibility; it needs throughput, and a large discrete GPU has an enormous amount of it.
>
> When screening finds something interesting, **that seed is re-run on the canonical CPU path**, which produces the real, shareable, reproducible result.

So a powerful Windows machine genuinely earns its keep — it explores the space of universes far faster than an M1 Pro can — without any of its output entering the reproducibility story. The two paths have different jobs and neither is a degraded version of the other.

Within a run, the GPU-suitable work is the bulk-tier reaction-diffusion stencil (embarrassingly parallel over ~10⁶ patches) and batch folding (thousands of independent anneals). Canonicalisation, RAF detection, and the event scheduler are branchy and serial, and stay on the CPU regardless of hardware.

### 13.6 Enforcement — the cross-platform golden matrix

None of the above is worth anything as an intention. It is worth something as a test.

CI runs the golden suite (§16) on **macOS/aarch64, Windows/x86-64 and Linux/x86-64**, and asserts the state hashes are identical across all three. A portability regression then fails the build on the commit that introduced it, rather than being discovered months later by the first person who tries to open a shared seed.

This matters more than it might appear. Every hazard in §13.4 is invisible locally — the code runs, the tests pass, the numbers look reasonable — and only manifests on a different machine. The matrix is the only mechanism that makes those failures visible at the moment they are cheap to fix.

---

## 14. Layer 7 — The viewer

**A native desktop application: `egui` for the interface, `wgpu` for the 3D scene beneath it.** One language, one build, and the viewer reads simulation state directly rather than serialising it across a boundary sixty times a second.

### 14.0 Why native, and why not a game engine

The original plan was a browser application over a local WebSocket. That was a reflexive choice rather than a reasoned one, and three things undermine it:

- The geometry is three-dimensional (§8.2), so the browser path needs WebGL regardless — it buys no rendering capability the native path lacks.
- `wgpu` is already required for GPU screening (§13.5), so it is in the dependency tree either way.
- A web viewer adds a second language, a second toolchain, and a wire protocol to design and version — for a program running on the same machine as the simulation.

There is also an audience argument, which matters here more than it usually would: a double-clickable application is meaningfully easier for a child than "start a server, then open localhost".

Going native does not close the browser door. **egui compiles to WebAssembly**, so a browser build remains available from the same source if it is ever wanted.

**Game engines** — Unity, Unreal, Godot — were considered and rejected for V1, but not for the obvious reason. They are excellent at what V5 needs. The problem is the FFI boundary: the simulation must stay Rust for determinism (§13.1) and throughput (§17), so a C# or C++ engine means marshalling state across a C ABI every frame. That is the same objection as the WebSocket, with a harder debugging story and a build system that wants to own the repository.

**Bevy** avoids that entirely, being Rust and built on wgpu — and it is the right thing to revisit for the Wonder skin at V5. It is wrong for V1 because of what V1's viewer actually is. Zoom levels, a depth dial, the Chronicle, population graphs, lineage trees, reaction-network diagrams, a molecule inspector: this is a *tool*, not a game. An engine's strengths — lighting, particles, terrain, physics, asset pipelines — go almost entirely unused, while its comparative weakness at dense two-dimensional interface work applies to nearly every screen.

The sequencing is what makes this safe. **Bevy is built on wgpu**, so renderer code written against wgpu now carries forward into Bevy later. Choosing Unity now would instead create a boundary that never goes away.

### 14.1 Zoom levels

Continuous navigation across four scales, each rendering the tier beneath it:

```
  REGION  ─────►  PATCH  ─────►  MOLECULAR  ─────►  INSPECTOR
  16 km            31 m          individual         one molecule:
  terrain,         local         molecules          graph, signature,
  temperature,     chemistry,    moving, binding    fold, reaction
  where the        gradients,    and reacting       history, lineage
  action is        compartments
```

### 14.2 Two skins, one data model

**Lab skin** ships first: schematic, high-contrast, precise. Molecules as clear diagrams, signatures as unfolded geodesic nets, reaction networks as graphs. This is the debugging tool, and it is what makes the simulation trustworthy — you can see exactly what the engine believes.

**Wonder skin** arrives once the simulation is believed: painterly, atmospheric, with depth and light and weather, and a microscopic view that feels like looking down a real instrument. Same data, different presentation. Lenia and Flow-Lenia (§2.3) are the bar.

This is the point at which a game engine genuinely earns its keep, and where **Bevy** should be revisited (§14.0). Because it is built on wgpu, the renderer written for the Lab skin carries forward rather than being thrown away.

### 14.3 The depth dial

One control, three registers, applied to vocabulary, chart density, and how much is explained inline:

| Register | Character |
|---|---|
| **Storybook** | Plain language, minimal numbers, colour and motion carry the meaning. *"This little shape just learned to make copies of itself."* |
| **Explorer** | Real vocabulary, always explained in context. Simple charts. Every claim expandable into why it happened. |
| **Scientist** | Full metrics, population graphs, reaction networks, lineage trees, arbitrary molecule inspection. |

This is what "build it to scale" means in practice: one interface that neither condescends to a nine-year-old nor expires when she turns fourteen. It also happens to be exactly what the developer needs.

### 14.4 The Chronicle

**The most important feature in the viewer.** An eight-hour run is unwatchable; a story about an eight-hour run is not.

The Chronicle is an automatically generated narrative log. Every significant event gets an entry with its timestamp, its location, a plain-language description at the current depth register, and a link that jumps the time machine straight to that moment.

Milestone events include: first polymer above threshold length; first stable fold; first enclosed cavity; first two-cavity enzyme; first membrane sheet; first closed vesicle; first RAF set; first compartment growth; first division; first lineage divergence.

Deaths are milestones too, and are reported with their computed cause (§9.6): the first compartment death, the first lineage to outlive its founder, the first repair of damage rather than mere replacement, and — if it happens — the first evidence of asymmetric damage segregation (§9.7). A chronicle of only the triumphs would misrepresent what is actually going on, which is overwhelmingly a story about things falling apart slightly slower than they are built.

The Chronicle turns "I left it running overnight" into "come and read what happened while we were asleep, and then let's go and watch the good bits."

### 14.5 The V0 shape inspector

The viewer proper is a V1 deliverable, but a shape-based chemistry cannot be developed blind. If you cannot see that two molecules interlock, you cannot distinguish a binding bug from a chemistry that is simply not very interesting — and those two failures look identical from the outside. So V0 ships a deliberately minimal precursor.

The Rust core emits **static SVG** directly:

- a molecule as a graph diagram, plus an orthographic projection of its 3D embedding
- its signature as an **unfolded geodesic net** — the icosahedron flattened into triangles, each sample direction coloured by extent and affinity. This is the map-projection trick, and it makes a spherical function legible on a flat page without any interaction
- a folded polymer in isometric projection with depth cueing, cavities highlighted, and buried monomers shaded by how enclosed they are
- a binding pair in the winning orientation, with per-direction contributions annotated so you can see *which* part of the fit carried the score

A contact sheet — a plain index page of the emitted SVGs — collects these for a run.

Three dimensions do make this harder, since a rotatable view is V1 work. But static projections chosen well are often *clearer* than interactive ones for the specific question "why did these two bind?", because a fixed viewpoint can be picked to show exactly the face that mattered.

Two reasons for SVG specifically:

- **It is testable.** Deterministic string output means renders can be golden-tested in CI, so a rendering regression is caught like any other regression. A GPU renderer would need image diffing and a windowing system to achieve the same thing.
- **It has no build step and no window.** An SVG opens in a browser, an editor, or a diff, and it can be produced headlessly on a CI runner. During V0 the render loop needs to be as short as the simulation loop.

**SVG does not become redundant when the interactive viewer arrives.** The two have different jobs and both keep them:

| | SVG (§14.5) | egui + wgpu (§14.0) |
|---|---|---|
| Runs headless in CI | yes | no |
| Golden-testable | yes — the output is a string | no, needs image diffing |
| Embeddable in docs and diffs | yes | no |
| Interactive, rotatable, live | no | yes |

The duplication is smaller than it appears, because **the geometry is computed once and rendered twice**. Projection, layout, cavity extraction and colour mapping live in `borbax-geometry`; SVG and wgpu are two backends consuming the same computed result. What SVG costs on top is a few hundred lines of string emission — a cheap price for having the renderers under test at all.

---

## 15. Success criteria

### 15.1 V1 is done when

Given a validated universe and a world seed, an 8-hour run on the target machine produces, unprompted:

1. A closed compartment
2. Containing a detected RAF set
3. That grows from its own internal chemistry
4. And divides
5. Whose descendants persist for ≥ 20 generations *against an active decay background* — persistence has to be continuously earned, not merely undisturbed (§9.4)
6. With measurable heritable variation between sibling lineages

All six, from chemistry alone, with no code path anywhere in the engine that mentions replication — and each verified against the neutral shadow (§15.3), not merely observed.

### 15.2 Instrumented metrics

Tracked continuously and graphable across a run. Note that **activity and novelty are deliberately separated**, per the Bedau class 3b failure mode in §2.7 — a system can generate unbounded activity while producing nothing new, and reading one as a proxy for the other is exactly how that goes unnoticed.

**Raw counts are not sufficient to detect it.** A species that appears once and dies counts the same as one that persists for a million years, so "distinct species" and "novel species rate" are proxies that miss the failure mode they exist to catch. The Bedau–Packard statistics are persistence-weighted by construction, and that weighting *is* the mechanism: per-species cumulative existence counters `aᵢ(t)`, diversity `D(t)`, cumulative activity `A(t) = Σaᵢ`, mean cumulative activity `Ā = A/D`, and new activity `A_new(t)`. Class 3b is exactly `A` unbounded with `D` bounded, which is only visible if both are computed.

| Metric | Family | What it tells us |
|---|---|---|
| Total reaction events | Activity | Is anything happening at all? |
| Distinct species count | Activity | Is chemical diversity growing or saturating? |
| **Novel species rate** | **Novelty** | Is the run producing genuinely unprecedented things? |
| **Shape-space novelty** | **Novelty** | Distance from the current shape-occupancy histogram to every prior one — adapted from ASAL's frame-similarity metric (Sakana AI). A direct runtime answer to "is this still interesting?" |
| Max polymer length | Complexity | Is complexity accumulating? |
| Reaction network diameter | Complexity | Is the chemistry becoming interconnected? |
| RAF set count and size | Organisation | Are self-sustaining systems appearing? |
| Compartment census | Organisation | Population dynamics of proto-life |
| **Turnover rate** | **Decay** | Mass destroyed per unit time as a fraction of mass created. The single clearest indicator that the system is a flow rather than a heap (§9.4) |
| **Mean molecular lifespan** | **Decay** | Per species class. Shows which shapes are surviving and lets us confirm the folding-protects-you prediction of §9.5 |
| **Compartment lifespan distribution** | **Decay** | Are protocells getting measurably better at staying alive? |
| **Damage load** | **Decay** | Unrepaired damage per compartment. This is the direct measurement of aging (§9.6) |
| **Cause-of-death census** | **Decay** | Which reaction stopped closing, and why. Often more informative about what the environment selects for than the survivors are |
| Simulated years per wall-second | Health | Governor performance |

Per Hintze's critique (§2.7), **every one of these is gameable in isolation.** None is trusted alone, and all are read against the shadow.

### 15.3 The neutral shadow

Adopted as a **V1 deliverable**, not an analysis tool added later.

Every run may be paired with a *shadow run*: the same universe, the same world, the same seeds, forked from the same keyframe — with selection switched off. Molecules still react, but differential persistence is disabled, so any structure appearing in the shadow is drift rather than selection.

Every claim of the form "look what evolved" is then read as **run minus shadow**. Without this, activity statistics are uninterpretable, and the entire success criterion in §15.1 reduces to a person looking at a graph and feeling optimistic.

**"Selection switched off" means decay rates are *equalised*, not disabled.** This distinction is easy to lose and gets the experiment wrong in a way that still produces plausible-looking output. Borbax has no fitness function — selection *is* differential persistence (§9.4) — so a shadow with decay turned off is not a neutral control, it is a different universe in which nothing dies. The shadow must apply a single common decay rate to every species, set so the total removal *flux* matches the focal run at the fork keyframe. Otherwise the two runs differ in mass balance and any diversity difference is explained by that rather than by selection.

Two honest limits. First, this is a **drift control on persistence only**: catalysis produces differential *formation* rates, which is also selection in Bedau's sense, and equalising it would destroy the chemistry we are trying to observe. Second, it does not answer the question that matters most about a RAF. For that, the standard control in this literature is different — **randomise which molecule catalyses which reaction**, holding catalysis density fixed, and ask how often a RAF still appears. That distinguishes "a RAF appeared because of the shape chemistry" from "any network this dense has one", and the persistence shadow cannot.

This is cheap to implement because the branching machinery (§13.3) already exists for the time machine, and it is the difference between a result and an anecdote.

### 15.4 Plateau detection, done properly

The most likely way this project fails is not a crash. It is becoming boring, silently, six hours in.

But per the boundedness illusion (§2.7), **a flattening curve is not evidence of a plateau.** Judged by eye, even genuinely unbounded processes look like they are levelling off; the same reasoning applied to real long-term evolution experiments turns out to be wrong once competing models are actually fitted.

So Borbax never declares stagnation from a graph. It fits several candidate models — saturating, linear, power-law — to the novelty trajectory, compares them by **AICc and BIC reported together with the fitted parameters**, and reports which fits best along with its confidence. A run flags itself as stalled only when a saturating model genuinely beats an unbounded one, and even then the finding is checked against the shadow.

Three constraints on the fit, each of which is a way the comparison silently becomes meaningless:

- **Fit every model on the same scale by the same method.** Fitting the power law on log–log axes and the saturating model on linear axes is not a comparison of models, it is a comparison of axes.
- **Fit increments, not cumulative totals.** A cumulative curve has near-perfectly autocorrelated residuals, and an information criterion computed on it under an i.i.d. error assumption will confidently "prove" whichever model was fitted.
- **Hold out the tail.** The boundedness illusion is a *projection* failure, so fit on the first half and score predictive error on the second. That tests the exact claim being made, and is more convincing than any information criterion.

---

## 16. Testing strategy

A simulation this stochastic is untestable by conventional means unless the strategy is designed in from the start.

**Property tests.** Mass and energy conservation across every reaction class, over randomly generated molecules. Canonical form stability — graph relabelling must not change the canonical hash. Signature stability — a small graph change must produce a small signature change, which is the property that makes evolution possible and would fail silently otherwise.

**Determinism tests.** Same seeds produce identical state hashes at fixed checkpoints. Run under varying thread counts to catch scheduler-order dependence, which is the bug class most likely to appear and least likely to be noticed.

**Golden runs, on every platform.** A handful of seeds with recorded state hashes at checkpoints, executed in CI on **macOS/aarch64, Windows/x86-64 and Linux/x86-64**, asserting the hashes are identical across all three (§13.6). Any unintended change to physics fails the build loudly; any *portability* regression fails it on the commit that caused it, rather than months later when someone opens a shared seed. Running the goldens on a single platform would pass happily while the cross-platform promise quietly rotted.

**Cross-tier agreement.** The critical validation of the whole multi-scale approach: run an identical scenario at bulk, stochastic, and molecular fidelity, and confirm the statistics agree within tolerance in the overlap regime. If they diverge, the promotion and demotion mapping is wrong and every result above it is suspect.

**Decay and turnover tests.** Cleavage must conserve mass — fragments sum to the parent — and turnover rate must respond monotonically to temperature. More interestingly, the beaker harness must be able to *demonstrate* the central claim of §9.4: with decay disabled, a run and its neutral shadow should become statistically indistinguishable, because without differential persistence there is no selection to detect. That is the direct experimental test of "no death, no life", it takes seconds to run, and it would catch a whole class of subtle errors in which decay is nominally implemented but not actually biting.

**Folding-map property tests.** Verify the three properties from §2.5 — neutral network connectivity, approximate shape-space covering, and plastogenetic congruence — as measurable assertions in the beaker harness. These are what make the chemistry evolvable, so they are tested, not assumed.

**The beaker harness.** A single well-mixed patch runnable in milliseconds. This is the tuning loop for chemistry parameters and it is worth building first — it converts a chemistry-design iteration from hours to seconds.

**Emergence regression suite.** Nightly: does universe X still reach autocatalysis by year Y? Guards against a well-meaning tuning change quietly destroying the emergent behaviour that is the entire point. Slow, allowed to be flaky in timing, judged on distributions rather than exact outcomes.

**Golden render tests.** SVG output (§14.5) is a deterministic string, so renders are golden-tested like any other output. This is cheap, it catches rendering regressions that would otherwise be discovered by a human noticing a diagram looks odd, and it means the diagrams stay trustworthy as the thing we debug the chemistry *with*.

**Fiction-guarantee checks.** CI enforcement of §5: no data files under chemistry crates, blocklist coverage verified, unit types enforced.

---

## 17. Performance budget

Target machine: Apple M1 Pro — 8 performance cores, 2 efficiency cores, 16-core GPU, 32 GB unified memory.

| Resource | Allocation |
|---|---|
| Worker threads | 8, on performance cores; efficiency cores left for OS and viewer |
| Bulk tier | ~10⁵ patch-updates/sec/core |
| Molecular tier | ~10⁵ molecule-steps/sec/core |
| Fold cache | 2 GB, target hit rate > 95% at steady state |
| Binding memo | Dense table over the most abundant species (§8.6). This, not the kernel, is what makes 3D binding affordable |
| Fold workspace | ~512 KB per worker thread for the FCC occupancy grid, allocated once and reset by trail rather than cleared |
| Bulk field state | ~256 MB |
| Total working set | < 8 GB, leaving headroom for viewer and keyframe buffers |
| Keyframe storage | < 20 GB per 8-hour run, compressed |

**GPU is explicitly out of scope for V1**, and when it arrives it arrives as a *screening* path rather than an acceleration of the canonical one — see §13.5 for why that distinction is the whole architecture rather than a caveat. The tool will be [wgpu](https://wgpu.rs), giving one WGSL source across Metal, DirectX 12 and Vulkan.

The reason to defer it is not difficulty. It is that a GPU screening path is an optimisation of a system that must first be correct and reproducible, and building it earlier would mean tuning chemistry against results we cannot reproduce.

**Borbax runs on macOS, Windows and Linux from V0**, on both aarch64 and x86-64, with identical results everywhere (§13.4). This is a portability requirement on the CPU path, enforced by the CI matrix in §13.6 — not a GPU concern.

---

## 18. Technology and repository layout

### 18.1 Toolchain

**proto is the sole toolchain manager**, with `.prototools` at the repo root pinning Rust. This matches the convention already established across the other repositories here, so there is nothing new to learn and no second mechanism to keep in sync.

There is also a reason specific to this project, and it is a strong one. **The toolchain pin is part of the reproducibility contract, not a developer convenience.** Determinism is a hard requirement (§13.1) and golden-run state hashes (§16) are only meaningful relative to a fixed compiler — a rustc bump can change floating-point codegen, autovectorisation, or intrinsic selection, any of which will silently change state hashes without changing a line of our code. Pinning makes that a deliberate, visible event: **a toolchain bump requires regenerating goldens, and the regeneration is reviewed like any other change to physics.** Without the pin, a background dependency update would look exactly like a physics regression, and we would waste a day finding that out.

```toml
# .prototools — toolchain pins. See §13.1: these are part of the
# determinism contract, not just developer convenience. Bumping any
# of them invalidates the golden-run hashes in §16 and requires
# deliberate regeneration.

# Current stable, released 2026-07-14.
# renovate: datasource=github-releases depName=rust-lang/rust
rust = "1.97.1"

[settings]
# Auto-install a pinned version when it isn't present yet, so a
# Renovate bump self-installs at job time rather than needing a
# manual step. Matches the convention in the other repos here.
auto-install = true
```

**Rust is the only entry, and that is a consequence worth noticing.** The native viewer decision in §14.0 removed Node, pnpm and TypeScript from the project entirely — one language, one toolchain, one dependency graph, one test command. A cross-language build was never going to be the hardest thing here, but it was pure overhead, and it is now simply absent.

### 18.2 Language choice — Rust, and why not C++

Rust for the core, and the honest comparison against C++ is worth recording because it is the kind of decision that gets re-argued later.

**The decisive argument is determinism under parallelism.** The hardest requirement in this project is bit-for-bit reproducibility across machines while running on eight cores (§13.4). In C++ a data race is undefined behaviour: it can silently produce a divergence that appears on one platform and not another, intermittently. The CI matrix (§13.6) would catch *that* something diverged, but hunting a race that only manifests on one architecture is genuinely miserable work. Rust makes that entire bug class structurally impossible. For a project whose central promise is reproducibility, that is not a preference.

Secondary but real: the invented-unit newtypes of G4 cost nothing at runtime and are enforced by the compiler; Cargo makes the toolchain pin — which is *part* of the determinism contract (§18.1) — a one-line fact rather than a CMake and package-manager project of its own; and `wgpu` (§13.5) gives one cross-platform GPU path where C++ would mean Vulkan and Metal maintained separately, or a translation layer.

**The fair concern is reinventing wheels, and it deserves a straight answer.** The honest position is that for this project the wheel mostly does not exist in either language:

| Need | C/C++ | Rust | Verdict |
|---|---|---|---|
| Graph canonicalisation | **nauty / Traces** — the gold standard | no mature equivalent | **C++ wins** — the one real gap |
| Cross-platform GPU | Vulkan + Metal, or MoltenVK | `wgpu`, one source | Rust |
| Immediate-mode UI | Dear ImGui, very mature | `egui` | roughly level |
| Deterministic parallelism | OpenMP / TBB, races are UB | `rayon`, races rejected at compile time | Rust, decisively |
| Build and dependencies | CMake + vcpkg/conan | Cargo | Rust, decisively |
| Serialisation | protobuf / cereal | `serde` | Rust |
| Generated periodic table | — | — | **nothing exists** |
| Shape-signature binding | — | — | **nothing exists** |
| Lattice folding on invented monomers | — | — | **nothing exists** |
| RAF detection | — | — | **nothing exists** |

The bottom four rows are most of the work, and they are novel *by construction* — inventing the chemistry is the premise (§1.1), so no library in any language implements it. What an ecosystem supplies here is the periphery: GPU, interface, parallelism, serialisation, build. Rust's periphery is at least as good as C++'s across that set and better in several places.

**On the one genuine gap**: nauty is C with a stable ABI and can be called from Rust directly, so even that does not require choosing a language. It is also unlikely to be needed — molecules are capped at twelve atoms with element and bond-order labels, and colour refinement alone discretises the large majority of such graphs, leaving the expensive search a rarely-taken path. If canonicalisation later proves a bottleneck or a correctness headache, **binding to nauty is the fix, not a rewrite** (`nauty-Traces-sys`, actively maintained; the petgraph-flavoured `graph-canon` wrapper is stale). That keeps the decision bounded and reversible.

`egui` and `wgpu` for the viewer (§14.0), which keeps the whole project in one language and lets the interface read simulation state directly. `wgpu` is needed regardless for GPU screening (§13.5), so the viewer's renderer and the screening path share a stack rather than duplicating one.

### 18.3 Repository layout

```
borbax/
├── .prototools              Toolchain pin — rust only (§18.1)
├── crates/
│   ├── borbax-units/        Newtype units, and det_math — see below
│   │                        (thermals, quanta, spans, world-years, fixed-point mass)
│   ├── borbax-rng/          Counter-based deterministic RNG streams
│   ├── borbax-universe/     Physics generation, periodic table, beaker battery
│   ├── borbax-molecule/     Graphs, canonicalisation, signatures, folding
│   ├── borbax-reaction/     Reaction classes, rates, energy, decay, RAF detection
│   ├── borbax-world/        Terrain, climate, abundances, ratchets
│   ├── borbax-sim/          Multi-scale engine, scheduler, governor, LOD
│   ├── borbax-metrics/      Activity/novelty metrics, shadow comparison, model fitting
│   ├── borbax-record/       Keyframes, event journal, branching
│   ├── borbax-geometry/     Projection, layout, colour mapping — shared by both renderers
│   ├── borbax-render/       Static SVG emitters, golden-testable (§14.5)
│   ├── borbax-ui/           egui interface + wgpu 3D scene (§14.0)
│   └── borbax-cli/          Headless runner, beaker harness, tooling
└── docs/
```

The crate split follows principle: each has one clear purpose, a defined interface, and can be understood and tested without reading the others. `borbax-molecule` in particular must be comprehensible in isolation — it is where the conceptual weight of the project sits, and where most of the iteration will happen.

**`det_math` belongs in `borbax-units`, at the bottom of the chain.** Every transcendental in the project routes through it (§13.1), and its callers include `borbax-rng`'s normal-variate draw and `borbax-universe`'s mass generation — both of which sit *below* `borbax-molecule`. Placing it any higher makes the one file that exists to guarantee portability unreachable from two of the crates that need it, and the only way out would be a dependency edge that inverts the crate order.

---

## 19. Roadmap

| Version | Delivers | Rough shape |
|---|---|---|
| **V0** | Universe generation, molecule model, shape and binding, folding, decay, beaker harness, folding-map property tests, and a static shape inspector (§14.5). No world, no live viewer. Proves interesting chemistry happens in a test tube — and lets us *see* that it does. **Exit criteria in §23.** | Foundation |
| **V1** | Focused region, multi-scale engine, determinism and keyframes, neutral shadow, Lab viewer, the Chronicle. **Target: protocells.** ← *this document* | The hard, novel part |
| **V2** | Template-copying polymers, mutation, true heredity. Darwinian evolution proper. | Life gets a genome |
| **V3** | Cells, ecology, and the expansion to the 1,000,000 km² planetary field with real biomes. | The world gets big |
| **V4** | Gene regulatory networks, morphogenesis, body plans, sensory organs. **Eyes.** | The original dream |
| **V5** | Wonder skin, intervention tools, branching timeline comparison. Revisit **Bevy** here (§14.0) — the wgpu renderer carries forward. | Making it hers |

V1 is deliberately the hardest and least certain. Everything above it is comparatively well-trodden — evolution, ecology, and morphogenesis all have substantial prior art. Getting non-living generated chemistry to produce a self-replicating compartment is the part nobody has done this way, and it is the part worth doing first.

**A caution on V3 and V4.** Per §2.1, the assumption that lower-level self-maintaining units compose into higher-order ones is the *least* well-supported claim in this field. When V2 succeeds, generating new niches is likely to need explicit machinery rather than emerging for free — POET and Minimal Criterion Coevolution (Wang, Lehman, Clune & Stanley; Brant & Stanley) are the strongest existing answers and should be evaluated then, not assumed away now.

---

## 20. Risks

| Risk | Severity | Mitigation |
|---|---|---|
| **Chemistry never produces anything interesting** | Critical | The beaker battery (§7.2) rejects dead universes in seconds rather than hours. Tuning happens against the battery, not against overnight runs. |
| **Decay rate mistuned** | Critical | The most sensitive parameter in the system. Too fast and nothing survives long enough to fold; too slow and there is no selection at all, only accumulation (§9.4). The beaker battery's decay band (§7.2) is the primary control, and it is the first thing to suspect whenever a run produces nothing. |
| **Runs plateau silently** | High | Novelty tracked separately from activity (§15.2); stagnation declared only by model fitting, never by eye (§15.4). |
| **Emergence is drift, not selection** | High | Neutral shadow runs (§15.3) as a V1 deliverable. Without a null model, every result is an anecdote. |
| **Time compression insufficient for deep time** | High | Event-driven scheduling, quasi-equilibrium jumping, rare-event sampling. If still short, the governor trades spatial extent for temporal reach. |
| **Multi-scale tiers disagree, invalidating results** | High | Cross-tier agreement testing (§16) is non-negotiable and gates V1 completion. |
| **Region too small for complexity to accumulate** | Medium | Complexity carrying capacity scales with space-time volume (§2.7). If a run stalls, suspect region size *before* tuning parameters — this is a physical budget, not a knob. |
| **Folding is the bottleneck** | Medium | Aggressive caching with high expected hit rates; fold approximation for the bulk tier; a fold budget per timestep with graceful degradation. |
| **Determinism breaks under parallelism** | Medium | Counter-based RNG and fixed-order reduction designed in from the start; tested under varying thread counts. |
| **Scope creep toward V2 features** | Medium | V1 succeeds or fails on §15.1 alone. Genomes are explicitly deferred, however tempting they become. |
| **It works but is boring to watch** | Medium | The Chronicle is a V1 deliverable, not a polish item. An unnarrated eight-hour run is a failure regardless of what the chemistry did. |

---

## 21. Explicit non-goals

Not in V1, and in several cases not ever:

- **Real chemistry, in any form** — permanent, per §5
- Genomes, template copying, or Darwinian selection — V2
- Multicellularity, morphogenesis, or organs — V4
- The 1,000,000 km² planetary field — V3
- Intervention and God-mode tools — V5
- Multiplayer, cloud execution, or any hosted service
- GPU screening via wgpu (§13.5) — revisit after V1 succeeds. Cross-platform *CPU* support is in from V0 and is not deferred.
- Any claim, anywhere in the product or its documentation, that Borbax results say something about real biology

---

## 22. Resolved design questions

These were open at first draft and are now settled. The reasoning is preserved because several will want revisiting once there is evidence.

### 22.1 Region choice — single volcanic coastal shelf

The shelf already contains several distinct micro-environments with real gradients between them: a vent field, a tidal flat that wets and dries, porous mineral rock, and shallow photic water. Migration between conditions is therefore already a factor without paying for multiple regions, and the boundaries between them are where the most interesting chemistry should sit.

Multi-region ecology is a V3 concern, when there is something alive to migrate. Revisit before then only if runs show chemistry saturating in a way that looks environment-limited rather than parameter-limited — and note that per §2.7 the more likely cause of saturation is total space-time volume, not variety.

### 22.2 Signature resolution — 42 directions, confirmed empirically in V0

The sample directions are vertices of a subdivided icosahedron, so the available resolutions are not free parameters but a discrete ladder: **12** (the icosahedron itself), **42** (one subdivision), **162** (two). Only these preserve the symmetry that makes rotation an exact permutation (§8.2), so intermediate values are not on the table.

42 is the prior. 12 is likely too coarse to distinguish shapes that ought to bind differently; 162 is likely to cost roughly four times as much for detail the chemistry cannot use. But this is an empirical question the beaker harness answers in minutes, so it gets answered rather than assumed.

V0 sweeps all three against the folding-map property tests (§16) and the chemistry-diversity criteria (§7.2), and picks the smallest that passes all of them — reporting the cost per binding call alongside the chemistry results, so "smallest that passes" is read against what each actually costs. **Locking this value is a V0 exit criterion**: changing it later invalidates every golden hash and every tuned parameter downstream, so it must not be left drifting.

### 22.3 Neutral shadow — spawned on demand, not run continuously

A continuous shadow doubles compute for a control we only actually consult when making a claim. Since branching (§13.3) already provides fork-from-keyframe, a shadow can be spawned *retroactively* over exactly the interval a claim covers.

Cheaper, and strictly more flexible — a single run can be interrogated at several different points rather than carrying one fixed control. The accepted trade-off is that shadow intervals are quantised to the keyframe spacing, so a claim finer-grained than that spacing cannot be tested without re-running from the previous keyframe.

### 22.4 Failed-run policy — stop and report, for now

While the engine is under development, a run that produces nothing is the single most informative event available, and discarding it to try another seed would be throwing away the best diagnostic we have.

So a null run halts and writes a full report: the metrics trajectory, which battery criteria the universe passed and by how much, the last milestone reached, and the state at the point progress stopped. Automatic retry becomes appropriate only once null runs are understood well enough to be boring — a V2 concern at the earliest.

### 22.5 Chronicle voice — sample entries written before the Chronicle is built

The voice is harder than the plumbing, and it constrains the event schema: what an entry needs to *say* determines what each event has to carry.

So roughly ten sample entries get drafted across all three depth registers, covering the milestone list in §14.4 including the deaths, and read aloud to the intended audience during V1 planning. Whatever survives that reading defines the event schema. Doing this after the Chronicle exists would mean discovering the schema is wrong at the point it is most expensive to change.

### 22.6 Decay granularity — global rate first, per-bond refinement second

The decay band is the most sensitive parameter in the system (§20), and searching for it in one dimension is dramatically easier than searching the full space of the bond-energy matrix.

So V0 locates the band with a single global cleavage rate. Once found, per-bond rates derived from the generated matrix — which is what principle 1 actually wants, and which produces the far richer outcome where some linkages are durable and others fragile — are introduced *with the global rate as their mean*, so the band established in the first phase stays valid and the refinement is a redistribution rather than a fresh search.

### 22.7 Damage representation — implicit, with a measurability gate

A damaged molecule is simply a different molecule. No damage flag, no parallel bookkeeping, consistent with principle 2 — the shape changed, so the identity changed, and that is the whole story.

Damage load (§15.2) is therefore measured as divergence between a compartment's actual composition and the composition its RAF set implies it should have. The risk is that this measurement is too noisy to act on, which would make the aging account in §9.6 unobservable even while it is happening.

So it gets a **V0 exit criterion**: inject a known quantity of damage into a beaker and confirm the metric recovers it above noise. If it fails, the fallback is lightweight provenance — tag each molecule with the reaction that produced it, so damage is identified by origin rather than inferred from composition.

The fallback is more expensive than it looks, and the cost should be weighed before taking it. Per-molecule provenance breaks the interchangeability that §9.5 depends on, which turns decay from a single propensity channel back into a sweep over every molecule in the simulation — a change of asymptotic class, not a constant factor. If the measurability gate fails, reworking the *measurement* is worth attempting before accepting that.

### 22.8 The binding search covers rotations but not reflections

The icosahedral **rotation** group has 60 elements. The full symmetry group, which also admits reflections, has 120. Binding searches the former. This has to be settled before the resolution sweep in §22.2, because it doubles the cost of what is being swept.

This is also why §8.3's `ANTI` lookup exists and why omitting it is so damaging: without it the kernel searches the 60 *improper* elements instead — reflections and nothing else — which is the precise opposite of the decision recorded here, arrived at by accident.

**Reflections are excluded.** Two reasons, and the second is the interesting one.

The cheap reason is that it halves the most expensive kernel in the engine.

The real reason is that excluding reflections means **Borbax chemistry is handed**. A molecule and its mirror image are genuinely different species, binding different partners, and neither can substitute for the other. In three dimensions this is a substantial claim rather than a bookkeeping choice: a left-handed helix genuinely cannot be rotated onto a right-handed one, no matter how it is turned. Handedness becomes a real structural property that a lineage can commit to.

Which makes **homochirality** — a chemistry coming to use predominantly one handedness — an *emergent result the simulation could produce* rather than something ruled out by construction. That is one of the genuinely open questions about the origin of life, and having it available as a possible headline event is worth more than the modelling convenience of ignoring it.

There is a natural extension deliberately not taken now. A molecule tumbling freely in solution is not obviously prevented from presenting either face, whereas one adsorbed onto a **mineral surface** demonstrably cannot flip. Permitting reflections in free solution and forbidding them on surfaces would make mineral surfaces chirality-selecting environments — which is one of the actual proposed mechanisms for how homochirality arose, and would connect it directly to a ratchet §11 already provides. A V1+ opportunity, noted so it is not lost.

---

## 23. V0 exit criteria

V0 has no world and no viewer. It exists to establish that the chemistry is worth building a world around, and to lock the parameters everything downstream depends on. It is complete when:

1. `universe_gen` produces universes passing the full beaker battery (§7.2), including the decay band
2. The folding map demonstrably exhibits neutral networks and approximate shape-space covering (§2.5, §16)
3. Signature dimensionality is swept and **locked** (§22.2)
4. The decay band is located with a global rate (§22.6)
5. Damage load is shown to be measurable above noise, or the provenance fallback is adopted (§22.7)
6. Catalysis is observed emerging — a folded polymer with two cavities measurably accelerating a reaction, with nothing in the code that knows what an enzyme is (§8.5)
7. With decay disabled, a run and its shadow become statistically indistinguishable (§16) — the direct test of "no death, no life". Note this is a *different* experiment from the neutral shadow of §15.3, which equalises decay rather than removing it; both are needed and conflating them produces plausible output from the wrong control
8. Molecules, signatures, folds, cavities and binding pairs can be rendered and visually inspected, with the renderers under golden-image test (§14.5)
9. The metric families of §15.2 are computed — including the persistence-weighted Bedau–Packard statistics — and a shadow can be forked from a beaker run and compared against it (§15.3)

Criterion 9 was added late, and the reason is worth recording. Criterion 7 requires a shadow comparison, but nothing in the plan built one: §15's novelty metrics, the shadow, and the plateau fitting were all specified and then never scheduled. The plan could not reach its own bar. The alternative was to move criterion 7 to V1, which would have removed the direct experimental test of "no death, no life" — the claim §9.4 rests on — from the version that establishes the chemistry. Building it was the cheaper mistake to make: on a well-mixed beaker the statistics are straightforward, and the shadow reuses machinery §13.3 already needs.

Criterion 8 is listed last but should be built early — several of the criteria above it are far easier to evaluate, and far easier to *trust*, once the shapes can actually be seen.

Criterion 6 is the one that matters. If shape-derived catalysis does not emerge in a beaker, no amount of world-building above it will help, and the design needs rethinking before another line is written.
