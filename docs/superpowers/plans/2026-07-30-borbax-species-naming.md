# Design note — species get two names, one derived and one drawn

**Status: scaffolding. Not V0, not scheduled, not costed.** Recorded because the
gap is real and the determinism constraint on it is easy to get wrong in a way
no test would catch.

## The gap

Elements are named. `naming::mint` draws a symbol and a name from a syllable
grammar and checks both against a real-name blocklist (§5, G2). Measured over
2000 universes: 165 usable symbols, a 120-element table filled in 480 attempts
at the median, zero fallbacks.

**Species have nothing.** Grepped across the spec and both plan files: no
molecule naming, no formula rendering, no `Display` for a species anywhere. A
species is a `SpeciesId(u32)` and a `CanonForm`. For a simulation whose entire
output is *what did this universe make*, the things it makes are unnameable.

This is not a rendering detail. The Chronicle (V1), the SVG output, the CLI's
`beaker` and `battery` subcommands and any future report all need a way to refer
to a species that a person can hold in their head, and `SpeciesId(4471)` is not
it.

## The real-world model, and the half of it that is forbidden

Real chemistry carries both kinds of name, and they do different jobs:

| | example | derived from structure? | job |
|---|---|---|---|
| systematic | `propane`, `cyclohexene` | **yes** — `prop-` = 3 carbons, `-ane` = all single | says what a thing *is* |
| common | `water`, `ammonia` | no — historical accident | says which thing you *mean* |

Both are wanted here, for the same two reasons.

**The morphemes are forbidden.** `meth-`, `eth-`, `prop-`, `-ane`, `-ene`,
`-yne`, `-ol` are real chemical nomenclature. Importing them would imply a
mapping onto real chemistry as surely as an element called `He` would — the
exact impression §5 exists to prevent. Borbax needs invented affixes doing the
same job.

**The structure is forbidden too, and this is the deeper point.** IUPAC
nomenclature is built on a *privileged element*: find the longest carbon chain,
name it, treat everything else as a substituent. That works because carbon is
special. Borbax has 60–120 elements drawn from packing geometry and **no
privileged one** — and inventing one would be precisely the hardcoded special
case §3 forbids. A "backbone element" chosen by the naming layer would be a
chemistry claim smuggled in through a display feature.

So the systematic name must be built from invariants that hold for any molecule
with no distinguished atom: composition, cyclomatic number, branching, and the
bond orders present.

## Why this depends on Task 6

A systematic name must be a function of the *species*, not of the labelling.
Two relabellings of one molecule that produced different names would be the same
defect as two relabellings interning as different species — and it is the same
function that prevents both. `canonicalise` is the prerequisite; the name is
computed from `CanonForm` and inherits its isomorphism-invariance for free.

## The systematic name

A pure function of `CanonForm`. No draw, no state, no per-universe morphology —
so the nomenclature is **learnable across universes** while the element symbols
inside it stay universe-specific.

Three parts, each a graph invariant:

1. **Composition.** For each distinct `ElementId` present, a count stem and the
   element's symbol. Ordered by `ElementId` ascending — which, because
   `generate_elements` builds `for units in 1..=n_elements`, is *lightest
   first*. That ordering is derived rather than chosen, unlike Hill notation's
   carbon-first convention, which is another privileged-element assumption
   Borbax cannot make.

2. **Topology.** From the cyclomatic number `c = edges − vertices + 1` and the
   maximum degree:
   - `c = 0`, max degree ≤ 2 → unbranched chain
   - `c = 0`, max degree ≥ 3 → branched
   - `c = 1` → one ring
   - `c ≥ 2` → polycyclic

3. **Bond character.** The highest bond order present.

Illustrative morphemes only — the *constraint* is what matters, not these
particular syllables:

| count | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| stem | ka | du | vor | tesh | pel | zan | hesh | nol | rix | dov | keth | zur |

| topology | chain | branched | one ring | polycyclic |
|---|---|---|---|---|
| suffix | -ath | -esk | -orn | -ux |

| max order | 1 | 2 | 3 | ≥4 |
|---|---|---|---|---|
| suffix | -im | -av | -uv | -ozh |

So a molecule of two `Vo` and three `Kel`, unbranched, with one double bond,
reads `du-Vo-vor-Kel-ath-av`. Decodable in both directions, like `propane`, and
mapping onto no real nomenclature.

**The hyphens are a fiction guarantee, not typography.** A systematic name is a
*composed* string, so its parts can concatenate into a real chemical or
biological word by accident — and unlike a drawn name there is nothing to redraw
if they do, because systematicity is the whole point. Enumerating the collisions
would be the enumeration-of-spellings failure this project has now recorded
three times. A mandatory separator makes the class structurally disjoint from
the blocklist instead: no entry in `REAL_ELEMENT_SYMBOLS` or `REAL_WORDS` — the
actual constants, an earlier draft of this line named one that does not exist —
contains a hyphen, so no composed name can equal one, whatever the morphemes are
later changed to. Checked mechanically by review, along with all 20 illustrative
morphemes and all 192 unhyphenated `stem+topology+bond` concatenations against
`naming::is_real`: zero collisions in either arm. State that as the property the check decides,
and `xtask` can assert it over the grammar rather than over the outputs.

## The common name

Drawn, memorable, undecodable — the `water` to the systematic `-ath-im`.

Reuse `naming::mint`'s grammar and its blocklist. Two things differ from element
naming, and both are load-bearing.

**Key it on `CanonForm`, never on `SpeciesId`.** Ids are assigned in first-seen
order, so an id-keyed name would differ between a run and its neutral shadow,
and between two world seeds in the same universe — the same molecule called two
things, which destroys the one property a name is for. Keyed on the canonical
form, a shape carries the same name in every run of that universe. This is the
interner's own argument applied one level up.

**Give it its own stream coordinate.** Element naming holds
`Stream::new(seed, Domain::Naming, 0)`. Species naming must not perturb it —
Task 5 established the pattern when it put the bond constants on their own
`Domain::Universe` index precisely so adding one could not shift `temp_min` and
with it the decay band. Concretely, `Stream::new(seed, Domain::Naming, digest)`
where `digest` is a 64-bit hash of the `CanonForm`: order-independent by
construction, and the probability of landing on index 0 is 2⁻⁶⁴.

Note that `Domain::Naming`'s doc comment already says "element **and compound**
names", so the domain was reserved for this. The index axis inside it was not,
and reserving it is free today for the same reason `packed_index` was free in
Task 3 — nothing else uses a non-zero index, so no golden moves.

## Nothing may branch on a name

The `catalytic_class` → `outer_fill_band` rename is the precedent, and its
lesson was that a doc comment saying "reserved" does not reserve anything: a
`pub` field with a chemical-sounding name on a `pub` struct is *available*.

A species name is a label. It carries no chemistry, and no rate, propensity,
reaction channel or selection rule may read it.

**A first version of this section got the enforcement wrong, and the correction
is the useful part.** It said: keep names out of `SpeciesRecord` and in a side
table owned by the display layer, so the rate path cannot see one. That closes
the `catalytic_class` route — a *stored field* — and leaves a wider one open.

Both names are **pure functions of `CanonForm`**, and the placement below put the
systematic one "beside `canonicalise`" (in `borbax-molecule`) and the common one
"at intern" (in `borbax-reaction`). Task 14 puts `Interner` and `rate()` in the
*same crate*, and `borbax-reaction` depends on `borbax-molecule` by the fixed
crate order — so `rate.rs` could call `systematic_name(&form)` directly and never
touch the side table at all. A callable pure function is strictly more available
than a stored field: no struct to change, no allocation, reachable from anywhere
in or below `borbax-reaction`. The remedy contradicted its own discriminator.

**The property that actually holds it: the name-producing functions must live in
a crate that `borbax-reaction` does not depend on.** Given the fixed order
(`units → rng → universe → molecule → reaction → beaker`, with `geometry →
render` and `cli` above), that means at or above `borbax-render`, taking
`CanonForm` as input. Then the enforcement is `cargo tree` rather than a grep of
the rate path — mechanical, which is what this note asked for and did not
deliver.

*Discriminator:* `cargo tree -p borbax-reaction` must not list the crate holding
the naming functions. A version that keeps `systematic_name` in
`borbax-molecule` and adds a doc comment saying "do not call this from a rate"
is the `catalytic_class` remedy that already failed once.

Keeping names out of `SpeciesRecord` remains right, and also keeps it free of two
`String`s per species — but that was never the argument and is not the guarantee.

## Scope

**Not V0.** V0 has no Chronicle and no report; `beaker` prints counts. Nothing
in §23's exit criteria requires a species to be nameable, which is why this is a
note and not a task.

**Its home when picked up is Task 14**, where a species first exists — the
systematic name beside `canonicalise`, the common name at intern, the side
table owned by `borbax-cli` and `borbax-render`.

**Mixtures are out of scope.** Real chemistry names some mixtures (brass,
brine) and they are not systematic. Borbax's beaker is a population of species,
not a named substance, and there is no evident reason for V1 to change that.

## What to measure first, when it is picked up

Three numbers, none of which should be written into the implementation until it
has been run — the standing rule after `FRONTIER_COEFF` orphaned five prose
sites:

0. **Species-against-element collisions, which the list below first omitted.**
   `naming::mint` takes `taken: &mut Vec<String>`, so a species-naming pass that
   starts a fresh list can mint a species with the same name as an element in
   that same universe. Not a G2 breach — both clear the blocklist — but it
   defeats the single property a common name has, which is naming one thing.
   The `taken` list must span both, or the two must be drawn from disjoint
   sub-grammars.

1. **Common-name collision rate.** Two species sharing a drawn name is
   confusing rather than corrupting (`SpeciesId` remains the identity), but the
   grammar has to be sized against the species count a real run reaches, and
   that count is currently unknown. Element naming's 165-symbol space is
   nowhere near enough; measure the reachable species count from a battery run
   first, then size the grammar to it and report the collision rate over
   universes.

2. **Systematic-name injectivity.** How often do two distinct species share a
   systematic name? It is lossy by construction — composition plus four
   topology classes plus four bond classes cannot separate all isomers, and
   real nomenclature has the same property before locants are added. The
   question is whether the collision rate is low enough to be useful, and if
   not, which invariant to add next. Candidates in order of cost: degree
   sequence, then a canonical locant set from `CanonForm`'s own vertex order.

3. **That the hyphen property actually holds.** Assert over the blocklist that
   no entry contains a separator, as a test rather than as the sentence above.
   The guard that fails loudly on a target it cannot find is the one that saved
   `extract_const_value`.
