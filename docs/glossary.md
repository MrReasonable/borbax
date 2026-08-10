# Glossary

Borbax's code and commit history use two kinds of vocabulary a newcomer won't
have both at once: **real chemistry and physics terms** (Slater, Z_eff,
lanthanide, and so on — used since issue #26 to approximate real atomic
structure) and **Borbax's own invented vocabulary** (`affinity`, shape
signature, `IDEAL_GAP` — the made-up chemistry mechanism the project is built
around). This page explains both, assuming no background in either.

If a term you hit isn't here, treat that as a bug in this file, not in you —
open an issue or ask.

## Real chemistry and physics

### Atom, element, atomic number

An **atom** is a nucleus (protons and neutrons) surrounded by **electrons**.
The number of protons is the atom's **atomic number**, conventionally written
`Z` — it's what makes an atom the **element** it is. Change the proton count
and you have a different element; add or remove electrons instead and it's
still the same element, just charged (an *ion*).

Borbax writes atomic number as `z` throughout the code (lowercase, since it
also uses `Z` for other things) and calls the count of protons+electrons
**`units`** — the two are the same number for a neutral atom, which is what
Borbax always generates.

### Electron shells, subshells, and the `1s`/`2p`/`4f` notation

Electrons don't orbit the nucleus in one big cloud — they sort into
**shells**, numbered `n = 1, 2, 3, …` outward from the nucleus, and each
shell splits further into **subshells**, labelled by a letter: `s`, `p`, `d`,
`f` (in that order, going outward *within* a shell). A subshell's electrons
are described as an *orbital* shape — `s` is roughly spherical, `p` has two
lobes, `d` and `f` are more complex and more "spread out" near the nucleus but
thinner right at the center — this shape is what "penetration" (below) is
about.

So `4f` means "shell 4, subshell f". Each subshell holds a fixed maximum
number of electrons: `s` holds 2, `p` holds 6, `d` holds 10, `f` holds 14 —
this project derives that count from the subshell's own quantum numbers
(`capacity_of`), never a hardcoded table.

The pair `(n, l)` that the code uses everywhere is this same idea in a more
compact form: `n` is the shell number, and `l` is the subshell letter turned
into a number (`s=0, p=1, d=2, f=3, …`).

### Aufbau principle and Madelung's rule (fill order)

As you build up an atom by adding electrons one at a time (Z = 1, then 2,
then 3, …), they fill subshells in a specific, mostly-predictable order —
this is called the **aufbau principle** (German for "building up"). The
concrete rule for *which* subshell fills next is **Madelung's rule**
(1936, independently found by Klechkovsky): subshells fill in order of
ascending `n + l`, and when two subshells tie on `n + l`, the one with the
smaller `n` goes first. This is why `4s` fills before `3d`, even though `3d`
is numerically a lower shell — `4s` has `n+l=4`, `3d` has `n+l=5`.

Borbax's `madelung_order()` computes exactly this rule and needs no
chemistry-specific input at all — it's pure arithmetic on `(n, l)` pairs.

### Screening (a.k.a. shielding) and effective nuclear charge (`Z_eff`)

A nucleus with `Z` protons pulls on every electron with a force proportional
to `Z` — but only if that electron were alone. In a real atom, the other
electrons sit *between* an outer electron and the nucleus and partially
cancel the pull, the same way standing behind a crowd blocks some of a
spotlight. This is **screening** (also called **shielding**): the inner
electrons *screen* an outer electron from the full nuclear charge.

The charge an electron actually *feels*, after screening, is called the
**effective nuclear charge**, written `Z_eff`. It's always somewhere between
1 (fully screened) and `Z` (not screened at all). A bigger `Z_eff` means the
electron is held more tightly and the atom is smaller in that direction;
Borbax's `radius` formula divides by `Z_eff` directly for exactly this
reason.

### Penetration

Not every electron screens equally well, and this is the physical reason
why. An `s` orbital's shape means its electron spends a small amount of time
very close to the nucleus, "inside" the other electrons — it **penetrates**
well, which makes it (a) hard for other electrons to screen, and (b) *itself*
a good screener of anything further out. A `d` or, especially, an `f` orbital
penetrates much more poorly — it stays further from the nucleus on average
and rarely gets inside its neighbours. Standard teaching states the ordering
plainly: penetration ability runs `s > p > d > f`, and — because it's the
same physical property either way — so does *how well an electron in that
subshell screens other electrons*.

This is the fact behind Borbax's `sigma_deep` constant: an `f` electron is a
worse screener of everything else specifically because it penetrates poorly,
not because of anything about the *candidate* electron it's screening.

### Slater's rules

In 1930, physicist John C. Slater published a simple, hand-calculable recipe
for estimating `Z_eff`: sort an atom's electrons into groups by shell and
subshell type, and give each group a fixed "how much does one electron in
this group screen a candidate in that group" number (0.35 for two electrons
in the exact same subshell, 0.85 for the shell just inside an `s`/`p`
candidate, 1.00 — full strength — for everything else). It's an
approximation, not exact physics, and it's specifically known to be weak for
`d`- and `f`-block elements, because it treats "everything else" as
screening at a flat 1.00 regardless of how poorly an `f` electron actually
penetrates (see above). Borbax's `orbital.rs` implements Slater's rule
faithfully, then adds `sigma_deep` on top specifically to correct that known
weak spot.

### Valence, groups, and periods

An atom's **valence** is (loosely) how many electrons it has available to
form bonds with other atoms — the property that decides most of an element's
chemistry. Elements are traditionally arranged in the **periodic table**,
where each row is a **period** (elements with the same outermost shell
number) and each column is a **group** (elements with similar outer-electron
structure, and therefore similar valence and similar chemistry).

### HOMO and LUMO

Standard chemistry uses **HOMO** (Highest Occupied Molecular Orbital) and
**LUMO** (Lowest Unoccupied Molecular Orbital) to describe a *molecule's*
frontier electrons — the ones most available to react. Borbax borrows the
same idea one level down, at the level of a single atom's subshells: an
atom's own "HOMO" is its most recently filled subshell, and its "LUMO" is
whatever subshell would fill *next*. The gap between the two — how much
energy it costs to promote an electron from one to the other — is what
Borbax's `gap` computes, and it's what decides whether an atom's valence can
increase past what its closed subshell alone would give it (this is why, in
this model, beryllium-like atoms can bond but noble-gas-like atoms can't).

### Lanthanides and the lanthanide contraction

The **lanthanides** are the block of 14 elements (real ones: Lanthanum
through Lutetium, abbreviated **La → Lu**) where the `4f` subshell is being
filled — they sit in their own row, usually drawn below the main periodic
table for space. Chemically they're famous for being unusually similar to
each other, because filling an *inner* subshell (`4f`, buried well below the
outermost `6s` electrons) barely changes an atom's outward-facing chemistry.

What it *does* change, slightly, is size — and in the wrong direction from
naive expectation. Each step across the series adds one proton (which should
pull everything in tighter) and one `4f` electron (which, if it screened
perfectly, would exactly cancel that extra pull, leaving atomic size
unchanged). Because `4f` electrons penetrate and therefore screen *poorly*
(see above), the cancellation isn't quite exact — each added proton wins by a
small margin, and the atoms genuinely shrink, one small step at a time,
across the whole 14-element series. That real, measured shrinkage is the
**lanthanide contraction**. It's the standard textbook illustration of
imperfect `f`-electron screening, and it's exactly the effect Borbax's
`sigma_deep` fix (issue #26, finding F3) was calibrated to reproduce, since
literal Slater screening (see above) predicts *zero* contraction here — a
known, documented failure of the plain rule.

### Shannon radii

**Shannon radii** are a standard, widely used table of measured ionic radii
(published by R. D. Shannon in 1976), used as the real-world reference
whenever a model's output needs checking against reality. Borbax's
`sigma_deep` was calibrated against the real La → Lu contraction using these
figures: 103.2 pm (picometres) for La³⁺ down to 86.1 pm for Lu³⁺, a 16.6%
shrink across the series — the target the model's own 14-element
lanthanide-analogue block was tuned to reproduce.

## Borbax's own vocabulary

Borbax invents a chemistry from a seed number (see the top-level README) and
gives that mechanism its own names. This section is the short version; the
PRD (`docs/superpowers/specs/2026-07-26-borbax-prd.md`) is the long one, with
every term's full reasoning under its own numbered section (§*n*), cited
throughout the codebase.

**`universe_seed`** — the one number that determines everything about a
generated universe's chemistry: its elements, valences, bond energies, and
folding rules. Two runs with the same seed produce bit-identical results,
anywhere.

**Shape signature** — the geometric fingerprint Borbax reduces a molecule to:
for each of 60 directions around a sphere, how far out the molecule's surface
reaches (its "shape") and how positive or negative that surface is (its
"character"). Binding is a comparison of two signatures, nothing more.

**`affinity(A, B)`** — the single function that decides whether two shapes
stick together: does molecule A's bump-and-hollow pattern fit molecule B's,
and do their surface charges attract. The README's own opening formula is
this function; every kind of "sticking together" in Borbax — two molecules
reacting, an enzyme recognising a substrate, a membrane holding together — is
this same call, just at a different scale.

**`IDEAL_GAP`** — the separation distance at which two touching shapes are
considered a perfect fit (a "bump" reaching exactly as far as its partner's
"hollow" is empty).

**The icosahedral rotation group (the "60 rotations")** — instead of trying
every possible orientation of one shape against another (infinitely many),
Borbax samples direction space at the vertices of a subdivided icosahedron,
which has exactly 60 rotational symmetries. Every geometric comparison in the
codebase searches this fixed, precomputed set of 60 poses, never a
continuous space.

**Antipodal mapping / `contact_perms`** — when two shapes touch, a bump
pointing outward from A meets a hollow pointing *inward* from B at the same
contact point — i.e. in the *opposite* direction. `contact_perms` is the
precomputed table that pairs each of the 60 rotations with its "look at the
opposite direction" partner, so the binding calculation always compares the
physically correct pair of directions.

**RAF closure (autocatalytic set)** — Borbax never declares that a
"protocell" or "living" thing exists; instead it *detects* a specific,
well-defined mathematical property called a **Reflexively Autocatalytic,
Food-generated (RAF) set**: a group of molecules where every reaction is
catalysed by something else already in the group, and everything traces back
to raw ingredients available from the environment. When that property shows
up in a simulated beaker, that's the closest thing to "life happened" the
project claims — and it always falls out of the geometry, never a
special-cased rule.

**Invented units — `Thermal`, `Quanta`, `Span`, `WorldYear`, `Mass`** —
Borbax refuses to use real-world units (Kelvin, Joules, metres, seconds) for
anything, on principle: this is an invented universe, and its temperature,
energy, distance and time are their own distinct types the Rust compiler
enforces — mixing a `Thermal` with a `Quanta` in an expression is a compile
error, not a runtime bug.

**`PhysicsVersion` (V1 vs V2)** — Borbax originally invented its chemistry
purely from the seed with no connection to real-world physics at all (`V1`).
Since issue #26, `V2` instead *approximates real atomic physics* (screened-
hydrogenic orbital filling, Slater-style screening, and so on — everything
in the "Real chemistry and physics" section above) and lets the seed perturb
the physical *constants* away from their real-world values, rather than
inventing the mechanism itself from scratch. The one seed whose constants are
left completely unperturbed reproduces the real periodic table exactly,
including real element names — every other seed produces a plausible,
fictional variation on real physics.

**Rung / perturbation strength** — how far a given universe's seed has
pushed its physical constants away from the real, unperturbed ("our
universe") values. Rung 0 is the identity configuration described above;
higher rungs perturb more.

**Canonical form / canonicalisation** — before two molecules can be compared,
they need to be described the same way regardless of which atom you happened
to start labelling from. *Canonicalisation* picks one unambiguous graph
labelling for every molecule, so "the same molecule described two different
ways" always produces identical output.

**FCC lattice, folding** — the geometric grid (face-centred cubic, the way
oranges naturally stack) that a molecular chain folds itself onto in 3D,
producing the buried cavities that this project's catalysis mechanism
depends on.

## Where to go from here

- `docs/superpowers/specs/2026-07-26-borbax-prd.md` — the design and the full
  reasoning behind every decision above, cited by its own numbered sections.
- `crates/borbax-universe/README.md` — the shell-packing model in more
  depth (currently describes the original `V1` mechanism; the `V2`
  real-physics work summarised above lives in `crates/borbax-universe/src/orbital.rs`'s
  own doc comments until it gets its own README section).
- The top-level [README](../README.md) — the project's premise and how the
  pieces fit together.
