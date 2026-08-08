# `borbax-universe` — inventing a periodic table from a seed

This crate takes a 64-bit number and produces a complete, self-consistent
chemistry: elements, how many bonds each can form, how strongly they bond, how
stable they are, and how common they are. Nothing is read from a data file and
nothing is copied from real chemistry (spec §5).

The trick is that almost everything is a **consequence of one geometric idea**,
rather than a separate number pulled out of a hat.

## The one idea: elements are packed spheres

Picture an element as a little cluster of identical units squeezed into a ball.
A cluster of 1 unit, 2 units, 3 units, … up to about 120. Everything else
follows from asking simple questions about that picture.

### How many units fit in each shell?

Real atoms have electron shells that hold 2, 8, 18, 32 … electrons. Borbax
invents its own rule instead, drawing one number `k` per universe:

```
shell_size(n) = k · n² + 2
```

so shell 1 holds `k + 2`, shell 2 holds `4k + 2`, shell 3 holds `9k + 2`, and so
on.

**Where the `n²` and the `+ 2` come from.** Not surface area — that heuristic
gets the leading term and nothing else, and it produces neither the `+ 2` nor the
fact that the linear term vanishes. The real source is **Euler's formula**: every
triangulated sphere satisfies `V − E + F = 2`. Carry that through the shell count
and the linear term cancels exactly, leaving `k·n² + 2` with no fudge. The same
identity forces `k = V − 2`, so `shell(1) = k + 2` is the number of units
touching the central one — which is why the coordination number is **not** a
second free parameter.

The `+ 2` is worth pinning down because the obvious reading is wrong and has
already been tried: it is *not* an antipodal pair. `packing.rs` records the
counterexample — a tetrahedron gives `k = 2`, shell 1 = 4, and there is no
antipodal pair anywhere in it.

Different universes draw different `k` and therefore get **genuinely different
periodic tables** — different period lengths, different points at which a shell
closes. That is spec §5's G3: the physics is deliberately *not* isomorphic to
ours.

**A coincidence you should know about, because it is live.** At `k = 10` the
shell law gives 12, 42, 92 — the shell populations of a real icosahedral packing
— and `k = 10` is **inside the drawn range** (6..14). `packing.rs` spends fifteen
lines on this: it is forced by Euler's formula with twelve five-fold vertices,
which is the same formula that produces this project's own sampling family
`D ∈ {12, 42, 162}`. **It is not corroboration and must never be read as any**,
and nothing may bias the draw toward `k = 10` on account of it.

(An earlier version of this README said the coincidence was about the number 2.
It is not, and 2 is outside the drawn range — so that sentence pointed away from
a hazard that can actually occur.)

### How many contacts does a cluster make?

Pack `N` units together and count how many touch each other: call it
`contacts(N)`. On average `contacts` grows *faster* than `N`, so contacts **per
unit** rises with cluster size — **but not monotonically, and the exception is
the interesting part.**

The unit that *opens* a new shell has almost nothing to nestle against, so
contacts-per-unit drops sharply just past every shell closure and climbs back as
the shell fills. Measured at `k = 6` it falls on **28 of 119 steps**, worst 7.3%
at `N = 9 → 10`.

That sawtooth is not noise. Its local maxima are exactly the closed shells, and
that is what gives the periodic families in §7.1 a *cause* rather than a
definition.

### How tightly bound is a cluster?

Two competing effects:

```
energy_per_unit(N)  =   ε · contacts(N) / N    ←  binding: touching is good
                      − σ · ∛(N²)              ←  strain: cramming is bad
```

- **The first term** rewards contacts. More neighbours, more binding.
- **The second term** punishes size. A lattice cannot tile a sphere, so each
  shell sits at a larger radius than the one below and its sites are stretched
  further apart than they want to be. Every site pays for that.

> **The exponent 2/3 is CHOSEN, not derived — and an earlier version of this
> README got that badly wrong.** It said strain "grows with the ball's surface".
> It does not: the term is a cost *per unit*, so **total** strain grows as
> `σ·N^(5/3)` — faster than the cluster's volume, let alone its surface. If
> total strain really grew as the surface (`N^(2/3)`), the per-unit term would be
> `σ·N^(−1/3)`, which *rises toward zero* and can never bend the binding curve
> down. Measured: under the shipped exponent the peak sits at `N = 55`; under the
> surface story it sits at the table edge and **there is no peak at all** — no
> instability gradient, no reason for heavy elements to decay.
>
> `element.rs` states the position directly and this README contradicted it:
> accumulating the strain shell-by-shell instead moves the peak in **48.1% of
> drawn universes**, and the file says in terms *"do not describe the exponent as
> derived."* It is a modelling choice with a recorded open item against it
> (Task 20), not a geometric necessity.

`ε` and `σ` are drawn per universe.

### The binding peak — Borbax's "iron"

Because binding rises and strain rises faster eventually, `energy_per_unit` goes
**up**, and in about **89%** of drawn universes peaks and comes back down. In the
other **~11%** the table is truncated before the peak, so the heaviest element is
the most tightly bound and there is no "past here it costs energy" regime at all.
Measured exhaustively over the drawn `(k, ε, σ, n_elements)` grid: 6,261 of
59,292 cells. Within a shell the curve saws — every closure is a local maximum.

The location of the peak is called `peak`, and it is an *output* — read off the
finished series, never an input.

Our universe has a peak of the same *shape* — in the iron-group region, and
strictly at nickel-62 rather than iron-56 once you measure binding per nucleon
rather than mass per nucleon. **That is where the analogy ends.** Nothing here is
calibrated against it, no constant is chosen to reproduce it, and the resemblance
is a consequence of "binding rises, strain rises faster" being a shape that
occurs in more than one place.

**`peak` must never be a parameter.** An earlier version passed it *in* and then
measured where the peak was — which is circular, and the code carries a long
comment about it. The test asserts it takes **at least** 6 distinct values across
24 drawn universes; the measured figure is **13**. (The 6 is the assertion's
floor, not a result — worth separating, because quoting a floor as a measurement
is how a figure goes stale with nobody noticing.)

### Everything else falls out

| Property | Where it comes from |
|---|---|
| `valence` — how many bonds an atom can form | Gaps in the partly-filled outer shell |
| `instability` — how far from stable | Relative distance *below* the binding peak: exactly 0 at the peak, rising toward both extremes **on average**, but sawing at every closure (non-monotone in 80% of grid cells below the peak, 43% above). Capped at 1.0 — and the cap fires on the monomer in **every** universe, a known open item routed to Task 15 |
| `abundance` — how common | Falls off with cluster size |
| `radius`, `mass` | Direct consequences of unit count and packing |

Note `instability` is the **nuclear-side** analogue — a dimensionless measure of
how far this element sits from the most tightly bound one. Deliberately *not* "a
probability that it falls apart": the field was renamed away from `decay_rate`
specifically because §7.1 and its consumer disagree about whether it is a
probability or an unbounded rate, and a likelihood gloss re-installs the reading
the rename exists to refuse. It is *not* the chemical decay of bonds breaking; see
the note in §7.1 of the spec, because the two were once confusingly both called
"decay".

## Bond energies

How strongly do two elements bond to each other? One line, and it is
deliberately not a drawn matrix:

```
E(a, b, order)  =  base · √( c(a) · c(b) ) · order^γ
```

- `c(a)` is element `a`'s **contact density** — contacts per unit, the same
  quantity from the packing model.
- `√(c(a) · c(b))` is the **geometric mean**. Two elements that both make many
  contacts bond strongly; if either is poor at contacts, the bond is weak.
  (Geometric rather than arithmetic mean because a bond needs *both* partners to
  contribute — an average would let a great partner rescue a hopeless one.)
- `order^γ` prices higher bond orders — up to **6**, not just double and triple;
  the attained valence ceiling is measured to be exactly `{4, 5, 6}`. `γ < 1` is
  drawn, which makes the series **diminishing**: each extra bond between two
  atoms is worth less than the one before.

**One consequence to know before building on this.** `base · √(c(a)·c(b))`
factorises into a per-element number, so the matrix is **exactly rank 1**. Every
element ranks its partners in the same order, and an `n`-element table stores
`n²` cells holding `n+1` independent numbers. There is no selective covalent bond
in this chemistry: partner preference comes from shape complementarity (§8.3) and
from nowhere else — deliberately, so it has one source rather than two.
`bonds.rs` flags this because the type is *called* a matrix, and a reader
building §9.1's energetics on top will otherwise assume pair structure exists.

Only **two dimensionless numbers** are drawn here (`base`'s multiplier and `γ`).
Everything with a unit attached comes from the element table, which is what
makes the bond matrix carry the table's energy *scale* rather than just its
ordering.

## Naming

Elements get invented names and symbols, checked against a blocklist of real
ones so no generated element can collide with a real element (spec §5, G2).
For every perturbed universe the blocklist exists purely to *exclude* — it is
never read as data. **G2 was revised 2026-08-07** (spec §5): the one
seed-equivalent configuration whose unperturbed constants reproduce the real
periodic table is expected to use the real names instead. That exception is
not implemented in this crate yet (see `naming.rs`'s module doc).

## Files

| File | What it holds |
|---|---|
| `packing.rs` | The shell law, contact counting, frontier geometry |
| `element.rs` | The element series and every derived property |
| `bonds.rs` | The pairwise bond-energy matrix and `BondOrder` |
| `naming.rs` | Name generation and the real-name blocklist |
| `lib.rs` | `Universe` — everything assembled, plus the drawn constants |

## How we know it still works: digests

Two tests hash every generated value into a single number and compare it against
a constant. If any physics changes, the number moves and the test fails loudly.

That is deliberate: **a moved digest is never a small thing.** It is either an
intentional physics change (say so in the commit message) or a bug. The tests
build the hash by *destructuring* each struct, so adding a field is a compile
error rather than a silently-unhashed value.
