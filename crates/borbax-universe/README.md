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
on. The `n²` is not arbitrary — the surface area of a sphere grows with the
square of its radius, so the `n`th shell out has room for roughly `n²` times
what the first one does.

Different universes draw different `k` and therefore get **genuinely different
periodic tables** — different period lengths, different points at which a shell
closes. That is spec §5's G3: the physics is deliberately *not* isomorphic to
ours.

A file called `packing.rs` warns in a comment that the real world's analogous
number happens to be 2. That is a coincidence, and the comment exists to stop
anyone reading it as corroboration.

### How many contacts does a cluster make?

Pack `N` units together and count how many touch each other: call it
`contacts(N)`. Adding a unit to a small cluster gives it few neighbours; adding
one to a large cluster lets it nestle against many. So `contacts` grows *faster*
than `N` — contacts **per unit** rises as the cluster gets bigger.

### How tightly bound is a cluster?

Two competing effects:

```
energy_per_unit(N)  =   ε · contacts(N) / N    ←  binding: touching is good
                      − σ · ∛(N²)              ←  strain: cramming is bad
```

- **The first term** rewards contacts. More neighbours, more binding.
- **The second term** punishes size. Squeezing `N` units into a ball costs
  effort that grows with the ball's *surface*, and surface grows as `N^(2/3)` —
  which is what `∛(N²)` means. (Double the units and the surface goes up by
  about 1.6×, not 2×.)

`ε` and `σ` are drawn per universe.

### The binding peak — Borbax's "iron"

Because binding rises and strain rises faster eventually, `energy_per_unit` goes
**up, peaks, and comes back down**. The location of that peak is called `peak`,
and it is an *output* — read off the finished series, never an input.

In our universe the equivalent peak sits at iron, and it is why stars stop
fusing there: past iron, combining nuclei costs energy rather than releasing it.
Borbax gets the same shape from its own invented arithmetic.

**`peak` must never be a parameter.** An earlier version passed it *in* and then
measured where the peak was — which is circular, and the code carries a long
comment about it. Currently it takes at least 6 distinct values across 24 drawn
universes, which is the test that it is genuinely an output.

### Everything else falls out

| Property | Where it comes from |
|---|---|
| `valence` — how many bonds an atom can form | Gaps in the partly-filled outer shell |
| `instability` — how far from stable | Distance *below* the binding peak, 0 at the peak, rising both ways |
| `abundance` — how common | Falls off with cluster size |
| `radius`, `mass` | Direct consequences of unit count and packing |

Note `instability` is the **nuclear-side** analogue — how likely this element is
to fall apart on its own. It is *not* the chemical decay of bonds breaking; see
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
- `order^γ` prices double and triple bonds. `γ < 1` is drawn, which makes the
  series **diminishing**: the second bond between two atoms is worth less than
  the first.

Only **two dimensionless numbers** are drawn here (`base`'s multiplier and `γ`).
Everything with a unit attached comes from the element table, which is what
makes the bond matrix carry the table's energy *scale* rather than just its
ordering.

## Naming

Elements get invented names and symbols, checked against a blocklist of real
ones so no generated element can collide with a real element (spec §5, G2). The
blocklist exists purely to *exclude* — it is never read as data.

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
