# Design note — a second axis, and stability from composition

**Status: scaffolding. Not V0, not scheduled, not costed.** Recorded because the
argument is sound and re-deriving it later would cost more than reading it.

## What this fixes, measured

Borbax today says the **light end is the most unstable**. `instability` is
distance below the binding peak, rising toward both extremes, so the light
quarter's mean is **7x** the heavy end's, and 86.8% of fission flux lands in the
lightest quintile. A lone atom is maximally unstable, tracing to
`energy_per_unit(1) = -sigma` — a cluster of one charged the strain of packing
with nothing to offset it.

That is backwards for anything resembling real matter, where small well-formed
things are stable and it is *imbalance* and *excess size* that destabilise. And
it is not merely cosmetic: the light end is where the interesting chemistry is
meant to happen.

The obvious repair — make instability one-sided, above the peak only — was
measured and **destroys the model**: the cascade stalls at dimers at every decay
strength, because the light end is where reservoir recycling happens.

So the light-side term is load-bearing *and* physically wrong. That combination
means the model is missing a degree of freedom, not a coefficient.

## What Borbax already has, stated precisely

Total binding is `eps * contacts(N)`, total cost is `sigma * N^(5/3)`. Measured,
doubling N:

| term | growth | |
|---|---|---|
| binding | **x2.33** | slightly faster than linear |
| cost | **x3.18** | **outgrows the binding** |
| *(a surface-like term would be)* | x1.59 | slower — Borbax has none |

So the cost term is **not** a surface term. It is a cost that outgrows the
binding, which is structurally the role repulsion plays in a real nucleus, and it
is what creates the binding peak. Borbax has a binding term and a
growth-limiting term.

**What it has no representation of at all is composition.** Every element is one
number. There is no ratio to be wrong, which is why `instability` has to be a
distance along size — the only axis available.

## The proposal

Give an element **two counts** rather than one. Nothing needs naming after
anything real; call them what the code calls them. Total size is their sum.

Then three terms, all with **drawn** coefficients:

```
binding      eps * contacts(N)              — already exists
size cost    sigma * N^(5/3)                — already exists
imbalance    kappa * (a - b)^2 / N          — NEW: the missing axis
```

The third term is what buys everything:

- **A valley of stability.** For each total size there is an optimal ratio.
  On-ratio elements are stable *regardless of size*, which is the fix: light
  elements stop being maximally unstable.
- **Isotopes fall out.** Same first count, different second — same "element",
  different mass, different stability. That is a second axis for neutral drift
  and for §9.5's mutation source.
- **Decay becomes a mechanism rather than a proxy.** An off-valley element moves
  toward the valley, releasing the difference. That is a real energy source
  (§9.5) instead of `instability` standing in for one.
- **Heavy instability keeps its cause.** The size cost still outgrows binding, so
  large elements remain unstable even on-ratio.
- **The binding peak stays emergent** — now from three competing terms instead of
  two, and `peak` remains an output.

## The §5 line, named precisely, because it was crossed once already

This session recorded a **G3 breach**: a channel-split argument derived from real
stellar nucleosynthesis, naming a real element and a real capture reaction, in a
block the plan says becomes code comments. It was retracted and rewritten from
Borbax's own arithmetic. The same discipline applies here and the line is sharp:

**Safe** — a functional *form*, with coefficients **drawn per universe**. That is
exactly what `eps` and `sigma` already are. "A binding term, a cost that outgrows
it, and a penalty for compositional imbalance" is a generic cluster-stability
shape, not an import.

**Not safe** — real coefficient values; real element or isotope names; any
calibration against real abundances or half-lives; and naming the terms after
real physics in code or comments. A term called `coulomb_term` fails G3 on its
face for the same reason "a rule that helps membranes form" fails: it names its
own referent.

**The test the codebase already states**, from `PhysicsVersion`: state the
physical process a law adds in terms that name no outcome outside the model. "A
cost that grows with compositional imbalance" passes. Anything naming a real
particle does not.

*And the exploration Ian wants is the §5-safe one by construction*: drawing the
coefficients and seeing what falls out is the opposite of calibrating to ours.
The interesting result is the **shape of the outcome space**, not a match.

## Scope, stated so this is not picked up by accident

This is a **Layer 0 redesign**, not a task:

- `Element` gains an axis; `ElementId` may stop being a single index.
- The shell law, valence, radius, mass and abundance are all functions of what is
  currently one number.
- Both pinned digests move, and every golden with them.
- §7.2's `f_w` reads valences that would now vary within an "element".
- `borbax-molecule` upward treats species identity as settled; isotopes change
  what "same element" means for canonicalisation.

V0 is chemistry in a test tube and its exit criteria (§23) do not ask for any of
this. Sequence it after V0 closes, or as a deliberate V1 physics change with the
golden regeneration that implies.

## What to do first, when it is picked up

**One measurement decides whether the third term earns its place**, and it can be
run before any of the above:

Add `kappa * (a - b)^2 / N` to a *probe* copy of the energy model, hold the
existing draws, and measure whether a valley of stability appears — i.e. whether
for each total size the most-bound composition is interior rather than at an
extreme, and in what fraction of drawn universes. If the valley is degenerate
(always at a - b = 0, or always at the edge) the term is decoration and the axis
is not paying for itself.

*Discriminator:* the fraction of universes with an interior valley, reported
alongside its spread. A term that produces the same ratio in every universe is a
constant wearing a mechanism's clothes — the defect this project has recorded
three times.
