# `borbax-rng` — where every random number comes from

**In one sentence:** this crate turns "which draw do you want?" into a number,
instead of keeping a running dice cup that has to be shaken in exactly the right
order.

## The problem it solves

A normal random number generator holds a hidden state and shuffles it forward
every time you ask for a number. That works fine until you want two guarantees
Borbax needs absolutely:

1. **The same seed must give the same universe** — on your Mac, on a Linux CI
   machine, on Windows, forever.
2. **Different parts of the simulation must not interfere.** If generating
   element names consumes one extra random number, the bond energies drawn
   afterwards must not all shift.

A running-state generator fails (2) badly. Everything downstream shares one
sequence, so an innocent change in one place silently rewrites everything after
it. This has already bitten the project once — see the note about deleting an
unused draw in the Task 5b plan, which would have renamed **1971 of 2000**
universes.

## The idea: counter-based generation

Instead of *shuffling* a state, a counter-based generator **computes** the
answer:

```
draw = F(key, counter)
```

`F` is a fixed scrambling function. Ask for counter 5 and you get the fifth
draw — directly, without generating draws 0 through 4 first. Ask again a week
later and you get the same number.

That makes each stream **independent by construction**. Borbax gives each
purpose its own key (`Domain::Universe`, `Domain::Naming`, …) and its own index,
so adding a draw to one cannot disturb another.

## The scrambling function: Philox 4x64-10

`F` is **Philox**, from Salmon, Moraes, Dror & Shaw, *"Parallel random numbers:
as easy as 1, 2, 3"*, SC'11, [doi:10.1145/2063384.2063405](https://doi.org/10.1145/2063384.2063405).

Here is the whole idea in miniature. Take two numbers. Multiply one by a fixed
constant, producing a result too big to fit — keep both the part that fits (the
"low" half) and the part that overflowed (the "high" half). Combine those with
the other number using XOR. Swap. Repeat **ten times**.

```
   L, R  ──►  multiply L by M  ──►  hi, lo
              R XOR hi XOR key  ──►  new L
              lo               ──►  new R          … ten rounds
```

The `key` is the part that matters for Borbax — it is what makes
`Domain::Universe` and `Domain::Naming` independent sequences rather than
offsets into one. (It is also bumped by a fixed constant each round, so the ten
rounds do not all mix in the same value.)

Two properties matter, and only the second needed a paper:

- **It is reversible.** For any fixed key, different counters *must* give
  different outputs — not "almost always", but always. Nothing can accidentally
  collide.
- **It scrambles thoroughly.** After ten rounds, flipping a single input bit
  changes each output bit with probability about one half. That is why counters
  0, 1, 2, 3 produce results with no visible relationship.

The multipliers `M0`/`M1` are the published ones, found by search. The
key-schedule constants `W0`/`W1` are **not** search results — they are
"nothing-up-my-sleeve" numbers you can recompute in one line:
`W0 = ⌊2⁶⁴/φ⌋` and `W1 = ⌊(√3 − 1)·2⁶⁴⌋`, both verified to the bit. That
distinction is worth keeping, because "recomputable" is what would catch a
transposed digit and "chosen by search" is not.

None of them may be adjusted. And if one is, **you will know**: `philox.rs`
carries the published known-answer vectors, so replacing `W0` with a
round-looking constant fails 6 of the crate's 46 tests — the first being
`matches_the_published_known_answer_vectors`, by name. (An earlier version of
this README said substitutes leave "most tests passing", which is arithmetically
true and reads as *you would not notice*. The one test whose entire purpose is to
notice fails first.)

## Why only whole-number arithmetic

Everything here is shifts, XOR, and multiplication that wraps around on
overflow. Those operations are **exactly specified by Rust's own semantics** —
two's complement, defined wrapping behaviour, no undefined cases — so they give
bit-identical answers on every processor. (Not by IEEE-754, which governs the
floating-point side and is `borbax-units`' subject, not this crate's.)

Anything involving decimals would not. That is why `sin`, `exp` and friends are
banned from this crate entirely — see [`borbax-units`](../borbax-units/README.md).

## What "wrapping" means

When a multiplication produces a number too large for 64 bits, Rust would
normally treat that as an error. Here it is the *point* — we deliberately keep
only the low 64 bits, which is what `wrapping_mul` does. The discarded high
bits are not lost; Philox captures them separately and mixes them back in.

## Files

| File | What it holds |
|---|---|
| `philox.rs` | The scrambling function and its constants |
| `lib.rs` | `Stream` and `Domain` — the interface everything else uses |

## The rules that go with it

- Never `rand::random`, never a thread-local generator, never `RandomState`.
  All three vary between runs or between machines.
- `Stream` is deliberately **not** `Copy`. Copying one would silently fork the
  sequence and replay the same draws — and duplicate "random" molecules would
  look like a chemistry result rather than a bug.
- Adding a draw shifts everything after it *in that stream only*. That is a
  deliberate physics change and needs the goldens regenerating.
