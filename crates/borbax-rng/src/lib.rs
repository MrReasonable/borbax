//! Counter-based deterministic randomness (spec §13.1).
//!
//! Every random draw in Borbax comes from here. Streams are *counter-based*
//! rather than stateful-shared: a stream's output is a pure function of
//! (seed, domain, index, counter), so two threads drawing from different
//! streams cannot influence each other and results never depend on
//! scheduling order. That property is what makes bit-for-bit reproducibility
//! survive parallelism.
//!
//! Never use `rand`, thread-local RNGs, or `HashMap`'s default hasher in
//! any code path that affects simulation results.
//!
//! **Why this is hand-rolled, stated accurately** — because the wrong reason
//! gets reused to reject a dependency that would have been correct. It is
//! *not* portability: `rand_chacha` and `rand_pcg` are pure integer
//! arithmetic and are bit-identical across our three targets. The two real
//! reasons are that `rand` guarantees value stability only across *patch*
//! releases and has already broken it in a minor bump (0.9.0's changelog
//! has a "Reproducibility-breaking optimisations" section listing Canon's
//! and Lemire's methods for `Uniform`, "breaks value stability"). Note the
//! break is *not* lockfile-silent — cargo will not cross a 0.x minor, so
//! `rand = "0.9"` stays on 0.9.x — which makes this a weaker analogy to the
//! `libm` hazard than it first appears; the policy is the reason, not the
//! drift. The second is that `rand_distr`'s normal calls `f64::exp` and
//! `f64::ln` — the exact paths `clippy.toml` denies, invisible to clippy
//! because they sit inside the dependency. What is left to outsource after
//! `next_f64`, `next_range` and `next_normal` are ours anyway is the
//! five-line `mix` finalizer.
//!
//! **Stream identity is a 256-bit counter, not a hash, and that is the whole
//! determinism argument.** Draws come from Philox 4x64-10 (see `philox`), a
//! bijection of its counter for any fixed key. Every coordinate gets its own
//! 64-bit field — `[domain, index, sub, block]` — so distinct coordinates are
//! **provably** distinct streams rather than improbably-colliding ones.
//!
//! That is a change of kind, not degree, and it is worth recording why. The
//! first version of this crate hashed `(seed, domain, index)` — 192 bits of
//! coordinate — into a single 64-bit key. Two triples whose keys collided
//! emitted identical sequences forever, the seed cancelled out of the
//! collision equation so a colliding pair collided in *every* universe, and a
//! complete search below index 2³⁰ found five such relations among the nine
//! domains, two of them pairing [`Domain::Decay`] with [`Domain::Shadow`].
//! Nothing crashed and no test saw it.
//!
//! **Five was a local fluctuation, not structure**, and the honest number
//! belongs here rather than only the alarming half of it. The expectation for
//! distinct relations is
//! `36 × C(2³⁰,2)/2⁶⁴ = 1.125`, so five is a 4.4× excess at p = 0.006; but
//! extending the search to 2³² gives 25 against 18 expected, p = 0.068, and
//! `Decay`/`Shadow` does not grow while two other pairs overtake it. The old
//! design's defect was compressing 192 bits into 64 — an argument that needs
//! no inflated anecdote. Salmon et al. state the rule this design now satisfies:
//! "(key, counter) tuples must never be improperly re-used, since an
//! inadvertent re-use of the same tuple will result in exactly the same random
//! number, with potentially dire consequences for simulation accuracy."
//!
//! The 128-bit key is half-used on purpose: word 0 is the universe seed, word
//! 1 is reserved for the world seed. **Two words, two seeds** — but the
//! reservation needs one qualifier it did not have, and four things had been
//! promised those two words at once. The budget, stated once:
//!
//! | Coordinate | Home | Status |
//! |---|---|---|
//! | universe seed | `key[0]` | settled |
//! | world seed | `key[1]`, **universe-scoped domains excepted** | see below |
//! | `config_hash` | not a coordinate — it is `@n` | settled |
//! | replicate | `index`, via [`Stream::packed_index`] | settled, see below |
//! | branch sub-seed (§13.3) | **unallocated** | V1's problem, named below |
//!
//! **The world seed cannot reach universe-scoped domains.** `universe_gen` is
//! generated once per universe seed and cached; two worlds under one universe
//! share a periodic table by construction. So keying [`Domain::Universe`],
//! [`Domain::Naming`], [`Domain::Molecule`], [`Domain::Fold`] and
//! [`Domain::Hash`] on `[universe, world]` would regenerate the chemistry for
//! every planet, which is the architecture inverted. `key[1]` is zero for
//! those domains. An unqualified "word 1 is the world seed" was the earlier
//! reservation and it is not safe to implement literally.
//!
//! **A replicate is not a key, and putting it in `key[1]` would be a defect
//! rather than a contested choice.** A plan note once said that knob "is a
//! seed and belongs in `key[1]`". It cannot: the key reaches *every* domain,
//! so replicate *n* would draw a different periodic table — and the whole
//! purpose of replicates is to vary the run while holding the chemistry fixed.
//! It goes in the `index`, packed, where `Domain::Universe` never sees it.
//!
//! **The branch sub-seed has no home yet, and this is the honest state.**
//! §13.3 forks a keyframe with a new sub-seed; the two available shapes are a
//! third `index` field (there is no spare counter word — all four are spoken
//! for) or displacing the world seed in `key[1]` for world-scoped domains
//! only. Both are V1 decisions that need the keyframe format, which does not
//! exist. What must not happen is a fifth consumer quietly assuming `key[1]`
//! is free.
//!
//! §13.1's third element,
//! `config_hash`, is deliberately *not* a stream coordinate at all.
//! Keying on it would make universe generation config-dependent, so changing a
//! keyframe interval would regenerate the periodic table and break §13.4's
//! promise that a shared `U-…/W-…` pair means the same planet. It belongs in
//! the *physics version* — §13.1's `config_hash` and §6's `@n` are the same
//! thing under two names, so it travels in the shared world address rather than
//! beside it, and a recipient with different physics has a visibly different
//! address. Runtime knobs (keyframe interval, output, thread count) may not
//! reach the arithmetic at all; a knob that would is not a knob.

#![no_std]

// `core` is all the library needs. What this does and does not buy, measured
// rather than asserted, because the first version of this comment overclaimed
// it:
//
//   DOES block path-based `std` — a planted `std::time::Instant::now()` and a
//   planted `HashMap` both give `E0433`. Those are real §13.1 hazards and this
//   stops them at the crate boundary rather than in review.
//
//   DOES NOT block `x.exp()`, `x.ln()`, `x.cos()`, `x.powf()` or
//   `x.mul_add()`. All five still compile here. They are *inherent* impls on
//   the `f64` primitive contributed by `std`, and `std` is still in the crate
//   graph because `borbax-units` links it; rustc collects primitive inherent
//   impls from every loaded crate regardless of this attribute. So the one
//   §13.1 hazard this crate is most able to trip — `ln`/`cos` in `normal_from`,
//   which route through `det_math` and so do not trip it today — is NOT
//   protected here. `clippy::disallowed_methods` remains the sole
//   authority, exactly as it was before.
//
// One inversion worth recording: `f64::sqrt` resolves here *because*
// `borbax-units` links `std`. `next_down`, `MIN_POSITIVE`, `abs` and `TAU` all
// come from `core`, but `sqrt` does not — so if `borbax-units` ever becomes
// `no_std` (plausible; `det_math` is `libm`-backed and `libm` is `no_std`),
// this crate stops compiling with `E0599`. That is the wrong direction for the
// crate order and should be fixed there, not worked around here.
#[cfg(test)]
extern crate std;

mod philox;

use borbax_units::det_math;
use philox::philox4x64_10;

/// Which subsystem a stream belongs to.
///
/// Adding a new random draw in one subsystem cannot shift the sequence any
/// other subsystem sees — which means a change to folding does not silently
/// perturb reaction outcomes. **That property is exact**: a stream is a pure
/// function of `(seed, domain, index)`, so no number of draws in one domain
/// can move another.
///
/// Separation is **exact**, not probabilistic. The domain is its own 64-bit
/// field in Philox's counter rather than one input to a hash, so two domains
/// cannot share a stream at any index. An earlier design hashed all three
/// coordinates into a 64-bit key and admitted permanent collisions; see the
/// module documentation for what that cost.
///
/// **Discriminants are explicit and permanent, and a new variant may only be
/// appended.** This is the seam that lets physics be extended later without
/// invalidating every world already generated, and it costs nothing to
/// establish now. [`Stream::new`] derives its key from `(seed, domain,
/// index)`, so a V3 law whose constants come from a new `Domain = 10` leaves
/// every constant an existing universe already generated bit-identical — the
/// periodic table, the bond matrix, the binding weights, all unchanged.
///
/// Without that discipline the alternative is grim: adding a single
/// `next_f64` call to universe generation would shift every subsequent draw
/// in that stream, so every existing seed would silently generate a
/// *different* periodic table. Worlds shared as `U-7F3A21C9` would stop
/// meaning what they meant (§6), and nothing would announce it.
///
/// Two rules follow, and both are cheap:
///
/// 1. **Append; never renumber, never reuse.** A removed subsystem's
///    discriminant is retired, not recycled. Reordering this enum is a
///    physics change.
/// 2. **New physics draws from a new `Domain`.** Reaching into an existing
///    one to add "just one more constant" is the mistake this enum exists to
///    prevent.
// `Ord` is not decoration. Without it `BTreeMap<Domain, _>` does not compile
// while `HashMap<Domain, _>` does, so the public API would be steering every
// downstream crate toward the container §13.1 bans in result-affecting paths.
// The order is the discriminant order the whole append-only argument rests on.
// Note *why* that is free: derived `Ord` on a fieldless enum compares the
// discriminant VALUE, not the declaration index, so a derive and the
// discriminant order cannot disagree — an earlier version of this comment
// claimed a reordered enum would break them apart, which is impossible. What
// `ordering_agrees_with_discriminant` actually guards is a hand-written `impl
// Ord`, which it does catch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u64)]
pub enum Domain {
    /// Universe generation — elements, valences, bond energies, folding rules.
    Universe = 1,
    /// Generated element and compound names (§5, G2).
    Naming = 2,
    /// Per-species molecular work: graph construction and embedding.
    Molecule = 3,
    /// Chain folding on the FCC lattice.
    Fold = 4,
    /// Reaction selection and outcome.
    Reaction = 5,
    /// Decay *outcomes* — which bond orbit of the species, which neighbouring
    /// *species* drawn from abundance (§9.1 cleave, §9.5 radiogenic; §9.4 is
    /// why decay is mandatory, not how it resolves).
    ///
    /// Not decay scheduling. §9.5 is emphatic that decay is not separate
    /// machinery: thermal cleavage is "one more channel in the same scheduler
    /// that handles every other reaction", so channel *selection* belongs to
    /// [`Domain::Beaker`] by definition. A domain named for a decay scheduler
    /// is how a per-molecule decay sweep gets built by someone reading this
    /// enum for permission.
    Decay = 6,
    /// The beaker's own scheduling — Gillespie draws and mixing.
    Beaker = 7,
    /// The neutral shadow run (§2.7, §15.3).
    ///
    /// **The shadow draws from its own streams rather than sharing the live
    /// run's. Decided — and not for either reason this comment gave before.**
    ///
    /// Common random numbers are the alternative, and the prize is not small:
    /// a common-reaction-path coupling measures a **variance ratio of 0.0007
    /// at T = 400** on *cumulative* statistics, which §15.2 makes the
    /// load-bearing ones. Any argument that dismisses CRN as speculative is
    /// wrong, and this comment made it twice — once from a
    /// **noninterruption** failure, which is a property of GSMP clocks and
    /// structurally inapplicable to CRP; once by filing `run − shadow` as
    /// Glasserman & Yao category II when it is category **I**.
    ///
    /// **The obstruction is that the pairing CRN needs does not exist here.**
    /// A new species mints new reaction channels, and CRN requires draw *k* of
    /// channel *c* in one run to correspond to draw *k* of channel *c* in the
    /// other. Selection on versus off is *precisely* a difference in which
    /// species arise and when — so a channel present in both runs was born at
    /// different times after different numbers of draws, and a channel present
    /// in only one has no counterpart at all. There is no canonical Poisson
    /// index to share. That is an architectural obstruction rather than a
    /// variance claim, and unlike the two it replaces it needs no literature
    /// and cannot be refuted by a better-chosen coupling.
    ///
    /// Task 20b carries what would reopen it, and given the measured 1400×
    /// it should be looked at: a channel correspondence stable across both
    /// runs — a shared species-ID space minted from a shared stream is the
    /// obvious candidate — after which the remaining cost is the partial-CRN
    /// hazard of missing a draw on the shadow path.
    ///
    /// The variant's *original* justification — "must not consume the live
    /// run's draws" — was vacuous and is retracted. A `Stream` is a value with
    /// its own counter, so a shadow drawing from its own stream cannot advance
    /// the live run's whatever domain it uses; every domain gets that from the
    /// counter-based construction.
    ///
    /// See [`Stream::packed_index`] for the replicate coordinate this run
    /// needs and does not yet use.
    Shadow = 8,
    /// Hashing only — cache keys and content digests. Separated from the live
    /// domains so that using a `Stream` as a hash function cannot couple to
    /// the subsystem it names: adding a draw in folding must not change fold
    /// cache keys.
    Hash = 9,
}

impl Domain {
    /// Every variant, in discriminant order.
    ///
    /// Appending a variant to the enum means appending it here too. The
    /// compiler does not force that — what it does force, via the exhaustive
    /// `match` in `every_domain_stream_is_pinned`, is that a new variant
    /// cannot be added without giving it a pinned stream.
    pub const ALL: [Self; 9] = [
        Self::Universe,
        Self::Naming,
        Self::Molecule,
        Self::Fold,
        Self::Reaction,
        Self::Decay,
        Self::Beaker,
        Self::Shadow,
        Self::Hash,
    ];

    /// This domain's permanent discriminant.
    ///
    /// The single place the `#[repr(u64)]` cast is written, so the one
    /// suppression it needs sits here rather than at each call site. A
    /// `From<Domain> for u64` impl would avoid the cast, but [`Stream::new`]
    /// is `const fn` and trait methods cannot be called in const context.
    #[must_use]
    #[expect(
        clippy::as_conversions,
        reason = "the enum is #[repr(u64)] with explicit small literal discriminants, so this \
                  widens a known value and cannot truncate or change sign"
    )]
    pub const fn discriminant(self) -> u64 {
        self as u64
    }
}

/// Box-Muller's transform, split out from the draws so its boundary can be
/// tested directly: `next_f64` returns exactly `0.0` with probability 2^-53,
/// so sampling will never reach the one input that matters.
///
/// **Total over the reachable domain, which is `u1 ∈ [0, 1)` and `u2` finite.**
/// `ln(0)` is `-inf`, so `u1` is clamped up before it reaches `ln`, and the
/// cosine is bounded by 1. Both arguments come from [`Stream::next_f64`], so
/// both preconditions hold at the only call site — and that is what makes
/// `Mass::from_f64_quantised`'s `Err` arm unreachable in Task 4's generator,
/// where `mass = base + next_normal() * spread` is the only realistic route to
/// a bad value.
///
/// It is **not** total for every `f64`, and an earlier version of this comment
/// claimed it was. `u1 > 1` gives `ln(u1) > 0`, hence `sqrt` of a negative,
/// hence NaN — measured NaN at `1.0 + EPSILON`, `1.5`, `2.0`, `f64::MAX` and
/// `inf`. A non-finite `u2` is NaN through `cos`. Only the `u1` *lower* bound
/// and NaN-`u1` cases are handled, and a NaN reaching a state hash is the one
/// hazard `det_math`'s own module doc says survives for `+ - * / sqrt`.
///
/// **The clamp target is `2⁻⁵³`, the smallest value `next_f64` can return —
/// not `f64::MIN_POSITIVE`, which was wrong by 292 orders of magnitude.** Both
/// make the function total; only this one makes it continuous over the
/// reachable domain. Clamping to `MIN_POSITIVE` put a single atom of
/// probability at **37.640σ** when every other reachable input tops out at
/// **8.5717σ** — a 29σ gap of empty support, and a far tail about 9× heavier
/// than a true normal's. Nothing would ever have drawn it (P = 2⁻⁵³), but
/// anyone reasoning "`|z| ≤ 8.6`, so `base/spread > 8.6` keeps mass positive"
/// would have been wrong by 4.4×.
///
/// **Clamping, not rejection.** Redrawing until `u1 > 0` would also make this
/// total, but its draw count would depend on the value drawn, so a caller
/// could not know a stream's position without simulating it — and stream
/// position is what §13.2 replay and §13.3 branching restore. Clamping keeps
/// the count *statically* two, which is what makes
/// `next_normal_consumes_exactly_two_draws` a meaningful test. (The weaker
/// argument, that rejection would shift subsequent draws, is true but proves
/// little: it would do so once in 2⁵³ draws, and `next_range` accepts exactly
/// that kind of value-dependent count as a matter of course.)
///
/// The clamp is written as a comparison rather than `f64::max`, which §13.1
/// bans for its non-deterministic `±0.0` behaviour. The comparison form also
/// folds a NaN `u1` to the floor instead of propagating it.
fn normal_from(u1: f64, u2: f64) -> f64 {
    // 2^-53 — exactly `next_f64`'s resolution, so the clamp lands on the edge
    // of the reachable support rather than 969 binades below it.
    const MIN_U1: f64 = 1.0 / 9_007_199_254_740_992.0;
    let u1 = if u1 > MIN_U1 { u1 } else { MIN_U1 };
    (-2.0 * det_math::ln(u1)).sqrt() * det_math::cos(core::f64::consts::TAU * u2)
}

/// The largest multiple of `n` at or below `u64::MAX`. Draws at or above it
/// are rejected, so the accepted region holds a whole number of residue
/// classes and `v % n` is unbiased.
///
/// Split out from [`Stream::next_range`] so a test can assert the divisibility
/// directly. It has to be: the bias from getting this wrong is on the order of
/// 2^-64, which no statistical test on this side of the heat death would ever
/// see, and a test that recomputed the same expression would pass for a wrong
/// one just as happily.
///
/// Computing from `u64::MAX` rather than 2^64 discards up to `n` extra values.
/// That is never a bias — only a marginally higher rejection rate — and it
/// avoids needing 65 bits.
const fn rejection_zone(n: u64) -> u64 {
    u64::MAX - (u64::MAX % n)
}

/// Whether a raw draw falls in the accepted region for `zone`.
///
/// One comparison, split out for the same reason as [`rejection_zone`] itself
/// and pinned by `the_accept_boundary_is_strict`: the boundary is the only
/// input that distinguishes it from `v <= zone`, and sampling reaches `zone`
/// with probability 2⁻⁶⁴. Getting it wrong puts `zone % n == 0` back into the
/// accepted set one time too many — a 6.5e-19 relative excess on residue 0 at
/// `n = 12`, which is unbiased for every practical purpose and wrong for
/// exactly the reason the rejection loop exists at all.
const fn accepts(v: u64, zone: u64) -> bool {
    v < zone
}

/// Place a uniform draw `u ∈ [0, 1)` into `[lo, hi)`.
///
/// Split out from [`Stream::next_f64_range`] for the same reason as
/// [`normal_from`] and [`rejection_zone`]: the interesting inputs are the two
/// ends of the 2⁵³ grid, which sampling reaches with probability 2⁻⁵³. A test
/// that recomputed this expression instead of calling it would pass for a
/// wrong one just as happily.
///
/// **`u` is NOT monotone in the result**, and an earlier version of this
/// comment claimed it was and used that to argue four endpoint checks covered
/// all 2⁵³ draws. Measured on this function: `(1000.0, 1000.0 + 1e-9)`
/// descends at 0.39% of consecutive grid steps, and `(180.0, 260.0)` — a range
/// the tests actually use — descends at 18.75% near the top of the grid. The
/// bound still holds everywhere it was checked; the *argument* for why a few
/// points sufficed did not. The endpoint test is a spot check and says so.
fn range_from(u: f64, lo: f64, hi: f64) -> f64 {
    // Written as `lo < hi` rather than a negated comparison so it reads as a
    // total ordering question it is not: a NaN bound compares false against
    // everything and lands in the `else`, which returns `lo`.
    if lo < hi {
        let v = lo * (1.0 - u) + hi * u;
        if v < hi { v } else { hi.next_down() }
    } else {
        lo
    }
}

/// A deterministic stream of draws, keyed by `(seed, domain, index)`.
///
/// **Deliberately not `Copy`** (spec §13.1). A silently duplicated stream
/// replays the same sequence twice, which is a determinism bug that produces
/// entirely plausible output. Requiring an explicit `.clone()` makes the
/// duplication visible at the call site, and dropping `Copy` surfaces every
/// existing accidental by-value duplication as a compile error — which is the
/// whole point of the change.
#[derive(Clone)]
pub struct Stream {
    /// `[universe_seed, reserved]`. The second word is deliberately unused and
    /// deliberately present: it is the world seed's home for world-scoped
    /// domains, and a 64-bit key would force the two seeds together — exactly
    /// the hashing that produced permanent stream collisions in the first
    /// design. It stays zero for universe-scoped domains, or one universe's
    /// chemistry would vary by planet. §13.1's third element, `config_hash`,
    /// is not a stream coordinate. The full key/counter budget, including the
    /// one coordinate that still has no home, is in the module doc.
    key: [u64; 2],
    /// `[domain, index, sub, block]` — one 64-bit field each, no packing.
    /// Because Philox is a bijection of the counter, distinct coordinates
    /// give distinct output *by construction*.
    counter: [u64; 4],
    /// The current block. Philox emits four `u64` at once; buffering them is
    /// what makes the per-draw cost a quarter of a block rather than a whole
    /// one.
    buf: [u64; 4],
    /// How much of `buf` is spent, `0..=4`. Starts at 4, so a fresh stream is
    /// pure field assignment and the first draw pays for the first block —
    /// `Stream::new` does no Philox work at all.
    spent: u32,
}

/// Redacts the key.
///
/// **The derived `Debug` printed `key[0]` — the universe seed — verbatim, and
/// that was a live breach of a barrier this crate documents as holding.**
/// Demonstrated by a reviewer: format a `&Stream`, parse the seed out of the
/// string, call [`Stream::sub`], and the drawn words match
/// `Stream::sub(SECRET, ..)` exactly. So a routine holding only a `&Stream`
/// could mint a per-molecule sub-stream — the §8.6 regression that
/// [`Stream::sub`]'s own doc says is blocked at every leaf frame, and that
/// Task 20b's keyframe note names as the reason a field-wise `Serialize` must
/// never be added.
///
/// `E0599` on `s.seed()` verifies only that no *accessor* exists. It says
/// nothing about `Debug`, which is why the barrier needs this impl rather than
/// a note.
///
/// The counter and position stay visible: they are what a diagnostic actually
/// wants, and neither can reconstruct the key.
impl core::fmt::Debug for Stream {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Stream")
            .field("key", &"<redacted>")
            .field("counter", &self.counter)
            .field("spent", &self.spent)
            .finish_non_exhaustive()
    }
}

impl Stream {
    /// Open the stream for `(seed, domain, index)`.
    ///
    /// The result is a pure function of those three, which is the whole
    /// determinism argument: no stream can be influenced by how many draws
    /// any other stream has taken, so drawing order across subsystems is not
    /// part of the contract and never has to be preserved.
    ///
    /// **Distinct *counters* are provably distinct streams**, and distinct
    /// `(domain, index)` pairs always give distinct counters — they occupy
    /// their own 64-bit fields of a 256-bit Philox counter, and Philox is a
    /// bijection of that counter.
    ///
    /// The guarantee is stated over counters rather than over the whole
    /// `(seed, domain, index)` triple because **the seed is the key, and
    /// bijectivity is of the counter for a *fixed* key**. By pigeonhole it
    /// could not be otherwise: 320 bits of input map to 256 bits of output.
    /// Across seeds the guarantee is a birthday bound over 256 bits — which is
    /// unassailable in practice (2²² seeds at a fixed counter give 2²²
    /// distinct blocks, and adjacent seed pairs share none of 14 400 words),
    /// but it is the "improbable" kind, and this file draws a sharp line
    /// between the two. [`Stream::sub`] documents the one exception of its
    /// own: `sub(.., u64::MAX)` aliases the root.
    ///
    /// This is still the property the first design did not have, and its
    /// absence was not theoretical: hashing the same coordinates into a
    /// 64-bit key produced permanent seed-independent collisions — five below
    /// index 2³⁰ against 1.125 expected, which the module doc records as a
    /// local fluctuation rather than structure.
    ///
    /// `Stream::new` does no Philox work — it is field assignment, and the
    /// first draw pays for the first block.
    #[inline]
    #[must_use]
    pub const fn new(seed: u64, domain: Domain, index: u64) -> Self {
        Self {
            key: [seed, 0],
            counter: [domain.discriminant(), index, 0, 0],
            buf: [0; 4],
            spent: 4,
        }
    }

    /// Open a sub-stream of `(seed, domain, index)` — one per fold attempt,
    /// one per species.
    ///
    /// **This is a constructor rather than a method on [`Stream`], and that is
    /// what makes nesting impossible instead of merely forbidden.** The
    /// earlier design was `fork(&self, index) -> Stream`, an operation on a
    /// stream returning a stream, so `base.fork(a).fork(b)` type-checked — and
    /// was measured bit-identical to `base.fork(b)` for every `a`, because the
    /// child *overwrote* its parent's sub-coordinate rather than extending it.
    /// That is a duplicate stream, precisely the hazard this design exists to
    /// remove, and it shipped with a test asserting it as the contract.
    ///
    /// No operation here *derives a different* stream from an existing one, so
    /// the second level cannot be written down. (`Clone` does take a `&Stream`
    /// and yield a `Stream` — an earlier version of this sentence said nothing
    /// does, which is false and checkable. It duplicates rather than extends,
    /// `!Copy` exists to make that duplication visible at the call site, and
    /// `xtask` enforces the real invariant on the AST: nothing in this crate
    /// takes a `Stream` and returns one — any receiver, any return position
    /// including `Option`, tuples, arrays and `impl Trait`, free functions,
    /// trait methods, `impl` blocks on other types, and items nested in any
    /// function body including a trait method's default body. Names are
    /// resolved through `type` aliases and renamed imports within the file.
    /// **What it cannot see**, stated because three earlier versions of this
    /// sentence claimed coverage the code did not have: anything produced by
    /// macro expansion, and any alias imported from another module. Both are
    /// *reported* as unanalysable rather than passed over in silence, so the
    /// gate still fails — but it fails saying "did not look", not "looked and
    /// it was clean".) The alternative —
    /// `Stream<const FORKED: bool>` — would have made every downstream
    /// function that merely *draws* generic over depth, since the ordinary
    /// pattern is to hand a sub-stream to a routine that anneals. Removing the
    /// operation costs one type parameter, one `debug_assert` and one test,
    /// and buys a stronger guarantee.
    ///
    /// **Not one per molecule at the bulk or stochastic tiers.** To call this
    /// per molecule you need a stable `sub` per molecule, and stable means
    /// stored — which is per-molecule state arriving disguised as "just an
    /// index", the shape §8.6's regression always takes. §9.5 makes all
    /// molecules of a species interchangeable precisely so none needs
    /// individual state, and a per-molecule index turns a keyframe from
    /// O(species) into O(molecules). V1's molecular tier (§12.1) is the one
    /// place a per-molecule stream is legitimate; V0 has no per-molecule
    /// representation at all.
    ///
    /// **One sub-level is enough for §13.1**, which asks for a stream per
    /// subsystem, per patch, per time index — `domain`, `index`, `sub`, with
    /// `block` left as the draw counter. A further level needs its own
    /// [`Domain`], or its two indices **packed injectively** before the call —
    /// a 32|32 shift, not a hash. [`Stream::packed_index`] is that shift.
    /// Measured, a 64-bit hash of two coordinates
    /// collides only about 4e-7 of the time at 4e6 pairs, so it is nothing
    /// like the original defect; but it silently moves the guarantee from
    /// provably distinct to improbably colliding, and this file's whole
    /// argument is that those are different.
    ///
    /// **A sub-index range is owned by exactly one call site.** This warning
    /// existed on `fork` and was deleted with it; the hazard was not. Two sites
    /// both drawing `Stream::sub(seed, Domain::Fold, 7, 0..k)` share a
    /// sequence — a duplicate stream arriving through the one door `!Copy`
    /// cannot guard, and the only such door left. Give each call site a
    /// disjoint sub-range, or its own [`Domain`].
    ///
    /// The refactor made this **more** reachable, not less, and that is the
    /// honest cost of it: `fork` at least required a `&Stream` in hand, so
    /// provenance was structural. A constructor lets anything holding the seed
    /// mint any sub-stream for any `(domain, index)`. Measured consolation:
    /// transposing `index` and `sub` is a bijection on the pair, so it fails
    /// loudly rather than aliasing two streams silently. Not *every* word
    /// moves — `(0, 0)` is a fixed point, so 3 of `GOLDEN_SUB_STREAMS`' 9 are
    /// unchanged and 6 move. The test still fails; the protection is real.
    ///
    /// Sub-coordinates are offset by one because the root occupies zero. That
    /// costs exactly one value: `sub(.., u64::MAX)` wraps onto the root, since
    /// 64 bits cannot hold 2⁶⁴ sub-streams *and* a distinct root. Pinned in
    /// `sub_of_u64_max_aliases_the_root`.
    #[inline]
    #[must_use]
    pub const fn sub(seed: u64, domain: Domain, index: u64, sub: u64) -> Self {
        Self {
            key: [seed, 0],
            counter: [domain.discriminant(), index, sub.wrapping_add(1), 0],
            buf: [0; 4],
            spent: 4,
        }
    }

    /// Combine two 32-bit coordinates into one `index`, injectively.
    ///
    /// This is the "32|32 shift, not a hash" [`Stream::sub`] asks for, written
    /// once so the two halves cannot drift apart between call sites. `major`
    /// takes the high word, `minor` the low; distinct pairs give distinct
    /// results by construction, with no collision probability to quote.
    ///
    /// **It exists now for a consumer that does not exist yet, and that is the
    /// point.** §15.3's neutral shadow needs a *replicate* coordinate: a
    /// shadow run drawn from `(seed, Domain::Shadow, patch, event)` is fully
    /// determined by its fork point, so N replicates of it are bit-identical.
    /// Measured consequence — zero variance across replicates, which makes
    /// every difference between run and shadow read as significant, on a
    /// metric whose whole job is to say which differences are not. The fix is
    /// an axis to vary, and `index = packed_index(patch, replicate)` is it.
    ///
    /// **`packed_index(0, 0) == 0`**, so adopting this convention in V0 — one
    /// beaker, one shadow, no patches — moves no golden and changes no
    /// number. That is the entire argument for reserving it today rather than
    /// when the second replicate appears: at that point every stream
    /// coordinate downstream of a patch has to move, and the goldens with
    /// them. The cost now is one function and one test.
    ///
    /// It is deliberately *not* wired into [`Stream::new`] or [`Stream::sub`]
    /// as a separate pair of parameters. Most domains index by a single thing
    /// — a species, an element, a fold attempt — and a two-word signature
    /// would invite callers to invent a meaningless second coordinate. Callers
    /// that genuinely have two say so at the call site.
    #[inline]
    #[must_use]
    #[expect(
        clippy::as_conversions,
        reason = "u32 -> u64 is a widening conversion and cannot lose a bit; that losslessness \
                  is exactly what makes the packing injective, and is asserted in \
                  packed_index_is_injective"
    )]
    pub const fn packed_index(major: u32, minor: u32) -> u64 {
        ((major as u64) << 32) | (minor as u64)
    }

    /// The next raw 64 bits.
    ///
    /// **Cost: no fixed figure, deliberately.** Three review rounds put three
    /// different per-draw numbers in this comment and all three were wrong, so
    /// what is recorded here is the bound and the variable rather than a
    /// fourth number.
    ///
    /// The hardware bound is the **multiply port**: 40 multiply µops per block
    /// over two pipes is 20 cycles/block ≈ **1.55 ns/draw** on an M1 Pro
    /// P-core. It is not a latency chain — successive blocks of one stream
    /// differ only in `counter[3]`, which no block's output feeds, so they
    /// overlap. Measured, making blocks independent is 17–20% faster, which
    /// only overlap can explain.
    ///
    /// The variable that dominated every disagreement is **whether this
    /// function inlines into the caller's loop**, worth about 2.5× (~1.5 ns
    /// inlined against ~3.8 out-of-line), and decided by how many call sites
    /// the *calling* crate has — LLVM declines to inline a 210-instruction
    /// body into many. So no number here can be right for callers that do not
    /// exist yet, and Task 18's Gillespie loop is exactly where it will
    /// matter. Measure in the hot loop that cares.
    ///
    /// Two rounds of wrong diagnosis are worth more than the numbers were: the
    /// first asserted figures from a harness whose noise floor was never
    /// measured, and the second explained them away with a *physical floor*
    /// that does not exist — which is the more expensive error, because
    /// "measured wrong" invites re-measurement while "physically impossible"
    /// closes the question.
    ///
    /// `const` is not for the sake of const contexts — a stream is drawn from
    /// at runtime. It is a cheap guardrail: a `const fn` cannot reach a clock,
    /// an allocator, or a thread-local, which is the exact list §13.1 bans
    /// from a result-affecting path.
    ///
    /// The block counter wraps rather than saturating. Exhausting it means
    /// 2⁶⁶ draws on one stream, which the simulation will not reach; a wrap
    /// repeats the sequence rather than panicking mid-run, which is the right
    /// failure for something that cannot happen.
    #[inline]
    #[expect(
        clippy::indexing_slicing,
        clippy::as_conversions,
        reason = "`spent` is invariantly 0..=4 — set to 4 at construction, reset to 0 on refill, \
                  and incremented only after the branch above has forced it below 4; the widening \
                  to usize is of a value that has just been bounded by 4"
    )]
    pub const fn next_u64(&mut self) -> u64 {
        if self.spent >= 4 {
            self.buf = philox4x64_10(self.counter, self.key);
            self.counter[3] = self.counter[3].wrapping_add(1);
            self.spent = 0;
        }
        let v = self.buf[self.spent as usize];
        self.spent += 1;
        v
    }

    /// Uniform on the grid `{ n / 2⁵³ : 0 ≤ n < 2⁵³ }`, using the top 53 bits
    /// — exactly the mantissa width of `f64`. Every value on that grid is
    /// reachable and none is favoured.
    ///
    /// **Not "every representable value in `[0, 1)`"**, which an earlier
    /// version of this sentence claimed and which is false by a factor of 512:
    /// `f64` has about 2⁶² representable values in the interval, this reaches
    /// 2⁵³ of them, and everything below 2⁻⁵³ — 94.8% of the interval by
    /// count — is unreachable. The distinction is not pedantry, because three
    /// downstream constants follow from the real resolution and none from the
    /// false one: `0.0` occurs with probability 2⁻⁵³, so `if next_f64() < p`
    /// fires at 1.11e-16 for *any* smaller `p`; **`-ln(1 - u)`** truncates at
    /// 36.7368, which caps a Gillespie waiting time; and Box-Muller truncates
    /// at 8.5717σ (see `normal_from`, whose clamp target is this same 2⁻⁵³).
    ///
    /// **`-ln(u)` does *not* truncate, and writing the Gillespie step that way
    /// is a live hazard.** `u == 0.0` is reachable, so `-ln(u)` is `+inf`: the
    /// clock jumps to infinity once in 2⁵³ draws and the run silently stops
    /// scheduling. The bound above is real only for `1 - u`, or for a clamped
    /// argument the way `normal_from` clamps `u1`. **Task 18** owns this — the
    /// well-mixed beaker, where the Gillespie clock lives — and it is written
    /// into that task rather than left here, because a doc comment in this
    /// crate is not read while implementing that one.
    ///
    /// Both steps are exact: the shifted integer is below 2⁵³ so it converts
    /// without rounding, and `SCALE` is a power of two so the multiply cannot
    /// round either. The result is `n / 2⁵³` for an integer `n`, exactly.
    #[expect(
        clippy::as_conversions,
        reason = "std offers no lossless u64 -> f64 conversion; the shift above is what makes \
                  this one lossless, and the exactness is asserted in \
                  floats_land_exactly_on_the_53_bit_grid"
    )]
    #[expect(
        clippy::cast_precision_loss,
        reason = "the shift leaves at most 53 significant bits, which is exactly f64's mantissa \
                  width, so no precision is lost for any input"
    )]
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        // The denominator is 2^53, so SCALE is 2^-53 — both exact in f64.
        // Spelled as a literal rather than `(1u64 << 53) as f64` so the scale
        // factor needs no cast of its own.
        const SCALE: f64 = 1.0 / 9_007_199_254_740_992.0;
        (self.next_u64() >> 11) as f64 * SCALE
    }

    /// Uniform in `[0, n)`, rejection-sampled so there is no modulo bias.
    /// Returns 0 when `n == 0`.
    ///
    /// The unbiasedness rests entirely on the private `rejection_zone`, which
    /// is where it is tested.
    ///
    /// Returning 0 for an empty range is a deliberate contract, not an
    /// oversight: the alternatives are a panic (`clippy::panic` is `deny`
    /// workspace-wide, and a panic inside a propensity loop kills an
    /// overnight run) or an `Option` that every call site would immediately
    /// unwrap. It is pinned in `range_of_zero_returns_zero`.
    #[inline]
    pub const fn next_range(&mut self, n: u64) -> u64 {
        if n == 0 {
            return 0;
        }
        let zone = rejection_zone(n);
        loop {
            let v = self.next_u64();
            if accepts(v, zone) {
                return v % n;
            }
        }
    }

    /// Uniform in `[lo, hi)` — half-open at both ends, for every finite
    /// `lo < hi` with `|lo| > 2⁻⁹⁷⁰`. Returns `lo` when `lo >= hi` or either
    /// is NaN.
    ///
    /// **The magnitude qualifier is not decoration.** Below about 2⁻⁹⁷⁰ the
    /// product `hi * u` underflows into the subnormals, where its rounding
    /// error is absolute (2⁻¹⁰⁷⁵) rather than relative and can exceed the
    /// margin `u · (hi − lo)`, so the result lands one ulp *below* `lo`.
    /// Reproducible at `lo = 4.450_147_717_014_403_75e-308`, `hi = lo.next_up()`,
    /// `u = 2⁻⁵³`. Nothing in Borbax reaches 1e-293 — the reason to write the
    /// bound down is that an unqualified guarantee becomes folklore, and the
    /// next caller asserts on it.
    ///
    /// **A non-finite `hi` is swallowed, deliberately and dangerously; a
    /// non-finite `lo` is not.** Measured, and the asymmetry is the part worth
    /// knowing: `hi = NaN` returns `lo`; `hi = +inf` returns `f64::MAX` for
    /// every draw; but **`lo = -inf` with finite `hi` returns `-inf` on every
    /// draw** — and *not* via the fallback, which is the mechanism an earlier
    /// version of this comment named. `-inf < hi` is true, so the guard
    /// **passes** and the convex combination produces it: `lo * (1.0 - u)` is
    /// `-inf` for every reachable `u`, since `1 - u > 0` always.
    /// (`lo = -inf, hi = +inf` is finite, at `f64::MAX` — the endpoint clamp
    /// catches the `NaN` that `-inf + inf` produces.) Tightening the guard to
    /// demand finite bounds does *not* fix this: the fallback still returns
    /// `lo`. Worse, it would **regress three cases that currently return
    /// finite in-range values** — `(0.0, +inf)` and `(1.0, +inf)` from
    /// `f64::MAX` to `lo`, and `(-inf, +inf)` from `f64::MAX` to `-inf` — and
    /// the suite has teeth here: applying it fails
    /// `non_finite_bounds_are_swallowed_not_propagated`. There is no honest
    /// value to return for a non-finite `lo`, so the behaviour is pinned
    /// rather than invented.
    /// **The consumer this once named does not exist.** Task 18's selector is
    /// `SegmentTree::sample(u)` with `u` from `next_f64`, and the propensity
    /// multiply happens inside `sample`; `next_f64_range` appears in the plans
    /// only with literal constant bounds. Measured, the real hazard also
    /// inverts: a single NaN rate makes `sample` pick the **last** channel
    /// forever, not channel 0. That belongs to Task 18 — `total()` must be
    /// asserted finite before `sample` — and is written there rather than
    /// here. Pinned in `non_finite_bounds_are_swallowed_not_propagated`
    /// so the behaviour is a decision rather than an accident; the caller is
    /// responsible for not building a non-finite bound.
    ///
    /// **The obvious form is wrong twice, and its own test could not see
    /// either.** `lo + u * (hi - lo)` rounds up to exactly `hi`, at a rate of
    /// about `2⁻⁵² · |lo| / (hi − lo)` — negligible for a wide range, but
    /// **measured at 50%** for `(1.0, 1.0 + EPSILON)` and affecting roughly
    /// 42% of all `(lo, hi)` pairs at some `u`. It also returns `inf` or NaN
    /// whenever `hi - lo` overflows: `next_f64_range(-1e308, 1e308)` was
    /// non-finite on **every** draw. The original test used `(-3.5, 2.25)`,
    /// which is provably immune on both counts — 0 failures in 2⁵³ — so it
    /// would have passed against an arbitrarily broken implementation. rand
    /// shipped the same overshoot and had to break value stability in a minor
    /// release to fix it.
    ///
    /// The convex combination cannot overflow, because a weighted average of
    /// two finite values is bounded by them. It is bit-identical to the affine
    /// form whenever `lo == 0.0`, which is the shape a future Gillespie
    /// channel selection over `[0, total_propensity)` would use.
    ///
    /// The final comparison is the whole half-open guarantee. For a range
    /// narrower than `ulp(lo)` there are only a handful of representable
    /// values in it, so *no* implementation can be uniform there; what this
    /// guarantees is that the result is in range, with the top grid point
    /// folded onto the last representable value below `hi`.
    ///
    /// The draw is taken before the guard so the draw count is exactly one for
    /// every argument, including the rejected ones. A draw count that depended
    /// on the arguments would make stream position a function of caller data.
    #[inline]
    pub fn next_f64_range(&mut self, lo: f64, hi: f64) -> f64 {
        range_from(self.next_f64(), lo, hi)
    }

    /// Standard normal via Box-Muller. Used for property noise in the
    /// periodic-table generator.
    ///
    /// `ln` and `cos` route through [`det_math`] — they are not
    /// correctly-rounded by IEEE-754 and genuinely differ between platform
    /// libm implementations (§13.1). This is not a "last digit" concern: the
    /// noise here sets element masses and affinities, which set signatures,
    /// which are compared by `total_cmp` to pick a winning rotation. One ULP
    /// flips that argmax, and the result is a different species — not a
    /// slightly different number.
    ///
    /// The second Box-Muller variate is deliberately discarded rather than
    /// cached. Caching it would halve the draws and move every golden hash.
    ///
    /// **Never returns a non-finite value** — see the private `normal_from`,
    /// which is where that guarantee is implemented and tested.
    ///
    /// **`#[inline]` here was a 12% regression against the five-line mixer and
    /// is a win against Philox** — +4.7% in dev (10/10 paired runs), and
    /// inside the noise floor in release, where thin LTO inlines it anyway.
    /// An earlier version claimed 3–5% "in both profiles"; only the dev figure
    /// is resolvable. The attribute is justified on that alone.
    ///
    /// The two draws are bound to locals rather than passed inline. Rust does
    /// specify left-to-right argument evaluation, so the inline form would be
    /// correct; but "this is deterministic because of an evaluation-order
    /// rule" is a claim a future reader has to go and check, and two `let`s
    /// cost nothing.
    #[inline]
    pub fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        normal_from(u1, u2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // `Vec` is not in the prelude under `#![no_std]`. The tests are the only
    // place in the crate that allocates, which is the point.
    use std::vec::Vec;

    #[test]
    fn same_inputs_give_same_sequence() {
        let mut a = Stream::new(42, Domain::Universe, 0);
        let mut b = Stream::new(42, Domain::Universe, 0);
        let sa: Vec<u64> = (0..64).map(|_| a.next_u64()).collect();
        let sb: Vec<u64> = (0..64).map(|_| b.next_u64()).collect();
        assert_eq!(sa, sb);
    }

    #[test]
    fn different_domains_are_independent() {
        let mut a = Stream::new(42, Domain::Universe, 0);
        let mut b = Stream::new(42, Domain::Fold, 0);
        let sa: Vec<u64> = (0..64).map(|_| a.next_u64()).collect();
        let sb: Vec<u64> = (0..64).map(|_| b.next_u64()).collect();
        assert_ne!(sa, sb);
    }

    #[test]
    fn different_indices_are_independent() {
        let mut a = Stream::new(42, Domain::Fold, 0);
        let mut b = Stream::new(42, Domain::Fold, 1);
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn floats_are_in_unit_interval() {
        let mut s = Stream::new(7, Domain::Reaction, 0);
        for _ in 0..10_000 {
            let v = s.next_f64();
            assert!((0.0..1.0).contains(&v), "out of range: {v}");
        }
    }

    #[test]
    fn range_is_bounded_and_covers() {
        let mut s = Stream::new(7, Domain::Reaction, 0);
        // A bitmask rather than an indexed array: no `as`, no panicking index,
        // and a failure prints which residues were missed rather than just
        // that one was.
        let mut seen = 0_u64;
        for _ in 0..10_000 {
            let v = s.next_range(6);
            assert!(v < 6, "out of range: {v}");
            seen |= 1 << v;
        }
        assert_eq!(seen, 0b11_1111, "some values in 0..6 were never produced");
    }

    /// `next_range(0)` has no value it could honestly return, and the
    /// alternatives — panic, or `Option` — are both worse for a generator
    /// called from inside a propensity loop. It returns 0, and that contract
    /// is pinned here so it cannot drift into a panic later.
    #[test]
    fn range_of_zero_returns_zero() {
        let mut s = Stream::new(7, Domain::Reaction, 0);
        assert_eq!(s.next_range(0), 0);
    }

    /// `next_range(1)` is the other degenerate case, and unlike zero it has a
    /// correct answer. Worth pinning because the rejection zone is computed
    /// from `u64::MAX` rather than `2^64`, so `n == 1` is the input most
    /// likely to expose an off-by-one in that arithmetic.
    #[test]
    fn range_of_one_is_always_zero() {
        let mut s = Stream::new(13, Domain::Reaction, 0);
        for _ in 0..1000 {
            assert_eq!(s.next_range(1), 0);
        }
    }

    /// The property that makes [`rejection_zone`] unbiased, asserted on the
    /// arithmetic because it cannot be asserted on the output: at u64 width
    /// the bias from a wrong zone is around 2^-64, so `range_is_bounded_and_covers`
    /// would pass against `zone = u64::MAX` and every other plausible error.
    ///
    /// The boundary values are the interesting ones. Above `u64::MAX / 2` the
    /// zone collapses to `n` itself, and `n == u64::MAX` is where an
    /// off-by-one would leave a zone of 0 and spin the loop forever.
    #[test]
    fn the_rejection_zone_is_a_whole_number_of_classes() {
        let boundaries = [
            u64::MAX / 3,
            u64::MAX / 2,
            u64::MAX / 2 + 1,
            u64::MAX - 2,
            u64::MAX - 1,
            u64::MAX,
        ];
        for n in (1..=2000_u64).chain(boundaries) {
            let zone = rejection_zone(n);
            assert_eq!(zone % n, 0, "zone for n = {n} is not a multiple of n");
            assert!(zone > 0, "zone for n = {n} rejects every draw");
        }
    }

    /// **`next_range` was the last public draw method pinned by nothing, and
    /// six mutations survived the whole suite** — including Lemire's method,
    /// which changes 91.7% of every integer draw and is *the* optimisation
    /// `rand` 0.9.0 shipped that broke value stability, cited by name in this
    /// crate's own module doc as the reason it is hand-rolled. Also surviving:
    /// `(v ^ 1) % n` (100% of draws), `n - 1 - (v % n)` (100%),
    /// `v.swap_bytes() % n` (75%), `(v >> 1) % n` (91.6%), and removing the
    /// rejection entirely.
    ///
    /// Every one stays uniform — all residues present, counts within 5σ over
    /// 200 000 draws — so no statistical test can see them.
    /// `the_rejection_zone_is_a_whole_number_of_classes` tests the *helper*'s
    /// arithmetic; nothing tested that `next_range` uses it, or how the
    /// accepted value is folded. The plans call this for FCC neighbour
    /// selection, element ids, valence, masses, radii and decay constants — a
    /// changed fold would be a different periodic table for every seed, with
    /// the suite green.
    #[test]
    fn ranged_draws_are_pinned() {
        let mut s = Stream::new(0x5EED, Domain::Reaction, 0);
        let got = [
            s.next_range(1),
            s.next_range(2),
            s.next_range(6),
            s.next_range(12),
            s.next_range(60),
            s.next_range(1000),
            s.next_range(u64::MAX / 3),
            s.next_range(u64::MAX),
        ];
        assert_eq!(got, GOLDEN_RANGE);
    }

    const GOLDEN_RANGE: [u64; 8] = [
        0,
        0,
        1,
        8,
        21,
        539,
        417_892_413_605_005_757,
        12_312_119_947_423_759_239,
    ];

    /// The property no golden can express: that the rejection loop is *live*.
    ///
    /// Deleting `if v < zone` changes no value for any small `n` — the zone
    /// covers all but ~2⁻⁶⁴ of the space — so a golden cannot distinguish
    /// rejection from no rejection. What distinguishes them is *draw count*:
    /// at `n` just above `u64::MAX / 2` the zone is `n` itself, so nearly half
    /// of all draws are rejected and the stream advances by more than one.
    ///
    /// **Asserted as "some call took at least three draws", not "some call
    /// took more than one".** The weaker form was the first draft, and it
    /// proves only that retrying happens — a *bounded* retry passes it.
    /// Measured survivors of the weak assertion: `if v < zone { v % n } else {
    /// self.next_u64() % n }` (one retry, then take whatever comes) and a
    /// two-attempt bound. Both reintroduce exactly the bias the loop exists to
    /// remove, and both are the natural shape someone reaches for when
    /// "unbounded loop in a simulation step" looks like a hazard worth
    /// capping.
    ///
    /// The number that makes it work: at `n = 2⁶³ + 1` the shipped code
    /// consumes 393 draws over 200 calls, and **51 of those 200 calls need
    /// three or more**. P(a given call needs ≥3) is 1/4, so P(never seeing one
    /// across 200 calls) is 0.75²⁰⁰ ≈ 1e-25.
    ///
    /// This test carries the loop alone, and that is not an oversight:
    /// `GOLDEN_RANGE` is *provably* blind to rejection, its eight rows having
    /// measured rejection probabilities between 5.4e-20 and 3.3e-17. Neither
    /// test is redundant with the other.
    #[test]
    fn next_range_actually_rejects() {
        // Zone is `n` itself here, so P(reject) is just under 1/2.
        let n = u64::MAX / 2 + 2;
        let mut probe = Stream::new(3, Domain::Reaction, 0);
        let mut plain = Stream::new(3, Domain::Reaction, 0);
        let mut most_rejections = 0_u32;
        for _ in 0..200 {
            let _ = probe.next_range(n);
            // Advance the reference by one and count how far probe has gone
            // past it by comparing the next value each would produce.
            let _ = plain.next_u64();
            let mut rejections = 0_u32;
            while plain.clone().next_u64() != probe.clone().next_u64() {
                let _ = plain.next_u64();
                rejections += 1;
                assert!(rejections < 10_000, "runaway resync");
            }
            most_rejections = most_rejections.max(rejections);
        }
        assert!(
            most_rejections >= 2,
            "no call to next_range took three or more draws in 200 tries at a ~50% rejection \
             rate — the retry is bounded rather than a loop, and modulo bias is only partly \
             guarded (saw at most {most_rejections} rejection(s) in one call)"
        );
    }

    /// The accept boundary is strict, pinned at the three inputs that decide it.
    ///
    /// `v <= zone` instead of `v < zone` is unreachable by sampling — `zone`
    /// comes up once in 2⁶⁴ draws — so no golden, no statistical test and no
    /// draw-count test can distinguish it. It is a real defect: `zone` is a
    /// multiple of `n`, so admitting it adds one extra hit to residue 0.
    #[test]
    fn the_accept_boundary_is_strict() {
        let zone = rejection_zone(12);
        assert!(accepts(zone - 1, zone), "rejected the last accepted value");
        assert!(!accepts(zone, zone), "accepted the zone boundary itself");
        assert!(!accepts(zone + 1, zone), "accepted a value above the zone");
    }

    /// The seed must not be recoverable from a `Stream` by any route.
    ///
    /// A derived `Debug` printed it verbatim, and a reviewer demonstrated the
    /// full attack: format, parse `key[0]`, call `Stream::sub`, get identical
    /// draws. The named probe at the time — `s.seed()` must not compile —
    /// passed while that route was open, which is why this asserts on the
    /// rendered output rather than on the type system.
    ///
    /// **Asserted as independence, not as absence of four substrings**, and
    /// the difference is five measured leaks. The first draft searched the
    /// rendering for the seed in decimal, lower hex, upper hex and `{:#x}`;
    /// `swap_bytes()`, `to_le_bytes()`, `rotate_left(1)`, `{:o}` and
    /// `& 0xFFFF_FFFF` all pass it while printing the key. `to_le_bytes()` is
    /// the one to picture: it renders `key: [52, 18, 254, 202, 239, 190, 173,
    /// 222]` — the seed's own bytes, in order, for a caller holding nothing
    /// but a `&Stream` to read off and feed back to [`Stream::sub`]. Two of
    /// the four probes were redundant anyway, since `{secret:#x}` contains
    /// `{secret:x}` as a substring.
    ///
    /// Independence covers every encoding at once: if the rendering is the
    /// same for three different seeds, it carries no information about which
    /// one produced it.
    #[test]
    fn debug_is_independent_of_the_seed() {
        // Same domain, same index, same draw position — so the *only* thing
        // that differs between these is the key.
        let reference = std::format!("{:?}", Stream::new(0, Domain::Fold, 7));
        for seed in [1, 0xDEAD_BEEF_CAFE_1234, u64::MAX] {
            let rendered = std::format!("{:?}", Stream::new(seed, Domain::Fold, 7));
            assert_eq!(
                rendered, reference,
                "Debug output varies with the seed, so it encodes it: seed {seed:#x} rendered \
                 {rendered}, seed 0 rendered {reference}"
            );
        }
        // The counter is deliberately still visible — it is what a diagnostic
        // wants, and it cannot reconstruct the key.
        assert!(reference.contains("counter"), "Debug lost its useful half");
    }

    /// The mean of many uniform draws should sit near 0.5. This is a smoke
    /// test for gross bias, not a serious statistical test — a mixer that
    /// passed this and nothing else would still be inadequate, but one that
    /// fails it is definitely broken.
    ///
    /// The `sum` here is a sequential fold over a range, so its accumulation
    /// order is fixed by the iterator and not by a scheduler. §13.1 bans the
    /// *parallel* reductions, which is a different construct.
    #[test]
    fn is_not_grossly_biased() {
        let mut s = Stream::new(99, Domain::Decay, 0);
        let n = 100_000;
        let mean: f64 = (0..n).map(|_| s.next_f64()).sum::<f64>() / f64::from(n);
        assert!((mean - 0.5).abs() < 0.01, "mean was {mean}");
    }

    /// **The first draft compared mismatched positions** — `f1.next_u64()`
    /// against `f3.next_u64()` *after* `f1` had already been advanced by the
    /// equality check, so it tested `f1`'s word 1 against `f3`'s word 0. That
    /// passes with the `sub` argument ignored entirely: measured, a mutation
    /// writing a constant into `counter[2]` (making every sub-stream the same
    /// stream) left it green, because `w1 != w0` holds with probability
    /// 1 - 2^-64 whatever the streams are. It was renamed from its `fork`
    /// ancestor rather than repaired. Compare position for position.
    #[test]
    fn sub_streams_are_deterministic_and_independent() {
        let mut f1 = Stream::sub(1, Domain::Fold, 0, 3);
        let mut f2 = Stream::sub(1, Domain::Fold, 0, 3);
        let mut f3 = Stream::sub(1, Domain::Fold, 0, 4);
        let a: Vec<u64> = (0..8).map(|_| f1.next_u64()).collect();
        let b: Vec<u64> = (0..8).map(|_| f2.next_u64()).collect();
        let c: Vec<u64> = (0..8).map(|_| f3.next_u64()).collect();
        assert_eq!(a, b, "the same sub-coordinate gave different streams");
        assert_ne!(a, c, "different sub-coordinates gave the same stream");
    }

    /// **`Stream::sub`'s encoding was pinned by nothing, and the gap was not
    /// theoretical.** Measured: replacing `sub.wrapping_add(1)` with `!sub`
    /// survives all 38 other tests, because `!` is a bijection (so injectivity
    /// holds) and `!u64::MAX == 0` (so the root alias holds) — the only two
    /// properties the `sub`-named tests check. It changes every sub-stream in
    /// every universe: 0 of the first 4 words survive. `Stream::new` carried
    /// two goldens and `Stream::sub` carried none, with nothing explaining the
    /// asymmetry.
    ///
    /// Fold attempts and per-species work draw from sub-streams, so a tidy-up
    /// of that one line would silently regenerate every fold trajectory ever
    /// produced — observable only as "the folds got worse".
    ///
    /// **Keep an odd `sub` in the list.** `s ^ 1 == s + 1` for every even `s`,
    /// so an all-even row set is *provably blind* to `sub ^ 1` as an encoding —
    /// computed, it produces byte-identical goldens. Row `(0, 3)` is the only
    /// one of the three that separates them, and it carries 3 of the 9 words.
    /// Same shape as the `(-3.5, 2.25)` note on `floats_in_a_range_are_bounded`
    /// and for the same reason: keep the number.
    ///
    /// A failure here means the same thing as any other moved golden.
    #[test]
    fn sub_streams_are_pinned() {
        // Collected rather than written through an iterator: the first draft
        // used `if let Some(slot) = w.next()`, so adding a fourth coordinate
        // silently discarded its draws and the test still passed — broadening
        // the golden's coverage would have broadened nothing.
        let got: Vec<u64> = [(0_u64, 0_u64), (0, 3), (7, 2)]
            .into_iter()
            .flat_map(|(index, sub)| {
                let mut s = Stream::sub(0x5EED, Domain::Fold, index, sub);
                core::iter::repeat_with(move || s.next_u64()).take(3)
            })
            .collect();
        assert_eq!(got, GOLDEN_SUB_STREAMS);
    }

    const GOLDEN_SUB_STREAMS: [u64; 9] = [
        14_164_468_725_724_755_671,
        12_873_852_386_610_874_040,
        7_495_883_468_757_775_053,
        9_018_786_775_305_268_180,
        3_890_255_439_183_115_723,
        11_443_634_209_256_710_128,
        9_659_346_300_237_333_348,
        3_948_479_596_175_888_705,
        2_411_997_774_858_210_089,
    ];

    /// A sub-stream writes `sub + 1` into the counter's third field, so
    /// distinct sub-coordinates are distinct counters and Philox's bijectivity
    /// does the rest. The old design needed a hand-chosen odd constant to stop
    /// `fork(0)`'s derived key colliding with the parent's first output; that
    /// whole class of concern is gone, because there is no key derivation left
    /// to collide.
    #[test]
    fn sub_streams_are_distinct_from_the_root_and_each_other() {
        let mut root = Stream::new(1, Domain::Fold, 0);
        let mut s0 = Stream::sub(1, Domain::Fold, 0, 0);
        let mut s1 = Stream::sub(1, Domain::Fold, 0, 1);
        let (p, a, b) = (root.next_u64(), s0.next_u64(), s1.next_u64());
        assert_ne!(a, p, "sub(0) collides with the root");
        assert_ne!(a, b, "sub(0) and sub(1) collide");
        // The name says "from the root and each other"; the first draft never
        // checked sub(1) against the root, so a sub-stream colliding with it at
        // any index but 0 would have passed.
        assert_ne!(b, p, "sub(1) collides with the root");
    }

    /// The counter-based property itself, asserted directly rather than
    /// through constants.
    ///
    /// **This is strictly stronger than lengthening a golden, and it is what a
    /// golden cannot do.** `GOLDEN_UNIVERSE_0` exercises the block-advance
    /// function only at the blocks it happens to span — at 8 words that was
    /// `g(0)` alone, which let a mutation giving every stream period 8 pass;
    /// at 12 it is `g(0)` and `g(1)`. A reviewer then found a survivor even at
    /// 12: `counter[3] += 1 + counter[2]`, identical to the real generator for
    /// the root and so invisible to every golden. Lengthening the constant
    /// only moves that boundary.
    ///
    /// The defining property of a counter-based generator is that block `k` is
    /// computable *directly*, without having produced blocks `0..k`. Asserting
    /// that pins the advance for every `k` at once, with no constants to
    /// regenerate and nothing to go stale.
    #[test]
    fn draws_are_the_direct_block_function_of_their_position() {
        for (seed, domain, index, sub) in [
            (0_u64, Domain::Universe, 0_u64, None),
            (0x5EED, Domain::Fold, 7, None),
            // **Keep this row.** It is the only one with an odd `sub`, and the
            // only thing that catches `counter[3] += 1 + counter[2]` — a block
            // advance that agrees with the real generator whenever `counter[2]`
            // is zero, which is every root stream.
            (0x5EED, Domain::Fold, 7, Some(3_u64)),
            // **Keep these two rows too, for the axis one below them.** Every
            // other row here, and every index anywhere else in this suite, sits
            // in `0..8`. Measured before they existed: `index as u32 as u64` in
            // `new`, the same in `sub`, `index & 0xFF`, `index & 0x7`, and
            // `counter[3] += 1 + (counter[1] >> 3)` — five mutations, all five
            // green against the whole suite. The sub axis was already probed at
            // `u64::MAX` by `sub_of_u64_max_aliases_the_root`; the index axis
            // had no equivalent, and `Stream::sub`'s own doc *recommends*
            // packing two indices into it, which produces exactly the large
            // values nothing else here reaches.
            (0x5EED, Domain::Molecule, (1_u64 << 32) + 9, None),
            (0x5EED, Domain::Molecule, u64::MAX, Some(2)),
        ] {
            let mut s = sub.map_or_else(
                || Stream::new(seed, domain, index),
                |v| Stream::sub(seed, domain, index, v),
            );
            let sub_field = sub.map_or(0, |v| v.wrapping_add(1));
            for k in 0..20_u64 {
                let want =
                    philox::philox4x64_10([domain.discriminant(), index, sub_field, k], [seed, 0]);
                for (j, expected) in want.iter().enumerate() {
                    assert_eq!(
                        s.next_u64(),
                        *expected,
                        "block {k} word {j} is not the direct block function's output"
                    );
                }
            }
        }
    }

    /// **Which coordinate lives in which counter field is load-bearing, and
    /// four mutations of it passed all 38 tests before this existed.**
    /// Measured, each silently catastrophic and each invisible:
    ///
    /// - the sub-coordinate written to field 3: a sub-stream *is* the root
    ///   advanced 4(j+1) draws — `sub 0` and `sub 1` shared 36 of their first
    ///   40 words.
    /// - `counter[1]` dropped: every index in a domain shares one sub-stream
    ///   family, 40 of 40 words identical.
    /// - `new` writing index into field 2: `new(s, Fold, 3)` *is*
    ///   `sub(s, Fold, 0, 2)` — the index and sub axes become one.
    /// - `new` writing index into field 3: every index stream is a bit-exact
    ///   *suffix* of index 0's.
    ///
    /// The old sub-stream tests all survived because they used `index = 0` and
    /// asserted a single-word `assert_ne!`. The consequence of the last one is
    /// the sharpest: `Domain::Molecule` indexed by species means species *k*
    /// draws numbers species 0 already used, signatures correlate along the
    /// index axis, and the periodic table's diversity quietly collapses —
    /// which reads as "this universe was a dud", with the decay band blamed
    /// first exactly as CLAUDE.md predicts.
    ///
    /// Sixteen words per stream is deliberately four blocks, so a
    /// block-coordinate shift cannot hide inside one. Accidental-collision
    /// probability is 1440²/2⁶⁵ ≈ 5.6e-14.
    #[test]
    fn distinct_coordinates_share_no_words() {
        let mut seen = std::collections::BTreeSet::new();
        let mut total = 0_usize;
        for domain in [Domain::Universe, Domain::Fold, Domain::Molecule] {
            for index in 0..6_u64 {
                for sub in 0..5_u64 {
                    let mut s = if sub == 0 {
                        Stream::new(0x5EED, domain, index)
                    } else {
                        Stream::sub(0x5EED, domain, index, sub - 1)
                    };
                    for _ in 0..16 {
                        seen.insert(s.next_u64());
                        total += 1;
                    }
                }
            }
        }
        assert_eq!(
            seen.len(),
            total,
            "{} of {total} words are shared between distinct (domain, index, sub) streams — \
             a counter coordinate has been transposed or dropped",
            total - seen.len()
        );
    }

    /// The index axis is 64 bits wide, and before this test nothing in the
    /// suite reached past 7.
    ///
    /// `distinct_coordinates_share_no_words` sweeps the domain and sub axes
    /// densely at small index; this sweeps index sparsely and far. Both are
    /// needed, and the asymmetry that motivated it is that the *sub* axis was
    /// already probed at `u64::MAX` while the index axis — the one
    /// [`Stream::sub`] tells callers to pack two coordinates into — was not.
    ///
    /// The consequence worth holding on to is the mask, not the truncation.
    /// Under `index & 0xFF`, `Stream::new(seed, Molecule, 256)` emits the
    /// identical sequence to index 0. [`Domain::Molecule`] is indexed by
    /// species, so a beaker reaching 300 species has species 256..300 drawing
    /// numbers species 0..44 already used: signatures correlate along the
    /// index axis and the periodic table's diversity collapses. That surfaces
    /// as "this universe was a dud", with the decay band blamed first, exactly
    /// as CLAUDE.md predicts — and `distinct_coordinates_share_no_words`
    /// names that same failure while sweeping `0..6` and being unable to see it.
    #[test]
    fn the_index_axis_is_sixty_four_bits_wide() {
        let mut seen = std::collections::BTreeSet::new();
        let mut total = 0_usize;
        for index in [
            0,
            1,
            7,
            8,
            255,
            256,
            1_u64 << 32,
            (1_u64 << 32) + 7,
            u64::MAX - 1,
            u64::MAX,
        ] {
            // Both constructors: the truncation was measured surviving in each
            // of them independently, so probing one would have caught half.
            for sub in [None, Some(0_u64)] {
                let mut s = sub.map_or_else(
                    || Stream::new(0x5EED, Domain::Molecule, index),
                    |v| Stream::sub(0x5EED, Domain::Molecule, index, v),
                );
                for _ in 0..8 {
                    seen.insert(s.next_u64());
                    total += 1;
                }
            }
        }
        assert_eq!(
            seen.len(),
            total,
            "{} of {total} words are shared between streams differing only in `index` — the \
             index axis is truncated or masked somewhere below 64 bits",
            total - seen.len()
        );
    }

    /// The packing is injective, and `(0, 0)` is 0.
    ///
    /// The second half is what makes adopting the convention in V0 free: a
    /// caller with one patch and one replicate passes `packed_index(0, 0)`,
    /// which is the `index = 0` it already passes, so no stream moves. If that
    /// ever stops holding, adopting it becomes a golden regeneration instead
    /// of a no-op — which is the whole reason it is reserved now.
    #[test]
    fn packed_index_is_injective() {
        assert_eq!(Stream::packed_index(0, 0), 0, "the V0 coordinate moved");
        let corners = [0_u32, 1, 2, 0xFFFF, 0xFFFF_FFFE, u32::MAX];
        let mut seen = std::collections::BTreeSet::new();
        for major in corners {
            for minor in corners {
                assert!(
                    seen.insert(Stream::packed_index(major, minor)),
                    "packed_index({major}, {minor}) collides with an earlier pair"
                );
            }
        }
        // The two arguments are not symmetric. This pins that they are
        // *distinguished*, not which one takes the high word — swapping the
        // halves is still injective and still sends `(0, 0)` to 0, so nothing
        // here or downstream can tell, and the doc claims no more than that.
        assert_ne!(Stream::packed_index(1, 0), Stream::packed_index(0, 1));
    }

    /// The one value the sub-coordinate offset costs, pinned rather than left
    /// to be discovered. Sixty-four bits cannot hold 2⁶⁴ sub-streams *and* a distinct
    /// root, so some aliasing is unavoidable; stating which is better than an
    /// unqualified claim. See [`Stream::sub`].
    #[test]
    fn sub_of_u64_max_aliases_the_root() {
        let mut aliased = Stream::sub(1, Domain::Fold, 0, u64::MAX);
        let mut root = Stream::new(1, Domain::Fold, 0);
        assert_eq!(aliased.next_u64(), root.next_u64());
    }

    /// Ranges chosen to *break* the affine form, not to pass. `(-3.5, 2.25)`
    /// — the only range the first draft tested — is **provably immune on all
    /// 2⁵³ grid points**, so it would have passed against an implementation
    /// that returned `hi` half the time. The narrow-at-large-offset rows are
    /// the discriminator; the wide row is kept only to show it still works.
    ///
    /// Keep that number. It is what stops the array being trimmed back to the
    /// immune case in a later tidy-up — and this doc was itself deleted once,
    /// in a commit that had nothing to do with it.
    #[test]
    fn floats_in_a_range_are_bounded() {
        const RANGES: &[(f64, f64)] = &[
            (1.0, 1.0 + f64::EPSILON),
            (1000.0, 1000.0 + 1e-9),
            (0.1, 0.1 + f64::EPSILON),
            (1e16, 1.000_000_000_000_000_2e16),
            (-1e308, 1e308),
            (0.0, f64::MAX),
            (-3.5, 2.25),
            (180.0, 260.0),
        ];
        for &(lo, hi) in RANGES {
            let mut s = Stream::new(5, Domain::Beaker, 0);
            for _ in 0..20_000 {
                let v = s.next_f64_range(lo, hi);
                assert!(v.is_finite(), "non-finite in [{lo}, {hi}): {v}");
                assert!(v >= lo && v < hi, "out of [{lo}, {hi}): {v}");
            }
        }
    }

    /// Sampling cannot reach either end of the grid, so the boundaries are
    /// checked at the values, through the real function.
    ///
    /// **This is a spot check, not a proof of coverage.** An earlier version
    /// claimed `u` was monotone in the result and that the extremes therefore
    /// bracketed all 2⁵³ draws. `range_from`'s own doc retracts that with
    /// measurements — `(180.0, 260.0)`, a range this very test uses, descends
    /// at 18.75% of consecutive grid steps near the top. The retraction landed
    /// in one copy and not this one, which is the third time that has happened
    /// in this file.
    #[test]
    fn the_range_endpoints_hold_at_the_extreme_draws() {
        // The largest and smallest values next_f64 can return.
        #[expect(
            clippy::as_conversions,
            clippy::cast_precision_loss,
            reason = "same 53-bit-exact conversion next_f64 itself performs, reproduced here to \
                      name the extreme draw rather than to sample it"
        )]
        let u_max = ((u64::MAX >> 11) as f64) / 9_007_199_254_740_992.0;
        for &(lo, hi) in &[
            (1.0_f64, 1.0 + f64::EPSILON),
            (1e300, 1.000_000_000_000_000_2e300),
            (1000.0, 1000.0 + 1e-9),
            (0.1, 0.1 + f64::EPSILON),
            (-1e308, 1e308),
            (0.0, f64::MAX),
            (-3.5, 2.25),
        ] {
            for u in [0.0, 1.0 / 9_007_199_254_740_992.0, 0.5, u_max] {
                let v = range_from(u, lo, hi);
                assert!(v.is_finite(), "non-finite at u = {u} in [{lo}, {hi}): {v}");
                assert!(v >= lo && v < hi, "u = {u} in [{lo}, {hi}) gave {v}");
            }
        }
    }

    /// The documented `|lo| > 2⁻⁹⁷⁰` qualifier, pinned at the value that
    /// violates it. Below that threshold `hi * u` underflows into the
    /// subnormals, its rounding error becomes absolute rather than relative,
    /// and the result lands one ulp below `lo`.
    ///
    /// Asserted as the *known* behaviour rather than as a bound that holds:
    /// this is out of range, deliberately documented, and unreachable in
    /// Borbax — nothing here goes near 1e-293. It is pinned so that a future
    /// change which silently fixes or worsens it is visible.
    #[test]
    fn the_range_bound_fails_below_the_documented_magnitude() {
        let lo = f64::from_bits(0x0020_0000_0000_0001);
        let hi = lo.next_up();
        let u = 1.0 / 9_007_199_254_740_992.0;
        let v = range_from(u, lo, hi);
        assert!(
            v < lo,
            "the subnormal underflow case no longer reproduces — if this is a \
             deliberate fix, the |lo| > 2^-970 qualifier on next_f64_range's \
             doc must be removed with it. got {v:e}, lo {lo:e}"
        );
    }

    /// The convex combination must not move the common case. `lo == 0.0` is
    /// the shape a Gillespie channel selection over `[0, total_propensity)`
    /// would use, and there the two formulations agree bit-for-bit.
    ///
    /// **The first draft of this test did not call `range_from` at all** — it
    /// asserted an algebraic identity between two inline expressions and would
    /// have passed with the function deleted. Confirmed: it survived every
    /// degenerate mutation, including "always return `lo`". Same shape as the
    /// `fork_zero_is_not_the_parents_first_draw` first draft this crate had to
    /// repair before `fork` was removed.
    ///
    /// The equality also needs `hi > f64::MIN_POSITIVE`, which the original
    /// claim omitted: below that the affine form returns `hi` and the clamp
    /// correctly fires, so the two *should* differ. Measured — at `hi = 1e-320`
    /// they differ on 3 000 000 of 3 000 000 draws.
    #[test]
    fn a_zero_lower_bound_matches_the_affine_form_bitwise() {
        let mut s = Stream::new(21, Domain::Beaker, 0);
        for &hi in &[12.5_f64, 1.0, f64::MAX, 1e-300, 2.0 * f64::MIN_POSITIVE] {
            for _ in 0..20_000 {
                let u = s.next_f64();
                let affine = 0.0 + u * (hi - 0.0);
                assert_eq!(
                    range_from(u, 0.0, hi).to_bits(),
                    affine.to_bits(),
                    "hi = {hi:e}, u = {u}"
                );
            }
        }
    }

    /// **The finding this crate's tests most needed.** The two bound tests pin
    /// only "finite and in `[lo, hi)`", and `lo` satisfies that — so a
    /// `range_from` that ignores `u` entirely passes the whole suite. Measured
    /// against six degenerate mutations (always `lo`, always `next_down(hi)`,
    /// always the midpoint, weights swapped, clamped affine, clamp to `lo`):
    /// **all six passed 28 of 28 tests.** The clamped-affine one is the
    /// plausible refactor, and through the real stream it yields *one distinct
    /// value in 20 000 draws* over `(-1e308, 1e308)`.
    ///
    /// `next_u64` has a golden and `next_f64` is fully specified against it;
    /// these two had nothing. A golden is the house style here for exactly this
    /// reason, and it is what makes the difference between a property that
    /// happens to hold and the function actually being the one intended.
    ///
    /// A failure here means the same thing as any other moved golden: the
    /// physics changed.
    #[test]
    fn ranged_and_normal_draws_are_pinned() {
        let mut s = Stream::new(0x5EED, Domain::Beaker, 0);
        let got = [
            s.next_f64_range(-3.5, 2.25).to_bits(),
            s.next_f64_range(0.0, 1.0).to_bits(),
            s.next_f64_range(180.0, 260.0).to_bits(),
            s.next_normal().to_bits(),
            s.next_normal().to_bits(),
        ];
        assert_eq!(got, GOLDEN_RANGE_AND_NORMAL);
    }

    const GOLDEN_RANGE_AND_NORMAL: [u64; 5] = [
        4_602_769_819_129_994_958,
        4_605_768_476_265_832_909,
        4_642_522_722_373_854_909,
        4_599_707_269_504_212_671,
        4_603_032_638_160_435_161,
    ];

    /// The golden pins the values; this pins that they are not all the *same*
    /// value. A degenerate implementation that happened to match one golden
    /// entry would still be caught here, and this is the assertion that
    /// separates the convex form from the clamped-affine one — which collapses
    /// a wide range to a single output.
    #[test]
    fn a_wide_range_does_not_collapse_to_one_value() {
        let mut s = Stream::new(3, Domain::Beaker, 0);
        let mut seen = std::collections::BTreeSet::new();
        for _ in 0..20_000 {
            seen.insert(s.next_f64_range(-1e308, 1e308).to_bits());
        }
        assert!(
            seen.len() >= 19_000,
            "only {} distinct values in 20000 draws — the range sampler has collapsed",
            seen.len()
        );
    }

    /// Non-finite bounds produce a plausible value rather than a loud failure.
    /// Pinned so it is a documented decision: a NaN `total_propensity` reaching
    /// a Gillespie channel selection would otherwise silently select channel 0
    /// for the rest of the run.
    ///
    /// **The `lo` rows were missing and a non-finite `lo` does *not* behave
    /// like a non-finite `hi`.** `lo = -inf` with a finite `hi` returns `-inf`
    /// on every draw, because the fallback returns `lo`. That is asserted here
    /// as known behaviour rather than as a bound that holds — there is no
    /// honest value to return — and it is why the doc on
    /// [`Stream::next_f64_range`] scopes its swallowing claim to `hi`.
    #[test]
    fn non_finite_bounds_are_swallowed_not_propagated() {
        // A non-finite `hi` is absorbed.
        assert_eq!(range_from(0.5, 0.0, f64::NAN).to_bits(), 0.0_f64.to_bits());
        for u in [0.0, 0.5, 1.0 - f64::EPSILON] {
            assert_eq!(
                range_from(u, 0.0, f64::INFINITY).to_bits(),
                f64::MAX.to_bits()
            );
        }
        // A non-finite `lo` is not: the fallback hands it straight back.
        assert_eq!(range_from(0.5, f64::NAN, 1.0).to_bits(), f64::NAN.to_bits());
        for u in [0.0, 0.5, 1.0 - f64::EPSILON] {
            // `.to_bits()`, not `is_infinite()`: the first draft accepted
            // `+inf` too, which would violate the function's own `[lo, hi)`
            // postcondition. A planted mutation returning `+inf` passed all
            // 38 tests. Pin what the doc claims, exactly.
            assert_eq!(
                range_from(u, f64::NEG_INFINITY, 1.0).to_bits(),
                f64::NEG_INFINITY.to_bits(),
                "if a non-finite lo now yields anything else, the doc on \
                 next_f64_range must lose its asymmetry caveat"
            );
        }
        // ...except against an infinite `hi`, where the endpoint clamp catches
        // the NaN that `-inf + inf` produces.
        assert_eq!(
            range_from(0.5, f64::NEG_INFINITY, f64::INFINITY).to_bits(),
            f64::MAX.to_bits()
        );
    }

    /// `next_normal` needs a distributional check, not just a support bound.
    /// A clamp floor of 0.5 collapses half of all draws and is caught today
    /// only by a bound assertion; nothing looks at the shape.
    #[test]
    fn normals_have_the_right_first_two_moments() {
        let mut s = Stream::new(0x00_1D, Domain::Universe, 0);
        let n = 400_000;
        let mut sum = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        for _ in 0..n {
            let z = s.next_normal();
            sum += z;
            sum_sq += z * z;
        }
        let mean = sum / f64::from(n);
        let var = sum_sq / f64::from(n) - mean * mean;
        assert!(mean.abs() < 0.01, "mean was {mean}");
        assert!((var - 1.0).abs() < 0.02, "variance was {var}");
    }

    /// The draw count must not depend on the arguments, or stream position
    /// becomes a function of caller data and §13.2 replay cannot restore it.
    #[test]
    fn next_f64_range_always_consumes_exactly_one_draw() {
        for &(lo, hi) in &[(0.0_f64, 1.0_f64), (5.0, 1.0), (2.0, 2.0), (f64::NAN, 1.0)] {
            let mut a = Stream::new(31, Domain::Beaker, 0);
            let _ = a.next_f64_range(lo, hi);
            let after = a.next_u64();
            let mut b = Stream::new(31, Domain::Beaker, 0);
            let _ = b.next_f64();
            assert_eq!(after, b.next_u64(), "draw count moved for [{lo}, {hi})");
        }
    }

    /// A reversed or empty range has no honest answer; it returns `lo`. Pinned
    /// so it cannot drift into returning something below `lo`, which the
    /// endpoint clamp would otherwise do for `lo == hi`.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "CLAUDE.md: tests may assert exactly. Returning `lo` bit-for-bit is the \
                  contract under test, so an epsilon would defeat it"
    )]
    fn a_reversed_or_empty_range_returns_lo() {
        let mut s = Stream::new(9, Domain::Beaker, 0);
        assert_eq!(s.next_f64_range(5.0, 1.0), 5.0);
        assert_eq!(s.next_f64_range(2.0, 2.0), 2.0);
        assert_eq!(s.next_f64_range(-1.5, -1.5), -1.5);
    }

    /// The property Task 4 depends on: `Mass::from_f64_quantised` rejects
    /// non-finite input, and `mass = base + next_normal() * spread` is the
    /// only realistic way it would ever see one. If this holds, that `Err`
    /// arm is unreachable in the generator.
    #[test]
    fn normals_are_always_finite() {
        let mut s = Stream::new(0xB0_1D, Domain::Universe, 0);
        for _ in 0..200_000 {
            let v = s.next_normal();
            assert!(v.is_finite(), "non-finite normal: {v}");
        }
    }

    /// The boundary the loop above cannot reach. `next_f64` returns 0.0 with
    /// probability 2^-53, so sampling will not find it this side of the heat
    /// death — the clamp has to be tested at the value directly.
    ///
    /// Deleting the clamp makes this fail and leaves every other test in this
    /// file passing: `ln(0)` is `-inf`, `sqrt(inf)` is `inf`, and `inf` times
    /// a cosine is `±inf` for every `u2` that does not land exactly on a
    /// zero of the cosine.
    #[test]
    fn the_normal_transform_is_finite_at_u1_zero() {
        for step in 0..64_u32 {
            let u2 = f64::from(step) / 64.0;
            let v = normal_from(0.0, u2);
            assert!(v.is_finite(), "non-finite at u1 = 0, u2 = {u2}: {v}");
        }
    }

    /// Finiteness alone does not pin the clamp — it passed against
    /// `f64::MIN_POSITIVE`, which put a lone atom of probability at 37.64σ
    /// while every reachable input topped out at 8.5717σ. What distinguishes a
    /// real clamp from a merely-total one is that `u1 == 0` is
    /// *indistinguishable from the edge of the reachable support*.
    ///
    /// **Discriminator, with its resolution.** This fails at any floor more
    /// than about 1e-14 relatively *below* 2⁻⁵³ — verified at
    /// `f64::MIN_POSITIVE`, `1e-300` and `2⁻⁵⁴`. Two blind spots, both
    /// measured, both covered elsewhere:
    ///
    /// - Floors **above** 2⁻⁵³ pass vacuously — at `f64::EPSILON`, `1e-3` and
    ///   `0.9` both arguments clamp to the same value and the equality is
    ///   trivial. `the_normal_support_is_bounded_at_the_documented_sigma`
    ///   catches those, which is why neither test is redundant.
    /// - The nearest **~90 representable floats** below 2⁻⁵³ also pass, because
    ///   `sqrt(-2 ln u)` compresses a relative change in `u` by roughly
    ///   1/(2·8.57): a 1-ulp move in the floor is ~2.6e-17 against an ulp of
    ///   1.78e-15 at 8.5717. Nothing plausible lands there.
    ///
    /// Two earlier versions of this sentence overstated it — first "any other
    /// floor", then "any floor below". Both were written to repair an
    /// overstatement of exactly this kind.
    ///
    /// Floors *above* 2⁻⁵³ are caught by
    /// `the_normal_support_is_bounded_at_the_documented_sigma` — a floor of 0.9
    /// collapses the normal to `|z| ≤ 0.459` and only that test sees it. The
    /// two together pin the constant; **neither is redundant**, which is worth
    /// saying because a reader trusting the old sentence would delete the
    /// second one.
    #[test]
    fn the_normal_clamp_lands_on_the_edge_of_the_reachable_support() {
        let smallest_reachable = 1.0 / 9_007_199_254_740_992.0;
        for step in 0..64_u32 {
            let u2 = f64::from(step) / 64.0;
            assert_eq!(
                normal_from(0.0, u2).to_bits(),
                normal_from(smallest_reachable, u2).to_bits(),
                "the clamp is not at 2^-53 — a gap in the support has opened at u2 = {u2}"
            );
        }
    }

    /// The support bound, pinned as a number so the doc cannot rot away from
    /// the code. Anything downstream reasoning about `base / spread` to keep a
    /// generated mass positive is reasoning against this constant.
    #[test]
    fn the_normal_support_is_bounded_at_the_documented_sigma() {
        let smallest_reachable = 1.0 / 9_007_199_254_740_992.0;
        let max_radius = normal_from(smallest_reachable, 0.0);
        assert!(
            (max_radius - 8.571_674_348_652_905).abs() < 1e-12,
            "max |z| moved: {max_radius}"
        );
        let mut s = Stream::new(0xB0_1D, Domain::Universe, 0);
        for _ in 0..200_000 {
            assert!(s.next_normal().abs() <= max_radius);
        }
    }

    /// Clamping and redrawing both make `next_normal` total, and they are not
    /// interchangeable: rejecting consumes further draws and shifts every
    /// subsequent value in the domain's stream, which regenerates the whole
    /// universe. Clamping was chosen. From this commit the draw count is part
    /// of the determinism contract, so it is pinned rather than assumed.
    #[test]
    fn next_normal_consumes_exactly_two_draws() {
        let mut normal = Stream::new(11, Domain::Universe, 0);
        let _ = normal.next_normal();
        let after_normal = normal.next_u64();

        let mut plain = Stream::new(11, Domain::Universe, 0);
        let _ = plain.next_u64();
        let _ = plain.next_u64();
        assert_eq!(
            after_normal,
            plain.next_u64(),
            "next_normal's draw count moved — every golden in the workspace moves with it"
        );
    }

    /// Golden values. If the mixer changes, every golden run hash in the
    /// project changes with it — so this is pinned deliberately and a
    /// failure here means "you just changed all the physics", not
    /// "the test is stale".
    #[test]
    fn golden_sequence_is_pinned() {
        let mut s = Stream::new(0, Domain::Universe, 0);
        let got: Vec<u64> = (0..12).map(|_| s.next_u64()).collect();
        assert_eq!(got, GOLDEN_UNIVERSE_0);
    }

    /// **The length was load-bearing, and it no longer is — the job moved.**
    /// Philox emits four `u64` per block, so a four-draw golden sits inside
    /// block 0 and pins nothing about how the block counter advances; the
    /// length went 4 → 8 → 12 chasing that, because eight words is two block
    /// *fills* but only one *transition*, and `counter[3] ^= 1` (period 8) and
    /// `counter[3] = 2c+1` both agree with the real generator through word 7.
    /// Twelve covers blocks 0, 1 and 2 — and a reviewer then found a survivor
    /// even at twelve, `counter[3] += 1 + counter[2]`, which no length of
    /// golden can reach because it is identical for every root stream.
    ///
    /// `draws_are_the_direct_block_function_of_their_position` now pins the
    /// advance for **every** `k`, with no constants at all, which is what
    /// ended that chase. Measured: with this test `#[ignore]`d, a 29-mutation
    /// battery has an **identical** catch set — this golden is the sole
    /// catcher of nothing.
    ///
    /// **Kept anyway, and for a reason that is not inertia.** It is a second,
    /// independently-recorded evaluation of the kernel: the direct-block test
    /// checks the crate against `philox4x64_10`, and if that function is wrong
    /// they agree with each other. These twelve words were reproduced from an
    /// independent Python re-derivation of the paper. That is the same role
    /// `matches_the_published_known_answer_vectors` plays one level down, and
    /// it costs twelve lines. Do not extend it; there is nothing left for
    /// length to buy.
    const GOLDEN_UNIVERSE_0: [u64; 12] = [
        213_000_021_201_967_259,
        4_455_796_210_202_625_458,
        2_055_444_239_878_205_049,
        10_411_612_076_246_414_556,
        6_312_158_256_571_094_726,
        10_634_814_581_434_429_480,
        1_598_446_939_479_630_672,
        11_723_492_092_571_950_057,
        16_006_495_525_962_256_124,
        11_263_462_740_154_766_164,
        18_401_285_607_914_570_996,
        16_002_093_844_216_690_389,
    ];

    /// Adding a `Domain` variant must not perturb any existing stream.
    ///
    /// These are golden values, and the point is that they can only change
    /// deliberately. If appending `Domain::Field = 10` moves any of them, the
    /// enum's discriminants have been renumbered or `Stream::new`'s mixing has
    /// changed — either of which regenerates the chemistry of every universe
    /// ever produced. A reviewer seeing this test fail should treat it exactly
    /// like a moved golden hash.
    ///
    /// **What the exhaustive `match` enforces, stated precisely.** Appending a
    /// variant makes it non-exhaustive, so the crate stops compiling until the
    /// new variant is given its own pinned triple. That is a compile-time
    /// gate, not a comment. What it does *not* enforce is that the new variant
    /// also reaches [`Domain::ALL`] — if it is omitted there, the new stream
    /// goes unexercised while all nine existing ones stay pinned, which is the
    /// property this test exists to protect.
    #[test]
    fn every_domain_stream_is_pinned() {
        for domain in Domain::ALL {
            let want: [u64; 3] = match domain {
                Domain::Universe => [
                    3_929_780_809_458_279_810,
                    14_808_739_883_066_402_667,
                    6_821_936_907_122_725_605,
                ],
                Domain::Naming => [
                    1_829_254_752_149_997_246,
                    5_728_333_709_296_349_463,
                    9_288_166_565_660_478_743,
                ],
                Domain::Molecule => [
                    4_049_216_183_878_968_044,
                    15_790_089_183_108_422_852,
                    16_750_542_454_782_617_201,
                ],
                Domain::Fold => [
                    2_306_140_887_577_208_375,
                    11_691_205_199_486_685_485,
                    2_911_623_241_027_970_606,
                ],
                Domain::Reaction => [
                    17_766_946_860_710_959_347,
                    6_927_189_939_051_538_356,
                    13_202_956_210_144_916_749,
                ],
                Domain::Decay => [
                    4_064_337_356_115_406_327,
                    2_355_560_016_391_046_765,
                    4_745_691_913_722_414_085,
                ],
                Domain::Beaker => [
                    12_864_929_427_389_045_100,
                    15_550_989_763_699_698_652,
                    13_012_296_228_043_260_465,
                ],
                Domain::Shadow => [
                    721_593_208_716_469_688,
                    16_408_411_480_340_345_767,
                    14_942_133_425_965_051_356,
                ],
                Domain::Hash => [
                    543_525_597_416_629_224,
                    15_581_591_749_131_607_509,
                    11_412_060_199_028_555_773,
                ],
            };
            let mut s = Stream::new(0x5EED, domain, 0);
            let got = [s.next_u64(), s.next_u64(), s.next_u64()];
            assert_eq!(got, want, "{domain:?} stream moved");
        }
    }

    /// The property the pinning test protects, stated directly: a stream is a
    /// pure function of (seed, domain, index), so no domain can influence
    /// another regardless of how many draws it takes.
    #[test]
    fn domains_are_independent_of_each_others_draw_counts() {
        let mut fold = Stream::new(0x5EED, Domain::Fold, 0);
        for _ in 0..1000 {
            let _ = fold.next_u64();
        }
        let mut a = Stream::new(0x5EED, Domain::Universe, 0);
        let mut b = Stream::new(0x5EED, Domain::Universe, 0);
        assert_eq!(a.next_u64(), b.next_u64());
    }

    /// "Append; never renumber, never reuse" is the rule the whole
    /// `Domain`-separation argument rests on, and this asserts *that* rule
    /// rather than a stronger one.
    ///
    /// **The first draft asserted contiguity from 1, which contradicts the
    /// rule it was protecting.** Retiring a subsystem leaves a gap — that is
    /// what "retired, not recycled" means — and a contiguity test fails on a
    /// gap with the message "or its discriminant was renumbered". The next
    /// person reads that, concludes the test is stale, and renumbers to close
    /// the gap, recycling a retired discriminant and silently regenerating
    /// every universe. The test would have caused the exact failure it was
    /// written to prevent.
    ///
    /// Strictly increasing is the real invariant: gaps are legal, reordering
    /// and reuse are not. (Reuse is additionally a compile error — `rustc`
    /// rejects duplicate discriminants with `E0081` — so what this test
    /// uniquely catches is the out-of-order case.)
    #[test]
    fn discriminants_are_strictly_increasing_in_declaration_order() {
        let mut previous = 0_u64;
        for domain in Domain::ALL {
            let d = domain.discriminant();
            assert!(
                d > previous,
                "{domain:?} has discriminant {d}, not greater than the previous {previous} — \
                 variants may only be appended, and a retired discriminant is never reused. \
                 A gap is legal; going backwards or repeating is not."
            );
            previous = d;
        }
    }

    /// `Ord` exists so `BTreeMap<Domain, _>` compiles — without it the API
    /// steers callers to the `HashMap` §13.1 bans. Deriving it is not enough:
    /// the order must agree with the discriminant order the append-only
    /// argument rests on, which a derive placed above a reordered enum would
    /// not.
    #[test]
    fn ordering_agrees_with_discriminant() {
        // `first() < last()` rather than `w[0] < w[1]`: same comparison on a
        // 2-window, without an indexing operation that could panic.
        assert!(Domain::ALL.windows(2).all(|w| w.first() < w.last()));
        let sorted_by_ord = {
            let mut v = Domain::ALL;
            v.sort_unstable();
            v
        };
        let sorted_by_discriminant = {
            let mut v = Domain::ALL;
            v.sort_unstable_by_key(|d| d.discriminant());
            v
        };
        assert_eq!(sorted_by_ord, sorted_by_discriminant);
    }

    /// The `#[expect]` on `next_f64` claims its conversion is lossless. The
    /// first draft cited `floats_are_in_unit_interval` as the proof, which
    /// asserts only the range — it would pass against a non-power-of-two
    /// `SCALE` that made every multiply inexact. This is the assertion that
    /// claim actually needs.
    /// **What actually fails to discriminate is the range check, and a
    /// `trunc()` check is not in that category.** Measured against
    /// `SCALE = 1/(2⁵³+2)` — the nearest inexact scale the f64 grid contains —
    /// over 200,000 draws: the range assertion passes on all of them, while
    /// `(v * 2⁵³).trunc() == v * 2⁵³` **fails on 75,034**, so a trunc check
    /// would catch that mutation inside the first block. Two earlier versions
    /// of this sentence said the opposite. The first cited `1/(2⁵³+1)`, and
    /// `2⁵³+1` is not representable — that literal *is* `2⁻⁵³`, so the
    /// mutation was a no-op and a reader reproducing it saw "survived" and
    /// would have concluded the test was weak. The second fixed the literal
    /// but kept the claim, which was false for the different reason above.
    /// It is written out rather than deleted because "a trunc check is
    /// inadequate here" is precisely the kind of sentence that outlives its
    /// evidence and gets cited to reject a perfectly adequate assertion
    /// somewhere else.
    ///
    /// The bit comparison is still the right test, for a reason that never
    /// depended on any of that: it computes the expected value by a *different
    /// operation* — exact division by a power of two — rather than re-deriving
    /// it from the value under test.
    #[test]
    fn floats_land_exactly_on_the_53_bit_grid() {
        let mut raws = Stream::new(77, Domain::Reaction, 0);
        let mut floats = raws.clone();
        for _ in 0..200_000 {
            let raw = raws.next_u64() >> 11;
            let got = floats.next_f64();
            #[expect(
                clippy::as_conversions,
                clippy::cast_precision_loss,
                reason = "raw is below 2^53 by the shift, so this conversion is exact — which is \
                          the property under test"
            )]
            let want = raw as f64 / 9_007_199_254_740_992.0;
            assert_eq!(got.to_bits(), want.to_bits(), "raw {raw} gave {got:e}");
        }
    }

    /// A sub-stream takes no parent, so the parent's draw position cannot
    /// influence it — the path-indexed property the old `fork` documented, now
    /// true by construction rather than by care. Kept as a regression test
    /// because a future `sub` that consulted a parent would reintroduce
    /// draw-order dependence, which §13.1 forbids.
    #[test]
    fn a_sub_stream_is_unaffected_by_how_far_the_root_has_advanced() {
        let early = Stream::sub(1, Domain::Fold, 0, 3).next_u64();
        let mut root = Stream::new(1, Domain::Fold, 0);
        for _ in 0..1000 {
            let _ = root.next_u64();
        }
        assert_eq!(Stream::sub(1, Domain::Fold, 0, 3).next_u64(), early);
    }
}
