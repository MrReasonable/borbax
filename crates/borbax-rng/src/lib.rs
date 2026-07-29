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
//! 64-bit field — `[domain, index, fork, block]` — so distinct coordinates are
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
//! **Five was a local fluctuation, not structure, and the honest number
//! belongs here** — an earlier draft of this paragraph dropped the correction
//! and kept only the alarming half. The expectation for distinct relations is
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
//! 1 is reserved for the world seed, so §13.1's `(universe_seed, world_seed,
//! config_hash)` tuple never has to be compressed into one word.

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
//   §13.1 hazard this crate actually trips — `ln`/`cos` in `normal_from` — is
//   NOT protected here. `clippy::disallowed_methods` remains the sole
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
/// meaning what they meant (§13.3), and nothing would announce it.
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
    /// **Its earlier justification — "must not consume the live run's draws" —
    /// was vacuous, and the correction matters.** A `Stream` is a value with
    /// its own counter and `next_u64` takes `&mut self`, so a shadow drawing
    /// from its own `Stream` is *arithmetically incapable* of advancing the
    /// live run's, whatever domain it uses. That guarantee is delivered by the
    /// counter-based construction for every domain and buys this variant
    /// nothing.
    ///
    /// What this variant actually does is make the shadow's draws *different*
    /// from the live run's, and whether that is wanted is an open question,
    /// not a settled one. §15.3 makes every emergence claim a `run − shadow`
    /// difference, and a paired difference estimated from one run and one
    /// shadow *can* have far lower variance under common random numbers than
    /// under independent streams — "can", because Glasserman and Yao (1992)
    /// find the class of systems for which CRN is *provably* advantageous
    /// "rather limited". Against that, common random numbers decouple in a
    /// Gillespie setting as the integrated intensities diverge; indexing draws
    /// per reaction channel (the common-reaction-path method) *reduces but
    /// does not remove* that, and Anderson (2012) shows it decouples too.
    ///
    /// **And the whole literature is for infinitesimal parameter
    /// perturbations.** `run − shadow` is selection on versus selection off, a
    /// large structural difference — the regime where every one of these
    /// couplings decays fastest. The pessimistic side of this is understated
    /// above, not overstated.
    ///
    /// §15.3 also defines a *second* control that this variant will be read as
    /// covering and does not: randomising which molecule catalyses which
    /// reaction at fixed catalysis density. That is a different experiment
    /// with the opposite requirement — it needs a permutation that does *not*
    /// move when the live run's draw count changes — and may want its own
    /// domain. Whatever Task 15 decides must then apply to *every* draw on the
    /// shadow path; mixing `Shadow` and `Beaker` draws gives partial CRN,
    /// silently, which is the worst of the options.
    ///
    /// **Decide it in Task 15/20; do not read this variant as having decided
    /// it.** What would settle it is cheap: run both indexings at fixed budget
    /// and compare the variance of the paired difference.
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
#[derive(Debug, Clone)]
pub struct Stream {
    /// `[universe_seed, reserved]`. The second word is deliberately unused
    /// and deliberately present: §13.1's reproducibility tuple is
    /// `(universe_seed, world_seed, config_hash)`, and a 64-bit key would
    /// force those to be compressed together — which is exactly the hashing
    /// that produced permanent stream collisions in the first design.
    key: [u64; 2],
    /// `[domain, index, fork, block]` — one 64-bit field each, no packing.
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
    /// between the two. [`Stream::fork`] has two documented exceptions of its
    /// own.
    ///
    /// This is still the property the first design did not have, and its
    /// absence was not theoretical: hashing the same coordinates into a
    /// 64-bit key produced five permanent seed-independent collisions below
    /// index 2³⁰.
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

    /// Derive a child stream. Used where a fixed number of sub-streams is
    /// needed — one per fold attempt, one per species — without threading a
    /// mutable parent through the call graph.
    ///
    /// **Not one per molecule at the bulk or stochastic tiers.** To call
    /// `fork(i)` per molecule you need a stable `i` per molecule, and stable
    /// means stored — which is per-molecule state arriving disguised as "just
    /// an index", the shape §8.6's regression always takes. §9.5 makes all
    /// molecules of a species interchangeable precisely so none needs
    /// individual state, and a per-molecule index turns a keyframe from
    /// O(species) into O(molecules). V1's molecular tier (§12.1) is the one
    /// place a per-molecule stream is legitimate; V0 has no per-molecule
    /// representation at all.
    ///
    /// **`fork` is a pure function of the parent's key, not of its state.**
    /// It deliberately ignores the parent's counter — the path-indexed style,
    /// which is the right choice for a §13.1 that forbids results depending on
    /// draw order. The cost is that `parent.fork(0)` called at two different
    /// points returns the *same* child, so two call sites forking the same
    /// parent with overlapping index ranges silently share a sequence. That is
    /// a duplicate stream arriving through the one door `!Copy` does not
    /// guard. Give each call site a disjoint index range, or a domain.
    ///
    /// **Forking is one level deep. That limit is inherent, and a second fork
    /// is now rejected rather than silently honoured.**
    ///
    /// The child writes `index` into the counter's third field, so a fork of a
    /// fork *replaces* its parent's fork coordinate instead of extending it —
    /// measured, `base.fork(a).fork(b)` was bit-identical to `base.fork(b)`
    /// for every `a`, which is a duplicate stream and precisely the hazard
    /// this whole design removes. It cannot be made to work: there is exactly
    /// one spare counter word, and hashing a fork *path* into it would
    /// reintroduce the collisions Philox was adopted to eliminate.
    ///
    /// One level is enough for §13.1, which asks for a stream per subsystem,
    /// per patch, per time index — that is `domain`, `index`, `fork`, with
    /// `block` left as the draw counter. A second level needs its own
    /// `Domain`, or its two indices composed before the first fork.
    ///
    /// **Type-state was considered and rejected on cost.** Making the second
    /// call a compile error means `Stream<const FORKED: bool>`, which makes
    /// every downstream function that merely *draws* generic over fork depth —
    /// viral across eight crates, to forbid a call nobody has yet written. The
    /// `debug_assert` catches it in any test run; `forking_a_fork_is_rejected`
    /// explains why that is the right trade and what it does not buy.
    ///
    /// **Fork coordinates are offset by one, because the root occupies zero.**
    /// The first draft wrote `index` directly, so `fork(0)` produced the
    /// parent's own counter — a duplicate stream, which is the precise hazard
    /// this whole design exists to remove. Its own test caught it immediately.
    /// The offset costs one value: `fork(u64::MAX)` wraps back onto the root.
    /// Sixty-four bits cannot hold 2⁶⁴ forks *and* a distinct root, so some
    /// aliasing is unavoidable; stating which is better than an unqualified
    /// claim. Pinned in `fork_of_u64_max_aliases_the_root`.
    #[inline]
    #[must_use]
    pub const fn fork(&self, index: u64) -> Self {
        // A root has `counter[2] == 0`; a fork never does (the offset below is
        // why). So this fires exactly when someone forks a fork — which cannot
        // be made to work, and until now failed by silently returning a
        // duplicate of an unrelated sibling.
        debug_assert!(
            self.counter[2] == 0,
            "fork is one level deep: forking a fork would overwrite the fork \
             coordinate rather than extend it, silently duplicating another \
             stream. Give the second level its own Domain, or compose the two \
             indices before the first fork."
        );
        Self {
            key: self.key,
            counter: [self.counter[0], self.counter[1], index.wrapping_add(1), 0],
            buf: [0; 4],
            spent: 4,
        }
    }

    /// The next raw 64 bits.
    ///
    /// **Cost, both shapes, because reporting only one of them misleads.**
    /// Against the mixer this replaced, release: 1.544 → 1.639 ns in a
    /// 1000-draw accumulator loop (**+6.2%**), and 79.2 → 96.7 µs over a
    /// 20 000-step fold anneal (**+22%**). Put a `black_box` on every single
    /// draw and it reads 1.37 → 3.52 ns (2.6×) — a correct measurement of a
    /// serialised pipeline that no caller doing work between draws will pay,
    /// and the number this commit originally reported on its own. The kernel
    /// sits at 94.5% of the hardware multiplier-port limit, so there is
    /// essentially nothing left to win here.
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
    /// argument the way `normal_from` clamps `u1`. Task 15 owns this.
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
        // 2^53, exact in f64. Spelled as a literal rather than `(1u64 << 53)
        // as f64` so the scale factor needs no cast of its own.
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
            if v < zone {
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
    /// **Non-finite bounds are swallowed, deliberately and dangerously.**
    /// `hi = NaN` returns `lo`; `hi = +inf` returns `f64::MAX` for every draw.
    /// This function is named below as the shape a Gillespie channel selection
    /// over `[0, total_propensity)` would use — so a single NaN rate anywhere
    /// in that sum makes this return `0.0` forever and the scheduler picks
    /// channel 0 for the rest of the run, with no crash and nothing in the
    /// state hash. Pinned in `non_finite_bounds_are_swallowed_not_propagated`
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
    /// is a 3–5% win against Philox** — measured in both profiles, twice each,
    /// after the rewrite. The attribute set is kept matching its evidence
    /// rather than its history.
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

    #[test]
    fn forking_is_deterministic_and_independent() {
        let base = Stream::new(1, Domain::Fold, 0);
        let mut f1 = base.fork(3);
        let mut f2 = base.fork(3);
        let mut f3 = base.fork(4);
        assert_eq!(f1.next_u64(), f2.next_u64());
        assert_ne!(f1.next_u64(), f3.next_u64());
    }

    /// Forking writes `index` into the counter's third field, so distinct
    /// fork indices are distinct counters and Philox's bijectivity does the
    /// rest. The old design needed a hand-chosen odd constant here to stop
    /// `fork(0)`'s key colliding with the parent's first output; that whole
    /// class of concern is gone, because there is no key derivation left to
    /// collide.
    #[test]
    fn forks_are_distinct_from_the_parent_and_from_each_other() {
        let base = Stream::new(1, Domain::Fold, 0);
        let mut parent = base.clone();
        let mut f0 = base.fork(0);
        let mut f1 = base.fork(1);
        let p = parent.next_u64();
        let a = f0.next_u64();
        let b = f1.next_u64();
        assert_ne!(a, p, "fork(0) collides with the parent");
        assert_ne!(a, b, "fork(0) and fork(1) collide");
    }

    /// **Which coordinate lives in which counter field is load-bearing, and
    /// four mutations of it passed all 38 tests before this existed.**
    /// Measured, each silently catastrophic and each invisible:
    ///
    /// - `fork` writing field 3: `fork(j)` *is* the parent advanced 4(j+1)
    ///   draws — `fork(0)` and `fork(1)` shared 36 of their first 40 words.
    /// - `fork` dropping `counter[1]`: every index in a domain shares one fork
    ///   family, 40 of 40 words identical.
    /// - `new` writing index into field 2: `new(s, Fold, 3)` *is*
    ///   `new(s, Fold, 0).fork(2)` — the index and fork axes become one.
    /// - `new` writing index into field 3: every index stream is a bit-exact
    ///   *suffix* of index 0's.
    ///
    /// The old fork tests all survived because they used `index = 0` and
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
                let base = Stream::new(0x5EED, domain, index);
                for fork in 0..5_u64 {
                    let mut s = if fork == 0 {
                        base.clone()
                    } else {
                        base.fork(fork - 1)
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
            "{} of {total} words are shared between distinct (domain, index, fork) streams — \
             a counter coordinate has been transposed or dropped",
            total - seen.len()
        );
    }

    /// The one value the fork offset costs, pinned rather than left to be
    /// discovered. Sixty-four bits cannot hold 2⁶⁴ forks *and* a distinct
    /// root, so some aliasing is unavoidable; stating which is better than an
    /// unqualified claim. See [`Stream::fork`].
    #[test]
    fn fork_of_u64_max_aliases_the_root() {
        let base = Stream::new(1, Domain::Fold, 0);
        let mut aliased = base.fork(u64::MAX);
        let mut root = base;
        assert_eq!(aliased.next_u64(), root.next_u64());
    }

    /// A second fork is a programming error and now says so.
    ///
    /// **The test this replaces asserted the bug as the contract.** It pinned
    /// `base.fork(3).fork(7) == base.fork(7)` and passed green, which converts
    /// a silent duplicate-stream defect into a documented guarantee — far
    /// harder to withdraw than code. Its deletion is the discriminator: a fix
    /// that leaves it passing is cosmetic.
    ///
    /// `debug_assert!` rather than a hard panic, because `clippy::panic` is
    /// `deny` in library code and a panic inside a propensity loop kills an
    /// overnight run. That is the right trade for a *caller* error — it cannot
    /// arise from data, so any test run reaches it, and CI runs the debug leg.
    /// A release build still collapses silently; the honest statement is that
    /// this catches the mistake during development rather than making it
    /// impossible.
    ///
    /// Gated on `debug_assertions` because CI also runs `--release`, where
    /// nothing panics and an ungated `should_panic` would fail.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "fork is one level deep")]
    fn forking_a_fork_is_rejected() {
        let base = Stream::new(1, Domain::Fold, 0);
        let _ = base.fork(3).fork(7);
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
    /// `fork_zero` first draft this crate already had to repair once.
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
    #[test]
    fn non_finite_bounds_are_swallowed_not_propagated() {
        assert_eq!(range_from(0.5, 0.0, f64::NAN).to_bits(), 0.0_f64.to_bits());
        assert_eq!(range_from(0.5, f64::NAN, 1.0).to_bits(), f64::NAN.to_bits());
        for u in [0.0, 0.5, 1.0 - f64::EPSILON] {
            assert_eq!(
                range_from(u, 0.0, f64::INFINITY).to_bits(),
                f64::MAX.to_bits()
            );
        }
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
    /// **Discriminator, stated accurately after measurement.** This fails at
    /// any floor *below* 2⁻⁵³ — verified at `f64::MIN_POSITIVE`, `1e-300` and
    /// `2⁻⁵⁴`. It does **not** catch floors above: at `f64::EPSILON` (2⁻⁵²),
    /// `1e-3` and `0.9` it passes vacuously, because for any floor ≥ 2⁻⁵³ both
    /// arguments clamp to the same value and the equality is trivial. An
    /// earlier version of this comment claimed it caught "any other floor",
    /// which is the same defect it was written to repair.
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
        let got: Vec<u64> = (0..8).map(|_| s.next_u64()).collect();
        assert_eq!(got, GOLDEN_UNIVERSE_0);
    }

    /// **The length is load-bearing.** Philox emits four `u64` per block, so
    /// a four-draw golden sits entirely inside block 0 and pins nothing about
    /// how the block counter advances. Measured: changing
    /// `counter[3].wrapping_add(1)` to `add(2)` — still a bijection, still
    /// statistically clean, and it changes every number in the project — was
    /// caught by exactly one test, and only because that test happened to take
    /// seven draws. Eight words spans two boundaries. Do not shorten this.
    const GOLDEN_UNIVERSE_0: [u64; 8] = [
        213_000_021_201_967_259,
        4_455_796_210_202_625_458,
        2_055_444_239_878_205_049,
        10_411_612_076_246_414_556,
        6_312_158_256_571_094_726,
        10_634_814_581_434_429_480,
        1_598_446_939_479_630_672,
        11_723_492_092_571_950_057,
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
    /// **The grid assertion alone is not enough, measured.** Setting
    /// `SCALE = 1/(2⁵³+1)` — which makes every multiply inexact — still leaves
    /// `f * 2⁵³` landing on an integer, so a `trunc()` check passes against
    /// the very defect it names. The discriminating assertion compares the
    /// produced bits against the raw draw divided by 2⁵³: division by a power
    /// of two is exact, so that is the value `next_f64` must produce, computed
    /// by a different operation than the one under test.
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

    /// `fork` ignores the parent's counter, so the same index gives the same
    /// child however far the parent has advanced. That is the correct
    /// path-indexed behaviour and also a duplicate-stream hazard, so it is
    /// pinned as a decision rather than left to be rediscovered.
    #[test]
    fn fork_is_a_function_of_the_key_not_the_parents_position() {
        let base = Stream::new(1, Domain::Fold, 0);
        let early = base.fork(3).next_u64();
        let mut advanced = base;
        for _ in 0..1000 {
            let _ = advanced.next_u64();
        }
        assert_eq!(advanced.fork(3).next_u64(), early);
    }
}
