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
//! **The key is 64 bits wide, and that is a bound worth knowing.** A stream's
//! whole output is determined by one `u64`, so distinct streams are
//! *probabilistically* rather than provably distinct — unlike Philox and
//! Threefry, whose 192-bit or wider `(key, counter)` tuple lets the
//! simulation's own indices be *embedded* rather than hashed down. (Philox4x32
//! is 64-bit key + 128-bit counter = exactly the 192 bits `Stream::new` hashes
//! into 64. An earlier version of this sentence said those families carry
//! 128-bit keys; they do not, and the correction matters because adopting one
//! would fix nothing unless the coordinates went in the *counter*.) Salmon et
//! al. put it directly: "(key, counter) tuples must never be improperly
//! re-used, since an inadvertent re-use of the same tuple will result in
//! exactly the same random number, with potentially dire consequences for
//! simulation accuracy." Expected colliding stream
//! pairs among `S` streams is about `S²/2⁶⁵`: negligible at 10⁶ streams
//! (3e-8), around 0.5 at 2³², and expected-many if anything ever opens a
//! stream per molecule per step. Two colliding streams emit *identical*
//! sequences, which nothing crashes on and no test sees. See
//! [`Stream::new`] for the concrete consequence.

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

use borbax_units::det_math;

/// Which subsystem a stream belongs to.
///
/// Adding a new random draw in one subsystem cannot shift the sequence any
/// other subsystem sees — which means a change to folding does not silently
/// perturb reaction outcomes. **That property is exact**: a stream is a pure
/// function of `(seed, domain, index)`, so no number of draws in one domain
/// can move another.
///
/// Domain *separation* is the weaker, probabilistic claim, and the difference
/// matters. See [`Stream::new`]: because all three coordinates are hashed into
/// one 64-bit key, distinct `(domain, index)` pairs can collide, and a
/// colliding pair collides in every universe ever generated.
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
    /// *species* drawn from abundance (§9.4).
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

/// `SplitMix64`'s finalizer — Stafford's Variant 13. No state, trivially
/// verifiable, and a bijection, so a stream's period is exactly 2^64 with no
/// value repeated inside it.
///
/// **Do not reduce [`Stream::next_u64`] to a single application of this.** The
/// obvious justification for the doubling — "strong avalanche" — is the wrong
/// one, and measurement says so: strict avalanche cannot tell one round from
/// two, both sitting at the noise floor over all 64×64 input/output bit pairs.
/// Avalanche is not the property doing the work.
///
/// What does the work shows up two ways. Second-order avalanche over all 2080
/// one- and two-bit input differences finds **three cells of total diffusion
/// failure** in the single round — flipping input bits 29 and 59 together
/// *never* flips output bit 0 — and **zero** in `mix(key ^ mix(c))`, whose
/// worst deviation over 133 120 cells is exactly what that many samples
/// predict. And empirically, a single finalizer over a sequential counter
/// **fails `PractRand` at 8 GB** (`BRank(12)`, p ≈ 1.7e-53), where the composed
/// form is clean to 128 GB. Measured here rather than cited: the published
/// material on this mixer family is grey literature (Evensen's *Mostly
/// Mangling*, Lemire's `testingRNG`) and reports a different setup, so it
/// corroborates the direction and not the number.
///
/// The second application costs about 0.75 ns per draw on aarch64 (M1 Pro),
/// roughly doubling the cost of `next_u64`. It is worth it, and the evidence
/// has to be stated against the right reduction or it will read as stale:
///
/// - With **stock `SplitMix64`** — one finalizer over a state advanced by the
///   golden-ratio constant — two seeds differing by `k·φ` are the same
///   sequence shifted by `k`: 4032 of 4032 constructed cases, against 0 of
///   4032 for the shipped form. That is what splitting in a loop produces.
/// - Deleting one `mix` from `next_u64` instead gives `mix(key ^ c)`, which
///   scores 0 of 4032 on *that* test and fails differently: 4032 of 4032
///   `XOR`-differing key pairs emit an exact permutation of each other's
///   draws.
///
/// Both reductions are unsafe; they are unsafe for different reasons, and a
/// reader who tries the second one against the first one's evidence will
/// wrongly conclude the comment is wrong.
#[inline]
#[must_use]
const fn mix(mut z: u64) -> u64 {
    z ^= z >> 30;
    z = z.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z ^= z >> 27;
    z = z.wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
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
    key: u64,
    counter: u64,
}

impl Stream {
    /// Open the stream for `(seed, domain, index)`.
    ///
    /// The result is a pure function of those three, which is the whole
    /// determinism argument: no stream can be influenced by how many draws
    /// any other stream has taken, so drawing order across subsystems is not
    /// part of the contract and never has to be preserved.
    ///
    /// **Distinct triples are not guaranteed distinct streams, and the seed
    /// does not help.** The key is `mix(seed ^ mix(D·φ ^ mix(index)))`, so two
    /// triples collide exactly when `mix(i₁) ^ mix(i₂) == (D₁·φ) ^ (D₂·φ)` —
    /// an equation the seed cancels out of. A colliding pair therefore
    /// collides *for every seed*, permanently, in every universe ever
    /// generated. A complete search below 2³⁰ finds five such index relations
    /// among the nine domains, two of them pairing `Decay` with `Shadow`.
    ///
    /// Five is a *local upward fluctuation*, not structure — and the arithmetic
    /// has to be stated carefully, because an earlier version of this comment
    /// called it "matching the birthday expectation" over numbers that do not
    /// match. For distinct unordered relations the expectation is
    /// `36 × C(2³⁰,2)/2⁶⁴ = 1.125`, against which five is a 4.4× excess at
    /// p = 0.006. (2.25 is the expectation for *ordered* index pairs, whose
    /// matching observation is 10, not 5.) Extending the search to 2³² gives 25
    /// relations against 18 expected, p = 0.068 — unremarkable, which is what
    /// says the mixer has no structure and the 2³⁰ window was just a bump.
    ///
    /// This is not a defect in the mixing — it is the price of hashing 192
    /// bits of coordinate into a 64-bit key, and it is why the module doc
    /// states the bound. At V0 scale it is unreachable. It stops being
    /// unreachable if `index` ever becomes a composite of patch and time step
    /// as §13.1 sketches, or a content digest as [`Domain::Hash`] invites,
    /// which is where the stream count reaches 10¹⁰–10¹². **Widening the key
    /// is a Task 13/20 decision; taking it accidentally is the thing to
    /// avoid.**
    #[inline]
    #[must_use]
    pub const fn new(seed: u64, domain: Domain, index: u64) -> Self {
        let key =
            mix(seed ^ mix(domain.discriminant().wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ mix(index)));
        Self { key, counter: 0 }
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
    /// The odd constant is real but weaker than "separation": `mix(0) == 0`,
    /// so without it `fork(0)`'s key would be bit-identical to the parent's
    /// first output. Since `mix` is a bijection, exactly one colliding index
    /// exists per counter either way — the constant *relocates* the collision
    /// to index `c − φ` ≈ 1.1e19, which no stream will reach, rather than
    /// eliminating it. Eliminating it would need the colliding index to depend
    /// on the key; swapping the constant for another would not achieve that.
    #[inline]
    #[must_use]
    pub const fn fork(&self, index: u64) -> Self {
        Self {
            key: mix(self.key ^ mix(index.wrapping_add(0x9e37_79b9_7f4a_7c15))),
            counter: 0,
        }
    }

    /// The next raw 64 bits.
    ///
    /// `const` is not for the sake of const contexts — a stream is drawn from
    /// at runtime. It is a cheap guardrail: a `const fn` cannot reach a clock,
    /// an allocator, or a thread-local, which is the exact list §13.1 bans
    /// from a result-affecting path.
    ///
    /// The counter wraps rather than saturating. Reaching 2^64 draws on one
    /// stream would take longer than the simulation will ever run, and a
    /// wrap repeats the sequence rather than panicking mid-run — which is the
    /// right failure for something that cannot happen.
    #[inline]
    pub const fn next_u64(&mut self) -> u64 {
        let v = mix(self.key ^ mix(self.counter));
        self.counter = self.counter.wrapping_add(1);
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
    /// The two draws are bound to locals rather than passed inline. Rust does
    /// specify left-to-right argument evaluation, so the inline form would be
    /// correct; but "this is deterministic because of an evaluation-order
    /// rule" is a claim a future reader has to go and check, and two `let`s
    /// cost nothing.
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

    /// `mix(0) == 0`, so without the odd constant inside [`Stream::fork`] the
    /// key of `fork(0)` is `mix(key)` — which is exactly the parent's own
    /// first output. This asserts the separation rather than describing it,
    /// because a comment saying "the constant is load-bearing" survives the
    /// constant being deleted and this does not.
    ///
    /// **It has to compare the child's key against the parent's output, not
    /// draw against draw.** The first draft compared first draws and passed
    /// with the constant removed: the child's first output is `mix(key')`, one
    /// further round than the parent's `mix(key)`, so the two never collide
    /// even when the coupling is fully present. The collision is one level up,
    /// and a test aimed one level below it proves nothing.
    #[test]
    fn fork_zero_is_not_the_parents_first_draw() {
        let base = Stream::new(1, Domain::Fold, 0);
        let child = base.fork(0);
        let mut parent = base;
        assert_ne!(
            child.key,
            parent.next_u64(),
            "fork(0)'s key is bit-identical to the parent's first output — the odd \
             constant in `fork` is what separates them"
        );
    }

    /// Ranges chosen to *break* the affine form, not to pass. `(-3.5, 2.25)`
    /// — the only range the first draft tested — is provably immune on all
    /// 2⁵³ grid points, so it would have passed against an implementation
    /// that returned `hi` half the time. The narrow-at-large-offset rows are
    /// the discriminator; the wide row is kept only to show it still works.
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
    /// checked at the values, through the real function. `u` is monotone in
    /// the result, so the two extreme draws bracket all 2⁵³.
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
        13_838_359_919_419_983_661,
        4_598_440_133_459_397_612,
        4_642_742_934_248_969_050,
        13_818_779_806_568_588_147,
        4_600_412_696_514_022_319,
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
        let got: Vec<u64> = (0..4).map(|_| s.next_u64()).collect();
        assert_eq!(got, GOLDEN_UNIVERSE_0);
    }

    const GOLDEN_UNIVERSE_0: [u64; 4] = [
        3_746_585_686_858_627_171,
        4_336_712_865_889_401_412,
        2_674_033_460_618_070_975,
        13_246_354_618_560_148_965,
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
                    9_457_681_971_551_545_159,
                    12_936_856_684_392_112_231,
                    6_926_625_805_525_389_621,
                ],
                Domain::Naming => [
                    6_644_770_736_882_075_221,
                    5_435_968_933_138_178_970,
                    12_450_347_572_679_043_309,
                ],
                Domain::Molecule => [
                    11_303_817_658_843_748_277,
                    14_443_300_859_520_213_470,
                    8_315_562_969_497_253_773,
                ],
                Domain::Fold => [
                    3_323_771_546_716_290_601,
                    2_332_186_771_234_600_492,
                    12_978_889_248_107_649_555,
                ],
                Domain::Reaction => [
                    16_203_379_500_996_458_828,
                    5_972_989_317_398_903_782,
                    7_171_348_452_277_123_823,
                ],
                Domain::Decay => [
                    6_069_625_794_531_361_978,
                    17_342_767_524_156_585_604,
                    4_292_406_413_987_601_818,
                ],
                Domain::Beaker => [
                    108_042_624_964_054_509,
                    4_882_957_866_487_493_608,
                    14_455_476_772_791_300_264,
                ],
                Domain::Shadow => [
                    12_942_870_876_989_012_825,
                    5_956_052_040_092_985_695,
                    10_235_981_688_369_823_609,
                ],
                Domain::Hash => [
                    5_891_565_877_906_749_992,
                    1_669_817_140_753_560_370,
                    17_555_392_361_777_784_652,
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
