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

use borbax_units::det_math;

/// Which subsystem a stream belongs to.
///
/// Streams in different domains are independent, so adding a new random draw
/// in one subsystem cannot shift the sequence any other subsystem sees —
/// which means a change to folding does not silently perturb reaction
/// outcomes.
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
    /// Decay channel selection (§9.4).
    Decay = 6,
    /// The beaker's own scheduling — Gillespie draws and mixing.
    Beaker = 7,
    /// The neutral shadow run, which must not consume the live run's draws
    /// (§2.7).
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

/// `SplitMix64`'s finalizer. Strong avalanche, no state, trivially verifiable.
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
/// **Total by construction, and that is the point.** `ln(0)` is `-inf`, so
/// `u1` is clamped up to the smallest positive normal before it reaches `ln`.
/// `ln(f64::MIN_POSITIVE)` is about `-708.4`, so the radius is at most `~37.7`
/// and the cosine is bounded by 1 — the result is finite for every input.
/// That is what makes `Mass::from_f64_quantised`'s `Err` arm unreachable in
/// Task 4's generator, where `mass = base + next_normal() * spread` is the
/// only realistic route to a bad value.
///
/// **Clamping, not rejection, and the choice is not free.** Redrawing until
/// `u1 > 0` would also make this total, but it would consume further draws
/// and shift every subsequent value in the domain's stream — regenerating the
/// entire universe. Clamping cannot. The draw count is pinned in
/// `next_normal_consumes_exactly_two_draws` from this commit onward.
///
/// The clamp is written as a comparison rather than `f64::max`, which §13.1
/// bans for its non-deterministic `±0.0` behaviour. The comparison form also
/// folds a NaN `u1` to `MIN_POSITIVE` instead of propagating it. That is
/// irrelevant to the only caller — `u1` comes from `next_f64` and cannot be
/// NaN — but it is why this is total for *every* input rather than only for
/// the reachable ones.
fn normal_from(u1: f64, u2: f64) -> f64 {
    let u1 = if u1 > f64::MIN_POSITIVE {
        u1
    } else {
        f64::MIN_POSITIVE
    };
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
    #[must_use]
    pub const fn new(seed: u64, domain: Domain, index: u64) -> Self {
        let key =
            mix(seed ^ mix(domain.discriminant().wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ mix(index)));
        Self { key, counter: 0 }
    }

    /// Derive a child stream. Used where a fixed number of sub-streams is
    /// needed — one per molecule, one per fold attempt — without threading
    /// a mutable parent through the call graph.
    ///
    /// The odd constant is domain separation, not decoration: `mix(0) == 0`,
    /// so without it `fork(0)`'s key would be bit-identical to the parent's
    /// first output. Harmless today, and exactly the kind of structural
    /// coupling that becomes a real correlation once someone forks in a loop.
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
    pub const fn next_u64(&mut self) -> u64 {
        let v = mix(self.key ^ mix(self.counter));
        self.counter = self.counter.wrapping_add(1);
        v
    }

    /// Uniform in `[0, 1)`. Uses the top 53 bits, which is exactly the
    /// mantissa width of `f64` — so every representable value in the
    /// interval is reachable and none is favoured.
    ///
    /// Both steps are exact: the shifted integer is below 2^53 so it converts
    /// without rounding, and `SCALE` is a power of two so the multiply cannot
    /// round either. The result is `n / 2^53` for an integer `n`, exactly.
    #[expect(
        clippy::as_conversions,
        reason = "std offers no lossless u64 -> f64 conversion; the shift above is what makes \
                  this one lossless, and the exactness is asserted in floats_are_in_unit_interval"
    )]
    #[expect(
        clippy::cast_precision_loss,
        reason = "the shift leaves at most 53 significant bits, which is exactly f64's mantissa \
                  width, so no precision is lost for any input"
    )]
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

    /// Uniform in `[lo, hi)`.
    pub fn next_f64_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.next_f64() * (hi - lo)
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

    #[test]
    fn floats_in_a_range_are_bounded() {
        let mut s = Stream::new(5, Domain::Beaker, 0);
        for _ in 0..10_000 {
            let v = s.next_f64_range(-3.5, 2.25);
            assert!((-3.5..2.25).contains(&v), "out of range: {v}");
        }
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
    /// `Domain`-separation argument rests on. Discriminants that are `1..=N`
    /// in declaration order are what make an append provably harmless, so the
    /// shape is asserted rather than left to the next person's care.
    #[test]
    fn discriminants_are_contiguous_from_one_in_declaration_order() {
        for (i, domain) in Domain::ALL.iter().enumerate() {
            let want = u64::try_from(i).unwrap_or(u64::MAX) + 1;
            assert_eq!(
                domain.discriminant(),
                want,
                "{domain:?} is out of order or its discriminant was renumbered"
            );
        }
    }
}
