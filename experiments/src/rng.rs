//! A deterministic stream, standing in for `borbax_rng::Stream`.
//!
//! `SplitMix64`: a fixed odd increment added to the state, then a strong
//! finalising mix. Steele, Lea & Flood, *"Fast splittable pseudorandom number
//! generators"*, OOPSLA 2014, doi:10.1145/2660193.2660195 — the algorithm behind
//! Java's `SplittableRandom`. The finaliser's constants are `MurmurHash3`'s,
//! found by search rather than derived. The mix is what makes seeds 0, 1, 2... behave as
//! independent streams rather than correlated ones, which matters here
//! because the harness derives one stream per trial index.
//!
//! Only shifts, xors, wrapping adds and wrapping multiplies — all exactly
//! specified, so the sequence is identical on every platform (§13.4).
//!
//! Deliberately not `Copy`, matching the doctrine for the real `Stream`: a
//! copy would silently fork the sequence and replay the same draws, and the
//! resulting duplicate "random" molecules would look like a chemistry result.

/// `SplitMix64`'s finalising mix, on its own.
///
/// `pub(crate)`, not `pub`: the only caller outside this file is
/// `molecule.rs`'s `form_hash`, and this module's own header says it stands in
/// for `borbax_rng::Stream` — which exports `Stream` and `Domain` and no bare
/// mixer. A `pub` mixer makes the stand-in diverge from the thing it stands in
/// for, in the direction that sends a reader looking for `borbax_rng::mix64`.
///
/// It exists because `molecule.rs` had grown a second, byte-identical copy of
/// these five lines, and two copies of a constant drift when one is edited and
/// the other is not — which is exactly how `bonds.rs` came to hold a `MAX` and
/// a `MAX_LEN` that disagreed, pricing bond order 7 at order 1's multiplier.
#[must_use]
pub(crate) const fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> MIX_SHIFT_1)).wrapping_mul(MIX_MULTIPLIER_1);
    z = (z ^ (z >> MIX_SHIFT_2)).wrapping_mul(MIX_MULTIPLIER_2);
    z ^ (z >> MIX_SHIFT_3)
}

/// A deterministic pseudo-random stream.
#[derive(Debug, Clone)]
pub struct Stream {
    state: u64,
}

/// The step added to the state on every draw: 2^64 divided by the golden ratio.
///
/// **Odd**, which is the whole requirement — an odd increment added repeatedly
/// to a 64-bit register visits all 2^64 states before repeating, so the period
/// is maximal for any seed. The golden-ratio value in particular spreads
/// successive states as evenly as possible around the register rather than
/// clustering them (a Weyl sequence).
const GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;

/// First multiplier of the finalising mix.
///
/// **The multipliers and the shifts must not be "tidied", but for two different
/// reasons — and an earlier version of this comment gave one reason for both,
/// wrongly.**
///
/// *The multipliers* are load-bearing for avalanche, and the claim was measured:
/// substituting `0x1000…0001`, `0xFFFF…FFFF` or small primes each collapses the
/// worst pairwise Hamming distance in
/// `adjacent_seeds_decorrelate_immediately` from 18 bits to 1 or 2, and that
/// test fails.
///
/// *The shifts* are **not** guarded by any test in this file, and are not
/// load-bearing for avalanche either. Measured: `(32,32,32)`, `(30,27,32)`,
/// `(30,27,1)` and `(1,1,1)` all leave the whole suite green, and a strict
/// avalanche measurement puts the published `(30,27,31)` at 0.0283 against
/// `(32,32,32)`'s 0.0293 — indistinguishable. What makes them load-bearing is
/// **reproducibility**: every figure in `docs/experiments/` was produced by this
/// exact sequence, so changing 30 to 32 silently re-bases published results
/// rather than degrading anything measurable.
///
/// They are the published MurmurHash3-derived finaliser `SplitMix64` uses.
const MIX_MULTIPLIER_1: u64 = 0xBF58_476D_1CE4_E5B9;

/// Second multiplier of the finalising mix. See [`MIX_MULTIPLIER_1`].
const MIX_MULTIPLIER_2: u64 = 0x94D0_49BB_1331_11EB;

/// Right-shift before the first multiply. See [`MIX_MULTIPLIER_1`].
const MIX_SHIFT_1: u32 = 30;

/// Right-shift before the second multiply. See [`MIX_MULTIPLIER_1`].
const MIX_SHIFT_2: u32 = 27;

/// Final right-shift, folding the high bits down over the low ones.
const MIX_SHIFT_3: u32 = 31;

/// Bits in an `f64` significand, counting the implicit leading one.
///
/// A `u64` shifted right by `64 - MANTISSA_BITS` leaves exactly this many
/// significant bits, so the integer-to-float conversion below is **exact** —
/// no value is rounded, and every representable float in `[0, 1)` at this
/// spacing is reachable.
const MANTISSA_BITS: u32 = 53;

/// How far to shift a 64-bit draw down to leave [`MANTISSA_BITS`] bits.
const UNIT_SHIFT: u32 = u64::BITS - MANTISSA_BITS;

impl Stream {
    /// A stream for `seed`. Any seed is valid, including zero.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next 64 bits.
    pub const fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(GAMMA);
        mix64(self.state)
    }

    /// A uniform value in `0..n`.
    ///
    /// Lemire's multiply-shift with rejection: exactly uniform, and the
    /// rejection branch is taken with probability under `n / 2^64`, so in
    /// practice never. Modulo would be biased toward small values, which for
    /// a harness that picks atoms and edit sites would quietly skew the
    /// molecules it generates.
    ///
    /// # Panics
    ///
    /// Never for `n > 0`. `n == 0` has no valid result and panics in debug
    /// via the assertion below.
    #[allow(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        reason = "the two truncations are the algorithm: Lemire splits the 128-bit product into \
                  its high half (the result) and its low half (the rejection test). `try_from` \
                  would be asserting the discarded half is zero, which is exactly what it is not"
    )]
    pub const fn below(&mut self, n: u64) -> u64 {
        debug_assert!(n > 0, "below(0) has no valid result");
        let threshold = n.wrapping_neg() % n;
        loop {
            let x = self.next_u64();
            let m = (x as u128) * (n as u128);
            if (m as u64) >= threshold {
                return (m >> 64) as u64;
            }
        }
    }

    /// A uniform index into a collection of length `n`, or `0` when empty.
    ///
    /// Exists so that no caller has to write its own `usize`/`u64` cast. Every
    /// conversion here is a `try_from`, so the whole harness stays free of
    /// `as` on integers — which is the point of the workspace's
    /// `as_conversions` lint, and worth one small function to keep.
    ///
    /// Empty collections return `0` rather than panicking. There is no valid
    /// index into nothing, and the callers that can hit it (a molecule with no
    /// leaves, a zero-atom edit) already guard the result.
    #[must_use]
    pub fn index(&mut self, n: usize) -> usize {
        let Ok(bound) = u64::try_from(n) else {
            return 0;
        };
        if bound == 0 {
            return 0;
        }
        usize::try_from(self.below(bound)).unwrap_or(0)
    }

    /// A uniform value in `[0, 1)`.
    ///
    /// Built from the top 53 bits — the exact width of an `f64` significand —
    /// so the conversion is exact rather than rounded: every result is `k/2^53`
    /// for an integer `k`, and every such value is representable.
    ///
    /// **Not "every representable value in the interval is reachable"**, which
    /// this said and which is wrong by a factor of 511.5: there are ~2^62
    /// distinct `f64` in `[0, 1)` because the exponent shrinks toward zero, and
    /// only 2^53 of them are on this grid. The qualifier "at that spacing" is
    /// load-bearing — `MANTISSA_BITS` has it and this did not.
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "the shift leaves 53 bits, so both conversions are exact — the precision loss \
                  the lint warns about is what the shift already removed"
    )]
    pub const fn f64_unit(&mut self) -> f64 {
        // Divide by 2^MANTISSA_BITS. Written as a division by a power of two
        // rather than a multiply by a decimal literal so the scaling is exact.
        (self.next_u64() >> UNIT_SHIFT) as f64 / (1u64 << MANTISSA_BITS) as f64
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds over sized collections; \
              casts are of values the surrounding assertion has just bounded"
)]
mod tests {
    use super::*;

    #[test]
    fn the_same_seed_gives_the_same_sequence() {
        let mut a = Stream::new(0xDEAD_BEEF);
        let mut b = Stream::new(0xDEAD_BEEF);
        let xs: Vec<u64> = (0..64).map(|_| a.next_u64()).collect();
        let ys: Vec<u64> = (0..64).map(|_| b.next_u64()).collect();
        assert_eq!(xs, ys);
    }

    #[test]
    fn different_seeds_give_different_sequences() {
        let mut a = Stream::new(1);
        let mut b = Stream::new(2);
        let xs: Vec<u64> = (0..64).map(|_| a.next_u64()).collect();
        let ys: Vec<u64> = (0..64).map(|_| b.next_u64()).collect();
        assert_ne!(xs, ys);
    }

    /// Adjacent seeds are the case a counter-based stream has to handle: the
    /// harness derives one stream per trial index, so seeds 0, 1, 2... must not
    /// produce correlated first draws. A raw LCG fails this badly.
    #[test]
    fn adjacent_seeds_decorrelate_immediately() {
        let firsts: Vec<u64> = (0..64u64).map(|s| Stream::new(s).next_u64()).collect();
        let distinct: std::collections::BTreeSet<u64> = firsts.iter().copied().collect();
        assert_eq!(
            distinct.len(),
            firsts.len(),
            "first draws collide across seeds"
        );
        // Every pair should differ in roughly half its bits. Anything under 16
        // means the seeds are leaking straight through into the output.
        for i in 0..firsts.len() {
            for j in (i + 1)..firsts.len() {
                let h = (firsts[i] ^ firsts[j]).count_ones();
                assert!(h >= 16, "seeds {i} and {j} differ in only {h} bits");
            }
        }
    }

    #[test]
    fn below_stays_in_range_and_covers_it() {
        let mut r = Stream::new(7);
        let mut seen = [false; 5];
        for _ in 0..1000 {
            let v = r.below(5);
            assert!(v < 5, "below(5) returned {v}");
            seen[v as usize] = true;
        }
        assert!(
            seen.iter().all(|&s| s),
            "below(5) never returned some value"
        );
    }

    /// `index` exists so that callers picking an atom or an edit site never
    /// write a cast of their own. It has to be total: a zero-length
    /// collection has no valid index, and returning one anyway would be an
    /// out-of-bounds panic at the call site instead of here.
    #[test]
    fn index_picks_within_a_collection_and_tolerates_an_empty_one() {
        let mut r = Stream::new(13);
        let mut seen = [false; 7];
        for _ in 0..1000 {
            let i = r.index(7);
            assert!(i < 7, "index(7) returned {i}");
            seen[i] = true;
        }
        assert!(
            seen.iter().all(|&s| s),
            "index(7) never returned some value"
        );
        assert_eq!(r.index(0), 0, "index(0) must not panic");
        assert_eq!(r.index(1), 0);
    }

    #[test]
    fn below_one_is_always_zero() {
        let mut r = Stream::new(11);
        for _ in 0..32 {
            assert_eq!(r.below(1), 0);
        }
    }

    #[test]
    fn unit_floats_are_in_the_half_open_unit_interval() {
        let mut r = Stream::new(99);
        for _ in 0..10_000 {
            let f = r.f64_unit();
            assert!((0.0..1.0).contains(&f), "out of range: {f}");
        }
    }
}
