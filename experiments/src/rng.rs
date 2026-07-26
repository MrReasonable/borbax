//! A deterministic stream, standing in for `borbax_rng::Stream`.
//!
//! `SplitMix64`: a fixed odd increment added to the state, then a strong
//! finalising mix. The mix is what makes seeds 0, 1, 2... behave as
//! independent streams rather than correlated ones, which matters here
//! because the harness derives one stream per trial index.
//!
//! Only shifts, xors, wrapping adds and wrapping multiplies — all exactly
//! specified, so the sequence is identical on every platform (§13.4).
//!
//! Deliberately not `Copy`, matching the doctrine for the real `Stream`: a
//! copy would silently fork the sequence and replay the same draws, and the
//! resulting duplicate "random" molecules would look like a chemistry result.

/// A deterministic pseudo-random stream.
#[derive(Debug, Clone)]
pub struct Stream {
    state: u64,
}

/// 2^64 divided by the golden ratio — odd, so adding it repeatedly visits
/// every 64-bit state before repeating.
const GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;

impl Stream {
    /// A stream for `seed`. Any seed is valid, including zero.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next 64 bits.
    pub const fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(GAMMA);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
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
    /// Built from the top 53 bits — the exact width of an f64 mantissa — so
    /// every representable value in the interval is reachable and the
    /// conversion is exact rather than rounded.
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "the shift leaves 53 bits, so both conversions are exact — the precision loss \
                  the lint warns about is what the shift already removed"
    )]
    pub const fn f64_unit(&mut self) -> f64 {
        // 2^-53, written as a division by a power of two so it is exact.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
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
