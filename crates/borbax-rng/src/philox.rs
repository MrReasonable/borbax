//! Philox 4x64-10, the counter-based bijection underneath every draw (§13.1).
//!
//! Salmon, Moraes, Dror and Shaw, *"Parallel Random Numbers: As Easy as 1, 2,
//! 3"*, SC'11, doi:10.1145/2063384.2063405. This is a transcription of the
//! reference `philox.h` from D. E. Shaw Research's Random123, not an
//! adaptation — the round function, the constants, the key schedule and the
//! round count are all as published.
//!
//! **Why a published algorithm rather than the hand-rolled mixer this crate
//! started with.** Two properties, and only the second needed a paper:
//!
//! 1. It is a **bijection of the counter for any fixed key**. Both multipliers
//!    are odd, so each round is invertible, so the ten-round composition is.
//!    That means distinct counters give distinct output — *by construction*,
//!    not below a birthday bound.
//! 2. That turns the stream-collision problem into a packing problem. The
//!    earlier design hashed `(seed, domain, index)` — 192 bits of coordinate —
//!    into a 64-bit key, which admits permanent collisions the seed cannot
//!    escape: five relations exist below index 2³⁰, two of them pairing
//!    `Decay` with `Shadow`. With a 256-bit counter every coordinate gets its
//!    own 64-bit field and nothing is compressed, so there is no key left to
//!    collide.
//!
//! **The 4x64 variant, not 4x32, and the reason is not speed.** 4x32 carries a
//! 64-bit key and a 128-bit counter — exactly full once domain, index and
//! block are placed, with nothing spare and no room for §13.1's
//! `(universe_seed, world_seed, config_hash)` tuple. 4x64 carries a 128-bit
//! key and a 256-bit counter. It also emits four `u64` per block rather than
//! two, which happens to halve the per-draw cost.
//!
//! **Portability.** Every operation here is `u64` wrapping arithmetic, `XOR`,
//! and a full 64×64→128 multiply via `u128`. All are exactly specified; there
//! is no float, no dispatch, no intrinsic and no `unsafe`. The one thing worth
//! naming is that `u128` multiplication lowers to `umulh`+`mul` on aarch64 and
//! `mulx` on x86-64 — different instructions computing the same exactly
//! defined product.
//!
//! **What makes this trustworthy is the test, not the citation.** The three
//! published known-answer vectors are asserted in this file. A crate reviewed
//! for this project shipped a `Threefry` that failed all three of its
//! reference vectors across nineteen releases, undetected, because its test
//! suite checked only self-consistency and statistics — which cannot tell a
//! correct implementation from a wrong one. That is the failure this module's
//! test exists to prevent, and it is why the vectors live here rather than in
//! a dependency's CI.

/// Philox's first 64-bit multiplier (`PHILOX_M4x64_0`).
const M0: u64 = 0xD2E7_470E_E14C_6C93;
/// Philox's second 64-bit multiplier (`PHILOX_M4x64_1`).
const M1: u64 = 0xCA5A_8263_9512_1157;
/// Key-schedule increment for word 0 — the golden ratio (`PHILOX_W64_0`).
const W0: u64 = 0x9E37_79B9_7F4A_7C15;
/// Key-schedule increment for word 1 — √3 − 1 (`PHILOX_W64_1`).
const W1: u64 = 0xBB67_AE85_84CA_A73B;

/// The published round count. Salmon et al. show 7 rounds already pass
/// `BigCrush`; 10 is the standard margin and is what the reference vectors
/// are computed at. Changing it changes every number this project will ever
/// produce.
const ROUNDS: usize = 10;

/// Full 64×64 → 128 multiply, returned as `(high, low)`.
///
/// Exactly specified, and identical on every target: `u128` multiplication is
/// integer arithmetic, not a widening float op.
#[inline]
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    reason = "both halves of an exact 128-bit product; the truncation to the low 64 bits is the \
              operation, and the shift supplies the high 64"
)]
const fn mulhilo(a: u64, b: u64) -> (u64, u64) {
    let p = (a as u128) * (b as u128);
    ((p >> 64) as u64, p as u64)
}

/// One Philox 4x64 round.
///
/// Note the crossed assignment — `hi1` pairs with `ctr[1]` and `hi0` with
/// `ctr[3]`. It is not a typo and it is load-bearing: it is what makes the
/// round a Feistel-like permutation rather than four independent lanes. The
/// reviewed crate that failed its reference vectors got a neighbouring detail
/// wrong by one term.
#[inline]
#[must_use]
const fn round(ctr: [u64; 4], key: [u64; 2]) -> [u64; 4] {
    let (hi0, lo0) = mulhilo(M0, ctr[0]);
    let (hi1, lo1) = mulhilo(M1, ctr[2]);
    [hi1 ^ ctr[1] ^ key[0], lo1, hi0 ^ ctr[3] ^ key[1], lo0]
}

/// Ten rounds of Philox 4x64 — four `u64` of output from a 256-bit counter
/// and a 128-bit key.
///
/// Stateless and pure: the whole point is that block `n` is computable
/// directly, without having produced blocks `0..n`.
#[must_use]
#[expect(
    clippy::redundant_pub_crate,
    reason = "clippy::redundant_pub_crate and rust's unreachable_pub directly contradict here — \
              the first wants `pub` because the module is private, the second rejects `pub` for \
              exactly that reason. The narrower visibility is the honest one, so the nursery \
              lint yields"
)]
pub(super) const fn philox4x64_10(ctr: [u64; 4], key: [u64; 2]) -> [u64; 4] {
    let mut c = ctr;
    let mut k = key;
    let mut i = 0;
    while i < ROUNDS {
        // The key is bumped *between* rounds, not before the first — the
        // reference applies round 0 with the caller's key unchanged.
        if i > 0 {
            k = [k[0].wrapping_add(W0), k[1].wrapping_add(W1)];
        }
        c = round(c, k);
        i += 1;
    }
    c
}

#[cfg(test)]
mod tests {
    use super::philox4x64_10;

    /// The three published `philox4x64 10` vectors from Random123's
    /// `tests/kat_vectors`, in the file's own order: all-zero, all-ones, and
    /// the digits-of-π vector.
    ///
    /// **This is the test the whole module exists for.** It is the difference
    /// between "we implemented Philox" and "we implemented something with
    /// Philox's statistical shape" — a distinction no amount of chi-squared
    /// testing can make, as the crate surveyed for this decision demonstrated
    /// by shipping a non-conforming `Threefry` through nineteen releases.
    ///
    /// A failure here does not mean the test is stale. It means the generator
    /// is no longer Philox.
    #[test]
    fn matches_the_published_known_answer_vectors() {
        const KAT: &[([u64; 4], [u64; 2], [u64; 4])] = &[
            (
                [0, 0, 0, 0],
                [0, 0],
                [
                    0x1655_4d9e_ca36_314c,
                    0xdb20_fe9d_672d_0fdc,
                    0xd7e7_72ce_e186_176b,
                    0x7e68_b68a_ec7b_a23b,
                ],
            ),
            (
                [u64::MAX, u64::MAX, u64::MAX, u64::MAX],
                [u64::MAX, u64::MAX],
                [
                    0x87b0_92c3_013f_e90b,
                    0x438c_3c67_be8d_0224,
                    0x9cc7_d7c6_9cd7_77b6,
                    0xa09c_aebf_594f_0ba0,
                ],
            ),
            (
                [
                    0x243f_6a88_85a3_08d3,
                    0x1319_8a2e_0370_7344,
                    0xa409_3822_299f_31d0,
                    0x082e_fa98_ec4e_6c89,
                ],
                [0x4528_21e6_38d0_1377, 0xbe54_66cf_34e9_0c6c],
                [
                    0xa528_f454_03e6_1d95,
                    0x38c7_2dbd_566e_9788,
                    0xa5a1_610e_72fd_18b5,
                    0x57bd_43b5_e52b_7fe6,
                ],
            ),
        ];
        for (ctr, key, want) in KAT {
            assert_eq!(
                &philox4x64_10(*ctr, *key),
                want,
                "philox4x64-10 no longer matches its reference vector for ctr {ctr:016x?}"
            );
        }
    }

    /// The property the whole design rests on: distinct counters give distinct
    /// output under a fixed key. Asserted directly over the coordinate layout
    /// `Stream` actually uses, because "provably distinct" is the claim that
    /// replaced a birthday bound and it should not rest on the word "provably".
    #[test]
    fn distinct_counters_give_distinct_blocks() {
        let mut seen = std::collections::BTreeSet::new();
        for domain in 1..=9_u64 {
            for index in 0..2000_u64 {
                seen.insert(philox4x64_10([domain, index, 0, 0], [0x5EED, 0]));
            }
        }
        assert_eq!(seen.len(), 9 * 2000, "counters collided");
    }

    /// Philox is a bijection of the counter, so it cannot map two counters to
    /// one block — but a transcription slip that broke invertibility would
    /// most likely show up as a *stuck* lane rather than as a collision. An
    /// all-zero counter with an all-zero key is the input most likely to
    /// expose one, and the reference vector above already pins it to a value
    /// with no zero word.
    #[test]
    fn the_all_zero_input_does_not_produce_a_zero_block() {
        assert!(philox4x64_10([0; 4], [0; 2]).iter().all(|&w| w != 0));
    }
}
