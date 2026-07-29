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
//! 1. It is a **bijection of the counter for any fixed key**, so distinct
//!    counters give distinct output *by construction* rather than below a
//!    birthday bound. Odd multipliers are necessary for that but **not
//!    sufficient**, and the sufficient part is the one worth writing down: a
//!    round emits `lo`, from which the input word is recoverable because
//!    `mullo(·, M)` is a bijection for odd `M`; the non-injective `hi` is only
//!    ever XOR-mixed into a lane recovered from a *different* low half. That
//!    Feistel pairing is what makes the round invert. Measured on a toy width
//!    with the same odd multiplier: emitting `hi` instead of `lo` is *not*
//!    bijective. "Any odd multiplier, any arrangement" is false, and the
//!    earlier wording licensed it.
//!
//!    Bijectivity is also **round-count independent** — one round is a
//!    bijection too — so it buys non-collision and nothing else. What buys
//!    *independence* between streams differing only in a low-entropy counter
//!    field is diffusion, which is bought entirely by rounds. Measured at one
//!    round: `Decay` and `Shadow` emit literally identical values in two of
//!    four words while distinctness still reads 40000/40000. Salmon et al.,
//!    footnote 11: "four-round Philox fails the avalanche criterion and is
//!    inadequate as a PRNG."
//! 2. That turns the stream-collision problem into a packing problem. The
//!    earlier design hashed `(seed, domain, index)` — 192 bits of coordinate —
//!    into a 64-bit key, which admits permanent collisions the seed cannot
//!    escape: five relations exist below index 2³⁰, against 1.125 expected —
//!    a local fluctuation rather than structure, as `borbax_rng`'s module doc
//!    records with the 2³² follow-up. With a 256-bit counter every coordinate
//!    gets its own 64-bit field and nothing is compressed, so there is no
//!    derived key left to collide.
//!
//! **The 4x64 variant, not 4x32, and the reason is not speed.** 4x32 carries a
//! 64-bit key and a 128-bit counter: the key holds one seed, so the world seed
//! would eat counter words, leaving 64 bits for domain, index, sub and block
//! together. Genuinely too tight. 4x64 carries a 128-bit key and a 256-bit
//! counter — two seeds and four coordinates, each in its own word. (An earlier
//! version of this paragraph counted §13.1's `config_hash` as a third key
//! element; it is not a stream coordinate. The conclusion is unaffected — two
//! seeds alone do not fit 4x32 — but the premise was overstated by one.)
//!
//! It also emits four `u64` per block rather than two, which halves the number
//! of Philox *invocations* per draw — not the
//! time: the paper's Table 2 puts 4x64-10 at 3.2 cpB against 2x64-10's 4.3,
//! which is 26% cheaper per byte, not 50%.
//!
//! **We inherit Philox; we do not inherit the paper's partitioning
//! validation.** §2.2.1 is explicit that stream partitioning needs its own
//! checking, and the partitions they sampled are affine in a single scalar
//! counter — where ours gives each coordinate its own 64-bit word, so the
//! per-draw coordinate sits in `ctr[3]` and enters a multiply only from round
//! 2 (9 rounds of diffusion against the 7 shown sufficient). Ample margin, and
//! measured: `PractRand` finds no anomaly at 1–2 GB on the raw stream, on nine
//! interleaved domains, on 64 interleaved sub-streams, on the top, middle and bottom
//! 32 bits, and on Box-Muller's `(u1, u2)` pairs. That is "no evidence of a
//! problem at this scale", not a replication of `Crush`-resistance.
//!
//! **Portability.** Every operation here is `u64` wrapping arithmetic, `XOR`,
//! and a full 64×64→128 multiply via `u128`. All are exactly specified; there
//! is no float, no dispatch, no intrinsic and no `unsafe`. The one thing worth
//! naming is that `u128` multiplication lowers to `mul`+`umulh` on aarch64 and
//! to the one-operand `mul` on x86-64 — different instructions computing the
//! same exactly specified product. (Not `mulx`: that is BMI2, and neither CI
//! leg enables it. Three reviewers measured `mulq` at our baseline
//! independently, and the digest is identical across the `mulq`/`mulx` split
//! anyway.)
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

/// The reference implementation's default (`PHILOX4x64_DEFAULT_ROUNDS`), and
/// the count the authors themselves recommend: "variants with additional
/// rounds as a safety margin … We favor use of the latter variants".
///
/// Salmon et al. verified *`Crush`-resistance* — `SmallCrush`, Crush and `BigCrush`
/// across their parallel-stream sampling — at **7** rounds for this variant.
/// The `philox.h` wording is weaker and dated ("As of September 2011, the
/// authors know of no statistical flaws with ROUNDS=7 or more"), and the
/// neighbouring 4x32 claim was later walked back to 8. The KAT file publishes
/// vectors at both 7 and 10; we assert the 10-round ones.
///
/// Changing this changes every number this project will ever produce.
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
///
/// **No `#[inline]`, deliberately — and the numbers that were here are
/// withdrawn.** A bare `#[inline]` produces a **byte-identical binary**
/// (sha256 verified, 63 019 lines of disassembly identical), because this
/// function is *already* fully inlined into [`super::Stream::next_u64`] at
/// every optimisation level the workspace uses. The attribute is inert, so the
/// "8.1% worse" once recorded here was measuring nothing.
///
/// `#[inline(always)]` is a real regression but a small one: **+3.35%**, not
/// the 17.1% once claimed, and the stated mechanism was backwards — it
/// produces 198 instructions against 210, so it is *smaller* and slower, which
/// makes the cost scheduling or register pressure rather than code size.
///
/// Keep the attribute off. Just on the +3.35%, and on the fact that the
/// compiler has already made the decision.
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

    /// The property the whole design rests on, tested by **inverting** the
    /// function rather than by counting outputs.
    ///
    /// **The first version of this test counted distinct blocks and asserted
    /// the count — which passes for any function with 256-bit output.**
    /// Measured: with `M0`'s low bit cleared — the single change that
    /// falsifies "both multipliers are odd, so each round is invertible" —
    /// it still reported 18000 distinct of 18000 and passed, while a
    /// hand-constructed collision existed. Its own doc said the claim
    /// "should not rest on the word provably", and it rested on the
    /// birthday bound it was written to replace.
    ///
    /// An explicit inverse cannot pass for a non-bijection.
    #[test]
    #[expect(
        clippy::as_conversions,
        clippy::indexing_slicing,
        reason = "the inverse mirrors the kernel's own 128-bit product and fixed-size [u64; 4] \
                  indexing; both are bounded by construction, as in `mulhilo` and `round`"
    )]
    fn the_kernel_is_invertible() {
        /// `m^-1 mod 2^64` by Newton iteration; exists iff `m` is odd.
        fn inverse_of(m: u64) -> u64 {
            let mut inv = m;
            for _ in 0..6 {
                inv = inv.wrapping_mul(2u64.wrapping_sub(m.wrapping_mul(inv)));
            }
            inv
        }
        let (i0, i1) = (inverse_of(super::M0), inverse_of(super::M1));
        assert_eq!(
            super::M0.wrapping_mul(i0),
            1,
            "M0 is not odd — no inverse exists"
        );
        assert_eq!(
            super::M1.wrapping_mul(i1),
            1,
            "M1 is not odd — no inverse exists"
        );

        // Invert one round: recover ctr[0] from out[3] and ctr[2] from out[1],
        // then the high halves, then un-XOR the remaining two lanes.
        let inv_round = |out: [u64; 4], key: [u64; 2]| -> [u64; 4] {
            let c0 = out[3].wrapping_mul(i0);
            let c2 = out[1].wrapping_mul(i1);
            let hi0 = ((u128::from(super::M0) * u128::from(c0)) >> 64) as u64;
            let hi1 = ((u128::from(super::M1) * u128::from(c2)) >> 64) as u64;
            [c0, out[0] ^ hi1 ^ key[0], c2, out[2] ^ hi0 ^ key[1]]
        };

        // A cheap in-test stirrer; its statistical quality is irrelevant, but
        // its *independence across words* is not. The first draft derived all
        // six input words from one `s` by rotation and XOR, so 20 000 cases
        // explored a one-parameter family of a 384-bit input space — an
        // invertibility test that never varies the words independently is a
        // much weaker test than its case count suggests.
        let mut s = 0x5EED_u64;
        let mut stir = || {
            s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            s ^ (s >> 31)
        };
        for _ in 0..20_000 {
            let ctr = [stir(), stir(), stir(), stir()];
            let key = [stir(), stir()];
            let out = philox4x64_10(ctr, key);
            // Walk the key schedule forward, then invert the rounds backward.
            let mut ks = [key; 10];
            for i in 1..10 {
                ks[i] = [
                    ks[i - 1][0].wrapping_add(super::W0),
                    ks[i - 1][1].wrapping_add(super::W1),
                ];
            }
            let mut c = out;
            for i in (0..10).rev() {
                c = inv_round(c, ks[i]);
            }
            assert_eq!(c, ctr, "philox4x64-10 failed to invert");
        }
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
