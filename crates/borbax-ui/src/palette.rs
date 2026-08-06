//! What colour an atom is drawn, as a function of what the generator made it.
//!
//! **This file imports nothing.** Not the engine, not the UI toolkit, not even
//! a `borbax-*` crate — it takes bare scalars and returns a bare triple. That is
//! deliberate twice over: the seam guard has nothing to say about a file with
//! no imports, and every test here runs as plain logic with no window and no
//! universe.
//!
//! # Why it still takes bare scalars, now that G5 is gone
//!
//! **§5's G5 was withdrawn on 2026-08-06 and the `xtask` guard that enforced
//! this file went with it.** That guard pinned [`colour`]'s signature and
//! banned lookups in it, on the grounds that a per-element colour table using
//! conventional colours is the "mapping table between Borbax entities and
//! real-world entities" G5 forbade. G5 no longer exists; see the spec §5 for
//! why, and issue #26 for the direction that made it untenable.
//!
//! The signature is unchanged anyway, because the reason for it was never only
//! G5. Taking bare scalars is what keeps this file in the seam's strictest tier
//! — it imports nothing, so every test here runs as plain logic with no window
//! and no universe. That is worth keeping on its own.
//!
//! **What was worth learning from the guard, since it cost five review rounds.**
//! No syntactic check can decide "closed-form formula" from "lookup": every
//! total function `u8 -> f64` is expressible as arithmetic, so the artefact was
//! rebuilt past each fix — in an `impl` block, in a helper, behind a `mut`
//! keyword, as a bit-packed integer, as a fitted polynomial. Each fix also cost
//! a false positive on honest code. If a structural guard is ever wanted here
//! again, that is the thing to weigh first.
//!
//! # The obligation that outlives the guard
//!
//! A guard could catch a colour constant. It could never catch a ramp
//! *anchored* so the picture happens to read like a real scheme, because no
//! forbidden word is written down and the breach lives entirely in the
//! justification. That obligation is unchanged and is now the only one:
//!
//! **Every channel assignment below names a generated quantity and a
//! measurement over it.** The executable form is to re-read each justification
//! with the seed changed — *"hue steps per group because families share a
//! group"* survives, because it is a claim about the generator; *"the heavy
//! atom is dark because heavy things read as dark"* does not.
//!
//! **It is a property, not a banned phrase, and the first version got that
//! wrong.** It forbade the words "looks right" and "reads as", which fired on
//! three legitimate lines in its own commit — all legibility arguments rather
//! than resemblance ones. An obligation that flags correct code in the commit
//! that introduces it is how guards get deleted.
//!
//! # Why each channel is the channel it is
//!
//! Measured over 500 universes on the demo molecule, before any of this was
//! written:
//!
//! - **Hue from `group`, stepped by the golden angle.** Families are elements
//!   sharing a group (§7.1), so a column sharing a hue is the payload §7.1
//!   names by name — *"wait, everything in this column makes rings"*. The
//!   stepping is what makes it legible rather than merely correct: the four
//!   leaves land at **consecutive** groups (1, 2, 3, 4) in every universe
//!   measured, so a smooth hue ramp gives them four near-identical colours.
//!   Measured, minimum pairwise separation across the five atoms: **0.1932**
//!   under this, and **0.6748** among the leaves alone. A hue ramp keyed on
//!   affinity instead collapses them to under **0.05** — three reviewers
//!   reproduced the two figures above to four digits and each got a different
//!   number for the ramp (0.041, 0.045, 0.046), because the counterfactual's
//!   construction was never written down. The order of magnitude is the claim;
//!   the digits are not, and the earlier "0.0239" implied a precision that
//!   nothing here can reproduce.
//! - **Saturation from `affinity`.** The one channel that adds information the
//!   picture does not already carry: radius order equals mass order in 500/500
//!   universes, so a mass-driven hue would merely restate the sphere sizes,
//!   while radius order equals affinity order in **0** of 500.
//! - **Zero saturation at `valence == 0`.** A closed outer shell has no contact
//!   left to make, so those elements bond with nothing and are drawn grey. This
//!   is read off the generated valence, never from a family name — the
//!   distinction Step 4b's plan step is explicit about.
//! - **Lightness from `mass`,** normalised over the universe's own attained
//!   range. It is the weakest channel deliberately: it duplicates what radius
//!   already shows, and it is here so two atoms sharing a group are still told
//!   apart.
//!
//! # Determinism
//!
//! **No transcendental appears here, and that is a property of the colour space
//! rather than a restraint.** HSL to RGB in this form is piecewise linear:
//! `abs`, `rem_euclid`, comparisons and the four arithmetic operations.
//!
//! **`rem_euclid` is not native, and an earlier version of this sentence said it
//! was.** Float `%` lowers to a **libm `fmod` call** — measured, one `bl _fmod`
//! on aarch64 and two `callq *fmod@GOTPCREL` on x86-64, with `_fmod` an
//! undefined import in both. The conclusion survives for the reason
//! `clippy.toml`'s header insists on: **the criterion is whether IEEE-754
//! specifies the result exactly, not whether the operation calls libm.** `fmod`
//! is exact — `x - trunc(x/y)*y` is representable in the same format, so there
//! is no rounding freedom to differ on. Verified rather than argued: Apple's
//! `fmod` against the vendored `libm::fmod` over 20M random finite inputs gives
//! **0 mismatches**, and `colour` hashed over 16.4M inputs built for both
//! architectures gives the **same hash**.
//!
//! The wrong reason mattered exactly as CLAUDE.md says it does: "`rem_euclid`
//! is native" is the sentence that would wave through the next `%`-adjacent
//! operation that is *not* exactly specified.
//!
//! The golden angle is a literal.
//!
//! A perceptual space (Oklab) would need `cbrt` and would be admissible —
//! `det_math::cbrt` exists — but it is not needed to place five colours apart,
//! and `cbrt` is the one call CLAUDE.md names as the only FMA site in the
//! release binaries. The cheaper space is also the one with less to argue
//! about.

/// The lowest `affinity` any element attains, over the universes surveyed.
///
/// Paired with [`AFFINITY_HI`] to map the channel onto `0..=1`. Both bounds are
/// **clamped against** rather than assumed, so an element outside the surveyed
/// range saturates instead of producing a colour outside the cube.
const AFFINITY_LO: f64 = 0.140_476;

/// The highest `affinity` any element attains.
const AFFINITY_HI: f64 = 1.0;

/// Degrees of hue per group.
///
/// **The golden angle, and the reason is a measurement rather than an
/// aesthetic.** Consecutive integers multiplied by it land maximally far apart
/// on the wheel and never repeat until the wheel is densely covered. The four
/// leaves of the demo molecule sit at consecutive groups in every universe
/// measured, so consecutive-group separation is exactly the property the
/// picture needs and exactly the one a linear ramp does not have.
const GOLDEN_ANGLE: f64 = 137.507_764_050_037_85;

/// The darkest an atom is drawn.
///
/// The band is narrow on purpose: lightness is the weakest channel, and a full
/// `0..=1` sweep would drive the lightest atoms to white and the heaviest to
/// black, in both cases destroying the hue that carries the family.
const LIGHT_LO: f64 = 0.34;

/// The lightest an atom is drawn.
const LIGHT_HI: f64 = 0.72;

/// The least saturated a *bonding* atom is drawn.
///
/// Deliberately well above zero, so a low-affinity atom is never confused with
/// the fully-desaturated grey that means "closed shell, bonds with nothing".
const SAT_LO: f64 = 0.45;

/// The most saturated an atom is drawn.
const SAT_HI: f64 = 0.82;

/// The colour an atom is drawn, as `[r, g, b]` in `0..=1`, sRGB.
///
/// Every argument is a **generated scalar** — no symbol, no name, no id, no
/// position in the table. Nothing enforces that any more (§5's G5 and its guard
/// were withdrawn on 2026-08-06), and it is kept because it is what holds this
/// file in the seam's strictest tier: with no identity type in the signature
/// there is nothing to import, so the whole file tests as plain logic.
///
/// **`group` is not an identity argument, and the difference is the whole
/// design.** It is the generated column index — the quantity that *makes* a
/// family under §7.1 — and many elements share one. An argument that named
/// *which* element this is would be injective; this one is deliberately not.
///
/// `mass_lo` and `mass_hi` are the attained mass range of the universe this
/// element belongs to, so lightness uses the spread that universe actually has
/// rather than a constant chosen from one of them. `mass_hi <= mass_lo` is not
/// an error — a universe with one distinct mass has no spread to show, and the
/// channel sits at its midpoint.
///
/// **One function, not a convenience wrapper beside it.** A second spelling
/// defaulting `group` would be a second surface for the signature guard to pin
/// and a second place for an identity argument to arrive.
#[must_use]
pub fn colour(
    affinity: f64,
    valence: u8,
    mass: f64,
    mass_lo: f64,
    mass_hi: f64,
    group: u8,
) -> [f64; 3] {
    let hue = f64::from(group) * GOLDEN_ANGLE;

    let reach = unit(affinity, AFFINITY_LO, AFFINITY_HI);
    // **Zero, not merely low.** A closed outer shell has no contact left to
    // make, so the atom bonds with nothing; grey says that without naming a
    // family. Reading it off the generated valence is what keeps it a
    // consequence rather than a declaration.
    let saturation = if valence == 0 {
        0.0
    } else {
        SAT_LO + reach * (SAT_HI - SAT_LO)
    };

    let heft = unit(mass, mass_lo, mass_hi);
    let lightness = LIGHT_LO + heft * (LIGHT_HI - LIGHT_LO);

    hsl_to_rgb(hue, saturation, lightness)
}

/// `x` mapped onto `0..=1` across `lo..=hi`, clamped at both ends.
///
/// A degenerate range gives the midpoint rather than a division by zero — the
/// honest answer for "this universe has one distinct mass" is "no spread to
/// show", not `NaN` painted onto every atom.
fn unit(x: f64, lo: f64, hi: f64) -> f64 {
    let span = hi - lo;
    if span <= 0.0 {
        return 0.5;
    }
    let t = (x - lo) / span;
    // **`clamp` propagates NaN, so the clamp is not the guard it looks like.**
    // A NaN input passes straight through it, out of `colour`, and through
    // `scene::srgb` into a NaN `Color` — an atom painted nothing at all, with
    // no error anywhere. The out-of-range test covered `f64::MIN` and
    // `f64::MAX` and not the one out-of-range value that is *not* clamped.
    //
    // The midpoint is the same answer a range with no spread gets: "no
    // information to show", rather than a hole in the picture.
    if t.is_nan() { 0.5 } else { t.clamp(0.0, 1.0) }
}

/// HSL to sRGB, with hue in degrees and both others in `0..=1`.
///
/// **Piecewise linear, which is the whole §13.1 argument.** `abs`, comparisons
/// and the four arithmetic operations are exactly specified by IEEE-754;
/// `rem_euclid` is a `fmod` libcall and is *also* exactly specified, which is
/// the criterion that matters. See the module header — the distinction is not
/// pedantry, it is the difference between a rule and a habit.
fn hsl_to_rgb(hue: f64, saturation: f64, lightness: f64) -> [f64; 3] {
    let chroma = (1.0 - (2.0 * lightness - 1.0).abs()) * saturation;
    // `rem_euclid` rather than `%`, so a hue below zero wraps onto the wheel
    // instead of reflecting through it. No caller passes one today; the
    // difference is invisible until one does.
    let sector = hue.rem_euclid(360.0) / 60.0;
    let ramp = chroma * (1.0 - (sector.rem_euclid(2.0) - 1.0).abs());

    // **`sector` is *not* always below 6.0, and the first version of this
    // comment said it was.** `(-1e-20_f64).rem_euclid(360.0)` returns exactly
    // 360.0 — the `r + 360.0` fixup rounds up for a tiny negative `r` — giving
    // `sector == 6.0` and taking the catch-all arm. The function is still
    // total, and correct, for a different reason than the one claimed: at
    // sector 6 the sixth arm's `ramp` is 0, so it returns `(chroma, 0, 0)`,
    // which is bit-identical to what sector 0 returns. Verified: hue 0, hue 360
    // and hue -1e-20 all give the same triple.
    let (red, green, blue) = if sector < 1.0 {
        (chroma, ramp, 0.0)
    } else if sector < 2.0 {
        (ramp, chroma, 0.0)
    } else if sector < 3.0 {
        (0.0, chroma, ramp)
    } else if sector < 4.0 {
        (0.0, ramp, chroma)
    } else if sector < 5.0 {
        (ramp, 0.0, chroma)
    } else {
        (chroma, 0.0, ramp)
    };

    let floor = lightness - chroma / 2.0;
    [red + floor, green + floor, blue + floor]
}

#[cfg(test)]
mod tests {
    use super::{AFFINITY_HI, AFFINITY_LO, GOLDEN_ANGLE, LIGHT_HI, LIGHT_LO, colour, hsl_to_rgb};

    /// How far apart two colours are, weighted so green counts most.
    ///
    /// **Not a perceptual metric and not claimed to be one.** It is a cheap
    /// monotone stand-in, used only to place a bar between two measured
    /// populations that are an order of magnitude apart — which is a job it can
    /// do and a CIE metric would do no better. Naming it honestly matters
    /// because the next reader will otherwise cite it as perceptual distance.
    fn separation(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dr = a[0] - b[0];
        let dg = a[1] - b[1];
        let db = a[2] - b[2];
        (2.0 * dr * dr + 4.0 * dg * dg + 3.0 * db * db).sqrt()
    }

    /// Every colour is inside the unit cube, over the whole input space.
    ///
    /// Swept rather than sampled: an HSL conversion that leaves the cube does
    /// so at the corners of the space, not in the middle of it.
    #[test]
    fn every_colour_is_inside_the_cube() {
        let mut checked = 0_u32;
        for group in 0..24_u8 {
            for valence in 0..7_u8 {
                for affinity_step in 0..11 {
                    let affinity =
                        AFFINITY_LO + f64::from(affinity_step) / 10.0 * (AFFINITY_HI - AFFINITY_LO);
                    for mass_step in 0..11 {
                        let mass = f64::from(mass_step);
                        let c = colour(affinity, valence, mass, 0.0, 10.0, group);
                        for channel in c {
                            assert!(
                                (0.0..=1.0).contains(&channel),
                                "colour {c:?} leaves the unit cube at group {group}, \
                                 valence {valence}, affinity {affinity}, mass {mass}"
                            );
                        }
                        checked += 1;
                    }
                }
            }
        }
        // The bar is the size of the sweep, so a loop that stopped selecting
        // anything is visible rather than silently green.
        assert_eq!(checked, 24 * 7 * 11 * 11, "the sweep did not run in full");
    }

    /// An input outside the surveyed range saturates rather than escaping.
    ///
    /// The clamps are the reason the sweep above can be believed for a universe
    /// that produces an affinity or a mass nobody surveyed.
    ///
    /// **`f64::NAN` is in the list because `clamp` propagates it**, so it was
    /// the one out-of-range value the clamp did *not* catch — it passed through
    /// `colour` and would have reached the engine as a `Color` of nothing.
    /// Found by a review; the earlier list covered `MIN`/`MAX`, which the clamp
    /// handles, and not the value it does not.
    #[test]
    fn an_out_of_range_input_is_clamped_rather_than_escaping() {
        for (affinity, mass) in [
            (-5.0, -5.0),
            (5.0, 500.0),
            (f64::MIN, f64::MAX),
            (f64::NAN, 5.0),
            (0.5, f64::NAN),
            (f64::NAN, f64::NAN),
        ] {
            let c = colour(affinity, 3, mass, 0.0, 10.0, 2);
            for channel in c {
                assert!(
                    (0.0..=1.0).contains(&channel),
                    "colour {c:?} left the cube for affinity {affinity}, mass {mass}"
                );
            }
        }
    }

    /// The largest `group` any universe attains.
    ///
    /// Measured over 500 universes. The sweeps below run the range the
    /// generator actually produces rather than a round number — an earlier
    /// version swept `0..24`, which is a third of it, and the interesting
    /// collisions are all outside that.
    const MAX_GROUP: u8 = 73;

    /// No two groups are drawn the same colour.
    ///
    /// **This is the test that pins the golden angle, and the one it replaced
    /// did the opposite.** `consecutive_groups_are_far_apart_in_colour`
    /// measured the separation between *neighbouring* groups — a quantity that
    /// is identical for every pair under a constant step, and is therefore
    /// maximised at 180°, which puts groups 0 and 2 on exactly the same colour.
    /// Measured in this metric over `0..24`: the golden angle scores 1.3455,
    /// **180° scores 1.8268 and 120° scores 1.3616** — both better, and both
    /// collapsing distinct families onto identical colours. A guard that ranks
    /// the two worst candidates above the shipped one is not a guard.
    ///
    /// The property that actually distinguishes the golden angle is that it is
    /// the most badly approximable rotation, so `{k·α}` never repeats and the
    /// minimum gap over the first N points stays as large as it can be **for
    /// every N at once** — which matters because the number of groups is drawn
    /// per universe and no constant can be tuned to it.
    ///
    /// Measured minimum over all pairs in `0..=MAX_GROUP`: **0.0420**, at groups
    /// 11 and 66. Every rejected candidate scores exactly **0.0000** somewhere
    /// in the same range, including a linear `360/24` ramp, which is exact at
    /// 24 groups and collides three times over the range a universe reaches.
    /// The bar sits between those two populations.
    ///
    /// **0.0420 is well below the 0.35 legibility bar**, and that is stated
    /// rather than hidden: this test says two families are never *the same*
    /// colour, not that every pair is *tellable apart*. The demo molecule's
    /// atoms are neighbours in group, which is what
    /// `neighbouring_groups_are_far_apart_in_colour` covers; a molecule
    /// spanning groups 11 and 66 would show two colours a person could not
    /// separate, and that belongs to whichever step first builds one.
    #[test]
    fn no_two_groups_are_drawn_the_same_colour() {
        let of = |group: u8| colour(0.6, 3, 5.0, 0.0, 10.0, group);
        let mut worst = f64::INFINITY;
        let mut worst_pair = (0_u8, 0_u8);
        let mut checked = 0_u32;
        for a in 0..=MAX_GROUP {
            for b in (a + 1)..=MAX_GROUP {
                let d = separation(of(a), of(b));
                if d < worst {
                    worst = d;
                    worst_pair = (a, b);
                }
                checked += 1;
            }
        }
        let pairs = u32::from(MAX_GROUP) * (u32::from(MAX_GROUP) + 1) / 2;
        assert_eq!(checked, pairs, "the sweep did not run in full");
        assert!(
            worst > 0.03,
            "groups {} and {} are drawn {worst:.4} apart — a step that repeats \
             puts two families on one colour, and every rejected candidate \
             (180°, 120°, a linear 360/24 ramp) scores exactly 0 somewhere in \
             this range",
            worst_pair.0,
            worst_pair.1
        );
    }

    /// Neighbouring groups are far apart in colour.
    ///
    /// **This does not pin the golden angle and no longer claims to** — see
    /// `no_two_groups_are_drawn_the_same_colour` for why the quantity it
    /// measures ranks 180° above the shipped constant. What it does cover is
    /// the case the demo molecule actually presents: its four leaves sit at
    /// *consecutive* groups in every universe measured, so neighbouring
    /// separation is what decides whether the molecule shows four leaves or one
    /// blur.
    ///
    /// The bar is 0.35 against a measured 1.3428 for the golden angle over the
    /// attained range, and 0.2153 for a linear `360/24` ramp. That is a margin
    /// of 1.6× below rather than the 6× an earlier version of this doc implied
    /// — it quoted 0.06, which is the figure for a ramp over ~86 groups, not 24.
    #[test]
    fn neighbouring_groups_are_far_apart_in_colour() {
        let of = |group: u8| colour(0.6, 3, 5.0, 0.0, 10.0, group);
        let mut worst = f64::INFINITY;
        let mut checked = 0_u32;
        for group in 0..MAX_GROUP {
            let d = separation(of(group), of(group + 1));
            if d < worst {
                worst = d;
            }
            checked += 1;
        }
        assert_eq!(
            checked,
            u32::from(MAX_GROUP),
            "the sweep did not run in full"
        );
        assert!(
            worst > 0.35,
            "neighbouring groups differ by only {worst:.4}, so a molecule whose \
             atoms sit at consecutive groups — which the demo molecule's leaves \
             do in every universe measured — is drawn as one colour"
        );
    }

    /// A closed outer shell is grey, and nothing else is.
    ///
    /// Fails on dropping the `valence == 0` arm, which is the edit that turns
    /// "bonds with nothing" into just another hue and takes Step 4b's inert
    /// story with it.
    #[test]
    fn only_a_closed_shell_is_grey() {
        for group in 0..12_u8 {
            let closed = colour(0.6, 0, 5.0, 0.0, 10.0, group);
            assert!(
                (closed[0] - closed[1]).abs() < 1e-12 && (closed[1] - closed[2]).abs() < 1e-12,
                "a closed shell at group {group} is {closed:?}, which is not grey"
            );
            for valence in 1..7_u8 {
                let bonding = colour(0.6, valence, 5.0, 0.0, 10.0, group);
                let grey = (bonding[0] - bonding[1]).abs() < 1e-12
                    && (bonding[1] - bonding[2]).abs() < 1e-12;
                assert!(
                    !grey,
                    "valence {valence} at group {group} is grey ({bonding:?}), so it \
                     is indistinguishable from a closed shell that bonds with nothing"
                );
            }
        }
    }

    /// The same properties give the same colour, always.
    ///
    /// The purity claim, asserted rather than assumed: `colour` reads
    /// nothing but its arguments, so two calls agree bit for bit.
    #[test]
    fn the_same_properties_give_the_same_colour() {
        let a = colour(0.61, 2, 7.5, 1.0, 9.0, 5);
        let b = colour(0.61, 2, 7.5, 1.0, 9.0, 5);
        assert_eq!(
            a.map(f64::to_bits),
            b.map(f64::to_bits),
            "the palette is not a pure function of its arguments"
        );
    }

    /// Every argument reaches the output.
    ///
    /// **An argument that changes nothing is a channel that carries nothing**,
    /// and the failure is silent: the colour still looks like a property map.
    /// Each arm moves exactly one argument and requires the colour to move.
    #[test]
    fn every_channel_reaches_the_colour() {
        let base = colour(0.5, 3, 5.0, 0.0, 10.0, 4);
        let moved = [
            ("affinity", colour(0.9, 3, 5.0, 0.0, 10.0, 4)),
            ("valence", colour(0.5, 0, 5.0, 0.0, 10.0, 4)),
            ("mass", colour(0.5, 3, 9.0, 0.0, 10.0, 4)),
            ("group", colour(0.5, 3, 5.0, 0.0, 10.0, 7)),
        ];
        for (name, c) in moved {
            assert_ne!(
                base.map(f64::to_bits),
                c.map(f64::to_bits),
                "moving {name} left the colour unchanged, so that channel carries \
                 nothing and the palette only looks like a property map"
            );
        }
    }

    /// A degenerate mass range gives the middle of the band, not `NaN`.
    ///
    /// A universe whose elements all share one mass has no spread to show. The
    /// failure this refuses is `0.0 / 0.0` painted onto every atom, which is
    /// not an error anywhere and renders as nothing at all.
    #[test]
    fn a_universe_with_no_mass_spread_paints_the_middle_of_the_band() {
        let c = colour(0.6, 3, 4.0, 4.0, 4.0, 3);
        for channel in c {
            assert!(channel.is_finite(), "a flat mass range gave {c:?}");
        }
        let midpoint = LIGHT_LO.midpoint(LIGHT_HI);
        // Lightness is recoverable from an HSL round trip as the midpoint of
        // the extreme channels.
        //
        // **Folded with explicit comparisons rather than `f64::max`/`f64::min`,
        // which §13.1 disallows** — they disagree about NaN between platforms,
        // and the deny reaches inside `#[cfg(test)]` because `clippy.toml` sets
        // no test exemption. A NaN here would sort to whichever end the
        // platform prefers and the assertion below would pass or fail by
        // architecture; the `is_finite` loop above is what makes that
        // impossible, and this spelling is what stops the next reader
        // reintroducing it.
        let mut hi = c[0];
        let mut lo = c[0];
        for channel in c {
            if channel > hi {
                hi = channel;
            }
            if channel < lo {
                lo = channel;
            }
        }
        assert!(
            (hi.midpoint(lo) - midpoint).abs() < 1e-9,
            "a flat mass range should sit at the middle of the lightness band"
        );
    }

    /// Hue wraps rather than reflecting, and the wheel is covered.
    ///
    /// Two arms, catching two different mutations, and **the second one had to
    /// be re-chosen after a probe**. `group` is a `u8` so no caller produces a
    /// negative hue today; this pins `hsl_to_rgb`'s own contract, against the
    /// day one does.
    ///
    /// - *past a full turn* (`390` against `30`) catches the modulo being
    ///   dropped altogether — without it, `390` lands in the sixth sector
    ///   rather than the first.
    /// - *below zero* (`-30` against `330`) catches `rem_euclid` being
    ///   swapped for a plain `%`, which follows the sign of the dividend.
    ///
    /// **The first version of the second arm compared `-330` against `30` and
    /// was silent under exactly that swap**, measured rather than reasoned:
    /// `-330 % 360` is `-330`, whose sector `-5.5` still selects the *first*
    /// arm of the colour chain, and the `ramp` term applies its own
    /// `rem_euclid` — so both spellings return the same colour and the test
    /// asserted nothing. At `-30` they diverge, because `330` selects the
    /// sixth sector and `-0.5` selects the first.
    #[test]
    fn the_hue_wheel_wraps() {
        let at = |h: f64| hsl_to_rgb(h, 0.7, 0.5);
        assert_eq!(
            at(30.0).map(f64::to_bits),
            at(390.0).map(f64::to_bits),
            "a hue past a full turn is not the same colour, so the wheel is not \
             being wrapped at all"
        );
        assert_eq!(
            at(330.0).map(f64::to_bits),
            at(-30.0).map(f64::to_bits),
            "a negative hue does not wrap onto the wheel — a plain `%` follows \
             the sign of the dividend, so it reflects instead"
        );
    }

    /// The golden angle is the constant it says it is.
    ///
    /// Pins the literal against a mistyped digit, which would weaken the
    /// separation above without failing it and would be invisible in review.
    #[test]
    fn the_golden_angle_is_the_golden_angle() {
        // 360 * (1 - 1/phi), phi = (1 + sqrt(5)) / 2. `sqrt` is exactly
        // specified by IEEE-754 and stays native under §13.1.
        let phi = f64::midpoint(1.0, 5.0_f64.sqrt());
        let want = 360.0 * (1.0 - 1.0 / phi);
        assert!(
            (GOLDEN_ANGLE - want).abs() < 1e-9,
            "the golden angle constant is {GOLDEN_ANGLE}, against {want}"
        );
    }
}
