//! Layer 0 — the generated chemistry of one universe (spec §7).
//!
//! A 64-bit seed becomes a complete periodic table. Nothing here is read from
//! a data file and nothing is modelled on real chemistry (§5, G1–G6).

pub mod bonds;
pub mod element;
pub mod naming;
/// **`pub(crate)`, not `pub`, and that is §8.6 enforced by the compiler rather
/// than by discipline.** Everything here runs at intern time and its results are
/// already materialised on [`element::Element`]; a step loop has no reason to
/// call any of it. `contacts_upto` costs 11.5 ns, so a redundant call from a hot
/// path would never show up in a profile — visibility is the only guard that
/// costs nothing.
#[expect(
    clippy::redundant_pub_crate,
    reason = "the items are `pub(crate)` *and* the module is, deliberately: if the \
              module is ever widened to `pub`, its contents must not silently become \
              public with it. The redundancy is the belt to the module's braces"
)]
pub(crate) mod packing;

// `NotABondOrder` travels with `BondOrder`: a public fallible conversion whose
// error type is not nameable from the same path forces a downstream `match` to
// reach into `bonds::` for one item.
pub use bonds::{BondEnergyMatrix, BondOrder, NotABondOrder};
pub use element::{Element, ElementId, PeriodicTable, ShellPattern};

use borbax_rng::{Domain, Stream};
use borbax_units::{Span, Thermal};

/// Universe-wide constants. Every one is generated, so two universes differ
/// in character rather than merely in their random draws.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct UniverseConsts {
    /// Target separation for a complementary fit in §8.3. Two surfaces
    /// interlock when their extents sum to about this.
    ///
    /// **`Span`, for the same reason `energy_per_unit` became `Quanta`** — it is
    /// a length crossing the public API, and [`Element::radius`] already
    /// establishes that lengths in this crate are `Span`. Its *value* is a
    /// separate open item (G1: it is on the wrong length scale), and typing it
    /// neither fixes nor hides that.
    pub ideal_gap: Span,
    /// Relative weight of geometric fit vs surface-character opposition.
    /// Dimensionless, so `f64` is the honest type rather than an omission.
    pub w_shape: f64,
    /// Weight of the surface-character term. Dimensionless.
    pub w_charge: f64,
    /// Pre-exponential factor in the Arrhenius rate law (§9.1).
    ///
    /// Bare `f64` and flagged rather than quietly accepted: this has dimension
    /// inverse-time, and there is no `Rate` newtype to hold it. Minting one is
    /// §9.1's business — it needs the rate law in front of it to know whether
    /// the unit is per [`borbax_units::WorldYear`] or per collision — so the
    /// decision routes to Task 14 (reaction rates) rather than being guessed
    /// here. Task 9 is the shape signature; §9.1 is Task 14's, and reading the
    /// spec section number as a task number is how §9 items get misrouted.
    pub rate_prefactor: f64,
    /// Global spontaneous-cleavage scale. V0 locates the decay band by
    /// moving this single number (spec §22.6); per-bond rates come later.
    /// Dimensionless multiplier on a rate, not a rate itself.
    pub decay_scale: f64,
    /// Which element acts as solvent. Its complementarity against a bond
    /// determines solvent attack (§9.5).
    pub solvent: ElementId,
    /// Coldest temperature this universe supports.
    pub temp_min: Thermal,
    /// Hottest temperature this universe supports.
    pub temp_max: Thermal,
}

/// The laws a universe was generated under.
///
/// **A universe is `(seed, physics)`, not `seed`.** §6 promises that
/// sharing `U-7F3A21C9` means the recipient gets the same planet; today the
/// real identity is `(seed, git SHA, rustc version)` and only the first third
/// is written down. This field writes down the rest of what matters.
///
/// The cost is stated honestly: every version here is a code path kept alive
/// forever, with its own goldens. That is the price of never breaking a world
/// somebody cared about, and it is bounded — versions are added when physics
/// changes, which should be roughly never.
///
/// **Gate on this, never on a zero weight.** The tempting shortcut for an
/// added term is to give it weight `0.0` for old universes and claim they are
/// unchanged. They are not: `0.0 * x` is `NaN` when `x` is infinite or NaN,
/// and a stray NaN's sign bit differs between the CI architectures (see
/// `det_math`). An inert term must be *skipped*:
///
/// ```text
/// if u.physics >= PhysicsVersion::V2 { score += w_field * field_term; }   // yes
/// score += u.consts.w_field * field_term;                                 // no
/// ```
///
/// New constants for a new version come from a **new `Domain`** (Task 3), so
/// adding them leaves every draw an existing universe already made untouched.
/// Those two rules together are what make a later law additive rather than
/// destructive.
///
/// # The admission criterion, which is not the same as the cost
///
/// This documents what a new version *costs* and offers no rule for when one
/// is *allowed*, and that gap is load-bearing. Before the seam existed, adding
/// a law carried a deterrent that was policing §3.1 by accident: it invalidated
/// every world ever generated, which is a visible, expensive, argued-about
/// event. Now it is a checklist — new `Domain`, new variant, regenerate
/// goldens. The seam is right (a world silently changing meaning is worse), but
/// it removed friction that a different invariant was relying on.
///
/// The `if physics >= V2 { score += w * term }` idiom makes it sharper: that is
/// monotone accretion inside `affinity()`, and two versions in, the scoring
/// function is a list of version-gated additive terms that all look alike.
/// "V2 physics has a term that helps membranes form" would be one more line
/// among them.
///
/// **So: a new `PhysicsVersion` must state the physical process it adds in
/// terms that name no biological structure.** "An environmental scalar enters
/// the fold energy" passes. "A rule that helps membranes form" fails on its
/// face, because it names its own outcome. That is §3.1's deletion test applied
/// to the justification rather than the code, and it is the thing that makes
/// this seam safe to have.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum PhysicsVersion {
    /// Shape and surface character, per spec §8.3. Everything in V0 and V1.
    V1 = 1,
}

/// The wire form of a version, for digests and for `U-7F3A21C9@1` addresses.
///
/// Spelled out rather than `as`, for the reason [`BondOrder`]'s inverse is:
/// `as` on a `#[non_exhaustive]` enum silently accepts a future variant, and a
/// version number reaching an address is the last place that should happen.
impl From<PhysicsVersion> for u8 {
    #[expect(
        clippy::as_conversions,
        reason = "the hazard this impl exists to avoid is `v as u8` on an opaque \
                  value, which silently accepts a future variant. Inside an \
                  exhaustive match arm the variant is already known, so writing \
                  the discriminant this way makes the wire form and the derived \
                  `Ord` one fact instead of two encodings that can drift — \
                  `V2 = 2` with a hand-written `=> 3` here compiles, orders \
                  correctly, and writes the wrong number into both the digest \
                  and the `U-7F3A21C9@N` address"
    )]
    fn from(v: PhysicsVersion) -> Self {
        match v {
            PhysicsVersion::V1 => PhysicsVersion::V1 as Self,
        }
    }
}

impl PhysicsVersion {
    /// The laws a *new* universe is generated under.
    ///
    /// Read only by [`Universe::generate`]. A universe being reloaded passes
    /// its own recorded version to [`Universe::generate_under`] instead — that
    /// asymmetry is the whole seam.
    pub const CURRENT: Self = Self::V1;
}

/// One generated universe: the laws of physics for one `universe_seed`.
///
/// Generated once, then immutable. Everything downstream — worlds, molecules,
/// reactions — reads these constants and never changes them.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Universe {
    /// The seed this universe was generated from.
    pub seed: u64,
    /// The laws this universe was generated under. Never inferred from the
    /// running binary — a universe loaded from disk keeps the physics it was
    /// born with, or it is a different universe wearing the same name.
    pub physics: PhysicsVersion,
    /// The elements, and the shell law they follow.
    pub table: PeriodicTable,
    /// Pairwise bond energies, derived from the table's binding energies.
    pub bonds: BondEnergyMatrix,
    /// The drawn constants of this universe.
    pub consts: UniverseConsts,
}

impl Universe {
    /// Generate a **new** universe, under the laws this binary implements.
    #[must_use]
    pub fn generate(seed: u64) -> Self {
        Self::generate_under(seed, PhysicsVersion::CURRENT)
    }

    /// Regenerate a universe under stated laws.
    ///
    /// **This is the half of the seam that is easy to leave out, and without it
    /// the `physics` field guards nothing.** A universe is *derived* from its
    /// seed rather than serialised, so "reloading" one means re-running
    /// generation — and with only [`Self::generate`] there is no argument that
    /// could say which laws to apply. A `V2` binary would stamp `V2` on every
    /// world ever shared, which is precisely the failure `physics` exists to
    /// prevent. The address parser feeds `@1` into this one.
    #[must_use]
    pub fn generate_under(seed: u64, physics: PhysicsVersion) -> Self {
        // Exhaustive *within this crate* despite `#[non_exhaustive]`, which is
        // the forcing function: adding `V2` breaks this match at compile time,
        // so a new law cannot be added without deciding what old worlds do.
        match physics {
            PhysicsVersion::V1 => Self::generate_v1(seed, physics),
        }
    }

    fn generate_v1(seed: u64, physics: PhysicsVersion) -> Self {
        // `generate_elements` returns the shell it actually used. An earlier
        // draft re-derived it by rebuilding the same stream and replaying —
        // which agreed by luck, and would silently stop agreeing the moment
        // anyone added a draw before the shell call inside `generate_elements`.
        // `Universe`'s pattern would then describe a table that does not exist,
        // with nothing failing. This is the class of bug dropping `Copy` on
        // `Stream` (Task 3) exists to surface.
        let table = element::generate_elements(seed);

        // A separate stream index so that adding a draw in element generation
        // cannot shift the constants, and vice versa.
        let mut rng = Stream::new(seed, Domain::Universe, 1);
        let bonds = BondEnergyMatrix::generate(&table, &mut rng);

        // Solvent: drawn uniformly from the lightest third of the table.
        //
        // **The comment this replaced said "a small, abundant, strongly-polar
        // element", and the code reads no polarity and no abundance.** It is a
        // uniform draw over a prefix — light and abundant follow from the
        // prefix, because `abundance` falls geometrically in `units`, but
        // "strongly-polar" was a property nothing here selects for. Saying so
        // matters more than it looks: §9.5 makes solvent attack a function of
        // the solvent's complementarity against a bond, and a reader who
        // believes polarity was already selected here will not think to check.
        // `try_from`/`unwrap_or` rather than `as`: both conversions are
        // infallible for a 60..=120 table on any supported target, and neither
        // fallback is reachable — but `clippy::as_conversions` is denied
        // workspace-wide precisely because a silent narrowing is the gotcha most
        // likely to produce a plausible wrong number.
        let solvent_pool = u64::try_from((table.len() / 3).max(1)).unwrap_or(1);
        let solvent = u8::try_from(rng.next_range(solvent_pool)).map_or(ElementId::ZERO, ElementId);

        let temp_min = Thermal(rng.next_f64_range(180.0, 260.0));
        let temp_max = temp_min + Thermal(rng.next_f64_range(200.0, 500.0));

        let consts = UniverseConsts {
            ideal_gap: Span(rng.next_f64_range(1.6, 2.8)),
            w_shape: rng.next_f64_range(0.6, 1.4),
            w_charge: rng.next_f64_range(0.6, 1.4),
            rate_prefactor: rng.next_f64_range(1e3, 1e5),
            decay_scale: rng.next_f64_range(1e-6, 1e-4),
            solvent,
            temp_min,
            temp_max,
        };

        Self {
            seed,
            // Carried in, never read from `PhysicsVersion::CURRENT` here. That
            // is the difference between a universe recording its laws and a
            // universe being restamped by whichever binary last opened it.
            physics,
            table,
            bonds,
            consts,
        }
    }

    /// The element with this id, or `None` if it is not from this universe.
    #[must_use]
    pub fn element(&self, id: ElementId) -> Option<&Element> {
        self.table.get(id)
    }

    /// This universe's solvent (§9.5).
    #[must_use]
    pub fn solvent(&self) -> Option<&Element> {
        self.element(self.consts.solvent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn universes_are_deterministic() {
        // Whole-value equality, not three sampled fields: a new field added to
        // `Universe` is covered the day it lands rather than the day someone
        // remembers to extend this list.
        assert_eq!(Universe::generate(77), Universe::generate(77));
    }

    /// A universe records the laws it was born under.
    ///
    /// This looks trivial and is not. The failure it guards against is a
    /// future `generate` that stamps whatever version the running binary
    /// happens to implement — at which point loading an old world under new
    /// laws silently reinterprets it, and `U-7F3A21C9` stops meaning what
    /// §6 promises. The real test arrives with `PhysicsVersion::V2`:
    /// `Universe::generate_under` must still stamp `V1` for a world being
    /// *reloaded*, and only new worlds get `CURRENT`.
    #[test]
    fn a_universe_records_its_physics() {
        assert_eq!(Universe::generate(77).physics, PhysicsVersion::V1);
    }

    /// The half of the seam that `generate` alone cannot provide. With only
    /// `generate(seed)` there is no argument that could say which laws to
    /// apply, so a `V2` binary would restamp every reloaded world.
    #[test]
    fn generate_under_stamps_what_it_is_told() {
        let u = Universe::generate_under(77, PhysicsVersion::V1);
        assert_eq!(u.physics, PhysicsVersion::V1);
        assert_eq!(u, Universe::generate(77));
    }

    #[test]
    fn element_lookup_round_trips() {
        let u = Universe::generate(3);
        for (id, e) in u.table.iter() {
            assert_eq!(
                u.element(id).map(|f| f.symbol.as_str()),
                Some(e.symbol.as_str())
            );
        }
    }

    /// The panic site decision 3 closes. Clamping an out-of-range id — the
    /// shape this replaced — returns a *different element's* properties, which
    /// is worse than a crash because nothing downstream can tell.
    #[test]
    fn an_out_of_range_id_is_none_rather_than_a_neighbour() {
        let u = Universe::generate(3);
        // `unwrap_or`, not `expect`: `clippy.toml` sets neither
        // `allow-expect-in-tests` nor `allow-unwrap-in-tests`, so the workspace
        // deny reaches inside `#[cfg(test)]`. The fallback is out of range too.
        let past_the_end = ElementId::from_index(u.table.len()).unwrap_or(ElementId(u8::MAX));
        assert!(u.element(past_the_end).is_none());
    }

    /// **Pins every value Task 5 added, because nothing else can see them.**
    ///
    /// `the_universe_digest_is_pinned` covers `ShellPattern` and `Element` and
    /// stops there. Task 5 added ten draws on `Domain::Universe` index 1 — four
    /// in [`BondEnergyMatrix::generate`] and six more for the constants — and
    /// none of them was in any digest. Three review lanes found the same gap.
    ///
    /// The hazard is not cross-platform, it is **cross-version**: `w_shape` and
    /// `w_charge` draw from the identical range, so transposing two fields in
    /// the `UniverseConsts` literal rewrites every universe while leaving every
    /// value inside its documented range, every type check passing and every
    /// range assertion green. A determinism lane verified exactly that — the
    /// transposition moved all eight probed seeds and **all six gate legs stayed
    /// green**. A seed shared before such an edit names a different universe
    /// after it, with `PhysicsVersion` still stamped `V1`, which is the failure
    /// §13.1's identity tuple exists to exclude.
    ///
    /// It is worse than it looks, because the literal's field order already
    /// disagrees with the draw order — `solvent` and the temperatures are bound
    /// above it. So "align the literal with the order things are drawn" reads as
    /// pure tidying and is a physics change.
    ///
    /// Raw `to_bits()`, deliberately, and not `canonical_bits()`: this is a
    /// **detector**, and canonicalisation is lossy by design — it would collapse
    /// a `-0.0` that differs by architecture onto `+0.0` and report "unchanged"
    /// while the value genuinely differed. Lossiness is right for a product
    /// hash and wrong here.
    #[test]
    fn the_assembled_universe_digest_is_pinned() {
        // FNV-1a: XOR-then-multiply, so the low k bits of the product depend
        // only on the low k bits of the operands — there is essentially no
        // avalanche downward. Measured: a mutation touching only exponent bits
        // of every Single cell left the low 32 bits identical. Compare the full
        // `u64` only — never a prefix, never truncated to a display address.
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |x: u64| {
            h ^= x;
            h = h.wrapping_mul(0x0100_0000_01b3);
        };
        for seed in 0..64 {
            let u = Universe::generate(seed);
            // **Destructured, not field-accessed, and that is the whole guard.**
            // A review added a ninth drawn constant three ways — appended to the
            // stream, and on a fresh `Domain::Universe` index as this enum's own
            // doc *recommends* — and all 54 tests stayed green each time. A
            // hand-written field list cannot see a field that is not in it.
            // `#[non_exhaustive]` does not bind in-crate, so destructuring makes
            // an added field `E0027` here instead of a value silently outside
            // the digest — the same forcing function as `generate_under`'s match.
            let Universe {
                seed: s,
                physics,
                table,
                bonds,
                consts,
            } = &u;
            let UniverseConsts {
                ideal_gap,
                w_shape,
                w_charge,
                rate_prefactor,
                decay_scale,
                solvent,
                temp_min,
                temp_max,
            } = consts;
            mix(*s);
            // `u8::from` on the discriminant rather than `as`: the enum is
            // `#[non_exhaustive]` and a future variant must not silently widen.
            mix(u64::from(u8::from(*physics)));
            // Fixed order, this test's own — so a reorder of the struct literal
            // cannot reorder the digest with it.
            mix(ideal_gap.0.to_bits());
            mix(w_shape.to_bits());
            mix(w_charge.to_bits());
            mix(rate_prefactor.to_bits());
            mix(decay_scale.to_bits());
            mix(u64::try_from(solvent.index()).unwrap_or(u64::MAX));
            mix(temp_min.0.to_bits());
            mix(temp_max.0.to_bits());
            // Every bond cell at every order, covering `e` and both order
            // multipliers.
            let ids: Vec<ElementId> = table.iter().map(|(id, _)| id).collect();
            for &a in &ids {
                for &b in &ids {
                    // Every representable order, not the first three: the cap
                    // is now `BondOrder::MAX`, derived from the valence ceiling.
                    for o in 1..=BondOrder::MAX {
                        if let Some(order) = BondOrder::new(o) {
                            mix(bonds.energy(a, b, order).get().to_bits());
                        }
                    }
                }
            }
        }
        assert_eq!(
            h, 0x59a6_b8d5_560f_ea5c,
            "the assembled-universe digest moved — say which of §18.1's three this \
             is, or the fourth: the digest's own seed range widened, which moves \
             the constant while moving no universe value. Recomputing over the \
             previous range distinguishes them, and did — this constant's \
             predecessor reproduces exactly at 0..16. Otherwise check whether a \
             draw was reordered or inserted"
        );
    }

    #[test]
    fn solvent_is_a_real_element_of_this_universe() {
        for seed in 0..20 {
            let u = Universe::generate(seed);
            assert!(u.solvent().is_some(), "seed {seed} has no solvent");
        }
    }
}
