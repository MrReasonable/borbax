//! Layer 0 — the generated chemistry of one universe (spec §7).
//!
//! A 64-bit seed becomes a complete periodic table. Nothing here is read from
//! a data file (§5, G1). **"Nothing here is modelled on real chemistry" no
//! longer holds without qualification** — `naming::REAL_TABLE` ships the
//! real periodic table verbatim, for the one seed-equivalent identity
//! configuration G2's 2026-08-07 revision carves out; every *perturbed*
//! universe is exactly as unmodelled as before. G5 (cited here as "G1–G6"
//! in an earlier version of this line) was withdrawn 2026-08-06; see §5 for
//! the current, revised set.
//!
//! # The migrated surface, as an executable table
//!
//! Task 5 renamed or retyped most of this crate's public API, and the plan files
//! carried hand-written "was / is now" tables to tell downstream tasks what
//! moved. Those tables went stale **in the commit that changed the type** —
//! `BondOrder::MAX` went 3 to 6 and both copies still spelled a
//! `BondOrder::Single` that no longer existed. They failed in the worst
//! available way: correct in form, wrong in content, so they looked maintained.
//!
//! A table is a copy of a fact and nothing compiles it. This is the same fact,
//! and `cargo test --locked --workspace` runs it — so a rename is a build
//! failure rather than a stale row. The plan's tables keep their "why" column,
//! which is genuinely a plan's job, and point here for "is now".
//!
//! ```
//! use borbax_universe::{
//!     Universe,
//!     bonds::BondOrder,
//!     element::{ElementId, PeriodicTable, generate_elements},
//! };
//! use borbax_units::{Quanta, Span};
//!
//! let u = Universe::generate(7);
//!
//! // `ElementId(3)` from another crate is a compile error — the field is
//! // `pub(crate)`. Construct through `from_index`, which is fallible.
//! let id = ElementId::from_index(3).expect("3 fits in a u8");
//!
//! // `u.elements` / `u.shell` are `u.table`, a `PeriodicTable`.
//! let table: &PeriodicTable = &u.table;
//! assert!(table.get(id).is_some());
//! assert!(table.iter().count() >= 1);
//! // `pattern()` is public; `eps` on it is not — it is the bond scale's input
//! // and lives inside the crate. `k`, `closures` and `peak` are the public part.
//! assert!(table.pattern().k >= 1);
//!
//! // `u.element(id)` answers `Option`, never a neighbouring element.
//! let el = u.element(id).expect("seed 7 has at least four elements");
//! let _epu: Quanta = el.energy_per_unit;
//! assert!(u.element(ElementId::from_index(255).expect("fits")).is_none());
//!
//! // `energy` takes ids and a validated `BondOrder`, not `&Element` and `u8`.
//! let _e: Quanta = u.bonds.energy(id, id, BondOrder::SINGLE);
//! assert!(BondOrder::new(0).is_none());
//! assert!(BondOrder::try_from(BondOrder::MAX).is_ok());
//! assert_eq!(BondOrder::ALL.len(), usize::from(BondOrder::MAX));
//! assert_eq!(BondOrder::PLANES, usize::from(BondOrder::MAX));
//! assert_eq!(BondOrder::SINGLE.plane_index(), 0);
//!
//! // G4: lengths crossing the boundary are `Span`, not `f64`.
//! let _gap: Span = u.consts.ideal_gap;
//! ```

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
/// **`pub(crate)`, not `pub` — narrow-API hygiene, not a fiction-guarantee
/// fix.** `Rung` and `Direction` have no consumer outside this crate
/// today (verified by grep across the whole workspace), and `Rung::draw`
/// being reachable let any crate compute `Rung::draw(seed).is_identity()`
/// without ever touching [`naming::IdentityWitness`]. That is worth closing
/// because it costs nothing, not because G2 is at stake: `Universe::seed` is
/// already `pub`, so the same predicate is already reconstructible today via
/// `borbax_rng::Domain::PerturbationRung` regardless of this module's own
/// visibility, and G2's 2026-08-07 revision does not treat the identity
/// configuration as a secret in the first place.
#[expect(
    clippy::redundant_pub_crate,
    reason = "same reasoning as packing, above: the items are `pub(crate)` *and* the module \
              is, deliberately, so a future widening of the module to `pub` cannot silently \
              carry its contents public with it"
)]
pub(crate) mod perturbation;
/// A typed phrase in, a seed out — and the phrase goes no further.
///
/// Deliberately **not** part of [`naming`], which is about *generated* element
/// names and the G2 blocklist that constrains them. Nothing a user types is
/// subject to G2.
pub mod phrase;

// `NotABondOrder` travels with `BondOrder`: a public fallible conversion whose
// error type is not nameable from the same path forces a downstream `match` to
// reach into `bonds::` for one item.
pub use bonds::{BondEnergyMatrix, BondOrder, NotABondOrder};
pub use element::{Element, ElementId, PeriodicTable, ShellPattern};
pub use phrase::seed_from_phrase;

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
    ///
    /// **An `ElementId` because Layer 0 has no molecules, not because a solvent
    /// is an element.** A universal solvent is far more plausibly a *molecule* —
    /// and by the time a beaker is filled, a pre-generated molecular inventory
    /// already exists (dust-cloud chemistry, plus whatever planetary heat
    /// drives), which is where a real solvent would come from. Making this a
    /// `SpeciesId` requires `borbax-molecule`, which sits above this crate, so
    /// it cannot be done here at any price.
    ///
    /// Recorded so the constraint is visible before §9.5's solvent attack is
    /// written against the element reading and inherits it. The decision belongs
    /// with Task 6, and Task 5b's scope block states it.
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
    fn from(v: PhysicsVersion) -> Self {
        v.discriminant()
    }
}

impl PhysicsVersion {
    /// The laws a *new* universe is generated under.
    ///
    /// Read only by [`Universe::generate`]. A universe being reloaded passes
    /// its own recorded version to [`Universe::generate_under`] instead — that
    /// asymmetry is the whole seam.
    pub const CURRENT: Self = Self::V1;

    /// Every variant, in discriminant order — the same shape as
    /// `borbax_rng::Domain::ALL`.
    ///
    /// **This is what lets `the_assembled_universe_digest_is_pinned` pin V1's
    /// digest by name instead of by `CURRENT`.** Before this existed, that
    /// test called [`Universe::generate`], which resolves through `CURRENT` —
    /// so the moment a second variant landed and `CURRENT` moved to it, the
    /// test would still fail, but for the wrong reason its own message names
    /// (a seed/order range widening, or a reordered draw), and "fixing" it by
    /// re-pinning the literal to the new digest would silently discard V1's
    /// own golden value forever, even though `V1` stays fully live and
    /// reachable via `generate_under` for every world generated before the
    /// bump. Iterating `ALL` and matching each variant to its own golden
    /// value, exhaustively, means a new variant with no arm is `E0004` at
    /// compile time — the same forcing function `Domain::ALL` already gives
    /// `every_domain_stream_is_pinned`.
    ///
    /// **What this does not enforce, stated precisely, the same caveat
    /// `Domain::ALL` carries — this *is* a hand-written array, and the
    /// compile-time forcing above is narrower than "checks against the
    /// enum's actual variants."** Appending a variant forces a golden arm to
    /// exist so the crate compiles — it does not force the new variant to
    /// also reach `ALL`. Omitted there, that arm still has to exist, but the
    /// loop never reaches it, so a wrong golden value goes unexercised until
    /// someone remembers to add the variant here too. `PhysicsVersion::CURRENT
    /// missing from ALL` is the one shape of that gap this module can close
    /// without a proc-macro dependency (`const _`, below) — `CURRENT` is the
    /// physics *every new universe* is generated under, so a universe with no
    /// pinned digest for its own laws is the failure that actually matters;
    /// "V2 exists, CURRENT is still V1, and V2 is missing from ALL" is not
    /// closed by anything here, and is not closed by `Domain::ALL` either.
    pub const ALL: [Self; 1] = [Self::V1];

    /// The single place the discriminant cast is written. `From<Self> for
    /// u8` delegates here rather than carrying its own second copy of the
    /// same match and the same `#[expect]` — two exhaustive matches over one
    /// enum, one hand-written cast each, is exactly the "two encodings that
    /// can drift" this fn's own `#[expect]` reason warns against, applied to
    /// itself.
    #[must_use]
    #[expect(
        clippy::as_conversions,
        reason = "the hazard this fn exists to avoid is `v as u8` on an opaque \
                  value, which silently accepts a future variant. Inside an \
                  exhaustive match arm the variant is already known, so writing \
                  the discriminant this way makes the wire form and the derived \
                  `Ord` one fact instead of two encodings that can drift — \
                  `V2 = 2` with a hand-written `=> 3` here compiles, orders \
                  correctly, and writes the wrong number into both the digest \
                  and the `U-7F3A21C9@N` address"
    )]
    const fn discriminant(self) -> u8 {
        match self {
            Self::V1 => Self::V1 as u8,
        }
    }
}

/// **`CURRENT` must be swept by `the_assembled_universe_digest_is_pinned`.**
/// Promoting a variant to `CURRENT` without adding it to `ALL` leaves the
/// physics every *new* universe is generated under with no pinned digest,
/// while every other test — including the digest test itself — stays green:
/// `ALL`'s own doc comment records that a golden arm is forced to *exist*,
/// not to *run*. A `const` block, not a runtime test, so the gap cannot
/// survive even a single `cargo test` invocation that happens to skip this
/// one function.
#[expect(
    clippy::indexing_slicing,
    reason = "a const context: `.get(i)` returns an `Option` this loop would have to \
              `unwrap()`, which is exactly as deniable, and the index is provably in \
              bounds — `i < PhysicsVersion::ALL.len()` is the loop condition, checked \
              before every indexing operation reached, never after"
)]
const _: () = {
    let mut i = 0;
    let mut found = false;
    while i < PhysicsVersion::ALL.len() {
        if PhysicsVersion::ALL[i].discriminant() == PhysicsVersion::CURRENT.discriminant() {
            found = true;
        }
        i += 1;
    }
    assert!(
        found,
        "PhysicsVersion::CURRENT is missing from PhysicsVersion::ALL — the physics \
         every new universe is generated under has no pinned digest"
    );
};

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
        let table = element::generate_elements(seed, physics);

        // **Three indices, not two, and the third was bought with a measured
        // defect.** The comment here used to read "a separate stream index so
        // that adding a draw in element generation cannot shift the constants,
        // and vice versa" — true of element generation, false of bonds, which
        // shared index 1 with the constants and sat *above* them on it.
        //
        // Retiring `weakest_share` and merging two order multipliers into one
        // exponent took the bond generator from four draws to two, and a review
        // measured the consequence: `solvent`, `temp_min`, `temp_max` and all
        // five `UniverseConsts` reals moved in **8 of 8** probed seeds. Seed 0's
        // `temp_min` went 206.09 -> 222.68. Every value stayed inside its
        // documented range, so nothing but a digest could see it.
        //
        // §13.1: the identity of a universe must not depend on how many numbers
        // an unrelated derivation happened to need. CLAUDE.md names the decay
        // band — downstream of these temperatures — as the most sensitive
        // parameter in the system, so the next person to add a bond constant
        // would have spent a week on a dead run.
        let bonds = BondEnergyMatrix::generate(&table, &mut Stream::new(seed, Domain::Universe, 2));
        let mut rng = Stream::new(seed, Domain::Universe, 1);

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
    ///
    /// **No comparison against `Universe::generate` here** — an earlier
    /// version of this test also asserted `u ==
    /// Universe::generate(77)`, which is exactly the coupling P2 (issue
    /// #26's prerequisites) removed from
    /// `the_assembled_universe_digest_is_pinned`: it agrees only because
    /// `CURRENT == V1` today, and at the V2 bump it fails naming no cause —
    /// the two natural repairs are deleting the assertion (correct) or
    /// changing `V1` to `CURRENT` (tautological, and it would silently
    /// retire the delegation check). `generate_follows_current` is that
    /// comparison's one home now.
    #[test]
    fn generate_under_stamps_what_it_is_told() {
        assert_eq!(
            Universe::generate_under(77, PhysicsVersion::V1).physics,
            PhysicsVersion::V1
        );
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
    ///
    /// **Pinned by [`PhysicsVersion::ALL`], one golden hash per variant, not by
    /// [`PhysicsVersion::CURRENT`].** An earlier version of this test called
    /// `Universe::generate(seed)`, which resolves through `CURRENT` — so the
    /// literal it pinned was never actually "V1's digest", it was "whatever
    /// `CURRENT` happens to point at today". P2 (issue #26's prerequisites)
    /// found the failure mode by mutation: add a byte-identical `V2` variant,
    /// move `CURRENT` to it, and the old test failed — but its own message
    /// named a seed-range widening, an order-range widening, or a reordered
    /// draw, none of which had happened, and "fixing" it the natural way (
    /// re-pinning the literal to the new hash) would have permanently
    /// discarded V1's golden value while `V1` stayed fully reachable via
    /// `generate_under` for every world already generated under it. Iterating
    /// `PhysicsVersion::ALL` and calling `generate_under` explicitly per
    /// version decouples "what does V1 digest to" from "what does `CURRENT`
    /// point at this week" — and the exhaustive `match` below means a new
    /// variant with no golden arm is `E0004` at compile time, not a value
    /// silently outside the pin.
    #[test]
    fn the_assembled_universe_digest_is_pinned() {
        for version in PhysicsVersion::ALL {
            // FNV-1a: XOR-then-multiply, so the low k bits of the product
            // depend only on the low k bits of the operands — there is
            // essentially no avalanche downward. Measured: a mutation
            // touching only exponent bits of every Single cell left the low
            // 32 bits identical. Compare the full `u64` only — never a
            // prefix, never truncated to a display address.
            let mut h: u64 = 0xcbf2_9ce4_8422_2325;
            let mut mix = |x: u64| {
                h ^= x;
                h = h.wrapping_mul(0x0100_0000_01b3);
            };
            let (mut mixed, mut pairs) = (0_usize, 0_usize);
            for seed in 0..64 {
                let u = Universe::generate_under(seed, version);
                // **Destructured, not field-accessed, and that is the whole
                // guard.** A review added a ninth drawn constant three ways —
                // appended to the stream, and on a fresh `Domain::Universe`
                // index as this enum's own doc *recommends* — and all 54
                // tests stayed green each time. A hand-written field list
                // cannot see a field that is not in it. `#[non_exhaustive]`
                // does not bind in-crate, so destructuring makes an added
                // field `E0027` here instead of a value silently outside the
                // digest — the same forcing function as `generate_under`'s
                // own match.
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
                // **The invariant P3's whole design rests on, pinned here
                // rather than merely relied upon.** `BondEnergyMatrix::generate`
                // dispatches on `table.physics()`, while the universe's own
                // address, digest and every downstream reader use `u.physics`
                // — `generate_v1` holds them equal by construction (one
                // `physics` binding feeds both), but nothing before this
                // assertion checked it. A future version whose derivation
                // touched some other subsystem and left this equality broken
                // would be undetectable by the digest alone.
                assert_eq!(
                    table.physics(),
                    *physics,
                    "version {version:?} seed {seed}: the table's laws disagree with \
                     the universe's own"
                );
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
                        //
                        // (This iterated `1..=MAX` through `BondOrder::new` and
                        // carried a paragraph about its `None` arm. It is
                        // `for order in ALL` now — there is no `Option` here.)
                        for order in BondOrder::ALL {
                            let v = bonds.energy(a, b, order).get();
                            // **A non-finite value must name itself here.** Both
                            // unreachable fallbacks in `bonds.rs` are `NAN`, and a
                            // review measured that NaN *propagation* through `*`,
                            // `/` and `sqrt` is bit-identical on aarch64 and x86-64
                            // — so a poisoned universe would move this hash the same
                            // way on all three CI legs, the §13.4 matrix would stay
                            // green, and "category one, regenerate" is the natural
                            // and wrong conclusion. Only *freshly generated* NaNs
                            // diverge by architecture.
                            assert!(
                                v.is_finite(),
                                "version {version:?} seed {seed}: {a:?}-{b:?} order \
                                 {order:?} is {v} — a non-finite value reached the \
                                 digest. This is NOT one of §18.1's three; an \
                                 unreachable fallback fired"
                            );
                            mix(v.to_bits());
                            mixed += 1;
                        }
                        pairs += 1;
                    }
                }
            }
            // Catches exactly one thing: a `continue` added inside the order
            // loop, which would silently hash fewer cells. `mixed` and
            // `pairs` exist only for this. `ALL`'s *contents* are guarded by
            // a `const _` in `bonds.rs`, at compile time.
            assert_eq!(
                mixed,
                pairs * BondOrder::ALL.len(),
                "version {version:?}: the digest skipped an order — `{mixed}` cells \
                 for `{pairs}` pairs"
            );
            let want = match version {
                PhysicsVersion::V1 => 0xef29_79c6_1879_20d8,
            };
            assert_eq!(
                h, want,
                "version {version:?}'s assembled-universe digest moved — say which \
                 of §18.1's three this is, or one of two widenings that move the \
                 constant while moving no universe value: the seed range (0..64 \
                 today) or the order range (1..=BondOrder::MAX today). Recompute \
                 over the previous range for whichever changed. Otherwise check \
                 whether a draw was reordered or inserted — and note that since \
                 the bond generator moved to its own stream index, its draw count \
                 can no longer do that silently, which it previously could and did"
            );
        }
    }

    /// **`Universe::generate`'s entire job is delegating to `generate_under`
    /// under `CURRENT`, and this test replaces the delegation check that
    /// used to live inside `generate_under_stamps_what_it_is_told`** — moved
    /// out, not newly added: that test's own `assert_eq!(u,
    /// Universe::generate(77))` was the same CURRENT-coupled comparison,
    /// under a name that didn't say so. [`the_assembled_universe_digest_is_pinned`]
    /// no longer calls `generate` at all — it iterates
    /// [`PhysicsVersion::ALL`] through `generate_under` explicitly, by
    /// design (see that test's own doc comment) — so this is now the one
    /// place in the file that would catch a `generate` body drifting from
    /// `generate_under(seed, CURRENT)`.
    ///
    /// **The comparison is `Universe`'s derived `PartialEq`, not a digest,
    /// and that has one live direction worth naming: a `NaN` on either side
    /// compares unequal to itself, so a poisoned universe (either of
    /// `bonds.rs`'s two unreachable `NAN` fallbacks firing) would fail this
    /// assertion with a message blaming the delegation, not the poison.**
    /// `-0.0`/`+0.0` is not a live hazard here the way it is for the digest
    /// test's `to_bits()`: that distinction exists to catch a value that
    /// differs *by architecture*, and both sides of this comparison run in
    /// the same process on the same architecture. Seeds `[0, 1, 63]` all sit
    /// inside `the_assembled_universe_digest_is_pinned`'s own `0..64`
    /// finiteness sweep, so this cannot fire today.
    #[test]
    fn generate_follows_current() {
        for seed in [0, 1, 63] {
            assert_eq!(
                Universe::generate(seed),
                Universe::generate_under(seed, PhysicsVersion::CURRENT),
                "seed {seed}: Universe::generate no longer matches \
                 generate_under(seed, CURRENT)"
            );
        }
    }

    #[test]
    fn solvent_is_a_real_element_of_this_universe() {
        for seed in 0..20 {
            let u = Universe::generate(seed);
            assert!(u.solvent().is_some(), "seed {seed} has no solvent");
        }
    }
}
