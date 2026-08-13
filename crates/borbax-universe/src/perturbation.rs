//! The shared perturbation mechanism (issue #26, Decision 10) — one integer
//! rung per universe, and a bounded per-constant direction, composed as
//! `base * (1.0 + strength * p)`.
//!
//! **`Rung` and `Direction` landed with zero in-crate callers; Task 26.1
//! gave them real ones, in two stages.** Steps 1-5:
//! [`MigratedConstant::ScreeningInner`](crate::perturbation::MigratedConstant::ScreeningInner),
//! `orbital.rs`'s screened-hydrogenic fill drawing its inner-shell screening
//! coefficient one-sided downward (see [`MigratedConstant`](crate::perturbation::MigratedConstant)'s own doc for
//! why only this one direction is legal for that particular constant). Steps
//! 6-9: `Universe.rung` itself, drawn for every physics version (costs V1
//! nothing — `Rung::draw` uses its own dedicated `Domain::PerturbationRung`,
//! entirely separate from `Domain::Universe`'s own indices), plus the rest
//! of the plan's migration table — `PromotionBudget`, `RadiusScale`,
//! `CapacityScale`, `BaseMass`, `ContactDefect`, `WShape`, `WCharge`,
//! `IdealGap`, `Scale`, `Gamma` — each landing only once its own real
//! consumer existed, never ahead of it (the "mechanism ahead of its
//! subject" shape this module's own earlier history, and P1's, and P4's,
//! already exist to avoid).
//!
//! **`Rung`'s own visibility, corrected in `/review-pr` — an earlier version
//! of this paragraph got both of its claims wrong.** It said `Rung` "has
//! to be" `pub` because `crate::naming::IdentityWitness::new` takes one as a
//! parameter — false: `new` is `pub(crate)`, so only in-crate callers ever
//! need to name the type, and this module (and the type) are `pub(crate)`
//! now. It also said "the trap is a future `Universe.rung` field, not
//! `Rung` itself" — also false, and demonstrated so in review:
//! `Rung::draw` is a pure function of the seed, so `Rung::draw(seed).is_identity()`
//! already recovered the identity predicate from `borbax-molecule` with no
//! field and no witness, while `Rung` was still `pub`. Narrowing this module
//! closes that route today, not merely at some future field. Now that
//! `Universe.rung` exists (`pub(crate)`, per this module's own guidance
//! below), it is held no more visible than `generate_elements_v2` and the
//! bond generator's V2 path — the code that actually needs to read it —
//! require.
//!
//! **This is narrow-API hygiene (CLAUDE.md's idiom tier), not a G2 fix, and
//! the distinction is worth keeping straight.** No fiction guarantee is at
//! stake either way: `Universe::seed` is already `pub`, so the identity
//! predicate is already reconstructible today via `borbax_rng::Domain::PerturbationRung`
//! regardless of this module's own visibility, and G2's 2026-08-07 revision
//! does not treat the identity configuration as a secret — narrowing here is
//! free (nothing outside this crate names `Rung`, `Direction` or `perturb`),
//! not load-bearing. Filing it as a G2 fix would imply the identity
//! configuration needs protecting, which contradicts the ruling that
//! unblocked issue #26 in the first place.
//!
//! **A blanket `#[expect(dead_code)]` still covers this module — narrower
//! again than before, not gone.** `Rung::new`/`draw`/`is_identity`/`strength`
//! and `Direction::get` are all live now, reached transitively from
//! `Universe::generate_v1`/`generate_v2` through
//! [`perturb`](crate::perturbation::perturb) and
//! [`draw_symmetric`](crate::perturbation::draw_symmetric)
//! — an earlier version of this doc listed them as still dead, which stopped
//! being true the moment those generators became real callers. What
//! genuinely remains dead outside this module's own tests:
//! `Rung::all_off_identity` (only `naming.rs`'s own exhaustive test walks
//! it) and `MigratedConstant::ALL`/`.axis()`/`.counterpart()` (nothing in
//! production branches on a constant's classification tags yet — Task
//! 26.2's independence test is the first thing that will). One blanket
//! `#[expect]` rather than four individually-reasoned ones, because the
//! reason is one fact ("these four are read only by tests today"), not
//! four — `cfg_attr`-gated to `not(test)` for the same reason every other
//! such attribute in this crate is: this module's own tests are callers,
//! and an unconditional `#[expect(dead_code)]` would fire
//! `unfulfilled_lint_expectations` on the test target.
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "Rung::all_off_identity and MigratedConstant::ALL/.axis()/.counterpart() have \
                  no in-crate caller outside this module's own tests — see this module's own \
                  doc for what used to be on this list and is now genuinely live"
    )
)]

use borbax_rng::{Domain, Stream};

/// Which physics axis a migrated constant belongs to (issue #26, P7's
/// round-4 correction) — electronic constants are Task 26.1's, nuclear
/// constants are Task 26.2's. `orbital.rs`'s independence from Task 26.2's
/// nuclear model depends on this being right: a constant tagged on the
/// wrong axis would let one task's perturbation silently affect the other's
/// physics.
///
/// **`Nuclear` is tagged by Task 26.2's coefficients (`eps`, `kappa`, `c`
/// — originally four including `sigma`, pre-registered 2026-08-12, `sigma`
/// deleted 2026-08-13 as structurally unjustified) — see
/// [`MigratedConstant`]'s own doc.** Still `#[allow(dead_code)]`, not
/// `#[expect(dead_code)]`: nothing in production branches on `.axis()`'s
/// classification yet (Step 3's independence test is what gives it a first
/// real reader, and Step 3 has not landed), so the variant is constructed
/// only by the generated `axis()` match arms below and read only by this
/// module's own tests — the same "flips between fulfilled and unfulfilled
/// across otherwise-identical clean rebuilds" property that ruled out
/// `#[expect]` for this enum originally still holds, since what changed is
/// which constants tag `Nuclear`, not whether anything outside a test reads
/// the tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Axis {
    Electronic,
    #[allow(
        dead_code,
        reason = "no production caller of .axis() exists yet -- Step 3's independence test is \
                  the first; see this enum's own doc for why #[allow] rather than #[expect]"
    )]
    Nuclear,
}

/// Whether a migrated constant's base value is pinned to its real physical
/// counterpart, or is a chosen constant with no real analogue to measure
/// against (issue #26, P7's round-4 correction, `signature.rs:141-144`'s
/// existing warning about weights with no physical counterpart).
///
/// **`HasReal` is tagged by Task 26.2's `kappa`/`c`** (`base_kappa = 23.7`,
/// `base_c = 0.015`, Rohlf's self-consistent semi-empirical mass formula
/// pair — see [`MigratedConstant`]'s own doc for the citation), pinned
/// 2026-08-12. See [`Axis`]'s own doc for why this is `#[allow]`, not
/// `#[expect]` — the same "no production reader yet" reasoning applies here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Counterpart {
    #[allow(
        dead_code,
        reason = "no production caller of .counterpart() exists yet -- Step 3's independence \
                  test is the first; see Axis's own doc for why #[allow] rather than #[expect]"
    )]
    HasReal,
    DefaultOnly,
}

/// One migrated constant's permanently-assigned stream index, plus its four
/// independent classification tags — generated from one macro-driven table
/// rather than four prose lists, per issue #26 P7's round-4 correction: a
/// constant added later to only some of the tags is a defect nothing
/// catches when the tags are separate lists, and is not representable at
/// all when they are arms of one exhaustive match on one enum.
///
/// **Grows one variant at a time, discriminants never renumbered — the same
/// discipline [`Domain`]'s own variants follow.**
macro_rules! migrated_constants {
    ($($(#[$doc:meta])* $variant:ident = $index:literal, $axis:expr, $counterpart:expr, $p_bound:expr;)+) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        #[repr(u64)]
        pub(crate) enum MigratedConstant {
            $($(#[$doc])* $variant = $index,)+
        }

        impl MigratedConstant {
            /// Every variant, in discriminant order — generated from the same
            /// token list that defines the enum itself, so (unlike
            /// `Domain::ALL` or `PhysicsVersion::ALL`, both hand-written) a
            /// variant cannot be added to the enum without also being added
            /// here: there is only one list to edit.
            pub(crate) const ALL: &'static [Self] = &[$(Self::$variant,)+];

            /// This constant's permanent index within
            /// `Domain::Perturbation` — see [`Domain::Perturbation`]'s own
            /// doc for why a shared `Domain` rather than one per constant.
            #[expect(
                clippy::as_conversions,
                reason = "the enum is #[repr(u64)] with explicit small literal discriminants, \
                          matching Domain::discriminant's own justification"
            )]
            pub(crate) const fn index(self) -> u64 {
                self as u64
            }

            pub(crate) const fn axis(self) -> Axis {
                match self { $(Self::$variant => $axis,)+ }
            }

            pub(crate) const fn counterpart(self) -> Counterpart {
                match self { $(Self::$variant => $counterpart,)+ }
            }

            /// The legal `|p|` excursion bound for this constant's symmetric
            /// draw (used by [`draw_symmetric`]; `ScreeningInner`'s bespoke
            /// one-sided draw does not read this). **Fourth tag, added when
            /// `w_shape`/`w_charge` gave it a real consumer** — every other
            /// variant uses the crate-wide default, [`Direction::P_MAX`].
            pub(crate) const fn p_bound(self) -> f64 {
                match self { $(Self::$variant => $p_bound,)+ }
            }
        }
    };
}

migrated_constants! {
    /// The screened-hydrogenic fold's inner-shell screening coefficient
    /// (`orbital.rs`, Decision 3). **Drawn one-sided downward — `p` is
    /// sampled from `[-Direction::P_MAX, 0.0]`, never the symmetric range
    /// most migrated constants use — at the draw call site in `orbital.rs`,
    /// not read from [`MigratedConstant::p_bound`].** The plan states this
    /// as a specific, one-off decision about this one constant ("draw the
    /// screening coefficients through P7's mechanism, one-sided downward on
    /// the inner-shell coefficient"), not a general dimension every
    /// migrated constant needs an opinion on. `p_bound` is still set here
    /// (to the crate-wide default) so the fourth tag stays total even for a
    /// variant whose own draw site ignores it.
    ScreeningInner = 0, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// Task 26.1 Steps 6-9's valence-promotion energy budget — `orbital.rs`'s
    /// [`crate::orbital::GapConsts`], not named in the plan's own P7
    /// migration table (same gap `ScreeningInner` was in until it gained a
    /// real consumer). A closed-subshell element promotes (gains a nonzero
    /// valence) iff its own HOMO-LUMO gap sits at or below this budget —
    /// see `orbital.rs`'s own doc for the numerical calibration that
    /// separates real alkaline-earth-like gaps (~1.1-1.6) from real
    /// noble-gas-like ones (~3.9+) by more than 2x. **No real physical
    /// counterpart** — this model's own promotion threshold, not a
    /// measured ionisation-energy cutoff, same tagging as `w_shape`/
    /// `w_charge`. Symmetric `p` bound (no one-sided restriction, unlike
    /// `ScreeningInner`): there is no structural reason to favour a
    /// downward-only excursion here.
    PromotionBudget = 1, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// Task 26.1 Steps 6-9's radius scale — `element.rs`'s V2 `radius(z) =
    /// radius_scale * outer_n^2 / Z_eff(outer_n, l=0)`, the Bohr-radius-like
    /// formula `orbital.rs`'s `ElectronicProperties::zeff_outer` feeds.
    /// **No real physical counterpart tagged** — a real Bohr radius
    /// (`a_0 ≈ 0.529` Å) exists, but this crate has no length scale to
    /// compare it against (G4: `radius` is `Span`, an invented unit with no
    /// real-world correspondence pinned anywhere else in this crate), so
    /// pinning `base` to `a_0`'s bare number would be a unit-free numerology
    /// match, not the kind of real correspondence `contact_defect`/
    /// `base_mass` have. Symmetric `p` bound, same reasoning as
    /// `PromotionBudget`.
    RadiusScale = 2, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// Task 26.1 Steps 6-9's bond-capacity scale — `bonds.rs`'s V2 path,
    /// `capacity(el) = capacity_scale + valence(el)`, replacing V1's
    /// packing-derived `contact_density`. **No real physical counterpart**
    /// — this model's own bonding-capacity floor, not a measured quantity,
    /// same tagging as `w_shape`/`w_charge`. Symmetric `p` bound.
    CapacityScale = 3, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// Task 26.1 Steps 6-9's V2 base-mass constant — `element.rs`'s
    /// `mass = units * base_mass_sub - defect_sub`, the same formula shape
    /// V1 uses with V2's own `paired_count` standing in for V1's
    /// `contacts`. **Tagged `DefaultOnly` for now, though the plan's own P7
    /// migration table calls this constant "has a real counterpart in
    /// principle" — this document does not pin a real value for it, and
    /// the plan's own text says to tag it `default only` explicitly rather
    /// than leave it silently untagged in that case.** Re-tag `HasReal`
    /// only once a real value is actually pinned. Re-quantised to the
    /// `1/1024` grid after perturbation — see `element.rs`'s own doc on why
    /// `base_mass` must stay dyadic.
    BaseMass = 4, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// Task 26.1 Steps 6-9's V2 contact-defect constant — same caveat as
    /// [`Self::BaseMass`]: the plan calls this "has a real counterpart in
    /// principle" but pins no value, so it stays `DefaultOnly` until one is.
    ContactDefect = 5, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// `UniverseConsts::w_shape`, §8.3's shape-complementarity weight.
    /// **The plan's own reason this pair (see [`Self::WCharge`]) gets the
    /// tighter `1/3` bound rather than the crate-wide default**: the
    /// prefilter's relative-gap function is a monotone Möbius transform of
    /// the shape/charge ratio `w_shape / w_charge`, verified over 200,000
    /// random draws with 0 monotonicity violations — a bound degenerating
    /// at its extremes needs its extremes kept away from, on pure
    /// numerical-soundness grounds, independent of what chemistry results.
    /// `|p| <= 1/3` gives a reachable ratio range of `[0.5, 2.0]` — tighter
    /// than V1's own independent-draw ratio (`[3/7, 7/3]`), not a
    /// reproduction of it. **No real physical counterpart** — a weight in
    /// this crate's own shape-complementarity kernel, not a measured
    /// quantity.
    WShape = 6, Axis::Electronic, Counterpart::DefaultOnly, 1.0 / 3.0;
    /// `UniverseConsts::w_charge` — see [`Self::WShape`] for the shared
    /// derivation of the `1/3` bound; only the ratio the prefilter reads is
    /// bounded, so both constants need the same tighter excursion.
    WCharge = 7, Axis::Electronic, Counterpart::DefaultOnly, 1.0 / 3.0;
    /// `UniverseConsts::ideal_gap` — **migrated for stream-layout
    /// uniformity, not because it is used.** Zero production consumers
    /// anywhere in the workspace (verified by search, matching the plan's
    /// own note). Migrating it costs one more permanently-assigned index
    /// and keeps every `UniverseConsts` field on one mechanism rather than
    /// one migrated and one not — not because anything downstream needs it
    /// identity-reachable. **No real physical counterpart** — G1 already
    /// flags its value as "on the wrong length scale", so nothing here
    /// changes that.
    IdealGap = 8, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// `bonds.rs`'s energy-scale multiple of `eps` (`BondEnergyMatrix`'s V2
    /// path) — V1's own `SCALE_RANGE = (28.0, 84.0)`, replaced by a
    /// perturbed base. `BASE_SCALE = 56.0` (V1's own range midpoint) at the
    /// crate-wide default `P_MAX = 0.5` reproduces exactly V1's own reachable
    /// range at the perturbed extremes (`56*(1-0.5) = 28`, `56*(1+0.5) = 84`)
    /// — a deliberate match, not a coincidence the constant depends on.
    /// **No real physical counterpart** — a chosen energy-scale multiplier,
    /// same tagging as `w_shape`/`w_charge`.
    Scale = 9, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// `bonds.rs`'s bond-order exponent `gamma` (`OrderScale`) — V1's own
    /// `GAMMA_RANGE = (0.68, 1.0)`. **`BASE_GAMMA = 0.6`, deliberately below
    /// V1's own range, not its midpoint** — `gamma <= 1` is "CHOSEN, not
    /// derived, and load-bearing" (this file's own header: it is what makes
    /// the bond-order series diminishing), and `OrderScale::new` performs
    /// **no clamping or validation** on its input, so a base chosen without
    /// checking the perturbed range against that ceiling could silently
    /// produce a non-diminishing series. At the crate-wide default
    /// `P_MAX = 0.5`, `0.6` keeps every reachable value in `[0.3, 0.9]`, a
    /// comfortable margin inside `(0, 1]` at both ends — not `0.84` (V1's
    /// own range midpoint), which would reach `1.26` at the worst-case
    /// corner and silently invert the series. **No real physical
    /// counterpart** — chosen to keep the series diminishing, not measured.
    Gamma = 10, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// `bonds.rs`'s ionic-excess scale (`BondEnergyMatrix`'s V2 path) —
    /// the rank-2 term the plan's scope section and Task 26.1 Step 6 both
    /// call for and V1 has no analogue of. **Pauling-motivated, found by
    /// direct numerical experiment, not stated in the plan**: real bond
    /// energy exceeds the geometric mean of the two homonuclear bond
    /// energies by an amount proportional to `(chi_A - chi_B)^2`
    /// (electronegativity difference squared) — the mechanism Pauling
    /// originally used to *derive* electronegativity from bond energies.
    /// The Borbax analogue adds `-ionic_base * affinity[a] * affinity[b]`
    /// to the existing rank-1 `base * sqrt(ca * cb)` term: a pure
    /// affinity-outer-product addition, which makes the *whole* matrix
    /// exactly rank 2 (verified: 2 nonzero singular values, the rest
    /// exactly 0, over 5 seeds' worth of V2 tables) rather than the
    /// higher-rank `(affinity[a]-affinity[b])^2` form, whose expansion
    /// contains two more rank-1 terms this project has no separate
    /// mechanism for. **Sign matches Principle 1** (`-x*y` is positive
    /// when `x`, `y` are opposite-signed): elements with complementary
    /// affinity bond *more* strongly, same-signed affinity bonds *less*
    /// strongly. Measured: a rank-1-removed residual of the resulting
    /// matrix correlates 0.86-0.89 (Pearson) with `(affinity[a] -
    /// affinity[b])^2` across those same 5 seeds — the discriminator Step
    /// 6 names.
    ///
    /// **`BASE_IONIC_SCALE = 2.0`, chosen against a proven worst-case
    /// corner, not observed over samples.** The term can only drive the
    /// bond energy negative when `affinity[a] * affinity[b]` is at its
    /// most positive (both elements' affinities near-saturated
    /// same-signed) simultaneously with `Scale` and `CapacityScale` both
    /// at their perturbed floor — `28 * 0.25 = 7.0` in eps-relative units
    /// (`Scale`'s base 56.0 and `CapacityScale`'s base 0.5, each at the
    /// crate-wide `P_MAX = 0.5` worst corner). Against that floor,
    /// `IonicScale`'s own worst-case ceiling (`2.0 * 1.5 = 3.0`) leaves a
    /// provable `4.0` margin — not the `1.0` a naive `BASE_IONIC_SCALE =
    /// 4.0` would leave, and not merely "unobserved over the seeds this
    /// was checked against": across seeds `[0, 1, 7, 42, 1000, 999_999]`
    /// and every reachable rung, the actual worst margin measured was
    /// `>16`, corroborating the analytic bound without depending on it.
    /// **No real physical counterpart** — Borbax's own invented ionic
    /// term, not a measured quantity, same tagging as `w_shape`/`w_charge`.
    /// Symmetric `p` bound: no structural reason to favour one direction.
    IonicScale = 11, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// The screened-hydrogenic fold's deep-screener coefficient
    /// (`orbital.rs`) — a ceiling every d/f-block electron's contribution
    /// to any *other* candidate's `near`/`far` screening sum is capped at,
    /// replacing Slater's idealised "everything closer screens at exactly
    /// 1.00" for d/f screeners specifically. Does **not** reach the
    /// same-shell `t(l)` term (`orbital.rs`'s own doc) — a same-shell peer
    /// is never routed through `near`/`far` at all, so "anyone's screening
    /// sum" would overstate the mechanism's reach. **No real physical
    /// counterpart tagged**, same reasoning as `ScreeningInner`: the
    /// *qualitative* principle (a d/f electron's poor nuclear penetration
    /// also makes it a poor screener of everything else, symmetric to it
    /// being *fully* screened itself — established atomic-physics teaching,
    /// not this codebase's invention) is real, but no citable Slater-style
    /// rule set turns it
    /// into a specific coefficient the way `ScreeningInner`'s `0.85` is
    /// Slater's own historical number; see `orbital.rs`'s own doc for the
    /// full accounting and the calibration this constant's base value
    /// comes from. **Drawn one-sided downward, same shape and same reason
    /// as `ScreeningInner`** — `p` sampled from `[-Direction::P_MAX, 0.0]`,
    /// at the draw call site in `orbital.rs`, not read from
    /// [`MigratedConstant::p_bound`]. `p_bound` is still set here (to the
    /// crate-wide default) so the fourth tag stays total.
    ScreeningDeep = 12, Axis::Electronic, Counterpart::DefaultOnly, Direction::P_MAX;
    /// Task 26.2's nuclear binding-energy scale — the model's `eps *
    /// contacts(total)` term. **Composite, tagged `DefaultOnly` rather than
    /// `HasReal` even though it approximates a real quantity**: it folds the
    /// real semi-empirical mass formula's volume *and* surface terms into
    /// one packing-derived contact count, so it has no single clean real
    /// coefficient to pin against — see the plan's own P7 pre-registration
    /// (2026-08-12) for the full accounting.
    ///
    /// **`BASE_EPS = 2.625`, corrected 2026-08-12 during Step 1's own gate —
    /// this entry originally pinned `1.0`, V1's own historical `eps`
    /// midpoint, which sits at an unrelated scale next to
    /// [`Self::Kappa`]/[`Self::C`]'s real MeV-pinned values.** Measured at
    /// the time, at `eps = 1.0` and the shipped `NUCLEAR_K = 10` on the
    /// four-coefficient model then shipped: 305 of 386 totals (79.0%) had
    /// negative on-valley per-unit binding energy — unbound over most of
    /// its domain, not merely peaked low, caught by a
    /// `physics-plausibility-reviewer` pass dispatched specifically because
    /// `kappa`/`c` carry real correspondence and `eps` did not. (An earlier
    /// draft of this correction quoted 288/386, 74.6% — measured against
    /// the superseded `NUCLEAR_K = 12` before that constant's own
    /// correction landed; and 305/386 was itself the four-coefficient
    /// model's figure — after `sigma`'s 2026-08-13 removal, the identical
    /// measurement at `eps = 1.0` on the current three-coefficient model
    /// gives 190/386, 49.2%, unbound over roughly half its domain, not
    /// most of it. Neither figure is stale *in its own historical context*;
    /// see `crates/borbax-universe/src/nuclear.rs`'s
    /// `the_identity_configuration_is_bound_over_all_but_two_totals` for the
    /// figure as a pinned test.) Corrected to `eps = 2*a_V/z`, `a_V = 15.75`
    /// `MeV` (Rohlf's real SEMF volume coefficient, the same citation `kappa`
    /// uses), `z = 12` (`crate::packing::shell_size(k, 1)` at the nuclear
    /// module's own fixed `k = 10`) — matching the real limit `BE/A -> a_V`
    /// as `A -> infinity`, a limit approached slowly and not attained
    /// within `T_MAX` (`eps * contacts(T_MAX)/T_MAX` is `13.19`, not
    /// `15.75`). See `crates/borbax-universe/src/nuclear.rs` for the landed
    /// constant and the full accounting.
    ///
    /// **Does not feed `energy_per_unit`** (V1's and V2's shared bond-pricing
    /// field) — see the plan's Round 5 correction for why that would
    /// re-price every bond by nuclear structure, ~10⁶× too strong.
    /// Symmetric `p` bound: no structural reason to favour one direction,
    /// and `eps` has zero effect on the composition argmax itself (a
    /// per-total additive constant that cancels out of the comparison
    /// across compositions at fixed total). **`eps`'s *direction* is
    /// independently drawn from `kappa`/`c`'s (each reads its own `Stream`
    /// index, per [`draw_symmetric`]'s doc) — "independent" here means
    /// independently *drawn*, not statistically independent: all three
    /// nuclear constants still share one rung, so their *magnitudes* move
    /// together. *Deriving* `eps` from `kappa`/`c` instead — the way
    /// `gamma` is derived from `kappa`/`c` — was considered and
    /// rejected**, not the independent draw itself: real SEMF
    /// volume/asymmetry coefficients are independently fitted parameters
    /// with no derivable relationship (`a_A/a_V` spans 1.348-1.505 across
    /// five published fits against `a_V`'s own ±5.7% spread), so deriving
    /// `eps` from `kappa` would assert a physical constraint real nuclear
    /// physics does not make.
    Eps = 13, Axis::Nuclear, Counterpart::DefaultOnly, Direction::P_MAX;
    // Index 14 is retired, not free: `Sigma` (the model's `sigma *
    // total^(5/3)` size-cost term) was deleted 2026-08-13 (see
    // `crates/borbax-universe/src/nuclear.rs`'s "Post-gate correction" —
    // a `physics-plausibility-reviewer` pass found it structurally
    // unjustified, not merely miscalibrated: `eps * contacts(total)`
    // already reproduces the real SEMF's volume-and-surface behaviour, so
    // a separate size-cost term duplicated a role already filled). This
    // index is never reused — see [`MigratedConstant`]'s own doc,
    // "discriminants never renumbered": reusing 14 for a future constant
    // would make that constant collide with every historical draw made
    // under the old `Sigma` meaning, silently. `Kappa`/`C` keep 15/16
    // unchanged, so every seed's `kappa`/`c` draw is bit-identical to
    // before this removal.
    /// Task 26.2's nuclear asymmetry coefficient — the model's `kappa *
    /// (a-b)^2 / total` term. **`HasReal`, pinned to Rohlf's semi-empirical
    /// mass formula asymmetry coefficient**: `BASE_KAPPA = 23.7` (`MeV`).
    /// Citation: J. W. Rohlf, *Modern Physics from α to Z⁰*, Wiley, New
    /// York (1994), ISBN 0-471-57270-5 — verified 2026-08-12 against an
    /// earlier draft that had silently mixed Rohlf's own `a_C` with a
    /// different fit's `a_sym = 23.2`; see the plan's citation correction
    /// for the full accounting. `kappa` is pure scale for the composition
    /// argmax (cancels out entirely once `gamma` is reparameterised as
    /// `2*kappa*c`, Decision 11) — its value matters for the model's energy
    /// *magnitude*, never for where the valley sits. Symmetric `p` bound.
    Kappa = 15, Axis::Nuclear, Counterpart::HasReal, Direction::P_MAX;
    /// Task 26.2's dimensionless Coulomb-to-asymmetry ratio, `c =
    /// gamma/(2*kappa)` — the one quantity that actually determines where
    /// the composition-argmax valley sits (Decision 11); `gamma` itself is
    /// **derived**, `2*kappa*c`, and is not an independent perturbation
    /// target. **`HasReal`**: `BASE_C = 0.015` exactly (`a_C / (2 *
    /// a_sym) = 0.711 / (2 * 23.7)`, Rohlf's self-consistent pair — see
    /// [`Self::Kappa`]'s own citation). **The excursion bound is the
    /// crate-wide default despite `c` setting a hard structural
    /// requirement** (coverage: no element's isotope distribution may be
    /// empty) **— the requirement is met by choosing `T_MAX`, not by
    /// tightening `p`.** At `T_MAX = 386` (the plan's own pre-registration,
    /// 2026-08-12), the worst-case reachable `c` (`BASE_C * 1.5`, at
    /// `Direction::P_MAX = 0.5`) still gives `a*(T_MAX) = 121.08`, a
    /// 1.5-unit buffer above the `119.5` bare-minimum coverage threshold
    /// for `n_elements = 120` — see the plan's own derivation, independently
    /// verified by `geometry-numerics-reviewer` (0/48,000 mismatches against
    /// brute-force integer argmax). `T_MAX` itself is not migrated through
    /// P7 — it is a fixed design constant of the model family, not a
    /// per-universe or per-constant quantity (Decision 11: using the drawn
    /// table's own size would make `c`'s legal range circular and couple
    /// the nuclear axis to Task 26.1's electronic one through `n_elements`,
    /// contradicting Step 3's independence requirement).
    ///
    /// **Correction, 2026-08-12, found during Step 1's own gate: this proof
    /// covers only the composition argmax's shape (does `a*(t)` reach every
    /// integer as `t` ranges to `T_MAX`), not whether binding stays
    /// positive out that far — and those are different questions once
    /// `eps` is perturbed too.** `contacts(total)` itself is **not**
    /// bounded (it grows roughly linearly in `total`) — what saturates is
    /// `contacts(total)/total`, toward `z/2 = 6`, so `eps *
    /// contacts(total)/total` (the volume term's own per-unit contribution)
    /// is bounded, and the on-valley per-unit energy can in principle turn
    /// negative past some total if `eps`/`kappa` are perturbed unfavourably
    /// — real physics too (real SEMF binding per nucleon crosses zero at
    /// `A ~ 3076` under the `Z(Z-1)` Coulomb convention this model itself
    /// uses, far past `T_MAX`, so the identity configuration correctly
    /// reproduces "bound everywhere in range"), but under perturbation that
    /// crossover can be pulled inside `1..=T_MAX`. `T_MAX`'s coverage proof
    /// is real and still holds for what it actually proves; it does not by
    /// itself guarantee the table stays bound that far. See
    /// `crates/borbax-universe/src/nuclear.rs`'s own gate for the
    /// measurement; how this bears on `n_elements`'s own draw is a question
    /// for Steps 2+, not resolved here.
    ///
    /// **Post-gate correction, 2026-08-13: this entry originally also
    /// carried a fourth coefficient, `sigma`, and roughly half this doc's
    /// own history discussed `sigma`'s outsized contribution to the
    /// shortfall above.** A `physics-plausibility-reviewer` pass found
    /// `sigma * total^(5/3)` structurally unjustified rather than merely
    /// miscalibrated — `eps * contacts(total)` already reproduces the real
    /// SEMF's volume-plus-surface term (worst-case error 2.63% at
    /// `total = 78`, mean absolute error 1.20%, within 1% from
    /// `total >= 266` onward), so the three surviving coefficients already
    /// form the real four-term SEMF at Rohlf's own values. `sigma` was
    /// deleted, along with this variant's own `Sigma` sibling
    /// (`MigratedConstant`'s index 14, now retired) — see
    /// `crates/borbax-universe/src/nuclear.rs`'s own module doc, "Post-gate
    /// correction," for the full account.
    ///
    /// **Measured over the same joint P7 sample, after the removal: the
    /// bound range — the last total any composition is still bound at, per
    /// seed — reaches p50 = 138, p10 = 130, p1 = 116, min = 38 elements,
    /// with 1.09% (419/38,416) of universes covering fewer than 120
    /// elements while still bound, down from 9.25% (3,554/38,416) at the
    /// four-coefficient model.** `composition_argmax` itself reaches every
    /// integer up to 120 given a large enough `total`, so this residual
    /// shortfall — smaller, not zero — is still a real, separate finding
    /// for Steps 2+, not something Step 1 itself is positioned to fix by
    /// picking a threshold here. Removing `sigma` was justified by the
    /// structural finding above, not by this improvement; the improvement
    /// is its consequence, not its reason.
    C = 16, Axis::Nuclear, Counterpart::HasReal, Direction::P_MAX;
}

/// How far a universe's constants sit from the identity configuration —
/// `0..=MAX`, drawn once per universe.
///
/// **`MAX = 99`, settling routed requirement 8 — a product decision, not a
/// tuning parameter, made explicitly rather than left implicit in `MAX`'s
/// value.** Ian, 2026-08-12, choosing `P(identity) = 1/(MAX + 1) = 1%`: the
/// search-cost case for a *common* identity seed is weak, since it only
/// needs to be found once and then ships as a named universe (Step 1b's own
/// precedent — typing a name like `emily` already reaches a specific seed
/// directly), so real chemistry stays a genuine rarity rather than an
/// ordinary rung a casual seed sweep turns up roughly one time in nine, as
/// the placeholder `MAX = 8` this replaces did. `MAX = 99` also gives every
/// migrated constant — Task 26.1's electronic axis included, not only Task
/// 26.2's new nuclear one — a finer 100-rung perturbation ladder than the
/// placeholder's 9 rungs. Nothing here depends on the exact number beyond
/// what `the_binomial_band_pins_p_identity_as_a_number` (this module's own
/// `#[cfg(test)]` tests, not linkable from public docs) checks against
/// whatever `MAX` currently is — but note that test's own corpus size and
/// band were recalibrated for this value, not merely left as they were: at
/// `P(identity) = 1%`, the gap to the `next_range(2·MAX + 1)` defect
/// reading (`1/199 ≈ 0.503%`) is under half a percentage point, an order of
/// magnitude tighter than the placeholder's `1/9` vs `1/17` gap, so the
/// same corpus size and band that safely discriminated the two at `MAX = 8`
/// do not at `MAX = 99` — see that test's own comment for the recomputed
/// figures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Rung(u8);

impl Rung {
    /// The top of the rung's legal range, `m` in the plan's notation. See
    /// this type's own doc for the 2026-08-12 pre-registration that pinned
    /// this value.
    pub(crate) const MAX: u8 = 99;

    /// `Rung(0)` — the identity configuration.
    pub(crate) const IDENTITY: Self = Self(0);

    /// A rung, or `None` outside `0..=MAX`.
    ///
    /// **Fallible, not a bare `u8` field, for the reason `BondOrder::new` and
    /// `Direction::new` both are: an unbounded rung reopens the excursion
    /// hazard `Direction`'s own bound exists to close.** At `rung > MAX` the
    /// perturbation factor's infimum is `+0.0` exactly once `rung = 2*MAX`,
    /// and past `4*MAX` a quarter of `p`'s range gives a **negative**
    /// factor — see [`perturb`]'s doc for the proof this constructor makes
    /// unreachable.
    #[must_use]
    pub(crate) const fn new(r: u8) -> Option<Self> {
        if r <= Self::MAX { Some(Self(r)) } else { None }
    }

    /// This universe's rung, from its own stream — never a shared index.
    ///
    /// **`next_range(MAX + 1)`, not `next_range(2*MAX + 1)`.** The latter
    /// gives `2*MAX + 1` outcomes and halves `P(identity)` silently; the
    /// former is the spelling that actually delivers `P(identity) =
    /// 1/(MAX + 1)`, pinned by `the_binomial_band_pins_p_identity_as_a_number`
    /// (this module's own tests) as a number, not an order-of-magnitude
    /// check that both readings would satisfy.
    ///
    /// Draws from `Stream::new(seed, Domain::PerturbationRung, 0)` — its own
    /// `Domain`, never an index shared with anything else. See
    /// `Domain::PerturbationRung`'s own doc for why a shared index is a live
    /// hazard under a discrete `p` ladder, not merely a style preference.
    ///
    /// **Routes through [`Self::new`], found bypassing it in `/review-pr`.**
    /// An earlier version built `Self(r as u8)` directly — the only
    /// production producer of a `Rung` skipping the one place `r <= MAX` is
    /// enforced. Under the exact defect named above, that bypass let an
    /// out-of-range `Rung` reach [`perturb`]'s annihilation proof, which
    /// assumes every `Rung` was checked. `new`'s own fallible check makes
    /// this the same shape of fix as `Direction::new`'s `bound <= P_MAX`
    /// clause: the invariant enforced at the one place a caller could skip
    /// it, not merely at every site that happens to construct one today.
    #[must_use]
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        reason = "next_range(u64::from(Self::MAX) + 1) returns a value in 0..=MAX, and MAX is a \
                  u8 literal, so the result fits u8 by construction — proven by the domain of \
                  next_range's input, not merely assumed"
    )]
    pub(crate) fn draw(seed: u64) -> Self {
        let mut stream = Stream::new(seed, Domain::PerturbationRung, 0);
        let r = stream.next_range(u64::from(Self::MAX) + 1);
        Self::new(r as u8).unwrap_or_else(|| unreachable!("next_range(MAX + 1) cannot exceed MAX"))
    }

    /// Every rung except the identity, `1..=MAX` — for Decision 9's
    /// exhaustive test (`IdentityWitness::new(r).is_none()` for every `r`
    /// here), so the corpus is derived from `MAX` rather than a hardcoded
    /// `1..=RUNGS` literal that could drift from it.
    pub(crate) fn all_off_identity() -> impl Iterator<Item = Self> {
        (1..=Self::MAX).map(Self)
    }

    /// `true` for `Rung(0)` — the identity configuration.
    #[must_use]
    pub(crate) const fn is_identity(self) -> bool {
        self.0 == 0
    }

    /// This rung's magnitude, `rung / MAX`, in `[0.0, 1.0]` — derived, never
    /// stored. "Store the integer rung, not the derived `f64` strength"
    /// (Decision 10): the integer is what every test and every digest
    /// should key on, so a strength recomputed here can never disagree with
    /// one stored somewhere else.
    #[must_use]
    pub(crate) fn strength(self) -> f64 {
        f64::from(self.0) / f64::from(Self::MAX)
    }
}

/// A bounded, finite perturbation direction — `p` in the plan's notation.
///
/// **The bound is a constructor parameter, not a shared constant.** A
/// single global `P_MAX` cannot express `w_shape`/`w_charge`'s tighter
/// `1/3` bound (issue #26's plan, P7), so `new` takes `bound` explicitly
/// rather than reading [`Self::P_MAX`] internally — that constant is the
/// *default* most constants use, not the only legal value.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct Direction(f64);

impl Direction {
    /// The default excursion bound most migrated constants use. `w_shape`
    /// and `w_charge` are tighter, at `1/3` — see issue #26's plan, P7, for
    /// the numerical-soundness argument (a monotone Möbius transform of
    /// their ratio, degenerating at its extremes) that bound rests on.
    pub(crate) const P_MAX: f64 = 0.5;

    /// A direction, or `None` if `v` is non-finite, `bound` exceeds
    /// [`Self::P_MAX`], or `v` exceeds `bound` in magnitude.
    ///
    /// **Finite-and-bounded by construction, not by a `debug_assert!`** —
    /// which is absent from the `--release` profile `goldens --emit` runs
    /// in, and this workspace has already shipped that exact debug/release
    /// asymmetry twice (`borbax-units`, Task 2). `borbax_rng`'s own
    /// `next_f64_range` documents swallowing a non-finite bound (`lo =
    /// -inf` with finite `hi` returns `-inf` on every draw), so an infinite
    /// `p` is producible by a typo'd range and nothing upstream of this
    /// constructor rejects it.
    ///
    /// **`bound` is itself checked against `P_MAX`, found missing in
    /// review.** [`perturb`]'s "annihilation is unreachable" proof holds
    /// "for any `Direction` constructed with `bound <= 0.5`" — a caller
    /// obligation the earlier version of this constructor never enforced,
    /// so `Direction::new(-2.0, 5.0)` would have succeeded and, fed through
    /// `perturb` at `Rung::MAX`, produced a negative multiplier the proof
    /// claims impossible. Checking `bound` here makes it a structural
    /// guarantee instead, at the one call site rather than at every future
    /// one Task 26.1 adds.
    ///
    /// **`-0.0` is canonicalised to `+0.0`, per Decision 10's own text**
    /// ("canonicalise `-0.0` in `p` before any digest hashes the
    /// perturbation vector"). Adding `+0.0` rather than comparing against
    /// zero: IEEE-754 addition leaves every finite non-zero value unchanged
    /// and turns negative zero into positive zero under the default
    /// rounding mode, so this moves exactly one bit pattern and needs no
    /// `float_cmp` suppression.
    ///
    /// **Not observable through [`perturb`] at `Rung::IDENTITY`, which is
    /// why this constructor is tested directly rather than through it —
    /// found claimed the other way round in `/review-pr`.** At
    /// `strength() == 0.0`, `0.0 * -0.0` and `0.0 * 0.0` are both `0.0`,
    /// `1.0 + 0.0` is exactly `1.0` regardless of which zero sign fed it,
    /// and `base * 1.0` preserves `base` bit for bit — so the *multiplier*
    /// is `1.0` either way and `perturb`'s result cannot distinguish the two
    /// input signs at that rung. An earlier version of this module claimed
    /// the opposite, that routing `-0.0` through `perturb` "verified" the
    /// canonicalisation; it verified only that identity-rung perturbation
    /// is a no-op, which was never in question. See `a_negative_zero_direction_is_canonicalised`
    /// (this module's own tests) for the assertion that actually observes
    /// this constructor's output.
    #[must_use]
    pub(crate) fn new(v: f64, bound: f64) -> Option<Self> {
        if v.is_finite() && bound <= Self::P_MAX && v.abs() <= bound {
            Some(Self(v + 0.0))
        } else {
            None
        }
    }

    /// The raw value.
    #[must_use]
    pub(crate) const fn get(self) -> f64 {
        self.0
    }
}

/// Apply Decision 10's mechanism: `base * (1.0 + strength * p)`, never the
/// additive form.
///
/// **Proof that annihilation (`base * 0.0`) is unreachable, restated from
/// the plan rather than merely cited.** `1.0 + strength * p` is exactly
/// `0.0` iff `strength * p == -1.0`. `rung.strength()` is in `[0.0, 1.0]` by
/// [`Rung::new`]'s domain, and every [`Direction`] is constructed with
/// `bound <= 0.5` — enforced by [`Direction::new`] itself, not merely a
/// caller convention, since a review found the earlier version of that
/// constructor left `bound` unchecked — so `strength * p` is in `[-0.5,
/// 0.5]` and `1.0 + strength * p` is in `[0.5, 1.5]`, which never reaches
/// `0.0` — proven at the worst-case corner by
/// `the_multiplier_never_reaches_zero_at_the_worst_case_corner` (this
/// module's own tests), not merely asserted.
///
/// **The multiplier staying in `[0.5, 1.5]` does not by itself mean `base *
/// multiplier` cannot underflow to `0.0` — found stated too broadly in
/// `/review-pr`.** For any `base` above the subnormal range this is moot,
/// but at `base` in the smallest few subnormals (`f64::from_bits(1)` being
/// the extreme case) `base * 0.5` itself underflows to `0.0` — not a live
/// hazard for anything this mechanism migrates today (every constant listed
/// for Task 26.1's migration is `O(1)` or larger), but worth stating
/// precisely rather than as an unqualified "unreachable".
///
/// **The expression shape is pinned, not merely preferred.** `base * (1.0 +
/// strength * p)` and the algebraically-equal `base + base * strength * p`
/// are different `f64` values — floating-point multiplication does not
/// distribute over addition — so redistributing this after Task 26.1 wires
/// the first migrated constant moves every golden downstream of it, the
/// same class of silent physics change §13.1 exists to catch. Commuting
/// either multiply (`p * strength`, or `(1.0 + ..) * base`) is bit-identical
/// and therefore safe; distributing is not.
#[must_use]
pub(crate) fn perturb(base: f64, rung: Rung, p: Direction) -> f64 {
    base * (1.0 + rung.strength() * p.get())
}

/// Draw and apply a migrated constant's symmetric excursion — `p` sampled
/// from `[-constant.p_bound(), constant.p_bound()]`, the shape every
/// migrated constant except [`MigratedConstant::ScreeningInner`] uses.
///
/// **Extracted once this shape reached its fifth call site** (`orbital.rs`'s
/// `GapConsts::draw`, `bonds.rs`'s V2 capacity draw, and `element.rs`'s
/// `radius_scale`/`base_mass`/`contact_defect` draws), so there is one
/// `Stream::new(seed, Domain::Perturbation, index)` → `Direction::new` →
/// [`perturb`] chain to keep correct rather than five copies that could
/// drift. `ScreeningInner`'s one-sided-downward draw stays bespoke at its
/// own call site — see [`MigratedConstant`]'s own doc for why that is a
/// draw-site decision, not a general shape.
///
/// **Reads the bound from `constant` itself, not a hardcoded
/// `Direction::P_MAX`** — added when `w_shape`/`w_charge` gave `p_bound`
/// its first non-default value.
#[must_use]
pub(crate) fn draw_symmetric(seed: u64, rung: Rung, constant: MigratedConstant, base: f64) -> f64 {
    let bound = constant.p_bound();
    let mut stream = Stream::new(seed, Domain::Perturbation, constant.index());
    let p = stream.next_f64_range(-bound, bound);
    let direction = Direction::new(p, bound)
        .unwrap_or_else(|| unreachable!("p is drawn within [-bound, bound] by construction"));
    perturb(base, rung, direction)
}

#[cfg(test)]
mod tests {
    use super::{Axis, Counterpart, Direction, MigratedConstant, Rung, draw_symmetric, perturb};

    #[test]
    fn migrated_constant_all_contains_every_variant_with_its_own_index() {
        // Generated from the same macro token list that defines the enum
        // itself (see MigratedConstant's own doc), so this is really a
        // regression test on the macro's own correctness, not on any one
        // variant.
        //
        // Pins explicit (variant, index) pairs, not "index == position in
        // ALL" -- the latter held only by coincidence while the enum had no
        // gaps, and Sigma's removal (2026-08-13, index 14 retired rather
        // than reused -- see MigratedConstant's own doc) created the first
        // one. Nothing in the workspace consumes a variant's *position*;
        // draw_symmetric and orbital.rs both consume .index() directly, so
        // the real property to pin is "ALL is exactly this list, in
        // strictly increasing index order, with no duplicate index" -- which
        // this is a strictly stronger check of than position-equality was
        // (it still catches a typo'd literal, and it survives every future
        // gap the same way this one had to be accommodated).
        let expected: &[(MigratedConstant, u64)] = &[
            (MigratedConstant::ScreeningInner, 0),
            (MigratedConstant::PromotionBudget, 1),
            (MigratedConstant::RadiusScale, 2),
            (MigratedConstant::CapacityScale, 3),
            (MigratedConstant::BaseMass, 4),
            (MigratedConstant::ContactDefect, 5),
            (MigratedConstant::WShape, 6),
            (MigratedConstant::WCharge, 7),
            (MigratedConstant::IdealGap, 8),
            (MigratedConstant::Scale, 9),
            (MigratedConstant::Gamma, 10),
            (MigratedConstant::IonicScale, 11),
            (MigratedConstant::ScreeningDeep, 12),
            (MigratedConstant::Eps, 13),
            // index 14 retired (formerly Sigma) -- not present, not reused.
            (MigratedConstant::Kappa, 15),
            (MigratedConstant::C, 16),
        ];
        assert_eq!(
            MigratedConstant::ALL.len(),
            expected.len(),
            "ALL's length should match the pinned (variant, index) list"
        );
        for (constant, (expected_variant, expected_index)) in
            MigratedConstant::ALL.iter().zip(expected)
        {
            assert_eq!(
                constant, expected_variant,
                "ALL's variant order should match the pinned list"
            );
            assert_eq!(
                constant.index(),
                *expected_index,
                "{constant:?}: index should match its pinned value"
            );
        }
        assert!(
            !MigratedConstant::ALL.iter().any(|c| c.index() == 14),
            "index 14 is retired (formerly Sigma, deleted 2026-08-13) and must never be reused \
             -- see MigratedConstant's own doc, 'discriminants never renumbered'. A future \
             constant drawn at this index would silently correlate its perturbation with every \
             historical seed's now-deleted sigma draw"
        );
        // Duplicate discriminants are already a compile error (rustc E0081) on this
        // #[repr(u64)] enum with explicit literals -- verified with a standalone repro.
        // Kept as a belt-and-braces read of ALL rather than of the enum, not because this
        // branch is reachable.
        let mut indices: Vec<u64> = MigratedConstant::ALL.iter().map(|c| c.index()).collect();
        indices.sort_unstable();
        indices.dedup();
        assert_eq!(
            indices.len(),
            MigratedConstant::ALL.len(),
            "every variant's index should be unique"
        );
        for pair in MigratedConstant::ALL.windows(2) {
            let [a, b] = pair else {
                unreachable!("windows(2) always yields a 2-element slice")
            };
            assert!(
                a.index() < b.index(),
                "ALL should be in strictly increasing index order: {a:?} ({}) then {b:?} ({})",
                a.index(),
                b.index()
            );
        }
    }

    /// The electronic/nuclear partition, shared by
    /// `electronic_constants_are_electronic_and_default_only`,
    /// `nuclear_constants_are_nuclear_and_only_kappa_and_c_have_real_counterparts`
    /// and `every_migrated_constant_is_electronic_or_nuclear_exactly_once` —
    /// one canonical pair of lists rather than three independent copies, so
    /// a constant added to the wrong half has exactly one place to hide,
    /// not three.
    const ELECTRONIC: &[MigratedConstant] = &[
        MigratedConstant::ScreeningInner,
        MigratedConstant::PromotionBudget,
        MigratedConstant::RadiusScale,
        MigratedConstant::CapacityScale,
        MigratedConstant::BaseMass,
        MigratedConstant::ContactDefect,
        MigratedConstant::WShape,
        MigratedConstant::WCharge,
        MigratedConstant::IdealGap,
        MigratedConstant::Scale,
        MigratedConstant::Gamma,
        MigratedConstant::IonicScale,
        MigratedConstant::ScreeningDeep,
    ];
    const NUCLEAR: &[MigratedConstant] = &[
        MigratedConstant::Eps,
        MigratedConstant::Kappa,
        MigratedConstant::C,
    ];

    /// **Task 26.1's electronic constants stay electronic, and only
    /// `Kappa`/`C` (of Task 26.2's three surviving nuclear additions) carry
    /// `HasReal`.** Superseded 2026-08-12: an earlier version of this test
    /// asserted every variant was `Axis::Electronic` and
    /// `Counterpart::DefaultOnly`, which stopped being true the moment
    /// Task 26.2's nuclear constants landed. Split by an explicit
    /// electronic/nuclear partition, derived from the same list the
    /// `ALL`-contents test pins, rather than iterating `ALL` and special-
    /// casing three variants inline — a constant added later to the wrong
    /// half of this partition is exactly the defect this test exists to
    /// catch, and an inline special case would be one place that defect
    /// could hide.
    #[test]
    fn electronic_constants_are_electronic_and_default_only() {
        for &constant in ELECTRONIC {
            assert_eq!(
                constant.axis(),
                Axis::Electronic,
                "{constant:?}: should be Axis::Electronic"
            );
            assert_eq!(
                constant.counterpart(),
                Counterpart::DefaultOnly,
                "{constant:?}: should be Counterpart::DefaultOnly"
            );
        }
    }

    /// The nuclear side of the same partition — see
    /// [`electronic_constants_are_electronic_and_default_only`] for why
    /// this is split rather than exhaustive over `ALL` with exceptions.
    #[test]
    fn nuclear_constants_are_nuclear_and_only_kappa_and_c_have_real_counterparts() {
        for &constant in NUCLEAR {
            assert_eq!(
                constant.axis(),
                Axis::Nuclear,
                "{constant:?}: should be Axis::Nuclear"
            );
        }
        assert_eq!(
            MigratedConstant::Eps.counterpart(),
            Counterpart::DefaultOnly,
            "Eps should be Counterpart::DefaultOnly (no single real coefficient)"
        );
        for &constant in &[MigratedConstant::Kappa, MigratedConstant::C] {
            assert_eq!(
                constant.counterpart(),
                Counterpart::HasReal,
                "{constant:?}: should be Counterpart::HasReal (Rohlf's self-consistent \
                 a_sym/a_C pair)"
            );
        }
    }

    /// **The partition above is exhaustive, not merely non-overlapping —
    /// found worth pinning separately after the split.** Two lists that
    /// each individually check out could still, together, omit a variant
    /// `ALL` contains (or double-count one) with neither per-list test
    /// noticing.
    ///
    /// **Corrected by `/review-pr` (`CodeRabbit`, corroborating an earlier
    /// informal `emergence-auditor` note): the original version of this
    /// test summed `axis() == Electronic` and `axis() == Nuclear` counts
    /// and compared the sum to `ALL.len()` — a tautology given `Axis` has
    /// exactly two variants, since every constant's `axis()` call returns
    /// one of the two by construction regardless of whether `ELECTRONIC`/
    /// `NUCLEAR` are correct. A third `Axis` variant added later without
    /// updating either list would still pass silently.** Checks per-
    /// constant membership against the two explicit lists instead: every
    /// `ALL` entry must appear in exactly one of `ELECTRONIC`/`NUCLEAR`,
    /// and its `axis()` must agree with which list it's in — the check
    /// this test's name has always claimed to make.
    #[test]
    fn every_migrated_constant_is_electronic_or_nuclear_exactly_once() {
        for &constant in MigratedConstant::ALL {
            let in_electronic = ELECTRONIC.contains(&constant);
            let in_nuclear = NUCLEAR.contains(&constant);
            assert!(
                in_electronic ^ in_nuclear,
                "{constant:?}: must appear in exactly one of ELECTRONIC/NUCLEAR -- found in \
                 ELECTRONIC={in_electronic}, in NUCLEAR={in_nuclear}"
            );
            let expected_axis = if in_nuclear {
                Axis::Nuclear
            } else {
                Axis::Electronic
            };
            assert_eq!(
                constant.axis(),
                expected_axis,
                "{constant:?}: axis() disagrees with which partition list it's in"
            );
        }
        assert_eq!(
            NUCLEAR.len(),
            3,
            "Task 26.2's nuclear model has three coefficients (eps, kappa, c) since the \
             2026-08-13 post-gate removal of sigma -- see nuclear.rs's own module doc"
        );
    }

    #[test]
    fn most_migrated_constants_use_the_default_p_bound() {
        for &constant in MigratedConstant::ALL {
            if matches!(
                constant,
                MigratedConstant::WShape | MigratedConstant::WCharge
            ) {
                continue;
            }
            assert_eq!(
                constant.p_bound().to_bits(),
                Direction::P_MAX.to_bits(),
                "{constant:?}: should use the crate-wide default p_bound"
            );
        }
    }

    /// The one place `p_bound` actually varies — checked directly, not just
    /// via the "everything else is default" sweep above, since a planted
    /// defect setting `WShape`/`WCharge` to the default would still pass
    /// that sweep (it explicitly skips them) but silently reopen exactly
    /// the size-comparison pathology the tighter bound exists to close.
    #[test]
    fn w_shape_and_w_charge_use_the_tighter_third_bound() {
        assert!((MigratedConstant::WShape.p_bound() - 1.0 / 3.0).abs() < f64::EPSILON);
        assert!((MigratedConstant::WCharge.p_bound() - 1.0 / 3.0).abs() < f64::EPSILON);
        assert!(
            MigratedConstant::WShape.p_bound() < Direction::P_MAX,
            "the whole point of this tag is that it's tighter than the default"
        );
    }

    /// `Gamma`'s base (0.6) was chosen specifically to keep the perturbed
    /// range inside `OrderScale`'s unvalidated `(0, 1]` requirement — this
    /// test is the actual guard against a future edit widening `BASE_GAMMA`
    /// (in `bonds.rs`'s `generate_v2`) back toward V1's own range midpoint
    /// without re-checking the worst-case corner, since `OrderScale::new`
    /// itself performs no clamping and would silently accept a broken
    /// value.
    #[test]
    fn gamma_stays_a_diminishing_series_at_every_reachable_rung() {
        const BASE_GAMMA: f64 = 0.6;
        for seed in [0_u64, 7, 42, 1000] {
            for r in 0..=Rung::MAX {
                let rung = Rung::new(r).unwrap_or_else(|| unreachable!("r <= Rung::MAX"));
                let gamma = draw_symmetric(seed, rung, MigratedConstant::Gamma, BASE_GAMMA);
                assert!(
                    gamma > 0.0 && gamma <= 1.0,
                    "seed {seed} rung {r}: gamma {gamma} left (0, 1], OrderScale::new performs \
                     no clamping so this would silently invert the bond-order series"
                );
            }
        }
    }

    #[test]
    fn draw_symmetric_round_trips_at_identity() {
        // At rung == 0, every migrated constant lands on exactly its base
        // value, regardless of which constant or which seed -- the same
        // property `OrbitalConsts::draw`/`GapConsts::draw`'s own round-trip
        // tests check for their bespoke draws.
        for seed in [0, 7, 42] {
            for constant in MigratedConstant::ALL {
                let base = 1.375;
                assert_eq!(
                    draw_symmetric(seed, Rung::IDENTITY, *constant, base).to_bits(),
                    base.to_bits(),
                    "seed {seed}, {constant:?}: should equal base exactly at rung == 0"
                );
            }
        }
    }

    #[test]
    fn a_rung_past_max_is_rejected() {
        assert!(Rung::new(Rung::MAX).is_some());
        assert!(Rung::new(Rung::MAX + 1).is_none());
        assert!(Rung::new(u8::MAX).is_none());
    }

    #[test]
    fn all_off_identity_is_exactly_one_through_max() {
        // `r.0`: `mod tests` is a child of `perturbation`, so the private
        // field is visible here — a round-trip through `strength()` would
        // test the wrong thing (that `strength` inverts cleanly), not that
        // `all_off_identity` yields the right integers.
        let got: Vec<u8> = Rung::all_off_identity().map(|r| r.0).collect();
        let want: Vec<u8> = (1..=Rung::MAX).collect();
        assert_eq!(got, want);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "CLAUDE.md: tests may assert exactly. 0/MAX and MAX/MAX are exact under \
                  IEEE-754 division (numerator equals denominator or is zero), so an epsilon \
                  would defeat the point of asserting on it"
    )]
    fn identity_has_zero_strength_and_every_other_rung_is_positive() {
        assert_eq!(Rung::IDENTITY.strength(), 0.0);
        assert!(Rung::IDENTITY.is_identity());
        for r in Rung::all_off_identity() {
            assert!(r.strength() > 0.0, "rung {r:?} has non-positive strength");
            assert!(!r.is_identity());
        }
        let max_rung = Rung::new(Rung::MAX)
            .unwrap_or_else(|| unreachable!("Rung::MAX is in 0..=Rung::MAX by construction"));
        assert_eq!(max_rung.strength(), 1.0);
    }

    #[test]
    fn a_direction_rejects_non_finite_and_out_of_bound_values() {
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                Direction::new(v, Direction::P_MAX).is_none(),
                "{v} should be rejected"
            );
        }
        assert!(Direction::new(0.6, Direction::P_MAX).is_none());
        assert!(Direction::new(-0.6, Direction::P_MAX).is_none());
        // The boundary itself is legal: `<=`, not `<`.
        assert!(Direction::new(Direction::P_MAX, Direction::P_MAX).is_some());
        assert!(Direction::new(-Direction::P_MAX, Direction::P_MAX).is_some());
        assert!(Direction::new(0.0, Direction::P_MAX).is_some());
    }

    /// **The `bound <= P_MAX` clause itself, found untested in `/review-pr`.**
    /// Deleting it left all 96 `borbax-universe` tests green — the clause is
    /// what makes [`perturb`]'s "annihilation is unreachable" proof
    /// structural rather than a caller convention. Without it,
    /// `Direction::new(-2.0, 5.0)` succeeds, and at `Rung::MAX` `perturb`'s
    /// multiplier is exactly `-1.0` — the sign flip the proof claims
    /// impossible for every migrated constant, once Task 26.1 wires a real
    /// call site.
    #[test]
    fn a_bound_wider_than_p_max_is_rejected() {
        assert!(Direction::new(-2.0, 5.0).is_none());
        assert!(Direction::new(0.0, f64::INFINITY).is_none());
        assert!(Direction::new(0.0, f64::NAN).is_none());
        // The property, not just the arm: every constructible pair keeps
        // perturb's multiplier inside [0.5, 1.5].
        for r in (0..=Rung::MAX).filter_map(Rung::new) {
            for raw in [-0.5, -0.25, 0.0, 0.5, -2.0, 3.0, 5.0] {
                for bound in [0.5, 1.0 / 3.0, 0.0, 1.0, 5.0] {
                    if let Some(d) = Direction::new(raw, bound) {
                        let m = perturb(1.0, r, d);
                        assert!((0.5..=1.5).contains(&m), "multiplier {m} escaped");
                    }
                }
            }
        }
    }

    /// Decision 10's own required test: at `rung == 0`, every base value is
    /// bit-identical, over the full finite range of `p` at both signs —
    /// **compared by `to_bits()`, not by value**, because an earlier
    /// version of this test used `assert_eq!` on the `f64`s directly, and
    /// `assert_eq!(-0.0, 0.0)` passes: the one input in the `base` list
    /// chosen to exercise sign-of-zero preservation was exactly the input
    /// the comparison could not see a mismatch in.
    ///
    /// **`-0.0` is in the `p` list for parity with Decision 10's text ("at
    /// both zero signs of `p`"), not because this test can see
    /// `Direction::new`'s canonicalisation — corrected in `/review-pr`.**
    /// An earlier version of this doc claimed the `-0.0` entry "now
    /// exercises that canonicalisation"; it cannot, by construction, at
    /// `Rung::IDENTITY` — see [`Direction::new`]'s own doc for why.
    /// `a_negative_zero_direction_is_canonicalised` (below) is the test
    /// that actually observes it.
    #[test]
    fn rung_zero_leaves_every_base_value_unchanged() {
        for base in [0.0, 1.0, -1.0, 1e-300, 1e300, 42.5, -0.0] {
            for p in [
                0.0,
                -0.0,
                Direction::P_MAX,
                -Direction::P_MAX,
                Direction::P_MAX / 3.0,
                -Direction::P_MAX / 3.0,
            ] {
                let d = Direction::new(p, Direction::P_MAX)
                    .unwrap_or_else(|| unreachable!("p is within [-P_MAX, P_MAX] by construction"));
                assert_eq!(
                    perturb(base, Rung::IDENTITY, d).to_bits(),
                    base.to_bits(),
                    "base {base}, p {p}: identity rung changed the value's bits"
                );
            }
        }
    }

    /// **The assertion `rung_zero_leaves_every_base_value_unchanged` cannot
    /// make, added in `/review-pr`.** That test's `p = -0.0` entry passes
    /// whether or not `Direction::new` canonicalises the sign, because
    /// `Rung::IDENTITY`'s multiplier is `1.0` regardless — see
    /// [`Direction::new`]'s own doc. This asserts on the constructor's
    /// output directly, where the canonicalisation is actually observable.
    #[test]
    fn a_negative_zero_direction_is_canonicalised() {
        let d = Direction::new(-0.0, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("-0.0 is within [-P_MAX, P_MAX] by construction"));
        assert_eq!(
            d.get().to_bits(),
            0.0_f64.to_bits(),
            "-0.0 survived into a Direction uncanonicalised"
        );
    }

    /// The worst-case corner of the excursion-bound proof in [`perturb`]'s
    /// own doc: `rung = MAX` (`strength = 1.0`) and `p = -P_MAX` exactly —
    /// the closest this mechanism can come to annihilating a constant.
    ///
    /// **Calls [`perturb`] itself, found recomputing its formula inline in
    /// `/review-pr`.** An earlier version wrote `1.0 + max_rung.strength() *
    /// worst_p.get()` by hand — pinning a copy of the expression rather than
    /// the function, so a change to `perturb`'s own shape (caught elsewhere
    /// by the pinned grid) would not necessarily be caught here too.
    /// `perturb(1.0, ..)` is exact against the hand-written form:
    /// multiplication by `1.0` is the identity, so the substitution changes
    /// nothing about what either assertion checks.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "CLAUDE.md: tests may assert exactly. Every step is error-free — 0.5 and \
                  -0.5 are exactly representable, strength(MAX) is exactly 1.0 (MAX/MAX), \
                  1.0 * -0.5 is exact (multiplication by 1.0 is the identity), and 1.0 - 0.5 \
                  is exact by Sterbenz's lemma — so an epsilon here (as an earlier version of \
                  this test used) tolerates up to 16 ulp of drift in strength() while the \
                  message claims exactness; same argument as \
                  identity_has_zero_strength_and_every_other_rung_is_positive, above"
    )]
    fn the_multiplier_never_reaches_zero_at_the_worst_case_corner() {
        let max_rung = Rung::new(Rung::MAX)
            .unwrap_or_else(|| unreachable!("Rung::MAX is in 0..=Rung::MAX by construction"));
        let worst_p = Direction::new(-Direction::P_MAX, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("-P_MAX is within [-P_MAX, P_MAX] by construction"));
        let multiplier = perturb(1.0, max_rung, worst_p);
        assert_eq!(multiplier, 0.5, "the worst-case corner is not exactly 0.5");
        assert!(
            multiplier > 0.0,
            "the multiplier reached zero or went negative"
        );
        // And the positive corner, for completeness: 1.0 + 1.0*0.5 = 1.5.
        let best_p = Direction::new(Direction::P_MAX, Direction::P_MAX)
            .unwrap_or_else(|| unreachable!("P_MAX is within [-P_MAX, P_MAX] by construction"));
        let upper = perturb(1.0, max_rung, best_p);
        assert_eq!(upper, 1.5, "the best-case corner is not exactly 1.5");
    }

    /// The tolerance on `P(identity)` and the corpus size it is measured
    /// over — shared by [`the_binomial_band_pins_p_identity_as_a_number`]
    /// and by [`the_defect_reading_is_distinguishable_from_the_correct_one`],
    /// so widening the band cannot silently outrun its own discriminability
    /// check. Found decoupled in `/review-pr`: the discriminability test
    /// previously asserted against a hardcoded `0.04` unrelated to this
    /// value — mutation-verified to let a widened band (`0.06`) and the
    /// named `2*MAX+1` defect both pass together, undetected by either
    /// test, catchable only by the exact pinned sequence below.
    ///
    /// **Recalibrated 2026-08-12 for `Rung::MAX = 99`, both values, not
    /// just the band — the placeholder's `N = 20_000` no longer has enough
    /// resolving power at the new `P(identity)`.** At `MAX = 8`,
    /// `P(identity) = 1/9 ≈ 11.1%` sat far enough from the `2*MAX+1` defect
    /// reading (`1/17 ≈ 5.9%`, a `~5.2` percentage-point gap) that a loose
    /// band and a modest corpus both worked. At `MAX = 99`, `P(identity) =
    /// 1%` sits much closer to its own defect reading (`1/199 ≈ 0.503%`, a
    /// `~0.497` percentage-point gap — an order of magnitude tighter) —
    /// smaller probabilities need a larger corpus to resolve to the same
    /// *relative* precision, and this document's own earlier `N = 20_000`
    /// draws left `BAND = 0.02` unable to distinguish the two at all
    /// (`0.02` alone exceeds the entire gap). `N = 200_000` brings the
    /// binomial standard error at `p ≈ 1%` to `≈0.02225` percentage
    /// points; `BAND = 0.15` percentage points is `≈6.74σ` on that scale —
    /// tight enough that `BAND + 6σ ≈ 0.283` percentage points still
    /// leaves `≈9.62σ` of margin against the `0.497`-point gap to the
    /// defect reading (see
    /// [`the_defect_reading_is_distinguishable_from_the_correct_one`] for
    /// the computed bound, not a restatement of this figure).
    const P_IDENTITY_BAND: f64 = 0.0015;

    /// See [`P_IDENTITY_BAND`].
    const P_IDENTITY_N: u32 = 200_000;

    /// Decision 10's other required test: `P(identity) = 1/(MAX + 1)`,
    /// pinned as a number over repeated draws — not an `O(m)`-order check,
    /// which both `next_range(m + 1)` and the `next_range(2*m + 1)` defect
    /// would satisfy identically.
    #[test]
    fn the_binomial_band_pins_p_identity_as_a_number() {
        let hits_count = (0..u64::from(P_IDENTITY_N))
            .filter(|&seed| Rung::draw(seed).is_identity())
            .count();
        let hits = u32::try_from(hits_count)
            .unwrap_or_else(|_| unreachable!("hits cannot exceed P_IDENTITY_N, which fits in u32"));
        let want_p = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        let observed_p = f64::from(hits) / f64::from(P_IDENTITY_N);
        // Binomial standard error at N=200,000: sqrt(p(1-p)/N) ~= 0.000222
        // at p ~= 0.01. P_IDENTITY_BAND = 0.0015 is ~6.74 sigma on that
        // scale — recalibrated 2026-08-12 for Rung::MAX=99 (was 0.02 at
        // N=20,000, ~9.0 sigma, for the placeholder MAX=8 — see
        // P_IDENTITY_BAND's own doc for why both the band and the corpus
        // size had to move together, not just the band). Generous enough
        // to survive sampling noise while tight enough that the 2*MAX+1
        // defect (p ~= 0.005025, a gap of ~0.004975 from the correct
        // reading, ~9.6-sigma beyond this band's own threshold) fails
        // loudly rather than marginally. See
        // `the_defect_reading_is_distinguishable_from_the_correct_one` for
        // the exact, computed discriminability bound, and note that an
        // off-by-one dropping the `+ 1` (`next_range(MAX)`, p ~= 0.0101,
        // a gap of only ~0.0001 from the correct reading) passes this band
        // easily — it is caught only by
        // `the_rung_sequence_over_0_to_64_is_pinned`, which is load-bearing
        // for that case, not belt-and-braces.
        assert!(
            (observed_p - want_p).abs() < P_IDENTITY_BAND,
            "observed P(identity) = {observed_p} over {P_IDENTITY_N} draws, want {want_p} +/- \
             {P_IDENTITY_BAND} ({hits} hits) — check next_range(MAX + 1) vs next_range(2*MAX + 1)"
        );
    }

    /// **A guard neither of `borbax-universe`'s two pinned goldens can
    /// provide, and the reason it exists rather than relying on them.**
    /// `every_domain_stream_is_pinned` (`borbax-rng`) pins three raw
    /// `next_u64` words per domain — it proves the *stream* moved, not
    /// that `next_range(MAX + 1)`'s rejection-sampling correctly turns
    /// those words into `0..=MAX`. The binomial test above is
    /// deliberately statistical, so it only catches a defect large enough
    /// to move `P(identity)` outside its band — a change to `m`, an
    /// off-by-one in the rejection zone, or a swapped comparison could
    /// leave `P(identity)` looking right by coincidence over 20,000 draws
    /// while every individual rung is wrong. An exact pin over a small,
    /// enumerable seed range closes that gap the way `the_universe_digest_is_pinned`
    /// closes the equivalent gap for element generation.
    #[test]
    fn the_rung_sequence_over_0_to_64_is_pinned() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |x: u64| {
            h ^= x;
            h = h.wrapping_mul(0x0100_0000_01b3);
        };
        for seed in 0..64 {
            mix(u64::from(Rung::draw(seed).0));
        }
        // Regenerated 2026-08-12: deliberate change to Rung::MAX (8 -> 99), see this
        // constant's own doc. Was 0xce60_cec5_2ce3_ee01 under MAX=8.
        assert_eq!(
            h, 0x48fe_64a7_e6b5_261b,
            "the rung sequence over seeds 0..64 moved — say whether this is a deliberate \
             change to Rung::MAX, to the draw's rejection-sampling, or a bug, and regenerate \
             deliberately if the first"
        );
    }

    /// A guard against the specific defect Decision 10 names by name:
    /// `next_range(2*MAX + 1)` gives `2*MAX + 1` outcomes and roughly
    /// halves `P(identity)`. This measures that the two readings are
    /// actually distinguishable at the corpus size
    /// [`the_binomial_band_pins_p_identity_as_a_number`] uses, not merely
    /// that they are different formulas.
    ///
    /// **Derived from [`P_IDENTITY_BAND`]/[`P_IDENTITY_N`], not a
    /// standalone literal — the bug this closes in `/review-pr`.** An
    /// earlier version asserted `(correct_p - defect_p).abs() > 0.04`, a
    /// constant with no relationship to the band it was meant to certify.
    /// Mutation-verified: widening the band to `0.06` and applying the
    /// named defect together left *both* this test and the band test
    /// green — the exact "the counter doesn't count what its message
    /// says" shape CLAUDE.md records from Task 8. The threshold here is
    /// the band plus six binomial standard errors at the shared corpus
    /// size: if the defect's miss does not clear that, the band test
    /// cannot reliably tell the two apart from sampling noise alone.
    #[test]
    fn the_defect_reading_is_distinguishable_from_the_correct_one() {
        let correct_p = 1.0 / f64::from(u32::from(Rung::MAX) + 1);
        let defect_p = 1.0 / f64::from(2 * u32::from(Rung::MAX) + 1);
        let se = (correct_p * (1.0 - correct_p) / f64::from(P_IDENTITY_N)).sqrt();
        let threshold = P_IDENTITY_BAND + 6.0 * se;
        assert!(
            (correct_p - defect_p).abs() > threshold,
            "the band ({P_IDENTITY_BAND}) is too wide to reliably see the 2*MAX+1 defect at \
             N={P_IDENTITY_N}: gap {} does not clear threshold {threshold}",
            (correct_p - defect_p).abs()
        );
    }

    /// The grid the two tests below share: every legal rung, ten base
    /// magnitudes, eighteen directions spanning the legal band.
    ///
    /// **Every rung, not just `0` and `MAX`.** Before these tests existed,
    /// `perturb` was called by exactly one test and only ever at
    /// `Rung::IDENTITY`, where `strength` is `0.0` and every
    /// algebraically-equal reformulation is bit-identical — so nothing in
    /// the workspace could see a change to the expression, and nothing at
    /// all could see a change to `Rung::strength`'s ladder between its two
    /// pinned endpoints.
    ///
    /// **Short-mantissa values cannot discriminate, so the grid avoids
    /// them.** With `base == 1.0`, `base * (1.0 + s*p)` and the distributed
    /// `base + base*s*p` are bit-identical for every `s` and `p` — a grid of
    /// round numbers pins nothing about the shape. `1.0` is kept as a
    /// readable anchor; the discrimination comes from the other nine.
    /// Likewise `k/17` is inexact for every `k` but `0` and `17`, which is
    /// why the directions are derived rather than hand-picked.
    fn pin_grid() -> impl Iterator<Item = (f64, Rung, Direction)> {
        const BASES: [f64; 10] = [
            1.0,
            -1.0,
            42.5,
            -13.7,
            0.1,
            6.022_140_76e23,
            1.380_649e-23,
            2.997_924_58e8,
            1e300,
            1e-300,
        ];
        (0..=Rung::MAX).flat_map(move |r| {
            let rung =
                Rung::new(r).unwrap_or_else(|| unreachable!("r is 0..=Rung::MAX by construction"));
            BASES.into_iter().flat_map(move |base| {
                (0..=17).map(move |k| {
                    let p = -Direction::P_MAX + f64::from(k) / 17.0;
                    let d = Direction::new(p, Direction::P_MAX).unwrap_or_else(|| {
                        unreachable!("p is within [-P_MAX, P_MAX] by construction")
                    });
                    (base, rung, d)
                })
            })
        })
    }

    /// **What [`perturb`]'s doc calls "pinned, not merely preferred", pinned
    /// — and `Rung::strength`'s ladder with it, which nothing else covers.**
    ///
    /// Measured before this test was written, on this commit: replacing
    /// `strength`'s body with a geometric ladder `(r/MAX)^2` — which
    /// preserves `strength(0) == 0.0`, `strength(MAX) == 1.0` and
    /// `strength(r) > 0.0` between, everything the rest of this module
    /// checks — passed all six gate legs, changing the perturbation
    /// magnitude in seven of the nine universes. Replacing `perturb`'s body
    /// with the distributed `base + base * strength * p` was caught only by
    /// `rung_zero_leaves_every_base_value_unchanged`, and only through a
    /// single sign-of-zero cell (`base == -0.0`, `p < 0.0`) — a value-neutral
    /// accident, not a test of the shape.
    ///
    /// The fused form `base * strength.mul_add(p, 1.0)` is *not* covered
    /// here and does not need to be: `clippy.toml`'s `disallowed_methods`
    /// and `xtask`'s `BANNED_CALLS` both reject `mul_add` outright, verified
    /// by mutation on this commit (clippy: "use of a disallowed method";
    /// xtask: a `§13.1: platform transcendental` hit). A test could not
    /// name it anyway — the textual scan reads test code too, and an
    /// `#[expect]` silences clippy while leaving `xtask` firing.
    #[test]
    fn the_perturb_expression_and_strength_ladder_are_pinned_bit_for_bit() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut cells = 0u32;
        for (base, rung, p) in pin_grid() {
            h ^= perturb(base, rung, p).to_bits();
            h = h.wrapping_mul(0x0100_0000_01b3);
            cells += 1;
        }
        // 18_000 = 100 rungs (0..=Rung::MAX=99) * 10 bases * 18 directions.
        // Was 1620 (9 rungs) under MAX=8 -- deliberate, see Rung::MAX's own doc.
        assert_eq!(cells, 18_000, "the pinned grid changed size");
        // Regenerated 2026-08-12 for Rung::MAX=99. Was 0x0839_1cd4_8bbc_19c1 under MAX=8.
        assert_eq!(
            h, 0xb5e2_f79e_a584_c7ce,
            "perturb's output moved. This is a physics change, not a tidy-up: say whether \
             the expression shape, Rung::strength's ladder, Rung::MAX or Direction::P_MAX \
             changed deliberately, and regenerate every downstream golden with it"
        );
    }

    /// **The bookkeeping check for the golden above**, in the same spirit as
    /// `the_defect_reading_is_distinguishable_from_the_correct_one`: that the
    /// grid can actually *see* the reformulation [`perturb`]'s doc names,
    /// rather than merely covering a lot of cells. A grid trimmed back to
    /// rungs `0` and `MAX`, or to round base values, would pass the golden
    /// and discriminate nothing.
    ///
    /// It is a readable second guard too: rewrite `perturb` to the
    /// distributed shape and this drops to zero, firing with a message that
    /// names the shape instead of with a moved hash.
    #[test]
    fn the_pinned_grid_still_sees_the_distributed_reshaping() {
        let differing = pin_grid()
            .filter(|&(base, rung, p)| {
                let (s, pv) = (rung.strength(), p.get());
                perturb(base, rung, p).to_bits() != (base + base * s * pv).to_bits()
            })
            .count();
        // 4063 of 18_000 cells when this was written (was 294 of 1620 under
        // Rung::MAX=8). The floor is 0, not 4063: a deliberate change to
        // `strength`'s ladder moves the count without weakening the grid,
        // and should fire the golden above alone.
        assert!(
            differing > 0,
            "the grid no longer distinguishes `base + base * strength * p` from \
             `base * (1.0 + strength * p)` — either perturb IS the distributed form now, \
             or the grid lost the full-mantissa base values that make the two differ \
             (with base == 1.0 they are bit-identical at every rung and direction)"
        );
    }
}
