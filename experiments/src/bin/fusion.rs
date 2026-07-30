//! Does the *rewritten* fusion scheme hold up — the one Task 4 will actually
//! instruct?
//!
//! `ptable.rs` measured the first draft and is kept as the record of what that
//! draft's claims survived. Three of its mechanisms were then condemned in the
//! same breath as they were reported, and this file measures their
//! replacements. Nothing here is shared with that file, deliberately: a probe
//! that reuses the code it is checking cannot find a defect in it.
//!
//! The four questions, each of which gates a sentence Task 4 would otherwise
//! assert without evidence:
//!
//! 1. **Does a drawn shell law give workable tables?** `ptable.rs` fixed the
//!    shell series at `10n^2 + 2`, which makes every universe's table the same
//!    shape. Here `k` is drawn and the shell holds `k*n^2 + 2`. Question: do
//!    the closure sets and table sizes stay inside the range §7.1 needs, or
//!    does some `k` produce a table of four elements or of nine hundred?
//! 2. **Does the contact-based valence reach 8, and is it 0 at closure without
//!    a branch?** The old `sqrt(min(outer, cap - outer))` was a fitted curve
//!    clamped at 6, and geometry measured that a *per-site* count cannot exceed
//!    `12 - 6 = 6` at all. The replacement counts unmade contacts summed over
//!    the frontier, which is a count *over* sites and has no such ceiling.
//! 3. **Is the binding peak an output?** `derive_table` took a `peak` argument
//!    and then reported the minimum it had been handed. Here the per-unit
//!    energy is contacts minus radial strain, and the peak is read off the
//!    series afterwards. The distinguishing test is structural and is in the
//!    signature: there is no `peak` parameter to pass.
//! 4. **Does the abundance-weighted valence distribution gel?** §7.2 as amended
//!    measures `f_w = <f^2>/<f>` and gels at extent `1/(f_w - 1)`. A valence
//!    ceiling of 8 is fine; a distribution centred near 8 is not. Abundance is
//!    the lever, so it has to be measured through the lever.
//!
//! Run: `cargo run --locked -p borbax-experiments --release --bin fusion`

// Same posture as `ptable.rs`: this is a spreadsheet, and every cast is a
// small count to `f64` for arithmetic. Nothing here reaches a simulation
// result.
#![expect(
    clippy::cast_precision_loss,
    clippy::as_conversions,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::too_many_lines,
    reason = "measurement probe: small counts to f64. A try_from on every line of a \
              spreadsheet would bury the arithmetic being checked"
)]

use borbax_experiments::rng::Stream;
use borbax_units::det_math;

/// Base seed for this experiment, `XOR`ed with a per-trial index.
///
/// Arbitrary but fixed.
///
/// **The `XOR` is not doing what an earlier version of this doc claimed.** It
/// said XOR was chosen over addition "so adjacent trial indices do not produce
/// adjacent seeds". For *this* constant that is vacuous: the low byte is `0x00`,
/// so `FUSION_SEED ^ s == FUSION_SEED + s` for every `s < 256`, verified
/// exhaustively — and every call site passes `s < 24`. It produces byte-identical
/// seeds to the addition it claimed to avoid.
///
/// Adjacent seeds are harmless anyway: `Stream`'s finaliser handles them, and
/// `adjacent_seeds_decorrelate_immediately` measures it. Keep the XOR because it
/// is the idiom, not because it buys anything here.
const FUSION_SEED: u64 = 0x00F0_5100;

/// Units in the `n`th packing shell: `k*n^2 + 2`.
///
/// **Only `k` is drawn; the rest is polyhedral combinatorics.** For any convex
/// simplicial polyhedron, shell `n` holds `V + E(n-1) + F(n-1)(n-2)/2`; apply
/// `3F = 2E` and the linear term vanishes, and `V - E + F = 2` leaves
/// `(F/2)n^2 + 2`. So `k = F/2 = E/3 = V - 2`, and **`shell(1) = k + 2 = V`**,
/// which forces `z = k + 2` in [`Consts`].
///
/// The `+ 2` is therefore the **Euler characteristic**, not "the antipodal pair
/// every shell carries" — an earlier version of this comment said the latter,
/// and the tetrahedron falsifies it outright: `k = 2`, shell 1 = 4 = 2*1^2 + 2,
/// and it has no antipodal vertex pair at all. The distinction never changes a
/// number here; it changes what the next person writes when they generalise to
/// a hemisphere or a torus.
///
/// `k` is how densely the base unit's own shape lets neighbours tile that
/// surface — the one thing a universe seed decides and no real table fixes.
///
/// Footnote for honesty: only `k` in {6, 7, 8, 10} admit a convex *deltahedron*,
/// i.e. a realisation from identical spheres. The rest are combinatorially fine
/// but not equal-sphere packings. In an invented chemistry the base unit need
/// not be a sphere, so this is a documentation point rather than a defect.
///
/// `k = 10` recovers `10n^2 + 2` = 12, 42, 92, the shell populations of the
/// Mackay icosahedral packing (A. L. Mackay, 1962). **The magic numbers seen
/// experimentally are the *cumulative* sizes 13, 55, 147 — not these.** An
/// earlier version attached the real-world observation to the shell
/// populations, silently swapping the two quantities the rest of this comment
/// correctly distinguishes.
///
/// **This is not corroboration and must never be read as any.** `10n^2 + 2` is
/// the vertex count of the nth icosahedral geodesic subdivision, forced by
/// Euler's formula with twelve five-fold vertices — the same formula that
/// generates this project's own sampling family `D` in {12, 42, 162}. The shell
/// law is re-deriving a consequence of `chi = 2` the codebase needed anyway,
/// not reaching toward anything real. Nothing may ever bias the draw toward
/// `k = 10` on account of the coincidence, and no claim here depends on the
/// experiment: delete the sentence naming it and everything above still holds
/// (G3; see the task preamble for the full position).
///
/// Defined for `n >= 1`. `shell_size(0)` would return 2 and is meaningless —
/// shell 0 is the single central unit, and every caller special-cases it.
const fn shell_size(k: usize, n: usize) -> usize {
    k * n * n + 2
}

/// Lateral coordination inside a shell.
///
/// Euler's formula, not a tunable. A triangulated sphere on `v` vertices has
/// exactly `3v - 6` edges, so its mean degree is `6 - 12/v` — approaching 6
/// from below for any shell large enough to matter. The twelve five-coordinate
/// sites Euler forces are the same twelve wherever they sit.
fn lateral_coordination(cap: usize) -> f64 {
    6.0 - 12.0 / cap as f64
}

/// Unmade lateral contacts on the frontier of a partly-filled shell.
///
/// The canonical filling order is a compact spherical cap growing from one
/// pole, so the frontier is a circle of latitude and its length scales as
/// `sqrt(cap*f(1-f))`. The **coefficient is fitted, not derived** — see
/// [`FRONTIER_COEFF`], and the `exact` module for the three errors the
/// derivation made. An earlier version of this comment carried that derivation
/// under the heading "Derived, and every factor has a reason", including the
/// "three of its six contacts outward" clause that measurement put at two.
///
/// **It is exactly 0 at a closure, and *two* independent factors put it there.**
/// The continuum's `f(1-f)` vanishes at both ends of a shell, and so does the
/// discrete bound's `min(outer, cap - outer)`. Probed: mutating either one
/// alone leaves `closed_shells_have_valence_zero_with_no_branch` passing,
/// because the other still zeroes it.
///
/// That makes the deletion test — `ptable.rs`'s argument for this property —
/// **uninformative here**. Continuity was the proposed replacement and it is
/// **also insufficient**: a review built a load-bearing declaration that passes
/// every test in this file and produces a bit-identical table, because both
/// neighbours of a closure are pinned at valence 1, so anything writing 0 sits
/// inside the bound. In one case it actively rewards the branch.
///
/// Continuity survives as a cheap behavioural backstop — it does catch a
/// declaration laid over a formula returning something *large*. The guarantee
/// is a source-level check that no branch reads the closure predicate, because
/// the claim is a property of the source and not of the output.
///
/// Note `compact_bound`'s own `a == 0` early return is such a branch, since
/// `a = min(outer, cap - outer)`. It is a legitimate guard against
/// `sqrt(12*0 - 3)`, and that is exactly why a source check has to name its
/// exemptions rather than scan a whole file.
fn unmade_lateral(cap: usize, outer: usize) -> f64 {
    // **There is deliberately no `if outer == 0 { return 0.0 }` here.** One
    // stood at the top of this function; a mutation probe showed the branch,
    // not the formula, was doing the work the formula was credited for — the
    // `if closed { 0 }` defect from `ptable.rs` in a different hat. It was also
    // redundant: both factors below already vanish at a closure.
    let f = outer as f64 / cap as f64;
    let continuum = continuum_at(cap, f);

    // **The continuum form overestimates tiny patches, and the bound has to be
    // the *compact* one.** A maximally compact patch of `a` sites exposes `6a`
    // less twice Harborth's internal-bond count, and any other arrangement of
    // `a` sites exposes *more* — so this is a floor on exposure, used here as a
    // cap on the continuum's estimate. An earlier version used a bare `6a`, the
    // bound for a patch whose sites touch nothing.
    //
    // **Where it actually binds, measured.** The discrete branch is taken iff
    // `compact_bound(a) < FRONTIER_COEFF*sqrt(a)`, which at 6.90 holds only at
    // `a = 1`: at `a = 2` it is 9.76 against 10, at `a = 3` 11.95 against 12,
    // and `a = 7` would need `cap > 254` where the table's caps run 8..130.
    // Confirmed over every `(cap, outer)` the table reaches — **51 of 1080**
    // cells (`k` in 6..=14, `units` 1..=120), all at
    // `min(outer, cap-outer) == 1`. The "54 of 1287" this replaces spliced two
    // different enumerations and neither produces it.
    //
    // So the justification this comment used to give — "it returned valence 3
    // where the exact count gives 2 at `outer = 3`" — is **stale by two
    // refits**. It was true at the naive 10.635; at 6.90 that case is
    // unreachable. Do not restore it as evidence, and do not delete the bound
    // on the strength of its being unreachable: it is what holds `outer == 1`.
    // By symmetry the same bound applies to the empty side.
    let smaller = if outer < cap - outer {
        outer
    } else {
        cap - outer
    };
    let discrete = compact_bound(smaller);
    if continuum < discrete {
        continuum
    } else {
        discrete
    }
}

/// The continuum frontier length as a function of the fill *fraction*, split out
/// so it can be evaluated off-lattice.
///
/// `cap as f64 * f * (1.0 - f)` associates left and must stay written that way:
/// the obvious tidy-up `cap as f64 * (f * (1.0 - f))` moves 315 of 2808 cells by
/// 1 ulp, and `peak` is a strict-`>` argmax.
///
/// `sqrt` is native on purpose: IEEE-754 specifies it exactly, so it is portable
/// without `det_math`. See the criterion in `clippy.toml`.
fn continuum_at(cap: usize, f: f64) -> f64 {
    FRONTIER_COEFF * (cap as f64 * f * (1.0 - f)).sqrt()
}

/// Unmade lateral contacts on the frontier, divided by the lateral
/// coordination — so the unit is "one absent neighbour's worth of surface" and
/// the number counts docking notches.
fn frontier_notches(cap: usize, outer: usize) -> f64 {
    unmade_lateral(cap, outer) / lateral_coordination(cap)
}

/// Continuum coefficient for the frontier length, **fitted to exact counts on
/// real triangulated spheres** rather than derived — see the `exact` module.
///
/// **The coefficient is cap-dependent and a single constant is a compromise.**
/// Least-squares against exact counts gives 6.075 at D = 12, 7.078 at D = 42
/// and 7.306 at D = 162: finite-size effects pull it down at small caps, where
/// Euler's twelve pentagons are a large fraction of the shell.
///
/// **6.90, not the asymptotic ~7.28, because of which caps this table uses.**
/// Over the drawn range `derive_table` uses caps of 8..128, and D = 12 and
/// D = 42 are real shells (`shell_size(10, 1)`, `shell_size(10, 2)`) while
/// D = 162 is not one of them. An earlier 7.30 was the D = 162 fit alone and
/// its single-level gate admitted -22%/+29% — a band over which the valence
/// series moves in up to 23 of 24 universes. 6.90 minimises the worst-case
/// error across all three levels and is best on D = 42.
///
/// `frontier_matches_exact_counts_on_a_real_sphere` gates all three levels at
/// 0.25. Bisected, the edges are [6.132767, 7.216878]; [6.2, 7.2] is that band
/// rounded inward. (An earlier "plus or minus 7%" was wrong: [6.2, 7.2] around
/// 6.90 is -10.1%/+4.3%.)
const FRONTIER_COEFF: f64 = 6.90;

/// **Minimum** unmade lateral contacts a patch of `a` sites can expose.
///
/// Harborth's result gives the **maximum** internal bonds over all arrangements
/// of `a` cells, `floor(3a - sqrt(12a - 3))`, attained by the compact ones — so
/// it is an upper bound on internal bonds and therefore a **lower** bound on
/// exposure, `6a - 2H(a)`, with equality only for maximally compact patches.
///
/// **The direction has been written three different ways in four places and
/// this is the settled one.** The title said "maximum" two lines above a body
/// saying "lower bound"; `unmade_lateral` said a patch "exposes at most `6a`
/// less twice Harborth's count", which is true only when the patch is
/// maximally compact; and the plan carried the "at least ... so it exposes at
/// most" form that this doc's own next paragraph names as the earlier error.
/// Falsified in this file's own ground truth: on `Geodesic<162>` in
/// `cap_order`, `a = 3` exposes **14** against `compact_bound(3) = 12`, because
/// a cap-order fill is not maximally compact — so "at most" is false for the
/// very fill order the derivation assumes.
///
/// No number moves today, because the bound is inert except at `a = 1` where
/// the fill is trivially compact (see the binding analysis in
/// `unmade_lateral`). It matters when the fill order is generalised, which the
/// strain accumulation contemplates.
///
/// Verified by brute force over every connected triangular-lattice animal for
/// `a = 1..10`. The result is provably positive (`>= 2*sqrt(12a-3)`), so no
/// clamp is written: a guard with no reachable case reads to a future
/// maintainer as evidence the case exists.
fn compact_bound(a: usize) -> f64 {
    if a == 0 {
        return 0.0;
    }
    let af = a as f64;
    let internal = (3.0 * af - (12.0 * af - 3.0).sqrt()).floor();
    6.0 * af - 2.0 * internal
}

/// Lateral contacts *made* inside a partly-filled shell.
///
/// Identity, not a model: summing lateral degree over the filled sites counts
/// every internal contact twice and every frontier contact once, so
/// `internal = (z_lat*a - unmade)/2`. At a full shell the frontier is empty and
/// this returns `(z_lat*cap)/2 = 3*cap - 6` exactly — Euler's edge count for a
/// triangulated sphere, which is the check that the two halves agree.
///
/// **"Identity, not a model" is true everywhere except `outer == 1`, and there
/// the clamp below is load-bearing.** The two halves use different coordination
/// conventions: `compact_bound(1) = 6` is a flat-lattice degree, while
/// `lateral_coordination(cap) = 6 - 12/cap` is the sphere's, which is strictly
/// less for every cap. So `raw` is negative at *every* `outer == 1` cell — 36
/// of 2808 across k = 6..14 and shells 1..4, worst −0.75 at cap = 8 — and the
/// clamp is what returns the correct 0, not the identity.
///
/// This is not cosmetic: `a_two_unit_cluster_has_exactly_one_contact` — the
/// test that caught the `z = 6 + k/2` defect — passes *only* because of the
/// clamp. Unclamped, `contacts_upto(c, 2)` gives 0.250 / 0.500 / 0.625 at
/// k = 6 / 10 / 14 against a tolerance of 1e-9. `the_clamp_fires_only_at_a_lone
/// _outer_site` pins the exemption rather than leaving it to be rediscovered.
///
/// The same mismatch is why `elements_one_past_a_closure_all_have_valence_one`
/// gets its result: its stated cause is "all six of its lateral contacts
/// unmade", but the shell offers 4.5–5.9, and the valence is 1 by rounding
/// `6/(6 - 12/cap)` ∈ [1.14, 1.33]. Margin to 1.5 is comfortable, so the family
/// holds — but the cause written down is not the cause.
fn lateral_made(cap: usize, outer: usize) -> f64 {
    let raw = lateral_made_raw(cap, outer);
    if raw > 0.0 { raw } else { 0.0 }
}

/// `lateral_made` before the clamp — exposed so a test can assert *where* the
/// clamp fires, instead of the clamp being an unnamed rescue inside a function
/// whose doc calls itself an identity.
fn lateral_made_raw(cap: usize, outer: usize) -> f64 {
    (lateral_coordination(cap) * outer as f64 - unmade_lateral(cap, outer)) * 0.5
}

/// Contacts made by a cluster of `units`, counted shell by shell.
///
/// **A continuum count does not work here and the first run of this probe
/// proved it.** The obvious form — `z*n/2` less a surface deficit growing as
/// `n^(2/3)` — has the deficit exceeding `n` itself for everything below
/// `n = 113`, so it collapsed to `z*n/4` across the whole table, made
/// contacts-per-unit a *constant*, and put the binding peak at `N = 1` in all
/// 24 drawn universes. The reported "peak lands on a closure in 24 of 24" was
/// an artefact of 1 being the first closure.
///
/// So contacts are accumulated from the packing itself. Each shell's sites make
/// inward contacts to the shell below and lateral contacts to their own
/// neighbours, and the inward budget is whatever outward capacity the shell
/// below has left after its own inward and lateral contacts are accounted for.
/// That recursion is what makes contacts-per-unit *rise* toward `z/2` instead
/// of sitting flat, which is what a binding peak needs.
fn contacts_upto(c: Consts, units: usize) -> f64 {
    if units <= 1 {
        return 0.0;
    }
    let mut total = 0.0_f64;
    let mut remaining = units - 1;
    // The core: one site, no lateral neighbours, no inward contacts, so its
    // whole coordination points outward.
    let (mut prev_cap, mut prev_z_lat, mut prev_inward) = (1_usize, 0.0_f64, 0.0_f64);
    let mut shell = 1_usize;
    while remaining > 0 {
        let cap = shell_size(c.k, shell);
        let take = if remaining < cap { remaining } else { cap };

        // Outward contacts the shell below still has to offer, shared out over
        // this shell's sites.
        // Divided by the shell's *capacity*, not by how many sites happen to
        // be placed: how many contacts a site makes inward is fixed by the
        // geometry, not by how full its shell currently is.
        //
        // A review found `z/cap != 1` for a lone unit on a bare core and
        // proposed dividing by `take` instead. That was measured against the
        // old `z = 6 + k/2`; with `z = k + 2` forced by the shell law,
        // `cap = shell_size(k, 1) = z`, so this is exactly 1.0 — and dividing
        // by `take` instead gives a lone unit the core's *entire* outward
        // budget, which is 8 contacts against a core that has one site.
        // Fixing `z` fixed this; the proposed repair would have broken it.
        let outward_budget = prev_cap as f64 * (c.z - prev_z_lat - prev_inward);
        let inward_per_site = outward_budget / cap as f64;

        total += inward_per_site * take as f64;
        total += lateral_made(cap, take);

        remaining -= take;
        prev_cap = cap;
        prev_z_lat = lateral_coordination(cap);
        prev_inward = inward_per_site;
        shell += 1;
    }
    total
}

/// A derived element. Every field is a function of `units` and the universe's
/// generated constants. Nothing is drawn per element.
#[derive(Debug, Clone, Copy)]
struct Derived {
    units: usize,
    shell: usize,
    outer: usize,
    mass: f64,
    radius: f64,
    valence: u8,
    affinity: f64,
    /// Binding energy per unit: made contacts, less radial strain.
    energy_per_unit: f64,
    /// Relative abundance, before normalisation.
    abundance: f64,
}

/// The generated constants of one universe's chemistry.
#[derive(Debug, Clone, Copy)]
struct Consts {
    /// Sites per unit area of shell — the drawn part of the shell law.
    k: usize,
    /// Bulk coordination of the base unit.
    z: f64,
    /// Mass of one base unit. **Dyadic, and load-bearing** — it is multiplied
    /// by an integer and never rounded, so it must land on the 1/1024 grid by
    /// itself. Probed: 0.3 steps instead of 0.25 fails exactness at N = 1.
    base_mass: f64,
    /// Mass carried away per made contact. **Need not be dyadic** — the
    /// defect is quantised explicitly at the point of use. Probed: changing
    /// the denominator from 512 to 500 leaves exactness passing, while
    /// deleting the quantisation fails it at N = 2. An earlier version of this
    /// line claimed it was load-bearing "for the same reason" as `base_mass`,
    /// which is the half of the pair that actually is.
    contact_defect: f64,
    base_radius: f64,
    /// Energy per made contact.
    eps: f64,
    /// Radial strain coefficient. The destabiliser.
    sigma: f64,
    /// How steeply abundance falls with unit count. The gel lever (§7.2).
    decay: f64,
    /// How many elements the table holds.
    ///
    /// **Drawn, and the probe must draw it too.** An earlier version fixed 120
    /// here while `generate_elements` drew 60..=120, so every figure this file
    /// reported was measured at a table size the scheme does not generate. That
    /// is how "the peak lands on a closure in 24 of 24" and "valence reaches 8
    /// in 19 of 24" both got into the plan: true at 120, false over the draw.
    n_elements: usize,
}

/// Draw one universe's constants.
///
/// `k` spans 6..=14. Below 6 a shell cannot triangulate a sphere at all; above
/// 14 the second shell exceeds 200 units and a 120-element table never leaves
/// it, so the table has one period and no families.
fn draw_consts(rng: &mut Stream) -> Consts {
    let k = 6 + rng.below(9) as usize;
    Consts {
        k,
        // **`z = k + 2` is forced, not chosen.** For any convex simplicial
        // polyhedron shell n holds `(F/2)n^2 + 2` with `F/2 = V - 2`, so
        // `shell(1) = k + 2 = V` — and shell 1 *is* the set of units touching
        // the core, which is what `z` means. An earlier `6 + k/2` agreed only
        // at k = 8, and put 16 sites in shell 1 around a core with 13
        // neighbours: three of them orbiting nothing.
        z: (k + 2) as f64,
        // Dyadic: 1, 1.25, 1.5 or 1.75. See the mass computation.
        base_mass: 1.0 + 0.25 * rng.below(4) as f64,
        // Dyadic: 1/512 .. 4/512.
        contact_defect: (1 + rng.below(4)) as f64 / 512.0,
        base_radius: 0.30 + 0.02 * rng.below(11) as f64,
        eps: 0.8 + 0.05 * rng.below(9) as f64,
        sigma: 0.02 + 0.01 * rng.below(12) as f64,
        decay: 0.04 + 0.01 * rng.below(13) as f64,
        n_elements: 60 + rng.below(61) as usize,
    }
}

/// Build a table by fusion: element `N` is `N` base units, packed.
///
/// **There is no `peak` parameter, and that is the point.** The first draft
/// took one, then reported the minimum it had been handed as an emergent
/// result — `0.013 at N = 55` was `((55-62)/62)^2`, the distance from a typed
/// constant. Here the per-unit energy is built from contacts and strain and the
/// peak is read off the finished series by `peak_of`.
fn derive_table(c: Consts, n_elements: usize) -> Vec<Derived> {
    let mut out = Vec::with_capacity(n_elements);
    for units in 1..=n_elements {
        // Which shell is filling, and how far into it.
        let (mut shell, mut filled) = (0_usize, 1_usize);
        while filled + shell_size(c.k, shell + 1) <= units {
            shell += 1;
            filled += shell_size(c.k, shell);
        }
        let outer = units - filled;
        let cap = shell_size(c.k, shell + 1);

        let contacts = contacts_upto(c, units);

        // Mass: whole base units, less the mass carried off as binding energy.
        //
        // **The two halves are guarded differently, and a mutation probe is
        // what established which guard does what.** `Mass` is fixed-point at
        // SCALE = 1024, so the series lands on its grid only if every term
        // does.
        //
        //  - `units * base_mass` is *never* rounded, so `base_mass` itself must
        //    be dyadic. Probed: drawing it in steps of 0.3 instead of 0.25
        //    fails `every_mass_is_exact_on_the_fixed_point_grid` at N = 1.
        //  - the defect is quantised explicitly, so `contact_defect` need
        //    **not** be dyadic. Probed: changing its denominator from 512 to
        //    500 leaves the test passing, while deleting the `round` below
        //    fails it at N = 2.
        //
        // An earlier version of this comment claimed both factors were
        // load-bearing. One is. That is the same shape as the `is_finite`
        // finding already on record: a comment telling a future reader that the
        // wrong half is what protects them.
        let raw_defect = contacts * c.contact_defect;
        // `round_ties_even`, matching `Mass::from_f64_quantised`. Both are
        // IEEE-754-exact, so this is not a portability fix — it is that two
        // rounding conventions for one grid let a later "consistency cleanup"
        // move a golden.
        //
        // An earlier version quoted "13 reachable ties at N = 2, 3, 4, 7".
        // Those were measured against the pre-fix frontier coefficient and did
        // not survive it. The convention still matters — one tie is enough and
        // the count was never the argument — but a quoted number no test holds
        // goes stale silently, which is this file's own recurring defect.
        let defect = (raw_defect * 1024.0).round_ties_even() / 1024.0;
        let mass = units as f64 * c.base_mass - defect;

        // Radius: the cube root of the unit count, times a bulge for the
        // partly-filled outer shell sitting proud of the packed core.
        //
        // G3: real atomic radius *decreases* across a period and jumps at a new
        // one. This does the opposite — grows through a shell, drops at
        // closure. Structurally anti-isomorphic, not a renamed trend.
        let fill = outer as f64 / cap as f64;
        // **Radius is set by which shell is occupied, and it cannot shrink.**
        // The enclosing radius of a superset of units is never smaller, and an
        // earlier `cbrt(units) * (1 + 0.55*fill)` collapsed it by 35% at every
        // closure (N=74 -> 75: 2.587 -> 1.687 on one unit added). For a shelled
        // cluster the radius is the shell index — each shell adds one lattice
        // spacing — so `shell + fill` is both monotone and the better model,
        // and it recovers the `N^(1/3)` scaling for free since N ~ k*shell^3/3.
        //
        // G3 still holds and by a cleaner route: real atomic radius *falls*
        // across a period and jumps at a new one. This rises monotonically and
        // never jumps. Structurally unlike, without being geometrically wrong.
        let radius = c.base_radius * (1.0 + shell as f64 + fill);

        // Valence: docking notches on the frontier. Zero at closure with no
        // branch — see `frontier_notches`.
        // `round_ties_even` for consistency with the mass quantiser above. An
        // earlier comment claimed exact 2.5 ties at k=8/outer=2 and outer=8;
        // measured, there are no exact half-integer notches anywhere, and even
        // pre-fix the indices were wrong. Consistency is the whole reason.
        let valence = frontier_notches(cap, outer).round_ties_even() as u8;

        // Affinity: exposed fraction of the cluster, mapped to [-1, 1].
        // Continuous across closures by construction; the first draft's version
        // jumped 1.86 at every closure and made the stable elements the
        // stickiest, which would have stopped anything accumulating.
        let shell_units = if shell == 0 {
            1
        } else {
            shell_size(c.k, shell)
        };
        let inner = filled - shell_units;
        let covered = shell_units as f64 * outer as f64 / cap as f64;
        let surface = units as f64 - (inner as f64 + covered);
        let affinity = 2.0 * (surface / units as f64) - 1.0;

        // Per-unit binding energy: made contacts, less radial strain.
        //
        // **The strain term is the invented destabiliser, and the obvious one
        // is a G3 breach.** The familiar way to bend a binding curve back down
        // is the semi-empirical mass formula's Coulomb term, and reaching for
        // it imports real nuclear physics with no data file in sight. This is a
        // different mechanism: a lattice cannot tile a sphere, so each shell is
        // stretched over a larger radius than the one below it and carries a
        // strain that grows as `n^2`. That competes with the contact term's
        // saturation and produces a maximum. (An earlier version quoted a closed
        // form `eps*z*coeff/(8*sigma)` with `coeff` undefined; solving for it
        // across the grid gives 0.80 to 5.00, so no such constant exists and the
        // claim was unfalsifiable as written.)
        //
        // **The exponent 2/3 is CHOSEN, not derived, and it is load-bearing.**
        // Executing the stated derivation does not produce it: per-*shell*
        // strain proportional to `n^2`, summed and divided by `N`, gives
        // `Sum n^2 / (k Sum n^2) = 1/k` — a *constant*, confirmed numerically
        // (0.154..0.161 against 1/6 at k = 6; 0.067..0.070 against 1/14 at
        // k = 14, approaching 1/k from below). Reaching `N^(2/3)` needs
        // per-*site* strain proportional to `n^2`, which is a different claim
        // and is not the one written above.
        //
        // **And the choice decides the headline result.** Accumulating strain
        // shell-by-shell — as `contacts_upto` does, for the identical reason —
        // moves the peak in **48.1% of drawn universes** (3170 of 6588 over
        // k x eps x sigma x n_elements): closure/edge/mid-shell goes from
        // 86.6/9.9/3.6% to 66.8/33.2/0.0%. So "the peak's position is a
        // function of the generated constants and of nothing typed in" — which
        // an earlier version of this comment asserted — is **false**, and the
        // `19 of 24` census below is a sample from the closed-form row only.
        //
        // This is the same defect `contacts_upto` was rewritten to fix, in its
        // direct competitor in the same subtraction: see that function's doc,
        // which records that a continuum count "put the binding peak at N = 1
        // in all 24 drawn universes". Deferred deliberately, not overlooked —
        // the remedy is shell-by-shell accumulation and it requotes every
        // number in this file, so it belongs with Task 20's battery. Precedence
        // 2 (correctness of physics), NOT a G3 breach: no data, no calibration,
        // no transfer function.
        let strain = c.sigma * det_math::cbrt(units as f64 * units as f64);
        let energy_per_unit = c.eps * contacts / units as f64 - strain;

        // Abundance: fusion builds heavy clusters from light ones, so each
        // extra unit costs a step and abundance falls geometrically. Constant
        // per-step survival compounded over N sequential additions gives
        // `r^N = exp(N ln r)`, and that is the whole law.
        //
        // Two honest caveats. The mechanism gives `r^(N-1)` and the code writes
        // `r^N`, which normalises away. And "constant per-step survival" is an
        // *assumption* stated as a consequence, in a file that computes a
        // per-element decay rate a few lines above.
        // **The energy tilt is deleted, and that is a finding rather than a
        // simplification.** It read `+ 0.6 * energy_per_unit`, described as
        // "tilted by how well the cluster binds". Three things were wrong with
        // it. It pushed the *wrong way* — removing it moves nominal p_c from
        // 0.358 to 0.381, better than the best the entire `decay` sweep
        // reaches. The stated compounding mechanism derives `E_total`, not
        // `E_per_unit`; implementing the stated story makes the *heaviest*
        // element the most abundant and collapses p_c to 0.112. And `0.6` was
        // typed in while every neighbouring constant is drawn, which is
        // structurally the same defect as `peak`-as-an-argument that this
        // scheme exists to remove.
        //
        // What remains is derived: sequential addition with a constant
        // per-step survival probability gives `r^N = exp(N ln r)`.
        let abundance = det_math::exp(-(units as f64) * c.decay);

        out.push(Derived {
            units,
            shell,
            outer,
            mass,
            radius,
            valence,
            affinity,
            energy_per_unit,
            abundance,
        });
    }
    out
}

/// Where the per-unit binding energy peaks — read off the series, not passed in.
fn peak_of(t: &[Derived]) -> usize {
    t.iter()
        .fold((0_usize, f64::MIN), |(bn, be), e| {
            if e.energy_per_unit > be {
                (e.units, e.energy_per_unit)
            } else {
                (bn, be)
            }
        })
        .0
}

/// Cumulative shell closures inside a table of `n` elements.
fn closures(k: usize, n: usize) -> Vec<usize> {
    let mut out = vec![1_usize];
    let (mut total, mut shell) = (1_usize, 0_usize);
    loop {
        shell += 1;
        total += shell_size(k, shell);
        if total > n {
            return out;
        }
        out.push(total);
    }
}

/// Weighted functionality `f_w = <f^2>/<f>` and the gel extent.
///
/// `p_c = 1/(f_w - 1)` is **Cohen/Erez/ben-Avraham/Havlin (2000)** and the
/// Flory-Stockmayer gel point, not Molloy-Reed — Molloy-Reed gives the
/// *existence* criterion `f_w > 2`, a yes/no test rather than a threshold in an
/// extent coordinate. An earlier version of this comment misattributed it, and
/// the correction reached the spec before it reached here.
///
/// **A high `p_c` is not by itself "the safe direction".** `p_c` moves under
/// any rescaling of the mean degree — bonds-per-node at gel is
/// `p_c*<f>/2 = <f>/(2(f_w-1))` — and roughly 84% of what the `decay` sweep
/// reports as a gain is that rescaling. Neither coordinate says where the
/// beaker actually sits; that needs an operating point, which is §7.2's
/// steady-state condition and Task 20's battery, not this file. Terminators (`f <= 1`) are counted in
/// the averages because they are what stops chains, which is the whole reason
/// the ceiling exists.
/// The three gel statistics, named.
///
/// **Not a `(f64, f64, f64)`.** Three same-typed floats read positionally at
/// four call sites is a transposed destructure away from silently retitling a
/// printed column, and the same expression is destined for Task 20's battery
/// where those numbers gate a universe. Fields make the transposition a compile
/// error instead of a wrong table.
#[derive(Debug, Clone, Copy)]
struct Gel {
    /// Weight-average functionality, `<f^2>/<f>`.
    f_w: f64,
    /// Critical conversion extent, `1/(f_w - 1)`.
    p_c: f64,
    /// Bonds per node at gel, `p_c*<f>/2` — the coordinate invariant under
    /// degree-independent thinning.
    bonds_per_node: f64,
}

fn gel_extent(weights: &[(u8, f64)]) -> Gel {
    let total: f64 = weights.iter().map(|&(_, w)| w).sum();
    let mean: f64 = weights.iter().map(|&(f, w)| f64::from(f) * w / total).sum();
    let mean_sq: f64 = weights
        .iter()
        .map(|&(f, w)| f64::from(f) * f64::from(f) * w / total)
        .sum();
    // `f_w = <f^2>/<f>` is *undefined* at zero mean degree, not zero. An
    // earlier version returned 0.0, which is a fabricated value for a quantity
    // that does not exist and would flow into a printed column looking like a
    // measurement. Unreachable over the drawn ranges (abundance bottoms out at
    // 4.6e-9 per element), so making it loud costs nothing.
    assert!(
        mean > 0.0,
        "gel_extent: mean degree {mean}, so f_w is undefined"
    );
    let f_w = mean_sq / mean;
    if f_w <= 1.0 {
        return Gel {
            f_w,
            p_c: f64::INFINITY,
            bonds_per_node: f64::INFINITY,
        };
    }
    let p_c = 1.0 / (f_w - 1.0);
    // Bonds per node at gel: `p_c * <f> / 2`. **The invariant coordinate** —
    // under degree-independent thinning at rate q, `f_w' = 1 + q(f_w - 1)` so
    // `p_c' = p_c/q`, and the q cancels here.
    //
    // It is **not** an operating point, and reporting it does not by itself fix
    // what §7.2 was missing: the ratio of where-the-beaker-sits to
    // where-it-gels is the same number in either coordinate. Where the beaker
    // sits needs the steady-state condition from §9.1/§9.4's kinetics, which
    // Task 4 has no rate constants for. Task 20's battery does.
    Gel {
        f_w,
        p_c,
        bonds_per_node: p_c * mean * 0.5,
    }
}

fn main() {
    println!("== the rewritten fusion scheme: does it hold up? ==\n");

    // -- Q1: does a drawn shell law give workable tables? -----------------
    println!("-- Q1: shell law k*n^2 + 2, k drawn --\n");
    println!(
        "{:>3}  {:>28}  {:>10}",
        "k", "closures within 120", "2nd shell"
    );
    for k in 6..=14 {
        println!(
            "{k:>3}  {:>28}  {:>10}",
            format!("{:?}", closures(k, 120)),
            shell_size(k, 2)
        );
    }

    // -- Q2 and Q3, over drawn universes ----------------------------------
    println!("\n-- Q2/Q3: valence reach and the peak, over 24 drawn universes --\n");
    println!(
        "{:>4} {:>3} {:>4} {:>7} {:>5} {:>5} {:>7} {:>8} {:>8}",
        "seed", "k", "n", "max val", "clos", "peak", "p@clos", "mass OK", "gel p_c"
    );
    let mut peaks_on_closure = 0_usize;
    let mut all_masses_exact = true;
    let mut max_valences = Vec::new();
    let mut gel_extents = Vec::new();
    for seed in 0..24_u64 {
        let mut rng = Stream::new(FUSION_SEED ^ seed);
        let c = draw_consts(&mut rng);
        let table = derive_table(c, c.n_elements);

        let max_val = table.iter().map(|e| e.valence).max().unwrap_or(0);
        max_valences.push(max_val);

        let peak = peak_of(&table);
        let clos = closures(c.k, c.n_elements);
        let on_closure = clos.contains(&peak);
        if on_closure {
            peaks_on_closure += 1;
        }

        // Mass must land on the 1/1024 grid for every element.
        let exact = table
            .iter()
            .all(|e| (e.mass * 1024.0 - (e.mass * 1024.0).round()).abs() < 1e-9);
        if !exact {
            all_masses_exact = false;
        }

        let weights: Vec<(u8, f64)> = table.iter().map(|e| (e.valence, e.abundance)).collect();
        let p_c = gel_extent(&weights).p_c;
        // Carries the seed so the sort below has a real tie-break (§13.1).
        gel_extents.push((p_c, seed));

        println!(
            "{seed:>4} {:>3} {:>4} {max_val:>7} {:>5} {peak:>5} {:>7} {:>8} {p_c:>8.3}",
            c.k,
            c.n_elements,
            clos.len(),
            if on_closure { "yes" } else { "no" },
            if exact { "yes" } else { "NO" }
        );
    }
    println!(
        "\npeak lands on a closure in {peaks_on_closure} of 24.\n  \
         The sawtooth IS structural -- a closed shell has no frontier, so it makes every \
         lateral\n  contact available to it and maximises contacts per unit, giving a local \
         maximum at\n  every closure. But the GLOBAL max of a series truncated mid-shell can \
         sit at the table\n  edge instead, and most misses here are exactly that. An earlier \
         version measured this\n  at a fixed 120 elements, read 24 of 24, and called it \
         structural without qualification.\n  What the constants choose is *which* closure \
         wins: {} distinct peaks over 24 draws.",
        {
            let mut ps: Vec<usize> = Vec::new();
            for seed in 0..24_u64 {
                let mut r = Stream::new(FUSION_SEED ^ seed);
                let cc = draw_consts(&mut r);
                let p = peak_of(&derive_table(cc, cc.n_elements));
                if !ps.contains(&p) {
                    ps.push(p);
                }
            }
            ps.len()
        }
    );
    println!(
        "every mass exact on the 1/1024 grid in all 24 universes: {}",
        if all_masses_exact { "yes" } else { "NO" }
    );
    let reach8 = max_valences.iter().filter(|&&v| v >= 8).count();
    println!("max valence >= 8 in {reach8} of 24 universes");

    // -- Q4: the gel lever, in detail on one universe ---------------------
    println!("\n-- Q4: does abundance keep the beaker off the gel point? --\n");
    let mut rng = Stream::new(FUSION_SEED);
    let c = draw_consts(&mut rng);
    let table = derive_table(c, c.n_elements);
    println!(
        "universe: k={} z={:.1} eps={:.2} sigma={:.3}",
        c.k, c.z, c.eps, c.sigma
    );
    println!(
        "closures: {:?}   peak: N={}\n",
        closures(c.k, c.n_elements),
        peak_of(&table)
    );

    println!(
        "{:>4} {:>5} {:>5} {:>8} {:>7} {:>4} {:>8} {:>8} {:>9}",
        "N", "shell", "outer", "mass", "radius", "val", "affin", "E/unit", "abundance"
    );
    let total_ab: f64 = table.iter().map(|e| e.abundance).sum();
    for e in table
        .iter()
        .filter(|e| e.outer == 0 || e.outer == 1 || e.units % 19 == 0)
    {
        println!(
            "{:>4} {:>5} {:>5} {:>8.3} {:>7.3} {:>4} {:>8.3} {:>8.3} {:>9.5}",
            e.units,
            e.shell,
            e.outer,
            e.mass,
            e.radius,
            e.valence,
            e.affinity,
            e.energy_per_unit,
            e.abundance / total_ab
        );
    }

    // Valence histogram, unweighted and abundance-weighted. The gap between
    // the two is the entire argument for abundance being the lever.
    let top = table.iter().map(|e| e.valence).max().unwrap_or(0);
    println!("\n{:>10} {:>8} {:>10}", "valence", "count", "ab.weight");
    for v in 0..=top {
        let count = table.iter().filter(|e| e.valence == v).count();
        let w: f64 = table
            .iter()
            .filter(|e| e.valence == v)
            .map(|e| e.abundance)
            .sum::<f64>()
            / total_ab;
        println!("{v:>10} {count:>8} {w:>10.4}");
    }

    let flat: Vec<(u8, f64)> = table.iter().map(|e| (e.valence, 1.0)).collect();
    let weighted: Vec<(u8, f64)> = table.iter().map(|e| (e.valence, e.abundance)).collect();
    let flat_gel = gel_extent(&flat);
    let (fw_flat, pc_flat, bc_flat) = (flat_gel.f_w, flat_gel.p_c, flat_gel.bonds_per_node);
    let w_gel = gel_extent(&weighted);
    let (fw_w, pc_w, bc_w) = (w_gel.f_w, w_gel.p_c, w_gel.bonds_per_node);
    println!(
        "\nunweighted (every element equally common):  f_w = {fw_flat:.3}  extent {pc_flat:.3}  bonds/node {bc_flat:.3}"
    );
    println!(
        "abundance-weighted:                          f_w = {fw_w:.3}  extent {pc_w:.3}  bonds/node {bc_w:.3}"
    );
    println!(
        "  -> extent moved {:+.1}%, bonds/node moved {:+.1}%. The gap between those two is\n     \
         mean-degree rescaling, not a change in how close the beaker is to gelling.",
        100.0 * (pc_w / pc_flat - 1.0),
        100.0 * (bc_w / bc_flat - 1.0)
    );
    // -- The gel lever, swept ---------------------------------------------
    //
    // §7.2 as amended by PR #8: a distribution *centred* near 8 gels at extent
    // 0.20, one dominated by f = 2 with a thin tail gels at 0.73. The valence
    // series here is fixed by the packing, so the only lever is how much
    // abundance sits on the low-valence end. This is the number Task 4's
    // battery has to hold, so the safe band has to be measured, not asserted.
    println!(
        "\n-- the gel lever: abundance decay swept, k = {} --\n",
        c.k
    );
    println!(
        "{:>7} {:>9} {:>9} {:>11} {:>19}",
        "decay", "f_w", "gel p_c", "bonds/node", "mass on val <= 2"
    );
    for step in 0..13_u64 {
        let decay = 0.04 + 0.01 * step as f64;
        let t = derive_table(Consts { decay, ..c }, c.n_elements);
        let tot: f64 = t.iter().map(|e| e.abundance).sum();
        let low: f64 = t
            .iter()
            .filter(|e| e.valence <= 2)
            .map(|e| e.abundance)
            .sum::<f64>()
            / tot;
        let w: Vec<(u8, f64)> = t.iter().map(|e| (e.valence, e.abundance)).collect();
        let g = gel_extent(&w);
        let (fw, pc, bc) = (g.f_w, g.p_c, g.bonds_per_node);
        println!("{decay:>7.2} {fw:>9.3} {pc:>9.3} {bc:>11.3} {low:>19.3}");
    }

    // Sorted once, then read at three positions, because `f64::min`/`max` are
    // disallowed (their tie-breaking differs by target).
    //
    // The seed is the ID tie-break §13.1 requires, and it is in the *key*
    // rather than in a comment. An earlier version relied on `sort_by` being
    // stable over seed-ordered input, which is true and invisible: swapping in
    // `sort_unstable_by` — the obvious cleanup for a `Vec<f64>` — would have
    // broken it silently, and nothing at this call site would have said so.
    // **The finiteness check is what makes the comparator below a total order,
    // and it is not decoration.** `partial_cmp(..).unwrap_or(Equal)` is not a
    // total order: with a NaN key it is intransitive — for (1.0,0), (NaN,1),
    // (0.5,2) it reports a<b, b<c, a>c — and `Equal.then(seed)` does not
    // rescue it, it makes std's "does not correctly implement a total order"
    // panic *more* likely (70.3% against 54.5% of NaN inputs, measured over
    // 360k fuzz cases). Trading a silent misorder for a loud panic is not the
    // same as having an order.
    //
    // NaN is unreachable today: `gel_extent` yields non-finite only from its
    // two early returns, and abundance cannot underflow to zero over the drawn
    // ranges (worst-case table total 5.763, smallest per-element term
    // 4.587e-9 — that needs `units*decay ~ 745` against a bound of 120*0.16).
    // So assert it, and let `total_cmp` be exactly what the library's own
    // `canonical_cmp` is: a total order made safe by a checked precondition.
    assert!(
        gel_extents.iter().all(|&(p, _)| p.is_finite()),
        "gel extents must be finite before sorting: {gel_extents:?}"
    );
    let mut sorted = gel_extents.clone();
    #[expect(
        clippy::disallowed_methods,
        reason = "§13.4: total_cmp is banned because a runtime NaN's sign bit \
                  is architecture-dependent. The assert above rules NaN out, \
                  which is the same precondition Unit::canonical_cmp relies on; \
                  the seed keeps the order total across equal p_c"
    )]
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    println!(
        "\ngel extents across the 24 universes: min {:.3}  median {:.3}  max {:.3}",
        sorted.first().map_or(0.0, |&(p, _)| p),
        sorted.get(sorted.len() / 2).map_or(0.0, |&(p, _)| p),
        sorted.last().map_or(0.0, |&(p, _)| p)
    );
}

#[cfg(test)]
// Scoped to the tests rather than the module: the binary itself indexes
// nothing, so a module-level expectation is *unfulfilled* there and
// `-D warnings` rejects it. Same trap as `[lints] workspace = true` — the
// attribute is present, believed, and in the wrong place.
#[expect(
    clippy::indexing_slicing,
    reason = "windows(2)/windows(3) indexing, where the slice length is the iterator's \
              own guarantee"
)]
mod tests {
    use super::*;

    /// Exact frontier counts on a real triangulated sphere.
    ///
    /// **This module is the specification, and the continuum formula is the
    /// approximation being checked against it.** The first version of this probe
    /// had the continuum formula and no ground truth, and it was wrong by 1.46x
    /// for two compounding reasons that no test could see:
    ///
    /// - *"each boundary site points three of its six contacts outward"* is wrong.
    ///   Measured here, unmade-contacts-per-boundary-site converges to **2.000**.
    /// - *site spacing `sqrt(4*pi/cap)`* treats the area per site as `s^2`. On a
    ///   triangular lattice it is `(sqrt(3)/2)*s^2`, so the spacing was 7.5% low.
    /// - and a third, which the first correction of this comment missed: a lattice
    ///   boundary **zig-zags**, passing through ~10% more sites than `L/s`.
    ///   Measured boundary-site density near the equator at D = 162 is 3.61-3.64
    ///   per `sqrt(cap*f(1-f))`, against the 3.299 the `L/s` route predicts.
    ///
    /// The first two alone compound to 1.61x and give 6.598, **not** the fitted
    /// value — so a reader re-deriving from them would "fix" the constant
    /// downward. The derivation that reproduces the measurement is Cauchy's
    /// line-crossing formula with bond density `2*sqrt(3)/s^2` over
    /// `L = 4*pi*sqrt(f(1-f))`: `16*sqrt(3)*sqrt(sqrt(3)/(8*pi))` = 7.2751. That
    /// is asymptotic; see [`FRONTIER_COEFF`] for why the shipped constant is lower.
    ///
    /// Both errors inflate the count, valence ran ~46% high everywhere, and the
    /// two headline chemistry conclusions were consequences of the error rather
    /// than of the design.
    ///
    /// Test-only: the ground truth exists to check the formula, not to be used by
    /// it. `#[cfg(test)]` rather than an `expect(dead_code)`, so that if the probe
    /// ever starts calling it in anger the gate says so.
    #[expect(
        clippy::indexing_slicing,
        reason = "fixed-size geodesic arrays indexed by their own 0..D loop bounds"
    )]
    mod exact {
        use borbax_experiments::geodesic::Geodesic;

        /// Adjacency of a geodesic sphere, by nearest-neighbour threshold.
        ///
        /// A geodesic's edges are its shortest vertex-vertex distances; the next
        /// shell of distances sits far enough above that a 1.3x threshold on the
        /// minimum separates them cleanly. Verified by the degree histogram: any
        /// triangulated sphere must have exactly twelve degree-5 vertices and the
        /// rest degree 6, and `degrees` asserts it.
        pub(super) fn adjacency<const D: usize>(g: &Geodesic<D>) -> Vec<Vec<usize>> {
            let mut min_d2 = f64::MAX;
            for i in 0..D {
                for j in (i + 1)..D {
                    let d2 = dist2(g.dirs[i], g.dirs[j]);
                    if d2 < min_d2 {
                        min_d2 = d2;
                    }
                }
            }
            let cut = min_d2 * 1.3 * 1.3;
            (0..D)
                .map(|i| {
                    (0..D)
                        .filter(|&j| j != i && dist2(g.dirs[i], g.dirs[j]) < cut)
                        .collect()
                })
                .collect()
        }

        fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
            let (dx, dy, dz) = (a[0] - b[0], a[1] - b[1], a[2] - b[2]);
            dx * dx + dy * dy + dz * dz
        }

        /// Degree histogram, as `(n_degree_5, n_degree_6, n_other)`.
        pub(super) fn degrees(adj: &[Vec<usize>]) -> (usize, usize, usize) {
            let mut out = (0, 0, 0);
            for a in adj {
                match a.len() {
                    5 => out.0 += 1,
                    6 => out.1 += 1,
                    _ => out.2 += 1,
                }
            }
            out
        }

        /// Fill order: a compact spherical cap growing from one pole. This is the
        /// *canonical filling order* the continuum derivation assumes, so the
        /// ground truth has to use it too or the comparison is meaningless.
        pub(super) fn cap_order<const D: usize>(g: &Geodesic<D>) -> Vec<usize> {
            let mut idx: Vec<usize> = (0..D).collect();
            // Sort by height along +z, with the index as tie-break (§13.4).
            //
            // Same reasoning as the gel-extent sort: `partial_cmp(..)
            // .unwrap_or(Equal)` is not a total order under NaN, so the
            // precondition has to be checked rather than assumed. Here the
            // components are coordinates of unit vectors on a geodesic, finite
            // by construction — asserted anyway, because "by construction" is
            // what the last three defects in this file were also confident of.
            assert!(
                g.dirs.iter().all(|d| d[2].is_finite()),
                "geodesic directions must be finite"
            );
            #[expect(
                clippy::disallowed_methods,
                reason = "§13.4: NaN is excluded by the assert above, and the \
                          index tie-break keeps the order total across equal \
                          heights — which are common on a geodesic"
            )]
            idx.sort_by(|&a, &b| g.dirs[b][2].total_cmp(&g.dirs[a][2]).then(a.cmp(&b)));
            idx
        }

        /// Unmade lateral contacts after filling the first `a` sites of `order`:
        /// the number of filled-empty adjacent pairs. Exact, by construction.
        pub(super) fn unmade(adj: &[Vec<usize>], order: &[usize], a: usize) -> usize {
            let mut filled = vec![false; adj.len()];
            for &i in order.iter().take(a) {
                filled[i] = true;
            }
            let mut n = 0;
            for (i, nbrs) in adj.iter().enumerate() {
                if filled[i] {
                    n += nbrs.iter().filter(|&&j| !filled[j]).count();
                }
            }
            n
        }
    }

    fn consts() -> Consts {
        let mut rng = Stream::new(FUSION_SEED);
        draw_consts(&mut rng)
    }

    /// The seed-0 fixture with `k` overridden — and `z` moved with it.
    ///
    /// **`consts_with_k(k)` is not a valid universe.** `draw_consts`
    /// derives `z = k + 2` from the shell law, so overriding `k` alone leaves
    /// `z` at the seed-0 draw and builds a state the generator cannot produce.
    /// The assertions that used it key on `k` only, so nothing passed falsely
    /// — but `contacts_upto` reads `c.z`, so `energy_per_unit` and `mass` in
    /// those tables were computed off an unreachable universe, and any later
    /// energy or mass assertion added to them would have measured it. Same
    /// defect as the hardcoded `n_elements` recorded above: a fixture drifting
    /// from the generator, invisible until something reads the drifted field.
    fn consts_with_k(k: usize) -> Consts {
        Consts {
            k,
            z: (k + 2) as f64,
            ..consts()
        }
    }

    /// The emergence claim, as a property rather than a deletion test.
    ///
    /// `ptable.rs` had an `if closed { 0 }` branch and earned its result by
    /// *deleting* it. Here there is no branch to delete: `f(1-f)` is zero at
    /// both ends of every shell, so the assertion is on the mechanism.
    #[test]
    fn closed_shells_have_valence_zero_with_no_branch() {
        let c = consts();
        for e in derive_table(c, 200) {
            if e.outer == 0 {
                assert_eq!(e.valence, 0, "N={} sits at a closure", e.units);
            }
        }
    }

    /// **The specification for `unmade_lateral`, and the test whose absence
    /// let a 1.46x error stand.**
    ///
    /// Builds a real geodesic sphere, derives its adjacency, fills it in
    /// spherical-cap order — the same canonical order the continuum derivation
    /// assumes — and counts filled-to-empty edges exactly. The formula must
    /// track that. Nothing else in this file could have caught the error,
    /// because every other test compares the formula against itself.
    #[test]
    fn frontier_matches_exact_counts_on_a_real_sphere() {
        use borbax_experiments::geodesic::Geodesic;

        // **All three levels, because the coefficient is cap-dependent.** An
        // earlier version gated D = 162 alone — the one level `derive_table`
        // never uses — and admitted -22%/+29% around the constant, a band over
        // which the valence series moves in up to 23 of 24 universes. D = 12
        // and D = 42 are real shells: `shell_size(10, 1)` and
        // `shell_size(10, 2)`. Three legs at 0.25 admit [6.2, 7.2].
        fn check<const D: usize>() -> (f64, usize) {
            let built = Geodesic::<D>::build();
            assert!(built.is_ok(), "geodesic D={D} failed to build");
            // Sentinel must FAIL the caller's `worst < 0.25`, not satisfy it.
            // `(0.0, 0)` was a *passing* value in the else-arm of the test whose
            // absence let a 1.46x error stand — unreachable only because of the
            // assert above it, so relaxing that assert would have converted this
            // into a silent pass.
            let Ok(g) = built else {
                return (f64::INFINITY, 0);
            };
            let adj = exact::adjacency(&g);

            // Euler: any triangulated sphere has exactly twelve degree-5
            // vertices. If this fails the adjacency threshold is wrong and
            // every count below is meaningless.
            let (d5, _, other) = exact::degrees(&adj);
            assert_eq!((d5, other), (12, 0), "D={D}: degree histogram is wrong");

            let order = exact::cap_order(&g);
            let (mut worst, mut worst_at) = (0.0_f64, 0);
            for a in 1..D {
                let truth = exact::unmade(&adj, &order, a) as f64;
                if truth > 0.0 {
                    let rel = (unmade_lateral(D, a) - truth).abs() / truth;
                    if rel > worst {
                        worst = rel;
                        worst_at = a;
                    }
                }
            }
            (worst, worst_at)
        }

        for (d, (worst, at)) in [
            (12, check::<12>()),
            (42, check::<42>()),
            (162, check::<162>()),
        ] {
            assert!(
                worst < 0.25,
                "D={d}: worst relative error {worst:.4} at outer={at} \
                 — the formula has drifted from the geometry it approximates"
            );
        }
    }

    /// The enclosing radius of a superset of units cannot be smaller. An
    /// earlier form collapsed it 35% at every closure, which — since G1
    /// measured binding at ~93% size and ~7% shape — would have read as
    /// closed-shell elements binding a different partner set entirely.
    #[test]
    fn radius_never_shrinks() {
        for k in 6..=14_usize {
            let c = consts_with_k(k);
            for pair in derive_table(c, 200).windows(2) {
                assert!(
                    pair[1].radius >= pair[0].radius,
                    "k={k}: radius fell from {} to {} at N={}",
                    pair[0].radius,
                    pair[1].radius,
                    pair[1].units
                );
            }
        }
    }

    /// Two units share exactly one contact. Not a modelling choice — and the
    /// earlier tolerance of 0.35 was exactly wide enough to span the defect it
    /// should have caught (`z = 6 + k/2` gave 0.81..1.13 across k).
    #[test]
    fn a_two_unit_cluster_has_exactly_one_contact() {
        for k in 6..=14_usize {
            let c = consts_with_k(k);
            let n = contacts_upto(c, 2);
            assert!(
                (n - 1.0).abs() < 1e-9,
                "k={k}: a 2-unit cluster made {n} contacts"
            );
        }
    }

    /// `z` is forced by the shell law: `shell_size(k, 1) == z`. Two hand-written
    /// copies of a relation are how they drift apart.
    #[test]
    fn coordination_agrees_with_the_first_shell() {
        for k in 6..=14_usize {
            let mut rng = Stream::new(FUSION_SEED ^ k as u64);
            let c = draw_consts(&mut rng);
            assert!(
                (shell_size(c.k, 1) as f64 - c.z).abs() < f64::EPSILON,
                "k={}: shell 1 holds {} sites but z is {}",
                c.k,
                shell_size(c.k, 1),
                c.z
            );
        }
    }

    /// A cheap behavioural backstop, **not** the guarantee.
    ///
    /// It catches a declaration laid over a formula returning something *large*
    /// at a closure — `notches + 3` fails it. It does not catch one writing 0
    /// over a formula returning ~1, because both neighbours of a closure are
    /// pinned at valence 1 and the declaration then sits inside the bound.
    /// Measured: a mutation putting every closure at valence 1 passes the whole
    /// suite. The guarantee is the source-level check.
    #[test]
    fn valence_is_continuous_across_closures() {
        let c = consts();
        for pair in derive_table(c, 200).windows(2) {
            let jump = i16::from(pair[1].valence) - i16::from(pair[0].valence);
            assert!(
                jump.abs() <= 2,
                "N={} -> {}: valence jumped {} ({} -> {}), closure={} -> {}",
                pair[0].units,
                pair[1].units,
                jump,
                pair[0].valence,
                pair[1].valence,
                pair[0].outer == 0,
                pair[1].outer == 0
            );
        }
    }

    /// §7.1 needs *families* — "everything in this column behaves alike" is the
    /// moment the table becomes learnable rather than watchable. Position in a
    /// drawn period was the old answer and it caused nothing.
    ///
    /// The fusion answer: the element one unit past any closure is a closed
    /// core plus a single adhering unit, and that unit has all six of its
    /// lateral contacts unmade wherever it sits. So it has exactly one docking
    /// notch, in every shell and every universe. That is a family with a cause.
    #[test]
    fn elements_one_past_a_closure_all_have_valence_one() {
        for k in 6..=14_usize {
            let c = consts_with_k(k);
            let t = derive_table(c, 200);
            for close in closures(k, 199) {
                let next = t.get(close).copied();
                assert!(
                    next.is_some_and(|e| e.valence == 1),
                    "k={k}: N={} follows closure {close} with valence {:?}, not 1",
                    close + 1,
                    next.map(|e| e.valence)
                );
            }
        }
    }

    /// **Valence 8 is not reachable in a table this scheme generates, and
    /// saying so is the finding.**
    ///
    /// A per-*site* count caps at `z - 6`, which is why the fitted curve this
    /// replaced could not exceed 6 even before its clamp. Summing over the
    /// frontier removes that ceiling in principle — measured, valence 9 appears
    /// once the fourth shell opens, past N = 200.
    ///
    /// But `n_elements` is drawn on 60..=120, and within that range the ceiling
    /// is **4 to 6**.
    ///
    /// This assertion has been wrong three times, which is the point of writing
    /// the history down: `5..=7` from a guess, which failed with a 4 in it;
    /// widened to `4..=7`, which was loose enough to survive the frontier
    /// coefficient being refitted from 7.30 to 6.90; now `4..=6`, measured. A
    /// range assertion that survives a change to the quantity it measures is
    /// not measuring it.
    ///
    /// An earlier version of the test asserted `max >= 8` and passed by
    /// measuring at 300 elements, outside anything the scheme builds — the same
    /// "measured at a size the plan does not generate" defect that put two
    /// false claims into the plan's prose.
    ///
    /// §7.2 wants a thin high-valence tail. What it gets is a ceiling of 6 at
    /// best and 4 at worst — and §7.2 asked for a *shape*, not a magnitude:
    /// abundance at valence >= 3 spans 0.077..0.566 across universes, so the
    /// shape it asked for is present in every one.
    /// Whether that is enough functionality is Task 5's battery to answer; it
    /// is not something this file may assert either way.
    /// **The gel band, asserted — because prose quoting it has drifted twice.**
    ///
    /// `FRONTIER_COEFF` was refitted 7.30 → 6.90 in this branch and the gel
    /// extents moved with it, from 0.424/0.649/0.892 to 0.468/0.741/0.914. Five
    /// separate places went on quoting the old triple — including §7.2 of the
    /// spec, *inside the parenthetical explaining that this had already
    /// happened once*. Editing five numbers by hand has the same half-life as
    /// the last time; this is the thing that fails.
    ///
    /// Bands are deliberately tight — ±0.02, about 3% — because the point is to
    /// catch a coefficient change, not to tolerate one. If this fails after a
    /// deliberate refit, requote every site listed in its message and *then*
    /// move the band. It is not a range to widen.
    #[test]
    fn gel_band_holds() {
        let mut extents: Vec<f64> = (0..24_u64)
            .map(|seed| {
                let mut rng = Stream::new(FUSION_SEED ^ seed);
                let c = draw_consts(&mut rng);
                let w: Vec<(u8, f64)> = derive_table(c, c.n_elements)
                    .iter()
                    .map(|e| (e.valence, e.abundance))
                    .collect();
                gel_extent(&w).p_c
            })
            .collect();
        assert!(
            extents.iter().all(|p| p.is_finite()),
            "non-finite gel extent"
        );
        #[expect(
            clippy::disallowed_methods,
            reason = "§13.4: finiteness asserted on the line above, which is the \
                      precondition that makes total_cmp safe"
        )]
        extents.sort_by(f64::total_cmp);
        // `extents[12]` is the upper order statistic of 24, not the mean of the
        // two central values. Labelled precisely because every figure in this
        // file gets requoted verbatim elsewhere.
        let (min, median, max) = (extents[0], extents[12], extents[23]);
        for (name, got, want) in [
            ("min", min, 0.468),
            ("median", median, 0.741),
            ("max", max, 0.914),
        ] {
            assert!(
                (got - want).abs() < 0.02,
                "gel {name} moved: {got:.3} vs {want:.3}. If this was a deliberate \
                 refit, requote it in prd.md §7.2, and in v0.md at the Q4 comment \
                 and the Step 13 commit message, BEFORE touching this band."
            );
        }
    }

    /// The clamp in `lateral_made` is an exemption, so pin its blast radius.
    ///
    /// It exists because the two halves of the "identity" use different
    /// coordination conventions at `outer == 1` (see that function's doc). If
    /// it ever fires anywhere else, the conventions have diverged somewhere new
    /// and the identity claim is wrong in a way no other test would show —
    /// `contacts_upto` would quietly gain spurious *made* contacts while every
    /// valence still read plausibly.
    #[test]
    fn the_clamp_fires_only_at_a_lone_outer_site() {
        for k in 6..=14_usize {
            for shell in 1..=4_usize {
                let cap = shell_size(k, shell);
                for outer in 0..=cap {
                    let fires = lateral_made_raw(cap, outer) < 0.0;
                    assert_eq!(
                        fires,
                        outer == 1,
                        "k={k} shell={shell} cap={cap} outer={outer}: clamp fired={fires}, \
                         raw={}",
                        lateral_made_raw(cap, outer)
                    );
                }
            }
        }
    }

    /// The clamp is load-bearing, and this is the mutation that proves it.
    ///
    /// `a_two_unit_cluster_has_exactly_one_contact` passes only because of it;
    /// unclamped, `contacts_upto(c, 2)` is 0.250 / 0.500 / 0.625 at
    /// k = 6 / 10 / 14 against a 1e-9 tolerance. This asserts that directly, so
    /// the dependency is recorded here rather than discovered by someone
    /// "simplifying" the clamp away.
    #[test]
    fn without_the_clamp_the_two_unit_cluster_test_would_fail() {
        for k in [6_usize, 10, 14] {
            let cap = shell_size(k, 1);
            assert!(
                lateral_made_raw(cap, 1) < 0.0,
                "k={k}: raw lateral_made at outer=1 should be negative, got {}",
                lateral_made_raw(cap, 1)
            );
        }
    }

    /// **Exhaustive, and asserting the attained set — both deliberate.**
    ///
    /// `valence` reads only `cap` and `outer`, and both are functions of `k`
    /// and `units` alone: the other seven `Consts` fields never enter. So the
    /// ceiling is a pure function of `(k, n_elements)` — 9 x 61 = 549 cells,
    /// about a millisecond. Sampling 48 seeds out of an enumerable space buys
    /// a probability where a certainty is free.
    ///
    /// It is not academic. `frontier_matches_exact_counts_on_a_real_sphere`
    /// admits `FRONTIER_COEFF` in [6.135, 7.215]; at 7.20 — *inside the band
    /// its own gate allows* — the true ceiling reaches 7 in 7 of the 549 cells,
    /// and 48 random draws miss all seven **54% of the time**. The previous
    /// sampled form was therefore still blind on the high side, which is the
    /// only side that matters: §7.1 puts valence at 0..6, so a 7 is a spec
    /// breach and this is the one test positioned to see it. Raising 48 to 480
    /// would not fix it — that shrinks a probability instead of removing it.
    ///
    /// `assert_eq!` on the set, not `.all(..)` containment: a refit that
    /// collapsed every universe to 5 satisfies containment and fails this.
    #[test]
    fn valence_ceiling_is_four_to_six_over_the_drawn_range() {
        let mut seen = std::collections::BTreeSet::new();
        for k in 6..=14_usize {
            let c = consts_with_k(k);
            for n in 60..=120_usize {
                if let Some(m) = derive_table(c, n).iter().map(|e| e.valence).max() {
                    seen.insert(m);
                }
            }
        }
        let expected: std::collections::BTreeSet<u8> = [4, 5, 6].into_iter().collect();
        assert_eq!(
            seen, expected,
            "the attained valence ceiling set moved over the full drawn grid"
        );
    }

    /// The ceiling *is* removed in principle — the mechanism has no per-site
    /// cap. This pins that separately from what a drawn table reaches, so the
    /// two claims cannot be confused again.
    #[test]
    fn the_frontier_count_has_no_per_site_ceiling() {
        let c = consts_with_k(14);
        let max = derive_table(c, 400).iter().map(|e| e.valence).max();
        assert!(
            max.is_some_and(|m| m >= 8),
            "a count over sites should exceed the per-site cap of 6 eventually: max {max:?}"
        );
    }

    /// §7.1: mass increases down the table. The first draft's series stopped
    /// increasing at N = 124, so the third closed shell weighed less than the
    /// element below it — checked against a 120-element table, that passed.
    #[test]
    fn mass_increases_with_unit_count() {
        let c = consts();
        for pair in derive_table(c, 300).windows(2) {
            assert!(
                pair[1].mass > pair[0].mass,
                "N={} mass {} >= N={} mass {}",
                pair[0].units,
                pair[0].mass,
                pair[1].units,
                pair[1].mass
            );
        }
    }

    /// The claim that makes `Mass`'s fixed-point representation the right home
    /// for this series. It was false for every element of the first draft,
    /// under a comment asserting it was true.
    #[test]
    fn every_mass_is_exact_on_the_fixed_point_grid() {
        for seed in 0..24_u64 {
            let mut rng = Stream::new(FUSION_SEED ^ seed);
            let c = draw_consts(&mut rng);
            for e in derive_table(c, c.n_elements) {
                let scaled = e.mass * 1024.0;
                assert!(
                    (scaled - scaled.round()).abs() < 1e-9,
                    "seed {seed} N={}: mass {} is not a multiple of 1/1024",
                    e.units,
                    e.mass
                );
            }
        }
    }

    /// Affinity is `2 * exposed_fraction - 1`, so it cannot leave `(-1, 1]`
    /// unless burial is miscounted. The first draft's version reached exactly
    /// +1.000 at every closure by reading the capacity of the *empty* shell.
    #[test]
    fn affinity_stays_in_range() {
        let c = consts();
        for e in derive_table(c, 200) {
            assert!(
                e.affinity > -1.0 && e.affinity <= 1.0,
                "N={}: affinity {} outside (-1, 1]",
                e.units,
                e.affinity
            );
        }
    }

    /// Shell closures must agree between the fill loop and the closed form,
    /// since valence, burial and the period index all key on them.
    #[test]
    fn closures_agree_with_the_shell_arithmetic() {
        for k in 6..=14_usize {
            let c = consts_with_k(k);
            let from_table: Vec<usize> = derive_table(c, 120)
                .iter()
                .filter(|e| e.outer == 0)
                .map(|e| e.units)
                .collect();
            assert_eq!(from_table, closures(k, 120), "k={k}");
        }
    }

    /// Different `k` must give different table *shapes*, or the drawn shell law
    /// buys nothing over the fixed `10n^2 + 2` it replaced.
    #[test]
    fn drawn_k_gives_distinct_closure_sets() {
        let sets: std::collections::BTreeSet<Vec<usize>> =
            (6..=14_usize).map(|k| closures(k, 120)).collect();
        assert!(sets.len() >= 8, "only {} distinct closure sets", sets.len());
    }

    /// The binding peak sits at a closure by *mechanism*, and this pins the
    /// mechanism rather than the coincidence.
    ///
    /// A closed shell has no frontier, so every lateral contact it can make is
    /// made, and contacts-per-unit is locally maximal there. If a later change
    /// makes the energy series smooth across closures, this fails — which is
    /// the point, because the first draft's minimum was the distance from a
    /// typed-in constant and looked identical from outside.
    #[test]
    fn energy_per_unit_is_locally_maximal_at_every_closure() {
        let c = consts();
        let t = derive_table(c, 200);
        for w in t.windows(3) {
            if w[1].outer == 0 && w[1].units > 1 {
                assert!(
                    w[1].energy_per_unit > w[0].energy_per_unit
                        && w[1].energy_per_unit > w[2].energy_per_unit,
                    "N={} is a closure but not a local energy maximum: {} / {} / {}",
                    w[1].units,
                    w[0].energy_per_unit,
                    w[1].energy_per_unit,
                    w[2].energy_per_unit
                );
            }
        }
    }

    /// The peak must be a function of the drawn constants, or `peak` has
    /// merely moved from an argument to a hard-coded location.
    #[test]
    fn the_peak_moves_with_the_generated_constants() {
        let mut peaks = std::collections::BTreeSet::new();
        for seed in 0..24_u64 {
            let mut rng = Stream::new(FUSION_SEED ^ seed);
            let c = draw_consts(&mut rng);
            peaks.insert(peak_of(&derive_table(c, c.n_elements)));
        }
        assert!(
            peaks.len() >= 6,
            "peak took only {} distinct values over 24 universes: {peaks:?}",
            peaks.len()
        );
    }
}
