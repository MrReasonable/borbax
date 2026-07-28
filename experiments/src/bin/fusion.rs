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
//! Run: `cargo run -p borbax-experiments --release --bin fusion`

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
/// `k = 10` recovers `10n^2 + 2` — the Mackay icosahedral numbers, which are
/// **real**, observed in rare-gas cluster mass spectra. That is one member of
/// the family and it is a coincidence of `k`, not a construction. Most drawn
/// `k` have no real counterpart at all, which is the G3 argument: the family is
/// invented and one member happens to coincide.
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
    // `sqrt` is native on purpose: IEEE-754 specifies it exactly, so it is
    // portable without `det_math`. See the criterion in `clippy.toml`.
    let continuum = FRONTIER_COEFF * (cap as f64 * f * (1.0 - f)).sqrt();

    // **The continuum form overestimates tiny patches, and the bound has to be
    // the *compact* one.** A patch of `a` sites laid out compactly — which is
    // what a spherical-cap fill order means — exposes at most `6a` less twice
    // Harborth's internal-bond count. An earlier version used a bare `6a`, the
    // bound for a patch whose sites touch *nothing*, and it returned valence 3
    // where the exact count gives 2 at `outer = 3`. By symmetry the same bound
    // applies to the empty side.
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
/// and 7.316 at D = 162: finite-size effects pull it down at small caps, where
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
/// 0.25, admitting [6.2, 7.2] — plus or minus 7%.
const FRONTIER_COEFF: f64 = 6.90;

/// Maximum unmade lateral contacts a *compact* patch of `a` sites can expose.
///
/// Harborth's result gives the **maximum** internal bonds over all arrangements
/// of `a` cells, `floor(3a - sqrt(12a - 3))`, attained by the compact ones — so
/// it is an upper bound on internal bonds and a lower bound on exposure. An
/// earlier version said "at least", which is true only for the compact case and
/// reads as a universal bound.
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
fn lateral_made(cap: usize, outer: usize) -> f64 {
    let raw = (lateral_coordination(cap) * outer as f64 - unmade_lateral(cap, outer)) * 0.5;
    if raw > 0.0 { raw } else { 0.0 }
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
        // strain that grows as `n^2`. Summed over shells that is `~N^(5/3)`
        // total, `~N^(2/3)` per unit — which competes with the contact term's
        // saturation and produces a maximum whose position is a
        // function of the generated constants and of nothing typed in. (An earlier
        // version quoted a closed form `eps*z*coeff/(8*sigma)` with `coeff`
        // undefined; solving for it across the grid gives 0.80 to 5.00, so no such
        // constant exists and the claim was unfalsifiable as written.)
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
fn gel_extent(weights: &[(u8, f64)]) -> (f64, f64, f64) {
    let total: f64 = weights.iter().map(|&(_, w)| w).sum();
    let mean: f64 = weights.iter().map(|&(f, w)| f64::from(f) * w / total).sum();
    let mean_sq: f64 = weights
        .iter()
        .map(|&(f, w)| f64::from(f) * f64::from(f) * w / total)
        .sum();
    if mean <= 0.0 {
        return (0.0, f64::INFINITY, f64::INFINITY);
    }
    let f_w = mean_sq / mean;
    if f_w <= 1.0 {
        return (f_w, f64::INFINITY, f64::INFINITY);
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
    (f_w, p_c, p_c * mean * 0.5)
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
        let mut rng = Stream::new(0xF0_51_00 ^ seed);
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
        let (_, p_c, _) = gel_extent(&weights);
        gel_extents.push(p_c);

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
                let mut r = Stream::new(0xF0_51_00 ^ seed);
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
    let mut rng = Stream::new(0xF0_51_00);
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
    let (fw_flat, pc_flat, bc_flat) = gel_extent(&flat);
    let (fw_w, pc_w, bc_w) = gel_extent(&weighted);
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
        let (fw, pc, bc) = gel_extent(&w);
        println!("{decay:>7.2} {fw:>9.3} {pc:>9.3} {bc:>11.3} {low:>19.3}");
    }

    // Sorted once, then read at three positions, because `f64::min`/`max` are
    // disallowed (their tie-breaking differs by target). This is a bare
    // `Vec<f64>` with **no** ID tie-break — an earlier version of this comment
    // claimed one, which would have been cited later as precedent for a sort
    // that genuinely needs it. It is deterministic for a different reason:
    // `sort_by` is stable and the input is already in seed order.
    let mut sorted = gel_extents.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    println!(
        "\ngel extents across the 24 universes: min {:.3}  median {:.3}  max {:.3}",
        sorted.first().copied().unwrap_or(0.0),
        sorted.get(sorted.len() / 2).copied().unwrap_or(0.0),
        sorted.last().copied().unwrap_or(0.0)
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
            idx.sort_by(|&a, &b| {
                g.dirs[b][2]
                    .partial_cmp(&g.dirs[a][2])
                    .unwrap_or(core::cmp::Ordering::Equal)
                    .then(a.cmp(&b))
            });
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
        let mut rng = Stream::new(0xF0_51_00);
        draw_consts(&mut rng)
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
            let Ok(g) = built else { return (0.0, 0) };
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
            let c = Consts { k, ..consts() };
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
            let c = Consts {
                k,
                z: (k + 2) as f64,
                ..consts()
            };
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
            let mut rng = Stream::new(0xF0_51_00 ^ k as u64);
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
            let c = Consts { k, ..consts() };
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
    /// is **4 to 7** — measured, after this assertion was first written as
    /// `5..=7` from a guess and failed with a 4 in it. An earlier version of this test asserted `max >= 8` and
    /// passed by measuring at 300 elements, outside anything the scheme builds
    /// — the same "measured at a size the plan does not generate" defect that
    /// put two false claims into the plan's prose.
    ///
    /// §7.2 wants a thin high-valence tail. What it gets is a ceiling of 6 at
    /// best and 4 at worst — and §7.2 asked for a *shape*, not a magnitude:
    /// abundance at valence >= 3 spans 0.077..0.566 across universes, so the
    /// shape it asked for is present in every one.
    /// Whether that is enough functionality is Task 5's battery to answer; it
    /// is not something this file may assert either way.
    #[test]
    fn valence_ceiling_is_four_to_six_over_the_drawn_range() {
        let mut seen = std::collections::BTreeSet::new();
        for seed in 0..48_u64 {
            let mut rng = Stream::new(0xF0_51_00 ^ seed);
            let c = draw_consts(&mut rng);
            if let Some(m) = derive_table(c, c.n_elements)
                .iter()
                .map(|e| e.valence)
                .max()
            {
                seen.insert(m);
            }
        }
        assert!(
            seen.iter().all(|&v| (4..=6).contains(&v)),
            "valence ceiling left 4..=6 over the drawn range: {seen:?}"
        );
    }

    /// The ceiling *is* removed in principle — the mechanism has no per-site
    /// cap. This pins that separately from what a drawn table reaches, so the
    /// two claims cannot be confused again.
    #[test]
    fn the_frontier_count_has_no_per_site_ceiling() {
        let c = Consts { k: 14, ..consts() };
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
            let mut rng = Stream::new(0xF0_51_00 ^ seed);
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
            let c = Consts { k, ..consts() };
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
            let mut rng = Stream::new(0xF0_51_00 ^ seed);
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
