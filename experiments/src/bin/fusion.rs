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
use core::f64::consts::PI;

/// Units in the `n`th packing shell: `k*n^2 + 2`.
///
/// **Only `k` is drawn; the rest is geometry.** A shell is a surface, and a
/// surface grows as the square of its radius, so `n^2` is not a choice. The
/// `+ 2` is the antipodal pair every shell carries. `k` is how densely the base
/// unit's own shape lets neighbours tile that surface, which is exactly the
/// kind of thing a universe seed should decide and the kind of thing no real
/// table fixes for us.
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

/// Unmade lateral contacts on the frontier of a partly-filled shell, divided
/// by the lateral coordination — so the unit is "one absent neighbour's worth
/// of surface", and the number is a count of docking notches.
///
/// **Derived, and every factor has a reason.** The canonical filling order is a
/// compact spherical cap growing from one pole, so the frontier is a circle of
/// latitude. A cap covering fraction `f` of a sphere subtends `cos(theta) =
/// 1 - 2f`, so its boundary circle has radius `sin(theta) = 2*sqrt(f(1-f))` in
/// units of the sphere radius. Site spacing on a shell of `cap` sites is
/// `sqrt(4*pi/cap)` in the same units, so the boundary carries
/// `2*pi*2*sqrt(f(1-f)) / sqrt(4*pi/cap) = 2*sqrt(pi*cap*f(1-f))` sites. Each
/// boundary site on a triangulated surface points three of its six contacts
/// outward into empty space. Divide the resulting `6*sqrt(pi*cap*f(1-f))` by
/// the lateral coordination and the sixes cancel.
///
/// **It is exactly 0 at a closure, and *two* independent factors put it there.**
/// The continuum's `f(1-f)` vanishes at both ends of a shell, and so does the
/// discrete bound's `min(outer, cap - outer)`. Probed: mutating either one
/// alone leaves `closed_shells_have_valence_zero_with_no_branch` passing,
/// because the other still zeroes it.
///
/// That makes the deletion test — `ptable.rs`'s argument for this property —
/// **uninformative here**, and saying "delete the branch and nothing changes"
/// would be claiming evidence this file cannot supply. What can be shown is
/// continuity: a hardcoded `if closed { 0 }` sitting on top of a formula that
/// returns something else at a closure produces a *jump*, and
/// `valence_is_continuous_across_closures` measures the jump. That is the test
/// that distinguishes a mechanism from a declaration, and it is the one Task 4
/// should carry.
fn unmade_lateral(cap: usize, outer: usize) -> f64 {
    // **There is deliberately no `if outer == 0 { return 0.0 }` here, and a
    // mutation probe is why.** One stood at the top of this function. With it
    // in place, rewriting the continuum below to drop its `(1 - f)` factor —
    // which is the whole reason a closed shell has no frontier — left
    // `closed_shells_have_valence_zero_with_no_branch` *passing*. The branch
    // was doing the work the formula was being credited for, which is the
    // `if closed { 0 }` defect from `ptable.rs` wearing a different hat.
    //
    // It was also redundant: at `outer == 0` the continuum is 0 and `smaller`
    // is 0, and at `outer == cap` both are 0 again. Deleting it changes no
    // value and makes the zero a property of the arithmetic.
    let f = outer as f64 / cap as f64;
    // `sqrt` is native on purpose: IEEE-754 specifies it exactly, so it is
    // portable without `det_math`. See the criterion in `clippy.toml`.
    let continuum = 6.0 * (PI * cap as f64 * f * (1.0 - f)).sqrt();

    // **The continuum form is wrong for tiny patches, and the first run of this
    // probe is what said so.** At `outer = 1` it returns 10.5 contacts for a
    // single site that has only six. A patch of `a` sites cannot expose more
    // than `6a` lateral contacts — that is the case where none of them touch
    // each other — and by symmetry the same bound applies to the `cap - a`
    // empty sites. Whichever side is smaller binds.
    let smaller = if outer < cap - outer {
        outer
    } else {
        cap - outer
    };
    let discrete = 6.0 * smaller as f64;
    if continuum < discrete {
        continuum
    } else {
        discrete
    }
}

/// Unmade lateral contacts on the frontier, divided by the lateral
/// coordination — so the unit is "one absent neighbour's worth of surface",
/// and the number counts docking notches.
fn frontier_notches(cap: usize, outer: usize) -> f64 {
    unmade_lateral(cap, outer) / lateral_coordination(cap)
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
    /// Mass of one base unit. Dyadic, so the mass series lands on `Mass`'s grid.
    base_mass: f64,
    /// Mass carried away per made contact. Dyadic for the same reason.
    contact_defect: f64,
    base_radius: f64,
    /// Energy per made contact.
    eps: f64,
    /// Radial strain coefficient. The destabiliser.
    sigma: f64,
    /// How steeply abundance falls with unit count. The gel lever (§7.2).
    decay: f64,
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
        // Coordination follows the tiling density: a denser shell means more
        // neighbours in contact. Not independent of `k`, and must not be.
        z: 6.0 + (k as f64) * 0.5,
        // Dyadic: 1, 1.25, 1.5 or 1.75. See the mass computation.
        base_mass: 1.0 + 0.25 * rng.below(4) as f64,
        // Dyadic: 1/512 .. 4/512.
        contact_defect: (1 + rng.below(4)) as f64 / 512.0,
        base_radius: 0.30 + 0.02 * rng.below(11) as f64,
        eps: 0.8 + 0.05 * rng.below(9) as f64,
        sigma: 0.02 + 0.01 * rng.below(12) as f64,
        decay: 0.04 + 0.01 * rng.below(13) as f64,
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
        let defect = (raw_defect * 1024.0).round() / 1024.0;
        let mass = units as f64 * c.base_mass - defect;

        // Radius: the cube root of the unit count, times a bulge for the
        // partly-filled outer shell sitting proud of the packed core.
        //
        // G3: real atomic radius *decreases* across a period and jumps at a new
        // one. This does the opposite — grows through a shell, drops at
        // closure. Structurally anti-isomorphic, not a renamed trend.
        let fill = outer as f64 / cap as f64;
        let radius = c.base_radius * det_math::cbrt(units as f64) * (1.0 + 0.55 * fill);

        // Valence: docking notches on the frontier. Zero at closure with no
        // branch — see `frontier_notches`.
        let valence = frontier_notches(cap, outer).round() as u8;

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
        // saturation and produces a maximum whose position is
        // `eps*z*coeff/(8*sigma)`, a function of the generated constants and of
        // nothing typed in.
        let strain = c.sigma * det_math::cbrt(units as f64 * units as f64);
        let energy_per_unit = c.eps * contacts / units as f64 - strain;

        // Abundance: fusion builds heavy clusters from light ones, so each
        // extra unit costs a step and abundance falls geometrically. Clusters
        // that bind better survive their surroundings longer, so the per-unit
        // energy tilts it. Both halves are processes already in the model.
        let abundance = det_math::exp(-(units as f64) * c.decay + energy_per_unit * 0.6);

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

/// Weighted functionality `f_w = <f^2>/<f>` and the Molloy-Reed gel extent.
///
/// §7.2 as amended: a giant component appears at bond-conversion extent
/// `p_c = 1/(f_w - 1)`. A *high* `p_c` is the safe direction — the beaker has
/// to be driven further before it gels. Terminators (`f <= 1`) are counted in
/// the averages because they are what stops chains, which is the whole reason
/// the ceiling exists.
fn gel_extent(weights: &[(u8, f64)]) -> (f64, f64) {
    let total: f64 = weights.iter().map(|&(_, w)| w).sum();
    let mean: f64 = weights.iter().map(|&(f, w)| f64::from(f) * w / total).sum();
    let mean_sq: f64 = weights
        .iter()
        .map(|&(f, w)| f64::from(f) * f64::from(f) * w / total)
        .sum();
    if mean <= 0.0 {
        return (0.0, f64::INFINITY);
    }
    let f_w = mean_sq / mean;
    if f_w <= 1.0 {
        return (f_w, f64::INFINITY);
    }
    (f_w, 1.0 / (f_w - 1.0))
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
        "{:>4} {:>3} {:>6} {:>8} {:>6} {:>6} {:>6} {:>8}",
        "seed", "k", "maxval", "closures", "peak", "p@clos", "massOK", "gel p_c"
    );
    let mut peaks_on_closure = 0_usize;
    let mut all_masses_exact = true;
    let mut max_valences = Vec::new();
    let mut gel_extents = Vec::new();
    for seed in 0..24_u64 {
        let mut rng = Stream::new(0xF0_51_00 ^ seed);
        let c = draw_consts(&mut rng);
        let table = derive_table(c, 120);

        let max_val = table.iter().map(|e| e.valence).max().unwrap_or(0);
        max_valences.push(max_val);

        let peak = peak_of(&table);
        let clos = closures(c.k, 120);
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
        let (_, p_c) = gel_extent(&weights);
        gel_extents.push(p_c);

        println!(
            "{seed:>4} {:>3} {max_val:>6} {:>8} {peak:>6} {:>6} {:>6} {p_c:>8.3}",
            c.k,
            clos.len(),
            if on_closure { "yes" } else { "no" },
            if exact { "yes" } else { "NO" }
        );
    }
    println!(
        "\npeak lands on a closure in {peaks_on_closure} of 24 -- and this is STRUCTURAL, \
         not luck.\n  A closed shell has no frontier, so it makes every lateral contact \
         available to it and\n  maximises contacts per unit; the energy series is a sawtooth \
         with a local maximum at\n  each closure. What the generated constants choose is \
         *which* closure wins, and that\n  does move: {} distinct peaks over the 24 draws.",
        {
            let mut ps: Vec<usize> = Vec::new();
            for seed in 0..24_u64 {
                let mut r = Stream::new(0xF0_51_00 ^ seed);
                let cc = draw_consts(&mut r);
                let p = peak_of(&derive_table(cc, 120));
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
    let table = derive_table(c, 120);
    println!(
        "universe: k={} z={:.1} eps={:.2} sigma={:.3}",
        c.k, c.z, c.eps, c.sigma
    );
    println!(
        "closures: {:?}   peak: N={}\n",
        closures(c.k, 120),
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
    let (fw_flat, pc_flat) = gel_extent(&flat);
    let (fw_w, pc_w) = gel_extent(&weighted);
    println!(
        "\nunweighted (every element equally common):  f_w = {fw_flat:.3}  gels at extent {pc_flat:.3}"
    );
    println!(
        "abundance-weighted:                          f_w = {fw_w:.3}  gels at extent {pc_w:.3}"
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
        "{:>7} {:>9} {:>9} {:>26}",
        "decay", "f_w", "gel p_c", "mass on valence <= 2"
    );
    for step in 0..13_u64 {
        let decay = 0.04 + 0.01 * step as f64;
        let t = derive_table(Consts { decay, ..c }, 120);
        let tot: f64 = t.iter().map(|e| e.abundance).sum();
        let low: f64 = t
            .iter()
            .filter(|e| e.valence <= 2)
            .map(|e| e.abundance)
            .sum::<f64>()
            / tot;
        let w: Vec<(u8, f64)> = t.iter().map(|e| (e.valence, e.abundance)).collect();
        let (fw, pc) = gel_extent(&w);
        println!("{decay:>7.2} {fw:>9.3} {pc:>9.3} {low:>25.3}");
    }

    // Sorted once, then read at three positions — `f64::min`/`max` are
    // disallowed (tie-breaking differs by target), and a sort with an explicit
    // tie-break is what §13.4 asks for anyway.
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

    /// The test that actually distinguishes a mechanism from a declaration.
    ///
    /// Zero-at-closure is held by two independent factors here, so mutating
    /// either leaves the zero intact and the deletion test proves nothing. A
    /// declaration does not hide from *continuity*, though: an `if closed { 0 }`
    /// laid over a formula that returns 5 at a closure makes valence jump by 5
    /// in one element, and this measures that.
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

    /// A per-*site* count caps at `z - 6`, which is why the first draft could
    /// not exceed 6. A count over sites has no such ceiling, and §7.2 needs the
    /// tail to reach 8 for a high-functionality species to exist at all.
    #[test]
    fn valence_reaches_eight_somewhere_in_a_large_shell() {
        // k = 14 gives a third shell of 128 sites; the frontier of a half-filled
        // shell that size carries well over eight notches.
        let c = Consts { k: 14, ..consts() };
        let max = derive_table(c, 300).iter().map(|e| e.valence).max();
        assert!(
            max.is_some_and(|m| m >= 8),
            "valence never reached 8: max {max:?}"
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
            for e in derive_table(c, 120) {
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
            peaks.insert(peak_of(&derive_table(c, 120)));
        }
        assert!(
            peaks.len() >= 6,
            "peak took only {} distinct values over 24 universes: {peaks:?}",
            peaks.len()
        );
    }
}
