//! How `N` base units pack into shells, and what that packing costs in
//! contacts (spec §7.1).
//!
//! Every element property in [`crate::element`] is a function of what this
//! module computes, which is the point: an element is not a row of drawn
//! numbers, it is a cluster of `N` units and its properties are consequences of
//! the cluster's shape.

/// Units in the `n`th shell: `k*n^2 + 2`. Only `k` is drawn.
///
/// **The `+ 2` is the Euler characteristic, not "the antipodal pair".** For any
/// convex simplicial polyhedron shell `n` holds `V + E(n-1) + F(n-1)(n-2)/2`;
/// apply `3F = 2E` and the linear term vanishes, and `V - E + F = 2` leaves
/// `(F/2)n^2 + 2`. So `k = F/2 = E/3 = V - 2`, and **`shell(1) = k + 2 = V`**,
/// which is what forces `z = k + 2` in [`PackingConsts`].
///
/// The tetrahedron falsifies the antipodal reading outright: `k = 2`,
/// shell 1 = 4 = 2*1^2 + 2, and it has no antipodal vertex pair at all. The
/// distinction changes no number here; it changes what the next person writes
/// when they generalise to a hemisphere or a torus.
///
/// `k` is how densely the base unit's own shape lets neighbours tile that
/// surface, and it is what makes two universes differ in table *shape*.
///
/// **Defined for `n >= 1`.** `shell_size(k, 0)` returns 2, which is
/// meaningless: shell 0 is the single central unit, and every caller
/// special-cases it.
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
/// **Do not reuse a geodesic's directions as site positions.** The shared
/// object is the icosahedron's face tiling and rotation group, not the point
/// set: geodesic is `10*4^L + 2` = 12, 42, 162, 642 against packing's
/// `k*n^2 + 2`. At `k = 10` they coincide only where `n = 2^L`, so 92 has no
/// geodesic rung at all and 162 is a different point set (mismatch 2.29e-2, 8%
/// of the neighbour spacing). Right at 12 and 42, silently wrong at 162.
#[must_use]
pub const fn shell_size(k: usize, n: usize) -> usize {
    k * n * n + 2
}

/// Lateral coordination inside a shell — Euler's formula, not a tunable.
///
/// A triangulated sphere on `v` vertices has exactly `3v - 6` edges, so its
/// mean degree is `6 - 12/v`. The twelve five-coordinate sites Euler forces are
/// the same twelve wherever they sit.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "`cap` is a shell capacity, bounded by 14*4^2+2 = 226 over the drawn range"
)]
pub fn lateral_coordination(cap: usize) -> f64 {
    6.0 - 12.0 / cap as f64
}

/// Unmade lateral contacts on the frontier of a partly-filled shell.
///
/// The canonical filling order is a compact spherical cap growing from one
/// pole, so the frontier is a circle of latitude. A cap covering fraction `f`
/// subtends `cos(theta) = 1 - 2f`, so its boundary circle has radius
/// `2*sqrt(f(1-f))`; site spacing on a shell of `cap` sites is
/// `sqrt(4*pi/cap)`; so the boundary carries `2*sqrt(pi*cap*f(1-f))` sites,
/// each pointing **two** of its six contacts outward — not three, which is one
/// of the three errors the fitted [`FRONTIER_COEFF`] exists to absorb. The
/// coefficient is fitted to exact counts on real geodesic spheres, not derived
/// from this sketch; see its own doc comment before changing either.
///
/// **The discrete bound is not tidying, and it binds in exactly one place.**
/// At `outer = 1` the continuum form returns 6.46 (cap 8) to 6.87 (cap 128) for
/// a site whose shell offers 4.5 to 5.9. (An earlier version said "10.5
/// contacts for a site that has six" — both numbers were measured at the naive
/// 10.635 coefficient and are stale by two refits.) By symmetry the same holds
/// for the `cap - a` empty sites, so the smaller side binds.
///
/// Measured, the discrete branch is taken iff `compact_bound(a) < COEFF*sqrt(a)`,
/// which at 6.90 is true only at `a = 1`: **51 of the 1080 `(cap, outer)` cells
/// the table reaches** (`k` in 6..=14, `units` 1..=120), all at
/// `min(outer, cap-outer) == 1`.
///
/// The figure this replaces, "54 of 1287", was a **splice of two different
/// enumerations** and neither produces it: 54 is the bind count over
/// `k x shells 1..3 x outer 0..=cap` (1341 cells), while 1287 is the cell count
/// at `units 1..=143`, where the count is 53. The `shells 1..4` sweep the tests
/// use is 72 of 2808. Quote an enumeration with every figure, or the next
/// reader re-splices it.
///
/// **It does not by itself stop `lateral_made` going negative** — an earlier
/// version of this comment implied it did. `compact_bound(1) = 6` exceeds
/// `lateral_coordination(cap)` for every cap in the table, so the raw value is
/// negative at every `outer == 1` cell and it is [`lateral_made`]'s own clamp
/// that returns 0. See that function.
///
/// **There is deliberately no `if outer == 0 { return 0.0 }` here.** It would
/// be redundant — both factors already vanish at a closure — and a redundant
/// special case on a chemically-loaded predicate is how a detector quietly
/// becomes a declaration two refactors later.
///
/// The property that no branch reads the closure predicate is checked at
/// *source* level by `xtask`, not behaviourally: deletion is uninformative here
/// because two independent factors hold the zero, and continuity was shown
/// defeatable. [`crate::element`]'s continuity assertion is a backstop, not the
/// guarantee.
///
/// **Which factor returns the zero at a closure is not a detail.** At
/// `outer == 0` and `outer == cap` both the continuum and the discrete bound are
/// `+0.0`, and `continuum < discrete` is false — so the value returned is
/// [`compact_bound`] at `a == 0`'s literal, the branch `xtask` exempts by name. The
/// continuum vanishes too, but it is never the value. Measured: poisoning
/// `compact_bound(0)` to return `NaN` makes `unmade_lateral(cap, 0)` return
/// `NaN` at every cap, while a continuum that does *not* vanish at `f = 0`
/// leaves it at exactly `0.0`. Anyone reasoning about the zero must reason about
/// `compact_bound`; `the_zero_at_a_closure_comes_from_the_discrete_bound` is
/// what stops the two being confused again.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "`cap` and `outer` are bounded shell counts, at most 226 over the drawn range"
)]
pub fn unmade_lateral(cap: usize, outer: usize) -> f64 {
    let f = outer as f64 / cap as f64;
    let continuum = continuum_at(cap, f);
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
/// so it can be evaluated off-lattice by
/// `the_continuum_frontier_vanishes_as_a_square_root`.
///
/// **`cap as f64 * f * (1.0 - f)` associates left and must stay written that
/// way.** Rust emits no fast-math flags and no FMA contraction, so no compiler
/// will move this — but a reader will. Measured: rewriting it as
/// `cap as f64 * (f * (1.0 - f))` moves **315 of the 2808 `(cap, outer)` cells**
/// the table reaches and **18 of 1080** `contacts_upto` values, by 1 ulp. `peak`
/// is a strict-`>` argmax over `energy_per_unit` ([`crate::element`]), so one ulp
/// can relocate it by a whole element and change every `decay_rate` below it.
/// `the_split_is_bit_exact_on_the_integer_domain` pins it.
///
/// `sqrt` stays native: IEEE-754 specifies it exactly, so it is portable without
/// `det_math`. See the criterion in `clippy.toml`.
#[must_use]
#[inline]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "`cap` is a bounded shell count, at most 226 over the drawn range"
)]
pub fn continuum_at(cap: usize, f: f64) -> f64 {
    FRONTIER_COEFF * (cap as f64 * f * (1.0 - f)).sqrt()
}

/// Unmade lateral contacts on the frontier, divided by the lateral
/// coordination — so the unit is "one absent neighbour's worth of surface" and
/// the number counts docking notches.
///
/// **This exists so `valence` has exactly one spelling.** The plan's Step 11
/// inlined the division at the call site while Step 9's exhaustive
/// `valence_ceiling_at` helper called this function, so the enumeration that is
/// supposed to be the whole input space for the valence ceiling would have been
/// measuring a *different expression* from the one the table uses. The probe
/// has always had this function and calls it at its own valence site.
#[must_use]
pub fn frontier_notches(cap: usize, outer: usize) -> f64 {
    unmade_lateral(cap, outer) / lateral_coordination(cap)
}

/// Frontier-length coefficient, **fitted to exact counts on real triangulated
/// spheres** rather than derived.
///
/// **Cap-dependent, so a single constant is a compromise.** Least-squares
/// against exact counts gives 6.075 at D = 12, 7.078 at D = 42 and 7.306 at
/// D = 162 — finite-size effects pull it down at small caps, where Euler's
/// twelve pentagons are a large fraction of the shell. This table's caps run
/// 8..130 (`shell_size(8, 4)`, reached at `units = 120`), and D = 12 and D = 42 are real shells (`shell_size(10,1)`,
/// `shell_size(10,2)`); D = 162 is not one of them. 6.90 minimises worst-case
/// error across all three.
///
/// The naive derivation gives `6*sqrt(pi) = 10.635` and is wrong by 1.46x, for
/// **three** reasons, not the two an earlier version listed: a boundary site
/// points *two* of six contacts outward, not three; the spacing must carry the
/// triangular lattice's `(sqrt(3)/2)s^2` area-per-site; and a lattice boundary
/// zig-zags through ~10% more sites than `L/s`. The first two alone give 6.598,
/// so a reader re-deriving from them would "fix" the constant downward. The
/// derivation that reproduces the measurement is Cauchy's line-crossing
/// formula: `16*sqrt(3)*sqrt(sqrt(3)/(8*pi))` = 7.2751, asymptotic.
///
/// `frontier_matches_exact_counts_on_a_real_sphere` holds this to the geometry
/// at **all three levels** at 0.25. The gate edges, bisected, are
/// **[6.132767, 7.216878]**; `[6.2, 7.2]` is that band rounded inward, and it is
/// the only form quoted anywhere so the two-numbers-for-one-quantity problem
/// does not arise. An earlier version
/// gated D = 162 alone and admitted -22%/+29% — a band over which the valence
/// series moves in up to 23 of 24 universes.
pub const FRONTIER_COEFF: f64 = 6.90;

/// **Minimum** unmade lateral contacts a patch of `a` sites can expose.
///
/// Harborth's result gives the **maximum** internal bonds over all arrangements
/// of `a` cells, `floor(3a - sqrt(12a - 3))`, attained by the compact ones — an
/// upper bound on internal bonds and therefore a **lower** bound on exposure,
/// with equality only for maximally compact patches. "At least ... so it
/// exposes at most" inverts it; a cap-order fill is not maximally compact, and
/// on a 162-vertex geodesic `a = 3` exposes 14 against this bound's 12.
/// Verified by brute force over every connected animal for `a = 1..10`.
///
/// A bare `6a` is the *dispersed* bound and contradicts the compactness the
/// continuum assumes. Its recorded symptom — "valence 3 where the exact count
/// gives 2 at `a = 3`" — was measured at the naive 10.635 coefficient and is
/// unreachable at 6.90.
///
/// The result is provably positive (`>= 2*sqrt(12a-3)`, since `floor(x) <= x`;
/// brute-forced to a = 20000, minimum 6.0 at a = 1), so **no clamp is
/// written** — a guard with no reachable case reads as evidence the case
/// exists. An earlier version of this block wrote one anyway, four lines under
/// this sentence.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "`a` is at most a shell capacity, 226 over the drawn range"
)]
pub fn compact_bound(a: usize) -> f64 {
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
/// every internal contact twice and every frontier contact once. At a full
/// shell the frontier is empty and this returns `3*cap - 6` exactly, which is
/// the check that the two halves agree.
///
/// **The identity holds everywhere except `outer == 1`, where the clamp below
/// is load-bearing — name it, do not delete it.** The two halves use different
/// coordination conventions: `compact_bound(1) = 6` is a flat-lattice degree,
/// `lateral_coordination(cap) = 6 - 12/cap` is the sphere's, and the second is
/// strictly smaller for every cap. So the raw value is negative at every
/// `outer == 1` cell (36 of 2808 over k = 6..14, shells 1..4; worst −0.75 at
/// cap = 8) and the clamp is what returns the correct 0.
///
/// `a_two_unit_cluster_has_exactly_one_contact` — the test that caught the
/// `z = 6 + k/2` defect — passes *only* because of it: unclamped,
/// `contacts_upto(c, 2)` is 0.250 / 0.500 / 0.625 at k = 6 / 10 / 14 against a
/// 1e-9 tolerance. `the_clamp_fires_only_at_a_lone_outer_site` pins that, so a
/// later "simplify the clamp away" fails a test instead of silently adding
/// spurious made contacts while every valence still reads plausibly.
#[must_use]
pub fn lateral_made(cap: usize, outer: usize) -> f64 {
    let raw = lateral_made_raw(cap, outer);
    if raw > 0.0 { raw } else { 0.0 }
}

/// [`lateral_made`] before the clamp, so a test can assert *where* the clamp
/// fires rather than leaving it an unnamed rescue inside a function whose doc
/// calls itself an identity.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "`outer` is at most a shell capacity, 226 over the drawn range"
)]
pub fn lateral_made_raw(cap: usize, outer: usize) -> f64 {
    (lateral_coordination(cap) * outer as f64 - unmade_lateral(cap, outer)) * 0.5
}

/// The packing constants of one universe.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackingConsts {
    k: usize,
    z: f64,
}

impl PackingConsts {
    /// The only constructor, and the fields are private, because `z` is not an
    /// independent parameter: `shell_size(k, 1) == k + 2` is the set of units
    /// touching the core, which is what `z` means. Written out by hand at three
    /// sites in an earlier draft, correct at all three by luck of copy-paste,
    /// with a doc comment saying "follows `k`, and must" above two `pub` fields
    /// that let `PackingConsts { k: 6, z: 40.0 }` compile.
    ///
    /// `z` feeds [`contacts_upto`], which sets the binding peak, which sets
    /// `decay_rate` and `abundance`. A fourth call site transposing the
    /// constant produces a plausible table.
    #[must_use]
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "`shell_size(k, 1)` is `k + 2`, at most 16 over the drawn range"
    )]
    pub const fn new(k: usize) -> Self {
        Self {
            k,
            z: shell_size(k, 1) as f64,
        }
    }

    /// Sites per unit area of shell — the drawn part of the shell law.
    #[must_use]
    pub const fn k(self) -> usize {
        self.k
    }

    /// Bulk coordination of the base unit.
    #[must_use]
    pub const fn z(self) -> f64 {
        self.z
    }
}

/// Contacts made by a cluster of `units`, accumulated shell by shell.
///
/// **A continuum count does not work and looks fine while failing.** The
/// obvious form — `z*n/2` less a surface deficit growing as `n^(2/3)` — has
/// the deficit exceeding `n` itself for everything below `n = 113`, so it
/// collapses to `z*n/4`, makes contacts-per-unit a *constant*, and removes the
/// binding peak entirely. Measured: the peak sat at `N = 1` in all 24 drawn
/// universes, reported as "lands on a closure in 24 of 24" because 1 is the
/// first closure.
///
/// So: each shell's sites make inward contacts to the shell below and lateral
/// contacts to their own neighbours, and the inward budget is whatever outward
/// capacity the shell below has left after its own inward and lateral contacts
/// are accounted for. That recursion is what makes contacts-per-unit *rise*
/// toward `z/2`, which is what a binding peak needs.
///
/// The inward budget is divided by the shell's *capacity*, not by how many
/// sites happen to be placed: how many contacts a site makes inward is fixed by
/// the geometry, not by how full its shell currently is. A review proposed
/// dividing by `take` instead — that was measured against the old
/// `z = 6 + k/2`; with `z = k + 2` forced by the shell law,
/// `cap = shell_size(k, 1) = z`, so this is exactly 1.0, and dividing by `take`
/// would give a lone unit the core's *entire* outward budget.
///
/// Accumulation is in shell order and must stay that way (§13.4).
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "`prev_cap`, `cap` and `take` are bounded shell counts; `units` is at most 120"
)]
pub fn contacts_upto(c: PackingConsts, units: usize) -> f64 {
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

/// Cumulative shell closures inside a table of `n` elements.
///
/// `1` is always the first: a single unit is a closed cluster with no frontier.
#[must_use]
pub fn closures(k: usize, n: usize) -> Vec<usize> {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Ground truth: a real triangulated sphere, built here and answerable to
    /// nothing but Euler.
    ///
    /// **Why this is a local copy rather than an import.** The probe's geodesic
    /// lives in `borbax-experiments` and Task 7's canonical one lands in
    /// `borbax-molecule` — which sits *downstream* of this crate in the crate
    /// order, so `borbax-universe` can never depend on it, now or later. The
    /// alternative was a dev-dependency on a probe crate that Step 8d already
    /// contemplates retiring, which would take this test with it.
    ///
    /// **A duplicate is partly self-certifying, and the boundary of "partly" is
    /// the thing to know.** Any triangulated sphere has exactly twelve degree-5
    /// vertices and no other degree, so [`degrees`] failing is proof the
    /// *triangulation* is wrong. That much needs no cross-check.
    ///
    /// It does **not** certify the fill order, and that gap was live in the
    /// first version of this module: without the canonical sort in [`geodesic`]
    /// the sphere is still correct and `degrees` still passes, while
    /// [`cap_order`]'s index tie-break silently follows the subdivision loop.
    /// Measured against the probe, that moved the admitted `FRONTIER_COEFF`
    /// band by about 1.1 in opposite directions at the two ends. With the sort
    /// in place the two copies agree to four decimals on the worst-case error
    /// at every level. See [`geodesic`].
    #[expect(
        clippy::indexing_slicing,
        reason = "fixed-size vertex arrays indexed by their own 0..len loop bounds"
    )]
    mod exact {
        type Vec3 = [f64; 3];

        fn dot(a: Vec3, b: Vec3) -> f64 {
            a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        }

        fn dist2(a: Vec3, b: Vec3) -> f64 {
            let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            dot(d, d)
        }

        fn normalise(a: Vec3) -> Vec3 {
            let n = dot(a, a).sqrt();
            if n > 0.0 {
                [a[0] / n, a[1] / n, a[2] / n]
            } else {
                a
            }
        }

        /// Midpoint of two directions, projected back onto the unit sphere.
        ///
        /// Added in the caller's order: IEEE addition is commutative, so the
        /// shared edge between two faces interns to a bit-identical vertex
        /// whichever face reaches it first.
        fn midpoint(a: Vec3, b: Vec3) -> Vec3 {
            normalise([
                f64::midpoint(a[0], b[0]),
                f64::midpoint(a[1], b[1]),
                f64::midpoint(a[2], b[2]),
            ])
        }

        /// The twelve icosahedron vertices, from the golden ratio — the cyclic
        /// permutations of `(0, ±1, ±φ)`.
        pub(super) fn icosahedron() -> Vec<Vec3> {
            #[expect(
                clippy::manual_midpoint,
                reason = "this is the closed form of the golden ratio; the rewrite hides \
                          the one thing a reader needs to recognise here"
            )]
            let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
            let mut v = Vec::with_capacity(12);
            for &s1 in &[1.0, -1.0] {
                for &s2 in &[1.0, -1.0] {
                    v.push([0.0, s1, s2 * phi]);
                    v.push([s1, s2 * phi, 0.0]);
                    v.push([s1 * phi, 0.0, s2]);
                }
            }
            v.into_iter().map(normalise).collect()
        }

        /// Squared length of the shortest edge — the adjacency threshold.
        ///
        /// A comparison rather than `f64::min`, which is banned under §13.1:
        /// std may return either input when they compare equal, measured
        /// differing by architecture.
        fn min_pair_dist2(v: &[Vec3]) -> f64 {
            let mut min_d2 = f64::MAX;
            for i in 0..v.len() {
                for j in (i + 1)..v.len() {
                    let d2 = dist2(v[i], v[j]);
                    if d2 < min_d2 {
                        min_d2 = d2;
                    }
                }
            }
            min_d2
        }

        fn faces_of(v: &[Vec3]) -> Vec<[usize; 3]> {
            let min_d2 = min_pair_dist2(v);
            let adjacent = |a: usize, b: usize| (dist2(v[a], v[b]) - min_d2).abs() < 1e-9;
            let mut out = Vec::new();
            for i in 0..v.len() {
                for j in (i + 1)..v.len() {
                    for k in (j + 1)..v.len() {
                        if adjacent(i, j) && adjacent(j, k) && adjacent(i, k) {
                            out.push([i, j, k]);
                        }
                    }
                }
            }
            out
        }

        /// Find `p` in `verts`, or append it. A linear scan avoids any hashing
        /// whose iteration order could vary.
        fn intern(verts: &mut Vec<Vec3>, p: Vec3) -> usize {
            for (i, q) in verts.iter().enumerate() {
                if dist2(*q, p) < 1e-18 {
                    return i;
                }
            }
            verts.push(p);
            verts.len() - 1
        }

        /// Subdivided-icosahedron directions: 12 at level 0, 42 at 1, 162 at 2.
        ///
        /// **The canonical `(z, y, x)` sort at the end is load-bearing, and
        /// omitting it is a defect this test cannot see.** [`cap_order`] breaks
        /// ties in z-height by *index*, and a geodesic has many vertices at
        /// equal z (whole rings of latitude). Without this sort that tie-break
        /// inherits the order the subdivision loop happened to intern vertices
        /// in — arbitrary — so the fill order, and therefore every exact count
        /// below, becomes a property of the construction rather than of the
        /// geometry.
        ///
        /// Measured: with the sort omitted, `FRONTIER_COEFF = 6.20` fails at
        /// D = 162 / outer = 80 while the probe passes it, and 7.30 passes here
        /// while the probe fails at D = 12 / outer = 6. Same spheres, different
        /// fills. `degrees` does **not** catch this — it certifies the
        /// triangulation, and an arbitrarily-ordered triangulation is still a
        /// correct one.
        pub(super) fn geodesic(level: u32) -> Vec<Vec3> {
            geodesic_from(&icosahedron(), level)
        }

        /// [`geodesic`] from an explicit base, so a test can permute it and
        /// check that the result depends on the point set and not on the order
        /// the subdivision loop happened to visit it in.
        pub(super) fn geodesic_from(base: &[Vec3], level: u32) -> Vec<Vec3> {
            let mut verts = base.to_vec();
            let mut faces = faces_of(base);
            for _ in 0..level {
                let mut next = Vec::with_capacity(faces.len() * 4);
                for f in &faces {
                    let (p0, p1, p2) = (verts[f[0]], verts[f[1]], verts[f[2]]);
                    let (m01, m12, m20) = (midpoint(p0, p1), midpoint(p1, p2), midpoint(p2, p0));
                    let m = [
                        intern(&mut verts, m01),
                        intern(&mut verts, m12),
                        intern(&mut verts, m20),
                    ];
                    next.push([f[0], m[0], m[2]]);
                    next.push([m[0], f[1], m[1]]);
                    next.push([m[2], m[1], f[2]]);
                    next.push([m[0], m[1], m[2]]);
                }
                faces = next;
            }
            assert!(
                verts.iter().all(|v| v.iter().all(|c| c.is_finite())),
                "geodesic coordinates must be finite"
            );
            #[expect(
                clippy::disallowed_methods,
                reason = "§13.4: finiteness asserted above, and every component is a \
                          normalised coordinate built from sums of squares, so no -0.0 \
                          arises; ordering all three coordinates makes the result a \
                          function of the point set rather than of the build loop"
            )]
            verts.sort_by(|a, b| {
                a[2].total_cmp(&b[2])
                    .then(a[1].total_cmp(&b[1]))
                    .then(a[0].total_cmp(&b[0]))
            });
            verts
        }

        /// Adjacency by nearest-neighbour threshold.
        ///
        /// A geodesic's edges are its shortest vertex-vertex distances; the next
        /// shell of distances sits far enough above that a 1.3x threshold on the
        /// minimum separates them cleanly. [`degrees`] is what verifies that.
        pub(super) fn adjacency(dirs: &[Vec3]) -> Vec<Vec<usize>> {
            let cut = min_pair_dist2(dirs) * 1.3 * 1.3;
            (0..dirs.len())
                .map(|i| {
                    (0..dirs.len())
                        .filter(|&j| j != i && dist2(dirs[i], dirs[j]) < cut)
                        .collect()
                })
                .collect()
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

        /// Fill order: a compact spherical cap growing from one pole. This is
        /// the *canonical filling order* the continuum derivation assumes, so
        /// the ground truth has to use it too or the comparison is meaningless.
        pub(super) fn cap_order(dirs: &[Vec3]) -> Vec<usize> {
            let mut idx: Vec<usize> = (0..dirs.len()).collect();
            // `partial_cmp(..).unwrap_or(Equal)` is not a total order under
            // NaN, so the precondition is checked rather than assumed. These
            // are coordinates of unit vectors, finite by construction —
            // asserted anyway, because "by construction" is what the last
            // several defects in this project were also confident of.
            assert!(
                dirs.iter().all(|d| d[2].is_finite()),
                "geodesic directions must be finite"
            );
            #[expect(
                clippy::disallowed_methods,
                reason = "§13.4: NaN is excluded by the assert above, and the index \
                          tie-break keeps the order total across equal heights — which \
                          are common on a geodesic"
            )]
            idx.sort_by(|&a, &b| dirs[b][2].total_cmp(&dirs[a][2]).then(a.cmp(&b)));
            idx
        }

        /// Unmade lateral contacts after filling the first `a` sites of
        /// `order`: the number of filled-empty adjacent pairs. Exact, by
        /// construction.
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

    /// **The specification for `unmade_lateral`, and the test whose absence let
    /// a 1.46x error stand for two rounds.**
    ///
    /// Builds a real geodesic, derives its adjacency, asserts the degree
    /// histogram is Euler's twelve pentagons (which is what proves the
    /// threshold right), fills in spherical-cap order — the same canonical
    /// order the continuum assumes — and counts filled-to-empty edges exactly.
    ///
    /// Gates **D = 12, 42 and 162** at 0.25. All three, because the coefficient
    /// is cap-dependent and D = 162 is the one level this table never uses.
    /// Nothing else in this file can catch a drift in [`FRONTIER_COEFF`],
    /// because every other test compares the formula against itself.
    ///
    /// **One margin worth knowing, measured rather than assumed.** A break that
    /// displaces every fill fraction by half a site — `f -> f + 0.5/cap` — scores
    /// 0.2449 against this 0.25 gate. It passes, by 2% of tolerance, on the only
    /// test in the crate that compares against real geometry.
    /// `the_continuum_frontier_vanishes_as_a_square_root` catches it outright, so
    /// the pair covers it; neither leg does alone, and this gate must not be
    /// tightened to try (the healthy continuum branch already sits at 0.1951).
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "exact contact counts are small integers being compared against a float formula"
    )]
    #[test]
    fn frontier_matches_exact_counts_on_a_real_sphere() {
        fn check(level: u32, d: usize) -> (f64, usize) {
            let dirs = exact::geodesic(level);
            assert_eq!(
                dirs.len(),
                d,
                "geodesic level {level} produced {} vertices",
                dirs.len()
            );
            let adj = exact::adjacency(&dirs);

            // Euler: any triangulated sphere has exactly twelve degree-5
            // vertices. If this fails the adjacency threshold is wrong and
            // every count below is meaningless.
            let (d5, _, other) = exact::degrees(&adj);
            assert_eq!((d5, other), (12, 0), "D={d}: degree histogram is wrong");

            let order = exact::cap_order(&dirs);
            let (mut worst, mut worst_at) = (0.0_f64, 0);
            for a in 1..d {
                let truth = exact::unmade(&adj, &order, a) as f64;
                if truth > 0.0 {
                    let rel = (unmade_lateral(d, a) - truth).abs() / truth;
                    if rel > worst {
                        worst = rel;
                        worst_at = a;
                    }
                }
            }
            (worst, worst_at)
        }

        for (d, (worst, at)) in [(12, check(0, 12)), (42, check(1, 42)), (162, check(2, 162))] {
            assert!(
                worst < 0.25,
                "D={d}: worst relative error {worst:.4} at outer={at} \
                 — the formula has drifted from the geometry it approximates"
            );
        }
    }

    /// The canonical sort in [`exact::geodesic`] is load-bearing and, until this
    /// test, was protected by nothing — neutralising the `sort_by` left all 28
    /// crate tests passing. That is CLAUDE.md's "present, believed, inert" shape,
    /// committed in the same change that documented why the sort mattered.
    ///
    /// **The guard is the sort's own stated purpose**, tested directly: the
    /// result must be a function of the point set, not of the build loop. So
    /// rebuild from a permuted base icosahedron and compare bit patterns.
    ///
    /// Deliberately independent of [`FRONTIER_COEFF`]. Tightening the 0.25 gate
    /// in `frontier_matches_exact_counts_on_a_real_sphere` would also catch the
    /// unsorted variant, but only by re-coupling this guarantee to the constant
    /// it exists to hold steady — and at the shipped 6.90 the unsorted variant
    /// passes that gate at all three levels anyway (0.2000 / 0.1401 / 0.1869),
    /// which is precisely why nothing caught it.
    #[expect(
        clippy::indexing_slicing,
        reason = "`k` ranges over 0..3 and both operands are `[f64; 3]`"
    )]
    #[test]
    fn the_canonical_sort_makes_the_sphere_a_function_of_the_point_set() {
        let base = exact::icosahedron();
        for level in 0..=2 {
            let reference = exact::geodesic_from(&base, level);
            for shift in 1..base.len() {
                let mut permuted = base.clone();
                permuted.rotate_left(shift);
                let got = exact::geodesic_from(&permuted, level);
                assert_eq!(
                    got.len(),
                    reference.len(),
                    "level {level} shift {shift}: vertex count"
                );
                for (i, (a, b)) in got.iter().zip(&reference).enumerate() {
                    for k in 0..3 {
                        assert_eq!(
                            a[k].to_bits(),
                            b[k].to_bits(),
                            "level {level} shift {shift}: vertex {i} coordinate {k} \
                             depends on the order the base was visited in"
                        );
                    }
                }
            }
        }
    }

    /// **What is observable about the zero at a closure — and what is not.**
    ///
    /// At `outer == 0` and `outer == cap` the continuum and the discrete bound
    /// are *both* `+0.0`, so `continuum < discrete` is false and the value
    /// returned comes from [`compact_bound`] at `a == 0` — the branch `xtask` exempts by
    /// name. Verified by poisoning: making `compact_bound(0)` return `NaN` makes
    /// `unmade_lateral(cap, 0)` return `NaN` at every cap, while a continuum that
    /// does *not* vanish at `f = 0` still leaves it exactly `0.0`.
    ///
    /// **That provenance is not testable from behaviour, and the first draft of
    /// this test pretended otherwise.** It asserted
    /// `unmade_lateral(cap, 0).to_bits() == compact_bound(0).to_bits()`, which
    /// holds tautologically — poison `compact_bound` and both sides move
    /// together, `NaN == NaN` bitwise, green. Two expressions that agree over the
    /// domain the test covers is the exact defect this project keeps recording.
    /// Because the two factors *coincide* at a closure, no behavioural assertion
    /// can separate them there. That is the whole reason the guarantee lives in
    /// a source check.
    ///
    /// So this asserts the three things that *are* observable and do
    /// discriminate: each factor vanishes exactly, and where the two genuinely
    /// differ — `outer == 1` — the discrete arm is the one selected.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "`cap` is a bounded shell count"
    )]
    #[test]
    fn each_factor_vanishes_exactly_and_the_discrete_arm_wins_at_a_lone_site() {
        for k in 6..=14 {
            for shell in 1..=4 {
                let cap = shell_size(k, shell);
                // Each factor vanishes on its own terms, bit-exactly.
                assert_eq!(
                    continuum_at(cap, 0.0).to_bits(),
                    0.0_f64.to_bits(),
                    "cap={cap}"
                );
                assert_eq!(
                    continuum_at(cap, 1.0).to_bits(),
                    0.0_f64.to_bits(),
                    "cap={cap}"
                );
                assert_eq!(compact_bound(0).to_bits(), 0.0_f64.to_bits());
                assert_eq!(
                    unmade_lateral(cap, 0).to_bits(),
                    0.0_f64.to_bits(),
                    "cap={cap}"
                );
                assert_eq!(
                    unmade_lateral(cap, cap).to_bits(),
                    0.0_f64.to_bits(),
                    "cap={cap}"
                );
                // Where they differ, the smaller must win and it must be the
                // discrete one — this is the cell that pins the `min`'s direction.
                let (c, d) = (continuum_at(cap, 1.0 / cap as f64), compact_bound(1));
                assert!(
                    c > d,
                    "cap={cap}: continuum {c} no longer exceeds the bound {d}"
                );
                assert_eq!(unmade_lateral(cap, 1).to_bits(), d.to_bits(), "cap={cap}");
            }
        }
    }

    /// The split of [`continuum_at`] out of [`unmade_lateral`] must not move a
    /// single bit on the integer domain.
    ///
    /// The reference is the pre-split expression written inline and verbatim.
    /// The hazard is not the compiler — Rust emits no fast-math and no FMA
    /// contraction — it is a reader tidying `cap * f * (1-f)` into
    /// `cap * (f * (1-f))`, which moves 315 of these 2808 cells by 1 ulp.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "`cap` and `outer` are bounded shell counts"
    )]
    #[test]
    fn the_split_is_bit_exact_on_the_integer_domain() {
        for k in 6..=14 {
            for shell in 1..=4 {
                let cap = shell_size(k, shell);
                for outer in 0..=cap {
                    let f = outer as f64 / cap as f64;
                    let reference = FRONTIER_COEFF * (cap as f64 * f * (1.0 - f)).sqrt();
                    assert_eq!(
                        continuum_at(cap, f).to_bits(),
                        reference.to_bits(),
                        "k={k} cap={cap} outer={outer}"
                    );
                }
            }
        }
    }

    /// The continuum's `f -> 0` asymptotics, and **not** a declaration guard.
    ///
    /// **Named for what it covers.** It was proposed as the replacement for the
    /// `xtask` source check, on the reasoning that a point-declaration fires at
    /// one lattice point and an off-lattice limit would step around it. Measured,
    /// that does not hold: `if outer == 0 { return 0.0 }` passes this, and so
    /// does `if f == 0.0 { return 0.0 }` — the sweep runs strictly between 0 and
    /// `1/cap` and neither branch fires inside it. A *redundant* declaration is
    /// observationally equivalent by definition, so no behavioural test can
    /// reach it; that is a theorem, not a shortcoming of this one.
    ///
    /// What it does cover is otherwise uncovered, which is why it is here: four
    /// arithmetic breaks of the continuum slip past
    /// `frontier_matches_exact_counts_on_a_real_sphere`'s 0.25 gate and fail
    /// here — an additive floor at +0.001 and +0.01, an exponent of 0.45 instead
    /// of 0.5, and a half-site shift in `f`. And that gate cannot be tightened to
    /// take over: the continuum branch already sits at 0.1951 against 0.25.
    ///
    /// The ratio is analytically `2*sqrt((1-f)/(1-f/4))`, free of both `cap` and
    /// [`FRONTIER_COEFF`], so sweeping caps repeats one check rather than
    /// widening coverage. Swept anyway, cheaply; claim no coverage from it.
    ///
    /// Note the sweep sits entirely below `outer = 1`, where the discrete bound
    /// is what `unmade_lateral` actually returns — so this tests an
    /// extrapolation of a branch the shipped function does not use there. That
    /// is deliberate and is the reason it is not the guard.
    #[test]
    fn the_continuum_frontier_vanishes_as_a_square_root() {
        for k in 6..=14 {
            for shell in 1..=4 {
                let cap = shell_size(k, shell);
                let mut f = 1e-2;
                for _ in 0..5 {
                    let ratio = continuum_at(cap, f) / continuum_at(cap, f / 4.0);
                    assert!(
                        (ratio - 2.0).abs() < 0.01,
                        "k={k} cap={cap} f={f:e}: ratio {ratio} is not the square-root law"
                    );
                    f /= 100.0;
                }
                // The tail must be tight, or an additive floor survives.
                let tail = continuum_at(cap, 1e-10) / continuum_at(cap, 2.5e-11);
                assert!(
                    (tail - 2.0).abs() < 5e-4,
                    "k={k} cap={cap}: tail ratio {tail} — the continuum does not vanish as sqrt(f)"
                );
            }
        }
    }

    /// **Vacuous, and kept only as a note.** `lateral_made(cap, cap)` reduces to
    /// `(z_lat*cap)/2 = 3cap - 6` for *any* formula whose unmade count vanishes
    /// at a full shell — which every candidate does. It tests one endpoint and
    /// cannot detect an error in the frontier formula, including the 1.46x one.
    /// The identity it names is sound; it is the *test* that is empty.
    /// `frontier_matches_exact_counts_on_a_real_sphere` is what this was trying
    /// to be.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "shell capacities are small bounded counts"
    )]
    #[test]
    fn a_complete_shell_makes_euler_many_lateral_contacts() {
        for k in 6..=14 {
            for n in 1..=3 {
                let cap = shell_size(k, n);
                let made = lateral_made(cap, cap);
                let euler = 3.0 * cap as f64 - 6.0;
                assert!(
                    (made - euler).abs() < 1e-9,
                    "k={k} n={n} cap={cap}: {made} lateral contacts, Euler says {euler}"
                );
            }
        }
    }

    /// The continuum boundary formula over-counts a lone site, and without the
    /// discrete bound the contact count goes negative — silently.
    ///
    /// (At `FRONTIER_COEFF = 6.90` the over-count is 6.46 at cap 8 to 6.87 at
    /// cap 128. An earlier version of this comment said "10.5 contacts for a
    /// site that has only six"; that was measured at the naive 10.635 and is
    /// stale by two refits. The *assertion* is unaffected.)
    #[test]
    fn a_lone_site_exposes_no_more_than_its_own_coordination() {
        for k in 6..=14 {
            let cap = shell_size(k, 2);
            assert!(
                unmade_lateral(cap, 1) <= 6.0 + 1e-9,
                "k={k}: a lone site exposed {} contacts",
                unmade_lateral(cap, 1)
            );
            assert!(lateral_made(cap, 1) >= 0.0);
        }
    }

    /// `lateral_made`'s clamp is an exemption; pin where it fires.
    ///
    /// The "identity" holds everywhere except `outer == 1`, where the two
    /// halves use different coordination conventions and the clamp supplies the
    /// correct 0. If it fires anywhere else the conventions have diverged
    /// somewhere new, and nothing else would show it — `contacts_upto` would
    /// gain spurious *made* contacts while every valence still read plausibly.
    #[test]
    fn the_clamp_fires_only_at_a_lone_outer_site() {
        for k in 6..=14 {
            for shell in 1..=4 {
                let cap = shell_size(k, shell);
                for outer in 0..=cap {
                    assert_eq!(
                        lateral_made_raw(cap, outer) < 0.0,
                        outer == 1,
                        "k={k} cap={cap} outer={outer}: raw={}",
                        lateral_made_raw(cap, outer)
                    );
                }
            }
        }
    }

    /// The output bound Step 8c's `compact_bound` exemption needs.
    ///
    /// Exempting a function by name in a source check bounds nothing about what
    /// it returns. `NaN` and negatives escape the `min` in `unmade_lateral`, and
    /// both saturate to **0** through `round_ties_even() as u8` — so
    /// `closed_shells_have_valence_zero` reads a correct zero while
    /// `contacts_upto` is poisoned. This is the assertion that makes the
    /// exemption's blast radius provable instead of argued.
    #[test]
    fn unmade_lateral_is_finite_and_non_negative_everywhere() {
        for k in 6..=14 {
            for shell in 1..=4 {
                let cap = shell_size(k, shell);
                for outer in 0..=cap {
                    let u = unmade_lateral(cap, outer);
                    assert!(
                        u.is_finite() && u >= 0.0,
                        "k={k} cap={cap} outer={outer}: unmade_lateral = {u}"
                    );
                }
            }
        }
    }

    /// Contacts per unit must *rise* toward `z/2`. A continuum estimate makes
    /// it constant, which removes the binding peak entirely — and the
    /// resulting table looks perfectly reasonable on a printout.
    #[test]
    fn contacts_per_unit_rises_with_cluster_size() {
        let c = PackingConsts::new(10);
        let (a, b, d) = (
            contacts_upto(c, 13) / 13.0,
            contacts_upto(c, 55) / 55.0,
            contacts_upto(c, 147) / 147.0,
        );
        assert!(
            a < b && b < d,
            "contacts per unit did not rise: {a} {b} {d}"
        );
        assert!(d < c.z() / 2.0, "contacts per unit passed the bulk bound");
    }

    /// Two units share exactly one contact, whatever the shell law says.
    #[test]
    fn a_two_unit_cluster_has_exactly_one_contact() {
        for k in 6..=14 {
            let c = PackingConsts::new(k);
            let n = contacts_upto(c, 2);
            // `1e-9`, not `0.35`. Two units share exactly one contact — that is
            // not an approximation. The old tolerance was exactly wide enough to
            // span the defect it should have caught: with the falsified
            // `z = 6 + k/2` this value runs 0.8125..1.1250 across k, all inside
            // 0.35. Built through `new`, it is 1.0 for every k.
            assert!(
                (n - 1.0).abs() < 1e-9,
                "k={k}: a 2-unit cluster made {n} contacts"
            );
        }
    }

    /// Different `k` must give different table *shapes*, or the drawn shell
    /// law buys nothing over the fixed `10n^2 + 2` it replaced.
    #[test]
    fn drawn_k_gives_distinct_closure_sets() {
        let sets: std::collections::BTreeSet<Vec<usize>> =
            (6..=14).map(|k| closures(k, 120)).collect();
        assert_eq!(sets.len(), 9, "closure sets collapsed: {sets:?}");
    }
}
