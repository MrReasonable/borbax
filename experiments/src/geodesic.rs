//! Sample directions on a sphere, and the rotations that permute them
//! (spec §8.2).
//!
//! The directions are vertices of a subdivided icosahedron. That choice is
//! not aesthetic: the icosahedral rotation group has 60 elements and every
//! one maps the vertex set exactly onto itself, so **a rotation is a
//! permutation** — a lookup, not an interpolation. With arbitrary sample
//! directions, rotating a signature would mean resampling, and the resulting
//! error would land squarely on the comparison that decides whether two
//! molecules bind.
//!
//! Built once and passed by reference. Deliberately not a global: global
//! state is one more thing that can differ between a test and a run.
//!
//! Nothing here uses a transcendental. `sqrt` and the four arithmetic
//! operations are exactly specified by IEEE-754, so this file is portable by
//! construction rather than by routing through `det_math` (§13.4).
//!
//! This is Task 7's content, written in `experiments` first because the G2
//! measurement needs it. It is intended to lift into `borbax-molecule`
//! unchanged.

#![allow(
    clippy::indexing_slicing,
    reason = "fixed-size arrays indexed by loop bounds proven against D and N_ROTATIONS; \
              a get().ok_or() on every coordinate access would obscure the geometry"
)]

/// Order of the icosahedral rotation group. Reflections are excluded — see
/// spec §22.8, which makes Borbax chemistry handed on purpose.
pub const N_ROTATIONS: usize = 60;

/// Why a geodesic could not be built. Both variants are bugs rather than bad
/// input in normal use; they exist so the failure is loud instead of silent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum GeoError {
    /// D must be one of the geodesic ladder values: 12, 42, 162.
    #[error("unsupported resolution {0}: D must be 12, 42 or 162")]
    UnsupportedResolution(usize),
    /// A rotation failed to map the vertex set onto itself. This would mean
    /// the symmetry assumption above is broken, which is a bug, not input.
    #[error("a rotation did not map the sample directions onto themselves")]
    RotationNotClosed,
}

/// Vertices produced by subdividing an icosahedron `level` times.
#[must_use]
pub const fn vertex_count(level: u32) -> usize {
    // 10 * 4^level + 2
    10 * 4usize.pow(level) + 2
}

const fn level_for(d: usize) -> Option<u32> {
    match d {
        12 => Some(0),
        42 => Some(1),
        162 => Some(2),
        _ => None,
    }
}

/// Sample directions at resolution `D`, with the rotation and antipode tables
/// that make comparing two signatures a table lookup.
#[derive(Debug, Clone)]
pub struct Geodesic<const D: usize> {
    /// Unit sample directions, in a canonical order.
    pub dirs: [[f64; 3]; D],
    /// `perms[r][i] = j` means rotation `r` carries direction `i` to `j`.
    pub perms: [[u8; D]; N_ROTATIONS],
    /// `anti[i] = j` where `dirs[j] == -dirs[i]`.
    ///
    /// **Required by the binding kernel** (§8.3). Two bodies in contact touch
    /// along opposite directions, so `affinity` indexes its partner through
    /// `anti[perms[r][i]]`. Omitting it converts the search into the 60
    /// *improper* elements of the icosahedral group — reflections only, the
    /// opposite of §22.8 — and nothing fails to signal it.
    ///
    /// The vertex set is antipodally closed at every subdivision level, so
    /// this is total; `build_anti` errors rather than approximating if it
    /// ever is not.
    pub anti: [u8; D],
}

type Vec3 = [f64; 3];

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn normalise(a: Vec3) -> Vec3 {
    let n = dot(a, a).sqrt();
    if n > 0.0 { scale(a, 1.0 / n) } else { a }
}

fn dist2(a: Vec3, b: Vec3) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}

/// Midpoint of two directions, projected back onto the unit sphere.
///
/// Written so the two endpoints are added in the caller's order: IEEE addition
/// is commutative, so the shared edge between two faces interns to a
/// bit-identical vertex whichever face reaches it first.
fn midpoint(a: Vec3, b: Vec3) -> Vec3 {
    // `f64::midpoint` is `(a + b) / 2` guarded against overflow. Every input
    // here is a component of a unit vector, so the guard never engages and
    // this is the exact halving of an exact sum — no rounding beyond the
    // addition itself, and no dependence on which face reaches the edge first.
    normalise([
        f64::midpoint(a[0], b[0]),
        f64::midpoint(a[1], b[1]),
        f64::midpoint(a[2], b[2]),
    ])
}

/// The twelve icosahedron vertices, from the golden ratio.
///
/// These are the cyclic permutations of `(0, ±1, ±φ)`. Generating all three
/// cyclic forms inside the sign loops keeps antipodal pairs exact: `-1.0` and
/// `-φ` are exact negations in IEEE, and `normalise` divides both members of a
/// pair by the same norm, so `dirs[anti[i]]` is the bit-exact negation of
/// `dirs[i]` rather than merely close to it.
fn icosahedron() -> Vec<Vec3> {
    #[allow(
        clippy::manual_midpoint,
        reason = "this is the closed form of the golden ratio, and both operands are constants \
                  near 1 — the overflow `f64::midpoint` guards against cannot arise, and the \
                  rewrite hides the one thing a reader needs to recognise here"
    )]
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
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

/// The squared length of the shortest edge, which is the adjacency threshold.
fn min_pair_dist2(v: &[Vec3]) -> f64 {
    // `f64::min` rather than a comparison would be the obvious spelling and is
    // banned (§13.1): std documents that when the inputs compare equal — which
    // `+0.0` and `-0.0` do — either may be returned non-deterministically, and
    // it was measured returning different signs on aarch64 and x86-64 *and*
    // under different codegen on one target. Squared distances make that
    // unreachable here, but this function lifts into `borbax-molecule` and the
    // comparison costs nothing.
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

/// Faces derived rather than hardcoded: on an icosahedron, a face is exactly
/// a triple of mutually-adjacent vertices, and adjacency is "at the minimum
/// pairwise distance". Deriving avoids a 20-row table that is easy to get
/// subtly wrong and impossible to spot by reading.
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

/// Find `p` in `verts`, or append it. Linear scan is fine at D <= 162 and
/// avoids any hashing whose iteration order could vary.
fn intern(verts: &mut Vec<Vec3>, p: Vec3) -> usize {
    for (i, q) in verts.iter().enumerate() {
        if dist2(*q, p) < 1e-18 {
            return i;
        }
    }
    verts.push(p);
    verts.len() - 1
}

/// Index of the sample direction equal to `p`, or `RotationNotClosed`.
///
/// `tol2` is a squared distance, and it is deliberately tight. The point is
/// not to find the nearest direction — it is to assert that `p` *is* one. A
/// loose match would let a broken symmetry assumption through as a slightly
/// wrong answer, which is the failure mode this whole module exists to make
/// impossible.
fn exact_index<const D: usize>(dirs: &[Vec3; D], p: Vec3, tol2: f64) -> Result<u8, GeoError> {
    // Scanned in index order with a strict `<`, so the first of any tie wins
    // and the result does not depend on iteration order. At the tolerances
    // below a tie cannot happen without the vertex set being degenerate, but
    // pinning it costs nothing and removes the question.
    let mut best = (usize::MAX, f64::MAX);
    for (j, dir) in dirs.iter().enumerate() {
        let d = dist2(p, *dir);
        if d < best.1 {
            best = (j, d);
        }
    }
    if best.1 > tol2 {
        return Err(GeoError::RotationNotClosed);
    }
    // D <= 162 by construction, so this cannot truncate — but `try_from`
    // states that rather than asserting it in a comment next to an `as`.
    u8::try_from(best.0).map_err(|_| GeoError::RotationNotClosed)
}

impl<const D: usize> Geodesic<D> {
    /// Build the direction set and both lookup tables.
    ///
    /// # Errors
    ///
    /// [`GeoError::UnsupportedResolution`] if `D` is not 12, 42 or 162.
    /// [`GeoError::RotationNotClosed`] if a rotation or an antipode fails to
    /// land exactly on a sample direction, which would mean the icosahedral
    /// symmetry assumption in §8.2 is false.
    pub fn build() -> Result<Self, GeoError> {
        let level = level_for(D).ok_or(GeoError::UnsupportedResolution(D))?;

        let base = icosahedron();
        let mut verts = base.clone();
        let mut faces = faces_of(&base);

        for _ in 0..level {
            let mut next = Vec::with_capacity(faces.len() * 4);
            for f in &faces {
                // Midpoints are computed into locals first. Writing them
                // inline as arguments to `intern(&mut verts, ...)` borrows
                // `verts` immutably inside a call that already borrows it
                // mutably — E0502, since two-phase borrows do not cover
                // explicit `&mut` arguments.
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

        if verts.len() != D {
            return Err(GeoError::UnsupportedResolution(D));
        }

        // Canonical ordering, and the *reason* matters more than the call.
        //
        // `total_cmp` is a pure function of the bits, so it is portable when
        // the bits are — not, as an earlier version of this comment claimed,
        // free of platform variation in general. It orders on the sign bit,
        // and a runtime NaN's sign bit differs between aarch64 and x86-64;
        // that is precisely the defect `borbax_units::Span::canonical_cmp`
        // exists to fix, and this comment asserting the opposite in the file
        // CLAUDE.md says lifts into `borbax-molecule` unchanged is how the
        // wrong lesson gets cited a year from now.
        //
        // It is safe *here* because these coordinates are provably free of
        // both hazards, measured at all three levels: 0 `-0.0` components,
        // 0 NaN, antipodes bit-exact, max |‖v‖-1| = 1.11e-16. Every value is
        // a normalised coordinate built from sums of squares. When this lifts,
        // that reasoning lifts with it or the call changes (spec §13.4).
        #[expect(
            clippy::disallowed_methods,
            reason = "§13.4: measured free of NaN and -0.0 at all three levels — see above"
        )]
        verts.sort_by(|a, b| {
            a[2].total_cmp(&b[2])
                .then(a[1].total_cmp(&b[1]))
                .then(a[0].total_cmp(&b[0]))
        });

        let mut dirs = [[0.0; 3]; D];
        dirs.copy_from_slice(&verts);

        let perms = Self::build_perms(&dirs)?;
        let anti = Self::build_anti(&dirs)?;
        Ok(Self { dirs, perms, anti })
    }

    /// Index of `-dirs[i]` for every `i`.
    ///
    /// The tolerance is `1e-24` on a squared distance — effectively exact.
    /// The vertex set is built so antipodal pairs are bit-exact negations
    /// (see [`icosahedron`]), so anything short of that means the set is not
    /// antipodally closed, and a nearest match would corrupt every binding
    /// comparison rather than failing visibly.
    fn build_anti(dirs: &[Vec3; D]) -> Result<[u8; D], GeoError> {
        let mut anti = [0u8; D];
        for i in 0..D {
            let target = [-dirs[i][0], -dirs[i][1], -dirs[i][2]];
            anti[i] = exact_index(dirs, target, 1e-24)?;
        }
        Ok(anti)
    }

    /// Build the 60 rotations as permutations of `dirs`.
    ///
    /// Construction: pick a reference icosahedron vertex and one of its
    /// neighbours. For every (vertex, neighbour) pair there is exactly one
    /// rotation carrying the reference pair onto it, giving 12 x 5 = 60.
    ///
    /// Every frame is built right-handed (`e2 = e0 × e1`), so each `R` maps
    /// one right-handed orthonormal frame onto another and is therefore a
    /// proper rotation. That is what keeps reflections out of the set, and
    /// `every_rotation_has_determinant_plus_one` measures it rather than
    /// trusting this paragraph.
    fn build_perms(dirs: &[Vec3; D]) -> Result<[[u8; D]; N_ROTATIONS], GeoError> {
        let mats = rotation_matrices()?;
        let mut perms = [[0u8; D]; N_ROTATIONS];
        for (r, m) in mats.iter().enumerate() {
            for i in 0..D {
                // Looser than `build_anti`'s tolerance because `apply_mat`
                // accumulates three products per coordinate; still far
                // tighter than the ~0.3 spacing between neighbouring
                // directions at D = 162, so a genuine miss cannot pass.
                perms[r][i] = exact_index(dirs, apply_mat(m, dirs[i]), 1e-12)?;
            }
        }
        Ok(perms)
    }
}

/// A 3x3 matrix, indexed `[row][column]`.
pub type Mat3 = [[f64; 3]; 3];

/// Apply a matrix to a vector.
#[must_use]
pub fn apply_mat(m: &Mat3, p: Vec3) -> Vec3 {
    [
        m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2],
        m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2],
        m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2],
    ]
}

/// The 60 rotations as matrices, in the same order as [`Geodesic::perms`].
///
/// `perms` is what the binding kernel needs — a rotation as a table lookup.
/// This is the same 60 rotations as actual geometry, and it exists so the two
/// can be checked against each other. Without it, `perms` is a table of
/// integers that no test can tie back to rotating anything: it would pass
/// every group-theoretic check while corresponding to the wrong rotations, or
/// to none.
///
/// Construction: pick a reference icosahedron vertex and one of its
/// neighbours. For every (vertex, neighbour) pair there is exactly one
/// rotation carrying the reference pair onto it, giving 12 x 5 = 60.
///
/// Every frame is built right-handed (`e2 = e0 × e1`), so each `R` maps one
/// right-handed orthonormal frame onto another and is therefore proper. That
/// is what keeps reflections out of the set, and
/// `every_rotation_has_determinant_plus_one` measures it rather than trusting
/// this paragraph.
///
/// # Errors
///
/// [`GeoError::RotationNotClosed`] if the icosahedron does not yield exactly
/// 60 directed edges, which would mean `icosahedron` is not one.
pub fn rotation_matrices() -> Result<Vec<Mat3>, GeoError> {
    let base = icosahedron();
    let min_d2 = min_pair_dist2(&base);
    let neighbours = |a: usize| -> Vec<usize> {
        (0..base.len())
            .filter(|&b| b != a && (dist2(base[a], base[b]) - min_d2).abs() < 1e-9)
            .collect()
    };

    // Orthonormal frame from an axis and a second direction.
    let frame = |v: Vec3, w: Vec3| -> [Vec3; 3] {
        let e0 = normalise(v);
        let e1 = normalise(sub(w, scale(e0, dot(w, e0))));
        let e2 = cross(e0, e1);
        [e0, e1, e2]
    };

    let first_neighbour = *neighbours(0).first().ok_or(GeoError::RotationNotClosed)?;
    let f0 = frame(base[0], base[first_neighbour]);

    let mut out = Vec::with_capacity(N_ROTATIONS);
    for v in 0..base.len() {
        for &w in &neighbours(v) {
            let f = frame(base[v], base[w]);
            // R = F^T F0 in the sense that R maps f0[m] onto f[m]:
            // expressing p in F0's basis and rebuilding it in F's gives
            // R[k][l] = sum_m f[m][k] * f0[m][l].
            let mut r = [[0.0; 3]; 3];
            for (k, row) in r.iter_mut().enumerate() {
                for (l, cell) in row.iter_mut().enumerate() {
                    *cell = f[0][k] * f0[0][l] + f[1][k] * f0[1][l] + f[2][k] * f0[2][l];
                }
            }
            out.push(r);
        }
    }

    if out.len() != N_ROTATIONS {
        return Err(GeoError::RotationNotClosed);
    }
    Ok(out)
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::as_conversions,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds over fixed-size arrays; \
              every cast is u8 -> usize, a widening that cannot lose information — `as_conversions` \
              is here to catch truncation, and `usize::from` on each index would bury the assertion"
)]
mod tests {
    use super::*;

    fn norm(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    /// Scalar triple product. Equal to the determinant of the matrix whose
    /// rows are the three vectors, so its sign is the handedness of the triple.
    fn det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
    }

    /// The most linearly independent triple of directions, by |det|. Chosen by
    /// search rather than by hand so the test does not depend on the canonical
    /// ordering, and taking the first maximum makes the choice deterministic.
    fn best_triple(dirs: &[[f64; 3]]) -> (usize, usize, usize) {
        let mut best = (0, 1, 2);
        let mut best_abs = 0.0;
        for i in 0..dirs.len() {
            for j in (i + 1)..dirs.len() {
                for k in (j + 1)..dirs.len() {
                    let d = det3(dirs[i], dirs[j], dirs[k]).abs();
                    if d > best_abs {
                        best_abs = d;
                        best = (i, j, k);
                    }
                }
            }
        }
        assert!(best_abs > 0.5, "no well-conditioned triple: {best_abs}");
        best
    }

    #[test]
    fn vertex_counts_are_the_geodesic_ladder() {
        assert_eq!(vertex_count(0), 12);
        assert_eq!(vertex_count(1), 42);
        assert_eq!(vertex_count(2), 162);
    }

    #[test]
    fn directions_are_unit_vectors() {
        let g = Geodesic::<42>::build().unwrap();
        for d in &g.dirs {
            assert!((norm(*d) - 1.0).abs() < 1e-12, "not unit: {d:?}");
        }
    }

    #[test]
    fn construction_is_deterministic() {
        let a = Geodesic::<42>::build().unwrap();
        let b = Geodesic::<42>::build().unwrap();
        assert_eq!(a.dirs, b.dirs);
        assert_eq!(a.perms, b.perms);
        assert_eq!(a.anti, b.anti);
    }

    #[test]
    fn every_rotation_is_a_bijection() {
        let g = Geodesic::<42>::build().unwrap();
        for (r, perm) in g.perms.iter().enumerate() {
            let mut seen = [false; 42];
            for &j in perm {
                assert!(!seen[j as usize], "rotation {r} is not injective");
                seen[j as usize] = true;
            }
        }
    }

    #[test]
    fn rotations_are_all_distinct() {
        let g = Geodesic::<42>::build().unwrap();
        for i in 0..N_ROTATIONS {
            for j in (i + 1)..N_ROTATIONS {
                assert_ne!(g.perms[i], g.perms[j], "rotations {i} and {j} coincide");
            }
        }
    }

    /// The strongest correctness check available: the permutations must be
    /// closed under composition. If they are, they genuinely form the
    /// icosahedral rotation group, and §8.2's claim that rotation is an exact
    /// permutation of the sample directions holds.
    #[test]
    fn permutations_form_a_group() {
        let g = Geodesic::<42>::build().unwrap();
        let set: std::collections::BTreeSet<Vec<u8>> = g.perms.iter().map(|p| p.to_vec()).collect();
        assert_eq!(set.len(), N_ROTATIONS);
        for a in &g.perms {
            for b in &g.perms {
                let composed: Vec<u8> = (0..42).map(|i| b[a[i] as usize]).collect();
                assert!(set.contains(&composed), "not closed under composition");
            }
        }
    }

    #[test]
    fn identity_is_present() {
        let g = Geodesic::<42>::build().unwrap();
        let identity: Vec<u8> = (0..42u8).collect();
        assert!(g.perms.iter().any(|p| p.to_vec() == identity));
    }

    /// Every element must be a *proper* rotation, and this measures it
    /// directly rather than inferring it. For any linear map `R`,
    /// `det[Ra, Rb, Rc] = det(R) * det[a, b, c]`, so the ratio of the image
    /// triple product to the source one **is** `det(R)`. A rotation has
    /// determinant `+1`; a reflection has `-1`.
    ///
    /// This is the per-element form of the check CLAUDE.md warns about: a
    /// binding kernel built over improper elements searches reflections only,
    /// and nothing else in the suite fails.
    #[test]
    fn every_rotation_has_determinant_plus_one() {
        let g = Geodesic::<42>::build().unwrap();
        let (a, b, c) = best_triple(&g.dirs);
        let source = det3(g.dirs[a], g.dirs[b], g.dirs[c]);
        for (r, p) in g.perms.iter().enumerate() {
            let image = det3(
                g.dirs[p[a] as usize],
                g.dirs[p[b] as usize],
                g.dirs[p[c] as usize],
            );
            let det = image / source;
            assert!(
                (det - 1.0).abs() < 1e-9,
                "rotation {r} has determinant {det}, not +1 — it is a reflection"
            );
        }
    }

    /// `-I` sends every direction to its antipode, so *as a permutation* it is
    /// exactly `anti`. In three dimensions `det(-I) = -1`, so it is improper
    /// and must not appear among the 60.
    ///
    /// Stated separately from the determinant test because this is the
    /// specific failure CLAUDE.md names: because `-I` is not in the rotation
    /// group, indexing a partner through `perm[r][i]` instead of
    /// `anti[perm[r][i]]` silently converts the binding search into the 60
    /// improper elements. Both sets have 60 members, both are closed under the
    /// composition test above, and every other property here still holds.
    #[test]
    fn minus_identity_is_not_in_the_group() {
        let g = Geodesic::<42>::build().unwrap();
        assert!(
            !g.perms.iter().any(|p| p[..] == g.anti[..]),
            "-I is in the rotation set: the search is over reflections, not rotations"
        );
    }

    /// The antipode table is what makes `affinity` search rotations rather
    /// than reflections (§8.3, §22.8). These three properties are what it
    /// means for it to be correct.
    #[test]
    fn the_antipode_table_is_an_involution_on_opposite_directions() {
        let g = Geodesic::<42>::build().unwrap();
        for i in 0..42 {
            let j = g.anti[i] as usize;
            assert_ne!(i, j, "direction {i} is its own antipode");
            assert_eq!(g.anti[j] as usize, i, "anti is not an involution at {i}");
            for k in 0..3 {
                assert!(
                    (g.dirs[j][k] + g.dirs[i][k]).abs() < 1e-12,
                    "dirs[{j}] != -dirs[{i}]"
                );
            }
        }
    }

    /// Rotations commute with `x -> -x`, so the two ways of composing them
    /// must agree exactly. If this fails, one of the two tables is built
    /// with the wrong handedness and the binding kernel silently inherits it.
    #[test]
    fn antipode_commutes_with_every_rotation() {
        let g = Geodesic::<42>::build().unwrap();
        for perm in &g.perms {
            for i in 0..42 {
                assert_eq!(g.anti[perm[i] as usize], perm[g.anti[i] as usize]);
            }
        }
    }

    #[test]
    fn all_three_resolutions_build() {
        assert!(Geodesic::<12>::build().is_ok());
        assert!(Geodesic::<42>::build().is_ok());
        assert!(Geodesic::<162>::build().is_ok());
    }

    /// The group structure is a property of the icosahedron, not of the
    /// subdivision level, so it must hold at every resolution the spec allows.
    /// Checking only D=42 would let a level-dependent indexing bug through at
    /// the one resolution §22.2 might lock to instead.
    #[test]
    fn the_group_properties_hold_at_every_resolution() {
        fn check<const D: usize>() {
            let g = Geodesic::<D>::build().unwrap();
            let set: std::collections::BTreeSet<Vec<u8>> =
                g.perms.iter().map(|p| p.to_vec()).collect();
            assert_eq!(set.len(), N_ROTATIONS, "D={D}: rotations are not distinct");
            for a in &g.perms {
                for b in &g.perms {
                    let composed: Vec<u8> = (0..D).map(|i| b[a[i] as usize]).collect();
                    assert!(set.contains(&composed), "D={D}: not closed");
                }
            }
            assert!(
                !g.perms.iter().any(|p| p[..] == g.anti[..]),
                "D={D}: -I is in the rotation set"
            );
            for i in 0..D {
                assert_eq!(g.anti[g.anti[i] as usize] as usize, i, "D={D}: anti at {i}");
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// Ties the integer table back to geometry. Every other test here treats
    /// `perms` as an abstract permutation group — and a table of the *wrong*
    /// 60 rotations, or of some unrelated group of order 60, would satisfy all
    /// of them. This is the only check that says the lookup means what the
    /// binding kernel will assume it means.
    #[test]
    fn the_permutation_table_agrees_with_the_rotation_matrices() {
        let g = Geodesic::<42>::build().unwrap();
        let mats = rotation_matrices().unwrap();
        assert_eq!(mats.len(), N_ROTATIONS);
        for (r, m) in mats.iter().enumerate() {
            for i in 0..42 {
                let rotated = apply_mat(m, g.dirs[i]);
                let expected = g.dirs[g.perms[r][i] as usize];
                for k in 0..3 {
                    assert!(
                        (rotated[k] - expected[k]).abs() < 1e-12,
                        "rotation {r} on direction {i}: matrix and table disagree"
                    );
                }
            }
        }
    }

    /// The matrices must be rotations in their own right, checked without
    /// reference to the direction set: orthonormal rows and determinant `+1`.
    #[test]
    fn rotation_matrices_are_orthonormal_with_determinant_one() {
        for (r, m) in rotation_matrices().unwrap().iter().enumerate() {
            for i in 0..3 {
                for j in 0..3 {
                    let d = m[i][0] * m[j][0] + m[i][1] * m[j][1] + m[i][2] * m[j][2];
                    let want = if i == j { 1.0 } else { 0.0 };
                    assert!(
                        (d - want).abs() < 1e-12,
                        "rotation {r}: rows {i},{j} give {d}"
                    );
                }
            }
            let det = det3(m[0], m[1], m[2]);
            assert!(
                (det - 1.0).abs() < 1e-12,
                "rotation {r} has determinant {det}"
            );
        }
    }

    #[test]
    fn wrong_d_is_rejected_not_silently_wrong() {
        assert!(matches!(
            Geodesic::<43>::build(),
            Err(GeoError::UnsupportedResolution(43))
        ));
    }
}
