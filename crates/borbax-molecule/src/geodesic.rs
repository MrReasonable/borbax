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
//! Nothing here uses a transcendental, and the operation inventory is short
//! enough to state in full rather than gesture at: `+ - * /`, `sqrt`, `abs`,
//! comparison and `total_cmp` — all either exactly specified by IEEE-754 or
//! pure functions of the bits — **plus [`f64::midpoint`], which is neither**.
//! std specifies nothing about its rounding. On the pinned rustc it is a plain
//! `const fn` in `core` with no `cfg(target_arch)`, computing `(a + b) * 0.5`
//! when `|a|` and `|b|` are both `≤ f64::MAX / 2` and `(a * 0.5) + (b * 0.5)`
//! otherwise — every input here is a unit-vector component, so only the first
//! branch is reachable. Portable across targets, exposed only to a toolchain
//! bump. That
//! matters more than it looks: `intern`'s edge dedup relies on
//! `midpoint(u, v)` and `midpoint(v, u)` being bit-identical, which holds
//! because IEEE addition commutes but would not for a future `a + (b - a) * 0.5`.
//!
//! Measured rather than argued: the full bit dump — every direction, every
//! permutation entry, `anti`, and all 60 matrices — is identical across
//! aarch64 and x86-64, at two optimisation levels, and with
//! `-C target-cpu=x86-64-v3` enabling FMA3 and AVX2. There are zero FMA-class
//! instructions in the emitted assembly, because rustc sets no `contract` flag
//! and so LLVM cannot fuse `dot`'s three-term sums (§13.4).
//!
//! `.prototools` already makes a toolchain bump a deliberate act;
//! `the_canonical_direction_ordering_is_pinned` is what makes it a *visible*
//! one.
//!
//! Written in `borbax-experiments` first, because the G2 locality measurement
//! needed a rotation table before Task 7 was due, and lifted here when it was.
//! `borbax-experiments` re-exports this module rather than holding a second
//! copy, so the table G2 measured and the table the binding kernel searches
//! cannot drift apart — a divergence there would invalidate the measurement
//! without failing anything.
//!
//! **That is an argument for one *shared* table, not for one copy of the
//! construction.** A third icosahedron-and-subdivision lives in
//! `borbax_universe::packing`, deliberately: it is an independent oracle,
//! answerable to nothing but Euler, and importing the implementation under
//! test would destroy the only property that makes it ground truth. It also
//! records what careless duplication costs — without the canonical sort, that
//! copy moved the admitted `FRONTIER_COEFF` band by about 1.1 in opposite
//! directions at the two ends. The right number of copies here is two, and
//! the distinction that matters is *shared table* versus *independent
//! oracle*, never copy count.
//!
//! The lift was **not** verbatim, and the differences are the interesting part.
//! Four comments asserted "this lifts into `borbax-molecule`" and had to become
//! statements about a file that has arrived; the test module traded `.unwrap()`
//! for this crate's `unwrap_or_else(|| unreachable!(…))` idiom, because
//! `clippy.toml` sets no `allow-unwrap-in-tests`; and one measured claim —
//! "antipodes bit-exact" — was wrong, in the benign direction, which nothing
//! had ever checked. See [`Geodesic::anti`].

#![allow(
    clippy::indexing_slicing,
    reason = "fixed-size arrays indexed by loop bounds proven against D and N_ROTATIONS, or \
              by a `Rotation`, which cannot be constructed out of range; a get().ok_or() on \
              every coordinate access would obscure the geometry"
)]

/// Order of the icosahedral rotation group. Reflections are excluded — see
/// spec §22.8, which makes Borbax chemistry handed on purpose.
pub const N_ROTATIONS: usize = 60;

/// An index naming one of the [`N_ROTATIONS`] elements, and nothing else.
///
/// The permutation accessors take this rather than a `usize` because they index
/// a fixed-size table, and a `usize` parameter makes them partial functions that
/// panic on out-of-range input. CLAUDE.md treats a library panic as a
/// correctness bug — "a simulation that dies eight hours into an overnight run"
/// — and `clippy::panic` does not catch a panic that arrives via indexing.
///
/// It is not defensive programming against a caller who cannot exist. Task 10's
/// `affinity_with_rotation` returns its winning rotation for §14.5's renderer to
/// apply, and that round-trip — table to caller and back — is exactly where an
/// unvalidated index arrives. This makes the round-trip total.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rotation(usize);

impl core::fmt::Display for Rotation {
    /// Prints the bare index. Assertion messages name rotations constantly, and
    /// `Rotation(7)` in a failure would be noise.
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Rotation {
    /// The identity rotation, index 0.
    ///
    /// Named rather than spelled `Rotation::new(0)` at call sites: §8.3's
    /// kernel seeds its argmax with it, and a seed that reads as "index zero"
    /// invites someone to change it to a sentinel.
    pub const IDENTITY: Self = Self(0);

    /// The rotation at `r`, or `None` when `r >= N_ROTATIONS`.
    #[must_use]
    pub const fn new(r: usize) -> Option<Self> {
        if r < N_ROTATIONS { Some(Self(r)) } else { None }
    }

    /// The index, for storing, printing, or feeding back through [`Rotation::new`].
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }

    /// All 60, in table order — the spelling every sweep should use.
    pub fn all() -> impl Iterator<Item = Self> {
        (0..N_ROTATIONS).map(Self)
    }
}

/// Why a geodesic could not be built.
///
/// Only [`GeoError::UnsupportedResolution`] is reachable from caller input; the
/// other two are broken construction invariants. All three exist so the failure
/// is loud instead of silent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum GeoError {
    /// D must be one of the geodesic ladder values: 12, 42, 162.
    #[error("unsupported resolution {0}: D must be 12, 42 or 162")]
    UnsupportedResolution(usize),
    /// A rotation failed to map the vertex set onto itself. This would mean
    /// the symmetry assumption above is broken, which is a bug, not input.
    #[error("a rotation did not map the sample directions onto themselves")]
    RotationNotClosed,
    /// The ladder promised `expected` vertices at this level and subdivision
    /// produced `got`.
    ///
    /// **Not a bad `D`** — `D` was already on the ladder when this fires, which
    /// is why it is not [`GeoError::UnsupportedResolution`]. That variant's
    /// message reads "D must be 12, 42 or 162", and printing it here would deny
    /// that anything is wrong at the one moment someone reads it. This is the
    /// guard that stops the `copy_from_slice` in [`Geodesic::build`] panicking
    /// on a length mismatch, so it fires exactly when subdivision has broken
    /// deep down.
    #[error("subdivision at level {level} produced {got} vertices, not {expected}")]
    SubdivisionCountMismatch {
        /// Subdivision level, from `level_for`.
        level: u32,
        /// Vertices the ladder promises at `level`.
        expected: usize,
        /// Vertices subdivision actually produced.
        got: usize,
    },
}

/// The geodesic ladder: `(subdivision level, vertex count)`.
///
/// **One encoding, read by both [`vertex_count`] and `level_for`**, so the two
/// cannot disagree. They used to be independent — a `10 * 4usize.pow(level) + 2`
/// formula and a hand-written `match` — tied together only by a test checking
/// three points in one direction.
///
/// Collapsing them also deletes an overflow rather than arguing it unreachable.
/// The formula wrapped at `level = 31` (not 32, as a first draft of the plan
/// claimed), and at `level = 32` it returned **`2`** in release: not an
/// obviously broken number, from a `pub` function anyone could call, in a
/// codebase where a plausible wrong answer is the thing most worth fearing.
/// A lookup cannot overflow.
const LADDER: [(u32, usize); 3] = [(0, 12), (1, 42), (2, 162)];

/// Vertices produced by subdividing an icosahedron `level` times, or `None`
/// when `level` is off the ladder.
///
/// Fallible because the ladder is finite: spec §8.2 admits D ∈ {12, 42, 162}
/// and nothing else, since intermediate resolutions destroy the property that
/// makes a rotation an exact permutation.
#[must_use]
pub const fn vertex_count(level: u32) -> Option<usize> {
    let mut i = 0;
    while i < LADDER.len() {
        if LADDER[i].0 == level {
            return Some(LADDER[i].1);
        }
        i += 1;
    }
    None
}

/// Subdivision level yielding `d` vertices, or `None` when `d` is off the
/// ladder. The inverse of [`vertex_count`] over [`LADDER`], and
/// `the_ladder_is_a_bijection` holds them to it.
const fn level_for(d: usize) -> Option<u32> {
    let mut i = 0;
    while i < LADDER.len() {
        if LADDER[i].1 == d {
            return Some(LADDER[i].0);
        }
        i += 1;
    }
    None
}

/// Sample directions at resolution `D`, with the rotation and antipode tables
/// that make comparing two signatures a table lookup.
///
/// **Fields are private and [`Geodesic::build`] is the only producer.** Three
/// invariants are established there and are unenforceable afterwards:
/// `rotation_perms`
/// is closed under composition, `anti` is an involution with
/// `dirs[anti[i]] == -dirs[i]`, and every index in either table is `< D`. With
/// the fields public a caller could desynchronise them — `g.dirs.reverse()` is
/// enough — and nothing here would notice, because every test but
/// `the_canonical_direction_ordering_is_pinned` is invariant under relabelling.
/// That is the same hazard the `-0.0` probe found, reopened at runtime.
///
/// Keeping the storage shape private is also what made
/// [`Geodesic::contact_perms`] cheap to add rather than a breaking change.
// No `Clone`. Nothing in the workspace clones a `Geodesic`, it is 23.5 KB at
// D=162, and the module doc's "built once and passed by reference" is a
// contract a derived `Clone` quietly contradicts. Adding a trait impl later is
// additive; removing one is breaking, so the asymmetry says omit it now.
#[derive(Debug)]
pub struct Geodesic<const D: usize> {
    /// Unit sample directions, in a canonical order.
    dirs: [Vec3; D],
    /// `rotation_perms[r][i] = j` — rotation `r` carries direction `i` to `j`.
    rotation_perms: [[u8; D]; N_ROTATIONS],
    /// `contact_perms[r][i] = anti[rotation_perms[r][i]]`.
    contact_perms: [[u8; D]; N_ROTATIONS],
    /// `anti[i] = j` where `dirs[j] == -dirs[i]`.
    anti: [u8; D],
}

impl<const D: usize> Geodesic<D> {
    /// Unit sample directions, in the canonical order [`Geodesic::build`] fixed.
    ///
    /// The index is the coordinate system for every signature (§8.2), which is
    /// why the order is pinned by a golden rather than left as an output.
    #[must_use]
    pub const fn dirs(&self) -> &[Vec3; D] {
        &self.dirs
    }

    /// Where rotation `r` sends each direction: `rotation_perms(r)[i] = j`.
    ///
    /// **For rotating a *single* signature** — Task 9's `canonicalise` sweeps
    /// the 60 elements looking for a canonical orientation, and no contact is
    /// involved. **Not for binding.** See [`Geodesic::contact_perms`].
    #[must_use]
    pub const fn rotation_perms(&self, r: Rotation) -> &[u8; D] {
        &self.rotation_perms[r.0]
    }

    /// The direction that *touches* direction `i` when the partner is presented
    /// at orientation `R_rᵀ`: `contact_perms(r)[i]` is the index of
    /// `−R_r · dirs[i]`, equivalently `anti[rotation_perms(r)[i]]`.
    ///
    /// **The transpose is not a typo, and an earlier version of this line got
    /// it backwards.** Pairing `a[i]` with `b[contact_perms(r)[i]]` asks for the
    /// partner's extent along `−R_r·u_i`; that lands on `−u_i`, where contact
    /// happens, only once the partner is turned by `R_rᵀ`. Measured on a
    /// support function at D=42 over all 60×42 pairs:
    /// `b[contact_perms(r)[i]] == h_{R_rᵀ·B}(−u_i)` to 1.8e-15, while the
    /// forward reading is off by 2.26 — not a rounding question.
    ///
    /// It costs nothing in `affinity`, which maximises over all 60 and the
    /// group is closed under inverse. It costs a great deal in **§14.5's
    /// renderer**, which is handed the winning `r` precisely so it does not
    /// re-derive the pose, and must apply `R_rᵀ`. Applying `R_r` draws a
    /// plausible, silent, wrong docking picture — in a project that ships SVG
    /// from day one because a shape-based chemistry cannot be developed blind.
    ///
    /// **§8.3's binding kernel indexes its partner through this and nothing
    /// else.** Two bodies in contact touch along opposite directions, so the
    /// antipode is not a refinement of the rotation — it is half of what
    /// "in contact" means. Reaching for [`Geodesic::rotation_perms`] here
    /// silently converts the search into the 60 *improper* elements of the
    /// icosahedral group — reflections only, the opposite of §22.8 — and
    /// nothing fails.
    ///
    /// The composition is precomputed for two reasons and the second is the
    /// one that mattered. It removes a dependent load from the kernel's
    /// innermost statement, worth a measured **−13.2%** at D = 42 against
    /// `anti[rotation_perms[r][i]]`, bit-identical. And it gives the correct
    /// operation the shorter name: a reviewer reading `rotation_perms` inside
    /// a function called `affinity` sees something wrong, where `perms` — the
    /// name this field used to carry — read like the obvious choice. Adding a
    /// table does not remove an affordance; renaming the other one does.
    #[must_use]
    pub const fn contact_perms(&self, r: Rotation) -> &[u8; D] {
        &self.contact_perms[r.0]
    }

    /// `anti()[i] = j` where `dirs()[j] == -dirs()[i]`.
    ///
    /// Exposed because Task 9 and the harnesses want the antipode directly.
    /// The binding kernel should not: [`Geodesic::contact_perms`] already
    /// carries the composition.
    ///
    /// The vertex set is antipodally closed at every subdivision level, so this
    /// is total; `build_anti` errors rather than approximating if it ever is not.
    #[must_use]
    pub const fn anti(&self) -> &[u8; D] {
        &self.anti
    }
}

/// A direction or point in three dimensions.
///
/// Public because [`apply_mat`] is: a `pub fn` signed with a private alias
/// renders in rustdoc as a bare `[f64; 3]` beside a linkable `Mat3`, so a
/// reader sees two spellings for one idea and cannot tell which means
/// "direction".
pub type Vec3 = [f64; 3];

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
/// pair by the same norm, so `dirs[anti[i]]` is the *exact* negation of
/// `dirs[i]` rather than merely close to it — `dirs[anti[i]][k] + dirs[i][k]`
/// is `0.0` in every component at every level, measured, not argued.
///
/// Exact numerically, not bitwise, and the distinction is worth keeping
/// because the two claims point at different hazards. The zero component of
/// each `(0, ±1, ±φ)` vertex is stored as `+0.0` and negates to `-0.0`: equal
/// under `==`, different under `to_bits`. Numeric equality is what `dist2` and
/// [`Geodesic::build_anti`]'s tolerance see. Bitwise, the fact that matters is
/// the *other* one — that no `-0.0` is ever stored — because that is what makes
/// the `total_cmp` sort in [`Geodesic::build`] safe.
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
    // unreachable here — they cannot be negative, so the two zeros cannot both
    // arise — but the comparison costs nothing and does not depend on that
    // argument staying true.
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
    /// [`GeoError::SubdivisionCountMismatch`] if subdivision at the ladder
    /// level for `D` yields a vertex count other than `D`.
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
            return Err(GeoError::SubdivisionCountMismatch {
                level,
                expected: D,
                got: verts.len(),
            });
        }

        // Canonical ordering, and the *reason* matters more than the call.
        //
        // `total_cmp` is a pure function of the bits, so it is portable when
        // the bits are — not, as an earlier version of this comment claimed,
        // free of platform variation in general. It orders on the sign bit,
        // and a runtime NaN's sign bit differs between aarch64 and x86-64;
        // that is precisely the defect `borbax_units::Span::canonical_cmp`
        // exists to fix, and a comment asserting the opposite is how the wrong
        // lesson gets cited a year from now.
        //
        // It is safe *here* because these coordinates are free of both hazards,
        // measured at all three levels: 0 `-0.0` components, 0 NaN, max
        // |‖v‖-1| = 1.11e-16 (one ulp below 1.0, at D=162; exact at 12 and
        // 42). Every
        // value is a normalised coordinate built from sums of squares.
        //
        // The predecessor of this comment ended "when this lifts, that
        // reasoning lifts with it or the call changes". It has lifted, so the
        // reasoning was re-measured rather than re-copied — which is how the
        // **third** item in that list, "antipodes bit-exact", turned out to be
        // false at every zero component and was corrected (see [`icosahedron`]).
        // The fourth, the norm bound, is retained above because it is true.
        // The two claims this exemption actually rests on are the first two,
        // and they are no longer only a comment:
        // `the_direction_set_carries_no_negative_zero_and_no_nan` fails if
        // either stops holding, at all three levels (spec §13.4).
        //
        // **What those two claims do NOT determine is the resulting order**,
        // and that is the sharper hazard. The sort key genuinely ties: 121 of
        // the 161 adjacent pairs at D=162 share a bit-identical `z` (29/41 at
        // 42, 7/11 at 12), because the icosahedron in the `(0, ±1, ±φ)`
        // orientation has 2-fold axes along all three coordinate axes and puts
        // whole families of directions at one latitude. They resolve today only
        // because mathematically equal `z` values round to the same bits.
        // Replacing `scale(a, 1.0 / n)` in `normalise` with componentwise
        // division — mathematically identical, *better* rounded, no behavioural
        // intent — splits two of the four directions at z = -0.85065… by one
        // ulp and **moves six of the 162 index positions** at D=162 (none at
        // 12 or 42): `dirs[12]` moves to a direction 0.743 away, and every
        // other test still passes in both profiles.
        // `the_canonical_direction_ordering_is_pinned` is what catches that.
        #[expect(
            clippy::disallowed_methods,
            reason = "§13.4: measured free of NaN and -0.0 at all three levels — see above"
        )]
        verts.sort_by(|a, b| {
            a[2].total_cmp(&b[2])
                .then(a[1].total_cmp(&b[1]))
                .then(a[0].total_cmp(&b[0]))
        });

        // **Nothing may transform between the sort and this copy.** Sorting is
        // a permutation and `copy_from_slice` is bitwise, so the multiset of
        // bit patterns in `dirs` is exactly the multiset in `verts` at sort
        // time — which is the only reason a test inspecting `dirs` can license
        // an exemption about the sort's *input*. Insert a re-normalisation pass
        // or a `dirs[i] = f(verts[i])` here and
        // `the_direction_set_carries_no_negative_zero_and_no_nan` silently
        // stops covering what it claims to cover.
        let mut dirs = [[0.0; 3]; D];
        dirs.copy_from_slice(&verts);

        let rotation_perms = Self::build_perms(&dirs)?;
        let anti = Self::build_anti(&dirs)?;

        // `contact_perms = anti ∘ rotation_perms`, composed once here so the
        // binding kernel does one load instead of two. Bit-identical to
        // composing at the call site by construction — same values, same
        // order, only the addressing differs — and
        // `contact_perms_is_anti_composed_with_rotation_perms` holds it to that.
        let mut contact_perms = [[0u8; D]; N_ROTATIONS];
        for (r, row) in contact_perms.iter_mut().enumerate() {
            for (i, cell) in row.iter_mut().enumerate() {
                *cell = anti[usize::from(rotation_perms[r][i])];
            }
        }

        Ok(Self {
            dirs,
            rotation_perms,
            contact_perms,
            anti,
        })
    }

    /// Index of `-dirs[i]` for every `i`.
    ///
    /// The tolerance is `1e-24` on a squared distance — effectively exact, and
    /// it can be, because the vertex set is built so antipodal pairs are exact
    /// negations (see [`icosahedron`]). The measured `dist2` is `0.0` at every
    /// index at all three levels.
    ///
    /// "Headroom over zero" would be meaningless, so state what the tolerance
    /// actually brackets: it must accept `0.0` and reject the nearest wrong
    /// candidate, which sits at `dist2 = 7.61e-2`. That is 22 orders above
    /// `1e-24`. What the tightness buys is rejecting a *drifted* construction —
    /// a one-ulp perturbation would land near `dist2 ≈ 1.2e-32` and still pass,
    /// so this is not infinitely strict, but anything that has genuinely lost
    /// antipodal closure is nowhere near. A nearest match instead of an error
    /// would corrupt every binding comparison rather than failing visibly.
    ///
    /// "Exact" is numeric, not bitwise: at a zero component `dirs` holds `+0.0`
    /// and the negation is `-0.0`. `dist2` subtracts them to `+0.0`, so the
    /// tolerance sees exactness; `to_bits` would not. An earlier version of the
    /// comment above said "bit-exact", which is wrong in exactly this way — 12,
    /// 24 and 48 components at the three levels — and nothing tested it either
    /// way. `the_antipode_table_is_an_involution_on_opposite_directions` now
    /// asserts the numeric equality that this tolerance actually depends on.
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
                // accumulates three products per coordinate.
                //
                // `1e-12` is a **squared** distance, so comparing it to the
                // ~0.28 linear spacing between neighbouring directions at
                // D = 162 — as an earlier version of this comment did — is a
                // dimensional mismatch. The two quantities that actually
                // bracket it, both measured: the worst true hit is
                // `dist2 = 5.55e-31`, and the nearest wrong direction is at
                // `dist2 = 7.61e-2`. So the tolerance sits 18 orders above
                // what it must accept and 11 below what it must reject.
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

/// Whether a group element is the identity, to within the error of building it.
///
/// Lives here rather than beside its caller because it is a fact about
/// [`rotation_matrices`]. Its only callers are the G2 harness and the
/// signature probes in `borbax-experiments`, one crate up; those call sites do
/// not come with the module, which is why
/// `is_identity_selects_exactly_one_of_the_sixty` was written during the lift
/// rather than leaving the predicate to ship here untested.
///
/// **The binding kernel does not need it.** Task 10's plan section consumes
/// `Geodesic` and its permutation tables and asks for neither this nor
/// [`rotation_matrices`]; the placement here rests on the harness use above,
/// which is real and current.
///
/// Not `m == IDENTITY`: the 60 elements are frame products, so the identity
/// arrives with entries at `1.0000000000000004` and `-5.6e-17`, and an exact
/// comparison finds it in **zero** of the 60.
///
/// **The tolerance is structural, and the basis-free form of the argument is
/// the one to trust.** For a rotation by angle θ, `‖R − I‖_F = 2√2·|sin(θ/2)|`,
/// and the largest entry of a 3×3 matrix is at least `‖·‖_F / 3`. The smallest
/// non-identity angle in the icosahedral group is 72°, so **no element can come
/// within `2√2·sin36°/3 = 0.554` of `I` entrywise, in any orientation of the
/// icosahedron.** That bound needs no re-measuring when anything nearby
/// changes, and 1e-9 sits 8.74 orders below it.
///
/// In this construction the figure is tighter and exact: the nearest other
/// element is 0.80901699437494742…, which is φ/2 = cos 36°, attained by all
/// twelve 72° rotations. **That number is a constant of the group only in the
/// `(0, ±1, ±φ)` frame.** Max entrywise deviation is not a similarity
/// invariant — conjugating the group by a 0.3 rad rotation about `z` gives
/// 0.717 instead — so an earlier version of this paragraph over-claimed by
/// calling it independent of everything. It *is* independent of `D`, of the
/// element radii, and of which `(vertex, neighbour)` pair seeds `f0`: since
/// `frame(gv, gw) = frame(v, w)·gᵀ`, the set of 60 matrices is the group
/// whichever pair is chosen. Reorienting the base icosahedron is the one
/// change that moves it, and the 0.554 bound above is what keeps 1e-9 safe if
/// it does.
///
/// Contrast the G2 harness's `ROTATION_MUST_MOVE_RAW`, whose margin *is*
/// proportional to molecule size and radius scale, and which documents that
/// exposure. Here the identity deviates from `I` by 4.44e-16 entrywise —
/// 2 ulp **above** 1.0 — so 1e-9 sits 6.35 orders above the construction noise
/// and 8.91 below the nearest real rotation, selecting exactly one matrix.
///
/// Every figure in the preceding paragraph is asserted by
/// `is_identity_selects_exactly_one_of_the_sixty` — both edges of the window
/// and the count. "Structural" is a claim about why the numbers cannot drift,
/// not a licence to leave them unchecked.
///
/// The predicate also never has to discriminate a *near*-identity. Every
/// caller runs [`Geodesic::build`] first, and `build_perms` rejects any matrix
/// that fails to carry the direction set onto itself, so the only elements
/// that reach here are exact icosahedral symmetries — the identity, or at
/// least 72°. The band where this predicate says "not the identity" while the
/// rotation is too small to move a molecule is θ ∈ (1e-9, 1.26e-2) rad, and
/// nothing can present an element in it.
///
/// Entrywise rather than by trace or by angle, deliberately. `det_math`'s own
/// documentation records that `(tr(R) − 1) / 2` exceeds 1 for 685 of the 3,600
/// compositions — the identity among them — so a trace predicate needs a
/// clamp; and `det_math::acos` loses angles below ~1.5e-8 rad entirely, which
/// is *above* this tolerance, so an angle predicate would compare inside its
/// own noise floor. Entrywise is also the only one of the three that does not
/// assume its input is orthogonal, which matters when part of the job is to
/// notice a construction that has broken.
#[must_use]
pub fn is_identity(m: &Mat3) -> bool {
    /// Entrywise tolerance. See [`is_identity`] for why it is structural.
    const TOL: f64 = 1e-9;
    const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    m.iter()
        .flatten()
        .zip(IDENTITY.iter().flatten())
        .all(|(a, b)| (a - b).abs() < TOL)
}

/// The 60 rotations as matrices, in the same order as
/// [`Geodesic::rotation_perms`].
///
/// [`Geodesic::contact_perms`] is what the binding kernel needs — a rotation
/// composed with the antipode, as a table lookup.
/// This is the same 60 rotations as actual geometry, and it exists so the two
/// can be checked against each other. Without it, the table is one of
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
    clippy::indexing_slicing,
    clippy::as_conversions,
    reason = "indices are loop bounds over fixed-size arrays; every cast is u8 -> usize, a \
              widening that cannot lose information — `as_conversions` is here to catch \
              truncation, and `usize::from` on each index would bury the assertion"
)]
mod tests {
    use super::*;

    /// Build a geodesic, or fail naming the resolution that broke.
    ///
    /// Not `.unwrap()`. `clippy.toml` sets neither `allow-unwrap-in-tests` nor
    /// `allow-expect-in-tests`, so the workspace `unwrap_used = "deny"` reaches
    /// inside `#[cfg(test)]`, and this crate carries zero bare unwraps. The
    /// `unreachable!` spelling is preferred anyway because it names its own
    /// precondition: `build` is infallible at the three resolutions the ladder
    /// admits, so a failure here means that stopped being true — which is what
    /// the message says, rather than `called Result::unwrap on an Err value`.
    fn geo<const D: usize>() -> Geodesic<D> {
        Geodesic::<D>::build().unwrap_or_else(|e| unreachable!("D={D} is on the ladder: {e}"))
    }

    /// The 60 rotations as matrices, or fail naming why.
    fn mats() -> Vec<Mat3> {
        rotation_matrices().unwrap_or_else(|e| unreachable!("the icosahedron has 60 edges: {e}"))
    }

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
        assert_eq!(vertex_count(0), Some(12));
        assert_eq!(vertex_count(1), Some(42));
        assert_eq!(vertex_count(2), Some(162));
        // Off the ladder returns `None` rather than a wrapped formula. The
        // predecessor computed `10 * 4usize.pow(level) + 2` and returned `2`
        // at level 32 in release.
        assert_eq!(vertex_count(3), None);
        assert_eq!(vertex_count(31), None);
        assert_eq!(vertex_count(32), None);
    }

    /// [`vertex_count`] and `level_for` read one table, so they are inverse by
    /// construction — this is what makes that checkable rather than asserted.
    /// It fails if anyone reintroduces a second encoding of the ladder.
    #[test]
    fn the_ladder_is_a_bijection() {
        for (level, count) in LADDER {
            assert_eq!(vertex_count(level), Some(count));
            assert_eq!(level_for(count), Some(level));
        }
        assert_eq!(level_for(43), None);
    }

    /// **The only check here that is not invariant under relabelling `dirs`.**
    ///
    /// Every other test in this module asserts a *property* — group closure,
    /// bijectivity, distinctness, det = +1, the `anti` involution, `-I` absent,
    /// matrix/table agreement, unit norms — and every one of them survives an
    /// arbitrary permutation of the direction set with the matching permutation
    /// of `rotation_perms` and `anti`, because that is still a perfectly
    /// self-consistent
    /// table. Two review lanes measured this independently: re-keying the sort
    /// `(z,y,x)` → `(x,y,z)` leaves the whole 325-test workspace green.
    ///
    /// The index is not an implementation detail. It is the coordinate system
    /// for every signature (§8.2), so a silent permutation changes every
    /// downstream result while nothing here fails.
    ///
    /// The realistic mutation is not a deliberate re-key, it is a tidy-up.
    /// Rewriting `normalise`'s `scale(a, 1.0 / n)` as componentwise division is
    /// mathematically identical, *better* rounded, unflagged by clippy, and is
    /// how `borbax_universe::packing`'s independent copy already spells it —
    /// and it moves six of the 162 index positions at D=162, because a one-ulp
    /// difference in the `z` key splits directions sitting at one latitude.
    ///
    /// A move here is a physics change: regenerate deliberately (§13.6), never
    /// because it moved. Because the three `test` legs run on three platforms,
    /// this is also a working §13.4 golden for this table today, rather than
    /// waiting for Task 20's `goldens --emit`.
    #[test]
    fn the_canonical_direction_ordering_is_pinned() {
        /// FNV-1a over `u64` lanes. Hand-rolled rather than `DefaultHasher`,
        /// whose output std documents as unstable across releases — the whole
        /// value of these constants is that they cannot move for a reason
        /// outside this file.
        fn fnv1a(vals: impl IntoIterator<Item = u64>) -> u64 {
            let mut h = 0xcbf2_9ce4_8422_2325_u64;
            for v in vals {
                h = (h ^ v).wrapping_mul(0x0000_0100_0000_01b3);
            }
            h
        }
        // Raw `to_bits()`, not `canonical_bits()`. This is a *detector*, and
        // canonicalisation is lossy by design: it would fold a `-0.0` that
        // differed by architecture onto `+0.0` and report "unchanged". The
        // sibling test proves there is no `-0.0` or NaN here to fold.
        fn check<const D: usize>(want_dirs: u64, want_perms: u64, want_anti: u64) {
            let g = geo::<D>();
            let dirs = fnv1a(g.dirs().iter().flatten().map(|c| c.to_bits()));
            let perms = fnv1a(
                Rotation::all().flat_map(|r| g.rotation_perms(r).iter().map(|&p| u64::from(p))),
            );
            let anti = fnv1a(g.anti().iter().map(|&a| u64::from(a)));
            assert_eq!(dirs, want_dirs, "D={D}: dirs moved");
            assert_eq!(perms, want_perms, "D={D}: rotation_perms moved");
            assert_eq!(anti, want_anti, "D={D}: anti moved");
            // `contact_perms` is derived, so it gets no constant of its own —
            // a constant here would be a second encoding of the same fact and
            // could drift from the composition it is supposed to be.
            // `contact_perms_is_anti_composed_with_rotation_perms` pins it
            // against the two tables above, which is stronger.
        }
        check::<12>(
            0x71d3_a418_1bf1_6135,
            0x615d_533d_d6ca_ad15,
            0xa033_525a_4dc8_63f9,
        );
        check::<42>(
            0x25b9_6c8b_ef07_03f5,
            0x729e_90ea_76d5_23dd,
            0xe680_99a5_c14e_0434,
        );
        check::<162>(
            0xa106_c746_df87_c901,
            0x23c8_7115_e46e_1edd,
            0xc5ef_782c_3a5a_3684,
        );

        // The 60 matrices, which are `D`-independent and so sit outside `check`.
        // Pinned because the module doc's cross-ISA bit-dump claim names them
        // and nothing else did: `the_permutation_table_agrees_with_the_rotation_matrices`
        // admits 1e-12, so a change to `frame`'s arithmetic could move these
        // bits with every other test in the file green.
        assert_eq!(
            fnv1a(mats().iter().flatten().flatten().map(|c| c.to_bits())),
            0xe318_1547_c216_9bfd,
            "rotation_matrices moved"
        );
    }

    /// Checked at every resolution rather than only at D=42: normalisation
    /// error is the one quantity here that grows with subdivision depth, so
    /// D=42 is the level least able to detect it.
    ///
    /// The tolerance is set from measurement, not from habit. The worst
    /// deviation across the whole ladder is 1.11e-16 — exactly one ulp *below*
    /// 1.0 (the ulp above is 2.22e-16), at D=162; D=12 and D=42 are exact. Every coordinate is built from sums,
    /// products and `sqrt`, all of which IEEE-754 specifies exactly, so that
    /// figure is the same on every target (§13.4) and 1e-15 sits ~9x above it
    /// rather than being a round number chosen for comfort.
    #[test]
    fn directions_are_unit_vectors() {
        fn check<const D: usize>() {
            for d in geo::<D>().dirs() {
                assert!((norm(*d) - 1.0).abs() < 1e-15, "D={D}: not unit: {d:?}");
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// The premise the `total_cmp` exemption in [`Geodesic::build`] rests on.
    ///
    /// `total_cmp` orders on the sign bit, so it separates `-0.0` from `+0.0`
    /// and it ranks NaN by a sign that differs between aarch64 and x86-64.
    /// Neither hazard can reach the sort here — but that is a *measured*
    /// property of the coordinates, and until this test existed it was only a
    /// sentence in a comment. It is the sentence that licenses the
    /// `#[expect(clippy::disallowed_methods)]`, so it needs something that
    /// fails when it stops being true.
    ///
    /// Zero components are real and expected: the icosahedron's vertices are
    /// the cyclic permutations of `(0, ±1, ±φ)`, giving 12, 24 and 48 exact
    /// zeros at the three levels. What matters is that every one is `+0.0`.
    #[test]
    fn the_direction_set_carries_no_negative_zero_and_no_nan() {
        fn check<const D: usize>() {
            for (i, v) in geo::<D>().dirs().iter().enumerate() {
                for (k, &c) in v.iter().enumerate() {
                    assert!(!c.is_nan(), "D={D}: NaN at dirs[{i}][{k}]");
                    assert!(
                        !(c == 0.0 && c.is_sign_negative()),
                        "D={D}: -0.0 at dirs[{i}][{k}] — the total_cmp sort is no longer safe"
                    );
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// Weak on its own — two builds in one process share every code path, so
    /// this cannot see a cross-platform or cross-toolchain difference. It is
    /// kept because it localises a nondeterminism *within* `build` (an
    /// address-dependent order, say) to this function rather than to whatever
    /// consumes the table. `the_canonical_direction_ordering_is_pinned` is the
    /// test that carries the §13.4 weight.
    #[test]
    fn construction_is_deterministic() {
        let a = geo::<42>();
        let b = geo::<42>();
        assert_eq!(a.dirs(), b.dirs());
        assert_eq!(a.anti(), b.anti());
        for r in Rotation::all() {
            assert_eq!(a.rotation_perms(r), b.rotation_perms(r));
            assert_eq!(a.contact_perms(r), b.contact_perms(r));
        }
    }

    /// `contact_perms` must be exactly `anti ∘ rotation_perms`, at every
    /// rotation and every direction and every resolution.
    ///
    /// This is the whole safety argument for precomposing. The kernel gets a
    /// one-load table because this holds; if it ever stopped holding, binding
    /// would search something that is neither the rotations nor the
    /// reflections, and every other test here would still pass — `contact_perms`
    /// participates in none of them.
    #[test]
    fn contact_perms_is_anti_composed_with_rotation_perms() {
        fn check<const D: usize>() {
            let g = geo::<D>();
            for r in Rotation::all() {
                for i in 0..D {
                    assert_eq!(
                        g.contact_perms(r)[i],
                        g.anti()[g.rotation_perms(r)[i] as usize],
                        "D={D}: contact_perms disagrees with anti∘rotation_perms at ({r},{i})"
                    );
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// Ties `contact_perms` to geometry **directly**, not transitively.
    ///
    /// `contact_perms_is_anti_composed_with_rotation_perms` pins it against the
    /// other two tables, which is the right invariant but says nothing about
    /// what the composition *means*. This says it: `contact_perms(r)[i]` is the
    /// index of `−R_r · dirs[i]`, measured against the matrices themselves.
    ///
    /// It exists because a cross-check found the docstring naming the wrong
    /// pose while the formula beside it was correct — the sort of error no
    /// table-versus-table assertion can see, and one §14.5's renderer would
    /// have inherited as a silently wrong docking picture.
    #[test]
    fn contact_perms_is_the_antipode_of_the_rotated_direction() {
        fn check<const D: usize>() {
            let g = geo::<D>();
            for (r, m) in Rotation::all().zip(mats().iter()) {
                for i in 0..D {
                    let want = apply_mat(m, g.dirs()[i]);
                    let got = g.dirs()[g.contact_perms(r)[i] as usize];
                    for k in 0..3 {
                        assert!(
                            (got[k] + want[k]).abs() < 1e-12,
                            "D={D}: contact_perms[{r}][{i}] is not the index of -R_{r}.dirs[{i}]"
                        );
                    }
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// The precomposed table is a *coset*, not the rotation group, and that is
    /// the point: `contact_perms` is `{-I · R}`, the 60 improper elements. If
    /// it ever coincided with `rotation_perms` the antipode step would have
    /// vanished — which is exactly the defect CLAUDE.md's ANTI bullet names,
    /// and the reason the two tables must never be interchangeable.
    #[test]
    fn contact_perms_is_never_the_rotation_table() {
        let g = geo::<42>();
        for r in Rotation::all() {
            for q in Rotation::all() {
                assert_ne!(
                    g.contact_perms(r)[..],
                    g.rotation_perms(q)[..],
                    "contact_perms[{r}] equals rotation_perms[{q}]: the antipode step is gone"
                );
            }
        }
    }

    /// Parameterised over `D` rather than pinned at 42, because the `[false; D]`
    /// scratch array is sized by the resolution: hardcoding `42` here while the
    /// geodesic is built at another rung is an out-of-bounds index, not a
    /// compile error. §22.2 locks `D` to one of three values and this test
    /// should not have an opinion about which.
    #[test]
    fn every_rotation_is_a_bijection() {
        fn check<const D: usize>() {
            let g = geo::<D>();
            for r in Rotation::all() {
                let mut seen = [false; D];
                for &j in g.rotation_perms(r) {
                    assert!(!seen[j as usize], "D={D}: rotation {r} is not injective");
                    seen[j as usize] = true;
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    #[test]
    fn rotations_are_all_distinct() {
        let g = geo::<42>();
        for i in Rotation::all() {
            for j in Rotation::all().filter(|j| j.index() > i.index()) {
                assert_ne!(
                    g.rotation_perms(i),
                    g.rotation_perms(j),
                    "rotations {i} and {j} coincide"
                );
            }
        }
    }

    /// The strongest correctness check available: the permutations must be
    /// closed under composition. If they are, they genuinely form the
    /// icosahedral rotation group, and §8.2's claim that rotation is an exact
    /// permutation of the sample directions holds.
    /// Parameterised for the same reason as the bijection test: a `(0..42)`
    /// composition against a geodesic built at another rung truncates silently
    /// rather than failing.
    #[test]
    fn permutations_form_a_group() {
        fn check<const D: usize>() {
            let g = geo::<D>();
            let set: std::collections::BTreeSet<Vec<u8>> = Rotation::all()
                .map(|r| g.rotation_perms(r).to_vec())
                .collect();
            assert_eq!(set.len(), N_ROTATIONS, "D={D}");
            for a in Rotation::all().map(|r| g.rotation_perms(r)) {
                for b in Rotation::all().map(|q| g.rotation_perms(q)) {
                    let composed: Vec<u8> = (0..D).map(|i| b[a[i] as usize]).collect();
                    assert!(
                        set.contains(&composed),
                        "D={D}: not closed under composition"
                    );
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    #[test]
    fn identity_is_present() {
        let g = geo::<42>();
        let identity: Vec<u8> = (0..42u8).collect();
        assert!(Rotation::all().any(|r| g.rotation_perms(r).to_vec() == identity));
        // **And it is at index 0, which `Rotation::IDENTITY` asserts by
        // construction and nothing checked.** §8.3's kernel seeds its argmax
        // with `IDENTITY`, and `Fit::pose` hands that index to §14.5's
        // renderer — so if index 0 were some other rotation, a pair whose best
        // score never beat the seed would be drawn in a pose nobody chose.
        assert_eq!(
            g.rotation_perms(Rotation::IDENTITY).to_vec(),
            identity,
            "Rotation::IDENTITY is not the identity permutation"
        );
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
        let g = geo::<42>();
        let (a, b, c) = best_triple(g.dirs());
        let source = det3(g.dirs()[a], g.dirs()[b], g.dirs()[c]);
        for (r, p) in Rotation::all()
            .enumerate()
            .map(|(n, r)| (n, g.rotation_perms(r)))
        {
            let image = det3(
                g.dirs()[p[a] as usize],
                g.dirs()[p[b] as usize],
                g.dirs()[p[c] as usize],
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
        let g = geo::<42>();
        assert!(
            !Rotation::all().any(|r| g.rotation_perms(r)[..] == g.anti()[..]),
            "-I is in the rotation set: the search is over reflections, not rotations"
        );
    }

    /// The antipode table is what makes `affinity` search rotations rather
    /// than reflections (§8.3, §22.8). These three properties are what it
    /// means for it to be correct.
    ///
    /// The negation is asserted as **exact equality**, not within a tolerance,
    /// and that is a statement about the construction rather than an act of
    /// optimism: [`icosahedron`] generates antipodal pairs as exact IEEE
    /// negations and `normalise` divides both members by the same norm, so
    /// `dirs[anti[i]] + dirs[i]` is `0.0` in every component at every level —
    /// measured, not hoped. A `< 1e-12` version of this assertion would pass
    /// for a construction that had merely drifted close, which is precisely
    /// what [`Geodesic::build_anti`]'s 1e-24 tolerance is there to forbid.
    ///
    /// Note it is *numeric* equality, not bit equality. `dirs` holds `+0.0` at
    /// the zero components (12, 24 and 48 of them at the three levels) while
    /// their negation is `-0.0`; the two compare equal and differ in bits. It
    /// is the numeric equality that `dist2` sees and that the tolerance means.
    #[test]
    #[allow(
        clippy::float_cmp,
        reason = "exact equality is the property under test — see the doc above"
    )]
    fn the_antipode_table_is_an_involution_on_opposite_directions() {
        fn check<const D: usize>() {
            let g = geo::<D>();
            for i in 0..D {
                let j = g.anti()[i] as usize;
                assert_ne!(i, j, "D={D}: direction {i} is its own antipode");
                assert_eq!(
                    g.anti()[j] as usize,
                    i,
                    "D={D}: anti is not an involution at {i}"
                );
                for k in 0..3 {
                    assert_eq!(
                        g.dirs()[j][k],
                        -g.dirs()[i][k],
                        "D={D}: dirs[{j}][{k}] is not the exact negation of dirs[{i}][{k}]"
                    );
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// Rotations commute with `x -> -x`, so the two ways of composing them
    /// must agree exactly. If this fails, one of the two tables is built
    /// with the wrong handedness and the binding kernel silently inherits it.
    #[test]
    fn antipode_commutes_with_every_rotation() {
        let g = geo::<42>();
        for r in Rotation::all() {
            let perm = g.rotation_perms(r);
            for i in 0..42 {
                assert_eq!(g.anti()[perm[i] as usize], perm[g.anti()[i] as usize]);
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
            let g = geo::<D>();
            let set: std::collections::BTreeSet<Vec<u8>> = Rotation::all()
                .map(|r| g.rotation_perms(r).to_vec())
                .collect();
            assert_eq!(set.len(), N_ROTATIONS, "D={D}: rotations are not distinct");
            for a in Rotation::all().map(|r| g.rotation_perms(r)) {
                for b in Rotation::all().map(|q| g.rotation_perms(q)) {
                    let composed: Vec<u8> = (0..D).map(|i| b[a[i] as usize]).collect();
                    assert!(set.contains(&composed), "D={D}: not closed");
                }
            }
            assert!(
                !Rotation::all().any(|r| g.rotation_perms(r)[..] == g.anti()[..]),
                "D={D}: -I is in the rotation set"
            );
            for i in 0..D {
                assert_eq!(
                    g.anti()[g.anti()[i] as usize] as usize,
                    i,
                    "D={D}: anti at {i}"
                );
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// Ties the integer table back to geometry. Every other test here treats
    /// `rotation_perms` as an abstract permutation group — and a table of the *wrong*
    /// 60 rotations, or of some unrelated group of order 60, would satisfy all
    /// of them. This is the only check that says the lookup means what the
    /// binding kernel will assume it means.
    ///
    /// Run at every resolution, for the same reason
    /// `the_group_properties_hold_at_every_resolution` is: this is the test
    /// carrying the most weight, and restricting the one test that touches
    /// geometry to a single `D` is how a level-dependent indexing error
    /// survives. The inverse mistake is already recorded in this project —
    /// a coefficient validated at D=162, the one level its table never used.
    #[test]
    fn the_permutation_table_agrees_with_the_rotation_matrices() {
        fn check<const D: usize>() {
            let g = geo::<D>();
            let mats = mats();
            assert_eq!(mats.len(), N_ROTATIONS);
            for (r, m) in Rotation::all().zip(mats.iter()) {
                for i in 0..D {
                    let rotated = apply_mat(m, g.dirs()[i]);
                    let expected = g.dirs()[g.rotation_perms(r)[i] as usize];
                    for k in 0..3 {
                        assert!(
                            (rotated[k] - expected[k]).abs() < 1e-12,
                            "D={D}: rotation {r} on direction {i}: matrix and table disagree"
                        );
                    }
                }
            }
        }
        check::<12>();
        check::<42>();
        check::<162>();
    }

    /// The matrices must be rotations in their own right, checked without
    /// reference to the direction set: orthonormal rows and determinant `+1`.
    #[test]
    fn rotation_matrices_are_orthonormal_with_determinant_one() {
        for (r, m) in mats().iter().enumerate() {
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

    /// [`is_identity`] documents a *structural* tolerance — one that needs no
    /// re-measuring when something nearby changes — and this is what makes
    /// that claim checkable rather than a paragraph.
    ///
    /// The claim has three parts and all three are asserted: exactly one of
    /// the 60 is the identity; it deviates from `I` by ~1 ulp; and the nearest
    /// other element is 0.809… away, which is φ/2 = cos 36°, a constant of the
    /// icosahedral group independent of `D`, of radius and of any fixture.
    /// `TOL = 1e-9` therefore sits in a window over fifteen orders wide, and
    /// this test fails if either edge of that window moves.
    ///
    /// Without it, lifting the module out of `borbax-experiments` would have
    /// left `is_identity` shipping with no test in its own crate: its only
    /// coverage was a call site in the G2 harness, which does not come along.
    #[test]
    fn is_identity_selects_exactly_one_of_the_sixty() {
        const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mats = mats();
        assert_eq!(mats.iter().filter(|m| is_identity(m)).count(), 1);

        // An explicit comparison rather than `f64::max`, which §13.1 bans for
        // the reason `min_pair_dist2` documents: for inputs that compare equal
        // — `+0.0` and `-0.0` — std may return either, and it was measured
        // returning different signs by architecture. Unreachable for these
        // absolute values, and still not worth spelling the banned way.
        let deviation = |m: &Mat3| {
            let mut worst = 0.0_f64;
            for i in 0..3 {
                for j in 0..3 {
                    let d = (m[i][j] - IDENTITY[i][j]).abs();
                    if d > worst {
                        worst = d;
                    }
                }
            }
            worst
        };
        let mut identity_dev = f64::MAX;
        let mut nearest_other = f64::MAX;
        for m in &mats {
            let d = deviation(m);
            if is_identity(m) {
                identity_dev = d;
            } else if d < nearest_other {
                nearest_other = d;
            }
        }
        assert!(
            identity_dev < 1e-15,
            "the identity now deviates by {identity_dev}, not ~1 ulp"
        );
        assert!(
            (nearest_other - 0.809_016_994_374_947_4).abs() < 1e-12,
            "nearest non-identity is {nearest_other}, not phi/2 = cos 36 degrees"
        );
    }

    #[test]
    fn wrong_d_is_rejected_not_silently_wrong() {
        assert!(matches!(
            Geodesic::<43>::build(),
            Err(GeoError::UnsupportedResolution(43))
        ));
    }
}
