//! The molecule the scene shows, built out of whatever elements this universe
//! turned out to have.
//!
//! **No engine type and no UI type.** This file is in the seam's strictest
//! tier, so every test here runs as plain logic with no window and no device.
//!
//! **It does arithmetic on chemistry values, and Step 4 is when that started.**
//! An earlier version of this line said it did not, which stopped being true
//! the moment colours and drawn radii were baked here. The honest rule is
//! narrower and still load-bearing: everything computed here is a **display**
//! quantity derived from chemistry values, and **none of it flows back** — no
//! `Embedding`, no `Element` and no `Universe` is altered, and nothing here is
//! read by anything that computes chemistry. The reason it lives here rather
//! than at paint time is §8.6: it is per-species work, done once when the
//! molecule is laid out, and no frame may recompute it.
//!
//! # Why the molecule is specified by valence and not by element
//!
//! Element ids mean nothing across universes — slot 40 of one seed and slot 40
//! of another are different elements — so a molecule written as a list of
//! [`ElementId`]s is a molecule that exists in one universe by luck. Measured
//! over 5000 seeds and reproduced over 500 here:
//! `ElementId::from_index(0)` has **valence 0 in every one of them**, so the
//! obvious "chain the first few elements" spelling is refused by
//! [`Mol12::add_bond`] on its first bond, 100% of the time.
//!
//! **Three corpus sizes appear in this file and they are three different
//! measurements, not one quoted inconsistently.** The shipped tests run over
//! **500** seeds, which is the number every bar here was set from. The **5000**
//! figures come from an independent run during design and are quoted where they
//! are the stronger evidence. The **2000** figures come from the Step 3 design
//! memo's own survey of `embed`'s output and are about shapes this file does not
//! build.
//!
//! # Why a centre and four leaves, and why the valences differ
//!
//! The shape has to be worth turning around, and two obvious ones are not:
//!
//! - a **chain** embeds to an exactly straight line — `axis2/axis1 = 0.000`
//!   median over 2000 seeds, which is correct output (a path graph's hop metric
//!   really is realisable in 1D) and nothing to look at;
//! - a **ring** is exactly planar, bare or substituted.
//!
//! A star is rank-3, and the *mixed* valence floor is what makes it read as one
//! body: drawing every atom from one narrow band gives five same-sized balls
//! floating apart (radius ratio 1.04, no bonded pair touching), while a
//! high-valence centre with low-valence leaves gives a centre about 2.1x its
//! leaves and **94.5% of bonded pairs interpenetrating**, with non-bonded pairs
//! never overlapping. So the molecule looks connected without Step 4's
//! cylinders, and an overlap on screen carries information rather than being an
//! artefact.
//!
//! # What this is not
//!
//! **The topology is hard-coded; the chemistry is not, and that line is the
//! whole §5 argument for this file.** A star of five is admissible because it is
//! justified by a measurement about [`embed`]'s output, and it would be a G3/G6
//! breach the moment it were justified by resembling a real molecule —
//! *including in a comment*. No guard can see this: the shipped G1 check reads
//! data files, and a molecule built in code is invisible to it. It is a review
//! obligation, written down because that is the only enforcement available.

use borbax_molecule::{CanonMol, Embedding, Mol12, canonicalise, embed};
use borbax_universe::{BondOrder, ElementId, Universe};

/// The centre needs somewhere to put four bonds.
const CENTRE_VALENCE: u8 = 4;

/// A leaf needs one, and no more is required of it.
///
/// **Deliberately not equal to [`CENTRE_VALENCE`].** A uniform floor selects
/// every atom from the same narrow band of the table, and the molecule comes out
/// as five balls of one size that do not touch.
const LEAF_VALENCE: u8 = 1;

/// How many leaves hang off the centre.
///
/// **`u8`, so the assertion below needs no cast.** It indexes nothing; the two
/// places that want a length widen it losslessly.
const LEAVES: u8 = 4;

// **The floors must differ, and this is a compile-time check rather than a test
// because the comparison is between two constants.** A uniform floor is the
// obvious simplification — one constant instead of two — and it is the one
// change that quietly turns the molecule back into five identical balls floating
// apart. `the_centre_is_visibly_larger_than_its_leaves` measures the
// consequence; this refuses the edit outright, the way `element.rs` pins
// `ElementId`'s size.
const _: () = assert!(
    LEAF_VALENCE < CENTRE_VALENCE,
    "the valence floors are equal, so every atom is drawn from one band of the \
     table and the molecule has nothing to look at"
);

// **The relation that actually empties the window**, and it had only a prose
// comment while the weaker, aesthetic one above got the compile-time check.
// Raise `LEAVES` past the centre's bonding slots and `add_bond` refuses the last
// bond in *every* universe — not some, every — so the molecule fails to build
// everywhere and the scene is blank. A review found the asymmetry: the
// load-bearing constant relation was guarded by a sentence.
const _: () = assert!(
    LEAVES <= CENTRE_VALENCE,
    "the centre has fewer bonding slots than it has leaves, so `add_bond` refuses \
     the last bond in every universe and the window opens empty"
);

/// A laid-out demo molecule, and which elements it was built from.
#[derive(Debug, Clone)]
pub struct Demo {
    /// Coordinates and radii, in **canonical** order.
    pub embedding: Embedding,
    /// The elements used, in **build** order: the centre, then its leaves.
    ///
    /// **Not positionally aligned with [`Embedding::coords`], and saying so is
    /// the point.** `canonicalise` relabels atoms — that is its entire job — so
    /// index `i` here and index `i` there are different atoms. This list is for
    /// naming what is on screen; anything positional must come off the
    /// embedding, which carries its own radii for exactly this reason.
    ///
    /// **Measured, so the hazard is a number rather than a caution:** build
    /// order and canonical order disagree in **500 of 500** universes, and the
    /// permutation is exactly the one that matters — canonical order runs
    /// `[leaf, leaf, leaf, leaf, centre]` against this list's
    /// `[centre, leaf, leaf, leaf, leaf]`. So zipping this against [`bonds`] or
    /// [`colours`] paints the centre's colour onto a leaf, in every universe.
    ///
    /// [`bonds`]: Self::bonds
    /// [`colours`]: Self::colours
    pub elements: Vec<ElementId>,
    /// Which pairs of atoms are bonded, as **canonical** indices.
    ///
    /// Each pair is ordered `[low, high]` and the list is ascending, so two
    /// runs produce the same list and a test may compare it as a sequence.
    ///
    /// **Taken from the molecule's own bond list, never inferred from how close
    /// two atoms ended up.** Proximity is not adjacency: the embedding places
    /// atoms to satisfy hop distances, so a non-bonded pair in a ring can sit
    /// closer than a bonded pair elsewhere, and a renderer that guessed would
    /// draw a bond the chemistry does not have.
    ///
    /// **No bond order travels with it, deliberately.** Every bond this file
    /// builds is [`BondOrder::SINGLE`], so an order field would be a value no
    /// consumer could tell apart from a constant — the shape this crate's
    /// guards exist to refuse. It arrives with the builder, when it can vary.
    pub bonds: Vec<[u8; 2]>,
    /// The colour each atom is drawn, in **canonical** order, aligned with
    /// [`Embedding::coords`].
    ///
    /// Computed here rather than at paint time because this is where the
    /// `&Universe` is already in hand, and because it makes the invalidation
    /// question answer itself: a colour can only change when the molecule does,
    /// which is the event the scene already rebuilds on.
    pub colours: Vec<[f64; 3]>,
    /// The radius each atom is drawn at in the sticks view, in **canonical**
    /// order, aligned with [`Embedding::coords`].
    ///
    /// **Derived from this molecule's own geometry, not from a constant.** A
    /// stick between two atoms is only visible if the two drawn spheres do not
    /// swallow it, which happens exactly when their drawn radii sum to less
    /// than the distance between them. That is a property of *this* molecule,
    /// so it is computed from it — see `stick_geometry`. A hard-coded scale
    /// would be correct on the molecules it was surveyed against and
    /// unjustified for any other, which is the whole reason this is a list
    /// rather than a number.
    ///
    /// Every radius is the true radius times **one** scale shared by all of
    /// them, so the relative sizes are exactly as the solver computed them. A
    /// per-atom adjustment would be the fabrication the scene refuses.
    pub stick_radii: Vec<f64>,
    /// How fat a bond is drawn, in chemistry-space units.
    ///
    /// A share of the **smallest drawn atom of this molecule**, so a stick
    /// always reads as thinner than everything it joins — again a ratio of a
    /// real atom rather than a constant chosen once.
    pub stick_width: f64,
}

/// Lay out this universe's demo molecule, or say there is none.
///
/// `None` means the universe has no element that could be a centre, or too few
/// that could be leaves, or the canonical search hit its cap. Measured over 5000
/// seeds during design and 500 in the shipped tests, **none of those happened**:
/// every seed built, none lacked a centre, none was short of leaves, no bond was
/// refused and nothing was capped. The branch is kept as data rather than
/// removed because "measured never" is not "cannot", and an [`Option`] costs a
/// caption line where a `panic!` would cost the window.
///
/// **This is the only [`embed`] and the only [`canonicalise`] call site in the
/// crate**, which `cargo xtask` checks rather than this comment asserting it.
/// Both are per-species work under §8.6 and must never reach a frame.
#[must_use]
pub fn build(universe: &Universe) -> Option<Demo> {
    let mut elements: Vec<ElementId> = Vec::with_capacity(usize::from(LEAVES) + 1);
    let centre_element = lowest_unused(universe, CENTRE_VALENCE, &elements)?;
    elements.push(centre_element);

    let mut mol = Mol12::new();
    let centre = mol.add_atom(centre_element)?;

    for _ in 0..LEAVES {
        let leaf_element = lowest_unused(universe, LEAF_VALENCE, &elements)?;
        elements.push(leaf_element);
        let leaf = mol.add_atom(leaf_element)?;
        // `SINGLE`, so the centre's four bonds cost exactly its four slots. A
        // higher order would spend the budget faster than the floor above
        // guarantees it, which is a refusal that depends on the universe.
        mol.add_bond(centre, leaf, BondOrder::SINGLE, &universe.table)
            .ok()?;
    }

    let (species, _search) = canonicalise(&mol).ok()?;
    let embedding = embed(&species, universe);
    let bonds = bonds_of(&species)?;
    let (stick_radii, stick_width) = stick_geometry(&embedding, &bonds);
    Some(Demo {
        colours: colours_of(&species, universe)?,
        embedding,
        bonds,
        stick_radii,
        stick_width,
        elements,
    })
}

/// How much of the separating bound the drawn atoms actually take up.
///
/// The bound itself is the largest scale at which the closest bonded pair still
/// touches; drawing *at* it leaves a gap of exactly zero. 0.6 spends 60% of it
/// on the atoms and leaves **40% of that pair's separation as visible gap**,
/// which is the quantity a person sees. Measured on the demo molecule, whose
/// worst bound across 500 universes is 0.872035: a scale of about 0.523 and a
/// minimum surface gap of about 0.47 span.
const STICK_CLEARANCE: f64 = 0.6;

/// How fat a bond is, as a share of the smallest atom it could join.
///
/// Below 1.0 by enough that a stick never reads as a neck between two lumps,
/// and far enough above zero to be visible while the molecule turns.
const STICK_SHARE: f64 = 0.38;

// A stick has to be thinner than the atoms it joins or it is not a stick. A
// compile-time check because the comparison is between two constants, and
// because the failure is a picture nobody would call wrong — just muddled.
const _: () = assert!(
    STICK_SHARE < 1.0,
    "a bond is drawn at least as fat as the thinnest atom it joins, so the      molecule reads as lumps connected by necks rather than as balls and sticks"
);

/// The drawn radii and bond width for the sticks view.
///
/// **The scale is the largest one at which every bonded pair still separates,
/// less a clearance.** Two drawn spheres of radii `s·rₐ` and `s·r_b` whose
/// centres are `d` apart stop overlapping exactly when `s < d / (rₐ + r_b)`, so
/// the binding constraint is the smallest such ratio over the molecule's bonds.
/// Taking it from the molecule rather than from a survey is what makes this
/// correct for a shape nobody has measured yet — which is every shape the
/// builder will produce.
///
/// A molecule with no bonds has no constraint and nothing to hide a stick
/// behind, so it is drawn at true size; the scale is never above 1.0, because
/// drawing an atom *larger* than it is would be a fabrication rather than a
/// display choice.
///
/// **Folded with explicit comparisons rather than `f64::min`/`f64::max`**,
/// which §13.1 disallows for disagreeing about NaN between platforms.
fn stick_geometry(embedding: &Embedding, bonds: &[[u8; 2]]) -> (Vec<f64>, f64) {
    let radii: Vec<f64> = embedding.radii().iter().map(|r| r.get()).collect();

    let mut bound = f64::INFINITY;
    for [a, b] in bonds {
        let (Some(ra), Some(rb)) = (
            radii.get(usize::from(*a)).copied(),
            radii.get(usize::from(*b)).copied(),
        ) else {
            continue;
        };
        let Some(distance) = embedding.distance(usize::from(*a), usize::from(*b)) else {
            continue;
        };
        let sum = ra + rb;
        if sum <= 0.0 {
            continue;
        }
        let ratio = distance.get() / sum;
        if ratio < bound {
            bound = ratio;
        }
    }

    // `bound` is infinite when the molecule has no bonds, which the clamp turns
    // into "draw it at true size" rather than into a non-finite scale.
    let mut scale = STICK_CLEARANCE * bound;
    // **`is_finite` first, and it is not defensive padding.** `bound` is
    // `f64::INFINITY` for a molecule with no bonds, and a `NaN` would slip past
    // a bare `scale > 1.0` — every comparison against `NaN` is false, so the
    // guard shaped to catch the bad case is the one the bad case passes
    // through. That exact shape is recorded in `scene.rs`'s viewport check.
    if !scale.is_finite() || scale > 1.0 {
        scale = 1.0;
    }

    let drawn: Vec<f64> = radii.iter().map(|r| r * scale).collect();

    let mut smallest = f64::INFINITY;
    for r in &drawn {
        if *r < smallest {
            smallest = *r;
        }
    }
    // An empty molecule leaves `smallest` infinite; a width of zero is the
    // honest answer for a molecule with no atoms to join.
    let width = if smallest.is_finite() {
        STICK_SHARE * smallest
    } else {
        0.0
    };

    (drawn, width)
}

/// Which pairs of atoms are bonded, as canonical indices.
///
/// Ordered `[low, high]` within each pair and ascending overall, so the list is
/// a function of the species rather than of the traversal.
///
/// `None` if the species holds more atoms than an index can name, which
/// [`Mol12`] makes unreachable — kept as data rather than a cast, because the
/// cast is what `clippy::as_conversions` denies and the fallible conversion is
/// what says why.
fn bonds_of(species: &CanonMol) -> Option<Vec<[u8; 2]>> {
    let n = u8::try_from(species.len()).ok()?;
    let mut bonds = Vec::new();
    for a in 0..n {
        for b in (a + 1)..n {
            // **`bond_order`, not a distance test.** See [`Demo::bonds`].
            if species.mol().bond_order(a, b).is_some() {
                bonds.push([a, b]);
            }
        }
    }
    Some(bonds)
}

/// The colour of each atom, in canonical order.
///
/// The mass range is the universe's own attained spread, so lightness uses the
/// contrast that universe actually has rather than one borrowed from another.
///
/// `None` if a canonical index names no atom or an atom names no element,
/// neither of which [`canonicalise`] can produce — the [`Option`] costs a
/// caption line where a `panic!` would cost the window.
fn colours_of(species: &CanonMol, universe: &Universe) -> Option<Vec<[f64; 3]>> {
    let (mass_lo, mass_hi) = mass_range(universe);
    let n = u8::try_from(species.len()).ok()?;
    let mut colours = Vec::with_capacity(usize::from(n));
    for i in 0..n {
        let id = species.mol().element(i)?;
        let element = universe.table.get(id)?;
        colours.push(crate::palette::colour(
            element.affinity,
            element.valence,
            element.mass.to_f64(),
            mass_lo,
            mass_hi,
            element.group,
        ));
    }
    Some(colours)
}

/// The lightest and heaviest masses this universe's table attains.
///
/// **Folded with explicit comparisons rather than `f64::min`/`f64::max`**,
/// which §13.1 disallows for disagreeing about NaN between platforms.
///
/// A table with no elements gives a degenerate range, which the palette reads
/// as "no spread to show" rather than dividing by zero.
fn mass_range(universe: &Universe) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for (_, element) in universe.table.iter() {
        let mass = element.mass.to_f64();
        if mass < lo {
            lo = mass;
        }
        if mass > hi {
            hi = mass;
        }
    }
    (lo, hi)
}

/// The lowest-indexed element with at least `floor` bonding slots that `used`
/// does not already hold.
///
/// **"Lowest" is a rule, not a search, and that is the point.** Anything that
/// ranked candidates — most bonds, biggest radius, strongest affinity — would be
/// the viewer choosing which chemistry is interesting, which is a judgement it
/// is not entitled to make. Taking the first one that fits takes no view at all.
///
/// Worth knowing before it looks like a bug: across 5000 seeds the centre lands
/// on just six distinct indices, never below 24. That is the shell law showing
/// through — the early elements have small clusters and few frontier contacts —
/// not a stuck iterator.
fn lowest_unused(universe: &Universe, floor: u8, used: &[ElementId]) -> Option<ElementId> {
    universe
        .table
        .iter()
        .find(|(id, e)| e.valence >= floor && !used.contains(id))
        .map(|(id, _)| id)
}

#[cfg(test)]
mod measure {
    use super::{CENTRE_VALENCE, LEAVES, build};
    use borbax_universe::Universe;

    /// Not an assertion — a reading, printed so the bars below are chosen from
    /// measurement rather than from taste.
    #[test]
    #[ignore = "measurement, not a check"]
    fn survey() {
        let mut built = 0usize;
        let mut ratio_lo = f64::INFINITY;
        let mut ratio_hi: f64 = 0.0;
        let mut centres = std::collections::BTreeSet::new();
        let mut distinct_ok = 0usize;
        let mut gyr_lo = f64::INFINITY;
        let mut gyr_hi: f64 = 0.0;
        let seeds = 500u64;
        for seed in 0..seeds {
            let u = Universe::generate(seed);
            let Some(demo) = build(&u) else { continue };
            built += 1;
            let centre = demo.elements.first().copied();
            centres.insert(u.table.iter().position(|(id, _)| Some(id) == centre));
            let mut set = demo.elements.clone();
            set.sort_unstable();
            set.dedup();
            if set.len() == usize::from(LEAVES) + 1 {
                distinct_ok += 1;
            }
            let radii: Vec<f64> = demo.embedding.radii().iter().map(|s| s.get()).collect();
            let mut lo = f64::INFINITY;
            let mut hi = 0.0f64;
            for r in &radii {
                if *r < lo {
                    lo = *r;
                }
                if *r > hi {
                    hi = *r;
                }
            }
            let r = hi / lo;
            if r < ratio_lo {
                ratio_lo = r;
            }
            if r > ratio_hi {
                ratio_hi = r;
            }
            let g = demo.embedding.radius_of_gyration().get();
            if g < gyr_lo {
                gyr_lo = g;
            }
            if g > gyr_hi {
                gyr_hi = g;
            }
        }
        println!("seeds={seeds} built={built} distinct_ok={distinct_ok}");
        println!("radius ratio hi/lo: {ratio_lo:.4} .. {ratio_hi:.4}");
        println!("gyration: {gyr_lo:.4} .. {gyr_hi:.4}");
        println!("centre element indices ({}): {centres:?}", centres.len());
        println!("centre valence floor = {CENTRE_VALENCE}");
    }
}

#[cfg(test)]
mod shapes {
    //! Bond extraction, tested on molecules the demo molecule cannot stand in
    //! for.
    //!
    //! **The demo molecule discriminates nothing here, and that is measured
    //! rather than suspected.** Its bond set is the *same constant* in 500 of
    //! 500 universes — `[(0,4), (1,4), (2,4), (3,4)]` — so the correct answer
    //! agrees with "every pair containing the last canonical index", "every
    //! pair containing the highest-degree atom", "every pair containing the
    //! heaviest atom", and a hard-coded literal, in every universe there is. No
    //! seed separates them.
    //!
    //! So these fixtures are built here: a path, a ring, a two-level tree and a
    //! graph with two hubs. Each is checked against an oracle that does not
    //! restate the implementation — the degree sequence and the profile of
    //! degree-pairs over edges, both invariant under the relabelling
    //! `canonicalise` performs — and each additionally asserts its own answer
    //! *differs* from the star rule's, so a fixture that quietly degenerated
    //! into a star could not leave this test passing vacuously.

    use super::{bonds_of, lowest_unused};
    use borbax_molecule::{CanonMol, Mol12, canonicalise};
    use borbax_universe::{BondOrder, ElementId, Universe};

    /// Build a molecule with `n` atoms and these edges, picking elements with
    /// enough bonding slots for the degree each atom needs.
    ///
    /// Returns `None` if the universe cannot supply them, which the caller
    /// treats as a failed fixture rather than a passing test.
    fn shape(universe: &Universe, n: u8, edges: &[[u8; 2]]) -> Option<CanonMol> {
        let mut used: Vec<ElementId> = Vec::new();
        let mut mol = Mol12::new();
        let mut atoms = Vec::new();
        for atom in 0..n {
            // Enough bonding slots for the degree this atom is about to need.
            let element = lowest_unused(universe, degree_of(atom, edges), &used)?;
            used.push(element);
            atoms.push(mol.add_atom(element)?);
        }
        for [a, b] in edges {
            // `get`, not `[]` — `clippy::indexing_slicing` is denied across the
            // workspace and reaches inside `#[cfg(test)]`. A fixture naming an
            // atom it never built should fail as a `None` here, not as a panic
            // with no fixture name in it.
            let from = atoms.get(usize::from(*a)).copied()?;
            let to = atoms.get(usize::from(*b)).copied()?;
            mol.add_bond(from, to, BondOrder::SINGLE, &universe.table)
                .ok()?;
        }
        canonicalise(&mol).ok().map(|(species, _)| species)
    }

    /// How many edges touch this atom.
    ///
    /// Counted by scanning rather than by accumulating into an indexed table,
    /// which is what keeps these helpers free of slice indexing.
    fn degree_of(atom: u8, edges: &[[u8; 2]]) -> u8 {
        let touching = edges
            .iter()
            .filter(|[a, b]| *a == atom || *b == atom)
            .count();
        u8::try_from(touching).unwrap_or(u8::MAX)
    }

    /// The degree of each atom, ascending — invariant under relabelling.
    fn degree_profile(n: u8, bonds: &[[u8; 2]]) -> Vec<u8> {
        let mut degrees: Vec<u8> = (0..n).map(|atom| degree_of(atom, bonds)).collect();
        degrees.sort_unstable();
        degrees
    }

    /// The `(lower degree, higher degree)` of each edge's endpoints, ascending.
    ///
    /// **Strictly stronger than the degree sequence alone**, which two
    /// non-isomorphic graphs can share. This separates shapes the degree
    /// sequence cannot.
    fn edge_profile(bonds: &[[u8; 2]]) -> Vec<(u8, u8)> {
        let mut profile: Vec<(u8, u8)> = bonds
            .iter()
            .map(|[a, b]| {
                let (x, y) = (degree_of(*a, bonds), degree_of(*b, bonds));
                if x <= y { (x, y) } else { (y, x) }
            })
            .collect();
        profile.sort_unstable();
        profile
    }

    /// What "every pair containing the highest canonical index" would return.
    ///
    /// The wrong implementation the demo molecule cannot rule out. Each fixture
    /// asserts the real answer differs from this one, which is what makes the
    /// fixture's own discriminating power visible rather than assumed.
    fn star_rule(n: u8) -> Vec<[u8; 2]> {
        (0..n.saturating_sub(1)).map(|a| [a, n - 1]).collect()
    }

    /// The bond list is the molecule's own bonds, on shapes that can tell.
    ///
    /// **Discriminator:** replacing `bonds_of`'s body with [`star_rule`] must
    /// fail here. It passes on the demo molecule in every universe, so no other
    /// test in this crate would notice.
    #[test]
    fn the_bond_list_is_the_molecules_own_bonds() {
        let universe = Universe::generate(1);

        let mut checked = 0_u32;
        for fixture in FIXTURES {
            let Fixture {
                name,
                atoms,
                edges,
                degrees,
            } = fixture;
            let species = shape(&universe, *atoms, edges)
                .unwrap_or_else(|| unreachable!("seed 1 supplies the {name} fixture"));
            let bonds =
                bonds_of(&species).unwrap_or_else(|| unreachable!("a fixture fits in a u8 index"));

            assert_eq!(
                bonds.len(),
                edges.len(),
                "the {name} fixture has {} bonds against {} built",
                bonds.len(),
                edges.len()
            );
            assert_eq!(
                degree_profile(*atoms, &bonds),
                *degrees,
                "the {name} fixture's degree sequence is not the one it was built with"
            );
            assert_eq!(
                edge_profile(&bonds),
                edge_profile(edges),
                "the {name} fixture's edges do not join the same kinds of atom"
            );

            // The corpus-filter check: a fixture that degenerated into a star
            // would satisfy every assertion above while proving nothing.
            assert_ne!(
                bonds,
                star_rule(*atoms),
                "the {name} fixture *is* the star rule's answer, so it cannot \
                 discriminate a real bond lookup from that mutation"
            );
            checked += 1;
        }
        assert_eq!(
            checked,
            u32::try_from(FIXTURES.len()).unwrap_or(u32::MAX),
            "not every shape was checked"
        );
    }

    /// One shape to build and what it should come back as.
    ///
    /// A named struct rather than a tuple, because
    /// `clippy::type_complexity` denies the four-field slice type and — more
    /// usefully — `degrees` and `atoms` are both integers that would otherwise
    /// be told apart only by position.
    struct Fixture {
        /// What to call it when the assertion fails.
        name: &'static str,
        /// How many atoms to build.
        atoms: u8,
        /// Which pairs to bond, as build indices.
        edges: &'static [[u8; 2]],
        /// The ascending degree sequence the result must have.
        degrees: &'static [u8],
    }

    /// Shapes the demo molecule cannot stand in for.
    ///
    /// Each has a degree sequence unlike a star's `[1, 1, 1, 3]`, which is what
    /// gives it the power to fail under the star-rule mutation.
    const FIXTURES: &[Fixture] = &[
        Fixture {
            name: "path",
            atoms: 4,
            edges: &[[0, 1], [1, 2], [2, 3]],
            degrees: &[1, 1, 2, 2],
        },
        Fixture {
            name: "ring",
            atoms: 4,
            edges: &[[0, 1], [1, 2], [2, 3], [3, 0]],
            degrees: &[2, 2, 2, 2],
        },
        Fixture {
            name: "two-level tree",
            atoms: 6,
            edges: &[[0, 1], [0, 2], [1, 3], [1, 4], [2, 5]],
            degrees: &[1, 1, 1, 2, 2, 3],
        },
        Fixture {
            name: "two hubs",
            atoms: 6,
            edges: &[[0, 1], [0, 2], [0, 3], [1, 4], [1, 5]],
            degrees: &[1, 1, 1, 1, 3, 3],
        },
    ];

    /// A bond list is ascending and each pair is ordered.
    ///
    /// Not cosmetic: the scene compares stick endpoints as a set keyed on these
    /// pairs, and an unordered pair would make `[a, b]` and `[b, a]` two
    /// different bonds.
    #[test]
    fn the_bond_list_is_ordered() {
        let universe = Universe::generate(1);
        let species = shape(&universe, 4, &[[2, 3], [0, 1], [1, 2]])
            .unwrap_or_else(|| unreachable!("seed 1 supplies a path"));
        let bonds =
            bonds_of(&species).unwrap_or_else(|| unreachable!("a fixture fits in a u8 index"));
        for [a, b] in &bonds {
            assert!(a < b, "the pair [{a}, {b}] is not ordered");
        }
        let mut sorted = bonds.clone();
        sorted.sort_unstable();
        assert_eq!(bonds, sorted, "the bond list is not ascending");
    }

    /// Bonds come from the graph, not from how close two atoms ended up.
    ///
    /// **A ring is the shape that separates them.** In a 4-ring every atom is
    /// bonded to two others and *not* bonded to the one across the diagonal —
    /// so a proximity rule with any threshold either misses a real bond or
    /// invents the diagonal. Asserting the diagonal is absent is what says the
    /// list is adjacency rather than nearness.
    #[test]
    fn a_bond_is_adjacency_and_not_nearness() {
        let universe = Universe::generate(1);
        let species = shape(&universe, 4, &[[0, 1], [1, 2], [2, 3], [3, 0]])
            .unwrap_or_else(|| unreachable!("seed 1 supplies a ring"));
        let bonds =
            bonds_of(&species).unwrap_or_else(|| unreachable!("a fixture fits in a u8 index"));
        assert_eq!(
            bonds.len(),
            4,
            "a 4-ring has four bonds, not {}",
            bonds.len()
        );
        // Every atom has exactly two neighbours; a diagonal would give one three.
        assert_eq!(
            degree_profile(4, &bonds),
            vec![2, 2, 2, 2],
            "a ring atom has two neighbours — a diagonal drawn from proximity \
             would give some atom three"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::{LEAVES, build};
    use borbax_universe::Universe;

    /// How many universes every claim below is made over.
    ///
    /// Enumerable spaces get enumerated; this one is not, so the corpus is
    /// stated once and every count is asserted as an **equality** against it. A
    /// `> 0` bar passes on a corpus that silently generated one seed, which is
    /// the shape that has hidden a dead test here before.
    const SEEDS: u64 = 500;

    /// The molecule builds in every universe it is offered.
    ///
    /// **The discriminator is a hard-coded `ElementId`.** Ids mean different
    /// elements in different universes, and slot 0 has valence 0 in every one
    /// measured — so a molecule written as a list of ids is refused by
    /// `add_bond` on its first bond. Selection by valence is what this pins.
    #[test]
    fn the_demo_molecule_is_buildable_in_every_universe_it_is_shown_in() {
        let mut built = 0u64;
        for seed in 0..SEEDS {
            if build(&Universe::generate(seed)).is_some() {
                built += 1;
            }
        }
        assert_eq!(
            built,
            SEEDS,
            "the demo molecule failed to build in {} of {SEEDS} universes, so the \
             window would open on an empty scene for anyone who typed one of them",
            SEEDS - built
        );
    }

    /// Five atoms, no more and no fewer.
    ///
    /// Fails on an off-by-one in the leaf loop, which is invisible on screen for
    /// a symmetric molecule and invisible to every other test here.
    #[test]
    fn the_molecule_is_a_centre_and_four_leaves() {
        let demo = build(&Universe::generate(1))
            .unwrap_or_else(|| unreachable!("seed 1 builds — the test above says so over 500"));
        assert_eq!(
            demo.elements.len(),
            usize::from(LEAVES) + 1,
            "wrong number of elements"
        );
        assert_eq!(
            demo.embedding.len(),
            usize::from(LEAVES) + 1,
            "the embedding does not hold one point per atom"
        );
    }

    /// The centre is substantially larger than its leaves, in every universe.
    ///
    /// **This is the test that pins the *mixed* valence floor**, and it is the
    /// whole reason the molecule reads as one object rather than five. Measured
    /// over 500 universes the largest-to-smallest radius ratio runs
    /// **2.0609 .. 2.8523**. Setting `LEAF_VALENCE` equal to `CENTRE_VALENCE` —
    /// a *uniform* floor, which is the obvious simplification — draws every atom
    /// from one narrow band of the table and collapses that ratio to about 1.04:
    /// five identical balls, floating apart, with nothing to look at.
    ///
    /// The bar is 1.5 because it sits between the two measured populations with
    /// room on both sides, not because 1.5 made the current code pass.
    #[test]
    fn the_centre_is_visibly_larger_than_its_leaves() {
        let mut checked = 0u64;
        let mut worst = f64::INFINITY;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds"));
            let mut lo = f64::INFINITY;
            let mut hi = 0.0f64;
            for radius in demo.embedding.radii() {
                let r = radius.get();
                assert!(
                    r > 0.0,
                    "an atom has radius {r} at seed {seed}, so it is a point"
                );
                if r < lo {
                    lo = r;
                }
                if r > hi {
                    hi = r;
                }
            }
            let ratio = hi / lo;
            if ratio < worst {
                worst = ratio;
            }
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this bar was measured on"
        );
        assert!(
            worst > 1.5,
            "the widest-to-narrowest radius ratio fell to {worst:.4}, against 2.0609 \
             measured across {SEEDS} universes. A uniform valence floor gives about \
             1.04 — five same-sized balls that do not touch"
        );
    }

    /// No element is used twice.
    ///
    /// Fails on dropping the "not already used" clause from the selection, which
    /// makes all four leaves the same element: still buildable, still a star,
    /// and every leaf identical.
    #[test]
    fn each_atom_is_a_different_element() {
        let mut checked = 0u64;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds"));
            let mut ids = demo.elements.clone();
            ids.sort_unstable();
            ids.dedup();
            assert_eq!(
                ids.len(),
                usize::from(LEAVES) + 1,
                "seed {seed} reused an element, so its leaves are indistinguishable"
            );
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
    }

    /// Every atom is coloured from the element that atom actually is.
    ///
    /// **This is the test for the defect the canonical relabelling guarantees.**
    /// `Demo::elements` is build order and the embedding is canonical order, and
    /// they disagree in 500 of 500 universes — canonical runs
    /// `[leaf, leaf, leaf, leaf, centre]` against build's
    /// `[centre, leaf, leaf, leaf, leaf]`. So the obvious implementation, zipping
    /// `elements` against `coords`, paints the centre's colour onto a leaf and a
    /// leaf's onto the big centre sphere, in every universe there is.
    ///
    /// **The oracle is the radius, and it is deliberately not the index.** A
    /// `want` list built by walking canonical indices the same way `colours_of`
    /// does would reproduce whatever `colours_of` did, including the defect —
    /// which is exactly how this class of test has shipped green here before.
    /// `Embedding::radii()[i]` *is* the element's own `radius`, and every demo
    /// radius is unique across the whole table in 500/500 universes, so the
    /// radius identifies the element on screen with no reference to any index at
    /// all. The uniqueness is asserted rather than trusted, so this test reports
    /// the oracle failing instead of silently matching the wrong element.
    #[test]
    fn every_atom_is_coloured_from_its_own_element() {
        let mut checked = 0_u64;
        let mut atoms_checked = 0_u64;
        for seed in 0..SEEDS {
            let universe = Universe::generate(seed);
            let demo =
                build(&universe).unwrap_or_else(|| unreachable!("every seed builds a molecule"));

            let (mass_lo, mass_hi) = super::mass_range(&universe);
            assert_eq!(
                demo.colours.len(),
                demo.embedding.len(),
                "seed {seed} has {} colours for {} atoms",
                demo.colours.len(),
                demo.embedding.len()
            );

            for (i, radius) in demo.embedding.radii().iter().enumerate() {
                // The oracle: which elements of the whole table have this radius?
                let matches: Vec<_> = universe
                    .table
                    .iter()
                    .filter(|(_, e)| e.radius.get().to_bits() == radius.get().to_bits())
                    .map(|(_, e)| e)
                    .collect();
                assert_eq!(
                    matches.len(),
                    1,
                    "seed {seed}: radius {} names {} elements, so it cannot \
                     identify the atom at canonical index {i} and this test's \
                     oracle has stopped working",
                    radius.get(),
                    matches.len()
                );
                let element = matches
                    .first()
                    .unwrap_or_else(|| unreachable!("exactly one match, asserted above"));

                let want = crate::palette::colour(
                    element.affinity,
                    element.valence,
                    element.mass.to_f64(),
                    mass_lo,
                    mass_hi,
                    element.group,
                );
                let got = demo
                    .colours
                    .get(i)
                    .unwrap_or_else(|| unreachable!("one colour per atom, asserted above"));
                assert_eq!(
                    got.map(f64::to_bits),
                    want.map(f64::to_bits),
                    "seed {seed}: the atom at canonical index {i} has radius {} \
                     — which is element {:?} — but is painted a colour that \
                     element does not earn. Build order and canonical order \
                     disagree in every universe, so this is what a `zip` of \
                     `elements` against `coords` produces",
                    radius.get(),
                    element.symbol
                );
                atoms_checked += 1;
            }
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
        // Five atoms per universe, so the count is an equality rather than a
        // floor: a loop that stopped selecting atoms would otherwise be silent.
        assert_eq!(
            atoms_checked,
            SEEDS * (u64::from(LEAVES) + 1),
            "the number of atoms examined is not five per universe"
        );
    }

    /// Every bond names two atoms that exist, and no bond is a self-loop.
    ///
    /// Cheap, and it is what stops an out-of-range index reaching the scene,
    /// where it would index a coordinate list and take the window with it.
    #[test]
    fn every_bond_names_two_atoms_of_this_molecule() {
        let mut checked = 0_u64;
        let mut bonds_checked = 0_u64;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            let n = u8::try_from(demo.embedding.len())
                .unwrap_or_else(|_| unreachable!("a molecule fits in a u8 index"));
            for [a, b] in &demo.bonds {
                assert!(
                    *a < n && *b < n,
                    "seed {seed}: bond [{a}, {b}] names an atom outside a \
                     molecule of {n}"
                );
                assert_ne!(a, b, "seed {seed}: bond [{a}, {b}] is a self-loop");
                bonds_checked += 1;
            }
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
        // The demo molecule is a centre and four leaves, so it has exactly
        // `LEAVES` bonds — an equality, not a floor.
        assert_eq!(
            bonds_checked,
            SEEDS * u64::from(LEAVES),
            "the number of bonds examined is not one per leaf per universe"
        );
    }

    /// The sticks view separates every bonded pair, in every universe.
    ///
    /// **This is the test that says the sticks view does its job**, and it is
    /// the one that would fail for the obvious implementation. Step 3 made
    /// interpenetration the way the molecule reads as one body — 94.7% of
    /// bonded pairs overlap at true size — so a stick drawn between two
    /// true-size atoms is buried inside them and shows on about 9% of bonds.
    /// Every structural test passes over that: the entity count is right, the
    /// transforms are right, and nothing is visible.
    ///
    /// So the assertion is the *geometric* one — two drawn spheres do not
    /// reach each other — rather than a count of entities.
    #[test]
    fn the_sticks_view_separates_every_bonded_pair() {
        let mut checked = 0_u64;
        let mut bonds_checked = 0_u64;
        let mut tightest = f64::INFINITY;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            for [a, b] in &demo.bonds {
                let ra = demo
                    .stick_radii
                    .get(usize::from(*a))
                    .copied()
                    .unwrap_or_else(|| unreachable!("one drawn radius per atom"));
                let rb = demo
                    .stick_radii
                    .get(usize::from(*b))
                    .copied()
                    .unwrap_or_else(|| unreachable!("one drawn radius per atom"));
                let distance = demo
                    .embedding
                    .distance(usize::from(*a), usize::from(*b))
                    .unwrap_or_else(|| unreachable!("a bond names two atoms"))
                    .get();
                let gap = distance - (ra + rb);
                assert!(
                    gap > 0.0,
                    "seed {seed}: the atoms of bond [{a}, {b}] still overlap by \
                     {:.6} in the sticks view, so the stick between them is \
                     buried inside them and the view shows nothing new",
                    -gap
                );
                if gap < tightest {
                    tightest = gap;
                }
                bonds_checked += 1;
            }
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
        assert_eq!(
            bonds_checked,
            SEEDS * u64::from(LEAVES),
            "the number of bonds examined is not one per leaf per universe"
        );
        // Reported rather than asserted against a chosen bar: the assertion
        // above is the requirement, and this says how much room it had.
        assert!(
            tightest.is_finite(),
            "no bond was measured, so the loop above asserted nothing"
        );
    }

    /// A molecule with no bonds is drawn at true size, not at infinity.
    ///
    /// **This test exists because a mutation probe found the clamp it guards to
    /// be vacuous without it.** The separating bound is a minimum over bonded
    /// pairs, so a molecule with *no* bonds leaves it at `f64::INFINITY` and the
    /// scale would be infinite. The clamp catches that — but the demo molecule
    /// can never reach it: its bound runs about 0.87 to 1.16 across 500
    /// universes, so the scale never approaches 1.0 and deleting the clamp
    /// outright left the whole suite green.
    ///
    /// Calling [`stick_geometry`] directly with an empty bond list is what makes
    /// the branch reachable, and it is a real case rather than a contrived one:
    /// the builder can produce a lone atom, and a `NaN` scale would size every
    /// sphere to nothing.
    #[test]
    fn a_molecule_with_no_bonds_is_drawn_at_true_size() {
        let demo = build(&Universe::generate(1))
            .unwrap_or_else(|| unreachable!("seed 1 builds a molecule"));
        let (drawn, width) = super::stick_geometry(&demo.embedding, &[]);

        assert_eq!(
            drawn.len(),
            demo.embedding.len(),
            "one drawn radius per atom, even with nothing to separate"
        );
        for (got, want) in drawn.iter().zip(demo.embedding.radii()) {
            assert!(
                got.is_finite(),
                "an unbonded atom is drawn at {got}, so the infinite bound \
                 reached the scale"
            );
            assert!(
                (got - want.get()).abs() < 1e-12,
                "an unbonded atom is drawn at {got:.6} against a true radius of \
                 {:.6} — with no bond to hide, there is nothing to shrink for",
                want.get()
            );
        }
        assert!(
            width.is_finite() && width > 0.0,
            "the bond width is {width}, which is not a size"
        );
    }

    /// A bond is drawn thinner than every atom it could join.
    ///
    /// Fails if the stick width stops being a share of the *smallest* drawn
    /// atom — a width taken from the average or the largest gives a stick fatter
    /// than the leaves, and the molecule reads as lumps joined by necks.
    #[test]
    fn a_stick_is_thinner_than_every_atom_it_joins() {
        let mut checked = 0_u64;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            for (i, radius) in demo.stick_radii.iter().enumerate() {
                assert!(
                    demo.stick_width < *radius,
                    "seed {seed}: a bond is drawn at {:.6} against atom {i} at \
                     {radius:.6}, so the stick is not thinner than the atoms",
                    demo.stick_width
                );
            }
            assert!(
                demo.stick_width > 0.0,
                "seed {seed}: bonds are drawn with no width at all"
            );
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
    }

    /// No atom is ever drawn larger than it is.
    ///
    /// The scale is a *display* choice and shrinking is one; enlarging would be
    /// a fabrication. Guards the clamp, which is also what stops a molecule
    /// with no bonds — where the separating bound is infinite — producing a
    /// non-finite scale.
    #[test]
    fn no_atom_is_drawn_larger_than_it_is() {
        let mut checked = 0_u64;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            assert_eq!(
                demo.stick_radii.len(),
                demo.embedding.len(),
                "seed {seed} has {} drawn radii for {} atoms",
                demo.stick_radii.len(),
                demo.embedding.len()
            );
            for (drawn, true_radius) in demo.stick_radii.iter().zip(demo.embedding.radii()) {
                assert!(
                    drawn.is_finite() && *drawn > 0.0,
                    "seed {seed}: an atom is drawn at {drawn}, which is not a size"
                );
                assert!(
                    *drawn <= true_radius.get(),
                    "seed {seed}: an atom is drawn at {drawn:.6} against a true \
                     radius of {:.6} — the view may shrink an atom and may not \
                     invent one",
                    true_radius.get()
                );
            }
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
    }

    /// The drawn radii are the true radii times **one** shared scale.
    ///
    /// **This is the guard against a per-atom fudge**, which is the fabrication
    /// the scene's own header refuses and which no separation test can see: an
    /// implementation that shrank only the two atoms of the tightest bond would
    /// satisfy every assertion above while destroying the relative sizes the
    /// solver computed.
    #[test]
    fn every_atom_is_shrunk_by_the_same_amount() {
        let mut checked = 0_u64;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds a molecule"));
            let mut ratios = demo
                .stick_radii
                .iter()
                .zip(demo.embedding.radii())
                .map(|(drawn, r)| drawn / r.get());
            let first = ratios
                .next()
                .unwrap_or_else(|| unreachable!("the molecule has atoms"));
            for ratio in ratios {
                assert!(
                    (ratio - first).abs() < 1e-12,
                    "seed {seed}: atoms are shrunk by {first:.9} and {ratio:.9}, \
                     so the relative sizes the solver computed have been changed"
                );
            }
            checked += 1;
        }
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this claim was made over"
        );
    }

    /// The molecule has a size the camera can frame.
    ///
    /// Guards the degenerate case the framing constant divides no attention to:
    /// a gyration radius of zero puts the camera at the origin, inside the
    /// scene. Measured range across 500 universes: **0.8047 .. 1.5764**.
    #[test]
    fn the_molecule_has_a_positive_extent() {
        let mut checked = 0_u64;
        for seed in 0..SEEDS {
            let demo = build(&Universe::generate(seed))
                .unwrap_or_else(|| unreachable!("every seed builds"));
            let gyration = demo.embedding.radius_of_gyration().get();
            assert!(
                gyration > 0.5,
                "seed {seed} has gyration radius {gyration:.4}, against a measured \
                 minimum of 0.8047 — a camera framed on that stands too close"
            );
            checked += 1;
        }
        // The equality this file requires of every count: a `> 0` bar passes on
        // a corpus that silently generated one seed. This was the one test here
        // not carrying it.
        assert_eq!(
            checked, SEEDS,
            "the corpus is not the one this bar was measured on"
        );
    }
}
