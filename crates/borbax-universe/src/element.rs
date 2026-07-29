//! The generated periodic table (spec §7.1).
//!
//! **Element `N` is `N` copies of one base unit, packed.** Every property below
//! is a function of that packing rather than a curve chosen because it looked
//! plausible — which is principle 2 applied to the file that defines what an
//! atom is. See Task 4's preamble in the plan for what is claimed and, more
//! importantly, what is not: there is no shape-diversity advantage, only the
//! measured finding that a fusion-derived radius series costs nothing.
//!
//! Periodicity is a *consequence* here. A shell closes when it fills, and the
//! elements at and around a closure form families with a cause: a closed shell
//! has no frontier and so no bonding slots, and the element one unit past any
//! closure has exactly one.

use crate::naming;
use crate::packing::{self, PackingConsts};
use borbax_rng::{Domain, Stream};
use borbax_units::{Mass, Span, det_math};

/// Index into a universe's element table. `u8` because tables are capped
/// well below 256, which keeps `Mol12`'s atom array one byte per slot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ElementId(pub u8);

/// One generated element: a cluster of [`Element::units`] base units, and the
/// properties that follow from how they pack.
#[derive(Debug, Clone, PartialEq)]
pub struct Element {
    /// Index into the universe's table.
    pub id: ElementId,
    /// Generated symbol, 1–2 characters (§5, G2).
    pub symbol: String,
    /// Generated name (§5, G2).
    pub name: String,
    /// Base units in this element's cluster. The element *is* this number.
    pub units: usize,
    /// Which shell is filling. The period index, derived rather than drawn.
    pub period: u8,
    /// Units in the incomplete outer shell. The group index.
    pub group: u8,
    /// Cluster mass: units carried, less the mass defect of made contacts.
    pub mass: Mass,
    /// Bonding slots: unmade lateral contacts on the frontier, over the
    /// lateral coordination. Zero at a closure because a closed shell has no
    /// frontier — not because a branch says so.
    pub valence: u8,
    /// Exposed fraction of the cluster, mapped to [-1, +1]. Drives
    /// complementarity in §8.3. **Not** electronegativity: different range,
    /// different meaning, different behaviour (G3).
    pub affinity: f64,
    /// Enclosing radius of the packed cluster. Rises monotonically with the
    /// occupied shell and never jumps — the opposite of a real atomic radius,
    /// which falls across a period (G3).
    pub radius: Span,
    /// Binding energy per unit: made contacts, less radial strain.
    ///
    /// **Task 5 must derive `bond_energies` from this**, not draw an
    /// independent matrix. Two encodings of "these elements bind well" in
    /// incommensurable units is the defect that made the predecessor's `peak`
    /// unfalsifiable, one level up.
    ///
    /// **G4 note, flagged deliberately rather than silently accepted.** This is
    /// an energy travelling as a bare `f64` into a task that mints it as
    /// `Quanta` at the far end, with no type-level relationship — the exact
    /// "bare `f64` carrying a quantity through three functions defeats this
    /// quietly" shape §5's G4 exists to stop. `borbax-units` provides
    /// `Div<f64> for Quanta`, so cluster-energy-per-unit is expressible as
    /// `Quanta` and this *could* be typed. It is left bare for now because
    /// `energy_per_unit` is also compared, summed and scaled inside this file
    /// where the newtype buys nothing, and changing it is Task 5's call when it
    /// writes the consumer. **Task 5 must decide explicitly**, not inherit the
    /// `f64` by default — that decision is what sets the precedent for every
    /// energy that follows.
    pub energy_per_unit: f64,
    /// Decay probability per world-year, from the distance below the binding
    /// peak. **Named for what it holds** — it is 0 at the peak and rises toward
    /// the extremes. An earlier draft called this `stability`, which inverted
    /// the sense of the field against its own description, in the quantity §9.4
    /// calls the most sensitive parameter in the system. Unstable elements are
    /// both an energy source and a mutation source (spec §9.5).
    pub decay_rate: f64,
    /// Relative abundance. **The gel lever** (§7.2) — see [`generate_elements`].
    pub abundance: f64,
    /// Which band of outer-shell occupancy this element exposes — `floor(fill ×
    /// n_catalytic)`.
    ///
    /// **Not derived, and this doc is what `cargo doc` renders.** An earlier
    /// version said "coordination class of the frontier sites"; no coordination
    /// mix is computed anywhere, and enumerated across all `k`, `n_catalytic`
    /// and `N` this takes exactly one value per `(outer, cap)` pair — zero
    /// information beyond `group`. Nothing downstream may branch on it as
    /// though it carried chemistry.
    pub catalytic_class: u8,
}

/// The derived shape of one universe's table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShellPattern {
    /// The drawn part of the shell law: shell `n` holds `k*n^2 + 2` units.
    pub k: usize,
    /// Cumulative unit counts at which a shell closes. Derived from `k`, not
    /// drawn — this replaces the predecessor's drawn table of period lengths.
    pub closures: Vec<usize>,
    /// Where the per-unit binding energy peaks. **Read off the finished
    /// series, never passed in** — that is the property distinguishing this
    /// from the predecessor, in which `peak` was an argument whose value was
    /// then reported back as an emergent minimum.
    pub peak: usize,
}

/// Convenience for tests and callers that need the universe-domain stream.
#[must_use]
pub const fn stream(seed: u64) -> Stream {
    Stream::new(seed, Domain::Universe, 0)
}

/// Build a universe's periodic table by fusion.
///
/// **There is no `peak` parameter and there must never be one.** The
/// predecessor took one, computed instability as `((N - peak)/peak)^2`, and
/// reported the resulting minimum as emergent — it was the distance from a
/// typed-in constant, and sweeping `peak` landed the minimum on a closure in
/// 13 of 24 draws. Here one binding energy per unit is derived from the
/// packing counts and both the mass defect and the decay rate are functions
/// of it.
///
/// Infallible: the mass is built in sub-units and handed to `Mass::from_raw`,
/// so there is no range check to fail and no error to propagate. See the task
/// preamble for why a draft made it fallible and why that was wrong.
#[must_use]
#[expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    reason = "every narrowing here is bounded by construction and the bound is the \
              point: `out.len()` and `units` are < `n_elements` <= 120; `shell` <= 3 \
              because the fourth shell opens past N = 189 for every drawn `k`; \
              `outer` < `cap` <= 226; `valence` is `frontier_notches` rounded, pinned \
              to 4..=6 by `valence_ceiling_is_four_to_six_over_the_drawn_range`; \
              `catalytic_class` < `n_catalytic` <= 6; the sub-unit products are \
              exact integers well under 2^31 (see \
              `every_mass_is_inside_the_range_from_raw_does_not_check`)"
)]
pub fn generate_elements(seed: u64) -> (ShellPattern, Vec<Element>) {
    let mut rng = stream(seed);

    // `k` spans 6..=14. Below 6 a shell cannot triangulate a sphere; above 14
    // the second shell exceeds 200 units, so a 120-element table never leaves
    // it and there are no families at all.
    let k = 6 + rng.next_range(9) as usize;
    // `z = k + 2` is forced by the shell law, not chosen: `shell(1) = k + 2` and
    // shell 1 *is* the units touching the core. An earlier `6 + k/2` agreed only
    // at k = 8 and put 16 sites in shell 1 around a core with 13 neighbours.
    // Construct through `PackingConsts::new` so the relation has one home.
    let pack = PackingConsts::new(k);

    // **`base_mass` must be dyadic and `contact_defect` need not be.** `Mass`
    // is fixed-point at SCALE = 1024; `units * base_mass` is never rounded, so
    // that factor has to land on the grid by itself, while the defect is
    // quantised explicitly below. Probed both ways: drawing `base_mass` in
    // steps of 0.3 fails exactness at N = 1, changing `contact_defect`'s
    // denominator from 512 to 500 changes nothing, and deleting the
    // quantisation fails at N = 2. An earlier comment called both
    // load-bearing, which would have sent a future reader to guard the wrong
    // half.
    let base_mass = 1.0 + 0.25 * rng.next_range(4) as f64;
    let contact_defect = (1 + rng.next_range(4)) as f64 / 512.0;
    let base_radius = 0.30 + 0.02 * rng.next_range(11) as f64;
    let eps = 0.8 + 0.05 * rng.next_range(9) as f64;
    let sigma = 0.02 + 0.01 * rng.next_range(12) as f64;
    let decay = 0.04 + 0.01 * rng.next_range(13) as f64;
    let n_elements = 60 + rng.next_range(61) as usize; // 60..=120
    let n_catalytic = 3 + rng.next_range(4) as u8;

    let mut naming_rng = Stream::new(seed, Domain::Naming, 0);
    let mut taken = Vec::new();
    let mut out: Vec<Element> = Vec::with_capacity(n_elements);

    for units in 1..=n_elements {
        // Which shell is filling, and how far into it.
        let (mut shell, mut filled) = (0_usize, 1_usize);
        while filled + packing::shell_size(k, shell + 1) <= units {
            shell += 1;
            filled += packing::shell_size(k, shell);
        }
        let outer = units - filled;
        let cap = packing::shell_size(k, shell + 1);
        let contacts = packing::contacts_upto(pack, units);

        let (symbol, name, _) = naming::mint(&mut naming_rng, &mut taken);

        // Mass, in sub-units of 1/1024, so the grid is a property of the
        // construction rather than a claim about it. `base_mass` is drawn on a
        // 0.25 step so `base_mass * 1024.0` is an exact integer; the defect is
        // quantised with `round_ties_even` to match `Mass`'s own convention —
        // ties are reachable (13 in this series, at N = 2, 3, 4, 7 among
        // others), and two conventions for one grid means a later "harmonise
        // the rounding" cleanup silently moves a golden.
        let base_sub = (base_mass * 1024.0).round_ties_even() as i64;
        let defect_sub = (contacts * contact_defect * 1024.0).round_ties_even() as i64;
        let mass = Mass::from_raw(units as i64 * base_sub - defect_sub);

        // Radius: set by which shell is occupied, and it cannot shrink.
        //
        // The enclosing radius of a superset of units is never smaller, and an
        // earlier `cbrt(units) * (1 + 0.55*fill)` collapsed it 35% at every
        // closure (N=74 -> 75: 2.587 -> 1.687 on one unit added). Since G1
        // measured binding at ~93% size and ~7% shape, that cliff would have
        // read as closed-shell elements binding a different partner set — an
        // arithmetic accident wearing the look of interesting periodicity.
        //
        // For a shelled cluster the radius *is* the shell index: each shell adds
        // one lattice spacing. So `shell + fill` is both monotone and the better
        // model, and it recovers `N^(1/3)` asymptotically since N ~ k*shell^3/3.
        // Over the range a table actually spans the fitted exponent is
        // 0.29..0.32 against a true cube root's 0.333 — the `1.0 +` offset drags
        // it down. Not a defect; an unqualified range claim would be.
        //
        // G3: real atomic radius falls across a period and jumps at a new one.
        // This rises monotonically and never jumps — structurally unlike,
        // without being geometrically wrong.
        let fill = outer as f64 / cap as f64;
        let radius = base_radius * (1.0 + shell as f64 + fill);

        // Valence: docking notches on the frontier. Zero at a closure because
        // there is no frontier there.
        //
        // `round_ties_even`, for the same reason as the mass defect: `notches`
        // lands on an exact 2.5 at k=8/outer=2 and k=8/outer=8, where ties-away
        // and ties-even differ by a whole bonding slot.
        let valence = packing::frontier_notches(cap, outer).round_ties_even() as u8;

        // Affinity: exposed fraction of the cluster. Continuous across every
        // closure by construction — at `outer == 0` it reduces to the
        // outermost complete shell, and at `outer == cap` to the same
        // expression for the next one.
        let shell_units = if shell == 0 {
            1
        } else {
            packing::shell_size(k, shell)
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
        // stretched over a larger radius than the one below and carries a
        // strain growing as `n^2`, which competes with the contact term's
        // saturation and produces a maximum.
        //
        // (The defence "no charge, no isospin, no pairing and no asymmetry
        // term" stood here and is deleted: it enumerates the terms that are
        // *absent* and says nothing about the one that is present, which is the
        // shape of the `libm` doctrine error. §7.2's prose forbids it.)
        //
        // **The exponent 2/3 is CHOSEN, not derived, and it is load-bearing.**
        // Per-*shell* strain proportional to `n^2`, summed and divided by `N`,
        // gives `1/k` — a constant, not `N^(2/3)`. Reaching 2/3 needs
        // per-*site* strain proportional to `n^2`, which is a different claim.
        // Accumulating shell-by-shell instead — as `contacts_upto` does, for
        // the identical reason — moves the peak in **48.1% of drawn
        // universes**, so the `19 of 24` closure census is a property of this
        // choice and not of the packing. Deferred to Task 20 deliberately:
        // precedence 2, not a G3 breach, remedy known and it requotes every
        // number. Do not describe the exponent as derived.
        let energy_per_unit =
            eps * contacts / units as f64 - sigma * det_math::cbrt((units * units) as f64);

        // Abundance: fusion builds heavy clusters from light ones, so each
        // extra unit costs a step and abundance falls geometrically. Sequential
        // addition with a constant per-step survival gives `r^N = exp(N ln r)`.
        //
        // **There is deliberately no `exp(beta * energy_per_unit)` tilt.** A
        // draft had one at beta = 0.6. Its stated mechanism compounds to
        // `E_total`, not `E_per_unit`; implementing the stated story makes the
        // heaviest element the most abundant and collapses p_c to 0.112. And
        // 0.6 was typed in where every neighbouring constant is drawn, which is
        // the same shape of defect as `peak`-as-an-argument.
        //
        // **This is the gel lever (§7.2) and Task 4 gates on nothing.** Measured
        // nominal p_c is 0.47..0.91 across 24 universes (median 0.74) — read it
        // from a probe run rather than from this comment. Do NOT justify that
        // with "burial makes realised functionality lower, so realised p_c is
        // higher" — §7.2 carried that argument and it is wrong: under
        // degree-independent thinning the beaker gels at the same bonds per atom
        // regardless, and only the coordinate was rescaled. Report nominal f_w
        // and bonds-per-node at gel as diagnostics; the gate is Task 20's
        // battery, on `p_ss / p_c`.
        let abundance = det_math::exp(-(units as f64) * decay);

        // Catalytic class: which band of outer-shell occupancy this element
        // exposes.
        //
        // **This is a rebinning of `group`, and the claim is narrowed to say
        // so.** A draft's comment here said Euler's twelve five-coordinate sites
        // mean a sparsely-filled shell exposes a different coordination mix from
        // a nearly-full one, "and that mix is what a mineral surface presents".
        // The code computed no such mix. Enumerated across all k, n_catalytic
        // and N, this takes exactly one value per `(outer, cap)` pair — zero
        // information beyond `group`. Two reviewers found it independently.
        //
        // So Task 4's headline is five §7.1 properties derived, plus
        // `abundance` — and this one honest rebinning, named as such.
        // Until it means something, nothing downstream may branch on it as
        // though it carried chemistry.
        //
        // `outer < cap` by loop construction and `cap >= 8`, so this lands in
        // `0..n_catalytic` with neither a `max(1)` nor a trailing modulo.
        let catalytic_class = ((outer * usize::from(n_catalytic)) / cap) as u8;

        out.push(Element {
            id: ElementId(out.len() as u8),
            symbol,
            name,
            units,
            period: shell as u8,
            group: outer as u8, // `outer < cap` and the table caps at 120
            mass,
            valence,
            affinity,
            radius: Span(radius),
            energy_per_unit,
            decay_rate: 0.0, // filled below, once the peak is known
            abundance,
            catalytic_class,
        });
    }

    // **The peak is read off here, after the series exists.** This ordering is
    // the mechanism: nothing above could have consulted it.
    //
    // Every closure is a *local* maximum by mechanism — a closed shell has no
    // frontier, so it makes every lateral contact available to it and maximises
    // contacts per unit, giving a sawtooth. That is structural: zero failures
    // across the full (k, eps, sigma) product.
    //
    // The **global** max is not. A series truncated mid-shell can peak at its
    // edge, and over the drawn table size it does.
    //
    // **Measured here, over `generate_elements` seeds 0..24:** on a closure in
    // 19 of 24, at the table edge in 2, mid-shell in 3; 13 distinct peak
    // positions; table sizes 61..120.
    //
    // The plan reports the same 19 of 24 and 13 distinct, but splits the
    // remaining five as 4 at the edge and 1 mid-shell. That is not a
    // contradiction and it is worth knowing why: those figures come from
    // `experiments/src/bin/fusion.rs`, which draws its constants from a
    // *different* stream, so its "24 seeds" are 24 different universes. Any
    // per-universe census quoted in this file has to be measured from this
    // function. The aggregate shape agrees; the individual tail does not, and
    // only the aggregate was ever a property of the packing.
    //
    // Both errors have been made here — an earlier version called this a
    // coincidence rate, and its correction called it structural without
    // qualification. Neither is right, which is why the census is reported
    // rather than asserted.
    let peak = out
        .iter()
        .fold((1_usize, f64::MIN), |(bn, be), e| {
            if e.energy_per_unit > be {
                (e.units, e.energy_per_unit)
            } else {
                (bn, be)
            }
        })
        .0;

    // Decay rate: distance below the peak, so the most tightly bound elements
    // persist and the extremes decay. One quantity, two consequences — which
    // is what the predecessor's two unlinked encodings could not give.
    let peak_energy = out.get(peak - 1).map_or(0.0, |e| e.energy_per_unit);
    for e in &mut out {
        // `f64::max` is disallowed — it returns either input on a tie and
        // measured `(+0.0).max(-0.0)` differs between aarch64 and x86-64.
        let scale = peak_energy.abs();
        let scale = if scale > f64::EPSILON {
            scale
        } else {
            f64::EPSILON
        };
        let deficit = (peak_energy - e.energy_per_unit) / scale;
        // **Cap at 1.0, and the 0.9 it replaces was never argued.** `decay_rate`
        // is a probability per world-year (§7.1), so 1.0 is the value the type
        // implies and needs no defence. 0.9 does need one, and had none —
        // measured over the 2916-universe grid it binds on 1.19% of elements
        // but on **at least one element in every universe**, and those elements
        // sit at median position N/n_elements = 0.011: it is the monomer, the
        // feedstock, the most abundant species. So a bare 0.9 systematically
        // slowed the decay of the one species whose persistence matters most,
        // by up to 26%, as a side effect of a guard nobody described. The raw
        // deficit does exceed 1 (max measured 1.218), so a cap is genuinely
        // needed — it is the value that was wrong, not the clamp.
        //
        // If a future measurement wants the monomer to persist longer, that is
        // a physics claim about feedstock and belongs in §9.4 with the
        // measurement attached, not in a magic number here.
        e.decay_rate = if deficit > 1.0 { 1.0 } else { deficit };
    }

    let shell = ShellPattern {
        k,
        closures: packing::closures(k, n_elements),
        peak,
    };
    (shell, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Named for brevity below; `generate_elements` is infallible because the
    /// mass is built in sub-units (see its doc), so there is nothing to unwrap.
    fn table(seed: u64) -> (ShellPattern, Vec<Element>) {
        generate_elements(seed)
    }

    /// The maximum valence of a table at an explicit `(k, n_elements)`.
    ///
    /// **This is the whole input space for that quantity.** `valence` is
    /// `frontier_notches(cap, outer).round_ties_even()`, and `cap`/`outer` are
    /// functions of `k` and `units` alone — none of `base_mass`,
    /// `contact_defect`, `decay`, `eps` or `sigma` enters. So enumerating
    /// `9 x 61` cells is exhaustive where sampling seeds is not.
    ///
    /// The shell walk below mirrors `generate_elements` exactly; keep them in
    /// step, or the enumeration stops measuring the thing it names. No
    /// `PackingConsts` is built because none is needed — `z` never enters a
    /// valence, and constructing one here would imply otherwise.
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "`frontier_notches` is non-negative and bounded well under 255 — \
                  `unmade_lateral_is_finite_and_non_negative_everywhere` pins the \
                  first and this test's own assertion pins the second"
    )]
    fn valence_ceiling_at(k: usize, n_elements: usize) -> u8 {
        (1..=n_elements)
            .map(|units| {
                let (mut shell, mut filled) = (0_usize, 1_usize);
                while filled + packing::shell_size(k, shell + 1) <= units {
                    shell += 1;
                    filled += packing::shell_size(k, shell);
                }
                let cap = packing::shell_size(k, shell + 1);
                packing::frontier_notches(cap, units - filled).round_ties_even() as u8
            })
            .max()
            .unwrap_or(0)
    }

    #[test]
    fn generation_is_deterministic() {
        // Purity, not determinism: `RandomState` is seeded once per process,
        // so two identically-built `HashMap`s iterate identically here and
        // would pass. Cross-platform reproducibility is the §13.6 golden
        // matrix's job and nothing in this file can stand in for it.
        assert_eq!(table(42), table(42));
    }

    /// Symbol and name must not be swapped on the way into `Element`.
    ///
    /// `let (symbol, name, _) = naming::mint(..)` feeds `Element { symbol,
    /// name, .. }` by field-init shorthand, so `let (name, symbol, _)` compiles
    /// and silently swaps them. `symbols_are_well_formed` tests `mint`
    /// directly and never sees it; nothing else checks either field at element
    /// level. Two `String`s in a tuple carry no type-level distinction between
    /// a 1–2 character symbol and a 4+ character name, so the assertion has to
    /// supply it.
    #[test]
    fn symbol_and_name_are_not_transposed() {
        for seed in 0..8 {
            for e in table(seed).1 {
                assert!(
                    (1..=2).contains(&e.symbol.len()),
                    "seed {seed}: symbol {:?} is not 1-2 chars",
                    e.symbol
                );
                assert!(
                    e.name.len() >= 4,
                    "seed {seed}: name {:?} is shorter than any minted name",
                    e.name
                );
            }
        }
    }

    #[test]
    fn different_seeds_give_different_tables() {
        assert_ne!(table(1), table(2));
    }

    /// **A cheap behavioural backstop, not the guarantee.**
    ///
    /// Zero-valence-at-closure is held by two independent factors, so the
    /// deletion test cannot tell a formula from an `if closed { 0 }`.
    /// Continuity was the proposed replacement and is **also insufficient**: it
    /// catches a declaration laid over a formula returning something *large*
    /// (`notches + 3` fails it) but not one writing 0 over a formula returning
    /// ~1, because both neighbours of a closure are pinned at valence 1 and the
    /// declaration sits inside the bound. Measured: a mutation putting every
    /// closure at valence 1 passes the entire suite.
    ///
    /// The guarantee is `xtask`'s source check; the output is guarded by
    /// `closed_shells_have_valence_zero`. This is the third line of defence,
    /// and the cheapest.
    #[expect(
        clippy::indexing_slicing,
        reason = "`windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn valence_is_continuous_across_closures() {
        for seed in 0..12 {
            let (_, els) = table(seed);
            for pair in els.windows(2) {
                let jump = i16::from(pair[1].valence) - i16::from(pair[0].valence);
                assert!(
                    jump.abs() <= 2,
                    "seed {seed} N={} -> {}: valence jumped {jump}",
                    pair[0].units,
                    pair[1].units
                );
            }
        }
    }

    /// The family §7.1 asks for, with a cause behind it. The element one unit
    /// past any closure is a closed core plus a single adhering unit, so it has
    /// exactly one docking notch, in every shell of every universe.
    ///
    /// **Its stated cause is not its cause**, and the correction is worth
    /// keeping: the valence is 1 by rounding `6/(6 - 12/cap)` in [1.14, 1.33],
    /// not because a lone site has "all six lateral contacts unmade" — the
    /// shell offers 4.5 to 5.9. Margin to 1.5 is comfortable, so the family
    /// holds.
    #[test]
    fn elements_one_past_a_closure_all_have_valence_one() {
        for seed in 0..12 {
            let (shell, els) = table(seed);
            for &c in &shell.closures {
                if let Some(e) = els.get(c) {
                    assert_eq!(
                        e.valence, 1,
                        "seed {seed}: N={} follows closure {c}",
                        e.units
                    );
                }
            }
        }
    }

    /// **Valence 8 is not reachable in a table this scheme generates, and
    /// saying so is the finding.**
    ///
    /// A per-*site* count caps at `z - 6`; summing over the frontier removes
    /// that ceiling in principle, and valence 9 does appear once a fourth shell
    /// opens past N = 200. But `n_elements` is drawn on 60..=120, and inside
    /// that range the ceiling is **4 to 6**.
    ///
    /// **Assert `4..=6`, not `4..=7`.** This assertion has been wrong three
    /// times: `5..=7` from a guess, which failed with a 4 in it; widened to
    /// `4..=7`, which was then loose enough to survive `FRONTIER_COEFF` being
    /// refitted from 7.30 to 6.90 — the refit moved the measured ceiling and
    /// the test did not notice. A range assertion that survives a change to the
    /// quantity it measures is not measuring it. Do not widen it back.
    ///
    /// A draft asserted `reached >= 12` of 24. The true rate is 273/549 = 49.7%,
    /// so the threshold sat on the population mean and the test would have
    /// failed on ~43% of seed sets — and its guidance sent the implementer to
    /// widen `k` when the controlling term is `n_elements`.
    ///
    /// **Enumerate the grid; do not sample seeds.** `valence` reads only `cap`
    /// and `outer`, both functions of `k` and `units` alone, so the ceiling is
    /// a pure function of `(k, n_elements)` — 9 x 61 = 549 cells, about a
    /// millisecond. The *third* draft of this assertion still sampled 48 seeds,
    /// and was still blind on the high side: at `FRONTIER_COEFF = 7.20`, inside
    /// the band `frontier_matches_exact_counts_on_a_real_sphere` admits, the
    /// ceiling reaches 7 in 7 of 549 cells and 48 draws miss all seven **54% of
    /// the time**. §7.1 puts valence at 0..6, so a 7 is a spec breach and this
    /// is the only test positioned to catch it. Raising 48 to 480 shrinks the
    /// probability instead of removing it and is not a fix.
    ///
    /// Assert the attained **set**, not containment — a refit collapsing every
    /// universe to 5 satisfies `.all(|v| (4..=6).contains(v))`.
    #[test]
    fn valence_ceiling_is_four_to_six_over_the_drawn_range() {
        let mut seen = std::collections::BTreeSet::new();
        for k in 6..=14 {
            for n in 60..=120 {
                seen.insert(valence_ceiling_at(k, n));
            }
        }
        let expected: std::collections::BTreeSet<u8> = [4, 5, 6].into_iter().collect();
        assert_eq!(
            seen, expected,
            "the attained valence ceiling set moved over the full drawn grid"
        );
    }

    /// §7.1: mass increases down the table. The fitted series this replaced
    /// stopped increasing at N = 124 — so the third closed shell weighed less
    /// than the element below it — and a 120-element table did not show it.
    #[expect(
        clippy::indexing_slicing,
        reason = "`windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn mass_increases_with_index() {
        for seed in 0..12 {
            let (_, els) = table(seed);
            for w in els.windows(2) {
                assert!(
                    w[1].mass > w[0].mass,
                    "seed {seed}: {} !> {}",
                    w[1].name,
                    w[0].name
                );
            }
        }
    }

    /// **`Mass::from_raw` validates nothing** — its own doc says so, and `impl
    /// Add` wraps in release on an out-of-range value. Since the construction is
    /// now integer, grid-exactness is a property of the type rather than a claim
    /// about the series; what still needs checking is that the drawn range
    /// cannot leave `MAX_MAGNITUDE`.
    ///
    /// Three drafts of this test were wrong in three ways. The first compared
    /// `Mass` to `Result<Mass, _>` and did not compile. The second round-tripped
    /// an already-quantised value — `to_f64` is an exact division by a power of
    /// two and `from_f64_quantised` an exact multiplication back, so it could
    /// not fail for anything the construction can produce. The third *described*
    /// the repair in this comment and left the round-trip in the body, which is
    /// how it reached review.
    ///
    /// **Assert the raw sub-unit count against the tight bound**, and the
    /// reason is *tightness* — nothing to do with overflow. A fourth draft
    /// justified it by overflow; the overflow needs `units >= 2^63/1792 ~ 5.2e15`
    /// against a loop bound of 120, unreachable by thirteen orders of magnitude.
    ///
    /// What the tight bound buys: a draw range that moves fails *here*, instead
    /// of silently spending 4993x of headroom. `MAX_MAGNITUDE` would not notice.
    #[test]
    fn every_mass_is_inside_the_range_from_raw_does_not_check() {
        // `const`, not `checked_mul(..).expect(..)`: const evaluation makes an
        // overflowing product a *compile* error (E0080), which is strictly
        // stronger than a runtime check, and avoids an `expect` that
        // `clippy::expect_used` denies — `clippy.toml` sets neither
        // `allow-expect-in-tests` nor `allow-unwrap-in-tests`, so the deny
        // reaches inside `#[cfg(test)]`.
        const MAX_UNITS: i64 = 120;
        const MAX_BASE_SUB: i64 = 1792;
        const CEILING: i64 = MAX_UNITS * MAX_BASE_SUB;
        for seed in 0..24 {
            let (_, els) = table(seed);
            for e in els {
                let raw = e.mass.raw();
                assert!(
                    (1..=CEILING).contains(&raw),
                    "seed {seed} N={}: raw {raw} outside 1..={CEILING} — a draw range moved",
                    e.units
                );
            }
        }
    }

    /// A closed shell has no frontier, so it has no unmade lateral contacts.
    /// This is the mechanism's headline output and what §7.2's "the most
    /// tightly bound element cannot bond" reading rests on.
    ///
    /// **It was in the probe and not carried across, and nothing else covers
    /// it.** Measured: a mutation breaking both zeroing factors puts every
    /// closure at valence 1 and the whole suite still passes — continuity
    /// trivially (1→1), the ceiling test, and `elements_one_past_a_closure`.
    ///
    /// It is *not* the test that distinguishes mechanism from declaration; two
    /// independent factors hold the zero, so a declaration can be laid over a
    /// broken formula and this still passes. That claim is a property of the
    /// source — see the `xtask` check.
    #[test]
    fn closed_shells_have_valence_zero() {
        for seed in 0..24 {
            let (sp, els) = table(seed);
            for n in packing::closures(sp.k, els.len()) {
                if let Some(e) = els.get(n - 1) {
                    assert_eq!(e.valence, 0, "seed {seed}: N={n} is a closure");
                }
            }
        }
    }

    /// Affinity is `2 * exposed_fraction - 1`, so it cannot leave `(-1, 1]`
    /// unless burial is miscounted. The predecessor reached exactly +1.000 at
    /// every closure by reading the capacity of the *empty* shell, making the
    /// stable elements the stickiest and stopping anything accumulating.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        clippy::indexing_slicing,
        reason = "`units` is at most 120; `windows(2)` yields slices of exactly length 2"
    )]
    #[test]
    fn affinity_is_in_range_and_continuous() {
        for seed in 0..12 {
            let (_, els) = table(seed);
            assert!(els.iter().all(|e| e.affinity > -1.0 && e.affinity <= 1.0));
            for w in els.windows(2) {
                let bound = 2.0 / w[0].units as f64;
                assert!(
                    (w[1].affinity - w[0].affinity).abs() <= bound,
                    "seed {seed} N={}: affinity jumped {}",
                    w[0].units,
                    (w[1].affinity - w[0].affinity).abs()
                );
            }
        }
    }

    /// The binding peak sits at a closure by *mechanism* — a closed shell has
    /// no frontier, so it maximises contacts per unit and the energy series is
    /// a sawtooth. If a change smooths the series across closures this fails,
    /// which is the point: the predecessor's minimum was the distance from a
    /// typed-in constant and looked identical from outside.
    #[test]
    fn energy_per_unit_is_locally_maximal_at_every_closure() {
        let (shell, els) = table(3);
        for &c in shell.closures.iter().filter(|&&c| c > 1) {
            let (lo, at, hi) = (els.get(c - 2), els.get(c - 1), els.get(c));
            if let (Some(lo), Some(at), Some(hi)) = (lo, at, hi) {
                assert!(
                    at.energy_per_unit > lo.energy_per_unit
                        && at.energy_per_unit > hi.energy_per_unit,
                    "closure {c} is not a local energy maximum"
                );
            }
        }
    }

    /// `peak` is an output. If it took only one value across universes it
    /// would merely have moved from an argument to a hard-coded location.
    #[test]
    fn the_binding_peak_moves_with_the_seed() {
        let peaks: std::collections::BTreeSet<usize> = (0..24).map(|s| table(s).0.peak).collect();
        assert!(
            peaks.len() >= 6,
            "peak took only {} values: {peaks:?}",
            peaks.len()
        );
    }

    /// Different universes must get different table *shapes*, not the same
    /// shape with different numbers in it.
    #[test]
    fn shell_patterns_vary_between_universes() {
        let patterns: std::collections::BTreeSet<Vec<usize>> =
            (0..30).map(|s| table(s).0.closures).collect();
        assert!(
            patterns.len() > 5,
            "shell patterns barely vary: {}",
            patterns.len()
        );
    }

    #[test]
    fn table_size_is_workable() {
        for seed in 0..30 {
            let n = table(seed).1.len();
            // The draw is 60..=120, so the lower half of a `40..=120`
            // assertion was a guard that could not fire — the defect class this
            // task is partly about, asserted 200 lines from a comment rejecting
            // it.
            assert!((60..=120).contains(&n), "seed {seed} produced {n} elements");
        }
    }
}
