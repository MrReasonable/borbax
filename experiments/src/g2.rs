//! The G2 measurement: does signature space have locality?
//!
//! # The question
//!
//! §8.4's evolvability argument assumes a one-atom edit moves a molecule a
//! *small* distance in signature space — that chemistry is explorable by
//! small steps rather than by teleportation. If that is false, the argument
//! has nothing under it, and that is a falsified design premise rather than a
//! bug to fix.
//!
//! # The protocol
//!
//! Each trial draws a parent molecule, applies one single-atom edit to get a
//! mutant, and draws an unrelated molecule independently. Every descriptor
//! then answers one question: **is the mutant nearer to its parent than a
//! stranger is?** Aggregated over trials that gives a concordance index —
//! `1.0` for perfect locality, `0.5` for none, and the earlier measurement's
//! "32% of mutants land further than an unrelated molecule" corresponds to
//! `0.68`.
//!
//! Reported in two regimes. Unrelated molecules are normally drawn at their
//! own size, so a descriptor can succeed just by tracking how big something
//! is; the **size-matched** regime holds atom count equal, leaving only shape.
//! Given that binding was separately measured as ~93% size and ~7% shape, the
//! gap between the two regimes is the interesting quantity, not either alone.
//!
//! # Why the controls are not optional
//!
//! A harness that reports a number it cannot calibrate is measuring itself.
//! Six earlier harnesses produced figures nobody can now check, and this one
//! refuses to report a G2 result at all unless a descriptor known to have
//! locality shows it and a descriptor known to have none does not. That gate
//! has its own tests, because a gate that cannot fail is the same mistake one
//! level up — the same principle as exit criterion 9.
#![expect(
    clippy::disallowed_methods,
    reason = "§13.4: concordance values are ratios of non-negative counts — never NaN, \
                  never -0.0"
)]

use crate::embed::{Point, embed};
use crate::geodesic::Geodesic;
use crate::geodesic::{apply_mat, is_identity, rotation_matrices};
use crate::molecule::Molecule;
use crate::rng::Stream;
use crate::signature::{D, Descriptor, d_group, d_raw, profile, signature};

/// The positive control must clear this, or the harness cannot see locality.
///
/// Not tuned to the observed value. A one-atom edit changes composition by at
/// most two by construction, while two independently drawn molecules differ by
/// far more, so any working instrument clears this comfortably. It is set well
/// below what a working run produces precisely so that passing it is weak
/// evidence of health and failing it is strong evidence of breakage.
pub const POSITIVE_CONTROL_MIN: f64 = 0.75;

/// How far the negative control may sit from `0.5` in either direction.
///
/// Two independently avalanched hashes make the comparison a coin flip, so the
/// standard error at 2000 trials is about `0.011`; `0.05` is roughly four of
/// those. Two-sided deliberately — a reading well *below* `0.5` would mean the
/// statistic is inverted somewhere, which would flip every conclusion while
/// looking like a strong result.
pub const NEGATIVE_CONTROL_TOLERANCE: f64 = 0.05;

/// The two controls must be this far apart for any reading between them to
/// carry information.
pub const MIN_CONTROL_SEPARATION: f64 = 0.25;

/// What a run measures over.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Config {
    /// Independent parent/mutant/stranger triples.
    pub trials: usize,
    /// Smallest parent molecule.
    pub min_atoms: usize,
    /// Largest parent molecule.
    pub max_atoms: usize,
    /// Root seed. Each trial derives its own stream from this and its index.
    pub seed: u64,
    /// Draw the unrelated molecule at the parent's exact atom count, removing
    /// size as a signal and leaving only shape.
    pub size_matched: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            trials: 2000,
            min_atoms: 8,
            max_atoms: 20,
            seed: 0x0B0B_BA05,
            size_matched: false,
        }
    }
}

/// One descriptor's result.
#[derive(Debug, Clone, PartialEq)]
pub struct Stats {
    /// Which descriptor.
    pub descriptor: Descriptor,
    /// Median parent-to-mutant distance.
    pub median_mutant: f64,
    /// Median parent-to-stranger distance.
    pub median_unrelated: f64,
    /// Fraction of trials where the mutant landed *further* than the stranger.
    /// This is the figure the earlier measurement reported as 32%.
    pub p_mutant_further: f64,
    /// Fraction of trials where the two distances were exactly equal.
    pub ties: f64,
    /// `P(mutant nearer) + 0.5 * P(tie)`. `1.0` is perfect locality, `0.5` none.
    pub concordance: f64,
}

/// Minimum distance `D_raw` must report between a molecule and a bodily
/// rotated copy of itself. Measured on a working pipeline it is ~1-4 at
/// D = 42; a dead one gives ~1e-15.
///
/// **This is an absolute distance in signature units, not a scale-free ratio,
/// and that is worth knowing before anything rescales the radii.** Signature
/// extents are O(1-10) at D = 42 under the current `ELEMENT_RADII`. Measured
/// over 400 random 12-atom molecules x 59 non-identity rotations, the smallest
/// `D_raw` observed is **1.40** (1.30 on an independent draw). The nearest
/// genuinely degenerate case is a *2-atom* molecule, whose floor is **7.8e-4**
/// and which this correctly rejects. So 0.1 sits ~14x below the working floor
/// and ~128x above the degenerate one — a wide gap, not a tuned value.
///
/// The exposure: that 13x margin is proportional to molecule size and to radius
/// scale. If a derived radius series replaces the ladder at a smaller scale, or
/// the fixture shrinks, the margin shrinks with it and this constant has to be
/// re-measured rather than assumed. [`ROTATION_MUST_NOT_MOVE_GROUP`] has no such
/// exposure — its margin is 1.65e8x and it is comparing against zero.
const ROTATION_MUST_MOVE_RAW: f64 = 0.1;

/// Maximum residual `D_group` may leave after recovering a rotation that is
/// itself one of the 60. Measured at ~1e-15 on a working pipeline.
const ROTATION_MUST_NOT_MOVE_GROUP: f64 = 1e-6;

/// Minimum `D_group` between a 12-path and a 12-star — same composition,
/// different topology. Measured on a working pipeline: 14.42 (7.76 at n = 8,
/// 21.53 at n = 16, so it scales with size rather than sitting near a cliff).
/// An `embed` that ignores the molecule gives exactly 0.0, so 1.0 sits 14x
/// below the working value and infinitely above the failure.
const TOPOLOGY_MUST_CHANGE_SHAPE: f64 = 1.0;

/// Atom count of the path/star fixtures the check above thresholds on.
const TOPOLOGY_FIXTURE_ATOMS: usize = 12;

/// Whether a measured quantity clears a floor it is required to exceed.
///
/// **Two separate guards, and it is worth knowing which does what**, because
/// the obvious reading credits the wrong one.
///
/// The first is the *shape* of the test. The gate used to read
/// `if got <= min { Err }`, and `NaN <= 0.1` is `false`, so a `NaN` was waved
/// through as healthy — the identical defect
/// [`ControlFailure::ControlNotFinite`] exists to catch one level up, and one
/// that `the_control_gate_rejects_controls_that_are_not_numbers` already
/// describes in its own doc comment: "every ordering comparison against `NaN`
/// is false — so a gate written as `positive < MIN` would wave it through as
/// healthy". Stating the *positive* condition and negating it at the call site
/// fixes that on its own, because `NaN > 0.1` is false too.
///
/// The `is_finite` call is the second guard and catches something else:
/// `+inf`, which genuinely does exceed any floor and would pass the rewritten
/// test on its own terms rather than by the `NaN` rule. Measured — deleting
/// `is_finite` from this function fails
/// `a_non_finite_measurement_clears_no_threshold` on the `inf` case and not on
/// the `NaN` one.
///
/// Neither is reachable today, and that was checked rather than assumed —
/// `normalise` guards `n > 0.0`, `guttman` skips non-positive distances, and
/// `signature` returns `NEG_INFINITY` only for empty coordinates, which the
/// gate's 12-atom fixtures cannot produce. It is here because a gate that
/// cannot fail is the mistake this module exists to avoid making.
fn clears_floor(got: f64, min: f64) -> bool {
    got.is_finite() && got > min
}

/// Whether a measured quantity stays below a ceiling it must not reach.
///
/// The mirror of [`clears_floor`], guarded twice for the same two reasons:
/// `NaN >= 1e-6` is false, so the old `if got >= max { Err }` waved it through;
/// and `-inf` is under every ceiling on its own terms.
fn under_ceiling(got: f64, max: f64) -> bool {
    got.is_finite() && got < max
}

/// Does the shape pipeline actually compute shape?
///
/// **The two controls this joins do not touch the geometry, and that was
/// measured rather than suspected.** `PositiveControl` reads only
/// `Profile::composition` and `NegativeControl` only `Profile::form_hash`;
/// neither goes through `embed`, `Geodesic`, `signature` or any descriptor.
/// Replacing `embed` with one returning every atom at the origin — no geometry
/// computed at all — left this harness printing `controls: PASS — the
/// descriptor rows may be read`, exiting 0, and passing its own calibration
/// test, while all four descriptor rows read 0.52 and the tie fraction was
/// 0.936. That is this module's own thesis failing one level up: a harness
/// that reports a number it cannot calibrate is measuring itself.
///
/// The check is a property, not a threshold on the reported rows. Rotate a
/// molecule bodily by one of the 60 group elements, and:
///
/// - `D_raw` compares direction-by-direction with no alignment, so it **must**
///   move. A dead pipeline reports ~0 and fails here.
/// - `D_group` searches all 60, one of which is exactly the rotation applied,
///   so it **must not** move. A broken permutation table or a signature that is
///   not a pure function of direction fails here.
///
/// **A dropped `ANTI` cannot fail either of those, and an earlier version of
/// this comment claimed it could.** `d_group` reads `g.perms` and never touches
/// `g.anti` — correctly, because `ANTI` exists for the binding kernel, where two
/// bodies in contact touch along *opposite* directions, and this compares one
/// shape against itself in two poses. That coverage has to come from Task 8's
/// kernel tests; advertising it here is worse than the gap, because it is the
/// sentence that stops someone writing the real test for the failure CLAUDE.md
/// names as the most dangerous in the project.
///
/// The rotation pair alone is also not enough, and that was measured rather
/// than reasoned. Replacing `embed` with one returning a **fixed configuration
/// keyed only on `mol.len()`** — no graph, no SMACOF, nothing from the molecule
/// — still moves `D_raw` and still recovers under `D_group`, so both checks
/// above pass while the harness reports a stable, plausible, entirely spurious
/// `D_frame` concordance of 0.88. So the third check compares two molecules of
/// **identical composition and different topology**: a 12-path against a
/// 12-star. It is the composition control's shape counterpart — `PositiveControl`
/// proves composition reaches the descriptor, this proves topology does.
///
/// Two molecules, one rotation: this is a liveness check, not a correctness
/// proof. The proofs are `rotating_the_molecule_permutes_the_signature` and
/// `the_permutation_table_agrees_with_the_rotation_matrices`, which are unit
/// tests and do not gate the binary.
fn check_geometry_is_live() -> Result<(), ControlFailure> {
    let g = Geodesic::<D>::build().map_err(|_| ControlFailure::ControlMissing("geodesic"))?;
    let mats = rotation_matrices().map_err(|_| ControlFailure::ControlMissing("rotations"))?;
    // *Which* rotation is immaterial — any non-identity element moves the
    // molecule. This was `mats[7]`, which is a fact about `rotation_matrices()`'
    // construction order standing in for a property nothing checked: reorder the
    // group so index 7 becomes the identity, and the gate fails on a healthy
    // pipeline with a diagnostic blaming `embed`. Naming the index in a constant
    // documents the requirement without enforcing it, and this repository has
    // four recorded defects of exactly that shape. So ask for the property.
    let mat = mats
        .iter()
        .find(|m| !is_identity(m))
        .ok_or(ControlFailure::ControlMissing("a non-identity rotation"))?;

    let mut rng = Stream::new(0x006E_01E7);
    let m = Molecule::random_tree(&mut rng, 12);
    let coords = embed(&m);
    let turned: Vec<Point> = coords.iter().map(|&p| apply_mat(mat, p)).collect();

    let a = signature(&g, &coords, &m.elements);
    let b = signature(&g, &turned, &m.elements);

    let moved = d_raw(&a, &b);
    if !clears_floor(moved, ROTATION_MUST_MOVE_RAW) {
        return Err(ControlFailure::GeometryIsDead {
            got: moved,
            min: ROTATION_MUST_MOVE_RAW,
        });
    }
    let recovered = d_group(&g, &a, &b);
    if !under_ceiling(recovered, ROTATION_MUST_NOT_MOVE_GROUP) {
        return Err(ControlFailure::RotationSearchIsBroken {
            got: recovered,
            max: ROTATION_MUST_NOT_MOVE_GROUP,
        });
    }

    // Same composition, different topology. Every atom is element 0 in both, so
    // nothing but the graph can separate them — which is exactly what an
    // `embed` keyed on `mol.len()` alone cannot do.
    let atoms = TOPOLOGY_FIXTURE_ATOMS;
    let path = Molecule::from_parts(vec![0; atoms], (1..atoms).map(|i| (i - 1, i)).collect());
    let star = Molecule::from_parts(vec![0; atoms], (1..atoms).map(|i| (0, i)).collect());
    let sig_path = signature(&g, &embed(&path), &path.elements);
    let sig_star = signature(&g, &embed(&star), &star.elements);

    let separation = d_group(&g, &sig_path, &sig_star);
    if !clears_floor(separation, TOPOLOGY_MUST_CHANGE_SHAPE) {
        return Err(ControlFailure::GeometryIgnoresMolecule {
            got: separation,
            min: TOPOLOGY_MUST_CHANGE_SHAPE,
            n: atoms,
        });
    }
    Ok(())
}

/// Why a run's numbers may not be believed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ControlFailure {
    /// The shape pipeline produced no shape: a molecule and a bodily rotated
    /// copy were indistinguishable to an unaligned descriptor.
    #[error(
        "the geometry is dead: a rotated copy moved D_raw by only {got:.3e}, expected > {min} \
         — embed() or signature() is not computing shape, so every descriptor row is noise"
    )]
    GeometryIsDead {
        /// Measured `D_raw` between a molecule and its rotated copy.
        got: f64,
        /// The floor it had to clear.
        min: f64,
    },
    /// The 60-rotation search failed to recover a rotation that is one of the
    /// 60 — so `perms` or the signature's direction-purity is broken.
    ///
    /// **Not `anti`.** `d_group` never reads it; see `check_geometry_is_live`.
    #[error(
        "the rotation search is broken: D_group left {got:.3e} between a molecule and a copy \
         rotated by a group element, expected < {max} — check perms and that the signature \
         is a pure function of direction (not anti; d_group does not read it)"
    )]
    RotationSearchIsBroken {
        /// Measured `D_group` residual.
        got: f64,
        /// The ceiling it had to stay under.
        max: f64,
    },
    /// The pipeline computes *a* shape, but not one that depends on the
    /// molecule: two graphs of identical composition came out the same.
    #[error(
        "the geometry ignores the molecule: a {n}-path and a {n}-star of identical composition \
         differ by only {got:.3e} under D_group, expected > {min} — embed() is producing a \
         configuration that does not depend on the graph, so every descriptor row is an \
         artefact of atom count alone"
    )]
    GeometryIgnoresMolecule {
        /// Measured `D_group` between the path and the star.
        got: f64,
        /// The floor it had to clear.
        min: f64,
        /// Atom count of both fixtures.
        n: usize,
    },
    /// A control was not measured at all.
    #[error("the {0} was not measured — the run cannot be calibrated")]
    ControlMissing(&'static str),
    /// A control produced a non-finite value, so no comparison against it
    /// means anything.
    ///
    /// Checked explicitly rather than left to the comparisons below, because
    /// every ordering comparison against `NaN` is false — so `positive < MIN`
    /// would *pass* a control that had not been measured. A gate that lets
    /// `NaN` through is the failure this whole module exists to prevent,
    /// arriving through the back door.
    #[error("the {which} is {got}, not a usable number — the run cannot be calibrated")]
    ControlNotFinite {
        /// Which control.
        which: &'static str,
        /// The value it produced.
        got: f64,
    },
    /// The positive control failed to show locality that is there by
    /// construction, so the harness cannot detect locality at all.
    #[error(
        "positive control concordance {got:.4} is below {min:.2}: the harness cannot see \
         locality that is present by construction, so no other number here means anything"
    )]
    PositiveControlBlind {
        /// Measured concordance.
        got: f64,
        /// Required minimum.
        min: f64,
    },
    /// The negative control showed structure it cannot have, so the statistic
    /// is manufacturing signal.
    #[error(
        "negative control concordance {got:.4} is not within {tol:.2} of 0.5: an avalanched \
         hash has no locality to find, so the statistic is manufacturing it"
    )]
    NegativeControlLeaks {
        /// Measured concordance.
        got: f64,
        /// Allowed deviation from 0.5.
        tol: f64,
    },
    /// The controls are too close together for a reading between them to mean
    /// anything.
    #[error("controls are only {got:.4} apart, under {min:.2}: the scale cannot resolve a result")]
    ControlsTooClose {
        /// Measured separation.
        got: f64,
        /// Required minimum.
        min: f64,
    },
}

/// Everything one run produced.
#[derive(Debug, Clone)]
pub struct Report {
    /// What was measured.
    pub config: Config,
    /// One entry per descriptor, in [`Descriptor::ALL`] order.
    pub stats: Vec<Stats>,
}

/// The middle value, averaging the two middle elements of an even sample.
///
/// Sorted with `total_cmp`, which is a total order with no tolerance, so the
/// result does not depend on how the input happened to be ordered.
#[must_use]
pub fn median(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NAN;
    }
    let mut v = xs.to_vec();
    v.sort_by(f64::total_cmp);
    let mid = v.len() / 2;
    if v.len() % 2 == 1 {
        v.get(mid).copied().unwrap_or(f64::NAN)
    } else {
        match (v.get(mid.wrapping_sub(1)), v.get(mid)) {
            (Some(&lo), Some(&hi)) => f64::midpoint(lo, hi),
            _ => f64::NAN,
        }
    }
}

/// `P(mutant nearer) + 0.5 * P(tie)` over `(mutant, unrelated)` pairs.
///
/// Ties count as half because a discrete descriptor genuinely produces them —
/// composition differs by the same amount for a mutant and a stranger more
/// often than a continuous distance ever would — and scoring them as failures
/// would penalise the positive control for being discrete.
#[must_use]
#[allow(
    clippy::float_cmp,
    reason = "exact equality is the definition of a tie. The descriptors that actually produce \
              ties are discrete — a composition L1 or a popcount — so their equal values are \
              bit-identical. An epsilon would reclassify near-misses of the continuous \
              descriptors as ties, which would change the statistic rather than stabilise it"
)]
pub fn concordance(pairs: &[(f64, f64)]) -> f64 {
    if pairs.is_empty() {
        return f64::NAN;
    }
    let mut score = 0.0;
    for &(mutant, unrelated) in pairs {
        if mutant < unrelated {
            score += 1.0;
        } else if mutant == unrelated {
            score += 0.5;
        }
    }
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "a trial count, far below 2^53"
    )]
    let n = pairs.len() as f64;
    score / n
}

/// Fraction of pairs where the mutant landed strictly further away.
#[must_use]
pub fn p_mutant_further(pairs: &[(f64, f64)]) -> f64 {
    fraction(pairs, |m, u| m > u)
}

/// Fraction of pairs where the two distances were exactly equal.
#[must_use]
#[allow(
    clippy::float_cmp,
    reason = "see `concordance`: exact equality is the definition being measured"
)]
pub fn tie_fraction(pairs: &[(f64, f64)]) -> f64 {
    fraction(pairs, |m, u| m == u)
}

fn fraction(pairs: &[(f64, f64)], predicate: impl Fn(f64, f64) -> bool) -> f64 {
    if pairs.is_empty() {
        return f64::NAN;
    }
    let hits = pairs.iter().filter(|&&(m, u)| predicate(m, u)).count();
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "counts of trials, far below 2^53"
    )]
    let out = hits as f64 / pairs.len() as f64;
    out
}

/// Run the measurement.
#[must_use]
pub fn run(cfg: &Config) -> Report {
    let Ok(geo) = Geodesic::<D>::build() else {
        return Report {
            config: cfg.clone(),
            stats: Vec::new(),
        };
    };

    // One bucket of (mutant, unrelated) distance pairs per descriptor.
    let mut pairs: Vec<Vec<(f64, f64)>> =
        vec![Vec::with_capacity(cfg.trials); Descriptor::ALL.len()];
    let span = cfg.max_atoms.saturating_sub(cfg.min_atoms) + 1;

    for trial in 0..cfg.trials {
        // One stream per trial index, so a trial's molecules depend on its own
        // index and nothing else. Trials could be reordered or run in any
        // grouping and the numbers would not move.
        let mut rng = Stream::new(cfg.seed.wrapping_add(u64::try_from(trial).unwrap_or(0)));

        let n = cfg.min_atoms + rng.index(span);
        let parent = Molecule::random_tree(&mut rng, n);
        let mutant = parent.mutate(&mut rng);
        let stranger_size = if cfg.size_matched {
            n
        } else {
            cfg.min_atoms + rng.index(span)
        };
        let stranger = Molecule::random_tree(&mut rng, stranger_size);

        let p = profile(&geo, &parent);
        let m = profile(&geo, &mutant);
        let s = profile(&geo, &stranger);

        for (slot, descriptor) in pairs.iter_mut().zip(Descriptor::ALL) {
            slot.push((
                descriptor.distance(&geo, &p, &m),
                descriptor.distance(&geo, &p, &s),
            ));
        }
    }

    let stats = Descriptor::ALL
        .iter()
        .zip(pairs.iter())
        .map(|(&descriptor, bucket)| {
            let mutants: Vec<f64> = bucket.iter().map(|&(m, _)| m).collect();
            let strangers: Vec<f64> = bucket.iter().map(|&(_, u)| u).collect();
            Stats {
                descriptor,
                median_mutant: median(&mutants),
                median_unrelated: median(&strangers),
                p_mutant_further: p_mutant_further(bucket),
                ties: tie_fraction(bucket),
                concordance: concordance(bucket),
            }
        })
        .collect();

    Report {
        config: cfg.clone(),
        stats,
    }
}

impl Report {
    /// The result for one descriptor, if it was measured.
    #[must_use]
    pub fn get(&self, d: Descriptor) -> Option<&Stats> {
        self.stats.iter().find(|s| s.descriptor == d)
    }

    /// Whether this run's numbers may be reported at all.
    ///
    /// # Errors
    ///
    /// Any [`ControlFailure`]. A failure here is not a G2 result — it is a
    /// statement that the instrument is broken and the run says nothing about
    /// signature space either way.
    pub fn check_controls(&self) -> Result<(), ControlFailure> {
        // First, because without it every check below can pass on a pipeline
        // that computed no geometry at all.
        check_geometry_is_live()?;

        let positive = self
            .get(Descriptor::PositiveControl)
            .ok_or(ControlFailure::ControlMissing("positive control"))?
            .concordance;
        let negative = self
            .get(Descriptor::NegativeControl)
            .ok_or(ControlFailure::ControlMissing("negative control"))?
            .concordance;

        if !positive.is_finite() {
            return Err(ControlFailure::ControlNotFinite {
                which: "positive control",
                got: positive,
            });
        }
        if !negative.is_finite() {
            return Err(ControlFailure::ControlNotFinite {
                which: "negative control",
                got: negative,
            });
        }
        if positive < POSITIVE_CONTROL_MIN {
            return Err(ControlFailure::PositiveControlBlind {
                got: positive,
                min: POSITIVE_CONTROL_MIN,
            });
        }
        if (negative - 0.5).abs() > NEGATIVE_CONTROL_TOLERANCE {
            return Err(ControlFailure::NegativeControlLeaks {
                got: negative,
                tol: NEGATIVE_CONTROL_TOLERANCE,
            });
        }
        let separation = positive - negative;
        if separation < MIN_CONTROL_SEPARATION {
            return Err(ControlFailure::ControlsTooClose {
                got: separation,
                min: MIN_CONTROL_SEPARATION,
            });
        }
        Ok(())
    }

    /// The report as a table, with the control verdict at the foot.
    #[must_use]
    pub fn render(&self) -> String {
        use std::fmt::Write as _;
        let regime = if self.config.size_matched {
            "size-matched (shape only)"
        } else {
            "free size"
        };
        let mut out = String::new();
        let _ = writeln!(
            out,
            "G2 — locality of signature space\n\
             {} trials, {}-{} atoms, seed {:#x}, regime: {regime}\n",
            self.config.trials, self.config.min_atoms, self.config.max_atoms, self.config.seed
        );
        let _ = writeln!(
            out,
            "{:<24} {:>10} {:>12} {:>10} {:>8} {:>12}",
            "descriptor", "med(mut)", "med(unrel)", "P(further)", "ties", "concordance"
        );
        let _ = writeln!(out, "{}", "-".repeat(80));
        for s in &self.stats {
            let _ = writeln!(
                out,
                "{:<24} {:>10.4} {:>12.4} {:>10.3} {:>8.3} {:>12.4}",
                s.descriptor.name(),
                s.median_mutant,
                s.median_unrelated,
                s.p_mutant_further,
                s.ties,
                s.concordance
            );
        }
        let _ = writeln!(out);
        match self.check_controls() {
            Ok(()) => {
                let _ = writeln!(out, "controls: PASS — the descriptor rows may be read.");
            }
            Err(e) => {
                let _ = writeln!(
                    out,
                    "controls: FAIL — {e}\n\
                     The descriptor rows above are NOT a G2 result and must not be quoted."
                );
            }
        }
        out
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::float_cmp,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds; the float comparisons \
              are against values this module produced exactly, not accumulated"
)]
mod tests {
    use super::*;

    fn stats_with(d: Descriptor, concordance: f64) -> Stats {
        Stats {
            descriptor: d,
            median_mutant: 1.0,
            median_unrelated: 2.0,
            p_mutant_further: 1.0 - concordance,
            ties: 0.0,
            concordance,
        }
    }

    /// A report whose controls behave: the positive one sees locality, the
    /// negative one does not.
    fn healthy() -> Report {
        Report {
            config: Config::default(),
            stats: vec![
                stats_with(Descriptor::Frame, 0.7),
                stats_with(Descriptor::PositiveControl, 0.95),
                stats_with(Descriptor::NegativeControl, 0.50),
            ],
        }
    }

    #[test]
    fn median_handles_odd_and_even_samples() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        assert_eq!(median(&[7.0]), 7.0);
        assert!(
            median(&[]).is_nan(),
            "the median of nothing is not a number"
        );
    }

    #[test]
    fn concordance_is_one_when_every_mutant_is_nearer() {
        let pairs = [(1.0, 5.0), (0.5, 9.0), (2.0, 2.5)];
        assert_eq!(concordance(&pairs), 1.0);
        assert_eq!(p_mutant_further(&pairs), 0.0);
    }

    #[test]
    fn concordance_is_zero_when_every_mutant_is_further() {
        let pairs = [(5.0, 1.0), (9.0, 0.5)];
        assert_eq!(concordance(&pairs), 0.0);
        assert_eq!(p_mutant_further(&pairs), 1.0);
    }

    /// A descriptor carrying no information puts the mutant nearer half the
    /// time. This is the value the negative control has to land on, so the
    /// statistic has to produce it exactly on a balanced sample.
    #[test]
    fn concordance_is_a_half_on_a_balanced_sample() {
        let pairs = [(1.0, 2.0), (2.0, 1.0), (3.0, 4.0), (4.0, 3.0)];
        assert_eq!(concordance(&pairs), 0.5);
    }

    /// Ties are half a point. Without this the positive control would be
    /// penalised for the discreteness of a count — a one-atom edit and an
    /// unrelated molecule genuinely can differ by the same amount.
    #[test]
    fn ties_count_as_half() {
        let pairs = [(1.0, 1.0), (1.0, 1.0)];
        assert_eq!(concordance(&pairs), 0.5);
        assert_eq!(tie_fraction(&pairs), 1.0);
        assert_eq!(concordance(&[(1.0, 1.0), (1.0, 2.0)]), 0.75);
    }

    #[test]
    fn the_control_gate_accepts_a_working_instrument() {
        assert!(healthy().check_controls().is_ok());
    }

    /// The property this whole file exists to have. An instrument that cannot
    /// report failure is not an instrument, and the previous six harnesses
    /// produced numbers nobody can now check precisely because nothing here
    /// was ever written down.
    #[test]
    fn the_control_gate_rejects_a_blind_positive_control() {
        let mut r = healthy();
        r.stats[1] = stats_with(Descriptor::PositiveControl, 0.55);
        assert!(matches!(
            r.check_controls(),
            Err(ControlFailure::PositiveControlBlind { .. })
        ));
    }

    #[test]
    fn the_control_gate_rejects_a_negative_control_that_shows_locality() {
        let mut r = healthy();
        r.stats[2] = stats_with(Descriptor::NegativeControl, 0.72);
        assert!(matches!(
            r.check_controls(),
            Err(ControlFailure::NegativeControlLeaks { .. })
        ));
    }

    /// A negative control can fail in the other direction too — reading well
    /// *below* a half means the statistic is inverted somewhere, which would
    /// silently flip every conclusion.
    #[test]
    fn the_control_gate_rejects_an_inverted_negative_control() {
        let mut r = healthy();
        r.stats[2] = stats_with(Descriptor::NegativeControl, 0.28);
        assert!(matches!(
            r.check_controls(),
            Err(ControlFailure::NegativeControlLeaks { .. })
        ));
    }

    /// A zero-trial run makes every statistic `NaN`, and every ordering
    /// comparison against `NaN` is false — so a gate written as
    /// `positive < MIN` would wave it through as healthy. The instrument must
    /// refuse a run it never actually made.
    #[test]
    fn the_control_gate_rejects_controls_that_are_not_numbers() {
        let empty = run(&Config {
            trials: 0,
            ..Config::default()
        });
        assert!(matches!(
            empty.check_controls(),
            Err(ControlFailure::ControlNotFinite { .. })
        ));

        let mut r = healthy();
        r.stats[2] = stats_with(Descriptor::NegativeControl, f64::NAN);
        assert!(matches!(
            r.check_controls(),
            Err(ControlFailure::ControlNotFinite { .. })
        ));
    }

    #[test]
    fn the_control_gate_rejects_a_report_missing_a_control() {
        let r = Report {
            config: Config::default(),
            stats: vec![stats_with(Descriptor::Frame, 0.7)],
        };
        assert!(matches!(
            r.check_controls(),
            Err(ControlFailure::ControlMissing(_))
        ));
    }

    #[test]
    fn the_report_is_deterministic() {
        let cfg = Config {
            trials: 40,
            ..Config::default()
        };
        let a = run(&cfg);
        let b = run(&cfg);
        for (x, y) in a.stats.iter().zip(b.stats.iter()) {
            assert_eq!(
                x.concordance,
                y.concordance,
                "{} moved",
                x.descriptor.name()
            );
            assert_eq!(x.median_mutant, y.median_mutant);
            assert_eq!(x.median_unrelated, y.median_unrelated);
        }
    }

    /// **The gate must fail when the geometry is broken, and before this it
    /// did not.** Measured: replacing `embed` with one returning every atom at
    /// the origin left the harness printing "controls: PASS — the descriptor
    /// rows may be read", exiting 0, and passing this very module's
    /// calibration test, because neither control touches `embed`. All four
    /// descriptor rows read 0.52 with a 0.936 tie fraction and nothing
    /// objected.
    ///
    /// This is the liveness property stated directly, so a future refactor of
    /// `check_geometry_is_live` has something to fail against.
    #[test]
    fn a_rotated_copy_moves_the_unaligned_descriptor_and_not_the_aligned_one() {
        let g = Geodesic::<D>::build().unwrap();
        let mats = rotation_matrices().unwrap();
        let mut rng = Stream::new(0x006E_01E7);
        let m = Molecule::random_tree(&mut rng, 12);
        let coords = embed(&m);
        let a = signature(&g, &coords, &m.elements);

        // Nothing else checks this at the matrix level: `identity_is_present`
        // and `rotations_are_all_distinct` both test `perms`, and the matrices
        // inherit it only transitively through
        // `the_permutation_table_agrees_with_the_rotation_matrices`.
        assert_eq!(
            mats.iter().filter(|m| is_identity(m)).count(),
            1,
            "the group must contain the identity exactly once"
        );

        // Every group element, not just the one the gate samples.
        for (r, mat) in mats.iter().enumerate() {
            let turned: Vec<Point> = coords.iter().map(|&p| apply_mat(mat, p)).collect();
            let b = signature(&g, &turned, &m.elements);
            assert!(
                d_group(&g, &a, &b) < ROTATION_MUST_NOT_MOVE_GROUP,
                "rotation {r}: the group search failed to recover a group element"
            );
            // Not `r != 0`: the identity's *position* is a fact about
            // construction order, so naming it would reintroduce the coupling
            // the gate above just removed — and if a reordering put some other
            // element at index 0 and the identity at index 5, this loop would
            // fail at `r == 5` with the "geometry is dead" misdiagnosis this
            // whole change exists to eliminate. Ask which element it is.
            if !is_identity(mat) {
                assert!(
                    d_raw(&a, &b) > ROTATION_MUST_MOVE_RAW,
                    "rotation {r}: an unaligned descriptor did not move — geometry is dead"
                );
            }
        }
    }

    /// The rotation pair above passes for an `embed` that ignores the molecule
    /// entirely — a fixed configuration keyed on atom count still moves under
    /// rotation and still recovers. So the property that mutation violates is
    /// stated here separately: **same composition, different graph, different
    /// shape.**
    ///
    /// Checked across sizes because the separation scales with atom count
    /// (7.8 at n = 8, 14.4 at n = 12, 21.5 at n = 16) and a threshold chosen
    /// at one size should visibly clear it at the others.
    #[test]
    fn identical_composition_with_different_topology_gives_a_different_shape() {
        let g = Geodesic::<D>::build().unwrap();
        for n in [8_usize, 12, 16] {
            let path = Molecule::from_parts(vec![0; n], (1..n).map(|i| (i - 1, i)).collect());
            let star = Molecule::from_parts(vec![0; n], (1..n).map(|i| (0, i)).collect());
            assert_eq!(
                path.elements, star.elements,
                "n={n}: the fixtures must differ only in topology, or this proves nothing"
            );
            let a = signature(&g, &embed(&path), &path.elements);
            let b = signature(&g, &embed(&star), &star.elements);
            assert!(
                d_group(&g, &a, &b) > TOPOLOGY_MUST_CHANGE_SHAPE,
                "n={n}: a path and a star of identical composition were not separated — \
                 embed() is not reading the graph"
            );
        }
    }

    /// And the gate itself passes on the real pipeline, so the thresholds are
    /// not so tight that a working run trips them.
    #[test]
    fn the_geometry_liveness_gate_passes_on_a_working_pipeline() {
        assert_eq!(check_geometry_is_live(), Ok(()));
    }

    /// A non-finite measurement must fail every threshold, in both directions.
    ///
    /// This is the seam the geometry gate did not have. Written as
    /// `got <= min`, all three of its comparisons waved `NaN` through, because
    /// every ordering comparison against `NaN` is false — the identical defect
    /// `the_control_gate_rejects_controls_that_are_not_numbers` guards one
    /// level up. Without extracting the comparison there is nothing to inject
    /// a `NaN` into: the quantities are computed from real geometry that cannot
    /// currently produce one.
    #[test]
    fn a_non_finite_measurement_clears_no_threshold() {
        // `+inf` is the one a bare `got > min` would let past on its own terms
        // rather than by the NaN rule — it genuinely does exceed the floor. It
        // is still not a measurement.
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(!clears_floor(bad, 0.1), "{bad} cleared a floor");
            assert!(!under_ceiling(bad, 1e-6), "{bad} stayed under a ceiling");
        }

        // And the guards did not break the ordinary cases they wrap.
        assert!(clears_floor(5.98, 0.1));
        assert!(!clears_floor(0.05, 0.1));
        assert!(!clears_floor(0.1, 0.1), "the floor is exclusive");
        assert!(under_ceiling(3.9e-15, 1e-6));
        assert!(!under_ceiling(0.5, 1e-6));
        assert!(!under_ceiling(1e-6, 1e-6), "the ceiling is exclusive");
    }

    /// The real calibration check, run against the actual molecule generator
    /// rather than synthetic numbers. Both regimes, because a control that
    /// only works when sizes differ is not controlling for shape.
    #[test]
    fn the_controls_calibrate_in_both_regimes() {
        for size_matched in [false, true] {
            let cfg = Config {
                trials: 400,
                size_matched,
                ..Config::default()
            };
            let report = run(&cfg);
            assert!(
                report.check_controls().is_ok(),
                "size_matched={size_matched}: controls failed:\n{}",
                report.render()
            );
        }
    }
}
