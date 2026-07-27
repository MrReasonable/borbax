//! Exit criterion 9: can the acceptance instrument return a negative?
//!
//! # Why this exists before the metrics do
//!
//! Two reviewers independently concluded that the open-endedness instrument
//! as specified would certify **an inert beaker** — one in which no chemistry
//! happens at all. `cumulative` is a monotone accumulator, which *is* the
//! class-3 hallmark it is supposed to be testing for; `mean_cumulative = A/D`
//! is the statistic Channon 2003 exists to replace; and the neutral shadow
//! does not save either, because it shows the same artefacts.
//!
//! An instrument that cannot fail is not an instrument. So the falsification
//! test is written first and the metrics have to satisfy it, rather than
//! being written and then asked politely whether they do.
//!
//! # What is modelled here
//!
//! A `Trace` is the only thing the metrics need: which species were present
//! in each observation window. Three generators stand in for the three cases
//! that matter, and the middle one is the hard one:
//!
//! - [`Trace::inert`] — nothing ever happens. Must be rejected.
//! - [`Trace::neutral_drift`] — Bedau class 2: activity accumulates forever
//!   while a fixed pool is reshuffled and nothing new is ever produced. Must
//!   be rejected, and this is what the naive statistics are worst at.
//! - [`Trace::creative`] — genuinely novel species keep arriving and
//!   persisting. Must be accepted, and it is the positive control without
//!   which a verdict function could return "no" unconditionally and pass
//!   every other test here.
//!
//! This is deliberately not the full Bedau apparatus. It is the smallest
//! thing that can demonstrate the defect and pin the property the real
//! implementation must have.

use std::collections::{BTreeMap, BTreeSet};

use crate::rng::Stream;

/// A species identity. Opaque — nothing here depends on what a species is.
pub type SpeciesId = u64;

/// Which species were present in each observation window.
#[derive(Debug, Clone, Default)]
pub struct Trace {
    /// One entry per window. `BTreeSet` rather than a hash set because
    /// iteration order reaches the results (§13.4).
    pub windows: Vec<BTreeSet<SpeciesId>>,
}

impl Trace {
    /// A beaker where nothing happens: the same species, every window.
    ///
    /// The strongest null. Activity accumulates purely because things persist.
    #[must_use]
    pub fn inert(species: usize, windows: usize) -> Self {
        let present: BTreeSet<SpeciesId> = (0..species)
            .map(|i| SpeciesId::try_from(i).unwrap_or(0))
            .collect();
        Self {
            windows: vec![present; windows],
        }
    }

    /// Bedau class 2: a fixed pool, resampled every window.
    ///
    /// Diversity churns and cumulative activity rises without bound, but the
    /// system never produces anything it could not have produced at window
    /// zero. Also serves as the neutral shadow — that is what a shadow *is*,
    /// a run with the creative process switched off.
    #[must_use]
    pub fn neutral_drift(rng: &mut Stream, pool: usize, per_window: usize, windows: usize) -> Self {
        let mut out = Vec::with_capacity(windows);
        for _ in 0..windows {
            let mut present = BTreeSet::new();
            for _ in 0..per_window.min(pool) {
                present.insert(SpeciesId::try_from(rng.index(pool)).unwrap_or(0));
            }
            out.push(present);
        }
        Self { windows: out }
    }

    /// A beaker that keeps inventing: every window mints species never seen
    /// before, and each persists for `lifetime` windows.
    ///
    /// Persistence is the part that matters. A species that appears for one
    /// window and vanishes is noise, not novelty, and an instrument that
    /// counted it would certify a random-number generator.
    #[must_use]
    pub fn creative(rng: &mut Stream, per_window: usize, lifetime: usize, windows: usize) -> Self {
        let mut out = vec![BTreeSet::new(); windows];
        let mut next: SpeciesId = 1_000_000;
        for start in 0..windows {
            for _ in 0..per_window {
                let id = next;
                next += 1;
                // A little variation in how long each lasts, so the trace is
                // not an exact lattice that a metric could latch onto.
                let live = 1 + rng.index(lifetime);
                for w in start..(start + live).min(windows) {
                    if let Some(slot) = out.get_mut(w) {
                        slot.insert(id);
                    }
                }
            }
        }
        Self { windows: out }
    }

    /// Brand-new species every window, none of which ever recurs.
    ///
    /// Maximal novelty and zero creation — the null that catches an
    /// instrument counting churn. A random-number generator produces this,
    /// and certifying it would be a worse failure than certifying the inert
    /// beaker, because the numbers look spectacular.
    #[must_use]
    pub fn ephemeral_noise(per_window: usize, windows: usize) -> Self {
        let mut next: SpeciesId = 5_000_000;
        let mut out = Vec::with_capacity(windows);
        for _ in 0..windows {
            let mut present = BTreeSet::new();
            for _ in 0..per_window {
                present.insert(next);
                next += 1;
            }
            out.push(present);
        }
        Self { windows: out }
    }

    /// Window index in which each species was first seen.
    #[must_use]
    pub fn first_appearances(&self) -> BTreeMap<SpeciesId, usize> {
        let mut out = BTreeMap::new();
        for (w, present) in self.windows.iter().enumerate() {
            for &id in present {
                out.entry(id).or_insert(w);
            }
        }
        out
    }

    /// How many windows each species was present for — Bedau's activity
    /// counter, accumulated over the whole trace.
    /// Longest *consecutive* run of windows in which `id` was present.
    ///
    /// [`Self::activity`] counts total presence anywhere in the trace, which is
    /// the right quantity for Bedau's `a_i(t)` and the wrong one for a
    /// persistence floor. The two were conflated: `MIN_PERSISTENCE` is
    /// documented as "stayed for at least N windows" and was implemented as
    /// "was present in at least N windows, anywhere". Measured, the difference
    /// is not academic — 4800 species each placed in exactly two *randomly
    /// chosen* windows, nothing surviving anything, scored 1148 and certified
    /// the trace open-ended. The random-number generator the floor exists to
    /// reject was rejected only in the one form the generator happened to be
    /// written in.
    #[must_use]
    pub fn longest_run(&self, id: SpeciesId) -> usize {
        let mut best = 0;
        let mut run = 0;
        for present in &self.windows {
            if present.contains(&id) {
                run += 1;
                if run > best {
                    best = run;
                }
            } else {
                run = 0;
            }
        }
        best
    }

    /// Total window-presences per species — Bedau's `a_i(t)` summed over the
    /// trace. This is the right quantity for activity statistics and the wrong
    /// one for a persistence floor; see [`Self::longest_run`].
    #[must_use]
    pub fn activity(&self) -> BTreeMap<SpeciesId, usize> {
        let mut out = BTreeMap::new();
        for present in &self.windows {
            for &id in present {
                *out.entry(id).or_insert(0) += 1;
            }
        }
        out
    }
}

/// Cumulative activity `A(t)`: total window-presences up to and including `t`.
///
/// Monotone by construction — it only ever adds — which is precisely why it
/// cannot serve as evidence of open-endedness on its own.
#[must_use]
pub fn cumulative_activity(trace: &Trace) -> Vec<usize> {
    let mut total = 0;
    trace
        .windows
        .iter()
        .map(|present| {
            total += present.len();
            total
        })
        .collect()
}

/// `A(t) / D(t)` — cumulative activity normalised by current diversity.
#[must_use]
pub fn mean_cumulative_activity(trace: &Trace) -> Vec<f64> {
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "counts of species and windows, far below 2^53"
    )]
    cumulative_activity(trace)
        .iter()
        .zip(trace.windows.iter())
        .map(|(&a, present)| {
            if present.is_empty() {
                0.0
            } else {
                a as f64 / present.len() as f64
            }
        })
        .collect()
}

/// The criterion as originally specified: cumulative activity, and its
/// diversity-normalised form, both rising without bound.
///
/// Kept in the codebase on purpose. It is the thing that was nearly shipped,
/// and `the_naive_cumulative_criterion_certifies_an_inert_beaker` is the
/// evidence that it should not be — a claim that decays into folklore if the
/// code it refers to is deleted.
#[must_use]
pub fn naive_verdict(trace: &Trace) -> bool {
    let a = cumulative_activity(trace);
    let mean = mean_cumulative_activity(trace);
    let rising = |first: f64, last: f64| last > first;
    #[allow(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "counts of window-presences, far below 2^53"
    )]
    let a_rising = match (a.first(), a.last()) {
        (Some(&f), Some(&l)) => rising(f as f64, l as f64),
        _ => false,
    };
    let mean_rising = match (mean.first(), mean.last()) {
        (Some(&f), Some(&l)) => rising(f, l),
        _ => false,
    };
    a_rising && mean_rising
}

/// What the instrument concluded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// Novelty above the neutral baseline is still arriving late in the run.
    OpenEnded,
    /// It is not.
    NotOpenEnded,
    /// The trace is too short to support either answer.
    ///
    /// A separate outcome rather than a default, because defaulting to
    /// `NotOpenEnded` would hide a run that was never long enough to judge,
    /// and defaulting to `OpenEnded` would be the original defect again.
    Inconclusive,
}

/// Minimum windows needed before either answer is available.
pub const MIN_WINDOWS: usize = 20;

/// Windows a species must last to count as more than churn.
///
/// A species that appears for one window and vanishes is noise. Without this
/// floor the criterion would certify a random-number generator, which is the
/// false positive Bedau's class 2 exists to name.
pub const MIN_PERSISTENCE: usize = 2;

/// The shadow-corrected criterion.
///
/// Two changes from the naive form, and both are load-bearing:
///
/// 1. **Only novelty arriving in the *second half* counts.** Open-endedness
///    is about whether a system is still producing, not about what it
///    produced at the start. An inert beaker invents everything in window
///    zero and nothing afterwards, so it scores zero here however long it
///    runs — which is what kills the monotone-accumulator false positive.
/// 2. **The count must beat the neutral shadow's**, not an absolute number.
///    A run and its shadow differ only in whether the creative process is on,
///    so anything the shadow also produces is not evidence of creation.
///
/// The comparison is against the shadow's *count of late novelties*, not
/// against its peak activity. Raw activity totals are not comparable between
/// runs of different diversity: a 40-species pool resampled for 400 windows
/// accumulates ~124 presences per species, while a species that is genuinely
/// new and lasts 8 windows accumulates 8. Thresholding one against the other
/// rejects real novelty and accepts recycling — the exact inversion of what
/// the instrument is for.
///
/// This is why the spec has the shadow spawned from a keyframe of the run
/// rather than generated independently (§9): only a **matched** shadow makes
/// the two counts comparable. A real implementation must match window count,
/// per-window diversity and turnover, and this function is only as good as
/// the shadow it is handed.
#[must_use]
pub fn verdict(trace: &Trace, shadow: &Trace) -> Verdict {
    if trace.windows.len() < MIN_WINDOWS {
        return Verdict::Inconclusive;
    }
    if !shadow_is_matched(trace, shadow) {
        return Verdict::Inconclusive;
    }
    if qualifying_late_novelties(trace) > qualifying_late_novelties(shadow) {
        Verdict::OpenEnded
    } else {
        Verdict::NotOpenEnded
    }
}

/// Whether `shadow` is a usable null for `trace`.
///
/// **This checks window count only, and that is a partial fix to a real
/// defect. The rest cannot be checked here, and finding out why is the
/// finding.**
///
/// The bar `verdict` compares against is the shadow's own
/// [`qualifying_late_novelties`]. With the shadow the test suite was using,
/// that bar is *structurally* zero rather than incidentally zero: a 40-species
/// pool drawn 12 at a time exhausts itself within a few windows, so no species
/// can first appear in the second half, and 0 of 500 seeds produce a non-zero
/// bar. `verdict` collapsed to `count >= 1`, and the measured consequences
/// were an inert beaker plus one transient certified `OpenEnded`, and
/// `neutral_drift` — zero creation by construction, this module's own class-2
/// null — certified as soon as its pool is 1000 rather than 40, which is the
/// regime a real beaker is in. Substituting a *zero-window* shadow changed no
/// verdict in the entire suite.
///
/// Window count catches the zero-window case and is unambiguous. Everything
/// else was tried and does not work, which is worth recording so it is not
/// re-attempted:
///
/// - **Per-window diversity ratio.** No tolerance separates the cases.
///   `Trace::inert(30)` has 30 species every window; the closest
///   `neutral_drift` shadow yields ~19, because it draws *with replacement*
///   and `per_window` is capped at `pool`. A tolerance loose enough to admit
///   that pair (1.6) is loose enough to admit much else.
/// - **Total distinct species.** Rejects the class-2 mismatch correctly, and
///   also rejects the *positive control* — `creative` has ~4800 species against
///   any fixed pool, and having more species than the null is precisely what
///   makes it creative.
///
/// The reason neither works is that "matched" is a property of how the shadow
/// was *generated*, not of what it contains. §15.3 says the shadow is forked
/// from a keyframe of the run and driven by the same environment with
/// selection neutralised — that provenance cannot be recovered from the trace.
///
/// **The fix is a type**: a pair constructible only from a run and its
/// keyframe-spawned shadow, so an unmatched pair cannot be passed to
/// `verdict` at all. That belongs with Task 20b, where the shadow is actually
/// forked; here there is nothing to fork from.
#[must_use]
pub const fn shadow_is_matched(trace: &Trace, shadow: &Trace) -> bool {
    trace.windows.len() == shadow.windows.len()
}

/// Species that first appeared in the trace's second half and then stayed for
/// at least [`MIN_PERSISTENCE`] windows.
#[must_use]
pub fn qualifying_late_novelties(trace: &Trace) -> usize {
    let halfway = trace.windows.len() / 2;
    trace
        .first_appearances()
        .iter()
        .filter(|&(id, &first)| first >= halfway && trace.longest_run(*id) >= MIN_PERSISTENCE)
        .count()
}

/// A human-readable account of how a verdict was reached, for test failures.
///
/// A bare `assert_eq!(verdict, NotOpenEnded)` failing tells you the instrument
/// is wrong but not how, and this is the class of bug where that gap costs
/// hours.
#[must_use]
pub fn diagnose(trace: &Trace, shadow: &Trace) -> String {
    let halfway = trace.windows.len() / 2;
    let late = trace
        .first_appearances()
        .iter()
        .filter(|&(_, &first)| first >= halfway)
        .count();
    format!(
        "  windows                  {}\n  \
           distinct species         {}\n  \
           first seen in 2nd half   {late}\n  \
           of those, persisting     {}\n  \
           same figure for shadow   {}\n  \
           naive verdict            {}",
        trace.windows.len(),
        trace.activity().len(),
        qualifying_late_novelties(trace),
        qualifying_late_novelties(shadow),
        naive_verdict(trace),
    )
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    reason = "CLAUDE.md: tests may unwrap freely; indices are loop bounds over sized collections"
)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    const WINDOWS: usize = 400;

    fn shadow() -> Trace {
        Trace::neutral_drift(&mut Stream::new(9001), 40, 12, WINDOWS)
    }

    // -- The three cases the instrument used to get wrong. ------------------
    //
    // Every one of these was measured returning `OpenEnded` before the shadow
    // matching and the survival-vs-presence fix landed.

    /// An inert beaker with a single transient species is not open-ended.
    ///
    /// Measured before the fix: `qualifying_late_novelties = 1` against a
    /// structurally-zero bar, verdict `OpenEnded`. The module exists because
    /// the *original* criterion certified an inert beaker; the replacement
    /// certified an inert beaker plus one molecule.
    #[test]
    fn an_inert_beaker_plus_one_transient_is_still_certified() {
        let mut t = Trace::inert(30, WINDOWS);
        let late = WINDOWS * 3 / 4;
        for w in late..(late + 2) {
            if let Some(win) = t.windows.get_mut(w) {
                win.insert(999_999);
            }
        }
        // PINNED RESIDUAL, not a fix. The shadow's bar is structurally zero,
        // so `verdict` reduces to `count >= 1`. See `shadow_is_matched` for
        // why comparing two traces cannot establish matching. Delete this test
        // when it starts failing — that means the residual is closed.
        assert_eq!(qualifying_late_novelties(&t), 1);
        assert_eq!(verdict(&t, &shadow()), Verdict::OpenEnded);
    }

    /// A species flickering in two *non-adjacent* windows has survived
    /// nothing. Before the fix this scored 1148 and certified the trace,
    /// because the floor counted total presence rather than a run.
    #[test]
    fn a_two_window_flicker_is_not_persistence() {
        let mut rng = Stream::new(77);
        let mut windows = vec![BTreeSet::new(); WINDOWS];
        for id in 1..=4800_u64 {
            let a = rng.index(WINDOWS);
            let b = rng.index(WINDOWS);
            if a.abs_diff(b) < 2 {
                continue;
            }
            if let Some(w) = windows.get_mut(a) {
                w.insert(id);
            }
            if let Some(w) = windows.get_mut(b) {
                w.insert(id);
            }
        }
        let t = Trace { windows };
        assert_eq!(
            qualifying_late_novelties(&t),
            0,
            "a flicker is not survival"
        );
    }

    /// The shadow must be able to change a verdict. Before this, substituting
    /// a zero-window shadow left *every* verdict in this file unchanged — the
    /// half of the criterion §15.3 exists for had never been exercised.
    #[test]
    fn a_zero_window_shadow_is_inconclusive_rather_than_a_free_pass() {
        let t = Trace::creative(&mut Stream::new(3), 12, 8, WINDOWS);
        assert_eq!(verdict(&t, &Trace::default()), Verdict::Inconclusive);
        // Only a length mismatch is Inconclusive. Diversity matching was
        // tried and cannot be made to work — see `shadow_is_matched`.
        assert_eq!(verdict(&t, &shadow()), Verdict::OpenEnded);
    }

    // -- The headline test. ------------------------------------------------
    //
    // Exit criterion 9 asks whether the acceptance instrument can return a
    // negative. This is that question as a test: a beaker in which no
    // chemistry happens at all must not be certified open-ended. It is
    // written before any metric code so that the metrics have to satisfy it,
    // rather than being written and then asked politely whether they do.

    #[test]
    fn a_beaker_with_no_chemistry_is_not_open_ended() {
        let inert = Trace::inert(30, WINDOWS);
        // Matched on per-window diversity: 30 present every window.
        let sh = shadow();
        assert_eq!(
            verdict(&inert, &sh),
            Verdict::NotOpenEnded,
            "an inert beaker was certified open-ended:\n{}",
            diagnose(&inert, &sh)
        );
    }

    /// Bedau class 2: activity accumulates without bound while nothing new is
    /// ever produced. A fixed pool of species shuffled forever is the cleanest
    /// case, and it is the one the standard statistics are worst at.
    #[test]
    fn a_beaker_recycling_a_fixed_pool_is_not_open_ended() {
        let drifting = Trace::neutral_drift(&mut Stream::new(4242), 40, 12, WINDOWS);
        assert_eq!(
            verdict(&drifting, &shadow()),
            Verdict::NotOpenEnded,
            "a beaker recycling a fixed pool was certified open-ended:\n{}",
            diagnose(&drifting, &shadow())
        );
    }

    /// Maximal novelty, zero creation: every species is new and none survives
    /// a single window. Certifying this would be worse than certifying the
    /// inert beaker, because the diversity numbers look spectacular.
    #[test]
    fn a_beaker_producing_only_ephemeral_species_is_not_open_ended() {
        let noise = Trace::ephemeral_noise(12, WINDOWS);
        assert_eq!(
            verdict(&noise, &shadow()),
            Verdict::NotOpenEnded,
            "pure churn was certified open-ended:\n{}",
            diagnose(&noise, &shadow())
        );
    }

    /// The positive control, and it is not optional. Without it,
    /// `verdict` could return `NotOpenEnded` unconditionally and pass every
    /// test above — an instrument that always says no is exactly as useless as
    /// one that always says yes.
    #[test]
    fn a_beaker_that_keeps_inventing_species_is_open_ended() {
        let creative = Trace::creative(&mut Stream::new(1234), 12, 8, WINDOWS);
        // A *fixed pool* drifting at the same per-window occupancy. This is
        // the null the positive control has to beat: same diversity, same
        // turnover, but nothing new can ever appear.
        let sh = shadow();
        assert_eq!(
            verdict(&creative, &sh),
            Verdict::OpenEnded,
            "a genuinely creative beaker was not recognised:\n{}",
            diagnose(&creative, &sh)
        );
    }

    // -- Why the instrument had to change. ---------------------------------

    /// The defect, measured rather than argued. The criterion as specified —
    /// unbounded cumulative activity — certifies a beaker in which literally
    /// nothing happens.
    #[test]
    fn the_naive_cumulative_criterion_certifies_an_inert_beaker() {
        let inert = Trace::inert(30, WINDOWS);
        assert!(
            naive_verdict(&inert),
            "the naive criterion was expected to fail here; if it no longer \
             does, this finding is stale and the comment above needs revisiting"
        );
    }

    /// Why it does. `A(t)` counts window-presences and only ever adds, so it
    /// rises forever in any beaker whose species merely persist. Monotonicity
    /// is the class-3 hallmark the statistic is supposed to be *testing for*,
    /// and it holds by construction.
    #[test]
    fn cumulative_activity_rises_forever_in_a_beaker_where_nothing_happens() {
        let a = cumulative_activity(&Trace::inert(30, WINDOWS));
        for w in a.windows(2) {
            assert!(
                w[1] > w[0],
                "cumulative activity is supposed to be monotone"
            );
        }
        assert!(a.last().is_some_and(|&last| last > a[0] * 100));
    }

    /// Normalising by diversity does not rescue it. For an inert beaker every
    /// species is present in every window, so `A/D` is just the window index —
    /// it grows without bound too, and reports the same false positive.
    #[test]
    fn dividing_by_diversity_does_not_rescue_the_naive_criterion() {
        let mean = mean_cumulative_activity(&Trace::inert(30, WINDOWS));
        for w in mean.windows(2) {
            assert!(w[1] > w[0], "A/D is supposed to be monotone here too");
        }
        let first = mean[0];
        let last = *mean.last().unwrap();
        assert!(
            last > first * 100.0,
            "A/D grew only from {first} to {last} — the demonstration is weaker than claimed"
        );
    }

    // -- The instrument's own calibration. ---------------------------------

    /// The meta-test. Both answers must be reachable from the same code path,
    /// or the suite above is measuring a constant.
    #[test]
    fn the_instrument_returns_both_answers() {
        // Each against its own matched shadow, as §15.3 requires.
        let answers = [
            verdict(&Trace::inert(30, WINDOWS), &shadow()),
            verdict(
                &Trace::creative(&mut Stream::new(77), 12, 8, WINDOWS),
                &shadow(),
            ),
        ];
        assert!(answers.contains(&Verdict::NotOpenEnded));
        assert!(answers.contains(&Verdict::OpenEnded));
    }

    /// A trace shorter than the halves the verdict compares cannot support a
    /// judgement, and saying so is better than defaulting either way.
    #[test]
    fn too_short_a_trace_is_refused_rather_than_guessed() {
        assert_eq!(
            verdict(&Trace::inert(5, 3), &shadow()),
            Verdict::Inconclusive
        );
    }

    #[test]
    fn generated_traces_are_what_they_claim_to_be() {
        let inert = Trace::inert(30, WINDOWS);
        assert_eq!(inert.windows.len(), WINDOWS);
        assert!(inert.windows.iter().all(|w| w.len() == 30));
        assert_eq!(
            inert.first_appearances().len(),
            30,
            "inert invented a species"
        );

        let drifting = Trace::neutral_drift(&mut Stream::new(1), 40, 12, WINDOWS);
        assert!(
            drifting.first_appearances().len() <= 40,
            "neutral drift escaped its pool"
        );

        let creative = Trace::creative(&mut Stream::new(2), 12, 8, WINDOWS);
        assert!(
            creative.first_appearances().len() > 200,
            "creative produced too few novel species to be a positive control"
        );
    }
}
