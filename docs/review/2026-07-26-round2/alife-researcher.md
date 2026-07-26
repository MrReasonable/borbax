# ALife literature review — round 2 (re-run with tools)

**Scope:** spec §2.1, §2.7, §9.3, §9.6, §15.2–15.4, §23; plan Tasks 14, 17, 18, 20, 20b, 21.

**Standing constraint honoured:** nothing below returns element data, reaction data,
structural data, sequence data, or any Borbax↔real mapping. Every result cited is a
property of maps, graphs, likelihoods or estimators. Where a source's original subject
matter was chemical or biological, only its *computational* content is reported.

**What changed from the previous pass.** That pass was dispatched without write tools and
its findings were hand-transcribed. This one has `Bash`, so the arithmetic is computed
rather than asserted, and the papers are downloaded and read rather than recalled. Two of
the previous pass's own findings turn out to be wrong and are withdrawn below
(`alr-115`, and the AlChemy figures in `alr-101` are now settled with a different answer
than *either* the spec's or my previous guess).

**Sources actually opened this pass** (as opposed to recalled):

| Source | How obtained | Confidence |
|---|---|---|
| Mathis, Patel, Weimer & Forrest, "Return to AlChemy", arXiv:2408.12137v2 | PDF downloaded; text extracted; **Figure 4(B) raster extracted and read** | Verified from source |
| Channon, "Improving and Still Passing the ALife Test", ALIFE VIII | PDF downloaded from author's site, full text | Verified from source |
| Dolson, Vostinar, Wiser & Ofria, "The MODES Toolbox", *Artificial Life* 25(1) | LaTeX source | Verified from source |
| Hordijk & Steel, "A Concise and Formal Definition of RAF Sets and the RAF Algorithm", arXiv:2303.01809 | PDF, full text | Verified from source |
| Packard, Bedau, Channon, Ikegami, Rasmussen, Stanley & Taylor, OEE-II editorial, arXiv:1909.04430 | PDF, full text | Verified from source |
| SakanaAI/asal `asal_metrics.py`, `foundation_models/clip.py` | GitHub API, verbatim | Verified from source |
| Wiser, Ribeck & Lenski 2013, *Science* 342:1364 | **Paywalled — 403.** Abstract + secondary descriptions only | **Could not verify** |
| Bedau, Snyder & Packard 1998 | **Not opened.** Relying on Channon 2003's verbatim reproduction of its Table 1 and footnote 1 | Second-hand, but from a source that quotes it directly |

---

## Summary table

| id | sev | category | where |
|---|---|---|---|
| alr-100 | **Critical** | numerical-error | plan 2997–3009 (`ActivityStats`) |
| alr-101 | **Critical** | misattribution | spec §2.1 L59; `.claude/agents/alife-researcher.md` L44–45 |
| alr-102 | High | numerical-error | plan 3004–3005 (`mean_cumulative`) |
| alr-103 | High | spec-divergence | plan 3002–3003 (`diversity`) |
| alr-104 | High | spec-divergence | plan 3006–3007 (`new_activity`) |
| alr-105 | High | scope-gap | plan 2983 vs 3036 (step ordering) |
| alr-106 | High | numerical-error | plan 3017–3034 (`shape_novelty`) |
| alr-107 | High | numerical-error | plan 3063–3069 (plateau RSS) |
| alr-108 | High | spec-divergence | plan 3063; spec §15.4 L870 |
| alr-109 | High | scope-gap | plan 1896–1902, 2886, 2893 |
| alr-110 | High | type-inconsistency | plan 1871, 2483, 2721 |
| alr-111 | Medium | numerical-error | plan 3069 (AICc `k`) |
| alr-112 | Medium | misattribution | plan 3011–3013 (MODES filter) |
| alr-113 | Medium | misattribution | plan 3018–3019 ("Lehman-Stanley") |
| alr-114 | Medium | misattribution | spec §2.7 L110; §20 L1042 |
| alr-115 | Medium | *correction to my own prior finding* | spec §2.7 L106 |
| alr-116 | Medium | numerical-error | plan 2948–2949 (neutral network) |
| alr-117 | Medium | scope-gap | spec §15.2 L839, §9.6 L511; plan 2452 |
| alr-118 | Medium | scope-gap | plan 3047–3054 |
| alr-119 | Medium | scope-gap | plan 3036–3041 |
| alr-120 | Low | type-inconsistency | plan 2440, 2288 |
| alr-121 | Low | scope-gap | spec §23 L1147 |
| alr-122 | Low | scope-gap | plan 2955–2958 |
| alr-123 | Low | numerical-error | plan 3063 |
| alr-124 | Low | numerical-error | plan 3065, 3077 |
| alr-125 | Low | numerical-error | plan 3069 |

---

## Verified correct — the previous round's fixes hold

### alr-000 · RAF closure fix is correct · *verified by execution*

I transcribed the plan's `find_raf` (plan 2452–2520) into a standalone Rust file and ran it.
Scratch file: `/private/tmp/claude-501/-Users-iandominey-projects-borbax/572831ef-ed3e-4266-af6e-48db0ed426dc/scratchpad/raf.rs`.

```
P1 bimolecular A+B->C, B unreachable      -> None   (correct)
P2 unimolecular encoded as [0,0], cat B   -> Some([0])  (correct)
```

The fix at plan 2483 (`missing = r.reactants.len()`, counting *all* reactants including
food ones) matches `Incidence::build`'s indexing of every reactant occurrence, so the
counter reaches zero exactly when all reactants are present. The duplicate-reactant case
works because `Incidence::build` writes two entries and the queue pops the species once,
decrementing twice.

Against source: this is Algorithm 1 of Hordijk & Steel (arXiv:2303.01809), verbatim —

> "Starting with the full set of reactions R′ = R, the algorithm repeatedly calculates the
> closure of the food set relative to the current reaction set R′, and then removes from
> R′ all reactions that have none of their catalysts or not all of their reactants in this
> closure. This is repeated until no more reactions can be removed. If upon termination of
> the algorithm R′ is non-empty, then R′ is the unique maximal RAF set (maxRAF) contained
> in Q."

The plan's `retain(reactants_ok && catalysed)` is the De Morgan dual of the paper's
removal condition. Correct. Note `ComputeClosure` (their Algorithm 2) deliberately does
*not* consult catalysts — that is what distinguishes RAF from CAF — and the plan's closure
correctly does not either.

### alr-001 · The four doc-comment claims check out

- **Maximality well-defined via closure under union** — confirmed: "the unique maximal RAF
  set (maxRAF) contained in Q, i.e., a RAF that contains every other RAF in Q as a subset."
- **RAF not CAF** — confirmed; arXiv:2303.01809 Algorithm 3 gives the CAF variant, which
  additionally requires a catalyst present during expansion.
- **Structural not dynamic** — confirmed; the CRS is `{X, R, C, F}` with no counts anywhere.
- **Inhibition is NP-hard** — the general RAF-with-inhibition decision problem is
  NP-complete (Mossel & Steel, *J. Theor. Biol.* 233(3):327–336, 2005, "Random biochemical
  networks: the probability of self-sustaining autocatalysis", DOI
  10.1016/j.jtbi.2004.10.011). I did **not** open Mossel & Steel this pass; this is from
  secondary descriptions plus recollection. It is however consistent with the tractability
  caveat in the recent literature (bounded inhibitor count remains tractable).

Complexity, from arXiv:2303.01809 verbatim: "A straightforward computational complexity
analysis of the RAF algorithm gives a worst-case running time of O(|X||R|³) … In practice,
the average running time on a simple polymer-based model of CRSs turns out to be
sub-quadratic." Worth putting in the doc comment — `O(|X||R|³)` is not free, and the plan
calls `find_raf` from the battery.

### alr-002 · `1.0 - u1` is correct (plan 2709–2712)

`next_f64` is uniform on [0,1). Using `1.0 - u1` moves the draw to (0,1]. The trade is
`τ = +∞` at `u1 == 0` for `τ = 0` at `u1` within 2⁻⁵³ of 1, which stalls the clock for one
event instead of ending the run. That reason — not the guard itself — should be the test's
assertion string.

### alr-003 · `n(n-1)/2` is correct (plan 2721)

Gillespie, D. T. (1977), "Exact Stochastic Simulation of Coupled Chemical Reactions",
*J. Phys. Chem.* 81(25):2340–2361, DOI 10.1021/j100540a008. Computed magnitudes of the
error the plan is guarding against:

```
n        n(n-1)/2        n^2/2      ratio
2             1.0          2.0     2.0000x
5            10.0         12.5     1.2500x
10           45.0         50.0     1.1111x
100        4950.0       5000.0     1.0101x
```

A factor of 2 at n=2 is exactly the regime where a RAF nucleates. The plan's note is right.

---

## Critical

### alr-100 · `ActivityStats.cumulative` is unbounded by construction

**Where:** plan `2997–3009`.
**Claim:** the plan describes `per_species` as a monotone cumulative counter and
`cumulative` as its sum over all species; in Bedau/Channon's definition `a_i(t)` is **zero
whenever component *i* is absent**, and `A_cum` sums only over *extant* components.
**Confidence:** verified from source (Channon 2003, equations 1, 2, 4, 12).

Channon 2003, verbatim:

> Activity increment (by presence). `Δᵢ(t) = 1 if component i exists at t, 0 otherwise` (1)
>
> Evolutionary Activity of a component.
> `aᵢ(t) = Σ_{τ=0}^{t} Δᵢ(τ) if component i exists at t, 0 otherwise` (2)
>
> `D(t) = #{i : aᵢ(t) > 0}` (3)   `A_cum(t) = Σᵢ aᵢ(t)` (4)

and for the normalised form, equation (12):

> `A^N_cum(t) = Σ_{i: component i exists in the real run at t} a^N_i(t)`

The `0 otherwise` clause is load-bearing. It makes `A_cum(t)` a sum over the *currently
present* components only, so it can and does fall.

**Why it matters.** Table 1 of Channon 2003 (reproducing Bedau, Snyder & Packard 1998 with
rows 3b and 3c added):

```
CLASS  EVOLUTIONARY DYNAMICS        D          A_new     Ā_cum
1      none                         bounded    zero      zero
2      bounded                      bounded    positive  bounded
3a     unbounded (D)                unbounded  positive  bounded
3b     unbounded (Ā_cum)            bounded    positive  unbounded
3c     unbounded (D & Ā_cum)        unbounded  positive  unbounded
```

Every open-ended class is defined by *unbounded* activity. If `per_species` is a monotone
counter and `cumulative` sums over all species ever interned, then `cumulative` is monotone
non-decreasing **by construction** — it cannot be bounded, so the statistic returns "class
3" for a beaker containing nothing but inert feedstock. The metric answers the question it
exists to ask, before the chemistry gets a vote. This is the same defect the plan correctly
diagnoses for raw species counts, one level up.

**Fix.** Keep the monotone history separately from the reported statistic, and gate on
presence:

```rust
/// Per-species activity counters, following Channon (2003) eqs. (1)–(5),
/// which reproduce Bedau, Snyder & Packard (1998) Table 1 with rows 3b/3c
/// made explicit.
///
/// `history[i]` is the monotone sum of activity increments Σ Δ_i(τ).
/// `a_i(t)` — the quantity the statistics are defined over — is
/// `history[i]` **gated by current presence**: it is zero while the species
/// is absent. Without the gate `cumulative` is monotone non-decreasing and
/// therefore unbounded by construction, which is precisely the class-3
/// hallmark the statistic exists to test for (spec §2.7, §15.2).
pub struct ActivityStats {
    /// Σ Δ_i(τ) over all τ. Monotone. **Not** a_i(t); not reported directly.
    history: Vec<f64>,
    /// a_i(t) = history[i] if present at t, else 0.0. Recomputed each snapshot.
    per_species: Vec<f64>,
    pub diversity: f64,
    pub cumulative: f64,
    pub mean_cumulative: f64,
    pub median_cumulative: f64,
    pub new_activity: f64,
}

impl ActivityStats {
    /// `present` is the extant species at this snapshot, sorted ascending.
    /// Sorted, not a HashSet: this feeds a result-affecting statistic.
    pub fn snapshot(&mut self, present: &[SpeciesId], counts: &[f64]) {
        for (i, a) in self.per_species.iter_mut().enumerate() {
            *a = 0.0;
            let _ = i;
        }
        for &s in present {
            let i = s.0 as usize;
            if counts[i] > 0.0 {
                self.history[i] += 1.0;          // Δ_i(t) = 1 (eq. 1)
                self.per_species[i] = self.history[i]; // a_i(t) (eq. 2)
            }
        }
        // eq. (3): D counts components with a_i > 0 — i.e. present ones.
        self.diversity = self.per_species.iter().filter(|a| **a > 0.0).count() as f64;
        // eq. (4): sum over extant components only.
        self.cumulative = kahan_ordered_sum(&self.per_species);
    }
}
```

**Test that would catch it:**

```rust
/// Channon (2003) eq. (2): a_i(t) is zero while component i is absent, so
/// A_cum(t) sums over extant components only and CAN fall. A monotone
/// accumulator makes A_cum unbounded by construction, which is the class-3
/// hallmark — the statistic would return its own answer.
#[test]
fn cumulative_activity_falls_when_species_go_extinct() {
    let mut s = ActivityStats::with_capacity(3);
    let counts = [1.0, 1.0, 1.0];
    for _ in 0..10 {
        s.snapshot(&[SpeciesId(0), SpeciesId(1), SpeciesId(2)], &counts);
    }
    let peak = s.cumulative;
    assert!(peak > 0.0);

    // species 1 and 2 die out
    let counts = [1.0, 0.0, 0.0];
    s.snapshot(&[SpeciesId(0)], &counts);
    assert!(
        s.cumulative < peak,
        "A_cum did not fall on extinction ({peak} -> {}): per-species activity \
         is being accumulated monotonically, so A_cum is unbounded by \
         construction and the class-3 test is vacuous",
        s.cumulative
    );
    assert_eq!(s.diversity, 1.0, "D counts extant components (eq. 3)");
}
```

---

### alr-101 · The AlChemy figures in the spec are wrong. Settled from the paper.

**Where:** spec §2.1 line 59; `.claude/agents/alife-researcher.md` lines 44–45.
**Claim:** the spec says organisations "coexist only about 16% of the time and mutually
destroy each other around 60% of the time." The paper's Figure 4(B) table reads
**Mutual Destruction ~68%, Domination ~27%, Coexistence ~5%, N = 455**.
**Confidence:** verified from source. The table is a raster inside the figure and does not
survive text extraction, so I extracted `Im4.jpg` from page 10 of the PDF and read it
directly (both polarities, to rule out a CMYK inversion artefact — identical).

This also **withdraws my own previous pass's guess** of ~13% / ~17% / ~70%. That was a
recollection of a render and it was wrong too. The numbers above are read off the figure.

**Citation:** Cole Mathis, Devansh Patel, Westley Weimer & Stephanie Forrest,
"Self Organization in Computation & Chemistry: Return to AlChemy", *Chaos* 34(9):093142
(2024), DOI 10.1063/5.0207358; arXiv:2408.12137v2.

**Three conditions on the result that the spec's one-sentence summary loses, and all three
matter for Borbax:**

1. **The experiment combined organisations from *different* simulations.** Verbatim: "We
   started with two L1 simulations (L and R) each with 1000 expressions remaining at the
   end of a series of over 10⁶ collisions and five perturbations, and combined them into a
   single simulation with the maximum number of objects set to 2500." The conclusion is
   stated with that condition attached: "These results show that the organizations produced
   by AlChemy can possibly coexist, but **it rarely occurs for organizations evolved in
   different simulations**."

2. **The other route to L2 was not attempted.** Verbatim: "In the original work, Fontana &
   Buss identified L2 organizations in two different ways, first by combining two
   independent L1 organizations found in different simulations, and then by identifying an
   L1 simulation in which two separable organizations emerged. **We did not attempt the
   latter.** … It is possible that the results we reported in the previous sections were in
   fact L2 organizations with separable internal structure."

   So this is a non-replication of *one* of the two original routes. The within-one-system
   route — which is the one Borbax's V2/V3 composition question actually resembles, since
   Borbax protocells would arise inside a single continuing chemistry — is **untested by
   this paper**. The countermeasure in §2.7 ("never design a version on the assumption that
   lower-level units will stack") remains correct and I would not weaken it. But the spec
   currently presents the evidence as stronger and more general than it is.

3. **The categories are defined by a similarity threshold**, not by a structural criterion.
   Verbatim: "(i) Dominance—the organization retains non-zero similarity to one input but
   not the other, (ii) Coexistence—the organization retains an **average similarity > 0.1**
   to both inputs, and (iii) Mutual Destruction—the organization retains **< 0.1** average
   similarity to both inputs." The similarity measure is the Jaccard index. 0.1 is a chosen
   cut; the 5%/68% split is a function of it.

**A separate error in the same sentence.** The spec says "Level-1 self-maintaining
organisations turn out to be *more* frequent and robust than originally reported." Two
halves, and only one is what the paper says.

- *More frequent*: yes. Abstract, verbatim: "complex, stable organizations emerge more
  frequently than previously expected".
- *More robust*: the abstract's robustness claim is specifically "robust against collapse
  into **trivial fixed-points**", and the experiment behind it (Figure 2D) is on **L0**
  organisations perturbed by replacing a fraction of expressions with the identity
  function. For **L1**, the paper's finding is the opposite of a blanket robustness claim.
  Verbatim: "These results show that **the stability of 'L1 organizations' varies widely**,
  both through time and in robustness to external perturbations. In some cases (as in
  organization 2) they are highly robust, while in others (organization 1), they are highly
  sensitive." Of organisation 1: "suggesting that the original system was not truly stable."

  This is arguably *encouraging* for Borbax for a different reason than the spec gives. The
  authors say so: "This makes them interesting targets for a Darwinian process, because we
  know, e.g., that different organizations will have different capacities to respond to
  selective conditions."

**Fix — replace spec §2.1 line 59 with:**

> **A warning from this literature.** "Return to AlChemy" (Mathis, Patel, Weimer & Forrest,
> *Chaos* 34:093142, 2024, DOI 10.1063/5.0207358) re-ran the founding experiment of the
> field and partially falsified it. Complex stable organisations emerge *more* frequently
> than originally reported, and L0 organisations are strikingly robust against collapse
> into the trivial fixed point — encouraging. L1 stability, though, **varies widely**: some
> organisations are highly robust to perturbation and others turn out not to have been
> stable at all. And level 2 does not replicate. Combining two L1 organisations from
> *different* runs (N = 455 pairs), the outcome was **mutual destruction ~68%, one input
> dominating ~27%, coexistence ~5%**. The paper explicitly did not attempt the other route
> Fontana & Buss used — finding two separable organisations within a single run — so the
> non-replication covers one of the two original routes, not both.
>
> The implication for Borbax is unchanged and is if anything sharper. **We must not assume
> protocells will spontaneously compose into higher-order structures.** The single
> most-cited claim in artificial chemistry is, on the one route that has been re-tested,
> the field's weakest result. V2 onward must treat composition as an open research
> question.

Also update `.claude/agents/alife-researcher.md` lines 44–45 — that file currently teaches
the wrong numbers to every future invocation of this agent, which is how a bad figure
becomes project folklore.

**Test that would catch it:** none — this is a documentation defect. The mechanism that
should have caught it is a rule that quoted figures carry the table or figure number they
came from. Recommend: any numeric claim attributed to a paper in the spec cites the
specific table/figure, so a later reader can check it in one step.

---

## High

### alr-102 · `mean_cumulative` is the statistic Channon 2003 exists to replace

**Where:** plan `3004–3005`.
**Claim:** `Ā(t) = A/D` is unbounded for free in a system that retains any single component
indefinitely; Channon's correction is to use the **median**, and a well-mixed beaker with
persistent feedstock species is the worst case for it.
**Confidence:** verified from source.

Channon 2003, verbatim:

> "My other criticism of the test was in its use of mean activity when looking for unbounded
> activity growth, especially when classifying a system as belonging to class 3b. **When
> diversity is bounded, retention (forever) of a single component results in unbounded mean
> activity.** The test should not be so influenced by such components, and should rather
> look for trends in typical components. So it is **median activity, not mean activity**,
> that should be measured, and required to be unbounded for a system to be classified as
> within class 3b."

And MODES (Dolson et al. 2019) endorses it, verbatim: "the presence of a single component
under stabilizing selection will trivially cause the mean evolutionary activity to increase
indefinitely; such a component will sit in the population, increasing the population's
activity counter despite being quickly lost from the shadow population. … Channon suggested
that we should look for unbounded growth in median (rather than mean) per-component
evolutionary activity. This adjustment is a drastic improvement."

**Why it matters for Borbax specifically.** Channon's pathological case is "retention of a
single component forever with bounded diversity". A V0 beaker is *made of* that case: the
food set is by definition maintained, feedstock species persist for the whole run, and
diversity in a single well-mixed volume is bounded by the interner. So `mean_cumulative`
will climb monotonically in every beaker run, including a run of nothing but inert
feedstock, and the plan's own docstring says class 3b is "`cumulative` unbounded with
`diversity` bounded — visible only if both are computed". Both *will* be exactly that, for
free. The metric would report class-3b open-endedness on a dead universe.

**Fix.** Compute both, report the median as the classifier, and say why in the comment:

```rust
/// Ā(t) = A_cum/D — Bedau, Snyder & Packard's original statistic. **Kept for
/// comparability with the older literature and NOT used to classify.**
/// Channon (2003) showed it is unbounded whenever a single component is
/// retained indefinitely with bounded diversity, which is the default state
/// of a well-mixed beaker with a maintained food set (spec §9.3).
pub mean_cumulative: f64,
/// Ã(t) = median over extant components of a_i(t). This is the statistic the
/// class-3b verdict is read from (Channon 2003, eq. 14).
pub median_cumulative: f64,
```

```rust
/// Median over extant components. Deterministic tie-break: the input is
/// sorted by (activity, SpeciesId) so equal activities order by id, and the
/// even-length case averages the two central elements in index order.
fn median_extant(per_species: &[f64], present: &[SpeciesId]) -> f64 {
    let mut v: Vec<(f64, u32)> = present
        .iter()
        .map(|s| (per_species[s.0 as usize], s.0))
        .filter(|(a, _)| *a > 0.0)
        .collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    match v.len() {
        0 => 0.0,
        n if n % 2 == 1 => v[n / 2].0,
        n => (v[n / 2 - 1].0 + v[n / 2].0) / 2.0,
    }
}
```

**Test that would catch it:**

```rust
/// Channon (2003): retention of a single component forever, with bounded
/// diversity, drives MEAN cumulative activity to infinity while the typical
/// component is going nowhere. A beaker with a maintained food set is exactly
/// this case, so the mean must never be the class-3b classifier.
#[test]
fn a_single_immortal_species_does_not_look_open_ended() {
    let mut s = ActivityStats::with_capacity(4);
    // species 0 is immortal feedstock; 1..4 churn every snapshot
    for t in 0..2_000 {
        let churn = SpeciesId(1 + (t % 3) as u32);
        s.snapshot(&[SpeciesId(0), churn], &counts_for(&[SpeciesId(0), churn]));
    }
    assert!(
        s.mean_cumulative > 100.0,
        "sanity: the mean is supposed to blow up here — that is the point"
    );
    assert!(
        s.median_cumulative < 10.0,
        "median activity {} also grew: the classifier is inheriting the \
         immortal-component artefact Channon (2003) identified",
        s.median_cumulative
    );
}
```

---

### alr-103 · `diversity` is defined wrongly

**Where:** plan `3002–3003`, doc comment "`D(t)`: count of species above the activity
threshold."
**Claim:** Channon 2003 eq. (3) is `D(t) = #{i : aᵢ(t) > 0}` — the count of components with
*any* activity, i.e. the count of extant components. The threshold `[a₀, a₁]` belongs to
`A_new` (eq. 6), not to `D`.
**Confidence:** verified from source.

**Why it matters.** `D` is the denominator of both `Ā_cum` and `A_new`. Channon is explicit
about which `D` goes in the denominator, verbatim: "`D_R` (not `D_N`) is the relevant value
to use when calculating `Ā^N_cum`, `Ã^N_cum` and `A^N_new`, because `D_R` is the number of
components that contribute to those statistics." Putting a threshold in `D` shrinks the
denominator and inflates both derived statistics, in a way that grows as the threshold
bites harder — which happens exactly as the run matures. The bias is in the direction of
declaring open-endedness.

**Fix:**

```rust
/// D(t) = #{i : a_i(t) > 0}  — Channon (2003) eq. (3). The count of components
/// with non-zero activity, i.e. of extant components. **No threshold here.**
/// The [a0, a1] band belongs to A_new (eq. 6). D is the denominator of both
/// Ā_cum and A_new, so a threshold applied here inflates both.
pub diversity: f64,
```

**Test:**

```rust
#[test]
fn diversity_counts_every_extant_component_not_only_active_ones() {
    let mut s = ActivityStats::with_capacity(4);
    s.set_activity_band(50.0, 500.0);          // a0, a1
    // one long-lived species well above a0, three that just appeared
    for _ in 0..100 { s.snapshot(&[SpeciesId(0)], &c(&[0])); }
    s.snapshot(&[SpeciesId(0), SpeciesId(1), SpeciesId(2), SpeciesId(3)], &c(&[0,1,2,3]));
    assert_eq!(s.diversity, 4.0, "D applied the A_new band (Channon 2003 eq. 3 vs 6)");
}
```

---

### alr-104 · `new_activity` is defined wrongly

**Where:** plan `3006–3007`, "`A_new(t)`: activity of species newly crossing the threshold."
**Claim:** Channon 2003 eq. (6) is
`A_new(t) = (1/D(t)) · Σ_{i : aᵢ(t) ∈ [a₀, a₁]} aᵢ(t)`.
Three differences from the plan: it is a **mean** (divided by `D`), the window has an
**upper** bound `a₁`, and membership is by *current* activity lying in the band, not by a
*crossing event*.
**Confidence:** verified from source.

Channon on why `a₁` exists, verbatim: "For `A_new` to be a good measure of new activity, the
range `[a₀, a₁]` should be chosen such that component activities within it can be considered
both adaptively significant **and not amongst the highest**."

**Why it matters.** Without `a₁`, `A_new` is contaminated by the long-lived high-activity
components — the same ones that break `Ā_cum` in `alr-102`. Without the `1/D`
normalisation, `A_new` scales with population size rather than being per-component, so it
is not comparable between a shadow and a focal run with different diversities, which is the
one comparison the plan needs it for. Without the "current activity in band" semantics, a
component that crosses `a₀` and then keeps growing is counted once and never again — which
happens to be MODES's *novelty* metric, a different (also useful) statistic. Verbatim from
MODES: "Our novelty metric is functionally equivalent to `A_new` in evolutionary activity
statistics" — but MODES also says "Once a component has been counted as novel, however, it
is part of the permanent history and will never be counted in the novelty metric again."
The two are equivalent in spirit, not in formula. Pick one and cite the one you picked.

**Fix:**

```rust
/// A_new(t) = (1/D(t)) · Σ_{i : a_i(t) ∈ [a0, a1]} a_i(t)   — Channon (2003) eq. (6).
///
/// A **mean over the whole population**, not a sum, so it is comparable
/// between the focal run and the shadow when their diversities differ — which
/// is the only reason this statistic is computed.
///
/// The band is closed at BOTH ends. a1 excludes the long-lived high-activity
/// components: Channon requires activities in the band to be "adaptively
/// significant and not amongst the highest". Both bounds come from the shadow
/// (see `shadow::activity_band`), never from a constant.
fn new_activity(per_species: &[f64], present: &[SpeciesId], band: (f64, f64), d: f64) -> f64 {
    if d == 0.0 {
        return 0.0;
    }
    let mut ids: Vec<u32> = present.iter().map(|s| s.0).collect();
    ids.sort_unstable();                       // fixed accumulation order
    let mut acc = 0.0;
    for id in ids {
        let a = per_species[id as usize];
        if a >= band.0 && a <= band.1 {
            acc += a;
        }
    }
    acc / d
}
```

**Test:**

```rust
/// Channon (2003) eq. (6): the band is closed at both ends and the sum is
/// normalised by D. A component whose activity has run past a1 must drop out.
#[test]
fn new_activity_excludes_the_highest_activity_components() {
    let per = vec![0.0, 60.0, 120.0, 5_000.0];   // ids 1,2 in band; id 3 far above
    let present = vec![SpeciesId(1), SpeciesId(2), SpeciesId(3)];
    let a = new_activity(&per, &present, (50.0, 500.0), 3.0);
    assert!((a - (60.0 + 120.0) / 3.0).abs() < 1e-12,
            "A_new = {a}: either a1 was ignored or the 1/D normalisation is missing");
}
```

---

### alr-105 · Task 20b Steps 1 and 3 are ordered backwards

**Where:** plan Step 1 at `2983`, Step 3 at `3036`.
**Claim:** the activity band `[a₀, a₁]` in `A_new` is *defined* as the level the shadow says
is adaptively significant. Step 1 cannot be completed before Step 3.
**Confidence:** verified from source.

Channon 2003, verbatim:

> "For artificial systems, a 'shadow' is run, mirroring the real run in every detail except
> that whenever selection (artificial or natural) operates in the real system, random
> selection is employed in the shadow. **The statistics from this shadow can then be used to
> determine a₀** and levels of total and mean activity that can be considered adaptively
> significant."

And MODES on what activity *means*, verbatim: "evolutionary activity has been measured as
the length of time that components exist in the population **beyond what would be expected
in the absence of selection**."

**Why it matters.** The failure is not that the numbers are wrong; it is that as written an
implementer completing Step 1 in isolation has no source for `a₀` and will hard-code a
constant. Hard-coding the threshold is the practice this entire method was written to
replace, and once a constant is in the code and the tests pass, nothing later forces it
out. The spec already knows this — §15.3: "Without this, activity statistics are
uninterpretable" — but Task 20b's step order silently permits it.

There is also a second, subtler dependency the plan does not mention. Channon's headline
contribution in the 2003 paper is *shadow resetting*, verbatim: "immediately after each
snapshot (when an entry is made in the component existence record), the shadow run has its
components reset to those of the real run. This allows us to compare inter-snapshot changes
in activity in the real run with the changes we would expect from random selection." His
reason: "the components that exist in the real population at any one time … are almost
certainly more densely clustered than those in the shadow. So the mutation of a real
component is more likely to produce another high-activity component … Once the real and
shadow populations have been allowed to evolve, we are no longer comparing the real run with
a true shadow."

Borbax's §22.3 already spawns shadows on demand from keyframes rather than running them
continuously. That is *closer* to Channon's resetting method than to the drifting-shadow
method he criticises, which is good and worth saying in the module header — but the reset
cadence should be a stated parameter (keyframe spacing), not an accident of §22.3.

**Fix — reorder and make the dependency a type:**

```rust
/// The activity band, which only a shadow can produce.
///
/// Channon (2003): "The statistics from this shadow can then be used to
/// determine a0 and levels of total and mean activity that can be considered
/// adaptively significant." There is no default and no `Default` impl, on
/// purpose: `ActivityStats` cannot be constructed without one, so the
/// shadow cannot be skipped and a constant cannot be substituted for it.
#[derive(Debug, Clone, Copy)]
pub struct ActivityBand {
    pub a0: f64,
    pub a1: f64,
}

impl ActivityBand {
    /// a0 = the quantile of the SHADOW's per-component activity distribution
    /// above which a real component's activity is not explained by drift.
    /// a1 = the upper cut excluding "amongst the highest" (Channon 2003 eq. 6).
    #[must_use]
    pub fn from_shadow(shadow: &[f64], q_lo: f64, q_hi: f64) -> Self {
        let mut v: Vec<f64> = shadow.iter().copied().filter(|a| *a > 0.0).collect();
        v.sort_by(f64::total_cmp);
        Self { a0: quantile_sorted(&v, q_lo), a1: quantile_sorted(&v, q_hi) }
    }
}

impl ActivityStats {
    /// No `Default`. The band is a required argument.
    #[must_use]
    pub fn new(capacity: usize, band: ActivityBand) -> Self { /* ... */ }
}
```

Then swap the plan's step order: **Step 1 becomes the shadow, Step 2 the activity band,
Step 3 the statistics.**

**Test:**

```rust
/// The band must be derived, not chosen. This test fails to COMPILE if
/// ActivityStats gains a Default or a band-free constructor.
#[test]
fn activity_stats_cannot_be_built_without_a_shadow_derived_band() {
    let shadow = run_shadow_activity(seed(9), 10_000);
    let band = ActivityBand::from_shadow(&shadow, 0.95, 0.999);
    let _s = ActivityStats::new(64, band);
    // There is deliberately no ActivityStats::new(64) and no Default.
}

/// A shadow-derived a0 is not a round number. If it comes out as exactly the
/// value someone would have typed, it was typed.
#[test]
fn the_derived_threshold_is_not_a_hardcoded_constant() {
    let a = ActivityBand::from_shadow(&run_shadow_activity(seed(1), 10_000), 0.95, 0.999);
    let b = ActivityBand::from_shadow(&run_shadow_activity(seed(2), 10_000), 0.95, 0.999);
    assert_ne!(a.a0.to_bits(), b.a0.to_bits(),
               "a0 identical across two shadows — it is not being derived from them");
}
```

**Also fix Step 5b (plan `3155–3160`) and exit criterion 9 (spec `1147`)** — see `alr-121`.

---

### alr-106 · The null-trajectory test fails as written; the bias is computed and there is an exact fix

**Where:** plan `3017–3034` (`shape_novelty` doc comment), `3034` ("must score ~0"),
`3080–3082`.
**Claim:** subsampling histograms to a constant `N` removes the sample-size dependence of
the plug-in distance estimator but **not** the dependence on `K`, the number of occupied
bins. Since `K` grows with species count over a run, novelty rises with diversity even when
the underlying distribution is unchanged — the exact activity/novelty conflation §15.2
exists to prevent.
**Confidence:** computed.

**The measurement.** Two independent samples of size `N` drawn from the *same* uniform
distribution over `K` bins. True distance is exactly zero. Plug-in estimate, mean of 300
replicates:

```
     N     K         TV   Hellinger         JS   sqrt(K/N)
   200     5    0.08007     0.06681    0.00502     0.1581
   200    50    0.27497     0.27142    0.06754     0.5000
   200   200    0.52213     0.63285    0.28773     1.0000
  1000    10    0.05335     0.04601    0.00223     0.1000
  1000   100    0.17655     0.16035    0.02545     0.3162
  1000   200    0.24779     0.23797    0.05336     0.4472
```

At fixed `N = 1000`, going from 10 occupied shape bins to 100 moves "novelty" from 0.053 to
0.177 — a **3.3× rise with zero change in the distribution**. The bias tracks `√(K/N)`.

**Why it matters.** In a beaker, `K` is the count of occupied shape bins, which is the
species count in disguise. The plan's own docstring says the subsampling exists so that
"novelty [does not] rise merely because species count rose" — and the subsampling does not
achieve that. Worse, Task 20b Step 6 asks for the null-trajectory score to be recorded in
the commit message as validation. It will not read ~0, the implementer will find a
threshold that makes it look like it does, and the metric ships calibrated to noise.

**Fix — an exactly unbiased estimator of squared L2 distance.** For multinomial counts
`a ~ Mult(N, p)` and `b ~ Mult(M, q)`, `E[Σ(p̂ᵢ − q̂ᵢ)²] = Σ(pᵢ − qᵢ)² + (1 − Σpᵢ²)/N +
(1 − Σqᵢ²)/M`, and `Σ aᵢ(aᵢ − 1) / (N(N−1))` is unbiased for `Σpᵢ²`. So:

```
D̂² = Σᵢ (aᵢ/N − bᵢ/M)²  −  (1 − Σᵢ aᵢ(aᵢ−1)/(N(N−1)))/N  −  (1 − Σᵢ bᵢ(bᵢ−1)/(M(M−1)))/M
```

Verified numerically (400 reps per cell, plug-in vs unbiased against the analytic truth):

```
    K      N   true L2^2      plug-in     unbiased
  200    500   0.000e+00    3.993e-03    1.345e-05
  200    500   1.250e-03    5.238e-03    1.261e-03
 1000    500   0.000e+00    4.007e-03    1.112e-05
 1000    500   2.500e-04    4.258e-03    2.625e-04
 1000   5000   2.500e-04    6.477e-04    2.482e-04
```

At `K = 1000, N = 500` the plug-in reads `4.0e-3` for two samples of the *same*
distribution — **16× larger than the true signal** of a genuine 50% redistribution of mass
(`2.5e-4`). The corrected estimator recovers the truth to three significant figures at
every `K` and `N` tested, and reads within `1.3e-5` of zero on the null.

```rust
/// Unbiased estimator of the squared L2 distance between the two multinomial
/// distributions two histograms are samples of.
///
/// The plug-in estimator Σ(p̂-q̂)² is positively biased by
/// (1-Σp²)/N + (1-Σq²)/M, and that bias grows with the number of occupied
/// bins. Measured: two samples of size 500 from the SAME uniform distribution
/// over 1000 bins give a plug-in distance of 4.0e-3 — sixteen times the true
/// distance produced by a genuine 50% redistribution of mass. Since occupied
/// bins track species count, an uncorrected novelty metric rises because
/// diversity rose, which is precisely the activity/novelty conflation §15.2
/// exists to prevent.
///
/// The correction uses E[Σ a_i(a_i-1)] = N(N-1)·Σp_i², so it costs one extra
/// pass and no tuning. The result may be slightly negative on the null; that
/// is correct behaviour for an unbiased estimator and must not be clamped —
/// clamping reintroduces a positive bias.
#[must_use]
pub fn l2_sq_unbiased(a: &[u64], b: &[u64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let n = a.iter().sum::<u64>() as f64;
    let m = b.iter().sum::<u64>() as f64;
    if n < 2.0 || m < 2.0 {
        return 0.0;
    }
    // Fixed accumulation order: index order. Load-bearing (§13.1).
    let mut plug = 0.0;
    let mut pa = 0.0;
    let mut pb = 0.0;
    for i in 0..a.len() {
        let (x, y) = (a[i] as f64, b[i] as f64);
        let d = x / n - y / m;
        plug += d * d;
        pa += x * (x - 1.0);
        pb += y * (y - 1.0);
    }
    let pur_a = pa / (n * (n - 1.0));
    let pur_b = pb / (m * (m - 1.0));
    plug - (1.0 - pur_a) / n - (1.0 - pur_b) / m
}
```

**Test that would catch it:**

```rust
/// Task 20b Step 6 requires the novelty metric to read ~0 on a null
/// trajectory. It must do so INDEPENDENTLY OF K, because K tracks species
/// count and a K-dependent floor is the activity/novelty conflation the
/// metric exists to prevent.
#[test]
fn null_trajectory_novelty_is_flat_in_the_number_of_occupied_bins() {
    let mut r = Stream::new(3, Domain::Metrics, 0);
    let mut scores = Vec::new();
    for k in [10usize, 100, 1_000] {
        let mut acc = 0.0;
        for _ in 0..200 {
            let a = multinomial_uniform(&mut r, 500, k);
            let b = multinomial_uniform(&mut r, 500, k);
            acc += l2_sq_unbiased(&a, &b);
        }
        scores.push(acc / 200.0);
    }
    for (k, s) in [10, 100, 1000].iter().zip(&scores) {
        assert!(s.abs() < 1e-4, "null novelty at K={k} reads {s}, not ~0");
    }
    // and the K-dependence itself, which is the actual defect:
    let spread = scores.iter().cloned().fold(f64::MIN, f64::max)
        - scores.iter().cloned().fold(f64::MAX, f64::min);
    assert!(spread < 1e-4,
            "null novelty varies by {spread} across K=10..1000: the estimator \
             is plug-in, so novelty will rise as diversity rises");
}
```

---

### alr-107 · Unweighted RSS on increments cannot separate bounded from unbounded — measured

**Where:** plan `3063–3069`; spec §15.4 line 875.
**Claim:** the plan is right that increments must be fitted rather than cumulative totals,
and right that all models must be fitted on the same scale by the same method. But it then
prescribes homoscedastic least squares (`rss`) on decaying counts, and **that fit is
dominated by the head of the curve, which carries no information about boundedness.**
**Confidence:** computed, on noiseless data.

**Measurement 1 — where the RSS lives.** For a bounded (hyperbolic) increment series over
2000 windows, perturbed by a constant *relative* 10% everywhere:

```
first   1% of windows carry:  87.02% of Gaussian RSS    49.88% of Poisson deviance
first   5% of windows carry:  99.51% of Gaussian RSS    83.84% of Poisson deviance
first  10% of windows carry:  99.92% of Gaussian RSS    91.63% of Poisson deviance
```

The plateau verdict on a 2000-window run is decided by the first 20 windows.

**Measurement 2 — what that does to the answer.** Fit a single power law `g = a·t^c` to
noiseless increments from each generating model. The cumulative diverges iff `c > −1`:

```
t in [1,2000]:
  hyperbolic b=20 (BOUNDED, true c=-2)   unweighted c= -0.718   log-space c= -1.772
  power law b=0.5 (UNBOUNDED, c=-0.5)    unweighted c= -0.500   log-space c= -0.500
  exp sat  b=0.05 (BOUNDED)              unweighted c= -0.715   log-space c=-25.356
```

Under unweighted least squares the two **bounded** generators come out at `c = −0.718` and
`c = −0.715`, versus `−0.500` for the genuinely unbounded one — and all three are declared
unbounded (`c > −1`). That is on data with **no noise at all**. Under a constant-relative-
error (equivalently, Poisson-deviance) fit the same data separates cleanly.

**Why it matters.** This is the boundedness illusion arriving through the fitting procedure
rather than through the eye — which is worse, because it wears the authority of a model
comparison. Borbax would report "power law beats saturating, therefore unbounded" on a run
that is genuinely saturating.

**Fix.** Novelty increments are **counts**, so the error model is Poisson, not Gaussian.
Fit every model by Poisson deviance and score every model the same way. This satisfies the
spec's own "same scale by the same method" rule better than RSS does, and it does not break
at zero counts the way a log-space fit would.

```rust
/// Poisson deviance, the residual for a count process.
///
/// Novelty increments are counts, so Var ≈ mean and homoscedastic least
/// squares is the wrong error model. Measured on noiseless data over 2000
/// windows: unweighted RSS puts 87% of its weight in the first 1% of windows,
/// and fits a bounded hyperbolic series and a bounded exponential series to
/// power-law exponents of -0.718 and -0.715 against -0.500 for a genuinely
/// unbounded one — declaring all three unbounded. The boundedness illusion
/// (spec §2.7, §15.4) arriving through the estimator.
///
/// D = 2 Σ [ y ln(y/mu) - (y - mu) ], with the y=0 term taken as 2*mu.
/// levenberg-marquardt minimises a sum of squares, so residuals are handed to
/// it as signed square roots of the per-point deviance.
#[must_use]
pub fn poisson_residual(y: f64, mu: f64) -> f64 {
    let mu = mu.max(1e-12);
    let d = if y > 0.0 {
        2.0 * (y * det_math::ln(y / mu) - (y - mu))
    } else {
        2.0 * mu
    };
    // signed sqrt: LM squares it back, recovering the deviance
    d.max(0.0).sqrt() * if y >= mu { 1.0 } else { -1.0 }
}

/// AIC for a Poisson fit: -2·loglik + 2K, K = number of regression parameters
/// only. There is no sigma^2 to estimate, so — unlike the Gaussian case
/// (see `alr-111`) — K does not gain one.
#[must_use]
pub fn poisson_aic(nll: f64, k_reg: usize) -> f64 {
    2.0 * nll + 2.0 * k_reg as f64
}
```

**Test that would catch it:**

```rust
/// Synthetic data with a KNOWN answer, and no noise, so a failure is the
/// estimator and nothing else. Task 21 Step 5b already asks for exactly this
/// ("distinguishes a saturating series from a power-law one on synthetic data
/// before it is trusted on real") — this is that test with the hard case in it.
#[test]
fn the_fitter_calls_a_hyperbolic_series_bounded() {
    let t: Vec<f64> = (1..=2000).map(f64::from).collect();
    // hyperbolic cumulative a*t/(t+b) -> increments a*b/(t+b)^2. BOUNDED.
    let y: Vec<f64> = t.iter().map(|&t| 100.0 * 20.0 / (t + 20.0).powi(2)).collect();
    let v = fit_models(&t, &y);
    assert!(
        v.best_is_bounded(),
        "declared {} on a hyperbolic series whose cumulative provably \
         converges — the fit is dominated by the head of the curve",
        v.best_name()
    );
}

#[test]
fn the_fitter_calls_a_power_law_unbounded() {
    let t: Vec<f64> = (1..=2000).map(f64::from).collect();
    let y: Vec<f64> = t.iter().map(|&t| 25.0 * det_math::powf(t, -0.5)).collect();
    let v = fit_models(&t, &y);
    assert!(!v.best_is_bounded(), "declared bounded on a divergent power law");
}
```

---

### alr-108 · The bounded competitor should be hyperbolic, not exponential-saturating

**Where:** plan `3063`; spec §15.4 line 870 ("saturating, linear, power-law").
**Claim:** the bounded model in the source Borbax cites is hyperbolic. Exponential
saturation is a straw man: its increments decay geometrically, so rejecting it is trivially
easy and proves nothing about boundedness.
**Confidence:** the *mathematics* is computed (below). The *functional forms in Wiser,
Ribeck & Lenski 2013* I **could not verify** — science.org returned 403. What I have is the
abstract and secondary descriptions, which agree that the two fitted models were
"hyperbolic" and "power law", that both decelerate, and that "only the hyperbolic model has
an upper limit or asymptote." I have **not** read the equations and will not reproduce them
from memory.

**The mathematics, which does not depend on the paper.** For cumulative `A(t)`:

| model | cumulative | increments | bounded? |
|---|---|---|---|
| exponential saturating | `a(1 − e^{−bt})` | `ab·e^{−bt}` | yes, geometrically |
| hyperbolic | `a·t/(t+b)` | `ab/(t+b)²` ~ `t^{−2}` | yes, as `1/t` |
| power law | `a·t^b`, `b>0` | `ab·t^{b−1}` ~ `t^{c}`, `c>−1` | no |

A hyperbolic bounded model **is a power law in increments**, with exponent `−2`. So the
boundedness question reduces to a single sharp one: *is the increment exponent above or
below `−1`?* An exponential-saturating model cannot express that question — it is not in
the power-law family at all, and a fit will reject it for reasons that have nothing to do
with boundedness. From the measurement in `alr-107`, over `t ∈ [50, 2000]` under a
correctly-weighted fit the exponential series comes out at `c = −32.5` while the hyperbolic
comes out at `c = −1.90`, close to its true `−2`. One of those is a real contest.

**Fix — four models, and the verdict is about `c`, not about which model wins:**

```rust
/// Candidate models for the novelty INCREMENT series.
///
/// The bounded competitor in the source Borbax cites (Wiser, Ribeck & Lenski
/// 2013, Science 342:1364, DOI 10.1126/science.1243357) is **hyperbolic**, not
/// exponential-saturating, and that choice is the whole point. A hyperbolic
/// cumulative a·t/(t+b) has increments ~ t^-2, which lives inside the
/// power-law family; distinguishing it from an unbounded power law is the
/// discrimination that actually tests the boundedness illusion. An
/// exponential-saturating model decays geometrically and is rejected for
/// reasons unrelated to boundedness — measured exponent -32.5 against -1.90
/// for hyperbolic over the same window. Keep it, but as a control, not as
/// the bounded case.
pub enum IncrementModel {
    /// g = a                       cumulative a·t          UNBOUNDED
    Constant,
    /// g = a·b·t^(b-1)             cumulative a·t^b        UNBOUNDED iff b > 0
    PowerLaw,
    /// g = a·b/(t+b)^2             cumulative a·t/(t+b)    BOUNDED  (asymptote a)
    Hyperbolic,
    /// g = a·b·exp(-b·t)           cumulative a(1-e^-bt)   BOUNDED  (control)
    ExpSaturating,
}

impl IncrementModel {
    /// A model's verdict is structural, not a fitted quantity — except for
    /// PowerLaw, where boundedness is exactly `b <= 0` (equivalently, fitted
    /// increment exponent <= -1).
    #[must_use]
    pub fn is_bounded(self, theta: &[f64]) -> bool {
        match self {
            Self::Constant => false,
            Self::PowerLaw => theta[1] <= 0.0,
            Self::Hyperbolic | Self::ExpSaturating => true,
        }
    }
}
```

**Test:** the two in `alr-107`, plus a control that the exponential model is *not* what
carries the verdict:

```rust
#[test]
fn removing_the_exponential_model_does_not_change_the_verdict() {
    let (t, y) = hyperbolic_series(2000);
    let with = fit_models(&t, &y);
    let without = fit_models_excluding(&t, &y, IncrementModel::ExpSaturating);
    assert_eq!(with.best_is_bounded(), without.best_is_bounded(),
               "the boundedness verdict is being carried by the straw-man model");
}
```

---

### alr-109 · `rate_prefactor` is never read; two battery tests cannot pass; volume `V` does not exist

**Where:** `rate()` at plan `1896–1902`; battery tests at plan `2886` and `2893`;
`rate_prefactor` generated at `2026-07-26-borbax-v0.md:2049`.
**Claim:** `rate()` computes `arrhenius * a * b * catalysis` and never multiplies by
`u.consts.rate_prefactor`. `rate_prefactor` is written by `Universe::generate` and read
nowhere in either plan file. Grepped both plans: the only three occurrences outside its
declaration are the two battery tests that set it.
**Confidence:** verified by grep.

```
$ grep -n "rate_prefactor" docs/superpowers/plans/*.md
v0.md:2000:    pub rate_prefactor: f64,
v0.md:2049:            rate_prefactor: rng.next_f64_range(1e3, 1e5),
v0-chemistry.md:2886:        u.consts.rate_prefactor = 1e-30; // nothing will ever react
v0-chemistry.md:2893:        u.consts.rate_prefactor = 1e30;
```

**Why it matters — three ways.**

1. `it_rejects_an_inert_universe` and `it_rejects_a_universe_that_burns` set a value that
   nothing reads. Both will run an ordinary universe and assert `!passes()`. They will fail
   or, worse, pass by luck of the seed and encode nothing.
2. Spec §9.1 line 446 gives the rate law as `A · exp(−Ea/T) · [X] · [Y] · catalysis_factor`.
   `A` is missing from the implementation of the law it cites.
3. In Gillespie's stochastic formulation the deterministic rate constant `k` and the
   stochastic one `c` are related by `c = k/V` for a bimolecular reaction and `c = 2k/V`
   for a same-species one. Borbax has no `V` anywhere. That is *defensible* for a
   well-mixed V0 — set `V = 1` and absorb it into the prefactor — but it must be a stated
   choice, because §2.7's Moreno/Ofria countermeasure ("treat region size as a complexity
   budget; if a run stalls, suspect this before touching parameters") is **not testable in
   V0 without a volume knob**. A universe rejected by the battery cannot be distinguished
   from a beaker too small to nucleate anything.

**Fix:**

```rust
/// Arrhenius-style rate (spec §9.1):  A · exp(-Ea/T) · [X] · [Y] · catalysis.
///
/// `prefactor` is `A` = `Universe::consts.rate_prefactor`. It was previously
/// absent, which made two battery tests (inert / burning universe)
/// unfalsifiable, since both manipulate a constant nothing read.
///
/// **On volume.** Gillespie's stochastic rate constant is c = k/V for a
/// bimolecular reaction and 2k/V for a same-species one (Gillespie 1977,
/// J. Phys. Chem. 81:2340). V0 fixes V = 1 by construction — one well-mixed
/// beaker, counts rather than concentrations — so V is absorbed into
/// `prefactor` and never appears. This is a decision, not an omission:
/// it means V0 cannot vary reaction volume, and therefore cannot separate
/// "this universe is dead" from "this beaker is too small" (spec §2.7,
/// complexity carrying capacity). V1 must reintroduce V explicitly when
/// compartments arrive, because a compartment IS a different volume.
#[must_use]
pub fn rate(r: &Reaction, t: Thermal, prefactor: f64, counts: &[f64], catalysis: f64) -> f64 {
    let temp = t.get().max(1e-6);
    let arrhenius = det_math::exp(-r.activation.get() / temp);
    prefactor * arrhenius * combinatorial_factor(r, counts) * catalysis
}
```

**Test:**

```rust
/// The two battery tests set rate_prefactor and assert on the outcome. If the
/// constant is not on the rate path they assert nothing. This test fails
/// loudly rather than silently.
#[test]
fn rate_scales_linearly_with_the_universe_prefactor() {
    let r = fixture_reaction();
    let c = [10.0, 10.0, 0.0, 0.0];
    let lo = rate(&r, Thermal(300.0), 1.0,   &c, 1.0);
    let hi = rate(&r, Thermal(300.0), 1e6,   &c, 1.0);
    assert!((hi / lo - 1e6).abs() < 1e-6,
            "rate_prefactor is not on the rate path: §9.1's `A` is missing and \
             the inert/burning battery tests are unfalsifiable");
}
```

---

### alr-110 · `[SpeciesId; 2]` has no arity, so unimolecular reactions have no correct encoding

**Where:** `Reaction.reactants` at plan `1871`; `missing` at plan `2483`; propensity rule at
plan `2721`; `rate()` at `1901`.
**Claim:** `ReactionKind::Rearrange` is unimolecular and `Cleave` is 1→2, but `reactants` is
a fixed `[SpeciesId; 2]` with no arity tag and no documented `NONE` convention. Every
available encoding is wrong somewhere.
**Confidence:** verified by execution.

I ran the three candidate encodings through the plan's own `find_raf` and `Incidence`:

```
P2 unimolecular encoded as duplicated reactant [A,A]  -> RAF found   (RAF ok)
P3 unimolecular encoded with a pad id INSIDE n_species -> None        (RAF broken)
P4 pad id ABOVE n_species:
     incidence entries for reaction 0 = 1
     missing (plan)  = [2]   <- unreachable, closure can never complete
     missing (fixed) = [1]   <- reachable
```

So:

- **Duplicate `[A, A]`** satisfies the RAF closure — but collides head-on with the plan's
  own propensity rule at line 2721. A reaction whose two reactant slots hold the same
  species is indistinguishable from a genuine `A + A → …`, so the beaker will score a
  *unimolecular* channel as `n(n−1)/2` instead of `n`. Computed error factor:

  ```
  n            n     n(n-1)/2    over-fast by
  2            2          1.0           0.5x
  10          10         45.0           4.5x
  100        100       4950.0          49.5x
  1000      1000     499500.0         499.5x
  ```

  The error is `(n−1)/2` — unbounded, and it is *fastest where counts are largest*, so a
  `Rearrange` channel on an abundant species runs hundreds of times too fast.
- **In-range pad** silently blocks the RAF forever (the pad is simply an absent reactant).
- **Out-of-range pad** breaks the closure permanently, because plan line 2483 computes
  `missing = r.reactants.len()` (always 2) while `Incidence::build` only writes entries for
  slots with `s.0 < n_species`. The two disagree and the counter never reaches zero. Worse,
  `n_species` at plan 2457–2463 is a max over reactant ids, so a sentinel of `u32::MAX`
  would allocate a 4-billion-entry `Vec`.

**Why it matters.** `Cleave` is the engine's primary decay path (§9.5) and §20 names decay
as the most sensitive parameter in the system. Task 15's `decay_channels` is correct — it is
linear in count — but *catalysed* cleavage and `Rearrange` go through `Reaction` and the
beaker's propensity rule. Getting them 50× too fast poisons the decay band search
(Task 21 Step 2, exit criterion 4) in a way that will read as "the universe burns".

**Fix — make arity explicit and derive everything from one predicate:**

```rust
/// Reactant multiset. Arity is explicit because a fixed [SpeciesId; 2] cannot
/// distinguish `A -> …` (propensity c·n) from `A + A -> …`
/// (propensity c·n(n-1)/2), and Rearrange and uncatalysed Cleave are both
/// unimolecular. Encoding a unimolecular reaction as [A, A] makes it run
/// (n-1)/2 times too fast — 49.5x at n=100 — and Cleave is the primary decay
/// path (spec §9.5), the most sensitive parameter in the system (§20).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reactants {
    One(SpeciesId),
    Two(SpeciesId, SpeciesId),
}

impl Reactants {
    /// The single source of truth for "which slots are real reactants".
    /// `Incidence::build` and the closure's `missing` counter must BOTH use
    /// this, or they disagree and the RAF closure silently never completes.
    #[must_use]
    pub fn as_slice(&self) -> &[SpeciesId] {
        match self {
            Self::One(a) => std::slice::from_ref(a),
            Self::Two(a, b) => std::slice::from_ref(a), // see note
        }
    }
    #[must_use]
    pub fn len(&self) -> usize {
        match self { Self::One(_) => 1, Self::Two(_, _) => 2 }
    }
    /// Gillespie 1977 combinatorial factor h_mu.
    #[must_use]
    pub fn propensity_factor(&self, counts: &[f64]) -> f64 {
        match *self {
            Self::One(a) => counts[a.0 as usize],
            Self::Two(a, b) if a == b => {
                let n = counts[a.0 as usize];
                n * (n - 1.0) / 2.0
            }
            Self::Two(a, b) => counts[a.0 as usize] * counts[b.0 as usize],
        }
    }
}
```

and in `find_raf`, derive `missing` from the same predicate rather than restating it:

```rust
// missing[ri] must equal the number of incidence entries actually written
// for ri, or the closure cannot complete. Deriving both from
// `Reactants::len()` makes disagreement impossible; restating it as
// `reactants.len()` on a fixed-arity array does not.
let missing: Vec<u32> = rxns.iter().map(|r| r.reactants.len() as u32).collect();
debug_assert_eq!(
    missing.iter().map(|&m| m as usize).sum::<usize>(),
    inc.total_entries(),
    "Incidence and the missing-reactant counter disagree: the closure will \
     never complete for at least one reaction"
);
```

**Tests:**

```rust
/// Gillespie (1977) h_mu forms. A unimolecular channel is linear in n; a
/// same-species bimolecular one is n(n-1)/2. At n=100 the difference is 49.5x.
#[test]
fn unimolecular_and_same_species_propensities_differ() {
    let counts = vec![0.0, 100.0];
    let uni = Reactants::One(SpeciesId(1)).propensity_factor(&counts);
    let bi  = Reactants::Two(SpeciesId(1), SpeciesId(1)).propensity_factor(&counts);
    assert_eq!(uni, 100.0);
    assert_eq!(bi, 4950.0);
    assert!((bi / uni - 49.5).abs() < 1e-9,
            "a unimolecular reaction encoded as [A,A] would run 49.5x too fast");
}

/// The incidence structure and the closure's counter must be derived from the
/// same predicate. This asserts the invariant directly rather than hoping.
#[test]
fn incidence_entries_match_the_missing_reactant_counters() {
    let rxns = vec![
        Reaction { reactants: Reactants::One(SpeciesId(0)), ..fixture() },
        Reaction { reactants: Reactants::Two(SpeciesId(0), SpeciesId(1)), ..fixture() },
        Reaction { reactants: Reactants::Two(SpeciesId(2), SpeciesId(2)), ..fixture() },
    ];
    let inc = Incidence::build(&rxns, 3);
    let total: usize = rxns.iter().map(|r| r.reactants.len()).sum();
    assert_eq!(inc.total_entries(), total,
               "closure would stall: some reaction's missing-counter can never reach 0");
}

/// A unimolecular reaction grounded in food must yield a RAF.
#[test]
fn a_unimolecular_reaction_can_close_a_raf() {
    let rxns = vec![Reaction {
        reactants: Reactants::One(SpeciesId(0)),
        products: [SpeciesId(1), SpeciesId(1)],
        catalyst: Some(SpeciesId(1)),
        ..fixture()
    }];
    let cat = vec![(SpeciesId(1), 0usize)];
    assert!(find_raf(&rxns, &cat, &food(&[0])).is_some(),
            "a unimolecular reaction was treated as bimolecular with a missing partner");
}
```

*Note: the `as_slice` sketch above has a deliberate hole in the `Two` arm — a slice cannot
be formed over two separately-stored fields. Store `Reactants::Two` as `[SpeciesId; 2]`
internally, or return a `SmallVec`/iterator. Flagging so the fix is not pasted as-is;
`rust-developer-expert` should own the final shape.*

---

## Medium

### alr-111 · AICc `K` must include `σ²` — computed

**Where:** plan `3069`.
**Claim:** for a least-squares fit, the number of estimated parameters `K` includes the
residual variance `σ̂²`. The plan's formula uses the regression-parameter count. The extra
term cancels in `AIC` and `BIC` differences between models of equal `K`, but **does not
cancel in `AICc`**, and the plan's model set does not have equal `K`: on increments, linear
has one parameter (a constant) while power-law, hyperbolic and saturating have two.
**Confidence:** computed.
**Citation:** Burnham & Anderson (2002), *Model Selection and Multimodel Inference*, 2nd ed.,
Springer, §2.2 — "the value of K … must include the estimated residual variance σ̂² as a
parameter."

Computed ΔAICc between a 2-regression-parameter and a 1-regression-parameter model at
identical RSS:

```
n=  10   plan=  3.2143   correct=  4.2857   error=+1.0714
n=  12   plan=  2.9333   correct=  3.6667   error=+0.7333
n=  20   plan=  2.4837   correct=  2.7941   error=+0.3105
n=  40   plan=  2.2306   correct=  2.3846   error=+0.1540
n= 100   plan=  2.0829   correct=  2.1263   error=+0.0434
```

The plan makes the two-parameter models look 1.07 AICc units better than they are at
`n = 10`. Burnham & Anderson's own rule of thumb is that `ΔAIC < 2` means substantial
support for the weaker model, so an error of 1.07 is inside the band where verdicts flip.
And `n = 10` is realistic: the plan holds out the tail (line 3077), so the fit sees half the
windows.

**Fix:**

```rust
/// Information criteria for a least-squares fit.
///
/// `k_reg` is the number of REGRESSION parameters. K = k_reg + 1, because the
/// residual variance sigma^2 is also estimated (Burnham & Anderson 2002 §2.2).
/// The +1 cancels in AIC and BIC differences between models of equal k_reg —
/// but NOT in AICc, and not between models of unequal k_reg, which is exactly
/// this model set (constant has one parameter, the rest have two).
/// Computed: at n=10 omitting it shifts ΔAICc by 1.07, inside the ΔAIC<2 band
/// where a verdict flips.
#[must_use]
pub fn info_criteria(n: usize, rss: f64, k_reg: usize) -> Option<InfoCriteria> {
    let k = (k_reg + 1) as f64;   // + sigma^2
    let n_f = n as f64;
    let denom = n_f - k - 1.0;
    if denom <= 0.0 {
        // AICc is undefined here and returns a *negative* correction that
        // would rank the most complex model best. Refuse rather than mislead.
        return None;
    }
    let aic = n_f * det_math::ln(rss / n_f) + 2.0 * k;
    Some(InfoCriteria {
        aic,
        aicc: aic + 2.0 * k * (k + 1.0) / denom,
        bic: n_f * det_math::ln(rss / n_f) + k * det_math::ln(n_f),
    })
}
```

**Test:**

```rust
/// Burnham & Anderson (2002) §2.2: K includes sigma^2. Checked against a
/// hand-computed value so a later "simplification" cannot silently drop it.
#[test]
fn aicc_counts_the_residual_variance_as_a_parameter() {
    let ic = info_criteria(10, 10.0, 2).unwrap();   // k_reg=2 -> K=3
    let k = 3.0;
    let expect_aic = 10.0 * (10.0f64 / 10.0).ln() + 2.0 * k;      // = 6.0
    let expect_aicc = expect_aic + 2.0 * k * (k + 1.0) / (10.0 - k - 1.0);
    assert!((ic.aicc - expect_aicc).abs() < 1e-12,
            "AICc {} != {expect_aicc}: K is probably k_reg, not k_reg+1", ic.aicc);
}

/// n - K - 1 <= 0 makes the AICc correction negative, which would rank the
/// most complex model best. Must refuse.
#[test]
fn aicc_refuses_when_the_sample_is_too_small() {
    assert!(info_criteria(4, 1.0, 2).is_none());
}
```

---

### alr-112 · The MODES persistence filter is a *lineage* filter and needs a phylogeny — but Channon's filter does not

**Where:** plan `3011–3013`.
**Claim:** the plan says "Add the **MODES persistence filter** alongside: discard components
not surviving a set interval before counting. It is the cheap, well-tested substitute for a
full shadow." Two errors: the filter is lineage-based, and the paper does not present it as
a substitute for a shadow.
**Confidence:** verified from source (MODES LaTeX).

Verbatim, on what the filter is:

> "we limit our analysis to those components whose **descendants** persist for a substantial
> number of generations. … We mark each organism with a **lineage ID** at a given time point
> A … The lineage IDs are passed on to **offspring** for the next t generations … At time
> point A+t, we determine which components from the population at A have descendants at
> A+t."

and on cost: "The next largest cost is imposed by needing to keep track of the **phylogeny**
over time."

Borbax V0 has no descent. Species arise by reaction, not by reproduction; there is no
parent-offspring relation to carry a lineage ID along. The filter is not implementable as
published.

Verbatim, on the relationship to a shadow — which is the opposite of the plan's framing:

> "All of the MODES metrics assume that some form of filtering has been applied before-hand
> … Here, we describe a persistence filter that we use in our experiments. **However, shadow
> runs are also a viable filter option**"
>
> "Whereas shadow runs filter out the effect of neutral processes, **the persistence filter
> does not entirely**. We view this reduced filtering primarily as an advantage — drift can
> be an important part of the evolutionary process — but there may also be situations where
> it is undesirable. Our metrics are unable to distinguish between class 1 and 2 dynamics or
> between class 3 and 4b dynamics"
>
> (Conclusions) "further investigation on the differences between using a shadow run as a
> filter and using the persistence filter described here would be worthwhile. Ultimately,
> **these two techniques capture sufficiently different information that it may be valuable
> to use each in turn.**"

They are alternatives with different blind spots, not a cheap substitute and an expensive
original. And Task 20b builds the shadow anyway (Step 3), so the "substitute" framing does
not even describe the plan's own scope.

**The better fix than dropping the citation: use Channon's filter, which fits Borbax
exactly.** Channon 2003, verbatim: "Because activity is intended as a measure of how much a
component both is used … and persists, I screen out (in each of the real and shadow
populations) **isolated occurrences: when a component occurs in the current snapshot but not
the previous one**." That is a time-persistence filter over consecutive snapshots. No
descent, no phylogeny, and it is the filter used by the very paper the activity statistics
come from. MODES itself notes this precedent: "prior open-ended evolution research has used
what is effectively a persistence filter with t = 1".

```rust
/// Isolated-occurrence screen (Channon 2003).
///
/// A species counts toward the statistics only if it was present in the
/// PREVIOUS snapshot as well as this one. Applied identically to the focal run
/// and the shadow, which is what makes their statistics comparable.
///
/// **Not** the MODES persistence filter (Dolson, Vostinar, Wiser & Ofria,
/// Artificial Life 25(1):50-73, 2019, DOI 10.1162/artl_a_00280). That one
/// requires a component's DESCENDANTS to persist for t generations and is
/// implemented over a phylogeny; Borbax species arise by reaction, not by
/// descent, so it has no analogue here. The MODES paper also treats the
/// persistence filter and shadow runs as ALTERNATIVE filters with different
/// blind spots — "shadow runs are also a viable filter option" — not as a
/// cheap substitute for one another.
pub struct IsolatedOccurrenceScreen {
    previous: BTreeSet<SpeciesId>,
}

impl IsolatedOccurrenceScreen {
    /// Returns the species present in both this snapshot and the last.
    pub fn admit(&mut self, current: &BTreeSet<SpeciesId>) -> BTreeSet<SpeciesId> {
        let admitted: BTreeSet<SpeciesId> =
            current.intersection(&self.previous).copied().collect();
        self.previous = current.clone();
        admitted
    }
}
```

**Test:**

```rust
/// Channon (2003): a component appearing in one snapshot and gone by the next
/// contributes nothing. This is what makes the statistics persistence-weighted
/// without a phylogeny.
#[test]
fn a_species_that_appears_for_one_snapshot_contributes_no_activity() {
    let mut s = ActivityStats::new(8, band());
    let mut screen = IsolatedOccurrenceScreen::default();
    for _ in 0..50 {
        s.snapshot(&screen.admit(&set(&[0, 1])), &c(&[0, 1]));
    }
    let before = s.cumulative;
    s.snapshot(&screen.admit(&set(&[0, 1, 7])), &c(&[0, 1, 7])); // 7 flickers
    s.snapshot(&screen.admit(&set(&[0, 1])), &c(&[0, 1]));
    assert_eq!(s.per_species_activity(SpeciesId(7)), 0.0,
               "a one-snapshot species accumulated activity");
    assert!(s.cumulative >= before);
}
```

---

### alr-113 · "the Lehman-Stanley archive measure" is a misattribution; ASAL's choice is verified, and it has a detail not to copy

**Where:** plan `3018–3019`.
**Claim:** the plan credits min-distance-to-archive jointly to ASAL and "the Lehman-Stanley
archive measure". Lehman & Stanley's novelty search uses the **mean distance to the k
nearest neighbours**, not the minimum — the mean is chosen precisely because a bare minimum
is too noisy. Min-distance-to-history is ASAL's choice alone.
**Confidence:** ASAL verified verbatim from source. Lehman & Stanley's `k` I did **not**
re-verify this pass; `k = 15` is from recollection and their published pseudocode, and
should be checked before being quoted in a comment.

**Citation for ASAL:** Kumar, Lu, Faldor, Nguyen, Ha & Tang, "Automating the Search for
Artificial Life with Foundation Models", Sakana AI, arXiv:2412.17799.
**Citation for novelty search:** Lehman & Stanley, "Abandoning Objectives: Evolution
Through the Search for Novelty Alone", *Evolutionary Computation* 19(2):189–223 (2011),
DOI 10.1162/EVCO_a_00025.

ASAL's actual metric, verbatim from `SakanaAI/asal` `asal_metrics.py`:

```python
def calc_open_endedness_score(z):
    """
    Calculates the open-endedness score from ASAL.
    The returned score should be minimized.
    """
    kernel = (z @ z.T) # T, T
    kernel = jnp.tril(kernel, k=-1)
    return kernel.max(axis=-1).mean()
```

So: **max similarity over strictly prior frames** (`tril(k=-1)` is the strict lower
triangle), then **mean over t**. Minimising max-similarity is maximising min-distance. The
plan's "min, not mean" is faithful at the inner reduction — but note the **outer** reduction
*is* a mean over time. The plan's `shape_novelty(current, history) -> f64` is a per-frame
quantity, which is fine, but the aggregate reported for a run must be the time-average of
per-frame minima, not something else.

**A detail from the source not to copy.** `jnp.tril(kernel, k=-1)` zero-fills the masked
region rather than filling with `-inf`. `foundation_models/clip.py` L27 confirms `z` is
L2-normalised, so `z @ z.T` is cosine similarity in `[-1, 1]` and negatives are possible: if
every prior frame has negative similarity, `max` returns `0` from the mask, not the true
maximum. ASAL's own `calc_illumination_score` does it correctly with `-jnp.inf`. In
Borbax's histogram version the analogous quantity is a distance (non-negative), so
`f64::INFINITY` is the right sentinel and the hazard does not arise — but say so, rather
than inheriting a mask value that only happens to be safe.

```rust
/// Novelty as *minimum* distance to any strictly prior state.
///
/// Following ASAL (Kumar, Lu, Faldor, Nguyen, Ha & Tang, arXiv:2412.17799),
/// whose `calc_open_endedness_score` takes the MAXIMUM similarity over
/// strictly prior frames and then the MEAN over time. Minimising max
/// similarity is maximising min distance; the outer reduction is a mean, so
/// the run-level figure is the time-average of these per-frame minima.
///
/// **Not** Lehman & Stanley's novelty-search measure (Evolutionary Computation
/// 19(2):189-223, 2011), which is the mean distance to the k nearest
/// neighbours in an archive — the mean is used there precisely because a bare
/// minimum is noisy.
///
/// The empty-history sentinel is `f64::INFINITY`, not zero. ASAL's reference
/// implementation zero-fills its mask, which is safe only because CLIP
/// similarities are almost always positive; a distance version zero-filled
/// would floor the score.
#[must_use]
pub fn shape_novelty(current: &Histogram, history: &[Histogram]) -> f64 {
    history
        .iter()
        .map(|h| l2_sq_unbiased(&current.bins, &h.bins))
        .fold(f64::INFINITY, f64::min)
}
```

**Test:**

```rust
#[test]
fn novelty_with_no_history_is_infinite_not_zero() {
    assert_eq!(shape_novelty(&h(&[1, 2, 3]), &[]), f64::INFINITY);
}

/// Min, not mean: adding a distant prior state must not raise the score.
#[test]
fn novelty_does_not_grow_with_run_length() {
    let cur = h(&[10, 0, 0]);
    let near = h(&[9, 1, 0]);
    let far = h(&[0, 0, 10]);
    let a = shape_novelty(&cur, &[near.clone()]);
    let b = shape_novelty(&cur, &[near, far]);
    assert!(b <= a, "novelty rose when history grew — this is a mean, not a min");
}
```

---

### alr-114 · "Complexity carrying capacity (Moreno & Ofria)" is misattributed, and it is a definition rather than a measurement

**Where:** spec §2.7 line 110; §20 line 1042; §22.1 line 1073.
**Claim:** the statement Borbax cites is from the OEE-II editorial, not from Moreno & Ofria,
and it is framed there as a *property OEE systems should have*, not as a measured scaling
law.
**Confidence:** verified from source.

Verbatim, from Packard, Bedau, Channon, Ikegami, Rasmussen, Stanley & Taylor, "An Overview
of Open-Ended Evolution: Editorial Introduction to the Open-Ended Evolution II Special
Issue", *Artificial Life* 25(2):93–103 (2019), DOI 10.1162/artl_a_00291, arXiv:1909.04430:

> "Moreno and Ofria consider the issue of indefinitely scalable complexity in their DISHTINY
> model. Channon's contribution to this special issue succinctly describes how indefinite
> scalability is linked to accumulation of adaptive success. Essentially, indefinitely
> scalable complexity implies that the 'complexity carrying capacity' of a chunk of space and
> time filled with interacting components **should, for systems that display OEE, be a
> thermodynamically extensive quantity**, which is to say it will scale with the space-time
> volume of the world. Indefinitely scalable complexity has become a hallmark of OEE systems,
> **essentially considered a necessary condition**."

Two things follow.

1. The sentence is the editors'. Moreno & Ofria's contribution is "Practical Steps Toward
   Indefinite Scalability: In Pursuit of Robust Computational Substrates for Open-Ended
   Evolution" — about substrate robustness. The evidence for indefinite scalability in the
   editorial is attributed to **Channon**'s Geb contribution. The spec's attribution should
   be corrected to Packard et al. 2019, or the claim restated.
2. **"should, for systems that display OEE"** is a conditional. It is a criterion for what
   counts as OEE, not a measured relationship. §20's mitigation — "If a run stalls, suspect
   region size *before* tuning parameters — this is a physical budget, not a knob" — reads
   the conditional backwards. It infers "the region is too small" from "the run stalled",
   which is valid only if Borbax already has OEE. Borbax having OEE is the thing being
   tested. Circular.

**Why it matters, mildly but concretely.** §20 tells a future debugger to suspect region
size before parameters, and §22.1 repeats it ("the more likely cause of saturation is total
space-time volume, not variety"). In V0 there is no region and no volume (`alr-109`), so the
advice is unactionable; in V1 it will send someone to enlarge the region when the honest
answer is "we do not know whether this system has OEE, so this diagnostic does not apply
yet."

**Fix — rewrite §2.7 line 110:**

> | **Complexity carrying capacity may be physical** — for systems that *do* display
> open-ended evolution, the complexity carrying capacity of a region of space-time is
> argued to be a thermodynamically extensive quantity, scaling with space-time volume
> (Packard, Bedau, Channon, Ikegami, Rasmussen, Stanley & Taylor, *Artificial Life*
> 25(2):93–103, 2019, DOI 10.1162/artl_a_00291, summarising Channon's and Moreno & Ofria's
> contributions to that issue). Note this is stated as a *necessary condition for* OEE, not
> as a measured law, so it cannot be used to explain a stall in a system whose OEE is
> unestablished. | Treat region size as a *candidate* explanation for a stall and test it by
> varying it, rather than assuming it (§20). |

**Test:** none — documentation. But the actionable half is a V1 experiment: make reaction
volume a swept parameter so "stalled because the budget is small" is falsifiable rather
than assumed. Recommend adding it to the V1 plan.

---

### alr-115 · Correction to my own previous finding: the 3b attribution is *partly* right

**Where:** spec §2.7 line 106.
**Previous finding (mine, now withdrawn):** "The 3b/3c subdivision is Channon 2003 (ALIFE
VIII), not Bedau/Snyder/Packard 1998. INTRODUCED BY THE FIX PASS."
**Corrected claim:** Channon made rows 3b and 3c explicit, but Bedau, Snyder & Packard 1998
*acknowledged* them. The spec's attribution is incomplete rather than wrong, and it uses a
numbering that a later paper has superseded.
**Confidence:** verified from source (Channon 2003 verbatim).

> "Table 1 in Bedau, Snyder & Packard (1998) only shows the first row (3a) for class 3, but
> **footnote 1 in that paper acknowledges the other rows (3b and 3c)**."

and the table caption: "based on table 1 from Bedau, Snyder & Packard (1998). **Rows 3b and
3c have been added to class 3** (see text)." Plus Channon's own footnote 2: "Note that Bedau
has since altered his class numbering scheme."

I did not open Bedau, Snyder & Packard 1998 itself. This is Channon quoting it directly,
which I consider adequate for a footnote-existence claim but not for the footnote's wording.

**The numbering has since moved, twice.** MODES, verbatim: "in order for a system to be
categorized among the most open-ended systems (**originally class 3, now class 4**) …
Channon suggested that class 3 open-ended dynamics should be broken up into three
subcategories … In parallel, **Skusa and Bedau refined the classification in a different
way**, inserting a new second class in which evolutionary activity was unbounded but no
novel components came into being. Such a situation would describe purely **ecological**
dynamics." So in MODES's merged numbering, Channon's 3b is **4b**, and there is a separate
class for "unbounded activity, no novelty".

**Why it matters, and it is more than pedantry.** The spec's §2.7 description — "unbounded
activity concentrated in a *bounded* diversity of components, so cumulative activity climbs
while **nothing new ever appears**" — merges two distinct classes. Channon's 3b has
`A_new` **positive** (see the table in `alr-100`); a system with `A_new = 0` is class 1 in
his scheme, verbatim: "Other possibilities exist with zero `A_new`, but these belong in class
1 (no evolutionary activity)." The "nothing new ever appears" case is Skusa & Bedau's
inserted class, which describes *ecology*, not stagnation.

Borbax should care about **both**, and they need different detectors: 3b/4b needs the median
activity statistic (`alr-102`), and the no-novelty class needs `A_new` read against the
shadow. The spec currently names one and describes the other.

**Fix — rewrite spec §2.7 line 106's parenthetical:**

> **Uncreative unbounded** — cumulative activity climbs while the system produces nothing
> genuinely new. The classification has been revised twice since Bedau, Snyder & Packard
> (1998), and two distinct classes matter here. **Channon's 3b** (Channon, "Improving and
> Still Passing the ALife Test: Component-Normalised Activity Statistics Classify Evolution
> in Geb as Unbounded", ALIFE VIII; rows 3b/3c made explicit from a footnote in Bedau,
> Snyder & Packard 1998) is unbounded *mean* activity with **bounded diversity** and
> positive `A_new`. **Skusa & Bedau's inserted class** is unbounded activity with **zero
> novelty**, which describes purely ecological dynamics. In the merged numbering used by
> Dolson et al. (2019) these are 4b and 3 respectively.
>
> Countermeasure: track novelty as a metric entirely separate from activity, never infer one
> from the other, use **median** rather than mean cumulative activity for the 3b/4b verdict
> (Channon's own correction — a single immortal component makes the mean unbounded for
> free), and read `A_new` against the neutral shadow.

---

### alr-116 · The neutral-network criterion needs the alphabet-dependent threshold, and Hamming-1 assumes fixed length

**Where:** plan `2948–2949`.
**Claim:** the plan asks for a largest-connected-component fraction with no reference point.
Random-graph theory on the generalised hypercube gives a percolation threshold that depends
on alphabet size, and Borbax generates its alphabet per universe. Separately, "Hamming-1
edges" presupposes fixed-length sequences, which Borbax polymers are not.
**Confidence:** threshold computed; the theory is second-hand (see below).

**Citation:** Reidys, Stadler & Schuster, "Generic properties of combinatory maps: neutral
networks of RNA secondary structures", *Bulletin of Mathematical Biology* 59(2):339–397
(1997), DOI 10.1007/BF02462007. I did **not** open this paper this pass. The threshold
formula and its derivation as a random-subgraph percolation result on `Q^n_κ` are from
recollection and from its wide secondary citation. The sanity check below is consistent with
it, but a reader who needs the exact statement of the theorem's hypotheses should open it.

Computed, `λ* = 1 − κ^(−1/(κ−1))`:

```
 kappa   lambda*    mean neutral neighbours at threshold
     2   0.50000     0.500 of  1
     4   0.37004     1.110 of  3     <- the textbook value, 1 - 4^(-1/3)
     8   0.25700     1.799 of  7
    12   0.20220     2.224 of 11
    32   0.10578     3.279 of 31
    64   0.06388     4.025 of 63
```

Two things fall out that are worth putting in the code. `λ*` is **decreasing** in `κ`, so a
larger generated alphabet percolates at lower neutrality. And the mean neutral degree at
threshold, `λ*·(κ−1)`, grows like `ln κ` — `ln 64 = 4.16` against the computed `4.025`. So
the criterion has a clean form: **a neutral network percolates once each sequence has on
average about `ln κ` neutral neighbours**, whatever the alphabet size.

**Why it matters.** A fixed neutrality bar rejects universes for having the wrong generated
alphabet size rather than the wrong folding map — and `universe_gen` chooses that alphabet.
The battery would reject workable universes systematically, and §7.2's stated purpose is to
reject *dead* ones.

**The second problem is larger and the plan does not mention it.** The theory is about a
fixed-length sequence space `Q^n_κ`. Borbax polymers grow and shrink by condensation and
cleavage, so the mutational neighbourhood includes length changes. Restricting the graph to
Hamming-1 edges undercounts the neighbourhood, which biases the measured neutrality
*downward* against a threshold derived for a graph that does not include those edges. The
comparison is not like-for-like. The honest options are (a) measure the LCC on the
fixed-length slice only and say so, or (b) include insertion/deletion edges and stop quoting
`λ*`, which no longer applies.

**Fix:**

```rust
/// Percolation reference for a neutral network on the generalised hypercube
/// Q^n_kappa: a random subgraph in which each vertex is neutral with
/// probability lambda percolates for lambda > lambda*, where
///
///     lambda* = 1 - kappa^(-1/(kappa-1))
///
/// (Reidys, Stadler & Schuster, Bull. Math. Biol. 59(2):339-397, 1997,
/// DOI 10.1007/BF02462007.)
///
/// This DEPENDS ON ALPHABET SIZE, and Borbax generates its monomer alphabet
/// per universe. Computed: lambda* = 0.370 at kappa=4 but 0.106 at kappa=32,
/// so a fixed neutrality bar rejects universes for their generated alphabet
/// size rather than for their folding map. Equivalently, and more memorably:
/// a neutral network percolates once each sequence has roughly ln(kappa)
/// neutral neighbours.
///
/// **Two caveats, both load-bearing.** The theorem is about a RANDOM subgraph;
/// a real folding map is not random, so this is a reference point and not a
/// pass/fail line — report the measured LCC fraction alongside it. And the
/// theorem is about FIXED-LENGTH sequences. Borbax polymers change length, so
/// the Hamming-1 graph omits insertion/deletion edges and understates
/// connectivity. Measure on the fixed-length slice and say so.
#[must_use]
pub fn percolation_threshold(kappa: usize) -> f64 {
    let k = kappa as f64;
    1.0 - det_math::powf(k, -1.0 / (k - 1.0))
}
```

**Test:**

```rust
/// kappa = 4 gives 1 - 4^(-1/3) = 0.3700394750525634, the value quoted
/// throughout this literature. If this moves, the formula was mistyped.
#[test]
fn percolation_threshold_matches_the_published_value_at_kappa_4() {
    assert!((percolation_threshold(4) - 0.370_039_475_052_563_4).abs() < 1e-15);
}

/// Decreasing in kappa. A larger alphabet needs LESS neutrality, so a fixed
/// bar over-rejects large-alphabet universes.
#[test]
fn percolation_threshold_falls_as_the_alphabet_grows() {
    for k in 3..64 {
        assert!(percolation_threshold(k) > percolation_threshold(k + 1),
                "threshold not monotone decreasing at kappa={k}");
    }
}

/// The battery bar must move with the universe's own alphabet.
#[test]
fn the_battery_bar_depends_on_the_generated_alphabet() {
    let small = Universe::with_alphabet_size(4);
    let large = Universe::with_alphabet_size(32);
    assert!(neutral_network_bar(&small) > neutral_network_bar(&large),
            "a fixed bar rejects universes for their alphabet, not their folding map");
}
```

---

### alr-117 · "RAF set count" is unreachable, and "the RAF stopped closing" is a blunt death criterion

**Where:** spec §15.2 line 839 ("RAF set count and size"); spec §9.6 line 511; plan `2452`
(`find_raf -> Option<RafSet>`).
**Claim:** `find_raf` returns the unique maxRAF, so "count" is 0 or 1. The literature's
answer is the count of **irreducible** RAFs, and finding *one* irrRAF is polynomial.
Separately, the maxRAF can shrink drastically while an irreducible core still closes, so
§9.6's death test fires late.
**Confidence:** verified by execution and from source.

From arXiv:2303.01809, verbatim: "note that this maxRAF contains two (nested) RAF subsets
(subRAFs) … The subRAF R′₂ is a so-called **irreducible RAF (iRAF)**, as it cannot be
reduced any further without losing the RAF property. Such subRAFs and iRAFs can be detected
with **repeated applications of the RAF algorithm to the maxRAF after removal of one or more
(random) reactions.**"

Executed, on the plan's own algorithm:

```
P5 two disjoint RAFs + one dead reaction    -> maxRAF = [0, 1]
P6 maxRAF = [0,1,2]; after losing the rim, still a RAF: [0]
   'the RAF stopped closing' is FALSE here even though 2/3 of it died.
```

**Why §9.6 matters more than §15.2 here.** §9.6 line 511: "death is emergent and precisely
defined: a compartment dies when its autocatalytic set no longer closes." P6 shows a
compartment can lose two-thirds of its RAF and still satisfy that criterion. The Chronicle
would report the compartment as alive right up to the moment the last irreducible core
fails — which is *defensible* (it is still, formally, self-maintaining) but it is not what
"damage load rises until the set stops closing" describes, and it makes the cause-of-death
census far less informative than §15.2 line 845 promises ("which reaction stopped being
catalysed, and what stopped catalysing it"). Tracking maxRAF *size* over time, not just its
existence, is what gives that report content.

**Fix:**

```rust
/// Find one irreducible RAF within the maxRAF.
///
/// "Such subRAFs and iRAFs can be detected with repeated applications of the
/// RAF algorithm to the maxRAF after removal of one or more (random)
/// reactions" — Hordijk & Steel, arXiv:2303.01809. Polynomial: |R| calls to
/// find_raf, each O(|X||R|^3) worst case.
///
/// **Removal order is ascending reaction index, never random** (§13.1). The
/// source says "random" because in that setting different irrRAFs are equally
/// interesting; here a reproducible one is required, and a fixed order gives
/// one. Note the irrRAF found is order-dependent — a different order finds a
/// different irrRAF — so this is "an irrRAF", not "the" one. Finding the
/// SMALLEST RAF is NP-hard and is not attempted.
#[must_use]
pub fn find_irr_raf(
    rxns: &[Reaction],
    catalysts: &[(SpeciesId, usize)],
    food: &BTreeSet<SpeciesId>,
) -> Option<RafSet> {
    let mut keep: Vec<usize> = find_raf(rxns, catalysts, food)?.reactions;
    let mut i = 0;
    while i < keep.len() {
        let candidate: Vec<usize> =
            keep.iter().copied().filter(|&r| r != keep[i]).collect();
        let sub: Vec<Reaction> = candidate.iter().map(|&r| rxns[r]).collect();
        let sub_cat: Vec<(SpeciesId, usize)> = catalysts
            .iter()
            .filter_map(|&(c, r)| candidate.iter().position(|&k| k == r).map(|j| (c, j)))
            .collect();
        if find_raf(&sub, &sub_cat, food).is_some_and(|s| s.reactions.len() == candidate.len()) {
            keep = candidate;      // reaction was removable; do not advance i
        } else {
            i += 1;
        }
    }
    let sub: Vec<Reaction> = keep.iter().map(|&r| rxns[r]).collect();
    let sub_cat: Vec<(SpeciesId, usize)> = catalysts
        .iter()
        .filter_map(|&(c, r)| keep.iter().position(|&k| k == r).map(|j| (c, j)))
        .collect();
    find_raf(&sub, &sub_cat, food).map(|s| RafSet {
        reactions: s.reactions.into_iter().map(|j| keep[j]).collect(),
        species: s.species,
    })
}
```

and for §15.2 line 839, replace "RAF set count and size" with:

> | maxRAF size (reactions, species) and irrRAF size | Organisation | Are self-sustaining systems appearing, and how much slack do they have? The count of RAF *sets* is always 0 or 1 — the maxRAF is unique. maxRAF size is the informative quantity, and the ratio maxRAF:irrRAF is the margin between the system and death (§9.6). |

**Tests:**

```rust
/// The maxRAF is unique and contains every other RAF (Hordijk & Steel), so
/// two disjoint RAFs come back as one set. "RAF set count" is 0 or 1.
#[test]
fn the_maximal_raf_is_the_union_of_disjoint_rafs() {
    let rxns = vec![uni(0, 1), uni(0, 2)];
    let cat = vec![(SpeciesId(1), 0), (SpeciesId(2), 1)];
    let raf = find_raf(&rxns, &cat, &food(&[0])).unwrap();
    assert_eq!(raf.reactions, vec![0, 1]);
}

/// An irrRAF is strictly smaller when the maxRAF has removable rim reactions.
#[test]
fn an_irreducible_raf_is_smaller_than_the_maximal_one() {
    let rxns = vec![uni(0, 1), uni(1, 2), uni(2, 3)];
    let cat = vec![(SpeciesId(1), 0), (SpeciesId(2), 1), (SpeciesId(3), 2)];
    let max = find_raf(&rxns, &cat, &food(&[0])).unwrap();
    let irr = find_irr_raf(&rxns, &cat, &food(&[0])).unwrap();
    assert_eq!(max.reactions.len(), 3);
    assert_eq!(irr.reactions.len(), 1, "irrRAF not reduced to its core");
}

/// §9.6's death criterion fires only when the LAST irreducible core fails.
/// This documents that, so nobody later reads a shrinking maxRAF as death.
#[test]
fn losing_most_of_the_max_raf_is_not_death() {
    let rxns = vec![uni(0, 1), uni(1, 2), uni(2, 3)];
    let cat = vec![(SpeciesId(1), 0), (SpeciesId(2), 1), (SpeciesId(3), 2)];
    assert_eq!(find_raf(&rxns, &cat, &food(&[0])).unwrap().reactions.len(), 3);
    // rim destroyed
    assert!(find_raf(&rxns[..1], &cat[..1], &food(&[0])).is_some(),
            "core stopped closing when only the rim was lost");
}
```

---

### alr-118 · The randomised-catalysis control needs a degree-preserving shuffle

**Where:** plan `3047–3054`; spec §15.3 line 860.
**Claim:** "Reassign which species catalyses which reaction uniformly at random, holding
catalysis density fixed" preserves the total number of catalytic assignments but destroys
the *degree distribution* on both sides. The result answers "shape chemistry versus a
uniform random network", which is not the question §15.3 poses.
**Confidence:** structural argument, not computed. Standard in network null-model
methodology (Maslov & Sneppen, "Specificity and stability in topology of protein networks",
*Science* 296:910–913, 2002, DOI 10.1126/science.1065103 — the edge-swap null model; the
result is a graph-theoretic one and is cited here for the algorithm only).

**Why it matters.** Shape-derived catalysis will not produce a uniform bipartite graph. A
species with a broad signature catalyses many reactions; a reaction with an unusual
transition geometry is catalysed by few. That heterogeneity is a *consequence* of the shape
chemistry, and a uniform reshuffle removes it along with everything else — so if a RAF
survives the real chemistry but not the uniform control, the conclusion available is only
"the degree distribution matters", which is weaker than and different from "the shape
chemistry matters". Hordijk & Steel's own random-catalysis analyses use a fixed mean
catalysis rate per molecule, which is a degree constraint.

**Fix — double-edge swap on the bipartite catalyst→reaction graph:**

```rust
/// Degree-preserving null model for the catalysis assignment.
///
/// Repeatedly pick two edges (c1, r1) and (c2, r2) and swap them to
/// (c1, r2) and (c2, r1). Every catalyst keeps its out-degree and every
/// reaction keeps its in-degree, so the only thing destroyed is WHICH
/// catalyst goes with WHICH reaction — which is exactly the shape-derived
/// content (spec §15.3, §8.5).
///
/// A uniform reassignment holding only total edge count fixed also destroys
/// the degree distribution, which shape chemistry produces and which is
/// therefore part of the result rather than part of the null. That control
/// answers "shape chemistry vs a uniform random network"; this one answers
/// "shape chemistry vs any network with the same degrees", which is the
/// question §15.3 asks.
///
/// Deterministic: edges are held in index order and both endpoints are drawn
/// from a `borbax_rng::Stream`. Rejected swaps (self-loops, duplicates) are
/// counted, not retried with a fresh draw, so the stream consumption is a
/// function of the swap count alone.
pub fn degree_preserving_shuffle(
    edges: &mut [(SpeciesId, usize)],
    rng: &mut Stream,
    swaps: usize,
) -> ShuffleReport {
    let mut seen: BTreeSet<(u32, usize)> = edges.iter().map(|&(c, r)| (c.0, r)).collect();
    let (mut done, mut rejected) = (0usize, 0usize);
    for _ in 0..swaps {
        let i = rng.next_range(edges.len() as u64) as usize;
        let j = rng.next_range(edges.len() as u64) as usize;
        if i == j {
            rejected += 1;
            continue;
        }
        let (c1, r1) = edges[i];
        let (c2, r2) = edges[j];
        if c1 == c2 || r1 == r2
            || seen.contains(&(c1.0, r2)) || seen.contains(&(c2.0, r1))
        {
            rejected += 1;
            continue;
        }
        seen.remove(&(c1.0, r1));
        seen.remove(&(c2.0, r2));
        seen.insert((c1.0, r2));
        seen.insert((c2.0, r1));
        edges[i] = (c1, r2);
        edges[j] = (c2, r1);
        done += 1;
    }
    ShuffleReport { done, rejected }
}
```

**Tests:**

```rust
/// The whole point: degrees on both sides are invariant.
#[test]
fn the_shuffle_preserves_both_degree_sequences() {
    let mut e = fixture_catalysis_edges();
    let before = (out_degrees(&e), in_degrees(&e));
    let mut r = Stream::new(1, Domain::Shadow, 0);
    degree_preserving_shuffle(&mut e, &mut r, 10_000);
    assert_eq!((out_degrees(&e), in_degrees(&e)), before);
}

/// And it must actually rewire — a null model that changes nothing is not one.
#[test]
fn the_shuffle_changes_the_assignment() {
    let orig = fixture_catalysis_edges();
    let mut e = orig.clone();
    let mut r = Stream::new(1, Domain::Shadow, 0);
    let report = degree_preserving_shuffle(&mut e, &mut r, 10_000);
    assert!(report.done > 1_000, "almost every swap was rejected: {report:?}");
    assert_ne!(e, orig);
}

/// A uniform reassignment does NOT preserve degrees. Asserted so the two
/// controls cannot be confused later.
#[test]
fn uniform_reassignment_destroys_the_degree_distribution() {
    let e = fixture_catalysis_edges();
    let mut r = Stream::new(1, Domain::Shadow, 0);
    let u = uniform_reassign(&e, &mut r);
    assert_ne!(out_degrees(&u), out_degrees(&e),
               "if these match, the fixture is too regular to test with");
}
```

---

### alr-119 · The shadow's equalised rate: constant or re-equalised must be a decision

**Where:** plan `3036–3041`; spec §15.3 line 858.
**Claim:** `k* = (Σ kᵢ nᵢ)/(Σ nᵢ)` matches the focal run's total removal flux **at the fork
instant only**. The plan does not say whether `k*` is then held constant or recomputed as
the shadow's composition drifts, and the two are materially different experiments.
**Confidence:** structural argument.

Held constant, the shadow is a clean neutral control: one fixed removal rate, matched to the
focal run at the fork, and any subsequent divergence in *flux* is itself a result. Re-
equalised, the shadow's rate depends on its own trajectory, which couples the control to the
thing being controlled and makes "run minus shadow" no longer a difference of one variable.

The spec's stated purpose (line 858) is "so the total removal *flux* matches the focal run at
the fork keyframe" — "at the fork keyframe" reads as constant-after-fork. Say so in the code
so it does not get "fixed" into re-equalisation by someone who notices the fluxes diverging.

```rust
/// The common decay rate applied to every species in the shadow.
///
/// k* = (Σ k_i n_i) / (Σ n_i), evaluated **once, at the fork keyframe**, and
/// then HELD CONSTANT for the shadow's lifetime.
///
/// Constant, not re-equalised: recomputing k* as the shadow's composition
/// drifts makes the control's parameter a function of the control's own
/// trajectory, so "run minus shadow" is no longer a difference in one
/// variable. The consequence — the shadow's removal flux will diverge from
/// the focal run's as the two compositions diverge — is a RESULT and is
/// reported, not corrected away (spec §15.3).
#[must_use]
pub fn equalised_rate(rates: &[f64], counts: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..rates.len() {          // index order; load-bearing (§13.1)
        num += rates[i] * counts[i];
        den += counts[i];
    }
    if den > 0.0 { num / den } else { 0.0 }
}
```

**Test:**

```rust
/// At the fork, total removal flux must match exactly.
#[test]
fn the_equalised_rate_matches_removal_flux_at_the_fork() {
    let rates = [0.1, 0.5, 2.0];
    let counts = [100.0, 10.0, 1.0];
    let k = equalised_rate(&rates, &counts);
    let focal: f64 = (0..3).map(|i| rates[i] * counts[i]).sum();
    let shadow: f64 = k * counts.iter().sum::<f64>();
    assert!((focal - shadow).abs() < 1e-12);
}

/// And it must not be recomputed afterwards.
#[test]
fn the_shadow_rate_is_frozen_after_the_fork() {
    let mut s = ShadowRun::fork_from(&keyframe(), &universe());
    let k0 = s.decay_rate();
    for _ in 0..10_000 { s.step(); }
    assert_eq!(s.decay_rate().to_bits(), k0.to_bits(),
               "shadow re-equalised: its parameter now depends on its own trajectory");
}
```

---

## Low

### alr-120 · `find_raf_present` and `rxn` are referenced and undefined

Plan `2440` — the doc comment says "Callers should therefore pass only species actually
present — see `find_raf_present`", and no such function exists in either plan file. Plan
`2288` and `2333` use a helper `rxn(a, b, c, d)` that is never defined; its four-argument
shape and the `[SpeciesId; 2]` reactant array interact with `alr-110`, so it is not merely
missing boilerplate. Either define both or reword the comment.

### alr-121 · Exit criterion 9 passes on a constant

Spec `1147` requires that "the metric families of §15.2 are computed — including the
persistence-weighted Bedau–Packard statistics — and a shadow can be forked from a beaker run
and compared against it." A hard-coded activity threshold satisfies every clause. Add the
provenance requirement:

> 9. The metric families of §15.2 are computed — including the persistence-weighted
>    Bedau–Packard statistics, using **median** cumulative activity for the class-3b verdict
>    (§2.7) — a shadow can be forked from a beaker run and compared against it, and the
>    activity band `[a₀, a₁]` is **demonstrably derived from that shadow's activity
>    distribution** rather than chosen. The shape-novelty metric reads within noise of zero
>    on a null trajectory **at three different occupied-bin counts**, which is what shows the
>    estimator is unbiased rather than merely tuned.

Same for plan Step 5b at `3155–3160`.

### alr-122 · Plastogenetic congruence measured from anneal endpoints is partly an annealer artefact

Plan `2955–2958` proposes folding under `k` diagnostic seeds and treating the distinct
results as the plastic repertoire. Those endpoints are a distribution over local minima of
one particular schedule, so measured congruence is confounded with the schedule. The
faithful variant costs no more: run a fixed-temperature Metropolis chain after burn-in and
sample the resulting Boltzmann ensemble at `T`, which makes the temperature an explicit
parameter of the measurement instead of a hidden one. Ancel & Fontana, "Plasticity,
evolvability, and modularity in RNA", *J. Exp. Zool.* 288(3):242–283 (2000),
DOI 10.1002/1097-010X(20001015)288:3<242::AID-JEZ5>3.0.CO;2-O. **I did not open this paper
this pass**; the ensemble-versus-endpoint point is a general statement about sampling and
does not depend on it, but the citation should be checked before it goes in a comment.

### alr-123 · The power-law increment is singular at `t = 0`

Plan `3063`. `a·b·t^(b−1)` diverges at `t = 0` for `b < 1`. Index windows from `t = 1`, and
assert it, so a later "off-by-one tidy-up" cannot reintroduce it:

```rust
#[test]
fn the_increment_index_starts_at_one() {
    let t = window_times(10);
    assert_eq!(t[0], 1.0, "power-law increment a·b·t^(b-1) is singular at t=0 for b<1");
}
```

### alr-124 · Information criteria must be computed on the training half only

Plan `3065` compares by AICc and BIC; plan `3077` holds out the tail. Both are right and
they must not be mixed: `n` and `rss` in the information criteria are the *training* `n` and
the *training* `rss`. Held-out error is reported alongside, never combined. Worth one line in
the `ModelVerdict` doc comment, because the obvious refactor is to compute the criteria on
everything.

### alr-125 · Guard `n − K − 1 ≤ 0`

Plan `3069`. At small `n` the AICc correction `2K(K+1)/(n−K−1)` goes negative and ranks the
most complex model best. Covered by the `Option` return in `alr-111`'s fix; noted separately
because it is a distinct failure and deserves its own assertion.

---

## What the plan gets right about this literature

More than the list above suggests, and several of these are choices that most projects in
this space get wrong.

- **RAF is RAF, not "RAF-like".** The implementation at plan 2452–2520 is Hordijk & Steel's
  Algorithm 1, verified line by line against arXiv:2303.01809, including the detail that the
  closure computation deliberately ignores catalysts (that is what separates RAF from CAF).
  The doc comment names the CAF distinction, the structural-not-dynamic limitation, and the
  NP-hardness of the inhibition variant. That is a more honest set of caveats than most
  published RAF implementations carry.
- **The bimolecular closure test** at plan 2314–2329 is exactly the test the previous four
  could not be: it is the first one using a two-reactant reaction, and it catches a bug that
  would have made `find_raf` return non-RAFs. I ran it; it works.
- **Increments, not cumulative totals** (§15.4 line 875). This is the single most common
  error in plateau analysis and the plan gets it right and explains why.
- **Held-out tail** (§15.4 line 876). This is the method that actually settled the LTEE
  question — truncating to early generations and testing which model predicted the rest —
  and the plan adopts it and states the reason ("the boundedness illusion is a *projection*
  failure").
- **Activity and novelty as separate families that may never proxy for one another**
  (§15.2 line 827). This is the correct reading of the failure mode, and the plan states it
  three times in three places because it is the kind of thing that erodes.
- **"Selection switched off means decay rates are *equalised*, not disabled"** (§15.3 line
  858), together with §23 criterion 7 explicitly naming decay-off as a *different*
  experiment. Conflating those two controls is a mistake that produces plausible output, and
  the spec has already spotted it and separated them.
- **The neutral-network criterion needs three numbers, not one** (plan 2941–2960). The plan
  correctly identifies that redundancy is necessary and nowhere near sufficient, and names
  connectivity, covering radius and the phenotype-frequency distribution. That is the actual
  content of the Fontana/Schuster work rather than a summary of it.
- **The plateau-fitting scale rule** (§15.4 line 874): "Fitting the power law on log–log
  axes and the saturating model on linear axes is not a comparison of models, it is a
  comparison of axes." Correct, and rarely stated.
- **`n(n−1)/2`, the `ln(0)` guard, and the segment-tree-over-Fenwick argument** are all
  correct and correctly reasoned in the comments.

---

## The single most important thing it still misses

**The activity statistics, as specified, cannot return a negative answer.**

`ActivityStats.cumulative` is a monotone accumulator (`alr-100`), so unbounded activity is
guaranteed by construction. `mean_cumulative` is `A/D`, which Channon showed goes unbounded
whenever a single component persists (`alr-102`) — and a beaker with a maintained food set
persists several by definition. The activity band that `A_new` depends on has no source in
the plan's step order and will be hard-coded (`alr-105`). The novelty metric's estimator has
a positive bias that grows with species count (`alr-106`), so novelty rises as diversity
rises. And the plateau fitter's error model puts 87% of its weight on the first 1% of the
trajectory (`alr-107`), so the tail — the only part that carries information about
boundedness — is effectively not fitted.

Every one of those defects points the same way: **toward declaring the run open-ended.** A
beaker containing nothing but inert feedstock would produce unbounded cumulative activity,
unbounded mean activity, rising shape-space novelty, and a power-law-beats-saturating
verdict. Task 21 Step 5b would tick criterion 9 and the run would be written up as a
success.

This is precisely Hintze's critique arriving from the inside. §2.7 line 109 quotes it — "a
trivial system can satisfy every published open-endedness criterion" — and the countermeasure
is "no single metric is trusted; every novelty measure is read against the neutral shadow."
That countermeasure does not help here, because the shadow would exhibit the same artefacts.
A monotone counter is monotone in the shadow too. `A/D` blows up in the shadow too. The
plug-in histogram bias is present in the shadow too. Subtracting one biased statistic from
another cancels the bias only if the two runs have the same diversity — and the entire point
of the comparison is that they do not.

**What to do about it, concretely.** Before Task 20b is written, add one test to
`activity.rs` that is not about any single statistic:

```rust
/// The whole metric suite, run on a beaker containing ONLY maintained
/// feedstock and no reactions at all. Nothing evolves; nothing can.
///
/// Every metric family must return its null answer. If any of them reports
/// open-endedness here, that metric is measuring its own implementation and
/// no amount of shadow subtraction will fix it — the shadow exhibits the same
/// artefact. This is Hintze's critique (spec §2.7) applied to our own code:
/// "a trivial system can satisfy every published open-endedness criterion."
#[test]
fn a_beaker_with_no_chemistry_is_not_open_ended() {
    let (u, g) = inert_fixture();            // food species only, zero reactions
    let run = run_beaker(&u, &g, 1_000_000);
    let shadow = ShadowRun::fork_from(&run.keyframe(0), &u).run(1_000_000);
    let band = ActivityBand::from_shadow(&shadow.per_species_activity(), 0.95, 0.999);
    let m = metrics(&run, band);

    assert!(m.activity.cumulative_is_bounded(),
            "A_cum unbounded on a beaker with no reactions");
    assert!(m.activity.median_cumulative_is_bounded(),
            "median activity unbounded on a beaker with no reactions");
    assert!(m.activity.new_activity < 1e-9,
            "A_new = {} on a beaker with no reactions", m.activity.new_activity);
    assert!(m.novelty.shape_novelty_mean.abs() < 1e-4,
            "shape novelty = {} on a beaker with no reactions",
            m.novelty.shape_novelty_mean);
    assert!(m.plateau.best_is_bounded(),
            "plateau fitter called {} on a beaker with no reactions",
            m.plateau.best_name());
}
```

That test is cheap, it runs in milliseconds, and it is the only thing in the plan that would
have caught all five defects at once. It should exist before any of the metrics do, and it
should be the first thing Task 21 Step 5b checks — ahead of the synthetic power-law-versus-
saturating discrimination, which tests the fitter but not the pipeline.

The plan's own §23 note explains why criterion 9 was added late: "The plan could not reach
its own bar." The observation this review adds is narrower and worse — as specified, the
metrics **would** reach the bar, on a universe with no chemistry in it.

---

## Where I could not verify, stated plainly

| Claim | What I actually have |
|---|---|
| Wiser, Ribeck & Lenski 2013 model forms and comparison statistic | **Paywalled, 403.** Abstract plus secondary descriptions. I know the two models were hyperbolic and power law, that both decelerate, that only the hyperbolic has an asymptote, and that the discriminating test was prediction from a truncated dataset. I have **not** read the equations and did not reproduce them. `alr-108`'s mathematics does not depend on the paper. |
| Bedau, Snyder & Packard 1998, Table 1 and footnote 1 | Not opened. Channon 2003 reproduces the table and quotes the footnote's existence directly; that is what `alr-115` rests on. The footnote's wording is not verified. |
| Mossel & Steel 2005 NP-hardness of RAF under inhibition | Not opened. Secondary descriptions plus recollection. The claim is uncontroversial and widely restated, but the doc comment at plan 2443–2445 asserts it and should carry the DOI. |
| Reidys, Stadler & Schuster 1997 percolation theorem | Not opened. The formula `λ* = 1 − κ^(−1/(κ−1))` is from recollection and wide secondary citation; the `κ = 4 → 0.370` check is consistent. The theorem's exact hypotheses (in particular how "randomly chosen neutral vertex set" is formalised) are not verified, and they matter for `alr-116`'s "this is a reference point, not a pass/fail line" caveat. |
| Lehman & Stanley `k = 15` | Recollection. The *shape* of the measure (k-NN mean, not minimum) I am confident in; the value of `k` should be checked before it appears in a comment. |
| Ancel & Fontana 2000 | Not opened. `alr-122`'s argument is about sampling and does not depend on it. |
| "Return to AlChemy" Figure 4(B) | **Verified** — but by reading a raster extracted from the PDF, not from a machine-readable table. I read it in both polarities and the values were identical. If anyone wants a second pair of eyes, the extracted image is at `…/scratchpad/fig4_right_inv.png`. |

## Referred to other reviewers

- **`rust-developer-expert`:** the `Reactants` enum in `alr-110` needs a real shape (my
  `as_slice` sketch cannot form a slice over the two-species arm); and plan line 1582's
  `Reaction` literal omits the `catalysis` and `catalyst` fields the struct at 1867
  declares, so that test does not compile.
- **`emergence-auditor`:** `alr-121`'s proposed exit-criterion wording adds a *derivation*
  requirement rather than a threshold, which I believe is on the right side of the "no
  hardcoded life" line — but it is their call. Also worth their eye: `alr-117` argues for
  tracking maxRAF *size* over time, which is a new metric, not a new mechanism, and should
  not trip the special-case rule.
- **`geometry-numerics-reviewer`:** `alr-116`'s second half — the Hamming-1 neighbourhood
  is the wrong graph for variable-length polymers — is a folding-map question as much as a
  literature one.
- **`determinism-auditor`:** `alr-117`'s `find_irr_raf` replaces the source's "random"
  removal order with ascending index; `alr-118`'s edge swap counts rejections rather than
  redrawing, so stream consumption is a function of the swap count alone. Both are
  deliberate determinism choices and should be checked.
