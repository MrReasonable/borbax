# ALife literature re-review — findings
(Transcribed by the parent: this agent has no write tools. Fix its `tools:` line.)

## VERIFIED CORRECT (fixes that landed)
- alr-000 RAF closure fix is correct. Walked A+B->C with A in F, B not: missing=2,
  popping A decrements to 1, 7 never queued, C never produced, prune drops r0,
  returns None. Duplicate-reactant case also correct. Matches Hordijk & Steel 2004.
- alr-001 All four doc-comment claims check out: maximality via closure-under-union,
  RAF-not-CAF, structural-not-dynamic, inhibition NP-hard (Mossel & Steel 2005).
- alr-002 `1.0 - u1` is correct. Trades tau=+inf at p~2^-53 for tau=0, which is
  harmless (clock stalls one event). That reason should be the test's assertion.
- alr-003 `n(n-1)/2` correct (Gillespie 1977).

## HIGH
- alr-010 [misattribution] §2.7. The 3b/3c subdivision is **Channon 2003** (ALIFE VIII),
  not Bedau/Snyder/Packard 1998. First three clauses verified correct; the fourth is
  Channon's extension. INTRODUCED BY THE FIX PASS.
- alr-011 [scope-gap] Task 18. **Volume V is not a parameter anywhere.** Gillespie's
  c = k/V (or 2k/V same-species), so beaker volume is silently baked into
  rate_prefactor. Consequence: §2.7's Moreno-Ofria "region size is a complexity
  budget" countermeasure is NOT implementable in V0 — a rejected universe cannot be
  distinguished from a beaker too small to nucleate. Fix before Task 18 lands.
- alr-012 [scope-gap] Task 20b Steps 1 and 3 are presented as independent and are not.
  The activity threshold must come from the shadow cross-over (Rechtsteiner & Bedau
  1999 — gave 50 vs the arbitrary 10 used previously). Step 1 cannot compute
  `diversity` or `new_activity` without Step 3's output. As written an implementer
  hard-codes a constant, which is the practice the paper was written to replace.
- alr-013 [numerical-error] Task 20b. **The null-trajectory test will fail as written.**
  Plug-in histogram-distance estimators are positively biased; subsampling to constant
  N removes the N-dependence but not the K-dependence (K = occupied bins, grows with
  species count). Change the criterion from "reads ~0" to "constant in time, matching
  the analytic null" — or bias-correct. Use the null as a subtracted baseline.
- alr-014 [numerical-error] Task 20b. Novelty increments are **counts**: variance scales
  with mean, so homoscedastic-Gaussian RSS is violated exactly where models differ most
  (saturating predicts small late increments, RSS weights them equally). Use
  Poisson/quasi-Poisson likelihood, or variance-stabilise (Anscombe/sqrt) first.
  Otherwise "fit increments" swaps one wrong error assumption for another.

## MEDIUM
- alr-020 [numerical-error] AICc `K` must include sigma^2 => K=3 for 2-param models,
  K=2 for linear. Cancels in AIC/BIC; does NOT cancel in AICc. At n=10 the ΔAICc error
  is ~1.07 and biases *toward* the 2-param models. (Burnham & Anderson 2002 §2.2.)
- alr-021 [scope-gap] The bounded model in the source is **hyperbolic**, not
  exponential-saturating (Wiser, Ribeck & Lenski 2013). Hyperbolic decays as t^-2 and is
  much harder to separate from a power law — it is the discrimination that actually
  tests the boundedness illusion. Costs one more model.
- alr-022 [misattribution] MODES persistence filter is **lineage** persistence on a
  **phylogeny**. Borbax has no phylogeny (species arise by reaction, not descent) so it
  is not implementable as published. Also the plan calls it "the substitute for a full
  shadow" while building the shadow — inverts the relationship. Keep a time-persistence
  filter, drop or qualify the citation.
- alr-023 [misattribution] "Lehman-Stanley archive measure" is wrong: theirs is k-NN
  mean with k~15, chosen because a pure minimum is too noisy. Min-distance is ASAL's
  choice. Verified against SakanaAI/asal `asal_metrics.py` — max similarity over
  strictly prior frames, then mean over t. Plan's "min not mean" is faithful to ASAL.
- alr-024 [scope-gap] Task 20b. Equalised rate k* = (sum k_i n_i)/(sum n_i) is correct
  at the fork ONLY. State whether it is held constant (cleanest) or re-equalised
  (couples the shadow to its own trajectory). Must be a decision, not an accident.
- alr-025 [scope-gap] Randomised-catalysis control: "holding density fixed" under-
  specified. Uniform reassignment also destroys the degree distribution, so it measures
  "shape chemistry vs uniform" not "vs same-density random". Add a **degree-preserving
  edge shuffle** on the bipartite catalyst->reaction graph. Only that answers "because
  of the shape chemistry".
- alr-026 [scope-gap] §15.2 asks for "RAF set count and size" but `find_raf` returns
  Option<RafSet> — count is 0 or 1. The literature's answer is the count of
  **irreducible** RAFs. One irrRAF is polynomial (~15 lines on find_raf); smallest is
  NP-hard. Also bears on §9.6: maxRAF can shrink drastically while a core irrRAF still
  closes, so "the RAF stopped closing" is a blunt death criterion.
- alr-027 [numerical-error] Task 20 battery. Neutral-network percolation threshold is
  lambda* = 1 - kappa^(-1/(kappa-1)), depending on **alphabet size kappa**, which Borbax
  generates per universe. A fixed LCC-fraction bar would reject universes for having a
  small alphabet. Pure random-graph theory (Reidys, Stadler & Schuster 1997).

## LOW
- alr-030 [type-inconsistency] `find_raf_present` referenced in the doc comment,
  defined nowhere.
- alr-031 [scope-gap] Exit criterion 9 should require the threshold was *derived from
  the shadow*, not merely that the statistics compute — else it passes on a constant.
- alr-032 Nothing computes the shadow's own activity distribution, which is the object
  Rechtsteiner & Bedau's method consumes. Missing ~5 lines connecting Steps 1 and 3.
- alr-033 Plastogenetic congruence: k-seed anneal endpoints are a distribution over
  local minima of a particular schedule, so measured congruence is partly annealer
  artefact. Faithful variant costs no more: fixed-temperature Metropolis chain after
  burn-in gives a Boltzmann ensemble at T (Ancel & Fontana 2000), making T explicit.
- alr-034 [could-not-verify] AlChemy 16%/60% figures could not be confirmed from the
  paper; one render gave ~13% coexistence / 17% mutual destruction / 70% dominance,
  and the paper's prose emphasises *dominance*. Fix against the published table or
  drop the numbers. Qualitative conclusion (never assume composition) is unaffected.
  Citation: Chaos 34(9):093142, 2024 — "Return to AlChemy" is the subtitle.
- alr-035 Power-law increment a*b*t^(b-1) is singular at t=0 for b<1. Index from t=1.
- alr-036 AICc/BIC must be computed on the training half only, reported alongside
  (not mixed with) held-out error.
