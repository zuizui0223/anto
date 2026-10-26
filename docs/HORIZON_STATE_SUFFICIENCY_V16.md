# ANTO v16 — horizon-specific sufficiency of evolutionary potential

2026-10-26. This is a **small exact synthetic counterexample**, not a measured return of ancestral flower colour, a newly discovered evolutionary law, a demographic field estimate, or evidence of a universal optimal waiting time.

## Research question

Can two plant-lineage histories with the same presently white flower, all three currently complete functional toolkit modules, the same immediate floral induction probability and the same descendant-reproduction rule be indistinguishable at one short endpoint yet differ in eventual trait accessibility?

**Yes, if the proposed state omits history-correlated future retention/transition information.** It does **not** follow that the past possesses an irreducible causal memory. When the complete present state includes the appropriate future transition law, adding a history label cannot improve the conditional prediction under the stipulated Markov closure.

This is a scoped ANTO theoretical complement to the EGWE criterion of endpoint-relevant state adequacy and to izu-core's separation of genetic accessibility from realised persistence. No outcome table or parameter is transferred from these repositories into the ANTO probability model.

## Two timescales and two endpoint contracts

Every modeled lineage begins with:
- The same white visible phenotype and all three required functional modules present.
- An induction opportunity in each generation with probability a = 1/2, independently between attempts.
- After an attempt that does not qualify, whole-toolkit competence survives the intergeneration interval with probability r_t. This retention process is a hidden transition law, not a measured gene deletion rate.
- No immigration, introgression, genetic mutation, developmental reconstitution, or retry cost.

The two separately evaluated endpoints are:

1. Visible: at least one expressed phenotype by model horizon T, regardless of later descendants.
2. Lineage: at least one expression attempt by T seeds a branch that *will eventually be nonextinct* under an independent Binomial(2, p) Galton–Watson process. This is retrospectively defined, does not imply fixation, and does not measure natural populations. Failed offspring branches do not consume later expression opportunities, so this is not the v13 one-shot preemption artefact.

The constructed daughter survival probability p = 2/3 gives per-expression ultimate branch nonextinction s = 3/4, whereas p <= 1/2 yields s = 0. Let b = a for visible events and b = a*s for descendant-lineage events. The exact qualified-event probability is

P(Y_T=1 | r_1,...,r_(T-1)) = sum over k=0,...,T-1 of [b * product over j=1,...,k of ((1-b)*r_j)].

For constant r the finite geometric expression is b*[1 - ((1-b)*r)^T]/[1 - (1-b)*r]. These are ordinary first-success/competence-retention identities, not new probability theorems.

## Result A — a one-generation sufficient-looking representation can fail at two generations

The source A and B worlds have r_t = 1/4 and r_t = 3/4, respectively, at all future intervals. All present-state variables listed above are matched, but the interval retention operator is not present in the *proposed coarse representation*.

| Endpoint | Horizon 1: A / B | Horizon 2: A / B |
|---|---|---|
| Any visible expression | 1/2 / 1/2 | 9/16 / 11/16 |
| First ultimately nonextinct branch | 3/8 / 3/8 | 111/256 / 141/256 |

Thus identical immediate expression and descendant-success forecasts are not evidence that the same current coarse state remains predictive for longer horizons.

For equally weighted historical cohorts, the exact Bayes-optimal Brier-score improvement from observing the history label (over a coarse model that cannot distinguish the two worlds) is (P_A - P_B)^2 / 4. At T=1 this is zero for both endpoints; at T=2 it is 1/256 for expression and 225/65536 for the lineage endpoint.

These are exact constructed expected errors, **not fitted cross-validation results** and not a sample-based claim about biological predictability.

## Result B — even two checked horizons need not qualify the third

The delayed-retention pair starts with an identical first interval r_1=1/2. The second interval is r_2=1/4 (A) or 3/4 (B), and all later intervals are identical. Immediate and two-generation outcomes match exactly; the three-generation contrasts are:

- Visible: P_B minus P_A = 1/32.
- Ultimately nonextinct lineage: P_B minus P_A = 75/2048.

In other words, apparent sufficiency at horizons one and two can fail at horizon three if a biologically relevant transition coordinate or age/developmental context is missing. An age-dependent schedule can be made Markov by including age/stage or the appropriate time-varying forcing in the full state; this is not a counterexample to a complete Markov model.

## What is and is not general

The formal transferable statement is conditional:

When an observable candidate state fails to capture a component of the future transition law that differs across histories, **predictive sufficiency can depend on the forecast horizon and specified endpoint**. For a transition-complete state with the same forcing law, the historical label has no additional predictive role in this generative example.

Important negative controls:
- Equal transition laws across histories yield zero residual historical prediction gain at all tested horizons.
- A subcritical daughter process has visible expression but no ultimate nonextinct descendant branch.
- More elapsed opportunities cannot decrease either first-qualification event in this no-cost, non-degrading-within-generation model.
- Exact fractions are checked against a separately implemented competent/success/lost state-tree oracle and the constant-hazard geometric solution.

The model does not estimate molecular regulatory recovery, historic dispersal, flowering frequency, ecological selection, selfing, introgression, lasting allele frequencies, or actual plant population survival. Here the three-module joint toolkit is stipulated to be complete in all histories; the v7 and v9 *joint module composition* counterexamples remain distinct. Costed multi-route re-expression is handled by v15 and does not inherit this model's hazard simplifications automatically.

## Practical tests before any biological claim

A genuinely decisive plant study should record, within independently replicated histories: baseline floral phenotype; same-individual, functional multi-organ ancestral-pathway assay; stage- and environment-dependent retention/activation; ancestry and gene-flow controls; seed survival and descendant reproduction over predeclared horizons. Evaluate held-out populations/clades at short and longer horizons, contrasting a visible-only model, a functional-state model, and a transition/context-complete model. Report residual history prediction only after stating assay limits; repeated photos and taxon-average expression cannot replace within-lineage follow-up.

ANTO owns this general calculation, EGWE owns its independent state-adequacy program, izu-core owns its own genetic and demographic experiments, and aza3 exclusively owns the empirical Cirsium mechanism/field analyses.

## Reproduce

- Main exact implementation: src/anto_horizon_state_sufficiency_v16.py
- Negative/closed-form/state-tree tests: tests/test_horizon_state_sufficiency_v16.py
- Deterministic fractions receipt: scripts/run_horizon_state_sufficiency_v16.py
- Dedicated CI: .github/workflows/anto-v16-horizon-sufficiency.yml

Commands:

    python -m unittest discover -s tests -p test_horizon_state_sufficiency_v16.py -v
    python scripts/run_horizon_state_sufficiency_v16.py --output /tmp/anto_v16.json

The inclusion gate requires exact synthetic replay and explicit nonclaim flags; no cross-repository numerical calibration is asserted.
