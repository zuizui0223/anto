# ANTO v17 — Exactly when does ecological memory become predictively relevant?

26 October 2026. **General synthetic theory only.** An exact second-generation result connects the [horizon-specific sufficiency counterexample in ANTO v16](HORIZON_STATE_SUFFICIENCY_V16.md) to the [shared-Markov-environment multitype reproduction model](SLK_ANTO_MARKOV_ECOLOGY_V6.md) integrated from #19. This is not a measured floral evolutionary reversal, novel branching-process theorem, source-derived environmental persistence parameter or direct SLK fitness estimate.

## Central claim, with its scope

Two prospective histories can have **the same initially nonexpressing but completely competent phenotype state, the same current one-generation offspring law, the same mean environment composition and the same immediate expression probability**, yet different multi-generation retention prospects because the proposed observed state omits the future environmental transition operator.

This is **not** intrinsic, unmeasurable evolutionary memory: once the full current developmental/genetic state and the shared environmental transition law (or known time-varying forcing) are supplied, this Markov model needs no separate historical label. The hypothetical label is a proxy for a missing p parameter.

## Biological endpoints must be separate

The model starts with one nonexpressing yet reproductively active **C** parent that retains all necessary joint functional modules. Competent offspring may express (**D**) or remain silent (**C**). Children losing competence exit the counted genealogy, not necessarily the living plant population.

- First-expression event: at least one D child appears from the initial C parent.
- Finite competent-lineage retention: at least one D or C descendant exists at generation H.
- Eventual nonextinction or genetic fixation: **not calculated here**. A positive finite-H probability is insufficient for either conclusion.

The model reuses v4/v6 without mutation, immigration, unknown reacquisition, externally calibrated selection or independent environmental histories assigned to different descendants.

## Exact two-generation theorem for the declared synthetic rates

In environment A, the chance of producing a viable child per reproductive slot is q_D=9/10 and q_C=2/5. Environment B reverses these values. Each child retains the full relevant toolkit with probability k=4/5. The environment is a symmetric stationary Markov chain shared by all descendants within the generation, with probability p of remaining in the same environment. Its marginal occupancy is always P(A)=P(B)=1/2 regardless of p.

The parameter w is the probability that a competent child switches expressing/silent state compared with the parent.

The one-generation probability of at least one competent offspring from a C parent is

\[
S_1=1-\frac{(1-kq_{C,A})^2+(1-kq_{C,B})^2}{2}
=\frac{456}{625}=72.96\%.
\]

This does not depend on p or w. First D-child expression depends on w but not p.

For environments E0=i and E1=j, the original C parent's probability of **zero competent grandchildren** is

\[
f_{ij}(w)=\left[1-kq_{C,i}+kq_{C,i}
 \left(w(1-kq_{D,j})^2+(1-w)(1-kq_{C,j})^2\right)\right]^2.
\]

Integrating the exact shared-environment path probabilities gives

\[
S_2(p,w)=1-\tfrac12\left[p(f_{AA}+f_{BB})+(1-p)(f_{AB}+f_{BA})\right].
\]

For this explicitly fixed scenario, the linear derivative is exactly

\[
\frac{\partial S_2}{\partial p}
=-\frac{361824}{9765625}(2w-1).
\]

Therefore environmental persistence increases two-generation capacity retention when w<0.5, has no effect when w=0.5, and decreases it when w>0.5. The reversal arises at generation two, not from the longer numerical enumeration or stochastic simulation.

| Expression switching w | S2 at p=0.1 | S2 at p=0.9 | Difference (p=0.9 minus 0.1) |
|---|---:|---:|---:|
| 0.1 | 57.437356% | 59.808606% | +2.371250 percentage points |
| 0.5 | 59.355306% | 59.355306% | 0 |
| 0.9 | 59.808606% | 57.437356% | −2.371250 percentage points |

The exact signed difference in the off-center cases is respectively ±5789184/244140625 = ±0.023712497664.

The larger previously integrated v6 14-generation contrast (16.16993% versus 34.50465%, with ranking reversed at w=0.9) is a **separate finite-horizon genealogy-presence result**. Neither endpoint measures natural genetic reversal, stable restoration of a flower trait, long-run extinction or fixation.

## Horizon-specific value of a proposed historical proxy

If two synthetic history labels have equal prior weight, identical currently recorded states but different unrecorded environmental persistence p, and their horizon-specific event probabilities are A and B, then the exact expected Bayes-optimal Brier-score gain of observing the label is

\[
\Delta\mathrm{Brier}=\frac{(A-B)^2}{4}.
\]

At horizon one it is **zero**. At horizon two it is **strictly positive for w unequal to 1/2** in the declared witness. At w=1/2 it vanishes at all tested horizons: offspring expression assignment does not preserve the parent's expression state, so environmental autocorrelation stops affecting this model's competent-lineage endpoint.

This specifically connects ANTO v16's horizon-dependent state-sufficiency criterion with v6's multitype reproduction under a shared random environment. No statistical or demographic result is transferred from EGWE or izu-core. Under a transition-complete Markov state, p is an observed or known model component and the label is redundant by construction.

## Falsifiers and prospective organismal tests

The transferable ecological question is not whether environmental persistence is always beneficial, but whether **regime duration interacts with the retention time of alternative functional expression states**.

A biological test requires independent tracked reproductive lineages, properly measured gene and regulatory *functional* competence, recurrent expression transitions, absolute whole-life reproduction of expressing and silent competent states, relevant partner-environment histories, and descendants assessed at predeclared multiple horizons. Immigration, introgression, genetic repair, stage, ploidy and induced plasticity must be separate.

- If expression-state inheritance is absent (the w=0.5 analogue) under the specified symmetry, the model predicts no environmental-autocorrelation effect on competent-lineage presence.
- If environment-specific reproductive returns are equal, this particular environmental-timescale mechanism disappears.
- If a proposed complete-state model still leaves a residual history effect on held-out lineages, test omitted molecular competence, age/stage, environmental feedback and migration, rather than treating the residual as evidence for a metaphysical memory.
- A flower-colour difference alone does not establish SLK multifunctionality or a recovered ancestral genetic route.

ANTO owns only the general theoretical test. Real Cirsium molecular/field/head-fitness evidence belongs exclusively to aza3; separate SLK, EGWE and izu-core studies are methodological comparators, not independent replications or calibrated parameter sources.

## Reproducibility contract

- Exact rational source: src/anto_environment_horizon_bridge_v17.py
- 12 independent and adverse tests: tests/test_environment_horizon_bridge_v17.py
- Frozen JSON receipt: scripts/run_environment_horizon_bridge_v17.py
- Dedicated CI: .github/workflows/anto-v17-environment-horizon.yml

The closed-form environment-pair formula is independently cross-checked against full rational path enumeration using the original v6 offspring PGF. Finite-horizon floats from the earlier v6 method provide another consistency check. CI fails if first-horizon equivalence, horizon-two reversal, w=0.5 null control, Brier-score identity or source/nonclaim contract is violated.
