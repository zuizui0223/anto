# ANTO v13 — Re-expression probability is not demographic establishment

**As of 23 October 2026.** This is a mathematical, synthetic extension of [ANTO v12](EVOLUTIONARY_AFTERLIFE_V12_RECOVERY_DIRECTION_AND_CAPACITY.md) and its alternative first-recovery model. **No mean offspring rate, floral reversion, genetic fixation, seed germination, or population persistence was estimated from real plants.** All original Cirsium observation/biochemical/field analyses stay in [aza3](https://github.com/zuizui0223/aza3).

## Main question

If historical loss makes the alternative phenotype more likely to be the **first** recovered state, does that alternative actually contribute more to long-term demographic persistence?

The answer is **not necessarily**. It depends separately on the chance of first expression and survival of the descendant lineage. The difference is fundamental because “first visible pigmentation”, “retention for three model snapshots”, “positive long-term population growth”, and “adaptive genetic fixation” are *different estimands*.

## Model: two deliberately separated filters

1. **Previously validated v12**: Following a phenotypically silent interval of duration t, irreversible loss of common C, ancestral-route A, and alternative-route B resources changes the relative probabilities that the *first* detectable re-expression is A or B. These are the predeclared v12 continuous-time competing risks.
2. **New v13 filter**: Conditional on the first expression, one independent individual founds a hypothetical route-specific **Poisson Galton–Watson process** with mean descendant number R_A or R_B, constant over time. Its probability of ultimate nonextinction is the positive root

   q(R) = 1 − exp[−R q(R)], with q(R)=0 for R≤1.

   The *one-attempt* probability of a retained lineage from route j is p(first route j by H | t) × q(R_j), **not** the v12 expression probability alone.

3. **Strong assumption**: the first phenotypic event prevents any later route switch even if its descendant lineage fails. This is a restrictive experimental semantics, not a general biological fact. It is tested against a counterfactual where B remains possible after an earlier A event. Neither model includes densities, pollen limitation, immigration, seed banks, environmentally changing fecundity, inheritable trait switching, or plant allele-frequency fixation.

Branch survival is not equivalent to fitness measured by viable seeds or ecological establishment in natural populations.

## Exact synthetic reference contrast

Keep **the same** v12 expressed phenotype activation and all resource-loss hazards (a_A=0.4, a_B=0.2, loss C=0.02, loss A=0.30, loss B=0.01). Give A a much larger hypothetical Poisson mean (R_A=3.0), while B is just above its critical threshold (R_B=1.05).

Conditional survival is q_A=0.9404797907, q_B=0.0937018371; hence A is far superior *once expressed*. All time units are arbitrary and unrelated to known plant years or generations.

| Silent time t | First A (%) | First B (%) | A branch nonextinction, per initial founder (%) | B branch nonextinction (%) | Either nonextinct (%) |
|---|---:|---:|---:|---:|---:|
| 0 | 43.608 | 49.556 | **41.013** | 4.643 | **45.656** |
| 5 | 8.922 | 67.661 | **8.391** | 6.340 | **14.731** |
| 7 | 4.728 | 66.773 | 4.447 | **6.257** | **10.704** |
| 20 | 0.076 | 47.672 | 0.072 | 4.467 | **4.538** |

**New model-bounded result:** The A-vs-B contributions to ultimate lineage nonextinction cross at **t = 5.8896093398868** despite A's 10-fold-higher conditional branch survival. This is a *single specified hypothetical parameterization*, not a measured evolutionary-age threshold or a new general statistical theorem.

Notice the first B route dominates at t=5, but its demographic contribution remains **smaller** than A's in that row. First expression and surviving reproductive lineage cannot be silently substituted.

Even when B eventually dominates, total probability of a surviving lineage is decreasing: the apparent evolutionary direction shift is **not evidence of rising population-wide evolvability**.

## Sharper ecological criterion: how much B fecundity is actually necessary?

For first-expression probabilities p_A(t), p_B(t), and conditional ancestral q_A, B overtakes A iff

    q_B > p_A(t) q_A / p_B(t).

When the right side lies strictly between 0 and 1, invert the Poisson branching equation:

    R_B* = −log(1 − q_B*) / q_B*.

A **strictly larger** R_B is required for greater alternative contribution. At the illustrative R_A=3:

| Silent time | Critical B mean offspring R_B* |
|---|---:|
| 0 | 2.124151 |
| 5 | 1.067663 |
| 7 | 1.034856 |

The actual example R_B=1.05 is below the threshold at t=5 but above it at t=7. This is a **conditional frontier**, not a field-derived fitness estimate.

When B has R_B≤1, it may be the first re-expressed phenotype often while its long-run Poisson branching nonextinction stays exactly zero. Thus neither first return nor repeated finite observations alone establish demographic persistence.

## Essential negative control: removing first-event preemption

Consider **another hypothetical** environment with subcritical A (R_A=0.95) and supercritical B (R_B=1.5), while keeping the same v12 expression-resource parameters. Under the one-shot experiment, only the first route gets a chance to reproduce. The total nonextinction probability now apparently improves after silence:
- t=0: **0.288817**
- t=5: **0.394338**

It would be incorrect to claim that “resource loss increases evolutionary resilience.” The intervention removed many A-first occurrences that had no chance of surviving, and thereby exposed more B-first occurrences.

If B can arise even after an earlier A, assuming that early A does not change B's resources, B's **ever-expression × survival** chance instead **decreases**:
- t=0: **0.506793**
- t=5: **0.436201**

Therefore the apparent one-shot demographic gain is **endpoint- and preemption-dependent**. It is not a general ecological benefit of losing ancestral resources, and it must not be advertised as an adaptive trade-off proven in nature.

## Scope of mathematical verification

- Existing v12 continuous-time and independent simulation checks are preserved.
- Twelve new v13 tests independently guard Poisson extinction threshold, first-event vs demographic chance, crossover root, finite observation horizon, mass partition, subcritical B control, no-preemption adverse control, valid inputs and a deterministic **540-settings grid**.
- The grid is a factorial collection of arbitrary time/offspring means, **not an inferred probability distribution of plant reproductive strategies**. Do not quote proportions of grid settings as natural frequencies.
- Code: [v13 conditional model](../src/anto_afterlife_ecological_filter_v13.py); [negative tests](../tests/test_afterlife_ecological_filter_v13.py); [replay](../scripts/run_afterlife_ecological_filter_v13.py); [CI](../.github/workflows/anto-v13-ecological-filter.yml).

## Cross-repository empirical discriminant (not yet measured)

- **anto** owns the general first-event and branching filter only. It must NOT relabel events from two Petunia plants or one Iochrominae radiation as empirical calibration of route-specific R.
- **aza3** recently reanalysed **105 real Cirsium undulatum mothers** with two *contemporaneous reproductive routes* (apical/lateral seed successes) and found pooled versus eight-stratum co-success signs differ, with exact two-sided p=0.09078. This is a real joint-route denominator/result, but it does **not** measure ancestral vs novel floral-colour regain, cross-generation descendants, or demographic Galton–Watson means. [Original mother-specific audit](https://github.com/zuizui0223/aza3/blob/main/docs/CIRSIUM_UNDULATUM_WITHIN_MOTHER_ROUTE_COVARIANCE_20261023.md).
- **egwe** distinguishes simulated first reappearance from 12-generation trait persistence under its own reproductive life cycle. It is a *different simulation*, not plant-data confirmation of v13.

A true plant experiment will need **verified genetic and floral biochemical route identity**, relevant pollen/seed/seedling/descendant fitness measures by lineage, whole-lineage genotyping and an explicit repeated-reexpression protocol. A new coloured flower cannot be assumed to found a self-sustaining lineage; separate ongoing re-expression opportunities can change the first-event model conclusion.

## Publication assessment

No Ecology Letters/NEE-ready empirical general law arises solely from adding Poisson branching to a known competing-risk result. It is a useful **quantitative experimental discriminant and failure-mode** for later cross-clade tests. Stronger originality would require genuinely independent source-qualified plants showing a specific resource-inventory × route × descendant-fitness interaction, preferably prospectively predicted by this model. Until then, preserve v13 as a general-theory extension and do not invent biological establishment data.
