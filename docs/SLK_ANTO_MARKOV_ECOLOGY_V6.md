# ANTO × SLK v6 — ecological memory selects which functional capacities remain possible

**20 October 2026. Theoretical comparison under a shared stationary Markov environment. Not empirical plant trait re-evolution, ecological fitness calibration, or a novel theorem about branching processes.**

## Biological question

Suppose a plant retains the developmental/genetic capacity to express a previously absent structure, but expressing or suppressing that structure changes reproductive performance depending on the ecological partner environment. **Can environmental persistence—not environmental mean, not the first re-expression event—change whether the capacity remains in reproducing descendants?**

ANTO v4 distinguished visibly expressing D from nonexpressing but reproducing *functionally competent* C. ANTO v5 showed a temporal effect in deliberately deterministic ecological cycles. V6 uses an **ergodic Markov environment shared among all descendants in each generation**. It asks about genuinely unpredictable environmental runs instead of only ABABAB and AAAABBBB cycles.

## Ecological comparison: identical averages, different memory

Environment A: per-reproductive-slot probabilities `q_D=0.9,q_C=0.4`. Environment B reverses these: `q_D=0.4,q_C=0.9`. Each biological child retains all required functional resources with probability `1-d=0.8`. Conditional on retention, D/C offspring change their expression state with either symmetric probability `w=0.1`, `w=0.5` or `w=0.9`. These are arbitrary **conditional model parameters**, not natural pollen, climate or gene-loss rates.

Let `p=P(E_(t+1)=E_t)`, with symmetric transitions between A and B. The chain begins at its invariant distribution `P(A)=P(B)=1/2`. Thus the marginal proportion of BOTH environments is 50% at every generation for all p. The lag-one environmental autocorrelation is `2p-1`; expected length of a run is `1/(1-p)`:

- `p=0.1`: frequent alternating ecology; runs average 1.11 generations
- `p=0.9`: persistent ecological conditions; runs average 10 generations

The stochastic environments are **shared by all offspring in a given generation**, not assigned independently to individuals. They are exogenous: no dispersal, density dependence, partner evolution or feedback of phenotype abundance on weather or pollinator communities.

## Exact finite-genealogy outcome: a clear ranking reversal

For a single initially **nonexpressing, competent C parent**, calculate the chance of at least one descendant with the complete functional toolkit at generation 14. For every environmental sequence of length 14, take the exact v4 two-type offspring-generating-function composition in its biological chronological order; average the resulting extinction probability under the exact Markov path probabilities. There are `2^14=16,384` distinct environmental histories; the path weights sum to 1, and the within-history PGF calculation uses floating-point arithmetic. An independent **rational** path-enumeration oracle checks the first four generations.

| D↔C inherited expression switching w | Frequent environmental alternation p=.1 | Long environmental runs p=.9 |
|---|---:|---:|
| .1 — strong expression-state retention | **16.16993%** | **34.50465%** |
| .5 — offspring expression has no parental memory | 23.90731% | 23.90731% |
| .9 — offspring expression tends to flip | **34.50465%** | **16.16993%** |

**This result does NOT say the trait has re-evolved or established.** It is a finite presence event for *competent descendants*, whether they currently express the phenotype or not. Extinction of competent descendants is not extinction of organisms; plants with an irrecoverable loss can reproduce outside the focal genealogy.

Within each w comparison the next-generation probability of at least one visibly expressing D child from a C parent is also exactly the same for p=.1 and p=.9, because the initial environment is stationary and identical. At w=.1 it is `10.0896%`; at w=.9 it is `68.4576%`. The differences in competent genealogy survival cannot be attributed to a different initial phenotype or one-off induction rate.

## Long-run growth: sampled typical history is not an average environment

A second, deliberately distinct endpoint is the **top pathwise growth** of the product of the v4 mean offspring matrices along long Markov environmental histories. We normalize the two-type mean vector after every generation, record log growth and exponentiate. For three deterministic independent random seeds of 100,000 generations each (2,500 burn-in), sample means are:

| w | p=.1 typical-growth estimate | p=.9 typical-growth estimate |
|---|---:|---:|
| .1 | 0.97698 | 1.18901 |
| .5 | 1.04000 | 1.04000 |
| .9 | 1.18943 | 0.97708 |

For this strictly positive two-type process, Lyapunov growth is a standard theoretical guide to long-term branching persistence under appropriate stationary-ergodic assumptions. These numbers are **Monte Carlo estimates, not rigorous error bounds or a certified survival classification**. The separately computed 14-generation probabilities above are actual finite-time event probabilities under the declared model (subject to numerical rounding).

The arithmetic mean offspring matrix has Perron growth **1.04** in all six comparisons. It cannot stand in for chronological multiplication in the stochastic environments. The distinction between *quenched* typical growth and *annealed* population expectations is established random-environment branching theory, not a new mathematical result.

Because the same toolkit loss `d` multiplies every generation's mean matrix by `1-d`, the estimated critical loss threshold is
`d*≈1-(1-d)/g_typical`.
For low w this is approximately **.1811 with alternating ecology versus .3272 with persistent ecology**. The displayed d=.20 lies between these estimates, illustrating different likely long-run retention regimes under the same average ecological occupancy. The threshold values inherit the Monte Carlo uncertainty; they are neither plant mutation rates nor formal numerical certificates.

## Biological predictions and counterexamples

**1. Temporal continuity can preserve dormant-looking evolutionary options.** When expression states are strongly retained, a long favorable episode allows lineages to remain in a productive state; frequently switching partner ecology can erode the competent genealogy even though both partner states are equally common.

**2. Temporal continuity does not always help.** When expression-state inheritance tends to flip, alternating ecology can match the inherited state changes, and long runs become detrimental. In this symmetric witness the ranking reverses. It is not a universal evolutionary optimum for w=.9, and no mutation-selection process of w itself is simulated.

**3. A valid negative control removes the effect.** At w=.5 the offspring's visible state is independent of the parent's type (conditional on retention). The run-duration effect on finite survival disappears under the stated symmetrical environment/offspring configuration; sampled pathwise growth is exactly 1.04. Changing the baseline symmetry can remove this equality.

**4. A photograph cannot identify this mechanism.** Same initial visible phenotype, first induced-expression chance, environmental mean and current functional-gene availability do not determine future capacity retention if expression-state inheritance and environmental autocorrelation differ. Same-plant functional assays across generations are necessary.

### Empirical biological design, not yet executed

The prospective unit is a **tracked lineage**, not a photo or a set of unrelated published species. For an SLK-compatible flower system with documented opposing functions of the **same architecture** and alternative integrated/differentiated forms, independently measure:

- Partner-environment sequences and their lagged autocorrelation/run lengths rather than just mean visitation or climate.
- Consecutive-generation expressed/nonexpressed architecture states, checking inheritance versus developmental plasticity and explicitly distinguishing them from seed dormancy.
- Validated co-retention of all necessary coding/regulatory **functional modules in the same lineages**, with leaf/flower assay scope, ancestry and introgression controls.
- Absolute reproductive success and persistence of the functional capacity, not merely pollen deposition, induced colour or higher relative fitness.
- Paired lineages under matched marginal ecology but differing duration histories, and out-of-sample predictive evaluation.

For the **Cirsium colour-regain system**, this cannot be called an SLK joint test without an independently established shared-function/division-of-labour component. All Cirsium specimen-level work belongs to [aza3](https://github.com/zuizui0223/aza3). A heterantherous floral comparison with retained mechanical support could be a more direct SLK candidate, if its genetic capacity is independently measurable.

## Source attribution and scientific ceiling

Environmental timescale × phenotype switching effects have prior art, notably Acar, Mettetal & van Oudenaarden (2008, *Nature Genetics*, DOI [10.1038/ng.110](https://doi.org/10.1038/ng.110)). The mathematical connection between multitype Markovian branching environments and Lyapunov exponents was studied by Hautphenne & Latouche (2016, DOI [10.1007/s10955-016-1474-3](https://doi.org/10.1007/s10955-016-1474-3)); generic branching/ergodic results are not a novelty claim. The contribution here is a **specific ANTO functional-toolkit retention/loss contrast**, using an unchanged source-registered v4 process, exact exhaustive finite environmental enumeration, control of marginal environment and initial expression, and an ecologically falsifiable capacity-retention interpretation.

No natural re-expression, molecular recovery, evolved optimal switching, organismal extinction, frequency-dependent SLK payoff calibration, demographic fixation, or population rescue is inferred. The Markov chain is stipulated, not fitted; environmental randomness is common within a generation. This v6 extension is intentionally a child draft of the v1–v5 ANTO branch, leaving stable ANTO main and the SLK Am Nat flagship untouched.

## Reproduce

- Model: `src/anto_slk_markov_ecology_v6.py`
- 13 focused tests, including independent exact rational environmental path/PGF oracle: `tests/test_slk_anto_markov_ecology_v6.py`
- Deterministic repeated-seed result generator: `scripts/run_slk_anto_markov_ecology_v6.py`
- CI: `.github/workflows/slk-anto-markov-ecology-v6.yml`
- Compare matched pairs `w=.1,.9` and `p=.1,.9` and the explicit w=.5 negative control.

The CI performs 13 focused tests, a full regression suite, then reruns deterministic summaries. Numerical finite-survival values are double-precision approximations of exhaustive weighted PGF sums; the Lyapunov-growth values are sampled quantities and must not be described as exact.
