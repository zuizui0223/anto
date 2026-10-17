# anto v2 synthetic research receipt — 2026-10-17

**Status: executed locally, source and replay shipped in this branch; NO real biological data admitted.** All rates and assay accuracies below are illustrative, not empirical thistle parameters. An independent GitHub Actions replay is specified by the accompanying workflow; do not claim it passed until the workflow status is checked.

## Reproduction
- Existing one-lineage continuous-time baseline: src/anto_model.py
- New finite-region/nonoverlapping-generation stochastic model: src/anto_colonization.py
- Synthetic mechanism identifiability: src/anto_identifiability.py
- Monte Carlo runner: scripts/run_synthetic_baselines.py
- Tests: tests/test_model.py and tests/test_colonization.py
- Commands:

    python -m pip install numpy scipy
    python -m unittest discover -s tests -v
    python scripts/run_synthetic_baselines.py --replicates 400 --horizon 8 --output results/synthetic_baselines_v2.json

## Findings conditional on stipulated assumptions

### 1. Geographic arrival order and latent pathway loss
Suppose a retained white pathway degrades each generation with probability delta_source=0.02 at its source and delta_island=0.20 after arriving on an island. Re-expression disabled during the first delay generations. Then both paths face eight island generations with per-competent-lineage effective establishment probability r=0.1, and eight independent original W+ founders. The closed-form chance that **at least one** re-expression event establishes is:

| Delay generations | Early migration | Late migration |
|---:|---:|---:|
| 0 | 0.9510 | 0.9510 |
| 5 | 0.5806 | 0.9309 |
| 10 | 0.2400 | 0.9068 |
| 20 | 0.0286 | 0.8478 |

Formula per competent founder after arrival: p_gain = r / (r+delta_island) * [1-(1-r-delta_island)^H]. Early vs late remaining competent probability: (1-delta_island)^T vs (1-delta_source)^T. Probability for n founders is 1-[1-q*p_gain]^n.

**Logical boundary:** these are mathematical consequences of user-specified rates, not an estimated evolutionary reversibility window, a universal migration effect or evidence about ancestral Cirsium. If delta_source equals delta_island, the ordering effect vanishes in this simplified model.

### 2. Geographically identical recolouring need not have a common history
In four separate 400-replicate, eight-generation *toy* scenarios, the finite-population simulator distinguishes offspring recruited by:
- intrinsic re-expression from W+ (12 initial competent-white island plants; r=.08; delta=.06);
- inability to intrinsically reactivate W- without repair (12 initial degraded-white island plants);
- already-pigmented immigrants from source through seed movement (seed migration=.06; r=0);
- pigmented offspring of white maternal ancestry after cross-region pollen introgression (pollen delivery=.5; introgression probability=.15; r=0).

Observed coloured occupancy by the final generation was 0.995, 0, 0.995 and 1.000 respectively. Route labels arise from simulated ancestry/events, **not** from a field photograph. These scenarios were not parameter-fitted to equal occupancy. Both seed immigration and pollen introgression yield island colour in models with **zero intrinsic re-expression**.

### 3. Exact visible-state non-identifiability, possible rescue with diagnostic information
Construct three histories with equal prior weights and the **same** Bernoulli(0.25) observed coloured endpoint in all histories. Therefore with visible colour alone, the true history has no identifying signal: balanced accuracy is exactly 1/3 in a balanced design.

Introduce two hypothetical independent binary laboratory markers:
A = foreign genome-wide ancestry; B = foreign colour-linked haplotype.
Specify truth profiles intrinsic=(0,0), coloured immigration=(1,1), pollen introgression=(0,1). Assign both assays sensitivity 0.90 and specificity 0.95. Estimate assay emissions on 300 labelled synthetic samples per history, then classify 1,500 untouched synthetic samples per history.

Observed test results (seed=20261017):
- colour-only: 1500 / 4500 = **33.33%**;
- two hypothetical diagnostics: 3986 / 4500 = **88.58%**.

Sensitivity analysis (same random seed):
- sens/spec 0.50/0.50: **33.73%** (near chance, as expected);
- sens/spec 0.65/0.65: **49.13%**;
- sens/spec 0.80/0.80: **68.91%**;
- sens/spec 0.90/0.95: **88.58%**;
- sens/spec 1.00/1.00: **100%**.

**This benchmark is intentionally favourable:** real local ancestry and pigment haplotypes may not be diagnostic, known training histories rarely exist, multiple histories may coexist, and gene flow may erase simple ancestry markers. Hence it demonstrates what *sort of information* would be needed; it does not claim an achievable experimental classifier, validate a causal marker or identify an empirical re-expression event.

### 4. Competing first-recolouring routes
In an exploratory mixture with seed-migration probability .005, within-W+ re-expression probability .005, pollen-delivery probability .4, introgressive-reacquisition probability .005 and island competence degradation .08, the island starts with 8 W+ plants and receives seeds/pollen from a source with 8 P plants. Fecundity is 1.6 and the horizon is 8 generations.

Among 400 synthetic seeds of randomness, first P establishment was assigned to **pigmented immigrants 163**, **intrinsic re-expression 119**, **pollen-associated introgressive reacquisition 64**, **multiple same-generation causes 20**, and **no appearance 34**. All these are oracle simulator route labels, not categories inferable from visible flowers. Simultaneous events cannot legitimately be forced to one cause. The parameter regime is exploratory and is not an empirical calibrated frequency distribution.

## Model scope and caveats
The simulator has finite regional carrying capacities, stochastic propagules, seed movement, local recruitment, catastrophe, optional pollen-origin genetic rescue, and separate observer colour vs oracle-only latent states. It is not a coalescent, diploid genotype model, phylogenetic diversification/DEC model, detailed Mendelian hybrid simulator, seed-bank model or real-world dispersal kernel. Each simulated parental lineage carries a simplistic maternal tag; genotype ancestry should **not** be reconstructed from that tag in nature. Wilson intervals over Monte Carlo draws quantify simulator randomness, not parameter/model uncertainty.

## Gates to avoid an empty "theory result"
1. New empirical same-individual, organ-scoped molecular data independently measure pigment pathways: needed for W+ vs W-.
2. Realistic founder/seed and pollen gene-flow operators: estimate from independent data, never assign from observed colour.
3. Independent dated lineage/region histories with uncertainty: distinguish old variants, new mutation and reticulation.
4. Discriminative information improves **whole-clade held-out** prediction over matched visible/geographic/ancestry baselines; quantify calibration.
5. If those gates fail, report mechanism not identifiable rather than claiming a history-dependent accessibility effect.

The non-commutation of movement and degradation is a *property of specified stochastic operators*, not a new mathematical theorem by itself. Its usefulness depends on empirical measurements that break otherwise observationally equivalent histories.
