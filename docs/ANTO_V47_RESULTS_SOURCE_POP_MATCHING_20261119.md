# ANTO v47 results — Phlox leaf trait predictive transfer: does deamii's three-population sample count explain the negative effect?

2026-11-19. Original 2024 Phlox source-author field dataset was already SHA-verified/reanalyzed in ANTO v46. This follow-up is a **post-outcome, retrospective influence/reference analysis on the exact original v46 20 source-population out-of-sample MSE receipt**, not a new trained model, original C++ run, independent natural study or valid permutation test. Source Git blob `f03b0eb0177f807317ee1f8a7c0e65ee408d903b`; 2024 author original CSV blob `4c780b0abe644ac69dae86a0450eb598c160626b`.

## Original population grains and the frozen v46 result

Same three original common gardens, 963 clonal ramets, 284 original genotypes, 20 real source populations: amoena 8, pilosa 9, deamii 3. The v46 primary fruit prediction source-population holdout comparison produced an equal-source population MSE improvement **+9.39%, +9.54%, +9.43%, +10.63%** for G/T/E/TE PC slope models, respectively, over nine original garden × taxon mean baselines. The within-deamii result was **negative in every model** and v46's preregistered all-taxon **PASS_SOURCE_POP_TRANSFER** rule failed across all four models. This failure remains unchanged in v47. No newly measured seed paternity, floral pigmentation, visitor or descendant genotype outcome.

## Exact size-three source-population reference

For each existing model and endpoint, subset the already-heldout **population-level MSE** records into all 3-population subsets from the eight original amoena populations (**56** overlapping subsets) and nine original pilosa populations (**84**), for a descriptive reference of **140**. Compute each triplet's fractional MSE gain against the same parent v46 baseline (sum equally over three source populations). Compare with the ONE original deamii 3-population group. Do not treat 140 overlapping subsets as independent experiments or a randomization null. As a secondary, different comparison, combine the 17 original non-deamii populations into C(17,3)=680 (mixed-taxon triples are allowed and not exchangeable with single-taxon deamii).

**Primary original fruit endpoint:**

| v46 predictor | Deamii mean reduction in MSE | Out of 140 non-deamii same-taxon 3-pop subsets scoring as poorly or worse | Negative-reference subsets, irrespective of magnitude | Full 2024 all-taxon transfer gate |
| --- | ---: | ---: | ---: | --- |
| G, shared leaf-PC slopes | **−9.83%** | **0/140** | 17/140 | **NOT ESTABLISHED** |
| T, taxon-specific slopes | **−3.99%** | **8/140** | 25/140 | **NOT ESTABLISHED** |
| E, garden-specific slopes | **−9.25%** | **0/140** | 15/140 | **NOT ESTABLISHED** |
| TE, garden × taxon slopes | **−2.80%** | **11/140** | 30/140 | **NOT ESTABLISHED** |

There are **three original deamii source populations**, IDs 1719/1720/1721. All three separately worsen heldout **fruit** MSE under G, T and E, and two of three worsen under TE. Each of the three two-population reaggregations formed by deleting one deamii population remains negative for **ALL FOUR** fruit models (arithmetic only, no refit). Thus no single deamii population uniquely generates the sign of the original within-taxon negative fruit result; that result still has unknown sampling uncertainty because n=3.

Secondary original **flowers_2019** endpoint remains negative for deamii: G −12.54%, T −1.13%, E −11.19%, TE −1.00%; same-taxon reference counts at or below deamii are respectively **0/140, 21/140, 0/140, 21/140**. These are not a second independent ecological experiment; they reuse the same source plants and fitted model families. Baseline log1p biomass remains UNMEASURED by v46 in this criterion because 13 raw signed negative biomass rows made the predeclared log1p endpoint inadmissible.

## Precise ecological interpretation and severe limits

The source-specific *population-count-only* explanation is insufficient as a complete description. Merely selecting three evaluation source populations from amoena or pilosa does **not** commonly recreate the negative deamii aggregate, particularly for shared G and garden E slopes. Therefore the deamii nontransfer pattern is **not just an unavoidable algebraic consequence of having three held-out populations**. But this does not isolate a true biological taxon interaction, nor does it rule out chance, unequal original trait range, genotype composition, garden representation, or population baseline loss. Within-source comparisons are correlated and selected after knowing deamii's v46 result, so **0/140 is NOT a p-value, scientific discovery threshold or independent replication**.

Training regimes are not sample-size matched; reference grouping only matches the number of evaluated source populations, not their traits, environment or training populations. When permitting taxon × garden slopes (TE), negative size-three original references occur more frequently and 11/140 are at or below deamii. Thus the strong descriptive contrast has **model-dependence**. This is consistent with ecological context limiting transfer of trait information, but there is not enough independent biological sampling to say precisely what causes that limitation.

## What is genuinely new

- Exact enumeration of *evaluation-source count matched* reference scores across the already-validated original-v46 predictions, retaining all unfavorable models and both previous outcomes.
- Independent confirmation that all three original deamii source populations contribute to the shared- and garden-model negative fruit prediction; no single original deamii population dominates the sign.
- A clear failure mode for extrapolating global average predictive gains to phylogenetic taxa: a model can improve most source populations and still worsen **all three original source populations of another taxon**.
- Explicit separation of informative **diagnostic** from **inferential test**: no new p-value, no cross-clade macroevolutionary rate, no phylogenetic species holdout, no flower-colour trait reversal, no biological mechanism identification.

## Reproduction

- Frozen posthoc protocol and limited scope: docs/ANTO_V47_THREE_POPULATION_POSTHOC_PROTOCOL_20261119.md
- Source identity registry: data/external/anto_v47_three_pop_posthoc_registry_20261119.json
- Standard-library reaggregator and strict v46 parent blob identity check: scripts/analyze_Phlox_three_pop_conditional_v47.py
- Adversarial 17-test suite: tests/test_Phlox_three_pop_conditional_v47.py
- Derived complete reference and all deamii 2-pop arithmetic omissions: results/anto_v47_three_pop_reference_20261119.json
- Read-only CI: .github/workflows/anto-v47-three-pop-posthoc-audit.yml

**Keep Draft / stacked on ANTO v46 PR #56; do not merge to main automatically.** The next genuinely new data should add external **deamii** populations and independent garden conditions, ideally with direct measured floral color/pollinator/seed-genotype transmission endpoints. Repeating reshuffles of these same 20 original source populations cannot turn a descriptive pattern into a universal rule.
