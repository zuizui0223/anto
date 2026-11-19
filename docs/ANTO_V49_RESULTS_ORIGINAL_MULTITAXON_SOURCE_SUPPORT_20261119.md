# ANTO v49 — Focal-population count matched, full v46 multi-taxon model: deamii nontransfer depends on sharing architecture

**2026-11-19. Targeted POST-v46/v47/v48 diagnostic, not independent hypothesis confirmation.** The original 2024 Phlox author data are the same authenticated source throughout these analyses: Goulet-Scott et al. (2024) Nature Communications DOI 10.1038/s41467-024-49075-6, original PhloxHopkins/PhloxFieldAdaptation commit d792fe58cce36721394adfd97e9f412396578e86, source CSV blob SHA1 4c780b0abe644ac69dae86a0450eb598c160626b. **Original 963 ramet rows, 284 genotypes, 20 donor populations, three gardens and three taxa.** Original source rows are not redistributed; no new biological field data, independent gardens, flowers' colour, pollinator behaviour or descendant genotype fitness have been introduced.

## Why the previous conclusion needed challenge

In v46, garden×taxon group-mean fruit baseline B0 was improved by shared flower-unrelated **leaf** PC1/PC2 slope models in overall 20-source-population holdouts but NOT in each taxon. For deamii, source-population-mean fruit MSE worsened under all four models, so v46's pre-outcome all-taxon transfer gate **FAILED**.

In v47, evaluation-source count only was matched. Deamii's full-v46 G/E source-population score was more negative than **all** 140 overlapping three-population amoena/pilosa subgroup aggregate references. This did not match actual training availability.

In v48, both evaluated and trained same-taxon source counts were matched (two original same-taxon training populations + one heldout) but the training removed other taxa and reduced the model from **nine garden×taxon intercept cells** to only **three garden cells**. Under this changed fitting domain, the deamii G/E negative-score ranking was no longer exceptional (34/140 and 44/140 references at or below deamii). This could not isolate population support from inclusion of other taxa in the coefficient fit.

## v49 new source-original model experiment

Re-fit **423 full source-population-heldout folds**, all original three gardens, retaining exact original v46 unmodified Python **fit_train** and **prediction** methods for B0/G/T/E/TE, nine garden×taxon training-only baselines, global train-only leaf-PC scaling, alpha=10 ridge penalty on slopes and original fruit/flower endpoints. For each chosen focal 3-population panel (C(8,3)=56 amoena, C(9,3)=84 pilosa, C(3,3)=1 deamii), train the model on:

1. **Two same-taxon source populations** from that focal 3-pop panel;
2. **All original source populations of the other two taxa**, with the third focal source population fully withheld across all three gardens.

This preserves other-taxon availability and all nine intercept groups and limits only the number of focal-taxon training populations compared to v46. However, removing other focal populations necessarily changes original training ramet counts and trait representation, not a unique causal "sample size only" contrast.

### Exact original deamii and baseline consistency checks

Each deamii heldout fold trains on the same **19 original source populations** as in v46; all **3 populations × 5 original models × 2 response endpoints = 30** exact original-source predicted-MSE values match v46's previous fixed receipt (1e-7 check).

An additional AFTER-result v48-v49 baseline parity audit verifies **846 of 846** three-population-heldout B0 score cells (423 folds × original fruit and original flower counts) equal to within **5e-10** absolute, when each model is trained from two focal-taxon source populations. Thus v48-v49 reported differences are **NOT due to changed evaluation source populations or baseline intercept fit**. The differences arise in how leaf-PC slopes are fitted when other taxa are included: that includes sample composition, pooling of slopes and training-wide feature scaling/implicit ridge strength. This parity check was added post-outcome and is not disguised as preregistration.

## Primary fruit-number source-population-heldout result (log1p fruits)

| Unmodified v46 model | Deamii relative MSE gain | Non-deamii **140 overlapping** same-taxon triplet panels with gain <= deamii | Number of negative reference panels |
| --- | ---: | ---: | ---: |
| G — common leaf-PC slopes | **−9.83%** | **0/140** | 14/140 |
| T — taxon-specific slopes | **−3.99%** | 34/140 | 53/140 |
| E — garden-specific slopes | **−9.25%** | **0/140** | 17/140 |
| TE — garden × taxon slopes | **−2.80%** | 47/140 | 58/140 |

Original v48 same-taxon-only G/E from the same two-focal-pop scheme: deamii **−4.53%, −3.90%**, with **34/140 and 44/140** references as negative or more; v49's G/E difference does not reflect baseline prediction differences. It reflects changed fitted slopes, pooling support and model scaling.

**Important architecture distinction:** G and E use information across taxa in common coefficients and deamii's negative outcome is unusual relative to the original other-taxon reference; T and TE permit taxon-specific coefficients and the negative relative ranking is far less unusual. This is compatible with "negative transfer" of pooled trait–fitness predictive information when a taxon is sparsely represented, **but** is not proof that slope biology, pollinators or genes causally differ. It is also possible that the taxon-specific models are more heavily regularized differently due to training-global standardization.

## Secondary flower *count* endpoint, not flower colour

v49 source mean deamii log1p flower count gains: G **−12.54%** (0/140 as poor), T **−1.13%** (29/140), E **−11.19%** (0/140), TE **−1.00%** (39/140). Same original plants/field study, not independent biological replication. The original predeclared log1p biomass remains BLOCKED because 13 original recorded biomass values were negative (minimum −233); none clipped or imputed.

## What this DOES NOT demonstrate

- 140 highly overlapping triplet panels and 423 refits are NOT 140/423 independent new biological experiments, valid permutation tests, inferential p-values or confidence intervals.
- The deamii taxon has only THREE source populations and all studies retain the same three original gardens. No externally sampled original deamii population or fourth garden has been added.
- v46 global heldout success does **not** become a PASS for universal 3-taxon prediction; the original primary all-taxon decision remains **NOT_ESTABLISHED**.
- Cross-taxon G/E performance differences reflect predictor estimation *architecture*, training composition and feature normalization as a bundle, not a uniquely identified causal biological mechanism.
- No experimental flower-colour manipulation, actual pollination, independent seed-genotype inheritance, ancestral pigment pathway functional restoration or natural multi-generation persistence were assessed.

## An operational general-theory implication

Cross-taxon predictors can be **apparently informative in aggregate but systematically harmful** in a weakly represented biological domain. Such negative transfer must be distinguished from an intrinsic evolutionary mechanism by controlling: source biological grain, support and genotype overlaps, training baseline, model architecture, and out-of-sample endpoints. This is a methodology/identifiability finding in **one** original study, not a new universal ecological law or cross-radiation meta-analytic effect.

Next internally testable alternative: keep v49's exact same training panels, same baseline and fixed penalty, but compare training-only standardization or raw-PC ridge penalty between the shared and taxon-specific slope models. Only new independent field gardens and populations can identify genuinely transferable ecological mechanisms.

## Exact evidence and CI

- Before-v49-readout posthoc frozen scope: `docs/ANTO_V49_MULTI_TAXON_FOCAL_SUPPORT_PROTOCOL_20261119.md`
- Original 2024 immutable source / frozen v46 receipt manifest: `data/external/anto_v49_multi_taxon_focal2_registry_20261119.json`
- SHA-verified source, imported unchanged v46 fit_train/prediction, full 423 population holdouts and original deamii anchor: `scripts/analyze_Phlox_2024_multitaxon_focal2_v49.py`
- Fail-closed unit tests and explicit post-outcome v48 846/846 baseline-parity audit: `tests/test_Phlox_2024_multitaxon_focal2_v49.py`
- **All 141 panels, 423 heldout population model losses for both endpoints and all five models, without original source CSV rows**: `results/anto_v49_original_Phlox_multitaxon_focal2_20261119.json`
- Source-only Actions: `.github/workflows/anto-v49-original-multitaxon-focal2.yml`

Keep **Draft / stacked on v48 PR #58**. Main is not changed or merged. Source-derived readout and independent arithmetic cross-check complete; dedicated full Python source CI status must be checked and must PASS before scientific merge claims.
