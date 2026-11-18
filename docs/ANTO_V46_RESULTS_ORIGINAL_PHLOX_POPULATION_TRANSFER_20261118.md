# ANTO v46 — Measured 2024 Phlox fitness: positive average leaf-trait information but no universal across-taxon transfer

2026-11-18; original-source reanalysis on a NEW independent 2024 common-garden cohort, NOT the 2014 PhloxSelection example input, NOT flower-colour experiments, and NOT a claim of heritable trait regain. All five models and the strict fruit decision were frozen in docs/ANTO_V46_EXTERNAL_PHLOX_POPULATION_HOLDOUT_PROTOCOL_20261118.md BEFORE inspecting prediction outcomes.

## Original source and biological unit

Goulet-Scott et al. (2024) Nature Communications, DOI 10.1038/s41467-024-49075-6; original authors' compiled field CSV at PhloxHopkins/PhloxFieldAdaptation pinned Git SHA d792fe58cce36721394adfd97e9f412396578e86, Git blob SHA-1 4c780b0abe644ac69dae86a0450eb598c160626b. The 963 RAMET rows belong to 284 genotypes, 20 source populations (amoena=8, pilosa=9, deamii=3), sampled in THREE original common gardens (321 rows each). All repeated clones for every test source population remain entirely withheld across all gardens. Population is the holdout unit and original common-garden site is not new.

The original authors reported local adaptation and differences in the trait axes shaping within- versus between-taxon fitness; those original published findings are NOT independent rediscoveries or new taxon-specific phenotypes from this ANTO analysis. The fitted predictors here are greenhouse-assayed LEAF PC1 and PC2, not flower colour, floral gene expression, pollinator visits or natural allele transmission.

## Frozen analysis and new heldout figures

All 20 leave-source-population-out folds were scored on log(1+fruits). Every model has nine separate training-only garden x taxon group mean baselines, with a fixed alpha=10 ridge correction on PC1/PC2: common slopes G, taxon slopes T, garden slopes E, garden-by-taxon slopes TE. There is NO response-based hyperparameter tuning. Baseline B0 equal-population MSE=2.20305469 and record-weighted MSE=2.07599011. All losses are from original nonmissing, unshifted fruit counts in 963 ramets.

| Model | Heldout equal-population fruit RMSE-independent MSE gain vs B0 | Heldout ramet-weighted MSE gain | Source populations with improvement | Amoena 8-pop gain | Pilosa 9-pop gain | Deamii 3-pop gain | Preregistered all-taxon decision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Common PC slopes G | **+9.39%** | **+4.87%** | 13/20 | +9.93% | +13.75% | **-9.83%** | **NOT ESTABLISHED** |
| Taxon-specific slopes T | **+9.54%** | **+5.23%** | 13/20 | +9.56% | +12.84% | **-3.99%** | **NOT ESTABLISHED** |
| Garden-specific slopes E | **+9.43%** | **+4.89%** | 13/20 | +8.98% | +14.28% | **-9.25%** | **NOT ESTABLISHED** |
| Full garden x taxon slopes TE | **+10.63%** | **+5.98%** | 13/20 | +8.99% | +14.92% | **-2.80%** | **NOT ESTABLISHED** |

The strict source-population transfer gate requires gains >2% in both weighting schemes, at least 12/20 improving heldout populations, and positive heldout equal-population gain within EACH of amoena, pilosa and deamii. All four tested PC models pass the first two requirements, and ALL FOUR fail the last one due to deamii. The deamii direction remains negative even when permitting 18 separate group-specific slopes; but three source populations are too few to infer a general negative evolutionary process or true taxon-level null.

We intentionally do NOT re-select the target model or original taxa and do not manufacture a favorable cross-taxon transfer claim. The strongest defensible positive is within-source heldout fruit information **for amoena and pilosa**; the overall mean predictive gain is not universally portable.

## Secondary original outcome and data-quality nonadmission

Prespecified secondary flowers_2019 log1p heldout MSE gains are G +11.65%, T +11.98%, E +11.90%, TE +11.53% on population-equal weighting. Deamii is negative in all four flower models, too. These are SECONDARY and do not rescue the frozen primary fruit criterion.

Prespecified log1p biomass was **blocked** after a source-domain check: 13/963 original biomass records have negative values (minimum **-233**), which are not admissible to log1p. No imputation, deletion, clipping, artificial offset or post hoc change of transform. The original file remains intact; our derived JSON explicitly records the numeric technical nonadmission and contains no biomass model scores. These source values alone do not demonstrate an author data error: corrected biomass measures may be signed, and the original paper analyzed its own endpoints.

## Integration into ANTO × SDEA × EGWE × izu-core

The original 2014 PhloxSelection source-C++ v44/v45 comparison was a fixed unoptimized model-likelihood sensitivity across FIVE spatial transects. In contrast, v46 is a different actual 2024 study: repeated clonal RAMETS, THREE gardens, source-population-heldout original fitness-related observations, greenhouse leaf PCs. Their numerical likelihood deltas and log-transformed MSE percentages cannot be pooled.

**Convergent limited inference:** selected trait-information layers can add observable predictive information while ecological/taxonomic context limits where the gain transports. This does not prove that genetic or floral mechanisms are universally environment-gated. Stage-wise compatibility, real pollination, same-genotype seedset and later-generation allele retention remain unlinked across studies. No whole-chain fitness rule or trait-regain probability has been estimated.

## Reproduction, source restrictions and review

- Frozen prior design: docs/ANTO_V46_EXTERNAL_PHLOX_POPULATION_HOLDOUT_PROTOCOL_20261118.md
- Pinned full original source manifest: data/external/anto_v46_phlox2024_source_registry_20261118.json
- Hash-verified original remote CSV fetch, 20 complete-population holdouts, nine group baselines, train-only feature standardization and ridge solver: scripts/analyze_original_Phlox_2024_population_holdout_v46.py
- Adversarial source/leakage/domain tests: tests/test_original_Phlox_2024_population_holdout_v46.py
- Derived-only fold record: results/anto_v46_phlox2024_population_heldout_20261118.json
- Dedicated read-only source rerun: .github/workflows/anto-v46-original-phlox-population-transfer.yml

Draft / stacked on v45 PR #55. Keep main intact pending original-source exact CI, scientific review and clearer independent garden/taxon transfer. Next actual biological test requires independently sampled gardens or populations from other clades and actual floral/pollinator/seed measurements, not more resplitting of 2024 data.