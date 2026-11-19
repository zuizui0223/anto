# ANTO v50 — Original Phlox 2024 source: coefficient sharing versus ridge-coordinate sensitivity

Date 2026-11-19. **POSTHOC** targeted methods analysis: the new question was selected after v46–v49 scores were examined. This document and the source registry are frozen BEFORE any v50 scenario readout. It is not a preregistered independent biological study, p-value, new field dataset, phylogenetic replication, or evidence of ancestral pigment regain.

## Why

v49 preserved the exact 2024 original Phlox full multispecies training structure and original v46 fit, while restricting the focal taxon to two training source populations. For 423 source-population holdouts within 141 overlapping three-population panels, the original garden×taxon baseline (B0) was unchanged. Deamii's original fruit MSE gain was G −9.83%, T −3.99%, E −9.25%, TE −2.80%, with 0/140 same-size non-deamii panels as poor or worse for shared G/E but 34/140 (T), 47/140 (TE). This can reflect differences in learned shared versus separate slopes, but fixed original alpha=10 under training-global PC standardization may shrink different parameter blocks differently. The full 2024 v46 primary universal transfer criterion remains FAILED (NOT_ESTABLISHED).

## Original biological source and immutable v49 anchor

Source: Goulet-Scott et al. 2024 Nature Communications, original author Git commit d792fe58cce36721394adfd97e9f412396578e86 in PhloxHopkins/PhloxFieldAdaptation, exact Dataset_fitness_leaftraits_distance_lat_long.csv, 175377 bytes, original Git blob SHA1 4c780b0abe644ac69dae86a0450eb598c160626b. Original 963 clonal ramet rows; 284 original genotypes; original 20 populations (amoena8, pilosa9, deamii3); three original gardens. ANTO v49 frozen output results/anto_v49_original_Phlox_multitaxon_focal2_20261119.json, Git blob 94947d03e85766a3c5c1f82fde2bb076f2ca42e8. All heldout source-population folds and original genotype clone grouping remain exactly v49; no additional biological units.

Every panel is 3 original source populations of ONE focal taxon (C8,3=56 amoena, C9,3=84 pilosa, C3,3=1 deamii). For each heldout population, fit on the other TWO focal-taxonomic populations PLUS **all original nonfocal taxa** and all three original gardens. Exactly 141×3 = 423 folds. All training input row order stays as in original v46. Preserve nine unpenalized original garden×taxon response group intercepts, training-only within-cell leaf PC centering, original two leaf PC predictors, and original transformed log1p fruit primary and log1p flower NUMBER secondary endpoint. Original 2019 biomass transform is technically nonadmitted because 13 original raw signed values <0, min −233; no clipping or invented endpoint.

## Frozen design of the extra model coordinates

Four coefficient-sharing models from original v46:
- G: two global leaf PC slopes shared across taxa/gardens.
- T: two slopes for each taxon (6 slopes).
- E: two slopes per garden (6 slopes).
- TE: two slopes for each garden×taxon (18 slopes).

No separate model intercept is fitted per individual: all intercepts are the same nine original v46 group training means, and slopes regress outcome residuals against training-group-centered PCs. For each scenario, train-only PC SD is calculated globally from ALL original training-group-centered rows (population SD denominator n, not n−1). Each model's block-specific coefficients are solved by independent 2×2 normal equations, including the same v46 ridge loss. Require positive determinant, finite coefficients and no heldout outcome leakage. Use exactly the same frozen v49 population panels, original CSV bytes and response observations.

The **predeclared 7 sensitivity scenarios**:
1. STD_R0 — centered leaf PC divided by training-only global SD; no ridge penalty (OLS; lambda 0).
2. STD_R1 — same standardized PCs; lambda=1.
3. STD_R10 — same standardized PCs; lambda=10, EXACT original v46/v49 target. 423×4×2 per-fold model MSE values MUST reproduce the frozen v49 receipt within 1e−7, all 9 B0 groups and 30 deamii source fit anchors.
4. STD_R100 — same standardized PCs; lambda=100.
5. RAW_R0 — centered PCs in original author numeric units (no global SD); lambda=0. Must be prediction equivalent to STD_R0 (absolute per-population MSE 1e−7), because ordinary least squares with a nonsingular diagonal predictor transformation is equivariant.
6. RAW_R10 — same original centered PC units and same scalar lambda=10 on raw coefficients. THIS CHANGES EFFECTIVE PENALTY RELATIVE TO STD_R10 and is an exploratory coordinate-sensitive penalty comparison.
7. RAW_EQUIV_R10 — original centered PC units, but ridge penalty **per PC coefficient** lambda×(original global training SD for that PC)^2. Must reproduce STD_R10 exactly for every original heldout cell at 1e−7, because b_std = SD × b_raw. This is a numerical invariance control, NOT a new biological scenario.

All scenarios apply to G/T/E/TE using identical original source folds and garden×taxon B0 means. No hyperparameter optimizations, fold selection based on target outcomes, post-outcome lambda changes or new model families. Original v46 B0 baseline MSEs do not depend on the scenario.

## All outcomes, reference, scope and interpretation

For each of the 141 original three-population panels, every 3 source-population test fold, both endpoints and all four coefficient models × seven scenarios, retain full numerical source-derived per-fold MSE *or* a cryptographically frozen full machine-readable record. Compute panel-equal source MSE gain against the SAME original B0 for each model/scenario. Full 140 overlapping amoena/pilosa same-taxon reference comparison with original deamii one-panel gain: counts <= deamii, number of negative original panels, median/min/max per scenario, endpoint, and model. Additionally report paired within-deamii contrasts (T minus G; TE minus E) and whether their sign persists over lambda values. Do not choose a favorable lambda or treat 140 overlapping panels as independently drawn evolutionary studies, exchangeable null replicates, calibrated uncertainty or valid p-values.

Falsification:
- If negative shared-model G/E deamii transfer persists for STD_R0, it cannot be explained wholly by a scalar ridge penalty under fixed original standardized coordinates.
- If the shared-model loss disappears at STD_R0 or varies strongly with penalty/coordinates, that weakens an interpretation of intrinsically nonportable shared trait information.
- If STD_R10 and RAW_EQUIV_R10 differ, the numerical/parameterization implementation FAILED; do not interpret biology.
- If RAW_R0 and STD_R0 differ, same nonadmission applies.
- If v49 original MSE reproduced only after changing source or hyperparameter, fail.

These tests cannot uniquely identify biological species-specific selection slopes, functional genes, pollinator coupling, or adaptive isolation; coefficient sharing co-varies with model flexibility, variance and sampling support. Even a robust shared-model failure is an original 2024 Phlox prediction-domain result, **not a universal ecological law** or independent cross-radiation validation.

## Reproducibility and release

Create a new Python standard-library script importing unchanged v46/v49 source methods for validation and original group means; pin v49 derived parent blob; add independent 2x2 ridge solver with diagonal and RAW_EQUIV penalty, 7 regimes, all 423 fold MSE and full panel references, strict OLS/ridge equivariance and frozen original v49 MSE. Add adverse tests, source-gated CI and immutable derived-only receipt (never original author CSV rows). Keep Draft stacked on v49 PR #59; do not merge main without exact HEAD CI and independent scientific review. New data need independently sampled deamii source populations/gardens and directly linked flower–pollinator–pollen–offspring endpoints.