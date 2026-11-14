# ANTO v31 — single-species holdout, phylogenetic isolation and pigment prediction

**2026-11-14. Exploratory continuation after v30; analytical choices below are fixed before viewing v31 distance/error outcomes.**

## Purpose and approved scope

Use the same original Wheeler et al. (2023), DOI 10.1098/rspb.2023.0275, processed floral RNA/HPLC measurements and authenticated dated Newick to ask whether a species' nearest-training-relative distance predicts its leave-one-species-out error. ANTO only. The primary panel contains all **59** measured Petunieae ingroup species; the secondary panel contains the existing **31** focal species shared with Pezzi et al. (2024). No further taxon selection, hyperparameter tuning, new chemical observations or ancestral-regain claim.

The user approved the species-LOO continuation of v30. Unlike complete-clade exclusions, one-species exclusions permit distances to vary across test species. They **do not** randomize phylogenetic distance or make sister species independent. Near-ultrametric trees can still create tied nearest-relative distances. v30's within-clade non-identification remains unchanged.

## Frozen source/model contract

Parent: `967dfb39c048301612f291d102309720f9ba31e0`, v30 Draft PR #40. Original v26–v30 sources and readouts remain immutable. Before accepting the new result, replay v30 from original pinned OSF and published-tree files and compare its entire corrected JSON exactly, thereby retaining the previous replay chain. Cache downloads only within the run. Source failures are technical failures, never ecological zeros.

Use the six existing v26 ridge models, lambda 10, training-only gene standardization and intercept, and two unchanged targets: log1p anthocyanidin total and log1p flavonol total. Primary molecular model remains F3primeH/F35H/FLS. Each test species is absent from all fitting, baseline calculations and neighbor selection. At 59 species, all six LOO RMSE/MAE values must reproduce the frozen v26 LOO results at their original rounding precision.

## Distance and comparator definitions

1. Nearest-training patristic distance uses original branch lengths before pruning or topology-only rerooting. Missing/invalid lengths fail. Independently reconstruct paths using Decimal representations of recorded edge lengths and compare all distances with the existing v30 path implementation. Ignore the root's unconnected stem. Units are source-Newick units, **not validated years**.
2. Preserve near-equal distances rather than interpreting rounding artifacts. Primary absolute tolerance is **1e-5 source branch units**; mandatory sensitivity tolerances are **1e-6 and 1e-4**, motivated by the ~1e-6 within-clade spread found in v30. Sorted complete-link bands have maximum-minus-minimum <= tolerance; no chain-merging or outcome-based grouping. Report both raw and tolerance-aware ranks.
3. Molecular novelty is the nearest training-species Euclidean distance in the three branch genes, standardized from training species only, as in v30.
4. Add a deliberately simple, fixed **phylogenetic-nearest-neighbor-only** predictor: the mean of the measured log-transformed target among all training species whose distance is within tolerance of the minimum. No test target is used; no arbitrary single neighbor is selected among ties. This is a secondary comparator, not PGLS or causal adjustment. Repeat it under all three tolerances.

## Pre-outcome fixed summaries

For each panel, keep every model's paired OOF RMSE and MAE, branch-model gain versus the training mean and the nearest-relative-only predictor, and the number of species whose squared error improves. For the primary branch model report:

- Spearman correlation between nearest-relative gap and squared error, training-mean squared error, and training-mean-minus-branch squared-error gain; both raw and tolerance-aware distance ranks. Report analogous molecular-novelty correlations.
- Distance-band summaries with all tied species kept together. Divide whole bands into near/middle/far using each band's midpoint in cumulative species count; never split a band by taxon alphabet or outcome to force equal sizes. Empty categories are not estimable. Report counts, ranges, null/branch/nearest-relative RMSE and gains.
- Equal-distance-band correlation between mean distance and mean squared error as a descriptive sensitivity, not an independent-clade test.
- Delete-one-**already-scored**-species influence ranges for distance–error rank correlation and branch-versus-neighbor RMSE gain. Do not advertise these as refitted CV, CIs or independent resampling.

The exploratory distance-only pattern is recorded separately for each endpoint: positive tie-aware rho **and** far-band RMSE greater than near-band RMSE. Joint generalization requires both endpoints in both panels; report failures and undefined summaries. Comparator superiority is assessed by lower whole-panel RMSE and MAE; any disagreement is retained. No p-values, posterior probabilities, causal distance slope, independent evolutionary event count or pooling across panels.

## Implementation plan and validation

Use a standalone standard-library module with lazy imports of the pinned parent analysis, a unit-test file, a dedicated workflow and aggregate JSON. Work on the isolated v31 branch, not main or v30. Sequence: (1) write failing tests for distance ties, target leakage, neighbor averaging and branch distances; (2) implement helpers and verify synthetic tests; (3) source-authenticated CI replays v30 and computes v31; (4) save the exact derived receipt and interpretation; (5) independently replay the frozen receipt and leave a Draft PR. Review missing labels/lengths, nonfinite data, duplicated IDs, zero denominators, empty bins and uninformative distances. Previous source outputs are never rewritten to make tests pass.

All original source RNA/HPLC rows and Newick bytes remain outside the repository; only derived metrics, source digests and fold/distance counts are retained. This is a single already-inspected original chemical study, not independent cross-study or future-fitness validation.
