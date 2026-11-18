# ANTO v46 — Prereadout protocol: cross-population portability of leaf-trait fitness information

Frozen 2026-11-18 before inspecting the population-heldout predictive metrics. Base ANTO v45, not main. This is a *different 2024 original Phlox cohort* with three other taxa, independent of the 2014 drummondii example input; it is NOT an independent phylogenetic radiation nor a flower-color experiment.

## Authenticated original source and hard gates
Original Goulet-Scott et al. 2024 Nature Communications DOI 10.1038/s41467-024-49075-6 and Dryad DOI 10.5061/dryad.gxd2547sx; original author repo PhloxHopkins/PhloxFieldAdaptation pinned commit d792fe58cce36721394adfd97e9f412396578e86, original author file Dataset_fitness_leaftraits_distance_lat_long.csv, Git blob SHA-1 4c780b0abe644ac69dae86a0450eb598c160626b, 175377 bytes. Source identity must be verified before processing; no raw author observations redistributed into ANTO.

Exactly 963 original plant/ramet rows, 26 CSV columns, 284 repeated genotype IDs, 20 source population IDs, 3 gardens of 321 rows each. Original taxa: amoena 366 rows and 8 populations, pilosa 375 rows and 9 populations, deamii 222 rows and 3 populations. Each garden-by-taxon contains 122 amoena, 125 pilosa or 74 deamii ramets. Population IDs 1729_amo and 1729_pil remain separate. Any changed original identity, missing numeric response, duplicate observation grain, incompatible garden/taxon/geno assignment or invalid cohort must FAIL CLOSED. The 2024 study measured greenhouse leaf PC axes and fitness-related common-garden traits; NOT flower colour, pollination channel, gene expression, adaptive mutations or descendant genotype persistence.

## Fixed outcomes, holdout and model families
Primary endpoint: log(1+fruits_2019). Secondary descriptive endpoints: log(1+flowers_2019) and log(1+biomass_2019). All outcomes and leaf PC1/PC2 must be finite and nonnegative where applicable. Exactly 20 source-population-heldout folds; withhold every ramet from each population in all 3 gardens. The 284 unique genotypes (each repeated across ramets) must be entirely confined to their source population, so no heldout clone leaks into training. These are *heldout populations within the same three original gardens*, not a leave-garden-out test or a replicated new study.

For each training fold, the 9 disjoint garden x taxon cells have separately fitted training-only outcome means (B0 reference). Center both leaf PCs and outcomes using the means of their respective training-only garden x taxon cells. Standardize the centered PC features using pooled training-only standard deviations (population variance, zero-variance error). Ridge alpha=10 fixed in advance, penalizing only slopes. Solve linear ridge normal equations by a pivoted linear system; never use heldout outcomes for fitting, normalization, penalty choice or group means.

Candidate models: B0 = 9 training-only garden x taxon means; G = B0 + one pair of shared PC1/PC2 slopes; T = B0 + slopes unique to each of 3 taxa (6 slopes); E = B0 + slopes unique to each of 3 gardens (6 slopes); TE = B0 + slopes unique to each of 9 garden x taxon cells (18 slopes). Evaluate all models regardless of outcome sign; no post hoc tuning.

Each model: store source-population heldout squared-error sums and counts, equal-population mean MSE and record-weighted mean MSE, number of populations with strict positive prediction gain vs B0. Gain = 1 - MSE_model/MSE_B0. Do not treat the 963 ramets, 284 cloned genotypes or 20 source populations as independent original evolutionary studies, and do not pool endpoints.

## Falsifiable predeclared primary rule
For FRUIT only, classify each model PASS_SOURCE_POP_TRANSFER if all hold:
- relative MSE gain > 2% on both equal-population and record-weighted heldout risk;
- positive MSE differences for at least 12 of 20 withheld source populations;
- strictly positive equal-population MSE gain within EACH of the three original taxa (8/9/3 populations).

Otherwise classify NOT_ESTABLISHED; preserve all negative folds/models and secondary endpoints without reselecting the winner. No nominal p-value, bootstrap confidence interval or universal ecological claim attached to this descriptive predictive gate. This is one 2024 original experiment, and the original authors previously reported multidimensional selection: any positive v46 result is incremental within-study predictive portability, not a new floral-color or across-family fitness law.

## Required reproducibility
Create a standard-library source SHA-verified script with hard cohort and train/test leakage checks, exact derived result JSON, adversarial tests, and a read-only GitHub Actions workflow with exact regeneration. Keep v46 Draft stacked on v45 and leave main intact until scientific review.