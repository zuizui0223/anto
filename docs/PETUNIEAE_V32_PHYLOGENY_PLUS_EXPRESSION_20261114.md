# ANTO v32 — does floral expression improve prediction conditional on a phylogenetic predictor?

**2026-11-14 | Exploratory follow-up after v31, design fixed before v32 results.**

## Decision and scope

V31 showed that the fixed three-gene branch module predicts current chemical totals better than a simple nearest-relative copying rule. However that copying rule was itself worse than the training mean. A stronger question is whether the SAME phylogeny-based model improves when the three expression measurements are added. Use only the original authenticated Wheeler 2023 processed species-mean RNA/HPLC and original 11-gene dated tree; do not call this a new experiment, genetic causality, ancestral trait regain or independent replication.

Parent commit: `45154b7fb4ad653e6c380afbc2a707df40c9a6fa` (v31 PR #41). No changes to v26–v31. The primary evaluation is all **59** species leave-one-species-out; the existing **31** species panel is secondary. Additional 59-species complete-clade exclusions use the original v27 **8/15** maximum-tip cuts. Report all four designs, including negative results. No new complete-case filtering.

## Fixed models and genuinely training-only selection

The original source tree defines Brownian-style covariance K(i,j) as the sum of branch lengths shared on root-to-tip paths. The root's unconnected stem is ignored. K is a modelling assumption, not evidence that these chemical traits actually evolve by Brownian motion. Validate K symmetry/positive semidefiniteness and Kii+Kjj−2Kij against the separately implemented v31 patristic distances. Tree information may include heldout taxa, but their chemical outcomes are never used.

For each outer training set, h is its mean diagonal of K; use C = K_train/h + eta*I and cross-covariance c = K_test,train/h. Fixed candidate eta values are **0.01, 0.1, 1, 10**, in ascending tie-break order.

- **P (phylogeny only):** fitted unpenalized intercept and conditional residual prediction, i.e. generalized least-squares intercept followed by c C^-1(y_train−intercept).
- **P+G:** identical covariance and residual prediction, plus original log1p F3primeH/F35H/FLS expression features standardized ONLY within that outer training set. Coefficients minimize (y−Xb)'C^-1(y−Xb) + **10** times the sum of squared gene coefficients; intercept unpenalized. Predicted mean is X_test b + c C^-1(y_train−X_train b).
- **M (training mean)** and **G (the existing fixed v26 three-gene ridge model)** remain comparators; not substituted for P.

Choose eta using **inner leave-one-training-species-out error of P only**, separately for each chemical endpoint and each OUTER training set. Every inner fit recomputes h and its intercept using only its own training species. Then use that same eta in outer P and P+G. This deliberately holds the phylogenetic nuisance setting fixed between models; it does not tune P+G on heldout errors. Gene penalty remains 10, inherited from the earlier fixed design. Report selected-eta counts, inner risk and all four fixed-eta outer results as sensitivity analyses. Do not choose a favorable fixed eta after seeing outcomes.

Linear solves use Cholesky factorization, not explicit matrix inversion. Independent tests verify the equations against a direct small-matrix reference. No test target is passed to fitting or inner selection. No full-dataset expression scaling or future trait leakage.

## Endpoints and reporting

Targets remain log1p(sum original anthocyanidin compounds) and log1p(sum original flavonols). All 59-species LOO G and M metrics must reproduce the frozen v31 results. Report species-weighted RMSE/MAE for M/G/P/P+G, P+G versus P relative RMSE improvement and signed MSE gain, number of species with lower squared error, and per-outer-fold signed mean gain. Clade designs additionally report equal-clade mean MSE gain. Delete-one-already-scored-species influence ranges are descriptive only: no refit, CI, p-value or independence assumption.

The primary *descriptive* conditional-improvement criterion requires P+G to beat P in both RMSE and MAE for BOTH chemical endpoints in the **59-species LOO** design. External portability requires separate evidence; the 31 subset and whole-clade holds are not independent biological replications. All disagreements across endpoints, fixed eta sensitivities or clade exclusions are retained. Even a positive result shows only extra predictive information within this chosen covariance/model family, not independence from every aspect of ancestry or causal enzyme function.

## Source and verification contract

Use v31's official OSF metadata/download transport and verify the SAME four original file IDs, project/path, size, MD5 and SHA256. Cache original bodies only on the runner; never commit them or upload them as artifacts. Replay the entire frozen v31 record (which replays v26–v30) before new v32 calculations. Missing sources, unknown IDs, nonfinite/negative branch lengths, mismatched taxa, nonpositive covariance scales, test/train overlap, altered parent metrics or failed algebra tests STOP the analysis.

Implementation: a new standalone module, tests, dedicated read-only CI and derived aggregate JSON. TDD: write and run tests first, implement and re-run, then authenticate real sources, freeze aggregate output, and replay. Test original covariance paths; constant genes => P+G=P; no shared covariance => intercept/ordinary ridge limit; training-only eta and scaling; response mutation cannot alter heldout prediction; correct nested train indices; negative/duplicate IDs. Full repository unittest suite is run. No shared/main branch update or merge.

## Methodological sources

Felsenstein (1985), *Phylogenies and the Comparative Method*, DOI 10.1086/284325, motivates accounting for phylogenetic nonindependence under a specified Brownian model, not treating species as independent evolutionary replicates. The shared-path covariance formulation is also given in Freckleton & Jetz (2009), *Space versus phylogeny*, DOI 10.1098/rspb.2008.0905. Numerical implementation uses SciPy official `scipy.linalg.cho_factor` and `cho_solve` documentation. The particular penalized conditional predictor, four-value nuisance grid and shared-nuisance comparison above are our explicit exploratory design, not an original-paper result.

**No v32 real-data result has been inspected at this design commit.**
