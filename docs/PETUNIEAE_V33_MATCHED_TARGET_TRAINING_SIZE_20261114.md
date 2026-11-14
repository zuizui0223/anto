# ANTO v33 — same target, same training count: clade exclusions versus random omissions

**2026-11-14. Exploratory follow-up after inspecting v32. Design fixed before new matched-control outcomes; not a new biological preregistration or independent chemical study.**

## Question and unchanged evidence

The v32 phylogeny-plus-expression increment was small/mixed in species LOO but positive for whole-clade exclusions. Clade exclusions also remove more training species. Test whether the larger increment survives comparison with **the very same evaluation species and the very same number of training species**, when omitted training taxa are chosen randomly rather than by clade membership.

Parent `9185e9f2766df62df1a94b668f79841756263a38`, Draft PR #42. All v26–v32 sources, model functions and previous scientific conclusions remain unchanged. In particular, v32's primary species-LOO success criterion failed. It must not be replaced by v33's outcome.

## Fixed design

Use all 59 original Wheeler et al. 2023 Petunieae floral RNA/HPLC species means, DOI 10.1098/rspb.2023.0275, exact authenticated dated tree, and the original v27/v32 maximal monophyletic subtree caps **8 and 15**. Every original target is retained. Chemical endpoints remain log1p anthocyanidin total and log1p flavonol total. No new source observations are added.

For target species i in original heldout clade C of size m:

- **Real clade arm:** predict i using all 59−m species outside C, exactly as in v32.
- **Matched control arm:** exclude i, then sample exactly 59−m training species without replacement from the other 58. Other members of C are eligible. There is no constraint using gene expression, pigment or model performance. All unused taxa are neither train nor test for this one target. This keeps target identity and training count fixed, but changes training composition.
- Run **31 deterministic control repetitions**, seed base **20261114**, with independent per-target/per-cap/per-repetition seeds from SHA256 of the literal key. Keep all 31 draws; do not select favorable draws. Singleton real folds have identical full remaining training sets in every control and must yield zero real-minus-control contrast.
- At each real/control training set, select covariance nugget among **0.01, 0.1, 1, 10** by the unchanged v32 phylogeny-only inner species LOO. P and P+G share that selected nugget. Fixed gene ridge penalty **10**, fixed genes **F3primeH/F35H/FLS**, original shared-path covariance, train-only standardization. Only original training outcomes enter fitting; outer target outcomes are used solely for scoring.
- Also retain train-mean and genes-only comparator forecasts. No model tuning or additional predictor selection.

## Primary estimand and checks

For each original target, let g = squared error(P) − squared error(P+G). Report **D_i = g_real_clade − mean_r(g_matched_control)**. The main quantity is the species-average D_i and the equal-original-clade average D_i, separately for each endpoint and cap. Positive D means the molecular increment is larger when the target's real clade is absent than when equally many arbitrary species are absent. It is **not** a gene's causal effect or evolutionary selection coefficient.

Report the decomposition D = [MSE(P_real)−MSE(P_control)] − [MSE(P+G_real)−MSE(P+G_control)] and verify its arithmetic. The control MSE pools all target×31 scored errors; its RMSE is sqrt(mean squared error), **not** the mean of replicate RMSEs. Keep per-repetition aggregate scores and all original-clade signed contrasts. Also report RMSE/MAE of mean, genes-only, P and P+G in each arm.

Directional diagnostic: all four endpoint×cap species-average contrasts positive, with a majority of original clades positive in each cell. If any fails, report mixed/failed diagnostic. This is descriptive and exploratory, not a permutation p-value or independent replication. Show delete-one-already-scored-clade influence ranges (no refitting, confidence interval or bootstrap claim).

Verify the matching directly: all training sets unique within a draw, target never used, exact original train count per target, exactly 59 scored targets per repetition, no outcome-dependent seeds, and real-clade predictions reproduce all v32 main RMSE/MAE values. Track the nearest-training patristic gap for real/control arms and how often original clade peers are present in controls. Do not interpret a change in gap as the only change: training chemistry, gene distribution and species composition change together.

## Provenance and validation plan

1. Separate branch; main and parent code untouched. Implement new source runner, adverse tests, read-only workflow and derived result/receipt.
2. Tests before implementation: deterministic outcome-blind control membership, exact training sizes, singleton identity, invalid/overlapping IDs, heldout-target response mutations, v32 prediction equivalence, pooled-risk arithmetic and paired-effect decomposition. Observe failed tests before implementing.
3. Re-fetch original sources only through their official metadata links, verify original sizes/MD5/SHA256, replay v32 and compare its complete output digest to the committed receipt. Network failures remain technical failures. Never substitute simulated responses or cached summary metrics for new fitting.
4. Run the full repository test suite and actual source analysis. Preserve all negative arms and numerical settings. Save aggregate results, source/run IDs and full SHA256; never publish original measurements or tree bytes.
5. Re-evaluate the frozen result and leave the PR Draft/unmerged pending scientific review.

## Alternatives and limits

A new random fold partition alone would change the target composition in each comparison; per-target controls are preferred because target and training count match exactly. A learning curve alone would confound both counts and composition. Neither replaces independent RNA/HPLC evidence.

31 controls are repeated training designs on one measured panel, not 31 studies or new biological replicates. Species and fits remain dependent. No inferential p-values, confidence interval or universal phylogenetic law. A positive D rules out training **count alone** for this controlled computational contrast; it does not isolate ancestry from every other property of which taxa were withheld. One topology and one original chemical study remain the evidence ceiling.
