# ANTO v30 — phylogenetic isolation, molecular predictor novelty, and heldout pigment-prediction loss

**14 November 2026 — protocol frozen prior to v30 source-run outcomes; explicitly exploratory after observing v26–v29.**

## Scientific question

Previous source-authenticated 2023 Petunieae data reanalysis found that a fixed three-gene branch expression module (F3primeH/F35H/FLS) predicts *current* species-mean anthocyanidin and flavonol totals. In v29 on the 31 focal species overlapping published 2024 external trees, the fine (2–5-tip) holdouts gave only **+0.22–0.89%** flavonol RMSE gain, despite >~8% anthocyanidin gain. But **clade size, heldout species composition, predictor novelty, and taxonomic separation co-vary**. This study probes which *descriptive* factors track the error without claiming a causal evolutionary distance mechanism.

**Target:** Does a heldout species become harder to predict when its nearest training species is more phylogenetically separated, or when its three-gene expression vector lies further from available training species? Is the answer endpoint-specific? We will use observed *raw source tree patristic path lengths*, not simply count topology edges, genetic variation, geographic distance, time since trait loss, or absolute geological time.

## Frozen sources, denominators, and candidate predictors

1. **Wheeler et al. (2023)** source-verified original processed three Petunieae flower-bud expression/HPLC CSVs and 1,481-byte dated-pruned Newick tree from OSF zg9cu, strict registered size+MD5+SHA256. Only current RNA/HPLC species means, not independent biological RNA quantification or independently measured descendant fitness.
2. **Pezzi et al. (2024)** published actual ASTRAL/IQ-TREE tree + author individual-to-species map, pinned to exact original public GitHub commit/three Git-blob SHA1s. This is used **only to reproduce v28/v29 and to define the 31 focal species overlap**, not as an alternative patristic branch-length scale: between-tree branch lengths are not directly comparable.
3. Original v26, v27, v28 and v29 derived source receipts must **all exactly replay** before v30 results are accepted. No prior source rows/results are replaced or tuned.
4. **Primary matched panel:** 31 original 2023 chemistry/RNA species also included among 2024 published *focal* species, HUTE-rooted true original 2023 subtree folds with pre-existing max heldout sizes **2,3,4,5,6,8,15**. **Secondary panel:** all 59 original Petunieae ingroup chemistry/RNA species, source 2023 BROW-rooted original dated tree, old v27 max heldout sizes **8,15**.
5. Each heldout fold's training set includes **every other measured panel species**. Exact same **6 v26 fixed ridge models**, training-only standardization/intercept, lambda=10, original anthocyanidin/flavonol chemical-total log1p targets. No tuning on heldout outcomes. Each species is scored exactly once per panel × cap × endpoint, keeping every negative case.
6. Predictor fixed **before** outcomes:
   - **Patristic nearest-train gap** for each heldout species: minimum source-Newick sum of original branch lengths along undirected paths to any other species *in that fold's training set*. Must preserve all source edge lengths before any v28 re-root/pruning (those procedures intentionally discard branch lengths). Newick root's own length, if present, is an unconnected stem and not part of distances. Reject any missing/nonfinite/negative child length and mismatched source-tip identity; do NOT silently substitute hop count.
   - **Molecular novelty:** minimum Euclidean distance in **three pre-existing v26 branch gene log1p expression coordinates**, each standardized with **training-only** means/SDs for that specific fold (same no-test-leakage convention as ridge). This is an observational expression-distribution distance, not a genetic linkage, allele or heritability estimate.
   - **Reference forecast errors:** source-original per species out-of-fold squared error for fixed train-only mean and branch-three model, plus signed *loss improvement* = squared error(mean) − squared error(branch-three). All predictions must exactly reconstruct v27/v28/v29 headline RMSE for every evaluated source panel/cap/endpoint/model, including the five comparator modules.
7. For each panel × cap × endpoint, report **all**:
   - species-weighted original model RMSE, training-only-null RMSE and signed relative RMSE improvement;
   - **Spearman rank correlations** between nearest-train patristic gap, expression novelty, and each model's per-species squared error (and null-minus-branch squared-error gain). No p-values or independent observation assumption; source species and folds share phylogeny/training data.
   - fixed **three equal-count tiers** of patristic distance, ranked by distance and deterministically tie-broken by *source species code alone* (without response data), with tier sizes and RMSE improvement for branch-three vs null in each tier. If distance ties cross boundaries, report tie count, not causal threshold.
   - **within-original-fold demeaned correlation** between patristic gap and branch-three squared error / loss improvement. Exclude singleton folds only from covariance (but not from overall RMSE/tiers), and report how many species/blocks and how much within-fold distance variance is available. This controls merely a fold-constant shift, *not* species phylogenetic dependence or measurement error.
   - signed block-level mean gain and its association with mean block gap (descriptive only). Do not use one tip as one independently evolved event.
8. **Exploratory directional diagnostics fixed before v30 outcomes**:
   - For **31 focal species, cap=3 and cap=5**, record whether **far-tier branch RMSE exceeds near-tier branch RMSE** for both chemical endpoints. Positive correlation/consistent tier differences are descriptive support, not confirmatory evidence. Explicitly preserve inversions and mixed cases.
   - For all 9 panel/cap combinations, report whether distance/gene-novelty rank associations have the same sign for anthocyanidins and flavonols and whether the effect is entirely fold-mediated (within-block association weak or sign reversed).
   - **No fit/threshold selected based on error patterns.** Any of these can fail without changing v26/v29 true results. No global pooled coefficient or meta-analytic effect-size across repeated caps.

## Negative controls, source and scientific gates

- Fail closed on wrong 2023 60-tip original tree or stale SHA, missing child branch length, zero/invalid all-pair distance, wrong 31/59 source species denominators, distance matrix non-symmetry or nonzero diagonal, missing training coverage, unpaired/test-leaked prediction, mismatched source model RMSE, mislabeled branch genes, or v26–v29 source receipt drift.
- A synthetic toy tree with known path lengths must reproduce exact shortest paths; a zero-length edge is allowed but cannot make all measured distances zero. A missing branch length is **not** silently a zero. Mutating test HPLC outcomes after fold predictions must NOT change training-fitted predictions, distance tier membership, or gene-expression novelty.
- No correlation p-values, bootstrap pseudo-replicated species, CI from overlapping folds, or universal "evolutionary predictability law".
- Source record is an original phylogeny and processed chemical observation, not dated ancestral colour-loss events or taxon-level independent evolutionary experiments. Neither pigment synthesis gene causal function, functional enzyme assays, mating system, pollinator interaction, offspring persistence nor genetic trait re-expression is measured.
- If the 2023 published Newick lacks the required **per-edge branch-length information**, stop and record **SOURCE_BRANCH_LENGTH_NOT_AVAILABLE** rather than claiming numeric patristic results. Do not fall back to topology-hop distances without a separately frozen new protocol.

## Deliverables

Versioned standalone calculation, adverse tests, and dedicated source-authenticated GitHub Actions, plus aggregate-only derived readout and scientific interpretation. No original processed RNA/HPLC data or original Newick published in ANTO. Draft stacked on v29 PR #39 and earlier PRs; do not merge before review.

**At this protocol commit no v30 distance/error outcomes had been evaluated.**
