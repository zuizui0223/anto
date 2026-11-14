# ANTO v29 — fine, discordance-bearing true-tree folds versus fold-size matched label controls

**2026-11-14; exploratory follow-up fixed AFTER inspecting v28's identical coarse-fold limitation, BEFORE evaluating any new fine-fold prediction.**

This is explicitly **not an independent preregistered confirmation** of v26–v28. Reuse precisely the 2023 original Wheeler et al. species-mean processed flower RNA/HPLC and original 2023 dated tree, the published 2024 Pezzi et al. ASTRAL tree / IQ-TREE / sample map, exact original-source bytes/pinned Git blob identities, and the immutable v26/v27/v28 original readout receipts.

## Scientific gap

ANTO v28's predeclared 31-focal-species predictions improved under multiple published tree definitions, but **the three largest heldout blocks were IDENTICAL** across all three trees at the 15-tip maximum, and 2023 and 2024 ASTRAL were also identical at the 8-tip maximum. Repeating predictions on the same folds does not test the predictive consequences of phylogenetic disagreement.

## Fixed exploratory tests

1. Same **31 real 2023 HPLC/RNA species** in the same published 2024 focal species overlap; source authors' original OUT anchor **HUTE**; same three topology sources: Wheeler 2023 dated 11-gene source, Pezzi 2024 ASTRAL species-assigned tree and the exact one-individual-per-species selection already frozen in v28 for Pezzi 2024 IQ-TREE. Same original models and lambda=10, both current measured chemical totals. NO new model selection, response-specific clade selection, updated tuning, or new HPLC cohort.
2. Five new fine clade maxima **2, 3, 4, 5, 6 original tips** (all five retained, however they perform). A clade is a complete rooted subtree after HUTE-outgroup rooting and pruning to the common 31 species, chosen solely by topology. Any unavoidably singleton folds remain in the denominator. Require each species predicted exactly once per fixed tree/cut/endpoint/model. Count and publish positive and negative folds.
3. Before looking at any new fine-cutout result, compute **fold-partition identity**, numbers of distinct partitions and pairwise species-pair co-holdout disagreements (among 31 choose 2 = 465 pairs) across the three trees. A purportedly distinct tree that generates the *same* heldout partitions supplies **no independent CV comparison**.
4. A **destructive topology-label negative reference** under the original 2023 tree only, at maximum-fold sizes **3** and **5**: generate **99 predetermined seeded uniform species-label permutations** across the original tree's fixed fold slots, preserving the original fold-size multiset, response and expression values and every original measured species, while breaking their phylogenetic co-assignment. No real clade or outcome is used to accept/reject permutations. Use deterministic Python random.Random seed **20261114 + cutoff**, generate 99 draws in order, do not pick a favorable draw. The random fold comparison is NOT a statistically valid species-exchangeable permutation test; phylogenetic dependence, selection of module in the inspected original sample and nonindependent folds prevent a confirmatory p-value.
5. For every fixed model/end-point/fold report original **species-weighted RMSE improvement (%) versus the training mean in exactly the same fold**, and block count / multi-species coverage. For the size-matched label reference report all 99 RMSE improvements per non-null module with median and 2.5/97.5 empirical interpolation quantiles, and fraction equal or higher than the original v29 tree score (descriptive rank only). Keep all results even if reference random partitions predict as well as the original clades.
6. Directional diagnostic: if there are NO discordance-bearing partitions at caps 2–6, verdict **NO_EVALUABLE_PUBLISHED_FOLD_DISCORDANCE**. Otherwise test whether branch-three predictions improve on the train-only mean for both chemical endpoints in **every distinct real published fold partition across all five maxima**; if any fails, verdict **MIXED_DIRECTION_WHEN_CLADES_DISCORD**, regardless of earlier positive v27/v28 results. This is an exploratory falsification diagnostic, not a new formal power analysis or external plant study.

## Main interpretation/limits

If original and arbitrary fixed-size partitions perform similarly, the original signal is consistent with a **general within-panel gene–pigment association**, not a specific predictive benefit unique to phylogenetic groupings. If some fine discordance-bearing clade partitions fail, retain those negatives: the original v27/v28 coarse-fold results remain true at their reported granularity but are NOT universally transportable across alternative clades.

The three tree inputs are **not** independent pigment experiments or posterior phylogeny samples; the 31 common species are a selected overlap of published RNA and HPLC sampling, not a global or cross-tribe reference distribution. Repeated splits and label permutations cannot be pooled as 99 independent natural evolutionary histories or combined into a valid meta-analysis. No visible floral colour, ancestral genetic toolkit function, intergenerational assay, pollinator behavior or descendant fitness is measured.

## Reproducibility/provenance

The v29 script must FIRST reproduce the entire original v26/v27/v28 source-authenticated, checksum-locked result, including published 2024 original trees, then evaluate the exact declared fine cuts and 99 synthetic relabelings. Original-source GitHub Actions CI emits an exact derived aggregate ledger, real source trees and original raw study data remain uncommitted. Review the new result after the source-run succeeds and do not change the predefined cuts or seed to maintain a desired biological conclusion.

**Status at protocol commit: no v29 fine-cut pigment prediction or random grouping result viewed.**
