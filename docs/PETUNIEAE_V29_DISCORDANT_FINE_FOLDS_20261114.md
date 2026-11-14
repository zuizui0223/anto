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


## Source-authenticated results (2026-11-14; AFTER the frozen exploratory design)

**Original completed source run:** https://github.com/zuizui0223/anto/actions/runs/43584570303 (**SUCCESS**); every original v26, v27 and v28 derived source receipt was recalculated and checked for exact equality before v29's new models. The 2023 original HPLC/RNA, 2023 published original tree and two genuine 2024 published phylogenomic methods were reused without replacing any source, outcome or molecular predictor.

### Real fine tree discordance and predictive consequences

The true published 2023 source tree and 2024 ASTRAL species-assigned tree still produce identical heldout clade membership at **every** new cap 2–6. The source-registered 2024 IQ-TREE first-sample-per-species reduction, however, generates a different real clade partition at each cap.

| Maximum heldout subtree | Distinct CV partitions across 3 real published source trees | Species pairs differently co-heldout by IQ-TREE vs the original, of 465 | Anthocyanidin branch3 RMSE gain, 2023 / IQTREE | Flavonol branch3 RMSE gain, 2023 / IQTREE |
|---|---:|---:|---:|---:|
| 2 species | 2 | 4 / 465 | +7.9257% / +8.1435% | +0.2197% / +0.1618% |
| 3 species | 2 | 8 / 465 | +8.4755% / +9.4416% | +0.7612% / +0.4374% |
| 4 species | 2 | 15 / 465 | +10.0303% / +10.0692% | +0.7614% / +1.0729% |
| 5 species | 2 | 29 / 465 | +10.2781% / +10.7941% | +0.8869% / +1.9901% |
| 6 species | 2 | 34 / 465 | +12.1232% / +11.0175% | +2.9231% / +1.9553% |

The predeclared **post-v28 exploratory directional criterion** remains positive across all 10 distinct source-tree partitions × both chemicals. Its exact verdict is **EXPLORATORY_ALL_DISTINCT_FINE_PARTITIONS_POSITIVE**. This is still **one original 2023 chemical dataset with overlapping training folds**, not 10 independent biological replications, and the actual discordance in co-heldout species pairs is limited (4–34/465), not an exhaustive posterior sample of plausible phylogenies.

**Critically, positive signs are not robust magnitudes.** Restricting holdouts to real clades of 2–5 species reduces the species-weighted **flavonol** branch3 RMSE improvement to **+0.22–0.89%** on the original source tree, while anthocyanidin improvement stays around **+7.93–10.28%**. The common 31-species panel with the coarse 8/15 cuts had roughly +9.82–10.40% flavonol gains. These within-source contrasts have interdependent test folds and must not be interpreted as causal changes in pigment pathways.

### Pre-fixed 99 same-fold-size taxon-label control partitions

| Endpoint | Fine cap | Authentic original tree branch3 RMSE gain | Median from 99 source-label-randomized foldings | 97.5th percentile | Random foldings ≥ authentic |
|---|---:|---:|---:|---:|---:|
| Anthocyanidin | 3 | +8.4755% | +8.5788% | +9.8684% | 54/99 |
| Anthocyanidin | 5 | +10.2781% | +8.6762% | +10.5789% | 7/99 |
| Flavonol | 3 | +0.7612% | +6.0162% | +8.5204% | 98/99 |
| Flavonol | 5 | +0.8869% | +6.1244% | +9.8539% | 97/99 |

All 99 assignments per cap retain exactly the original source's 31 taxa, gene predictors, true chemical endpoints, exact number and sizes of folds; **only group membership is shuffled**. These are **descriptive destructive-phylogeny controls**, NOT a phylogenetic exchangeability null, valid p-values or new external replicates. The original v26 was inspected before the new exploratory v29 cut sizes were set; no confirmatory multiple-testing error control can be inferred.

### Implications and negative model outcomes

1. **Both chemicals contain a branch-module association**, but phylogenetically heldout **flavonol amount** is less well extrapolated between fine original source clades than between random partitions. This may indicate clade-specific response baselines or other structured chemical/expression confounding; it does NOT establish that the three branch genes lack function or that phylogenetic ancestry is a cause of the error.
2. At very fine caps, the previously positive 6-gene module is not consistently better than the 3-gene module, especially for flavonols, and other genes can outperform branch-three for anthocyanidin. All five alternative modules including CHI, early-three, late-two and joint-six are retained in the immutable result; the branch-module choice was carried forward from v26 rather than chosen using v29 outcomes.
3. **Reframe paper language:** evidence of **within-source target- and holdout-scale-specific predictive sufficiency**, not a general law of molecular evolutionary potential, gene regulatory restoration or validated cross-phylogeny portability. It would be incorrect to call the identical 2023/ASTRAL test splits independent confirmations.
4. The next credible empirical strengthening is a truly *independent original HPLC/RNA* dataset and/or functionally validated multi-species biochemical assays, with a prospective frozen model and taxon-level lineage holdouts. Gene-tree alternative estimates alone, on the same 2023 chemical rows, cannot supply that replication.

The full source-locked, **negative-case-preserving** machine-readable ledger is results/petunieae_fine_discordance_label_reference_v29_20261114.json. Original source files are not recommitted to ANTO. No natural ancestral flower-colour re-evolution, genotype–pollinator–fitness chain or cross-taxonomic quantitative meta-analysis was measured or estimated.
