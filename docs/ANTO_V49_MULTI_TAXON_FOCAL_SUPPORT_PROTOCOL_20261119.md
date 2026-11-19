# ANTO v49 — source-focal training-population restriction with v46 multispecies model held fixed

**2026-11-19; protocol frozen before examining v49 outputs.** This is a *post-v46/v47/v48 hypothesis development* analysis, not independent confirmatory science. ANTO v46's pre-outcome all-taxon fruit transfer criterion was NOT_ESTABLISHED; retain that failure without changing any earlier receipt.

## Precise motivation

v47 tested only evaluation-population count, and showed that the full v46 deamii source-average negative fruit MSE gains in shared (G) and garden (E) slopes were below 140/140 non-deamii 3-population reaggregations. However, v46 had **2 deamii same-taxon training populations** for each deamii holdout versus 7 amoena and 8 pilosa for theirs, while sharing all remaining taxa. v48 matched **2 same-taxon train + 1 heldout source population**, but also *removed other taxa* entirely and changed v46 nine garden×taxon intercepts into three garden-only intercepts. In v48, 34/140 (G) and 44/140 (E) non-deamii panels scored as badly or worse than deamii: this cannot identify the population-support effect separately from changed model domain.

## New, more controlled scope — original 2024 source, original v46 predictive architecture

Source remains Goulet-Scott et al. 2024 original **PhloxHopkins/PhloxFieldAdaptation**, commit d792fe58cce36721394adfd97e9f412396578e86, file Dataset_fitness_leaftraits_distance_lat_long.csv, Git blob SHA-1 4c780b0abe644ac69dae86a0450eb598c160626b, exactly 175,377 source bytes. Reuse v46's authenticated source parser, garden × taxon intercepts, 2 leaf PC axes, training-only standardization, fixed ridge lambda 10, and its *actual unchanged Python fit/prediction functions*. Do not change original v46 or 2014 PhloxSelection model. Original: 963 ramet rows, 284 cloned genotypes, 20 source populations, three common gardens, three taxa (amoena 8, pilosa 9, deamii 3). No original raw observations redistributed, new field samples or independent biological replicates.

Enumerate exactly **C(8,3)+C(9,3)+C(3,3)=56+84+1=141** original *same-taxon focal* triplet panels. For every triple, **each** one-source-population holdout (423 complete folds):
- Include **all original source populations of every NONFOCAL taxon** in the training set, unchanged for the focal taxon.
- Include exactly **two other source populations of the FOCAL taxon**, with the third completely withheld as the target; do not train on other focal populations.
- Retain all three original gardens for training and heldout evaluation; all original cloned genotype ramets remain within their source-population split.
- Fit the original v46 B0/G/T/E/TE models with **nine garden×taxon intercept cells**, fixed original source model ridge lambda 10, training-only centering and standardization. Nothing is recomputed from heldout response; model selection is not optimized using targets.
- The deamii three-member panel has the other 17 nonfocal source populations always included, so each deamii fold training set equals the original v46 full 19-population training set **exactly**. Hard numeric gate: deamii test MSEs for B0, G, T, E, TE and both fruit/flower endpoints MUST match all three corresponding original v46 frozen population MSEs within 1e-7 (parent stored eight decimals). Fail if this anchor breaks.
- Non-deamii panels retain all other original species in train, but drop their other 5/6 focal source populations. This changes **focal training support only** relative to v46 within each focal taxon, not the model family or availability of nonfocal original taxa. However, original row counts, source genetic diversity, focal trait distributions, and relative between-taxon sample weights all change together: isolate no unique causal effect of population count.

## Frozen outcome and complete reporting

Primary: original log1p fruit counts, secondary: original log1p *flower number* (NOT floral colour). v46 biomass log1p remains blocked by 13 negative biomass observations, minimum −233, and must never be reranked or imputed. Each triple panel has three heldout population MSEs per B0/G/T/E/TE for fruit and flowers, plus training/heldout ramet counts and all 2 focal + all nonfocal original training population keys. Full 141×3 fold-by-model derived JSON, no raw CSV source.

For each triple and model m, compute population-equal relative heldout MSE gain = 1 − sum(m fold MSE)/sum(B0 fold MSE). Compare deamii's **one** 3-source-pop panel with 56 amoena panels plus 84 pilosa panels, record #negative and #at-or-below-deamii for all four original model variants, both endpoints, plus taxon-specific min/median/max and pooled reference median. No user-selected favorable subset. All 140 source-triplet panels *overlap*, hence not a valid randomization p-value, independent study replication, calibrated predictive interval or species effect test.

A further numeric transparency checkpoint records original v46 gains, original v47 140-subset reference counts, and v48 two-train-only benchmark as *existing* historical baselines; does not rewrite them or pool quantities as independent effects.

**Decision framing prereadout:** We do NOT define a new hypothesis-test pass/fail, because this posthoc sequence has undergone several related model/sampling comparisons. If deamii stays below all/most 140 references when the nine-cell, multi-taxon predictor is preserved, report a *source/fitting-domain sensitivity*: v48's domain exclusion likely materially affects ranking. If many 140 panels are as poor as deamii, v47's exceptional-tail interpretation is NOT robust even when nonfocal taxa and model family are retained. Neither finding establishes a biological taxon interaction or shows that sample size alone is causal. Preserve negative outcomes and model dependence.

## Reproducibility, limits, review

Commit protocol and source registry BEFORE any v49 outcome. Build standard-library original-source script importing v46 parsing/model methods. Include fail-closed input source blob, population taxon classification, genotype leakage, two-focal + all-nonfocal train rule, nine full training cells, three original gardens, and original v46 deamii fold MSE exact anchor. Retain every original sampling interval and method caveat. Complete source-derived frozen JSON, adverse tests, CI that downloads pinned original author CSV and exactly regenerates full receipt. Draft PR stacked on v48 #58, not main; no new independent external biological validation until new source populations/gardens/flower–visitor–seed data are collected.
