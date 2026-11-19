# ANTO v47 — post-outcome source-population-count matching diagnostic

**2026-11-19 | EXPLORATORY / POSTHOC**. This is NOT a preregistered biological hypothesis test. v46 already revealed a global positive leaf-trait → fruit prediction improvement and negative within-deamii source-population mean. The scientific question now is narrow: *is the exceptional negative deamii three-population score ordinary among arbitrary groups of three source populations from either other original taxon?*

## Immutable evidence provenance and admissible units

Read the exact committed ANTO v46 source-derived aggregate (NOT raw original 2024 observations):
- repo: zuizui0223/anto, branch analysis/anto-v46-external-phlox-population-transfer-20261118
- derived path: results/anto_v46_phlox2024_population_heldout_20261118.json
- pinned Git blob SHA-1: f03b0eb0177f807317ee1f8a7c0e65ee408d903b
- original author's source pinned commit d792fe58cce36721394adfd97e9f412396578e86, CSV Git blob 4c780b0abe644ac69dae86a0450eb598c160626b
- 963 original clonal ramet records, 284 original genotypes, 20 source populations from three original Phlox taxa in three common gardens; no new original studies or generations.

Original v46 20 holdout folds all exclude each source population and all of its clones across three gardens. This v47 reaggregates **already held-out fold scores** only; it never refits slopes, resplits ramets, uses heldout targets to train, invents synthetic new field observations, or replaces original baseline. The v46's full prereadout **fruit** transfer criterion FAILS for all four leaf-PC model variants. Preserve that failure exactly.

## Posthoc comparison, not a null test

For endpoint e and model m in G/T/E/TE, and any subset U of three original SOURCE POPULATIONS, define:

    relative_gain(U,e,m) = 1 - sum_{p in U} MSE(p,e,m) / sum_{p in U} MSE(p,e,B0).

This weights each distinct source population once irrespective of its different number of original ramet observations. Exactly one deamii triplet exists (all its 3 original populations), compared to **all** C(8,3)=56 amoena-only and C(9,3)=84 pilosa-only triplets. The 140 reference triplets are overlapping combinatorial subsets of original v46 scores, not 140 independently collected studies, random samples, or refit CV. Also report an expressly SECONDARY **pooled** C(17,3)=680 mixed- or same-taxon triplet reference. No pooled reference is an exchangeable sampling null without additional assumptions; **no p-values, confidence intervals, null rejections, posterior odds or general causal species effects**.

Per model/end point, freeze output for: deamii three-population gain; population-level sign and number of positive folds; each deamii leave-one-population *arithmetic* two-population residual gain (NOT refit); reference counts n with gain <= deamii; counts with negative gain; min, median, max; and checks that parent 20-pop source summary and source-denominator integrity remain unchanged. For the reference, use the 8-decimal parent per-population MSEs directly and do not pretend 8-digit rounding is full original high-precision per-ramet output. A tie within 1e-10 for reference comparison is treated as equal.

The primary comparison domain is original fruit count (log1p); original flower count (log1p) is an already-observed SECONDARY diagnostic. Source original biomass has 13 negatives (min -233), hence v46's preregistered log1p biomass remains a technical NONADMISSION. No new scoring.

## Interpretation gate

If the deamii gain is below all 140 same-taxon size-three references, say ONLY: *within this original 2024 study's fixed v46 20-population predictions, equal evaluation-population count alone does not reproduce the deamii aggregate for that model*. This is descriptive evidence against the idea that merely selecting 3 out of the other taxa commonly produces an equally poor score. The analysis cannot exclude chance across **independent original populations**, explain the causal basis (genetic, pollinator, site, leaf morphology, sample variance) or establish a universal taxon-specific slope, because:
- deamii has only three biologically original source populations;
- both reference lists reuse original populations and overlap substantially;
- all original 20-fold models still train against the same three gardens and, for most designs, many other taxa;
- the reference matches **evaluation-source count** only, not training sample size, population baseline risk, habitat, trait range, evolutionary distance, or phylogenetic independence.

Do not choose the most favorable model/endpoint post hoc. Four original model families and both original endpoints must be retained. No need to claim ANTO × EGWE × izu-core mechanism convergence from this single leaf-fitness cohort.

## Required reproducibility and review

Implement a hash-locked Python 3 standard-library reaggregator, exact deterministic receipt with all model/endpoint reference outcomes, adversarial numeric/source/denominator/grain tests, read-only GitHub Actions, and an explicit result report. v47 is a **Draft stacked on v46 PR #56**, no main merge without external scientific review. The next stronger biological data step is independent deamii source populations and genuinely new gardens/experimental treatment combinations, not additional resampling of this same CSV.
