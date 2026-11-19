# ANTO v51 — Original Ulex six sampling sites: environmental-analogue transfer fails despite positive relative incremental floral information

2026-11-19. The original-author source MD5 and source-grain controls were applied to Castellanos et al. 2023 Proceedings B original Zenodo record 7761289. ANTO v51 is a NEW predictive contrast compared with the 2024 Phlox source, but is NOT a new Ulex original study relative to already-authenticated SDEA v15. All v51 methods were frozen before the first readout; the scientific question was retrospective after prior results.

## Mandatory biological-grain correction (original author)

The 2023 original authors explicitly wrote in their field-sampling section that **the six Ulex sites fall within one genetically continuous regional population** with little differentiation, so six sites MUST NOT be called six genetically distinct evolutionary populations. Original source: 240 tagged individual plants, 212 usable phenotype–fruit records, original relatedness matrix 225 genotypes, **202 SAME individual genotype×flower trait×fruit rows**, six distinct sampling locations (AresM32, Cheste37, Chiva33, Montserrat37, Simat37, Sot.Chera26). Original 569 patch census pollinator observations do not share individual plant IDs, so no flower-visit-to-individual-genotype/fruit link was invented. Source original five files MD5 verified; no raw author rows redistributed.

The original paper already reported stabilizing selection on standard petal area and relatively high heritability for flower dry mass but not current selection on that trait. ANTO v51 does **not** rediscover these original author results or measure new causal pollination.

## Source-frozen original model replication

Model M0 = original SDEA v15 ridge(5) elevation-only; M1 = M0 plus original standard-petal area and flower dry weight and their standardized quadratic terms, training-only variable scaling, predictions clipped to [0,1]. Original current-season outcome = fruits divided by twig flower scars. ALL six heldout-site M0/M1 scores reproduce the independent SDEA v15 source-derived receipt within **max absolute MSE error 1.214306433183765e−17**. The source is the SAME original Ulex 2023 biological experiment already in SDEA, not a second species effect in any meta-analysis.

## New site environmental analogue comparison

For each original site, compare ALL five other source training sites (FULL5) versus two training sites selected by closest original site-mean elevation (NEAR2) versus two with farthest site-mean elevation (FAR2). All possible two-site training pairs, 10 per heldout site = 60 total paired-model exercises, are recorded with original source-site and heldout individual counts. Original heldout site mean elevations are available as environmental metadata, not heldout fruit outcomes: this is transductive, not prospective blind-site prediction.

| Source train support | Equal-site M0 MSE | Equal-site M1 MSE | Relative floral trait increment 1 - M1/M0 | Site M1 improves vs M0 |
| --- | ---: | ---: | ---: | ---: |
| FULL5 | 0.008734302507 | **0.007713521166** | +11.69% | 5/6 |
| NEAR2 | 0.040194630939 | 0.031636807505 | +21.29% | 5/6 |
| FAR2 | 0.449083092642 | 0.449083092642 | 0% | 0/6 |

**New crucial contrast:** selecting the two supposedly most analogous source sites raises the actual M1 prediction MSE by a factor of **4.1015** relative to full five-site training, despite the near2 trait *increment* being a larger fraction of its own substantially worsened M0 baseline. The near2 M1 total prediction risk is better than FULL5 in only **1/6 original heldout sites**. FAR2 is dramatically worse (≈58.22× FULL5 M1 risk), and equal FAR2 M0/M1 scores must NOT be misinterpreted as proof of zero flower biology: training outside the target environmental range and prediction clipping can cause degenerate loss equality.

The nearest-elevation chosen pair's M1 risk rank among the original ten different two-site training choices (1=best, 10=worst), by heldout original site is:
- AresM: rank 3/10; near2 M1=0.003055 vs FULL5 0.003723 (only near2 improvement);
- Cheste: rank 5/10; near2 M1=0.011430 vs FULL5 0.010701;
- Chiva: rank 4/10; near2 M1=0.002520 vs FULL5 0.002109;
- Montserrat: rank 5/10; near2 M1=0.015810 vs FULL5 0.015325;
- Simat: rank **10/10**; near2 M1=0.141155 vs FULL5 0.009311;
- Sot.Chera: rank **9/10**; near2 M1=0.015851 vs FULL5 0.005112.

Thus “environmentally closest site is always best training source” **fails** on this original six-site dataset. The original site-mean elevation alone does not ensure matched flowering trait distributions, sample support or other environmental covariates; FULL5 retains more original training individuals. The positive FULL5 floral trait MSE increment is **not** evidence of heritable adaptation or a new pollinator-function result.

## Why relative trait gain alone misleads

The relative quantity 1 - MSE_M1/MSE_M0 is *within a training regime*. It can become more positive merely because that regime has a substantially worse ecological M0 predictor, even while the complete M1 absolute risk gets much worse. Compare absolute heldout loss and same biological unit before calling a model or mechanism superior.

## Strong limitations / ANTO cross-study interpretation

Ulex 2023 (one Fabaceae species, one continuous regional population, 6 field sites) and Phlox 2024 (three Polemoniaceae taxa, 20 origin populations, 3 common gardens) are distinct original source studies but their endpoints are **fruit fraction vs log1p fruit count**, and their heldout units are sampling **sites vs source populations**. MSE magnitudes, evidence units, mechanism and test population distributions cannot be numerically pooled or called independent natural evolutionary mechanism replications. Original Ulex pollinator census is not individually joinable; no direct genotype–flower–visitor–stigma–compatible pollen–next-generation lineage was identified.

An appropriately bounded generalizable research issue is the danger of interpreting **positive incremental information** or a convenient ecological distance axis as adequate support for **actual out-of-sample portability**. The next data-driven step must explore the environmental analogue failure using independent sampling sites or a different original field experimental system; an isolated field study cannot establish a universal macroecological boundary.

## Reproduction

- Frozen v51 source-first protocol: docs/ANTO_V51_ULEX_SITE_TRANSFER_PROTOCOL_20261119.md
- Exact original Zenodo 5 MD5 and source grain manifest: data/external/anto_v51_ulex_original_site_registry_20261119.json
- Immutable derived benchmark against SDEA v15 original: data/external/anto_v51_sdea_v15_derived_anchor_20261119.json
- Python original-source downloader and all original six-site and 60-pair comparisons: scripts/analyze_original_Ulex_elevation_analogues_v51.py
- Adversarial source/biological-grain tests: tests/test_original_Ulex_elevation_analogues_v51.py
- Full source-derived six-site original environmental pairing, train sizes and losses: results/anto_v51_ulex_six_site_elevation_analogue_20261119.json
- Read-only source-md5 rerun: .github/workflows/anto-v51-ulex-site-analogue.yml

**GitHub Actions source-first run #44527971383 completed SUCCESS** with exact five original source MD5, original n=202, source SDEA v15 parity and 60-pair original readout. Keep Draft / stacked on v50 PR #60. The full stored-result exact-head replay is a separate CI gate and must pass before merge; main untouched.
