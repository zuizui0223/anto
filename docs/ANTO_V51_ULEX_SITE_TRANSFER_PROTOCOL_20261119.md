# ANTO v51 — Original Ulex six-site elevation-analogue transfer (posthoc)

2026-11-19. Study design fixed before v51 results were read, but research question is retrospective AFTER ANTO v46–v50 and SDEA v15. This is one independent original *Fabaceae* author study relative to the 2024 *Phlox* source; it is the SAME original Ulex archive already used in SDEA v15 and therefore adds no independently collected Ulex samples.

## Crucial source grain
Castellanos et al. 2023, Proc. R. Soc. B DOI 10.1098/rspb.2023.0141; archived data DOI 10.5061/dryad.9zw3r22jp; canonical Zenodo record 7761289. Authors explicitly state the six sampling sites belong to **one genetically continuous regional Ulex population**, NOT six differentiated genetic populations. Do not count six sites as evolutionary replicates. The original 240 tagged plants include 212 eligible fruit/flower trait plants, of which 202 overlap with source 225-genotype SNP relatedness labels. Source six site counts on the SAME 202: AresM 32, Cheste 37, Chiva 33, Montserrat 37, Simat 37, Sot.Chera 26. Original pollinator census (569 windows) is patch-level with NO individual key, so it is quarantined from original same-plant reproductive outcome data. Flower area and dry flower mass are morphological traits, not color regain or fitness-causal pollinator experiments.

## Fixed original source admission
Fetch ORIGINAL canonical Zenodo source metadata record 7761289 and require all five files and MD5 values registered in data/external/anto_v51_ulex_original_site_registry_20261119.json (fruits.and.mean.floral.traits.csv, floral.traits.csv, relatedness.matrix.csv, pollinator.census.csv, README.md). Never redistribute original per-plant source rows, SNP matrix or pollinator census. Parse 240 source plants, skip only explicitly missing/author #VALUE! trait cells, maintain original fruit scar counts identity scars=fruits+nofruit and 212 eligible phenotype-fruit rows. Relatedness matrix used to identify only the 202 source-genotyped matched plants, not to fit a genetic fitness effect. Require source original six sites, numeric provenance, full 202 identity and original SDEA v15 source-derived M0/M1 heldout MSE reference.

## Fixed whole-site forecast and distance analogues
Original response: individual twig fruit fraction fruits/scars in [0,1]. Original fixed predictor models from SDEA v15:
M0: intercept + training-standardized elevation with ridge alpha=5 on coefficient (unpenalized intercept).
M1: M0 plus training-standardized area, weight and the squares of their standardized values, all ridge alpha=5 except intercept.
Predictions clipped to [0,1]. All normalization is training-only; target fruit outcome never influences fit.

For each of the six **original sampling-site** holdouts (no within-site train/test mixing):
- FULL5: fit on ALL five remaining sites; anchor exact original SDEA v15 M0/M1 results.
- NEAR2: fit on two of the five training sites whose SOURCE SITE MEAN elevation is closest to heldout site's SOURCE SITE MEAN elevation, based only on environmental covariates; tie-break alphabetical site ID. Target elevation covariates are allowed for site-analogue selection; this is transductive environmental matching, not unknown-site prospective prediction.
- FAR2: fit on two most distant original source site mean elevations.
- ALL2: evaluate each of 10 possible unordered training-site pairs from five remaining. EXACTLY six times ten = 60 pair experiments for each M0/M1, with NEAR2 and FAR2 included rather than extra independence.

Site-mean elevations computed from original 202 eligible plant-level measured elevation; no response-dependent pair choice. Each site heldout model predicts only its eligible original plants once within its own regime. FULL5 and NEAR2/FAR2 models use identical features and ridge/response, differing in site-support set only (also changes the number/distribution of training plants and characteristics). No hyperparameter search, site reclassification, clustering or phylogenetic effect inference.

## Frozen outcomes, limits and falsification
Retain original six-site FULL5 M0/M1 MSE values and equal-site mean MSE gain, exactly 6 NEAR2 and FAR2 values and all 60 pair MSEs. For each compare M0→M1 trait contribution per site; also compare M1 total predictive risk of FULL5 vs NEAR2 vs FAR2 while explicitly reporting M0 baseline shifts. Report source-site-equal and original-individual-weighted metrics and number of improving sites, including negative. For each heldout site, report NEAR2 and FAR2 rank of M1 test MSE among all ten pairs, no p-values or confidence intervals. No favorable subset cherry-picking.

A positive NEAR2 transfer does NOT prove that geographic or pollinator environment causally gates floral selection. A negative NEAR2 transfer falsifies the naive “nearest elevation always improves transfer” within this one source. The six reference sites are not distinct genetic populations, and Phlox log1p fruit COUNT and Ulex fruit FRACTION are noncommensurable endpoints: do NOT pool numerical effect sizes, or call this general adaptive trait regain.

## Release/repro
One standard-library Python source-MD5 audit and training-only ridge, unit tests for source identity/cohort original source site grain, six FULL5 parent benchmark parity, 60 ALL2 and sorted-near/far ranks, no pollinator-plant join, no genotype reuse as independent observations, source-derived complete JSON receipt. GitHub Actions reruns original MD5-gated source and compares JSON after frozen. Keep Draft stacked on ANTO v50 PR #60. No main merge or general biological claim before source CI and independent review.
