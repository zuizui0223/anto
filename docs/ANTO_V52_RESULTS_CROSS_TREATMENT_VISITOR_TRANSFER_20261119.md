# ANTO v52 — Sinapis original-source rainfall-treatment out-of-domain visitor-information transfer

2026-11-19. This is a retrospective, **posthoc scientific question** motivated by prior SDEA v17 shelter holdouts and ANTO v46–v51 prediction-domain failures. The v52 protocol and fixed model/placebo readiness contract were committed **before its original-data readout**. This is NOT independently collected ecological evidence relative to SDEA v17, a new trait-evolution experiment, a causal test of visitation or a general macroevolutionary transfer theorem.

## Exact original study provenance and methodological distinction

Höfer et al. 2022 Functional Ecology DOI 10.1111/1365-2435.14233, original author Dryad DOI 10.5061/dryad.66t1g1k5f and Zenodo record 7335086. The original study experimentally varied water availability (Control, Drought, PrecipRed) under five rain shelters in ONE 2020 field trial. It originally reported that lower water supply reduces plant seed production predominantly through reduced flower abundance, rather than a large induced shift in flower visitation.

The source-locked SDEA v17 original parser was imported at pinned user-owned Git commit a8990f5255b1b27ed5fec85587bf09ef6390882a, Git blob SHA1 5f797c1d891c8f8de711909a56ae26530841a2ba. Original SDEA v17 derived study receipt Git blob SHA1 fa2419d32b74f22f67ea9cbfd79123109db4fdf0. **ALL FOUR canonical Zenodo author files passed original MD5**, with 53 source plant list rows, 53 floral trait rows, 427 repeated visitor-group records and exactly 52 joined complete individual plant samples. The original five shelter counts remain (11, 9, 8, 12, 12). Original source *observetime* is counted once per plant, not 427 times; original visitor NrInteraction summed over original groups. Each same-plant motherID is a maternal family identifier, NOT a sequenced genotype.

**Crucial reproducibility gate:** Exact original SDEA v17 source code replayed all 5 old original SHELTER-heldout folds × 4 original SDEA models with **maximum absolute MSE difference exactly 0** from previously frozen source. These old results are not new independent source evidence.

## New 3-rainfall-TREATMENT heldout prediction target and model

In v52 each rainfall treatment is withheld across ALL five shelters; train on the other two treatment conditions. Original complete-case n: Control 16 (train 36), Drought 17 (train 35), PrecipRed 19 (train 33); every original eligible plant appears in exactly one test fold.

All 4 v52 models omit treatment indicators on purpose, since an entire treatment level absent from training cannot be assigned an estimable treatment-dummy coefficient. Thus the models are **treatment-blind** and not directly interchangeable with original SDEA v17 treatment-adjusted shelter-holdout regressions. Predict original source log1p estimated same-season seedsperplant with ridge lambda 4 on slopes and training-only covariate normalization:
- U0: log1p(total flower supply) only.
- U1: U0 plus same-plant original mean floral display size and petal width.
- U2: U1 plus same-plant original log1p(number of recorded visitor interactions divided by its own observation effort once).
- U3: U0 plus visitor rate directly, with no additional shape predictors.

All predictor observations come from the same season as seed estimates: this tests the transport of an existing **conditional association** into an unseen experimental treatment, NOT a deployable predrought prediction or causal mediation.

## Original source-derived results: independent holdout axis (all measured plants retained)

**Primary comparison U1 minus U2 MSE:**

| Entire original rainfall treatment withheld | original n plants | U1, flower supply + shape MSE | U2, plus direct observed visit intensity MSE | U1−U2 MSE gain | relative MSE reduction |
|---|---:|---:|---:|---:|---:|
| Control | 16 | 0.414627221 | 0.342116004 | **+0.072511217** | **+17.49%** |
| Drought | 17 | 0.191035055 | 0.087503948 | **+0.103531108** | **+54.19%** |
| PrecipRed | 19 | 0.230948716 | 0.171285362 | **+0.059663354** | **+25.83%** |

- **Positive in all 3 of 3 withheld rainfall treatments**; equal-treatment MSE gain **+0.07856855933109674**, original-plant-weighted gain **+0.07795792338930461**.
- Mean equal-treatment U0 MSE **0.231895982**, U1 **0.278870331**, U2 **0.200301771**, U3 **0.174178772**.
- Floral morphology added to flower supply **worsened U0→U1 in all 3** treatments: equal-treatment absolute gain **−0.04697434866912429**. Original U0→U3 direct-visitor increment was positive in **all 3** treatments (+0.057717210150544006 equal-treatment).
- U3, the parsimonious supply+visits alternative, has the lowest **whole-study equal-treatment mean MSE** (0.174178772); in the Drought fold U2 (0.087504) is better than U3 (0.095945). Do not omit this directional heterogeneity or reselect U3 to rebrand the prereadout U2 primary comparison.

### Original source plant-visit label placebo, no randomized causal p-value

Fixed 200 source stratum-preserving shuffles of the *observed plant visit-rate predictor* within each original shelter×treatment combination, RNG seed 20261119, preserving each plant's other traits, treatment assignment and same-season seed estimate, then refitting the same U2. Observed equal-treatment U1−U2 gain +0.07856856; placebo distribution median **−0.00263180**, q97.5 **+0.02203237**, maximum **+0.03519084**; **0/200** shuffled label configurations met or exceeded the original gain.

The frozen three-part **descriptive** gate (positive equal/weighted mean, >=2/3 positive treatment folds, observed > placebo q975) **PASSES**: DESCRIPTIVE_WITHIN_SINGLE_TRIAL_CROSS_TREATMENT_VISITOR_INFORMATION. It is NOT a formal p-value or an experimental manipulation of insect visitation: original rainfall was randomized, but the visit predictor and reproduction were observational and contemporaneous. Maternal background, flowering phenology, local source shelter differences and unobserved insect quality remain possible confounds.

## What this adds beyond the original SDEA v17 and ANTO Phlox/Ulex results

- Prior SDEA v17 *held out whole rain shelters*: visit increment after morphology was positive for **4/5 shelters**, equal-shelter +0.063121875; those models included rainfall treatment indicator variables and were evaluated against new shelters.
- **New ANTO v52 target** holds out an entire original *rainfall treatment* across existing five shelters. Visitor intensity retains descriptive predictive information in **3/3** treatment groups after the fixed U0/U1 morphology controls, even though the original morphology-only layer is harmful under this model.
- ANTO v51 Ulex 2023 (one continuous population, six sampling sites) was a separate Fabaceae source study without same-individual pollinator visits; ANTO v46–v50 Phlox 2024 original leaf-PC→fruit count was Polemoniaceae without direct visitor records. Their heterogeneous responses and measurement scales **cannot** be treated as three identical independent cross-species effect estimates or quantitatively meta-analyzed.

**Bounded ecological interpretation:** The plant's observed interaction-rate information can transfer in this one Sinapis experiment between drought-manipulation treatments, but no general proof exists that species-specific pollination services or behavioural learning mechanisms explain ecological rescue. In Sinapis the actual year-2020 trial original-author estimated seed number, flowering supply and visitor encounters were measured in the same season; visitor count is NOT pollen deposition or offspring viability, and future trait expression/ancestral flower colour has not been measured. The original drought causal biology remains as the authors reported it.

## Reproduction and merge gate

- Frozen posthoc protocol before outcome: docs/ANTO_V52_SINAPIS_TREATMENT_HOLDOUT_PROTOCOL_20261119.md
- Exact original four source MD5 + SDEA v17 pinned code/source/receipt: data/external/anto_v52_sinapis_treatment_original_source_registry_20261119.json
- Frozen derived SDEA v17 prior shelter-fold MSE and source cohort: data/external/anto_v52_sdea_v17_derived_original_anchor_20261119.json
- Original source parsing + treatment-heldout ridge + 200 plant-identity placebos: scripts/analyze_original_Sinapis_treatment_holdout_v52.py
- Adverse original source, treatment grain, leakage and noncausal placebo tests: tests/test_original_Sinapis_treatment_holdout_v52.py
- Exact derived-only source outcome: results/anto_v52_sinapis_treatment_heldout_original_source_20261119.json
- Read-only original source Actions: .github/workflows/anto-v52-sinapis-rain-treatment.yml

Initial v52 CI job failure was a **coding manifest-field typo** in accessing the prior SDEA receipt (original_source_zenodo_record vs source_zenodo_record), corrected without modifying source, model, training design, permutations, or outcome. After this repair the original-author four-MD5 and prior SDEA 20 MSE replay CI run **#44533054601 SUCCESS**. Exact current HEAD re-download/reproduction and scientific review remain required before any main merge. Keep v52 **Draft**, stacked on ANTO v51 PR #61. No prospective or independent biological reviewer approval claimed.
