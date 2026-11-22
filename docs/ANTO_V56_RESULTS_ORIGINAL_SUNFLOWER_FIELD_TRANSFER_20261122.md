# ANTO v56 — Original sunflower field-level bee visitation predicts heldout seed set, with source failures retained

**2026-11-22.** This is one original Mota et al. (2022) sunflower crop field experiment (original Dryad DOI 10.5061/dryad.h18931zpg, Zenodo 6671199), re-analysed **retrospectively**, not a newly collected study and not a preregistered prospective ecology experiment. The v56 source-inclusion, models, physical-field grouping, 200 placebos and descriptive pass/fail rule were frozen in `docs/ANTO_V56_PREOUTCOME_SUNFLOWER_REPAIRED_FRACTION_PROTOCOL_20261122.md` **before v56 predictive output**. The reason to study the source and the three invalid rows were known from v53–v55 earlier audits.

## Source admission: exact original author records and the prior STOPs

All three original Mota author files match original registered MD5. Exact v53 source-key ancestry: 80 original matched field-years (2017: 41, 2018: 39) and **51 UNIQUE physical sunflower fields**, some measured in both years. Both survey and harvest author tables have the same original physical field name, region, field management category and field-year key, but **do not identify the same observed and harvested individual sunflower heads**. The original sources have 2,875 harvested source rows and 11,683 original one-minute visitor survey rows, of which 2,217 and 11,682 rows, respectively, belong to joined field-years.

- **v54 SOURCE STOP RETAINED:** 1,056/2,217 original source abortion-count cells were missing and 4,087/11,682 original survey rows had `No_Interactions_TOTAL != No_Interactions_APIS + No_Interactions_WILD`; the original v54 stricter three-count arithmetic criterion was not met. No v54 fitted effect existed.
- **v55 SOURCE STOP RETAINED:** directly observed `No_Seeds_Total` and `No_rep_units` were available for every original joined harvest row, but three 2017 source rows had **No_Seeds_Total greater than No_rep_units**. v55 strictly rejected these original ratios; its registered model did not run.
- **v56 EXPLICIT NEW DATA DEFINITION:** exclude exactly the three impossible seed-ratio original rows identified in the independent prior v55 source quality audit; **2,214/2,217 (99.865%)** harvest observations remain. No missing abortion counts are imputed, no seed fractions are clipped during source input, no alternative favorable source responses are substituted, and all 11,682 source bee surveys remain. The original 80 field-years/51 physical field groups and all original 3–4 aligned distance bands remain available after excluding three lines.
- Report original direct source TOTAL bee interactions per monitored sunflower unit as PRIMARY visitation. The disagreeing APIS+WILD reconstructed total is used only as a distinct, predeclared **measurement-definition sensitivity**, not a taxon-specific pollination effect or a forced count identity.

The registered original source quality ledger is frozen at Git blob `f9d6c01720ae0657977d09d0a92043b0eff0d7b5`; old v54 stop at `55eb9b942e7a0ab7f4e32139ac5305325fb74a83`; v53 original fieldyear identity ledger at `7a4fb27c7f26abfa8400d23b041e1a3d3e3424c9`. Old v55 fitting script `scripts/analyze_original_Mota_reported_total_bees_v55.py` was **reused unchanged** at Git blob `3f4b936101ef51e56e1fc661089d49b0efe8810b`. None of these earlier analyses are rewritten as successful.

## New original spatial heldout predictions — entire physical fields, including both years

Construct original y_field-year = equal mean across matched source distance bands of pooled full-seed / original total-reproductive-unit ratio. Original x_field-year = equal band mean of pooled recorded total bee interactions / monitored sunflower unit during one-minute census, then log1p. Each response remains a **field-year observed seed fraction**, NOT a same-individual bee–harvest effect.

All 51 original physical fields were held out one at a time, **both their years together**, using fixed ridge alpha=8 on slopes, unpenalized intercept, field context Year+Region+Field_Type, training-only visitor normalization and bounded predicted seed fraction [0,1]. No heldout seed outcome informed feature selection, transformation, intercept or training. There are 80 field-year predictions from 51 full source-field holdouts.

| Frozen original model | Equal-physical-field MSE (smaller is better) | Field-year-weighted MSE | Relative improvement over Z0, equal-field |
| --- | ---: | ---: | ---: |
| Z0: Region + Year + Field_Type only | **0.013559594** | 0.012809129 | reference |
| Z1: Z0 + source-reported TOTAL bee interaction rate | **0.011316545** | 0.010701690 | **+16.54%** |
| Z2: Z0 + separately reconstructed APIS+WILD rate, measurement sensitivity only | **0.011263140** | 0.010626575 | **+16.94%** |

**The original predeclared Z0→Z1 primary mean improvement is +0.0022430493** seed-fraction-squared MSE on equal-physical-field weighting, and **+0.0021074385** on fieldyear weighting, both positive. It improves **34/51** physical fields (66.7%) and worsens the remaining 17/51; all are preserved.

The Z0→Z1 observational fieldyear mean MSE gain is positive within **Burgos (+0.002544856)** and **Cuenca (+0.001711680)**, as well as within all three original habitat/field management categories **NonV (+0.003135554), SNV (+0.000851555), WFS (+0.002894468)**. No source raw predicted score required clipping to zero or one in any model (Z0/Z1/Z2 all 0/80 outside [0,1]).

The secondary Z2 separately reconstructed exposure gives source equal-physical-field gain +0.002296454, with **34/51** original fields improving and 17/51 worsening. This similar direction across two data definitions is a *measurement-description robustness check* only. Original total and reconstructed APIS+WILD are empirically discrepant in 4,087 survey rows; do not call these independent biological confirmations or estimate honeybee-versus-wild-bee mechanisms from this result.

## Original visitor identity placebo — NOT causal inference

Exactly 200 fixed RNG seed 20261119 shuffles of the paired (reported total, reconstructed Apis+wild) **field-year-level** observed visitor rates, permuted only within original Region × Year × Field_Type strata. Same original field-year seed outcomes, field identity, folds, context and training policy were maintained; Z1 was re-fitted each shuffle.

- Observed original source equal-field Z0→Z1 gain: **+0.0022430493**.
- Original label-shuffle distribution median **+0.000041679**; 97.5th percentile **+0.001031873**; largest sampled gain **+0.002080661**.
- **0/200** shuffled exposure assignments equaled or exceeded the observed gain.
- The original fixed descriptive readout gate (both weighting gains positive, > half original fields better, both regions and >=2 management classes positive, observed above 200-shuffle q97.5) PASSES, recorded as `DESCRIPTIVE_SINGLE_ORIGINAL_CROP_STUDY_FIELD_TRANSFER`.

This is NOT a valid randomization-derived effect p-value: original bee visits were NOT experimentally assigned, field-year shuffling is not independent pollination treatment randomization, physical fields can share landscapes, and the seed/visitor records are observational. Original Mota source author paper already considered plant habitat, bees and crop productivity: v56 is a *heldout-predictive* rather than new causal pollination effect question.

## Scientific comparison with Sinapis v52: no false same-plant replication

Sinapis (Höfer et al. 2022, Brassicaceae, 52 individual plants, five original rain shelters) reported observed original same-plant visits and contemporary seed estimates. ANTO v52 showed direct plant-level visitation added **positive heldout seed-prediction information across 3/3 entirely omitted rain treatments**, but original visitors were not randomized. Mota sunflower (Asteraceae) differs: original bee census and harvested seed-set data meet only at field-year level. **Both** separate original studies show positive predictive information from measured insect interactions in distinct holdout domains, but the measurement units, outcomes, experimental setting, biological mechanism and MSE scales are noncommensurable and must **NOT** be pooled as one causal interaction effect.

This gives a bounded cross-study research insight: **visitor encounter information can retain out-of-sample predictive value beyond simple contextual baselines in some distinct botanical settings**, even while the quality and scale of the predictor–outcome join can strongly restrict inference. It does not imply a universal pollinator effect, taxon-independent causal adaptation or evidence of ancestor petal-colour reappearance. Negative original source audits (v54/v55) and negative fields (17/51) are inseparable parts of the evidence.

## Complete source-reproducibility and merge gate

- Freeze-before-v56-output but post-v55 source STOP protocol: `docs/ANTO_V56_PREOUTCOME_SUNFLOWER_REPAIRED_FRACTION_PROTOCOL_20261122.md`
- Exact prior original MD5, author data grain, 3 invalid source rows, model blob and 51-field policy: `data/external/anto_v56_valid_seed_fraction_registry_20261122.json`
- Exact original Zenodo source retained external; source-ratio-validated 2,214-row aggregator + immutable v55 original models: `scripts/analyze_original_Mota_valid_seed_fraction_v56.py`
- Adversarial source and physical-field/whole-year leakage tests: `tests/test_original_Mota_valid_fraction_v56.py`
- Full 51-field/80 field-year original derived MSEs, original source counts, signed placebos, numeric clipping and verdict, **without original raw row values**: `results/anto_v56_sunflower_valid_fraction_field_transfer_20261122.json`
- Dedicated read-only source replay: `.github/workflows/anto-v56-original-valid-fraction.yml`

Initial MD5-gated original-source Github Actions #45012856622 completed **SUCCESS** with all model/quality tests, source 2214/2217 registration, 80 fieldyears and 51 full-field holdouts, and 200 observational placebos. A separate exact-current-HEAD replay after storing result JSON remains mandatory before treating this as reproducible for review.

**DRAFT stacked on PR #65**, main untouched. Independent biological peer review needed before merge and broader confirmation needs additional ORIGINAL studies with aligned same-plant timing, pollinator deposition and later seed-genetic/recruitment outcomes; this source cannot supply them.
