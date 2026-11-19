# ANTO v55 — Original Mota sunflower source reported bee-rate predictor, exact physical field heldout

**2026-11-19 after v54 strict source STOP; BEFORE any v55 seed–visit predictive readout.** This is a SEPARATE retrospective analysis, not a silent amendment of v54 or a new original study. v54's original failure is retained exactly in Draft PR #64 and docs/ANTO_V54_SOURCE_QUALITY_STOP_AND_NEW_ESTIMAND_20261119.md. The new question is whether the *directly recorded overall bee-visit count*, without using irreconcilable guild subtotals, adds heldout field-level information about seed set.

## Source provenance and the reason the prior version stopped
Mota et al. 2022 original Dryad DOI 10.5061/dryad.h18931zpg, Zenodo record 6671199, original author file MD5s in registry. ANTO v53 independently authenticated 80 common field-years, 51 true physical sunflower fields, 2,217 source harvest records and 11,682 one-minute visitor records, using exact original field names across files. ANTO v54 exact original-source integrity receipt Git blob SHA1 **55eb9b942e7a0ab7f4e32139ac5305325fb74a83** recorded:
- original No_Seeds_Total and No_rep_units are present/valid in **2,217/2,217** matched harvest observations;
- original No_Aborts_Total is MISSING in **1,056/2,217**, so the registered three-column identity guard of v54 exceeded its tolerance and halted, even though the recorded independent numerator and denominator exist. No imputing missing Abort values;
- original No_Interactions_TOTAL, No_Interactions_APIS, No_Interactions_WILD and Total_no_monitored_units are present in **11,682/11,682** matched visitor surveys, but TOTAL differs from APIS+WILD in **4,087/11,682**, largest absolute discrepancy **9**. Do NOT silently force the three to add, replace mismatches with 0 or discard contradictory surveys as favorable outliers.

This post-stop SOURCE QUALITY audit examined original numerical source field consistency but did NOT fit any seed/visitor predictive model; thus v55's model, response and readiness gate are frozen before any new predictive result. No source harvesting file rehosted.

## New tighter estimand, exact source admission

Same original Mota source field/year/region/type/distance controls as v53/v54; 80 field-years (2017 n41, 2018 n39), 51 physical fields (region and original name); source types NonV10/SNV32/WFS38; regions Burgos38/Cuenca42; 79/80 field-years with four aligned original distance bands, one with three.

For each original Field_code×Year×Distance band 0/15/30/60:
- y_d = original sum(No_Seeds_Total)/sum(No_rep_units), using ONLY the TWO original fields with complete source values, both nonnegative integers and No_Seeds_Total ≤ No_rep_units. No_Aborts_Total used only to AUDIT original source consistency when present; source blanks not assumed zero and not a v55 exclusion. All original 2,217 matched source harvest rows must remain eligible with positive denominator. Original source outcomes are harvested seed-set fraction, not seedling survival.
- x_total_d = original sum(No_Interactions_TOTAL)/sum(Total_no_monitored_units), where the latter is the observed sunflower units during each one-minute visitor census; each 11,682 original visitor row must be nonnegative and have positive exposure. The source author describes this as the direct total bee-interaction count but guild-specific totals disagree in many records: **do not infer visit taxonomy from the total**.
- x_reconstituted_d = sum(No_Interactions_APIS + No_Interactions_WILD)/sum(Total_no_monitored_units) as a SEPARATE **SECONDARY measurement-definition sensitivity**, *not* a valid duplicate measurement of TOTAL or a species-specific ecological effect. Include all 11,682 original visitor surveys without excluding incongruent rows.
- y_fieldyear, x_total_fieldyear, x_reconstituted_fieldyear are equal means across 3/4 original source-valid shared distance bands. This matches original location bands without weighting by 1-minute survey counts versus harvested seed numbers. One field-year is one modeled outcome/visit-rate unit; original plant rows/bee observations are NOT independent model sample units. Require all **80 fieldyears / 51 physical fields**, 2,217 harvest rows and 11,682 surveys, exact original source MD5s and quality receipts, fixed source taxonomy/inconsistency counts before model fit; any new source drift BLOCKS.

## Fixed exact leave-physical-FIELD-out prediction

For EACH of 51 original physical fields, hold out ALL of its year-2017/year-2018 field-year records together, fit on the other original fields; 80 source heldout field-years each scored exactly once. Every model predicts original seed fraction, linear fixed ridge alpha 8 on slopes, unpenalized intercept, clipped to [0,1] for MSE; report all raw prediction boundary clip counts.

- Z0: original region (Cuenca indicator), year (2018 indicator), and field management category (SNV and WFS relative to NonV) only; these are observational field contexts, not random treatment effects.
- Z1 PRIMARY: Z0 + training-only standardized log1p(original **reported TOTAL** bee interactions per monitored sunflower-minute).
- Z2 SECONDARY MEASUREMENT SENSITIVITY: Z0 + training-only standardized log1p(reconstituted APIS+WILD interactions per monitored sunflower-minute); no joint coefficient fitting or source guild additive biological assumption.

No respondent seed data from heldout fields used in model fitting, predictor normalization, choosing guild formula, hyperparameter or source field-year representation. Every original physical field must have ≥1 heldout field-year, all years absent together. Report equal-physical-field and fieldyear-weighted MSE, direction by physical field, original region/type summaries, absolute differences, and Z1 vs Z2 signed differences *without* postselecting best response predictor.

## Frozen visitor-placebo and interpretation boundary

200 deterministic visits-label permutations (seed 20261119) among original field-years only within source (Region×Year×Field_Type). Shuffle the Z1 and Z2 rates **jointly** to preserve original source reporting discrepancies, leave original field grouping, field metadata and measured seed set intact. Refit the same Z1 for each permutation, keep the complete set of physical-field holdouts and baseline Z0, report equal-physical-field gain reference median/q025/q975/max and fraction ≥ original; NOT a valid randomized ecological causal p value.

A positive field-transfer information verdict must meet ALL:
1. Z0−Z1 MSE gain >0 for both physical-field-equal and fieldyear-weighted means;
2. > half of 51 original physical fields improve;
3. original region-specific fieldyear mean gains positive in BOTH regions and source field-type mean gains positive for at least TWO of THREE original types;
4. source observed gain strictly above the 97.5% quantile of 200 stratified identity shuffles.

Otherwise use **VISITOR_INFORMATION_FIELD_TRANSFER_NOT_ESTABLISHED**, preserving all sources and negative results. A pass earns only **DESCRIPTIVE_SINGLE_ORIGINAL_CROP_STUDY_FIELD_TRANSFER**, not evidence of causal pollen deposition/realized outcrossing, independent multi-species universal law or individual sunflower identification.

## Limits, future tests and release

The source 51 original physical fields belong to one Mota 2017–2018 sunflower experiment, with some same fields repeated in the second year; exact field leaveout is not separate original published study replication. Both visits and seed set are observed under source management/region/year settings; mismatched monitoring and harvested sunflower individuals prevent one-plant inference, and breeding/genetic recombination/next-generation fitness are unmeasured. Mota authors already studied visitation and crop productivity: v55 contributes a different strict physical-field-heldout scoring estimand, not rediscovery of their original ecological effect. Sinapis seeds per plant and Ulex fruit fraction and Phlox log fruit counts cannot be numerically pooled.

Implement new versioned source-data parser with v53/v54 original exact MD5/fieldyear/quality receipt locks; no alteration to v54 prereadout code. Source-derived field-year/fold MSE only, 200 permutation results, adverse tests and reproducible read-only Actions. Draft stacked on PR #64, main unmodified, peer review required.
