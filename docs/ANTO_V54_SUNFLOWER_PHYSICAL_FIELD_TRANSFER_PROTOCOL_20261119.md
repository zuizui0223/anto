# ANTO v54 — Original sunflower crop: field-heldout bee visitation versus seed-set prediction

2026-11-19, prereadout prediction contract frozen before opening any Mota original visitation or seed outcome value. Scientific QUESTION is posthoc and follows ANTO v52 Sinapis visitor results and v53 source-grain gate. This is one independent original Mota et al. 2022 field study relative to Sinapis, Phlox and Ulex, NOT 51 independent original published studies or an intervention on bees.

## Immutable original source and grain

Mota et al. (2022), Flower strips and remnant semi-natural vegetation have different impacts on pollination and productivity of sunflower crops, original data DOI 10.5061/dryad.h18931zpg, archived Zenodo 6671199. Three source author files, exact original MD5:
- Database_Productivity.csv = aeea30a21e495498295241bcfd8c693d
- Database_Visitation_rates.csv = e13131ad60b133cc6ab25e20bbd3554b
- README_file.txt = 479b26d3cebb72c8491b0ec00609c08c

The SOURCE-FIRST v53 gate obtained and verified these all, without reading numerical outcomes. Frozen derived v53 original physical-field/key receipt at results/anto_v53_mota_fieldyear_source_key_audit_20261119.json, Git blob SHA1 7a4fb27c7f26abfa8400d23b041e1a3d3e3424c9; original 80 shared field-years (2017 n41, 2018 n39) and **51 physically unique sunflower fields** from region + original name (many physical fields surveyed both years). All 80 shared keys (Year,Field_code) have exact agreement between productivity Field_code_original and visitation Field_ID, agree on Region and Field_Type, and share 4 original distance-to-field-edge bands for 79 field-years and 3 bands for one. Source eligibility: 102 productivity field-years (2875 original harvested units), 81 visitation field-years (11683 one-minute survey rows); 80 shared field-years encompass 2217 productivity units and 11682 surveys. Original field types WFS 38 / SNV 32 / NonV 10, original regions Cuenca 42 / Burgos 38.

NO individual sunflower plant ID is common across visitor survey and harvested seed tables. Every predictive claim here is **FIELD-YEAR** ecology, not same-plant visitation/seed association. Explicitly exclude original 22 productivity field-years without visitation and 1 visitation field-year without source harvest; these are not no-visit/no-seed biological zeros. Never silently join the files by row number or treat 11682 census rows as individual independent plant replicates.

## Predeclared original band-aligned ecological variables

Valid original rows are in exact matched Field_code×Year×Distance bands 0,15,30,60 metres from field edge. Source original README defines:
- productivity No_Seeds_Total = full seeds, No_Aborts_Total = empty seeds and florets, No_rep_units = full+empty+unfertilized florets, and Seed_set = full seeds / No_rep_units;
- one-minute visitor census Total_no_monitored_units is the number of sunflower flower units observed for one minute, No_Interactions_TOTAL = Apis+wilder bee interaction contacts, separate No_Interactions_APIS and No_Interactions_WILD.

Hard numeric/data checks:
- columns/exact original names, field/year/region/type/site-name and distances fixed, no strings outside 2017,2018 / Burgos,Cuenca / NonV,SNV,WFS or distance 0,15,30,60;
- all present numeric source counts finite and nonnegative, seed/abort/reproductive-unit and Apis+wild arithmetic consistent, positive exposure denominator;
- explicit source cells blank, NA or dot classified missing and excluded only from that numerator/denominator pair, never imputed as zero; unknown nonnumeric, negative or contradictory count is a fail-closed source quality stop; report skipped counts and do not exceed 5% exclusions of the 2217 harvest and 11682 visitor matched source rows;
- for each original field-year and each common original distance band compute y_d=sum(No_Seeds_Total)/sum(No_rep_units) of source harvested units with valid positive denominators, and x_total_d=sum(No_Interactions_TOTAL)/sum(Total_no_monitored_units) of valid one-minute survey flower exposures, similarly x_apis_d, x_wild_d;
- average each **distance band's** y_d and visit rates equally across the 3 or 4 verified common bands in that field-year (don't weight by harvested seed units or survey count across bands). Require >=3 valid shared bands per retained field-year. This creates exactly one FIELD-YEAR y, x, and region/year/field-type record, with NO individual flowers or pollinator rows as independent prediction units.
- target field-year y is the original harvested seed-set fraction in [0,1], not log seed counts, seed viability, next generation survival or genetic evolution; predictors observed during same original field-year, not known before season.
- admission requires >=70 original matched field-years and >=45 unique physical fields, each of the two regions, two years and three field types represented.

## Primary model, holdout and outcome BEFORE opening original responses

Keep **all field-years of a given physical field (Region,Field_code_original) together** in one heldout group; absolutely no 2017/2018 leakage of the same physical field. Leave ONE of the >=45 physical fields entirely out, train on all remaining physical fields and predict each of the heldout field's eligible field-years. This is a single original field study, with repeated field-year measurements; report both physical-field-equal and field-year-equal MSE, source group sizes and sign, but do not call 51 independent experimental studies.

Source-normalized fixed bounded-risk models (all slopes fixed ridge lambda=8, intercept unpenalized):
- B0: original Region, Year and Field_Type categorical main-effect indicators, NO visitor information; treatment reference NonV, region Burgos, year 2017;
- B1: B0 + training-only-standardized log1p(original source total recorded bee interactions per monitored flower-minute), exactly one source field-year-level x computed by averaging aligned distance bands;
- B2 secondary: B0 + separately training-standardized log1p(Apis-only interactions per monitored flower-minute) and log1p(wild-bee interactions per monitored flower-minute), not selected posthoc over B1.

No field/year/region treatment code from a heldout field's seed outcomes enters training, and training-only standardization of visit predictors is re-fit for each physical-field holdout. All models use source original physical field-year response as a fraction (not per-seed Bernoulli pseudo observations), raw squared-error after predictions clipped into [0,1]. Track numeric raw-score clipping 0/1 counts separately for all B0/B1/B2. The primary visitor information contrast is MSE(B0)−MSE(B1) for each heldout physical field, averaged equally across physical fields. Secondary B0→B2 and B1 vs B2, original site/year/type heterogeneity reported without reselecting the primary endpoint.

## Fixed source-placebo contrast (not randomized pollinator intervention)

Use 200 deterministic seed=20261119 permutations of the original field-year total bee visitation predictor **jointly WITH Apis and wild** across other field-years within the same original Region×Year×Field_Type stratum, preserving physical field holdout groups, source seed outcomes, year, landscape category, source distances and all other covariates. If a stratum has <=1 field-year, its visitation rates stay unchanged. Refit B1 on 200 permuted datasets under the SAME physical field holdout. Summarize source equal-physical-field MSE increment distribution (median, 2.5% /97.5%, maximum and n >= observed); this is observational within-stratum identity placebo, NOT random-assignment causal p-value or valid independent field-year randomization inference.

**Predeclared descriptive cross-field information gate:**
(1) B0→B1 mean MSE reduction strictly positive in both equal-physical-field and equal-fieldyear weights;
(2) B1 better than B0 for MORE THAN HALF of the original physical fields;
(3) positive mean source group MSE gain within both regions and at least two of three original Field_Type categories;
(4) observed mean physical-field MSE gain exceeds the fixed 200-placebo q97.5.
If any fail: VISITOR_INFORMATION_FIELD_TRANSFER_NOT_ESTABLISHED. If pass: DESCRIPTIVE_SINGLE_ORIGINAL_CROP_STUDY_FIELD_TRANSFER. Never call a p-value or causal seed-setting mechanism.

## Scientific boundaries

The original author paper already studied managed flower strips, seminatural vegetation, insect visitation and crop productivity; this reanalysis asks a narrower field-heldout predictive question. Bee interaction frequencies are not pollen deposition or seed paternity; season-end harvested sunflowers are not necessarily the same plants watched in visitor censuses. Location/year/field management confounding can remain even after B0. The 51 original physical fields and two years arise from one original study in two regions, not 51 research publications nor naturally independent evolutionary radiations. Outcomes are seed set FRACTION, Sinapis v52 had estimated seed NUMBER and a different measured plant-level exposure, Ulex fruit fraction and Phlox fruit counts: DO NOT POOL numerical MSE effects.

Source original first MD5 and grain auth, fail closed on numeric domain, retain all positive and negative original complete field holdouts, source-derived aggregated-only JSON, exact-head CI, DRAFT PR stacked on v53 PR #63. Do not touch main or change registered models/outcomes to chase positives.
