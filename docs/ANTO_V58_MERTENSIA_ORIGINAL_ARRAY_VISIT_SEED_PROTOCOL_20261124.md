# ANTO v58 — Original Mertensia same-plant flower/visitor/seed prediction under whole-array holdouts

**2026-11-24, frozen BEFORE reading numeric original seed output or visitation rates.** This is a post-v57, retrospective study-design question, informed by existing original author Gallagher & Campbell (2021) report where plant seed output responded to water availability and a visitor–seed observational association was not compelling. A failed predictive visitor-gain test is scientifically admissible, and author findings are not independent ANTO replications. No new raw original individual data are committed.

## Precisely authenticated original biological source

Gallagher & Campbell (2021), Frontiers in Ecology and Evolution DOI 10.3389/fevo.2021.641693, original author Dryad DOI **10.7280/D16D7Z** (official dataset ID 59594, published version 4). The individual official Dryad source API routes returned HTTP401, but original Zenodo archived record **4731053** held the exact same original author's two files as confirmed by BOTH original original-archive checksums:
- Original xlsx `gallagher&campbell_frontiersinecology&evolution_data.xlsx`, **339816 bytes**, Dryad SHA256 **152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d** = SHA256 of Zenodo source, and original Zenodo MD5 **db86f8286007b6176bf957334eba9ca4**.
- Original metadata docx 47265 bytes, Dryad SHA256 b41dd78899e1b451731917e923e739bbf1d0cbbd734867b1faa7d04713a8d968 and Zenodo MD5 ed346bf60b6f0a2092a6cdc5efed4218.
Do NOT substitute a title-only similar file, GitHub copies, edited author workbooks or a synthetic test cohort for original data.

Parent ANTO v57 source-grain author result (Git blob SHA1 **TO_PIN_BLOB_BEFORE_OUTCOME**) audited **115** distinct `week×array×treat×plantID` source records from `exp3_meansForAnalysis`, **20 original week×array groups** of 5–7 original plants, weeks (1..4) and water treatment A/C/D. Candidate fields `flrCnt` (115 present), `mcorW` / `mcorL` (106 each), `mVisRate` (115), `mseed` (114). Exactly **105** plant records have all five proposed predictors+reproductive field PRESENT. Presence does not guarantee numerical validity: a raw formula error, nonfinite response, negative flower/visitor/seed count or contradictory sheet must **TECHNICALLY BLOCK** rather than be silently excluded. The paper reported ~114 seed-related plants, NOT the same analytic N=105.

## Original source observation grain, fixed cohort gate, transformations

Use ONLY 105 complete numeric plant-level rows selected by all source fields `flrCnt`, `mcorW`, `mcorL`, `mVisRate`, `mseed` plus original categorical `week`, `array`, `treat`, `plantID`. Source numeric source entries must be finite and nonnegative, and `mcorW`/`mcorL` positive. If any candidate record fails numerical admission, reject instead of deleting to improve a positive score. The original `mseed` variable is analyzed as **source-coded seed-production proxy**; do not describe it as genetically inherited fitness or claim per-flower units unless the original metadata confirms. Response `log1p(mseed)`. Flower supply feature `log1p(flrCnt)`, original corolla widths/lengths `mcorW`, `mcorL`, visitation feature `log1p(mVisRate)`. These are all recorded contemporaneously; the visitor predictor is NOT experimental visitor manipulation or a prospectively deployable seed-set forecast.

Use source nominal water treatment levels A,C,D as labels (do NOT guess biological meaning from letters), original weeks 1,2,3,4, original arrays 1..5. Group holdout ID = **(week,array)**, all original plants of that week-array held out together regardless of treatment or shared plantID strings. Require original 20 groups remain represented after complete-case selection and each has ≥3 complete plants; require source 105 complete plants in exactly 20 unique groups and total original 115 candidate grain. Original plant key duplicates technically reject. Do not interpret 20 groups as 20 independent published studies, independent genetic populations or independent evolutionary radiations.

## Fixed predictive model

Ridge alpha=**5.0**, unpenalized intercept, source original log1p seed proxy prediction with raw MSE (nonnegative response no clipping because response transformed log-positive; track raw negative predictions but don't select score based on clipping). Nominal factors as reference dummy indicators (source `C`,`D` vs original A; weeks 2/3/4 vs 1), never fold-dependent target labels. Numeric features standardized with training-only mean and population-SD, rejecting degenerate SD; no across-fold visitor imputation and no response encoding in predictors.

- `Q0` baseline: intercept + week and water-treatment dummies + `log1p(flrCnt)` only.
- `Q1` source floral morphology: Q0 + original `mcorW` and `mcorL` (two slope additions).
- `Q2` PRIMARY visitor increment after morphology: Q1 + independently recorded source plant `log1p(mVisRate)`.
- `Q3` SECONDARY: Q0 + same plant visitor rate, without corolla predictors; do not promote after viewing results.

Exactly **20 whole-original-week×array heldouts**, each of the 105 complete original plant records predicted exactly once. Per fold: original group n, original week and water treatment counts, training n, full Q0/Q1/Q2/Q3 MSE and signed Q1−Q2, Q0−Q1, Q0−Q3. Report equal-20-group MSE and 105-plant-weighted MSE. No row-level random splits, in-sample AIC tuning, hyperparameter selection, and no filtering of negative groups.

## Fixed observational visitor identity placebo

200 RNG seed **20261124** permutations of the observed plant `log1p(mVisRate)` predictor within each original (week,treat) stratum, preserving water/wk factors, flowers, corolla traits, seed-proxy responses and all 20 whole-array holdouts. Refit Q2 under each placebo; original Q1 remains the fixed baseline. Report actual equal-group Q1−Q2 gain vs 200 placebo median/q025/q975/max and number ≥ actual. This is a noncausal source plant-label comparison, not experimental randomization of pollinator visitation or a p value.

**Pre-outcome descriptive information-transport gate**: Q1→Q2 equal-group and plant-weighted **absolute MSE gain both >0**, ≥12 of 20 original heldout array groups improved, ≥3 of 4 week-specific group-equal mean gains >0, and observed equal-group gain exceeds placebo q97.5. If any fail: `THIRD_ORIGINAL_MERTENSIA_VISITOR_TRANSFER_NOT_ESTABLISHED`. Passing means ONLY `DESCRIPTIVE_ORIGINAL_MERTENSIA_ARRAY_HELDOUT_VISITOR_INFORMATION`, not causal pollination benefit or universal reproduction mechanism. Opposite results relative to Sinapis/Mota are not erased or pooled; they may indicate context-specific information or missing field-grain comparability.

## Scientific limits / reproducibility

Original author's watering treatment may be manipulated, but source plant visitation is observational and seed production was not a randomized bee receipt intervention. Experimental array assignments share field garden and source donor cohort; they are not independent original studies. Flower traits and visit rate might be linked through unmeasured phenotype/plant vigor; same-source seed proxy and visit features may share time windows. Original plant identity, original metadata and year must not be invented. Some flower species have heterogeneous pollinator guilds: `mVisRate` total source rate is NOT automatically stigma contact, compatible pollen, outcross paternity, seed recruitment or multi-generation hereditary fitness. Inference is limited to original site/arrays, not ecological universality or ancestral pigment trait return.

Script must retrieve original author workbook through original Zenodo 4731053, check both Dryad SHA256 and Zenodo MD5, parse original Excel worksheet XML `sheet11.xml` directly and detect numerical/formula errors before any Q0..Q3. Import neither v57 source-only outcome results as model effect nor v52/Sinapis coefficients. Freeze full derived-only JSON and exact model parity in read-only original-source GitHub Actions, plus ≥12 adversarial tests for clone/key/group leakage, identity, missingness, model feature source, placebo, and negative results.

**Draft stacked on ANTO v57 PR #67; main untouched** pending original author-source CI and independent scientific review.