# ANTO v59 — disentangle within-season array transfer, temporal extrapolation, and measured soil moisture

2026-11-26. **Retrospective diagnostic of the same 2016 Mertensia experiment, not new biological replication.** Contract fixed before v59 refit/readout; v58 original outcome already known. Preserve every v58 result, including negative week 3. Original article DOI 10.3389/fevo.2021.641693 describes a key design limitation: 12 plants ceased flowering during week 4 and the remaining flowering plants were rearranged into three new arrays. Their contexts were shared across the nominal original array labels. Thus a nominal week×array split does not guarantee disjoint environmental context for week 4, although no actual duplicated response leakage has yet been shown.

## Source and source grain

Reuse precisely the 105 complete plants and original log1p(mseed) response of v58; source author Dryad version 4 SHA256 + Zenodo 4731053 MD5, original v58 script and receipt Git blobs frozen in v59 registry. First independently regenerate all 80 v58 model-fold losses and compare with the frozen receipt (absolute tolerance 1e-10). Check original raw plantID counts across weeks, arrays and treatments (no actual tokens exported). If raw plantID labels occur in different weeks, stop new analysis pending source semantics; if same ID occurs in different arrays in weeks 1–3, stop safe16; week 4 is always contained. For the treatment split, the raw plantID must not occur in multiple treatment groups. A shared array/environment does not imply genetic independence.

Re-read the original author data dictionary docx (exact original double archive checksums) and collect only bounded field-definition snippets for mseed, mVisRate, mVWC. Do not invent units. Article defines seed set per flower and per-flower visitation; retain v58 source-coded outcome until dictionary corroborates field semantics. Match mVWC from the SAME original sheet row by exact composite key for every v58 eligible plant. All 105 mVWC values must be finite and between 0 and 100; no imputation, outcome-driven exclusion or cohort change. Use original water treatment codes A/C/D without assigning labels unless source verified.

## New fixed comparisons

Base predictors B: unpenalized intercept, original treatment dummies C,D (omitted ONLY in whole-water-treatment splits), train-standardized log1p(flrCnt), mcorW, mcorL, plus the specified contextual terms. V = B plus train-standardized original log1p(mVisRate). All non-intercept coefficients ridge alpha=5; raw log-response squared errors, no clipping. No outcome tuning.

For each of these contexts, score all frozen 105 plants once:
- no_week (secondary): no week predictors.
- linear_week (primary): numeric week 1–4, train-standardized. It can extrapolate to a withheld week, unlike an unseen categorical dummy; no target-week outcome enters training.
- linear_week_vwc (secondary sensitivity): same linear week + independently measured original mVWC, standardized on training only. This is contemporaneous resource context, not a causal adjustment or prospective pre-season forecast.

Each context uses three holdout policies:
- array20: nominal 20 original week×array folds, diagnostic including possible week-4 environmental overlap.
- safe16: hold out each of 15 original arrays in weeks 1–3, plus the entire fourth week as a single fold, all its 23 eligible plants together.
- week4: ALL FOUR whole phenology-week folds (not merely week number 4), each week completely excluded from training. Primary comparison B vs V in linear_week.
- water3 (secondary): whole original water-code holdouts, three folds; fit linear_week and linear_week_vwc with NO water-treatment dummies because the category is unobserved in training. Other plants from the same arrays can occur in training; this is treatment-context transfer, not spatial independence.

Report absolute losses and gains with equal-fold and plant-weighted averages, and mean of four per-week plant-weighted losses. For each week and treatment, report n, B and V MSE, signed gain, and count test plants outside each training predictor range. Do not compare unequal fold-weighting without displaying the plant/weekly common weights. Training-fitted coefficients/response corrections are descriptive, not causal selection slopes.

## Diagnostic and falsification

Replay v58 Q1/Q2 in every original fold and compute exact arithmetic loss decomposition per week: e=y-p_B, d=p_V-p_B, gain=2 mean(e*d)-mean(d^2). The alignment term and squared correction separate a wrong/weak directional prediction from an excessive correction. Same diagnostics in the new folds. This is not a mechanistic partition of reproductive biology. Store all signs, no selecting favorable week or model.

Primary directional gate: whole-week linear_week both week-equal and plant-weighted gains >0 and at least 3/4 weeks positive. Separately compute 200 fixed-seed (20261126) within original week×treatment plant-visit identity shuffles, refit V in the four whole-week folds, unchanged B, record actual gain versus q97.5 and number equal/exceed. For stronger descriptive transfer require both directional gate and actual strictly above q97.5. Neither shuffle nor four correlated folds constitute a randomized visitation intervention, valid inferential p value, or a confidence interval. Failed gate is retained even if a secondary model performs better.

Measured VWC may change a predictive comparison; it does not identify plant resource limitation causally. Seed set, effective pollen deposition, compatible pollen, genetic inheritance and future recruitment are different quantities; the same study cannot establish universal pollination effects. No raw author rows or plant ID tokens committed. Dedicated read-only source CI, derived full results and Draft stacked on PR #68, main untouched.
