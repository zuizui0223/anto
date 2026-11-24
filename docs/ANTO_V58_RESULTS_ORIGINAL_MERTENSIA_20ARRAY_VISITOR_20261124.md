# ANTO v58 — Third original independent plant system: same-plant visitor information adds array-heldout seed-proxy prediction, but NOT every week

**2026-11-24; one new original published study relative to Sinapis and Mota, NOT 20 independent papers, no new independently collected data and NO causal visitation experiment.** The scientific hypothesis was chosen retrospectively after positive ANTO v52/v56. However this v58's fixed predictors, nonadmission gates, exact original sheet/cohort, 20 group holdouts, 200 source-identity placebos and descriptive pass/fail rule were committed **before reading original numerical visitor or seed outcome values**.

## Authoritative original source and actual biological unit

Gallagher & Campbell (2021), *Frontiers in Ecology and Evolution* DOI **10.3389/fevo.2021.641693**, original data DOI **10.7280/D16D7Z**. Original Dryad dataset 59594 version 4, files authenticated from original archived Zenodo record 4731053 by **both** Dryad original SHA256 and Zenodo MD5 (not a title-only mirror). Original exact Excel `gallagher&campbell_frontiersinecology&evolution_data.xlsx`, size **339816** bytes, original SHA256 **152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d**. The original v57 author-source audit is frozen at Git blob `5659545f1376aeb715aefa4308a0fe1e8673da39`.

Single author original `exp3_meansForAnalysis` sheet of the original 2016 experiment:
- 115 original distinct `week×array×treat×plantID` records (0 duplicate source keys);
- **105 same-plant COMPLETE numerical records** for original flower count `flrCnt`, corolla width `mcorW`, corolla length `mcorL`, plant visitation `mVisRate`, source-coded seed output `mseed`;
- **20 whole experimental week×array clusters**, 4 original weeks, 5 arrays per week; **4–6 eligible plants per array** (original full-sheet groups had 5–7).
- Water-treatment levels retained as original source codes A, C, D, without invented wet/dry labels.
- The source-coded `mseed` is used as a seed-production proxy transformed `log1p`. Do NOT claim genetic selection, offspring viability, seeds-per-flower physical units or per-visit pollen transmission from this variable alone.

Original author publication reported water availability as a stronger determinant of seed production and a not-clearly-significant observed visitation–seed relation. We do NOT reinterpret that published statistical result as false: the new ANTO estimand is out-of-array conditional prediction across the existing source cohort.

## Frozen predictive comparison, whole original array out

20 group folds exclude **all plant rows of one original (week,array) experimental cluster**, fit on the other 19 clusters. All 105 eligible source plant observations predicted once, no original plant ID or outcome leaks across training, all numeric features standardized using ONLY train rows and ridge alpha=5 on slopes (unpenalized intercept). Original factor controls are water-code A/C/D and week1–4, plus log1p original flower count.

| Source model | Original predictors beyond water/week/flower count | Original 20-group-equal MSE (lower is better) | 105-plant weighted MSE |
| --- | --- | ---: | ---: |
| Q0 | None | **0.169035110** | 0.158917300 |
| Q1 | Original corolla width and length | **0.175770351** | 0.165542350 |
| **Q2 (primary visitor increment)** | Q1 plus source plant visitation rate | **0.156447157** | 0.146388569 |
| Q3 (secondary alternative) | Original visitation rate without corolla widths | **0.151260952** | 0.141307934 |

- Adding corolla morphology Q0→Q1 **worsens mean source-heldout MSE by −0.006735241** (negative gain), source arrays 10 improve and 10 worsen; this is not evidence that corolla morphology has no evolutionary function.
- Adding independently recorded original same-plant visitor rate **Q1→Q2 improves heldout MSE by +0.019323194 in equal-array weighting**, an **11.0% reduction relative to Q1's original MSE**; +0.019153781 on original eligible plant-weighted loss. **16/20 original arrays improve, 4/20 worsen**.
- The simpler Q3 has the lowest overall array mean MSE, but was predeclared SECONDARY; do **not** change the primary readout to Q0→Q3 after results. Original Q0→Q3 absolute gain is +0.017774158 equal-array, positive 15/20 and negative 5/20.

## Source week heterogeneity and negative source groups retained

Primary Q1−Q2 equal-array original 5-group mean per original week:
- Week 1: **+0.044686008** (5/5 arrays positive)
- Week 2: **+0.003249530** (3/5 positive)
- Week 3: **−0.005270504** (**3/5 positive but mean negative** due to negative losses)
- Week 4: **+0.034627742** (5/5 positive)

This is **3/4 original weeks with positive mean**, not a uniform visitor-transfer law. In 4 individual original arrays, adding visits worsens risk. The negative week/arrays remain in the frozen result without selective omission.

Separately marked **AFTER-OUTCOME arithmetic influence** diagnostic (no new source model refit, no inferential CI): delete exactly one original group from the frozen 20 group MSE score sums and recompute mean Q1−Q2 risk gain. **All 20/20 arithmetic group deletions keep the aggregate mean positive**, minimum **+0.0159986384**, maximum **+0.0222578920**. This shows the observed positive mean is not caused by exactly one array, but cannot establish true independence across phenology weeks or new ecology studies.

## Visitor plant-ID placebo (not experimental intervention)

Source visit-rate predictor `log1p(mVisRate)` shuffled 200 times with fixed seed 20261124 **within original week×water code** among eligible plants, all arrays/seed responses/corolla measures/flower counts and holdout groups kept unchanged. Refit Q2 after each shuffle, compare to the unchanged original Q1 baseline.

- Observed equal-array Q1−Q2 gain: **+0.019323194**.
- 200 source-label shuffled median **−0.001572113**, q97.5 **+0.006553548**, largest **+0.016280857**.
- **0/200** source visitor-label relabelings produced an equal or greater gain. This is a descriptive predictor–response identity reference, **NOT** a causally valid pollinator randomization p value (original water/phenology treatment does not randomize actual visits).

The **preoutcome v58 descriptive gate** (positive both weighting measures, at least 12/20 original groups better, at least 3/4 source weeks positive mean, observed > source-placebo q97.5) **PASSES** exactly as frozen: `DESCRIPTIVE_ORIGINAL_MERTENSIA_ARRAY_HELDOUT_VISITOR_INFORMATION`.

## Interpretation across ANTO v52 Sinapis / v56 Mota / v58 Mertensia

Three **distinct published original studies** contribute narrow, related support for the empirical question of **observed insect interaction information** versus static morphology in heldout reproductive prediction:
- Sinapis 2022, 52 same-plant samples, three unseen experimental rain-treatment targets, source **estimated seeds per plant** (original v52 +0.07857 absolute MSE conditional visitor gain, positive 3/3);
- Mota sunflower 2022, **51 original field-group** targets and 80 field-years, 2214 source harvest units + 11682 bee censuses aggregated across matched source distance bands, original source **seed-set fraction** (original v56 +0.002243 field-equal MSE, positive 34/51), **NOT individual plant-level visit/seed link**;
- Gallagher Mertensia 2021, **105 linked original plants** and 20 original week×array experimental groups, original **log1p(source-coded mseed)** seed proxy (v58 +0.019323 equal-array MSE, positive 16/20, but week3 source mean NEGATIVE).

These data units, holdout domains, model baselines, source endpoints, trait pools and field settings differ. **DO NOT pool raw MSE deltas or present 3/3 studies as three independent identical randomized ecological effects.** The Mota source required exactly three author-impossible seed-ratio exclusions after earlier v54/v55 source-quality STOPS, which remain documented; no plant-visitor causal effect or pollen deposition established by the Mota data.

The strongest common interpretation is **conditional out-of-domain predictive information from directly observed interaction rates can persist where floral morphology alone does not improve risk**, but is heterogeneous across environmental weeks and source groups. There is no reliable estimate of universal natural selection, no experimental randomized bee mediation, no offspring paternity, pollen compatibility measurement, lineage persistence or ancestral petal colour reacquisition. A truly general ecological evolutionary claim requires independent additional original systems and **aligned intervention-to-compatible-pollen-to-descendant endpoints**.

## Reproducibility and review

- Before outcome readout frozen protocol: `docs/ANTO_V58_MERTENSIA_ORIGINAL_ARRAY_VISIT_SEED_PROTOCOL_20261124.md`
- Immutable original SHA/MD5, original 105 source cohort and fixed model gate: `data/external/anto_v58_mertensia_original_source_model_registry_20261124.json`
- Original author XLSX downloaded through DOI-verified, **double** original archive digest, exact source sheet XML parse, original 20 whole-array folds, Q0–Q3 and 200 placebos: `scripts/analyze_original_Mertensia_20array_visitation_v58.py`
- Original group/key/leakage and source-only tests: `tests/test_original_Mertensia_20array_v58.py`
- Complete 20 original per-array original-derived MSEs, all response/predictor provenance and placebo, **NO original author seed/visitor raw plant rows**: `results/anto_v58_mertensia_array_heldout_visitation_20261124.json`
- Read-only original-author source repeat: `.github/workflows/anto-v58-original-mertensia-20array.yml`

First actual original-author source CI run `45380516060` completed SUCCESS with original double-digest/cohort tests and full 20 original array model first readout. **Latest exact-head stored-result CI remains a further gate**. Keep Draft stacked on ANTO v57 PR #67, main unchanged and no scientific peer review presumed.
