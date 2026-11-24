# ANTO v57 — Third independent same-plant source admission, no outcome scoring

2026-11-24; this is an explicitly POST-v56 **SOURCE-ONLY** source-author eligibility gate. No original response/visit numeric values will be fitted, no p values, no biological replication inferred from source headers. Design a v58 frozen predictor+holdout after passing original source-grain checks and before outcome modelling.

## Established prior results and boundary
- ANTO v52: Sinapis original 2022 (Brassicaceae), same-plant visits+same-season seed estimates, visitor data improved source-treatment holdout 3/3. Observational visitation, not randomized effect.
- ANTO v56: Mota sunflower original 2022 (Asteraceae), source-field-year join only, new 51 physical-field heldout validation improved 34/51 with 16.54% relative MSE, BUT original visited sunflowers were NOT linked by plant-ID to harvested heads. V54/v55 source admission STOPS retained.
- Neither is a causal pollination mechanism or general across-species ecological law; new original plant-study systems and same-plant/whole-block holdouts are required.

## Original author candidate A — Physalis
Davis et al. 2026 Oikos, original paper DOI **10.1002/oik.11650**, Dryad original dataset DOI **10.5061/dryad.xpnvx0kvp** (published Feb 2026), Solanaceae *Physalis philadelphica* tomatillo. Original author file names from the public Dryad codebook:
`Oikos_2020_cross_season_abundance.csv` (plant block/treat/id, total flowers, observation dates, visitors),
`Oikos_2020_total_plant_yield.csv` (plant block/treat/id, fruit yield grams),
`Oikos_2021_fruit_set_poll_obs.csv` (single bee visit, plant_id/block/treatment, later fruit_set),
`Oikos_2021_seeds_fruit_sugar.csv`,
`Oikos_2021_poll_importance.csv` (derived pollinator importance incorporating fruit weight and fruit set),
`README.md`.
Study's 2020 plant-season and 2021 flower/single-visit are separate original cohorts **within ONE original publication**, not two independently published sources. Never use original `poll_imp`, `prop_fs`, `avg_fruit_weight` or fruit outcome as visitation predictors: they are calculated partly from target outcomes and leak the answer.

## Original author candidate B — Mertensia
Gallagher & Campbell (2021) *Frontiers in Ecology and Evolution*, DOI **10.3389/fevo.2021.641693**, original Dryad DOI **10.7280/D16D7Z** (Mar 2021, updated Apr 2021), Boraginaceae *Mertensia ciliata*. The original study manipulated water availability and flowering start, documented corolla/nectar traits, same potted plant visits and subsequent seed production. Original publication reports visitation–seed association not statistically decisive and stronger water effect on seed production; these are EXISTING author results, not ANTO v57 findings. Dataset original `gallagher&campbell_frontiersinecology&evolution_data.xlsx` and `gallagher&campbell_frontiersinecology&evolution_metaDataFile.docx`. Original paper seed N=114 is NOT assumed to be source complete-case or linked-plant N.

## Source author and exact original cohort gate
Use official Dryad public REST v2 documented DOI/versions/files and original download endpoints, or **official Zenodo archive whose DOI/title metadata unambiguously matches exact original author data**. No guessed github mirrors, source name substitutions, webpage scrapers, or unauthenticated access workarounds. Bounded original file downloads <=2 MB each, document official dataset version/file IDs, precise original file SHA256/MD5 and sizes. If download denied (403), metadata absent, or version/file linkage not resolvable, classify `ORIGINAL_SOURCE_ACCESS_STOP`, **not** an ecological null result.

**v57 only inspects original field names, original row counts and ID cardinalities, never numeric seed/visit outcome values.** Physalis 2020 direct original (block,treat,id) join between visitation and fruit-yield tables: verify source unique IDs, exact overlap, treatment/block consistency and original dates/exposure. Physalis 2021 source plant_ID, treatment, block, fruit_set and bee identity/behavior source headers and source block×plant group counts; avoid duplicating flower rows as independent plants. Mertensia XLSX first row labels and sheet inventory using ZIP/XML standard library without calculating response values; metadata docx source verified but not interpreted as equivalent data table.

Minimum admission for follow-on original study: >=30 directly same-plant linked visitor and subsequent reproduction IDs, >=3 independent original field blocks/arrays with at least 5 linked original plants each; effort denominator or matched visit duration; original outcome not derived using visit counts; group identifiers and original time semantics verified; no group leakage. If these do not hold, do not fit/claim v58. Even if pass, v58 preoutcome protocol must freeze models and holdout units and then run original-model CI. Source individual rows never committed.

Draft stacked on ANTO v56 PR #66; main unchanged. No quantitative meta-analysis of the existing Phlox, Ulex, Sinapis and sunflower MSEs because outcomes/sampling units differ.
