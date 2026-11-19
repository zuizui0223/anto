# ANTO v53 — Two independent original source audits reveal distinct limits of plant-level replication and an eligible sunflower field-year dataset

**2026-11-19. This is a source-identity and research-unit gate ONLY**, not a fitted ecological model, a negative pollination effect, or an independent verification of Sinapis visitor→seed predictive strength. The question was chosen after knowing Sinapis v52 successes. Exact source and source identity rules were committed before any new numerical outcomes were opened.

## Two verified original archives

Both studies were downloaded with the original publication-linked Zenodo record IDs and **all ten original files (MD5 authenticated)**. Only CSV headers, row counts, and sampling ID values/overlaps were inspected. No original seed production, visitor count, flower morphology or environmental effect coefficients were examined or fit. Original individual-level observations and files were NOT rehosted in ANTO.

| Original author experiment | Original source | Source files | Observed tables/grain | Direct same-individual visit→seed join? |
| --- | --- | ---: | --- | --- |
| Pardee et al. 2019 montane snowmelt/frost, original Dryad DOI 10.5061/dryad.kd1r166, archived Zenodo 4970401 | Original four wildflowers *Claytonia, Delphinium, Mertensia, Potentilla*, not four independent studies | 7/7 MD5 | 320 plot×species×year visitation records; 319 phenology records; per-plant seed set from 231+229+243+260 source rows | **NO**: visitation source has plot/species/year but NO plant ID |
| Mota et al. 2022 sunflower flower strips vs semi-natural habitat, original Dryad DOI 10.5061/dryad.h18931zpg, archived Zenodo 6671199 | Original sunflower crop field data from 2017 and 2018 | 3/3 MD5 | 2,875 reproductive productivity observations; 11,683 one-minute flower-visitor survey rows | **NO**: original visitor surveys at field/flower-unit level and distinct sunflower heads harvested; no same individual plant key |

**Do not claim** that these data replicate the same-plant Sinapis v52 visitor association, or that the four Pardee species represent four independent randomized study replications. In Pardee, plant-number tags in species-level seed-set files cannot be matched to individual plants in the species/plot/year visitation table. Aggregating visits to plot is an entirely different biological estimand, and the seed files do not expose an obvious year to support an unambiguous plot×year reproductive join.

## New strong admission for a different valid spatial prediction question: Mota original sunflower fields

The original Mota README specifically defines `Field_code`, `Field_ID` / `Field_code_original`, `Region`, `Year`, `Field_Type` and `Distance` in both visitor and productivity archives. It specifies one-minute visitor censuses and original total interactions / monitored sunflower units, plus reproductive seed set (full seeds vs full+empty+unfertilized florets).

**Original source-grain-only exact join counts:**
- Original productivity file: **102 field-year units** (52 in 2017, 50 in 2018) over 2,875 observation rows.
- Original visitor file: **81 field-year units** (41 in 2017, 40 in 2018) over 11,683 survey rows.
- Exact matching `(Year, Field_code)`: **80 field-year units** (41 in 2017, 39 in 2018), containing 2,217 original productivity rows and 11,682 original visitor surveys.
- All 80/80 source-matched field-years have internally unique and **identical original physical field name** (`Field_code_original` vs `Field_ID`) and agree exactly on original region and field type.
- **51 unique physical fields** identified from `(Region, Field_code_original)` across these 80 matched field-years. They are not 80 independent physical fields, and repeated field observations in different years must remain in the same holdout fold.
- The 80 matched field-years include original treatments `SNV` 32, `WFS` 38, `NonV` 10, and original regions Burgos 38, Cuenca 42. 79/80 have all four source distance-to-edge sampling bands in BOTH archives, and one has three common bands.

This evidence is sufficient to set up an independently motivated **FIELD / FIELD-YEAR spatial forecast experiment** in the next version, with explicit original biological-grain and inference limits. It does NOT yet demonstrate a join across individual flowers or causal bee pollination. Note that visitors are counted during 1-minute flower censuses, while seeds reflect later harvested crop reproductive units, which might not be identical sampled heads.

## Prior access stops preserved

The *Brassica rapa* six-generation Traine et al. (New Phytologist 2026) dataset appears structurally ideal for individually measured flower, visitor and mature seed information, but the user's EGWE repository already documented the exact official Dryad file stream **HTTP 403**, with no CSV row or header examined. The v53 selection registry preserves that source-specific denial and does **not** repeat uncredentialled access workarounds. Treat it as access stop, not a negative biological effect.

## Exact registered artifacts

- Source-first no-outcome protocol: `docs/ANTO_V53_INDEPENDENT_SOURCE_GRAIN_PROTOCOL_20261119.md`
- Original two-study 10-file expected-MD5 inventory and prior Traine stop: `data/external/anto_v53_independent_source_registry_20261119.json`
- Original 7+3 file header/key source pass: `scripts/audit_independent_pollinator_reproduction_sources_v53.py`
- Source-only original biological grain and field/year codebook audit: `scripts/audit_mota_fieldyear_source_codebook_v53.py`
- Two complete original-source outcome-blind receipts: `results/anto_v53_original_source_grain_audit_20261119.json`; `results/anto_v53_mota_fieldyear_source_key_audit_20261119.json`
- Negative tests against false same-plant joins: `tests/test_independent_original_source_grain_v53.py`
- Exact MD5+key receipt CI: `.github/workflows/anto-v53-original-source-key-scout.yml`

Source-only original Zenodo and codebook checks executed in GitHub Actions, source original values NOT modelled. **Draft stacked on ANTO v52 PR #62, main unchanged.** A new model must register its eligible field-year outcome and visit exposure, holdout whole physical fields (all years), conditioning, response scale and first falsification rule BEFORE examining numeric visitation/seed outcomes.
