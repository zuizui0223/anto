# ANTO v57 — Original third-study evidence search: authenticated Mertensia 115 plants and 20 intact array groups

**Source identity / eligibility audit, 2026-11-24. No third-study seed–visitor predictor fitted in this version.** This is a retrospective candidate selection after ANTO v52/v56, not independent outcome support. All raw author records remain in the original published archives, with source-derived hashes and group counts committed only.

## Published original sources, not inferred access or pseudo-replications

**Physalis philadelphica**: Davis et al. 2026 Oikos original paper DOI 10.1002/oik.11650, author Dryad dataset DOI 10.5061/dryad.xpnvx0kvp, official dataset ID 176194, published version **5**, exactly **10** original files. Official public DOI/versions/file metadata accessible. But all ten directly advertised official original file streams `/api/v2/files/<ID>/download` returned **HTTP 401**, and the separately documented **first-party original whole-dataset ZIP** endpoint also returned HTTP 401. No original author plant records, raw visitor values, plant-yield records or single-visit fruit-set values downloaded or fitted; **ORIGINAL_SOURCE_ACCESS_STOP** (not a biological negative). The published 2020 plant-season visit/yield records and distinct 2021 single-visit records must not be merged into one independent experiment, nor should originally defined `poll_imp` (which embeds fruit output) be used as a predictive visit variable. Source access receipts: `results/anto_v57_author_source_identity_grain_20261124.json` and `results/anto_v57_physalis_official_dataset_package_20261124.json`.

**Mertensia ciliata**: Gallagher & Campbell (2021), Frontiers in Ecology and Evolution DOI 10.3389/fevo.2021.641693, original Dryad DOI 10.7280/D16D7Z, official dataset ID 59594, published original version **4**. Official individual Dryad streams also returned HTTP 401, but official original Zenodo archived record **4731053** provided the original two author files. Both were accepted ONLY after double original-source checksum verification:
- `gallagher&campbell_frontiersinecology&evolution_data.xlsx`: original **339816 bytes**, original Dryad v4 SHA-256 **152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d**, Zenodo MD5 **db86f8286007b6176bf957334eba9ca4**;
- `gallagher&campbell_frontiersinecology&evolution_metaDataFile.docx`: original **47265 bytes**, Dryad SHA-256 **b41dd78899e1b451731917e923e739bbf1d0cbbd734867b1faa7d04713a8d968**, Zenodo MD5 **ed346bf60b6f0a2092a6cdc5efed4218**.
The original title/source version, both file lengths and two independent provider digests all match exactly. Thus original workbook is legitimately accessible despite the first-party per-file HTTP401. Source authenticity receipt `results/anto_v57_mertensia_dual_archive_author_grain_20261124.json`.

## Outcome-blind original manuscript-table grain and completeness

Original workbook has **11 original author sheets**. One source analysis sheet, `exp3_meansForAnalysis`, contains categorical `week`, `array`, `treat`, `plantID` and floral/visitor/seed fields. A source-only Excel XML parser checked **only nominal group/plant identifiers and presence/absence flags** for outcome/trait fields (not their magnitudes). The authoritative source result `results/anto_v57_mertensia_exp3_original_plant_grain_20261124.json` records:
- **115 original analysis-sheet plant records**, **115 distinct** `week×array×treat×plantID` keys, **0 duplicates** and **0 missing categorical keys**.
- **20 original week×array experimental clusters**, each **5–7 original plant rows**, across four weeks (28, 28, 30, 29) and three source water treatments A/C/D (39,39,37 rows). These 20 clusters are the candidate *heldout* sampling groups, NOT 20 independently published studies, neither phylogenetic replicates nor 20 independent evolutionary radiations.
- Original candidate column presence: `flrCnt` 115, `mcorW` and `mcorL` 106 each, `mVisRate` 115, `mseed` 114, `sCnt` 92. **105 records** have the proposed flower-count/corolla-width/corolla-length/visitation/seed fields present concurrently. The paper's ~114 reproductive-individual figure is not the 105-plant predictive complete-case cohort.
- Availability of water/phenology context and same-plant visitor plus reproductive measurements is a **source-level eligibility result**. Original published paper already argued that water availability affected seed production more strongly than visitation; no ANTO independent confirmation or contradiction is asserted.

## What is scientifically admissible next

A separate v58 prereadout contract can now specify `exp3_meansForAnalysis` same original 105 candidate plants, response `mseed` pending original author field semantics/numeric validation, flower number and corolla morphology, observed same-plant visit rate, and **intact week×array group holdout**. It must fix its original cohort, transformations, predictor leakage exclusion, family/treatment covariates, ridge penalty, baseline and failure rule BEFORE opening these 2021-author numeric observations. A whole-array experiment with same original 20 groups, independent of Sinapis and sunflower original datasets, would be a third botanical ORIGINAL STUDY, but still observational insect encounters and not an experimental causal pollen-receipt demonstration.

No original author raw XLSX, source plant IDs, seed or visitor counts were rehosted. Original author's seed-count biological meaning, visit-rate effort and missingness must be interpreted cautiously; original pollen deposition, seed paternity, multi-generation genetic persistence and ancestral flower colour regain are NOT measured.

## Reproduction / review
- `docs/ANTO_V57_INDEPENDENT_SAMEPLANT_ORIGINAL_SOURCE_GATE_20261124.md`
- `data/external/anto_v57_original_thirdstudy_source_registry_20261124.json`
- `data/external/anto_v57_verified_dryad_file_ids_20261124.json` (exact author version/file ID/name/bytes)
- `scripts/discover_original_Physalis_Mertensia_v57.py`
- `scripts/audit_original_Physalis_Mertensia_file_grain_v57.py`
- `scripts/audit_original_Mertensia_Zenodo_equivalence_v57.py`
- `scripts/audit_original_Mertensia_exp3_sameplant_grain_v57.py`
- `scripts/audit_original_Physalis_official_Dryad_dataset_ZIP_v57.py`
- `results/anto_v57_*20261124.json` (source-only immutable receipts)
- dedicated source-only read-only GitHub Actions (five distinct v57 workflows).

**Draft stacked on ANTO PR #66, no main merge.** v58's outcomes and scientific decision are not yet observed or tested here.
