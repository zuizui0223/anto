# ANTO v57 supplement — Mertensia original author bytes through DOI-authenticated Zenodo archive

2026-11-24. This is a **source-access recovery**, not a new original study or a new fit, and no predictor or holdout model is accepted before original author cohort validation. The original Dryad metadata is public for original Gallagher & Campbell 2021 *Mertensia ciliata* dataset DOI **10.7280/D16D7Z** (Dryad dataset numeric ID **59594**, published version **4**). Its TWO original files have exact original author names/sizes:
- `gallagher&campbell_frontiersinecology&evolution_data.xlsx`, **339,816 bytes**;
- `gallagher&campbell_frontiersinecology&evolution_metaDataFile.docx`, **47,265 bytes**.
Their canonical individual file download API URLs produce **HTTP 401** in GitHub Actions for both files; this is documented and remains an official Dryad-source access STOP.

A Zenodo repository search returned official record **4731053** with exactly the original paper/data title and the TWO original file names. Its search hit did not itself contain the Dryad DOI, so **a title match alone is INSUFFICIENT to authenticate**. Before accepting anything:
1. Require full official Dryad v2 original DOI metadata and exact publication version, 2 filenames, byte lengths, original published digest fields. No fabricated version/file IDs.
2. Require official Zenodo REST `https://zenodo.org/api/records/4731053`, exact full original publication title, same 2 original filenames and exact sizes. Verify publisher/creator context where available.
3. Download the two files only through the official Zenodo record's own HTTPS content URLs (not a GitHub mirror or third-party scraping), with 1 MB per-file cap. Validate source bytes against **BOTH** each original Dryad published digest (digestType-aware MD5/SHA256, if supported) AND Zenodo record's own checksum. No original data admitted with mismatched digests or unverified source link.
4. For the original Excel, report only workbook sheet names/first row field labels and worksheet approximate row counts; no original seed, visitation, water or flowering quantitative values read/analyzed. The docx is checked as an intact original office ZIP, not used to infer plant ID.
5. The original article's 114 seed-bearing plants do NOT imply 114 complete matched visit–seed source observations. Next separate prereadout protocol must verify same-plant keys, original blocks, missingness and effort before model fit.
6. If byte equivalence cannot be proven, fail closed as `SOURCE_EQUIVALENCE_STOP`, and explicitly preserve the previous official Dryad HTTP401. No interpretive biological null.

Physalis (2026) still has HTTP401 on ALL 10 original Dryad file streams and no verified DOI-identical Zenodo archive; cannot be scored. v57 is source-only; any individual visit/seed prediction requires a new frozen v58 model and new original author data CI. Retain all positive/negative outcomes from prior ANTO PR #66 as previously scoped and no merger into main.
