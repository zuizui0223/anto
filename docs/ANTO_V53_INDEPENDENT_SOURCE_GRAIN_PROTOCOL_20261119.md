# ANTO v53 — Independent original pollination–reproduction source-grain gate

2026-11-19. After ANTO v52 original *Sinapis arvensis* whole-rainfall-treatment heldout predictive success, the next task is to identify **truly independent original source studies** containing linked floral visitation and reproductive outcomes, not endlessly resplit Sinapis. This is a **source-only, outcome-blind eligibility audit**. It does NOT estimate any plant-effect coefficient, fit a visitor–seed model, claim positive or negative biological replication, nor alter the prior ANTO PR #62.

## Previously known source and exclusions
- Sinapis v52: one original Hoefer et al. 2022 study, 52 eligible plants and 3 entirely heldout rain treatments, visitor information prediction improved 3/3, observational even though rainfall was manipulated.
- The apparently ideal Traine, Rusman & Schiestl New Phytologist 2026 *Brassica rapa* six-generation experimental evolution archive (Dryad DOI 10.5061/dryad.2ngf1vj15) was ALREADY officially audited under zuizui0223/egwe (source gate recorded in artifacts/empirical/traine_plant_assay_source_20261008.json and docs/TRAINE_PARTIAL_OPERATOR_SOURCE_AUDIT_2026-10-08.md): latest source version 429932, file ID 4644044, 288058 bytes, source archive digest SHA256 4c9162866b4d90afd31c847e19e9de39514ab168088de4ba57e7f3c70d37af77. Exact identified Dryad file-stream returned HTTP 403 in GitHub Actions run 37788113341. **No original source row/header access and no biological test.** Do not repeat unauthenticated URL workarounds or call it new.
- A separate recently found Dryad Oikos herbivory/bee single-visit study, https://doi.org/10.5061/dryad.xpnvx0kvp, advertises original same-plant visits and fruit-set records but has no VERIFIED source stream currently. It is a candidate, not an admitted biological replication.

## Registered source-first new audits (original archive, never outcome modeling)
**A. Pardee, Jensen, Inouye & Irwin (2019)** experimental snowmelt/frost × pollination and seed-set across four subalpine wildflower species, original data DOI 10.5061/dryad.kd1r166, archived Zenodo record 4970401. Author source files: AllSpeciesVisitation.csv, AllSpeciesPhenology.csv, four {Claytonia,Delphinium,Mertensia,Potentilla}_SeedSet.csv and MetaData.xlsx. Original file names and exact published MD5 are in v53 registry. This is ONE original experiment, not four independent papers. A shared plot/species/time key is NOT a same-plant key unless one source ID is verified to refer to the same plant in both visitation and reproduction.

**B. Mota et al. (2022)** flower strips / semi-natural habitat and sunflower visitor/service/productivity measurements, original data DOI 10.5061/dryad.h18931zpg, archived Zenodo record 6671199. Original files Database_Visitation_rates.csv, Database_Productivity.csv and README_file.txt; fixed original MD5 in v53 registry. Do NOT join distinct flowers, sampling surveys or field sites as if same plant identities.

Read exact registered source Zenodo record IDs via official https://zenodo.org/api/records/{id}, require every pinned original file + exact MD5 checksum, file key and bounded size. If any file cannot be accessed or verified, classify ACCESS_OR_SOURCE_STOP instead of a null effect. Do not substitute source releases.

## Before any outcome rows
For each CSV: report only original column names, file byte count, original row count, a source file MD5, and candidate identity-field names/cardinalities (plant_id, individual, plot, site, field, transect, season/year, day/date, taxon/species). Numeric visitation frequencies, fruit/seed outcome values, correlations, effects and treatment differences MUST NOT be examined, returned or committed. Source record must not be treated as independent biological replications across files.

Eligibility thresholds for next-stage direct plant-level predictive refit:
1. at least 30 original distinct individual plant IDs linking direct pollinator visitation and measured fruit/seed outcome within the SAME authored experiment (verified by ID semantics, not matching a generic numeric sequence);
2. at least 3 source-independent experimental plots/blocks that can be held out intact, with at least 5 linked distinct plants in each; all same individual repeated records stay within their plant and block;
3. original visitor measurement effort or per-plant exposure is supplied, otherwise any visit-rate comparison must be changed to clearly labeled raw count and independent exposure source;
4. reproductive outcome is subsequently measured and not itself a derived multiplication involving visitor counts; source README must support timing and grain;
5. no source leakage across holdout plots, no invented genotype/parent identity or effective pollination event; exclude incomparable taxa if plant ID semantics differ.

If both sources fail same-individual eligibility, do not invent successful plant-level replication. A separate prereadout spatial site-level test may be designed only with valid shared plot identifiers, biological model and target frozen before outcomes are opened. Never pool distinct source response scales.

## Output
A source SHA/MD5-locked audit script with bounded original Zenodo access and no numerical outcome fitting; machine-readable column/count/key provenance receipt and adversarial tests. Run on read-only GitHub Actions, keep Draft PR stacked on ANTO #62, no main merge. All legitimate source-access stops and unidentified joins are valid scientific outcomes to report.
