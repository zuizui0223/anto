# ANTO v34 — nearest-available-relative matched-set positivity audit

**2026-11-15 | Analysis frozen prior to any v34 original-tree feasibility results.** Exploratory sequel to ANTO v33, not a biological replication or a significance test.

## Question and why it matters

v33 preserved target identity and training-species count, but controls almost always retained a target's own heldout clade peers. Their nearest available phylogenetic relatives were therefore much closer than in true clade exclusions. This confounds *which species are retained* with *how far away the nearest retained relative is*. Before fitting another model, verify whether a comparator can exist that fixes (a) original target, (b) training count, and (c) nearest-training patristic gap.

**Predeclared estimand:** *structural existence of alternative training sets*, not a model error or a trait association. If no such alternative exists under exact matching, report **NOT IDENTIFIABLE UNDER SAME TARGET / COUNT / GAP**. Do not simulate synthetic substitutes or define an empty control as zero effect.

## Frozen sources, comparisons, tolerances

- Original Wheeler et al. 2023 1,481-byte OSF original dated pruned 11-gene Newick `phyloCCA/11genestre_dated_pruned.tre` (file ID `63ae2453c8da860249df343a`, SHA256 `95b4a688d3d71417b712b37a2b04cdc22a9431be3172d6509def4435f5fd8614`) **only**, fetched from authenticated official public OSF metadata and verified against the exact v26 file length/MD5/SHA256. Do not require new pigment/source datasets to answer structural support.
- Original 59-species ingroup codes, original `BROW` external outgroup pruned without changing numeric source-to-source paths; original rooted maximal monophyletic holdout caps **8** and **15**, as v33. Verify all 59 test species occur once per cap, and source clade counts and sizes exactly match prior v33 receipts (14 and 6). Assert v33 original-source full-record digest `a61d5c7c7e0552a45356f0f966a2dbdbecd91e47a3980a06bc06d078ed312042`; do not change v33 results.
- Distances are complete-tree original Newick patristic sums, calculated independently in v31 with Decimal and compared to v30's original parser; source units are **not verified calendar years**. No gene, HPLC or other phenotype values may be read for selecting controls.
- Compare four **fixed tolerances** for nearest-training gap matching: **absolute ±0.00001 source branch units**, and **relative ±5%, ±10%, ±20%** of each species' original clade-exclusion nearest gap. The latter are exploration of loss of exact matching, NOT additional independent biological studies.
- For each target i in a real clade C of m tips (actual train size T=59−m), define g0=minimum patristic distance from i to the 59−m original training species. A candidate training set S must omit i, have exactly T distinct source measured species, and obey min_j∈S d(i,j) in [g0−tol,g0+tol]. Reconstruct this by partitioning the 58 candidates into **closer** (distance below lower), **boundary** (within both bounds), and **farther** (above upper). All closer must be excluded and ≥1 boundary retained. Number of admissible sets is exactly `C(boundary+farther,T)−C(farther,T)` (negative or impossible counts fail, not coerced to one). Because original training is valid, the count must be ≥1. All admissible sets that are truly different from original count as distinct counterfactual options; `N=1` is **not** a replication.
- Report all 59 targets × 2 caps × 4 tolerance definitions structurally, but no original raw tree or individual molecular observations in repository or artifact. Aggregate by cap/tolerance: count with **only one possible set**, >1, impossible, median forced exclusions, minimum/maximum alternatives, fraction with >1; count how many *original clades* have at least one supported target, and original training count/nearest gap summaries. Preserve clade-specific support signatures (source tip IDs only, no genes/pigment).

## Interpretation guardrails

The strict claim `same target, same number of retained species, and exactly the same nearest usable evolutionary relative distance` may be mathematically inseparable from the same original full-clade exclusion in an ultrametric tree. The audit will expose this **positivity failure** and show how broad a gap-band relaxation is necessary for any alternate training composition. A relaxed interval does **not** isolate the causal effect of phylogenetic distance; additional training-set molecular and HPLC compositions change.

No biological model is re-fitted in v34. The v33 D>0 exploratory finding stays true under the v33 comparator, and v32 LOO primary gate remains false. Do NOT claim v34's structural support implies an incremental gene benefit. If matched controls are identifiable under a relaxed tolerance, a subsequent separately labeled study may run source-authenticated P and P+G comparisons conditional on that tolerance.

## Adversarial tests and workflow

Pure-source structural functions: reject duplicate/missing tips, invalid or negative edge lengths, wrong original species denominators, impossible original folds, wrong root/tip set, source digest drift, alternate-set combinatorics with no boundary, equal-distance polytomies and exactly one alternative. An independently calculated small ultrametric rooted toy proves that exact clade-exclusion nearest-gap matching has one admissible training set. One non-ultrametric toy proves genuinely distinct conditional sets can be found. No original phenotype outcomes are read. CI executes unit tests and original authenticated Newick source audit, writes derived-only `results/petunieae_v34_gap_matching_support_20261115.json`, then an exact frozen read-only replay. PR is Draft, stacked on v33, not merged.

**At design freeze: v34's original-tree structural support counts were not viewed.**
