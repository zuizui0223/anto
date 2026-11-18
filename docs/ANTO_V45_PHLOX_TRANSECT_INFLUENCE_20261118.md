# ANTO v45 — Phlox same-source transect influence, component signs and reaggregation

2026-11-18. POSTHOC DESCRIPTIVE ANALYSIS of immutable ANTO v44 original-C++ forward-model derived receipt. Not a prospective model-selection test, new simulation, fitted migration/selection effect, paper-reported full-cohort reproduction or an independent wild population replication.

## Precise question

In v44, removing phenotype-level epistasis and halving the pollen dispersal kernel each produced a negative aggregate adult log-likelihood change but a small positive aggregate mother-offspring component change. Does this reflect opposite signs WITHIN original geographic transects, and does the aggregate direction survive omitting any one of the five source transects?

The correct biological observation unit is not 5 independent studies: these are five spatial transects from ONE original PhloxSelection example-input source, with 32 populations, 11,379 adults, 371 mother-family records and 2,496 offspring. The parent 2014 paper stated 12,024 / 382 / 2,581; the missing 645 / 11 / 85 observations have NOT been reconciled. All scenario changes are around ONE unfitted 13-parameter START vector.

## Source-preserving new calculation

Read parent Git blob SHA-1 92d940d60051f57bad023e290e0dd5fbff5b9822. Freeze all 7 source-case statuses. For each of FOUR admitted nonbaseline scenarios, collect source-computed per-transect changes in adult and conditional offspring-family log-likelihood; independently verify each scenario's all-five source aggregate is the sum of its exact five components. For every source transect, describe whether adult and offspring delta signs are opposite. Finally, OMIT each existing source transect in turn and subtract its original contribution from the full source total. These omissions are arithmetic influence diagnostics WITHOUT refitting original source parameters or rerunning its C++ equilibrium. They are not new cross-validation folds, uncertainty intervals or leave-region-out model skill.

## New numerical findings

| Fixed intervention | Whole-sample adult ΔlogL | Whole-sample offspring ΔlogL | Opposite adult/offspring sign within transects | Which omitted transect reverses aggregate adult / offspring / total sign |
| --- | ---: | ---: | ---: | --- |
| Fitness neutral | -9655.426 | -994.213 | 0 / 5 | none / none / none |
| Epistasis zero | -261.291 | +4.474 | 2 / 5 (3 and 5) | none / 2 / none |
| Dispersal half | -432.213 | +8.674 | 0 / 5 | 1 / 2 / 1 |
| Dispersal double | -437.574 | -89.802 | 0 / 5 | none / none / none |

**Important qualitative change from a naive v44 reading:**

1. Half dispersal shows **zero local transects where adult and offspring likelihood signs conflict**: transects 1, 4, 5 are adult-negative AND offspring-negative; transects 2, 3 are adult-positive AND offspring-positive. The full-sample adult-negative/offspring-positive contrast is an aggregation of different geographic contributions, not an observed within-transect sign reversal.
2. Removing transect 1 from that original fixed-start half-dispersal comparison changes remaining 4-transect adult ΔlogL to **+27.569182** and total to **+45.253809**. Thus original full-sample negative adult/total direction is not robust to single-transect exclusion.
3. Removing transect 2 changes the half-dispersal offspring ΔlogL from +8.674 to **−17.064269**. Removing transect 3 leaves a marginally positive offspring ΔlogL of **+0.162803** (not zero).
4. Epistasis-zero has within-transect adult-vs-offspring sign opposition only in transects **3 and 5**. Removing transect 2 reverses the global offspring delta to **−7.452721**; adult and total signs remain negative under all five single deletions.
5. Fitness-neutral and doubled-dispersal original fixed-start total likelihood losses are negative in each original transect and after every single-transect deletion. This only demonstrates descriptive source START sensitivity, not fitness estimation or optimized null rejection.

## Why this matters to ANTO × EGWE × izu-core × SDEA

The correct claim is **stage-specific predictive/likelihood information can differ, and geography changes the apparent aggregate**. A coarse system-level mismatch is NOT automatically a local developmental/reproductive tradeoff. If spatial blocks weight different mechanisms differently, adult-only observations cannot identify cross-generation transmission, even when the whole original-model likelihood difference has a stable sign. This is an identification boundary, not a general natural evolutionary law. Do not pool these five transects with independent Ulex, Sinapis or Ipomopsis original studies, and do not claim that flower-colour regained function was measured.

## Reproduction and review

Source-only standard-library script: scripts/analyze_Phlox_transect_influence_v45.py.
Adversarial tests: tests/test_Phlox_transect_influence_v45.py.
Frozen derived-only analysis: results/anto_v45_phlox_transect_influence_20261118.json.
GitHub Actions: .github/workflows/anto-v45-phlox-transect-influence.yml.

Run: python -m unittest discover -s tests -p test_Phlox_transect_influence_v45.py -v
Then: python scripts/analyze_Phlox_transect_influence_v45.py --check

All positive and negative v44 observations retained; v45 does NOT rerun or alter v44 source author code, parent receipt, original study or the main branch. Keep stacked Draft until reproducibility and scientific review. A subsequent genuine inferential step needs predeclared parameter refitting AND the full original paper cohort, with appropriate geographical heldouts; neither is supplied here.
