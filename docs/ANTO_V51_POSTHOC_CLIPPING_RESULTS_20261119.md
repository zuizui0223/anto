# ANTO v51 — Original Ulex prediction saturation audit: FAR2 zero incremental effect is an artifact of [0,1] clipping

**2026-11-19; POST-OUTCOME diagnostic**, designed only AFTER the original v51 first readout revealed exactly equal FAR2 M0 and M1 MSE at all six original sampling sites. It must not be represented as preregistered or independent field evidence.

## Exact original author source

Same Castellanos et al. 2023 original Ulex Zenodo record 7761289, five original MD5 files, 240 source tagged plants / 212 admissible fruit-trait records / 202 genotyped flower+fruit individuals in SIX sampling sites of ONE genetically continuous regional population. Same old SDEA v15 elevation-only M0 and floral-area/flower-weight quadratic M1, ridge alpha=5, site-fold train membership and exact original source outcome. The parent v51 first-source derived JSON Git blob SHA1 is **8f1d5435d2d14c2b2c999671ff9801cffab01a71**. Nothing is refitted adaptively using outcomes, no original author plant observations are redistributed.

## Posthoc raw-score-before-clipping counts (202 original source individuals each regime)

| Frozen training-site support | M0 raw linear score <0 | M0 >1 | M1 raw linear score <0 | M1 >1 | Same clipped M0/M1 prediction at the same 0 or 1 bound |
|---|---:|---:|---:|---:|---:|
| FULL5, all five other sites | **0/202** | 0/202 | **0/202** | 0/202 | **0/202** |
| NEAR2, two closest mean-elevation sites | **59/202** | 0/202 | **27/202** | 0/202 | **27/202** |
| FAR2, two farthest mean-elevation sites | **111/202** | **91/202** | **111/202** | **91/202** | **202/202** |

In every one of six FAR2 heldout sites, **every individual** had a raw M0 and M1 prediction outside the probability unit interval, and both model predictions clipped to the exact SAME boundary for every individual: at the low-elevation original sites Cheste, Montserrat and Simat all source predictions clipped to 0; at high-elevation AresM, Chiva and Sot.Chera all clipped to 1. Site counts corresponding to those two groups are **111** and **91**. This exactly explains why original FAR2 M0 and M1 squared-error scores were identical across all six source sites, **despite the pre-clipping predictions being numerically different**.

The result is a source-verified **statistical extrapolation / bound-saturation** failure, not a natural zero floral-trait selection coefficient or evidence that pollinator-mediated function was absent in distant environments. The original Ulex field reproductive endpoint is a twig fruit fraction, not measured physiological genotype fitness or descendant viability.

For NEAR2, the boundary counts differ between M0 and M1 (59→27 raw negative scores). The original relative MSE gain M0→M1 is partly dependent on how the model gets back inside [0,1] as well as on same-season floral trait information. Relative within-regime information must not be compared without its absolute MSE and baseline counterpart. FULL5 has no pre-clipping predictions beyond the [0,1] interval and has much smaller absolute heldout MSE.

**v51 primary conclusion unchanged:** FULL5 original M1 site-equal MSE 0.007713521166, NEAR2 M1 MSE 0.031636807505 (**4.10× higher**), FAR2 0.449083092642 (**58.22× higher**); NEAR2 beats FULL5 in only one of six original sampling sites. There are NO six genetically independent field populations. Original future trait evolution, parent-offspring genetic effects, individual pollinator visits and causal selection mechanisms remain unidentified.

## Reproducibility

- Original prereadout v51 primary protocol and original result: docs/ANTO_V51_ULEX_SITE_TRANSFER_PROTOCOL_20261119.md and results/anto_v51_ulex_six_site_elevation_analogue_20261119.json
- Explicit post-outcome extra diagnostics design: docs/ANTO_V51_POSTHOC_CLIPPING_PROTOCOL_20261119.md
- Source-identity and no-new-independent-study lock: data/external/anto_v51_posthoc_clipping_registry_20261119.json
- Original 202-plant raw-versus-clipped prediction audit: scripts/audit_original_Ulex_v51_clipping_posthoc.py
- Adversarial source and source-parent tests: tests/test_Ulex_v51_clipping_posthoc.py
- Source-only frozen derived counts, no individual plant values: results/anto_v51_posthoc_Ulex_clipping_audit_20261119.json
- Both primary and posthoc source checks in .github/workflows/anto-v51-ulex-site-analogue.yml

The source-first original-author GitHub Actions run #44529280701 finished its original MD5, source-parity, 60-pair primary and posthoc clipping steps successfully. A fresh exact-head source-replay check after committing this result file is the separate final gate. Maintain PR #61 as Draft and main unchanged.
