# ANTO v50 — Same-source coefficient sharing versus ridge coordinate sensitivity in Phlox

**2026-11-19. Explicitly POSTHOC extension of ANTO v46–v49.** The 7 scenarios below were frozen before the v50 first readout in docs/ANTO_V50_SHARING_PENALTY_PROTOCOL_20261119.md; the scientific question was chosen after reading the older Phlox outcomes. These results are not independent study replication, a randomization test, estimates of natural flower-colour selection or a proof of a general macroevolutionary law.

## Original field source, biological grain and actual analysis

Original Goulet-Scott et al. (2024) Nature Communications, DOI 10.1038/s41467-024-49075-6; original author source PhloxHopkins/PhloxFieldAdaptation, pinned commit d792fe58cce36721394adfd97e9f412396578e86; original author CSV exact Git blob SHA1 4c780b0abe644ac69dae86a0450eb598c160626b. Source values: **963 ramets of 284 original genotypes, 20 original source populations across three Phlox taxa and three original common gardens**. No additional source studies, gardens, pollen observations, genetic transmission measurements or functional floral gene assays added in this analysis.

v50 reuses the **same original v49 141 three-population panels, 423 complete heldout-population folds, 2 same-taxon + all remaining nonfocal training-population policy, 9 original garden×taxon B0 mean cells**, and both original response endpoints. Each original genotype clone is held out with its complete source population across all gardens. Predictors are the two original author-assayed GREENHOUSE LEAF PC axes; not colour, corolla shape, hybrid fitness, pollinator recognition or genetic toolkit competence. All response values and training-only centering are original v46; no select-and-refit on heldout outcome.

After the v49 frozen source result (Git blob 94947d03e85766a3c5c1f82fde2bb076f2ca42e8), v50 recomputes original source panel risk under four original models:
- G: one pooled two-leaf-PC slope across taxa and gardens.
- T: two slopes unique to each original taxon.
- E: two slopes unique to each original garden.
- TE: two slopes unique to each garden×taxon group.

Seven prespecified original-data model coordinate/scaling settings: standardized PC with λ=0,1,10,100; raw centered PC with λ=0,10; and raw PC with analytically equivalent diagonal penalty λ×SD_PC² to reproduce the original standardized λ=10. B0 retains the same original nine training-only group intercepts in EVERY setting.

## Independent numerical readout — main fruit endpoint

Values are relative **source-population-equal MSE gain compared with the same original B0**, deamii one 3-population source panel. **Negative means predictive performance gets worse** after adding leaf PC information; it is not negative organismal selection or an estimated genotype-fitness coefficient.

| Source-locked coordinate / penalty | G: pooled | T: by taxon | E: by garden | TE: garden×taxon |
|---|---:|---:|---:|---:|
| Original standardized PC, λ=0 (OLS) | **−9.98%** | −5.35% | **−9.72%** | −6.28% |
| Original standardized PC, λ=1 | **−9.97%** | −5.18% | **−9.67%** | −5.70% |
| Original v49 standardized PC, λ=10 | **−9.83%** | −3.99% | **−9.25%** | −2.80% |
| Original standardized PC, λ=100 | **−8.58%** | −0.80% | **−6.40%** | −0.05% |
| Original raw centered PC, λ=0 (OLS) | **−9.98%** | −5.35% | **−9.72%** | −6.28% |
| Raw PC, scalar λ=10 | **−9.88%** | −3.75% | **−9.41%** | −2.19% |
| Raw PC, mathematically equivalent penalty | **−9.83%** | −3.99% | **−9.25%** | −2.80% |

In **all** four standardized λ=0/1/10/100 scenarios, the G and E original deamii source-panel MSE gain stays negative, and **0 of 140** overlapping same-taxon amoena/pilosa reference triplet panels score at or below deamii under either G or E. By contrast the taxon-specific T and TE models have many comparison panels equal or more negative (e.g., at standardized λ=10: T 34/140, TE 47/140; at λ=0: T 37/140, TE 53/140). Therefore the strong G/E source-specific negative-transfer ranking does NOT vanish when scalar ridge shrinkage is removed.

The regularization affects magnitude and some predictions considerably: deamii TE fruit MSE deterioration shrinks from −6.28% at λ=0 to −0.05% at λ=100; under fixed global PC scaling that shift cannot be interpreted as trait biological rescue, since strong shrinkage suppresses predictor contributions.

## Strong mathematical and source-replication controls

For each **423 heldout source-population fold × four original models × two original endpoints**, demand:
1. Standardized λ=10 MSE reproduces unchanged v49 original model within absolute 1e−7;
2. Raw PC OLS λ=0 MSE reproduces standardized PC OLS λ=0 (parameterization invariance);
3. Raw PC with **PC-specific coefficient penalty λ×global SD_PC²** reproduces standardized λ=10 (same physical ridge prior in different units).

Independent original source calculations report maximum absolute MSE deviations from frozen v49 **about 5.0e−10**, raw-vs-standardized OLS **about 2.7e−15**, and ridge-equivalent raw-vs-standardized **about 2.7e−15**. These are numerical controls, not biological replication. Original v46 all-taxon fruit prediction gate remains FAILED. The entire source-derived receipt must be regenerated by exact-head original-author SHA-verified Python GitHub Actions before labeling v50 verified.

The *counterfactual* raw centered PCs with a scalar λ=10 **does change** model shrinkage relative to standardized λ=10. It alters magnitudes but leaves deamii G/E negative. That is a legitimate fixed-coordinate sensitivity, not evidence that scale never matters.

## Secondary flower-number endpoint

Same original 2024 plants, **count of flowers**, not flower colour, colour saturation, nectar guide or natural pollination. Standardized λ=0 deamii G **−12.75%**, T −2.17%, E **−11.78%**, TE −4.52%; standardized λ=10 G **−12.54%**, T −1.13%, E **−11.19%**, TE −1.00%. At λ=100, T and TE both become weakly positive (+0.69%, +0.63%), while G/E remain negative (−10.93%, −7.62%). They do not constitute an independent natural evolutionary experiment. The original v46 predeclared log1p biomass still fails source-domain admission due to **13 negative signed source biomass entries, minimum −233**: no clipping, shifted log, deletion or posthoc positive selection.

## Inference boundary and integration

The strongest defensible result is **data- and model-specific negative transfer from shared leaf–fitness predictive coefficients under asymmetric phylogenetic/taxonomic support**. It persists under OLS with exactly matched response/grain/baselines, so a *positive scalar ridge penalty is not required to generate it*. It is not uniquely isolated as a biological taxon-selection mechanism: parameter sharing is also a model constraint, original taxa occupy potentially different predictor distributions, and the original 2024 study had unequal population counts. These are 140 overlapping reference panels from ONE original study, not 140 independent studies, calibration data for a macroevolutionary prior, p-values or uncertainty intervals.

In ANTO × SDEA × izu-core × EGWE terms, this is an **identifiability/transferability caution**: predictive effects inferred across ecological domains can fail when sharing coefficients despite true additional within-domain information. It does not assert a universal law of ancestral floral trait regain, or supply same-genotype gene→pollinator→compatible-pollen→seed→multigeneration evidence.

## Files and review gate

- Frozen *before v50 readout* design: docs/ANTO_V50_SHARING_PENALTY_PROTOCOL_20261119.md
- Source SHA, model and nonclaim manifest: data/external/anto_v50_ridge_sharing_original_source_registry_20261119.json
- Standard-library exact original source fitting and 2×2 coordinate control: scripts/analyze_Phlox_original_scale_penalty_v50.py
- Adversarial parent/source/coordinate checks: tests/test_Phlox_original_scale_penalty_v50.py
- Full original derived 141-panel × seven-scenario × four-model × two-endpoint gains and exact source fold-loss digest: results/anto_v50_sharing_penalty_source_sensitivity_20261119.json
- Read-only original author GitHub Actions: .github/workflows/anto-v50-original-scale-penalty.yml

**Keep Draft / stacked on PR #59.** Main unchanged. No merge or universal biological claim before exact result source CI and independent scientific review. Future substantive test requires separate original Phlox/other clade gardens and biological pollen/seed/offspring responses, not another split of this 2024 source alone.
