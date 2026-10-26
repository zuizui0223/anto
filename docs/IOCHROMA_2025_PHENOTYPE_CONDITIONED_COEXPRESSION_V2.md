# Iochroma 2025 original flower RNA — pooled versus within-colour transcript correlation

**26 October 2026. New independent source-data *reanalysis* of original author-provided files, not a new sequencing experiment, natural trait regain or functionally verified regulatory network.** The original dataset was published by Wheeler, Larter & Smith (2025), *G3*, DOI 10.1093/g3journal/jkaf230; source https://osf.io/j5m8f/.

## Main original-data result

The source-verified 24 original flower-bud RNA libraries belong to **one backcross** of blue and red flower phenotypes, **12 individuals per flower-colour class**. Previously the complete, MD5-verified matrix of **19,184 transcript entries** showed that **16 transcript entries** strictly separate all blue from red individuals in the author's selected sample, including an authors-nominated MYB12-like row. This selection was not external validation or confirmation of 16 independent causative genes.

The new audit tests **the other 15 strict blue-above-red transcript entries**, with the original nominated MYB12-like transcript `TRINITY_DN3000_c0_g1_i10.p1` as reference.

| Descriptive result, same 15 other transcripts | Entries |
|---|---:|
| Positive Pearson r on all 24 individuals combined | **15 / 15** |
| Negative Pearson r on all 24 individuals combined | 0 / 15 |
| Positive Pearson r after within-colour mean-centering | 6 / 15 |
| Negative Pearson r after within-colour mean-centering | **9 / 15** |
| Undefined because a centered vector has zero variance | 0 / 15 |

**Result:** nine of the 15 correlations change from positive to negative when the blue/red group mean of *each* transcript is subtracted. This is a demonstrable case in which the pooled transcript coexpression picture is sensitive to a strong observed phenotype-group contrast, not evidence that MYB12 actively represses nine genes.

**Important meaning of "within-colour":** the centered Pearson coefficient is calculated **across all 24 group-centered residuals**. It is *not* the assertion that a negative correlation was separately observed in **both** the blue group and the red group. No transcript-wise independent intervention, gene-regulatory direction, or external predictive generality follows.

## Exactly reproducible calculation

For original individual i in blue (B) or red (R) group g and the original untransformed TPM value of transcript j, construct

[
X^*_{ij}=X_{ij}-overline X_{j,g}.
]

For MYB12 transcript m and each other pre-selected transcript j, calculate two ordinary Pearson coefficients:

[
r_{\text{pooled},j}=\mathrm{cor}(X_{im},X_{ij}),\qquad
r_{\text{colour-centred},j}=\mathrm{cor}(X^*_{im},X^*_{ij}).
]

Each original transcript is centered **separately** in the 12 B and 12 R individuals, then pooled for the residual correlation. **TPM values are not log transformed**, and no arbitrary pseudocount is added. An exactly constant centered transcript returns undefined (null), not r=0. The coefficient has no hypothesis-test p-value attached: all 15 original entries were selected post hoc from the *same* complete-separation dataset, and are not 15 independent genomic loci or independent plant lineages.

### All 15 original-source coefficients

| Original transcript entry | Pooled Pearson r | Flower-colour-centered r |
|---|---:|---:|
| `TRINITY_DN11681_c0_g3_i5.p1` | +0.632 | +0.189 |
| `TRINITY_DN12748_c0_g2_i2.p1` | +0.581 | +0.083 |
| `TRINITY_DN1372_c0_g1_i30.p1` | +0.568 | -0.283 |
| `TRINITY_DN2718_c1_g1_i1.p1` | +0.727 | +0.222 |
| `TRINITY_DN2718_c1_g2_i1.p1` | +0.789 | +0.400 |
| `TRINITY_DN3164_c0_g1_i15.p1` | +0.523 | -0.305 |
| `TRINITY_DN3164_c0_g1_i6.p1` | +0.549 | -0.321 |
| `TRINITY_DN5710_c0_g1_i14.p1` | +0.551 | -0.132 |
| `TRINITY_DN5810_c0_g1_i18.p1` | +0.629 | -0.254 |
| `TRINITY_DN6389_c0_g1_i4.p2` | +0.912 | +0.794 |
| `TRINITY_DN6574_c0_g1_i12.p1` | +0.534 | -0.331 |
| `TRINITY_DN6994_c1_g1_i1.p1` | +0.475 | -0.357 |
| `TRINITY_DN753_c0_g1_i2.p1` | +0.681 | -0.015 |
| `TRINITY_DN901_c0_g1_i2.p1` | +0.902 | +0.752 |
| `TRINITY_DN9157_c0_g1_i2.p1` | +0.384 | -0.015 |

The above values are derived from **live OSF original source data**, with original provider byte lengths and MD5 hashes independently verified before computation. The stable detailed receipt is [results/iochroma_2025_color_conditioned_coexpression_v2_20261026.json](../results/iochroma_2025_color_conditioned_coexpression_v2_20261026.json). The original large data matrix is **not republished** into this repository.

## Why the finding matters, and what it does not show

This is a genuine additional *diagnostic result* beyond merely reproducing the authors' MYB12 expression/blue-red association: the 15 other strictly discriminating transcripts do not all display positive *within-group-adjusted* coexpression. Strong phenotypic group differences can generate a pooled positive association despite opposite within-group residual covariation.

It **does not** show gene repression or an ancestral pigment regulatory memory state. The original study is **blue versus red within a backcross**, not loss of colour followed by regained ancestral pigment expression. Neither genotype-level complete function of a floral ancestral toolkit, long-term lineage persistence, heredity of a reversible state nor actual seed/offspring fitness was measured by this analysis.

The independently published 2019 Iochrominae comparison (Larter et al., DOI 10.1002/dvdy.82) distinguished white morphs from fixed-white species, but those original Dryad files were not reanalysed by ANTO (the repository's prior retrieval was blocked HTTP401). The current blue-red RNA measurements must not be misrepresented as a replication of that white-flower paper.

For ANTO v9–v18, the lesson is about **measurement sufficiency**: pooled transcript signatures of a visible state are not necessarily reliable functional-state proxies. More original lineage-level evidence is required before testing historic competence or inferring parameters for v17/v18's environment-memory and multigeneration branching models. No mathematical demographic rate in those models has been fitted from these 24 individuals.

## Locked evidence pipeline

- Source script and group-stratified correlation method: [scripts/audit_iochroma_2025_osf_myb12_backcross.py](../scripts/audit_iochroma_2025_osf_myb12_backcross.py).
- Frozen all-15 coefficient receipt with 24-individual original-source and biological-scope flags: [results/iochroma_2025_color_conditioned_coexpression_v2_20261026.json](../results/iochroma_2025_color_conditioned_coexpression_v2_20261026.json).
- Source MD5 verification and exact frozen-receipt comparison in [GitHub CI](../.github/workflows/iochroma-2025-backcross-source-audit.yml).
- Positive/negative synthetic probes for color confounding, undefined centered variance and source drift: [tests/test_iochroma_2025_osf_myb12_backcross.py](../tests/test_iochroma_2025_osf_myb12_backcross.py).

Reproduce the externally sourced diagnostic with:

```bash
python scripts/audit_iochroma_2025_osf_myb12_backcross.py \
  --output /tmp/iochroma_original_tpm_audit.json \
  --expected results/iochroma_2025_24_individual_myb12_16_transcript_source_audit_20261026.json \
  --conditioned-expected results/iochroma_2025_color_conditioned_coexpression_v2_20261026.json
```

**Claim ceiling:** source-derived descriptive coexpression sensitivity in one original blue/red Iochroma backcross. It is not a new functional-regain demonstration, a genetic causal effect, a cross-clade meta-analysis, or empirical confirmation of ANTO's branching parameters. Original Cirsium specimens/fitness/genetic data remain in [aza3](https://github.com/zuizui0223/aza3).
