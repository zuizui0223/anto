# ANTO v48 — when training-source availability is matched, deamii is no longer the exceptional negative transfer case

**2026-11-19. Post-outcome source-population coverage challenge.** This refits the same authentic Goulet-Scott et al. (2024) original Phlox common-garden field dataset (Nature Communications 15:4661; DOI 10.1038/s41467-024-49075-6), already authenticated in ANTO v46, after the v46/v47 deamii prediction failures were known. Source original author commit d792fe58cce36721394adfd97e9f412396578e86; source CSV 175377-byte Git blob 4c780b0abe644ac69dae86a0450eb598c160626b. Nothing here is a new original study, new independent garden, independent plant radiation, source cohort expansion, or flower colour / pollinator / molecular-functional experiment.

## Motivating scientific issue

v46 fruit model held out an entire source population, fitting the remaining **19 source populations** including different numbers from the same taxon (deamii **2**, amoena **7**, pilosa **8**) in three common gardens. Source-population-equal fruit MSE improved in aggregate (+9.39% G, +9.43% E) yet worsened in deamii. v47 then matched only **evaluation-group count**, comparing deamii's three heldout population losses with C(8,3)=56 amoena and C(9,3)=84 pilosa subset aggregates (without any new refit); for G and E, deamii was worse than 140/140 reference subsets. This **cannot** separate taxon-specific biology from asymmetric train-population representation.

## New original-source re-fits, complete design

All original 963 ramet rows / 284 genotypes / 20 source populations / 3 common gardens were source hash-verified again. Every panel is three ORIGINAL source populations chosen from a single taxon. For each panel, three source-population heldout fits use **two other same-taxon populations for training**, all three original gardens, and no other taxa. Thus all train/test panels have exact source-population count (2 training, 1 evaluation). Original genotype clones never straddle train/test.

Exactly **56+84+1 = 141 triple panels**, **423 source-population-heldout refits**, each with a training-only garden intercept B0 plus G pooled-two-leaf-PC slopes or E garden-specific-two-leaf-PC slopes, fixed ridge alpha=10, training-only PC scaling, and the same original log1p fruit outcome. With one focal taxon per fit, the original v46 T model equals G and v46 TE equals E in parameter layout, so do not pretend to test four distinct nested candidate models. Every original source-population observation is stored only as a derived MSE/fit receipt; original source CSV is not redistributed.

## Primary original fruit endpoint — new source-derived numeric conclusion

| Predictor | deamii equal-pop MSE gain (3 source populations) | Number of 140 other original three-pop panels equal or worse | Number of 140 references with negative MSE gain | Median reference gain |
|---|---:|---:|---:|---:|
| G: shared two leaf slopes across the original three gardens | **−4.53%** | **34/140** | 55/140 | +2.64% |
| E: garden-specific leaf slopes | **−3.90%** | **44/140** | 60/140 | +1.96% |

The negative deamii sign remains, but the earlier v47 "unusually negative 0/140" conditional comparison **DOES NOT SURVIVE** this change in comparable training information. Other taxa often have a three-original-population panel that predicts as badly or worse when trained using only two within-taxon original source populations. The source-specific mean positive transfer from v46 also cannot be assumed under this reduced-data target domain.

Original deamii population-level fruit results in v48 are:
- G: original heldout populations 1719, 1720, 1721 have prediction MSE gain −5.822%, −8.332%, −0.551%, respectively.
- E: −6.167%, −6.517%, +0.001% (last approximately zero; not evidence of an actual beneficial biological effect).

This demonstrates remaining **model and population heterogeneity**, not a uniform biological sign. Per-test training RAMET counts also differ: 114, 156, 174, even though the number of source training populations is exactly two in each fold.

## Secondary floral-count outcome (NOT floral morphology)

The original flower-COUNT target (log1p flowers_2019) shows deamii G **−1.56%** with **29/140** other panels as poorly predicted, and E **−2.07%**, **39/140**. It is the same original 2024 plants, not independent biological replication, nor phenotype pigmentation, corolla size, pollinator attraction or stable flower allele selection. The original 13 negative signed biomass observations (minimum −233) still block v46's prereadout log1p biomass endpoint; no alternative scale or posthoc wins are manufactured.

## Critical model-family and causal limitations

This test better matches **within-taxon training-source population support**, but it also changes the model estimation domain relative to v46/v47:
- Original v46 G/E fit slopes using other TAXA as well, while v48 G/E is fitted **within one taxon only**;
- thus source-population count, other-taxon pooling, training row count, genotype spread and PC distribution all change;
- the 140 reference panels heavily overlap, and many panel folds recycle the same original populations, so **34/140 and 44/140 are NOT p-values, 140 independent population experiments or a valid randomization null**;
- the three original gardens are unchanged, with zero independent new common gardens; the total source count remains 20;
- the 2024 authors had already found that different trait axes predict between- and within-taxon fitness differences. ANTO v48 does not rediscover that selection biology or isolate a genetic/fecundity mechanism.

Therefore the careful result is: **the unusual deamii ranking was sensitive to the training-population support/composition and fitted scope**, whereas the limited-leaf-information transport gate remains unestablished. No independent mechanism is identified. It would be incorrect to say the number of training populations *alone* explains deamii's trait–fitness relationship, or that deamii truly shares all fitness slopes with amoena/pilosa.

## Theory connection and next falsification

This supports a methodological principle for ANTO × SDEA × EGWE × izu-core: before interpreting apparent out-of-domain predictive failure as biological mechanism divergence, audit (i) heldout observational grain, (ii) training-population support, (iii) model-family domain, and (iv) independently linked stage-specific biological endpoints. In this one 2024 Phlox source, the apparent diagnostic tail changed from v47 0/140 to v48 34/140 or 44/140 when training representation was altered. This is a demonstration of sensitivity, NOT a universal evolutionary theory proof.

The next controlled analysis would retain v46's multi-taxon nine-intercept G/T/E/TE model architecture and change *only* focal-taxon training population representation while keeping other-taxon training as comparable as possible, ideally in an independently sampled dataset. To establish biology rather than a failure mode of prediction, use independent deamii origins and additional gardens alongside directly measured floral phenotype/pollinator/mating/seed outcomes.

## Reproduction and status

- Freeze-before-v48-readout posthoc protocol: docs/ANTO_V48_TWO_TRAIN_POP_MATCHED_PROTOCOL_20261119.md
- Immutable 2024 author-source and population/model scope: data/external/anto_v48_original2024_train2_registry_20261119.json
- SHA-verified original CSV input, train-only fit and 423 refits: scripts/analyze_Phlox_2024_train2_v48.py
- Source and synthetic adverse tests: tests/test_Phlox_2024_train2_v48.py
- **Complete derived receipt: all 141 source triple panels and each of their 423 source-population-heldout B0/G/E MSE rows for fruit and flower:** results/anto_v48_phlox2024_train2_matched_source_20261119.json
- Original-source Actions: .github/workflows/anto-v48-two-train-pop-matched.yml

Dedicated original-source repeat run [GitHub Actions #44431670331](https://github.com/zuizui0223/anto/actions/runs/44431670331) **SUCCESS** with exact frozen receipt replay; general theory run on same head also **SUCCESS**. Draft stacked on v47, **no main merge** without independent scientific review.
