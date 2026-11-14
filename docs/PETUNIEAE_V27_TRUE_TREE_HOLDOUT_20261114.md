# ANTO v27 — 59 real Petunieae species, original-tree-defined clade holdouts

**2026-11-14 | Frozen scientific design before any v27 tree-based result is observed.**

## Why

ANTO v26 identified predictive gains for Petunieae species-level anthocyanidin/flavonol total HPLC abundance using F3′H/F3′5′H/FLS RNA expression. However, its 19 blocks formed by **the first two characters of a source species code** are not established taxonomic genera, monophyletic clades or independent lineages.

The next direct challenge uses the **original author dated pruned Newick tree, `phyloCCA/11genestre_dated_pruned.tre`** from official OSF zg9cu, file ID `63ae2453c8da860249df343a`, exact length **1,481 bytes**, MD5 `f270feef2909dec19a67eefe15efcb2f` and SHA-256 `95b4a688d3d71417b712b37a2b04cdc22a9431be3172d6509def4435f5fd8614`. The original independently authenticated 3 CSVs and all v26 results remain frozen and must pass exact reproduction. Only the original file download is allowed. Tree authenticity is byte verified, not inferred from a filename.

**Primary question:** Do a fixed 3-gene branch module and the 6-gene combination improve heldout predictions over a train-only mean when **entire genuine topological clades** are excluded, rather than one species or unverified code prefixes?

**Response:** Current original author *processed*, species-mean chemical anthocyanidin and flavonol totals `log1p(sum)`, not eye-visible colour, heritable pathway competence, ancestral regain, or descendant fitness.

## Before-outcome fixed analysis

1. Parse Newick, reject nonfinite lengths, unbalanced punctuation, duplicate/unmatched tips, mismatched species identifiers. No fuzzy mapping based on colour or taxon-code prefixes. Remove `BROW` only when source tree contains it. All 59 original Petunieae species must match exactly. If a source tree's labels are on another taxonomic convention, STOP and establish an independently sourced species-key mapping **before** any model result.
2. For each **fixed maximum clade size 8 and 15**, traverse the source-tree rooted topology from its root: every complete subtree containing at most the specified count is a heldout block; larger subtrees are split recursively into their direct descendants. This produces deterministic nonoverlapping, monophyletic-in-the-source-tree blocks covering **every original species exactly once**. An inevitable singleton tip is flagged separately, not called a multi-species replicate. Training has at least ten original species.
3. Run the exact predeclared six v26 ridge models (`train_mean_baseline`, `CHI_only`, `early_three`, `late_two`, `branch_three`, `joint_six`), fixed lambda=10, original train-only scaling and intercept, with exactly the same raw-source-derived response and gene definitions. Each species receives exactly one out-of-fold prediction for each model and outcome; no response is used in selecting clades.
4. Report **species-weighted** RMSE relative to **the same-clade holdout mean baseline**, equal-block signed MSE gain, number of heldout multi-species clades, fraction of tested species in multi-species blocks, number of positive blocks, and minimum gain upon dropping one existing block's score. Do not treat that last score as a newly refitted holdout or CI.
5. Perform both 8-tip and 15-tip sensitivity cuts whether the results look favourable or not. Both are topological convenience thresholds, **not biological radiations or formal lineage independence**. If phylogeny has poor coverage or little multi-species coverage, the claim is blocked. No post-outcome optimisation.
6. Keep source biological interpretation strictly bounded to one 2023 Petunieae original study and **one published topology**. Not a proper phylogenetic regression, independent phylogeny bootstrap, ancestral-state reconstruction, flower-colour phenotyping, plant-level genotype, genetic functional assay, ecological selection, or intergenerational persistence.

## Predeclared interpretation gates

- **Descriptive within-source tree transfer**: module beats the same-fold train mean on **both endpoints in both cut granularities**, with at least 3 multi-species clade folds per cut and at least half the tested species in non-singleton folds. Any failure => `NO_BROAD_TREE_TRANSFER_CLAIM`. The numerical improvements alone are not independent evolutionary replication.
- **Critical failure mode**: source original tree tips are not exactly matched. Analysis must STOP; no clade label derived by code prefix or from source colour response.
- **Sensitivity**: retain all six models and both full result matrices; do not drop a negative clade, select one convenient clade size, add a tuned lambda, or call transient colour chemistry inherited capacity.

## Reproducibility

```bash
python -m unittest discover -s tests -p test_original_petunieae_phylogeny_v27.py -v
python scripts/analyze_original_petunieae_phylogeny_v27.py \
  --expected-v26 results/petunieae_2023_original_species_pigment_prediction_v26_20261028.json \
  --output /tmp/anto_v27_real_tree_holdout.json
```

CI uses byte-verified original-author tree and processed HPLC + expression files in memory, emits only derived species-code clade membership and aggregates, and retains raw OSF source bytes **off GitHub**.

**Evidence ceiling:** This is retrospective phylogeny-structured within-Petunieae chemical prediction. Genetic access, local pollinator interaction and multigeneration fitness outcomes from ANTO × EGWE × izu-core × SDEA are not tested jointly in these sources, and original Cirsium experimental data stay in aza3.
