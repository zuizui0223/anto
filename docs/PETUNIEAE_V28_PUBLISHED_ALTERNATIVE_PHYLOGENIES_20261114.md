# ANTO v28 — published alternative-phylogeny challenge (design frozen before outcomes)

**2026-11-14 | Prereadout scientific design.** The previous ANTO v27 59-species result is a within-source current-chemical prediction, not a natural trait regain.

## Question

Does the pre-existing branch expression module (F3primeH/F35H/FLS) retain lower heldout RMSE versus the train-only mean when its training folds come from genuinely different **published phylogenomic trees** rather than the 2023 original dated tree?

Pezzi, Wheeler, Freitas & Smith (2024), Molecular Phylogenetics and Evolution 198:108136, doi:10.1016/j.ympev.2024.108136, report gene-tree discordance, ILS and hybridization in Petunia–Calibrachoa–Fabiana. The 2024 study reuses related Wheeler RNA resources; this is a **tree-method sensitivity**, NOT independent experimental HPLC/RNA validation.

## Source identities fixed before prediction

1. **2023 original chemical/RNA source:** Wheeler et al. (2023), doi:10.1098/rspb.2023.0275; OSF zg9cu. Inherit existing authenticated expression, HPLC and original Newick source file ids, bytes, MD5/SHA256, and the complete frozen v26/v27 aggregate readouts. Source rows not redistributed.
2. **2024 published tree source:** Authors' original GitHub repository https://github.com/pedrohpezzi/Petunia-Calibrachoa-Fabiana_TreeDiscordance, fixed exact commit **0b33d5e992c39e191393bc0b86671c5c6a378e1e**.
   - ASTRAL species-assigned tree: Data/ASTRAL/SpeciesAssigned/ASTRAL_SpeciesAssigned_pp_scored.treefile, **1850 bytes**, original Git blob SHA-1 **04b4333e7a2e505e87dd321df1f4be88562ea4c4**; expected 43 species-level tips.
   - IQ-TREE individual consensus: Data/IQTree/SuperMatrix_PFCOut.partitions_iqtree.txt.contree, **4271 bytes**, original Git blob SHA-1 **dbbe0ff0f352c6fc163f6c1a3117612bf0982984**; expected 117 individual-level tips.
   - Original sample mapping: Data/HyDe/speciesmap.txt, **1269 bytes**, original Git blob SHA-1 **30b670968ab1033cee6cb3eb5f74904529d9794b**; expected 117 original sample-ID mappings.
   - Fetch only three pinned small external author files; verify body byte lengths and computed Git-blob hashes. Any download, identity or mapping mismatch is a **technical STOP**, never a biological zero or a filled-in dummy tree.
3. The published 2024 species tree contains 37 focal Petunia/Calibrachoa/Fabiana plus 6 reference outgroups. Among the 2023 chemistry's 59 ingroup codes, exact overlap must be **31 focal species**, plus **six 2024 source-designated outgroups** (37 combined). All species and OUT status are derived from original author mapping; no response-driven species selection.

## Before-outcome fixed contrasts

- ASTRAL uses the actual published 2024 species tree.
- IQ-TREE's published tree has individuals, not species. Check every individual code against the original author mapping and corresponding ASTRAL species ID. Choose the **lexicographically first individual per species** as primary tree extraction, with **middle and last** representatives as two mandatory secondary sensitivities; do not imply these are independent published species-tree estimates.
- Root each tree at the pendant edge of published OUT code **HUTE**, then prune to the same **31** focal overlapping original measured species. Alternative OUT **PLNY** is a prespecified secondary root, and retaining all **37** overlapping codes is an independently labeled secondary panel. All trees are re-rooted with the same chosen anchor before pruning; branch lengths are not model inputs.
- Primary trees: (i) 2023 original tree, (ii) 2024 ASTRAL species tree, (iii) 2024 IQ-TREE first-individual reduced tree. For each, apply exactly the predeclared v27 maximal monophyletic subtree size caps **8** and **15**. Keep folds strictly species-disjoint with >=10 training species.
- Same **six** v26 ridge models, pre-fixed lambda=10, no test-feature/outcome tuning, both endpoints separately: log1p(original anthocyanidin total) and log1p(original flavonol total). All targets reflect **present** species-mean chemical abundance, not visible flower hue, ancestral capacity, selective advantage or heritable regain.
- **Primary descriptive criterion**: for all **3 trees × 2 caps × 2 endpoints**, branch-three RMSE improvement against the *same-fold train-only mean* must be positive. Each cut must have >=3 multi-species clades containing >=50% of the **31** focal taxa. Otherwise verdict **NOT_UNIFORMLY_PORTABLE_ACROSS_PUBLISHED_TREE_DEFINITIONS**. Preserve failures and all six modules; no post-result changes.
- All three published source families' outcomes are correlated because the chemical observations are identical. Do not pool 12 contrasts, treat tree alternatives as independent pigment replicates, or calculate lineage-wide evolutionary fitness probabilities.
- The two additional IQ-tree representative rules, alternate root PLNY and the 37-species outgroup-inclusive panel must be reported independently even if their results are negative. These are robustness checks, **not** separate HPLC samples or preregistered confirmations.

## Mandatory integrity and falsification gates

Check original 2023 v26 and v27 readouts exactly; 2024 source SHA and lengths; original author 43 species tree and 117 IQTREE tips/mapping IDs; unique IDs; exact 31/37 overlap; 6 source-classified OUT codes; no taxon-name guessing; real original 2023 tips and published external tree origin; correct HUTE/PLNY rooting and exclusive test folds; negative controls for missing tips, duplicate samples, genotype-to-species incorrect mapping, missing samples and perturbed test outcomes.

All original source trees, source processed HPLC/RNA observations and unredacted per-individual genotype source data are kept off ANTO. Store only reproducible fold-defined aggregate metrics, source digests and rejection states.

**Claim ceiling:** a retrospective published-phylogeny sensitivity of species-level current flower pigment chemistry. No direct causal gene function, phenotype-state regain, 2024 independent chemical study, future floral persistence, or joint ANTO×EGWE×izu-core natural mechanism is identified. The real Cirsium organismal study stays in aza3.

**Pre-readout state: no new v28 model outcomes yet examined.**


## Authenticated v28 results and important limitations (2026-11-14)

First original-source GitHub Actions success: https://github.com/zuizui0223/anto/actions/runs/43583999875 . All three Wheeler 2023 original processed CSVs and original Newick are hash-verified, the **previous v26 AND v27 results are reconstructed exactly**, and the Pezzi 2024 ASTRAL tree, IQ-TREE and sample mapping satisfy all author-file byte-length and Git-blob SHA-1 source pins. These are not arbitrary in-code trees. The declared **31 focal species** and **37 including author-OUT species** were recovered exactly. No rejected source data or taxon aliases were silently replaced.

### Main fixed 31 focal-species panel: branch-three versus same-fold training mean

| Published tree definition | Max heldout tips | Heldout blocks / multi-species | Anthocyanidin RMSE reduction | Flavonol RMSE reduction |
| --- | ---: | ---: | ---: | ---: |
| Wheeler 2023 original 11-gene tree | 8 | 7 / 6 | **12.6208%** | **9.8195%** |
| Pezzi 2024 ASTRAL species-assigned | 8 | 7 / 6 | **12.6208%** | **9.8195%** |
| Pezzi 2024 IQ-TREE first sample/species | 8 | 7 / 6 | **12.5520%** | **10.3590%** |
| All three tree definitions | 15 | 3 / 3 | **11.1933%** | **10.4002%** |

All **12** predeclared directional comparisons improve, so the narrowly worded **descriptive numerical gate passes**. Alternative original OUT rooting with PLNY, last/middle original IQTREE samples, and all-37-taxon secondary panels also retain positive branch-three RMSE improvement; they are **methodological sensitivities on the exact same original pigment data**, not independent biological confirmations. Keep all other five predictors and all negative cases in the frozen result JSON.

### Crucial post-result identity audit: numerical agreement is largely because the partitions are IDENTICAL

Exact fold membership was compared **after** the original runs, not used to choose cutoffs or models.

- For the **15-tip limit**, all three trees produce **the very same three test blocks** for the 31 shared focal species. Hence the three identical numerical RMSE values are **literally one reused holdout design**, not three independent topology challenges.
- For the **8-tip limit**, Wheeler 2023 and Pezzi 2024 ASTRAL have **identical seven test blocks**. Only the IQTREE first-sample tree changes some assignments; its predictively modest differences are not enough to establish that uncertain shallow topology is harmless.
- The alternate published tree samples and roots are not posterior phylogenetic uncertainty draws, and the overlapping focal source is smaller than the original v27 59-species panel. The apparent reduction from v27's 16–18% to v28's 10–13% mixes **changed species sampling and changed heldout partitions**, not a simple decline attributable to phylogeny.

**Therefore a scientifically stronger claim that branch-three pigment predictions are robust across genuinely discordant test clades is NOT justified by this v28 alone.** The positive empirical conclusion is narrower: within these 31 shared focal species and the realized coarse holdout partitions, current chemical totals remain partly predictable from branch-pathway gene expression.

This audit mandates an explicitly **separate, post-v28 exploratory follow-up** at finer tree-defined partitions that actually differ, without silently promoting it into the v28 predeclared test.

Full original-source-derived frozen numerical receipt (all modules, folds, reference digests, secondary negatives): ../results/petunieae_alternative_published_phylogenies_v28_20261114.json . Dedicated workflow: ../.github/workflows/anto-v28-published-phylogenies.yml .
