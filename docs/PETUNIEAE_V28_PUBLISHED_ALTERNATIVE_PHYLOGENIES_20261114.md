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
