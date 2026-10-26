# A lost floral phenotype does not imply loss of its recorded coding machinery: an exact 22-taxon Iochrominae source join

**2026-10-19.** This is an **original-source biological crosswalk**, not newly collected plant data, an independently estimated molecular decay rate, or a verified re-expression mechanism. The original molecular evolutionary conclusions belong to Ho & Smith (2016).

## Registered analysis unit and source provenance

The original figure 1 of Ho & Smith (2016, [BMC Evolutionary Biology](https://doi.org/10.1186/s12862-016-0675-3)) identifies **7 anthocyanin-unpigmented and 15 pigmented Iochrominae taxa**. Their article separately registered GenBank source accessions **KT898394–KT898451**. Our [NCBI retrieval CI](https://github.com/zuizui0223/anto/actions/runs/39232695589) downloaded all **58 original entries** and reconstructed 22 taxon names.

A subsequent taxonomic nomenclature screen reconciled five names that have been reassigned or synonymized since publication. In particular:
- `Acnistus arborescens` → `Iochroma arborescens` (Kew [POWO](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77186334-1)).
- `Iochroma australe` → `Eriolarynx australis` (Kew [POWO](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77186335-1)).
- `Iochroma grandiflorum` → `Trozelia grandiflora` (Kew [POWO](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77186340-1)).
- `Iochroma lehmannii` ↔ `Iochroma squamosum` (World Flora Online [taxon concept](https://www.worldfloraonline.org/taxon/wfo-0001021914)).
- `Dunalia solanacea` ↔ `Saracha nigribaccata` (both linked by [Kew](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77217353-1) under `Saracha ferruginea`; GBIF currently [treats the latter as an accepted species name](https://www.gbif.org/taxon/79MBC)).

This supports a **22-record taxon-to-source match**, not a new voucher-level historical identity adjudication. Any taxonomic disagreement in original accession names must retain the historical concept and source aliases, not silently recompute the result with newer taxon names.

## Completed biological status × record-level genotype availability

| 2016 flower phenotype (per named taxon) | Taxa | CHI deposited in KT | F3H deposited in KT | DFR deposited in KT | DFR with CDS feature | DFR present in KT without CDS feature |
|---|---:|---:|---:|---:|---:|---:|
| **Unpigmented** | **7** | 7 | 7 | 5 | 3 | 2 |
| **Pigmented** | **15** | 15 | 15 | 9 | 9 | 0 |

**What is actually independently verified:** the correspondence of source GenBank accession, original species phenotype classification, source annotation presence and exact denominators. These are *not* per-plant simultaneous genomic/genetic-functional assays. Original article sequence alignments, branch ancestral polarity and the paper's dN/dS fit have **not** been independently rerun.

**Very important counterexample to careless macroanalyses:** In these genuine source records, the only **two deposited DFR sequences lacking a CDS annotation** belong to unpigmented species (*I. confertiflorum* and *I. ellipticum*). Coding-feature absence would appear enriched in the unpigmented group if an analyst replaced unknown/missing with a gene-loss zero. But the original Ho–Smith study reported intact full coding regions after phenotype loss; unannotated CDS records are not identified pseudogenes. No test of biological gene erosion should treat these 2/5 unannotated deposits or the 8 missing KT deposit slots as knockout evidence.

**This is a source ascertainment caution, not evidence that GenBank deliberately censors genes by flower colour.** It is also not a new estimated selection coefficient. The sample is small and the discrepancy was discovered after the original data were read.

## Iochroma edule: a concrete independent study linkage, not a regain proof

Source NCBI record accessions for the *same named taxon* are:
- CHI = **KT898423**;
- F3H = **KT898401**;
- DFR = **KT898445**.

Ho & Smith's source tree marks *I. edule* pigmented within a largely historically unpigmented radiation. Larter et al. (2018; [DOI:10.1093/molbev/msy117](https://doi.org/10.1093/molbev/msy117)) discuss reappearance of pigmentation and gene-expression evolution in this clade. Thus there is a source-backed **taxon-level intersection of structural pathway transcripts and phylogenetically inferred phenotypic recurrence**.

However, the historic phylogenetic regain remains an inference in the original research. Our crosswalk is not functional validation of CHI/F3H/DFR catalysis in the same propagule, nor a proof of regulatory reactivation of a retained ancestral toolkit, nor a measured lineage persistence window.

## Relation to ANTO, AZA3 and observation-series repos

**ANTO v10/v10b** separates visual trait loss, latent functional resource availability, first induction and repeated/heritable expression. **AZA3 P04** has 25 public young-leaf transcriptome runs but no Ryukyu flower-stage genetic or biochemical competence results. The Iochrominae source audit does not fill that gap, but shows **published independent clade evidence worth following up** and exact pitfalls for inferring genetic loss from database annotations.

At a methodological level, the record-to-biological-state issue connects to **REC** (missing support/denominator), **TNOA** (retaining unresolved annotations rather than collapsing them to zeros), and **CED** (a scientific report cannot be promoted without evidence for the required distinction). This is a conceptual connection; our join does **not** constitute four independent empirical validations of one observation theory.

## What is next, and what is not

- Supported current result: a 22-taxon source trait × sequence deposit availability table with 7 original unpigmented and 15 pigmented species; within both status groups recorded structural sequence sources are common and gene-absence cannot be inferred from missing GenBank annotation.
- Not supported: a new estimate of how long evolutionary potential persists after colour loss, a new rate of molecular pseudogenization, proof of a functioning whole pigment toolkit, causal evidence that pigmentation returned in *I. edule*, or generality beyond one Iochrominae radiation.
- The next **genuine biological** comparative test needs stage-specific expression/metabolites, organ-specific function and ancestry/hybridization resolution in two or more independent clades, with controls where source coding persists but floral induction fails.

## Exact reproduction

- Original [2016 Figure-1 phenotype name-to-taxonomy source ledger](../data/iochrominae_2016_fig1_source_trait_taxon_crosswalk_v1.json).
- Real [58-accession source gene-coverage ledger](../data/iochrominae_ho_smith_original_KT_taxon_gene_crosswalk_20261019.json).
- [GenBank CDS annotation uncertainty ledger](../data/iochrominae_ho_smith_2016_original_genbank_CDS_annotation_audit_20261019.json).
- [Frozen 22 species joint matrix](../results/iochrominae_2016_trait_KT_source_join_20261019.json).
- [Validation script](../scripts/audit_iochrominae_2016_pigment_source_join_v1.py) and [tests](../tests/test_iochrominae_2016_pigment_source_join_v1.py).

Run with Python standard library:

    python -m unittest discover -s tests -p test_iochrominae_2016_pigment_source_join_v1.py -v
    python scripts/audit_iochrominae_2016_pigment_source_join_v1.py --expected results/iochrominae_2016_trait_KT_source_join_20261019.json --output /tmp/iochrominae_source_join.json

This exact join reproduces an original published dataset composition and gene-entry accounting. It deliberately does not claim discovery of a new trait reversal.
