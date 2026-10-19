# Iochrominae original GenBank sequence re-audit: retained coding annotation is not recovered floral competence

**2026-10-19 JST.** Outcome: 58/58 original Ho–Smith (2016) GenBank accessions independently retrieved, exact 22-taxon × gene-source crosswalk materialized, deposited CDS annotations assessed. **This is a direct original-source audit, not a discovery of functional reacquisition or a fitted evolutionary-clock study.**

### Actual original data now retrieved, independently of Dryad

- Source study: Ho & Smith (2016), *Molecular evolution of anthocyanin pigmentation genes following losses of flower color*, [DOI:10.1186/s12862-016-0675-3](https://doi.org/10.1186/s12862-016-0675-3).
- Original data registration in paper: **GenBank KT898394–KT898451**.
- Authoritative [GitHub Actions run 39232695589](https://github.com/zuizui0223/anto/actions/runs/39232695589) actually fetched and identified all **58 records** through NCBI nuccore; local workflow used Biopython and exact accession-identity gates.
- [58-accession/22-taxon source crosswalk](../data/iochrominae_ho_smith_original_KT_taxon_gene_crosswalk_20261019.json), [source script](../scripts/audit_iochrominae_original_genbank_kt_accessions_v1.py).
- Compact real original-record [coding annotation audit](../data/iochrominae_ho_smith_2016_original_genbank_CDS_annotation_audit_20261019.json).

| Original GenBank group | Exact deposited records | With CDS annotation | Raw record lengths (nt) |
|---|---:|---:|---|
| CHI | 22 | 22 | 705 in every record |
| F3H | 22 | 22 | 1089 in every record |
| DFR | 14 | 12 | 1146 ×11, 1145 ×1, 1143 ×1, 1170 ×1 |
| **Total** | **58** | **56** | Mixed record lengths not direct presence/absence phenotypes |

The 2016 article discusses coding alignments of CHI 717 nt, F3H 1113 nt and DFR 1179 nt. Those quantities need not be identical to raw GenBank record lengths. They refer to different sequence representations and/or processing; we did not audit the author alignments here. **The differences do not establish biologically meaningful deletions or frame-disrupting mutations.**

Specific deposited DFR entries without a CDS feature are **KT898443 / Iochroma confertiflorum** (1145 nt) and **KT898446 / I. ellipticum** (1170 nt). The lack of the GenBank feature is an **annotation uncertainty**, not proof of an inactive enzyme or evolutionary loss. Before inferring an open reading frame, an independent source-aware reading-frame check is necessary.

Eight other Iochrominae have a CHI and F3H accession within this original registered KT interval but no DFR accession *within the same interval*. This does **not** demonstrate absent DFR loci: published DFR sequences also exist in antecedent research such as Smith, Wang & Rausher (2013) ([DOI](https://doi.org/10.1093/molbev/mss255)), including functionally examined blue/red Iochroma lineages. A source-registration gap and a biological gene deletion are different estimands.

### Mechanistic synthesis with an actually independently documented candidate

**Iochroma edule** has all three original registered sequence records in the NCBI source: CHI `KT898423`, F3H `KT898401` and DFR `KT898445`. In a separate published comparative study (Larter et al. 2018, [DOI:10.1093/molbev/msy117](https://doi.org/10.1093/molbev/msy117)), this taxon was interpreted as a possible return to floral anthocyanin synthesis within an ancestrally unpigmented subgroup.

This gives a **cross-study species-level alignment of structural sequence availability and a phylogenetically inferred phenotype recurrence**, not same-individual genetic co-retention, a validated ancestral gene-restoration mutation, a controlled induction experiment or an independently demonstrated persistence window. The original two studies are not two independent evolutionary regain origins.

The paper's original empirical conclusion is that coding sequences remained intact through repeated anthocyanin losses, despite elevated dN/dS after loss at CHI and F3H, with most ratios still under purifying selection. These are **the original authors' findings**, not new dN/dS tests in this repository.

### What this does and does not say about ANTO's hypothesis

The cross-repository question is **whether phenotypic absence can coexist with retained functional resources and, if so, under what conditions those resources enable repeated re-expression before becoming inaccessible**.

This source audit establishes a *provenance-backed opportunity* to independently recheck part of the structural sequence story in an independent plant clade. It has **not** established the functional three-module state required by ANTO v10. Annotated structural coding does not equal expression in a flower, pathway activity, gene regulation, regenerated pigmentation, heritable persistence or mature seed fitness.

The decision gate is to resolve the biologically meaningful difference among:
1. structural sequence remains annotated while floral expression is absent;
2. floral regulatory or biochemical activity is retained/recoverable;
3. colour is induced transiently;
4. true heritable floral expression persists and affects fitness;
5. a colored phenotype returns only through external introgression/migration rather than intrinsic reactivation.

Only the first level has partial original-source evidence here. Keep all downstream levels **UNMEASURED** in ANTO, as in the [aza3 P04 source audit](https://github.com/zuizui0223/aza3/blob/main/docs/P04_AFTERLIFE_SOURCE_READINESS_AND_MEASUREMENT_GAP_20261019.md).

### Next paper requirement

Do **not** submit a new Evolution Letters/NEE universal biological claim from 58 deposited sequences, since the original paper already established structural persistence and the genetic mechanism of any actual return remains unresolved here.

A novel cross-clade result would need a prequalified comparison in multiple independent flowering-plant clades of *gene function + floral regulatory competence + ancestry-conditioned regain + repeated ecological success*, with model predictions tested on whole held-out clades/populations. At present, the source matrix is eligible for such future work but cannot alone test a rate of closure.
