# anto v5 — first real Asteraceae molecular references and a fail-closed Ryukyu CDS screen

**Date:** 2026-10-17  
**Decision:** INDEPENDENT REFERENCE CDS VERIFIED; RYUKYU FUNCTIONAL RE-EXPRESSION NOT TESTED.

## What v5 actually establishes

We retrieved and validated **three real GenBank annotated coding sequences** as independent positive/reference controls, not as observed alleles in Ryukyu plants.

| Verified GenBank CDS | Biological material | Locus and actual CDS length | Relevant independent evidence |
|---|---|---|---|
| EF600682.1, positions 61–1089 | *Saussurea medusa* | DFR, 1,029 bp | DFR enzymatic activity validated in yeast; flower-vs-leaf differential expression documented in source publication |
| AY547342.1, positions 29–1099 | *Saussurea medusa* | ANS, 1,071 bp | Anthocyanidin-synthesis candidate, with cell-line pigmentation and expression analysis in published experiments |
| AJ007709.1, positions 310–1911 | *Gerbera hybrida* | GMYC1 bHLH regulator, 1,602 bp | Tissue/region-specific regulatory context for gerbera DFR; GMYC1 is NOT itself a DFR enzyme |

All three CDSs start with ATG, end with a canonical stop codon, have lengths divisible by 3 and contain no internal stop codon under the standard genetic code. The real sequences are in [data/molecular/asteraceae_reference_cds_v5.fasta](../data/molecular/asteraceae_reference_cds_v5.fasta), with [accessions, exact annotated CDS coordinates and literature sources](../data/molecular/asteraceae_reference_provenance_v5.json). This validates **reference extraction**, not homology or biological function in *Cirsium*.

**Literature anchors:**
- SmDFR EF600682: Li et al., DOI [10.1007/s11033-011-1061-2](https://doi.org/10.1007/s11033-011-1061-2); PMID 21701830.
- SmANS AY547342: Cheng et al., DOI [10.1007/s11240-007-9211-x](https://doi.org/10.1007/s11240-007-9211-x).
- Gerbera GMYC1 AJ007709: Elomaa et al., DOI [10.1046/j.1365-313x.1998.00273.x](https://doi.org/10.1046/j.1365-313x.1998.00273.x); PMID 9807831.

The gerbera finding is a particularly useful warning: organ- and position-specific anthocyanin regulatory change is plausible in composite flower heads. However GMYC1 is a bHLH regulator and must not be substituted for an anthocyanin-specific SG6 MYB. The DFR/ANS/GMYC1 target set is **a deliberate small pilot panel**, not an exhaustive pigmentation pathway inventory.

## Status of actual Ryukyu sequence analysis

The 2026 *Cirsium* [source publication](https://doi.org/10.1186/s12870-026-08097-6) points to **raw young-leaf reads** from 25 newly sequenced accessions in [PRJNA1311153](https://www.ncbi.nlm.nih.gov/bioproject/1311153). We queried ENA's project-scoped secondary-analysis and sequence registries: both returned **zero entries** on 2026-10-17. This does not demonstrate that processed alignments are absent at other repositories or in author files; it means **we have not independently retrieved target-specific, validated ortholog sequences for any of the six focal Ryukyu samples**.

The frozen [source availability audit](../data/molecular/source_availability_audit_v5.json) records exactly which registries were queried and the evidence boundary.

In the [18-row register](../data/molecular/ryukyu_gene_screen_register_v5.csv), six already-identified biological specimens are crossed with DFR, ANS and GMYC1. All 18 have status **NOT_TESTED**, and all six samples have no expression/anthocyanin confirmation in their flower tissue. This is intentional, not a negative biological gene test.

## First test once Cirsium candidate coding sequences become accessible

The [offline candidate coding screen](../src/anto_molecular_screen_v5.py) accepts optional externally generated candidate CDS FASTA with exact provenance-constrained headers `>SRR35152725|DFR|contig_identifier`. It verifies:

1. The sample and target were in the predeclared run×locus register — otherwise reject.
2. A/C/G/T/N sequence integrity, reading frame, start/terminal codon, internal stops and ambiguity.
3. Multiple contig/paralog candidates are **retained and flagged**, never collapsed into an invented single orthologue.
4. Apparent premature stops are coded **POTENTIAL CANDIDATE REQUIRING GENOMIC VALIDATION**, never confirmed gene loss.
5. The output never claims independently validated orthology, floral pathway function, phenotype causality, introgression or local re-expression.

A negative leaf-RNA result cannot exclude a functional gene, because transcription is tissue- and time-specific; 80.66 GB of RNA reads should not be downloaded indiscriminately just to fill a speculative matrix. A full, intact coding RNA still does not establish intact regulatory circuits in florets or the ability to undergo re-expression.

## Minimum study design for a decisive Ryukyu history

| Competing account | Molecular evidence required | Unresolved even if a DFR/ANS-like leaf transcript is seen |
|---|---|---|
| Coloured ancestor, southern colour retained | Rooted species/allele histories, outgroups, uncertainty and flower-stage pigment assays | Branch on which white expression disappeared |
| White ancestor, coloured southern lineage re-expresses | Independently polarized ancestral loss, retained regulatory/coding route, floral restoration experiment | Whether alternative pathway or novel allele was involved |
| Introgression of colour-associated alleles | Dense nuclear ancestry and locus genealogies with coalescent/ILS null, ploidy-aware callset | Age/source of donated haplotype |
| De novo reconstruction | Derived coding/regulatory variants with ancestral age exclusion, validation of flower-stage causality | Rate and likelihood of reconstructed state |
| Environmental plasticity | Controlled shared-growth conditions, stage-matched pigment/expression | Genetically heritable change |

The **first cost-effective molecular step** is to obtain the published orthogroup/coding alignments from the original author if available, or reconstruct a rigorously bounded subset of gene candidates after an explicit resource budget is approved. Target-read subsets cannot reliably be identified from raw compressed FASTQ by family name alone; reference-guided alignment or assembly is needed, with controls for paralogy and missingness.

The source report's lineage split at approximately 0.93 Ma is not the date of any anthocyanin reacquisition. The same two-tip white/coloured contrast remains non-polarizing under the symmetric ancestral model in v4.

## Run, tests and claim ceiling

```bash
python -m pip install numpy scipy
python -m unittest discover -s tests -v
python scripts/run_molecular_preflight_v5.py --output results/ryukyu_molecular_preflight_v5.json
```

This v5 analysis is **source-backed biological reference preparation and method validation**, not observed comparative DFR/ANS coding divergence in *Cirsium*. There is no basis yet to assert that white *C. brevicaule* populations have intact floral anthocyanin competence, that *C. irumtiense* reacquired colour, or that a universal re-evolutionary time window exists across angiosperms. Anecdotal leaf homologous read detections remain candidate retention signals only until independent QC.
