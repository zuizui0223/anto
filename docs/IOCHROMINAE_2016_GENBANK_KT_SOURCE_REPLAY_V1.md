# Iochrominae 2016 original-author GenBank source — alternative to restricted Dryad

**19 October 2026. Source-first original-record verification, not an estimate of regaining colour.**

Ho & Smith (2016; https://doi.org/10.1186/s12862-016-0675-3) explicitly state that their deposited primary nucleotide records are **GenBank KT898394–KT898451**. This inclusive accession interval contains 58 entries. It is different from the Dryad six-file archive, whose original ZIP and individually registered file-download endpoints returned HTTP 401 in [ANTO's previous exact-source audit](IOCHROMINAE_2016_ORIGINAL_SEQUENCE_SOURCE_GATE_V1.md).

**The 58 accessions are an author-registered collection of deposited sequence records**, not 58 species, 58 independent gene losses, or 58 recoveries. The paper's main study used 22 Iochrominae ingroup species (7 without floral anthocyanins, 15 with) and 6 outgroup taxa. The other outgroup reference sequences can predate the study and need not fall in the KT range. The paper reports coding regions Chi (717 bp), F3h (1113 bp) and Dfr (1179 bp), with low divergence and no obvious truncations in its ingroup sequences. Those are the **original authors' findings**, not new results until verified against accessions.

## Execution

The [source audit](../scripts/audit_iochrominae_original_genbank_kt_accessions_v1.py) retrieves precisely the registered KT interval directly from the NCBI nuccore GenBank API in bounded 10-accession batches, checks all accession identities, and stores only per-record taxon, gene annotations, sequence length and nucleotide SHA256. Raw sequences remain with NCBI; ANTO does not rehost them.

A successful source audit provides an independent accession registry and *provenance*. It does **not** verify enzymatic function from the CDS field, certify historic colour-loss polarity, reproduce the published branch dN/dS likelihood models, or estimate the waiting time before floral colour can no longer reappear. Source access errors must remain acquisition failures; they are not biological null results.

## Next inference gate

Before contrasting genes among pigmented and unpigmented lineages, source-authorized links are needed between recorded accession, taxon, original article colour state, biological sample and original rooted gene tree. The published article's supplemental Table S1 includes GenBank taxon crosswalks; this must be independently checked before recording an ancestor-state or pigment-loss branch. Selection and phenotype claims also need appropriate outgroup and paralog controls.

The new generality question remains **whether conserved pathways after colour loss support actual re-expression and persistence**, which the author's no-frameshift observation alone cannot establish.
