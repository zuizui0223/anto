# Independent Iochrominae sequence source — 58 published INSDC accessions

**Source:** Ho & Smith (2016), [DOI 10.1186/s12862-016-0675-3](https://doi.org/10.1186/s12862-016-0675-3), original GenBank accession range **KT898394–KT898451** (58 records inclusive). This is an independent plant clade relative to Cirsium. The dataset is still **one historical study**, not 58 biological regain origins.

## Executed original-source evidence — 19 October 2026

ENA Browser API retrieval of all **58/58** original GenBank records succeeded in [source Actions run 39229422231](https://github.com/zuizui0223/anto/actions/runs/39229422231) (7/7 offline gate tests). The source artifact preserves each accession's nucleotide sequence checksum, EMBL record checksum, submitted taxon and CDS annotation; no raw sequence is redistributed.

Actual source accession coverage: **22 F3h**, **22 Chi**, **14 Dfr-range deposited accession records** across **22 submitted taxa**. Only **12/14** Dfr-range records carry a CDS/product annotation; KT898443 (*Iochroma confertiflorum*) and KT898446 (*Iochroma ellipticum*) have no such annotation. **14/22 taxa** have new accessions in all three gene blocks, whereas **8/22** have only the new F3h/Chi accession pair in this particular published block. These counts are accessions/annotations, not gene presence/absence or loss events.

Other actual metadata: **56/58** records have CDS annotations, **23/58** records contain ambiguous nucleotide codes (52 individual ambiguous positions); new deposited sequence lengths are Chi 705 nt, F3h 1089 nt, Dfr-range 1143–1170 nt. The original article separately describes coding alignments of Chi 717, F3h 1113 and Dfr 1179 bp and analysed pre-existing accessions as appropriate. Differences between the newly deposited *individual record lengths* and the author's *analysis-alignment lengths* must be investigated rather than interpreted as deletions or pseudogenization.

[Source-locked 22-taxon matrix](../data/external/iochrominae_2016_genbank_58_source_coverage_v1.json) and [replay evaluator](../scripts/check_iochrominae_ena_coverage_v1.py) derive these counts from the actual SHA-verified 58-record Actions artifact. The evaluator explicitly fails if the original file universe or gene annotations differ.

**Crucial publication boundary:** Ho & Smith (2016) already concluded structural genes were conserved and selection was relaxed at some loci following floral pigment loss. The present independent computation confirms *source-record availability and coverage*, not those original biological inferential results, functional enzyme competence, historic regain, or evolutionary persistence. The 8 taxa without *new deposited* Dfr accessions must **not** be coded as Dfr-null.

## Why the second source path matters

The original Dryad 2016 archive and registered individual-file downloads both returned HTTP 401 in public GitHub Actions (the [explicit stop receipt](../data/external/iochrominae_2016_dryad_access_gate_20261019.json)). The original authors also deposited nucleotide records at GenBank, accessible via the European Nucleotide Archive's INSDC mirror.

The [ENA Browser API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html) has a documented accession-specific EMBL flatfile endpoint. The [source audit code](../scripts/audit_iochrominae_ena_genbank_58_v1.py) retrieves each of the **58 exact accession IDs**, compares submitted accession IDs, records sequence size, checksums, organism annotations, coding features and gene names, and emits only a source receipt. No raw sequences are rehosted.

## What a successful source replay would establish

Actual retrieval confirms that the archived nucleotide **records exist and are readable** at specified versioned IDs. A CDS annotation or amino-acid translation is metadata/sequence evidence, NOT validated gene function, expression in white flowers, independently counted pigment losses, a restored allele, or viable progeny.

Ho & Smith already reported intact coding regions at Chi/F3h/Dfr and changing constraint in unpigmented Iochrominae. Repeating their author conclusions without an independent codon test is not a new biological discovery. A new cross-clade test must also match flower phenotype and genealogical origin at compatible taxa and assess the functionally complete joint pigment pathway, not an isolated nuclear gene.

## Next strictly delimited biological test

If all accession records can be retrieved, the **first genuinely empirical independent step** is a specimen/taxon-verified source crosswalk for the original Chi/F3h/Dfr genes, identify which lineages were independently inferred by the original authors to be pigment-lost, reconstruct comparable lineage-level coding integrity and note that functional competence remains unknown. Use full published source status; do not infer ancestry from species names or individual sequence presence.

For a truly new answer to the ANTO evolutionary-afterlife question, compare more than one independent non-Cirsium lineage with same-individual flower-stage pigment enzyme and regulatory evidence, dated trait-loss/re-expression state, historical gene flow, and heritable offspring reproduction, with clade-level holdout. This source audit **does not accomplish that final test**.
