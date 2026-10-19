# Iochrominae independent public molecular-source gate — ANTO

**19 October 2026.** This source-first verification must not be confused with a new re-evolution, genetic restoration, flower colour selection, or 12-generation persistence finding.

The original independent plant study is Ho & Smith (2016), [BMC Evolutionary Biology](https://doi.org/10.1186/s12862-016-0675-3), and [Dryad dataset 10.5061/dryad.49rs8](https://doi.org/10.5061/dryad.49rs8). The authors examined *Chi*, *F3h* and *Dfr* from flower-bud cDNA in 22 Iochrominae species plus six Solanaceae outgroups. They reported intact coding regions in unpigmented lineages as well as changes in selection constraint; those are **their** results, not newly discovered in ANTO.

## Executed access decision on 19 October 2026

The dedicated public GitHub Actions run first retrieved the dataset metadata, but the official whole-dataset /download endpoint returned **HTTP 401**. A second attempt followed the **registered version-file metadata and exact official individual /api/v2/files/{id}/download links**, and these sequence-byte endpoints also returned **HTTP 401**. These are two legitimate documented retrieval routes, not source-data negative results. No sequence bytes, codon alignments, coding integrity measurements or per-species regained-colour estimates were independently produced.

**Current machine decision: SOURCE_FILE_ACCESS_RESTRICTED_HTTP_401 / SOURCE_ONLY_METADATA_FILES_UNAVAILABLE.** The workflow now records a source-file metadata-only receipt and stops biological promotion while keeping the actual authenticated-byte verifier available for future authorized replay. Passing metadata/fixture tests must not be interpreted as having reproduced Ho–Smith sequence comparisons.

## Actual-source next work

The runnable [audit script](../scripts/audit_iochrominae_ho_smith_original_dryad_v1.py) uses the [official Dryad API](https://github.com/datadryad/dryad-app/blob/main/documentation/apis/README.md) and original dataset DOI to retrieve only the public 2016 archive (no sequence files are added to this repository). It checks exact three alignment file names, bytes and SHA256, detecting FASTA/PHYLIP/PAML/NEXUS-format first headers where possible, and emits a compact original-data receipt through dedicated GitHub Actions.

The workflow is allowed to say **source download succeeded**. It is *not* allowed to say the 2016 codon-model selection tests were independently recreated, that each species' flower pigmentation was validated, or that a complete jointly functional toolkit has been established. No original taxon-colour states are yet matched to rows.

If the source is inaccessible, report technical SOURCE_UNAVAILABLE, not biological NO_CAPACITY. The source file format will determine whether the later biological analysis is feasible; do not guess PHYLIP interleaving or treat the author-reported 22 species as 22 independent regain events.

## Why the source is independent but not yet a generality test

Iochrominae is independent of the Ryukyu *Cirsium* field system and the originally documented Petunia regulatory regain. Yet Iochrominae remains **one comparative clade** across the 2016 Ho–Smith coding and 2018 Larter et al. RNA/pigment work. The very fact that an intact gene can remain under a visibly unpigmented flower is a **necessary but insufficient** source clue for ANTO's latent competence variable.

Next admissible analysis requires a full species/source and pigment classification crosswalk, codon-frame integrity assessment, uncertainty-aware phylogenetic history, expression and biochemical data at matched biological grain. Even then, the outcome may be limited to comparative retained structural sequences and not experimentally demonstrated flower regain.

**Scientific hard limit:** do not use the frozen ANTO v10/v10b synthetic event rates as estimates for this species set.
