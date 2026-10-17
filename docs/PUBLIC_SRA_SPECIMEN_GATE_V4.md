# anto v4 — public SRA specimen crosswalk and phylogenetic polarity gate
Status: 2026-10-17, publicly verified metadata, not molecular re-expression evidence.

## Actual public study intake

The [ENA PRJNA1311153 read-run report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1311153&result=read_run&fields=run_accession,experiment_accession,sample_accession,experiment_title,library_source,library_strategy,library_layout,fastq_bytes,study_accession&format=json) supplies **25** distinct SRR/SRX/SAMN registrations from Chang et al. (2026). Titles disclose species and specimen codes. These source experiment titles are not independently verified herbarium vouchers; original Table S1 remains to be inspected.

| Published taxon label | ENA run | BioSample | Title specimen code | Tissue |
|---|---|---|---|---|
| C. brevicaule | SRR35152725 | SAMN50798036 | OKN-4295 | young leaf |
| C. brevicaule | SRR35152729 | SAMN50798033 | YRN-4166 | young leaf |
| C. brevicaule | SRR35152730 | SAMN50798032 | AMM-4163 | young leaf |
| C. irumtiense | SRR35152724 | SAMN50798037 | YGN-4296 | young leaf |
| C. irumtiense | SRR35152731 | SAMN50798031 | ISG-4111 | young leaf |
| C. irumtiense | SRR35152732 | SAMN50798030 | MMY-4078 | young leaf |

Total public taxon-title labels (8): C. brevicaule 3, C. irumtiense 3, C. japonicum var. albescens 2, var. australe 3, var. fukienense 4, var. japonicum 2, var. takaoense 6, C. morii 2. The 2026 phylogeny incorporated another 12 older accessions to analyze 37 total samples. The six Ryukyu records have **no individually measured corolla pigment data** here; source-study taxon descriptions white/purple are not flower-stage measurements in those six specimens.

### Takaoense molecular-phenotype link

The separately maintained [aza3 metadata contract](https://github.com/zuizui0223/aza3/blob/main/data/contracts/takaoense_phaseb_prefield_contract_v1.json) links public codes by SRR and BioSample to colour-morph labels:

| Colour morph | Specimen | Run |
|---|---|---|
| white | WY-3560 | SRR35152717 |
| white | FB-3629 | SRR35152738 |
| white | LT-3839 | SRR35152734 |
| bluish-purple | FC-3559 | SRR35152718 |
| bluish-purple | TJ-3807 | SRR35152736 |
| bluish-purple | NH-3835 | SRR35152735 |

These are externally curated flower-morph labels, **not molecular expressions inferred from leaf RNA**. Each sample represents a different locality and morph is confounded with geography/elevation. This is not an experiment for pigmentation causation, functional restoration or adaptation.

### Sequencing-platform discrepancy (unresolved)

The [published Methods](https://doi.org/10.1186/s12870-026-08097-6) state **Illumina NovaSeq, 150-bp paired end**, while the ENA run metadata list **Illumina MiniSeq in 25/25 experiments**. This is a provenance discrepancy, **not proof of invalid sequencing or which instrument actually produced the reads**. The original ENA metadata, file byte sizes and provider-declared MD5 hashes are retained in [platform fingerprints](../data/public_sra/prjna1311153_ena_platform_fingerprints_20261017.json). ENA advertises approximately 80.66 decimal GB of paired FASTQs; we **did not download raw reads**.

## Important phylogenetic insight: colour contrast does not determine direction

Chang et al. report two Ryukyu evolutionary units, **C. brevicaule (white)** and **C. irumtiense (bluish-purple)**, whose *lineages* diverged around 0.93 Ma (95% HPD 0.71–1.33 Ma). This is not a floral pigment gain/loss date.

With the two taxa coded as observed white/pigmented tips, a symmetric equal-rate two-state CTMC on equal-time branches and a 50:50 root prior gives exactly **P(ancestral pigmented | white/pigmented tips)=0.5**: likelihoods of an ancestral white and ancestral pigmented lineage are equal by construction. Six within-species transcriptome specimens are not six independent gains or losses.

Under a deliberately hypothetical unequal-rate sensitivity at t=0.93 Myr:
- white-to-pigmented 0.1/Myr, pigmented-to-white 0.5/Myr: P(root pigmented)=0.7760;
- equal rates 0.5/Myr: P(root pigmented)=0.5000;
- white-to-pigmented 2.5/Myr, pigmented-to-white 0.5/Myr: P(root pigmented)=0.4364.

Different rate assumptions produce different ancestral conclusions from the same visible tips. These numbers are **toy mathematical consequences, not posterior inference from measured rates**. The source's C. morii samples are currently placed with Nipponocirsium, not the immediate sister to the Ryukyu Arenicola pair; do not use it as an arbitrarily close outgroup. See [the reproducible polarity diagnostic](../src/anto_ryukyu_polarity.py).

## Next biological sequence of work

1. Recover supplementary author specimen-to-herbarium-voucher crosswalk and independently validate locality metadata; no population coordinates inferred from specimen-code initials.
2. Seek existing deposited orthogroup alignments or affordable target-gene sequence subsets before requesting or downloading ~80GB of reads. Validate paralogy, coverage and cytotype.
3. Reconstruct pigment-locus genealogies against ancestry nulls including incomplete lineage sorting, historical introgression and different dated branching histories.
4. In *the same flowering individual*, measure floral pigment chemistry, flower-stage expression, nonfloral use of the pathway and inherited/genetic competence. Homologous DFR/ANS reads from leaf are not sufficient.
5. Distinguish inherited re-expression of retained pathways from new genetic reconstruction, coloured-seed immigration and pollen introgression using independently verified molecular mechanisms. Keep independence of different plant clades (e.g. Petunia) explicit.

If pigment-stage measurements, ancestry diagnostics and historical orientation cannot be independently validated, the output remains **MECHANISM_NOT_IDENTIFIABLE**, not "re-evolution confirmed".

## Reproducibility

Public sources: [ENA report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1311153&result=read_run&fields=run_accession,experiment_accession,sample_accession,experiment_title,library_source,library_strategy,library_layout,fastq_bytes,study_accession&format=json), [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/1311153), [primary article](https://doi.org/10.1186/s12870-026-08097-6), [aza3 source contract](https://github.com/zuizui0223/aza3/blob/main/data/contracts/takaoense_phaseb_prefield_contract_v1.json).

Run the exact offline source validation and polarity sensitivity using:

    python -m pip install numpy scipy
    python -m unittest discover -s tests -v
    python scripts/run_public_sra_intake_v4.py --out results/public_sra_audit_v4.json
    python src/anto_ryukyu_polarity.py

Run IDs, titles, sample codes and provider-supplied file sizes are provenance facts; no phenotype measurements, RNA abundance, gene-loss mutations, fitness effects or evolutionary accessibility probabilities have been inferred from these sources. No field sampling/permissions are conferred by this work.
