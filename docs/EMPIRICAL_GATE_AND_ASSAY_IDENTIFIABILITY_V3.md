# anto v3 — empirical readiness and history-identifiability results

**Research question:** When plant lineages colonize new regions, can they carry and reuse a reproductive programme no longer visible in the flower? How can a recoloured population be distinguished from a coloured immigrant, an introgression product, or a new mutation?

**Status (17 October 2026):** Primary-source metadata audit + executable exact calculations of hypothetical observation models. No real *Cirsium* regulatory regain, mutation rate, colonization direction, re-expression ability, or universal evolutionary window has been demonstrated.

## 1. Verified Ryukyu source and exact limitations

Chang et al. (2026): [original publication](https://doi.org/10.1186/s12870-026-08097-6), [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13020037/).

- **37 analysed accessions** total (33 *Cirsium*, four non-*Cirsium*); **25 new young-leaf RNA-seq accessions** plus 12 reused comparative accessions. [NCBI BioProject PRJNA1311153](https://www.ncbi.nlm.nih.gov/bioproject/1311153) currently lists **25 SRA experiments / BioSamples**, NOT 37 new ones. Registration reports approx. 57 GB sequence data. Do not download all reads just for the source audit.
- The authors report reciprocally supported *C. brevicaule* and *C. irumtiense* Ryukyu lineages separated near the Miyako Strait, white vs bluish-purple corollas, and **lineage divergence** near **0.93 Ma (95% HPD 0.71–1.33 Ma)**. This is **not** a pigment-loss event date, a re-expression date, or a directional dispersal reconstruction. The names are treated differently by older/synonymizing taxonomies.
- **Crucial source limitation:** the registered samples are *young-leaf* RNA, not same-individual flowering-stage expression data. The article's conclusion mentions anthocyanin expression/pollinator preference for *takaoense*, but accessible young-leaf methods and morphology alone do not validate a flower-specific expression or causal fitness experiment. Any such additional claim requires matching verifiable experimental materials.
- The [aza3 *takaoense* prefield registry](https://github.com/zuizui0223/aza3/blob/main/docs/TAKAOENSE_PHASEB_PREFIELD_CANDIDATE_REDUCTION_V1.md) links **3 white + 3 bluish-purple** morph labels to public young-leaf libraries, one locality per plant. Geography and elevation are confounded. Existing DFR/ANS homologous reads support plausible coding-pathway retention **only**, not demonstrated floral reactivation.

**Mechanistic positive control:** [Esfeld et al. 2018](https://www.sciencedirect.com/science/article/pii/S0960982218313484) established AN2 functional reading-frame restoration in *Petunia secreta*; [Berardi et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8364234/) linked DPL to regain of anthocyanin synthesis in *P. exserta*, including functional tests. These are closely related within one genus, not two independent cross-angiosperm colonization replicates. [chun](https://github.com/zuizui0223/chun/blob/main/docs/FLOWER_COLOUR_VARIATION_TEMPORAL_PROGRAM.md) independently establishes recurrent but **non-universal** memory of visible flower colour; visible similarity is NOT a molecular competence assay.

[Structured literature-gate registry](../data/empirical_sources_v3.json) is checked by [src/anto_evidence_gate.py](../src/anto_evidence_gate.py):

| Evidence system | Admissible tier | Integrated geography × competence claim |
|---|---|---|
| Ryukyu *Cirsium* pair | VISIBLE_LINEAGE_HISTORY_ONLY | HOLD |
| *Cirsium* var. *takaoense* white/blue-purple leaf samples | VISIBLE_LINEAGE_HISTORY_ONLY | HOLD |
| *Petunia* independently documented molecular regain | MOLECULAR_MECHANISM_ONLY | HOLD for dispersal history |
| chun visible-colour evolutionary time | VISIBLE_COLOUR_ONLY | HOLD for competence |

**Result:** **0 of 4 registered evidence systems** currently meets the full joint history/competence standard. This says nothing about whether the biological hypothesis is true; the necessary measurements are not yet integrated.

## 2. Exact rank criterion and why a new method is needed

Let H have K latent evolutionary histories and Y be an observed assay profile. If the observation emission probabilities E[y,h] = Pr(Y=y | H=h) are **known**, a mixture of unknown history frequencies π produces q = Eπ.

The history-mixture weights are identifiable from the population observation distribution q **if and only if E has full column rank K**. If rank(E) < K, distinct mixtures produce exactly the same observable distribution, however many photos one collects. Even full rank is not enough when emission parameters are unknown, assays are correlated, marker histories are wrong, or phylogenetic observations are not independent. Full rank of mixtures does not prove the causal history of a single individual.

The [three-history design](../src/anto_assay_design.py) computes exact probabilities under made-up independent ancestry assays (90% sensitivity, 95% specificity; equal class frequencies):

| Hypothetical observations | Rank (K=3) | Bayes-optimal known-emission accuracy |
|---|---:|---:|
| Colour only | 1 | 33.33% |
| Genome-wide ancestry alone | 2 | 61.67% |
| Pigment-linked haplotype ancestry alone | 2 | 61.67% |
| Both calibrated and independent | 3 | 88.58% |
| Both, but 50% independent missingness | 3 | 61.31% |
| Both uninformative (50% sensitivity and specificity) | 1 | 33.33% |
| Shared ancestral haplotypes mimic intrinsic gain | 2 | 63.08% |

The [four-history design](../src/anto_four_history_design.py) explicitly adds **retained ancestral pathway reuse vs local de novo mutation**. These are identical under genome background and simple foreign colour-haplotype observations.

| Synthetic observables | Rank (K=4) | Exact best-case accuracy |
|---|---:|---:|
| Visible colour only | 1 | 25.00% |
| Genome background + colour-locus ancestry | 3 | 66.44% |
| Above + independently calibrated ancestral-reuse evidence | 4 | 86.57% |
| The extra history proxy fails to distinguish old/new variants | 3 | 66.44% |
| Three weak proxies (65% sens/spec) | 4, low information | 42.25% |

**A hypothetical "ancestral reuse" marker is NOT presently a validated laboratory test.** Establishing variant origin needs multiple gene trees or haplotype-age distributions, orthology and ploidy checks, ancestral polymorphism/ILS and reticulation nulls, and flower-specific functional verification. The 86.57% number is an exact feature of the stipulated probability kernel, not a prediction of genotyping accuracy in real plants.

## 3. Prospective discriminators and stop rules

1. Same-plant voucher-linked corolla phenotype, pigmentation onset stage, VIS/UV optical measurements and pigment chemistry (LC-MS/HPLC). Floral white by camera is not biochemical absence.
2. Comparative corolla and nonfloral gene expression/functional pathway evidence. A preserved CHS/DFR/ANS coding read in leaf cannot classify W+.
3. Population- and cytotype-aware whole-genome ancestry + local pigment-gene genealogies. Compare with alternative ancestral sorting and historical pollen introgression; simple species-tip shared states cannot date pigment events.
4. Independent restoration tests (common garden, segregation, crosses, specific regulatory perturbation) before labeling a white plant 'competent' or a new coloured plant 're-expression'.
5. Independent dated dispersal histories with uncertainty, ideally replicate other genera selected before inspecting outcomes; rate claims must distinguish population colonization from phenotype transmission.
6. Rank/precision/held-out clade tests with missingness and weak-assay controls. If no fully validated assays, return NOT_IDENTIFIABLE; do not insert latent competence imputed from images into empirical results.

**Relation to other repositories:** aza3 is the molecular/genomic biological authority; [izu-core](https://github.com/zuizui0223/izu-core) provides model-specific path-dependence insights; [egwe](https://github.com/zuizui0223/egwe) demands state sufficiency and movement-operator validation; chun offers evolutionary visible-colour memory with nonuniversal biochemical generalization. anto uniquely focuses on **plant lineage/allele movement and historical reachability of an unexpressed reproductive phenotype**.

## 4. Reproduce without downloading raw reads

\`\`\`bash
python -m pip install numpy scipy
python -m unittest discover -s tests -v
python scripts/run_empirical_gate_v3.py --output results/empirical_design_gate_v3.json
\`\`\`

The evidence gate is an audit of **source studies**, not independent natural replicates or RNA-seq measurement. Every row declares its highest justified level and missing measurements. The exact assay-rank calculations use assumed known independent binary observation models; no empirical rate, allele origin, fitness advantage or dispersal mechanism is estimated here.
