# anto — The biogeography of floral evolutionary potential

**Core question:** When plants colonize new regions, to what extent can lineages retain, transport, and reuse reproductive traits that were previously lost from visible expression?

**Origin:** The contrasting white and coloured flower phenotypes of the Ryukyu *Cirsium brevicaule* / *irumtiense* complex motivated this programme through [aza3](https://github.com/zuizui0223/aza3). A 2026 phylotranscriptomic study supports distinct Ryukyu lineages (doi:10.1186/s12870-026-08097-6), while POWO retains a broader synonym treatment. Neither the ancestral pigment state nor any loss/regain event is established. This remains the hypothesis-generating origin, not proof of reversal.

**Scope:** Cross-angiosperm evolutionary biogeography, led by plant lineage dispersal, colonization, extinction, introgression, and the movement of latent developmental/genetic capacities. Pollinators enter as ecological filters and selection agents, **not** the default moving protagonist.

## Distinguish five outcomes
1. Ancestral retention of an expressed trait.
2. Genuine re-expression within a lineage after documented loss of expression.
3. Immigration of an already-expressed lineage (geographic reappearance, not re-evolution).
4. Introgression or hybridization restoring alleles or regulation (genetic reacquisition with external input).
5. Environmental plasticity, measurement error, or alternative pigment pathways.

A white flower is **not** proof that anthocyanin synthesis is genetically impossible. A preserved coding gene is **not** proof of floral re-expression competence. Do not infer evolutionary reversal from present-day geographic ordering.

## State and geography
Let x=(region, lineage ancestry, observed floral phenotype, latent pathway competence, other-organ pathway use, cytotype). The model compares lineage dispersal/colonization, local extinction, expression transitions, competence degradation/repair, and introgression. Separate:
- **visible recolonization:** coloured flowers present again in a region;
- **lineage re-expression:** a formerly non-expressing lineage produces and establishes a coloured morph without importing the coloured lineage;
- **trait reacquisition through gene flow:** a local ancestry background regains expression using externally acquired genetic material.

Define accessibility as Pr(a specified trait establishes within H generations | present full state, evolutionary history, stated migration and selection regime). Never label it a measured probability without calibration.

## Research programme
- **WP0 — Origin and alternatives:** test the Ryukyu colour contrast against retention, re-expression, introgression, plasticity and taxonomy uncertainty.
- **WP1 — Theory:** derive identifiable predictions and counterexamples for visible-state-only vs hidden-competence vs migration/introgression models. Compare held-out predictions under matched information.
- **WP2 — Cross-clade history:** prospectively choose independent angiosperm clades with dated trees, geography, observed colour and molecular evidence; account for sampling, uncertainty and phylogenetic dependence.
- **WP3 — Molecular validation:** same-individual floral colour/UV, pigment chemistry, flower and other-organ expression, ancestry, cytotype and controlled growth or functional perturbation when feasible.
- **WP4 — Generalization:** only after colour validation, compare other reproductive traits (nectar guides, nectar secretion, orientation or architecture) with trait-specific definitions of latent capacity.

## Relationship to other projects
- [aza3](https://github.com/zuizui0223/aza3): molecular mechanism, standing variation, introgression and same-individual *Cirsium* sampling; anto does not claim or overwrite its empirical results.
- [izu-core](https://github.com/zuizui0223/izu-core): mechanistic models show reproductive and genetic evolutionary sequences can differ; model-specific results are not empirical proof of trait reacquisition.
- [egwe](https://github.com/zuizui0223/egwe): test whether visible states are sufficient to predict future outcomes; avoid assuming that hidden state always adds predictive value.
- [chun](https://github.com/zuizui0223/chun): visible flower-colour evolutionary memory is recurrent but non-universal; historical similarity cannot be equated with pathway competence.
- [fcp](https://github.com/zuizui0223/fcp): photographic spatial colour structure, not a direct genotype or pigment assay.

## Start here
- [Research charter and testable predictions](docs/RESEARCH_CHARTER_V1.md)
- [Model, identifiability, and falsification contract](docs/MODEL_AND_IDENTIFIABILITY_V1.md)

**Evidence status (2026-10-17):** This repository contains a new hypothesis and study design, **not** a demonstrated reversal, calibrated re-expression rate, dated accessibility window, or universal law. Existing projects remain independently governed. No fieldwork, sampling, sequencing, or ethics authorization is implied.

## Working v2 — plant dispersal × latent competence × genetic rescue

The [finite-region colonization simulator](src/anto_colonization.py) now explicitly models discrete plant generations, seed-based migration/recruitment, regional capacity, extinction shocks, inherited pigment expression capacity and pollen-mediated genetic rescue. It logs **intrinsic re-expression**, **already-pigmented seed immigration**, and **introgressive reacquisition** as *different histories*. Its observer output contains only visible pigmentation; hidden states live in a separately identified simulator-oracle field.

The [identifiability benchmark](src/anto_identifiability.py) demonstrates that three different histories can produce exactly the same visible-colour endpoint distribution. Synthetic ancestry/pigment-haplotype assays help only when genuinely diagnostic. [Replay experiments](scripts/run_synthetic_baselines.py) and [tests](tests/test_colonization.py) are included. This is model development, not Cirsium parameter estimation or evidence of real re-evolution.

### New research navigation

- [Ryukyu 2026 source audit, taxonomy conflict, BioProject and ancestry hypotheses](docs/RYUKYU_ORIGIN_EVIDENCE_V2.md)
- [Cross-angiosperm selection panel and falsifiers](docs/CROSS_CLADE_PROSPECTIVE_PANEL_V2.md)
- [Reproducible synthetic results and stringent inference limits](docs/SYNTHETIC_RESULTS_V2.md)

Commands: `python -m pip install numpy scipy`; `python -m unittest discover -s tests -v`; `python scripts/run_synthetic_baselines.py --replicates 400 --output results/synthetic_baselines_v2.json`.

**Claim ceiling:** no biological specimens were newly admitted, no genuine colour re-evolution has been inferred for the Ryukyu pair, and no universal evolutionary window has been demonstrated.

## v3 — empirical source intake and identifiability gate (2026-10-17)

**New:** The [Ryukyu/Petunia/chun source-level evidence registry](data/empirical_sources_v3.json) is validated by a [fail-closed gate](src/anto_evidence_gate.py); *no* registered system yet combines experimentally validated floral re-expression competence, dated colonization direction, dense pigment-locus ancestry and independent replication. The new source audit verifies that the 2026 Ryukyu *Cirsium* study used **37 analysed accessions** but the NCBI [PRJNA1311153](https://www.ncbi.nlm.nih.gov/bioproject/1311153) project contains **25 new young-leaf transcriptome samples**. These do not measure flower-stage re-expression.

The [exact three-history observation model](src/anto_assay_design.py) calculates why colour-only records are structurally non-identifying and how *hypothetical*, calibrated ancestry and pigment haplotype assays can help. A [four-history extension](src/anto_four_history_design.py) adds the most important anto contrast: **retained ancestral pathway reuse versus local de novo construction**. Two simple ancestry assays remain rank deficient for four histories; the third ancestral-reuse signal is hypothetical and requires independent molecular/dated validation.

- [Evidence audit and full inference boundaries](docs/EMPIRICAL_GATE_AND_ASSAY_IDENTIFIABILITY_V3.md)
- [Source-level evidence record](data/empirical_sources_v3.json)
- [One-command validation runner](scripts/run_empirical_gate_v3.py)
- [Regression tests](tests/test_empirical_v3.py)

Run \`python scripts/run_empirical_gate_v3.py --output results/empirical_design_gate_v3.json\` after installing numpy/scipy. The exact assay-rank values are **synthetic design results**, not observed molecular classification accuracy. The repository has not established actual *Cirsium* re-expression or a universal rule across flowering plants.

## v4 — verified public Ryukyu specimen crosswalk (2026-10-17)

**Actual public data, not synthetic:** [25 ENA leaf RNA-seq runs](data/public_sra/prjna1311153_run_manifest_20261017.csv) now have exact SRR experiment, SAMN BioSample, published taxon title and specimen code; among them **3 C. brevicaule**, **3 C. irumtiense**, and **6 var. takaoense**. All 6 takaoense runs independently match [aza3's three-white/three-blue-purple labels](data/public_sra/takaoense_aza3_morph_crosswalk_20261017.csv). The Ryukyu 3+3 colour contrast remains a taxon-level description, not direct individual pigment measurement. See [full evidence and data provenance audit](docs/PUBLIC_SRA_SPECIMEN_GATE_V4.md).

The [offline intake validator](src/anto_sra_intake.py) rejects invented individual floral colours, wrong tissue, swapped BioSamples, duplicate runs and unverified morph labels; the [sequencing platform audit](data/public_sra/prjna1311153_ena_platform_fingerprints_20261017.json) flags unresolved **MiniSeq (ENA) vs NovaSeq (paper)** metadata disagreement. No raw reads have been downloaded; ENA's declared paired FASTQ files total approximately 80.66 GB.

The [Ryukyu phylogenetic-polarity diagnostic](src/anto_ryukyu_polarity.py) demonstrates a substantive limitation: with symmetric transition rates, equal temporal branches and a 50:50 prior, one white and one coloured sister lineage yield **no ancestral direction** (posterior remains 0.5 regardless of elapsed time). A 0.93 Ma **lineage** split is not a pigment transition date. Next biological experiments must resolve ancestral polymorphism, introgression, pathway retention and flower-stage function before claiming re-expression.

Commands: `python -m unittest discover -s tests -v`; `python scripts/run_public_sra_intake_v4.py`; `python src/anto_ryukyu_polarity.py`.

**Status:** provenance and model-based limitations established; no demonstrated Ryukyu floral reactivation or cross-angiosperm evolutionary window.
