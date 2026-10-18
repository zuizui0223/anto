# anto v6 — Floral loss is not pathway loss: within-genus organ evidence and a modular competence model

**Status: source-linked empirical organ evidence + independently labelled theoretical counterexample, 2026-10-17. No Ryukyu regained colour demonstrated.**

## 1. New empirical observation from within Cirsium

The article [Park et al. (2020), *Integrated Analysis of Transcriptome and Metabolome in Cirsium japonicum Fisch ex DC*](https://doi.org/10.1021/acsomega.0c04001) provides real multi-organ qRT-PCR and metabolite contrasts. The source identifies 104,890 unigene candidates and 29 phenylpropanoid-pathway candidate orthologues. Author-reported comparisons in the text include:

| Assay/target | Numerator / denominator | Reported fold |
|---|---|---:|
| CjCHS1 transcript | flower / leaf | 500 |
| CjCHS2 transcript | flower / leaf | 112 |
| CjF3'H transcript | flower / leaf | 16.6 |
| CjCHI transcript | flower / leaf | 13 |
| CjFLS transcript | flower / leaf | 5.6 |
| CjFSII transcript | flower / leaf | 3.3 |
| CjDFR isoform A transcript (exact accession unverified) | root / leaf | 22 |
| CjDFR isoform B transcript (exact accession unverified) | root / leaf | 2.6 |
| Total flavonols | flower / root | 21 |
| Total flavones | leaf / flower | 58 |
| Total flavanols | leaf / root | 19 |

**Interpretation discipline:** the two DFR contrasts are ROOT/LEAF, **not** root/flower; do not say DFR is expressed more in roots than in flowers. CHS, FLS and FSII belong to broader flavonoid pathways; their expression cannot be substituted for an anthocyanin-specific floral pigment measurement. The metabolite classes are not measurements of flower anthocyanin. Distinct transcripts/metabolites with different denominators must not be averaged as if they were exchangeable independent replicates. These results are from *C. japonicum* sensu the 2020 source, **not** the six Ryukyu *C. brevicaule* and *C. irumtiense* transcriptomes, and show no trait regain or island colonization.

What this *does* justify is a prospective prediction: selection on secondary-metabolite functions outside the flower could maintain pieces of the biosynthetic toolkit even if the corolla becomes white. But organ-specific regulation could decay independently, so **coding retention does not guarantee floral reacquisition**.

Machine-readable [11 source-linked reported contrasts](../data/molecular/cirsium_japonicum_organ_evidence_v6.json) are validated by [organ evidence gate](../src/anto_cirsium_organ_evidence_v6.py), and source-fold denominators and original limitations are locked with tests. No raw qPCR replicate data were obtained; hence no new significance tests, pooled meta-statistic or causal effect estimate are authorized.

## 2. Biological mechanism: two hidden modules, three white types

The previous anto toy model distinguished P, W+ and W−. But this conflates structural coding integrity with flower-specific regulation. v6 instead defines two separable modules:

- **C:** structural pigment-production toolkit functionally intact in some tissue (not merely a transcript match);
- **F:** developmental/regulatory module capable of activating the pigment in the flower;
- **U:** independently observed nonfloral use that may impose maintenance selection on C.

| Latent state | Coding toolkit C | Flower activation F | Visible phenotype |
|---|---:|---:|---|
| C1F1 | 1 | 1 | Pigmented |
| C1F0 | 1 | 0 | White — coding intact, flower switch off |
| C0F1 | 0 | 1 | White — coding unavailable, activation present |
| C0F0 | 0 | 0 | White — both missing |

These are **abstract pathway states**, not directly equivalent to a single DFR transcript or single MYB gene. Multiple genes, paralogs, biochemical fluxes and environmental conditions determine real flower pigmentation. In particular, there is no evidence any Ryukyu white plant has been directly assigned to C1F0.

Why this matters for dispersal: propagules transmit alleles/genetic background that can retain C and F to different extents, so visually indistinguishable white colonists can carry radically different future phenotype possibilities. This provides a more testable meaning of "plants transport evolutionary capability" than counting white/coloured island occurrences alone.

## 3. A prespecified quantitative mechanism and negative controls

Let the **structural degradation hazard** depend on nonfloral pathway usage \(U\):

\[
\delta_C(U)=\delta_0\exp(-\beta U)
\]

This is a **testable hypothesis** that regulatory/coding reuse under nonfloral selection protects C when \(\beta>0\). It is NOT a new theorem and Park et al. did not estimate beta. If beta=0, nonfloral expression gives no structural-retention advantage. Flower-regulatory decay/recovery is independently parameterized.

A deliberately restricted special case starts with white C1F0, prohibits reactivation during T generations of exposure, then permits regulatory reactivation at hazard \(r\) during the next H generations, competing with structural decay:

\[
A(T,H\mid U)=e^{-\delta_C(U)T}\frac{r}{r+\delta_C(U)}(1-e^{-(r+\delta_C(U))H}).
\]

This is one-founder event probability, **not** population establishment in a changing metapopulation, excludes de novo structural reconstruction, gene flow, drift, changing environment and founder extinction. Independent n-founder probability is \(1-(1-A)^n\), under the same unrealistic independence conditions.

For a **synthetic illustration only**, set coding loss hazards \(\delta_C(0)=0.12\), \(\delta_C(1)=0.02\), and regulatory reactivation \(r=0.05\), with H=10:

| T since loss | U=0 reactivation chance | U=1 reactivation chance |
|---:|---:|---:|
| 0 | 0.2404 | 0.3596 |
| 5 | 0.1319 | 0.3254 |
| 10 | 0.0724 | 0.2944 |
| 20 | 0.0218 | 0.2410 |
| 40 | 0.0020 | 0.1616 |

These numbers are forced by arbitrary **hypothetical rates** and should never be called Cirsium probabilities. If \(\beta=0\), the two columns coincide; if the flower-specific switch never reactivates, neither group regains floral pigment under this restricted model even if C is fully preserved. Structural rebuilding or introgression may reopen accessibility.

The full [four-state continuous-time Markov generator](../src/anto_modular_competence_v6.py) distinguishes intact coding versus floral activation, shows identical white visible states can have different first-hit probabilities for coloured flowers, and separately acknowledges that new structural reconstruction is **not** ancestral regulatory re-expression.

## 4. What to measure to test the mechanism rather than assume it

A decisive contrast needs flower colour and leaf/root/involucre pigment pathways from the **same genetic individual** and developmental stage, plus ancestry-aware population replication. For each candidate homologous locus and tissue:

1. Locate functional coding orthologs with valid paralog/ploidy-aware genealogies; don't equate leaf RNA nondetection with gene deletion.
2. Measure flower and nonfloral transcription, pathway metabolites and calibrated visible/UV optical phenotypes.
3. Determine heritability and causal competence from controlled expression restoration/segregation or perturbation where feasible.
4. Test whether measured nonfloral use predicts structural retention, after lineage age, Ne, ancestry and introgression controls.
5. Independently resolve geographic dispersal/colonization before asking whether this hidden capacity predicts regional phenotypic recolonization better than colour-only models.

A potentially novel discovery would be **a geographic lag between loss of visible pigmentation and loss of genomic/developmental reachability**, conditional on an independently verified nonfloral function. It remains a falsifiable hypothesis.

## 5. Supplementary author material and targeted retrieval

The 2020 paper advertises [annotated unigene TXT](https://pmc.ncbi.nlm.nih.gov/articles/instance/7675961/bin/ao0c04001_si_002.txt), [unigene XLSX](https://pmc.ncbi.nlm.nih.gov/articles/instance/7675961/bin/ao0c04001_si_003.xlsx) and [PCR/figure PDF](https://pmc.ncbi.nlm.nih.gov/articles/instance/7675961/bin/ao0c04001_si_001.pdf). They are distinct from the 2026 young-leaf Ryukyu BioProject PRJNA1311153. The author files' **content has not yet been independently inspected and no bona fide Cirsium japonicum unigene sequence is incorporated in the reference panel on the basis of a supplement title alone.**

The [bounded supplementary ZIP probe](../scripts/probe_cirsium_japonicum_supplement_v6.py) can query Europe PMC's public supplemental archive with a 50 MB compressed limit and 180 MB extracted limit. It records available file types and candidate annotation strings, but never promotes text matches to functionally validated orthologs. A blocked URL means **NOT_RETRIEVED**, not "gene absent" or "data not published". The opt-in CI job records whether this route works on a GitHub runner.

## 6. Reproduce

```bash
python -m pip install numpy scipy
python -m unittest discover -s tests -v
python scripts/run_organ_history_v6.py --output results/cirsium_organ_history_v6.json
python scripts/probe_cirsium_japonicum_supplement_v6.py
# Optional public fetch, size-limited:
python scripts/probe_cirsium_japonicum_supplement_v6.py --fetch-public --output results/supplement_probe_v6.json
```

**No Ryukyu molecular comparison completed; no anthocyanin re-expression experimentally verified; no historical colonization rate or universal window estimated.** This version provides an actual Cirsium within-genus organ mechanism lead, a discriminable module model and a source-aware route to existing unigene supplements without fetching ~80 GB of RNA reads.
