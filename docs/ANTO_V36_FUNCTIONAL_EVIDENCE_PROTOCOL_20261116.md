# ANTO v36 — tissue-specific pathway capacity, regulatory inheritance, and experimental repression

**16 November 2026 — retrospective source-authenticated published-evidence extension of v35.** Not new genotype/chemical measurements or formal meta-analysis. Claims about original authors' results are separated from our source-audit calculations.

## Driving question

Does a visibly white flower imply the anthocyanin pathway is destroyed in the whole plant? Distinguish (i) flower appearance and chemistry, (ii) transcription and organ-dependent biochemical production, and (iii) experimental perturbation of an inhibitory regulator. Function of a regulator that SUPPRESSES anthocyanin is not experimental restoration of pigment in the native white-flowered species.

## Original evidence sources (one original paper = one program)

**Duncan & Rausher (2020)**, DOI 10.1371/journal.pone.0231263, https://pmc.ncbi.nlm.nih.gov/articles/PMC7153891/ ; original Tables 1 and 3, Figs 2–3:
- Original HPLC: white Ipomoea lacunosa n=3 plants, corolla anthocyanidins undetected but stems anthocyanidins detected; purple n=3, both organs positive. This demonstrates floral-specific loss despite whole-organism production. It does not demonstrate a previously white flower naturally regained pigment.
- Separate field qPCR: white n=4, purple n=4 flower-throats. Author-reported CHS ~100-fold, DFR >100-fold, MYB1 ~100-fold lower in white. Values are approximate descriptions, not re-estimated from raw qPCR source.
- Table3 genotyped crosses: purple MYB genotype T_ and purple phenotype 23; tt and white 30; 0 discordants across four families, total 53. Families 1–3 = 41 B2 and family 4 = 12 F3. Dfr alleles did NOT co-segregate with flower color. A separate S2 phenotyping experiment had 732 offspring (530 purple vs 202 white throat) and is not added to the 53 genotype-table denominator. Original authors' allele-specific expression in three F1 heterozygotes supports cis-regulation but no site-directed editing is performed.

**Gates, Olson, Clemente & Smith (2018, online 2017)**, DOI 10.1111/nph.14830, original author-hosted PDF https://www.colorado.edu/smithlab/sites/default/files/attached-files/gatesetal2017.pdf :
- Original Iochroma cross white loxense × blue cyaneum: backcross to blue 35 blue + 38 pale; backcross to white 8 white + 11 pale. One crossing program, not 92 independent losses.
- Novel R3 MYB inhibitor MYBL1 is associated with reduced anthocyanin structural-gene expression in white and pale pools.
- Original HETEROLOGOUS intervention: 20 transgenic Nicotiana tabacum plants/lines expressing I. loxense MYBL1, all white; authors report floral pigment 530-nm absorbance 11–20-fold lower across lines and 14-fold lower on average relative to wild-type. This demonstrates MYBL1 suppresses pigmentation when expressed in tobacco, NOT endogenous rescue or the exact causal nucleotide in Iochroma. Absorbance is not HPLC concentration.
- The native MYBL1 cis-regulatory explanation is a mechanism proposed from the expression/genetic results, not an independently verified genome-editing result.

**Meanchaipiboon, Kobayashi & Nakatsuka (2020)**, DOI 10.2503/hortj.UTD-142, https://www.jstage.jst.go.jp/article/hortj/89/3/89_UTD-142/_html/-char/en :
- White Hirado azalea cultivars show flower transcripts F3'H, F3'5'H, DFR and ANS. Author-reported white-petal flavonol and a trace anthocyanidin precursor in named exceptions.
- HPLC 6 wild samples + 25 cultivars, qPCR 6 wild samples + 18 cultivars. Cultivars are not 31 independent species; no verified matched individual RNA/HPLC pairs or functional transgene result. FLS is NOT among the four assayed genes. This is the same source already in v35, so cannot be counted as a new fourth original program.

## Evidence ladder; independent STUDY units, not probabilities

T1 = appearance + chemical assay, T2 = source-assayed expression, T3 = genotype/flower-colour segregation and allele-specific expression, T4 = heterologous regulatory perturbation with measured absorbance, T5 = DIRECT rescue in original white-flowered native species, T6 = NATURAL ancestral-state re-expression and descendant persistence.

Exactly 3 original paper programs, 3 families. T3 is documented in Ipomoea, T4 documented in Iochroma, and an organ-specific white-flower/stem biochemical contrast in Ipomoea. **T5 and T6 are not supplied by any of these programs.** Gene-expression detection does not imply enzyme functional capacity or future endogenous reactivation. Assay units, time, organ and response differ, making numerical effect pooling improper.

Exact quantitative count denominators are retained only for audit and must not be treated as independent evolutionary events or meta-analysis weights. Never infer that MYBL1 transformation into tobacco restores Iochroma, or that white-flower Ipomoea stems predict a return to purple flowers.

## Implementation / review

Freeze study IDs and DOI, organ and genotype/phenotype count tables, transformation and absorbance summaries as original-published descriptions. Create source registry, source-checked standard-library audit, adverse tests, deterministic derived-only result, comprehensive repository CI and Draft stacked on v35. No original source files, raw sequences or HPLC data are copied into the repository. If an original paper has no verified direct measure, the corresponding field MUST be absent/false rather than inferred.

**This is a new structured literature- and experiment-evidence audit, NOT independent raw-data reanalysis of identical gene-expression predictors in unseen species.**
