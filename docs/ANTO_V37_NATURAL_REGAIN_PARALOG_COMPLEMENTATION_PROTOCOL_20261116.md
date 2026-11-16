# ANTO v37: gene resurrection, paralog recruitment, experimental native-flower rescue and pollinator choice

**2026-11-16. Retrospective original-evidence audit after v36's three-program non-regain results.** All original papers/quantitative denominators are recorded as AUTHOR-PUBLISHED, not new sequence, HPLC, or insect behavior raw-data reanalysis. This is an independent next evidence domain; does not retroactively change v36's original three-study eligibility.

## Main question

When an ancestral flower colour is lost, can floral pigmentation reappear, and must the original lost gene be repaired? Disentangle four testable, different phenomena:
1. Natural historical regain of visible flower color inferred from published Petunia phylogeny/ancestral sequences.
2. Distinct molecular paths of regain (same inactive regulator resurrected versus paralog redeployed).
3. Experimental restoration of pigmentation in an originally white genetic background.
4. Ecological pollinator behavior after restoring floral color.

A natural species' continued existence is not a measured per-generation gain or persistence probability, and experimental petals in a greenhouse are not historical natural evolutionary reversions.

## Four original sources, original-study units

1. **Esfeld et al. (2018)**, *Current Biology*, DOI **10.1016/j.cub.2018.10.019**, title Pseudogenization and Resurrection of a Speciation Gene. Original authors infer *Petunia secreta*'s restored purple flower from a nonfunctional ancestral AN2; a **2-base-pair compensatory deletion** restores the reading frame. This is a natural lineages' inferred genetic reacquisition route, not a new experimentally watched ancestral-to-descendant time series or native CRISPR edit.
2. **Berardi et al. (2021)**, *Plant Cell*, DOI **10.1093/plcell/koab114**, title Complex evolution of novel red floral color in Petunia, original full text https://pmc.ncbi.nlm.nih.gov/articles/PMC8364234/ and https://academic.oup.com/plcell/article/33/7/2273/6237947 :
   - Red-flowered *P. exserta* retains nonfunctional AN2; **DPL** paralog spatially recruited to petal limbs, plus coordinated modifications to F3′H, F3′5′H and AAT/MYB-FL; the red chemical phenotype uses cyanidin and delphinidin, *not* simply red pelargonidin.
   - Native-species functional validation: DPL knockdown in red *P. exserta* produces white sectors. **Nine independent 35S:DPL stable lines were generated in naturally white *P. axillaris***, two lines selected for chemical analysis; HPLC petal-limb anthocyanidins measured in **3 nontransgenic controls and 6 transformed plants (3 per line)**. DPL restores violet/purple anthocyanin in *P. axillaris* corollas despite dysfunctional AN2. Count 9 generated lines and 2 assayed lines as different denominators. Laboratory rescue != historical ancestor reversal.
3. **Hoballah et al. (2007)**, *Plant Cell*, DOI **10.1105/tpc.106.048694**, title Single Gene–Mediated Shift in Pollinator Attraction in Petunia, original https://pmc.ncbi.nlm.nih.gov/articles/PMC1867374/ :
   - Restored a functional *P. integrifolia* AN2 by transformation into naturally white *P. axillaris*, changed limb coloration/reflectance.
   - Greenhouse pollinator choice: naïve **Manduca sexta** first chose **wild white P. axillaris 30** versus **AN2-transformed coloured 6** (n=36); author-reported binomial P<0.0001. **Bombus terrestris** preferred transformed-coloured flowers (n=12 observational trials, author Wilcoxon P=0.017). Morphology/scent/nectar measurements reported unchanged in the specific assessed conditions. This demonstrates *pollinator preference*, not field pollen delivery, seed production, or fitness.
   - This is the SAME white Petunia host species and ancestral AN2 pathway as Esfeld/Berardi: do not count each 2007/2018/2021 paper as independent macroevolutionary radiations.
4. **Fraser et al. (2013)**, *BMC Genomics*, DOI **10.1186/1471-2164-14-28**, title An R2R3 MYB transcription factor determines red petal colour in an Actinidia (kiwifruit) hybrid population, original https://pmc.ncbi.nlm.nih.gov/articles/PMC3618344/ :
   - The white-petalled Actinidia hybrid F2 backcross genotype **21-05-09e** was transformed at petal edges *transiently* with **35S::MYB110a**, restoring visible anthocyanin in cells near the infiltration site. Negative 35S::GFP control did not elicit pigment. This is a genuine **white petal functional complementation** in a separate plant FAMILY **Actinidiaceae**. A transient/local response in a hybrid, **not** a stable heritable transformation nor natural macroevolutionary reacquisition.
   - The original R2R3-MYB marker segregated with red/white petal appearance in four related hybrid families and the factor transactivated relevant promoters; do not invent per-family individual sample denominators from captions alone.

## Data model and fixed outcome categories

- **Natural regain channel:** `AN2_FRAME_RESTORATION` in *P. secreta* versus `PARALOG_DPL_REDEPLOYMENT_WITH_OTHER_PATHWAY_CHANGES` in *P. exserta*. Both are inferred within the **one shared Petunia long-tube radiation**, hence two molecular alternatives, not independent random draws of all plant evolution.
- **Lab native-white-petal restorations:** AN2 transgene in *P. axillaris*, DPL transgene in the SAME white *P. axillaris*, MYB110a transient complementation in white *Actinidia hybrid* petals. Exactly three distinct original perturbation programs, across exactly two flowering-plant families. One is transient/local; two are stable Petunia transgenics. Do not label a pigment-loss repressor (v36 MYBL1) as a functional rescue.
- **Ecological behavior:** only the 2007 AN2 experiment, controlled greenhouse moth and bumblebee *preferences*, not natural seed-set or long-term survival. Keep author-tested sample sizes and outcomes separate.
- **Mechanism specificity:** AN2 coding repair ≠ DPL subfunctional shift ≠ MYB110a activation; AN2 inactivation does not imply all anthocyanin enzyme loci broken. FLS expression is NOT necessarily different in Berardi's DPL overexpression assay; F3′5′H may remain low. Do not generalize that all genes are restored or same RNA module from ANTO v33 was independently quantitatively replicated.
- **No pseudo-replication:** Original publications=4; taxonomic families=2; independent Petunia radiation=1; natural regained flower-color tips=2 **from the one radiation**; native white petal functional intervention publications=3, but same Petunia host reused. Not 4 separate evolutionary transitions or 3 independent macro-biogeographic clades.
- **No fake probabilities:** Natural multi-generation regain rate and probability of retention NOT measurable from these studies. No root-resolved population sampling denominator for reversal frequency, no formal cross-study RMSE or meta-analysis. No direct field seed-set effect of AN2 restoration.

## Validation / implementation

Use fixed paper DOI + taxon + gene + readout + observation unit and strict counts in an original study registry; fail closed if natural and experimental categories are interchanged, if the 9 transgenic lines become 9 HPLC biological replicates, if individual moth choices are equated with multiple independent pollination events, or if *P. secreta* and *P. exserta* are counted as different family-level radiations. Emit derived-only, source-aware evidence ledger with negative and scope markers. Run new adversarial tests plus full ANTO repository CI. Keep original v26–v36 files unchanged, Draft stacked on v36.

## Source inference ceiling

This v37 directly documents **published naturally evolved regained flower colour and experimental restoration in white corollas**, stronger than v36's three original study set. It is still a source audit, not a new wet-lab result; the 2021 DPL and 2018 AN2 mechanisms share a clade's ancestral-state model. Do not claim that every white plant retains a restorable intact pathway, that restoration is always single-step or that pollinator choice automatically becomes realized offspring fitness.
