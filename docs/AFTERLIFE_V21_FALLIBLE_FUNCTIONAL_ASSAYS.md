# ANTO v21 — When fallible floral assays invalidate apparent evolutionary-capacity contrasts

**2026-10-26. Prospective measurement-validity mathematics, not an experiment on plants.** This extends the [real-source admissibility audit in v19](AFTERLIFE_V19_SOURCE_OBSERVABILITY_AND_ADMISSIBILITY.md) and the [22-taxon finite-panel identifiability work in v20](AFTERLIFE_V20_FUNCTIONAL_ASSAY_IDENTIFIABILITY.md). The original Ho & Smith (2016; DOI [10.1186/s12862-016-0675-3](https://doi.org/10.1186/s12862-016-0675-3)) crosswalk is 7 anthocyanin-unpigmented plus 15 pigmented taxa; all 22 lack qualifying newly measured jointly functional floral-toolkit assays in ANTO.

## Main new result

v20 correctly stated that with **hypothetical, perfectly valid, binary, taxon-level floral function tests**, 7 uniformly positive unpigmented taxa plus 1 negative pigmented taxon is a best-case eight-taxon witness for a strictly positive difference in this fixed 22-taxon panel. This has never happened as a biological experiment, and eight is not a field sample-size/power recommendation.

But the witness is fragile to even one possibly incorrect taxon-level assay call. An exact adversarial error model changes the best possible strict-sign identification thresholds as follows.

| *Hypothetical* maximum incorrect reported functional calls | Smallest **best-case** number of tested taxa for a strict finite-panel contrast | Best-case allocation |
|---|---:|---|
| No calls can be incorrect in either group | **8** | 7 unpigmented, 1 pigmented |
| Up to 1 incorrect call in pigmented group only | **9** | 7 unpigmented, 2 pigmented |
| Up to 1 incorrect call in unpigmented group only | **10** | 7 unpigmented, 3 pigmented |
| Up to 1 incorrect call **in each group** | **11** | 7 unpigmented, 4 pigmented |
| Up to 2 incorrect unpigmented and 1 pigmented | **13** | 7 unpigmented, 6 pigmented |

These are **necessary-and-attainable mathematical minima for a favorable pattern of hypothetical calls** within the declared finite panel and error limits, not experimentally expected numbers or guarantees that any observed sign will be certified. The mirrored favorable opposite-contrast design has the same respective minima.

In the previous eight-taxon v20 scenario, where all seven unpigmented taxa are reported positive and one pigmented taxon negative, **at most one wrong call across both groups (wherever it occurred) is already enough to destroy strict positive-sign certification**. With at most one wrong call in each group, the exact contrast interval is \([-1/7,1]\). If ideal reports instead include all seven unpigmented positive and four pigmented negative, the same group-specific one-call error caps produce a strictly positive lower bound of \(2/35\). No known actual floral assays are being assigned these calls.

## Biological target and limits

For only the **named 22 original source taxa**, let the binary *true* outcome be whether the complete set of functionally necessary floral developmental/biochemical resources was jointly validated in an appropriately controlled flower organ/stage. This is **not** the mere presence of a deposited CHI/F3H/DFR accession or a CDS feature. A source taxon is not an independent evolutionary loss or regain. Assigning a single binary status to a species further requires adequate multiple-individual, developmental-stage and population replication; this is assumed for the mathematical exercise, not accomplished by archived data.

At present there are **zero** admitted actual function outcomes. Thus the biologically relevant contrast in current complete floral competence still has sharp logical range \([-1,+1]\). It is scientifically incorrect to read an assay error model as a new plant observation.

## Derivation: deterministic maximum number of wrong binary calls

For a group with \(n\) original taxa, suppose \(t\) taxa have been **hypothetically** assayed under a validated complete-functional protocol, and \(a\) are reported positive while \(b=t-a\) are reported negative. Suppose an independent validation supplies a **deterministic upper bound \(e\)** on wrong binary calls *among these t reported results*, not merely an expected error rate.

The untested \(n-t\) statuses remain arbitrary. Let \(K\) be the true number of functionally positive taxa. The following finite-population bounds are sharp:

\[
\boxed{\max(0,a-e)\ \le K \le\
 a+(n-t)+\min(e,b).}
\]

The lower endpoint permits up to \(e\) reported positives to be false. The upper permits up to \(e\) reported negatives to be false, plus all untested taxa to be positive. Every integer endpoint is attained under an allowed assignment; explicit enumeration of all latent binary configurations independently confirms this for small populations in CI.

Apply separately to U (\(n=7\)) and P (\(n=15\)). With lower and upper true-positive counts \(K_U^L,K_U^H,K_P^L,K_P^H\), the fixed-panel prevalence difference \(D=p_U-p_P\) obeys

\[
\boxed{
D\in\left[
\frac{K_U^L}{7}-\frac{K_P^H}{15},\
\frac{K_U^H}{7}-\frac{K_P^L}{15}
\right].}
\]

Only if the lower bound is strictly positive (or upper strictly negative) does the *hypothetical* validated binary test evidence rule out the opposite sign. This is not phylogenetic statistical significance, a causal effect, ancestral regain probability, or long-run population survival.

The most optimistic one-call error robust configuration—7 U reported positive and 4 P reported negative, with \(e_U=e_P=1\)—has \(D\ge 6/7-12/15=2/35\), positive. Conversely the v20 eight-taxon configuration \(7U^+ /1P^-\) has \(D\ge 6/7-1=-1/7\), and no strict conclusion.

## An important distinction: sensitivity/specificity are not hard error-count caps

Suppose instead the per-call assay has published or independently estimated sensitivity and specificity between zero and one, for instance *hypothetical* 95% sensitivity and 95% specificity. Under the conventional full-support independent misclassification likelihood, **every possible underlying true binary state vector has a nonzero probability of producing any given observed set of positive and negative calls**.

Therefore those nominal accuracy probabilities alone cannot certify "at most one of 22 calls is wrong" as a logical bound. Without additional deterministic validation or a declared probabilistic statistical model and coverage target, the strict logical support for the full-panel contrast remains \([-1,+1]\), even after all 22 taxa have reported assay calls. A Bayesian posterior credible interval or frequentist confidence interval could narrow beliefs/coverage **under its own explicit assumptions**, but is not an all-worlds mathematical guarantee and is outside this v21 analysis.

This is also why repeatedly counting the same flower, gene transcript or colour-labelled photograph does not directly certify a functional species state. Actual replicate errors can be correlated through age, tissue, genotype, collection location and biochemical assay performance. The full-support observation model here is a logical counterexample, not an estimated laboratory error distribution.

## Why this matters for ANTO, EGWE and izu-core

The scientific bridge is **measurement-conditioned state adequacy**:

- ANTO v7–v18 separates visible floral phenotype, functional resource combinations, re-expression and descendant population persistence.
- ANTO v19–v20 documents why deposited sequence records do not independently measure the complete functional state.
- **v21** now shows that even an attempt to directly assay the relevant resource can fail to identify a *directional group difference* if classification error is not independently bounded or modeled.

A detected phenotype contrast, a pooled floral expression signature and a partly validated gene inventory cannot automatically be promoted to an observed inherited evolutionary potential. EGWE's state-sufficiency methods and izu-core's demographic/genetic models remain separate; no data or parameters from them were used to calibrate this source-specific construction. The same-clade Iochrominae 22 species are not 22 independent historical regain events.

## Preregistered assay-quality gates before empirical use

1. **Define the target in an individual, organ and time:** same-voucher flower pigment chemistry, complete pathway enzyme and regulator function, and measured developmental stage, not nominal molecular homology.
2. **Assess measurement error independently:** blinded reference controls, error rates by pigment expression group, taxon and tissue, replicate concordance, and stage sensitivity. A claimed deterministic wrong-call cap needs a real *design-backed* basis; "95% sensitivity" is not that basis.
3. **Resolve within-taxon heterogeneity:** independent genotype/population replicates and an explicit rule translating a sample's function to a taxon's outcome.
4. **Measure historical re-expression separately:** phylogenetic polarity, cross/recombination or native-repair versus donor introgression, direct phenotypic reappearance, descendant heritable reproduction and seed fitness are independent gates.
5. **Evaluate whole-clade holdouts:** do not infer a universal ability-retention rate from these 22 correlated taxa or from Petunia's related but distinct mechanisms. A real M2 functional-state forecast must beat visible/historical M0/M1 baselines in independent evolutionary lineages.

Actual Ryukyu Cirsium voucher assays, flower developmental/biochemical tests and population reproduction remain exclusively in [aza3](https://github.com/zuizui0223/aza3); nothing in this design authorizes a new sample collection or genetic manipulation.

## Reproduction

- Source-locked exact error-bound code: [scripts/anto_fallible_function_assays_v21.py](../scripts/anto_fallible_function_assays_v21.py).
- Negative and brute-force binary-truth tests: [tests/test_fallible_function_assays_v21.py](../tests/test_fallible_function_assays_v21.py).
- Machine-readable real-source no-assay statement and counterfactual witness fractions: [results/afterlife_v21_fallible_assay_design_20261026.json](../results/afterlife_v21_fallible_assay_design_20261026.json).
- [Dedicated CI](../.github/workflows/anto-v21-fallible-functional-assays.yml) checks the original v19/v20 evidence gates, all 4,896 original v20 label-count combinations, the new independent truth-vector oracle, and the absence of claimed new biological functional measurements.

**Publication ceiling:** exact source-anchored mathematical sensitivity of a currently **unmeasured** floral functional-competence group contrast to potential assay errors. No actual molecular functionality, re-evolution, adaptation, retention half-life, or organismal fitness was estimated.
