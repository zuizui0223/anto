# Penstemon original experimental counterexample: structural decline versus donor-mediated recovery

**Original empirical study:** Wessinger & Rausher (2014), *Evolution* 68:1058–1070, [10.1111/evo.12340](https://doi.org/10.1111/evo.12340). Published original data at [Dryad 10.5061/dryad.sp800](https://doi.org/10.5061/dryad.sp800). **This is a source-based re-expression hypothesis control, not a new field experiment or an independent evolutionary fitness study by ANTO.**

## Why this system is essential to the afterlife question

Unlike Cirsium's unresolved white versus coloured contrast and Iochrominae's conserved coding following loss of *floral* anthocyanin, the original Penstemon comparison concerns **loss of the ancestral blue-purple/delphinidin pathway branch while the flowers remain red and anthocyanin-pigmented**.

The original paper's endogenous F3′5′H enzyme is structurally inactivated in *P. barbatus*. A stop codon at amino acid 407 removes the downstream heme-binding domain; in vitro enzyme activity is undetectable. Reversing **that stop alone** failed to restore activity, supporting at least one other loss-of-function change. This is genuine **multiple-hit biochemical inactivation**, not merely suppressed expression. But the genes enabling blue-purple colour in the other parental species can also be introduced in the interspecies cross.

**This establishes three distinct objects:** the red visible phenotype, the damaged endogenous genetic pathway, and the accessible donor pathway under an artificial cross. A locus can have severely impaired intrinsic reversibility without making the phenotype inaccessible under external gene flow. The latter was **tested in crosses, not observed as naturally occurring introgression or re-evolution in the wild**.

## Original Table 1: 161 source-authored F2 individuals

Alleles F/f identify the original donor species: F is from blue-purple *Penstemon neomexicanus* (functionally intact F3′5′h), f from red *P. barbatus* (damaged F3′5′h).

| Source F2 flower colour | FF | Ff | ff | Total |
|---|---:|---:|---:|---:|
| Blue-purple | **57** | **65** | **0** | **122** |
| Red | **0** | **15** | **24** | **39** |
| **Total** | **57** | **80** | **24** | **161** |

The original phenotype response is **NOT deterministic given this single locus alone**, because among **80 Ff heterozygotes, 65 are blue-purple (81.25%) and 15 remain red (18.75%)**. Source authors interpreted the discrepancy using genotype-dependent additional background factors. It is not evidence of an ANTO-calibrated parameter, nor is the F2 offspring frequency the probability of evolved colour restoration after a natural loss.

[Machine original source transcription](../data/external/penstemon_2014_f3p5h_source_genetic_accessibility_v1.json), [exact count validator](../scripts/audit_penstemon_2014_afterlife_contrast_v1.py), [negative controls](../tests/test_penstemon_2014_afterlife_contrast_v1.py).

## Cross-clade *mechanism types*, not pseudo-independent evolutionary probabilities

| Clade | Paper-derived supported biological phenomenon | Crucial unverified target |
|---|---|---|
| *Petunia*, Solanaceae | Documented natural colour regain with distinct regulatory/coding paths within one genus (Esfeld 2018; Berardi 2021) | General loss→regain probability, quantitative persistence |
| Iochrominae, Solanaceae | Original Ho–Smith 2016 reported Chi/F3h/Dfr coding conservation despite loss of floral anthocyanins; independently retrieved 58 public accession records with 22 submitted taxa | Direct functional inducibility and inherited return in comparable full individual datasets |
| *Penstemon*, Plantaginaceae | Wessinger & Rausher 2014 source authors documented multiple inactivating variants at F3′5′h, no activity after single-stop correction, and recovery of blue-purple colour among interspecific crosses via donor alleles | Does NOT demonstrate endogenous backmutation restoration, natural ongoing introgression, or natural long-run persistence |
| *Cirsium*, Asteraceae | Public leaf transcript and reference family tree outcomes are unresolved for floral function (aza3 owns empirical evidence) | Functional joint pigment pathway, floral activation, ancestry direction, heritable return, viable achenes |

Different traits are contrasted: white versus anthocyanin-positive in Iochrominae, red versus blue-purple branch in Penstemon, more diverse Petunia regain mechanisms. These are biologically related **pathway-level resource states**, but **not one homogeneous binary flower-colour loss outcome**. A general model must represent which function is lost and which is retained; otherwise the comparisons are pseudoreplication across incomparable targets.

## New prospective, falsifiable comparative test (not carried out)

**Hypothesis:** When the original target pigment pathway becomes silent, future accessible states depend more strongly on the *combined molecular mechanism* (intact coding, regulatory repression, multi-hit loss, external allele availability) than on the current visible hue alone.

Cross-system estimands should be indexed by target pathway/anthocyanidin and comparable assay scope. Candidate models: visible pigment alone; visible + complete functional coding and regulation; visible + ancestry/cytotype and donor gene-flow; full model + experimental expression/viable-offspring data. Evaluate independent clade/population holdouts; include damaged structural negatives and still-competent white controls.

**No universal time-to-irreversibility claim is identifiable today.** Penstemon's documented multiple-hit structural damage is an existence case, not a hazard rate, and PETUNIA and Iochrominae are in the same larger family. ANTO v10/v10b parameters are synthetic and must not be fitted from the source Table 1 as if 161 observations were 161 independent historical loss/reversal episodes.

## Reproduce

    python -m unittest discover -s tests -p test_penstemon_2014_afterlife_contrast_v1.py -v
    python scripts/audit_penstemon_2014_afterlife_contrast_v1.py --output /tmp/penstemon-source.json

The original source is article Table 1, deliberately pinned in a clear provenance file. The five original Dryad spreadsheet files were **not downloaded or fitted** in this analysis; any later use requires a separate source-verified snapshot.
