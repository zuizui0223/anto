# anto v9 — Does evolutionary history still matter when the present resource state is known?

Date: 2026-10-19  
**A new falsifiable theoretical framing, NOT a demonstrated plant evolutionary law.**

## Central distinction: historical dependence versus hidden-state mediation

The claim "a flower remembers its evolutionary history" can be misleading. For a Markov-complete present biological state \(X_t\), a future trait \(Y_{t+H}\), and past dispersal/ancestry \(H_{<t}\),

\[
P(Y_{t+H}\mid X_t,H_{<t})=P(Y_{t+H}\mid X_t).
\]

The past can still have important predictive value when the measured present state is incomplete, because past dispersal, founder sorting, ploidy changes, gene flow, tissue-specific repression or structural degradation altered a *hidden* joint functional resource inventory.

The two compatible hypotheses we must separate are:

- **H-history-proxy:** including the dispersal/lineage history improves a future trait prediction only because the present visible flower hides necessary developmental/genetic resources.
- **H-residual:** history predicts future outcomes even when an apparently complete state was measured. This means our proposed state may have omitted allele age, molecular quality, epigenetic or ontogenetic state, population demography, gene exchange or relevant environmental exposures. It **does not automatically prove a non-Markov evolutionary law**.

The target is a state-sufficiency test informed by egwe, but explicitly about **plant propagule/allele movement and future floral phenotypic accessibility**, not pollinator-driven shifts as the primary moving protagonist.

## Exact three-module generative example

Choose two equally weighted hypothetical plant histories (history A, history B). They produce identical visible WHITE flowers and, even in the latent functional modules \(M_1,M_2,M_3\), both have:

\[
P(M_i=1)=1/2,\qquad P(M_i=M_j=1)=1/4\ (i\ne j).
\]

A and B differ only in their higher-order resource composition:

- A: even parity genotypes 000/011/101/110, each 1/4; **P(111)=0**.
- B: odd parity genotypes 001/010/100/111, each 1/4; **P(111)=1/4**.

All current flowers are WHITE regardless of their latent inventory. A future composite phenotype can be established with probability \(\alpha=0.4\) **conditional on a single lineage retaining ALL THREE required resources**; it is impossible in this toy model otherwise. Therefore:

\[
P(Y=1\mid A)=0,\quad P(Y=1\mid B)=0.1,\quad P(Y=1)=0.05.
\]

These probabilities are not measured flower re-expression or known biological establishment rates. Both histories are constructed prior distributions, not dated reconstructions of any taxon or locale.

### Exact Bayes-optimal prediction errors

Expected Brier error for the true conditional event probability \(p\) is \(p(1-p)\). Computing the population conditional forecasts exactly, not fitting simulated data, yields:

| Information available at prediction | Exact expected Brier error |
|---|---:|
| Visible white-flower state alone | 0.047500 |
| Visible colour + history label | 0.045000 |
| Perfect direct three-module functional resource inventory | 0.030000 |
| Perfect resources + history | 0.030000 |

The gain from history is fully mediated by \(G_{\mathrm{joint}}\): knowing the COMPLETE current resource state makes the ancient history label redundant under this generative model. Thus an empirical result "history seems important" should not automatically be interpreted as an irreducible, unmeasurable evolutionary memory effect.

## Molecular assays are imperfect: a measurement-design prediction

Assume each functional module is tested by an independent noisy binary assay with known sensitivity and specificity. Perfect functional genomic assays do not exist merely because gene-family transcripts were sequenced.

With sensitivity 90% and specificity 95% stipulated (NOT validated in Cirsium):

| Available information | Bayes-optimal Brier error |
|---|---:|
| White flower only | 0.047500 |
| White flower + noisy three-module assays | 0.037183 |
| White flower + history label + noisy assays | 0.033187 |
| White flower + fully observed functional inventory | 0.030000 |

**Counterintuitive implication:** imperfect molecular measurement can leave substantial independent value in historical ancestry/colonization information, even if a fully measured biological state would mathematically make history redundant. Increasing the number of photographs cannot fix that measurement gap.

If assay sensitivity and specificity are both 50%, the assay adds no information; if both are 100%, history adds precisely no predictive improvement beyond fully measured resources. In the null where both historical groups have the *same* current joint inventory, history provides no prediction advantage at all. If structural modules coexist but developmental assembly probability is zero, no floral emergence is predicted despite apparent resource retention.

These predictions are exact and falsifiable under the stipulated model; they are not validation of a real species' functional molecular markers. The known-assay assumptions matter.

## A pathway to a genuinely biological cross-clade result

**Prediction A — hidden geography.** Regions with the same visible white-flower composition and equal per-module retention rates may still differ in the presence of jointly retained *functional* developmental pathways. Such differences should be related to independently reconstructed plant colonization/founder histories, not simply assigned from a white-flower map.

**Prediction B — state sufficiency.** If validated same-individual multi-organ functional module assays are added to a history-aware model, the held-out predictive value of ancestral dispersal history should shrink when the measured modules account for the mechanism. A residual should trigger explicit tests of missing genetic quality/regulatory state or non-Markov environment, not automatic invention of a "memory" gene.

**Prediction C — mechanism-specific exceptions.** Complete structure is still not sufficient when flower-specific activation, developmental compatibility, seed establishment and local ecology are absent. Introgression or de novo reconstruction can reopen pathways that the toy's no-repair transition matrix calls "closed."

**Operational empirical validation** requires independent lineages with time-aware colonization/pigment histories, functional genetics in appropriate floral and nonfloral tissues, proper ancestry and ploidy control, and held-out whole populations/clades. At least two additional independent plant clades besides Cirsium and Petunia should be sought before using a broad cross-angiosperm generalization; no fixed numerical "universal loss rate" should be assumed.

Current known studies have complementary but nonidentical evidence levels: Petunia is a documented molecular re-acquisition control, aza3 Cirsium is an unresolved historical colour contrast and prospective head LEGO mechanism, and other published angiosperm colour histories are not automatically measurements of latent floral competence. Missing values must remain missing.

## Reproducibility and claim ceilings

- Exact probability kernel: [src/anto_history_state_sufficiency_v9.py](../src/anto_history_state_sufficiency_v9.py)
- Negative controls and no-assay/no-success tests: [tests/test_history_state_sufficiency_v9.py](../tests/test_history_state_sufficiency_v9.py)
- Frozen numerical scenarios: [scripts/run_history_state_sufficiency_v9.py](../scripts/run_history_state_sufficiency_v9.py)

Run:

    python -m unittest discover -s tests -p 'test_history_state_sufficiency_v9.py' -v
    python scripts/run_history_state_sufficiency_v9.py

This is synthetic *conditional probability theory*, not empirical Cirsium sequencing, molecular orthology inference, field phenology, natural population history, trait regain or a causal fitness finding. **aza3 alone** is the active owner of Cirsium empirical data and downstream genetic LEGO experiments. The anto v9 test suite requires no Cirsium files.
