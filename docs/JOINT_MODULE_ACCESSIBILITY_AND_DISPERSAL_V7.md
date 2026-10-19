# anto v7 — The hidden joint inventory of floral evolutionary potential

**Date:** 2026-10-19  
**Scope:** Exact probability bounds and explicit synthetic dispersal counterfactual, independent of any particular angiosperm genus. It is **not** a new probability theorem or evidence that real thistle florets have re-evolved.

## Main contribution: first- and second-order module frequencies do not determine complex evolutionary accessibility

Let three necessary developmental/genetic resources be \(M_1,M_2,M_3\), where a 1 means **functionally verified capacity is retained**, and not simply that a homologous RNA read was found. Let \(G=(M_1,M_2,M_3)\) be the **joint resource inventory carried by a propagule**.

Suppose in each geographic region or plant lineage we know with infinite precision:

\[
P(M_i=1)=1/2,\quad P(M_i=1,M_j=1)=1/4\quad(i\ne j).
\]

These are the same marginals and pairwise joint retention frequencies as independent fair binary modules. But the probability that *all three necessary resources coexist in the same lineage* can range over the sharp interval:

\[
\boxed{0\le P(M_1=M_2=M_3=1)\le 1/4}.
\]

This is the Fréchet feasibility result for three Bernoulli moments, **not** a newly discovered theorem. The important biological application is that **even perfect pairwise association analyses cannot identify whether the required three-part floral toolkit is present**. This is not resolved by enlarging a flower-photo dataset; the same plant's joint, functional modules must be measured.

### Four exact constructed inventory regimes

| Joint inventory model | Marginal retention of each module | Pairwise co-retention of each pair | All three jointly retained | Upper bound E[2^K] |
|---|---|---|---|---|
| Even parity (000, 011, 101, 110; each 1/4) | 1/2 | 1/4 | **0** | 13/4 = 3.25 |
| Independent (all eight masks each 1/8) | 1/2 | 1/4 | **1/8 = 12.5%** | 27/8 = 3.375 |
| Odd parity (001, 010, 100, 111; each 1/4) | 1/2 | 1/4 | **1/4 = 25%** | 7/2 = 3.50 |
| Complete positive linkage (000/111 each 1/2) | 1/2 | 1/2 | **1/2 = 50%** | 9/2 = 4.50 |

The first three distributions have **identical first- and second-order statistics**; even parity and odd parity even have the same joint genotype entropy (2 bits). But they differ maximally in a three-component target's **necessary resource availability**.

The upper bound E[2^K] is the expected number of nominal ON/OFF subsets if every retained module could be freely toggled. It is **not measured head phenotype richness**: developmental incompatibility, linkage, dosage constraints, epistasis, stage and ecological selection can eliminate many nominal configurations.

The corresponding exact third-order centered moment is -1/8, 0 and +1/8 for even, independent and odd parity, respectively. It must not be replaced with a proxy derived from pairwise correlations.

### General sharp bounds (arbitrary marginal and pairwise data)

For \(m_i=P(M_i=1)\), \(p_{ij}=P(M_i=M_j=1)\) and \(t=P(111)\):

\[
\begin{aligned}
t_{\min}&=\max(0,p_{12}+p_{13}-m_1,\ p_{12}+p_{23}-m_2,\ p_{13}+p_{23}-m_3),\\
t_{\max}&=\min(p_{12},p_{13},p_{23},\ 1-m_1-m_2-m_3+p_{12}+p_{13}+p_{23}).
\end{aligned}
\]

Every \(t\) within that interval gives a feasible full joint distribution; outside it, the candidate marginals/pairs are inconsistent. These bounds generalize the toy counterexample. They do not infer genotype accessibility from floral images.

## Correlation change is not automatically adaptive modular recombination

Introduce an abstract mixing operator \(R_\rho(P)=(1-\rho)P+\rho \prod_i P(M_i)\). The word *mixing* means a population-level redistribution of joint inventory toward marginal independence. It is **not a measured Mendelian recombination fraction**, and may not even be biologically reachable with limited fertility or genotype assortment.

With \(\rho=1/2\), exactly the same operator changes complete-resource availability as follows:

- Even-parity population: **0 → 1/16 = 6.25%**.
- Odd-parity population: **1/4 → 3/16 = 18.75%**.
- Positive-linkage population: **1/2 → 5/16 = 31.25%**.

Thus decoupling modules can create a previously inaccessible joint configuration **or destroy the frequency of an already favourable linked inventory**. This is a conditional statement about a specified probability operation, not evidence that real genetic recombination has these effects.

## Why moving plants, not just moving pollinators, matters

Suppose a source and a newly colonized region have identical observed white flowers and the same initial joint inventory. Let the structural module loss hazard per generation differ among regions. With no repair, no gene flow and no new synthesis, the probability that all \(k\) required modules survive a period \(T\) is

\[
q_{\mathrm{complete}}(T\mid r)=q_{\mathrm{complete}}(0)
\exp\left[-T\sum_{i=1}^k\delta_{i,r}\right].
\]

If one lineage moves to the new region **before** the \(T\)-generation waiting period and another moves **after** it, their relative chance to retain all resources is

\[
\boxed{\frac{q_{\mathrm{early}}}{q_{\mathrm{late}}}
=\exp\left[-T\sum_{i=1}^k(\delta_{i,\mathrm{new}}-\delta_{i,\mathrm{source}})\right]}.
\]

When the destination-specific hazard exceeds the source by the same positive amount for each module, the **relative** disadvantage accumulates with the number of necessary components \(k\). The **absolute** difference need not increase monotonically, and the sign reverses if the destination protects resources better. With no complete inventory initially (even parity), neither movement order permits a complete configuration absent rebuilding or genetic mixing.

Illustrative numbers only: \(T=10\), independent initial module retention 1/2, structural hazards 0.02 (source) and 0.12 (destination), no repair. Retained complete-inventory probabilities are:

| Required components k | Early migration | Late migration | Early/late ratio |
|---:|---:|---:|---:|
| 1 | 0.15060 | 0.40937 | 0.36788 |
| 2 | 0.02268 | 0.16758 | 0.13534 |
| 3 | 0.00342 | 0.06860 | 0.04979 |

For \(k=3\), eight **independent biological founders**, hypothetical developmental assembly probability 0.7 and phenotype establishment probability 0.5, the probabilities of **at least one successful realized trait** are about **0.95%** early versus **17.67%** late. These are *not estimates of thistle traits* and rely on unrealistically simple independence and time-invariant selection.

## The anto ↔ aza3 hypothesis stack

**anto theorizes resource availability under plant history:**

\[
P(G_{\mathrm{joint}}\mid\text{plant dispersal history, drift, protection, recombination, genealogy}).
\]

**aza3 measures the downstream LEGO assembly:**

\[
P(\text{capitulum configuration }Z\mid G_{\mathrm{joint}},\text{development},E)
\]

and finally a separate ecological fitness contrast in viable achenes. Same visible white flowers may conceal distinct inventories; even a complete inventory may fail to produce a trait because floral regulators remain off or because module combinations are incompatible.

This is not a scientific identity between "trait regain" and "LEGO". Re-expression of a retained pathway is one route; de novo mutation, migrant coloured propagules, introgression, and non-genetic plasticity remain distinct. A functional module cannot be classified from leaf transcript similarity alone.

### Repository role freeze

- **anto** remains the cross-angiosperm theoretical and prospective-comparison home. No new *Cirsium* empirical sample source is interpreted or developed here after v6.
- **[aza3](https://github.com/zuizui0223/aza3)** is the canonical **Cirsium empirical** owner. This includes Ryukyu white/coloured candidates, *takaoense* 3+3 source morphs, the *C. sieboldii* genetic LEGO test, all head modules, pigment/regulatory assays, phylogeny and fitness. See [aza3 empirical bridge](https://github.com/zuizui0223/aza3/blob/main/docs/ANTO_ACCESSIBILITY_TO_AZA3_HEAD_LEGO_EMPIRICAL_MASTER_V1.md).
- **izucore** informs history/order/founder sensitivity but its simulated effect sizes are not Cirsium parameters.
- **egwe** provides the measurement-sufficiency criterion; here even perfect first- and second-order measurement is demonstrably insufficient.
- **chun** gives visible-colour evolutionary memory across radiations, not measured floral pathway competence. It warns against a universal coarse flower-colour scale.

## What is and is not demonstrated

**Proved under the explicit mathematical assumptions**:
1. Sharp third-order inventory bounds given all first-/pairwise retention probabilities.
2. A shared pairwise distribution and shared entropy can hide distinct necessary composite resources.
3. One abstract decorrelation operation can have opposite effects on complete resource availability.
4. Geographic migration order has a relative effect that is multiplicative across independently degraded necessary modules when regional loss hazards differ.

**Not established in living plants**:
- Actual functional inventory correlations in any flowering-plant clade.
- Heritable reactivation of floral colour in the Ryukyu *Cirsium* sister lineages.
- Recurrent genetic LEGO assembly from old alleles or an ecological fitness benefit.
- A universal evolutionary time window, globally calibrated decay rate or a common direction of colonization effect.

The most valuable *biological* tests will come from independent, same-individual, multi-organ, multi-module genetic studies and whole-clade replication beyond Cirsium. A photos-only macroanalysis may locate visible colour boundaries but cannot identify these latent mechanisms.

## Reproduce

\`\`\`bash
python -m pip install numpy scipy
python -m unittest discover -s tests -v
python scripts/run_joint_accessibility_v7.py --output results/joint_accessibility_v7.json
\`\`\`

The v7 source and unit tests are deliberately **independent of aza3 source files**, so mathematical generality cannot be inflated by reuse of one biological cohort.
