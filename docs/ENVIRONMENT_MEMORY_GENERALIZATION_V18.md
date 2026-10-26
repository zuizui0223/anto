# ANTO v18 — Limits of the environmental-memory × latent-capacity reversal

26 October 2026. **Exact conditional synthetic branching model only.** This extends [v17](ENVIRONMENT_HORIZON_SUFFICIENCY_V17.md) while preserving the two-reproductive-slot expressing/silent-competent kernel from [v4](SLK_ANTO_RECURRENT_CAPACITY_V4.md) and [v6](SLK_ANTO_MARKOV_ECOLOGY_V6.md). It does not establish a general law of ecological adaptation, measured genetic retention, natural ancestral phenotype regain, plant reproductive fitness, organism extinction or genetic fixation.

## Main question

The original ANTO v17 constructed a symmetric 50:50 pair of opposite ecological conditions in which increasing environmental persistence raises the chance of retaining toolkit-competent descendants at generation two when expressing/silent switching w is less than 0.5, lowers it when w is greater than 0.5, and has no effect at w=0.5. **Are those sign and null claims robust to unequal phase frequencies, non-reciprocal D/C reproductive responses and reproductive-output saturation?**

The answer is model-specific: the second-generation effect is affine in environmental lag-one autocorrelation at fixed stationary occupancy, but its sign has **three separable factors**. The w=0.5 null is special to the v17 witness. Furthermore, two-generation equality need not persist to generation three.

## The biological unit is a genealogy of retained competence

- One initial **nonexpressing yet reproductively active and completely toolkit-competent C parent** begins each trajectory; this is not a dormant seed.
- Descendants may express the relevant candidate trait (D) or remain silent with complete toolkit competence (C). Children losing competence are excluded from the tracked *competent genealogy*, not assumed to be biologically dead.
- Endpoint S_H: **at least one complete-toolkit descendant at exactly generation H**. It is not first visible expression, an event ever occurring before H, eventual branching nonextinction, organism population persistence or genetic fixation.
- In environmental phase e each D/C parent has two independent reproductive slots, producing a viable child with probability q_D,e or q_C,e. The child retains the toolkit with common probability k and switches expressing/silent state relative to its parent with probability w.
- A or B environment is shared by **every descendant in a given generation**. An exogenous stationary two-state Markov chain drives it; independent weather per offspring, ecological feedback, migration, mutation, introgression and pathway reconstruction are excluded.

## Exact non-symmetric environment theorem for the declared kernel

Let π_A=π and π_B=1−π with 0<π<1. Define lag-one environment autocorrelation ρ. For strictly positive ergodic phase transitions,

\[
T_{ij}=(1-\rho)\pi_j+\rho\,\mathbf{1}_{i=j}.
\]

Changing ρ holds each phase's stationary occupancy constant. Negative ρ is permitted where all probabilities remain positive. The symmetric v17 persistence p is (1+ρ)/2 when π=1/2.

Define

\[
u_i=kq_{C,i},\qquad
h_j=w[1-(1-kq_{D,j})^2]+(1-w)[1-(1-kq_{C,j})^2].
\]

Here u_i is a single C-parent reproductive slot's probability of creating a competent child. The conditional mixture h_j is that child's chance of producing at least one further competent descendant in phase j. The original C parent's probability of **zero competent grandchildren** given environments (i,j) is

\[
f_{ij}=(1-u_i h_j)^2.
\]

Using the **shared** environment-path probabilities,

\[
S_2(\rho)=1-\sum_{i,j}\pi_i T_{ij}f_{ij},
\]

and the exact derivative at fixed π is

\[
\boxed{\frac{\partial S_2}{\partial \rho}=
\pi(1-\pi)(u_A-u_B)(h_A-h_B)
[2-(u_A+u_B)(h_A+h_B)]}.
\]

The identity is derived from the classical two-slot branching PGF. It is an exact consequence within this parameterized model, not a novel universal probability theorem.

The three factors have distinct ecological meanings: (1) environment-specific **initial C reproduction**; (2) environment-dependent competent offspring reproduction conditional on D/C state; and (3) a **nonlinear at-least-one-offspring saturation** term. Positive initial and later reproductive contrasts do not necessarily imply a positive persistence effect if the saturation term is negative. This is a binary descendant-event effect and does not prove that reducing toolkit retention is beneficial in real organisms.

The model computes all 2^H shared-environment paths with exact fractions for H up to 4, checks the independent two-generation environment-pair expression, and validates it against explicit enumeration of two initial reproductive slots, each potentially producing no child, an expressed child or a silent competent child. All results are exact model predictions, not fitted parameter estimates.

## Witness 1: v17 recovered without changes

Set π=1/2, k=4/5, q_D,A=9/10, q_C,A=2/5, and reverse q_D and q_C in B. The offspring-state switching null is w*=1/2 and S1=456/625=72.96%. Raising ρ from −0.8 to +0.8 (v17 p=0.1 to p=0.9) changes S2 by **+5789184/244140625**, or +2.3712497664 percentage points, at w=0.1. Its sign reverses at w=0.9 and vanishes at w=0.5. This is regression recovery of the same constructed witness, **not an independent biological replicate**.

## Witness 2: the critical switching value shifts and a two-horizon null fails at three

Use nonreciprocal phase-specific outputs (q_D,A,q_C,A)=(0.85,0.35) and (q_D,B,q_C,B)=(0.45,0.75), k=0.8, with unequal stationary phase occupancy π_A=0.30. Write g_D,e=1−(1−kq_D,e)^2 and g_C,e=1−(1−kq_C,e)^2, and define Δg_D=g_D,A−g_D,B and Δg_C=g_C,A−g_C,B.

Only when biologically admissible, the interior cancellation point for the h_A−h_B factor is

\[
w^*=\frac{\Delta g_C}{\Delta g_C-\Delta g_D}.
\]

For these particular outputs, **w*=7/13≈0.53846**, rather than 0.5. At exactly w=7/13, the two synthetic environmental histories ρ=0 and ρ=0.8 give the **same** first- and second-generation competent-descendant probabilities. Yet they diverge at the third:

\[
S_3(\rho=0.8)-S_3(\rho=0)
=-\frac{492874414932892640256}{1513399183750152587890625}
\approx-0.000325674.
\]

The difference is only **−0.03257 percentage points**. Two horizons of apparent sufficiency do not certify the third even in this exact Markov model. With equal hypothetical prior weights, the history-label Bayes-Brier gain is exactly zero at horizons one and two but becomes strictly positive at three. The magnitude is extremely small: it is **not** an assurance of practical field statistical power.

There is no general interior w* when D/C environment contrasts do not oppose one another; cancellation can instead arise from identical founding C output or a zero saturation factor. The precise numerical value w*=7/13 is not a transferable plant trait threshold.

## Witness 3: even w=0.5 may preserve an environmental-memory effect

Set q_D,A=q_C,A=0.9, q_D,B=q_C,B=0.4, k=0.8, π=1/2 and **w=0.5**. In each phase, D and C have identical reproductive outputs, so expressing/silent switching does not mediate relative fecundity. Nevertheless, raising ρ from 0 to 0.8 increases S2 by **+1.482031104 percentage points**.

This is a **generic shared-environment demographic correlation effect**, not evidence for ancestral latent trait or phenotype-state inheritance. It is a crucial confounder: observing that environment persistence predicts competent-genealogy presence does not identify the D/C-memory mechanism unless phase-specific D/C reproduction and competence are measured separately.

## Witness 4: the nonlinear saturation factor can reverse the sign

Hold the same phase outputs in two counterfactual models, (q_D,A,q_C,A)=(0.98,0.98) and (q_D,B,q_C,B)=(0.85,0.85), with w=0.5 and π=1/2. The early and later reproductive contrasts are both positive. Changing only common toolkit retention k changes the sign of the exact two-generation environmental-autocorrelation response.

| Per-child toolkit retention k | Two-slot saturation factor | Change in S2 when ρ goes 0→0.8 |
|---|---|---|
| 0.35 | positive | +0.07383 percentage points |
| 0.95 | negative | −0.11202 percentage points |

This is the nonlinear saturation of the **binary at-least-one competent descendant** event, not measured physiological exhaustion, genetic reversion or evidence that deleting a pigment gene is adaptive. It is a separate mechanism from ANTO v15's explicit failed-expression shared-capacity costs.

## What the comparisons mean scientifically

The most defensible general statement is that **environment autocorrelation can affect future retained competence, but the direction and explanatory mechanism depend on founder reproduction, the future expression-state mixture, and nonlinear endpoint saturation**. Thus the original v17 w=0.5 null is not a universal biological invariant.

The state-sufficiency interpretation is also conditional. If histories with a matching observable present state differ only in an unmeasured future environmental transition operator, history labels can predict longer-horizon outcomes as proxies. A complete current biological state together with that transition law removes the proxy under the Markov model; residual history after incomplete molecular assays is not evidence of an irreducible inherited memory.

## Prospective biological falsification

Independent tracked lineages—not repeated flowers, heads, photographs or independent species averages—must be followed across predeclared horizons. Measure: complete **functional** coding and regulatory toolkit co-retention in the same organisms, the expressing/nonexpressing state transition of that same floral structure, partner-environment marginal frequencies **and autocorrelation**, absolute reproduction of both D and C within each ecological regime, and the genotype/functional competence of descendants. Include ancestry, introgression, ontogenetic stage, seasonality and plasticity controls.

For held-out populations, compare an averaged-environment baseline, a **phase-transition-only demographic model** with no D/C relative-fecundity contrast, and the full **D/C × phase** model. The difference between the last two is necessary before attributing a temporal effect specifically to ancestral floral functional memory. A real phenotype return must also be distinguished from donor-gene reacquisition and immediate but nonheritable plasticity.

ANTO owns only this generic probability model. Empirical Cirsium genetics, fieldwork and reproductive endpoints belong to [aza3](https://github.com/zuizui0223/aza3); SLK has a separate structural architecture payoff definition; EGWE has independent state-validity tests; izu-core has distinct genetic-accessibility and demographic simulations. None of their numerical results were used to fit these synthetic probabilities. The original v4/v6 documentation acknowledges prior literature on switching and multitype branching in fluctuating environments; this v18 extension does **not** claim new general mathematical priority.

## Reproduce

- Exact generic model: src/anto_environment_memory_generalization_v18.py
- 12 independent slot-enumeration, identity and negative-control tests: tests/test_environment_memory_generalization_v18.py
- Five-witness exact fraction receipt: scripts/run_environment_memory_generalization_v18.py
- Dedicated CI: .github/workflows/anto-v18-memory-generalization.yml

The GitHub Actions gate rejects non-matching H2 analytic identity, false H2→H3 sufficiency, a universal w=0.5 null, incorrect saturation-sign claims and unsupported natural plant/fitness/fixation extrapolations.
