# anto model and identifiability v1

## Explicit states
P = expressed pigmented flower.
W+ = non-expressing flower with experimentally supported capacity for re-expression under a specified intervention/condition.
W? = white flower of unknown competence (default for photo-only data).
W- = non-expressing flower with strong evidence of reduced competence under a specified assay; **not** absolute irreversibility.

Separate lineage ancestry, gene flow, and regional occupancy from these phenotype/competence states. Environmental plasticity is an observation/emission process, not automatically an evolutionary transition.

## Baselines
M0: geography and visible phenotype only.
M1: dispersal, local extinction and ordinary observed phenotype transitions.
M2: M1 plus hidden competence retention, degradation and reacquisition.
M3: M2 plus introgression, hybrid ancestry and alternative regulatory routes.

Do not rank M2/M3 by fit to visible data alone if parameters are non-identifiable. Use parametric simulation, posterior predictive checks, prior sensitivity, phylogenetic uncertainty, and held-out independent clades.

## Accessibility estimand
A(x,t,H) = Pr(locally derived re-expressed phenotype establishes within H generations | full state x at time t, specified environment, migration/gene flow regime).
G(x,t,H) = Pr(coloured phenotype occupies a geographic region within H generations | x).
Typically G and A are different, even when both are observable as coloured flowers at a later date.

## Minimal competing-risk illustration
For one W+ lineage with constant rates r (successful re-expression) and delta (irreversible competence degradation), the conditional probability of re-expression over H among populations still white at time t after initially becoming W+ is:
R(t,H) = [r exp(-(r+delta)t) / (delta+r exp(-(r+delta)t))] * [1-exp(-(r+delta)H)].
This requires constant rates, no repair, no gene flow, no branching and conditioning on white persistence. It is a teaching baseline, **not** an empirical thistle estimate or novel theorem. Note that r is an effective establishment-event rate, not per-allele mutation probability.

## Key theoretical question
Do dispersal and degradation operators commute? If the degradation rate depends on region or other-organ use, applying dispersal then degradation need not yield the same distribution as degradation then dispersal. Quantify the noncommutativity only for biologically specified operators; do not call a generic nonzero matrix commutator a biological discovery.

## Simulation experiment (preregister before empirical fitting)
Simulate replicated dated trees and geographic colonizations under M0–M3, crossing competence degradation, migration, introgression, sampling, observation error and independent clade counts. Evaluate (a) mechanism classification confusion matrix; (b) calibration of A and G; (c) held-out log score; (d) minimum molecular samples required; (e) robustness to tree and taxonomy uncertainty. Negative controls: no competence decay, no gene flow, plastic colour shifts, and pure recolonization without local re-expression.

## Empirical claim ceilings
Photo records: visible colour and geography only.
Dated tree + colour: probabilistic historical transitions conditional on model.
Same-individual floral chemistry/expression: pathway activity, not necessarily heritable competence.
Common garden, crosses, genetic perturbation and ancestry-aware haplotypes: stronger evidence of reacquisition mechanism.
No transition date or quantitative evolutionary window without validated dated history.
