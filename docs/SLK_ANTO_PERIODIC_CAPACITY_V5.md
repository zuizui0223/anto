# SLK x ANTO v5: ecological run length can determine whether recoverability persists

2026-10-20. **Conditional theoretical result; no plant observations or fitted rates.**
This extends the unchanged v4 two-state offspring kernel to a known periodic
sequence of environments. It does not add an unobserved priming mechanism,
perfect expression inheritance, a finite census target, or a positive SLK
payoff-to-fecundity mapping.

## Biological question

V4's illustrative reversal assumed that the nonexpressing competent state C
reproduced better than the expressing state D in a constant environment.
Does retention still matter when neither state is always favored, and the
same two environments recur with the same frequencies?

**Model result:** yes. The persistence boundary of the recoverable genealogy
can change when only the consecutive duration of the environments changes.
The direction is not universal: it reverses with the probability of offspring
switching their expression state. This is an explicit application to ongoing
functional-toolkit loss, not the first theory of phenotypic switching.

## What is held fixed

D expresses a focal candidate structure. C reproduces without expressing it
but retains all functionally necessary resources. Each parent has two
independent reproductive slots, exactly as in `anto_slk_recurrent_capacity_v4.py`.
In a slot a child is produced with probability q_D or q_C; it retains the toolkit
with probability 1-d. Conditional on retention, an expressing parent's child
silences with probability s, and a silent parent's child expresses with
probability a. Toolkit-lost children leave the modeled genealogy: they need
not die as organisms. There is no gene flow or reconstruction after loss.

The illustrative two environments reverse the reproductive advantage:

| Environment | q_D | q_C | Which expression state produces more offspring? |
|---|---:|---:|---|
| A | 0.90 | 0.40 | D |
| B | 0.40 | 0.90 | C |

Both environments occupy exactly half of all generations. The independent
per-child toolkit loss is d=0.20 throughout. These rates are constructed,
not natural mutation rates or measured flower fitness. In particular, the
20% value must not be described as a measured genetic degradation rate.
The initial state is one C parent immediately before an A generation.

Compare a fast environmental sequence `ABABABAB` to `AAAABBBB`, repeated
indefinitely. Both have four A and four B generations per eight-generation
supercycle. The only difference WITHIN each pair is environmental run length.
Period 2 is the minimal representation of the first sequence; repeating it
four times leaves its long-run boundary and starting-phase survival unchanged.

## Positive result and counterexample

| Offspring switching a=s | Ecological schedule | Growth per generation of competent genealogy | Indefinite competent-genealogy survival from the same initial C |
|---|---|---:|---:|
| 0.10 | ABABABAB | 0.968925 | **0%** |
| 0.10 | AAAABBBB | 1.016400 | **2.08593%** |
| 0.90 | ABABABAB | 1.308580 | **10.35677%** |
| 0.90 | AAAABBBB | 0.977826 | **0%** |

The zero probabilities refer to ultimate loss of all competent descendants
under these permanently repeated, closed environments. They do NOT imply
extinction of the organisms or fixation of a nonexpressing phenotype.
Nonzero survival does not mean every lineage survives or that every
individual expresses the trait. Founder phase and starting type matter to
these probabilities, even when the survival boundary is unchanged.

For a=s=0.10, the two schedules also have the SAME next-generation chance
of at least one expressed child from initial C: **6.2976%**. Therefore the
contrast is not an artifact of unequal first-induction success. Initial
capacity, first environment, and the one-step offspring distribution agree.

The low-switching run-length sweep is 0.968925, 0.978386, 0.994363 and
1.016400 for runs of length 1, 2, 3 and 4, respectively. This is a numerical
boundary in a declared parameter slice, NOT a four-generation rule for plants.
At high switching, run-length responses need not be monotonic: length 3 is
supercritical even though lengths 2 and 4 are subcritical in this example.
Perfectly periodic resonance must not be generalized to random environments.

## Why averaging the environment gives the wrong answer here

With parent rows and offspring columns, the one-generation mean matrix is

```
M_e=(1-d)*[[2*q_D,e*(1-s), 2*q_D,e*s],
           [2*q_C,e*a,     2*q_C,e*(1-a)]].
```

For a period e_0,...,e_(L-1), the mean cycle is the CHRONOLOGICAL product
`P=M_e0 M_e1 ... M_e(L-1)`. Its normalized Perron root is
`g=rho(P)^(1/L)`. In this strictly positive, nonsingular finite-type process,
indefinite competent-genealogy survival is possible iff `rho(P)>1`.

The arithmetic average matrix is identical for the two schedules. In the
symmetric example it has growth 1.04 for both switching probabilities.
That average would predict growth even for the schedules with actual g<1.
The mismatch arises because descendants have a history-dependent distribution
of expression states when the environment changes. Multiplying the phase
matrices preserves that distribution; averaging them does not.

A permanently D lineage with switching disabled and a permanently C lineage
with switching disabled both have geometric growth
`2*(1-d)*sqrt(q_high*q_low)=0.96` under equal environmental exposure.
Neither can retain capacity indefinitely on its own in this construction.
The coupled reproducing genealogy can: it traverses expressing and silently
competent states. C is not a dormant seed or a permanently nonreproductive
reservoir.

## A sharp ecological effect on the loss threshold

When the same per-child loss d applies in EVERY phase, it factors from each
matrix. If P_0 is the period matrix at zero loss, then

```
P(d)=(1-d)^L P_0,
d_critical=1-rho(P_0)^(-1/L).
```

The boundary is signed; if d_critical<=0, even zero loss cannot produce a
strictly supercritical process. For switching 0.10, the positive threshold
is **0.174343** for ABABABAB and **0.212908** for AAAABBBB. Thus any d strictly
between those numbers gives opposite classifications, an open parameter
interval rather than a random-seed accident. At the illustrated d=0.20,
only the block schedule permits competent-genealogy survival.

This threshold concerns functional competence retained through reproducing
lineages, not the lifetime of a gene in one individual or the proportion
of white flowers. Positive initial competence frequency affects the chance
of founding a surviving genealogy, not this periodic critical surface.

## Controls, interpretation, and connection to SLK

Rotating the starting phase of a cycle preserves its Perron root. In
particular, `AB` versus `BA` cannot be advertised as a different long-run
threshold: those are phase shifts of the same periodic regime. They can
have different survival probabilities from a single specified founder.
The tests check both the invariant boundary and the correct phase-linked PGFs.

At a=s=0.50 in the symmetric example, the offspring's expression assignment
is independent of the parent's expression state. The run-length dependence
vanishes (g=1.04). This is a specific no-memory control, not a general claim
that histories cannot matter after any noisy assay.

SLK contributes the need to compare the functional returns of ALL viable
architectures under the relevant ecological partner context. V5 makes D and
C alternative reproducing states whose ranking actually reverses between
contexts; it does not assume D beats C because it beats some original S.
ANTO contributes inherited competence, irreversible loss and recurrent
expression. The remaining same-organism bridge is substantive: the retained
resource must belong to the SAME multifunctional structure whose division or
reorganization SLK compares. Flower recolouring alone does not establish that.

A biological prediction worth testing is that two populations with similar
current phenotypes, environmental means and retained pathways can differ in
future functional re-expression because their ecological regimes have
different durations relative to expression-state inheritance. Paired exposure
histories would need matched ancestry, measured full functional resources,
stage-specific expression, and absolute reproductive output. No experiment
or macroecological validation of this full prediction was performed here.

## Prior art and scientific ceiling

Acar, Mettetal & van Oudenaarden (2008), *Nature Genetics* 40:471-475,
DOI 10.1038/ng.110, experimentally demonstrated in engineered yeast that
fast versus slow phenotype switching has environment-frequency-dependent
growth consequences. The broad switching/time-scale principle is prior art.
Hufton, Lin & Galla (2018), DOI 10.1088/1742-5468/aaa78e,
arXiv:1706.07789, explicitly distinguish periodic and stochastic environments.
Carja & Plotkin (2019), DOI 10.1534/genetics.118.301758, examine partly
heritable phenotypic variability and evolutionary rescue. These are not
new independent observations of floral trait regain by this project.

The new repository contribution is the explicit competence-loss boundary,
matched environmental compositions, invariant/counterexample checks and
phase-specific long-run survival in the existing ANTO resource process.
No new general branching theorem, universal optimal switching strategy,
measured plant loss rate, genetic fixation, population rescue or natural
re-evolution is claimed. The environment is prescribed, periodic and common
to all descendants; random environmental noise, density dependence and
feedback of trait frequencies on partners remain outside this model.
All Cirsium empirical work stays in aza3. Neither draft PR is merged, and
SLK's main submission manuscript is not changed by this extension.

## Reproduction

- Model: `src/anto_slk_periodic_capacity_v5.py` (imports unchanged v4).
- Tests: `tests/test_slk_anto_periodic_capacity_v5.py` (15 focused tests).
- Runner: `scripts/run_slk_anto_periodic_capacity_v5.py`.
- Selected results: `results/slk_anto_periodic_capacity_v5_summary.json`.
- Independent grid: `results/slk_anto_periodic_capacity_v5_validation.json`.
- CI: `.github/workflows/slk-anto-periodic-capacity-v5.yml`.

The cycle mean/classification use rational arithmetic. Survival uses
70-digit Decimal PGF composition with converged lower/upper iterations;
these are numerical brackets, not outward-rounded interval certificates.
Independent forward offspring-count enumeration checks PGF order and means.
364 independent NumPy product/eigenvalue comparisons gave maximum normalized
growth error 1.12e-15; 357 non-boundary classifications agreed, and seven
numerically near-critical grid points were not assigned a class by floating
point. An exact critical unit test remains separate. Local replay ran all
29 tests present in the v4+v5 source bundle. Full repository verification is
performed separately in GitHub CI and must be reported from its actual run.
