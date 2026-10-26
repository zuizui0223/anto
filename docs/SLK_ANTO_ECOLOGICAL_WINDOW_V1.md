# SLK × ANTO: when an ecological opportunity outlives evolutionary competence

**19 October 2026. Theory-only cross-repository bridge v1.** This is a mathematically exact **constructed example** using an existing SLK comparison and an existing ANTO hidden-resource/phenotype transition kernel. It is neither a new mathematical theorem about Markov chains nor a fit to plant data. It does not estimate re-evolution, fitness, establishment or fixation in any flowering plant.

## Biological question and falsifiable candidate hypothesis

**Does the timing of ecological conditions favorable for renewed functional differentiation overlap with the period when a lineage can still express the required developmental architecture?**

An advantageous phenotype can be absent either because its ecological return is not yet favorable (SLK), or because the necessary joint functional modules are not available or not expressed (ANTO). Their interaction adds a third outcome beyond each repository alone: *a favorable ecological window may pass without any eligible repeatedly expressed phenotype*, or phenotypes can briefly reappear after a window but never overlap the favorable conditions.

This is **not** the assertion that two mathematically specified probabilities can be multiplied. It is a pathwise overlap experiment in which the same phenotype trajectory is crossed with distinct environmental histories, and a three-snapshot event is recorded only when **all** qualifying ecology signs hold during all three snapshots.

## Exact interface and its limitations

### SLK side (registered quadratic witness, not arbitrary biology)

Consider only a **candidate paired comparison**: integrated shared-function architecture S versus a specified differentiated architecture D, with a proposed development/re-expression route. The existing SLK witness uses degree x in [0,1] and

```text
R(x) = x + x^2
K(x,E) = k(E) x
Phi(E) = R(1)-K(1,E) = 2-k(E)
g0(E) = R'(0)-k(E) = 1-k(E)
Delta_R(E) = Phi(E)-eta(E)
```

Under the registered canonical positive-frequency feedback `Delta(p,E)=Phi(E)+eta(E)(2p-1)`, `Delta_R` is selection against an S-resident **at infinitesimal D frequency**. The three nested *strict* conditions are:

```text
VALUE: Phi>0
PATH: Phi>0 and g0>0
ALL:  Phi>0 and g0>0 and Delta_R>0
```

This last condition **excludes three early barriers**, but does not establish a globally feasible developmental/mutational path, absolute population growth, heritable phenotype, fixation or persistence. It is a deliberately conservative gate. Even `g0>0` is meaningful only when the chosen differentiated form is related to the specified SLK local release path. **A flower-colour return by itself is not SLK functional division of labour.**

### ANTO side (actual existing transition kernel, not an ecological population model)

Reuse `Rates(activation, repeat, degradation, acquisition)` and `transitions` in `src/anto_afterlife_frontier_v10.py`; no redefinition of ANTO's biological states:

- `theta` is the initial probability that **all necessary joint functional modules** exist within a lineage, not mean gene presence.
- `a` is expression activation, `r` is repeat-expression probability, `d` is irreversible toolkit-loss probability for that transition, and `g` abstracts reacquisition/reassembly.
- `g>0` does **not** specifically identify gene flow, native pathway reactivation, de novo mutation or viable seed establishment.
- Every lineage begins visibly negative; none is an already visibly positive immigrant.

The bridge leaves expression and resource state transitions **unchanged** when ecology changes. It marks whether a particular expressed snapshot is within an SLK-eligible generation. The target event is three **consecutive** expressed snapshots, all passing the selected gate, at least once within 14 model generations. Once this *first-event occurrence* is recorded, it remains counted even if resources are later lost. Therefore it is emphatically **not a probability of demographic establishment, heritable restoration, or long-run adaptation**.

## Exact illustrative results

Fixed settings: `theta=1/8`, `a=2/5`, `r=3/4`, `d=1/10`, `g=0`, horizon 14 generations.

### Result 1 — equal total favorable time does not imply equal opportunity

Both ecological schedules have exactly **five** all-gate-positive generations:

| Equal-duration schedule | Strictly favorable generations | Qualified three-expression event by generation 14 |
|---|---|---:|
| Favorable early | 1–5 | **5.090625%** |
| Favorable late | 10–14 | **2.291109%** |

The probability of **ever expressing at least once** (irrespective of ecology) is **identical** for the two schedules, because the ANTO transition kernel is identical. So is the raw ANTO three-expression event. The 2.22-fold difference is caused specifically by the temporal placement of five eligible SLK generations against the fixed ANTO expression/loss dynamics. An ecological window late enough to miss much retained competence can contain the same total favorable time and still have much lower qualified-event probability.

This is a constructive contrast **under the declared toy rates and starting state**, not a universal monotone effect of late colonization or a measured effect of climate change.

### Result 1b — even the direction of the timing effect can reverse

It would be misleading to generalize Result 1 into "earlier favorable conditions always give more evolutionary opportunity." With *the same loss rate* `d=1/20`, the *same initial competence* `theta=1/8`, and the *same two schedules of five favorable generations*, changing expression activation and continuation reverses the ranking:

| Expression regime, both with d=0.05 | Activation a | Repetition r | Favorable early | Favorable late | Direction |
|---|---:|---:|---:|---:|---|
| Faster activation, moderate continuation | 0.40 | 0.75 | **5.435156%** | 3.901086% | Early > late |
| Slow activation, high continuation | 0.05 | 0.95 | 1.528609% | **3.194658%** | Late > early |

In the second constructed regime, the unfavorable pre-window generations can allow latent lineages to *begin expressing and continue expressing* before the ecological opportunity opens. With slow activation, an early five-generation window is often missed; late-arriving favorable conditions can coincide with an already-expressing pool even though a positive resource-loss hazard operates throughout.

This **temporal preconditioning effect** is model-specific, not evidence that plants can freely sustain a selected-against phenotype before a future ecological change. If expression itself is costly under unfavorable ecology, or if expression/decay rates depend strongly on the environment, the ranking may change. Crucially, merely knowing the duration of the future selective window and the loss hazard is not enough to rank histories without activation and continuation.

**Biological implication to test:** ancestral latent trait retention can create either historical inertia or delayed re-expression opportunity; neither earlier colonization nor faster first expression is universally the better predictor of repeated functional expression. Timing, hidden resource inventory and expression kinetics jointly matter.

### Result 2 — selective profitability, access and invasion open at different times

On a monotonically decreasing cost example `k(t)=max(0,12/5-t/5)`, `eta=3/2` (cost saturates at zero in generations 12–14), the gates become strictly positive at different generations:

| Qualifying condition | First positive generation | Eligible generations | Qualified three-expression event by gen 14 |
|---|---:|---:|---:|
| Positive architecture value only | 3 | 12 | **5.983948%** |
| Value and positive local release | 8 | 7 | **3.262062%** |
| Value, path and positive rare-invasion margin | 10 | 5 | **2.291109%** |

The ordering is exact because the events are nested on identical ANTO sample paths. The particular calendar thresholds are chosen numbers, **not estimates of ecological time-to-regain**. In particular, `Phi>0` at generation 3 does not imply that a three-generation phenotype event could establish there.

### Result 3 — genomic summary statistics may not identify the window

For three hypothetical functional modules, `P(M_i=1)=1/2` and `P(M_i=M_j=1)=1/4` for every pair permit `theta=P(111)` anywhere from **0 to 1/4**. Under the late schedule with **no reacquisition**, the exact conditional qualified-event interval is **[0, 4.582218%]**. This interval is sharp for the specified two-initial-state ANTO closure and equal exogenous transition rates across the unknown toolkit states; it is **not a confidence interval**. With reacquisition, a positive event can occur even from `theta=0`, but that route must not be called ancestral endogenous reuse.

## What this adds, and what it does not

**Positive model deductions:** A stationary visible phenotype can conceal changing opportunities for reorganization, *and* an ever-expressed phenotype can conceal different overlap with eligible ecological states. The joint timing of genetic competence and ecological selection is a distinguishable constraint, not reducible to current phenotype, initial gene counts, or the total duration of favorable conditions. The resource-state effect and the SLK invasion effect remain causally distinct.

**Not demonstrated:** common real-world loss rates, a globally unimodal re-evolution window, any observed independent Cirsium regain, an experimental SLK cost `k` or frequency feedback `eta`, or positive effects on reproductive fitness. SLK gate positivity does not itself supply an establishment probability; ANTO's three-snapshot persistence proxy is not ecological fixation. Phenotypic lack of return might reflect no capacity, insufficient induction, ecological exclusion, source sorting, or observation error.

**Falsifiable biological extension:** In genuinely independently replicated lineages with known histories, separately quantify (i) functional joint components in the same individual, (ii) stage-matched induction and persistence, (iii) fitness and rare-variant performance under the relevant environment, and (iv) the timing/order of ecological change. The key outcome is whether **history-by-environment timing predicts functional reorganization beyond contemporary ecology and visible phenotype**, after controls for ancestry, introgression, founder filtering and retained resources. Replicated flowers from the same plant are not independent evolutionary regain events.

The *Cirsium* practical phenotype/genomic/fitness pipeline belongs exclusively to [aza3](https://github.com/zuizui0223/aza3). Floral colour (ANTO) can support the competence question there, but it is not automatically the multifunctional-architecture conflict required by SLK. A valid joint organismal application needs an empirically established SLK shared-function conflict as well as an ANTO-style retained/degraded **same-architecture** developmental resource.

## Repository contracts and reproducibility

- ANTO general-model owner: `src/anto_slk_ecological_window_v1.py` (imports **unchanged** ANTO v10 transition kernel and v7 exact feasibility bounds).
- Executable tests: `tests/test_slk_anto_ecological_window_v1.py`.
- Exact JSON model receipt: `scripts/run_slk_anto_ecological_window_v1.py`.
- GitHub CI: `.github/workflows/slk-anto-ecological-window-v1.yml`.
- SLK theoretical reference: [SLK unified threshold atlas](https://github.com/zuizui0223/slk/blob/biology/scaffold-bypass-test-20261014/theory/UNIFIED_THRESHOLD_ATLAS_V1.md) on its still-unmerged biology branch; **no assertion** that the current SLK main already contains this merged cross-repository code.
- Primary existing ANTO references: `docs/EVOLUTIONARY_AFTERLIFE_V10_EXACT_FRONTIER_AND_UNIVERSALITY.md`; `docs/JOINT_MODULE_ACCESSIBILITY_AND_DISPERSAL_V7.md`; `docs/AFTERLIFE_V10B_INDUCTION_IDENTIFIABILITY_AND_EMPIRICAL_GATES.md`.

Run: `python -m unittest discover -s tests -p test_slk_anto_ecological_window_v1.py -v`; `python scripts/run_slk_anto_ecological_window_v1.py`.

**Manuscript decision:** treat this as a cross-repository *conditional ecological prediction* and potential joint falsifier, not a novel universal mathematical theorem or evidence of adaptive trait resurrection. Keep the SLK and ANTO flagship questions independent until an organismal dataset demonstrates a genuinely shared comparison.
