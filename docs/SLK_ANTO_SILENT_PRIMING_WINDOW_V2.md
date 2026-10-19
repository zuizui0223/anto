# SLK × ANTO v2 — Invisible developmental priming can create a transient evolutionary opportunity

**19 October 2026. Exact conditional synthetic mechanism; no measured plant priming, lineage return, fitness effect or genetic restoration.** This extends [v1](SLK_ANTO_ECOLOGICAL_WINDOW_V1.md) to address an explicit weakness: v1's late-window advantage sometimes required **visible trait expression before the SLK ecological gate became favorable**. Positive selection for such expression before the favorable period was not established. Here NO trait is visibly expressed before the favorable three-generation interval; only a *hypothetical hidden regulatory-readiness state* changes.

## Central biological question

Can ecologically neutral/nonqualifying history silently prepare a retained but unexpressed trait for later reactivation, so that a delayed favorable selection window outperforms an immediate one *despite ongoing irreversible toolkit loss*?

**Testable answer in this declared model:** yes, provided silent regulatory priming outruns resource erosion. The effect reverses when priming is absent, activation gains are absent, or resource degradation wins. This is not a general law of phenotypic resurrection, nor a new probability theorem.

## Exact mechanistic states (same organism / retained ancestral toolkit)

- **U (unprimed):** fully joint-competent functional resource inventory is retained, expression off and not regulatory primed; phenotype absent.
- **P (silently primed):** SAME complete inventory, expression still OFF, but activation after a later relevant cue is more probable.
- **D (resource missing):** required joint toolkit has been irreversibly lost in the closed model.
- **X (expressed):** visible/functional trait has appeared under the later favorable environment; continuing X for two further snapshots counts an *observational three-expression streak*, not adaptive fixation.

Unlike introgression or de novo pathway reconstruction, U→P **does not add missing coding function** and does not cause a visible phenotype. Both U and P would be white for a colour example. The existence of an unexpressed yet primable molecular state is a **research hypothesis**: its genetics, epigenetics, ontogenetic timing, measurement and heritability are not validated in Cirsium.

## Coupling to SLK is explicit but deliberately narrow

Use the same candidate-specific registered SLK architecture witness from [SLK draft #117](https://github.com/zuizui0223/slk/pull/117) and the [v1 bridge](SLK_ANTO_ECOLOGICAL_WINDOW_V1.md):

```text
R(x)=x+x²; K(x,E)=k(E)x; Phi=2-k; g0=1-k; Delta_R=2-k-eta.

Ecologically NONQUALIFYING history: the same candidate fails at least one
  of {Phi>0, g0>0, Delta_R>0}; no visible ANTO trait is produced then.

Ecologically QUALIFYING opportunity: all three strict signs are positive
  in each of EXACTLY three subsequent generations.
```

The gate does NOT determine the genetic state, supply of structural mutants, demographic growth, selection on hidden readiness, field pollen reward/transfer or long-run persistence. Candidate-specific adaptive benefit is a condition for opportunity, not a numerical establishment probability. A pigmented flower alone is **not** necessarily SLK functional division of labour. To apply both frameworks to one organism, the latent component must be identified as part of the **same multifunctional-to-differentiated architectural comparison**.

## Model transitions

Start with initial probability `theta` of being completely competent but unprimed U, and `1-theta` resource missing D.

During each of `m` nonqualifying generations (no visible expression at all):

```text
U -> P with b                 [silent regulatory priming]
U -> D with d                 [irreversible resource loss]
U -> U with 1-b-d
P -> D with d
P -> P with 1-d
D -> D with 1
```

This is a **different explicitly assumed ecological-transition kernel**, not just the v1 gating mask. During the subsequent three qualifying generations, U activates with probability `a_U`, P with probability `a_P`, and each successful expressed state continues with probability `r` into the next qualifying generation. The probabilities `a_U+d`, `a_P+d`, `r+d` and `b+d` are each at most one, preserving competing-event validity.

Neither ANTO v10's acquisition route nor immigration is present. No latent-priming state exists in the original v10 kernel; this extension is intentionally separate rather than silently changing that model.

### Exact result (not a fitted biological formula)

Let `S=1-d`, `T=1-b-d`. After the nonqualifying history:

```text
P(U_m)=theta*T^m
P(P_m)=theta*(S^m-T^m)
P(D_m)=1-theta*S^m
P(visible positive during bad history)=0
```

The opportunity event is **three consecutive positive snapshots within the only three SLK-qualifying generations**. Thus its exact probability is:

```text
Q(m)=theta*r^2 * [a_U*T^m + a_P*(S^m-T^m)].
```

This is a conditional phenotype *event*, not measured natural adaptation. An independent iterative U/P/D Markov oracle checks the closed form.

**One-step delay threshold:**

```text
Q(1)-Q(0)=theta*r^2*[b*(a_P-a_U)-a_U*d].
```

For positive initial toolkit probability and continuation, a one-generation delay strictly **improves** the opportunity if and only if

```text
b*(a_P-a_U)>a_U*d.
```

Here `b*(a_P-a_U)` is the new *silent readiness premium*, while `a_U*d` is the baseline activation probability lost through genetic toolkit erosion. This is a falsifiable **parameter inequality**, not a universal ecological law.

**Finite interior maximum:** If `d>0`, `b>0`, `a_P>a_U>0` and the one-step inequality is strict, then `Q(m)` initially increases and eventually tends to zero. In fact:

```text
Q(m+1)-Q(m)=theta*r²*[(a_P-a_U)*(b+d)*T^m-a_P*d*S^m].
```

Since `0<=T<S<1`, the increment changes sign at most once from positive to negative. Hence the integer delay profile has a single peak (or at most an adjacent tie), not a permanently improving late-arrival rule. If `T>0`, the continuous increment-crossing is at

```text
m_c=log(a_P*d/[(a_P-a_U)*(b+d)]) / log(T/S),
```

with the first discrete optimum at `ceil(m_c)` when the crossing is positive and interior; equality permits adjacent tied delays. If no priming (`b=0`) or priming does not improve activation (`a_P=a_U`), then `Q(m)=theta*r²*a_U*(1-d)^m` and delay cannot help. If `d=0`, no finite closure is guaranteed; priming can accumulate indefinitely in this toy setting.

## Exact illustrative results and an independently testable control

Set `theta=1/8`, `a_U=0.10`, `a_P=0.60`, `b=0.15`, `d=0.04`, `r=0.80`.

| Wait before the 3 favorable generations | Qualified three-snapshot probability |
|---|---:|
| 0 | **0.8000%** |
| 1 | 1.3680% |
| 2 | 1.7993% |
| 6 | 2.6275% |
| 8 | 2.7215% |
| **9 (peak)** | **2.7238%** |
| 10 | 2.7049% |
| 15 | 2.4325% |
| 25 | 1.7093% |

The first-delay boundary equals `0.071>0`, and the best wait in the registered 0–50 horizon is exactly **9 generations**. The optimum produces approximately **3.40 times** as many eligible three-snapshot events as an immediate favorable window; these are *synthetic* percentages.

**More important same-present-state diagnostic:** After nine unfavorable generations, both a `b=0.15` priming history and a `b=0` nonpriming history have:

- exactly the **same complete functional toolkit fraction** `theta*(1-d)^9=8.656675%`;
- exactly **zero visibly expressed trait**;
- the **same future SLK environment** and same potential `a_U,a_P,r` values;
- different hidden shares of U and P, and therefore qualified three-expression probabilities of **2.723785% versus 0.554027%** (about **4.92-fold**).

This does **not** violate present-state sufficiency: it says a **coding/inventory-only assay** is an incomplete present-state assay. A truly complete assay including the distinct molecular priming state would distinguish the cases. The biological question is whether there exists a measurable flower-specific latent readiness variable that carries some information beyond sequence integrity and visible morphology.

## Falsifiers and organismal research programme

In replicated whole-lineage cohorts with identified ancestral capacity (not photos treated as independent origin events):

1. **Verify silent status** throughout the unfavorable phase by pigment, organ-appropriate expression, and developmental-stage assays. A visible precursor is not a silent-priming demonstration.
2. **Measure competence and priming separately** on the SAME individuals or tracked genetically defined lines, including complete joint functional modules and independently chosen regulatory/expression markers. No invented transcript "priming" labels.
3. **Pre-expose versus immediate-induce**, holding future favorable induction conditions identical. Quantify activation after `m` pre-exposure generations and measure resource damage/loss over time. Include `b=0`-type cue controls, alternative histories, and matched population ancestry/cytotype.
4. **Distinguish physiological memory from evolution.** A transient within-plant priming effect without inherited effect is evidence for developmental/environmental memory, not ancestral genomic reacquisition. Measure multigeneration offspring and heritable allele/regulation differences if claiming evolution.
5. **Close the SLK contrast independently.** Test the same focal structure's opposing functions, its matched full reproductive fitness under integrated versus differentiated architecture, initial selective release and rare-variant establishment. Do not substitute pigment presence, pollen deposition or flower colour for SLK `Phi`.
6. **Include negative lineages** whose functional toolkit is demonstrably absent, and donor-gene-flow controls. Recovery after introgression is a separate path, not latent native regulation.

The reported *Cirsium japonicum* multi-organ contrasts in [ANTO v6](CIRSIUM_ORGAN_RETENTION_AND_HIDDEN_COMPETENCE_V6.md) motivate possible nonfloral conservation, **not** evidence of floral priming or the parameter `b`. All empirical Cirsium field, molecular, head LEGO and reproductive-fitness work remains with [aza3](https://github.com/zuizui0223/aza3).

## Code, checks and publication scope

- Exact model: `src/anto_slk_silent_priming_window_v2.py`.
- Independent hidden-state iterative oracle and negative controls: `tests/test_slk_anto_silent_priming_window_v2.py`.
- Machine-readable fraction-accurate receipt: `scripts/run_slk_anto_silent_priming_v2.py`.
- CI: `.github/workflows/slk-anto-silent-priming-v2.yml`.
- Complementary v1: `docs/SLK_ANTO_ECOLOGICAL_WINDOW_V1.md`.
- SLK companion: [draft PR #119](https://github.com/zuizui0223/slk/pull/119).

The point of this interface is **not** to suggest that "reversal is possible" or "selection plus hidden history matters" are new established ideas. The scientific proposal is the **specific hidden-readiness versus irreversible-resource-loss test** and its conditional single-peaked timing signature. Without real same-system priming, ancestry and SLK fitness data, it remains a model-based falsifier and must not be sold as an experimentally supported evolutionary rescue law.
