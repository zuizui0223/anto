# ANTO × SLK v3 — regained expression is not evolutionary rescue

**19 October 2026. Theory-only conditional demographic comparison. No new plant observations, calibrated natural loss/regain rates, heritable recovered thistle, or evidence of population rescue.** This extends [v2 silent competence](SLK_ANTO_SILENT_PRIMING_WINDOW_V2.md) by separating the *supply of a re-expressed phenotype* from whether a resulting rare phenotype lineage can produce descendants under an ecological fitness regime. It does not replace SLK's population game or anto's molecular evidence.

## Biological question

If an ancestral trait's developmental machinery survives its visibly silent period, can renewed selection after an environmental reversal actually produce a lasting biological response? Why can the *same apparent advantage when rare* still fail to restore a population?

**Central model-conditional finding:** The same latent toolkit frequencies, preparation history, candidate-specific positive SLK signs, and **identical expected relative offspring ratios** can give sharply different numbers of descendant-bearing lineages, because a selection advantage is not an absolute demographic growth rate. The finite chance of leaving descendants and the genuinely long-run nonextinction probability are themselves different endpoints.

This is the classic distinction between selection and demographic viability, applied to the **specific conjunction of ANTO hidden capacity and SLK candidate-architecture thresholds**. Prior art on branching evolutionary rescue and population-size/demographic feedback is acknowledged: [A branching process model of evolutionary rescue](https://doi.org/10.1016/j.mbs.2021.108708), [Demographic feedbacks during evolutionary rescue can slow or speed adaptive evolution](https://pmc.ncbi.nlm.nih.gov/articles/PMC10865011/). The branching PGF mathematics and absolute/relative fitness distinction are **not** claimed as new.

## Separate biological stages: what is actually modeled

1. **ANO/ANTO hidden preparation**: before the good environment, a lineage remains phenotypically negative, while retained joint functional material transitions between unprimed U, primed P, and irreplaceably missing D, exactly as in ANTO v2. Initial complete-toolkit probability is theta; the negative period lasts m generations.
2. **Specified SLK advantage**: the ecological state after the switch has `Phi>0`, `g0>0`, and `Delta_R>0` for one *explicit differentiated architecture*. The numerical example is the registered witness `R(x)=x+x²`, `K(x)=kx`, `k=0.4`, `eta=1.5`, so `Phi=1.6`, `g0=0.6`, `Delta_R=0.1`. These are hypothetical payoff units. They are NOT a measured fitness advantage of a pigmented flower.
3. **An explicit NEW heritability assumption**: each founder who first activates when the environment switches produces one individual of a phenotypically differentiated D lineage that **transmits its differentiated state to descendants**. This step is *not justified by the earlier ANTO three-expression observation model*; a one-off induced pigment change need not be inherited.
4. **Absolute demographic outcome**: in each good generation, each D descendant has two independent reproductive opportunities, each successful with probability `q_D`; hence one descendant's offspring count is `Binomial(2,q_D)`. For comparison the unmodeled resident S has `Binomial(2,q_S)` offspring. The condition `q_D>q_S` is an additional, explicitly declared demographic mapping of SLK's positive rare-selection direction, **not automatically its numerical payoff-to-fecundity conversion**. The generation-by-generation S–D competition, population density, introgression, and feedback of the new lineage on pollinators are not modeled.
5. **Prespecified finite census**: at generation 14 since initial trait absence, ask whether **at least one out of 40 independent hypothetical initial lineages** has at least **four D descendants**. The 40 lineages are not 40 flower images or 40 within-species plants from aza3. This 4-descendant threshold is a transparent *transient census criterion*, NOT persistence, fixation, population rescue, viable flower-colour evolution or lifetime fitness.

## Mathematical connection, with its identifying boundaries

Let the exact ANTO v2 activation probability after m unfavorable generations be

```text
A(m) = theta*[a_U*(1-b-d)^m
               + a_P*((1-d)^m-(1-b-d)^m)].
```

This is the fraction of initially sampled lineages able to generate one D founder at the switch, **conditional on experimental heritable re-expression being possible**. Unlike v2's three-consecutive-expression event, `A(m)` has no factor `r²`.

The independent Galton–Watson offspring PGF is `f_q(z)=(1-q+qz)^2`. For a single D founder and T good reproductive generations, `G_T(z)=f_q^{∘T}(z)`. With a target of at least `k` descendants,

```text
B_q(T,k) = P(N_T>=k | N_0=1)
         = 1 - sum_{j=0}^{k-1} [z^j] G_T(z).
```

The exact polynomial recurrence is implemented with 65-digit decimal arithmetic (finite numerical precision) and checked against **independently enumerated full offspring-count distributions** for small T. No approximation by `q_D-q_S` or `Delta_R` is made when computing `B_q`.

For n independent initial propagules, the **finite descendant-census event** is

```text
C(m;H,k,n) = 1 - {1-A(m)*B_q(H-m,k)}^n.
```

Because C is monotone in the single-founder joint event, the exact discrete decision between waiting m versus m+1 (when denominators are positive) is

```text
C(m+1)>C(m)
iff
A(m+1)/A(m) > B_q(H-m,k)/B_q(H-m-1,k).
```

**Biological interpretation:** The relative gain in preparedness from one more generation of history must exceed the cost of **one fewer reproductive generation to build an established-looking descendant lineage**. Even if `A(m)` peaks late, the maximum `C(m)` can occur earlier. This "two competing clocks" statement is conditional on the declared model, not a new theorem about branching processes.

The actual **infinite-environment** probability that an activated lineage avoids eventual extinction in this Binomial(2,q) model is

```text
E(q) = 0                        if q <= 1/2
       1 - [(1-q)/q]^2          if q > 1/2.
```

This is the classical Galton–Watson nonextinction result and assumes a permanently stable demographic environment. It is categorically NOT the finite-census `B_q(T,k)`, or genetic fixation in a bounded population.

## Constructed contrast: same relative selection, dramatically different demographic fate

Hold `theta=1/8`, `a_U=0.1`, `a_P=0.6`, `b=0.15`, `d=0.04` fixed. The priming and activation probabilities are IDENTICAL in both worlds, including the same all-white pre-exposure phenotype.

The ecological SLK `Phi`, local slope, and rare margin are held fixed, while the new absolute demographic mapping is changed:

| | Absolute-supercritical world | Absolute-subcritical world |
|---|---:|---:|
| S resident offspring success `q_S` | 0.55 | 0.33 |
| D candidate offspring success `q_D` | 0.65 | 0.39 |
| D/S expected offspring ratio | **13/11** | **13/11** |
| Candidate absolute mean offspring | **1.30** | **0.78** |
| Infinite-horizon D nonextinction, if ecology stays favorable | **120/169 = 71.01%** | **0%** |
| Best finite waiting period before 14-gen census | **6 generations** | **9 generations** |
| Maximum 40-founder finite census event | **64.68%** | **1.62%** |

The finite-census percentages are not population persistence probabilities. Particularly in the **subcritical** world, observing four D descendants temporarily is compatible with **certain ultimate extinction** in the stipulated persistent environment.

### Fixed-delay comparison (more controlled than comparing each world's best)

| Unfavorable pre-exposure (m) | D activation fraction `A(m)` shared by both | Supercritical 40-founder census | Subcritical 40-founder census |
|---|---:|---:|---:|
| 0 | 1.250% | 29.64% | 0.082% |
| 6 | 4.105% | **64.68%** | 1.022% |
| 9 | 4.256% | 55.08% | **1.615%** |
| 12 | 4.097% | 11.64% | 0.575% |

**The more informative result is the timing disagreement:** the underlying genetic/developmental founder supply `A(m)` peaks at 9 in the example. Yet in a supercritical environment, the fixed generation-14 threshold for four descendants peaks **earlier, at 6**, because later activation leaves too little time for expansion. Under the subcritical process, the finite 4-descendant transient can peak at 9, but eventual persistence remains zero.

*Neither the precise 6/9 generations nor the numerical percentages are estimates of any plant species.* The values require a one-founder activation event, inherited D, independent whole lineages, no density dependence or partner turnover, and a fixed finite census. Conditioning on the existence of flowering and viable offspring is NOT resolved by these numbers.

## What is actually nontrivial for organismal evolutionary ecology?

1. **Regained expression is not re-evolutionary success.** An induced floral pigment and the genetic capacity to induce it can be mechanistically real while the resulting genetic lineage cannot persist. A three-snapshot ANTO event and SLK relative invasion sign do not solve absolute demographic viability.
2. **Phenotypic memory, selective opportunity, and growth are separately timed.** Under a shared environmental cycle, heritable priming can lag; a chance of generating a variant can peak late; the best deadline for accumulating descendants can peak earlier; and future population persistence can still be impossible.
3. **The same selection gradient does not fix the probability of rescue.** Need realized absolute births, deaths, competitive context and full-cycle fecundity, not only pollinator visitation, anther pollen release, pigment status, or rare-relative fitness.

## Biological evidence and falsifiable tests

- Original *Iochrominae* research reports several phylogenetic floral anthocyanin losses with maintained coding regions; see [Ho & Smith 2016](https://doi.org/10.1186/s12862-016-0675-3). ANTO's source-controlled 22-taxon colour/GenBank join currently has seven unpigmented taxa and 15 pigmented taxa, but **no joint functional induction or absolute reproductive outcome**. The 58 KT accessions are sequence records, not 58 independent evolutionary events; missing Dfr KT deposits or CDS annotations are not evidence of gene loss.
- The [2018 Iochrominae expression analysis](https://doi.org/10.1093/molbev/msy117) motivates tissue-specific expression control but does not provide a quantitative absolute demographic re-establishment probability.
- SLK's documented *Solanum*, *Chamaecrista* and *Bredia* pollen-mechanics contrasts motivate **multifunctional-anatomy** tests, not colour restoration. A pigmentation reversal per se does NOT identify a functional division-of-labor comparison.
- A definitive joint living-system study would measure a common plant genotype/lineage's hidden functional resource inventory, stage-matched induced trait and offspring heritability, the **same architectural candidate's** matched reproductive-function benefit and cost, its rare advantage **and absolute survival/fecundity** under multiple partner environments, and population history. Ancestral allele source and introgression must be resolved before labeling restoration intrinsic.

**Mandatory zero-result categories:** a transcript without verified function is not "competence"; an induced transient without heritable evidence is not "restoration"; a finite descendant increase with eventual extinction is not "rescue"; unavailable data are unmeasured, not zero.

**Repository boundaries:** This is a generic ANTO-theory extension linked to SLK's candidate-specific selection gates. All *Cirsium* empirical work remains exclusively at [aza3](https://github.com/zuizui0223/aza3). Do not replace SLK's Am Nat flagship empirical ceiling or promote it to an observed cross-clade evolutionary law.

## Reproduction and validation

- `src/anto_slk_absolute_demography_v3.py`
- `tests/test_slk_anto_absolute_demography_v3.py`: full-count independent branching oracle, threshold and input edge cases, long-run sign differences, conditional fitness mapping, and biological claim firewalls
- `scripts/run_slk_anto_absolute_demography_v3.py`: JSON receipt with 18-decimal-place output of 65-digit calculated numbers
- `.github/workflows/slk-anto-absolute-demography-v3.yml`: CI runs v2 and v3 tests and checks the receipt
- [SLK companion draft PR #119](https://github.com/zuizui0223/slk/pull/119) contains the biological interpretation, not a duplicated ANTO transition kernel.

No frequency-dependent competition, immigration, ecological partner evolution, or direct estimation of selection intensity is implemented here. A transition to such biology requires separate measurements and a fully specified population model, rather than multiplying this toy's q values into existing SLK fixation identities.
