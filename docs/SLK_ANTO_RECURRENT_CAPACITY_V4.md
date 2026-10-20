# SLK × ANTO v4 — recurrent capacity can persist without continuous expression

2026-10-20. **Conditional theoretical result; no new plant observations.**

## Biological conclusion

A lineage's capacity to recreate a phenotype is not the same as uninterrupted
inheritance of that phenotype. When nonexpressing but competent descendants can
reproduce, they can keep the capacity alive and repeatedly supply expressing
descendants. Conversely, directing more of those descendants into an expressing
state can improve the immediate visible return while making the entire
competent lineage subcritical. This reversal requires a demographic advantage
of the nonexpressing state in the specified environment; it is not universal.

This is a cleaner endpoint than the preceding v3 census: **no generation-14
cutoff, four-descendant target, 40-founder sample size, or silent-priming state
is needed.** The outcome is the infinite-horizon nonextinction probability of
all descendants retaining the necessary toolkit in a constant environment.
It is NOT genetic fixation, whole-organism population survival, or an observed
plant re-evolution event.

## What changed from v3

The prior [absolute-demography model](SLK_ANTO_ABSOLUTE_DEMOGRAPHY_V3.md), pinned
at `195f87be7b8599adff999542d89e2a5e014c3d25`, assumes that an activated founder's
phenotype is inherited thereafter. v4 permits both silencing and renewed
expression among descendants and lets the required toolkit be lost AFTER
re-expression. v1–v3 are retained unchanged as explicitly narrower models.

The two states are D, expressing the focal trait, and C, not expressing it but
retaining the complete functional resource. **C is reproductively active, not
a dormant seed.** Offspring losing the toolkit exit this two-type process;
they may remain viable organisms, so toolkit extinction is not plant extinction.
Loss is irreversible in the present model: immigration, introgression and new
construction could reopen routes but are not included.

## Defined reproduction and transmission

Each parent has two independent reproductive slots. A type-i slot produces a
viable offspring with probability q_i. Conditional on birth, that offspring
retains the necessary toolkit with probability 1-d. Conditional on retention,
a D parent's child remains D with probability 1-s and becomes C with probability
s; a C parent's child becomes D with probability a and stays C otherwise.
Switching is an assumed heritable or transgenerational state process, not an
experimentally established pigment-regulatory mechanism.

The mean competent-offspring matrix (rows parents D,C; columns offspring D,C) is

```
M = (1-d) * [[2*q_D*(1-s), 2*q_D*s],
             [2*q_C*a,     2*q_C*(1-a)]].
```

For the declared nonsingular irreducible domain, a positive probability of
competent-lineage survival exists exactly when the largest eigenvalue rho(M)
exceeds one. At or below one, extinction of the competent lineage is certain.
This is a classical multitype branching-process result, not a new theorem.
Because both switching directions are positive, surviving competent genealogies
produce expressing descendants at arbitrarily late generations almost surely;
this does not imply that expression is continuous or fixed throughout the tree.

With slot PGFs `F_i(z)=(1-sum(M_i)/2 + M_i*z/2)^2`, extinction is the smallest
fixed point of F. Two-sided iteration is evaluated with 50-digit Decimal
arithmetic. Returned numerical brackets are not outward-rounded interval
certificates. Criticality is classified with exact rational arithmetic, so
slow convergence near the boundary cannot turn a critical case into apparent
positive permanence.

## Result 1 — more immediate re-expression, less indefinite recoverability

Hold `q_D=0.65`, `q_C=0.90`, `d=0.30`, and `s=0.25` fixed. These illustrative
rates are not estimated in any species. Change only activation a:

| Conditional activation a | At least one D child next generation, from one competent C | rho(M) | Indefinite competent-lineage survival, from one competent C |
|---:|---:|---:|---:|
| 0.30 | 34.2279% | 1.092000 | 35.6158% |
| 0.50 | 53.0775% | 1.035742 | 16.2928% |
| 0.90 | 81.2511% | 0.983394 | **0%** |

The first column of outcomes is `1-[1-q_C*(1-d)*a]^2`; it is an immediate
expression event, not the conditional activation parameter itself. Long-run
survival is from the SAME single competent silent founder in all scenarios.
No arbitrary large-founder multiplier creates this difference.

Here C has higher viable offspring output than D. Increasing a sends a greater
share of C descendants into the less productive D state, progressively reducing
the route that maintains capacity. The important ecological hypothesis is not
that more expression is harmful everywhere, but that **current expression and
retention of a reproductive lineage capable of future expression can oppose
one another when the silent state has a demographic advantage**.

If D instead has higher competent-offspring output than C, the direction of the
spectral effect reverses. With equal output it disappears. The model does not
simulate evolutionary change in a and does not identify an evolutionarily
optimal switching rate.

## Result 2 — an exact boundary, not one chosen simulation

Write `v_D=2*q_D*(1-d)` and `v_C=2*q_C*(1-d)`. In the case `v_D<1<v_C`, an
interior crossing occurs at

```
a* = [1-v_D*(1-s)]*(v_C-1) / [v_C*(1-v_D)].
```

When `0<a*<1`, positive competent-lineage survival holds for `0<a<a*` and
fails for `a>=a*`. For the stated example, `a*=1651/2268 = 0.7279541446...`.
This percentage is a boundary of the declared synthetic slice, not a biological
constant. The exact identity follows from

```
det(I-M) = [1-v_D*(1-s)]*(1-v_C) + v_C*a*(1-v_D).
```

More generally the derivative of the largest eigenvalue with respect to a is
`-v_C*(rho-v_D)/(2*rho-trace(M))`. Row-sum bounds put rho strictly between
v_D and v_C when they differ and switching is irreducible. Thus higher
activation reduces spectral growth when v_C>v_D, increases it when v_D>v_C,
and leaves it unchanged when the two row sums coincide.

At a=0.30 the diagonal entries are `M_DD=0.6825` and `M_CC=0.882`, BOTH below
one, yet `rho(M)=1.092>1`. Each state considered only through its same-state
children is subcritical, but cross-state contributions sustain the joint
competent genealogy. Removing cross-type offspring is a bookkeeping control,
not automatically a feasible biological knockout preserving all other rates.

Common toolkit loss gives another exact boundary: if M0 is the zero-loss
matrix, `rho(M)=(1-d)*rho(M0)`, so survival requires `d<1-1/rho(M0)`.
At a=0.30 this critical loss is approximately 0.35897436. If rho(M0)<=1,
there is no nonnegative loss rate allowing survival.

## Result 3 — initial hidden resources and long-run dynamics have different roles

If only a fraction theta of founders retains the toolkit and all competent
founders initially are C, the overall nonextinction probability for one sampled
founder is `theta*survival_C`. For positive theta, initial resource abundance
changes that probability but not the spectral boundary. No finite number of
initial founders can make a genuinely subcritical closed process persist
forever. This is specific to the unbounded branching approximation with constant
rates, not a claim about density-dependent populations or recurrent gene flow.

The v3 perfect-inheritance baseline is recovered as a control: with equal
q_D=q_C=0.65 and zero toolkit loss, total competent offspring are Binomial(2,0.65)
regardless of expression switching, giving survival `120/169`. Adding loss d=0.25
makes the effective per-slot probability 0.4875 and survival zero. Loss of
recoverable capacity after expression therefore cannot be silently omitted.

## Actual connection to SLK, and what is deliberately not connected

SLK's [scaffold-preserving alternative](https://github.com/zuizui0223/slk/blob/3d6dea545b3389a73e1ceb369279a01af15a5d30/theory/FLORAL_SCAFFOLD_BYPASS_V1.md)
shows why beating one integrated reference is not equivalent to beating every
viable alternative. C is an additional viable, nonexpressing state here. The
example makes C demographically better than D, so it would be inconsistent to
also call D favored over C merely because a previous D-versus-S SLK comparison
was positive. **No SLK value/path/invasion mask is imported as if it calibrated
these offspring rates.** The whole candidate comparison must be remeasured.

ANTO supplies the retained-versus-lost functional toolkit question; SLK supplies
the requirement to compare credible alternative organizations and ecological
returns. This supplement identifies a candidate-set and lineage-persistence
hypothesis. It is not a calibrated joint SLK–ANTO population model or a fourth
SLK selection threshold. A colour return alone is not anatomical division of
labour. All Cirsium empirical ownership remains with `zuizui0223/aza3`.

## Prior art and the remaining biological test

Carja & Plotkin (2019), *Genetics*, DOI
[10.1534/genetics.118.301758](https://doi.org/10.1534/genetics.118.301758), already
study partly heritable phenotypic variability and nonmonotonic persistence.
Blath, Hermann & Slowik (2021), *Journal of Mathematical Biology*, DOI
[10.1007/s00285-021-01639-6](https://doi.org/10.1007/s00285-021-01639-6), use
active/dormant two-type branching processes and resource-aware comparisons.
Switching-assisted persistence, the spectral criterion and nonmonotonic
switching effects are therefore not claimed as discoveries of this project.
The difference between their dormancy state and the reproducing C state here
must not be obscured.

A biological test would have to distinguish complete pathway retention from
expression, measure inherited transitions and full-cycle reproduction of both
states under the same partner environment, and demonstrate the same trait's
functional comparison rather than assign q from photographs. The decisive
comparison is whether greater immediate re-expression coincides with lower
long-term persistence of recoverable capacity after accounting for the
nonexpressing state. No such plant result has been obtained in this update.

## Reproduction and numerical verification

- Model: `src/anto_slk_recurrent_capacity_v4.py`.
- Tests: `tests/test_slk_anto_recurrent_capacity_v4.py` (14 tests).
- Runner: `scripts/run_slk_anto_recurrent_capacity_v4.py`.
- Summary: `results/slk_anto_recurrent_capacity_v4_summary.json`.
- Workflow: `.github/workflows/slk-anto-recurrent-capacity-v4.yml`.

The short-generation oracle enumerates the joint distribution of D and C
counts independently of the PGF recursion. A separate local NumPy eigensolver
checked 324 parameter combinations; all classifications agreed and all
50-digit survival brackets converged to the requested tolerance. A near-critical
regression reproduces a binary64 plateau and checks against an independently
computed 70-digit root. Full repository CI is reported separately, not inferred
from this local focused suite. Existing v1–v3 source is unchanged.
