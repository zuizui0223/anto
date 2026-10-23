"""Post-expression demographic *filter* on top of v12 first expression.

Models two INDEPENDENT components:
  (1) frozen v12 C,A,B first-recovery competing-risk probabilities;
  (2) a Poisson Galton-Watson descendant branching process following the
      first expression, under a STRICT one-shot route selection assumption.

Biological caveats:
- "nonextinction" = survival probability of the abstract unbounded, unregulated
  branch; NOT genetic fixation, ecological establishment, fitness measured in
  plants, or real seed production.
- A failed first-expressing lineage is not allowed another phenotype route
  under one-shot semantics. When no preemption is imposed, the claimed
  apparent B-benefit may vanish: this adverse control is indispensable.
- Real Cirsium maternal apical/lateral seed-route data are NOT estimates of
  ancestral-vs-alternative restored pigment reproduction parameters.
"""
from __future__ import annotations

from itertools import product
from math import expm1, isfinite, log1p
from typing import Any

from anto_afterlife_route_window_v12 import Rates, Loss, after_delay, alternative_ever_without_preemption


def _nonnegative_finite(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite nonnegative number")
    try:
        answer=float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite nonnegative number") from exc
    if not isfinite(answer) or answer<0:
        raise ValueError(f"{name} must be finite nonnegative number")
    return answer


def poisson_lineage_nonextinction(mean_offspring: float) -> float:
    """P(nonextinction) for a Poisson-offspring *independent* lineage.

    Solves q = 1 - exp(-R q) for positive root, using numerically stable
    bisection. For R<=1 extinction is certain, including the critical case.
    This classical branching-process formula is NOT a new theorem.
    """
    mu=_nonnegative_finite(mean_offspring,"mean_offspring")
    if mu<=1.:
        return 0.
    lo,hi=0.,1.
    for _ in range(160):
        mid=(lo+hi)/2
        f=-expm1(-mu*mid)-mid
        if f>0:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def one_attempt_ecological_filter(
    rates: Rates, pre_loss: Loss, age: float,
    mean_ancestral: float, mean_alternative: float,
    horizon: float | None = None,
) -> dict:
    """One v12 first recovery followed by a conditional unbounded branch.

    At most one successful *first phenotype choice* per founder: even if that
    expressing lineage later dies, the founder is not offered a new route.
    Each route's GW process starts with exactly one individual and has
    independent Poisson-mean offspring. No back-mutation, flow, density,
    environment fluctuations, migration or continued inter-route competition.
    """
    ra=_nonnegative_finite(mean_ancestral,"mean_ancestral")
    rb=_nonnegative_finite(mean_alternative,"mean_alternative")
    first=after_delay(rates,pre_loss,age,horizon)
    qa=poisson_lineage_nonextinction(ra)
    qb=poisson_lineage_nonextinction(rb)
    fa,fb=first["ancestral_first"], first["alternative_first"]
    a_live=fa*qa
    b_live=fb*qb
    partition={
        "ancestral_nonextinct":a_live,
        "alternative_nonextinct":b_live,
        "ancestral_first_then_extinct":fa*(1-qa),
        "alternative_first_then_extinct":fb*(1-qb),
        "no_first_expression":1-first["any_first"],
    }
    if abs(sum(partition.values())-1)>5e-12:
        raise ArithmeticError("Partition does not conserve starting founder mass")
    return {
        "silent_age_model_units":float(age),
        "first":first,
        "assumed_Poisson_mean_offspring":{"ancestral":ra,"alternative":rb},
        "conditional_lineage_nonextinction":{"ancestral":qa,"alternative":qb},
        "one_attempt_nonextinction":{"ancestral":a_live,"alternative":b_live,
                                     "any":a_live+b_live},
        "partition":partition,
        "abstract_branching_model_only":True,
        "true_genetic_fixation_or_real_plant_fitness_measured":False,
        "first_expression_preempts_other_route_after_failure":True,
    }


def fitness_route_crossover(rates:Rates, pre_loss:Loss,
                            mean_ancestral:float, mean_alternative:float,
                            bracket:tuple[float,float]) -> dict:
    """Locate a nonextinction-dominance crossover within a supplied interval.

    A sign-changing bracket is required. The result is a root of the
    *specified* one-shot model, not a biological evolutionary transition time
    or a theorem proving global uniqueness over all model parameters.
    """
    if not isinstance(bracket,(list,tuple)) or len(bracket)!=2:
        raise ValueError("bracket must have two ordered nonnegative ages")
    left=_nonnegative_finite(bracket[0],"left age")
    right=_nonnegative_finite(bracket[1],"right age")
    if right<=left:
        raise ValueError("bracket must have increasing ages")
    qa=poisson_lineage_nonextinction(mean_ancestral)
    qb=poisson_lineage_nonextinction(mean_alternative)
    def difference(t):
        p=after_delay(rates,pre_loss,t)
        return p["ancestral_first"]*qa - p["alternative_first"]*qb
    fl,fr=difference(left),difference(right)
    if fl==0:
        return {"kind":"endpoint","delay":left}
    if fr==0:
        return {"kind":"endpoint","delay":right}
    if fl*fr>0:
        return {"kind":"not_bracketed","delay":None,
                "endpoint_differences":[fl,fr]}
    for _ in range(110):
        middle=(left+right)/2
        fm=difference(middle)
        if fm==0:
            left=right=middle
            break
        if (fm>0)==(fl>0):
            left=middle
            fl=fm
        else:
            right=middle
    root=(left+right)/2
    return {"kind":"interior","delay":root,
            "balance_difference":difference(root),
            "interpretation":"equal first-expression × branch-nonextinction across the two routes"}


def minimum_alternative_offspring_mean_for_dominance(
    rates:Rates, pre_loss:Loss, age:float, mean_ancestral:float,
) -> dict:
    """Exact sharp *conditional* R_B threshold for competing first expression.

    Let p_A and p_B be v12 first-expression probabilities and q_A be the
    Poisson nonextinction chance at chosen R_A. Alternative persistence
    exceeds ancestral persistence if q_B > p_A*q_A/p_B.

    If target q_B in (0,1), invert q_B=1-exp(-R_B*q_B):
        R_B* = -log(1-q_B)/q_B.
    No finite mean B can overtake if target >=1. This is an elementary
    branching model inversion, NOT a new general eco-evolution theorem.
    """
    qa=poisson_lineage_nonextinction(mean_ancestral)
    first=after_delay(rates,pre_loss,age)
    pa,pb=first["ancestral_first"],first["alternative_first"]
    if pb<=0:
        return {"kind":"no_alternative_route","target_qB":None,
                "critical_mean_offspring":None}
    qtarget=pa*qa/pb
    if qtarget>=1:
        return {"kind":"not_attainable_for_finite_Poisson_R_B",
                "target_qB":qtarget,"critical_mean_offspring":None}
    if qtarget<=0:
        return {"kind":"any_supercritical_B_suffices",
                "target_qB":0.,"critical_mean_offspring":1.}
    critical=-log1p(-qtarget)/qtarget
    return {"kind":"finite","target_qB":qtarget,
            "critical_mean_offspring":critical,
            "strict_requirement":"R_B must be greater than this threshold for B to dominate"}


def nonpreemptive_alternative_probe(rates:Rates, pre_loss:Loss,
                                    age:float, mean_alternative:float,
                                    horizon:float|None=None)->dict:
    """Adverse control: B remains eligible after an earlier A event.

    In this special construction the resource and expression hazard for B are
    unaffected by prior A, and the branch attached to B is independent.
    Consequently B's event probability is its noncompeting 'ever' probability,
    monotonically declining with shared/B-resource age.
    """
    qb=poisson_lineage_nonextinction(mean_alternative)
    pb=alternative_ever_without_preemption(rates,pre_loss,age,horizon)
    return {
        "alternative_expression_ever":pb,
        "conditional_nonextinction":qb,
        "alternative_nonextinction":pb*qb,
        "adverse_control":"prior A expression does not censor B; no fitness competition modeled",
        "true_regain_and_fitness_in_plants_unmeasured":True,
    }


def sensitivity_grid(rates:Rates, pre_loss:Loss)->list[dict]:
    """540 deterministic toy combinations, not prevalence or fitted plant rates."""
    ages=(0.,2.,5.,7.,10.,20.)
    ra_values=(.8,.95,1.,1.05,1.1,1.5,2.,3.,5.)
    rb_values=(.8,.95,1.,1.05,1.1,1.2,1.5,2.,3.,5.)
    grid=[]
    for age,ra,rb in product(ages,ra_values,rb_values):
        value=one_attempt_ecological_filter(rates,pre_loss,age,ra,rb)
        fa,fb=value["first"]["ancestral_first"],value["first"]["alternative_first"]
        pa,pb=value["one_attempt_nonextinction"]["ancestral"],value["one_attempt_nonextinction"]["alternative"]
        grid.append({
            "silent_age_model_units":age,
            "Poisson_R_ancestral":ra,"Poisson_R_alternative":rb,
            "first_ancestral":fa,"first_alternative":fb,
            "nonextinction_ancestral":pa,"nonextinction_alternative":pb,
            "nonextinction_any":pa+pb,
            "first_alternative_more_but_ancestral_nonextinct_more":fb>fa and pa>pb,
            "first_alternative_more_and_alternative_nonextinct_more":fb>fa and pb>pa,
            "alternative_first_positive_but_alt_branching_zero":fb>0 and poisson_lineage_nonextinction(rb)==0,
            "model_only":True,
        })
    return grid
