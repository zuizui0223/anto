"""Conditional extension of ANTO v14: repeated expression attempts can be costly.

Resource scaffold C is an abstract common physiological/developmental
reproductive capacity. A/B are route-specific competences. Resource loss
is represented as irreversible for the *modeled opportunity episode*,
NOT as literal DNA deletion, germline mutation or a measured plant gene loss.

Raw attempts of route j occur at a_j. Each starts an independent abstract
Poisson branching descendant, ultimately nonextinct with q_j.
A failed attempt consumes:
  - its own route's capacity with probability route_loss_j;
  - common capacity C with probability common_loss_j;
  - neither with the remaining probability.
The two cost events are mutually exclusive; source event rates are thinned
by independent Bernoulli draws. The 'no consequence' failed attempts can
retry indefinitely while capacity remains. An A-success ends the "first
ultimately nonextinct route" analysis; B-ever controls ignore that censor.
All events occur in constant rates during one hypothetical opportunity period.

This is classical competing-risk/Poisson-thinning mathematics, not fitted
flower biology and not proof of ecological or genetic advantage of trait loss.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import exp, isfinite, log
from typing import Any

import numpy as np
from scipy.linalg import expm

from anto_afterlife_route_window_v12 import Rates, Loss, aged_inventory
from anto_afterlife_ecological_filter_v13 import poisson_lineage_nonextinction


def _unit_fraction(value:Any)->Fraction:
    if isinstance(value,bool):
        raise ValueError("Cost probability must lie in [0,1]")
    try:
        v=Fraction(str(value))
        if not 0 <= v <= 1 or not isfinite(float(v)):
            raise ValueError
    except (ValueError,TypeError,OverflowError,ZeroDivisionError) as e:
        raise ValueError("Cost probability must lie in [0,1]") from e
    return v


@dataclass(frozen=True)
class RetryCost:
    route_loss: Fraction
    common_loss: Fraction

    def __post_init__(self):
        for k in ("route_loss","common_loss"):
            object.__setattr__(self,k,_unit_fraction(getattr(self,k)))
        if self.route_loss+self.common_loss>1:
            raise ValueError("A failed attempt cannot destroy disjoint capacities with total probability >1")


@dataclass(frozen=True)
class CostlyRoutes:
    rates: Rates
    ancestral_failure: RetryCost
    alternative_failure: RetryCost

    def __post_init__(self):
        if not isinstance(self.rates,Rates) or not isinstance(self.ancestral_failure,RetryCost) or not isinstance(self.alternative_failure,RetryCost):
            raise ValueError("Expected v12 Rates and two declared RetryCost objects")


def _hazards(model:CostlyRoutes, qa:float, qb:float)->dict:
    if not isinstance(model,CostlyRoutes):
        raise ValueError("Invalid model instance")
    for q in (qa,qb):
        if not isfinite(float(q)) or q<0 or q>1:
            raise ValueError("Nonextinction probabilities must be [0,1]")
    r=model.rates
    a,b=float(r.ancestral),float(r.alternative)
    c,x,y=float(r.loss.common),float(r.loss.ancestral),float(r.loss.alternative)
    fa=a*(1-qa)
    fb=b*(1-qb)
    return {
        "success_A":a*qa,"success_B":b*qb,
        "cost_A_route":fa*float(model.ancestral_failure.route_loss),
        "cost_B_route":fb*float(model.alternative_failure.route_loss),
        "cost_A_common":fa*float(model.ancestral_failure.common_loss),
        "cost_B_common":fb*float(model.alternative_failure.common_loss),
        "baseline_common":c,"baseline_A":x,"baseline_B":y,
    }


def _div(a:float,b:float)->float:
    return a/b if b>0 else 0.


def exact_eventual(model:CostlyRoutes, qa:float, qb:float)->dict:
    """First-pass probabilities from all three transient inventory states."""
    h=_hazards(model,qa,qb)
    sA,sB=h["success_A"],h["success_B"]
    kA,kB=h["cost_A_route"],h["cost_B_route"]
    gA,gB=h["cost_A_common"],h["cost_B_common"]
    c,x,y=h["baseline_common"],h["baseline_A"],h["baseline_B"]
    from_A_only=_div(sA,sA+c+x+kA+gA)
    from_B_only=_div(sB,sB+c+y+kB+gB)
    total=sA+sB+c+x+y+kA+kB+gA+gB
    from_both_A=_div(sA+(y+kB)*from_A_only,total)
    from_both_B=_div(sB+(x+kA)*from_B_only,total)
    return {
        "from_A_only_A":from_A_only,
        "from_B_only_B":from_B_only,
        "from_both_A":from_both_A,
        "from_both_B":from_both_B,
        "hazards":h,
    }


def generator(model:CostlyRoutes,qa:float,qb:float)->np.ndarray:
    """Six-state CTMC in frozen v12 order [both,Aonly,Bonly,Asuccess,Bsuccess,lost]."""
    h=_hazards(model,qa,qb)
    sA,sB=h["success_A"],h["success_B"]
    kA,kB=h["cost_A_route"],h["cost_B_route"]
    gA,gB=h["cost_A_common"],h["cost_B_common"]
    c,x,y=h["baseline_common"],h["baseline_A"],h["baseline_B"]
    q=np.zeros((6,6),dtype=float)
    transitions=((0,1,y+kB),(0,2,x+kA),(0,3,sA),(0,4,sB),(0,5,c+gA+gB),
                 (1,3,sA),(1,5,c+x+kA+gA),
                 (2,4,sB),(2,5,c+y+kB+gB))
    for i,j,hazard in transitions:
        q[i,j]+=hazard
        q[i,i]-=hazard
    if not np.allclose(q.sum(axis=1),0,atol=1e-12):
        raise ArithmeticError("Transient generator does not preserve mass")
    return q


def after_silence(
    model:CostlyRoutes,age:float,
    mean_ancestral:float,mean_alternative:float,
    horizon:float|None=None
)->dict:
    """Probability of first ultimately surviving route after silent age t.

    Finite horizon bounds the date of the branch-initiating attempt, NOT
    the hypothetical infinite descendant time needed to adjudicate survival.
    """
    qa=poisson_lineage_nonextinction(mean_ancestral)
    qb=poisson_lineage_nonextinction(mean_alternative)
    p0=aged_inventory(model.rates.loss,age)
    if horizon is None:
        e=exact_eventual(model,qa,qb)
        fA=float(p0[0])*e["from_both_A"]+float(p0[1])*e["from_A_only_A"]
        fB=float(p0[0])*e["from_both_B"]+float(p0[2])*e["from_B_only_B"]
    else:
        if isinstance(horizon,bool) or not isfinite(float(horizon)) or float(horizon)<0:
            raise ValueError("Finite followup horizon must be finite nonnegative")
        p=p0 @ expm(generator(model,qa,qb)*float(horizon))
        fA,fB=float(p[3]),float(p[4])
    any_success=fA+fB
    return {
        "silent_age_model_units":float(age),
        "ancestral_first_nonextinct":fA,
        "alternative_first_nonextinct":fB,
        "any_nonextinct":any_success,
        "no_nonextinct_branch_by_endpoint":max(0.,1-any_success),
        "conditional_Poisson_nonextinction":{"A":qa,"B":qb},
        "resource_scope":"abstract shared reproductive/developmental capacity, not gene allele",
        "synthetic_only":True,
        "real_plant_recovery_demonstrated":False,
        "full_pathway_function_measured":False,
        "true_genetic_fixation_estimated":False,
        "claim_endpoint":"First ultimately nonextinct branch; not first observed pigment or actual plant establishment",
    }


def alternative_ever_without_success_preemption(
    model:CostlyRoutes,age:float,
    mean_ancestral:float,mean_alternative:float,
)->float:
    """Negative control: A may establish first but does not censor B.

    A successes leave resource inventory unchanged, A failed attempts still
    incur shared costs. This "coexistence" is an artificial comparison, not
    an agent-based ecological coexistence or true multiplex colonization.
    """
    qa=poisson_lineage_nonextinction(mean_ancestral)
    qb=poisson_lineage_nonextinction(mean_alternative)
    h=_hazards(model,qa,qb)
    sB=h["success_B"]
    kA,kB=h["cost_A_route"],h["cost_B_route"]
    gA,gB=h["cost_A_common"],h["cost_B_common"]
    c,x,y=h["baseline_common"],h["baseline_A"],h["baseline_B"]
    from_B_only=_div(sB,sB+c+y+kB+gB)
    from_both=_div(sB+(x+kA)*from_B_only,
                   sB+x+kA+y+kB+c+gA+gB)
    p=aged_inventory(model.rates.loss,age)
    return float(p[0])*from_both+float(p[2])*from_B_only


def subcritical_ancestral_peak(model:CostlyRoutes,mean_ancestral:float,
                               mean_alternative:float)->dict:
    """Sharp interior-max criterion in a restricted A-subcritical case.

    Assumptions: q_A=0, route-specific failure costs=0, B shared cost=0,
    no regeneration. Let k=c+y>0, x=baseline A loss>0,
    beta=b*q_B>0, h=a*g_A where g_A is A failed common cost.
      P(any eventually surviving route | t)
         = exp(-k*t) * u * [1 - h/(beta+x+k+h)*exp(-x*t)]
      u=beta/(beta+k).
    The curve peaks at positive t iff x*h > k*(beta+x+k).
    If g_A <= k*(beta+x+k)/(x*a), no positive peak exists.
    This is standard competing hazards with a harmful futile route, not
    a proven real biological 'loss benefit'.
    """
    qa=poisson_lineage_nonextinction(mean_ancestral)
    qb=poisson_lineage_nonextinction(mean_alternative)
    h=_hazards(model,qa,qb)
    if qa>0 or h["cost_A_route"] or h["cost_B_route"] or h["cost_B_common"]:
        raise ValueError("Peak formula requires subcritical A and no other failure costs")
    a=float(model.rates.ancestral)
    x=h["baseline_A"]
    k=h["baseline_common"]+h["baseline_B"]
    beta=h["success_B"]
    if a<=0 or x<=0 or k<=0 or beta<=0:
        raise ValueError("Peak formula requires positive opportunity, x, k, beta")
    d=h["cost_A_common"]
    threshold=k*(beta+x+k)/(x*a)
    result={"critical_common_cost_probability":threshold,
            "source_assumptions":"subcritical A; no A/B route-specific failure costs; no B shared failure cost",
            "equivalent_hazard_inequality":"x * a * gA > (c+y) * (beta+x+c+y)",
            "infinite_branch_retries_no_offspring_costs_measured":True}
    if x*d>k*(beta+x+k):
        t=log((k+x)*d/(k*(beta+x+k+d)))/x
        return {**result,"kind":"interior","peak_age":t}
    return {**result,"kind":"boundary_at_zero","peak_age":0.}


def cost_grid(rates:Rates,pre_loss:Loss)->list[dict]:
    """210 synthetic settings; no parameter prior or natural species frequency."""
    out=[]
    ages=(0.,1.,3.,5.,10.,20.)
    ras=(.95,1.01,1.05,1.2,3.)
    gammas=(0,.05,.1,.15,.2,.3,.6)
    for mean_a,gamma in product(ras,gammas):
        model=CostlyRoutes(rates,RetryCost(0,gamma),RetryCost(0,0))
        age_zero=after_silence(model,0,mean_a,1.5)["any_nonextinct"]
        for age in ages:
            result=after_silence(model,age,mean_a,1.5)
            value=result["any_nonextinct"]
            out.append({
                "mean_offspring_ancestral":mean_a,
                "mean_offspring_alternative":1.5,
                "failed_A_common_capacity_cost_probability":gamma,
                "silent_age":age,
                "any_nonextinct":value,
                "age_zero_same_model":age_zero,
                "cost_early_delay_raises_rescue":value>age_zero+1e-12,
                "cost_early_delay_lowers_rescue":value<age_zero-1e-12,
                "synthetic_only":True,
                "max_observed_flower_pigment_measurements":0,
            })
    return out
