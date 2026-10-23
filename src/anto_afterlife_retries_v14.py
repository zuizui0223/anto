"""Conditional retry-allowed recovery after phenotypic silence.

An exact finite-state *retrospective* model for the FIRST offspring branch
that will ultimately not become extinct after a phenotype re-expression.
All founder lineages begin phenotypically silent with common resource C
and route-specific resources A and B. Before opportunity returns, each
resource is irreversibly lost at a source-declared exponential rate.

Unlike v13's ONE-SHOT filter, each route's raw expression events may
recur while its functional resources remain. Each independent expression
attempt seeds an independent Poisson Galton-Watson branch with eventual
nonextinction probability q_j. Failed branches have NO cost, interaction,
effect on resources, or occupation of opportunity in this toy model.

By standard Poisson thinning, successful branch-founding attempts have
hazards a*q_A, b*q_B. The original v12 competitive resource-loss model
can therefore be reused EXACTLY without a new arbitrary simulation. The
winner is "first retrospectively successful branch", not the first visible
pigment event, fixation, fitness, plant seed success or a prospective
time-known establishment decision.

DO NOT fit actual Cirsium data: field maternal seed heads have neither
ancestral/alternative regained pigmentation nor heritable adult descendants.
"""
from __future__ import annotations

from itertools import product
from math import isfinite
from typing import Any

from anto_afterlife_route_window_v12 import (
    Rates, Loss, after_delay, aged_inventory, alternative_ever_without_preemption,
)
from anto_afterlife_ecological_filter_v13 import (
    poisson_lineage_nonextinction, one_attempt_ecological_filter,
)

def _valid_mean(mean: Any) -> float:
    if isinstance(mean,bool):
        raise ValueError("Mean offspring must be nonnegative finite")
    try:
        value=float(mean)
    except (TypeError, ValueError, OverflowError) as e:
        raise ValueError("Mean offspring must be nonnegative finite") from e
    if value < 0 or not isfinite(value):
        raise ValueError("Mean offspring must be nonnegative finite")
    return value

def effective_surviving_attempt_rates(
    raw_expression_rates: Rates, mean_ancestral: float, mean_alternative: float
) -> Rates:
    """Thinned clocks for events spawning ultimately surviving branches."""
    if not isinstance(raw_expression_rates,Rates):
        raise ValueError("Expected a valid v12 Rates instance")
    q_a=poisson_lineage_nonextinction(_valid_mean(mean_ancestral))
    q_b=poisson_lineage_nonextinction(_valid_mean(mean_alternative))
    # A successful branch is selected retrospectively from repeated
    # expression opportunities; all failed attempts leave toolkit unchanged.
    return Rates(float(raw_expression_rates.ancestral)*q_a,
                 float(raw_expression_rates.alternative)*q_b,
                 raw_expression_rates.loss)

def first_nonextinct_route_with_retries(
    raw_expression_rates: Rates, pre_opportunity_loss: Loss,
    silent_age: float, mean_ancestral: float, mean_alternative: float,
    horizon: float | None = None,
) -> dict:
    """P(first branch that is eventually nonextinct originates via A or B).

    This is *not* the probability of first visible expression. If horizon
    is finite, it limits the TIME of the ultimately surviving branch's
    initiating expression attempt, not the lifetime of its descendants.
    """
    effective=effective_surviving_attempt_rates(
        raw_expression_rates,mean_ancestral,mean_alternative)
    first=after_delay(effective,pre_opportunity_loss,silent_age,horizon)
    return {
        "silent_age_model_units":float(silent_age),
        "ancestral_first_nonextinct":first["ancestral_first"],
        "alternative_first_nonextinct":first["alternative_first"],
        "any_nonextinct":first["any_first"],
        "none_by_endpoint":first["no_recovery_by_endpoint"],
        "alternative_fraction_among_nonextinct_first":first["alternative_fraction_among_first_recoveries"],
        "effective_surviving_branch_hazards":{
            "ancestral":float(effective.ancestral),
            "alternative":float(effective.alternative)
        },
        "conditional_branch_nonextinction":{
            "ancestral":poisson_lineage_nonextinction(mean_ancestral),
            "alternative":poisson_lineage_nonextinction(mean_alternative)
        },
        "biological_model_semantics":"Retrospective first eventually surviving offspring branch, failed expression may recur without cost",
        "synthetic_only":True,
        "no_observed_plant_fitness_or_heritable_reexpression":True,
    }

def eventual_alternative_nonextinct_without_preemption(
    rates: Rates, pre_loss: Loss, silent_age: float,
    mean_alternative: float, horizon: float | None = None
) -> float:
    """P(any B-expressed branch eventually persists; A cannot censor B).

    Even an earlier successful A does not change the subsequent B hazard,
    C/B resource lifetime or the distribution of B opportunities.
    As age increases, this quantity NEVER increases in the declared model.
    """
    effective=effective_surviving_attempt_rates(
        rates,mean_ancestral=0,mean_alternative=mean_alternative)
    return alternative_ever_without_preemption(
        effective,pre_loss,silent_age,horizon)

def crossing_age_for_persistent_route(
    rates: Rates, pre_loss: Loss, mean_ancestral: float,
    mean_alternative: float, bracket: tuple[float,float]
) -> dict:
    """Find a sign-changing comparison of A-first vs B-first successful branch.

    Does NOT assert global uniqueness over arbitrary rate regimes or
    population growth histories. Returns not_bracketed when unavailable.
    """
    if not isinstance(bracket,(list,tuple)) or len(bracket)!=2:
        raise ValueError("Provide two model ages in a bracket")
    a,b=map(float,bracket)
    if (not isfinite(a) or not isfinite(b) or a<0 or b<=a):
        raise ValueError("Expected finite ordered nonnegative ages")
    eff=effective_surviving_attempt_rates(rates,mean_ancestral,mean_alternative)
    def diff(t):
        r=after_delay(eff,pre_loss,t)
        return r["ancestral_first"]-r["alternative_first"]
    left,right=diff(a),diff(b)
    if left==0:
        return {"kind":"endpoint","delay":a}
    if right==0:
        return {"kind":"endpoint","delay":b}
    if left*right>0:
        return {"kind":"not_bracketed","delay":None,
                "endpoint_differences":[left,right]}
    for _ in range(100):
        mid=(a+b)/2
        dm=diff(mid)
        if dm==0:
            a=b=mid
            break
        if (left>0)==(dm>0):
            a=mid
            left=dm
        else:
            b=mid
    root=(a+b)/2
    return {"kind":"interior","delay":root,
            "difference_at_root":diff(root),
            "claim_ceiling":"Crossover only for first eventually surviving lineage in the specified model"}

def coupling_grid(rates: Rates, pre_loss: Loss) -> list[dict]:
    """Deterministic 540-setting check, not plant prevalence or a prior."""
    ages=(0.,2.,5.,7.,10.,20.)
    ra_grid=(.8,.95,1.,1.05,1.1,1.5,2.,3.,5.)
    rb_grid=(.8,.95,1.,1.05,1.1,1.2,1.5,2.,3.,5.)
    rows=[]
    for ra,rb in product(ra_grid,rb_grid):
        shot0=one_attempt_ecological_filter(rates,pre_loss,0,ra,rb)["one_attempt_nonextinction"]["any"]
        retry0=first_nonextinct_route_with_retries(rates,pre_loss,0,ra,rb)["any_nonextinct"]
        for age in ages:
            newer=first_nonextinct_route_with_retries(rates,pre_loss,age,ra,rb)
            old=one_attempt_ecological_filter(rates,pre_loss,age,ra,rb)
            shot=old["one_attempt_nonextinction"]["any"]
            pre=aged_inventory(pre_loss,age)
            rows.append({
                "silent_age":age,"mean_offspring_ancestral":ra,
                "mean_offspring_alternative":rb,
                "retry_ancestral":newer["ancestral_first_nonextinct"],
                "retry_alternative":newer["alternative_first_nonextinct"],
                "retry_any":newer["any_nonextinct"],
                "one_shot_any":shot,
                "probability_any_initial_capacity":float(pre[0]+pre[1]+pre[2]),
                "retry_any_at_age0_for_same_R":retry0,
                "one_shot_any_at_age0_for_same_R":shot0,
                "one_shot_total_increase_but_retry_total_decrease":(
                    shot > shot0+1e-12 and newer["any_nonextinct"] < retry0-1e-12),
                "alternative_first_success_increased_over_age0":(
                    newer["alternative_first_nonextinct"] >
                    first_nonextinct_route_with_retries(
                        rates,pre_loss,0,ra,rb)["alternative_first_nonextinct"]+1e-12),
                "synthetic_only":True,
            })
    return rows
