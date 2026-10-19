"""Exact, genus-agnostic history/hidden-state sufficiency experiment.

Question: does dispersal history *directly* determine later phenotype
appearance, or does its predictive information reside in hidden resource
states that plant propagules carried and retained?

Construct histories A (even-parity module inventory) and B (odd-parity).
Each has IDENTICAL one-module frequencies and pairwise co-retention.
All current visible flowers are stipulated white. Yet the probability of
joint availability of 3 necessary resources differs (0 vs 1/4).

A future phenotypic event occurs with probability alpha iff all three
resources co-occur in a propagule. Thus its event probability by history
is 0 vs alpha/4. Once the FULL retained module inventory is measured,
history adds no further prediction under the stipulated Markov/sufficiency
assumptions. This is a mathematical mechanism-identifiability counterexample,
NOT a measured plant rate, direct trait reversal or novel Markov theorem.

We calculate exact Bayes-optimal expected Brier errors, not simulated fits.
Independent binary marker assays have explicit sensitivity/specificity;
perfect measurement is NOT assumed for real floral genomic functions.
"""
from fractions import Fraction
from itertools import product
from typing import Mapping

from anto_joint_accessibility_v7 import MASKS,pairwise_constrained_family

HISTORIES=("history_even","history_odd")
OBS=tuple(MASKS)
HALF=Fraction(1,2)
ZERO=Fraction(0)
ONE=Fraction(1)
DEFAULT_ALPHA=Fraction(2,5)


def _prob(x,name):
    try:q=Fraction(x)
    except (TypeError,ValueError,ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a finite rational probability") from exc
    if q<0 or q>1:raise ValueError(f"{name} must be in [0,1]")
    return q


def _history_inventories():
    return {
        "history_even":pairwise_constrained_family(ZERO),
        "history_odd":pairwise_constrained_family(Fraction(1,4)),
    }


def _likelihood(observed,true,sensitivity,specificity):
    p=ONE
    for measured,actual in zip(observed,true):
        if actual=="1":
            p*=sensitivity if measured=="1" else ONE-sensitivity
        else:
            p*=ONE-specificity if measured=="1" else specificity
    return p


def _outcomes(success_if_complete=DEFAULT_ALPHA,
              sensitivity=Fraction(9,10),
              specificity=Fraction(19,20),
              history_inventories=None):
    success_if_complete=_prob(success_if_complete,"success_if_complete")
    sensitivity=_prob(sensitivity,"sensitivity")
    specificity=_prob(specificity,"specificity")
    d=_history_inventories() if history_inventories is None else history_inventories
    if set(d)!=set(HISTORIES):
        raise ValueError("Must supply both named independent history scenarios")

    rows=[]
    for history in HISTORIES:
        p=d[history]
        if set(p)!=set(MASKS) or sum(p.values())!=1:
            raise ValueError("Each history requires an 8-state normalized inventory distribution")
        for inventory,prob in p.items():
            if prob<0:raise ValueError("Inventory probability cannot be negative")
            mass=HALF*Fraction(prob)
            if not mass:continue
            event=success_if_complete if inventory=="111" else ZERO
            for z in OBS:
                assay=_likelihood(z,inventory,sensitivity,specificity)
                if assay:
                    rows.append({
                        "history":history,"inventory":inventory,
                        "measured_bits":z,"mass":mass*assay,
                        "event_probability":event,
                    })
    if sum(r["mass"] for r in rows)!=ONE:
        raise ArithmeticError("Generative probabilities do not sum to one")
    return rows


def _optimal_brier(rows,features):
    """Population risk of Bayes forecast E[Y | selected measured features]."""
    groups={}
    for r in rows:
        key=tuple(r[field] for field in features)
        mass,events=groups.get(key,(ZERO,ZERO))
        groups[key]=(mass+r["mass"],events+r["mass"]*r["event_probability"])
    out=ZERO
    for mass,events in groups.values():
        if mass:
            q=events/mass
            out+=mass*q*(ONE-q)
    return out


def compare_models(success_if_complete=DEFAULT_ALPHA,
                   sensitivity=Fraction(9,10),
                   specificity=Fraction(19,20),
                   histories_differ=True):
    if histories_differ:
        d=_history_inventories()
    else:
        common=pairwise_constrained_family(Fraction(1,8))
        d={h:common for h in HISTORIES}
    rows=_outcomes(success_if_complete,sensitivity,specificity,
                   history_inventories=d)
    models={
        "visible_flower_only":(),
        "visible_plus_history":("history",),
        "visible_plus_noisy_module_assays":("measured_bits",),
        "visible_plus_history_plus_noisy_assays":("history","measured_bits"),
        "visible_plus_complete_functional_inventory":("inventory",),
        "visible_plus_complete_inventory_plus_history":("inventory","history")
    }
    scores={name:_optimal_brier(rows,features) for name,features in models.items()}
    if scores["visible_plus_complete_functional_inventory"] != scores["visible_plus_complete_inventory_plus_history"]:
        raise ArithmeticError("History unexpectedly predicts after sufficient complete state")
    # A more informative variable cannot increase Bayes-optimal Brier error.
    if not (scores["visible_flower_only"]>=scores["visible_plus_history"]
            >=scores["visible_plus_history_plus_noisy_assays"]
            >=scores["visible_plus_complete_functional_inventory"]):
        raise ArithmeticError("Information monotonicity failed")
    return {
        "definition":"Exact Bayes-optimal forecast; all currently visible flowers white",
        "synthetic":True,
        "no_plant_data_fitted":True,
        "event_requires_all_three_modules":True,
        "history_classes_same_single_and_pairwise_module_moments":True,
        "success_if_all_resources":str(_prob(success_if_complete,"success_if_complete")),
        "assay_sensitivity":str(_prob(sensitivity,"sensitivity")),
        "assay_specificity":str(_prob(specificity,"specificity")),
        "histories_differ_in_third_order_joint_inventory":bool(histories_differ),
        "history_specific_future_event_probability":{
            h:str(sum(r["mass"]*r["event_probability"] for r in rows if r["history"]==h)/HALF)
            for h in HISTORIES
        },
        "bayes_optimal_brier_exact":{k:str(v) for k,v in scores.items()},
        "bayes_optimal_brier_decimal":{k:float(v) for k,v in scores.items()},
        "history_predictive_gain_before_complete_state":
            float(scores["visible_flower_only"]-scores["visible_plus_history"]),
        "state_measurement_predictive_gain_after_history":
            float(scores["visible_plus_history"]-scores["visible_plus_complete_functional_inventory"]),
        "residual_history_gain_after_complete_state":
            float(scores["visible_plus_complete_functional_inventory"]-
                  scores["visible_plus_complete_inventory_plus_history"]),
        "qualified_interpretation":"History predicts only through stipulated hidden resource inventory. If history still predicts after a putatively complete present state is measured, test omitted regulatory, epigenetic, age, gene-flow or environmental factors; no literal evidence of causal history beyond an adequate Markov state."
    }


def comparative_sensitivity_grid():
    cases=[
        ("uninformative_assay",Fraction(1,2),Fraction(1,2)),
        ("weak_assay",Fraction(13,20),Fraction(13,20)),
        ("moderate_assay",Fraction(4,5),Fraction(4,5)),
        ("strong_assay",Fraction(9,10),Fraction(19,20)),
        ("perfect_assay",ONE,ONE),
    ]
    out={}
    for name,s,sp in cases:
        out[name]=compare_models(sensitivity=s,specificity=sp)
    out["matched_history_null"]=compare_models(histories_differ=False)
    out["no_success_null"]=compare_models(success_if_complete=ZERO)
    return {
        "version":"anto_history_state_sufficiency_v9",
        "central_claim":"Historical prediction can be mediated entirely by present hidden resources; a full sufficient state removes residual past dependence by construction.",
        "empirical_generalization_demonstrated":False,
        "models":out
    }


if __name__=="__main__":
    import json
    print(json.dumps(comparative_sensitivity_grid(),indent=2,ensure_ascii=False))
