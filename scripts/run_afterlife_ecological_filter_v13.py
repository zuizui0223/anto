"""Source-pinned reproduction of ANTO v13 synthetic demographic filtering.

The measured endpoint is NOT plant fitness/establishment: one-off first
pigment-expression route × nonextinction of a hypothetical Poisson offspring
Galton–Watson descendant branch. Does not sample real plant data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_afterlife_ecological_filter_v13 import (
    one_attempt_ecological_filter,
    poisson_lineage_nonextinction,
    fitness_route_crossover,
    nonpreemptive_alternative_probe,
    minimum_alternative_offspring_mean_for_dominance,
    sensitivity_grid,
)
from anto_afterlife_route_window_v12 import Rates,Loss


def report():
    loss=Loss(".02",".30",".01")
    rates=Rates(".40",".20",loss)
    mean_a=3.0
    mean_b=1.05
    ages=(0.,5.,7.,20.)
    data=[]
    for age in ages:
        d=one_attempt_ecological_filter(rates,loss,age,mean_a,mean_b)
        d["critical_R_B_to_outcompete_ancestral"]=minimum_alternative_offspring_mean_for_dominance(rates,loss,age,mean_a)
        data.append(d)

    crossing=fitness_route_crossover(rates,loss,mean_a,mean_b,(0,20))
    if crossing["kind"]!="interior":
        raise AssertionError("Original demonstration is not bracketed")
    primary=one_attempt_ecological_filter(rates,loss,5,mean_a,mean_b)
    finite=one_attempt_ecological_filter(rates,loss,5,mean_a,mean_b,horizon=10)

    # Strongly deliberately adverse case: the A phenotype appears first,
    # but its descendant mean is subcritical. Treating A as terminal preempts
    # the B lineage. B nevertheless remains biologically available in
    # a counterfactual with no such preemption.
    temp0=one_attempt_ecological_filter(rates,loss,0,.95,1.5)
    temp5=one_attempt_ecological_filter(rates,loss,5,.95,1.5)
    free0=nonpreemptive_alternative_probe(rates,loss,0,1.5)
    free5=nonpreemptive_alternative_probe(rates,loss,5,1.5)
    spurious_gain=(
        temp5["one_attempt_nonextinction"]["any"]>temp0["one_attempt_nonextinction"]["any"] and
        free5["alternative_nonextinction"]<free0["alternative_nonextinction"]
    )
    if not spurious_gain:
        raise AssertionError("Essential adverse control failed")

    grid=sensitivity_grid(rates,loss)
    n=Counter()
    for g in grid:
        for label in ("first_alternative_more_but_ancestral_nonextinct_more",
                      "first_alternative_more_and_alternative_nonextinct_more",
                      "alternative_first_positive_but_alt_branching_zero"):
            n[label]+=int(g[label])

    source=ROOT/"src/anto_afterlife_ecological_filter_v13.py"
    return {
        "schema":"anto_afterlife_ecological_filter_v13",
        "analysis_date":"2026-10-23",
        "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
        "v12_source_model":"src/anto_afterlife_route_window_v12.py",
        "synthetic_only":True,
        "real_plant_mean_offspring_estimated":False,
        "genetic_fixation_or_true_ecological_establishment_demonstrated":False,
        "first_phenotype_event_preempts_subsequent_routes_even_after_branch_extinction":True,
        "biological_caveat":"Poisson Galton-Watson branch without density, migration, ecological seasonality, mate limitation, genetic composition or re-expression retries",
        "rates_v12":{"activation_ancestral":0.4,"activation_alternative":0.2,
                     "common_loss":0.02,"ancestral_loss":0.30,"alternative_loss":0.01},
        "conditional_mean_offspring":{"ancestral":mean_a,"alternative":mean_b},
        "conditional_branching_nonextinction":{"ancestral":poisson_lineage_nonextinction(mean_a),
                                              "alternative":poisson_lineage_nonextinction(mean_b)},
        "ancestral_vs_alternative_crossover_delay":crossing["delay"],
        "crossover_root":crossing,
        "primary_at_age_5":primary,
        "finite_first_event_horizon_10_at_age_5":finite,
        "reference_cases":data,
        "sensitivity_grid":{"settings":len(grid),"counts":dict(n),
                            "grid_has_no_natural_frequency_interpretation":True},
        "no_preemption_adverse_control":{
            "spurious_gain_removed":spurious_gain,
            "one_shot_with_subcritical_ancestral":{
                "mean_ancestral":.95,"mean_alternative":1.5,
                "time0_nonextinction_any":temp0["one_attempt_nonextinction"]["any"],
                "time5_nonextinction_any":temp5["one_attempt_nonextinction"]["any"]},
            "nonpreemptive_B_ever_branch_survival":{
                "time0":free0["alternative_nonextinction"],
                "time5":free5["alternative_nonextinction"]},
            "interpretation":"One-shot selection creates a B-filter benefit absent if A-first does not preclude later B expression"
        },
        "original_study_outcome_control":{
            "aza3_C_undulatum_105_maternal_apical_lateral_routes_are_not_fitted_as_R_A_R_B":True,
            "same_mother_head_seed_yield_is_not_an_ancient_flower_pigment_regain_event":True,
            "egwe_simulation_is_not_independent_empirical_confirmation":True,
        },
        "software_environment":{"python":platform.python_version()},
        "next_real_plant_gate":"Independently measured stage-specific flower pigment, actual functional coding and regulation, descendant survival/reproductive maturity, ancestry, and explicit re-expression retries"
    },grid


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--grid-csv",type=Path)
    args=parser.parse_args()
    data,grid=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    if args.grid_csv:
        args.grid_csv.parent.mkdir(parents=True,exist_ok=True)
        with args.grid_csv.open("w",newline="",encoding="utf-8") as f:
            fields=list(grid[0].keys())
            writer=csv.DictWriter(f,fieldnames=fields)
            writer.writeheader()
            writer.writerows(grid)
    print("ANTO_V13_CONDITIONAL_ECOLOGICAL_FILTER",
          json.dumps({"crossover":data["ancestral_vs_alternative_crossover_delay"],
                      "grid_settings":data["sensitivity_grid"]["settings"],
                      "grid_conditions":data["sensitivity_grid"]["counts"],
                      "no_preemption_adverse_control":data["no_preemption_adverse_control"]["spurious_gain_removed"],
                      "natural_plants_calibrated":False},sort_keys=True))


if __name__=="__main__":
    main()
