"""Reproduce retry-allowed conditional ANTO v14 first nonextinct route analysis."""
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
from anto_afterlife_route_window_v12 import Rates,Loss
from anto_afterlife_ecological_filter_v13 import one_attempt_ecological_filter
from anto_afterlife_retries_v14 import (
    first_nonextinct_route_with_retries, crossing_age_for_persistent_route,
    eventual_alternative_nonextinct_without_preemption, coupling_grid,
)

def report():
    rates=Rates(".4",".2",Loss(".02",".30",".01"))
    loss=rates.loss
    ages=(0,5,7,20)
    cases=[]
    for age in ages:
        new=first_nonextinct_route_with_retries(rates,loss,age,3,1.05)
        one=one_attempt_ecological_filter(rates,loss,age,3,1.05)
        new["v13_one_shot_surviving_lineages"]=one["one_attempt_nonextinction"]
        new["alternative_ever_surviving_if_A_cannot_censor_B"]=(
            eventual_alternative_nonextinct_without_preemption(
                rates,loss,age,1.05))
        cases.append(new)
    crossing=crossing_age_for_persistent_route(rates,loss,3,1.05,(0,20))
    adverse={
        "one_shot_at_age0":one_attempt_ecological_filter(
            rates,loss,0,.95,1.5)["one_attempt_nonextinction"]["any"],
        "one_shot_at_age5":one_attempt_ecological_filter(
            rates,loss,5,.95,1.5)["one_attempt_nonextinction"]["any"],
        "retry_at_age0":first_nonextinct_route_with_retries(
            rates,loss,0,.95,1.5)["any_nonextinct"],
        "retry_at_age5":first_nonextinct_route_with_retries(
            rates,loss,5,.95,1.5)["any_nonextinct"],
        "alternative_without_A_preemption_age0":
            eventual_alternative_nonextinct_without_preemption(
                rates,loss,0,1.5),
        "alternative_without_A_preemption_age5":
            eventual_alternative_nonextinct_without_preemption(
                rates,loss,5,1.5),
    }
    adverse["one_shot_spurious_gain_disappears_when_failed_expression_does_not_preempt"]=(
        adverse["one_shot_at_age5"]>adverse["one_shot_at_age0"] and
        adverse["retry_at_age5"]<adverse["retry_at_age0"] and
        adverse["alternative_without_A_preemption_age5"]<
        adverse["alternative_without_A_preemption_age0"]
    )
    if not adverse["one_shot_spurious_gain_disappears_when_failed_expression_does_not_preempt"]:
        raise AssertionError("Important adverse control not recovered")
    grid=coupling_grid(rates,loss)
    counts=Counter()
    for r in grid:
        for field in ("one_shot_total_increase_but_retry_total_decrease",
                      "alternative_first_success_increased_over_age0"):
            counts[field]+=int(r[field])
        if not (r["retry_any"]>=r["one_shot_any"]-1e-12 and
                r["retry_any"]<=r["probability_any_initial_capacity"]+1e-12 and
                r["retry_any"]<=r["retry_any_at_age0_for_same_R"]+1e-12):
            raise AssertionError("Mathematically impossible retry bound in grid")
    p={
        "schema":"anto_afterlife_retries_v14",
        "analysis_date":"2026-10-23",
        "synthetic_only":True,
        "natural_phylogeny_or_reversion_rate_estimated":False,
        "observed_plant_reexpression_events":0,
        "population_genetics_and_real_fitness_measured":False,
        "model_endpoint":"Retrospective first ever-offspring branch that would remain nonextinct in an unlimited Poisson Galton-Watson process",
        "attempt_semantics":"Expression may retry after failed independent branch; failure does not deplete genetic resources or exclude alternative route",
        "limit":"No ecological feedback from temporarily extant failed branches, no genetic switching, no density dependence or environmental stochasticity",
        "source_code_sha256":hashlib.sha256(
            (ROOT/"src/anto_afterlife_retries_v14.py").read_bytes()).hexdigest(),
        "earlier_commit_reference":"a1f804dd1bc989b735fdeb790e421fad6c8d6232",
        "original_v12_rates":{"A":.4,"B":.2,"loss_C":.02,"loss_A":.30,"loss_B":.01},
        "hypothetical_mean_offspring":{"A":3.,"B":1.05},
        "comparison_to_v13_crossover_delay":5.8896093398868,
        "crossover_delay":crossing["delay"],
        "crossover":crossing,
        "reference_cases":cases,
        "adverse_control":adverse,
        "sensitivity_grid":{
            "settings":len(grid),
            "counts":dict(counts),
            "grid_is_not_a_natural_species_frequency_or_parameter_prior":True},
        "existing_biological_records_are_not_fitted_branching_inputs":{
            "Cirsium_undulatum_maternal_apical_lateral_seed_route_105_plants":True,
            "Ryukyu_Cirsium_young_leaf_RNA_reads_not_mature_plant_descendant_lineages":True,
        },
        "software":{"python":platform.python_version()},
        "next_admissible_empirical_test":
            "Repeated developmental expression and true descendant survival, same genetic lineages, adjudicated pathway identity and source-confounded conditions",
    }
    if abs(p["crossover_delay"]-2.17240158406441)>1e-7:
        raise AssertionError("Unexpected conditional model crossover")
    return p,grid

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--grid",type=Path)
    args=parser.parse_args()
    result,grid=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    if args.grid:
        args.grid.parent.mkdir(parents=True,exist_ok=True)
        with args.grid.open("w",encoding="utf-8",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=list(grid[0]))
            wr.writeheader()
            wr.writerows(grid)
    print("ANTO_V14_RETRIES_VERIFIED",json.dumps({
        "crossover":result["crossover_delay"],
        "grid":result["sensitivity_grid"],
        "adverse_control":result["adverse_control"],
        "no_natural_plant_rate_fitted":True
    },sort_keys=True))

if __name__=="__main__":
    main()
