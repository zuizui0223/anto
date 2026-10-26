"""Deterministically reproduce ANTO v15 costly retry model-only results.

This emits a machine-verifiable *synthetic* comparison to the zero-cost
v14 boundary. Original Cirsium seed data and Petunia perturbation papers
are cited as context, NOT imported as trait-reversion/failure-cost rates.
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
from anto_afterlife_route_window_v12 import Loss,Rates
from anto_afterlife_costly_retries_v15 import (
    RetryCost,CostlyRoutes,after_silence,
    alternative_ever_without_success_preemption,subcritical_ancestral_peak,
    cost_grid,
)

def report():
    rates=Rates(".4",".2",Loss(".02",".3",".01"))
    no=CostlyRoutes(rates,RetryCost(0,0),RetryCost(0,0))
    costly=CostlyRoutes(rates,RetryCost(0,".3"),RetryCost(0,0))
    below=CostlyRoutes(rates,RetryCost(0,".1"),RetryCost(0,0))
    peak=subcritical_ancestral_peak(costly,".95","1.5")
    if peak["kind"]!="interior":
        raise AssertionError("Reference cost window absent")
    reference={}
    for label,mdl,t in [
        ("no_cost_age_zero",no,0.),
        ("no_cost_age_peak",no,peak["peak_age"]),
        ("costly_age_zero",costly,0.),
        ("costly_age_peak",costly,peak["peak_age"]),
        ("costly_age_five",costly,5.),
        ("costly_age_twenty",costly,20.),
        ("below_threshold_age_zero",below,0.),
        ("below_threshold_age_peak",below,peak["peak_age"]),
    ]:
        reference[label]=after_silence(mdl,t,.95,1.5)
        reference[label]["alternative_ever_not_preempted"]=(
            alternative_ever_without_success_preemption(mdl,t,.95,1.5))
    if not (reference["costly_age_peak"]["any_nonextinct"]>
            reference["costly_age_zero"]["any_nonextinct"] and
            reference["no_cost_age_peak"]["any_nonextinct"]<
            reference["no_cost_age_zero"]["any_nonextinct"]):
        raise AssertionError("Cost and zero-cost tests disagree")
    grid=cost_grid(rates,rates.loss)
    counts=Counter()
    for x in grid:
        for label in ("cost_early_delay_raises_rescue","cost_early_delay_lowers_rescue"):
            counts[label]+=int(x[label])
    return ({
        "schema":"anto_afterlife_costly_retries_v15",
        "analysis_date":"2026-10-23",
        "synthetic_only":True,
        "plant_reversion_hazard_measured":False,
        "underlying_R_A_R_B_fitted_from_plant_data":False,
        "cost_rate_units":"dimensionless probability of modeled shared functional capacity destruction per failed expression attempt",
        "resource_status":"shared physiological/developmental reproductive capacity; not heritable pigment gene deletion",
        "previous_model":"anto v14 cost-free retry, with no effect of expression failures on capacity",
        "expression_attempt_hazards":{"A":.4,"B":.2},
        "pre_opportunity_resource_loss":{"C":.02,"A":.30,"B":.01},
        "offspring_means_Poisson_illustrative":{"A":.95,"B":1.5},
        "failure_cost_proposal":{"A_failure_route_cost":0.,
                                 "A_failure_shared_cost":.3,
                                 "B_failure_route_cost":0.,
                                 "B_failure_shared_cost":0.},
        "cost_window":peak,
        "reference":reference,
        "sensitivity_grid":{"settings":len(grid),"counts":dict(counts),
            "parameter_frequencies_not_estimated_from_nature":True},
        "source_sha256":hashlib.sha256(
            (ROOT/"src/anto_afterlife_costly_retries_v15.py").read_bytes()).hexdigest(),
        "source_biological_conclusions":{
            "new_Cirsium_flower_pigment_expression_assays":0,
            "Cirsium_arvense_2022_615_head_seed_stage_data_in_aza3":
                "Descriptive seed-conversion contrast; not a failed colour-expression functional cost",
            "Petunia_DFR_Cas9_2026_DOI":"10.1007/s00299-026-03903-9",
            "Petunia_pathway_expression_trait_evo_2023_url":
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC10320354/",
            "external_2026_paper_is_not_rerun_or_model_calibration":True,
            "sources_from_different_species_do_not_constitute_cross_clade_model_validation":True,
        },
        "important_adverse_controls":{
            "zero_cost_recovers_v14":True,
            "subthreshold_shared_cost_has_no_peak":True,
            "costly_retries_never_more_effective_at_same_silent_age_than_zero_cost":True,
            "selected_parameter_window_does_not_imply_real_evolvability_increase":True,
            "natural_genetic_rescue_fitness_allele_fixation_or_reversal_hazard_unmeasured":True,
        },
        "environment":{"python":platform.python_version()},
        "next_experiment":"Measure functional floral tissue state, failed-expression-associated physiological resource loss and repeated opportunities with lineage-level descendants; avoid treating seed yield as genotype survival."
    },grid)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--grid",type=Path)
    args=parser.parse_args()
    output,grid=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    if args.grid:
        args.grid.parent.mkdir(parents=True,exist_ok=True)
        with args.grid.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(grid[0]))
            w.writeheader()
            w.writerows(grid)
    print("ANTO_V15_FAILED_ATTEMPT_COST_MODEL",
          json.dumps({"peak":output["cost_window"],
                      "source_results":{name:round(row["any_nonextinct"],10)
                            for name,row in output["reference"].items()},
                      "grid":output["sensitivity_grid"],
                      "plant_regain_measured":False},sort_keys=True))


if __name__=="__main__":
    main()
