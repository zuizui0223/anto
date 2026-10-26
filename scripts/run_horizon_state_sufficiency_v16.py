"""Reproduce exact ANTO v16 horizon/endpoint sufficiency contrasts.

All numeric inputs are stipulated synthetic rational probabilities.
No species-specific afterlife rate, genetic regain, or demographic
persistence rate is estimated or calibrated from field data.
"""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_horizon_state_sufficiency_v16 import (
    exact_reference_panel, binomial_two_nonextinction,
)


def serialize(obj):
    if isinstance(obj,Fraction):
        return str(obj)
    if isinstance(obj,dict):
        return {str(k):serialize(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):
        return [serialize(x) for x in obj]
    return obj


def receipt():
    return serialize({
        "version":"v16",
        "study_type":"exact synthetic finite-horizon state-sufficiency and branching filter",
        "same_present_visible_phenotype":"white in both constructed histories",
        "same_present_functional_inventory":"all three required modules stipulated complete",
        "conditional_branching_nonextinction":binomial_two_nonextinction(Fraction(2,3)),
        "full_future_transition_kernel_is_part_of_sufficient_state":True,
        "short_horizon_acceptance_establishes_long_horizon_sufficiency":False,
        "no_empirical_plant_parameter_calibrated":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "wild_trait_regain_observed":False,
        "permanent_genetic_fixation_inferred":False,
        "ecological_population_rescue_measured":False,
        "results":exact_reference_panel(),
    })


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    s=json.dumps(receipt(),sort_keys=True,indent=2,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding="utf-8")
    else:
        print(s)


if __name__=="__main__":main()
