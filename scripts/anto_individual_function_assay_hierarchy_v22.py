"""ANTO v22: biological-individual versus technical-repeat assay inference.

The Iochrominae Ho–Smith 2016 7-unpigmented / 15-pigmented taxon
original-source registry contains ZERO independently validated complete
floral functional-competence assay results in ANTO v19–v21.

All numerical inputs here are DECLARED SYNTHETIC calibration assumptions.
They must not be represented as estimates of pollen selection, retained
ancestral pigment function, genotype frequencies, trait reacquisition or
transgenerational fixation.

Within ONE hypothetical sampled taxon/population:
    Z_i ~ Bernoulli(theta), independently across n biological individuals.
Each individual gets r technical measurements. A single shared handling
failure C_i ~ Bernoulli(c) makes all r results falsely positive if Z_i=0,
or falsely negative if Z_i=1. Otherwise r technical errors are conditionally
independent with sensitivity se and specificity sp. Across independently
handled individuals (Z_i,C_i) pairs are independent.

The common-mode assumption is a deliberate sensitivity contrast, not a
calibrated laboratory failure model. An assay batch shared ACROSS plants,
family relatedness, site/genotype correlations, developmental plasticity,
nonrepresentative specimen selection, and incomplete gold-standard pathways
all violate this toy independence assumption.

A taxon is said to have a *positive screen* if at least one individual has
at least k positive technical calls out of r. That positive screen does
NOT establish species-level fixation, ancestral allele function, re-expression
after historical loss, or long-term population viability.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb
from pathlib import Path

from anto_fallible_function_assays_v21 import source_panel_no_assays

ROOT=Path(__file__).resolve().parents[1]


def _rate(x,name,*,allow_zero=True,allow_one=False):
    if type(x) is bool:
        raise ValueError(f"{name} cannot be a boolean")
    try:
        value=F(str(x)) if isinstance(x,float) else F(x)
    except (ValueError,TypeError,ZeroDivisionError,OverflowError) as exc:
        raise ValueError(f"{name} must be a finite numeric rational probability") from exc
    if not ((0 <= value if allow_zero else 0 < value) and
            (value <= 1 if allow_one else value < 1)):
        raise ValueError(f"{name} outside allowed probability range")
    return value


def _positive_int(v,name):
    if type(v) is not int or v<1:
        raise ValueError(f"{name} must be a positive integer")
    return v


@dataclass(frozen=True)
class HypotheticalAssay:
    sensitivity:F
    specificity:F
    shared_plant_handling_failure:F=F(0)

    def __post_init__(self):
        object.__setattr__(self,"sensitivity",
                           _rate(self.sensitivity,"sensitivity",allow_zero=False))
        object.__setattr__(self,"specificity",
                           _rate(self.specificity,"specificity",allow_zero=False))
        object.__setattr__(self,"shared_plant_handling_failure",
                           _rate(self.shared_plant_handling_failure,"shared_plant_handling_failure"))

    @property
    def synthetic_and_not_validated(self):
        return True


def binomial_tail(repetitions:int,minimum_positives:int,probability)->F:
    """Independent technical assay positives, conditional on genotype."""
    r=_positive_int(repetitions,"repetitions")
    k=_positive_int(minimum_positives,"minimum_positives")
    if k>r:
        raise ValueError("minimum_positives exceeds technical repetitions")
    p=F(probability)
    if p<0 or p>1:
        raise ValueError("Binomial success probability must lie in [0,1]")
    return sum(F(comb(r,i))*p**i*(1-p)**(r-i) for i in range(k,r+1))


def individual_call_probabilities(assay:HypotheticalAssay,repetitions:int,minimum_positives:int):
    """Conditional screen-positive probabilities given Z=1 or Z=0."""
    r=_positive_int(repetitions,"repetitions")
    k=_positive_int(minimum_positives,"minimum_positives")
    if k>r:
        raise ValueError("Positive-call threshold exceeds repetitions")
    c=assay.shared_plant_handling_failure
    # Shared plant-wise failure: a competent specimen's calls are all 0;
    # an incompetent specimen's calls are all 1. Otherwise independent.
    tpr=(1-c)*binomial_tail(r,k,assay.sensitivity)
    fpr=c+(1-c)*binomial_tail(r,k,1-assay.specificity)
    return tpr,fpr


def one_taxon_screen(n_individuals:int,repetitions:int,minimum_positives:int,
                     within_taxon_competence_fraction,assay:HypotheticalAssay)->dict:
    """Exact finite synthetic probability of a taxon screen and its errors.

    The first-order units are n independent biological individuals, NOT n*r
    independently inherited gene lineages.
    """
    n=_positive_int(n_individuals,"n_individuals")
    r=_positive_int(repetitions,"repetitions")
    k=_positive_int(minimum_positives,"minimum_positives")
    if k>r:
        raise ValueError("Positive-call threshold exceeds repetitions")
    theta=_rate(within_taxon_competence_fraction,"within_taxon_competence_fraction",allow_one=True)
    tpr,fpr=individual_call_probabilities(assay,r,k)
    marginal=theta*tpr+(1-theta)*fpr
    any_screen=1-(1-marginal)**n
    true_positive_screen=1-(1-theta*tpr)**n
    zero_competent=(1-theta)**n
    false_screen_if_none=1-(1-fpr)**n
    return {
        "n_biologically_distinct_plants":n,
        "technical_repeats_per_plant":r,
        "minimum_positive_repeats_for_individual_call":k,
        "total_assay_calls":n*r,
        "hypothetical_within_taxon_competence_fraction":theta,
        "plant_positive_call_probability_if_competent":tpr,
        "plant_positive_call_probability_if_incompetent":fpr,
        "chance_at_least_one_plant_screen_positive":any_screen,
        "chance_at_least_one_TRULY_competent_plant_correctly_called":true_positive_screen,
        "chance_zero_competent_plants_sampled":zero_competent,
        "false_taxon_screen_if_all_plants_incompetent":false_screen_if_none,
        "not_taxon_genetic_fixation_or_multi_generation_regain":True,
    }


def same_marginal_nonidentification():
    """One raw call per plant cannot identify prevalence AND test error.

    Both distinct hypothetical triples yield the exact same P(Y=1)=.22.
    This is a strict negative identifiability example, not fitting to
    22 Iochrominae species or to 2025 backcross expression records.
    """
    a=(F(1,5),F(9,10),F(19,20))
    b=(F(1,4),F(17,20),F(99,100))
    def q(t):
        theta,se,sp=t
        return theta*se+(1-theta)*(1-sp)
    if q(a)!=F(11,50) or q(b)!=F(11,50):
        raise ArithmeticError("Identifiability witness lost")
    return {
        "hypothetical_world_A":a,
        "hypothetical_world_B":b,
        "same_observable_single_binary_call_probability":q(a),
        "different_true_within_taxon_prevalence":True,
        "not_observed_Iochrominae_assay_calibration":True,
    }


def finite_budget_search(total_call_budget:int,theta,assay:HypotheticalAssay,
                         maximum_false_taxon_screen: F)->dict:
    """Enumerate uniform (n biological plants, r technical repeats, k cutoff).

    Objective: maximize P(at least one TRUE competent plant correctly
    screened) among designs satisfying a declared FPR constraint at theta=0.
    Uses <= budget technical assay calls. This is NOT a statistical power
    analysis, expected information gain, phylogenetic design or biological
    sampling recommendation; theta,se,sp,c are all stipulated synthetic.
    """
    budget=_positive_int(total_call_budget,"total_call_budget")
    theta=_rate(theta,"theta",allow_one=True)
    ceiling=F(maximum_false_taxon_screen)
    if not 0<=ceiling<=1:
        raise ValueError("False screen ceiling must lie in [0,1]")
    eligible=[]
    all_count=0
    for n in range(1,budget+1):
        for r in range(1,budget//n+1):
            for k in range(1,r+1):
                all_count+=1
                row=one_taxon_screen(n,r,k,theta,assay)
                if row["false_taxon_screen_if_all_plants_incompetent"]<=ceiling:
                    eligible.append(row)
    eligible.sort(key=lambda x:(
        -x["chance_at_least_one_TRULY_competent_plant_correctly_called"],
        x["false_taxon_screen_if_all_plants_incompetent"],
        x["total_assay_calls"],
        x["n_biologically_distinct_plants"],
        x["technical_repeats_per_plant"],
        x["minimum_positive_repeats_for_individual_call"],
    ))
    return {
        "total_call_budget":budget,
        "tested_synthetic_candidate_designs":all_count,
        "designs_passing_declared_false_screen_ceiling":len(eligible),
        "minimum_true_competence_frequency_assumed":None,
        "assumed_within_taxon_competence_frequency":theta,
        "assumed_false_taxon_screen_ceiling":ceiling,
        "assumed_plant_handling_shared_failure":assay.shared_plant_handling_failure,
        "optimal_under_THIS_STIPULATED_likelihood_only":eligible[0] if eligible else None,
        "not_real_sample_size_power_or_estimated_field_error":True,
    }


def report():
    source=source_panel_no_assays()
    assay=HypotheticalAssay(F(9,10),F(19,20))
    common=HypotheticalAssay(F(9,10),F(19,20),F(1,50))
    high_common=HypotheticalAssay(F(9,10),F(19,20),F(1,20))
    theta=F(1,5)
    benchmarks={
        "four_independent_plants_one_call":one_taxon_screen(4,1,1,theta,assay),
        "one_plant_four_technical_repeats_any_positive":one_taxon_screen(1,4,1,theta,assay),
        "one_plant_four_technical_repeats_two_positive":one_taxon_screen(1,4,2,theta,assay),
        "six_plants_two_technical_repeats_both_positive":one_taxon_screen(6,2,2,theta,assay),
    }
    designs={
        "independent_replicate_error_assumption":finite_budget_search(
            12,theta,assay,F(1,20)),
        "two_percent_per_plant_shared_handling_failure":finite_budget_search(
            12,theta,common,F(1,20)),
        "five_percent_per_plant_shared_handling_failure":finite_budget_search(
            12,theta,high_common,F(1,20)),
    }
    def spec(row):
        if row is None: return None
        return [row["n_biologically_distinct_plants"],
                row["technical_repeats_per_plant"],
                row["minimum_positive_repeats_for_individual_call"]]
    if (spec(designs["independent_replicate_error_assumption"]
             ["optimal_under_THIS_STIPULATED_likelihood_only"])!=[6,2,2]
        or spec(designs["two_percent_per_plant_shared_handling_failure"]
             ["optimal_under_THIS_STIPULATED_likelihood_only"])!=[2,6,3]
        or designs["five_percent_per_plant_shared_handling_failure"]
             ["optimal_under_THIS_STIPULATED_likelihood_only"] is not None):
        raise ArithmeticError("Declared synthetic assay-design sensitivity witness changed")
    # Same call-count and any-one-positive rule: biological rather than
    # technical replication alters TRUE-plant detection.
    a=benchmarks["four_independent_plants_one_call"]
    b=benchmarks["one_plant_four_technical_repeats_any_positive"]
    if (a["false_taxon_screen_if_all_plants_incompetent"]!=
        b["false_taxon_screen_if_all_plants_incompetent"] or
        a["chance_at_least_one_TRULY_competent_plant_correctly_called"]<=
        b["chance_at_least_one_TRULY_competent_plant_correctly_called"]):
        raise ArithmeticError("Original four-test distinct-individual witness lost")
    return {
        "schema":"anto_v22_within_taxon_biological_vs_technical_replication",
        "date":"2026-10-26",
        "empirical_source_frame":"Ho & Smith 2016 22 original Iochrominae taxa via ANTO v19-v21",
        "observed_original_source_taxa":len(source),
        "observed_original_source_unpigmented_taxa":sum(
            x["flower_pigment_category"]=="UNPIGMENTED" for x in source),
        "observed_original_source_pigmented_taxa":sum(
            x["flower_pigment_category"]=="PIGMENTED" for x in source),
        "actual_qualified_complete_flower_function_assays":0,
        "source_DFR_CDS_annotation_is_not_a_function_assay":True,
        "input_rates_are_hypothetical_not_calibrated":True,
        "within_taxon_frequency_across_source_clades_is_not_estimated":True,
        "one_sampled_population_not_an_entire_species_or_ancestral_transition":True,
        "technical_repeat_independence_requires_external_assay_validation":True,
        "real_field_shared_batch_and_plant_relatedness_not_covered":True,
        "no_natural_floral_trait_regain_heritability_or_fitness_measurement":True,
        "mechanistic_single_call_nonidentifiability":same_marginal_nonidentification(),
        "exact_four_call_biological_vs_technical_benchmarks":benchmarks,
        "hypothetical_fixed_12_call_budget_search_under_declared_assumptions":designs,
        "Cirsium_organismal_empirical_owner":"zuizui0223/aza3",
    }


def serialize(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {str(k):serialize(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [serialize(v) for v in value]
    return value


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    data=serialize(report())
    if args.expected and read_expected(args.expected)!=data:
        raise ValueError("Frozen original-source v22 no-assay experiment-design receipt mismatch")
    s=json.dumps(data,sort_keys=True,indent=2,ensure_ascii=False)+"\n"
    if args.output is None:
        print(s)
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding="utf-8")


def read_expected(path):
    return json.loads(path.read_text(encoding="utf-8"))


if __name__=="__main__":
    main()
