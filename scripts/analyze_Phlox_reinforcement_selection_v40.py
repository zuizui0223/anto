"""ANTO v40: exact published Phlox reinforcing selection equation + conditional costs.

Original Hopkins & Rausher (2012) aggregate hybridization model, and separate
(2014) allopatric visitor preference (original author bootstrap CIs).
Only deterministic mathematical counterfactuals are newly calculated.
No raw HPLC, pollen paternity, genotype selection time series or 2014 Dryad
individual data are remeasured or downloaded here.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/external/anto_v40_original_Phlox_selection_2012_2014_author_registry_20261116.json"
SOURCE_IDS={
    "phlox_2012_sympatric_reinforcement":"10.1126/science.1215198",
    "phlox_2014_allopatric_pollinator_cost":"10.1086/675495",
}
Q_GRID=(0,.25,.5,.75,.9,1)
COST_GRID=(0,.05,.1,.15,.2,.25,.3)


def require(ok,code):
    if not ok:
        raise ValueError(code)


def prob(value):
    if (isinstance(value,bool) or not isinstance(value,(float,int))
            or not math.isfinite(value) or not 0<=value<=1):
        raise ValueError("INVALID_PARAMETER_OUTSIDE_UNIT_INTERVAL")
    return float(value)


def viable_fraction(h,q):
    return 1-prob(q)*prob(h)


def hybrid_advantage(h_light,h_dark,q):
    w_light=viable_fraction(h_light,q)
    w_dark=viable_fraction(h_dark,q)
    if w_light<=0:
        raise ValueError("ZERO_REFERENCE_VIABLE_FITNESS")
    return w_dark/w_light-1


def net_selection(h_light,h_dark,q,nonhybrid_dark_cost):
    cost=prob(nonhybrid_dark_cost)
    w_light=viable_fraction(h_light,q)
    w_dark=viable_fraction(h_dark,q)
    if w_light<=0:
        raise ValueError("ZERO_REFERENCE_VIABLE_FITNESS")
    return (1-cost)*w_dark/w_light-1


def break_even_nonhybrid_cost(h_light,h_dark,q):
    w_light=viable_fraction(h_light,q)
    w_dark=viable_fraction(h_dark,q)
    if w_dark<=0:
        raise ValueError("ZERO_DARK_VIABLE_FITNESS")
    if w_light<=0:
        raise ValueError("ZERO_REFERENCE_VIABLE_FITNESS")
    return 1-w_light/w_dark


def required_sterility_for_cost(h_light,h_dark,nonhybrid_dark_cost):
    hL,hD,c=prob(h_light),prob(h_dark),prob(nonhybrid_dark_cost)
    if c==0:
        return 0.0
    denom=hL-(1-c)*hD
    if denom<=0:
        return None
    return c/denom


def _rounded(x):
    return None if x is None else round(x,12)


def validate_registry(source):
    require(source.get("schema")=="anto_v40_Phlox_original_reinforcement_assumptions_and_allopatry_source_registry",
            "SOURCE_SCHEMA_DRIFT")
    papers=source.get("source_papers")
    require(isinstance(papers,list) and len(papers)==2 and
            len({x.get("id") for x in papers})==2,
            "INDEPENDENT_ORIGINAL_SOURCE_SET_DRIFT")
    d={x["id"]:x for x in papers}
    require({k:x["doi"] for k,x in d.items()}==SOURCE_IDS,
            "SOURCE_ORIGINAL_DOI_DRIFT")
    require(source.get("independent_original_publication_count")==2 and
            source.get("one_shared_Phlox_drumondii_phylogenetic_radiation") is True
            and all(x["species"]=="Phlox drummondii" and
                    x["congener"]=="Phlox cuspidata" for x in papers),
            "ORIGINAL_SHARED_RADIATION_SOURCE_DRIFT")
    require(source.get("fixed_scenario_q_grid")==list(Q_GRID)
            and source.get("fixed_scenario_nonhybrid_dark_cost_grid")==list(COST_GRID),
            "SCENARIO_SOURCE_PREDECLARATION_DRIFT")
    require(all(source.get(k) is True for k in (
        "no_uncertainty_intervals_computed",
        "cost_is_hypothetical_not_measured_dark_allele_penalty",
        "no_independent_RNA_HPLC_or_natural_regain_measured",
        "no_new_2014_original_Dryad_raw_data_acquired",
        "no_numerical_meta_analysis_between_the_two_studies",
        "no_2012_and_2014_independent_radiation_claim"
    )),"SCIENTIFIC_GUARD_DRIFT")
    p=d["phlox_2012_sympatric_reinforcement"]
    a=p["original_measurements"]
    require(p["year"]==2012 and p["journal"]=="Science"
            and p["experimental_context"]==
                "FOCAL_PATERNITY_BLOCKS_WITH_P_CUSPIDATA_COEXISTENCE"
            and p["original_author_PDF_url"]==
                "https://hopkins-lab.org/wp-content/uploads/2014/12/Hopkins-and-Rausher_2012.pdf"
            and p["original_Science_supporting_material"]=="science.1215198.DC1",
            "SOURCE_2012_BIBLIOGRAPHY_DRIFT")
    require(a["light_reference_hybrid_paternity_fraction"]==.43
            and a["dark_focal_hybrid_paternity_fraction"]==.21
            and a["hybrid_sterility_assumed_from_earlier_studies"]==.9
            and a["original_selection_coefficient_rounded"]==.32
            and a["reference_hybridization_four_blocks_range"]==[.28,.44],
            "SOURCE_2012_PARAMETRIC_INPUT_DRIFT")
    require(a["paternity_figure_seeds_genotyped"]==931
            and a["original_field_paternity_blocks"]==4
            and a["separate_pollinator_individuals_observed"]==181
            and a["separate_interplant_pollinator_transitions"]==2301
            and a["same_observed_visitation_and_paternity_individuals_joined"] is False
            and a["first_paper_s_is_author_model_estimate_not_observed_allele_frequency_change"] is True,
            "SOURCE_2012_POLLEN_AND_VISITOR_GRAIN_DRIFT")
    require(a["original_author_relative_seed_viability_formula"]=="W=1-q*h"
            and a["original_author_s_definition"]=="(1-q*h_dark)/(1-q*h_light)-1"
            and a["dark_intensity_locus"]=="R2R3-MYB_UPREGULATED_INTENSITY"
            and a["observed_dark_reduces_interspecific_hybridization_more_than_half"] is True
            and a["source_maternal_fitness_in_congener_absent_common_garden_difference_detected"] is False,
            "SOURCE_2012_REGULATORY_INTENSITY_AND_MECHANISM_DRIFT")
    require(a["hue_red_F3prime5primeH_locus_hybridization_reduction_detected"] is False,
            "HUE_LOCUS_NULL_GUARD_DRIFT")

    s=d["phlox_2014_allopatric_pollinator_cost"]
    b=s["original_measurements"]
    require(s["year"]==2014 and s["journal"]=="The American Naturalist"
            and s["original_Dryad_DOI"]=="10.5061/dryad.f710n"
            and s["original_author_PDF_url"]==
                "https://hopkins-lab.org/wp-content/uploads/2015/03/HopkinsRausher_2014.pdf"
            and s["experimental_context"]==
                "ALLOPATRIC_NATIVE_P_DRUMMONDII_FLOWER_COLOUR_2010_2011_FIELD_AND_POLLINATORS",
            "SOURCE_2014_BIBLIOGRAPHY_DRIFT")
    require(b["original_field_years"]==[2010,2011]
            and b["total_observed_pollinators_all_arrays"]==278
            and b["total_observed_visits_all_arrays"]==1000
            and b["array_n_sum_can_differ_from_original_total_278"] is True
            and b["all_three_pollinator_comparisons_ancestral_light_blue_preferred"] is True
            and b["year2010_2011_maternal_fitness_point_estimates_and_weather_differ"] is True
            and b["year2011_derived_allele_net_invasion_depends_on_fecundity_definition_and_male_assumptions"] is True
            and b["original_field_paternity_sire_genotypes_not_directly_measured"] is True
            and b["original_author_CI_are_clustered_bootstrap_not_simple_binomial"] is True
            and b["not_same_2012_selection_model_or_original_2012_field_population"] is True
            and b["source_underlying_2014_Dryad_rows_reanalyzed_in_v40"] is False,
            "SOURCE_2014_OBSERVATION_GRAIN_DRIFT")
    expected=[
        ("light_red",84,.619,[.571,.666]),
        ("dark_blue",118,.603,[.553,.653]),
        ("dark_red",74,.682,[.634,.732]),
    ]
    arr=b["arrays"]
    require(isinstance(arr,list) and len(arr)==3 and
            all([r["focal_derived_colour"],r["original_pollinators_in_table4"],
                 r["share_visits_to_ancestral_light_blue"],r["author_bootstrap_95ci"]]
                ==list(record) for r,record in zip(arr,expected)),
            "ORIGINAL_2014_TABLE4_DRIFT")
    return a,b


def review(source):
    a,b=validate_registry(source)
    hL=a["light_reference_hybrid_paternity_fraction"]
    hD=a["dark_focal_hybrid_paternity_fraction"]
    original_q=a["hybrid_sterility_assumed_from_earlier_studies"]
    W_L=viable_fraction(hL,original_q)
    W_D=viable_fraction(hD,original_q)
    s=hybrid_advantage(hL,hD,original_q)
    require(round(s,2)==a["original_selection_coefficient_rounded"],
            "ORIGINAL_2012_S_SELECTION_REPLAY_DRIFT")
    cstar=break_even_nonhybrid_cost(hL,hD,original_q)

    q_results=[]
    scenarios=[]
    for q in Q_GRID:
        q_results.append({
            "hybrid_sterility_fraction":q,
            "modelled_light_viability_proxy":viable_fraction(hL,q),
            "modelled_dark_viability_proxy":viable_fraction(hD,q),
            "modelled_selection":hybrid_advantage(hL,hD,q),
            "original_study_estimate_only_if_q_is_0_90":q==original_q,
            "not_a_new_original_study_observation":True,
        })
        for c in COST_GRID:
            value=net_selection(hL,hD,q,c)
            scenarios.append({
                "hybrid_sterility_fraction":q,
                "hypothetical_dark_nonhybrid_cost_fraction":c,
                "modelled_light_viability_proxy":viable_fraction(hL,q),
                "modelled_dark_viability_after_hypothetical_cost":
                    (1-c)*viable_fraction(hD,q),
                "relative_net_selection":value,
                "dark_intensity_allele_modelled_advantage_sign":
                    "POSITIVE" if value>1e-12 else
                    "NEGATIVE" if value< -1e-12 else "ZERO",
                "not_empirical_confidence_interval":True,
                "not_an_independent_genotype_experiment":True,
            })

    cost_thresholds=[]
    for c in COST_GRID:
        q0=required_sterility_for_cost(hL,hD,c)
        cost_thresholds.append({
            "hypothetical_dark_nonhybrid_cost_fraction":c,
            "required_q_for_positive_dark_net_selection":q0,
            "threshold_admissible_q_LE_1":q0 is not None and q0<=1,
            "fixed_original_q_0_90_net_dark_allele_favored":
                net_selection(hL,hD,original_q,c)>1e-12,
            "not_a_measured_cost_or_empirical_site_sterility":True,
        })
    allopatry=[]
    for row in b["arrays"]:
        share=row["share_visits_to_ancestral_light_blue"]
        ci=row["author_bootstrap_95ci"]
        allopatry.append({
            "derived_flower_colour":row["focal_derived_colour"],
            "author_original_array_N":row["original_pollinators_in_table4"],
            "light_blue_visit_share":share,
            "other_derived_colour_visit_share_complement":1-share,
            "light_blue_vs_focal_visit_share_ratio":share/(1-share),
            "author_bootstrap_95CI":ci,
            "author_reported_bootstrap_CI_excludes_50pct":ci[0]>.5,
            "NOT_same_as_paternal_fitness_or_paternity_hybridization_rate":True,
            "source_individual_visits_not_all_independent":True,
        })
    return {
        "schema":"anto_v40_original_Phlox_conditional_reinforcement_and_nonhybrid_cost_envelope",
        "original_source_registry":"data/external/anto_v40_original_Phlox_selection_2012_2014_author_registry_20261116.json",
        "n_original_published_papers":source["independent_original_publication_count"],
        "n_independent_phylogenetic_radiations":1,
        "two_papers_counted_as_independent_phylogenetic_reacquisitions":False,
        "source_DOIs":["10.1126/science.1215198","10.1086/675495"],
        "original_paternity_figure_seeds_genotyped":931,
        "original_2012_pollinator_individuals_separate_from_paternity":181,
        "original_2012_interplant_transitions_separate_from_paternity":2301,
        "reproduced_original_2012_author_selection":{
            "reference_light_hybrid_paternity_fraction":hL,
            "dark_hybrid_paternity_fraction":hD,
            "hybrid_sterility_fraction_original_assumption":original_q,
            "light_viable_fitness_proxy":W_L,
            "dark_viable_fitness_proxy":W_D,
            "source_relative_selection_coefficient":s,
            "source_selection_rounded_by_original_authors":a["original_selection_coefficient_rounded"],
            "hybridization_reduction_fraction":(hL-hD)/hL,
            "actual_net_population_allele_frequency_increase_measured":False,
            "model_not_directly_observed_32pct_seed_number_gain":True
        },
        "source_parameter_q_sensitivity":q_results,
        "q_by_hypothetical_dark_nonhybrid_cost_grid":scenarios,
        "hypothetical_break_even_dark_nonhybrid_cost_for_q_0_90":{
            "cost_fraction":cstar,
            "relative_cost_percent":100*cstar,
            "net_dark_fitness_advantage_at_cost":net_selection(hL,hD,original_q,cstar),
            "not_a_measured_actual_metabolic_or_mating_cost":True
        },
        "q_required_by_cost_scenarios":cost_thresholds,
        "independent_2014_author_allopatry_pollinator_preferences":{
            "original_DOI":"10.1086/675495",
            "total_original_pollinators_reported":b["total_observed_pollinators_all_arrays"],
            "total_original_visits_reported":b["total_observed_visits_all_arrays"],
            "sum_of_three_table4_array_N":sum(x["original_pollinators_in_table4"]
                                               for x in b["arrays"]),
            "author_CIs_not_recalculated_from_independent_binomial_visits":True,
            "arrays":allopatry,
            "all_three_original_table4_light_blue_bootstrap_intervals_above_0_5":
                all(x["author_reported_bootstrap_CI_excludes_50pct"] for x in allopatry),
            "2010_and_2011_maternal_component_estimates_and_fecundity_model_matter":True,
            "same_2012_paternity_blocks_not_used":True,
        },
        "original_2012_hue_red_locus_no_detected_hybridization_reduction":True,
        "original_2012_derived_s_transferred_to_red_hue_allele":False,
        "allopatry_pollinator_visit_share_is_hybridization_rate":False,
        "2014_allopatry_direct_dark_paternal_siring_success_measured":False,
        "2014_allopatry_net_colour_selection_same_as_2012_s":False,
        "2014_2010_2011_net_predictions_depend_on_fecundity_model":True,
        "2014_Dryad_original_individual_rows_reanalyzed":False,
        "2012_paternity_931_seed_rows_reanalyzed":False,
        "numerical_sensitivity_scenarios_not_independent_observed_genotypes":True,
        "original_cross_paper_quantitative_meta_selection_effect":None,
        "original_2014_allopatric_dark_allele_cost_empirically_measured":None,
        "observed_allele_fixation_rate":None,
        "no_new_species_phylogenetic_replication_or_direct_colour_regain":True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--expected",type=Path)
    args=ap.parse_args()
    result=review(json.loads(SOURCE.read_text(encoding="utf-8")))
    if args.expected and json.loads(args.expected.read_text(encoding="utf-8"))!=result:
        raise ValueError("V40_FROZEN_ORIGINAL_SOURCE_AND_MATHEMATICAL_SCENARIO_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    important={
        "author_formula_s":result["reproduced_original_2012_author_selection"]["source_relative_selection_coefficient"],
        "break_even_nonhybrid_cost":result["hypothetical_break_even_dark_nonhybrid_cost_for_q_0_90"]["cost_fraction"],
        "s_by_sterility":{str(x["hybrid_sterility_fraction"]):x["modelled_selection"]
                         for x in result["source_parameter_q_sensitivity"]},
        "q_required_at_costs":{str(x["hypothetical_dark_nonhybrid_cost_fraction"]):
                               x["required_q_for_positive_dark_net_selection"]
                               for x in result["q_required_by_cost_scenarios"]},
        "original_allopatry_table4_light_blue_shares":{
            x["derived_flower_colour"]:x["light_blue_visit_share"]
            for x in result["independent_2014_author_allopatry_pollinator_preferences"]["arrays"]},
        "no_empirical_dark_cost_or_fixation_rate":True,
    }
    print("ANTO_V40_SOURCE_CONDITIONAL_FITNESS_RESULT",json.dumps(important,sort_keys=True))
    packed=json.dumps(result,sort_keys=True,separators=(",",":"))
    stride=6000
    n=(len(packed)+stride-1)//stride
    for i in range(n):
        print("ANTO_V40_DERIVED_CHUNK",i+1,n,packed[i*stride:(i+1)*stride])


if __name__=="__main__":main()
