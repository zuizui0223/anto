"""ANTO v37: original-published natural regain, native flower rescue, pollinator choice.

No genomic sequences, HPLC rows, experimentally transformed plants or insect
trials are downloaded or re-analyzed in this audit. Distinguishes original
publication, shared radiation, native white petal restoration, and fitness.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/external/anto_v37_natural_regain_and_native_reactivation_original_sources_20261116.json"
DOIS={
    "Esfeld2018_Psecreta":"10.1016/j.cub.2018.10.019",
    "Berardi2021_Pexserta":"10.1093/plcell/koab114",
    "Hoballah2007_AN2_pollinators":"10.1105/tpc.106.048694",
    "Fraser2013_Actinidia_MYB110a":"10.1186/1471-2164-14-28",
}
RADIATION={
    "Esfeld2018_Psecreta":"Petunia_long_tube_clade",
    "Berardi2021_Pexserta":"Petunia_long_tube_clade",
    "Hoballah2007_AN2_pollinators":"Petunia_long_tube_clade",
    "Fraser2013_Actinidia_MYB110a":"Actinidia_chinensis_eriantha_hybrid_backcross",
}
NATURAL_ROUTES={
    "Petunia secreta":"AN2_FRAME_RESTORATION",
    "Petunia exserta":"PARALOG_DPL_REDEPLOYMENT_WITH_OTHER_PATHWAY_CHANGES",
}


def require(condition,code):
    if not condition:
        raise ValueError(code)


def review(data):
    require(data.get("schema")=="anto_v37_original_petunia_regain_actinidia_reactivation_ecological_evidence_registry",
            "SOURCE_SCHEMA_DRIFT")
    studies=data.get("original_study_programs")
    require(isinstance(studies,list) and len(studies)==4 and
            len({x.get("id") for x in studies})==4,
            "ORIGINAL_STUDY_SET_DRIFT")
    by={x["id"]:x for x in studies}
    require({key:row["doi"] for key,row in by.items()}==DOIS,
            "ORIGINAL_DOI_DRIFT")
    require({key:row["shared_radiation"] for key,row in by.items()}==RADIATION,
            "SOURCE_SHARED_PHYLOGENETIC_RADIATION_DRIFT")
    require({x["family"] for x in studies}=={"Solanaceae","Actinidiaceae"}
            and all(x["one_independent_original_study"] is True for x in studies)
            and all(x["family"]=="Solanaceae" for x in studies[:3])
            and studies[3]["family"]=="Actinidiaceae",
            "ORIGINAL_FAMILY_GRAIN_DRIFT")
    require(all(data[k] is True for k in (
        "do_not_count_independent_Petunia_studies_as_independent_macro_radiations",
        "no_crossstudy_numeric_RMSE_or_regain_rate",
        "no_4_study_independent_evolutionary_event_claim",
        "no_native_ancestral_gene_resurrection_rate_estimated",
        "pollinator_visits_not_equivalent_to_fitness",
    )),"SOURCE_INFERENCE_GUARD_DRIFT")

    es=by["Esfeld2018_Psecreta"]
    an=es["natural_colour_regain"]
    require(an is not None and an["species"]=="Petunia secreta"
            and an["visible_colour"]=="PURPLE"
            and an["ancestor_reconstructed"]=="WHITE_NONFUNCTIONAL_AN2"
            and an["molecular_channel"]=="AN2_FRAME_RESTORATION"
            and an["function_restoration"]=="AN2_CODING_READING_FRAME"
            and an["indel"]=="2_BASE_PAIR_COMPENSATORY_DELETION",
            "SOURCE_AN2_REPAIR_MECHANISM_DRIFT")
    require(an["direct_historical_transition_observed_in_real_time"] is False
            and es["original_white_petal_transformation"] is None
            and es["controlled_pollinator_preference"] is None,
            "HISTORICAL_DIRECT_OBSERVATION_FORGED")

    be=by["Berardi2021_Pexserta"]
    d=be["natural_colour_regain"]
    require(d is not None and d["species"]=="Petunia exserta" and
            d["visible_colour"]=="RED"
            and d["ancestor_reconstructed"]=="WHITE_NONFUNCTIONAL_AN2"
            and d["molecular_channel"]==
            "PARALOG_DPL_REDEPLOYMENT_WITH_OTHER_PATHWAY_CHANGES"
            and d["used_regulator"]=="DPL" and d["AN2_still_nonfunctional"] is True
            and d["regulatory_spatial_shift"]=="DPL_IN_PETAL_LIMBS",
            "DPL_AND_AN2_STATE_DRIFT")
    require(set(d["additional_required_pathway_changes"])==
            {"MYB-FL","F3'H/HT1","F3'5'H/HF1","F3'5'H/HF2","AAT"}
            and d["original_red_chemical_components_not_pelargonidin"] is True
            and d["petal_red_anthocyanidins"]==["cyanidin","delphinidin"]
            and d["DPL_native_VIGS_white_sectors"] is True,
            "SOURCE_MULTIGENE_RED_PIGMENT_DRIFT")
    require(d["direct_historical_transition_observed_in_real_time"] is False,
            "HISTORICAL_DIRECT_OBSERVATION_FORGED")
    tr=be["original_white_petal_transformation"]
    require(tr is not None and tr["native_white_species"]=="Petunia axillaris"
            and tr["regulator"]=="DPL"
            and tr["design"]=="STABLE_35S_OVEREXPRESSION"
            and tr["generated_independent_transgenic_lines"]==9
            and tr["lines_selected_for_chemical_assay"]==2
            and tr["HPLC_nontransgenic_samples"]==3
            and tr["HPLC_transgenic_samples"]==6
            and tr["HPLC_transgenics_per_assayed_line"]==3
            and tr["HPLC_transgenic_samples"]==
                tr["lines_selected_for_chemical_assay"]*
                tr["HPLC_transgenics_per_assayed_line"],
            "DPL_TRANSFORMATION_UNIT_DRIFT")
    require(tr["white_to_coloured_corolla_restored"] is True
            and tr["anthocyanidin_increase_original_assay"] is True
            and tr["natural_ancestral_reversal"] is False
            and tr["other_plant_traits_assumed_unchanged"] is False
            and be["controlled_pollinator_preference"] is None,
            "DPL_NOT_A_WILD_REGAIN_ASSAY")

    ho=by["Hoballah2007_AN2_pollinators"]
    require(ho["natural_colour_regain"] is None,"FORGED_THIRD_NATURAL_REGAIN")
    an2=ho["original_white_petal_transformation"]
    require(an2 is not None and an2["native_white_species"]=="Petunia axillaris"
            and an2["donor_of_functional_regulator"]=="Petunia integrifolia"
            and an2["regulator"]=="AN2"
            and an2["design"]=="STABLE_TRANSFORMED_FUNCTIONAL_AN2_CDNA"
            and an2["white_to_coloured_corolla_restored"] is True
            and an2["pigment_functional_complementation"] is True
            and an2["not_spontaneous_pseudogene_resurrection"] is True,
            "AN2_NATIVE_WHITE_TRANSFORMATION_DRIFT")
    pol=ho["controlled_pollinator_preference"]
    require(pol is not None and pol["study_organisms"]==
            ["Manduca sexta","Bombus terrestris"] and
            pol["setup"]=="CONTROLLED_GREENHOUSE" and
            pol["moth_first_choice_white_P_axillaris"]==30 and
            pol["moth_first_choice_ANT2_restored_P_axillaris"]==6 and
            pol["moth_first_choice_total"]==36 and
            pol["moth_first_choice_white_P_axillaris"]+
            pol["moth_first_choice_ANT2_restored_P_axillaris"]==
            pol["moth_first_choice_total"] and
            pol["bee_n_trials"]==12 and pol["bee_direction"]==
            "PREFERENCE_AN2_RESTORED_COLOURED_P_AXILLARIS"
            and pol["moth_direction"]=="PREFERENCE_WILD_WHITE_P_AXILLARIS"
            and pol["moth_source_p_upper_bound"]==0.0001
            and pol["bee_author_reported_p"]==0.017,
            "POLLINATOR_SOURCE_DENOMINATOR_DRIFT")
    require(pol["field_pollination_or_seed_set_measured"] is False and
            pol["other_flower_traits_experimentally_proved_identical_everywhere"] is False
            and an2["pollinator_assay_same_plant_genotype_contrast"] is True,
            "PREFERENCE_CANNOT_BE_TREATED_AS_FITNESS")

    fr=by["Fraser2013_Actinidia_MYB110a"]
    require(fr["natural_colour_regain"] is None
            and fr["controlled_pollinator_preference"] is None,
            "ACTINIDIA_NOT_NATURAL_ANCESTRAL_GAIN")
    ax=fr["original_white_petal_transformation"]
    require(ax is not None
            and ax["native_white_material"].endswith("21-05-09e")
            and ax["regulator"]=="MYB110a"
            and ax["design"]=="TRANSIENT_AGROBACTERIUM_35S_EXPRESSION_IN_WHITE_PETAL_EDGE"
            and ax["white_to_coloured_corolla_restored"] is True
            and ax["pigmented_cells_local_to_infiltration"] is True
            and ax["control"]=="35S_GFP_NEGATIVE_NO_PIGMENT"
            and ax["long_term_stable_transgenic_plant_established"] is False
            and ax["natural_ancestral_reversal"] is False
            and ax["number_of_independent_Petunia_natural_color_regains_added"]==0,
            "ACTINIDIA_LOCAL_TRANSIENT_ONLY")

    natural=[s for s in studies if s["natural_colour_regain"] is not None]
    rescue=[s for s in studies if s["original_white_petal_transformation"] is not None]
    stable=[s for s in rescue if s["original_white_petal_transformation"]["design"].startswith("STABLE")]
    transient=[s for s in rescue if s not in stable]
    poll=[s for s in studies if s["controlled_pollinator_preference"] is not None]
    require(len(natural)==2 and len(rescue)==3 and len(stable)==2
            and len(transient)==1 and len(poll)==1,
            "SOURCE_EXPERIMENTAL_STUDY_COUNT_DRIFT")
    mapping={s["natural_colour_regain"]["species"]:s["natural_colour_regain"]["molecular_channel"]
             for s in natural}
    require(mapping==NATURAL_ROUTES,"ORIGINAL_REGAIN_CHANNEL_DRIFT")
    # Independent, strictly binomial arithmetic on the authors' already
    # published greenhouse FIRST-CHOICE counts only. Not a genetic/fecundity P.
    moth_n=pol["moth_first_choice_total"]
    moth_extreme=pol["moth_first_choice_white_P_axillaris"]
    binomial_p=2*math.fsum(math.comb(moth_n,k) for k in range(moth_extreme,moth_n+1))/(2**moth_n)
    require(binomial_p < pol["moth_source_p_upper_bound"],
            "ORIGINAL_BINOMIAL_BOUND_NOT_REPRODUCED")
    return {
        "schema":"anto_v37_original_published_natural_regain_white_petal_rescue_pollinator_evidence",
        "original_source_registry":"data/external/anto_v37_natural_regain_and_native_reactivation_original_sources_20261116.json",
        "n_original_publications":len(studies),
        "registered_original_DOIs":sorted(s["doi"] for s in studies),
        "n_independent_botanical_families":len({s["family"] for s in studies}),
        "n_distinct_petunia_radiations":len({s["shared_radiation"] for s in studies if s["family"]=="Solanaceae"}),
        "related_Petunia_studies_not_independent_macro_evolutionary_replicates":True,
        "natural_regain_source_species":sorted(mapping),
        "natural_regain_source_study_count":len(natural),
        "natural_regain_mechanism_count":len(set(mapping.values())),
        "natural_regain_route_by_species":mapping,
        "natural_AN2_repair":{"gene":"AN2","indel":an["indel"],
                              "state":"READING_FRAME_RESTORED","species":"Petunia secreta"},
        "natural_DPL_paralog_redeployment":{
            "species":"Petunia exserta","AN2_still_nonfunctional":d["AN2_still_nonfunctional"],
            "regulator":d["used_regulator"],
            "spatial_expression_change":d["regulatory_spatial_shift"],
            "other_pathway_changes":d["additional_required_pathway_changes"],
            "native_DPL_silencing_white_sectors":d["DPL_native_VIGS_white_sectors"],
            "red_colour_assumed_pelargonidin":False},
        "observed_in_real_time_ancestral_to_descendant_reversion":False,
        "n_original_native_white_petal_complementation_studies":len(rescue),
        "n_families_with_native_white_petal_complementation":len({x["family"] for x in rescue}),
        "n_stable_native_white_petal_complementation_studies":len(stable),
        "n_transient_local_native_white_petal_complementation_studies":len(transient),
        "native_white_petal_restorations":[
            {"original_study":s["doi"],"family":s["family"],
             "original_white_host":s["original_white_petal_transformation"].get(
                 "native_white_species",s["original_white_petal_transformation"].get("native_white_material")),
             "regulator":s["original_white_petal_transformation"]["regulator"],
             "experiment_type":s["original_white_petal_transformation"]["design"],
             "natural_historical_reversion":False}
            for s in rescue],
        "Petunia_axillaris_DPL_petal_rescue":{
            "source":be["doi"],"independent_lines_created":tr["generated_independent_transgenic_lines"],
            "lines_assayed_for_chemistry":tr["lines_selected_for_chemical_assay"],
            "HPLC_transgenic_plants":tr["HPLC_transgenic_samples"],
            "HPLC_control_plants":tr["HPLC_nontransgenic_samples"],
            "transgenic_plants_per_assayed_line":tr["HPLC_transgenics_per_assayed_line"],
            "corolla_pigmentation_restored":True,
            "not_longitudinal_ancestral_reversal":True},
        "Petunia_axillaris_AN2_petal_rescue":{
            "source":ho["doi"],
            "original_white_host":an2["native_white_species"],
            "gene_donor":an2["donor_of_functional_regulator"],
            "regulator":an2["regulator"],
            "spontaneous_ancestral_pseudogene_repair":False},
        "Actinidia_white_hybrid_MYB110a_petal_rescue":{
            "source":fr["doi"],
            "host_white_genotype":"21-05-09e",
            "introduced_gene":ax["regulator"],
            "local_cells_pigmented":ax["pigmented_cells_local_to_infiltration"],
            "GFP_control_negative":True,
            "stable_intergenerational_inheritance_established":False},
        "original_pollinator_preference":{
            "source":ho["doi"],"original_study_count":len(poll),
            "experimental_scope":"controlled greenhouse first-choice/visit behavior, not wild fitness",
            "Manduca_first_choose_white":pol["moth_first_choice_white_P_axillaris"],
            "Manduca_first_choose_coloured":pol["moth_first_choice_ANT2_restored_P_axillaris"],
            "Manduca_n":pol["moth_first_choice_total"],
            "Manduca_author_p_less_than":pol["moth_source_p_upper_bound"],
            "Manduca_exact_two_sided_binomial_p":binomial_p,
            "Manduca_exact_binomial_null":"equal_first_choice_probability_0.5",
            "not_phylogenetic_or_fitness_p_value":True,
            "Bombus_n_trials":pol["bee_n_trials"],
            "Bombus_direction":"AN2_COLOURED_P_AXILLARIS",
            "Bombus_author_p":pol["bee_author_reported_p"],
            "field_seed_set_or_offspring_fitness_measured":False},
        "measured_natural_outcrossing_advantage":False,
        "pollinator_choices_equivalent_to_seed_fitness":False,
        "combined_causal_genotype_to_wild_fitness_across_programs":False,
        "all_white_flowering_plants_functionally_restorable":False,
        "natural_regain_probability":None,
        "cross_source_pooled_regain_effect":None,
        "independent_original_raw_sequence_reanalysis_count":0,
        "exact_v33_F3Hprime_F3H5prime_FLS_RMSE_externally_replayed":False,
        "original_assays_and_evolutionary_history_are_distinct_grains":True,
        "not_independent_phylogenetic_radiation_level_replicates":True,
        "source_values_author_published_not_remeasured":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    result=review(json.loads(SOURCE.read_text(encoding="utf-8")))
    if args.expected and json.loads(args.expected.read_text(encoding="utf-8"))!=result:
        raise ValueError("V37_FROZEN_AUTHOR_ORIGINAL_EVIDENCE_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ANTO_V37_ORIGINAL_REGAIN_AND_FUNCTIONAL_RESCUE",
          json.dumps({"n_papers":result["n_original_publications"],
                      "n_families":result["n_independent_botanical_families"],
                      "n_Petunia_radiations":result["n_distinct_petunia_radiations"],
                      "natural_regain_mechanisms":result["natural_regain_route_by_species"],
                      "native_white_colla_restoration_programs":result["n_original_native_white_petal_complementation_studies"],
                      "controlled_moth_first_choices":[result["original_pollinator_preference"]["Manduca_first_choose_white"],
                                                       result["original_pollinator_preference"]["Manduca_first_choose_coloured"]],
                      "no_wild_seedset":not result["original_pollinator_preference"]["field_seed_set_or_offspring_fitness_measured"]},
                     sort_keys=True))
if __name__=="__main__":main()
