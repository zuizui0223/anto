"""ANTO v38: independently published pollinator-channel and service-to-seed gates.

All quantitative rows are author-reported ORIGINAL study denominators; only
exact binomial first-choice arithmetic is newly recomputed. Not a meta-analysis,
raw Petunia gene/seed assay refit, AN2 colour→fruit causal experiment or
multi-generation natural evolutionary trajectory.
"""
from __future__ import annotations

import argparse
import json
from math import comb, fsum, isclose
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/external/anto_v38_petunia_original_pollinator_seed_set_programs_20261116.json"
PARENT=ROOT/"results/anto_v37_original_natural_regain_native_rescue_20261116.json"
ORIGINAL={
    "hoballah2007_AN2":"10.1105/tpc.106.048694",
    "dellolivo2013_AN2":"10.1111/evo.12178",
    "brandenburg2012_nectar":"10.1016/j.cub.2012.06.058",
    "rodrigues2018_secreta":"10.1093/aobpla/ply057",
}


def require(ok,message):
    if not ok: raise ValueError(message)


def exact_two_sided_binomial(first,second):
    """Two-sided simple p=.5 binomial tail, not authors' chi-square tests."""
    if any(isinstance(k,bool) or not isinstance(k,int) or k<0 for k in (first,second)):
        raise ValueError("FIRST_CHOICE_COUNTS_NOT_ADMITTED")
    n=first+second
    if n<2: raise ValueError("INSUFFICIENT_FIRST_CHOICES")
    tail=fsum(comb(n,k) for k in range(max(first,second),n+1))/(2**n)
    return min(1.0,2*tail)


def review(source,parent):
    require(source.get("schema")=="anto_v38_original_petunia_pollinator_channel_effectiveness_and_seed_source_registry",
            "V38_SOURCE_SCHEMA_DRIFT")
    programs=source.get("programs")
    require(isinstance(programs,list) and len(programs)==4 and
            len({p.get("id") for p in programs})==4,
            "PAPER_STUDY_SET_DRIFT")
    rows={p["id"]:p for p in programs}
    require({k:v["doi"] for k,v in rows.items()}==ORIGINAL,
            "ORIGINAL_STUDY_DOI_DRIFT")
    require(all(p["family"]=="Solanaceae" and p["shared_radiation"]=="Petunia"
                and p["n_original_studies"]==1 for p in programs)
            and source["all_four_papers_one_shared_Petunia_radiation"] is True
            and source["duplicate_study_policy"]==
            "HOBALLAH_2007_ALREADY_IN_V37_DO_NOT_COUNT_AS_INDEPENDENT_NEW_STUDY",
            "FAMILY_RADIATION_REPLICATION_OR_DUPLICATE_DRIFT")
    require(all(source[k] is True for k in (
        "no_original_raw_pollen_or_seed_numbers_reconstructed",
        "no_equivalence_claim_from_non_significant_P_values",
        "no_direct_AN2_colour_change_to_field_seedset_link",
        "no_source_has_natural_genetic_trait_frequency_time_series",
        "no_multi_study_numeric_effect_pooling",
    )),"SOURCE_INTERPRETATION_GUARD_DRIFT")
    require(parent.get("schema")==
            "anto_v37_original_published_natural_regain_white_petal_rescue_pollinator_evidence",
            "V37_PARENT_SCHEMA_DRIFT")
    old=parent["original_pollinator_preference"]
    require(old["Manduca_first_choose_white"]==30 and
            old["Manduca_first_choose_coloured"]==6 and
            old["Manduca_n"]==36,
            "V37_SOURCE_COUNT_DRIFT")

    h=rows["hoballah2007_AN2"]["experimental_substudies"]
    t=h["transgenic_colour_choice"]
    require(t["species"]=="Petunia axillaris" and
            t["intervention"]=="STABLE_FUNCTIONAL_AN2_ADD_TO_WHITE" and
            t["readout"]=="GREENHOUSE_FIRST_CHOICE" and
            t["moth"]=="Manduca sexta" and
            t["moth_first_choice_wild_white"]==30 and
            t["moth_first_choice_transgenic_coloured"]==6 and
            t["n_first_choices"]==36 and
            t["moth_first_choice_wild_white"]+
            t["moth_first_choice_transgenic_coloured"]==t["n_first_choices"]
            and t["author_p_upper"]==.0001 and
            t["study_specific_colour_contrast"] is True,
            "HOBALLAH_COLOUR_TRANSFORMED_CHOICE_DENOMINATOR_DRIFT")
    require(t["same_treatment_stigma_pollen_deposition_measured"] is False
            and t["same_treatment_capsule_or_seedset_measured"] is False
            and t["same_treatment_heritable_natural_population_frequency_measured"] is False,
            "AN2_FIRST_CHOICE_NOT_FECUNDITY")
    p2007=exact_two_sided_binomial(30,6)
    require(isclose(p2007,old["Manduca_exact_two_sided_binomial_p"],abs_tol=1e-16,rel_tol=0)
            and p2007<t["author_p_upper"],
            "V37_EXACT_BINOMIAL_ARITHMETIC_DRIFT")

    wild=h["independent_wild_native_pollinator_exclusion"]
    require(wild["organism_genotype_equivalent_to_transgenic_color_choice"] is False and
            wild["genetic_intervention"]=="NONE_IN_NATURAL_POPULATIONS"
            and wild["population_taxa"]==["Petunia axillaris","Petunia integrifolia"]
            and wild["measurement"]=="CAPSULE_FORMATION_NOT_SEED_NUMBER"
            and wild["separate_wild_genotypes_not_edited_AN2_colour_contrast"] is True,
            "FORGED_COLOUR_TO_WILD_CAPSULE_LINK")
    site1=wild["p_axillaris_Jose_Ignacio"]
    site2=wild["p_axillaris_Carmelo"]
    require(site1["day_exclusion_and_night_exclusion_both_made_capsules"] is True
            and site2["day_exclusion_and_night_exclusion_both_made_capsules"] is True
            and site1["original_Mann_Whitney_Z"]==-.711
            and site2["original_Mann_Whitney_Z"]==-1.38
            and site1["original_Mann_Whitney_p"]==.477
            and site2["original_Mann_Whitney_p"]==.167
            and site1["reported_original_replicates"]==4
            and site2["reported_original_replicates"]==4,
            "WILD_CAPSULE_EXCLUSION_AUTHOR_RESULTS_DRIFT")
    wild_int=wild["p_integrifolia_night_only"]
    require(wild_int["capsule_formation_in_night_only"] is False
            and wild_int["reported_original_replicates"]==4
            and wild_int["original_sites"]==["Puerto Viejo","Rivera"],
            "SPECIES_SPECIFIC_DAY_VS_NIGHT_CAPSULE_DRIFT")

    d=rows["dellolivo2013_AN2"]["experimental_substudies"]
    ga=d["axillaris_colour_gain"]
    lo=d["inflata_colour_loss"]
    require(ga["species"]=="Petunia axillaris"
            and ga["intervention"]=="NIL_FUNCTIONAL_AN2_GAIN_WHITE_TO_COLOURED"
            and ga["moth_first_choice_An2plus_coloured"]==15
            and ga["moth_first_choice_original_white"]==14
            and ga["moth_n_first_choices"]==29
            and ga["author_chi_square_p"]==.85268
            and ga["moth_established_guild_difference_detected"] is False
            and ga["moth_equivalence_proved"] is False,
            "AN2_GAIN_MOTH_FIRST_CHOICES_AUTHOR_DRIFT")
    require(ga["greenhouse_bumblebee_genotype_preference"]=="AN2PLUS_COLOURED"
            and ga["greenhouse_bumblebee_original_hives"]==4
            and ga["greenhouse_bumblebee_author_p"]==.001
            and ga["native_solitary_bee_genotype_preference"]=="AN2PLUS_COLOURED"
            and ga["field_bee_original_observation_days"]==7
            and ga["field_bee_author_p"]==.001,
            "AN2_GAIN_BEE_VISIT_CHANNEL_DRIFT")
    require(lo["species"]=="Petunia inflata"
            and lo["intervention"]=="NIL_NONFUNCTIONAL_AN2_LOSS_PURPLE_TO_WHITE"
            and lo["moth_first_choice_An2minus_white"]==15
            and lo["moth_first_choice_original_purple"]==3
            and lo["moth_n_first_choices"]==18
            and lo["author_chi_square_p"]==.005
            and lo["moth_novel_guild_increased_visitation"] is True,
            "AN2_LOSS_MOTH_FIRST_CHOICES_AUTHOR_DRIFT")
    require(lo["greenhouse_bumblebee_genotype_difference_detected"] is False
            and lo["greenhouse_bumblebee_original_hives"]==4
            and lo["greenhouse_bumblebee_author_p"]==.525
            and lo["native_solitary_bee_genotype_difference_detected"] is False
            and lo["field_bee_original_observation_days"]==7
            and lo["field_bee_author_p"]==.277
            and lo["additional_flower_projected_area_change"]==
            "AN2MINUS_INFLATA_PROJECTED_AREA_28_PERCENT_SMALLER_THAN_PARENT",
            "AN2_LOSS_EXISTING_BEE_GUILD_CAVEATS_DRIFT")
    require(all(item["stigma_pollen_deposition_measured"] is False and
                item["same_NIL_seedset_measured"] is False
                for item in (ga,lo))
            and lo["wild_natural_evolution_observed_over_time"] is False
            and d["no_direct_evidence_of_new_guild_reproductive_assurance_in_same_NIL"] is True,
            "NIL_POLLEN_FITNESS_OUTCOME_NOT_MEASURED")

    n=rows["brandenburg2012_nectar"]["experimental_substudies"]["low_nectar_introgression"]
    require(n["species"]=="Petunia axillaris"
            and n["intervention"]=="LOW_NECTAR_INTROGRESSION_NOT_AN2_COLOUR"
            and n["low_nectar_relative_wildtype"]=="APPROXIMATELY_ONE_THIRD_VOLUME"
            and n["initial_attraction_discrimination_detected"] is False
            and n["probing_duration_low_nectar"]=="SHORTER"
            and n["seed_set_low_nectar_under_moth_service"]=="LOWER_THAN_NORMAL_NECTAR"
            and n["low_nectar_manual_hand_pollination_reproductive_potential"]==
            "HIGHER_THAN_MOTH_SERVICE"
            and n["supplementation_to_parental_nectar_level"]==
            "RESTORED_SEED_SET_NEAR_HAND_POLLINATION"
            and n["main_original_conclusion"]==
            "PROBING_DURATION_MEDIATES_REPRODUCTIVE_OUTPUT",
            "NECTAR_SEEDSET_SERVICE_SOURCE_DRIFT")
    require(n["original_paper_raw_seed_set_counts_recomputed_here"] is False
            and n["novel_flower_colour_or_AN2_perturbation"] is False,
            "NECTAR_NOT_AN2_COLOUR")

    s=rows["rodrigues2018_secreta"]["experimental_substudies"]["natural_pink_Psecreta"]
    require(s["species"]=="Petunia secreta"
            and s["field_study_years"]==[2014,2015]
            and s["original_bee_genus"]=="Pseudagapostemon"
            and s["total_observed_bee_visits"]==38
            and s["bee_visits_with_documented_reproductive_structure_contact"]==22
            and s["pollen_collection"] is True
            and s["other_visitors_birds_and_bees_occur"] is True
            and s["source_replication_unit"]=="VISITS_NOT_DISTINCT_FLOWER_GENOTYPES",
            "PSECRETA_WILD_OBSERVATIONS_DRIFT")
    require(s["viable_pollen_deposition_determined_per_visit"] is False
            and s["measured_capsule_set_from_same_visited_flowers"] is False
            and s["measured_seed_set_from_same_visited_flowers"] is False
            and s["AN2_causal_gene_editing_in_the_wild_field_study"] is False
            and s["naturally_regained_pigment_phylogenetic_inference_from_different_Esfeld2018_study"] is True,
            "PSECRETA_CONTACT_NOT_SEED_OR_POLLEN_FITNESS")

    stage={
      "Hoballah2007_AN2_transformed_colour":{
         "gene_colour_manipulated":True,"AN2_colour_changed":True,
         "pollinator_first_choice_measured":True,
         "same_genotype_pollen_deposition_count_measured":False,
         "same_genotype_field_capsule_or_seed_number_measured":False,
         "same_genotype_wild_population_genetic_persistence_measured":False},
      "Hoballah2007_wild_pollinator_exclusion":{
         "gene_colour_manipulated":False,"pollinator_first_choice_measured":False,
         "wild_day_and_night_both_yielded_capsules":True,
         "outcome":"WILD_CAPSULE_FORMATION_NOT_AN2_TRANSGENE"},
      "DellOlivo2013_axillaris_AN2_gain":{
         "gene_colour_manipulated":True,"AN2_colour_changed":True,
         "pollinator_first_choice_measured":True,
         "same_genotype_pollen_deposition_count_measured":False,
         "same_genotype_field_capsule_or_seed_number_measured":False},
      "DellOlivo2013_inflata_AN2_loss":{
         "gene_colour_manipulated":True,"AN2_colour_changed":True,
         "pollinator_first_choice_measured":True,
         "same_genotype_pollen_deposition_count_measured":False,
         "same_genotype_field_capsule_or_seed_number_measured":False},
      "Brandenburg2012_low_nectar_not_colour":{
         "gene_colour_manipulated":False,"low_nectar_introgression":True,
         "pollinator_probing_time_measured":True,"same_trait_seed_set_measured":True,
         "nectar_supplementation_test":True,"not_AN2_colour_seed_set":True},
      "Rodrigues2018_wild_Psecreta_pollen_contacts":{
         "gene_colour_manipulated":False,
         "wild_reproductive_structure_contact_observed":True,
         "same_visit_pollen_deposition_quantified":False,
         "same_flowers_seed_set_measured":False}
    }
    return {
       "schema":"anto_v38_original_Petunia_channel_addition_effectiveness_and_seedset_gap",
       "original_primary_publication_count":len(programs),
       "n_original_publications_in_v38_registry":len(programs),
       "n_botanical_families":1,
       "n_independent_petunia_radiations":1,
       "n_new_genetically_independent_colour_evolution_experiments":0,
       "registered_original_DOIs_sorted":sorted(row["doi"] for row in programs),
       "hoballah2007_existing_v37_study_deduplicated":True,
       "parent_v37_first_choice_reproduced_unchanged":True,
       "author_published_count_rechecks":{
           "hoballah2007_moth_white":30,"hoballah2007_moth_coloured":6,
           "hoballah2007_moth_n":36,
           "hoballah2007_exact_binomial_two_sided_p":p2007,
           "2013_axillaris_coloured_first":ga["moth_first_choice_An2plus_coloured"],
           "2013_axillaris_white_first":ga["moth_first_choice_original_white"],
           "2013_axillaris_total_first":ga["moth_n_first_choices"],
           "2013_axillaris_author_chi_square_p":ga["author_chi_square_p"],
           "2013_axillaris_exact_binomial_two_sided_p":exact_two_sided_binomial(15,14),
           "2013_inflata_white_first":lo["moth_first_choice_An2minus_white"],
           "2013_inflata_purple_first":lo["moth_first_choice_original_purple"],
           "2013_inflata_total_first":lo["moth_n_first_choices"],
           "2013_inflata_author_chi_square_p":lo["author_chi_square_p"],
           "2013_inflata_exact_binomial_two_sided_p":exact_two_sided_binomial(15,3),
       },
       "source_chi_square_p_does_not_need_equal_exact_binomial_p":True,
       "AN2_NIL_secondary_guild_sources":{
           "axillaris_bumblebee_hives":ga["greenhouse_bumblebee_original_hives"],
           "axillaris_field_bee_days":ga["field_bee_original_observation_days"],
           "inflata_bumblebee_hives":lo["greenhouse_bumblebee_original_hives"],
           "inflata_field_bee_days":lo["field_bee_original_observation_days"],
           "new_guild_can_be_attracted_without_detectable_old_guild_loss":True,
           "non_significance_proves_equal_original_guild_visits":False},
       "wild_field_capsule_published_source":{
           "source_DOI":rows["hoballah2007_AN2"]["doi"],
           "n_source_localities_P_axillaris":2,
           "capsules_day_visitor_access_observed":True,
           "capsules_night_visitor_access_observed":True,
           "jose_ignacio_original_Mann_Whitney_p":site1["original_Mann_Whitney_p"],
           "carmelo_original_Mann_Whitney_p":site2["original_Mann_Whitney_p"],
           "n_source_reported_replicates_per_locality":4,
           "transgenic_AN2_colour_was_the_field_exclusion_treatment":False,
           "day_night_equal_effectiveness_established":False,
           "measured_seed_counts":False,
           "integrifolia_capsules_without_day_visits":False,
           "species_specific_pollinator_channel_availability":True},
       "nectar_service_efficiency_source":{
           "source_DOI":rows["brandenburg2012_nectar"]["doi"],
           "relative_nectar":n["low_nectar_relative_wildtype"],
           "initial_attraction_difference_not_detected":True,
           "probing_duration_reduced":True,
           "low_nectar_seed_production_reduced":True,
           "nectar_supplementation_recovered_seed_set":True,
           "this_is_AN2_colour_perturbation":False,
           "seed_set_numeric_risk_ratio":None,
           "seed_set_model_p_value":None,
           "original_seed_count_rows_recomputed":0},
       "wild_Psecreta_pollen_reproductive_contact":{
           "source_DOI":rows["rodrigues2018_secreta"]["doi"],
           "total_observed_pollinator_visits":s["total_observed_bee_visits"],
           "observed_reproductive_contact_visits":s["bee_visits_with_documented_reproductive_structure_contact"],
           "descriptive_fraction_of_observed_visits_with_contact":22/38,
           "successful_pollen_deposition_per_visit_confirmed":False,
           "corresponding_visited_flowers_seedset_measured":False,
           "ancestral_AN2_reactivation_experiment_on_field_plants":False},
       "experimental_contrast_observation_matrix":stage,
       "n_original_programs_with_real_capsule_or_seed_outcome":2,
       "n_original_programs_with_colour_AN2_first_choice_manipulation":2,
       "n_same_colour_AN2_intervention_with_direct_reproductive_output":0,
       "joint_intergenerational_genetic_colour_to_fitness_chains":0,
       "fully_linked_AN2_molecular_to_long_term_fitness_verdict":
           "NOT_IDENTIFIED_SAME_TREATMENT_END_TO_END",
       "same_AN2_colour_perturbation_field_capsule_or_seedset_measured":False,
       "cannot_infer_causal_same_gene_reproductive_assurance":True,
       "non_significant_old_guild_preference_is_equivalence_evidence":False,
       "natural_ancestral_allele_recovery_caused_contact_verified_in_same_individuals":False,
       "new_externally_replicated_v33_species_RNA_HPLC_model":False,
       "independent_genetic_fitness_effects_for_meta_analysis":0,
       "no_meta_analytic_effect_size_or_pseudoreplicate_reversal_rate":True,
       "source_reported_events_not_raw_assays_or_new_wetlab":True
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--expected",type=Path)
    opts=ap.parse_args()
    result=review(json.loads(SOURCE.read_text()),json.loads(PARENT.read_text()))
    if opts.expected and json.loads(opts.expected.read_text())!=result:
        raise ValueError("V38_FROZEN_AUTHOR_REPORTED_STAGE_EVIDENCE_DRIFT")
    opts.output.parent.mkdir(parents=True,exist_ok=True)
    opts.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("ANTO_V38_SOURCE_LINKAGE_GATE",json.dumps({
        "n_papers":result["n_original_publications_in_v38_registry"],
        "new_guild_an2_gain_and_loss":True,
        "moth_2013_counts":[[15,14],[15,3]],
        "moth_2013_exact_binomial":[
            result["author_published_count_rechecks"]["2013_axillaris_exact_binomial_two_sided_p"],
            result["author_published_count_rechecks"]["2013_inflata_exact_binomial_two_sided_p"]],
        "wild_axillaris_both_day_night_capsule":True,
        "nectar_behavior_drives_seed":True,
        "same_AN2_colour_to_seed_measured":result["same_AN2_colour_perturbation_field_capsule_or_seedset_measured"],
        "full_causal_chain":result["fully_linked_AN2_molecular_to_long_term_fitness_verdict"]
    },sort_keys=True))

if __name__=="__main__":main()
