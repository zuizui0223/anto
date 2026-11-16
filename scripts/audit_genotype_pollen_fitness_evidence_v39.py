"""ANTO v39: original floral colour–pollen–offspring evidence chain.

Four original paper programs across three current plant families, with original
source-grain guards. No new raw pollen/fruit effects after Dryad 401, and no
numeric pooling of incomparable experiments or pseudoreplicated petals.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STUDIES=ROOT/"data/external/anto_v39_original_genotype_pollen_fruit_studies_20261116.json"
ACCESS=ROOT/"results/anto_v39_phlox2026_dryad_original_access_status_20261116.json"

DOIS={
    "antirrhinum_2001":"10.2307/2657109",
    "phlox_2012_reinforcement":"10.1126/science.1215198",
    "wahlenbergia_2013":"10.1111/1365-2435.12131",
    "phlox_2026_si":"10.1073/pnas.2525299123",
}
FAMILIES={
    "antirrhinum_2001":"Plantaginaceae",
    "phlox_2012_reinforcement":"Polemoniaceae",
    "wahlenbergia_2013":"Campanulaceae",
    "phlox_2026_si":"Polemoniaceae",
}

def require(x,code):
    if not x: raise ValueError(code)


def review(source,access):
    require(source.get("schema")=="anto_v39_original_cross_family_colour_pollination_reproduction_evidence_registry",
            "SOURCE_REGISTRY_SCHEMA_DRIFT")
    studies=source.get("studies")
    require(isinstance(studies,list) and len(studies)==4 and
            len({s["id"] for s in studies})==4,
            "ORIGINAL_STUDY_SET_DRIFT")
    rows={r["id"]:r for r in studies}
    require({k:v["doi"] for k,v in rows.items()}==DOIS,
            "ORIGINAL_SOURCE_DOI_DRIFT")
    require({k:v["family_current"] for k,v in rows.items()}==FAMILIES
            and rows["antirrhinum_2001"]["family_historical_original"]=="Scrophulariaceae",
            "ORIGINAL_CURRENT_FAMILY_DRIFT")
    require(source["original_paper_program_count"]==4
            and source["family_replicate_count"]==3
            and source["original_phlox_papers_are_one_shared_genus_not_independent_radiations"] is True
            and source["no_pooled_effect_size_assumption"] is True
            and source["same_species_controlled_color_gene_to_free_pollinator_to_seedset_chain_measured"] is False
            and source["prior_v38_original_source_untouched"] is True,
            "INDEPENDENCE_OR_MEASUREMENT_FIREWALL_MISSING")

    a=rows["antirrhinum_2001"]
    am=a["measurements"]
    require(a["genus"]=="Antirrhinum"
            and all(am[x] is True for x in (
                "heritable_flower_colour","pollinator_visit_number","visit_duration",
                "pollen_receipt","pollen_export","offspring_allele_genotype_frequencies"))
            and am["genetically_engineered_flower_colour"] is False,
            "ANTIRRHINUM_CAUSAL_EDIT_NOT_PRESENT")
    require(am["longitudinal_ancestral_reversal"] is False
            and am["flower_fruit_seed_count_directly_reported_in_abstract"] is False
            and a["source_measurement_bridge"]==
                "1997_within_array_visits_and_offspring;1994_1995_separate_pollen_receipt_export_arrays"
            and a["direct_flower_gene_edit_to_fitness"] is False
            and a["raw_individual_data_acquired_by_ANTO_v39"] is False,
            "ANTIRRHINUM_VISIT_POLLEN_OFFSPRING_GRAIN_DRIFT")
    require(a["reported_mechanism"]==
            "pollen receipt depended most on visit count; pollen export best predicted by time inside flower; assortative visits aligned with offspring allele/genotype frequency departures",
            "ANTIRRHINUM_AUTHOR_MECHANISM_DRIFT")
    expected_a={
        "1994_Barr":"visit_counts_durations_pollen_receipt_export_and_dye;RED_AND_YELLOW_COMMERCIAL",
        "1995_Avery":"visit_counts_durations_pollen_receipt_export_and_dye;WHITE_AND_YELLOW_BACKCROSS",
        "1997_Rosy":"same_array_visit_behavior_and_offspring_genotypes;WHITE_YELLOW_NO_NECTAR_AUGMENTATION",
        "1997_Barr":"same_array_visit_behavior_and_offspring_genotypes;WHITE_YELLOW_YELLOW_NECTAR_SUPPLEMENTED",
        "pollen_receipt_export_and_1997_offspring_not_exactly_same_array_and_flower":True,
    }
    require(a["original_measurement_cohort_mapping"]==expected_a
            and "not individually isolated" in a["original_interpretation_limit"],
            "ANTIRRHINUM_FIELD_YEAR_COHORT_LINK_DRIFT")
    outcomes=a["original_author_1997_offspring_frequencies"]
    exp={"Rosy":(.5,.44,.5,.41,False),"Barr":(.5,.52,.5,.43,True)}
    require(set(outcomes)==set(exp),"ANTIRRHINUM_SOURCE_ARRAY_SET_DRIFT")
    for key,(p0,p1,h0,h1,nectar) in exp.items():
        o=outcomes[key]
        require(o["initial_sulfurea_allele_frequency"]==p0
                and o["offspring_sulfurea_allele_frequency"]==p1
                and o["initial_heterozygote_fraction"]==h0
                and o["offspring_heterozygote_fraction"]==h1
                and o["nectar_supplemented_to_yellow"] is nectar
                and o["original_field_array_year"]==1997,
                "ANTIRRHINUM_AUTHOR_1997_OFFSPRING_FREQUENCY_DRIFT")

    p=rows["phlox_2012_reinforcement"]
    pm=p["measurements"]
    require(p["genus"]=="Phlox" and p["flower_colour_genetics"]==
            ["F3'5'H hue","R2R3-MYB intensity"]
            and all(pm[x] is True for x in (
                "flower_colour_genotype","pollinator_inter_species_movement",
                "hybridization","selection_field_context"))
            and pm["floral_gene_edit_no"] is False
            and pm["independent_multigeneration_colour_revival"] is False
            and p["allopatric_reinforcing_selection_detected"] is False
            and p["cannot_infer_absolute_maternal_seed_advantage_from_selection_metric"] is True
            and p["do_not_use_different_Dryad_8792d_as_fitness_data"] is True
            and p["raw_individual_data_acquired_by_ANTO_v39"] is False,
            "PHLOX_REINFORCEMENT_SOURCE_DRIFT")

    pub=p["original_author_2012_field_denominators"]
    require(pub["common_garden_seeds_planted"]==2720
            and pub["survival_outcome_figure_n"]==1909
            and pub["fruit_production_figure_n"]==1240
            and pub["hybridization_paternity_figure_n"]==931
            and pub["female_fitness_genotype_direct_significant_difference_detected"] is False
            and pub["paternity_reference_hybridization_min"]==.28
            and pub["paternity_reference_hybridization_max"]==.44
            and pub["paternity_mean_reference_hybridization_used_in_selection_model"]==.43
            and pub["interspecific_hybrids_assumed_sterility_for_model"]==.90
            and pub["dark_intensity_allele_reported_hybridization_reduction"]=="GREATER_THAN_50_PERCENT"
            and pub["author_derived_selection_coefficient_dark_allele"]==.32
            and pub["hue_red_allele_hybridization_effect_detected"] is False
            and p["conditional_selection_cannot_be_merged_with_red_hue_effect"] is True,
            "PHLOX_DARK_INTENSITY_NOT_RED_HUE_SELECTION_SOURCE_DRIFT")
    poll=p["original_author_2012_pollinator_array_denominators"]
    require(poll["observed_pollinators"]==181
            and poll["observed_total_transitions_between_plants"]==2301
            and poll["Battus_philenor_pollinators"]==108
            and poll["skippers_Hesperiidae_pollinators"]==73
            and poll["Battus_philenor_pollinators"]+poll["skippers_Hesperiidae_pollinators"]==
                poll["observed_pollinators"]
            and poll["visitation_tracking_arrangement"]==
               "THREE_OBSERVATIONAL_ARRAYS_SEPARATE_FROM_FOUR_PATERNITY_GARDEN_BLOCKS"
            and poll["pollinator_movement_events_not_one_to_one_linked_to_genotyped_seeds"] is True,
            "PHLOX_2012_POLLINATOR_EVENT_GRAIN_DRIFT")

    w=rows["wahlenbergia_2013"]
    wm=w["measurements"]
    require(w["genus"]=="Wahlenbergia"
            and all(wm[x] is True for x in (
                "visible_colour_quantification","pollen_hand_supplementation",
                "pollen_reduction","per_plant_seed_production"))
            and wm["floral_colour_gene_manipulation"] is False
            and wm["pollinator_direct_pollen_grain_receipt_assayed"] is False
            and wm["ancestral_trait_regain"] is False
            and w["original_author_reported_added_pollen_seed_increase_fraction"]==.43
            and w["selection_favors_whiter_flowers"] is True
            and w["selection_white_stronger_under_added_pollen"] is True
            and w["resource_mechanism_proven"] is False
            and w["source_mechanism"]==
               "pollen_limitation_observed_but_colour_selection_not_explained_by_pollen_receipt"
            and w["raw_individual_data_acquired_by_ANTO_v39"] is False,
            "WAHLENBERGIA_NONPOLLEN_MECHANISM_SOURCE_DRIFT")

    b=rows["phlox_2026_si"]
    bm=b["measurements"]
    require(b["genus"]=="Phlox"
            and b["publication_preprint_superseded"]=="10.1101/2025.06.17.660178"
            and b["data_doi"]=="10.5061/dryad.zgmsbccr3"
            and all(bm[x] is True for x in (
                "targeted_PdPIRK_silencing","DFR_white_colour_reporter",
                "hand_crossed_pollen_allele_genotype","stigma_pollen_tubes","fruit_binary"))
            and bm["free_flying_pollinator_choice"] is False
            and bm["observed_plant_multigeneration_colour_regain"] is False
            and b["gene_colour_reporter_is_not_same_as_pistil_recognition_gene"] is True
            and b["original_gene_causality"]==
                "PdPIRK allele-specific pollen recognition gene, not DFR flower-color gene",
            "PDPIRK_DFR_MECHANISM_CONFUSED")
    required_csv={"vigs_fruitset.csv","vigs_pollentubes.csv","vigs_qpcr.csv",
                  "diallelcross_fruitdata.csv","diallelcross_pollentubedata.csv"}
    require(set(b["source_files_to_authenticate"])==required_csv
            and b["source_raw_reanalysis_status"]=="TO_BE_COMPUTED_IF_OFFICIAL_DRYAD_BYTES_ACCESSIBLE",
            "DRYAD_ORIGINAL_SMALL_SOURCE_LIST_DRIFT")
    require(access.get("schema")=="anto_v39_original_phlox_PNAS2026_dryad_source_access"
            and access.get("dataset_DOI")==b["data_doi"]
            and access.get("latest_published_paper_DOI")==b["doi"]
            and access.get("original_Dryad_official_API_response")=="SOURCE_ACCESS_RESTRICTED_HTTP_401"
            and access.get("row_level_fruit_pollen_tubes_reanalysis_completed") is False
            and access.get("verified_original_CSV_downloads")==0
            and access.get("new_PdPIRK_gene_allele_fruit_effect_estimated") is False
            and access.get("technical_access_failure_not_biological_negative") is True
            and access.get("no_raw_data_substitution_or_permission_bypass") is True,
            "DRYAD_ACCESS_CLAIM_DRIFT")

    return {
        "schema":"anto_v39_cross_family_colour_pollinator_mating_fitness_chain_source_audit",
        "n_original_published_study_programs":len(studies),
        "n_distinct_current_plant_families":len(set(FAMILIES.values())),
        "n_original_Phlox_studies_same_genus":2,
        "original_papers_not_independent_phylogenetic_transition_events":True,
        "original_source_DOIs_sorted":sorted(DOIS.values()),
        "same_study_confirmed_processes":{
            "Antirrhinum_array_visit_to_pollen_to_offspring_allele":{
                "doi":a["doi"],
                "flower_colour_locus_inherited":True,
                "visit_time_and_pollen_receipt_export_measured":True,
                "pollen_measurements_1994_1995_separate_from_offspring_1997":True,
                "offspring_genotypes_measured":True,
                "genic_nucleotide_identity_or_transgenic_phenotype_not_tested":True,
                "long_term_allele_fixation_not_observed":True,
            },
            "Phlox_SI_receptor_allele_to_controlled_pollen_to_fruit_original_authors":{
                "doi":b["doi"],
                "allele_specific_PdPIRK_pistil_RNA_silencing":True,
                "DFR_petal_bleaching_is_reporter_not_Fruit_gene":True,
                "pollen_tubes_and_fruit_original_author_measured":True,
                "natural_free_foraging_visitation_NOT_measured":True,
                "original_source_rows_authenticated_for_new_fit":False,
                "new_same_individual_gene_colour_seed_effect_computed":False
            }
        },
        "source_stage_matrix":{
            "antirrhinum_2001":{
                "natural_pollinator_visits":True,
                "pollen_receipt_and_export":True,
                "offspring_allele_genotype_measured":True,
                "original_within_array_visits_pollen_offspring_chain":False,
                "original_within_1997_arrays_visitation_and_offspring_genotypes":True,
                "pollen_receipt_export_from_1994_1995_other_arrays":True,
                "original_flower_colour_gene_edited":False,
                "natural_multigeneration_fixation_measured":False,
            },
            "phlox_2012_reinforcement":{
                "field_selection_and_hybridization_reported":True,
                "pollinator_interspecific_movement":True,
                "author_reported_dark_allele_selection_coefficient":pub["author_derived_selection_coefficient_dark_allele"],
                "inferred_coefficient_assumes_hybrid_sterility":pub["interspecific_hybrids_assumed_sterility_for_model"],
                "dark_allele_hybridization_reduction_more_than_half":True,
                "red_hue_allele_hybridization_reduction_detected":False,
                "original_paternity_genotyped_seeds_n":pub["hybridization_paternity_figure_n"],
                "source_pollinator_individuals":poll["observed_pollinators"],
                "source_interplant_transitions":poll["observed_total_transitions_between_plants"],
                "paternity_and_visitor_arrays_same_plants":False,
                "allopatric_reinforcing_selection_detected":False,
                "natural_ancestral_colour_regain_observed":False,
            },
            "wahlenbergia_2013":{
                "pollen_limited_seed_production":True,
                "published_mean_pollen_supplementation_seed_increase_percent":43,
                "selection_on_whiteness_stronger_with_excess_pollen":True,
                "floral_colour_selection_proved_pollen_receipt_mediated":False,
                "genotyped_flower_colour_locus":False,
            },
            "phlox_2026_si":{
                "allele_specific_receptor_intervention_reported":True,
                "hand_cross_pollen_genotype_and_pollen_tubes":True,
                "controlled_cross_fruit_binary_outcome":True,
                "DFR_petal_bleaching_is_visible_silencing_reporter":True,
                "DFR_flower_colour_causally_changes_compatibility_demonstrated":False,
                "free_foraging_pollinator_choice_study":False,
            }
        },
        "Antirrhinum_1997_original_author_array_allele_shifts":{
            key:{
                "field_year":1997,
                "nectar_supplemented_on_yellow":item["nectar_supplemented_to_yellow"],
                "sulfurea_allele_initial":item["initial_sulfurea_allele_frequency"],
                "sulfurea_allele_offspring":item["offspring_sulfurea_allele_frequency"],
                "sulfurea_delta_offspring_minus_initial":round(
                    item["offspring_sulfurea_allele_frequency"]-
                    item["initial_sulfurea_allele_frequency"],8),
                "heterozygote_initial":item["initial_heterozygote_fraction"],
                "heterozygote_offspring":item["offspring_heterozygote_fraction"],
                "heterozygote_delta_offspring_minus_initial":round(
                    item["offspring_heterozygote_fraction"]-
                    item["initial_heterozygote_fraction"],8),
            } for key,item in sorted(outcomes.items())
        },
        "Antirrhinum_pollen_arrays_1994_1995_not_same_year_as_1997_offspring":True,
        "Antirrhinum_1997_neither_colour_only_causal_effect_nor_multigeneration_selection_estimated":True,
        "Phlox_2012_author_derived_dark_intensity_selection":{
            "source_DOI":p["doi"],
            "dark_intensity_allele_s_author_inferred":.32,
            "hybrid_sterility_assumption_for_s":.90,
            "paternity_reference_hybridization_range":[.28,.44],
            "mean_reference_hybridization_in_author_model":.43,
            "original_paternity_figure_seeds_genotyped_n":931,
            "dark_intensity_hybridization_lower_by_more_than_half":True,
            "hue_red_allele_hybridization_effect_detected":False,
            "direct_female_fecundity_genotype_significant_difference_detected":False,
            "raw_individual_seed_genotype_data_reanalyzed":False,
            "original_2012_field_common_garden_seeds_planted":2720,
            "different_field_arrays_used_for_pollen_movement_and_seed_paternity":True
        },
        "Phlox_2012_pollinator_source_count":{
            "total_pollinator_individuals":181,
            "Battus_philenor":108,
            "skippers":73,
            "interplant_transitions":2301,
            "not_same_unit_as_paternity_genotyped_seeds":True
        },
        "sympatry_specific_fitness_selection_not_universal_colour_benefit":True,
        "source_wahlenbergia_original_response":"SEED_OUTPUT",
        "author_published_43_percent":43,
        "different_Phlox_microsatellite_Dryad_not_fitness_reanalysis":True,
        "original_Phlox2026_Dryad_source_status":"SOURCE_ACCESS_RESTRICTED_HTTP_401",
        "raw_Phlox2026_fruit_pollen_tube_rows_reanalyzed":0,
        "new_original_PdPIRK_allele_specific_fruit_risk_ratio":None,
        "technical_source_failure_not_biological_zero_effect":True,
        "2025_bioRxiv_is_duplicate_not_new_paper":True,
        "published_2026_PNAS_DOI":"10.1073/pnas.2525299123",
        "original_Antirrhinum_colour_gene_causal_nucleotide_identified_v39":False,
        "one_locus_mendelian_colour_state_in_original_Antirrhinum":True,
        "cross_study_genetic_pollen_fitness_pooled_beta":None,
        "identical_gene_colour_treatment_replicated_across_all_studies":False,
        "n_same_AN2_colour_gene_to_wild_seedset_new_replicates":0,
        "all_five_stages_same_gene_plant_cohort_demonstrated":False,
        "subprocesses_identified_but_not_one_connected_experiment":True,
        "prior_Petunia_v38_separate_and_not_overridden":True,
        "confirmed_natural_ancestral_colour_regain_rate":False,
        "author_published_studies_not_original_assay_bytes_remeasured":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    result=review(json.loads(STUDIES.read_text()),json.loads(ACCESS.read_text()))
    if args.expected and json.loads(args.expected.read_text())!=result:
        raise ValueError("V39_AUTHOR_ORIGINAL_STAGE_AUDIT_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ANTO_V39_REPRODUCTIVE_GATE_EVIDENCE",json.dumps({
        "programs":result["n_original_published_study_programs"],
        "families":result["n_distinct_current_plant_families"],
        "Antirrhinum_visit_to_offspring":True,
        "Phlox_PdPIRK_original_gene_to_fruit":True,
        "Wahlenbergia_seed_increase_on_hand_pollen_pct":43,
        "Wahlenbergia_colour_selection_pollen_receipt_mediated":False,
        "original_Dryad_reanalysis_admitted":False,
        "new_completed_full_AN2_colour_gene_wild_fitness_chain":False,
    },sort_keys=True))


if __name__=="__main__":
    main()
