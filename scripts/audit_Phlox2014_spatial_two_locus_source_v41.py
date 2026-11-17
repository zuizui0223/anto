"""ANTO v41: published original two-locus Phlox spatial fitness audit.

2014 author model/CI reproduced from Table 1, original supplementary Table S1
counts provisionally transcribed from a text-extracted PDF (NOT image/byte
validated): 11,379 vs original main article's 12,024 adults UNRECONCILED.
This module does NOT refit the paper's original MCMC, deduce observed alleles
from dominant flower phenotypes, or claim multi-generation field tracking.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"data/external/anto_v41_Phlox2014_original_spatial_study_tableS1_provisional_20261117.json"
DOI="10.1016/j.cub.2014.07.027"
COLOURS=("LB","DB","LR","DR")
EXPECTED_MAIN_ADULTS=12024
EXPECTED_SOURCE_TRANSCRIPTION_TOTAL=11379
EXPECTED_TRANSCRIPTS_BY_COLOUR={"LB":4292,"DB":1396,"LR":1084,"DR":4607}
EXPECTED_FITNESS={
    "allopatry":{"LB":1,"DB":.86,"LR":.89,"DR":.90},
    "sympatry":{"LB":.17,"DB":.53,"LR":.86,"DR":1},
}
EXPECTED_LOCUS={
    "allopatry":{
       "red_hue_allele_h":{"estimate":-.11,"author_interval_95":[-.14,-.08]},
       "dark_intensity_allele_I":{"estimate":-.14,"author_interval_95":[-.17,-.10]},
    },
    "sympatry":{
       "ancestral_blue_allele_H":{"estimate":-.14,"author_interval_95":[-.19,-.09]},
       "ancestral_light_intensity_allele_i":{"estimate":-.47,"author_interval_95":[-.69,-.30]},
    },
}
EXPECTED_EPI={
    "allopatry":{"estimate":.15,"author_interval_95":[.11,.19]},
    "sympatry":{"estimate":-.22,"author_interval_95":[-.37,0]},
}
# The 32 original Table S1 source population IDs, fixed before script runs.
EXPECTED_TRANSECTS={
    1:[301,302,303,304,305,306,307,308,314],
    2:[325,326,327,333,334,328,329],
    3:[330,311,312,310,309,332],
    4:[335,317,316,315,319,318],
    5:[324,323,322,320],
}


def require(ok,why):
    if not ok: raise ValueError(why)


def is_finite_nonnegative_int(n):
    return isinstance(n,int) and not isinstance(n,bool) and n>=0


def source_s1_summary(data):
    trans=data["source_transcription"]
    require(trans["method"]==
            "MANUAL_FROM_PUBLISHED_SUPPLEMENT_TABLE_S1_PDF_TEXT_EXTRACTION_NOT_IMAGE_VERIFIED"
            and trans["source_original_page_image_independently_verified"] is False
            and trans["original_file_byte_checksum_authenticated"] is False
            and trans["source_count_denominator_reconciliation_status"]==
                "TABLE_S1_TEXT_DIGITIZATION_UNRECONCILED"
            and trans["preliminary_text_extraction_count_mismatch_known_before_script"] is True
            and trans["transcribed_rows_are_not_complete_verified_original_mcmc_input"] is True
            and trans["main_paper_sample_total_do_not_replace_with_transcription"] is True,
            "SOURCE_S1_TEXT_TRANSCRIPTION_SCOPE_DRIFT")
    require(trans["columns"]==
            ["transect","source_population_id","LB","DB","LR","DR"],
            "DIGITIZED_S1_COLUMN_GRAIN_DRIFT")
    rows=trans.get("population_phenotype_counts")
    require(isinstance(rows,list) and len(rows)==32,
            "ORIGINAL_32_SOURCE_POPULATIONS_DRIFT")
    by=Counter()
    per_t=Counter()
    color=Counter()
    sites_4=0
    mixed=0
    for raw in rows:
        require(isinstance(raw,list) and len(raw)==6 and
                all(is_finite_nonnegative_int(n) for n in raw),
                "INVALID_DIGITIZED_SOURCE_COUNT")
        transect,pop,*values=raw
        require(transect in EXPECTED_TRANSECTS and pop in EXPECTED_TRANSECTS[transect],
                "ORIGINAL_32_SOURCE_POPULATIONS_DRIFT")
        require(pop not in by, "ORIGINAL_32_SOURCE_POPULATIONS_DRIFT")
        require(sum(values)>0, "EMPTY_SOURCE_POPULATION_NOT_ADMITTED")
        by[pop]=transect
        per_t[transect]+=sum(values)
        color.update(dict(zip(COLOURS,values)))
        sites_4+=int(all(x>0 for x in values))
        mixed+=int(values[0]>0 and values[3]>0)
    require(list(by)==[p for ts in EXPECTED_TRANSECTS.values() for p in ts],
            "ORIGINAL_32_SOURCE_POPULATIONS_DRIFT")
    require(sum(color.values())==EXPECTED_SOURCE_TRANSCRIPTION_TOTAL
            and dict(color)==EXPECTED_TRANSCRIPTS_BY_COLOUR,
            "SOURCE_S1_DIGITIZATION_SUM_DRIFT")
    require(sites_4==17 and mixed==18,
            "SOURCE_S1_PHENOTYPE_CROSSCHECK_DRIFT")
    require([sum(1 for v in by.values() if v==t) for t in range(1,6)]==
            [9,7,6,6,4],
            "ORIGINAL_TRANSECT_SAMPLE_COUNT_DRIFT")
    return {
        "source_status":"TEXT_EXTRACTED_PROVISIONAL_TABLE_S1_NOT_IMAGE_VERIFIED",
        "n_populations_digitized":len(rows),
        "transect_populations":{str(t):sum(1 for v in by.values() if v==t) for t in range(1,6)},
        "transect_colour_sample_totals":{str(t):per_t[t] for t in range(1,6)},
        "by_colour":dict(color),
        "total_digitized_flower_count_entries":sum(color.values()),
        "recombinant_morph_count":color["DB"]+color["LR"],
        "recombinant_fraction":(color["DB"]+color["LR"])/sum(color.values()),
        "sites_with_all_four_colours":sites_4,
        "sites_with_both_LB_and_DR":mixed,
        "first_population_code":rows[0][1],
        "last_population_code":rows[-1][1],
        "not_equivalent_to_genotyped_haplotypes":True,
    }


def reported_epistasis(fitness):
    # Main source paper's definitions, based on RELATIVE phenotype fitness.
    a=fitness["allopatry"]
    s=fitness["sympatry"]
    return {
        "allopatry":1-a["DB"]+a["DR"]-a["LR"],
        "sympatry":1-s["DB"]+s["LB"]-s["LR"],
    }


def review(source):
    require(source.get("schema")==
            "anto_v41_original2014_Phlox_spatial_population_selection_source_registry",
            "ORIGINAL_SOURCE_SCHEMA_DRIFT")
    s=source.get("source_original",{})
    require(s.get("doi")==DOI and s.get("main_article_claim_adult_survey_total")==
            EXPECTED_MAIN_ADULTS and s.get("main_article_populations")==32 and
            s.get("main_article_transects")==5 and
            s.get("main_article_offspring_colored_count")==2581 and
            s.get("main_article_maternal_field_plants")==382 and
            s.get("natural_multigeneration_allele_frequency_time_series_measured") is False,
            "ORIGINAL_SOURCE_DOI_DRIFT" if s.get("doi")!=DOI else
            "MULTIGENERATION_TRACKING_NOT_IN_SOURCE" if s.get("natural_multigeneration_allele_frequency_time_series_measured") is not False else
            "AUTHOR_MAIN_SAMPLE_DENOMINATOR_DRIFT")
    require(s["original_published_light_cline_width_km"]==4.
            and s["original_published_red_cline_width_km"]==3.7
            and s["source_supplement_pollen_kernel_variance_km2"]==.83
            and s["exact_original_MCMC_model_reimplemented_in_v41"] is False
            and s["replay_of_full_offspring_Table_S2_completed"] is False,
            "ORIGINAL_SPATIAL_MODEL_EVIDENCE_GRANULARITY_DRIFT")
    m=source.get("original_study_published_model",{})
    f=m.get("original_relative_fitness_per_environment_normalized_separately")
    require(f==EXPECTED_FITNESS, "AUTHOR_2014_FITNESS_TABLE_DRIFT")
    require(m.get("original_table1_locus_selection")==EXPECTED_LOCUS,
            "ORIGINAL_2014_LOCUS_TABLE1_CI_DRIFT")
    require(m.get("epistasis_original_table1")==EXPECTED_EPI,
            "AUTHOR_2014_EPISTASIS_DRIFT")
    require(all(m.get(k) is True for k in (
         "original_authors_spatial_model_selective_cause_not_uniquely_identified",
         "original_authors_different_selection_and_migration_models_validated",
         "original_2012_field_red_hue_allele_hybridization_reduction_not_detected",
         "original_2014_model_ancestral_blue_hue_selection_significant",
         "two_source_hue_effects_distinct_study_and_estimand")),
        "SOURCE_INFERENCE_FIREWALL_DRIFT")
    require(abs(m["not_2012_s_hybrid_proxy"]-.32300163132137016)<1e-13,
            "V40_DIFFERENT_SOURCE_ESTIMAND_CONFUSED")
    require(all(source.get(k) is True for k in (
         "no_molecular_resurrection_probability_in_this_source",
         "one_original_P_drumondii_radiation_not_new_35_evolutionary_events",
         "no_empirical_coloured_AN2_restoration_rate_measured")),
        "FICTITIOUS_INDEPENDENT_EVOLUTION_EVENTS")
    typed=source_s1_summary(source)
    eps=reported_epistasis(f)
    for env in ("allopatry","sympatry"):
        require(math.isclose(eps[env],EXPECTED_EPI[env]["estimate"],abs_tol=1e-13),
                "AUTHOR_2014_EPISTASIS_ARITHMETIC_MISMATCH")
    n=typed["total_digitized_flower_count_entries"]
    difference=EXPECTED_MAIN_ADULTS-n
    # The mismatch is a NON-RESOLUTION, not a newly created negative selection result.
    require(difference==645 and n!=EXPECTED_MAIN_ADULTS,
            "SOURCE_S1_DIGITIZATION_SUM_DRIFT")
    return {
       "schema":"anto_v41_source_reported_spatial_two_locus_fitness_and_population_digitization_audit",
       "original_publication_doi":DOI,
       "n_original_studies":1,
       "source_primary_2014_original_main_and_supplement_urls":[s["original_author_pdf"],
                s["original_supplement_text_mirror"]],
       "original_main_article_observation_denominators":{
           "adult_phenotypes_main_paper_reported":EXPECTED_MAIN_ADULTS,
           "field_populations":32,"transects":5,
           "offspring_phenotypes_author_reported":2581,
           "field_maternal_plants_reported":382,
       },
       "source_spatial_model_metadata":{
           "demes_per_transect":36,"n_transects":5,
           "light_cline_km":s["original_published_light_cline_width_km"],
           "red_cline_km":s["original_published_red_cline_width_km"],
           "pollen_dispersal_variance_km2":s["source_supplement_pollen_kernel_variance_km2"],
           "weak_selection_one_locus_reconstructed":False,
           "gene_flow_parameter_newly_estimated_from_S1":False,
       },
       "author_published_fitness_per_region":f,
       "author_original_table1_locus_selection":m["original_table1_locus_selection"],
       "author_original_table1_epistasis":m["epistasis_original_table1"],
       "independent_arithmetic_on_original_published_epistasis":{
            "allopatry":round(eps["allopatry"],12),
            "sympatry":round(eps["sympatry"],12),
            "sympatry_minus_allopatry":round(eps["sympatry"]-eps["allopatry"],12),
            "not_a_new_MCMC_posterior_estimate":True,
       },
       "fitness_epistasis_changes_sign_across_regions":eps["allopatry"]>0>eps["sympatry"],
       "dark_red_relative_to_light_blue_by_environment":{
            k: f[k]["DR"]/f[k]["LB"] for k in ("allopatry","sympatry")
       },
       "regionally_normalized_fitness_as_common_absolute_seed_numbers":False,
       "different_2012_and_2014_hue_estimands":True,
       "2012_hue_field_hybridization_effect_not_detected":True,
       "2014_spatial_hue_allele_selection_author_inferred":True,
       "source_2014_intensity_coefficient_not_2012_viability_proxy_s":True,
       "provisional_original_S1_text_transcription":typed,
       "source_denominator_reconciliation":{
            "state":"SOURCE_ADULT_DENOMINATOR_NOT_RECONCILED",
            "main_paper_claim_adults":EXPECTED_MAIN_ADULTS,
            "S1_digitized_count":n,
            "main_minus_S1_digitized":difference,
            "full_original_2014_survey_reproduced":False,
            "original_S1_PDF_page_image_independently_checked":False,
            "not_claimed_original_publication_error":True,
            "requires_official_supplement_image_or_original_csv_before_inferential_refit":True,
       },
       "cross_sectional_population_and_one_offspring_generation_source":True,
       "new_field_evolutionary_time_series_fitted":False,
       "independent_full_2014_MCMC_reestimated":False,
       "independent_offspring_tableS2_reparsed":False,
       "new_parent_offspring_likelihood_fit":False,
       "new_original_selection_coefficient_estimate":None,
       "not_independent_repeat_of_v40_2012_hybrid_viability_selection":True,
       "only_original_author_reported_2014_model_inference":True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--expected",type=Path)
    a=ap.parse_args()
    j=review(json.loads(SOURCE.read_text(encoding="utf-8")))
    if a.expected and json.loads(a.expected.read_text(encoding="utf-8"))!=j:
        raise ValueError("V41_FROZEN_SOURCE_REPLAY_DRIFT")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(j,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("ANTO_V41_2014_SOURCE_SPATIAL_SELECTION_AUDIT",json.dumps({
        "n_populations":j["provisional_original_S1_text_transcription"]["n_populations_digitized"],
        "n_published_main_adults":j["source_denominator_reconciliation"]["main_paper_claim_adults"],
        "n_tableS1_digitized":j["source_denominator_reconciliation"]["S1_digitized_count"],
        "unreconciled":j["source_denominator_reconciliation"]["main_minus_S1_digitized"],
        "epistasis":j["independent_arithmetic_on_original_published_epistasis"],
        "fitness":j["author_published_fitness_per_region"],
        "no_new_MCMC_fit":True,
    },sort_keys=True))


if __name__=="__main__":
    main()
