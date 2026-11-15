"""ANTO v35 published original-source evidence-state comparison.

No raw HPLC/qPCR retrieved, no fabricated common quantitative gene metric.
Findings belong to the ORIGINAL published experiments; our new result is
a logical measurement-identifiability audit of source-backed claims.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"data/external/anto_v35_white_flower_cross_family_published_source_registry_20261115.json"
DOIS={
    "hirado_azalea_2020":"10.2503/hortj.UTD-142",
    "ipomoea_lacunosa_2020":"10.1371/journal.pone.0231263",
}
STATES={
    "hirado_azalea_2020":"F3HPRIME_F3H5PRIME_DFR_ANS_TRANSCRIPTS_DETECTED_IN_WHITE",
    "ipomoea_lacunosa_2020":"CHS_DFR_MYB1_STRONGLY_DOWNREGULATED_IN_WHITE",
}


def review(source):
    if source.get("schema")!="anto_v35_published_original_white_flower_expression_hplc_source_registry":
        raise ValueError("SOURCE_REGISTRY_SCHEMA_DRIFT")
    rows=source.get("external_independent_primary_study_programs")
    if not isinstance(rows,list) or len(rows)!=2 or len({r.get("id") for r in rows})!=2:
        raise ValueError("INDEPENDENT_PRIMARY_SOURCE_DENOMINATOR_DRIFT")
    if {r["id"]:r["doi"] for r in rows}!=DOIS:
        raise ValueError("SOURCE_DOI_DRIFT")
    if any(r["visible_flower_class"]!="WHITE" for r in rows):
        raise ValueError("VISIBLE_PHENOTYPE_NOT_MATCHED")
    if {r["id"]:r["transcription_state"] for r in rows}!=STATES:
        raise ValueError("ORIGINAL_TRANSCRIPTION_STATE_DRIFT")
    if {r["plant_family"] for r in rows}!={"Ericaceae","Convolvulaceae"}:
        raise ValueError("NOT_INDEPENDENT_FAMILY_GROUPS")
    az=next(r for r in rows if r["id"]=="hirado_azalea_2020")
    ip=next(r for r in rows if r["id"]=="ipomoea_lacunosa_2020")
    if (az["measured_genes"]!=["F3'H","F3'5'H","DFR","ANS"]
            or az["absent_from_measured_gene_set"]!=["FLS"]):
        raise ValueError("UNMEASURED_GENE_CANNOT_BE_PROMOTED")
    if "FLS" in az["measured_genes"]:
        raise ValueError("UNMEASURED_GENE")
    if (az["reported_hplc_wild_samples"]!=6 or
            az["reported_hplc_cultivar_samples"]!=25 or
            az["reported_qpcr_wild_samples"]!=6 or
            az["reported_qpcr_cultivar_samples"]!=18 or
            len(az["white_hirado_cultivars_author_named"])!=4):
        raise ValueError("SOURCE_DENOMINATOR_DRIFT")
    if az["white_anthocyanidin_aglycone_exception"]!="SMALL_DETECTABLE_PEAKS_REPORTED_FOR_SHIRO_RYUKYU":
        raise ValueError("ANTHOCYANIDIN_TRACE_MUST_NOT_BE_CALLED_ABSENT")
    if (ip["measured_genes"]!=["CHS","DFR","MYB1"]
            or ip["reported_field_qpcr_white_individuals"]!=4
            or ip["reported_field_qpcr_purple_individuals"]!=4
            or ip["reported_hplc_white_I_lacunosa_individuals"]!=3
            or ip["reported_hplc_purple_I_lacunosa_individuals"]!=3
            or ip["reported_hplc_purple_I_cordatotriloba_individuals"]!=2):
        raise ValueError("SOURCE_IPOMOEA_SAMPLE_OR_GENE_DRIFT")
    DFR_STATES={
        "hirado_azalea_2020":
            "WHITE_DFR_DETECTED_SIMILAR_TO_COLORED_PUBLISHED_DESCRIPTION",
        "ipomoea_lacunosa_2020":
            "WHITE_DFR_OVER_100_FOLD_DOWN_VS_PURPLE_PUBLISHED_DESCRIPTION",
    }
    if {x["id"]:x.get("common_DFR_relative_to_colored_reference") for x in rows}!=DFR_STATES:
        raise ValueError("SOURCE_COMMON_DFR_RELATIVE_CONTRAST_DRIFT")
    if (any("DFR" not in x["measured_genes"] for x in rows)
            or source.get("scientific_ceilings",{}).get(
                "same_DFR_gene_reported_but_assays_and_reference_taxa_not_harmonized") is not True):
        raise ValueError("COMMON_DFR_ASSAY_ADMISSIBILITY_DRIFT")
    if ip["approx_expression_change_white_relative_purple"]!={
            "CHS":"APPROXIMATELY_100_FOLD_LOWER",
            "DFR":"GREATER_THAN_100_FOLD_LOWER",
            "MYB1":"APPROXIMATELY_100_FOLD_LOWER"}:
        raise ValueError("SOURCE_APPROXIMATE_EXPRESSION_REPORTING_DRIFT")
    if any(r["hplc_rna_same_biological_individuals_confirmed"] for r in rows):
        raise ValueError("FORGED_JOINT_MEASUREMENT")
    if any(r["enzyme_function_or_ancestral_reversal_measured"] for r in rows):
        raise ValueError("UNSUPPORTED_FUNCTIONAL_ASSAY")
    if any(r["new_species_holdout_prediction_from_author_raw_data"] for r in rows):
        raise ValueError("FORGED_INDEPENDENT_OUT_OF_SAMPLE_MODEL")
    ctx=source.get("context_not_new_independent_primary_qpcr_hplc_reanalysis")
    if not isinstance(ctx,list) or len(ctx)!=2 or {x["id"] for x in ctx}!={"iochrominae_2018_2019","petunieae_2023"}:
        raise ValueError("PRIOR_CONTEXT_SOURCE_SET_DRIFT")
    if any(x["raw_qpcr_HPLC_reparsed_in_v35"] for x in ctx):
        raise ValueError("FORGED_PREVIOUS_SOURCE_RAW_REANALYSIS")
    dryad=next(x for x in ctx if x["id"]=="iochrominae_2018_2019")
    if dryad["original_dryad_RAR_access"]!="PREVIOUS_V23_HTTP_401_403":
        raise ValueError("DRYAD_ACCESS_CLAIM_DRIFT")
    if not all(x["no_extra_independent_program_count"] for x in ctx):
        raise ValueError("STUDY_SAMPLE_PSEUDOREPLICATION")
    ceilings=source.get("scientific_ceilings",{})
    if not ceilings or any(x is not True for x in ceilings.values()):
        raise ValueError("SOURCE_CLAIM_FIREWALL_DRIFT")
    return {
        "schema":"anto_v35_original_published_white_flower_state_nonuniqueness_audit",
        "original_source_input_schema":source["schema"],
        "study_observation_grains":{
            "hirado_azalea_2020":{
                "HPLC_samples":az["reported_hplc_wild_samples"]+az["reported_hplc_cultivar_samples"],
                "qPCR_samples":az["reported_qpcr_wild_samples"]+az["reported_qpcr_cultivar_samples"],
                "n_author_named_white_Hirado_cultivars":len(az["white_hirado_cultivars_author_named"]),
                "distinct_species_count_not_31":True,
                "individual_HPLC_RNA_join_not_verified":True,
            },
            "ipomoea_lacunosa_2020":{
                "field_qpcr_individuals":ip["reported_field_qpcr_white_individuals"]+ip["reported_field_qpcr_purple_individuals"],
                "I_lacunosa_HPLC_individuals":ip["reported_hplc_white_I_lacunosa_individuals"]+ip["reported_hplc_purple_I_lacunosa_individuals"],
                "additional_I_cordatotriloba_HPLC_individuals":ip["reported_hplc_purple_I_cordatotriloba_individuals"],
                "individual_HPLC_RNA_join_not_verified":True,
            }
        },
        "n_external_independent_study_programs":2,
        "n_independent_original_family_groups":len({r["plant_family"] for r in rows}),
        "independent_family_names":sorted({r["plant_family"] for r in rows}),
        "author_published_dois":sorted(r["doi"] for r in rows),
        "n_context_prior_studies_not_new_reanalyses":len(ctx),
        "all_primary_programs_share_visible_white_state":all(r["visible_flower_class"]=="WHITE" for r in rows),
        "distinct_documented_white_transcription_state_classes":len({r["transcription_state"] for r in rows}),
        "white_colour_does_not_uniquely_identify_regulatory_state":
            len({r["transcription_state"] for r in rows})>1,
        "same_visible_white_class_different_expression_reporting":True,
        "observed_pairwise_evidence_states":[{
            "doi":r["doi"],
            "family":r["plant_family"],
            "visible_state":r["visible_flower_class"],
            "gene_expression_reporting":r["transcription_state"],
            "published_assay_genes":r["measured_genes"],
            "pigment_reporting":r["white_anthocyanin_observation"],
            "not_functional_enzyme_validation":True,
        } for r in rows],
        "common_DFR_is_measured_in_both_original_programs":True,
        "source_reported_DFR_relative_contrasts_differ":True,
        "common_DFR_within_study_descriptions":DFR_STATES,
        "absolute_DFR_expression_is_cross_study_comparable":False,
        "common_DFR_relative_contrast_does_not_prove_nonoverlap_of_absolute_expression":True,
        "original_ancestral_gene_function_demonstrated":False,
        "intergenerational_natural_regain_demonstrated":False,
        "same_individual_HPLC_RNA_matched_comparison_computed":False,
        "raw_joint_individual_HPLC_qPCR_data_unavailable_in_registry":True,
        "gene_expression_units_and_pigment_methods_not_poolable":True,
        "independent_quantitative_meta_analysis_possible":False,
        "numeric_cross_study_effect_size":None,
        "new_source_raw_expression_rows_fitted":0,
        "ANTOv33_exact_F3Hprime_F3H5prime_FLS_model_externally_replayed":False,
        "white_state_nonuniqueness_is_observational_across_studies_not_ancestral_regain":True,
        "source_provenance_original_2020_articles_not_original_assay_files":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    output=review(json.loads(INPUT.read_text(encoding="utf-8")))
    if args.expected and json.loads(args.expected.read_text(encoding="utf-8"))!=output:
        raise ValueError("V35_FROZEN_EVIDENCE_STATE_RESULT_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ANTO_V35_SOURCE_BOUNDED_EXTERNAL_STATE_AUDIT",
          json.dumps({"n_external_original_programs":output["n_external_independent_study_programs"],
                      "families":output["independent_family_names"],
                      "white_regulatory_states":output["distinct_documented_white_transcription_state_classes"],
                      "visible_state_does_not_uniquely_specify_transcripts":output["white_colour_does_not_uniquely_identify_regulatory_state"],
                      "independent_quantitative_HPLC_expression_model":False,
                      "no_ancestral_regain":True},sort_keys=True))


if __name__=="__main__":
    main()
