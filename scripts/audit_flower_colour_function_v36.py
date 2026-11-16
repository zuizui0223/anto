"""Audit published ORIGINAL experiments behind white floral phenotypes.

Independent source programs (N=3), not independent genotype events.
Numeric sample denominators represent original author experiments, not
new raw-data fits, or prospective natural ancestral phenotype reversions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v36_independent_original_white_functional_programs_20261116.json"
DOI={
    "ipomoea_2020":"10.1371/journal.pone.0231263",
    "iochroma_MYBL1_2018":"10.1111/nph.14830",
    "azalea_2020":"10.2503/hortj.UTD-142",
}
FAMILY={
    "ipomoea_2020":"Convolvulaceae",
    "iochroma_MYBL1_2018":"Solanaceae",
    "azalea_2020":"Ericaceae",
}

def fail_if(condition, message):
    if condition:
        raise ValueError(message)


def review(source):
    fail_if(source.get("schema")!="anto_v36_independent_original_white_flower_biochemical_genetic_intervention_registry",
            "SOURCE_ORIGINAL_SCHEMA_DRIFT")
    programs=source.get("original_programs")
    fail_if(not isinstance(programs,list) or len(programs)!=3
            or len({p.get("id") for p in programs})!=3,
            "ORIGINAL_PROGRAM_SET_DRIFT")
    observed={p["id"]:p for p in programs}
    fail_if({k:p["doi"] for k,p in observed.items()}!=DOI,
            "SOURCE_DOI_IDENTITY_DRIFT")
    fail_if({k:p["family"] for k,p in observed.items()}!=FAMILY
            or any(p.get("independent_paper_program_count")!=1 for p in programs),
            "SOURCE_ORIGINAL_FAMILY_GRAIN_DRIFT")
    fail_if(not source["no_current_original_species_white_to_pigmented_rescue_in_any_program"]
            or not source["no_natural_ancestral_loss_then_regain_trajectory_measured"]
            or not source["do_not_mix_Nicotiana_transgenic_lines_with_original_species_reversions"]
            or not source["do_not_pool_crossstudy_gene_effects_or_absorbance_with_HPLC"],
            "ORIGINAL_STUDY_INFERENCE_FIREWALL_MISSING")

    ip=observed["ipomoea_2020"]
    t=ip["chemistry"]
    fail_if(t["white_individuals"]!=3 or t["purple_individuals"]!=3 or
            t["white_corolla_anthocyanidins_detected"] is not False or
            t["white_stem_anthocyanidins_detected"] is not True or
            t["purple_corolla_anthocyanidins_detected"] is not True or
            t["purple_stem_anthocyanidins_detected"] is not True or
            t["white_corolla_failure_does_not_equal_whole_plant_loss"] is not True,
            "TISSUE_BIOCHEMISTRY_DRIFT")
    e=ip["expression"]
    fail_if(e["field_white_individuals"]!=4 or e["field_purple_individuals"]!=4
            or e["genes"]!=["CHS","DFR","MYB1"]
            or e["white_vs_purple_source_verbal_effects"]!={"CHS":"ABOUT_100_FOLD_LOWER",
                 "DFR":"MORE_THAN_100_FOLD_LOWER","MYB1":"ABOUT_100_FOLD_LOWER"}
            or e["same_individual_as_HPLC_not_demonstrated"] is not True
            or e["direct_molecular_gene_editing"] is not False,
            "IPOMOEA_EXPRESSION_SAMPLE_OR_CLAIM_DRIFT")
    c=ip["crossing"]
    x=c["published_individually_genotyped_MYB_offspring"]
    fail_if(x["family_1_3_B2"]!=41 or x["family_4_F3"]!=12
            or x["purple_phenotype_with_purple_MYB_marker"]!=23
            or x["white_phenotype_with_white_MYB_marker"]!=30
            or x["discordant_marker_phenotype_pairs"]!=0
            or x["DFR_marker_did_not_cosegregate"] is not True
            or x["family_1_3_B2"]+x["family_4_F3"]!=
            x["purple_phenotype_with_purple_MYB_marker"]+
            x["white_phenotype_with_white_MYB_marker"]+
            x["discordant_marker_phenotype_pairs"],
            "ORIGINAL_TABLE3_GENETIC_DENOMINATOR_DRIFT")
    s2=c["separate_original_S2_throat_phenotypes"]
    fail_if(s2["total"]!=732 or s2["purple"]!=530 or s2["white"]!=202
            or s2["total"]!=s2["purple"]+s2["white"]
            or c["allele_specific_expression_F1_individuals"]!=3
            or c["direct_native_white_to_purple_restoration_experiment"] is not False
            or c["MYB_cis_regulatory_interpretation"]!="ORIGINAL_AUTHORS_SUPPORTED_NOT_DIRECT_SITE_EDIT"
            or ip["heterologous_regulator_perturbation"] is not None,
            "IPOMOEA_CROSSING_ORIGINAL_DATA_DRIFT")

    io=observed["iochroma_MYBL1_2018"]
    b=io["crossing"]["backcross_to_blue"]
    w=io["crossing"]["backcross_to_white"]
    fail_if(b!={"blue":35,"pale_F1like":38} or w!={"white":8,"pale_F1like":11}
            or io["crossing"]["original_floral_color_segregation"] is not True
            or io["crossing"]["native_single_causal_MYBL1_variant_genetically_edited"] is not False
            or io["crossing"]["pooled_genotype_association_not_individual_site_rescue"] is not True,
            "IOCHROMA_BACKCROSS_DENOMINATOR_DRIFT")
    fail_if(io["expression"]["upregulated_repressor"]!="MYBL1"
            or io["expression"]["downstream_anthocyanin_genes_suppressed_in_white_and_pale_pools"] is not True
            or io["expression"]["not_full_individual_gene_to_HPLC_matrix"] is not True
            or io["expression"]["native_causal_nucleotide_directly_identified"] is not False,
            "IOCHROMA_EXPRESSION_CAUSALITY_DRIFT")
    h=io["heterologous_regulator_perturbation"]
    fail_if(h is None or h["donor"]!="Iochroma loxense"
            or h["transformed_host"]!="Nicotiana tabacum"
            or h["transgene"]!="MYBL1 R3-MYB repressor"
            or h["transformed_lines_grown"]!=20
            or h["transformed_lines_white_author_reported"] is not True
            or h["original_host_wildtype_comparator"] is not True
            or h["author_reported_floral_absorbance_fold_lower_mean"]!=14
            or h["author_reported_floral_absorbance_fold_lower_min"]!=11
            or h["author_reported_floral_absorbance_fold_lower_max"]!=20,
            "Nicotiana_TRANSFORMATION_OR_ASSAY_DRIFT")
    fail_if(h["not_endogenous_Iochroma_rescue"] is not True
            or h["no_evidence_of_restored_ancestral_blue_flower"] is not True,
            "FORGED_NATIVE_RESCUE")
    fail_if(h["assay"]!="530-nm light absorbance per gram fresh corolla tissue"
            or h["functional_conclusion"]!="HETEROLOGOUS_SUFFICIENCY_FOR_PIGMENT_REPRESSION",
            "TRANSFORMATION_READOUT_MISCLASSIFIED")

    az=observed["azalea_2020"]
    a=az["chemistry"]
    q=az["expression"]
    fail_if(a["wild_origin_samples"]!=6 or a["cultivar_samples"]!=25
            or len(a["named_white_cultivars"])!=4
            or a["not_31_independent_species"] is not True
            or a["white_cultivars_flavonol_present"] is not True
            or a["reported_white_anthocyanidin_precursor_exception"]!="SHIRO_RYUKYU_AND_WHITE_HIRADO_REPORTED",
            "AZALEA_HPLC_DENOMINATOR_OR_PIGMENT_DRIFT")
    fail_if(q["wild_origin_samples"]!=6 or q["cultivar_samples"]!=18
            or q["white_detected_genes"]!=["F3'H","F3'5'H","DFR","ANS"]
            or q["FLS_assayed"] is not False
            or q["RNA_HPLC_individual_pairing_not_verified"] is not True,
            "AZALEA_QPCR_DENOMINATOR_OR_GENE_DRIFT")
    fail_if(q["not_original_enzyme_functional_test"] is not True,
            "PCR_IS_NOT_FUNCTIONAL_ASSAY")
    fail_if(az["heterologous_regulator_perturbation"] is not None or az["crossing"] is not None,
            "AZALEA_NO_FUNCTIONAL_PERTURBATION")
    fail_if(az["already_in_v35"] is not True,
            "DUPLICATE_AZALEA_PROGRAM_NOT_ADMITTED")

    return {
        "schema":"anto_v36_original_study_white_flower_functional_evidence_ladder",
        "source_registry":"data/external/anto_v36_independent_original_white_functional_programs_20261116.json",
        "original_study_programs":len(programs),
        "distinct_plant_families":len({p["family"] for p in programs}),
        "original_dois_sorted":sorted(p["doi"] for p in programs),
        "quantitatively_poolable_study_effects":False,
        "no_source_original_raw_expression_HPLC_refit":True,
        "cross_tissue_white_flower_stem_anthocyanidin":{
            "source":ip["doi"],
            "white_plants_assayed":t["white_individuals"],
            "white_flower_detected":t["white_corolla_anthocyanidins_detected"],
            "white_stem_detected":t["white_stem_anthocyanidins_detected"],
            "purple_plants_assayed":t["purple_individuals"],
            "purple_flower_detected":t["purple_corolla_anthocyanidins_detected"],
            "purple_stem_detected":t["purple_stem_anthocyanidins_detected"],
            "native_white_flower_reactivation_not_tested":True,
        },
        "ipomoea_MYB_genotype_phenotype_original_table3":{
            "source":ip["doi"],
            "B2_genotyped":x["family_1_3_B2"],
            "F3_genotyped":x["family_4_F3"],
            "total_genotyped":x["family_1_3_B2"]+x["family_4_F3"],
            "concordant":x["purple_phenotype_with_purple_MYB_marker"]+
                         x["white_phenotype_with_white_MYB_marker"],
            "discordant":x["discordant_marker_phenotype_pairs"],
            "purple_marker_and_flower":x["purple_phenotype_with_purple_MYB_marker"],
            "white_marker_and_flower":x["white_phenotype_with_white_MYB_marker"],
            "DFR_marker_did_not_cosegregate":x["DFR_marker_did_not_cosegregate"],
            "not_an_editing_experiment":True,
        },
        "separate_Ipomoea_S2_throat_phenotyping":{
            "offspring":s2["total"],
            "purple":s2["purple"],
            "white":s2["white"],
            "not_genotyped_Table3_53":True,
        },
        "Iochroma_published_crossing_design":{
            "source":io["doi"],
            "backcross_blue_parent_offspring":b["blue"]+b["pale_F1like"],
            "backcross_white_parent_offspring":w["white"]+w["pale_F1like"],
            "all_crossing_offspring":sum(b.values())+sum(w.values()),
            "not_independent_mutational_origins":True,
        },
        "functional_interventions":{
            "Iochroma_MYBL1_heterologous_in_Nicotiana":{
                "source":io["doi"],
                "donor":h["donor"],"transgenic_host":h["transformed_host"],
                "source_transformed_lines_grown":h["transformed_lines_grown"],
                "author_reported_all_transgenic_lines_white":True,
                "mean_fold_lower_reported":h["author_reported_floral_absorbance_fold_lower_mean"],
                "line_fold_lower_reported_range":[
                    h["author_reported_floral_absorbance_fold_lower_min"],
                    h["author_reported_floral_absorbance_fold_lower_max"]],
                "absorbance_not_HPLC_concentration":True,
                "target_process":"PIGMENT_REPRESSION",
                "endogenous_Iochroma_anthocyanin_restoration":False,
                "functionally_active_regulator_not_direct_restoration":True
            }
        },
        "azalea_original_grain":{
            "source":az["doi"],
            "HPLC_wild_plus_cultivar":a["wild_origin_samples"]+a["cultivar_samples"],
            "qPCR_wild_plus_cultivar":q["wild_origin_samples"]+q["cultivar_samples"],
            "cultivars_not_independent_species":True,
            "native_gene_function_verified":False,
            "FLS_gene_assayed":False,
        },
        "source_evidence_levels":{
            "original_flower_chemical_assay":sorted([ip["id"],io["id"],az["id"]]),
            "original_flower_transcript_assay":sorted([ip["id"],io["id"],az["id"]]),
            "cross_tissue_white_flower_stem_biochemistry":[ip["id"]],
            "original_genotype_flower_cosegregation":[ip["id"]],
            "native_original_floral_colour_segmentation_backcrosses":sorted([ip["id"],io["id"]]),
            "heterologous_gene_pigment_suppression":[io["id"]],
            "native_white_flower_restoration":[],
            "natural_ancestral_reacquisition_with_persistence":[],
        },
        "native_corolla_reexpression_experiment_count":0,
        "direct_native_gene_edit_rescue_confirmed":False,
        "regain_probability_estimated":False,
        "cross_family_pooled_molecular_effect":None,
        "exact_v33_three_gene_FLS_model_replicated":False,
        "no_direct_cross_study_identical_gene_HPLC_prediction":True,
        "genotyped_offspring_not_independent_reversion_events":True,
        "assays_and_original_published_counts_not_independent_phylogenetic_replications":True,
        "white_phenotype_not_equivalent_to_whole_plant_pathway_loss":True,
        "heterologous_pigment_repression_not_endogenous_ancestral_regain":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--expected",type=Path)
    opts=parser.parse_args()
    result=review(json.loads(REGISTRY.read_text(encoding="utf-8")))
    if opts.expected and json.loads(opts.expected.read_text(encoding="utf-8"))!=result:
        raise ValueError("V36_FROZEN_ORIGINAL_AUTHOR_STUDY_EVIDENCE_DRIFT")
    opts.output.parent.mkdir(parents=True,exist_ok=True)
    opts.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ANTO_V36_ORIGINAL_EXPERIMENTAL_EVIDENCE",json.dumps({
        "original_papers":result["original_study_programs"],
        "white_flower_stem_functional_product":result["cross_tissue_white_flower_stem_anthocyanidin"]["white_stem_detected"],
        "genotype_color_concordant":result["ipomoea_MYB_genotype_phenotype_original_table3"]["concordant"],
        "genotype_color_discordant":result["ipomoea_MYB_genotype_phenotype_original_table3"]["discordant"],
        "heterologous_MYBL1_repressor_n":result["functional_interventions"]["Iochroma_MYBL1_heterologous_in_Nicotiana"]["source_transformed_lines_grown"],
        "heterologous_fold_decrease_reported_mean":result["functional_interventions"]["Iochroma_MYBL1_heterologous_in_Nicotiana"]["mean_fold_lower_reported"],
        "native_white_to_colour_restoration_performed":result["direct_native_gene_edit_rescue_confirmed"],
    },sort_keys=True))


if __name__=="__main__":
    main()
