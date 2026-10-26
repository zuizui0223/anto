"""Reconcile original Penstemon 2014 reported cross outcomes with ANTO afterlife concepts.

Wessinger & Rausher (2014), Evolution 68:1058-1070, Table 1, DOI
10.1111/evo.12340. Counts were published by the original authors.

This source is a BLUE-PURPLE pathway component loss in a still-RED
anthocyanin-containing lineage, NOT a white/achromatic flower transition.
Hybrid donor restoration is different from endogenous genetic backmutation.

This code does not infer natural gene flow, viable seed fitness, historical
flower regain rate, true cross-clade universality, or ANTO a/r/d/g.
"""
from fractions import Fraction as F
from copy import deepcopy
from pathlib import Path
import argparse
import json


def audit(source):
    if source.get("schema")!="anto_penstemon_barbatus_F3p5H_genetic_accessibility_source_table_v1":
        raise ValueError("Wrong paper-level measurement source")
    if source.get("source",{}).get("paper_doi")!="10.1111/evo.12340":
        raise ValueError("Mismatched external original paper")
    if "not a white flower" not in source.get("trait_target",""):
        raise ValueError("Do not misclassify ancestral blue lost as total floral pigment loss")
    expected={"FF","Ff","ff"}
    table=source.get("observed_F2_table",{})
    if set(table)!={"purple","red"}:
        raise ValueError("Two original flower phenotypes required")
    for state in table:
        if set(table[state])!=expected:
            raise ValueError("Original genotype classes missing")
        if any(type(v)!=int or v<0 for v in table[state].values()):
            raise ValueError("Nonnegative integer original source counts required")
    purple=sum(table["purple"].values())
    red=sum(table["red"].values())
    n=purple+red
    if (n,purple,red)!=(161,122,39):
        raise ValueError("Original published F2 denominator and phenotype totals drift")
    pff=table["purple"]["FF"]+table["red"]["FF"]
    pFf=table["purple"]["Ff"]+table["red"]["Ff"]
    pff_small=table["purple"]["ff"]+table["red"]["ff"]
    if (pff,pFf,pff_small)!=(57,80,24):
        raise ValueError("Original three genotype denominators drift")
    if (table["purple"]["FF"],table["purple"]["Ff"],table["purple"]["ff"])!=(57,65,0):
        raise ValueError("Published purple genotype counts drift")
    if (table["red"]["FF"],table["red"]["Ff"],table["red"]["ff"])!=(0,15,24):
        raise ValueError("Published red genotype counts drift")
    mechanistic=source.get("reported_source_molecular_results",{})
    if (mechanistic.get("P_barbatus_F3p5H_stop_codon_AA")!=407
        or mechanistic.get("stop_repair_only_restored_activity") is not False
        or mechanistic.get("multiple_function_disabling_changes_supported_by_source") is not True
        or mechanistic.get("unmeasured_actual_wild_gene_flow") is not True
        or mechanistic.get("no_direct_natural_pbarbatus_backmutation_reversal_observed") is not True):
        raise ValueError("Original source molecular nonclaims mutated")
    if source.get("user_discovered_new_empirical_evolutionary_transition") is not False:
        raise ValueError("The reported source article is not the user's new experiment")
    conditional={"FF":F(table["purple"]["FF"],pff),
                 "Ff":F(table["purple"]["Ff"],pFf),
                 "ff":F(table["purple"]["ff"],pff_small)}
    return {
        "schema":"anto_penstemon_2014_table1_exact_afterlife_distinction_receipt_v1",
        "original_paper_DOI":source["source"]["paper_doi"],
        "original_Dryad_DOI":source["source"]["data_doi"],
        "observed_F2_individuals":n,
        "phenotypes_source_still_red_or_blue_purple_not_all_white":True,
        "source_genotype_denominators":{"FF":pff,"Ff":pFf,"ff":pff_small},
        "source_purple_by_genotype":{"FF":57,"Ff":65,"ff":0},
        "source_red_by_genotype":{"FF":0,"Ff":15,"ff":24},
        "blue_purple_probability_conditional_on_genotype_fraction":
            {g:str(v) for g,v in conditional.items()},
        "heterozygote_red_probability":str(F(15,80)),
        "blue_purple_in_cross_offspring_total":122,
        "functional_donor_F_allele_but_still_red_Ff":15,
        "source_internal_F3p5H_repair_of_single_stop_failed_to_restore_in_vitro_activity":True,
        "external_functional_donor_allele_can_produce_blue_purple_in_hybrid_offspring":True,
        "difference_native_repair_vs_hybrid_inheritance_is_mechanistically_necessary":True,
        "plant_phenotype_after_regain_fitness_measured_here":False,
        "natural_introgression_frequency_or_permanence_estimated":False,
        "synthetic_ANTO_model_parameter_g_inferred":False,
        "species_wide_evolutionary_regain_probability_inferred":False,
        "original_source_results_only":True,
    }


def main():
    root=Path(__file__).resolve().parents[1]
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=root/"data/external/penstemon_2014_f3p5h_source_genetic_accessibility_v1.json")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    source=json.loads(args.source.read_text(encoding="utf-8"))
    report=audit(source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("PENSTEMON_ORIGINAL_TABLE_SOURCE_AUDIT_PASS",
          json.dumps({"F2_n":report["observed_F2_individuals"],
                      "heterozygote_red_fraction":report["heterozygote_red_probability"],
                      "no_natural_reversal_inferred":True},sort_keys=True))

if __name__=="__main__":
    main()
