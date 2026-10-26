"""ANTO v20: exact prospective function-assay partial-identification design.

No actual new functional assays are present. Ho & Smith 2016's 22
Iochrominae taxa are a SOURCE PHENOTYPE and deposited-accession panel,
not independently replicated floral regain events. The source CD S
annotation and sequence presence are *not* functional capacity assays.

Question: after genuinely validated taxon-level binary FUNCTIONAL-toolkit
assays are performed (hypothetically), at what observed counts would the
finite-panel difference in competence prevalence be sign-identified?

This is exact finite-population partial identification, not a causal
pollinator effect, heritability, fixation, a dated transition or a real
re-expression probability. Future actual specimens and permits are outside
this ANTO theoretical design; any taxon-level phenotype endpoint must have
independent provenance and an appropriate within-taxon sampling protocol.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from anto_afterlife_empirical_admissibility_v19 import (
    SOURCE_JOIN, STATUS, contingency, read, record_status
)

GROUPS=("UNPIGMENTED","PIGMENTED")
GROUP_SIZE={"UNPIGMENTED":7,"PIGMENTED":15}
SOURCE_STATUS={
    STATUS[0]:"DFR_SOURCE_CDS_FEATURE_PRESENT_NOT_FUNCTION",
    STATUS[1]:"DFR_SOURCE_RECORD_NO_CDS_FEATURE_NOT_GENE_LOSS",
    STATUS[2]:"DFR_NOT_IN_NEW_KT_SOURCE_RANGE_NOT_GENE_LOSS",
}


def _integer(x,name):
    if type(x) is not int or x<0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return x


def function_bound(
    unpigmented_tested:int, unpigmented_positive:int,
    pigmented_tested:int, pigmented_positive:int,
    *, n_unpigmented:int=7, n_pigmented:int=15
)->dict:
    """Sharp binary finite-panel bounds IF independent functional calls exist.

    Not a sampling confidence interval. Unknown taxon statuses can freely
    vary 0/1, without imputing 0 for missing gene/CDS source annotation.
    Counts are HYPOTHETICAL until validated comparable organ/stage and
    appropriately replicated taxon-level functional assays are actually done.
    """
    variables={
        "unpigmented_tested":unpigmented_tested,
        "unpigmented_positive":unpigmented_positive,
        "pigmented_tested":pigmented_tested,
        "pigmented_positive":pigmented_positive,
        "n_unpigmented":n_unpigmented,
        "n_pigmented":n_pigmented,
    }
    for name,x in variables.items():
        _integer(x,name)
    if n_unpigmented<1 or n_pigmented<1:
        raise ValueError("Both fixed groups require at least one taxon")
    if not (unpigmented_positive<=unpigmented_tested<=n_unpigmented
            and pigmented_positive<=pigmented_tested<=n_pigmented):
        raise ValueError("Positive functional assays cannot exceed examined taxa")
    # All unmeasured taxa may be positive or negative, independently.
    u_low=F(unpigmented_positive,n_unpigmented)
    u_high=F(unpigmented_positive+n_unpigmented-unpigmented_tested,
             n_unpigmented)
    p_low=F(pigmented_positive,n_pigmented)
    p_high=F(pigmented_positive+n_pigmented-pigmented_tested,n_pigmented)
    difference_low=u_low-p_high
    difference_high=u_high-p_low
    if difference_low>0:
        result="UNPIGMENTED_GT_PIGMENTED_CERTIFIED"
    elif difference_high<0:
        result="UNPIGMENTED_LT_PIGMENTED_CERTIFIED"
    elif difference_low==difference_high==0:
        result="EQUALITY_CERTIFIED_FOR_THIS_FINITE_PANEL"
    else:
        result="SIGN_NOT_IDENTIFIED"
    return {
        "group_denominators":[n_unpigmented,n_pigmented],
        "taxa_measured":[unpigmented_tested,pigmented_tested],
        "measured_positive":[unpigmented_positive,pigmented_positive],
        "unpigmented_prevalence_bound":[str(u_low),str(u_high)],
        "pigmented_prevalence_bound":[str(p_low),str(p_high)],
        "unpigmented_minus_pigmented_bound":[str(difference_low),str(difference_high)],
        "finite_panel_sign_decision":result,
        "sharp_given_binary_validated_taxon_level_assay":True,
        "does_not_identify_ancestral_regain_heritability_fitness_or_phylogenetic_causality":True,
    }


def earliest_possible_sign_certificate(n_unpigmented=7,n_pigmented=15):
    """Best-case *theoretical* minimum assays for a STRICT directional sign.

    All assayed nonpigmented taxa are hypothetically competent (positive);
    all assayed pigmented taxa are hypothetically not competent (negative).
    An opposite-sign certificate is obtained by swapping all outcomes.
    This is an existence witness, not a guarantee for unknown data and not
    an optimal real taxon-selection plan in a phylogenetically correlated clade.
    """
    _integer(n_unpigmented,"n_unpigmented")
    _integer(n_pigmented,"n_pigmented")
    if not n_unpigmented or not n_pigmented:
        raise ValueError("Need two nonempty groups")
    candidates=[]
    for u in range(n_unpigmented+1):
        for p in range(n_pigmented+1):
            decision=function_bound(u,u,p,0,n_unpigmented=n_unpigmented,
                                    n_pigmented=n_pigmented)
            if decision["finite_panel_sign_decision"]=="UNPIGMENTED_GT_PIGMENTED_CERTIFIED":
                candidates.append((u+p,u,p))
    if not candidates:
        raise ArithmeticError("Strict directional certification should be feasible")
    best=min(k for k,_,_ in candidates)
    allocation=sorted((u,p) for k,u,p in candidates if k==best)
    return {
        "minimum_hypothetical_taxon_functional_calls":best,
        "best_case_allocations_unpigmented_and_pigmented":[list(t) for t in allocation],
        "strict_certificate_requires_both_groups":True,
        "not_a_guaranteed_sample_size_or_power_analysis":True,
        "assumes_perfectly_qualified_binary_functional_taxon_assays":True,
    }


def build_taxon_ledger(source:dict):
    """Preserve all 22 source taxa; fill no missing function with a zero."""
    tab=contingency(source)  # Reuse v19 frozen source membership checks.
    if {g:sum(tab[g].values()) for g in GROUPS}!=GROUP_SIZE:
        raise ValueError("Original 7/15 Ho–Smith source membership changed")
    result=[]
    for row in sorted(source["taxon_gene_records"],
                      key=lambda r:r["current_taxon_name"]):
        if row.get("functional_joint_toolkit")!="UNMEASURED":
            raise ValueError("Source CDS cannot be promoted to plant functional competence")
        if row.get("direct_flower_reexpression")!="UNMEASURED":
            raise ValueError("Source phenotype cannot be promoted to direct regain")
        result.append({
            "taxon":row["current_taxon_name"],
            "original_2016_taxon":row["original_2016_figure_taxon"],
            "original_flower_pigment_category":row["trait_state_original_fig1"],
            "source_DFR_annotation_class":SOURCE_STATUS[record_status(row)],
            "DFR_new_KT_accession":row["DFR_accession_KT_registered"],
            "taxonomy_source_join":row["taxonomic_join_type"],
            "native_flower_joint_function":"UNMEASURED",
            "floral_induction_response":"UNMEASURED",
            "native_ancestral_regain_vs_donor_origin":"UNRESOLVED",
            "transgenerational_stability":"UNMEASURED",
            "trait_specific_reproduction":"UNMEASURED",
            "no_new_taxon_function_assay_admitted":True,
        })
    if len(result)!=22 or len({x["taxon"] for x in result})!=22:
        raise ArithmeticError("Taxon-level ledger must preserve 22 independent names")
    return result


def design_receipt(source:dict)->dict:
    ledger=build_taxon_ledger(source)
    count={g:{x:0 for x in SOURCE_STATUS.values()} for g in GROUPS}
    for record in ledger:
        count[record["original_flower_pigment_category"]][record["source_DFR_annotation_class"]]+=1
    if sorted(count["UNPIGMENTED"].values())!=[2,2,3] or sorted(count["PIGMENTED"].values())!=[0,6,9]:
        raise ArithmeticError("Source metadata coverage has unexpectedly changed")
    null=function_bound(0,0,0,0)
    best=earliest_possible_sign_certificate()
    # The following finite-population witness outcomes are artificial.
    scenarios={
        "actual_no_new_functional_assays":null,
        "synthetic_best_case_7U_positive_1P_negative":function_bound(7,7,1,0),
        "synthetic_best_case_7U_negative_1P_positive":function_bound(7,0,1,1),
        "synthetic_incomplete_no_sign":function_bound(4,2,4,2),
        "synthetic_full_panel_positive_group_difference":function_bound(7,5,15,6),
        "synthetic_full_panel_all_nonfunctional":function_bound(7,0,15,0),
    }
    if null["unpigmented_minus_pigmented_bound"]!=["-1","1"]:
        raise ArithmeticError("Zero biological assays cannot identify function")
    if best["minimum_hypothetical_taxon_functional_calls"]!=8:
        raise ArithmeticError("Unexpected mathematical best-case certificate")
    return {
        "schema":"anto_v20_functional_assay_identifiability_and_taxon_enrollment",
        "date":"2026-10-26",
        "source":"Original Ho and Smith 2016 fig.1 taxon/GenBank crosswalk, verified in ANTO v19",
        "source_reference":"results/iochrominae_2016_trait_KT_source_join_20261019.json",
        "original_source_taxa":22,
        "original_source_taxa_by_flower_pigment":dict(GROUP_SIZE),
        "source_annotation_coverage_not_function":count,
        "new_functional_assays_or_regain_events_added":0,
        "new_multigeneration_reexpression_fitness_experiments":0,
        "minimal_functional_assay_definition":"Functionally validated jointly required ancestral/potential pathway modules measured at appropriate flower organ/stage with a prespecified taxon-level sampling and comparison protocol; a CDS feature alone does not qualify",
        "target_estimand_if_future_assays_exist":"Finite 22-taxon panel prevalence difference in validated flower functional competence, NOT the frequency/rate of evolutionary regained traits",
        "entire_clade_phylogenetic_nonindependence_not_resolved":True,
        "genetic_retention_vs_current_flower_expression_vs_heritable_regain_are_separate_endpoints":True,
        "functional_assay_status":"NO_QUALIFYING_REAL_ASSAYS_ADMITTED",
        "identification_if_future_direct_taxon_assays_exist":scenarios,
        "minimum_best_case_design_bound":best,
        "taxon_enrollment":ledger,
        "no_ancestral_state_polarity_measured":True,
        "no_direct_function_or_reexpression_measured_in_this_KT_audit":True,
        "no_cross_angiosperm_rate_or_predictive_clade_holdout_estimated":True,
        "Cirsium_organismal_empirical_owner":"zuizui0223/aza3",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    parser.add_argument("--expected",type=Path)
    args=parser.parse_args()
    receipt=design_receipt(read(SOURCE_JOIN))
    if args.expected is not None:
        if read(args.expected)!=receipt:
            raise RuntimeError("v20 functional-assay design frozen receipt mismatch")
    s=json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n"
    if args.output is None:
        print(s)
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding="utf-8")


if __name__=="__main__":
    main()
