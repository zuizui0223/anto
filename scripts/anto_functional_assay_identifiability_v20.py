"""ANTO v20: source-anchored prospective FUNCTION assay identifiability.

No actual complete-floral-toolkit assays are present in the 22 original
Iochrominae source taxa. Neither archived gene accessions nor CDS features
constitute joint floral enzymatic and regulatory FUNCTION measurements.

The calculations below are exact partial-identification bounds for a
HYPOTHETICAL binary, independently gold-standard-assayed capacity status
for each of the 22 named taxa. Positives/negatives in example scenarios are
deliberately invented mathematical observations, NEVER experimental results.
We do not estimate real trait re-expression, ancestral allele repair,
population persistence, pollinator adaptation, or a genus-wide effect.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

from anto_afterlife_empirical_admissibility_v19 import (
    SOURCE_JOIN, read, contingency, record_status, STATUS,
)

ROOT = Path(__file__).resolve().parents[1]
N_WHITE=7
N_PIG=15
PHENOTYPES=("UNPIGMENTED","PIGMENTED")
PROTOCOL = (
    "functional coding AND regulatory components assayed jointly in the same plant",
    "flower-specific organ and matched developmental stage independently recorded",
    "biochemical pigment pathway and enzyme/regulator function tested, not homology only",
    "genotype, cytotype, ancestry and introgression controls specified",
    "induction/sham and offspring inheritance separately assessed when claiming re-expression",
    "source, biological replicate, population and voucher identifiers verified",
)


def panel_from_original_source(join):
    """Read the 22 real original-source rows, preserving ALL 22 missing assays."""
    source_table=contingency(join)
    rows=[]
    for row in join["taxon_gene_records"]:
        rows.append({
            "taxon": row["current_taxon_name"],
            "original_taxon":row["original_2016_figure_taxon"],
            "flower_pigment_category":row["trait_state_original_fig1"],
            "source_DFR_metadata":record_status(row),
            "original_DFR_KT_accession":row["DFR_accession_KT_registered"],
            "flower_functional_joint_toolkit_assay":None,
            "assay_admitted":False,
            "historical_natural_trait_regain_in_same_lineage_assayed":False,
            "multi_generation_functional_retention_assayed":False,
            "source_observation_unit":"original named taxon, not an independent evolution event",
        })
    rows.sort(key=lambda x:x["taxon"])
    if len(rows)!=22 or len(set(x["taxon"] for x in rows))!=22:
        raise ValueError("Source taxon inventory changed")
    if Counter(x["flower_pigment_category"] for x in rows)!={"UNPIGMENTED":7,"PIGMENTED":15}:
        raise ValueError("Source phenotype stratification changed")
    for category in PHENOTYPES:
        for metadata in STATUS:
            if sum(x["flower_pigment_category"]==category and
                   x["source_DFR_metadata"]==metadata for x in rows)!=source_table[category][metadata]:
                raise ValueError("Source phenotype/CDS censoring table inconsistent")
    return rows


def _source_and_label_gate(panel, labels, *, synthetic):
    """Hard separation between a source observation and an invented oracle label.

    Labels must be an explicit dict taxon -> bool; even a complete set of
    labels is NOT an actual botanical assay. Use only synthetic=True.
    """
    if synthetic is not True:
        raise ValueError("Functional labels cannot be promoted from source metadata")
    if not isinstance(labels,dict):
        raise ValueError("Hypothetical oracle labels must be a dict")
    source_ids=[x["taxon"] for x in panel]
    if len(panel)!=22 or len(set(source_ids))!=22:
        raise ValueError("Full fixed 22 taxon denominator is required")
    if set(labels)-set(source_ids):
        raise ValueError("Functional label belongs to no registered source taxon")
    if any(type(x) is not bool for x in labels.values()):
        raise ValueError("Functional label must be a hypothetical Boolean, not an accession flag")
    if any(x["flower_functional_joint_toolkit_assay"] is not None or
           x["assay_admitted"] is not False for x in panel):
        raise ValueError("Original source record contains NO admitted functional outcome")
    return labels


def sharp_bounds(panel, labels=None, *, synthetic=True):
    """Sharp finite-22-taxon proportion bounds under unrestricted missingness.

    The only assumption in hypothetical labelled scenarios is that oracle
    labels accurately identify joint floral competence in the full taxon
    scope. With no labels, true functional state of EVERY taxon is unknown.
    Missing outcomes may independently take either value, so endpoints
    of the interval are attainable by construction. No exchangeability,
    phylogenetic independence or missing-at-random assumption is made.
    """
    labels={} if labels is None else labels
    _source_and_label_gate(panel,labels,synthetic=synthetic)
    groups={}
    for state,n in (("UNPIGMENTED",N_WHITE),("PIGMENTED",N_PIG)):
        taxa=[x["taxon"] for x in panel if x["flower_pigment_category"]==state]
        if len(taxa)!=n:
            raise ValueError("Source pigment group denominator changed")
        positive=sum(x in labels and labels[x] is True for x in taxa)
        negative=sum(x in labels and labels[x] is False for x in taxa)
        unknown=n-positive-negative
        groups[state]={
            "source_taxa":n,
            "hypothetically_assayed_positive":positive,
            "hypothetically_assayed_negative":negative,
            "not_functionally_assayed":unknown,
            "functional_proportion_lower":F(positive,n),
            "functional_proportion_upper":F(positive+unknown,n),
        }
    w,p=groups["UNPIGMENTED"],groups["PIGMENTED"]
    low=w["functional_proportion_lower"]-p["functional_proportion_upper"]
    high=w["functional_proportion_upper"]-p["functional_proportion_lower"]
    width=high-low
    assert width==F(w["not_functionally_assayed"],N_WHITE)+F(p["not_functionally_assayed"],N_PIG)
    decision=("UNPIGMENTED_HIGHER" if low>0 else
              "PIGMENTED_HIGHER" if high<0 else
              "SIGNED_DIFFERENCE_NOT_IDENTIFIED")
    return {
        "groups":groups,
        "difference_unpigmented_minus_pigmented_lower":low,
        "difference_unpigmented_minus_pigmented_upper":high,
        "difference_identification_width":width,
        "strict_direction_identified":decision,
        "case_is_hypothetical_or_unassayed_not_empirical_validation":True,
    }


def metadata_strata(panel):
    """Nonempty (flower pigment × source annotation) cells, NOT genotype states."""
    result=[]
    for state in PHENOTYPES:
        for mask in STATUS:
            ids=sorted(x["taxon"] for x in panel if x["flower_pigment_category"]==state and
                       x["source_DFR_metadata"]==mask)
            if ids:
                result.append({
                    "flower_pigment_category":state,
                    "new_KT_source_metadata":mask,
                    "original_taxa":len(ids),
                    "source_only_first_taxon_for_possible_protocol_planning":ids[0],
                    "functional_assays_currently_qualified":0,
                })
    if len(result)!=5 or sum(x["original_taxa"] for x in result)!=22:
        raise ValueError("Source metadata strata unexpectedly changed")
    return result


def hypothetical_annotations_all_functional(panel):
    """Fiction: each of 12 CDS-annotated source taxa gets a positive gold-standard assay.

    This is intentionally NOT the assertion that CDS annotation proves
    function. It helps distinguish experimental observation from metadata
    while estimating how many additional independent full assays would
    be necessary to resolve one finite-panel rate contrast.
    """
    return {x["taxon"]:True for x in panel
            if x["source_DFR_metadata"]==STATUS[0]}


def hypothetical_certification_worlds(panel):
    """Two incompatible purely fictional futures starting from same 12 positives."""
    base=hypothetical_annotations_all_functional(panel)
    unknown_w=sorted(x["taxon"] for x in panel if x["flower_pigment_category"]=="UNPIGMENTED" and
                     x["taxon"] not in base)
    unknown_p=sorted(x["taxon"] for x in panel if x["flower_pigment_category"]=="PIGMENTED" and
                     x["taxon"] not in base)
    if len(unknown_w)!=4 or len(unknown_p)!=6:
        raise ValueError("2016 non-annotated-source denominator changed")
    # Positive difference: all four unknown white taxa function, and
    # one initially unknown pigmented taxon DOES NOT function.
    positive={**base,**{x:True for x in unknown_w},unknown_p[0]:False}
    # Negative difference: three unknown white taxa do not function;
    # their one still-unmeasured white taxon cannot reverse the sign.
    negative={**base,**{x:False for x in unknown_w[:3]}}
    return {"assume_12_annotated_have_confirmed_function":base,
            "additional_5_best_case_unpigmented_higher":positive,
            "additional_3_best_case_pigmented_higher":negative}


def minimal_best_case_certification(panel, start):
    """Exact oracle enumeration over all admissible extra +/− counts.

    This is BEST-CASE attainable evidence after seeing idealized assay
    outcomes, not expected power, worst-case guaranteed certification,
    actual experimental observation or optimal adaptive sampling.
    """
    summary=sharp_bounds(panel,start)
    g=summary["groups"]
    remaining_w=g["UNPIGMENTED"]["not_functionally_assayed"]
    remaining_p=g["PIGMENTED"]["not_functionally_assayed"]
    existing_pw=g["UNPIGMENTED"]["hypothetically_assayed_positive"]
    existing_nw=g["UNPIGMENTED"]["hypothetically_assayed_negative"]
    existing_pp=g["PIGMENTED"]["hypothetically_assayed_positive"]
    existing_np=g["PIGMENTED"]["hypothetically_assayed_negative"]
    out={"UNPIGMENTED_HIGHER":None,"PIGMENTED_HIGHER":None}
    for wp in range(remaining_w+1):
        for wn in range(remaining_w-wp+1):
            for pp in range(remaining_p+1):
                for pn in range(remaining_p-pp+1):
                    wplus,wminus=existing_pw+wp,existing_nw+wn
                    pplus,pminus=existing_pp+pp,existing_np+pn
                    remainw=N_WHITE-wplus-wminus
                    remainp=N_PIG-pplus-pminus
                    lower=F(wplus,N_WHITE)-F(pplus+remainp,N_PIG)
                    upper=F(wplus+remainw,N_WHITE)-F(pplus,N_PIG)
                    classification=("UNPIGMENTED_HIGHER" if lower>0 else
                                    "PIGMENTED_HIGHER" if upper<0 else None)
                    if classification:
                        cost=wp+wn+pp+pn
                        candidate=(cost,wp,wn,pp,pn)
                        if out[classification] is None or candidate<out[classification]:
                            out[classification]=candidate
    def encode(c):
        if c is None:
            return {"achievable":False,"minimum_additional_assays":None}
        cost,wp,wn,pp,pn=c
        return {"achievable":True,"minimum_additional_assays":cost,
                "illustrative_additional_positive_unpigmented":wp,
                "illustrative_additional_negative_unpigmented":wn,
                "illustrative_additional_positive_pigmented":pp,
                "illustrative_additional_negative_pigmented":pn}
    return {k:encode(v) for k,v in out.items()}


def minimax_width_budget(budget):
    """Exact outcome-independent interval-width minimization for a fixed panel.

    With no existing assays, each *perfectly valid* new white assay shrinks
    the sharp contrast interval by 1/7, and each pigmented assay by 1/15,
    REGARDLESS of the measured outcome. Therefore white-first strictly
    minimizes width at any fixed count, but does not ensure the SIGN is
    identifiable, does not establish phylogenetic portability and is not
    an ecological priority ranking for field sample collection.
    """
    if type(budget) is not int or not 0<=budget<=22:
        raise ValueError("Budget must be integer 0..22")
    w=min(budget,N_WHITE)
    p=budget-w
    return {
        "assay_budget_assuming_one_qualified_binary_taxon_assay_per_taxon":budget,
        "width_only_optimal_unpigmented_taxa":w,
        "width_only_optimal_pigmented_taxa":p,
        "possible_sharp_difference_interval_width":F(N_WHITE-w,N_WHITE)+F(N_PIG-p,N_PIG),
        "not_guaranteed_sign_certification_or_phylogenetic_predictive_validation":True,
    }


def report(join):
    panel=panel_from_original_source(join)
    scenarios=hypothetical_certification_worlds(panel)
    summaries={k:sharp_bounds(panel,v) for k,v in scenarios.items()}
    empty=sharp_bounds(panel)
    budget=[minimax_width_budget(b) for b in (0,1,5,7,8,12,17,22)]
    minimal_unassayed=minimal_best_case_certification(panel,{})
    minimal_after_12=minimal_best_case_certification(
        panel,scenarios["assume_12_annotated_have_confirmed_function"])
    if (empty["difference_unpigmented_minus_pigmented_lower"]!=-1 or
        empty["difference_unpigmented_minus_pigmented_upper"]!=1 or
        minimal_unassayed["UNPIGMENTED_HIGHER"]["minimum_additional_assays"]!=8 or
        minimal_unassayed["PIGMENTED_HIGHER"]["minimum_additional_assays"]!=8 or
        minimal_after_12["UNPIGMENTED_HIGHER"]["minimum_additional_assays"]!=5 or
        minimal_after_12["PIGMENTED_HIGHER"]["minimum_additional_assays"]!=3):
        raise ArithmeticError("Original 7×15 sharp-bound witness changed")
    return {
        "schema":"anto_v20_functional_flower_assay_sharp_bounds_prospective",
        "date":"2026-10-26",
        "source":"Ho & Smith 2016 original flower phenotype × KT source metadata (22 taxa)",
        "original_author_study_DOI":"10.1186/s12862-016-0675-3",
        "original_source_join":"results/iochrominae_2016_trait_KT_source_join_20261019.json",
        "all_source_taxon_functional_assays_currently_unmeasured":True,
        "ancestral_reexpression_rate_estimated":False,
        "validated_flower_toolkit_function_measured_by_ANTO":False,
        "independent_macroevolutionary_regain_clades_replicated_here":False,
        "same_taxon_samples_not_independent_loss_events":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "prospective_assay_definition_not_yet_implemented_in_field":list(PROTOCOL),
        "source_taxon_panel":panel,
        "source_annotation_strata":metadata_strata(panel),
        "observed_source_only_identification":empty,
        "source_annotated_CDS_fictions_not_observed":summaries,
        "unassayed_best_case_first_strict_sign_certificates":minimal_unassayed,
        "after_fictional_12_positives_best_case_additional_certificates":minimal_after_12,
        "width_only_budget_illustrations":budget,
        "no_worst_case_nonzero_sign_guarantee_even_after_all_22_if_true_rates_equal":True,
    }


def serializable(x):
    if isinstance(x,F):
        return str(x)
    if isinstance(x,dict):
        return {k:serializable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [serializable(v) for v in x]
    return x


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    p.add_argument("--expected",type=Path)
    args=p.parse_args()
    value=serializable(report(read(SOURCE_JOIN)))
    if args.expected and value!=json.loads(args.expected.read_text(encoding="utf-8")):
        raise RuntimeError("Source-locked v20 prospective assay identification receipt drift")
    out=json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(out,encoding="utf-8")
    else:
        print(out)


if __name__=="__main__":
    main()
