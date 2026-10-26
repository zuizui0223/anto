"""ANTO v21: bounded diagnostic misclassification and floral competence.

Extends v20's exact 7-unpigmented/15-pigmented SOURCE taxon frame. No new
functional enzyme, regulatory, re-expression, ancestry or fitness observation
is admitted. GenBank accession presence/CDS annotations must NEVER be treated
as binary functional-test outcomes. All positive/negative assays below are
explicit mathematical counterfactuals.

An explicit deterministic upper bound e_g on the NUMBER OF MISCLASSIFIED
taxon-level calls in each group makes sharp adversarial partial-identification
possible. Nominal sensitivity and specificity strictly between 0 and 1 DO
NOT imply any such deterministic finite-sample guarantee: all true binary
competence vectors have positive likelihood for any observed binary calls in
the usual full-support independent-assay likelihood.

The resulting finite-panel calculations are not confidence intervals,
probabilistic error models, phylogenetic tests, evolutionary regain rates,
real field sampling targets, or recommendations to assume independent taxa.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from anto_afterlife_empirical_admissibility_v19 import SOURCE_JOIN, read
from anto_functional_assay_identifiability_v20 import panel_from_original_source

ROOT=Path(__file__).resolve().parents[1]
FROZEN_V20=ROOT/"results/anto_v20_functional_assay_identification_20261026.json"
UNPIGMENTED=7
PIGMENTED=15


def _nonnegative_integer(x, label):
    if type(x) is not int or x < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return x


def false_call_bound(n:int, tested:int, reported_positive:int, max_miscalls:int)->dict:
    """Sharp true-positive range given <=e incorrect *tested* binary calls.

    Unknown/unmeasured taxa have arbitrary binary status. A positive call
    can be false, and a negative call can be false. No wrong-call pattern
    with more than e mismatches is admitted.
    """
    for name,v in (("n",n),("tested",tested),
                   ("reported_positive",reported_positive),
                   ("max_miscalls",max_miscalls)):
        _nonnegative_integer(v,name)
    if n < 1 or tested > n or reported_positive > tested:
        raise ValueError("Original group size, tested counts or positive calls invalid")
    if max_miscalls > tested:
        raise ValueError("Cannot cap miscalls beyond number of tested taxa")
    negative=tested-reported_positive
    unknown=n-tested
    minimum=reported_positive-min(max_miscalls,reported_positive)
    maximum=reported_positive+unknown+min(max_miscalls,negative)
    assert 0 <= minimum <= maximum <= n
    return {
        "source_taxa":n,
        "hypothetically_tested":tested,
        "reported_positive":reported_positive,
        "reported_negative":negative,
        "not_tested":unknown,
        "assumed_deterministic_miscall_cap":max_miscalls,
        "true_positive_minimum":minimum,
        "true_positive_maximum":maximum,
        "true_prevalence_lower":F(minimum,n),
        "true_prevalence_upper":F(maximum,n),
    }


def contrast_bounds(u_tested:int,u_positive:int,p_tested:int,p_positive:int,
                    e_unpigmented:int=0,e_pigmented:int=0)->dict:
    """Exact finite-panel contrast, 0-error reproduces ANTO v20.

    A group-level error allowance is an extra assumed deterministic
    validation result, not supplied by current sequencing/RNA records.
    """
    white=false_call_bound(UNPIGMENTED,u_tested,u_positive,e_unpigmented)
    colour=false_call_bound(PIGMENTED,p_tested,p_positive,e_pigmented)
    low=white["true_prevalence_lower"]-colour["true_prevalence_upper"]
    high=white["true_prevalence_upper"]-colour["true_prevalence_lower"]
    decision=("UNPIGMENTED_HIGHER" if low>0 else
              "PIGMENTED_HIGHER" if high<0 else
              "EXACT_EQUALITY" if low==high==0 else
              "STRICT_DIRECTION_UNIDENTIFIED")
    return {
        "hypothetical_binary_calls_only":True,
        "source_annotation_is_not_function":True,
        "unpigmented":white,
        "pigmented":colour,
        "unpigmented_minus_pigmented_lower":low,
        "unpigmented_minus_pigmented_upper":high,
        "finite_panel_direction":decision,
        "not_sampling_confidence_interval_or_phylogenetic_rate":True,
    }


def max_one_error_any_group(u_tested,u_positive,p_tested,p_positive):
    """Allow <=1 total miscalled taxon, without knowing which group.

    This is the union of 0- or 1-error possibilities across the 2 groups,
    not a cap of 1 error PER group. Ignore a location where no call exists.
    """
    allocations=[(0,0)]
    if u_tested>0:
        allocations.append((1,0))
    if p_tested>0:
        allocations.append((0,1))
    panels=[contrast_bounds(u_tested,u_positive,p_tested,p_positive,*caps)
            for caps in allocations]
    low=min(p["unpigmented_minus_pigmented_lower"] for p in panels)
    high=max(p["unpigmented_minus_pigmented_upper"] for p in panels)
    return {
        "at_most_one_total_miscall_hypothesis_only":True,
        "possible_allocation_caps":[list(z) for z in allocations],
        "difference_lower":low,
        "difference_upper":high,
        "positive_sign_certified":low>0,
        "negative_sign_certified":high<0,
    }


def best_case_minimum_per_group_caps(e_unpigmented:int,e_pigmented:int):
    """Minimal number of calls that CAN prove strict sign, if ideal outcomes.

    U-positive/P-negative is most favorable to U>P and its mirror to P>U.
    Completely prospective counterfactual *existence* witness, NOT power,
    sampling prescription, error-rate estimate or a guarantee.
    """
    _nonnegative_integer(e_unpigmented,"e_unpigmented")
    _nonnegative_integer(e_pigmented,"e_pigmented")
    if e_unpigmented>UNPIGMENTED or e_pigmented>PIGMENTED:
        raise ValueError("Error cap exceeds source group size")
    out={}
    for direction in ("UNPIGMENTED_HIGHER","PIGMENTED_HIGHER"):
        possible=[]
        for u in range(UNPIGMENTED+1):
            for p in range(PIGMENTED+1):
                if e_unpigmented>u or e_pigmented>p:
                    continue
                r=contrast_bounds(
                    u,u if direction=="UNPIGMENTED_HIGHER" else 0,
                    p,0 if direction=="UNPIGMENTED_HIGHER" else p,
                    e_unpigmented,e_pigmented,
                )
                if r["finite_panel_direction"]==direction:
                    possible.append((u+p,u,p))
        if possible:
            shortest=min(z[0] for z in possible)
            choices=sorted([list(z[1:]) for z in possible if z[0]==shortest])
            out[direction]={
                "possible_under_ideal_hypothetical_calls":True,
                "minimum_tested_taxa":shortest,
                "minimum_allocations_unpigmented_pigmented":choices,
            }
        else:
            out[direction]={
                "possible_under_ideal_hypothetical_calls":False,
                "minimum_tested_taxa":None,
                "minimum_allocations_unpigmented_pigmented":[],
            }
    return out


def stochastic_assay_support(sensitivity,specificity):
    """Nonzero random misclassification gives no logical sign guarantee.

    A Bernoulli assay with both sensitivity and specificity strictly in
    (0,1) gives positive likelihood to each observed result for each true
    state. Thus even a completely assayed finite panel cannot *logically*
    exclude any competence configuration from the observations alone.
    Bayesian posterior or confidence coverage may be useful, but require
    separately validated statistical assumptions and do not turn into
    the hard count cap used by contrast_bounds().
    """
    def rate(x,name):
        if type(x) is bool:
            raise ValueError(f"{name} must be a numeric probability")
        try:
            v=F(str(x))
        except (ValueError,TypeError,ZeroDivisionError) as e:
            raise ValueError(f"{name} must be rational") from e
        if not 0<v<1:
            raise ValueError(f"{name} must lie strictly inside (0,1)")
        return v
    se,sp=rate(sensitivity,"sensitivity"),rate(specificity,"specificity")
    return {
        "assumed_independent_full_support_binary_error_likelihood":True,
        "sensitivity":se,
        "specificity":sp,
        "all_22_functional_state_vectors_compatible_with_any_binary_call_vector":True,
        "full_panel_true_difference_logical_support":[F(-1),F(1)],
        "deterministic_strict_sign_can_be_proven_from_accuracy_rates_alone":False,
        "hard_miscall_cap_derived_from_sensitivity_specificity":False,
    }


def source_panel_no_assays():
    """Verify underlying frozen source; do not promote CDS to function."""
    panel=panel_from_original_source(read(SOURCE_JOIN))
    locked=read(FROZEN_V20)
    if panel!=locked["source_taxon_panel"]:
        raise ValueError("Original 22-taxon v20 source enrollment changed")
    if not locked["all_source_taxon_functional_assays_currently_unmeasured"]:
        raise ValueError("Previously unmeasured source was promoted incorrectly")
    if len(panel)!=22 or any(x["flower_functional_joint_toolkit_assay"] is not None
                             or x["assay_admitted"] is not False for x in panel):
        raise ValueError("v21 requires no admitted real functional outcome")
    return panel


def panel_receipt():
    panel=source_panel_no_assays()
    cap_cases={}
    for e_u,e_p in ((0,0),(0,1),(1,0),(1,1),(2,1)):
        cap_cases[f"U{e_u}_P{e_p}"]=best_case_minimum_per_group_caps(e_u,e_p)
    original=contrast_bounds(7,7,1,0)
    one_total=max_one_error_any_group(7,7,1,0)
    two_individual_caps=contrast_bounds(7,7,1,0,1,1)
    robust=contrast_bounds(7,7,4,0,1,1)
    if (original["unpigmented_minus_pigmented_lower"]!=F(1,15)
            or one_total["positive_sign_certified"]
            or two_individual_caps["finite_panel_direction"]!="STRICT_DIRECTION_UNIDENTIFIED"
            or robust["unpigmented_minus_pigmented_lower"]!=F(2,35)
            or cap_cases["U1_P1"]["UNPIGMENTED_HIGHER"]["minimum_tested_taxa"]!=11):
        raise ArithmeticError("v21 error-budget witness drift")
    return {
        "schema":"anto_v21_fallible_joint_flower_assay_identifiability",
        "date":"2026-10-26",
        "source":"Ho and Smith 2016 source 22-taxon pigment and CDS registry via ANTO v19/v20",
        "source_ledger_sha_not_a_functional_assay":"results/anto_v20_functional_assay_identification_20261026.json",
        "n_original_taxa":len(panel),
        "n_original_unpigmented":UNPIGMENTED,
        "n_original_pigmented":PIGMENTED,
        "source_taxa_with_new_joint_flower_function_assays":0,
        "native_ancestral_trait_regain_or_transgenerational_establishment_measured":False,
        "phylogenetically_independent_regain_events_added":0,
        "zero_misclassification_reproduces_v20":original,
        "v20_eight_taxon_certificate_after_up_to_one_total_miscall":one_total,
        "v20_eight_taxon_certificate_if_each_group_may_have_one_miscall":two_individual_caps,
        "ideal_eleven_taxon_certificate_if_each_group_may_have_one_miscall":robust,
        "ideal_counterfactual_minima_for_error_caps":cap_cases,
        "full_support_random_assay_example_not_a_validation_rate":stochastic_assay_support(
            F(19,20),F(19,20)),
        "assay_claim_boundary":(
            "A deterministic error cap requires validation independent of genomic/CDS "
            "annotation; a nominal 95% assay accuracy is NOT a maximum of one wrong "
            "taxon in 22, and no inferred natural functional positivity is reported."
        ),
        "beyond_current_scope":[
            "stage/tissue-matched same-organism coding and regulation assays",
            "intraspecific replicate validation and assay-calibration specimens",
            "ancestral polarity, introgression, donor complementation and plasticity",
            "lineage-specific re-expression and transgenerational stability",
            "trait-specific viable progeny and fitness",
            "new independent clades for holdout external biological validation",
        ],
        "Cirsium_organismal_empirical_owner":"zuizui0223/aza3",
    }


def serialize(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {k:serialize(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [serialize(x) for x in value]
    return value


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    parser.add_argument("--expected",type=Path)
    args=parser.parse_args()
    data=serialize(panel_receipt())
    if args.expected and data!=read(args.expected):
        raise ValueError("Source-locked v21 misclassification sensitivity receipt drift")
    payload=json.dumps(data,sort_keys=True,ensure_ascii=False,indent=2)+"\n"
    if args.output is None:
        print(payload)
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload,encoding="utf-8")


if __name__=="__main__":
    main()
