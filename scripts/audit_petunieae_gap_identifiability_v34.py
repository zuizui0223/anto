"""ANTO v34: can count/gap matched phylogenetic training controls exist?

Pure original-tree structural audit. NO pigment, gene-expression, predictor
fitting, random reassignments, hypothesis p-values, or new field data.
Source Newick must pass existing Wheeler-2023 original SHA256+MD5 gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from audit_petunieae_2023_processed_files_v26 import FILES, verify_original
from analyze_original_petunieae_phylogeny_v27 import (
    parse_newick, prune_outgroup, exact_monophyletic_partitions
)
from analyze_petunieae_loo_distance_v31 import tree_distances

CAPS=(8,15)
TOLERANCES={
    "abs_1e-5":("absolute",1e-5),
    "rel_05pct":("relative",0.05),
    "rel_10pct":("relative",0.10),
    "rel_20pct":("relative",0.20),
}
EXPECTED_INGROUP=59
EXPECTED_FOLD_SIGNATURE_SHA={
    8:"aedf4c5fc6790f6966421fa2e7b5cb642df4b96914b1dfdc3dc43ff78c2ea24a",
    15:"00c49ce4f95a703de674c84410719287b82d7d4caec5d43d583a326fc5a07020",
}
EXPECTED_V33_FULL_SHA="a61d5c7c7e0552a45356f0f966a2dbdbecd91e47a3980a06bc06d078ed312042"


def _real_nonnegative(x):
    if isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) or x<0:
        raise ValueError("INVALID_SOURCE_DISTANCE_OR_TOLERANCE")
    return float(x)


def conditional_training_set_count(boundary,farther,training_size):
    """Exact # T-subsets from eligible species containing >=1 boundary taxon."""
    if any(isinstance(x,bool) or not isinstance(x,int) or x<0
           for x in (boundary,farther,training_size)):
        raise ValueError("INVALID_CONDITIONAL_COMBINATORIC_ARGUMENT")
    from math import comb
    total=comb(boundary+farther,training_size) if training_size<=boundary+farther else 0
    all_far=comb(farther,training_size) if training_size<=farther else 0
    return total-all_far


def effective_tolerance(original_gap,mode):
    gap=_real_nonnegative(original_gap)
    if not isinstance(mode,(tuple,list)) or len(mode)!=2:
        raise ValueError("UNDECLARED_MATCHING_TOLERANCE")
    name,value=mode
    _real_nonnegative(value)
    if name=="absolute":
        return float(value)
    if name=="relative":
        if value>1:
            raise ValueError("UNREALISTIC_RELATIVE_TOLERANCE")
        return gap*float(value)
    raise ValueError("UNDECLARED_MATCHING_TOLERANCE")


def matching_structure(ids,target,heldout,pairwise,tolerance):
    """Counterfactual support for fixed target, train count and nearest-gap band.

    All closer candidates must be omitted; >=1 boundary candidate must be
    trained. Uses exact integer combinatorics, never a biological effect.
    """
    species=sorted(ids)
    if (not species or len(species)!=len(set(species)) or
            len(heldout)!=len(set(heldout)) or not heldout or
            target not in heldout or not set(heldout).issubset(species)):
        raise ValueError("TARGET_OR_FOLD_SOURCE_IDENTITY")
    tol=_real_nonnegative(tolerance)
    if tol==0:
        raise ValueError("NONPOSITIVE_MATCHING_TOLERANCE")
    train=sorted(set(species)-set(heldout))
    if not train:
        raise ValueError("EMPTY_ORIGINAL_TRAINING_GROUP")
    if target not in pairwise or target not in pairwise[target]:
        raise ValueError("MISSING_TARGET_DISTANCE_ROW")
    if not math.isclose(_real_nonnegative(pairwise[target][target]),0,abs_tol=1e-12):
        raise ValueError("DISTANCE_DIAGONAL_NOT_ZERO")
    for other in species:
        try:
            d=_real_nonnegative(pairwise[target][other])
            back=_real_nonnegative(pairwise[other][target])
        except KeyError as exc:
            raise ValueError("MISSING_ORIGINAL_DISTANCE_PAIR") from exc
        if not math.isclose(d,back,abs_tol=1e-9,rel_tol=1e-12):
            raise ValueError("NON_SYMMETRIC_SOURCE_DISTANCE_MATRIX")
    original_nearest=min(pairwise[target][j] for j in train)
    lower=max(0,original_nearest-tol)
    upper=original_nearest+tol
    closer=[]
    boundary=[]
    farther=[]
    for j in species:
        if j==target:
            continue
        d=pairwise[target][j]
        if d<lower:
            closer.append(j)
        elif d<=upper:
            boundary.append(j)
        else:
            farther.append(j)
    if (set(train)&set(closer) or not set(train)&set(boundary)
            or len(closer)+len(boundary)+len(farther)!=len(species)-1):
        raise ValueError("ORIGINAL_TRAINING_NOT_IN_PRESPECIFIED_GAP_BAND")
    train_size=len(train)
    admissible=conditional_training_set_count(len(boundary),len(farther),train_size)
    if admissible<1:
        raise ValueError("STRUCTURAL_SOURCE_MATCHING_NOT_ADMITTED")
    return {
        "original_training_species_count":train_size,
        "original_heldout_species_count":len(heldout),
        "original_nearest_training_gap":round(original_nearest,9),
        "allowed_nearest_gap_lower":round(lower,9),
        "allowed_nearest_gap_upper":round(upper,9),
        "n_forced_closer_exclusions":len(closer),
        "n_boundary_candidates":len(boundary),
        "n_farther_candidates":len(farther),
        "n_extra_exclusion_slots_after_forced_nearer_removal":len(heldout)-1-len(closer),
        "n_admissible_training_sets_exact":admissible,
        "n_nonidentical_training_sets_exact":admissible-1,
        "independent_counterfactual_set_exists":admissible>1,
        "all_training_sets_have_fixed_size_and_one_boundary_taxon":True,
    }


def _summary(values):
    if not values:
        return None
    v=sorted(values)
    return [round(v[0],8),round(v[len(v)//2],8),round(v[-1],8)]


def audit_folds(ids,folds,pairs):
    species=sorted(ids)
    if (not species or len(species)!=len(set(species)) or
            sorted(x for f in folds for x in f)!=species or any(not f for f in folds)):
        raise ValueError("FOLD_PARTITION_ORIGINAL_SOURCE_DRIFT")
    output={}
    per_target=[]
    for mode,setting in TOLERANCES.items():
        entries=[]
        for ci,fold in enumerate(folds):
            for target in sorted(fold):
                gap=min(pairs[target][other] for other in species if other not in fold)
                tol=effective_tolerance(gap,setting)
                row=matching_structure(species,target,fold,pairs,tol)
                row["target_species_source_code"]=target
                row["original_clade_number"]=ci
                entries.append(row)
        if len(entries)!=len(species) or len({e["target_species_source_code"] for e in entries})!=len(species):
            raise ValueError("ORIGINAL_SPECIES_COVERAGE_DRIFT")
        supported=[x for x in entries if x["independent_counterfactual_set_exists"]]
        alts=[e["n_nonidentical_training_sets_exact"] for e in entries]
        output[mode]={
            "n_targets":len(entries),
            "n_original_clades":len(folds),
            "n_targets_with_alternative":len(supported),
            "n_targets_without_alternative":len(entries)-len(supported),
            "n_only_original_training_set":sum(e["n_admissible_training_sets_exact"]==1 for e in entries),
            "n_not_admitted":sum(e["n_admissible_training_sets_exact"]<1 for e in entries),
            "fraction_targets_with_alternative":round(len(supported)/len(entries),8),
            "n_original_clades_with_any_alternative":len({e["original_clade_number"] for e in supported}),
            "forced_nearer_exclusions_min_median_max":_summary([e["n_forced_closer_exclusions"] for e in entries]),
            "remaining_exclusion_slots_min_median_max":_summary(
                [e["n_extra_exclusion_slots_after_forced_nearer_removal"] for e in entries]),
            "nearest_actual_gap_min_median_max":_summary(
                [e["original_nearest_training_gap"] for e in entries]),
            "n_alternative_sets_min_max_exact_as_strings":[str(min(alts)),str(max(alts))],
            "multi_clade_and_species_independence_not_established":True,
        }
        for e in entries:
            per_target.append({"mode":mode,"target":e["target_species_source_code"],
                "clade_index":e["original_clade_number"],
                "train_species":e["original_training_species_count"],
                "original_gap":e["original_nearest_training_gap"],
                "closer_forced":e["n_forced_closer_exclusions"],
                "boundary":e["n_boundary_candidates"],
                "other_than_original_admissible":str(e["n_nonidentical_training_sets_exact"])})
    return {
        "species_scored_once":len(species),
        "n_original_clades":len(folds),
        "clade_sizes_sorted":sorted(map(len,folds)),
        "tolerance_results":output,
        "primary_exact_gap_conditional_comparison_identifiable":
            output["abs_1e-5"]["n_targets_with_alternative"]>0,
        "per_target_structural_feasibility":per_target,
        "no_model_fit_or_pigment_expression_values_used":True,
    }


def _fingerprint(folds):
    ordered=sorted(sorted(f) for f in folds)
    return hashlib.sha256(
        json.dumps(ordered,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def report(tree_bytes,receipt):
    original=FILES["species_tree_pruned"]
    verify_original(tree_bytes,original)
    if (receipt.get("full_result_sha256")!=EXPECTED_V33_FULL_SHA
            or receipt.get("all_original_sources_verified") is not True
            or receipt.get("all_four_directional_composition_checks") is not True
            or len(receipt.get("comparisons",[]))!=4):
        raise ValueError("V33_SOURCE_OR_FROZEN_RECEIPT_DRIFT")
    tree=parse_newick(tree_bytes.decode("utf-8-sig"))
    ingroup=prune_outgroup(tree)
    ids=sorted(ingroup.tips())
    if (len(ids)!=EXPECTED_INGROUP or len(tree.tips())!=60
            or "BROW" not in tree.tips() or "BROW" in ids):
        raise ValueError("ORIGINAL_TREE_SOURCE_TIP_COUNT_DRIFT")
    distances=tree_distances(tree)
    caps={}
    for cap in CAPS:
        folds=exact_monophyletic_partitions(ingroup,cap)
        if _fingerprint(folds)!=EXPECTED_FOLD_SIGNATURE_SHA[cap]:
            raise ValueError("V33_ORIGINAL_MONOPHYLETIC_CLADE_MEMBERSHIP_DRIFT")
        expected=14 if cap==8 else 6
        if len(folds)!=expected:
            raise ValueError("V33_SOURCE_CLADE_NUMBER_DRIFT")
        cells=[c for c in receipt["comparisons"] if c["max_clade_size"]==cap]
        if len(cells)!=2 or any(c["total_clades"]!=expected for c in cells):
            raise ValueError("V33_PRIOR_COMPARISON_GATE_DRIFT")
        caps["max_clade_"+str(cap)]=audit_folds(ids,folds,distances)
    exact=[caps["max_clade_"+str(cap)]["tolerance_results"]["abs_1e-5"] for cap in CAPS]
    return {
        "schema":"anto_v34_original_tree_matched_training_count_and_nearest_gap_structural_support",
        "design":"POST_V33_EXPLORATORY_TOPOLOGICAL_POSITIVITY_AUDIT",
        "original_study_DOI":"10.1098/rspb.2023.0275",
        "original_tree_source_osf":"zg9cu",
        "original_tree_sha256":original["sha256"],
        "original_newick_md5_sha256_and_size_verified":True,
        "previous_v33_full_derived_result_sha256_authenticated_by_compact_receipt":EXPECTED_V33_FULL_SHA,
        "previous_v33_result_recomputed_in_this_run":False,
        "old_v33_original_clade_membership_sha256_matched":True,
        "original_source_species_n":EXPECTED_INGROUP,
        "pre_frozen_absolute_and_relative_matching_tolerances":TOLERANCES,
        "per_clade_cap":caps,
        "strict_identifiability_verdict":(
            "NO_SOURCE_TREE_TARGET_HAS_ALTERNATIVE_EXACT_GAP_COUNT_MATCH"
            if sum(e["n_targets_with_alternative"] for e in exact)==0 else
            "SOME_ORIGINAL_TARGETS_HAVE_STRUCTURALLY_ADMISSIBLE_ALTERNATIVES"),
        "no_numeric_predictive_model_was_refitted":True,
        "no_inference_about_gene_effect_or_ancestral_regain":True,
        "under_band_relaxation_distance_and_composition_still_change":True,
        "n_statistically_independent_empirical_studies_added":0,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--verified-original-tree",required=True,type=Path)
    parser.add_argument("--v33-receipt",type=Path,default=Path(
        "results/petunieae_target_size_control_v33_20261114.receipt.json"))
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--expected",type=Path)
    args=parser.parse_args()
    result=report(args.verified_original_tree.read_bytes(),json.loads(args.v33_receipt.read_text()))
    if args.expected:
        frozen=json.loads(args.expected.read_text())
        if frozen!=result:
            def first_diff(a,b,where="$"):
                if type(a)!=type(b):
                    return where + " type " + type(a).__name__ + "/" + type(b).__name__
                if isinstance(a,dict):
                    if set(a)!=set(b):
                        return where + " keys mismatch " + repr(sorted(set(a)^set(b))[:6])
                    for key in sorted(a):
                        if a[key]!=b[key]:
                            return first_diff(a[key],b[key],where+"."+str(key))
                elif isinstance(a,list):
                    if len(a)!=len(b):
                        return where+" list length "+str(len(a))+"/"+str(len(b))
                    for i,(x,y) in enumerate(zip(a,b)):
                        if x!=y:
                            return first_diff(x,y,where+"["+str(i)+"]")
                else:
                    return where+" frozen="+repr(a)+" computed="+repr(b)
                return where+" unknown drift"
            raise ValueError("V34_SOURCE_DERIVED_FROZEN_RESULT_DRIFT: "+first_diff(frozen,result))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    small={}
    for cap,info in result["per_clade_cap"].items():
        small[cap]={k:{
            "alt_targets":v["n_targets_with_alternative"],
            "alt_clades":v["n_original_clades_with_any_alternative"],
            "unique_targets":v["n_only_original_training_set"],
            "max_alternatives":v["n_alternative_sets_min_max_exact_as_strings"][1],
            "free_exclusion_slots":v["remaining_exclusion_slots_min_median_max"],
        } for k,v in info["tolerance_results"].items()}
    print("ANTO_V34_SOURCE_AUTHENTICATED_STRUCTURAL_RESULT",
          json.dumps({"verdict":result["strict_identifiability_verdict"],
                      "caps":small,"not_biological_effect":True},sort_keys=True))
    packed=json.dumps(result,sort_keys=True,separators=(",",":"))
    chunks=(len(packed)+5999)//6000
    for i in range(chunks):
        print("ANTO_V34_DERIVED_CHUNK",i+1,chunks,packed[i*6000:(i+1)*6000])


if __name__=="__main__":
    main()
