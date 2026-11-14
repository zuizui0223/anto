"""ANTO v29 exploratory fine clade discordance and size-matched relabeling controls.

Requires source-authenticated v26/v27/v28 results to replay EXACTLY first.
No new chemical observations; no numeric pooling across tree definitions.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from collections import defaultdict
from pathlib import Path

import analyze_original_petunieae_alternative_trees_v28 as v28
from analyze_original_petunieae_2023_species_prediction_v26 import (
    MODELS, source_real_data
)
from analyze_original_petunieae_phylogeny_v27 import (
    TreeGateError, evaluate_folds, exact_monophyletic_partitions,
    parse_newick
)
from audit_petunieae_2023_processed_files_v26 import (
    FILES, get_official_file, registered_location, verify_original
)

FINE_CUTS=(2,3,4,5,6)
NULL_CUTS=(3,5)
NULL_REPS=99
NULL_SEED_BASE=20261114
ENDPOINTS=("anthocyanidins_log1p_total","flavonols_log1p_total")
TREES=(v28.ORIGINAL_TREE,v28.ASTRAL_TREE,v28.IQ_FIRST_TREE)
MODELS_PLUS=tuple(k for k in MODELS if k!="train_mean_baseline")
EXPECTED_N=31
EXPECTED_PAIRWISE=31*30//2


def topology_folds(tree, taxa, cap, anchor=v28.ROOT_PRIMARY):
    rooted=v28.reroot_outgroup_edge(tree,anchor)
    sub=v28.prune_to_taxa(rooted,set(taxa))
    if sub is None or sorted(sub.tips())!=sorted(taxa):
        raise TreeGateError("V29_FINE_ROOT_PRUNE_EXACT_SPECIES_DRIFT")
    folds=exact_monophyletic_partitions(sub,cap)
    if sorted(t for f in folds for t in f)!=sorted(taxa):
        raise TreeGateError("V29_FINE_CLADE_PARTITION_COVERAGE_DRIFT")
    return folds


def signature(folds):
    """Exact unlabeled partition, irrespective of tree left-right order."""
    cells=[tuple(sorted(b)) for b in folds]
    if any(not cell for cell in cells):
        raise TreeGateError("EMPTY_CLADE")
    if sum(len(c) for c in cells)!=len(set(t for c in cells for t in c)):
        raise TreeGateError("DUPLICATE_TAXON_IN_PARTITION")
    return tuple(sorted(cells))


def coheldout_disagreements(folds_a, folds_b, taxa):
    """Number of species pairs treated jointly in only one of two foldings."""
    a=signature(folds_a)
    b=signature(folds_b)
    aa={taxon:i for i,cell in enumerate(a) for taxon in cell}
    bb={taxon:i for i,cell in enumerate(b) for taxon in cell}
    if set(aa)!=set(taxa) or set(bb)!=set(taxa):
        raise TreeGateError("COHELDOUT_PAIR_SPECIES_DRIFT")
    return sum((aa[x]==aa[y]) != (bb[x]==bb[y])
               for x,y in itertools.combinations(sorted(taxa),2))


def pct_gain(v,endpoint,model):
    row=v["endpoints"][endpoint]["models"][model]
    return 100*(1-row["RMSE_relative_to_same_fold_training_mean"])


def quantile(xs,q):
    if not xs or not 0<=q<=1 or any(not math.isfinite(x) for x in xs):
        raise TreeGateError("BAD_FIXED_NULL_QUANTILE")
    z=sorted(xs)
    pos=(len(z)-1)*q
    a=math.floor(pos)
    b=math.ceil(pos)
    return z[a]*(b-pos)+z[b]*(pos-a) if a!=b else z[a]


def same_size_label_reference(original_folds, species, seed, reps):
    """Uniform random assignments to fixed original cut group slots.

    THIS IS AN OBSERVATIONAL LABEL CONTROL, NOT A VALID PHYLOGENETIC PERMUTATION P.
    """
    if len(species)!=EXPECTED_N or reps!=NULL_REPS:
        raise TreeGateError("NULL_REFERENCE_DENOMINATOR_OR_REPLICATE_DRIFT")
    ids=sorted(species)
    sizes=[len(b) for b in original_folds]
    if sum(sizes)!=EXPECTED_N:
        raise TreeGateError("NULL_REFERENCE_FOLD_SIZE_DRIFT")
    rng=random.Random(seed)
    for _ in range(reps):
        shuffled=ids.copy()
        rng.shuffle(shuffled)
        cursor=0
        folds=[]
        for count in sizes:
            folds.append(sorted(shuffled[cursor:cursor+count]))
            cursor+=count
        if cursor!=EXPECTED_N or sorted(t for f in folds for t in f)!=ids:
            raise TreeGateError("NULL_REFERENCE_MISSING_SPECIES")
        yield folds


def size_reference_summary(samples,source_folds,observed_result,cap):
    scores={e:{m:[] for m in MODELS_PLUS} for e in ENDPOINTS}
    for i,folds in enumerate(same_size_label_reference(
            source_folds, sorted(samples),NULL_SEED_BASE+cap,NULL_REPS)):
        res=evaluate_folds(samples,folds,"size_random_control_"+str(cap)+"_"+str(i))
        for endpoint in ENDPOINTS:
            for model in MODELS_PLUS:
                scores[endpoint][model].append(round(pct_gain(res,endpoint,model),8))
    return {
        "fixed_seed":NULL_SEED_BASE+cap,
        "same_fold_size_multiset_sorted":sorted(map(len,source_folds)),
        "draws_count":NULL_REPS,
        "all_species_and_outcomes_preserved":True,
        "responses_never_used_to_select_random_folds":True,
        "no_valid_exchangeability_or_permutation_p_value":True,
        "models":{
            e:{
                m:{
                    "observed_original_tree_improvement_pct":round(
                        pct_gain(observed_result,e,m),8),
                    "null_all_99_improvements_pct":xs,
                    "null_median_pct":round(quantile(xs,.5),8),
                    "null_q025_pct":round(quantile(xs,.025),8),
                    "null_q975_pct":round(quantile(xs,.975),8),
                    "n_random_partitions_at_least_observed":sum(
                        x>=pct_gain(observed_result,e,m) for x in xs),
                    "denominator":NULL_REPS,
                } for m,xs in mod.items()
            } for e,mod in scores.items()
        },
    }


def fine_gate(fine,structural):
    has_discordance=any(v["n_distinct_holdout_partitions"]>1 for v in structural.values())
    fails=[]
    unique_checked=0
    for cap in FINE_CUTS:
        seen=set()
        for name in TREES:
            row=fine[name]["max_clade_"+str(cap)]
            sig=signature([f["source_tip_ids"] for f in row["fold_membership"]])
            if sig in seen:
                continue
            seen.add(sig)
            unique_checked+=1
            for endpoint in ENDPOINTS:
                val=pct_gain(row,endpoint,"branch_three")
                if val<=0:
                    fails.append([cap,name,endpoint,round(val,6)])
    verdict=("NO_EVALUABLE_PUBLISHED_FOLD_DISCORDANCE" if not has_discordance
             else "MIXED_DIRECTION_WHEN_CLADES_DISCORD" if fails
             else "EXPLORATORY_ALL_DISTINCT_FINE_PARTITIONS_POSITIVE")
    return {
        "verdict":verdict,
        "unique_real_fold_designs_across_caps":unique_checked,
        "some_actual_published_tree_partitions_distinct":has_discordance,
        "branch_three_negative_distinct_designs":fails,
        "not_an_independent_preregistration_or_field_replication":True,
    }


def analyze(*,expected_v26,expected_v27,expected_v28,
            fetch_original=get_official_file,fetch_2024=v28.fetch_published):
    cache_o={}
    cache_p={}
    def old(url):
        if url not in cache_o:
            cache_o[url]=fetch_original(url)
        return cache_o[url]
    def recent(url):
        if url not in cache_p:
            cache_p[url]=fetch_2024(url)
        return cache_p[url]
    frozen=json.loads(Path(expected_v28).read_text(encoding="utf-8"))
    recreated=v28.report(fetch_original=old,fetch_2024=recent,
                          expected_v26=expected_v26,expected_v27=expected_v27)
    if frozen!=recreated:
        raise TreeGateError("V29_UPSTREAM_V28_ORIGINAL_SOURCE_DRIFT")
    samples,_=source_real_data(fetch=old)
    source_tree_bytes=old(registered_location(FILES["species_tree_pruned"]))
    verify_original(source_tree_bytes,FILES["species_tree_pruned"])
    source_tree=parse_newick(source_tree_bytes.decode("utf-8-sig"))
    psource=v28.load_published(fetch=recent)
    astral=parse_newick(psource["astral_species_tree"])
    iqtree=parse_newick(psource["iqtree_individual_consensus"])
    mapping,outgroups=v28.sample_map(psource["published_speciesmap"])
    v28.published_tree_admission(astral,iqtree,mapping,outgroups)
    trees,overlap=v28.panels_from_sources(samples,source_tree,astral,iqtree,
                                          mapping,outgroups)
    names=overlap["primary_common_focal_codes_sorted"]
    if len(names)!=EXPECTED_N:
        raise TreeGateError("V29_ORIGINAL_SOURCE_SHARED_FOCAL_TAXA_DRIFT")
    shared={k:samples[k] for k in names}
    all_folds={}
    for cap in FINE_CUTS:
        all_folds[cap]={t:topology_folds(trees[t],names,cap) for t in TREES}
    structural={}
    for cap in FINE_CUTS:
        groups=all_folds[cap]
        unique={signature(f) for f in groups.values()}
        pairs={}
        for t1,t2 in itertools.combinations(TREES,2):
            count=coheldout_disagreements(groups[t1],groups[t2],names)
            pairs[t1+"_vs_"+t2]={
                "discordant_species_pairs":count,
                "of_total_species_pairs":EXPECTED_PAIRWISE,
                "exact_same_CV_partition":count==0,
            }
        structural["max_clade_"+str(cap)]={
            "n_distinct_holdout_partitions":len(unique),
            "n_evaluated_published_trees":len(TREES),
            "pairwise_cv_membership_disagreements":pairs,
        }
    # All topology checks are completed before reading any heldout performance.
    evaluated={t:{} for t in TREES}
    for cap in FINE_CUTS:
        for name in TREES:
            evaluated[name]["max_clade_"+str(cap)]=evaluate_folds(
                shared,all_folds[cap][name],name+"_max_clade_"+str(cap))
    controls={}
    for cap in NULL_CUTS:
        controls["max_clade_"+str(cap)]=size_reference_summary(
            shared,all_folds[cap][v28.ORIGINAL_TREE],
            evaluated[v28.ORIGINAL_TREE]["max_clade_"+str(cap)],cap)
    return {
        "schema":"anto_v29_post_v28_fine_published_phylogeny_and_label_reference",
        "design_status":"EXPLORATORY_AFTER_V28_WITH_FULL_NEGATIVE_REPORTING",
        "original_source_v26_v27_v28_exact_replay_passed":True,
        "same_2023_processed_RNA_HPLC_data_no_new_natural_dataset":True,
        "original_2024_ASTAL_IQtree_samplemap_hashes_verified":True,
        "source_DOIs":[v28.ORIGINAL_STUDY_DOI,v28.LATER_TREE_STUDY_DOI],
        "2024_external_source_commit":v28.LATER_PIN,
        "species_n":len(names),
        "source_taxa_exact_codes":names,
        "source_gene_models_no_retuning":list(MODELS),
        "pre_frozen_fine_maximum_clade_size_caps":list(FINE_CUTS),
        "tree_definitions":list(TREES),
        "rooting_original_OUT":"HUTE",
        "structural_folds_without_outcome":structural,
        "fine_CV_results":evaluated,
        "post_v28_exploratory_directional_gate":fine_gate(evaluated,structural),
        "size_matched_label_random_controls":controls,
        "same_size_random_reference_is_not_a_valid_phylogenetic_permutation_test":True,
        "same_2023_RNA_and_HPLC_not_2024_independent_pigment_validation":True,
        "does_not_prove_flower_colour_regain_or_ancestral_function_or_selection":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--expected-v26",type=Path,required=True)
    p.add_argument("--expected-v27",type=Path,required=True)
    p.add_argument("--expected-v28",type=Path,required=True)
    p.add_argument("--expected-v29",type=Path)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    r=analyze(expected_v26=args.expected_v26,
              expected_v27=args.expected_v27,expected_v28=args.expected_v28)
    if args.expected_v29 and json.loads(args.expected_v29.read_text(encoding="utf-8"))!=r:
        raise TreeGateError("V29_FROZEN_SOURCE_RESULT_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    z={}
    for cap in FINE_CUTS:
        z["max_clade_"+str(cap)]={
            "n_distinct_holdout_partitions":r["structural_folds_without_outcome"][
                "max_clade_"+str(cap)]["n_distinct_holdout_partitions"],
            "tree_directional_improvements":{
                tree:{
                    endpoint:round(pct_gain(r["fine_CV_results"][tree]["max_clade_"+str(cap)],
                                            endpoint,"branch_three"),4)
                    for endpoint in ENDPOINTS
                } for tree in TREES
            },
        }
    print("ANTO_V29_NEW_FINE_GENUINE_TOPOLOGY_RESULTS",json.dumps({
        "gate":r["post_v28_exploratory_directional_gate"],
        "fine":z,
        "random_label_controls":{
            key:{
                e: {
                    "observed_branch3_pct":val["branch_three"][
                        "observed_original_tree_improvement_pct"],
                    "median_random_pct":val["branch_three"]["null_median_pct"],
                    "q975_random_pct":val["branch_three"]["null_q975_pct"],
                    "null_at_least_observed":val["branch_three"][
                        "n_random_partitions_at_least_observed"],
                } for e,val in ctl["models"].items()
            } for key,ctl in r["size_matched_label_random_controls"].items()
        },
        "not_meta_analysis_or_ancestral_reacquisition":True,
    },sort_keys=True,separators=(",",":")))
    packed=json.dumps(r,sort_keys=True,separators=(",",":"))
    step=6000
    count=(len(packed)+step-1)//step
    for j in range(count):
        print("ANTO_V29_DERIVED_CHUNK",j+1,count,packed[j*step:(j+1)*step])


if __name__=="__main__":
    main()
