"""ANTO v30: nearest-training patristic gap and molecular novelty vs heldout pigment error.

Same *original* 2023 source RNA/HPLC+tree and exactly replayed v26-v29 receipts.
All inference here is descriptive repeated-fold evidence, not an independent
source/taxon experiment, a phylogenetic independent-significance test or a
causal molecular mechanism. No original author source rows are emitted.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import analyze_original_petunieae_fine_discordance_v29 as v29
import analyze_original_petunieae_alternative_trees_v28 as v28
from analyze_original_petunieae_2023_species_prediction_v26 import (
    MODELS, RIDGE_LAMBDA, source_real_data, ridge_predict
)
from analyze_original_petunieae_phylogeny_v27 import (
    TreeGateError, parse_newick, prune_outgroup, exact_monophyletic_partitions
)
from audit_petunieae_2023_processed_files_v26 import (
    FILES, registered_location, get_official_file, verify_original
)

FINE_AND_COARSE_31 = (2,3,4,5,6,8,15)
COARSE_59 = (8,15)
ENDPOINTS = ("anthocyanidins_log1p_total","flavonols_log1p_total")
BRANCH_MODEL = "branch_three"
NULL_MODEL = "train_mean_baseline"
BRANCH_GENES = ("F3primeH","F35H","FLS")
EXPECTED_ORIGINAL_SOURCE_SPECIES = 59
EXPECTED_SHARED_FOCAL_SPECIES = 31
PAT_UNIT = "published_2023_Newick_branch_length_units_not_calibrated_years"


def patristic_pairs(tree):
    """Read original undirected branch-weighted tree, BEFORE pruning or rerooting.

    A child-length is the parent's connecting edge length; the original root
    length, if any, is an unconnected stem. Missing child lengths are never
    interpreted as zero.
    """
    adjacent={}
    leaf_nodes={}
    number=0

    def descend(node,parent=None):
        nonlocal number
        index=number
        number+=1
        adjacent[index]=[]
        if parent is not None:
            length=node.length
            if length is None:
                raise TreeGateError("SOURCE_BRANCH_LENGTH_NOT_AVAILABLE")
            if not isinstance(length,(float,int)) or not math.isfinite(length) or length<0:
                raise TreeGateError("SOURCE_BRANCH_LENGTH_INVALID")
            adjacent[index].append((parent,float(length)))
            adjacent[parent].append((index,float(length)))
        if node.children:
            for child in node.children:
                descend(child,index)
        else:
            if not node.label or node.label in leaf_nodes:
                raise TreeGateError("SOURCE_TREE_TIP_LABEL_DUPLICATE_OR_EMPTY")
            leaf_nodes[node.label]=index
    descend(tree)
    if len(leaf_nodes)<20:
        raise TreeGateError("TOO_FEW_PUBLISHED_ORIGINAL_TIPS")
    pairs={}
    for source,src_idx in leaf_nodes.items():
        found={}
        def travel(n,previous,length):
            if n in found:
                raise TreeGateError("PATRISTIC_NON_TREE_CYCLE")
            found[n]=length
            for neighbor,edge in adjacent[n]:
                if neighbor!=previous:
                    travel(neighbor,n,length+edge)
        travel(src_idx,None,0.0)
        if len(found)!=number:
            raise TreeGateError("PATRISTIC_DISCONNECTED_SOURCE_TREE")
        pairs[source]={dest:found[idx] for dest,idx in leaf_nodes.items()}
    for a in leaf_nodes:
        if abs(pairs[a][a])>1e-12:
            raise TreeGateError("PATRISTIC_NONZERO_DIAGONAL")
        for b in leaf_nodes:
            if not math.isclose(pairs[a][b],pairs[b][a],rel_tol=1e-12,abs_tol=1e-12):
                raise TreeGateError("PATRISTIC_ASYMMETRIC")
    if not any(pairs[a][b]>0 for a,b in itertools.combinations(sorted(leaf_nodes),2)):
        raise TreeGateError("PATRISTIC_ALL_ZERO_SOURCE_TREE")
    return pairs


def training_only_nearest_expression(test,samples,train,genes=BRANCH_GENES):
    if tuple(genes)!=BRANCH_GENES or len(train)!=len(set(train)) or test in train:
        raise TreeGateError("MOLECULAR_MODULE_CHANGED")
    if len(train)<10:
        raise TreeGateError("MOLECULAR_TRAINING_SPECIES_TOO_FEW")
    x=[[samples[k]["genes"][gene] for gene in genes] for k in train]
    means=[math.fsum(v[j] for v in x)/len(x) for j in range(len(genes))]
    std=[math.sqrt(math.fsum((v[j]-means[j])**2 for v in x)/len(x))
         for j in range(len(genes))]
    scale=[s if s>1e-10 else 1.0 for s in std]
    tx=[(samples[test]["genes"][gene]-means[j])/scale[j]
        for j,gene in enumerate(genes)]
    distances=[math.sqrt(math.fsum(
        (tx[j]-(samples[t]["genes"][gene]-means[j])/scale[j])**2
        for j,gene in enumerate(genes))) for t in train]
    if any(not math.isfinite(z) or z<0 for z in distances):
        raise TreeGateError("MOLECULAR_DISTANCE_NONFINITE")
    return min(distances)


def pearson(x,y):
    if len(x)!=len(y) or len(x)<2:
        return None
    if any(not math.isfinite(a) for a in x+y):
        raise TreeGateError("NONFINITE_SOURCE_CORRELATION_INPUT")
    mx=math.fsum(x)/len(x); my=math.fsum(y)/len(y)
    xx=math.fsum((a-mx)**2 for a in x)
    yy=math.fsum((b-my)**2 for b in y)
    if xx<1e-20 or yy<1e-20:
        return None
    return math.fsum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(xx*yy)


def midrank(xs):
    if any(not math.isfinite(x) for x in xs):
        raise TreeGateError("NONFINITE_SOURCE_RANK")
    order=sorted(range(len(xs)),key=lambda i:(xs[i],i))
    ranks=[0.0]*len(xs)
    start=0
    while start<len(order):
        end=start+1
        while end<len(order) and xs[order[end]]==xs[order[start]]:
            end+=1
        val=(start+1+end)/2
        for j in range(start,end):
            ranks[order[j]]=val
        start=end
    return ranks


def spearman(xs,ys):
    return pearson(midrank(xs),midrank(ys))


def round_or_none(x,digits=7):
    return None if x is None else round(x,digits)


def within_block_corr(records, xkey, ykey):
    by_fold=defaultdict(list)
    for r in records:
        by_fold[r["fold_id"]].append(r)
    dx=[];dy=[]
    used_blocks=0
    for group in by_fold.values():
        if len(group)<2:
            continue
        used_blocks+=1
        mx=math.fsum(r[xkey] for r in group)/len(group)
        my=math.fsum(r[ykey] for r in group)/len(group)
        dx.extend(r[xkey]-mx for r in group)
        dy.extend(r[ykey]-my for r in group)
    return {
        "multi_species_blocks_used":used_blocks,
        "species_in_non_singleton_blocks_used":len(dx),
        "within_block_x_squared_deviation":round(math.fsum(x*x for x in dx),9),
        "demeaned_pearson":round_or_none(pearson(dx,dy)),
        "not_an_independent_cluster_adjusted_regression":True,
    }


def source_distances_in_folds(samples, folds, pairwise):
    ids=sorted(samples)
    if sorted(t for f in folds for t in f)!=ids:
        raise TreeGateError("PREDICTION_DENOMINATOR_DRIFT")
    if not set(ids).issubset(pairwise):
        raise TreeGateError("DISTANCE_SOURCE_TIPS_NOT_SAME_MODEL_SPECIES")
    records=[]
    for fold_id, group in enumerate(folds):
        training=sorted(set(ids)-set(group))
        if len(training)<10:
            raise TreeGateError("PREDICTION_TRAINING_TOO_SMALL")
        source_training={t:samples[t] for t in training}
        for species in group:
            phylo=min(pairwise[species][t] for t in training)
            nov=training_only_nearest_expression(species,samples,training)
            if not math.isfinite(phylo) or phylo<0:
                raise TreeGateError("PATRISTIC_NEAREST_SOURCE_INVALID")
            d={"species":species,"fold_id":fold_id,"block_n":len(group),
               "patristic":phylo,"expression_novelty":nov,"responses":{}}
            for endpoint in ENDPOINTS:
                observed=samples[species]["responses"][endpoint]
                errors={}
                for model,genes in MODELS.items():
                    forecast=ridge_predict(source_training,samples[species],endpoint,genes)
                    if not math.isfinite(forecast):
                        raise TreeGateError("MODEL_NONFINITE_SOURCE_FORECAST")
                    errors[model]=observed-forecast
                d["responses"][endpoint]={
                    "residuals":errors,
                    "sq_branch":errors[BRANCH_MODEL]**2,
                    "sq_baseline":errors[NULL_MODEL]**2,
                    "gain_sq":errors[NULL_MODEL]**2-errors[BRANCH_MODEL]**2,
                }
            records.append(d)
    if len(records)!=len(ids) or len({r["species"] for r in records})!=len(ids):
        raise TreeGateError("ORIGINAL_SPECIES_EXACTLY_ONCE_GATE")
    return sorted(records,key=lambda r:r["species"])


def check_existing_model_metrics(records,old):
    for e in ENDPOINTS:
        for model in MODELS:
            err=[r["responses"][e]["residuals"][model] for r in records]
            score=round(math.sqrt(math.fsum(x*x for x in err)/len(err)),8)
            expected=old["endpoints"][e]["models"][model][
                "RMSE_log1p_current_chemical_total"]
            if score != expected:
                raise TreeGateError("V30_UPSTREAM_CV_SAME_RMSE_METRIC_DRIFT_"+e+"_"+model)


def tier_assignments(records):
    ranks=sorted(records,key=lambda r:(r["patristic"],r["species"]))
    n=len(records)
    if n not in (EXPECTED_SHARED_FOCAL_SPECIES,EXPECTED_ORIGINAL_SOURCE_SPECIES):
        raise TreeGateError("DISTANCE_TIER_SPECIES_DENOMINATOR_DRIFT")
    tiers=([],[],[])
    for i,r in enumerate(ranks):
        tiers[min(2,3*i//n)].append(r)
    if sorted(len(t) for t in tiers) not in ([10,10,11],[19,20,20]):
        raise TreeGateError("DISTANCE_TERTILE_SIZE_DRIFT")
    boundaries=(tiers[0][-1]["patristic"],tiers[1][-1]["patristic"])
    crosses=sum(
        r["patristic"]==border
        for border in boundaries
        for r in records
    )  # tie count at boundary values, not an inference statistic.
    return tiers,crosses


def summarize_endpoint(records,endpoint):
    gaps=[r["patristic"] for r in records]
    novelty=[r["expression_novelty"] for r in records]
    values=[r["responses"][endpoint] for r in records]
    b=[r["sq_branch"] for r in values]
    a=[r["sq_baseline"] for r in values]
    g=[r["gain_sq"] for r in values]
    tiers,ties=tier_assignments(records)
    rmse_null=math.sqrt(math.fsum(a)/len(a))
    rmse_branch=math.sqrt(math.fsum(b)/len(b))
    if rmse_null<=0:
        raise TreeGateError("NULL_RMSE_ZERO")
    strata=[]
    for label,group in zip(("near","middle","far"),tiers):
        se_b=[r["responses"][endpoint]["sq_branch"] for r in group]
        se_a=[r["responses"][endpoint]["sq_baseline"] for r in group]
        rb=math.sqrt(math.fsum(se_b)/len(se_b))
        ra=math.sqrt(math.fsum(se_a)/len(se_a))
        strata.append({
            "source_distance_rank_bin":label,
            "n_original_test_species":len(group),
            "patristic_min":round(min(x["patristic"] for x in group),8),
            "patristic_max":round(max(x["patristic"] for x in group),8),
            "RMSE_null":round(ra,8),
            "RMSE_branch_three":round(rb,8),
            "RMSE_improvement_pct":round_or_none(
                100*(1-rb/ra) if ra>0 else None,5),
            "mean_branch_squared_error":round(math.fsum(se_b)/len(se_b),8),
            "mean_null_minus_branch_squared_error":round(
                math.fsum(aa-bb for aa,bb in zip(se_a,se_b))/len(group),8),
        })
    fold_groups=defaultdict(list)
    for r in records:
        fold_groups[r["fold_id"]].append(r)
    block_d=[];block_gain=[]
    block_pos=0
    for group in fold_groups.values():
        dg=math.fsum(r["patristic"] for r in group)/len(group)
        gg=math.fsum(r["responses"][endpoint]["gain_sq"] for r in group)/len(group)
        block_d.append(dg)
        block_gain.append(gg)
        block_pos+=gg>0
    with_rows=[
        {"fold_id":r["fold_id"],"patristic":r["patristic"],
         "sq_branch":r["responses"][endpoint]["sq_branch"],
         "sq_baseline":r["responses"][endpoint]["sq_baseline"],
         "gain_sq":r["responses"][endpoint]["gain_sq"]}
        for r in records
    ]
    return {
        "n_original_source_species_scored_once":len(records),
        "n_clade_holdout_blocks":len(fold_groups),
        "RMSE_null":round(rmse_null,8),
        "RMSE_branch_three":round(rmse_branch,8),
        "RMSE_improvement_pct":round(100*(1-rmse_branch/rmse_null),5),
        "rank_correlations_no_significance":{
            "phylo_gap_vs_branch_sq_error":round_or_none(spearman(gaps,b)),
            "phylo_gap_vs_train_mean_sq_error":round_or_none(spearman(gaps,a)),
            "phylo_gap_vs_null_minus_branch_sq_gain":round_or_none(spearman(gaps,g)),
            "expression_novelty_vs_branch_sq_error":round_or_none(spearman(novelty,b)),
            "expression_novelty_vs_null_minus_branch_sq_gain":round_or_none(spearman(novelty,g)),
            "phylo_gap_vs_expression_novelty":round_or_none(spearman(gaps,novelty)),
            "mean_fold_phylo_gap_vs_mean_fold_gain":round_or_none(
                spearman(block_d,block_gain)),
            "n_blocks_with_positive_mean_sq_gain":block_pos,
            "n_total_blocks":len(fold_groups),
        },
        "within_fold_demeaned":{
            "gap_vs_branch_sq_error":within_block_corr(with_rows,"patristic","sq_branch"),
            "gap_vs_null_minus_branch_sq_gain":within_block_corr(
                with_rows,"patristic","gain_sq"),
        },
        "distance_only_tertiles":strata,
        "distance_ties_at_rank_boundary_count":ties,
        "not_independent_species_folds_or_95pct_intervals":True,
    }


def summarize_fold(records, expected):
    check_existing_model_metrics(records,expected)
    gaps=[r["patristic"] for r in records]
    expressions=[r["expression_novelty"] for r in records]
    out={
        "n_source_species":len(records),
        "n_folds":len({r["fold_id"] for r in records}),
        "n_test_species_in_multi_species_folds":sum(
            r["block_n"]>1 for r in records),
        "distance_units":PAT_UNIT,
        "nearest_train_patristic_summary":{
            "min":round(min(gaps),8),
            "max":round(max(gaps),8),
            "median":round(sorted(gaps)[len(gaps)//2],8),
            "n_nonzero":sum(x>0 for x in gaps),
        },
        "training_only_expression_novelty":{
            "min":round(min(expressions),8),
            "max":round(max(expressions),8),
            "median":round(sorted(expressions)[len(expressions)//2],8),
        },
        "endpoints":{e:summarize_endpoint(records,e) for e in ENDPOINTS},
    }
    return out


def diagnostic(focal):
    out=[]
    for cap in (3,5):
        entry=focal["max_clade_"+str(cap)]
        for ep in ENDPOINTS:
            bins=entry["endpoints"][ep]["distance_only_tertiles"]
            near=bins[0]["RMSE_branch_three"]
            far=bins[2]["RMSE_branch_three"]
            corr=entry["endpoints"][ep]["rank_correlations_no_significance"][
                "phylo_gap_vs_branch_sq_error"]
            out.append({
                "max_clade_size":cap,
                "chemical_endpoint":ep,
                "far_branch_RMSE_gt_near_branch_RMSE":far>near,
                "phylo_distance_error_spearman_positive":corr is not None and corr>0,
                "near_branch_RMSE":near,
                "far_branch_RMSE":far,
                "spearman_gap_to_branch_sq_error":corr,
            })
    return {
        "hypotheses_declared_before_v30_original_data_readout":True,
        "checks":out,
        "all_four_far_RMSE_greater_near":all(
            a["far_branch_RMSE_gt_near_branch_RMSE"] for a in out),
        "all_four_Spearman_gap_error_positive":all(
            a["phylo_distance_error_spearman_positive"] for a in out),
        "not_a_causal_distance_effect_or_p_value":True,
    }


def analyze(*,expected_v26,expected_v27,expected_v28,expected_v29,
            fetch_original=get_official_file,fetch_published=v28.fetch_published):
    if not all((expected_v26,expected_v27,expected_v28,expected_v29)):
        raise TreeGateError("REQUIRED_FROZEN_SOURCE_RECEIPTS_MISSING")
    cache_o={};cache_p={}
    def source(url):
        if url not in cache_o:
            cache_o[url]=fetch_original(url)
        return cache_o[url]
    def published(url):
        if url not in cache_p:
            cache_p[url]=fetch_published(url)
        return cache_p[url]
    expected=json.loads(Path(expected_v29).read_text(encoding="utf-8"))
    prior=v29.analyze(
        expected_v26=expected_v26,
        expected_v27=expected_v27,
        expected_v28=expected_v28,
        fetch_original=source,
        fetch_2024=published)
    if prior!=expected:
        raise TreeGateError("V30_V29_ORIGINAL_SOURCE_REPLAY_OR_RESULT_DRIFT")
    samples,_=source_real_data(fetch=source)
    original_bytes=source(registered_location(FILES["species_tree_pruned"]))
    verify_original(original_bytes,FILES["species_tree_pruned"])
    tree=parse_newick(original_bytes.decode("utf-8-sig"))
    pairwise=patristic_pairs(tree)
    ids59=sorted(samples)
    if len(ids59)!=EXPECTED_ORIGINAL_SOURCE_SPECIES:
        raise TreeGateError("V30_2023_59_SOURCE_SPECIES_DRIFT")
    if set(tree.tips())!=set(ids59)|{"BROW"}:
        raise TreeGateError("V30_SOURCE_TREE_SPECIES_SET_DRIFT")
    shared=prior["source_taxa_exact_codes"]
    if len(shared)!=EXPECTED_SHARED_FOCAL_SPECIES or not set(shared).issubset(samples):
        raise TreeGateError("V30_SOURCE_31_FOCAL_PANEL_DRIFT")
    all31={k:samples[k] for k in shared}
    eval31={}
    for cap in FINE_AND_COARSE_31:
        folds=v29.topology_folds(tree,shared,cap)
        if cap in v29.FINE_CUTS:
            saved=prior["fine_CV_results"][v28.ORIGINAL_TREE]["max_clade_"+str(cap)]
        else:
            old28=json.loads(Path(expected_v28).read_text(encoding="utf-8"))
            saved=old28["primary_original_vs_two_published_alternative_definitions"][
                v28.ORIGINAL_TREE]["max_clade_"+str(cap)]
        records=source_distances_in_folds(all31,folds,pairwise)
        eval31["max_clade_"+str(cap)]=summarize_fold(records,saved)
    eval59={}
    source_ingroup=prune_outgroup(tree)
    v27=json.loads(Path(expected_v27).read_text(encoding="utf-8"))
    for cap in COARSE_59:
        folds=exact_monophyletic_partitions(source_ingroup,cap)
        records=source_distances_in_folds(samples,folds,pairwise)
        saved=v27["evaluation"]["max_clade_"+str(cap)]
        eval59["max_clade_"+str(cap)]=summarize_fold(records,saved)
    return {
        "schema":"anto_v30_original_2023_source_phylogenetic_gap_vs_heldout_pigment_errors",
        "design_status":"EXPLORATORY_POST_V29_WITH_PRE_OUTCOME_PUBLISHED_PROTOCOL",
        "source_2023_DOI":v28.ORIGINAL_STUDY_DOI,
        "external_2024_tree_DOI":v28.LATER_TREE_STUDY_DOI,
        "original_2023_tree_sha256":FILES["species_tree_pruned"]["sha256"],
        "original_2024_tree_source_commit":v28.LATER_PIN,
        "original_source_v26_to_v29_all_receipts_exactly_replayed":True,
        "source_only_real_original_RNA_HPLC_and_Newick":True,
        "original_source_tip_count":len(pairwise),
        "source_tree_path_distance_units":PAT_UNIT,
        "source_only_branch_lengths_not_true_calendar_divergence_time":True,
        "covariate_definitions":{
            "phylogenetic_gap":"minimum 2023 original complete-tree patristic path to any training species; 2023 original branch units",
            "expression_novelty":"nearest standardized F3primeH/F35H/FLS log1p gene-expression Euclidean training species with within-training scaling",
            "response":"current 2023 original HPLC species log1p anthocyanidin / flavonol totals",
            "prediction":"same six fixed v26 ridge modules and each fold train-only intercept and scaling",
        },
        "n_original_2023_ingroup_species":len(samples),
        "n_original_2023_and_2024_focal_overlap":len(shared),
        "panels":{
            "31_focal_original_tree_maxclades_2_3_4_5_6_8_15":eval31,
            "59_all_original_ingroup_maxclades_8_15":eval59,
        },
        "exploratory_near_far_error_prediction_gate_31_focal":diagnostic(eval31),
        "species_and_repeated_clade_folds_not_independent_genetic_events":True,
        "multiple_cutoff_spearman_and_tertiles_no_inferential_p_values":True,
        "clade_composition_and_chemical_scale_not_causally_separated":True,
        "not_visible_flower_colour_or_biochemical_functional_reexpression_or_fitness":True,
        "no_cross_study_HPLC_or_intergenerational_species_replication":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("26","27","28","29"):
        p.add_argument("--expected-v"+n,required=True,type=Path)
    p.add_argument("--expected-v30",type=Path)
    p.add_argument("--output",required=True,type=Path)
    x=p.parse_args()
    result=analyze(
        expected_v26=x.expected_v26,expected_v27=x.expected_v27,
        expected_v28=x.expected_v28,expected_v29=x.expected_v29)
    if x.expected_v30 and json.loads(x.expected_v30.read_text(encoding="utf-8"))!=result:
        raise TreeGateError("V30_ORIGINAL_PHYLO_DISTANCE_FROZEN_RECEIPT_DRIFT")
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    brief={}
    for panel,rows in result["panels"].items():
        brief[panel]={}
        for cut,obj in rows.items():
            brief[panel][cut]={
                "patristic_range":obj["nearest_train_patristic_summary"],
                "endpoints":{name:{
                    "branch_rmse_gain_pct":v["RMSE_improvement_pct"],
                    "rho_phylo_vs_branch_sq_error":v["rank_correlations_no_significance"][
                        "phylo_gap_vs_branch_sq_error"],
                    "rho_expression_novelty_vs_branch_sq_error":v["rank_correlations_no_significance"][
                        "expression_novelty_vs_branch_sq_error"],
                    "within_fold_corr_gap_vs_branch_sq_error":v["within_fold_demeaned"][
                        "gap_vs_branch_sq_error"]["demeaned_pearson"],
                    "near_mid_far_branch_RMSE":[z["RMSE_branch_three"]
                                               for z in v["distance_only_tertiles"]],
                    "near_mid_far_rmse_gain_pct":[z["RMSE_improvement_pct"]
                                                for z in v["distance_only_tertiles"]],
                } for name,v in obj["endpoints"].items()}
            }
    print("ANTO_V30_PATRISTIC_SOURCE_AUTHENTICATED_RESULT",json.dumps({
        "directional_diagnostic":result["exploratory_near_far_error_prediction_gate_31_focal"],
        "panels":brief,
        "original_branch_units_not_calendar":True,
        "not_causal_molecular_evolution":True,
    },sort_keys=True,separators=(",",":")))
    packed=json.dumps(result,sort_keys=True,separators=(",",":"),ensure_ascii=False)
    step=6000
    total=(len(packed)+step-1)//step
    for i in range(total):
        print("ANTO_V30_DERIVED_RECEIPT_CHUNK",i+1,total,packed[i*step:(i+1)*step])


if __name__=="__main__":
    main()
