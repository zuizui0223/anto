"""ANTO v28: 2023 pigment prediction under *real published* 2024 alternative trees.

Strict source gate: original Wheeler 2023 OSF CSV+tree v26/v27 AND three
author-Pezzi 2024 GitHub repository blobs pinned to an exact commit/blob hash.
No original RNA/HPLC rows, phylogenies or individual sequences are redistributed.
The newer tree study reuses related RNA resources: NOT independent pigment validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

from audit_petunieae_2023_processed_files_v26 import (
    FILES, registered_location, get_official_file, verify_original
)
from analyze_original_petunieae_2023_species_prediction_v26 import (
    MODELS, source_real_data
)
from analyze_original_petunieae_phylogeny_v27 import (
    Node, TreeGateError, parse_newick, prune_outgroup,
    exact_monophyletic_partitions, evaluate_folds, analyze_original_phylogeny
)

ORIGINAL_STUDY_DOI = "10.1098/rspb.2023.0275"
LATER_TREE_STUDY_DOI = "10.1016/j.ympev.2024.108136"
LATER_REPO = "pedrohpezzi/Petunia-Calibrachoa-Fabiana_TreeDiscordance"
LATER_PIN = "0b33d5e992c39e191393bc0b86671c5c6a378e1e"
ORIGINAL_TREE = "2023_11gene_original"
ASTRAL_TREE = "2024_ASTRAL_species_assigned"
IQ_FIRST_TREE = "2024_IQtree_first_sample"
ROOT_PRIMARY = "HUTE"
ROOT_SECONDARY = "PLNY"
MAX_CLADE_SIZES = (8, 15)
OUTGROUP_CODES = frozenset(("BOER", "BRCU", "HUTE", "LEAC", "NIRI", "PLNY"))
EXPECTED_PUBLISHED_SPECIES = 43
EXPECTED_PUBLISHED_INDIVIDUALS = 117
EXPECTED_SHARED_FOCAL = 31
EXPECTED_SHARED_ALL = 37

EXTERNAL = {
    "astral_species_tree": {
        "path": "Data/ASTRAL/SpeciesAssigned/ASTRAL_SpeciesAssigned_pp_scored.treefile",
        "bytes": 1850,
        "git_blob_sha1": "04b4333e7a2e505e87dd321df1f4be88562ea4c4",
    },
    "iqtree_individual_consensus": {
        "path": "Data/IQTree/SuperMatrix_PFCOut.partitions_iqtree.txt.contree",
        "bytes": 4271,
        "git_blob_sha1": "dbbe0ff0f352c6fc163f6c1a3117612bf0982984",
    },
    "published_speciesmap": {
        "path": "Data/HyDe/speciesmap.txt",
        "bytes": 1269,
        "git_blob_sha1": "30b670968ab1033cee6cb3eb5f74904529d9794b",
    },
}


def git_blob_sha(data):
    if not isinstance(data, bytes):
        raise TreeGateError("BLOB_BYTES_REQUIRED")
    head = b"blob " + str(len(data)).encode("ascii") + b"\x00"
    return hashlib.sha1(head + data).hexdigest()


def pinned_url(record):
    return ("https://raw.githubusercontent.com/" + LATER_REPO + "/" +
            LATER_PIN + "/" + record["path"])


def fetch_published(url, max_bytes=10000):
    prefix = "https://raw.githubusercontent.com/" + LATER_REPO + "/" + LATER_PIN + "/"
    if not isinstance(url, str) or not url.startswith(prefix) or "?" in url or ".." in url:
        raise TreeGateError("UNAUTHORIZED_2024_SOURCE_URL")
    req = urllib.request.Request(
        url, headers={"User-Agent": "ANTO-v28-original-published-phylogeny-audit/1.0"}
    )
    with urllib.request.urlopen(req, timeout=25) as resp:
        data = resp.read(max_bytes + 1)
    if not data or len(data) > max_bytes:
        raise TreeGateError("EXTERNAL_PUBLISHED_SOURCE_NOT_ADMITTED")
    return data


def load_published(fetch=fetch_published):
    result = {}
    for key, meta in EXTERNAL.items():
        raw = fetch(pinned_url(meta))
        if not isinstance(raw, bytes) or len(raw) != meta["bytes"]:
            raise TreeGateError("EXTERNAL_PUBLISHED_BYTE_LENGTH_MISMATCH:" + key)
        if git_blob_sha(raw) != meta["git_blob_sha1"]:
            raise TreeGateError("EXTERNAL_PUBLISHED_BLOB_SHA_MISMATCH:" + key)
        result[key] = raw.decode("utf-8-sig")
    return result


def sample_map(text):
    """Only actual author rows; never infer phenotype-based taxon synonyms."""
    mapping = {}
    outgroups = set()
    for line in text.splitlines():
        if not line.strip():
            continue
        cells = line.strip().split("\t")
        if len(cells) != 2:
            raise TreeGateError("SOURCE_SAMPLE_MAP_BAD_ROW")
        key, source_target = cells
        if not re.fullmatch(r"[A-Z]{4}[1-9]", key):
            raise TreeGateError("SOURCE_SAMPLE_MAP_ID_NOT_AUTHOR_FORMAT")
        if key in mapping:
            raise TreeGateError("SOURCE_SAMPLE_MAP_DUPLICATE")
        species = key[:4]
        if source_target == "OUT":
            outgroups.add(species)
        elif source_target != species:
            raise TreeGateError("SOURCE_SAMPLE_MAP_SPECIES_DISAGREEMENT")
        mapping[key] = species
    if (len(mapping) != EXPECTED_PUBLISHED_INDIVIDUALS or
            outgroups != OUTGROUP_CODES or
            len(set(mapping.values())) != EXPECTED_PUBLISHED_SPECIES):
        raise TreeGateError("SOURCE_SAMPLE_MAP_DENOMINATOR_OR_OUTGROUP_DRIFT")
    return mapping, outgroups


def published_tree_admission(astral, iqtree, mapping, outgroups):
    ext_species = set(astral.tips())
    if len(ext_species) != EXPECTED_PUBLISHED_SPECIES:
        raise TreeGateError("ASTRAL_PUBLISHED_SPECIES_DENOMINATOR_DRIFT")
    if ext_species != set(mapping.values()):
        raise TreeGateError("ASTRAL_SAMPLE_MAP_SPECIES_MISMATCH")
    if set(iqtree.tips()) != set(mapping):
        raise TreeGateError("IQTREE_INDIVIDUAL_SAMPLE_MAP_MISMATCH")
    if outgroups != OUTGROUP_CODES or not {ROOT_PRIMARY, ROOT_SECONDARY}.issubset(ext_species):
        raise TreeGateError("PUBLISHED_OUTGROUP_ANCHOR_DRIFT")
    if not set(t[:4] for t in iqtree.tips()) == ext_species:
        raise TreeGateError("IQTREE_ASTRAL_TIP_PREFIX_NOT_SOURCE_VALIDATED")
    return {
        "original_2024_source_tree_species":len(ext_species),
        "original_2024_IQtree_individual_tips":len(iqtree.tips()),
        "original_2024_author_sample_map_rows":len(mapping),
        "source_designated_OUT_codes":sorted(outgroups),
        "map_join_strict_all_samples":True,
    }


def prune_to_taxa(node, retained):
    """Prune original Newick without imputing missing branch information."""
    if not node.children:
        return Node(label=node.label) if node.label in retained else None
    children = tuple(z for child in node.children
                     if (z := prune_to_taxa(child, retained)) is not None)
    if not children:
        return None
    if len(children) == 1:
        return children[0]
    return Node(children=children)


def collapse_to_one_sample_per_species(tree, mapping, rank):
    """One predeclared original sequencing sample per original author species."""
    if rank not in ("first", "middle", "last"):
        raise TreeGateError("UNDECLARED_INDIVIDUAL_SELECTION")
    samples = defaultdict(list)
    for sample, species in mapping.items():
        samples[species].append(sample)
    chosen={}
    for species, seq in samples.items():
        ordered=sorted(seq)
        idx = (0 if rank == "first" else len(ordered)//2 if rank=="middle" else -1)
        chosen[ordered[idx]]=species
    def visit(node):
        if not node.children:
            if node.label not in chosen:
                return None
            return Node(label=chosen[node.label])
        children=tuple(c for ch in node.children if (c:=visit(ch)) is not None)
        if not children:
            return None
        if len(children)==1:
            return children[0]
        return Node(children=children)
    reduced=visit(tree)
    if (reduced is None or len(reduced.tips())!=EXPECTED_PUBLISHED_SPECIES or
            set(reduced.tips())!=set(mapping.values()) or
            len(chosen)!=EXPECTED_PUBLISHED_SPECIES):
        raise TreeGateError("IQTREE_SPECIES_REPRESENTATIVE_PRUNE_DRIFT")
    return reduced, chosen


def reroot_outgroup_edge(tree, anchor):
    """Reorient a tree as undirected edges and root at one declared tip edge.

    No branch lengths are used in rooted clade selection. Internal labels are
    lost; no bootstrap/posterior-supported rooting is inferred here.
    """
    if anchor not in tree.tips():
        raise TreeGateError("ROOT_OUTGROUP_NOT_PRESENT")
    adjacency = {}
    labels = {}
    serial=0
    def fill(node,parent=None):
        nonlocal serial
        idx=serial
        serial+=1
        adjacency[idx]=[]
        labels[idx]=node.label if not node.children else None
        if parent is not None:
            adjacency[idx].append(parent)
            adjacency[parent].append(idx)
        for child in node.children:
            fill(child,idx)
    fill(tree)
    hits=[n for n,label in labels.items() if label==anchor]
    if len(hits)!=1 or len(adjacency[hits[0]])!=1:
        raise TreeGateError("ROOT_ANCHOR_NOT_UNIQUE_TERMINAL")
    anchor_node=hits[0]
    origin=adjacency[anchor_node][0]
    def orient(idx,parent):
        descend=tuple(orient(z,idx) for z in adjacency[idx] if z!=parent)
        if not descend:
            if not labels[idx]:
                raise TreeGateError("ROOT_CREATED_INTERNAL_TIP")
            return Node(label=labels[idx])
        if len(descend)==1:
            return descend[0]
        return Node(children=descend)
    rerooted=Node(children=(Node(label=anchor),orient(origin,anchor_node)))
    if set(rerooted.tips())!=set(tree.tips()) or len(rerooted.tips())!=len(tree.tips()):
        raise TreeGateError("REROOT_CHANGED_SOURCE_TAXON_IDS")
    return rerooted


def panels_from_sources(samples, source_2023, astral, iqtree, mapping, outgroups):
    measured=set(samples)
    external=set(astral.tips())
    shared_all=measured&external
    shared_focal=shared_all-outgroups
    if (len(measured)!=59 or len(shared_focal)!=EXPECTED_SHARED_FOCAL or
            len(shared_all)!=EXPECTED_SHARED_ALL or
            shared_all-shared_focal!=outgroups):
        raise TreeGateError("REGISTERED_31_FOCAL_37_SHARED_PANELS_DRIFT")
    if set(source_2023.tips()) not in (measured, measured|{"BROW"}):
        raise TreeGateError("SOURCE_2023_TREE_SPECIES_DRIFT")
    source_tree = prune_outgroup(source_2023)
    if source_tree is None or set(source_tree.tips())!=measured:
        raise TreeGateError("SOURCE_2023_INGROUP_TREE_DRIFT")
    iqt={}
    for rank in ("first","middle","last"):
        iqt[rank],_ = collapse_to_one_sample_per_species(iqtree,mapping,rank)
    trees={
        ORIGINAL_TREE:source_tree,
        ASTRAL_TREE:astral,
        IQ_FIRST_TREE:iqt["first"],
        "2024_IQtree_middle_sample":iqt["middle"],
        "2024_IQtree_last_sample":iqt["last"],
    }
    for name,t in trees.items():
        if not shared_all.issubset(set(t.tips())):
            raise TreeGateError("PUBLISHED_TREE_MISSING_ALL_SHARED_TAXA:" + name)
    return trees, {
        "primary_common_focal_taxa_count":len(shared_focal),
        "primary_common_focal_codes_sorted":sorted(shared_focal),
        "secondary_shared_including_OUT_count":len(shared_all),
        "secondary_common_codes_sorted":sorted(shared_all),
        "2024_OUT_codes_matched_to_2023_chemistry":sorted(outgroups),
        "no_source_taxon_phenotype_based_matching":True,
    }


def evaluate_tree_panel(samples, tree, retained, anchor, label):
    rooted=reroot_outgroup_edge(tree,anchor)
    subset=prune_to_taxa(rooted,set(retained))
    if subset is None or set(subset.tips())!=set(retained):
        raise TreeGateError("ROOTED_TREE_PANEL_EXACT_COVERAGE_DRIFT")
    sub={k:samples[k] for k in retained}
    out={}
    for limit in MAX_CLADE_SIZES:
        folds=exact_monophyletic_partitions(subset,limit)
        result=evaluate_folds(sub,folds,label+"__maxclade"+str(limit))
        if result["n_species_tested_once"] != len(retained):
            raise TreeGateError("PUBLISHED_TREE_CV_DENOMINATOR_DRIFT")
        out["max_clade_"+str(limit)]=result
    return out


def positive_primary_gate(results):
    failures=[]
    for name in (ORIGINAL_TREE,ASTRAL_TREE,IQ_FIRST_TREE):
        for cap in MAX_CLADE_SIZES:
            dat=results[name]["max_clade_"+str(cap)]
            if not (dat["n_multispecies_clade_blocks"]>=3 and
                    dat["n_species_in_multispecies_clade_blocks"]>=math.ceil(
                        dat["n_species_tested_once"]/2)):
                failures.append([name,cap,"insufficient_multispecies_partition"])
            for endpoint in ("anthocyanidins_log1p_total","flavonols_log1p_total"):
                rel = dat["endpoints"][endpoint]["models"]["branch_three"][
                    "RMSE_relative_to_same_fold_training_mean"]
                if not rel<1:
                    failures.append([name,cap,endpoint,"branch_three_not_better_than_null"])
    return {
        "verdict":("DESCRIPTIVE_PORTABILITY_ACROSS_THREE_PUBLISHED_TREE_DEFINITIONS"
                   if not failures
                   else "NOT_UNIFORMLY_PORTABLE_ACROSS_PUBLISHED_TREE_DEFINITIONS"),
        "predeclared_total_direction_checks":12,
        "failed_predeclared_gates":failures,
        "does_not_establish_phylogenetic_independence_or_empirical_trait_regain":True,
    }


def report(*, fetch_original=get_official_file, fetch_2024=fetch_published,
           expected_v26, expected_v27):
    if not expected_v26 or not expected_v27:
        raise TreeGateError("FROZEN_SOURCE_BASELINE_REQUIRED")
    cache={}
    def cache_original(url):
        if url not in cache:
            cache[url]=fetch_original(url)
        return cache[url]
    # Exact source and old-result lock before the new 31-species contrast.
    prior=analyze_original_phylogeny(fetch=cache_original,expected_v26=expected_v26)
    if prior != json.loads(Path(expected_v27).read_text(encoding="utf-8")):
        raise TreeGateError("V27_PREVIOUS_FROZEN_REAL_SOURCE_OR_OUTCOME_DRIFT")
    samples,reconciled=source_real_data(fetch=cache_original)
    orig_bytes=cache_original(registered_location(FILES["species_tree_pruned"]))
    verify_original(orig_bytes,FILES["species_tree_pruned"])
    source_2023=parse_newick(orig_bytes.decode("utf-8-sig"))
    source_2024=load_published(fetch=fetch_2024)
    astral=parse_newick(source_2024["astral_species_tree"])
    iqtree=parse_newick(source_2024["iqtree_individual_consensus"])
    mapping,outgroups=sample_map(source_2024["published_speciesmap"])
    external_admission=published_tree_admission(astral,iqtree,mapping,outgroups)
    trees, overlap=panels_from_sources(samples,source_2023,astral,iqtree,mapping,outgroups)
    focal=overlap["primary_common_focal_codes_sorted"]
    allshared=overlap["secondary_common_codes_sorted"]
    primary={}
    for name in (ORIGINAL_TREE,ASTRAL_TREE,IQ_FIRST_TREE):
        primary[name]=evaluate_tree_panel(samples,trees[name],focal,ROOT_PRIMARY,name)
    controls={
        "IQtree_other_individual_representatives_same_31_focal":{},
        "alternate_root_PLNY_same_31_focal":{},
        "all_37_shared_including_source_OUT":{},
    }
    for name in ("2024_IQtree_middle_sample","2024_IQtree_last_sample"):
        controls["IQtree_other_individual_representatives_same_31_focal"][name]=(
            evaluate_tree_panel(samples,trees[name],focal,ROOT_PRIMARY,name))
    for name in (ORIGINAL_TREE,ASTRAL_TREE,IQ_FIRST_TREE):
        controls["alternate_root_PLNY_same_31_focal"][name]=(
            evaluate_tree_panel(samples,trees[name],focal,ROOT_SECONDARY,name))
        controls["all_37_shared_including_source_OUT"][name]=(
            evaluate_tree_panel(samples,trees[name],allshared,ROOT_PRIMARY,name))
    return {
        "schema":"anto_v28_published_2023_2024_common_species_alternative_topology_cv",
        "original_chemical_source_DOI":ORIGINAL_STUDY_DOI,
        "original_alternative_phylogenomic_source_DOI":LATER_TREE_STUDY_DOI,
        "external_github_source_repository":LATER_REPO,
        "external_source_pinned_commit":LATER_PIN,
        "external_source_integrity":{
            k:{"original_byte_length":m["bytes"],
               "git_blob_sha1_verified":m["git_blob_sha1"]}
            for k,m in EXTERNAL.items()
        },
        "source_2023_v26_and_v27_readouts_exactly_reproduced":True,
        "source_2023_original_tree_md5_sha256_verified":True,
        "original_2023_species_rows":reconciled["Petunieae_ingroup_species_rows"],
        "2024_published_tree_taxonomy":external_admission,
        "source_intersection":overlap,
        "primary_root_outgroup":ROOT_PRIMARY,
        "secondary_root_outgroup":ROOT_SECONDARY,
        "fixed_prediction_models":list(MODELS),
        "fixed_ridge_lambda":10.0,
        "cutoffs":list(MAX_CLADE_SIZES),
        "primary_original_vs_two_published_alternative_definitions":primary,
        "predeclared_primary_gate":positive_primary_gate(primary),
        "mandatory_secondary_sensitivity":controls,
        "one_original_2023_RNA_HPLC_study_not_independent_pigment_replication":True,
        "2024_alternative_phylogenies_reuse_related_RNA_transcriptomes":True,
        "IQ_tree_individual_to_species_reduction_not_independent_published_species_tree":True,
        "models_not_reselected_to_favour_external_phylogenies":True,
        "current_chemical_amount_not_visible_colour_or_ancestral_regain":True,
        "not_cross_clade_HPLC_or_multigeneration_fitness_validation":True,
        "Cirsium_organismal_empirical_owner":"zuizui0223/aza3",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--expected-v26",type=Path,required=True)
    p.add_argument("--expected-v27",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected-v28",type=Path)
    args=p.parse_args()
    result=report(expected_v26=args.expected_v26,expected_v27=args.expected_v27)
    if args.expected_v28 and result!=json.loads(args.expected_v28.read_text(encoding="utf-8")):
        raise TreeGateError("V28_FROZEN_REAL_ALTERNATIVE_SOURCE_METRIC_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    compact={}
    for group_name,group in (
        ("primary", result["primary_original_vs_two_published_alternative_definitions"]),
        ("other_representatives",result["mandatory_secondary_sensitivity"][
            "IQtree_other_individual_representatives_same_31_focal"]),
        ("alternative_root",result["mandatory_secondary_sensitivity"][
            "alternate_root_PLNY_same_31_focal"]),
        ("outgroup_panel",result["mandatory_secondary_sensitivity"][
            "all_37_shared_including_source_OUT"]),
    ):
        compact[group_name]={}
        for name,cuts in group.items():
            compact[group_name][name]={}
            for cap,v in cuts.items():
                compact[group_name][name][cap]={
                    "tested_species":v["n_species_tested_once"],
                    "n_blocks":v["n_phylogeny_defined_blocks"],
                    "n_multispecies_blocks":v["n_multispecies_clade_blocks"],
                    "branch_three_RMSE_improvement_pct":{
                        e:round(100*(1-p["models"]["branch_three"][
                            "RMSE_relative_to_same_fold_training_mean"]),4)
                        for e,p in v["endpoints"].items()
                    },
                }
    print("ANTO_V28_VERIFIED_ALTERNATIVE_PHYLOGENY_SUMMARY",json.dumps({
        "gate":result["predeclared_primary_gate"],
        "overlap_focal":len(result["source_intersection"]["primary_common_focal_codes_sorted"]),
        "overlap_all":len(result["source_intersection"]["secondary_common_codes_sorted"]),
        "source_integrity":result["external_source_integrity"],
        "evaluations":compact,
    },sort_keys=True,separators=(",",":")))
    # Entire derived-only record, chunked to avoid truncation of long CI log lines.
    packed=json.dumps(result,sort_keys=True,separators=(",",":"),ensure_ascii=False)
    length=6000
    n=(len(packed)+length-1)//length
    for i in range(n):
        print("ANTO_V28_DERIVED_CHUNK",i+1,n,packed[i*length:(i+1)*length])


if __name__ == "__main__":
    main()
