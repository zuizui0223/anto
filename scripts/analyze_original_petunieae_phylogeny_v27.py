"""ANTO v27: original Petunieae phylogenetic-clade holdout of molecular pigment prediction.

Read ONLY five already registered, MD5/SHA256/size-pinned original author
processed files (three CSVs plus the original dated pruned Newick tree).
Retains v26's six models, lambda=10 and the same two chemical endpoints.
Clade folds are selected from tree topology BEFORE reading any outcomes;
the tree is not a claim of 59 phylogenetically independent evolution events.
No original source-row data or tree is redistributed in output.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from audit_petunieae_2023_processed_files_v26 import (
    FILES, get_official_file, registered_location, verify_original
)
from analyze_original_petunieae_2023_species_prediction_v26 import (
    MODELS, RIDGE_LAMBDA, OUTGROUP, source_real_data, evaluate, ridge_predict
)

ENDPOINTS = ("anthocyanidins_log1p_total", "flavonols_log1p_total")
MAX_CLADE_SIZES = (8, 15)
MIN_TRAIN = 10


class TreeGateError(ValueError):
    """Invalid source topology or inability to map original real species."""


class Node:
    def __init__(self, label=None, children=(), length=None):
        self.label = label
        self.children = tuple(children)
        self.length = length

    def tips(self):
        if not self.children:
            return [self.label]
        return [tip for child in self.children for tip in child.tips()]


class NewickReader:
    """Fail-closed parser of conventional quoted/unquoted Newick plus comments."""

    def __init__(self, text):
        self.text = text
        self.i = 0

    def _space(self):
        s = self.text
        while self.i < len(s):
            if s[self.i].isspace():
                self.i += 1
            elif s[self.i] == "[":
                end = s.find("]", self.i + 1)
                if end < 0:
                    raise TreeGateError("NEWICK_UNTERMINATED_COMMENT")
                self.i = end + 1
            else:
                break

    def _peek(self):
        self._space()
        return self.text[self.i] if self.i < len(self.text) else ""

    def _read_label(self):
        self._space()
        s = self.text
        if self.i >= len(s):
            return None
        if s[self.i] in ",():;":
            return None
        if s[self.i] in "'\"":
            quote = s[self.i]
            self.i += 1
            out = []
            while self.i < len(s):
                ch = s[self.i]
                self.i += 1
                if ch == quote:
                    if self.i < len(s) and s[self.i] == quote:
                        out.append(quote)
                        self.i += 1
                        continue
                    return "".join(out)
                out.append(ch)
            raise TreeGateError("NEWICK_UNTERMINATED_QUOTE")
        start = self.i
        while self.i < len(s) and not s[self.i].isspace() and s[self.i] not in ",():;[]":
            self.i += 1
        if self.i == start:
            raise TreeGateError("NEWICK_BAD_LABEL")
        return s[start:self.i]

    def _read_length(self):
        if self._peek() != ":":
            return None
        self.i += 1
        text = self._read_label()
        try:
            num = float(text)
        except (TypeError, ValueError) as exc:
            raise TreeGateError("NEWICK_INVALID_LENGTH") from exc
        if not math.isfinite(num) or num < 0:
            raise TreeGateError("NEWICK_NEGATIVE_OR_NONFINITE_LENGTH")
        return num

    def _node(self):
        if self._peek() == "(":
            self.i += 1
            children = [self._node()]
            while self._peek() == ",":
                self.i += 1
                children.append(self._node())
            if self._peek() != ")" or len(children) < 2:
                raise TreeGateError("NEWICK_UNBALANCED_OR_UNARY")
            self.i += 1
            label = self._read_label()
            return Node(label, children, self._read_length())
        label = self._read_label()
        if not label:
            raise TreeGateError("NEWICK_UNNAMED_TIP")
        return Node(label, (), self._read_length())

    def parse(self):
        tree = self._node()
        if self._peek() != ";":
            raise TreeGateError("NEWICK_MISSING_TERMINATOR")
        self.i += 1
        if self._peek() != "":
            raise TreeGateError("NEWICK_TRAILING_CONTENT")
        leaves = tree.tips()
        if not leaves or len(set(leaves)) != len(leaves):
            raise TreeGateError("NEWICK_EMPTY_OR_DUPLICATE_TIPS")
        if len(leaves) < 20 or len(leaves) > 300:
            raise TreeGateError("NEWICK_UNEXPECTED_TAXON_NUMBER")
        return tree


def parse_newick(text):
    if not isinstance(text, str) or not text.strip():
        raise TreeGateError("NEWICK_EMPTY_TEXT")
    return NewickReader(text).parse()


def original_tip_gate(tree, samples):
    """No prefix guessing, name normalization or dropping unmatched endpoints."""
    tree_labels = set(tree.tips())
    species = set(samples)
    expected_sets = (species, species | {OUTGROUP})
    if tree_labels not in expected_sets:
        raise TreeGateError(
            "ORIGINAL_TREE_SPECIES_ID_MISMATCH "
            +json.dumps({
                "tree_tip_count":len(tree_labels),
                "trait_species_count":len(species),
                "missing_from_tree": sorted(species-tree_labels),
                "extra_tree_tips": sorted(tree_labels-(species|{OUTGROUP})),
                "outgroup_in_tree": OUTGROUP in tree_labels,
                "first_original_tree_tips": sorted(tree_labels)[:12],
            },sort_keys=True))
    if len(species) != 59 or OUTGROUP in species:
        raise TreeGateError("ORIGINAL_59_INGROUP_SPECIES_DRIFT")
    return {
        "source_tree_tips_before_ingroup_prune":len(tree_labels),
        "source_tree_ingroups_reconciled":59,
        "BROW_outgroup_present_in_original_tree":OUTGROUP in tree_labels,
        "exact_original_species_code_match_no_prefix_guessing":True,
    }


def prune_outgroup(tree):
    if not tree.children:
        return None if tree.label == OUTGROUP else tree
    keep = tuple(x for child in tree.children if (x:=prune_outgroup(child)) is not None)
    if not keep:
        return None
    if len(keep) == 1:
        # Preserve monophyly after dropping outgroup (unary branch carries no split).
        return keep[0]
    return Node(tree.label, keep, tree.length)


def exact_monophyletic_partitions(tree, max_tip_count, min_train=MIN_TRAIN):
    """Deterministic maximal monophyletic components no larger than max_tip_count.

    No topology pruning uses expression, pigment, species code prefixes, or
    outcomes. Some singleton fold leaves are unavoidable in unbalanced trees
    and are identified separately, not claimed as multispecies clades.
    """
    if max_tip_count < 2:
        raise TreeGateError("CLADE_SIZE_TOO_SMALL")
    tip_list = tree.tips()
    ids = set(tip_list)
    n = len(ids)
    if len(tip_list) != n:
        raise TreeGateError("DUPLICATE_TREE_SPECIES")
    if max_tip_count > n - min_train:
        raise TreeGateError("CLADE_SIZE_LEAVES_TOO_FEW_TRAIN")
    blocks = []
    def visit(node):
        if len(node.tips()) <= max_tip_count:
            blocks.append(sorted(node.tips()))
        else:
            if not node.children:
                raise TreeGateError("TREE_NONTRIVIAL_SPLIT_UNAVAILABLE")
            for child in node.children:
                visit(child)
    visit(tree)
    flat = [tip for group in blocks for tip in group]
    if set(flat) != ids or len(flat) != n or not blocks:
        raise TreeGateError("CLADE_PARTITION_COVERAGE_OR_DISJOINTNESS")
    if any(len(b) > max_tip_count or n-len(b)<min_train for b in blocks):
        raise TreeGateError("CLADE_TRAIN_OR_MAXSIZE_VIOLATION")
    return blocks


def evaluate_folds(samples, folds, label):
    """Train-only ridge, same 6 preselected modules and 2 measured endpoints."""
    ids = sorted(samples)
    if sorted(t for b in folds for t in b) != ids:
        raise TreeGateError("PREDICTION_SPECIES_NOT_EXACTLY_ONCE")
    out = {}
    for endpoint in ENDPOINTS:
        model_results = {}
        truth = {t:samples[t]["responses"][endpoint] for t in ids}
        for model, predictors in MODELS.items():
            predictions = {}
            block_errors = {}
            for group in folds:
                train_ids = sorted(set(ids)-set(group))
                if len(train_ids)<MIN_TRAIN:
                    raise TreeGateError("INSUFFICIENT_TRAINING_SPECIES")
                train = {t:samples[t] for t in train_ids}
                block_predictions = {
                    t:ridge_predict(train, samples[t], endpoint,predictors) for t in group
                }
                for t in group:
                    if t in predictions:
                        raise TreeGateError("DUPLICATE_TEST_PREDICTION")
                    predictions[t] = block_predictions[t]
                block_errors["/".join(group)] = [truth[t]-block_predictions[t] for t in group]
            if set(predictions) != set(ids):
                raise TreeGateError("MISSING_OUT_OF_FOLD_PREDICTION")
            errors = [truth[t]-predictions[t] for t in ids]
            se = math.fsum(e*e for e in errors)
            baseline = None
            model_results[model] = {
                "RMSE_log1p_current_chemical_total": round(math.sqrt(se/len(ids)), 8),
                "MAE_log1p_current_chemical_total": round(math.fsum(abs(e) for e in errors)/len(ids),8),
                "test_species_exactly_once":len(ids),
                "training_only_parameters":True,
            }
            if model == "train_mean_baseline":
                baseline_predictions = predictions
                baseline_se = se
                baseline_rmse = math.sqrt(se/len(ids))
                baseline_block_mse = {
                    key:math.fsum(e*e for e in err)/len(err)
                    for key,err in block_errors.items()
                }
            else:
                if not baseline_rmse>0:
                    raise TreeGateError("ZERO_BASELINE_ERROR")
                model_results[model]["RMSE_relative_to_same_fold_training_mean"] = round(
                    math.sqrt(se/len(ids))/baseline_rmse, 8)
                gains = []
                multi_gains = []
                for group in folds:
                    key="/".join(group)
                    model_mse=math.fsum(e*e for e in block_errors[key])/len(group)
                    gain=baseline_block_mse[key]-model_mse
                    gains.append(gain)
                    if len(group)>1:
                        multi_gains.append(gain)
                model_results[model]["n_blocks_with_lower_MSE"] = sum(g>0 for g in gains)
                model_results[model]["n_multispecies_blocks_with_lower_MSE"] = sum(
                    g>0 for g in multi_gains)
                model_results[model]["equal_clade_block_MSE_gain"] = round(
                    math.fsum(gains)/len(gains), 8)
                model_results[model]["equal_multispecies_block_MSE_gain"] = (
                    round(math.fsum(multi_gains)/len(multi_gains), 8)
                    if multi_gains else None)
                model_results[model]["min_omit_existing_block_mean_MSE_gain"] = (
                    round(min((math.fsum(gains)-x)/(len(gains)-1) for x in gains),8)
                    if len(gains)>1 else None)
        model_results["train_mean_baseline"]["RMSE_relative_to_same_fold_training_mean"] = 1.0
        out[endpoint] = {
            "models":model_results,
            "improvement_percent_in_RMSE_vs_same_fold_mean":{
                k:round(100*(1-m["RMSE_relative_to_same_fold_training_mean"]),4)
                for k,m in model_results.items()
            },
        }
    sizes=sorted(len(x) for x in folds)
    return {
        "label":label,
        "n_species_tested_once":len(ids),
        "n_phylogeny_defined_blocks":len(folds),
        "n_multispecies_clade_blocks":sum(n>1 for n in sizes),
        "n_singleton_tips":sum(n==1 for n in sizes),
        "n_species_in_multispecies_clade_blocks":sum(n for n in sizes if n>1),
        "fold_sizes_sorted":sizes,
        "fold_membership": [{"source_tip_ids":x,"clade_tip_count":len(x)}
                            for x in folds],
        "all_folds_defined_from_original_tree_without_looking_at_pigment":True,
        "species_groups_are_mathematically_monophyletic_in_this_source_tree":True,
        "endpoints":out,
    }


def analyze_original_phylogeny(*, fetch=get_official_file, expected_v26=None):
    samples, source_reconciliation = source_real_data(fetch=fetch)
    original_tree_record = FILES["species_tree_pruned"]
    tree_body = fetch(registered_location(original_tree_record))
    verify_original(tree_body, original_tree_record)
    tree = parse_newick(tree_body.decode("utf-8-sig"))
    audit = original_tip_gate(tree,samples)
    ingroup_tree = prune_outgroup(tree)
    if ingroup_tree is None or set(ingroup_tree.tips()) != set(samples):
        raise TreeGateError("INGROUP_TREE_PRUNE_MISMATCH")
    if expected_v26 is not None:
        old = json.loads(Path(expected_v26).read_text(encoding="utf-8"))
        if old.get("schema") != "anto_v26_original_2023_petunieae_species_pigment_expression_prediction":
            raise TreeGateError("V26_EXPECTED_SCHEMA_DRIFT")
        if old.get("source_reconciled") != source_reconciliation:
            raise TreeGateError("V26_SOURCE_RECONCILIATION_DRIFT")
        for kind in ("species_LOO", "first_two_source_code_characters_NOT_accepted_genus"):
            if evaluate(samples,kind) != old["validation"][kind]:
                raise TreeGateError("V26_OLD_ORIGINAL_HOLDOUT_DRIFT")
    results = {}
    for limit in MAX_CLADE_SIZES:
        folds=exact_monophyletic_partitions(ingroup_tree,limit)
        results[f"max_clade_{limit}"]=evaluate_folds(samples,folds,f"max_tree_clade_size_{limit}")
    return {
        "schema":"anto_v27_original_petunieae_tree_defined_multispecies_blocked_prediction",
        "source_paper_doi":"10.1098/rspb.2023.0275",
        "source_osf_project":"zg9cu",
        "original_source_tree_sha256":original_tree_record["sha256"],
        "original_tree_hash_and_byte_length_verified":True,
        "source_tree_tip_audit":audit,
        "original_processed_2023_species_data":source_reconciliation,
        "model":"v26 prespecified 6 train-only ridge models, lambda 10, identical chemical endpoints",
        "phylogenetic_sensitivity_cut_sizes_predeclared":list(MAX_CLADE_SIZES),
        "v26_original_species_and_prefix_metrics_unchanged":expected_v26 is not None,
        "evaluation":results,
        "source_tree_one_topology_only_not_gene_tree_uncertainty":True,
        "derived_blocks_can_include_singleton_species_not_new_evolutionary_replication":True,
        "not_a_new_independent_plant_species_dataset":True,
        "not_a_phylogenetic_PGLS_or_ancestral_trait_state_reconstruction":True,
        "not_visible_colour_or_ancestral_reexpression_or_fitness":True,
        "no_direct_causal_genotype_function_inference":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--expected-v26",type=Path)
    p.add_argument("--expected-v27",type=Path)
    args=p.parse_args()
    r=analyze_original_phylogeny(expected_v26=args.expected_v26)
    if args.expected_v27 and json.loads(args.expected_v27.read_text(encoding="utf-8")) != r:
        raise TreeGateError("V27_FROZEN_RECEIPT_DRIFT")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("ANTO_V27_TREE_HOLDOUT_RESULT",json.dumps({
        "tip_audit":r["source_tree_tip_audit"],
        "cut_sensitivity":{
            key:{
                "n_clades":v["n_phylogeny_defined_blocks"],
                "n_multispecies_clades":v["n_multispecies_clade_blocks"],
                "singletons":v["n_singleton_tips"],
                "n_species_covered":v["n_species_tested_once"],
                "anthocyanidin_improvement_pct":{
                    k:round(100*(1-z["RMSE_relative_to_same_fold_training_mean"]),4)
                    for k,z in v["endpoints"]["anthocyanidins_log1p_total"]["models"].items()},
                "flavonol_improvement_pct":{
                    k:round(100*(1-z["RMSE_relative_to_same_fold_training_mean"]),4)
                    for k,z in v["endpoints"]["flavonols_log1p_total"]["models"].items()},
            } for key,v in r["evaluation"].items()
        },
        "trait_regain":False,"biological_model_refit":True
    },sort_keys=True))


if __name__=="__main__":
    main()
