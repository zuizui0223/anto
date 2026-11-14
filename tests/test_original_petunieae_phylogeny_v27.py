"""V27 actual-tree parsing, exact-tip and no-leakage adversarial controls."""
import copy
import hashlib
import io
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

from analyze_original_petunieae_phylogeny_v27 import (
    MAX_CLADE_SIZES, MIN_TRAIN, MODELS, Node, NewickReader, TreeGateError,
    parse_newick, original_tip_gate, prune_outgroup,
    exact_monophyletic_partitions, evaluate_folds
)


def synthetic_panel():
    out={}
    for i in range(59):
        code = f"SP{i:02d}"
        expr = {g:math.log1p(i/9 + j/7+1)
                for j,g in enumerate(sorted({gene for xs in MODELS.values() for gene in xs}))}
        out[code] = {"genes":expr, "responses":{
            "anthocyanidins_log1p_total":math.log1p((i%7)+1) + 0.15*expr["FLS"],
            "flavonols_log1p_total":math.log1p((i%9)+1) + 0.13*expr["FLS"],
        }}
    return out


def tree_for_samples(samples):
    tips=sorted(samples)
    def t(nodes):
        if len(nodes)==1:
            return nodes[0]+":0.1"
        mid=len(nodes)//2
        return "("+t(nodes[:mid])+","+t(nodes[mid:]) + "):0.2"
    return parse_newick("(BROW:0.5,"+t(tips)+"):0.3;")


class PhylogeneticSourceV27Tests(unittest.TestCase):
    def test_original_tree_gate_with_exact_59_and_real_outgroup_label(self):
        s=synthetic_panel()
        tree=tree_for_samples(s)
        self.assertEqual(len(tree.tips()),60)
        a=original_tip_gate(tree,s)
        self.assertEqual(a["source_tree_ingroups_reconciled"],59)
        self.assertTrue(a["BROW_outgroup_present_in_original_tree"])
        self.assertEqual(set(prune_outgroup(tree).tips()),set(s))
        self.assertEqual(MAX_CLADE_SIZES,(8,15))
        self.assertEqual(MIN_TRAIN,10)

    def test_tree_grouping_monophyletic_and_predefined_without_outcomes(self):
        s=synthetic_panel()
        tree=prune_outgroup(tree_for_samples(s))
        for limit in MAX_CLADE_SIZES:
            groups=exact_monophyletic_partitions(tree,limit)
            self.assertEqual(set(t for b in groups for t in b),set(s))
            self.assertEqual(sum(map(len,groups)),59)
            self.assertTrue(all(1<=len(g)<=limit for g in groups))
            copy_s=copy.deepcopy(s)
            for row in copy_s.values():
                row["responses"]["anthocyanidins_log1p_total"]=1e8
            self.assertEqual(groups,exact_monophyletic_partitions(tree,limit))
            self.assertLess(max(map(len,groups)),59)

    def test_real_original_species_code_not_in_tree_must_fail(self):
        s=synthetic_panel()
        tree=tree_for_samples(s)
        d=copy.deepcopy(s)
        d["WRONG"]=d.pop("SP00")
        with self.assertRaisesRegex(TreeGateError,"ORIGINAL_TREE_SPECIES_ID_MISMATCH"):
            original_tip_gate(tree,d)

    def test_mismatched_original_outgroup_or_duplicates_must_fail(self):
        s=synthetic_panel()
        bad=tree_for_samples(s)
        with self.assertRaisesRegex(TreeGateError,"ORIGINAL_TREE_SPECIES_ID_MISMATCH"):
            original_tip_gate(parse_newick(
                "(NOT_BROW:0.5,"+
                _fixture_ingroup(s)+"):0.1;"),s)
        with self.assertRaisesRegex(TreeGateError,"NEWICK_EMPTY_OR_DUPLICATE_TIPS"):
            parse_newick("(" + ",".join("DUP" for _ in range(21)) + ");")

    def test_newick_parser_rejects_invalid_input(self):
        for tree in (
            "(A:1,B:1)", "(A:1,B:1));", "(A:-1,B:1);",
            "(A:nan,B:1);", "(A:inf,B:1);", "(A:0,B:1",
            "(A:0);", "(A:0,B:1)[unclosed;",
        ):
            with self.assertRaises(TreeGateError,msg=tree):
                parse_newick(tree)
        with self.assertRaises(TreeGateError):
            parse_newick("")

    def test_quoted_tip_names_and_comments(self):
        taxa=["'Tip "+str(i)+"':0.1" for i in range(21)]
        tree=parse_newick("[&R]("+",".join(taxa)+")Root:0.2;")
        self.assertEqual(len(tree.tips()),21)
        self.assertEqual(tree.label,"Root")
        self.assertIn("Tip 19",tree.tips())

    def test_59_folds_predict_all_models_endpoints_and_keep_frozen_lambda(self):
        s=synthetic_panel()
        tree=prune_outgroup(tree_for_samples(s))
        fold=exact_monophyletic_partitions(tree,8)
        result=evaluate_folds(s,fold,"artificial tree only")
        self.assertEqual(result["n_species_tested_once"],59)
        self.assertEqual(result["n_phylogeny_defined_blocks"],len(fold))
        for endpoint in result["endpoints"].values():
            self.assertEqual(len(endpoint["models"]),6)
            self.assertEqual(endpoint["models"]["train_mean_baseline"][
                "RMSE_relative_to_same_fold_training_mean"],1.)
            self.assertTrue(all(x["test_species_exactly_once"]==59
                for x in endpoint["models"].values()))
            self.assertTrue(all(math.isfinite(
                x["RMSE_log1p_current_chemical_total"])
                for x in endpoint["models"].values()))
        self.assertTrue(result["all_folds_defined_from_original_tree_without_looking_at_pigment"])

    def test_duplicate_or_dropped_test_species_blocked(self):
        s=synthetic_panel()
        b=exact_monophyletic_partitions(prune_outgroup(tree_for_samples(s)),8)
        with self.assertRaisesRegex(TreeGateError,"PREDICTION_SPECIES_NOT_EXACTLY_ONCE"):
            evaluate_folds(s,b+[b[0]],"duplicate")
        with self.assertRaisesRegex(TreeGateError,"PREDICTION_SPECIES_NOT_EXACTLY_ONCE"):
            evaluate_folds(s,b[:-1],"dropped")
        with self.assertRaisesRegex(TreeGateError,"CLADE_SIZE_LEAVES_TOO_FEW_TRAIN"):
            exact_monophyletic_partitions(prune_outgroup(tree_for_samples(s)),50)

    def test_root_source_assumptions_are_not_independence(self):
        s=synthetic_panel()
        b=exact_monophyletic_partitions(prune_outgroup(tree_for_samples(s)),15)
        result=evaluate_folds(s,b,"synthetic")
        self.assertLessEqual(result["n_multispecies_clade_blocks"],
                             result["n_phylogeny_defined_blocks"])
        self.assertEqual(result["n_species_in_multispecies_clade_blocks"]+
                         result["n_singleton_tips"],59)
        self.assertNotIn("regain",result["endpoints"]["flavonols_log1p_total"])


def _fixture_ingroup(s):
    # One balanced source-like synthetic clade; no source outcomes used.
    t=tree_for_samples(s)
    return _to_newick(prune_outgroup(t))


def _to_newick(node):
    if not node.children:
        return node.label + (":"+str(node.length) if node.length is not None else "")
    return "("+",".join(_to_newick(x) for x in node.children)+")"+(
        ":"+str(node.length) if node.length is not None else "")


if __name__=="__main__":
    unittest.main()
