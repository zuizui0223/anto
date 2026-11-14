"""Fail-closed tests of v29 topology-only fine-fold and destructive-label controls."""
from __future__ import annotations

import copy
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"tests"))
from test_original_petunieae_alternative_trees_v28 import fake_source
import analyze_original_petunieae_fine_discordance_v29 as m
import analyze_original_petunieae_alternative_trees_v28 as v28
from analyze_original_petunieae_phylogeny_v27 import TreeGateError, evaluate_folds


class V29DiscordanceAndReference(unittest.TestCase):
    def test_frozen_cutset_seed_and_data_boundaries(self):
        self.assertEqual(m.FINE_CUTS,(2,3,4,5,6))
        self.assertEqual(m.NULL_CUTS,(3,5))
        self.assertEqual(m.NULL_REPS,99)
        self.assertEqual(m.NULL_SEED_BASE,20261114)
        self.assertEqual(m.EXPECTED_N,31)
        self.assertEqual(m.EXPECTED_PAIRWISE,465)
        self.assertEqual(len(m.MODELS_PLUS),5)

    def test_signature_order_invariant_and_malformed_rejected(self):
        a=[["A","B"],["C"],["D","E"]]
        b=[["E","D"],["C"],["B","A"]]
        self.assertEqual(m.signature(a),m.signature(b))
        self.assertEqual(m.coheldout_disagreements(a,b,list("ABCDE")),0)
        with self.assertRaisesRegex(TreeGateError,"DUPLICATE_TAXON"):
            m.signature([["A","B"],["A","C"]])
        with self.assertRaisesRegex(TreeGateError,"EMPTY_CLADE"):
            m.signature([["A"],[]])

    def test_coheldout_pairwise_disagreement_is_not_tree_count(self):
        a=[["A","B"],["C","D"],["E"]]
        b=[["A","C"],["B","D"],["E"]]
        self.assertEqual(m.coheldout_disagreements(a,b,list("ABCDE")),4)
        with self.assertRaisesRegex(TreeGateError,"COHELDOUT_PAIR_SPECIES_DRIFT"):
            m.coheldout_disagreements(a,b,list("ABCDEF"))

    def test_all_five_fine_topology_caps_cover_exact_real_grain(self):
        astral,iq,mapping,out,original,data=fake_source()
        trees,overlap=v28.panels_from_sources(data,original,astral,iq,mapping,out)
        taxa=overlap["primary_common_focal_codes_sorted"]
        self.assertEqual(len(taxa),31)
        for cap in m.FINE_CUTS:
            for name in m.TREES:
                folds=m.topology_folds(trees[name],taxa,cap)
                self.assertEqual(sorted(t for cell in folds for t in cell),taxa)
                self.assertTrue(all(1<=len(b)<=cap for b in folds))
                self.assertLessEqual(max(map(len,folds)),cap)

    def test_fixed_random_folds_preserve_species_exactly_once_and_sizes(self):
        astral,iq,mapping,out,original,data=fake_source()
        trees,overlap=v28.panels_from_sources(data,original,astral,iq,mapping,out)
        taxa=overlap["primary_common_focal_codes_sorted"]
        base=m.topology_folds(trees[v28.ORIGINAL_TREE],taxa,3)
        draws=list(m.same_size_label_reference(base,taxa,m.NULL_SEED_BASE+3,m.NULL_REPS))
        self.assertEqual(len(draws),99)
        sizes=sorted(map(len,base))
        for blocks in draws:
            self.assertEqual(sorted(map(len,blocks)),sizes)
            self.assertEqual(sorted(t for b in blocks for t in b),taxa)
        first_run=draws[0]
        second_run=next(m.same_size_label_reference(
            base,taxa,m.NULL_SEED_BASE+3,m.NULL_REPS))
        self.assertEqual(first_run,second_run)
        with self.assertRaisesRegex(TreeGateError,"NULL_REFERENCE_DENOMINATOR_OR_REPLICATE_DRIFT"):
            list(m.same_size_label_reference(base,taxa,m.NULL_SEED_BASE+3,98))

    def test_fixed_null_percentile_interpolation_and_unfit_values(self):
        self.assertEqual(m.quantile([1,2,3,4,5],0.5),3)
        self.assertAlmostEqual(m.quantile([1,2,3,4,5],.025),1.1)
        self.assertAlmostEqual(m.quantile([1,2,3,4,5],.975),4.9)
        for d in ([float("nan")],[],[float("inf")]):
            with self.assertRaisesRegex(TreeGateError,"BAD_FIXED_NULL_QUANTILE"):
                m.quantile(d,.5)

    def test_primary_fine_gate_cannot_hide_negative_result(self):
        astral,iq,mapping,out,original,data=fake_source()
        trees,overlap=v28.panels_from_sources(data,original,astral,iq,mapping,out)
        taxa=overlap["primary_common_focal_codes_sorted"]
        structural={}
        fine={t:{} for t in m.TREES}
        for cap in m.FINE_CUTS:
            d={t:m.topology_folds(trees[t],taxa,cap) for t in m.TREES}
            structural["max_clade_"+str(cap)]={
                "n_distinct_holdout_partitions":len({m.signature(x) for x in d.values()})
            }
            for name in m.TREES:
                fine[name]["max_clade_"+str(cap)]=evaluate_folds(
                    {t:data[t] for t in taxa},d[name],"synthetic")
        a=m.fine_gate(fine,structural)
        self.assertIn(a["verdict"],(
            "NO_EVALUABLE_PUBLISHED_FOLD_DISCORDANCE",
            "MIXED_DIRECTION_WHEN_CLADES_DISCORD",
            "EXPLORATORY_ALL_DISTINCT_FINE_PARTITIONS_POSITIVE"))
        forced=copy.deepcopy(fine)
        forced[v28.ORIGINAL_TREE]["max_clade_2"]["endpoints"][
            "anthocyanidins_log1p_total"]["models"]["branch_three"][
                "RMSE_relative_to_same_fold_training_mean"]=2
        b=m.fine_gate(forced,{k:{"n_distinct_holdout_partitions":2} for k in structural})
        self.assertEqual(b["verdict"],"MIXED_DIRECTION_WHEN_CLADES_DISCORD")
        self.assertTrue(b["branch_three_negative_distinct_designs"])

    def test_chosen_predictor_model_is_fixed_before_fine_outcomes(self):
        self.assertEqual(set(m.MODELS_PLUS),set(m.MODELS)-{"train_mean_baseline"})
        self.assertIn("branch_three",m.MODELS_PLUS)
        self.assertIn("CHI_only",m.MODELS_PLUS)
        self.assertIn("early_three",m.MODELS_PLUS)
        self.assertIn("late_two",m.MODELS_PLUS)
        self.assertIn("joint_six",m.MODELS_PLUS)


if __name__=="__main__":
    unittest.main()
