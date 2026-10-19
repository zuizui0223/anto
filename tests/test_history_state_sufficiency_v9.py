"""Exact counterfactual tests: history is useful until current modules are known."""
import sys
from fractions import Fraction
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_history_state_sufficiency_v9 import (
    compare_models,comparative_sensitivity_grid,_outcomes,_optimal_brier
)
from anto_joint_accessibility_v7 import pairwise_constrained_family,marginals,pairwise_joint


class HistoryStateSufficiencyTests(unittest.TestCase):
    def test_identical_first_and_pairwise_moments_hide_history_dependent_triples(self):
        a=pairwise_constrained_family(Fraction(0))
        b=pairwise_constrained_family(Fraction(1,4))
        self.assertEqual(marginals(a),marginals(b))
        self.assertEqual(pairwise_joint(a),pairwise_joint(b))
        self.assertEqual(a["111"],0)
        self.assertEqual(b["111"],Fraction(1,4))

    def test_exact_brier_risks_and_history_erasure_with_full_state(self):
        r=compare_models()
        b=r["bayes_optimal_brier_exact"]
        self.assertEqual(b["visible_flower_only"],"19/400")
        self.assertEqual(b["visible_plus_history"],"9/200")
        self.assertEqual(b["visible_plus_complete_functional_inventory"],"3/100")
        self.assertEqual(b["visible_plus_complete_inventory_plus_history"],"3/100")
        self.assertEqual(r["history_specific_future_event_probability"],{
            "history_even":"0","history_odd":"1/10"})
        self.assertEqual(r["residual_history_gain_after_complete_state"],0)

    def test_more_information_cannot_increase_bayes_risk(self):
        for r in comparative_sensitivity_grid()["models"].values():
            b=r["bayes_optimal_brier_decimal"]
            self.assertGreaterEqual(b["visible_flower_only"]+1e-13,
                                    b["visible_plus_history"])
            self.assertGreaterEqual(b["visible_plus_history"]+1e-13,
                                    b["visible_plus_history_plus_noisy_assays"])
            self.assertGreaterEqual(b["visible_plus_history_plus_noisy_assays"]+1e-13,
                                    b["visible_plus_complete_functional_inventory"])

    def test_uninformative_assays_do_not_improve_over_visible_only(self):
        r=compare_models(sensitivity=Fraction(1,2),specificity=Fraction(1,2))
        b=r["bayes_optimal_brier_exact"]
        self.assertEqual(b["visible_flower_only"],
                         b["visible_plus_noisy_module_assays"])
        self.assertEqual(b["visible_plus_history"],
                         b["visible_plus_history_plus_noisy_assays"])

    def test_strong_assays_approach_complete_state(self):
        r=compare_models(sensitivity=Fraction(9,10),specificity=Fraction(19,20))
        b=r["bayes_optimal_brier_decimal"]
        self.assertAlmostEqual(b["visible_plus_noisy_module_assays"],
                               0.037182992562845,places=11)
        self.assertAlmostEqual(b["visible_plus_history_plus_noisy_assays"],
                               0.033187058125219,places=11)
        self.assertLess(b["visible_plus_complete_functional_inventory"],
                        b["visible_plus_history_plus_noisy_assays"])

    def test_perfect_assay_makes_history_irrelevant(self):
        r=compare_models(sensitivity=1,specificity=1)
        b=r["bayes_optimal_brier_exact"]
        self.assertEqual(b["visible_plus_complete_functional_inventory"],
                         b["visible_plus_noisy_module_assays"])
        self.assertEqual(b["visible_plus_complete_functional_inventory"],
                         b["visible_plus_history_plus_noisy_assays"])

    def test_same_hidden_state_distribution_removes_history_association(self):
        r=compare_models(histories_differ=False)
        self.assertFalse(r["histories_differ_in_third_order_joint_inventory"])
        self.assertEqual(r["bayes_optimal_brier_exact"]["visible_flower_only"],
                         r["bayes_optimal_brier_exact"]["visible_plus_history"])

    def test_no_developmental_success_no_predictive_signal(self):
        r=compare_models(success_if_complete=0)
        self.assertTrue(all(v==0 for v in r["bayes_optimal_brier_decimal"].values()))

    def test_complete_present_resource_and_history_label_are_not_same_variable(self):
        rows=_outcomes()
        self.assertEqual(sum(x["mass"] for x in rows),1)
        self.assertEqual(_optimal_brier(rows,("inventory",)),
                         _optimal_brier(rows,("history","inventory")))
        self.assertLess(_optimal_brier(rows,("inventory",)),
                        _optimal_brier(rows,("history",)))

    def test_invalid_assays_and_alpha(self):
        for kwargs in (
            {"sensitivity":Fraction(5,4)},
            {"specificity":Fraction(-1,5)},
            {"success_if_complete":Fraction(8,7)}
        ):
            with self.assertRaises(ValueError):compare_models(**kwargs)


if __name__=="__main__":
    unittest.main()
