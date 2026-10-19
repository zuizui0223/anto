"""Exact model identifiability gates: one induction is not a return forecast."""
import sys
from pathlib import Path
from fractions import Fraction as F
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_afterlife_induction_identifiability_v10b import (
    bayes_error_equal_prior_binomial, first_induction_counterexample,
    known_initial_toolkit_still_not_sufficient,
    equal_pairwise_hidden_capacity, panel
)


class InductionIdentifiabilityTests(unittest.TestCase):
    def test_indistinguishable_one_shot_for_any_independent_sample_size(self):
        r=first_induction_counterexample()
        self.assertEqual(r["first_pulse_visible_positive_both"],"1/10")
        self.assertEqual(r["first_pulse_bayes_error_any_n"],"1/2")
        for n in (1,10,40,100):
            self.assertEqual(bayes_error_equal_prior_binomial(F(1,10),F(1,10),n), F(1,2))

    def test_three_expression_event_is_distinguishable(self):
        r=first_induction_counterexample()
        p=F(r["three_in_a_row_by_horizon_A"])
        q=F(r["three_in_a_row_by_horizon_B"])
        self.assertGreater(q,p)
        self.assertGreater(q, F(8,100))
        self.assertLess(p,F(1,200))
        self.assertEqual(q-p,F(r["recurrent_event_gap_B_minus_A"]))

    def test_hypothetical_power_is_not_observed_plant_sample_size(self):
        r=first_induction_counterexample()
        risks=[r["equal_prior_model_discrimination_error"][str(n)] for n in (10,20,40,80)]
        self.assertTrue(all(x>y for x,y in zip(risks,risks[1:])))
        self.assertAlmostEqual(risks[2],0.0756957247794568,places=13)
        self.assertTrue(r["not_real_plant_sample_power"])

    def test_exact_toolkit_and_one_exposure_still_do_not_identify_continuation(self):
        r=known_initial_toolkit_still_not_sufficient()
        self.assertEqual(r["same_full_initial_toolkit_fraction"],"1/8")
        self.assertEqual(r["same_first_generation_expression"],"1/20")
        d=r["different_three_streak_after_horizon"]
        self.assertLess(F(d["low_continuation"]),F(d["high_continuation"]))

    def test_arbitrarily_long_uninduced_white_history_does_not_measure_capacity(self):
        for steps in (1,5,20,60):
            r=equal_pairwise_hidden_capacity(steps)
            self.assertEqual(r["uninduced_phenotype_trace"],"all_white")
            self.assertEqual(r["joint_toolkit_A"],"0")
            self.assertEqual(r["joint_toolkit_B"],"1/4")
            self.assertEqual(r["induced_positive_A"],"0")
            self.assertEqual(r["induced_positive_B"],"1/8")
            self.assertTrue(r["not_an_empirical_ancestry_loss_reconstruction"])

    def test_binomial_edge_cases(self):
        self.assertEqual(bayes_error_equal_prior_binomial(0,1,20),0)
        self.assertEqual(bayes_error_equal_prior_binomial(0,0,20),F(1,2))
        with self.assertRaises(ValueError):bayes_error_equal_prior_binomial(F(1,4),F(1,2),0)
        with self.assertRaises(ValueError):bayes_error_equal_prior_binomial(F(1,4),F(1,2),True)
        with self.assertRaises(ValueError):bayes_error_equal_prior_binomial(F(-1,2),F(1,2),5)
        with self.assertRaises(ValueError):first_induction_counterexample(2)
        with self.assertRaises(ValueError):equal_pairwise_hidden_capacity(0)

    def test_claim_firewall(self):
        p=panel()
        self.assertTrue(p["conditional_exact_synthetic"])
        self.assertTrue(p["independent_lineage_counts_are_only_hypothetical"])
        self.assertFalse(p["ecological_fitness_or_true_trait_regain_measured"])
        self.assertEqual(p["Cirsium_empirical_owner"],"zuizui0223/aza3")


if __name__=="__main__":
    unittest.main()
