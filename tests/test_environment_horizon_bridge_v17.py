"""Analytic gates for v16 × v6 environmental-memory sufficiency."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from anto_environment_horizon_bridge_v17 import (
    exact_competent_horizon, direct_two_generation_survival, persistence_slope,
    visible_child_from_silent, equal_history_brier_gain, bridge_panel, probability
)
from anto_slk_markov_ecology_v6 import symmetric_markov


class EnvironmentHorizonBridgeV17(unittest.TestCase):
    def test_all_scenarios_identical_first_horizon(self):
        for p in (F(1,10),F(1,2),F(9,10)):
            for w in (F(1,10),F(1,2),F(9,10)):
                self.assertEqual(exact_competent_horizon(p,w,1),F(456,625))

    def test_first_visible_expression_not_environment_autocorrelation(self):
        for w in (F(1,10),F(1,2),F(9,10)):
            self.assertEqual(visible_child_from_silent(w),
                             symmetric_markov(F(1,10),w).first_expression_probability_from_C())
            self.assertEqual(visible_child_from_silent(w),
                             symmetric_markov(F(9,10),w).first_expression_probability_from_C())
        self.assertEqual(visible_child_from_silent(F(1,10)),F(3153,31250))
        self.assertEqual(visible_child_from_silent(F(9,10)),F(21393,31250))

    def test_direct_second_generation_formula_equals_independent_pgf_paths(self):
        for p in (F(1,10),F(1,2),F(9,10)):
            for w in (F(1,10),F(1,2),F(9,10)):
                self.assertEqual(direct_two_generation_survival(p,w),
                                 exact_competent_horizon(p,w,2))

    def test_analytical_persistence_slope_matches_exact_parameter_contrast(self):
        for w in (F(1,10),F(1,4),F(1,2),F(3,4),F(9,10)):
            left=direct_two_generation_survival(F(1,10),w)
            right=direct_two_generation_survival(F(9,10),w)
            self.assertEqual(right-left,F(4,5)*persistence_slope(w))
            self.assertEqual(persistence_slope(w),
                             -F(361824,9765625)*(2*w-1))

    def test_horizon_two_switching_reverses_environment_memory_sign(self):
        delta=F(5789184,244140625)
        a=direct_two_generation_survival(F(9,10),F(1,10))-direct_two_generation_survival(F(1,10),F(1,10))
        b=direct_two_generation_survival(F(9,10),F(9,10))-direct_two_generation_survival(F(1,10),F(9,10))
        self.assertEqual(a,delta)
        self.assertEqual(b,-delta)

    def test_half_switching_removes_memory_at_horizons_one_through_four(self):
        for horizon in (1,2,3,4):
            x=exact_competent_horizon(F(1,10),F(1,2),horizon)
            y=exact_competent_horizon(F(9,10),F(1,2),horizon)
            self.assertEqual(x,y)

    def test_first_environment_mean_alone_misses_brier_gain_later(self):
        for w in (F(1,10),F(9,10)):
            h1a=exact_competent_horizon(F(1,10),w,1)
            h1b=exact_competent_horizon(F(9,10),w,1)
            h2a=exact_competent_horizon(F(1,10),w,2)
            h2b=exact_competent_horizon(F(9,10),w,2)
            self.assertEqual(equal_history_brier_gain(h1a,h1b)["history_information_gain"],0)
            self.assertGreater(equal_history_brier_gain(h2a,h2b)["history_information_gain"],0)
        a=exact_competent_horizon(F(1,10),F(1,2),2)
        b=exact_competent_horizon(F(9,10),F(1,2),2)
        self.assertEqual(equal_history_brier_gain(a,b)["history_information_gain"],0)

    def test_v6_existing_float_oracle_matches_exact_horizons(self):
        for p in (F(1,10),F(9,10)):
            for w in (F(1,10),F(1,2),F(9,10)):
                model=symmetric_markov(p,w)
                for h in (1,2,3,4):
                    self.assertAlmostEqual(float(exact_competent_horizon(p,w,h)),
                        model.finite_capacity_survival(h)["survival_by_initial_type"][1],places=12)

    def test_predicted_probability_not_monotonic_in_environment_persistence(self):
        self.assertGreater(persistence_slope(F(1,10)),0)
        self.assertEqual(persistence_slope(F(1,2)),0)
        self.assertLess(persistence_slope(F(9,10)),0)

    def test_probability_input_validation(self):
        for x in (-.1,1.1,"NaN",True):
            with self.assertRaises(ValueError):
                probability(x,"sample")
        for bad in (0,1,-1,True):
            with self.assertRaises(ValueError):
                exact_competent_horizon(bad,F(1,10),2)
            with self.assertRaises(ValueError):
                direct_two_generation_survival(F(1,10),bad)
        for bad in (0,6,True,1.5):
            with self.assertRaises(ValueError):
                exact_competent_horizon(F(1,10),F(1,10),bad)

    def test_history_brier_identity_on_declared_binary_endpoint(self):
        d=equal_history_brier_gain(F(1,4),F(3,4))
        self.assertEqual(d["history_information_gain"],F(1,16))
        self.assertEqual(d["coarse_brier"]-d["history_brier"],F(1,16))

    def test_frozen_receipt_contains_no_empirical_promotions(self):
        r=bridge_panel()
        self.assertEqual(r["same_first_generation_probability"],F(456,625))
        self.assertFalse(r["ancestral_trait_regain_measured"])
        self.assertFalse(r["SLK_payoff_calibrated"])
        self.assertTrue(r["no_natural_plant_parameters_estimated"])
        self.assertTrue(r["no_permanent_genealogy_nonextinction_claim"])
        self.assertTrue(r["known_markov_transition_law_closes_history_proxy_by_construction"])
        self.assertEqual(r["Cirsium_empirical_owner"],"zuizui0223/aza3")
        a=r["scenarios"]["1/10"]
        b=r["scenarios"]["9/10"]
        c=r["scenarios"]["1/2"]
        self.assertEqual(a["difference_p9_minus_p1_at_h2"],F(5789184,244140625))
        self.assertEqual(b["difference_p9_minus_p1_at_h2"],-F(5789184,244140625))
        self.assertEqual(c["difference_p9_minus_p1_at_h2"],0)
        self.assertEqual(a["brier_gain_at_h1"],0)
        self.assertEqual(a["brier_gain_at_h2"],(F(5789184,244140625)**2)/4)


if __name__=="__main__":unittest.main()
