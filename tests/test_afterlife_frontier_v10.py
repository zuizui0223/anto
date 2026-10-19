"""Internal, exactly-rational tests of synthetic trait-afterlife model."""
import sys
from fractions import Fraction as F
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_afterlife_frontier_v10 import (Rates, transitions, run,
    exact_infinite_three_streak_no_acquisition, sharp_moment_bounds,
    independent_founder_opportunity, closed_resource_success_score,
    eventual_three_streak_any_acquisition, panel)


class AfterlifeFrontierV10Tests(unittest.TestCase):
    def test_all_transition_rows_normalized(self):
        x=Rates("1/5", "3/5", "1/10", "1/7")
        for state in ("lost","latent_native","latent_acquired", "expr1_native",
                      "expr2_native","expr1_acquired","expr2_acquired",
                      "expr3_native","expr3_acquired"):
            row=transitions(state,x)
            self.assertEqual(sum((p for _,p in row), F(0)),1)
            self.assertTrue(all(p>=0 for _,p in row))

    def test_monotone_ever_and_streak_probabilities_with_12_step_conservation(self):
        path=run(F(1,8), Rates("1/5", "3/4", "2/25", "1/20"),12)
        self.assertEqual(path[0]["visibly_white_now"],1)
        self.assertEqual(path[0]["present_expression"],0)
        self.assertEqual(path[0]["latent_toolkit_unexpressed_ever"],F(1,8))
        for before,after in zip(path,path[1:]):
            self.assertLessEqual(before["first_expression_ever"],after["first_expression_ever"])
            self.assertLessEqual(before["three_consecutive_expression_ever"],after["three_consecutive_expression_ever"])
            self.assertLessEqual(after["three_consecutive_expression_ever"],after["first_expression_ever"])
            self.assertEqual(after["three_streak_native"]+after["three_streak_acquired"],
                             after["three_consecutive_expression_ever"])

    def test_independent_short_path_oracle_and_visibility_partition(self):
        a,r,d,g=F(3,5),F(3,4),F(1,10),F(1,3)
        rates=Rates(a,r,d,g)
        # Independent three-step path: L -> E1 -> E2 -> E3.
        self.assertEqual(run(1,rates,3)[-1]["three_streak_native"],a*r*r)
        # Four steps from no toolkit: D -> newly acquired L -> E1 -> E2 -> E3.
        self.assertEqual(run(0,rates,4)[-1]["three_streak_acquired"],g*a*r*r)
        path=run(F(1,8),rates,12)
        self.assertEqual(path[1]["first_expression_ever"],F(1,8)*a)
        for item in path:
            self.assertEqual(item["present_expression"]+item["visibly_white_now"],1)

    def test_no_initial_resources_no_acquisition_absolutely_no_return(self):
        for t in run(0,Rates("3/5","3/4","1/10",0),12):
            self.assertEqual(t["first_expression_ever"],0)
            self.assertEqual(t["three_consecutive_expression_ever"],0)
            self.assertEqual(t["whole_toolkit_missing_now"],1)

    def test_acquisition_can_restore_capacity_but_not_immediate_expression(self):
        rate=Rates("3/5","3/4","1/10","1/3")
        path=run(0,rate,5)
        self.assertEqual(path[1]["first_expression_ever"],0)
        self.assertEqual(path[2]["three_consecutive_expression_ever"],0)
        self.assertGreater(path[2]["first_expression_ever"],0)
        self.assertGreater(path[4]["three_streak_acquired"],0)
        self.assertEqual(path[-1]["three_streak_native"],0)

    def test_separate_flash_and_repeated_retention_with_reversed_ranking(self):
        fast=run(F(1,8),Rates("3/5","3/20","2/25"),12)
        slow=run(F(1,8),Rates("1/5","3/4","2/25"),12)
        self.assertGreater(fast[1]["first_expression_ever"],slow[1]["first_expression_ever"])
        self.assertLess(fast[-1]["three_consecutive_expression_ever"],
                        slow[-1]["three_consecutive_expression_ever"])
        self.assertLess(closed_resource_success_score(Rates("3/5","3/20","2/25")),
                        closed_resource_success_score(Rates("1/5","3/4","2/25")))

    def test_infinite_risk_closed_formula_and_finite_convergence(self):
        x=Rates("1/5","3/4","2/25")
        exact=exact_infinite_three_streak_no_acquisition(F(1,8),x)
        self.assertLess(run(F(1,8),x,60)[-1]["three_consecutive_expression_ever"],exact)
        self.assertLess(float(exact-run(F(1,8),x,200)[-1]["three_consecutive_expression_ever"]),1e-10)
        self.assertEqual(exact_infinite_three_streak_no_acquisition(F(1,8),Rates(0,"3/4",0)),0)
        self.assertEqual(exact_infinite_three_streak_no_acquisition(F(1,8),Rates("1/5",0,0)),0)
        self.assertEqual(exact_infinite_three_streak_no_acquisition(F(1,8),Rates("1/5","3/4",0)),F(1,8))
        with self.assertRaises(ValueError):
            exact_infinite_three_streak_no_acquisition(F(1,8),Rates("1/5","3/4",0,"1/20"))

    def test_sharp_moment_bounds_and_rescue_confounding(self):
        no=sharp_moment_bounds(Rates("1/5","3/4","2/25"),12)
        yes=sharp_moment_bounds(Rates("1/5","3/4","2/25","1/20"),12)
        self.assertEqual((no["initial_joint_toolkit_lower"],no["initial_joint_toolkit_upper"]),(F(0),F(1,4)))
        self.assertEqual(no["forecast_lower"],0)
        self.assertGreater(no["forecast_upper"],0)
        self.assertFalse(no["initial_no_complete_toolkit_can_recover"])
        self.assertGreater(yes["forecast_lower"],0)
        self.assertTrue(yes["initial_no_complete_toolkit_can_recover"])
        self.assertEqual(no["forecast_width"],no["forecast_upper"]-no["forecast_lower"])

    def test_asymptotic_relief_is_model_conditional(self):
        with_repair=Rates("1/5","3/4","2/25","1/20")
        self.assertEqual(eventual_three_streak_any_acquisition(0,with_repair),1)
        self.assertEqual(eventual_three_streak_any_acquisition(1,with_repair),1)
        self.assertEqual(eventual_three_streak_any_acquisition(0,Rates("1/5","3/4","2/25")),0)
        self.assertEqual(eventual_three_streak_any_acquisition(0,Rates("1/5",0,"2/25","1/20")),0)

    def test_valid_moment_point_family_and_restricted_endpoint(self):
        rates=Rates("1/3","1/2","1/4")
        point=sharp_moment_bounds(rates,7,marginals=(F(1),F(1),F(1)),
            pairs={(0,1):F(1),(0,2):F(1),(1,2):F(1)})
        self.assertEqual(point["initial_joint_toolkit_lower"],1)
        self.assertEqual(point["forecast_width"],0)
        with self.assertRaises(ValueError):
            sharp_moment_bounds(rates,4,endpoint="unknown")

    def test_region_order_changes_success_at_identical_duration(self):
        p=panel()
        self.assertNotEqual(p["same_exposure_duration_order_test"]["order_difference"],"0")

    def test_independent_founder_opportunity_is_conditional(self):
        self.assertEqual(independent_founder_opportunity(F(1,4),3),F(37,64))
        for wrong in (0,-1,1.5,True):
            with self.assertRaises(ValueError):independent_founder_opportunity(F(1,4),wrong)

    def test_invalid_inputs_and_not_confusing_event_with_fitness(self):
        with self.assertRaises(ValueError):Rates("3/4","3/4","1/2")
        with self.assertRaises(ValueError):Rates("1/4","3/4","1/3")
        with self.assertRaises(ValueError):Rates("nan","1/2","1/10")
        with self.assertRaises(ValueError):run(F(5,4),Rates("1/2","1/2",0),2)
        with self.assertRaises(ValueError):run(0,[Rates("1/2","1/2",0)],3)
        r=panel()
        self.assertTrue(r["synthetic"])
        self.assertTrue(r["three_streak_not_fitness_establishment"])
        self.assertTrue(r["external_acquisition_not_verified_introgression"])


if __name__=="__main__":unittest.main()
