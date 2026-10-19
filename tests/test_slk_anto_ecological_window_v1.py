"""Exact checks for the theoretical SLK × ANTO ecological timing bridge."""
import sys
import unittest
from pathlib import Path
from fractions import Fraction as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))

from anto_afterlife_frontier_v10 import Rates
from anto_slk_ecological_window_v1 import (
    SLKSlice, opportunity_trace, sharp_opportunity_bounds, gradient_slices, panel
)


class SLKAntoWindowTests(unittest.TestCase):
    def setUp(self):
        self.rate = Rates(F(2,5), F(3,4), F(1,10))
        self.on = SLKSlice(F(2,5), F(3,2))
        self.off = SLKSlice(F(12,5), F(3,2))

    def test_strict_separated_slk_surfaces(self):
        # R(x)=x+x^2; K(x)=k*x; Delta_rare=2-k-eta.
        self.assertEqual(SLKSlice(2,F(3,2)).phi,0)
        self.assertFalse(SLKSlice(2,F(3,2)).allows("value"))
        self.assertTrue(SLKSlice(F(9,5),F(3,2)).allows("value"))
        self.assertFalse(SLKSlice(F(9,5),F(3,2)).allows("path"))
        self.assertFalse(SLKSlice(1,F(3,2)).allows("path"))
        self.assertTrue(SLKSlice(F(4,5),F(3,2)).allows("path"))
        self.assertFalse(SLKSlice(F(4,5),F(3,2)).allows("all"))
        self.assertEqual(SLKSlice(F(1,2),F(3,2)).rare_margin,0)
        self.assertFalse(SLKSlice(F(1,2),F(3,2)).allows("all"))
        self.assertTrue(self.on.allows("all"))
        self.assertLess(self.off.phi,0)

    def test_raw_anto_transition_reused_exactly_when_all_gates_open(self):
        path=opportunity_trace(1,self.rate,[self.on]*6)
        for t,row in enumerate(path):
            self.assertEqual(row["qualified_native"],row["gate_qualified_three_streak_ever"])
            self.assertEqual(row["qualified_acquired"],0)
            self.assertEqual(row["any_expression_ever"], path[t]["any_expression_ever"])
            self.assertEqual(row["gate_qualified_three_streak_ever"], row["unqualified_three_streak_ever"])
        # Exactly L -> E1 -> E2 -> E3; only one possible route in 3 steps.
        self.assertEqual(path[3]["gate_qualified_three_streak_ever"],
                         F(2,5)*F(3,4)**2)

    def test_exclusion_and_no_initial_toolkit(self):
        a=opportunity_trace(F(1,8),self.rate,[self.off]*9)
        self.assertEqual(a[-1]["gate_qualified_three_streak_ever"],0)
        self.assertGreater(a[-1]["unqualified_three_streak_ever"],0)
        none=opportunity_trace(0,self.rate,[self.on]*9)
        self.assertEqual(none[-1]["any_expression_ever"],0)
        self.assertEqual(none[-1]["gate_qualified_three_streak_ever"],0)

    def test_equal_duration_timing_reversal_under_fixed_anto_rates(self):
        horizon=14
        early=[self.on]*5+[self.off]*9
        late=[self.off]*9+[self.on]*5
        self.assertEqual(sum(c.allows("all") for c in early),5)
        self.assertEqual(sum(c.allows("all") for c in late),5)
        a=opportunity_trace(F(1,8),self.rate,early)[-1]
        b=opportunity_trace(F(1,8),self.rate,late)[-1]
        self.assertEqual(a["any_expression_ever"],b["any_expression_ever"])
        self.assertEqual(a["unqualified_three_streak_ever"],b["unqualified_three_streak_ever"])
        self.assertEqual(a["gate_qualified_three_streak_ever"], F("0.05090625"))
        self.assertAlmostEqual(float(b["gate_qualified_three_streak_ever"]),
                               0.022911090658500367, places=14)
        self.assertGreater(a["gate_qualified_three_streak_ever"],2*b["gate_qualified_three_streak_ever"])

    def test_timing_ranking_reverses_under_constant_loss(self):
        early=[self.on]*5+[self.off]*9
        late=[self.off]*9+[self.on]*5
        high=Rates(F(2,5),F(3,4),F(1,20))
        slow=Rates(F(1,20),F(19,20),F(1,20))
        f1=opportunity_trace(F(1,8),high,early)[-1]["gate_qualified_three_streak_ever"]
        f2=opportunity_trace(F(1,8),high,late)[-1]["gate_qualified_three_streak_ever"]
        s1=opportunity_trace(F(1,8),slow,early)[-1]["gate_qualified_three_streak_ever"]
        s2=opportunity_trace(F(1,8),slow,late)[-1]["gate_qualified_three_streak_ever"]
        self.assertGreater(f1,f2)
        self.assertLess(s1,s2)
        self.assertAlmostEqual(float(f1),0.0543515625,places=12)
        self.assertAlmostEqual(float(f2),0.03901085888547753,places=12)
        self.assertAlmostEqual(float(s1),0.01528609375,places=12)
        self.assertAlmostEqual(float(s2),0.03194657665380759,places=12)

    def test_shifted_value_local_and_invasion_crossings(self):
        contexts=gradient_slices()
        self.assertEqual([sum(c.allows(g) for c in contexts)
                          for g in ("value","path","all")],[12,7,5])
        first=[next(i for i,c in enumerate(contexts,1) if c.allows(g))
               for g in ("value","path","all")]
        self.assertEqual(first,[3,8,10])
        results=[opportunity_trace(F(1,8),self.rate,contexts,g)[-1]
                 ["gate_qualified_three_streak_ever"]
                 for g in ("value","path","all")]
        self.assertTrue(results[0]>results[1]>results[2]>0)
        self.assertAlmostEqual(float(results[0]),0.059839483592109374,places=14)
        self.assertAlmostEqual(float(results[1]),0.03262062077939575,places=14)
        self.assertAlmostEqual(float(results[2]),0.022911090658500367,places=14)

    def test_sharp_joint_inventory_bounds_are_affine_not_CI(self):
        env=[self.off]*9+[self.on]*5
        lohi=sharp_opportunity_bounds(self.rate,env)
        self.assertEqual((lohi["joint_toolkit_lower"],lohi["joint_toolkit_upper"]),(0,F(1,4)))
        self.assertEqual(lohi["forecast_lower"],0)
        self.assertEqual(lohi["forecast_upper"],
                         2*opportunity_trace(F(1,8),self.rate,env)[-1]["gate_qualified_three_streak_ever"])
        self.assertTrue(lohi["model_conditional_not_confidence_interval"])
        # The propagation is exactly linear in initial mixture theta.
        z=opportunity_trace(0,self.rate,env)[-1]["gate_qualified_three_streak_ever"]
        o=opportunity_trace(1,self.rate,env)[-1]["gate_qualified_three_streak_ever"]
        t=opportunity_trace(F(1,8),self.rate,env)[-1]["gate_qualified_three_streak_ever"]
        self.assertEqual(t,F(1,8)*o+F(7,8)*z)

    def test_acquisition_route_does_not_pretend_native_regain(self):
        rate=Rates(F(2,5),F(3,4),F(1,10),F(1,5))
        path=opportunity_trace(0,rate,[self.on]*7)
        self.assertGreater(path[-1]["qualified_acquired"],0)
        self.assertEqual(path[-1]["qualified_native"],0)
        self.assertEqual(path[-1]["qualified_acquired"],path[-1]["gate_qualified_three_streak_ever"])
        for item in path:
            self.assertGreaterEqual(item["unqualified_three_streak_ever"],item["gate_qualified_three_streak_ever"])

    def test_input_and_semantic_firewalls(self):
        with self.assertRaises(ValueError):
            SLKSlice(-1,0)
        with self.assertRaises(ValueError):
            self.on.allows("invalid")
        with self.assertRaises(ValueError):
            opportunity_trace(F(1,8),self.rate,[])
        with self.assertRaises(ValueError):
            opportunity_trace(F(1,8),self.rate,[self.on],"invalid")
        with self.assertRaises(ValueError):
            opportunity_trace(-F(1,8),self.rate,[self.on])
        with self.assertRaises(ValueError):
            opportunity_trace(F(1,8),[],[self.on])
        p=panel()
        self.assertTrue(p["synthetic_conditional_model"])
        self.assertTrue(p["qualified_expression_is_not_establishment"])
        self.assertFalse(p["actual_evolutionary_regain_or_fixation_observed"])
        self.assertEqual(p["Cirsium_empirical_owner"],"zuizui0223/aza3")


if __name__=="__main__":
    unittest.main()
