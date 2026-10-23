"""Evolving first-expression direction is not synonymous with demographic rescue.

These are unit/negative controls for a hypothetical single-attempt branching
process, NOT estimates of real species' offspring means or fixation.
"""
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import anto_afterlife_route_window_v12 as previous
from anto_afterlife_ecological_filter_v13 import (
    poisson_lineage_nonextinction, one_attempt_ecological_filter,
    fitness_route_crossover, nonpreemptive_alternative_probe,
    sensitivity_grid,
)

class BranchingGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.loss=previous.Loss(".02",".30",".01")
        cls.rates=previous.Rates(".40",".20",cls.loss)

    def test_branching_extinction_threshold_and_reference_roots(self):
        for mean in (0, .5, .95, 1):
            self.assertEqual(poisson_lineage_nonextinction(mean),0)
        self.assertAlmostEqual(poisson_lineage_nonextinction(1.05),.0937018370729,places=10)
        self.assertAlmostEqual(poisson_lineage_nonextinction(3),.94047979070736,places=10)
        for mean in (1.001,1.05,1.1,1.5,3,8):
            q=poisson_lineage_nonextinction(mean)
            self.assertGreater(q,0)
            self.assertLessEqual(q,1)
            self.assertAlmostEqual(q,1-math.exp(-mean*q),places=11)

    def test_first_route_is_not_branching_nonextinction(self):
        d=one_attempt_ecological_filter(self.rates,self.loss,5,3,1.05)
        self.assertAlmostEqual(d["first"]["alternative_first"],.676613795817674,places=12)
        self.assertAlmostEqual(d["first"]["ancestral_first"],.089219701320716,places=12)
        self.assertGreater(d["first"]["alternative_first"],d["first"]["ancestral_first"])
        self.assertGreater(d["one_attempt_nonextinction"]["ancestral"],d["one_attempt_nonextinction"]["alternative"])
        self.assertAlmostEqual(d["one_attempt_nonextinction"]["ancestral"],.08390932602508,places=10)
        self.assertAlmostEqual(d["one_attempt_nonextinction"]["alternative"],.06339995565699,places=10)

    def test_crossover_despite_ancestral_route_higher_lineage_persistence(self):
        crossing=fitness_route_crossover(self.rates,self.loss,3,1.05,(0,20))
        self.assertEqual(crossing["kind"],"interior")
        self.assertAlmostEqual(crossing["delay"],5.8896093398868,places=7)
        younger=one_attempt_ecological_filter(self.rates,self.loss,5,3,1.05)
        older=one_attempt_ecological_filter(self.rates,self.loss,7,3,1.05)
        self.assertGreater(younger["one_attempt_nonextinction"]["ancestral"],younger["one_attempt_nonextinction"]["alternative"])
        self.assertGreater(older["one_attempt_nonextinction"]["alternative"],older["one_attempt_nonextinction"]["ancestral"])
        self.assertLess(older["one_attempt_nonextinction"]["any"],younger["one_attempt_nonextinction"]["any"])

    def test_no_resource_loss_is_not_demographic_improvement(self):
        x=one_attempt_ecological_filter(self.rates,self.loss,0,3,1.05)
        y=one_attempt_ecological_filter(self.rates,self.loss,5,3,1.05)
        self.assertGreater(x["one_attempt_nonextinction"]["any"],y["one_attempt_nonextinction"]["any"])
        self.assertAlmostEqual(x["one_attempt_nonextinction"]["any"],.4565603536784469,places=9)
        self.assertAlmostEqual(y["one_attempt_nonextinction"]["any"],.1473092816820666,places=9)

    def test_subcritical_alternative_may_appear_but_cannot_persist(self):
        d=one_attempt_ecological_filter(self.rates,self.loss,5,3,.95)
        self.assertGreater(d["first"]["alternative_first"],.65)
        self.assertEqual(d["one_attempt_nonextinction"]["alternative"],0)
        self.assertGreater(d["one_attempt_nonextinction"]["ancestral"],0)

    def test_fail_or_no_first_mass_is_explicit(self):
        for t in (0,5,7,20):
            x=one_attempt_ecological_filter(self.rates,self.loss,t,3,1.05)
            full=x["partition"]
            self.assertAlmostEqual(sum(full.values()),1,places=12)
            self.assertTrue(all(-1e-12 <= v <=1+1e-12 for v in full.values()))
            self.assertAlmostEqual(
                full["ancestral_nonextinct"]+full["alternative_nonextinct"],
                x["one_attempt_nonextinction"]["any"],places=12)
            self.assertAlmostEqual(full["no_first_expression"],1-x["first"]["any_first"],places=12)

    def test_preemption_can_create_spurious_demographic_benefit(self):
        # Under a one-shot, "failed A ends the experiment" rule, loss of
        # subcritical A creates extra first-B events; removing A censoring
        # abolishes the gain without making B any more intrinsically viable.
        a=one_attempt_ecological_filter(self.rates,self.loss,0,.95,1.5)
        b=one_attempt_ecological_filter(self.rates,self.loss,5,.95,1.5)
        self.assertGreater(b["one_attempt_nonextinction"]["any"],a["one_attempt_nonextinction"]["any"])
        u0=nonpreemptive_alternative_probe(self.rates,self.loss,0,1.5)
        u5=nonpreemptive_alternative_probe(self.rates,self.loss,5,1.5)
        self.assertLess(u5["alternative_nonextinction"],u0["alternative_nonextinction"])
        self.assertGreater(u0["alternative_nonextinction"],a["one_attempt_nonextinction"]["alternative"])
        self.assertGreater(u5["alternative_nonextinction"],b["one_attempt_nonextinction"]["alternative"])

    def test_identical_postexpression_means_remain_first_route_direction(self):
        a=one_attempt_ecological_filter(self.rates,self.loss,5,3,3)
        qa=a["conditional_lineage_nonextinction"]["ancestral"]
        qb=a["conditional_lineage_nonextinction"]["alternative"]
        self.assertAlmostEqual(qa,qb)
        self.assertGreater(a["one_attempt_nonextinction"]["alternative"],a["one_attempt_nonextinction"]["ancestral"])

    def test_zero_first_expression_cannot_reach_branching_gate(self):
        r=previous.Rates(0,0,self.loss)
        z=one_attempt_ecological_filter(r,self.loss,12,3,3)
        self.assertEqual(z["one_attempt_nonextinction"]["any"],0)
        self.assertEqual(z["partition"]["no_first_expression"],1)

    def test_finite_observation_horizon_is_distinct_from_eventual(self):
        e=one_attempt_ecological_filter(self.rates,self.loss,5,3,1.05)
        f=one_attempt_ecological_filter(self.rates,self.loss,5,3,1.05,horizon=10)
        self.assertLess(f["one_attempt_nonextinction"]["any"],e["one_attempt_nonextinction"]["any"])

    def test_invalid_parameter_rejection(self):
        for bad in (-1,float("nan"),float("inf"),True):
            with self.assertRaises(ValueError):
                poisson_lineage_nonextinction(bad)
            with self.assertRaises(ValueError):
                one_attempt_ecological_filter(self.rates,self.loss,0,bad,2)
        with self.assertRaises(ValueError):
            fitness_route_crossover(self.rates,self.loss,3,1.05,(5,0))
        with self.assertRaises(ValueError):
            one_attempt_ecological_filter(self.rates,self.loss,-1,3,1.05)

    def test_sensitivity_grid_no_population_estimates(self):
        grid=sensitivity_grid(self.rates,self.loss)
        self.assertEqual(len(grid),540)
        self.assertEqual(sum(r["nonextinction_any"]<0 for r in grid),0)
        self.assertTrue(any(r["first_alternative_more_but_ancestral_nonextinct_more"] for r in grid))
        self.assertTrue(any(r["first_alternative_more_and_alternative_nonextinct_more"] for r in grid))
        self.assertTrue(any(r["alternative_first_positive_but_alt_branching_zero"] for r in grid))
        self.assertTrue(all(r["model_only"] for r in grid))

if __name__ == "__main__":
    unittest.main()
