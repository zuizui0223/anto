"""ANTO v14: first *surviving* route when failed expressions do not preempt retry.

No empirical plant offspring means or allele fixation are inferred.
"""
import math
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import anto_afterlife_route_window_v12 as v12
import anto_afterlife_ecological_filter_v13 as v13
from anto_afterlife_retries_v14 import (
    effective_surviving_attempt_rates,
    first_nonextinct_route_with_retries,
    eventual_alternative_nonextinct_without_preemption,
    crossing_age_for_persistent_route,
    coupling_grid,
)

class ReexpressionRetriesV14(unittest.TestCase):
    def setUp(self):
        self.loss=v12.Loss(".02",".30",".01")
        self.rates=v12.Rates(".40",".20",self.loss)

    def test_poisson_branch_gate_thins_attempt_clocks(self):
        active=effective_surviving_attempt_rates(self.rates,3,1.05)
        self.assertAlmostEqual(float(active.ancestral),.4*.9404797907073595,places=12)
        self.assertAlmostEqual(float(active.alternative),.2*.09370183707290118,places=12)
        self.assertEqual(active.loss,self.loss)

    def test_stated_four_ages_and_first_success_partition(self):
        expected={
           0:(.5263877622003944,.18496692749843047),
           5:(.10641340226760972,.29261774005853564),
           7:(.05613899134573933,.2918593055544138),
           20:(.000878830686789013,.21074322275299479)
        }
        for age,(a,b) in expected.items():
            r=first_nonextinct_route_with_retries(self.rates,self.loss,age,3,1.05)
            self.assertAlmostEqual(r["ancestral_first_nonextinct"],a,places=10)
            self.assertAlmostEqual(r["alternative_first_nonextinct"],b,places=10)
            self.assertAlmostEqual(r["any_nonextinct"],a+b,places=10)
            self.assertAlmostEqual(r["none_by_endpoint"]+r["any_nonextinct"],1,places=12)

    def test_route_dominance_crosses_earlier_than_one_shot_v13(self):
        r=crossing_age_for_persistent_route(self.rates,self.loss,3,1.05,(0,20))
        self.assertEqual(r["kind"],"interior")
        self.assertAlmostEqual(r["delay"],2.17240158406441,places=8)
        self.assertLess(r["delay"],5.8896093398868)
        before=first_nonextinct_route_with_retries(self.rates,self.loss,0,3,1.05)
        after=first_nonextinct_route_with_retries(self.rates,self.loss,5,3,1.05)
        self.assertGreater(before["ancestral_first_nonextinct"],before["alternative_first_nonextinct"])
        self.assertLess(after["ancestral_first_nonextinct"],after["alternative_first_nonextinct"])

    def test_failed_expression_retry_must_not_preempt_b_route(self):
        # A is below demographic replacement even though it activates rapidly.
        a=v13.one_attempt_ecological_filter(self.rates,self.loss,0,.95,1.5)
        b=v13.one_attempt_ecological_filter(self.rates,self.loss,5,.95,1.5)
        self.assertGreater(b["one_attempt_nonextinction"]["any"],a["one_attempt_nonextinction"]["any"])
        new0=first_nonextinct_route_with_retries(self.rates,self.loss,0,.95,1.5)
        new5=first_nonextinct_route_with_retries(self.rates,self.loss,5,.95,1.5)
        self.assertAlmostEqual(new0["any_nonextinct"],.7953089293058951,places=11)
        self.assertAlmostEqual(new5["any_nonextinct"],.6845287391756563,places=11)
        self.assertLess(new5["any_nonextinct"],new0["any_nonextinct"])
        self.assertEqual(new0["ancestral_first_nonextinct"],0)

    def test_total_rescue_cannot_improve_from_delaying_reexpression_in_this_model(self):
        # Coupling: source inventories at t+dt are subsets of those at t.
        for mean_a,mean_b in ((.95,1.5),(3,1.05),(1.05,1.05),(5,5),(0.9,0.8)):
            last=2
            for age in range(41):
                r=first_nonextinct_route_with_retries(
                    self.rates,self.loss,age*.5,mean_a,mean_b)
                self.assertLessEqual(r["any_nonextinct"],last+1e-12)
                last=r["any_nonextinct"]

    def test_alternative_ever_success_decreases_even_when_B_first_increases(self):
        for age in (0,5,20):
            r=first_nonextinct_route_with_retries(self.rates,self.loss,age,3,1.05)
            B=eventual_alternative_nonextinct_without_preemption(
                self.rates,self.loss,age,1.05)
            self.assertGreaterEqual(B+1e-12,r["alternative_first_nonextinct"])
        r0=first_nonextinct_route_with_retries(self.rates,self.loss,0,3,1.05)
        r5=first_nonextinct_route_with_retries(self.rates,self.loss,5,3,1.05)
        self.assertGreater(r5["alternative_first_nonextinct"],r0["alternative_first_nonextinct"])
        self.assertGreater(
            eventual_alternative_nonextinct_without_preemption(self.rates,self.loss,0,1.05),
            eventual_alternative_nonextinct_without_preemption(self.rates,self.loss,5,1.05))
        self.assertAlmostEqual(
            eventual_alternative_nonextinct_without_preemption(self.rates,self.loss,0,1.05),
            .38449376581790445,places=10)

    def test_finite_followup_has_less_nonextinction_started_than_eventual(self):
        a=first_nonextinct_route_with_retries(self.rates,self.loss,5,3,1.05)
        b=first_nonextinct_route_with_retries(self.rates,self.loss,5,3,1.05,horizon=10)
        self.assertLess(b["any_nonextinct"],a["any_nonextinct"])

    def test_first_success_matches_independent_absorbing_chain(self):
        active=effective_surviving_attempt_rates(self.rates,3,1.05)
        q=v12.generator(active)
        assert q.shape==(6,6)
        for age in (0,2,5,15):
            p0=v12.aged_inventory(self.loss,age)
            m=np.linalg.solve(-q[:3,:3],q[:3,3:5])
            res=p0[:3] @ m
            actual=first_nonextinct_route_with_retries(self.rates,self.loss,age,3,1.05)
            np.testing.assert_allclose(
                res,[actual["ancestral_first_nonextinct"],
                     actual["alternative_first_nonextinct"]],atol=1e-12)

    def test_independent_poisson_thinning_lifetime_simulation(self):
        # Independently sample pre-aged C/A/B inventories and exponential
        # successful expression event clocks. This does not simulate
        # real reproductive output; it validates exact rates.
        rng=np.random.default_rng(20261023)
        n=130000
        loss=(.02,.3,.01)
        qa=v13.poisson_lineage_nonextinction(3)
        qb=v13.poisson_lineage_nonextinction(1.05)
        for age in (0,5,12):
            live=[rng.exponential(1/k,n)>age for k in loss]
            has_c,has_a,has_b=live
            t_c=rng.exponential(1/loss[0],n)
            t_a_loss=rng.exponential(1/loss[1],n)
            t_b_loss=rng.exponential(1/loss[2],n)
            t_a=rng.exponential(1/(.4*qa),n)
            t_b=rng.exponential(1/(.2*qb),n)
            t_a[~(has_c&has_a&(t_a<t_c)&(t_a<t_a_loss))]=np.inf
            t_b[~(has_c&has_b&(t_b<t_c)&(t_b<t_b_loss))]=np.inf
            outcome=((t_a<t_b)&np.isfinite(t_a)).mean(),((t_b<t_a)&np.isfinite(t_b)).mean()
            exact=first_nonextinct_route_with_retries(self.rates,self.loss,age,3,1.05)
            for empirical,predicted in zip(outcome,(
                exact["ancestral_first_nonextinct"],exact["alternative_first_nonextinct"])):
                self.assertLess(abs(empirical-predicted),
                                6*math.sqrt(max(1e-12,predicted*(1-predicted))/n)+1/n)

    def test_zero_fitness_and_zero_resource_loss_limits(self):
        r=first_nonextinct_route_with_retries(self.rates,self.loss,0,0.8,0.9)
        self.assertEqual(r["any_nonextinct"],0)
        no_loss=v12.Loss(0,0,0)
        rates=v12.Rates(.4,.2,no_loss)
        value=first_nonextinct_route_with_retries(rates,no_loss,1000,3,1.05)
        self.assertAlmostEqual(value["any_nonextinct"],1)
        self.assertAlmostEqual(value["ancestral_first_nonextinct"]+value["alternative_first_nonextinct"],1)
        self.assertEqual(
            eventual_alternative_nonextinct_without_preemption(rates,no_loss,1000,1.05),1)

    def test_bad_parameters_rejected(self):
        for invalid in (-1,float("inf"),float("nan"),True):
            with self.assertRaises(ValueError):
                first_nonextinct_route_with_retries(self.rates,self.loss,0,invalid,1.1)
            with self.assertRaises(ValueError):
                eventual_alternative_nonextinct_without_preemption(
                    self.rates,self.loss,0,invalid)
        with self.assertRaises(ValueError):
            crossing_age_for_persistent_route(self.rates,self.loss,3,1.05,(12,2))

    def test_sensitivity_grid_retains_model_scope(self):
        result=coupling_grid(self.rates,self.loss)
        self.assertEqual(len(result),540)
        self.assertTrue(all(x["synthetic_only"] for x in result))
        self.assertTrue(any(x["one_shot_total_increase_but_retry_total_decrease"]
                            for x in result))
        self.assertTrue(all(0<=x["retry_any"]<=1 for x in result))
        self.assertTrue(all(x["retry_any"]+1e-12 >= x["one_shot_any"] for x in result))
        self.assertTrue(all(x["retry_any"]<=x["probability_any_initial_capacity"]+1e-12
                            for x in result))
if __name__=="__main__":
    unittest.main()
