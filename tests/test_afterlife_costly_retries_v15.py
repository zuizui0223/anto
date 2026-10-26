"""Source-free exact synthetic tests for a retry cost's effect on trait recovery.

RED tests written before v15 implementation. Do not derive natural gene loss
or fitness from the modeled destruction of a shared *physiological* budget.
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
import anto_afterlife_retries_v14 as v14
from anto_afterlife_costly_retries_v15 import (
    RetryCost, CostlyRoutes, generator, after_silence, exact_eventual,
    alternative_ever_without_success_preemption, subcritical_ancestral_peak,
    cost_grid,
)

class CostlyV15(unittest.TestCase):
    def setUp(self):
        self.loss=v12.Loss(".02",".30",".01")
        self.rates=v12.Rates(".4",".2",self.loss)
        self.empty=RetryCost(0,0)
        self.no=CostlyRoutes(self.rates, self.empty, self.empty)

    def test_zero_cost_exactly_recovers_prior_v14_in_all_stages(self):
        for ra,rb in ((3,1.05),(.95,1.5),(.8,.9)):
            for age in (0,1,5,12):
                reference=v14.first_nonextinct_route_with_retries(self.rates,self.loss,age,ra,rb)
                new=after_silence(self.no,age,ra,rb)
                self.assertAlmostEqual(new["ancestral_first_nonextinct"],
                    reference["ancestral_first_nonextinct"],places=12)
                self.assertAlmostEqual(new["alternative_first_nonextinct"],
                    reference["alternative_first_nonextinct"],places=12)

    def test_shared_attempt_cost_creates_recovery_window_not_a_new_gene(self):
        cfg=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        outcomes=[after_silence(cfg,t,.95,1.5) for t in (0,3,5,20)]
        self.assertEqual(outcomes[0]["ancestral_first_nonextinct"],0)
        self.assertAlmostEqual(outcomes[0]["any_nonextinct"],.6268595519472424,places=11)
        self.assertAlmostEqual(outcomes[1]["any_nonextinct"],.664265770173676,places=11)
        self.assertAlmostEqual(outcomes[2]["any_nonextinct"],.6521780516467042,places=11)
        self.assertGreater(outcomes[1]["any_nonextinct"],outcomes[0]["any_nonextinct"])
        self.assertGreater(outcomes[1]["any_nonextinct"],outcomes[3]["any_nonextinct"])
        self.assertTrue(all(x["synthetic_only"] for x in outcomes))
        self.assertFalse(outcomes[1]["real_plant_recovery_demonstrated"])

    def test_calculated_interior_peak_condition_and_threshold(self):
        cfg=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        peak=subcritical_ancestral_peak(cfg,.95,1.5)
        self.assertEqual(peak["kind"],"interior")
        self.assertAlmostEqual(peak["critical_common_cost_probability"],.11164058219329057,places=11)
        self.assertAlmostEqual(peak["peak_age"],2.8193330561728502,places=11)
        p=after_silence(cfg,peak["peak_age"],.95,1.5)["any_nonextinct"]
        self.assertAlmostEqual(p,.6643708203694473,places=11)
        self.assertGreater(p,after_silence(cfg,0,.95,1.5)["any_nonextinct"])
        for delta in (-.01,.01):
            self.assertLess(after_silence(cfg,peak["peak_age"]+delta,.95,1.5)["any_nonextinct"],p)

    def test_no_cost_and_below_threshold_have_no_interior_peak(self):
        for gamma in (0,.05,.1):
            cfg=CostlyRoutes(self.rates,RetryCost(0,gamma),self.empty)
            peak=subcritical_ancestral_peak(cfg,.95,1.5)
            self.assertEqual(peak["kind"],"boundary_at_zero")
            p0=after_silence(cfg,0,.95,1.5)["any_nonextinct"]
            p3=after_silence(cfg,3,.95,1.5)["any_nonextinct"]
            self.assertGreater(p0,p3)

    def test_supercritical_ancestral_can_still_show_small_cost_window(self):
        cfg=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        p0=after_silence(cfg,0,1.05,1.5)["any_nonextinct"]
        p2=after_silence(cfg,2,1.05,1.5)["any_nonextinct"]
        self.assertGreater(p2,p0)
        self.assertGreater(p2,0.68)
        self.assertLess(p2,0.70)
        with self.assertRaisesRegex(ValueError,"subcritical"):
            subcritical_ancestral_peak(cfg,1.05,1.5)

    def test_original_2026_plant_genetic_retention_is_unmeasured(self):
        cfg=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        result=after_silence(cfg,3,.95,1.5)
        self.assertEqual(result["resource_scope"],"abstract shared reproductive/developmental capacity, not gene allele")
        self.assertFalse(result["real_plant_recovery_demonstrated"])
        self.assertFalse(result["full_pathway_function_measured"])
        self.assertFalse(result["true_genetic_fixation_estimated"])

    def test_costly_generator_and_finite_matrix_exponential_independently_agree(self):
        cfg=CostlyRoutes(self.rates,RetryCost(.2,.3),RetryCost(.1,.05))
        qa=v13.poisson_lineage_nonextinction(1.05)
        qb=v13.poisson_lineage_nonextinction(1.5)
        Q=generator(cfg,qa,qb)
        np.testing.assert_allclose(Q.sum(axis=1),0,atol=1e-13)
        self.assertTrue(np.all(Q-np.diag(np.diag(Q)) >=0))
        for age in (0,1,5,8):
            p0=v12.aged_inventory(self.loss,age)
            x=after_silence(cfg,age,1.05,1.5,10)
            finite=p0 @ expm(Q*10)
            self.assertAlmostEqual(x["ancestral_first_nonextinct"],float(finite[3]),places=12)
            self.assertAlmostEqual(x["alternative_first_nonextinct"],float(finite[4]),places=12)
            self.assertLessEqual(x["any_nonextinct"],after_silence(cfg,age,1.05,1.5)["any_nonextinct"])

    def test_eventual_solution_equals_first_step_independent_linear_system(self):
        cfg=CostlyRoutes(self.rates,RetryCost(.2,.3),RetryCost(.1,.05))
        qa=v13.poisson_lineage_nonextinction(1.05)
        qb=v13.poisson_lineage_nonextinction(1.5)
        Q=generator(cfg,qa,qb)
        exact=exact_eventual(cfg,qa,qb)
        independent=np.linalg.solve(-Q[:3,:3],Q[:3,3:5])
        np.testing.assert_allclose(independent[0],
            [exact["from_both_A"],exact["from_both_B"]],atol=1e-12)
        np.testing.assert_allclose(independent[1],
            [exact["from_A_only_A"],0],atol=1e-12)
        np.testing.assert_allclose(independent[2],
            [0,exact["from_B_only_B"]],atol=1e-12)

    def test_nonpreemptive_B_ever_higher_than_B_first_and_changes_with_A_fail_cost(self):
        for gamma in (0,.3,.8):
            cfg=CostlyRoutes(self.rates,RetryCost(0,gamma),self.empty)
            for age in (0,3,5):
                result=after_silence(cfg,age,3,1.05)
                no_preempt=alternative_ever_without_success_preemption(cfg,age,3,1.05)
                self.assertGreaterEqual(no_preempt+1e-12,result["alternative_first_nonextinct"])
        cfg0=CostlyRoutes(self.rates,self.empty,self.empty)
        cfg1=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        self.assertLess(alternative_ever_without_success_preemption(cfg1,0,3,1.05),
                        alternative_ever_without_success_preemption(cfg0,0,3,1.05))

    def test_extra_route_failure_cost_never_help_at_fixed_silence_age(self):
        low=CostlyRoutes(self.rates,RetryCost(0,0),self.empty)
        high=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        for age in (0,3,5,10):
            l=after_silence(low,age,.95,1.5)["any_nonextinct"]
            h=after_silence(high,age,.95,1.5)["any_nonextinct"]
            self.assertLessEqual(h,l+1e-12)
        # Comparing age at fixed *different* failed attempt risk is not monotone.
        self.assertGreater(after_silence(high,3,.95,1.5)["any_nonextinct"],
                           after_silence(high,0,.95,1.5)["any_nonextinct"])

    def test_if_cost_is_only_route_specific_the_total_monotonicity_remains(self):
        cfg=CostlyRoutes(self.rates,RetryCost(.7,0),RetryCost(.5,0))
        previous=2.
        for t in range(41):
            val=after_silence(cfg,t/2,.95,1.5)["any_nonextinct"]
            self.assertLessEqual(val,previous+1e-12)
            previous=val

    def test_invalid_costs_and_malformed_inputs_fail_closed(self):
        for bad in (-.1,1.1,float("nan"),float("inf"),True):
            with self.assertRaises(ValueError):
                RetryCost(bad,0)
        with self.assertRaises(ValueError):
            RetryCost(.6,.5)
        cfg=CostlyRoutes(self.rates,self.empty,self.empty)
        for bad in (-1,float("nan"),float("inf"),True):
            with self.assertRaises(ValueError):
                after_silence(cfg,bad,.95,1.5)
        with self.assertRaises(ValueError):
            after_silence(cfg,0,.95,1.5,-5)

    def test_independent_stochastic_event_race_for_costly_recovery(self):
        # Competing exponential clocks with independently drawn pre-aging
        # resource lifetimes, rather than using matrix exponential or source
        # exact_eventual probabilities. Here R_A=.95 implies q_A=0.
        rng=np.random.default_rng(20261023)
        beta=.2*v13.poisson_lineage_nonextinction(1.5)
        c,x,y=.02,.3,.01
        failed_A_common=.4*.3
        hazard_both=beta+x+c+y+failed_A_common
        chance_bonly=beta/(beta+c+y)
        n=180000
        model=CostlyRoutes(self.rates,RetryCost(0,.3),self.empty)
        for age in (0.,2.8193330561728502,5.):
            # The common/A/B resources survive the silent phase independently.
            alive_c=rng.exponential(1/c,n)>age
            alive_a=rng.exponential(1/x,n)>age
            alive_b=rng.exponential(1/y,n)>age
            both=alive_c&alive_a&alive_b
            bonly=alive_c&(~alive_a)&alive_b
            count=np.count_nonzero(bonly)
            result=np.count_nonzero(rng.random(count)<chance_bonly)
            if np.any(both):
                # While both exist, a B-success races A-background-loss
                # against common, B and futile-A common-capacity depletion.
                u=rng.random(np.count_nonzero(both))
                direct_B=(u<beta/hazard_both)
                A_is_lost=(u>=beta/hazard_both)&(
                    u<(beta+x)/hazard_both)
                result+=np.count_nonzero(direct_B)
                result+=np.count_nonzero(
                    rng.random(np.count_nonzero(A_is_lost))<chance_bonly)
            estimate=result/n
            predicted=after_silence(model,age,.95,1.5)["any_nonextinct"]
            sigma=math.sqrt(predicted*(1-predicted)/n)
            self.assertLess(abs(estimate-predicted),6*sigma+1/n)

    def test_declared_parameter_grid_and_missing_outcomes(self):
        grid=cost_grid(self.rates,self.loss)
        self.assertEqual(len(grid),210)
        self.assertTrue(all(x["synthetic_only"] for x in grid))
        self.assertTrue(any(x["cost_early_delay_raises_rescue"] for x in grid))
        self.assertTrue(any(x["cost_early_delay_lowers_rescue"] for x in grid))
        self.assertTrue(all(x["max_observed_flower_pigment_measurements"]==0 for x in grid))

if __name__=="__main__":unittest.main()
