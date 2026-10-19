"""Exact small-generation oracle and claim boundaries for ANTO x SLK demography."""
import sys
import unittest
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from math import comb

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from anto_slk_absolute_demography_v3 import (
    AbsoluteDemography, one_founder_at_least_k,
    first_activation_after_silent_history, census_probability,
    infinite_horizon_nonextinction, censused_optimum, panel
)
from anto_slk_silent_priming_window_v2 import (
    SilentPriming, hidden_inventory, qualified_pulse
)
from anto_slk_ecological_window_v1 import SLKSlice


def brute_force_full_count_distribution(q,generations):
    """Independent full offspring histogram, not truncated PGF recurrence."""
    states={1:F(1)}
    for _ in range(generations):
        nxt={}
        for n,p in states.items():
            for k in range(2*n+1):
                z=p*comb(2*n,k)*q**k*(1-q)**(2*n-k)
                nxt[k]=nxt.get(k,F(0))+z
        if sum(nxt.values())!=1:
            raise ArithmeticError("Full distribution oracle lost mass")
        states=nxt
    return states


class ANTOAbsoluteDemographyV3Tests(unittest.TestCase):
    def setUp(self):
        self.theta=F(1,8)
        self.priming=SilentPriming(F(1,10),F(3,5),F(3,20),F(1,25),F(4,5))
        self.high=AbsoluteDemography(F(13,20),F(11,20))
        self.low=AbsoluteDemography(F(39,100),F(33,100))

    def test_small_generation_independent_full_distribution_oracle(self):
        for q in (F(2,5),F(1,2),F(13,20),F(4,5)):
            for t in range(5):
                dist=brute_force_full_count_distribution(q,t)
                for k in (1,2,3,4,5):
                    exact=sum((v for n,v in dist.items() if n>=k),F(0))
                    pgf=one_founder_at_least_k(q,t,k)
                    self.assertAlmostEqual(float(pgf),float(exact),places=14)

    def test_all_descendant_pgf_boundaries_and_smooth_precision(self):
        self.assertEqual(one_founder_at_least_k(1,0,1),1)
        self.assertEqual(one_founder_at_least_k(1,0,4),0)
        self.assertEqual(one_founder_at_least_k(1,2,4),1)
        self.assertEqual(one_founder_at_least_k(1,2,5),0)
        self.assertEqual(one_founder_at_least_k(0,1,1),0)
        a=one_founder_at_least_k(F(13,20),14,4,digits=50)
        b=one_founder_at_least_k(F(13,20),14,4,digits=95)
        self.assertAlmostEqual(float(a),float(b),places=14)
        self.assertAlmostEqual(float(b),0.7001113900921574,places=13)

    def test_relative_advantage_same_absolute_growth_different(self):
        self.assertTrue(self.high.slk.allows("all"))
        self.assertTrue(self.low.slk.allows("all"))
        self.assertEqual(self.high.slk.rare_margin,self.low.slk.rare_margin)
        self.assertEqual(self.high.candidate_success/self.high.resident_success,F(13,11))
        self.assertEqual(self.low.candidate_success/self.low.resident_success,F(13,11))
        self.assertGreater(self.high.candidate_success*2,1)
        self.assertLess(self.low.candidate_success*2,1)
        self.assertEqual(infinite_horizon_nonextinction(self.high.candidate_success),
                         F(120,169))
        self.assertEqual(infinite_horizon_nonextinction(self.low.candidate_success),0)

    def test_original_anto_hidden_state_controls_founder_supply(self):
        p=self.priming
        for wait in (0,1,5,9,14):
            stocks=hidden_inventory(self.theta,p,wait)
            activation=first_activation_after_silent_history(self.theta,p,wait)
            self.assertEqual(activation,stocks["unprimed_complete"]*p.activation_unprimed+
                                         stocks["primed_complete"]*p.activation_primed)
            self.assertEqual(activation*p.continuation**2,
                             qualified_pulse(self.theta,p,wait))
        self.assertEqual(first_activation_after_silent_history(self.theta,p,0),F(1,80))

    def test_fixed_deadline_shifts_optimum_before_readiness_peak(self):
        high=censused_optimum(self.theta,self.priming,self.high)
        self.assertEqual(high["earliest_best_pre_exposure"],6)
        at_6=high["all_waits"][6]["at_least_one_of_n_independent_lineages"]
        at_9=high["all_waits"][9]["at_least_one_of_n_independent_lineages"]
        self.assertGreater(at_6,at_9)
        self.assertAlmostEqual(float(at_6),0.6467768821574177,places=13)
        self.assertAlmostEqual(float(at_9),0.5507944786106482,places=12)
        self.assertEqual(high["all_waits"][14]["at_least_one_of_n_independent_lineages"],0)

    def test_same_relative_selection_subcritical_control(self):
        high=censused_optimum(self.theta,self.priming,self.high)
        low=censused_optimum(self.theta,self.priming,self.low)
        self.assertEqual(low["earliest_best_pre_exposure"],9)
        best_low=low["best_census_event"]
        self.assertAlmostEqual(float(best_low),0.016153014030613713,places=13)
        self.assertGreater(high["best_census_event"],35*best_low)
        for wait in (0,3,6,9,12):
            a=high["all_waits"][wait]
            b=low["all_waits"][wait]
            self.assertEqual(a["first_heritable_variant_activation_hypothetical"],
                             b["first_heritable_variant_activation_hypothetical"])
            self.assertEqual(a["relative_expected_offspring_ratio"],
                             b["relative_expected_offspring_ratio"])
            self.assertGreater(a["at_least_one_of_n_independent_lineages"],
                               b["at_least_one_of_n_independent_lineages"])

    def test_zero_capacity_zero_demographic_outcome(self):
        for model in (self.high,self.low):
            for wait in (0,6,14):
                result=census_probability(0,self.priming,model,wait)
                self.assertEqual(result["at_least_one_of_n_independent_lineages"],0)

    def test_target_more_demanding_less_probable(self):
        for target in (1,2,4,6):
            p=one_founder_at_least_k(F(13,20),8,target)
            q=one_founder_at_least_k(F(13,20),8,target+1)
            self.assertGreaterEqual(p,q)
        self.assertLess(one_founder_at_least_k(F(13,20),8,4),
                        one_founder_at_least_k(F(13,20),8,1))

    def test_independent_founder_poisson_binomial_formula(self):
        a=census_probability(self.theta,self.priming,self.high,6)
        one=a["one_original_founder_joint_event"]
        anyn=a["at_least_one_of_n_independent_lineages"]
        self.assertAlmostEqual(float(anyn),1-(1-float(one))**40,places=13)
        self.assertGreaterEqual(anyn,one)
        self.assertLess(anyn,40*one)

    def test_invalid_demography_and_ambiguous_selection_fail_closed(self):
        with self.assertRaises(ValueError):
            AbsoluteDemography(F(11,20),F(13,20))
        with self.assertRaises(ValueError):
            AbsoluteDemography(F(2,5),0)
        with self.assertRaises(ValueError):
            AbsoluteDemography(F(13,20),F(11,20),
                               slk=SLKSlice(F(4,5),F(3,2)))
        with self.assertRaises(ValueError):
            AbsoluteDemography(F(13,20),F(11,20),independent_founders=True)
        with self.assertRaises(ValueError):
            AbsoluteDemography(F(13,20),F(11,20),descendant_target=0)
        with self.assertRaises(ValueError):
            census_probability(self.theta,self.priming,self.high,15)
        with self.assertRaises(ValueError):
            census_probability(self.theta,self.priming,self.high,-1)
        with self.assertRaises(ValueError):
            census_probability(self.theta,self.priming,self.high,True)
        with self.assertRaises(ValueError):
            one_founder_at_least_k(F(1,2),-1,3)
        with self.assertRaises(ValueError):
            one_founder_at_least_k(F(1,2),1,0)
        with self.assertRaises(ValueError):
            one_founder_at_least_k(F(1,2),1,3,digits=20)

    def test_conditionally_heritable_does_not_claim_measured_reversal(self):
        out=panel()
        self.assertTrue(out["synthetic_conditional_no_field_rates"])
        self.assertTrue(out["activation_heritable_is_added_unvalidated_assumption"])
        self.assertTrue(out["relative_advantage_not_population_rescue"])
        self.assertTrue(out["finite_threshold_not_establishment_or_fixation"])
        self.assertEqual(out["Cirsium_empirical_owner"],"zuizui0223/aza3")
        self.assertEqual(out["scenarios"]["absolute_supercritical"]
                         ["earliest_best_pre_exposure"],6)
        self.assertEqual(out["scenarios"]["absolute_subcritical"]
                         ["earliest_best_pre_exposure"],9)


if __name__=="__main__":
    unittest.main()
