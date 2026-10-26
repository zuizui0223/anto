"""Exact negative controls for ANTO v18 stationary-environment generalization."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from anto_environment_memory_generalization_v18 import (
    CapacityRates, StationaryEcology, prob, first_competent_child,
    first_visible_child, direct_horizon_two, exact_horizon_presence,
    slope_factors, threshold_for_switching, brier_history_gain,
    witness_catalog, panel,
)


def independent_h2_probability(rates, eco):
    """Explicit first-generation offspring states and grandchildren failures.

    Does not call source PGF, path extinction or analytic slope function.
    Two reproductive slots are enumerated across three outcomes each:
    none, expressing D or nonexpressing C. Offspring of each type have two
    independent competent birth slots in the next shared environment.
    """
    k,w=rates.toolkit_retention,rates.switching
    survive=F(0)
    for e0 in (0,1):
        qc0=rates.q(e0)[1]
        slot=(1-k*qc0,k*qc0*w,k*qc0*(1-w))
        for e1 in (0,1):
            qd1,qc1=rates.q(e1)
            no_grandchild=(F(1),(1-k*qd1)**2,(1-k*qc1)**2)
            failure=F(0)
            for x,y in product(range(3),repeat=2):
                failure += slot[x]*slot[y]*no_grandchild[x]*no_grandchild[y]
            survive += eco.stationary[e0]*eco.transition(e0,e1)*(1-failure)
    return survive


class EnvironmentMemoryGeneralizationV18(unittest.TestCase):
    def setUp(self):
        self.cases=witness_catalog()

    def test_fixed_occupancy_even_with_asymmetric_environment_and_negative_rho(self):
        for pi,rho in ((F(3,10),F(-1,4)), (F(3,10),F(4,5)),
                       (F(1,2),F(-4,5)),(F(7,10),F(2,3))):
            env=StationaryEcology(pi,rho)
            self.assertEqual(sum(env.stationary),1)
            for i in (0,1):
                self.assertEqual(sum(env.transition(i,j) for j in (0,1)),1)
                self.assertTrue(all(0<env.transition(i,j)<1 for j in (0,1)))
            for j in (0,1):
                self.assertEqual(sum(env.stationary[i]*env.transition(i,j)
                                     for i in (0,1)),env.stationary[j])
            self.assertEqual(env.transition(0,0)+env.transition(1,1)-1,rho)

    def test_independently_enumerated_slot_tree_matches_both_h2_formulations(self):
        for rates,pi,lo,hi in self.cases.values():
            for rho in (lo,hi):
                env=StationaryEcology(pi,rho)
                explicit=independent_h2_probability(rates,env)
                self.assertEqual(explicit,direct_horizon_two(rates,env))
                self.assertEqual(explicit,exact_horizon_presence(rates,env,2))

    def test_general_derivative_exact_across_all_five_witnesses(self):
        for rates,pi,lo,hi in self.cases.values():
            a,b=StationaryEcology(pi,lo),StationaryEcology(pi,hi)
            observed=direct_horizon_two(rates,b)-direct_horizon_two(rates,a)
            expected=(hi-lo)*slope_factors(rates,a)['rho_slope']
            self.assertEqual(observed,expected)
            self.assertEqual(slope_factors(rates,a),slope_factors(rates,b))

    def test_recover_integrated_v17_original_symmetric_analytic_reversal(self):
        r,pi,lo,hi=self.cases['v17_symmetric']
        a,b=StationaryEcology(pi,lo),StationaryEcology(pi,hi)
        self.assertEqual(first_competent_child(r,a),F(456,625))
        self.assertEqual(first_competent_child(r,b),F(456,625))
        self.assertEqual(direct_horizon_two(r,b)-direct_horizon_two(r,a),
                         F(5789184,244140625))
        self.assertEqual(slope_factors(r,a)['rho_slope'],F(723648,48828125))
        self.assertEqual(threshold_for_switching(r),F(1,2))
        self.assertEqual(first_visible_child(r,a), first_visible_child(r,b))

    def test_threshold_moves_when_reproductive_environment_is_asymmetric(self):
        rates,pi,lo,hi=self.cases['asymmetric_transition']
        self.assertEqual(threshold_for_switching(rates),F(7,13))
        a,b=StationaryEcology(pi,lo),StationaryEcology(pi,hi)
        self.assertEqual(slope_factors(rates,a)['rho_slope'],0)
        self.assertEqual(direct_horizon_two(rates,a),direct_horizon_two(rates,b))
        self.assertEqual(first_competent_child(rates,a),first_competent_child(rates,b))

    def test_two_horizon_null_can_hide_three_horizon_predictive_divergence(self):
        rates,pi,lo,hi=self.cases['asymmetric_transition']
        a,b=StationaryEcology(pi,lo),StationaryEcology(pi,hi)
        for h in (1,2):
            self.assertEqual(exact_horizon_presence(rates,a,h),
                             exact_horizon_presence(rates,b,h))
        h3a=exact_horizon_presence(rates,a,3)
        h3b=exact_horizon_presence(rates,b,3)
        self.assertEqual(h3b-h3a,F(-492874414932892640256,
                                    1513399183750152587890625))
        self.assertLess(h3b,h3a)
        self.assertEqual(brier_history_gain(h3a,h3b),(h3a-h3b)**2/4)
        self.assertEqual(brier_history_gain(exact_horizon_presence(rates,a,2),
                                            exact_horizon_presence(rates,b,2)),0)

    def test_half_switching_not_universal_null(self):
        rates,pi,lo,hi=self.cases['nonphenotype_environment']
        self.assertEqual(rates.switching,F(1,2))
        self.assertEqual(rates.q_d_a,rates.q_c_a)
        self.assertEqual(rates.q_d_b,rates.q_c_b)
        self.assertIsNone(threshold_for_switching(rates))
        a,b=StationaryEcology(pi,lo),StationaryEcology(pi,hi)
        self.assertEqual(first_visible_child(rates,a),first_visible_child(rates,b))
        self.assertEqual(first_competent_child(rates,a),first_competent_child(rates,b))
        self.assertGreater(direct_horizon_two(rates,b),direct_horizon_two(rates,a))
        self.assertEqual(direct_horizon_two(rates,b)-direct_horizon_two(rates,a),
                         F(92626944,6250000000))

    def test_saturation_can_reverse_effect_at_fixed_environmental_outputs(self):
        high,pi,lo,hi=self.cases['high_retention_saturation']
        low,pi2,lo2,hi2=self.cases['low_retention_unsaturated']
        self.assertEqual((pi,lo,hi),(pi2,lo2,hi2))
        self.assertEqual(high.q(0),low.q(0))
        self.assertEqual(high.q(1),low.q(1))
        self.assertEqual(high.switching,low.switching)
        a=StationaryEcology(pi,lo)
        top,bottom=slope_factors(high,a),slope_factors(low,a)
        self.assertGreater(top['early_C_output_contrast'],0)
        self.assertGreater(top['later_competent_output_contrast'],0)
        self.assertGreater(bottom['early_C_output_contrast'],0)
        self.assertGreater(bottom['later_competent_output_contrast'],0)
        self.assertLess(top['two_slot_saturation_factor'],0)
        self.assertGreater(bottom['two_slot_saturation_factor'],0)
        self.assertLess(top['rho_slope'],0)
        self.assertGreater(bottom['rho_slope'],0)

    def test_flat_C_fecundity_removes_horizon_two_env_memory_even_if_D_differs(self):
        rates=CapacityRates(F(9,10),F(1,2),F(2,5),F(1,2),F(4,5),F(3,10))
        a,b=StationaryEcology(F(3,10),0),StationaryEcology(F(3,10),F(4,5))
        self.assertEqual(slope_factors(rates,a)['early_C_output_contrast'],0)
        self.assertEqual(direct_horizon_two(rates,a),direct_horizon_two(rates,b))

    def test_environmentally_equal_fecundities_removes_effect_for_every_switch(self):
        for w in (F(1,10),F(1,2),F(9,10)):
            r=CapacityRates(F(1,2),F(3,5),F(1,2),F(3,5),F(4,5),w)
            a,b=StationaryEcology(F(3,10),0),StationaryEcology(F(3,10),F(4,5))
            self.assertEqual(slope_factors(r,a)['rho_slope'],0)
            for h in (1,2,3):
                self.assertEqual(exact_horizon_presence(r,a,h),
                                 exact_horizon_presence(r,b,h))

    def test_adverse_inputs_are_rejected(self):
        base=CapacityRates(F(9,10),F(2,5),F(2,5),F(9,10),F(4,5),F(1,10))
        for bad in (False,True,-.01,1.01,'NaN'):
            with self.assertRaises(ValueError):
                CapacityRates(bad,F(2,5),F(2,5),F(9,10),F(4,5),F(1,10))
        for bad in (0,1,-F(1,8),F(5,4),True):
            with self.assertRaises(ValueError):
                StationaryEcology(bad,F(1,2))
        for rho in (F(-1,2),F(1),True,'NaN'):
            with self.assertRaises(ValueError):
                StationaryEcology(F(3,10),rho)
        e=StationaryEcology(F(1,2),0)
        for h in (0,True,-1,5,1.1):
            with self.assertRaises(ValueError):
                exact_horizon_presence(base,e,h)
        with self.assertRaises(ValueError):
            e.transition(5,0)
        with self.assertRaises(ValueError):
            brier_history_gain(-F(1,2),F(1,2))

    def test_audit_receipt_preserves_source_and_nonclaim_contract(self):
        receipt=panel()
        self.assertTrue(receipt['no_plant_calibration'])
        self.assertTrue(receipt['no_natural_trait_regain'])
        self.assertTrue(receipt['no_demographic_rescue_or_genetic_fixation'])
        self.assertFalse(receipt['no_SLK_ecological_payoff_measurement'] is False)
        self.assertEqual(receipt['Cirsium_empirical_owner'],'zuizui0223/aza3')
        self.assertEqual(set(receipt['witnesses']),set(self.cases))
        for key,row in receipt['witnesses'].items():
            self.assertTrue(row['first_visible_equal'])
            self.assertEqual(row['h1_low'],row['h1_high'])
            self.assertEqual(row['h1_history_brier_gain'],0)
            self.assertEqual(row['h2_delta'],row['analytic_h2_delta'])
        self.assertEqual(receipt['witnesses']['asymmetric_transition']['h2_delta'],0)
        self.assertNotEqual(receipt['witnesses']['asymmetric_transition']['h3_delta'],0)


if __name__=='__main__':
    unittest.main()
