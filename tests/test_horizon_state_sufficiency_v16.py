"""Negative controls for ANTO v16 horizon-specific state sufficiency."""
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_horizon_state_sufficiency_v16 import (
    World, binomial_two_nonextinction, qualifying_probability,
    compare_worlds, exact_brier_contrast, exact_reference_panel,
    reference_worlds,
)


def enumerated_oracle(world, horizon, endpoint):
    """Independent small-horizon explicit competent/success state tree."""
    q = world.induction
    if endpoint == "lineage":
        q *= binomial_two_nonextinction(world.daughter_viability)
    states = {"competent":F(1), "success":F(0), "lost":F(0)}
    for generation in range(horizon):
        nxt={"competent":F(0),"success":states["success"],"lost":states["lost"]}
        nxt["success"] += states["competent"]*q
        remaining=states["competent"]*(1-q)
        if generation < horizon-1:
            keep=world.retention[generation]
            nxt["competent"]+=remaining*keep
            nxt["lost"]+=remaining*(1-keep)
        else:
            nxt["lost"]+=remaining
        states=nxt
    assert sum(states.values())==1
    return states["success"]


class FutureStateSufficiencyV16(unittest.TestCase):
    def setUp(self):
        self.examples=reference_worlds()
        self.panel=exact_reference_panel()

    def test_exact_binomial_descendant_filter_supercritical_and_subcritical(self):
        self.assertEqual(binomial_two_nonextinction(F(2,3)),F(3,4))
        self.assertEqual(binomial_two_nonextinction(F(1,2)),F(0))
        self.assertEqual(binomial_two_nonextinction(F(2,5)),F(0))
        p=F(2,3); extinction=1-binomial_two_nonextinction(p)
        self.assertEqual(extinction,(1-p+p*extinction)**2)

    def test_same_immediate_visible_and_lineage_for_all_histories(self):
        for a,b in self.examples.values():
            self.assertEqual(qualifying_probability(a,1,"visible"),F(1,2))
            self.assertEqual(qualifying_probability(b,1,"visible"),F(1,2))
            self.assertEqual(qualifying_probability(a,1,"lineage"),F(3,8))
            self.assertEqual(qualifying_probability(b,1,"lineage"),F(3,8))
            for endpoint in ("visible","lineage"):
                self.assertEqual(compare_worlds(a,b,1,endpoint)
                                 ["residual_history_gain"],0)

    def test_constant_kernel_diverges_at_second_generation_for_two_endpoints(self):
        a,b=self.examples["constant_hidden_retention"]
        self.assertEqual(qualifying_probability(a,2,"visible"),F(9,16))
        self.assertEqual(qualifying_probability(b,2,"visible"),F(11,16))
        self.assertEqual(qualifying_probability(a,2,"lineage"),F(111,256))
        self.assertEqual(qualifying_probability(b,2,"lineage"),F(141,256))
        self.assertEqual(compare_worlds(a,b,2,"visible")
                         ["residual_history_gain"],F(1,256))
        self.assertEqual(compare_worlds(a,b,2,"lineage")
                         ["residual_history_gain"],F(225,65536))

    def test_delayed_kernel_passes_two_horizons_then_fails_third(self):
        a,b=self.examples["delayed_hidden_retention"]
        for endpoint in ("visible","lineage"):
            for t in (1,2):
                self.assertEqual(compare_worlds(a,b,t,endpoint)
                                 ["residual_history_gain"],0)
            self.assertGreater(compare_worlds(a,b,3,endpoint)
                               ["residual_history_gain"],0)
        self.assertEqual(compare_worlds(a,b,3,"visible")
                         ["difference_B_minus_A"],F(1,32))
        self.assertEqual(compare_worlds(a,b,3,"lineage")
                         ["difference_B_minus_A"],F(75,2048))

    def test_independent_path_tree_matches_exact_formula_at_all_small_horizons(self):
        for a,b in self.examples.values():
            for w in (a,b):
                for endpoint in ("visible","lineage"):
                    for t in (1,2,3,4,5):
                        self.assertEqual(qualifying_probability(w,t,endpoint),
                                         enumerated_oracle(w,t,endpoint))

    def test_constant_hazard_geometric_closed_form(self):
        for a in self.examples["constant_hidden_retention"]:
            for endpoint in ("visible","lineage"):
                rate=a.induction
                if endpoint=="lineage":
                    rate*=binomial_two_nonextinction(a.daughter_viability)
                r=a.retention[0]
                for t in (1,3,8):
                    x=(1-rate)*r
                    expected=rate*(1-x**t)/(1-x)
                    self.assertEqual(qualifying_probability(a,t,endpoint),expected)

    def test_no_history_gain_after_matched_transition_law(self):
        a,_=self.examples["constant_hidden_retention"]
        alias=World(a.induction,a.daughter_viability,a.retention)
        for t in (1,2,3,5,8):
            for endpoint in ("visible","lineage"):
                diff=compare_worlds(a,alias,t,endpoint)
                self.assertEqual(diff["difference_B_minus_A"],0)
                self.assertFalse(diff["horizon_coarse_state_rejected"])

    def test_horizon_monotonicity_and_lineage_is_not_expression(self):
        for a,b in self.examples.values():
            for w in (a,b):
                v=[qualifying_probability(w,t,"visible") for t in range(1,9)]
                s=[qualifying_probability(w,t,"lineage") for t in range(1,9)]
                self.assertTrue(all(F(0)<=si<=vi<=1 for vi,si in zip(v,s)))
                self.assertTrue(all(x<=y for x,y in zip(v,v[1:])))
                self.assertTrue(all(x<=y for x,y in zip(s,s[1:])))

    def test_nonextinction_is_zero_in_subcritical_env_despite_visible_return(self):
        w=World(F(1,2),F(2,5),(F(3,4),)*4)
        self.assertGreater(qualifying_probability(w,5,"visible"),0)
        self.assertEqual(qualifying_probability(w,5,"lineage"),0)

    def test_exact_brier_identity_and_no_extrapolation_to_empirical_sample(self):
        d=exact_brier_contrast(F(1,4),F(3,4))
        self.assertEqual(d["residual_history_gain"],F(1,16))
        self.assertEqual(d["coarse_brier"]-d["history_brier"],F(1,16))
        self.assertEqual(len(self.panel),2)
        self.assertEqual(set(self.panel["delayed_hidden_retention"]),{"1","2","3","5"})

    def test_fail_closed_for_mismatched_current_induction_or_demography(self):
        a,b=self.examples["constant_hidden_retention"]
        with self.assertRaisesRegex(ValueError,"matched"):
            compare_worlds(a,World(F(1,3),b.daughter_viability,b.retention),2,"visible")
        with self.assertRaisesRegex(ValueError,"matched"):
            compare_worlds(a,World(b.induction,F(3,5),b.retention),2,"lineage")

    def test_bad_probability_horizon_and_endpoint_rejected(self):
        for wrong in (-F(1,8),F(9,8),"NaN",True):
            with self.assertRaises(ValueError):World(F(1,2),F(2,3),(wrong,))
        a=self.examples["constant_hidden_retention"][0]
        for invalid in (0,-1,1.5,True,9):
            with self.assertRaises(ValueError):
                qualifying_probability(a,invalid,"visible")
        with self.assertRaisesRegex(ValueError,"endpoint"):
            qualifying_probability(a,2,"fixation")
        with self.assertRaisesRegex(ValueError,"max_horizon"):
            reference_worlds(2)


if __name__=="__main__":unittest.main()
