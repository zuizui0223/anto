"""Independent exact audits for the ANTO-SLK silent priming hypothesis."""
import sys
import unittest
from pathlib import Path
from fractions import Fraction as F

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from anto_slk_ecological_window_v1 import SLKSlice
from anto_slk_silent_priming_window_v2 import (
    SilentPriming, hidden_inventory, qualified_pulse,
    first_delay_boundary, exact_delay_increment, best_delay, identify_priming_from_three_pulses, panel
)


class SilentPrimingExactTests(unittest.TestCase):
    def setUp(self):
        self.theta=F(1,8)
        self.params=SilentPriming(F(1,10),F(3,5),F(3,20),F(1,25),F(4,5))
        self.off=SLKSlice(F(12,5),F(3,2))
        self.on=SLKSlice(F(2,5),F(3,2))

    def independent_hidden_distribution(self,wait):
        # This is iterative and deliberately independent of the closed form.
        s={"U":self.theta,"P":F(0),"D":1-self.theta}
        b,d=self.params.silent_priming,self.params.resource_loss
        for _ in range(wait):
            s={"U":s["U"]*(1-b-d),
               "P":s["P"]*(1-d)+s["U"]*b,
               "D":s["D"]+d*(s["U"]+s["P"])}
            self.assertEqual(sum(s.values()),1)
        return s

    def test_hidden_state_distribution_matches_separate_recursive_oracle(self):
        for wait in (0,1,2,5,9,30,70):
            a=hidden_inventory(self.theta,self.params,wait)
            b=self.independent_hidden_distribution(wait)
            self.assertEqual(a["unprimed_complete"],b["U"])
            self.assertEqual(a["primed_complete"],b["P"])
            self.assertEqual(a["missing_joint_toolkit"],b["D"])
            self.assertEqual(a["complete_total"],self.theta*(1-self.params.resource_loss)**wait)
            self.assertEqual(a["visible_trait_positive"],0)

    def test_pulse_requires_three_good_steps_not_just_one_expression(self):
        # Exact independent path enumeration for precisely 3 good steps.
        p=self.params
        for m in (0,1,2,7,9,15):
            stocks=self.independent_hidden_distribution(m)
            e1=stocks["U"]*p.activation_unprimed+stocks["P"]*p.activation_primed
            e2=e1*p.continuation
            e3=e2*p.continuation
            self.assertEqual(qualified_pulse(self.theta,p,m,self.off,self.on),e3)

    def test_initial_silent_priming_advantage_exact_sign_and_difference(self):
        p=self.params
        boundary=first_delay_boundary(p)
        self.assertEqual(boundary,F(71,1000))
        earliest=qualified_pulse(self.theta,p,0)
        after_one=qualified_pulse(self.theta,p,1)
        self.assertEqual(after_one-earliest,self.theta*p.continuation**2*boundary)
        self.assertEqual(earliest,F(1,125))
        self.assertEqual(after_one,F(171,12500))
        self.assertTrue(after_one>earliest)

    def test_pulse_has_delayed_peak_and_later_closes(self):
        result=best_delay(self.theta,self.params,50)
        self.assertEqual(result["best_wait"],9)
        self.assertGreater(result["best_probability"],3*result["earliest_probability"])
        self.assertGreater(result["best_probability"],result["last_probability"])
        self.assertEqual(result["best_probability"],qualified_pulse(self.theta,self.params,9))
        self.assertGreater(exact_delay_increment(self.theta,self.params,8),0)
        self.assertLess(exact_delay_increment(self.theta,self.params,9),0)
        for m in range(0,49):
            left=qualified_pulse(self.theta,self.params,m)
            right=qualified_pulse(self.theta,self.params,m+1)
            self.assertEqual(right-left,exact_delay_increment(self.theta,self.params,m))

    def test_same_visible_colour_and_same_complete_toolkit_hide_readiness(self):
        p=self.params
        q=SilentPriming(p.activation_unprimed,p.activation_primed,0,
                        p.resource_loss,p.continuation)
        a=hidden_inventory(self.theta,p,9)
        b=hidden_inventory(self.theta,q,9)
        self.assertEqual(a["complete_total"],b["complete_total"])
        self.assertEqual(a["visible_trait_positive"],b["visible_trait_positive"])
        self.assertGreater(a["primed_complete"],0)
        self.assertEqual(b["primed_complete"],0)
        self.assertGreater(qualified_pulse(self.theta,p,9),
                           4*qualified_pulse(self.theta,q,9))

    def test_no_priming_control_drops_monotonically_with_resource_loss(self):
        p=self.params
        null=SilentPriming(p.activation_unprimed,p.activation_primed,0,
                          p.resource_loss,p.continuation)
        self.assertLess(first_delay_boundary(null),0)
        ys=[qualified_pulse(self.theta,null,m) for m in range(30)]
        self.assertTrue(all(x>y for x,y in zip(ys,ys[1:])))
        for m in (0,3,9,25):
            self.assertEqual(ys[m],self.theta*p.activation_unprimed*
                             p.continuation**2*(1-p.resource_loss)**m)

    def test_equal_initial_delay_boundary_not_a_general_early_advantage(self):
        # b*(ap-au) == au*d: gen-0 and gen-1 are equal; after that loss wins.
        p=SilentPriming(F(1,5),F(2,5),F(1,10),F(1,10),F(4,5))
        self.assertEqual(first_delay_boundary(p),0)
        self.assertEqual(qualified_pulse(self.theta,p,0),
                         qualified_pulse(self.theta,p,1))
        self.assertLess(qualified_pulse(self.theta,p,2),
                        qualified_pulse(self.theta,p,1))

    def test_no_toolkit_no_phenotype_and_no_gene_rescue(self):
        for m in (0,1,9,50):
            self.assertEqual(qualified_pulse(0,self.params,m),0)
            self.assertEqual(hidden_inventory(0,self.params,m)["missing_joint_toolkit"],1)

    def test_slk_value_and_local_path_positivity_not_enough(self):
        path_only=SLKSlice(F(4,5),F(3,2))  # Phi>0 and g0>0, rare D margin<0
        self.assertTrue(path_only.allows("path"))
        self.assertFalse(path_only.allows("all"))
        self.assertEqual(qualified_pulse(self.theta,self.params,9,self.off,path_only),0)
        with self.assertRaises(ValueError):
            qualified_pulse(self.theta,self.params,9,self.on,self.on)

    def test_alternate_hazards_and_no_silent_readiness_gain(self):
        a=SilentPriming(F(1,4),F(1,4),F(1,5),F(1,20),F(4,5))
        self.assertLess(first_delay_boundary(a),0)
        self.assertEqual(qualified_pulse(self.theta,a,7),
                         self.theta*a.activation_unprimed*a.continuation**2*
                         (1-a.resource_loss)**7)
        b=SilentPriming(F(1,10),F(3,5),F(3,20),0,F(4,5))
        self.assertGreater(first_delay_boundary(b),0)
        values=[qualified_pulse(self.theta,b,m) for m in (0,1,5,20)]
        self.assertTrue(all(x<y for x,y in zip(values,values[1:])))

    def test_invalid_rates_and_input_fail_closed(self):
        for bad in (
            (F(1,2),F(3,5),F(1,10),F(1,2),F(1,2)),
            (F(1,10),F(9,10),F(1,10),F(1,5),F(1,10)),
            (F(1,10),F(3,5),F(9,10),F(1,5),F(1,10)),
            (-1,F(1,10),F(1,10),F(1,10),F(1,10)),
        ):
            with self.assertRaises(ValueError):
                SilentPriming(*bad)
        with self.assertRaises(ValueError):
            hidden_inventory(self.theta,self.params,-1)
        with self.assertRaises(ValueError):
            hidden_inventory(self.theta,self.params,True)
        with self.assertRaises(ValueError):
            hidden_inventory(F(6,5),self.params,2)
        with self.assertRaises(ValueError):
            best_delay(self.theta,self.params,-1)
        with self.assertRaises(TypeError):
            qualified_pulse(self.theta,self.params,2,favorable="good")

    def test_three_arm_readiness_identification_in_exact_model(self):
        p=self.params
        q0,q1,q2=(qualified_pulse(self.theta,p,m) for m in (0,1,2))
        identified=identify_priming_from_three_pulses(q0,q1,q2,1-p.resource_loss)
        self.assertTrue(identified["identified"])
        self.assertEqual(identified["silent_priming"],p.silent_priming)
        self.assertEqual(identified["activation_primed_over_unprimed"],
                         p.activation_primed/p.activation_unprimed)
        self.assertEqual(identified["predicted_wait3_event_for_holdout"],
                         qualified_pulse(self.theta,p,3))

    def test_two_worlds_same_first_two_cohorts_but_different_third(self):
        p=self.params
        q=SilentPriming(F(1,10),F(2,5),F(1,4),F(1,25),F(4,5))
        for m in (0,1):
            self.assertEqual(qualified_pulse(self.theta,p,m),
                             qualified_pulse(self.theta,q,m))
        self.assertNotEqual(qualified_pulse(self.theta,p,2),
                            qualified_pulse(self.theta,q,2))
        identified=identify_priming_from_three_pulses(
            *(qualified_pulse(self.theta,q,m) for m in (0,1,2)),
            1-q.resource_loss)
        self.assertEqual(identified["silent_priming"],F(1,4))
        self.assertEqual(identified["activation_primed_over_unprimed"],4)

    def test_identification_does_not_require_known_theta_or_repeat(self):
        p=self.params
        alternative=SilentPriming(
            p.activation_unprimed,p.activation_primed,p.silent_priming,
            p.resource_loss,F(3,5))
        a=identify_priming_from_three_pulses(
            *(qualified_pulse(F(1,8),p,m) for m in (0,1,2)),
            1-p.resource_loss)
        b=identify_priming_from_three_pulses(
            *(qualified_pulse(F(1,4),alternative,m) for m in (0,1,2)),
            1-p.resource_loss)
        self.assertEqual(a["silent_priming"],b["silent_priming"])
        self.assertEqual(a["activation_primed_over_unprimed"],
                         b["activation_primed_over_unprimed"])

    def test_nonidentification_and_incompatible_rate_sequences(self):
        p=self.params
        no=SilentPriming(p.activation_unprimed,p.activation_primed,0,
                         p.resource_loss,p.continuation)
        vals=[qualified_pulse(self.theta,no,m) for m in (0,1,2)]
        n=identify_priming_from_three_pulses(*vals,1-p.resource_loss)
        self.assertFalse(n["identified"])
        self.assertEqual(n["status"],"PRIMING_NOT_IDENTIFIED")
        n=identify_priming_from_three_pulses(
            vals[0],vals[1],vals[2]+F(1,1000),1-p.resource_loss)
        self.assertFalse(n["identified"])
        self.assertEqual(n["status"],"INCOMPATIBLE_WITH_CONSTANT_RATE_MODEL")
        n=identify_priming_from_three_pulses(0,0,0,1-p.resource_loss)
        self.assertFalse(n["identified"])
        self.assertEqual(n["status"],"NO_BASELINE_EVENT")

    def test_receipt_firewall(self):
        result=panel()
        self.assertTrue(result["no_new_specimens_or_empirical_rates"])
        self.assertTrue(result["priming_is_unverified_mechanism"])
        self.assertTrue(result["three_expression_snapshots_not_fitness_or_fixation"])
        self.assertTrue(result["SLK_positive_gate_is_not_population_establishment"])
        self.assertEqual(result["Cirsium_empirical_owner"],"zuizui0223/aza3")
        self.assertEqual(result["finite_optimum"]["best_wait"],9)


if __name__=="__main__":
    unittest.main()
