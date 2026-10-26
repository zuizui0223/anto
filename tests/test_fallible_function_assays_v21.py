"""Adversarial exact verification of ANTO v21 diagnostic-error constraints.

The 22 source taxa have ZERO measured qualifying complete floral functional
assays. Mathematical calls in tests are counterfactual, not plant measurements.
"""
import itertools
import json
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from anto_fallible_function_assays_v21 import (
    false_call_bound, contrast_bounds, max_one_error_any_group,
    best_case_minimum_per_group_caps, stochastic_assay_support,
    source_panel_no_assays, panel_receipt, serialize,
)
from anto_functional_assay_bounds_v20 import function_bound


def independent_all_truth_assignments(n,tested,positive,error_budget):
    """Enumerate TRUE latent taxonomy and compare to binary reported calls."""
    calls=(1,)*positive+(0,)*(tested-positive)
    possible=[]
    for truth in itertools.product((0,1),repeat=n):
        observed_errors=sum(truth[i]!=calls[i] for i in range(tested))
        if observed_errors<=error_budget:
            possible.append(sum(truth))
    assert possible
    return min(possible),max(possible)


class V21FallibleAssayBounds(unittest.TestCase):
    def test_real_original_22_taxa_have_no_functional_tests(self):
        rows=source_panel_no_assays()
        self.assertEqual(len(rows),22)
        self.assertEqual(len({r["taxon"] for r in rows}),22)
        self.assertEqual(sum(r["flower_pigment_category"]=="UNPIGMENTED" for r in rows),7)
        self.assertEqual(sum(r["flower_pigment_category"]=="PIGMENTED" for r in rows),15)
        self.assertTrue(all(r["flower_functional_joint_toolkit_assay"] is None
                            and r["assay_admitted"] is False for r in rows))

    def test_no_assays_leave_full_difference_range(self):
        z=contrast_bounds(0,0,0,0,0,0)
        self.assertEqual((z["unpigmented_minus_pigmented_lower"],
                          z["unpigmented_minus_pigmented_upper"]),(-1,1))
        self.assertEqual(z["finite_panel_direction"],"STRICT_DIRECTION_UNIDENTIFIED")

    def test_every_small_truth_vector_agrees_with_sharp_single_group_bounds(self):
        checked=0
        for n in (2,3,4):
            for tested in range(n+1):
                for positive in range(tested+1):
                    for error in range(tested+1):
                        oracle=independent_all_truth_assignments(n,tested,positive,error)
                        predicted=false_call_bound(n,tested,positive,error)
                        self.assertEqual(
                            (predicted["true_positive_minimum"],
                             predicted["true_positive_maximum"]),oracle)
                        self.assertEqual(predicted["true_prevalence_lower"],F(oracle[0],n))
                        self.assertEqual(predicted["true_prevalence_upper"],F(oracle[1],n))
                        checked+=1
        self.assertGreaterEqual(checked,90)

    def test_equality_to_old_v20_for_all_4896_original_count_scenarios(self):
        count=0
        for ut in range(8):
            for up in range(ut+1):
                for pt in range(16):
                    for pp in range(pt+1):
                        old=function_bound(ut,up,pt,pp)
                        new=contrast_bounds(ut,up,pt,pp,0,0)
                        self.assertEqual(
                            [str(new["unpigmented_minus_pigmented_lower"]),
                             str(new["unpigmented_minus_pigmented_upper"])],
                            old["unpigmented_minus_pigmented_bound"])
                        count+=1
        self.assertEqual(count,4896)

    def test_original_eight_taxon_best_case_is_not_error_robust(self):
        original=contrast_bounds(7,7,1,0,0,0)
        with_one_in_each=contrast_bounds(7,7,1,0,1,1)
        self.assertEqual(original["unpigmented_minus_pigmented_lower"],F(1,15))
        self.assertEqual(with_one_in_each["unpigmented_minus_pigmented_lower"],-F(1,7))
        self.assertEqual(with_one_in_each["finite_panel_direction"],
                         "STRICT_DIRECTION_UNIDENTIFIED")
        one_total=max_one_error_any_group(7,7,1,0)
        self.assertEqual(one_total["difference_lower"],-F(8,105))
        self.assertFalse(one_total["positive_sign_certified"])

    def test_one_error_only_in_pigmented_requires_nine_ideal_calls(self):
        fail=contrast_bounds(7,7,1,0,0,1)
        passed=contrast_bounds(7,7,2,0,0,1)
        self.assertEqual(fail["unpigmented_minus_pigmented_lower"],0)
        self.assertEqual(passed["unpigmented_minus_pigmented_lower"],F(1,15))
        self.assertEqual(passed["finite_panel_direction"],"UNPIGMENTED_HIGHER")

    def test_one_error_only_in_unpigmented_requires_ten_ideal_calls(self):
        fail=contrast_bounds(7,7,2,0,1,0)
        passed=contrast_bounds(7,7,3,0,1,0)
        self.assertLessEqual(fail["unpigmented_minus_pigmented_lower"],0)
        self.assertEqual(passed["unpigmented_minus_pigmented_lower"],F(2,35))
        self.assertEqual(passed["finite_panel_direction"],"UNPIGMENTED_HIGHER")
        single_total=max_one_error_any_group(7,7,3,0)
        self.assertTrue(single_total["positive_sign_certified"])
        self.assertEqual(single_total["difference_lower"],F(2,35))

    def test_one_per_group_requires_eleven_not_eight_ideal_calls(self):
        fail=contrast_bounds(7,7,3,0,1,1)
        pass_=contrast_bounds(7,7,4,0,1,1)
        self.assertLessEqual(fail["unpigmented_minus_pigmented_lower"],0)
        self.assertEqual(pass_["unpigmented_minus_pigmented_lower"],F(2,35))
        self.assertEqual(pass_["finite_panel_direction"],"UNPIGMENTED_HIGHER")
        inverse=contrast_bounds(7,0,4,4,1,1)
        self.assertEqual(inverse["unpigmented_minus_pigmented_upper"],-F(2,35))
        self.assertEqual(inverse["finite_panel_direction"],"PIGMENTED_HIGHER")

    def test_full_exact_counterfactual_minimum_is_eight_nine_ten_eleven(self):
        expected={(0,0):8,(0,1):9,(1,0):10,(1,1):11,(2,1):13}
        for pair,n in expected.items():
            ans=best_case_minimum_per_group_caps(*pair)
            for direction in ("UNPIGMENTED_HIGHER","PIGMENTED_HIGHER"):
                self.assertEqual(ans[direction]["minimum_tested_taxa"],n)
                self.assertTrue(ans[direction]["possible_under_ideal_hypothetical_calls"])
                self.assertEqual(ans[direction]["minimum_allocations_unpigmented_pigmented"],
                                 [[7,n-7]])

    def test_hard_error_budgets_only_widen_sharp_identification(self):
        for ut,up,pt,pp in ((7,7,1,0),(5,3,8,2),(7,0,15,15),(0,0,0,0)):
            base=contrast_bounds(ut,up,pt,pp)
            for eU in range(min(ut,2)+1):
                for eP in range(min(pt,2)+1):
                    r=contrast_bounds(ut,up,pt,pp,eU,eP)
                    self.assertLessEqual(r["unpigmented_minus_pigmented_lower"],
                                         base["unpigmented_minus_pigmented_lower"])
                    self.assertGreaterEqual(r["unpigmented_minus_pigmented_upper"],
                                            base["unpigmented_minus_pigmented_upper"])

    def test_fully_tested_ideally_opposing_groups_can_remain_distinct_even_if_fallible(self):
        x=contrast_bounds(7,7,15,0,1,1)
        self.assertEqual(x["finite_panel_direction"],"UNPIGMENTED_HIGHER")
        self.assertEqual(x["unpigmented_minus_pigmented_lower"],F(83,105))
        self.assertEqual(x["unpigmented_minus_pigmented_upper"],1)

    def test_a_nonzero_misclassification_probability_has_no_logical_error_cap(self):
        for v in (F(19,20),F(999,1000),F(4,5)):
            model=stochastic_assay_support(v,v)
            self.assertEqual(model["full_panel_true_difference_logical_support"],
                             [-F(1),F(1)])
            self.assertFalse(model["deterministic_strict_sign_can_be_proven_from_accuracy_rates_alone"])
            self.assertFalse(model["hard_miscall_cap_derived_from_sensitivity_specificity"])
            self.assertTrue(model["all_22_functional_state_vectors_compatible_with_any_binary_call_vector"])

    def test_full_support_assumption_fail_closed_when_extremal_rates(self):
        for bad in (0,1,True,"NaN",-0.1,1.01):
            with self.assertRaises(ValueError):
                stochastic_assay_support(bad,F(19,20))
        with self.assertRaises(ValueError):
            stochastic_assay_support(F(19,20),0)

    def test_invalid_evidence_counts_and_boolean_int_confusion(self):
        invalid=((7,8,0,0,0,0),(8,7,0,0,0,0),
                 (7,7,1,0,8,0),(0,0,0,0,1,0),
                 (True,0,0,0,0,0),(-1,0,0,0,0,0))
        for x in invalid:
            with self.subTest(x=x),self.assertRaises(ValueError):
                contrast_bounds(*x)

    def test_receipt_nonclaim_and_frozen_v20_source_boundary(self):
        receipt=panel_receipt()
        self.assertEqual(receipt["source_taxa_with_new_joint_flower_function_assays"],0)
        self.assertFalse(receipt["native_ancestral_trait_regain_or_transgenerational_establishment_measured"])
        self.assertEqual(receipt["phylogenetically_independent_regain_events_added"],0)
        self.assertEqual(receipt["Cirsium_organismal_empirical_owner"],"zuizui0223/aza3")
        self.assertEqual(receipt["n_original_taxa"],22)
        self.assertEqual(receipt["zero_misclassification_reproduces_v20"]
                         ["unpigmented_minus_pigmented_lower"],F(1,15))
        self.assertEqual(receipt["ideal_eleven_taxon_certificate_if_each_group_may_have_one_miscall"]
                         ["unpigmented_minus_pigmented_lower"],F(2,35))
        encoded=serialize(receipt)
        self.assertEqual(encoded["zero_misclassification_reproduces_v20"]
                         ["unpigmented_minus_pigmented_lower"],"1/15")
        self.assertEqual(encoded["ideal_counterfactual_minima_for_error_caps"]
                         ["U1_P1"]["UNPIGMENTED_HIGHER"]["minimum_tested_taxa"],11)
        json.dumps(encoded,ensure_ascii=False)

if __name__=="__main__":
    unittest.main()
