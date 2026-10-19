"""Exact three-module accessibility and migration-history counterexamples."""
import math
from fractions import Fraction
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_joint_accessibility_v7 import (
    MASKS,co_retained_distribution,complete_retention_probability,
    decorrelate,expected_expressible_subsets_upper_bound,
    feasible_triple_bounds,founders_at_least_one_complete,
    independent_marginal_distribution,joint_from_moments,
    joint_inventory_entropy_bits,marginals,migration_order_counterfactual,
    observed_modules_report,pairwise_constrained_family,pairwise_joint,
    probability_required_modules_after_decay,third_order_central_moment,
    validate_joint,exact_report,
)

P_HALF=(Fraction(1,2),)*3
P_PAIRS={(0,1):Fraction(1,4),(0,2):Fraction(1,4),(1,2):Fraction(1,4)}


class TripleAccessibiltyTest(unittest.TestCase):
    def test_sharp_bounds_given_identical_first_and_second_moments(self):
        self.assertEqual(feasible_triple_bounds(P_HALF,P_PAIRS),
                         (Fraction(0),Fraction(1,4)))
        for theta in (Fraction(0),Fraction(1,32),Fraction(1,8),Fraction(1,4)):
            p=joint_from_moments(P_HALF,P_PAIRS,theta)
            self.assertEqual(complete_retention_probability(p),theta)
            self.assertEqual(marginals(p),P_HALF)
            self.assertEqual(pairwise_joint(p),P_PAIRS)

    def test_zero_to_quarter_identical_lower_order_observations(self):
        a=pairwise_constrained_family(0)
        b=pairwise_constrained_family(Fraction(1,8))
        c=pairwise_constrained_family(Fraction(1,4))
        self.assertEqual(marginals(a),marginals(b))
        self.assertEqual(marginals(a),marginals(c))
        self.assertEqual(pairwise_joint(a),pairwise_joint(b))
        self.assertEqual(pairwise_joint(a),pairwise_joint(c))
        self.assertEqual(tuple(complete_retention_probability(x) for x in (a,b,c)),
                         (Fraction(0),Fraction(1,8),Fraction(1,4)))

    def test_third_order_interaction_not_predictable_from_pairs(self):
        self.assertEqual(third_order_central_moment(pairwise_constrained_family(0)),
                         -Fraction(1,8))
        self.assertEqual(third_order_central_moment(pairwise_constrained_family(Fraction(1,8))),
                         Fraction(0))
        self.assertEqual(third_order_central_moment(pairwise_constrained_family(Fraction(1,4))),
                         Fraction(1,8))

    def test_even_and_odd_parity_identical_entropy_but_different_accessibility(self):
        even=pairwise_constrained_family(0)
        odd=pairwise_constrained_family(Fraction(1,4))
        self.assertAlmostEqual(joint_inventory_entropy_bits(even),2)
        self.assertAlmostEqual(joint_inventory_entropy_bits(odd),2)
        self.assertNotEqual(complete_retention_probability(even),
                            complete_retention_probability(odd))

    def test_genotype_expression_combo_upper_bound_not_empirical_flower_richness(self):
        self.assertEqual(expected_expressible_subsets_upper_bound(
            pairwise_constrained_family(0)),Fraction(13,4))
        self.assertEqual(expected_expressible_subsets_upper_bound(
            pairwise_constrained_family(Fraction(1,8))),Fraction(27,8))
        self.assertEqual(expected_expressible_subsets_upper_bound(
            pairwise_constrained_family(Fraction(1,4))),Fraction(7,2))
        self.assertEqual(expected_expressible_subsets_upper_bound(
            co_retained_distribution()),Fraction(9,2))

    def test_full_coupling_can_raise_or_lower_compound_accessibility(self):
        positive=co_retained_distribution()
        even=pairwise_constrained_family(0)
        odd=pairwise_constrained_family(Fraction(1,4))
        half=Fraction(1,2)
        self.assertEqual(complete_retention_probability(decorrelate(positive,half)),
                         Fraction(5,16))
        self.assertEqual(complete_retention_probability(decorrelate(even,half)),
                         Fraction(1,16))
        self.assertEqual(complete_retention_probability(decorrelate(odd,half)),
                         Fraction(3,16))

    def test_decorrelation_preserves_marginals_and_pairwise_for_parity_cases(self):
        for theta in (0,Fraction(1,4)):
            before=pairwise_constrained_family(theta)
            after=decorrelate(before,Fraction(1,2))
            self.assertEqual(marginals(before),marginals(after))
            self.assertEqual(pairwise_joint(before),pairwise_joint(after))
            self.assertNotEqual(complete_retention_probability(before),
                                complete_retention_probability(after))

    def test_all_white_visible_emission_would_not_disclose_resource_inventory(self):
        # No expression is allowed in a white plant here. All inventories
        # yield identical observable white flower even though they differ.
        p=pairwise_constrained_family(0)
        q=pairwise_constrained_family(Fraction(1,4))
        visible=lambda dist:{"white":sum(dist.values()),"pigmented":Fraction(0)}
        self.assertEqual(visible(p),visible(q))
        self.assertNotEqual(complete_retention_probability(p),
                            complete_retention_probability(q))

    def test_reconstruct_general_nonsymmetric_distribution(self):
        p={
            "000":Fraction(1,10),"001":Fraction(1,10),
            "010":Fraction(1,20),"011":Fraction(1,20),
            "100":Fraction(1,10),"101":Fraction(1,10),
            "110":Fraction(1,5),"111":Fraction(3,10)
        }
        p=validate_joint(p)
        lower,upper=feasible_triple_bounds(marginals(p),pairwise_joint(p))
        self.assertTrue(lower<=p["111"]<=upper)
        self.assertEqual(joint_from_moments(marginals(p),pairwise_joint(p),p["111"]),p)

    def test_impossible_joint_moments_rejected(self):
        with self.assertRaises(ValueError):
            feasible_triple_bounds(
                (Fraction(1,10),)*3,
                {(0,1):Fraction(1,5),(0,2):Fraction(1,5),(1,2):Fraction(1,5)}
            )
        with self.assertRaises(ValueError):
            pairwise_constrained_family(Fraction(3,10))
        with self.assertRaises(ValueError):
            joint_from_moments(P_HALF,P_PAIRS,Fraction(1,3))

    def test_incomplete_or_invalid_inventory_rejected(self):
        with self.assertRaises(ValueError):
            validate_joint({"111":Fraction(1)})
        with self.assertRaises(ValueError):
            complete_retention_probability(co_retained_distribution(),required="000")
        with self.assertRaises(ValueError):
            decorrelate(co_retained_distribution(),Fraction(3,2))

    def test_subset_requirement_can_differ_from_full(self):
        even=pairwise_constrained_family(0)
        self.assertEqual(complete_retention_probability(even,"110"),Fraction(1,4))
        self.assertEqual(complete_retention_probability(even,"111"),Fraction(0))


class GeographicOrderingTest(unittest.TestCase):
    def setUp(self):
        self.p=pairwise_constrained_family(Fraction(1,8))

    def test_equal_geographic_decay_rates_remove_order_effect(self):
        r=migration_order_counterfactual(
            self.p,[.02]*3,[.02]*3,10,n_founders=8)
        self.assertAlmostEqual(r["complete_resource_early"],r["complete_resource_late"])
        self.assertAlmostEqual(r["realized_early"],r["realized_late"])
        self.assertAlmostEqual(r["ratio_early_over_late"],1)

    def test_target_higher_degradation_gives_stronger_relative_contrast_for_more_modules(self):
        r=migration_order_counterfactual(
            self.p,[.02]*3,[.12]*3,10,
            n_founders=8,assembly_success=.7,establishment_success=.5)
        self.assertAlmostEqual(r["ratio_early_over_late"],math.exp(-3))
        self.assertLess(r["realized_early"],r["realized_late"])
        self.assertGreater(r["realized_late"],0)
        result=exact_report()["migration_example"]["complexity_grid"]
        self.assertAlmostEqual(result["1"]["relative_early_over_late"],math.exp(-1))
        self.assertAlmostEqual(result["2"]["relative_early_over_late"],math.exp(-2))
        self.assertAlmostEqual(result["3"]["relative_early_over_late"],math.exp(-3))

    def test_no_initial_complete_resources_means_no_regain_without_recombination(self):
        even=pairwise_constrained_family(0)
        r=migration_order_counterfactual(even,[.02]*3,[.12]*3,10)
        self.assertEqual(r["complete_resource_early"],0)
        self.assertEqual(r["complete_resource_late"],0)
        self.assertEqual(r["realized_early"],0)
        self.assertEqual(r["realized_late"],0)
        self.assertIsNone(r["ratio_early_over_late"])

    def test_zero_assembly_or_establishment_eliminates_realization(self):
        p=.125
        self.assertEqual(founders_at_least_one_complete(p,20,assembly_success=0),0)
        self.assertEqual(founders_at_least_one_complete(p,20,establishment_success=0),0)
        self.assertGreater(founders_at_least_one_complete(p,20),0)

    def test_genetic_loss_on_required_subset(self):
        survival=probability_required_modules_after_decay(
            self.p,[.02,.05,.5],10,required="110")
        self.assertAlmostEqual(survival,.25*math.exp(-.07*10))

    def test_initial_distribution_can_have_same_marginals_but_different_migration_outcomes(self):
        odd=pairwise_constrained_family(Fraction(1,4))
        even=pairwise_constrained_family(0)
        self.assertEqual(marginals(odd),marginals(even))
        self.assertEqual(pairwise_joint(odd),pairwise_joint(even))
        aa=migration_order_counterfactual(odd,[.02]*3,[.12]*3,10)
        bb=migration_order_counterfactual(even,[.02]*3,[.12]*3,10)
        self.assertGreater(aa["realized_early"],bb["realized_early"])

    def test_invalid_rate_or_founder_count_rejected(self):
        with self.assertRaises(ValueError):
            migration_order_counterfactual(self.p,[-.02]*3,[.12]*3,10)
        with self.assertRaises(ValueError):
            migration_order_counterfactual(self.p,[.02]*2,[.12]*3,10)
        with self.assertRaises(ValueError):
            founders_at_least_one_complete(.2,0)
        with self.assertRaises(ValueError):
            founders_at_least_one_complete(1.1,10)

    def test_report_explicitly_synthetic(self):
        r=exact_report()
        self.assertFalse(r["empirical_plant_data_admitted"])
        self.assertEqual(r["exact_sharp_bound_p111_given_marginals_half_and_pairs_quarter"],
                         ["0","1/4"])
        self.assertEqual(set(r["distributions"]),
                         {"even_parity","independent","odd_parity","positive_linkage"})


if __name__=="__main__":
    unittest.main()
