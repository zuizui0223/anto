"""Independent exhaustive oracles for ANTO v22 biological vs technical assays.

All prevalence, test specificity/sensitivity and shared handling-error
numbers are declared synthetic. No real Iochrominae functional assays.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from anto_individual_function_assay_hierarchy_v22 import (
    HypotheticalAssay, binomial_tail, individual_call_probabilities,
    one_taxon_screen, finite_budget_search, same_marginal_nonidentification,
    report, serialize
)


def raw_assay_outcome_oracle(r,k,theta,assay):
    """Enumerate latent competence, handling state and each binary test."""
    prob_call=F(0)
    prob_true_call=F(0)
    for competent in (0,1):
        z_weight=theta if competent else 1-theta
        for failure in (0,1):
            h_weight=(assay.shared_plant_handling_failure if failure else
                      1-assay.shared_plant_handling_failure)
            for ys in product((0,1),repeat=r):
                if failure:
                    # Shared handling mistake reverses all replicates
                    if any(y != (0 if competent else 1) for y in ys):
                        continue
                    cond=F(1)
                else:
                    p1=(assay.sensitivity if competent else 1-assay.specificity)
                    cond=F(1)
                    for y in ys:
                        cond*=p1 if y else (1-p1)
                wt=z_weight*h_weight*cond
                if sum(ys)>=k:
                    prob_call+=wt
                    if competent:
                        prob_true_call+=wt
    return prob_call,prob_true_call


class V22BiologicalVersusTechnical(unittest.TestCase):
    def setUp(self):
        self.i=HypotheticalAssay(F(9,10),F(19,20))
        self.c2=HypotheticalAssay(F(9,10),F(19,20),F(1,50))

    def test_binomial_tail_boundary_and_exact_negative_controls(self):
        self.assertEqual(binomial_tail(4,2,F(9,10)),F(9963,10000))
        self.assertEqual(binomial_tail(4,2,F(1,20)),F(11215,800000))
        self.assertEqual(binomial_tail(3,3,F(1,20)),F(1,8000))
        self.assertEqual(binomial_tail(3,1,F(0)),0)
        self.assertEqual(binomial_tail(3,1,F(1)),1)
        for wrong in (0,-1,True,1.1):
            with self.assertRaises(ValueError):
                binomial_tail(wrong,1,F(1,2))
        with self.assertRaisesRegex(ValueError,"exceeds"):
            binomial_tail(2,3,F(1,2))

    def test_explicit_latent_individual_and_test_enumeration_matches_formula(self):
        examined=0
        for assay in (self.i,self.c2):
            for theta in (F(0),F(1,5),F(1,2),F(1)):
                for r in (1,2,3):
                    for k in range(1,r+1):
                        pc,pt=raw_assay_outcome_oracle(r,k,theta,assay)
                        for n in (1,2,3):
                            target=one_taxon_screen(n,r,k,theta,assay)
                            self.assertEqual(
                                target["chance_at_least_one_plant_screen_positive"],
                                1-(1-pc)**n)
                            self.assertEqual(
                                target["chance_at_least_one_TRULY_competent_plant_correctly_called"],
                                1-(1-pt)**n)
                            examined+=1
        self.assertEqual(examined,144)

    def test_four_calls_are_not_four_biological_lineages(self):
        a=one_taxon_screen(4,1,1,F(1,5),self.i)
        b=one_taxon_screen(1,4,1,F(1,5),self.i)
        c=one_taxon_screen(1,4,2,F(1,5),self.i)
        self.assertEqual(a["total_assay_calls"],b["total_assay_calls"])
        self.assertEqual(a["total_assay_calls"],c["total_assay_calls"])
        self.assertEqual(a["false_taxon_screen_if_all_plants_incompetent"],
                         b["false_taxon_screen_if_all_plants_incompetent"])
        self.assertEqual(a["false_taxon_screen_if_all_plants_incompetent"],
                         F(29679,160000))
        self.assertEqual(a["chance_at_least_one_TRULY_competent_plant_correctly_called"],
                         F("0.54787824"))
        self.assertEqual(b["chance_at_least_one_TRULY_competent_plant_correctly_called"],
                         F("0.19998"))
        self.assertEqual(c["chance_at_least_one_TRULY_competent_plant_correctly_called"],
                         F("0.19926"))
        self.assertEqual(c["false_taxon_screen_if_all_plants_incompetent"],
                         F(11215,800000))
        self.assertLess(c["false_taxon_screen_if_all_plants_incompetent"],
                        a["false_taxon_screen_if_all_plants_incompetent"])

    def test_false_taxon_screens_not_diagnosis_of_species_function(self):
        no_function=one_taxon_screen(6,2,2,F(0),self.i)
        self.assertEqual(no_function["chance_at_least_one_TRULY_competent_plant_correctly_called"],0)
        self.assertEqual(no_function["chance_at_least_one_plant_screen_positive"],
                         no_function["false_taxon_screen_if_all_plants_incompetent"])
        all_function=one_taxon_screen(6,2,2,F(1),self.i)
        self.assertEqual(all_function["chance_zero_competent_plants_sampled"],0)
        self.assertGreater(all_function["chance_at_least_one_TRULY_competent_plant_correctly_called"],0)

    def test_common_handling_failure_cannot_be_repaired_by_technical_repeats(self):
        for r in (1,2,4,12):
            for k in (1,r):
                _,fp=individual_call_probabilities(self.c2,r,k)
                self.assertGreaterEqual(fp,F(1,50))
                tp,_=individual_call_probabilities(self.c2,r,k)
                self.assertLessEqual(tp,F(49,50))

    def test_one_raw_call_prevalence_is_confounded_with_unknown_assay_accuracy(self):
        d=same_marginal_nonidentification()
        self.assertEqual(d["same_observable_single_binary_call_probability"],F(11,50))
        a,b=d["hypothetical_world_A"],d["hypothetical_world_B"]
        self.assertEqual(a,(F(1,5),F(9,10),F(19,20)))
        self.assertEqual(b,(F(1,4),F(17,20),F(99,100)))
        self.assertNotEqual(a[0],b[0])
        self.assertTrue(d["not_observed_Iochrominae_assay_calibration"])

    def test_nominal_twelve_call_optimizer_gives_six_biological_plants(self):
        ans=finite_budget_search(12,F(1,5),self.i,F(1,20))
        chosen=ans["optimal_under_THIS_STIPULATED_likelihood_only"]
        self.assertEqual((chosen["n_biologically_distinct_plants"],
                          chosen["technical_repeats_per_plant"],
                          chosen["minimum_positive_repeats_for_individual_call"]), (6,2,2))
        self.assertEqual(chosen["total_assay_calls"],12)
        self.assertLessEqual(chosen["false_taxon_screen_if_all_plants_incompetent"],F(1,20))
        self.assertTrue(ans["not_real_sample_size_power_or_estimated_field_error"])

    def test_two_percent_shared_handling_changes_optimal_synthetic_allocation(self):
        result=finite_budget_search(12,F(1,5),self.c2,F(1,20))
        chosen=result["optimal_under_THIS_STIPULATED_likelihood_only"]
        self.assertEqual((chosen["n_biologically_distinct_plants"],
                          chosen["technical_repeats_per_plant"],
                          chosen["minimum_positive_repeats_for_individual_call"]), (2,6,3))
        self.assertEqual(chosen["total_assay_calls"],12)
        self.assertLessEqual(chosen["false_taxon_screen_if_all_plants_incompetent"],F(1,20))
        base=finite_budget_search(12,F(1,5),self.i,F(1,20))
        self.assertLess(
            chosen["chance_at_least_one_TRULY_competent_plant_correctly_called"],
            base["optimal_under_THIS_STIPULATED_likelihood_only"]
                ["chance_at_least_one_TRULY_competent_plant_correctly_called"])

    def test_five_percent_common_failure_eliminates_all_strict_five_percent_screen_designs(self):
        invalid=HypotheticalAssay(F(9,10),F(19,20),F(1,20))
        found=finite_budget_search(12,F(1,5),invalid,F(1,20))
        self.assertEqual(found["designs_passing_declared_false_screen_ceiling"],0)
        self.assertIsNone(found["optimal_under_THIS_STIPULATED_likelihood_only"])

    def test_small_budget_optimizer_matches_independent_enumeration(self):
        for assay in (self.i,self.c2):
            for ceiling in (F(1,20),F(1,4)):
                b=6
                expected=[]
                for n in range(1,b+1):
                    for r in range(1,b//n+1):
                        for k in range(1,r+1):
                            call,truth=raw_assay_outcome_oracle(r,k,F(1,5),assay)
                            fake=raw_assay_outcome_oracle(r,k,F(0),assay)[0]
                            false=1-(1-fake)**n
                            true=1-(1-truth)**n
                            if false<=ceiling:
                                expected.append((true,false,n*r,n,r,k))
                expected.sort(key=lambda x:(-x[0],x[1],x[2],x[3],x[4],x[5]))
                got=finite_budget_search(b,F(1,5),assay,ceiling)
                self.assertEqual(got["designs_passing_declared_false_screen_ceiling"],
                                 len(expected))
                if expected:
                    c=got["optimal_under_THIS_STIPULATED_likelihood_only"]
                    self.assertEqual(
                        (c["chance_at_least_one_TRULY_competent_plant_correctly_called"],
                         c["false_taxon_screen_if_all_plants_incompetent"],
                         c["total_assay_calls"],
                         c["n_biologically_distinct_plants"],
                         c["technical_repeats_per_plant"],
                         c["minimum_positive_repeats_for_individual_call"]),
                        expected[0])
                else:
                    self.assertIsNone(got["optimal_under_THIS_STIPULATED_likelihood_only"])

    def test_error_rate_or_technical_replicate_input_rejections(self):
        for wrong in (False,True,-F(1,10),F(11,10),"NaN"):
            with self.assertRaises(ValueError):
                HypotheticalAssay(wrong,F(19,20))
        for wrong in (True,-F(1,10),F(1),F(11,10),"NaN"):
            with self.assertRaises(ValueError):
                HypotheticalAssay(F(9,10),F(19,20),wrong)
        for wrong in (True,-1,0,1.4):
            with self.assertRaises(ValueError):
                one_taxon_screen(wrong,2,1,F(1,5),self.i)
        with self.assertRaises(ValueError):
            one_taxon_screen(3,2,3,F(1,5),self.i)
        with self.assertRaises(ValueError):
            one_taxon_screen(3,2,1,F(11,10),self.i)
        with self.assertRaises(ValueError):
            finite_budget_search(12,F(1,5),self.i,F(11,10))

    def test_true_rate_monotonic_and_false_positive_not_improved_by_more_plants(self):
        for theta in (F(0),F(1,20),F(1,5),F(1,2),F(1)):
            old=F(-1)
            for n in range(1,9):
                r=one_taxon_screen(n,2,2,theta,self.i)
                self.assertGreaterEqual(r["chance_at_least_one_TRULY_competent_plant_correctly_called"],old)
                old=r["chance_at_least_one_TRULY_competent_plant_correctly_called"]
                if n>1:
                    prev=one_taxon_screen(n-1,2,2,theta,self.i)
                    self.assertGreater(
                       r["false_taxon_screen_if_all_plants_incompetent"],
                       prev["false_taxon_screen_if_all_plants_incompetent"])

    def test_original_biological_evidence_grain_not_overridden(self):
        d=report()
        self.assertEqual(d["observed_original_source_taxa"],22)
        self.assertEqual((d["observed_original_source_unpigmented_taxa"],
                          d["observed_original_source_pigmented_taxa"]),(7,15))
        self.assertEqual(d["actual_qualified_complete_flower_function_assays"],0)
        self.assertTrue(d["source_DFR_CDS_annotation_is_not_a_function_assay"])
        self.assertTrue(d["input_rates_are_hypothetical_not_calibrated"])
        self.assertTrue(d["one_sampled_population_not_an_entire_species_or_ancestral_transition"])
        self.assertTrue(d["technical_repeat_independence_requires_external_assay_validation"])
        self.assertTrue(d["real_field_shared_batch_and_plant_relatedness_not_covered"])
        self.assertTrue(d["no_natural_floral_trait_regain_heritability_or_fitness_measurement"])
        self.assertEqual(d["Cirsium_organismal_empirical_owner"],"zuizui0223/aza3")

    def test_full_serialized_receipt_and_optional_source_lock(self):
        raw=serialize(report())
        json.dumps(raw,ensure_ascii=False)
        self.assertEqual(raw["mechanistic_single_call_nonidentifiability"]
                         ["same_observable_single_binary_call_probability"],"11/50")
        self.assertEqual(raw["observed_original_source_taxa"],22)
        frozen=ROOT/"results/anto_v22_individual_function_assay_20261026.json"
        if frozen.is_file():
            self.assertEqual(raw,json.loads(frozen.read_text(encoding="utf-8")))

if __name__=="__main__":
    unittest.main()
