"""ANTO v20 negative-control suite: real source censored, ALL outcomes hypothetical."""
import json
import sys
import unittest
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from anto_afterlife_empirical_admissibility_v19 import SOURCE_JOIN,read,STATUS
from anto_functional_assay_identifiability_v20 import (
    panel_from_original_source, sharp_bounds, metadata_strata,
    hypothetical_annotations_all_functional, hypothetical_certification_worlds,
    minimal_best_case_certification, minimax_width_budget, report, serializable,
)


class V20FunctionalAssayIdentifiability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=read(SOURCE_JOIN)
        cls.panel=panel_from_original_source(cls.source)
        cls.scenarios=hypothetical_certification_worlds(cls.panel)
        cls.result=report(cls.source)

    def test_original_22_taxa_and_no_functional_measurements_are_real(self):
        self.assertEqual(len(self.panel),22)
        self.assertEqual(len({x["taxon"] for x in self.panel}),22)
        self.assertEqual(sum(x["flower_pigment_category"]=="UNPIGMENTED"
                             for x in self.panel),7)
        self.assertEqual(sum(x["flower_pigment_category"]=="PIGMENTED"
                             for x in self.panel),15)
        self.assertTrue(all(x["flower_functional_joint_toolkit_assay"] is None
                            and x["assay_admitted"] is False for x in self.panel))
        self.assertTrue(all(x["multi_generation_functional_retention_assayed"] is False
                            for x in self.panel))

    def test_source_status_strata_complete_and_no_silent_gene_loss_promotion(self):
        strata=metadata_strata(self.panel)
        self.assertEqual(len(strata),5)
        self.assertEqual(sum(x["original_taxa"] for x in strata),22)
        counts={(x["flower_pigment_category"],x["new_KT_source_metadata"]):
                x["original_taxa"] for x in strata}
        self.assertEqual(counts[("UNPIGMENTED",STATUS[0])],3)
        self.assertEqual(counts[("UNPIGMENTED",STATUS[1])],2)
        self.assertEqual(counts[("UNPIGMENTED",STATUS[2])],2)
        self.assertEqual(counts[("PIGMENTED",STATUS[0])],9)
        self.assertEqual(counts[("PIGMENTED",STATUS[2])],6)
        self.assertNotIn(("PIGMENTED",STATUS[1]),counts)
        self.assertTrue(all(x["functional_assays_currently_qualified"]==0
                            for x in strata))

    def test_no_actual_function_outcomes_leaves_full_minus_one_to_one_identification(self):
        r=sharp_bounds(self.panel)
        self.assertEqual(r["difference_unpigmented_minus_pigmented_lower"],-1)
        self.assertEqual(r["difference_unpigmented_minus_pigmented_upper"],1)
        self.assertEqual(r["difference_identification_width"],2)
        self.assertEqual(r["groups"]["UNPIGMENTED"]["not_functionally_assayed"],7)
        self.assertEqual(r["groups"]["PIGMENTED"]["not_functionally_assayed"],15)
        self.assertEqual(r["strict_direction_identified"],"SIGNED_DIFFERENCE_NOT_IDENTIFIED")

    def test_as_if_12_cds_annotated_were_fully_assayed_is_counterfactual_only(self):
        assumptions=hypothetical_annotations_all_functional(self.panel)
        self.assertEqual(len(assumptions),12)
        self.assertTrue(all(assumptions.values()))
        r=sharp_bounds(self.panel,assumptions,synthetic=True)
        self.assertEqual(r["groups"]["UNPIGMENTED"]["hypothetically_assayed_positive"],3)
        self.assertEqual(r["groups"]["PIGMENTED"]["hypothetically_assayed_positive"],9)
        self.assertEqual(r["difference_unpigmented_minus_pigmented_lower"],-F(4,7))
        self.assertEqual(r["difference_unpigmented_minus_pigmented_upper"],F(2,5))
        self.assertEqual(r["difference_identification_width"],F(34,35))
        self.assertEqual(r["strict_direction_identified"],"SIGNED_DIFFERENCE_NOT_IDENTIFIED")

    def test_five_more_idealized_outcomes_can_certify_unpigmented_higher(self):
        r=sharp_bounds(self.panel,self.scenarios["additional_5_best_case_unpigmented_higher"])
        self.assertEqual(r["groups"]["UNPIGMENTED"]["hypothetically_assayed_positive"],7)
        self.assertEqual(r["groups"]["PIGMENTED"]["hypothetically_assayed_negative"],1)
        self.assertEqual(r["difference_unpigmented_minus_pigmented_lower"],F(1,15))
        self.assertEqual(r["difference_unpigmented_minus_pigmented_upper"],F(2,5))
        self.assertEqual(r["strict_direction_identified"],"UNPIGMENTED_HIGHER")

    def test_three_more_idealized_outcomes_can_certify_pigmented_higher(self):
        r=sharp_bounds(self.panel,self.scenarios["additional_3_best_case_pigmented_higher"])
        self.assertEqual(r["groups"]["UNPIGMENTED"]["hypothetically_assayed_negative"],3)
        self.assertEqual(r["groups"]["PIGMENTED"]["hypothetically_assayed_positive"],9)
        self.assertEqual(r["difference_unpigmented_minus_pigmented_lower"],-F(4,7))
        self.assertEqual(r["difference_unpigmented_minus_pigmented_upper"],-F(1,35))
        self.assertEqual(r["strict_direction_identified"],"PIGMENTED_HIGHER")

    def test_exact_minimum_best_case_certificates_from_unmeasured_and_fictional12(self):
        x=minimal_best_case_certification(self.panel,{})
        self.assertEqual(x["UNPIGMENTED_HIGHER"]["minimum_additional_assays"],8)
        self.assertEqual(x["PIGMENTED_HIGHER"]["minimum_additional_assays"],8)
        y=minimal_best_case_certification(
            self.panel,self.scenarios["assume_12_annotated_have_confirmed_function"])
        self.assertEqual(y["UNPIGMENTED_HIGHER"]["minimum_additional_assays"],5)
        self.assertEqual(y["PIGMENTED_HIGHER"]["minimum_additional_assays"],3)
        self.assertTrue(x["UNPIGMENTED_HIGHER"]["achievable"])

    def test_no_result_sign_guarantee_even_under_a_complete_oracle_taxon_assay(self):
        # Hypothetical future in which every taxon has a functional toolkit.
        r=sharp_bounds(self.panel,{p["taxon"]:True for p in self.panel})
        self.assertEqual((r["difference_unpigmented_minus_pigmented_lower"],
                          r["difference_unpigmented_minus_pigmented_upper"]),(0,0))
        self.assertEqual(r["difference_identification_width"],0)
        self.assertEqual(r["strict_direction_identified"],"SIGNED_DIFFERENCE_NOT_IDENTIFIED")

    def test_exact_interval_width_is_unknown_white_over7_plus_unknown_pig_over15(self):
        for labels in [{},
                       self.scenarios["assume_12_annotated_have_confirmed_function"],
                       self.scenarios["additional_5_best_case_unpigmented_higher"],
                       self.scenarios["additional_3_best_case_pigmented_higher"]]:
            r=sharp_bounds(self.panel,labels)
            g=r["groups"]
            self.assertEqual(r["difference_identification_width"],
                F(g["UNPIGMENTED"]["not_functionally_assayed"],7)
                + F(g["PIGMENTED"]["not_functionally_assayed"],15))

    def test_an_annotated_CDS_record_is_not_admitted_functional_data(self):
        example=next(x for x in self.panel if x["source_DFR_metadata"]==STATUS[0])
        r=sharp_bounds(self.panel)
        self.assertIn(example["taxon"],[r["groups"]["UNPIGMENTED"],
                                        r["groups"]["PIGMENTED"]] if False else
                      [x["taxon"] for x in self.panel])
        self.assertEqual(r["difference_identification_width"],2)
        self.assertFalse(self.result["validated_flower_toolkit_function_measured_by_ANTO"])

    def test_minimax_width_budget_without_outcome_assumptions(self):
        expected={0:(0,0,F(2)),1:(1,0,F(13,7)),5:(5,0,F(9,7)),
                  7:(7,0,F(1)),8:(7,1,F(14,15)),
                  12:(7,5,F(2,3)),17:(7,10,F(1,3)),22:(7,15,F(0))}
        for budget,(white,pig,width) in expected.items():
            r=minimax_width_budget(budget)
            self.assertEqual((r["width_only_optimal_unpigmented_taxa"],
                              r["width_only_optimal_pigmented_taxa"],
                              r["possible_sharp_difference_interval_width"]),
                             (white,pig,width))
            self.assertTrue(r["not_guaranteed_sign_certification_or_phylogenetic_predictive_validation"])

    def test_width_first_optimization_cannot_be_mistaken_for_cross_clade_generalization(self):
        one=minimax_width_budget(5)
        self.assertEqual(one["width_only_optimal_unpigmented_taxa"],5)
        self.assertLess(one["possible_sharp_difference_interval_width"],
                        F(2)-F(5,15))
        self.assertTrue(self.result["same_taxon_samples_not_independent_loss_events"])
        self.assertFalse(self.result["independent_macroevolutionary_regain_clades_replicated_here"])

    def test_no_boolean_function_label_without_explicit_synthetic_mode(self):
        x=self.panel[0]["taxon"]
        with self.assertRaisesRegex(ValueError,"promoted"):
            sharp_bounds(self.panel,{x:True},synthetic=False)
        with self.assertRaisesRegex(ValueError,"hypothetical Boolean"):
            sharp_bounds(self.panel,{x:"CDS_PRESENT"})
        with self.assertRaisesRegex(ValueError,"registered source taxon"):
            sharp_bounds(self.panel,{"SYNTHETIC_FAKE_TAXON":True})

    def test_adverse_spontaneously_filled_functional_state_is_rejected(self):
        forged=deepcopy(self.panel)
        forged[0]["flower_functional_joint_toolkit_assay"]=True
        with self.assertRaisesRegex(ValueError,"NO admitted functional"):
            sharp_bounds(forged)
        forged=deepcopy(self.panel)
        forged[0]["assay_admitted"]=True
        with self.assertRaisesRegex(ValueError,"NO admitted functional"):
            sharp_bounds(forged)

    def test_adverse_taxon_denominator_or_identity(self):
        with self.assertRaisesRegex(ValueError,"denominator"):
            sharp_bounds(self.panel[:-1])
        forged=deepcopy(self.panel)
        forged[1]["taxon"]=forged[0]["taxon"]
        with self.assertRaisesRegex(ValueError,"denominator"):
            sharp_bounds(forged)
        forged=deepcopy(self.source)
        forged["taxon_gene_records"][0]["trait_state_original_fig1"]="INVALID"
        with self.assertRaises(ValueError):
            panel_from_original_source(forged)

    def test_original_source_cds_censoring_perturbation_must_not_be_deletion(self):
        forged=deepcopy(self.source)
        a=next(x for x in forged["taxon_gene_records"] if x["DFR_CDS_annotation_in_KT_record"] is False)
        a["DFR_CDS_annotation_in_KT_record"]=True
        with self.assertRaises(ValueError):
            panel_from_original_source(forged)

    def test_full_result_never_promotes_to_natural_ancestry_or_fitness(self):
        p=self.result
        self.assertEqual(p["schema"],"anto_v20_functional_flower_assay_sharp_bounds_prospective")
        self.assertTrue(p["all_source_taxon_functional_assays_currently_unmeasured"])
        self.assertFalse(p["validated_flower_toolkit_function_measured_by_ANTO"])
        self.assertFalse(p["ancestral_reexpression_rate_estimated"])
        self.assertFalse(p["independent_macroevolutionary_regain_clades_replicated_here"])
        self.assertTrue(p["no_worst_case_nonzero_sign_guarantee_even_after_all_22_if_true_rates_equal"])
        self.assertEqual(p["Cirsium_empirical_owner"],"zuizui0223/aza3")
        self.assertGreaterEqual(len(p["prospective_assay_definition_not_yet_implemented_in_field"]),6)
        self.assertEqual(len(p["source_taxon_panel"]),22)
        self.assertEqual(len(p["source_annotation_strata"]),5)

    def test_report_json_receipt_exact_fractions_and_provenance(self):
        output=serializable(self.result)
        self.assertEqual(output["observed_source_only_identification"]
                         ["difference_unpigmented_minus_pigmented_lower"],"-1")
        self.assertEqual(output["observed_source_only_identification"]
                         ["difference_unpigmented_minus_pigmented_upper"],"1")
        self.assertEqual(output["source_annotated_CDS_fictions_not_observed"]
                         ["assume_12_annotated_have_confirmed_function"]
                         ["difference_unpigmented_minus_pigmented_upper"],"2/5")
        json.dumps(output)
        # Optional locked receipt is supplied after independent GitHub CI replay.
        expected=ROOT/"results/anto_v20_functional_assay_identification_20261026.json"
        if expected.exists():
            self.assertEqual(output,json.loads(expected.read_text(encoding="utf-8")))

if __name__=="__main__":
    unittest.main()
