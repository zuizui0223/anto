"""Source-grain and numerical adversarial gates for posthoc v47 reference."""
from __future__ import annotations

import copy
import itertools
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_Phlox_three_pop_conditional_v47 as v


class Phlox2024PosthocThreeSourcePopulationAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = v.read_registry()
        cls.parent = v.read_parent(cls.registry)
        cls.output = v.analyze(cls.parent, cls.registry)

    def test_sha1_git_blob_is_real_git_format(self):
        self.assertEqual(v.git_blob_sha1(b""), "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391")
        self.assertEqual(v.git_blob_sha1(v.PARENT.read_bytes()),
                         self.registry["v46_parent_git_blob_sha1"])

    def test_original_source_registry_grain_and_evolutionary_nonclaims(self):
        self.assertEqual(self.registry["source_2024_original_populations"], 20)
        self.assertEqual(self.registry["source_2024_taxon_populations"],
                         {"amoena": 8, "deamii": 3, "pilosa": 9})
        self.assertEqual(self.registry["new_2014_Phlox_selection_refits"], 0)
        self.assertEqual(self.registry["independent_new_biological_studies"], 0)
        self.assertFalse(self.registry["has_p_value"])

    def test_original_v46_parent_continues_to_fail_primary(self):
        self.assertEqual(self.parent["primary_fruit_transfer_pass_models"], [])
        self.assertEqual(self.output["primary_v46_all_taxa_transfer_pass_models"], [])
        for m in v.MODELS:
            self.assertEqual(
                self.parent["endpoints"]["fruits_2019"]["model_results"][m]
                    ["predeclared_primary_fruit_decision"],
                "NOT_ESTABLISHED")

    def test_exact_56_84_140_680_combinations(self):
        self.assertEqual(len(list(itertools.combinations(range(8), 3))), 56)
        self.assertEqual(len(list(itertools.combinations(range(9), 3))), 84)
        self.assertEqual(len(list(itertools.combinations(range(17), 3))), 680)
        self.assertEqual(self.output["n_same_taxon_triplets_total"], 140)

    def test_no_new_models_or_original_field_observations(self):
        self.assertEqual(self.output["new_prediction_fits"], 0)
        self.assertEqual(self.output["original_source_studies_added_in_v47"], 0)
        self.assertEqual(self.output["new_original_field_populations"], 0)
        self.assertFalse(self.output["new_hypothesis_test_or_p_values"])

    def test_deamii_primary_common_and_garden_are_more_negative_than_all_140(self):
        obj = self.output["endpoints"]["fruits_2019"]["models"]
        for m in ("G", "E"):
            self.assertEqual(
                obj[m]["primary_same_taxon_choose3_140_reference"]
                    ["number_at_or_below_deamii"], 0)
            self.assertLess(obj[m]["deamii_3_pop_relative_MSE_gain"], -0.09)

    def test_deamii_taxon_specific_models_are_not_universally_extreme(self):
        obj = self.output["endpoints"]["fruits_2019"]["models"]
        self.assertEqual(
            obj["T"]["primary_same_taxon_choose3_140_reference"]
                ["number_at_or_below_deamii"], 8)
        self.assertEqual(
            obj["TE"]["primary_same_taxon_choose3_140_reference"]
                ["number_at_or_below_deamii"], 11)

    def test_deamii_common_slope_fruit_all_three_negative(self):
        d = self.output["endpoints"]["fruits_2019"]["models"]["G"]
        self.assertEqual(d["deamii_positive_original_population_folds"], 0)
        self.assertEqual(d["deamii_3_original_population_ids"],
                         ["1719", "1720", "1721"])
        self.assertTrue(all(q["remaining_relative_MSE_gain"] < 0
                            for q in d["deamii_two_pop_delete_one_arithmetic"]))

    def test_original_flower_secondary_refuses_primary_promotion(self):
        d = self.output["endpoints"]["flowers_2019"]
        self.assertEqual(d["outcome_role"], "SECONDARY_DESCRIPTIVE_REAGGREGATION")
        self.assertEqual(d["models"]["G"]
                         ["primary_same_taxon_choose3_140_reference"]
                         ["number_at_or_below_deamii"], 0)

    def test_biomass_negative_source_gating_untouched(self):
        self.assertTrue(self.output["technical_biomass_nonadmission_retained"])
        self.assertNotIn("biomass_2019", self.output["endpoints"])
        self.assertEqual(self.parent["source_original_grain"]
                         ["biomass_original_negative_n"], 13)

    def test_exact_new_frozen_receipt_after_reaggregation(self):
        self.assertEqual(self.output, json.loads(v.RESULT.read_text()))

    def test_fake_original_study_claim_rejected(self):
        wrong = copy.deepcopy(self.parent)
        wrong["independent_evolutionary_radiations_added"] = 10
        with self.assertRaisesRegex(v.V47Error, "SOURCE_OR_BIOLOGICAL_SCOPE"):
            v.validate_parent(wrong, self.registry)

    def test_forged_20_pop_outcome_rejected(self):
        wrong = copy.deepcopy(self.parent)
        wrong["endpoints"]["fruits_2019"]["per_population_holdout"][1]["population"] = "1719"
        with self.assertRaisesRegex(v.V47Error, "ORIGINAL_POPULATION_GRAIN"):
            v.validate_parent(wrong, self.registry)

    def test_forged_outcome_drift_rejected(self):
        wrong = copy.deepcopy(self.parent)
        wrong["endpoints"]["fruits_2019"]["per_population_holdout"][0]["mse_by_model"]["G"] += 1.0
        with self.assertRaisesRegex(v.V47Error, "LOSS_AGGREGATE_DRIFT"):
            v.validate_parent(wrong, self.registry)

    def test_nonfinite_or_nonpositive_fake_squared_error_fails(self):
        for bad in (0, float("nan"), float("inf"), -1):
            wrong = copy.deepcopy(self.parent)
            wrong["endpoints"]["flowers_2019"]["per_population_holdout"][0]["mse_by_model"]["TE"] = bad
            with self.assertRaisesRegex(v.V47Error, "INVALID_HOLDOUT_MSE"):
                v.validate_parent(wrong, self.registry)

    def test_rounded_or_mixed_reference_not_valid_randomization_pvalue(self):
        for e in v.ENDPOINTS:
            for m in v.MODELS:
                r = self.output["endpoints"][e]["models"][m]
                for key in ("amoena_same_taxon_choose3_reference",
                            "pilosa_same_taxon_choose3_reference",
                            "primary_same_taxon_choose3_140_reference",
                            "secondary_mixed_taxon_choose3_680_reference"):
                    self.assertFalse(r[key]["is_valid_randomization_p_value"])
                    self.assertFalse(r[key]["is_independent_study_replication"])

    def test_combinatorial_reference_is_not_independent_taxa(self):
        for e in v.ENDPOINTS:
            for m in v.MODELS:
                r = self.output["endpoints"][e]["models"][m]
                self.assertEqual(
                    r["primary_same_taxon_choose3_140_reference"]
                        ["n_overlapping_reference_combinations"], 140)
                self.assertEqual(
                    r["secondary_mixed_taxon_choose3_680_reference"]
                        ["n_overlapping_reference_combinations"], 680)


if __name__ == "__main__":
    unittest.main()
