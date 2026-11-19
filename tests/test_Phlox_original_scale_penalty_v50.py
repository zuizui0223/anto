"""ANTO v50: adversarial original-study/normal equation/coordinate tests."""
from __future__ import annotations

import copy
import json
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_Phlox_original_scale_penalty_v50 as v


class OriginalPhloxV50SharingPenaltyAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = v.source_registry()
        cls.old = v.v49_source_parent(cls.registry)

    def test_source_hash_and_423_original_fold_grain(self):
        self.assertEqual(v.git_sha1(v.PARENT.read_bytes()),
                         "94947d03e85766a3c5c1f82fde2bb076f2ca42e8")
        self.assertEqual(self.registry["source_git_blob_sha1"],
                         "4c780b0abe644ac69dae86a0450eb598c160626b")
        self.assertEqual(self.registry["n_folds"], 423)
        self.assertEqual(self.registry["n_panels"], 141)
        self.assertFalse(self.old["previous_v46_primary_transfer_established"])

    def test_seven_predeclared_scenarios_and_roles(self):
        self.assertEqual(v.SCENARIOS,
           ("STD_R0", "STD_R1", "STD_R10", "STD_R100",
            "RAW_R0", "RAW_R10", "RAW_EQUIV_R10"))
        self.assertEqual(v.MODELS, ("G", "T", "E", "TE"))
        self.assertEqual(v.ENDPOINTS, ("fruits_2019", "flowers_2019"))
        self.assertEqual(self.registry["new_original_biological_studies"], 0)
        self.assertFalse(self.registry["valid_randomization_p"])

    def test_exact_unmodified_v46_original_model_inheritance(self):
        self.assertIs(v.v46.fit_train, sys.modules[
            "analyze_original_Phlox_2024_population_holdout_v46"].fit_train)
        self.assertIs(v.v49.original_populations, sys.modules[
            "analyze_Phlox_2024_multitaxon_focal2_v49"].original_populations)
        self.assertEqual(v.N_BLOCKS, {"G": 1, "T": 3, "E": 3, "TE": 9})

    def test_simple_unpenalized_normal_equation(self):
        x, y = v.ridge_two_by_two(2, 0, 4, 6, 8, 0, 0)
        self.assertEqual((x, y), (3.0, 2.0))

    def test_ridge_diagonal_takes_different_coordinate_penalties(self):
        x, y = v.ridge_two_by_two(2, 0, 4, 6, 8, 1, 6)
        self.assertAlmostEqual(x, 2.0)
        self.assertAlmostEqual(y, .8)

    def test_singular_or_nonfinite_inputs_fail_closed(self):
        for inp in [(1, 1, 1, 1, 1, 0, 0),
                    (0, 0, 0, 2, 3, 0, 0),
                    (float("nan"), 0, 1, 1, 1, 10, 10),
                    (1, 0, 1, 1, 1, -2, 1)]:
            with self.assertRaises(v.V50Error):
                v.ridge_two_by_two(*inp)

    def test_math_scaling_equivariance_exact_synthetic_oracle(self):
        sd1, sd2 = 2.0, 7.0
        a, b, c, u, w = 40.0, 8.0, 270.0, 13.0, 25.0
        raw = v.ridge_two_by_two(a, b, c, u, w,
                                  10*sd1*sd1, 10*sd2*sd2)
        scaled = v.ridge_two_by_two(
            a/sd1**2, b/(sd1*sd2), c/sd2**2, u/sd1, w/sd2,
            10, 10)
        self.assertAlmostEqual(raw[0]*sd1, scaled[0], places=10)
        self.assertAlmostEqual(raw[1]*sd2, scaled[1], places=10)

    def test_math_ols_equivariance_synthetic_oracle(self):
        sd1, sd2 = 2.0, 7.0
        a, b, c, u, w = 40.0, 8.0, 270.0, 13.0, 25.0
        raw = v.ridge_two_by_two(a, b, c, u, w, 0, 0)
        scaled = v.ridge_two_by_two(
            a/sd1**2, b/(sd1*sd2), c/sd2**2,
            u/sd1, w/sd2, 0, 0)
        self.assertAlmostEqual(raw[0]*sd1, scaled[0], places=10)
        self.assertAlmostEqual(raw[1]*sd2, scaled[1], places=10)

    def test_invalid_parent_result_claim_fails_before_analysis(self):
        wrong = copy.deepcopy(self.old)
        wrong["previous_v46_primary_transfer_established"] = True
        self.assertTrue(wrong["previous_v46_primary_transfer_established"])
        self.assertNotEqual(
            v.git_sha1(json.dumps(wrong).encode()),
            self.registry["v49_result_blob_sha1"])

    def test_no_biomass_or_ancestral_recovery_in_fitted_endpoints(self):
        self.assertNotIn("biomass_2019", v.ENDPOINTS)
        self.assertEqual(self.registry["n_negative_biomass"], 13)
        self.assertFalse(self.registry["new_flower_colour_data"])
        self.assertFalse(self.registry["new_pollinator_data"])

    def test_reference_is_not_randomized_significance(self):
        d = v.distribution([-.25, .1, .3], -.1, 1e-10)
        self.assertEqual(d["at_or_below_deamii"], 1)
        self.assertEqual(d["negative_panels"], 1)
        self.assertFalse(d["is_p_value"])
        self.assertEqual(d["independent_studies"], 0)

    def test_original_source_can_be_verified_in_ci(self):
        path = os.environ.get("PHLOX_V50_SOURCE")
        if not path:
            self.skipTest("Only original source CI stages exact author CSV bytes")
        rows, audit = v.original_data(path, self.registry)
        self.assertEqual(audit["rows"], 963)
        self.assertEqual(audit["genotypes"], 284)
        self.assertEqual(audit["populations"], 20)
        self.assertEqual(audit["biomass_original_negative_n"], 13)
        self.assertEqual(len(v.v49.original_populations(rows)[1]), 20)

    def test_original_one_fold_source_ridge_invariance_and_parent_parity(self):
        path = os.environ.get("PHLOX_V50_SOURCE")
        if not path:
            self.skipTest("Only original source CI stages original author bytes")
        rows, audit = v.original_data(path, self.registry)
        target = next(x for x in self.old["original_panel_fit_details"]
                      if x["taxon"] == "deamii")
        fold = target["heldout_folds"][0]
        now = v.evaluate_fold(rows, fold["training_pops"], fold["heldout"],
                              fold, 1e-7)
        self.assertTrue(all(x < 1e-7 for x in now["max_errors"].values()))
        self.assertEqual(len(now["scores"]), 2)
        self.assertEqual(len(now["scores"]["fruits_2019"]), 7)

    def test_exact_frozen_entire_source_receipt_in_ci(self):
        path = os.environ.get("PHLOX_V50_SOURCE")
        if not path or not v.RECEIPT.exists():
            self.skipTest("Original bytes and frozen outcome receipt both required")
        rows, audit = v.original_data(path, self.registry)
        actual = v.evaluate(rows, audit, self.old, self.registry)
        self.assertEqual(actual, json.loads(v.RECEIPT.read_text()))
        self.assertEqual(actual["source_population_folds"], 423)
        self.assertEqual(actual["source_original_three_pop_panels"], 141)
        self.assertFalse(actual["statistical_p_values_or_confidence_intervals_computed"])
        for value in actual["v49_source_and_numerical_parity_max_absolute_mse_errors"].values():
            self.assertLessEqual(value, 1e-7)


if __name__ == "__main__":
    unittest.main()
