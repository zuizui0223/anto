"""ANTO v49: adversarial multi-taxon training support / source-grain tests."""
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
import analyze_Phlox_2024_multitaxon_focal2_v49 as v


class V49MultiTaxonOriginal2024Holdouts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = v.registry()
        cls.old = v.frozen_v46(cls.r)

    def test_locked_source_2024_and_v46_parent_git_blob(self):
        self.assertEqual(v.git_blob_sha1(v.PARENT.read_bytes()),
                         "f03b0eb0177f807317ee1f8a7c0e65ee408d903b")
        self.assertEqual(self.r["source_CSV_git_blob_sha1"],
                         "4c780b0abe644ac69dae86a0450eb598c160626b")
        self.assertEqual(self.r["n_original_ramets"], 963)
        self.assertEqual(self.r["n_original_clonal_genotypes"], 284)

    def test_registry_does_not_claim_new_wild_flower_or_pollinators(self):
        self.assertEqual(self.r["n_additional_original_studies"], 0)
        self.assertFalse(self.r["flower_pigment_measured"])
        self.assertFalse(self.r["pollinator_measured"])
        self.assertFalse(self.r["genetic_fitness_measured"])
        self.assertTrue(self.r["is_posthoc"])
        self.assertFalse(self.r["valid_randomization_p"])

    def test_original_v46_primary_failure_retained(self):
        self.assertEqual(self.old["primary_fruit_transfer_pass_models"], [])
        for m in v.MODELS[1:]:
            self.assertEqual(
                self.old["endpoints"]["fruits_2019"]["model_results"][m]
                    ["predeclared_primary_fruit_decision"],
                "NOT_ESTABLISHED")

    def test_multitaxon_model_implementation_is_exact_v46(self):
        self.assertEqual(v.MODELS, ("B0", "G", "T", "E", "TE"))
        self.assertIs(v.v46.fit_train, sys.modules[
            "analyze_original_Phlox_2024_population_holdout_v46"].fit_train)
        self.assertIs(v.v46.prediction, sys.modules[
            "analyze_original_Phlox_2024_population_holdout_v46"].prediction)
        self.assertEqual(self.r["original_v46_ridge_penalty"], 10)

    def test_training_policy_two_focal_plus_all_nonfocal(self):
        p = {"a1": "amoena", "a2": "amoena", "a3": "amoena",
             "d1": "deamii", "d2": "deamii", "d3": "deamii",
             "p1": "pilosa", "p2": "pilosa"}
        actual = v.training_populations("amoena", ("a1", "a2", "a3"), "a1", p)
        self.assertEqual(set(actual),
                         {"a2", "a3", "d1", "d2", "d3", "p1", "p2"})
        self.assertNotIn("a1", actual)

    def test_reject_heldout_not_inside_focal_panel(self):
        p = {"a1": "amoena", "a2": "amoena", "a3": "amoena",
             "d1": "deamii", "d2": "deamii"}
        with self.assertRaisesRegex(v.V49Error, "INVALID_FOCAL_TRIPLET"):
            v.training_populations("amoena", ("a1", "a2", "a3"), "d1", p)

    def test_reject_cross_taxon_triple(self):
        p = {"a1": "amoena", "a2": "amoena", "d1": "deamii"}
        with self.assertRaisesRegex(v.V49Error, "INVALID_FOCAL_TRIPLET"):
            v.training_populations("amoena", ("a1", "a2", "d1"), "a1", p)

    def test_reject_repeated_population_in_triple(self):
        p = {"a1": "amoena", "a2": "amoena"}
        with self.assertRaisesRegex(v.V49Error, "INVALID_FOCAL_TRIPLET"):
            v.training_populations("amoena", ("a1", "a2", "a2"), "a1", p)

    def test_deamii_v46_original_anchor_rejects_fake_loss(self):
        anchor = self.old["endpoints"]["fruits_2019"]["per_population_holdout"]
        x = next(z for z in anchor if z["population"] == "1719")
        y = next(z for z in
                 self.old["endpoints"]["flowers_2019"]["per_population_holdout"]
                 if z["population"] == "1719")
        fold = {"heldout": "1719", "training_pops": [str(i) for i in range(19)],
                "two_focal_training_populations": ["1720", "1721"],
                "n_test_ramets": x["n_ramets"],
                "loss": {"fruits_2019": copy.deepcopy(x["mse_by_model"]),
                         "flowers_2019": copy.deepcopy(y["mse_by_model"])}}
        v.check_deamii_against_v46(fold, self.old, 1e-7)
        fold["loss"]["fruits_2019"]["G"] += 0.001
        with self.assertRaisesRegex(v.V49Error,
                                    "SAME_MODEL_DEAMII_ANCHOR_DIFFERENT"):
            v.check_deamii_against_v46(fold, self.old, 1e-7)

    def test_parent_source_claim_change_rejected(self):
        wrong = copy.deepcopy(self.old)
        wrong["primary_fruit_transfer_pass_models"] = ["G"]
        self.assertEqual(wrong["primary_fruit_transfer_pass_models"], ["G"])
        # Parent binding SHA is checked before the parent can be read by the analyzer.
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "fake.json"
            p.write_text(json.dumps(wrong))
            self.assertNotEqual(v.git_blob_sha1(p.read_bytes()),
                                self.r["v46_derived_result_git_blob_sha1"])

    def test_count_exact_141_panels_423_folds(self):
        self.assertEqual(self.r["panels_by_taxon"],
                         {"amoena": 56, "deamii": 1, "pilosa": 84})
        self.assertEqual(self.r["n_source_triplet_panels"], 141)
        self.assertEqual(self.r["n_holdout_folds"], 423)

    def test_reference_is_not_randomization_significance(self):
        r = v.describe_reference([-.2, -.03, .05, .2], -.05, 1e-10)
        self.assertEqual(r["n_at_or_below_deamii"], 1)
        self.assertEqual(r["n_negative"], 2)
        self.assertFalse(r["randomization_p_value"])
        self.assertEqual(r["independent_field_studies"], 0)

    def test_reference_inadmissible_zero_or_nonfinite(self):
        with self.assertRaisesRegex(v.V49Error, "NONFINITE_SOURCE_OUTPUT"):
            v.num(float("nan"))
        with self.assertRaisesRegex(v.V49Error, "NONFINITE_SOURCE_OUTPUT"):
            v.num(float("inf"))

    def test_panel_gain_preserves_original_b0_denominator(self):
        folds = [{"loss": {"fruits_2019": {"B0": 2.0, "G": 1.0}}},
                 {"loss": {"fruits_2019": {"B0": 4.0, "G": 3.0}}},
                 {"loss": {"fruits_2019": {"B0": 4.0, "G": 2.0}}}]
        self.assertAlmostEqual(v.panel_gain(folds, "fruits_2019", "G"), 0.4)
        with self.assertRaisesRegex(v.V49Error, "NOT_THREE_FOLD_PANEL"):
            v.panel_gain(folds[:2], "fruits_2019", "G")

    def test_fruit_and_flower_are_separate_original_endpoints(self):
        self.assertEqual(v.ENDPOINTS, ("fruits_2019", "flowers_2019"))
        self.assertEqual(self.old["nonadmitted_secondary_biomass"]["status"],
                         "BLOCKED_PREDECLARED_LOG1P_DOMAIN")

    def test_replay_original_source_if_staged_in_ci(self):
        source = os.environ.get("PHLOX_V49_SOURCE")
        if not source:
            self.skipTest("Original source CSV only staged in dedicated CI")
        rows, audit = v.author_source(source, self.r)
        bypop, poptax, taxon = v.original_populations(rows)
        self.assertEqual(len(poptax), 20)
        self.assertEqual({k: len(vv) for k, vv in taxon.items()},
                         {"amoena": 8, "deamii": 3, "pilosa": 9})
        self.assertEqual(audit["genotypes"], 284)
        self.assertEqual(audit["biomass_original_negative_n"], 13)
        for target in ("1719", "1720", "1721"):
            f = v.score_one_fold(rows, bypop, poptax, "deamii",
                                 ("1719", "1720", "1721"), target)
            v.check_deamii_against_v46(f, self.old, 1e-7)

    def test_full_423_source_fits_receipt_if_available(self):
        source = os.environ.get("PHLOX_V49_SOURCE")
        if not source or not v.RECEIPT.exists():
            self.skipTest("Frozen v49 receipt and exact source both required")
        rows, audit = v.author_source(source, self.r)
        output = v.summarize(rows, audit, self.old, self.r)
        self.assertEqual(output, json.loads(v.RECEIPT.read_text()))
        self.assertEqual(output["source_population_model_fits"], 423)
        self.assertTrue(output["all_deamii_per_population_loss_reproduced_from_v46"])
        self.assertFalse(output["new_p_values_or_confidence_intervals"])


    def test_posthoc_v48_v49_baseline_parity_846_folds(self):
        """After-result diagnostic: same B0, but slopes use other taxa in v49."""
        earlier = ROOT / "results/anto_v48_phlox2024_train2_matched_source_20261119.json"
        newer = v.RECEIPT
        if not newer.exists():
            self.skipTest("v49 receipt not frozen for baseline-parity diagnostic")
        self.assertEqual(
            v.git_blob_sha1(earlier.read_bytes()),
            "e07fb43c92a6ef3fd5c790720daab8047f937803")
        old = json.loads(earlier.read_text())
        now = json.loads(newer.read_text())
        self.assertEqual(old["new_source_model_fits"], 423)
        self.assertEqual(now["source_population_model_fits"], 423)
        index = {(p["taxon"], tuple(p["three_source_populations"])): p
                 for p in old["all_panels"]}
        self.assertEqual(len(index), 141)
        seen = 0
        for panel in now["original_panel_fit_details"]:
            key = (panel["taxon"], tuple(panel["three_focal_population_ids"]))
            self.assertIn(key, index)
            orig = index[key]
            for fold in panel["heldout_folds"]:
                prior = next(q for q in orig["three_folds_one_source_pop_out"]
                             if q["heldout_source_population"] == fold["heldout"])
                self.assertEqual(fold["n_test_ramets"], prior["n_heldout_ramets"])
                for endpoint in v.ENDPOINTS:
                    self.assertAlmostEqual(
                        fold["loss"][endpoint]["B0"],
                        prior["loss"][endpoint]["B0"], delta=1e-8)
                    seen += 1
        self.assertEqual(seen, 846)


if __name__ == "__main__":
    unittest.main()
