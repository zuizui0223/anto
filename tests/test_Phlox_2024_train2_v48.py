"""Adversarial source-population and fit scope tests for ANTO v48."""
from __future__ import annotations

import copy
import json
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_Phlox_2024_train2_v48 as v


class OriginalPhloxTrainTwoSourcePopulations(unittest.TestCase):
    def test_frozen_source_and_141_423_design(self):
        r = v.fixed_registry()
        self.assertEqual(r["original_ramets"], 963)
        self.assertEqual(r["total_panels"], 141)
        self.assertEqual(r["total_population_holdout_fits"], 423)
        self.assertEqual(r["ridge_lambda"], 10)
        self.assertTrue(r["posthoc_after_v46_and_v47"])

    def test_two_original_pop_training_each_fold_enforced(self):
        train = self.synthetic_rows(["p1", "p2"])
        test = self.synthetic_rows(["p3"])
        self.assertEqual(len(train), 12)
        self.assertEqual(len(test), 6)
        with self.assertRaisesRegex(v.V48Error, "WRONG_TAXON_OR_TWO_TRAIN_POP"):
            v.score_fold(train, test, ["p1", "p1"], "p3", "amoena")
        with self.assertRaisesRegex(v.V48Error, "WRONG_TAXON_OR_TWO_TRAIN_POP"):
            v.score_fold(train, test, ["p1", "p2"], "p2", "amoena")

    @staticmethod
    def synthetic_rows(pop_ids):
        rows = []
        for pi, pop in enumerate(pop_ids):
            for j, garden in enumerate(v.GARDENS):
                for k in range(2):
                    rows.append({
                        "garden": garden, "taxon": "amoena",
                        "population": pop,
                        "genotype": "%s_%s" % (pop, k),
                        "leaf_traits_pc1": -1.0 + 0.9 * k + 0.8 * pi,
                        "leaf_traits_pc2": 0.2 + 0.5 * k + 0.2 * j + 0.17 * pi,
                        "fruits_2019": float(1 + k + pi + j),
                        "flowers_2019": float(2 + 2*k + pi + j),
                    })
        return rows

    def test_test_response_cannot_influence_prediction(self):
        fit = v.build_train(self.synthetic_rows(["p1", "p2"]))
        row = copy.deepcopy(self.synthetic_rows(["p3"])[0])
        for model in v.MODELS:
            before = v.predict(fit, row, model, "fruits_2019")
            row["fruits_2019"] += 100000
            self.assertEqual(before, v.predict(fit, row, model,
                                               "fruits_2019"))
            row["fruits_2019"] -= 100000

    def test_all_three_gardens_have_training_only_intercepts(self):
        fit = v.build_train(self.synthetic_rows(["p1", "p2"]))
        self.assertEqual(set(fit["garden_means"]), set(v.GARDENS))
        self.assertEqual(len(fit["beta"]["G"]["fruits_2019"]), 2)
        self.assertEqual(len(fit["beta"]["E"]["fruits_2019"]), 6)
        self.assertGreater(fit["sd1"], 0)
        self.assertGreater(fit["sd2"], 0)

    def test_genotype_clones_never_cross_original_source_pop(self):
        train = self.synthetic_rows(["p1", "p2"])
        test = self.synthetic_rows(["p3"])
        test[0]["genotype"] = train[0]["genotype"]
        with self.assertRaisesRegex(v.V48Error, "CLONE_LEAKS"):
            v.score_fold(train, test, ["p1", "p2"], "p3", "amoena")

    def test_only_original_fruit_flower_not_signed_biomass(self):
        self.assertEqual(v.ENDPOINTS, ("fruits_2019", "flowers_2019"))
        self.assertEqual(v.MODELS, ("B0", "G", "E"))

    def test_model_predictions_finite_and_loss_positive(self):
        train = self.synthetic_rows(["p1", "p2"])
        test = self.synthetic_rows(["p3"])
        result = v.score_fold(train, test, ["p1", "p2"], "p3", "amoena")
        self.assertEqual(set(result), set(v.ENDPOINTS))
        for by_model in result.values():
            self.assertEqual(set(by_model), set(v.MODELS))
            self.assertTrue(all(x >= 0 for x in by_model.values()))

    def test_test_taxon_mismatch_rejected(self):
        train = self.synthetic_rows(["p1", "p2"])
        test = self.synthetic_rows(["p3"])
        test[0]["taxon"] = "pilosa"
        with self.assertRaisesRegex(v.V48Error, "WRONG_TAXON"):
            v.score_fold(train, test, ["p1", "p2"], "p3", "amoena")

    def test_missing_original_test_garden_rejected(self):
        train = self.synthetic_rows(["p1", "p2"])
        test = [x for x in self.synthetic_rows(["p3"])
                if x["garden"] != v.GARDENS[-1]]
        with self.assertRaisesRegex(v.V48Error, "NOT_ALL_THREE"):
            v.score_fold(train, test, ["p1", "p2"], "p3", "amoena")

    def test_panel_gain_algebra_and_limits(self):
        folds = [
            {"loss": {"fruits_2019": {"B0": 10.0, "G": 8.0}}},
            {"loss": {"fruits_2019": {"B0": 20.0, "G": 10.0}}}
        ]
        self.assertAlmostEqual(v.gain(folds, "fruits_2019", "G"), 0.4)
        with self.assertRaisesRegex(v.V48Error, "INVALID_GAIN_INPUT"):
            v.gain(folds, "fruits_2019", "B0")

    def test_no_valid_p_value_for_overlapping_references(self):
        x = v.summarize([-.10, 0.02, 0.20], -.05)
        self.assertEqual(x["n_at_or_below_deamii"], 1)
        self.assertFalse(x["valid_randomization_p_value"])
        self.assertEqual(x["independent_biological_experiments"], 0)

    def test_actual_original_author_source_lock_if_staged(self):
        staged = os.environ.get("PHLOX_V48_SOURCE")
        if not staged:
            self.skipTest("Original 2024 CSV is staged only in source CI")
        rows, audit = v.load_author_source(staged, v.fixed_registry())
        self.assertEqual(audit["rows"], 963)
        self.assertEqual(audit["genotypes"], 284)
        index, pop = v.source_panel_index(rows)
        self.assertEqual({k: len(q) for k, q in pop.items()},
                         {"amoena": 8, "deamii": 3, "pilosa": 9})
        self.assertEqual(len(index), 20)

    def test_complete_original_source_replay_if_receipt_frozen(self):
        staged = os.environ.get("PHLOX_V48_SOURCE")
        if not staged or not v.RECEIPT.exists():
            self.skipTest("Both pinned original 2024 bytes and frozen receipt required")
        reg = v.fixed_registry()
        rows, audit = v.load_author_source(staged, reg)
        result = v.analyze(rows, audit, reg)
        self.assertEqual(result, json.loads(v.RECEIPT.read_text()))
        self.assertEqual(result["new_source_model_fits"], 423)
        self.assertFalse(result["valid_randomization_p_values"])
        self.assertTrue(result["original_biomass_log1p_technical_nonadmission"])


if __name__ == "__main__":
    unittest.main()
