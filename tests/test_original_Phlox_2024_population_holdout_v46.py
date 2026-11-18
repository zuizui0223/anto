"""Adversarial source, inference-grain and leakage tests for ANTO v46."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_original_Phlox_2024_population_holdout_v46 as v


class Original2024PhloxPopulationTransfer(unittest.TestCase):
    def test_frozen_original_repo_and_source_registry(self):
        r = v.source_registry()
        self.assertEqual(r["source_git_blob_sha1"],
                         "4c780b0abe644ac69dae86a0450eb598c160626b")
        self.assertEqual(r["source_size"], 175377)
        self.assertEqual(r["ridge_alpha"], 10)
        self.assertFalse(r["flower_colour_assay"])

    def test_blob_sha_bytes_sensitive(self):
        self.assertNotEqual(v.git_blob_sha1(b"original"),
                            v.git_blob_sha1(b"original\n"))
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "fake.csv"
            p.write_bytes(b"not original")
            with self.assertRaisesRegex(v.V46Error, "BYTE_COUNT_DRIFT"):
                v.acquire_raw(v.source_registry(), p)

    def test_normal_equation_positive_ridge_exact(self):
        beta = v.solve([[12.0, 0.0], [0.0, 10.0]], [6.0, 0.0])
        self.assertAlmostEqual(beta[0], 0.5)
        self.assertAlmostEqual(beta[1], 0.0)
        with self.assertRaisesRegex(v.V46Error, "SINGULAR_RIDGE"):
            v.solve([[0.0]], [1.0])

    def synthetic_train(self):
        rows = []
        for gi, garden in enumerate(v.GARDENS):
            for ti, taxon in enumerate(v.TAXA):
                for k in range(3):
                    rows.append({"garden": garden, "taxon": taxon,
                        "population": str(gi * 10 + ti),
                        "genotype": str(gi * 10 + ti) + "_" + str(k),
                        "leaf_traits_pc1": -1.1 + 1.17 * k + 0.3 * ti,
                        "leaf_traits_pc2": -0.8 + 0.71 * k + 0.2 * gi,
                        "fruits_2019": float(3 + k + gi + ti),
                        "flowers_2019": float(7 + 2 * k + gi + ti)})
        return rows

    def test_required_nine_source_cells_and_group_means(self):
        fitted = v.fit_train(self.synthetic_train())
        self.assertEqual(len(fitted["group_mean"]), 9)
        self.assertEqual(set(fitted["coefficients"]), set(v.MODELS[1:]))
        for m, co in fitted["coefficients"].items():
            self.assertEqual(len(co["fruits_2019"]),
                {"G": 2, "T": 6, "E": 6, "TE": 18}[m])

    def test_no_heldout_outcome_used_by_predictor(self):
        fitted = v.fit_train(self.synthetic_train())
        heldout = dict(self.synthetic_train()[0])
        for m in v.MODELS:
            p1 = v.prediction(fitted, heldout, m, "fruits_2019")
            heldout["fruits_2019"] += 100000
            p2 = v.prediction(fitted, heldout, m, "fruits_2019")
            self.assertEqual(p1, p2)
            heldout["fruits_2019"] -= 100000

    def test_split_requires_twenty_distinct_source_pops(self):
        with self.assertRaisesRegex(v.V46Error, "TWENTY_SOURCE_POPULATIONS"):
            v.cohort_folds(self.synthetic_train())

    def test_original_column_order_is_not_fuzzy(self):
        self.assertEqual(len(v.KEYS), 26)
        self.assertEqual(v.KEYS[20:23],
                         ("leaf_traits_pc1", "leaf_traits_pc2",
                          "leaf_traits_pc3"))

    def test_no_biomass_proxy_imputation(self):
        self.assertNotIn(v.BLOCKED_ENDPOINT, v.ACCEPTED_ENDPOINTS)
        self.assertEqual(v.BLOCKED_ENDPOINT, "biomass_2019")

    def test_original_clones_not_independent_studies(self):
        self.assertEqual(v.TAXA, ("amoena", "deamii", "pilosa"))
        self.assertEqual(len(v.GARDENS), 3)
        self.assertEqual(v.MODELS, ("B0", "G", "T", "E", "TE"))

    def test_slope_index_order_and_nonoverlap(self):
        row = {"garden": v.GARDENS[2], "taxon": v.TAXA[1]}
        self.assertEqual([v.slope_offset(m, row)
                          for m in v.MODELS[1:]], [0, 2, 4, 14])
        with self.assertRaisesRegex(v.V46Error, "UNKNOWN_MODEL"):
            v.slope_offset("something", row)

    def test_original_source_replay_when_author_bytes_are_available(self):
        p = os.environ.get("PHLOX_V46_SOURCE")
        if not p:
            self.skipTest("Original author CSV not staged for offline run")
        reg = v.source_registry()
        rows, audit = v.parse_source(v.acquire_raw(reg, p), reg)
        self.assertEqual(audit["rows"], 963)
        self.assertEqual(audit["genotypes"], 284)
        self.assertEqual(audit["populations"], 20)
        self.assertEqual(audit["biomass_original_negative_n"], 13)
        self.assertEqual(audit["biomass_original_minimum"], -233.0)
        folds = v.cohort_folds(rows)
        self.assertEqual(len(folds), 20)
        self.assertEqual(sum(len(te) for _, tr, te in folds), 963)
        for _, tr, te in folds:
            self.assertTrue(set(r["genotype"] for r in tr).isdisjoint(
                set(r["genotype"] for r in te)))

    def test_frozen_receipt_replay_when_author_bytes_available(self):
        p = os.environ.get("PHLOX_V46_SOURCE")
        if not p or not v.RECEIPT.exists():
            self.skipTest("Source and derived receipt both required")
        reg = v.source_registry()
        rows, audit = v.parse_source(v.acquire_raw(reg, p), reg)
        output = v.report(rows, audit, reg)
        self.assertEqual(output, json.loads(v.RECEIPT.read_text()))
        self.assertFalse(output["primary_fruit_transfer_pass_models"])
        self.assertEqual(output["nonadmitted_secondary_biomass"]["status"],
                         "BLOCKED_PREDECLARED_LOG1P_DOMAIN")
        self.assertEqual(set(output["endpoints"]), set(v.ACCEPTED_ENDPOINTS))


if __name__ == "__main__":
    unittest.main()
