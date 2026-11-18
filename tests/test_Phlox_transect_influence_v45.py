"""Adversarial tests for Phlox v45 cross-transect diagnostic (not new CV)."""
import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_Phlox_transect_influence_v45 as v


class V45TransectAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = v.read_source()
        cls.summary = v.report(cls.data)

    def test_exact_parent_git_blob_and_denominators(self):
        self.assertEqual(len(v.SOURCE_GIT_BLOB_SHA1), 40)
        self.assertEqual(
            v.git_blob_sha1(v.SOURCE.read_bytes()), v.SOURCE_GIT_BLOB_SHA1)
        self.assertEqual(self.summary["n_adult_records_example"], 11379)
        self.assertEqual(self.summary["n_offspring_example"], 2496)
        self.assertFalse(self.summary["original_2014_full_paper_cohort_reconciled"])

    def test_receipt_recomputed_from_frozen_parent(self):
        self.assertEqual(self.summary, json.loads(v.FROZEN.read_text()))

    def test_one_corrupted_parent_byte_is_rejected(self):
        raw = v.SOURCE.read_bytes()
        with tempfile.TemporaryDirectory() as d:
            altered = Path(d) / "fake.json"
            altered.write_bytes(raw + b"\n")
            with self.assertRaisesRegex(v.AuditError, "GIT_BLOB_CHANGED"):
                v.read_source(altered)

    def test_mother_and_adult_components_not_mixed(self):
        for rec in self.summary["scenarios"].values():
            self.assertEqual(len(rec["source_per_transect"]), 5)
            for t in rec["source_per_transect"]:
                self.assertAlmostEqual(
                    t["adult"] + t["offspring"], t["total"], places=6)

    def test_aggregate_discordance_not_within_transect_replication(self):
        half = self.summary["scenarios"]["dispersal_half"]
        epi = self.summary["scenarios"]["epistasis_zero"]
        self.assertTrue(half["global_adult_offspring_opposite"])
        self.assertEqual(half["within_transect_opposite_count"], 0)
        self.assertTrue(epi["global_adult_offspring_opposite"])
        self.assertEqual(epi["within_transect_opposite_count"], 2)
        self.assertEqual(
            [t["transect"] for t in epi["source_per_transect"]
             if t["adult_offspring_opposite"]], [3, 5])

    def test_half_dispersal_global_negative_not_stable_to_transect_one(self):
        half = self.summary["scenarios"]["dispersal_half"]
        self.assertEqual(
            half["sign_reversals_under_one_transect_deletion"],
            {"adult": [1], "offspring": [2], "total": [1]})
        cut = half["delete_one_transect_remaining_delta_logL"][0]
        self.assertGreater(cut["adult"], 0)
        self.assertGreater(cut["total"], 0)
        self.assertAlmostEqual(cut["adult"], 27.569182, places=5)

    def test_epistasis_offspring_improvement_not_delete_one_stable(self):
        epi = self.summary["scenarios"]["epistasis_zero"]
        self.assertEqual(
            epi["sign_reversals_under_one_transect_deletion"],
            {"adult": [], "offspring": [2], "total": []})
        self.assertAlmostEqual(
            epi["delete_one_transect_remaining_delta_logL"][1]["offspring"],
            -7.452721, places=5)

    def test_fitness_neutral_and_doubled_dispersal_signs_all_negative(self):
        for name in ("fitness_neutral", "dispersal_double"):
            s = self.summary["scenarios"][name]
            self.assertEqual(s["within_transect_opposite_count"], 0)
            self.assertEqual(s["sign_reversals_under_one_transect_deletion"],
                             {"adult": [], "offspring": [], "total": []})
            self.assertTrue(all(row["total"] < 0
                                for row in s["source_per_transect"]))

    def test_blocked_source_cases_never_become_zero_effects(self):
        self.assertEqual(set(self.summary["source_technical_nonadmissions"]),
                         set(v.BLOCKED))
        self.assertTrue(set(self.summary["scenarios"]).isdisjoint(v.BLOCKED))
        corrupted = copy.deepcopy(self.data)
        corrupted["scenarios"]["assortative_25pct"]["comparison"] = {}
        with self.assertRaisesRegex(v.AuditError, "TECHNICAL_BLOCK"):
            v.validate(corrupted)

    def test_missing_or_duplicate_transect_must_fail(self):
        wrong = copy.deepcopy(self.data)
        a = wrong["scenarios"]["epistasis_zero"]["comparison"]
        a["per_transect_delta_logL"][4]["transect"] = 4
        with self.assertRaisesRegex(v.AuditError, "TRANSECT_GRAIN_CHANGED"):
            v.validate(wrong)

    def test_single_source_cell_mutation_must_fail(self):
        wrong = copy.deepcopy(self.data)
        c = wrong["scenarios"]["dispersal_half"]["comparison"]
        c["per_transect_delta_logL"][0]["Delta_logL_offspring_vs_original_START"] += 20
        with self.assertRaisesRegex(v.AuditError, "TRANSECT_SUM_MISMATCH"):
            v.validate(wrong)

    def test_no_fake_source_independent_replication(self):
        wrong = copy.deepcopy(self.data)
        wrong["n_independent_original_phylogenetic_studies_added"] = 5
        with self.assertRaisesRegex(v.AuditError, "CLAIM_OR_DENOMINATOR_DRIFT"):
            v.validate(wrong)

    def test_no_original_mcmc_or_published_full_cohort(self):
        for key in ("original_2014_75000_step_MCMC_performed",
                    "original_2014_published_full_12024_382_2581_used"):
            wrong = copy.deepcopy(self.data)
            wrong[key] = True
            with self.assertRaisesRegex(v.AuditError, "CLAIM_OR_DENOMINATOR_DRIFT"):
                v.validate(wrong)

    def test_sign_threshold_fails_closed_for_nonfinite(self):
        for x in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaisesRegex(v.AuditError, "NONFINITE_COMPONENT"):
                v.sign(x)


if __name__ == "__main__":
    unittest.main()
