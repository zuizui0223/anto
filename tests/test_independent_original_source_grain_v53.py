"""Strict source-first original-study eligibility and no-outcome-reading tests."""
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_independent_pollinator_reproduction_sources_v53 as v


class OriginalSourceIdentityGate(unittest.TestCase):
    def test_two_independent_original_candidate_archives(self):
        x=v.registry()
        self.assertEqual([z["zenodo_record_id"] for z in x["candidates"]],
                         [4970401,6671199])
        self.assertEqual(sum(len(z["files_md5"]) for z in x["candidates"]),10)
        self.assertFalse(x["fit_any_model"])
        self.assertFalse(x["may_read_reproductive_outcome_values"])

    def test_prior_traine_dryad_403_stop_not_retried(self):
        x=v.registry()
        z=x["prior_brassica_traine_egwe_access_stop"]
        self.assertEqual(z["original_dryad_file_id"],4644044)
        self.assertEqual(z["exact_file_stream_status"],403)
        self.assertFalse(z["original_csv_accessed"])
        self.assertFalse(z["should_retry_unauthenticated"])

    def test_source_names_registered_without_guessing(self):
        x=v.registry()
        snow=x["candidates"][0]["files_md5"]
        self.assertEqual(snow["AllSpeciesVisitation.csv"],
                         "5b7ec80a6183546b203a4830311157df")
        sun=x["candidates"][1]["files_md5"]
        self.assertEqual(sun["Database_Productivity.csv"],
                         "aeea30a21e495498295241bcfd8c693d")

    def test_only_zenodo_https_permitted(self):
        self.assertTrue(v.zenodo_url("https://zenodo.org/api/records/4970401"))
        self.assertFalse(v.zenodo_url("https://evil.example/"))
        self.assertFalse(v.zenodo_url("http://zenodo.org/api"))
        self.assertFalse(v.zenodo_url("https://bad@zenodo.org/abc"))

    def test_column_audit_uses_no_seed_or_visit_values(self):
        sample=(b"Plant_ID,Plot,Visits,Seed_Set\n"
                b"A1,North,9999,100000000\n"
                b"A2,South,0,NaN\n")
        got=v.table_schema("test.csv",sample)
        self.assertEqual(got["n_original_source_rows"],2)
        self.assertIn("Plant_ID",got["candidate_plant_fields"])
        self.assertIn("Plot",got["potential_identity_or_block_field_cardinality"])
        self.assertNotIn("Seed_Set",got["potential_identity_or_block_field_cardinality"])
        self.assertNotIn("Visits",got["potential_identity_or_block_field_cardinality"])
        self.assertTrue(got["no_numeric_flower_visit_or_fruit_seed_outcomes_examined"])

    def test_do_not_auto_admit_matching_number_tokens(self):
        c={"code":"X","original_data_doi":"test","zenodo_record_id":1}
        payload={"visits.csv":b"Plant_ID,Visitation\n1,4\n2,3\n",
                 "seed.csv":b"Plant_ID,Seed_Count\n1,999\n2,0\n"}
        got=v.source_audit(c,payload)
        self.assertEqual(got["n_source_files_exact_MD5_verified"],2)
        self.assertEqual(got["same_token_plant_key_candidate_comparisons"][0]
                         ["n_exact_shared_original_plant_id_tokens"],2)
        self.assertFalse(got["direct_individual_plant_visitation_to_seed_set_join_ADMITTED"])
        self.assertFalse(got["heldout_block_and_temporal_semantics_verified"])

    def test_bad_csv_row_shape_fails_closed(self):
        with self.assertRaisesRegex(v.SourceGateError,"CSV_INVALID_ROW_SHAPE"):
            v.table_schema("bad.csv",b"plant_id,plot\nx,y,z\n")

    def test_binary_metadata_xlsx_only_authenticated_not_parsed(self):
        x=v.table_schema("MetaData.xlsx",b"not a CSV document")
        self.assertFalse(x["parsed_as_csv"])
        self.assertFalse(x["original_binary_values_parsed"])

    def test_result_receipt_after_source_freeze(self):
        if not v.RECEIPT.exists():
            self.skipTest("Original source-gate receipt not committed")
        out=json.loads(v.RECEIPT.read_text())
        self.assertFalse(out["no_original_reproductive_or_visitor_values_parsed_for_effect"] is False)
        self.assertTrue(out["no_independent_original_fitted_predictive_studies_added"])
        self.assertEqual(len(out["records"]),2)
        for item in out["records"]:
            self.assertFalse(item["direct_individual_plant_visitation_to_seed_set_join_ADMITTED"])
            self.assertFalse(item["source_visit_or_reproductive_values_modeled"])


if __name__=="__main__":
    unittest.main()
