"""OSF original source manifest negative controls; no fictional plant assays."""
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_petunieae_2023_osf_source_v26 import (
    PROJECT, PAPER_DOI, FILES, normalize_api_url,
    inspect_source,
)


def mock_api(duplicate=False):
    def fetch(url):
        if url.endswith("/nodes/zg9cu/"):
            return json.dumps({"data":{"id":"zg9cu","attributes":{
                "title":"Petunieae original analytical source repository","public":True}}}).encode()
        if url==FILES:
            return json.dumps({"data":[
                {"id":"file-111","attributes":{"kind":"file","name":"species_expression.csv",
                 "size":1000,"materialized_path":"/processed/species_expression.csv",
                 "extra":{"hashes":{"md5":"fake_illustrative_not_real_checksum"}}},
                 "links":{"download":"https://files.osf.io/FAKE_EXAMPLE"}},
                {"id":"folder-222","attributes":{"kind":"folder","name":"Processed"},
                 "relationships":{"files":{"links":{"related":{
                    "href":FILES+"path/processed/"}}}}}
            ],"links":{"next":None}}).encode()
        if url==FILES+"path/processed/":
            return json.dumps({"data":[
                {"id":"file-333","attributes":{"kind":"file","name":"pigment.csv",
                 "size":950,"materialized_path":(
                     "/processed/species_expression.csv" if duplicate else "/processed/pigment.csv")},
                 "links":{"download":"https://files.osf.io/FAKE_EXAMPLE"}}
            ],"links":{"next":None}}).encode()
        raise ValueError("unrecognized synthetic source endpoint")
    return fetch


class Petunieae2023SourceGate(unittest.TestCase):
    def test_original_paper_species_panel_not_2019_Iochrominae(self):
        self.assertEqual(PROJECT,"zg9cu")
        self.assertEqual(PAPER_DOI,"10.1098/rspb.2023.0275")

    def test_only_official_public_https_osf_api_traversal(self):
        for unsafe in ("http://api.osf.io/v2/nodes/zg9cu",
                       "https://evil.example/v2/nodes/zg9cu",
                       "https://api.osf.io.evil.example/v2",
                       "https://username:password@api.osf.io/v2",
                       "file:///etc/passwd",FILES.replace("https","ftp")):
            with self.subTest(unsafe=unsafe),self.assertRaises(ValueError):
                normalize_api_url(unsafe)

    def test_complete_nested_manifest_has_no_false_original_reads(self):
        result=inspect_source(mock_api())
        self.assertTrue(result["node_metadata_available"])
        self.assertTrue(result["folder_inventory_complete"])
        self.assertTrue(result["source_manifest_admitted"])
        self.assertEqual(result["status"],"PUBLIC_OSF_MANIFEST_ADMITTED_METADATA_ONLY")
        self.assertEqual(result["source_file_count_metadata_only"],2)
        self.assertEqual(result["small_tabular_candidate_count_metadata_only"],2)
        self.assertFalse(result["original_2023_processed_expression_or_HPLC_bytes_downloaded"])
        self.assertEqual(result["original_species_level_joint_expression_and_pigment_rows_parsed"],0)
        self.assertEqual(result["independent_evolutionary_regain_events_newly_measured"],0)
        self.assertTrue(all(row["file_body_retrieved"] is False
                            for row in result["source_file_records"]))
        self.assertEqual(result["Cirsium_empirical_owner"],"zuizui0223/aza3")
        self.assertTrue(result["source_only_not_cross_clade_predictive_validation"])

    def test_original_file_name_collision_aborts_completeness(self):
        result=inspect_source(mock_api(duplicate=True))
        self.assertFalse(result["folder_inventory_complete"])
        self.assertFalse(result["source_manifest_admitted"])
        self.assertTrue(result["status"].startswith("SOURCE_TREE_"))

    def test_bad_root_identity_fails_and_prevents_assay_claims(self):
        def bad_root(url):
            return json.dumps({"data":{"id":"not_petunieae","attributes":{
                "title":"Unrelated dataset","public":True}}}).encode()
        result=inspect_source(bad_root)
        self.assertFalse(result["node_metadata_available"])
        self.assertFalse(result["source_manifest_admitted"])
        self.assertFalse(result["original_2023_processed_expression_or_HPLC_bytes_downloaded"])

    def test_private_root_does_not_prove_no_molecular_function(self):
        def private_root(url):
            return json.dumps({"data":{"id":"zg9cu","attributes":{
                "title":"Private","public":False}}}).encode()
        result=inspect_source(private_root)
        self.assertFalse(result["source_manifest_admitted"])
        self.assertEqual(result["new_provenance_matched_plant_functional_assays"],0)
        self.assertFalse(result["original_2019_Iochrominae_white_morph_RAR_reanalyzed"])

    def test_technical_access_error_is_not_biological_negative(self):
        def offline(url):
            raise OSError("synthetic testing: source API unavailable")
        data=inspect_source(offline)
        self.assertFalse(data["source_manifest_admitted"])
        self.assertEqual(data["original_species_level_joint_expression_and_pigment_rows_parsed"],0)
        self.assertEqual(data["independent_evolutionary_regain_events_newly_measured"],0)
        self.assertTrue(data["source_only_not_cross_clade_predictive_validation"])


if __name__=="__main__":
    unittest.main()
