"""Do not confuse published source metadata with downloaded data or assays."""
import io
import subprocess
import tempfile
import json
import sys
import unittest
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_larter_2019_original_dryad_v23 import (
    DATA_DOI,REGISTERED_BYTES,REGISTERED_MD5,REGISTERED_FILE_ID,
    original_bytes_valid, registered_api_details, attempt_original_source,
    _find_exact_rar_in_bundle,receipt,
)


class OriginalLarter2019DryadGateV23(unittest.TestCase):
    def setUp(self):
        self.metadata={"dataset_metadata_retrieved":True,
                       "version_manifest_retrieved":True,
                       "original_rar_metadata_match":True,
                       "file_metadata_status":"VERIFIED_OFFICIAL_MANIFEST"}

    def test_source_exact_identity_never_pretend_to_resequence(self):
        self.assertEqual(DATA_DOI,"10.5061/dryad.p5dq84v")
        self.assertEqual(REGISTERED_FILE_ID,108456)
        self.assertEqual(REGISTERED_BYTES,3054999)
        self.assertEqual(REGISTERED_MD5,"76b46e384fb7c9bf1ef3fbd7d1e5d2f0")
        for junk in (b"",b"Rar!\x1a\x07\x00",b"not the original",b"<html>Validating...</html>",b"PK\x03\x04"):
            self.assertFalse(original_bytes_valid(junk))

    def test_harmless_nonadmitted_html_or_wrong_payloads_are_not_source(self):
        def fake(url):
            if "download" in url: return b"<html>Validating</html>"
            return b"invalid/not original"
        raw,mode,attempts=attempt_original_source(fake)
        self.assertIsNone(raw)
        self.assertIsNone(mode)
        self.assertEqual(len(attempts),3)
        self.assertFalse(any(x["original_md5_matched"] for x in attempts))

    def test_blocked_requests_remain_technical_not_negative_biology(self):
        def offline(url):
            raise OSError("blocked source in synthetic negative control")
        raw,mode,attempts=attempt_original_source(offline)
        d=receipt(self.metadata,raw,mode,attempts)
        self.assertFalse(d["original_RAR_bytes_verified"])
        self.assertTrue(d["technical_blocker_not_biological_negative"])
        self.assertFalse(d["white_morph_versus_fixed_white_new_statistical_result"])
        self.assertFalse(d["trait_regain_or_multigeneration_retention_estimated"])
        self.assertEqual(d["raw_gene_expression_sample_rows_parsed"],0)
        self.assertTrue(d["not_Wheeler2025_red_blue_backcross"])

    def test_fake_bundle_with_wrong_archive_bytes_never_passes(self):
        mem=io.BytesIO()
        with zipfile.ZipFile(mem,"w") as z:
            z.writestr("original-name.rar",b"Rar!\x1a\x07\x00"+"fake".encode())
        self.assertIsNone(_find_exact_rar_in_bundle(mem.getvalue()))

    def test_true_unknown_original_cannot_be_fabricated_using_metadata(self):
        d=receipt(self.metadata,None,None,[])
        self.assertFalse(d["biological_source_admitted_for_future_raw_reanalysis"])
        self.assertFalse(d["original_RAR_bytes_verified"])
        self.assertFalse(d["biological_2019_expression_reanalysed_in_this_gate"])
        self.assertEqual(d["Cirsium_empirical_owner"],"zuizui0223/aza3")

    def test_official_metadata_requires_exact_registered_doi_and_file(self):
        def wrong(url):
            if url.endswith("p5dq84v"):
                return json.dumps({"identifier":"doi:wrong","_links":{
                    "stash:version":{"href":"/api/v2/versions/28856"}}}).encode()
            return b"{}"
        x=registered_api_details(wrong)
        self.assertFalse(x["dataset_metadata_retrieved"])
        self.assertFalse(x["original_rar_metadata_match"])

    def test_invalid_local_archive_is_rejected_without_network(self):
        """A fake file never qualifies even through the supplied-original route."""
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/"source.rar"
            file.write_bytes(b"Rar!\x1a\x07\x00"+b"fake")
            receipt_path=Path(d)/"result.json"
            command=[sys.executable,str(ROOT/"scripts/audit_larter_2019_original_dryad_v23.py"),
                     "--from-local-rar",str(file),"--output",str(receipt_path)]
            result=subprocess.run(command,capture_output=True,text=True,timeout=10)
            self.assertNotEqual(result.returncode,0)
            self.assertIn("original registered size",result.stderr)
            self.assertFalse(receipt_path.exists())

    def test_frozen_access_denial_is_not_an_empirical_negative_result(self):
        path=ROOT/"results/iochroma_2019_source_access_v23_20261028.json"
        x=json.loads(path.read_text(encoding="utf-8"))
        self.assertTrue(x["source_metadata"]["original_rar_metadata_match"])
        self.assertFalse(x["original_RAR_bytes_verified"])
        self.assertFalse(x["biological_source_admitted_for_future_raw_reanalysis"])
        self.assertFalse(x["biological_2019_expression_reanalysed_in_this_gate"])
        self.assertFalse(x["white_morph_versus_fixed_white_new_statistical_result"])
        self.assertEqual([z["retrieval_status"] for z in x["source_acquisition_attempts"]],
                         ["HTTP_401","HTTP_403","HTTP_401"])
        self.assertTrue(x["technical_blocker_not_biological_negative"])
        self.assertFalse(x["source_file_raw_bytes_rehosted_in_git"])

    def test_study_identity_and_assay_gate_are_separate(self):
        x=receipt({},None,None,[{"route":"file_API",
                                 "retrieval_status":"HTTP_401",
                                 "original_md5_matched":False}])
        self.assertEqual(x["paper_observational_design"],
                         "Larter2019 white morphs and fixed white taxa")
        self.assertFalse(x["source_assay_functionally_validated"])
        self.assertFalse(x["source_file_raw_bytes_rehosted_in_git"])
        self.assertTrue(x["technical_blocker_not_biological_negative"])


if __name__=="__main__":
    unittest.main()
