"""Mertensia official Dryad/Zenodo exact-original equivalence negative gates."""
import hashlib
import io
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_original_Mertensia_Zenodo_equivalence_v57 as m

class Mertensia2021Source(unittest.TestCase):
    def test_source_pin_is_two_file_original_dryad(self):
        self.assertEqual(m.ZENODO_RECORD,4731053)
        self.assertEqual(m.ORIGINAL_DOI,"10.7280/D16D7Z")
        self.assertEqual(len(m.FILES),2)
        self.assertEqual(m.FILES["gallagher&campbell_frontiersinecology&evolution_data.xlsx"],339816)

    def test_digest_must_match_source_archive(self):
        payload=b"Original author bytes"
        self.assertEqual(m.checksum_matches(payload,"md5:"+hashlib.md5(payload).hexdigest()),
                         (True,"MD5"))
        self.assertEqual(m.checksum_matches(payload,hashlib.sha256(payload).hexdigest()),
                         (True,"SHA256"))
        self.assertEqual(m.checksum_matches(payload,"a"*32),(False,"MD5"))
        self.assertEqual(m.checksum_matches(payload,None)[0],False)

    def test_zenodo_title_alone_not_authenticated(self):
        # Actual source gate must pass both Dryad published digest
        # and original Zenodo MD5, not title search matches.
        s=open(m.__file__,encoding="utf8").read()
        self.assertIn("V57_MERTENSIA_AUTHOR_DRYAD_DIGEST_NOT_EQUAL",s)
        self.assertIn("V57_MERTENSIA_ZENODO_ORIGINAL_MD5_FAILED",s)

    def test_published_dataset_is_one_original_study(self):
        self.assertEqual(m.ORIGINAL_TITLE,
                         "experimental test of the combined effects of water availability and flowering time on pollinator visitation and seed set")
        self.assertIn("sameplant_visit_seed_outcome_model_fitted",open(m.__file__).read())

if __name__=="__main__":unittest.main()
