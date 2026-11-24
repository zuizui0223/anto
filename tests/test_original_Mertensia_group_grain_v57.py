"""Mertensia source-only Excel identifiers, no original seed/visit values."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_original_Mertensia_exp3_sameplant_grain_v57 as p

class Original2016MertensiaXLSX(unittest.TestCase):
    def test_dual_source_digest_frozen_before_composite_IDs(self):
        r=p.check_receipt()
        self.assertEqual(p.author.gitblob(p.SOURCE.read_bytes()),
                         "1071510aa742cb3f17b45c1fd6d8e0cc5412a746")
        self.assertEqual(r["original_dryad_doi"],"10.7280/D16D7Z")
        self.assertTrue(all(x["byte_identity_validated_by_both_archives"]
                            for x in r["dual_digest_identical_original_files"]))

    def test_composite_source_identifier_distinct_from_individual(self):
        self.assertEqual(p.CATEGORICAL,("week","array","treat","plantID"))
        self.assertIn("mVisRate",p.MEASURED)
        self.assertIn("mseed",p.MEASURED)
        self.assertNotIn("plantID",p.MEASURED)

    def test_original_excel_column_number_is_deterministic(self):
        self.assertEqual(p.colnumber("A"),0)
        self.assertEqual(p.colnumber("Z"),25)
        self.assertEqual(p.colnumber("AA"),26)
        with self.assertRaisesRegex(p.OriginalGrainError,"INVALID_XLSX_COLUMN_LABEL"):
            p.colnumber("A1")

    def test_reproductive_magnitudes_not_copied_into_source_receipt(self):
        content=open(p.__file__,encoding="utf8").read()
        self.assertIn('"source_original_measurement_magnitudes_examined":False',content)
        self.assertIn('"source_original_visitation_seed_association_fitted":False',content)
        self.assertIn('"original_plant_ID_tokens_or_outcomes_exported":False',content)

if __name__=="__main__":unittest.main()
