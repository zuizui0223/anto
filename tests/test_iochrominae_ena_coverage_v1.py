"""Offline source-coverage admission; original actual file tested in CI after GET."""
import json
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import check_iochrominae_ena_coverage_v1 as c


class SourceCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected=json.loads((ROOT/"data/external/iochrominae_2016_genbank_58_source_coverage_v1.json").read_text(encoding="utf-8"))

    def test_gene_blocks_are_registered_but_not_functionality(self):
        self.assertEqual(c.classify_number("KT898394"),"F3h")
        self.assertEqual(c.classify_number("KT898415"),"F3h")
        self.assertEqual(c.classify_number("KT898416"),"Chi")
        self.assertEqual(c.classify_number("KT898437"),"Chi")
        self.assertEqual(c.classify_number("KT898438"),"Dfr accession block")
        self.assertEqual(c.classify_number("KT898451"),"Dfr accession block")

    def test_outside_paper_accessions_fail(self):
        for acc in ("KT898393","KT898452","XX898400"):
            with self.subTest(acc=acc),self.assertRaises(ValueError):
                c.classify_number(acc)

    def test_missing_new_Dfr_accession_does_not_imply_gene_absence(self):
        a=self.expected["per_taxon"]
        self.assertEqual(len(a),22)
        self.assertEqual(sum(x["has_new_Dfr_accession_in_this_span"] for x in a),14)
        self.assertEqual(sum(not x["has_new_Dfr_accession_in_this_span"] for x in a),8)
        self.assertTrue(all(x["no_Dfr_in_this_new_accession_span_is_gene_loss"] is False for x in a))

    def test_two_unannotated_records_are_not_record_deletion(self):
        a=self.expected["observed"]
        self.assertEqual(a["successful_accessions"],58)
        self.assertEqual(a["annotated_CDS_features"],56)
        self.assertEqual(a["records_without_CDS_feature"],2)
        self.assertEqual({r["accession"] for r in a["source_no_CDS_feature_accessions"]},
                         {"KT898443","KT898446"})
        self.assertFalse(self.expected["independent_real_reexpression_rate_estimate"])

    def test_source_artifact_reference_is_exact_not_synthetic(self):
        self.assertEqual(len(c.SOURCE_RECEIPT_SHA256),64)
        self.assertEqual(self.expected["provenance"]["github_workflow_run_id"],39229422231)

if __name__=="__main__":
    unittest.main()
