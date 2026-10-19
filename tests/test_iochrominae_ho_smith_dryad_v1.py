"""Synthetic negative controls for an official Dryad source-only intake."""
import io
import json
from pathlib import Path
import sys
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_iochrominae_ho_smith_original_dryad_v1 as d

def fixture(dfr_header="28 1179",missing=None):
    f={"EVOL_CHI_SEQFILE.txt":"28 717\n"+"ATGC"*350,
       "EVOL_DFR_SEQFILE.txt":dfr_header+"\n"+"ATGC"*350,
       "F3H_SEQFILE.txt":"28 1113\n"+"ATGC"*350,
       "README_for_EVOL_CHI_SEQFILE.txt":"original source legend"}
    if missing:f.pop(missing)
    bio=io.BytesIO()
    with zipfile.ZipFile(bio,"w") as z:
        for p,s in f.items():z.writestr(p,s)
    return bio.getvalue()


class SourceIntakeTests(unittest.TestCase):
    def test_exact_three_original_gene_file_names_and_headers(self):
        report=d.inspect_zip(fixture())
        self.assertEqual(set(report["target_gene_file_metadata"]),set(d.EXPECTED_FILES))
        self.assertEqual(report["target_gene_file_metadata"]["EVOL_CHI_SEQFILE.txt"]["declared_alignment_sites"],717)
        self.assertEqual(report["target_gene_file_metadata"]["EVOL_DFR_SEQFILE.txt"]["declared_sequences"],28)
        self.assertEqual(report["target_gene_file_metadata"]["F3H_SEQFILE.txt"]["declared_alignment_sites"],1113)

    def test_missing_one_gene_alignment_is_not_success(self):
        with self.assertRaisesRegex(ValueError,"missing"):d.inspect_zip(fixture(missing="F3H_SEQFILE.txt"))

    def test_source_identity_and_limits(self):
        m={"identifier":"doi:10.5061/dryad.49rs8","versionStatus":"published"}
        r=d.evaluate(m,fixture())
        self.assertEqual(r["biological_inference_level"],"SOURCE_FILE_RETRIEVAL_AND_HEADER_ONLY")
        self.assertEqual(r["outcomes_causally_measured_here"],[])
        self.assertFalse(r["original_sequences_rehosted"])
        self.assertTrue(r["no_ancestral_regain_or_phenotypic_persistence_inferred"])

    def test_wrong_doi_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"dataset"):d.evaluate({"identifier":"doi:10.5061/dryad.WRONG"},fixture())

    def test_url_restriction(self):
        with self.assertRaisesRegex(ValueError,"official"):d.get_bytes("https://example.com/other.txt")

    def test_unknown_format_must_remain_unknown(self):
        x=d.read_alignment_header(b"not-an-alignment\nACGT\n")
        self.assertEqual(x["format"],"UNCLASSIFIED")

    def test_fasta_is_not_silently_interpreted_as_coding(self):
        x=d.read_alignment_header(b">species_A\nATGC\n>species_B\nATGC\n")
        self.assertEqual(x["format"],"FASTA")
        self.assertIsNone(x["declared_alignment_sites"])

if __name__=="__main__":unittest.main()
