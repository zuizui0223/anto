"""Source-bounded 2020 Cirsium transcript sequence extraction regression tests."""
import io
import json
from pathlib import Path
import sys
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from extract_cirsium_2020_26_sequences_v6 import (
    load_registry,extract_text_stream,extract_archive
)

class SourceSequenceExtractionTests(unittest.TestCase):
    def test_real_registry_has_26_precisely_registered_identifiers(self):
        d=load_registry()
        self.assertEqual(d["selected_rows"],26)
        self.assertEqual(len(d["records"]),26)
        self.assertEqual(len({r["unigene_id"] for r in d["records"]}),26)
        self.assertEqual(d["source_zip_sha256"],
          "12238a97225d34a8eaf66bf1983030fa7aa0be3c4d7723f88bf44251ae634df8")
        self.assertTrue(all(r["sequence_extracted"] is False for r in d["records"]))

    def test_exact_full_transcript_fasta_join(self):
        fasta=b"""Publisher metadata before FASTA
>TR1|c0_g1_i1 len=9 path=[bla]
ATGCCC
TAA
>TR2|c0_g2_i1 len=12 path=[bla]
ATGTTCGGGTAA
>TR3|c0_g1_i1 len=9 path=[bla]
ATGAAATAG
"""
        found,n=extract_text_stream(io.BytesIO(fasta),{
            "TR1|c0_g1_i1","TR3|c0_g1_i1"}, expected_total_headers=3)
        self.assertEqual(n,3)
        self.assertEqual(found["TR1|c0_g1_i1"]["sequence"],"ATGCCCTAA")
        self.assertEqual(found["TR3|c0_g1_i1"]["source_fasta_length"],9)

    def test_unmatched_target_rejected(self):
        with self.assertRaises(ValueError):
            extract_text_stream(io.BytesIO(
                b">TR1|c0_g1_i1 len=6\nATGTAA\n"),{
                "TR2|c0_g1_i1"},expected_total_headers=1)

    def test_transcript_length_inflation_rejected(self):
        with self.assertRaises(ValueError):
            extract_text_stream(io.BytesIO(
                b">TR1|c0_g1_i1 len=100\nATGTAA\n"),{
                "TR1|c0_g1_i1"},expected_total_headers=1)

    def test_duplicate_transcript_rejected(self):
        with self.assertRaises(ValueError):
            extract_text_stream(io.BytesIO(
                b">TR1|c0_g1_i1 len=6\nATGTAA\n>TR1|c0_g1_i1 len=6\nATGTAA\n"),{
                "TR1|c0_g1_i1"},expected_total_headers=2)

    def test_source_header_count_gate(self):
        with self.assertRaises(ValueError):
            extract_text_stream(io.BytesIO(
                b">TR1|c0_g1_i1 len=6\nATGTAA\n"),{
                "TR1|c0_g1_i1"},expected_total_headers=104890)

    def test_bogus_archive_fails_uncompressed_source_gate(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w",compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr("ao0c04001_si_002.txt",">TR1|c0_g1_i1 len=6\nATGTAA\n")
        with self.assertRaisesRegex(ValueError,"missing"):
            extract_archive(b.getvalue(),load_registry())


if __name__=="__main__":
    unittest.main()
