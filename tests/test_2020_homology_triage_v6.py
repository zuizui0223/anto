"""Strict tests for protein-similarity triage using real published source sequences."""
from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_2020_homology_triage_v6 import (
    ref_proteins,six_frame_segments,alignment_summary,triage_sequence,
    run_study,read_fasta,
)
from Bio.Seq import Seq

REF=ROOT/"data/molecular/asteraceae_reference_cds_v5.fasta"
SOURCE=ROOT/"data/molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta"
META=ROOT/"data/molecular/cirsium_japonicum_2020_26_sequence_provenance_v6.json"


class ProteinSimilarityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prot=ref_proteins(REF)
        cls.refnt=read_fasta(REF)

    def test_validated_source_protein_sizes(self):
        self.assertEqual(len(self.prot["DFR"]),342)
        self.assertEqual(len(self.prot["ANS"]),356)

    def test_dfr_reference_scores_as_dfr_not_ans(self):
        nt=next(x for k,x in self.refnt.items() if k.startswith("EF600682.1|"))
        score=triage_sequence(nt,self.prot)
        self.assertEqual(score["best_reference_family"],"DFR")
        self.assertEqual(score["family_screen"],"STRONG_REFERENCE_SIMILARITY_CANDIDATE")
        self.assertAlmostEqual(score["highest_score_features"]["aligned_identity"],1)
        self.assertAlmostEqual(score["highest_score_features"]["reference_coverage"],1)

    def test_ans_reference_scores_as_ans_not_dfr(self):
        nt=next(x for k,x in self.refnt.items() if k.startswith("AY547342.1|"))
        score=triage_sequence(nt,self.prot)
        self.assertEqual(score["best_reference_family"],"ANS")
        self.assertEqual(score["family_screen"],"STRONG_REFERENCE_SIMILARITY_CANDIDATE")
        self.assertAlmostEqual(score["highest_score_features"]["aligned_identity"],1)

    def test_reverse_complement_detected(self):
        nt=next(x for k,x in self.refnt.items() if k.startswith("EF600682.1|"))
        rev=str(Seq(nt).reverse_complement())
        x=triage_sequence(rev,self.prot)
        self.assertEqual(x["best_reference_family"],"DFR")
        self.assertEqual(x["highest_score_features"]["orientation"],"-")

    def test_no_long_orf_does_not_report_gene_absent(self):
        x=triage_sequence("ATGTAA",self.prot)
        self.assertEqual(x["family_screen"],"NO_ORF_LONG_ENOUGH")
        self.assertFalse(x["orthology_confirmed"])

    def test_random_alanine_sequence_not_certificate(self):
        x=triage_sequence("GCT"*450+"TAA",self.prot)
        self.assertNotEqual(x["family_screen"],"STRONG_REFERENCE_SIMILARITY_CANDIDATE")

    def test_local_alignment_matches_known_protein(self):
        x=alignment_summary(self.prot["DFR"],self.prot["DFR"])
        self.assertEqual(x["exact_matches"],342)
        self.assertEqual(x["aligned_aa"],342)

    def test_real_26_transcripts_stay_non_focal(self):
        x=run_study(SOURCE,META,REF)
        self.assertEqual(x["n_transcripts"],26)
        self.assertFalse(x["Ryukyu_2026_sequences_included"])
        self.assertFalse(x["orthology_resolved"])
        self.assertFalse(x["floral_reexpression_proven"])
        self.assertTrue(all(not r["orthology_confirmed"] for r in x["records"]))

    def test_swapped_real_source_sequence_rejected(self):
        orig=SOURCE.read_text(encoding="utf-8")
        i=orig.find("\n",orig.find(">"))+1
        replacement = ("C" if orig[i]!="C" else "A")
        with TemporaryDirectory() as d:
            p=Path(d)/"altered.fasta"
            p.write_text(orig[:i]+replacement+orig[i+1:])
            with self.assertRaisesRegex(ValueError,"corrupted"):
                run_study(p,META,REF)


if __name__=="__main__":
    unittest.main()
