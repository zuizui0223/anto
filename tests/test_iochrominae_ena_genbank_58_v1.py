"""Strict, synthetic unit checks of original-genbank record ingestion."""
import io
from pathlib import Path
import sys
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation,SeqFeature
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_iochrominae_ena_genbank_58_v1 as a


def example(acc="KT898394",gene="DFR"):
    r=SeqRecord(Seq("ATGGCTTAA"),id=acc+".1",name=acc,description="Iochrominae example coding",
        annotations={"molecule_type":"DNA","organism":"Iochroma edule"})
    r.features=[SeqFeature(FeatureLocation(0,9),type="CDS",qualifiers={
        "gene":[gene],"product":["annotated protein"],"translation":["MA"]})]
    b=io.StringIO()
    SeqIO.write([r],b,"embl")
    return b.getvalue().encode()


class ENARecordTests(unittest.TestCase):
    def test_exhaustive_source_span(self):
        self.assertEqual(len(a.ACCESSIONS),58)
        self.assertEqual(a.ACCESSIONS[0],"KT898394")
        self.assertEqual(a.ACCESSIONS[-1],"KT898451")
        self.assertEqual(len(set(a.ACCESSIONS)),58)

    def test_parse_realistic_embl_record_without_phenotype_filling(self):
        r=a.parse_record("KT898394",example())
        self.assertEqual(r["organism_as_submitted"],"Iochroma edule")
        self.assertEqual(r["sequence_length_nt"],9)
        self.assertEqual(r["cds_annotations"][0]["gene_annotations"],["DFR"])
        self.assertTrue(r["cds_annotations"][0]["translation_annotated"])
        self.assertFalse(r["functional_enzyme_activity_measured"])
        self.assertFalse(r["real_regained_flower_colour_measured"])

    def test_wrong_accession_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"mismatches"):
            a.parse_record("KT898395",example())

    def test_empty_sequence_fail_closed(self):
        with self.assertRaises(ValueError):
            a.parse_record("KT898394",b"")

    def test_unregistered_accession_never_fetched(self):
        with self.assertRaisesRegex(ValueError,"Unregistered"):
            a.get_record("KT898452")

    def test_summary_does_not_promote_a_single_gene_to_phenotype_repair(self):
        p=a.parse_record("KT898394",example())
        errs=[{"accession":acc,"source_failure_type":"TEST_NOT_FETCHED",
               "source_failure_message":"fixture"} for acc in a.ACCESSIONS[1:]]
        s=a.summarize([p],errs)
        self.assertEqual((s["expected_original_accession_count"],s["retrieved_count"],s["not_retrieved_count"]),(58,1,57))
        self.assertFalse(s["all_records_retrieved"])
        self.assertEqual(s["reported_CDS_gene_annotations"],{"DFR":1})
        self.assertFalse(s["flower_trait_regain_and_persistence_data_collected"])
        self.assertTrue(s["original_author_coding_integrity_result_not_independently_retested"])
        self.assertTrue(s["no_false_independent_clade_or_transition_replication"])

    def test_fake_duplicate_accession_in_summary_is_rejected(self):
        p=a.parse_record("KT898394",example())
        errs=[{"accession":acc,"source_failure_type":"not-retrieved",
               "source_failure_message":"fixture"} for acc in a.ACCESSIONS]
        with self.assertRaisesRegex(ValueError,"overlap"):
            a.summarize([p],errs)

if __name__=="__main__":
    unittest.main()
