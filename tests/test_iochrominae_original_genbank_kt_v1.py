"""Offline proof that GenBank source audit does not manufacture biology."""
import sys
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio import SeqIO
from unittest.mock import patch
from tempfile import TemporaryDirectory
import unittest
import io

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_iochrominae_original_genbank_kt_accessions_v1 as g

def fake(accession,organism="Iochroma candidate",gene="chalcone isomerase"):
    rec=SeqRecord(Seq("ATG"*300),id=accession+".1",name=accession,
                  description=gene)
    rec.annotations["molecule_type"]="DNA"
    rec.annotations["organism"]=organism
    rec.features=[SeqFeature(FeatureLocation(0,900),type="CDS",
                   qualifiers={"gene":[gene],"product":[gene]})]
    return rec

class OriginalAccessionTests(unittest.TestCase):
    def test_exact_original_author_accession_range(self):
        self.assertEqual(g.ACCESSIONS[0],"KT898394")
        self.assertEqual(g.ACCESSIONS[-1],"KT898451")
        self.assertEqual(len(g.ACCESSIONS),58)
        self.assertEqual(len(set(g.ACCESSIONS)),58)

    def test_only_registered_accessions_get_requested(self):
        url=g.accession_query(g.ACCESSIONS[:10])
        self.assertIn("KT898394,KT898395",url)
        self.assertIn("db=nuccore",url)
        with self.assertRaises(ValueError):g.accession_query(["AF123456"])

    def test_annotated_genes_have_distinct_names(self):
        self.assertEqual(g.gene_label(fake("KT898394",gene="chalcone isomerase")),"CHI")
        self.assertEqual(g.gene_label(fake("KT898394",gene="flavanone 3-hydroxylase")),"F3H")
        self.assertEqual(g.gene_label(fake("KT898394",gene="dihydroflavonol 4-reductase")),"DFR")
        self.assertEqual(g.gene_label(fake("KT898394",gene="unknown gene")),"UNRESOLVED_GENE_LABEL")

    def test_source_cds_annotation_does_not_prove_function(self):
        r=g.compact_record(fake("KT898394"))
        self.assertEqual(r["CDS_annotations"],1)
        self.assertFalse(r["full_functional_anthocyanin_pathway_verified"])
        self.assertEqual(r["phenotypic_state_from_this_record"],"UNKNOWN")
        self.assertEqual(r["trait_reversal_event_from_this_record"],"NOT_IDENTIFIABLE")
        self.assertEqual(len(r["source_nucleotide_sequence_sha256"]),64)

    def test_complete_original_source_census_is_a_metadata_census_not_reversals(self):
        allrecords=[fake(ac,organism=f"Iochroma taxon {i%22}",
                         gene=("chalcone isomerase","flavanone 3-hydroxylase",
                               "dihydroflavonol 4-reductase")[i%3])
                    for i,ac in enumerate(g.ACCESSIONS)]
        p=g.audit(allrecords)
        self.assertEqual(p["registered_accession_count"],58)
        self.assertEqual(p["source_records_recovered"],58)
        self.assertEqual(p["unique_organism_names_in_deposited_records"],22)
        self.assertTrue(p["not_independent_evolutionary_trait_reversal_events"])
        inv=p["taxon_gene_accession_inventory"]
        self.assertEqual(len(inv),22)
        self.assertEqual(sum(len(v["CHI"]) for v in inv.values()),20)
        self.assertEqual(sum(sum(len(a) for a in v.values()) for v in inv.values()),58)
        self.assertEqual(sum(len(v["DFR"]) for v in inv.values()),19)

        self.assertFalse(p["phenotypic_color_crosswalk_independently_qualified"])
        self.assertFalse(p["coding_integrity_verified_by_functional_experiment"])
        self.assertFalse(p["regain_or_perpetuation_probability_estimated"])

    def test_missing_accession_refuses_success(self):
        allrecords=[fake(ac) for ac in g.ACCESSIONS[:-1]]
        with self.assertRaisesRegex(ValueError,"not exhaustively"):
            g.audit(allrecords)

    def test_batch_genbank_roundtrip_matches_ids_and_count(self):
        expected=g.ACCESSIONS[:2]
        out=io.StringIO()
        SeqIO.write([fake(x) for x in expected],out,"genbank")
        class Response:
            def __init__(self,text):self.text=text
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,lim):return self.text.encode()
        result=g.fetch_batch(expected,opener=lambda req,timeout:Response(out.getvalue()))
        self.assertEqual([x.id.split(".")[0] for x in result],list(expected))

    def test_refuse_wrong_records(self):
        requested=g.ACCESSIONS[:2]
        out=io.StringIO()
        SeqIO.write([fake(requested[0]),fake("KT898499")],out,"genbank")
        class Response:
            def __enter__(self):return self
            def __exit__(self,*args):return False
            def read(self,lim):return out.getvalue().encode()
        with self.assertRaisesRegex(RuntimeError,"ID_MISMATCH"):
            g.fetch_batch(requested,opener=lambda req,timeout:Response())
if __name__=="__main__":unittest.main()
