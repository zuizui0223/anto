"""Molecular source acceptance and negative controls; synthetic candidates are NOT observations."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_molecular_screen_v5 import (
    read_fasta, cds_quality, validate_references, load_register,
    inspect_optional_candidates, screen,
)

FASTA=ROOT/"data/molecular/asteraceae_reference_cds_v5.fasta"
PROV=ROOT/"data/molecular/asteraceae_reference_provenance_v5.json"
REGISTER=ROOT/"data/molecular/ryukyu_gene_screen_register_v5.csv"
RUNS=ROOT/"data/public_sra/prjna1311153_run_manifest_20261017.csv"


class MolecularReferenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs=read_fasta(FASTA)
        cls.r=load_register(REGISTER,RUNS)

    def test_three_real_published_complete_cds_sequences(self):
        result=validate_references(FASTA,PROV)
        self.assertEqual(set(result),{"DFR","ANS","GMYC1"})
        self.assertEqual(result["DFR"]["nucleotide_length"],1029)
        self.assertEqual(result["ANS"]["nucleotide_length"],1071)
        self.assertEqual(result["GMYC1"]["nucleotide_length"],1602)
        self.assertTrue(all(v["reference_only_no_focal_sequence"] for v in result.values()))
        self.assertTrue(all(v["sequence_qc_status"]=="OPEN_READING_FRAME_CANDIDATE_ONLY"
                            for v in result.values()))

    def test_six_ryukyu_runs_predeclare_eighteen_unmeasured_loci(self):
        self.assertEqual(len(self.r),18)
        runs={r["run_accession"] for r in self.r}
        self.assertEqual(len(runs),6)
        self.assertEqual({r["status"] for r in self.r},{"NOT_TESTED"})
        self.assertTrue(all(r["floral_expression_measured"]=="false" for r in self.r))

    def test_no_focal_coding_sequences_yet(self):
        record=screen(FASTA,PROV,REGISTER,RUNS)
        self.assertEqual(record["n_candidate_sequences_supplied"],0)
        self.assertEqual(record["focal_gene_state"],"UNMEASURED")
        self.assertFalse(record["ancestral_regain_demonstrated"])
        self.assertFalse(record["flower_specific_competence_demonstrated"])

    def test_short_or_frameshifted_candidate_not_gene_loss(self):
        self.assertEqual(cds_quality("AT")["sequence_qc_status"],"INCOMPLETE_OR_UNCERTAIN_CDS")
        complete=next(v for k,v in self.refs.items() if "|DFR|" in k)
        self.assertEqual(cds_quality(complete[:-1])["sequence_qc_status"],"INCOMPLETE_OR_UNCERTAIN_CDS")
        self.assertEqual(cds_quality(complete[:-4])["sequence_qc_status"],"INCOMPLETE_OR_UNCERTAIN_CDS")

    def test_synthetic_premature_stop_is_unvalidated_candidate(self):
        c=next(v for k,v in self.refs.items() if "|ANS|" in k)
        mutant=c[:30]+"TAA"+c[33:]
        self.assertEqual(len(mutant),len(c))
        v=cds_quality(mutant)
        self.assertEqual(v["sequence_qc_status"],"PREMATURE_STOP_CANDIDATE_REQUIRES_GENOMIC_VALIDATION")
        self.assertGreaterEqual(v["internal_stop_codon_count"],1)

    def test_many_ambiguous_bases_is_not_callable(self):
        c=next(v for k,v in self.refs.items() if "|DFR|" in k)
        mutant=c[:60]+"N"*45+c[105:]
        self.assertEqual(cds_quality(mutant)["sequence_qc_status"],"AMBIGUOUS_CDS_NOT_CALLABLE")

    def test_candidate_fasta_remains_sequence_qc_only(self):
        cds=next(v for k,v in self.refs.items() if "|DFR|" in k)
        with TemporaryDirectory() as folder:
            path=Path(folder)/"proposed.fa"
            path.write_text(">SRR35152725|DFR|synthetic_copy\n"+cds+"\n",encoding="utf-8")
            results=inspect_optional_candidates(path,self.r)
        self.assertEqual(len(results),1)
        self.assertEqual(results[0]["sequence_qc_status"],"OPEN_READING_FRAME_CANDIDATE_ONLY")
        self.assertFalse(results[0]["orthology_verified"])
        self.assertFalse(results[0]["floral_pathway_competence_established"])

    def test_multiple_candidate_paralogs_not_forced_to_single_gene(self):
        cds=next(v for k,v in self.refs.items() if "|DFR|" in k)
        with TemporaryDirectory() as folder:
            path=Path(folder)/"candidate.fa"
            path.write_text(">SRR35152725|DFR|contigA\n"+cds+
                            "\n>SRR35152725|DFR|contigB\n"+cds+"\n",encoding="utf-8")
            result=inspect_optional_candidates(path,self.r)
        self.assertEqual(len(result),2)
        self.assertTrue(all(v["multiple_paralogs_or_assembly_candidates"] for v in result))

    def test_unregistered_sample_or_wrong_locus_rejected(self):
        cds=next(v for k,v in self.refs.items() if "|DFR|" in k)
        for header in (">SRR000000000|DFR|fake",">SRR35152725|AN2|fake",
                       ">SRR35152725|DFR"):
            with TemporaryDirectory() as folder:
                path=Path(folder)/"candidate.fa"
                path.write_text(header+"\n"+cds+"\n",encoding="utf-8")
                with self.assertRaises(ValueError):inspect_optional_candidates(path,self.r)

    def test_fake_flower_expression_flag_rejected(self):
        with REGISTER.open(encoding="utf-8",newline="") as handle:
            rows=list(csv.DictReader(handle))
        rows[0]["floral_expression_measured"]="true"
        with TemporaryDirectory() as folder:
            dest=Path(folder)/"register.csv"
            with dest.open("w",encoding="utf-8",newline="") as handle:
                w=csv.DictWriter(handle,fieldnames=list(rows[0]))
                w.writeheader();w.writerows(rows)
            with self.assertRaises(ValueError):load_register(dest,RUNS)

    def test_gmyc1_reference_is_regulator_not_dfr(self):
        result=validate_references(FASTA,PROV)
        self.assertNotEqual(result["GMYC1"]["reference_accession"],result["DFR"]["reference_accession"])

    def test_invalid_fasta_rejected(self):
        with TemporaryDirectory() as folder:
            f=Path(folder)/"bad.fa"
            f.write_text(">duplicate\nATGTAA\n>duplicate\nATGTAA\n",encoding="utf-8")
            with self.assertRaises(ValueError):read_fasta(f)


if __name__=="__main__":
    unittest.main()
