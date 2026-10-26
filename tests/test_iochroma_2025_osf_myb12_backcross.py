"""Source-locked audit for original Iochroma 2025 12 red/12 blue individuals.

Only original author transcript TPM and specimen phenotypic label counts:
NO experiments with loss-and-regain, no causal gene fine mapping, no phylo
replication, and no claim original paper discoveries are newly ours.
"""
import hashlib
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_iochroma_2025_osf_myb12_backcross import (
    SOURCE, MYB_ID, parse_manifest, audit_tpm_table, validate_payload,
)

class IochromaSourceTests(unittest.TestCase):
    def test_fixed_osf_original_file_ids(self):
        self.assertEqual(SOURCE["tpm"]["id"],"6622e385ea2b630b23532fe5")
        self.assertEqual(SOURCE["sample"]["id"],"6622e067387f9d08bf7ac834")
        self.assertEqual(SOURCE["tpm"]["bytes"],4792342)
        self.assertEqual(SOURCE["sample"]["bytes"],955)
        self.assertEqual(MYB_ID,"TRINITY_DN3000_c0_g1_i10.p1")

    def test_small_tpm_table_can_find_exact_individual_level_separator(self):
        manifest="ID,condition,path,location\nB1,B,,\nB2,B,,\nR1,R,,\nR2,R,,\n"
        mat=("gene_id,R2,B1,R1,B2\n"
             "MYB,0.1,8,0.3,9\n"
             "OTHER,3,2,4,6\n"
             "RED_HIGH,9,0,8,0\n")
        x=audit_tpm_table(mat,parse_manifest(manifest),selected_gene="MYB")
        self.assertEqual((x["n_samples_blue"],x["n_samples_red"]),(2,2))
        self.assertEqual(x["total_filtered_transcripts"],3)
        self.assertEqual(x["perfect_blue_over_red"],1)
        self.assertEqual(x["perfect_red_over_blue"],1)
        self.assertEqual(x["selected_gene"]["blue_min_tpm"],8)
        self.assertEqual(x["selected_gene"]["red_max_tpm"],.3)
        self.assertTrue(x["selected_gene"]["strict_rank_discrimination"])
        self.assertEqual(x["selected_gene"]["pairwise_blue_over_red_AUC"],1.0)

    def test_frozen_study_not_white_morph_vs_fixed_white_case(self):
        manifest="ID,condition,path,location\nB1,B,,\nR1,R,,\n"
        mat="gene_id,R1,B1\nMYB,0.1,4.0\n"
        x=audit_tpm_table(mat,parse_manifest(manifest),selected_gene="MYB")
        self.assertFalse(x["trait_loss_then_historical_regain_measured"])
        self.assertFalse(x["causal_MYB_gene_demonstrated_by_reanalysis"])
        self.assertFalse(x["white_morph_vs_fixed_white_species_measured"])
        self.assertTrue(x["single_original_backcross_not_independent_macro_clades"])

    def test_no_fake_perfect_separation_from_zero_or_ties(self):
        manifest="ID,condition,path,location\nB1,B,,\nB2,B,,\nR1,R,,\nR2,R,,\n"
        mat="gene_id,B1,R1,B2,R2\nMYB,2,0,2,2\nSAME,0,0,0,0\n"
        z=audit_tpm_table(mat,parse_manifest(manifest),selected_gene="MYB")
        self.assertEqual(z["perfect_blue_over_red"],0)
        self.assertEqual(z["perfect_red_over_blue"],0)
        self.assertEqual(z["selected_gene"]["pairwise_blue_over_red_AUC"],.75)

    def test_duplicate_sample_or_missing_assignment_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"duplicate"):
            parse_manifest("ID,condition,path,location\nB1,B,,\nB1,R,,\n")
        x=parse_manifest("ID,condition,path,location\nB1,B,,\nR1,R,,\n")
        with self.assertRaisesRegex(ValueError,"sample"):
            audit_tpm_table("gene_id,B1,XXX\nMYB,1,0\n",x,selected_gene="MYB")

    def test_negative_tpm_missing_gene_and_duplicate_transcript_fail_closed(self):
        m=parse_manifest("ID,condition,path,location\nB1,B,,\nR1,R,,\n")
        for mat,word in (
            ("gene_id,B1,R1\nMYB,-1,0\n","TPM"),
            ("gene_id,B1,R1\nX,1,0\n","MYB"),
            ("gene_id,B1,R1\nMYB,1,0\nMYB,1,0\n","duplicate"),
        ):
            with self.assertRaisesRegex(ValueError,word):
                audit_tpm_table(mat,m,selected_gene="MYB")

    def test_source_byte_digest_or_count_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"source"):
            validate_payload(b"not the original bytes",SOURCE["sample"])
        with self.assertRaisesRegex(ValueError,"source"):
            validate_payload(b"",SOURCE["tpm"])

if __name__=="__main__":
    unittest.main()
