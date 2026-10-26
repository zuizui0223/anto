"""Source-locked audit for original Iochroma 2025 12 red/12 blue individuals.

Only original author transcript TPM and specimen phenotypic label counts:
NO experiments with loss-and-regain, no causal gene fine mapping, no phylo
replication, and no claim original paper discoveries are newly ours.
"""
import hashlib
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_iochroma_2025_osf_myb12_backcross import (
    SOURCE, MYB_ID, parse_manifest, audit_tpm_table, validate_payload,
    assert_expected_receipt, coexpression_of_selected_perfect,
    assert_expected_coexpression,
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

    def test_frozen_original_source_receipt_preserves_exact_separators(self):
        path=ROOT/"results/iochroma_2025_24_individual_myb12_16_transcript_source_audit_20261026.json"
        ref=json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(ref["original_n_transcript_entries"],19184)
        self.assertEqual(len(ref["our_descriptive_audit"]["strict_perfect_source_entries"]),16)
        self.assertEqual(ref["our_descriptive_audit"]["n_distinct_TRINITY_DN_assembly_component_labels"],14)
        fake={
            "n_samples_blue":12,"n_samples_red":12,
            "total_filtered_transcripts":19184,
            "perfect_blue_over_red":16,"perfect_red_over_blue":0,
            "perfect_blue_transcript_IDs":[x["id"] for x in ref["our_descriptive_audit"]["strict_perfect_source_entries"]],
            "perfect_red_transcript_IDs":[],
            "selected_gene":{
                "blue_mean_tpm":ref["our_descriptive_audit"]["mean_blue_TPM"],
                "red_mean_tpm":ref["our_descriptive_audit"]["mean_red_TPM"],
                "blue_min_tpm":ref["our_descriptive_audit"]["min_blue_TPM"],
                "red_max_tpm":ref["our_descriptive_audit"]["max_red_TPM"],
                "pairwise_blue_over_red_AUC":1.,
                "individuals":[{"source_original_ID":v["sample_id"],
                                "paper_phenotype":"blue" if v["source_flower_colour"]=="B" else "red",
                                "TPM":v["TPM"]}
                               for v in ref["our_descriptive_audit"]["original_Myb12_TPM_by_sample"]],
            },
        }
        assert_expected_receipt(fake,ref)
        perturbed=copy.deepcopy(fake)
        perturbed["perfect_blue_transcript_IDs"][0]="SYNTHETIC_ID"
        with self.assertRaisesRegex(ValueError,"receipt"):
            assert_expected_receipt(perturbed,ref)
        perturbed=copy.deepcopy(fake)
        perturbed["selected_gene"]["individuals"][0]["TPM"]+=1
        with self.assertRaisesRegex(ValueError,"receipt"):
            assert_expected_receipt(perturbed,ref)

    def test_pooled_correlation_is_not_within_color_coregulation(self):
        assigned={"B1":"B","B2":"B","R1":"R","R2":"R"}
        ids=["B1","R1","B2","R2"]
        vectors={
            "MYB":[8,0,9,1],
            "ANTI":[20,3,10,1],
            "SAME":[12,0,14,1],
        }
        diagnostic=coexpression_of_selected_perfect(vectors,"MYB",ids,assigned)
        self.assertEqual(diagnostic["other_blue_perfect_transcripts"],2)
        self.assertEqual(diagnostic["pooled_positive"],2)
        self.assertEqual(diagnostic["color_centered_negative"],1)
        self.assertEqual(diagnostic["color_centered_positive"],1)
        self.assertEqual(diagnostic["color_centered_undefined"],0)
        rows={v["transcript_id"]:v for v in diagnostic["per_transcript"]}
        self.assertGreater(rows["ANTI"]["pooled_r"],0)
        self.assertLess(rows["ANTI"]["within_color_centered_r"],0)
        self.assertGreater(rows["SAME"]["within_color_centered_r"],0)

    def test_centered_correlation_handles_constant_within_groups_without_false_sign(self):
        assigned={"B1":"B","B2":"B","R1":"R","R2":"R"}
        ids=["B1","B2","R1","R2"]
        vectors={
            "MYB":[8,9,0,1],
            "STRONG_BUT_FLAT":[15,15,0,0],
        }
        diagnostic=coexpression_of_selected_perfect(vectors,"MYB",ids,assigned)
        self.assertEqual(diagnostic["pooled_positive"],1)
        self.assertEqual(diagnostic["color_centered_undefined"],1)
        self.assertIsNone(diagnostic["per_transcript"][0]["within_color_centered_r"])

    def test_original_24_individual_15_comparator_frozen_source_correlations(self):
        receipt=json.loads((ROOT/"results/iochroma_2025_color_conditioned_coexpression_v2_20261026.json")
                           .read_text(encoding="utf-8"))
        d=receipt["source_verified_descriptive_result"]
        self.assertEqual(d["original_individuals"],24)
        self.assertEqual((d["blue_individuals"],d["red_individuals"]),(12,12))
        self.assertEqual(d["reference_transcript"],MYB_ID)
        self.assertEqual(d["other_blue_perfect_transcripts"],15)
        self.assertEqual(d["pooled_positive"],15)
        self.assertEqual(d["pooled_negative"],0)
        self.assertEqual(d["color_centered_negative"],9)
        self.assertEqual(d["color_centered_positive"],6)
        self.assertEqual(d["color_centered_undefined"],0)
        self.assertEqual(len({x["transcript_id"] for x in d["per_transcript"]}),15)
        self.assertTrue(all(x["pooled_r"]>0 for x in d["per_transcript"]))
        self.assertEqual(sum(x["within_color_centered_r"]<0
                             for x in d["per_transcript"]),9)
        self.assertTrue(receipt["no_multiple_testing_corrected_discovery_or_causal_coregulation"])
        self.assertTrue(receipt["not_ancestral_trait_regain_fitness_or_multigeneration_retention"])
        self.assertTrue(receipt["do_not_merge_transcript_rows_into_independent_genetic_loci"])
        assert_expected_coexpression(d,receipt)

        # Counts and all 15 actual transcript coefficients are source-locked.
        for mode in ("count","positive","negative","id","rho","source_digest"):
            bad=copy.deepcopy(receipt)
            entries=bad["source_verified_descriptive_result"]
            if mode=="count": entries["color_centered_negative"]=8
            if mode=="positive": entries["pooled_positive"]=14
            if mode=="negative": entries["color_centered_positive"]=5
            if mode=="id": entries["per_transcript"][0]["transcript_id"]="SYNTHETIC"
            if mode=="rho": entries["per_transcript"][0]["pooled_r"]+=.01
            if mode=="source_digest": bad["source_original_tpm_md5"]="WRONG"
            with self.assertRaisesRegex(ValueError,"Frozen"):
                assert_expected_coexpression(d,bad)

    def test_source_byte_digest_or_count_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"source"):
            validate_payload(b"not the original bytes",SOURCE["sample"])
        with self.assertRaisesRegex(ValueError,"source"):
            validate_payload(b"",SOURCE["tpm"])

if __name__=="__main__":
    unittest.main()
