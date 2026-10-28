"""Independent validation of exact ORIGINAL Petunieae 2023 source-body gates."""
import hashlib
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_petunieae_2023_processed_files_v26 import (
    PAPER_DOI,OSF_PROJECT,FILES,registered_location,verify_original,
    source_summary_from_bytes,preflight,
)


class PetunieaeOriginalProcessedFilesV26(unittest.TestCase):
    def test_exact_original_source_file_manifest_has_five_verifiable_IDs(self):
        self.assertEqual(PAPER_DOI,"10.1098/rspb.2023.0275")
        self.assertEqual(OSF_PROJECT,"zg9cu")
        self.assertEqual(len(FILES),5)
        self.assertEqual(FILES["expression_pigment_matched"]["size"],26917)
        self.assertEqual(FILES["expression_gene_mean"]["size"],30410)
        self.assertEqual(FILES["pigment_mean"]["size"],7039)
        self.assertEqual(FILES["species_tree_pruned"]["size"],1481)
        self.assertEqual(FILES["pigment_cluster_source"]["size"],551)
        self.assertEqual(len({z["file_id"] for z in FILES.values()}),5)
        for row in FILES.values():
            self.assertEqual(len(row["md5"]),32)
            self.assertEqual(len(row["sha256"]),64)
            self.assertTrue(registered_location(row).startswith("https://osf.io/download/"))

    def test_correct_original_body_authentication_is_not_filename_authentication(self):
        raw=b"study\tmeasured\nplant-001\t0.52\n"
        record={"size":len(raw),"md5":hashlib.md5(raw).hexdigest(),
                "sha256":hashlib.sha256(raw).hexdigest()}
        self.assertTrue(verify_original(raw,record))
        for content in (raw+b"more", b"<html>Validating...</html>", raw.replace(b"0.52",b"0.51")):
            with self.assertRaisesRegex(ValueError,"mismatch"):
                verify_original(content,record)
        with self.assertRaises(ValueError):
            verify_original(raw.decode(),record)

    def test_verified_csv_structure_still_not_a_validated_species_phenotype_join(self):
        data=b"species,CHS,CHI,flavonol\nSp_a,1.2,2.7,0.02\nSp_b,5,0,\n"
        d=source_summary_from_bytes("expression_pigment_matched",data)
        self.assertEqual(d["source_row_count"],2)
        self.assertEqual(d["source_column_count"],4)
        self.assertEqual(d["empty_cells"],1)
        self.assertEqual(d["distinct_first_column_values"],2)
        self.assertTrue(d["original_individual_to_species_grain_NOT_YET_ADMITTED"])
        self.assertTrue(d["no_new_heldout_or_causal_prediction_calculated"])

    def test_invalid_rectangular_table_and_duplicate_headers_rejected(self):
        for content in (b"species,CHS,CHS\nx,1,2\n",b"species,CHI\nx\n",
                        b"species,CHI\n",b"species,CHI\nx,1,wrong\n"):
            with self.assertRaises(ValueError):
                source_summary_from_bytes("pigment_mean",content)

    def test_tree_is_not_multiple_ancestral_regains(self):
        x=source_summary_from_bytes("species_tree_pruned",b"(Petunia_a:1,Petunia_b:1);\n")
        self.assertEqual(x["type"],"NEWICK_TREE")
        self.assertTrue(x["source_tip_labels_not_yet_reconciled"])
        with self.assertRaises(ValueError):
            source_summary_from_bytes("species_tree_pruned",b"fake_no_tree")

    def test_unavailable_original_processed_source_yields_no_biological_numbers(self):
        def unavailable(url):
            raise OSError("network unavailable in negative unit test")
        x=preflight(unavailable)
        self.assertFalse(x["complete_joint_expression_pigment_source_verified"])
        self.assertFalse(x["independent_species_holdout_analysis_PERFORMED"])
        self.assertFalse(x["full_functional_pathway_or_ancestral_regain_measured"])
        self.assertFalse(x["natural_selection_heritable_return_or_fitness_estimated"])
        self.assertEqual(x["Cirsium_empirical_owner"],"zuizui0223/aza3")
        self.assertTrue(all(not y["real_source_file_body_hash_verified"] for y in x["source_file_statuses"].values()))

    def test_fake_osf_body_with_real_original_source_metadata_is_not_admitted(self):
        def fake(url):
            return b"species,genotype,flavonol\nFAKE,0,1\n"
        x=preflight(fake)
        self.assertFalse(x["complete_joint_expression_pigment_source_verified"])
        self.assertTrue(all("SOURCE_TECHNICAL" in z["status"]
                            for z in x["source_file_statuses"].values()))
        self.assertTrue(x["source_rows_are_authors_processed_not_independent_RNA_HPLC_requantification"])


if __name__=="__main__":
    unittest.main()
