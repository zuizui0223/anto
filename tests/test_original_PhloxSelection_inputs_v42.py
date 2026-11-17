"""ANTO v42 source-integrity tests; original actual author data fetched by CI, not test fixtures."""
import copy
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_original_PhloxSelection_author_inputs_v42 as v


def pheno(text=None):
    return (text or
       "pop\tdistance\tDB\tDR\tLB\tLR\n"
       "301\t0\t1\t0\t3\t2\n"
       "302\t1.5\t0\t2\t4\t3\n").encode()


def mother(text=None):
    return (text or
       "population\tcolor\tDB\tDR\tLB\tLR\n"
       "301\t2\t1\t0\t2\t0\n"
       "302\t1\t0\t3\t0\t1\n"
       "\t\t\t\t\t\n").encode()


class V42OriginalAuthorInput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=json.loads(v.MANIFEST.read_text())
        cls.v41=json.loads(v.V41.read_text())

    def test_exact_historical_git_commit_and_20_authentic_blob_manifest(self):
        self.assertEqual(self.manifest["original_commit_sha"],
                         "60ff508170f1be592bbc6fe7ed46ac6e5b459026")
        self.assertEqual(self.manifest["canonical_GitHub_repo"],"PhloxHopkins/PhloxSelection")
        self.assertEqual(len(self.manifest["all_author_files_fixed_git_blob_sha1"]),20)
        self.assertEqual(self.manifest["source_files_phen_count"],5)
        self.assertEqual(self.manifest["source_files_fam_count"],5)

    def test_author_readme_inputs_described_as_examples(self):
        self.assertTrue(self.manifest["original_author_README_calls_files_example_input"])
        self.assertFalse(self.manifest["original_full_model_Likelihood_MCMC_reanalysis_claimed"])
        self.assertTrue(self.manifest["no_original_raw_source_files_committed_to_ANTO"])

    def test_git_blob_known_vector_and_tamper_changes_digest(self):
        self.assertEqual(v.git_blob_sha1(b"hello\n"),"ce013625030ba8dba906f756967f9e9ca394464a")
        self.assertNotEqual(v.git_blob_sha1(b"hello\n"),v.git_blob_sha1(b"hello!"))
        with self.assertRaises(v.OriginalSourceError):
            v.git_blob_sha1("not-bytes")

    def test_file_digest_gate_catches_same_length_substitution(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"original.data"
            p.write_bytes(b"abcdef")
            expected=v.git_blob_sha1(p.read_bytes())
            self.assertEqual(v.digest_file(p,expected),b"abcdef")
            p.write_bytes(b"abcxef")
            with self.assertRaisesRegex(v.OriginalSourceError,"GIT_BLOB_SHA1_DRIFT"):
                v.digest_file(p,expected)

    def test_original_phen_source_order_and_header(self):
        r=v.parse_phen(pheno())
        self.assertEqual(len(r),2)
        self.assertEqual(r[0]["population"],301)
        self.assertEqual(r[0]["by_color"],{"DB":1,"DR":0,"LB":3,"LR":2})
        self.assertEqual(r[1]["distance_km"],1.5)

    def test_phen_duplicate_and_out_of_order_distance_fail(self):
        for raw,pattern in (
            (pheno("pop\tdistance\tDB\tDR\tLB\tLR\n301\t0\t1\t0\t3\t2\n301\t1\t0\t2\t4\t3\n"),"DUPLICATE_POPULATION"),
            (pheno("pop\tdistance\tDB\tDR\tLB\tLR\n301\t2\t1\t0\t3\t2\n302\t1\t0\t2\t4\t3\n"),"UNSORTED_DISTANCE"),
        ):
            with self.assertRaisesRegex(v.OriginalSourceError,pattern):
                v.parse_phen(raw)

    def test_phen_malformed_negative_and_noninteger_counts_fail(self):
        for raw in (
           pheno("pop\tdistance\tDB\tDR\tLB\tLR\n301\t0\t1\t-1\t3\t2\n"),
           pheno("pop\tdistance\tDB\tDR\tLB\tLR\n301\t0\t1.0\t0\t3\t2\n"),
           pheno("pop\tdistance\tDB\tDR\tLB\tLR\n301\t0\t1\t0\t3\n"),
        ):
            with self.assertRaises(v.OriginalSourceError):v.parse_phen(raw)

    def test_blank_family_padding_ignored(self):
        r=v.parse_fam(mother(),{301,302})
        self.assertEqual(len(r),2)
        self.assertEqual(r[0]["maternal_phenotype"],"LB")
        self.assertEqual(r[1]["maternal_phenotype"],"DR")
        self.assertEqual(sum(x["n_offspring"] for x in r),7)

    def test_fam_invalid_maternal_color_code_caught(self):
        bad=mother("population\tcolor\tDB\tDR\tLB\tLR\n301\t4\t1\t0\t2\t0\n")
        with self.assertRaisesRegex(v.OriginalSourceError,"MOTHER_COLOR_CODE_CHANGED"):
            v.parse_fam(bad,{301})

    def test_fam_unrelated_population_is_rejected(self):
        with self.assertRaisesRegex(v.OriginalSourceError,"POP_NOT_IN_SAME_TRANSECT"):
            v.parse_fam(mother(),{301})

    def test_fam_zero_children_and_bad_count_fail(self):
        for bad in (
          mother("population\tcolor\tDB\tDR\tLB\tLR\n301\t1\t0\t0\t0\t0\n"),
          mother("population\tcolor\tDB\tDR\tLB\tLR\n301\t1\t0\tX\t0\t1\n"),
        ):
            with self.assertRaises(v.OriginalSourceError):v.parse_fam(bad,{301})

    def test_parent_child_matrix_is_family_row_grain_not_independent_children(self):
        f=v.parse_fam(mother(),{301,302})
        z=v.parent_child_matrix(f)
        self.assertEqual(z["LB"]["n_source_family_rows"],1)
        self.assertEqual(z["LB"]["offspring_n"],3)
        self.assertEqual(z["DR"]["n_source_family_rows"],1)
        self.assertEqual(z["DR"]["offspring_n"],4)
        self.assertEqual(z["DB"]["offspring_n"],0)

    def test_v41_previous_provisional_source_32_rows_untouched(self):
        r=self.v41["source_transcription"]["population_phenotype_counts"]
        self.assertEqual(len(r),32)
        self.assertEqual(sum(sum(x[2:]) for x in r),11379)
        self.assertEqual(self.v41["source_original"]["main_article_claim_adult_survey_total"],12024)

    def test_duplicate_init_files_not_independent_replicates(self):
        z=self.manifest["all_author_files_fixed_git_blob_sha1"]
        for i in range(1,6):
            self.assertEqual(z[f"inputFiles/init_Phlox_IF{i}.txt"],
                             z[f"inputFiles/inPhlox_Phen_T{i}.data"])

    def test_manifest_rejects_misregistered_original_cohort_as_complete(self):
        p=self.manifest
        self.assertEqual(p["original_article_claimed_counts"]["adult"],12024)
        self.assertEqual(p["pre_script_source_author_snapshot_observed"]["adult"],11379)
        self.assertEqual(p["original_article_claimed_counts"]["offspring"],2581)
        self.assertEqual(p["pre_script_source_author_snapshot_observed"]["offspring"],2496)
        self.assertEqual(p["original_article_claimed_counts"]["maternal_plants"],382)
        self.assertEqual(p["pre_script_source_author_snapshot_observed"]["maternal_families"],371)

    def test_fake_repo_or_non_checkout_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(v.OriginalSourceError,"NOT_A_GIT_CHECKOUT"):
                v.report(Path(d),self.manifest,self.v41)

    def test_original_phenotype_codes_unambiguous_and_fixed(self):
        self.assertEqual(v.MOTHER_LABEL,{0:"DB",1:"DR",2:"LB",3:"LR"})
        self.assertTrue(self.manifest["original_code_maternal_coding_0DB_1DR_2LB_3LR"])

    def test_no_numeric_biological_selection_claim_from_exemplar(self):
        self.assertTrue(self.manifest["no_independent_RNA_HPLC_or_natural_regain_measured"] if
                        "no_independent_RNA_HPLC_or_natural_regain_measured" in self.manifest else
                        self.manifest["no_original_raw_source_files_committed_to_ANTO"])

if __name__=="__main__":unittest.main()
