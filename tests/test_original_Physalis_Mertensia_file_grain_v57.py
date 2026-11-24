"""V57 original Physalis/Mertensia author-byte and plant key negative gates."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_original_Physalis_Mertensia_file_grain_v57 as f

class V57OriginalCohort(unittest.TestCase):
    def test_author_source_pin_is_a_real_committed_git_blob(self):
        r,p=f.check_gate()
        self.assertEqual(f.gitblob(f.PIN.read_bytes()),f.PIN_GIT_BLOB)
        self.assertEqual([len(x["author_files"]) for x in p["source_study_versions"]],[10,2])
        self.assertEqual([x["published_version"] for x in p["source_study_versions"]],[5,4])
        self.assertTrue(r["no_ecological_effect_fit"])

    def test_exact_official_paths_not_guessed_community_mirrors(self):
        self.assertEqual(f.exact_official_download("/api/v2/files/4573992/download"),
                         "https://datadryad.org/api/v2/files/4573992/download")
        with self.assertRaisesRegex(f.V57FileError,"NOT_IN_ORIGINAL_DRYAD_API"):
            f.exact_official_download("http://evil.example/4573992")

    def test_original_csv_source_plant_identifiers_only(self):
        x=f.rows_csv(b"block,treat,id,number_visitor,fruit_total\nA,C,1,57,800\nB,H,2,33,0\n")
        self.assertEqual(x["rows"],2)
        self.assertEqual(x["nominal_key_cardinality"],{"block":2,"treat":2,"id":2})
        self.assertNotIn("number_visitor",x["nominal_key_cardinality"])
        self.assertNotIn("fruit_total",x["nominal_key_cardinality"])

    def test_sameplant_composite_key_not_plant_number_only(self):
        a=f.rows_csv(b"block,treat,id,total_flow\nA,C,1,22\nB,C,1,18\n")
        b=f.rows_csv(b"block,treat,id,harvest_yield\nA,C,1,20\nB,C,2,21\n")
        d=f.cohort_2020(a,b)
        self.assertEqual(d["original_same_plant_key_overlap"],1)
        self.assertEqual(d["original_visitation_unique_plantkeys"],2)
        self.assertFalse(d["linked_plant_reproduction_or_visit_values_analysed"])

    def test_repeat_observations_not_distinct_2021_plants(self):
        a=f.rows_csv(b"block,treatment,plant_id,morphospecies,fruit_set\nA,c,1,Bombus,1\nA,c,1,Apis,0\nB,h,1,Bombus,1\n")
        r=f.cohort_2021(a)
        self.assertEqual(r["original_single_visit_rows"],3)
        self.assertEqual(r["n_distinct_original_plants_by_block_treatment_id"],2)
        self.assertEqual(r["repeated_original_plant_visit_records"],1)

    def test_illegal_fruit_outcome_as_visit_predictor_forbidden(self):
        import json
        r=json.loads((ROOT/"data/external/anto_v57_original_thirdstudy_source_registry_20261124.json").read_text())
        for name in ("poll_imp","prop_fs","avg_fruit_weight","fruit_set","total_yield_grams"):
            self.assertIn(name,r["physalis"]["forbidden_outcome_leaking_predictors"])

    def test_provider_sha_verified_not_decoration(self):
        source=b"abc"
        ok=f.original_digest_check("x.csv",source,{"digestType":"md5","digest":"900150983cd24fb0d6963f7d28e17f72"})
        self.assertEqual(ok["digest_comparison_status"],"provider_md5_verified")
        with self.assertRaisesRegex(f.V57FileError,"ORIGINAL_MD5_MISMATCH"):
            f.original_digest_check("x.csv",source,{"digestType":"md5","digest":"a00150983cd24fb0d6963f7d28e17f72"})

    def test_xlsx_invalid_original_zip_blocks(self):
        with self.assertRaises(Exception):
            f.inspect_xlsx(b"bad xlsx not original")
        
    def test_scope_original_site_not_six_populations(self):
        r,p=f.check_gate()
        self.assertEqual(r["minimum_original_sameplant_keys"],30)
        self.assertEqual(r["minimum_original_blocks"],3)
        self.assertEqual(r["new_independent_original_study_results_admitted"],0)

if __name__=="__main__":unittest.main()
