"""Original Mota v55: strict source-grain, no invented bee-sum, holdout tests."""
import json
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mota_reported_total_bees_v55 as v

class SourceVerifiedSunflowerTotal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reg=v.registry()
        cls.prior=v.require_prior_quality(cls.reg)

    def test_immutable_prior_quality_and_registrations(self):
        r=self.reg
        self.assertEqual(v.git_blob(v.QUALITY.read_bytes()),
                         "55eb9b942e7a0ab7f4e32139ac5305325fb74a83")
        self.assertEqual(r["original_v53_fieldyear_blob_sha1"],
                         "7a4fb27c7f26abfa8400d23b041e1a3d3e3424c9")
        self.assertEqual(r["n_source_missing_aborts"],1056)
        self.assertEqual(r["n_source_visits_total_not_apis_plus_wild"],4087)
        self.assertEqual(r["n_original_physical_fields"],51)

    def test_original_v54_stop_does_not_migrate_to_biological_null(self):
        self.assertFalse(self.prior["predictive_effect_or_visitor_seed_association_computed"])
        self.assertFalse(self.prior["original_v54_preregistered_source_identity_guard_modified"])
        self.assertTrue(self.reg["v54_source_stop_retained"])
        self.assertFalse(self.reg["source_model_fitted_before_protocol"])
        self.assertFalse(self.reg["visitor_group_component_counts_assumed_equal_total"])

    def test_no_same_individual_or_causal_inference(self):
        self.assertFalse(self.reg["same_individual_plant_visit_and_harvest"])
        self.assertFalse(self.reg["causal_randomization_p_value"])
        self.assertEqual(self.reg["new_original_biological_studies_vs_v54"],0)

    def test_three_original_bounded_prediction_models(self):
        self.assertEqual(v.MODELS,("Z0","Z1","Z2"))
        self.assertEqual(v.LAMBDA,8.0)
        self.assertEqual(v.BANDS,("0","15","30","60"))

    @staticmethod
    def fixture():
        rows=[]
        for i in range(80):
            idx=i%51
            region="Burgos" if idx%2==0 else "Cuenca"
            kind=v.TYPES[idx%3]
            field=(region,"plot"+str(idx))
            rows.append({
                "year":"2017" if i<51 else "2018",
                "region":region,
                "field_type":kind,
                "physical_field":field,
                "physical_field_hash":"hash"+str(idx),
                "seed_fraction":.20+.45*(i%13)/13,
                "log_reported_total_rate":.02+.03*(i%15),
                "log_reconstituted_rate":.03+.04*(i%11)
            })
        return rows

    def test_test_seed_outcome_does_not_appear_in_prediction(self):
        rows=self.fixture()
        train=rows[1:]
        test=dict(rows[0])
        for k in v.MODELS:
            f=v.fit(train,k)
            val=v.predict(f,test)
            test["seed_fraction"]=1.0
            self.assertEqual(val,v.predict(f,test))
            test["seed_fraction"]=0.0

    def test_exact_model_dimensions_and_bounded_predictions(self):
        rows=self.fixture()
        self.assertEqual(len(v.fit(rows,"Z0")["beta"]),5)
        self.assertEqual(len(v.fit(rows,"Z1")["beta"]),6)
        self.assertEqual(len(v.fit(rows,"Z2")["beta"]),6)
        for m in v.MODELS:
            val,raw=v.predict(v.fit(rows,m),rows[0])
            self.assertTrue(0<=val<=1)
            self.assertTrue(math.isfinite(raw))

    def test_full_two_year_physical_field_is_always_held_out_together(self):
        folds=v.holdouts(self.fixture())
        self.assertEqual(len(folds),51)
        self.assertEqual(sum(f["n_original_fieldyears"] for f in folds),80)
        twice=[f for f in folds if f["n_original_fieldyears"]==2]
        self.assertEqual(len(twice),29)
        self.assertEqual(v.summary(folds)["n_original_fieldyear_predictions"],80)
        self.assertEqual(set(f["region"] for f in folds),set(v.REGIONS))

    def test_inference_does_not_use_numeric_visitor_counts_as_plant_replicates(self):
        rows=self.fixture()
        out=v.summary(v.holdouts(rows))
        self.assertEqual(out["n_original_physical_field_holdouts"],51)
        for x in ("Z0","Z1","Z2"):
            self.assertIn("equal_physical_field_mse",out["model_risks"][x])
            self.assertIn("raw_predictions_below_zero",out["model_risks"][x])

    def test_unknown_model_and_bad_train_n_fails(self):
        with self.assertRaisesRegex(v.V55Error,"MODEL_OR_TRAIN_N"):
            v.fit(self.fixture(),"nonexistent")
        with self.assertRaisesRegex(v.V55Error,"MODEL_OR_TRAIN_N"):
            v.fit(self.fixture()[:2],"Z1")

    def test_source_eligibility_is_still_80_fieldyears_51_fields(self):
        self.assertEqual(self.reg["n_expected_original_fieldyear_predictions"],80)
        self.assertEqual(self.reg["n_expected_original_field_holdouts"],51)
        self.assertEqual(self.reg["source_original_fieldyear_band_counts"],
                         {"3":1,"4":79})
        self.assertEqual(self.reg["n_source_visit_rows_matched"],11682)
        self.assertEqual(self.reg["n_source_harvest_rows_matched"],2217)

    def test_full_result_if_frozen(self):
        if not v.RECEIPT.exists():
            self.skipTest("Only source-complete result available after MD5 first run")
        output=json.loads(v.RECEIPT.read_text())
        self.assertEqual(output["summary"]["n_original_physical_field_holdouts"],51)
        self.assertEqual(output["summary"]["n_original_fieldyear_predictions"],80)
        self.assertEqual(output["source_audit"]["original_abort_source_rows_missing_quarantined_not_imputed"],1056)
        self.assertEqual(output["source_audit"]["original_apis_plus_wild_vs_reported_total_mismatch_count"],4087)
        self.assertTrue(output["does_not_link_same_sunflower_plant_visitor_and_seed"])
        self.assertTrue(output["no_valid_randomization_p_value_computed"])
        self.assertEqual(output["visitor_source_label_placebo"]["n_original_fieldyear_stratum_permutations"],200)

if __name__=="__main__":
    unittest.main()
