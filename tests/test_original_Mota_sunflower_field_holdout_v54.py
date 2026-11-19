"""ANTO v54 source-locked sunflower physical-field split and no-causality tests."""
import copy
import json
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mota_sunflower_field_holdout_v54 as v


class OriginalFieldYearSunflowerBenchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reg=v.read_registry()
        cls.parent=v.verified_anchor(cls.reg)

    def test_original_one_study_three_md5_and_80_fieldyears(self):
        r=self.reg
        self.assertEqual(r["zenodo_record_id"],6671199)
        self.assertEqual(len(r["original_files_md5"]),3)
        self.assertEqual(r["original_matched_fieldyears"],80)
        self.assertEqual(r["original_matched_physical_fields"],51)
        self.assertEqual(self.parent["n_shared_field_year_units"],80)
        self.assertEqual(self.parent["physical_field_metadata_join_audit"]
                         ["shared_fieldyear_same_original_physical_field_name"],80)

    def test_original_51_fields_not_80_independent_fields(self):
        self.assertEqual(self.reg["original_matched_fieldyears_by_year"],
                         {"2017":41,"2018":39})
        self.assertEqual(self.reg["original_matched_fieldyears_by_region"],
                         {"Burgos":38,"Cuenca":42})
        self.assertEqual(self.reg["original_matched_fieldyears_by_field_type"],
                         {"NonV":10,"SNV":32,"WFS":38})

    def test_no_same_plant_causal_or_genetic_claim(self):
        self.assertFalse(self.reg["original_individual_sunflower_visitation_linked_to_harvest"])
        self.assertTrue(self.reg["study_is_posthoc"])
        self.assertFalse(self.reg["valid_randomization_p_value"])
        self.assertEqual(self.reg["new_original_independent_studies_vs_v53_source_gate"],0)

    def test_registered_three_models_and_original_ridge(self):
        self.assertEqual(v.MODELS,("B0","B1","B2"))
        self.assertEqual(v.RADIUS,8)
        self.assertEqual(v.BANDS,("0","15","30","60"))
        self.assertEqual(self.reg["permutations"],200)

    def test_original_missing_does_not_become_zero(self):
        self.assertIsNone(v.parse_number("NA"))
        self.assertIsNone(v.parse_number(""))
        self.assertEqual(v.parse_number("0"),0)
        with self.assertRaisesRegex(v.V54Error,"INVALID_COUNT"):
            v.parse_number("-2")
        with self.assertRaisesRegex(v.V54Error,"UNKNOWN_SOURCE_NUMERIC"):
            v.parse_number("??")

    def test_seed_accounting_source_semantic_fixture(self):
        row={"Year":"2017","Field_code":"fA","Region":"Burgos",
             "Field_Type":"SNV","Field_code_original":"physical1",
             "Distance":"15"}
        self.assertEqual(v.fieldyear_key(row),("2017","fA"))
        self.assertEqual(v.distance_key(row),"15")
        self.assertEqual(v.sanitize_metadata(row,"harvest"),
                         ("Burgos","SNV","physical1"))
        row["Distance"]="17"
        with self.assertRaisesRegex(v.V54Error,"SOURCE_DISTANCE"):
            v.distance_key(row)

    @staticmethod
    def synthetic(n=65):
        rows=[]
        for i in range(n):
            row={"year":"2017" if i%2==0 else "2018",
                 "region":"Burgos" if i%3 else "Cuenca",
                 "field_type":v.TYPES[i%3],
                 "physical_field":("Burgos","p"+str(i)),
                 "physical_field_sha256":str(i),
                 "n_common_valid_distance_bands":4,
                 "seed_fraction":.10+.7*(i%17)/17,
                 "log_total_visit_rate":.1+.05*(i%23),
                 "log_apis_visit_rate":.07+.03*(i%19),
                 "log_wild_visit_rate":.05+.02*(i%13)}
            rows.append(row)
        return rows

    def test_fit_does_not_need_heldout_seed_fraction(self):
        data=self.synthetic()
        train=data[1:]
        test=dict(data[0])
        for name in v.MODELS:
            model=v.fit(train,name)
            original=v.predict(model,test)
            test["seed_fraction"]=0
            self.assertEqual(original,v.predict(model,test))
            test["seed_fraction"]=.9
            self.assertEqual(original,v.predict(model,test))

    def test_fixed_models_dimensions(self):
        data=self.synthetic()
        self.assertEqual(len(v.fit(data,"B0")["coef"]),5)
        self.assertEqual(len(v.fit(data,"B1")["coef"]),6)
        self.assertEqual(len(v.fit(data,"B2")["coef"]),7)
        for m in v.MODELS:
            self.assertTrue(0<=v.predict(v.fit(data,m),data[0])[0]<=1)

    def test_original_two_year_same_field_is_heldout_together(self):
        data=self.synthetic()
        copyrow=dict(data[0])
        copyrow["year"]="2018"
        copyrow["seed_fraction"]=0
        data.append(copyrow)
        folds=v.leave_field_out(data)
        field=next(x for x in folds if x["physical_field_sha256"]=="0")
        self.assertEqual(field["n_fieldyears"],2)
        self.assertEqual(field["n_training_other_physical_fields"],64)
        self.assertEqual(field["heldout_fieldyear_years"],["2017","2018"])

    def test_prediction_saturation_metadata_retained(self):
        data=self.synthetic(75)
        scored=v.leave_field_out(data)
        s=v.score_summary(scored)
        self.assertEqual(s["original_heldout_physical_fields"],75)
        self.assertEqual(s["original_heldout_fieldyears"],75)
        self.assertIn("total_raw_below_zero",s["models"]["B1"])

    def test_incorrect_training_source_denominator_blocks(self):
        with self.assertRaisesRegex(v.V54Error,"INVALID_TRAINING"):
            v.fit(self.synthetic()[:5],"B1")
        with self.assertRaisesRegex(v.V54Error,"INVALID_TRAINING"):
            v.fit(self.synthetic(),"invented")

    def test_same_original_source_aggregate_adversarial_group(self):
        x=self.synthetic()
        x[1]["physical_field"]=x[0]["physical_field"]
        x[1]["physical_field_sha256"]=x[0]["physical_field_sha256"]
        # In this fixture, one same field appears twice and must be
        # contained in a single heldout prediction (never divided).
        folds=v.leave_field_out(x)
        field=next(z for z in folds if z["physical_field_sha256"]=="0")
        self.assertEqual(field["n_fieldyears"],2)

    def test_original_source_full_receipt_when_available(self):
        if not v.RECEIPT.exists():
            self.skipTest("Only current MD5 original source completed receipt")
        report=json.loads(v.RECEIPT.read_text())
        self.assertFalse(report["fitted_individual_sunflower_visitor_to_harvest_effect"])
        self.assertFalse(report["causal_visitor_intervention_or_genetic_inheritance_estimated"])
        self.assertEqual(report["source_audit"]["source_exact_files_md5"],
                         self.reg["original_files_md5"])
        self.assertEqual(report["original_v53_source_fieldyear_receipt_git_blob_sha1"],
                         self.reg["prior_original_v53_fieldyear_receipt_git_blob_sha1"])
        self.assertEqual(report["within_original_study_visit_identity_placebo"]
                         ["replicate_count"],200)
        self.assertFalse(report["independent_randomization_p_value_calculated"])


if __name__=="__main__":
    unittest.main()
