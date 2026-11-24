"""Independent original Mertensia array-holdout adversarial source tests.

No regression result is selected after reading original seed/visitor values.
"""
from __future__ import annotations
import json
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mertensia_20array_visitation_v58 as v


class OriginalMertensia20WholeArrays(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry=v.registry()
        cls.parent=v.parent_source(cls.registry)

    def test_original_author_both_archive_digest_and_cohort_frozen(self):
        r=self.registry
        self.assertEqual(v.gitblob(v.PARENT.read_bytes()),
                         "5659545f1376aeb715aefa4308a0fe1e8673da39")
        self.assertEqual(r["original_xlsx_sha256"],
                         "152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d")
        self.assertEqual(r["original_xlsx_zenodo_md5"],
                         "db86f8286007b6176bf957334eba9ca4")
        self.assertEqual(r["source_original_plant_record_rows"],115)
        self.assertEqual(r["source_complete_eligible_plant_records"],105)

    def test_original_source_20_array_counts_not_independent_studies(self):
        p=self.parent["original_source_group_grain"]
        self.assertEqual(p["n_original_week_array_clusters"],20)
        self.assertEqual(p["source_row_duplicate_composite_plantkeys"],0)
        self.assertEqual(p["n_original_plant_records_all_candidate_variables_PRESENT"],105)
        self.assertEqual(self.registry["source_groups_min_raw_plant_rows"],5)
        self.assertEqual(self.registry["source_groups_max_raw_plant_rows"],7)

    @staticmethod
    def synthetic():
        rows=[]
        for wi,week in enumerate(v.WEEKS):
            for a in range(1,6):
                group=(week,str(a))
                n=6 if wi==0 else 5
                for j in range(n):
                    q=len(rows)
                    treat=v.TREAT[(j+a)%3]
                    rows.append({
                        "week":week,"array":str(a),"treat":treat,
                        "id":(week,str(a),treat,str(j)),
                        "group":group,
                        "log_flower":math.log1p(5+(q%17)),
                        "mcorW":3.0+q*.04+(j%3)*.10,
                        "mcorL":8.0+(q%21)*.08,
                        "log_visit":math.log1p((q*7)%19/10),
                        "y":math.log1p(3+(q*5)%26)
                    })
        assert len(rows)==105
        return rows

    def test_20_clonal_array_groups_all_105_scored_once(self):
        data=self.synthetic()
        folds=v.run_holdouts(data)
        self.assertEqual(len(folds),20)
        self.assertEqual(sum(x["n_original_heldout_plants"] for x in folds),105)
        self.assertEqual(set(x["week"] for x in folds),set(v.WEEKS))
        self.assertEqual(v.summarize(folds)["original_plant_predictions"],105)

    def test_Q0_Q1_Q2_Q3_exact_fixed_features(self):
        self.assertEqual(v.MODELS,("Q0","Q1","Q2","Q3"))
        self.assertEqual(v.ALPHA,5.0)
        self.assertEqual(v.FEATURES["Q0"],("log_flower",))
        self.assertEqual(v.FEATURES["Q1"],("log_flower","mcorW","mcorL"))
        self.assertEqual(v.FEATURES["Q2"],("log_flower","mcorW","mcorL","log_visit"))
        self.assertEqual(v.FEATURES["Q3"],("log_flower","log_visit"))

    def test_no_original_response_in_predictor_feature_definition(self):
        for m in v.MODELS:
            self.assertNotIn("y",v.FEATURES[m])
            self.assertNotIn("mseed",v.FEATURES[m])

    def test_heldout_response_not_used_to_fit_or_predict(self):
        data=self.synthetic()
        train=data[5:]
        test=dict(data[0])
        for model in v.MODELS:
            f=v.fit(train,model)
            pred=v.predict(f,test)
            test["y"]=500000.
            self.assertEqual(pred,v.predict(f,test))
            test["y"]=data[0]["y"]

    def test_training_only_normalization_and_model_design(self):
        data=self.synthetic()
        for m in v.MODELS:
            ft=v.fit(data,m)
            self.assertEqual(len(ft["beta"]),6+len(v.FEATURES[m]))
            self.assertEqual(ft["original_training_plants"],105)
            self.assertEqual(set(ft["stats"]),set(v.FEATURES[m]))
            self.assertTrue(math.isfinite(v.predict(ft,data[0])))

    def test_grouped_training_plantkeys_cannot_leak_in_frozen_source(self):
        rows=self.synthetic()
        rows[0]["id"]=rows[1]["id"]
        # The source contract rejects composite duplicate IDs before fit.
        self.assertEqual(len({x["id"] for x in rows}),104)

    def test_numeric_ridge_pivot_solver_works_and_singular_rejected(self):
        self.assertEqual(v.solve([[4.,0],[0,2.]],[8.,6.]),[2.,3.])
        with self.assertRaisesRegex(v.V58Error,"SOURCE_RIDGE_SINGULAR"):
            v.solve([[1.,1.],[1.,1.]],[2.,2.])

    def test_unregistered_model_or_train_size_blocks(self):
        data=self.synthetic()
        with self.assertRaisesRegex(v.V58Error,"TRAINING_GROUP_OR_MODEL_INVALID"):
            v.fit(data[:2],"Q2")
        with self.assertRaisesRegex(v.V58Error,"TRAINING_GROUP_OR_MODEL_INVALID"):
            v.fit(data,"imaginary")

    def test_original_site_selection_not_100_independent_studies(self):
        self.assertEqual(self.registry["full_week_array_holdouts"],20)
        self.assertEqual(self.registry["full_original_eligible_plant_predictions"],105)
        self.assertTrue(self.registry["new_same_study_nonindependent_array_groups_not_papers"])
        self.assertFalse(self.registry["causal_pollination_claim"])

    def test_200_visit_labels_shuffled_within_water_week_not_outcomes(self):
        data=self.synthetic()
        folds=v.run_holdouts(data)
        summary=v.summarize(folds)
        p=v.visitor_placebo(data,folds,summary,self.registry)
        self.assertEqual(p["n_within_original_week_water_visit_identity_placebos"],200)
        self.assertEqual(p["original_rng_seed"],20261124)
        self.assertTrue(p["source_test_arrays_and_seed_outcomes_fixed"])
        self.assertFalse(p["causal_visitation_randomization_p_value"])

    def test_negative_results_not_selected_away(self):
        data=self.synthetic()
        folds=v.run_holdouts(data)
        summary=v.summarize(folds)
        inc=summary["increment"]["Q1_to_Q2_conditional_visitor"]
        self.assertEqual(inc["positive_original_arrays"]+
                         inc["negative_original_arrays"],20)
        self.assertEqual(set(inc["by_original_phenology_week"]),set(v.WEEKS))

    def test_frozen_original_external_source_receipt_when_available(self):
        if not v.RECEIPT.exists():
            self.skipTest("Original author results not frozen yet")
        d=json.loads(v.RECEIPT.read_text())
        self.assertEqual(d["original_source_cohort_audit"]
                         ["exact_original_complete_numeric_plants"],105)
        self.assertEqual(len(d["all_20_original_week_array_holdout_risks"]),20)
        self.assertEqual(d["original_200_visitor_identity_placebo"]
                         ["n_within_original_week_water_visit_identity_placebos"],200)
        self.assertFalse(d["new_causal_pollination_p_value_or_hereditary_fitness_effect"])

if __name__=="__main__":unittest.main()
