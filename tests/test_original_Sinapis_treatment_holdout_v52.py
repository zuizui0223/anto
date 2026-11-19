"""Adverse source, numerical, biological-grain and no-causal-claims checks."""
from __future__ import annotations

import copy
import json
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Sinapis_treatment_holdout_v52 as v


class OriginalSinapisTreatmentTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=v.read_registry()
        cls.prior=v.prior_source()

    def test_original_MD5_and_author_plant_grain(self):
        self.assertEqual(self.r["zenodo_record_id"],7335086)
        self.assertEqual(self.r["author_source_MD5"]["Int_data.csv"],
                         "d67c32cc5aa357fe6e6d952a801ea38f")
        self.assertEqual(self.r["original_total_plant_records"],53)
        self.assertEqual(self.r["original_complete_case_plants"],52)
        self.assertEqual(self.r["original_visitor_group_rows"],427)
        self.assertEqual(self.r["original_rain_shelters"],5)

    def test_frozen_parent_is_prior_not_new_independent_study(self):
        self.assertEqual(self.prior["original_file_MD5"],
                         self.r["author_source_MD5"])
        self.assertEqual(self.prior["parent_git_blob_sha1"],
                         "fa2419d32b74f22f67ea9cbfd79123109db4fdf0")
        self.assertEqual(len(self.prior["original_shelter_heldout_folds"]),5)
        self.assertTrue(self.prior["this_is_prior_replay_not_new_evidence"])
        self.assertEqual(self.r["new_original_biological_studies_vs_sdea_v17"],0)

    def test_treatment_dummies_must_be_excluded_in_unseen_treatment(self):
        self.assertTrue(self.r["model_treatment_indicators_excluded"])
        self.assertEqual(v.MODELS,("U0","U1","U2","U3"))
        self.assertEqual(v.FEATURES["U1"],
                         ("log_flower","display","petal_width"))
        self.assertEqual(v.FEATURES["U2"],
                         ("log_flower","display","petal_width","log_visit_rate"))
        self.assertTrue(all("treatment" not in f
                            for fields in v.FEATURES.values() for f in fields))

    @staticmethod
    def example_plant(treat,shelter,i):
        return {"treatment":treat,"shelter":str(shelter),
                "log_flower":1.0+.18*i+.12*shelter,
                "display":5+.3*i+.4*shelter,
                "petal_width":.9+.07*(i%4)+.01*shelter,
                "log_visit_rate":.1+.16*(i%7)+.04*shelter,
                "log_seed":2+.2*i+.11*shelter+.28*(i%7)}

    def synthetic_52_rows(self):
        rows=[]
        for si,s in enumerate(v.SHELTERS):
            for ti,t in enumerate(v.TREATMENTS):
                n=4 if si*3+ti<7 else 3
                for j in range(n):
                    k=len(rows)
                    rows.append(self.example_plant(t,int(s),k))
        self.assertEqual(len(rows),52)
        return rows

    def test_numerical_ridge_with_unpenalized_intercept(self):
        sol=v.solve([[2.,0.],[0.,3.]],[10.,6.])
        self.assertAlmostEqual(sol[0],5)
        self.assertAlmostEqual(sol[1],2)
        with self.assertRaisesRegex(v.V52Error,"SINGULAR"):
            v.solve([[1.,1.],[1.,1.]],[1.,2.])

    def test_all_four_new_models_use_same_original_y(self):
        rows=self.synthetic_52_rows()
        for m in v.MODELS:
            tr=[r for r in rows if r["treatment"]!="Control"]
            fitted=v.fit(tr,m)
            self.assertEqual(fitted["model"],m)
            self.assertEqual(fitted["n_train"],len(tr))
            self.assertNotIn("treatment",fitted["stats"])
            self.assertEqual(len(fitted["beta"]),len(v.FEATURES[m])+1)

    def test_heldout_log_seed_outcome_never_used_in_prediction(self):
        rows=self.synthetic_52_rows()
        train=[r for r in rows if r["treatment"]!="Drought"]
        observed=dict(next(r for r in rows if r["treatment"]=="Drought"))
        for m in v.MODELS:
            fit=v.fit(train,m)
            start=v.predict(fit,observed)
            observed["log_seed"]+=100000
            self.assertEqual(start,v.predict(fit,observed))
            observed["log_seed"]-=100000

    def test_three_whole_treatment_folds_cover_each_plant_exactly_once(self):
        rows=self.synthetic_52_rows()
        got=v.treatment_folds(rows)
        self.assertEqual([x["heldout_treatment"] for x in got],
                         list(v.TREATMENTS))
        self.assertEqual(sum(x["n_heldout_original_plants"] for x in got),52)
        self.assertEqual(set(got[0]["mse"]),set(v.MODELS))
        for x in got:
            self.assertNotIn(x["heldout_treatment"],
                             x["original_train_treatments"])
            self.assertEqual(x["n_heldout_original_plants"]+
                             x["n_training_original_plants"],52)

    def test_unknown_model_or_insufficient_training_fails_closed(self):
        rows=self.synthetic_52_rows()
        with self.assertRaisesRegex(v.V52Error,"UNREGISTERED_MODEL"):
            v.fit(rows,"fake")
        with self.assertRaisesRegex(v.V52Error,"SMALL_TRAINING"):
            v.fit(rows[:3],"U2")

    def test_prediction_clarity_uses_seed_output_not_reproduction_causality(self):
        self.assertTrue(self.r["no_molecular_genotype"])
        self.assertTrue(self.r["no_pollen_deposition_measured"])
        self.assertTrue(self.r["no_ancestral_flower_regain_measured"])
        self.assertFalse(self.r["randomization_causal_p_value"])
        self.assertEqual(self.r["primary_comparison"],
                         "U1_MSE minus U2_MSE")

    def test_quantile_interpolation(self):
        self.assertEqual(v.quantile([0,1,2,3,4],.5),2)
        self.assertAlmostEqual(v.quantile([0.,1.],.975),.975)

    def test_source_identity_bytes_mutation_fails_before_import(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"sdea_mutated.py"
            p.write_bytes(b"# counterfeited code\n")
            with self.assertRaisesRegex(v.V52Error,"SOURCE_IMPLEMENTATION_CHANGED"):
                v.load_frozen_sdea(p,self.r)

    def test_placebo_does_not_invent_causal_randomization_p(self):
        rows=self.synthetic_52_rows()
        observed=v.treatment_folds(rows)
        placebo=v.visit_identity_placebo(rows,observed,self.r)
        self.assertEqual(placebo["permutation_replicates"],200)
        self.assertEqual(placebo["seed"],20261119)
        self.assertFalse(placebo["valid_randomized_causal_p_value"])
        self.assertTrue(0<=placebo["n_placebo_gains_at_or_above_observed"]<=200)

    def test_strict_source_integrity_with_immutable_SDEA_module_when_staged(self):
        module_file=os.environ.get("SDEA_V52_PINNED_SOURCE_CODE")
        source_root=os.environ.get("SINAPIS_V52_SOURCE_ROOT")
        if not module_file or not source_root:
            self.skipTest("Pinned user-owned SDEA v17 code and original Zenodo files required")
        mod=v.load_frozen_sdea(module_file,self.r)
        raw=v.source_payload(mod,self.r,source_root=source_root)
        rows,audit=v.original_cohort(mod,raw,self.r,self.prior)
        self.assertEqual(len(rows),52)
        self.assertEqual(audit["original_visitor_group_rows"],427)
        self.assertEqual(audit["n_original_rain_shelters"],5)
        self.assertLessEqual(audit["original_sdea_v17_5_shelter_20_mse_max_absolute_error"],1e-8)
        self.assertEqual(sum(audit["treatment_eligible_plant_counts"].values()),52)

    def test_full_immutable_52_plant_three_treatment_replay_when_available(self):
        module_file=os.environ.get("SDEA_V52_PINNED_SOURCE_CODE")
        source_root=os.environ.get("SINAPIS_V52_SOURCE_ROOT")
        if not module_file or not source_root or not v.RECEIPT.exists():
            self.skipTest("Pinned source, data and frozen 200-placebo receipt required")
        mod=v.load_frozen_sdea(module_file,self.r)
        rows,audit=v.original_cohort(
            mod,v.source_payload(mod,self.r,source_root),self.r,self.prior)
        output=v.report(rows,audit,self.prior,self.r)
        self.assertEqual(output,json.loads(v.RECEIPT.read_text()))
        self.assertEqual(len(output["original_three_complete_treatment_folds"]),3)
        self.assertFalse(output["new_causal_p_values_or_confidence_intervals"])


if __name__=="__main__":
    unittest.main()
