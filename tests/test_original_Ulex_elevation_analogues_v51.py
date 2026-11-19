"""ANTO v51 original Ulex six sites are one continuous population, not six populations."""
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
import analyze_original_Ulex_elevation_analogues_v51 as v


class OriginalUlexSiteTransport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=v.source_registry()
        cls.parent=v.load_anchor(cls.r)

    def test_source_zenodo_and_five_original_md5(self):
        self.assertEqual(self.r["zenodo_record_id"],7761289)
        self.assertEqual(set(self.r["source_md5"]),
                         {v.FRUIT,v.TRAIT,v.KIN,v.POLL,v.READ})
        self.assertEqual(self.r["source_md5"][v.FRUIT],
                         "8ca72b7d9721f6d6279e793a1526c43e")
        self.assertEqual(self.r["source_matched_individuals"],202)
        self.assertEqual(self.r["source_original_tagged_individuals"],240)

    def test_original_six_sampling_sites_one_genetically_continuous_population(self):
        self.assertTrue(self.r["sites_within_one_continuous_population"])
        self.assertFalse(self.r["six_sites_are_separate_genetic_populations"])
        self.assertEqual(self.r["n_site_holdouts"],6)
        self.assertEqual(sum(self.r["source_sites"].values()),202)
        self.assertEqual(self.r["n_independent_original_studies_added_vs_sdea_v15"],0)

    def test_original_individuals_must_not_be_pollinator_patch_rows(self):
        self.assertEqual(self.r["source_pollinator_census_n"],569)
        self.assertEqual(self.r["source_eligible_fruit_trait_individuals"],212)
        self.assertEqual(self.r["source_relatedness_matrix_genotypes"],225)

    def test_unchanged_existing_sdea_six_site_baseline_is_authoritative(self):
        self.assertEqual(self.parent["sdea_parent_original_receipt_blob_sha1"],
                         "1fb22c1028388e720f388f2345c168270b065af4")
        self.assertEqual({x["site"] for x in self.parent["original_site_folds"]},
                         set(v.SITES))
        self.assertEqual(self.parent["source_matched_individuals"],202)

    def test_md5_changed_source_is_hard_failure(self):
        p={name:b"not a real original source" for name in self.r["source_md5"]}
        with self.assertRaisesRegex(v.V51Error,"SOURCE_BYTES_UNVERIFIED"):
            v.parse_original(p,self.r)

    def test_untrusted_zenodo_url_rejected(self):
        self.assertTrue(v.allowed_zenodo("https://zenodo.org/api/records/7761289"))
        self.assertFalse(v.allowed_zenodo("https://evil.test/record"))
        self.assertFalse(v.allowed_zenodo("https://foo@zenodo.org/api/records/7761289"))
        self.assertFalse(v.allowed_zenodo("http://zenodo.org/api/records/7761289"))

    def test_source_formula_error_is_missing_not_zero(self):
        self.assertIsNone(v.number("#VALUE!"))
        self.assertIsNone(v.number("NA"))
        self.assertIsNone(v.number("."))
        self.assertEqual(v.number("0.0",nonnegative=True),0)
        with self.assertRaisesRegex(v.V51Error,"ORIGINAL_SOURCE_BAD_NUMBER"):
            v.number("-1",nonnegative=True)
        with self.assertRaisesRegex(v.V51Error,"ORIGINAL_SOURCE_BAD_NUMBER"):
            v.number("0",positive=True)
        with self.assertRaisesRegex(v.V51Error,"SOURCE_NONNUMERIC"):
            v.number("INVALID")

    def synthetic_rows(self,n=20):
        values=[]
        for i in range(n):
            values.append({"site":"TRAIN","ind":"k"+str(i),
                "elev":200+i*3.2,
                "weight":5.0+i*.3,
                "area":9+i*.17+(i%3)*.7,
                "fruit_set":max(0,min(1,.3+i*.01))})
        return values

    def test_training_only_ridge_can_fit_both_frozen_model_families(self):
        a=self.synthetic_rows(35)
        m0=v.fit(a,False)
        m1=v.fit(a,True)
        self.assertEqual(len(m0["beta"]),2)
        self.assertEqual(len(m1["beta"]),6)
        self.assertTrue(0<=v.predict(m0,a[0])<=1)
        self.assertTrue(0<=v.predict(m1,a[0])<=1)

    def test_training_fit_does_not_use_heldout_fruit_response(self):
        train=self.synthetic_rows(36)
        test={"elev":400.0,"area":15.4,"weight":8.2,"fruit_set":.8}
        models=[v.fit(train,False),v.fit(train,True)]
        original=[v.predict(f,test) for f in models]
        test["fruit_set"]=0.0
        self.assertEqual(original,[v.predict(f,test) for f in models])

    def test_exact_simple_linear_system(self):
        a=v.solve([[3,1],[1,2]],[7,5])
        self.assertAlmostEqual(a[0],1.8)
        self.assertAlmostEqual(a[1],1.6)
        with self.assertRaisesRegex(v.V51Error,"SINGULAR"):
            v.solve([[1.,1.],[1.,1.]],[1.,2.])

    def test_test_site_cannot_overlap_training(self):
        train=self.synthetic_rows(35)
        test=[dict(train[0])]
        with self.assertRaisesRegex(v.V51Error,"SITE_OR_INDIVIDUAL_LEAKAGE"):
            v.mse_model(train,test)

    def test_fixed_models_and_ten_pair_combinatorics(self):
        self.assertEqual((v.M0,v.M1),("M0","M1"))
        self.assertEqual(v.RIDGE,5)
        self.assertEqual(self.r["n_training_pairs_per_heldout"],10)
        self.assertEqual(self.r["n_all_pairs"],60)
        self.assertEqual(self.r["site_train_regimes"],
                         ["FULL5","NEAR2","FAR2","ALL2"])

    def test_original_source_full_six_site_replay_if_staged(self):
        folder=os.environ.get("ULEX_V51_SOURCE_ROOT")
        if not folder:
            self.skipTest("Official original Zenodo source is staged only in CI")
        payload=v.download_source(self.r,source_root=folder)
        matched,audit=v.parse_original(payload,self.r)
        self.assertEqual(len(matched),202)
        self.assertEqual(audit["original_patch_census_rows_quarantined"],569)
        self.assertFalse(audit["sampling_sites_are_distinct_genetic_populations"])
        output=v.report(matched,audit,self.parent,self.r)
        self.assertEqual(output["n_exact_train_two_site_pair_fits"],60)
        self.assertLessEqual(
            output["original_SDEA_v15_full5_max_absolute_MSE_parity_error"],
            1e-7)

    def test_all_source_results_match_frozen_receipt_if_available(self):
        folder=os.environ.get("ULEX_V51_SOURCE_ROOT")
        if not folder or not v.RECEIPT.exists():
            self.skipTest("Official original source AND frozen receipt required")
        payload=v.download_source(self.r,source_root=folder)
        matched,audit=v.parse_original(payload,self.r)
        output=v.report(matched,audit,self.parent,self.r)
        self.assertEqual(output,json.loads(v.RECEIPT.read_text()))
        self.assertEqual(len(output["per_original_site_complete_scores"]),6)
        self.assertTrue(output["six_sites_within_one_continuous_original_population"])
        self.assertFalse(output["new_causal_genetic_selection_or_ancestral_flower_reexpression"])


if __name__=="__main__":
    unittest.main()
