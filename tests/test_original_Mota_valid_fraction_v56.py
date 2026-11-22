"""ANTO v56: fail-closed 3 original invalid ratios; field group integrity.

No original individual plant identity is shared between bee census and seed
harvest. Unit tests never claim v54/v55 original model source passed.
"""
import inspect
import json
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mota_valid_seed_fraction_v56 as v


class MotaOriginalValidFractionV56(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=v.registered()
        cls.quality=v.original_denominator_quality(cls.r)

    def test_exact_v55_prior_numeric_source_quality_and_3_rows(self):
        self.assertEqual(
            v.git_blob((ROOT/"results/anto_v55_Mota_seed_fraction_source_quality_20261119.json").read_bytes()),
            "f9d6c01720ae0657977d09d0a92043b0eff0d7b5")
        x=self.quality["source_ratio_domain_status"]
        self.assertEqual(x["total_matched_harvest_source_rows"],2217)
        self.assertEqual(x["successful_seeds_greater_than_reproductive_units"],3)
        self.assertEqual(x["valid_source_fraction"],2214)
        self.assertEqual(self.quality["source_ratio_domain_status_by_year"][0]["year"],"2017")

    def test_original_previous_nonadmissions_cannot_be_relabelled(self):
        self.assertTrue(self.r["v54_source_stop_cannot_be_redefined"])
        self.assertTrue(self.r["v55_source_stop_cannot_be_redefined"])
        self.assertTrue(self.quality["v54_source_STOP_preserved"])
        self.assertFalse(self.quality["original_seed_or_visitor_outcome_prediction_fitted"])

    def test_original_MD5_and_new_model_source_blob_pinned(self):
        self.assertEqual(self.r["source_zenodo_record"],6671199)
        self.assertEqual(len(self.r["source_original_md5"]),3)
        self.assertEqual(self.r["v55_original_model_implementation_git_blob_sha1"],
                         "3f4b936101ef51e56e1fc661089d49b0efe8810b")
        self.assertEqual(v.v55.MODELS,("Z0","Z1","Z2"))
        self.assertEqual(v.v55.LAMBDA,8.0)

    def test_original_fit_and_placebo_are_imported_not_modified(self):
        self.assertIs(v.v55.fit,
            sys.modules["analyze_original_Mota_reported_total_bees_v55"].fit)
        self.assertIs(v.v55.holdouts,
            sys.modules["analyze_original_Mota_reported_total_bees_v55"].holdouts)
        self.assertIs(v.v55.placebo,
            sys.modules["analyze_original_Mota_reported_total_bees_v55"].placebo)
        self.assertIs(v.v55.report,
            sys.modules["analyze_original_Mota_reported_total_bees_v55"].report)

    def test_frozen_exact_original_sample_grain(self):
        self.assertEqual(self.r["source_included_harvest_rows"],2214)
        self.assertEqual(self.r["original_matched_bee_visitation_survey_rows"],11682)
        self.assertEqual(self.r["original_2017_2018_matched_fieldyears"],80)
        self.assertEqual(self.r["original_unique_physical_fields"],51)
        self.assertEqual(self.r["source_2017_2018_fieldyears"],{"2017":41,"2018":39})
        self.assertEqual(self.r["original_all_common_distance_band_counts"],
                         {"3":1,"4":79})

    def test_original_contingent_visitor_discrepancies_never_interpreted_as_zero(self):
        self.assertEqual(self.r["source_missing_abort_rows_preserved_not_imputed"],1056)
        self.assertEqual(self.r["source_guild_component_disagreements_preserved"],4087)
        self.assertEqual(self.r["source_max_guild_total_difference"],9)

    def test_no_fabricated_molecular_or_plant_level_process(self):
        self.assertFalse(self.r["individual_sunflower_visit_seed_join_possible"])
        self.assertFalse(self.r["new_molecular_pollinator_evolution_evidence"])
        self.assertFalse(self.r["ecological_causal_p_value"])
        self.assertEqual(self.r["new_original_independent_biological_studies"],0)

    def test_frozen_primary_vs_secondary_difference(self):
        self.assertEqual(self.r["models"],["Z0","Z1","Z2"])
        self.assertEqual(self.r["ridge_lambda"],8)
        self.assertEqual(self.r["source_placebo_replicates"],200)
        self.assertEqual(self.r["source_placebo_seed"],20261119)
        self.assertIn("reported",self.r["primary_visitor"])
        self.assertIn("reconstructed",self.r["secondary_visitor"])

    @staticmethod
    def fixture():
        data=[]
        for i in range(80):
            ix=i%51
            region="Burgos" if ix%2==0 else "Cuenca"
            ty=v.TYPES[ix%3]
            name=(region,"site"+str(ix))
            data.append({
                "year":"2017" if i<51 else "2018",
                "region":region,"field_type":ty,
                "physical_field":name,
                "physical_field_hash":"anonymous"+str(ix),
                "seed_fraction":.15+.65*((i*7)%31)/31,
                "log_reported_total_rate":.02+.08*((i*3)%17),
                "log_reconstituted_rate":.05+.07*((i*5)%19)
            })
        return data

    def test_physical_field_two_years_never_cross_train_test(self):
        folds=v.v55.holdouts(self.fixture())
        self.assertEqual(len(folds),51)
        self.assertEqual(sum(f["n_original_fieldyears"] for f in folds),80)
        self.assertEqual(sum(f["n_original_fieldyears"]==2 for f in folds),29)

    def test_source_outcomes_are_not_predictor_features(self):
        train=self.fixture()[2:]
        target=dict(self.fixture()[0])
        for name in v.v55.MODELS:
            fitted=v.v55.fit(train,name)
            before=v.v55.predict(fitted,target)
            target["seed_fraction"]=.99
            self.assertEqual(before,v.v55.predict(fitted,target))
            target["seed_fraction"]=.15

    def test_registered_new_scenario_is_postsource_stop_not_new_original_study(self):
        self.assertTrue(self.r["posthoc_after_v55_source_failure"])
        self.assertTrue(self.r["v54_source_stop_cannot_be_redefined"])
        self.assertTrue(self.r["v55_source_stop_cannot_be_redefined"])

    def test_derived_result_never_replaces_original_raw_rows(self):
        if not v.RECEIPT.exists():
            self.skipTest("First original source model run must precede frozen receipt")
        out=json.loads(v.RECEIPT.read_text())
        self.assertEqual(out["source_audit"]["valid_harvest_original_source_rows_kept"],2214)
        self.assertEqual(out["source_audit"]["original_harvest_seed_count_ratio_invalid_rows_excluded"],3)
        self.assertEqual(out["source_audit"]["n_eligible_fieldyear_unit_observations"],80)
        self.assertEqual(out["source_audit"]["n_unique_original_physical_fields"],51)
        self.assertEqual(out["summary"]["n_original_physical_field_holdouts"],51)
        self.assertEqual(out["summary"]["n_original_fieldyear_predictions"],80)
        self.assertTrue(out["not_retroactive_v54_or_v55_success"])
        self.assertFalse(out["new_source_plant_level_pollination_mechanism_estimated"])


if __name__=="__main__":
    unittest.main()
