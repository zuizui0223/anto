"""A source title/DOI is not a row, fit, or ecological independent replication."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import discover_original_Physalis_Mertensia_v57 as m

class SourceOnlyThirdStudy(unittest.TestCase):
    def test_two_original_studies_and_not_four_independent_cohorts(self):
        r=m.read_registry()
        self.assertEqual(r["physalis"]["original_dataset_doi"],"10.5061/dryad.xpnvx0kvp")
        self.assertEqual(r["mertensia"]["original_dataset_doi"],"10.7280/D16D7Z")
        self.assertTrue(r["physalis_2020_and_2021_not_new_independent_papers"])

    def test_official_origins_not_fallback_arbitrary_hosts(self):
        self.assertTrue(m.allowed("https://datadryad.org/api/v2/datasets"))
        self.assertTrue(m.allowed("https://zenodo.org/api/records"))
        self.assertFalse(m.allowed("https://example.org/2021.csv"))
        self.assertFalse(m.allowed("http://datadryad.org/data"))
        self.assertFalse(m.allowed("https://user@zenodo.org/example"))

    def test_original_doi_is_pct_encoded_for_official_api(self):
        import urllib.parse
        self.assertEqual(urllib.parse.quote("doi:10.5061/dryad.xpnvx0kvp",safe=""),
                        "doi%3A10.5061%2Fdryad.xpnvx0kvp")

    def test_original_causal_nonclaims(self):
        r=m.read_registry()
        self.assertTrue(r["no_ecological_effect_fit"])
        self.assertTrue(r["no_original_raw_data_rehosted"])
        self.assertEqual(r["new_independent_original_study_results_admitted"],0)
        self.assertFalse(r["valid_causal_p_value_computed"])

    def test_derived_pollinator_importance_not_as_predictor(self):
        r=m.read_registry()
        bad=set(r["physalis"]["forbidden_outcome_leaking_predictors"])
        self.assertTrue({"poll_imp","prop_fs","avg_fruit_weight"}.issubset(bad))

    def test_source_reproduction_grain_not_imputed_from_article(self):
        r=m.read_registry()
        self.assertEqual(r["mertensia"]["original_paper_seed_n_not_assumed_source_cohort"],114)
        self.assertEqual(r["minimum_original_sameplant_keys"],30)
        self.assertEqual(r["minimum_original_blocks"],3)

    def test_prior_v54_v55_source_stops_retained(self):
        self.assertTrue(m.read_registry()["prior_anto_v54_v55_source_STOPS_preserved"])

if __name__=="__main__":unittest.main()
