"""Adverse original-source-grain tests for Phlox 2014 published spatial study."""
import copy
import json
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_Phlox2014_spatial_two_locus_source_v41 as m


class Phlox2014SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads((ROOT/"data/external/anto_v41_Phlox2014_original_spatial_study_tableS1_provisional_20261117.json").read_text())
        cls.out=m.review(cls.source)

    def test_study_identity_and_no_historical_time_series_forgery(self):
        a=self.out
        self.assertEqual(a["original_publication_doi"],"10.1016/j.cub.2014.07.027")
        self.assertFalse(a["new_field_evolutionary_time_series_fitted"])
        self.assertFalse(a["independent_full_2014_MCMC_reestimated"])
        self.assertTrue(a["cross_sectional_population_and_one_offspring_generation_source"])

    def test_author_study_12024_32pop_and_2581_382_mothers(self):
        x=self.out["original_main_article_observation_denominators"]
        self.assertEqual(x["adult_phenotypes_main_paper_reported"],12024)
        self.assertEqual(x["field_populations"],32)
        self.assertEqual(x["transects"],5)
        self.assertEqual(x["offspring_phenotypes_author_reported"],2581)
        self.assertEqual(x["field_maternal_plants_reported"],382)

    def test_original_source_phenotypic_fitness_two_zone_values(self):
        f=self.out["author_published_fitness_per_region"]
        self.assertEqual(f["allopatry"],{"LB":1.,"DB":.86,"LR":.89,"DR":.90})
        self.assertEqual(f["sympatry"],{"LB":.17,"DB":.53,"LR":.86,"DR":1.})
        self.assertAlmostEqual(f["sympatry"]["DR"]/f["sympatry"]["LB"],1/.17)

    def test_original_epistasis_exact_from_rounded_fitness(self):
        x=self.out["independent_arithmetic_on_original_published_epistasis"]
        self.assertAlmostEqual(x["allopatry"],.15,places=12)
        self.assertAlmostEqual(x["sympatry"],-.22,places=12)
        self.assertAlmostEqual(x["sympatry_minus_allopatry"],-.37,places=12)
        self.assertTrue(self.out["fitness_epistasis_changes_sign_across_regions"])

    def test_source_locus_selection_and_original_CI_not_new_p_values(self):
        x=self.out["author_original_table1_locus_selection"]
        self.assertEqual(x["allopatry"]["dark_intensity_allele_I"]["estimate"],-.14)
        self.assertEqual(x["sympatry"]["ancestral_light_intensity_allele_i"]["estimate"],-.47)
        self.assertEqual(x["sympatry"]["ancestral_light_intensity_allele_i"]["author_interval_95"],[-.69,-.3])
        self.assertTrue(self.out["source_2014_intensity_coefficient_not_2012_viability_proxy_s"])
        self.assertIsNone(self.out["new_original_selection_coefficient_estimate"])

    def test_hue_locus_2012_field_null_cannot_override_2014_spatial_result(self):
        self.assertTrue(self.out["different_2012_and_2014_hue_estimands"])
        self.assertTrue(self.out["2012_hue_field_hybridization_effect_not_detected"])
        self.assertTrue(self.out["2014_spatial_hue_allele_selection_author_inferred"])

    def test_32_rows_and_five_original_transect_counts(self):
        a=self.out["provisional_original_S1_text_transcription"]
        self.assertEqual(a["n_populations_digitized"],32)
        self.assertEqual(a["transect_populations"],{"1":9,"2":7,"3":6,"4":6,"5":4})
        self.assertEqual(a["first_population_code"],301)
        self.assertEqual(a["last_population_code"],320)

    def test_sum_of_phenotype_count_entries_and_recombinants(self):
        a=self.out["provisional_original_S1_text_transcription"]
        self.assertEqual(a["total_digitized_flower_count_entries"],11379)
        self.assertEqual(a["by_colour"],{"LB":4292,"DB":1396,"LR":1084,"DR":4607})
        self.assertEqual(a["recombinant_morph_count"],2480)
        self.assertAlmostEqual(a["recombinant_fraction"],2480/11379)
        self.assertEqual(a["sites_with_all_four_colours"],17)
        self.assertEqual(a["sites_with_both_LB_and_DR"],18)

    def test_published_main_vs_manual_source_table_difference_must_be_retained(self):
        x=self.out["source_denominator_reconciliation"]
        self.assertEqual(x["main_paper_claim_adults"],12024)
        self.assertEqual(x["S1_digitized_count"],11379)
        self.assertEqual(x["main_minus_S1_digitized"],645)
        self.assertEqual(x["state"],"SOURCE_ADULT_DENOMINATOR_NOT_RECONCILED")
        self.assertFalse(x["full_original_2014_survey_reproduced"])
        self.assertFalse(x["original_S1_PDF_page_image_independently_checked"])
        self.assertTrue(x["not_claimed_original_publication_error"])

    def test_different_base_fitness_normalization_retained(self):
        f=self.out["author_published_fitness_per_region"]
        self.assertEqual(f["allopatry"]["LB"],1)
        self.assertEqual(f["sympatry"]["DR"],1)
        self.assertAlmostEqual(self.out["dark_red_relative_to_light_blue_by_environment"]["allopatry"],.9)
        self.assertAlmostEqual(self.out["dark_red_relative_to_light_blue_by_environment"]["sympatry"],1/.17)
        self.assertFalse(self.out["regionally_normalized_fitness_as_common_absolute_seed_numbers"])

    def test_invalid_negative_count_fails_closed(self):
        s=copy.deepcopy(self.source)
        s["source_transcription"]["population_phenotype_counts"][0][2]=-1
        with self.assertRaisesRegex(ValueError,"INVALID_DIGITIZED_SOURCE_COUNT"):
            m.review(s)

    def test_omitted_population_and_duplicate_tips_fail(self):
        s=copy.deepcopy(self.source)
        s["source_transcription"]["population_phenotype_counts"].pop()
        with self.assertRaisesRegex(ValueError,"ORIGINAL_32_SOURCE_POPULATIONS_DRIFT"):
            m.review(s)
        s=copy.deepcopy(self.source)
        s["source_transcription"]["population_phenotype_counts"][1][1]=301
        with self.assertRaisesRegex(ValueError,"ORIGINAL_32_SOURCE_POPULATIONS_DRIFT"):
            m.review(s)

    def test_tampered_transcribed_counts_must_not_silently_fill_to_main_sample(self):
        s=copy.deepcopy(self.source)
        s["source_transcription"]["population_phenotype_counts"][0][2]+=645
        with self.assertRaisesRegex(ValueError,"SOURCE_S1_DIGITIZATION_SUM_DRIFT"):
            m.review(s)

    def test_source_new_native_ancestral_regain_experiment_fake_fails(self):
        s=copy.deepcopy(self.source)
        s["source_original"]["natural_multigeneration_allele_frequency_time_series_measured"]=True
        with self.assertRaisesRegex(ValueError,"MULTIGENERATION_TRACKING_NOT_IN_SOURCE"):
            m.review(s)

    def test_mismatched_fitness_or_epistasis_sign_fails(self):
        s=copy.deepcopy(self.source)
        s["original_study_published_model"]["original_relative_fitness_per_environment_normalized_separately"]["sympatry"]["LB"]=.4
        with self.assertRaisesRegex(ValueError,"AUTHOR_2014_FITNESS_TABLE_DRIFT"):
            m.review(s)
        s=copy.deepcopy(self.source)
        s["original_study_published_model"]["epistasis_original_table1"]["sympatry"]["estimate"]=+.22
        with self.assertRaisesRegex(ValueError,"AUTHOR_2014_EPISTASIS_DRIFT"):
            m.review(s)

    def test_source_offspring_2581_not_read_as_individually_parsed_maternal_rows(self):
        self.assertFalse(self.out["independent_offspring_tableS2_reparsed"])
        self.assertFalse(self.out["new_parent_offspring_likelihood_fit"])

    def test_cline_width_and_dispersal_are_author_not_new_weak_selection_formula(self):
        x=self.out["source_spatial_model_metadata"]
        self.assertEqual(x["light_cline_km"],4)
        self.assertEqual(x["red_cline_km"],3.7)
        self.assertEqual(x["pollen_dispersal_variance_km2"],.83)
        self.assertEqual(x["demes_per_transect"],36)
        self.assertFalse(x["weak_selection_one_locus_reconstructed"])

    def test_bad_original_DOI_fails(self):
        s=copy.deepcopy(self.source)
        s["source_original"]["doi"]="10.0000/WRONG"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_SOURCE_DOI_DRIFT"):
            m.review(s)

if __name__=="__main__":
    unittest.main()
