"""Adverse/original-author publication-grain tests for ANTO v24.

Six Table 1 species are original-published real source information; survey
percentages are ESTIMATES from herbaria/natural field, not random binomial
counts. Purposefully recruited experimental flower-morph samples are NOT
representative prevalence, molecular functional assays or new regain events.
"""
import json
import sys
import unittest
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_larter2019_published_morph_sampling_v24 import (
    ORIGINAL_TABLE,SOURCE2016,SOURCE2019_GATE,FROZEN_ROWS,
    source_rows,group_summary,rounding_delta,audit,read,serialize
)


class Larter2019PublishedMorphAscertainmentV24(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table=read(ORIGINAL_TABLE)
        cls.gene=read(SOURCE2016)
        cls.gate=read(SOURCE2019_GATE)
        cls.original=source_rows(cls.table)
        cls.result=audit(cls.table,cls.gene,cls.gate)

    def test_six_original_species_source_exact_morphs_and_survey_values(self):
        self.assertEqual(len(self.original),6)
        self.assertEqual([
          (x["published_taxon"],x["color_labels"],x["pigmented_sampled"],
           x["acyanic_sampled"],x["estimated_acyanic_frequency_percent"],
           x["surveyed_individuals"])
          for x,_,_,_ in self.original
        ],list(FROZEN_ROWS))
        self.assertEqual(sum(n for _,_,_,n in self.original),30)
        self.assertEqual(sum(x["acyanic_sampled"] for x,_,_,_ in self.original),12)
        self.assertEqual(sum(x["surveyed_individuals"] for x,_,_,_ in self.original),1025)
        self.assertFalse(self.table["raw_2019_individual_HPLC_qPCR_retrieved"])

    def test_rare_four_species_recruitment_fraction_differs_from_published_estimate(self):
        rare=self.result["rare_below_2percent_published_frequency_4_species"]
        self.assertEqual(rare["intentionally_recruited_total"],19)
        self.assertEqual(rare["intentionally_recruited_acyanic"],7)
        self.assertEqual(rare["intentionally_recruited_pigmented"],12)
        self.assertEqual(rare["original_mixed_field_herbarium_surveyed_total"],675)
        self.assertEqual(rare["experimental_acyanic_sample_fraction"],F(7,19))
        self.assertEqual(rare["survey_n_weighted_original_PUBLISHED_estimated_frequency"],F(1,75))
        self.assertEqual(rare["relative_acyanic_enrichment_of_experimental_composition"],F(525,19))
        self.assertGreater(rare["relative_acyanic_enrichment_rounding_interval_not_sampling_CI"][0],26)
        self.assertLess(rare["relative_acyanic_enrichment_rounding_interval_not_sampling_CI"][1],29)
        self.assertTrue(rare["not_a_census_or_sampling_based_population_prevalence"])

    def test_full_six_morphs_are_fourty_percent_of_targeted_samples_not_population(self):
        all_six=self.result["all_six_published_polymorphic_species"]
        self.assertEqual(all_six["intentionally_recruited_total"],30)
        self.assertEqual(all_six["intentionally_recruited_acyanic"],12)
        self.assertEqual(all_six["experimental_acyanic_sample_fraction"],F(2,5))
        self.assertEqual(all_six["survey_n_weighted_original_PUBLISHED_estimated_frequency"],F(69,410))
        self.assertEqual(all_six["relative_acyanic_enrichment_of_experimental_composition"],F(164,69))
        self.assertEqual(self.result["common_at_least_30percent_published_frequency_2_species"]
                         ["intentionally_recruited_acyanic"],5)
        self.assertTrue(self.result["original_survey_mixed_herbarium_and_field_not_unbiased_iid"])

    def test_rare_vs_common_source_frequency_ranges_are_strongly_separated(self):
        source=self.result["published_species"]
        p=[x["published_morph_frequency_as_fraction"] for x in source]
        self.assertEqual(p[:4],[F(9,500),F(9,1000),F(13,1000),F(13,1000)])
        self.assertEqual(p[-2:],[F(57,100),F(33,100)])
        for r in source[:4]:
            self.assertGreater(r["rounding_relative_enrichment_min"],18)
            self.assertGreater(r["experimental_fraction_to_survey_frequency_ratio"],19)
        for r in source[-2:]:
            self.assertGreater(r["experimental_fraction_to_survey_frequency_ratio"],F(4,5))
            self.assertLess(r["experimental_fraction_to_survey_frequency_ratio"],F(5,4))

    def test_precisely_two_white_morph_singletons_not_six_independent_transitions(self):
        d=[r for r in self.result["published_species"] if
           r["original_morph_recruited_acyanic"]==1]
        self.assertEqual([r["published_2019_taxon"] for r in d],
                         ["Iochroma calycinum","Iochroma parvifolium"])
        self.assertTrue(self.result["single_acyanic_morph_recruited_in_two_original_species"])
        self.assertTrue(self.result["no_new_estimated_genetic_ancestral_regain_or_functional_retention_rate"])

    def test_original_Iochroma_parvifolium_is_yellow_acyanic_not_literally_white(self):
        yellow=[r for r in self.result["published_species"] if
                r["acyanic_morph_is_yellow_not_white"]]
        self.assertEqual([r["published_2019_taxon"] for r in yellow],
                         ["Iochroma parvifolium"])
        self.assertTrue(self.result["purple_yellow_acyanic_counted_but_not_literally_white"])

    def test_source_specific_overlap_2016_is_five_species_not_six(self):
        rows=self.result["published_species"]
        self.assertEqual(self.result["2016_original_22taxon_new_KT_panel_overlap_species"],5)
        self.assertEqual(self.result["2016_unmatched_source_taxon"],"Saracha punctata")
        self.assertFalse(rows[-1]["overlap_original_2016_GenBank_taxon_name"])
        self.assertIsNone(rows[-1]["2016_original_DFR_new_KT_accession_if_same_species"])
        self.assertTrue(all(row["2016_registrations_not_same_plant_or_2019_flower_function"]
                            for row in rows))
        self.assertTrue(self.result[
            "not_independent_2016_genomic_plus_2019_flower_assay_same_individual"])

    def test_published_rounding_sensitivity_is_not_sampling_CI(self):
        self.assertEqual(rounding_delta("1.8"),F(1,2000))
        self.assertEqual(rounding_delta("0.9"),F(1,2000))
        self.assertEqual(rounding_delta("57"),F(1,200))
        all_six=self.result["all_six_published_polymorphic_species"]
        low,high=all_six["relative_enrichment_rounding_interval_not_sampling_CI"]
        exact=all_six["relative_acyanic_enrichment_of_experimental_composition"]
        self.assertLess(low,exact)
        self.assertLess(exact,high)
        self.assertTrue(self.result["source_percentages_are_rounded_estimates_not_integer_raw_white_counts"])

    def test_random_sample_detection_probabilities_are_counterfactual_only(self):
        for r in self.result["published_species"][:4]:
            p=r["published_morph_frequency_as_fraction"]
            n=r["experimental_morph_total"]
            expected=1-(1-p)**n
            self.assertEqual(r[
              "prob_random_at_least_one_acyanic_under_ARTIFICIAL_IID_sampling"],expected)
            self.assertGreater(expected,F(3,100))
            self.assertLess(expected,F(11,100))
        self.assertTrue(self.result[
            "morph_presence_in_four_to_six_iid_random_samples_is_counterfactual_only"])

    def test_source_survey_percent_does_not_become_integer_acyanic_count(self):
        for r in self.result["published_species"]:
            self.assertNotIn("survey_acyanic_raw_count",r)
            self.assertNotIn("census_population_prevalence",r)
            self.assertNotIn("true_gene_function_positive",r)
            self.assertNotIn("independent_pigment_loss_origin_count",r)
        self.assertTrue(self.table["no_original_survey_acyanic_counts_inferred_from_rounded_percentages"])
        self.assertTrue(self.table["no_binomial_convenience_interval_from_mixed_survey"])

    def test_bad_source_phenotype_or_sample_counts_rejected(self):
        for index,field,bad in ((0,"estimated_acyanic_frequency_percent","2.8"),
                                (2,"surveyed_individuals",151),
                                (3,"color_labels","purple/white"),
                                (5,"published_taxon","Saracha quitensis"),
                                (1,"acyanic_sampled",3),
                                (0,"pigmented_sampled",4)):
            altered=deepcopy(self.table)
            altered["rows"][index][field]=bad
            with self.subTest(index=index,field=field),self.assertRaisesRegex(
                    ValueError,"Original 2019 Table 1"):
                source_rows(altered)

    def test_original_2019_raw_archive_cannot_be_claimed_without_access(self):
        changed=deepcopy(self.gate)
        changed["original_RAR_bytes_verified"]=True
        with self.assertRaisesRegex(ValueError,"raw reanalysis"):
            audit(self.table,self.gene,changed)
        changed=deepcopy(self.gate)
        changed["raw_gene_expression_sample_rows_parsed"]=10
        with self.assertRaisesRegex(ValueError,"raw reanalysis"):
            audit(self.table,self.gene,changed)

    def test_hypothetical_later_source_gene_function_cannot_be_promoted(self):
        different=deepcopy(self.gene)
        row=next(r for r in different["taxon_gene_records"]
                 if r["current_taxon_name"]=="Iochroma cyaneum")
        row["functional_joint_toolkit"]="MEASURED_FUNCTION"
        with self.assertRaisesRegex(ValueError,"2016 recorded color or competence"):
            audit(self.table,different,self.gate)

    def test_frozen_fraction_receipt_and_external_Cirsium_boundary(self):
        x=serialize(self.result)
        self.assertEqual(x["rare_below_2percent_published_frequency_4_species"]
                         ["relative_acyanic_enrichment_of_experimental_composition"],"525/19")
        self.assertEqual(x["all_six_published_polymorphic_species"]
                         ["survey_n_weighted_original_PUBLISHED_estimated_frequency"],"69/410")
        self.assertFalse(x["individually_measured_2019_gene_expression_or_HPLC_obtained"])
        self.assertTrue(x["rar_original_source_gate_currently_blocked"])
        self.assertTrue(x["no_actual_2019_expression_regression_recomputed"])
        self.assertEqual(x["Cirsium_empirical_owner"],"zuizui0223/aza3")
        json.dumps(x,sort_keys=True)
        expected=ROOT/"results/larter2019_morph_ascertainment_v24_20261028.json"
        if expected.exists():
            self.assertEqual(x,read(expected))


if __name__=="__main__":
    unittest.main()
