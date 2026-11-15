"""Adversarial evidence-grain tests for v35 published external comparisons."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_external_white_flower_v35 as v


class V35ExternalEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/"data/external/anto_v35_white_flower_cross_family_published_source_registry_20261115.json").read_text())

    def test_two_actual_primary_doi_registered_and_separate_families(self):
        result=v.review(self.data)
        self.assertEqual(result["n_external_independent_study_programs"],2)
        self.assertEqual(result["n_independent_original_family_groups"],2)
        self.assertEqual(result["independent_family_names"],["Convolvulaceae","Ericaceae"])
        self.assertEqual(result["author_published_dois"],[
            "10.1371/journal.pone.0231263","10.2503/hortj.UTD-142"])

    def test_white_state_does_not_uniquely_identify_expression_state(self):
        result=v.review(self.data)
        self.assertTrue(result["all_primary_programs_share_visible_white_state"])
        self.assertTrue(result["distinct_documented_white_transcription_state_classes"]==2)
        self.assertTrue(result["white_colour_does_not_uniquely_identify_regulatory_state"])
        self.assertFalse(result["original_ancestral_gene_function_demonstrated"])

    def test_published_gene_sets_not_misrepresented_as_same_three_gene_model(self):
        data=self.data
        az=next(x for x in data["external_independent_primary_study_programs"] if "azalea" in x["id"])
        self.assertEqual(az["measured_genes"],["F3'H","F3'5'H","DFR","ANS"])
        self.assertNotIn("FLS",az["measured_genes"])
        result=v.review(data)
        self.assertFalse(result["ANTOv33_exact_F3Hprime_F3H5prime_FLS_model_externally_replayed"])

    def test_distinguish_original_rhododendron_accessions_from_species(self):
        result=v.review(self.data)
        self.assertEqual(result["study_observation_grains"]["hirado_azalea_2020"]["HPLC_samples"],31)
        self.assertEqual(result["study_observation_grains"]["hirado_azalea_2020"]["qPCR_samples"],24)
        self.assertTrue(result["raw_joint_individual_HPLC_qPCR_data_unavailable_in_registry"])

    def test_distinguish_independent_original_Ipomoea_qpcr_HPLC_sample_sets(self):
        r=v.review(self.data)
        z=r["study_observation_grains"]["ipomoea_lacunosa_2020"]
        self.assertEqual(z["field_qpcr_individuals"],8)
        self.assertEqual(z["I_lacunosa_HPLC_individuals"],6)
        self.assertEqual(z["additional_I_cordatotriloba_HPLC_individuals"],2)
        self.assertFalse(r["same_individual_HPLC_RNA_matched_comparison_computed"])

    def test_external_context_studies_not_counted_as_new_replication(self):
        r=v.review(self.data)
        self.assertEqual(r["n_context_prior_studies_not_new_reanalyses"],2)
        self.assertEqual(r["n_external_independent_study_programs"],2)
        self.assertFalse(r["independent_quantitative_meta_analysis_possible"])

    def test_access_block_gate_cannot_be_silently_promoted(self):
        d=copy.deepcopy(self.data)
        d["context_not_new_independent_primary_qpcr_hplc_reanalysis"][0]["original_dryad_RAR_access"]="SOURCE_VERIFIED"
        with self.assertRaisesRegex(ValueError,"DRYAD_ACCESS_CLAIM_DRIFT"):
            v.review(d)

    def test_fake_source_doi_fails(self):
        d=copy.deepcopy(self.data)
        d["external_independent_primary_study_programs"][0]["doi"]="10.0000/fake"
        with self.assertRaisesRegex(ValueError,"SOURCE_DOI_DRIFT"):
            v.review(d)

    def test_fls_fabrication_fails(self):
        d=copy.deepcopy(self.data)
        d["external_independent_primary_study_programs"][0]["measured_genes"].append("FLS")
        with self.assertRaisesRegex(ValueError,"UNMEASURED_GENE"):
            v.review(d)

    def test_cultivar_count_not_species_count(self):
        d=copy.deepcopy(self.data)
        d["external_independent_primary_study_programs"][0]["reported_hplc_cultivar_samples"]=31
        with self.assertRaisesRegex(ValueError,"SOURCE_DENOMINATOR_DRIFT"):
            v.review(d)

    def test_forged_joint_HPLC_qpcr_same_individuals_fails(self):
        d=copy.deepcopy(self.data)
        d["external_independent_primary_study_programs"][0]["hplc_rna_same_biological_individuals_confirmed"]=True
        with self.assertRaisesRegex(ValueError,"FORGED_JOINT_MEASUREMENT"):
            v.review(d)

    def test_gene_expression_cannot_be_promoted_to_verified_enzyme_function(self):
        d=copy.deepcopy(self.data)
        d["external_independent_primary_study_programs"][0]["enzyme_function_or_ancestral_reversal_measured"]=True
        with self.assertRaisesRegex(ValueError,"UNSUPPORTED_FUNCTIONAL_ASSAY"):
            v.review(d)

    def test_azalea_anthocyanidin_trace_not_silently_called_zero(self):
        d=copy.deepcopy(self.data)
        d["external_independent_primary_study_programs"][0]["white_anthocyanidin_aglycone_exception"]="ABSENT"
        with self.assertRaisesRegex(ValueError,"ANTHOCYANIDIN_TRACE"):
            v.review(d)

    def test_shared_DFR_relative_to_colored_comparison_is_explicit_but_not_cross_calibrated(self):
        result=v.review(self.data)
        self.assertTrue(result["common_DFR_is_measured_in_both_original_programs"])
        self.assertTrue(result["source_reported_DFR_relative_contrasts_differ"])
        self.assertFalse(result["absolute_DFR_expression_is_cross_study_comparable"])
        descriptions=result["common_DFR_within_study_descriptions"]
        self.assertIn("WHITE_DFR_DETECTED_SIMILAR_TO_COLORED_PUBLISHED_DESCRIPTION",
                      list(descriptions.values()))
        self.assertIn("WHITE_DFR_OVER_100_FOLD_DOWN_VS_PURPLE_PUBLISHED_DESCRIPTION",
                      list(descriptions.values()))
        wrong=copy.deepcopy(self.data)
        wrong["external_independent_primary_study_programs"][0][
            "common_DFR_relative_to_colored_reference"]="SAME_AS_IPOMOEA"
        with self.assertRaisesRegex(ValueError,"SOURCE_COMMON_DFR_RELATIVE_CONTRAST_DRIFT"):
            v.review(wrong)

    def test_intergenerational_regain_not_claimed(self):
        result=v.review(self.data)
        self.assertFalse(result["intergenerational_natural_regain_demonstrated"])
        self.assertTrue(result["same_visible_white_class_different_expression_reporting"])
        self.assertTrue(result["gene_expression_units_and_pigment_methods_not_poolable"])

    def test_no_numeric_RMSE_or_pvalue_fabrication(self):
        result=v.review(self.data)
        self.assertNotIn("p_value",result)
        self.assertNotIn("pooled_RMSE",result)
        self.assertEqual(result["numeric_cross_study_effect_size"],None)
        self.assertEqual(result["new_source_raw_expression_rows_fitted"],0)


if __name__=="__main__":
    unittest.main()
