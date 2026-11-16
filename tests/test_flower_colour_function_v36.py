"""Adversarial original-study evidence checks, not new chemistry fitting."""
import copy
import json
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_flower_colour_function_v36 as v


class OriginalFunctionalEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads((ROOT/"data/external/anto_v36_independent_original_white_functional_programs_20261116.json").read_text())

    def test_original_three_programs_three_families_not_new_meta_replicates(self):
        r=v.review(self.source)
        self.assertEqual(r["original_study_programs"],3)
        self.assertEqual(r["distinct_plant_families"],3)
        self.assertFalse(r["quantitatively_poolable_study_effects"])
        self.assertTrue(r["no_source_original_raw_expression_HPLC_refit"])

    def test_white_plant_with_stem_anthocyanidins_is_a_real_cross_tissue_observation(self):
        r=v.review(self.source)
        self.assertEqual(r["cross_tissue_white_flower_stem_anthocyanidin"]["white_plants_assayed"],3)
        self.assertEqual(r["cross_tissue_white_flower_stem_anthocyanidin"]["white_flower_detected"],False)
        self.assertEqual(r["cross_tissue_white_flower_stem_anthocyanidin"]["white_stem_detected"],True)
        self.assertEqual(r["native_corolla_reexpression_experiment_count"],0)

    def test_myb_cross_segregation_original_table_counts_and_not_separate_evolution_events(self):
        r=v.review(self.source)
        a=r["ipomoea_MYB_genotype_phenotype_original_table3"]
        self.assertEqual(a["B2_genotyped"],41)
        self.assertEqual(a["F3_genotyped"],12)
        self.assertEqual(a["total_genotyped"],53)
        self.assertEqual(a["concordant"],53)
        self.assertEqual(a["discordant"],0)
        self.assertTrue(a["DFR_marker_did_not_cosegregate"])
        self.assertTrue(r["genotyped_offspring_not_independent_reversion_events"])

    def test_S2_segregation_is_not_additional_allele_genotyped_offspring(self):
        r=v.review(self.source)
        self.assertEqual(r["separate_Ipomoea_S2_throat_phenotyping"]["offspring"],732)
        self.assertEqual(r["separate_Ipomoea_S2_throat_phenotyping"]["purple"],530)
        self.assertEqual(r["separate_Ipomoea_S2_throat_phenotyping"]["white"],202)
        self.assertEqual(r["ipomoea_MYB_genotype_phenotype_original_table3"]["total_genotyped"],53)

    def test_heterologous_regulator_sufficiency_with_original_20_transgenics(self):
        r=v.review(self.source)
        t=r["functional_interventions"]["Iochroma_MYBL1_heterologous_in_Nicotiana"]
        self.assertEqual(t["source_transformed_lines_grown"],20)
        self.assertEqual(t["mean_fold_lower_reported"],14)
        self.assertEqual(t["line_fold_lower_reported_range"],[11,20])
        self.assertEqual(t["target_process"],"PIGMENT_REPRESSION")
        self.assertFalse(t["endogenous_Iochroma_anthocyanin_restoration"])

    def test_Iochroma_cross_counts_not_confused_with_transgenic_numbers(self):
        r=v.review(self.source)
        c=r["Iochroma_published_crossing_design"]
        self.assertEqual(c["backcross_blue_parent_offspring"],73)
        self.assertEqual(c["backcross_white_parent_offspring"],19)
        self.assertEqual(c["all_crossing_offspring"],92)
        self.assertEqual(r["functional_interventions"]["Iochroma_MYBL1_heterologous_in_Nicotiana"]["source_transformed_lines_grown"],20)

    def test_azalea_assayed_cultivars_are_not_31_species_or_functional_enzyme(self):
        r=v.review(self.source)
        a=r["azalea_original_grain"]
        self.assertEqual(a["HPLC_wild_plus_cultivar"],31)
        self.assertEqual(a["qPCR_wild_plus_cultivar"],24)
        self.assertTrue(a["cultivars_not_independent_species"])
        self.assertFalse(a["native_gene_function_verified"])
        self.assertFalse(r["exact_v33_three_gene_FLS_model_replicated"])

    def test_evidence_levels_count_separate_original_experiment_types(self):
        r=v.review(self.source)
        levels=r["source_evidence_levels"]
        self.assertEqual(levels["cross_tissue_white_flower_stem_biochemistry"],["ipomoea_2020"])
        self.assertEqual(levels["original_genotype_flower_cosegregation"],["ipomoea_2020"])
        self.assertEqual(levels["heterologous_gene_pigment_suppression"],["iochroma_MYBL1_2018"])
        self.assertEqual(levels["native_white_flower_restoration"],[])
        self.assertEqual(levels["natural_ancestral_reacquisition_with_persistence"],[])

    def test_unverified_native_rescue_must_fail(self):
        d=copy.deepcopy(self.source)
        d["original_programs"][1]["heterologous_regulator_perturbation"]["not_endogenous_Iochroma_rescue"]=False
        with self.assertRaisesRegex(ValueError,"FORGED_NATIVE_RESCUE"):
            v.review(d)

    def test_fake_stem_loss_or_white_corolla_anthocyanidins_fails(self):
        for field,val in (("white_stem_anthocyanidins_detected",False),
                          ("white_corolla_anthocyanidins_detected",True)):
            d=copy.deepcopy(self.source)
            d["original_programs"][0]["chemistry"][field]=val
            with self.assertRaisesRegex(ValueError,"TISSUE_BIOCHEMISTRY_DRIFT"):
                v.review(d)

    def test_fake_cosegregation_or_genotyped_total_mismatch_fails(self):
        for key in ("discordant_marker_phenotype_pairs","purple_phenotype_with_purple_MYB_marker"):
            d=copy.deepcopy(self.source)
            d["original_programs"][0]["crossing"]["published_individually_genotyped_MYB_offspring"][key]+=1
            with self.assertRaisesRegex(ValueError,"ORIGINAL_TABLE3_GENETIC_DENOMINATOR_DRIFT"):
                v.review(d)

    def test_fake_transformation_fold_or_line_counts_fails(self):
        for key in ("transformed_lines_grown","author_reported_floral_absorbance_fold_lower_mean"):
            d=copy.deepcopy(self.source)
            d["original_programs"][1]["heterologous_regulator_perturbation"][key]+=1
            with self.assertRaisesRegex(ValueError,"Nicotiana_TRANSFORMATION_OR_ASSAY_DRIFT"):
                v.review(d)

    def test_fake_functional_azalea_transgene_must_fail(self):
        d=copy.deepcopy(self.source)
        d["original_programs"][2]["heterologous_regulator_perturbation"]={"donor":"Rhododendron"}
        with self.assertRaisesRegex(ValueError,"AZALEA_NO_FUNCTIONAL_PERTURBATION"):
            v.review(d)

    def test_original_program_duplicate_or_new_doi_fails(self):
        d=copy.deepcopy(self.source)
        d["original_programs"][1]["doi"]="10.0000/fake"
        with self.assertRaisesRegex(ValueError,"SOURCE_DOI_IDENTITY_DRIFT"):
            v.review(d)
        d=copy.deepcopy(self.source)
        d["original_programs"].append(d["original_programs"][0])
        with self.assertRaisesRegex(ValueError,"ORIGINAL_PROGRAM_SET_DRIFT"):
            v.review(d)

    def test_dont_claim_original_2019_Dryad_or_normalized_cross_study_effect(self):
        r=v.review(self.source)
        self.assertIsNone(r["cross_family_pooled_molecular_effect"])
        self.assertFalse(r["regain_probability_estimated"])
        self.assertTrue(r["no_direct_cross_study_identical_gene_HPLC_prediction"])
        self.assertFalse(r["direct_native_gene_edit_rescue_confirmed"])

    def test_posted_transcript_detection_is_not_functional_enzyme_recovery(self):
        d=copy.deepcopy(self.source)
        d["original_programs"][2]["expression"]["not_original_enzyme_functional_test"]=False
        with self.assertRaisesRegex(ValueError,"PCR_IS_NOT_FUNCTIONAL_ASSAY"):
            v.review(d)


if __name__=="__main__":
    unittest.main()
