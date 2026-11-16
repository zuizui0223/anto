"""ANTO v38: hard-fail source-linked flower-colour channel vs fitness audit."""
import copy
import json
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_petunia_pollinator_channel_v38 as v

class OriginalPollinatorStageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads((ROOT/"data/external/anto_v38_petunia_original_pollinator_seed_set_programs_20261116.json").read_text())
        cls.parent=json.loads((ROOT/"results/anto_v37_original_natural_regain_native_rescue_20261116.json").read_text())
        cls.output=v.review(cls.source,cls.parent)

    def test_original_four_papers_only_one_shared_petunia_radiation(self):
        x=self.output
        self.assertEqual(x["n_original_publications_in_v38_registry"],4)
        self.assertEqual(x["n_botanical_families"],1)
        self.assertEqual(x["n_independent_petunia_radiations"],1)
        self.assertEqual(x["n_new_genetically_independent_colour_evolution_experiments"],0)
        self.assertTrue(x["hoballah2007_existing_v37_study_deduplicated"])

    def test_existing_v37_moth_30_6_preserved_no_gene_colour_fitness(self):
        x=self.output["author_published_count_rechecks"]
        self.assertEqual(x["hoballah2007_moth_white"],30)
        self.assertEqual(x["hoballah2007_moth_coloured"],6)
        self.assertEqual(x["hoballah2007_moth_n"],36)
        self.assertAlmostEqual(x["hoballah2007_exact_binomial_two_sided_p"],0.00006960169412195683,places=16)
        self.assertFalse(self.output["same_AN2_colour_perturbation_field_capsule_or_seedset_measured"])

    def test_2013_two_AN2_backgrounds_first_choice_and_exact_binomial_distinct_author_chi2(self):
        x=self.output["author_published_count_rechecks"]
        self.assertEqual(x["2013_axillaris_coloured_first"],15)
        self.assertEqual(x["2013_axillaris_white_first"],14)
        self.assertEqual(x["2013_axillaris_total_first"],29)
        self.assertAlmostEqual(x["2013_axillaris_exact_binomial_two_sided_p"],1)
        self.assertEqual(x["2013_inflata_white_first"],15)
        self.assertEqual(x["2013_inflata_purple_first"],3)
        self.assertEqual(x["2013_inflata_total_first"],18)
        self.assertAlmostEqual(x["2013_inflata_exact_binomial_two_sided_p"],.007537841796875)
        self.assertEqual(x["2013_inflata_author_chi_square_p"],.005)
        self.assertTrue(self.output["source_chi_square_p_does_not_need_equal_exact_binomial_p"])

    def test_bee_data_hives_and_field_days_are_not_species(self):
        counts=self.output["AN2_NIL_secondary_guild_sources"]
        self.assertEqual(counts["axillaris_bumblebee_hives"],4)
        self.assertEqual(counts["axillaris_field_bee_days"],7)
        self.assertEqual(counts["inflata_bumblebee_hives"],4)
        self.assertEqual(counts["inflata_field_bee_days"],7)
        self.assertTrue(counts["new_guild_can_be_attracted_without_detectable_old_guild_loss"])
        self.assertFalse(counts["non_significance_proves_equal_original_guild_visits"])

    def test_original_AN2_intervention_visitation_not_pollen_or_seed_set(self):
        x=self.output["experimental_contrast_observation_matrix"]
        for key in ("Hoballah2007_AN2_transformed_colour","DellOlivo2013_axillaris_AN2_gain",
                    "DellOlivo2013_inflata_AN2_loss"):
            self.assertTrue(x[key]["AN2_colour_changed"])
            self.assertTrue(x[key]["pollinator_first_choice_measured"])
            self.assertFalse(x[key]["same_genotype_field_capsule_or_seed_number_measured"])
            self.assertFalse(x[key]["same_genotype_pollen_deposition_count_measured"])

    def test_2007_wild_axillaris_day_and_night_capsules_real_but_not_edited_gene(self):
        x=self.output["wild_field_capsule_published_source"]
        self.assertEqual(x["n_source_localities_P_axillaris"],2)
        self.assertTrue(x["capsules_day_visitor_access_observed"])
        self.assertTrue(x["capsules_night_visitor_access_observed"])
        self.assertEqual(x["jose_ignacio_original_Mann_Whitney_p"],.477)
        self.assertEqual(x["carmelo_original_Mann_Whitney_p"],.167)
        self.assertEqual(x["n_source_reported_replicates_per_locality"],4)
        self.assertFalse(x["transgenic_AN2_colour_was_the_field_exclusion_treatment"])
        self.assertFalse(x["day_night_equal_effectiveness_established"])
        self.assertFalse(x["measured_seed_counts"])

    def test_wild_integrifolia_no_night_capsules_not_universal_day_dependence(self):
        x=self.output["wild_field_capsule_published_source"]
        self.assertFalse(x["integrifolia_capsules_without_day_visits"])
        self.assertTrue(x["species_specific_pollinator_channel_availability"])

    def test_no_detection_does_not_imply_equivalence(self):
        r=self.output
        self.assertTrue(r["cannot_infer_causal_same_gene_reproductive_assurance"])
        self.assertFalse(r["non_significant_old_guild_preference_is_equivalence_evidence"])

    def test_2012_nectar_low_same_choice_different_seed_outcome(self):
        r=self.output["nectar_service_efficiency_source"]
        self.assertEqual(r["relative_nectar"],"APPROXIMATELY_ONE_THIRD_VOLUME")
        self.assertTrue(r["initial_attraction_difference_not_detected"])
        self.assertTrue(r["probing_duration_reduced"])
        self.assertTrue(r["low_nectar_seed_production_reduced"])
        self.assertTrue(r["nectar_supplementation_recovered_seed_set"])
        self.assertFalse(r["this_is_AN2_colour_perturbation"])

    def test_2012_seed_amount_not_fabricated_when_abstract_only(self):
        r=self.output["nectar_service_efficiency_source"]
        self.assertIsNone(r["seed_set_numeric_risk_ratio"])
        self.assertIsNone(r["seed_set_model_p_value"])
        self.assertEqual(r["original_seed_count_rows_recomputed"],0)

    def test_2018_Psecreta_contacts_not_actual_stigma_deposition_measure(self):
        r=self.output["wild_Psecreta_pollen_reproductive_contact"]
        self.assertEqual(r["total_observed_pollinator_visits"],38)
        self.assertEqual(r["observed_reproductive_contact_visits"],22)
        self.assertAlmostEqual(r["descriptive_fraction_of_observed_visits_with_contact"],22/38)
        self.assertFalse(r["successful_pollen_deposition_per_visit_confirmed"])
        self.assertFalse(r["corresponding_visited_flowers_seedset_measured"])
        self.assertFalse(r["ancestral_AN2_reactivation_experiment_on_field_plants"])

    def test_stage_overlap_paper_count_does_not_identify_full_AN2_colour_chain(self):
        r=self.output
        self.assertEqual(r["n_original_programs_with_real_capsule_or_seed_outcome"],2)
        self.assertEqual(r["n_original_programs_with_colour_AN2_first_choice_manipulation"],2)
        self.assertEqual(r["n_same_colour_AN2_intervention_with_direct_reproductive_output"],0)
        self.assertEqual(r["joint_intergenerational_genetic_colour_to_fitness_chains"],0)
        self.assertEqual(r["fully_linked_AN2_molecular_to_long_term_fitness_verdict"],
                         "NOT_IDENTIFIED_SAME_TREATMENT_END_TO_END")

    def test_fake_2007_edited_colour_field_capsule_claim_fails(self):
        d=copy.deepcopy(self.source)
        d["programs"][0]["experimental_substudies"]["independent_wild_native_pollinator_exclusion"][
            "organism_genotype_equivalent_to_transgenic_color_choice"]=True
        with self.assertRaisesRegex(ValueError,"FORGED_COLOUR_TO_WILD_CAPSULE_LINK"):
            v.review(d,self.parent)

    def test_fake_2013_same_NIL_seedset_data_fails(self):
        d=copy.deepcopy(self.source)
        d["programs"][1]["experimental_substudies"]["axillaris_colour_gain"]["same_NIL_seedset_measured"]=True
        with self.assertRaisesRegex(ValueError,"NIL_POLLEN_FITNESS_OUTCOME_NOT_MEASURED"):
            v.review(d,self.parent)

    def test_fake_2012_colour_change_causal_seedset_fails(self):
        d=copy.deepcopy(self.source)
        d["programs"][2]["experimental_substudies"]["low_nectar_introgression"]["novel_flower_colour_or_AN2_perturbation"]=True
        with self.assertRaisesRegex(ValueError,"NECTAR_NOT_AN2_COLOUR"):
            v.review(d,self.parent)

    def test_fake_2018_pollen_deposition_or_seedcounts_fails(self):
        d=copy.deepcopy(self.source)
        d["programs"][3]["experimental_substudies"]["natural_pink_Psecreta"]["viable_pollen_deposition_determined_per_visit"]=True
        with self.assertRaisesRegex(ValueError,"PSECRETA_CONTACT_NOT_SEED_OR_POLLEN_FITNESS"):
            v.review(d,self.parent)

    def test_genetic_ancestor_loss_not_corollary_of_2018_pollen_field_study(self):
        r=self.output
        self.assertFalse(r["natural_ancestral_allele_recovery_caused_contact_verified_in_same_individuals"])
        self.assertFalse(r["new_externally_replicated_v33_species_RNA_HPLC_model"])

    def test_missing_or_duplicate_program_prevents_false_replicates(self):
        d=copy.deepcopy(self.source)
        d["programs"].append(copy.deepcopy(d["programs"][0]))
        with self.assertRaisesRegex(ValueError,"PAPER_STUDY_SET_DRIFT"):
            v.review(d,self.parent)
        d=copy.deepcopy(self.source)
        d["programs"][0]["doi"]="10.0000/fake"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_STUDY_DOI_DRIFT"):
            v.review(d,self.parent)

    def test_parent_v37_source_record_must_not_silently_change(self):
        d=copy.deepcopy(self.parent)
        d["original_pollinator_preference"]["Manduca_first_choose_white"]=29
        with self.assertRaisesRegex(ValueError,"V37_SOURCE_COUNT_DRIFT"):
            v.review(self.source,d)

if __name__=="__main__":
    unittest.main()
