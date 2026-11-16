"""Fail-closed audit of original pollination-genetics-fitness studies and source."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_genotype_pollen_fitness_evidence_v39 as m

class StageEvidenceV39(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry=json.loads((ROOT/"data/external/anto_v39_original_genotype_pollen_fruit_studies_20261116.json").read_text())
        cls.access=json.loads((ROOT/"results/anto_v39_phlox2026_dryad_original_access_status_20261116.json").read_text())
        cls.out=m.review(cls.registry,cls.access)

    def test_four_papers_but_only_three_families_with_two_phlox_programs(self):
        self.assertEqual(self.out["n_original_published_study_programs"],4)
        self.assertEqual(self.out["n_distinct_current_plant_families"],3)
        self.assertEqual(self.out["n_original_Phlox_studies_same_genus"],2)
        self.assertTrue(self.out["original_papers_not_independent_phylogenetic_transition_events"])

    def test_antirrhinum_visitation_pollen_to_offspring_genetic_change_same_array(self):
        x=self.out["source_stage_matrix"]["antirrhinum_2001"]
        self.assertTrue(x["natural_pollinator_visits"])
        self.assertTrue(x["pollen_receipt_and_export"])
        self.assertTrue(x["offspring_allele_genotype_measured"])
        self.assertFalse(x["original_within_array_visits_pollen_offspring_chain"])
        self.assertTrue(x["original_within_1997_arrays_visitation_and_offspring_genotypes"])
        self.assertTrue(x["pollen_receipt_export_from_1994_1995_other_arrays"])
        self.assertFalse(x["original_flower_colour_gene_edited"])
        self.assertFalse(x["natural_multigeneration_fixation_measured"])

    def test_1997_offspring_allele_and_heterozygote_changes_from_author_arrays(self):
        a=self.out["Antirrhinum_1997_original_author_array_allele_shifts"]
        self.assertAlmostEqual(a["Rosy"]["sulfurea_delta_offspring_minus_initial"],-.06)
        self.assertAlmostEqual(a["Rosy"]["heterozygote_delta_offspring_minus_initial"],-.09)
        self.assertAlmostEqual(a["Barr"]["sulfurea_delta_offspring_minus_initial"],.02)
        self.assertAlmostEqual(a["Barr"]["heterozygote_delta_offspring_minus_initial"],-.07)
        self.assertFalse(a["Rosy"]["nectar_supplemented_on_yellow"])
        self.assertTrue(a["Barr"]["nectar_supplemented_on_yellow"])
        self.assertTrue(self.out["Antirrhinum_pollen_arrays_1994_1995_not_same_year_as_1997_offspring"])

    def test_forged_same_cohort_pollen_offspring_and_frequency_fails(self):
        d=copy.deepcopy(self.registry)
        d["studies"][0]["original_measurement_cohort_mapping"][
          "pollen_receipt_export_and_1997_offspring_not_exactly_same_array_and_flower"]=False
        with self.assertRaisesRegex(ValueError,"ANTIRRHINUM_FIELD_YEAR_COHORT_LINK_DRIFT"):
            m.review(d,self.access)
        d=copy.deepcopy(self.registry)
        d["studies"][0]["original_author_1997_offspring_frequencies"]["Barr"][
            "offspring_sulfurea_allele_frequency"]=0.72
        with self.assertRaisesRegex(ValueError,"ANTIRRHINUM_AUTHOR_1997_OFFSPRING_FREQUENCY_DRIFT"):
            m.review(d,self.access)

    def test_phlox_field_context_selection_not_regain(self):
        x=self.out["source_stage_matrix"]["phlox_2012_reinforcement"]
        self.assertTrue(x["field_selection_and_hybridization_reported"])
        self.assertTrue(x["pollinator_interspecific_movement"])
        self.assertFalse(x["allopatric_reinforcing_selection_detected"])
        self.assertFalse(x["natural_ancestral_colour_regain_observed"])
        self.assertTrue(self.out["sympatry_specific_fitness_selection_not_universal_colour_benefit"])

    def test_phlox_original_931_genotyped_seeds_and_dark_not_hue_selective_component(self):
        x=self.out["Phlox_2012_author_derived_dark_intensity_selection"]
        self.assertEqual(x["original_paternity_figure_seeds_genotyped_n"],931)
        self.assertAlmostEqual(x["dark_intensity_allele_s_author_inferred"],.32)
        self.assertEqual(x["hybrid_sterility_assumption_for_s"],.90)
        self.assertEqual(x["paternity_reference_hybridization_range"],[.28,.44])
        self.assertTrue(x["dark_intensity_hybridization_lower_by_more_than_half"])
        self.assertFalse(x["hue_red_allele_hybridization_effect_detected"])
        self.assertFalse(x["direct_female_fecundity_genotype_significant_difference_detected"])
        self.assertTrue(x["different_field_arrays_used_for_pollen_movement_and_seed_paternity"])

    def test_phlox_original_181_pollinators_2301_transitions_are_not_931_paternity_seeds(self):
        p=self.out["Phlox_2012_pollinator_source_count"]
        self.assertEqual(p["total_pollinator_individuals"],181)
        self.assertEqual(p["Battus_philenor"],108)
        self.assertEqual(p["skippers"],73)
        self.assertEqual(p["interplant_transitions"],2301)
        self.assertTrue(p["not_same_unit_as_paternity_genotyped_seeds"])

    def test_fake_red_hue_selection_or_doubled_paternity_sample_fails(self):
        for change in ("hue_red_allele_hybridization_effect_detected",
                       "hybridization_paternity_figure_n",
                       "author_derived_selection_coefficient_dark_allele"):
            d=copy.deepcopy(self.registry)
            pub=d["studies"][1]["original_author_2012_field_denominators"]
            pub[change]=True if change=="hue_red_allele_hybridization_effect_detected" else 1862 if change=="hybridization_paternity_figure_n" else .64
            with self.assertRaisesRegex(ValueError,"PHLOX_DARK_INTENSITY_NOT_RED_HUE_SELECTION_SOURCE_DRIFT"):
                m.review(d,self.access)

    def test_wahlenbergia_pollen_43pct_but_color_selection_nonpollen(self):
        x=self.out["source_stage_matrix"]["wahlenbergia_2013"]
        self.assertEqual(x["published_mean_pollen_supplementation_seed_increase_percent"],43)
        self.assertTrue(x["pollen_limited_seed_production"])
        self.assertTrue(x["selection_on_whiteness_stronger_with_excess_pollen"])
        self.assertFalse(x["floral_colour_selection_proved_pollen_receipt_mediated"])
        self.assertFalse(x["genotyped_flower_colour_locus"])

    def test_Phlox2026_PdPIRK_gene_pollen_tubes_fruit_but_DFR_only_color_reporter(self):
        x=self.out["source_stage_matrix"]["phlox_2026_si"]
        self.assertTrue(x["allele_specific_receptor_intervention_reported"])
        self.assertTrue(x["hand_cross_pollen_genotype_and_pollen_tubes"])
        self.assertTrue(x["controlled_cross_fruit_binary_outcome"])
        self.assertTrue(x["DFR_petal_bleaching_is_visible_silencing_reporter"])
        self.assertFalse(x["DFR_flower_colour_causally_changes_compatibility_demonstrated"])
        self.assertFalse(x["free_foraging_pollinator_choice_study"])

    def test_2026_PNAS_not_counted_again_as_2025_preprint(self):
        self.assertEqual(self.out["published_2026_PNAS_DOI"],"10.1073/pnas.2525299123")
        self.assertEqual(self.out["2025_bioRxiv_is_duplicate_not_new_paper"],True)
        self.assertEqual(self.out["n_original_published_study_programs"],4)

    def test_Dryad_401_no_raw_aggregate_effects_claimed(self):
        self.assertEqual(self.out["original_Phlox2026_Dryad_source_status"],"SOURCE_ACCESS_RESTRICTED_HTTP_401")
        self.assertEqual(self.out["raw_Phlox2026_fruit_pollen_tube_rows_reanalyzed"],0)
        self.assertIsNone(self.out["new_original_PdPIRK_allele_specific_fruit_risk_ratio"])
        self.assertTrue(self.out["technical_source_failure_not_biological_zero_effect"])

    def test_cross_study_numeric_effect_pool_prohibited(self):
        self.assertIsNone(self.out["cross_study_genetic_pollen_fitness_pooled_beta"])
        self.assertFalse(self.out["identical_gene_colour_treatment_replicated_across_all_studies"])
        self.assertEqual(self.out["n_same_AN2_colour_gene_to_wild_seedset_new_replicates"],0)

    def test_genetic_marker_not_assume_named_single_causal_nucleotide(self):
        self.assertFalse(self.out["original_Antirrhinum_colour_gene_causal_nucleotide_identified_v39"])
        self.assertTrue(self.out["one_locus_mendelian_colour_state_in_original_Antirrhinum"])

    def test_forged_colour_gene_edit_in_antirrhinum_fails(self):
        d=copy.deepcopy(self.registry)
        d["studies"][0]["measurements"]["genetically_engineered_flower_colour"]=True
        with self.assertRaisesRegex(ValueError,"ANTIRRHINUM_CAUSAL_EDIT_NOT_PRESENT"):
            m.review(d,self.access)

    def test_forged_pollen_mediated_colour_selection_in_wahlenbergia_fails(self):
        d=copy.deepcopy(self.registry)
        d["studies"][2]["selection_white_stronger_under_added_pollen"]=False
        with self.assertRaisesRegex(ValueError,"WAHLENBERGIA_NONPOLLEN_MECHANISM_SOURCE_DRIFT"):
            m.review(d,self.access)

    def test_forged_fruit_readout_for_DFR_colour_gene_fails(self):
        d=copy.deepcopy(self.registry)
        d["studies"][3]["original_gene_causality"]="DFR_CAUSES_FRUIT_SET"
        with self.assertRaisesRegex(ValueError,"PDPIRK_DFR_MECHANISM_CONFUSED"):
            m.review(d,self.access)

    def test_forged_numeric_reanalysis_through_401_fails(self):
        a=copy.deepcopy(self.access)
        a["row_level_fruit_pollen_tubes_reanalysis_completed"]=True
        with self.assertRaisesRegex(ValueError,"DRYAD_ACCESS_CLAIM_DRIFT"):
            m.review(self.registry,a)

    def test_wrong_peer_reviewed_PNAS_DOI_fails(self):
        d=copy.deepcopy(self.registry)
        d["studies"][3]["doi"]="10.1101/2025.06.17.660178"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_SOURCE_DOI_DRIFT"):
            m.review(d,self.access)

    def test_nonduplicate_Pholx2012_Dryad_microsat_not_fitness_data(self):
        self.assertTrue(self.out["different_Phlox_microsatellite_Dryad_not_fitness_reanalysis"])

    def test_no_same_original_genotype_across_multiple_papers_bridge_claim(self):
        self.assertFalse(self.out["all_five_stages_same_gene_plant_cohort_demonstrated"])
        self.assertTrue(self.out["subprocesses_identified_but_not_one_connected_experiment"])

    def test_author_published_43_percent_not_43_percent_pollen_increase(self):
        self.assertEqual(self.out["source_wahlenbergia_original_response"],"SEED_OUTPUT")
        self.assertEqual(self.out["author_published_43_percent"],43)

    def test_access_status_mutation_rejected(self):
        d=copy.deepcopy(self.access)
        d["original_Dryad_official_API_response"]="SOURCE_VERIFIED"
        with self.assertRaisesRegex(ValueError,"DRYAD_ACCESS_CLAIM_DRIFT"):
            m.review(self.registry,d)

    def test_parent_v38_Petunia_results_not_overridden(self):
        self.assertTrue(self.out["prior_Petunia_v38_separate_and_not_overridden"])
        self.assertFalse(self.out["confirmed_natural_ancestral_colour_regain_rate"])

if __name__=="__main__": unittest.main()
