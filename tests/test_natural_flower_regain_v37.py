"""v37 original-published functional/phylogenetic/behavioural evidence gates."""
from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_natural_flower_regain_v37 as v


class OriginalPublishedRegainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads((ROOT/"data/external/anto_v37_natural_regain_and_native_reactivation_original_sources_20261116.json").read_text())

    def test_four_original_DOIs_two_families_one_petunia_radiation(self):
        x=v.review(self.source)
        self.assertEqual(x["n_original_publications"],4)
        self.assertEqual(x["n_independent_botanical_families"],2)
        self.assertEqual(x["n_distinct_petunia_radiations"],1)
        self.assertEqual(len(x["registered_original_DOIs"]),4)
        self.assertTrue(x["related_Petunia_studies_not_independent_macro_evolutionary_replicates"])

    def test_two_distinct_natural_regain_routes_in_one_original_radiation(self):
        x=v.review(self.source)
        self.assertEqual(x["natural_regain_source_species"],["Petunia exserta","Petunia secreta"])
        self.assertEqual(x["natural_regain_mechanism_count"],2)
        self.assertEqual(x["natural_regain_source_study_count"],2)
        self.assertEqual(x["natural_regain_route_by_species"]["Petunia secreta"],"AN2_FRAME_RESTORATION")
        self.assertEqual(x["natural_regain_route_by_species"]["Petunia exserta"],
                         "PARALOG_DPL_REDEPLOYMENT_WITH_OTHER_PATHWAY_CHANGES")
        self.assertFalse(x["observed_in_real_time_ancestral_to_descendant_reversion"])

    def test_secreta_two_base_pair_restores_AN2_not_DPL(self):
        x=v.review(self.source)
        self.assertEqual(x["natural_AN2_repair"]["indel"],"2_BASE_PAIR_COMPENSATORY_DELETION")
        self.assertEqual(x["natural_AN2_repair"]["gene"],"AN2")
        self.assertEqual(x["natural_AN2_repair"]["state"],"READING_FRAME_RESTORED")

    def test_exserta_retains_nonfunctional_AN2_dpl_and_additional_hydroxylation(self):
        x=v.review(self.source)
        self.assertTrue(x["natural_DPL_paralog_redeployment"]["AN2_still_nonfunctional"])
        self.assertEqual(x["natural_DPL_paralog_redeployment"]["regulator"],"DPL")
        self.assertIn("AAT",x["natural_DPL_paralog_redeployment"]["other_pathway_changes"])
        self.assertFalse(x["natural_DPL_paralog_redeployment"]["red_colour_assumed_pelargonidin"])
        self.assertTrue(x["natural_DPL_paralog_redeployment"]["native_DPL_silencing_white_sectors"])

    def test_three_original_native_white_petal_complementation_programs_two_families(self):
        x=v.review(self.source)
        self.assertEqual(x["n_original_native_white_petal_complementation_studies"],3)
        self.assertEqual(x["n_families_with_native_white_petal_complementation"],2)
        self.assertEqual(x["n_stable_native_white_petal_complementation_studies"],2)
        self.assertEqual(x["n_transient_local_native_white_petal_complementation_studies"],1)
        self.assertEqual(len(x["native_white_petal_restorations"]),3)

    def test_DPL_9_transgenic_lines_not_9_HPLC_plant_replicates(self):
        x=v.review(self.source)
        d=x["Petunia_axillaris_DPL_petal_rescue"]
        self.assertEqual(d["independent_lines_created"],9)
        self.assertEqual(d["lines_assayed_for_chemistry"],2)
        self.assertEqual(d["HPLC_transgenic_plants"],6)
        self.assertEqual(d["HPLC_control_plants"],3)
        self.assertEqual(d["transgenic_plants_per_assayed_line"],3)
        self.assertTrue(d["corolla_pigmentation_restored"])

    def test_AN2_transgenic_host_P_axillaris_not_historical_secreta(self):
        x=v.review(self.source)
        t=x["Petunia_axillaris_AN2_petal_rescue"]
        self.assertEqual(t["original_white_host"],"Petunia axillaris")
        self.assertEqual(t["gene_donor"],"Petunia integrifolia")
        self.assertEqual(t["regulator"],"AN2")
        self.assertFalse(t["spontaneous_ancestral_pseudogene_repair"])

    def test_Actinidia_transient_recovery_in_white_petal_local_with_GFP_negative(self):
        x=v.review(self.source)
        a=x["Actinidia_white_hybrid_MYB110a_petal_rescue"]
        self.assertEqual(a["host_white_genotype"],"21-05-09e")
        self.assertEqual(a["introduced_gene"],"MYB110a")
        self.assertTrue(a["local_cells_pigmented"])
        self.assertTrue(a["GFP_control_negative"])
        self.assertFalse(a["stable_intergenerational_inheritance_established"])

    def test_greenhouse_moth_and_bee_preferences_opposite_and_no_seed_set(self):
        x=v.review(self.source)
        p=x["original_pollinator_preference"]
        self.assertEqual(p["original_study_count"],1)
        self.assertEqual(p["Manduca_first_choose_white"],30)
        self.assertEqual(p["Manduca_first_choose_coloured"],6)
        self.assertEqual(p["Manduca_n"],36)
        self.assertEqual(p["Bombus_n_trials"],12)
        self.assertEqual(p["Bombus_direction"],"AN2_COLOURED_P_AXILLARIS")
        self.assertFalse(p["field_seed_set_or_offspring_fitness_measured"])

    def test_original_moth_30_6_not_replaced_by_pooled_trials(self):
        d=copy.deepcopy(self.source)
        p=next(s["controlled_pollinator_preference"] for s in d["original_study_programs"]
               if s["id"]=="Hoballah2007_AN2_pollinators")
        p["moth_first_choice_white_P_axillaris"]=29
        with self.assertRaisesRegex(ValueError,"POLLINATOR_SOURCE_DENOMINATOR_DRIFT"):
            v.review(d)

    def test_bombus_preference_not_equated_with_moth_pollen_transfer(self):
        x=v.review(self.source)
        self.assertFalse(x["measured_natural_outcrossing_advantage"])
        self.assertFalse(x["pollinator_choices_equivalent_to_seed_fitness"])
        self.assertFalse(x["combined_causal_genotype_to_wild_fitness_across_programs"])

    def test_forged_restored_exserta_AN2_fails(self):
        d=copy.deepcopy(self.source)
        d["original_study_programs"][1]["natural_colour_regain"]["AN2_still_nonfunctional"]=False
        with self.assertRaisesRegex(ValueError,"DPL_AND_AN2_STATE_DRIFT"):
            v.review(d)

    def test_fake_natural_live_reversion_fails(self):
        d=copy.deepcopy(self.source)
        d["original_study_programs"][0]["natural_colour_regain"]["direct_historical_transition_observed_in_real_time"]=True
        with self.assertRaisesRegex(ValueError,"HISTORICAL_DIRECT_OBSERVATION_FORGED"):
            v.review(d)

    def test_fake_DPL_HPLC_denominator_fails(self):
        d=copy.deepcopy(self.source)
        d["original_study_programs"][1]["original_white_petal_transformation"]["HPLC_transgenic_samples"]=9
        with self.assertRaisesRegex(ValueError,"DPL_TRANSFORMATION_UNIT_DRIFT"):
            v.review(d)

    def test_fake_Actinidia_stable_F1_genetic_transformation_fails(self):
        d=copy.deepcopy(self.source)
        d["original_study_programs"][3]["original_white_petal_transformation"]["long_term_stable_transgenic_plant_established"]=True
        with self.assertRaisesRegex(ValueError,"ACTINIDIA_LOCAL_TRANSIENT_ONLY"):
            v.review(d)

    def test_fake_four_evolutionary_independent_radiations_fails(self):
        d=copy.deepcopy(self.source)
        d["original_programs_not_independent_macro_radiations"]=False
        # No corresponding declaration: the original frozen constraint remains
        self.assertTrue(v.review(self.source)["related_Petunia_studies_not_independent_macro_evolutionary_replicates"])
        d=copy.deepcopy(self.source)
        d["original_study_programs"][1]["shared_radiation"]="InventedDifferentPetuniaClade"
        with self.assertRaisesRegex(ValueError,"SOURCE_SHARED_PHYLOGENETIC_RADIATION_DRIFT"):
            v.review(d)

    def test_regain_rates_not_invented(self):
        x=v.review(self.source)
        self.assertIsNone(x["cross_source_pooled_regain_effect"])
        self.assertIsNone(x["natural_regain_probability"])
        self.assertEqual(x["independent_original_raw_sequence_reanalysis_count"],0)
        self.assertFalse(x["all_white_flowering_plants_functionally_restorable"])

    def test_do_not_forge_full_ANTO_v33_three_genes_external_model(self):
        x=v.review(self.source)
        self.assertFalse(x["exact_v33_F3Hprime_F3H5prime_FLS_RMSE_externally_replayed"])
        self.assertTrue(x["original_assays_and_evolutionary_history_are_distinct_grains"])

if __name__=="__main__":
    unittest.main()
