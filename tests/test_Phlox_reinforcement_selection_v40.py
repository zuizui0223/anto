"""ANTO v40: original-source selection equation and counterfactual-cost guards."""
import copy
import json
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_Phlox_reinforcement_selection_v40 as m


class ConditionalSelectionV40(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry=json.loads((ROOT/"data/external/anto_v40_original_Phlox_selection_2012_2014_author_registry_20261116.json").read_text())
        cls.result=m.review(cls.registry)

    def test_independent_original_papers_but_only_one_radiation(self):
        self.assertEqual(self.result["n_original_published_papers"],2)
        self.assertEqual(self.result["n_independent_phylogenetic_radiations"],1)
        self.assertFalse(self.result["two_papers_counted_as_independent_phylogenetic_reacquisitions"])

    def test_source_original_0_43_0_21_0_90_exact_reconstruction(self):
        a=self.result["reproduced_original_2012_author_selection"]
        self.assertAlmostEqual(a["light_viable_fitness_proxy"],.613)
        self.assertAlmostEqual(a["dark_viable_fitness_proxy"],.811)
        self.assertAlmostEqual(a["source_relative_selection_coefficient"],.32300163132137016)
        self.assertAlmostEqual(a["hybridization_reduction_fraction"],.5116279069767442)
        self.assertAlmostEqual(round(a["source_relative_selection_coefficient"],2),.32)
        self.assertFalse(a["actual_net_population_allele_frequency_increase_measured"])

    def test_no_hybrid_sterility_no_reinforcement_gain(self):
        self.assertEqual(m.hybrid_advantage(.43,.21,0),0)
        self.assertAlmostEqual(m.hybrid_advantage(.43,.21,.9),.32300163132137016)

    def test_no_hybridization_difference_no_reinforcement_gain(self):
        self.assertEqual(m.hybrid_advantage(.43,.43,.9),0)
        self.assertEqual(m.net_selection(.43,.43,.9,0),0)
        self.assertLess(m.net_selection(.43,.43,.9,.10),0)

    def test_focal_has_higher_hybridization_advantage_reverses(self):
        self.assertLess(m.hybrid_advantage(.21,.43,.9),0)

    def test_invalid_probabilities_and_zero_denominator_fail(self):
        for value in (-.01,1.01,float("nan"),float("inf")):
            with self.assertRaises(ValueError):m.hybrid_advantage(value,.21,.9)
            with self.assertRaises(ValueError):m.hybrid_advantage(.43,.21,value)
            with self.assertRaises(ValueError):m.net_selection(.43,.21,.9,value)
        with self.assertRaisesRegex(ValueError,"ZERO_REFERENCE_VIABLE_FITNESS"):
            m.hybrid_advantage(1,.2,1)

    def test_2012_q_sensitivity_at_0_half_90_and_one(self):
        grid={x["hybrid_sterility_fraction"]:x["modelled_selection"] for x in
              self.result["source_parameter_q_sensitivity"]}
        self.assertAlmostEqual(grid[0],0)
        self.assertAlmostEqual(grid[.5],.14012738853503184)
        self.assertAlmostEqual(grid[.9],.32300163132137016)
        self.assertAlmostEqual(grid[1],.38596491228070173)
        self.assertEqual(list(sorted(grid)),[0,.25,.5,.75,.9,1])

    def test_break_even_cost_formula(self):
        c=m.break_even_nonhybrid_cost(.43,.21,.9)
        self.assertAlmostEqual(c,.198/.811)
        self.assertAlmostEqual(m.net_selection(.43,.21,.9,c),0,places=13)
        self.assertGreater(m.net_selection(.43,.21,.9,.2),0)
        self.assertLess(m.net_selection(.43,.21,.9,.25),0)
        self.assertEqual(m.break_even_nonhybrid_cost(.43,.21,0),0)

    def test_required_q_for_hypothetical_cost(self):
        self.assertAlmostEqual(m.required_sterility_for_cost(.43,.21,.1),.1/.241)
        self.assertAlmostEqual(m.required_sterility_for_cost(.43,.21,.2),.2/.262)
        self.assertGreater(m.required_sterility_for_cost(.43,.21,.3),1)
        self.assertEqual(m.required_sterility_for_cost(.43,.21,0),0)

    def test_countervailing_cost_grid_all_cells_and_original_model_limit(self):
        table=self.result["q_by_hypothetical_dark_nonhybrid_cost_grid"]
        self.assertEqual(len(table),6*7)
        rec=lambda q,c:next(x for x in table if
            x["hybrid_sterility_fraction"]==q and x["hypothetical_dark_nonhybrid_cost_fraction"]==c)
        self.assertAlmostEqual(rec(.9,0)["relative_net_selection"],.32300163132137016)
        self.assertGreater(rec(.9,.20)["relative_net_selection"],0)
        self.assertLess(rec(.9,.25)["relative_net_selection"],0)
        self.assertAlmostEqual(rec(0,.2)["relative_net_selection"],-.2)
        self.assertTrue(all(x["not_empirical_confidence_interval"] for x in table))

    def test_allopatric_original_arrays_bootstrap_source_not_binomial(self):
        a=self.result["independent_2014_author_allopatry_pollinator_preferences"]
        self.assertEqual(a["total_original_pollinators_reported"],278)
        self.assertEqual(a["total_original_visits_reported"],1000)
        self.assertEqual(a["sum_of_three_table4_array_N"],276)
        self.assertTrue(a["author_CIs_not_recalculated_from_independent_binomial_visits"])
        self.assertEqual([round(x["light_blue_visit_share"],3) for x in a["arrays"]],
                         [.619,.603,.682])
        self.assertEqual([x["author_original_array_N"] for x in a["arrays"]],[84,118,74])
        self.assertTrue(all(x["author_reported_bootstrap_CI_excludes_50pct"]
                            for x in a["arrays"]))

    def test_2014_preference_not_converted_to_2012_hybrid_paternity(self):
        self.assertFalse(self.result["allopatry_pollinator_visit_share_is_hybridization_rate"])
        self.assertFalse(self.result["2014_allopatry_direct_dark_paternal_siring_success_measured"])
        self.assertFalse(self.result["2014_allopatry_net_colour_selection_same_as_2012_s"])
        self.assertTrue(self.result["2014_2010_2011_net_predictions_depend_on_fecundity_model"])

    def test_negative_hue_locus_cannot_be_assigned_dark_selection(self):
        self.assertTrue(self.result["original_2012_hue_red_locus_no_detected_hybridization_reduction"])
        self.assertFalse(self.result["original_2012_derived_s_transferred_to_red_hue_allele"])
        d=copy.deepcopy(self.registry)
        d["source_papers"][0]["original_measurements"]["hue_red_F3prime5primeH_locus_hybridization_reduction_detected"]=True
        with self.assertRaisesRegex(ValueError,"HUE_LOCUS_NULL_GUARD_DRIFT"):
            m.review(d)

    def test_forged_2012_hybrid_fraction_or_original_assumption_fails(self):
        for key,val in (("light_reference_hybrid_paternity_fraction",.35),
                        ("dark_focal_hybrid_paternity_fraction",.29),
                        ("hybrid_sterility_assumed_from_earlier_studies",.99)):
            d=copy.deepcopy(self.registry)
            d["source_papers"][0]["original_measurements"][key]=val
            with self.assertRaisesRegex(ValueError,"SOURCE_2012_PARAMETRIC_INPUT_DRIFT"):
                m.review(d)

    def test_forged_original_2014_bootstrap_interval_fails(self):
        d=copy.deepcopy(self.registry)
        d["source_papers"][1]["original_measurements"]["arrays"][1]["author_bootstrap_95ci"]=[.7,.9]
        with self.assertRaisesRegex(ValueError,"ORIGINAL_2014_TABLE4_DRIFT"):
            m.review(d)

    def test_missing_source_and_duplicate_author_program_rejected(self):
        d=copy.deepcopy(self.registry)
        d["source_papers"].append(d["source_papers"][0])
        with self.assertRaisesRegex(ValueError,"INDEPENDENT_ORIGINAL_SOURCE_SET_DRIFT"):
            m.review(d)
        d=copy.deepcopy(self.registry)
        d["source_papers"][1]["doi"]="10.0000/fake"
        with self.assertRaisesRegex(ValueError,"SOURCE_ORIGINAL_DOI_DRIFT"):
            m.review(d)

    def test_no_reanalysis_of_2014_Dryad_or_seed_rows_and_no_numeric_meta(self):
        a=self.result
        self.assertFalse(a["2014_Dryad_original_individual_rows_reanalyzed"])
        self.assertFalse(a["2012_paternity_931_seed_rows_reanalyzed"])
        self.assertTrue(a["numerical_sensitivity_scenarios_not_independent_observed_genotypes"])
        self.assertIsNone(a["original_cross_paper_quantitative_meta_selection_effect"])
        self.assertIsNone(a["original_2014_allopatric_dark_allele_cost_empirically_measured"])
        self.assertIsNone(a["observed_allele_fixation_rate"])

if __name__=="__main__":unittest.main()
