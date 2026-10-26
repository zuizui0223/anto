"""Evidence boundary tests for the original 161 Penstemon F2 observations."""
import sys
import json
import unittest
from pathlib import Path
from copy import deepcopy
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_penstemon_2014_afterlife_contrast_v1 import audit

class PenstemonMechanisticContrastTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=json.loads((ROOT/"data/external/penstemon_2014_f3p5h_source_genetic_accessibility_v1.json").read_text())
    def test_full_original_cross_contingency(self):
        r=audit(self.source)
        self.assertEqual(r["observed_F2_individuals"],161)
        self.assertEqual(r["source_genotype_denominators"],{"FF":57,"Ff":80,"ff":24})
        self.assertEqual(r["source_purple_by_genotype"],{"FF":57,"Ff":65,"ff":0})
        self.assertEqual(r["source_red_by_genotype"],{"FF":0,"Ff":15,"ff":24})
        self.assertEqual({k:F(v) for k,v in r["blue_purple_probability_conditional_on_genotype_fraction"].items()},{"FF":F(1),"Ff":F(13,16),"ff":F(0)})
    def test_hybrid_donor_not_sufficient_among_heterozygotes(self):
        r=audit(self.source)
        self.assertEqual(F(r["heterozygote_red_probability"]),F(3,16))
        self.assertEqual(r["functional_donor_F_allele_but_still_red_Ff"],15)
        self.assertTrue(r["difference_native_repair_vs_hybrid_inheritance_is_mechanistically_necessary"])
    def test_structural_native_gene_failure_and_external_donor_are_separate(self):
        r=audit(self.source)
        self.assertTrue(r["source_internal_F3p5H_repair_of_single_stop_failed_to_restore_in_vitro_activity"])
        self.assertTrue(r["external_functional_donor_allele_can_produce_blue_purple_in_hybrid_offspring"])
        self.assertFalse(r["natural_introgression_frequency_or_permanence_estimated"])
        self.assertFalse(r["synthetic_ANTO_model_parameter_g_inferred"])
        self.assertFalse(r["plant_phenotype_after_regain_fitness_measured_here"])
    def test_one_person_total_mismatch_is_rejected(self):
        x=deepcopy(self.source);x["observed_F2_table"]["purple"]["FF"]=56
        with self.assertRaisesRegex(ValueError,"denominator"):audit(x)
    def test_one_transposed_cell_even_total_preserved_is_rejected(self):
        x=deepcopy(self.source)
        x["observed_F2_table"]["purple"]["FF"]=56
        x["observed_F2_table"]["purple"]["Ff"]=66
        with self.assertRaises(ValueError):audit(x)
    def test_phenotypic_loss_of_blue_not_complete_flower_depigmentation(self):
        x=deepcopy(self.source)
        x["trait_target"]="total flower color absent"
        with self.assertRaisesRegex(ValueError,"total floral pigment loss"):audit(x)
    def test_no_false_natural_gene_flow_from_cross(self):
        x=deepcopy(self.source)
        x["reported_source_molecular_results"]["unmeasured_actual_wild_gene_flow"]=False
        with self.assertRaises(ValueError):audit(x)
    def test_original_author_is_not_user_discovered_source(self):
        x=deepcopy(self.source)
        x["user_discovered_new_empirical_evolutionary_transition"]=True
        with self.assertRaisesRegex(ValueError,"user's new experiment"):audit(x)

if __name__=="__main__":unittest.main()
