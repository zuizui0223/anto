"""Guard published-research cases against fake replication or regain-rate claims."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from validate_afterlife_v11_empirical_regain_routes import validate,SOURCE

class RegainRouteEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(SOURCE.read_text(encoding="utf-8"))

    def test_original_source_contract(self):
        r=validate(self.data)
        self.assertEqual((r["n_cases"],r["n_distinct_clades"]),(6,5))
        self.assertEqual((r["ancestral_hue_restoration_cases"],r["altered_hue_regain_cases"]),(1,2))
        self.assertEqual((r["causal_regain_examples_same_Petunia_genus"],
                          r["causal_regain_independent_genus_count"]),(2,1))
        self.assertFalse(r["legitimate_cross_clade_frequency_inference"])

    def test_two_petunia_routes_are_not_two_independent_angiosperm_clades(self):
        r=validate(self.data)
        self.assertEqual(r["per_clade_case_counts"]["Petunia"],2)
        self.assertEqual(r["causal_regain_independent_genus_count"],1)
        self.assertTrue(r["source_cases_are_not_6_independent_regains"])

    def test_duplicate_case_fails(self):
        d=copy.deepcopy(self.data);d["cases"][0]["case_id"]=d["cases"][1]["case_id"]
        with self.assertRaisesRegex(ValueError,"uniquely"):validate(d)

    def test_iochroma_ancestry_not_promoted_to_causal_repair(self):
        d=copy.deepcopy(self.data)
        for c in d["cases"]:
            if c["taxon"]=="Iochroma edule":
                c["molecular_causal_evidence"]="SOURCE_SUPPORTED_FUNCTIONAL_MUTATION"
        with self.assertRaisesRegex(ValueError,"single causal-recovery clade"):validate(d)

    def test_antirrhinum_regulation_is_not_ancient_whole_flower_return(self):
        d=copy.deepcopy(self.data)
        for c in d["cases"]:
            if c["clade"]=="Antirrhinum":c["case_role"]="MOLECULARLY_SUPPORTED_PIGMENT_REGAIN"
        with self.assertRaises(ValueError):validate(d)

    def test_cirsium_leaf_source_not_phenotype_recovery(self):
        d=copy.deepcopy(self.data)
        for c in d["cases"]:
            if c["clade"]=="Cirsium":c["ancestral_hue_restoration"]="SOURCE_SUPPORTED"
        with self.assertRaisesRegex(ValueError,"Leaf-only"):validate(d)

    def test_pigment_restoration_not_same_as_any_pigment(self):
        d=copy.deepcopy(self.data)
        for c in d["cases"]:
            if c["taxon"]=="Petunia exserta":c["ancestral_hue_restoration"]="SOURCE_SUPPORTED"
        with self.assertRaisesRegex(ValueError,"Original hue"):validate(d)

    def test_disallow_published_cases_as_independent_fitness_effect(self):
        d=copy.deepcopy(self.data);d["cases"][0]["trait_specific_viable_seed_fitness"]="DIRECT_MEASURED"
        with self.assertRaisesRegex(ValueError,"fitness"):validate(d)

    def test_no_invented_longitudinal_trait_recovery_rate(self):
        d=copy.deepcopy(self.data);d["cases"][0]["quantitative_lineage_reexpression_rate"]=0.9
        with self.assertRaisesRegex(ValueError,"cannot estimate"):validate(d)

    def test_no_unsupported_literature_case(self):
        d=copy.deepcopy(self.data);d["cases"][0]["source_dois"]=[]
        with self.assertRaisesRegex(ValueError,"Missing source"):validate(d)

if __name__=="__main__":
    unittest.main()
