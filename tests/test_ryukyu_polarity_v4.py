import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from anto_ryukyu_polarity import transition_matrix, posterior_root_pigmented, sensitivity


class RyukyuPolarityTests(unittest.TestCase):
    def test_transition_probs_valid(self):
        for g,l in ((.1,.2),(.5,.5),(2.5,.5),(0,1)):
            P=transition_matrix(g,l,1)
            for row in P:
                self.assertAlmostEqual(sum(row),1.0,places=12)
                self.assertTrue(all(0<=p<=1 for p in row))

    def test_ancestral_polarity_not_identified_by_symmetric_two_tips(self):
        # 3 + 3 individuals sampled within two species do NOT equal
        # three independent white-to-pigmented origin events.
        for t in (.71,.93,1.33,100):
            out=posterior_root_pigmented(.5,.5,t,t)
            self.assertAlmostEqual(out["posterior_ancestral_pigmented"],.5,places=12)
            self.assertAlmostEqual(out["likelihood_if_white_ancestor"],
                                   out["likelihood_if_pigmented_ancestor"],places=12)

    def test_posterior_equals_prior_if_likelihoods_symmetric(self):
        out=posterior_root_pigmented(.5,.5,.93,.93,prior_p=.2)
        self.assertAlmostEqual(out["posterior_ancestral_pigmented"],.2,places=12)

    def test_prior_rate_changes_polarity(self):
        a=posterior_root_pigmented(.1,.5,.93,.93)
        b=posterior_root_pigmented(2.5,.5,.93,.93)
        self.assertNotAlmostEqual(a["posterior_ancestral_pigmented"],
                                  b["posterior_ancestral_pigmented"],places=5)

    def test_non_ultrametric_terminal_branches_can_break_symmetry(self):
        r=posterior_root_pigmented(.5,.5,.25,1.0)
        self.assertNotAlmostEqual(r["posterior_ancestral_pigmented"],.5,places=5)

    def test_zero_mutation_contrast_is_impossible(self):
        with self.assertRaises(ValueError): posterior_root_pigmented(0,0,.93,.93)
        with self.assertRaises(ValueError): posterior_root_pigmented(.2,.5,-1,.93)

    def test_grid_is_explicitly_hypothetical(self):
        r=sensitivity()
        self.assertEqual(len(r["cases"]),12)
        self.assertTrue(r["not_phenotype_transition_date"])
        self.assertTrue(r["rates_hypothetical"])


if __name__=="__main__":
    unittest.main()
