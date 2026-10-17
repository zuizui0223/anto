import unittest
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from anto_model import generator, distribution, observable_colour, conditional_reexpression

class TestModel(unittest.TestCase):
    def test_conservation(self):
        q=generator(["A","B"],{("A","B"):0.4},0.1,0.2,{"A":0.1,"B":0.8})
        self.assertTrue(np.allclose(q.sum(axis=1),0))
        p=distribution(q,[0,1,0,0,0,0],12)
        self.assertAlmostEqual(p.sum(),1,places=10)
        self.assertGreaterEqual(min(p),-1e-12)
    def test_visible_cannot_identify_competence(self):
        q=generator(["A"],{},0,0.3,0.2)
        a=distribution(q,[0,1,0],0)
        b=distribution(q,[0,0,1],0)
        self.assertEqual(observable_colour(a,["A"])["A"]["white"],
                         observable_colour(b,["A"])["A"]["white"])
        self.assertNotEqual(observable_colour(a,["A"])["A"]["competent_white"],
                            observable_colour(b,["A"])["A"]["competent_white"])
    def test_window_decreases(self):
        self.assertGreater(conditional_reexpression(0,10,0.1,0.2),
                           conditional_reexpression(20,10,0.1,0.2))
        self.assertAlmostEqual(conditional_reexpression(0,10,0.1,0),
                               conditional_reexpression(20,10,0.1,0))
    def test_noncommuting_regional_degradation(self):
        d=generator(["A","B"],{("A","B"):0.4},0,0,0)
        g=generator(["A","B"],{},0,0,{"A":0.01,"B":0.5})
        self.assertGreater(np.max(np.abs(d@g-g@d)),0.01)

if __name__=="__main__": unittest.main()
