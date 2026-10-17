import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from anto_colonization import Parameters, Plant, simulate, visible, first_pigmented, order_experiment_probability
from anto_identifiability import benchmark, generate

class TestDispersalReexpression(unittest.TestCase):
    def test_invalid_probabilities_and_capacity(self):
        with self.assertRaises(ValueError):
            Parameters(('main', 'island'), {'main': 10, 'island': 10},
                       {('main', 'island'): 1.1})
        with self.assertRaises(ValueError):
            Parameters(('main',), {'main': 10}, reexpression=.8,
                       degradation={'main': .4})
        with self.assertRaises(ValueError):
            simulate(Parameters(('main',), {'main': 1}),
                     {'main': [Plant('A', 'P'), Plant('B', 'W+')]}, 3)

    def test_reexpression_requires_competence(self):
        p = Parameters(('island',), {'island': 40}, fecundity=5, reexpression=1)
        a = simulate(p, {'island': [Plant('A', 'W+')]}, 1)
        b = simulate(p, {'island': [Plant('A', 'W-')]}, 1)
        self.assertGreater(a['histories'][1]['regions']['island']['pigmented'], 0)
        self.assertEqual(b['histories'][1]['regions']['island']['pigmented'], 0)
        self.assertGreater(a['events'][0]['regions']['island'].get('intrinsic_reexpression', 0), 0)
        self.assertNotIn('intrinsic_reexpression', b['events'][0]['regions']['island'])

    def test_immigration_without_intrinsic_reexpression(self):
        p = Parameters(('source', 'island'), {'source': 40, 'island': 40},
                       seed_migration={('source', 'island'): 1},
                       fecundity=6, reexpression=0)
        result = simulate(p, {'source': [Plant('S', 'P')],
                              'island': [Plant('I', 'W-')]}, 1)
        self.assertEqual(first_pigmented(result, 'island'), 1)
        self.assertGreater(result['events'][0]['regions']['island'].get('pigmented_immigrant', 0), 0)
        self.assertEqual(result['events'][0]['regions']['island'].get('intrinsic_reexpression', 0), 0)

    def test_pollen_introgression_without_seed_migration(self):
        p = Parameters(('source', 'island'), {'source': 50, 'island': 50},
                       pollen_delivery={('source', 'island'): 1},
                       fecundity=10, introgression=1, reexpression=0)
        result = simulate(p, {'source': [Plant('S', 'P')],
                              'island': [Plant('I', 'W-')]}, 1)
        self.assertGreater(result['events'][0]['regions']['island'].get('introgressive_reacquisition', 0), 0)
        self.assertEqual(result['events'][0]['regions']['island'].get('seed_immigrant', 0), 0)
        self.assertTrue(any(p['introgressed'] and p['state'] == 'P'
                            for p in result['final']['island']))

    def test_empty_island_can_be_colonized(self):
        p = Parameters(('source', 'island'), {'source': 100, 'island': 100},
                       seed_migration={('source', 'island'): 1}, fecundity=20)
        result = simulate(p, {'source': [Plant('S', 'W+')], 'island': []}, 1)
        self.assertGreater(result['events'][0]['regions']['island']['seed_recolonist'], 0)
        self.assertEqual(result['histories'][1]['regions']['island']['pigmented'], 0)

    def test_capacity(self):
        p = Parameters(('A',), {'A': 5}, fecundity=12)
        result = simulate(p, {'A': [Plant('A', 'P')]*5}, 5)
        self.assertTrue(all(sum(h['regions']['A'][k] for k in ('pigmented','white')) <= 5
                            for h in result['histories']))

    def test_same_visible_different_future(self):
        self.assertEqual(set(visible({'A':[Plant('A','W+')]})['A']), {'pigmented', 'white'})
        self.assertEqual(visible({'A': [Plant('A', 'W+')]} )['A']['white'],
                         visible({'A': [Plant('A', 'W-')]} )['A']['white'])
        p = Parameters(('A',), {'A': 100}, fecundity=10, reexpression=1)
        a = simulate(p, {'A': [Plant('A', 'W+')]}, 1)
        b = simulate(p, {'A': [Plant('A', 'W-')]}, 1)
        self.assertGreater(a['histories'][-1]['regions']['A']['pigmented'], 0)
        self.assertEqual(b['histories'][-1]['regions']['A']['pigmented'], 0)

    def test_disperse_then_degrade_differs(self):
        p = order_experiment_probability(.02, .2, delay=12, horizon=8,
                                         r_island=.1, n_founders=8)
        self.assertGreater(p['late'], p['early'])
        p = order_experiment_probability(.2, .2, 12, 8, .1, 8)
        self.assertAlmostEqual(p['late'], p['early'])

    def test_no_reexpression_if_no_ability_or_time(self):
        self.assertEqual(order_experiment_probability(.2, .2, 10, 0, .1)['early'], 0)
        self.assertEqual(order_experiment_probability(.2, .2, 10, 5, 0)['early'], 0)

class TestIdentifiability(unittest.TestCase):
    def test_equivalent_visible_mechanism(self):
        rows = generate(15)
        self.assertEqual(len(rows), 45)
        self.assertEqual({r[0] for r in rows}, {'intrinsic','immigration','introgression'})

    def test_assays_improve_synthetic_discrimination(self):
        b = benchmark()
        self.assertAlmostEqual(b['predictions']['visible_only']['accuracy'], 1/3)
        self.assertGreater(b['predictions']['two_molecular_markers']['accuracy'], .8)

if __name__ == '__main__':
    unittest.main()
