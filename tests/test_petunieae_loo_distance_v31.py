"""Source-grain and numerical adversarial checks; no real outcomes in fixtures."""
import copy
import math
import sys
import unittest
from dataclasses import dataclass
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import analyze_petunieae_loo_distance_v31 as m
except ModuleNotFoundError:
    m = None

@dataclass
class Node:
    label: str = None
    children: tuple = ()
    length: float = None


def source(n=12):
    samples = {}
    for i in range(n):
        samples[f'S{i:02}'] = {'genes': {'F3primeH': i/10, 'F35H': (i%3)/4, 'FLS': (i%5)/3},
                              'responses': {'anthocyanidins_log1p_total': 1+i/9,
                                            'flavonols_log1p_total': 2+i/8}}
    distances = {a: {b: 0.0 if a == b else float(1+abs(int(a[1:])-int(b[1:])))
                      for b in samples} for a in samples}
    return samples, distances


def prediction(train, test, endpoint, genes):
    assert 'responses' not in test, 'heldout target was exposed'
    avg = sum(r['responses'][endpoint] for r in train.values()) / len(train)
    return avg if not genes else avg + 0.1*test['genes']['FLS']

class V31(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m, 'v31 implementation not yet present')

    def test_declared_scope(self):
        self.assertEqual(m.TOLERANCES, (1e-6, 1e-5, 1e-4))
        self.assertEqual(m.PRIMARY_TOL, 1e-5)
        self.assertEqual(m.BRANCH_GENES, ('F3primeH', 'F35H', 'FLS'))

    def test_decimal_tree_path_and_root_stem(self):
        tree = Node(children=(Node(children=(Node('A',length=.1), Node('B',length=.2)),length=.3),
                              Node('C',length=.4)), length=999)
        d = m.tree_distances(tree)
        self.assertEqual(d['A']['B'], .3)
        self.assertEqual(d['A']['C'], .8)
        self.assertEqual(d['B']['C'], .9)
        self.assertEqual(d['C']['A'], .8)
        self.assertEqual(d['A']['A'], 0)

    def test_bad_edges_and_duplicate_tips(self):
        for bad in (None, -1, float('nan'), float('inf')):
            t = Node(children=(Node('A',length=bad), Node('B',length=1)))
            with self.assertRaises(ValueError): m.tree_distances(t)
        with self.assertRaises(ValueError):
            m.tree_distances(Node(children=(Node('A',length=1), Node('A',length=2))))

    def test_zero_length_edge_allowed_but_all_zero_not(self):
        d = m.tree_distances(Node(children=(Node('A',length=0), Node('B',length=1))))
        self.assertEqual(d['A']['B'],1)
        with self.assertRaises(ValueError):
            m.tree_distances(Node(children=(Node('A',length=0), Node('B',length=0))))

    def test_complete_link_bands_not_transitive_chaining(self):
        b = m.distance_bands({'A':1, 'B':1.000006, 'C':1.000012},1e-5)
        self.assertEqual(b, [['A','B'], ['C']])
        self.assertEqual(m.band_ranks({'A':1, 'B':1.000006, 'C':1.000012},1e-5),
                         {'A':1.5,'B':1.5,'C':3.0})

    def test_nonfinite_distance_and_invalid_tolerance_rejected(self):
        for val in (-1, float('nan'), float('inf')):
            with self.assertRaises(ValueError): m.distance_bands({'A':val},1e-5)
        for tol in (0,-1,float('nan')):
            with self.assertRaises(ValueError): m.distance_bands({'A':1},tol)

    def test_tiers_do_not_split_equal_distance_species(self):
        gaps={'A':1,'B':1,'C':1,'D':2,'E':3,'F':4}
        tiers=m.tier_ids(gaps,1e-5)
        positions={x:k for k,v in tiers.items() for x in v}
        self.assertEqual(len({positions[x] for x in ('A','B','C')}),1)
        self.assertEqual(set(positions),set(gaps))
        t=m.tier_ids({x:1 for x in gaps},1e-5)
        self.assertEqual(t['near'],[])
        self.assertEqual(t['far'],[])
        self.assertEqual(len(t['middle']),6)

    def test_midrank_and_constant_correlation(self):
        self.assertEqual(m.midrank([1,1,3]),[1.5,1.5,3])
        self.assertAlmostEqual(m.spearman([1,2,3],[3,2,1]),-1)
        self.assertIsNone(m.spearman([1,1,1],[2,3,4]))

    def test_neighbor_average_uses_all_source_precision_ties(self):
        train={'A':{'responses':{'e':2}}, 'B':{'responses':{'e':8}}, 'C':{'responses':{'e':100}}}
        value, ids=m.neighbor_prediction(train,{'A':1,'B':1.000002,'C':1.01},'e',1e-5)
        self.assertEqual(value,5)
        self.assertEqual(ids,['A','B'])
        value,ids=m.neighbor_prediction(train,{'A':1,'B':1.000002,'C':1.01},'e',1e-6)
        self.assertEqual(value,2)
        self.assertEqual(ids,['A'])

    def test_neighbor_mapping_cannot_drop_a_training_species(self):
        with self.assertRaises(ValueError):
            m.neighbor_prediction({'A':{'responses':{'e':1}}},{'B':1},'e',1e-5)

    def test_loo_scores_once_and_no_test_y_exposure(self):
        samples,pairs=source()
        models={'train_mean_baseline':(), 'branch_three':m.BRANCH_GENES}
        rows=m.loo_rows(samples,pairs,models,prediction)
        self.assertEqual(len(rows),12)
        self.assertEqual(len({r['id'] for r in rows}),12)
        self.assertTrue(all(r['n_train']==11 for r in rows))
        self.assertTrue(all(r['neighbor_counts']['1e-05']>=1 for r in rows))
        self.assertTrue(all(math.isfinite(r['gap']) for r in rows))

    def test_mutating_target_preserves_own_predictions_and_neighborhood(self):
        samples,pairs=source()
        models={'train_mean_baseline':(), 'branch_three':m.BRANCH_GENES}
        a=m.loo_rows(samples,pairs,models,prediction)
        changed=copy.deepcopy(samples)
        changed['S00']['responses']['anthocyanidins_log1p_total']=100000
        b=m.loo_rows(changed,pairs,models,prediction)
        self.assertEqual(a[0]['pred'], b[0]['pred'])
        self.assertEqual(a[0]['novelty'], b[0]['novelty'])
        self.assertEqual(a[0]['neighbor_counts'], b[0]['neighbor_counts'])

    def test_train_only_novelty(self):
        train={'A':{'genes':dict(zip(m.BRANCH_GENES,(0,0,0)))},
               'B':{'genes':dict(zip(m.BRANCH_GENES,(2,2,2)))}}
        x=dict(zip(m.BRANCH_GENES,(1,1,1)))
        self.assertAlmostEqual(m.expression_novelty(train,x),math.sqrt(3))

    def test_pairwise_identity_symmetry_validation(self):
        s,d=source()
        d['S00']['S01']=999
        with self.assertRaises(ValueError): m.loo_rows(s,d,{'train_mean_baseline':()},prediction)

    def test_mean_and_pair_metric_arithmetic(self):
        rows=[{'y':{'e':1},'pred':{'e':{'a':0,'b':1}}},
              {'y':{'e':3},'pred':{'e':{'a':1,'b':3}}}]
        z=m.risk(rows,'e','a')
        self.assertAlmostEqual(z['rmse'],math.sqrt(2.5),places=7)
        self.assertEqual(z['mae'],1.5)
        self.assertEqual(m.risk(rows,'e','b')['rmse'],0)

    def test_summary_non_estimable_constant_gap_and_no_fake_p_value(self):
        s,d=source()
        d={a:{b:0 if a==b else 1 for b in s} for a in s}
        models={'train_mean_baseline':(), 'branch_three':m.BRANCH_GENES}
        rows=m.loo_rows(s,d,models,prediction)
        out=m.summarize(rows,models)
        e=out['endpoints']['anthocyanidins_log1p_total']
        p=e['tolerance_sensitivity']['1e-05']
        self.assertIsNone(p['rho_gap_branch_squared_error'])
        self.assertIsNone(p['distance_pattern_support'])
        self.assertIsNone(p['tiers']['near'])
        self.assertNotIn('p_value',out)

if __name__=='__main__': unittest.main()
