"""Matched-target/count controls; synthetic algebra fixtures only."""
import copy, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import analyze_target_size_controls_v33 as m
except ModuleNotFoundError:
    m=None

class MatchedTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m,'v33 not implemented')
    def test_frozen_constants(self):
        self.assertEqual(m.CAPS,(8,15)); self.assertEqual(m.REPS,31)
        self.assertEqual(m.SEED,20261114)
    def test_training_exact_size_and_target_absent(self):
        ids=[f'S{i:02}' for i in range(20)]
        for target in ids:
            for r in range(m.REPS):
                train=m.control_train(ids,target,4,8,r)
                self.assertEqual(len(train),16)
                self.assertEqual(len(set(train)),16)
                self.assertNotIn(target,train)
    def test_deterministic_order_independent(self):
        ids=[f'S{i:02}' for i in range(20)]
        self.assertEqual(m.control_train(ids,'S00',5,8,3),m.control_train(ids[::-1],'S00',5,8,3))
        self.assertNotEqual(m.control_train(ids,'S00',5,8,3),m.control_train(ids,'S00',5,8,4))
    def test_singleton_control_is_same_training_set(self):
        ids=[f'S{i:02}' for i in range(20)]
        for r in range(m.REPS):
            self.assertEqual(m.control_train(ids,'S00',1,8,r),ids[1:])
    def test_invalid_design_rejected(self):
        ids=list('ABCDEFGHIJ')
        for target,size,cap,rep in [('Z',2,8,0),('A',0,8,0),('A',7,8,0),('A',2,2,0),('A',2,8,31)]:
            with self.assertRaises(ValueError):m.control_train(ids,target,size,cap,rep)
        with self.assertRaises(ValueError):m.control_train(ids+['A'],'A',2,8,0)
    def test_risk_pools_squared_errors_not_rmse(self):
        a=m.risk_from_losses(np.array([1.,9.]),np.array([1.,3.]))
        self.assertAlmostEqual(a['rmse'],np.sqrt(5),7)
        self.assertEqual(a['mae'],2.)
    def test_pair_decomposition(self):
        real={'phylo':np.array([4.,9.,1.]),'phylo_genes':np.array([1.,4.,1.])}
        controls={'phylo':np.array([[1.,4.,1.],[4.,4.,1.]]),
                  'phylo_genes':np.array([[1.,1.,1.],[1.,1.,1.]])}
        x=m.paired_effect(real,controls,[[0,1],[2]])
        self.assertAlmostEqual(x['mean_increment_difference'],(1.5+2+0)/3,7)
        self.assertAlmostEqual(x['mean_increment_difference'],x['phylo_loss_penalty']-x['phylo_genes_loss_penalty'],7)
        self.assertEqual(x['n_positive_original_clades'],1)
    def test_all_control_draws_retained(self):
        real={'phylo':np.array([4.,4.]),'phylo_genes':np.array([1.,1.])}
        c={'phylo':np.array([[4.,4.],[100.,100.]]),'phylo_genes':np.array([[1.,1.],[1.,1.]])}
        z=m.paired_effect(real,c,[[0],[1]])
        self.assertLess(z['mean_increment_difference'],0)
        self.assertEqual(z['n_control_repetitions'],2)
    def test_bad_loss_geometry_rejected(self):
        real={'phylo':np.ones(3),'phylo_genes':np.ones(3)}
        c={'phylo':np.ones((2,2)),'phylo_genes':np.ones((2,2))}
        with self.assertRaises(ValueError):m.paired_effect(real,c,[[0,1],[2]])
        c={key:np.ones((2,3)) for key in real}
        with self.assertRaises(ValueError):m.paired_effect(real,c,[[0,1],[1,2]])
    def test_constant_response_prediction(self):
        k=np.eye(9)+.2; y=np.ones((9,2))*3; g=np.arange(27).reshape(9,3)/20
        out=m.predict_subset(k,y,g,list(range(7)),[7])
        for val in out['pred'].values():np.testing.assert_allclose(val,3,atol=1e-12)
    def test_test_y_and_unused_y_do_not_enter_fit(self):
        rng=np.random.default_rng(13); a=rng.normal(size=(9,3)); k=a@a.T+np.eye(9)
        y=rng.normal(size=(9,2)); g=rng.normal(size=(9,3)); tr=list(range(7))
        first=m.predict_subset(k,y,g,tr,[7])
        changed=y.copy(); changed[7:]=99999
        second=m.predict_subset(k,changed,g,tr,[7])
        self.assertEqual(first['eta'],second['eta'])
        for name in first['pred']:np.testing.assert_array_equal(first['pred'][name],second['pred'][name])
    def test_predict_subset_matches_parent_whole_fold(self):
        import analyze_phylo_expression_v32 as old
        rng=np.random.default_rng(17); a=rng.normal(size=(10,3)); k=a@a.T+np.eye(10)
        y=rng.normal(size=(10,2)); g=rng.normal(size=(10,3)); tr=list(range(8));te=[8,9]
        ref=old.outer_fold(k,y,g,tr,te); new=m.predict_subset(k,y,g,tr,te)
        self.assertEqual(ref['eta'],new['eta'])
        for name in m.MODELS:np.testing.assert_allclose(ref['pred'][name],new['pred'][name],rtol=1e-12,atol=1e-12)
    def test_overlap_and_duplicate_index_rejected(self):
        k=np.eye(8); y=np.ones((8,2));g=np.ones((8,3))
        for tr,te in [([0,1,2,3],[3]),([0,1,1,2,3],[4]),([0,1,2,3],[8])]:
            with self.assertRaises(ValueError):m.predict_subset(k,y,g,tr,te)
    def test_null_matched_effect_exactly_zero(self):
        real={'phylo':np.array([1.,2.]),'phylo_genes':np.array([.5,1.])}
        c={key:np.tile(x,(31,1)) for key,x in real.items()}
        x=m.paired_effect(real,c,[[0],[1]])
        self.assertEqual(x['mean_increment_difference'],0)
    def test_primary_gate_preserves_any_negative(self):
        good={'mean_increment_difference':1,'n_positive_original_clades':4,'n_original_clades':6}
        self.assertTrue(m.directional_gate([good]*4))
        bad=dict(good,mean_increment_difference=-1)
        self.assertFalse(m.directional_gate([good,good,bad,good]))
        bad=dict(good,n_positive_original_clades=3)
        self.assertFalse(m.directional_gate([good,good,bad,good]))

if __name__=='__main__':unittest.main()
