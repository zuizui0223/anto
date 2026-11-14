"""Offline algebra, source-tree and nested holdout integrity checks."""
import sys, unittest
from pathlib import Path
from dataclasses import dataclass
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import analyze_phylo_expression_v32 as m
except ModuleNotFoundError:
    m=None
@dataclass
class Node:
    label: str=None
    children: tuple=()
    length: float=None

def panel(n=9):
    x=np.arange(n*3,dtype=float).reshape(n,3)/7
    x[:,1]=np.sin(x[:,1]); x[:,2]=np.cos(x[:,2])
    k=np.exp(-np.abs(np.arange(n)[:,None]-np.arange(n)[None,:])/3)
    y=np.column_stack((np.sin(np.arange(n))+x[:,0],np.cos(np.arange(n))))
    return k,y,x

class V32(unittest.TestCase):
    def setUp(self):self.assertIsNotNone(m,'v32 implementation absent: expected RED')
    def test_fixed_spec(self):
        self.assertEqual(m.ETAS,(.01,.1,1.,10.));self.assertEqual(m.GENE_LAMBDA,10.)
    def test_shared_path_covariance(self):
        t=Node(children=(Node(children=(Node('A',length=1),Node('B',length=2)),length=3),Node('C',length=4)),length=999)
        ids,k=m.shared_covariance(t)
        self.assertEqual(ids,['A','B','C']); np.testing.assert_allclose(k,[[4,3,0],[3,5,0],[0,0,4]])
        self.assertEqual(k[0,0]+k[1,1]-2*k[0,1],3)
    def test_bad_tree_rejected(self):
        for v in (None,-1,float('nan')):
            with self.assertRaises(ValueError):m.shared_covariance(Node(children=(Node('A',length=v),Node('B',length=1))))
        with self.assertRaises(ValueError):m.shared_covariance(Node(children=(Node('A',length=1),Node('A',length=2))))
    def test_phylo_zero_crosscov_mean_limit(self):
        out=m.predict(np.eye(5),np.zeros((2,5)),np.arange(5.),.1)
        np.testing.assert_allclose(out,[2,2])
    def test_constant_genes_equal_phylogeny_only(self):
        k,y,x=panel();tr=np.arange(8);te=[8]
        a=m.predict(k[:8,:8],k[te,:8],y[:8,0],.1)
        b=m.predict(k[:8,:8],k[te,:8],y[:8,0],.1,np.ones((8,3))*4,np.ones((1,3))*4)
        np.testing.assert_allclose(a,b,atol=1e-12)
    def test_constant_response(self):
        k,y,x=panel()
        out=m.predict(k[:8,:8],k[[8],:8],np.ones(8)*5,1.,x[:8],x[[8]])
        np.testing.assert_allclose(out,[5],atol=1e-10)
    def test_direct_gls_ridge_reference(self):
        k,y,x=panel();yt=y[:8,0];h=np.mean(np.diag(k[:8,:8]));c=k[:8,:8]/h+.1*np.eye(8)
        mu=x[:8].mean(0);sd=x[:8].std(0);a=np.column_stack((np.ones(8),(x[:8]-mu)/sd));b=np.column_stack((np.ones(1),(x[[8]]-mu)/sd))
        inv=np.linalg.inv(c);beta=np.linalg.solve(a.T@inv@a+np.diag([0,10,10,10]),a.T@inv@yt)
        expected=b@beta+k[[8],:8]/h@inv@(yt-a@beta)
        np.testing.assert_allclose(m.predict(k[:8,:8],k[[8],:8],yt,.1,x[:8],x[[8]]),expected,atol=1e-11)
    def test_train_feature_affine_invariance(self):
        k,y,x=panel();a=m.predict(k[:8,:8],k[[8],:8],y[:8,0],1.,x[:8],x[[8]])
        b=m.predict(k[:8,:8],k[[8],:8],y[:8,0],1.,x[:8]*3+100,x[[8]]*3+100)
        np.testing.assert_allclose(a,b,atol=1e-10)
    def test_nested_eta_matches_direct_inner_loo(self):
        k,y,x=panel(7);eta,score=m.select_eta(k,y[:,0])
        manual=[]
        for e in m.ETAS:
            errors=[]
            for j in range(7):
                tr=[i for i in range(7) if i!=j];kt=k[np.ix_(tr,tr)];cc=k[np.ix_([j],tr)];h=np.diag(kt).mean();iv=np.linalg.inv(kt/h+e*np.eye(6));v=y[tr,0];o=np.ones(6);mu=(o@iv@v)/(o@iv@o);pr=mu+(cc/h@iv@(v-mu))[0];errors.append((pr-y[j,0])**2)
            manual.append(np.mean(errors))
        self.assertEqual(eta,m.ETAS[int(np.argmin(manual))]);np.testing.assert_allclose(score,manual,rtol=1e-10)
    def test_heldout_outcome_mutation_cannot_change_prediction_or_selection(self):
        k,y,x=panel();a=m.outer_fold(k,y,x,list(range(7)),[7,8]);yy=y.copy();yy[7:]=1e8
        b=m.outer_fold(k,yy,x,list(range(7)),[7,8]);self.assertEqual(a['eta'],b['eta'])
        for key in a['pred']:np.testing.assert_array_equal(a['pred'][key],b['pred'][key])
    def test_train_test_overlap_rejected(self):
        k,y,x=panel()
        with self.assertRaises(ValueError):m.outer_fold(k,y,x,list(range(8)),[7,8])
        with self.assertRaises(ValueError):m.outer_fold(k,y,x,[0,0,1,2,3,4,5],[6,7,8])
    def test_invalid_covariance_or_shape(self):
        with self.assertRaises(ValueError):m.predict(np.zeros((4,4)),np.zeros((1,4)),np.arange(4.),1.)
        with self.assertRaises(ValueError):m.predict(np.eye(4),np.zeros((1,3)),np.arange(4.),1.)
        k,y,x=panel();k[0,1]=99
        with self.assertRaises(ValueError):m.predict(k,k[[0]],y[:,0],1.)
    def test_both_metrics_required(self):
        self.assertFalse(m.improves({'rmse':1.,'mae':2.},{'rmse':2.,'mae':1.}))
        self.assertTrue(m.improves({'rmse':1.,'mae':1.},{'rmse':2.,'mae':2.}))

if __name__=='__main__':unittest.main()
