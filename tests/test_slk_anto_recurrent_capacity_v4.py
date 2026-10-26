"""Independent finite-population enumeration and exact boundary checks."""
import importlib.util
import itertools
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
if importlib.util.find_spec('anto_slk_recurrent_capacity_v4'):
    import anto_slk_recurrent_capacity_v4 as model
else:
    model = None


def full_oracle(p, generations, start):
    """Enumerate all outcomes of each of 2*(D+C) reproductive slots."""
    states = {(1, 0) if start == 0 else (0, 1): F(1)}
    for _ in range(generations):
        nxt = {}
        for (d, c), weight in states.items():
            individuals = [0] * d + [1] * c
            conditional = {(0, 0): F(1)}
            for parent in individuals:
                q = p.q_expressed if parent == 0 else p.q_silent
                mix = (1 - p.silencing, p.silencing) if parent == 0 else (p.activation, 1 - p.activation)
                slots = ((0, 0, 1-q*(1-p.loss)),
                         (1, 0, q*(1-p.loss)*mix[0]),
                         (0, 1, q*(1-p.loss)*mix[1]))
                for _ in range(2):
                    new = {}
                    for (u,v), mass in conditional.items():
                        for ud,uc,prob in slots:
                            key = (u+ud, v+uc)
                            new[key] = new.get(key,F(0)) + mass*prob
                    conditional = new
            for key, prob in conditional.items():
                nxt[key] = nxt.get(key,F(0)) + weight*prob
        assert sum(nxt.values()) == 1
        states = nxt
    return states.get((0,0), F(0))


class RecurrentCapacityTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(model, 'New recurrent-capacity implementation is missing')
        self.p = model.CapacityProcess('13/20','9/10','3/10','1/4','3/10')

    def test_mean_matrix_is_exact_and_coupling_is_essential(self):
        m = self.p.mean_matrix()
        self.assertEqual(m, ((F(273,400),F(91,400)),(F(189,500),F(441,500))))
        self.assertTrue(m[0][0] < 1 and m[1][1] < 1)
        self.assertEqual(self.p.regime(), 'supercritical')
        self.assertAlmostEqual(self.p.spectral_radius(),1.092,places=14)

    def test_generating_function_has_mass_one(self):
        self.assertEqual(self.p.pgf((F(1),F(1))), (F(1),F(1)))
        for row in self.p.mean_matrix():
            self.assertTrue(sum(row) < 2)
        self.assertEqual(self.p.pgf((F(0),F(0))),
                         ((1-self.p.q_expressed*(1-self.p.loss))**2,
                          (1-self.p.q_silent*(1-self.p.loss))**2))

    def test_independent_joint_count_oracle(self):
        for t in range(4):
            result = model.finite_extinction(self.p,t)
            for start in (0,1):
                self.assertEqual(result[start],full_oracle(self.p,t,start))

    def test_equal_demography_recovers_v3_scalar_limit(self):
        p = model.CapacityProcess('13/20','13/20',0,'1/4','3/10')
        r = model.survival(p)
        self.assertTrue(r['converged'])
        for lo,hi in zip(r['survival_lower'],r['survival_upper']):
            self.assertLessEqual(lo, float(F(120,169))+1e-13)
            self.assertGreaterEqual(hi, float(F(120,169))-1e-13)

    def test_ongoing_loss_reverses_perfect_inheritance_conclusion(self):
        p = model.CapacityProcess('13/20','13/20','1/4','1/4','3/10')
        self.assertEqual(p.regime(),'subcritical')
        self.assertEqual(model.survival(p)['survival_upper'],(0.0,0.0))

    def test_more_first_expression_but_no_long_run_capacity(self):
        p = model.CapacityProcess('13/20','9/10','3/10','1/4','9/10')
        self.assertEqual(self.p.first_expression_from_silent(),F(342279,1000000))
        self.assertEqual(p.first_expression_from_silent(),F(812511,1000000))
        self.assertGreater(p.first_expression_from_silent(),self.p.first_expression_from_silent())
        self.assertEqual(p.regime(),'subcritical')
        self.assertEqual(model.survival(p)['survival_upper'],(0.0,0.0))
        good=model.survival(self.p)
        self.assertGreater(good['survival_lower'][1],0.35)
        self.assertLess(good['survival_upper'][1],0.36)

    def test_exact_critical_activation_and_nearby_signs(self):
        ac=model.activation_boundary(self.p)
        self.assertEqual(ac,F(1651,2268))
        args=(self.p.q_expressed,self.p.q_silent,self.p.loss,self.p.silencing)
        self.assertEqual(model.CapacityProcess(*args,ac).regime(),'critical')
        self.assertEqual(model.survival(model.CapacityProcess(*args,ac))['survival_upper'],(0.0,0.0))
        for a,regime in ((ac-F(1,10**8),'supercritical'),(ac+F(1,10**8),'subcritical')):
            self.assertEqual(model.CapacityProcess(*args,a).regime(),regime)

    def test_boundary_only_returned_in_declared_regime(self):
        p=model.CapacityProcess('13/20','13/20',0,'1/4','3/10')
        with self.assertRaises(ValueError):model.activation_boundary(p)

    def test_loss_boundary_and_initial_state_are_different_objects(self):
        dc=self.p.loss_boundary()
        low=model.CapacityProcess('13/20','9/10',dc-0.001,'1/4','3/10')
        high=model.CapacityProcess('13/20','9/10',dc+0.001,'1/4','3/10')
        self.assertEqual(low.regime(),'supercritical')
        self.assertEqual(high.regime(),'subcritical')
        self.assertEqual(model.initial_capacity_survival(self.p,0)['upper'],0)
        single=model.initial_capacity_survival(self.p,F(1,8))
        doubled=model.initial_capacity_survival(self.p,F(1,4))
        self.assertAlmostEqual(doubled['lower'],2*single['lower'],places=13)

    def test_finite_nonextinction_is_an_upper_bound_not_permanence(self):
        p=model.CapacityProcess('13/20','9/10','3/10','1/4','9/10')
        self.assertGreater(1-model.finite_extinction(p,3)[1],0)
        self.assertEqual(model.survival(p)['survival_upper'][1],0)

    def test_solver_brackets_narrow_and_residual_small(self):
        r=model.survival(self.p,tolerance=1e-12)
        self.assertTrue(r['converged'])
        self.assertLessEqual(r['width'],1e-12)
        self.assertLess(r['residual'],1e-12)
        self.assertLess(r['iterations'],2000)
        for lo,hi in zip(r['survival_lower'],r['survival_upper']):
            self.assertTrue(0 <= lo <= hi <= 1)
        limited=model.survival(self.p,max_iterations=1)
        self.assertFalse(limited['converged'])
        self.assertGreater(limited['width'],1e-3)

    def test_near_critical_brackets_do_not_stall_at_float_precision(self):
        p=model.CapacityProcess('13/20','9/10','2/5','1/10','1/10')
        r=model.survival(p,tolerance=1e-12)
        self.assertTrue(r['converged'])
        expected=(0.00027405774194360315,0.0010472818792832478)
        for value,lo,hi in zip(expected,r['survival_lower'],r['survival_upper']):
            self.assertLessEqual(lo,value+1e-15)
            self.assertGreaterEqual(hi,value-1e-15)
        self.assertLessEqual(r['width'],1e-12)

    def test_invalid_inputs_rejected(self):
        base=['13/20','9/10','3/10','1/4','3/10']
        for index,bad in ((0,0),(1,1),(2,1),(3,0),(4,0),(0,True),(4,'nan')):
            x=base.copy();x[index]=bad
            with self.assertRaises(ValueError):model.CapacityProcess(*x)
        for t in (-1,True,11):
            with self.assertRaises(ValueError):model.finite_extinction(self.p,t)
        with self.assertRaises(ValueError):model.survival(self.p,max_iterations=0)
        with self.assertRaises(ValueError):model.survival(self.p,tolerance=0)

    def test_report_keeps_biology_and_process_limits(self):
        out=model.panel()
        self.assertFalse(out['natural_regain_demonstrated'])
        self.assertFalse(out['organism_extinction_inferred'])
        self.assertFalse(out['fixation_inferred'])
        self.assertFalse(out['slk_payoff_calibrated'])
        self.assertEqual(out['empirical_owner'],'zuizui0223/aza3')
        self.assertEqual(len(out['activation_sweep']),19)

if __name__ == '__main__':unittest.main()
