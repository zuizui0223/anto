"""Independent checks of route-specific first-recovery risks (synthetic only)."""
import importlib.util
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
SPEC = importlib.util.find_spec("anto_afterlife_route_window_v12")
if SPEC is not None:
    import anto_afterlife_route_window_v12 as m
else:
    m = None


class ModuleExists(unittest.TestCase):
    def test_route_resolved_model_exists(self):
        self.assertIsNotNone(SPEC, "Route-specific recovery model has not been implemented")


@unittest.skipIf(m is None, "Implementation not yet present")
class RouteWindowTests(unittest.TestCase):
    def setUp(self):
        self.loss = m.Loss("0.02", "0.30", "0.01")
        self.rates = m.Rates("0.40", "0.20", self.loss)

    def test_exact_first_step_solution(self):
        r = m.exact_from_both(self.rates)
        a, b, c, x, y = F(2,5), F(1,5), F(1,50), F(3,10), F(1,100)
        z = a+b+c+x+y
        self.assertEqual(r["ancestral"], (a+y*a/(a+c+x))/z)
        self.assertEqual(r["alternative"], (b+x*b/(b+c+y))/z)
        # Independent linear system for the three transient resource states.
        q = m.generator(self.rates)
        solved = np.linalg.solve(-q[:3,:3], q[:3,3:5])
        np.testing.assert_allclose(solved[0], [float(r["ancestral"]),float(r["alternative"])], atol=1e-13)

    def test_generator_conserves_mass_and_has_no_return_from_logged_first_event(self):
        q = m.generator(self.rates)
        np.testing.assert_allclose(q.sum(axis=1), 0, atol=1e-14)
        self.assertTrue(np.all(q-np.diag(np.diag(q)) >= 0))
        self.assertTrue(np.all(q[3:,:] == 0))
        for t in (0, 1, 5, 20):
            p = m.aged_inventory(self.loss, t)
            self.assertAlmostEqual(p.sum(), 1)
            self.assertTrue(np.all(p >= 0))
            self.assertEqual(p[3]+p[4], 0)

    def test_finite_alternative_probability_matches_independent_quadrature(self):
        h=7.0
        a,b=float(self.rates.ancestral),float(self.rates.alternative)
        c,x,y=map(float,(self.loss.common,self.loss.ancestral,self.loss.alternative))
        lam=a+b+c+x+y; d=b+c+y
        # Direct B event, plus A loss at s followed by B event by h.
        independent=quad(lambda s: b*np.exp(-lam*s),0,h)[0]
        independent+=quad(lambda s: x*np.exp(-lam*s)*b/d*(-np.expm1(-d*(h-s))),0,h)[0]
        actual=m.after_delay(self.rates,self.loss,0,h)["alternative_first"]
        self.assertAlmostEqual(actual,independent,places=12)

    def test_aging_formula_matches_matrix_exponential(self):
        for t in (0, 1, 5, 15):
            analytic=m.after_delay(self.rates,self.loss,t)
            finite=m.after_delay(self.rates,self.loss,t,2000)
            self.assertAlmostEqual(analytic["alternative_first"],finite["alternative_first"],places=12)
            self.assertAlmostEqual(analytic["ancestral_first"],finite["ancestral_first"],places=12)
            self.assertAlmostEqual(analytic["any_first"],finite["any_first"],places=12)

    def test_interior_window_increases_absolute_alternative_but_reduces_total(self):
        peak=m.alternative_peak(self.rates,self.loss)
        self.assertEqual(peak["kind"],"interior")
        t=peak["age"]
        start=m.after_delay(self.rates,self.loss,0)
        top=m.after_delay(self.rates,self.loss,t)
        self.assertGreater(top["alternative_first"],start["alternative_first"])
        self.assertLess(top["any_first"],start["any_first"])
        self.assertLess(top["ancestral_first"],start["ancestral_first"])
        for dt in (-.001,.001):
            self.assertLess(m.after_delay(self.rates,self.loss,t+dt)["alternative_first"],top["alternative_first"])
        self.assertLess(m.after_delay(self.rates,self.loss,1000)["alternative_first"],1e-10)

    def test_more_ancestral_resource_loss_can_increase_alternative_first_not_total(self):
        last_b=-1.; last_total=2.
        for x in (.001,.01,.1,.3,1.,10.):
            r=m.exact_from_both(m.Rates(.4,.2,m.Loss(.02,x,.01)))
            self.assertGreater(float(r["alternative"]),last_b)
            self.assertLess(float(r["any"]),last_total)
            last_b=float(r["alternative"]);last_total=float(r["any"])

    def test_no_first_event_censoring_removes_the_increase(self):
        values=[m.alternative_ever_without_preemption(self.rates,self.loss,t) for t in (0,1,5,10,30)]
        self.assertTrue(all(a>b for a,b in zip(values,values[1:])))
        other=m.Rates(100,.2,m.Loss(.02,100,.01))
        self.assertEqual(m.alternative_ever_without_preemption(other,self.loss,5),values[2])

    def test_symmetric_routes_have_equal_probabilities_and_no_peak(self):
        loss=m.Loss(.02,.1,.1); rates=m.Rates(.2,.2,loss)
        self.assertEqual(m.alternative_peak(rates,loss)["kind"],"boundary_at_zero")
        previous=1.
        for t in (0,1,5,20):
            r=m.after_delay(rates,loss,t)
            self.assertAlmostEqual(r["ancestral_first"],r["alternative_first"])
            self.assertLess(r["any_first"],previous)
            previous=r["any_first"]

    def test_permanently_preserved_alternative_has_no_finite_optimal_age(self):
        loss=m.Loss(0,.3,0); rates=m.Rates(.4,.2,loss)
        self.assertEqual(m.alternative_peak(rates,loss)["kind"],"increasing_to_limit")
        self.assertAlmostEqual(m.after_delay(rates,loss,100)["any_first"],1.)

    def test_zero_rates_never_create_expression_or_divide_by_zero(self):
        zero=m.Loss(0,0,0); rates=m.Rates(0,0,zero)
        r=m.after_delay(rates,zero,10)
        self.assertEqual(r["any_first"],0.)
        self.assertIsNone(r["alternative_fraction_among_first_recoveries"])
        self.assertEqual(m.after_delay(rates,zero,10,0)["any_first"],0.)
        self.assertEqual(m.alternative_peak(rates,zero)["kind"],"constant_zero")

    def test_reject_nonfinite_negative_time_and_rate(self):
        for bad in (-.1,float("inf"),float("nan"),True):
            with self.assertRaises(ValueError):m.Loss(bad,0,0)
            with self.assertRaises(ValueError):m.Rates(bad,0,self.loss)
        for bad in (-1,float("nan"),float("inf")):
            with self.assertRaises(ValueError):m.after_delay(self.rates,self.loss,bad)
            with self.assertRaises(ValueError):m.after_delay(self.rates,self.loss,0,bad)

    def test_independent_sampled_lifetime_oracle(self):
        # Independent construction: resource lifetimes and two latent clocks,
        # not simulation of generator rows or the analytic recovery formula.
        rng=np.random.default_rng(20261020)
        n=100000; horizon=10.
        c,x,y=.02,.3,.01; a,b=.4,.2
        for age in (0.,5.,20.):
            c_live=rng.exponential(1/c,n)>age
            a_live=rng.exponential(1/x,n)>age
            b_live=rng.exponential(1/y,n)>age
            tc=rng.exponential(1/c,n)
            ta_loss=rng.exponential(1/x,n)
            tb_loss=rng.exponential(1/y,n)
            ta=rng.exponential(1/a,n)
            tb=rng.exponential(1/b,n)
            ta[~(c_live & a_live & (ta<tc) & (ta<ta_loss))]=np.inf
            tb[~(c_live & b_live & (tb<tc) & (tb<tb_loss))]=np.inf
            observed=((ta<tb)&(ta<=horizon)).mean(),((tb<ta)&(tb<=horizon)).mean()
            result=m.after_delay(self.rates,self.loss,age,horizon)
            expected=result["ancestral_first"],result["alternative_first"]
            for value,prob in zip(observed,expected):
                self.assertLess(abs(value-prob),6*np.sqrt(prob*(1-prob)/n)+1/n)

    def test_continuous_neighborhood_and_deterministic_grid(self):
        grid=m.parameter_grid()
        self.assertEqual(len(grid),1125)
        interior=[r for r in grid if r["peak_kind"]=="interior"]
        self.assertTrue(interior)
        for r in interior:
            self.assertGreater(r["alternative_at_peak"],r["alternative_at_zero"])
            self.assertLessEqual(r["any_at_peak"],r["any_at_zero"]+1e-12)
        for factor in (.9,1.,1.1):
            rates=m.Rates(.4*factor,.2,self.loss)
            self.assertEqual(m.alternative_peak(rates,self.loss)["kind"],"interior")

if __name__ == "__main__": unittest.main()
