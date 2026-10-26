"""Tests of environmental-run persistence using the unchanged v4 offspring model."""
import importlib
import importlib.util
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from anto_slk_recurrent_capacity_v4 import CapacityProcess, survival


def joint_count_oracle(phases, initial):
    """Forward enumeration of EACH reproductive slot, not a PGF composition."""
    states = {(1, 0) if initial == 0 else (0, 1): F(1)}
    for phase in phases:
        matrix = phase.mean_matrix()
        following = {}
        for (nd, nc), mass in states.items():
            offspring = {(0, 0): F(1)}
            for parent in [0]*nd + [1]*nc:
                pd, pc = (x/2 for x in matrix[parent])
                for _ in range(2):
                    nxt = {}
                    for (xd, xc), prob in offspring.items():
                        for key, weight in (((xd, xc), 1-pd-pc),
                                            ((xd+1, xc), pd), ((xd, xc+1), pc)):
                            nxt[key] = nxt.get(key, F(0)) + prob*weight
                    offspring = nxt
            for key, prob in offspring.items():
                following[key] = following.get(key, F(0)) + mass*prob
        states = following
    assert sum(states.values()) == 1
    return states


class PeriodicCapacityTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('anto_slk_periodic_capacity_v5'),
                             'Periodic-capacity implementation is not yet present')
        self.m = importlib.import_module('anto_slk_periodic_capacity_v5')
        self.A = CapacityProcess('9/10','2/5','1/5','1/10','1/10')
        self.B = CapacityProcess('2/5','9/10','1/5','1/10','1/10')

    def test_same_composition_different_persistence(self):
        fast = self.m.PeriodicCapacity((self.A,self.B)*4)
        block = self.m.PeriodicCapacity((self.A,)*4+(self.B,)*4)
        self.assertEqual(fast.average_mean(), block.average_mean())
        self.assertEqual(fast.regime(), 'subcritical')
        self.assertEqual(block.regime(), 'supercritical')
        self.assertAlmostEqual(fast.growth_per_generation(), .9689254303117699, 13)
        self.assertAlmostEqual(block.growth_per_generation(), 1.0164000507618296, 13)
        self.assertEqual(fast.phases[0].first_expression_from_silent(), F(984,15625))
        self.assertEqual(fast.phases[0], block.phases[0])

    def test_chronological_matrix_product_and_pgf_are_consistent(self):
        p = self.m.PeriodicCapacity((self.A,self.B))
        a,b = self.A.mean_matrix(), self.B.mean_matrix()
        expected = tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))
                               for j in range(2)) for i in range(2))
        self.assertEqual(p.mean_cycle(), expected)
        z = (F(1,3),F(2,5))
        self.assertEqual(p.exact_pgf(z), self.A.pgf(self.B.pgf(z)))
        self.assertNotEqual(p.exact_pgf(z), self.B.pgf(self.A.pgf(z)))
        self.assertEqual(p.exact_pgf((F(1),F(1))), (F(1),F(1)))

    def test_independent_joint_count_oracle(self):
        p = self.m.PeriodicCapacity((self.A,self.B))
        for initial in (0,1):
            dist = joint_count_oracle(p.phases, initial)
            self.assertEqual(p.exact_pgf((F(0),F(0)))[initial], dist.get((0,0),0))
            for child in (0,1):
                mean = sum(key[child]*prob for key,prob in dist.items())
                self.assertEqual(p.mean_cycle()[initial][child], mean)
            z=(F(1,3),F(2,5))
            value=sum(prob*z[0]**nd*z[1]**nc for (nd,nc),prob in dist.items())
            self.assertEqual(p.exact_pgf(z)[initial], value)

    def test_constant_environment_is_v4_not_a_new_model(self):
        p = CapacityProcess('13/20','9/10','3/10','1/4','3/10')
        ref = survival(p)
        for length in (1,2,4):
            period=self.m.PeriodicCapacity((p,)*length)
            result=period.survival()
            self.assertTrue(result['converged'])
            self.assertAlmostEqual(period.growth_per_generation(),p.spectral_radius(),13)
            for x,y in zip(result['survival_lower'],ref['survival_lower']):
                self.assertAlmostEqual(x,y,10)

    def test_cycle_rotation_preserves_boundary_not_founder_phase(self):
        phases=(self.A,)*4+(self.B,)*4
        first=self.m.PeriodicCapacity(phases)
        rotated=self.m.PeriodicCapacity(phases[1:]+phases[:1])
        self.assertAlmostEqual(first.growth_per_generation(),rotated.growth_per_generation(),13)
        x,y=first.survival(),rotated.survival()
        self.assertTrue(x['converged'] and y['converged'])
        qy=tuple(F(str(1-z)) for z in y['survival_lower'])
        mapped=self.A.pgf(qy)
        for a,b in zip(x['survival_lower'],mapped):
            self.assertAlmostEqual(a,1-float(b),10)
        self.assertNotAlmostEqual(x['survival_lower'][1],y['survival_lower'][1],6)

    def test_repeated_cycle_preserves_survival(self):
        phases=(self.A,)*4+(self.B,)*4
        once=self.m.PeriodicCapacity(phases)
        twice=self.m.PeriodicCapacity(phases*2)
        self.assertAlmostEqual(once.growth_per_generation(),twice.growth_per_generation(),13)
        for x,y in zip(once.survival()['survival_lower'],twice.survival()['survival_lower']):
            self.assertAlmostEqual(x,y,10)

    def test_fast_phenotype_switching_reverses_environmental_ranking(self):
        a=CapacityProcess('9/10','2/5','1/5','9/10','9/10')
        b=CapacityProcess('2/5','9/10','1/5','9/10','9/10')
        fast=self.m.PeriodicCapacity((a,b)*4)
        block=self.m.PeriodicCapacity((a,)*4+(b,)*4)
        self.assertEqual(fast.regime(),'supercritical')
        self.assertEqual(block.regime(),'subcritical')
        self.assertAlmostEqual(fast.growth_per_generation(),1.3085801927102407,13)
        self.assertAlmostEqual(block.growth_per_generation(),.9778264981191912,12)

    def test_common_retention_scales_growth_not_initial_supply(self):
        for k in (1,4):
            p=self.m.symmetric_period(k, F(1,10), F(1,5))
            no_loss=self.m.symmetric_period(k,F(1,10),0)
            self.assertAlmostEqual(p.growth_per_generation(),.8*no_loss.growth_per_generation(),13)
            self.assertAlmostEqual(p.common_loss_boundary(),1-1/no_loss.growth_per_generation(),13)
        threshold=self.m.symmetric_period(4,F(1,10),F(1,5)).common_loss_boundary()
        self.assertGreater(threshold,.21)
        self.assertLess(threshold,.22)

    def test_same_state_only_cannot_persist_in_the_witness(self):
        # Pure D and pure C, with switching DISABLED, each grow geometrically
        # at 2*(1-loss)*sqrt(q_hi*q_lo)=.96 over equal time in A and B.
        pure_growth=2*.8*(.9*.4)**.5
        self.assertAlmostEqual(pure_growth,.96,14)
        self.assertGreater(self.m.symmetric_period(4,F(1,10),F(1,5)).growth_per_generation(),1)

    def test_no_memory_half_switching_removes_run_length_effect(self):
        for k in (1,2,4,8):
            p=self.m.symmetric_period(k,F(1,2),F(1,5))
            self.assertAlmostEqual(p.growth_per_generation(),1.04,13)

    def test_exact_critical_case_and_zero_survival(self):
        p=CapacityProcess('5/8','5/8','1/5','1/5','2/5')
        period=self.m.PeriodicCapacity((p,)*3)
        self.assertEqual(period.regime(),'critical')
        self.assertEqual(period.survival()['survival_upper'],(0.,0.))
        fast=self.m.symmetric_period(1,F(1,10),F(1,5))
        self.assertEqual(fast.survival()['survival_upper'],(0.,0.))

    def test_convergence_and_numerical_brackets(self):
        p=self.m.symmetric_period(4,F(1,10),F(1,5))
        r=p.survival()
        self.assertTrue(r['converged'])
        self.assertLess(r['width'],1e-11)
        self.assertLess(r['residual'],1e-11)
        self.assertGreater(r['survival_lower'][1],0)
        q=p.survival(max_iterations=1)
        self.assertFalse(q['converged'])
        self.assertTrue(r['numerical_brackets_not_interval_certificates'])

    def test_reject_unsupported_limits_and_bad_inputs(self):
        for bad in ((),('A',), (self.A,)*65):
            with self.assertRaises((ValueError,TypeError)):
                self.m.PeriodicCapacity(bad)
        boundary=CapacityProcess('9/10','2/5','1/5',1,1)
        with self.assertRaises(ValueError): self.m.PeriodicCapacity((boundary,))
        for k in (0,True,-1):
            with self.assertRaises(ValueError): self.m.symmetric_period(k,F(1,10),F(1,5))
        p=self.m.symmetric_period(1,F(1,10),F(1,5))
        for tolerance in (0,True,float('nan')):
            with self.assertRaises(ValueError): p.survival(tolerance=tolerance)
        with self.assertRaises(ValueError): p.survival(max_iterations=0)
        with self.assertRaises(ValueError): p.exact_pgf((F(-1),F(0)))
        mixed=self.m.PeriodicCapacity((self.A,CapacityProcess('.4','.9','.1','.1','.1')))
        with self.assertRaises(ValueError): mixed.common_loss_boundary()

    def test_tiny_per_generation_growth_does_not_underflow_cycle_product(self):
        loss=1-F(1,10**20)
        p=self.m.symmetric_period(32,F(1,10),loss)
        result=p.growth_per_generation()
        self.assertGreater(result,0)
        self.assertLess(result,2e-20)
        self.assertEqual(p.regime(),'subcritical')

    def test_report_claims_are_conditional(self):
        report=self.m.panel()
        self.assertFalse(report['natural_regain_demonstrated'])
        self.assertFalse(report['organism_extinction_inferred'])
        self.assertFalse(report['fixation_or_optimal_strategy_inferred'])
        self.assertFalse(report['slk_payoff_calibrated'])
        self.assertEqual(report['empirical_owner'],'zuizui0223/aza3')
        self.assertEqual(report['environment_process'],'deterministic periodic, shared across descendants')


if __name__=='__main__': unittest.main()
