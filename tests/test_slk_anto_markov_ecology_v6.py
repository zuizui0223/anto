"""Markov environmental persistence checks for the SLK x ANTO capacity bridge.

Tests exercise conditional finite extinction, pathwise growth, and distinct
mean-versus-quenched objects. One shared environment is experienced by ALL
descendants in each generation. It is not independent per individual.
"""
import itertools
import math
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from anto_slk_markov_ecology_v6 import (
    MarkovCapacity, symmetric_markov, panel
)


def independent_finite_oracle(process,horizon):
    """Separate rational enumeration over every Markov environment history."""
    s=F(0)
    total=[F(0),F(0)]
    p=process.persistence
    for seq in itertools.product(range(2),repeat=horizon):
        weight=F(1,2)
        for x,y in zip(seq,seq[1:]):
            weight *= p if x==y else 1-p
        z=(F(0),F(0))
        for env in reversed(seq):
            z=process.phases[env].pgf(z)
        s+=weight
        for i in range(2):
            total[i]+=weight*z[i]
    if s!=1:
        raise ArithmeticError("Exact environment path weights do not sum to one")
    return tuple(1-x for x in total)


class MarkovCapacityV6Tests(unittest.TestCase):
    def setUp(self):
        self.low_fast=symmetric_markov(F(1,10),F(1,10))
        self.low_slow=symmetric_markov(F(9,10),F(1,10))
        self.high_fast=symmetric_markov(F(1,10),F(9,10))
        self.high_slow=symmetric_markov(F(9,10),F(9,10))

    def test_stationary_environment_fraction_and_run_length(self):
        for p in (self.low_fast,self.low_slow):
            self.assertEqual(p.stationary_distribution(),(F(1,2),F(1,2)))
            self.assertEqual([sum(x) for x in p.transition_matrix()],[1,1])
            self.assertEqual(p.environment_autocorrelation(),2*p.persistence-1)
            self.assertEqual(p.mean_environment_run_length(),1/(1-p.persistence))
            self.assertEqual(p.stationary_next(),(F(1,2),F(1,2)))

    def test_first_reexpression_is_identical_across_environment_memory(self):
        for x,y in ((self.low_fast,self.low_slow),
                    (self.high_fast,self.high_slow)):
            self.assertEqual(x.first_expression_probability_from_C(),
                             y.first_expression_probability_from_C())
        self.assertEqual(self.low_fast.first_expression_probability_from_C(),
                         F(3153,31250))
        self.assertEqual(self.high_fast.first_expression_probability_from_C(),
                         F(21393,31250))

    def test_environment_averaged_reproduction_is_identical_for_each_pair(self):
        for fast,slow in ((self.low_fast,self.low_slow),
                          (self.high_fast,self.high_slow)):
            self.assertEqual(fast.arithmetic_mean_matrix(),
                             slow.arithmetic_mean_matrix())
            self.assertTrue(all(sum(row)==F(26,25)
                                for row in fast.arithmetic_mean_matrix()))

    def test_finite_shared_environment_exact_rational_oracle(self):
        for process in (self.low_fast,self.low_slow,self.high_fast,self.high_slow):
            for horizon in (1,2,3,4):
                expected=independent_finite_oracle(process,horizon)
                observed=process.finite_capacity_survival(horizon)
                for a,b in zip(observed["survival_by_initial_type"],expected):
                    self.assertAlmostEqual(a,float(b),places=13)
                self.assertAlmostEqual(observed["environment_mass"],1,places=13)
                self.assertEqual(observed["enumerated_histories"],2**horizon)

    def test_first_generation_capacity_survival_depends_not_on_memory_or_switching(self):
        for process in (self.low_fast,self.low_slow,self.high_fast,self.high_slow):
            r=process.finite_capacity_survival(1)
            self.assertAlmostEqual(r["survival_by_initial_type"][0],
                r["survival_by_initial_type"][1],places=13)
        self.assertAlmostEqual(
            self.low_fast.finite_capacity_survival(1)["survival_by_initial_type"][1],
            self.high_slow.finite_capacity_survival(1)["survival_by_initial_type"][1],
            places=13)

    def test_stochastic_environment_reverses_14_generation_capacity_survival(self):
        low_a=self.low_fast.finite_capacity_survival(14)["survival_by_initial_type"][1]
        low_b=self.low_slow.finite_capacity_survival(14)["survival_by_initial_type"][1]
        high_a=self.high_fast.finite_capacity_survival(14)["survival_by_initial_type"][1]
        high_b=self.high_slow.finite_capacity_survival(14)["survival_by_initial_type"][1]
        self.assertAlmostEqual(low_a,0.16169930199163,places=10)
        self.assertAlmostEqual(low_b,0.34504653033817,places=10)
        self.assertAlmostEqual(high_a,low_b,places=10)
        self.assertAlmostEqual(high_b,low_a,places=10)
        self.assertGreater(low_b,2*low_a)
        self.assertGreater(high_a,2*high_b)

    def test_stationary_environment_phenotype_half_switching_control(self):
        vals=[]
        for memory in (F(1,10),F(1,2),F(9,10)):
            process=symmetric_markov(memory,F(1,2))
            vals.append(process.finite_capacity_survival(12)["survival_by_initial_type"][1])
            r=process.sample_typical_growth(steps=8000,burn_in=100,seed=20261020)
            self.assertAlmostEqual(r["growth_per_generation"],1.04,places=12)
        self.assertLess(max(vals)-min(vals),1e-11)

    def test_pathwise_log_growth_reversal_but_not_exact_survival_certificate(self):
        # Margins deliberately wide enough not to depend on one sampled path.
        for sw,fast,slow in ((F(1,10),self.low_fast,self.low_slow),
                             (F(9,10),self.high_fast,self.high_slow)):
            a=fast.sample_typical_growth(steps=90000,burn_in=2000,seed=20261020)
            b=slow.sample_typical_growth(steps=90000,burn_in=2000,seed=20261020)
            if sw==F(1,10):
                self.assertLess(a["growth_per_generation"],0.99)
                self.assertGreater(b["growth_per_generation"],1.15)
            else:
                self.assertGreater(a["growth_per_generation"],1.15)
                self.assertLess(b["growth_per_generation"],0.99)
            for out,expected in ((a,float(fast.persistence)),
                                  (b,float(slow.persistence))):
                self.assertTrue(out["sampled_only_not_a_survival_proof"])
                self.assertAlmostEqual(out["environment_A_fraction"],.5,delta=.03)
                self.assertAlmostEqual(out["environment_stay_frequency"],
                                       expected,delta=.015)

    def test_common_loss_scales_every_sample_path_growth(self):
        for p in (F(1,10),F(9,10)):
            with_loss=symmetric_markov(p,F(1,10),F(1,5))
            without=symmetric_markov(p,F(1,10),0)
            x=with_loss.sample_typical_growth(steps=10000,burn_in=100,seed=8)
            y=without.sample_typical_growth(steps=10000,burn_in=100,seed=8)
            self.assertAlmostEqual(x["growth_per_generation"],
                                   F(4,5)*y["growth_per_generation"],places=12)
            self.assertAlmostEqual(
                with_loss.estimated_common_loss_boundary(x["growth_per_generation"]),
                1-1/y["growth_per_generation"],places=12)

    def test_finite_survival_does_not_prove_eventual_survival(self):
        r=self.low_fast.finite_capacity_survival(14)
        self.assertGreater(r["survival_by_initial_type"][1],0)
        self.assertTrue(r["finite_horizon_not_permanent_persistence"])
        self.assertFalse(r["empirical_plant_evolution_estimated"])

    def test_independent_environment_chain_paths_are_reproducible(self):
        p=self.low_fast
        x=p.sample_typical_growth(steps=5000,burn_in=100,seed=1234)
        y=p.sample_typical_growth(steps=5000,burn_in=100,seed=1234)
        self.assertEqual(x,y)
        z=p.sample_typical_growth(steps=5000,burn_in=100,seed=1235)
        self.assertNotEqual(x["empirical_log_growth"],z["empirical_log_growth"])

    def test_invalid_inputs_fail_closed(self):
        from anto_slk_recurrent_capacity_v4 import CapacityProcess
        for bad in (0,1,True,-1,F(3,2)):
            with self.assertRaises(ValueError):
                symmetric_markov(bad,F(1,10))
        for bad in (0,1,True,-1):
            with self.assertRaises(ValueError):
                symmetric_markov(F(1,10),bad)
        for h in (0,True,-1,17):
            with self.assertRaises(ValueError):
                self.low_fast.finite_capacity_survival(h)
        for args in ((0,100),(100,0),(100,100)):
            with self.assertRaises(ValueError):
                self.low_fast.sample_typical_growth(steps=args[0],burn_in=args[1])
        with self.assertRaises(TypeError):
            MarkovCapacity((self.low_fast.phases[0],"invalid"),F(1,2))
        with self.assertRaises(ValueError):
            MarkovCapacity(self.low_fast.phases,F(1,2),initial_environment=(1,0))

    def test_replay_reproduces_environment_memory_contrast_and_nonclaims(self):
        r=panel()
        self.assertEqual(r["schema"],"anto_slk_markov_environment_v6")
        self.assertEqual(r["environment_process"],"shared stationary two-state Markov")
        self.assertFalse(r["observed_flower_trait_regain"])
        self.assertFalse(r["field_measured_environment_process"])
        self.assertFalse(r["formal_lyapunov_sign_certificate"])
        self.assertFalse(r["slk_selection_fitness_calibrated"])
        self.assertEqual(r["Cirsium_empirical_owner"],"zuizui0223/aza3")
        index={(v["trait_switch"],v["environment_persistence"]):v for v in r["scenarios"]}
        self.assertEqual(len(index),6)
        self.assertTrue(all(v["finite_horizon"]==14 for v in index.values()))
        self.assertTrue(all(v["finite_genealogy_survival_from_C"]>0 for v in index.values()))
        self.assertGreater(index["1/10","9/10"]["finite_genealogy_survival_from_C"],
                           2*index["1/10","1/10"]["finite_genealogy_survival_from_C"])
        self.assertGreater(index["9/10","1/10"]["finite_genealogy_survival_from_C"],
                           2*index["9/10","9/10"]["finite_genealogy_survival_from_C"])


if __name__=="__main__":
    unittest.main()
