"""ANTO × SLK v6: a SHARED stochastic Markov environment.

This does not create a third floral phenotype or modify the v4 reproductive
kernel. D is visibly expressing, C is a reproducing nonexpressing but fully
competent state. Toolkit-lost offspring leave the modeled COMPETENT genealogy,
not necessarily the population of organisms.

A stationary two-state Markov chain drives ALL descendants at the same time.
The first environmental state is sampled from its exact (1/2,1/2) stationary
law. Environmental persistence alters autocorrelation without changing the
stationary time share of either environment. Do not confuse finite genealogical
presence with long-run survival or pathwise expected growth with fixation.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import exp, fsum, log
import random

from anto_slk_recurrent_capacity_v4 import CapacityProcess, probability

HALF = F(1, 2)


@dataclass(frozen=True)
class MarkovCapacity:
    phases: tuple
    persistence: F
    initial_environment: tuple = (HALF, HALF)

    def __post_init__(self):
        phases = tuple(self.phases)
        if len(phases) != 2 or any(not isinstance(p, CapacityProcess) for p in phases):
            raise TypeError("Require two v4 CapacityProcess regimes")
        p = probability(self.persistence)
        if not 0 < p < 1:
            raise ValueError("Restrict Markov environmental persistence to (0,1)")
        z = tuple(probability(x) for x in self.initial_environment)
        if z != (HALF, HALF):
            raise ValueError("Require stationary (1/2,1/2) initialization")
        if any(not 0 < c.activation < 1 or not 0 < c.silencing < 1
               for c in phases):
            raise ValueError("Require strictly positive capacity mean matrices")
        object.__setattr__(self, "phases", phases)
        object.__setattr__(self, "persistence", p)
        object.__setattr__(self, "initial_environment", z)

    def transition_matrix(self):
        p = self.persistence
        return ((p, 1-p), (1-p, p))

    def stationary_distribution(self):
        return self.initial_environment

    def stationary_next(self):
        P = self.transition_matrix()
        return tuple(sum(self.initial_environment[i]*P[i][j]
                         for i in range(2)) for j in range(2))

    def environment_autocorrelation(self):
        """For +/-1 coding, lag-one correlation is exactly 2*p_stay-1."""
        return 2*self.persistence-1

    def mean_environment_run_length(self):
        """Mean geometric run duration; not a measured ecological turnover."""
        return 1/(1-self.persistence)

    def arithmetic_mean_matrix(self):
        A, B = (x.mean_matrix() for x in self.phases)
        return tuple(tuple((A[i][j]+B[i][j])/2 for j in range(2))
                     for i in range(2))

    def first_expression_probability_from_C(self):
        """One C founder, stationary starting environment, one generation."""
        return sum(self.initial_environment[i]*
                   self.phases[i].first_expression_from_silent()
                   for i in range(2))

    def finite_capacity_survival(self, horizon):
        """All Markov histories; conditional offspring extinction via v4 PGFs.

        Finite-horizon probability that AT LEAST ONE toolkit-competent
        descendant is present at generation H, starting from one D or C.
        A single environmental history is SHARED by every offspring that
        generation. Averaging two PGFs before propagation is incorrect.

        This exhausts 2**H environment histories with stationarity exactly
        specified; only the PGF arithmetic is float, not rational.
        For independent validation tests, a separate short-horizon Fraction
        enumeration is used. Limit H to prevent exponential CPU/memory work.
        """
        if type(horizon) is not int or not 1 <= horizon <= 16:
            raise ValueError("H must be a positive integer no greater than 16")
        stay = float(self.persistence)
        matrices = [tuple(tuple(float(v)/2 for v in row)
                          for row in p.mean_matrix()) for p in self.phases]
        prob_terms = []
        extinction_d = []
        extinction_c = []
        for envs in product((0, 1), repeat=horizon):
            weight = .5
            for i in range(1,horizon):
                weight *= stay if envs[i] == envs[i-1] else 1-stay
            x=y=0.0
            for e in reversed(envs):
                (ad, ac), (cd, cc) = matrices[e]
                xx = (1-ad-ac+ad*x+ac*y)**2
                yy = (1-cd-cc+cd*x+cc*y)**2
                x,y=xx,yy
            prob_terms.append(weight)
            extinction_d.append(weight*x)
            extinction_c.append(weight*y)
        total=fsum(prob_terms)
        if abs(total-1)>1e-10:
            raise ArithmeticError("Environmental path weights lost mass")
        results = (1-fsum(extinction_d),1-fsum(extinction_c))
        if any(v < -1e-10 or v > 1+1e-10 for v in results):
            raise ArithmeticError("Invalid finite genealogical survival")
        return {
            "generations": horizon,
            "enumerated_histories": 2**horizon,
            "environment_mass": total,
            "survival_by_initial_type": results,
            "finite_horizon_not_permanent_persistence": True,
            "empirical_plant_evolution_estimated": False,
        }

    def sample_typical_growth(self, steps=100000, burn_in=2000, seed=20261020):
        """One *pathwise* normalized mean-matrix product, sampled.

        State composition starts at (1/2,1/2). Repeatedly multiply row
        v by the current environment's v4 reproduction mean matrix, normalize
        its sum, and accumulate the log normalization. This estimates the
        maximal (quenched) Lyapunov growth, not E[N(t)] or an extinction
        probability. Independent seeds gauge Monte Carlo variation only;
        they are not a formal error certificate for the Lyapunov exponent.
        """
        if (type(steps) is not int or type(burn_in) is not int
                or steps <= 0 or not 1 <= burn_in < steps):
            raise ValueError("Require integer steps>burn_in>=1")
        if type(seed) is not int:
            raise ValueError("Seed must be integer")
        stay=float(self.persistence)
        matrices = [tuple(tuple(float(z) for z in row)
                          for row in p.mean_matrix()) for p in self.phases]
        rng=random.Random(seed)
        e=rng.randrange(2)
        x=y=.5
        log_sum=0.0
        times_A=0
        same=0
        for t in range(steps):
            (a,b),(c,d)=matrices[e]
            xx=x*a+y*c
            yy=x*b+y*d
            norm=xx+yy
            if norm <= 0:
                raise ArithmeticError("Nonpositive normalization")
            x,y=xx/norm,yy/norm
            if t>=burn_in:
                log_sum+=log(norm)
                times_A+=(e==0)
            if t<steps-1:
                next_e = e if rng.random()<stay else 1-e
                if t>=burn_in:
                    same+=(next_e==e)
                e=next_e
        sample_length=steps-burn_in
        mean_log=log_sum/sample_length
        return {
            "steps":steps,
            "burn_in":burn_in,
            "seed":seed,
            "growth_per_generation":exp(mean_log),
            "empirical_log_growth":mean_log,
            "environment_A_fraction":times_A/sample_length,
            "environment_stay_frequency":same/(sample_length-1),
            "sampled_only_not_a_survival_proof": True,
        }

    def estimated_common_loss_boundary(self, sampled_growth):
        """d* ~ 1-(1-d)/exp(sampled_lambda), for common loss in both phases.

        The exact scaling of each path is known, but the input growth is
        sampled, so this threshold is NOT an exact statistical or theorem
        certificate. Do not treat its digits as an estimate for plants.
        """
        if len({p.loss for p in self.phases})!=1:
            raise ValueError("The loss threshold requires common loss")
        if not isinstance(sampled_growth,(float,int)) or sampled_growth<=0:
            raise ValueError("Supply a positive sampled growth multiplier")
        return 1-float(1-self.phases[0].loss)/sampled_growth


def symmetric_markov(persistence, trait_switch, loss=F(1,5)):
    """Matching ANTO v5 biological parameters; only E autocorrelation varies."""
    p=probability(persistence)
    w=probability(trait_switch)
    d=probability(loss)
    if not 0<p<1 or not 0<w<1 or not 0<=d<1:
        raise ValueError("Require ergodic ecology, interior trait switch and d<1")
    A=CapacityProcess(F(9,10),F(2,5),d,w,w)
    B=CapacityProcess(F(2,5),F(9,10),d,w,w)
    return MarkovCapacity((A,B),p)


def panel():
    """Source-locked, finite-horizon and sampled theoretical comparison."""
    rows=[]
    seeds=(20261020,20261021,20261022)
    for w in (F(1,10),F(1,2),F(9,10)):
        for stay in (F(1,10),F(9,10)):
            process=symmetric_markov(stay,w)
            path_samples=[process.sample_typical_growth(
                steps=100000,burn_in=2500,seed=s) for s in seeds]
            values=[x["growth_per_generation"] for x in path_samples]
            finite=process.finite_capacity_survival(14)
            if all(x>1.005 for x in values):
                direction="sampled_positive_log_growth"
            elif all(x<.995 for x in values):
                direction="sampled_negative_log_growth"
            else:
                direction="unresolved_near_zero_log_growth"
            mean=sum(values)/len(values)
            rows.append({
                "trait_switch":str(w),
                "environment_persistence":str(stay),
                "environment_A_stationary_fraction":.5,
                "environment_B_stationary_fraction":.5,
                "lag1_environment_autocorrelation":float(
                    process.environment_autocorrelation()),
                "mean_environment_run_length":float(
                    process.mean_environment_run_length()),
                "first_expression_from_C":float(
                    process.first_expression_probability_from_C()),
                "finite_horizon":14,
                "finite_genealogy_survival_from_C":
                    finite["survival_by_initial_type"][1],
                "finite_genealogy_survival_from_D":
                    finite["survival_by_initial_type"][0],
                "finite_path_count":finite["enumerated_histories"],
                "sampled_typical_growth_mean":mean,
                "sampled_typical_growth_range":[min(values),max(values)],
                "sampled_growth_direction_not_certified":direction,
                "estimated_loss_threshold_not_certified":
                    process.estimated_common_loss_boundary(mean),
                "sampling_steps_per_seed":100000,
                "sampling_seed_count":len(seeds),
                "sampled_paths_not_confidence_intervals":True,
            })
    return {
        "schema":"anto_slk_markov_environment_v6",
        "date":"2026-10-20",
        "environment_process":"shared stationary two-state Markov",
        "trait_states":"expressed_D or silently_competent_C",
        "no_new_functional_resource_after_loss":True,
        "observed_flower_trait_regain":False,
        "field_measured_environment_process":False,
        "formal_lyapunov_sign_certificate":False,
        "slk_selection_fitness_calibrated":False,
        "population_fixation_inferred":False,
        "organism_extinction_inferred":False,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "reproductive_slot_success_A":{"D":"9/10","C":"2/5"},
        "reproductive_slot_success_B":{"D":"2/5","C":"9/10"},
        "per_child_toolkit_loss":"1/5",
        "stationary_environment_occupancy":[.5,.5],
        "mean_arithmetic_matrix_growth":1.04,
        "initial_capacity_type_for_primary_event":"one C founder",
        "scenarios":rows,
    }
