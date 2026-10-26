"""Exact two-type competent-genealogy response to non-symmetric ecological memory.

The only biological endpoint is at least one D/C toolkit-competent descendant
*at an exact finite generation*. Descendants without the requisite toolkit
may be viable organisms; this is not organismal extinction, natural ancestral
trait regain, or fixation. No empirical model parameters are estimated.

Extends ANTO v17's two-environment Galton--Watson witness with a stationary but
unequally occupied environment, unequal q_D/q_C responses, and rational
lag-one autocorrelation rho. The same environmental state is shared by all
individuals within any one generation.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product


def prob(value, label):
    if isinstance(value, bool):
        raise ValueError(f"{label} must be numeric")
    try:
        ans = F(str(value)) if isinstance(value, float) else F(value)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f"{label} must be a rational probability") from exc
    if not 0 <= ans <= 1:
        raise ValueError(f"{label} must be in [0,1]")
    return ans


@dataclass(frozen=True)
class CapacityRates:
    """Phase-specific slot success and shared phenotype-switching mechanism."""
    q_d_a: F
    q_c_a: F
    q_d_b: F
    q_c_b: F
    toolkit_retention: F
    switching: F

    def __post_init__(self):
        for name in ("q_d_a", "q_c_a", "q_d_b", "q_c_b", "toolkit_retention", "switching"):
            val = prob(getattr(self, name), name)
            if val <= 0 or (name != "toolkit_retention" and val >= 1):
                raise ValueError(f"{name} must be in (0,1), with toolkit_retention <=1")
            object.__setattr__(self, name, val)

    def q(self, environment):
        if environment == 0:
            return self.q_d_a, self.q_c_a
        if environment == 1:
            return self.q_d_b, self.q_c_b
        raise ValueError("environment index must be 0 or 1")

    def pgf(self, environment, extinction):
        """Unchanged ANTO v4/v6 two-independent-offspring-slot PGF."""
        qd, qc = self.q(environment)
        k, w = self.toolkit_retention, self.switching
        d, c = extinction
        return (
            (1-k*qd + k*qd*((1-w)*d+w*c))**2,
            (1-k*qc + k*qc*(w*d+(1-w)*c))**2,
        )


@dataclass(frozen=True)
class StationaryEcology:
    """Two-state shared Markov environment with fixed occupancy and rho.

    P(i,j) = rho * 1{i=j} + (1-rho) * pi_j.
    rho has the exact lag-one autocorrelation and is allowed to be negative
    when every transition remains strictly positive.
    """
    pi_a: F
    rho: F

    def __post_init__(self):
        pi = prob(self.pi_a, "pi_a")
        if not 0 < pi < 1:
            raise ValueError("pi_a must be in (0,1)")
        try:
            rho = F(str(self.rho)) if isinstance(self.rho, float) else F(self.rho)
        except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
            raise ValueError("rho must be a finite rational") from exc
        pi_b = 1-pi
        lower = -min(pi/pi_b, pi_b/pi)
        if not lower < rho < 1:
            raise ValueError("rho is outside strictly positive ergodic transitions")
        object.__setattr__(self, "pi_a", pi)
        object.__setattr__(self, "rho", rho)

    @property
    def stationary(self):
        return self.pi_a, 1-self.pi_a

    def transition(self, i, j):
        if i not in (0,1) or j not in (0,1):
            raise ValueError("environment indices must be 0 or 1")
        return (1-self.rho)*self.stationary[j]+(self.rho if i == j else 0)


def first_competent_child(rates, ecology):
    """A positive competent descendant at generation one from a silent C."""
    k = rates.toolkit_retention
    return sum(pi*(1-(1-k*rates.q(e)[1])**2)
               for e, pi in enumerate(ecology.stationary))


def first_visible_child(rates, ecology):
    """A D child from one C (not an established ancestral trait)."""
    k, w = rates.toolkit_retention, rates.switching
    return sum(pi*(1-(1-k*rates.q(e)[1]*w)**2)
               for e, pi in enumerate(ecology.stationary))


def _u_h(rates, e):
    """u: one C slot produces competent offspring;
    h: that offspring produces at least one competent child next generation,
    averaged over its D/C state assignment under an environment e.
    """
    k, w = rates.toolkit_retention, rates.switching
    qd, qc = rates.q(e)
    u = k*qc
    h = 1-w*(1-k*qd)**2-(1-w)*(1-k*qc)**2
    return u, h


def two_generation_path_extinction(rates, early_env, later_env):
    u = _u_h(rates, early_env)[0]
    h = _u_h(rates, later_env)[1]
    return (1-u*h)**2


def exact_horizon_presence(rates, ecology, horizon):
    """Exact rational P(any competent descendant at exactly generation H).

    Enumerates all 2**H shared-environment trajectories and computes the
    offspring PGF composition in chronological biological order. Deliberately
    restricted to H<=4 to control exact Fraction denominator growth.
    """
    if type(horizon) is not int or not 1 <= horizon <= 4:
        raise ValueError("horizon must be an integer from 1 to 4")
    loss = F(0)
    mass = F(0)
    for envs in product((0,1), repeat=horizon):
        wt = ecology.stationary[envs[0]]
        for a,b in zip(envs,envs[1:]):
            wt *= ecology.transition(a,b)
        extinction = (F(0),F(0))
        for e in reversed(envs):
            extinction = rates.pgf(e, extinction)
        mass += wt
        loss += wt * extinction[1]
    if mass != 1:
        raise ArithmeticError("Markov path weights do not sum to one")
    return 1-loss


def direct_horizon_two(rates, ecology):
    """Environment-pair formula independent of the multigeneration PGF."""
    loss = sum(ecology.stationary[i] * ecology.transition(i,j)
               * two_generation_path_extinction(rates,i,j)
               for i in (0,1) for j in (0,1))
    return 1-loss


def slope_factors(rates, ecology):
    """Exactly d S2 / d rho at fixed stationary occupancy.

    The response is a product of C founder output contrast, its descendants'
    environment-sensitive reproductive success contrast, and a nonlinear
    saturation correction from two independent offspring slots.
    """
    ua, ha = _u_h(rates,0)
    ub, hb = _u_h(rates,1)
    early = ua-ub
    later = ha-hb
    saturation = 2-(ua+ub)*(ha+hb)
    result = ecology.pi_a*(1-ecology.pi_a)*early*later*saturation
    return {"early_C_output_contrast":early,
            "later_competent_output_contrast":later,
            "two_slot_saturation_factor":saturation,
            "rho_slope":result}


def threshold_for_switching(rates):
    """Interior switching w* that cancels generation-two environment effect.

    Only applies to the h_A-h_B contrast, not to saturation or the initial
    C contrast; no interior threshold is guaranteed or universally w*=1/2.
    """
    k = rates.toolkit_retention
    g = lambda q: 1-(1-k*q)**2
    dd = g(rates.q_d_a)-g(rates.q_d_b)
    dc = g(rates.q_c_a)-g(rates.q_c_b)
    if dd == dc:
        return None
    root = dc/(dc-dd)
    return root if 0 < root < 1 else None


def brier_history_gain(p_a, p_b):
    """Exact equal-prior Bayes-Brier gain of an unobserved history label."""
    a, b = prob(p_a,"p_a"), prob(p_b,"p_b")
    mean = (a+b)/2
    coarse = mean*(1-mean)
    labelled = (a*(1-a)+b*(1-b))/2
    gain = (a-b)**2/4
    if coarse-labelled != gain:
        raise ArithmeticError("Brier decomposition failed")
    return gain


def witness_catalog():
    """Fixed non-calibrated witnesses and stringent counterexamples."""
    return {
        "v17_symmetric":(CapacityRates(F(9,10),F(2,5),F(2,5),F(9,10),F(4,5),F(1,10)),
                         F(1,2),F(-4,5),F(4,5)),
        "asymmetric_transition":(CapacityRates(F(85,100),F(35,100),F(45,100),F(75,100),F(4,5),F(7,13)),
                                 F(3,10),F(0),F(4,5)),
        "nonphenotype_environment":(CapacityRates(F(9,10),F(9,10),F(2,5),F(2,5),F(4,5),F(1,2)),
                                    F(1,2),F(0),F(4,5)),
        "high_retention_saturation":(CapacityRates(F(49,50),F(49,50),F(17,20),F(17,20),F(19,20),F(1,2)),
                                     F(1,2),F(0),F(4,5)),
        "low_retention_unsaturated":(CapacityRates(F(49,50),F(49,50),F(17,20),F(17,20),F(7,20),F(1,2)),
                                    F(1,2),F(0),F(4,5)),
    }


def panel():
    items={}
    for name,(rates,pi,rho_lo,rho_hi) in witness_catalog().items():
        lo,hi=StationaryEcology(pi,rho_lo),StationaryEcology(pi,rho_hi)
        slo,s_hi=slope_factors(rates,lo),slope_factors(rates,hi)
        assert slo==s_hi
        h1a,h1b=first_competent_child(rates,lo),first_competent_child(rates,hi)
        h2a,h2b=direct_horizon_two(rates,lo),direct_horizon_two(rates,hi)
        h3a,h3b=exact_horizon_presence(rates,lo,3),exact_horizon_presence(rates,hi,3)
        items[name]={
            "stationary_pi_a":pi,
            "rho_low":rho_lo,
            "rho_high":rho_hi,
            "first_visible_equal":first_visible_child(rates,lo)==first_visible_child(rates,hi),
            "h1_low":h1a,"h1_high":h1b,
            "h2_low":h2a,"h2_high":h2b,
            "h3_low":h3a,"h3_high":h3b,
            "h2_delta":h2b-h2a,
            "h3_delta":h3b-h3a,
            "analytic_h2_delta":(rho_hi-rho_lo)*slo["rho_slope"],
            "h1_history_brier_gain":brier_history_gain(h1a,h1b),
            "h2_history_brier_gain":brier_history_gain(h2a,h2b),
            "h3_history_brier_gain":brier_history_gain(h3a,h3b),
            "switching_null_threshold":threshold_for_switching(rates),
            **slo,
        }
    return {
        "schema":"anto_environment_memory_generalization_v18",
        "scope":"synthetic fixed-occupancy two-environment shared branching model",
        "endpoint":"presence of a complete-toolkit descendant at exactly generation H",
        "no_plant_calibration":True,
        "no_natural_trait_regain":True,
        "no_demographic_rescue_or_genetic_fixation":True,
        "no_SLK_ecological_payoff_measurement":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "w_half_is_not_a_universal_null":True,
        "horizon_two_null_not_a_horizon_three_null":True,
        "witnesses":items,
    }
