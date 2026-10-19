"""SLK × ANTO v3: latent readiness and ABSOLUTE demographic fate.

Constructed, uncalibrated ecological model. Reuse v2's U/P/toolkit-loss
kernel for a lineage that remains visibly negative before the switch.
At the switch, a capable lineage activates with a_U or a_P. ONLY as a
new explicit hypothesis, each activation seeds ONE *heritable* variant
founder. Each descendant has two independent Bernoulli reproduction
opportunities per favorable generation (Binomial(2,q_D) offspring).
A fixed resident reference has q_S and is NOT dynamically simulated.

Absolute growth (2*q_D) is distinct from the focal positive SLK
candidate-versus-resident sign. Threshold event at a finite deadline:
at least k DESCENDANTS of at least one activated founder among n
independent founding lineages. It is NOT genetic fixation, ecological
equilibrium, field evolutionary rescue or evidence of trait reacquisition.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as F

from anto_slk_ecological_window_v1 import SLKSlice
from anto_slk_silent_priming_window_v2 import (
    SilentPriming, hidden_inventory
)

ZERO, ONE = F(0), F(1)
DEFAULT_SLK = SLKSlice(F(2, 5), F(3, 2))


def _fraction(value):
    try:
        x = F(str(value))
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("Expected a rational offspring probability") from exc
    if not ZERO <= x <= ONE:
        raise ValueError("Offspring probabilities must lie in [0,1]")
    return x


def _decimal(value):
    if isinstance(value, F):
        return Decimal(value.numerator) / Decimal(value.denominator)
    return Decimal(str(value))


@dataclass(frozen=True)
class AbsoluteDemography:
    """One marked D founder versus an unsimulated S reference.

    SLK's declared Phi, local-release and rare-invasion signs must all
    be positive. Positivity only orders focal relative payoff; the
    offspring q_D,q_S are an EXTRA demographic calibration, not implied
    numerically by the payoff. q_D>q_S is necessary for a consistent
    directional comparison in this conditional model.

    An equality of q_D/q_S between scenarios gives an equal expected
    relative offspring ratio 2*q_D/(2*q_S), yet does not determine
    whether 2*q_D exceeds 1 (absolute demographic viability).
    """
    candidate_success: F
    resident_success: F
    independent_founders: int = 40
    descendant_target: int = 4
    fixed_census_generation: int = 14
    slk: SLKSlice = DEFAULT_SLK

    def __post_init__(self):
        q = _fraction(self.candidate_success)
        s = _fraction(self.resident_success)
        object.__setattr__(self, "candidate_success", q)
        object.__setattr__(self, "resident_success", s)
        if not isinstance(self.slk, SLKSlice) or not self.slk.allows("all"):
            raise ValueError("SLK value/path/rare-invasion signs must all be positive")
        if s == ZERO or q <= s:
            raise ValueError("Require positive resident reference and q_D>q_S")
        for name in ("independent_founders", "descendant_target",
                     "fixed_census_generation"):
            val = getattr(self, name)
            if type(val) is not int or val < (0 if name == "fixed_census_generation" else 1):
                raise ValueError("Counts/horizon must be valid integers")


def truncated_poly_product(a, b, limit):
    """Convolution retaining only the coefficients below z**limit."""
    c = [Decimal(0)]*limit
    for i, x in enumerate(a[:limit]):
        if not x:
            continue
        for j, y in enumerate(b[:limit-i]):
            if y:
                c[i+j] += x*y
    return c


def one_founder_at_least_k(q, generations, k, digits=65):
    """PGF exact polynomial recurrence, numerically evaluated in Decimal.

    Start G_0(z)=z, offspring PGF f(z)=(1-q+q*z)^2, and
    G_{t+1}(z)=f(G_t(z)). To find P(N_t>=k), coefficients of degrees
    0...k-1 suffice. The recurrence is an EXACT finite branching-model
    identity; Decimal rounding means numerical answers are NOT exact
    fractions. A separate full-distribution oracle tests short horizons.
    """
    q = _fraction(q)
    if type(generations) is not int or generations < 0:
        raise ValueError("Generations must be nonnegative integer")
    if type(k) is not int or k < 1:
        raise ValueError("Descendant threshold must be positive integer")
    if type(digits) is not int or digits < 32:
        raise ValueError("Decimal precision must be >=32")
    with localcontext() as ctx:
        ctx.prec = digits
        p = _decimal(q)
        coef = [Decimal(0)]*k
        if k > 1:
            coef[1] = Decimal(1)
        for _ in range(generations):
            b = [p*x for x in coef]
            b[0] += 1-p
            coef = truncated_poly_product(b, b, k)
        tail = Decimal(1)-sum(coef)
        if tail < -Decimal("1e-35") or tail > Decimal(1)+Decimal("1e-35"):
            raise ArithmeticError("Branching PGF lost probability mass")
        return max(Decimal(0), min(Decimal(1), +tail))


def first_activation_after_silent_history(theta, priming, wait):
    """Conditional probability of one expressed founder at favorable switch.

    This is ANTO v2's positive first step BEFORE multiplying by its
    additional r**2 three-expression observation endpoint.
    Assumes an activation can create a heritable D founder, which is a
    new *hypothesis* and never follows from one pigment photograph.
    """
    if not isinstance(priming, SilentPriming):
        raise TypeError("Expected SilentPriming model")
    p = hidden_inventory(theta, priming, wait)
    return (priming.activation_unprimed*p["unprimed_complete"]+
            priming.activation_primed*p["primed_complete"])


def census_probability(theta, priming, model, wait, digits=65):
    """Probability at least one of n independent founders has >=k D descendants.

    Each independent initial lineage can activate ONCE when the ecology
    becomes favorable. It then runs a separate branching process for
    remaining time (fixed census generation minus delay). Multiple D
    founders do not cooperate or compete, and S is only a reference.
    """
    if not isinstance(model, AbsoluteDemography):
        raise TypeError("Expected AbsoluteDemography model")
    if type(wait) is not int or not 0 <= wait <= model.fixed_census_generation:
        raise ValueError("Wait must be within the predetermined census horizon")
    a = first_activation_after_silent_history(theta, priming, wait)
    remaining = model.fixed_census_generation-wait
    with localcontext() as ctx:
        ctx.prec = digits
        tail = one_founder_at_least_k(
            model.candidate_success, remaining, model.descendant_target,
            digits=digits)
        per_founder = _decimal(a)*tail
        observed = Decimal(1)-(Decimal(1)-per_founder)**model.independent_founders
        return {
            "pre_exposure_generations": wait,
            "available_growth_generations": remaining,
            "first_heritable_variant_activation_hypothetical": a,
            "one_activated_founder_descendant_threshold": tail,
            "one_original_founder_joint_event": +per_founder,
            "at_least_one_of_n_independent_lineages": +observed,
            "relative_expected_offspring_ratio": (
                model.candidate_success/model.resident_success),
            "candidate_absolute_mean_offspring": 2*model.candidate_success,
            "resident_absolute_mean_offspring": 2*model.resident_success,
        }


def infinite_horizon_nonextinction(q):
    """Exact classic Galton-Watson survival in constant environment.

    This infinite-time limit is NOT what the finite census measures.
    Environmental permanence and inherited functional phenotype are assumed.
    """
    q = _fraction(q)
    if q <= F(1, 2):
        return ZERO
    return ONE-((ONE-q)/q)**2


def censused_optimum(theta, priming, model, digits=65):
    """Earliest max of a finite, prespecified census among all wait choices."""
    series = [
        census_probability(theta, priming, model, i, digits=digits)
        for i in range(model.fixed_census_generation+1)
    ]
    best = max(range(len(series)),
               key=lambda i: series[i]["at_least_one_of_n_independent_lineages"])
    return {
        "earliest_best_pre_exposure":best,
        "best_census_event":series[best]["at_least_one_of_n_independent_lineages"],
        "all_waits":series,
        "best_event_is_not_fixation_or_rescue":True,
    }


def panel():
    theta = F(1, 8)
    priming = SilentPriming(F(1,10), F(3,5), F(3,20),
                           F(1,25), F(4,5))
    # The two ecologies have identical relative offspring ratios 13/11
    # but different ABSOLUTE demography.
    supercritical = AbsoluteDemography(F(13,20),F(11,20))
    subcritical = AbsoluteDemography(F(39,100),F(33,100))
    if (supercritical.candidate_success/supercritical.resident_success !=
        subcritical.candidate_success/subcritical.resident_success):
        raise ArithmeticError("Relative selection comparator changed")
    cases = {
        "absolute_supercritical":censused_optimum(theta,priming,supercritical),
        "absolute_subcritical":censused_optimum(theta,priming,subcritical),
    }
    return {
        "schema":"anto_slk_absolute_demography_v3",
        "date":"2026-10-19",
        "synthetic_conditional_no_field_rates":True,
        "activation_heritable_is_added_unvalidated_assumption":True,
        "relative_advantage_not_population_rescue":True,
        "finite_threshold_not_establishment_or_fixation":True,
        "no_density_dependence_or_interacting_founders":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "theta":theta,
        "silent_priming":{
            "a_U":priming.activation_unprimed,
            "a_P":priming.activation_primed,
            "b":priming.silent_priming,
            "d":priming.resource_loss,
            "r_observation_only_not_used":priming.continuation,
        },
        "favorable_slk":{"phi":supercritical.slk.phi,
                         "local_release":supercritical.slk.local_release,
                         "rare_margin":supercritical.slk.rare_margin},
        "same_relative_expected_offspring_ratio":F(13,11),
        "population_reference_conditions":{
            "absolute_supercritical":{"qD":supercritical.candidate_success,
                                      "qS":supercritical.resident_success,
                                      "candidate_nonextinction_under_permanent_environment":
                                          infinite_horizon_nonextinction(supercritical.candidate_success)},
            "absolute_subcritical":{"qD":subcritical.candidate_success,
                                    "qS":subcritical.resident_success,
                                    "candidate_nonextinction_under_permanent_environment":
                                        infinite_horizon_nonextinction(subcritical.candidate_success)}
        },
        "founders_per_scenario":40,
        "descendant_target":4,
        "fixed_census_generation":14,
        "scenarios":cases
    }
