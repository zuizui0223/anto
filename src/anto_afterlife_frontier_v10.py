"""Exact finite afterlife audit for a *hypothetical* lost biological trait.

This is a conditional model, not a novelty claim about generic Markov chains,
not a population-genetic simulation and not fitted plant/pollinator evidence.
Initial visible state is absent for every lineage. A lineage may possess a
complete joint functional resource inventory (latent) or lack it (lost).
Each generation can show transient trait expression, a three-snapshot streak,
resource loss, or external acquisition/reassembly of missing resources.

The three-expression streak is a short observation endpoint, NOT reproductive
establishment, fixation, ecological function or permanent restoration.
"""
from dataclasses import dataclass
from fractions import Fraction

from anto_joint_accessibility_v7 import feasible_triple_bounds

ZERO = Fraction(0)
ONE = Fraction(1)


@dataclass(frozen=True)
class Rates:
    activation: Fraction
    repeat: Fraction
    degradation: Fraction
    acquisition: Fraction = ZERO

    def __post_init__(self):
        for name in ("activation", "repeat", "degradation", "acquisition"):
            try:
                value = Fraction(str(getattr(self, name)))
            except (TypeError, ValueError, ZeroDivisionError) as exc:
                raise ValueError(f"{name} must be a finite rational number") from exc
            if not ZERO <= value <= ONE:
                raise ValueError(f"{name} must be a probability in [0,1]")
            object.__setattr__(self, name, value)
        if self.activation + self.degradation > ONE:
            raise ValueError("Activation and degradation are competing latent-state events")
        if self.repeat + self.degradation > ONE:
            raise ValueError("Repeat expression and degradation are competing expressed-state events")


def transitions(state: str, rate: Rates):
    """One declared discrete-generation transition. Probabilities sum to one."""
    if state == "lost":
        return (("latent_acquired", rate.acquisition),
                ("lost", ONE-rate.acquisition))
    if state.startswith("latent_"):
        tag = state.removeprefix("latent_")
        return (("expr1_" + tag, rate.activation),
                ("lost", rate.degradation),
                (state, ONE-rate.activation-rate.degradation))
    if state.startswith(("expr1_", "expr2_", "expr3_")):
        two = state.startswith("expr2_")
        three = state.startswith("expr3_")
        tag = state.split("_", 1)[1]
        return ((("expr3_" if (two or three) else "expr2_")+tag, rate.repeat),
                ("lost", rate.degradation),
                ("latent_"+tag, ONE-rate.repeat-rate.degradation))
    raise ValueError("Invalid trait/resource state")


def _mass(states, *, which):
    return sum((p for (name, seen, route), p in states.items() if which(name, seen, route)), ZERO)


def run(theta, rates, horizon: int):
    """Exact joint state + ever-expressed propagation, including route provenance.

    theta = fraction of independent lineages with all necessary functional
    resources already co-retained. 1-theta starts without a complete toolkit.
    acquisition represents whole-toolkit recovery/reassembly from *outside*
    the stipulated original inventory. It is not a direct immigrant phenotype
    and is not a calibrated allele migration or molecular mutation model.
    """
    if type(horizon) is not int or horizon < 0:
        raise ValueError("horizon must be a nonnegative integer")
    try:
        theta = Fraction(str(theta))
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("theta must be rational") from exc
    if not ZERO <= theta <= ONE:
        raise ValueError("theta must be a probability")
    if isinstance(rates, Rates):
        schedule = (rates,)*horizon
    else:
        schedule = tuple(rates)
        if len(schedule) != horizon or not all(isinstance(r, Rates) for r in schedule):
            raise ValueError("Supply one Rates or a Rates sequence of exactly horizon entries")

    states = {("latent_native", False, ""): theta, ("lost", False, ""): ONE-theta}
    out = []
    for gen in range(horizon+1):
        if sum(states.values(), ZERO) != ONE:
            raise ArithmeticError("Mass not conserved")
        report = {
            "generation": gen,
            "first_expression_ever": _mass(states, which=lambda s,e,h: e),
            "present_expression": _mass(states, which=lambda s,e,h:
                s.startswith("expr")),
            "three_consecutive_expression_ever": _mass(states, which=lambda s,e,h: bool(h)),
            "three_streak_native": _mass(states, which=lambda s,e,h: h=="native"),
            "three_streak_acquired": _mass(states, which=lambda s,e,h: h=="acquired"),
            "whole_toolkit_missing_now": _mass(states, which=lambda s,e,h: s=="lost"),
            "latent_toolkit_unexpressed_ever": _mass(states, which=lambda s,e,h: s.startswith("latent") and not e),
            "visibly_white_now": _mass(states, which=lambda s,e,h: s=="lost" or s.startswith("latent")),
            "never_expressed_and_white_now": _mass(states, which=lambda s,e,h: not e),
        }
        out.append(report)
        if gen == horizon:
            break
        moved = {}
        for (state, seen, route), p in states.items():
            if not p:
                continue
            for next_state, w in transitions(state, schedule[gen]):
                newseen = seen or next_state.startswith("expr")
                next_route = route
                if not next_route and state.startswith("expr2_") and next_state.startswith("expr3_"):
                    next_route = state.split("_",1)[1]
                key = (next_state, newseen, next_route)
                moved[key] = moved.get(key, ZERO)+p*w
        states = moved
    return out


def exact_infinite_three_streak_no_acquisition(theta, rates: Rates):
    """P(ever seeing 3 consecutive positive snapshots) before closure.

    Homogeneous rates, without acquisition. For one initially competent
    lineage the elementary first-step equations give
      a*r**2 / (a*r**2 + d*(1+a*(1+r))).
    This is a standard absorbing competing-risk calculation, NOT a new
    probabilistic theorem or the establishment probability of real plants.
    """
    if rates.acquisition:
        raise ValueError("Closed-resource infinite formula requires acquisition=0")
    theta = Fraction(str(theta))
    if not ZERO <= theta <= ONE:
        raise ValueError("theta must be within [0,1]")
    a, r, d = rates.activation, rates.repeat, rates.degradation
    if not a or not r or not theta:
        return ZERO
    denominator = a*r*r + d*(ONE + a*(ONE+r))
    return theta*a*r*r/denominator


def closed_resource_success_score(rates: Rates):
    """Monotone ratio of activation/repeat to degradation at fixed d > 0.

    For two closed and homogeneous models with equal degradation d>0 and
    identical initial toolkit frequency, comparing this score compares
    eventual three-expression probabilities. It does not measure fitness.
    """
    if rates.acquisition:
        raise ValueError("A closed-resource score is undefined with acquisition")
    a,r=rates.activation,rates.repeat
    return a*r*r/(ONE+a*(ONE+r))


def eventual_three_streak_any_acquisition(theta, rates: Rates):
    """Asymptotic event probability, assuming constant non-exhausting inputs.

    With positive acquisition, activation and repeat, the finite chain reaches
    a first three-positive streak almost surely as T tends to infinity.
    No actual population extinction, limited external allele supply, source
    dependence or biological establishment is represented.
    """
    if not rates.acquisition:
        return exact_infinite_three_streak_no_acquisition(theta,rates)
    p=Fraction(str(theta))
    if not ZERO<=p<=ONE:
        raise ValueError("theta must lie in [0,1]")
    return ONE if rates.activation>0 and rates.repeat>0 else ZERO


def sharp_moment_bounds(rates, horizon: int,
                        marginals=(Fraction(1,2),)*3,
                        pairs=None, endpoint="three_consecutive_expression_ever"):
    """Sharp forecast interval given only first/pairwise functional moments.

    The interval is sharp under this finite model because the admissible
    inventory polytope permits every theta between its Fréchet endpoints,
    and the forecast probability is affine in theta. It is not a confidence
    interval or an estimate from a natural population.
    """
    if pairs is None:
        pairs = {(0,1): Fraction(1,4),(0,2): Fraction(1,4),(1,2): Fraction(1,4)}
    lo, hi = feasible_triple_bounds(marginals, pairs)
    left = run(lo, rates, horizon)[-1]
    right = run(hi, rates, horizon)[-1]
    if endpoint not in left or endpoint == "generation":
        raise ValueError("endpoint must be a named probability endpoint")
    xs = (left[endpoint], right[endpoint])
    return {"initial_joint_toolkit_lower": lo,
            "initial_joint_toolkit_upper": hi,
            "endpoint": endpoint,
            "forecast_lower": min(xs),
            "forecast_upper": max(xs),
            "forecast_width": max(xs)-min(xs),
            "slope_initial_joint_toolkit": (right[endpoint]-left[endpoint])/(hi-lo) if hi>lo else None,
            "initial_no_complete_toolkit_can_recover": left[endpoint]>0}


def independent_founder_opportunity(single_lineage_event, n_founders: int):
    """P(>=1 event) if whole lineage histories are independent, not photos."""
    if type(n_founders) is not int or n_founders<1:
        raise ValueError("n_founders must be a positive number of independent lineages")
    p = Fraction(str(single_lineage_event))
    if p<ZERO or p>ONE:
        raise ValueError("single_lineage_event must be a probability")
    return ONE-(ONE-p)**n_founders


def panel():
    """Four declared illustrative scenarios; no parameter inferred from plants."""
    theta = Fraction(1,8)
    models = {
        "fast_flash_weak_retention": Rates("3/5", "3/20", "2/25"),
        "slow_activation_strong_retention": Rates("1/5", "3/4", "2/25"),
        "slow_activation_with_acquisition": Rates("1/5", "3/4", "2/25", "1/20"),
        "silent_loss_no_recovery": Rates(0, "3/4", "2/25"),
    }
    results = {}
    for name, params in models.items():
        path = run(theta, params, 12)
        bounds = sharp_moment_bounds(params, 12)
        results[name] = {
            "assumptions": {k:str(getattr(params,k)) for k in ("activation","repeat","degradation","acquisition")},
            "initial_joint_fraction": str(theta),
            "trajectory": [{k:(v if k=="generation" else str(v)) for k,v in row.items()} for row in path],
            "sharp_bounds": {k:(str(v) if isinstance(v,Fraction) else v) for k,v in bounds.items()},
        }
        results[name]["infinite_horizon_three_streak_probability"] = str(
            eventual_three_streak_any_acquisition(theta, params))
        if not params.acquisition:
            results[name]["closed_resource_success_score"] = str(closed_resource_success_score(params))
    early = Rates("1/2", "4/5", "1/10")
    late = Rates("1/20", "3/20", "7/20")
    order = {
        "favorable_first": run(theta, [early]*6+[late]*6, 12)[-1]["three_consecutive_expression_ever"],
        "favorable_last": run(theta, [late]*6+[early]*6, 12)[-1]["three_consecutive_expression_ever"],
    }
    return {"synthetic":True, "no_real_plant_fitted_rates":True,
            "three_streak_not_fitness_establishment":True,
            "external_acquisition_not_verified_introgression":True,
            "current_visible_phenotype_absent_for_all_starting_lineages":True,
            "model_family":"finite discrete-generation resource / phenotype / consecutive-reappearance chain",
            "comparison_horizon_generations":12,
            "scenarios":results,
            "same_exposure_duration_order_test":{"favorable_first":str(order["favorable_first"]),
                "favorable_last":str(order["favorable_last"]),
                "order_difference":str(order["favorable_first"]-order["favorable_last"])} }
