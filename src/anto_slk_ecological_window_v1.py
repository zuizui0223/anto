"""SLK × ANTO conditional ecological-opportunity window (theory-only).

An ANTO lineage may be latent, visibly expressing, or resource-lost.
An SLK *candidate-specific* comparison supplies three distinct strict signs:
net endpoint value, local release, and rarity fitness. A successful event here
is ONLY three consecutive expressed snapshots during generations when the
selected signs hold. It is NOT establishment, fixation, heritable trait regain,
or any claim that an actual plant's structural and floral phenotypes coincide.
"""
from dataclasses import dataclass
from fractions import Fraction as F

from anto_afterlife_frontier_v10 import Rates, transitions, run as afterlife_run
from anto_joint_accessibility_v7 import feasible_triple_bounds


def _fraction(value):
    try:
        return F(str(value))
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("Expected a finite rational quantity") from exc


@dataclass(frozen=True)
class SLKSlice:
    """The REGISTERED witness R(x)=x+x^2, K(x)=k*x, x in [0,1].

    Phi=2-k, local release gradient g0=1-k, and rare-D selection
    Delta_R=Phi-eta under the canonical symmetric frequency map.
    k is nonnegative; eta can have either sign.
    """
    k: F
    eta: F

    def __post_init__(self):
        k, eta = _fraction(self.k), _fraction(self.eta)
        if k < 0:
            raise ValueError("Architecture debit k must be nonnegative")
        object.__setattr__(self, "k", k)
        object.__setattr__(self, "eta", eta)

    @property
    def phi(self):
        return F(2)-self.k

    @property
    def local_release(self):
        return F(1)-self.k

    @property
    def rare_margin(self):
        return self.phi-self.eta

    def allows(self, gate="all"):
        if gate not in ("value", "path", "all"):
            raise ValueError("gate must be value, path, or all")
        return self.phi > 0 and (
            gate == "value" or
            (self.local_release > 0 and
             (gate == "path" or self.rare_margin > 0))
        )


def _schedule(rates, length):
    if isinstance(rates, Rates):
        return (rates,)*length
    seq = tuple(rates)
    if len(seq) != length or any(not isinstance(x, Rates) for x in seq):
        raise ValueError("Provide Rates or one Rates per ecological generation")
    return seq


def opportunity_trace(theta, rates, contexts, gate="all"):
    """Exact event probabilities, with fractions and ANTO's state transitions.

    A qualifying event requires three successive visibly expressed snapshots,
    each occurring while the declared SLK gates have strictly positive signs.
    Ecology changes the observation/eligibility mask ONLY. It does not change
    ANTO expression, genetic loss, or transition rates by fiat.

    First qualified streak provenance is native versus acquired toolkit;
    every previously qualified event stays counted even if resources later
    disappear. This is an occurrence event, NOT durable evolutionary success.
    """
    contexts = tuple(contexts)
    if not contexts or any(not isinstance(c, SLKSlice) for c in contexts):
        raise ValueError("Supply at least one SLKSlice per generation")
    if gate not in ("value", "path", "all"):
        raise ValueError("gate must be value, path, or all")
    length = len(contexts)
    schedule = _schedule(rates, length)
    plain = afterlife_run(theta, schedule, length)  # validates theta and Rates
    theta = _fraction(theta)
    states = {
        ("latent_native", 0, False, ""): theta,
        ("lost", 0, False, ""): F(1)-theta,
    }

    def snapshot(t, context):
        total = sum(states.values(), F(0))
        if total != F(1):
            raise ArithmeticError("Transition probability mass not conserved")
        return {
            "generation": t,
            "ecological_gate_open": None if context is None else context.allows(gate),
            "visible_expression_now": sum((
                p for (state, _, _, _), p in states.items()
                if state.startswith("expr")
            ), F(0)),
            "any_expression_ever": sum((
                p for (_, _, seen, _), p in states.items() if seen
            ), F(0)),
            "unqualified_three_streak_ever": plain[t]["three_consecutive_expression_ever"],
            "gate_qualified_three_streak_ever": sum((
                p for (_, _, _, route), p in states.items() if route
            ), F(0)),
            "qualified_native": sum((
                p for (_, _, _, route), p in states.items() if route == "native"
            ), F(0)),
            "qualified_acquired": sum((
                p for (_, _, _, route), p in states.items() if route == "acquired"
            ), F(0)),
        }

    out = [snapshot(0, None)]
    for t, (context, rate) in enumerate(zip(contexts, schedule), 1):
        moved = {}
        open_now = context.allows(gate)
        for (state, streak, seen, route), mass in states.items():
            if mass == 0:
                continue
            for next_state, weight in transitions(state, rate):
                if weight == 0:
                    continue
                expressed = next_state.startswith("expr")
                next_streak = min(3, streak+1) if (expressed and open_now) else 0
                next_route = route
                if not route and next_streak == 3:
                    next_route = next_state.split("_", 1)[1]
                key = (next_state, next_streak, seen or expressed, next_route)
                moved[key] = moved.get(key, F(0)) + mass*weight
        states = moved
        current = snapshot(t, context)
        if current["gate_qualified_three_streak_ever"] > current["unqualified_three_streak_ever"]:
            raise ArithmeticError("An ecological subset exceeded raw ANTO expression")
        if current["qualified_native"]+current["qualified_acquired"] != current["gate_qualified_three_streak_ever"]:
            raise ArithmeticError("First-event origin partition failed")
        out.append(current)
    return out


def sharp_opportunity_bounds(rates, contexts, gate="all", marginals=None, pairs=None):
    """Sharp model-conditional range when ONLY initial P(111) is unknown.

    Each incomplete initial joint inventory is grouped into the existing
    ANTO lost state with the SAME exogenous transition/acquisition kernel.
    The first/pairwise constraints alone yield theta in a Frechet interval.
    The event probability is affine in theta because the Markov dynamics are
    fixed and founder lineages do not alter one another's environment.
    This is NOT a sampling CI or an empirical recovery-rate estimate.
    """
    if marginals is None:
        marginals = (F(1,2),)*3
    if pairs is None:
        pairs = {(0,1):F(1,4), (0,2):F(1,4), (1,2):F(1,4)}
    contexts = tuple(contexts)
    low, high = feasible_triple_bounds(marginals, pairs)
    name = "gate_qualified_three_streak_ever"
    left = opportunity_trace(low, rates, contexts, gate)[-1][name]
    right = opportunity_trace(high, rates, contexts, gate)[-1][name]
    return {
        "joint_toolkit_lower": low,
        "joint_toolkit_upper": high,
        "forecast_lower": min(left, right),
        "forecast_upper": max(left, right),
        "forecast_at_lower_theta": left,
        "forecast_at_upper_theta": right,
        "affine_slope": (right-left)/(high-low) if high>low else None,
        "model_conditional_not_confidence_interval": True,
    }


def gradient_slices(horizon=14):
    """Cost-lowering example; cost saturates at zero after generation 12."""
    if type(horizon) is not int or horizon < 1:
        raise ValueError("horizon must be a positive integer")
    return tuple(SLKSlice(max(F(0), F(12,5)-F(t,5)), F(3,2))
                 for t in range(1, horizon+1))


def panel():
    """Exactly reproducible, explicitly synthetic SLK × ANTO examples."""
    horizon = 14
    theta = F(1,8)
    rates = Rates(F(2,5), F(3,4), F(1,10), F(0))
    favorable = SLKSlice(F(2,5), F(3,2))
    unfavorable = SLKSlice(F(12,5), F(3,2))
    early = (favorable,)*5+(unfavorable,)*9
    late = (unfavorable,)*9+(favorable,)*5
    gradient = gradient_slices(horizon)
    outcome = "gate_qualified_three_streak_ever"
    e = opportunity_trace(theta, rates, early)
    l = opportunity_trace(theta, rates, late)
    v = opportunity_trace(theta, rates, gradient, "value")
    p = opportunity_trace(theta, rates, gradient, "path")
    a = opportunity_trace(theta, rates, gradient, "all")
    return {
        "schema": "anto_slk_ecological_opportunity_window_v1",
        "synthetic_conditional_model": True,
        "actual_evolutionary_regain_or_fixation_observed": False,
        "qualified_expression_is_not_establishment": True,
        "Cirsium_empirical_owner": "zuizui0223/aza3",
        "theta": theta,
        "rates": {"activation": rates.activation, "repeat": rates.repeat,
                  "degradation": rates.degradation, "acquisition": rates.acquisition},
        "horizon_generations": horizon,
        "equal_favorable_generations": 5,
        "early_qualified": e[-1][outcome],
        "late_qualified": l[-1][outcome],
        "early_late_ratio": e[-1][outcome]/l[-1][outcome],
        "early_any_expression": e[-1]["any_expression_ever"],
        "late_any_expression": l[-1]["any_expression_ever"],
        "gradient_gate_counts": {
            gate: sum(ctx.allows(gate) for ctx in gradient)
            for gate in ("value", "path", "all")
        },
        "gradient_probabilities": {
            "value": v[-1][outcome],
            "value_and_local_path": p[-1][outcome],
            "value_path_and_rare_invader": a[-1][outcome],
            "raw_no_ecological_gate": a[-1]["unqualified_three_streak_ever"],
        },
        "gradient_first_positive_generation": {
            gate: next((i for i, c in enumerate(gradient, 1) if c.allows(gate)), None)
            for gate in ("value", "path", "all")
        },
        "late_sharp_partial_identification": sharp_opportunity_bounds(
            rates, late, "all"
        ),
    }
