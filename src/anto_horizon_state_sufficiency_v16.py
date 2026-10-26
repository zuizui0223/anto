"""Exact horizon-specific state-sufficiency counterexamples for floral afterlife.

This is a minimal one-route restriction of ANTO's more general afterlife
models. At t=0 the model plants all have the *same* white visible phenotype
and the same complete three-component potential toolkit. There is no
mutation, immigration, introgression, or induction cost.

An expression opportunity is tried every model generation with probability
a while competence remains. Following an opportunity that does not yield a
qualifying outcome, the shared complete toolkit remains competent until
the next generation with the specified interval-dependent probability r_t.
The first-visible endpoint treats any expression as qualifying; the
independently assessed demographic endpoint qualifies only when an
expressed descendant lineage is ultimately nonextinct in an explicit
Binomial(2, p) Galton-Watson process. A failed descendant attempt costs
nothing, so retries are permitted, as in ANTO v14.

If r_t differs between histories but was not included in a proposed
observable state, immediate endpoint predictions can be identical even
when later horizon-specific predictions diverge. Conditioning additionally
on the full transition schedule closes that history proxy by construction;
the conclusion does not challenge Markov closure or infer real plant rates.
"""
from dataclasses import dataclass
from fractions import Fraction as F


def _prob(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} is not a boolean")
    try:
        out = F(str(value)) if isinstance(value, float) else F(value)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite rational probability") from exc
    if out < 0 or out > 1:
        raise ValueError(f"{name} must be in [0,1]")
    return out


@dataclass(frozen=True)
class World:
    """All worlds share a white flower and initially complete toolkit.

    retention is a schedule of probabilities for keeping the still-silent
    full toolkit between consecutive modeled generations. This is not
    measured gene retention or a fitted environmental survival law.
    """
    induction: F
    daughter_viability: F
    retention: tuple

    def __post_init__(self):
        object.__setattr__(self, "induction", _prob(self.induction, "induction"))
        object.__setattr__(self, "daughter_viability",
                           _prob(self.daughter_viability, "daughter_viability"))
        try:
            rates = tuple(_prob(x, "retention") for x in self.retention)
        except TypeError as exc:
            raise ValueError("retention must be an iterable of probabilities") from exc
        object.__setattr__(self, "retention", rates)


def binomial_two_nonextinction(daughter_viability):
    """Ultimate nonextinction for offspring count Binomial(2, p).

    The smallest fixed point of G(z)=(1-p+p*z)^2 is extinction probability
    1 for p <= 1/2, or ((1-p)/p)^2 for p > 1/2. Infinite genealogy is a
    mathematical endpoint, not measured population persistence or fixation.
    """
    p = _prob(daughter_viability, "daughter_viability")
    return F(0) if p <= F(1, 2) else 1 - ((1-p)/p)**2


def qualifying_probability(world, horizon, endpoint):
    """Exact event by the first horizon model generations.

    Endpoint names: visible (any visible expression) and lineage
    (at least one expression by horizon whose independent daughter
    lineage will not go extinct, retrospectively identified).
    At an unqualified trial the full toolkit is unchanged except for
    subsequent background interval retention. No first-event preemption.
    """
    if endpoint not in ("visible", "lineage"):
        raise ValueError("endpoint must be visible or lineage")
    if isinstance(horizon, bool) or not isinstance(horizon, int) or horizon < 1:
        raise ValueError("horizon must be a positive integer")
    if horizon > len(world.retention) + 1:
        raise ValueError("retention schedule does not cover this horizon")
    attempt = world.induction
    if endpoint == "lineage":
        attempt *= binomial_two_nonextinction(world.daughter_viability)
    reach = F(1)
    result = F(0)
    for step in range(horizon):
        result += reach * attempt
        if step < horizon - 1:
            reach *= (1-attempt) * world.retention[step]
    return result


def exact_brier_contrast(prob_A, prob_B):
    """Expected Bayes Brier risks under two equally weighted histories.

    The coarse state records the same visible flower, complete toolkit,
    induction and offspring process but hides interval competence
    retention. There is no fitted prediction model or empirical holdout.
    """
    a = _prob(prob_A, "prob_A")
    b = _prob(prob_B, "prob_B")
    mean = (a+b)/2
    coarse = mean*(1-mean)
    with_history = (a*(1-a)+b*(1-b))/2
    gain = coarse-with_history
    assert gain == (a-b)**2/4
    return {
        "coarse_brier": coarse,
        "history_brier": with_history,
        "residual_history_gain": gain,
        "horizon_coarse_state_rejected": gain > 0,
    }


def compare_worlds(world_A, world_B, horizon, endpoint):
    """Return a history/forecast contrast for a specified future endpoint.

    Worlds must share immediate reproductive and phenotype parameters.
    Differences in the retention transition operator are deliberately
    hidden from the proposed current snapshot representation.
    """
    if (world_A.induction != world_B.induction or
            world_A.daughter_viability != world_B.daughter_viability):
        raise ValueError("comparison requires matched induction and offspring law")
    pa = qualifying_probability(world_A, horizon, endpoint)
    pb = qualifying_probability(world_B, horizon, endpoint)
    return {
        "prob_A": pa,
        "prob_B": pb,
        "difference_B_minus_A": pb-pa,
        **exact_brier_contrast(pa, pb),
    }


def reference_worlds(max_horizon=8):
    """Two preregistered mathematical witnesses, not natural cohorts.

    Constant contrast differs in each retention interval and therefore
    diverges at horizon 2. Delayed contrast has an identical FIRST
    retention interval; it first diverges at horizon 3.
    """
    if isinstance(max_horizon, bool) or not isinstance(max_horizon, int) or max_horizon < 3:
        raise ValueError("max_horizon must be at least three")
    slots = max_horizon-1
    half, quarter, three_quarters = F(1,2), F(1,4), F(3,4)
    return {
        "constant_hidden_retention": (
            World(half,F(2,3),(quarter,)*slots),
            World(half,F(2,3),(three_quarters,)*slots),
        ),
        "delayed_hidden_retention": (
            World(half,F(2,3),(half,quarter)+(half,)*(slots-2)),
            World(half,F(2,3),(half,three_quarters)+(half,)*(slots-2)),
        ),
    }


def exact_reference_panel(horizons=(1,2,3,5)):
    """Exact fraction-valued model comparisons, no empirical parameters."""
    horizons = tuple(horizons)
    if not horizons or min(horizons) < 1:
        raise ValueError("horizons must be nonempty and positive")
    worlds = reference_worlds(max(horizons))
    return {
        name: {
            str(h): {endpoint: compare_worlds(a,b,h,endpoint)
                     for endpoint in ("visible", "lineage")}
            for h in horizons
        } for name,(a,b) in worlds.items()
    }
