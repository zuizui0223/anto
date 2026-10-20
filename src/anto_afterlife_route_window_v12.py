"""Two recovery routes after a phenotypically silent period (synthetic only).

A common necessary resource C supports ancestral A and alternative B routes.
Resource loss is irreversible; each block may represent several independently
lost, jointly necessary components with their hazards summed. A and B are
alternative *first re-expression events*, not measured fitness or fixation.
The absorbing labels record which event happened first; they do not assert
that biological recovery is permanent or that later changes are impossible.

This applies standard exponential-loss / competing-risk mathematics to the
v11 route distinction. It is NOT a new general probability theorem, a fitted
plant model, or evidence that natural loss causes biological innovation.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import exp, isfinite, log
from typing import Any

import numpy as np
from scipy.linalg import expm

# Resource states and first-event records; all start visibly absent.
STATES = ("both_routes", "ancestral_only", "alternative_only",
          "ancestral_first", "alternative_first", "capacity_lost")


def _rate(value: Any) -> Fraction:
    if isinstance(value, bool):
        raise ValueError("Boolean is not a rate")
    try:
        f = Fraction(str(value))
        valid = f >= 0 and isfinite(float(f))
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError("Rates must be finite nonnegative numbers") from exc
    if not valid:
        raise ValueError("Rates must be finite nonnegative numbers")
    return f


def _time(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Time must be finite and nonnegative") from exc
    if isinstance(value, bool) or not isfinite(result) or result < 0:
        raise ValueError("Time must be finite and nonnegative")
    return result


@dataclass(frozen=True)
class Loss:
    common: Fraction
    ancestral: Fraction
    alternative: Fraction

    def __post_init__(self):
        for key in ("common", "ancestral", "alternative"):
            object.__setattr__(self, key, _rate(getattr(self, key)))


@dataclass(frozen=True)
class Rates:
    ancestral: Fraction
    alternative: Fraction
    loss: Loss

    def __post_init__(self):
        for key in ("ancestral", "alternative"):
            object.__setattr__(self, key, _rate(getattr(self, key)))
        if not isinstance(self.loss, Loss):
            raise ValueError("Specify a Loss object")


def _divide(numerator, denominator):
    return numerator / denominator if denominator else Fraction(0)


def exact_from_both(rates: Rates) -> dict[str, Fraction]:
    """Exact eventual first-route probabilities, starting with C,A,B intact."""
    a, b = rates.ancestral, rates.alternative
    c, x, y = rates.loss.common, rates.loss.ancestral, rates.loss.alternative
    ua, ub = _divide(a, a+c+x), _divide(b, b+c+y)
    total_rate = a+b+c+x+y
    fa = _divide(a+y*ua, total_rate)
    fb = _divide(b+x*ub, total_rate)
    return {"ancestral": fa, "alternative": fb, "any": fa+fb,
            "only_ancestral": ua, "only_alternative": ub}


def generator(rates: Rates) -> np.ndarray:
    """Six-state row generator; p(t) = p(0) @ expm(Q*t)."""
    a,b = map(float, (rates.ancestral,rates.alternative))
    c,x,y = map(float, (rates.loss.common,rates.loss.ancestral,rates.loss.alternative))
    q = np.zeros((6,6), dtype=float)
    for i,j,r in ((0,1,y),(0,2,x),(0,3,a),(0,4,b),(0,5,c),
                  (1,3,a),(1,5,c+x),(2,4,b),(2,5,c+y)):
        q[i,j] += r
        q[i,i] -= r
    return q


def aged_inventory(loss: Loss, age: float) -> np.ndarray:
    """Independent resource loss before opportunity; no pre-window expression.

    Initially C,A,B are all intact; the age is model time, not observed years
    or an ancestor-state estimate. Any loss of C disables both routes.
    """
    t = _time(age)
    pc,pa,pb = (exp(-float(r)*t) for r in (loss.common,loss.ancestral,loss.alternative))
    both=pc*pa*pb; onlya=pc*pa*(1-pb); onlyb=pc*(1-pa)*pb
    return np.array([both,onlya,onlyb,0.,0.,max(0.,1-both-onlya-onlyb)])


def after_delay(rates: Rates, pre_loss: Loss, age: float,
                horizon: float | None = None) -> dict[str, float | None]:
    """First re-expression outcomes after silence, for finite or eventual follow-up."""
    p0=aged_inventory(pre_loss, age)
    if horizon is None:
        values=exact_from_both(rates)
        fa=p0[0]*float(values["ancestral"])+p0[1]*float(values["only_ancestral"])
        fb=p0[0]*float(values["alternative"])+p0[2]*float(values["only_alternative"])
    else:
        final=p0 @ expm(generator(rates)*_time(horizon))
        fa,fb=float(final[3]),float(final[4])
    any_event=fa+fb
    return {"ancestral_first":fa, "alternative_first":fb, "any_first":any_event,
            "no_recovery_by_endpoint":max(0.,1-any_event),
            "alternative_fraction_among_first_recoveries":fb/any_event if any_event>0 else None}


def alternative_peak(rates: Rates, pre_loss: Loss) -> dict[str, float | str | None]:
    """Peak of the eventual ALTERNATIVE-FIRST probability versus pre-silence age.

    B(t) = exp(-k*t) * [u - v*exp(-x*t)],
    k=pre_C+pre_B; x=pre_A; u=B-only success; v=u-B-from-both.
    An interior maximum exists iff x*v > k*(u-v), with k,x positive.
    This is a conditional first-event effect, not faster intrinsic B activation.
    """
    exact=exact_from_both(rates)
    u=float(exact["only_alternative"])
    v=float(exact["only_alternative"]-exact["alternative"])
    k=float(pre_loss.common+pre_loss.alternative)
    x=float(pre_loss.ancestral)
    if u==0:
        return {"kind":"constant_zero", "age":None}
    if k==0:
        return {"kind":"increasing_to_limit" if x*v>0 else "constant", "age":None}
    if x>0 and x*v>k*(u-v):
        return {"kind":"interior", "age":log((k+x)*v/(k*u))/x}
    return {"kind":"boundary_at_zero", "age":0.}


def alternative_ever_without_preemption(rates: Rates, pre_loss: Loss,
                                        age: float, horizon: float | None = None) -> float:
    """Negative control: allow B after an earlier A event, with unchanged B resources.

    A events do not alter C/B loss or activation. B-ever no longer depends
    on the ancestral activation/loss rates; the first-event bump disappears.
    This is NOT the same endpoint as the route-specific cumulative incidence.
    """
    survival=exp(-float(pre_loss.common+pre_loss.alternative)*_time(age))
    b=rates.alternative; total=b+rates.loss.common+rates.loss.alternative
    result=survival*float(_divide(b,total))
    if horizon is not None:
        h=_time(horizon)
        result *= -np.expm1(-float(total)*h)
    return float(result)


def parameter_grid() -> list[dict]:
    """A deterministic 1125-setting sensitivity grid, NOT a parameter prior."""
    out=[]
    for a,b,c,x,y in product((.05,.1,.2,.4,.8),(.05,.1,.2,.4,.8),
                            (.005,.02,.08),(.01,.05,.1,.3,.6),(.005,.02,.08)):
        loss=Loss(c,x,y); rates=Rates(a,b,loss)
        peak=alternative_peak(rates,loss)
        p0=after_delay(rates,loss,0)
        at=float(peak["age"]) if peak["age"] is not None else 0.
        pt=after_delay(rates,loss,at)
        out.append({"activation_A":a,"activation_B":b,"loss_C":c,"loss_A":x,"loss_B":y,
                    "peak_kind":peak["kind"],"peak_age":peak["age"],
                    "alternative_at_zero":p0["alternative_first"],
                    "alternative_at_peak":pt["alternative_first"],
                    "any_at_zero":p0["any_first"],"any_at_peak":pt["any_first"]})
    return out
