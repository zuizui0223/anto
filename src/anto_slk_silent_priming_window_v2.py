"""Exact ANTO-SLK conditional *silent priming* ecological-window extension.

Theory-only. This introduces one NEW HYPOTHETICAL hidden regulator-readiness
state within a functionally complete, visibly trait-negative lineage.
An ecologically nonqualifying interval may shift an unprimed toolkit to a
primed but still invisible state. Resource loss independently eliminates it.
A qualifying SLK window then permits expression, with three-generation
continuation as a deliberately narrow event proxy (NOT establishment).

No data were used to estimate priming, loss, reactivation, or SLK signs.
"""
from dataclasses import dataclass
from fractions import Fraction as F

from anto_slk_ecological_window_v1 import SLKSlice

ZERO, ONE = F(0), F(1)


def _prob(value):
    try:
        x=F(str(value))
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("Each input must be a rational probability") from exc
    if not ZERO <= x <= ONE:
        raise ValueError("Probability must lie in [0,1]")
    return x


@dataclass(frozen=True)
class SilentPriming:
    """Conditional two-ecology one-toolkit model.

    Before the selective window, no visible phenotype is expressed:
      U -> P with rate b (silent readiness change);
      U or P -> D with rate d (permanent joint-toolkit loss);
      U -> U with 1-b-d; P -> P with 1-d.

    Within the 3-generation all-SLK-gates-positive window, first expression
    occurs with a_U from U or a_P from P; an expressing lineage continues
    to the next snapshot with unconditional probability r (competing
    expressed-state resource loss, deactivation, etc. are not separately
    identified here). For consistency with ANTO's competing hazards,
    a_U+d, a_P+d, and r+d must each be <= 1.

    The conditional readiness operator is an added MODEL HYPOTHESIS and is
    not an experimentally documented Cirsium or Petunia transition.
    """
    activation_unprimed: F
    activation_primed: F
    silent_priming: F
    resource_loss: F
    continuation: F

    def __post_init__(self):
        for name in ("activation_unprimed","activation_primed",
                     "silent_priming","resource_loss","continuation"):
            object.__setattr__(self,name,_prob(getattr(self,name)))
        d=self.resource_loss
        if self.silent_priming+d > ONE:
            raise ValueError("Pre-window priming and loss are competing")
        for name in ("activation_unprimed","activation_primed","continuation"):
            if getattr(self,name)+d > ONE:
                raise ValueError("Favorable expression and loss are competing")


DEFAULT_UNFAVORABLE=SLKSlice(F(12,5),F(3,2))  # Phi < 0
DEFAULT_FAVORABLE=SLKSlice(F(2,5),F(3,2))     # Phi,g0,Delta_R all > 0


def hidden_inventory(theta,params,wait):
    """Exact distribution after wait *nonqualifying*, nonexpressing steps."""
    theta=_prob(theta)
    if not isinstance(params,SilentPriming):
        raise TypeError("Expected a SilentPriming rate contract")
    if type(wait) is not int or wait < 0:
        raise ValueError("wait must be a nonnegative integer")
    s=ONE-params.resource_loss
    t=ONE-params.silent_priming-params.resource_loss
    u=theta*t**wait
    p=theta*(s**wait-t**wait)
    lost=ONE-theta*s**wait
    if min(u,p,lost)<ZERO or u+p+lost!=ONE:
        raise ArithmeticError("Latent capacity distribution not conserved")
    return {"unprimed_complete":u, "primed_complete":p,
            "missing_joint_toolkit":lost, "complete_total":u+p,
            "visible_trait_positive":ZERO}


def qualified_pulse(theta,params,wait,unfavorable=DEFAULT_UNFAVORABLE,
                    favorable=DEFAULT_FAVORABLE):
    """At least one 3-consecutive positive run in EXACTLY 3 favorable steps.

    No prior expression is allowed. All favorable steps have the SAME
    context (no time-varying rates or independent genetic restoration).
    The first good step must activate, and two further good steps must
    continue expression. Thus this is also an exact first-window event.
    The positive SLK signs are a candidate-specific *eligibility gate*
    rather than a demographic establishment probability.
    """
    if not isinstance(unfavorable,SLKSlice) or not isinstance(favorable,SLKSlice):
        raise TypeError("Contexts must be SLKSlice")
    if unfavorable.allows("all"):
        raise ValueError("Pre-window context must lack an SLK early gate")
    inventory=hidden_inventory(theta,params,wait)
    if not favorable.allows("all"):
        return ZERO
    return params.continuation**2*(
        params.activation_unprimed*inventory["unprimed_complete"]
       +params.activation_primed*inventory["primed_complete"])


def first_delay_boundary(params):
    """Exact one-step timing sign, apart from positive theta*r**2.

    F(1)-F(0) = theta*r^2*[b*(aP-aU)-aU*d].
    If d>0, b>0 and aP>aU, the discrete delay series is unimodal,
    because F(m+1)-F(m) changes sign at most once from + to -.
    No empirical threshold is identified.
    """
    if not isinstance(params,SilentPriming):
        raise TypeError("Expected SilentPriming")
    return (params.silent_priming*
            (params.activation_primed-params.activation_unprimed)
            -params.activation_unprimed*params.resource_loss)


def exact_delay_increment(theta,params,wait):
    """Direct closed-form first difference, independently testable."""
    hidden_inventory(theta,params,wait)
    s=ONE-params.resource_loss
    t=s-params.silent_priming
    d=params.resource_loss
    ap=params.activation_primed
    au=params.activation_unprimed
    return (_prob(theta)*params.continuation**2*
            ((ap-au)*(params.silent_priming+d)*t**wait
              - ap*d*s**wait))


def best_delay(theta,params,max_wait):
    """Earliest maximizing integer delay in a PREDECLARED finite horizon."""
    if type(max_wait) is not int or max_wait < 0:
        raise ValueError("max_wait must be a nonnegative integer")
    values=[qualified_pulse(theta,params,m) for m in range(max_wait+1)]
    idx=max(range(max_wait+1),key=lambda k:values[k])
    return {"best_wait":idx,"best_probability":values[idx],
            "earliest_probability":values[0],
            "last_probability":values[-1],
            "horizon_limited":True}



def identify_priming_from_three_pulses(p0,p1,p2,retention_ratio):
    """Identify *conditional* b and aP/aU using 0-,1-,2-wait cohorts.

    Assume shared initial toolkit theta, expression continuation r, and
    activation probabilities across separately allocated cohorts; identical
    pre-window loss hazard d with independently measured survival S=1-d.
    Normalize pulses y_m=Q(m)/Q(0), so theta and r cancel:
      v1=y1-S=b*(aP/aU-1)
      v2=y2-S**2=v1*(2*S-b)
      b=2*S-v2/v1; aP/aU=1+v1/b.
    When v1=0, no priming or no activation premium is observationally
    inseparable by these two extra pulse arms. Assay noise requires a
    statistical model; this function is structural algebra, NOT an estimator
    of real rates or confidence intervals.
    """
    p0,p1,p2,s=(_prob(x) for x in (p0,p1,p2,retention_ratio))
    base={"conditional_structural_design":True,
          "no_real_population_rates_inferred":True,
          "requires_independent_toolkit_survival":True}
    if p0==ZERO:
        return {**base,"identified":False,"status":"NO_BASELINE_EVENT"}
    v1=p1/p0-s
    v2=p2/p0-s*s
    if v1==ZERO:
        return {**base,"identified":False,"status":(
            "PRIMING_NOT_IDENTIFIED" if v2==ZERO
            else "INCOMPATIBLE_WITH_CONSTANT_RATE_MODEL"),
                "normalized_first_pulse_effect":v1}
    b=2*s-v2/v1
    if not ZERO < b <= s:
        return {**base,"identified":False,
                "status":"INCOMPATIBLE_PRIMING_RATE"}
    ratio=ONE+v1/b
    if ratio<ZERO:
        return {**base,"identified":False,
                "status":"INCOMPATIBLE_ACTIVATION_RATIO"}
    expected_p3=p0*(ratio*s**3-(ratio-ONE)*(s-b)**3)
    return {**base,"identified":True,
            "status":"STRUCTURALLY_IDENTIFIED_IF_ASSUMPTIONS_HOLD",
            "silent_priming":b,"activation_primed_over_unprimed":ratio,
            "normalized_first_pulse_effect":v1,
            "predicted_wait3_event_for_holdout":expected_p3}


def panel():
    theta=F(1,8)
    params=SilentPriming(
        activation_unprimed=F(1,10),
        activation_primed=F(3,5),
        silent_priming=F(3,20),
        resource_loss=F(1,25),
        continuation=F(4,5),
    )
    no_priming=SilentPriming(
        activation_unprimed=params.activation_unprimed,
        activation_primed=params.activation_primed,
        silent_priming=ZERO,
        resource_loss=params.resource_loss,
        continuation=params.continuation,
    )
    optimum=best_delay(theta,params,50)
    at_nine=qualified_pulse(theta,params,9)
    control=qualified_pulse(theta,no_priming,9)
    q0,q1,q2=(qualified_pulse(theta,params,m) for m in (0,1,2))
    ident=identify_priming_from_three_pulses(
        q0,q1,q2,ONE-params.resource_loss)
    # Same q0 and q1: higher priming rate versus smaller activation advantage.
    alt=SilentPriming(
        params.activation_unprimed,F(2,5),F(1,4),
        params.resource_loss,params.continuation)
    alt0,alt1,alt2=(qualified_pulse(theta,alt,m) for m in (0,1,2))
    if q0!=alt0 or q1!=alt1 or q2==alt2:
        raise AssertionError("Prespecified nonidentifiability example changed")
    return {
        "three_preexposure_arm_identification":ident,
        "same_pulse0_pulse1_different_pulse2":{
            "baseline_p0_both":q0,
            "one_delay_p1_both":q1,
            "world_A_priming":params.silent_priming,
            "world_A_primed_activation":params.activation_primed,
            "world_B_priming":alt.silent_priming,
            "world_B_primed_activation":alt.activation_primed,
            "world_A_p2":q2,
            "world_B_p2":alt2,
            "no_measured_biological_parameters":True,
        },
        "schema":"anto_slk_silent_priming_opportunity_v2",
        "model_type":"exact_conditional_synthetic_hypothesis",
        "no_new_specimens_or_empirical_rates":True,
        "priming_is_unverified_mechanism":True,
        "three_expression_snapshots_not_fitness_or_fixation":True,
        "SLK_positive_gate_is_not_population_establishment":True,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "initial_complete_toolkit":theta,
        "parameters":{
            "activation_unprimed":params.activation_unprimed,
            "activation_primed":params.activation_primed,
            "silent_priming":params.silent_priming,
            "resource_loss":params.resource_loss,
            "continuation":params.continuation,
        },
        "one_step_delay_boundary":first_delay_boundary(params),
        "one_step_delayed_better":first_delay_boundary(params)>0,
        "sampled_delay_event_probabilities":{
            str(m):qualified_pulse(theta,params,m)
            for m in (0,1,2,6,8,9,10,15,25,50)
        },
        "finite_optimum":optimum,
        "same_all_white_and_same_toolkit_at_delay9":{
            "complete_toolkit_both":hidden_inventory(theta,params,9)["complete_total"],
            "visible_positive_both":ZERO,
            "with_silent_priming":at_nine,
            "without_priming":control,
            "priming_divided_by_control":at_nine/control,
        },
    }
