"""Exact environment-memory × forecast-horizon bridge for ANTO.

This analytic test uses the offspring mechanism already declared in the
SLK–ANTO v4/v6 modules. An ergodic two-state environment is shared by all
offspring of a generation. Stationary occupancy is 1/2 for each state
regardless of environmental persistence p.

Initially one nonexpressing but completely competent C parent. D denotes a
visibly expressing toolkit-competent parent. Toolkit-lost offspring leave
the counted genealogy but need not die as whole organisms.

Distinguish finite-horizon presence, immediate visible expression, and
eventual nonextinction: the first two are studied here, and neither
establishes fixation or permanent survival. All parameters are synthetic.
"""
from fractions import Fraction as F
from itertools import product
from anto_slk_markov_ecology_v6 import symmetric_markov

LOW=F(2,5)
HIGH=F(9,10)
KEEP=F(4,5)
HALF=F(1,2)


def probability(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a numeric probability")
    try:
        result=F(str(value)) if isinstance(value, float) else F(value)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite rational probability") from exc
    if not 0<=result<=1:
        raise ValueError(f"{name} must lie in [0,1]")
    return result


def _params(persistence, switching):
    p=probability(persistence,"persistence")
    w=probability(switching,"switching")
    if not 0<p<1 or not 0<w<1:
        raise ValueError("v6 bridge requires interior environmental persistence and switching")
    return p,w


def visible_child_from_silent(switching):
    """First D offspring event, averaging the identical stationary initial E."""
    w=probability(switching,"switching")
    if not 0<w<1:
        raise ValueError("Require interior switching")
    return 1-HALF*((1-KEEP*LOW*w)**2 + (1-KEEP*HIGH*w)**2)


def exact_competent_horizon(persistence, switching, horizon):
    """Exact finite-generation presence of competent C ancestry.

    Exhaustive rational enumeration of the original v6 offspring PGFs is
    bounded to H<=5. Not an infinite-time survival calculation.
    """
    p,w=_params(persistence,switching)
    if type(horizon) is not int or not 1<=horizon<=5:
        raise ValueError("horizon must be an integer from 1 to 5")
    process=symmetric_markov(p,w)
    extinct=F(0)
    mass=F(0)
    for sequence in product((0,1),repeat=horizon):
        wt=HALF
        for i,j in zip(sequence,sequence[1:]):
            wt*=p if i==j else 1-p
        z=(F(0),F(0))
        for state in reversed(sequence):
            z=process.phases[state].pgf(z)
        extinct+=wt*z[1]
        mass+=wt
    if mass!=1:
        raise ArithmeticError("stationary Markov environment mass must equal one")
    return 1-extinct


def direct_two_generation_survival(persistence, switching):
    """Closed two-environment-path formula, independent of v6 classes.

    For E0=i and E1=j, one C ancestor has zero competent descendants at
    generation 2 with probability:
      [1-KEEP*q_C,i + KEEP*q_C,i*(w*z_D,j +(1-w)*z_C,j)]**2,
    where z_{D/C,j}=[1-KEEP*q_{D/C,j}]**2.
    """
    p,w=_params(persistence,switching)

    def fail(i,j):
        ci=LOW if i==0 else HIGH
        dj=HIGH if j==0 else LOW
        cj=LOW if j==0 else HIGH
        zd=(1-KEEP*dj)**2
        zc=(1-KEEP*cj)**2
        slot=1-KEEP*ci+KEEP*ci*(w*zd+(1-w)*zc)
        return slot**2

    within=fail(0,0)+fail(1,1)
    cross=fail(0,1)+fail(1,0)
    return 1-HALF*(p*within+(1-p)*cross)


def persistence_slope(switching):
    """Exact dP(competent by H=2)/dp in the fixed A/B v6 witness.

    This slope factors by (1-2*w), so the environmental-persistence
    effect changes sign at w=1/2 under the declared symmetrical model.
    """
    w=probability(switching,"switching")
    if not 0<w<1:
        raise ValueError("Require interior switching")
    zl=(1-KEEP*LOW)**2
    zh=(1-KEEP*HIGH)**2
    return (HALF*KEEP*(HIGH-LOW)*(1-2*w)*(zl-zh)*
            (2-KEEP*(HIGH+LOW)*(2-zl-zh)))


def equal_history_brier_gain(prob_a,prob_b):
    """Exact history-information gain for a binary finite-horizon endpoint.

    Equal hypothetical history weights; Bayes risks, not fitted validation.
    """
    a=probability(prob_a,"prob_a")
    b=probability(prob_b,"prob_b")
    pooled=(a+b)/2
    coarse=pooled*(1-pooled)
    labelled=(a*(1-a)+b*(1-b))/2
    gain=(a-b)**2/4
    if coarse-labelled!=gain:
        raise ArithmeticError("Bayes risk identity failed")
    return {"coarse_brier":coarse,"history_brier":labelled,
            "history_information_gain":gain}


def bridge_panel():
    """Six rational scenarios: p=.1/.9, w=.1/.5/.9."""
    scenarios={}
    for w in (F(1,10),F(1,2),F(9,10)):
        p1,p9=F(1,10),F(9,10)
        outcomes={}
        for p in (p1,p9):
            outcomes[str(p)]={
                "first_visible_D_child":visible_child_from_silent(w),
                "competent_h1":exact_competent_horizon(p,w,1),
                "competent_h2":exact_competent_horizon(p,w,2),
                "competent_h3":exact_competent_horizon(p,w,3),
            }
        s1=outcomes[str(p1)]["competent_h2"]
        s9=outcomes[str(p9)]["competent_h2"]
        scenarios[str(w)]={
            "outcomes":outcomes,
            "difference_p9_minus_p1_at_h2":s9-s1,
            "exact_slope":persistence_slope(w),
            "brier_gain_at_h1":equal_history_brier_gain(
                outcomes[str(p1)]["competent_h1"],
                outcomes[str(p9)]["competent_h1"])["history_information_gain"],
            "brier_gain_at_h2":equal_history_brier_gain(s1,s9)["history_information_gain"],
            "common_stationary_environment_distribution":[HALF,HALF],
        }

    return {
        "schema":"anto_environment_horizon_identifiability_v17",
        "reference_model":"unchanged SLK-ANTO v4/v6 competent-genealogy branching PGFs",
        "forecast_endpoint":"finite-generation presence of at least one toolkit-competent descendant",
        "parent_type":"initial silent but reproductively competent C",
        "same_first_generation_probability":F(456,625),
        "no_natural_plant_parameters_estimated":True,
        "ancestral_trait_regain_measured":False,
        "no_permanent_genealogy_nonextinction_claim":True,
        "no_organism_population_extinction_claim":True,
        "SLK_payoff_calibrated":False,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
        "known_markov_transition_law_closes_history_proxy_by_construction":True,
        "scenarios":scenarios,
    }
