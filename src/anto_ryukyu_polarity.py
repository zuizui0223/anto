"""Sensitivity of ancestral floral-pigment polarity in a two-tip Ryukyu pair.

States: W (white), P (pigmented). Observed *taxon-level* contrast:
brevicaule=W, irumtiense=P as reported by Chang et al. 2026.
No individual floral pigment chemistry or new ancestral trait is observed.

This two-tip continuous-time Markov model is a POLARITY DIAGNOSTIC, not a
fitted phylogenetic model. Positive branch times are illustrative; the
0.71-1.33 Ma uncertainty belongs to the LINEAGE split, NOT a colour event.
The model does not include outgroups, polymorphism, gene flow, sampling bias,
uncertainty in taxonomy or separate rate heterogeneity (all future work).
"""
import math


def transition_matrix(gain,loss,time):
    """Row-stochastic 2x2 CTMC probabilities for [W,P].

    W -> P = gain, P -> W = loss, rates per time unit.
    """
    for x in (gain,loss,time):
        if not math.isfinite(x) or x < 0: raise ValueError("Rates/time must be nonnegative")
    rate=gain+loss
    if rate==0: return ((1.0,0.0),(0.0,1.0))
    dec=math.exp(-rate*time)
    pi_p=gain/rate
    pi_w=loss/rate
    pwp=pi_p*(1-dec)
    ppw=pi_w*(1-dec)
    return ((1-pwp,pwp),(ppw,1-ppw))


def posterior_root_pigmented(gain,loss,t_brevicaule,t_irumtiense,prior_p=.5):
    """Root P posterior for observed tips (W in brev, P in irum)."""
    if not 0<prior_p<1 or not math.isfinite(prior_p):
        raise ValueError("Prior root P must be strictly between 0 and 1")
    a=transition_matrix(gain,loss,t_brevicaule)
    b=transition_matrix(gain,loss,t_irumtiense)
    likelihood_w=a[0][0]*b[0][1]
    likelihood_p=a[1][0]*b[1][1]
    total=(1-prior_p)*likelihood_w+prior_p*likelihood_p
    if total==0:
        raise ValueError("Observed tip contrast has zero likelihood under this rate model")
    return {"posterior_ancestral_pigmented":prior_p*likelihood_p/total,
            "likelihood_if_white_ancestor":likelihood_w,
            "likelihood_if_pigmented_ancestor":likelihood_p,
            "posterior_ancestral_white":(1-prior_p)*likelihood_w/total,
            "assumptions":"Independent 2-tip equal-rate process, no outgroups/reticulation"}


def sensitivity():
    """Declare a PRIOR-RATE grid, not a data-selected model."""
    estimates=[]
    for t in (.71,.93,1.33):
        for gain,loss in ((.1,.5),(.5,.5),(2.5,.5),(.5,2.5)):
            result=posterior_root_pigmented(gain,loss,t,t)
            estimates.append({"terminal_branch_time_ma":t,
                              "rate_white_to_pigmented_per_ma":gain,
                              "rate_pigmented_to_white_per_ma":loss,
                              **result})
    return {"evidence":"taxon-level two-tip white/pigmented contrast only",
            "source_split_ma":.93,"source_split_95hpd_ma":[.71,1.33],
            "not_phenotype_transition_date":True,
            "rates_hypothetical":True,
            "cases":estimates}


if __name__=="__main__":
    import json
    print(json.dumps(sensitivity(),indent=2,ensure_ascii=False))
