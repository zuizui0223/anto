"""Organ-aware pathway capability model (illustrative, not a Cirsium fit).

Phenotype and potential are *not* the same molecular state.

States encode two abstract binary modules:
    C: functionally intact structural pigment toolkit (not just a gene hit)
    F: effective flower-specific activation circuit
Index 0: C1 F1  = floral pigment expressed
Index 1: C1 F0  = white, structural toolkit intact, flower regulation OFF
Index 2: C0 F1  = white, structure unavailable but floral regulation ON
Index 3: C0 F0  = white, both modules unavailable

An external measured NONFLORAL USE indicator U in [0,1] affects only
structural degradation, via a prespecified hypothesized pleiotropic-constraint
coefficient beta. U is NOT estimated from Park et al. 2020 qPCR and cannot
be inferred for the focal Ryukyu species from those organ contrasts.

This is a basic mechanistic counterexample to one-dimensional W+/W-.
It omits individual loci, gene flow, seed dispersal, allele-frequency drift,
demography, tissue chemistry, cytotype and direct adaptive performance.
"""
from dataclasses import dataclass
from math import exp, isfinite
import numpy as np
from scipy.linalg import expm

STATES = ("C1_F1", "C1_F0", "C0_F1", "C0_F0")
PIGMENTED_INDEX = 0


@dataclass(frozen=True)
class Rates:
    structural_degradation: float
    structural_reconstruction: float
    floral_switch_off: float
    floral_switch_on: float
    nonfloral_constraint_strength: float = 0.0

    def __post_init__(self):
        for field in (self.structural_degradation,self.structural_reconstruction,
                      self.floral_switch_off,self.floral_switch_on,
                      self.nonfloral_constraint_strength):
            if not isfinite(field) or field < 0:
                raise ValueError("All continuous-time rates and constraint strengths must be finite and nonnegative")


def effective_structural_decay(rates, nonfloral_use):
    """Predictive hypothesis: maintaining a nonfloral role protects coding toolkit."""
    if not isfinite(nonfloral_use) or not 0 <= nonfloral_use <= 1:
        raise ValueError("Nonfloral-use indicator must be in [0,1]")
    return rates.structural_degradation * exp(-rates.nonfloral_constraint_strength*nonfloral_use)


def generator(rates, nonfloral_use=0.0):
    """Continuous-time row generator, valid for either retained or lost coding."""
    decay = effective_structural_decay(rates,nonfloral_use)
    Q=np.zeros((4,4),dtype=float)
    def add(i,j,rate):
        Q[i,j] += rate
        Q[i,i] -= rate
    # coding loss or (hypothetical de novo) rebuilding, F held constant
    add(0,2,decay)
    add(1,3,decay)
    add(2,0,rates.structural_reconstruction)
    add(3,1,rates.structural_reconstruction)
    # flower regulatory activation/off, C held constant
    add(0,1,rates.floral_switch_off)
    add(2,3,rates.floral_switch_off)
    add(1,0,rates.floral_switch_on)
    add(3,2,rates.floral_switch_on)
    assert np.allclose(Q.sum(axis=1),0)
    return Q


def evolved_distribution(rates, nonfloral_use, initial_state, time):
    """Distribution over 4 modules after time, not a genealogy or phenotype rate."""
    if initial_state not in STATES or not isfinite(time) or time<0:
        raise ValueError("Known state and nonnegative time required")
    row=expm(generator(rates,nonfloral_use)*time)[STATES.index(initial_state)]
    return {state:float(max(0,value)) for state,value in zip(STATES,row)}


def observed_flower_distribution(modular_state_distribution):
    """Visible-colour observation collapses THREE molecularly distinct white states."""
    return {
        "pigmented":float(modular_state_distribution["C1_F1"]),
        "white":float(sum(v for k,v in modular_state_distribution.items() if k != "C1_F1"))
    }


def probability_first_reexpression(rates, nonfloral_use, start_white_state, horizon):
    """Probability of ever attaining C1_F1 within a horizon from a white state.

    This is a first-hitting probability under this exact four-state CTMC.
    An event reached via C0->C1 may reflect reconstruction, not intrinsic
    regulatory re-expression. Hence *do not* interpret this result as the
    frequency of one uniquely identified 'ancestral re-expression' mechanism.
    """
    if start_white_state not in STATES[1:] or not isfinite(horizon) or horizon<0:
        raise ValueError("Start in a recognized white state; horizon >=0")
    Qa=generator(rates,nonfloral_use)
    Qa[0,:]=0.0   # pigmented becomes an absorbing first-arrival state
    result=float(expm(Qa*horizon)[STATES.index(start_white_state),0])
    return max(0.0,min(1.0,result))


def closed_form_retention_then_regulatory_reactivation(
        time_since_loss, horizon, n_founders, rates, nonfloral_use):
    """Special case: all founders start C1_F0; no flowering reactivation
    during delay; then single reactivation hazard competes with coding decay.

    Structural reconstruction and flower-regulatory shutdown are explicitly
    absent from this special case, as are gene flow and population death.
    Returns probability at least one of n original founders makes a pigment
    reactivation event in the specified horizon. NOT estimated from plants.
    """
    for value in (time_since_loss,horizon):
        if not isfinite(value) or value<0:
            raise ValueError("Time/horizon must be nonnegative")
    if type(n_founders) is not int or n_founders<1:
        raise ValueError("Number of independent founders must be positive integer")
    if rates.structural_reconstruction!=0:
        raise ValueError("Simple form forbids structural rebuilding")
    decay=effective_structural_decay(rates,nonfloral_use)
    r=rates.floral_switch_on
    opportunity = (
        exp(-decay*time_since_loss) *
        (r/(r+decay))*(1-exp(-(r+decay)*horizon))
        if r+decay>0 else 0.0
    )
    return float(1-(1-opportunity)**n_founders)


def synthetic_organ_protection_grid():
    """Hypothesis demonstration, not biologically calibrated rates."""
    rates=Rates(
        structural_degradation=.12,
        structural_reconstruction=0,
        floral_switch_off=.1,
        floral_switch_on=.05,
        nonfloral_constraint_strength=np.log(6.0)
    )
    rows=[]
    for age in (0,5,10,20,40):
        row={"time_since_floral_loss":age}
        for nonfloral_use in (0.0,1.0):
            tag="nonfloral_use" if nonfloral_use else "no_nonfloral_use"
            row[tag]=closed_form_retention_then_regulatory_reactivation(
                age,10,1,rates,nonfloral_use)
        rows.append(row)
    return {
        "synthetic":True,
        "not_Cirsium_rates":True,
        "parameter_values":{
            "coding_decay_without_nonfloral_use":effective_structural_decay(rates,0),
            "coding_decay_with_nonfloral_use":effective_structural_decay(rates,1),
            "reactivation_rate_per_unit_time":rates.floral_switch_on,
            "nonfloral_constraint":rates.nonfloral_constraint_strength,
            "founders":1,"horizon":10,
            "regulatory_reactivation_disabled_during_loss_age":True,
            "structural_reconstruction":0
        },
        "scenario_comparison":rows,
        "interpretation":"Retained nonfloral usage protects the structural pathway only if beta>0, but flower-specific regulatory competence may still decay. This is an explicit hypothesis, not a property inferred from organ-specific qPCR."
    }


if __name__=="__main__":
    import json
    print(json.dumps(synthetic_organ_protection_grid(),indent=2))
