"""Exact observable-equivalence witnesses for lost-floral-trait induction.

Two biologically different latent worlds can have exactly the same first
induction outcome distribution, yet different repeated-expression forecasts.
No novel theorem about hidden Markov chains is claimed. The goal is an
empirical measurement design: measure functional inventory and re-expression
*continuity* rather than counting one colour change as restored evolvability.

Every first-pulse trial is an independent lineage Bernoulli observation in this
toy model. Population kinship, flowering phases, recurrence dependence and
fitness are NOT represented or calibrated in real plants.
"""
from fractions import Fraction as F
from math import comb
from anto_afterlife_frontier_v10 import Rates, run
from anto_joint_accessibility_v7 import pairwise_constrained_family, marginals, pairwise_joint

ZERO = F(0)
ONE = F(1)


def _fraction(value):
    p = F(str(value))
    if not ZERO <= p <= ONE:
        raise ValueError("Expected a finite [0, 1] probability")
    return p


def bayes_error_equal_prior_binomial(p, q, n_independent_lineages: int):
    """Exact minimum equal-prior 0/1 misclassification risk of two Bin(n,p) laws.

    Standard total-variation identity, not a new result. The effective sample
    size is independent *lineages*, not plant images, flowers, or markers.
    """
    p, q = _fraction(p), _fraction(q)
    if type(n_independent_lineages) is not int or n_independent_lineages < 1:
        raise ValueError("Need positive integer independent biological lineages")
    n = n_independent_lineages
    overlap = sum((
        min(F(comb(n,k))*p**k*(ONE-p)**(n-k),
            F(comb(n,k))*q**k*(ONE-q)**(n-k))
        for k in range(n+1)
    ), ZERO)
    return overlap/2


def first_induction_counterexample(horizon=12):
    """Equal complete first-pulse outcomes, sharply different recurrent events."""
    if type(horizon) is not int or horizon < 3:
        raise ValueError("Three or more generations needed")
    A = {"theta":F(1,4), "rates": Rates(F(2,5), F(1,10), F(1,10))}
    B = {"theta":F(1,8), "rates": Rates(F(4,5), F(4,5), F(1,10))}
    pa = run(A["theta"], A["rates"], horizon)
    pb = run(B["theta"], B["rates"], horizon)
    if pa[0]["present_expression"] != pb[0]["present_expression"]:
        raise ArithmeticError("Initial visible states unexpectedly differ")
    if pa[1]["present_expression"] != pb[1]["present_expression"]:
        raise ArithmeticError("The deliberately matched first-induction probabilities differ")
    if pa[1]["first_expression_ever"] != pb[1]["first_expression_ever"]:
        raise ArithmeticError("The first-induction observation distributions differ")
    if pa[-1]["three_consecutive_expression_ever"] == pb[-1]["three_consecutive_expression_ever"]:
        raise ArithmeticError("The later event should separate the worlds")
    p, q = (pa[-1]["three_consecutive_expression_ever"],
            pb[-1]["three_consecutive_expression_ever"])
    first = pa[1]["present_expression"]
    return {
        "first_pulse_visible_positive_both": str(first),
        "first_pulse_bayes_error_any_n": "1/2",
        "initial_joint_functional_toolkit_A": str(A["theta"]),
        "initial_joint_functional_toolkit_B": str(B["theta"]),
        "activation_rate_A": str(A["rates"].activation),
        "activation_rate_B": str(B["rates"].activation),
        "continuation_rate_A": str(A["rates"].repeat),
        "continuation_rate_B": str(B["rates"].repeat),
        "same_resource_loss_rate": str(A["rates"].degradation),
        "three_in_a_row_by_horizon_A": str(p),
        "three_in_a_row_by_horizon_B": str(q),
        "recurrent_event_gap_B_minus_A": str(q-p),
        "equal_prior_model_discrimination_error": {
            str(n): float(bayes_error_equal_prior_binomial(p,q,n))
            for n in (10,20,40,80)
        },
        "not_real_plant_sample_power": True,
    }


def known_initial_toolkit_still_not_sufficient(horizon=12):
    """Even an exact resource census + one induction can leave r unidentified."""
    theta = F(1,8)
    A = run(theta, Rates(F(2,5),F(1,10),F(1,10)),horizon)
    B = run(theta, Rates(F(2,5),F(4,5),F(1,10)),horizon)
    assert A[1]["present_expression"] == B[1]["present_expression"]
    return {
        "same_full_initial_toolkit_fraction":str(theta),
        "same_first_generation_expression":str(A[1]["present_expression"]),
        "different_three_streak_after_horizon":{
            "low_continuation":str(A[-1]["three_consecutive_expression_ever"]),
            "high_continuation":str(B[-1]["three_consecutive_expression_ever"])
        },
        "inference_ceiling":"Perfect baseline resource measurement and one induction do not identify the later continuation rate in the declared model.",
    }


def equal_pairwise_hidden_capacity(history_steps=8):
    """Identical visible white trajectory and pairwise modules, different rescue.

    Preintervention a=0,d=0. At final step apply identical induction a=1/2.
    In the two latent inventory laws all single and pairwise module
    co-retention moments match, but P(111) differs by exactly 1/4.
    """
    if type(history_steps) is not int or history_steps < 1:
        raise ValueError("Need positive number of uninduced steps")
    even = pairwise_constrained_family(F(0))
    odd = pairwise_constrained_family(F(1,4))
    if marginals(even)!=marginals(odd) or pairwise_joint(even)!=pairwise_joint(odd):
        raise ArithmeticError("Pairwise moment parity witness invalid")
    no_induction = Rates(0,F(1,2),0)
    induction = Rates(F(1,2),F(1,2),0)
    pre = [no_induction]*history_steps
    a = run(F(0),pre+[induction],history_steps+1)
    b = run(F(1,4),pre+[induction],history_steps+1)
    for x,y in zip(a[:-1],b[:-1]):
        if x["present_expression"] or y["present_expression"]:
            raise ArithmeticError("Non-intervention baseline was not all-white")
    return {
        "uninduced_phenotype_trace":"all_white",
        "uninduced_generations":history_steps,
        "first_module_marginals_each":"1/2",
        "pairwise_module_co_retention_each":"1/4",
        "joint_toolkit_A":"0",
        "joint_toolkit_B":"1/4",
        "identical_final_induction_a":"1/2",
        "induced_positive_A":str(a[-1]["present_expression"]),
        "induced_positive_B":str(b[-1]["present_expression"]),
        "not_an_empirical_ancestry_loss_reconstruction":True,
    }


def panel():
    return {
        "schema":"anto_afterlife_induction_identifiability_v10b",
        "conditional_exact_synthetic":True,
        "independent_lineage_counts_are_only_hypothetical":True,
        "first_induction_counterexample":first_induction_counterexample(),
        "measured_toolkit_but_unknown_continuation":known_initial_toolkit_still_not_sufficient(),
        "matched_marginal_pairwise_and_visible_history":equal_pairwise_hidden_capacity(),
        "ecological_fitness_or_true_trait_regain_measured":False,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }
