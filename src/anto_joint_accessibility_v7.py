"""Exact bounds on composite floral-trait accessibility under hidden-module linkage.

Aim: test when identical marginal and pairwise retention frequencies are
INSUFFICIENT to predict whether all 3 required developmental resources
coexist in the same propagule. This is a general probabilistic counterexample,
not a new probability theorem and not empirical angiosperm evidence.

A "module" is an independently specified *necessary genetic/developmental
resource*, NOT automatically a locus or a visible head phenotype. Co-presence
is an upper bound on accessibility unless compatibility, floral activation,
selection and establishment are also modeled.

For X1,X2,X3 binary, the sharp family with
    P(Xi=1)=1/2; P(Xi=Xj=1)=1/4 for all i!=j
has P(X1=X2=X3=1)=theta, any theta in [0,1/4].
Even when ALL marginals and pairwise measurements are identical,
three-module accessibility can range from none to 25 percent.

A decorrelation operation is a distributional mixing counterfactual,
NOT a calibrated Mendelian recombination or a claim that gene exchange
breaks physical linkage at the assumed rate.
"""
from fractions import Fraction
from itertools import combinations, product
from math import exp, isfinite, prod, log
from typing import Mapping

MODULES=("M1","M2","M3")
MASKS=tuple("".join(bits) for bits in product("01",repeat=3))
ZERO=Fraction(0,1)
ONE=Fraction(1,1)


def validate_joint(p: Mapping[str,Fraction]):
    if set(p)!=set(MASKS):
        raise ValueError("Exactly eight 3-module inventory states required")
    distribution={s:Fraction(p[s]) for s in MASKS}
    if any(x < 0 for x in distribution.values()) or sum(distribution.values())!=ONE:
        raise ValueError("Probabilities must be nonnegative and sum to 1 exactly")
    return distribution


def pairwise_constrained_family(theta:Fraction):
    """Construct every feasible 3-module exchangeable parity-like family.

    P(000) = P(011) = P(101) = P(110) = 1/4 - theta;
    P(001) = P(010) = P(100) = P(111) = theta.

    All such distributions have P(Mi)=1/2 and P(Mi Mj)=1/4.
    theta = 0: even parity (no complete 111)
    theta = 1/8: independent Bernoulli(1/2) for each module
    theta = 1/4: odd parity (111 occurs 1/4 of the time).
    """
    theta=Fraction(theta)
    if not ZERO <= theta <= Fraction(1,4):
        raise ValueError("Feasible theta is between 0 and 1/4")
    return validate_joint({
        mask:(theta if (mask.count("1")%2) else Fraction(1,4)-theta)
        for mask in MASKS
    })


def co_retained_distribution():
    """Positive complete linkage: 000 and 111, each probability 1/2."""
    return validate_joint({mask:Fraction(1,2) if mask in ("000","111") else ZERO
                           for mask in MASKS})


def marginals(p):
    p=validate_joint(p)
    return tuple(sum(v for mask,v in p.items() if mask[i]=="1")
                 for i in range(3))


def pairwise_joint(p):
    p=validate_joint(p)
    return {(i,j):sum(v for mask,v in p.items()
                       if mask[i]==mask[j]=="1")
            for i,j in combinations(range(3),2)}


def complete_retention_probability(p, required="111"):
    """Requires all specified 1 bits; unspecified bits can be either state."""
    if len(required)!=3 or any(c not in "01" for c in required) or "1" not in required:
        raise ValueError("A nonempty required resource mask is necessary")
    p=validate_joint(p)
    return sum(v for state,v in p.items()
               if all(state[i]=="1" or required[i]=="0" for i in range(3)))


def third_order_central_moment(p):
    """Exactly E[(X1-mean1)(X2-mean2)(X3-mean3)]."""
    p=validate_joint(p)
    mu=marginals(p)
    return sum(v*prod(Fraction(int(mask[i]))-mu[i] for i in range(3))
               for mask,v in p.items())


def expected_expressible_subsets_upper_bound(p):
    """E[2^(number of retained modules)], assuming perfect free ON/OFF control.

    This is at best an upper bound; developmental incompatibility, dosage,
    regulatory epistasis or fitness filtering could collapse real expression.
    """
    p=validate_joint(p)
    return sum(v*(2**mask.count("1")) for mask,v in p.items())


def independent_marginal_distribution(p):
    """Product of component marginal retention probabilities."""
    mu=marginals(p)
    return validate_joint({
        mask:prod(mu[i] if mask[i]=="1" else ONE-mu[i] for i in range(3))
        for mask in MASKS
    })


def decorrelate(p,rho=Fraction(1,2)):
    """Mixture with independent-marginal null, *not* literal gene recombination."""
    rho=Fraction(rho)
    if not ZERO<=rho<=ONE:
        raise ValueError("Mixing fraction must be 0..1")
    p=validate_joint(p)
    baseline=independent_marginal_distribution(p)
    return validate_joint({mask:(ONE-rho)*p[mask]+rho*baseline[mask]
                           for mask in MASKS})


def observed_modules_report(p):
    p=validate_joint(p)
    marg=marginals(p)
    pair=pairwise_joint(p)
    joint=complete_retention_probability(p)
    return {
        "marginals":[str(v) for v in marg],
        "pairwise_joint":{"M1_M2":str(pair[0,1]),"M1_M3":str(pair[0,2]),
                          "M2_M3":str(pair[1,2])},
        "all_three_joint":str(joint),
        "third_order_interaction":str(third_order_central_moment(p)),
        "expected_potential_expression_combinations_upper_bound":
            str(expected_expressible_subsets_upper_bound(p)),
        "n_support_genotypes":sum(v>0 for v in p.values()),
        "probability_mass":{mask:str(v) for mask,v in p.items() if v>0}
    }


def _check_decay(rates, time):
    if len(rates)!=3:
        raise ValueError("Specify degradation hazard for each of three modules")
    if not isfinite(time) or time<0 or any(not isfinite(x) or x<0 for x in rates):
        raise ValueError("All hazards and durations must be finite nonnegative")


def probability_required_modules_after_decay(p, rates, generations,
                                             required="111"):
    """Probability required subset is retained after independent losses.

    This function assumes that module loss is independent conditional on
    initial genotype, exponential per-module hazard, no repair/introgression,
    and no post-dispersal reshuffling of founder lineages.
    """
    _check_decay(rates,generations)
    baseline=complete_retention_probability(p,required)
    survive=exp(-generations*sum(rates[i] for i in range(3)
                                 if required[i]=="1"))
    return float(baseline)*survive


def founders_at_least_one_complete(p_complete, n_founders,
                                  assembly_success=1.0,
                                  establishment_success=1.0):
    """Chance at least one independent founder produces a realized trait.

    A complete inventory is only a *necessary* resource; assembly and
    establishment are stipulated independent conditional probabilities.
    Counts are biological founders, NOT images, transcript isoforms or cells.
    """
    if type(n_founders) is not int or n_founders<1:
        raise ValueError("Positive integer independent founder count required")
    for x in (p_complete, assembly_success, establishment_success):
        if not isfinite(x) or not 0<=x<=1:
            raise ValueError("All probabilities must be in [0,1]")
    q=p_complete*assembly_success*establishment_success
    return -(__import__("math").expm1(n_founders*__import__("math").log1p(-q))) if q<1 else 1.0


def migration_order_counterfactual(p, rate_source, rate_target,
                                   time_exposed, n_founders=1,
                                   assembly_success=1.0,
                                   establishment_success=1.0):
    """Early = spend T in target; late = spend T in source then migrate.

    Both share identical initial source genotype distribution and subsequent
    assembly/establishment conditional probabilities. True seed exchange,
    selective colonization, pollinator changes, and recombination are absent.

    The relative advantage scales as exp(-T*sum(delta_target-delta_source))
    for per-propagule *complete-resource* probability, when source p111>0;
    this is NOT a universal positive absolute difference in accessibility.
    """
    _check_decay(rate_source,time_exposed)
    _check_decay(rate_target,time_exposed)
    early=probability_required_modules_after_decay(p,rate_target,time_exposed)
    late=probability_required_modules_after_decay(p,rate_source,time_exposed)
    ratio=early/late if late>0 else None
    return {
        "complete_resource_early":early,"complete_resource_late":late,
        "ratio_early_over_late":ratio,
        "realized_early":founders_at_least_one_complete(
            early,n_founders,assembly_success,establishment_success),
        "realized_late":founders_at_least_one_complete(
            late,n_founders,assembly_success,establishment_success),
        "rate_contrast_sum":sum(t-s for s,t in zip(rate_source,rate_target)),
        "assumption":"independent genotype loss, identical founder prior, no post-arrival reconstruction, no demographic selection"
    }


def exact_report():
    labels={
        "even_parity":pairwise_constrained_family(Fraction(0)),
        "independent":pairwise_constrained_family(Fraction(1,8)),
        "odd_parity":pairwise_constrained_family(Fraction(1,4)),
        "positive_linkage":co_retained_distribution()
    }
    reports={}
    for name,p in labels.items():
        exact=observed_modules_report(p)
        after=decorrelate(p,Fraction(1,2))
        exact["half_decorrelated_complete"]=str(complete_retention_probability(after))
        exact["half_decorrelated_expected_subsets_upper_bound"]=str(
            expected_expressible_subsets_upper_bound(after))
        reports[name]=exact

    # Independent module rates; hypothetical time units, no empirical plants.
    source=[.02,.02,.02]
    island=[.12,.12,.12]
    exemplar=labels["independent"]
    travel={
        str(k):{
            "sum_source_hazard":round(.02*k,6),
            "sum_island_hazard":round(.12*k,6),
            "relative_early_over_late":exp(-10*.10*k),
            "earlier_migration_resource_probability":
                (Fraction(1,2)**k).__float__()*exp(-10*.12*k),
            "later_migration_resource_probability":
                (Fraction(1,2)**k).__float__()*exp(-10*.02*k)
        } for k in (1,2,3)
    }
    return {
        "project":"anto","version":"v7",
        "empirical_plant_data_admitted":False,
        "claim":"Identical module marginals AND pairwise joint probabilities do not identify the availability of a 3-component target.",
        "exact_sharp_bound_p111_given_marginals_half_and_pairs_quarter":
            ["0","1/4"],
        "distributions":reports,
        "migration_example":{
            "duration":10,"source_rate_per_module":.02,"target_rate_per_module":.12,
            "source_to_target_hypothesis":migration_order_counterfactual(
                exemplar,source,island,10,n_founders=8,
                assembly_success=.7,establishment_success=.5),
            "complexity_grid":travel
        },
        "interpretation_boundary":"Necessary joint resource availability only; realized trait requires measured regulatory/developmental compatibility and ecological establishment. This is a constructed probabilistic demonstration, not novel probability theory or universal biological evidence.",
    }


if __name__=="__main__":
    import json
    print(json.dumps(exact_report(),indent=2,ensure_ascii=False))
