"""Exact observation-design identifiability for synthetic floral re-expression histories.

This is a DESIGN tool, not empirical assignment of real evolutionary histories.

Three causal hypotheses can produce identical observed flower colour:
    intrinsic      W+ locally re-expresses; no foreign coloured haplotype
    immigration    already-coloured plant establishes; foreign background and haplotype
    introgression  colour haplotype crosses into local maternal background

Genomic assay proxies are imperfect, and the assumed truth profiles may fail
under incomplete lineage sorting, shared ancestral variation, and backcrossing.
Under KNOWN emission probabilities, linear independence of observation columns
is necessary and sufficient to identify mixture weights. It does not identify
historical paths of specific individuals without additional data/assumptions.
"""
from itertools import product
from math import log
import numpy as np

HISTORIES = ('intrinsic', 'immigration', 'introgression')
PROFILES = {
    'intrinsic': {'foreign_genome': 0, 'foreign_pigment_haplotype': 0},
    'immigration': {'foreign_genome': 1, 'foreign_pigment_haplotype': 1},
    'introgression': {'foreign_genome': 0, 'foreign_pigment_haplotype': 1},
}


def _validate_probability(x):
    if not np.isfinite(x) or not 0 <= x <= 1:
        raise ValueError('probabilities must lie between 0 and 1')


def emission_matrix(markers=(), sensitivity=.90, specificity=.95,
                    missing_probability=0.0, profiles=None):
    """Return (observation outcomes, emission matrix Y x H).

    Marker observation outcomes: 0, 1, or -1 for missing. Missing-at-random
    is assumed, independent of history and marker status. Conditional assay
    independence is also assumed; correlated observations require a joint
    emission model not implemented here.

    profiles can specify different TRUE marker patterns by history to model
    failure of proposed lineage/haplotype diagnostics.
    """
    profiles = PROFILES if profiles is None else profiles
    for p in (sensitivity, specificity, missing_probability):
        _validate_probability(p)
    markers = tuple(markers)
    if len(set(markers)) != len(markers):
        raise ValueError('duplicate assay markers')
    if any(marker not in PROFILES[HISTORIES[0]] for marker in markers):
        raise ValueError('unknown marker')
    for h in HISTORIES:
        if h not in profiles or any(profiles[h].get(m) not in (0, 1) for m in markers):
            raise ValueError('every history needs a binary true marker profile')
    alphabet = (-1, 0, 1) if missing_probability > 0 else (0, 1)
    outcomes = list(product(alphabet, repeat=len(markers)))
    mat = np.zeros((len(outcomes), len(HISTORIES)), dtype=float)
    for i, obs in enumerate(outcomes):
        for j, history in enumerate(HISTORIES):
            prob = 1.0
            for name, value in zip(markers, obs):
                if value == -1:
                    prob *= missing_probability
                else:
                    true = profiles[history][name]
                    q = sensitivity if true else 1-specificity
                    prob *= (1-missing_probability)*(q if value == 1 else 1-q)
            mat[i, j] = prob
    if not np.allclose(mat.sum(axis=0), 1.0, atol=1e-12):
        raise ArithmeticError('emissions must sum to 1 within each history')
    return outcomes, mat


def diagnose(markers=(), sensitivity=.9, specificity=.95,
             missing_probability=0.0, priors=(1/3, 1/3, 1/3), profiles=None):
    """Exact diagnostic capacity under specified known likelihoods.

    rank 3 -> mixture-weight-identifiable WITH KNOWN emission functions.
    expected_optimal_accuracy is the Bayes upper benchmark from hypothetically
    validated history-dependent measurements, not a fitted classifier result.
    mutual_information_bits quantifies expected history information per sample.
    """
    priors = np.asarray(priors, dtype=float)
    if priors.shape != (3,) or np.any(~np.isfinite(priors)) or np.any(priors <= 0) or not np.isclose(priors.sum(), 1):
        raise ValueError('positive priors for 3 histories must sum to 1')
    outcomes, mat = emission_matrix(markers, sensitivity, specificity,
                                    missing_probability, profiles=profiles)
    singular = np.linalg.svd(mat, compute_uv=False)
    rank = int(np.linalg.matrix_rank(mat, tol=1e-10))
    weighted = mat * priors[None, :]
    event_probabilities = weighted.sum(axis=1)
    optimal_accuracy = float(weighted.max(axis=1).sum())
    mi = 0.0
    for i in range(len(outcomes)):
        for j in range(len(HISTORIES)):
            if weighted[i,j] > 0 and event_probabilities[i] > 0:
                mi += weighted[i,j] * log(mat[i,j] / event_probabilities[i], 2)
    return {
        'synthetic_design': True,
        'assays': list(markers),
        'rank_known_emissions': rank,
        'mixture_weights_identifiable_if_emissions_known': rank == 3,
        'singular_values': [float(v) for v in singular],
        'expected_bayes_accuracy': optimal_accuracy,
        'chance_accuracy': float(priors.max()),
        'mutual_information_bits': float(mi),
        'n_possible_observations': len(outcomes),
        'assumed_independent_assays': True,
        'assumed_missing_at_random': True,
    }


def run_design_grid():
    """Predefined scenario comparison; never infer a field marker from this."""
    both = ('foreign_genome', 'foreign_pigment_haplotype')
    collapse = {
        'intrinsic': {'foreign_genome': 0, 'foreign_pigment_haplotype': 0},
        'immigration': {'foreign_genome': 1, 'foreign_pigment_haplotype': 1},
        'introgression': {'foreign_genome': 0, 'foreign_pigment_haplotype': 0},
    }
    scenarios = [
        ('colour_only', (), .90, .95, 0, None),
        ('genome_background_only', both[:1], .90, .95, 0, None),
        ('pigment_haplotype_only', both[1:], .90, .95, 0, None),
        ('two_reliable_diagnostics', both, .90, .95, 0, None),
        ('two_noisy_diagnostics', both, .65, .65, 0, None),
        ('two_uninformative_diagnostics', both, .50, .50, 0, None),
        ('two_reliable_but_half_missing', both, .90, .95, .50, None),
        ('shared_ancestral_haplotype', both, .90, .95, 0, collapse),
    ]
    return {name: diagnose(markers, sens, spec, missing, profiles=profiles)
            for name, markers, sens, spec, missing, profiles in scenarios}


if __name__ == '__main__':
    import json
    print(json.dumps(run_design_grid(), ensure_ascii=False, indent=2))
