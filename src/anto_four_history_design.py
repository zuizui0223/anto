"""Exact synthetic assay-selection study for four floral history mechanisms.

The hardest anto distinction is NOT merely pigment vs white, but:
 (1) expression of an already retained ancestral pathway (reuse),
 (2) construction by a new local causal mutation (de novo),
 (3) arrival of an expressed plant via seed migration (immigrant),
 (4) restoration via donor pollen/haplotype introgression.

All 'assays' below are *idealized hypothetical* proxies. Notably, inferring that
an allele predates a split is an uncertain genealogical inference, not a
single experimentally measurable yes/no marker. This code therefore describes
necessary diagnostic information conditional on strong identification assumptions,
not a method to categorize natural samples from a single SNP table.
"""
from itertools import product
import math
import numpy as np

HISTORIES = ('retained_reexpression','de_novo_reconstruction',
             'coloured_immigration','introgressive_reacquisition')
MARKERS = ('foreign_background_ancestry', 'foreign_colour_locus',
           'independent_ancestral_reuse_evidence')
# hypothetical, not measured: most profiles rely on independent history evidence
TRUTH = {
 'retained_reexpression': (0,0,1),
 'de_novo_reconstruction': (0,0,0),
 'coloured_immigration': (1,1,0),
 'introgressive_reacquisition': (0,1,0),
}


def experiment(markers=(), sensitivity=.9, specificity=.95,
               prior=None, collapse_reuse_with_denovo=False):
    for p in (sensitivity,specificity):
        if not math.isfinite(p) or not 0 <= p <= 1:
            raise ValueError('Assay accuracy must be within [0,1]')
    markers=tuple(markers)
    if len(set(markers)) != len(markers) or any(m not in MARKERS for m in markers):
        raise ValueError('Unique registered markers only')
    prior = np.full(4,.25) if prior is None else np.asarray(prior,float)
    if prior.shape != (4,) or np.any(~np.isfinite(prior)) or np.any(prior<=0) or not np.isclose(prior.sum(),1):
        raise ValueError('Four strictly positive prior probabilities required')
    profiles = dict(TRUTH)
    if collapse_reuse_with_denovo:
        # Shared ancestral allele information is not diagnostic (ILS/orthology etc.)
        profiles['retained_reexpression']=profiles['de_novo_reconstruction']
    y=list(product((0,1),repeat=len(markers)))
    kernel=np.ones((len(y),4))
    for i,observed in enumerate(y):
        for j,h in enumerate(HISTORIES):
            for m,val in zip(markers,observed):
                q = sensitivity if profiles[h][MARKERS.index(m)] else 1-specificity
                kernel[i,j] *= q if val else (1-q)
    if not np.allclose(kernel.sum(axis=0),1):
        raise ArithmeticError('Invalid observation probabilities')
    joint=kernel*prior[None,:]
    bayes_accuracy=float(joint.max(axis=1).sum())
    rank=int(np.linalg.matrix_rank(kernel,tol=1e-10))
    return {
        'markers':list(markers), 'hypotheses':list(HISTORIES),
        'rank_known_emissions':rank,
        'all_four_mixture_weights_identifiable_if_emissions_known':rank==4,
        'bayes_accuracy_known_calibrated_assays':bayes_accuracy,
        'naive_chance_accuracy':float(prior.max()),
        'reusability_assay_assumption':'independent and genealogically validated, not directly observable',
        'n_observation_patterns':len(y),
        'all_results_synthetic':True,
    }


def comparison():
    two=MARKERS[:2]
    all_three=MARKERS
    return {
        'visible_colour_only':experiment(),
        'two_ancestry_markers':experiment(two),
        'plus_ancestral_reuse_evidence':experiment(all_three),
        'reuse_evidence_uninformative':experiment(all_three,collapse_reuse_with_denovo=True),
        'three_weak_markers':experiment(all_three,sensitivity=.65,specificity=.65),
    }

if __name__=='__main__':
    import json
    print(json.dumps(comparison(),ensure_ascii=False,indent=2))
