"""Negative-control experiment: visible-colour equivalence of distinct histories.

These are stylized *synthetic* causal histories (intrinsic re-expression,
coloured-seed immigration, pollen introgression). Their visible endpoint has
exactly the same Bernoulli probability by construction, so colour alone has
zero identifying information about the mechanism. Diagnostic genomic assays
are synthetic noisy measurements, not real validated markers.
"""
import math
import random
from collections import Counter

HISTORIES = ("intrinsic", "immigration", "introgression")
# (nonlocal genome-wide ancestry, nonlocal colour-linked haplotype)
DIAGNOSTIC_TRUTHS = {
    "intrinsic": (0, 0),
    "immigration": (1, 1),
    "introgression": (0, 1),
}


def generate(n_per_class=600, colour_probability=0.25,
             sensitivity=0.9, specificity=0.95, seed=20261017):
    for p in (colour_probability, sensitivity, specificity):
        if not 0 <= p <= 1:
            raise ValueError("Probabilities must lie in [0,1]")
    if type(n_per_class) is not int or n_per_class < 1:
        raise ValueError("n_per_class must be positive integer")
    rng = random.Random(seed)
    records = []
    for history in HISTORIES:
        truth = DIAGNOSTIC_TRUTHS[history]
        for _ in range(n_per_class):
            coloured = int(rng.random() < colour_probability)
            assays = []
            for bit in truth:
                p_positive = sensitivity if bit else (1-specificity)
                assays.append(int(rng.random() < p_positive))
            records.append((history, coloured, tuple(assays)))
    rng.shuffle(records)
    return records


def _fit(train):
    """Estimate molecular assay emission probabilities from labelled training.

    In nature such labels require independent known histories; here oracle
    synthetic labels exist solely for a negative-control benchmark.
    """
    fitted = {}
    for h in HISTORIES:
        values = [assays for history, _, assays in train if history == h]
        if not values:
            raise ValueError("Training requires labelled observations for each history")
        fitted[h] = [(1 + sum(a[j] for a in values)) / (2 + len(values))
                     for j in range(2)]
    return fitted


def _predict(assays, probabilities=None):
    if probabilities is None:
        return HISTORIES[0]   # visible-only tie, all histories equiprobable
    scores = {}
    for h in HISTORIES:
        scores[h] = sum(math.log(p if x else (1-p))
                        for x, p in zip(assays, probabilities[h]))
    return max(HISTORIES, key=lambda h: scores[h])


def benchmark(seed=20261017, train_per_class=300, test_per_class=1500,
              sensitivity=0.9, specificity=0.95):
    data = generate(train_per_class + test_per_class, seed=seed,
                    sensitivity=sensitivity, specificity=specificity)
    training, testing = [], []
    for h in HISTORIES:
        subset = [row for row in data if row[0] == h]
        training.extend(subset[:train_per_class])
        testing.extend(subset[train_per_class:])
    emission = _fit(training)
    predictions = {}
    for method in ("visible_only", "two_molecular_markers"):
        counts = Counter()
        correct = 0
        for h, colour, assays in testing:
            predicted = _predict(assays, None if method == "visible_only" else emission)
            counts[f"{h}->{predicted}"] += 1
            correct += (predicted == h)
        predictions[method] = {
            "accuracy": correct / len(testing),
            "correct": correct, "n": len(testing),
            "confusion": dict(sorted(counts.items())),
        }
    return {"seed": seed, "synthetic": True,
            "mechanism_prevalence": "balanced",
            "visible_endpoint_probability_by_history": {
                h: 0.25 for h in HISTORIES},
            "trained_assay_positive_rates": emission,
            "predictions": predictions,
            "warning": "A proof-of-need for history-specific measurements; not empirical inference."}


if __name__ == "__main__":
    import json
    print(json.dumps(benchmark(), indent=2, ensure_ascii=False))
