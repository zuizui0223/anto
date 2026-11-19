#!/usr/bin/env python3
"""ANTO v47: posthoc 3-original-source-population matching, NO new model refit.

Input is the SHA-locked v46 derived JSON (20 original population-heldout MSEs);
these are not independent studies and no p-value is defined.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data/external/anto_v47_three_pop_posthoc_registry_20261119.json"
PARENT = ROOT / "results/anto_v46_phlox2024_population_heldout_20261118.json"
RESULT = ROOT / "results/anto_v47_three_pop_reference_20261119.json"
MODELS = ("G", "T", "E", "TE")
ENDPOINTS = ("fruits_2019", "flowers_2019")
TAXA = ("amoena", "deamii", "pilosa")


class V47Error(ValueError):
    pass


def require(ok, label):
    if not ok:
        raise V47Error(label)


def git_blob_sha1(raw):
    return hashlib.sha1(("blob %d\0" % len(raw)).encode("ascii") + raw).hexdigest()


def read_registry():
    info = json.loads(REGISTRY.read_text())
    require(
        info["analysis"] == "POSTHOC_ORIGINAL_2024_PHLOX_THREE_SOURCE_POPULATION_MATCHING"
        and info["v46_parent_git_blob_sha1"] ==
            "f03b0eb0177f807317ee1f8a7c0e65ee408d903b"
        and info["source_2024_CSV_blob_sha"] ==
            "4c780b0abe644ac69dae86a0450eb598c160626b"
        and info["source_2024_author_git_sha"] ==
            "d792fe58cce36721394adfd97e9f412396578e86"
        and info["source_2024_rows"] == 963
        and info["source_2024_genotypes"] == 284
        and info["source_2024_gardens"] == 3
        and info["source_2024_original_populations"] == 20
        and info["source_2024_taxon_populations"] ==
            {"amoena": 8, "deamii": 3, "pilosa": 9}
        and info["endpoints"] == list(ENDPOINTS)
        and info["fixed_original_models"] == list(MODELS)
        and info["reference_methods"]["amoena_triplets"] == 56
        and info["reference_methods"]["pilosa_triplets"] == 84
        and info["reference_methods"]["pooled_triplets_primary"] == 140
        and info["reference_methods"]["pooled_triplets_secondary"] == 680
        and info["independent_new_biological_studies"] == 0
        and info["new_2014_Phlox_selection_refits"] == 0
        and info["has_p_value"] is False
        and info["regime"] == "posthoc_descriptive_no_new_fits",
        "V47_REGISTRY_OR_CLAIM_DRIFT",
    )
    return info


def read_parent(registry, path=PARENT):
    raw = Path(path).read_bytes()
    require(git_blob_sha1(raw) == registry["v46_parent_git_blob_sha1"],
            "V47_V46_PARENT_SOURCE_BLOB_CHANGED")
    data = json.loads(raw)
    validate_parent(data, registry)
    return data


def validate_parent(data, registry):
    require(data["analysis_type"] ==
            "ORIGINAL_2024_SOURCE_POPULATION_HELDOUT_PREDICTION"
            and data["source_git_blob_sha1"] ==
            registry["source_2024_CSV_blob_sha"]
            and data["source_commit"] == registry["source_2024_author_git_sha"]
            and data["distinct_original_programs_added"] == 1
            and data["independent_evolutionary_radiations_added"] == 0
            and data["original_2014_Phlox_selection_source_used_for_model_fit"] is False
            and data["measured_floral_color_or_pollinator_or_reproductive_alleles"] is False
            and data["primary_fruit_transfer_pass_models"] == []
            and data["fixed_models"] == ["B0"] + list(MODELS),
            "V47_PARENT_SOURCE_OR_BIOLOGICAL_SCOPE_DRIFT")
    grain = data["source_original_grain"]
    require(
        grain["rows"] == 963 and grain["genotypes"] == 284
        and grain["populations"] == 20
        and grain["populations_by_taxon"] ==
            {"amoena": 8, "deamii": 3, "pilosa": 9}
        and grain["biomass_original_negative_n"] == 13
        and grain["biomass_original_minimum"] == -233,
        "V47_ORIGINAL_SOURCE_DENOMINATORS_CHANGED")
    blocked = data["nonadmitted_secondary_biomass"]
    require(blocked["status"] == "BLOCKED_PREDECLARED_LOG1P_DOMAIN"
            and blocked["numerical_scores_calculated"] is False,
            "V47_TECHNICAL_BIOMASS_NONADMISSION_CHANGED")
    require(set(data["endpoints"]) == set(ENDPOINTS),
            "V47_ENDPOINT_SCOPE_CHANGED")
    for e in ENDPOINTS:
        obj = data["endpoints"][e]
        require(set(obj["model_results"]) == set(("B0",) + MODELS)
                and len(obj["per_population_holdout"]) == 20,
                "V47_HOLDOUT_MODEL_OR_POPULATION_SCOPE_CHANGED")
        rows = obj["per_population_holdout"]
        require(len({r["population"] for r in rows}) == 20
                and sum(r["n_ramets"] for r in rows) == 963
                and Counter(r["taxon"] for r in rows) ==
                    Counter({"amoena": 8, "deamii": 3, "pilosa": 9}),
                "V47_HOLDOUT_ORIGINAL_POPULATION_GRAIN_CHANGED")
        for r in rows:
            require(r["n_ramets"] > 0 and
                    set(r["mse_by_model"]) == set(("B0",) + MODELS),
                    "V47_FOLD_MODEL_OR_COUNTS_DRIFT")
            require(all(isinstance(x, (int, float))
                        and math.isfinite(x) and x > 0
                        for x in r["mse_by_model"].values()),
                    "V47_NONFINITE_OR_INVALID_HOLDOUT_MSE")
        for m in ("B0",) + MODELS:
            mr = obj["model_results"][m]
            # v46 stores per-population prediction losses at 8 decimals;
            # numeric parent summary was generated before that rounding.
            avg = sum(r["mse_by_model"][m] for r in rows) / 20
            weighted = sum(r["n_ramets"] * r["mse_by_model"][m]
                           for r in rows) / 963
            require(math.isclose(avg, mr["population_equal_mse"],
                                 rel_tol=0, abs_tol=2e-7)
                    and math.isclose(weighted,
                                     mr["ramet_record_weighted_mse"],
                                     rel_tol=0, abs_tol=2e-7),
                    "V47_PARENT_LOSS_AGGREGATE_DRIFT")
            positives = sum(
                r["mse_by_model"][m] < r["mse_by_model"]["B0"] - 1e-12
                for r in rows)
            require(positives == mr["positive_population_folds"],
                    "V47_PARENT_SIGN_COUNT_DRIFT")
            if e == "fruits_2019":
                require(mr["predeclared_primary_fruit_decision"] ==
                        "NOT_ESTABLISHED",
                        "V47_PARENT_FAILED_PRIMARY_GATE_CHANGED")


def gain(rows, model):
    require(rows, "V47_EMPTY_SOURCE_BLOCK")
    base = sum(x["mse_by_model"]["B0"] for x in rows)
    other = sum(x["mse_by_model"][model] for x in rows)
    require(base > 0 and math.isfinite(base) and math.isfinite(other),
            "V47_BAD_REAGGREGATION_DENOMINATOR")
    return 1 - other / base


def metric(v):
    require(math.isfinite(v), "V47_NONFINITE_DIAGNOSTIC")
    return round(v, 10)


def summarize(reference, deamii, equal_tol):
    vals = sorted(reference)
    require(vals and 0 <= equal_tol <= 1e-8,
            "V47_INVALID_REFERENCE_OR_TOLERANCE")
    n = len(vals)
    return {
        "n_overlapping_reference_combinations": n,
        "number_negative": sum(v < 0 for v in vals),
        "number_at_or_below_deamii": sum(v <= deamii + equal_tol for v in vals),
        "minimum": metric(vals[0]),
        "median": metric((vals[(n - 1) // 2] + vals[n // 2]) / 2),
        "maximum": metric(vals[-1]),
        "is_independent_study_replication": False,
        "is_valid_randomization_p_value": False,
    }


def analyze(source, registry):
    validate_parent(source, registry)
    output = {
        "analysis_type": "POSTHOC_FROZEN_V46_SOURCE_POPULATION_COUNT_MATCHING",
        "original_input_v46_git_blob_sha1": registry["v46_parent_git_blob_sha1"],
        "original_CSV_source_git_blob_sha1": registry["source_2024_CSV_blob_sha"],
        "original_source_studies_added_in_v47": 0,
        "new_prediction_fits": 0,
        "new_original_field_populations": 0,
        "primary_v46_all_taxa_transfer_pass_models": [],
        "technical_biomass_nonadmission_retained": True,
        "reference_is_no_refit_overlapping_combinatorial_conditioning": True,
        "new_hypothesis_test_or_p_values": False,
        "n_same_taxon_amoena_triplets": 56,
        "n_same_taxon_pilosa_triplets": 84,
        "n_same_taxon_triplets_total": 140,
        "n_secondary_pooled_17_other_taxon_triplets": 680,
        "endpoints": {},
    }
    for e in ENDPOINTS:
        original = source["endpoints"][e]["per_population_holdout"]
        groups = {
            t: sorted((r for r in original if r["taxon"] == t),
                      key=lambda r: r["population"])
            for t in TAXA
        }
        summary = {
            "outcome_role": "PRIMARY_REAGGREGATION" if e == "fruits_2019"
                            else "SECONDARY_DESCRIPTIVE_REAGGREGATION",
            "models": {},
        }
        for m in MODELS:
            d = gain(groups["deamii"], m)
            am = [gain(c, m) for c in
                  itertools.combinations(groups["amoena"], 3)]
            pi = [gain(c, m) for c in
                  itertools.combinations(groups["pilosa"], 3)]
            within = am + pi
            pooled = [gain(c, m) for c in
                      itertools.combinations(
                          groups["amoena"] + groups["pilosa"], 3)]
            require(len(am) == 56 and len(pi) == 84 and
                    len(within) == 140 and len(pooled) == 680,
                    "V47_INCORRECT_EXACT_COMBINATORIAL_COUNTS")
            pairs = [
                {"deleted_deamii_population": drop["population"],
                 "remaining_deamii_population_count": 2,
                 "remaining_relative_MSE_gain": metric(gain(
                     [x for x in groups["deamii"]
                      if x["population"] != drop["population"]], m))}
                for drop in groups["deamii"]
            ]
            summary["models"][m] = {
                "deamii_3_original_population_ids": [
                    r["population"] for r in groups["deamii"]],
                "deamii_3_pop_relative_MSE_gain": metric(d),
                "deamii_positive_original_population_folds": sum(
                    r["mse_by_model"][m] <
                    r["mse_by_model"]["B0"] - 1e-12
                    for r in groups["deamii"]),
                "deamii_two_pop_delete_one_arithmetic": pairs,
                "amoena_same_taxon_choose3_reference": summarize(
                    am, d, registry["comparison_tolerance"]),
                "pilosa_same_taxon_choose3_reference": summarize(
                    pi, d, registry["comparison_tolerance"]),
                "primary_same_taxon_choose3_140_reference": summarize(
                    within, d, registry["comparison_tolerance"]),
                "secondary_mixed_taxon_choose3_680_reference": summarize(
                    pooled, d, registry["comparison_tolerance"]),
            }
        output["endpoints"][e] = summary
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    registry = read_registry()
    value = analyze(read_parent(registry), registry)
    if args.check:
        require(value == json.loads(RESULT.read_text()),
                "V47_FROZEN_RECEIPT_REGENERATION_MISMATCH")
        print("V47_SOURCE_BLOB_AND_ALL_140_PLUS_680_COMBINATIONS_EXACT_PASS")
    else:
        print(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
