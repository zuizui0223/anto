#!/usr/bin/env python3
"""ANTO v49: original 2024 Phlox focal-2/pop heldout using v46 model unchanged.

POSTHOC diagnostic after v46-v48. Biological source is the same 2024
experiment; source-population panels overlap and are NOT independent trials.
No author raw observations are redistributed in the receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import statistics
from collections import Counter
from pathlib import Path

import analyze_original_Phlox_2024_population_holdout_v46 as v46

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "data/external/anto_v49_multi_taxon_focal2_registry_20261119.json"
PARENT = ROOT / "results/anto_v46_phlox2024_population_heldout_20261118.json"
RECEIPT = ROOT / "results/anto_v49_original_Phlox_multitaxon_focal2_20261119.json"
TAXA = ("amoena", "deamii", "pilosa")
ENDPOINTS = ("fruits_2019", "flowers_2019")
MODELS = ("B0", "G", "T", "E", "TE")


class V49Error(ValueError):
    pass


def require(ok, code):
    if not ok:
        raise V49Error(code)


def registry():
    r = json.loads(REGISTRY_PATH.read_text())
    require(
        r["analysis_type"] == "POSTHOC_MULTI_TAXON_FOCAL_2_TRAIN_SOURCE_POP_REFIT"
        and r["source_repo"] == "PhloxHopkins/PhloxFieldAdaptation"
        and r["source_git_commit"] ==
            "d792fe58cce36721394adfd97e9f412396578e86"
        and r["source_CSV_path"] == "Dataset_fitness_leaftraits_distance_lat_long.csv"
        and r["source_CSV_git_blob_sha1"] ==
            "4c780b0abe644ac69dae86a0450eb598c160626b"
        and r["source_CSV_bytes"] == 175377
        and r["v46_derived_result_git_blob_sha1"] ==
            "f03b0eb0177f807317ee1f8a7c0e65ee408d903b"
        and r["n_original_ramets"] == 963
        and r["n_original_clonal_genotypes"] == 284
        and r["n_original_source_populations"] == 20
        and r["original_taxon_source_population_counts"] ==
            {"amoena": 8, "deamii": 3, "pilosa": 9}
        and r["original_gardens"] == list(v46.GARDENS)
        and r["models"] == list(MODELS)
        and r["original_v46_ridge_penalty"] == 10
        and r["model_implementation"] ==
            "import v46.fit_train and v46.prediction without modification"
        and r["target_taxon_training_pops"] == 2
        and r["heldout_focal_source_populations"] == 1
        and r["nonfocal_source_pops_training"] ==
            "ALL original nonfocal source populations"
        and r["panels_by_taxon"] == {"amoena": 56, "deamii": 1, "pilosa": 84}
        and r["n_source_triplet_panels"] == 141
        and r["n_holdout_folds"] == 423
        and r["n_reference_panels"] == 140
        and r["primary_endpoint"] == "fruits_2019"
        and r["secondary_endpoint"] == "flowers_2019"
        and r["biomass_technical_nonadmission"] is True
        and r["negative_biomass_n"] == 13
        and r["negative_biomass_min"] == -233
        and r["deamii_v46_mse_replay_tolerance"] == 1e-7
        and r["reference_comparison_tolerance"] == 1e-10
        and r["n_additional_original_studies"] == 0
        and r["flower_pigment_measured"] is False
        and r["pollinator_measured"] is False
        and r["genetic_fitness_measured"] is False
        and r["primary_v46_all_taxon_transfer_established"] is False
        and r["is_posthoc"] is True
        and r["valid_randomization_p"] is False,
        "V49_PROTOCOL_SOURCE_OR_CLAIM_DRIFT",
    )
    return r


def author_source(path, r):
    inherited = v46.source_registry()
    require(inherited["source_commit"] == r["source_git_commit"]
            and inherited["source_git_blob_sha1"] ==
                r["source_CSV_git_blob_sha1"]
            and inherited["source_size"] == r["source_CSV_bytes"],
            "V49_V46_SOURCE_REGISTRY_CHANGED")
    rows, audit = v46.parse_source(v46.acquire_raw(inherited, path), inherited)
    require(audit["rows"] == 963 and audit["genotypes"] == 284
            and audit["populations"] == 20
            and audit["biomass_original_negative_n"] == 13
            and audit["biomass_original_minimum"] == -233,
            "V49_SOURCE_COHORT_CHANGED")
    return rows, audit


def git_blob_sha1(raw):
    return hashlib.sha1(("blob " + str(len(raw)) + "\0").encode("ascii") + raw).hexdigest()


def frozen_v46(r):
    content = PARENT.read_bytes()
    require(git_blob_sha1(content) == r["v46_derived_result_git_blob_sha1"],
            "V49_FROZEN_V46_RECEIPT_CHANGED")
    parent = json.loads(content)
    require(parent["analysis_type"] ==
                "ORIGINAL_2024_SOURCE_POPULATION_HELDOUT_PREDICTION"
            and parent["source_git_blob_sha1"] == r["source_CSV_git_blob_sha1"]
            and parent["primary_fruit_transfer_pass_models"] == []
            and parent["fixed_models"] == list(MODELS)
            and parent["source_original_grain"]["rows"] == 963,
            "V49_PREVIOUS_FAILED_GATE_OR_SOURCE_REINTERPRETED")
    for e in ENDPOINTS:
        require(len(parent["endpoints"][e]["per_population_holdout"]) == 20
                and set(parent["endpoints"][e]["model_results"]) == set(MODELS),
                "V49_INHERITED_V46_FOLDS_CHANGED")
    return parent


def original_populations(rows):
    poptax = {}
    bypop = {}
    clone_owner = {}
    for row in rows:
        pop, taxon, geno = row["population"], row["taxon"], row["genotype"]
        require(pop not in poptax or poptax[pop] == taxon,
                "V49_ORIGINAL_POP_MULTITAXON")
        poptax[pop] = taxon
        bypop.setdefault(pop, []).append(row)
        require(geno not in clone_owner or clone_owner[geno] == pop,
                "V49_ORIGINAL_CLONE_REASSIGNED")
        clone_owner[geno] = pop
    keys = {t: sorted(p for p in poptax if poptax[p] == t) for t in TAXA}
    require({t: len(keys[t]) for t in TAXA} ==
            {"amoena": 8, "deamii": 3, "pilosa": 9}
            and len(poptax) == 20 and len(clone_owner) == 284,
            "V49_ORIGINAL_POPULATION_COUNTS_CHANGED")
    return bypop, poptax, keys


def training_populations(focal, panel, heldout, poptax):
    require(focal in TAXA and len(set(panel)) == 3
            and set(panel).issubset(
                {p for p in poptax if poptax[p] == focal})
            and heldout in panel,
            "V49_INVALID_FOCAL_TRIPLET")
    chosen_focal = [p for p in panel if p != heldout]
    nonfocal = sorted(p for p in poptax if poptax[p] != focal)
    chosen = sorted(chosen_focal + nonfocal)
    require(len(chosen_focal) == 2 and heldout not in chosen
            and set(chosen) == (
                set(chosen_focal) | {p for p in poptax if poptax[p] != focal})
            and set(p for p in chosen if poptax[p] == focal) == set(chosen_focal),
            "V49_NOT_TWO_FOCAL_ALL_NONFOCAL_TRAIN")
    return chosen


def score_one_fold(rows, bypop, poptax, focal, panel, heldout):
    train_pops = training_populations(focal, panel, heldout, poptax)
    chosen_set = set(train_pops)
    # Preserve the original v46 source input row sequence (not population grouping)
    train = [r for r in rows if r["population"] in chosen_set]
    test = [r for r in rows if r["population"] == heldout]
    require(len(test) == len(bypop[heldout]) and test
            and len(train) == sum(len(bypop[p]) for p in train_pops)
            and not ({r["genotype"] for r in train} &
                     {r["genotype"] for r in test})
            and set(r["garden"] for r in test) == set(v46.GARDENS),
            "V49_TRAIN_TEST_SOURCE_OR_GENOTYPE_LEAKAGE")
    require(set(v46.group(r) for r in train) ==
            {(g, t) for g in v46.GARDENS for t in TAXA},
            "V49_ALL_NINE_ORIGINAL_V46_GROUP_CELLS_REQUIRED")
    fit = v46.fit_train(train)
    loss = {}
    for e in ENDPOINTS:
        loss[e] = {}
        for model in MODELS:
            sse = sum((math.log1p(row[e]) -
                       v46.prediction(fit, row, model, e)) ** 2
                      for row in test)
            val = sse / len(test)
            require(math.isfinite(val) and val >= 0,
                    "V49_INVALID_OUT_OF_SAMPLE_MSE")
            loss[e][model] = val
    return {"heldout": heldout, "training_pops": train_pops,
            "two_focal_training_populations": [p for p in train_pops
                                                if poptax[p] == focal],
            "n_train_ramets": len(train), "n_test_ramets": len(test),
            "loss": loss}


def check_deamii_against_v46(fold, parent, tolerance):
    target = fold["heldout"]
    require(target in {"1719", "1720", "1721"}
            and len(fold["training_pops"]) == 19
            and len(fold["two_focal_training_populations"]) == 2,
            "V49_DEAMII_ORIGINAL_FULL_TRAINING_NOT_RESTORED")
    for e in ENDPOINTS:
        old = {x["population"]: x for x in
               parent["endpoints"][e]["per_population_holdout"]}
        require(target in old
                and fold["n_test_ramets"] == old[target]["n_ramets"],
                "V49_V46_DEAMII_ORIGINAL_FOLD_MISSING")
        for m in MODELS:
            require(
                math.isclose(fold["loss"][e][m],
                             old[target]["mse_by_model"][m],
                             rel_tol=0, abs_tol=tolerance),
                "V49_V46_SAME_MODEL_DEAMII_ANCHOR_DIFFERENT_" + e + "_" + m,
            )


def panel_gain(folds, e, model):
    require(model != "B0" and len(folds) == 3,
            "V49_NOT_THREE_FOLD_PANEL")
    original = sum(f["loss"][e]["B0"] for f in folds)
    candidate = sum(f["loss"][e][model] for f in folds)
    require(original > 0 and math.isfinite(candidate),
            "V49_PANEL_BASELINE_RISK_INVALID")
    return 1 - candidate / original


def num(x):
    require(math.isfinite(x), "V49_NONFINITE_SOURCE_OUTPUT")
    return round(x, 9)


def describe_reference(scores, deamii, tol):
    scores = sorted(scores)
    n = len(scores)
    require(n > 0, "V49_EMPTY_COMPARISON_REFERENCE")
    return {
        "n_overlapping_panels": n,
        "n_negative": sum(x < 0 for x in scores),
        "n_at_or_below_deamii": sum(x <= deamii + tol for x in scores),
        "min": num(scores[0]), "median": num(statistics.median(scores)),
        "max": num(scores[-1]),
        "independent_field_studies": 0, "randomization_p_value": False
    }


def summarize(rows, audit, parent, r):
    require(len(rows) == 963 and audit["rows"] == 963,
            "V49_BAD_ORIGINAL_INPUT_DENOMINATOR")
    bypop, poptax, bytaxon = original_populations(rows)
    panels = []
    values = {e: {m: {t: [] for t in TAXA} for m in MODELS[1:]}
              for e in ENDPOINTS}
    n_folds = 0
    for taxon in TAXA:
        for triplet in itertools.combinations(bytaxon[taxon], 3):
            folds = [score_one_fold(rows, bypop, poptax, taxon, triplet, heldout)
                     for heldout in triplet]
            if taxon == "deamii":
                for f in folds:
                    check_deamii_against_v46(
                        f, parent, r["deamii_v46_mse_replay_tolerance"])
            n_folds += len(folds)
            gains = {}
            for e in ENDPOINTS:
                gains[e] = {}
                for m in MODELS[1:]:
                    raw = panel_gain(folds, e, m)
                    values[e][m][taxon].append(raw)
                    gains[e][m] = num(raw)
            panels.append({"taxon": taxon, "three_focal_population_ids":
                           list(triplet), "heldout_folds": [
                {**f, "loss": {
                    e: {m: num(f["loss"][e][m]) for m in MODELS}
                    for e in ENDPOINTS}} for f in folds
                ], "equal_population_panel_gains": gains})
    require(len(panels) == 141 and n_folds == 423
            and Counter(p["taxon"] for p in panels) ==
                Counter({"amoena": 56, "deamii": 1, "pilosa": 84}),
            "V49_EXPECT_141_PANELS_423_FULL_POP_HOLDOUTS")
    comparisons = {}
    for e in ENDPOINTS:
        comparisons[e] = {}
        for m in MODELS[1:]:
            d = values[e][m]["deamii"][0]
            a = values[e][m]["amoena"]
            p = values[e][m]["pilosa"]
            require(len(a) == 56 and len(p) == 84,
                    "V49_REFERENCES_NOT_56_PLUS_84")
            comparisons[e][m] = {
                "deamii_original_v46_full_train_population_mean_gain": num(d),
                "amoena_56_reference": describe_reference(
                    a, d, r["reference_comparison_tolerance"]),
                "pilosa_84_reference": describe_reference(
                    p, d, r["reference_comparison_tolerance"]),
                "nondeamii_140_reference": describe_reference(
                    a+p, d, r["reference_comparison_tolerance"]),
            }
            old = parent["endpoints"][e]["model_results"][m]["by_taxon"]["deamii"]["relative_equal_population_mse_gain"]
            require(math.isclose(d, old, rel_tol=0, abs_tol=1e-7),
                    "V49_ORIGINAL_V46_DEAMII_MEAN_RISK_DRIFT")
    return {
        "analysis_type": r["analysis_type"],
        "source_original_CSV_git_blob_sha1": r["source_CSV_git_blob_sha1"],
        "v46_parent_derived_git_blob_sha1": r["v46_derived_result_git_blob_sha1"],
        "original_grain": audit,
        "previous_v46_primary_transfer_established": False,
        "new_original_field_studies": 0,
        "new_original_common_gardens": 0,
        "source_population_model_fits": n_folds,
        "n_same_taxon_panels": len(panels),
        "all_nonfocal_source_populations_preserved_in_training": True,
        "v46_unmodified_five_model_ridge_and_nine_intercepts_reused": True,
        "all_deamii_per_population_loss_reproduced_from_v46": True,
        "biomass_log1p_technical_nonadmission_retained": True,
        "new_p_values_or_confidence_intervals": False,
        "original_panel_fit_details": panels,
        "comparisons": comparisons,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", default=None)
    p.add_argument("--check", action="store_true")
    args = p.parse_args()
    r = registry()
    parent = frozen_v46(r)
    rows, audit = author_source(args.source, r)
    output = summarize(rows, audit, parent, r)
    if args.check:
        require(output == json.loads(RECEIPT.read_text()),
                "V49_ORIGINAL_SOURCE_FULL_DERIVED_RECEIPT_MISMATCH")
        print("V49_ORIGINAL_AUTHOR_SHA_423_FITS_V46_DEAMII_ANCHOR_EXACT_PASS")
    else:
        raw = json.dumps(output, separators=(",", ":"), sort_keys=True,
                         ensure_ascii=True)
        n_chunks = (len(raw) + 3999) // 4000
        for i in range(n_chunks):
            print("V49_SOURCE_READOUT_CHUNK_%03d_OF_%03d=" %
                  (i+1, n_chunks) + raw[4000*i:4000*(i+1)])


if __name__ == "__main__":
    main()
