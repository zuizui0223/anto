#!/usr/bin/env python3
"""ANTO v50: source-locked Phlox 2024 ridge coordinates/sharing diagnostics.

Posthoc to v46-v49. Not a phylogenetic replication, flower color,
natural selection coefficient, or a new original study.
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
import analyze_Phlox_2024_multitaxon_focal2_v49 as v49

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data/external/anto_v50_ridge_sharing_original_source_registry_20261119.json"
PARENT = ROOT / "results/anto_v49_original_Phlox_multitaxon_focal2_20261119.json"
RECEIPT = ROOT / "results/anto_v50_sharing_penalty_source_sensitivity_20261119.json"
MODELS = ("G", "T", "E", "TE")
ENDPOINTS = ("fruits_2019", "flowers_2019")
TAXA = ("amoena", "deamii", "pilosa")
SCENARIOS = ("STD_R0", "STD_R1", "STD_R10", "STD_R100",
             "RAW_R0", "RAW_R10", "RAW_EQUIV_R10")
N_BLOCKS = {"G": 1, "T": 3, "E": 3, "TE": 9}


class V50Error(ValueError):
    pass


def require(ok, reason):
    if not ok:
        raise V50Error(reason)


def git_sha1(raw):
    return hashlib.sha1(("blob %s\0" % len(raw)).encode() + raw).hexdigest()


def source_registry():
    r = json.loads(REGISTRY.read_text())
    require(
        r["scope"] ==
            "POSTHOC_ORIGINAL_2024_PHLOX_FROZEN_SOURCE_RIDGE_AND_SHARING_DECOMPOSITION"
        and r["source_repo"] == "PhloxHopkins/PhloxFieldAdaptation"
        and r["source_commit"] ==
            "d792fe58cce36721394adfd97e9f412396578e86"
        and r["source_git_blob_sha1"] ==
            "4c780b0abe644ac69dae86a0450eb598c160626b"
        and r["source_size_bytes"] == 175377
        and r["v49_result_blob_sha1"] ==
            "94947d03e85766a3c5c1f82fde2bb076f2ca42e8"
        and r["n_ramets"] == 963 and r["n_genotypes"] == 284
        and r["n_source_populations"] == 20
        and r["n_source_populations_by_taxon"] ==
            {"amoena": 8, "deamii": 3, "pilosa": 9}
        and r["panels_by_taxon"] ==
            {"amoena": 56, "deamii": 1, "pilosa": 84}
        and r["n_panels"] == 141 and r["n_folds"] == 423
        and r["original_v46_models"] == ["B0", "G", "T", "E", "TE"]
        and r["primary_endpoint"] == "fruits_2019"
        and r["secondary_endpoint"] == "flowers_2019"
        and [x["name"] for x in r["scenarios"]] == list(SCENARIOS)
        and [x["lambda"] for x in r["scenarios"]] ==
            [0, 1, 10, 100, 0, 10, 10]
        and [x["coordinate"] for x in r["scenarios"]] ==
            ["standardized"] * 4 + ["raw", "raw", "raw_equivalent"]
        and r["reference_comparison_tol"] == 1e-10
        and r["source_mse_parity_tolerance"] == 1e-7
        and r["n_negative_biomass"] == 13
        and r["min_biomass"] == -233
        and r["new_original_biological_studies"] == 0
        and r["new_flower_colour_data"] is False
        and r["new_pollinator_data"] is False
        and r["new_genetic_persistence"] is False
        and r["original_v46_all_taxa_transfer_falsified"] is True
        and r["is_posthoc"] is True
        and r["valid_randomization_p"] is False,
        "V50_SOURCE_METHOD_OR_ADMISSIBILITY_REGISTRY_DRIFT",
    )
    return r


def original_data(path, registry):
    old = v46.source_registry()
    require(old["source_commit"] == registry["source_commit"]
            and old["source_git_blob_sha1"] ==
                 registry["source_git_blob_sha1"]
            and old["source_size"] == registry["source_size_bytes"],
            "V50_V46_ORIGINAL_CSV_SOURCE_CHANGED")
    rows, audit = v46.parse_source(v46.acquire_raw(old, path), old)
    require(audit["rows"] == 963 and audit["genotypes"] == 284
            and audit["populations"] == 20
            and audit["biomass_original_negative_n"] == 13
            and audit["biomass_original_minimum"] == -233,
            "V50_SOURCE_COHORT_CHANGED")
    return rows, audit


def v49_source_parent(registry):
    content = PARENT.read_bytes()
    require(git_sha1(content) == registry["v49_result_blob_sha1"],
            "V50_PARENT_SOURCE_V49_GIT_SHA_CHANGED")
    parent = json.loads(content)
    require(parent["analysis_type"] ==
            "POSTHOC_MULTI_TAXON_FOCAL_2_TRAIN_SOURCE_POP_REFIT"
            and parent["source_original_CSV_git_blob_sha1"] ==
                registry["source_git_blob_sha1"]
            and parent["source_population_model_fits"] == 423
            and parent["n_same_taxon_panels"] == 141
            and parent["previous_v46_primary_transfer_established"] is False
            and parent["all_deamii_per_population_loss_reproduced_from_v46"]
            and parent["new_original_field_studies"] == 0,
            "V50_PARENT_V49_CLAIMS_CHANGED")
    return parent


def stable_number(x):
    require(math.isfinite(x), "V50_NONFINITE_DERIVED_SCORE")
    return round(x, 9)


def ridge_two_by_two(a, b, c, u, v, penalty_a, penalty_b):
    require(all(math.isfinite(z) for z in
                (a, b, c, u, v, penalty_a, penalty_b))
            and a >= 0 and c >= 0 and penalty_a >= 0
            and penalty_b >= 0, "V50_INVALID_RIDGE_NORMAL_EQUATION")
    aa = a + penalty_a
    cc = c + penalty_b
    det = aa * cc - b * b
    require(det > max(1e-12, 1e-12 * (aa * cc)),
            "V50_NONIDENTIFIABLE_OR_SINGULAR_2X2")
    x = (cc * u - b * v) / det
    y = (aa * v - b * u) / det
    require(math.isfinite(x) and math.isfinite(y),
            "V50_NONFINITE_FITTED_RIDGE_SLOPE")
    return (x, y)


def block_id(model, row):
    require(model in MODELS, "V50_UNREGISTERED_MODEL")
    return v46.slope_offset(model, row) // 2


def frozen_panel_index(parent):
    idx = {}
    for p in parent["original_panel_fit_details"]:
        key = (p["taxon"], tuple(p["three_focal_population_ids"]))
        require(key not in idx, "V50_DUPLICATE_V49_FROZEN_PANEL")
        idx[key] = p
    require(len(idx) == 141, "V50_PARENT_NOT_141_UNIQUE_PANELS")
    return idx


def fit_one_training(rows):
    """Training-only nine B0 means and standardized/raw 2-PC normal equations."""
    # Explicitly replay v46's exact frozen lambda=10 baseline and features.
    original = v46.fit_train(rows)
    sd1 = original["sd_pc1"]
    sd2 = original["sd_pc2"]
    require(sd1 > 1e-10 and sd2 > 1e-10,
            "V50_TRAIN_GLOBAL_PC_SCALE_DEGENERATE")
    b = {}
    for model in MODELS:
        blocks = []
        for _ in range(N_BLOCKS[model]):
            blocks.append({"rr11": 0.0, "rr12": 0.0, "rr22": 0.0,
                           "rhs": {e: [0.0, 0.0] for e in ENDPOINTS}})
        for row in rows:
            z = original["group_mean"][v46.group(row)]
            x1 = row["leaf_traits_pc1"] - z["pc1"]
            x2 = row["leaf_traits_pc2"] - z["pc2"]
            block = blocks[block_id(model, row)]
            block["rr11"] += x1*x1
            block["rr12"] += x1*x2
            block["rr22"] += x2*x2
            for e in ENDPOINTS:
                residual = math.log1p(row[e]) - z["y"][e]
                block["rhs"][e][0] += x1 * residual
                block["rhs"][e][1] += x2 * residual
        b[model] = blocks
    all_beta = {}
    for scenario in SCENARIOS:
        std = scenario.startswith("STD_")
        raw_equiv = scenario == "RAW_EQUIV_R10"
        alpha = 0.0 if scenario.endswith("R0") else (
            1.0 if scenario == "STD_R1" else
            100.0 if scenario == "STD_R100" else 10.0)
        fit_model = {}
        for model in MODELS:
            fit_model[model] = []
            for block in b[model]:
                if std:
                    a = block["rr11"] / (sd1*sd1)
                    ab = block["rr12"] / (sd1*sd2)
                    c = block["rr22"] / (sd2*sd2)
                    p1 = alpha
                    p2 = alpha
                else:
                    a = block["rr11"]
                    ab = block["rr12"]
                    c = block["rr22"]
                    p1 = alpha*sd1*sd1 if raw_equiv else alpha
                    p2 = alpha*sd2*sd2 if raw_equiv else alpha
                found = {}
                for e in ENDPOINTS:
                    rhs = block["rhs"][e]
                    u = rhs[0]/sd1 if std else rhs[0]
                    v = rhs[1]/sd2 if std else rhs[1]
                    found[e] = ridge_two_by_two(a, ab, c, u, v, p1, p2)
                fit_model[model].append(found)
        all_beta[scenario] = fit_model
    return original, all_beta


def predict(original, betas, scenario, model, row, e):
    s = original["group_mean"][v46.group(row)]
    x1 = row["leaf_traits_pc1"] - s["pc1"]
    x2 = row["leaf_traits_pc2"] - s["pc2"]
    if scenario.startswith("STD_"):
        x1 /= original["sd_pc1"]
        x2 /= original["sd_pc2"]
    v1, v2 = betas[scenario][model][block_id(model, row)][e]
    return s["y"][e] + x1*v1 + x2*v2


def evaluate_fold(rows, chosen_pops, heldout, parent_fold, tolerance):
    selected = set(chosen_pops)
    require(heldout not in selected
            and len(set(chosen_pops)) == len(chosen_pops),
            "V50_SOURCE_FOLD_SELF_TRAINING_OR_DUPLICATE")
    train = [r for r in rows if r["population"] in selected]
    test = [r for r in rows if r["population"] == heldout]
    require(train and test and len(test) == parent_fold["n_test_ramets"]
            and len(train) == parent_fold["n_train_ramets"]
            and chosen_pops == parent_fold["training_pops"]
            and {r["genotype"] for r in train}.isdisjoint(
                {r["genotype"] for r in test}),
            "V50_SOURCE_POPULATION_OR_CLONE_FOLD_DRIFT")
    original, beta = fit_one_training(train)
    scores = {e: {} for e in ENDPOINTS}
    baseline = {}
    max_err = {"v49_original_std10": 0.0,
               "raw_ols_vs_std_ols": 0.0,
               "raw_equiv_vs_std_ridge": 0.0}
    for e in ENDPOINTS:
        b0 = sum((math.log1p(row[e]) -
                  v46.prediction(original, row, "B0", e))**2
                 for row in test) / len(test)
        baseline[e] = b0
        require(math.isclose(b0, parent_fold["loss"][e]["B0"],
                             rel_tol=0, abs_tol=tolerance),
                "V50_PARENT_V49_B0_NOT_REPRODUCED")
        for scenario in SCENARIOS:
            scores[e][scenario] = {}
            for model in MODELS:
                loss = sum((math.log1p(row[e]) -
                            predict(original, beta, scenario, model,
                                    row, e))**2
                           for row in test) / len(test)
                scores[e][scenario][model] = loss
        for model in MODELS:
            for key, left, right, target in (
                ("v49_original_std10", scores[e]["STD_R10"][model],
                 parent_fold["loss"][e][model], tolerance),
                ("raw_ols_vs_std_ols", scores[e]["RAW_R0"][model],
                 scores[e]["STD_R0"][model], 1e-7),
                ("raw_equiv_vs_std_ridge",
                 scores[e]["RAW_EQUIV_R10"][model],
                 scores[e]["STD_R10"][model], 1e-7)
            ):
                diff = abs(left-right)
                max_err[key] = max(max_err[key], diff)
                require(diff <= target,
                        "V50_MANDATORY_NUMERIC_INVARIANCE_FAILED_" +
                        key + "_" + model + "_" + e)
    return {"n_training_ramets": len(train),
            "n_heldout_ramets": len(test),
            "b0": baseline, "scores": scores,
            "max_errors": max_err}


def panel_gain(folds, e, scenario, model):
    base = sum(x["b0"][e] for x in folds)
    predicted = sum(x["scores"][e][scenario][model] for x in folds)
    require(base > 0 and math.isfinite(predicted),
            "V50_BAD_PANEL_BASELINE_DENOMINATOR")
    return 1-predicted/base


def distribution(a, d, tolerance):
    vals = sorted(a)
    require(len(vals) > 0, "V50_EMPTY_REFERENCE")
    return {"overlapping_panel_count": len(vals),
            "negative_panels": sum(x < 0 for x in vals),
            "at_or_below_deamii": sum(x <= d+tolerance for x in vals),
            "minimum": stable_number(vals[0]),
            "median": stable_number(statistics.median(vals)),
            "maximum": stable_number(vals[-1]),
            "independent_studies": 0, "is_p_value": False}


def evaluate(rows, audit, parent, source):
    pp = frozen_panel_index(parent)
    bypop, poptax, popkeys = v49.original_populations(rows)
    scores = {e: {s: {m: {taxon: [] for taxon in TAXA}
                      for m in MODELS} for s in SCENARIOS}
              for e in ENDPOINTS}
    panel_receipts = []
    fold_count = 0
    max_parity = {
        "v49_original_std10": 0.0,
        "raw_ols_vs_std_ols": 0.0,
        "raw_equiv_vs_std_ridge": 0.0,
    }
    hash_fold_scores = hashlib.sha256()
    deamii_fold_details = []
    for taxon in TAXA:
        for panel in itertools.combinations(popkeys[taxon], 3):
            original_panel = pp.get((taxon, tuple(panel)))
            require(original_panel is not None, "V50_PARENT_V49_PANEL_MISSING")
            folds = []
            compact = {"taxon": taxon, "source_populations": list(panel),
                       "gains": {}}
            for target in panel:
                source_fold = next((x for x in original_panel["heldout_folds"]
                                    if x["heldout"] == target), None)
                require(source_fold is not None,
                        "V50_PARENT_V49_FOLD_MISSING")
                next_fold = evaluate_fold(rows, source_fold["training_pops"],
                                          target, source_fold,
                                          source["source_mse_parity_tolerance"])
                for k in max_parity:
                    max_parity[k] = max(max_parity[k],
                                        next_fold["max_errors"][k])
                # A complete canonical digest covers every original
                # scenario×model×endpoint fold loss before rounding.
                payload = {
                    "taxon": taxon, "panel": panel, "target": target,
                    "n_train": next_fold["n_training_ramets"],
                    "n_test": next_fold["n_heldout_ramets"],
                    "b0": next_fold["b0"],
                    "losses": next_fold["scores"],
                }
                hash_fold_scores.update(
                    (json.dumps(payload, sort_keys=True,
                                separators=(",", ":"), allow_nan=False)
                     + "\n").encode("utf-8"))
                if taxon == "deamii":
                    deamii_fold_details.append({
                        "population": target,
                        "scores": {
                            e: {s: {m: stable_number(
                                next_fold["scores"][e][s][m])
                                for m in MODELS} for s in SCENARIOS}
                            for e in ENDPOINTS}})
                folds.append(next_fold)
                fold_count += 1
            for e in ENDPOINTS:
                compact["gains"][e] = {}
                for s in SCENARIOS:
                    compact["gains"][e][s] = {}
                    for model in MODELS:
                        gain = panel_gain(folds, e, s, model)
                        scores[e][s][model][taxon].append(gain)
                        compact["gains"][e][s][model] = stable_number(gain)
            panel_receipts.append(compact)
    require(fold_count == 423 and len(panel_receipts) == 141
            and Counter(x["taxon"] for x in panel_receipts) ==
                Counter({"amoena": 56, "deamii": 1, "pilosa": 84}),
            "V50_ORIGINAL_141_PANELS_423_FOLDS_NOT_COMPLETE")
    comparisons = {}
    for e in ENDPOINTS:
        comparisons[e] = {}
        for s in SCENARIOS:
            comparisons[e][s] = {}
            for model in MODELS:
                d = scores[e][s][model]["deamii"][0]
                amoena = scores[e][s][model]["amoena"]
                pilosa = scores[e][s][model]["pilosa"]
                require(len(amoena)==56 and len(pilosa)==84,
                        "V50_NONDEAMII_REFERENCE_NOT_140")
                comparisons[e][s][model] = {
                    "deamii_original_population_gain": stable_number(d),
                    "nondeamii_140_reference": distribution(
                        amoena+pilosa,d,source["reference_comparison_tol"]),
                    "amoena_56": distribution(
                        amoena,d,source["reference_comparison_tol"]),
                    "pilosa_84": distribution(
                        pilosa,d,source["reference_comparison_tol"]),
                }
    for e in ENDPOINTS:
        for model in MODELS:
            now = comparisons[e]["STD_R10"][model]
            old = parent["comparisons"][e][model]
            require(math.isclose(
                now["deamii_original_population_gain"],
                old["deamii_original_v46_full_train_population_mean_gain"],
                rel_tol=0,abs_tol=1e-7)
                and (now["nondeamii_140_reference"]["at_or_below_deamii"] ==
                     old["nondeamii_140_reference"]["n_at_or_below_deamii"]),
                "V50_V49_REFERENCE_SUMMARY_DRIFT")
    return {
        "analysis_type": source["scope"],
        "source_csv_git_blob_sha1": source["source_git_blob_sha1"],
        "v49_parent_result_git_blob_sha1": source["v49_result_blob_sha1"],
        "original_source_grain": audit,
        "new_original_biological_studies": 0,
        "n_original_source_populations": 20,
        "source_population_folds": fold_count,
        "source_original_three_pop_panels": len(panel_receipts),
        "no_original_source_raw_rows_redistributed": True,
        "v46_all_taxa_fruit_transfer_still_NOT_ESTABLISHED": True,
        "original_biomass_technical_nonadmission_retained": True,
        "is_posthoc": True,
        "statistical_p_values_or_confidence_intervals_computed": False,
        "v49_source_and_numerical_parity_max_absolute_mse_errors": max_parity,
        "full_unrounded_fold_loss_sha256": hash_fold_scores.hexdigest(),
        "full_original_three_pop_panel_source_derived_gains": panel_receipts,
        "original_deamii_three_fold_diagnostic_scores": deamii_fold_details,
        "comparisons": comparisons
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",default=None)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    fixed=source_registry()
    parent=v49_source_parent(fixed)
    rows,audit=original_data(args.source,fixed)
    output=evaluate(rows,audit,parent,fixed)
    if args.check:
        require(output == json.loads(RECEIPT.read_text()),
                "V50_SOURCE_DERIVED_EXACT_STORED_RECEIPT_DRIFT")
        print("V50_ORIGINAL_SHA_423_FOLDS_SEVEN_CASES_NUMERIC_INVARIANCE_AND_STORED_RECEIPT_PASS")
    else:
        raw=json.dumps(output,sort_keys=True,separators=(",",":"),
                       ensure_ascii=True)
        nchunks=(len(raw)+3999)//4000
        for i in range(nchunks):
            print("V50_SOURCE_READOUT_CHUNK_%03d_OF_%03d="%
                  (i+1,nchunks)+raw[4000*i:4000*(i+1)])


if __name__=="__main__":
    main()
