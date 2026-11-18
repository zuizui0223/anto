#!/usr/bin/env python3
"""ANTO v46: original 2024 Phlox population-heldout fitness prediction.

This is a single original cross-taxon leaf-trait study, NOT a flower-colour
assay, 2014 PhloxSelection refit, causal model, or new independent garden.
All original rows are accessed only on a temporary runner/source path.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data/external/anto_v46_phlox2024_source_registry_20261118.json"
RECEIPT = ROOT / "results/anto_v46_phlox2024_population_heldout_20261118.json"
TAXA = ("amoena", "deamii", "pilosa")
GARDENS = ("1-Amoena-PB", "2-Deamii-JC", "3-Pilosa-LBL")
MODELS = ("B0", "G", "T", "E", "TE")
ACCEPTED_ENDPOINTS = ("fruits_2019", "flowers_2019")
BLOCKED_ENDPOINT = "biomass_2019"
KEYS = ("label", "garden", "bed", "taxon", "HvA", "population",
        "genotype", "longitude", "latitude", "herbivory_2019",
        "flowers_2019", "fruits_2019", "biomass_2019", "survival_2020",
        "length", "width", "length_width_ratio", "area",
        "chlorophyll", "specific_leaf_area", "leaf_traits_pc1",
        "leaf_traits_pc2", "leaf_traits_pc3", "geographic_distance",
        "genetic_distance", "environmental_distance")


class V46Error(ValueError):
    pass


def require(condition, code):
    if not condition:
        raise V46Error(code)


def git_blob_sha1(raw):
    return hashlib.sha1(("blob %d\0" % len(raw)).encode() + raw).hexdigest()


def source_registry():
    r = json.loads(REGISTRY.read_text())
    require(r["source_repo"] == "PhloxHopkins/PhloxFieldAdaptation"
            and r["source_commit"] ==
            "d792fe58cce36721394adfd97e9f412396578e86"
            and r["source_file"] ==
            "Dataset_fitness_leaftraits_distance_lat_long.csv"
            and r["source_git_blob_sha1"] ==
            "4c780b0abe644ac69dae86a0450eb598c160626b"
            and r["source_size"] == 175377,
            "V46_AUTHOR_SOURCE_REGISTRY_DRIFT")
    require(r["rows"] == 963 and r["columns"] == 26
            and r["genotypes"] == 284 and r["populations"] == 20
            and r["models"] == list(MODELS)
            and r["primary"] == "fruits_2019"
            and r["secondary"] == ["flowers_2019", "biomass_2019"]
            and r["ridge_alpha"] == 10,
            "V46_FROZEN_ANALYSIS_DRIFT")
    require(r["source_programs_added"] if "source_programs_added" in r
            else r["new_original_study_count"] == 1,
            "V46_ORIGINAL_STUDY_UNIT_CHANGED")
    require(r["flower_colour_assay"] is False
            and r["trait_ancestral_regain_assay"] is False
            and r["independent_original_2014_sample"] is False,
            "V46_FALSE_FLORAL_MECHANISM_CLAIM")
    return r


def acquire_raw(registry, source=None):
    if source is None:
        url = ("https://raw.githubusercontent.com/" + registry["source_repo"]
               + "/" + registry["source_commit"] + "/" + registry["source_file"])
        with urllib.request.urlopen(url, timeout=90) as response:
            raw = response.read()
    else:
        raw = Path(source).read_bytes()
    require(len(raw) == registry["source_size"],
            "V46_AUTHOR_ORIGINAL_CSV_BYTE_COUNT_DRIFT")
    require(git_blob_sha1(raw) == registry["source_git_blob_sha1"],
            "V46_AUTHOR_ORIGINAL_CSV_GIT_BLOB_DRIFT")
    return raw


def parse_source(raw, registry):
    require(len(raw) == registry["source_size"]
            and git_blob_sha1(raw) == registry["source_git_blob_sha1"],
            "V46_ORIGINAL_CSV_SOURCE_IDENTITY_UNVERIFIED")
    try:
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline=""))
        require(reader.fieldnames == list(KEYS), "V46_SOURCE_COLUMNS_CHANGED")
        rows = list(reader)
    except UnicodeDecodeError as exc:
        raise V46Error("V46_ORIGINAL_UTF8_UNREADABLE") from exc
    require(len(rows) == registry["rows"], "V46_ORIGINAL_ROW_DENOMINATOR_DRIFT")
    require(all(None not in row and all(x is not None for x in row.values())
                for row in rows), "V46_INCOMPLETE_SOURCE_ROW")
    labels, genotypes, original_pop, pop_tax, gbytax = set(), {}, set(), {}, Counter()
    count_tax, count_garden = Counter(), Counter()
    negative_biomass = []
    for i, row in enumerate(rows, 1):
        require(row["label"] not in labels, "V46_RAMET_LABEL_DUPLICATED")
        labels.add(row["label"])
        garden, taxon, pop, geno = (row[k] for k in
                                    ("garden", "taxon", "population", "genotype"))
        require(garden in GARDENS and taxon in TAXA,
                "V46_SOURCE_TAXON_GARDEN_CHANGED")
        require((row["HvA"] == "home") ==
                ((garden == "1-Amoena-PB" and taxon == "amoena")
                 or (garden == "2-Deamii-JC" and taxon == "deamii")
                 or (garden == "3-Pilosa-LBL" and taxon == "pilosa")),
                "V46_SOURCE_HOME_AWAY_DRIFT")
        original_pop.add(pop)
        count_tax[taxon] += 1
        count_garden[garden] += 1
        gbytax[(garden, taxon)] += 1
        require(pop not in pop_tax or pop_tax[pop] == taxon,
                "V46_SOURCE_POPULATION_SPANS_TAXA")
        pop_tax[pop] = taxon
        d = {}
        for key in ("leaf_traits_pc1", "leaf_traits_pc2",
                    "fruits_2019", "flowers_2019", "biomass_2019"):
            try:
                v = float(row[key])
            except (TypeError, ValueError) as exc:
                raise V46Error("V46_BAD_ORIGINAL_NUMERIC_" + key) from exc
            require(math.isfinite(v), "V46_NONFINITE_ORIGINAL_" + key)
            d[key] = v
        for key in ("fruits_2019", "flowers_2019"):
            require(d[key] >= 0
                    and math.isclose(d[key], round(d[key]), abs_tol=1e-8),
                    "V46_INVALID_PRIMARY_OR_SECONDARY_COUNT_" + key)
        if d["biomass_2019"] < 0:
            negative_biomass.append(d["biomass_2019"])
        info = genotypes.setdefault(geno, {"pop": pop, "taxon": taxon,
                    "pc": (d["leaf_traits_pc1"], d["leaf_traits_pc2"]),
                    "garden_count": Counter(), "n": 0})
        require(info["pop"] == pop and info["taxon"] == taxon,
                "V46_CLONE_SPANS_ORIGINAL_SOURCE_POPULATIONS")
        require(info["pc"] == (d["leaf_traits_pc1"], d["leaf_traits_pc2"]),
                "V46_CLONE_LEAF_PCS_INCONSISTENT")
        info["garden_count"][garden] += 1
        info["n"] += 1
        row.update(d)
    require(set(pop_tax) == original_pop and len(original_pop) == 20,
            "V46_SOURCE_POPULATION_COUNT_CHANGED")
    require(len(genotypes) == registry["genotypes"],
            "V46_ORIGINAL_GENOTYPE_COUNT_CHANGED")
    require(count_tax == Counter(registry["taxon_counts"])
            and count_garden == Counter(registry["garden_counts"]),
            "V46_ORIGINAL_GROUP_SAMPLE_SIZES_CHANGED")
    for t in TAXA:
        require(sum(1 for p in pop_tax if pop_tax[p] == t)
                == registry["populations_by_taxon"][t],
                "V46_NESTED_SOURCE_TAXON_DENOMINATOR_CHANGED")
        for g in GARDENS:
            require(gbytax[(g, t)] ==
                    registry["garden_taxon_count_per_garden"][t],
                    "V46_NINE_SOURCE_CELL_DENOMINATOR_CHANGED")
    require(all(set(info["garden_count"]) == set(GARDENS)
                and info["n"] in (3, 6) for info in genotypes.values()),
            "V46_CLONE_GARDEN_REPLICATION_DRIFT")
    require(len(negative_biomass) == 13
            and min(negative_biomass) == -233.0,
            "V46_ORIGINAL_BIOMASS_NONADMISSION_COUNT_DRIFT")
    # NONADMISSION: preregistered log1p biomass is undefined for
    # negative source observations; do not clip, drop, shift or refit.
    return rows, {"rows": len(rows), "genotypes": len(genotypes),
            "populations": len(original_pop),
            "populations_by_taxon": registry["populations_by_taxon"],
            "gardens": registry["garden_counts"],
            "biomass_original_negative_n": len(negative_biomass),
            "biomass_original_minimum": min(negative_biomass)}


def group(row):
    return (row["garden"], row["taxon"])


def fit_train(train, endpoints=ACCEPTED_ENDPOINTS, lam=10.0):
    require(len(train) > 0 and lam == 10, "V46_CHANGED_TRAIN_PENALTY")
    stats = {}
    for r in train:
        cell = group(r)
        s = stats.setdefault(cell, {"n": 0, "pc1": 0.0, "pc2": 0.0,
                                    "y": {e: 0.0 for e in endpoints}})
        s["n"] += 1
        s["pc1"] += r["leaf_traits_pc1"]
        s["pc2"] += r["leaf_traits_pc2"]
        for e in endpoints:
            s["y"][e] += math.log1p(r[e])
    require(set(stats) == {(g, t) for g in GARDENS for t in TAXA},
            "V46_SOURCE_CELL_ABSENT_IN_TRAIN")
    for s in stats.values():
        s["pc1"] /= s["n"]
        s["pc2"] /= s["n"]
        for e in endpoints:
            s["y"][e] /= s["n"]
    sd1 = math.sqrt(sum((r["leaf_traits_pc1"] - stats[group(r)]["pc1"]) ** 2
                        for r in train) / len(train))
    sd2 = math.sqrt(sum((r["leaf_traits_pc2"] - stats[group(r)]["pc2"]) ** 2
                        for r in train) / len(train))
    require(sd1 > 1e-10 and sd2 > 1e-10, "V46_TRAIN_FEATURE_UNIDENTIFIABLE")
    betas = {}
    for model in MODELS[1:]:
        ncols = {"G": 2, "T": 6, "E": 6, "TE": 18}[model]
        lhs = [[lam if i == j else 0.0 for j in range(ncols)]
               for i in range(ncols)]
        rhs = {e: [0.0] * ncols for e in endpoints}
        for r in train:
            cell = group(r)
            s = stats[cell]
            z1 = (r["leaf_traits_pc1"] - s["pc1"]) / sd1
            z2 = (r["leaf_traits_pc2"] - s["pc2"]) / sd2
            p = slope_offset(model, r)
            lhs[p][p] += z1 * z1
            lhs[p][p + 1] += z1 * z2
            lhs[p + 1][p] += z2 * z1
            lhs[p + 1][p + 1] += z2 * z2
            for e in endpoints:
                residual = math.log1p(r[e]) - s["y"][e]
                rhs[e][p] += z1 * residual
                rhs[e][p + 1] += z2 * residual
        betas[model] = {e: solve(lhs, rhs[e]) for e in endpoints}
    return {"group_mean": stats, "sd_pc1": sd1, "sd_pc2": sd2,
            "coefficients": betas}


def slope_offset(model, row):
    if model == "G":
        return 0
    if model == "T":
        return 2 * TAXA.index(row["taxon"])
    if model == "E":
        return 2 * GARDENS.index(row["garden"])
    if model == "TE":
        return 2 * (len(TAXA) * GARDENS.index(row["garden"])
                    + TAXA.index(row["taxon"]))
    raise V46Error("V46_UNKNOWN_MODEL")


def solve(lhs, rhs):
    n = len(rhs)
    require(len(lhs) == n and all(len(a) == n for a in lhs),
            "V46_BAD_RIDGE_MATRIX_SHAPE")
    a = [list(lhs[i]) + [rhs[i]] for i in range(n)]
    for j in range(n):
        piv = max(range(j, n), key=lambda i: abs(a[i][j]))
        require(abs(a[piv][j]) > 1e-12, "V46_SINGULAR_RIDGE")
        a[j], a[piv] = a[piv], a[j]
        den = a[j][j]
        for k in range(j, n + 1):
            a[j][k] /= den
        for i in range(n):
            if i != j:
                factor = a[i][j]
                for k in range(j, n + 1):
                    a[i][k] -= factor * a[j][k]
    result = [a[i][n] for i in range(n)]
    require(all(math.isfinite(z) for z in result),
            "V46_NONFINITE_COEFFICIENT")
    return result


def prediction(fitted, row, model, endpoint):
    s = fitted["group_mean"][group(row)]
    base = s["y"][endpoint]
    if model == "B0":
        return base
    z1 = (row["leaf_traits_pc1"] - s["pc1"]) / fitted["sd_pc1"]
    z2 = (row["leaf_traits_pc2"] - s["pc2"]) / fitted["sd_pc2"]
    p = slope_offset(model, row)
    b = fitted["coefficients"][model][endpoint]
    return base + z1 * b[p] + z2 * b[p + 1]


def cohort_folds(rows):
    pops = sorted(set(r["population"] for r in rows))
    require(len(pops) == 20, "V46_EXPECT_TWENTY_SOURCE_POPULATIONS")
    folds = []
    for pop in pops:
        train = [r for r in rows if r["population"] != pop]
        test = [r for r in rows if r["population"] == pop]
        require(train and test
                and not set(r["genotype"] for r in train).intersection(
                    r["genotype"] for r in test),
                "V46_ORIGINAL_GENOTYPE_LEAKAGE_ACROSS_FOLD")
        require(len(train) + len(test) == len(rows)
                and len({r["population"] for r in train}) == 19
                and len({r["population"] for r in test}) == 1,
                "V46_ORIGINAL_SOURCE_POPULATION_FOLD_DRIFT")
        folds.append((pop, train, test))
    return folds


def metric(value):
    require(math.isfinite(value), "V46_NONFINITE_RESULT")
    return round(value, 8)


def report(rows, audit, registry):
    folds = cohort_folds(rows)
    results = {e: [] for e in ACCEPTED_ENDPOINTS}
    for pop, train, test in folds:
        fitted = fit_train(train)
        taxon = test[0]["taxon"]
        require(all(r["taxon"] == taxon for r in test),
                "V46_ONE_SOURCE_POPULATION_SPANS_TAXA")
        for e in ACCEPTED_ENDPOINTS:
            mse = {}
            for model in MODELS:
                sse = sum((math.log1p(row[e]) -
                           prediction(fitted, row, model, e)) ** 2
                          for row in test)
                mse[model] = sse / len(test)
            results[e].append((pop, taxon, len(test), mse))
    output = {
        "analysis_type": "ORIGINAL_2024_SOURCE_POPULATION_HELDOUT_PREDICTION",
        "source_git_blob_sha1": registry["source_git_blob_sha1"],
        "source_commit": registry["source_commit"],
        "distinct_original_programs_added": 1,
        "independent_evolutionary_radiations_added": 0,
        "original_2014_Phlox_selection_source_used_for_model_fit": False,
        "measured_floral_color_or_pollinator_or_reproductive_alleles": False,
        "source_original_grain": audit,
        "holdout": "20 source populations, all genotype clones held out in three original gardens",
        "fixed_models": list(MODELS),
        "fixed_ridge_lambda": 10.0,
        "primary_fruit_transfer_pass_models": [],
        "nonadmitted_secondary_biomass": {
            "status": "BLOCKED_PREDECLARED_LOG1P_DOMAIN",
            "reason": "13 original signed biomass observations below zero; no clipping, deletion, shift or alternative transformation",
            "negative_source_rows": audit["biomass_original_negative_n"],
            "minimum_original_biomass": audit["biomass_original_minimum"],
            "numerical_scores_calculated": False
        },
        "endpoints": {}
    }
    for e, folds in results.items():
        baseline_eq = sum(m["B0"] for _, _, _, m in folds) / len(folds)
        baseline_wt = sum(n * m["B0"] for _, _, n, m in folds) / sum(
            n for _, _, n, _ in folds)
        summary = {"role": "PRIMARY" if e == "fruits_2019"
                   else "SECONDARY_EXPLORATORY",
                   "scale": "log1p_original_count", "model_results": {},
                   "per_population_holdout": []}
        for pop, taxon, n, model_mse in folds:
            summary["per_population_holdout"].append({
                "population": pop, "taxon": taxon, "n_ramets": n,
                "mse_by_model": {k: metric(model_mse[k]) for k in MODELS}})
        for m in MODELS:
            mse_eq = sum(x[3][m] for x in folds) / len(folds)
            mse_wt = sum(x[2] * x[3][m] for x in folds) / sum(
                x[2] for x in folds)
            bytax = {}
            for t in TAXA:
                subset = [x for x in folds if x[1] == t]
                subbase = sum(x[3]["B0"] for x in subset) / len(subset)
                submse = sum(x[3][m] for x in subset) / len(subset)
                bytax[t] = {
                    "n_source_populations": len(subset),
                    "relative_equal_population_mse_gain": metric(
                        1 - submse / subbase),
                    "positive_population_folds": sum(
                        1 for x in subset if
                        x[3][m] < x[3]["B0"] - 1e-12)}
            gain_eq = 1 - mse_eq / baseline_eq
            gain_wt = 1 - mse_wt / baseline_wt
            positives = sum(1 for x in folds if
                            x[3][m] < x[3]["B0"] - 1e-12)
            primary_pass = (e == "fruits_2019" and m != "B0"
                and gain_eq > registry["min_gain_equal_population"]
                and gain_wt > registry["min_gain_record_weighted"]
                and positives >= registry["min_positive_source_population_folds"]
                and all(bytax[t]["relative_equal_population_mse_gain"] > 0
                        for t in TAXA))
            summary["model_results"][m] = {
                "population_equal_mse": metric(mse_eq),
                "ramet_record_weighted_mse": metric(mse_wt),
                "relative_population_equal_gain": metric(gain_eq),
                "relative_record_weighted_gain": metric(gain_wt),
                "positive_population_folds": positives,
                "n_source_population_folds": len(folds),
                "by_taxon": bytax,
                "predeclared_primary_fruit_decision":
                    ("PASS_SOURCE_POP_TRANSFER" if primary_pass else
                     "NOT_ESTABLISHED") if e == "fruits_2019" else
                    "NOT_PRIMARY_ENDPOINT"}
            if primary_pass:
                output["primary_fruit_transfer_pass_models"].append(m)
        output["endpoints"][e] = summary
    return output


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--source", default=os.environ.get("PHLOX_V46_SOURCE"))
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    registry = source_registry()
    rows, source_audit = parse_source(acquire_raw(registry, args.source),
                                      registry)
    result = report(rows, source_audit, registry)
    if args.check:
        expected = json.loads(RECEIPT.read_text())
        require(expected == result, "V46_DERIVED_SOURCE_RESULT_RECEIPT_DRIFT")
        print("V46 ORIGINAL 2024 CSV HASH + ALL HOLDOUT PREDICTIONS + RECEIPT PASS")
        for e, summary in result["endpoints"].items():
            for m, v in summary["model_results"].items():
                print(e, m, "population_gain",
                    v["relative_population_equal_gain"],
                    "ramet_gain", v["relative_record_weighted_gain"],
                    "positive_folds", v["positive_population_folds"],
                    "decision", v["predeclared_primary_fruit_decision"])
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
