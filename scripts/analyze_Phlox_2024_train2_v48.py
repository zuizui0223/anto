#!/usr/bin/env python3
"""Original Phlox 2024 two-source-population training benchmark (v48).

This is a POST-v46 diagnostic, not independent biological replication or
a preregistered causal selection test. All original author source records
are fetched ephemerally and their Git SHA-1 is checked by ANTO v46.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import analyze_original_Phlox_2024_population_holdout_v46 as v46

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data/external/anto_v48_original2024_train2_registry_20261119.json"
RECEIPT = ROOT / "results/anto_v48_phlox2024_train2_matched_source_20261119.json"
TAXA = ("amoena", "deamii", "pilosa")
GARDENS = ("1-Amoena-PB", "2-Deamii-JC", "3-Pilosa-LBL")
ENDPOINTS = ("fruits_2019", "flowers_2019")
MODELS = ("B0", "G", "E")


class V48Error(ValueError):
    pass


def require(ok, code):
    if not ok:
        raise V48Error(code)


def fixed_registry():
    d = json.loads(REGISTRY.read_text())
    require(
        d["analysis_type"] ==
            "POSTHOC_TWO_TRAIN_ONE_TEST_ORIGINAL_SOURCE_POPULATION_MATCHING"
        and d["original_author_repo"] == "PhloxHopkins/PhloxFieldAdaptation"
        and d["original_author_commit"] ==
            "d792fe58cce36721394adfd97e9f412396578e86"
        and d["original_csv_blob_sha1"] ==
            "4c780b0abe644ac69dae86a0450eb598c160626b"
        and d["original_csv_bytes"] == 175377
        and d["original_ramets"] == 963
        and d["original_genotypes"] == 284
        and d["original_populations"] == 20
        and d["original_gardens"] == list(GARDENS)
        and d["original_taxon_population_counts"] ==
            {"amoena": 8, "deamii": 3, "pilosa": 9}
        and d["models"] == list(MODELS)
        and d["primary_outcome"] == "fruits_2019"
        and d["secondary_outcome"] == "flowers_2019"
        and d["biomass_log1p_nonadmitted"] is True
        and d["negative_biomass_rows"] == 13
        and d["negative_biomass_minimum"] == -233
        and d["ridge_lambda"] == 10
        and d["train_populations_per_refit"] == 2
        and d["test_populations_per_refit"] == 1
        and d["panels_per_taxon"] ==
            {"amoena": 56, "deamii": 1, "pilosa": 84}
        and d["total_panels"] == 141
        and d["total_population_holdout_fits"] == 423
        and d["reference_panels"] == 140
        and d["source_family_replicates_added"] == 0
        and d["real_floral_colour_or_pollinator_data_admitted"] is False
        and d["biological_mechanism_inferred"] is False
        and d["valid_randomization_p_value"] is False
        and d["posthoc_after_v46_and_v47"] is True,
        "V48_FIXED_POSTHOC_REGISTRY_DRIFT",
    )
    return d


def load_author_source(source_path, registry):
    original_registry = v46.source_registry()
    require(
        original_registry["source_commit"] == registry["original_author_commit"]
        and original_registry["source_git_blob_sha1"] ==
            registry["original_csv_blob_sha1"]
        and original_registry["source_size"] == registry["original_csv_bytes"],
        "V48_PARENT_V46_AUTHOR_SOURCE_DRIFT",
    )
    rows, audit = v46.parse_source(
        v46.acquire_raw(original_registry, source_path), original_registry)
    require(
        audit["rows"] == 963 and audit["genotypes"] == 284
        and audit["populations"] == 20
        and audit["biomass_original_negative_n"] == 13
        and audit["biomass_original_minimum"] == -233,
        "V48_ORIGINAL_AUTHOR_COHORT_CHANGED",
    )
    return rows, audit


def source_panel_index(rows):
    index = defaultdict(list)
    genome_owner = {}
    for row in rows:
        k = (row["taxon"], row["population"])
        index[k].append(row)
        genotype = row["genotype"]
        if genotype in genome_owner:
            require(genome_owner[genotype] == k, "V48_CLONE_TRAIN_TEST_LEAKAGE")
        genome_owner[genotype] = k
    original = {t: sorted(p for k, p in index if k == t) for t in TAXA}
    require({t: len(original[t]) for t in TAXA} ==
            {"amoena": 8, "deamii": 3, "pilosa": 9},
            "V48_SOURCE_POPULATION_KEYS_CHANGED")
    return index, original


def build_train(train):
    require(train, "V48_EMPTY_TWO_SOURCE_POP_TRAINING")
    garden_means = {}
    for row in train:
        g = row["garden"]
        s = garden_means.setdefault(
            g, {"n": 0, "pc1": 0.0, "pc2": 0.0,
                "y": {e: 0.0 for e in ENDPOINTS}})
        s["n"] += 1
        s["pc1"] += row["leaf_traits_pc1"]
        s["pc2"] += row["leaf_traits_pc2"]
        for endpoint in ENDPOINTS:
            s["y"][endpoint] += math.log1p(row[endpoint])
    require(set(garden_means) == set(GARDENS),
            "V48_MISSING_ORIGINAL_GARDEN_FROM_TRAIN")
    for g in GARDENS:
        s = garden_means[g]
        require(s["n"] > 0, "V48_EMPTY_SOURCE_GARDEN")
        s["pc1"] /= s["n"]
        s["pc2"] /= s["n"]
        for endpoint in ENDPOINTS:
            s["y"][endpoint] /= s["n"]
    n = len(train)
    sd1 = math.sqrt(sum((r["leaf_traits_pc1"] -
                         garden_means[r["garden"]]["pc1"]) ** 2
                        for r in train) / n)
    sd2 = math.sqrt(sum((r["leaf_traits_pc2"] -
                         garden_means[r["garden"]]["pc2"]) ** 2
                        for r in train) / n)
    require(sd1 > 1e-10 and sd2 > 1e-10,
            "V48_TWO_POP_TRAIN_FEATURES_NOT_IDENTIFIABLE")
    fit = {"garden_means": garden_means, "sd1": sd1, "sd2": sd2,
           "beta": {}}
    for model, p in (("G", 2), ("E", 6)):
        lhs = [[10.0 if i == j else 0.0 for j in range(p)]
               for i in range(p)]
        rhs = {e: [0.0] * p for e in ENDPOINTS}
        for row in train:
            m = garden_means[row["garden"]]
            z1 = (row["leaf_traits_pc1"] - m["pc1"]) / sd1
            z2 = (row["leaf_traits_pc2"] - m["pc2"]) / sd2
            off = 0 if model == "G" else 2 * GARDENS.index(row["garden"])
            lhs[off][off] += z1 * z1
            lhs[off][off + 1] += z1 * z2
            lhs[off + 1][off] += z2 * z1
            lhs[off + 1][off + 1] += z2 * z2
            for endpoint in ENDPOINTS:
                resid = math.log1p(row[endpoint]) - m["y"][endpoint]
                rhs[endpoint][off] += z1 * resid
                rhs[endpoint][off + 1] += z2 * resid
        fit["beta"][model] = {e: v46.solve(lhs, rhs[e])
                              for e in ENDPOINTS}
    return fit


def predict(fit, row, model, endpoint):
    m = fit["garden_means"][row["garden"]]
    base = m["y"][endpoint]
    if model == "B0":
        return base
    require(model in ("G", "E"), "V48_UNREGISTERED_SLOPE_MODEL")
    z1 = (row["leaf_traits_pc1"] - m["pc1"]) / fit["sd1"]
    z2 = (row["leaf_traits_pc2"] - m["pc2"]) / fit["sd2"]
    off = 0 if model == "G" else 2 * GARDENS.index(row["garden"])
    beta = fit["beta"][model][endpoint]
    return base + z1 * beta[off] + z2 * beta[off + 1]


def score_fold(train, test, train_pop_ids, heldout_pop_id, taxon):
    require(len(set(train_pop_ids)) == 2
            and set(r["population"] for r in train) == set(train_pop_ids)
            and set(r["population"] for r in test) == {heldout_pop_id}
            and heldout_pop_id not in train_pop_ids
            and all(r["taxon"] == taxon for r in train + test),
            "V48_WRONG_TAXON_OR_TWO_TRAIN_POP_HOLDOUT")
    tgen = {r["genotype"] for r in train}
    qgen = {r["genotype"] for r in test}
    require(tgen.isdisjoint(qgen), "V48_CLONE_LEAKS_ACROSS_TRAIN_TEST")
    require(set(r["garden"] for r in train) == set(GARDENS)
            and set(r["garden"] for r in test) == set(GARDENS),
            "V48_NOT_ALL_THREE_COMMON_GARDENS_TESTED")
    fit = build_train(train)
    scores = {}
    for e in ENDPOINTS:
        scores[e] = {}
        for m in MODELS:
            err = sum((math.log1p(r[e]) -
                       predict(fit, r, m, e)) ** 2 for r in test)
            scores[e][m] = err / len(test)
            require(math.isfinite(scores[e][m]) and scores[e][m] >= 0,
                    "V48_BAD_PREDICTION_LOSS")
    return scores


def gain(folds, endpoint, model):
    require(folds and model != "B0", "V48_INVALID_GAIN_INPUT")
    base = sum(f["loss"][endpoint]["B0"] for f in folds)
    model_loss = sum(f["loss"][endpoint][model] for f in folds)
    require(base > 0 and math.isfinite(model_loss),
            "V48_INVALID_BASELINE_LOSS")
    return 1 - model_loss / base


def num(v):
    require(math.isfinite(v), "V48_NONFINITE_RECEIPT_RESULT")
    return round(v, 9)


def summarize(vals, deamii):
    v = sorted(vals)
    require(v, "V48_EMPTY_REFERENCE")
    n = len(v)
    return {
        "n_overlapping_three_pop_source_panels": n,
        "n_negative_panel_gains": sum(x < 0 for x in v),
        "n_at_or_below_deamii": sum(x <= deamii + 1e-10 for x in v),
        "min_gain": num(v[0]),
        "median_gain": num(statistics.median(v)),
        "max_gain": num(v[-1]),
        "valid_randomization_p_value": False,
        "independent_biological_experiments": 0,
    }


def analyze(rows, audit, registry):
    index, pops = source_panel_index(rows)
    panels = {t: [] for t in TAXA}
    raw_gains = {e: {m: {t: [] for t in TAXA}
                     for m in ("G", "E")} for e in ENDPOINTS}
    counter = 0
    for taxon in TAXA:
        for panel_pops in itertools.combinations(pops[taxon], 3):
            panel = {
                "taxon": taxon,
                "three_source_populations": list(panel_pops),
                "three_folds_one_source_pop_out": [],
                "equal_source_pop_panel_relative_gain": {}
            }
            for heldout in panel_pops:
                train_pops = [p for p in panel_pops if p != heldout]
                train = [r for p in train_pops for r in index[(taxon, p)]]
                test = index[(taxon, heldout)]
                outcome = score_fold(train, test, train_pops, heldout, taxon)
                panel["three_folds_one_source_pop_out"].append({
                    "training_source_populations": train_pops,
                    "heldout_source_population": heldout,
                    "n_training_ramets": len(train),
                    "n_heldout_ramets": len(test),
                    "loss": outcome,
                })
                counter += 1
            for e in ENDPOINTS:
                panel["equal_source_pop_panel_relative_gain"][e] = {}
                for m in ("G", "E"):
                    raw = gain(panel["three_folds_one_source_pop_out"], e, m)
                    raw_gains[e][m][taxon].append(raw)
                    panel["equal_source_pop_panel_relative_gain"][e][m] = num(raw)
            panels[taxon].append(panel)
    require(counter == 423
            and {t: len(panels[t]) for t in TAXA} ==
                 {"amoena": 56, "deamii": 1, "pilosa": 84},
            "V48_EXPECT_141_PANELS_423_COMPLETE_FITS")
    comparison = {}
    for e in ENDPOINTS:
        comparison[e] = {}
        for m in ("G", "E"):
            d = raw_gains[e][m]["deamii"][0]
            am = raw_gains[e][m]["amoena"]
            pi = raw_gains[e][m]["pilosa"]
            comparison[e][m] = {
                "deamii_one_panel_relative_gain": num(d),
                "deamii_original_source_populations": pops["deamii"],
                "amoena_same_train2_reference": summarize(am, d),
                "pilosa_same_train2_reference": summarize(pi, d),
                "combined_140_overlapping_reference": summarize(am + pi, d),
            }
    return {
        "analysis_type": registry["analysis_type"],
        "source_original_author_commit": registry["original_author_commit"],
        "source_original_CSV_git_blob_sha1": registry["original_csv_blob_sha1"],
        "original_source_grain": audit,
        "models": list(MODELS),
        "endpoints": list(ENDPOINTS),
        "fit_scope": "two same-taxon training source populations, one complete heldout source population, all three original gardens",
        "new_source_model_fits": 423,
        "all_panels": [x for t in TAXA for x in panels[t]],
        "same_train_two_pop_comparisons": comparison,
        "original_biomass_log1p_technical_nonadmission": True,
        "original_v46_all_taxa_primary_transfer_remains_NOT_ESTABLISHED": True,
        "new_independent_field_populations_collected": 0,
        "new_original_studies": 0,
        "new_flower_colour_or_pollinator_experiments": 0,
        "valid_randomization_p_values": False,
        "is_posthoc_model_refit": True
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", default=None)
    p.add_argument("--check", action="store_true")
    args = p.parse_args()
    registry = fixed_registry()
    rows, audit = load_author_source(args.source, registry)
    output = analyze(rows, audit, registry)
    if args.check:
        require(json.loads(RECEIPT.read_text()) == output,
                "V48_EXACT_SOURCE_REGENERATED_RECEIPT_DIFFERS")
        print("V48_ORIGINAL_AUTHOR_SHA_AND_ALL_423_FITS_EXACT_PASS")
    else:
        source_json = json.dumps(output, sort_keys=True,
                                 separators=(",", ":"), ensure_ascii=True)
        n_chunks = (len(source_json) + 3999) // 4000
        for k in range(n_chunks):
            part = source_json[k*4000:(k+1)*4000]
            print("V48_SOURCE_READOUT_CHUNK_%03d_OF_%03d=" %
                  (k+1, n_chunks) + part)


if __name__ == "__main__":
    main()
