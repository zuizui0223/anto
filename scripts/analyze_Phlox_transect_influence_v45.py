#!/usr/bin/env python3
"""Source-frozen v45 descriptive Phlox transect influence audit.

This reads the committed ANTO v44 source-C++ receipt; it does not run
new author MCMC, optimize parameters, or claim biological replication.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/anto_v44_original_CPP_seven_scenarios_likelihood_20261117.json"
FROZEN = ROOT / "results/anto_v45_phlox_transect_influence_20261118.json"
SOURCE_GIT_BLOB_SHA1 = "92d940d60051f57bad023e290e0dd5fbff5b9822"
ADMITTED = ("baseline", "fitness_neutral", "epistasis_zero",
            "dispersal_half", "dispersal_double")
BLOCKED = ("fitness_uniform_zone0", "assortative_25pct")
THRESHOLD = 1e-7


class AuditError(ValueError):
    pass


def require(ok: bool, code: str) -> None:
    if not ok:
        raise AuditError(code)


def sign(value: float) -> int:
    require(math.isfinite(value), "NONFINITE_COMPONENT")
    return 1 if value > THRESHOLD else (-1 if value < -THRESHOLD else 0)


def git_blob_sha1(content: bytes) -> str:
    prefix = ("blob " + str(len(content)) + "\0").encode("ascii")
    return hashlib.sha1(prefix + content).hexdigest()


def read_source(path: Path = SOURCE) -> dict:
    raw = path.read_bytes()
    require(git_blob_sha1(raw) == SOURCE_GIT_BLOB_SHA1,
            "V44_FROZEN_PARENT_GIT_BLOB_CHANGED")
    source = json.loads(raw)
    validate(source)
    return source


def validate(source: dict) -> None:
    require(source.get("analysis_type") ==
            "EXPLORATORY_FIXED_START_SOURCE_CPP_SENSITIVITY_NOT_MODEL_SELECTION",
            "V44_ANALYSIS_TYPE_CHANGED")
    for key, expected in (
        ("n_independent_original_phylogenetic_studies_added", 0),
        ("n_source_scenarios_admitted", 5),
        ("n_source_scenarios_technically_blocked", 2),
        ("original_2014_75000_step_MCMC_performed", False),
        ("original_2014_published_full_12024_382_2581_used", False),
        ("new_ancestral_flower_colour_regain_or_allele_frequency_estimated", False),
        ("original_source_32p_11379adults_371families_2496offspring", True),
        ("no_original_source_raw_rows_redistributed", True),
    ):
        require(source.get(key) == expected, "V44_CLAIM_OR_DENOMINATOR_DRIFT_" + key)
    expected = set(ADMITTED + BLOCKED)
    require(set(source.get("scenarios", {})) == expected,
            "V44_SEVEN_SCENARIOS_CHANGED")
    require(set(source.get("scenario_names_admitted", [])) == set(ADMITTED)
            and set(source.get("scenario_names_blocked", [])) == set(BLOCKED),
            "V44_ADMISSION_LIST_DRIFT")
    for name in BLOCKED:
        row = source["scenarios"][name]
        require(row.get("status") == "BLOCKED_NUMERICAL_OR_ORIGINAL_CODE"
                and "comparison" not in row
                and "failure_category" in row
                and row.get("numerical_failure_not_biological_zero") is True,
                "V44_TECHNICAL_BLOCK_TREATED_AS_BIOLOGY")
    baseline = source["scenarios"]["baseline"]["output"]["all_transects"]
    require((baseline["n_adults"], baseline["n_mothers"],
             baseline["n_offspring"], baseline["n_populations"]) ==
            (11379, 371, 2496, 32), "V44_ORIGINAL_EXAMPLE_DENOMINATORS_DRIFT")
    for name in ADMITTED:
        row = source["scenarios"][name]
        require(row.get("status") == "ADMITTED_ORIGINAL_FORWARD_SOURCE"
                and row.get("source_true_change_tolerance_verified") is True,
                "V44_SCENARIO_CONVERGENCE_NOT_ADMITTED")
        comparison = row.get("comparison", {})
        totals = comparison.get("all_transects_delta_logL", {})
        transects = comparison.get("per_transect_delta_logL", [])
        require(set(totals) == {"adult", "offspring", "total"}
                and len(transects) == 5
                and [v.get("transect") for v in transects] == [1, 2, 3, 4, 5],
                "V44_TRANSECT_GRAIN_CHANGED")
        for key in ("adult", "offspring", "total"):
            require(math.isfinite(totals[key]), "V44_NONFINITE_AGGREGATE")
            col = "Delta_logL_" + key + "_vs_original_START"
            vals = [v[col] for v in transects]
            require(all(math.isfinite(v) for v in vals) and
                    math.isclose(sum(vals), totals[key], rel_tol=0, abs_tol=1e-6),
                    "V44_TRANSECT_SUM_MISMATCH_" + key)
            require(sum(sign(v) > 0 for v in vals) ==
                    comparison["n_transects_delta_positive_by_component"][key]
                    and sum(sign(v) < 0 for v in vals) ==
                    comparison["n_transects_delta_negative_by_component"][key],
                    "V44_COMPONENT_SIGN_COUNT_MISMATCH_" + key)
        require(math.isclose(totals["adult"] + totals["offspring"],
                             totals["total"], rel_tol=0, abs_tol=1e-6),
                "V44_LIKELIHOOD_COMPONENTS_FAIL_TO_ADD")
        for v in transects:
            require(math.isclose(v["Delta_logL_adult_vs_original_START"] +
                                 v["Delta_logL_offspring_vs_original_START"],
                                 v["Delta_logL_total_vs_original_START"],
                                 rel_tol=0, abs_tol=1e-6),
                    "V44_TRANSECT_COMPONENTS_FAIL_TO_ADD")
        if name == "baseline":
            require(all(sign(totals[k]) == 0
                        for k in ("adult", "offspring", "total")),
                    "V44_BASELINE_NOT_ZERO")


def report(source: dict) -> dict:
    validate(source)
    result = {
        "analysis_type": "POSTHOC_SOURCE_FROZEN_TRANSECT_INFLUENCE_NOT_INFERENCE",
        "v44_parent_git_blob_sha1": SOURCE_GIT_BLOB_SHA1,
        "original_programs_added": 0,
        "original_2014_full_paper_cohort_reconciled": False,
        "new_forward_model_executions": 0,
        "new_optimization_or_MCMC": False,
        "n_transects": 5,
        "n_adult_records_example": 11379,
        "n_mother_families_example": 371,
        "n_offspring_example": 2496,
        "source_technical_nonadmissions": {
            name: source["scenarios"][name]["failure_category"]
            for name in BLOCKED
        },
        "scenarios": {}
    }
    for name in ADMITTED[1:]:
        original = source["scenarios"][name]["comparison"]
        sums = original["all_transects_delta_logL"]
        rows = [
            {
                "transect": v["transect"],
                "adult": v["Delta_logL_adult_vs_original_START"],
                "offspring": v["Delta_logL_offspring_vs_original_START"],
                "total": v["Delta_logL_total_vs_original_START"],
                "adult_offspring_opposite": (
                    sign(v["Delta_logL_adult_vs_original_START"]) *
                    sign(v["Delta_logL_offspring_vs_original_START"]) == -1
                ),
            }
            for v in original["per_transect_delta_logL"]
        ]
        omitted = []
        for v in rows:
            omitted.append({
                "omitted_transect": v["transect"],
                "adult": sums["adult"] - v["adult"],
                "offspring": sums["offspring"] - v["offspring"],
                "total": sums["total"] - v["total"],
            })
        sign_changes = {
            key: [v["omitted_transect"] for v in omitted
                  if sign(v[key]) != sign(sums[key])]
            for key in ("adult", "offspring", "total")
        }
        result["scenarios"][name] = {
            "global_delta_logL": sums,
            "global_adult_offspring_opposite": (
                sign(sums["adult"]) * sign(sums["offspring"]) == -1
            ),
            "source_per_transect": rows,
            "within_transect_opposite_count": sum(
                x["adult_offspring_opposite"] for x in rows),
            "delete_one_transect_remaining_delta_logL": omitted,
            "sign_reversals_under_one_transect_deletion": sign_changes,
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = report(read_source())
    if args.check:
        require(json.loads(FROZEN.read_text()) == output,
                "V45_STORED_RECEIPT_REGENERATION_MISMATCH")
        print("V45_SOURCE_BLOB_AND_STORED_REAGGREGATION_PASS")
    else:
        print(json.dumps(output, indent=2, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
