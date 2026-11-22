#!/usr/bin/env python3
"""v55 original Mota sunflower source stop: successful seed / total denominator.

Original Zenodo MD5 and v53/v54 source ledger pinned. No visitor→seed
model, correlations, coefficients or per-original-field rows disclosed.
"""
from __future__ import annotations
import json
import sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mota_reported_total_bees_v55 as v55
import analyze_original_Mota_sunflower_field_holdout_v54 as v54

def examine():
    r=v55.registry()
    v55.require_prior_quality(r)
    raw=v55.original_source(r)
    h=v54.source_csv(raw["Database_Productivity.csv"],(
        "Year","Field_code","No_Seeds_Total","No_rep_units","No_Aborts_Total"))
    b=v54.source_csv(raw["Database_Visitation_rates.csv"],(
        "Year","Field_code"))
    hkeys={v54.fieldyear_key(z) for z in h}
    bkeys={v54.fieldyear_key(z) for z in b}
    common=hkeys&bkeys
    c=Counter()
    category_by_year=Counter()
    # Classification for 2217 matched harvest units only. Source successful
    # seeds and original denominator are never exported, only reasons.
    for z in h:
        key=v54.fieldyear_key(z)
        if key not in common: continue
        c["total_matched_harvest_source_rows"]+=1
        a=v54.parse_number(z["No_Seeds_Total"])
        d=v54.parse_number(z["No_rep_units"])
        if a is None or d is None:
            reason="missing_numerator_or_denominator"
        elif d==0:
            reason="zero_denominator"
        elif a>d:
            reason="successful_seeds_greater_than_reproductive_units"
        else:
            reason="valid_source_fraction"
        c[reason]+=1
        category_by_year[(key[0],reason)]+=1
    assert c["total_matched_harvest_source_rows"]==2217
    return {
        "analysis_type":"POST_V55_SOURCE_SEED_RATIO_QUALITY_ONLY_NOT_PREDICTION",
        "source_original_record":6671199,
        "source_original_file_md5":r["original_file_md5"],
        "v54_source_STOP_preserved":True,
        "new_v55_source_invalid_seed_fraction_not_imputed":True,
        "source_original_harvest_rows":2217,
        "source_ratio_domain_status":dict(sorted(c.items())),
        "source_ratio_domain_status_by_year":[
            {"year":year,"status":label,"n":n}
            for (year,label),n in sorted(category_by_year.items())],
        "original_seed_or_visitor_outcome_prediction_fitted":False,
        "no_source_plant_or_field_identifier_exported":True,
    }

if __name__=="__main__":
    print("V55_AUTHOR_SOURCE_SEED_DENOMINATOR_QUALITY="+json.dumps(
        examine(),sort_keys=True,separators=(",",":")))
