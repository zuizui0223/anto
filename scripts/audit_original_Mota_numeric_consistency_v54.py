#!/usr/bin/env python3
"""Original Mota 2022 sunflower source numeric consistency gate only.

Post-failed-v54 source inspection. Reports quality counts, never original
individual sunflower, field name, seed outcome, or visitor treatment effect.
Does NOT run prediction models, not a biological null.
"""
from __future__ import annotations
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_Mota_sunflower_field_holdout_v54 as v54

def fieldpair(rows):
    return {(r["Year"].strip(),r["Field_code"].strip()) for r in rows}

def parse(raw,key):
    r=v54.source_csv(raw,key)
    return r

def number(x):
    s=(x or "").strip()
    if s in ("","NA","NaN","nan",".","#VALUE!"):
        return None,"MISSING"
    try:
        v=float(s)
    except ValueError:
        return None,"UNKNOWN"
    if not math.isfinite(v):
        return None,"NONFINITE"
    if v<0:
        return None,"NEGATIVE"
    if abs(v-round(v))>1e-6:
        return None,"NONINTEGER"
    return v,"OK"

def examine():
    reg=v54.read_registry()
    v54.verified_anchor(reg)
    raw=v54.source_bytes(reg)
    p=parse(raw["Database_Productivity.csv"],(
        "Year","Field_code","No_Seeds_Total","No_Aborts_Total","No_rep_units"))
    v=parse(raw["Database_Visitation_rates.csv"],(
        "Year","Field_code","No_Interactions_TOTAL","No_Interactions_APIS",
        "No_Interactions_WILD","Total_no_monitored_units"))
    shared=fieldpair(p)&fieldpair(v)
    counters=Counter()
    states=defaultdict(Counter)
    max_abs_seed=max_abs_bee=0.
    for name,rows,fields in (
        ("harvest",p,("No_Seeds_Total","No_Aborts_Total","No_rep_units")),
        ("visits",v,("No_Interactions_TOTAL","No_Interactions_APIS",
                   "No_Interactions_WILD","Total_no_monitored_units"))):
        for r in rows:
            if (r["Year"].strip(),r["Field_code"].strip()) not in shared:
                continue
            counters[name+"_matched_rows"]+=1
            values=[]
            reasons=[]
            for f in fields:
                value,status=number(r[f])
                states[name+"_"+f][status]+=1
                values.append(value)
                reasons.append(status)
            if any(s!="OK" for s in reasons):
                counters[name+"_rows_any_numeric_not_admitted"]+=1
                continue
            if name=="harvest":
                d=values[0]+values[1]-values[2]
                max_abs_seed=max(max_abs_seed,abs(d))
                if abs(d)>1e-5:
                    counters["harvest_rows_seed_units_identity_disagree"]+=1
                if values[2]==0:
                    counters["harvest_rows_zero_reproductive_denominator"]+=1
            else:
                d=values[1]+values[2]-values[0]
                max_abs_bee=max(max_abs_bee,abs(d))
                if abs(d)>1e-5:
                    counters["visits_rows_total_vs_apis_wild_disagree"]+=1
                    counters["visits_discrepancy_direction_" +
                             ("APIS_WILD_EXCEEDS_TOTAL" if d>0 else
                              "TOTAL_EXCEEDS_APIS_WILD")]+=1
                if values[3]==0:
                    counters["visits_rows_zero_monitored_units"]+=1
    return {
        "analysis_type":"POST_GATE_ORIGINAL_MOTA_NUMERIC_SOURCE_INTEGRITY_DIAGNOSTIC_ONLY",
        "source":6671199,"original_file_md5":reg["original_files_md5"],
        "original_matched_fieldyears":len(shared),
        "original_matched_harvest_rows":len([1 for r in p if
                    (r["Year"].strip(),r["Field_code"].strip()) in shared]),
        "original_matched_visit_rows":len([1 for r in v if
                    (r["Year"].strip(),r["Field_code"].strip()) in shared]),
        "source_count_categories":dict(sorted(counters.items())),
        "per_source_numeric_column_statuses":{
            name:dict(sorted(c.items())) for name,c in sorted(states.items())},
        "maximum_absolute_author_seed_identity_difference":max_abs_seed,
        "maximum_absolute_author_visitation_identity_difference":max_abs_bee,
        "original_row_values_or_identifiers_exported":False,
        "predictive_effect_or_visitor_seed_association_computed":False,
        "original_v54_preregistered_source_identity_guard_modified":False,
        "explanation":"This count-only source audit follows v54 fail-closed source gate; not evidence of a zero biological effect."}

if __name__=="__main__":
    print("V54_ORIGINAL_MOTA_SOURCE_NUMERIC_GATE="+json.dumps(
        examine(),sort_keys=True,separators=(",",":")))
