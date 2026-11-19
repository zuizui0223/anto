#!/usr/bin/env python3
"""Mota original sunflower archive: field×year identity / codebook ONLY.

MD5-authenticated source, no visitation/seed numerical outcomes examined.
"""
from __future__ import annotations
import csv
import io
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_independent_pollinator_reproduction_sources_v53 as base

def read_csv(raw):
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    return list(reader)

def field_keys(rows):
    keys=Counter()
    miss=0
    for r in rows:
        year=r["Year"].strip()
        field=r["Field_code"].strip()
        if not year or not field:
            miss+=1
            continue
        keys[(year,field)]+=1
    return keys,miss

def audit():
    reg=base.registry()
    entry=next(x for x in reg["candidates"] if x["code"]=="Mota2022")
    files=base.fetch_original(entry,reg["max_bytes_per_file"])
    product=read_csv(files["Database_Productivity.csv"])
    visits=read_csv(files["Database_Visitation_rates.csv"])
    kp,missp=field_keys(product)
    kv,missv=field_keys(visits)
    both=set(kp)&set(kv)
    def profile(values):
        counts=sorted(values.values())
        return {"n_original_distinct_field_years":len(counts),
                "total_source_rows":sum(counts),
                "minimum_records_per_field_year":counts[0] if counts else None,
                "median_records_per_field_year":statistics.median(counts) if counts else None,
                "maximum_records_per_field_year":counts[-1] if counts else None}
    candidates={}
    for label,values in (("productivity",kp),("visitation",kv)):
        candidates[label]={"original_field_year_profile":profile(values),
                           "year_counts_original_field_years":dict(sorted(
                                Counter(key[0] for key in values).items()))}
    text=files["README_file.txt"].decode("utf-8-sig",errors="replace")
    # Source README is used only to check field and year semantics; the full
    # author note is not redistributed with any source records in GitHub.
    print("V53_MOTA_README_START")
    print(text)
    print("V53_MOTA_README_END")
    return {"zenodo_record":6671199,
        "source_original_files_MD5":entry["files_md5"],
        "source_cross_file_common_column_names":["Year","Field_code"],
        "source_identity_units":"FIELD×YEAR ONLY, not same individual plants",
        "productivity":candidates["productivity"],
        "visitation":candidates["visitation"],
        "n_productivity_blank_field_year":missp,
        "n_visitation_blank_field_year":missv,
        "n_shared_field_year_units":len(both),
        "n_shared_field_year_units_by_year":dict(sorted(Counter(
            x[0] for x in both).items())),
        "n_productivity_records_on_shared_field_year":
            sum(kp[k] for k in both),
        "n_visitation_records_on_shared_field_year":
            sum(kv[k] for k in both),
        "same_plant_id_join_identifiable":False,
        "source_original_seed_yield_or_visit_count_values_opened":False,
        "no_predictive_model_fitted":True
    }

if __name__=="__main__":
    print("V53_MOTA_SOURCE_KEY_AUDIT="+json.dumps(audit(),sort_keys=True))
