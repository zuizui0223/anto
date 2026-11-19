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
    text=raw.decode("utf-8-sig")
    try:
        delim=csv.Sniffer().sniff(text[:16000],delimiters=",;").delimiter
    except csv.Error:
        delim=","
    reader=csv.DictReader(io.StringIO(text),delimiter=delim)
    assert reader.fieldnames is not None and "Year" in reader.fieldnames and "Field_code" in reader.fieldnames, "V53_SOURCE_FIELD_YEAR_COLUMNS_MISSING"
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
    # Independent source field-year and physical field identity gate.
    from collections import defaultdict
    meta_p=defaultdict(lambda:{"names":set(),"region":set(),"field_type":set(),
                               "distance":set()})
    meta_v=defaultdict(lambda:{"names":set(),"region":set(),"field_type":set(),
                               "distance":set()})
    for source,metas,namecol in (
        (product,meta_p,"Field_code_original"),
        (visits,meta_v,"Field_ID")):
        for row in source:
            k=(row["Year"].strip(),row["Field_code"].strip())
            q=metas[k]
            q["names"].add(row[namecol].strip())
            q["region"].add(row["Region"].strip())
            q["field_type"].add(row["Field_Type"].strip())
            q["distance"].add(row["Distance"].strip())
    by_type=Counter()
    by_region=Counter()
    n_same_orig_name=0
    n_exact_single_name=0
    n_single_meta=0
    band_counts=Counter()
    distance_intersection=Counter()
    physical_fields=set()
    for k in both:
        p=meta_p[k];w=meta_v[k]
        n_exact_single_name+=(len(p["names"])==len(w["names"])==1)
        n_same_orig_name+=(len(p["names"])==len(w["names"])==1
                           and p["names"]==w["names"])
        all_same=(len(p["region"])==len(w["region"])==1 and
                  len(p["field_type"])==len(w["field_type"])==1 and
                  p["region"]==w["region"] and
                  p["field_type"]==w["field_type"])
        n_single_meta+=bool(all_same)
        if all_same:
            by_type[next(iter(p["field_type"]))]+=1
            by_region[next(iter(p["region"]))]+=1
        inter=p["distance"]&w["distance"]
        distance_intersection[len(inter)]+=1
        band_counts[len(p["distance"])]+=1
        if len(p["names"])==1 and len(p["region"])==1:
            physical_fields.add((next(iter(p["region"])),next(iter(p["names"]))))
    identity_summary={
        "shared_fieldyear_original_name_singleton_in_both":n_exact_single_name,
        "shared_fieldyear_same_original_physical_field_name":n_same_orig_name,
        "shared_fieldyear_region_fieldtype_consistent_and_singleton":n_single_meta,
        "shared_fieldyear_by_original_field_type":dict(sorted(by_type.items())),
        "shared_fieldyear_by_original_region":dict(sorted(by_region.items())),
        "shared_fieldyear_shared_distance_band_count":dict(sorted(distance_intersection.items())),
        "productivity_fieldyear_distance_band_count":dict(sorted(band_counts.items())),
        "distinct_original_region_plus_field_names_in_shared_fieldyears":len(physical_fields),
        "seed_or_visitation_response_data_examined":False,
    }
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
        "physical_field_metadata_join_audit":identity_summary,
        "same_plant_id_join_identifiable":False,
        "source_original_seed_yield_or_visit_count_values_opened":False,
        "no_predictive_model_fitted":True
    }

if __name__=="__main__":
    print("V53_MOTA_SOURCE_KEY_AUDIT="+json.dumps(audit(),sort_keys=True))
