#!/usr/bin/env python3
"""ANTO v56 — source-locked Mota sunflower 2214 valid harvest rows.

POST-v55 source STOP design, preregistered before v56 predicted outcomes.
Only exactly three 2017 author rows with full seeds > total reproduction
units are excluded. Reuses unchanged v55 whole-field Z0/Z1/Z2 ridge
and the 200 observational visit-label comparisons.
Not same-plant visitor–seed, not experimental pollination, no p values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import analyze_original_Mota_reported_total_bees_v55 as v55

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v56_valid_seed_fraction_registry_20261122.json"
RECEIPT=ROOT/"results/anto_v56_sunflower_valid_fraction_field_transfer_20261122.json"
PRIOR_V55_CODE=ROOT/"scripts/analyze_original_Mota_reported_total_bees_v55.py"
BANDS=("0","15","30","60")
YEARS=("2017","2018")
REGIONS=("Burgos","Cuenca")
TYPES=("NonV","SNV","WFS")


class V56Error(ValueError):
    pass


def require(ok,code):
    if not ok:
        raise V56Error(code)


def git_blob(raw):
    return hashlib.sha1(("blob %d\0"%len(raw)).encode()+raw).hexdigest()


def registered():
    r=json.loads(REGISTRY.read_text())
    require(
        r["analysis_type"]==
           "POST_V55_ORIGINAL_SUNFLOWER_EXACT_THREE_INVALID_ROWS_EXCLUDED_FIELD_TRANSFER"
        and r["source_zenodo_record"]==6671199
        and r["source_original_md5"]=={
           "Database_Productivity.csv":"aeea30a21e495498295241bcfd8c693d",
           "Database_Visitation_rates.csv":"e13131ad60b133cc6ab25e20bbd3554b",
           "README_file.txt":"479b26d3cebb72c8491b0ec00609c08c"}
        and r["v53_fieldyear_source_blob_sha1"]==
           "7a4fb27c7f26abfa8400d23b041e1a3d3e3424c9"
        and r["v54_original_quality_blob_sha1"]==
           "55eb9b942e7a0ab7f4e32139ac5305325fb74a83"
        and r["v55_invalid_fraction_quality_blob_sha1"]==
           "f9d6c01720ae0657977d09d0a92043b0eff0d7b5"
        and r["v55_original_model_implementation_git_blob_sha1"]==
           "3f4b936101ef51e56e1fc661089d49b0efe8810b"
        and r["original_2017_2018_matched_fieldyears"]==80
        and r["original_unique_physical_fields"]==51
        and r["original_matched_harvest_source_rows"]==2217
        and r["original_matched_bee_visitation_survey_rows"]==11682
        and r["source_invalid_seed_fraction_above_one_rows_excluded"]==3
        and r["source_invalid_above_one_rows_year"]=={"2017":3}
        and r["source_included_harvest_rows"]==2214
        and r["source_missing_abort_rows_preserved_not_imputed"]==1056
        and r["source_guild_component_disagreements_preserved"]==4087
        and r["source_max_guild_total_difference"]==9
        and r["source_original_years"]==list(YEARS)
        and r["source_original_regions"]==list(REGIONS)
        and r["source_original_fieldtypes"]==list(TYPES)
        and r["source_2017_2018_fieldyears"]=={"2017":41,"2018":39}
        and r["source_fieldyears_by_region"]=={"Burgos":38,"Cuenca":42}
        and r["source_fieldyears_by_type"]=={"NonV":10,"SNV":32,"WFS":38}
        and r["source_original_distance_bands"]==list(BANDS)
        and r["original_all_common_distance_band_counts"]=={"3":1,"4":79}
        and r["models"]==["Z0","Z1","Z2"]
        and r["ridge_lambda"]==8
        and r["n_physical_field_holdouts"]==51
        and r["n_fieldyear_heldouts"]==80
        and r["source_placebo_seed"]==20261119
        and r["source_placebo_replicates"]==200
        and r["v54_source_stop_cannot_be_redefined"] is True
        and r["v55_source_stop_cannot_be_redefined"] is True
        and r["new_original_independent_biological_studies"]==0
        and r["new_molecular_pollinator_evolution_evidence"] is False
        and r["individual_sunflower_visit_seed_join_possible"] is False
        and r["ecological_causal_p_value"] is False
        and r["posthoc_after_v55_source_failure"] is True,
        "V56_PREOUTCOME_SOURCE_PROTOCOL_OR_CLAIM_CHANGED")
    require(git_blob(PRIOR_V55_CODE.read_bytes())==
            r["v55_original_model_implementation_git_blob_sha1"],
            "V56_UNCHANGED_REGISTERED_V55_MODEL_SOURCE_MOVED")
    return r


def original_denominator_quality(r):
    path=ROOT/"results/anto_v55_Mota_seed_fraction_source_quality_20261119.json"
    raw=path.read_bytes()
    require(git_blob(raw)==r["v55_invalid_fraction_quality_blob_sha1"],
            "V56_V55_DENOMINATOR_AUDIT_PARENT_BYTES_CHANGED")
    quality=json.loads(raw)
    x=quality["source_ratio_domain_status"]
    require(quality["source_original_record"]==6671199
            and quality["original_source_file_md5"]==r["source_original_md5"]
            and quality["v54_source_STOP_preserved"] is True
            and quality["original_seed_or_visitor_outcome_prediction_fitted"] is False
            and x["total_matched_harvest_source_rows"]==2217
            and x["valid_source_fraction"]==2214
            and x["successful_seeds_greater_than_reproductive_units"]==3
            and not set(x).difference({
              "total_matched_harvest_source_rows",
              "valid_source_fraction",
              "successful_seeds_greater_than_reproductive_units"}),
            "V56_HISTORIC_V55_SOURCE_QUALITY_REINTERPRETED")
    return quality


def validate_source_archives(r):
    old=v55.registry()
    require(old["original_file_md5"]==r["source_original_md5"]
            and old["original_v53_fieldyear_blob_sha1"]==
                r["v53_fieldyear_source_blob_sha1"]
            and old["original_v54_quality_blob_sha1"]==
                r["v54_original_quality_blob_sha1"],
            "V56_PARENT_V55_RECORDED_SOURCE_MOVED")
    v55.require_prior_quality(old)
    return v55.original_source(old)


def aggregate_original(raw,reg):
    src=v55.old
    harvest=src.source_csv(raw["Database_Productivity.csv"],(
        "Region","Year","Field_code","Field_Type","Field_code_original",
        "Distance","No_Seeds_Total","No_rep_units","No_Aborts_Total"))
    visits=src.source_csv(raw["Database_Visitation_rates.csv"],(
        "Region","Year","Field_code","Field_Type","Field_ID",
        "Distance","No_Interactions_TOTAL","No_Interactions_APIS",
        "No_Interactions_WILD","Total_no_monitored_units"))
    require(len(harvest)==2875 and len(visits)==11683,
            "V56_ORIGINAL_AUTHOR_TABLE_LENGTH_DRIFT")
    hk=Counter(src.fieldyear_key(z) for z in harvest)
    vk=Counter(src.fieldyear_key(z) for z in visits)
    shared=set(hk)&set(vk)
    require(len(shared)==80
            and sum(hk[k] for k in shared)==2217
            and sum(vk[k] for k in shared)==11682,
            "V56_V53_ORIGINAL_JOIN_DRIFT")

    meta_h=defaultdict(lambda:{"region":set(),"type":set(),
                               "physical":set(),"band":set()})
    meta_v=defaultdict(lambda:{"region":set(),"type":set(),
                               "physical":set(),"band":set()})
    seed=defaultdict(lambda:{"numerator":0.0,"denominator":0.0,"n":0})
    bee=defaultdict(lambda:{"reported":0.0,"reconstructed":0.0,
                            "exposure":0.0,"n":0})
    invalid=Counter()
    abort_missing=triads=bee_disagree=0
    max_bee_disagreement=0.0
    seed_count=visit_count=0
    for z in harvest:
        fieldyear=src.fieldyear_key(z)
        if fieldyear not in shared:continue
        region,kind,field=src.sanitize_metadata(z,"harvest")
        band=src.distance_key(z)
        m=meta_h[fieldyear]
        m["region"].add(region);m["type"].add(kind)
        m["physical"].add(field);m["band"].add(band)
        full=src.parse_number(z["No_Seeds_Total"])
        total=src.parse_number(z["No_rep_units"])
        abort=src.parse_number(z["No_Aborts_Total"])
        require(full is not None and total is not None and total>0,
                "V56_NEW_UNEXPECTED_AUTHOR_SEED_COUNT_MISSING_OR_ZERO")
        if abort is None:
            abort_missing+=1
        else:
            triads+=1
            require(abs(full+abort-total)<=1e-6,
                    "V56_ORIGINAL_COMPLETE_ABORT_TRIAD_MISMATCH")
        if full>total:
            invalid[fieldyear[0]]+=1
            # This source-invalid row is explicitly excluded BEFORE
            # source outcome aggregation, never replaced by zero/one.
            continue
        require(0<=full<=total,"V56_UNRECOGNIZED_SEED_FRACTION_DOMAIN")
        s=seed[(fieldyear,band)]
        s["numerator"]+=full;s["denominator"]+=total;s["n"]+=1
        seed_count+=1
    for z in visits:
        fieldyear=src.fieldyear_key(z)
        if fieldyear not in shared:continue
        region,kind,field=src.sanitize_metadata(z,"visits")
        band=src.distance_key(z)
        m=meta_v[fieldyear]
        m["region"].add(region);m["type"].add(kind)
        m["physical"].add(field);m["band"].add(band)
        numbers=[src.parse_number(z[c]) for c in (
            "No_Interactions_TOTAL","No_Interactions_APIS",
            "No_Interactions_WILD","Total_no_monitored_units")]
        require(all(n is not None for n in numbers)
                and numbers[3]>0,
                "V56_NEW_UNEXPECTED_VISIT_EXPOSURE_MISSING_OR_ZERO")
        direct,a,w,exposure=numbers
        diff=abs(direct-(a+w))
        if diff>1e-5:bee_disagree+=1
        max_bee_disagreement=max(max_bee_disagreement,diff)
        s=bee[(fieldyear,band)]
        s["reported"]+=direct
        s["reconstructed"]+=a+w
        s["exposure"]+=exposure
        s["n"]+=1
        visit_count+=1

    require(dict(invalid)=={"2017":3} and seed_count==2214
            and visit_count==11682
            and abort_missing==1056 and triads==1161
            and bee_disagree==4087 and max_bee_disagreement==9,
            "V56_HISTORIC_SOURCE_INVALID_ROWS_OR_VISITS_COUNT_DRIFT")
    data=[]
    original_bands=Counter()
    for key in sorted(shared):
        h,v=meta_h[key],meta_v[key]
        require(all(len(h[k])==len(v[k])==1 and h[k]==v[k]
                    for k in ("region","type","physical")),
                "V56_ORIGINAL_HARVEST_VISIT_FIELD_IDENTITY_CHANGED")
        common=h["band"]&v["band"]
        original_bands[len(common)]+=1
        observations=[]
        for band in BANDS:
            if band not in common:continue
            s=seed[(key,band)]
            b=bee[(key,band)]
            require(s["n"]>0 and s["denominator"]>0
                    and b["n"]>0 and b["exposure"]>0,
                    "V56_ORIGINAL_DISTANCE_BAND_BECAME_UNOBSERVABLE")
            y=s["numerator"]/s["denominator"]
            reported=b["reported"]/b["exposure"]
            reconstructed=b["reconstructed"]/b["exposure"]
            require(0<=y<=1 and min(reported,reconstructed)>=0,
                    "V56_INVALID_AGGREGATED_SOURCE_RATIO")
            observations.append((y,reported,reconstructed))
        require(len(observations)==len(common)
                and 3<=len(common)<=4,
                "V56_3_4_BAND_REQUIREMENT_CHANGED")
        n=len(observations)
        region=next(iter(h["region"]))
        source_type=next(iter(h["type"]))
        physical=next(iter(h["physical"]))
        physical_key=(region,physical)
        data.append({
            "year":key[0],
            "region":region,"field_type":source_type,
            "physical_field":physical_key,
            "physical_field_hash":hashlib.sha256(
               ("ANTO_v55_Mota_source_physical_field:"+
                "\x1f".join(physical_key)).encode()).hexdigest(),
            "n_original_matched_distance_bands":n,
            "seed_fraction":sum(z[0] for z in observations)/n,
            "log_reported_total_rate":math.log1p(
                sum(z[1] for z in observations)/n),
            "log_reconstituted_rate":math.log1p(
                sum(z[2] for z in observations)/n),
        })
    audit={
        "original_archive_record_id":6671199,
        "original_MD5_source_files":reg["source_original_md5"],
        "original_source_harvest_rows":len(harvest),
        "original_source_visit_survey_rows":len(visits),
        "original_common_fieldyears":len(shared),
        "matched_harvest_rows_pre_numeric_gate":2217,
        "original_visitation_surveys_on_common_fieldyears":visit_count,
        "original_harvest_seed_count_ratio_invalid_rows_excluded":sum(invalid.values()),
        "excluded_harvest_seed_ratio_invalid_by_year":dict(sorted(invalid.items())),
        "valid_harvest_original_source_rows_kept":seed_count,
        "missing_Abort_entries_never_imputed":abort_missing,
        "complete_Abort_numerator_denominator_triads_audited":triads,
        "source_reported_total_vs_group_sum_disagreement_rows":bee_disagree,
        "maximum_source_reported_vs_group_sum_difference":max_bee_disagreement,
        "original_matched_distance_bands_per_fieldyear":
            {str(k):v for k,v in sorted(original_bands.items())},
        "n_eligible_fieldyear_unit_observations":len(data),
        "n_unique_original_physical_fields":len(
            {r["physical_field"] for r in data}),
        "fieldyears_by_year":dict(sorted(Counter(
            x["year"] for x in data).items())),
        "fieldyears_by_original_region":dict(sorted(Counter(
            x["region"] for x in data).items())),
        "fieldyears_by_original_field_type":dict(sorted(Counter(
            x["field_type"] for x in data).items())),
        "no_link_between_same_individual_sunflower_visits_and_harvest":True,
        "original_v54_and_v55_source_numerical_STOPS_retained":True,
    }
    require(audit["n_eligible_fieldyear_unit_observations"]==80
            and audit["n_unique_original_physical_fields"]==51
            and audit["original_matched_distance_bands_per_fieldyear"]==
                {"3":1,"4":79}
            and audit["fieldyears_by_year"]=={"2017":41,"2018":39}
            and audit["fieldyears_by_original_region"]=={
                "Burgos":38,"Cuenca":42}
            and audit["fieldyears_by_original_field_type"]=={
                "NonV":10,"SNV":32,"WFS":38},
            "V56_ORIGINAL_FIELDYEAR_OR_PHYSICAL_FIELD_TARGET_CHANGED")
    return data,audit


def report(data,audit,reg):
    original_fit_reg={
        "analysis_type":reg["analysis_type"],
        "original_file_md5":reg["source_original_md5"],
        "original_v53_fieldyear_blob_sha1":
            reg["v53_fieldyear_source_blob_sha1"],
        "original_v54_quality_blob_sha1":
            reg["v54_original_quality_blob_sha1"],
        "placebo_seed":reg["source_placebo_seed"],
        "placebo_repeats":reg["source_placebo_replicates"],
    }
    out=v55.report(data,audit,original_fit_reg)
    require(out["summary"]["n_original_physical_field_holdouts"]==51
            and out["summary"]["n_original_fieldyear_predictions"]==80
            and out["visitor_source_label_placebo"][
                "n_original_fieldyear_stratum_permutations"]==200,
            "V56_NEW_SOURCE_MODEL_GROUP_OR_PLACEBO_CHANGED")
    out["v55_source_ratio_quality_parent_blob_sha1"]=reg[
        "v55_invalid_fraction_quality_blob_sha1"]
    out["unchanged_original_v55_model_blob_sha1"]=reg[
        "v55_original_model_implementation_git_blob_sha1"]
    out["not_retroactive_v54_or_v55_success"]=True
    out["excluded_source_rows_reason_only"]={
        "No_Seeds_Total_gt_No_rep_units":3}
    out["original_independent_biological_studies_added"]=0
    out["new_source_plant_level_pollination_mechanism_estimated"]=False
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    r=registered()
    original_denominator_quality(r)
    payload=validate_source_archives(r)
    data,audit=aggregate_original(payload,r)
    result=report(data,audit,r)
    if args.check:
        require(result==json.loads(RECEIPT.read_text()),
                "V56_ORIGINAL_SOURCE_FROZEN_RESULT_NOT_REPRODUCED")
        print("V56_SOURCE_2214_VALID_SEED_ROWS_80_FIELDYEARS_51_PHYSICAL_FIELDS_EXACT_PASS")
    else:
        s=json.dumps(result,sort_keys=True,separators=(",",":"),
                     ensure_ascii=True,allow_nan=False)
        n=(len(s)+3999)//4000
        for i in range(n):
            print("V56_ORIGINAL_REPORT_CHUNK_%03d_OF_%03d="%
                  (i+1,n)+s[4000*i:4000*(i+1)])


if __name__=="__main__":
    main()
