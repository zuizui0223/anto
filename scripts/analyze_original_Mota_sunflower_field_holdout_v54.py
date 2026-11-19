#!/usr/bin/env python3
"""Original sunflower one-study physical-field holdout (Mota et al. 2022).

Source is MD5-verified original Zenodo 6671199 and source fields/years
are pinned to ANTO v53. No original crop rows rehosted. This is a
POSTHOC predictive information test, not observed same-plant bee→seed
effects or an ecological causal pollination intervention.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import audit_independent_pollinator_reproduction_sources_v53 as v53

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v54_sunflower_field_source_registry_20261119.json"
SOURCE_ANCHOR=ROOT/"results/anto_v53_mota_fieldyear_source_key_audit_20261119.json"
RECEIPT=ROOT/"results/anto_v54_original_sunflower_field_heldout_20261119.json"
YEARS=("2017","2018")
REGIONS=("Burgos","Cuenca")
TYPES=("NonV","SNV","WFS")
BANDS=("0","15","30","60")
MODELS=("B0","B1","B2")
RADIUS=8.0


class V54Error(ValueError): pass

def need(ok,msg):
    if not ok: raise V54Error(msg)

def git_sha1(raw):
    return hashlib.sha1(("blob %d\0"%len(raw)).encode()+raw).hexdigest()

def read_registry():
    r=json.loads(REGISTRY.read_text())
    need(
        r["analysis_type"]=="ORIGINAL_MOTA_SUNFLOWER_LEAVE_PHYSICAL_FIELD_OUT_VISIT_VS_SEED_SET"
        and r["zenodo_record_id"]==6671199
        and r["original_files_md5"]=={
            "Database_Productivity.csv":"aeea30a21e495498295241bcfd8c693d",
            "Database_Visitation_rates.csv":"e13131ad60b133cc6ab25e20bbd3554b",
            "README_file.txt":"479b26d3cebb72c8491b0ec00609c08c"}
        and r["prior_original_v53_fieldyear_receipt_git_blob_sha1"]==
            "7a4fb27c7f26abfa8400d23b041e1a3d3e3424c9"
        and r["original_source_productivity_rows"]==2875
        and r["original_source_visitation_census_rows"]==11683
        and r["original_matched_fieldyears"]==80
        and r["original_matched_productivity_rows"]==2217
        and r["original_matched_visitation_rows"]==11682
        and r["original_matched_physical_fields"]==51
        and r["original_matched_fieldyears_by_year"]==
            {"2017":41,"2018":39}
        and r["original_matched_fieldyears_by_region"]==
            {"Burgos":38,"Cuenca":42}
        and r["original_matched_fieldyears_by_field_type"]==
            {"NonV":10,"SNV":32,"WFS":38}
        and r["original_common_distance_band_counts"]==
            {"3":1,"4":79}
        and r["original_distance_bands"]==list(BANDS)
        and r["original_field_type_categories"]==list(TYPES)
        and r["original_regions"]==list(REGIONS)
        and r["years"]==list(YEARS)
        and r["source_seed_set_numerator"]=="No_Seeds_Total"
        and r["source_seed_set_denominator"]=="No_rep_units"
        and r["visitor_per_one_minute_source_total"]=="No_Interactions_TOTAL"
        and r["visitor_source_apis"]=="No_Interactions_APIS"
        and r["visitor_source_wild"]=="No_Interactions_WILD"
        and r["source_visitation_denominator"]=="Total_no_monitored_units"
        and r["min_common_distance_bands"]==3
        and r["min_fieldyears_valid"]==70
        and r["min_physical_fields_valid"]==45
        and r["source_missing_records_max_fraction"]==.05
        and r["models"]==list(MODELS)
        and r["ridge_lambda"]==8
        and r["permutations"]==200
        and r["permutation_seed"]==20261119
        and r["original_individual_sunflower_visitation_linked_to_harvest"] is False
        and r["new_original_independent_studies_vs_v53_source_gate"]==0
        and r["study_is_posthoc"] is True
        and r["valid_randomization_p_value"] is False
        and r["numerical_original_seed_or_visitation_outcomes_inspected_before_freeze"] is False,
        "V54_ORIGINAL_SOURCE_OR_REGISTERED_PREDICTION_CONTRACT_CHANGED")
    return r


def verified_anchor(reg):
    raw=SOURCE_ANCHOR.read_bytes()
    need(git_sha1(raw)==reg["prior_original_v53_fieldyear_receipt_git_blob_sha1"],
         "V54_ORIGINAL_V53_FIELD_IDENTITY_ANCHOR_CHANGED")
    a=json.loads(raw)
    q=a["physical_field_metadata_join_audit"]
    need(
        a["n_shared_field_year_units"]==80
        and a["n_productivity_records_on_shared_field_year"]==2217
        and a["n_visitation_records_on_shared_field_year"]==11682
        and q["distinct_original_region_plus_field_names_in_shared_fieldyears"]==51
        and q["shared_fieldyear_same_original_physical_field_name"]==80
        and q["shared_fieldyear_region_fieldtype_consistent_and_singleton"]==80
        and q["shared_fieldyear_shared_distance_band_count"]=={"3":1,"4":79}
        and q["shared_fieldyear_by_original_field_type"]==
            {"NonV":10,"SNV":32,"WFS":38}
        and q["seed_or_visitation_response_data_examined"] is False,
        "V54_V53_SOURCE_ORIGINAL_IDENTITY_GRAIN_DRIFT")
    return a


def source_bytes(reg):
    old=v53.registry()
    cand=next(c for c in old["candidates"] if c["code"]=="Mota2022")
    need(cand["zenodo_record_id"]==reg["zenodo_record_id"]
         and cand["files_md5"]==reg["original_files_md5"],
         "V54_V53_SOURCE_ARCHIVE_CHANGED")
    return v53.fetch_original(cand,old["max_bytes_per_file"])


def source_csv(raw,required):
    try:
        text=raw.decode("utf-8-sig")
        dialect=csv.Sniffer().sniff(text[:16000],delimiters=",;")
        reader=csv.DictReader(io.StringIO(text),delimiter=dialect.delimiter)
        headers=tuple(reader.fieldnames or ())
        need(set(required).issubset(headers) and len(set(headers))==len(headers),
             "V54_SOURCE_CSV_HEADERS_CHANGED")
        out=list(reader)
        need(all(None not in row and all(z is not None for z in row.values())
                 for row in out),"V54_SOURCE_INVALID_ROW_SHAPE")
        return out
    except (UnicodeDecodeError,csv.Error) as exc:
        raise V54Error("V54_SOURCE_TEXT_CSV_PARSE") from exc


def parse_number(token,integer=True,positive=False):
    value=(token or "").strip()
    if value in ("","NA","NaN","nan",".","#VALUE!"):
        return None
    try:
        z=float(value)
    except (TypeError,ValueError) as exc:
        raise V54Error("V54_UNKNOWN_SOURCE_NUMERIC_TOKEN") from exc
    need(math.isfinite(z) and z>=0 and
         (not integer or abs(z-round(z))<=1e-6) and
         (not positive or z>0),
         "V54_SOURCE_NEGATIVE_NONFINITE_OR_INVALID_COUNT")
    return z


def fieldyear_key(row):
    y=row["Year"].strip()
    code=row["Field_code"].strip()
    need(y in YEARS and code,"V54_UNREGISTERED_YEAR_OR_FIELD_CODE")
    return (y,code)


def distance_key(row):
    token=row["Distance"].strip()
    try:
        value=float(token)
    except ValueError as exc:
        raise V54Error("V54_SOURCE_UNKNOWN_DISTANCE") from exc
    need(math.isfinite(value) and value.is_integer()
         and str(int(value)) in BANDS,
         "V54_SOURCE_DISTANCE_NOT_AN_ORIGINAL_FROZEN_BAND")
    return str(int(value))


def sanitize_metadata(row,which):
    region=row["Region"].strip()
    typ=row["Field_Type"].strip()
    physical=row["Field_code_original" if which=="harvest" else "Field_ID"].strip()
    need(region in REGIONS and typ in TYPES and physical,
         "V54_ORIGINAL_REGION_FIELD_TYPE_OR_PHYSICAL_NAME_CHANGED")
    return region,typ,physical


def aggregate_original_source(raw,reg,anchor):
    harvest=source_csv(raw["Database_Productivity.csv"],(
        "Region","Year","Field_code_original","Field_code","Field_Type",
        "Distance","No_Seeds_Total","No_Aborts_Total","No_rep_units"))
    bees=source_csv(raw["Database_Visitation_rates.csv"],(
        "Region","Year","Field_ID","Field_code","Field_Type","Distance",
        "No_Interactions_TOTAL","No_Interactions_APIS",
        "No_Interactions_WILD","Total_no_monitored_units"))
    need(len(harvest)==2875 and len(bees)==11683,
         "V54_ORIGINAL_MOTA_SOURCE_ROW_COUNT_DRIFT")
    hkeys=Counter(fieldyear_key(row) for row in harvest)
    bkeys=Counter(fieldyear_key(row) for row in bees)
    both=set(hkeys)&set(bkeys)
    need(len(both)==80 and sum(hkeys[k] for k in both)==2217
         and sum(bkeys[k] for k in both)==11682,
         "V54_ORIGINAL_COMMON_FIELDYEAR_KEYS_CHANGED")

    meta_h=defaultdict(lambda:{"region":set(),"type":set(),"physical":set(),"bands":set()})
    meta_b=defaultdict(lambda:{"region":set(),"type":set(),"physical":set(),"bands":set()})
    harvested=defaultdict(lambda:{"seeds":0.0,"units":0.0,"valid_rows":0})
    visits=defaultdict(lambda:{"total":0.0,"apis":0.0,"wild":0.0,
                               "exposure":0.0,"valid_rows":0})
    skipped_h=skipped_b=0
    for rows,which,meta in ((harvest,"harvest",meta_h),(bees,"visits",meta_b)):
        for row in rows:
            key=fieldyear_key(row)
            if key not in both:
                continue
            region,typ,field=sanitize_metadata(row,which)
            dist=distance_key(row)
            info=meta[key]
            info["region"].add(region)
            info["type"].add(typ)
            info["physical"].add(field)
            info["bands"].add(dist)
            k=(key,dist)
            if which=="harvest":
                nums=[parse_number(row[c]) for c in
                      ("No_Seeds_Total","No_Aborts_Total","No_rep_units")]
                if any(z is None for z in nums):
                    skipped_h+=1
                    continue
                full,empty,total=nums
                need(abs((full+empty)-total)<=1e-5,
                     "V54_HARVEST_ORIGINAL_SEED_ACCOUNTING_FAILED")
                if total<=0:
                    skipped_h+=1
                    continue
                item=harvested[k]
                item["seeds"]+=full
                item["units"]+=total
                item["valid_rows"]+=1
            else:
                nums=[parse_number(row[c]) for c in
                     ("No_Interactions_TOTAL","No_Interactions_APIS",
                      "No_Interactions_WILD","Total_no_monitored_units")]
                if any(z is None for z in nums):
                    skipped_b+=1
                    continue
                total,apis,wild,exposure=nums
                need(abs((apis+wild)-total)<=1e-5,
                     "V54_BEE_ORIGINAL_APIS_PLUS_WILD_COUNT_IDENTITY_FAILED")
                if exposure<=0:
                    skipped_b+=1
                    continue
                item=visits[k]
                item["total"]+=total;item["apis"]+=apis;item["wild"]+=wild
                item["exposure"]+=exposure
                item["valid_rows"]+=1
    need(skipped_h<=.05*2217 and skipped_b<=.05*11682,
         "V54_TOO_MANY_ORIGINAL_NUMERIC_OR_ZERO_DENOMINATOR_ROWS")
    result=[]
    bandcount=Counter()
    for key in sorted(both):
        h,b=meta_h[key],meta_b[key]
        need(all(len(x[k])==1 for x in (h,b) for k in
                 ("region","type","physical"))
             and all(h[k]==b[k] for k in
                     ("region","type","physical")),
             "V54_MISMATCHED_OR_MULTIPLE_ORIGINAL_FIELD_YEAR_SOURCE_KEYS")
        common=h["bands"]&b["bands"]
        bandcount[len(common)]+=1
        valid=[]
        for dist in BANDS:
            if dist not in common:continue
            H=harvested[(key,dist)]
            B=visits[(key,dist)]
            if H["units"]<=0 or B["exposure"]<=0:
                continue
            seed=H["seeds"]/H["units"]
            total=B["total"]/B["exposure"]
            apis=B["apis"]/B["exposure"]
            wild=B["wild"]/B["exposure"]
            need(0<=seed<=1 and min(total,apis,wild)>=0
                 and abs(total-apis-wild)<=1e-7,
                 "V54_BAND_AGGREGATION_OUT_OF_SUPPORT")
            valid.append({"distance":dist,"seed":seed,
                          "visitors":total,"apis":apis,"wild":wild})
        if len(valid)<reg["min_common_distance_bands"]:
            continue
        region=next(iter(h["region"]))
        typ=next(iter(h["type"]))
        physical=next(iter(h["physical"]))
        group=(region,physical)
        no=len(valid)
        y=sum(x["seed"] for x in valid)/no
        v=sum(x["visitors"] for x in valid)/no
        a=sum(x["apis"] for x in valid)/no
        w=sum(x["wild"] for x in valid)/no
        result.append({
            "year":key[0],"field_code":key[1],"region":region,
            "field_type":typ,"physical_field":group,
            "physical_field_sha256":hashlib.sha256(
                ("ANTO_v54_physical_field_source_20261119:"+
                 "\x1f".join(group)).encode("utf8")).hexdigest(),
            "n_common_valid_distance_bands":no,
            "seed_fraction":y,
            "log_total_visit_rate":math.log1p(v),
            "log_apis_visit_rate":math.log1p(a),
            "log_wild_visit_rate":math.log1p(w)
        })
    physical={x["physical_field"] for x in result}
    need(len(result)>=reg["min_fieldyears_valid"]
         and len(physical)>=reg["min_physical_fields_valid"]
         and Counter(x["year"] for x in result).keys()==set(YEARS)
         and Counter(x["region"] for x in result).keys()==set(REGIONS)
         and Counter(x["field_type"] for x in result).keys()==set(TYPES),
         "V54_INSUFFICIENT_ORIGINAL_FIELDYEAR_OR_PHYSICAL_FIELD_COVERAGE")
    # An eligible predictor may NOT be a single harvest or visit observation.
    return result,{
        "source_exact_files_md5":reg["original_files_md5"],
        "original_harvest_rows_total":len(harvest),
        "original_visitation_survey_rows_total":len(bees),
        "original_common_fieldyears_before_numeric_gates":len(both),
        "original_harvest_rows_on_common_fieldyears":2217,
        "original_visitation_rows_on_common_fieldyears":11682,
        "numeric_missing_or_zero_denominator_harvest_rows":skipped_h,
        "numeric_missing_or_zero_exposure_visitor_surveys":skipped_b,
        "original_common_band_count_before_numeric_gates":dict(sorted(
            bandcount.items())),
        "eligible_source_fieldyears":len(result),
        "eligible_unique_physical_fields":len(physical),
        "eligible_fieldyears_by_year":dict(sorted(Counter(
            x["year"] for x in result).items())),
        "eligible_fieldyears_by_original_region":dict(sorted(Counter(
            x["region"] for x in result).items())),
        "eligible_fieldyears_by_original_field_type":dict(sorted(Counter(
            x["field_type"] for x in result).items())),
        "source_individual_sunflower_visitation_to_harvest_key":False,
        "source_one_minute_census_flower_exposure_used_once":True,
        "source_equal_distance_band_aggregation":True,
    }


def design(row,model,stats=None):
    x=[1.,float(row["year"]=="2018"),
       float(row["region"]=="Cuenca"),
       float(row["field_type"]=="SNV"),
       float(row["field_type"]=="WFS")]
    if model in ("B1","B2"):
        keys=("log_total_visit_rate",) if model=="B1" else (
             "log_apis_visit_rate","log_wild_visit_rate")
        for k in keys:
            mean,sd=stats[k]
            x.append((row[k]-mean)/sd)
    return x


def fit(rows,model):
    need(model in MODELS and len(rows)>=40,
         "V54_INVALID_TRAINING_COHORT_OR_MODEL")
    stats={}
    if model!="B0":
        features=("log_total_visit_rate",) if model=="B1" else (
                  "log_apis_visit_rate","log_wild_visit_rate")
        for k in features:
            xs=[row[k] for row in rows]
            mu=sum(xs)/len(xs)
            sd=math.sqrt(sum((a-mu)**2 for a in xs)/len(xs))
            need(math.isfinite(sd) and sd>1e-12,
                 "V54_VISIT_PREDICTOR_UNIDENTIFIABLE_IN_TRAIN")
            stats[k]=(mu,sd)
    X=[design(r,model,stats) for r in rows]
    Y=[r["seed_fraction"] for r in rows]
    p=len(X[0])
    a=[[0.]*p for _ in range(p)]
    b=[0.]*p
    for x,y in zip(X,Y):
        for i in range(p):
            b[i]+=x[i]*y
            for j in range(p):a[i][j]+=x[i]*x[j]
    for i in range(1,p):a[i][i]+=RADIUS
    coef=solve(a,b)
    return {"model":model,"stats":stats,"coef":coef,
            "n_original_training_fieldyears":len(rows)}


def solve(A,b):
    n=len(b)
    need(n>0 and len(A)==n and all(len(x)==n for x in A),
         "V54_INVALID_ORIGINAL_RIDGE_MATRIX")
    aug=[list(A[i])+[b[i]] for i in range(n)]
    for j in range(n):
        k=max(range(j,n),key=lambda x:abs(aug[x][j]))
        need(abs(aug[k][j])>1e-11,
             "V54_SINGULAR_ORIGINAL_FIELD_RIDGE")
        aug[j],aug[k]=aug[k],aug[j]
        den=aug[j][j]
        aug[j]=[v/den for v in aug[j]]
        for i in range(n):
            if i==j:continue
            mult=aug[i][j]
            aug[i]=[v-mult*w for v,w in zip(aug[i],aug[j])]
    out=[x[-1] for x in aug]
    need(all(math.isfinite(x) for x in out),
         "V54_NONFINITE_SOURCE_COEFFICIENT")
    return out


def predict(model,row):
    raw=sum(c*x for c,x in zip(model["coef"],
                             design(row,model["model"],model["stats"])))
    need(math.isfinite(raw),"V54_NONFINITE_RAW_PREDICTION")
    return (max(0.,min(1.,raw)),raw)


def physical_groups(rows):
    grouped=defaultdict(list)
    for row in rows:grouped[row["physical_field"]].append(row)
    return grouped


def leave_field_out(rows,models=MODELS):
    grouped=physical_groups(rows)
    need(len(grouped)>=45 and sum(map(len,grouped.values()))==len(rows),
         "V54_INVALID_PHYSICAL_FIELD_SCORES")
    folders=[]
    for key in sorted(grouped):
        train=[r for r in rows if r["physical_field"]!=key]
        test=grouped[key]
        need(train and test and all(r["physical_field"]!=key for r in train)
             and set(r["physical_field_sha256"] for r in test)==
                 {test[0]["physical_field_sha256"]},
             "V54_HELDOUT_ORIGINAL_PHYSICAL_FIELD_LEAKAGE")
        learned={name:fit(train,name) for name in models}
        losses={}
        clip={}
        for name in models:
            preds=[predict(learned[name],r) for r in test]
            losses[name]=sum((p[0]-r["seed_fraction"])**2
                             for p,r in zip(preds,test))/len(test)
            clip[name]={"raw_below_zero":sum(p[1]<0 for p in preds),
                        "raw_above_one":sum(p[1]>1 for p in preds)}
        folders.append({
            "physical_field_sha256":test[0]["physical_field_sha256"],
            "region":test[0]["region"],
            "n_fieldyears":len(test),
            "heldout_fieldyear_types":[x["field_type"] for x in test],
            "heldout_fieldyear_years":[x["year"] for x in test],
            "n_training_other_physical_fields":len(grouped)-1,
            "mse":losses,"clip_counts":clip,
            "original_fieldyear_scores":[{
                "year":r["year"],"field_type":r["field_type"],
                "mse_by_model":{name:(predict(learned[name],r)[0]-
                         r["seed_fraction"])**2 for name in models}}
                for r in test]})
    return folders


def score_summary(folds):
    nfields=len(folds)
    nfieldyears=sum(f["n_fieldyears"] for f in folds)
    need(nfields>=45 and nfieldyears>=70,
         "V54_INCOMPLETE_FIELD_HELDOUT_SCORE")
    stats={}
    for m in MODELS:
        stats[m]={
          "equal_physical_field_mse":sum(f["mse"][m] for f in folds)/nfields,
          "fieldyear_weighted_mse":sum(f["mse"][m]*f["n_fieldyears"]
                       for f in folds)/nfieldyears,
          "total_raw_below_zero":sum(f["clip_counts"][m]["raw_below_zero"]
                                    for f in folds),
          "total_raw_above_one":sum(f["clip_counts"][m]["raw_above_one"]
                                   for f in folds)}
    gains={}
    for before,after,label in (("B0","B1","B0_to_B1_total_bees"),
                               ("B0","B2","B0_to_B2_apis_and_wild")):
        delta=[f["mse"][before]-f["mse"][after] for f in folds]
        by_fieldyear=[
            (f["region"],y["field_type"],y["mse_by_model"][before]-
             y["mse_by_model"][after])
            for f in folds for y in f["original_fieldyear_scores"]]
        gmean=sum(delta)/len(delta)
        gwt=sum(f["n_fieldyears"]*d for f,d in zip(folds,delta))/nfieldyears
        byreg={name:(sum(x for region,typ,x in by_fieldyear if region==name)/
                        sum(1 for region,typ,x in by_fieldyear if region==name))
               for name in REGIONS}
        bytype={name:(sum(x for region,typ,x in by_fieldyear if typ==name)/
                         sum(1 for region,typ,x in by_fieldyear if typ==name))
                for name in TYPES}
        gains[label]={
            "equal_physical_field_absolute_mse_reduction":gmean,
            "fieldyear_weighted_absolute_mse_reduction":gwt,
            "positive_physical_fields":sum(x>0 for x in delta),
            "negative_physical_fields":sum(x<0 for x in delta),
            "per_region_mean_fieldyear_mse_gain":byreg,
            "per_original_field_type_mean_fieldyear_mse_gain":bytype}
    return {"original_heldout_physical_fields":nfields,
            "original_heldout_fieldyears":nfieldyears,
            "models":stats,"gains":gains}


def quantile(a,p):
    values=sorted(a)
    idx=(len(values)-1)*p
    lo=int(idx)
    return values[lo]*(1-(idx-lo))+values[min(lo+1,len(values)-1)]*(idx-lo)


def placebo(rows,original_folds,original,reg):
    observed=original["gains"]["B0_to_B1_total_bees"][
        "equal_physical_field_absolute_mse_reduction"]
    baseline={f["physical_field_sha256"]:f["mse"]["B0"]
              for f in original_folds}
    strata=defaultdict(list)
    for i,r in enumerate(rows):
        strata[(r["region"],r["year"],r["field_type"])].append(i)
    need(set(r["region"] for r in rows)==set(REGIONS)
         and set(r["field_type"] for r in rows)==set(TYPES),
         "V54_PLACEBO_INVALID_ORIGINAL_GRAIN")
    randomizer=random.Random(reg["permutation_seed"])
    vals=[]
    for _ in range(reg["permutations"]):
        permuted=[dict(r) for r in rows]
        for key in sorted(strata):
            indices=strata[key]
            # Keep each original ecology-stratum and the original
            # joint bee-count signal together during relabelling.
            bees=[(rows[i]["log_total_visit_rate"],
                   rows[i]["log_apis_visit_rate"],
                   rows[i]["log_wild_visit_rate"]) for i in indices]
            randomizer.shuffle(bees)
            for i,b in zip(indices,bees):
                permuted[i]["log_total_visit_rate"]=b[0]
                permuted[i]["log_apis_visit_rate"]=b[1]
                permuted[i]["log_wild_visit_rate"]=b[2]
        simulated=leave_field_out(permuted,models=("B1",))
        benefits=[baseline[f["physical_field_sha256"]]-f["mse"]["B1"]
                  for f in simulated]
        need(len(benefits)==len(baseline),
             "V54_PSEUDO_PLACEBO_NOT_ALL_ORIGINAL_FIELDS")
        vals.append(sum(benefits)/len(benefits))
    return {
        "replicate_count":len(vals),"seed":reg["permutation_seed"],
        "original_shuffle_stratum":"Region×Year×Field_Type at source fieldyear level",
        "source_original_seed_fraction_and_field_holdout_unchanged":True,
        "observed_equal_physical_field_gain":observed,
        "shuffle_median":quantile(vals,.5),
        "shuffle_q025":quantile(vals,.025),
        "shuffle_q975":quantile(vals,.975),
        "shuffle_min":min(vals),"shuffle_max":max(vals),
        "n_shuffles_at_or_above_observed":sum(x>=observed for x in vals),
        "valid_causal_randomization_p_value":False}


def report(rows,audit,reg):
    folds=leave_field_out(rows)
    summ=score_summary(folds)
    sensitivity=placebo(rows,folds,summ,reg)
    g=summ["gains"]["B0_to_B1_total_bees"]
    passing=(
        g["equal_physical_field_absolute_mse_reduction"]>0 and
        g["fieldyear_weighted_absolute_mse_reduction"]>0 and
        g["positive_physical_fields"]>summ["original_heldout_physical_fields"]/2 and
        all(v>0 for v in g["per_region_mean_fieldyear_mse_gain"].values()) and
        sum(v>0 for v in
            g["per_original_field_type_mean_fieldyear_mse_gain"].values())>=2 and
        g["equal_physical_field_absolute_mse_reduction"]>
            sensitivity["shuffle_q975"])
    out={
        "analysis_type":reg["analysis_type"],
        "zenodo_original_record":6671199,
        "original_three_file_md5":reg["original_files_md5"],
        "original_v53_source_fieldyear_receipt_git_blob_sha1":
            reg["prior_original_v53_fieldyear_receipt_git_blob_sha1"],
        "source_audit":audit,
        "full_physical_field_holdout_scores":[{
            "physical_field_sha256":f["physical_field_sha256"],
            "region":f["region"],"n_fieldyears":f["n_fieldyears"],
            "heldout_years":f["heldout_fieldyear_years"],
            "heldout_field_types":f["heldout_fieldyear_types"],
            "MSE_by_original_model":f["mse"],
            "original_fieldyear_derived_squared_losses":
                f["original_fieldyear_scores"],
            "model_raw_clip_counts":f["clip_counts"]}
            for f in folds],
        "summary":summ,
        "within_original_study_visit_identity_placebo":sensitivity,
        "registered_descriptive_field_transfer_verdict":(
            "DESCRIPTIVE_SINGLE_ORIGINAL_CROP_STUDY_FIELD_TRANSFER"
            if passing else "VISITOR_INFORMATION_FIELD_TRANSFER_NOT_ESTABLISHED"),
        "n_additional_independent_original_studies_vs_v53":0,
        "fitted_individual_sunflower_visitor_to_harvest_effect":False,
        "causal_visitor_intervention_or_genetic_inheritance_estimated":False,
        "this_is_a_posthoc_source_locked_predictive_study":True,
        "independent_randomization_p_value_calculated":False,
        "original_source_raw_rows_redistributed":False}
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    reg=read_registry()
    verified_anchor(reg)
    raw=source_bytes(reg)
    rows,audit=aggregate_original_source(raw,reg,None)
    outcome=report(rows,audit,reg)
    if args.check:
        need(outcome==json.loads(RECEIPT.read_text()),
             "V54_ORIGINAL_SOURCE_MODEL_RECEIPT_NOT_EXACT")
        print("V54_ORIGINAL_MOTA_MD5_51_FIELD_LORO_REPRO_AND_PLACEBO_PASS")
    else:
        report_string=json.dumps(outcome,sort_keys=True,
                                 separators=(",",":"),ensure_ascii=True)
        n=(len(report_string)+3999)//4000
        for i in range(n):
            print("V54_ORIGINAL_CROP_READOUT_%03d_OF_%03d="%
                  (i+1,n)+report_string[4000*i:4000*(i+1)])


if __name__=="__main__":
    main()
