#!/usr/bin/env python3
"""ANTO v55: source-locked Mota sunflower field-year bee/seed forecast.

Distinct POST-v54-source-STOP hypothesis, frozen before outcome modeling.
Uses original directly reported total interactions; guild reconstructed
counts are secondary measurement sensitivity, NOT additive by assumption.
Not the same original individual plants in pollinator and seed tables.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path
import analyze_original_Mota_sunflower_field_holdout_v54 as old

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v55_mota_reported_total_source_registry_20261119.json"
QUALITY=ROOT/"results/anto_v54_original_Mota_source_integrity_20261119.json"
RECEIPT=ROOT/"results/anto_v55_sunflower_total_bee_field_holdout_20261119.json"
MODELS=("Z0","Z1","Z2")
YEARS=("2017","2018")
REGIONS=("Burgos","Cuenca")
TYPES=("NonV","SNV","WFS")
BANDS=("0","15","30","60")
LAMBDA=8.0

class V55Error(ValueError):
    pass

def require(ok,msg):
    if not ok: raise V55Error(msg)

def git_blob(data):
    return hashlib.sha1(("blob %d\0"%len(data)).encode()+data).hexdigest()

def registry():
    r=json.loads(REGISTRY.read_text())
    require(
        r["analysis_type"]=="POST_V54_STOP_MOTA_REPORTED_TOTAL_BEE_FIELD_TRANSFER"
        and r["zenodo_record"]==6671199
        and r["original_file_md5"]=={
            "Database_Productivity.csv":"aeea30a21e495498295241bcfd8c693d",
            "Database_Visitation_rates.csv":"e13131ad60b133cc6ab25e20bbd3554b",
            "README_file.txt":"479b26d3cebb72c8491b0ec00609c08c"}
        and r["original_v53_fieldyear_blob_sha1"]==
            "7a4fb27c7f26abfa8400d23b041e1a3d3e3424c9"
        and r["original_v54_quality_blob_sha1"]==
            "55eb9b942e7a0ab7f4e32139ac5305325fb74a83"
        and r["n_source_harvest_rows_matched"]==2217
        and r["n_source_visit_rows_matched"]==11682
        and r["n_source_fieldyears"]==80
        and r["n_original_physical_fields"]==51
        and r["n_source_missing_aborts"]==1056
        and r["n_source_visits_total_not_apis_plus_wild"]==4087
        and r["max_absolute_visit_count_discrepancy"]==9
        and r["original_distance_bands"]==list(BANDS)
        and r["source_original_fieldyear_band_counts"]=={"3":1,"4":79}
        and r["original_years"]==list(YEARS)
        and r["original_regions"]==list(REGIONS)
        and r["original_field_types"]==list(TYPES)
        and r["source_original_fieldyear_by_year"]=={"2017":41,"2018":39}
        and r["source_original_fieldyear_by_region"]=={"Burgos":38,"Cuenca":42}
        and r["source_original_fieldyear_by_type"]==
            {"NonV":10,"SNV":32,"WFS":38}
        and r["models"]==list(MODELS)
        and r["ridge_alpha"]==8
        and r["n_expected_original_field_holdouts"]==51
        and r["n_expected_original_fieldyear_predictions"]==80
        and r["placebo_seed"]==20261119
        and r["placebo_repeats"]==200
        and r["source_model_fitted_before_protocol"] is False
        and r["v54_source_stop_retained"] is True
        and r["new_original_biological_studies_vs_v54"]==0
        and r["source_study_is_posthoc"] is True
        and r["same_individual_plant_visit_and_harvest"] is False
        and r["visitor_group_component_counts_assumed_equal_total"] is False
        and r["causal_randomization_p_value"] is False,
        "V55_REGISTERED_SOURCE_OR_PRE_OUTCOME_MODEL_DRIFT")
    return r

def require_prior_quality(reg):
    raw=QUALITY.read_bytes()
    require(git_blob(raw)==reg["original_v54_quality_blob_sha1"],
            "V55_ORIGINAL_V54_SOURCE_STOP_RECEIPT_REWRITTEN")
    z=json.loads(raw)
    require(
        z["source"]==6671199
        and z["original_file_md5"]==reg["original_file_md5"]
        and z["original_matched_fieldyears"]==80
        and z["original_matched_harvest_rows"]==2217
        and z["original_matched_visit_rows"]==11682
        and z["source_count_categories"]["visits_rows_total_vs_apis_wild_disagree"]==4087
        and z["per_source_numeric_column_statuses"]
            ["harvest_No_Aborts_Total"]["MISSING"]==1056
        and z["maximum_absolute_author_visitation_identity_difference"]==9
        and z["predictive_effect_or_visitor_seed_association_computed"] is False
        and z["original_v54_preregistered_source_identity_guard_modified"] is False,
        "V55_ORIGINAL_V54_NONADMISSION_REINTERPRETED")
    return z

def original_source(reg):
    prior=old.read_registry()
    require(prior["original_files_md5"]==reg["original_file_md5"]
            and prior["prior_original_v53_fieldyear_receipt_git_blob_sha1"]==
                reg["original_v53_fieldyear_blob_sha1"],
            "V55_ORIGINAL_ARCHIVE_IDENTITY_CHANGED")
    old.verified_anchor(prior)
    return old.source_bytes(prior)

def parse_original(raw,reg,prior):
    harvest=old.source_csv(raw["Database_Productivity.csv"],(
        "Year","Field_code","Region","Field_Type","Field_code_original",
        "Distance","No_Seeds_Total","No_rep_units","No_Aborts_Total"))
    visit=old.source_csv(raw["Database_Visitation_rates.csv"],(
        "Year","Field_code","Region","Field_Type","Field_ID","Distance",
        "No_Interactions_TOTAL","No_Interactions_APIS","No_Interactions_WILD",
        "Total_no_monitored_units"))
    require(len(harvest)==2875 and len(visit)==11683,
            "V55_ORIGINAL_AUTHOR_TOTAL_SOURCE_ROWS_CHANGED")
    hkeys=Counter(old.fieldyear_key(z) for z in harvest)
    vkeys=Counter(old.fieldyear_key(z) for z in visit)
    both=set(hkeys)&set(vkeys)
    require(len(both)==80
            and sum(hkeys[k] for k in both)==2217
            and sum(vkeys[k] for k in both)==11682,
            "V55_FIELDYEAR_JOIN_NOT_EXACT_ORIGINAL")
    meta={kind:defaultdict(lambda:{
        "region":set(),"type":set(),"physical":set(),"bands":set()})
        for kind in ("harvest","visit")}
    harvest_band=defaultdict(lambda:{"seeds":0.0,"units":0.0,"n":0})
    bee_band=defaultdict(lambda:{"reported_total":0.0,
                                 "reconstructed_total":0.0,"flower_minutes":0.0,"n":0})
    n_aborts_missing=n_visits_disagree=0
    max_visits_diff=0.0
    n_aborts_triads=0
    for which,records in (("harvest",harvest),("visit",visit)):
        for z in records:
            k=old.fieldyear_key(z)
            if k not in both:continue
            region,kind,physical=old.sanitize_metadata(
                z,"harvest" if which=="harvest" else "visits")
            dist=old.distance_key(z)
            m=meta[which][k]
            m["region"].add(region)
            m["type"].add(kind)
            m["physical"].add(physical)
            m["bands"].add(dist)
            key=(k,dist)
            if which=="harvest":
                full=old.parse_number(z["No_Seeds_Total"])
                units=old.parse_number(z["No_rep_units"])
                aborts=old.parse_number(z["No_Aborts_Total"])
                require(full is not None and units is not None and
                        0<=full<=units and units>0,
                        "V55_ORIGINAL_SEED_FRACTION_SOURCE_INVALID_OR_MISSING")
                if aborts is None:
                    n_aborts_missing+=1
                else:
                    n_aborts_triads+=1
                    require(abs(full+aborts-units)<=1e-5,
                            "V55_ORIGINAL_COMPLETE_SEED_TRIAD_CHANGED")
                h=harvest_band[key]
                h["seeds"]+=full
                h["units"]+=units
                h["n"]+=1
            else:
                nums=[old.parse_number(z[c]) for c in
                      ("No_Interactions_TOTAL","No_Interactions_APIS",
                       "No_Interactions_WILD","Total_no_monitored_units")]
                total,apis,wild,exposure=nums
                require(all(x is not None for x in nums) and exposure>0,
                        "V55_REPORTED_TOTAL_OR_BEE_COMPONENT_SOURCE_MISSING")
                gap=abs(total-apis-wild)
                max_visits_diff=max(max_visits_diff,gap)
                if gap>1e-5:n_visits_disagree+=1
                b=bee_band[key]
                b["reported_total"]+=total
                b["reconstructed_total"]+=apis+wild
                b["flower_minutes"]+=exposure
                b["n"]+=1
    require(n_aborts_missing==1056 and n_aborts_triads==1161
            and n_visits_disagree==4087 and max_visits_diff==9,
            "V55_SOURCE_QUALITY_V54_COUNTS_CHANGED")
    allrows=[]
    band_counts=Counter()
    for fieldyear in sorted(both):
        h=meta["harvest"][fieldyear]
        b=meta["visit"][fieldyear]
        require(all(len(h[z])==len(b[z])==1
                    and h[z]==b[z]
                    for z in ("region","type","physical")),
                "V55_V53_PHYSICAL_FIELDYEAR_METADATA_DISAGREE")
        common=h["bands"]&b["bands"]
        band_counts[len(common)]+=1
        observed=[]
        for band in BANDS:
            if band not in common:continue
            H=harvest_band[(fieldyear,band)]
            B=bee_band[(fieldyear,band)]
            require(H["units"]>0 and B["flower_minutes"]>0
                    and H["n"]>0 and B["n"]>0,
                    "V55_INCOMPLETE_ORIGINAL_SOURCE_DISTANCE_BAND")
            frac=H["seeds"]/H["units"]
            v1=B["reported_total"]/B["flower_minutes"]
            v2=B["reconstructed_total"]/B["flower_minutes"]
            require(0<=frac<=1 and v1>=0 and v2>=0,
                    "V55_BAND_SEED_FRACTION_OR_VISITS_OUT_OF_RANGE")
            observed.append((frac,v1,v2))
        require(len(observed)==len(common) and len(observed)>=3,
                "V55_ORIGINAL_MATCHED_SOURCE_BANDS_NOT_ALL_RETAINED")
        n=len(observed)
        region=next(iter(h["region"]))
        ftype=next(iter(h["type"]))
        physical=next(iter(h["physical"]))
        fieldid=(region,physical)
        allrows.append({
            "year":fieldyear[0],"region":region,"field_type":ftype,
            "physical_field":fieldid,
            "physical_field_hash":hashlib.sha256(
                ("ANTO_v55_Mota_source_physical_field:"+"\x1f".join(fieldid))
                .encode()).hexdigest(),
            "n_shared_original_distance_bands":n,
            "seed_fraction":sum(z[0] for z in observed)/n,
            "log_reported_total_rate":math.log1p(
                sum(z[1] for z in observed)/n),
            "log_reconstituted_rate":math.log1p(
                sum(z[2] for z in observed)/n),
        })
    audit={"n_original_2022_productivity_rows":len(harvest),
           "n_original_one_minute_bee_census_rows":len(visit),
           "n_exact_source_fieldyears_before_outcome":80,
           "n_source_seed_units_aggregated":2217,
           "n_source_visitation_rows_aggregated":11682,
           "original_abort_source_rows_missing_quarantined_not_imputed":
               n_aborts_missing,
           "original_apis_plus_wild_vs_reported_total_mismatch_count":
               n_visits_disagree,
           "original_visit_count_max_discrepancy":max_visits_diff,
           "original_common_distance_band_counts":
               dict(sorted(band_counts.items())),
           "n_eligible_original_fieldyears":len(allrows),
           "n_eligible_original_physical_fields":
               len({x["physical_field"] for x in allrows}),
           "eligible_by_year":dict(sorted(
               Counter(x["year"] for x in allrows).items())),
           "eligible_by_region":dict(sorted(
               Counter(x["region"] for x in allrows).items())),
           "eligible_by_original_field_type":dict(sorted(
               Counter(x["field_type"] for x in allrows).items())),
           "original_reproductive_outcome_is_fieldyear_not_individual":True,
           "original_visit_rate_is_census_not_same_harvest_plant":True,
           "previous_v54_source_nonadmission_preserved":True}
    require(len(allrows)==80 and audit["n_eligible_original_physical_fields"]==51
            and audit["eligible_by_year"]==reg["source_original_fieldyear_by_year"]
            and audit["eligible_by_region"]==reg["source_original_fieldyear_by_region"]
            and audit["eligible_by_original_field_type"]==
                reg["source_original_fieldyear_by_type"]
            and audit["original_common_distance_band_counts"]==
                reg["source_original_fieldyear_band_counts"],
            "V55_ORIGINAL_FIELDYEAR_OR_BAND_COHORT_DRIFT")
    return allrows,audit

def features(row,model,norm=None):
    require(model in MODELS,"V55_UNKNOWN_PREDECLARED_MODEL")
    x=[1.0,float(row["year"]=="2018"),
       float(row["region"]=="Cuenca"),
       float(row["field_type"]=="SNV"),
       float(row["field_type"]=="WFS")]
    if model!="Z0":
        k=("log_reported_total_rate" if model=="Z1" else
           "log_reconstituted_rate")
        mean,sd=norm[k]
        x.append((row[k]-mean)/sd)
    return x

def fit(train,model):
    require(model in MODELS and len(train)>=40,"V55_MODEL_OR_TRAIN_N_INVALID")
    norm={}
    if model!="Z0":
        k=("log_reported_total_rate" if model=="Z1" else
           "log_reconstituted_rate")
        vals=[r[k] for r in train]
        mean=sum(vals)/len(vals)
        sd=math.sqrt(sum((v-mean)**2 for v in vals)/len(vals))
        require(math.isfinite(sd) and sd>1e-12,
                "V55_TRAIN_EXPOSURE_ZERO_VARIANCE")
        norm[k]=(mean,sd)
    p=len(features(train[0],model,norm))
    mat=[[0.]*p for _ in range(p)]
    rhs=[0.]*p
    for row in train:
        x=features(row,model,norm)
        y=row["seed_fraction"]
        for i in range(p):
            rhs[i]+=x[i]*y
            for j in range(p):mat[i][j]+=x[i]*x[j]
    for i in range(1,p):mat[i][i]+=LAMBDA
    beta=old.solve(mat,rhs)
    return {"model":model,"norm":norm,"beta":beta,"ntrain":len(train)}

def predict(fitted,row):
    val=sum(a*b for a,b in zip(fitted["beta"],
                              features(row,fitted["model"],fitted["norm"])))
    require(math.isfinite(val),"V55_NONFINITE_LINEAR_PREDICTION")
    return max(0.0,min(1.0,val)),val

def holdouts(rows,models=MODELS):
    grouped=defaultdict(list)
    for x in rows:grouped[x["physical_field"]].append(x)
    require(len(grouped)==51 and sum(map(len,grouped.values()))==80,
            "V55_NOT_51_PHYSICAL_FIELDS_80_FIELDYEARS")
    out=[]
    for group in sorted(grouped):
        target=grouped[group]
        train=[r for r in rows if r["physical_field"]!=group]
        require(len(train)+len(target)==80
                and not ({r["physical_field"] for r in train}&
                         {r["physical_field"] for r in target}),
                "V55_ORIGINAL_FIELD_YEAR_LEAKAGE")
        fitted={name:fit(train,name) for name in models}
        losses={}
        clipping={}
        peryear=[]
        for m in models:
            pred=[predict(fitted[m],x) for x in target]
            loss=[(p[0]-x["seed_fraction"])**2
                  for p,x in zip(pred,target)]
            losses[m]=sum(loss)/len(loss)
            clipping[m]={
                "n_before_clip_below_zero":sum(p[1]<0 for p in pred),
                "n_before_clip_above_one":sum(p[1]>1 for p in pred)}
        for row in target:
            peryear.append({
                "year":row["year"],"field_type":row["field_type"],
                "squared_loss_by_model":{
                    m:(predict(fitted[m],row)[0]-row["seed_fraction"])**2
                    for m in models}})
        out.append({
            "physical_field_hash":target[0]["physical_field_hash"],
            "region":target[0]["region"],
            "n_original_fieldyears":len(target),
            "years":[x["year"] for x in target],
            "mse":losses,
            "raw_prediction_clipping":clipping,
            "heldout_original_fieldyear_loss":peryear})
    return out

def summary(folds):
    require(len(folds)==51 and
            sum(f["n_original_fieldyears"] for f in folds)==80,
            "V55_INCOMPLETE_PHYSICAL_FIELD_HOLDOUTS")
    models={}
    for m in MODELS:
        models[m]={
            "equal_physical_field_mse":sum(f["mse"][m] for f in folds)/51,
            "fieldyear_weighted_mse":sum(
                f["mse"][m]*f["n_original_fieldyears"] for f in folds)/80,
            "raw_predictions_below_zero":sum(
                f["raw_prediction_clipping"][m]["n_before_clip_below_zero"]
                for f in folds),
            "raw_predictions_above_one":sum(
                f["raw_prediction_clipping"][m]["n_before_clip_above_one"]
                for f in folds)}
    increments={}
    for m in ("Z1","Z2"):
        field_gains=[f["mse"]["Z0"]-f["mse"][m] for f in folds]
        fieldyear_gains=[
            (f["region"],t["field_type"],
             t["squared_loss_by_model"]["Z0"]-
             t["squared_loss_by_model"][m])
            for f in folds for t in f["heldout_original_fieldyear_loss"]]
        regmean={k:sum(g for rr,tt,g in fieldyear_gains if rr==k)/
                 sum(rr==k for rr,tt,g in fieldyear_gains)
                 for k in REGIONS}
        typemean={k:sum(g for rr,tt,g in fieldyear_gains if tt==k)/
                  sum(tt==k for rr,tt,g in fieldyear_gains)
                  for k in TYPES}
        increments["Z0_to_"+m]={
            "equal_physical_field_mse_gain":sum(field_gains)/51,
            "fieldyear_weighted_mse_gain":sum(
                f["n_original_fieldyears"]*g for f,g in zip(folds,field_gains))/80,
            "positive_physical_fields":sum(g>0 for g in field_gains),
            "negative_physical_fields":sum(g<0 for g in field_gains),
            "by_original_region_fieldyear_mean_gain":regmean,
            "by_original_field_type_fieldyear_mean_gain":typemean}
    return {"n_original_physical_field_holdouts":51,
            "n_original_fieldyear_predictions":80,
            "model_risks":models,"increments":increments}

def quantile(xs,p):
    ordered=sorted(xs)
    idx=(len(ordered)-1)*p
    j=int(idx)
    return ordered[j]*(1-(idx-j))+ordered[min(j+1,len(ordered)-1)]*(idx-j)

def placebo(rows,actual,results,reg):
    observed=results["increments"]["Z0_to_Z1"]["equal_physical_field_mse_gain"]
    known={f["physical_field_hash"]:f["mse"]["Z0"] for f in actual}
    strata=defaultdict(list)
    for i,r in enumerate(rows):
        strata[(r["region"],r["year"],r["field_type"])].append(i)
    rng=random.Random(reg["placebo_seed"])
    vals=[]
    for _ in range(reg["placebo_repeats"]):
        moved=[dict(r) for r in rows]
        for k in sorted(strata):
            indices=strata[k]
            source=[(rows[i]["log_reported_total_rate"],
                     rows[i]["log_reconstituted_rate"]) for i in indices]
            rng.shuffle(source)
            for i,pair in zip(indices,source):
                moved[i]["log_reported_total_rate"]=pair[0]
                moved[i]["log_reconstituted_rate"]=pair[1]
        shuffled=holdouts(moved,models=("Z1",))
        gains=[known[f["physical_field_hash"]]-f["mse"]["Z1"]
               for f in shuffled]
        require(len(gains)==51,
                "V55_PLACEBO_LOST_ORIGINAL_PHYSICAL_FIELD")
        vals.append(sum(gains)/len(gains))
    return {
        "n_original_fieldyear_stratum_permutations":len(vals),
        "seed":reg["placebo_seed"],
        "stratum":"Region×Year×Field_Type",
        "original_reported_and_reconstructed_rates_shuffled_jointly":True,
        "original_field_seed_outcomes_and_holdout_membership_fixed":True,
        "observed_equal_physical_field_gain":observed,
        "shuffle_median":quantile(vals,.5),
        "shuffle_q025":quantile(vals,.025),
        "shuffle_q975":quantile(vals,.975),
        "shuffle_max":max(vals),"shuffle_min":min(vals),
        "n_shuffle_gains_ge_observed":sum(x>=observed for x in vals),
        "is_causal_randomization_p_value":False}

def report(rows,audit,reg):
    folds=holdouts(rows)
    scores=summary(folds)
    null=placebo(rows,folds,scores,reg)
    positive=scores["increments"]["Z0_to_Z1"]
    passed=(positive["equal_physical_field_mse_gain"]>0
        and positive["fieldyear_weighted_mse_gain"]>0
        and positive["positive_physical_fields"]>51/2
        and all(x>0 for x in
                positive["by_original_region_fieldyear_mean_gain"].values())
        and sum(x>0 for x in positive[
            "by_original_field_type_fieldyear_mean_gain"].values())>=2
        and positive["equal_physical_field_mse_gain"]>
            null["shuffle_q975"])
    return {
        "analysis_type":reg["analysis_type"],
        "source_original_archive_record":6671199,
        "original_three_author_file_md5":reg["original_file_md5"],
        "v53_prior_fieldyear_source_blob_sha1":
            reg["original_v53_fieldyear_blob_sha1"],
        "v54_registered_quality_nonadmission_blob_sha1":
            reg["original_v54_quality_blob_sha1"],
        "source_audit":audit,
        "source_original_sunflower_total_bees_not_assumed_equal_guild_sum":True,
        "source_original_Abort_missing_never_zero_imputed":True,
        "original_physical_field_scores":[{
            "physical_field_hash":x["physical_field_hash"],
            "region":x["region"],
            "n_fieldyears":x["n_original_fieldyears"],
            "years":x["years"],
            "mse":x["mse"],
            "clipped_prediction_count":x["raw_prediction_clipping"],
            "original_fieldyear_derived_losses":x["heldout_original_fieldyear_loss"]}
            for x in folds],
        "summary":scores,
        "visitor_source_label_placebo":null,
        "registered_primary_field_transfer_verdict":(
            "DESCRIPTIVE_SINGLE_ORIGINAL_CROP_STUDY_FIELD_TRANSFER"
            if passed else "VISITOR_INFORMATION_FIELD_TRANSFER_NOT_ESTABLISHED"),
        "n_new_original_studies_vs_v54":0,
        "does_not_link_same_sunflower_plant_visitor_and_seed":True,
        "no_causal_pollination_or_genetic_evolution_estimated":True,
        "no_valid_randomization_p_value_computed":True,
        "source_original_individual_observations_not_exported":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    r=registry()
    require_prior_quality(r)
    raw=original_source(r)
    rows,audit=parse_original(raw,r,None)
    out=report(rows,audit,r)
    if args.check:
        require(out==json.loads(RECEIPT.read_text()),
                "V55_ORIGINAL_SOURCE_FULL_REGENERATION_MISMATCH")
        print("V55_AUTHOR_MD5_80_FIELDYEARS_51_PHYSICAL_FIELDS_200_PLACEBO_RECEIPT_PASS")
    else:
        content=json.dumps(out,sort_keys=True,separators=(",",":"),
                           ensure_ascii=True,allow_nan=False)
        n=(len(content)+3999)//4000
        for i in range(n):
            print("V55_ORIGINAL_SOURCE_READOUT_%03d_OF_%03d="%
                  (i+1,n)+content[4000*i:4000*(i+1)])

if __name__=="__main__":
    main()
