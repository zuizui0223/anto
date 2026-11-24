#!/usr/bin/env python3
"""ANTO v58: one independent original Mertensia plant experiment, array holdout.

Exactly the verified original Gallagher & Campbell 2021 workbook, source
Dryad SHA256 and Zenodo MD5, 115 original plant rows, 105 predeclared
complete plants, 20 (week,array) intact withheld arrays. No raw author
plant data are committed. Observed visitation is NOT randomized.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import re
import statistics
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import audit_original_Mertensia_Zenodo_equivalence_v57 as orig
import audit_original_Mertensia_exp3_sameplant_grain_v57 as grain

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v58_mertensia_original_source_model_registry_20261124.json"
PARENT=ROOT/"results/anto_v57_mertensia_exp3_original_plant_grain_20261124.json"
RECEIPT=ROOT/"results/anto_v58_mertensia_array_heldout_visitation_20261124.json"
MODELS=("Q0","Q1","Q2","Q3")
NUMERIC=("flrCnt","mcorW","mcorL","mVisRate","mseed")
ID_FIELDS=("week","array","treat","plantID")
WEEKS=("1","2","3","4")
TREAT=("A","C","D")
ALPHA=5.0
SOURCE_XLSX="gallagher&campbell_frontiersinecology&evolution_data.xlsx"
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
FEATURES={
    "Q0":("log_flower",),
    "Q1":("log_flower","mcorW","mcorL"),
    "Q2":("log_flower","mcorW","mcorL","log_visit"),
    "Q3":("log_flower","log_visit")
}


class V58Error(ValueError):
    pass


def need(ok,reason):
    if not ok:
        raise V58Error(reason)


def gitblob(raw):
    return hashlib.sha1(("blob %d\0"%len(raw)).encode()+raw).hexdigest()


def registry():
    r=json.loads(REGISTRY.read_text())
    need(
        r["analysis_type"]==
            "POSTHOC_AUTHENTIC_MERTENSIA_ARRAY_HELDOUT_FLOWER_VISITOR_SEED_INFORMATION"
        and r["source_dryad_doi"]=="10.7280/D16D7Z"
        and r["source_original_dryad_id"]==59594
        and r["source_dryad_version"]==4
        and r["source_zenodo_record"]==4731053
        and r["original_xlsx_name"]==SOURCE_XLSX
        and r["original_xlsx_size"]==339816
        and r["original_xlsx_sha256"]==
           "152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d"
        and r["original_xlsx_zenodo_md5"]==
           "db86f8286007b6176bf957334eba9ca4"
        and r["original_author_sheet_name"]=="exp3_meansForAnalysis"
        and r["original_xlsx_worksheet"]=="xl/worksheets/sheet11.xml"
        and r["parent_v57_source_grain_blob_sha1"]==
           "5659545f1376aeb715aefa4308a0fe1e8673da39"
        and r["source_original_plant_record_rows"]==115
        and r["source_distinct_week_array_treatment_plant_keys"]==115
        and r["source_complete_eligible_plant_records"]==105
        and r["source_original_week_array_groups"]==20
        and r["source_groups_min_raw_plant_rows"]==5
        and r["source_groups_max_raw_plant_rows"]==7
        and r["source_weeks"]==list(WEEKS)
        and r["source_water_treatments"]==list(TREAT)
        and r["source_candidate_fields"]==list(NUMERIC)
        and r["models"]==list(MODELS)
        and r["ridge_alpha"]==ALPHA
        and r["full_week_array_holdouts"]==20
        and r["full_original_eligible_plant_predictions"]==105
        and r["min_complete_plants_per_group"]==3
        and r["permutation_seed"]==20261124
        and r["permutation_replicates"]==200
        and r["no_outcome_derived_predictor"] is True
        and r["no_inherited_phenotype_gene_effect"] is True
        and r["original_source_not_rehosted"] is True
        and r["research_question_is_posthoc"] is True,
        "V58_SOURCE_FIT_OR_INFERENCE_CONTRACT_CHANGED"
    )
    need(r["Q0"] and r["Q1"] and r["Q2"] and r["Q3"],
         "V58_MODEL_REGISTER_FIELDS_ABSENT")
    return r


def parent_source(r):
    raw=PARENT.read_bytes()
    need(gitblob(raw)==r["parent_v57_source_grain_blob_sha1"],
         "V58_PREOUTCOME_V57_SOURCE_GRAIN_RECEIPT_CHANGED")
    p=json.loads(raw)
    g=p["original_source_group_grain"]
    need(p["source_original_author_workbook_authentic"] is True
         and p["original_workbook_sha256"]==r["original_xlsx_sha256"]
         and g["original_distinct_composite_week_array_treat_plantID_keys"]==115
         and g["n_original_plant_records_all_candidate_variables_PRESENT"]==105
         and g["n_original_week_array_clusters"]==20
         and g["source_row_duplicate_composite_plantkeys"]==0
         and g["source_only_prerequisite_three_original_group_cohort_ready"] is True
         and g["source_original_measurement_magnitudes_examined"] is False,
         "V58_PREOUTCOME_MERTENSIA_SOURCE_COHORT_ALTERED")
    return p


def get_original_excel(r):
    original=orig.get_original_official_metadata()
    archival=orig.verify_zenodo()
    raw=orig.fetch_zenodo_exact(
        archival[SOURCE_XLSX],SOURCE_XLSX,r["original_xlsx_size"])
    author_ok,original_type=orig.checksum_matches(
        raw,original[SOURCE_XLSX]["digest"])
    zenodo_ok,archival_type=orig.checksum_matches(
        raw,archival[SOURCE_XLSX]["checksum"])
    need(author_ok and zenodo_ok
         and original_type=="SHA256" and archival_type=="MD5"
         and hashlib.sha256(raw).hexdigest()==r["original_xlsx_sha256"],
         "V58_ORIGINAL_AUTHOR_DOUBLE_SOURCE_DIGEST_NOT_ACCEPTED")
    return raw


def interpret_cell(cell,shared):
    # Formulas are allowed only when they contain a source cached numeric
    # result; no formula expression is evaluated and no missing is imputed.
    typ=cell.attrib.get("t")
    if typ=="inlineStr":
        x=cell.find("m:is",NS)
        return ("".join(t.text or "" for t in x.findall(".//m:t",NS))
                if x is not None else "")
    v=cell.find("m:v",NS)
    if v is None:
        return ""
    if typ=="s":
        i=int(v.text or "-1")
        need(0<=i<len(shared),"V58_SOURCE_SHARED_STRING_INDEX_INVALID")
        return shared[i]
    if typ=="e":
        return "#ERROR:"+str(v.text or "")
    return v.text or ""


def load_original_plant_rows(raw,r,parent):
    with zipfile.ZipFile(io.BytesIO(raw)) as book:
        need(book.testzip() is None and
             r["original_xlsx_worksheet"] in book.namelist(),
             "V58_ORIGINAL_AUTHOR_ANALYSIS_WORKSHEET_MISSING")
        shared=[]
        if "xl/sharedStrings.xml" in book.namelist():
            root=ET.fromstring(book.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si",NS):
                shared.append("".join(t.text or "" for t in si.findall(".//m:t",NS)))
        sheet=ET.fromstring(book.read(r["original_xlsx_worksheet"]))
        lines=sheet.findall(".//m:sheetData/m:row",NS)
        need(len(lines)>=116,"V58_ORIGINAL_ANALYSIS_WORKSHEET_ROW_COUNT_DRIFT")
        def selected(row):
            d={}
            for cell in row.findall("m:c",NS):
                ref=cell.attrib.get("r","")
                m=re.match(r"^([A-Z]+)(\d+)$",ref)
                if not m:continue
                col=grain.colnumber(m.group(1))
                d[col]=interpret_cell(cell,shared)
            return d
        first=selected(lines[0])
        headers={i:v.strip() for i,v in first.items() if v.strip()}
        reverse={val:key for key,val in headers.items()}
        need(set(ID_FIELDS+NUMERIC).issubset(reverse),
             "V58_ORIGINAL_SOURCE_SHEET_FIELD_NAMES_CHANGED")
        if "mseed" not in headers.values() or "mVisRate" not in headers.values():
            raise V58Error("V58_SOURCE_REPRO_VISITS_NOT_SAME_AUTHOR_SHEET")
        originals=[]
        found=Counter()
        missing=Counter()
        groups=Counter()
        duplicates=0
        fullkeys=set()
        for line in lines[1:]:
            data=selected(line)
            original={k:str(data.get(reverse[k],"")).strip()
                      for k in ID_FIELDS}
            if not any(original.values()):
                continue
            need(all(original.values()),
                 "V58_ORIGINAL_EXPERIMENT_PLANT_KEY_MISSING")
            key=tuple(original[k] for k in ID_FIELDS)
            if key in fullkeys:duplicates+=1
            fullkeys.add(key)
            originals.append((original,data))
            groups[(original["week"],original["array"])]+=1
            for name in NUMERIC:
                cell=str(data.get(reverse[name],"")).strip()
                if cell.lower() in ("","na","nan","n/a",".","#n/a","#value!"):
                    missing[name]+=1
                else:
                    found[name]+=1
        need(len(originals)==115 and len(fullkeys)==115 and duplicates==0
             and len(groups)==20 and min(groups.values())>=5
             and max(groups.values())<=7,
             "V58_ORIGINAL_SAMEPLANT_ROW_GROUP_GRAIN_CHANGED")
        verified=parent["original_source_group_grain"]
        need(found["flrCnt"]==verified["candidate_var_present_count"]["flrCnt"]
             and found["mseed"]==verified["candidate_var_present_count"]["mseed"]
             and found["mcorW"]==verified["candidate_var_present_count"]["mcorW"]
             and found["mcorL"]==verified["candidate_var_present_count"]["mcorL"]
             and found["mVisRate"]==verified["candidate_var_present_count"]["mVisRate"],
             "V58_ORIGINAL_CANDIDATE_FIELD_PRESENCE_COUNTS_DRIFT")
        data=[]
        for labels,record in originals:
            text={k:str(record.get(reverse[k],"")).strip() for k in NUMERIC}
            if any(v.lower() in ("","na","nan","n/a",".","#n/a","#value!")
                   for v in text.values()):
                continue
            values={}
            for k,v in text.items():
                try:z=float(v)
                except ValueError as exc:
                    raise V58Error("V58_NONNUMERIC_ORIGINAL_CANDIDATE_CELL_"+k) from exc
                need(math.isfinite(z) and z>=0,
                     "V58_SOURCE_NEGATIVE_OR_NONFINITE_VALUE_"+k)
                if k in ("mcorW","mcorL"):
                    need(z>0,"V58_NONPOSITIVE_ORIGINAL_COROLLA_SIZE")
                values[k]=z
            wk,arr,treat=labels["week"],labels["array"],labels["treat"]
            need(wk in WEEKS and arr in ("1","2","3","4","5")
                 and treat in TREAT,
                 "V58_ORIGINAL_SOURCE_WEEK_ARRAY_TREATMENT_CHANGED")
            data.append({
                "week":wk,"array":arr,"treat":treat,
                "id":tuple(labels[k] for k in ID_FIELDS),
                "group":(wk,arr),
                "log_flower":math.log1p(values["flrCnt"]),
                "mcorW":values["mcorW"],"mcorL":values["mcorL"],
                "log_visit":math.log1p(values["mVisRate"]),
                "y":math.log1p(values["mseed"])
            })
    counts=Counter(z["group"] for z in data)
    need(len(data)==105 and len(counts)==20
         and min(counts.values())>=r["min_complete_plants_per_group"]
         and set(z["week"] for z in data)==set(WEEKS)
         and set(z["treat"] for z in data)==set(TREAT),
         "V58_PREDECLARED_ORIGINAL_COMPLETE_CASE_OR_HOLDOUT_GRAIN_CHANGED")
    report={"original_author_plant_sheet_rows":115,
            "exact_original_complete_numeric_plants":len(data),
            "exact_original_composite_plantkeys":len(fullkeys),
            "original_20_week_array_groups":len(counts),
            "original_complete_plants_per_group":{
                "|".join(k):v for k,v in sorted(counts.items())},
            "original_water_treatment_names_used_as_source_codes":list(TREAT),
            "original_phenology_weeks":list(WEEKS),
            "original_source_mseed_used_as_original_coded_proxy_no_unit_claim":True,
            "source_original_author_sheet":r["original_author_sheet_name"],
            "original_source_raw_rows_exported":False}
    return data,report


def vector(r,model,stats):
    need(model in MODELS,"V58_UNKNOWN_MODEL")
    x=[1.,
       float(r["treat"]=="C"),float(r["treat"]=="D"),
       float(r["week"]=="2"),float(r["week"]=="3"),
       float(r["week"]=="4")]
    for f in FEATURES[model]:
        mean,sd=stats[f]
        x.append((r[f]-mean)/sd)
    return x


def solve(A,y):
    n=len(y)
    need(n>0 and len(A)==n and all(len(x)==n for x in A),
         "V58_BAD_RIDGE_NORMAL_EQUATION_SHAPE")
    m=[list(A[i])+[y[i]] for i in range(n)]
    for j in range(n):
        p=max(range(j,n),key=lambda x:abs(m[x][j]))
        need(abs(m[p][j])>1e-12,"V58_SOURCE_RIDGE_SINGULAR")
        m[j],m[p]=m[p],m[j]
        div=m[j][j]
        m[j]=[x/div for x in m[j]]
        for i in range(n):
            if i==j:continue
            factor=m[i][j]
            m[i]=[a-factor*b for a,b in zip(m[i],m[j])]
    result=[row[-1] for row in m]
    need(all(math.isfinite(x) for x in result),"V58_SOURCE_RIDGE_NONFINITE")
    return result


def fit(train,model):
    need(model in MODELS and len(train)>=90,"V58_TRAINING_GROUP_OR_MODEL_INVALID")
    fields=FEATURES[model]
    stats={}
    for f in fields:
        vals=[z[f] for z in train]
        mean=sum(vals)/len(vals)
        sd=math.sqrt(sum((a-mean)**2 for a in vals)/len(vals))
        need(sd>1e-12 and math.isfinite(sd),
             "V58_TRAIN_NUMERIC_FEATURE_DEGENERATE_"+f)
        stats[f]=(mean,sd)
    X=[vector(z,model,stats) for z in train]
    n=len(X[0])
    A=[[0.]*n for _ in range(n)]
    b=[0.]*n
    for row,item in zip(X,train):
        for i in range(n):
            b[i]+=row[i]*item["y"]
            for j in range(n):
                A[i][j]+=row[i]*row[j]
    for i in range(1,n):
        A[i][i]+=ALPHA
    return {"model":model,"beta":solve(A,b),"stats":stats,
            "original_training_plants":len(train)}


def predict(fitted,r):
    vec=vector(r,fitted["model"],fitted["stats"])
    return sum(a*b for a,b in zip(vec,fitted["beta"]))


def run_holdouts(data,models=MODELS):
    group=defaultdict(list)
    for z in data:group[z["group"]].append(z)
    need(len(group)==20 and sum(len(x) for x in group.values())==105,
         "V58_NOT_TWENTY_SOURCE_ARRAY_GROUPS")
    outcome=[]
    for grp in sorted(group):
        train=[z for z in data if z["group"]!=grp]
        test=group[grp]
        need(set(z["id"] for z in train).isdisjoint(z["id"] for z in test)
             and len(test)>=3 and len(train)+len(test)==105,
             "V58_SOURCE_ARRAY_FOLD_OR_PLANT_LEAKAGE")
        learned={m:fit(train,m) for m in models}
        scores={}
        for model in models:
            err=sum((z["y"]-predict(learned[model],z))**2 for z in test)
            scores[model]=err/len(test)
        outcome.append({"week":grp[0],"array":grp[1],
                        "n_original_heldout_plants":len(test),
                        "n_original_training_plants":len(train),
                        "water_treatments_in_test":dict(sorted(
                            Counter(z["treat"] for z in test).items())),
                        "original_test_mse":scores})
    return outcome


def summarize(folds):
    need(len(folds)==20 and sum(x["n_original_heldout_plants"]
                                 for x in folds)==105,
         "V58_INCOMPLETE_SOURCE_20_ARRAY_RISK")
    models={}
    for m in MODELS:
        models[m]={
            "equal_array_mse":sum(f["original_test_mse"][m] for f in folds)/20,
            "plant_weighted_mse":sum(f["original_test_mse"][m]*
                       f["n_original_heldout_plants"] for f in folds)/105}
    diffs={}
    for name,before,after in [
        ("Q0_to_Q1_morphology","Q0","Q1"),
        ("Q1_to_Q2_conditional_visitor","Q1","Q2"),
        ("Q0_to_Q3_visitor_without_morphology","Q0","Q3")]:
        ds=[x["original_test_mse"][before]-x["original_test_mse"][after]
            for x in folds]
        perweek={}
        for week in WEEKS:
            y=[d for f,d in zip(folds,ds) if f["week"]==week]
            need(len(y)==5,"V58_FIVE_SOURCE_ARRAYS_PER_ORIGINAL_WEEK_CHANGED")
            perweek[week]={"n_original_arrays":len(y),
                           "equal_array_mse_gain":sum(y)/len(y),
                           "positive_original_arrays":sum(v>0 for v in y)}
        diffs[name]={
            "equal_array_absolute_mse_gain":sum(ds)/20,
            "plant_weighted_absolute_mse_gain":sum(
                f["n_original_heldout_plants"]*d for f,d in zip(folds,ds))/105,
            "positive_original_arrays":sum(v>0 for v in ds),
            "negative_original_arrays":sum(v<0 for v in ds),
            "by_original_phenology_week":perweek}
    return {"original_heldout_arrays":20,
            "original_plant_predictions":105,
            "model_risk":models,
            "increment":diffs}


def quantile(arr,p):
    x=sorted(arr)
    pos=(len(x)-1)*p
    j=int(pos)
    return x[j]*(1-(pos-j))+x[min(j+1,len(x)-1)]*(pos-j)


def visitor_placebo(data,folds,summary,registered):
    baseline={x["week"]+"|"+x["array"]:x["original_test_mse"]["Q1"]
              for x in folds}
    obs=summary["increment"]["Q1_to_Q2_conditional_visitor"][
        "equal_array_absolute_mse_gain"]
    strata=defaultdict(list)
    for i,z in enumerate(data):
        strata[(z["week"],z["treat"])].append(i)
    need(len(strata)==12 and min(len(x) for x in strata.values())>=6,
         "V58_ORIGINAL_WEEK_WATER_PLACEBO_SOURCE_GROUP_CHANGED")
    rng=random.Random(registered["permutation_seed"])
    vals=[]
    for _ in range(registered["permutation_replicates"]):
        randomdata=[dict(z) for z in data]
        for group in sorted(strata):
            ix=strata[group]
            visitors=[data[i]["log_visit"] for i in ix]
            rng.shuffle(visitors)
            for i,value in zip(ix,visitors):
                randomdata[i]["log_visit"]=value
        permuted=run_holdouts(randomdata,models=("Q2",))
        gains=[baseline[x["week"]+"|"+x["array"]]-
               x["original_test_mse"]["Q2"] for x in permuted]
        vals.append(sum(gains)/20)
    need(len(vals)==200 and all(math.isfinite(v) for v in vals),
         "V58_VISITOR_PLACEBO_INCOMPLETE")
    return {
        "n_within_original_week_water_visit_identity_placebos":200,
        "original_rng_seed":registered["permutation_seed"],
        "source_test_arrays_and_seed_outcomes_fixed":True,
        "original_visitor_rate_shuffled_within_week_treat":True,
        "observed_equal_array_visit_mse_gain":obs,
        "placebo_median":quantile(vals,.5),
        "placebo_q025":quantile(vals,.025),
        "placebo_q975":quantile(vals,.975),
        "placebo_max":max(vals),
        "placebo_min":min(vals),
        "n_source_shuffles_at_or_above_observed":
            sum(v>=obs for v in vals),
        "causal_visitation_randomization_p_value":False}


def report(data,audit,r):
    folds=run_holdouts(data)
    summary=summarize(folds)
    placebo=visitor_placebo(data,folds,summary,r)
    inc=summary["increment"]["Q1_to_Q2_conditional_visitor"]
    positive=(
        inc["equal_array_absolute_mse_gain"]>0
        and inc["plant_weighted_absolute_mse_gain"]>0
        and inc["positive_original_arrays"]>=12
        and sum(x["equal_array_mse_gain"]>0 for x in
                inc["by_original_phenology_week"].values())>=3
        and inc["equal_array_absolute_mse_gain"]>placebo["placebo_q975"]
    )
    return {
        "analysis_type":r["analysis_type"],
        "original_author_dryad_doi":r["source_dryad_doi"],
        "original_author_data_version":4,
        "original_zenodo_archival_record":4731053,
        "original_workbook_sha256":r["original_xlsx_sha256"],
        "prior_original_source_grain_git_blob_sha1":
            r["parent_v57_source_grain_blob_sha1"],
        "original_experiment":"Mertensia original Gallagher Campbell 2016 original 20 week-array groups within one 2021 paper",
        "original_coded_response":"log1p(source mseed) without inferred trait genetic fitness units",
        "original_coded_visit_feature":"log1p(source mVisRate) observed plant visit rate",
        "original_coded_flower_features":["log1p(flrCnt)","mcorW","mcorL"],
        "water_levels_preserved_as_original_source_codes":list(TREAT),
        "original_source_cohort_audit":audit,
        "model_features":{m:list(FEATURES[m]) for m in MODELS},
        "fixed_ridge_alpha":5,
        "all_20_original_week_array_holdout_risks":folds,
        "summary":summary,
        "original_200_visitor_identity_placebo":placebo,
        "preoutcome_descriptive_transfer_gate_verdict":(
            "DESCRIPTIVE_ORIGINAL_MERTENSIA_ARRAY_HELDOUT_VISITOR_INFORMATION"
            if positive else
            "THIRD_ORIGINAL_MERTENSIA_VISITOR_TRANSFER_NOT_ESTABLISHED"),
        "new_original_plant_studies_collected_in_ANTO":0,
        "source_original_paper_independent_of_2022_sinapis_and_mota":True,
        "plant_visit_rate_randomly_manipulated":False,
        "new_causal_pollination_p_value_or_hereditary_fitness_effect":False,
        "original_author_plant_raw_records_rehosted":False
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    r=registry()
    parent=parent_source(r)
    original=get_original_excel(r)
    plants,audit=load_original_plant_rows(original,r,parent)
    result=report(plants,audit,r)
    if args.check:
        need(result==json.loads(RECEIPT.read_text()),
             "V58_ORIGINAL_MERTENSIA_AUTHOR_SOURCE_RECEIPT_REPLAY_FAILED")
        print("V58_ORIGINAL_DRYAD_ZENODO_105_PLANTS_20_ARRAYS_200_PLACEBOS_EXACT_PASS")
    else:
        j=json.dumps(result,sort_keys=True,separators=(",",":"),
                     allow_nan=False,ensure_ascii=True)
        n=(len(j)+3999)//4000
        for i in range(n):
            print("V58_AUTHENTIC_MERTENSIA_SOURCE_%03d_OF_%03d="%
                  (i+1,n)+j[4000*i:4000*(i+1)])

if __name__=="__main__":main()
