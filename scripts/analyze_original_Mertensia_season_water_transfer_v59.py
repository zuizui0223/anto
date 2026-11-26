#!/usr/bin/env python3
"""Post-v58 Mertensia context-transfer diagnostic. Same 105 original plants.

No new ecological study, no causal visit effect. Whole-week and week4-safe
holdouts address author-described late-week array rearrangement. Verify all
v58 losses and original source identity before any new prediction.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
import random
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import analyze_original_Mertensia_20array_visitation_v58 as v58
import audit_original_Mertensia_Zenodo_equivalence_v57 as source
import audit_original_Mertensia_exp3_sameplant_grain_v57 as grain

ROOT=Path(__file__).resolve().parents[1]
REG=ROOT/"data/external/anto_v59_mertensia_season_water_registry_20261126.json"
RECEIPT=ROOT/"results/anto_v59_mertensia_season_water_transfer_20261126.json"
FROZEN_V58=ROOT/"results/anto_v58_mertensia_array_heldout_visitation_20261124.json"
CONTEXTS=("no_week","linear_week","linear_week_vwc")
POLICIES=("array20","safe16","week4","water3")
BASE=("log_flower","mcorW","mcorL")
NS=v58.NS
WEEKS=v58.WEEKS
TREAT=v58.TREAT
DOCX="gallagher&campbell_frontiersinecology&evolution_metaDataFile.docx"
WNS={"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


class V59Error(ValueError):pass
def need(ok,reason):
    if not ok:raise V59Error(reason)

def registry():
    r=json.loads(REG.read_text())
    need(r["version"]=="v59" and r["n_plants"]==105
         and r["n_original_groups"]==20 and r["n_weeks"]==4
         and r["n_treatments"]==3 and r["ridge_alpha"]==5
         and r["base_features"]==list(BASE)
         and r["primary_holdouts"]==["array20","safe16","week4"]
         and r["secondary_holdout"]=="water3"
         and r["primary_models"]==["linear_week","linear_week_vwc"]
         and r["secondary_model"]=="no_week"
         and r["placebo_n"]==200 and r["placebo_seed"]==20261126
         and r["original_week4_rearranged"] is True
         and r["new_original_studies"]==0 and r["causal_claim"] is False,
         "V59_REGISTERED_CONTEXT_OR_GRAIN_CHANGED")
    need(v58.gitblob(FROZEN_V58.read_bytes())==r["v58_result_git_blob"]
         and v58.gitblob(Path(v58.__file__).read_bytes())==r["v58_script_git_blob"],
         "V59_ORIGINAL_V58_MODEL_OR_RECEIPT_CHANGED")
    need(r["source_sha256"]==v58.registry()["original_xlsx_sha256"],
         "V59_ORIGINAL_WORKBOOK_SOURCE_DRIFT")
    return r

def source_dictionary(r):
    original=source.get_original_official_metadata()
    archived=source.verify_zenodo()
    b=source.fetch_zenodo_exact(archived[DOCX],DOCX,source.FILES[DOCX])
    a,at=source.checksum_matches(b,original[DOCX]["digest"])
    z,zt=source.checksum_matches(b,archived[DOCX]["checksum"])
    need(a and z and at=="SHA256" and zt=="MD5"
         and hashlib.sha256(b).hexdigest()==r["metadata_sha256"],
         "V59_ORIGINAL_AUTHOR_DICTIONARY_DIGEST_MISMATCH")
    with zipfile.ZipFile(io.BytesIO(b)) as f:
        xml=ET.fromstring(f.read("word/document.xml"))
    # Source dictionary definition snippets only; never output raw plant rows.
    entries=[]
    for row in xml.findall(".//w:tr",WNS):
        txt=" ".join(t.text or "" for t in row.findall(".//w:t",WNS)).strip()
        if txt:entries.append(txt)
    for para in xml.findall(".//w:p",WNS):
        txt=" ".join(t.text or "" for t in para.findall(".//w:t",WNS)).strip()
        if txt:entries.append(txt)
    ans={}
    for field in ("mseed","mVisRate","mVWC","plantID","array"):
        matches=[]
        for txt in entries:
            if re.search(r"\b"+field+r"\b",txt,re.I):
                s=" ".join(txt.split())
                if s not in matches:matches.append(s[:240])
        ans[field]=matches[:3]
    return {"original_dictionary_sha256":r["metadata_sha256"],
            "definitions_bounded_to_registered_field_names":ans}

def original_extra_rows(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as f:
        shared=[]
        if "xl/sharedStrings.xml" in f.namelist():
            e=ET.fromstring(f.read("xl/sharedStrings.xml"))
            for node in e.findall("m:si",NS):
                shared.append("".join(t.text or "" for t in node.findall(".//m:t",NS)))
        sheet=ET.fromstring(f.read("xl/worksheets/sheet11.xml"))
        rows=sheet.findall(".//m:sheetData/m:row",NS)
        def read(row):
            out={}
            for cell in row.findall("m:c",NS):
                m=re.match(r"^([A-Z]+)\d+$",cell.attrib.get("r",""))
                if m:out[grain.colnumber(m.group(1))]=v58.interpret_cell(cell,shared).strip()
            return out
        cols={value:k for k,value in read(rows[0]).items() if value}
        need(all(k in cols for k in v58.ID_FIELDS+("mVWC",)),
             "V59_ORIGINAL_EXTRA_SOURCE_FIELDS_MISSING")
        extra={}
        for row in rows[1:]:
            d=read(row)
            key=tuple(d.get(cols[k],"").strip() for k in v58.ID_FIELDS)
            if not any(key):continue
            need(all(key) and key not in extra,"V59_ORIGINAL_EXTRA_SOURCE_ID_DUPLICATE")
            value=d.get(cols["mVWC"],"")
            extra[key]=value
        need(len(extra)==115,"V59_ORIGINAL_115_SOURCE_IDENTITIES_CHANGED")
        return extra

def identity_audit(data,extra):
    # Source author used distinct potted cohorts across weeks. Raw plantID
    # labels reused across weeks or treatment groups cannot be assumed distinct.
    # Count tokens only; never export plant IDs.
    def info(keys):
        p=defaultdict(list)
        for key in keys:p[key[3]].append(key)
        return {
            "n_distinct_raw_plantID_tokens":len(p),
            "n_tokens_repeated":sum(len(x)>1 for x in p.values()),
            "n_tokens_in_multiple_weeks":sum(len({z[0] for z in x})>1 for x in p.values()),
            "n_tokens_in_multiple_treatments":sum(len({z[2] for z in x})>1 for x in p.values()),
            "n_tokens_in_multiple_arrays_within_weeks1to3":sum(
                any(len({z[1] for z in x if z[0]==w})>1 for w in ("1","2","3"))
                for x in p.values()),
            "n_tokens_in_multiple_week4_arrays":sum(
                len({z[1] for z in x if z[0]=="4"})>1 for x in p.values())}
    a=info(extra)
    b=info([x["id"] for x in data])
    need(b["n_tokens_in_multiple_weeks"]==0,
         "V59_SAME_RAW_PLANTID_RECURS_ACROSS_WEEKS_SOURCE_SEMANTICS_REQUIRED")
    need(b["n_tokens_in_multiple_arrays_within_weeks1to3"]==0,
         "V59_SAME_RAW_PLANTID_CROSSES_EARLY_ARRAYS")
    need(b["n_tokens_in_multiple_treatments"]==0,
         "V59_SAME_RAW_PLANTID_CROSSES_WATER_TREATMENTS")
    return {"all_115_source_keys":a,"eligible_105_source_keys":b,
            "original_week4_environment_rearrangement_reported_in_article":True}

def hydrate(data,extra):
    out=[]
    for item in data:
        row=dict(item)
        need(row["id"] in extra,"V59_SOURCE_VWC_SAMEPLANT_JOIN_MISSING")
        try:v=float(extra[row["id"]])
        except (TypeError,ValueError) as exc:
            raise V59Error("V59_SOURCE_VWC_NOT_NUMERIC") from exc
        need(math.isfinite(v) and 0<=v<=100,"V59_VWC_OUTSIDE_REGISTERED_RANGE")
        row["mVWC"]=v
        row["week_num"]=float(row["week"])
        out.append(row)
    need(len(out)==105,"V59_COHORT_105_NOT_RETAINED")
    return out

def group_key(row,policy):
    if policy=="array20":return row["week"]+"|"+row["array"]
    if policy=="safe16":
        return "4|ALL" if row["week"]=="4" else row["week"]+"|"+row["array"]
    if policy=="week4":return "week|"+row["week"]
    if policy=="water3":return "water|"+row["treat"]
    raise V59Error("V59_UNREGISTERED_SPLIT")

def splits(data,policy):
    groups=defaultdict(list)
    for i,x in enumerate(data):groups[group_key(x,policy)].append(i)
    counts={"array20":20,"safe16":16,"week4":4,"water3":3}
    need(len(groups)==counts[policy] and len(data)==105,"V59_FOLD_GRAIN_DRIFT")
    out=[]
    for k,ix in sorted(groups.items()):
        testset=set(ix)
        tr=[i for i in range(len(data)) if i not in testset]
        need(len(tr)+len(ix)==105 and len(ix)>=4,"V59_INVALID_FOLD_SUPPORT")
        if policy!="array20":
            need({data[i]["id"][3] for i in tr}.isdisjoint({data[i]["id"][3] for i in ix}),
                 "V59_POTENTIAL_SAME_PHYSICAL_PLANT_LEAKAGE")
        if policy=="week4":
            need(set(data[i]["week"] for i in tr).isdisjoint(data[i]["week"] for i in ix),
                 "V59_TEST_WEEK_LEAKAGE")
        if policy=="water3":
            need(set(data[i]["treat"] for i in tr).isdisjoint(data[i]["treat"] for i in ix),
                 "V59_TEST_WATER_CODE_LEAKAGE")
        out.append((k,tr,ix))
    need(sorted(i for k,tr,ix in out for i in ix)==list(range(105)),
         "V59_INCOMPLETE_OR_DUPLICATED_TEST_COVERAGE")
    return out

def feature_names(context,visit):
    fs=list(BASE)
    if context!="no_week":fs.append("week_num")
    if context=="linear_week_vwc":fs.append("mVWC")
    if visit:fs.append("log_visit")
    return fs

def design(row,fs,stats,policy):
    x=[1.]
    if policy!="water3":x.extend([float(row["treat"]=="C"),float(row["treat"]=="D")])
    x.extend((row[f]-stats[f][0])/stats[f][1] for f in fs)
    return x

def fit(train,context,visit,policy):
    fs=feature_names(context,visit)
    stats={}
    for f in fs:
        a=[x[f] for x in train]
        mean=sum(a)/len(a)
        sd=math.sqrt(sum((z-mean)**2 for z in a)/len(a))
        need(math.isfinite(sd) and sd>1e-12,"V59_ZERO_TRAIN_VARIANCE_"+f)
        stats[f]=(mean,sd)
    X=[design(row,fs,stats,policy) for row in train]
    p=len(X[0]);A=[[0.]*p for _ in range(p)];b=[0.]*p
    for row,z in zip(X,train):
        for i in range(p):
            b[i]+=row[i]*z["y"]
            for j in range(p):A[i][j]+=row[i]*row[j]
    for i in range(1,p):A[i][i]+=5.
    beta=v58.solve(A,b)
    return {"fs":fs,"stats":stats,"beta":beta,"policy":policy}

def predict(fit,row):
    return sum(a*b for a,b in zip(fit["beta"],design(row,fit["fs"],fit["stats"],fit["policy"])))

def errors(row,b,v):
    e=row["y"]-b
    d=v-b
    return {"base_sq":e*e,"visit_sq":(e-d)**2,
            "alignment":2*e*d,"correction_sq":d*d,
            "delta_pred":d,"base_residual":e}

def agg(rows):
    n=len(rows)
    need(n>0,"V59_EMPTY_SCORE_GROUP")
    m={k:sum(row[k] for row in rows)/n
       for k in ("base_sq","visit_sq","alignment","correction_sq","delta_pred","base_residual")}
    gain=m["base_sq"]-m["visit_sq"]
    need(abs(gain-(m["alignment"]-m["correction_sq"]))<1e-10,
         "V59_CORRECTION_ALGEBRA_MISMATCH")
    return {"n":n,"B_mse":m["base_sq"],"V_mse":m["visit_sq"],
            "gain":gain,"alignment_2ed":m["alignment"],
            "correction_d2":m["correction_sq"],
            "mean_prediction_correction":m["delta_pred"],
            "mean_baseline_residual":m["base_residual"],
            "n_plants_improved":sum(x["base_sq"]>x["visit_sq"] for x in rows)}

def outside(train,test,fields):
    out={}
    for f in fields:
        lo=min(x[f] for x in train);hi=max(x[f] for x in train)
        out[f]=sum(x[f]<lo or x[f]>hi for x in test)
    return out

def evaluate(data,policy,context,return_private=False):
    grouped=defaultdict(list);folds=[];allrows=[]
    for label,tri,tei in splits(data,policy):
        tr=[data[i] for i in tri];te=[data[i] for i in tei]
        b=fit(tr,context,False,policy);v=fit(tr,context,True,policy)
        out=[]
        for x in te:
            score=errors(x,predict(b,x),predict(v,x))
            score.update({"week":x["week"],"treat":x["treat"]})
            out.append(score);allrows.append(score)
        e=agg(out)
        e.update({"fold":label,"n_train":len(tr),
                  "heldout_weeks":dict(sorted(Counter(x["week"] for x in te).items())),
                  "heldout_treatments":dict(sorted(Counter(x["treat"] for x in te).items())),
                  "train_visit_slope_per_log_rate_unit":v["beta"][-1]/v["stats"]["log_visit"][1],
                  "n_test_outside_training_range":outside(tr,te,BASE+("log_visit","mVWC","week_num"))})
        folds.append(e)
    byweek={w:agg([x for x in allrows if x["week"]==w]) for w in WEEKS}
    bytreat={t:agg([x for x in allrows if x["treat"]==t]) for t in TREAT}
    equal_fold={k:sum(f[k] for f in folds)/len(folds) for k in ("B_mse","V_mse","gain")}
    equal_week={k:sum(byweek[w][k] for w in WEEKS)/4 for k in ("B_mse","V_mse","gain")}
    result={"policy":policy,"context":context,"n_folds":len(folds),"n_plants":len(allrows),
            "fold_equal":equal_fold,"plant_weighted":agg(allrows),"week_equal":equal_week,
            "n_positive_folds":sum(f["gain"]>0 for f in folds),
            "n_positive_weeks":sum(f["gain"]>0 for f in byweek.values()),
            "by_week":byweek,"by_water_code":bytreat,"folds":folds}
    return (result,allrows) if return_private else result

def replay_v58(data):
    frozen=json.loads(FROZEN_V58.read_text())
    expected={(x["week"],x["array"]):x for x in frozen["all_20_original_week_array_holdout_risks"]}
    folds=v58.run_holdouts(data)
    ncheck=0;maxerr=0.
    for x in folds:
        orig=expected[(x["week"],x["array"])]
        need(x["n_original_heldout_plants"]==orig["n_original_heldout_plants"],
             "V59_V58_FROZEN_GROUP_N_MISMATCH")
        for m in v58.MODELS:
            delta=abs(x["original_test_mse"][m]-orig["original_test_mse"][m])
            maxerr=max(maxerr,delta);ncheck+=1
            need(delta<=1e-10,"V59_FROZEN_V58_RISK_REPLAY_FAILED")
    scores=[]
    for k,tri,tei in splits(data,"array20"):
        tr=[data[i] for i in tri];te=[data[i] for i in tei]
        b=v58.fit(tr,"Q1");v=v58.fit(tr,"Q2")
        for x in te:
            y=errors(x,v58.predict(b,x),v58.predict(v,x))
            y.update({"week":x["week"],"treat":x["treat"]})
            scores.append(y)
    return {"all_original_v58_80_fold_model_losses_reproduced":ncheck==80,
            "maximum_absolute_mse_error":maxerr,
            "original_plant_weighted_correction_decomposition_by_week":{
                w:agg([x for x in scores if x["week"]==w]) for w in WEEKS},
            "original_plant_weighted_correction_decomposition_all":agg(scores)}

def placebo(data,observed,r):
    # Same test groups and original baseline predictions; within source
    # week x water label shuffle, deliberately descriptive, not causal.
    folds=splits(data,"week4")
    strata=defaultdict(list)
    for i,x in enumerate(data):strata[(x["week"],x["treat"])].append(i)
    need(len(strata)==12,"V59_ORIGINAL_VISITOR_PLACEBO_STRATA_CHANGED")
    baseline={f["fold"]:f["B_mse"] for f in observed["folds"]}
    rng=random.Random(r["placebo_seed"])
    values=[]
    for _ in range(r["placebo_n"]):
        perm=[dict(x) for x in data]
        for key,ids in sorted(strata.items()):
            v=[data[i]["log_visit"] for i in ids]
            rng.shuffle(v)
            for i,val in zip(ids,v):perm[i]["log_visit"]=val
        gains=[]
        for key,tri,tei in folds:
            tr=[perm[i] for i in tri];te=[perm[i] for i in tei]
            v=fit(tr,"linear_week",True,"week4")
            mse=sum((x["y"]-predict(v,x))**2 for x in te)/len(te)
            gains.append(baseline[key]-mse)
        values.append(sum(gains)/len(gains))
    actual=observed["week_equal"]["gain"]
    q=v58.quantile(values,.975)
    return {"n":len(values),"observed_week_equal_gain":actual,
            "median":v58.quantile(values,.5),"q025":v58.quantile(values,.025),
            "q975":q,"min":min(values),"max":max(values),
            "n_at_or_above":sum(x>=actual for x in values),
            "above_q975":actual>q,
            "reference_is_not_a_valid_causal_p_value":True}

def run():
    r=registry()
    parent=v58.parent_source(v58.registry())
    raw=v58.get_original_excel(v58.registry())
    old,audit=v58.load_original_plant_rows(raw,v58.registry(),parent)
    extra=original_extra_rows(raw)
    ids=identity_audit(old,extra)
    data=hydrate(old,extra)
    replay=replay_v58(data)
    dictionary=source_dictionary(r)
    cases={}
    for policy in POLICIES:
        for ctx in CONTEXTS:
            if policy=="water3" and ctx=="no_week":continue
            cases[policy+"__"+ctx]=evaluate(data,policy,ctx)
    primary=cases["week4__linear_week"]
    ref=placebo(data,primary,r)
    directional=(primary["week_equal"]["gain"]>0
                 and primary["plant_weighted"]["gain"]>0
                 and primary["n_positive_weeks"]>=3)
    return {"analysis":"V59_RETROSPECTIVE_MERTENSIA_SAME_STUDY_CONTEXT_TRANSFER",
            "original_data_sha256":r["source_sha256"],
            "frozen_v58_result_git_blob":r["v58_result_git_blob"],
            "n_new_independent_studies":0,"n_plants":len(data),
            "source_identity_audit":ids,"author_dictionary":dictionary,
            "v58_replay_and_week_decomposition":replay,
            "original_vwc_by_water_code":{
                t:{"n":len([x for x in data if x["treat"]==t]),
                   "mean":sum(x["mVWC"] for x in data if x["treat"]==t)/
                          sum(x["treat"]==t for x in data),
                   "min":min(x["mVWC"] for x in data if x["treat"]==t),
                   "max":max(x["mVWC"] for x in data if x["treat"]==t)}
                for t in TREAT},
            "fixed_context_policy_results":cases,
            "primary_whole_week_directional_gate_pass":directional,
            "primary_whole_week_placebo_reference":ref,
            "primary_strong_descriptive_gate_pass":directional and ref["above_q975"],
            "primary_verdict":("DESCRIPTIVE_WHOLE_WEEK_VISITOR_INFORMATION" if directional and ref["above_q975"] else
                               "WHOLE_WEEK_VISITOR_TRANSFER_NOT_ESTABLISHED"),
            "new_causal_pollination_inference":False,
            "raw_plant_ids_or_raw_measurements_exported":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=run()
    if args.check:
        need(out==json.loads(RECEIPT.read_text()),
             "V59_ORIGINAL_CONTEXT_TRANSFER_RECEIPT_CHANGED")
        print("V59_ORIGINAL_SOURCE_REPLAY_ALL_CONTEXT_RESULTS_EXACT_PASS")
    else:
        raw=json.dumps(out,sort_keys=True,separators=(",",":"),ensure_ascii=True,allow_nan=False)
        n=(len(raw)+2999)//3000
        for i in range(n):
            print("V59_CONTEXT_READOUT_%03d_OF_%03d="%(i+1,n)+raw[3000*i:3000*(i+1)])
if __name__=="__main__":main()
