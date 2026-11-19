#!/usr/bin/env python3
"""ANTO v52: original Sinapis rainfall-treatment-heldout visitor prediction.

POSTHOC new holdout domain, not new ecological specimens or causal mediation.
Imports SHA-authenticated SDEA v17 original-source join + historical comparator.
No original author's plant or visitor records are committed to ANTO.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT/"data/external/anto_v52_sinapis_treatment_original_source_registry_20261119.json"
PRIOR = ROOT/"data/external/anto_v52_sdea_v17_derived_original_anchor_20261119.json"
RECEIPT = ROOT/"results/anto_v52_sinapis_treatment_heldout_original_source_20261119.json"
TREATMENTS = ("Control","Drought","PrecipRed")
SHELTERS = ("1","2","3","4","5")
MODELS = ("U0","U1","U2","U3")
FEATURES = {
    "U0":("log_flower",),
    "U1":("log_flower","display","petal_width"),
    "U2":("log_flower","display","petal_width","log_visit_rate"),
    "U3":("log_flower","log_visit_rate"),
}
RIDGE = 4.0


class V52Error(ValueError):
    pass


def require(valid, code):
    if not valid:
        raise V52Error(code)


def git_blob_sha1(data:bytes)->str:
    return hashlib.sha1(("blob %d\0"%len(data)).encode("ascii")+data).hexdigest()


def read_registry():
    r=json.loads(REGISTRY.read_text())
    require(
        r["analysis_type"]==
            "POSTHOC_ORIGINAL_SINAPIS_THREE_TREATMENT_HELDOUT_VISITOR_INFORMATION"
        and r["original_publication_doi"]=="10.1111/1365-2435.14233"
        and r["zenodo_record_id"]==7335086
        and r["author_source_MD5"]=={
            "Floral_traits.csv":"24f1f01435fa22faa832c415bbd0b77e",
            "Int_data.csv":"d67c32cc5aa357fe6e6d952a801ea38f",
            "Plantlist.csv":"630eb7d3d7c871f65c2bd62890e1597c",
            "README.md":"68bf11faaa9786fbb39f9e32006b336d"}
        and r["sdea_repo"]=="zuizui0223/SDEA"
        and r["sdea_original_source_commit"]==
            "a8990f5255b1b27ed5fec85587bf09ef6390882a"
        and r["sdea_original_parser_git_blob_sha1"]==
            "5f797c1d891c8f8de711909a56ae26530841a2ba"
        and r["sdea_v17_prior_derived_receipt_git_blob_sha1"]==
            "fa2419d32b74f22f67ea9cbfd79123109db4fdf0"
        and r["original_total_plant_records"]==53
        and r["original_visitor_group_rows"]==427
        and r["original_complete_case_plants"]==52
        and r["original_shelter_counts"]==
            {"1":11,"2":9,"3":8,"4":12,"5":12}
        and r["original_treatments"]==list(TREATMENTS)
        and r["original_rain_shelters"]==5
        and r["treatment_holdout_folds"]==3
        and r["model_names"]==list(MODELS)
        and all(r[m]==list(FEATURES[m]) for m in MODELS)
        and r["model_treatment_indicators_excluded"] is True
        and r["training_ridge_penalty"]==4
        and r["training_only_standardization"] is True
        and r["original_sdea_shelter_20_cell_mse_abs_tolerance"]==1e-8
        and r["placebo_seed"]==20261119
        and r["placebo_replicates"]==200
        and r["pass_requirements"]["improved_treatments_at_least"]==2
        and r["new_original_biological_studies_vs_sdea_v17"]==0
        and r["randomization_causal_p_value"] is False
        and r["source_is_posthoc"] is True
        and r["no_molecular_genotype"] is True
        and r["no_pollen_deposition_measured"] is True
        and r["no_ancestral_flower_regain_measured"] is True,
        "V52_SOURCE_PROTOCOL_OR_BIOLOGICAL_SCOPE_DRIFT")
    return r


def load_frozen_sdea(sdea_script,registry):
    require(sdea_script is not None,"V52_PINNED_SDEA_SOURCE_PARSER_REQUIRED")
    raw=Path(sdea_script).read_bytes()
    require(git_blob_sha1(raw)==registry["sdea_original_parser_git_blob_sha1"],
            "V52_ORIGINAL_SDEA_SOURCE_IMPLEMENTATION_CHANGED")
    spec=importlib.util.spec_from_file_location("sdea_v17_immutable_source",str(sdea_script))
    require(spec is not None and spec.loader is not None,
            "V52_SDEA_SOURCE_MODULE_UNREADABLE")
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    require(mod.RECORD_ID==7335086
            and mod.ORIGINAL_MD5==registry["author_source_MD5"]
            and tuple(mod.TREATMENTS)==TREATMENTS
            and tuple(mod.SHELTERS)==SHELTERS
            and mod.RIDGE==RIDGE
            and mod.N_ORIGINAL_PLANTS==53
            and mod.N_ORIGINAL_VISIT_GROUP_ROWS==427,
            "V52_SDEA_SOURCE_ORIGINAL_CONSTANTS_CHANGED")
    return mod


def prior_source():
    p=json.loads(PRIOR.read_text())
    r=read_registry()
    require(p["parent_git_blob_sha1"]==
            r["sdea_v17_prior_derived_receipt_git_blob_sha1"]
            and p["source_zenodo_record"]==7335086
            and p["original_file_MD5"]==r["author_source_MD5"]
            and p["source_n_plants"]==53
            and p["source_n_visitor_group_rows"]==427
            and p["source_n_complete_plant_rows"]==52
            and p["source_shelters"]==r["original_shelter_counts"]
            and len(p["original_shelter_heldout_folds"])==5
            and p["previous_original_process_readiness"]==
                "PROCESS_ADDS_PORTABLE_INFORMATION_IN_FIVE_SHELTERS"
            and p["this_is_prior_replay_not_new_evidence"] is True
            and p["individual_original_source_records_included"] is False,
            "V52_FROZEN_SDEA_V17_PARENT_RECORD_CHANGED")
    return p


def source_payload(mod,r,source_root=None):
    if source_root is None:
        raw=mod.get_sources()
    else:
        raw={name:(Path(source_root)/name).read_bytes()
             for name in r["author_source_MD5"]}
    require(set(raw)==set(r["author_source_MD5"]),
            "V52_ORIGINAL_SOURCE_FOUR_FILES_CHANGED")
    mod.validate_bytes(raw)
    for name,digest in r["author_source_MD5"].items():
        require(hashlib.md5(raw[name]).hexdigest()==digest,
                "V52_ORIGINAL_AUTHOR_FILE_MD5_CHANGED_"+name)
    return raw


def original_cohort(mod,payload,r,p):
    mod.validate_bytes(payload)
    rows,audit=mod.source_join(payload)
    require(len(rows)==52
            and audit["plantlist_source_rows"]==53
            and audit["floral_source_rows"]==53
            and audit["visitor_group_source_rows"]==427
            and audit["individuals_with_all_model_fields"]==52
            and audit["observed_shelter_counts"]==
                r["original_shelter_counts"]
            and audit["plant_genetic_alleles_or_ancestors_measured"] is False
            and audit["visitor_observation_time_reused_once_per_plant"] is True,
            "V52_ORIGINAL_COMPLETE_PLANT_COHORT_DRIFT")
    treat=Counter(row["treatment"] for row in rows)
    strata=Counter((row["shelter"],row["treatment"]) for row in rows)
    require(set(treat)==set(TREATMENTS)
            and sum(treat.values())==52
            and min(treat.values())>=10
            and set(strata)=={(s,t) for s in SHELTERS for t in TREATMENTS},
            "V52_ORIGINAL_TREATMENT_SUPPORT_OR_SHELTER_CELLS_CHANGED")
    # Replay EXACT old SDEA v17 5 shelter-fold model architecture,
    # including original rain-treatment dummy indicators.
    old=mod.score_shelters(rows)
    require(old["total_heldout_individuals_once_each"]==52,
            "V52_PREVIOUS_SHELTER_HOLDOUT_GRAIN_CHANGED")
    by={str(q["shelter"]):q for q in p["original_shelter_heldout_folds"]}
    largest=0.
    for entry in old["folds"]:
        reference=by[entry["shelter"]]
        require(entry["heldout_plants"]==reference["n"]
                and set(entry["mse_models"])=={"M0","M1","M2","M3"},
                "V52_FROZEN_SDEA_SHELTER_FOLDS_CHANGED")
        for model in ("M0","M1","M2","M3"):
            diff=abs(entry["mse_models"][model]-
                     reference["models"][model])
            largest=max(largest,diff)
            require(diff<=r["original_sdea_shelter_20_cell_mse_abs_tolerance"],
                    "V52_ORIGINAL_SDEA_V17_SHELTER_MODEL_MSE_DRIFT")
    audit_out={
        "original_author_plant_rows":53,
        "original_visitor_group_rows":427,
        "original_complete_plant_rows":52,
        "n_original_rain_shelters":5,
        "shelter_eligible_plant_counts":dict(sorted(
            (str(s),n) for s,n in
            Counter(row["shelter"] for row in rows).items())),
        "treatment_eligible_plant_counts":dict(sorted(treat.items())),
        "shelter_treatment_eligible_cell_sizes":{
            s:{t:strata[(s,t)] for t in TREATMENTS} for s in SHELTERS},
        "original_sdea_v17_5_shelter_20_mse_max_absolute_error":largest,
        "same_original_biological_study_not_independent_replicate":True,
        "original_source_visitor_effort_once_per_plant":True,
        "original_source_molecular_genotype_measured":False,
    }
    return rows,audit_out


def solve(a,b):
    n=len(b)
    require(n>0 and len(a)==n and all(len(row)==n for row in a),
            "V52_INVALID_RIDGE_NORMAL_EQUATION")
    m=[list(a[i])+[b[i]] for i in range(n)]
    for j in range(n):
        pivot=max(range(j,n),key=lambda k:abs(m[k][j]))
        require(abs(m[pivot][j])>=1e-11,
                "V52_SINGULAR_TRAINING_RIDGE")
        m[j],m[pivot]=m[pivot],m[j]
        div=m[j][j]
        m[j]=[x/div for x in m[j]]
        for i in range(n):
            if i==j:continue
            mult=m[i][j]
            m[i]=[x-mult*y for x,y in zip(m[i],m[j])]
    x=[row[-1] for row in m]
    require(all(math.isfinite(v) for v in x),
            "V52_NONFINITE_RIDGE_ESTIMATE")
    return x


def fit(train,model):
    require(model in MODELS and len(train)>=20,
            "V52_UNREGISTERED_MODEL_OR_SMALL_TRAINING")
    features=FEATURES[model]
    transformed={}
    for f in features:
        vals=[float(row[f]) for row in train]
        mean=sum(vals)/len(vals)
        sd=math.sqrt(sum((x-mean)**2 for x in vals)/len(vals))
        require(math.isfinite(sd),
                "V52_NONFINITE_TRAIN_PC_VARIANCE")
        transformed[f]=(mean,max(1e-12,sd))
    def vector(row):
        return [1.]+[(row[f]-transformed[f][0])/transformed[f][1]
                     for f in features]
    xx=[vector(row) for row in train]
    yy=[float(row["log_seed"]) for row in train]
    dim=len(xx[0])
    A=[[0.]*dim for _ in range(dim)]
    b=[0.]*dim
    for x,y in zip(xx,yy):
        for i in range(dim):
            b[i]+=x[i]*y
            for j in range(dim):
                A[i][j]+=x[i]*x[j]
    for i in range(1,dim):
        A[i][i]+=RIDGE
    return {"model":model,"stats":transformed,"beta":solve(A,b),
            "n_train":len(train),"treatment_indicators_included":False}


def predict(fitted,row):
    features=FEATURES[fitted["model"]]
    x=[1.]+[(row[f]-fitted["stats"][f][0])/fitted["stats"][f][1]
             for f in features]
    return sum(a*b for a,b in zip(x,fitted["beta"]))


def treatment_folds(rows, models=MODELS):
    original=Counter(r["treatment"] for r in rows)
    require(len(rows)==52
            and set(original)==set(TREATMENTS)
            and sum(original.values())==52,
            "V52_ORIGINAL_TREATMENT_STUDY_GRAIN_CHANGED")
    out=[]
    for hold in TREATMENTS:
        train=[r for r in rows if r["treatment"]!=hold]
        test=[r for r in rows if r["treatment"]==hold]
        require(set(r["treatment"] for r in train)==set(TREATMENTS)-{hold}
                and len(train)+len(test)==52
                and len(train)>=20 and len(test)>=10,
                "V52_UNSEEN_TREATMENT_FOLD_TRAIN_TEST_LEAKAGE")
        fitted={m:fit(train,m) for m in models}
        loss={}
        for m in models:
            loss[m]=sum((predict(fitted[m],r)-r["log_seed"])**2
                         for r in test)/len(test)
            require(math.isfinite(loss[m]) and loss[m]>=0,
                    "V52_UNSEEN_TREATMENT_INVALID_MSE")
        out.append({"heldout_treatment":hold,
                    "n_heldout_original_plants":len(test),
                    "n_training_original_plants":len(train),
                    "original_train_treatments":sorted(
                        set(r["treatment"] for r in train)),
                    "original_heldout_shelters":dict(sorted(
                        Counter(r["shelter"] for r in test).items())),
                    "mse":loss})
    return out


def summarise(folds):
    require(len(folds)==3
            and {x["heldout_treatment"] for x in folds}==set(TREATMENTS)
            and sum(x["n_heldout_original_plants"] for x in folds)==52,
            "V52_MISSING_THREE_FULL_TREATMENT_FOLDS")
    a={}
    for m in MODELS:
        eq=sum(f["mse"][m] for f in folds)/3
        wt=sum(f["n_heldout_original_plants"]*f["mse"][m]
                for f in folds)/52
        a[m]={"equal_treatment_mse":eq,
              "plant_weighted_mse":wt}
    diffs={}
    for name,a_model,b_model in (
        ("morphology_U0_to_U1","U0","U1"),
        ("visitation_conditional_on_morphology_U1_to_U2","U1","U2"),
        ("visitation_without_morphology_U0_to_U3","U0","U3")
    ):
        gains=[f["mse"][a_model]-f["mse"][b_model] for f in folds]
        n=[f["n_heldout_original_plants"] for f in folds]
        diffs[name]={
            "equal_treatment_absolute_mse_gain":sum(gains)/3,
            "plant_weighted_absolute_mse_gain":
                sum(w*z for w,z in zip(n,gains))/52,
            "positive_heldout_treatments":sum(z>0 for z in gains),
            "per_treatment_signed_mse_gains":dict(zip(TREATMENTS,gains)),
        }
    return {"models":a,"increments":diffs}


def quantile(xs,prob):
    values=sorted(xs)
    position=(len(values)-1)*prob
    i=int(position)
    return values[i]*(1-(position-i))+values[min(i+1,len(values)-1)]*(position-i)


def visit_identity_placebo(rows, real_folds, r):
    require(r["placebo_replicates"]==200 and
            r["placebo_seed"]==20261119,
            "V52_VISITOR_PLACEBO_SETTINGS_CHANGED")
    original_summary=summarise(real_folds)
    goal=original_summary["increments"][
        "visitation_conditional_on_morphology_U1_to_U2"][
        "equal_treatment_absolute_mse_gain"]
    u1_reference={f["heldout_treatment"]:f["mse"]["U1"]
                  for f in real_folds}
    strata=defaultdict(list)
    for i,row in enumerate(rows):
        strata[(row["shelter"],row["treatment"])].append(i)
    require(set(strata)=={(s,t) for s in SHELTERS for t in TREATMENTS},
            "V52_ORIGINAL_PLACEBO_SHELTER_TREATMENT_STRATA_CHANGED")
    rng=random.Random(20261119)
    repeats=[]
    for _ in range(200):
        shuffled=[dict(z) for z in rows]
        for key in sorted(strata):
            idx=strata[key]
            original_values=[rows[i]["log_visit_rate"] for i in idx]
            rng.shuffle(original_values)
            for i,val in zip(idx,original_values):
                shuffled[i]["log_visit_rate"]=val
        scores=treatment_folds(shuffled,models=("U2",))
        gains=[u1_reference[f["heldout_treatment"]]-f["mse"]["U2"]
               for f in scores]
        repeats.append(sum(gains)/3)
    require(len(repeats)==200
            and all(math.isfinite(x) for x in repeats),
            "V52_INCOMPLETE_OR_BAD_SOURCE_PLACEBO")
    return {
        "permutation_replicates":200,
        "seed":20261119,
        "unit":"original plant visit rate shuffled only within original shelter×treatment source cell",
        "same_plant_seed_and_floral_data_unchanged":True,
        "observed_equal_treatment_gain":goal,
        "placebo_min":min(repeats),
        "placebo_median":quantile(repeats,.5),
        "placebo_q025":quantile(repeats,.025),
        "placebo_q975":quantile(repeats,.975),
        "placebo_max":max(repeats),
        "n_placebo_gains_at_or_above_observed":
            sum(x>=goal for x in repeats),
        "valid_randomized_causal_p_value":False,
    }


def report(rows, audit, parent, reg):
    result=treatment_folds(rows)
    summary=summarise(result)
    placebo=visit_identity_placebo(rows,result,reg)
    inc=summary["increments"][
        "visitation_conditional_on_morphology_U1_to_U2"]
    positive=(inc["equal_treatment_absolute_mse_gain"]>0
              and inc["plant_weighted_absolute_mse_gain"]>0
              and inc["positive_heldout_treatments"]>=2
              and inc["equal_treatment_absolute_mse_gain"]>
                   placebo["placebo_q975"])
    return {
        "analysis_type":reg["analysis_type"],
        "source_original_zenodo_record":7335086,
        "original_file_md5":reg["author_source_MD5"],
        "source_sdea_v17_parser_git_blob_sha1":
            reg["sdea_original_parser_git_blob_sha1"],
        "source_sdea_v17_derived_git_blob_sha1":
            reg["sdea_v17_prior_derived_receipt_git_blob_sha1"],
        "source_audit":audit,
        "previous_sdea_v17_shelter_holdout_replayed":True,
        "previous_original_sdea_v17_result_is_new_study":False,
        "n_new_original_biological_studies_vs_sdea_v17":0,
        "original_study_year":"2020",
        "model_treatment_labels_excluded_because_level_heldout":True,
        "model_training_ridge_alpha":4,
        "models":{k:list(FEATURES[k]) for k in MODELS},
        "new_holdout_axis":"original rainfall treatment fully withheld across all 5 shelters",
        "original_three_complete_treatment_folds":result,
        "heldout_rainfall_treatment_summary":summary,
        "source_visit_identity_placebo":placebo,
        "fixed_exploratory_process_transfer_verdict":
            ("DESCRIPTIVE_WITHIN_SINGLE_TRIAL_CROSS_TREATMENT_VISITOR_INFORMATION"
             if positive else
             "CROSS_TREATMENT_VISITOR_TRANSFER_NOT_ESTABLISHED"),
        "full_gene_pollinator_fitness_ancestral_regain_mechanism_identified":False,
        "new_causal_p_values_or_confidence_intervals":False,
        "source_individual_plant_records_redistributed":False,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--sdea-script",type=Path,required=True)
    p.add_argument("--source-root",type=Path)
    p.add_argument("--fetch-source-to",type=Path)
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    r=read_registry()
    mod=load_frozen_sdea(args.sdea_script,r)
    if args.fetch_source_to is not None:
        require(args.source_root is None and not args.check,
                "V52_FETCH_STAGING_NOT_AN_ANALYSIS_RUN")
        payload=source_payload(mod,r)
        args.fetch_source_to.mkdir(parents=True,exist_ok=True)
        for name,raw in payload.items():
            (args.fetch_source_to/name).write_bytes(raw)
        print("V52_STAGED_FOUR_EXACT_ZENODO_MD5_FILES")
        return
    pprior=prior_source()
    raw=source_payload(mod,r,source_root=args.source_root)
    rows,audit=original_cohort(mod,raw,r,pprior)
    out=report(rows,audit,pprior,r)
    if args.check:
        require(out==json.loads(RECEIPT.read_text()),
                "V52_ORIGINAL_SINAPIS_TREATMENT_RECEIPT_REPLAY_CHANGED")
        print("V52_ORIGINAL_MD5_SDEA_20_MSE_AND_TREATMENT_PLACEBO_RECEIPT_PASS")
    else:
        print("V52_ORIGINAL_SOURCE_READOUT="+json.dumps(
            out,sort_keys=True,separators=(",",":"),ensure_ascii=True))


if __name__=="__main__":
    main()
