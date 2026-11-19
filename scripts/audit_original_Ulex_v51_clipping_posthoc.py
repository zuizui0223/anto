#!/usr/bin/env python3
"""Post-outcome Ulex v51 clipping audit: not a new selection result."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import analyze_original_Ulex_elevation_analogues_v51 as v

ROOT=Path(__file__).resolve().parents[1]
REG=ROOT/"data/external/anto_v51_posthoc_clipping_registry_20261119.json"
PARENT=ROOT/"results/anto_v51_ulex_six_site_elevation_analogue_20261119.json"
RECEIPT=ROOT/"results/anto_v51_posthoc_Ulex_clipping_audit_20261119.json"
REGIMES=("FULL5","NEAR2","FAR2")
MODELS=("M0","M1")

class ClippingError(ValueError): pass
def guard(x,msg):
    if not x: raise ClippingError(msg)

def git_blob_sha1(raw):
    return hashlib.sha1(("blob %d\0"%len(raw)).encode()+raw).hexdigest()

def lock():
    x=json.loads(REG.read_text())
    guard(x["analysis_type"]=="POST_OUTCOME_ULEX_V51_FAR_PAIR_PREDICTION_CLIP_DIAGNOSTIC"
          and x["source_original_v51_result_git_blob_sha1"]==
              "8f1d5435d2d14c2b2c999671ff9801cffab01a71"
          and x["original_zenodo_record_id"]==7761289
          and x["original_source_individuals"]==202
          and x["original_sampling_sites"]==6
          and x["regimes"]==list(REGIMES)
          and x["models"]==list(MODELS)
          and x["original_ridge_alpha"]==5
          and x["full_parent_mse_parity_abs_tolerance"]==1e-10
          and x["retrospective_posthoc"] is True
          and x["independently_sampled_original_studies_added"]==0
          and x["new_original_genetic_populations_added"]==0
          and x["external_pollinator_individual_joined"] is False
          and x["valid_randomization_p_value"] is False,
          "V51_CLIPPING_POSTHOC_SOURCE_OR_SCOPE_DRIFT")
    return x

def read_parent(r):
    raw=PARENT.read_bytes()
    guard(git_blob_sha1(raw)==r["source_original_v51_result_git_blob_sha1"],
          "V51_CLIPPING_PARENT_RESULT_BYTES_CHANGED")
    x=json.loads(raw)
    guard(x["n_exact_train_two_site_pair_fits"]==60
          and x["original_source_audit"]["joined_original_plant_individuals"]==202
          and x["six_sites_within_one_continuous_original_population"] is True
          and len(x["per_original_site_complete_scores"])==6,
          "V51_CLIPPING_PARENT_BIOLOGICAL_GRAIN_CHANGED")
    return x

def preclip(model,row):
    return sum(a*b for a,b in zip(model["beta"],
                                v.vec(row,model["stats"],model["traits"])))

def source_audit(rows, parent, r):
    index={site:[row for row in rows if row["site"]==site] for site in v.SITES}
    sites={x["heldout_site"]:x for x in parent["per_original_site_complete_scores"]}
    output={"analysis_type":r["analysis_type"],
            "source_original_v51_result_git_blob_sha1":
                r["source_original_v51_result_git_blob_sha1"],
            "n_original_individuals":202,
            "original_sites_one_continuous_population":True,
            "source_fitted_models_only":list(MODELS),
            "source_train_regimes_only":list(REGIMES),
            "same_original_folds_reused":True,
            "n_additional_independent_biological_studies":0,
            "derived_original_site_diagnostics":[],
            "aggregated_all_sites":{},
            "posthoc_not_primary_confirmation":True}
    for site in v.SITES:
        rr=sites[site]
        for regime in REGIMES:
            source=rr[regime]
            train_sites=source["training_sites"]
            guard(site not in train_sites and
                  len(train_sites)==(5 if regime=="FULL5" else 2),
                  "V51_CLIPPING_FOLD_SOURCE_SUPPORT_DRIFT")
            train=[row for s in train_sites for row in index[s]]
            test=index[site]
            guard(len(test)==rr["n_original_heldout_genotyped_plants"] and
                  len(train)==source["training_n_original_individuals"],
                  "V51_CLIPPING_ORIGINAL_INDIVIDUAL_DENOMINATOR_DRIFT")
            stats={}
            values={}
            for model in MODELS:
                fit=v.fit(train,model=="M1")
                before=[preclip(fit,row) for row in test]
                after=[min(1,max(0,x)) for x in before]
                loss=sum((out-row["fruit_set"])**2
                         for out,row in zip(after,test))/len(test)
                guard(abs(loss-source[model])<=r["full_parent_mse_parity_abs_tolerance"],
                      "V51_CLIPPING_SOURCE_LOSS_MUST_REPLAY")
                stats[model]={"n_original_heldout_individuals":len(test),
                              "raw_below_zero_n":sum(x<0 for x in before),
                              "raw_above_one_n":sum(x>1 for x in before),
                              "raw_inside_closed_unit_interval_n":
                                  sum(0<=x<=1 for x in before),
                              "raw_min":round(min(before),12),
                              "raw_max":round(max(before),12),
                              "clipped_original_fruit_fraction_mse":round(loss,12)}
                values[model]=after
            equal=sum(abs(a-b)<=1e-12 for a,b in zip(values["M0"],values["M1"]))
            both_at_same_edge=sum((a==0 and b==0) or (a==1 and b==1)
                                  for a,b in zip(values["M0"],values["M1"]))
            output["derived_original_site_diagnostics"].append({
                "site":site,"regime":regime,
                "original_training_site_count":len(train_sites),
                "original_heldout_individual_count":len(test),
                "model_clip_counts":stats,
                "same_clipped_prediction_M0_M1_n":equal,
                "both_models_at_same_zero_or_one_boundary_n":both_at_same_edge})
    for regime in REGIMES:
        group=[d for d in output["derived_original_site_diagnostics"]
               if d["regime"]==regime]
        guard(len(group)==6 and
              sum(d["original_heldout_individual_count"] for d in group)==202,
              "V51_CLIPPING_NOT_SIX_SITES_202_INDIVIDUALS")
        output["aggregated_all_sites"][regime]={
            "n_original_heldout_individuals_once_each":202,
            "by_model":{
                m:{
                    "n_raw_predictions_below_zero":
                        sum(d["model_clip_counts"][m]["raw_below_zero_n"] for d in group),
                    "n_raw_predictions_above_one":
                        sum(d["model_clip_counts"][m]["raw_above_one_n"] for d in group),
                    "n_raw_predictions_within_unit":
                        sum(d["model_clip_counts"][m]["raw_inside_closed_unit_interval_n"]
                            for d in group)
                } for m in MODELS
            },
            "n_same_after_clipping_M0_and_M1":
                sum(d["same_clipped_prediction_M0_M1_n"] for d in group),
            "n_both_models_at_same_clipping_boundary":
                sum(d["both_models_at_same_zero_or_one_boundary_n"] for d in group)
        }
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",default=None)
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    x=lock()
    source=v.source_registry()
    matched,audit=v.parse_original(v.download_source(
        source,source_root=args.source_root),source)
    guard(len(matched)==202 and audit["six_sites_within_one_continuous_original_population"],
          "V51_CLIPPING_WRONG_ORIGINAL_SOURCE_GRAIN")
    output=source_audit(matched,read_parent(x),x)
    if args.check:
        guard(output==json.loads(RECEIPT.read_text()),
              "V51_CLIPPING_SOURCE_RESULT_REGENERATION_CHANGED")
        print("V51_POSTHOC_CLIPPING_6_SITE_SOURCE_REPLAY_EXACT_PASS")
    else:
        print("V51_CLIPPING_SOURCE_READOUT="+json.dumps(
            output,sort_keys=True,separators=(",",":"),ensure_ascii=True))

if __name__=="__main__":
    main()
