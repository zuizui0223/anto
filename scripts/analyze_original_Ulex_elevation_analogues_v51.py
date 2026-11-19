#!/usr/bin/env python3
"""ANTO v51 — original 2023 Ulex six-site environment-analogue prediction.

SIX ORIGINAL SITES ARE WITHIN ONE GENETICALLY CONTINUOUS REGIONAL POPULATION,
NOT six genetically independent populations. This is posthoc and neither
independent site experiments nor causal selection or future trait regain.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v51_ulex_original_site_registry_20261119.json"
ANCHOR=ROOT/"data/external/anto_v51_sdea_v15_derived_anchor_20261119.json"
RECEIPT=ROOT/"results/anto_v51_ulex_six_site_elevation_analogue_20261119.json"
FRUIT="fruits.and.mean.floral.traits.csv"
TRAIT="floral.traits.csv"
KIN="relatedness.matrix.csv"
POLL="pollinator.census.csv"
READ="README.md"
SITES=("AresM","Cheste","Chiva","Montserrat","Simat","Sot.Chera")
M0="M0"
M1="M1"
RIDGE=5.0


class V51Error(ValueError):
    pass


def require(ok, msg):
    if not ok:
        raise V51Error(msg)


def md5(data):
    return hashlib.md5(data).hexdigest()


def source_registry():
    r=json.loads(REGISTRY.read_text())
    require(
        r["analysis_type"]==
            "POSTHOC_ULEX_ORIGINAL_SIX_SAMPLING_SITES_ELEVATION_ANALOGUE_TRANSFER"
        and r["zenodo_record_id"]==7761289
        and r["zenodo_url"]=="https://zenodo.org/api/records/7761289"
        and r["source_md5"]=={
            FRUIT:"8ca72b7d9721f6d6279e793a1526c43e",
            TRAIT:"14216bb46af39bfd0b959c1aa092485b",
            KIN:"ecaa79680ffa25cb45ebdb73c11fe407",
            POLL:"78494af547de11c1800af50a74c5642f",
            READ:"9a6d62c907a93fe7b67268922f756995"}
        and r["sdea_v15_result_git_blob_sha1"]==
            "1fb22c1028388e720f388f2345c168270b065af4"
        and r["source_original_tagged_individuals"]==240
        and r["source_eligible_fruit_trait_individuals"]==212
        and r["source_relatedness_matrix_genotypes"]==225
        and r["source_matched_individuals"]==202
        and r["source_sites"]==
            {"AresM":32,"Cheste":37,"Chiva":33,
             "Montserrat":37,"Simat":37,"Sot.Chera":26}
        and r["source_pollinator_census_n"]==569
        and r["six_sites_are_separate_genetic_populations"] is False
        and r["sites_within_one_continuous_population"] is True
        and r["models"]==["M0","M1"]
        and r["ridge_alpha"]==5
        and r["site_train_regimes"]==["FULL5","NEAR2","FAR2","ALL2"]
        and r["n_all_pairs"]==60
        and r["n_site_holdouts"]==6
        and r["reference_parity_abs_tolerance"]==1e-7
        and r["n_independent_original_studies_added_vs_sdea_v15"]==0
        and r["source_is_posthoc"] is True
        and r["valid_randomization_p_value"] is False,
        "V51_REGISTERED_SOURCE_OR_INFERENCE_GRAIN_CHANGED")
    return r


def load_anchor(r):
    x=json.loads(ANCHOR.read_text())
    require(x["sdea_parent_original_receipt_blob_sha1"]==
            r["sdea_v15_result_git_blob_sha1"]
            and x["source_zenodo_record"]==7761289
            and x["source_md5"]==r["source_md5"]
            and x["source_matched_individuals"]==202
            and x["original_six_sampling_sites_not_genetically_independent"] is True
            and {q["site"] for q in x["original_site_folds"]}==set(SITES)
            and {q["site"]:q["n"] for q in x["original_site_folds"]}==
                r["source_sites"],"V51_ORIGINAL_SDEA_REFERENCE_CHANGED")
    return x


def allowed_zenodo(url):
    u=urlsplit(url)
    return (u.scheme=="https" and u.hostname=="zenodo.org"
            and u.username is None and u.password is None)


def download_source(r, *, source_root=None, timeout=90):
    if source_root is not None:
        folder=Path(source_root)
        out={name:(folder/name).read_bytes() for name in r["source_md5"]}
    else:
        try:
            with urlopen(Request(r["zenodo_url"],
                                 headers={"User-Agent":"ANTO-v51-original-Ulex-source-audit/1.0"}),
                         timeout=timeout) as resp:
                meta=json.load(resp)
        except Exception as exc:
            raise V51Error("V51_OFFICIAL_ZENODO_METADATA_UNAVAILABLE") from exc
        require(meta.get("id")==7761289,"V51_WRONG_ORIGINAL_ZENODO_RECORD")
        files=meta.get("files",[])
        if isinstance(files,dict):
            files=files.get("entries",files)
        if isinstance(files,dict):
            files=[dict(value,key=name) for name,value in files.items()]
        require(isinstance(files,list),"V51_UNEXPECTED_ZENODO_FILE_INVENTORY")
        inventory={str(f.get("key",f.get("filename"))):f for f in files}
        out={}
        for name,want in r["source_md5"].items():
            require(name in inventory,"V51_SOURCE_REQUIRED_FILE_MISSING_"+name)
            record=inventory[name]
            reported=str(record.get("checksum","")).lower().split(":")[-1]
            require(not reported or reported==want,
                    "V51_METADATA_MD5_DRIFT_"+name)
            links=record.get("links",{})
            url=links.get("content") or links.get("download") or links.get("self")
            require(isinstance(url,str) and allowed_zenodo(url),
                    "V51_UNTRUSTED_SOURCE_FILE_URL_"+name)
            max_bytes=1_000_000 if name==KIN else 150_000
            try:
                with urlopen(Request(url,headers={
                        "User-Agent":"ANTO-v51-original-Ulex-source-audit/1.0",
                        "Accept":"text/csv,text/plain,*/*"}),timeout=timeout) as resp:
                    raw=resp.read(max_bytes+1)
            except Exception as exc:
                raise V51Error("V51_OFFICIAL_SOURCE_DOWNLOAD_FAILED_"+name) from exc
            require(len(raw)<=max_bytes,"V51_ORIGINAL_SOURCE_TOO_LARGE_"+name)
            out[name]=raw
    require(set(out)==set(r["source_md5"]),"V51_SOURCE_FILE_SET_CHANGED")
    for name,expected in r["source_md5"].items():
        require(md5(out[name])==expected,
                "V51_ORIGINAL_SOURCE_BYTE_MD5_MISMATCH_"+name)
    return out


def csvrows(raw,header=None):
    try:
        reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
        names=tuple(reader.fieldnames or ())
        if header is not None:
            require(names==header,"V51_ORIGINAL_CSV_HEADER_DRIFT")
        records=list(reader)
    except (UnicodeDecodeError,csv.Error) as exc:
        raise V51Error("V51_ORIGINAL_SOURCE_CSV_PARSE") from exc
    require(all(None not in row and all(v is not None for v in row.values())
                for row in records),"V51_ORIGINAL_CSV_ROW_SHAPE_CHANGED")
    return names,records


def number(v,positive=False,nonnegative=False):
    if v is None or v.strip() in ("","NA",".","#VALUE!"):
        return None
    try:
        x=float(v)
    except ValueError as exc:
        raise V51Error("V51_SOURCE_NONNUMERIC_UNKNOWN_TOKEN") from exc
    require(math.isfinite(x) and (not positive or x>0)
            and (not nonnegative or x>=0),"V51_ORIGINAL_SOURCE_BAD_NUMBER")
    return x


def parse_original(payloads,r):
    for name,expected in r["source_md5"].items():
        require(name in payloads and md5(payloads[name])==expected,
                "V51_ORIGINAL_SOURCE_BYTES_UNVERIFIED")
    _,fruit=csvrows(payloads[FRUIT],
                   ("site","ind","elev","weight","area",
                    "scars","fruits","nofruit"))
    _,flowers=csvrows(payloads[TRAIT],
                     ("site","ind","flower","weight","area","elev"))
    poll_names,poll=csvrows(payloads[POLL],
                     ("year","locality","flowers","tot.visits",
                      "Apis","Bombus","other","visits.flower","observations"))
    require("ind" not in poll_names and "plant" not in poll_names
            and len(poll)==569,
            "V51_PATCH_POLLINATOR_CENSUS_INDIVIDUAL_JOIN_FORBIDDEN")
    cols,genome=csvrows(payloads[KIN])
    require(len(cols)==226 and cols[0]=="ind"
            and len(set(cols))==226 and len(genome)==225,
            "V51_ORIGINAL_GENOTYPED_MATRIX_SHAPE_CHANGED")
    k={row["ind"].strip() for row in genome}
    require(len(k)==225 and k==set(cols[1:]),
            "V51_GENETIC_SOURCE_INDIVIDUAL_IDENTIFIERS_CHANGED")
    require(len(fruit)==240 and len(flowers)==1125,
            "V51_ORIGINAL_PLANT_OR_FLOWER_DENOMINATOR_DRIFT")
    flowerkeys={(x["site"].strip(),x["ind"].strip()) for x in flowers}
    require(len({(x["site"],x["ind"]) for x in fruit})==240,
            "V51_DUPLICATE_ORIGINAL_FRUIT_INDIVIDUAL")
    eligible=[]
    source_key_sites=defaultdict(set)
    formula_error_n=0
    for row in fruit:
        site=row["site"].strip();ind=row["ind"].strip()
        require(site and ind,"V51_EMPTY_ORIGINAL_SITE_INDIVIDUAL_ID")
        source_key_sites[ind].add(site)
        formula_error_n+=sum(row[x].strip()=="#VALUE!" for x in
                             ("elev","weight","area","scars","fruits","nofruit"))
        vals={c:number(row[c],positive=c in ("area","weight"),
                       nonnegative=c in ("scars","fruits","nofruit"))
              for c in ("elev","weight","area","scars","fruits","nofruit")}
        if any(val is None for val in vals.values()):
            continue
        a=vals["scars"];b=vals["fruits"];c=vals["nofruit"]
        require(all(x.is_integer() for x in (a,b,c))
                and b<=a and math.isclose(a-b,c,abs_tol=1e-9),
                "V51_ORIGINAL_FRUIT_SCARS_ACCOUNTING_DRIFT")
        if a==0:
            continue
        eligible.append({"site":site,"ind":ind,"elev":vals["elev"],
                         "area":vals["area"],"weight":vals["weight"],
                         "fruit_set":b/a,"scars":int(a),"fruits":int(b)})
    require(all(len(v)==1 for v in source_key_sites.values()),
            "V51_ORIGINAL_GENOTYPE_REUSED_ACROSS_SITES")
    matched=[x for x in eligible if x["ind"] in k]
    counts=Counter(x["site"] for x in matched)
    require(len(eligible)==212 and len(matched)==202
            and dict(counts)==r["source_sites"]
            and formula_error_n==1,
            "V51_SDEA_V15_SOURCE_COHORT_OR_SITE_COUNTS_CHANGED")
    return matched, {
        "tagged_original_plants":240,
        "usable_original_fruit_trait_plants":len(eligible),
        "source_genotyped_individual_ids":len(k),
        "joined_original_plant_individuals":len(matched),
        "site_counts":dict(sorted(counts.items())),
        "original_repeated_flower_rows":len(flowers),
        "original_patch_census_rows_quarantined":len(poll),
        "source_spreadsheet_formula_error_cells":formula_error_n,
        "sampling_sites_are_distinct_genetic_populations":False,
        "six_sites_within_one_continuous_original_population":True,
        "new_biological_original_studies_vs_sdea_v15":0
    }


def solve(a,b):
    n=len(b)
    require(n>0 and len(a)==n and all(len(z)==n for z in a),
            "V51_BAD_NORMAL_EQUATION_SHAPE")
    aug=[list(a[i])+[b[i]] for i in range(n)]
    for col in range(n):
        pivot=max(range(col,n),key=lambda k:abs(aug[k][col]))
        require(abs(aug[pivot][col])>=1e-12,
                "V51_SINGULAR_ORIGINAL_RIDGE")
        aug[col],aug[pivot]=aug[pivot],aug[col]
        div=aug[col][col]
        aug[col]=[v/div for v in aug[col]]
        for row in range(n):
            if row==col:
                continue
            factor=aug[row][col]
            aug[row]=[v-factor*w for v,w in zip(aug[row],aug[col])]
    return [z[-1] for z in aug]


def stats(train,key):
    values=[r[key] for r in train]
    mu=sum(values)/len(values)
    sd=math.sqrt(sum((x-mu)**2 for x in values)/len(values))
    return mu,max(sd,1e-10)


def vec(row,info,traits):
    e=(row["elev"]-info["elev"][0])/info["elev"][1]
    x=[1.,e]
    if traits:
        a=(row["area"]-info["area"][0])/info["area"][1]
        w=(row["weight"]-info["weight"][0])/info["weight"][1]
        x.extend([a,w,a*a,w*w])
    return x


def fit(train,traits):
    require(len(train)>=30,"V51_NOT_ENOUGH_ORIGINAL_TRAIN_INDIVIDUALS")
    info={key:stats(train,key) for key in ("elev","area","weight")}
    vectors=[vec(r,info,traits) for r in train]
    y=[r["fruit_set"] for r in train]
    n=len(vectors[0])
    mat=[[0.]*n for _ in range(n)]
    rhs=[0.]*n
    for x,target in zip(vectors,y):
        for i in range(n):
            rhs[i]+=x[i]*target
            for j in range(n):
                mat[i][j]+=x[i]*x[j]
    for i in range(1,n):
        mat[i][i]+=RIDGE
    return {"beta":solve(mat,rhs),"stats":info,"traits":traits}


def predict(model,row):
    raw=sum(v*w for v,w in zip(
        model["beta"],vec(row,model["stats"],model["traits"])))
    return max(0.,min(1.,raw))


def mse_model(train,test):
    require(train and test
            and {r["site"] for r in train}.isdisjoint(
                {r["site"] for r in test})
            and {r["ind"] for r in train}.isdisjoint(
                {r["ind"] for r in test}),
            "V51_ORIGINAL_SITE_OR_INDIVIDUAL_LEAKAGE")
    return {
        "M0":sum((predict(fit_model0,r)-r["fruit_set"])**2
                 for r in test)/len(test),
        "M1":sum((predict(fit_model1,r)-r["fruit_set"])**2
                 for r in test)/len(test)
    } if False else _mse_original(train,test)


def _mse_original(train,test):
    f0=fit(train,False)
    f1=fit(train,True)
    return {"M0":sum((predict(f0,r)-r["fruit_set"])**2 for r in test)/len(test),
            "M1":sum((predict(f1,r)-r["fruit_set"])**2 for r in test)/len(test)}


def number_out(x):
    require(math.isfinite(x),"V51_NONFINITE_SOURCE_DERIVED_OUTPUT")
    return round(x,12)


def report(rows,audit,original_anchor,reg):
    site_rows={site:[r for r in rows if r["site"]==site] for site in SITES}
    site_mean={site:sum(r["elev"] for r in values)/len(values)
               for site,values in site_rows.items()}
    expected={x["site"]:x for x in original_anchor["original_site_folds"]}
    per_site=[]
    all_pairs=0
    max_anchor_error=0.0
    for heldout in SITES:
        train_sites=[s for s in SITES if s!=heldout]
        test=site_rows[heldout]
        fulltrain=[r for site in train_sites for r in site_rows[site]]
        # The published v15 original estimator is preserved in full.
        full=_mse_original(fulltrain,test)
        reference=expected[heldout]
        for key,col in (("M0","m0"),("M1","m1")):
            delta=abs(full[key]-reference[col])
            max_anchor_error=max(max_anchor_error,delta)
            require(delta<=reg["reference_parity_abs_tolerance"],
                    "V51_SDEA_V15_ORIGINAL_FIVE_SITE_PARITY_FAILED_"+heldout+"_"+key)
        by_distance=sorted(train_sites,key=lambda s:
            (abs(site_mean[s]-site_mean[heldout]),s))
        near=tuple(sorted(by_distance[:2]))
        far=tuple(sorted(sorted(train_sites,key=lambda s:
            (-abs(site_mean[s]-site_mean[heldout]),s))[:2]))
        pairdata=[]
        for pair in itertools.combinations(sorted(train_sites),2):
            selected=[r for site in pair for r in site_rows[site]]
            mse=_mse_original(selected,test)
            pairdata.append({
                "training_sites":list(pair),
                "training_n_original_individuals":len(selected),
                "elevation_mean_gap_average":number_out(
                    sum(abs(site_mean[s]-site_mean[heldout]) for s in pair)/2),
                "M0":number_out(mse["M0"]),
                "M1":number_out(mse["M1"]),
                "trait_mse_increment_M0_minus_M1":number_out(mse["M0"]-mse["M1"]),
            })
            all_pairs+=1
        require(len(pairdata)==10,"V51_MISSING_TEN_ORIGINAL_SITE_PAIRS")
        def choose(pair):
            hit=[x for x in pairdata if tuple(x["training_sites"])==pair]
            require(len(hit)==1,"V51_NONUNIQUE_ENVIRONMENT_PAIR")
            return hit[0]
        best=choose(near);worst=choose(far)
        sorted_risks=sorted(x["M1"] for x in pairdata)
        rank=1+sum(x<best["M1"]-1e-12 for x in sorted_risks)
        far_rank=1+sum(x<worst["M1"]-1e-12 for x in sorted_risks)
        per_site.append({
            "heldout_site":heldout,
            "n_original_heldout_genotyped_plants":len(test),
            "heldout_site_mean_elevation":number_out(site_mean[heldout]),
            "FULL5":{"training_sites":sorted(train_sites),
                     "training_n_original_individuals":len(fulltrain),
                     "M0":number_out(full["M0"]),"M1":number_out(full["M1"]),
                     "trait_mse_increment_M0_minus_M1":number_out(full["M0"]-full["M1"])},
            "NEAR2":dict(best,rank_among_ten_two_site_pairs_M1=rank),
            "FAR2":dict(worst,rank_among_ten_two_site_pairs_M1=far_rank),
            "all_ten_source_train_pair_scores":pairdata,
        })
    require(len(per_site)==6 and all_pairs==60,
            "V51_EXPECT_SIX_SITES_AND_SIXTY_ORIGINAL_TRAIN_PAIRS")
    main={}
    for regime in ("FULL5","NEAR2","FAR2"):
        a=[(q["n_original_heldout_genotyped_plants"],q[regime]) for q in per_site]
        n_tot=sum(n for n,_ in a)
        m0=sum(x["M0"] for _,x in a)/6
        m1=sum(x["M1"] for _,x in a)/6
        w0=sum(n*x["M0"] for n,x in a)/n_tot
        w1=sum(n*x["M1"] for n,x in a)/n_tot
        require(m0>0 and w0>0,"V51_BAD_SITE_BASELINE_DENOMINATOR")
        main[regime]={
            "n_heldout_original_sampling_sites":6,
            "n_unique_original_individuals_heldout_once_in_regime":n_tot,
            "equal_site_M0_MSE":number_out(m0),
            "equal_site_M1_MSE":number_out(m1),
            "equal_site_absolute_trait_MSE_gain":number_out(m0-m1),
            "equal_site_relative_trait_MSE_gain":number_out(1-m1/m0),
            "individual_weighted_relative_trait_MSE_gain":number_out(1-w1/w0),
            "sites_with_M1_better_than_M0":sum(x["M1"]<x["M0"] for _,x in a),
        }
    return {
        "analysis_type":reg["analysis_type"],
        "source_original_zenodo_record":7761289,
        "source_original_file_MD5":reg["source_md5"],
        "parent_SDEA_v15_derived_result_git_blob_sha1":
            reg["sdea_v15_result_git_blob_sha1"],
        "n_independent_Ulex_original_studies_added_vs_SDEA_v15":0,
        "n_distinct_original_population_genetic_replicates_in_2023_Ulex":1,
        "six_sites_within_one_continuous_original_population":True,
        "original_source_audit":audit,
        "target_site_elevation_is_used_transductively_without_target_fruit_outcome":True,
        "all_site_group_heldout_no_individual_replicate_pseudoreplication":True,
        "pollinator_patch_level_data_never_plant_joined":True,
        "all_models_original_fixed_ridge_penalty":5,
        "original_SDEA_v15_full5_max_absolute_MSE_parity_error":max_anchor_error,
        "n_original_site_holdouts":6,
        "n_exact_train_two_site_pair_fits":60,
        "site_mean_elevation_by_original_site":{site:number_out(site_mean[site])
                                             for site in SITES},
        "per_original_site_complete_scores":per_site,
        "summary_by_training_support":main,
        "near2_vs_full5_M1_sites_better":
            sum(q["NEAR2"]["M1"]<q["FULL5"]["M1"] for q in per_site),
        "far2_vs_full5_M1_sites_better":
            sum(q["FAR2"]["M1"]<q["FULL5"]["M1"] for q in per_site),
        "new_causal_genetic_selection_or_ancestral_flower_reexpression":False,
        "randomization_p_values_or_confidence_intervals_calculated":False,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-root",default=None)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    r=source_registry()
    payload=download_source(r,source_root=args.source_root)
    rows,audit=parse_original(payload,r)
    result=report(rows,audit,load_anchor(r),r)
    if args.check:
        require(result==json.loads(RECEIPT.read_text()),
                "V51_SOURCE_REPLAY_RECEIPT_DIFFERENT")
        print("V51_ORIGINAL_ZENODO_SOURCE_AND_SDEA_V15_PARITY_EXACT_PASS")
    else:
        raw=json.dumps(result,sort_keys=True,separators=(",",":"),
                       ensure_ascii=True)
        n=(len(raw)+3999)//4000
        for i in range(n):
            print("V51_SOURCE_READOUT_CHUNK_%03d_OF_%03d="%(
                i+1,n)+raw[i*4000:(i+1)*4000])


if __name__=="__main__":
    main()
