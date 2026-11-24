#!/usr/bin/env python3
"""Official Dryad documented dataset ZIP route, Physalis author byte gate.

Previous per-file Dryad API URLs returned HTTP401. This checks a SECOND
official public API route, not a mirror. If STOP, do not fabricate evidence.
"""
from __future__ import annotations
import hashlib
import io
import json
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
import discover_original_Physalis_Mertensia_v57 as dryad
import audit_original_Physalis_Mertensia_file_grain_v57 as original

RESULT=Path(__file__).resolve().parents[1]/"results/anto_v57_physalis_official_dataset_package_20261124.json"
DOI="10.5061/dryad.xpnvx0kvp"
class PhysalisPackageStop(ValueError):pass
def need(ok,why):
    if not ok:raise PhysalisPackageStop(why)

def inspect():
    _,p=original.check_gate()
    record=p["source_study_versions"][0]
    need(record["dataset_doi"]==DOI and record["dryad_dataset_id"]==176194
         and record["published_version"]==5 and len(record["author_files"])==10,
         "V57_ORIGINAL_PHYSALIS_REPO_AUTHORITY_DRIFT")
    key=urllib.parse.quote("doi:"+DOI,safe="")
    url="https://datadryad.org/api/v2/datasets/"+key+"/download"
    req=urllib.request.Request(url,headers={
        "User-Agent":"ANTO-v57-official-original-dataset-ZIP/1.0",
        "Accept":"application/zip,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=95) as resp:
        body=resp.read(1_000_001)
    need(0<len(body)<=1_000_000,"V57_ORIGINAL_DATASET_PACKAGE_TOO_LARGE")
    with zipfile.ZipFile(io.BytesIO(body)) as z:
        need(z.testzip() is None,"V57_ORIGINAL_ZIP_FAILED_INTEGRITY")
        found={}
        for info in z.infolist():
            if info.is_dir():continue
            name=Path(info.filename).name
            need(name and name not in found,"V57_DUPLICATE_OR_EMPTY_ZIP_FILENAMES")
            found[name]=z.read(info)
    wanted={r["file_name"]:r for r in record["author_files"]}
    need(set(found)==set(wanted),"V57_DRYAD_PACKAGE_DOES_NOT_MATCH_SOURCE_FILE_INVENTORY")
    # Official current published version's metadata provides *author digest*
    # for each original file; ZIP package is admitted only if all agree.
    official=original.d.fetch_json("https://datadryad.org/api/v2/datasets/"+
                                  key+"/versions")
    versions=official.get("_embedded",{}).get("stash:versions") or []
    need(len(versions)==1 and versions[0]["versionNumber"]==5,
         "V57_PHYSALIS_ORIGINAL_VERSION_CHANGED")
    href=(versions[0].get("_links",{}).get("stash:files") or {}).get("href")
    need(href.startswith("/api/v2/"),"V57_DRYAD_PACKAGE_FILELIST_MISSING")
    meta=original.d.fetch_json("https://datadryad.org"+href)
    filemeta={z["path"]:z for z in meta.get("_embedded",{}).get("stash:files",[])}
    need(set(filemeta)==set(wanted),"V57_ORIGINAL_METADATA_ZIP_COHORT_DIFFERS")
    file_receipts=[];csv={}
    for name,record in wanted.items():
        blob=found[name]
        need(len(blob)==record["size_bytes"]
             and len(blob)==filemeta[name]["size"],
             "V57_ORIGINAL_PACKAGE_INDIVIDUAL_FILE_BYTES_CHANGED")
        receipt=original.original_digest_check(name,blob,filemeta[name])
        need(receipt["digest_comparison_status"].startswith("provider_")
             and receipt["digest_comparison_status"].endswith("_verified"),
             "V57_ORIGINAL_PACKAGE_AUTHOR_FILE_DIGEST_UNVERIFIED")
        file_receipts.append(receipt)
        if name.lower().endswith(".csv"):
            csv[name]=original.rows_csv(blob)
    abund="Oikos_2020_cross_season_abundance.csv"
    yieldfile="Oikos_2020_total_plant_yield.csv"
    trial="Oikos_2021_fruit_set_poll_obs.csv"
    need(all(k in csv for k in (abund,yieldfile,trial)),
         "V57_PHYSALIS_2020_2021_ORIGINAL_SOURCE_TABLES_MISSING")
    return {
        "original_study":"Davis et al 2026 Oikos one published original Physalis study",
        "original_dryad_dataset_doi":DOI,
        "original_author_dataset_numeric_id":176194,
        "original_author_dataset_version":5,
        "documented_official_firstparty_Dryad_dataset_ZIP_used":True,
        "original_individual_file_api_previous_status":"HTTP401",
        "original_raw_plant_visits_seed_rows_not_exported":True,
        "original_reproduction_or_visitation_values_modelled":False,
        "n_verified_original_author_files":len(file_receipts),
        "original_author_file_digests":file_receipts,
        "original_2020_sameplant_nominal_id_cohort":
            original.cohort_2020(csv[abund],csv[yieldfile]),
        "original_2021_singlevisit_nominal_id_cohort":
            original.cohort_2021(csv[trial])
    }

def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    try:
        out=inspect()
    except Exception as exc:
        out={"original_dryad_dataset_doi":DOI,
             "status":"ORIGINAL_SOURCE_ACCESS_OR_EQUIVALENCE_STOP",
             "error_type":type(exc).__name__,
             "source_error_code":str(exc)[:160],
             "original_reproduction_or_visitation_values_modelled":False,
             "new_independent_ecological_effect_admitted":False}
    if args.check:
        need(json.loads(RESULT.read_text())==out,
             "V57_PHYSALIS_OFFICIAL_DATASET_ZIP_SOURCE_RECEIPT_DRIFT")
        print("V57_ORIGINAL_PHYSALIS_OFFICIAL_ZIP_SOURCE_EXACT_PASS")
    else:
        print("V57_PHYSALIS_OFFICIAL_ZIP_GATE="+json.dumps(out,sort_keys=True,separators=(",",":")))
if __name__=="__main__":main()
