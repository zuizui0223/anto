#!/usr/bin/env python3
"""Source-only discovery of two original independent plant-visitor studies.

No original reproductive, visitation or numerical outcome values read or fit.
Uses only the official Dryad DOI-indexed REST endpoints and metadata-search
on official Zenodo. No guessed source files or unauthorized mirrors.
"""
from __future__ import annotations
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REG=ROOT/"data/external/anto_v57_original_thirdstudy_source_registry_20261124.json"
RESULT=ROOT/"results/anto_v57_thirdstudy_original_metadata_gate_20261124.json"
BASE="https://datadryad.org/api/v2/"
ZENODO="https://zenodo.org/api/records"

class V57Error(ValueError):pass
def need(ok,code):
    if not ok:raise V57Error(code)
def read_registry():
    r=json.loads(REG.read_text())
    need(r["analysis_type"]=="V57_INDEPENDENT_ORIGINAL_SAME_PLANT_SOURCE_ELIGIBILITY_NO_OUTCOME_SCORING"
         and r["physalis"]["original_dataset_doi"]=="10.5061/dryad.xpnvx0kvp"
         and r["mertensia"]["original_dataset_doi"]=="10.7280/D16D7Z"
         and r["minimum_original_sameplant_keys"]==30
         and r["minimum_original_blocks"]==3
         and r["no_ecological_effect_fit"] is True
         and r["no_original_raw_data_rehosted"] is True
         and r["new_independent_original_study_results_admitted"]==0,
         "V57_ORIGINAL_SOURCE_ONLY_SCOPE_CHANGED")
    return r

def allowed(url):
    p=urllib.parse.urlsplit(url)
    return p.scheme=="https" and p.hostname in ("datadryad.org","zenodo.org") and p.username is None and p.password is None

def fetch_json(url,limit=2_000_000):
    need(allowed(url),"V57_UNAUTHORIZED_SOURCE_METADATA_URL")
    r=urllib.request.Request(url,headers={"User-Agent":"ANTO-v57-20261124-original-source-gate/1.0",
                                         "Accept":"application/json"})
    with urllib.request.urlopen(r,timeout=35) as h:
        body=h.read(limit+1)
    need(len(body)<=limit,"V57_SOURCE_METADATA_EXCEEDS_BOUND")
    return json.loads(body)

def status(exc):
    if isinstance(exc,urllib.error.HTTPError):
        return {"status":"HTTP_SOURCE_ACCESS_STOP","http_status":exc.code}
    if isinstance(exc,urllib.error.URLError):
        return {"status":"NETWORK_SOURCE_ACCESS_STOP","error_type":type(exc.reason).__name__}
    return {"status":"SCHEMA_OR_SOURCE_ACCESS_STOP","error_type":type(exc).__name__}

def record_search(doi,title):
    # Do not accept apparent title-only matches as original-author provenance.
    results=[]
    qs=[doi, '"'+title+'"']
    for q in qs:
        url=ZENODO+"/?"+urllib.parse.urlencode({"q":q,"size":8})
        try:
            o=fetch_json(url)
            hits=o.get("hits",{}).get("hits",[])
            for record in hits:
                metadata=record.get("metadata") or {}
                identifiers=[str(y.get("identifier","")).lower() for y in
                             metadata.get("related_identifiers",[]) if isinstance(y,dict)]
                has_source=doi.lower() in identifiers or doi.lower() in str(metadata.get("description","")).lower()
                results.append({"id":record.get("id"),
                    "title":str(metadata.get("title",""))[:150],
                    "explicit_source_doi_in_metadata":has_source,
                    "file_key_inventory":[str(f.get("key",f.get("filename",""))) for f in
                                          (record.get("files") if isinstance(record.get("files"),list) else [])][:30]})
        except Exception as exc:
            results.append({"search":q,**status(exc)})
    # Do not infer a source equivalent without exact DOI.
    return {"search_hits":results[:20],
            "exact_original_Dryad_DOI_matches":[x for x in results if x.get("explicit_source_doi_in_metadata")]}

def inspect_dryad(doi):
    encoded=urllib.parse.quote("doi:"+doi,safe="")
    url=BASE+"datasets/"+encoded
    evidence={"original_dataset_doi":doi,
              "original_source_endpoint":"official_dryad_v2"}
    try:
        dataset=fetch_json(url)
        evidence["dataset_metadata_status"]="ORIGINAL_METADATA_RETRIEVED"
        evidence["dataset_keys"]=sorted(dataset.keys())
        evidence["dataset_id"]=dataset.get("id")
        evidence["published_version_key"]=dataset.get("versionNumber",dataset.get("latestVersion"))
        evidence["dataset_title"]=str(dataset.get("title",""))[:180]
        evidence["dataset_links"]=sorted((dataset.get("_links") or {}).keys())
    except Exception as exc:
        evidence["dataset_metadata_status"]=status(exc)
    try:
        versions=fetch_json(url+"/versions")
        evidence["versions_metadata_status"]="ORIGINAL_VERSION_LIST_RETRIEVED"
        evidence["version_response_keys"]=sorted(versions.keys()) if isinstance(versions,dict) else ["list"]
        potential=versions if isinstance(versions,list) else (
            versions.get("_embedded",{}).get("stash:versions")
            or versions.get("_embedded",{}).get("versions")
            or versions.get("versions") or [])
        evidence["n_versions"]=len(potential)
        evidence["versions_sanitized"]=[{"keys":sorted(x.keys()),
                 "id":x.get("id"),"versionNumber":x.get("versionNumber"),
                 "links":sorted((x.get("_links") or {}).keys())} for x in potential[:12]]
        for x in potential[-3:]:
            ident=x.get("id")
            if not isinstance(ident,int):continue
            try:
                files=fetch_json(BASE+"versions/"+str(ident)+"/files")
                farr=files if isinstance(files,list) else (
                    files.get("_embedded",{}).get("stash:files")
                    or files.get("_embedded",{}).get("files")
                    or files.get("files") or [])
                evidence.setdefault("version_file_inventory",[]).append({
                    "version_id":ident,
                    "file_response_keys":sorted(files.keys()) if isinstance(files,dict) else ["list"],
                    "n_files":len(farr),
                    "original_files":[{"keys":sorted(f.keys()),
                        "id":f.get("id"),"fileName":f.get("fileName",f.get("path")),
                        "size":f.get("size"),"sha256":f.get("sha256"),"md5":f.get("md5"),
                        "links":sorted((f.get("_links") or {}).keys())}
                        for f in farr[:25]]})
            except Exception as e:
                evidence.setdefault("version_file_inventory",[]).append({
                    "version_id":ident,"status":status(e)})
    except Exception as exc:
        evidence["versions_metadata_status"]=status(exc)
    return evidence

def main():
    r=read_registry()
    candidates=[
        ("Physalis 2026 Davis et al",r["physalis"]["original_dataset_doi"],"Plant responses to herbivory shift the contributions of different pollinator taxa to plant reproduction"),
        ("Mertensia 2021 Gallagher Campbell",r["mertensia"]["original_dataset_doi"],"Experimental test of the combined effects of water availability and flowering time on pollinator visitation and seed set")]
    out=[]
    for code,doi,title in candidates:
        dryad=inspect_dryad(doi)
        # Lookup archive only if direct original API cannot supply author file versions.
        record={"candidate":code,"doi":doi,"dryad":dryad,
                "files_opened_or_outcomes_modelled":False}
        if not dryad.get("version_file_inventory"):
            record["zenodo_original_doi_metadata_search"]=record_search(doi,title)
        out.append(record)
    result={"analysis_type":r["analysis_type"],
            "n_original_studies_checked":2,
            "source_discovery_only":True,
            "no_original_numerical_trait_visitor_seed_values_read":True,
            "no_predictive_model_fitted":True,
            "no_new_original_study_outcome_evidence_claimed":True,
            "records":out}
    raw=json.dumps(result,sort_keys=True,separators=(",",":"),ensure_ascii=True)
    print("V57_OFFICIAL_ORIGINAL_SOURCE_METADATA="+raw)

if __name__=="__main__":main()
