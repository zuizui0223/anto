#!/usr/bin/env python3
"""Recover exact original 2021 Mertensia workbook using dual archive checksums.

Source-only: no plant responses/pollinator counts read or modeled.
The official Dryad per-file endpoints returned HTTP401 and remain stopped.
Zenodo 4731053 is admissible ONLY if byte-for-byte original Dryad digest
and its own original archived MD5 match with exact size and source titles.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import discover_original_Physalis_Mertensia_v57 as d
import audit_original_Physalis_Mertensia_file_grain_v57 as author

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v57_verified_dryad_file_ids_20261124.json"
RECEIPT=ROOT/"results/anto_v57_mertensia_dual_archive_author_grain_20261124.json"
ZENODO_RECORD=4731053
ORIGINAL_DOI="10.7280/D16D7Z"
ORIGINAL_TITLE=("experimental test of the combined effects of water availability "
                "and flowering time on pollinator visitation and seed set")
FILES={
 "gallagher&campbell_frontiersinecology&evolution_data.xlsx":339816,
 "gallagher&campbell_frontiersinecology&evolution_metaDataFile.docx":47265
}

class SourceIdentityError(ValueError): pass
def need(cond,label):
    if not cond:raise SourceIdentityError(label)

def get_original_official_metadata():
    r,p=author.check_gate()
    expected=p["source_study_versions"][1]
    need(expected["dataset_doi"]==ORIGINAL_DOI
         and expected["published_version"]==4
         and expected["dryad_dataset_id"]==59594
         and {x["file_name"]:x["size_bytes"] for x in expected["author_files"]}==FILES,
         "V57_MERTENSIA_ORIGINAL_DRYAD_PRIOR_VERSION_DRIFT")
    key=urllib.parse.quote("doi:"+ORIGINAL_DOI,safe="")
    head="https://datadryad.org/api/v2/datasets/"+key
    doc=d.fetch_json(head)
    need(doc["id"]==59594
         and doc.get("versionNumber")==4
         and str(doc.get("title","")).strip().lower()==ORIGINAL_TITLE,
         "V57_MERTENSIA_CURRENT_DRYAD_DOI_TITLE_VERSION_CHANGED")
    versions=d.fetch_json(head+"/versions")
    items=versions.get("_embedded",{}).get("stash:versions") or []
    need(len(items)==1 and items[0]["versionNumber"]==4,
         "V57_MERTENSIA_ORIGINAL_VERSION_NOT_UNIQUE")
    files_link=(items[0].get("_links",{}).get("stash:files") or {}).get("href")
    need(files_link and files_link.startswith("/api/v2/"),
         "V57_MERTENSIA_SOURCE_FILE_LINK_CHANGED")
    resp=d.fetch_json("https://datadryad.org"+files_link)
    files=resp.get("_embedded",{}).get("stash:files") or []
    source={z["path"]:z for z in files}
    need(set(source)==set(FILES),
         "V57_MERTENSIA_PUBLISHED_AUTHOR_FILES_DIFFER")
    for f,size in FILES.items():
        need(source[f].get("size")==size
             and source[f].get("_links",{}).get("stash:download",{}).get("href")==
             next(x["official_download_href"] for x in expected["author_files"]
                  if x["file_name"]==f),
             "V57_MERTENSIA_AUTHOR_SOURCE_NAME_SIZE_OR_DRYAD_FILE_DRIFT")
    return source

def checksum_matches(data,hashstr):
    source=(hashstr or "").lower().strip().split(":")[-1]
    if len(source)==32:
        return source==hashlib.md5(data).hexdigest(),"MD5"
    if len(source)==40:
        return source==hashlib.sha1(data).hexdigest(),"SHA1"
    if len(source)==64:
        return source==hashlib.sha256(data).hexdigest(),"SHA256"
    return False,"ORIGINAL_DIGEST_UNRECOGNIZED"

def verify_zenodo():
    r=d.fetch_json("https://zenodo.org/api/records/"+str(ZENODO_RECORD))
    need(r.get("id")==ZENODO_RECORD,
         "V57_WRONG_OFFICIAL_ZENODO_RECORD")
    meta=r.get("metadata") or {}
    title=" ".join(str(meta.get("title","")).strip().lower().split())
    need(title==ORIGINAL_TITLE,
         "V57_MERTENSIA_ZENODO_TITLE_IS_NOT_ORIGINAL_PAPER")
    files=r.get("files") or []
    if isinstance(files,dict):
        files=[dict(item,key=name) for name,item in files.items()]
    need(isinstance(files,list),
         "V57_MERTENSIA_ZENODO_FILE_LIST_MISSING")
    known={str(z.get("key",z.get("filename",""))):z for z in files}
    need(set(known)==set(FILES),
         "V57_MERTENSIA_ZENODO_ARCHIVED_FILE_NAMES_CHANGED")
    for f,size in FILES.items():
        need(known[f]["size"]==size,
             "V57_MERTENSIA_ZENODO_AUTHOR_FILE_SIZE_CHANGED")
    return known

def fetch_zenodo_exact(filemeta,filename,size):
    link=(filemeta.get("links") or {}).get("content")
    need(isinstance(link,str),"V57_ZENODO_OFFICIAL_CONTENT_URL_MISSING")
    u=urllib.parse.urlsplit(link)
    need(u.scheme=="https" and u.hostname=="zenodo.org"
         and u.username is None and u.password is None,
         "V57_ZENODO_UNTRUSTED_FILE_STREAM")
    req=urllib.request.Request(link,headers={
        "User-Agent":"ANTO-v57-original-archive-digest-verify/1.0"})
    with urllib.request.urlopen(req,timeout=90) as resp:
        raw=resp.read(size+1)
    need(len(raw)==size,
         "V57_ZENODO_ARCHIVED_SOURCE_BYTE_COUNT_CHANGED")
    return raw

def assess_original_source():
    source=get_original_official_metadata()
    archived=verify_zenodo()
    verified=[]
    headers=None
    for fname,size in FILES.items():
        original=source[fname]
        mirror=archived[fname]
        b=fetch_zenodo_exact(mirror,fname,size)
        ok_zenodo,checksumtype=checksum_matches(b,mirror.get("checksum"))
        need(ok_zenodo and checksumtype=="MD5",
             "V57_MERTENSIA_ZENODO_ORIGINAL_MD5_FAILED_"+fname)
        ok_dryad,dryadtype=checksum_matches(b,original.get("digest"))
        need(ok_dryad and dryadtype in ("MD5","SHA1","SHA256"),
             "V57_MERTENSIA_AUTHOR_DRYAD_DIGEST_NOT_EQUAL_"+fname)
        verified.append({
            "author_file":fname,"bytes":len(b),
            "sha256":hashlib.sha256(b).hexdigest(),
            "zenodo_md5":hashlib.md5(b).hexdigest(),
            "original_dryad_digest_type":original.get("digestType"),
            "byte_identity_validated_by_both_archives":True})
        if fname.endswith(".xlsx"):
            headers=author.inspect_xlsx(b)
        elif fname.endswith(".docx"):
            with zipfile.ZipFile(io.BytesIO(b)) as z:
                need(z.testzip() is None
                     and "word/document.xml" in z.namelist(),
                     "V57_ORIGINAL_DOCX_INVALID")
    need(len(verified)==2 and headers,
         "V57_MERTENSIA_COMPLETE_ORIGINAL_SOURCE_NOT_ADMITTED")
    return {"source_identity":"exact_original_Dryad_v4_files_via_original_Zenodo_archival_record",
            "original_dryad_doi":ORIGINAL_DOI,"dryad_id":59594,
            "dryad_original_version":4,
            "original_zenodo_archival_record":4731053,
            "official_dryad_individual_download_401_remains_a_stop":True,
            "dual_digest_identical_original_files":verified,
            "original_workbook_header_grain":headers,
            "sameplant_visit_seed_outcome_model_fitted":False,
            "original_individual_seed_visit_numeric_values_inspected":False,
            "new_independent_original_study_effect_confirmed":False,
            "sameplant_outcome_model_ready_without_future_cohort_check":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    try:
        result=assess_original_source()
    except Exception as err:
        result={"source_identity":"SOURCE_EQUIVALENCE_STOP",
                "original_dryad_doi":ORIGINAL_DOI,
                "original_zenodo_archival_record":ZENODO_RECORD,
                "error_type":type(err).__name__,
                "source_error_code":str(err)[:140],
                "sameplant_visit_seed_outcome_model_fitted":False,
                "new_independent_original_study_effect_confirmed":False}
    if args.check:
        need(result==json.loads(RECEIPT.read_text()),
             "V57_ORIGINAL_DOUBLE_ARCHIVE_RECEIPT_REGENERATION_FAILED")
        print("V57_MERTENSIA_ORIGINAL_DOUBLE_ARCHIVE_SOURCE_GATE_EXACT_PASS")
    else:
        raw=json.dumps(result,sort_keys=True,separators=(",",":"),ensure_ascii=True)
        n=(len(raw)+3999)//4000
        for i in range(n):
            print("V57_MERTENSIA_DOUBLE_ARCHIVE_%03d_OF_%03d="%
                 (i+1,n)+raw[4000*i:4000*(i+1)])
if __name__=="__main__":main()
