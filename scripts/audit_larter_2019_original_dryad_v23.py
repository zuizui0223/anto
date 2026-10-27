"""Source-exact 2019 Iochrominae white-morph versus fixed-white Dryad gate.

The study's original authors, Larter, Dunbar-Wallis, Berardi & Smith
(2019), Developmental Dynamics, DOI 10.1002/dvdy.82, deposited at Dryad
doi:10.5061/dryad.p5dq84v. This is an attempt to recover the original
published RAR archive without falsely claiming observed new biological
re-expression or reconstructed ancestral state. This study must NOT be
merged with Wheeler et al. 2025's 24 blue/red hybrid RNA specimens.

Download via official Dryad endpoints only. Admit only if bytes match
the registered original provider MD5 AND declared length AND RAR magic.
An HTTP 401/403, HTML bot challenge, ZIP lacking the true RAR, incorrect
digest, or incomplete download is a non-admission, NOT an observed
biological negative outcome.

This source access gate does not claim an independent reanalysis. Original
RAR bytes must remain ephemeral on the Actions runner; no raw data rehost.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

DATA_DOI="10.5061/dryad.p5dq84v"
PUB_DOI="10.1002/dvdy.82"
REGISTERED_FILE_ID=108456
README_ID=108458
REGISTERED_BYTES=3054999
REGISTERED_MD5="76b46e384fb7c9bf1ef3fbd7d1e5d2f0"
REGISTERED_DATASET_VERSION_ID=28856
RAR4_MAGIC=b"Rar!\x1a\x07\x00"
RAR5_MAGIC=b"Rar!\x1a\x07\x01\x00"
DOI_TOKEN=urllib.parse.quote("doi:"+DATA_DOI,safe="")
API_BASE="https://datadryad.org/api/v2/"
DATASET_API=API_BASE+"datasets/"+DOI_TOKEN
FILE_API=API_BASE+f"files/{REGISTERED_FILE_ID}/download"
VERSION_API=API_BASE+f"versions/{REGISTERED_DATASET_VERSION_ID}/files"
STREAM_URL=f"https://datadryad.org/downloads/file_stream/{REGISTERED_FILE_ID}"
BUNDLE_URL=DATASET_API+"/download"
MAX_REPLY=6_000_000


def original_bytes_valid(data:bytes,md5=REGISTERED_MD5,length=REGISTERED_BYTES):
    """A deterministic strict gate; never infer raw data from HTML or JSON."""
    if type(data) is not bytes or len(data)!=length:
        return False
    if not (data.startswith(RAR4_MAGIC) or data.startswith(RAR5_MAGIC)):
        return False
    return hashlib.md5(data).hexdigest()==md5


def _safe_status(exc):
    if isinstance(exc,urllib.error.HTTPError):
        return f"HTTP_{exc.code}"
    if isinstance(exc,TimeoutError):
        return "TIMEOUT"
    if isinstance(exc,urllib.error.URLError):
        return "NETWORK_UNAVAILABLE"
    return "INVALID_SOURCE_OR_UNEXPECTED_RESPONSE"


def fetch_limited(url,timeout=18,limit=MAX_REPLY):
    if not isinstance(url,str) or not (
            url.startswith("https://datadryad.org/api/v2/") or
            url.startswith("https://datadryad.org/downloads/file_stream/")):
        raise ValueError("Only official Dryad public source paths allowed")
    req=urllib.request.Request(url,headers={
        "Accept":"application/octet-stream,application/x-rar-compressed,"
                 "application/zip,application/json;q=0.9,*/*;q=0.5",
        "User-Agent":"ANTO-2019-Original-Dryad-Source-Integrity-Audit/1.0",
    })
    with urllib.request.urlopen(req,timeout=timeout) as resp:
        if getattr(resp,"status",200)!=200:
            raise ValueError("Original source did not return status 200")
        data=resp.read(limit+1)
    if not data or len(data)>limit:
        raise ValueError("Original source empty or exceeds safe response cap")
    return data


def registered_api_details(fetch=fetch_limited):
    """Check official file metadata WITHOUT taking its statements as an assay."""
    out={"dataset_metadata_retrieved":False,
         "version_manifest_retrieved":False,
         "original_rar_metadata_match":False,
         "file_metadata_status":"NOT_CHECKED"}
    try:
        dataset=json.loads(fetch(DATASET_API).decode("utf-8"))
        if dataset.get("identifier")!="doi:"+DATA_DOI:
            raise ValueError("Dryad DOI identity mismatch")
        out["dataset_metadata_retrieved"]=True
        out["dryad_dataset_id"]=dataset.get("id")
        v=dataset.get("_links",{}).get("stash:version",{}).get("href")
        if v!=f"/api/v2/versions/{REGISTERED_DATASET_VERSION_ID}":
            raise ValueError("Dryad original version pointer drift")
        manifest=json.loads(fetch(VERSION_API).decode("utf-8"))
        rows=manifest.get("_embedded",{}).get("stash:files",[])
        matches=[r for r in rows if (str(r.get("id"))==str(REGISTERED_FILE_ID)
                 or r.get("_links",{}).get("self",{}).get("href")==
                     f"/api/v2/files/{REGISTERED_FILE_ID}")]
        if len(matches)!=1:
            raise ValueError("Original Dryad file ID missing or ambiguous")
        f=matches[0]
        if (f.get("size")!=REGISTERED_BYTES or
                f.get("digest","").lower()!=REGISTERED_MD5 or
                f.get("digestType","").lower()!="md5" or
                not f.get("path","").lower().endswith(".rar")):
            raise ValueError("Dryad original RAR source digest/size/name drift")
        out.update({"version_manifest_retrieved":True,
                    "original_rar_metadata_match":True,
                    "file_metadata_status":"VERIFIED_OFFICIAL_MANIFEST",
                    "file_name":Path(f["path"]).name})
    except Exception as exc:
        out["file_metadata_status"]=_safe_status(exc)
        out["metadata_failure_type"]=type(exc).__name__
    return out


def _find_exact_rar_in_bundle(raw:bytes):
    """Search original provider zip for exactly the MD5-verified RAR, no rehost."""
    if len(raw)>MAX_REPLY:
        return None
    try:
        archive=zipfile.ZipFile(io.BytesIO(raw))
        members=archive.infolist()
        if len(members)>100:
            return None
        for info in members:
            if info.is_dir() or info.file_size!=REGISTERED_BYTES:
                continue
            name=Path(info.filename)
            if name.is_absolute() or ".." in name.parts:
                continue
            data=archive.read(info)
            if original_bytes_valid(data):
                return data
    except (zipfile.BadZipFile,RuntimeError,ValueError):
        return None
    return None


def attempt_original_source(fetch=fetch_limited):
    """Live official Dryad public acquisition, honest status per candidate."""
    candidates=[
        ("file_API",FILE_API),
        ("file_stream",STREAM_URL),
        ("original_dataset_ZIP",BUNDLE_URL),
    ]
    attempts=[]
    admitted=None
    retrieval=None
    for label,url in candidates:
        try:
            data=fetch(url)
            passed=original_bytes_valid(data)
            if not passed and label=="original_dataset_ZIP":
                extracted=_find_exact_rar_in_bundle(data)
                if extracted is not None:
                    data=extracted
                    passed=True
            attempts.append({
                "route":label,
                "retrieval_status":"ADMITTED_ORIGINAL_MD5_RAR" if passed
                    else "REJECTED_NOT_ORIGINAL_RAR",
                "returned_bytes":len(data),
                "original_md5_matched":passed,
            })
            if passed:
                admitted=data
                retrieval=label
                break
        except Exception as exc:
            attempts.append({
                "route":label,
                "retrieval_status":_safe_status(exc),
                "original_md5_matched":False,
            })
    return admitted,retrieval,attempts


def receipt(metadata,acquired,route,attempts):
    """No phenotype, enzyme, expression or inheritance result is created."""
    accepted=acquired is not None and original_bytes_valid(acquired)
    if acquired is not None and not accepted:
        raise ValueError("Attempt to admit unverified original dataset")
    return {
        "schema":"anto_v23_larter2019_white_morph_fixed_white_original_access_gate",
        "audit_date":"2026-10-28",
        "original_pub_DOI":PUB_DOI,
        "original_public_dataset_DOI":DATA_DOI,
        "original_Dryad_file_ID":REGISTERED_FILE_ID,
        "registered_source_bytes":REGISTERED_BYTES,
        "registered_source_MD5":REGISTERED_MD5,
        "source_metadata":metadata,
        "source_acquisition_attempts":attempts,
        "original_RAR_bytes_verified":accepted,
        "successful_route":route if accepted else None,
        "biological_source_admitted_for_future_raw_reanalysis":accepted,
        "biological_2019_expression_reanalysed_in_this_gate":False,
        "raw_gene_expression_sample_rows_parsed":0,
        "white_morph_versus_fixed_white_new_statistical_result":False,
        "source_assay_functionally_validated":False,
        "trait_regain_or_multigeneration_retention_estimated":False,
        "paper_observational_design":"Larter2019 white morphs and fixed white taxa",
        "not_Wheeler2025_red_blue_backcross":True,
        "source_file_raw_bytes_rehosted_in_git":False,
        "technical_blocker_not_biological_negative":not accepted,
        "Cirsium_empirical_owner":"zuizui0223/aza3",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--original-rar",type=Path)
    parser.add_argument("--from-local-rar",type=Path,
                        help="Optional independently supplied original RAR; always check original MD5")
    parser.add_argument("--expected",type=Path,
                        help="Fail closed if the frozen source acquisition gate has changed")
    args=parser.parse_args()
    metadata=registered_api_details()
    if args.from_local_rar is not None:
        local=args.from_local_rar
        if not local.is_file() or local.stat().st_size!=REGISTERED_BYTES:
            raise ValueError("Locally supplied archive does not have original registered size")
        data=local.read_bytes()
        if not original_bytes_valid(data):
            raise ValueError("Locally supplied archive failed original MD5 or RAR signature")
        route="locally_supplied_exact_registered_RAR"
        attempts=[{"route":route,"retrieval_status":"ADMITTED_ORIGINAL_MD5_RAR",
                   "returned_bytes":len(data),"original_md5_matched":True}]
    else:
        data,route,attempts=attempt_original_source()
    result=receipt(metadata,data,route,attempts)
    if args.expected is not None:
        expected=json.loads(args.expected.read_text(encoding="utf-8"))
        # A new successful retrieval, source version change, or previously
        # inaccessible endpoint becoming accessible triggers explicit review.
        keys=("original_public_dataset_DOI","registered_source_bytes",
              "registered_source_MD5","source_metadata",
              "source_acquisition_attempts","original_RAR_bytes_verified",
              "biological_source_admitted_for_future_raw_reanalysis",
              "technical_blocker_not_biological_negative")
        if any(result.get(k)!=expected.get(k) for k in keys):
            raise RuntimeError(
                "DRYAD_2019_GATE_TRANSITION_REQUIRES_MANUAL_AUDIT")
    if data is not None and args.original_rar:
        args.original_rar.parent.mkdir(parents=True,exist_ok=True)
        args.original_rar.write_bytes(data)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("ANTO_V23_ORIGINAL_SOURCE_GATE",json.dumps({
        "source_admitted":result["original_RAR_bytes_verified"],
        "metadata":result["source_metadata"]["file_metadata_status"],
        "attempts":result["source_acquisition_attempts"],
        "2019_expression_reanalysed":False
    },sort_keys=True))


if __name__=="__main__":
    main()
