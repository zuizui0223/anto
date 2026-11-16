"""Source-only original Dryad small CSV fetch; never download 1.74-GB archive.

Official Dryad DOI 10.5061/dryad.zgmsbccr3; metadata, author version
and file-ID / digest validation. No original raw rows committed to GitHub.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import quote, urlsplit

DOI="10.5061/dryad.zgmsbccr3"
ROOT="https://datadryad.org"
META=ROOT+"/api/v2/datasets/"+quote("doi:"+DOI,safe="")
NAMES=(
    "vigs_fruitset.csv","vigs_pollentubes.csv","vigs_qpcr.csv",
    "diallelcross_fruitdata.csv","diallelcross_pollentubedata.csv",
)
MAX_FILE_BYTES=150_000
MAX_METADATA_BYTES=300_000

class SourceGateError(ValueError):
    pass


def admitted_api_url(url):
    p=urlsplit(url)
    if (p.scheme!="https" or p.hostname!="datadryad.org" or p.port is not None
            or p.username or p.password or p.query or p.fragment
            or not p.path.startswith("/api/v2/")
            or ".." in p.path):
        raise SourceGateError("SOURCE_URL_NOT_OFFICIAL_DRYAD_API")
    return url


def fetch_bytes(url, limit):
    admitted_api_url(url)
    req=urllib.request.Request(
        url,headers={"User-Agent":"ANTO-v39-source-authenticated-Phlox-SI/1.0",
                     "Accept":"application/json,text/csv,*/*;q=0.6"})
    with urllib.request.urlopen(req,timeout=40) as resp:
        effective=resp.geturl()
        # Dryad's official file endpoint may redirect to object storage; only
        # source API URL, original metadata IDs and digest authenticate bytes.
        if url.endswith("/files") and urlsplit(effective).hostname!="datadryad.org":
            raise SourceGateError("METADATA_UNTRUSTED_REDIRECT")
        body=resp.read(limit+1)
    if not body or len(body)>limit:
        raise SourceGateError("ORIGINAL_SOURCE_RESPONSE_SIZE_UNADMITTED")
    return body


def retry(url,limit,attempts=3):
    last=None
    for i in range(attempts):
        try:
            return fetch_bytes(url,limit)
        except urllib.error.HTTPError as e:
            if e.code in (401,403):
                raise SourceGateError("SOURCE_ACCESS_RESTRICTED_HTTP_"+str(e.code)) from e
            if e.code not in (429,500,502,503,504):
                raise
            last=e
        except (urllib.error.URLError,TimeoutError,ConnectionError) as e:
            last=e
        if i<attempts-1:
            time.sleep(2)
    raise SourceGateError("SOURCE_TEMPORARY_TRANSPORT_FAILURE_"+type(last).__name__)


def official_version_ref(meta):
    if meta.get("identifier")!="doi:"+DOI:
        raise SourceGateError("DRYAD_DATASET_DOI_IDENTITY_DRIFT")
    path=meta.get("_links",{}).get("stash:version",{}).get("href")
    if not isinstance(path,str) or not re.fullmatch(r"/api/v2/versions/\d+",path):
        raise SourceGateError("DRYAD_OFFICIAL_VERSION_LINK_UNADMITTED")
    return ROOT+path+"/files"


def original_file_inventory(version_listing):
    files=version_listing.get("_embedded",{}).get("stash:files",[])
    total=version_listing.get("total",len(files))
    if not isinstance(files,list) or not files or total>len(files):
        raise SourceGateError("DRYAD_AUTHOR_FILE_PAGE_INCOMPLETE")
    names={}
    for item in files:
        filename=Path(str(item.get("path",""))).name
        if filename not in NAMES:
            continue
        if filename in names:
            raise SourceGateError("DRYAD_DUPLICATE_REGISTERED_FILENAME")
        source_link=item.get("_links",{}).get("stash:download",{}).get("href")
        if not isinstance(source_link,str) or not re.fullmatch(r"/api/v2/files/\d+/download",source_link):
            raise SourceGateError("DRYAD_FILE_ORIGINAL_ID_UNADMITTED")
        raw_size=item.get("size")
        if isinstance(raw_size,bool) or not isinstance(raw_size,int) or not 50<=raw_size<=MAX_FILE_BYTES:
            raise SourceGateError("DRYAD_FILE_DECLARED_SIZE_UNADMITTED")
        kind=str(item.get("digestType","")).lower()
        digest=item.get("digest")
        if kind not in ("md5","sha256","sha-256","") or (digest and not re.fullmatch(r"[a-fA-F0-9]{32}|[a-fA-F0-9]{64}",str(digest))):
            raise SourceGateError("DRYAD_DIGEST_MANIFEST_UNRECOGNIZED")
        names[filename]={
            "file_name":filename,"original_path":item["path"],
            "download_url":ROOT+source_link,
            "official_file_id":int(source_link.split("/")[4]),
            "declared_bytes":raw_size,
            "original_manifest_digest_type":kind,
            "original_manifest_digest":digest,
        }
    if set(names)!=set(NAMES):
        raise SourceGateError("DRYAD_SOURCE_CSV_SET_INCOMPLETE")
    return names


def admit_original_file(raw,meta):
    if not isinstance(raw,bytes) or len(raw)!=meta["declared_bytes"]:
        raise SourceGateError("DRYAD_ORIGINAL_FILE_SIZE_DRIFT")
    declared=meta["original_manifest_digest"]
    kind=meta["original_manifest_digest_type"]
    if not declared:
        raise SourceGateError("NO_ORIGINAL_PROVIDER_DIGEST_FOR_AUTHENTICATION")
    if kind=="md5":
        computed=hashlib.md5(raw).hexdigest()
    elif kind in ("sha256","sha-256"):
        computed=hashlib.sha256(raw).hexdigest()
    else:
        raise SourceGateError("NO_ADMISSIBLE_ORIGINAL_DIGEST_TYPE")
    if computed.lower()!=declared.lower():
        raise SourceGateError("DRYAD_OFFICIAL_FILE_DIGEST_MISMATCH")
    return {
        "source_file_id":meta["official_file_id"],
        "original_filename":meta["file_name"],
        "declared_bytes":meta["declared_bytes"],
        "real_bytes_verified":len(raw),
        "original_provider_digest_type":kind,
        "original_provider_digest":declared,
        "sha256_computed_from_original_bytes":hashlib.sha256(raw).hexdigest(),
        "original_provider_digest_matched":True,
    }


def run(cache):
    doc=json.loads(retry(META,MAX_METADATA_BYTES))
    listing=json.loads(retry(official_version_ref(doc),MAX_METADATA_BYTES))
    candidates=original_file_inventory(listing)
    cache.mkdir(parents=True,exist_ok=True)
    receipt={}
    for key in NAMES:
        meta=candidates[key]
        raw=retry(meta["download_url"],MAX_FILE_BYTES)
        receipt[key]=admit_original_file(raw,meta)
        (cache/key).write_bytes(raw)
    return {
        "schema":"anto_v39_original_dryad_Phlox_PdPIRK_five_CSV_source_integrity",
        "dataset_DOI":DOI,"official_metadata_url":META,
        "source_version_from_original_Dryad_metadata":True,
        "all_five_original_files_source_digest_verified":True,
        "downloaded_total_original_bytes":sum(r["real_bytes_verified"] for r in receipt.values()),
        "no_large_transcriptome_or_images_downloaded":True,
        "source_files":receipt,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--status",type=Path,required=True)
    opts=p.parse_args()
    try:
        response=run(opts.cache)
    except Exception as exc:
        # Technical inability is explicitly not a negative functional result.
        response={
            "schema":"anto_v39_original_Dryad_source_access_gate",
            "dataset_DOI":DOI,
            "source_files_admitted":False,
            "technical_source_status":"NOT_ADMITTED",
            "error_type":type(exc).__name__,
            "reason":str(exc) if isinstance(exc,(SourceGateError,urllib.error.HTTPError)) else type(exc).__name__,
            "original_biological_result_recomputed":False,
            "source_substitution_or_permission_bypass":False,
        }
        opts.status.write_text(json.dumps(response,indent=2,sort_keys=True)+"\n")
        print("ANTO_V39_SOURCE_ACCESS_BLOCK",json.dumps(response,sort_keys=True))
        raise SystemExit(2)
    opts.status.write_text(json.dumps(response,indent=2,sort_keys=True)+"\n")
    print("ANTO_V39_AUTHENTICATED_ORIGINAL_DRYAD_SOURCE",json.dumps(response,sort_keys=True))


if __name__=="__main__":
    main()
