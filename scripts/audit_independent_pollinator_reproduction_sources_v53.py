#!/usr/bin/env python3
"""ANTO v53: two independent Zenodo original-source ID/column audits only.

NO seed/fruit outcome column values are numerically inspected. This is a
post-Sinapis eligibility screen, NOT an original model or biological null.
Original author files are NOT committed to GitHub.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit

ROOT=Path(__file__).resolve().parents[1]
REGISTRY=ROOT/"data/external/anto_v53_independent_source_registry_20261119.json"
RECEIPT=ROOT/"results/anto_v53_original_source_grain_audit_20261119.json"
MAX_METADATA=2_000_000

class SourceGateError(ValueError): pass

def require(ok,code):
    if not ok: raise SourceGateError(code)

def registry():
    x=json.loads(REGISTRY.read_text())
    require(x["analysis_type"]=="V53_INDEPENDENT_SOURCE_IDENTITY_AND_GRAIN_ONLY_NO_OUTCOME_READOUT"
            and [a["code"] for a in x["candidates"]]==["Pardee2019","Mota2022"]
            and [a["zenodo_record_id"] for a in x["candidates"]]==[4970401,6671199]
            and [a["original_data_doi"] for a in x["candidates"]]==[
                "10.5061/dryad.kd1r166","10.5061/dryad.h18931zpg"]
            and x["prior_brassica_traine_egwe_access_stop"]["exact_file_stream_status"]==403
            and x["prior_brassica_traine_egwe_access_stop"]["should_retry_unauthenticated"] is False
            and x["max_bytes_per_file"]==2_000_000
            and x["per_plant_same_id_join_min"]==30
            and x["min_independent_heldout_blocks"]==3
            and x["min_linked_plants_per_block"]==5
            and x["only_metadata_columns_counts_identifiers"] is True
            and x["may_read_reproductive_outcome_values"] is False
            and x["fit_any_model"] is False
            and x["this_is_not_biological_replication"] is True,
            "V53_SOURCE_REGISTRY_OR_SCIENTIFIC_SCOPE_DRIFT")
    return x


def md5(data):
    return hashlib.md5(data).hexdigest()


def zenodo_url(url):
    parsed=urlsplit(url or "")
    return (parsed.scheme=="https" and parsed.hostname=="zenodo.org"
            and parsed.username is None and parsed.password is None)


def get_bytes(url,max_size):
    require(zenodo_url(url),"V53_UNTRUSTED_SOURCE_HOST_OR_URL")
    req=urllib.request.Request(url,headers={
        "User-Agent":"ANTO-v53-original-source-grain-20261119/1.0",
        "Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=75) as response:
        raw=response.read(max_size+1)
    require(len(raw)<=max_size,"V53_EXCESSIVE_ORIGINAL_SOURCE_SIZE")
    return raw


def fetch_original(candidate, cap):
    rid=candidate["zenodo_record_id"]
    url="https://zenodo.org/api/records/%d"%rid
    meta=json.loads(get_bytes(url,MAX_METADATA))
    require(str(meta.get("id"))==str(rid),"V53_ZENODO_ORIGINAL_RECORD_ID_MISMATCH")
    files=meta.get("files",[])
    if isinstance(files,dict):files=files.get("entries",files)
    if isinstance(files,dict):
        files=[dict(item,key=name) for name,item in files.items()]
    require(isinstance(files,list),"V53_UNKNOWN_ORIGINAL_RECORD_FILES_SCHEMA")
    found={str(k.get("key",k.get("filename"))):k
           for k in files if isinstance(k,dict)}
    require(set(candidate["files_md5"]).issubset(found),
            "V53_ZENODO_REGISTERED_SOURCE_FILE_MISSING")
    authenticated={}
    for name,expected in candidate["files_md5"].items():
        entry=found[name]
        digest=str(entry.get("checksum","")).lower().split(":")[-1]
        require(not digest or digest==expected,
                "V53_ORIGINAL_ARCHIVE_METADATA_MD5_CHANGED_"+name)
        links=entry.get("links") or {}
        url=links.get("content") or links.get("download") or links.get("self")
        require(isinstance(url,str) and zenodo_url(url),
                "V53_ORIGINAL_ARCHIVE_FILE_URL_UNTRUSTED_"+name)
        raw=get_bytes(url,cap)
        require(md5(raw)==expected,"V53_MD5_ORIGINAL_AUTHOR_BYTES_CHANGED_"+name)
        authenticated[name]=raw
    return authenticated


IDENTITY_HINTS=(
    "plant","indiv","tree","genot","clon","plot","site",
    "field","farm","orchard","block","rep","transect","quadrat",
    "location","species","taxon","year","season","date","day",
    "treatment","sample","tag","record","id"
)

def key_field(name):
    word=re.sub(r"[^a-z0-9]+","_",name.strip().lower())
    return (any(s in word for s in IDENTITY_HINTS)
            and not any(s in word for s in (
                "flower_visited","fruitset","seed_set","fruit_set",
                "visit_rate","visitation_rate","number","count",
                "interaction_id","species_richness","abundance","density")))


def is_plant_field(field):
    x=re.sub(r"[^a-z0-9]+","_",field.lower())
    return x in ("plant","plantid","plant_id","id_plant","plant_code",
                 "individual","individual_id","individualid","tag",
                 "plantnumber","plant_no","plantid_number")


def table_schema(name,raw):
    if not name.lower().endswith(".csv"):return {
        "file":name,"original_file_md5":md5(raw),
        "bytes":len(raw),"parsed_as_csv":False,
        "original_binary_values_parsed":False,
    }
    try:
        s=raw.decode("utf-8-sig")
        start=s[:16000]
        try:
            dialect=csv.Sniffer().sniff(start,delimiters=",;\t")
            delimiter=dialect.delimiter
        except csv.Error:
            delimiter=","
        reader=csv.DictReader(io.StringIO(s),delimiter=delimiter)
        names=reader.fieldnames
        require(names and len(set(names))==len(names),
                "V53_UNKNOWN_OR_DUPLICATE_ORIGINAL_CSV_COLUMNS_"+name)
        total=0
        ids={field:set() for field in names if key_field(field)}
        for row in reader:
            require(None not in row and all(x is not None for x in row.values()),
                    "V53_ORIGINAL_CSV_INVALID_ROW_SHAPE_"+name)
            total+=1
            # Access and hash only original identity/group key fields.
            # No visit count, seed number, reproductive outcome or floral
            # trait magnitude is parsed or used in this stage.
            for field in ids:
                value=(row[field] or "").strip()
                if value:
                    ids[field].add(value)
    except UnicodeDecodeError as exc:
        raise SourceGateError("V53_ORIGINAL_CSV_ENCODING_FAILED_"+name) from exc
    key_counts={field:len(v) for field,v in ids.items()}
    return {
        "file":name,"original_file_md5":md5(raw),
        "bytes":len(raw),"parsed_as_csv":True,
        "delimiter":repr(delimiter),
        "n_original_source_rows":total,
        "original_column_names":names,
        "potential_identity_or_block_field_cardinality":key_counts,
        "candidate_plant_fields":[field for field in ids if is_plant_field(field)],
        "no_numeric_flower_visit_or_fruit_seed_outcomes_examined":True,
        "_id_sets_for_internal_join_only":ids,
    }


def source_audit(candidate,files):
    summary=[table_schema(name,body) for name,body in files.items()]
    # Hash/count individual key intersections inside the process only:
    # never save or transmit any original individual IDs.
    relationships=[]
    candidate_plant_files=[x for x in summary if x.get("candidate_plant_fields")]
    for i,a in enumerate(candidate_plant_files):
        for b in candidate_plant_files[i+1:]:
            for ca in a["candidate_plant_fields"]:
                for cb in b["candidate_plant_fields"]:
                    va=a["_id_sets_for_internal_join_only"][ca]
                    vb=b["_id_sets_for_internal_join_only"][cb]
                    relationships.append({
                        "first_file":a["file"],"first_plant_field":ca,
                        "second_file":b["file"],"second_plant_field":cb,
                        "n_exact_shared_original_plant_id_tokens":len(va&vb),
                        "n_first_unique_original_tokens":len(va),
                        "n_second_unique_original_tokens":len(vb),
                        "is_biologically_same_plant_ID_confirmed":False})
    # No automatic admission: overlapping integers and plot labels can
    # spuriously match across species; source README must establish semantics.
    for x in summary:x.pop("_id_sets_for_internal_join_only",None)
    n_high=sum(x["n_exact_shared_original_plant_id_tokens"]>=30
               for x in relationships)
    return {"candidate":candidate["code"],
        "original_data_doi":candidate["original_data_doi"],
        "zenodo_record_id":candidate["zenodo_record_id"],
        "n_source_files_exact_MD5_verified":len(files),
        "source_files_schema":summary,
        "same_token_plant_key_candidate_comparisons":relationships,
        "n_unverified_raw_token_overlaps_at_least_30":n_high,
        "direct_individual_plant_visitation_to_seed_set_join_ADMITTED":False,
        "heldout_block_and_temporal_semantics_verified":False,
        "source_visit_or_reproductive_values_modeled":False,
        "original_data_used_for_numerical_model":False,
        "source_access_status":"ALL_ORIGINAL_SOURCE_MD5_AUTHENTICATED_HEADERS_AND_KEYS_ONLY"}


def evaluate(x):
    reports=[]
    for c in x["candidates"]:
        try:
            source=fetch_original(c,x["max_bytes_per_file"])
            reports.append(source_audit(c,source))
        except Exception as exc:
            reports.append({
                "candidate":c["code"],"zenodo_record_id":c["zenodo_record_id"],
                "source_access_status":"ACCESS_OR_SOURCE_STOP",
                "error_type":type(exc).__name__,
                "source_error_code":str(exc)[:130],
                "direct_individual_plant_visitation_to_seed_set_join_ADMITTED":False,
                "source_visit_or_reproductive_values_modeled":False})
    require(len(reports)==2,"V53_CANDIDATE_REGISTRY_LOST")
    return {
        "analysis_type":x["analysis_type"],
        "registered_candidates_count":2,
        "previous_traine_access_stop_respected":True,
        "no_original_reproductive_or_visitor_values_parsed_for_effect":True,
        "no_independent_original_fitted_predictive_studies_added":True,
        "no_plant_visitor_reproductive_same_ID_automatically_admitted":True,
        "records":reports,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=evaluate(registry())
    if args.check:
        require(result==json.loads(RECEIPT.read_text()),
                "V53_SOURCE_SCHEMA_RECEIPT_REGENERATION_MISMATCH")
        print("V53_ORIGINAL_SOURCE_IDENTITY_RECEIPT_PASS")
    else:
        s=json.dumps(result,sort_keys=True,separators=(",",":"),ensure_ascii=True)
        chunks=(len(s)+3999)//4000
        for i in range(chunks):
            print("V53_ORIGINAL_SOURCE_GRANULARITY_%03d_OF_%03d="%
                  (i+1,chunks)+s[4000*i:4000*(i+1)])


if __name__=="__main__":
    main()
