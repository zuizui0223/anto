#!/usr/bin/env python3
"""ANTO v57 original-author file-byte, cohort identity and group gate.

Fetch ONLY canonical original Dryad file URLs from a previously frozen
official v2 DOI/version/file ID inventory. Never fit predictor or examine
seed/fruit/visit numeric outcomes; never commit original author source rows.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ORIG=ROOT/"data/external/anto_v57_original_thirdstudy_source_registry_20261124.json"
PIN=ROOT/"data/external/anto_v57_verified_dryad_file_ids_20261124.json"
RESULT=ROOT/"results/anto_v57_author_source_identity_grain_20261124.json"
PIN_GIT_BLOB="66ea1296a3a5f3957d20f777dd49faf815f4902c"
DRYAD="https://datadryad.org"
MAX_FILE=2_000_000

class V57FileError(ValueError): pass
def need(test,label):
    if not test:raise V57FileError(label)
def gitblob(raw):
    return hashlib.sha1(("blob %d\0"%len(raw)).encode()+raw).hexdigest()
def check_gate():
    raw=PIN.read_bytes()
    need(gitblob(raw)==PIN_GIT_BLOB,"V57_FROZEN_ORIGINAL_DRYAD_SOURCE_FILE_INVENTORY_DRIFT")
    pin=json.loads(raw)
    reg=json.loads(ORIG.read_text())
    need(reg["no_ecological_effect_fit"] is True
         and reg["new_independent_original_study_results_admitted"]==0
         and len(pin["source_study_versions"])==2
         and pin["source_study_versions"][0]["published_version"]==5
         and pin["source_study_versions"][1]["published_version"]==4
         and pin["source_study_versions"][0]["dryad_dataset_id"]==176194
         and pin["source_study_versions"][1]["dryad_dataset_id"]==59594
         and {z["dataset_doi"] for z in pin["source_study_versions"]}=={
            "10.5061/dryad.xpnvx0kvp","10.7280/D16D7Z"}
         and [len(z["author_files"]) for z in pin["source_study_versions"]]==[10,2],
         "V57_ORIGINAL_DOI_VERSION_FILE_COUNTS_CHANGED")
    return reg,pin

def exact_official_download(href):
    need(re.match(r"^/api/v2/files/[0-9]+/download$",href or "") is not None,
         "V57_FILE_STREAM_NOT_IN_ORIGINAL_DRYAD_API")
    return DRYAD+href

def bytes_from_official(href,expected_size):
    url=exact_official_download(href)
    need(0<expected_size<=MAX_FILE,"V57_SOURCE_FILE_SIZE_NOT_ADMITTED")
    req=urllib.request.Request(url,headers={
        "User-Agent":"ANTO-v57-authentic-original-20261124/1.0",
        "Accept":"application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=80) as response:
        payload=response.read(expected_size+1)
    need(len(payload)==expected_size,"V57_ORIGINAL_DRYAD_BYTE_COUNT_CHANGED")
    return payload

def original_digest_check(filename,data,source_record):
    actual=hashlib.sha256(data).hexdigest()
    dg=str(source_record.get("digest","") or "").strip().lower()
    typ=str(source_record.get("digestType","") or "").strip().lower().replace("-","")
    digeststatus="original_sha256_computed_no_comparable_provider_digest"
    if len(dg)==64 and all(c in "0123456789abcdef" for c in dg):
        need(actual==dg,"V57_ORIGINAL_SHA256_MISMATCH")
        digeststatus="provider_sha256_verified"
    elif len(dg)==32 and all(c in "0123456789abcdef" for c in dg):
        need(hashlib.md5(data).hexdigest()==dg,"V57_ORIGINAL_MD5_MISMATCH")
        digeststatus="provider_md5_verified"
    elif dg:
        digeststatus="original_digest_unrecognized_logged_not_verified"
    return {"file":filename,"bytes":len(data),"sha256":actual,
            "original_provider_digest_type":typ,
            "original_provider_digest_value":dg,
            "digest_comparison_status":digeststatus}

def rows_csv(data):
    try:
        raw=data.decode("utf-8-sig")
        sniff=csv.Sniffer()
        delimiter=sniff.sniff(raw[:8000],delimiters=",;\t").delimiter
        reader=csv.DictReader(io.StringIO(raw),delimiter=delimiter)
        fields=reader.fieldnames or []
        need(len(fields)>0 and len(set(fields))==len(fields),
             "V57_CSV_NO_VALID_AUTHOR_HEADERS")
        n=0
        keys={}
        wanted=("block","treat","treatment","id","plant_id")
        keynames=[x for x in fields if x.strip().lower() in wanted]
        unique={x:set() for x in keynames}
        identity=[]
        n_missing_key=0
        for row in reader:
            need(None not in row and all(v is not None for v in row.values()),
                 "V57_BAD_ORIGINAL_CSV_RECORD_SHAPE")
            n+=1
            values=[]
            for k in keynames:
                token=str(row[k]).strip()
                if token and token.lower() not in ("na","nan","null"):
                    unique[k].add(token)
                values.append(token)
            # Collect ONLY nominal identifiers, never numerical outcomes.
            identity.append({k:str(row[k]).strip() for k in keynames})
            n_missing_key+=any(not x for x in values)
        return {"source_headers":fields,"rows":n,
                "nominal_key_cardinality":{k:len(v) for k,v in unique.items()},
                "identity_missing_rows":n_missing_key,
                "delimiter":delimiter,
                "_only_identifier_rows":identity}
    except (UnicodeDecodeError,csv.Error) as exc:
        raise V57FileError("V57_SOURCE_CSV_UNDECODABLE_OR_BAD_DIALECT") from exc

def identity_key(r,columns):
    return tuple(str(r.get(k,"")).strip() for k in columns)

def cohort_2020(cross,fruit):
    def columns(table,poss):
        names=set(table["source_headers"])
        need(set(poss).issubset(names),"V57_PHYSALIS_2020_REQUIRED_PLANT_ID_HEADERS_ABSENT")
        return poss
    want=("block","treat","id")
    columns(cross,want);columns(fruit,want)
    a=[identity_key(r,want) for r in cross["_only_identifier_rows"]]
    b=[identity_key(r,want) for r in fruit["_only_identifier_rows"]]
    need(all(all(x for x in k) for k in a+b),
         "V57_PHYSALIS_2020_ORIGINAL_PLANT_KEY_BLANK")
    ac=Counter(a);bc=Counter(b)
    joined=set(ac)&set(bc)
    by_block=Counter(k[0] for k in joined)
    by_treatment=Counter(k[1] for k in joined)
    return {"original_2020_crossseason_records":len(a),
            "original_2020_fruit_yield_records":len(b),
            "original_plant_composite_key":"block×treat×id",
            "original_visitation_unique_plantkeys":len(ac),
            "original_yield_unique_plantkeys":len(bc),
            "original_same_plant_key_overlap":len(joined),
            "n_source_keys_repeated_in_visits":sum(q>1 for q in ac.values()),
            "n_source_keys_repeated_in_yield":sum(q>1 for q in bc.values()),
            "linked_plantkeys_by_block":dict(sorted(by_block.items())),
            "linked_plantkeys_by_treatment":dict(sorted(by_treatment.items())),
            "source_identity_is_experimentally_validated_by_readme":False,
            "linked_plant_reproduction_or_visit_values_analysed":False}

def cohort_2021(table):
    required=("fruit_set","plant_id","block","treatment","morphospecies")
    need(set(required).issubset(table["source_headers"]),
         "V57_PHYSALIS_2021_SINGLE_VISIT_SOURCE_FIELDS_CHANGED")
    keys=["plant_id","block","treatment"]
    rows=[identity_key(r,keys) for r in table["_only_identifier_rows"]]
    need(all(all(v for v in q) for q in rows),
         "V57_PHYSALIS_2021_BLANK_SOURCE_PLANT_KEY")
    n=Counter(rows)
    by_block=Counter(q[1] for q in n)
    return {"original_single_visit_rows":table["rows"],
            "n_distinct_original_plants_by_block_treatment_id":len(n),
            "n_plants_by_original_block":dict(sorted(by_block.items())),
            "repeated_original_plant_visit_records":sum(z-1 for z in n.values()),
            "direct_original_fruit_set_field_exists":True,
            "source_reproduction_values_analysed":False,
            "original_2021_independent_new_publications":0}

def inspect_xlsx(data):
    names={}
    ns={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        paths=set(z.namelist())
        need("xl/workbook.xml" in paths and z.testzip() is None,
             "V57_AUTHOR_MERTENSIA_XLSX_NOT_A_VALID_ZIP")
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=[n.attrib.get("name") for n in wb.findall(".//m:sheet",ns)]
        shared=[]
        if "xl/sharedStrings.xml" in paths:
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for node in root.findall("m:si",ns):
                shared.append("".join((t.text or "") for t in node.findall(".//m:t",ns)))
        outputs=[]
        for f in sorted(p for p in paths if re.match(r"xl/worksheets/sheet\d+\.xml$",p)):
            first_text=[]
            seen_rows=0
            # Read only the FIRST row strings to identify source labels;
            # numerical values in all sheets are NEVER collected/scored.
            elem=ET.fromstring(z.read(f))
            for row in elem.findall(".//m:sheetData/m:row",ns):
                seen_rows+=1
                if seen_rows>2:break
                if seen_rows!=1:continue
                for cell in row.findall("m:c",ns):
                    val=cell.find("m:v",ns)
                    if val is None:continue
                    if cell.attrib.get("t")=="s":
                        idx=int(val.text or "0")
                        if 0<=idx<len(shared):
                            first_text.append(shared[idx][:150])
                    elif cell.attrib.get("t") in ("inlineStr","str"):
                        txt=val.text or ""
                        if txt:first_text.append(txt[:150])
                if len(first_text)>80:first_text=first_text[:80]
            outputs.append({"worksheet_path":f,
                            "first_row_text_labels_only":first_text})
        need(sheets and outputs,"V57_MERTENSIA_XLSX_NO_VALID_SHEETS")
        return {"original_workbook_sheets":sheets,"n_sheets":len(sheets),
                "first_row_only_header_candidates":outputs,
                "source_seed_or_visitation_numeric_cell_values_parsed":False}

def inspect_author_files(candidate,inventory):
    original=urllib.parse.quote("doi:"+candidate["dataset_doi"],safe="")
    version_url="https://datadryad.org/api/v2/datasets/"+original+"/versions"
    # Cross-check source file list from LIVE official metadata.
    import discover_original_Physalis_Mertensia_v57 as d
    versions=d.fetch_json(version_url)
    items=versions.get("_embedded",{}).get("stash:versions") or []
    need(len(items)==1 and items[0].get("versionNumber")==candidate["published_version"],
         "V57_CURRENT_ORIGINAL_DRYAD_VERSION_DRIFT")
    href=(items[0].get("_links",{}).get("stash:files") or {}).get("href")
    need(isinstance(href,str) and href.startswith("/api/v2/"),
         "V57_ORIGINAL_DRYAD_FILES_METADATA_LINK_DRIFT")
    resp=d.fetch_json("https://datadryad.org"+href)
    farr=resp.get("_embedded",{}).get("stash:files") or []
    live={z.get("path"):z for z in farr}
    expected={a["file_name"]:a for a in candidate["author_files"]}
    need(set(live)==set(expected),"V57_ORIGINAL_FILE_NAMES_CHANGED")
    details=[];tables={};statuses=[]
    for filename,source in expected.items():
        meta=live[filename]
        b=source["size_bytes"]
        link=(meta.get("_links",{}).get("stash:download") or {}).get("href")
        need(int(meta.get("size"))==b and link==source["official_download_href"],
             "V57_ORIGINAL_FILE_LINK_OR_SIZE_CHANGED")
        try:
            raw=bytes_from_official(link,b)
            meta_info=original_digest_check(filename,raw,meta)
            details.append(meta_info)
            if filename.lower().endswith(".csv"):
                table=rows_csv(raw)
                tables[filename]=table
            elif filename.lower().endswith(".xlsx"):
                tables[filename]=inspect_xlsx(raw)
            elif filename.lower().endswith(".docx"):
                with zipfile.ZipFile(io.BytesIO(raw)) as zipf:
                    need(zipf.testzip() is None,
                         "V57_ORIGINAL_DOCX_FAILED_CONTAINER_INTEGRITY")
            elif filename.endswith("README.md"):
                # The original readme is source authorization, not a predicted outcome.
                pass
        except Exception as err:
            statuses.append({"file":filename,"error_type":type(err).__name__,
                             "status":("HTTP_"+str(err.code) if isinstance(err,urllib.error.HTTPError)
                                       else "ORIGINAL_FILE_SOURCE_STOP")})
            # No claim that a study is eligible unless all required files
            # and source identities are admitted.
    out={"candidate":candidate["candidate"],
         "dataset_doi":candidate["dataset_doi"],
         "n_expected_original_files":len(expected),
         "n_downloaded_verified_original_files":len(details),
         "files":details,
         "download_errors":statuses,
         "author_original_raw_observations_exported":False,
         "any_predictive_model_fitted":False}
    if candidate["dataset_doi"]=="10.5061/dryad.xpnvx0kvp":
        need(len(details)==10 or statuses,
             "V57_PHYSALIS_ORIGINAL_FILE_SET_INCOMPLETE")
        if (a:="Oikos_2020_cross_season_abundance.csv") in tables and (
            b:="Oikos_2020_total_plant_yield.csv") in tables:
            out["source_2020_sameplant_cohort"]=cohort_2020(tables[a],tables[b])
        if (c:="Oikos_2021_fruit_set_poll_obs.csv") in tables:
            out["source_2021_singlevisit_cohort"]=cohort_2021(tables[c])
    else:
        a="gallagher&campbell_frontiersinecology&evolution_data.xlsx"
        if a in tables:out["source_mertensia_original_xlsx_headers"]=tables[a]
    out["original_source_data_admitted"]=(len(details)==len(expected)
                                          and not statuses)
    return out

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    reg,pin=check_gate()
    out=[]
    for candidate in pin["source_study_versions"]:
        try:
            out.append(inspect_author_files(candidate,pin))
        except Exception as err:
            out.append({"candidate":candidate["candidate"],
                        "dataset_doi":candidate["dataset_doi"],
                        "status":"ORIGINAL_SOURCE_ACCESS_OR_SCHEMA_STOP",
                        "error_type":type(err).__name__,
                        "original_source_data_admitted":False,
                        "original_outcome_models_fitted":False})
    result={"analysis_type":reg["analysis_type"],
            "frozen_official_Dryad_file_inventory_Git_blob":PIN_GIT_BLOB,
            "two_original_published_studies_checked":2,
            "author_source_grain_not_fitted_outcome":True,
            "all_source_numeric_outcomes_held_unanalysed":True,
            "source_records":out}
    if args.check:
        need(json.loads(RESULT.read_text())==result,
             "V57_SOURCE_GRAIN_RECEIPT_CHANGED")
        print("V57_ORIGINAL_SOURCES_IDENTITY_RECEIPT_EXACT_PASS")
    else:
        dump=json.dumps(result,sort_keys=True,separators=(",",":"),ensure_ascii=True)
        n=(len(dump)+3999)//4000
        for i in range(n):
            print("V57_AUTHOR_BYTE_GRAIN_%03d_OF_%03d="%(i+1,n)+dump[4000*i:4000*(i+1)])
if __name__=="__main__":main()
