#!/usr/bin/env python3
"""Source-only original Mertensia plant-ID, week-array and missingness audit.

The original author Excel was authenticated as identical bytes by
Dryad version-4 SHA256 and Zenodo 4731053 MD5. This script reads only
nominal plant/group labels and whether candidate outcome/trait cells
are PRESENT; it does NOT read their magnitudes, fit models or infer
pollination or seed mechanisms.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from pathlib import Path

import audit_original_Mertensia_Zenodo_equivalence_v57 as orig
import audit_original_Physalis_Mertensia_file_grain_v57 as author

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/anto_v57_mertensia_dual_archive_author_grain_20261124.json"
SOURCE_GIT_BLOB="1071510aa742cb3f17b45c1fd6d8e0cc5412a746"
RESULT=ROOT/"results/anto_v57_mertensia_exp3_original_plant_grain_20261124.json"
XLSX="gallagher&campbell_frontiersinecology&evolution_data.xlsx"
EXPECTED_SHA256="152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d"
CATEGORICAL=("week","array","treat","plantID")
MEASURED=("mVWC","flrCnt","mcorW","mcorL","mNecVol","mNecCon",
          "mVisRate","sBombVis","sNonBombVis","sVis","sCnt","mseed")
ROW_FIELDS=CATEGORICAL+MEASURED
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

class OriginalGrainError(ValueError): pass
def need(ok,msg):
    if not ok: raise OriginalGrainError(msg)

def check_receipt():
    raw=SOURCE.read_bytes()
    need(author.gitblob(raw)==SOURCE_GIT_BLOB,
         "V57_DRYAD_ZENODO_ORIGINAL_AUTHOR_SOURCE_RECEIPT_CHANGED")
    z=json.loads(raw)
    need(z["source_identity"]==
         "exact_original_Dryad_v4_files_via_original_Zenodo_archival_record"
         and z["original_dryad_doi"]=="10.7280/D16D7Z"
         and z["dryad_original_version"]==4
         and z["original_zenodo_archival_record"]==4731053
         and len(z["dual_digest_identical_original_files"])==2,
         "V57_UNVERIFIED_MERTENSIA_AUTHOR_COHORT")
    f=next(x for x in z["dual_digest_identical_original_files"]
           if x["author_file"]==XLSX)
    need(f["sha256"]==EXPECTED_SHA256
         and f["byte_identity_validated_by_both_archives"] is True,
         "V57_MERTENSIA_AUTHOR_WORKBOOK_DIGEST_DRIFT")
    return z

def download_authenticated_book():
    base=orig.get_original_official_metadata()
    repo=orig.verify_zenodo()
    raw=orig.fetch_zenodo_exact(repo[XLSX],XLSX,orig.FILES[XLSX])
    a,dtype=orig.checksum_matches(raw,base[XLSX]["digest"])
    b,zt=orig.checksum_matches(raw,repo[XLSX]["checksum"])
    need(a and b and dtype in ("SHA256","MD5","SHA1") and zt=="MD5"
         and hashlib.sha256(raw).hexdigest()==EXPECTED_SHA256,
         "V57_MERTENSIA_ORIGINAL_AUTHOR_WORKBOOK_AUTH_FAILED")
    return raw

def colnumber(name):
    # Excel row reference, letters A through ZZ.
    ans=0
    for c in name.upper():
        need("A"<=c<="Z","V57_INVALID_XLSX_COLUMN_LABEL")
        ans=26*ans+ord(c)-64
    return ans-1

def text_values(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        need(archive.testzip() is None and
             "xl/worksheets/sheet11.xml" in archive.namelist(),
             "V57_ORIGINAL_EXCEL_MASTER_SHEET_MISSING")
        shared=[]
        if "xl/sharedStrings.xml" in archive.namelist():
            parent=ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for node in parent.findall("m:si",NS):
                shared.append("".join(v.text or "" for v in node.findall(".//m:t",NS)))
        xml=ET.fromstring(archive.read("xl/worksheets/sheet11.xml"))
        rows=xml.findall(".//m:sheetData/m:row",NS)
        need(rows,"V57_ORIGINAL_EXCEL_ANALYSIS_ROWS_EMPTY")
        def only_one_cell(row, wanted):
            for c in row.findall("m:c",NS):
                ref=c.attrib.get("r","")
                match=re.match(r"^([A-Z]+)\d+$",ref)
                if not match: continue
                col=colnumber(match.group(1))
                if col not in wanted:continue
                typ=c.attrib.get("t")
                if typ=="inlineStr":
                    item=c.find("m:is",NS)
                    text="".join(x.text or "" for x in item.findall(".//m:t",NS)) if item is not None else ""
                    yield col,text.strip()
                    continue
                val=c.find("m:v",NS)
                if val is None:continue
                if typ=="s":
                    idx=int(val.text or "0")
                    need(0<=idx<len(shared),"V57_SHARED_STRING_INDEX_INVALID")
                    text=shared[idx]
                elif typ in ("str","e"):
                    text=val.text or ""
                else:
                    # This is allowed only for (week,array,treat,plantID)
                    # categorical source IDs. Numeric measurement cell
                    # contents are NEVER returned or used, only presence.
                    text=val.text or ""
                yield col,text.strip()
        header_map=dict(only_one_cell(rows[0],set(range(80))))
        cols={val:key for key,val in header_map.items() if val}
        need(all(k in cols for k in ROW_FIELDS),
             "V57_MERTENSIA_ORIGINAL_MASTER_REQUIRED_COLUMN_NAMES_CHANGED")
        nominal_idx={cols[k] for k in CATEGORICAL}
        measure_idx={cols[k] for k in MEASURED}
        used_indices=nominal_idx|measure_idx
        def classify(cell):
            if not cell: return False
            if cell in ("NA","NaN","N/A","#VALUE!","#N/A",".","NA/NaN"):return False
            return True
        visited=set()
        duplicate=0
        counted=0
        omitted_nominal=0
        measured_present=Counter()
        measured_missing=Counter()
        arrays=Counter()
        treatments=Counter()
        weeks=Counter()
        cross=Counter()
        plant_keys_in_cluster=defaultdict(set)
        complete=0
        for row in rows[1:]:
            # We retrieve values only for source categorical identifiers;
            # for measured fields we inspect cell *presence*, not number.
            nominal={}
            present=set()
            for c in row.findall("m:c",NS):
                m=re.match(r"^([A-Z]+)\d+$",c.attrib.get("r",""))
                if not m:continue
                col=colnumber(m.group(1))
                if col not in used_indices:continue
                v=c.find("m:v",NS)
                if col in nominal_idx:
                    if c.attrib.get("t")=="s":
                        idx=int(v.text or "0") if v is not None else -1
                        text=shared[idx] if 0<=idx<len(shared) else ""
                    elif c.attrib.get("t")=="inlineStr":
                        text="".join((n.text or "") for n in c.findall(".//m:t",NS))
                    else:
                        text=v.text if v is not None else ""
                    nominal[next(k for k in CATEGORICAL if cols[k]==col)]=str(text or "").strip()
                else:
                    # Do not examine numerical outcome magnitude. Exclude
                    # explicitly missing formula-error/NA string cells.
                    if v is None and c.attrib.get("t")!="inlineStr":
                        continue
                    if c.attrib.get("t")=="s":
                        idx=int(v.text or "0")
                        rawtext=shared[idx] if 0<=idx<len(shared) else ""
                        active=classify(rawtext)
                    elif c.attrib.get("t")=="e":
                        active=False
                    else:
                        active=True
                    if active:
                        present.add(next(k for k in MEASURED if cols[k]==col))
            if not any(nominal.values()):
                continue
            counted+=1
            if any(not nominal.get(k) for k in CATEGORICAL):
                omitted_nominal+=1
                continue
            key=tuple(nominal[k] for k in CATEGORICAL)
            if key in visited:duplicate+=1
            visited.add(key)
            week,array,treat,plant=key
            wkarr=(week,array)
            arrays[wkarr]+=1
            weeks[week]+=1
            treatments[treat]+=1
            cross[(week,treat)]+=1
            plant_keys_in_cluster[wkarr].add(plant)
            for name in MEASURED:
                (measured_present if name in present else measured_missing)[name]+=1
            complete+=all(x in present for x in
                          ("flrCnt","mcorW","mcorL","mVisRate","mseed"))
        need(counted>0 and len(visited)>0,
             "V57_MERTENSIA_ORIGINAL_PLANT_LABELS_NOT_FOUND")
        groups={("|".join(k)):n for k,n in sorted(arrays.items())}
        eligible=[k for k,n in arrays.items() if n>=5]
        readiness=(len(visited)>=30 and len(eligible)>=3
                   and complete>=30 and duplicate==0 and omitted_nominal==0)
        return {
            "original_workbook_sheet":"exp3_meansForAnalysis",
            "source_sheet_xml_path":"xl/worksheets/sheet11.xml",
            "source_original_header_names":
                [header_map[k] for k in sorted(header_map)],
            "original_nonempty_plant_record_rows":counted,
            "original_distinct_composite_week_array_treat_plantID_keys":len(visited),
            "source_row_duplicate_composite_plantkeys":duplicate,
            "source_rows_missing_any_categorical_identifier":omitted_nominal,
            "nominal_week_categories":dict(sorted(weeks.items())),
            "nominal_original_treatment_categories":dict(sorted(treatments.items())),
            "nominal_week_treatment_record_counts":[
                {"week":w,"treat":t,"n":n} for (w,t),n in sorted(cross.items())],
            "original_week_by_array_clusters":groups,
            "n_original_week_array_clusters":len(arrays),
            "n_clusters_with_at_least_5_source_plant_records":len(eligible),
            "n_original_plant_records_all_candidate_variables_PRESENT":complete,
            "candidate_var_present_count":dict(sorted(measured_present.items())),
            "candidate_var_missing_count":dict(sorted(measured_missing.items())),
            "source_only_prerequisite_three_original_group_cohort_ready":readiness,
            "source_original_measurement_magnitudes_examined":False,
            "source_original_visitation_seed_association_fitted":False,
            "original_plant_ID_tokens_or_outcomes_exported":False}

def examine():
    evidence=check_receipt()
    book=download_authenticated_book()
    result=text_values(book)
    return {
        "analysis_type":"V57_ORIGINAL_MERTENSIA_EXPERIMENT3_SAMEPLANT_ORIGINAL_GROUP_GRAIN_ONLY",
        "original_dryad_doi":"10.7280/D16D7Z",
        "original_dryad_version":4,
        "original_zenodo_source_record":4731053,
        "original_workbook_sha256":EXPECTED_SHA256,
        "source_original_author_workbook_authentic":True,
        "new_observed_fitness_effect_fitted":False,
        "evolutionary_pollinator_mechanism_inferred":False,
        "original_source_group_grain":result}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    try:
        obj=examine()
    except Exception as exc:
        obj={"analysis_type":
             "V57_ORIGINAL_MERTENSIA_EXPERIMENT3_SAMEPLANT_ORIGINAL_GROUP_GRAIN_ONLY",
             "source_access_or_ID_schema_stop":True,
             "error_type":type(exc).__name__,
             "error_code":str(exc)[:160],
             "source_original_visitation_seed_association_fitted":False}
    if args.check:
        need(obj==json.loads(RESULT.read_text()),
             "V57_ORIGINAL_MERTENSIA_PLANT_GRAIN_RECEIPT_DRIFT")
        print("V57_MERTENSIA_ORIGINAL_2016_EXPERIMENT3_GROUP_SOURCE_ID_PASS")
    else:
        raw=json.dumps(obj,sort_keys=True,separators=(",",":"),
                       ensure_ascii=True)
        n=(len(raw)+3999)//4000
        for i in range(n):
            print("V57_MERTENSIA_SOURCE_GROUP_%03d_OF_%03d="%
                  (i+1,n)+raw[4000*i:4000*(i+1)])
if __name__=="__main__":main()
