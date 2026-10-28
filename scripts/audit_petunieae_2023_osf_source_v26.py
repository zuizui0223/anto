"""Source-only original processed-data discovery for Wheeler et al. 2023 Petunieae.

Author-published study: Wheeler, Dunbar-Wallis, Schutz & Smith 2023,
Proceedings of the Royal Society B, DOI 10.1098/rspb.2023.0275.
OSF project publicly declared in the original article: https://osf.io/zg9cu/
SRA raw RNA BioProject PRJNA746328, not read by this inventory.

Only requests official public OSF *metadata*, not original gene data, and
never invents measured species, expression, gene function or retention.
Not a biological replication of Iochrominae: different Petunieae clade,
but both in Solanaceae and same laboratory; not independent 59 regain events.

If source API is unavailable, emit an honest technical nonadmission.
No real processed-data regression without original downloaded,
checksum-verified, specimen-identified, response and predictor matched data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PROJECT="zg9cu"
PAPER_DOI="10.1098/rspb.2023.0275"
SRA_BIOPROJECT="PRJNA746328"
API="https://api.osf.io/v2"
NODE=f"{API}/nodes/{PROJECT}/"
FILES=f"{NODE}files/osfstorage/"
MAX_PAGES=80
MAX_FILES=2000
MAX_DEPTH=6
MAX_BYTES=3_000_000
EXTENSIONS=(".csv",".tsv",".xlsx",".txt",".json",".rds",".rdata",".rda",".nwk",".tree",".tre")


def normalize_api_url(url):
    """Never follow non-public or untrusted API locators for traversal."""
    if not isinstance(url,str) or not url.startswith(API+"/"):
        raise ValueError("Source traversal must stay within official OSF v2 API")
    parsed=urllib.parse.urlsplit(url)
    if parsed.scheme!="https" or parsed.netloc!="api.osf.io" or parsed.username or parsed.password:
        raise ValueError("Source traversal URL not official public OSF")
    return url


def get_bytes(url, *, timeout=20):
    url=normalize_api_url(url)
    req=urllib.request.Request(url,headers={"Accept":"application/vnd.api+json",
        "User-Agent":"ANTO-original-Petunieae-2023-dataset-source-audit/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as res:
        value=res.read(MAX_BYTES+1)
    if len(value)>MAX_BYTES or not value:
        raise ValueError("Unexpected OSF API metadata response size or empty")
    return value


def error_code(exc):
    if isinstance(exc,urllib.error.HTTPError):
        return f"HTTP_{exc.code}"
    if isinstance(exc,(urllib.error.URLError,TimeoutError)):
        return "NETWORK_UNAVAILABLE"
    return f"METADATA_ERROR_{type(exc).__name__}"


def _safe_filename(x):
    text=x if isinstance(x,str) else ""
    # Keep only a bounded display string; original filenames are source
    # metadata, not commands, paths to write, or biological observations.
    return text[:250].replace("\n"," ").replace("\r"," ")


def inspect_source(fetch=get_bytes):
    """Discover bounded official public file metadata, including nested folders."""
    out={
       "schema":"anto_v26_wheeler2023_original_petunieae_osf_qualified_source_inventory",
       "audit_date":"2026-10-28",
       "original_paper_DOI":PAPER_DOI,
       "original_2023_OSF_project":PROJECT,
       "original_author_declared_SRA_raw_RNA":SRA_BIOPROJECT,
       "source_API_root":FILES,
       "node_metadata_available":False,
       "folder_inventory_complete":False,
       "source_manifest_admitted":False,
       "original_species_level_joint_expression_and_pigment_rows_parsed":0,
       "original_2023_processed_expression_or_HPLC_bytes_downloaded":False,
       "new_provenance_matched_plant_functional_assays":0,
       "independent_evolutionary_regain_events_newly_measured":0,
       "original_2019_Iochrominae_white_morph_RAR_reanalyzed":False,
       "source_only_not_cross_clade_predictive_validation":True,
       "same_Solanaceae_family_not_independent_clade_loss_replication":True,
       "Cirsium_empirical_owner":"zuizui0223/aza3",
       "status":"NOT_CHECKED",
       "source_file_records":[],
       "tabular_candidate_files":[],
       "warnings":[],
    }
    try:
        metadata=json.loads(fetch(NODE))
        data=metadata.get("data",{})
        if data.get("id")!=PROJECT:
            raise ValueError("Expected original 2023 OSF project node identity")
        attrs=data.get("attributes",{})
        out["node_metadata_available"]=True
        out["source_node_title"]=_safe_filename(attrs.get("title"))
        out["source_node_public"]=bool(attrs.get("public",False))
        if not out["source_node_public"]:
            raise ValueError("Original 2023 OSF project is not publicly readable")
    except Exception as exc:
        out["status"]="SOURCE_NODE_METADATA_"+error_code(exc)
        out["warnings"].append(out["status"])
        return out

    queue=[(FILES,0)]
    seen_urls=set()
    seen_filenames=set()
    pages=0
    try:
        while queue:
            url,depth=queue.pop(0)
            if depth>MAX_DEPTH:
                raise ValueError("Maximum registered source-folder nesting exceeded")
            url=normalize_api_url(url)
            if url in seen_urls:
                continue
            seen_urls.add(url)
            pages+=1
            if pages>MAX_PAGES:
                raise ValueError("OSF public source pagination exceeded")
            response=json.loads(fetch(url))
            records=response.get("data")
            if not isinstance(records,list):
                raise ValueError("OSF source file response lacks JSONAPI list")
            for item in records:
                attrs=item.get("attributes",{})
                kind=attrs.get("kind")
                name=_safe_filename(attrs.get("name"))
                fid=item.get("id")
                size=attrs.get("size")
                if not name or not fid or kind not in ("folder","file"):
                    raise ValueError("OSF source file identity/kind missing")
                if kind=="folder":
                    related=item.get("relationships",{}).get("files",{}).get("links",{}).get("related",{})
                    child=(related.get("href") if isinstance(related,dict) else related)
                    if child is None:
                        child=item.get("links",{}).get("move")
                    if not child:
                        raise ValueError("OSF source folder has no children locator")
                    queue.append((normalize_api_url(child),depth+1))
                    continue
                full_path=_safe_filename(attrs.get("materialized_path")) or name
                if full_path in seen_filenames:
                    raise ValueError("Duplicated original source file materialized path")
                seen_filenames.add(full_path)
                if len(seen_filenames)>MAX_FILES:
                    raise ValueError("More source files than bounded audit supports")
                item_info={
                    "osf_file_id":str(fid),
                    "source_name":name,
                    "source_materialized_path":full_path,
                    "source_declared_bytes":size if isinstance(size,int) else None,
                    "source_declared_hashes":attrs.get("extra",{}).get("hashes"),
                    "source_original_file_download_url":item.get("links",{}).get("download"),
                    "file_body_retrieved":False,
                }
                out["source_file_records"].append(item_info)
                low=name.lower()
                if low.endswith(EXTENSIONS):
                    out["tabular_candidate_files"].append(item_info)
            nxt=response.get("links",{}).get("next")
            if nxt:
                queue.append((normalize_api_url(nxt),depth))
        out["source_manifest_admitted"]=True
        out["folder_inventory_complete"]=True
        out["status"]="PUBLIC_OSF_MANIFEST_ADMITTED_METADATA_ONLY"
        out["source_file_records"].sort(key=lambda z:z["source_materialized_path"])
        out["tabular_candidate_files"].sort(key=lambda z:z["source_materialized_path"])
    except Exception as exc:
        out["status"]="SOURCE_TREE_"+error_code(exc)
        out["warnings"].append(out["status"])
        out["source_manifest_admitted"]=False
        out["folder_inventory_complete"]=False
        # Partial file list is useful for debugging, but not an exhaustive
        # biological source panel and must never become an assay denominator.
    out["source_file_count_metadata_only"]=len(out["source_file_records"])
    out["small_tabular_candidate_count_metadata_only"]=len(out["tabular_candidate_files"])
    out["no_function_gene_ancestry_or_regain_claims_from_metadata"]=True
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=inspect_source()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    subset={
        "status":result["status"],
        "node_metadata_available":result["node_metadata_available"],
        "folder_inventory_complete":result["folder_inventory_complete"],
        "source_file_count":result.get("source_file_count_metadata_only",0),
        "tabular_candidates":result.get("small_tabular_candidate_count_metadata_only",0),
        "first_80_tabular_paths":[{"path":x["source_materialized_path"],
                                   "size":x["source_declared_bytes"],
                                   "file_id":x["osf_file_id"]}
                                  for x in result["tabular_candidate_files"][:80]],
        "warnings":result["warnings"],
    }
    print("ANTO_V26_PETUNIEAE_OSF_SOURCE_INVENTORY",json.dumps(subset,sort_keys=True))
    compact=json.dumps(result,ensure_ascii=False,sort_keys=True,separators=(',',':'))
    for j in range(0,len(compact),2500):
        print("ANTO_V26_SOURCE_RECEIPT_CHUNK",j//2500,compact[j:j+2500])


if __name__=="__main__":
    main()
