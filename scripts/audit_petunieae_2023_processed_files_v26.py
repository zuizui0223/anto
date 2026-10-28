"""Fetch and MD5/SHA256-verify five original 2023 Petunieae *processed* files.

All source identities (original OSF project zg9cu, version IDs, byte counts
and digests) were read from official OSF JSONAPI source-file metadata in
ANTO v26 source-discovery CI. Source: Wheeler et al. (2023), Proc. R.
Soc. B 290, 20230275, DOI 10.1098/rspb.2023.0275.

The original source bytes are retained only on the local runner. Do NOT
republish original data in ANTO or call this a genome-functional, causal,
regain-history or holdout analysis. A source-verified processed table
is not the raw RNA-seq and does not independently validate the authors'
quantification from raw sequencing/HPLC chromatograms.

A body with the wrong original file length, digest, or unusable CSV
is rejected. Network/authorization/format errors are *technical source
nonadmissions*, not empty gene expression or biological negative results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import urllib.error
import urllib.request
from pathlib import Path

PAPER_DOI="10.1098/rspb.2023.0275"
OSF_PROJECT="zg9cu"
FILES={
 "expression_pigment_matched":{
  "file_id":"63ae244c202f170475a6c560",
  "source_path":"/phyloCCA/tpm10k-mgg-combined-with-flavs.csv",
  "size":26917,
  "md5":"13a786471791dc1175849d66b2d4a2ef",
  "sha256":"5843d4cd4eb253046f97349fa6bd285ca77e43e7a9c3aaa0e78fae3e8e391edd",
 },
 "expression_gene_mean":{
  "file_id":"63ae217d3312f10280687d53",
  "source_path":"/salmon_expression_quantification/target_gene_expression_mean_handle_multiplicity-tpm10k.csv",
  "size":30410,
  "md5":"9f2ea379ecd7a8e2a56e739b8e56a3bd",
  "sha256":"6f62a4e385d7c9d12e63fd07f6d33a2c79e469af63223f01cec6e8268e68f5c5",
 },
 "pigment_mean":{
  "file_id":"63af0bf70b719802e57c3a08",
  "source_path":"/dataset_summary_statistics/pigment_distributions/combined-antho-flav-mgg-mean-4letter-codes.csv",
  "size":7039,
  "md5":"407400bde406e266456ab0b6589f2707",
  "sha256":"74b4d20a383d478266b67bb7b31c20dd4ec7cc949f2ffb612b19fd3a0e7b2d73",
 },
 "species_tree_pruned":{
  "file_id":"63ae2453c8da860249df343a",
  "source_path":"/phyloCCA/11genestre_dated_pruned.tre",
  "size":1481,
  "md5":"f270feef2909dec19a67eefe15efcb2f",
  "sha256":"95b4a688d3d71417b712b37a2b04cdc22a9431be3172d6509def4435f5fd8614",
 },
 "pigment_cluster_source":{
  "file_id":"63af2c55202f170510a6cfdf",
  "source_path":"/stochastic_mapping/phyloPCA-pigment-kmeansclust-IDs.csv",
  "size":551,
  "md5":"83be371a2cc535948886e8d38f95db9b",
  "sha256":"f30ea096c66a0cdbad159d42edaa818144b13c6d8d63553cd1da51c4331f46ca",
 },
}


def registered_location(record):
    fid=record["file_id"]
    if not (isinstance(fid,str) and len(fid)==24 and all(c in "0123456789abcdef" for c in fid)):
        raise ValueError("OSF original file id not a 24-character hex id")
    return f"https://osf.io/download/{fid}/"


def verify_original(data,record):
    if type(data) is not bytes:
        raise ValueError("Source body must be bytes")
    if len(data)!=record["size"]:
        raise ValueError("Original OSF source byte length mismatch")
    if hashlib.md5(data).hexdigest()!=record["md5"]:
        raise ValueError("Original OSF source MD5 mismatch")
    if hashlib.sha256(data).hexdigest()!=record["sha256"]:
        raise ValueError("Original OSF source SHA256 mismatch")
    return True


def get_official_file(url,max_bytes=1_000_000,timeout=25):
    if not url.startswith("https://osf.io/download/") or not url.endswith("/"):
        raise ValueError("Only source-registered OSF file downloads permitted")
    req=urllib.request.Request(url,headers={
        "User-Agent":"ANTO-original-Petunieae-2023-verified-source-reader/1.0",
        "Accept":"text/csv,text/plain,*/*;q=0.5"})
    with urllib.request.urlopen(req,timeout=timeout) as response:
        value=response.read(max_bytes+1)
    if len(value)>max_bytes or not value:
        raise ValueError("Unexpected original OSF response size")
    return value


def source_summary_from_bytes(label,original):
    """Read *verified* contents but preserve data row/column identifiers only."""
    if label=="species_tree_pruned":
        content=original.decode("utf-8-sig")
        if not content.strip().endswith(";") or "(" not in content:
            raise ValueError("Original source tree Newick structure absent")
        return {"type":"NEWICK_TREE","characters":len(content),
                "ends_in_semicolon":True,"source_tip_labels_not_yet_reconciled":True}
    text=original.decode("utf-8-sig")
    reader=csv.reader(io.StringIO(text))
    header=next(reader,None)
    if not header or len(set(header))!=len(header):
        raise ValueError("Original processed CSV header missing or duplicated")
    rows=list(reader)
    if not rows:
        raise ValueError("Original processed CSV has no data rows")
    if any(len(row)!=len(header) for row in rows):
        raise ValueError("Original processed CSV nonrectangular")
    if len(rows)>10000 or len(header)>600:
        raise ValueError("Unexpected original source table geometry")
    nonempty=[sum(bool(s.strip()) for s in row) for row in rows]
    firstcol=[z[0] for z in rows]
    return {"type":"CSV",
            "source_row_count":len(rows),
            "source_column_count":len(header),
            "original_column_names":header,
            "first_column_name":header[0],
            "distinct_first_column_values":len(set(firstcol)),
            "nonempty_cells":sum(nonempty),
            "empty_cells":len(header)*len(rows)-sum(nonempty),
            "first_6_source_first_column_values":firstcol[:6],
            "first_3_original_row_cells_per_column_first_12":[z[:12] for z in rows[:3]],
            "original_individual_to_species_grain_NOT_YET_ADMITTED":True,
            "no_new_heldout_or_causal_prediction_calculated":True,
           }


def preflight(fetch=get_official_file):
    """Bytes are not saved, unless later separately consented original analysis."""
    r={
      "schema":"anto_v26_original_2023_petunieae_processed_source_body_preflight",
      "audit_date":"2026-10-28",
      "original_paper_DOI":PAPER_DOI,
      "original_OSF_project":OSF_PROJECT,
      "source_rows_are_authors_processed_not_independent_RNA_HPLC_requantification":True,
      "original_processed_file_count_expected":len(FILES),
      "source_file_statuses":{},
      "complete_joint_expression_pigment_source_verified":False,
      "independent_species_holdout_analysis_PERFORMED":False,
      "full_functional_pathway_or_ancestral_regain_measured":False,
      "natural_selection_heritable_return_or_fitness_estimated":False,
      "original_2019_Iochrominae_RAR_source_restored":False,
      "Cirsium_empirical_owner":"zuizui0223/aza3",
    }
    for label,record in FILES.items():
        entry={"osf_id":record["file_id"],
               "source_materialized_path":record["source_path"],
               "source_metadata_registered_md5":record["md5"],
               "source_metadata_registered_sha256":record["sha256"],
               "source_declared_bytes":record["size"],
               "real_source_file_body_hash_verified":False}
        try:
            data=fetch(registered_location(record))
            verify_original(data,record)
            entry["real_source_file_body_hash_verified"]=True
            entry["source_structure"]=source_summary_from_bytes(label,data)
            entry["status"]="ORIGINAL_PROCESSED_BYTES_SHA256_AND_MD5_VERIFIED"
        except Exception as exc:
            entry["status"]=("SOURCE_DOWNLOAD_HTTP_"+str(exc.code) if isinstance(exc,urllib.error.HTTPError)
                             else "SOURCE_TECHNICAL_"+type(exc).__name__)
        r["source_file_statuses"][label]=entry
    r["complete_joint_expression_pigment_source_verified"]=all(
        r["source_file_statuses"][label]["real_source_file_body_hash_verified"]
        for label in ("expression_pigment_matched","pigment_mean","expression_gene_mean","species_tree_pruned"))
    return r


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    out=preflight()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ANTO_V26_PETUNIEAE_PROCESSED_SOURCE_BODY_PREFLIGHT",
          json.dumps(out,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":
    main()
