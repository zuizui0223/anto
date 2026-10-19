"""Source-only read of original published Dryad Iochrominae coding alignments.

Dataset: Ho & Smith (2016), 10.5061/dryad.49rs8.
This is an *independent replay of externally collected data inventory*, NOT a
new regain mechanism, a reanalysis of omega/dN/dS, verified enzyme function,
historical regain rates, or an independent plant species sample from ANTO.

Only dataset metadata and aggregate file identities are written. The source
sequences remain with Dryad; this repository does not copy the raw alignments.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request
import urllib.error
import zipfile

DOI="10.5061/dryad.49rs8"
DOI_ENC=urllib.parse.quote("doi:"+DOI,safe="")
METADATA_URL=f"https://datadryad.org/api/v2/datasets/{DOI_ENC}"
DOWNLOAD_URL=METADATA_URL+"/download"
EXPECTED_FILES=("EVOL_CHI_SEQFILE.txt","EVOL_DFR_SEQFILE.txt","F3H_SEQFILE.txt")
MAX_BYTES=5_000_000

def sha256(b):
    return hashlib.sha256(b).hexdigest()

def get_bytes(url, timeout=40, limit=MAX_BYTES):
    if not url.startswith("https://datadryad.org/api/v2/"):
        raise ValueError("Only the official, registered public Dryad API is permitted")
    headers={"Accept":"application/zip,application/json;q=0.9,*/*;q=0.5",
             "User-Agent":"ANTO-source-qualification/1.0"}
    req=urllib.request.Request(url,headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as res:
        result=res.read(limit+1)
    if len(result)>limit or not result:
        raise ValueError("Unexpected source download size or empty response")
    return result

def read_alignment_header(b):
    """Classify but DO NOT silently parse/translate unfamiliar alignments."""
    text=b.decode("utf-8-sig",errors="replace")
    first=next((line.strip() for line in text.splitlines() if line.strip()),"")
    if first.startswith(">"):
        return {"format":"FASTA","header":first[:160], "declared_sequences":None,
                "declared_alignment_sites":None}
    if first.upper().startswith("#NEXUS"):
        return {"format":"NEXUS","header":first[:160],"declared_sequences":None,
                "declared_alignment_sites":None}
    m=re.match(r"^(\d+)\s+(\d+)(?:\s+.*)?$",first)
    if m:
        return {"format":"PHYLIP_OR_PAML_HEADER","header":first[:160],
                "declared_sequences":int(m.group(1)),
                "declared_alignment_sites":int(m.group(2))}
    return {"format":"UNCLASSIFIED","header":first[:160],
            "declared_sequences":None,"declared_alignment_sites":None}

def inspect_zip(archive):
    if len(archive)>MAX_BYTES:
        raise ValueError("Oversized source bundle")
    z=zipfile.ZipFile(io.BytesIO(archive))
    members={Path(item.filename).name:item for item in z.infolist()
             if not item.is_dir() and not item.filename.startswith("__MACOSX")}
    absent=set(EXPECTED_FILES)-set(members)
    if absent:
        raise ValueError(f"Original Dryad 3-gene alignment files missing: {sorted(absent)}")
    output={}
    for key in EXPECTED_FILES:
        item=members[key]
        if item.file_size>MAX_BYTES:
            raise ValueError(f"Unexpected decompressed alignment size for {key}")
        b=z.read(item)
        if not (1000<=len(b)<=MAX_BYTES):
            raise ValueError(f"Bad alignment size for {key}")
        info=read_alignment_header(b)
        info.update({"filename":key,"bytes":len(b),"sha256":sha256(b)})
        output[key]=info
    return {"file_count":len(members),"all_file_names":sorted(members),
            "target_gene_file_metadata":output}

def download_registered_files(metadata):
    """Use official public version/file links when the dataset ZIP returns 401.

    Individual Dryad file GET is documented separately from dataset ZIP GET.
    Packaging the retrieved files into an in-memory ZIP lets the same file
    identity audit run, but this is NOT the original dataset ZIP checksum.
    """
    links=metadata.get("_links",{})
    version=links.get("stash:version",{}).get("href")
    if not version or not re.match(r"^/api/v2/versions/\d+$",version):
        raise RuntimeError("SOURCE_VERSION_METADATA_UNAVAILABLE")
    url="https://datadryad.org"+version+"/files"
    payload=json.loads(get_bytes(url))
    files=payload.get("_embedded",{}).get("stash:files",[])
    if not files or payload.get("total",len(files))>len(files):
        raise RuntimeError("SOURCE_FILE_PAGINATION_OR_EMPTY_LIST")
    lookup={}
    for row in files:
        name=Path(row.get("path","")).name
        if name in EXPECTED_FILES:
            target=row.get("_links",{}).get("stash:download",{}).get("href")
            if not target or not re.match(r"^/api/v2/files/\d+/download$",target):
                raise RuntimeError("SOURCE_FILE_DOWNLOAD_LINK_UNAVAILABLE")
            try:
                contents=get_bytes("https://datadryad.org"+target)
            except urllib.error.HTTPError as ex:
                raise RuntimeError(f"SOURCE_FILE_ACCESS_RESTRICTED_HTTP_{ex.code}") from ex
            declared_digest=row.get("digest")
            digest_type=row.get("digestType","").lower()
            if declared_digest and digest_type in ("sha256","sha-256"):
                if sha256(contents)!=declared_digest:
                    raise ValueError(f"Source file SHA256 mismatch: {name}")
            elif declared_digest and digest_type=="md5":
                if hashlib.md5(contents).hexdigest()!=declared_digest:
                    raise ValueError(f"Source file MD5 mismatch: {name}")
            lookup[name]=contents
    if set(lookup)!=set(EXPECTED_FILES):
        raise RuntimeError("SOURCE_REGISTERED_GENE_ALIGNMENT_MISSING")
    buff=io.BytesIO()
    with zipfile.ZipFile(buff,"w") as z:
        for key in sorted(lookup):z.writestr(key,lookup[key])
    return buff.getvalue()


def official_original_file_manifest(metadata):
    """Return only public API metadata when sequence bytes are inaccessible.

    The registered Dryad digest and size are *author/source statements*, not
    locally computed checksums. All claims about intact coding/trait reversal
    remain excluded until original bytes pass an independent source check.
    """
    links=metadata.get("_links",{})
    version=links.get("stash:version",{}).get("href")
    if not version or not re.match(r"^/api/v2/versions/\d+$",version):
        raise RuntimeError("SOURCE_VERSION_METADATA_UNAVAILABLE")
    response=json.loads(get_bytes("https://datadryad.org"+version+"/files"))
    records=response.get("_embedded",{}).get("stash:files",[])
    if not records or response.get("total",len(records))>len(records):
        raise RuntimeError("SOURCE_FILE_PAGINATION_OR_EMPTY_LIST")
    rows={}
    for row in records:
        name=Path(row.get("path","")).name
        if name in EXPECTED_FILES:
            link=row.get("_links",{}).get("stash:download",{}).get("href")
            if not link or not re.match(r"^/api/v2/files/\d+/download$",link):
                raise RuntimeError("SOURCE_FILE_ID_MISSING")
            rows[name]={"source_path":row.get("path"),
                        "source_declared_bytes":row.get("size"),
                        "source_declared_digest":row.get("digest"),
                        "source_declared_digest_type":row.get("digestType"),
                        "api_file_url":"https://datadryad.org"+link,
                        "raw_bytes_accessed":False}
    if set(rows)!=set(EXPECTED_FILES):
        raise RuntimeError("SOURCE_ALIGNMENT_MANIFEST_INCOMPLETE")
    return rows


def blocked_source_receipt(metadata,reason):
    return {"schema":"anto_iochrominae_2016_ho_smith_public_coding_source_v1",
            "source_status":reason,
            "scientific_admission":"SOURCE_ONLY_METADATA_FILES_UNAVAILABLE",
            "registered_dataset_DOI":DOI,
            "dataset_metadata_URL":METADATA_URL,
            "dataset_title":metadata.get("title",""),
            "original_files_metadata_only":official_original_file_manifest(metadata),
            "source_bytes_downloaded":False,
            "per_file_sha256_computed":False,
            "coding_integrity_independently_replayed":False,
            "published_original_author_claim_not_retested":
                "Ho and Smith (2016) reported intact coding regions in sampled unpigmented lineages.",
            "phenotypic_regain_or_fitness_estimated":False,
            "no_natural_trait_reexpression_event_count_assigned":True,
            "remediation":"Requires legitimate access to original raw Dryad files; rerun raw-byte verifier without changing source DOI."}


def public_original_material(metadata):
    try:
        return get_bytes(DOWNLOAD_URL),"original_dataset_zip"
    except urllib.error.HTTPError as err:
        if err.code not in (401,403):
            raise
        return download_registered_files(metadata),"official_registered_file_API_fallback"


def evaluate(metadata,archive,method="original_dataset_zip"):
    if metadata.get("identifier")!="doi:"+DOI:
        raise ValueError("Metadata did not identify declared Dryad dataset")
    if metadata.get("versionStatus","").lower() not in ("published","public",""):
        raise ValueError("Dataset is not confirmed as publicly published")
    info=inspect_zip(archive)
    return {"schema":"anto_iochrominae_2016_ho_smith_public_coding_source_v1",
            "registered_dataset_DOI":DOI,
            "data_owner":"Ho and Smith (2016) original authors",
            "dataset_title":metadata.get("title",""),
            "version":metadata.get("versionNumber"),
            "dataset_metadata_URL":METADATA_URL,
            "original_archive_URL":DOWNLOAD_URL,
            "archive_mode":method,
            "original_zip_sha256":sha256(archive) if method=="original_dataset_zip" else None,
            "repackaged_file_bundle_sha256":sha256(archive) if method!="original_dataset_zip" else None,
            "source_files":info,
            "published_flower_bud_cDNA_sequences_22_ingroups_and_6_outgroups":"AUTHOR_SOURCE_STATEMENT_NOT_RECOUNTED",
            "coding_integrity_result_already_published":"Ho and Smith 2016 report intact Chi/F3h/Dfr coding in their unpigmented lineages",
            "original_sequences_rehosted":False,
            "user_generated_cirsium_results":False,
            "outcomes_causally_measured_here":[],
            "biological_inference_level":"SOURCE_FILE_RETRIEVAL_AND_HEADER_ONLY",
            "no_ancestral_regain_or_phenotypic_persistence_inferred":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    meta=json.loads(get_bytes(METADATA_URL))
    try:
        archive,method=public_original_material(meta)
    except RuntimeError as ex:
        if str(ex) not in ("SOURCE_FILE_ACCESS_RESTRICTED_HTTP_401",
                           "SOURCE_FILE_ACCESS_RESTRICTED_HTTP_403"):
            raise
        result=blocked_source_receipt(meta,str(ex))
        status="IOCHROMINAE_DRYAD_SOURCE_DOWNLOAD_BLOCKED"
    else:
        result=evaluate(meta,archive,method)
        result["source_status"]="SOURCE_BYTES_VERIFIED"
        result["source_bytes_downloaded"]=True
        status="IOCHROMINAE_DRYAD_SOURCE_BYTES_VERIFIED"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(status,json.dumps({
        "source_status":result["source_status"],
        "original_bytes_available":result["source_bytes_downloaded"],
        "no_biological_return_claim":True,
        "file_names":list(result.get("source_files",{}).get("target_gene_file_metadata",
                         result.get("original_files_metadata_only",{})))
    },sort_keys=True))

if __name__=="__main__":
    main()
