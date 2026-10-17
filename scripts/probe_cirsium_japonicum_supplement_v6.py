"""Bounded retrieval probe for *Cirsium japonicum* 2020 open supplement ZIP.

The 2020 Park et al. ACS Omega source exposes supplementary annotated-unigene
text and XLSX through PMC/Europe PMC. This optional network probe examines
FILE ACCESSIBILITY and candidate annotation text only; it NEVER produces
Cirsium brevicaule/irumtiense genotypes or gene-function claims.

Default mode: no network, reports NOT_REQUESTED.
Public mode: --fetch-public (one allowlisted HTTPS endpoint, max 50MB).
Local mode: --zip-file ... for a legitimately downloaded supplementary ZIP.
All result states are bounded and do not imply evidence from missing files.
"""
import argparse
import io
import json
import re
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

SOURCE = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7675961/supplementaryFiles"
EXPECTED = {
    "ao0c04001_si_002.txt": "annotated_unigene_text",
    "ao0c04001_si_003.xlsx": "unigene_annotation_spreadsheet",
    "ao0c04001_si_001.pdf": "published_expression_figures",
}
PATTERN = re.compile(r"\b(?:DFR|ANS|LDOX|dihydroflavonol.{0,12}reductase|anthocyanidin.{0,12}synthase)\b",re.I)


def inspect_zip(data, max_uncompressed=180_000_000):
    """Inspect ZIP central directory and a bounded text scan; never execute files."""
    with zipfile.ZipFile(io.BytesIO(data), "r") as z:
        members=z.infolist()
        if len(members)>50:
            raise ValueError("Unexpectedly many supplemental archive members")
        sizes=[m.file_size for m in members]
        if sum(sizes)>max_uncompressed:
            raise ValueError("Archive decompressed bytes exceed resource budget")
        files=[]
        for m in members:
            name=Path(m.filename).name
            if not name or m.is_dir():
                continue
            files.append({"name":name,"size_bytes":m.file_size,
                          "matching_expected_supplement":name in EXPECTED})
        report={"archive_is_zip":True,"members":files,
                "compressed_bytes":len(data),
                "files_found":sorted({Path(m.filename).name for m in members}
                                     & set(EXPECTED)),
                "annotation_text_scan_status":"NOT_FOUND",
                "potential_gene_annotation_text_hits":0,
                "annotation_examples":[]}
        target=next((m for m in members if Path(m.filename).name == "ao0c04001_si_002.txt"),None)
        if target is not None:
            n=0
            seen=0
            header_mode=None
            samples=[]
            with z.open(target,"r") as raw:
                for b in raw:
                    n+=1
                    if n > 250000:
                        report["annotation_text_scan_status"]="TRUNCATED_AT_250000_LINES"
                        break
                    s=b.decode("utf-8",errors="replace").strip()
                    if header_mode is None and s:
                        header_mode="FASTA_FORMAT" if s.startswith(">") else "TABULAR_OR_OTHER_TEXT"
                    if PATTERN.search(s):
                        seen+=1
                        if len(samples)<6:
                            samples.append(s[:200])
            else:
                report["annotation_text_scan_status"]="SCANNED_ALL_LINES"
            report.update(annotation_lines_scanned=n,
                          annotation_text_structure=header_mode,
                          potential_gene_annotation_text_hits=seen,
                          annotation_examples=samples)
        report["biological_sequence_assignment"]="NOT_VALIDATED"
        report["focal_ryukyu_evidence"]="NOT_TESTED"
        return report


def retrieve_source(max_bytes=50_000_000,timeout=25):
    """One bounded GET; on ordinary blocks returns a recorded retrieval failure."""
    request=urllib.request.Request(SOURCE,headers={
        "User-Agent":"anto-open-literature-provenance/1.0 (academic metadata audit)",
        "Accept":"application/zip, application/octet-stream"
    })
    try:
        with urllib.request.urlopen(request,timeout=timeout) as response:
            mime=response.headers.get("Content-Type","")
            length=response.headers.get("Content-Length")
            if length and int(length)>max_bytes:
                return {"retrieval_status":"TOO_LARGE","reported_bytes":int(length),
                        "max_bytes":max_bytes,"url":SOURCE}
            b=response.read(max_bytes+1)
            if len(b)>max_bytes:
                return {"retrieval_status":"TOO_LARGE","bytes_read":len(b),
                        "max_bytes":max_bytes,"url":SOURCE}
            if not zipfile.is_zipfile(io.BytesIO(b)):
                return {"retrieval_status":"RESPONSE_NOT_ZIP",
                        "bytes_read":len(b),"content_type":mime,
                        "url":SOURCE,"message":"Blocked or incompatible content, NOT absence of published supplementary data"}
            return {"retrieval_status":"ARCHIVE_OBTAINED",
                    "url":SOURCE,"content_type":mime,
                    "inspection":inspect_zip(b)}
    except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,
            OSError,ValueError) as exc:
        return {"retrieval_status":"UNAVAILABLE_IN_THIS_RUN",
                "url":SOURCE,
                "error_type":type(exc).__name__,
                "message":str(exc)[:200],
                "not_evidence_of_absence":True}


def run(fetch_public=False,zip_file=None):
    if fetch_public and zip_file is not None:
        raise ValueError("One input source at a time")
    if zip_file is not None:
        b=Path(zip_file).read_bytes()
        if len(b)>50_000_000:
            raise ValueError("Local supplement exceeds bounded archive budget")
        return {"retrieval_status":"LOCAL_PUBLIC_SUPPLEMENT",
                "inspection":inspect_zip(b),"source_accession":"PMC7675961"}
    if fetch_public:
        return retrieve_source()
    return {"retrieval_status":"NOT_REQUESTED",
            "note":"Run --fetch-public to perform an explicit bounded network probe.",
            "source_accession":"PMC7675961"}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--fetch-public",action="store_true")
    p.add_argument("--zip-file",type=Path)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    summary={"version":"anto_cirsium_supplement_probe_v6",
             "not_focal_cirsium_regain_evidence":True,**run(args.fetch_public,args.zip_file)}
    s=json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding="utf-8")
    print(s)
